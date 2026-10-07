"""Gate di validazione G01-G15.

I gate sono i controlli che precedono l'esecuzione: servono a fermare un'analisi
prima che venga allocato qualunque calcolo, non a metà di un'elaborazione che
dura ore.

Il modulo implementa l'insieme completo dei quindici gate previsti (G01-G15):
dalla coerenza interna della configurazione (G15), alla verifica di esistenza e
integrità degli ingressi e delle tabelle ISA-Tab (G01-G06), fino al controllo di
layout single-end (G07), plausibilità dei lotti (G08), compatibilità delle
lunghezze e assenza primer con controllo positivo del motivo conservato (G09,
G10), classificazione controlli (G11), checksum del database tassonomico (G12),
grammatica dei file FASTQ (G13) e risorse hardware disponibili (G14).

G15 intercetta una classe di errori che la validazione parametro per parametro
non può cogliere: combinazioni di valori singolarmente validi ma insensati messi
insieme. Una soglia di avviso più alta della soglia di arresto è fatta di due
numeri plausibili e di una coppia che non lo è.

**Rapporto con il catalogo.** I codici ``E-G15-*`` e ``E-S0-*`` non sono un
insieme a parte: sono voci di :mod:`amplicon16s.errors.catalog` come tutte le
altre, e da lì prendono messaggio operativo e categoria di gestione. Qui resta
la sola informazione che è propria dei gate: quali parametri ciascun controllo
sorveglia, la logica di verifica e la raccolta delle violazioni.

**Rapporto con lo schema.** Quattro dei controlli di coerenza sono già
garantiti dai vincoli dichiarati in :mod:`amplicon16s.config.schema`. G15 non li
riscrive: li esegue attraverso la validazione dello schema e si limita ad
attribuire a ciascuno il proprio codice di errore. Il registro :data:`CONTROLLI`
dice per ognuno dove è implementato, così la divisione del lavoro è leggibile
nel codice e non solo nelle note.
"""

from __future__ import annotations

import csv
import hashlib
import os
import re
import shutil
import statistics
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Final

from amplicon16s.config.resolve import ConfigRisolta, risolvi
from amplicon16s.config.schema import (
    INTESTAZIONE_MANCANTI,
    Config,
    ErroreConfigurazione,
    processori_disponibili,
    valida,
)
from amplicon16s.errors.catalog import CATALOGO, Categoria, VoceCatalogo, voce
from amplicon16s.errors.exceptions import ErrorePipeline
from amplicon16s.io_layer.reads import StatisticheFile, espandi_iupac, scansiona
from amplicon16s.metadata.crosswalk import Analisi, analizza
from amplicon16s.metadata.models import CLASSI_CONTROLLATE, ClasseCampione, Inventario
from amplicon16s.metadata.tabelle import intestazione as _intestazione
from amplicon16s.metadata.tabelle import (
    nomi_in_collisione,
    tabella_delle_cellule,
    tabella_di_studio,
    valori_non_tabellari,
)

__all__ = [
    "Avviso",
    "CONTROLLI",
    "Contesto",
    "Controllo",
    "ErroreGate",
    "GATE_METADATI",
    "Violazione",
    "esegui_g15",
    "esegui_gate_metadati",
    "rifiuto_di_g15",
]


@dataclass(frozen=True)
class Controllo:
    """Un controllo di coerenza, identificato da un codice del catalogo.

    Messaggio e categoria non sono ripetuti qui: vivono nel catalogo degli
    errori, che è la sola fonte di verità dei codici. Questa classe aggiunge
    ciò che il catalogo non sa: su quali parametri il controllo vigila e chi
    lo esegue materialmente.
    """

    codice: str
    #: Parametri la cui combinazione il controllo sorveglia.
    parametri: tuple[str, ...]
    #: ``"schema"`` se il vincolo è già dichiarato nello schema e G15 lo riusa,
    #: ``"gate"`` se è G15 a verificarlo.
    implementato_da: str

    @property
    def voce(self) -> VoceCatalogo:
        """La voce del catalogo corrispondente al codice del controllo."""
        return voce(self.codice)

    @property
    def descrizione(self) -> str:
        """La sintesi del controllo, presa dal catalogo."""
        return self.voce.sintesi

    @property
    def categoria(self) -> Categoria:
        """La categoria di gestione del controllo, presa dal catalogo."""
        return self.voce.categoria


CONTROLLI: Final[tuple[Controllo, ...]] = (
    Controllo(
        "E-G15-02",
        ("decontam.threshold",),
        "gate",
    ),
    Controllo(
        "E-G15-03",
        ("asv.len_min", "asv.len_max", "asv.len_tol", "filter.truncLen"),
        "gate",
    ),
    Controllo(
        "E-G15-04",
        ("prev.min_fraction",),
        "schema",
    ),
    Controllo(
        "E-G15-05",
        ("qc.warn_frac_chimeric", "qc.stop_frac_chimeric"),
        "schema",
    ),
    Controllo(
        "E-G15-06",
        ("retry.max_attempts",),
        "schema",
    ),
    Controllo(
        "E-G15-08",
        ("asv.len_min", "asv.len_max"),
        "schema",
    ),
    Controllo(
        "E-G15-09",
        ("retry.whitelist",),
        "gate",
    ),
    Controllo(
        "E-G15-10",
        (),
        "schema",
    ),
    Controllo(
        "E-G15-11",
        ("io.batch_table", "out.batch_columns", "err.batch_column", "decontam.batch_column",
         "decontam.mode", "meta.batch_key_column", "meta.batch_module_column"),
        "gate",
    ),
    Controllo(
        "E-G15-12",
        ("ctrl.positive_values", "katharoseq.target_taxon", "katharoseq.cell_count_column"),
        "gate",
    ),
    Controllo(
        "E-G15-13",
        ("out.study_columns", "out.batch_columns"),
        "gate",
    ),
    Controllo(
        "E-G15-14",
        ("ctrl.blank_override_values", "ctrl.blank_override_column",
         "meta.study_sample_id_column", "io.study_table",
         "meta.module_regex", "meta.module_column"),
        "gate",
    ),
    Controllo(
        "E-G15-99",
        (),
        "schema",
    ),
    # --- S0: crosswalk, inventario e classificazione dei controlli ---------
    Controllo(
        "E-S0-03",
        ("io.assay_table", "io.study_table", "meta.sample_id_column"),
        "gate",
    ),
    Controllo(
        "E-S0-04",
        ("io.accession_regex", "io.fastq_glob", "meta.accession_column"),
        "gate",
    ),
    Controllo(
        "E-S0-05",
        ("io.accession_regex", "io.fastq_dir", "io.assay_table"),
        "gate",
    ),
    Controllo(
        "E-S0-06",
        ("io.fastq_dir", "io.fastq_glob", "io.assay_table"),
        "gate",
    ),
    Controllo(
        "E-S0-11",
        ("ctrl.column", "ctrl.biological_values", "ctrl.positive_values",
         "ctrl.blank_values"),
        "gate",
    ),
)

_PER_CODICE: Final[dict[str, Controllo]] = {c.codice: c for c in CONTROLLI}


@dataclass(frozen=True)
class Violazione:
    """Un controllo fallito, con il dettaglio di cosa è in conflitto con cosa."""

    codice: str
    dettaglio: str

    @property
    def voce(self) -> VoceCatalogo:
        """La voce del catalogo corrispondente al codice violato."""
        return voce(self.codice)

    @property
    def azione(self) -> str:
        """Che cosa fare, secondo il catalogo."""
        return self.voce.azione

    def __str__(self) -> str:
        return f"[{self.codice}] {self.dettaglio}"


class ErroreGate(ErrorePipeline):
    """Un gate ha respinto la configurazione.

    Raccoglie tutte le violazioni invece di fermarsi alla prima: chi corregge
    una configurazione vuole l'elenco completo.

    Deriva da :class:`ErrorePipeline` perché un gate respinto è un errore della
    pipeline come gli altri e va trattato allo stesso modo da chi lo cattura e
    da chi lo registra; quel che aggiunge è l'insieme delle violazioni, che un
    errore a codice singolo non potrebbe portare.
    """

    def __init__(self, gate: str, violazioni: Sequence[Violazione]) -> None:
        self.gate = gate
        self.violazioni = tuple(violazioni)
        # Il codice con cui l'esecuzione termina e' quello della prima
        # violazione; le altre restano nell'elenco e nel log.
        primo = self.violazioni[0].codice if self.violazioni else "E-G15-99"
        super().__init__(primo, gate=gate, violazioni=[v.codice for v in self.violazioni])

    def _testo(self) -> str:
        """Il messaggio con il nome del gate e l'elenco delle violazioni, una per riga.
        """
        return "{} ha respinto la configurazione: {} violazion{}\n{}".format(
            self.gate,
            len(self.violazioni),
            "e" if len(self.violazioni) == 1 else "i",
            "\n".join(f"  - {v}" for v in self.violazioni),
        )


# --------------------------------------------------------------------------- #
# Attribuzione dei codici ai vincoli già dichiarati nello schema               #
# --------------------------------------------------------------------------- #
# Lo schema respinge da sé queste combinazioni. La tabella serve solo a dire
# quale controllo di G15 ciascun rifiuto rappresenta, così il codice di errore
# è quello del controllo e non uno generico.

_ATTRIBUZIONI: Final[tuple[tuple[tuple[str, ...], frozenset[str], str], ...]] = (
    (
        ("prev", "min_fraction"),
        frozenset({"greater_than_equal", "less_than_equal", "greater_than", "less_than"}),
        "E-G15-04",
    ),
    (("qc",), frozenset({"value_error"}), "E-G15-05"),
    (
        ("retry", "max_attempts"),
        frozenset({"greater_than", "greater_than_equal"}),
        "E-G15-06",
    ),
)


def _codice_per_problema(problema: str) -> str:
    """Attribuisce un codice di controllo a un problema segnalato dallo schema."""
    if problema.startswith(INTESTAZIONE_MANCANTI):
        return "E-G15-10"
    if "parametro derivato" in problema:
        return "E-G15-08"
    # La regola di riclassificazione senza la sua colonna: lo schema la respinge
    # da se', e il rifiuto e' quello del controllo sulle dipendenze fra parametri.
    if problema.startswith("ctrl:") and "blank_override_column e' nullo" in problema:
        return "E-G15-14"
    for percorso, _tipi, codice in _ATTRIBUZIONI:
        if problema.startswith(".".join(percorso) + ":"):
            return codice
    return "E-G15-99"


def _violazione_da_schema(problema: str) -> Violazione:
    """Trasforma un problema dello schema nella violazione del controllo relativo.

    Il messaggio di pydantic nomina il solo parametro che ha sbagliato. Il
    controllo, invece, sorveglia una combinazione: il contesto viene aggiunto
    perche' chi legge sappia quali parametri sono in conflitto fra loro e non
    solo quale valore e' stato rifiutato.
    """
    codice = _codice_per_problema(problema)
    controllo = _PER_CODICE[codice]

    if not controllo.parametri or codice == "E-G15-14":
        return Violazione(codice, problema)

    return Violazione(
        codice,
        "{} : controllo {} su {}: {}".format(
            problema,
            codice,
            ", ".join(controllo.parametri),
            controllo.descrizione,
        ),
    )


# --------------------------------------------------------------------------- #
# Controlli propri di G15                                                      #
# --------------------------------------------------------------------------- #


def _controlla_coerenza(risolta: ConfigRisolta) -> list[Violazione]:
    """Esegue i controlli che lo schema non copre, sulla configurazione risolta."""
    config = risolta.config
    derivati = risolta.derivati
    violazioni: list[Violazione] = []

    # E-G15-02 : lo schema ammette [0, 1]; agli estremi la soglia e' degenere.
    soglia = config.decontam.threshold
    if not (0.0 < soglia < 1.0):
        violazioni.append(
            Violazione(
                "E-G15-02",
                f"decontam.threshold ({soglia}) deve essere strettamente compreso fra "
                f"0 e 1: a 0 nessuna variante verrebbe mai classificata come "
                f"contaminante, a 1 lo sarebbero tutte",
            )
        )

    # E-G15-03 : l'intervallo delle lunghezze ammesse deve essere sensato.
    if derivati.asv_len_min > derivati.asv_len_max:
        violazioni.append(
            Violazione(
                "E-G15-03",
                f"asv.len_min ({derivati.asv_len_min}) supera asv.len_max "
                f"({derivati.asv_len_max}): l'intervallo delle lunghezze ammesse "
                f"sarebbe vuoto",
            )
        )
    elif derivati.asv_len_min < 1:
        violazioni.append(
            Violazione(
                "E-G15-03",
                f"asv.len_min ({derivati.asv_len_min}) non e' positivo: asv.len_tol "
                f"({config.asv.len_tol}) e' troppo grande rispetto a filter.truncLen "
                f"({config.filter.truncLen}) meno filter.trimLeft "
                f"({config.filter.trimLeft})",
            )
        )

    # E-G15-09 : la whitelist e' l'autorita' su quali codici si ritentano, ma
    # puo' solo restringere l'elenco del catalogo: un codice che il catalogo
    # non ammette al retry verrebbe ritentato cambiando un'assunzione.
    for codice in config.retry.whitelist:
        if codice not in CATALOGO:
            violazioni.append(
                Violazione(
                    "E-G15-09",
                    f"retry.whitelist contiene {codice}, che non esiste nel catalogo",
                )
            )
        elif not CATALOGO[codice].ammette_retry:
            violazioni.append(
                Violazione(
                    "E-G15-09",
                    f"retry.whitelist contiene {codice}, classificato nel catalogo "
                    f"come {CATALOGO[codice].categoria.value}: non ammette il retry",
                )
            )

    # E-G15-11 : le colonne del file del lotto hanno senso solo se il file c'e'.
    # Senza, una colonna dichiarata resterebbe vuota per ogni campione, e lo si
    # scoprirebbe a calcolo avviato (S3, S10) o mai.
    if config.io.batch_table is None:
        dichiarate = [
            nome for nome, valore in (
                ("out.batch_columns", config.out.batch_columns),
                ("err.batch_column", config.err.batch_column),
                ("decontam.batch_column", config.decontam.batch_column),
                ("meta.batch_key_column", config.meta.batch_key_column),
                ("meta.batch_module_column", config.meta.batch_module_column),
            ) if valore
        ]
        if dichiarate:
            violazioni.append(
                Violazione(
                    "E-G15-11",
                    "io.batch_table e' nullo ma {} {} un valore: senza il file del "
                    "lotto vanno dichiarati vuoti o nulli".format(
                        ", ".join(dichiarate), "ha" if len(dichiarate) == 1 else "hanno"
                    ),
                )
            )
    elif config.meta.batch_key_column is None:
        violazioni.append(
            Violazione(
                "E-G15-11",
                "io.batch_table e' indicato ma meta.batch_key_column e' nullo: senza la "
                "colonna con la chiave del campione le righe del file non si agganciano",
            )
        )

    # La decontaminazione per piastra confronta ogni campione con i negativi
    # della sua piastra: senza la colonna della piastra non c'e' raggruppamento.
    if config.decontam.mode == "batch" and config.decontam.batch_column is None:
        violazioni.append(
            Violazione(
                "E-G15-11",
                "decontam.mode e' batch ma decontam.batch_column e' nullo: senza la "
                "colonna della piastra la decontaminazione per piastra non ha gruppi; "
                "dichiara decontam.mode aggregate",
            )
        )

    # E-G15-12 : controlli positivi dichiarati, ma non valutabili.
    if config.ctrl.positive_values:
        nulli = [
            nome for nome, valore in (
                ("katharoseq.target_taxon", config.katharoseq.target_taxon),
                ("katharoseq.cell_count_column", config.katharoseq.cell_count_column),
            ) if valore is None
        ]
        if nulli:
            violazioni.append(
                Violazione(
                    "E-G15-12",
                    "ctrl.positive_values elenca {} ma {} nullo: i controlli positivi "
                    "non sarebbero valutabili".format(
                        config.ctrl.positive_values,
                        " e ".join(nulli) + (" e'" if len(nulli) == 1 else " sono"),
                    ),
                )
            )

    # E-G15-14 : un parametro che ha effetto solo insieme a un altro. Dichiarato
    # da solo verrebbe ignorato in silenzio, e chi lo ha scritto crederebbe il
    # contrario. (La regola di riclassificazione senza colonna la respinge gia'
    # lo schema, con lo stesso codice.)
    if config.meta.study_sample_id_column is not None and config.io.study_table is None:
        violazioni.append(
            Violazione(
                "E-G15-14",
                "meta.study_sample_id_column e' dichiarato "
                f"({config.meta.study_sample_id_column!r}) ma io.study_table e' nullo: "
                "senza tabella di studio il nome del campione si legge dalla tabella di "
                "assay, in meta.sample_id_column, e il parametro non avrebbe effetto",
            )
        )
    if config.meta.module_regex is not None and config.meta.module_column is None:
        violazioni.append(
            Violazione(
                "E-G15-14",
                f"meta.module_regex e' dichiarata ({config.meta.module_regex!r}) ma "
                "meta.module_column e' nullo: l'espressione si applica alla posizione "
                "del campione, e senza la colonna che la dichiara non deriverebbe alcun "
                "modulo",
            )
        )

    # E-G15-13 : due colonne richieste non possono avere lo stesso nome
    # nell'oggetto, ne' quello di una colonna dell'inventario.
    for originale, nome in nomi_in_collisione(
        [*config.out.study_columns, *config.out.batch_columns]
    ):
        violazioni.append(
            Violazione(
                "E-G15-13",
                f"la colonna {originale!r} (out.study_columns o out.batch_columns) "
                f"diventerebbe {nome!r} nell'oggetto, nome gia' usato da un'altra colonna",
            )
        )

    return violazioni


# --------------------------------------------------------------------------- #
# Esecuzione del gate                                                          #
# --------------------------------------------------------------------------- #


def rifiuto_di_g15(errore: ErroreConfigurazione) -> ErroreGate:
    """Il rifiuto di G15 per una configurazione che lo schema non accetta: ogni
    problema diventa la violazione del controllo che lo sorveglia, con il suo
    codice. I parametri obbligatori non dichiarati sono una sola violazione,
    ``E-G15-10``, che li elenca tutti.
    """
    return ErroreGate("G15", [_violazione_da_schema(p) for p in errore.problemi])


def esegui_g15(dati: Mapping[str, Any]) -> ConfigRisolta:
    """Esegue G15 su una configurazione non ancora validata.

    Restituisce la configurazione risolta nella sua parte statica. Solleva
    :class:`ErroreGate` se un controllo fallisce, con il codice del controllo
    e il dettaglio di quali parametri sono in conflitto.
    """
    try:
        config = valida(dict(dati))
    except ErroreConfigurazione as errore:
        raise rifiuto_di_g15(errore) from errore

    risolta = risolvi(config)

    violazioni = _controlla_coerenza(risolta)
    if violazioni:
        raise ErroreGate("G15", violazioni)

    return risolta


# =========================================================================== #
# Gate sui metadati: G04, G05, G06, G03, G11                                  #
# =========================================================================== #
# Sorvegliano la costruzione dell'inventario dei campioni. L'ordine non e'
# arbitrario: senza un accession estraibile non si puo' verificarne
# l'univocita', senza accession univoci non si possono confrontare gli
# insiemi, e cosi' via. Ciascun gate si ferma appena trova violazioni, perche'
# proseguire su dati gia' noti come incoerenti produrrebbe errori derivati che
# confonderebbero la diagnosi.

_MAX_ELENCATI = 10


def _elenca(voci: Sequence[str]) -> str:
    """Elenca alcune voci, dicendo quante ne restano.

    Un elenco di centinaia di accession renderebbe illeggibile il messaggio;
    nasconderne il numero renderebbe impossibile capire la portata del
    problema.
    """
    mostrate = list(voci[:_MAX_ELENCATI])
    resto = len(voci) - len(mostrate)
    testo = ", ".join(mostrate)
    return f"{testo} e altri {resto}" if resto > 0 else testo


def _g04_accession_estraibile(analisi: Analisi) -> list[Violazione]:
    """L'accession dev'essere estraibile da ogni nome, e non ambiguo."""
    violazioni = []
    for nome, motivo in analisi.file_senza_accession:
        violazioni.append(
            Violazione(
                "E-S0-04",
                f"dal nome del file {nome!r} non si ricava un accession: {motivo}. "
                f"io.accession_regex non corrisponde a come sono nominati i file",
            )
        )
    for riferimento, motivo in analisi.assay_senza_accession:
        violazioni.append(
            Violazione(
                "E-S0-04",
                f"dal valore {riferimento!r} di meta.accession_column non si ricava "
                f"un accession: {motivo}",
            )
        )

    if analisi.arricchimento_senza_colonna is not None:
        violazioni.append(
            Violazione(
                "E-S0-04",
                "il file indicato da io.batch_table non ha la colonna "
                f"{analisi.arricchimento_senza_colonna!r} dichiarata in "
                "meta.batch_key_column. Quella colonna deve contenere l'accession "
                "dell'esperimento, perche' e' la chiave con cui la pipeline aggancia "
                "ogni riga: il nome del campione si ripete fra repliche e aggancerebbe "
                "la replica sbagliata. Se il file non ce l'ha, ricavala dalla tabella "
                "di assay indicata da io.assay_table, associando il nome del campione "
                "in meta.sample_id_column e leggendo l'accession da "
                "meta.accession_column; oppure togli io.batch_table e la pipeline "
                "procedera' senza lotto",
            )
        )

    for grezzo, motivo in analisi.arricchimento_senza_accession:
        violazioni.append(
            Violazione(
                "E-S0-04",
                f"nel file indicato da io.batch_table, dal valore {grezzo!r} della "
                f"colonna meta.batch_key_column non si ricava un accession: {motivo}",
            )
        )

    return violazioni


def _g05_accession_univoci(analisi: Analisi) -> list[Violazione]:
    """Lo stesso accession non puo' appartenere a due file o a due campioni."""
    violazioni = []
    for accession, file in sorted(analisi.accession_ripetuti_nei_file.items()):
        violazioni.append(
            Violazione(
                "E-S0-05",
                f"l'accession {accession} compare in {len(file)} file: "
                f"{_elenca(file)}. Se io.accession_regex ha un gruppo di cattura la "
                f"chiave e' il primo gruppo: un gruppo che serve solo a raggruppare va "
                f"scritto senza cattura, (?:...)",
            )
        )
    for accession, nomi in sorted(analisi.accession_ripetuti_nell_assay.items()):
        violazioni.append(
            Violazione(
                "E-S0-05",
                f"l'accession {accession} compare in {len(nomi)} righe della "
                f"tabella di assay, per i campioni {_elenca(nomi)}",
            )
        )
    return violazioni


def _g06_insiemi_simmetrici(analisi: Analisi) -> list[Violazione]:
    """I file e le righe dell'assay devono coprire gli stessi accession."""
    violazioni = []
    if analisi.solo_nei_file:
        violazioni.append(
            Violazione(
                "E-S0-06",
                f"{len(analisi.solo_nei_file)} accession compaiono fra i file ma non "
                f"nella tabella di assay: {_elenca(analisi.solo_nei_file)}",
            )
        )
    if analisi.solo_nell_assay:
        violazioni.append(
            Violazione(
                "E-S0-06",
                f"{len(analisi.solo_nell_assay)} accession compaiono nella tabella di "
                f"assay ma non fra i file: {_elenca(analisi.solo_nell_assay)}",
            )
        )
    return violazioni


def _g03_join_ristretto(analisi: Analisi) -> list[Violazione]:
    """Il join verso la tabella di studio non deve alterare l'insieme.

    La restrizione e' garantita per costruzione (si itera sulle righe
    dell'assay), ma resta un modo in cui puo' rompersi: un nome di campione
    che nella tabella di studio compare piu' volte, perche' riusato da un
    altro assay. In quel caso il join moltiplicherebbe le righe dell'assay, e
    quale delle due sia quella giusta non e' deducibile.
    """
    violazioni = []
    for nome, quante in sorted(analisi.nomi_ambigui_nello_studio.items()):
        violazioni.append(
            Violazione(
                "E-S0-03",
                f"il campione {nome!r} corrisponde a {quante} righe della tabella "
                f"campioni di studio: il join non sarebbe piu' ristretto a una riga "
                f"per campione dell'assay",
            )
        )
    if analisi.nomi_assenti_nello_studio:
        violazioni.append(
            Violazione(
                "E-S0-03",
                "{} campioni dell'assay non hanno una riga nella tabella campioni di "
                "studio, quindi restano senza classe: {}".format(
                    len(analisi.nomi_assenti_nello_studio),
                    _elenca(analisi.nomi_assenti_nello_studio),
                ),
            )
        )
    return violazioni


def _g11_classi_complete(analisi: Analisi) -> list[Violazione]:
    """Ogni campione ricade in una e una sola classe dichiarata."""
    violazioni = []
    for materiale, nomi in sorted(analisi.materiali_non_mappati.items()):
        violazioni.append(
            Violazione(
                "E-S0-11",
                "il valore {!r} di ctrl.column non compare in ctrl.biological_values, "
                "ctrl.positive_values ne' ctrl.blank_values, e riguarda {} campioni: "
                "{}".format(materiale or "<vuoto>", len(nomi), _elenca(nomi)),
            )
        )
    return violazioni


#: I gate sui metadati, nell'ordine in cui vanno eseguiti.
GATE_METADATI: Final[tuple[tuple[str, Any], ...]] = (
    ("G04", _g04_accession_estraibile),
    ("G05", _g05_accession_univoci),
    ("G06", _g06_insiemi_simmetrici),
    ("G03", _g03_join_ristretto),
    ("G11", _g11_classi_complete),
)


def esegui_gate_metadati(config: Config) -> Inventario:
    """Costruisce l'inventario dei campioni facendolo passare dai cinque gate.

    Solleva :class:`ErroreGate` al primo gate che trova violazioni, con il
    nome del gate e il codice del catalogo corrispondente.
    """
    analisi = analizza(config)

    for nome, controllo in GATE_METADATI:
        violazioni = controllo(analisi)
        if violazioni:
            raise ErroreGate(nome, violazioni)

    return analisi.inventario()


# =========================================================================== #
# Contesto condiviso fra i gate                                               #
# =========================================================================== #


@dataclass(frozen=True)
class Avviso:
    """Una segnalazione che non ferma l'esecuzione.

    Serve ai controlli di plausibilita': una composizione inconsueta merita di
    essere vista, ma non e' un errore e bloccare sarebbe sbagliato.
    """

    codice: str
    dettaglio: str

    def __str__(self) -> str:
        return f"[{self.codice}] {self.dettaglio}"


class Contesto:
    """Stato condiviso dai gate, calcolato una volta sola.

    Leggere i metadati e ispezionare le letture costa: ogni gate che ne ha
    bisogno attinge da qui invece di rifare il lavoro. Le proprieta' sono
    calcolate alla prima richiesta, cosi' un gate che fallisce presto non
    paga il costo di cio' che non serve piu'.
    """

    def __init__(self, config: Config) -> None:
        self.config = config
        self._analisi: Analisi | None = None
        self._scansione: dict[str, Any] | None = None

    @property
    def analisi(self) -> Analisi:
        """L'analisi dei metadati (crosswalk), calcolata alla prima richiesta."""
        if self._analisi is None:
            self._analisi = analizza(self.config)
        return self._analisi

    @property
    def inventario(self) -> Inventario:
        """L'inventario dei campioni, ricavato dall'analisi dei metadati."""
        return self.analisi.inventario()

    @property
    def scansione_gia_fatta(self) -> bool:
        """Se la scansione e' gia' stata calcolata da qualche gate."""
        return self._scansione is not None

    @property
    def scansione(self) -> dict[str, StatisticheFile]:
        """Statistiche delle prime letture di ogni file, una passata sola."""
        if self._scansione is None:
            # Con filter.trimLeft maggiore di zero l'inizio delle letture verra'
            # tolto dal filtro: il primer in testa non e' piu' un difetto e non si
            # cerca, e il motivo conservato si cerca dove iniziera' la lettura
            # filtrata.
            taglio = self.config.filter.trimLeft
            self._scansione = scansiona(
                self.analisi.file_per_accession,
                head_reads=self.config.qc.head_reads,
                primer=espandi_iupac(self.config.qc.primer_sequence) if taglio == 0 else None,
                motivo=self.config.qc.conserved_motif,
                inizio_motivo=taglio,
                lunghezza_richiesta=self.config.filter.truncLen,
            )
        return self._scansione


# =========================================================================== #
# Gate di ingresso e di risorse                                               #
# =========================================================================== #


def _g01_ingressi_leggibili(contesto: Contesto) -> tuple[list[Violazione], list[Avviso]]:
    """Ogni ingresso di S0 dichiarato esiste ed e' leggibile.

    Il riferimento tassonomico non e' un ingresso di S0: lo verifica G12, che
    si esegue a ogni avvio come precondizione e ne controlla anche il checksum.
    """
    config = contesto.config
    violazioni: list[Violazione] = []

    attesi: list[tuple[str, Path, bool]] = [
        ("io.fastq_dir", Path(config.io.fastq_dir), True),
        ("io.assay_table", Path(config.io.assay_table), False),
    ]
    if config.io.study_table is not None:
        attesi.append(("io.study_table", Path(config.io.study_table), False))
    if config.io.batch_table is not None:
        attesi.append(("io.batch_table", Path(config.io.batch_table), False))

    for parametro, percorso, e_cartella in attesi:
        if not percorso.exists():
            violazioni.append(
                Violazione("E-S0-01", f"{parametro}: {percorso} non esiste")
            )
            continue
        if e_cartella and not percorso.is_dir():
            violazioni.append(
                Violazione("E-S0-01", f"{parametro}: {percorso} non e' una cartella")
            )
            continue
        if not e_cartella and not percorso.is_file():
            violazioni.append(
                Violazione("E-S0-01", f"{parametro}: {percorso} non e' un file")
            )
            continue
        if not os.access(percorso, os.R_OK):
            violazioni.append(
                Violazione("E-S0-01", f"{parametro}: {percorso} non e' leggibile "
                                      f"dall'utente che esegue la pipeline")
            )

    cartella = Path(config.io.fastq_dir)
    if cartella.is_dir() and not any(cartella.glob(config.io.fastq_glob)):
        violazioni.append(
            Violazione(
                "E-S0-01",
                f"io.fastq_dir ({cartella}) non contiene alcun file che "
                f"corrisponda a io.fastq_glob ({config.io.fastq_glob!r})",
            )
        )
    return violazioni, []


def _colonne_mancanti(
    intestazione: Sequence[str], colonne: Sequence[tuple[str, str]], tabella: str, codice: str
) -> list[Violazione]:
    """Una violazione per ogni colonna dichiarata che non compare, o compare piu'
    volte, nell'intestazione. ``colonne`` sono coppie (parametro, nome).
    """
    violazioni = []
    for dichiarata, nome in colonne:
        occorrenze = list(intestazione).count(nome)
        if occorrenze == 1:
            continue
        stato = (
            "non compare nell'intestazione" if occorrenze == 0
            else f"compare {occorrenze} volte nell'intestazione"
        )
        violazioni.append(
            Violazione(codice, f"la colonna {nome!r}, dichiarata in {dichiarata}, {stato} di {tabella}")
        )
    return violazioni


def _valori_spezzati(
    percorso: Path, colonne: Sequence[tuple[str, str]], tabella: str, codice: str
) -> list[Violazione]:
    """Una violazione se una colonna dichiarata ha valori con una tabulazione o
    un a capo dentro il campo: le tabelle scritte a valle (crosswalk, metadati
    dell'oggetto) non hanno virgolette, e S10 se ne accorgerebbe dopo tutto il
    calcolo.
    """
    trovati = valori_non_tabellari(percorso, [nome for _, nome in colonne])
    if not trovati:
        return []
    return [
        Violazione(
            codice,
            "{} valori di {} contengono una tabulazione o un a capo dentro il campo: {}. "
            "Le tabelle prodotte dalla pipeline non li possono riportare: sostituiscili "
            "con uno spazio nel file dei metadati".format(
                len(trovati), tabella,
                _elenca([f"riga {riga}, colonna {nome!r}" for riga, nome in trovati]),
            ),
        )
    ]


def _g02_tabelle_apribili(contesto: Contesto) -> tuple[list[Violazione], list[Avviso]]:
    """Le tabelle di metadati si aprono e hanno ogni colonna dichiarata.

    Si verifica qui, prima di qualunque calcolo, ogni colonna che la
    configurazione chiede alle tabelle di assay e di studio: quelle con cui si
    costruisce l'inventario e quelle che S10 portera' nell'oggetto
    (``out.study_columns``). Una colonna assente scoperta in S10 costerebbe
    tutte le fasi di calcolo che la precedono.
    """
    config = contesto.config
    violazioni: list[Violazione] = []
    studio, colonna_id = tabella_di_studio(config)

    assay = [("meta.sample_id_column", config.meta.sample_id_column),
             ("meta.accession_column", config.meta.accession_column)]
    di_studio = [
        ("meta.study_sample_id_column" if config.meta.study_sample_id_column
         and config.io.study_table is not None else "meta.sample_id_column", colonna_id),
        ("ctrl.column", config.ctrl.column),
    ]
    if config.meta.module_column is not None:
        di_studio.append(("meta.module_column", config.meta.module_column))
    if config.ctrl.blank_override_column is not None:
        di_studio.append(("ctrl.blank_override_column", config.ctrl.blank_override_column))
    di_studio += [("out.study_columns", c) for c in config.out.study_columns]
    # I livelli dei controlli positivi: se la colonna e' fra quelle richieste
    # e' gia' verificata (qui o in G08); altrimenti si cerca nel file di
    # arricchimento e nella tabella di studio, e S10 la porta nell'oggetto da
    # dove la trova. Assente da entrambe, S11 non avrebbe i livelli.
    cellule = config.katharoseq.cell_count_column
    if cellule is not None and cellule not in (
        *config.out.batch_columns, *config.out.study_columns
    ):
        try:
            origine = tabella_delle_cellule(config)
        except (OSError, UnicodeDecodeError, csv.Error):
            origine = "illeggibile"  # la tabella che non si apre e' gia' una violazione
        if origine == "studio":
            di_studio.append(("katharoseq.cell_count_column", cellule))
        elif origine is None:
            violazioni.append(Violazione(
                "E-S0-02",
                f"la colonna {cellule!r}, dichiarata in katharoseq.cell_count_column, non "
                f"compare ne' nel file di arricchimento (io.batch_table) ne' in {studio.name}",
            ))

    if config.io.study_table is None:
        attese = [("io.assay_table", Path(config.io.assay_table), assay + di_studio[1:])]
    else:
        attese = [("io.assay_table", Path(config.io.assay_table), assay),
                  ("io.study_table", studio, di_studio)]

    for parametro, percorso, colonne in attese:
        try:
            intestazione = _intestazione(percorso)
        except (OSError, UnicodeDecodeError, csv.Error) as guasto:
            violazioni.append(
                Violazione("E-S0-02", f"{parametro}: {percorso} non si apre: {guasto}")
            )
            continue
        if not intestazione:
            violazioni.append(
                Violazione("E-S0-02", f"{parametro}: {percorso} non ha intestazione")
            )
            continue
        # Lo stesso nome dichiarato in due parametri si segnala una volta sola.
        uniche = list(dict((nome, (dichiarata, nome)) for dichiarata, nome in colonne).values())
        mancanti = _colonne_mancanti(
            intestazione, uniche, f"{percorso.name} ({parametro})", "E-S0-02"
        )
        violazioni += mancanti
        if not mancanti:
            try:
                violazioni += _valori_spezzati(
                    percorso, uniche, f"{percorso.name} ({parametro})", "E-S0-02"
                )
            except (OSError, UnicodeDecodeError, csv.Error) as guasto:
                violazioni.append(
                    Violazione("E-S0-02", f"{parametro}: {percorso} non si legge: {guasto}")
                )
    return violazioni, []


#: Il secondo file di una coppia, riconosciuto dal marcatore che precede
#: l'estensione: ``_R2``, ``.R2``, ``_2`` o ``.2``, con un eventuale numero di
#: blocco (``_001``). Ancorato alla fine del nome: un campione che si chiama
#: ``LAB_R2D2`` non e' una lettura inversa. Sono convenzioni diffuse, non una
#: proprieta' di un dataset.
_LETTURA_INVERSA: Final = re.compile(
    r"(?:[._]R2|[._]2)(?:_[0-9]+)?\.(?:fastq|fq)(?:\.gz)?$", re.IGNORECASE
)


def _g07_layout_single_end(contesto: Contesto) -> tuple[list[Violazione], list[Avviso]]:
    """Un solo file per campione, nessun file di lettura inversa, e nessun file
    con le due letture di ogni coppia.

    Il terzo controllo guarda dentro i file, nelle letture gia' ispezionate:
    alcuni archivi pubblici distribuiscono un dataset paired-end con un solo
    file per corsa, che contiene le letture forward e le reverse. Il nome non
    lo dice, e trattato come single-end darebbe varianti di due regioni
    diverse con il doppio delle letture.
    """
    config = contesto.config
    violazioni: list[Violazione] = []

    for accession, nomi in sorted(
        contesto.analisi.accession_ripetuti_nei_file.items()
    ):
        violazioni.append(
            Violazione(
                "E-S0-07",
                f"l'accession {accession} ha {len(nomi)} file: un layout "
                f"single-end ne prevede uno solo per campione",
            )
        )

    inverse = sorted(
        percorso.name
        for percorso in Path(config.io.fastq_dir).glob(config.io.fastq_glob)
        if _LETTURA_INVERSA.search(percorso.name)
    )
    if inverse:
        violazioni.append(
            Violazione(
                "E-S0-07",
                "{} file sembrano letture inverse: {}. Questa pipeline tratta "
                "solo dati single-end".format(len(inverse), _elenca(inverse)),
            )
        )

    intercalati = sorted(
        f"{s.nome} ({s.coppie_nello_stesso_file})"
        for s in contesto.scansione.values()
        if s.valido and s.coppie_nello_stesso_file
    )
    if intercalati:
        violazioni.append(
            Violazione(
                "E-S0-07",
                "{} file contengono le due letture di ogni coppia: {}. Estrai le sole "
                "letture forward (la prima di ogni coppia) in un file per campione e "
                "indica quelli in io.fastq_dir".format(len(intercalati), _elenca(intercalati)),
            )
        )
    return violazioni, []


def _g12_riferimento_verificato(contesto: Contesto) -> tuple[list[Violazione], list[Avviso]]:
    """Il database tassonomico esiste e il suo checksum e' quello dichiarato.

    Se ``tax.ref_bad_taxa`` indica l'elenco dei taxa difettosi del riferimento,
    anche quel file deve esistere: S8 lo legge.
    """
    config = contesto.config
    elenco = config.tax.ref_bad_taxa
    if elenco is not None and not Path(elenco).is_file():
        return [Violazione("E-S0-12", f"tax.ref_bad_taxa: {elenco} non esiste")], []
    percorso = Path(config.tax.ref_fasta)

    if not percorso.is_file():
        return [Violazione("E-S0-12", f"tax.ref_fasta: {percorso} non esiste")], []

    impronta = hashlib.md5()
    try:
        with open(percorso, "rb") as file:
            while blocco := file.read(1024 * 1024):
                impronta.update(blocco)
    except OSError as e:
        return [Violazione("E-S0-12", f"tax.ref_fasta: {percorso} non e' leggibile: {e}")], []
    ottenuto = impronta.hexdigest()

    if ottenuto.lower() != config.tax.ref_md5.lower():
        return [
            Violazione(
                "E-S0-12",
                f"il checksum di {percorso.name} e' {ottenuto}, mentre tax.ref_md5 "
                f"dichiara {config.tax.ref_md5}",
            )
        ], []
    return [], []


def _g14_risorse_disponibili(contesto: Contesto) -> tuple[list[Violazione], list[Avviso]]:
    """Le CPU richieste e lo spazio su disco sono disponibili.

    Con ``run.threads`` nullo (automatico) le fasi usano i processori
    utilizzabili, qualunque sia il loro numero: non c'e' una richiesta da
    confrontare, e il controllo riguarda il solo spazio su disco.
    """
    config = contesto.config
    violazioni: list[Violazione] = []

    disponibili = processori_disponibili()
    if config.run.threads is not None and config.run.threads > disponibili:
        violazioni.append(
            Violazione(
                "E-S0-14",
                f"run.threads vale {config.run.threads} ma sono utilizzabili "
                f"{disponibili} CPU",
            )
        )

    ingresso = sum(
        percorso.stat().st_size
        for percorso in Path(config.io.fastq_dir).glob(config.io.fastq_glob)
    )
    radice = Path(config.io.out_root)
    esistente = radice
    while not esistente.exists() and esistente != esistente.parent:
        esistente = esistente.parent
    libero = shutil.disk_usage(esistente).free

    # Soglia senza parametri: almeno quanto pesano i dati di ingresso. Gli
    # artefatti intermedi di una pipeline 16S sono dello stesso ordine, e una
    # stima piu' fine richiederebbe un parametro da indovinare.
    if libero < ingresso:
        violazioni.append(
            Violazione(
                "E-S0-14",
                f"su {esistente} sono liberi {libero / 1024**3:.1f} GB, meno dei "
                f"{ingresso / 1024**3:.1f} GB occupati dai dati di ingresso",
            )
        )
    return violazioni, []


# =========================================================================== #
# Gate che ispezionano le letture                                             #
# =========================================================================== #
# Tutti e tre attingono alla stessa scansione: le prime qc.head_reads letture
# di ogni file, lette una volta sola.


def _g13_archivi_validi(contesto: Contesto) -> tuple[list[Violazione], list[Avviso]]:
    """Ogni file e' un archivio valido con quattro righe per record."""
    violazioni = []
    for accession in sorted(contesto.scansione):
        statistiche = contesto.scansione[accession]
        if not statistiche.valido:
            violazioni.append(
                Violazione(
                    "E-S0-13",
                    f"{statistiche.nome} ({accession}): {statistiche.errore}",
                )
            )
    return violazioni, []


def _g09_troncamento_compatibile(contesto: Contesto) -> tuple[list[Violazione], list[Avviso]]:
    """Le letture piu' corte di filter.truncLen sono poche.

    Le letture piu' corte del troncamento vengono scartate, non accorciate: un
    valore troppo alto azzera interi campioni senza che nulla lo segnali. Una
    sola lettura corta, pero', non e' un motivo per fermarsi: in un dataset a
    lunghezza variabile ce n'e' quasi sempre qualcuna, e il filtro la toglie
    senza danno. Il gate misura la **frazione** di letture piu' corte **per
    classe**, separatamente per i campioni biologici e per i controlli positivi
    (i controlli negativi, che amplificano poco e male, non contano), e si
    ferma se una delle due supera ``qc.max_frac_short_reads``. Le due classi
    non si riuniscono in una sola frazione pesata sulle letture: la classe
    meno numerosa (di norma i positivi) potrebbe perdere tutte le sue letture
    senza spostare la frazione complessiva.

    La misura e' sulle prime qc.head_reads letture di ogni file, non su tutto
    il file: il limite e' chiuso da S1, che legge tutte le letture e ripete la
    stessa verifica (E-S1-02).
    """
    config = contesto.config
    troncamento = config.filter.truncLen
    massima = config.qc.max_frac_short_reads
    classi = {c.accession: c.classe for c in contesto.inventario}

    violazioni: list[Violazione] = []
    for classe in CLASSI_CONTROLLATE:
        esaminate = corte = 0
        con_corte: list[StatisticheFile] = []
        for accession, s in contesto.scansione.items():
            if not s.valido or classi.get(accession) is not classe:
                continue
            esaminate += s.letture_esaminate
            corte += s.piu_corte
            if s.piu_corte:
                con_corte.append(s)
        frazione = corte / esaminate if esaminate else 0.0
        if frazione <= massima:
            continue
        peggiori = sorted(con_corte, key=lambda s: -s.frazione_corte)
        dettaglio = ", ".join(
            f"{s.nome} ({s.frazione_corte:.1%}, minima {s.lunghezza_minima} bp)"
            for s in peggiori[:5]
        )
        violazioni.append(
            Violazione(
                "E-S0-09",
                f"filter.truncLen vale {troncamento} e il {frazione:.1%} delle prime "
                f"letture della classe {classe.value} e' piu' corto ({corte} su "
                f"{esaminate}), oltre il {massima:.1%} di qc.max_frac_short_reads; "
                f"{len(con_corte)} file ne hanno: {dettaglio}. Le letture piu' corte "
                f"del troncamento vengono scartate, non accorciate",
            )
        )
    if violazioni:
        return violazioni, []

    # Lo scarto opposto (troncare molto sotto la lettura piu' corta) e'
    # E-S1-01, e non si valuta qui: su una stima del minimo sarebbe sbagliato
    # nei due versi. Lo valuta S1, che legge tutte le letture.
    return [], []


def _g10_primer_assente(contesto: Contesto) -> tuple[list[Violazione], list[Avviso]]:
    """Il primer non dev'essere in testa alle letture, e il segnale dev'esserci.

    Due controlli in uno, perche' il primo da solo non basta: un file privo di
    segnale supererebbe la verifica di assenza del primer proprio perche' non
    contiene nulla. Il secondo e' il controllo positivo.

    **Che cosa misura il controllo positivo.** Misura che le letture
    contengano la regione amplificata dichiarata, non che il campione sia
    buono. Il motivo compare anche nei controlli negativi, quanto nei
    biologici: i bianchi amplificano contaminanti, e i contaminanti sono
    batteri con lo stesso 16S.
    Il motivo non distingue quindi il segnale dalla contaminazione, e chi lo
    usasse per quello leggerebbe una risposta a una domanda diversa.

    **Come e' costruito.** Su due scelte, ciascuna con la propria ragione:

    * il primer in testa si giudica file per file, ma sulle stesse classi: in
      un controllo negativo con una manciata di letture una sola che comincia
      per caso come il primer supera qualunque soglia in frazione, e non dice
      nulla di come sono state prodotte le letture;
    * riguarda le sole classi da cui ci si puo' attendere il segnale del
      bersaglio (:data:`~amplicon16s.metadata.models.CLASSI_CONTROLLATE`),
      quindi esclude i controlli negativi. E' un principio, non una
      constatazione: da un bianco non ci si aspetta il materiale in studio,
      anche dove di fatto il motivo ce l'ha per via dei contaminanti;
    * guarda la **mediana** fra quei campioni e non ciascuno di essi: un
      campione biologico legittimo e profondo puo' stare sotto la soglia (per
      esempio se e' dominato da 16S mitocondriale, che non porta il motivo), e
      una verifica file per file lo farebbe fallire.

    Cosi' costruito, il controllo fallisce solo se il segnale manca nel
    complesso - file sbagliati, regione diversa da quella dichiarata - che e'
    la condizione per cui esiste.
    """
    config = contesto.config
    inventario = contesto.inventario
    violazioni: list[Violazione] = []
    avvisi: list[Avviso] = []
    taglio = config.filter.trimLeft

    # Il primer in testa e' un difetto solo se il filtro non lo togliera': con
    # filter.trimLeft maggiore di zero non si cerca (la scansione non lo conta).
    # Si giudica sui campioni biologici e sui controlli positivi, come il
    # motivo: i controlli negativi hanno troppo poche letture perche' una
    # frazione vi abbia significato.
    controllati = {c.accession for c in inventario if c.classe in CLASSI_CONTROLLATE}
    col_primer = {
        accession: s
        for accession, s in contesto.scansione.items()
        if accession in controllati
        and s.valido and s.frazione_primer > config.qc.max_primer_hit_frac
    }
    if col_primer and taglio == 0:
        peggiori = sorted(
            col_primer.items(), key=lambda voce: -voce[1].frazione_primer
        )
        dettaglio = ", ".join(
            f"{s.nome} ({s.frazione_primer:.1%})" for _, s in peggiori[:5]
        )
        violazioni.append(
            Violazione(
                "E-S0-10",
                f"{len(col_primer)} file di campioni biologici o di controlli positivi "
                f"hanno piu' di {config.qc.max_primer_hit_frac:.0%} di letture che iniziano con il "
                f"primer {config.qc.primer_sequence}: {dettaglio}. Imposta "
                f"filter.trimLeft a {len(config.qc.primer_sequence)}, la lunghezza "
                f"del primer",
            )
        )
        # Con il primer in testa il motivo conservato non puo' stare all'inizio
        # della lettura: giudicarlo adesso darebbe una seconda violazione che e'
        # solo la conseguenza della prima. Si giudichera' alla posizione
        # filter.trimLeft, una volta corretto.
        return violazioni, avvisi

    if config.qc.conserved_motif is None:
        avvisi.append(
            Avviso(
                "E-S0-18",
                "qc.conserved_motif e' nullo: verificata la sola assenza del primer in "
                "testa alle letture, non la presenza della regione amplificata",
            )
        )
        return violazioni, avvisi

    con_segnale = [
        contesto.scansione[c.accession].frazione_motivo
        for c in inventario
        if c.classe in CLASSI_CONTROLLATE
        and c.accession in contesto.scansione
        and contesto.scansione[c.accession].valido
    ]
    if not con_segnale:
        violazioni.append(
            Violazione(
                "E-S0-16",
                "nessun campione di una classe da cui attendersi il segnale e' stato "
                "ispezionato: la presenza della regione amplificata non e' verificabile",
            )
        )
        return violazioni, avvisi

    mediana = statistics.median(con_segnale)
    if mediana < config.qc.min_motif_frac:
        dove = "all'inizio delle letture" if taglio == 0 else (
            f"alla posizione {taglio} delle letture (filter.trimLeft)"
        )
        violazioni.append(
            Violazione(
                "E-S0-16",
                f"il motivo conservato {config.qc.conserved_motif!r}, cercato {dove}, "
                f"compare in una mediana del {mediana:.1%} delle letture dei "
                f"{len(con_segnale)} campioni da cui ci si attende il segnale del "
                f"bersaglio, sotto il {config.qc.min_motif_frac:.0%} di "
                f"qc.min_motif_frac: i file potrebbero non contenere la regione "
                f"dichiarata, oppure filter.trimLeft non corrisponde a cio' che precede "
                f"il motivo",
            )
        )
    return violazioni, avvisi


# =========================================================================== #
# Gate sul lotto e sulla composizione delle piastre                           #
# =========================================================================== #


def _g08_lotto_coerente(contesto: Contesto) -> tuple[list[Violazione], list[Avviso]]:
    """Il lotto dichiarato e' leggibile, e ogni piastra e' plausibile.

    La plausibilita' non si misura confrontando una piastra con le altre: una
    composizione inconsueta non e' per questo sbagliata. Si misura su un
    criterio di metodo - una piastra e' plausibile se porta abbastanza
    controlli negativi da sostenere la decontaminazione per piastra, cioe'
    almeno decontam.min_blanks. Le piastre che non lo raggiungono vengono
    segnalate senza bloccare: la decontaminazione su di esse sara' meno
    fondata, e chi legge deve saperlo.
    """
    config = contesto.config
    violazioni: list[Violazione] = []
    avvisi: list[Avviso] = []

    if config.io.batch_table is not None:
        percorso = Path(config.io.batch_table)
        try:
            intestazione = _intestazione(percorso)
        except (OSError, UnicodeDecodeError, csv.Error) as guasto:
            return [
                Violazione("E-S0-08", f"io.batch_table: {percorso} non si apre: {guasto}")
            ], []

        colonne = [
            (dichiarata, nome) for dichiarata, nome in (
                ("decontam.batch_column", config.decontam.batch_column),
                ("err.batch_column", config.err.batch_column),
                ("meta.batch_module_column", config.meta.batch_module_column),
            ) if nome is not None
        ]
        colonne += [("out.batch_columns", c) for c in config.out.batch_columns]
        if tabella_delle_cellule(config) == "lotto":
            colonne.append(("katharoseq.cell_count_column", config.katharoseq.cell_count_column))
        uniche = list(dict((nome, (dichiarata, nome)) for dichiarata, nome in colonne).values())
        violazioni += _colonne_mancanti(
            intestazione, uniche, f"{percorso.name} (io.batch_table)", "E-S0-08"
        )
        if violazioni:
            return violazioni, []
        chiave = [("meta.batch_key_column", config.meta.batch_key_column)]
        violazioni += _valori_spezzati(
            percorso, [c for c in chiave + uniche if c[1] is not None],
            f"{percorso.name} (io.batch_table)", "E-S0-08",
        )
        if violazioni:
            return violazioni, []

        # Un file incompleto o ambiguo: ogni campione deve avervi una e una
        # sola riga, e la corsa se il modello d'errore e' per corsa. Scoprirlo
        # qui costa secondi; in S3 costerebbe il filtro di tutte le letture.
        analisi = contesto.analisi
        if analisi.senza_riga_di_arricchimento:
            violazioni.append(
                Violazione(
                    "E-S0-08",
                    "{} campioni non hanno una riga in {}: {}".format(
                        len(analisi.senza_riga_di_arricchimento), percorso.name,
                        _elenca(sorted(analisi.senza_riga_di_arricchimento)),
                    ),
                )
            )
        if analisi.arricchimento_ambiguo:
            violazioni.append(
                Violazione(
                    "E-S0-08",
                    "{} campioni hanno piu' di una riga in {}: {}".format(
                        len(analisi.arricchimento_ambiguo), percorso.name,
                        _elenca([f"{a} ({n} righe)" for a, n in
                                 sorted(analisi.arricchimento_ambiguo.items())]),
                    ),
                )
            )
        if analisi.senza_piastra:
            violazioni.append(
                Violazione(
                    "E-S0-08",
                    "{} campioni hanno la colonna {!r} (decontam.batch_column) vuota in {}: "
                    "{}. Un campione senza piastra non ha una soglia di profondita' ne' un "
                    "confronto di decontaminazione propri".format(
                        len(analisi.senza_piastra), config.decontam.batch_column, percorso.name,
                        _elenca(sorted(analisi.senza_piastra)),
                    ),
                )
            )
        if analisi.senza_corsa:
            violazioni.append(
                Violazione(
                    "E-S0-08",
                    "{} campioni hanno la colonna {!r} (err.batch_column) vuota in {}: "
                    "{}. Un campione senza corsa non ha un modello d'errore".format(
                        len(analisi.senza_corsa), config.err.batch_column, percorso.name,
                        _elenca(sorted(analisi.senza_corsa)),
                    ),
                )
            )
        if violazioni:
            return violazioni, []

        # Il verso opposto: righe del file che non corrispondono ad alcun
        # campione. Non fermano (un file del lotto puo' coprire piu' assay, o
        # campioni poi esclusi) ma vanno dette: sono anche il sintomo di una
        # chiave scritta in modo diverso da quella dei campioni.
        if analisi.righe_lotto_senza_campione:
            avvisi.append(
                Avviso(
                    "E-S0-19",
                    "{} righe di {} non corrispondono ad alcun campione e sono "
                    "ignorate: {}".format(
                        sum(analisi.righe_lotto_senza_campione.values()), percorso.name,
                        _elenca(sorted(analisi.righe_lotto_senza_campione)),
                    ),
                )
            )

    inventario = contesto.inventario
    per_piastra: dict[str, int] = {}
    for campione in inventario:
        if campione.piastra is None:
            continue
        per_piastra.setdefault(campione.piastra, 0)
        if campione.classe is ClasseCampione.CONTROLLO_NEGATIVO:
            per_piastra[campione.piastra] += 1

    for piastra in sorted(per_piastra, key=lambda p: (len(p), p)):
        negativi = per_piastra[piastra]
        if negativi < config.decontam.min_blanks:
            avvisi.append(
                Avviso(
                    # Codice distinto da quello delle violazioni di G08: questo
                    # avvisa e prosegue, quello ferma. Condividere il codice
                    # significherebbe condividere la categoria di gestione, e
                    # una sola delle due puo' essere giusta.
                    "E-S0-15",
                    f"la piastra {piastra} ha {negativi} controlli negativi, meno dei "
                    f"{config.decontam.min_blanks} di decontam.min_blanks: la "
                    f"decontaminazione su questa piastra sara' meno fondata",
                )
            )
    return violazioni, avvisi


def _g15_coerenza_configurazione(contesto: Contesto) -> tuple[list[Violazione], list[Avviso]]:
    """Adattatore di G15 al contesto: la configurazione e' gia' validata."""
    return _controlla_coerenza(risolvi(contesto.config)), []


def _g11_classi_e_controlli(contesto: Contesto) -> tuple[list[Violazione], list[Avviso]]:
    """G11 nel registro: ogni campione ha una classe, e i controlli si contano.

    Alla verifica delle classi (:func:`_g11_classi_complete`) aggiunge la
    dichiarazione di cio' che manca al dataset: nessun controllo positivo, o
    meno controlli negativi di ``decontam.min_blanks``. Non sono errori (un
    dataset puo' non averne) ma cambiano cio' che le fasi a valle possono fare,
    e vanno detti prima del calcolo.
    """
    violazioni = _g11_classi_complete(contesto.analisi)
    if violazioni:
        return violazioni, []
    conteggi = contesto.inventario.conteggi()
    minimo = contesto.config.decontam.min_blanks
    avvisi: list[Avviso] = []
    negativi = conteggi[ClasseCampione.CONTROLLO_NEGATIVO]
    if negativi < minimo:
        dichiarati = (
            "ctrl.blank_values e' vuoto" if not contesto.config.ctrl.blank_values
            and not contesto.config.ctrl.blank_override_values
            else "le etichette dichiarate ne riconoscono cosi' pochi"
        )
        avvisi.append(
            Avviso(
                "E-S0-17",
                f"il dataset ha {negativi} controlli negativi, meno dei {minimo} di "
                f"decontam.min_blanks ({dichiarati}): i contaminanti non sono "
                f"stimabili dai controlli",
            )
        )
    if conteggi[ClasseCampione.CONTROLLO_POSITIVO] == 0:
        dichiarati = (
            "ctrl.positive_values e' vuoto" if not contesto.config.ctrl.positive_values
            else "nessun campione porta le etichette di ctrl.positive_values"
        )
        avvisi.append(
            Avviso(
                "E-S0-17",
                f"il dataset non ha controlli positivi ({dichiarati}): la soglia di "
                f"profondita' non puo' essere derivata da una curva",
            )
        )
    return [], avvisi


def _adatta(controllo: Any) -> Any:
    """Adatta al contesto i gate sui metadati, che lavorano sull'analisi."""

    def eseguito(contesto: Contesto) -> tuple[list[Violazione], list[Avviso]]:
        return controllo(contesto.analisi), []

    eseguito.__name__ = getattr(controllo, "__name__", "gate")
    eseguito.__doc__ = controllo.__doc__
    return eseguito
