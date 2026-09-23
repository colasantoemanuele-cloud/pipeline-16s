"""Gate di validazione G01–G15.

I gate sono i controlli che precedono l'esecuzione: servono a fermare un'analisi
prima che venga allocato qualunque calcolo, non a metà di un'elaborazione che
dura ore.

Di quindici gate previsti è realizzato il solo **G15 — coerenza interna della
configurazione**. Gli altri non sono ancora implementati.

G15 intercetta una classe di errori che la validazione parametro per parametro
non può cogliere: combinazioni di valori singolarmente validi ma insensati messi
insieme. Una soglia di avviso più alta della soglia di arresto è fatta di due
numeri plausibili e di una coppia che non lo è.

**Rapporto con il catalogo.** I codici ``E-G15-*`` non sono un insieme a parte:
sono voci di :mod:`amplicon16s.errors.catalog` come tutte le altre, e da lì
prendono messaggio operativo e categoria di gestione. Qui resta la sola
informazione che è propria del gate e non avrebbe senso nel catalogo: quali
parametri ciascun controllo sorveglia e dove è implementato.

**Rapporto con lo schema.** Quattro dei sette controlli sono già garantiti dai
vincoli dichiarati in :mod:`amplicon16s.config.schema`. G15 non li riscrive: li
esegue attraverso la validazione dello schema e si limita ad attribuire a
ciascuno il proprio codice di errore. Il registro :data:`CONTROLLI` dice per
ognuno dove è implementato, così la divisione del lavoro è leggibile nel codice
e non solo nelle note.
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
    Config,
    ErroreConfigurazione,
    valida,
)
from amplicon16s.errors.catalog import Categoria, VoceCatalogo, voce
from amplicon16s.errors.exceptions import ErrorePipeline
from amplicon16s.io_layer.reads import StatisticheFile, espandi_iupac, scansiona
from amplicon16s.metadata.crosswalk import Analisi, analizza
from amplicon16s.metadata.models import ClasseCampione, Inventario

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
]


@dataclass(frozen=True)
class Controllo:
    """Un controllo di coerenza, identificato da un codice del catalogo.

    Messaggio e categoria non sono ripetuti qui: vivono nel catalogo degli
    errori, che è la sola fonte di verità dei codici. Questa classe aggiunge
    ciò che il catalogo non sa — su quali parametri il controllo vigila e chi
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
        return voce(self.codice)

    @property
    def descrizione(self) -> str:
        return self.voce.sintesi

    @property
    def categoria(self) -> Categoria:
        return self.voce.categoria


CONTROLLI: Final[tuple[Controllo, ...]] = (
    Controllo(
        "E-G15-01",
        ("filter.minLen", "filter.truncLen"),
        "gate",
    ),
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
        "E-G15-07",
        ("decontam.method", "ctrl.blank_values"),
        "schema",
    ),
    Controllo(
        "E-G15-08",
        ("filter.minLen", "asv.len_min", "asv.len_max", "prev.min_samples"),
        "schema",
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
    (("ctrl", "blank_values"), frozenset({"too_short"}), "E-G15-07"),
)


def _codice_per_problema(problema: str) -> str:
    """Attribuisce un codice di controllo a un problema segnalato dallo schema."""
    if "parametro derivato" in problema:
        return "E-G15-08"
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

    if not controllo.parametri:
        return Violazione(codice, problema)

    return Violazione(
        codice,
        "{} — controllo {} su {}: {}".format(
            problema,
            codice,
            ", ".join(controllo.parametri),
            controllo.descrizione,
        ),
    )


# --------------------------------------------------------------------------- #
# Controlli propri di G15                                                      #
# --------------------------------------------------------------------------- #


def _decontaminazione_attiva(config: Config) -> bool:
    """Se la decontaminazione verrà eseguita.

    Oggi non esiste un parametro che la disattivi, quindi è sempre attiva e la
    condizione è costante. Resta esplicita perché, se un interruttore verrà
    aggiunto, il controllo E-G15-07 debba cambiare in un punto solo.
    """
    del config  # nessun parametro la governa al momento
    return True


def _controlla_coerenza(risolta: ConfigRisolta) -> list[Violazione]:
    """Esegue i controlli che lo schema non copre, sulla configurazione risolta."""
    config = risolta.config
    derivati = risolta.derivati
    violazioni: list[Violazione] = []

    # E-G15-01 — invariante sulla derivazione: minLen discende da truncLen,
    # quindi la disuguaglianza non puo' essere violata da una configurazione.
    # Il controllo resta perche' sorveglia la regola di derivazione, non
    # l'utente: se quella regola cambiasse in modo incoerente, fallirebbe qui.
    if derivati.filter_minLen > config.filter.truncLen:
        violazioni.append(
            Violazione(
                "E-G15-01",
                f"filter.minLen ({derivati.filter_minLen}) supera filter.truncLen "
                f"({config.filter.truncLen}): la lunghezza minima accettata non puo' "
                f"eccedere quella di troncamento",
            )
        )

    # E-G15-02 — lo schema ammette [0, 1]; agli estremi la soglia e' degenere.
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

    # E-G15-03 — l'intervallo delle lunghezze ammesse deve essere sensato.
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

    # E-G15-07 — lo schema impone gia' che l'elenco non sia vuoto; qui resta la
    # sola parte condizionale, che oggi e' sempre vera.
    if _decontaminazione_attiva(config) and not config.ctrl.blank_values:
        violazioni.append(
            Violazione(
                "E-G15-07",
                "la decontaminazione e' attiva (decontam.method="
                f"{config.decontam.method!r}) ma ctrl.blank_values e' vuoto: senza "
                "controlli negativi non e' possibile individuare i contaminanti",
            )
        )

    return violazioni


# --------------------------------------------------------------------------- #
# Esecuzione del gate                                                          #
# --------------------------------------------------------------------------- #


def esegui_g15(dati: Mapping[str, Any]) -> ConfigRisolta:
    """Esegue G15 su una configurazione non ancora validata.

    Restituisce la configurazione risolta nella sua parte statica. Solleva
    :class:`ErroreGate` se un controllo fallisce, con il codice del controllo
    e il dettaglio di quali parametri sono in conflitto.
    """
    try:
        config = valida(dict(dati))
    except ErroreConfigurazione as errore:
        violazioni = [_violazione_da_schema(p) for p in errore.problemi]
        raise ErroreGate("G15", violazioni) from errore

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
                f"{_elenca(file)}",
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

    La restrizione e' garantita per costruzione — si itera sulle righe
    dell'assay — ma resta un modo in cui puo' rompersi: un nome di campione
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
        if self._analisi is None:
            self._analisi = analizza(self.config)
        return self._analisi

    @property
    def inventario(self) -> Inventario:
        return self.analisi.inventario()

    @property
    def scansione_gia_fatta(self) -> bool:
        """Se la scansione e' gia' stata calcolata da qualche gate."""
        return self._scansione is not None

    @property
    def scansione(self) -> dict[str, StatisticheFile]:
        """Statistiche delle prime letture di ogni file, una passata sola."""
        if self._scansione is None:
            self._scansione = scansiona(
                self.analisi.file_per_accession,
                head_reads=self.config.qc.head_reads,
                primer=espandi_iupac(self.config.qc.primer_sequence),
                motivo=self.config.qc.conserved_motif,
            )
        return self._scansione


# =========================================================================== #
# Gate di ingresso e di risorse                                               #
# =========================================================================== #


def _g01_ingressi_leggibili(contesto: Contesto) -> tuple[list[Violazione], list[Avviso]]:
    """Ogni ingresso dichiarato esiste ed e' leggibile."""
    config = contesto.config
    violazioni: list[Violazione] = []

    attesi: list[tuple[str, Path, bool]] = [
        ("io.fastq_dir", Path(config.io.fastq_dir), True),
        ("io.assay_table", Path(config.io.assay_table), False),
        ("io.study_table", Path(config.io.study_table), False),
        ("tax.ref_fasta", Path(config.tax.ref_fasta), False),
    ]
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


def _intestazione(percorso: Path) -> list[str]:
    with open(percorso, encoding="utf-8", newline="") as file:
        lettore = csv.reader(file, delimiter="\t")
        return [c.strip().strip('"') for c in next(lettore, [])]


def _g02_tabelle_apribili(contesto: Contesto) -> tuple[list[Violazione], list[Avviso]]:
    """Le tabelle di metadati si aprono e hanno le colonne dichiarate."""
    config = contesto.config
    violazioni: list[Violazione] = []

    attese = [
        ("io.assay_table", Path(config.io.assay_table),
         {"meta.sample_id_column": config.meta.sample_id_column,
          "meta.accession_column": config.meta.accession_column}),
        ("io.study_table", Path(config.io.study_table),
         {"meta.sample_id_column": config.meta.sample_id_column,
          "ctrl.column": config.ctrl.column,
          "meta.module_column": config.meta.module_column}),
    ]

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
        for dichiarata, nome in colonne.items():
            if nome not in intestazione:
                violazioni.append(
                    Violazione(
                        "E-S0-02",
                        f"{parametro}: la colonna {nome!r}, dichiarata in "
                        f"{dichiarata}, non compare nell'intestazione di "
                        f"{percorso.name}",
                    )
                )
    return violazioni, []


#: Marcatori convenzionali del secondo file di una coppia. Sono convenzioni
#: diffuse, non una proprieta' di questo dataset.
_MARCATORI_INVERSA: Final[tuple[str, ...]] = ("_R2", "_R2.", ".R2.", "_2.fastq", "_2.fq")


def _g07_layout_single_end(contesto: Contesto) -> tuple[list[Violazione], list[Avviso]]:
    """Un solo file per campione e nessun file di lettura inversa."""
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
        if any(marcatore in percorso.name for marcatore in _MARCATORI_INVERSA)
    )
    if inverse:
        violazioni.append(
            Violazione(
                "E-S0-07",
                "{} file sembrano letture inverse: {}. Questa pipeline tratta "
                "solo dati single-end".format(len(inverse), _elenca(inverse)),
            )
        )
    return violazioni, []


def _g12_riferimento_verificato(contesto: Contesto) -> tuple[list[Violazione], list[Avviso]]:
    """Il database tassonomico esiste e il suo checksum e' quello dichiarato."""
    config = contesto.config
    percorso = Path(config.tax.ref_fasta)

    if not percorso.is_file():
        return [Violazione("E-S0-12", f"tax.ref_fasta: {percorso} non esiste")], []

    impronta = hashlib.md5()
    with open(percorso, "rb") as file:
        while blocco := file.read(1024 * 1024):
            impronta.update(blocco)
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
    """Le CPU richieste e lo spazio su disco sono disponibili."""
    config = contesto.config
    violazioni: list[Violazione] = []

    try:
        disponibili = len(os.sched_getaffinity(0))
    except AttributeError:  # piattaforme senza affinita' di processo
        disponibili = os.cpu_count() or 1
    if config.run.threads > disponibili:
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
    """filter.truncLen non puo' superare la lunghezza minima osservata.

    Le letture piu' corte del troncamento vengono scartate, non accorciate: un
    valore troppo alto azzera interi campioni senza che nulla lo segnali.

    La lunghezza minima e' quella delle prime qc.head_reads letture, non di
    tutto il file. E' una stima per eccesso del vero minimo, e il limite e'
    dichiarato nell'avviso che accompagna i casi al margine.
    """
    config = contesto.config
    troncamento = config.filter.truncLen

    corti = {
        accession: s
        for accession, s in contesto.scansione.items()
        if s.valido and s.lunghezza_minima is not None and s.lunghezza_minima < troncamento
    }
    if corti:
        peggiori = sorted(corti.items(), key=lambda voce: voce[1].lunghezza_minima)
        dettaglio = ", ".join(
            f"{s.nome} ({s.lunghezza_minima} bp)" for _, s in peggiori[:5]
        )
        return [
            Violazione(
                "E-S0-09",
                f"filter.truncLen vale {troncamento} ma {len(corti)} file hanno "
                f"letture piu' corte: {dettaglio}. Le letture piu' corte del "
                f"troncamento vengono scartate, non accorciate",
            )
        ], []

    minimi = [
        s.lunghezza_minima
        for s in contesto.scansione.values()
        if s.valido and s.lunghezza_minima is not None
    ]
    if not minimi:
        return [], []

    minimo = min(minimi)
    scarto = minimo - troncamento
    avvisi: list[Avviso] = []
    # Uno scarto ampio non fa perdere campioni, fa perdere basi: si tronca a
    # una lunghezza molto inferiore alla piu' corta delle letture, e ogni
    # lettura cede basi utili senza motivo.
    if scarto > config.filter.truncLen_shortfall_warn:
        avvisi.append(
            Avviso(
                "E-S1-01",
                f"la lettura piu' corta osservata e' di {minimo} bp mentre "
                f"filter.truncLen vale {troncamento}: uno scarto di {scarto} bp, "
                f"oltre i {config.filter.truncLen_shortfall_warn} di "
                f"filter.truncLen_shortfall_warn. Ogni lettura cede {scarto} bp che "
                f"si potrebbero conservare alzando filter.truncLen",
            )
        )
    return [], avvisi


#: Classi di campione da cui ci si puo' attendere il segnale del bersaglio.
#: I controlli negativi ne sono esclusi per principio e non per misura: un
#: bianco non contiene il materiale che si sta studiando, quindi pretendere da
#: lui il segnale sarebbe chiedergli cio' che per definizione non ha. Sul
#: dataset di riferimento i bianchi il motivo ce l'hanno - fra il 79% e il 95%
#: delle letture - ma per via dei contaminanti che amplificano, non del
#: bersaglio: e' la conferma che il motivo misura la presenza della regione
#: amplificata e non la qualita' del campione.
_CLASSI_CON_SEGNALE: Final[tuple[ClasseCampione, ...]] = (
    ClasseCampione.BIOLOGICO,
    ClasseCampione.CONTROLLO_POSITIVO,
)


def _g10_primer_assente(contesto: Contesto) -> tuple[list[Violazione], list[Avviso]]:
    """Il primer non dev'essere in testa alle letture, e il segnale dev'esserci.

    Due controlli in uno, perche' il primo da solo non basta: un file privo di
    segnale supererebbe la verifica di assenza del primer proprio perche' non
    contiene nulla. Il secondo e' il controllo positivo.

    **Che cosa misura il controllo positivo.** Misura che le letture
    contengano la regione amplificata dichiarata, non che il campione sia
    buono. Sul dataset di riferimento il motivo compare nei controlli negativi
    fra il 79% e il 95% delle letture, quanto nei biologici: i bianchi
    amplificano contaminanti, e i contaminanti sono batteri con lo stesso 16S.
    Il motivo non distingue quindi il segnale dalla contaminazione, e chi lo
    usasse per quello leggerebbe una risposta a una domanda diversa.

    **Come e' costruito.** Su due scelte, ciascuna con la propria ragione:

    * riguarda le sole classi da cui ci si puo' attendere il segnale del
      bersaglio, quindi esclude i controlli negativi. E' un principio, non una
      constatazione: da un bianco non ci si aspetta il materiale in studio,
      anche dove di fatto il motivo ce l'ha per via dei contaminanti;
    * guarda la **mediana** fra quei campioni e non ciascuno di essi. Sul
      dataset di riferimento sei campioni biologici stanno sotto la soglia,
      uno al 5,8%, e una verifica file per file li avrebbe fatti fallire pur
      essendo campioni legittimi e profondi.

    Cosi' costruito, il controllo fallisce solo se il segnale manca nel
    complesso - file sbagliati, regione diversa da quella dichiarata - che e'
    la condizione per cui esiste.
    """
    config = contesto.config
    inventario = contesto.inventario
    violazioni: list[Violazione] = []

    col_primer = {
        accession: s
        for accession, s in contesto.scansione.items()
        if s.valido and s.frazione_primer > config.qc.max_primer_hit_frac
    }
    if col_primer:
        peggiori = sorted(
            col_primer.items(), key=lambda voce: -voce[1].frazione_primer
        )
        dettaglio = ", ".join(
            f"{s.nome} ({s.frazione_primer:.1%})" for _, s in peggiori[:5]
        )
        violazioni.append(
            Violazione(
                "E-S0-10",
                f"{len(col_primer)} file hanno piu' di "
                f"{config.qc.max_primer_hit_frac:.0%} di letture che iniziano con il "
                f"primer {config.qc.primer_sequence}: {dettaglio}. Imposta "
                f"filter.trimLeft a {len(config.qc.primer_sequence)}, la lunghezza "
                f"del primer",
            )
        )

    con_segnale = [
        contesto.scansione[c.accession].frazione_motivo
        for c in inventario
        if c.classe in _CLASSI_CON_SEGNALE
        and c.accession in contesto.scansione
        and contesto.scansione[c.accession].valido
    ]
    if not con_segnale:
        violazioni.append(
            Violazione(
                "E-S0-10",
                "nessun campione di una classe da cui attendersi il segnale e' stato "
                "ispezionato: il controllo positivo non e' verificabile",
            )
        )
        return violazioni, []

    mediana = statistics.median(con_segnale)
    if mediana < config.qc.min_motif_frac:
        violazioni.append(
            Violazione(
                "E-S0-10",
                f"il motivo conservato {config.qc.conserved_motif!r} compare in una "
                f"mediana del {mediana:.1%} delle letture dei {len(con_segnale)} "
                f"campioni da cui ci si attende il segnale del bersaglio, sotto il "
                f"{config.qc.min_motif_frac:.0%} di qc.min_motif_frac. L'assenza del "
                f"primer non basta a dire che i file siano corretti: potrebbero non "
                f"contenere la regione dichiarata",
            )
        )
    return violazioni, []


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

        for dichiarata, nome in (
            ("decontam.batch_column", config.decontam.batch_column),
            ("err.batch_column", config.err.batch_column),
        ):
            if nome not in intestazione:
                violazioni.append(
                    Violazione(
                        "E-S0-08",
                        f"la colonna {nome!r}, dichiarata in {dichiarata}, non compare "
                        f"in {percorso.name}: il lotto resterebbe nullo per ogni "
                        f"campione senza che nulla lo segnali",
                    )
                )
        if violazioni:
            return violazioni, []

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


def _adatta(controllo: Any) -> Any:
    """Adatta al contesto i gate sui metadati, che lavorano sull'analisi."""

    def eseguito(contesto: Contesto) -> tuple[list[Violazione], list[Avviso]]:
        return controllo(contesto.analisi), []

    eseguito.__name__ = getattr(controllo, "__name__", "gate")
    eseguito.__doc__ = controllo.__doc__
    return eseguito
