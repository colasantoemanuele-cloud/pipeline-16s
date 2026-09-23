"""Costruzione dell'inventario: file, metadati e classi messi in corrispondenza.

Tre sorgenti, tre ruoli distinti:

* i **file di letture** dicono quali dati esistono davvero su disco;
* la **tabella di assay dell'amplicone** definisce l'insieme dei campioni. È
  lei l'insieme di riferimento, non la tabella campioni di studio: quest'ultima
  è condivisa fra più assay dello stesso studio e contiene righe che a questo
  assay non appartengono. Il join parte dall'assay e la arricchisce, mai il
  contrario;
* la **tabella campioni di studio** porta la colonna della classe. Il join
  verso di essa non è un arricchimento facoltativo: è l'unico modo di sapere
  se un campione è materiale biologico o un controllo.

Un quarto ingresso, il **file di arricchimento del lotto**, è facoltativo.
Quando manca, piastra e corsa restano nulle e la pipeline procede in modalità
a corsa singola: serve a non legare la pipeline a un dataset che per fortuna
dispone di quell'informazione.

Questo modulo **legge e diagnostica, non decide**: raccoglie tutto ciò che non
torna e lo consegna ai gate, che stabiliscono quali problemi fermano
l'esecuzione e con quale codice.
"""

from __future__ import annotations

import csv
import re
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from amplicon16s.config.schema import Config
from amplicon16s.metadata.controls_map import MappaControlli
from amplicon16s.metadata.models import Campione, ClasseCampione, Inventario

__all__ = ["Analisi", "analizza", "estrai_accession"]


class AccessionNonEstraibile(ValueError):
    """Il nome non contiene esattamente un accession."""


def estrai_accession(testo: str, espressione: re.Pattern[str]) -> str:
    """Estrae l'unico accession contenuto nel testo.

    Zero corrispondenze e più di una sono entrambe errori, e per la stessa
    ragione: in nessuno dei due casi si sa a quale campione il file
    appartenga. Prendere la prima corrispondenza sarebbe la scelta comoda e
    silenziosamente sbagliata.
    """
    # finditer e group(0) invece di findall: findall restituisce il contenuto
    # dei gruppi di cattura quando l'espressione ne ha, e l'accession
    # predefinito ne ha uno - "(E|S|D)RX[0-9]{4,}" darebbe "E" al posto
    # dell'accession intero. La forma dell'espressione e' scelta dall'utente,
    # quindi l'estrazione non deve dipendere da come l'ha raggruppata.
    trovati = [corrispondenza.group(0) for corrispondenza in espressione.finditer(testo)]
    if not trovati:
        raise AccessionNonEstraibile("nessun accession riconosciuto")
    if len(trovati) > 1:
        raise AccessionNonEstraibile(
            "accession ambiguo, trovati {}: {}".format(len(trovati), ", ".join(trovati))
        )
    return trovati[0]


def _pulisci(valore: str | None) -> str:
    """Toglie spazi e virgolette: le tabelle ISA le usano in modo incostante."""
    return (valore or "").strip().strip('"').strip()


def _leggi_tsv(percorso: Path) -> list[dict[str, str]]:
    with open(percorso, encoding="utf-8", newline="") as file:
        return [
            {chiave: _pulisci(valore) for chiave, valore in riga.items() if chiave}
            for riga in csv.DictReader(file, delimiter="\t")
        ]


@dataclass(frozen=True)
class RigaAssay:
    """Una riga della tabella di assay, con la chiave estratta."""

    accession: str
    nome: str
    riferimento: str


@dataclass
class Analisi:
    """Esito della lettura delle sorgenti, con tutto ciò che non torna."""

    # --- G04: estrazione dell'accession -------------------------------------
    file_senza_accession: list[tuple[str, str]] = field(default_factory=list)
    assay_senza_accession: list[tuple[str, str]] = field(default_factory=list)

    # --- G05: univocità ------------------------------------------------------
    accession_ripetuti_nei_file: dict[str, list[str]] = field(default_factory=dict)
    accession_ripetuti_nell_assay: dict[str, list[str]] = field(default_factory=dict)

    # --- G06: simmetria ------------------------------------------------------
    solo_nei_file: tuple[str, ...] = ()
    solo_nell_assay: tuple[str, ...] = ()

    # --- G03: restrizione del join ------------------------------------------
    nomi_ambigui_nello_studio: dict[str, int] = field(default_factory=dict)
    nomi_assenti_nello_studio: list[str] = field(default_factory=list)
    righe_studio_estranee: int = 0

    # --- G11: classificazione ------------------------------------------------
    materiali_non_mappati: dict[str, list[str]] = field(default_factory=dict)

    # --- arricchimento del lotto (facoltativo) --------------------------------
    arricchimento_presente: bool = False
    #: Nome della colonna con l'accession, quando il file non ce l'ha. La sua
    #: assenza e' un errore di G04: senza accession il file non e' agganciabile.
    arricchimento_senza_colonna: str | None = None
    arricchimento_senza_accession: list[tuple[str, str]] = field(default_factory=list)
    senza_riga_di_arricchimento: list[str] = field(default_factory=list)
    arricchimento_ambiguo: dict[str, int] = field(default_factory=dict)

    # --- risultato ------------------------------------------------------------
    #: File di letture trovati, per accession. Serve ai gate che ispezionano
    #: le sequenze senza doverli cercare una seconda volta.
    file_per_accession: dict[str, Path] = field(default_factory=dict)
    campioni: list[Campione] = field(default_factory=list)

    def inventario(self) -> Inventario:
        return Inventario(campioni=tuple(self.campioni))


def _accessioni_dei_file(config: Config, analisi: Analisi) -> dict[str, Path]:
    espressione = re.compile(config.io.accession_regex)
    per_accession: dict[str, list[Path]] = defaultdict(list)

    for percorso in sorted(Path(config.io.fastq_dir).glob(config.io.fastq_glob)):
        try:
            accession = estrai_accession(percorso.name, espressione)
        except AccessionNonEstraibile as guasto:
            analisi.file_senza_accession.append((percorso.name, str(guasto)))
            continue
        per_accession[accession].append(percorso)

    analisi.accession_ripetuti_nei_file = {
        accession: [p.name for p in percorsi]
        for accession, percorsi in per_accession.items()
        if len(percorsi) > 1
    }
    return {accession: percorsi[0] for accession, percorsi in per_accession.items()}


def _righe_assay(config: Config, analisi: Analisi) -> dict[str, RigaAssay]:
    espressione = re.compile(config.io.accession_regex)
    per_accession: dict[str, list[RigaAssay]] = defaultdict(list)

    for riga in _leggi_tsv(Path(config.io.assay_table)):
        nome = riga.get(config.meta.sample_id_column, "")
        riferimento = riga.get(config.meta.accession_column, "")
        try:
            accession = estrai_accession(riferimento, espressione)
        except AccessionNonEstraibile as guasto:
            analisi.assay_senza_accession.append(
                (riferimento or f"<riga del campione {nome!r}>", str(guasto))
            )
            continue
        per_accession[accession].append(RigaAssay(accession, nome, riferimento))

    analisi.accession_ripetuti_nell_assay = {
        accession: [r.nome for r in righe]
        for accession, righe in per_accession.items()
        if len(righe) > 1
    }
    return {accession: righe[0] for accession, righe in per_accession.items()}


def _righe_studio(config: Config) -> dict[str, list[dict[str, str]]]:
    per_nome: dict[str, list[dict[str, str]]] = defaultdict(list)
    for riga in _leggi_tsv(Path(config.io.study_table)):
        per_nome[riga.get(config.meta.sample_id_column, "")].append(riga)
    return per_nome


#: Valori che le tabelle di metadati usano per dire "non applicabile". Sono
#: convenzioni diffuse, non una proprieta' di questo dataset: trattarli come
#: un valore vero produrrebbe un modulo chiamato "not applicable".
_VALORI_ASSENTI: frozenset[str] = frozenset(
    {"", "na", "n/a", "nan", "none", "null", "not applicable", "notapplicable", "-", "--"}
)


def _valore_o_assente(grezzo: str | None) -> str | None:
    return None if (grezzo or "").strip().casefold() in _VALORI_ASSENTI else grezzo.strip()


def _e_superficie(posizione: str | None, non_superfici: frozenset[str]) -> bool:
    """Se la posizione dichiarata è una superficie campionabile.

    Una posizione assente non basta a dire che lo sia, quindi conta come non
    superficie: in mancanza del dato non si attribuisce un modulo.
    """
    if posizione is None:
        return False
    return posizione.strip().casefold() not in non_superfici


def _righe_arricchimento(
    config: Config, analisi: Analisi
) -> dict[str, list[dict[str, str]]]:
    """Righe del file facoltativo, indicizzate per accession.

    La chiave e' l'accession come ovunque nella pipeline, non il nome del
    campione: il nome si ripete fra repliche, e un file che identifichi i
    campioni per nome aggancia la replica sbagliata senza che nulla se ne
    accorga.
    """
    per_accession: dict[str, list[dict[str, str]]] = defaultdict(list)
    percorso = config.io.batch_table
    if percorso is None:
        return per_accession

    espressione = re.compile(config.io.accession_regex)
    righe = _leggi_tsv(Path(percorso))

    if righe and config.meta.batch_key_column not in righe[0]:
        analisi.arricchimento_senza_colonna = config.meta.batch_key_column
        return per_accession

    for riga in righe:
        grezzo = riga.get(config.meta.batch_key_column, "")
        try:
            accession = estrai_accession(grezzo, espressione)
        except AccessionNonEstraibile as guasto:
            analisi.arricchimento_senza_accession.append((grezzo, str(guasto)))
            continue
        per_accession[accession].append(riga)
    return per_accession


def analizza(config: Config) -> Analisi:
    """Legge le sorgenti, costruisce i campioni e raccoglie le diagnosi.

    Non solleva eccezioni per i dati incoerenti: le registra. Decidere quali
    fermino l'esecuzione spetta ai gate.
    """
    analisi = Analisi()

    file_per_accession = _accessioni_dei_file(config, analisi)
    analisi.file_per_accession = file_per_accession
    assay_per_accession = _righe_assay(config, analisi)

    # G06 — i due insiemi devono coincidere, senza orfani da nessuna parte.
    nei_file = set(file_per_accession)
    nell_assay = set(assay_per_accession)
    analisi.solo_nei_file = tuple(sorted(nei_file - nell_assay))
    analisi.solo_nell_assay = tuple(sorted(nell_assay - nei_file))

    studio_per_nome = _righe_studio(config)
    arricchimento_per_accession = _righe_arricchimento(config, analisi)
    analisi.arricchimento_presente = config.io.batch_table is not None

    # G03 — la restrizione si misura qui: si itera sulle righe dell'assay e la
    # tabella di studio viene solo consultata. Le sue righe in piu' restano
    # fuori dall'inventario per costruzione, non per un filtro applicato dopo.
    nomi_assay = {r.nome for r in assay_per_accession.values()}
    analisi.righe_studio_estranee = sum(
        len(righe) for nome, righe in studio_per_nome.items() if nome not in nomi_assay
    )

    mappa = MappaControlli.da_configurazione(config.ctrl)
    espressione_modulo = (
        re.compile(config.meta.module_regex) if config.meta.derive_module else None
    )
    non_superfici = frozenset(
        v.strip().casefold() for v in config.meta.non_surface_positions
    )
    non_mappati: dict[str, list[str]] = defaultdict(list)

    for accession in sorted(assay_per_accession):
        riga_assay = assay_per_accession[accession]
        nome = riga_assay.nome

        righe_studio = studio_per_nome.get(nome, [])
        if not righe_studio:
            analisi.nomi_assenti_nello_studio.append(nome)
            continue
        if len(righe_studio) > 1:
            # Piu' righe per lo stesso nome: il join moltiplicherebbe le righe
            # dell'assay, e non e' piu' ristretto.
            analisi.nomi_ambigui_nello_studio[nome] = len(righe_studio)
            continue
        riga_studio = righe_studio[0]

        materiale = riga_studio.get(config.ctrl.column, "")
        classe = mappa.classifica(materiale)
        if classe is None:
            non_mappati[materiale].append(nome)
            continue

        posizione = (
            _valore_o_assente(riga_studio.get(config.meta.module_column, ""))
            if config.meta.module_column
            else None
        )
        # La regola vale in entrambe le vie di attribuzione: un campione che
        # non sta su una superficie non appartiene a un modulo, e non vi va
        # forzato ne' dalla regex ne' dal dato dichiarato. Applicarla a una
        # sola delle due farebbe dipendere il raggruppamento dalla presenza
        # del file di arricchimento, e lo stesso dataset darebbe due risultati.
        superficie = _e_superficie(posizione, non_superfici)

        piastra = corsa = None
        riga_lotto: dict[str, str] | None = None
        if analisi.arricchimento_presente:
            righe = arricchimento_per_accession.get(accession, [])
            if len(righe) == 1:
                riga_lotto = righe[0]
                piastra = _valore_o_assente(riga_lotto.get(config.decontam.batch_column))
                corsa = _valore_o_assente(riga_lotto.get(config.err.batch_column))
            elif not righe:
                analisi.senza_riga_di_arricchimento.append(accession)
            else:
                # Piu' righe per lo stesso accession: quale sia quella giusta
                # non e' deducibile, e attribuirne una a caso metterebbe il
                # campione nel lotto sbagliato. Resta senza.
                analisi.arricchimento_ambiguo[accession] = len(righe)

        # Il modulo viene dal file di arricchimento quando c'e': e' il dato
        # dichiarato, con il nome originale. Derivarlo per espressione
        # regolare e' una ricostruzione, e va usata solo quando il dato
        # dichiarato non esiste - un valore esplicito di "non applicabile" e'
        # un'informazione, non un'assenza da colmare con la regex.
        modulo = None
        if not superficie:
            # Aria, tubi non aperti, posizioni non dichiarate: restano una
            # categoria a parte, in entrambe le modalita'.
            pass
        elif (
            riga_lotto is not None
            and config.meta.batch_module_column
            and config.meta.batch_module_column in riga_lotto
        ):
            modulo = _valore_o_assente(riga_lotto.get(config.meta.batch_module_column))
        elif espressione_modulo is not None:
            trovato = espressione_modulo.match(posizione)
            if trovato is not None:
                # Il gruppo 1 se l'espressione ne ha uno, altrimenti l'intera
                # corrispondenza: anche qui la forma e' scelta dall'utente.
                modulo = trovato.group(1) if trovato.groups() else trovato.group(0)

        analisi.campioni.append(
            Campione(
                accession=accession,
                nome=nome,
                classe=classe,
                materiale=materiale,
                file=file_per_accession.get(accession),
                posizione=posizione,
                modulo=modulo,
                piastra=piastra,
                corsa=corsa,
            )
        )

    analisi.materiali_non_mappati = dict(non_mappati)
    return analisi
