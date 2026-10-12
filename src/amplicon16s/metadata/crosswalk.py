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

from amplicon16s.config.schema import Config
from amplicon16s.io_layer.reads import abbina_file, compagni_per_nome
from amplicon16s.metadata.controls_map import MappaControlli
from amplicon16s.metadata.models import Campione, ClasseCampione, Inventario
from amplicon16s.metadata.tabelle import leggi_tsv, tabella_di_studio

__all__ = ["Analisi", "analizza", "elementi_di_una_cella", "estrai_accession"]


class AccessionNonEstraibile(ValueError):
    """Il nome non contiene esattamente un accession."""


def estrai_accession(testo: str, espressione: re.Pattern[str]) -> str:
    """Estrae l'unica chiave del campione contenuta nel testo.

    La chiave e' la corrispondenza intera dell'espressione oppure, se
    l'espressione ha un gruppo di cattura, il primo gruppo: cosi' un file
    nominato per campione o per corsa si descrive delimitando la chiave con
    cio' che la circonda, per esempio ``^(CAMP[0-9]+)_R1``. Un gruppo che serve
    solo a raggruppare va scritto senza cattura, ``(?:...)``.

    Zero corrispondenze e più di una sono entrambe errori, e per la stessa
    ragione: in nessuno dei due casi si sa a quale campione il file
    appartenga. Prendere la prima corrispondenza sarebbe la scelta comoda e
    silenziosamente sbagliata.
    """
    gruppo = 1 if espressione.groups else 0
    trovati = [c.group(gruppo) for c in espressione.finditer(testo)]
    if not trovati or any(not chiave for chiave in trovati):
        raise AccessionNonEstraibile("nessun accession riconosciuto")
    if len(trovati) > 1:
        raise AccessionNonEstraibile(
            "accession ambiguo, trovati {}: {}".format(len(trovati), ", ".join(trovati))
        )
    return trovati[0]


def elementi_di_una_cella(valore: str) -> list[str]:
    """Gli elementi di una cella che ne elenca piu' d'uno, separati da virgola.

    In ISA-Tab la cella del file di un campione appaiato riporta i due file
    della coppia: ``..._R1.fastq.gz, ..._R2.fastq.gz``.
    """
    return [elemento.strip() for elemento in valore.split(",") if elemento.strip()]


def chiave_dalla_tabella(valore: str, espressione: re.Pattern[str]) -> str:
    """La chiave del campione ricavata dal valore di una colonna dei metadati.

    Il valore puo' essere il nome di un file, e allora la chiave vi si estrae
    come dai nomi dei file; oppure puo' essere gia' la chiave, e allora la si
    usa cosi' com'e' quando l'espressione, pensata per i nomi dei file, non vi
    trova nulla. Una cella che elenca piu' file separati da virgola (i due file
    di una coppia) si legge come elenco: la chiave di ciascun elemento deve
    coincidere, e allora non e' ambigua. Un valore vuoto, o con chiavi diverse,
    resta un errore.
    """
    if not valore:
        raise AccessionNonEstraibile("valore vuoto")
    try:
        return estrai_accession(valore, espressione)
    except AccessionNonEstraibile as guasto:
        if "ambiguo" not in str(guasto):
            return valore
        elementi = elementi_di_una_cella(valore)
        if len(elementi) < 2:
            raise
        chiavi = set()
        for elemento in elementi:
            try:
                chiavi.add(estrai_accession(elemento, espressione))
            except AccessionNonEstraibile:
                raise guasto from None
        if len(chiavi) != 1:
            raise
        return chiavi.pop()


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

    #: Le righe della tabella di assay escluse dalla regola dichiarata in
    #: meta.assay_keep_column: (nome del campione, valore della colonna).
    righe_assay_escluse: list[tuple[str, str]] = field(default_factory=list)

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
    #: Perche' il file non si legge (codifica diversa da UTF-8, struttura non
    #: tabellare), o ``None``. Lo dichiara G08: l'analisi prosegue senza il
    #: file, perche' i gate che la precedono non ne dipendono.
    arricchimento_illeggibile: str | None = None
    arricchimento_senza_accession: list[tuple[str, str]] = field(default_factory=list)
    senza_riga_di_arricchimento: list[str] = field(default_factory=list)
    arricchimento_ambiguo: dict[str, int] = field(default_factory=dict)
    #: Righe del file che non corrispondono ad alcun campione dell'assay, per
    #: chiave: il file ne dichiara piu' di quanti campioni la pipeline tratti.
    righe_lotto_senza_campione: dict[str, int] = field(default_factory=dict)
    #: Campioni con la riga del lotto ma senza il valore di una colonna
    #: dichiarata: la piastra (decontam.batch_column), la corsa (err.batch_column).
    senza_piastra: list[str] = field(default_factory=list)
    senza_corsa: list[str] = field(default_factory=list)

    # --- coppie di file --------------------------------------------------------
    #: In modalita' single-end, gli accession i cui file ripetuti sono le due
    #: letture di una coppia, riconosciute dal marcatore nel nome: accession ->
    #: (marcatori, file).
    coppie_in_single_end: dict[str, tuple[tuple[str, str], list[str]]] = field(
        default_factory=dict)
    #: In modalita' appaiata, i file senza compagno e quelli il cui nome porta
    #: i due marcatori in modo ambiguo.
    file_senza_compagno: tuple[str, ...] = ()
    marcatori_di_file: tuple[str, str] | None = None
    file_ambigui: tuple[str, ...] = ()
    #: In modalita' appaiata, il file reverse di ogni accession.
    reverse_per_accession: dict[str, Path] = field(default_factory=dict)

    # --- risultato ------------------------------------------------------------
    #: File di letture trovati, per accession. Serve ai gate che ispezionano
    #: le sequenze senza doverli cercare una seconda volta.
    file_per_accession: dict[str, Path] = field(default_factory=dict)
    campioni: list[Campione] = field(default_factory=list)

    def inventario(self) -> Inventario:
        """L'inventario dei campioni ricavati dal crosswalk."""
        return Inventario(campioni=tuple(self.campioni))


def _accessioni_dei_file(config: Config, analisi: Analisi) -> dict[str, Path]:
    """Associa a ogni accession il file di letture; registra nell'analisi i nomi senza
    accession e le ripetizioni.
    """
    espressione = re.compile(config.io.accession_regex)
    per_accession: dict[str, list[Path]] = defaultdict(list)

    trovati = sorted(Path(config.io.fastq_dir).glob(config.io.fastq_glob))
    reverse_di: dict[str, Path] = {}
    if config.appaiata:
        # I file si abbinano per nome prima di estrarre la chiave, che si
        # ricava dal solo file forward: il reverse e' il suo compagno.
        per_nome = {p.name: p for p in trovati}
        abbinamento = abbina_file(per_nome, config.pair.file_markers)
        analisi.file_senza_compagno = abbinamento.senza_compagno
        analisi.marcatori_di_file = tuple(config.pair.file_markers)
        analisi.file_ambigui = abbinamento.ambigui
        trovati = [per_nome[nome] for nome in sorted(abbinamento.coppie)]
        reverse_di = {f: per_nome[r] for f, r in abbinamento.coppie.items()}

    for percorso in trovati:
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
    if not config.appaiata:
        for accession, nomi in analisi.accession_ripetuti_nei_file.items():
            marcatori = compagni_per_nome(nomi)
            if marcatori is not None:
                analisi.coppie_in_single_end[accession] = (marcatori, sorted(nomi))
    analisi.reverse_per_accession = {
        accession: reverse_di[percorsi[0].name]
        for accession, percorsi in per_accession.items() if percorsi[0].name in reverse_di
    }
    return {accession: percorsi[0] for accession, percorsi in per_accession.items()}


def _righe_assay(config: Config, analisi: Analisi) -> dict[str, RigaAssay]:
    """Associa a ogni accession la riga della tabella di assay; registra nell'analisi le
    righe senza accession e le ripetizioni.
    """
    espressione = re.compile(config.io.accession_regex)
    per_accession: dict[str, list[RigaAssay]] = defaultdict(list)

    # Le righe da tenere, se la configurazione le dichiara: il confronto e'
    # quello delle etichette dei controlli, senza maiuscole ne' spazi attorno.
    colonna_da_tenere = config.meta.assay_keep_column
    tenuti = (
        None if colonna_da_tenere is None
        else frozenset(v.strip().casefold() for v in config.meta.assay_keep_values or ()))

    for riga in leggi_tsv(Path(config.io.assay_table)):
        nome = riga.get(config.meta.sample_id_column, "")
        riferimento = riga.get(config.meta.accession_column, "")
        if tenuti is not None and (
                riga.get(colonna_da_tenere or "", "").strip().casefold() not in tenuti):
            analisi.righe_assay_escluse.append(
                (nome or riferimento, riga.get(colonna_da_tenere or "", "")))
            continue
        try:
            accession = chiave_dalla_tabella(riferimento, espressione)
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
    """Le righe della tabella di studio raggruppate per nome del campione.

    Senza ``io.study_table`` la tabella di studio e' quella di assay stessa.
    """
    percorso, colonna = tabella_di_studio(config)
    per_nome: dict[str, list[dict[str, str]]] = defaultdict(list)
    for riga in leggi_tsv(percorso):
        per_nome[riga.get(colonna, "")].append(riga)
    return per_nome


#: Valori che le tabelle di metadati usano per dire "non applicabile". Sono
#: convenzioni diffuse, non una proprieta' di questo dataset: trattarli come
#: un valore vero produrrebbe un modulo chiamato "not applicable".
_VALORI_ASSENTI: frozenset[str] = frozenset(
    {"", "na", "n/a", "nan", "none", "null", "not applicable", "notapplicable", "-", "--"}
)


def _valore_o_assente(grezzo: str | None) -> str | None:
    """Il valore ripulito, o ``None`` se è una delle convenzioni per "non applicabile".
    """
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
    try:
        righe = leggi_tsv(Path(percorso))
    except (UnicodeDecodeError, csv.Error) as guasto:
        analisi.arricchimento_illeggibile = str(guasto)
        return per_accession

    colonna = config.meta.batch_key_column
    if colonna is None or (righe and colonna not in righe[0]):
        analisi.arricchimento_senza_colonna = colonna or "<meta.batch_key_column nullo>"
        return per_accession

    for riga in righe:
        grezzo = riga.get(colonna, "")
        try:
            accession = chiave_dalla_tabella(grezzo, espressione)
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

    # G06 : i due insiemi devono coincidere, senza orfani da nessuna parte.
    nei_file = set(file_per_accession)
    nell_assay = set(assay_per_accession)
    analisi.solo_nei_file = tuple(sorted(nei_file - nell_assay))
    analisi.solo_nell_assay = tuple(sorted(nell_assay - nei_file))

    studio_per_nome = _righe_studio(config)
    arricchimento_per_accession = _righe_arricchimento(config, analisi)
    analisi.arricchimento_presente = config.io.batch_table is not None
    analisi.righe_lotto_senza_campione = {
        accession: len(righe)
        for accession, righe in arricchimento_per_accession.items()
        if accession not in nell_assay
    }

    # G03 : la restrizione si misura qui: si itera sulle righe dell'assay e la
    # tabella di studio viene solo consultata. Le sue righe in piu' restano
    # fuori dall'inventario per costruzione, non per un filtro applicato dopo.
    nomi_assay = {r.nome for r in assay_per_accession.values()}
    analisi.righe_studio_estranee = sum(
        len(righe) for nome, righe in studio_per_nome.items() if nome not in nomi_assay
    )

    mappa = MappaControlli.da_configurazione(config.ctrl)
    espressione_modulo = (
        re.compile(config.meta.module_regex)
        if config.meta.derive_module and config.meta.module_regex is not None
        else None
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

        if config.ctrl.column is None:
            # Il dataset dichiara di non avere controlli ne' una colonna che
            # distingua i campioni: sono tutti biologici, senza un materiale.
            materiale, classe = "", ClasseCampione.BIOLOGICO
        else:
            materiale = riga_studio.get(config.ctrl.column, "")
            classe = mappa.classifica(materiale)
        if classe is None:
            non_mappati[materiale].append(nome)
            continue
        # La regola di riclassificazione cambia la classe, non il materiale,
        # che resta il valore originale.
        classe = mappa.riclassifica(classe, riga_studio)

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
        #
        # Senza meta.module_column il dataset non dichiara alcuna posizione, e
        # la regola non ha nulla da giudicare: il modulo, se c'e', e' quello
        # che il file del lotto dichiara. Con la colonna indicata, invece, una
        # posizione vuota resta una non superficie: il dato era atteso e manca.
        superficie = (
            True if config.meta.module_column is None
            else _e_superficie(posizione, non_superfici)
        )

        piastra = corsa = None
        riga_lotto: dict[str, str] | None = None
        if analisi.arricchimento_presente:
            righe = arricchimento_per_accession.get(accession, [])
            if len(righe) == 1:
                riga_lotto = righe[0]
                if config.decontam.batch_column is not None:
                    piastra = _valore_o_assente(riga_lotto.get(config.decontam.batch_column))
                    if piastra is None:
                        analisi.senza_piastra.append(accession)
                if config.err.batch_column is not None:
                    corsa = _valore_o_assente(riga_lotto.get(config.err.batch_column))
                    if corsa is None:
                        analisi.senza_corsa.append(accession)
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
            # Aria, tubi non aperti, posizioni attese e non compilate: restano
            # una categoria a parte, in entrambe le modalita'.
            pass
        elif (
            riga_lotto is not None
            and config.meta.batch_module_column
            and config.meta.batch_module_column in riga_lotto
        ):
            modulo = _valore_o_assente(riga_lotto.get(config.meta.batch_module_column))
        elif espressione_modulo is not None and posizione is not None:
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
                file_reverse=analisi.reverse_per_accession.get(accession),
            )
        )

    analisi.materiali_non_mappati = dict(non_mappati)
    return analisi
