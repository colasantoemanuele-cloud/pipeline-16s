"""Il sottocomando ``init``: dai dati a una configurazione di partenza.

Uso, dalla radice di un clone e nel container::

    python3 scripts/esegui.py init --reads CARTELLA --metadata TABELLA_O_CARTELLA \\
        --output configurazione.yaml [--batch-table FILE] [--reference FILE] \\
        [--out-root CARTELLA]

``--metadata`` e' una tabella dei campioni separata da tabulazioni oppure una
cartella ISA-Tab (una tabella di assay ``a_*.txt`` e una di studio
``s_*.txt``). Il comando legge e non modifica nulla: non scarica, non tocca i
file di ingresso, e non scrive su un file di configurazione che esiste gia'.

Che cosa fa, in ordine: elenca i file di letture; se risultano appaiati (dal
nome o dal contenuto, con le regole di G07 e S1) si ferma senza scrivere e
dice di usare le sole forward; mette in corrispondenza file e tabella; cerca
nelle letture un primer di un elenco di primer 16S comuni; calcola il profilo
di qualita' con lo script di S1 e ne ricava il troncamento che il report
suggerirebbe; propone colonne ed etichette; scrive il file e un riepilogo.

**Letture appaiate.** Se ogni file della cartella ha il suo compagno, con una
delle convenzioni dei nomi (``_R1`` e ``_R2``, ``_1`` e ``_2``), ``init``
rileva la modalita' appaiata e propone anche i parametri del gruppo ``pair``:
il primer e il motivo della lettura reverse come per la forward, la lunghezza
dell'inserto misurata sulle coppie, i due troncamenti dal profilo di qualita'
entro il vincolo che la fusione richiede, l'intervallo delle lunghezze fuse.

**Codici di uscita**: ``0`` file scritto; ``3`` rifiuto (file gia' esistente,
file con le due letture di ogni coppia, ingressi non trovati o non leggibili),
nulla e' stato scritto.
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import statistics
import tempfile
from collections import Counter
from collections.abc import Callable, Mapping
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Any, Final

from amplicon16s.config import defaults
from amplicon16s.config.schema import Config, processori_disponibili
from amplicon16s.errors.exceptions import ErrorePipeline
from amplicon16s.init import rilevamento as r
from amplicon16s.init import vocabolario as v
from amplicon16s.init.scrittura import SENZA_VALORE, Specie, Voce, testo_yaml, valore_yaml
from amplicon16s.io_layer.artifacts import AlberoOutput, Fase, scrivi_atomico
from amplicon16s.io_layer.reads import (
    MARCATORI_DI_FILE,
    RicercaPrimer,
    abbina_file,
    espandi_iupac,
)
from amplicon16s.metadata import ampliconi
from amplicon16s.metadata.crosswalk import elementi_di_una_cella
from amplicon16s.metadata.models import ClasseCampione
from amplicon16s.rbridge.runner import cartella_r, esegui_script
from amplicon16s.soglie.tabella import (
    MARGINE_FRAZIONE,
    SOGLIE,
    non_eseguito,
    valore_esplorativo,
)
from amplicon16s.report.builder import (
    FINESTRA_DI_QUALITA,
    QUALITA_MINIMA_SUGGERITA,
    troncamento_suggerito,
    ultima_posizione_di_qualita,
)
from amplicon16s.steps.s01_profile import (
    NOME_LUNGHEZZE,
    NOME_QUALITA,
    NOME_SCRIPT,
    NOME_VALORI_QUALITA,
    QUALITA_RAGGRUPPATE,
)

__all__ = ["USCITA_RIFIUTO", "Rifiuto", "accetta", "aggiungi_argomenti", "esegui", "prepara"]

USCITA_RIFIUTO: Final = 3

#: La quota della memoria della macchina che la tabella delle sequenze puo'
#: occupare nel valore di partenza di ``qc.max_asv_count``, e i byte di una
#: cella (un intero di R).
_QUOTA_DI_MEMORIA: Final = 0.1
_BYTE_PER_CELLA: Final = 4

#: Il valore di partenza di qc.head_reads: le letture con cui la frazione
#: misurata su un file (primer, motivo, letture corte) ha, nel caso peggiore di
#: una frazione di 0,5, un intervallo di confidenza al 95% di semiampiezza non
#: superiore a PRECISIONE_DELLE_FRAZIONI. Un punto percentuale e' un quinto del
#: margine con cui il sottocomando soglie propone le soglie in frazione:
#: l'errore di stima resta piccolo rispetto alla tolleranza.
QUANTILE_NORMALE_95: Final = 1.96
#: Le basi di sovrapposizione che i troncamenti proposti lasciano oltre la
#: sovrapposizione minima della fusione: con il predefinito di 12 fanno le 20
#: basi che il tutorial di dada2 chiede oltre la variabilita' di lunghezza
#: della regione ("truncLen must be large enough to maintain 20 +
#: biological.length.variation nucleotides of overlap").
MARGINE_DI_SOVRAPPOSIZIONE: Final = 8

#: Di quanto la lunghezza misurata dell'inserto puo' scostarsi da quella attesa
#: dai primer prima che init lo dica: la stessa quota con cui la fusione
#: confronta la mediana delle lunghezze fuse con l'attesa.
SCARTO_DALL_ATTESA_DEI_PRIMER: Final = 0.05

#: La qualita' mediana sotto la quale init non allunga un troncamento, in
#: modalita' appaiata, per rispettare il vincolo della fusione. Il troncamento
#: di una lettura si propone dove la qualita' mediana regge a 30 (un errore
#: atteso ogni mille basi); per fondere le coppie lo si puo' allungare solo
#: dentro posizioni in cui regge almeno a 20 (un errore ogni cento), con la
#: stessa statistica e la stessa finestra. Oltre, la sovrapposizione si
#: reggerebbe su basi di qualita' bassa e il filtro scarterebbe molte letture:
#: i troncamenti restano da scegliere, e il riepilogo lo dice.
QUALITA_MINIMA_PER_FONDERE: Final = 20

#: L'ampiezza, rispetto alla mediana, sotto la quale l'intervallo misurato delle
#: lunghezze dell'inserto non e' proposto come pair.len_min e pair.len_max: le
#: regioni ipervariabili del 16S variano di lunghezza fra i taxa di qualche
#: punto percentuale (per la V4 il tutorial di dada2 tiene le sequenze fra 250
#: e 256 basi, il 2,4% di 253), e un intervallo piu' stretto descrive i taxa di
#: quei campioni, non la regione: scarterebbe le varianti vere appena piu'
#: lunghe o piu' corte.
AMPIEZZA_MINIMA_DELL_INTERVALLO: Final = 0.02

PRECISIONE_DELLE_FRAZIONI: Final = 0.01

#: La quota di letture piu' corte che il troncamento proposto ammette in
#: ciascuna classe: cinque punti percentuali, lo stesso margine con cui il
#: sottocomando soglie propone le soglie in frazione. Con zero basterebbe una
#: lettura accorciata prima del deposito a trascinare la proposta alla sua
#: lunghezza, e con una quota di un punto bastano poche letture su cento; e'
#: una regola della proposta, non una soglia dell'analisi (quella e'
#: qc.max_frac_short_reads, che si sceglie).
QUOTA_DI_LETTURE_CORTE: Final = float(MARGINE_FRAZIONE)

#: Quanto della coda di una lettura il commento del troncamento descrive, in
#: basi, e ogni quante posizioni.
CODA_DESCRITTA: Final = 100
PASSO_DELLA_CODA: Final = 10
LETTURE_PER_UNA_FRAZIONE: Final = math.ceil(
    (QUANTILE_NORMALE_95 * 0.5 / PRECISIONE_DELLE_FRAZIONI) ** 2)

RILEVATO, INDICATO = Specie.RILEVATO, Specie.INDICATO
PROPOSTO, SCELTA = Specie.PROPOSTO, Specie.DA_SCEGLIERE


#: Il prefisso con cui le tabelle ISA-Tab dichiarano le variabili dello studio.
PREFISSO_DEI_FATTORI: Final = "Factor Value["


class Rifiuto(Exception):
    """``init`` non scrive alcuna configurazione, e dice perche'."""


def _intero_positivo(testo: str) -> int:
    """Un argomento intero di almeno uno."""
    try:
        valore = int(testo)
    except ValueError:
        raise argparse.ArgumentTypeError(f"{testo!r} non e' un numero intero") from None
    if valore < 1:
        raise argparse.ArgumentTypeError(f"{valore} non e' ammesso: serve almeno 1")
    return valore


def _sequenza_iupac(testo: str) -> str:
    """Una sequenza nucleotidica in codici IUPAC, in maiuscolo."""
    sequenza = testo.strip().upper()
    try:
        espandi_iupac(sequenza)
    except ValueError as guasto:
        raise argparse.ArgumentTypeError(str(guasto)) from None
    if not sequenza:
        raise argparse.ArgumentTypeError("sequenza vuota")
    return sequenza


def aggiungi_argomenti(parser: argparse.ArgumentParser) -> None:
    """Gli argomenti del sottocomando."""
    parser.add_argument("--reads", required=True, type=Path, metavar="CARTELLA",
                        help="la cartella con i file delle letture, uno per campione")
    parser.add_argument("--metadata", required=True, type=Path, metavar="FILE_O_CARTELLA",
                        help="la tabella dei campioni (TSV) oppure una cartella ISA-Tab")
    parser.add_argument("--assay", metavar="NOME",
                        help="in una cartella ISA-Tab con piu' tabelle di assay, il nome "
                             "del file da usare")
    parser.add_argument("--study-table", type=Path, metavar="FILE",
                        help="la tabella di studio, se --metadata e' una tabella di assay")
    parser.add_argument("--batch-table", type=Path, metavar="FILE",
                        help="il file del lotto (piastra e corsa per campione), se c'e'")
    parser.add_argument("--primer", type=_sequenza_iupac, metavar="SEQUENZA",
                        help="il primer forward usato per l'amplificazione, in codici IUPAC, "
                             "se non e' fra quelli che init riconosce da se'")
    parser.add_argument("--primer-reverse", type=_sequenza_iupac, metavar="SEQUENZA",
                        help="il primer reverse usato per l'amplificazione, in codici IUPAC "
                             "e scritto da 5' a 3', se non e' fra quelli che init riconosce da "
                             "se': in modalita' appaiata si cerca in testa alle letture reverse, "
                             "in single-end in coda alle letture piu' lunghe della regione")
    parser.add_argument("--reference", type=Path, metavar="FILE",
                        help="il riferimento tassonomico, nel formato di dada2")
    parser.add_argument("--out-root", type=Path, metavar="CARTELLA",
                        help="la cartella di uscita prevista per l'esecuzione")
    parser.add_argument("--output", required=True, type=Path, metavar="FILE",
                        help="il file di configurazione da scrivere: non deve esistere")
    parser.add_argument("--head-reads", type=_intero_positivo, metavar="N",
                        help="ispeziona le sole prime N letture di ogni file per le coppie "
                             "e il primer (senza: tutte)")
    parser.add_argument("--skip-profile", action="store_true",
                        help="non calcola il profilo di qualita': nessuna proposta di "
                             "troncamento ne' di funzione d'errore")
    parser.add_argument("--threads", type=_intero_positivo, metavar="N",
                        help="processi per l'ispezione (senza: quelli utilizzabili)")
    parser.add_argument("--accetta-proposte", action="store_true",
                        help="scrive come righe attive le proposte, i valori di partenza e, "
                             "per le soglie di qualita' senza partenza, il valore che non "
                             "ferma: restano da scrivere le sole scelte senza valore")


# --------------------------------------------------------------------------- #
# Ingressi                                                                     #
# --------------------------------------------------------------------------- #


def _tabelle(args: argparse.Namespace) -> tuple[Path, Path | None]:
    """La tabella di assay e quella di studio, dagli argomenti."""
    metadati: Path = args.metadata
    if not metadati.exists():
        raise Rifiuto(f"--metadata: {metadati} non esiste")
    if metadati.is_file():
        if args.study_table is not None and not args.study_table.is_file():
            raise Rifiuto(f"--study-table: {args.study_table} non esiste")
        return metadati, args.study_table
    assay = sorted(metadati.glob("a_*.txt"))
    if args.assay is not None:
        assay = [p for p in assay if p.name == args.assay]
        if not assay:
            raise Rifiuto(f"--assay: {args.assay} non e' fra le tabelle di assay di {metadati}")
    if not assay:
        raise Rifiuto(
            f"{metadati} non contiene una tabella di assay ISA-Tab (a_*.txt): indica con "
            "--metadata la tabella dei campioni")
    if len(assay) > 1:
        raise Rifiuto(
            f"{metadati} contiene {len(assay)} tabelle di assay: indica con --assay quella "
            "del sequenziamento 16S, fra " + ", ".join(repr(p.name) for p in assay))
    studio = args.study_table
    if studio is None:
        trovate = sorted(metadati.glob("s_*.txt"))
        if len(trovate) > 1:
            raise Rifiuto(
                f"{metadati} contiene {len(trovate)} tabelle di studio: indica quella da "
                "usare con --study-table")
        studio = trovate[0] if trovate else None
    return assay[0], studio


def _memoria_della_macchina() -> int | None:
    """I byte di memoria della macchina, o ``None`` se non si possono leggere."""
    try:
        return os.sysconf("SC_PHYS_PAGES") * os.sysconf("SC_PAGE_SIZE")
    except (ValueError, OSError, AttributeError):
        return None


def _due_cifre(valore: float) -> int:
    """Il valore arrotondato per difetto a due cifre significative."""
    intero = int(valore)
    if intero < 100:
        return max(intero, 1)
    scala = 10 ** (len(str(intero)) - 2)
    return intero // scala * scala


# --------------------------------------------------------------------------- #
# Profilo                                                                      #
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class TaglioCercato:
    """Il taglio dei primer cercati che ``init`` propone: i primer trovati
    lungo le letture (forward, e reverse in modalita' appaiata) e la regola di
    ricerca di partenza.
    """

    forward: r.PrimerSfasato
    reverse: r.PrimerSfasato | None

    @property
    def discordanze(self) -> int:
        """Le discordanze di partenza, sul primer piu' corto."""
        return r.discordanze_di_partenza(
            *(p.primer.sequenza for p in (self.forward, self.reverse) if p is not None))

    @property
    def sfasamento(self) -> int:
        """Lo sfasamento massimo proposto: il piu' grande osservato."""
        return max(p.sfasamento_massimo for p in (self.forward, self.reverse) if p is not None)

    def ricerca(self, reverse: bool = False) -> RicercaPrimer:
        """La ricerca del primer di una delle due letture, con i valori proposti."""
        primer = self.reverse if reverse else self.forward
        assert primer is not None
        return RicercaPrimer(primer.primer.sequenza, self.discordanze, self.sfasamento)

    def per_lo_script(self, reverse: bool = False) -> dict[str, Any]:
        """Il parametro ``taglio_primer`` di ``R/01_profile.R``."""
        ricerca = self.ricerca(reverse)
        return {"primer": ricerca.primer, "discordanze_massime": ricerca.discordanze_massime,
                "sfasamento_massimo": ricerca.sfasamento_massimo, "scarta": True}


def _descrivi_gli_sfasamenti(primer: r.PrimerSfasato) -> str:
    """La distribuzione degli sfasamenti di un primer, per il commento."""
    totale = sum(primer.sfasamenti.values())
    quote = ", ".join(
        f"{sfasamento}: {letture / totale:.0%}"
        for sfasamento, letture in primer.sfasamenti.items()
        if letture >= r.QUOTA_DI_UNO_SFASAMENTO * totale)
    return (f"basi che precedono il primer {primer.primer.nome} e quota delle letture che lo "
            f"portano ({quote})")


def _profilo(
    file_per_chiave: dict[str, Path], negativi: set[str], processi: int,
    taglio: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Il profilo delle letture con lo script di S1, e cio' che se ne ricava.

    Lo script e' quello della fase, eseguito in una cartella temporanea: le
    tabelle sono le stesse che S1 scrivera', e il troncamento suggerito e'
    quello che il report calcolerebbe da esse. La frazione ammessa di letture
    corte e' zero: nessuna tolleranza e' ancora stata scelta.
    """
    with tempfile.TemporaryDirectory(prefix="amplicon16s-init-") as temporanea:
        albero = AlberoOutput(temporanea)
        esegui_script(
            cartella_r() / NOME_SCRIPT,
            # Percorsi assoluti: lo script lavora in una cartella temporanea, e
            # una cartella delle letture indicata in forma relativa non vi si
            # troverebbe.
            {"campioni": {c: str(Path(p).resolve()) for c, p in file_per_chiave.items()},
             "processi": processi,
             # Con il primer cercato il profilo e' quello delle letture tagliate.
             **({"taglio_primer": taglio} if taglio is not None else {})},
            albero, Fase.QC_PROFILES, passo="S1",
        )
        cartella = albero.cartella(Fase.QC_PROFILES)

        def righe(nome: str) -> list[dict[str, str]]:
            with open(cartella / nome, encoding="utf-8", newline="") as file:
                return list(csv.DictReader(file, delimiter="\t"))

        lunghezze, qualita = righe(NOME_LUNGHEZZE), righe(NOME_QUALITA)
        valori = [int(x["qualita"]) for x in righe(NOME_VALORI_QUALITA)]
    classi = {
        c: (ClasseCampione.CONTROLLO_NEGATIVO if c in negativi else ClasseCampione.BIOLOGICO).value
        for c in file_per_chiave
    }
    suggerito = troncamento_suggerito(lunghezze, qualita, classi, QUOTA_DI_LETTURE_CORTE)
    distribuzione: Counter[int] = Counter()
    for riga in lunghezze:
        if riga["campione"] not in negativi:
            distribuzione[int(riga["lunghezza"])] += int(float(riga["letture"]))
    totale = sum(distribuzione.values())
    moda = max(distribuzione, key=lambda x: (distribuzione[x], x)) if distribuzione else None
    # La qualita' per posizione dei campioni non di controllo negativo, per il
    # commento: la mediana, fra i campioni, della media e del primo quartile.
    per_posizione: dict[int, tuple[list[float], list[float]]] = {}
    for riga in qualita:
        if riga["campione"] not in negativi:
            medie, quartili = per_posizione.setdefault(int(riga["posizione"]), ([], []))
            medie.append(float(riga["media"]))
            quartili.append(float(riga["p25"]))
    # La qualita' mediana per posizione, come la legge la regola del
    # troncamento: la mediana, fra i campioni, della mediana di ciascuno.
    mediane: dict[int, list[float]] = {}
    for riga in qualita:
        if riga["campione"] not in negativi:
            mediane.setdefault(int(riga["posizione"]), []).append(float(riga["mediana"]))
    mediana = {posizione: statistics.median(x) for posizione, x in mediane.items()}
    return {
        "suggerito": suggerito,
        "mediana": mediana,
        "per_fondere": ultima_posizione_di_qualita(mediana, QUALITA_MINIMA_PER_FONDERE),
        "valori_di_qualita": sorted(valori),
        "minima": min(distribuzione) if distribuzione else None,
        "massima": max(distribuzione) if distribuzione else None,
        "moda": moda,
        "distribuzione": dict(distribuzione),
        "totale": totale,
        "qualita": {
            posizione: (statistics.median(medie), statistics.median(quartili))
            for posizione, (medie, quartili) in per_posizione.items()},
    }


def _percento(frazione: float, decimali: int = 2) -> str:
    """Una frazione come percentuale, con l'articolo che la precede: «lo
    0.08%», «l'1%», «l'81.9%», «il 5%».
    """
    numero = f"{frazione:.{decimali}%}"
    intero = numero.split(".")[0].rstrip("%")
    if intero == "0":
        return f"lo {numero}"
    if intero in ("1", "11") or intero.startswith("8"):
        return f"l'{numero}"
    return f"il {numero}"


def _quota_piu_corte(profilo: Mapping[str, Any], lunghezza: int) -> float:
    """La quota di letture del profilo piu' corte di ``lunghezza``."""
    totale = profilo["totale"]
    return sum(
        n for osservata, n in profilo["distribuzione"].items() if osservata < lunghezza
    ) / totale if totale else 0.0


def _descrivi_le_lunghezze(profilo: Mapping[str, Any], troncamento: int | None) -> str:
    """Le lunghezze osservate, per il commento di un troncamento: gli estremi,
    la lunghezza piu' frequente con la sua quota, e la quota di letture piu'
    corte del troncamento proposto e di 10 e 20 basi sotto la piu' frequente.
    """
    moda = profilo["moda"]
    quota_della_moda = profilo["distribuzione"][moda] / profilo["totale"]
    corte = [
        f"di {lunghezza} ({che_cosa}) {_percento(_quota_piu_corte(profilo, lunghezza))}"
        for lunghezza, che_cosa in (
            (troncamento, "il troncamento proposto"),
            (moda - 10, "10 basi sotto la piu' frequente"),
            (moda - 20, "20 basi sotto la piu' frequente"),
        ) if lunghezza is not None and lunghezza > 0
    ]
    return (
        f"lunghezze osservate: da {profilo['minima']} a {profilo['massima']}; la piu' "
        f"frequente e' {moda}, con {_percento(quota_della_moda, 1)} delle letture; letture piu' corte "
        + ", ".join(corte))


def _descrivi_la_qualita(profilo: Mapping[str, Any]) -> str:
    """La qualita' in coda alle letture, per il commento di un troncamento:
    media e primo quartile ogni :data:`PASSO_DELLA_CODA` posizioni nelle ultime
    :data:`CODA_DESCRITTA`, e all'ultima.
    """
    qualita = profilo["qualita"]
    ultima = max(qualita)
    posizioni = [
        p for p in sorted(qualita)
        if p > ultima - CODA_DESCRITTA and (p % PASSO_DELLA_CODA == 0 or p == ultima)]
    return (
        "qualita' per posizione in coda (media e primo quartile; fra i campioni, la "
        "mediana): " + ", ".join(
            f"{p}: {qualita[p][0]:.0f} e {qualita[p][1]:.0f}" for p in posizioni))


def _regola_del_troncamento(suggerito: Mapping[str, Any]) -> str:
    """La regola con cui e' proposto un troncamento, con i suoi due termini."""
    return (
        f"il minore fra la lunghezza ammessa ({suggerito['per_lunghezza']}: la piu' lunga che "
        f"in ogni classe scarta non oltre {_percento(QUOTA_DI_LETTURE_CORTE, 0)} delle letture) e "
        f"l'ultima posizione fino a cui regge la qualita' mediana ({suggerito['per_qualita']}: "
        f"almeno {QUALITA_MINIMA_SUGGERITA} nella mediana di {FINESTRA_DI_QUALITA} posizioni "
        "consecutive)")


def _perdita(suggerito: Mapping[str, Any]) -> float:
    """La quota di letture piu' corte del troncamento suggerito, nella classe
    che ne perde di piu'.
    """
    return max(suggerito["perdita_per_classe"].values(), default=0.0)


def _voce_delle_letture_corte(metti: Callable[..., None], perdita: float, quali: str) -> None:
    """La voce di ``qc.max_frac_short_reads`` coerente con i troncamenti
    proposti: zero se non scartano letture, altrimenti da fissare sopra la
    quota che scartano.
    """
    if perdita == 0:
        metti("qc.max_frac_short_reads", SCELTA, 0.0,
              f"partenza: zero, perche' {quali} non scarta letture; va alzata solo se si "
              "sceglie un troncamento piu' lungo")
        return
    metti("qc.max_frac_short_reads", SCELTA, SENZA_VALORE,
          f"soglia di qualita' da fissare sui propri dati: {quali} scarta {_percento(perdita)} "
          "delle letture, nella classe che ne perde di piu', e la soglia deve stare sopra "
          "quella quota (altrimenti validate si ferma: E-S0-09)",
          "misura da guardare: le lunghezze delle letture, in 02_qc_profiles/lunghezze.tsv; "
          "dopo la prima esecuzione la propone il sottocomando soglie")


# --------------------------------------------------------------------------- #
# Primer cercato lungo la lettura                                              #
# --------------------------------------------------------------------------- #

def _voci_del_taglio_cercato(
    metti: Callable[..., None], note: list[str], cercato: TaglioCercato, appaiata: bool
) -> None:
    """Le voci del taglio dei primer cercati, quando le letture portano il
    primer preceduto da un numero variabile di basi: il modo, il primer
    forward, il taglio fisso a zero e la regola di ricerca di partenza.
    """
    PROPOSTO, SCELTA = Specie.PROPOSTO, Specie.DA_SCEGLIERE
    primer = cercato.forward
    entrambi = [p for p in (cercato.forward, cercato.reverse) if p is not None]
    variabile = max(entrambi, key=lambda p: p.fuori_posizione)
    metti("filter.primer_removal", PROPOSTO, defaults.PRIMER_CERCATO,
          f"il primer non sta a una posizione fissa: il {variabile.fuori_posizione:.0%} delle "
          f"letture che portano il primer {variabile.primer.nome} lo ha a uno sfasamento "
          f"diverso dal piu' frequente, oltre il {r.QUOTA_FUORI_POSIZIONE:.0%} che distingue "
          "un protocollo con distanziatori dagli errori di sequenza",
          *(_descrivi_gli_sfasamenti(p) for p in entrambi),
          "un taglio a posizione fissa (filter.trimLeft) lascerebbe le letture sfasate, e "
          "lo stesso taxon diviso in varianti diverse")
    metti("qc.primer_sequence", PROPOSTO, primer.primer.sequenza,
          f"primer {primer.primer.nome} ({primer.primer.fonte}), trovato lungo le letture: "
          f"nella mediana dei file il {primer.mediana:.0%} delle letture lo porta entro "
          f"{r.SFASAMENTO_CERCATO} basi dall'inizio (dal {primer.minima:.0%} al "
          f"{primer.massima:.0%}), su non piu' di {r.LETTURE_PER_IL_MOTIVO} letture per file")
    metti("filter.trimLeft", PROPOSTO, 0,
          "zero: con filter.primer_removal search il primer lo toglie la ricerca, e "
          "filter.truncLen e' la lunghezza della lettura dopo il taglio")
    metti("filter.primer_max_mismatch", SCELTA, cercato.discordanze,
          f"partenza: un decimo della lunghezza del primer piu' corto, la frazione di errori "
          f"predefinita di cutadapt ({r.FRAZIONE_DI_ERRORI}), arrotondata per difetto",
          "si contano le sole sostituzioni; una N della lettura e' una discordanza")
    metti("filter.primer_max_offset", PROPOSTO, cercato.sfasamento,
          "il piu' grande numero di basi che precedono il primer, fra quelli che portano "
          f"almeno il {r.QUOTA_DI_UNO_SFASAMENTO:.0%} delle letture con il primer",
          "se il protocollo dichiara distanziatori piu' lunghi, scrivi la loro lunghezza massima")
    metti("filter.primer_untrimmed", SCELTA, defaults.SENZA_PRIMER_SCARTATA,
          "partenza: una lettura senza il primer si scarta"
          + (", con la sua compagna" if appaiata else "")
          + " (come cutadapt con --discard-untrimmed): in una libreria che porta il primer "
          "e' fuori bersaglio o danneggiata, e tenuta intera comincerebbe in un punto "
          "diverso dalle altre")
    metti("qc.max_frac_untrimmed", SCELTA, SENZA_VALORE,
          "soglia di qualita' da fissare sui propri dati; misura da guardare: la quota di "
          "letture senza il primer, in 03_filtered/taglio_primer.tsv")
    note.append(
        f"Il primer {variabile.primer.nome} e' preceduto da un numero variabile di basi: "
        "proposto il taglio cercato lungo la lettura (filter.primer_removal: \"search\").")


# --------------------------------------------------------------------------- #
# Modalita' appaiata                                                           #
# --------------------------------------------------------------------------- #


def _voci_delle_coppie(
    args: argparse.Namespace,
    metti: Callable[..., None],
    voci: dict[str, Voce],
    note: list[str],
    non_indovinati: list[str],
    avanzamento: Callable[[str], None],
    *,
    file: list[Path],
    nomi: list[str],
    reverse_di: dict[str, Path],
    chiave: r.Chiave | None,
    negativi: set[str],
    ispezioni: dict[str, r.Ispezione],
    cercati: tuple[v.Primer, ...],
    inserti: list[list[int | None]],
    quante: str,
    processi: int,
    tabelle: list[Any],
    taglio_forward: int,
    primer_forward: v.Primer | None,
    profilo_forward: dict[str, Any] | None,
    cercato: TaglioCercato | None = None,
) -> None:
    """Le voci del gruppo ``pair`` per un dataset appaiato, e i troncamenti
    delle due letture entro il vincolo della fusione.

    Segue per la lettura reverse le stesse regole della forward (primer
    dall'elenco, motivo dal consenso, troncamento dal profilo), e in piu'
    misura sulle coppie la lunghezza dell'inserto: da quella vengono il vincolo
    sui troncamenti e l'intervallo delle lunghezze fuse.
    """
    PROPOSTO, RILEVATO, SCELTA, INDICATO = (
        Specie.PROPOSTO, Specie.RILEVATO, Specie.DA_SCEGLIERE, Specie.INDICATO)
    da_considerare = {
        n: i for n, i in ispezioni.items()
        if chiave is None or _chiave_del_file(chiave, n) not in negativi
    } or dict(ispezioni)

    # ---- Primer reverse e taglio ---------------------------------------------
    osservato, parziali = r.primer_nelle_letture(ispezioni, cercati)
    dai_metadati: list[str] = []
    con_etichetta_forward: list[str] = []
    _, nei_metadati, _ = r.primer_dai_metadati(tabelle, dai_metadati, con_etichetta_forward)
    # Candidate a primer reverse: le sequenze che i metadati dichiarano tali, e
    # quelle che coincidono con un primer reverse dell'elenco.
    candidati_reverse = list(dict.fromkeys(
        [*dai_metadati,
         *(s for s in nei_metadati if s not in con_etichetta_forward
           and r.primer_dell_elenco(s, v.PRIMER_16S_REVERSE) is not None)]))
    primer_noto: str | None = None
    taglio = 0
    # Il primer reverse riconosciuto nelle letture, in testa o lungo la lettura.
    primer_reverse = osservato.primer if osservato is not None else None
    if cercato is not None:
        assert cercato.reverse is not None
        primer_reverse = cercato.reverse.primer
        primer_noto = primer_reverse.sequenza
        metti("pair.primer_sequence_reverse", PROPOSTO, primer_noto,
              f"primer {primer_reverse.nome} ({primer_reverse.fonte}), trovato lungo le "
              f"letture reverse: nella mediana dei file il {cercato.reverse.mediana:.0%} delle "
              "letture lo porta", _descrivi_gli_sfasamenti(cercato.reverse))
        metti("pair.trimLeft_reverse", PROPOSTO, 0,
              "zero: con filter.primer_removal search il primer reverse lo toglie la ricerca")
    elif osservato is not None:
        trovato = osservato.primer
        primer_noto, taglio = trovato.sequenza, len(trovato.sequenza)
        if osservato.in_ogni_file and args.head_reads is None:
            specie, limite = RILEVATO, ""
            dove = (f"in testa alle letture reverse di ogni file: dal {osservato.minima:.0%} al "
                    f"{osservato.massima:.0%} delle letture, su {quante}")
        elif not osservato.in_ogni_file:
            sotto = osservato.sotto_soglia
            specie = PROPOSTO
            dove = (f"in testa alle letture reverse nella mediana dei file "
                    f"({osservato.mediana:.0%}), ma non in tutti: dal {osservato.minima:.0%} al "
                    f"{osservato.massima:.0%}")
            limite = (
                f"attenzione: {len(sotto)} file su {len(ispezioni)} non lo portano in testa "
                f"({', '.join(sotto[:5])}{'...' if len(sotto) > 5 else ''}): un taglio fisso "
                "toglierebbe loro basi vere. Vanno capiti prima di accettare")
            note.append(
                f"Il primer reverse {trovato.nome} e' in testa alle letture reverse della "
                f"maggior parte dei file ma non di {len(sotto)}: il primer e' proposto, il taglio "
                "resta da scegliere.")
        else:
            specie = PROPOSTO
            dove = f"in testa alle letture reverse di ogni file, su {quante}"
            limite = ("verificato sulle sole prime letture (--head-reads): senza l'opzione "
                      "init lo verifica su tutte")
        metti("pair.primer_sequence_reverse", specie, trovato.sequenza,
              f"primer {trovato.nome} ({trovato.fonte}), {dove}", limite)
        if osservato.in_ogni_file:
            metti("pair.trimLeft_reverse", specie, taglio,
                  f"la lunghezza del primer {trovato.nome}, che le letture reverse portano in "
                  "testa", limite)
        else:
            metti("pair.trimLeft_reverse", SCELTA, SENZA_VALORE,
                  f"non proposto: il taglio sarebbe {taglio}, la lunghezza del primer "
                  f"{trovato.nome}, ma non tutti i file reverse lo portano in testa", limite,
                  f"scrivi {taglio} dopo aver capito quei file; motivo delle reverse, "
                  "lunghezze fuse e troncamenti di questo file presuppongono quel taglio")
            non_indovinati.append(
                f"pair.trimLeft_reverse: il primer reverse {trovato.nome} non e' in testa alle "
                "letture reverse di ogni file.")
    else:
        variabile = "; ".join(f"{p.nome} nel {f:.1%}" for p, f in parziali)
        metti("pair.trimLeft_reverse", PROPOSTO, 0,
              f"nessun primer reverse dell'elenco e' in testa alle letture reverse (su "
              f"{quante}); un primer fuori elenco va indicato a mano, con la sua lunghezza qui",
              "primer cercati: " + ", ".join(p.nome for p in cercati),
              (f"attenzione: una parte delle letture reverse inizia con un primer "
               f"({variabile}): e' in posizione variabile, e un taglio fisso non lo toglie")
              if variabile else "")
        if args.primer_reverse is not None:
            primer_noto = args.primer_reverse
            metti("pair.primer_sequence_reverse", INDICATO, args.primer_reverse,
                  "--primer-reverse; non e' in testa alle letture reverse")
        elif len(candidati_reverse) == 1:
            primer_noto = candidati_reverse[0]
            dell_elenco = r.primer_dell_elenco(primer_noto, v.PRIMER_16S_REVERSE)
            metti("pair.primer_sequence_reverse", PROPOSTO, primer_noto,
                  "sequenza che i metadati dichiarano come primer reverse; non e' nelle "
                  "letture reverse" if primer_noto in dai_metadati else
                  f"sequenza dei metadati che coincide con il primer reverse {dell_elenco.nome} "
                  f"dell'elenco ({dell_elenco.fonte}): il verso e' quello del primer, anche se "
                  "i metadati non lo dicono; non e' nelle letture reverse")
        else:
            nominati = r.primer_nominati(tabelle, v.PRIMER_16S_REVERSE)
            metti("pair.primer_sequence_reverse", SCELTA, None,
                  "partenza: nullo, perche' il primer reverse non e' nelle letture reverse e "
                  "i metadati non ne riportano una sola sequenza. Con il primer nullo la "
                  "validazione non lo cerca. Se e' noto, rilancia init con --primer-reverse "
                  "SEQUENZA: init propone allora anche il motivo delle letture reverse",
                  ("piu' sequenze dei metadati sono candidate a primer reverse (dichiarate "
                   "tali, o coincidenti con un primer reverse dell'elenco): "
                   + ", ".join(candidati_reverse)) if candidati_reverse else "",
                  ("i metadati nominano un primer reverse senza riportarne la sequenza; le "
                   "forme che init conosce con quel nome: "
                   + "; ".join(f"{p.nome} {p.sequenza} ({p.fonte})" for p in nominati))
                  if nominati else "")
            non_indovinati.append(
                "pair.primer_sequence_reverse: il primer reverse non e' nelle letture reverse, "
                "e i metadati non ne riportano una sola sequenza.")

    # ---- Motivo della reverse --------------------------------------------------
    motivo = None
    dissentono: list[str] = []
    if primer_noto is not None:
        del_consenso = da_considerare
        if cercato is not None and cercato.reverse is not None:
            del_consenso = {
                n: i.dopo_il_primer(cercato.reverse.indice) for n, i in da_considerare.items()}
        dissentono = r.file_che_dissentono(del_consenso, taglio, len(primer_noto))
        motivo = r.motivo_dal_consenso(
            [i for n, i in del_consenso.items() if n not in dissentono],
            taglio, len(primer_noto))
    if dissentono:
        elenco = ", ".join(reverse_di[n].name for n in dissentono[:5]) + (
            "..." if len(dissentono) > 5 else "")
        note.append(
            f"{len(dissentono)} file reverse non iniziano come gli altri e sono esclusi dal "
            f"consenso del motivo delle letture reverse ({elenco}): da guardare (una corsa con "
            "i due file scambiati, un'altra regione, sequenze fuori bersaglio).")
    if motivo is not None:
        metti("pair.conserved_motif_reverse", PROPOSTO, motivo,
              f"il consenso delle {len(primer_noto or '')} basi "
              + ("che seguono il primer reverse" if taglio or cercato is not None else
                 "con cui iniziano le letture reverse")
              + f", su non piu' di {r.LETTURE_PER_IL_MOTIVO} letture per file: una base o "
              f"una coppia di basi dove copre almeno il {r.SOGLIA_BASE_DEL_MOTIVO:.0%} delle "
              "letture",
              "serve a riconoscere una corsa con i due file scambiati e una libreria con i "
              "due versi mescolati",
              (f"attenzione: {len(dissentono)} file reverse non iniziano come gli altri e "
               "sono esclusi dal consenso: da guardare prima di eseguire") if dissentono else "")
    else:
        metti("pair.conserved_motif_reverse", SCELTA, None,
              "partenza: nullo. "
              + ("Il consenso delle letture reverse non da' un motivo abbastanza definito"
                 if primer_noto else
                 "Senza conoscere il primer reverse init non propone il motivo con cui "
                 "iniziano le letture reverse")
              + ". Con motivo e primer reverse nulli l'orientamento delle letture reverse "
              "non si verifica")

    # ---- Lunghezza dell'inserto ---------------------------------------------------
    inserto = r.stima_dell_inserto(inserti)
    dei_primer = None
    if primer_forward is not None and primer_reverse is not None:
        if primer_forward.posizione is not None and primer_reverse.posizione is not None:
            dei_primer = primer_reverse.posizione - primer_forward.posizione + 1
            if cercato is not None:
                # L'inserto e' misurato senza i primer: lo e' anche l'atteso.
                dei_primer -= len(primer_forward.sequenza) + len(primer_reverse.sequenza)
    # Con meno di meta' delle coppie sovrapposte la misura vale solo se la
    # lunghezza attesa dai primer riconosciuti cade nell'intervallo misurato:
    # in letture lunghe le code rovinate nascondono la sovrapposizione di una
    # parte delle coppie, e il riscontro dei primer dice che le altre sono il
    # bersaglio.
    con_riscontro = False
    if inserto is None and dei_primer is not None:
        parziale = r.stima_dell_inserto(inserti, r.QUOTA_DI_COPPIE_CON_IL_RISCONTRO_DEI_PRIMER)
        if parziale is not None and (
                parziale.minima * (1 - SCARTO_DALL_ATTESA_DEI_PRIMER) <= dei_primer
                <= parziale.massima * (1 + SCARTO_DALL_ATTESA_DEI_PRIMER)):
            inserto, con_riscontro = parziale, True
    tagli = taglio_forward + taglio
    campionate = sum(len(x) for x in inserti)
    sovrapposte = sum(1 for x in inserti for y in x if y is not None)
    #: La lunghezza piu' grande da coprire con le due letture, nella stessa
    #: forma dei troncamenti (con i primer, se le letture li portano).
    da_coprire: int | None = None
    if inserto is not None:
        da_coprire = inserto.massima
        confronto = ""
        if dei_primer is not None:
            concordano = abs(inserto.mediana - dei_primer) <= (
                SCARTO_DALL_ATTESA_DEI_PRIMER * dei_primer)
            confronto = (
                f"dai primer riconosciuti nelle letture ({primer_forward.nome} e "
                f"{primer_reverse.nome}) su Escherichia coli se ne attendono {dei_primer}: "
                + ("concordano" if concordano else
                   "la mediana misurata se ne scosta di oltre il "
                   f"{SCARTO_DALL_ATTESA_DEI_PRIMER:.0%}, ma l'atteso e' nell'intervallo "
                   "misurato: la regione varia fra i taxa" if con_riscontro else
                   "NON concordano, e vale la lunghezza misurata sulle letture"))
            if not concordano and not con_riscontro:
                note.append(
                    f"La lunghezza dell'inserto misurata sulle letture (mediana "
                    f"{inserto.mediana}) non concorda con quella attesa dai primer riconosciuti "
                    f"({dei_primer}): vale la misurata. Verifica la regione amplificata.")
        misura = (
            f"misurata su {inserto.sovrapposte} coppie con le letture sovrapposte, su "
            f"{inserto.coppie} campionate (le prime {r.COPPIE_PER_L_INSERTO} di ogni "
            f"campione): mediana {inserto.mediana}, dal {inserto.minima} al {inserto.massima} "
            "(percentili 0,5 e 99,5)"
            + ("; comprende i primer, che le letture portano" if tagli else "")
            + ("; e' la lunghezza senza i primer e senza le basi che li precedono, tolte a "
               "ogni coppia come le toglie il taglio" if cercato is not None else ""))
        poche = (
            f"attenzione: la sovrapposizione si conferma in {inserto.sovrapposte} coppie su "
            f"{inserto.coppie}, meno di meta'. Nelle altre le code delle letture non hanno "
            f"almeno {r.SOVRAPPOSIZIONE_MINIMA_STIMATA} basi di qualita' almeno "
            f"{r.QUALITA_DI_UNA_BASE_AFFIDABILE} da confrontare: la misura vale perche' la "
            "lunghezza attesa dai primer riconosciuti cade nell'intervallo misurato, ma la "
            "quota di coppie fuse potra' essere bassa") if con_riscontro else ""
        if con_riscontro:
            note.append(
                f"La lunghezza dell'inserto e' misurata su {inserto.sovrapposte} coppie su "
                f"{inserto.coppie}: nelle altre le code delle letture sono troppo rovinate "
                "per confermare la sovrapposizione. Dopo la prima esecuzione guarda la quota "
                "di coppie fuse.")
        metti("pair.expected_merged_len", PROPOSTO, inserto.mediana - tagli,
              "la mediana della lunghezza dell'inserto, tolte le basi dei due tagli", misura,
              confronto, poche)
        stretto = (inserto.massima - inserto.minima
                   < AMPIEZZA_MINIMA_DELL_INTERVALLO * inserto.mediana)
        if stretto:
            # Un intervallo cosi' stretto descrive i taxa di questi campioni,
            # non la regione: non e' una proposta sicura.
            fusa = inserto.mediana - tagli
            largo = (max(1, math.floor(fusa * (1 - SCARTO_DALL_ATTESA_DEI_PRIMER))),
                     math.ceil(fusa * (1 + SCARTO_DALL_ATTESA_DEI_PRIMER)))
            perche = (
                f"l'intervallo misurato delle lunghezze, da {max(1, inserto.minima - tagli)} a "
                f"{max(1, inserto.massima - tagli)} (percentili 0,5 e 99,5, tolte le basi dei "
                f"due tagli), e' largo meno del {AMPIEZZA_MINIMA_DELL_INTERVALLO:.0%} della "
                "mediana: descrive i taxa di questi campioni, non la regione amplificata, che "
                "fra i taxa varia di piu'")
            metti("pair.len_min", SCELTA, SENZA_VALORE,
                  f"non proposto: {perche}", misura,
                  f"raccomandato: un intervallo piu' largo attorno alla mediana, per esempio da "
                  f"{largo[0]} a {largo[1]} (la mediana meno e piu' il "
                  f"{SCARTO_DALL_ATTESA_DEI_PRIMER:.0%}), da stringere dopo la prima esecuzione "
                  "guardando 05_pair_merge/lunghezze_fuse.tsv; una sequenza fusa fuori "
                  "dall'intervallo esce dalla tabella")
            metti("pair.len_max", SCELTA, SENZA_VALORE, "si sceglie insieme a pair.len_min")
            note.append(
                "Intervallo delle lunghezze fuse non proposto: " + perche + ". pair.len_min e "
                f"pair.len_max restano da scegliere (raccomandato: da {largo[0]} a {largo[1]}, "
                "da stringere dopo la prima esecuzione).")
            non_indovinati.append(
                "pair.len_min e pair.len_max: l'intervallo misurato e' troppo stretto per "
                "descrivere la regione.")
        else:
            metti("pair.len_min", PROPOSTO, max(1, inserto.minima - tagli),
                  "il percentile 0,5 della lunghezza dell'inserto, tolte le basi dei due tagli",
                  "una sequenza fusa piu' corta esce dalla tabella: allarga l'intervallo se la "
                  "regione varia di piu' fra i taxa attesi",
                  "l'intervallo descrive le coppie misurate, non la regione: se molto sotto la "
                  "lunghezza attesa c'e' una popolazione a se' (amplificazioni fuori bersaglio), "
                  "il percentile basso la comprende. Dopo la prima esecuzione guarda 05_pair_merge/lunghezze_fuse.tsv "
                  "e stringi l'intervallo attorno alla regione")
            metti("pair.len_max", PROPOSTO, max(1, inserto.massima - tagli),
                  "il percentile 99,5 della lunghezza dell'inserto, tolte le basi dei due tagli")
    elif dei_primer is not None:
        # La sovrapposizione non si misura, ma i due primer sono riconosciuti:
        # l'intervallo si propone dalla loro distanza, dichiarandolo.
        minima = math.floor(dei_primer * (1 - SCARTO_DALL_ATTESA_DEI_PRIMER))
        massima = math.ceil(dei_primer * (1 + SCARTO_DALL_ATTESA_DEI_PRIMER))
        da_coprire = massima
        origine = (
            f"dalla distanza fra i primer riconosciuti nelle letture ({primer_forward.nome} e "
            f"{primer_reverse.nome}): su Escherichia coli {dei_primer} basi"
            + ("" if cercato is not None else ", primer compresi")
            + f". Non e' una misura: solo {sovrapposte} coppie su {campionate} campionate "
            "hanno le due letture sovrapposte in modo verificabile")
        metti("pair.expected_merged_len", PROPOSTO, dei_primer - tagli,
              f"la lunghezza attesa {origine}", "tolte le basi dei due tagli")
        metti("pair.len_min", PROPOSTO, max(1, minima - tagli),
              f"la lunghezza attesa meno il {SCARTO_DALL_ATTESA_DEI_PRIMER:.0%}, {origine}",
              "la regione varia fra i taxa: dopo la prima esecuzione guarda "
              "05_pair_merge/lunghezze_fuse.tsv e correggi l'intervallo attorno alla regione")
        metti("pair.len_max", PROPOSTO, max(1, massima - tagli),
              f"la lunghezza attesa piu' il {SCARTO_DALL_ATTESA_DEI_PRIMER:.0%}, dalla "
              "distanza fra i primer riconosciuti")
        note.append(
            "La lunghezza dell'inserto non si misura dalle letture (solo "
            f"{sovrapposte} coppie su {campionate} con le letture sovrapposte): lunghezza "
            "attesa e intervallo delle sequenze fuse sono proposti dalla distanza fra i primer "
            f"riconosciuti ({primer_forward.nome} e {primer_reverse.nome}), da verificare dopo "
            "la prima esecuzione sulla quota di coppie fuse e su lunghezze_fuse.tsv.")
    else:
        perche = (
            f"la lunghezza dell'inserto non si stima dalle letture: solo {sovrapposte} coppie "
            f"su {campionate} campionate hanno le due letture sovrapposte")
        metti("pair.expected_merged_len", SCELTA, None,
              f"partenza: nullo, perche' {perche}. Senza una lunghezza attesa il confronto "
              "della mediana delle sequenze fuse non si fa")
        metti("pair.len_min", SCELTA, SENZA_VALORE,
              f"{perche}: l'intervallo delle lunghezze fuse va dichiarato conoscendo la "
              "regione amplificata")
        metti("pair.len_max", SCELTA, SENZA_VALORE, "si sceglie insieme a pair.len_min")
        note.append(
            perche[0].upper() + perche[1:] + ". Le letture non coprono l'intera regione "
            "amplificata, e le coppie non si potrebbero fondere: come alternativa, la "
            "modalita' single-end con le sole letture forward (io.layout: \"single-end\", "
            "senza il gruppo pair, e i soli file forward nella cartella).")
        non_indovinati.append(
            "pair.len_min, pair.len_max e i due troncamenti: la lunghezza dell'inserto non si "
            "stima dalle letture.")

    # ---- Profilo della reverse e troncamenti ------------------------------------
    profilo: dict[str, Any] | None = None
    if not args.skip_profile and profilo_forward is not None:
        avanzamento("profilo di qualita' delle letture reverse...")
        leggibili = {
            (_chiave_del_file(chiave, n) if chiave else n): reverse_di[n]
            for n in nomi if ispezioni[n].statistiche.valido
        }
        try:
            profilo = _profilo(
                leggibili, negativi, processi,
                cercato.per_lo_script(reverse=True) if cercato is not None else None)
        except ErrorePipeline as guasto:
            note.append(
                f"Profilo di qualita' delle letture reverse non calcolato ({guasto.codice}: "
                f"{guasto.dettaglio or guasto.voce.sintesi}).")
    sovrapposizione = defaults.PAIR_MIN_OVERLAP + MARGINE_DI_SOVRAPPOSIZIONE
    vincolo = (
        "vincolo della fusione: (filter.truncLen - filter.trimLeft) + (pair.truncLen_reverse "
        "- pair.trimLeft_reverse) deve superare la lunghezza della regione senza i primer di "
        f"almeno pair.min_overlap piu' {MARGINE_DI_SOVRAPPOSIZIONE} basi di margine")
    if cercato is not None:
        vincolo = (
            "vincolo della fusione: con il primer cercato i troncamenti sono lunghezze dopo il "
            "taglio, e filter.truncLen + pair.truncLen_reverse deve superare la lunghezza "
            "dell'inserto senza i primer di almeno pair.min_overlap piu' "
            f"{MARGINE_DI_SOVRAPPOSIZIONE} basi di margine")
    forward = profilo_forward["suggerito"] if profilo_forward else None
    della_reverse = profilo["suggerito"] if profilo else None
    proposti: tuple[int, int] | None = None
    senza_qualita = False
    if (da_coprire is not None and forward and della_reverse
            and forward["per_lunghezza"] is not None
            and della_reverse["per_lunghezza"] is not None):
        assert profilo_forward is not None and profilo is not None
        necessarie = da_coprire + sovrapposizione
        # Mai oltre la lunghezza ammessa: quella che in ogni classe scarta non
        # piu' della quota di letture corte della proposta.
        massimo_forward = int(forward["per_lunghezza"])
        massimo_reverse = int(della_reverse["per_lunghezza"])
        # Il troncamento di ciascuna lettura dalla sua qualita'; zero se la
        # qualita' non regge nemmeno all'inizio, e allora decide il vincolo.
        t_forward = min(int(forward["per_qualita"] or 0), massimo_forward)
        t_reverse = min(int(della_reverse["per_qualita"] or 0), massimo_reverse)
        ammessa_forward, ammessa_reverse = massimo_forward, massimo_reverse
        limite_forward, limite_reverse = t_forward, t_reverse
        allungati = False
        if t_forward + t_reverse < necessarie:
            # Prima la forward, che di norma ha la qualita' migliore in coda;
            # poi la reverse; mai oltre la lunghezza ammessa (che la gran parte
            # delle letture raggiunge), ne' oltre l'ultima posizione in cui la
            # qualita' mediana regge alla soglia per fondere: una
            # sovrapposizione fatta di basi di qualita' bassa non e' una
            # proposta.
            allungati = True
            massimo_forward = min(
                massimo_forward, max(limite_forward, int(profilo_forward["per_fondere"] or 0)))
            massimo_reverse = min(
                massimo_reverse, max(limite_reverse, int(profilo["per_fondere"] or 0)))
            t_forward = min(massimo_forward, max(t_forward, necessarie - t_reverse))
            t_reverse = min(massimo_reverse, max(t_reverse, necessarie - t_forward))
        if t_forward + t_reverse >= necessarie:
            proposti = (t_forward, t_reverse)
            conto = (
                f"{t_forward} + {t_reverse} = {t_forward + t_reverse} basi, e ne servono "
                f"{necessarie}: l'inserto piu' lungo "
                + ("misurato" if inserto is not None else "atteso dai primer")
                + f" ({da_coprire}) piu' {sovrapposizione} di sovrapposizione"
                + ("" if allungati or t_forward + t_reverse == necessarie else
                   ". Sono piu' basi di quelle che servono: accorciare un troncamento fin "
                   "dove la sua qualita' e' migliore (la riga della qualita' in coda, qui "
                   "sotto) riduce le coppie che il filtro scarta, e la fusione resta "
                   f"possibile finche' la somma dei due non scende sotto {necessarie}"))
            def allungato(del_profilo: Mapping[str, Any], limite: int, troncamento: int) -> str:
                """Di quanto un troncamento supera il limite di qualita' della
                sua lettura, e la qualita' mediana piu' bassa che raggiunge."""
                if troncamento <= limite:
                    return ""
                piu_bassa = min(
                    del_profilo["mediana"][p] for p in range(limite + 1, troncamento + 1)
                    if p in del_profilo["mediana"])
                return (
                    f"; allungato di {troncamento - limite} basi oltre il calo della qualita' "
                    "quanto serve a fondere le coppie, dentro posizioni in cui la qualita' "
                    f"mediana e' almeno {QUALITA_MINIMA_PER_FONDERE} (la piu' bassa: "
                    f"{piu_bassa:.0f})")

            for parametro, quale, troncamento, esito, del_profilo, in_testa, limite in (
                ("filter.truncLen", "forward", t_forward, forward, profilo_forward,
                 "filter.trimLeft", limite_forward),
                ("pair.truncLen_reverse", "reverse", t_reverse, della_reverse, profilo,
                 "pair.trimLeft_reverse", limite_reverse),
            ):
                perdita = _quota_piu_corte(del_profilo, troncamento)
                metti(parametro, PROPOSTO, troncamento,
                      f"troncamento della lettura {quale}: {_regola_del_troncamento(esito)}"
                      + allungato(del_profilo, limite, troncamento)
                      + ("; non scarta alcuna lettura" if perdita == 0 else
                                     f"; scarta {_percento(perdita)} delle letture {quale}"),
                      conto, vincolo,
                      _descrivi_le_lunghezze(del_profilo, troncamento),
                      _descrivi_la_qualita(del_profilo),
                      "e' la lunghezza dopo il taglio del primer" if cercato is not None else
                      f"comprende le basi di {in_testa}")
            for parametro, troncamento, esito in (
                ("filter.truncLen_shortfall_warn", t_forward, forward),
                ("pair.truncLen_shortfall_warn_reverse", t_reverse, della_reverse),
            ):
                metti(parametro, SCELTA,
                      min(max(int(esito["senza_perdite"]) - troncamento, 0), troncamento - 1),
                      "partenza: lo scarto fra la lettura piu' corta e il troncamento "
                      "proposto; con un valore minore il profilo segnalerebbe il troncamento "
                      "proposto (E-S1-01)")
            _voce_delle_letture_corte(
                metti,
                max(_quota_piu_corte(profilo_forward, t_forward),
                    _quota_piu_corte(profilo, t_reverse)),
                "la coppia di troncamenti proposti")
        elif ammessa_forward + ammessa_reverse >= necessarie:
            senza_qualita = True
            note.append(
                "Le coppie non si fondono senza basi di qualita' bassa: la qualita' mediana "
                f"regge fino alla posizione {limite_forward} nelle letture forward (lunghezza "
                f"piu' frequente {profilo_forward['moda']}) e fino alla {limite_reverse} nelle "
                f"reverse (lunghezza piu' frequente {profilo['moda']}); anche allungando ciascun "
                "troncamento fin dove la qualita' mediana e' almeno "
                f"{QUALITA_MINIMA_PER_FONDERE} le due letture danno {massimo_forward} + "
                f"{massimo_reverse} = "
                f"{massimo_forward + massimo_reverse} basi, e ne servono {necessarie} (l'inserto "
                f"piu' lungo, {da_coprire}, piu' {sovrapposizione} di sovrapposizione). I due "
                "troncamenti restano da scegliere. Si consiglia di usare le sole letture "
                "forward: prepara una cartella con i soli file forward (GUIDA.md, passo 2) e "
                "rilancia init con --reads su quella cartella e un altro --output.")
        else:
            note.append(
                f"Nessuna coppia di troncamenti entro la lunghezza ammessa permette la fusione: "
                f"le due letture danno al piu' {ammessa_forward} + {ammessa_reverse} = "
                f"{ammessa_forward + ammessa_reverse} basi, e ne servono {necessarie} "
                f"(l'inserto piu' lungo, {da_coprire}, piu' {sovrapposizione} di "
                "sovrapposizione). I due troncamenti restano da scegliere; come alternativa, "
                "la modalita' single-end con le sole letture forward (io.layout: "
                "\"single-end\", senza il gruppo pair, e i soli file forward nella cartella).")
    if proposti is None:
        perche = (
            "la lunghezza dell'inserto non e' stimata" if da_coprire is None else
            "il profilo di qualita' non e' stato calcolato" if not (forward and della_reverse)
            else "le coppie non si fondono senza basi di qualita' bassa (mediana sotto "
                 f"{QUALITA_MINIMA_PER_FONDERE}): si consiglia di usare le sole letture forward"
                 if senza_qualita
            else "nessuna coppia di troncamenti entro la lunghezza ammessa rispetta il vincolo")
        for parametro in ("filter.truncLen", "pair.truncLen_reverse"):
            metti(parametro, SCELTA, SENZA_VALORE,
                  f"da scegliere insieme all'altro troncamento: {perche}", vincolo,
                  *((_descrivi_le_lunghezze(del_profilo, None), _descrivi_la_qualita(del_profilo))
                    if (del_profilo := (profilo_forward if parametro == "filter.truncLen"
                                        else profilo)) is not None and del_profilo["totale"]
                    and del_profilo["qualita"] else ()))
        for parametro in ("filter.truncLen_shortfall_warn",
                          "pair.truncLen_shortfall_warn_reverse"):
            metti(parametro, SCELTA, SENZA_VALORE, "si sceglie insieme al suo troncamento")
        metti("qc.max_frac_short_reads", SCELTA, SENZA_VALORE,
              "si sceglie insieme ai troncamenti")

    # ---- Soglie della modalita' appaiata -------------------------------------------
    for parametro, misura in (
        ("pair.max_primer_hit_frac_reverse",
         "frazione di letture reverse con il primer reverse in testa, in "
         "01_input_validation/letture_ispezionate_reverse.tsv"),
        ("pair.min_merged_frac_warn", "quota di coppie fuse, in 05_pair_merge/fusione.json"),
        ("pair.min_merged_frac_stop", "quota di coppie fuse, in 05_pair_merge/fusione.json"),
    ):
        metti(parametro, SCELTA, SENZA_VALORE,
              f"soglia di qualita' da fissare sui propri dati; misura da guardare: {misura}")
    metti("pair.min_merged_frac_sample_warn", SCELTA, SENZA_VALORE,
          "soglia di qualita' facoltativa, per singolo campione: lasciata in commento non "
          "da' alcun controllo; misura da guardare: quota di coppie fuse di ogni campione, "
          "in 05_pair_merge/fusione.tsv")
    mancanti = [c for c in defaults.PARAMETRI_DELLE_COPPIE if c not in voci]
    assert not mancanti, mancanti


# --------------------------------------------------------------------------- #
# Dalle osservazioni alle voci                                                 #
# --------------------------------------------------------------------------- #


def prepara(args: argparse.Namespace, avanzamento: Callable[[str], None] = lambda _: None) -> tuple[
    dict[str, Voce], list[str], list[str], list[str]
]:
    """Esamina i dati e restituisce le voci della configurazione, le righe
    d'intestazione, cio' che ``init`` non stabilisce e le note del riepilogo.

    Solleva :class:`Rifiuto` se non c'e' una configurazione da proporre.
    """
    voci: dict[str, Voce] = {}
    note: list[str] = []
    non_indovinati: list[str] = []

    def metti(chiave: str, specie: Specie, valore: Any, origine: str, *altre: str) -> None:
        voci[chiave] = Voce(chiave, specie, valore, origine, tuple(a for a in altre if a))

    # ---- Letture ---------------------------------------------------------
    cartella: Path = args.reads
    if not cartella.is_dir():
        raise Rifiuto(f"--reads: {cartella} non e' una cartella")
    per_estensione = r.file_di_letture(cartella)
    if not per_estensione:
        raise Rifiuto(
            f"{cartella} non contiene file di letture ({', '.join(v.ESTENSIONI_FASTQ)})")
    file = [p for lista in per_estensione.values() for p in lista]
    nomi = [p.name for p in file]
    # Se ogni file ha il suo compagno il dataset e' appaiato: da qui in poi
    # ``file`` e ``nomi`` sono i soli file forward, e i reverse li seguono.
    appaiamento = r.coppie_di_file(nomi)
    reverse_di: dict[str, Path] = {}
    if appaiamento is not None:
        marcatori, coppie_di_file = appaiamento
        per_nome = {p.name: p for p in file}
        reverse_di = {f: per_nome[x] for f, x in coppie_di_file.items()}
        file = [per_nome[f] for f in sorted(coppie_di_file)]
        nomi = [p.name for p in file]
        for estensione in list(per_estensione):
            per_estensione[estensione] = [
                p for p in per_estensione[estensione] if p.name in coppie_di_file]
    else:
        inverse = r.secondi_di_coppia(nomi)
        if inverse:
            raise Rifiuto(_senza_compagno(nomi, inverse))

    avanzamento(f"letture: {len(file)} file in {cartella}"
                + (f", ciascuno con il suo file reverse ({marcatori[0]} e {marcatori[1]})"
                   if reverse_di else ""))
    # Il primer indicato da chi lancia il comando si cerca nelle letture come
    # quelli dell'elenco, e viene per primo.
    cercati = v.PRIMER_16S
    if args.primer is not None:
        cercati = (v.Primer("indicato", args.primer, "--primer"),
                  *(p for p in v.PRIMER_16S if p.sequenza != args.primer))
    primer = tuple(espandi_iupac(p.sequenza) for p in cercati)
    processi = args.threads or processori_disponibili()
    # I primer reverse: in modalita' appaiata si cercano in testa alle letture
    # reverse; in single-end, il loro complemento inverso in coda alle letture.
    cercati_reverse = v.PRIMER_16S_REVERSE
    if args.primer_reverse is not None:
        cercati_reverse = (
            v.Primer("indicato", args.primer_reverse, "--primer-reverse"),
            *(p for p in v.PRIMER_16S_REVERSE if p.sequenza != args.primer_reverse))
    with ProcessPoolExecutor(max_workers=processi) as gruppo:
        ispezioni = dict(zip(nomi, gruppo.map(
            r.ispeziona,
            [(str(p), args.head_reads, primer, tuple(x.sequenza for x in cercati),
              () if reverse_di else tuple(x.sequenza for x in cercati_reverse))
             for p in file])))
        # In modalita' appaiata: le letture reverse, con i primer reverse, e la
        # lunghezza dell'inserto su un campione di coppie di ogni campione.
        ispezioni_reverse: dict[str, r.Ispezione] = {}
        inserti: list[list[int | None]] = []
        if reverse_di:
            primer_reverse = tuple(espandi_iupac(p.sequenza) for p in cercati_reverse)
            ispezioni_reverse = dict(zip(nomi, gruppo.map(
                r.ispeziona,
                [(str(reverse_di[n]), args.head_reads, primer_reverse,
                  tuple(x.sequenza for x in cercati_reverse)) for n in nomi])))
        # Il primer preceduto da un numero variabile di basi: se lo e' in una
        # delle due letture il taglio cercato vale per entrambe, purche' la
        # ricerca trovi il primer di ciascuna.
        cercato: TaglioCercato | None = None
        sfasato = r.primer_sfasato(ispezioni, cercati)
        sfasato_reverse = r.primer_sfasato(ispezioni_reverse, cercati_reverse)
        if sfasato is not None and (not reverse_di or sfasato_reverse is not None) and (
                sfasato.da_cercare or (sfasato_reverse is not None and sfasato_reverse.da_cercare)):
            cercato = TaglioCercato(sfasato, sfasato_reverse if reverse_di else None)
        if reverse_di:
            # Con il taglio cercato l'inserto si misura senza i primer.
            ricerche = () if cercato is None else (cercato.ricerca(), cercato.ricerca(True))
            inserti = list(gruppo.map(
                r.inserti_di_una_coppia_di_file,
                [(str(p), str(reverse_di[p.name]), *ricerche) for p in file]))
    if reverse_di:
        # Un file con i segni di entrambe le letture contiene le coppie intere;
        # il marcatore dell'altra lettura dice una corsa con i file scambiati,
        # che init segnala e la validazione ferma.
        intercalati, invertite = [], []
        for n in nomi:
            for ispezione, attesa in ((ispezioni[n], "1"), (ispezioni_reverse[n], "2")):
                esito = (ispezione.statistiche.segni_di_coppia_inattesi(attesa)
                         if ispezione.statistiche.valido else None)
                if esito is not None:
                    (intercalati if esito[0] == "intercalato" else invertite).append(
                        f"{ispezione.statistiche.nome} ({esito[1]})")
        intercalati.sort()
        if invertite:
            note.append(
                f"{len(invertite)} file portano nelle intestazioni il marcatore dell'altra "
                f"lettura ({'; '.join(sorted(invertite)[:3])}"
                f"{'...' if len(invertite) > 3 else ''}): i due file di quelle corse sono "
                "scambiati, e validate si fermera' (E-S0-26) finche' non se ne scambiano i nomi.")
    else:
        intercalati = sorted(
            f"{n} ({i.statistiche.coppie_nello_stesso_file})" for n, i in ispezioni.items()
            if i.statistiche.valido and i.statistiche.coppie_nello_stesso_file)
    if intercalati:
        raise Rifiuto(_letture_appaiate(
            f"{len(intercalati)} file non contengono le sole prime letture delle coppie: "
            + "; ".join(intercalati[:3]) + ("..." if len(intercalati) > 3 else "")
        ) if not reverse_di else (
            f"{len(intercalati)} file contengono le due letture di ogni coppia nello stesso "
            "file: " + "; ".join(intercalati[:3]) + ("..." if len(intercalati) > 3 else "")
            + ". In modalita' appaiata ogni file riporta una sola lettura di ogni coppia: "
            "separa le due letture in due file per campione"))
    guasti = sorted(n for n, i in ispezioni.items()
                    if not i.statistiche.valido and not i.statistiche.vuoto)
    vuoti = sorted(n for n, i in ispezioni.items() if i.statistiche.vuoto)
    if guasti:
        note.append(
            f"{len(guasti)} file non si leggono come FASTQ ({', '.join(guasti[:3])}"
            f"{'...' if len(guasti) > 3 else ''}): validate li segnalera' (E-S0-13).")
    if vuoti:
        note.append(
            f"{len(vuoti)} file non contengono letture ({', '.join(vuoti[:3])}"
            f"{'...' if len(vuoti) > 3 else ''}): ammesso solo per i controlli negativi.")
    quante = "tutte le letture" if args.head_reads is None else f"le prime {args.head_reads} letture"

    metti("io.fastq_dir", INDICATO, str(cartella), "--reads")
    if reverse_di:
        metti("io.layout", RILEVATO, defaults.LAYOUT_APPAIATO,
              f"ognuno dei {2 * len(file)} file della cartella ha il suo compagno: "
              f"{len(file)} coppie di file che si distinguono per {marcatori[0]} e {marcatori[1]}",
              "per usare le sole letture forward: \"single-end\", togliendo il gruppo pair e "
              "lasciando nella cartella (o in io.fastq_glob) i soli file forward")
        metti("pair.file_markers", RILEVATO, list(marcatori),
              "i due marcatori con cui i nomi dei file di ogni coppia si distinguono, prima "
              "quello del file forward",
              "abbinamento, per esempio: " + "; ".join(
                  f"{n} e {reverse_di[n].name}" for n in nomi[:2]))
    else:
        metti("io.layout", RILEVATO, defaults.LAYOUT_SINGOLO,
              "nessun file della cartella ha un compagno che se ne distingua per un "
              "marcatore di coppia (R1 e R2, 1 e 2): un file per campione")
    if len(per_estensione) == 1:
        (estensione,) = per_estensione
        metti("io.fastq_glob", RILEVATO, f"*{estensione}",
              f"i {len(file)} file di letture della cartella hanno tutti questa estensione")
    else:
        metti("io.fastq_glob", SCELTA, SENZA_VALORE,
              "la cartella contiene file di letture con estensioni diverse ("
              + ", ".join(f"{len(l)} {e}" for e, l in per_estensione.items())
              + "): indica il modello che seleziona quelli da usare")

    # ---- Tabelle ----------------------------------------------------------
    assay, studio = _tabelle(args)
    try:
        colonne_assay, righe_assay = r.leggi_tabella(assay)
        colonne_studio, righe_studio = r.leggi_tabella(studio) if studio else ([], None)
    except ValueError as guasto:
        raise Rifiuto(str(guasto)) from guasto
    metti("io.assay_table", INDICATO, str(assay), "--metadata")
    metti("io.study_table", INDICATO, str(studio) if studio else None,
          "--study-table, o la tabella s_*.txt della cartella ISA-Tab" if studio else
          "nessuna tabella di studio: tipo e variabili si leggono dalla tabella dei campioni")
    if args.out_root is not None:
        metti("io.out_root", INDICATO, str(args.out_root), "--out-root")

    # ---- Chiave del campione ---------------------------------------------
    candidate = r.chiave_del_campione(nomi, righe_assay, colonne_assay)
    # Prima le colonne con il nome di un file per cella, poi quelle che ne
    # elencano piu' d'uno, poi le altre.
    di_file = [c for c in candidate if c.di_file] or [c for c in candidate if c.elenca_file]
    migliori = di_file or candidate
    chiave: r.Chiave | None = migliori[0] if len(migliori) == 1 else None

    # ---- Righe e file: che cosa la tabella dichiara, che cosa le letture portano
    # Dove la tabella dichiara per riga l'amplicone o il primer, lo si confronta
    # con il primer in testa alle letture del file collegato, su tutte le righe:
    # una riga di un altro amplicone collegata a un file di questa analisi va
    # vista prima di proporre quali righe tenere.
    tutte_le_righe = righe_assay
    colonna_amplicone = (
        ampliconi.COLONNA_AMPLICONE if colonne_assay.count(ampliconi.COLONNA_AMPLICONE) == 1
        else None)
    colonna_primer = (
        ampliconi.COLONNA_PRIMER if colonne_assay.count(ampliconi.COLONNA_PRIMER) == 1 else None)
    presenti = {p.name: p for p in (*file, *reverse_di.values())}
    # La colonna che collega le righe ai file: quella della chiave, o la prima
    # i cui valori nominano file della cartella.
    colonna_dei_file = chiave.colonna if chiave is not None else next(
        (c for c in colonne_assay
         if any(ampliconi.file_della_cella(x.get(c, ""), presenti)[0] for x in righe_assay)),
        None)
    nomi_provvisori = r.colonna_del_campione(
        righe_assay, colonne_assay, righe_studio, colonne_studio)
    colonna_dei_nomi = nomi_provvisori[0][0] if len(nomi_provvisori) == 1 else colonne_assay[0]
    incoerenze: list[ampliconi.Incoerenza] = []
    senza_file: list[tuple[str, str, ampliconi.Dichiarazione]] = []
    if colonna_dei_file is not None and (colonna_amplicone or colonna_primer):
        della_chiave = chiave

        def per_chiave(valore: str) -> Path | None:
            if della_chiave is None:
                return None
            nome_del_file = della_chiave.file_per_chiave.get(
                _chiave_della_riga(della_chiave, {della_chiave.colonna: valore}) or "")
            return presenti.get(nome_del_file or "")

        collegate, senza_file = ampliconi.collega(
            righe_assay, colonna_dei_nomi, colonna_dei_file, presenti,
            colonna_amplicone, colonna_primer, per_chiave)
        avanzamento("confronto fra cio' che i metadati dichiarano e il primer in testa alle "
                    "letture del file collegato...")
        incoerenze, giudicate = ampliconi.confronta(collegate)
        dichiarate = " e ".join(
            repr(c) for c in (colonna_amplicone, colonna_primer) if c is not None)
        if incoerenze:
            elenco = "; ".join(str(i) for i in incoerenze[:8]) + (
                f"; e altre {len(incoerenze) - 8}" if len(incoerenze) > 8 else "")
            note.append(
                f"{len(incoerenze)} righe della tabella sono collegate a un file le cui "
                f"letture portano il primer di un altro amplicone: {elenco}. La chiave del "
                "campione non e' proposta: correggi una copia della tabella (mai "
                "l'originale), collegando ogni campione al suo file, e rilancia init sulla "
                "copia. Anche le righe non elencate possono essere collegate al file di un "
                "altro campione dello stesso amplicone, che dal contenuto non si riconosce: "
                "l'abbinamento va verificato per intero. Con questa tabella validate si "
                "fermerebbe (E-S0-30).")
            chiave, migliori = None, []
        elif not giudicate:
            note.append(
                f"I metadati dichiarano per riga l'amplicone o il primer ({dichiarate}), ma il "
                "confronto con le letture non si e' potuto fare: "
                + ("nessuna riga e' collegata a un file della cartella" if not collegate else
                   "nessuna dichiarazione si interpreta senza ambiguita'"
                   if not any(c.dichiarazione.interpretabile for c in collegate) else
                   "le letture non portano in testa un primer riconosciuto (di solito perche' "
                   "e' gia' stato tolto)")
                + ". Che ogni riga sia collegata al suo file resta da verificare.")
        elif giudicate < len(collegate):
            # Una riga non giudicata puo' essere collegata al file di un altro
            # amplicone senza che lo si veda: va detto anche quando le altre
            # righe risultano coerenti.
            note.append(
                f"{len(collegate) - giudicate} righe su {len(collegate)} collegate a un file "
                "non si sono potute confrontare con le letture (la riga non dichiara "
                "amplicone o primer in modo interpretabile, o il file non porta in testa un "
                "primer riconosciuto): che siano collegate al loro file resta da verificare.")
        if senza_file:
            elenco = "; ".join(
                f"{nome} ({nome_file}"
                + (f", {d.amplicone}" if d.amplicone is not None else "") + ")"
                for nome, nome_file, d in senza_file[:8]) + (
                f"; e altre {len(senza_file) - 8}" if len(senza_file) > 8 else "")
            note.append(
                f"{len(senza_file)} righe della tabella nominano un file che non e' nella "
                f"cartella delle letture: {elenco}.")

    # ---- Righe di questa analisi ----------------------------------------------
    # Se una colonna dichiara l'amplicone e la tabella ne descrive piu' d'uno,
    # si propone di tenere le sole righe 16S; la chiave e' gia' stata cercata fra
    # le sole righe che nominano un file.
    tenute: list[int] | None = None
    if colonna_amplicone is not None:
        valori_amplicone = [x.get(colonna_amplicone, "") for x in righe_assay]
        del_16s = sorted({x for x in valori_amplicone
                          if ampliconi.amplicone_dichiarato(x) == "16S"})
        altri = sorted({x for x in valori_amplicone if x not in del_16s})
        if del_16s and altri:
            tenute = [i for i, x in enumerate(valori_amplicone) if x in del_16s]
            metti("meta.assay_keep_column", PROPOSTO, colonna_amplicone,
                  f"la colonna dichiara l'amplicone di ogni riga, e la tabella ne descrive "
                  f"piu' d'uno: {len(tenute)} righe {', '.join(map(repr, del_16s))} e "
                  f"{len(righe_assay) - len(tenute)} righe "
                  + ", ".join(repr(x or '<vuoto>') for x in altri),
                  "le righe che non portano uno dei valori di meta.assay_keep_values sono "
                  "escluse dalla validazione, con un avviso che le elenca (E-S0-31)")
            metti("meta.assay_keep_values", PROPOSTO, del_16s,
                  "i valori che dichiarano il 16S: le sole righe che questa analisi tratta")
            degli_esclusi = sorted({
                nome_file for i, riga in enumerate(righe_assay) if i not in tenute
                for nome_file in [ampliconi.file_della_cella(
                    riga.get(colonna_dei_file or "", ""), presenti)[0]] if nome_file})
            if degli_esclusi and not incoerenze:
                note.append(
                    f"{len(degli_esclusi)} file della cartella appartengono alle righe che la "
                    f"regola esclude ({', '.join(degli_esclusi[:5])}"
                    f"{'...' if len(degli_esclusi) > 5 else ''}): vanno tolti dalla cartella "
                    "delle letture, o lasciati fuori da io.fastq_glob, altrimenti validate si "
                    "ferma (E-S0-06: un file senza riga).")
            note.append(
                f"La tabella descrive anche campioni di un altro amplicone "
                f"({', '.join(repr(x or '<vuoto>') for x in altri)}: "
                f"{len(righe_assay) - len(tenute)} righe): proposto di tenere le sole "
                f"{len(tenute)} righe 16S (meta.assay_keep_column, meta.assay_keep_values).")
    if chiave is not None and chiave.righe is not None:
        con_file = set(chiave.righe)
        fuori = [i for i in range(len(righe_assay))
                 if i not in con_file and (tenute is None or i in tenute)]
        if fuori:
            elenco = ", ".join(
                repr(righe_assay[i].get(colonna_dei_nomi, "")) for i in fuori[:8]) + (
                f" e altre {len(fuori) - 8}" if len(fuori) > 8 else "")
            note.append(
                f"{len(fuori)} righe della tabella non nominano alcun file della cartella "
                f"({elenco}): la chiave e' cercata fra le altre {len(con_file)}. validate si "
                "ferma su quelle righe (E-S0-06): toglile da una copia della tabella, oppure, "
                "se una colonna le distingue, dichiara meta.assay_keep_column e "
                "meta.assay_keep_values.")
        righe_assay = [righe_assay[i] for i in sorted(con_file)]
    elif tenute is not None:
        righe_assay = [righe_assay[i] for i in tenute]

    # ---- File del lotto ----------------------------------------------------
    lotto: Path | None = args.batch_table
    colonne_lotto: list[str] = []
    righe_lotto: list[dict[str, str]] = []
    if lotto is not None:
        if not lotto.is_file():
            raise Rifiuto(f"--batch-table: {lotto} non esiste")
        try:
            colonne_lotto, righe_lotto = r.leggi_tabella(lotto)
        except ValueError as guasto:
            raise Rifiuto(str(guasto)) from guasto
    metti("io.batch_table", INDICATO, str(lotto) if lotto else None,
          "--batch-table" if lotto else "nessun file del lotto indicato: corsa singola, nessun lotto")

    colonne_chiave_lotto: list[str] = []
    if chiave is not None and lotto is not None:
        # Fra le espressioni verificate si preferisce quella con cui anche il
        # file del lotto si aggancia: la chiave deve essere la stessa ovunque.
        tutte = (chiave.espressione, *chiave.alternative)
        for espressione in tutte:
            chiavi = r.chiavi_dei_file(espressione, nomi)
            trovate = r.chiave_del_lotto(chiavi.values(), espressione, righe_lotto, colonne_lotto)
            if trovate:
                chiave = r.Chiave(
                    chiave.colonna, espressione, chiave.di_file,
                    {c: n for n, c in chiavi.items()},
                    tuple(e for e in tutte if e != espressione))
                colonne_chiave_lotto = trovate
                break

    if chiave is not None:
        # Una chiave di sole cifre o cortissima puo' corrispondere ai file per
        # coincidenza: l'abbinamento e' allora una proposta, da guardare.
        specie_chiave = SCELTA if chiave.debole else RILEVATO
        esempi = "abbinamento, per esempio: " + "; ".join(
            f"{nome} -> {c}" for c, nome in list(chiave.file_per_chiave.items())[:3])
        avvertenza = (
            "non proposta: la chiave e' di sole cifre o molto corta, e la corrispondenza "
            "con i file puo' essere una coincidenza. Controlla l'abbinamento: se e' giusto, "
            f"scrivi meta.accession_column: {valore_yaml(chiave.colonna)} e io.accession_regex: "
            f"{valore_yaml(chiave.espressione)}"
        ) if chiave.debole else ""
        metti("meta.accession_column", specie_chiave,
              SENZA_VALORE if chiave.debole else chiave.colonna,
              "la sola colonna "
              + ("con nomi di file di letture " if chiave.di_file or chiave.elenca_file else "")
              + f"i cui valori corrispondono uno a uno ai file ({len(nomi)})"
              + (f", fra le {len(chiave.righe)} righe su {len(tutte_le_righe)} che nominano "
                 "un file della cartella" if chiave.righe is not None else ""),
              esempi, avvertenza)
        metti("io.accession_regex", specie_chiave,
              SENZA_VALORE if chiave.debole else chiave.espressione,
              f"estrae da ogni file una chiave unica, la stessa che si ricava dalla "
              f"colonna {chiave.colonna!r}: verificato su tutti i file e le righe ({len(nomi)})",
              ("altre espressioni verificate: " + ", ".join(chiave.alternative[:4]))
              if chiave.alternative else "",
              # La chiave e' il nome del campione in ogni uscita: se e' il nome
              # intero del file, lo si dice, con il modo di accorciarla.
              ("la chiave e' il nome intero del file, e diventa l'identificativo di ogni "
               "campione nelle uscite: per una chiave piu' corta si puo' scrivere una delle "
               "altre espressioni verificate, se da' ancora un valore diverso per campione")
              if chiave.alternative and all(
                  c == r._senza_estensione(n) for c, n in chiave.file_per_chiave.items())
              else "", avvertenza)
        if chiave.debole:
            note.append(
                "Chiave del campione non proposta: la sola colonna che corrisponde ai file "
                f"({chiave.colonna!r}) da' una chiave di sole cifre o molto corta, e la "
                "corrispondenza puo' essere una coincidenza. meta.accession_column e "
                "io.accession_regex restano da scegliere: il commento riporta l'abbinamento "
                "trovato e i valori da scrivere se e' giusto.")
            non_indovinati.append(
                "meta.accession_column e io.accession_regex: la chiave trovata e' di sole "
                "cifre o molto corta.")
    elif migliori and not incoerenze:
        # Piu' colonne corrispondono ai file: init non sceglie, ma scrive per
        # ciascuna la coppia di righe pronta (colonna ed espressione), cosi'
        # chi sceglie toglie un commento e non scrive un'espressione regolare.
        # Si consiglia la colonna che riporta per intero il nome dei file
        # della cartella, se e' una sola: e' il collegamento dichiarato fra
        # riga e file, non una parte del nome in comune.
        dei_file = set(nomi) | {p.name for p in reverse_di.values()}
        con_il_nome = [
            c for c in migliori
            if all(any(e in dei_file for e in elementi_di_una_cella(x.get(c.colonna, "")))
                   for i, x in enumerate(tutte_le_righe) if c.righe is None or i in c.righe)]
        consigliata = con_il_nome[0] if len(con_il_nome) == 1 else None

        def etichetta(numero: int, candidata: r.Chiave) -> str:
            return f"coppia {numero}" + (
                " (consigliata: la colonna riporta il nome dei file)"
                if candidata is consigliata else "") + (
                " (chiave di sole cifre o molto corta: la corrispondenza puo' essere una "
                "coincidenza)" if candidata.debole else "")

        come_scegliere = (
            "togli il commento a una sola delle righe qui sotto, e alla riga con lo stesso "
            "numero di coppia sotto {altra}: le due righe di una coppia vanno insieme")
        metti("meta.accession_column", SCELTA, SENZA_VALORE,
              f"{len(migliori)} colonne corrispondono uno a uno ai file, e init non sceglie "
              "fra loro: " + come_scegliere.format(altra="io.accession_regex (gruppo io)"))
        metti("io.accession_regex", SCELTA, SENZA_VALORE,
              "si sceglie insieme a meta.accession_column (gruppo meta): "
              + come_scegliere.format(altra="meta.accession_column"))
        voci["meta.accession_column"] = replace(voci["meta.accession_column"], pronte=tuple(
            (etichetta(n, c) + ": i valori sono come "
             + repr(tutte_le_righe[c.righe[0] if c.righe else 0].get(c.colonna, "")), c.colonna)
            for n, c in enumerate(migliori, 1)))
        voci["io.accession_regex"] = replace(voci["io.accession_regex"], pronte=tuple(
            (etichetta(n, c) + f", con la colonna {c.colonna!r}: abbinamento, per esempio: "
             + "; ".join(f"{nome} -> {k}" for k, nome in list(c.file_per_chiave.items())[:2]),
             c.espressione)
            for n, c in enumerate(migliori, 1)))
        note.append(
            f"Chiave del campione da scegliere: {len(migliori)} colonne corrispondono uno a "
            "uno ai file ("
            + "; ".join(f"coppia {n}: {c.colonna!r}" for n, c in enumerate(migliori, 1))
            + "). Nel file le coppie sono righe pronte sotto meta.accession_column e "
            "io.accession_regex: si toglie il commento alle due righe della coppia scelta"
            + (f"; consigliata la coppia {migliori.index(consigliata) + 1}, la colonna che "
               "riporta il nome dei file" if consigliata is not None else
               "; nessuna e' consigliata: guarda l'abbinamento scritto accanto a ciascuna")
            + ".")
    else:
        if incoerenze:
            origine = (
                f"non proposta: {len(incoerenze)} righe della tabella sono collegate a un "
                "file le cui letture portano il primer di un altro amplicone (elencate nel "
                "riepilogo): si corregge una copia della tabella, poi si sceglie la chiave")
        else:
            origine = (
                "nessuna colonna della tabella corrisponde uno a uno ai file: i valori "
                f"devono contenere una parte del nome del file unica per campione ("
                f"{len(nomi)} file, {len(tutte_le_righe)} righe)")
        metti("meta.accession_column", SCELTA, SENZA_VALORE, origine)
        metti("io.accession_regex", SCELTA, SENZA_VALORE,
              "si sceglie insieme a meta.accession_column: la parte del nome del file "
              "che e' unica per campione")
        if not incoerenze:
            note.append("Chiave del campione non rilevata: " + origine + ".")
        else:
            non_indovinati.append(
                "meta.accession_column e io.accession_regex: la tabella collega alcune righe "
                "al file di un altro amplicone.")

    # ---- Nome del campione --------------------------------------------------
    agganci = r.colonna_del_campione(righe_assay, colonne_assay, righe_studio, colonne_studio)
    colonna_campione: str | None = None
    colonna_studio: str | None = None
    if len(agganci) == 1:
        colonna_campione, colonna_studio = agganci[0]
        metti("meta.sample_id_column", RILEVATO, colonna_campione,
              "la sola colonna che aggancia ogni riga della tabella di assay a una e una "
              "sola riga della tabella di studio")
        if colonna_studio != colonna_campione:
            metti("meta.study_sample_id_column", RILEVATO, colonna_studio,
                  "la colonna della tabella di studio che corrisponde a meta.sample_id_column")
    elif righe_studio is not None:
        metti("meta.sample_id_column", SCELTA, SENZA_VALORE,
              (f"{len(agganci)} coppie di colonne agganciano le due tabelle: "
               + "; ".join(f"{a!r} con {s!r}" for a, s in agganci[:6])) if agganci else
              "nessuna colonna aggancia ogni riga della tabella di assay a una sola riga "
              "della tabella di studio")
    else:
        uniche = [c for c in colonne_assay
                  if c != (chiave.colonna if chiave else None)
                  and len({x.get(c, "") for x in righe_assay} - {""}) == len(righe_assay)]
        if uniche:
            colonna_campione = uniche[0]
            metti("meta.sample_id_column", PROPOSTO, colonna_campione,
                  "la prima colonna con un valore diverso per ogni riga",
                  ("altre colonne con valori tutti distinti: " + ", ".join(map(repr, uniche[1:5])))
                  if uniche[1:] else "")
        else:
            metti("meta.sample_id_column", SCELTA, SENZA_VALORE,
                  "nessuna colonna, oltre a quella della chiave, ha un valore diverso per "
                  "ogni riga")

    # ---- Tipo di campione ----------------------------------------------------
    righe_classe = righe_assay
    colonne_classe = colonne_assay
    if righe_studio is not None:
        colonne_classe = colonne_studio
        if colonna_campione is not None:
            nomi_assay = {x.get(colonna_campione, "") for x in righe_assay}
            righe_classe = [x for x in righe_studio if x.get(colonna_studio or "", "") in nomi_assay]
        else:
            righe_classe = righe_studio
    escluse = {chiave.colonna if chiave else None, colonna_campione, colonna_studio}
    tipi = r.colonne_del_tipo(righe_classe, [c for c in colonne_classe if c not in escluse])
    negativi_proposti: set[str] = set()
    con_positivi = False
    if tipi and tipi[0].da == "parole" and not tipi[0].biologici:
        # Ogni valore della colonna porta una parola dei controlli: accettarla
        # classificherebbe come controllo ogni campione, e non ne resterebbe
        # nessuno da analizzare. Succede quando i campioni sono essi stessi lo
        # standard (una comunita' a composizione nota descritta come «controllo
        # positivo»): se lo siano lo sa solo chi conosce lo studio.
        tipo = tipi[0]
        metti("ctrl.column", SCELTA, SENZA_VALORE,
              f"non proposta: con la colonna {tipo.colonna!r} tutti i campioni "
              "risulterebbero controlli (ogni suo valore ha una parola dei controlli: "
              + ", ".join(map(repr, (tipo.negativi + tipo.positivi + tipo.non_attribuiti)[:5]))
              + ") e non ne resterebbe nessuno da analizzare",
              "se i campioni sono lo standard da analizzare, e non i controlli di altri "
              "campioni, scrivi ctrl.column: null (nessun controllo, ogni campione e' "
              "biologico) e lascia in commento le tre etichette",
              "se invece una colonna distingue i campioni dai controlli, scrivi qui il suo "
              "nome e nelle tre etichette i suoi valori")
        for nome in ("blank_values", "positive_values", "biological_values"):
            metti(f"ctrl.{nome}", SCELTA, SENZA_VALORE, "si sceglie insieme a ctrl.column")
        non_indovinati.append(
            "ctrl.column: la colonna trovata darebbe tutti i campioni per controlli.")
    elif tipi:
        tipo = tipi[0]
        altre = ("altre colonne candidate: " + ", ".join(repr(t.colonna) for t in tipi[1:5])
                 ) if tipi[1:] else ""
        come = ("ha valori con le parole dei controlli" if tipo.da == "parole" else
                "il nome dice il tipo del campione; nessun valore ha le parole dei controlli")
        colonna_dei_nomi_di_classe = colonna_studio or colonna_campione
        con_positivi = bool(tipo.positivi)
        nome_del_negativo: dict[str, str] = {}
        ricchi: list[str] = []
        mediana_biologici = 0.0
        letture_di: dict[str | None, int] = {}
        if chiave is not None and colonna_campione is not None:
            tipo_per_nome = {x.get(colonna_studio or colonna_campione, ""): x.get(tipo.colonna, "")
                             for x in righe_classe}
            for riga in righe_assay:
                if tipo_per_nome.get(riga.get(colonna_campione, "")) in tipo.negativi:
                    della_riga = _chiave_della_riga(chiave, riga)
                    if della_riga is not None:
                        negativi_proposti.add(della_riga)
                        nome_del_negativo[della_riga] = riga.get(colonna_campione, "")
            # Un controllo negativo con piu' letture della mediana degli altri
            # campioni si nomina, senza togliere la proposta: con campioni a
            # bassa biomassa i bianchi hanno spesso tante letture quante i
            # campioni, e la colonna e' comunque quella giusta.
            if negativi_proposti and args.head_reads is None:
                letture_di = {
                    _chiave_del_file(chiave, n): i.statistiche.letture_esaminate
                    for n, i in ispezioni.items() if i.statistiche.valido}
                dei_biologici = [
                    n for c, n in letture_di.items() if c not in negativi_proposti]
                if dei_biologici:
                    mediana_biologici = statistics.median(dei_biologici)
                    ricchi = sorted(
                        c for c in negativi_proposti
                        if c is not None and letture_di.get(c, 0) > mediana_biologici)
        # La proposta e' sicura solo se nulla la contraddice: le altre candidate
        # danno ai campioni le stesse classi e ogni valore e' attribuito.
        # Altrimenti la colonna resta da scegliere, con cio' che serve a
        # scegliere.
        # Una colonna con un valore diverso per ogni riga e' un nome, non una
        # classificazione: una parola dei controlli in un nome non contraddice
        # la colonna del tipo.
        discordi = [
            t.colonna for t in tipi[1:]
            if t.biologici
            and len({x.get(t.colonna, "") for x in righe_classe}) < len(righe_classe)
            and _classi_dei_campioni(righe_classe, t) != _classi_dei_campioni(
                righe_classe, tipo)]
        dubbi = [
            *([f"{len(discordi)} altre colonne candidate classificano i campioni in modo "
               f"diverso ({', '.join(map(repr, discordi[:4]))})"] if discordi else []),
            *(["alcuni valori hanno la sola parola «control» e non si attribuiscono a una "
               "classe (" + ", ".join(map(repr, tipo.non_attribuiti)) + ")"]
              if tipo.non_attribuiti else []),
        ]
        dei_ricchi = (
            ("1 campione che ne risulta controllo negativo ha " if len(ricchi) == 1 else
             f"{len(ricchi)} campioni che ne risultano controlli negativi hanno ")
            + f"piu' letture della mediana degli altri campioni ({mediana_biologici:.0f}): "
            + ", ".join(f"{nome_del_negativo.get(c) or c} ({letture_di[c]})" for c in ricchi[:8])
            + (f" e altri {len(ricchi) - 8}" if len(ricchi) > 8 else "")
            + ". Con campioni a bassa biomassa e' normale; altrimenti controlla che siano "
            "davvero bianchi") if ricchi else ""
        if ricchi:
            note.append(f"Controlli negativi con molte letture (colonna {tipo.colonna!r}): "
                        + dei_ricchi + ".")
        per_valore = (
            f"valori di {tipo.colonna!r} e campioni che li portano: "
            + _campioni_per_valore(righe_classe, tipo.colonna, colonna_dei_nomi_di_classe),
            altre,
            *(f"candidata {t.colonna!r}, valori e campioni: "
              + _campioni_per_valore(righe_classe, t.colonna, colonna_dei_nomi_di_classe)
              for t in tipi[1:5]))
        etichette = (
            ("blank_values", list(tipo.negativi),
             "i valori con una parola dei controlli negativi ("
             + ", ".join(v.PAROLE_NEGATIVO) + ")" if tipo.negativi else
             "nessun valore ha una parola dei controlli negativi", ""),
            ("positive_values", list(tipo.positivi),
             "i valori con una parola dei controlli positivi ("
             + ", ".join(v.PAROLE_POSITIVO) + ")" if tipo.positivi else
             "nessun valore ha una parola dei controlli positivi", ""),
            ("biological_values", list(tipo.biologici), "i valori senza parole dei controlli",
             ("valori con la sola parola «control», non attribuiti: da assegnare a mano a "
              "una delle tre etichette: " + ", ".join(map(repr, tipo.non_attribuiti)))
             if tipo.non_attribuiti else ""),
        )
        if not dubbi:
            metti("ctrl.column", PROPOSTO, tipo.colonna, f"colonna proposta perche' {come}",
                  "la proposta guarda le parole dei valori e puo' sbagliare: controlla quali "
                  "campioni porta ogni valore, qui e nelle altre candidate", dei_ricchi,
                  *per_valore)
            for nome, valori, origine, altro in etichette:
                metti(f"ctrl.{nome}", PROPOSTO, valori, origine, altro)
        else:
            perche = "; ".join(dubbi)
            metti("ctrl.column", SCELTA, SENZA_VALORE,
                  f"non proposta: la colonna piu' plausibile e' {tipo.colonna!r} ({come}), ma "
                  f"{perche}",
                  f"se e' la colonna giusta scrivi {tipo.colonna!r} qui e, nelle tre "
                  "etichette, "
                  + "; ".join(f"ctrl.{nome}: {valore_yaml(valori)}" for nome, valori, *_ in etichette)
                  + "; altrimenti scegli fra le candidate guardando quali campioni porta "
                  "ogni valore",
                  "motivo conservato, troncamento e note di questo file sono calcolati "
                  f"escludendo i campioni che {tipo.colonna!r} darebbe per controlli negativi: "
                  "con un'altra colonna conviene rilanciare init dopo aver corretto la tabella, "
                  "o verificare quelle righe a mano", dei_ricchi, *per_valore)
            for nome, valori, origine, altro in etichette:
                metti(f"ctrl.{nome}", SCELTA, SENZA_VALORE,
                      f"si sceglie insieme a ctrl.column; con {tipo.colonna!r} sarebbe "
                      f"{valore_yaml(valori)}: {origine}", altro)
            note.append(
                f"Colonna dei controlli non proposta: la piu' plausibile e' {tipo.colonna!r}, ma "
                f"{perche}. ctrl.column e le tre etichette restano da scegliere: il commento "
                "elenca i campioni di ogni valore, per questa colonna e per le altre candidate.")
            non_indovinati.append(
                "ctrl.column e le tre etichette: la colonna trovata non e' una proposta sicura.")
    else:
        metti("ctrl.column", PROPOSTO, None,
              "nessuna colonna ha valori con le parole dei controlli ne' un nome che dica "
              "il tipo del campione: nullo vuol dire nessun controllo, ogni campione e' "
              "biologico")
        for nome in ("blank_values", "positive_values", "biological_values"):
            metti(f"ctrl.{nome}", PROPOSTO, [], "senza la colonna del tipo le etichette sono vuote")
    metti("ctrl.blank_override_column", SCELTA, None,
          "partenza: nessuna riclassificazione. Campioni dichiarati di un tipo ma da "
          "trattare come controlli negativi (per esempio materiale mai usato) si "
          "riconoscono solo conoscendo lo studio: init non li puo' rilevare")
    metti("ctrl.blank_override_values", SCELTA, [], "partenza: nessuna riclassificazione")

    # ---- Modulo, colonne da portare nell'oggetto ---------------------------------
    metti("meta.derive_module", SCELTA, False,
          "partenza: nessun raggruppamento di posizioni da derivare")
    metti("meta.module_regex", SCELTA, None,
          "partenza: nessuna; serve solo se le posizioni hanno un prefisso che dice il gruppo")
    metti("meta.module_column", SCELTA, None,
          "partenza: nessuna colonna della posizione di prelievo",
          "colonne della tabella: " + ", ".join(map(repr, colonne_classe[:20])))
    metti("meta.non_surface_positions", SCELTA, [],
          "partenza: nessuna; ha senso solo con meta.module_column")
    # Nelle tabelle ISA-Tab le variabili dello studio sono dichiarate come
    # tali: le colonne «Factor Value[...]». Portarle nell'oggetto finale e'
    # cio' che permette alle analisi a valle di confrontare i gruppi; una
    # colonna gia' usata per il tipo del campione non si ripete.
    gia_usate = {voci["ctrl.column"].valore} if "ctrl.column" in voci else set()
    fattori = [c for c in colonne_classe
               if c.startswith(PREFISSO_DEI_FATTORI) and colonne_classe.count(c) == 1
               and c not in gia_usate]
    if fattori:
        metti("out.study_columns", PROPOSTO, fattori,
              f"le colonne «{PREFISSO_DEI_FATTORI}...]» della tabella: in ISA-Tab sono le "
              "variabili dello studio, quelle che le analisi a valle confrontano",
              "altre colonne disponibili: " + ", ".join(
                  repr(c) for c in colonne_classe[:20] if c not in fattori))
    else:
        metti("out.study_columns", SCELTA, [],
              "partenza: nessuna colonna oltre a quelle dell'inventario; elenca qui le "
              "variabili dello studio che le analisi a valle useranno",
              "colonne disponibili: " + ", ".join(map(repr, colonne_classe[:20])))
        note.append(
            "out.study_columns e' vuoto: nessuna variabile dello studio entra nell'oggetto "
            "finale, e le analisi ecologiche non avranno gruppi da confrontare. Se la tabella "
            "ha colonne con le variabili dello studio, scrivile in quella riga (il commento "
            "elenca le colonne disponibili).")

    # ---- File del lotto: colonne ---------------------------------------------------
    if lotto is None:
        for nome in ("meta.batch_key_column", "meta.batch_module_column", "err.batch_column",
                     "decontam.batch_column"):
            metti(nome, RILEVATO, None,
                  "nessun file del lotto indicato: e' il solo valore ammesso (E-G15-11)")
        metti("out.batch_columns", RILEVATO, [],
              "nessun file del lotto indicato: e' il solo valore ammesso (E-G15-11)")
        metti("katharoseq.cell_count_column", SCELTA, None,
              "partenza: nessuna; senza file del lotto la colonna delle cellule dei "
              "controlli positivi andrebbe cercata nella tabella dei campioni")
    else:
        if len(colonne_chiave_lotto) == 1:
            metti("meta.batch_key_column", RILEVATO, colonne_chiave_lotto[0],
                  "la sola colonna del file del lotto che da' per ogni campione la sua "
                  "chiave, in una e una sola riga")
        else:
            metti("meta.batch_key_column", SCELTA, SENZA_VALORE,
                  (f"{len(colonne_chiave_lotto)} colonne del file del lotto danno la chiave "
                   "di ogni campione: " + ", ".join(map(repr, colonne_chiave_lotto)))
                  if colonne_chiave_lotto else
                  "nessuna colonna del file del lotto da' per ogni campione la chiave "
                  "estratta dai file" + ("" if chiave else " (la chiave non e' stata rilevata)"))
        del_lotto = r.colonne_del_lotto(
            colonne_lotto, righe_lotto,
            colonne_chiave_lotto[0] if len(colonne_chiave_lotto) == 1 else None)
        for nome, tipo_colonna, che_cosa in (
            ("decontam.batch_column", "piastra", "la piastra di estrazione"),
            ("err.batch_column", "corsa", "la corsa di sequenziamento"),
            ("katharoseq.cell_count_column", "cellule", "le cellule dei controlli positivi"),
        ):
            trovate = del_lotto[tipo_colonna]
            # Piu' colonne candidate non sono un dubbio se raggruppano le righe
            # allo stesso modo (un numero e una sigla della stessa piastra).
            se_equivalenti = len({_raggruppamento(righe_lotto, c) for c in trovate}) == 1
            if trovate and se_equivalenti:
                metti(nome, PROPOSTO, trovate[0],
                      f"il nome della colonna dice {che_cosa}",
                      ("altre colonne candidate, che raggruppano i campioni allo stesso modo: "
                       + ", ".join(map(repr, trovate[1:]))) if trovate[1:] else "")
            elif trovate:
                metti(nome, SCELTA, SENZA_VALORE,
                      f"non proposta: {len(trovate)} colonne del file del lotto "
                      f"hanno un nome che dice {che_cosa}, raggruppano i campioni in modo "
                      "diverso e init non sceglie fra loro: "
                      + ", ".join(map(repr, trovate)),
                      f"raccomandata la prima, {trovate[0]!r}, se e' quella che descrive "
                      f"{che_cosa} di ogni campione: guarda i valori delle candidate; nullo "
                      "se nessuna lo e'")
                note.append(
                    f"{nome} non proposto: {len(trovate)} colonne del file del lotto dicono "
                    f"{che_cosa} e raggruppano i campioni in modo diverso "
                    f"({', '.join(map(repr, trovate))}).")
            else:
                metti(nome, SCELTA, None,
                      f"partenza: nessuna; nessuna colonna del file del lotto ha un nome "
                      f"che dica {che_cosa}",
                      "colonne del file: " + ", ".join(map(repr, colonne_lotto)))
        metti("meta.batch_module_column", SCELTA, None,
              "partenza: nessuna colonna del modulo nel file del lotto")
        metti("out.batch_columns", SCELTA, [],
              "partenza: nessuna colonna del file del lotto oltre a piastra e corsa",
              "colonne disponibili: " + ", ".join(map(repr, colonne_lotto)))
    con_cellule = con_positivi and voci["katharoseq.cell_count_column"].valore not in (None, SENZA_VALORE)

    # ---- Primer e motivo -----------------------------------------------------------
    osservato, parziali = r.primer_nelle_letture(ispezioni, cercati)
    da_considerare = [
        i for n, i in ispezioni.items()
        if chiave is None or _chiave_del_file(chiave, n) not in negativi_proposti
    ] or list(ispezioni.values())
    esclusi = len(ispezioni) - len(da_considerare)
    senza_negativi = (
        f"esclusi {esclusi} file dei campioni proposti come controlli negativi" if esclusi else "")
    con_etichetta_reverse: list[str] = []
    con_etichetta_forward: list[str] = []
    _, sequenze, _ = r.primer_dai_metadati(
        [righe_assay] + ([righe_studio] if righe_studio is not None else []),
        con_etichetta_reverse, con_etichetta_forward)
    # Le sequenze dei metadati candidate a primer forward: quelle che i metadati
    # dichiarano tali, e quelle che coincidono con un primer forward
    # dell'elenco, il cui verso e' noto anche senza etichetta.
    candidati_forward = list(dict.fromkeys(
        [*con_etichetta_forward,
         *(s for s in sequenze if s not in con_etichetta_reverse
           and r.primer_dell_elenco(s, v.PRIMER_16S) is not None)]))
    primer_noto: str | None = None
    taglio = 0
    if cercato is not None:
        trovato = cercato.forward.primer
        primer_noto = trovato.sequenza
        _voci_del_taglio_cercato(metti, note, cercato, bool(reverse_di))
    elif osservato is not None:
        trovato = osservato.primer
        primer_noto, taglio = trovato.sequenza, len(trovato.sequenza)
        dove = (f"in testa alle letture di ogni file: dal {osservato.minima:.0%} al "
                f"{osservato.massima:.0%} delle letture, su {quante}")
        # E' un fatto solo se vale per ogni file e su tutte le letture: un
        # taglio fisso toglierebbe basi vere ai file che non portano il primer.
        if osservato.in_ogni_file and args.head_reads is None:
            specie_primer, limite = RILEVATO, ""
        elif not osservato.in_ogni_file:
            sotto = osservato.sotto_soglia
            specie_primer = PROPOSTO
            dove = (f"in testa alle letture nella mediana dei file ({osservato.mediana:.0%}), "
                    f"ma non in tutti: dal {osservato.minima:.0%} al {osservato.massima:.0%}")
            limite = (
                f"attenzione: {len(sotto)} file su {len(ispezioni)} non lo portano in testa "
                f"({', '.join(sotto[:5])}{'...' if len(sotto) > 5 else ''}): un taglio fisso "
                "toglierebbe loro basi vere. Vanno capiti prima di accettare")
            note.append(
                f"Il primer {trovato.nome} e' in testa alle letture della maggior parte dei "
                f"file ma non di {len(sotto)}: il primer e' proposto, il taglio resta da scegliere.")
        else:
            specie_primer = PROPOSTO
            limite = ("verificato sulle sole prime letture (--head-reads): senza l'opzione "
                      "init lo verifica su tutte")
            note.append(
                f"Primer {trovato.nome} e taglio sono verificati sulle sole prime "
                f"{args.head_reads} letture di ogni file (--head-reads): restano proposte, e "
                "non fatti. Per verificarli su tutte le letture si rilancia init senza "
                "l'opzione.")
        metti("qc.primer_sequence", specie_primer, trovato.sequenza,
              f"primer {trovato.nome} ({trovato.fonte}), {dove}", limite)
        if osservato.in_ogni_file:
            metti("filter.trimLeft", specie_primer, taglio,
                  f"la lunghezza del primer {trovato.nome}, che le letture portano in testa",
                  limite)
        else:
            # Il primer e' riconosciuto, ma un taglio fisso toglierebbe basi
            # vere ai file che non lo portano: non e' una proposta sicura.
            metti("filter.trimLeft", SCELTA, SENZA_VALORE,
                  f"non proposto: il taglio sarebbe {taglio}, la lunghezza del primer "
                  f"{trovato.nome}, ma non tutti i file lo portano in testa", limite,
                  f"scrivi {taglio} dopo aver capito quei file (e averli tolti o sostituiti, "
                  "se sono di un'altra regione o contengono letture inverse); motivo "
                  "conservato e troncamento di questo file presuppongono quel taglio")
            non_indovinati.append(
                f"filter.trimLeft: il primer {trovato.nome} non e' in testa alle letture di "
                "ogni file.")
    else:
        variabile = "; ".join(f"{p.nome} nel {f:.1%}" for p, f in parziali)
        metti("filter.trimLeft", PROPOSTO, 0,
              f"nessun primer dell'elenco e' in testa alle letture (su {quante}); un primer "
              "fuori elenco va indicato a mano, con la sua lunghezza qui",
              "primer cercati: " + ", ".join(p.nome for p in cercati),
              (f"attenzione: una parte delle letture inizia con un primer ({variabile}): e' "
               "in posizione variabile, e un taglio fisso non lo toglie") if variabile else "")
        if variabile:
            note.append(
                f"Una parte delle letture inizia con un primer dell'elenco ({variabile}), "
                "sotto la soglia per dirlo in testa: e' probabilmente preceduto da basi di "
                "lunghezza variabile e va rimosso per sequenza prima della pipeline.")
        nominati = r.primer_nominati(
            [righe_assay] + ([righe_studio] if righe_studio is not None else []))
        if args.primer is not None:
            primer_noto = args.primer
            metti("qc.primer_sequence", INDICATO, args.primer,
                  "--primer; non e' in testa alle letture")
        elif len(candidati_forward) == 1:
            dichiarato = candidati_forward[0]
            primer_noto = dichiarato
            dell_elenco = r.primer_dell_elenco(dichiarato, v.PRIMER_16S)
            metti("qc.primer_sequence", PROPOSTO, dichiarato,
                  "sequenza che i metadati dichiarano come primer forward; non e' nelle letture"
                  if dichiarato in con_etichetta_forward else
                  f"sequenza dei metadati che coincide con il primer forward {dell_elenco.nome} "
                  f"dell'elenco ({dell_elenco.fonte}): il verso e' quello del primer, anche se i "
                  "metadati non lo dicono; non e' nelle letture",
                  ("sequenze trovate nei metadati: " + ", ".join(sequenze)) if len(sequenze) > 1 else "")
        elif candidati_forward:
            # Due sequenze candidate allo stesso verso: init non sceglie.
            metti("qc.primer_sequence", SCELTA, None,
                  f"partenza: nullo. {len(candidati_forward)} sequenze dei metadati sono "
                  "candidate a primer forward (dichiarate tali, o coincidenti con un primer "
                  "forward dell'elenco) e init non sceglie fra loro: "
                  + ", ".join(candidati_forward),
                  "scrivi qui quella usata, o rilancia init con --primer SEQUENZA: init "
                  "propone allora anche il motivo conservato")
            note.append(
                f"Primer forward non proposto: {len(candidati_forward)} sequenze dei metadati "
                f"sono candidate ({', '.join(candidati_forward)}). qc.primer_sequence parte da "
                "nullo: rilancia init con --primer.")
            non_indovinati.append(
                "qc.primer_sequence: piu' sequenze dei metadati sono candidate a primer forward.")
        else:
            metti("qc.primer_sequence", SCELTA, None,
                  "partenza: nullo, perche' il primer forward non e' nelle letture e i "
                  "metadati non ne riportano la sequenza. Con il primer nullo la validazione "
                  "non cerca alcun primer in testa alle letture. Se il primer e' noto, "
                  "rilancia init con --primer SEQUENZA: init propone allora anche il motivo "
                  "conservato",
                  ("i metadati riportano sequenze senza dire quale sia il forward, e nessuna "
                   "coincide con un primer forward dell'elenco: " + ", ".join(sequenze)
                   + ". Se una e' il primer usato, rilancia init con --primer") if sequenze else "",
                  ("i metadati nominano un primer senza riportarne la sequenza; le forme che "
                   "init conosce con quel nome: "
                   + "; ".join(f"{p.nome} {p.sequenza} ({p.fonte})" for p in nominati)
                   + ". Come scegliere fra piu' forme: RIFERIMENTO.md, 6.3")
                  if nominati else "")
            if sequenze:
                note.append(
                    "Primer forward non proposto: i metadati riportano sequenze ("
                    + ", ".join(sequenze[:4]) + ") che non dichiarano come primer forward e che "
                    "non coincidono con un primer forward dell'elenco. qc.primer_sequence parte "
                    "da nullo: se una e' il primer usato, rilancia init con --primer.")
            non_indovinati.append(
                "qc.primer_sequence: il primer forward non e' nelle letture, e i metadati "
                + ("riportano sequenze di cui non si conosce il verso." if sequenze else
                   "lo nominano senza riportarne la sequenza." if nominati else
                   "non ne riportano la sequenza."))
    motivo = None
    dissentono: list[str] = []
    if primer_noto is not None:
        # I file che non iniziano come la maggioranza escono dal consenso: un
        # file di letture inverse lo allargherebbe fino a comprendere se stesso.
        del_consenso = {
            n: i for n, i in ispezioni.items() if any(i is c for c in da_considerare)}
        if cercato is not None:
            # Con il primer cercato il consenso e' quello delle basi che
            # seguono il primer trovato in ciascuna lettura.
            del_consenso = {
                n: i.dopo_il_primer(cercato.forward.indice) for n, i in del_consenso.items()}
        dissentono = r.file_che_dissentono(del_consenso, taglio, len(primer_noto))
        motivo = r.motivo_dal_consenso(
            [i for n, i in del_consenso.items() if n not in dissentono],
            taglio, len(primer_noto))
    elenco_diversi = ", ".join(dissentono[:5]) + ("..." if len(dissentono) > 5 else "")
    quanti_diversi = (
        "1 file non inizia come gli altri ed e' escluso" if len(dissentono) == 1 else
        f"{len(dissentono)} file non iniziano come gli altri e sono esclusi")
    diversi = (
        f"attenzione: {quanti_diversi} dal consenso ({elenco_diversi}): da guardare prima "
        "di eseguire. Puo' contenere le letture inverse di una corsa appaiata, un'altra "
        "regione, o un campione dominato da sequenze fuori bersaglio"
    ) if dissentono else ""
    if dissentono:
        note.append(
            f"{quanti_diversi} dal consenso del motivo ({elenco_diversi}): da guardare "
            "(letture inverse, un'altra regione, sequenze fuori bersaglio).")
    if motivo is not None:
        metti("qc.conserved_motif", PROPOSTO, motivo,
              f"il consenso delle {len(primer_noto or '')} basi "
              + ("che seguono il primer" if taglio or cercato is not None
                 else "con cui iniziano le letture")
              + f", su non piu' di {r.LETTURE_PER_IL_MOTIVO} letture per file: una base o "
              f"una coppia di basi dove copre almeno il {r.SOGLIA_BASE_DEL_MOTIVO:.0%} delle letture",
              "validate ne misura la frequenza con qc.min_motif_frac", senza_negativi, diversi)
    else:
        metti("qc.conserved_motif", SCELTA, None,
              "partenza: nullo (si verifica la sola assenza del primer). "
              + ("Il consenso delle letture non da' un motivo abbastanza definito"
                 if primer_noto else
                 "Senza conoscere il primer init non propone il motivo della regione amplificata"))
        if primer_noto is None:
            non_indovinati.append(
                "qc.conserved_motif: senza il primer la regione amplificata non e' nota "
                "(nullo e' ammesso).")

    # ---- Profilo: troncamento e funzione d'errore --------------------------------
    profilo: dict[str, Any] | None = None
    if args.skip_profile:
        motivo_senza_profilo = "il profilo di qualita' non e' stato calcolato (--skip-profile)"
    else:
        avanzamento("profilo di qualita' con lo script della fase S1...")
        leggibili = {
            (_chiave_del_file(chiave, n) if chiave else n): p
            for n, p in zip(nomi, file) if ispezioni[n].statistiche.valido
        }
        try:
            profilo = _profilo(
                leggibili, negativi_proposti, processi,
                cercato.per_lo_script() if cercato is not None else None)
            motivo_senza_profilo = ""
        except ErrorePipeline as guasto:
            motivo_senza_profilo = (
                f"il profilo di qualita' non e' stato calcolato ({guasto.codice}: "
                f"{guasto.dettaglio or guasto.voce.sintesi})")
            note.append(
                "Profilo di qualita' non calcolato: " + motivo_senza_profilo
                + ". init va eseguito nel container, come gli altri comandi.")
    suggerito = profilo["suggerito"]["suggerito"] if profilo else None
    if profilo and suggerito is not None:
        s = profilo["suggerito"]
        troncamento = int(suggerito)
        # Letture piu' lunghe della regione amplificata portano in coda il
        # primer dell'altra estremita' e l'adattatore: in single-end nessuna
        # fusione li toglie, e il troncamento deve fermarsi prima.
        coda = None if reverse_di else r.primer_in_coda(ispezioni, cercati_reverse)
        in_coda = ""
        if coda is not None and not coda.chiude_l_amplicone:
            # Lo stesso sito sta anche dentro ampliconi piu' lunghi: senza un
            # adattatore dopo, non si sa se la regione finisce li'.
            note.append(
                f"Le letture contengono il sito del primer {coda.primer.nome} (complemento "
                f"inverso, dopo {coda.prima_posizione} basi nelle letture in cui comincia "
                "prima), ma non e' seguito da un adattatore: se e' il primer reverse usato "
                "per l'amplificazione, il troncamento deve fermarsi prima di quelle basi; se "
                "la regione amplificata prosegue oltre, non conta.")
        elif coda is not None:
            # Con il primer cercato le lunghezze sono quelle dopo il taglio: si
            # toglie cio' che il taglio toglie al piu'.
            prima_del_primer = coda.prima_posizione - (
                cercato.sfasamento + len(cercato.forward.primer.sequenza)
                if cercato is not None else 0)
            in_coda = (
                f"le letture sono piu' lunghe della regione amplificata: nella mediana dei "
                f"file il {coda.mediana:.0%} porta in coda il complemento inverso del primer "
                f"{coda.primer.nome} dell'altra estremita', seguito dall'adattatore, che nelle "
                f"letture in cui comincia prima e' preceduto da {coda.prima_posizione} basi "
                f"(percentile 0,5, su non piu' di {r.LETTURE_PER_IL_MOTIVO} letture per file)")
            if 0 < prima_del_primer < troncamento:
                troncamento = prima_del_primer
                in_coda += (
                    f": il troncamento proposto si ferma li' ({troncamento}), perche' il "
                    "primer e l'adattatore che lo segue non entrino nelle varianti")
            else:
                in_coda += "; il troncamento proposto si ferma gia' prima"
            note.append(
                f"Le letture sono piu' lunghe della regione amplificata (in coda portano il "
                f"primer {coda.primer.nome} dell'altra estremita' e l'adattatore): il "
                "troncamento proposto si ferma prima di quel primer.")
        perdita = _perdita(s) if troncamento == int(suggerito) else _quota_piu_corte(
            profilo, troncamento)
        metti("filter.truncLen", PROPOSTO, troncamento,
              _regola_del_troncamento(s)
              + (": non scarta alcuna lettura" if perdita == 0 else
                 f": scarta {_percento(perdita)} delle letture, nella classe che ne perde di piu'"),
              in_coda, _descrivi_le_lunghezze(profilo, troncamento),
              _descrivi_la_qualita(profilo),
              ("e' la lunghezza della lettura dopo il taglio del primer: le lunghezze sono "
               "misurate sulle letture tagliate" if cercato is not None else
               "comprende le basi di filter.trimLeft"), senza_negativi)
        scarto = max(int(s["senza_perdite"]) - troncamento, 0)
        metti("filter.truncLen_shortfall_warn", SCELTA, min(scarto, troncamento - 1),
              "partenza: lo scarto fra la lettura piu' corta e il troncamento proposto; con "
              "un valore minore il profilo segnalerebbe il troncamento proposto (E-S1-01)")
        _voce_delle_letture_corte(metti, perdita, "il troncamento proposto")
    else:
        metti("filter.truncLen", SCELTA, SENZA_VALORE,
              (motivo_senza_profilo or "le letture non bastano a suggerire un troncamento")
              + ": si sceglie dal profilo delle lunghezze e della qualita'")
        metti("filter.truncLen_shortfall_warn", SCELTA, SENZA_VALORE,
              "si sceglie insieme a filter.truncLen")
        metti("qc.max_frac_short_reads", SCELTA, SENZA_VALORE,
              "si sceglie insieme a filter.truncLen")
    if profilo and len(profilo["valori_di_qualita"]) <= QUALITA_RAGGRUPPATE:
        metti("err.error_function", PROPOSTO, "loess_monotono",
              f"le letture hanno {len(profilo['valori_di_qualita'])} valori di qualita' "
              f"distinti ({', '.join(map(str, profilo['valori_di_qualita']))}): qualita' "
              "raggruppate; il predefinito e' loess (dada2::learnErrors)")

    # ---- Modalita' appaiata: la lettura reverse e la fusione ----------------------
    if reverse_di:
        _voci_delle_coppie(
            args, metti, voci, note, non_indovinati, avanzamento,
            file=file, nomi=nomi, reverse_di=reverse_di, chiave=chiave,
            negativi=negativi_proposti, ispezioni=ispezioni_reverse, cercati=cercati_reverse,
            inserti=inserti, quante=quante, processi=processi,
            tabelle=[righe_assay] + ([righe_studio] if righe_studio is not None else []),
            taglio_forward=taglio,
            primer_forward=(cercato.forward.primer if cercato is not None
                            else osservato.primer if osservato else None),
            profilo_forward=profilo, cercato=cercato,
        )

    # ---- Riferimento ---------------------------------------------------------------
    riferimento: Path | None = args.reference
    if riferimento is not None:
        if not riferimento.is_file():
            raise Rifiuto(f"--reference: {riferimento} non esiste")
        metti("tax.ref_fasta", INDICATO, str(riferimento), "--reference")
        metti("tax.ref_md5", RILEVATO, r.md5_di(riferimento), "calcolato sul file indicato")
        nome, versione = r.nome_e_versione_del_riferimento(riferimento.name)
        for parametro, valore, che_cosa in (("tax.ref_name", nome, "il nome"),
                                            ("tax.ref_version", versione, "la versione")):
            if valore is not None:
                metti(parametro, PROPOSTO, valore, f"letto dal nome del file ({riferimento.name})")
            else:
                metti(parametro, SCELTA, SENZA_VALORE,
                      f"{che_cosa} del riferimento non si legge dal nome del file "
                      f"({riferimento.name})")
    else:
        for parametro in ("tax.ref_name", "tax.ref_version"):
            metti(parametro, SCELTA, SENZA_VALORE, "nessun riferimento indicato (--reference)")

    # ---- Scelte di analisi: valori di partenza -----------------------------------
    letture = [i.statistiche.letture_esaminate for i in ispezioni.values() if i.statistiche.valido]
    mediana = int(statistics.median(letture)) if letture else 0
    # Con ogni file letto per intero il file piu' grande limita il valore: uno
    # maggiore leggerebbe comunque tutte le letture di ogni file.
    intere = args.head_reads is None and bool(letture)
    ispezionate = min(LETTURE_PER_UNA_FRAZIONE, max(letture)) if intere else LETTURE_PER_UNA_FRAZIONE
    metti("qc.head_reads", SCELTA, ispezionate,
          "partenza: "
          + (f"le letture del file piu' grande ({max(letture)}), meno delle "
             f"{LETTURE_PER_UNA_FRAZIONE} con cui" if intere and ispezionate < LETTURE_PER_UNA_FRAZIONE
             else "le letture con cui")
          + " una frazione misurata su un file ha un intervallo di confidenza al "
          f"95% di non oltre {PRECISIONE_DELLE_FRAZIONI:.0%} in piu' o in meno; non e' una scelta "
          "scientifica: con piu' letture la validazione dura di piu' e le frazioni sono "
          "piu' precise",
          (f"il file mediano ha {mediana} letture: un valore pari o superiore legge quel "
           "file per intero") if args.head_reads is None else
          f"init ha ispezionato le prime {args.head_reads} letture di ogni file")
    metti("err.randomize", SCELTA, False,
          "partenza: il predefinito di dada2::learnErrors (randomize = FALSE); true se "
          "l'ordine dei file segue una struttura dello studio")
    metti("dada.pool", SCELTA, False,
          "partenza: il predefinito di dada2::dada (pool = FALSE); \"pseudo\" recupera le "
          "varianti rare presenti in piu' campioni")
    metti("chimera.min_fold_parent_over_abundance", SCELTA, 1.5,
          "partenza: il predefinito del metodo di consenso (dada2::isBimeraDenovoTable, "
          "minFoldParentOverAbundance = 1.5)")
    metti("tax.try_rc", SCELTA, False,
          "partenza: il predefinito di dada2::assignTaxonomy (tryRC = FALSE); true se le "
          "letture possono essere nell'orientamento opposto al riferimento")
    metti("filt.remove_na_phylum", SCELTA, SENZA_VALORE,
          "dipende da che cosa lo studio considera fuori bersaglio")
    metti("filt.exclude_taxa", SCELTA, SENZA_VALORE,
          "i nomi dipendono dallo studio e da come li scrive il riferimento; [] per non "
          "togliere nulla")
    metti("phylo.enabled", SCELTA, False,
          "partenza: disattivata; la filogenesi e' una fase facoltativa, che serve alle "
          "distanze filogenetiche fra campioni")
    metti("phylo.max_seqs", SCELTA, SENZA_VALORE,
          "con phylo.enabled false non ha effetto; altrimenti si sceglie dopo aver visto "
          "nel report quante varianti restano")
    for parametro, origine in (
        ("prev.apply", "se filtrare le varianti presenti in pochi campioni"),
        ("prev.min_fraction", "dipende da quanto sono eterogenei i campioni"),
        ("prev.min_count", "letture minime perche' una variante conti come presente"),
    ):
        metti(parametro, SCELTA, SENZA_VALORE, origine)

    # Controlli positivi e calibrazione.
    metti("ctrl.min_positives", SCELTA, 3,
          "partenza: un punto piu' dei due parametri della curva di calibrazione, il "
          "minimo perche' l'adattamento abbia un residuo"
          + ("" if con_positivi else "; senza controlli positivi non ha effetto"))
    metti("ctrl.min_positive_pass_frac", SCELTA, SENZA_VALORE,
          "quanti controlli positivi devono essere conformi al loro livello"
          + ("" if con_positivi else "; senza controlli positivi non ha effetto"))
    metti("ctrl.positive_gate", SCELTA, False,
          "partenza: un avviso e non un arresto, finche' la conformita' dei controlli "
          "non e' nota dal report della prima esecuzione")
    if con_positivi:
        metti("katharoseq.target_taxon", SCELTA, SENZA_VALORE,
              "il taxon atteso nei controlli positivi: lo sa solo chi li ha allestiti")
        non_indovinati.append(
            "katharoseq.target_taxon: il taxon dei controlli positivi.")
    else:
        metti("katharoseq.target_taxon", SCELTA, None,
              "partenza: nullo, coerente con l'assenza di controlli positivi nella proposta")
    metti("katharoseq.collapse_rank", SCELTA, "Genus",
          "partenza: il rango piu' fine che la pipeline assegna")
    metti("katharoseq.target_sensitivity", SCELTA, SENZA_VALORE,
          "la fedelta' a cui si legge la profondita' minima: una scelta su quanta "
          "contaminazione si accetta nei campioni a bassa biomassa"
          + ("" if con_positivi else "; senza controlli positivi non ha effetto"))
    metti("katharoseq.min_r2", SCELTA, SENZA_VALORE,
          "la bonta' minima di una curva valida"
          + ("" if con_positivi else "; senza controlli positivi non ha effetto"))
    metti("qc.min_reads_mode", SCELTA, "katharoseq_if_available" if con_cellule else "none",
          "partenza: le soglie dei controlli positivi, che la proposta riconosce con la "
          "colonna delle cellule" if con_cellule else
          "partenza: nessuna soglia di profondita', perche' la proposta non ha controlli "
          "positivi con la colonna delle cellule")

    # Decontaminazione.
    metti("decontam.threshold", SCELTA, SENZA_VALORE,
          "la soglia di decontaminazione dipende dalla biomassa: init non la propone",
          "il predefinito di decontam::isContaminant e' 0.1; un valore piu' alto toglie di piu'")
    metti("decontam.min_blanks", SCELTA, SENZA_VALORE,
          "quanti controlli negativi con letture servono per stimare i contaminanti: "
          "dipende dal disegno dei controlli")
    metti("decontam.mode", SCELTA, "aggregate",
          "partenza: tutti i biologici contro tutti i negativi; l'altra modalita' si "
          "calcola comunque come diagnostica")
    if voci["ctrl.blank_values"].valore:
        non_indovinati.append("decontam.threshold: la soglia di decontaminazione.")
    non_indovinati += [
        "Ogni scelta legata alla biomassa o al disegno dello studio: le soglie di qualita' "
        "(qc.*), la prevalenza (prev.*), i taxa da escludere (filt.*), la fedelta' della "
        "calibrazione, la riclassificazione di campioni in controlli negativi.",
    ]

    # Soglie di qualita'.
    memoria = _memoria_della_macchina()
    if memoria is not None:
        massimo = memoria * _QUOTA_DI_MEMORIA / (_BYTE_PER_CELLA * max(len(file), 1))
        if args.head_reads is None and letture:
            massimo = min(massimo, sum(letture))
        metti("qc.max_asv_count", SCELTA, _due_cifre(massimo),
              f"partenza: la tabella delle sequenze piu' grande che occupa il "
              f"{_QUOTA_DI_MEMORIA:.0%} della memoria di questa macchina "
              f"({memoria / 2**30:.0f} GB), con {len(file)} campioni e {_BYTE_PER_CELLA} byte "
              "per cella" + (", non oltre il numero di letture" if args.head_reads is None
                             else "; le letture non sono state contate (--head-reads)"))
    else:
        metti("qc.max_asv_count", SCELTA, SENZA_VALORE, "secondo la memoria disponibile")
    metti("qc.max_zeroed_samples", SCELTA, 0,
          "partenza: nessun campione biologico perso nel filtro senza che lo si veda")
    for parametro, misura in (
        ("qc.max_primer_hit_frac", "frazione di letture con il primer in testa, in "
                                   "01_input_validation/letture_ispezionate.tsv"),
        ("qc.min_motif_frac", "frazione di letture con il motivo, nello stesso file"),
        ("qc.max_frac_lost_filter", "perdita media nel filtro, nel tracciamento delle letture del report"),
        ("qc.warn_frac_chimeric", "frazione di letture chimeriche, nel report"),
        ("qc.stop_frac_chimeric", "frazione di letture chimeriche, nel report"),
        ("qc.min_frac_phylum", "frazione di varianti con il phylum assegnato, nel report"),
        ("qc.max_frac_contaminant", "frazione di letture rimosse dalla decontaminazione, nel report"),
        ("qc.min_reads_final", "letture per campione nell'oggetto finale, nel report"),
        ("qc.min_frac_reads_retained", "letture trattenute dall'insieme, nel report"),
    ):
        metti(parametro, SCELTA, SENZA_VALORE,
              f"soglia di qualita' da fissare sui propri dati; misura da guardare: {misura}")

    mancanti = [c for c in (*defaults.OBBLIGATORI, *defaults.CONDIZIONALI) if c not in voci]
    assert not mancanti, mancanti

    # I condizionali la cui condizione e' falsa con i valori rilevati, proposti
    # o di partenza: non si chiede di sceglierli. Un parametro senza valore
    # (da scegliere, o gia' non richiesto) vale come la sua assenza, e la
    # condizione che ne dipende e' falsa: validate chiedera' il condizionale se
    # l'utente la rende vera. L'ordine dell'elenco mette prima le etichette, da
    # cui dipendono i parametri dei controlli positivi.
    for chiave in defaults.CONDIZIONALI:
        valori = {c: v.valore for c, v in voci.items() if v.valore is not SENZA_VALORE}
        if not defaults.condizione_vera(chiave, valori):
            voci[chiave] = Voce(
                chiave, Specie.CONDIZIONALE, SENZA_VALORE,
                f"serve solo se {defaults.testo_condizione(chiave)}: con i valori di "
                "questo file non e' richiesto e si puo' lasciare in commento")

    intestazione = [
        (f"Letture appaiate: {len(file)} coppie di file in {cartella}, ispezionate {quante} "
         "di ogni file." if reverse_di else
         f"Letture: {len(file)} file in {cartella}, ispezionate {quante} di ogni file."),
        f"Metadati: {assay.name}" + (f" e {studio.name}" if studio else "")
        + (f"; file del lotto {lotto.name}" if lotto else "; nessun file del lotto") + ".",
    ]
    return voci, intestazione, non_indovinati, note


def accetta(voci: Mapping[str, Voce]) -> dict[str, Voce]:
    """Le voci con proposte e valori di partenza accettati in blocco.

    E' il gesto di ``--accetta-proposte``: ogni proposta diventa una riga
    attiva, ogni scelta con un valore di partenza lo prende (``[PARTENZA]``), e
    ogni soglia che ferma o avvisa rimasta senza valore prende quello che non
    ferma (``[ESPLORATIVA]``), se ne ha uno. Nessun valore e' inventato: una
    scelta senza valore di partenza resta da scegliere, un condizionale non
    pertinente resta omesso.
    """
    accettate = {}
    valori = {c: x.valore for c, x in voci.items() if x.valore is not SENZA_VALORE}
    esplorative = {s.chiave: s for s in SOGLIE}
    for chiave, voce in voci.items():
        if voce.specie is PROPOSTO:
            voce = replace(voce, accettata=True,
                           origine=voce.origine + "; accettato con --accetta-proposte")
        elif voce.specie is SCELTA and voce.valore is not SENZA_VALORE:
            voce = replace(voce, specie=Specie.PARTENZA,
                           origine=voce.origine.removeprefix("partenza: ")
                           + "; accettato con --accetta-proposte")
        elif voce.specie is SCELTA and chiave in esplorative:
            soglia = esplorative[chiave]
            valore = valore_esplorativo(soglia)
            if valore is not None:
                senza_effetto = non_eseguito(soglia, {**valori, chiave: valore})
                voce = replace(
                    voce, specie=Specie.ESPLORATIVA, valore=valore,
                    origine=(
                        "il valore che non ferma: il controllo "
                        f"({', '.join(soglia.codici)}) e' spento. Serve alla prima "
                        "esecuzione; poi il sottocomando soglie propone il valore dalla "
                        "misura" if senza_effetto is None else
                        f"il valore che non ferma; {senza_effetto}"),
                    note=(voce.origine, *voce.note))
        accettate[chiave] = voce
    return accettate


#: Quanti campioni si nominano per ogni valore di una colonna, nel commento.
_CAMPIONI_NOMINATI: Final = 6


def _campioni_per_valore(
    righe: list[dict[str, str]], colonna: str, colonna_dei_nomi: str | None
) -> str:
    """Per ogni valore della colonna, i campioni che lo portano: i primi
    :data:`_CAMPIONI_NOMINATI` per nome, e quanti sono gli altri. Chi legge il
    file vede cosi' quali campioni ogni etichetta raccoglierebbe.
    """
    per_valore: dict[str, list[str]] = {}
    for riga in righe:
        per_valore.setdefault(riga.get(colonna, ""), []).append(
            riga.get(colonna_dei_nomi, "") if colonna_dei_nomi else "")
    parti = []
    for valore in sorted(per_valore):
        campioni = sorted(per_valore[valore])
        if not colonna_dei_nomi:
            parti.append(f"{valore or '<vuoto>'!r}: {len(campioni)} campioni")
            continue
        resto = len(campioni) - _CAMPIONI_NOMINATI
        parti.append(
            f"{valore or '<vuoto>'!r}: {', '.join(campioni[:_CAMPIONI_NOMINATI])}"
            + (f" e altri {resto}" if resto > 0 else ""))
    return "; ".join(parti)


def _raggruppamento(righe: Sequence[dict[str, str]], colonna: str) -> tuple[int, ...]:
    """Come una colonna raggruppa le righe, senza i nomi dei gruppi: due colonne
    con lo stesso raggruppamento dividono i campioni allo stesso modo.
    """
    indice: dict[str, int] = {}
    return tuple(indice.setdefault(riga.get(colonna, ""), len(indice)) for riga in righe)


def _classi_dei_campioni(
    righe: list[dict[str, str]], tipo: r.TipoDiCampione
) -> list[str]:
    """La classe che una colonna candidata da' a ogni riga, nell'ordine delle
    righe: due colonne che danno le stesse classi dicono la stessa cosa, anche
    con valori diversi.
    """
    def classe(valore: str) -> str:
        return ("negativo" if valore in tipo.negativi else
                "positivo" if valore in tipo.positivi else
                "non attribuito" if valore in tipo.non_attribuiti else "biologico")
    return [classe(riga.get(tipo.colonna, "")) for riga in righe]


def _chiave_del_file(chiave: r.Chiave, nome: str) -> str | None:
    """La chiave del campione a cui appartiene un file."""
    for c, n in chiave.file_per_chiave.items():
        if n == nome:
            return c
    return None


def _chiave_della_riga(chiave: r.Chiave, riga: dict[str, str]) -> str | None:
    """La chiave del campione di una riga della tabella di assay."""
    import re

    from amplicon16s.metadata.crosswalk import AccessionNonEstraibile, chiave_dalla_tabella

    try:
        return chiave_dalla_tabella(riga.get(chiave.colonna, ""), re.compile(chiave.espressione))
    except AccessionNonEstraibile:
        return None


def _letture_appaiate(dettaglio: str) -> str:
    """Il messaggio con cui ``init`` rifiuta una cartella in cui alcuni file
    contengono seconde letture, ma i file non si abbinano a coppie per nome.
    """
    return (
        f"le letture risultano appaiate: {dettaglio}. I file non si abbinano pero' a coppie "
        "per nome: le due letture di un campione si riconoscono solo se i due nomi "
        "differiscono per il marcatore in fondo al nome, prima dell'estensione (_R1 e _R2, _1 "
        "e _2, anche seguito da _001), e non se dopo il marcatore il nome continua (per "
        "esempio _R1_raw). Usa le sole letture forward: una cartella a parte con un file per "
        "campione, quello della prima lettura di ogni coppia, indicata con --reads (GUIDA.md, "
        "passo 2)"
    )


def _senza_compagno(nomi: list[str], inverse: list[str]) -> str:
    """Il messaggio con cui ``init`` rifiuta una cartella in cui alcuni file
    hanno il nome di una seconda lettura ma non ogni file ha il suo compagno:
    nomina i file rimasti soli e le due strade possibili.
    """
    soli: list[str] = []
    for marcatori in MARCATORI_DI_FILE:
        abbinamento = abbina_file(nomi, marcatori)
        if abbinamento.coppie:
            soli = [*abbinamento.senza_compagno, *abbinamento.ambigui]
            break
    elenco = soli or inverse
    return (
        f"{len(inverse)} file hanno il nome della seconda lettura di una coppia, ma non ogni "
        f"file della cartella ha il suo compagno: {', '.join(elenco[:5])}"
        f"{'...' if len(elenco) > 5 else ''}"
        + (" (senza compagno)" if soli else "")
        + ". Per la modalita' appaiata aggiungi i file che mancano o togli dalla cartella "
        "quelli che non appartengono a una coppia; per usare le sole letture forward tieni "
        "nella cartella indicata con --reads un file per campione, con la prima lettura di "
        "ogni coppia (RIFERIMENTO.md, 1.2)"
    )


def esegui(args: argparse.Namespace, stampa: Callable[[str], None] = print) -> int:
    """Esegue ``init`` e restituisce il codice di uscita."""
    uscita: Path = args.output
    if uscita.exists():
        stampa(f"init: {uscita} esiste gia': init non scrive su una configurazione "
               "esistente. Indica con --output un file nuovo")
        return USCITA_RIFIUTO
    if not uscita.parent.is_dir():
        stampa(f"init: la cartella {uscita.parent} non esiste")
        return USCITA_RIFIUTO
    if not os.access(uscita.parent, os.W_OK | os.X_OK):
        # Si verifica prima di leggere le letture: una cartella montata nel
        # container e creata da Docker appartiene a un altro utente.
        stampa(f"init: la cartella {uscita.parent} non e' scrivibile dall'utente che esegue: "
               "creala prima di lanciare il comando, oppure indica un altro --output. "
               "Nulla e' stato scritto")
        return USCITA_RIFIUTO
    try:
        voci, intestazione, non_indovinati, note = prepara(args, lambda t: stampa(f"init: {t}"))
    except Rifiuto as rifiuto:
        stampa(f"init: {rifiuto}. Nulla e' stato scritto")
        return USCITA_RIFIUTO
    accettate = bool(getattr(args, "accetta_proposte", False))
    if accettate:
        voci = accetta(voci)
    scrivi_atomico(uscita, testo_yaml(voci, intestazione, non_indovinati, accettate) + "\n")

    da_schema = [
        f"{gruppo}.{nome}" for gruppo, campo in Config.model_fields.items()
        for nome in campo.annotation.model_fields
    ]
    voci = {chiave: voci[chiave] for chiave in da_schema if chiave in voci}
    per_specie = {s: [x for x in voci.values() if x.specie is s] for s in Specie}
    stampa(f"Configurazione di partenza scritta in {uscita}")
    for riga in intestazione:
        stampa(f"  {riga}")
    stampa(f"  rilevati: {len(per_specie[Specie.RILEVATO])}, indicati: "
           f"{len(per_specie[Specie.INDICATO])}, proposti: {len(per_specie[Specie.PROPOSTO])}, "
           f"da scegliere: {len(per_specie[Specie.DA_SCEGLIERE])}, "
           f"non richiesti: {len(per_specie[Specie.CONDIZIONALE])}")
    stampa("Rilevati sui file (righe attive):")
    for voce in per_specie[Specie.RILEVATO]:
        stampa(f"  {voce.chiave}: {_breve(voce.valore)}")
    stampa("Proposti (righe attive, accettati con --accetta-proposte: controllali):"
           if accettate else
           "Proposti (righe in commento: togli il commento per accettare):")
    for voce in per_specie[Specie.PROPOSTO]:
        stampa(f"  {voce.chiave}: {_breve(voce.valore)}")
    if accettate:
        stampa("Valori di partenza accettati (righe attive):")
        for voce in per_specie[Specie.PARTENZA]:
            stampa(f"  {voce.chiave}: {_breve(voce.valore)}")
        assenze = [x.chiave for x in per_specie[Specie.PARTENZA] if x.valore in (None, [])]
        if assenze:
            stampa("  Fra questi dichiarano un'assenza (nullo o elenco vuoto), da rivedere se "
                   "lo studio ha quell'informazione: " + ", ".join(assenze))
        stampa("Soglie esplorative, al valore che non ferma (righe attive):")
        for voce in per_specie[Specie.ESPLORATIVA]:
            stampa(f"  {voce.chiave}: {_breve(voce.valore)}")
    stampa("Restano da scegliere (righe in commento, da scrivere):" if accettate else
           "Da scegliere:")
    for voce in per_specie[Specie.DA_SCEGLIERE]:
        partenza = "" if voce.valore is SENZA_VALORE else f" (partenza: {_breve(voce.valore)})"
        stampa(f"  {voce.chiave}{partenza}")
    if per_specie[Specie.CONDIZIONALE]:
        stampa("Non richiesti con le proposte attuali (condizionali):")
        for voce in per_specie[Specie.CONDIZIONALE]:
            stampa(f"  {voce.chiave}: {voce.origine.split(':')[0]}")
    if non_indovinati:
        stampa("Che cosa init non stabilisce:")
        for testo in non_indovinati:
            stampa(f"  - {testo}")
    for testo in note:
        stampa(f"Nota: {testo}")
    stampa(f"Poi: rivedi {uscita}, e lancia validate: elenca i parametri ancora da "
           "dichiarare (E-G15-10).")
    if accettate:
        stampa("Dopo la prima esecuzione, il sottocomando soglie propone dalle misure le "
               "soglie ora esplorative.")
    return 0


def _breve(valore: Any) -> str:
    """Un valore su una riga, accorciato se lungo."""
    from amplicon16s.init.scrittura import valore_yaml

    testo = valore_yaml(valore)
    return testo if len(testo) <= 90 else testo[:87] + "..."
