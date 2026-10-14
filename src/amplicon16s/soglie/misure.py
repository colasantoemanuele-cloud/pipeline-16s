"""Le misure di un'esecuzione conclusa che i controlli confrontano con le soglie.

Per ogni soglia di :data:`~amplicon16s.soglie.tabella.SOGLIE` si legge, dagli
artefatti che l'esecuzione ha lasciato nella cartella di uscita, **la stessa
statistica che il controllo confronta con la soglia**: il valore peggiore fra i
file dove il controllo giudica file per file, la mediana dove giudica la
mediana, la frazione d'insieme dove giudica l'insieme. La cartella non viene
modificata, e nulla viene ricalcolato dalle letture, con una sola eccezione: le
letture piu' corte del troncamento fra le prime lette dalla validazione, che
nessun artefatto riporta quando il file non e' stato letto per intero, si
contano rileggendo quelle prime letture dal file FASTQ, con la stessa funzione
della validazione.

Dove l'artefatto riporta la statistica arrotondata e il controllo confronta il
valore non arrotondato, la misura porta una **tolleranza di lettura**: meta'
dell'ultima cifra scritta. Chi propone una soglia peggiora la misura di quella
tolleranza, cosi' che il suggerimento non sia mai piu' stretto del valore vero.

Una misura puo' mancare per due ragioni distinte, che chi legge deve poter
distinguere: il controllo non ha prodotto una misura in questa esecuzione
(``assente``), oppure la misura esiste ma gli artefatti non bastano a
ricostruirla (``non_leggibile``).

**Chi lo usa.** Il sottocomando ``soglie``
(:mod:`amplicon16s.soglie.comando`); la riga di comando, per il numero dei
campioni controllati con cui riconosce i controlli spenti.
"""

from __future__ import annotations

import csv
import json
import statistics
from collections import defaultdict
from collections.abc import Mapping
from dataclasses import dataclass
from decimal import Decimal
from fractions import Fraction
from pathlib import Path
from typing import Any, Final

from amplicon16s.io_layer.artifacts import Fase, nome_manifesto_passo
from amplicon16s.io_layer.reads import RicercaPrimer, scansiona_file
from amplicon16s.metadata.lettura_inventario import NOME_CROSSWALK
from amplicon16s.metadata.models import CLASSI_CONTROLLATE
from amplicon16s.metadata.tabelle import COLONNE_INVENTARIO
from amplicon16s.soglie.tabella import non_eseguito, soglia
from amplicon16s.steps.s00_validate import NOME_INVENTARIO, NOME_LETTURE, NOME_LETTURE_REVERSE
from amplicon16s.steps.s01_profile import NOME_LUNGHEZZE, nome_reverse
from amplicon16s.steps.s01_profile import NOME_RIEPILOGO as RIEPILOGO_PROFILO
from amplicon16s.steps.s04m_merge import NOME_FUSIONE
from amplicon16s.steps.s04m_merge import NOME_RIEPILOGO as RIEPILOGO_FUSIONE
from amplicon16s.steps.s05_seqtab import NOME_RIEPILOGO as RIEPILOGO_TABELLA
from amplicon16s.steps.s06_chimera import NOME_RIEPILOGO as RIEPILOGO_CHIMERE
from amplicon16s.steps.s08_taxonomy import NOME_RIEPILOGO as RIEPILOGO_TASSONOMIA
from amplicon16s.steps.s11_controls import NOME_RIEPILOGO as RIEPILOGO_CONTROLLI
from amplicon16s.steps.s12_decontam import NOME_RIEPILOGO as RIEPILOGO_DECONTAM
from amplicon16s.steps.s13_filtri import NOME_RIEPILOGO as RIEPILOGO_FILTRI

__all__ = ["ArtefattoNonLeggibile", "LettureFinali", "Misura", "campioni_controllati",
           "gruppi_dei_campioni_finali", "leggi", "letture_finali"]

_CONTROLLATE: Final = tuple(c.value for c in CLASSI_CONTROLLATE)
_INTERMEDI: Final = "intermedi"


class ArtefattoNonLeggibile(Exception):
    """Un artefatto da cui si legge una misura manca o non ha la forma attesa."""


@dataclass(frozen=True)
class Misura:
    """La statistica che un controllo confronta con la sua soglia."""

    #: Il valore letto, o ``None`` se la misura manca.
    valore: Fraction | None
    #: Che cosa e' stato misurato, a parole.
    statistica: str
    #: L'artefatto da cui e' letta, relativo alla cartella di uscita.
    origine: str
    #: Meta' dell'ultima cifra scritta, se il controllo confronta il valore
    #: non arrotondato; zero se il valore letto e' quello confrontato.
    tolleranza: Fraction = Fraction(0)
    #: Perche' il controllo non ha una misura in questa esecuzione.
    assente: str | None = None
    #: Perche' la misura esiste ma non si ricostruisce dagli artefatti.
    non_leggibile: str | None = None


def _json(radice: Path, relativo: str) -> Any:
    """Il contenuto di un artefatto JSON, con i numeri reali come decimali esatti."""
    try:
        return json.loads((radice / relativo).read_text(encoding="utf-8"), parse_float=Decimal)
    except (OSError, ValueError) as guasto:
        raise ArtefattoNonLeggibile(f"{relativo}: {guasto}") from guasto


def _tsv(radice: Path, relativo: str) -> list[dict[str, str]]:
    """Le righe di un artefatto tabellare."""
    try:
        with open(radice / relativo, encoding="utf-8", newline="") as file:
            return list(csv.DictReader(file, delimiter="\t"))
    except OSError as guasto:
        raise ArtefattoNonLeggibile(f"{relativo}: {guasto}") from guasto


def _campo(documento: Any, relativo: str, *percorso: str) -> Any:
    """Un campo annidato di un documento, o un errore che dice quale manca."""
    valore = documento
    for chiave in percorso:
        if not isinstance(valore, dict) or chiave not in valore:
            raise ArtefattoNonLeggibile(f"{relativo}: manca il campo {'.'.join(percorso)}")
        valore = valore[chiave]
    return valore


def _frazione(valore: Any) -> Fraction:
    """Un numero letto da un artefatto, come frazione esatta."""
    return Fraction(Decimal(str(valore)))


def _frazione_delle_letture(scritta: str, letture: int) -> tuple[Fraction, Fraction]:
    """La frazione di un file e la sua tolleranza di lettura.

    L'artefatto la scrive con sei decimali, ma e' un conteggio diviso per le
    letture esaminate: se un solo conteggio da' quella scrittura la frazione
    si ricostruisce esatta, altrimenti resta quella scritta, con la tolleranza
    di meta' dell'ultima cifra.
    """
    frazione = _frazione(scritta)
    conteggio = round(frazione * letture)
    if letture and f"{conteggio / letture:.6f}" == scritta.strip():
        vicini = [c for c in (conteggio - 1, conteggio + 1)
                  if 0 <= c <= letture and f"{c / letture:.6f}" == scritta.strip()]
        if not vicini:
            return Fraction(conteggio, letture), Fraction(0)
    return frazione, Fraction(1, 2_000_000)


def _classi(radice: Path) -> dict[str, str]:
    """La classe di ogni campione, dalla tabella di S0."""
    relativo = f"{Fase.INPUT_VALIDATION.value}/{NOME_CROSSWALK}"
    try:
        return {r["accession"]: r["classe"] for r in _tsv(radice, relativo)}
    except KeyError as guasto:
        raise ArtefattoNonLeggibile(f"{relativo}: manca la colonna {guasto}") from guasto


@dataclass(frozen=True)
class LettureFinali:
    """Le letture di un campione dell'oggetto finale."""

    accession: str
    nome: str
    letture: int


#: Le tabelle della fase dei filtri finali: letture per campione ed esclusi.
_LETTURE_FINALI: Final = "letture_finali.tsv"
_ESCLUSIONI: Final = "esclusioni.tsv"


def letture_finali(radice: Path | str) -> list[LettureFinali] | None:
    """Le letture finali dei campioni rimasti nell'oggetto finale, nell'ordine
    della tabella; ``None`` se le tabelle dei filtri finali non ci sono.

    Sono i campioni biologici della tabella delle letture finali che la
    tabella delle esclusioni non nomina: i soli su cui ha senso chiedersi se
    qualcuno ha molte meno letture degli altri. I controlli compaiono nella
    tabella con zero letture finali, perche' l'oggetto finale non li contiene:
    contarli farebbe credere che l'oggetto abbia campioni vuoti.
    """
    radice = Path(radice)
    cartella = f"{Fase.FINAL.value}/{_INTERMEDI}"
    if not all((radice / cartella / n).is_file() for n in (_LETTURE_FINALI, _ESCLUSIONI)):
        return None
    try:
        esclusi = {r["accession"] for r in _tsv(radice, f"{cartella}/{_ESCLUSIONI}")}
        relativo = f"{Fase.INPUT_VALIDATION.value}/{NOME_CROSSWALK}"
        inventario = _tsv(radice, relativo)
        nomi = {r["accession"]: r.get("sample_name", "") for r in inventario}
        controlli = {r["accession"] for r in inventario
                     if r.get("classe", "biologico") != "biologico"}
        return [
            LettureFinali(r["campione"], nomi.get(r["campione"], ""), int(r["letture"]))
            for r in _tsv(radice, f"{cartella}/{_LETTURE_FINALI}")
            if r["campione"] not in esclusi and r["campione"] not in controlli
        ]
    except (KeyError, ValueError, ArtefattoNonLeggibile):
        return None


def gruppi_dei_campioni_finali(radice: Path | str) -> dict[str, dict[str, str]]:
    """Per ogni colonna dello studio portata nell'oggetto finale, il valore di
    ogni campione: dalla tabella dei campioni consegnata in testo. Vuoto se la
    tabella non c'e' o non porta colonne oltre a quelle dell'inventario.
    """
    relativo = f"{Fase.FINAL.value}/metadati.tsv"
    if not (Path(radice) / relativo).is_file():
        return {}
    try:
        righe = _tsv(Path(radice), relativo)
    except ArtefattoNonLeggibile:
        return {}
    colonne = [c for c in (righe[0] if righe else {}) if c not in COLONNE_INVENTARIO]
    return {c: {r.get("accession", ""): r.get(c, "") for r in righe} for c in colonne}


def campioni_controllati(radice: Path | str) -> int | None:
    """Quanti campioni biologici e controlli positivi ha l'inventario di S0, o
    ``None`` se S0 non l'ha ancora scritto.
    """
    percorso = Path(radice) / Fase.INPUT_VALIDATION.value / NOME_INVENTARIO
    try:
        conteggi = json.loads(percorso.read_text(encoding="utf-8"))["conteggi"]
        return sum(int(conteggi.get(classe, 0)) for classe in _CONTROLLATE)
    except (OSError, ValueError, KeyError, TypeError, AttributeError):
        return None


def _letture_ispezionate(
    radice: Path, nome: str = NOME_LETTURE
) -> tuple[str, list[dict[str, str]]]:
    """Le righe dei file leggibili dei campioni controllati, come le giudicano
    i gate che leggono le sequenze. Con il nome della tabella dei file reverse,
    quelle della lettura reverse.
    """
    relativo = f"{Fase.INPUT_VALIDATION.value}/{nome}"
    classi = _classi(radice)
    righe = [
        r for r in _tsv(radice, relativo)
        if classi.get(r.get("accession", "")) in _CONTROLLATE and not r.get("errore")
    ]
    return relativo, righe


def _primer(radice: Path, valori: Mapping[str, Any], reverse: bool = False) -> Misura:
    """La misura per la soglia sulla frazione di letture che iniziano con il
    primer: della lettura forward, o della reverse con ``reverse``.

    Passi: legge le righe dei file ispezionati delle classi controllate; se il
    controllo non e' stato eseguito con questi valori, o nessun file e' stato
    ispezionato, la misura e' assente con il motivo; altrimenti e' il massimo
    delle frazioni per file. Il massimo perche' la soglia e' un tetto che il
    controllo applica file per file; la tolleranza di lettura e' la piu' larga
    fra quelle dei file. Una tabella che non si interpreta solleva
    :class:`ArtefattoNonLeggibile`.
    """
    statistica = ("massimo, fra i file dei campioni biologici e dei controlli positivi, della "
                  "frazione delle prime letture che iniziano con il primer")
    if reverse:
        statistica = statistica.replace("fra i file", "fra i file reverse").replace(
            "il primer", "il primer reverse")
    relativo, righe = _letture_ispezionate(
        radice, NOME_LETTURE_REVERSE if reverse else NOME_LETTURE)
    motivo = non_eseguito(
        soglia("pair.max_primer_hit_frac_reverse" if reverse else "qc.max_primer_hit_frac"),
        valori)
    if motivo is not None:
        return Misura(None, statistica, relativo, assente=motivo)
    if not righe:
        return Misura(None, statistica, relativo,
                      assente="nessun file di un campione biologico o di un controllo positivo "
                              "e' stato ispezionato")
    try:
        lette = [_frazione_delle_letture(r["frazione_primer"], int(r["letture_esaminate"]))
                 for r in righe]
    except (KeyError, ValueError, ArithmeticError) as guasto:
        raise ArtefattoNonLeggibile(f"{relativo}: frazione_primer: {guasto}") from guasto
    return Misura(max(f for f, _ in lette), statistica, relativo,
                  tolleranza=max(t for _, t in lette))


def _motivo(radice: Path, valori: Mapping[str, Any]) -> Misura:
    """La misura per la soglia sulla frazione di letture con il motivo
    conservato.

    Stessi passi della misura del primer, ma la statistica e' la mediana fra
    i file delle classi controllate: e' la mediana che il controllo confronta
    con la soglia. La tolleranza di lettura e' la piu' larga fra quelle dei
    file.
    """
    statistica = ("mediana, fra i file dei campioni biologici e dei controlli positivi, della "
                  "frazione delle prime letture con il motivo conservato")
    relativo, righe = _letture_ispezionate(radice)
    motivo = non_eseguito(soglia("qc.min_motif_frac"), valori)
    if motivo is not None:
        return Misura(None, statistica, relativo, assente=motivo)
    if not righe:
        return Misura(None, statistica, relativo,
                      assente="nessun file di un campione biologico o di un controllo positivo "
                              "e' stato ispezionato")
    try:
        lette = [_frazione_delle_letture(r["frazione_motivo"], int(r["letture_esaminate"]))
                 for r in righe]
    except (KeyError, ValueError, ArithmeticError) as guasto:
        raise ArtefattoNonLeggibile(f"{relativo}: frazione_motivo: {guasto}") from guasto
    return Misura(Fraction(statistics.median(f for f, _ in lette)), statistica, relativo,
                  tolleranza=max(t for _, t in lette))


def _letture_corte(radice: Path, valori: Mapping[str, Any]) -> Misura:
    """La frazione di letture piu' corte del troncamento, per classe: il
    massimo fra quella di S1 (tutte le letture) e quella della validazione
    (le prime letture), che non e' scritta e si ricostruisce solo se ogni file
    e' senza letture corte fra le prime o e' stato letto per intero.
    """
    statistica = ("massimo, fra le classi biologico e controllo positivo, della frazione di "
                  "letture piu' corte di filter.truncLen (su tutte le letture, e sulle prime "
                  "lette dalla validazione)")
    manifesto = f"{Fase.QC_PROFILES.value}/{nome_manifesto_passo('S1')}"
    metriche = _campo(_json(radice, manifesto), manifesto, "metriche")
    forward = _corte_di_una_lettura(
        radice, manifesto, _campo(metriche, manifesto, "letture_piu_corte_del_troncamento"),
        int(valori["filter.truncLen"]), NOME_LETTURE, NOME_LUNGHEZZE, statistica, valori)
    if valori.get("io.layout") != "paired-end" or forward.non_leggibile is not None:
        return forward
    # In modalita' appaiata la stessa soglia giudica anche le letture reverse,
    # con il loro troncamento: la misura e' la peggiore fra le due letture.
    statistica = statistica.replace(
        "filter.truncLen", "filter.truncLen o, per le reverse, di pair.truncLen_reverse")
    reverse = _corte_di_una_lettura(
        radice, manifesto,
        _campo(metriche, manifesto, "reverse", "letture_piu_corte_del_troncamento"),
        int(valori["pair.truncLen_reverse"]), NOME_LETTURE_REVERSE, nome_reverse(NOME_LUNGHEZZE),
        statistica, valori, reverse=True)
    origine = f"{forward.origine}; {Fase.INPUT_VALIDATION.value}/{NOME_LETTURE_REVERSE}"
    if reverse.non_leggibile is not None:
        return Misura(None, statistica, origine, non_leggibile=reverse.non_leggibile)
    presenti = [m.valore for m in (forward, reverse) if m.valore is not None]
    if not presenti:
        return Misura(None, statistica, origine, assente=forward.assente)
    return Misura(max(presenti), statistica, origine)


def _corte_fra_le_prime(
    valori: Mapping[str, Any], nome_file: str, esaminate: int, troncamento: int, reverse: bool
) -> int | None:
    """Quante delle prime ``esaminate`` letture di un file sono piu' corte del
    troncamento, contate sul file FASTQ con la funzione della validazione;
    ``None`` se il file non e' piu' quello ispezionato (non si trova, non si
    legge, o le sue prime letture non sono tante quante allora).

    La tabella della validazione riporta la lettura piu' corta fra le prime ma
    non quante sono corte: quando il file non e' stato letto per intero il
    numero non e' in alcun artefatto, e si ricava solo dalle letture.
    """
    ricerca = None
    if valori.get("filter.primer_removal") == "search":
        ricerca = RicercaPrimer(
            valori["pair.primer_sequence_reverse" if reverse else "qc.primer_sequence"],
            valori["filter.primer_max_mismatch"], valori["filter.primer_max_offset"],
            scarta_senza_primer=valori["filter.primer_untrimmed"] == "discard")
    statistiche = scansiona_file(
        Path(str(valori["io.fastq_dir"])) / nome_file, esaminate, None, None,
        lunghezza_richiesta=troncamento, ricerca=ricerca)
    if not statistiche.valido or statistiche.letture_esaminate != esaminate:
        return None
    return statistiche.piu_corte


def _corte_di_una_lettura(
    radice: Path, manifesto: str, corte: Any, troncamento: int, nome_letture: str,
    nome_lunghezze: str, statistica: str, valori: Mapping[str, Any], reverse: bool = False,
) -> Misura:
    """La frazione di letture corte di una delle due letture (o della sola
    lettura, in modalita' single-end): dal manifesto di S1 e dalla tabella
    delle letture ispezionate dalla validazione; per un file con letture corte
    fra le prime e non letto per intero, dal file delle letture.
    """
    del_profilo = [
        _frazione(_campo(corte, manifesto, classe, "frazione"))
        for classe in _CONTROLLATE if classe in corte
    ]
    origine = f"{manifesto}; {Fase.INPUT_VALIDATION.value}/{nome_letture}"

    relativo, righe = _letture_ispezionate(radice, nome_letture)
    classi = _classi(radice)
    lunghezze: dict[str, int] | None = None
    dai_file = False
    per_classe: dict[str, list[int]] = defaultdict(lambda: [0, 0])
    for riga in righe:
        try:
            esaminate = int(riga["letture_esaminate"])
            minima = int(riga["lunghezza_minima"]) if riga["lunghezza_minima"] else None
            esaurito = riga["file_esaurito"] == "si"
        except (KeyError, ValueError) as guasto:
            raise ArtefattoNonLeggibile(f"{relativo}: {guasto}") from guasto
        if minima is None or minima >= troncamento:
            corte_del_file = 0
        elif esaurito:
            # Il file e' stato letto per intero: le prime letture sono tutte le
            # letture, e quelle corte sono nella tabella delle lunghezze di S1.
            if lunghezze is None:
                lunghezze = _corte_per_campione(radice, troncamento, nome_lunghezze)
            corte_del_file = lunghezze.get(riga["accession"], 0)
        else:
            # Quante siano non e' scritto in alcun artefatto: si contano
            # rileggendo le stesse prime letture dal file.
            contate = _corte_fra_le_prime(
                valori, riga.get("file", ""), esaminate, troncamento, reverse)
            if contate is None:
                return Misura(
                    None, statistica, origine,
                    non_leggibile=(
                        f"il file {riga.get('file', riga['accession'])} ha letture piu' corte "
                        f"del troncamento fra le prime {esaminate} e non e' stato letto per "
                        "intero dalla validazione: quante siano non e' scritto in alcun "
                        "artefatto, e il file non e' piu' leggibile in io.fastq_dir com'era "
                        "alla validazione"))
            corte_del_file = contate
            dai_file = True
        voce = per_classe[classi[riga["accession"]]]
        voce[0] += corte_del_file
        voce[1] += esaminate
    della_validazione = [Fraction(c, n) for c, n in per_classe.values() if n]
    tutte = del_profilo + della_validazione
    if dai_file:
        origine += "; le prime letture dei file in io.fastq_dir"
    if not tutte:
        return Misura(None, statistica, origine,
                      assente="nessuna lettura di un campione biologico o di un controllo positivo")
    return Misura(max(tutte), statistica, origine)


def _corte_per_campione(
    radice: Path, troncamento: int, nome: str = NOME_LUNGHEZZE
) -> dict[str, int]:
    """Per campione, le letture piu' corte del troncamento secondo il profilo di S1."""
    relativo = f"{Fase.QC_PROFILES.value}/{nome}"
    corte: dict[str, int] = defaultdict(int)
    try:
        for riga in _tsv(radice, relativo):
            if int(riga["lunghezza"]) < troncamento:
                corte[riga["campione"]] += int(float(riga["letture"]))
    except (KeyError, ValueError) as guasto:
        raise ArtefattoNonLeggibile(f"{relativo}: {guasto}") from guasto
    return dict(corte)


def _scarto(radice: Path, valori: Mapping[str, Any], reverse: bool = False) -> Misura:
    """Le basi fra la lettura piu' corta, su tutte le letture, e il troncamento
    dichiarato (zero se il troncamento e' piu' lungo); con ``reverse`` per la
    lettura reverse. Dal riepilogo dei profili."""
    riepilogo = nome_reverse(RIEPILOGO_PROFILO) if reverse else RIEPILOGO_PROFILO
    parametro = "pair.truncLen_reverse" if reverse else "filter.truncLen"
    relativo = f"{Fase.QC_PROFILES.value}/{riepilogo}"
    minima = int(_campo(_json(radice, relativo), relativo, "lunghezza_minima"))
    return Misura(
        Fraction(max(minima - int(valori[parametro]), 0)),
        f"basi fra la lettura{' reverse' if reverse else ''} piu' corta, su tutte le letture, "
        f"e {parametro}", relativo)


def _fusione(radice: Path) -> Misura:
    """La quota di coppie fuse sulle coppie delle classi controllate, dal
    riepilogo della fusione; assente se nessuna coppia e' entrata nella fusione."""
    relativo = f"{Fase.PAIR_MERGE.value}/{RIEPILOGO_FUSIONE}"
    statistica = ("quota di coppie fuse, sulle coppie dei campioni biologici e dei controlli "
                  "positivi entrate nella fusione")
    quota = _campo(_json(radice, relativo), relativo, "quota_controllata")
    if quota is None:
        return Misura(None, statistica, relativo,
                      assente="nessuna coppia di un campione biologico o di un controllo positivo")
    return Misura(_frazione(quota), statistica, relativo)


def _fusione_per_campione(radice: Path) -> Misura:
    """La quota di coppie fuse piu' bassa fra i campioni biologici e i
    controlli positivi con coppie, dai conteggi per campione della fusione.
    """
    relativo = f"{Fase.PAIR_MERGE.value}/{NOME_FUSIONE}"
    statistica = ("minimo, fra i campioni biologici e i controlli positivi con coppie, della "
                  "quota di coppie fuse del campione")
    classi = _classi(radice)
    try:
        quote = [
            Fraction(int(float(r["fuse"])), int(float(r["coppie"])))
            for r in _tsv(radice, relativo)
            if classi.get(r["campione"]) in _CONTROLLATE and int(float(r["coppie"]))
        ]
    except (KeyError, ValueError) as guasto:
        raise ArtefattoNonLeggibile(f"{relativo}: {guasto}") from guasto
    if not quote:
        return Misura(None, statistica, relativo,
                      assente="nessun campione biologico o controllo positivo con coppie")
    return Misura(min(quote), statistica, relativo)


def _filtro(radice: Path) -> tuple[Misura, Misura]:
    """Le due misure del filtro, dal manifesto della fase: i campioni delle
    classi controllate rimasti senza letture, e la perdita media di letture.

    La perdita e' assente se nessun campione controllato aveva letture. La sua
    tolleranza e' meta' dell'ultima cifra con cui il manifesto la riporta.
    """
    relativo = f"{Fase.FILTERED.value}/{nome_manifesto_passo('S2')}"
    metriche = _campo(_json(radice, relativo), relativo, "metriche")
    azzerati = _campo(metriche, relativo, "azzerati")
    quanti = sum(len(azzerati.get(classe, [])) for classe in _CONTROLLATE)
    conta = Misura(
        Fraction(quanti),
        "campioni biologici e controlli positivi senza letture dopo il filtro", relativo)
    statistica = ("perdita media di letture nel filtro, sui campioni biologici e sui controlli "
                  "positivi")
    if "perdita_media_controllata" not in metriche:
        return conta, Misura(None, statistica, relativo,
                             assente="nessun campione biologico o controllo positivo con letture")
    return conta, Misura(_frazione(metriche["perdita_media_controllata"]), statistica, relativo,
                         tolleranza=Fraction(1, 20_000))


def _senza_primer(radice: Path) -> Misura:
    """La quota di letture senza il primer cercato, sull'insieme delle classi
    controllate, dal manifesto del filtro; assente se non c'erano letture.
    Calcolata dai due conteggi e non dalla quota scritta, che e' arrotondata."""
    relativo = f"{Fase.FILTERED.value}/{nome_manifesto_passo('S2')}"
    metriche = _campo(_json(radice, relativo), relativo, "metriche")
    statistica = ("quota di letture senza il primer cercato, sull'insieme dei campioni "
                  "biologici e dei controlli positivi")
    taglio = metriche.get("taglio_primer") or {}
    if not taglio.get("letture_controllate"):
        return Misura(None, statistica, relativo,
                      assente="nessun campione biologico o controllo positivo con letture")
    # Dai due conteggi, non dalla quota arrotondata: la misura e' esatta.
    return Misura(
        Fraction(int(taglio["senza_primer_controllate"]), int(taglio["letture_controllate"])),
        statistica, relativo)


def _varianti(radice: Path) -> Misura:
    """Il numero di varianti distinte, dal riepilogo della tabella delle sequenze."""
    relativo = f"{Fase.SEQTAB.value}/{RIEPILOGO_TABELLA}"
    return Misura(Fraction(int(_campo(_json(radice, relativo), relativo, "varianti"))),
                  "varianti distinte nella tabella delle sequenze", relativo)


def _chimere(radice: Path) -> Misura:
    """La frazione chimerica delle letture delle classi controllate, dal
    riepilogo della rimozione delle chimere; assente se non c'erano letture."""
    relativo = f"{Fase.CHIMERA.value}/{RIEPILOGO_CHIMERE}"
    statistica = ("frazione chimerica delle letture dei campioni biologici e dei controlli "
                  "positivi")
    frazione = _campo(_json(radice, relativo), relativo, "frazione_controllata")
    if frazione is None:
        return Misura(None, statistica, relativo,
                      assente="nessuna lettura di un campione biologico o di un controllo positivo")
    return Misura(_frazione(frazione), statistica, relativo)


def _phylum(radice: Path) -> Misura:
    """Il minimo, fra le classi controllate, della frazione di varianti con il
    phylum assegnato, dal riepilogo della tassonomia.

    Il minimo perche' il controllo giudica ogni classe da sola. Una classe
    senza varianti non entra; assente se non ne entra nessuna.
    """
    relativo = f"{Fase.TAXONOMY.value}/{RIEPILOGO_TASSONOMIA}"
    statistica = ("minimo, fra le classi biologico e controllo positivo, della frazione di "
                  "varianti con il phylum assegnato")
    per_classe = _campo(_json(radice, relativo), relativo, "frazione_con_phylum")
    frazioni = [
        _frazione(per_classe[classe]["varianti"]) for classe in _CONTROLLATE
        if isinstance(per_classe.get(classe), dict)
        and per_classe[classe].get("varianti") is not None
    ]
    if not frazioni:
        return Misura(None, statistica, relativo,
                      assente="nessuna variante di un campione biologico o di un controllo "
                              "positivo")
    return Misura(min(frazioni), statistica, relativo)


def _contaminanti(radice: Path, valori: Mapping[str, Any]) -> Misura:
    """La frazione delle letture dei campioni biologici rimossa come
    contaminante, dal riepilogo della decontaminazione.

    Passi: se il controllo non e' stato eseguito con questi valori la misura e'
    assente; lo e' anche se i controlli negativi con letture sono meno del
    minimo dichiarato, perche' allora la decontaminazione non c'e' stata;
    altrimenti si legge la frazione rimossa dalla modalita' dichiarata (per
    piastra o aggregata), che e' quella che il controllo confronta con la
    soglia, e non dall'altra, riportata solo come diagnostica. Assente anche
    senza letture di campioni biologici nel confronto.
    """
    relativo = f"{Fase.CONTROLS.value}/{RIEPILOGO_DECONTAM}"
    statistica = ("frazione delle letture dei campioni biologici rimossa come contaminante "
                  "dalla modalita' dichiarata")
    motivo = non_eseguito(soglia("qc.max_frac_contaminant"), valori)
    if motivo is not None:
        return Misura(None, statistica, relativo, assente=motivo)
    riepilogo = _json(radice, relativo)
    minimo = _campo(riepilogo, relativo, "min_negativi")
    negativi = _campo(riepilogo, relativo, "confronto", "negativi")
    if minimo is None or negativi < minimo:
        return Misura(None, statistica, relativo,
                      assente=f"nessuna decontaminazione: {negativi} controlli negativi con "
                              "letture, meno di decontam.min_blanks")
    modalita = _campo(riepilogo, relativo, "modalita")
    dichiarata = "per_piastra" if _campo(riepilogo, relativo, "modalita_dichiarata") == "batch" \
        else "aggregata"
    rimossa = _campo(modalita.get(dichiarata) or {}, relativo, "letture_rimosse", "biologico")
    if rimossa is None:
        return Misura(None, statistica, relativo,
                      assente="nessuna lettura di un campione biologico nel confronto")
    return Misura(_frazione(rimossa), statistica, relativo)


def _conformi(radice: Path, valori: Mapping[str, Any]) -> Misura:
    """La frazione dei controlli positivi valutabili conformi al loro livello,
    dal riepilogo dei controlli; assente se il controllo non e' stato eseguito
    o nessun controllo positivo e' valutabile."""
    relativo = f"{Fase.CONTROLS.value}/{RIEPILOGO_CONTROLLI}"
    statistica = "frazione dei controlli positivi valutabili conformi al loro livello"
    motivo = non_eseguito(soglia("ctrl.min_positive_pass_frac"), valori)
    if motivo is not None:
        return Misura(None, statistica, relativo, assente=motivo)
    frazione = _campo(_json(radice, relativo), relativo, "frazione_conformi")
    if frazione is None:
        return Misura(None, statistica, relativo,
                      assente="nessun controllo positivo valutabile")
    return Misura(_frazione(frazione), statistica, relativo)


def _trattenute(radice: Path) -> Misura:
    """Le letture finali sulle letture senza chimere, nell'insieme dei campioni
    finali, dal riepilogo dei filtri finali; zero senza letture di partenza,
    come confronta la fase."""
    relativo = f"{Fase.FINAL.value}/{_INTERMEDI}/{RIEPILOGO_FILTRI}"
    letture = _campo(_json(radice, relativo), relativo, "letture")
    finali = int(_campo(letture, relativo, "finali"))
    senza_chimere = int(_campo(letture, relativo, "nonchimeric_finali"))
    # Senza letture di partenza la fase confronta zero con la soglia.
    frazione = Fraction(finali, senza_chimere) if senza_chimere else Fraction(0)
    return Misura(frazione,
                  "letture finali su letture senza chimere, nell'insieme dei campioni finali",
                  relativo)


def leggi(radice: Path | str, valori: Mapping[str, Any]) -> dict[str, Misura]:
    """La misura di ogni soglia, dalla cartella di un'esecuzione conclusa.

    ``valori`` sono i parametri della configurazione con cui l'esecuzione e'
    stata fatta, nella forma ``gruppo.parametro``. Solleva
    :class:`ArtefattoNonLeggibile` se un artefatto manca o non ha la forma
    attesa.
    """
    try:
        return _leggi(Path(radice), valori)
    except (TypeError, ValueError, AttributeError, KeyError, ArithmeticError) as guasto:
        # Un artefatto che si apre ma non ha i tipi attesi (un campo nullo, un
        # elenco al posto di una mappa) non e' quello di un'esecuzione integra.
        raise ArtefattoNonLeggibile(
            f"un artefatto non ha la forma attesa ({type(guasto).__name__}: {guasto})"
        ) from guasto


def _leggi(radice: Path, valori: Mapping[str, Any]) -> dict[str, Misura]:
    """La misura di ogni soglia, come mappa dalla chiave del parametro.

    Passi: le misure che servono a due soglie (filtro, chimere, fusione) si
    leggono una volta sola; le soglie delle coppie entrano solo in modalita'
    appaiata e quella sulle letture senza primer solo con il taglio cercato,
    perche' negli altri casi i loro artefatti non esistono; le altre si
    leggono sempre. L'avviso e l'arresto delle chimere, e quelli della
    fusione, condividono la stessa misura: cambia solo la soglia con cui si
    confronta.
    """
    azzerati, perdita = _filtro(radice)
    chimere = _chimere(radice)
    delle_coppie: dict[str, Misura] = {}
    if valori.get("io.layout") == "paired-end":
        fusione = _fusione(radice)
        delle_coppie = {
            "pair.max_primer_hit_frac_reverse": _primer(radice, valori, reverse=True),
            "pair.truncLen_shortfall_warn_reverse": _scarto(radice, valori, reverse=True),
            "pair.min_merged_frac_warn": fusione,
            "pair.min_merged_frac_stop": fusione,
            "pair.min_merged_frac_sample_warn": _fusione_per_campione(radice),
        }
    del_taglio: dict[str, Misura] = {}
    if valori.get("filter.primer_removal") == "search":
        del_taglio = {"qc.max_frac_untrimmed": _senza_primer(radice)}
    return {
        **delle_coppie,
        **del_taglio,
        "qc.max_primer_hit_frac": _primer(radice, valori),
        "qc.min_motif_frac": _motivo(radice, valori),
        "qc.max_frac_short_reads": _letture_corte(radice, valori),
        "filter.truncLen_shortfall_warn": _scarto(radice, valori),
        "qc.max_zeroed_samples": azzerati,
        "qc.max_frac_lost_filter": perdita,
        "qc.max_asv_count": _varianti(radice),
        "qc.warn_frac_chimeric": chimere,
        "qc.stop_frac_chimeric": chimere,
        "qc.min_frac_phylum": _phylum(radice),
        "ctrl.min_positive_pass_frac": _conformi(radice, valori),
        "qc.max_frac_contaminant": _contaminanti(radice, valori),
        "qc.min_frac_reads_retained": _trattenute(radice),
    }
