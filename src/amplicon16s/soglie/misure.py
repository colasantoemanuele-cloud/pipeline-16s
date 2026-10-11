"""Le misure di un'esecuzione conclusa che i controlli confrontano con le soglie.

Per ogni soglia di :data:`~amplicon16s.soglie.tabella.SOGLIE` si legge, dagli
artefatti che l'esecuzione ha lasciato nella cartella di uscita, **la stessa
statistica che il controllo confronta con la soglia**: il valore peggiore fra i
file dove il controllo giudica file per file, la mediana dove giudica la
mediana, la frazione d'insieme dove giudica l'insieme. Nulla viene ricalcolato
dalle letture, e la cartella non viene modificata.

Dove l'artefatto riporta la statistica arrotondata e il controllo confronta il
valore non arrotondato, la misura porta una **tolleranza di lettura**: meta'
dell'ultima cifra scritta. Chi propone una soglia peggiora la misura di quella
tolleranza, cosi' che il suggerimento non sia mai piu' stretto del valore vero.

Una misura puo' mancare per due ragioni distinte, che chi legge deve poter
distinguere: il controllo non ha prodotto una misura in questa esecuzione
(``assente``), oppure la misura esiste ma gli artefatti non bastano a
ricostruirla (``non_leggibile``).
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
from amplicon16s.metadata.lettura_inventario import NOME_CROSSWALK
from amplicon16s.metadata.models import CLASSI_CONTROLLATE
from amplicon16s.soglie.tabella import non_eseguito, soglia
from amplicon16s.steps.s00_validate import NOME_INVENTARIO, NOME_LETTURE
from amplicon16s.steps.s01_profile import NOME_LUNGHEZZE
from amplicon16s.steps.s01_profile import NOME_RIEPILOGO as RIEPILOGO_PROFILO
from amplicon16s.steps.s05_seqtab import NOME_RIEPILOGO as RIEPILOGO_TABELLA
from amplicon16s.steps.s06_chimera import NOME_RIEPILOGO as RIEPILOGO_CHIMERE
from amplicon16s.steps.s08_taxonomy import NOME_RIEPILOGO as RIEPILOGO_TASSONOMIA
from amplicon16s.steps.s11_controls import NOME_RIEPILOGO as RIEPILOGO_CONTROLLI
from amplicon16s.steps.s12_decontam import NOME_RIEPILOGO as RIEPILOGO_DECONTAM
from amplicon16s.steps.s13_filtri import NOME_RIEPILOGO as RIEPILOGO_FILTRI

__all__ = ["ArtefattoNonLeggibile", "Misura", "campioni_controllati", "leggi"]

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


def _letture_ispezionate(radice: Path) -> tuple[str, list[dict[str, str]]]:
    """Le righe dei file leggibili dei campioni controllati, come le giudicano
    i gate che leggono le sequenze.
    """
    relativo = f"{Fase.INPUT_VALIDATION.value}/{NOME_LETTURE}"
    classi = _classi(radice)
    righe = [
        r for r in _tsv(radice, relativo)
        if classi.get(r.get("accession", "")) in _CONTROLLATE and not r.get("errore")
    ]
    return relativo, righe


def _primer(radice: Path, valori: Mapping[str, Any]) -> Misura:
    statistica = ("massimo, fra i file dei campioni biologici e dei controlli positivi, della "
                  "frazione delle prime letture che iniziano con il primer")
    relativo, righe = _letture_ispezionate(radice)
    motivo = non_eseguito(soglia("qc.max_primer_hit_frac"), valori)
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
    corte = _campo(_json(radice, manifesto), manifesto, "metriche",
                   "letture_piu_corte_del_troncamento")
    del_profilo = [
        _frazione(_campo(corte, manifesto, classe, "frazione"))
        for classe in _CONTROLLATE if classe in corte
    ]
    origine = f"{manifesto}; {Fase.INPUT_VALIDATION.value}/{NOME_LETTURE}"
    troncamento = int(valori["filter.truncLen"])

    relativo, righe = _letture_ispezionate(radice)
    classi = _classi(radice)
    lunghezze: dict[str, int] | None = None
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
                lunghezze = _corte_per_campione(radice, troncamento)
            corte_del_file = lunghezze.get(riga["accession"], 0)
        else:
            return Misura(
                None, statistica, origine,
                non_leggibile=(
                    f"il file {riga.get('file', riga['accession'])} ha letture piu' corte di "
                    f"filter.truncLen fra le prime {esaminate} e non e' stato letto per intero "
                    "dalla validazione: quante siano non e' scritto in alcun artefatto"))
        voce = per_classe[classi[riga["accession"]]]
        voce[0] += corte_del_file
        voce[1] += esaminate
    della_validazione = [Fraction(c, n) for c, n in per_classe.values() if n]
    tutte = del_profilo + della_validazione
    if not tutte:
        return Misura(None, statistica, origine,
                      assente="nessuna lettura di un campione biologico o di un controllo positivo")
    return Misura(max(tutte), statistica, origine)


def _corte_per_campione(radice: Path, troncamento: int) -> dict[str, int]:
    """Per campione, le letture piu' corte del troncamento secondo il profilo di S1."""
    relativo = f"{Fase.QC_PROFILES.value}/{NOME_LUNGHEZZE}"
    corte: dict[str, int] = defaultdict(int)
    try:
        for riga in _tsv(radice, relativo):
            if int(riga["lunghezza"]) < troncamento:
                corte[riga["campione"]] += int(float(riga["letture"]))
    except (KeyError, ValueError) as guasto:
        raise ArtefattoNonLeggibile(f"{relativo}: {guasto}") from guasto
    return dict(corte)


def _scarto(radice: Path, valori: Mapping[str, Any]) -> Misura:
    relativo = f"{Fase.QC_PROFILES.value}/{RIEPILOGO_PROFILO}"
    minima = int(_campo(_json(radice, relativo), relativo, "lunghezza_minima"))
    return Misura(
        Fraction(max(minima - int(valori["filter.truncLen"]), 0)),
        "basi fra la lettura piu' corta, su tutte le letture, e filter.truncLen", relativo)


def _filtro(radice: Path) -> tuple[Misura, Misura]:
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


def _varianti(radice: Path) -> Misura:
    relativo = f"{Fase.SEQTAB.value}/{RIEPILOGO_TABELLA}"
    return Misura(Fraction(int(_campo(_json(radice, relativo), relativo, "varianti"))),
                  "varianti distinte nella tabella delle sequenze", relativo)


def _chimere(radice: Path) -> Misura:
    relativo = f"{Fase.CHIMERA.value}/{RIEPILOGO_CHIMERE}"
    statistica = ("frazione chimerica delle letture dei campioni biologici e dei controlli "
                  "positivi")
    frazione = _campo(_json(radice, relativo), relativo, "frazione_controllata")
    if frazione is None:
        return Misura(None, statistica, relativo,
                      assente="nessuna lettura di un campione biologico o di un controllo positivo")
    return Misura(_frazione(frazione), statistica, relativo)


def _phylum(radice: Path) -> Misura:
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
    azzerati, perdita = _filtro(radice)
    chimere = _chimere(radice)
    return {
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
