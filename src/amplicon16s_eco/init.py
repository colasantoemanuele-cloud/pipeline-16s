"""Configurazione assistita delle analisi ecologiche: ``eco init``.

Dall'oggetto finale della pipeline scrive una configurazione di partenza delle
analisi, con le tre specie di valore della pipeline (:mod:`amplicon16s.init`):

* **rilevato**: cio' che l'oggetto contiene (ranghi tassonomici, albero e sua
  radice). Non sono parametri: stanno in testa al file, come fatti;
* **proposto**: l'esito di un'euristica, in commento con il valore gia' scritto
  (la variabile biologica, le variabili tecniche, il rango);
* **da scegliere**: le scelte di analisi, in commento, con un valore di
  partenza solo dove c'e' il predefinito di un pacchetto.

Il file non e' eseguibile cosi' com'e': ``validate`` nomina i parametri che
mancano. L'oggetto si legge da una copia e non viene modificato.

Due avvisi guardano i campioni uno per uno: la profondita' di rarefazione
molto piu' bassa della mediana, con i campioni che la abbassano, e le colonne
candidate a variabile biologica che sono confuse fra loro. I dati per campione
si leggono dagli export in testo che la pipeline consegna accanto all'oggetto
(conteggi e tabella dei campioni), dopo aver verificato con il file dei
checksum che siano quelli di questo oggetto; senza, i due avvisi non si
valutano e il riepilogo lo dice.

Con ``--accetta-proposte`` le proposte e i valori di partenza sono righe
attive, con la loro origine nel commento. Restano in commento, con il loro
valore, i parametri condizionali: sono richiesti da un metodo e respinti senza,
e quali metodi usare e' una scelta che ``init`` non fa.
"""

from __future__ import annotations

import csv
import json
import math
import statistics
import tempfile
from dataclasses import dataclass, replace
import textwrap
from pathlib import Path
from typing import Any, Final

from amplicon16s.init.scrittura import SENZA_VALORE, Specie, Voce, valore_yaml
from amplicon16s.io_layer.artifacts import scrivi_atomico
from amplicon16s.io_layer.checksums import checksum_file
from amplicon16s.metadata.crosswalk import _VALORI_ASSENTI
from amplicon16s.metadata.tabelle import COLONNE_INVENTARIO
from amplicon16s.soglie.tabella import LETTURE_FINALI_POCHE

from amplicon16s_eco.config import CONDIZIONALI, FACOLTATIVI, STANDARD_DEL_METODO, ConfigEco
from amplicon16s_eco.esecuzione import _copia, _impronta_oggetto
from amplicon16s_eco.ponte import esegui_r

__all__ = ["QUOTA_DELLA_MEDIANA", "QUOTA_MINIMA_ASSEGNATA", "Esportati", "accetta", "confuse",
           "descrivi", "esportati", "prepara", "rarefazione_bassa", "scrivi", "testo_yaml"]

#: La quota di varianti assegnate sotto la quale un rango non si propone: sotto
#: la meta', la categoria «non assegnato» sarebbe la piu' grande del grafico.
QUOTA_MINIMA_ASSEGNATA: Final = 0.5

#: Le colonne dell'inventario della pipeline che sono variabili tecniche (la
#: piastra di estrazione e la corsa di sequenziamento, le ultime due), e quelle
#: che identificano il campione o ne dicono la classe (le prime tre), che non
#: sono variabili da confrontare. I nomi sono quelli del contratto
#: dell'oggetto, non scritti qui.
_TECNICHE: Final = COLONNE_INVENTARIO[-2:]
_NON_VARIABILI: Final = COLONNE_INVENTARIO[:3]

PROPOSTO, SCELTA = Specie.PROPOSTO, Specie.DA_SCEGLIERE

#: Sotto questa quota della mediana delle letture per campione la profondita'
#: di rarefazione (la minima fra i campioni) e' molto piu' bassa del resto: e'
#: la stessa quota con cui la pipeline riconosce un campione con molte meno
#: letture degli altri.
QUOTA_DELLA_MEDIANA: Final = 0.1
#: Le letture che bastano a descrivere un campione: sotto, toglierlo e' il
#: rimedio; sopra, una profondita' bassa rispetto agli altri si puo' tenere.
#: E' lo stesso numero con cui la pipeline propone le letture minime.
LETTURE_SUFFICIENTI: Final = LETTURE_FINALI_POCHE

#: I file che la pipeline consegna accanto all'oggetto finale: l'elenco dei
#: checksum, i conteggi (varianti sulle righe, un campione per colonna) e la
#: tabella dei campioni. Sono i nomi del contratto di consegna.
_CHECKSUM: Final = "checksum.sha256"
_CONTEGGI: Final = "conteggi.tsv"
_CAMPIONI: Final = "metadati.tsv"
#: Quanti campioni o gruppi si nominano in un avviso.
_NOMINATI: Final = 8


@dataclass(frozen=True)
class Esportati:
    """Che cosa gli export in testo dicono dei campioni dell'oggetto."""

    #: Le letture di ogni campione.
    letture: dict[str, int]
    #: La tabella dei campioni: una riga per campione, colonna -> valore.
    campioni: list[dict[str, str]]


def esportati(oggetto: Path | str) -> Esportati | None:
    """Gli export in testo consegnati accanto all'oggetto, se sono i suoi.

    Razionale sistemistico: i conteggi in testo valgono per l'oggetto solo se
    il file dei checksum della stessa cartella elenca l'oggetto con la sua
    impronta di ora e i due export con le loro. Un oggetto spostato, o
    sostituito dopo la consegna, non ha export su cui contare: ``None``.
    """
    oggetto = Path(oggetto)
    cartella = oggetto.parent
    try:
        attesi = {}
        for riga in (cartella / _CHECKSUM).read_text(encoding="utf-8").splitlines():
            impronta, _, nome = riga.strip().partition("  ")
            attesi[nome] = impronta
        for nome in (oggetto.name, _CONTEGGI, _CAMPIONI):
            trovata = checksum_file(cartella / nome).removeprefix("sha256:")
            if attesi.get(nome, "").removeprefix("sha256:") != trovata:
                return None
        with open(cartella / _CONTEGGI, encoding="utf-8", newline="") as file:
            righe = csv.reader(file, delimiter="\t")
            nomi = next(righe)[1:]
            letture = dict.fromkeys(nomi, 0)
            for riga in righe:
                for nome, valore in zip(nomi, riga[1:]):
                    letture[nome] += int(valore)
        with open(cartella / _CAMPIONI, encoding="utf-8", newline="") as file:
            campioni = list(csv.DictReader(file, delimiter="\t"))
    except (OSError, ValueError, StopIteration, UnicodeDecodeError):
        return None
    if not letture or len(campioni) != len(letture):
        return None
    return Esportati(letture, campioni)


def rarefazione_bassa(letture: dict[str, int], nomi: dict[str, str] | None = None) -> str | None:
    """L'avviso sulla profondita' di rarefazione, o ``None`` se non serve.

    Senza ``alpha.rarefy_depth`` le analisi rarefanno alla profondita' del
    campione piu' povero: se e' sotto un decimo della mediana, quasi tutta
    l'informazione degli altri campioni va persa. L'avviso nomina i campioni
    sotto quella quota, dal piu' povero, e dice come rimediare.
    """
    if len(letture) < 2:
        return None
    mediana = statistics.median(letture.values())
    minima = min(letture.values())
    if minima >= QUOTA_DELLA_MEDIANA * mediana:
        return None
    limite = math.ceil(QUOTA_DELLA_MEDIANA * mediana)
    poveri = sorted((n, c) for c, n in letture.items() if n < limite)
    nomi = nomi or {}
    elenco = ", ".join(f"{nomi.get(c) or c} ({n})" for n, c in poveri[:_NOMINATI]) + (
        f" e altri {len(poveri) - _NOMINATI}" if len(poveri) > _NOMINATI else "")
    testa = (
        f"Rarefazione: la profondita' sarebbe {minima} letture, quella del campione piu' "
        f"povero, meno di un decimo della mediana ({mediana:.0f}): l'alfa diversita' di ogni "
        f"campione sarebbe calcolata su {minima} letture. La abbassano: {elenco}. ")
    if minima >= LETTURE_SUFFICIENTI:
        # Abbastanza letture per descrivere un campione: toglierlo costa piu'
        # di quanto renda, soprattutto se e' di un gruppo con pochi campioni.
        return testa + (
            f"Con almeno {LETTURE_SUFFICIENTI} letture per campione di solito si lascia "
            "cosi'. Per una profondita' piu' alta si tengono fuori quei campioni con "
            "design.subset, dopo aver controllato che non siano un intero gruppo dello studio.")
    return testa + (
        f"Rimedio: nella configurazione della pipeline alzare qc.min_reads_final a "
        f"{LETTURE_SUFFICIENTI} (escono i campioni con meno letture) e rilanciare resume, "
        "poi eco init sul nuovo oggetto; oppure tenerli fuori dalle analisi con "
        "design.subset. Prima va controllato che non siano un intero gruppo dello studio.")


def _presente(valore: str | None) -> bool:
    """Vero se il valore non e' vuoto ne' uno dei modi di scrivere un dato assente."""
    return (valore or "").strip().casefold() not in _VALORI_ASSENTI


def confuse(campioni: list[dict[str, str]], colonne: list[str]) -> list[str]:
    """Gli avvisi sulle colonne candidate che sono confuse fra loro.

    Razionale sperimentale: se ogni gruppo di una colonna cade in un solo
    gruppo dell'altra, un confronto fra i gruppi della seconda non si
    distingue da uno fra i gruppi della prima. Lo stesso vale, in parte,
    quando un gruppo intero di una colonna cade in un solo gruppo dell'altra
    (tutti i campioni di una condizione sono di un solo tessuto): la
    differenza che i test trovano puo' venire dall'una o dall'altra. Due
    colonne che dividono i campioni negli stessi gruppi sono la stessa
    variabile, e si confrontano con le altre una volta sola. Contano i gruppi
    di almeno due campioni, fra i campioni con un valore in entrambe.
    """
    def partizione(colonna: str) -> tuple[int, ...]:
        """Per ogni campione, l'indice del gruppo in cui la colonna lo mette; -1
        se il campione non ha un valore."""
        indice: dict[str, int] = {}
        return tuple(indice.setdefault(c.get(colonna, ""), len(indice)) if _presente(
            c.get(colonna)) else -1 for c in campioni)

    avvisi = []
    distinte: list[str] = []
    for colonna in colonne:
        uguale = next((d for d in distinte if partizione(d) == partizione(colonna)), None)
        if uguale is None:
            distinte.append(colonna)
        else:
            avvisi.append(
                f"{uguale} e {colonna} dividono i campioni negli stessi gruppi: per le "
                "analisi sono la stessa variabile, e basta sceglierne una.")

    def dentro(a: str, b: str) -> tuple[dict[str, str], int]:
        """I gruppi di ``a`` che cadono in un solo gruppo di ``b``, e quanti
        gruppi di almeno due campioni ha ``a``."""
        gruppi: dict[str, set[str]] = {}
        quanti: dict[str, int] = {}
        for c in campioni:
            if _presente(c.get(a)) and _presente(c.get(b)):
                gruppi.setdefault(c[a], set()).add(c[b])
                quanti[c[a]] = quanti.get(c[a], 0) + 1
        ripetuti = [g for g in gruppi if quanti[g] >= 2]
        return {g: next(iter(gruppi[g])) for g in ripetuti if len(gruppi[g]) == 1}, len(ripetuti)

    for i, a in enumerate(distinte):
        for b in distinte[i + 1:]:
            a_in_b, gruppi_a = dentro(a, b)
            b_in_a, gruppi_b = dentro(b, a)
            if min(gruppi_a, gruppi_b) < 2 or not (a_in_b or b_in_a):
                continue
            rimedio = ("Un confronto fra i gruppi dell'una misura anche la differenza fra i "
                       "gruppi dell'altra: scegli come design.variable quella che lo studio "
                       "confronta sapendo che l'altra vi e' confusa, oppure restringi con "
                       "design.subset ai campioni di un solo gruppo dell'altra.")
            if len(a_in_b) == gruppi_a or len(b_in_a) == gruppi_b:
                interna, esterna = (a, b) if len(a_in_b) == gruppi_a else (b, a)
                avvisi.append(
                    f"{a} e {b} sono confuse: ogni gruppo di {interna} cade in un solo gruppo "
                    f"di {esterna}. " + rimedio)
                continue
            parti = [f"i campioni {g!r} di {x} sono tutti {v!r} di {y}"
                     for x, y, trovati in ((a, b, a_in_b), (b, a, b_in_a))
                     for g, v in sorted(trovati.items())]
            avvisi.append(
                f"{a} e {b} sono confuse in parte: " + "; ".join(parti[:_NOMINATI])
                + (f"; e altri {len(parti) - _NOMINATI} gruppi" if len(parti) > _NOMINATI else "")
                + ". " + rimedio)
    return avvisi


def descrivi(oggetto: Path | str) -> dict[str, Any]:
    """Che cosa contiene l'oggetto: campioni, ranghi, albero, colonne.

    Solleva :class:`~amplicon16s_eco.catalogo.ErroreEco` se l'oggetto non e'
    utilizzabile (``E-ECO-04``) o se R non si esegue (``E-ECO-90``).
    """
    oggetto = Path(oggetto)
    impronta = _impronta_oggetto(oggetto)
    with tempfile.TemporaryDirectory(prefix="amplicon16s-eco-") as temporanea:
        lavoro = Path(temporanea)
        copia = _copia(oggetto, lavoro, impronta)
        esegui_r("descrivi.R", {"oggetto": str(copia),
                                "valori_mancanti": sorted(_VALORI_ASSENTI)}, lavoro, "descrizione")
        return json.loads((lavoro / "descrizione.json").read_text(encoding="utf-8"))


def prepara(
    descrizione: dict[str, Any], dei_campioni: Esportati | None = None
) -> tuple[dict[str, Voce], list[str]]:
    """Le voci della configurazione e i fatti rilevati sull'oggetto.

    ``dei_campioni`` sono gli export in testo dell'oggetto, se ci sono: con
    quelli si valutano la profondita' di rarefazione e le colonne confuse.
    """
    voci: dict[str, Voce] = {}

    def metti(chiave: str, specie: Specie, valore: Any, origine: str, *note: str) -> None:
        """Registra la voce di un parametro: specie, valore, origine e note non vuote."""
        voci[chiave] = Voce(chiave, specie, valore, origine, tuple(n for n in note if n))

    colonne = descrizione["colonne"]
    ranghi = descrizione["ranghi"]
    albero = descrizione["albero"]
    fatti = [
        f"Oggetto: {descrizione['campioni']} campioni, {descrizione['varianti']} varianti.",
        # Tutte le colonne, non solo le candidate: una variabile continua o il
        # soggetto di un disegno a misure ripetute si scelgono fra queste.
        "Colonne della tabella dei campioni (valori distinti; «numerica» se tutti i valori "
        "sono numeri): " + ", ".join(
            f"{c['nome']} ({c['distinti']}{', numerica' if c['numerica'] else ''})"
            for c in colonne) + ".",
        "Ranghi tassonomici (quota di varianti assegnate): "
        + ", ".join(f"{r['rango']} {r['quota_assegnata']:.2f}" for r in ranghi) + ".",
        "Albero filogenetico: " + (
            ("presente, " + ("radicato" if albero["radicato"] else "senza radice")
             + ("" if albero["con_lunghezze"] else ", senza lunghezze dei rami"))
            if albero["presente"] else "assente (le distanze UniFrac non sono disponibili)") + ".",
    ]

    # Variabile biologica: proposta solo se una sola colonna ha gruppi ripetuti.
    # Il numero di gruppi non dice quale colonna lo studio confronta: una
    # coordinata di laboratorio si ripete quanto un fattore sperimentale.
    # Una colonna che sembra una misura non e' candidata a variabile
    # categorica: si propone come variabile continua.
    continue_ = sorted(
        c["nome"] for c in colonne
        if c["nome"] not in _TECNICHE + _NON_VARIABILI and _sembra_continua(c))
    candidate = sorted(
        (c for c in colonne
         if c["nome"] not in _TECNICHE + _NON_VARIABILI and c["gruppi_ripetuti"] >= 2
         and c["nome"] not in continue_),
        key=lambda c: c["nome"])
    if continue_:
        metti("design.continuous_variables", PROPOSTO, continue_,
              "le colonne con valori tutti numerici, almeno 10 distinti e distinti in "
              "almeno l'80% dei campioni con un valore: misure, non gruppi. Per ciascuna "
              "l'analisi calcola la correlazione con l'alfa diversita' e una PERMANOVA con "
              "un termine numerico")
    if len(candidate) == 1:
        sola = candidate[0]
        metti("design.variable", PROPOSTO, sola["nome"],
              f"la sola colonna con almeno due gruppi di almeno due campioni "
              f"({sola['gruppi_ripetuti']} su {sola['distinti']} valori)",
              *(["attenzione: i valori sono tutti numerici e quasi tutti diversi; le analisi "
                 "la tratterebbero come categorica"] if _sembra_continua(sola) else []))
    elif candidate:
        metti("design.variable", SCELTA, SENZA_VALORE,
              "piu' colonne hanno almeno due gruppi di almeno due campioni, e init non sa "
              "quale lo studio confronta: il numero di gruppi non lo dice",
              "colonne con gruppi ripetuti (numero di gruppi): "
              + ", ".join(f"{c['nome']} ({c['gruppi_ripetuti']})" for c in candidate),
              *(["attenzione: " + ", ".join(c["nome"] for c in candidate if _sembra_continua(c))
                 + " ha valori tutti numerici e quasi tutti diversi; le analisi la "
                 "tratterebbero come categorica"]
                if any(_sembra_continua(c) for c in candidate) else []))
    else:
        metti("design.variable", SCELTA, SENZA_VALORE,
              "nessuna colonna ha almeno due gruppi di almeno due campioni: i test fra "
              "gruppi non si potranno eseguire, restano le uscite descrittive",
              "colonne dell'oggetto: " + ", ".join(c["nome"] for c in colonne))
    if dei_campioni is None:
        fatti.append(
            "Profondita' di rarefazione e variabili confuse: non valutate, perche' accanto "
            "all'oggetto non ci sono gli export in testo della pipeline con il loro file dei "
            f"checksum ({_CONTEGGI}, {_CAMPIONI}, {_CHECKSUM}), o non sono quelli di questo "
            "oggetto.")
    else:
        nomi = {c.get(COLONNE_INVENTARIO[0], ""): c.get(COLONNE_INVENTARIO[1], "")
                for c in dei_campioni.campioni}
        bassa = rarefazione_bassa(dei_campioni.letture, nomi)
        if bassa is not None:
            fatti.append(bassa)
        # I gruppi di ogni candidata con quanti campioni hanno: servono a
        # scegliere la variabile e la dimensione minima dei gruppi.
        dei_gruppi = []
        for c in candidate:
            quanti: dict[str, int] = {}
            for riga in dei_campioni.campioni:
                if _presente(riga.get(c["nome"])):
                    quanti[riga[c["nome"]]] = quanti.get(riga[c["nome"]], 0) + 1
            if 0 < len(quanti) <= _NOMINATI:
                dei_gruppi.append(c["nome"] + ": " + ", ".join(
                    f"{valore!r} {n}" for valore, n in sorted(quanti.items())))
        if dei_gruppi:
            fatti.append("Gruppi delle colonne candidate a variabile da confrontare "
                         "(campioni per gruppo): " + "; ".join(dei_gruppi) + ".")
        delle_confuse = confuse(dei_campioni.campioni, [c["nome"] for c in candidate])
        for avviso in delle_confuse:
            fatti.append("Variabili confuse: " + avviso)
        if delle_confuse and "design.variable" in voci:
            voce = voci["design.variable"]
            voci["design.variable"] = replace(voce, note=(
                *voce.note, *("attenzione: " + a for a in delle_confuse)))
    tecniche = sorted(
        (c for c in colonne if c["nome"] in _TECNICHE and c["distinti"] > 1),
        key=lambda c: c["distinti"])
    metti("design.technical_variables", PROPOSTO, [c["nome"] for c in tecniche],
          ("le colonne di corsa e di piastra dell'inventario che hanno piu' di un valore, "
           "dalla piu' grossolana alla piu' annidata (meno valori prima)") if tecniche else
          "ne' la corsa ne' la piastra hanno piu' di un valore fra i campioni")

    # Dal rango piu' generale in giu', finche' la quota regge: un rango fine
    # non si propone attraverso uno piu' generale che non la raggiunge.
    proponibili = []
    for rango in ranghi:
        if rango["quota_assegnata"] < QUOTA_MINIMA_ASSEGNATA:
            break
        proponibili.append(rango)
    if proponibili:
        rango = proponibili[-1]
        metti("comp.rank", PROPOSTO, rango["rango"],
              f"il rango piu' fine a cui almeno il {QUOTA_MINIMA_ASSEGNATA:.0%} delle "
              f"varianti ha un'assegnazione ({rango['quota_assegnata']:.2f} delle varianti; "
              "la quota e' sul numero delle varianti, non sulle letture)")
    else:
        metti("comp.rank", SCELTA, SENZA_VALORE,
              f"nessun rango ha almeno il {QUOTA_MINIMA_ASSEGNATA:.0%} delle varianti assegnate")
    if albero["presente"] and albero["radicato"] and albero["con_lunghezze"]:
        metti("beta.unifrac_root", PROPOSTO, "existing",
              "l'albero dell'oggetto e' gia' radicato; serve solo se beta.distances "
              "comprende una UniFrac, ed e' respinto altrimenti")

    metti("run.seed", SCELTA, SENZA_VALORE,
          "un intero qualunque, da dichiarare: fissa la rarefazione e le permutazioni")
    metti("comp.top_n", SCELTA, SENZA_VALORE,
          "quanti taxa mostrare nei grafici prima di riunire gli altri")
    metti("beta.distances", SCELTA, SENZA_VALORE,
          "fra bray, jaccard, aitchison" + (
              ", unifrac_weighted, unifrac_unweighted" if albero["presente"] else
              " (le due UniFrac richiedono un albero, che l'oggetto non ha)"))
    metti("beta.clr_pseudocount", SCELTA, SENZA_VALORE,
          "solo con aitchison: lo pseudoconteggio sommato ai conteggi prima del CLR")
    metti("ord.methods", SCELTA, SENZA_VALORE, "fra pcoa e nmds; [] per nessuna ordinazione")
    metti("ord.distances", SCELTA, SENZA_VALORE, "fra quelle di beta.distances")
    metti("ord.pcoa_correction", SCELTA, "none",
          "solo con pcoa. Partenza: il predefinito di ape::pcoa (correction = \"none\"); "
          "cailliez o lingoes correggono gli autovalori negativi delle distanze non euclidee")
    metti("ord.nmds_trymax", SCELTA, 20,
          "solo con nmds. Partenza: il predefinito di vegan::metaMDS (trymax = 20)")
    metti("stat.min_group_size", SCELTA, SENZA_VALORE,
          "quanti campioni deve avere un gruppo per entrare nei test (almeno 2)")
    metti("stat.significance_level", SCELTA, SENZA_VALORE,
          "il livello sotto il quale le dispersioni dei gruppi si dichiarano diverse")
    return voci, fatti


def accetta(voci: dict[str, Voce]) -> dict[str, Voce]:
    """Le voci con proposte e valori di partenza accettati in blocco.

    Un condizionale non si attiva: e' respinto se il metodo che lo richiede
    non e' fra quelli scelti, e gli elenchi dei metodi non hanno un valore di
    partenza. Resta in commento con il suo valore, e lo dice.
    """
    accettate = {}
    for chiave, voce in voci.items():
        if voce.valore is SENZA_VALORE:
            pass
        elif chiave in CONDIZIONALI:
            elenco, metodi = CONDIZIONALI[chiave]
            voce = replace(voce, note=(
                *voce.note,
                f"non attivato da --accetta-proposte: va dichiarato solo se {elenco} "
                f"comprende {' o '.join(metodi)}, ed e' respinto altrimenti"))
        elif voce.specie is PROPOSTO:
            voce = replace(voce, accettata=True,
                           origine=voce.origine + "; accettato con --accetta-proposte")
        else:
            voce = replace(voce, specie=Specie.PARTENZA,
                           origine=voce.origine + "; accettato con --accetta-proposte")
        accettate[chiave] = voce
    return accettate


def _sembra_continua(colonna: dict[str, Any]) -> bool:
    """La stessa regola dell'avviso sulle variabili che sembrano continue."""
    return bool(colonna["numerica"]) and colonna["distinti"] >= 10 and (
        colonna["distinti"] >= 0.8 * colonna["con_valore"])


def testo_yaml(
    voci: dict[str, Voce], fatti: list[str], oggetto: str, accettate: bool = False
) -> str:
    """Il testo del file di configurazione delle analisi."""
    da_scegliere = [v for v in voci.values() if v.specie is SCELTA]
    proposti = [v for v in voci.values() if v.specie is PROPOSTO]
    attive = [v for v in voci.values() if v.attiva]

    def commento(testo: str, rientro: str = "  ##") -> list[str]:
        """Il testo, a capo entro la larghezza del file, in righe con il prefisso dato."""
        # Dentro i gruppi le spiegazioni portano ##, le righe dei parametri in
        # commento un solo #: non si scambiano fra loro.
        return [f"{rientro} {r}" for r in textwrap.wrap(
            testo, 77 - len(rientro), break_long_words=False, break_on_hyphens=False)]

    righe = [
        "# " + "=" * 76,
        "# Configurazione di partenza delle analisi ecologiche, scritta da eco init",
        "# " + "=" * 76,
        "#",
        f"# Oggetto esaminato: {oggetto}",
        "# [RILEVATO] sull'oggetto:",
        *[riga for fatto in fatti for riga in commento(fatto, "#  ")],
        "#",
        *([
            "#   [PROPOSTO]      l'esito di un'euristica: riga in commento con il valore",
            "#                   gia' scritto. Per accettarlo, togli il commento.",
            "#   [DA SCEGLIERE]  una scelta di analisi: riga in commento, con un valore di",
            "#                   partenza solo dove c'e' il predefinito di un pacchetto.",
        ] if not accettate else [
            "#   [PROPOSTO]      l'esito di un'euristica: riga attiva, accettata con",
            "#                   --accetta-proposte. Controllala comunque.",
            "#   [PARTENZA]      il valore di partenza di una scelta, accettato con",
            "#                   --accetta-proposte: riga attiva.",
            "#   [DA SCEGLIERE]  una scelta di analisi: riga in commento, da scrivere. Un",
            "#                   condizionale resta in commento anche se ha un valore.",
        ]),
        "#   [CONDIZIONALE]  va dichiarato solo se il metodo indicato e' richiesto.",
        "#",
        f"# Proposti: {len(proposti)}; da scegliere: {len(da_scegliere)}."
        + (f" Righe attive, accettate con --accetta-proposte: {len(attive)}." if accettate
           else ""),
        "# Da scegliere:",
        *commento(", ".join(v.chiave for v in da_scegliere), "#"),
        "# Le spiegazioni cominciano con ##, le righe dei parametri in commento con #.",
        "",
    ]
    for gruppo, campo_gruppo in ConfigEco.model_fields.items():
        modello = campo_gruppo.annotation
        righe.append(f"# {(modello.__doc__ or '').strip().splitlines()[0]}")
        righe.append(f"{gruppo}:")
        for nome, campo in modello.model_fields.items():
            chiave = f"{gruppo}.{nome}"
            righe += commento(campo.description or "")
            voce = voci.get(chiave)
            if chiave in CONDIZIONALI:
                elenco, metodi = CONDIZIONALI[chiave]
                righe += commento(f"[CONDIZIONALE] solo se {elenco} comprende "
                                  + " o ".join(metodi))
            if voce is None and chiave in CONDIZIONALI:
                # Senza una proposta il condizionale non ha un valore da
                # mostrare: il nullo dello schema e' la sua assenza, e scritto
                # sembrerebbe una scelta.
                righe.append(f"  # {nome}:")
                continue
            if voce is None:
                predefinito = campo.get_default(call_default_factory=True)
                if chiave in STANDARD_DEL_METODO:
                    righe += commento(f"[STANDARD] {STANDARD_DEL_METODO[chiave][1]}")
                elif chiave in FACOLTATIVI:
                    righe += commento("[FACOLTATIVO] il predefinito e' l'assenza")
                valore = {} if isinstance(predefinito, dict) else (
                    list(predefinito) if isinstance(predefinito, tuple) else predefinito)
                scritto = "{}" if valore == {} else valore_yaml(valore)
                righe.append(f"  # {nome}: {scritto}")
                continue
            righe += commento(f"[{voce.specie.value}] {voce.origine}")
            for nota in voce.note:
                righe += commento(nota)
            scritto = "" if voce.valore is SENZA_VALORE else f" {valore_yaml(voce.valore)}"
            righe.append(f"  {'' if voce.attiva else '# '}{nome}:{scritto}")
        righe.append("")
    return "\n".join(righe)


def scrivi(
    oggetto: Path | str, uscita: Path | str, accettate: bool = False
) -> tuple[dict[str, Voce], list[str]]:
    """Esamina l'oggetto e scrive la configurazione di partenza.

    Con ``accettate`` le proposte e i valori di partenza sono righe attive
    (:func:`accetta`). Solleva ``FileExistsError`` se il file esiste gia': non
    lo sovrascrive.
    """
    uscita = Path(uscita)
    if uscita.exists():
        raise FileExistsError(uscita)
    voci, fatti = prepara(descrivi(oggetto), esportati(oggetto))
    if accettate:
        voci = accetta(voci)
    scrivi_atomico(uscita, testo_yaml(voci, fatti, Path(oggetto).name, accettate) + "\n")
    return voci, fatti
