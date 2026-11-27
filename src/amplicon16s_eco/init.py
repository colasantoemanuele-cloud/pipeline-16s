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

Sopra il gruppo ``design`` il file elenca le **divisioni suggerite dai
metadati**: il catalogo delle divisioni (:mod:`amplicon16s_eco.divisioni`)
applicato alla tabella dei campioni, piu' le divisioni generiche (una colonna
con gruppi ripetuti, il prefisso dei valori di una colonna), senza doppioni e
solo quelle che funzionano sul dataset. Per ciascuna una riga, con il nome e
una spiegazione; poi come sceglierne una (il nome nei campi che accettano una
colonna) e una riga con i ranghi tassonomici, che non sono divisioni. I gruppi
di ogni divisione, con i campioni, non stanno nel file di configurazione ma in
un file a parte accanto a esso (:func:`nome_del_tsv`), una riga per gruppo.

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

import yaml

from amplicon16s.init.scrittura import SENZA_VALORE, Specie, Voce, valore_yaml
from amplicon16s.io_layer.artifacts import scrivi_atomico
from amplicon16s.io_layer.checksums import checksum_file
from amplicon16s.metadata.crosswalk import _VALORI_ASSENTI
from amplicon16s.metadata.tabelle import COLONNE_INVENTARIO
from amplicon16s.soglie.tabella import LETTURE_FINALI_POCHE

from amplicon16s_eco.config import (
    CONDIZIONALI,
    FACOLTATIVI,
    STANDARD_DEI_CONDIZIONALI,
    STANDARD_DEL_METODO,
    ConfigEco,
    condizione,
)
from amplicon16s_eco.divisioni import (DIMENSIONE_DELLE_PROPOSTE, PREFISSO, Proposte, Tabella,
                                       proponi)
from amplicon16s_eco.esecuzione import _copia, _impronta_oggetto, descrizione_della_copia

__all__ = ["PREFISSO", "QUOTA_DELLA_MEDIANA", "QUOTA_MINIMA_ASSEGNATA", "Esportati", "accetta",
           "blocco_delle_divisioni", "confuse", "derivazioni_della_pipeline", "descrivi",
           "divisioni_confuse", "esportati", "nome_del_tsv", "prepara", "proposte_di",
           "rarefazione_bassa", "scrivi", "testo_yaml", "tsv_delle_divisioni"]

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

#: Le dimensioni minime dei gruppi per cui init conta i gruppi di ogni colonna:
#: stat.min_group_size non e' ancora scelto, e questi sono i valori di solito
#: considerati.
SOGLIE_DEI_GRUPPI: Final = (2, 5, 10)

#: La colonna che la pipeline deriva da una regola (il modulo, con
#: meta.derive_module e meta.module_regex): il nome e' quello del contratto
#: dell'oggetto.
_DERIVATE: Final = (COLONNE_INVENTARIO[5],)

_RIMANDO: Final = ("si sceglie per nome fra le divisioni suggerite elencate in testa al file, "
                   "sopra design; oppure una colonna dell'oggetto o un gruppo derivato di "
                   "design.groupings")

#: Le proposte che --accetta-proposte non attiva mai: una regola sui dati da
#: controllare prima di usarla.
_MAI_ACCETTATE: Final = ("design.groupings",)

#: Note fisse dei parametri che init non propone (facoltativi e standard): le
#: alternative e le fonti per sceglierli, scritte sopra la riga in commento.
NOTE_DEI_PARAMETRI: Final[dict[str, tuple[str, ...]]] = {
    "design.other_variables": (_RIMANDO,),
    "ord.color_by": (
        "per le figure descrittive, per esempio un raggruppamento visto dopo aver "
        "guardato i dati: nessun test, e il titolo lo dice; " + _RIMANDO,),
    "run.analysis_rank": (
        "nullo (il predefinito): ogni calcolo sulle varianti, la pratica comune. Con un "
        "rango, alfa diversita', distanze, ordinazioni, test e metodo vincolato lavorano "
        "sui taxa agglomerati con phyloseq::tax_glom(NArm = FALSE): si perde la "
        "risoluzione delle varianti, si guadagnano meno taxa rari e meno zeri; UniFrac "
        "usa l'albero potato. E' una scelta distinta dal rango dei grafici (comp.rank), "
        "che non puo' essere piu' fine (RAZIONALE.md, sezione 4.3)",),
    "comp.facet_wrap": (
        "barre per campione, un pannello per gruppo; " + _RIMANDO,),
    "comp.facet_grid": (
        "barre medie per gruppo di design.variable in una griglia {rows: <nome>, cols: "
        "<nome>}, oppure il nome di un incrocio fra le divisioni suggerite (\"A per B\"); "
        + _RIMANDO,),
    "design.groupings": (
        "per una divisione personalizzata, che il catalogo non conosce: {<nome>: {from: "
        "<colonna>, groups: {<gruppo>: [<valori>]}}} oppure {<nome>: {from: <colonna>, "
        "pattern: <espressione con un gruppo di cattura>}}; le divisioni suggerite non "
        "servono qui, si nominano direttamente",),
    "constrained.method": (
        "ordinazione vincolata: rda sui conteggi trasformati, dbrda sulle distanze; "
        "le variabili si scelgono fra quelle del disegno e i gruppi derivati, le "
        "condizioni (Condition()) di solito fra le variabili tecniche",),
    "comp.group_by": (
        "le barre della composizione per gruppo; vuoto: le variabili del disegno; "
        + _RIMANDO,),
    "plot.palette": (
        "un nome di grDevices::hcl.pals() (\"Dark 3\" e' quella di sempre), \"Okabe-Ito\" "
        "(adatta ai daltonici), un nome di RColorBrewer (\"Set2\", \"Dark2\"...) o di "
        "viridisLite (\"viridis\", \"cividis\"...), oppure un elenco di colori "
        "esadecimali; con piu' gruppi che colori la palette si estende per "
        "interpolazione",),
    "filter.min_relative_abundance": (
        "le varianti rare: un'alternativa o un complemento alla rarefazione. Nullo (il "
        "predefinito): nessun filtro. La soglia fissa e' una scelta dello studio: "
        "Bokulich et al. 2013 (Nat Methods 10:57) propongono lo 0,005% delle letture "
        "totali (filter.scope total), pensato per le OTU; Cao et al. 2021 (Front "
        "Microbiol 11:607325) e Nikodemova et al. 2023 (Front Cell Infect Microbiol "
        "13:1165295) misurano l'effetto del filtro: forte su ricchezza e Chao1, piccolo "
        "su Shannon e Bray-Curtis. Si applica dopo design.subset e prima di ogni "
        "analisi (RAZIONALE.md, sezione 4.2)",
    ),
}

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


def _gruppi_della_colonna(colonna: dict[str, Any]) -> dict[str, int]:
    """I valori di una colonna con i loro campioni, dalla descrizione."""
    return {v["valore"]: int(v["campioni"]) for v in colonna.get("valori") or []}


def _non_selezionabile(colonna: dict[str, Any]) -> str | None:
    """Perche' una colonna non e' selezionabile come gruppo, o ``None``."""
    if colonna["nome"] in _NON_VARIABILI:
        return "identifica il campione o ne dice la classe"
    if colonna["nome"] in _TECNICHE:
        return "variabile tecnica: si dichiara in design.technical_variables"
    if colonna["con_valore"] == 0:
        return "nessun valore"
    if colonna["distinti"] == 1:
        return "un solo valore"
    if _sembra_continua(colonna):
        return "continua (valori numerici quasi tutti diversi): design.continuous_variables"
    if colonna["distinti"] >= colonna["con_valore"]:
        return "tutti i valori diversi: un gruppo per campione"
    return None


def derivazioni_della_pipeline(oggetto: Path | str) -> dict[str, str] | None:
    """Le colonne che la pipeline ha derivato da una regola, con la regola.

    Razionale sistemistico: la regola sta nei file dell'esecuzione accanto
    all'oggetto: ``10_phyloseq/colonne_metadati.tsv`` (verificato sul
    manifesto di S10) descrive come e' nata ogni colonna, e
    ``00_config/resolved.yaml`` porta i valori dei parametri che l'hanno
    prodotta. Senza quei file (un oggetto spostato) la regola non si conosce:
    ``None``.
    """
    esecuzione = Path(oggetto).resolve().parent.parent
    tabella = esecuzione / "10_phyloseq" / "colonne_metadati.tsv"
    try:
        manifesto = json.loads((tabella.parent / "manifest.json").read_text(encoding="utf-8"))
        attesa = next(a["checksum"] for a in manifesto["artefatti"]
                      if a["nome"] == tabella.name)
        if checksum_file(tabella).removeprefix("sha256:") != attesa.removeprefix("sha256:"):
            return None
        with open(tabella, encoding="utf-8", newline="") as file:
            righe = {r["colonna"]: r for r in csv.DictReader(file, delimiter="\t")}
    except (OSError, ValueError, KeyError, StopIteration, TypeError):
        return None
    try:
        documenti = list(yaml.safe_load_all(
            (esecuzione / "00_config" / "resolved.yaml").read_text(encoding="utf-8")))
        # Con una ripresa che cambia i valori il file ne registra un'altra: vale
        # l'ultima.
        meta = (documenti[-1].get("parametri") or {}).get("meta") or {}
    except (OSError, yaml.YAMLError, AttributeError, IndexError):
        meta = {}
    derivazioni = {}
    for nome in _DERIVATE:
        if nome not in righe:
            continue
        regola = righe[nome].get("valore", "")
        if meta:
            regola += (f"; in questa esecuzione meta.derive_module = {meta.get('derive_module')}, "
                       f"meta.module_regex = {meta.get('module_regex')!r}, meta.module_column = "
                       f"{meta.get('module_column')!r}")
        derivazioni[nome] = regola
    return derivazioni


def proposte_di(descrizione: dict[str, Any]) -> Proposte:
    """Le divisioni che funzionano sui metadati dell'oggetto: il catalogo e le
    divisioni generiche delle colonne selezionabili (:func:`_non_selezionabile`).

    Fra le generiche le colonne dello studio vengono prima di quelle
    dell'inventario della pipeline, ciascun gruppo nell'ordine dell'oggetto:
    fra una colonna dello studio e una dell'inventario con gli stessi gruppi
    (la colonna dei controlli copiata in ``materiale``) resta quella dello
    studio, con il nome che ha nei metadati.
    """
    tabella = Tabella.dalla_descrizione(descrizione)
    selezionabili = [c["nome"] for c in descrizione["colonne"]
                     if not _non_selezionabile(c) and c["nome"] in tabella.valori]
    selezionabili.sort(key=lambda nome: nome in COLONNE_INVENTARIO)
    return proponi(tabella, selezionabili)


def divisioni_confuse(proposte: Proposte) -> list[str]:
    """Gli avvisi sulle divisioni proposte confuse fra loro (:func:`confuse`),
    coppia per coppia nell'ordine delle proposte.

    Gli incroci non contano, e nemmeno le coppie di due colonne cosi' come
    sono: quelle le avvisa gia' il controllo delle colonne candidate.
    """
    semplici = [p for p in proposte.proposte if p.applicazione.componenti is None]
    righe = [{p.nome: g or "" for p, g in zip(semplici, gruppi)}
             for gruppi in zip(*(p.applicazione.gruppi_dei_campioni for p in semplici))]

    def colonna(p) -> bool:
        return p.origine == "generica" and p.applicazione.espansione == {"colonna": p.nome}

    avvisi = []
    for i, a in enumerate(semplici):
        for b in semplici[i + 1:]:
            if not (colonna(a) and colonna(b)):
                avvisi += confuse(righe, [a.nome, b.nome])
    return avvisi


def nome_del_tsv(uscita: Path | str) -> Path:
    """Il file dei gruppi delle divisioni proposte, accanto alla configurazione:
    ``<nome della configurazione>_divisioni_proposte.tsv``."""
    uscita = Path(uscita)
    return uscita.with_name(f"{uscita.stem}_divisioni_proposte.tsv")


#: Le colonne del file dei gruppi delle divisioni proposte.
COLONNE_DEL_TSV: Final = ("divisione", "regola", "gruppo", "campioni", "campioni_senza_gruppo")


def tsv_delle_divisioni(proposte: Proposte) -> str:
    """Il file dei gruppi delle divisioni proposte: una riga per gruppo, con la
    divisione, la regola, il gruppo, i suoi campioni e i campioni che la
    divisione lascia senza gruppo; nell'ordine delle proposte e dei gruppi."""
    righe = ["\t".join(COLONNE_DEL_TSV)]
    for proposta in proposte.proposte:
        app = proposta.applicazione
        for gruppo, campioni in app.gruppi:
            righe.append("\t".join(
                " ".join(str(x).split()) for x in
                (app.nome, app.regola, gruppo, campioni, app.senza_gruppo)))
    return "\n".join(righe) + "\n"


def blocco_delle_divisioni(proposte: Proposte, ranghi: list[dict[str, Any]],
                           tsv: str | None) -> list[str]:
    """Le righe di commento che stanno sopra ``design``: quante divisioni sono
    suggerite, una riga per ciascuna (nome e spiegazione), come sceglierne una,
    dove sono i gruppi (``tsv``), e i ranghi tassonomici."""
    def paragrafo(testo: str) -> list[str]:
        return [f"# {r}" for r in textwrap.wrap(testo, 76, break_long_words=False,
                                                 break_on_hyphens=False)]

    elenco = proposte.proposte
    larghezza = max((len(p.nome) for p in elenco), default=0)
    righe = [f"# Divisioni suggerite dai metadati: {len(elenco)}"]
    righe += [f"#   {p.nome.ljust(larghezza)}  {p.spiegazione}" for p in elenco]
    if elenco:
        righe += paragrafo(
            "Per sceglierne una scrivi il suo nome dove va una colonna: design.variable, "
            "design.other_variables, comp.group_by, ord.color_by, comp.facet_wrap, le righe "
            "e le colonne di comp.facet_grid; un incrocio («A per B») e' il valore di "
            "comp.facet_grid. Lo strumento ricava la regola dal suo catalogo delle divisioni. "
            "Una colonna dell'oggetto con valori vince su una divisione con lo stesso nome; "
            "una divisione personalizzata si dichiara in design.groupings.")
        if tsv:
            righe += paragrafo(f"I gruppi di ogni divisione, con i campioni: {tsv}.")
    else:
        righe += paragrafo(
            "Nessuna divisione da' da 2 a 12 gruppi, almeno due con "
            f"{DIMENSIONE_DELLE_PROPOSTE} campioni, per almeno meta' dei campioni: si nomina una "
            "colonna dell'oggetto, o si dichiara un gruppo derivato in design.groupings.")
    righe += paragrafo(
        "Ranghi tassonomici (comp.rank, run.analysis_rank; il livello tassonomico non e' "
        "una divisione dei metadati): " + (", ".join(r["rango"] for r in ranghi) or "nessuno")
        + ".")
    return righe


def descrivi(oggetto: Path | str) -> dict[str, Any]:
    """Che cosa contiene l'oggetto: campioni, ranghi, albero, colonne.

    Solleva :class:`~amplicon16s_eco.catalogo.ErroreEco` se l'oggetto non e'
    utilizzabile (``E-ECO-04``) o se R non si esegue (``E-ECO-90``).
    """
    oggetto = Path(oggetto)
    impronta = _impronta_oggetto(oggetto)
    with tempfile.TemporaryDirectory(prefix="amplicon16s-eco-") as temporanea:
        lavoro = Path(temporanea)
        return descrizione_della_copia(_copia(oggetto, lavoro, impronta), lavoro)


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
    # Una sola divisione suggerita (gli incroci non contano), ricavata dalla
    # sola colonna candidata, e' la proposta: e' la sola divisione dei metadati
    # che funziona sul dataset, e nessun'altra colonna e' in gioco. Altrimenti
    # vale la regola sulle colonne.
    delle_divisioni = proposte_di(descrizione)
    semplici = [p for p in delle_divisioni.proposte if p.applicazione.componenti is None]
    if len(semplici) == 1 and all(c["nome"] == semplici[0].applicazione.colonna
                                  for c in candidate):
        metti("design.variable", PROPOSTO, semplici[0].nome,
              f"la sola divisione suggerita dai metadati ({semplici[0].spiegazione})",
              *(["colonne con gruppi ripetuti (numero di gruppi): " + ", ".join(
                  f"{c['nome']} ({c['gruppi_ripetuti']})" for c in candidate)]
                if candidate else []))
    elif len(candidate) == 1:
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
        # I gruppi di ogni divisione, con i loro campioni, sono nel file dei
        # gruppi delle divisioni proposte, non qui.
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
    # Qualunque rango si puo' scegliere: il piu' fine proposto puo' essere troppo
    # fine per lo studio.
    voce = voci["comp.rank"]
    voci["comp.rank"] = replace(voce, note=(*voce.note, (
        "si puo' scegliere qualunque rango dell'oggetto, anche meno specifico se il "
        "proposto e' troppo fine, oppure un elenco di ranghi (per esempio i due piu' "
        "fini), uno per sottocartella di composizione/. Ranghi "
        "dell'oggetto (quota di varianti assegnate): "
        + ", ".join(f"{x['rango']} {x['quota_assegnata']:.2f}" for x in ranghi))))
    if albero["presente"] and albero["radicato"] and albero["con_lunghezze"]:
        metti("beta.unifrac_root", PROPOSTO, "existing",
              "l'albero dell'oggetto e' gia' radicato; serve solo se beta.distances "
              "comprende una UniFrac, ed e' respinto altrimenti")

    # Le divisioni confuse fra loro: si guardano sulla tabella della
    # descrizione, che c'e' sempre; le coppie di colonne le ha gia' viste il
    # controllo delle colonne candidate.
    for avviso in divisioni_confuse(delle_divisioni):
        fatti.append("Divisioni confuse: " + avviso)
    if "design.variable" in voci:
        voce = voci["design.variable"]
        voci["design.variable"] = replace(voce, note=(*voce.note, _RIMANDO))
    metti("alpha.method", PROPOSTO, "rarefy",
          "la rarefazione e' il modo di controllare la profondita' che Schloss 2024 "
          "(mSphere 9:e00354-23) trova piu' affidabile per alfa e beta diversita'",
          "le due strade: rarefy rarefa i conteggi (alpha.rarefy_depth, "
          "alpha.rarefy_iterations: con piu' rarefazioni ogni indice e' la media); none "
          "calcola gli indici sui conteggi non rarefatti, di norma con un filtro delle "
          "varianti rare (filter.min_relative_abundance), e riporta in "
          "alfa/profondita.tsv se gli indici dipendono ancora dalle letture. McMurdie e "
          "Holmes 2014 (PLoS Comput Biol 10:e1003531) giudicano la rarefazione "
          "inammissibile per l'abbondanza differenziale; Weiss et al. 2017 (Microbiome "
          "5:27) concludono che dipende dai dati",
          "con none, se la profondita' e' molto diversa fra i campioni, controlla "
          "l'avviso E-ECO-39")
    metti("run.seed", SCELTA, SENZA_VALORE,
          "un intero qualunque, da dichiarare: fissa la rarefazione e le permutazioni")
    metti("comp.top_n", SCELTA, SENZA_VALORE,
          "quanti taxa mostrare nei grafici prima di riunire gli altri")
    metti("beta.distances", SCELTA, SENZA_VALORE,
          "fra bray, jaccard, aitchison" + (
              ", unifrac_weighted, unifrac_unweighted" if albero["presente"] else
              " (le due UniFrac richiedono un albero, che l'oggetto non ha)"),
          "bray (Bray-Curtis): le differenze di abbondanza relativa, pesate dai taxa "
          "abbondanti",
          "jaccard: la sola presenza o assenza dei taxa",
          "aitchison: la distanza euclidea dei log-rapporti (dati composizionali); chiede "
          "beta.clr_pseudocount",
          *(["unifrac_weighted (UniFrac pesato): la parentela filogenetica dei taxa, con le "
             "loro abbondanze; chiede beta.unifrac_root",
             "unifrac_unweighted (UniFrac non pesato): la parentela filogenetica dei taxa "
             "presenti, senza abbondanze; chiede beta.unifrac_root"]
            if albero["presente"] else []))
    metti("beta.clr_pseudocount", SCELTA, SENZA_VALORE,
          "solo con aitchison: lo pseudoconteggio sommato ai conteggi prima del CLR")
    metti("ord.methods", SCELTA, SENZA_VALORE,
          "fra pcoa, nmds e tsne; [] per nessuna ordinazione. Il t-SNE e' una figura "
          "descrittiva: conserva i vicini, non le distanze")
    metti("ord.distances", SCELTA, SENZA_VALORE, "fra quelle di beta.distances")
    metti("ord.pcoa_correction", SCELTA, "none",
          "solo con pcoa. Partenza: il predefinito di ape::pcoa (correction = \"none\"); "
          "cailliez o lingoes correggono gli autovalori negativi delle distanze non euclidee")
    metti("ord.nmds_trymax", SCELTA, 20,
          "solo con nmds. Partenza: il predefinito di vegan::metaMDS (trymax = 20)")
    metti("ord.tsne_perplexity", SCELTA, SENZA_VALORE,
          "solo con tsne: circa il numero di vicini efficaci di ogni campione, di solito da "
          "5 a 50; Rtsne richiede 3 * perplessita' <= campioni - 1: con "
          f"{descrizione['campioni']} campioni (prima di design.subset e del filtro) al "
          f"massimo {(descrizione['campioni'] - 1) / 3:.1f}")
    metti("constrained.max_df_fraction", SCELTA, SENZA_VALORE,
          "solo con constrained.method: oltre questa frazione dei campioni i gradi di "
          "liberta' del modello danno un avviso, perche' il modello si avvicina a quello "
          "saturo e spiega molto per costruzione. Una variabile con un livello per pochi "
          "campioni (una posizione) ne usa tanti: una frazione bassa, per esempio 0.1, "
          "avvisa presto; la scelta e' dello studio")
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
        if voce.valore is SENZA_VALORE or chiave in _MAI_ACCETTATE:
            pass
        elif chiave in CONDIZIONALI:
            voce = replace(voce, note=(
                *voce.note,
                f"non attivato da --accetta-proposte: va dichiarato solo se "
                f"{condizione(chiave)}, ed e' respinto altrimenti"))
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
    voci: dict[str, Voce], fatti: list[str], oggetto: str, accettate: bool = False,
    divisioni: list[str] | None = None,
) -> str:
    """Il testo del file di configurazione delle analisi; ``divisioni`` e' il
    blocco delle divisioni suggerite (:func:`blocco_delle_divisioni`), righe di
    commento gia' composte, scritto sopra ``design``."""
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
    if divisioni:
        righe += ["# " + "-" * 76, *divisioni, "# " + "-" * 76, ""]
    for gruppo, campo_gruppo in ConfigEco.model_fields.items():
        modello = campo_gruppo.annotation
        righe.append(f"# {(modello.__doc__ or '').strip().splitlines()[0]}")
        righe.append(f"{gruppo}:")
        for nome, campo in modello.model_fields.items():
            chiave = f"{gruppo}.{nome}"
            righe += commento(campo.description or "")
            voce = voci.get(chiave)
            if chiave in CONDIZIONALI:
                righe += commento(f"[CONDIZIONALE] solo se {condizione(chiave)}")
                if chiave in STANDARD_DEI_CONDIZIONALI:
                    valore, fonte = STANDARD_DEI_CONDIZIONALI[chiave]
                    righe += commento(f"[STANDARD {valore_yaml(valore)}] {fonte}; "
                                      "vale se il parametro e' omesso")
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
                for nota in NOTE_DEI_PARAMETRI.get(chiave, ()):
                    righe += commento(nota)
                valore = {} if isinstance(predefinito, dict) else (
                    list(predefinito) if isinstance(predefinito, tuple) else predefinito)
                scritto = "{}" if valore == {} else valore_yaml(valore)
                righe.append(f"  # {nome}: {scritto}")
                continue
            righe += commento(f"[{voce.specie.value}] {voce.origine}")
            for nota in voce.note:
                righe += commento(nota)
            if voce.valore is SENZA_VALORE:
                scritto = ""
            elif isinstance(voce.valore, dict):
                # Una mappa su una riga, in forma JSON: e' YAML valido.
                scritto = " " + json.dumps(voce.valore, ensure_ascii=False)
            else:
                scritto = f" {valore_yaml(voce.valore)}"
            righe.append(f"  {'' if voce.attiva else '# '}{nome}:{scritto}")
        righe.append("")
    return "\n".join(righe)


def _fatto_delle_derivate(descrizione: dict[str, Any],
                          derivazioni: dict[str, str] | None) -> str:
    """Il fatto sulle colonne che la pipeline ha derivato da una regola."""
    if derivazioni is None:
        return ("Colonne derivate dalla pipeline: non note (accanto all'oggetto non ci sono i "
                "file dell'esecuzione: 10_phyloseq/colonne_metadati.tsv).")
    parti = []
    for nome, regola in derivazioni.items():
        vuota = next((c["con_valore"] == 0 for c in descrizione["colonne"]
                      if c["nome"] == nome), False)
        parti.append(f"{nome}: {regola}" + ("; vuota in questo oggetto" if vuota else ""))
    return "Colonne derivate dalla pipeline, con la regola: " + ("; ".join(parti) or "nessuna") + "."


def scrivi(
    oggetto: Path | str, uscita: Path | str, accettate: bool = False
) -> tuple[dict[str, Voce], list[str], Proposte]:
    """Esamina l'oggetto e scrive la configurazione di partenza e, accanto, il
    file dei gruppi delle divisioni proposte (:func:`nome_del_tsv`).

    Con ``accettate`` le proposte e i valori di partenza sono righe attive
    (:func:`accetta`). Solleva ``FileExistsError`` se uno dei due file esiste
    gia': non li sovrascrive.
    """
    uscita = Path(uscita)
    tsv = nome_del_tsv(uscita)
    for percorso in (uscita, tsv):
        if percorso.exists():
            raise FileExistsError(percorso)
    descrizione = descrivi(oggetto)
    voci, fatti = prepara(descrizione, esportati(oggetto))
    fatti.append(_fatto_delle_derivate(descrizione, derivazioni_della_pipeline(oggetto)))
    if accettate:
        voci = accetta(voci)
    proposte = proposte_di(descrizione)
    blocco = blocco_delle_divisioni(proposte, descrizione["ranghi"], tsv.name)
    scrivi_atomico(tsv, tsv_delle_divisioni(proposte))
    scrivi_atomico(uscita, testo_yaml(voci, fatti, Path(oggetto).name, accettate, blocco)
                   + "\n")
    return voci, fatti, proposte
