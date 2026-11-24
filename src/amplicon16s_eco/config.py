"""Configurazione delle analisi ecologiche: schema, caricamento e digest.

Il file e' proprio del pacchetto, distinto da quello della pipeline, e ha
gruppi ``design``, ``run``, ``filter``, ``alpha``, ``comp``, ``beta``, ``ord``,
``stat``, ``constrained`` e ``plot``. Le regole sono quelle della pipeline:

* cio' che dipende dallo studio e' **obbligatorio** e non ha un valore
  predefinito (:data:`OBBLIGATORI`), a cominciare dalla variabile biologica;
* un predefinito esiste solo se e' quello del pacchetto R che realizza il
  metodo, con la fonte (:data:`STANDARD_DEL_METODO`);
* i parametri facoltativi hanno per predefinito l'assenza (:data:`FACOLTATIVI`);
* un parametro che serve a un solo metodo e' obbligatorio quando quel metodo e'
  richiesto e respinto quando non lo e' (:data:`CONDIZIONALI`): un valore
  dichiarato e mai usato entrerebbe nel digest senza effetto. Un condizionale
  puo' avere un predefinito, quello del pacchetto R del metodo
  (:data:`STANDARD_DEI_CONDIZIONALI`): vale solo quando il metodo e' richiesto;
* le chiavi sconosciute sono respinte.

La validazione qui riguarda il solo file. Cio' che dipende dall'oggetto
(colonne, gruppi, ranghi, albero) si verifica in R, prima di ogni calcolo
(:mod:`amplicon16s_eco.esecuzione`).

**Chi lo usa.** La riga di comando (``validate``, ``run``, ``schema``),
l'esecuzione, che passa a R i parametri e registra il digest, ed
``eco init``. Legge il file YAML; non scrive.
"""

from __future__ import annotations

import hashlib
import json
import re
from collections.abc import Callable
from pathlib import Path
from typing import Annotated, Any, Final, Literal

import yaml
from pydantic import BaseModel, ConfigDict, Field, ValidationError, model_validator

from amplicon16s_eco.catalogo import ErroreEco

__all__ = [
    "ASSENTI_SE_NON_DICHIARATI",
    "AMBITI_DEL_FILTRO",
    "CONDIZIONALI",
    "CONDIZIONALI_ANCHE",
    "CORREZIONI_PCOA",
    "DISTANZE",
    "FACOLTATIVI",
    "ORDINAZIONI",
    "OBBLIGATORI",
    "RADICI_UNIFRAC",
    "STANDARD_DEL_METODO",
    "STANDARD_DEI_CONDIZIONALI",
    "ConfigEco",
    "attivo",
    "carica",
    "condizione",
    "digest",
    "risolta",
    "valida",
]

#: Le distanze realizzate. Le due UniFrac richiedono l'albero nell'oggetto.
DISTANZE: Final = ("bray", "jaccard", "aitchison", "unifrac_weighted", "unifrac_unweighted")

Distanza = Literal["bray", "jaccard", "aitchison", "unifrac_weighted", "unifrac_unweighted"]

#: Come si ottiene la radice dell'albero per le due UniFrac: al punto medio
#: (phangorn::midpoint), oppure quella che l'albero dell'oggetto ha gia'.
RADICI_UNIFRAC: Final = ("midpoint", "existing")

#: I metodi di ordinazione realizzati.
ORDINAZIONI: Final = ("pcoa", "nmds", "tsne")

Ordinazione = Literal["pcoa", "nmds", "tsne"]

#: Le correzioni degli autovalori negativi della PCoA (ape::pcoa).
CORREZIONI_PCOA: Final = ("none", "cailliez", "lingoes")

#: Come si misura la soglia del filtro di abbondanza relativa: in ogni campione,
#: sul totale delle letture dei campioni analizzati, o sulla media delle
#: abbondanze relative fra i campioni.
AMBITI_DEL_FILTRO: Final = ("sample", "total", "mean")

#: Come si calcola l'alfa diversita': sui conteggi rarefatti o su quelli non
#: rarefatti (dopo l'eventuale filtro).
METODI_ALFA: Final = ("rarefy", "none")

#: Parametri senza predefinito: dipendono dallo studio o sono scelte di analisi.
OBBLIGATORI: Final[tuple[str, ...]] = (
    "design.variable",
    "run.seed",
    "alpha.method",
    "comp.rank",
    "comp.top_n",
    "beta.distances",
    "ord.methods",
    "ord.distances",
    "stat.min_group_size",
    "stat.significance_level",
)

#: Parametri obbligatori quando uno dei metodi indicati e' richiesto, respinti
#: altrimenti: chiave -> (parametro che indica i metodi, metodi). Il parametro
#: dei metodi e' un elenco o un solo valore; metodi ``None`` vuol dire «quando
#: il parametro e' dichiarato» (non nullo).
CONDIZIONALI: Final[dict[str, tuple[str, tuple[str, ...] | None]]] = {
    "alpha.rarefy_depth": ("alpha.method", ("rarefy",)),
    "alpha.rarefy_iterations": ("alpha.method", ("rarefy",)),
    "filter.scope": ("filter.min_relative_abundance", None),
    "beta.clr_pseudocount": ("beta.distances", ("aitchison",)),
    "beta.unifrac_root": ("beta.distances", ("unifrac_weighted", "unifrac_unweighted")),
    "ord.pcoa_correction": ("ord.methods", ("pcoa",)),
    "ord.nmds_trymax": ("ord.methods", ("nmds",)),
    "ord.tsne_perplexity": ("ord.methods", ("tsne",)),
    "ord.tsne_max_iter": ("ord.methods", ("tsne",)),
    "ord.tsne_theta": ("ord.methods", ("tsne",)),
    "constrained.variables": ("constrained.method", None),
    "constrained.max_df_fraction": ("constrained.method", None),
    "constrained.transform": ("constrained.method", ("rda",)),
    "constrained.distances": ("constrained.method", ("dbrda",)),
}

#: Altri parametri che rendono obbligatorio un condizionale: chiave ->
#: condizioni (parametro, metodi) nella forma di :data:`CONDIZIONALI`. Il
#: condizionale e' richiesto se almeno una condizione e' vera, respinto se
#: nessuna lo e'.
CONDIZIONALI_ANCHE: Final[dict[str, tuple[tuple[str, tuple[str, ...] | None], ...]]] = {
    # Lo pseudoconteggio serve anche alla RDA sui log-rapporti.
    "beta.clr_pseudocount": (("constrained.transform", ("clr",)),),
}

#: Condizionali con un predefinito, il valore del pacchetto R del metodo: vale
#: quando il metodo e' richiesto, e il parametro puo' allora essere omesso.
#: chiave -> (valore, fonte).
STANDARD_DEI_CONDIZIONALI: Final[dict[str, tuple[Any, str]]] = {
    "alpha.rarefy_depth": (
        None,
        "phyloseq::rarefy_even_depth, sample.size = min(sample_sums(physeq)): "
        "nullo significa la minima profondita' fra i campioni analizzati",
    ),
    "alpha.rarefy_iterations": (
        1, "phyloseq::rarefy_even_depth: una sola rarefazione per campione"),
    "ord.tsne_max_iter": (1000, "Rtsne::Rtsne, max_iter = 1000"),
    "ord.tsne_theta": (0.5, "Rtsne::Rtsne, theta = 0.5 (0: t-SNE esatto)"),
}

#: La palette di partenza dei grafici: quella con cui lo strumento li ha sempre
#: disegnati (colori_gruppi in R/lib/uscite.R).
PALETTE_PREDEFINITA: Final = "Dark 3"

#: Predefiniti che sono il valore del pacchetto R che realizza il metodo:
#: chiave -> (valore, fonte).
STANDARD_DEL_METODO: Final[dict[str, tuple[Any, str]]] = {
    "stat.permanova_permutations": (999, "vegan::adonis2, permutations = 999"),
    "plot.palette": (
        PALETTE_PREDEFINITA,
        "grDevices::hcl.colors(n, palette = \"Dark 3\"): la palette con cui lo strumento "
        "ha sempre colorato i gruppi (colori_gruppi, R/lib/uscite.R)"),
}

#: Parametri facoltativi il cui predefinito e' l'assenza.
FACOLTATIVI: Final[tuple[str, ...]] = (
    "design.other_variables",
    "design.technical_variables",
    "design.subset",
    "design.continuous_variables",
    "design.subject",
    "design.groupings",
    "filter.min_relative_abundance",
    "ord.color_by",
    "comp.group_by",
    "comp.facet_wrap",
    "comp.facet_grid",
    "run.analysis_rank",
    "constrained.method",
    "constrained.conditions",
    "stat.permanova_strata",
)

#: I facoltativi che, non dichiarati, non compaiono nella configurazione
#: risolta: ne' nel digest, ne' in ``configurazione.json``, ne' fra i parametri
#: che arrivano a R. Una configurazione che non li dichiara resta cosi' la
#: stessa configurazione, con lo stesso digest e gli stessi file, qualunque
#: analisi facoltativa il pacchetto realizzi.
ASSENTI_SE_NON_DICHIARATI: Final[tuple[str, ...]] = (
    "design.continuous_variables",
    "design.subject",
    "design.groupings",
    "ord.color_by",
    "comp.group_by",
    "comp.facet_wrap",
    "comp.facet_grid",
    "run.analysis_rank",
    "filter.min_relative_abundance",
    "filter.scope",
    "ord.tsne_perplexity",
    "constrained.method",
    "constrained.variables",
    "constrained.conditions",
    "constrained.transform",
    "constrained.distances",
    "constrained.max_df_fraction",
)

#: I gruppi che, senza alcun parametro dichiarato, non compaiono nella
#: configurazione risolta.
_GRUPPI_FACOLTATIVI: Final = ("filter", "constrained")

_Nome = Annotated[str, Field(min_length=1, strict=True)]

#: Il massimo intero che R rappresenta: oltre, un parametro arriverebbe come NA.
_MAX_INTERO: Final = 2**31 - 1
_Positivo = Annotated[int, Field(ge=1, le=_MAX_INTERO, strict=True)]


class _Gruppo(BaseModel):
    """Base dei gruppi: immutabile, chiavi sconosciute respinte."""

    model_config = ConfigDict(extra="forbid", frozen=True)


class Raggruppamento(_Gruppo):
    """Un gruppo derivato: una nuova variabile categorica ricavata da una
    colonna, con un elenco di valori per gruppo oppure con un'espressione
    regolare il cui primo gruppo di cattura da' il nome del gruppo."""

    da: _Nome = Field(alias="from", description=(
        "La colonna della tabella dei campioni da cui si ricava il gruppo."))
    groups: dict[_Nome, tuple[_Nome, ...]] | None = Field(None, description=(
        "Per ogni gruppo i valori della colonna che vi rientrano."))
    pattern: _Nome | None = Field(None, description=(
        "Espressione regolare (PCRE) sul valore della colonna: il primo gruppo di "
        "cattura e' il nome del gruppo."))

    @model_validator(mode="after")
    def _una_regola(self) -> "Raggruppamento":
        """Respinge una regola assente o doppia, un gruppo senza valori, un valore
        in due gruppi, un'espressione non valida o senza gruppo di cattura."""
        if (self.groups is None) == (self.pattern is None):
            raise ValueError("dichiara uno solo fra groups e pattern")
        if self.groups is not None:
            if not self.groups:
                raise ValueError("groups non elenca alcun gruppo")
            vuoti = sorted(g for g, valori in self.groups.items() if not valori)
            if vuoti:
                raise ValueError(f"gruppi senza valori: {vuoti}")
            visti: dict[str, str] = {}
            doppi = []
            for gruppo, valori in self.groups.items():
                for valore in valori:
                    if valore in visti and visti[valore] != gruppo:
                        doppi.append(f"{valore!r} in {visti[valore]} e {gruppo}")
                    visti[valore] = gruppo
            if doppi:
                raise ValueError(f"valori in due gruppi: {', '.join(doppi)}")
        if self.pattern is not None:
            try:
                compilata = re.compile(self.pattern)
            except re.error as e:
                raise ValueError(
                    f"pattern {self.pattern!r} non e' un'espressione valida ({e})") from None
            if compilata.groups < 1:
                raise ValueError(
                    f"pattern {self.pattern!r} non ha un gruppo di cattura: il primo gruppo "
                    "fra parentesi da' il nome del gruppo")
        return self


class Design(_Gruppo):
    """Le colonne della tabella dei campioni che descrivono lo studio."""

    #: La variabile biologica: una sola colonna, senza predefinito.
    variable: _Nome = Field(description=(
        "La variabile biologica: una colonna della tabella dei campioni "
        "dell'oggetto, trattata come categorica."))
    #: Altre variabili biologiche, analizzate una alla volta.
    other_variables: tuple[_Nome, ...] = Field((), description=(
        "Altre variabili biologiche, analizzate una alla volta."))
    #: Variabili tecniche (lotto): primi termini dei modelli.
    technical_variables: tuple[_Nome, ...] = Field((), description=(
        "Variabili tecniche (per esempio il lotto di estrazione), nell'ordine in cui "
        "entrano nei modelli: prima la piu' grossolana."))
    #: Sottoinsieme: colonna -> valori ammessi, applicato prima di tutto.
    subset: dict[_Nome, tuple[_Nome, ...]] = Field(default_factory=dict, description=(
        "Sottoinsieme dei campioni: per ogni colonna i valori ammessi. Si applica "
        "prima di tutto il resto."))
    #: Variabili continue: correlazione con l'alfa diversita' e PERMANOVA con
    #: un termine numerico.
    continuous_variables: tuple[_Nome, ...] = Field((), description=(
        "Variabili continue (misure numeriche): per ciascuna la correlazione di "
        "Spearman con l'alfa diversita' e una PERMANOVA con la variabile come termine "
        "numerico. I valori devono essere numeri."))
    #: Colonna del soggetto in un disegno a misure ripetute.
    subject: _Nome | None = Field(None, description=(
        "Colonna del soggetto in un disegno a misure ripetute: con il soggetto ogni "
        "variabile biologica categorica e' analizzata anche entro soggetto (Wilcoxon "
        "dei ranghi con segno o Friedman, PERMANOVA con le permutazioni entro il "
        "soggetto)."))
    #: Gruppi derivati da una colonna: nome -> regola.
    groupings: dict[_Nome, Raggruppamento] = Field(default_factory=dict, description=(
        "Nuove variabili categoriche ricavate da una colonna, con un elenco di valori per "
        "gruppo (groups) oppure con un'espressione regolare (pattern) il cui primo gruppo "
        "di cattura da' il nome del gruppo. Si usano come le colonne."))

    @model_validator(mode="after")
    def _senza_ripetizioni(self) -> "Design":
        """Respinge un disegno in cui una colonna ha piu' di un ruolo.

        Controlli, nell'ordine: nessuna variabile biologica ripetuta; nessuna
        tecnica ripetuta; nessuna colonna insieme biologica e tecnica; nessuna
        continua ripetuta; nessuna continua che sia anche categorica o tecnica;
        il soggetto non e' una variabile dell'analisi (identifica chi e' stato
        campionato piu' volte, non si confronta); ogni colonna del sottoinsieme
        ha almeno un valore ammesso. Una colonna con due ruoli entrerebbe due
        volte nello stesso modello.
        """
        biologiche = (self.variable, *self.other_variables)
        ripetute = sorted({n for n in biologiche if biologiche.count(n) > 1})
        if ripetute:
            raise ValueError(f"variabili biologiche ripetute: {ripetute}")
        tecniche = self.technical_variables
        ripetute = sorted({n for n in tecniche if tecniche.count(n) > 1})
        if ripetute:
            raise ValueError(f"variabili tecniche ripetute: {ripetute}")
        comuni = sorted(set(biologiche) & set(tecniche))
        if comuni:
            raise ValueError(
                f"una colonna non puo' essere insieme biologica e tecnica: {comuni}"
            )
        continue_ = self.continuous_variables
        ripetute = sorted({n for n in continue_ if continue_.count(n) > 1})
        if ripetute:
            raise ValueError(f"variabili continue ripetute: {ripetute}")
        comuni = sorted(set(continue_) & (set(biologiche) | set(tecniche)))
        if comuni:
            raise ValueError(
                "una colonna non puo' essere insieme continua e categorica o tecnica: "
                f"{comuni}; una variabile continua si dichiara solo in "
                "design.continuous_variables"
            )
        if self.subject is not None and self.subject in (*biologiche, *tecniche, *continue_):
            raise ValueError(
                f"design.subject ({self.subject!r}) e' anche una variabile dell'analisi: il "
                "soggetto identifica chi e' stato campionato piu' volte, non si confronta"
            )
        vuoti = sorted(c for c, valori in self.subset.items() if not valori)
        if vuoti:
            raise ValueError(f"design.subset senza valori ammessi per: {vuoti}")
        di_se = sorted(n for n, g in self.groupings.items() if g.da in self.groupings)
        if di_se:
            raise ValueError(
                "design.groupings: un gruppo derivato si ricava da una colonna dell'oggetto, "
                f"non da un altro gruppo derivato: {di_se}")
        return self


class Constrained(_Gruppo):
    """Ordinazione vincolata (supervisionata): RDA e db-RDA."""

    #: Il metodo; nullo: nessuna ordinazione vincolata.
    method: Literal["rda", "dbrda"] | None = Field(None, description=(
        "rda: RDA (vegan::rda) sui conteggi trasformati; dbrda: RDA basata sulle distanze "
        "(vegan::dbrda) sulle distanze dichiarate. Nullo: nessuna ordinazione vincolata."))
    #: Le variabili esplicative.
    variables: Annotated[tuple[_Nome, ...], Field(min_length=1)] | None = Field(None, description=(
        "Variabili esplicative: categoriche, gruppi derivati o continue del disegno. "
        "Obbligatorie con method, respinte senza."))
    #: Le condizioni, tolte prima (Condition()).
    conditions: tuple[_Nome, ...] = Field((), description=(
        "Variabili il cui effetto si toglie prima (Condition()), per esempio quelle "
        "tecniche."))
    #: Trasformazione dei conteggi: solo con rda.
    transform: Literal["hellinger", "clr"] | None = Field(None, description=(
        "Trasformazione dei conteggi per la RDA: hellinger (Legendre e Gallagher 2001) o "
        "clr (con beta.clr_pseudocount). Obbligatoria con rda, respinta altrimenti."))
    #: Distanze: solo con dbrda.
    distances: Annotated[tuple[Distanza, ...], Field(min_length=1)] | None = Field(
        None, description=(
            "Distanze della db-RDA, fra quelle di beta.distances (UniFrac comprese). "
            "Obbligatorie con dbrda, respinte altrimenti."))
    #: Frazione dei campioni oltre cui i gradi di liberta' del modello avvisano.
    max_df_fraction: Annotated[float, Field(gt=0, le=1, strict=True, allow_inf_nan=False)] | None = (
        Field(None, description=(
            "Frazione dei campioni oltre la quale i gradi di liberta' del modello danno un "
            "avviso: il modello si avvicina a quello saturo. Obbligatoria con method, senza "
            "predefinito.")))


class Plot(_Gruppo):
    """Grafici."""

    #: Palette dei gruppi: un nome o un elenco di colori.
    palette: _Nome | Annotated[tuple[_Nome, ...], Field(min_length=1)] = Field(
        PALETTE_PREDEFINITA, description=(
            "Palette dei gruppi in tutti i grafici: un nome di grDevices::hcl.colors, "
            "\"Okabe-Ito\" (grDevices::palette.colors, adatta ai daltonici), un nome di "
            "RColorBrewer o di viridisLite, oppure un elenco di colori esadecimali. Con piu' "
            "gruppi che colori la palette si estende per interpolazione."))

    @model_validator(mode="after")
    def _colori(self) -> "Plot":
        """Un elenco di colori e' fatto di colori esadecimali (#RRGGBB o
        #RRGGBBAA); un nome si verifica in R, che conosce le palette."""
        if isinstance(self.palette, tuple):
            sbagliati = [c for c in self.palette
                         if not re.fullmatch(r"#[0-9A-Fa-f]{6}([0-9A-Fa-f]{2})?", c)]
            if sbagliati:
                raise ValueError(f"colori non esadecimali (#RRGGBB): {sbagliati}")
        return self


class Run(_Gruppo):
    """Parametri dell'esecuzione."""

    #: Seme di rarefazione e di ogni altra operazione casuale: va dichiarato.
    seed: Annotated[int, Field(ge=0, le=_MAX_INTERO, strict=True)] = Field(description=(
        "Seme della rarefazione e di ogni altra operazione casuale."))
    #: Rango dei calcoli; nullo: le varianti.
    analysis_rank: _Nome | None = Field(None, description=(
        "Rango tassonomico a cui agglomerare i conteggi prima di ogni calcolo (alfa, "
        "distanze, ordinazioni, test, metodo vincolato), con phyloseq::tax_glom e "
        "NArm = FALSE. Nullo: le varianti."))


class Filter(_Gruppo):
    """Filtro delle varianti rare per abbondanza relativa."""

    #: Soglia in percentuale; nullo: nessun filtro.
    min_relative_abundance: Annotated[
        float, Field(gt=0, lt=100, strict=True, allow_inf_nan=False)] | None = Field(
        None, description=(
            "Soglia di abbondanza relativa, in percentuale (0.1 = 0,1%). Nullo: nessun "
            "filtro. Si applica dopo design.subset e prima di ogni analisi."))
    #: Su che cosa si misura la soglia.
    scope: Literal["sample", "total", "mean"] | None = Field(None, description=(
        "Su che cosa si misura la soglia: sample azzera in ogni campione le varianti "
        "sotto la soglia delle sue letture; total toglie le varianti la cui somma e' "
        "sotto la soglia delle letture totali (Bokulich 2013); mean toglie le varianti "
        "la cui abbondanza relativa media fra i campioni e' sotto la soglia. "
        "Obbligatorio con min_relative_abundance, respinto altrimenti."))


class Alpha(_Gruppo):
    """Alfa diversita'."""

    #: Rarefazione o conteggi non rarefatti: senza predefinito.
    method: Literal["rarefy", "none"] = Field(description=(
        "rarefy: alfa diversita' su conteggi rarefatti; none: sui conteggi non "
        "rarefatti, dopo l'eventuale filtro, con la correlazione fra letture e indici "
        "in alfa/profondita.tsv."))
    #: Profondita' di rarefazione. Nullo: la minima fra i campioni analizzati
    #: (phyloseq::rarefy_even_depth, sample.size = min(sample_sums(physeq))).
    rarefy_depth: _Positivo | None = Field(None, description=(
        "Profondita' di rarefazione. Nullo: la minima fra i campioni analizzati. Se "
        "dichiarata, i campioni con meno letture escono dall'alfa diversita'. Solo con "
        "alpha.method rarefy."))
    #: Rarefazioni ripetute: con piu' di una, ogni indice e' la media.
    rarefy_iterations: _Positivo | None = Field(None, description=(
        "Numero di rarefazioni (predefinito 1): con piu' di una, ogni indice e' la media "
        "sulle ripetizioni (rarefazione ripetuta). Solo con alpha.method rarefy."))


class Griglia(_Gruppo):
    """Le due colonne di una griglia di pannelli: righe e colonne."""

    rows: _Nome = Field(description="Colonna o gruppo derivato delle righe della griglia.")
    cols: _Nome = Field(description="Colonna o gruppo derivato delle colonne della griglia.")

    @model_validator(mode="after")
    def _distinte(self) -> "Griglia":
        """Righe e colonne sono due colonne diverse."""
        if self.rows == self.cols:
            raise ValueError("comp.facet_grid: rows e cols sono la stessa colonna")
        return self


class Comp(_Gruppo):
    """Composizione tassonomica."""

    rank: _Nome | Annotated[tuple[_Nome, ...], Field(min_length=1)] = Field(description=(
        "Rango tassonomico a cui agglomerare i grafici della composizione, fra quelli "
        "dell'oggetto, oppure un elenco di ranghi: uno per sottocartella di "
        "composizione/."))
    top_n: _Positivo = Field(description=(
        "Numero di taxa mostrati nei grafici; i rimanenti in «altri»."))
    #: Variabili delle barre della composizione; vuoto: quelle del disegno.
    group_by: tuple[_Nome, ...] = Field((), description=(
        "Colonne o gruppi derivati per cui disegnare le barre della composizione (media "
        "per gruppo). Vuoto: le variabili del disegno."))
    #: Barre per campione, un pannello per gruppo.
    facet_wrap: _Nome | None = Field(None, description=(
        "Colonna o gruppo derivato: barre per campione, un pannello per gruppo "
        "(ggplot2::facet_wrap, scale libere sull'asse x)."))
    #: Barre medie per gruppo della variabile, in una griglia.
    facet_grid: Griglia | None = Field(None, description=(
        "Barre medie per gruppo della variabile del disegno in una griglia di pannelli: "
        "righe e colonne sono due colonne o gruppi derivati (ggplot2::facet_grid)."))

    @model_validator(mode="after")
    def _senza_ripetizioni(self) -> "Comp":
        """Respinge una variabile ripetuta in group_by o un rango ripetuto."""
        ripetute = sorted({v for v in self.group_by if self.group_by.count(v) > 1})
        if ripetute:
            raise ValueError(f"comp.group_by con valori ripetuti: {ripetute}")
        if isinstance(self.rank, tuple):
            ripetuti = sorted({v for v in self.rank if self.rank.count(v) > 1})
            if ripetuti:
                raise ValueError(f"comp.rank con ranghi ripetuti: {ripetuti}")
        return self


class Beta(_Gruppo):
    """Beta diversita'."""

    distances: Annotated[tuple[Distanza, ...], Field(min_length=1)] = Field(description=(
        "Distanze da calcolare. Le due UniFrac richiedono un albero nell'oggetto."))
    #: Pseudoconteggio del CLR: solo con la distanza di Aitchison.
    clr_pseudocount: Annotated[float, Field(gt=0, strict=True, allow_inf_nan=False)] | None = Field(
        None, description=(
            "Pseudoconteggio sommato ai conteggi prima del CLR: obbligatorio con la "
            "distanza di Aitchison, respinto altrimenti."))
    #: Radice dell'albero: solo con una delle due UniFrac.
    unifrac_root: Literal["midpoint", "existing"] | None = Field(None, description=(
        "Radice dell'albero per UniFrac: midpoint lo radica al punto medio, existing "
        "usa la radice che ha gia'. Obbligatorio con una UniFrac, respinto altrimenti."))

    @model_validator(mode="after")
    def _senza_ripetizioni(self) -> "Beta":
        """Respinge un elenco di distanze con un valore ripetuto."""
        ripetute = sorted({d for d in self.distances if self.distances.count(d) > 1})
        if ripetute:
            raise ValueError(f"distanze ripetute: {ripetute}")
        return self


class Ord(_Gruppo):
    """Ordinazione."""

    #: Metodi di ordinazione; un elenco vuoto non ne produce.
    methods: tuple[Ordinazione, ...] = Field(description=(
        "Metodi di ordinazione. Un elenco vuoto non produce ordinazioni."))
    #: Distanze da ordinare, fra quelle di ``beta.distances``.
    distances: tuple[Distanza, ...] = Field(description=(
        "Distanze da ordinare, fra quelle di beta.distances."))
    #: Correzione degli autovalori negativi: solo con la PCoA.
    pcoa_correction: Literal["none", "cailliez", "lingoes"] | None = Field(None, description=(
        "Correzione degli autovalori negativi della PCoA (ape::pcoa): obbligatoria "
        "con pcoa, respinta altrimenti."))
    #: Numero di avvii casuali: solo con la NMDS.
    nmds_trymax: _Positivo | None = Field(None, description=(
        "Numero di avvii casuali della NMDS (vegan::metaMDS): obbligatorio con nmds, "
        "respinto altrimenti."))
    #: Perplessita' del t-SNE: solo con tsne, senza predefinito.
    tsne_perplexity: Annotated[float, Field(gt=0, strict=True, allow_inf_nan=False)] | None = Field(
        None, description=(
            "Perplessita' del t-SNE (Rtsne::Rtsne): circa il numero di vicini efficaci di "
            "ogni campione. Obbligatoria con tsne, respinta altrimenti; 3 * perplessita' non "
            "puo' superare il numero di campioni meno uno."))
    #: Iterazioni del t-SNE (Rtsne::Rtsne, max_iter = 1000).
    tsne_max_iter: _Positivo | None = Field(None, description=(
        "Iterazioni del t-SNE (predefinito 1000, Rtsne::Rtsne). Solo con tsne."))
    #: Approssimazione di Barnes-Hut (Rtsne::Rtsne, theta = 0.5).
    tsne_theta: Annotated[float, Field(ge=0, le=1, strict=True, allow_inf_nan=False)] | None = Field(
        None, description=(
            "Compromesso fra velocita' e precisione di Barnes-Hut (predefinito 0.5, "
            "Rtsne::Rtsne); 0 per il t-SNE esatto. Solo con tsne."))
    #: Colonne o gruppi derivati per cui disegnare le ordinazioni senza test.
    color_by: tuple[_Nome, ...] = Field((), description=(
        "Colonne o gruppi derivati per cui disegnare le ordinazioni colorate, senza "
        "alcun test: figure descrittive."))

    @model_validator(mode="after")
    def _coerente(self) -> "Ord":
        """Respinge metodi o distanze ripetuti, e un elenco dichiarato senza
        l'altro: un'ordinazione chiede insieme un metodo e una distanza."""
        for nome in ("methods", "distances"):
            valori = getattr(self, nome)
            ripetuti = sorted({v for v in valori if valori.count(v) > 1})
            if ripetuti:
                raise ValueError(f"ord.{nome} con valori ripetuti: {ripetuti}")
        ripetute = sorted({v for v in self.color_by if self.color_by.count(v) > 1})
        if ripetute:
            raise ValueError(f"ord.color_by con valori ripetuti: {ripetute}")
        if bool(self.methods) != bool(self.distances):
            raise ValueError(
                "ord.methods e ord.distances vanno dichiarati insieme: entrambi "
                "con almeno un valore, oppure entrambi vuoti"
            )
        return self


class Stat(_Gruppo):
    """Gruppi e test."""

    #: Dimensione minima di un gruppo per entrare nei test.
    min_group_size: Annotated[int, Field(ge=2, le=_MAX_INTERO, strict=True)] = Field(
        description="Numero minimo di campioni perche' un gruppo entri nei test.")
    #: Livello sotto il quale le dispersioni si dichiarano diverse.
    significance_level: Annotated[float, Field(gt=0, lt=1, strict=True)] = Field(description=(
        "Livello di significativita' sotto il quale le dispersioni dei gruppi si "
        "dichiarano diverse, con un avviso sulla PERMANOVA."))
    #: Permutazioni di PERMANOVA e dispersioni (vegan::adonis2, permutations = 999).
    permanova_permutations: _Positivo = Field(999, description=(
        "Permutazioni della PERMANOVA e del test delle dispersioni (vegan::adonis2)."))
    #: Colonna degli strati entro cui permutare nella PERMANOVA.
    permanova_strata: _Nome | None = Field(None, description=(
        "Colonna degli strati entro cui permutare nella PERMANOVA e nel test delle "
        "dispersioni."))


class ConfigEco(_Gruppo):
    """La configurazione di un'analisi."""

    design: Design
    run: Run
    filter: Filter = Field(default_factory=Filter)
    alpha: Alpha
    comp: Comp
    beta: Beta
    ord: Ord
    stat: Stat
    constrained: Constrained = Field(default_factory=Constrained)
    plot: Plot = Field(default_factory=Plot)

    @model_validator(mode="after")
    def _distanze_ordinate(self) -> "ConfigEco":
        """Respinge un'ordinazione su una distanza che il gruppo delle distanze non
        calcola."""
        estranee = [d for d in self.ord.distances if d not in self.beta.distances]
        if estranee:
            raise ValueError(
                f"ord.distances nomina distanze che beta.distances non calcola: {estranee}"
            )
        vincolata = self.constrained
        if vincolata.method is None and vincolata.conditions:
            raise ValueError("constrained: conditions dichiarate senza constrained.method")
        if vincolata.variables:
            ammesse = {self.design.variable, *self.design.other_variables,
                       *self.design.continuous_variables, *self.design.groupings}
            estranee = [v for v in vincolata.variables if v not in ammesse]
            if estranee:
                raise ValueError(
                    "constrained: le variabili devono essere variabili del disegno "
                    "(design.variable, other_variables, continuous_variables) o gruppi "
                    f"derivati: {estranee}")
            ripetute = sorted({v for v in (*vincolata.variables, *vincolata.conditions)
                               if (*vincolata.variables, *vincolata.conditions).count(v) > 1})
            if ripetute:
                raise ValueError(
                    f"constrained: variabili ripetute, o insieme variabile e condizione: {ripetute}")
        estranee = [d for d in vincolata.distances or () if d not in self.beta.distances]
        if estranee:
            raise ValueError(
                f"constrained: distances nomina distanze che beta.distances non calcola: {estranee}")
        biologiche = {self.design.variable, *self.design.other_variables}
        doppie = sorted(biologiche & set(self.ord.color_by))
        if doppie:
            raise ValueError(
                "ord.color_by nomina variabili del disegno, che hanno gia' i loro grafici: "
                f"{doppie}")
        return self


# --------------------------------------------------------------------------- #
# Dai problemi dello schema ai codici del catalogo                             #
# --------------------------------------------------------------------------- #

#: Parametri il cui valore e' il nome di un metodo: un valore non ammesso e' un
#: metodo non supportato (E-ECO-03), non un errore generico.
_PARAMETRI_DI_METODO: Final = (
    "alpha.method", "filter.scope", "beta.distances", "beta.unifrac_root", "ord.methods",
    "ord.distances", "ord.pcoa_correction", "constrained.method", "constrained.transform",
    "constrained.distances",
)

_MESSAGGI: Final[dict[str, str]] = {
    "extra_forbidden": "parametro sconosciuto (refuso o parametro non previsto)",
    "int_type": "deve essere un numero intero",
    "float_type": "deve essere un numero",
    "finite_number": "deve essere un numero finito",
    "invalid_key": "le chiavi devono essere testo (fra virgolette, se sono numeri)",
    "string_type": "deve essere una stringa (fra virgolette, se e' un numero)",
    "tuple_type": "deve essere un elenco",
    "dict_type": "deve essere una mappa",
    "model_type": "deve essere una mappa di parametri",
    "greater_than": "deve essere maggiore di {gt}",
    "less_than": "deve essere minore di {lt}",
    "greater_than_equal": "deve essere maggiore o uguale a {ge}",
    "less_than_equal": "deve essere minore o uguale a {le}",
    "too_short": "deve contenere almeno {min_length} elementi",
    "string_too_short": "non puo' essere vuoto",
}


def _chiave(loc: tuple[Any, ...]) -> str:
    """Il nome del parametro con il punto, senza gli indici degli elenchi."""
    return ".".join(str(p) for p in loc if not isinstance(p, int))


def _rifiuto(problema: dict[str, Any]) -> tuple[str, str]:
    """Il codice del catalogo e il dettaglio per un problema dello schema."""
    chiave = _chiave(problema["loc"])
    tipo = problema["type"]
    contesto = problema.get("ctx") or {}
    if chiave.startswith("design.groupings") or (
            tipo == "value_error" and "design.groupings:" in str(contesto.get("error", ""))):
        # Ogni problema di un gruppo derivato ha il suo codice.
        if tipo == "missing":
            testo = "parametro mancante"
        elif tipo == "value_error":
            testo = str(contesto.get("error", problema["msg"]))
        elif tipo in _MESSAGGI:
            testo = _MESSAGGI[tipo].format(**contesto)
        else:
            testo = problema["msg"]
        return "E-ECO-40", f"{chiave}: {testo}."
    if tipo == "value_error" and str(contesto.get("error", "")).startswith("constrained:"):
        return "E-ECO-43", str(contesto["error"]).removeprefix("constrained:").strip() + "."
    if chiave == "plot" or chiave.startswith("plot.palette"):
        testo = (str(contesto.get("error", problema["msg"])) if tipo == "value_error"
                 else "un nome di palette o un elenco di colori esadecimali")
        return "E-ECO-41", f"plot.palette: {testo}."
    if tipo == "missing":
        return "E-ECO-02", f"Parametro mancante: {chiave}."
    if tipo == "literal_error" and chiave in _PARAMETRI_DI_METODO:
        return "E-ECO-03", (
            f"{chiave}: {problema['input']!r} non e' ammesso; ammessi: "
            f"{contesto.get('expected', '?')}."
        )
    if tipo == "value_error":
        testo = str(contesto.get("error", problema["msg"]))
    elif tipo in _MESSAGGI:
        testo = _MESSAGGI[tipo].format(**contesto)
    else:
        testo = problema["msg"]
    dove = f"{chiave}: " if chiave else ""
    return "E-ECO-01", f"{dove}{testo}."


def condizioni(chiave: str) -> tuple[tuple[str, tuple[str, ...] | None], ...]:
    """Le condizioni che rendono richiesto il condizionale ``chiave``."""
    return (CONDIZIONALI[chiave], *CONDIZIONALI_ANCHE.get(chiave, ()))


def condizione(chiave: str) -> str:
    """Quando il condizionale ``chiave`` e' richiesto, in parole: per i commenti
    di ``eco init`` e per i messaggi."""
    parti = []
    for parametro, metodi in condizioni(chiave):
        if metodi is None:
            parti.append(f"{parametro} e' dichiarato")
        else:
            parti.append(f"{parametro} comprende {' o '.join(metodi)}")
    return " oppure ".join(parti)


def _richiesti(parametro_metodi: tuple[str, tuple[str, ...] | None], dichiarati: Any) -> list[str] | None:
    """I metodi della condizione presenti fra quelli dichiarati; ``None`` se la
    dichiarazione non si puo' giudicare (assente o di un tipo che lo schema
    respinge)."""
    _, metodi = parametro_metodi
    if metodi is None:
        if dichiarati is None or dichiarati == [] or dichiarati == ():
            return []
        return ["(dichiarato)"]
    if isinstance(dichiarati, str):
        dichiarati = [dichiarati]
    if dichiarati is None:
        return []
    if not isinstance(dichiarati, (list, tuple)):
        return None
    return [m for m in metodi if m in dichiarati]


def attivo(chiave: str, valore_di: Callable[[str], Any]) -> bool:
    """Vero se almeno una condizione del condizionale ``chiave`` e' vera."""
    return any(_richiesti(c, valore_di(c[0])) for c in condizioni(chiave))


def _dichiarato(valore: Any) -> bool:
    """Un valore conta come dichiarato se non e' nullo ne' un elenco o una mappa
    vuoti."""
    return valore is not None and valore != [] and valore != () and valore != {}


def _condizionali(
    valore_di: Callable[[str], Any], esclusi: frozenset[str] = frozenset()
) -> list[tuple[str, str]]:
    """I rifiuti dei parametri condizionali: mancanti con il metodo richiesto,
    oppure dichiarati senza il metodo che li usa.

    ``valore_di`` da' il valore di un parametro ``gruppo.nome``. Un parametro
    dei metodi non dichiarato, o non leggibile, non richiede nulla: il suo
    rifiuto viene dallo schema. Un condizionale non si giudica se lui o uno dei
    parametri delle sue condizioni sono fra gli ``esclusi``. Un condizionale
    con un predefinito (:data:`STANDARD_DEI_CONDIZIONALI`) non e' mai mancante.
    """
    rifiuti: list[tuple[str, str]] = []
    for chiave in CONDIZIONALI:
        if chiave in esclusi or any(c[0] in esclusi for c in condizioni(chiave)):
            continue
        valore = valore_di(chiave)
        giudizi = [(c, _richiesti(c, valore_di(c[0]))) for c in condizioni(chiave)]
        if any(r is None for _, r in giudizi):
            continue
        veri = [(c, r) for c, r in giudizi if r]
        if veri and not _dichiarato(valore) and chiave not in STANDARD_DEI_CONDIZIONALI:
            (parametro, metodi), richiesti = veri[0]
            motivo = (f"{parametro} e' dichiarato" if metodi is None
                      else f"{parametro} comprende {richiesti[0]}")
            rifiuti.append((
                "E-ECO-02",
                f"Parametro mancante: {chiave}, obbligatorio perche' {motivo}.",
            ))
        if not veri and _dichiarato(valore):
            rifiuti.append((
                "E-ECO-01",
                f"{chiave}: dichiarato, ma non serve: e' ammesso solo se "
                f"{condizione(chiave)}; toglilo, oppure richiedi il metodo che lo usa.",
            ))
    return rifiuti


def valida(dati: Any) -> ConfigEco:
    """Valida i dati letti dal file e restituisce la configurazione.

    Solleva :class:`~amplicon16s_eco.catalogo.ErroreEco` con tutti i problemi
    trovati, ciascuno con il suo codice, non solo con il primo.
    """
    if not isinstance(dati, dict):
        raise ErroreEco([("E-ECO-01", "Il file non contiene una mappa di gruppi di parametri.")])
    # Un gruppo con tutte le chiavi in commento arriva come nullo: vale come un
    # gruppo vuoto, cosi' il rifiuto nomina i parametri che gli mancano.
    dati = {gruppo: {} if valori is None else valori for gruppo, valori in dati.items()}
    try:
        config = ConfigEco.model_validate(dati)
    except ValidationError as e:
        problemi = e.errors()
        # Un elenco i cui elementi sono tutti respinti risulta anche vuoto: il
        # secondo problema e' una conseguenza del primo e non si riporta.
        respinti = {_chiave(p["loc"]) for p in problemi if p["type"] != "too_short"}
        # I condizionali che mancano ai metodi gia' dichiarati si chiedono nello
        # stesso giro degli altri parametri mancanti, sui dati letti; quelli di
        # un elenco di metodi non ancora dichiarato, o respinto dallo schema,
        # non si possono conoscere. Un condizionale dichiarato senza il suo
        # metodo si giudica solo sulla configurazione valida: qui sarebbe la
        # conseguenza di un elenco respinto.
        # Gia' respinti dallo schema: il condizionale stesso, o il suo elenco.
        gia_detti = {_chiave(p["loc"]) for p in problemi}

        def dai_dati(chiave: str) -> Any:
            """Il valore dichiarato nel file per la chiave ``gruppo.parametro``, o ``None``."""
            gruppo, nome = chiave.split(".")
            valori = dati.get(gruppo)
            return valori.get(nome) if isinstance(valori, dict) else None

        # Un valore che lo schema ammette in piu' forme (la palette) puo' dare lo
        # stesso rifiuto una volta per forma: si riporta una volta sola.
        raise ErroreEco(list(dict.fromkeys([
            _rifiuto(p) for p in problemi
            if not (p["type"] == "too_short" and _chiave(p["loc"]) in respinti)
        ] + [r for r in _condizionali(dai_dati, frozenset(gia_detti))
             if r[0] == "E-ECO-02"]))) from None

    def dal_modello(chiave: str) -> Any:
        """Il valore della chiave ``gruppo.parametro`` nella configurazione validata."""
        gruppo, nome = chiave.split(".")
        return getattr(getattr(config, gruppo), nome)

    rifiuti = _condizionali(dal_modello)
    if rifiuti:
        raise ErroreEco(rifiuti)
    return config


class _SenzaChiaviRipetute(yaml.SafeLoader):
    """Lettore YAML che respinge una chiave dichiarata due volte nella stessa
    mappa: il lettore normale terrebbe l'ultima in silenzio.
    """


def _mappa_senza_ripetizioni(lettore: yaml.SafeLoader, nodo: yaml.MappingNode) -> dict[Any, Any]:
    """Costruisce una mappa YAML dopo avere respinto le chiavi dichiarate piu'
    volte, che il lettore di base terrebbe in silenzio con l'ultimo valore."""
    chiavi = [lettore.construct_object(chiave, deep=True) for chiave, _ in nodo.value]
    ripetute = sorted({str(c) for c in chiavi if chiavi.count(c) > 1})
    if ripetute:
        raise yaml.YAMLError(f"chiavi dichiarate piu' volte: {', '.join(ripetute)}")
    return lettore.construct_mapping(nodo, deep=True)


_SenzaChiaviRipetute.add_constructor(
    yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG, _mappa_senza_ripetizioni
)


def carica(percorso: Path | str) -> ConfigEco:
    """Legge e valida il file di configurazione dell'analisi."""
    percorso = Path(percorso)
    try:
        dati = yaml.load(percorso.read_text(encoding="utf-8"), Loader=_SenzaChiaviRipetute)  # noqa: S506
    except (OSError, UnicodeDecodeError, yaml.YAMLError) as e:
        raise ErroreEco([("E-ECO-01", f"File non leggibile: {percorso} ({e}).")]) from None
    return valida(dati)


def risolta(config: ConfigEco) -> dict[str, Any]:
    """La configurazione risolta: ogni parametro con il valore che vale,
    predefiniti compresi, in tipi JSON. I facoltativi di
    :data:`ASSENTI_SE_NON_DICHIARATI` vi compaiono solo se dichiarati; un
    condizionale con un predefinito vi compare, con il predefinito, solo se il
    metodo che lo usa e' richiesto; un gruppo facoltativo senza parametri non vi
    compare.
    """
    dati = config.model_dump(mode="json", by_alias=True)

    def valore_di(chiave: str) -> Any:
        """Il valore della chiave ``gruppo.parametro`` nella configurazione."""
        gruppo, nome = chiave.split(".")
        return dati.get(gruppo, {}).get(nome)

    for chiave, (predefinito, _) in STANDARD_DEI_CONDIZIONALI.items():
        gruppo, nome = chiave.split(".")
        if not attivo(chiave, valore_di):
            del dati[gruppo][nome]
        elif dati[gruppo][nome] is None:
            dati[gruppo][nome] = predefinito
    for chiave in ASSENTI_SE_NON_DICHIARATI:
        gruppo, nome = chiave.split(".")
        if nome in dati[gruppo] and not _dichiarato(dati[gruppo][nome]):
            del dati[gruppo][nome]
    for gruppo in _GRUPPI_FACOLTATIVI:
        if not dati[gruppo]:
            del dati[gruppo]
    return dati


def digest(config: ConfigEco) -> str:
    """SHA-256 della configurazione risolta in forma canonica (chiavi
    ordinate): dipende dai valori, non dai byte del file.
    """
    canonico = json.dumps(risolta(config), sort_keys=True, ensure_ascii=False,
                          separators=(",", ":"))
    return hashlib.sha256(canonico.encode("utf-8")).hexdigest()
