"""Configurazione delle analisi ecologiche: schema, caricamento e digest.

Il file e' proprio del pacchetto, distinto da quello della pipeline, e ha
gruppi ``design``, ``run``, ``alpha``, ``comp``, ``beta``, ``ord`` e ``stat``. Le regole
sono quelle della pipeline:

* cio' che dipende dallo studio e' **obbligatorio** e non ha un valore
  predefinito (:data:`OBBLIGATORI`), a cominciare dalla variabile biologica;
* un predefinito esiste solo se e' quello del pacchetto R che realizza il
  metodo, con la fonte (:data:`STANDARD_DEL_METODO`);
* i parametri facoltativi hanno per predefinito l'assenza (:data:`FACOLTATIVI`);
* un parametro che serve a un solo metodo e' obbligatorio quando quel metodo e'
  richiesto e respinto quando non lo e' (:data:`CONDIZIONALI`): un valore
  dichiarato e mai usato entrerebbe nel digest senza effetto;
* le chiavi sconosciute sono respinte.

La validazione qui riguarda il solo file. Cio' che dipende dall'oggetto
(colonne, gruppi, ranghi, albero) si verifica in R, prima di ogni calcolo
(:mod:`amplicon16s_eco.esecuzione`).
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Annotated, Any, Final, Literal

import yaml
from pydantic import BaseModel, ConfigDict, Field, ValidationError, model_validator

from amplicon16s_eco.catalogo import ErroreEco

__all__ = [
    "CONDIZIONALI",
    "CORREZIONI_PCOA",
    "DISTANZE",
    "FACOLTATIVI",
    "ORDINAZIONI",
    "OBBLIGATORI",
    "RADICI_UNIFRAC",
    "STANDARD_DEL_METODO",
    "ConfigEco",
    "carica",
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
ORDINAZIONI: Final = ("pcoa", "nmds")

Ordinazione = Literal["pcoa", "nmds"]

#: Le correzioni degli autovalori negativi della PCoA (ape::pcoa).
CORREZIONI_PCOA: Final = ("none", "cailliez", "lingoes")

#: Parametri senza predefinito: dipendono dallo studio o sono scelte di analisi.
OBBLIGATORI: Final[tuple[str, ...]] = (
    "design.variable",
    "run.seed",
    "comp.rank",
    "comp.top_n",
    "beta.distances",
    "ord.methods",
    "ord.distances",
    "stat.min_group_size",
    "stat.significance_level",
)

#: Parametri obbligatori quando uno dei metodi indicati e' richiesto, respinti
#: altrimenti: chiave -> (parametro che elenca i metodi, metodi).
CONDIZIONALI: Final[dict[str, tuple[str, tuple[str, ...]]]] = {
    "beta.clr_pseudocount": ("beta.distances", ("aitchison",)),
    "beta.unifrac_root": ("beta.distances", ("unifrac_weighted", "unifrac_unweighted")),
    "ord.pcoa_correction": ("ord.methods", ("pcoa",)),
    "ord.nmds_trymax": ("ord.methods", ("nmds",)),
}

#: Predefiniti che sono il valore del pacchetto R che realizza il metodo:
#: chiave -> (valore, fonte).
STANDARD_DEL_METODO: Final[dict[str, tuple[Any, str]]] = {
    "alpha.rarefy_depth": (
        None,
        "phyloseq::rarefy_even_depth, sample.size = min(sample_sums(physeq)): "
        "nullo significa la minima profondita' fra i campioni analizzati",
    ),
    "stat.permanova_permutations": (999, "vegan::adonis2, permutations = 999"),
}

#: Parametri facoltativi il cui predefinito e' l'assenza.
FACOLTATIVI: Final[tuple[str, ...]] = (
    "design.other_variables",
    "design.technical_variables",
    "design.subset",
    "stat.permanova_strata",
)

_Nome = Annotated[str, Field(min_length=1, strict=True)]

#: Il massimo intero che R rappresenta: oltre, un parametro arriverebbe come NA.
_MAX_INTERO: Final = 2**31 - 1
_Positivo = Annotated[int, Field(ge=1, le=_MAX_INTERO, strict=True)]


class _Gruppo(BaseModel):
    """Base dei gruppi: immutabile, chiavi sconosciute respinte."""

    model_config = ConfigDict(extra="forbid", frozen=True)


class Design(_Gruppo):
    """Le colonne della tabella dei campioni che descrivono lo studio."""

    #: La variabile biologica: una sola colonna, senza predefinito.
    variable: _Nome
    #: Altre variabili biologiche, analizzate una alla volta.
    other_variables: tuple[_Nome, ...] = ()
    #: Variabili tecniche (lotto): primi termini dei modelli.
    technical_variables: tuple[_Nome, ...] = ()
    #: Sottoinsieme: colonna -> valori ammessi, applicato prima di tutto.
    subset: dict[_Nome, tuple[_Nome, ...]] = Field(default_factory=dict)

    @model_validator(mode="after")
    def _senza_ripetizioni(self) -> "Design":
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
        vuoti = sorted(c for c, valori in self.subset.items() if not valori)
        if vuoti:
            raise ValueError(f"design.subset senza valori ammessi per: {vuoti}")
        return self


class Run(_Gruppo):
    """Parametri dell'esecuzione."""

    #: Seme di rarefazione e di ogni altra operazione casuale: va dichiarato.
    seed: Annotated[int, Field(ge=0, le=_MAX_INTERO, strict=True)]


class Alpha(_Gruppo):
    """Alfa diversita'."""

    #: Profondita' di rarefazione. Nullo: la minima fra i campioni analizzati
    #: (phyloseq::rarefy_even_depth, sample.size = min(sample_sums(physeq))).
    rarefy_depth: _Positivo | None = None


class Comp(_Gruppo):
    """Composizione tassonomica."""

    rank: _Nome
    top_n: _Positivo


class Beta(_Gruppo):
    """Beta diversita'."""

    distances: Annotated[tuple[Distanza, ...], Field(min_length=1)]
    #: Pseudoconteggio del CLR: solo con la distanza di Aitchison.
    clr_pseudocount: Annotated[float, Field(gt=0, strict=True, allow_inf_nan=False)] | None = None
    #: Radice dell'albero: solo con una delle due UniFrac.
    unifrac_root: Literal["midpoint", "existing"] | None = None

    @model_validator(mode="after")
    def _senza_ripetizioni(self) -> "Beta":
        ripetute = sorted({d for d in self.distances if self.distances.count(d) > 1})
        if ripetute:
            raise ValueError(f"distanze ripetute: {ripetute}")
        return self


class Ord(_Gruppo):
    """Ordinazione."""

    #: Metodi di ordinazione; un elenco vuoto non ne produce.
    methods: tuple[Ordinazione, ...]
    #: Distanze da ordinare, fra quelle di ``beta.distances``.
    distances: tuple[Distanza, ...]
    #: Correzione degli autovalori negativi: solo con la PCoA.
    pcoa_correction: Literal["none", "cailliez", "lingoes"] | None = None
    #: Numero di avvii casuali: solo con la NMDS.
    nmds_trymax: _Positivo | None = None

    @model_validator(mode="after")
    def _coerente(self) -> "Ord":
        for nome in ("methods", "distances"):
            valori = getattr(self, nome)
            ripetuti = sorted({v for v in valori if valori.count(v) > 1})
            if ripetuti:
                raise ValueError(f"ord.{nome} con valori ripetuti: {ripetuti}")
        if bool(self.methods) != bool(self.distances):
            raise ValueError(
                "ord.methods e ord.distances vanno dichiarati insieme: entrambi "
                "con almeno un valore, oppure entrambi vuoti"
            )
        return self


class Stat(_Gruppo):
    """Gruppi e test."""

    #: Dimensione minima di un gruppo per entrare nei test.
    min_group_size: Annotated[int, Field(ge=2, le=_MAX_INTERO, strict=True)]
    #: Livello sotto il quale le dispersioni si dichiarano diverse.
    significance_level: Annotated[float, Field(gt=0, lt=1, strict=True)]
    #: Permutazioni di PERMANOVA e dispersioni (vegan::adonis2, permutations = 999).
    permanova_permutations: _Positivo = 999
    #: Colonna degli strati entro cui permutare nella PERMANOVA.
    permanova_strata: _Nome | None = None


class ConfigEco(_Gruppo):
    """La configurazione di un'analisi."""

    design: Design
    run: Run
    alpha: Alpha = Field(default_factory=Alpha)
    comp: Comp
    beta: Beta
    ord: Ord
    stat: Stat

    @model_validator(mode="after")
    def _distanze_ordinate(self) -> "ConfigEco":
        estranee = [d for d in self.ord.distances if d not in self.beta.distances]
        if estranee:
            raise ValueError(
                f"ord.distances nomina distanze che beta.distances non calcola: {estranee}"
            )
        return self


# --------------------------------------------------------------------------- #
# Dai problemi dello schema ai codici del catalogo                             #
# --------------------------------------------------------------------------- #

#: Parametri il cui valore e' il nome di un metodo: un valore non ammesso e' un
#: metodo non supportato (E-ECO-03), non un errore generico.
_PARAMETRI_DI_METODO: Final = (
    "beta.distances", "beta.unifrac_root", "ord.methods", "ord.distances",
    "ord.pcoa_correction",
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


def _condizionali(config: ConfigEco) -> list[tuple[str, str]]:
    """I rifiuti dei parametri condizionali: mancanti con il metodo richiesto,
    oppure dichiarati senza il metodo che li usa.
    """
    rifiuti: list[tuple[str, str]] = []
    for chiave, (elenco, metodi) in CONDIZIONALI.items():
        gruppo, nome = chiave.split(".")
        gruppo_elenco, nome_elenco = elenco.split(".")
        valore = getattr(getattr(config, gruppo), nome)
        dichiarati = getattr(getattr(config, gruppo_elenco), nome_elenco)
        richiesti = [m for m in metodi if m in dichiarati]
        if richiesti and valore is None:
            rifiuti.append((
                "E-ECO-02",
                f"Parametro mancante: {chiave}, obbligatorio perche' {elenco} "
                f"comprende {richiesti[0]}.",
            ))
        if not richiesti and valore is not None:
            alternative = " o ".join(metodi)
            rifiuti.append((
                "E-ECO-01",
                f"{chiave}: dichiarato, ma {elenco} non comprende {alternative}; "
                f"toglilo, oppure richiedi {alternative}.",
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
        raise ErroreEco([
            _rifiuto(p) for p in problemi
            if not (p["type"] == "too_short" and _chiave(p["loc"]) in respinti)
        ]) from None
    rifiuti = _condizionali(config)
    if rifiuti:
        raise ErroreEco(rifiuti)
    return config


class _SenzaChiaviRipetute(yaml.SafeLoader):
    """Lettore YAML che respinge una chiave dichiarata due volte nella stessa
    mappa: il lettore normale terrebbe l'ultima in silenzio.
    """


def _mappa_senza_ripetizioni(lettore: yaml.SafeLoader, nodo: yaml.MappingNode) -> dict[Any, Any]:
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
    predefiniti compresi, in tipi JSON.
    """
    return config.model_dump(mode="json")


def digest(config: ConfigEco) -> str:
    """SHA-256 della configurazione risolta in forma canonica (chiavi
    ordinate): dipende dai valori, non dai byte del file.
    """
    canonico = json.dumps(risolta(config), sort_keys=True, ensure_ascii=False,
                          separators=(",", ":"))
    return hashlib.sha256(canonico.encode("utf-8")).hexdigest()
