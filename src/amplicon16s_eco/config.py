"""Configurazione delle analisi ecologiche: schema, caricamento e digest.

Il file e' proprio del pacchetto, distinto da quello della pipeline, e ha
gruppi ``design``, ``run``, ``alpha``, ``comp``, ``beta`` e ``stat``. Le regole
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
    "DISTANZE",
    "FACOLTATIVI",
    "OBBLIGATORI",
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

#: Parametri senza predefinito: dipendono dallo studio o sono scelte di analisi.
OBBLIGATORI: Final[tuple[str, ...]] = (
    "design.variable",
    "run.seed",
    "comp.rank",
    "comp.top_n",
    "beta.distances",
    "stat.min_group_size",
)

#: Parametri obbligatori quando il metodo indicato e' richiesto, respinti
#: altrimenti: chiave -> (parametro che elenca i metodi, metodo).
CONDIZIONALI: Final[dict[str, tuple[str, str]]] = {
    "beta.clr_pseudocount": ("beta.distances", "aitchison"),
}

#: Predefiniti che sono il valore del pacchetto R che realizza il metodo:
#: chiave -> (valore, fonte).
STANDARD_DEL_METODO: Final[dict[str, tuple[Any, str]]] = {
    "alpha.rarefy_depth": (
        None,
        "phyloseq::rarefy_even_depth, sample.size = min(sample_sums(physeq)): "
        "nullo significa la minima profondita' fra i campioni analizzati",
    ),
}

#: Parametri facoltativi il cui predefinito e' l'assenza.
FACOLTATIVI: Final[tuple[str, ...]] = (
    "design.other_variables",
    "design.technical_variables",
    "design.subset",
)

_Nome = Annotated[str, Field(min_length=1, strict=True)]


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
    seed: Annotated[int, Field(ge=0, le=2**31 - 1, strict=True)]


class Alpha(_Gruppo):
    """Alfa diversita'."""

    #: Profondita' di rarefazione. Nullo: la minima fra i campioni analizzati
    #: (phyloseq::rarefy_even_depth, sample.size = min(sample_sums(physeq))).
    rarefy_depth: Annotated[int, Field(ge=1, strict=True)] | None = None


class Comp(_Gruppo):
    """Composizione tassonomica."""

    rank: _Nome
    top_n: Annotated[int, Field(ge=1, strict=True)]


class Beta(_Gruppo):
    """Beta diversita'."""

    distances: Annotated[tuple[Distanza, ...], Field(min_length=1)]
    #: Pseudoconteggio del CLR: solo con la distanza di Aitchison.
    clr_pseudocount: Annotated[float, Field(gt=0)] | None = None

    @model_validator(mode="after")
    def _senza_ripetizioni(self) -> "Beta":
        ripetute = sorted({d for d in self.distances if self.distances.count(d) > 1})
        if ripetute:
            raise ValueError(f"distanze ripetute: {ripetute}")
        return self


class Stat(_Gruppo):
    """Gruppi e test."""

    #: Dimensione minima di un gruppo per entrare nei test.
    min_group_size: Annotated[int, Field(ge=2, strict=True)]


class ConfigEco(_Gruppo):
    """La configurazione di un'analisi."""

    design: Design
    run: Run
    alpha: Alpha = Field(default_factory=Alpha)
    comp: Comp
    beta: Beta
    stat: Stat


# --------------------------------------------------------------------------- #
# Dai problemi dello schema ai codici del catalogo                             #
# --------------------------------------------------------------------------- #

#: Parametri il cui valore e' il nome di un metodo: un valore non ammesso e' un
#: metodo non supportato (E-ECO-03), non un errore generico.
_PARAMETRI_DI_METODO: Final = ("beta.distances",)

_MESSAGGI: Final[dict[str, str]] = {
    "extra_forbidden": "parametro sconosciuto (refuso o parametro non previsto)",
    "int_type": "deve essere un numero intero",
    "float_type": "deve essere un numero",
    "string_type": "deve essere una stringa (fra virgolette, se e' un numero)",
    "tuple_type": "deve essere un elenco",
    "dict_type": "deve essere una mappa",
    "model_type": "deve essere una mappa di parametri",
    "greater_than": "deve essere maggiore di {gt}",
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
    for chiave, (elenco, metodo) in CONDIZIONALI.items():
        gruppo, nome = chiave.split(".")
        gruppo_elenco, nome_elenco = elenco.split(".")
        valore = getattr(getattr(config, gruppo), nome)
        richiesto = metodo in getattr(getattr(config, gruppo_elenco), nome_elenco)
        if richiesto and valore is None:
            rifiuti.append((
                "E-ECO-02",
                f"Parametro mancante: {chiave}, obbligatorio perche' {elenco} "
                f"comprende {metodo}.",
            ))
        if not richiesto and valore is not None:
            rifiuti.append((
                "E-ECO-01",
                f"{chiave}: dichiarato, ma {elenco} non comprende {metodo}; "
                f"toglilo, oppure richiedi {metodo}.",
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


def carica(percorso: Path | str) -> ConfigEco:
    """Legge e valida il file di configurazione dell'analisi."""
    percorso = Path(percorso)
    try:
        dati = yaml.safe_load(percorso.read_text(encoding="utf-8"))
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
