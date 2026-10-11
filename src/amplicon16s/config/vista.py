"""La configurazione come la vede una fase: solo i parametri che dichiara.

Ogni fase dichiara i parametri da cui dipende, e la validita' di un suo
risultato si giudica su quelli: cambiare un parametro che non dichiara non la
rende da rifare. Una dipendenza dimenticata produrrebbe quindi un risultato
obsoleto senza alcun errore. Per questo la dichiarazione e' vincolante per
costruzione: il codice di una fase vede la configurazione attraverso una
:class:`VistaConfig`, e leggere un parametro non dichiarato solleva
:class:`ParametroNonDichiarato`. L'errore si presenta al primo test che
esercita la fase, non in un risultato sbagliato mesi dopo.

Oltre ai parametri dichiarati la vista lascia leggere quelli che per
costruzione non incidono sui risultati (``PARAMETRI_SENZA_EFFETTO``): quanti
processori usare, quanti campioni elaborare insieme, dove scrivere, se
conservare le letture filtrate, se e quante volte ritentare. Non entrano
nell'impronta di nessuna fase, e dichiararli sarebbe un errore.

Un parametro condizionale omesso perche' non pertinente non ha un valore: la
vista solleva :class:`CondizionaleNonPertinente` se una fase lo legge senza
aver prima verificato la condizione che lo rende pertinente. Fanno eccezione i
condizionali che una fase riporta in un riepilogo
(``defaults.CONDIZIONALI_RIPORTATI``): omessi si leggono nulli, perche' la fase
li scrive tali, e sono gli script di calcolo a fermarsi se un nullo servisse.

``ParametroNonDichiarato`` non deriva da ``AttributeError`` di proposito: un
``getattr(vista, nome, predefinito)`` o un ``hasattr`` lo inghiottirebbero, e
la dipendenza dimenticata tornerebbe silenziosa.
"""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import fields
from typing import Any

from amplicon16s.config.defaults import CONDIZIONALI_SENZA_VALORE, testo_condizione
from amplicon16s.config.resolve import PARAMETRI_SENZA_EFFETTO, ConfigRisolta, Derivati
from amplicon16s.config.schema import Config

__all__ = [
    "CondizionaleNonPertinente", "ParametroNonDichiarato", "VistaConfig", "risolta_ristretta",
]


class ParametroNonDichiarato(LookupError):
    """Una fase ha letto un parametro che non ha dichiarato.

    Razionale sistemistico: eredita da LookupError anziche' da AttributeError
    affinche' costrutti come getattr(vista, nome, default) o hasattr(vista, nome)
    non catturino l'eccezione restituendo silenziosamente un valore predefinito.
    """


class CondizionaleNonPertinente(LookupError):
    """Una fase ha letto un condizionale omesso senza guardarne la condizione.

    Un condizionale non pertinente non ha un valore: restituire un nullo
    lascerebbe passare in silenzio una lettura che la fase non dovrebbe fare.
    Come ``ParametroNonDichiarato``, non deriva da ``AttributeError``.
    """


def _ammessi(parametri: Iterable[str]) -> tuple[frozenset[str], frozenset[str]]:
    """Gruppi dichiarati per intero, e chiavi ``gruppo.parametro`` ammesse una a una."""
    gruppi = {p for p in parametri if "." not in p}
    # La modalita' di lettura si legge sempre: decide quali parametri la fase
    # dichiara (PipelineStep.parametri_per), e una fase che non la dichiara e'
    # una fase che in modalita' single-end non ne dipende.
    chiavi = {p for p in parametri if "." in p} | set(PARAMETRI_SENZA_EFFETTO) | {"io.layout"}
    return frozenset(gruppi), frozenset(chiavi)


class _VistaGruppo:
    """Un gruppo di parametri visto da una fase: le sole chiavi dichiarate, o
    tutte (``chiavi`` nullo) se la fase ha dichiarato il gruppo per intero."""

    def __init__(self, nome: str, gruppo: Any, chiavi: frozenset[str] | None) -> None:
        object.__setattr__(self, "_nome", nome)
        object.__setattr__(self, "_gruppo", gruppo)
        object.__setattr__(self, "_chiavi", chiavi)

    def __getattr__(self, chiave: str) -> Any:
        if chiave.startswith("_"):
            raise AttributeError(chiave)
        if self._chiavi is not None and chiave not in self._chiavi:
            raise ParametroNonDichiarato(
                f"{self._nome}.{chiave} non e' fra i parametri dichiarati dalla fase"
            )
        valore = getattr(self._gruppo, chiave)
        nome = f"{self._nome}.{chiave}"
        if valore is None and nome in CONDIZIONALI_SENZA_VALORE:
            raise CondizionaleNonPertinente(
                f"{nome} non e' dichiarato perche' non pertinente: si legge solo dopo "
                f"aver verificato che {testo_condizione(nome)}"
            )
        return valore

    def __setattr__(self, nome: str, valore: Any) -> None:
        raise AttributeError("la configurazione di una fase non si modifica")


class VistaConfig:
    """La configurazione ristretta ai parametri dichiarati da una fase."""

    def __init__(self, config: Config, parametri: Iterable[str]) -> None:
        gruppi, chiavi = _ammessi(parametri)
        object.__setattr__(self, "_config", config)
        object.__setattr__(self, "_gruppi", gruppi)
        object.__setattr__(self, "_chiavi", chiavi)

    def __getattr__(self, nome: str) -> Any:
        if nome.startswith("_"):
            raise AttributeError(nome)
        if nome == "appaiata":
            return self._config.appaiata
        if nome not in Config.model_fields:
            # Anche model_dump e simili: l'intera configurazione non si legge.
            raise ParametroNonDichiarato(
                f"{nome}: una fase legge solo i parametri che dichiara, non "
                "l'intera configurazione"
            )
        gruppo = getattr(self._config, nome)
        if nome in self._gruppi:
            # Dichiarato per intero: ogni chiave si legge, sotto la sola guardia
            # dei condizionali senza valore.
            return _VistaGruppo(nome, gruppo, None)
        # Anche un gruppo non dichiarato si attraversa: l'errore arriva sulla
        # chiave letta, e ne porta il nome completo.
        chiavi = frozenset(c.split(".", 1)[1] for c in self._chiavi if c.split(".", 1)[0] == nome)
        return _VistaGruppo(nome, gruppo, chiavi)

    def __setattr__(self, nome: str, valore: Any) -> None:
        raise AttributeError("la configurazione di una fase non si modifica")


#: A quale parametro corrisponde ciascun derivato: si legge se quello si legge.
_DERIVATI = {
    "asv_len_min": "asv.len_min",
    "asv_len_max": "asv.len_max",
}


class _VistaDerivati:
    """I parametri derivati, sotto lo stesso vincolo dei dichiarati."""

    def __init__(self, derivati: Derivati, parametri: Iterable[str]) -> None:
        gruppi, chiavi = _ammessi(parametri)
        object.__setattr__(self, "_derivati", derivati)
        object.__setattr__(self, "_gruppi", gruppi)
        object.__setattr__(self, "_chiavi", chiavi)

    def __getattr__(self, nome: str) -> Any:
        if nome.startswith("_"):
            raise AttributeError(nome)
        chiave = _DERIVATI.get(nome)
        if chiave is None:
            raise ParametroNonDichiarato(f"{nome}: derivato sconosciuto o accesso non ammesso")
        if chiave.split(".")[0] not in self._gruppi and chiave not in self._chiavi:
            raise ParametroNonDichiarato(
                f"{chiave} non e' fra i parametri dichiarati dalla fase"
            )
        return getattr(self._derivati, nome)


def risolta_ristretta(risolta: ConfigRisolta, parametri: Iterable[str]) -> ConfigRisolta:
    """La configurazione risolta come la vede una fase.

    Configurazione e derivati passano per la vista; il digest e la mappa
    completa, che leggerebbero tutto, sollevano ``ParametroNonDichiarato``.
    """
    # Razionale sistemistico: la sostituzione di config e derivati con le
    # rispettive viste ristrette mantiene intatti i valori gia' risolti (inclusi
    # gli aggiustamenti applicati dal retry in memoria), ma impedisce alla fase
    # di calcolare il digest globale o serializzare l'intera configurazione.
    parametri = tuple(parametri)
    assert {f.name for f in fields(Derivati)} == set(_DERIVATI)
    return ConfigRisolta(
        config=VistaConfig(risolta.config, parametri),  # type: ignore[arg-type]
        derivati=_VistaDerivati(risolta.derivati, parametri),  # type: ignore[arg-type]
    )
