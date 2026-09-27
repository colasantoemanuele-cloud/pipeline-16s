"""La politica dei tentativi ripetuti e le azioni correttive.

**Il retry è un elenco chiuso.** Si ritenta solo un codice che sta in
``retry.whitelist`` e che il catalogo ammette al retry: il criterio è che
l'azione correttiva non modifichi alcuna assunzione metodologica. La
whitelist della configurazione è l'autorità, ma può solo restringere l'elenco
del catalogo; G15 respinge una whitelist che lo allarghi (``E-G15-09``) prima
che giri qualunque fase.

**``retry.max_attempts`` conta i tentativi totali, compreso il primo.** Con il
valore 2 un codice ammesso viene ritentato una volta sola; con 1 non viene
ritentato mai, come con ``retry.enabled`` falso.

**Ritentare identico non serve.** Ogni tentativo successivo applica l'azione
correttiva che il catalogo dichiara per quel codice: una fase la dichiara con
:class:`Aggiustamento`, per ciascun codice ripetibile che può sollevare. Se
per un codice la fase non dichiara un aggiustamento, o l'aggiustamento non
può più cambiare il valore, non si ritenta. L'unica eccezione è dichiarata,
non presunta: per un errore transitorio, come un errore di lettura, la fase
puo' dichiarare :func:`senza_modifiche`, e il nuovo tentativo identico viene
registrato come ogni altro aggiustamento.

**L'aggiustamento non tocca la configurazione dell'utente.** Produce una
configurazione nuova, in memoria, con il solo parametro cambiato e
rivalidata dallo schema. La fase la usa per il calcolo; il manifesto registra
il valore dichiarato e quello usato, con il codice che ha causato il cambio,
mentre la validità della fase alla ripresa si giudica sulla configurazione
dichiarata. È coerente solo se il parametro aggiustato non incide sui
risultati: è l'assunzione su cui si regge la whitelist, e per
``run.batch_size`` non è ancora verificata (vedi :data:`PARAMETRI_AGGIUSTABILI`).
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass
from enum import StrEnum
from typing import Any, Final

from amplicon16s.config.schema import Config, valida
from amplicon16s.errors.catalog import Categoria, voce

__all__ = [
    "PARAMETRI_AGGIUSTABILI",
    "Aggiustamento",
    "Decisione",
    "Motivo",
    "PoliticaRetry",
    "RITENTARE_INUTILE",
    "applica",
    "dimezza",
    "raddoppia",
    "senza_modifiche",
    "spiega_arresto",
    "valore",
]

#: I parametri che un'azione correttiva può cambiare. Un aggiustamento su un
#: parametro fuori elenco è un difetto della fase.
#:
#: Sono quelli che i messaggi del catalogo indicano per i codici ripetibili:
#: ``run.batch_size`` per E-S4-02 ed E-S5-01, ``err.nbases`` per E-S3-01.
#: E-S2-03, un errore di lettura, si ritenta invece senza modifiche
#: (:func:`senza_modifiche`): se era transitorio basta rileggere. **Tensione aperta**: entrambi restano nell'impronta dei risultati,
#: perché non è dimostrato che non incidano sui risultati. Con
#: ``dada.pool: pseudo``, il predefinito, i priori della seconda passata di
#: dada2 si costruiscono dai campioni elaborati insieme: se l'inferenza
#: procede per lotti di ``run.batch_size`` campioni, ridurre il lotto può
#: cambiare le varianti inferite, e allora ridurlo in un tentativo automatico
#: violerebbe il criterio stesso della whitelist. ``err.nbases`` cambia per
#: definizione la stima del modello d'errore. La decisione va presa quando S2-S5
#: esisteranno, misurando l'effetto sul dataset di riferimento; fino ad allora
#: ogni aggiustamento resta registrato nel manifesto con il valore usato.
PARAMETRI_AGGIUSTABILI: Final[tuple[str, ...]] = ("run.batch_size", "err.nbases")


def valore(config: Config, parametro: str) -> Any:
    """Il valore di un parametro ``gruppo.chiave``."""
    gruppo, chiave = parametro.split(".")
    return getattr(getattr(config, gruppo), chiave)


def applica(config: Config, parametro: str, nuovo: Any) -> Config:
    """Una configurazione nuova con un solo parametro cambiato, rivalidata.

    La configurazione di partenza, e il file da cui viene, restano intatti.
    """
    gruppo, chiave = parametro.split(".")
    dati = config.model_dump(mode="python")
    dati[gruppo][chiave] = nuovo
    return valida(dati, origine=f"aggiustamento di {parametro}")


@dataclass(frozen=True)
class Aggiustamento:
    """L'azione correttiva di un codice ripetibile: come cambiare un parametro.

    Senza parametro, l'aggiustamento dichiara un nuovo tentativo **senza
    modifiche** (:func:`senza_modifiche`): e' l'azione giusta solo per un
    errore transitorio, e va dichiarata esplicitamente, codice per codice, con
    il suo motivo. In assenza di una dichiarazione la regola resta quella
    generale: non si ritenta identico.
    """

    parametro: str | None
    regola: Callable[[Any], Any] | None
    descrizione: str

    def __post_init__(self) -> None:
        if self.parametro is None:
            if self.regola is not None:
                raise ValueError("un nuovo tentativo senza modifiche non ha una regola")
            return
        if self.parametro not in PARAMETRI_AGGIUSTABILI:
            raise ValueError(
                f"{self.parametro} non e' un parametro aggiustabile: ammessi "
                f"{', '.join(PARAMETRI_AGGIUSTABILI)}"
            )

    @property
    def invariato(self) -> bool:
        """Se il nuovo tentativo avviene senza cambiare alcun parametro."""
        return self.parametro is None

    def prossimo(self, config: Config) -> Any:
        """Il valore da usare nel tentativo successivo."""
        if self.parametro is None or self.regola is None:
            raise ValueError("un nuovo tentativo senza modifiche non cambia valori")
        return self.regola(valore(config, self.parametro))


def senza_modifiche(motivo: str) -> Aggiustamento:
    """Nuovo tentativo identico, per un errore dichiaratamente transitorio.

    ``motivo`` spiega perche' ritentare identico puo' riuscire, e finisce nel
    manifesto e nel log accanto al codice. Un errore che non e' transitorio,
    come un file davvero corrotto, fallira' di nuovo e si fermera' dopo i
    tentativi ammessi.
    """
    if not motivo.strip():
        raise ValueError("un nuovo tentativo senza modifiche va motivato")
    return Aggiustamento(None, None, f"nuovo tentativo senza modifiche: {motivo}")


def dimezza(parametro: str, minimo: int = 1) -> Aggiustamento:
    """Dimezza un parametro intero, senza scendere sotto ``minimo``."""
    return Aggiustamento(
        parametro, lambda v: max(minimo, v // 2), f"{parametro} dimezzato"
    )


def raddoppia(parametro: str) -> Aggiustamento:
    """Raddoppia un parametro."""
    return Aggiustamento(parametro, lambda v: v * 2, f"{parametro} raddoppiato")


class Motivo(StrEnum):
    """Perché si ritenta, o perché no."""

    RITENTA = "ritenta"
    REVISIONE_UMANA = "revisione_umana"
    RETRY_DISATTIVATO = "retry_disattivato"
    FUORI_WHITELIST = "fuori_whitelist"
    TENTATIVI_ESAURITI = "tentativi_esauriti"
    NESSUNA_AZIONE = "nessuna_azione_correttiva"
    AZIONE_ESAURITA = "azione_correttiva_esaurita"
    #: La fase che ha sollevato l'errore sa che l'azione correttiva non
    #: cambierebbe l'esito, e lo dichiara nel contesto dell'errore.
    AZIONE_INUTILE = "azione_correttiva_inutile"


#: Chiave del contesto di un errore con cui la fase dichiara che l'azione
#: correttiva non cambierebbe l'esito, e perche'. Solo la fase lo sa: per
#: esempio err.nbases raddoppiato non serve se la stima usava gia' tutte le
#: basi disponibili.
RITENTARE_INUTILE: Final = "ritentare_inutile"


@dataclass(frozen=True)
class Decisione:
    motivo: Motivo
    #: Per :attr:`Motivo.RITENTA`, il valore da usare.
    nuovo_valore: Any = None

    @property
    def ritenta(self) -> bool:
        return self.motivo is Motivo.RITENTA


@dataclass(frozen=True)
class PoliticaRetry:
    """Quando si ritenta: ``retry.enabled``, ``retry.max_attempts``, ``retry.whitelist``."""

    abilitata: bool
    #: Tentativi totali, compreso il primo.
    tentativi_massimi: int
    whitelist: frozenset[str]

    @classmethod
    def da_config(cls, config: Config) -> PoliticaRetry:
        return cls(
            abilitata=config.retry.enabled,
            tentativi_massimi=config.retry.max_attempts,
            whitelist=frozenset(config.retry.whitelist),
        )

    def decidi(
        self,
        codice: str,
        tentativo: int,
        aggiustamento: Aggiustamento | None,
        config: Config,
        inutile: str | None = None,
    ) -> Decisione:
        """Se, dopo il fallimento del tentativo ``tentativo`` (da 1), si ritenta."""
        if not voce(codice).ammette_retry:
            return Decisione(Motivo.REVISIONE_UMANA)
        if not self.abilitata:
            return Decisione(Motivo.RETRY_DISATTIVATO)
        if codice not in self.whitelist:
            return Decisione(Motivo.FUORI_WHITELIST)
        if tentativo >= self.tentativi_massimi:
            return Decisione(Motivo.TENTATIVI_ESAURITI)
        if aggiustamento is None:
            return Decisione(Motivo.NESSUNA_AZIONE)
        if inutile:
            return Decisione(Motivo.AZIONE_INUTILE)
        if aggiustamento.invariato:
            return Decisione(Motivo.RITENTA)
        assert aggiustamento.parametro is not None
        attuale = valore(config, aggiustamento.parametro)
        nuovo = aggiustamento.prossimo(config)
        if nuovo == attuale:
            return Decisione(Motivo.AZIONE_ESAURITA)
        return Decisione(Motivo.RITENTA, nuovo)


def spiega_arresto(
    codice: str, motivo: Motivo, tentativi: int, massimi: int, inutile: str | None = None
) -> str:
    """Perché l'esecuzione si è fermata, secondo la categoria del codice."""
    categoria = voce(codice).categoria
    if motivo is Motivo.AZIONE_INUTILE:
        return (
            "Il codice ammetterebbe il retry, ma l'azione correttiva non cambierebbe "
            f"l'esito: {inutile}. Nessun nuovo tentativo."
        )
    if motivo is Motivo.REVISIONE_UMANA:
        return "Il codice richiede la revisione umana: nessun tentativo automatico."
    if motivo is Motivo.RETRY_DISATTIVATO:
        return "Il codice ammetterebbe il retry, ma retry.enabled e' false."
    if motivo is Motivo.FUORI_WHITELIST:
        return "Il codice ammetterebbe il retry, ma non e' in retry.whitelist."
    if motivo is Motivo.NESSUNA_AZIONE:
        return (
            "La fase non dichiara un'azione correttiva per questo codice, e "
            "ritentare identico darebbe lo stesso esito."
        )
    if motivo is Motivo.AZIONE_ESAURITA:
        return (
            "L'azione correttiva non puo' piu' cambiare il parametro: e' gia' al "
            "suo limite."
        )
    if categoria is Categoria.RETRY_POI_REVISIONE:
        return (
            f"Tentativi automatici esauriti ({tentativi} su {massimi}). Per questo "
            "codice, esauriti i tentativi, serve la revisione umana prima di "
            "riprendere: la causa va capita, non ritentata."
        )
    return (
        f"Tentativi automatici esauriti ({tentativi} su {massimi}). Il codice "
        "ammette il retry: rimossa la causa indicata qui sotto, la ripresa "
        "ricomincia i tentativi da capo."
    )
