"""La classe base delle fasi: lo scheletro comune a S0-S14.

Ogni fase segue gli stessi quattro momenti, e :meth:`PipelineStep.esegui` li
fissa una volta sola:

1. **verifica dei prerequisiti**: le fasi da cui dipende sono concluse;
2. **esecuzione**: il calcolo proprio della fase, :meth:`PipelineStep.calcola`,
   l'unica parte che una fase concreta deve scrivere;
3. **validazione degli artefatti**: ciò che la fase dichiara di aver prodotto
   esiste, sta nella sua cartella ed è integro;
4. **registrazione dell'esito**: il manifesto della fase, con gli artefatti e
   ciò su cui è stata calcolata.

**Su che cosa è stata calcolata una fase.** Artefatti integri non bastano a
dire che una fase è conclusa: se una fase a monte viene ricalcolata, o cambia
un parametro, i file a valle restano integri ma sono stati calcolati su
ingressi che non esistono piu'. Il manifesto registra quindi, in
``calcolata_su``, tre cose:

* ``configurazione``: l'impronta dei parametri da cui la fase dipende — per
  difetto l'intera configurazione risolta, meno i parametri che non incidono
  sui risultati, come il numero di thread o la cartella di output;
* ``a_monte``: l'impronta del manifesto di ogni fase di cui consuma gli
  artefatti. Ogni esecuzione produce un'impronta nuova, quindi ricalcolare una
  fase invalida tutte quelle che l'avevano registrata, e a cascata le loro;
* ``dati_esterni``: per le fasi che leggono dati fuori dall'albero di output,
  un'impronta di quei dati.

Una fase è conclusa se il suo manifesto esiste, i suoi artefatti sono integri
e ``calcolata_su`` coincide con ciò che si otterrebbe adesso. Il confronto lo
fa :class:`~amplicon16s.runner.project.ProjectRun`, con lo stesso metodo
:meth:`PipelineStep.calcolata_su` che la fase usa per registrarlo.
"""

from __future__ import annotations

import hashlib
import json
import logging
import time
from abc import ABC, abstractmethod
from collections.abc import Mapping
from dataclasses import dataclass, field
from enum import StrEnum
from typing import Any, ClassVar

from amplicon16s.config.resolve import PARAMETRI_SENZA_EFFETTO, ConfigRisolta
from amplicon16s.config.schema import Config
from amplicon16s.errors.exceptions import errore
from amplicon16s.io_layer.artifacts import AlberoOutput, Artefatto, Fase
from amplicon16s.metadata.models import Inventario
from amplicon16s.runner.graph import GRAFO, Grafo, Nodo, Passo

__all__ = [
    "Esito",
    "PipelineStep",
    "Produzione",
    "StepContext",
    "StepResult",
    "impronta_parametri",
]


def _impronta(valore: Any) -> str:
    canonico = json.dumps(
        valore, sort_keys=True, ensure_ascii=False, separators=(",", ":"), default=str
    )
    return "sha256:" + hashlib.sha256(canonico.encode("utf-8")).hexdigest()


def impronta_parametri(
    risolta: ConfigRisolta, parametri: tuple[str, ...] | None
) -> str:
    """Impronta dei parametri da cui una fase dipende.

    ``None`` indica l'intera configurazione, meno i parametri che non
    incidono sui risultati (``PARAMETRI_SENZA_EFFETTO`` in
    :mod:`amplicon16s.config.resolve`). Altrimenti ``parametri`` elenca gruppi
    (``"tax"``) o chiavi (``"filter.truncLen"``), e anche da un gruppo quei
    parametri restano fuori. Un nome inesistente è un errore, perché una fase
    che dichiara male i propri parametri non vedrebbe cambiare quelli che usa
    davvero; e lo è nominarne uno escluso, che non entrerebbe nell'impronta.
    """
    if parametri is None:
        return risolta.impronta_risultati

    mappa = risolta.come_mappa()
    for escluso in PARAMETRI_SENZA_EFFETTO:
        gruppo, chiave = escluso.split(".")
        del mappa[gruppo][chiave]

    scelti: dict[str, Any] = {}
    for nome in parametri:
        if nome in PARAMETRI_SENZA_EFFETTO:
            raise ValueError(f"{nome} non incide sui risultati ed e' escluso dall'impronta")
        gruppo, _, chiave = nome.partition(".")
        if gruppo not in mappa or (chiave and chiave not in mappa[gruppo]):
            raise KeyError(f"parametro inesistente: {nome}")
        scelti[nome] = mappa[gruppo][chiave] if chiave else mappa[gruppo]
    return _impronta(scelti)


@dataclass(frozen=True)
class StepContext:
    """Ciò che una fase riceve per essere eseguita."""

    risolta: ConfigRisolta
    albero: AlberoOutput
    logger: logging.Logger
    #: ``None`` prima che S0 lo abbia costruito.
    inventario: Inventario | None
    #: Impronta del manifesto di ogni dipendenza attiva, ``None`` se quella
    #: dipendenza non risulta conclusa.
    a_monte: Mapping[Passo, str | None] = field(default_factory=dict)

    @property
    def config(self) -> Config:
        return self.risolta.config


class Esito(StrEnum):
    """Come si è conclusa l'esecuzione di una fase."""

    #: La fase ha prodotto i suoi artefatti ed è registrata come conclusa.
    COMPLETATA = "completata"
    #: La fase ha girato ma il suo esito non la conclude: gli artefatti
    #: prodotti restano per diagnosi, ma nessun manifesto li dichiara validi.
    NON_SUPERATA = "non_superata"


@dataclass(frozen=True)
class Produzione:
    """Ciò che :meth:`PipelineStep.calcola` restituisce."""

    artefatti: tuple[Artefatto, ...]
    metriche: Mapping[str, Any] = field(default_factory=dict)
    superata: bool = True
    #: Risultato proprio della fase, per chi la invoca direttamente.
    dettaglio: Any = None


@dataclass(frozen=True)
class StepResult:
    """Esito di una fase: artefatti con i loro checksum, metriche, esito."""

    passo: Passo
    esito: Esito
    artefatti: tuple[Artefatto, ...]
    metriche: Mapping[str, Any]
    secondi: float
    #: Impronta del manifesto scritto; ``None`` se la fase non è conclusa.
    impronta: str | None = None
    dettaglio: Any = None

    @property
    def completata(self) -> bool:
        return self.esito is Esito.COMPLETATA


class PipelineStep(ABC):
    """Una fase della pipeline.

    Una fase concreta dichiara quale ``passo`` realizza, e se serve da quali
    ``parametri`` dipende, e implementa :meth:`calcola`. Il resto — ordine,
    prerequisiti, validazione, registrazione — è comune.
    """

    #: La fase del grafo che questa classe realizza.
    passo: ClassVar[Passo]
    #: Parametri da cui il risultato dipende; ``None`` è l'intera
    #: configurazione. Restringerli evita di ricalcolare la fase quando cambia
    #: un parametro che non usa, ma un parametro dimenticato qui renderebbe
    #: invisibile un cambiamento che la riguarda: nel dubbio, ``None``.
    parametri: ClassVar[tuple[str, ...] | None] = None
    #: Il grafo di riferimento.
    grafo: ClassVar[Grafo] = GRAFO

    @property
    def nodo(self) -> Nodo:
        return self.grafo.nodo(self.passo)

    @property
    def cartella(self) -> Fase:
        return self.nodo.cartella

    # ----------------------------------------------------------------- #
    # Cio' che una fase concreta scrive                                  #
    # ----------------------------------------------------------------- #

    @abstractmethod
    def calcola(self, contesto: StepContext) -> Produzione:
        """Il calcolo proprio della fase."""

    def impronta_dati_esterni(self, config: Config) -> str | None:
        """Impronta dei dati letti fuori dall'albero di output, se ce ne sono.

        Le fasi che leggono solo artefatti di altre fasi non ne hanno: i loro
        ingressi sono gia' identificati dalle impronte a monte.
        """
        return None

    # ----------------------------------------------------------------- #
    # Lo scheletro comune                                                #
    # ----------------------------------------------------------------- #

    def calcolata_su(
        self, risolta: ConfigRisolta, a_monte: Mapping[Passo, str | None]
    ) -> dict[str, Any]:
        """Su che cosa la fase viene calcolata, nella forma del manifesto."""
        return {
            "configurazione": {
                "parametri": list(self.parametri) if self.parametri is not None else None,
                "impronta": impronta_parametri(risolta, self.parametri),
            },
            "a_monte": {str(p): a_monte[p] for p in sorted(a_monte, key=_ordine)},
            "dati_esterni": self.impronta_dati_esterni(risolta.config),
        }

    def verifica_prerequisiti(self, contesto: StepContext) -> None:
        """Le dipendenze attive devono essere concluse.

        Eseguire una fase prima di una da cui dipende produrrebbe un risultato
        calcolato su artefatti assenti o superati. Se la dipendenza mancante
        e' una precedenza obbligatoria del grafo, l'errore porta il codice di
        quella precedenza (``E-S13-01`` per S12 prima di S13); altrimenti e'
        ``E-GRAFO-01``.
        """
        attese = self.grafo.dipendenze_attive(self.passo, contesto.config)
        mancanti = [p for p in attese if contesto.a_monte.get(p) is None]
        if not mancanti:
            return

        codice = "E-GRAFO-01"
        for prima, dopo, codice_precedenza in self.grafo.precedenze:
            if dopo is self.passo and prima in mancanti:
                codice = codice_precedenza
                break
        raise errore(
            codice,
            f"{self.passo} richiede fasi non concluse: {', '.join(mancanti)}",
            passo=str(self.passo),
            mancanti=[str(p) for p in mancanti],
        )

    def valida_artefatti(
        self, contesto: StepContext, artefatti: tuple[Artefatto, ...]
    ) -> None:
        """Gli artefatti dichiarati stanno nella cartella della fase e sono integri.

        Un artefatto assente o alterato subito dopo il calcolo è un difetto
        della fase, non una condizione dei dati. Una fase concreta può
        aggiungere le proprie verifiche chiamando questa per prima.
        """
        for artefatto in artefatti:
            if artefatto.fase is not self.cartella:
                raise RuntimeError(
                    f"{self.passo} ha dichiarato {artefatto.nome} in "
                    f"{artefatto.fase.value}, ma scrive in {self.cartella.value}"
                )
            if not artefatto.integro:
                raise RuntimeError(
                    f"{self.passo}: l'artefatto {artefatto.nome} manca o non "
                    "corrisponde al checksum registrato"
                )

    def esegui(self, contesto: StepContext) -> StepResult:
        """Prerequisiti, calcolo, validazione, registrazione."""
        log = contesto.logger
        self.verifica_prerequisiti(contesto)

        # Si fissa prima del calcolo su che cosa la fase viene calcolata: e'
        # cio' che le e' stato dato, anche se durante il calcolo cambiasse.
        calcolata_su = self.calcolata_su(contesto.risolta, dict(contesto.a_monte))

        # Da qui la fase non risulta piu' conclusa: se il calcolo si
        # interrompe, la ripresa la rifa' invece di fidarsi del manifesto
        # precedente.
        contesto.albero.rimuovi_manifesto_passo(self.passo, self.cartella)

        inizio = time.perf_counter()
        produzione = self.calcola(contesto)
        secondi = round(time.perf_counter() - inizio, 3)

        self.valida_artefatti(contesto, produzione.artefatti)

        impronta = None
        esito = Esito.NON_SUPERATA
        if produzione.superata:
            manifesto = contesto.albero.concludi_passo(
                self.passo,
                self.cartella,
                produzione.artefatti,
                calcolata_su,
                produzione.metriche,
            )
            impronta = manifesto.impronta
            esito = Esito.COMPLETATA

        log.log(
            logging.INFO if produzione.superata else logging.WARNING,
            f"{self.passo} {esito.value}",
            extra={
                "passo": str(self.passo),
                "esito": esito.value,
                "secondi": secondi,
                "artefatti": [a.nome for a in produzione.artefatti],
                "metriche": dict(produzione.metriche),
            },
        )
        return StepResult(
            passo=self.passo,
            esito=esito,
            artefatti=produzione.artefatti,
            metriche=dict(produzione.metriche),
            secondi=secondi,
            impronta=impronta,
            dettaglio=produzione.dettaglio,
        )


def _ordine(passo: Passo) -> int:
    return list(Passo).index(passo)
