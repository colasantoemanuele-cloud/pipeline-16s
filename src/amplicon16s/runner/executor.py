"""L'esecutore: percorre il grafo, esegue le fasi da fare, ritenta o si ferma.

A ogni avvio, prima di qualunque fase:

1. **coerenza della configurazione** (G15), compresa quella fra
   ``retry.whitelist`` e il catalogo: un errore qui è di configurazione e
   nessuna fase parte;
2. **risorse della macchina** (G14): processori e spazio su disco. G14 fa
   parte di S0, ma verifica la macchina in uso e non i dati, e ``run.threads``
   e ``io.out_root`` non entrano nell'impronta dei risultati: una ripresa con
   S0 già conclusa non ripasserebbe da S0, e una macchina con meno processori
   o un volume piu' piccolo non verrebbero controllati. Qui G14 gira da solo,
   senza rieseguire S0.

Superati i controlli, **registra la configurazione** in ``00_config`` senza
mai sovrascrivere (:func:`~amplicon16s.config.resolve.registra_risolta`): la
prima esecuzione scrive ``resolved.yaml``; una ripresa con lo stesso digest
non scrive nulla; una ripresa con un digest diverso — anche solo per
``run.threads`` o un parametro di retry, fuori dall'impronta dei risultati ma
non dal digest — scrive la nuova versione accanto alle precedenti. Il log
dice a ogni avvio con quale versione si esegue.

Poi esegue, una alla volta, le fasi che la :class:`Valutazione` dà da
eseguire. L'esito di una fase fallita dipende dalla categoria del suo codice:

* **revisione umana**: ci si ferma subito;
* **retry automatico** e **retry poi revisione umana**: si ritenta con
  l'azione correttiva dichiarata dalla fase, entro ``retry.max_attempts``, se
  il codice è in ``retry.whitelist`` e ``retry.enabled`` è vero; poi ci si
  ferma, con un messaggio che distingue le due categorie;
* **degradazione automatica**: non ferma. La fase registra il ripiego con
  :meth:`~amplicon16s.steps.base.StepContext.degrada`, o lo dichiara con
  :meth:`~amplicon16s.steps.base.PipelineStep.ripiega`, e la degradazione
  finisce nel log e nel manifesto; l'esecuzione prosegue.

Quando si ferma, l'esecutore dichiara il **punto di ripresa**: la fase, il
codice, il messaggio operativo del catalogo, i tentativi fatti e il comando
per ripartire. Lo stampa chi lo invoca, e l'esecutore lo scrive in
``99_logs/punto_di_ripresa.json`` e ``.txt``, dove resta dopo la chiusura del
terminale. Un'esecuzione che arriva in fondo lo rimuove: non descriverebbe
piu' nulla.
"""

from __future__ import annotations

import json
import logging
from dataclasses import dataclass, field
from enum import StrEnum
from pathlib import Path
from typing import Any, Final

from amplicon16s.config.resolve import Registrazione, registra_risolta, risolvi
from amplicon16s.config.schema import Config
from amplicon16s.errors.catalog import voce
from amplicon16s.errors.exceptions import ErrorePipeline
from amplicon16s.gates.g01_g15 import (
    Contesto,
    ErroreGate,
    _controlla_coerenza,
    _g14_risorse_disponibili,
)
from amplicon16s.io_layer.artifacts import Fase
from amplicon16s.logging.logger import registra_errore
from amplicon16s.runner.graph import Passo
from amplicon16s.runner.project import ProjectRun, StatoPasso, Valutazione
from amplicon16s.runner.retry import Motivo, PoliticaRetry, applica, spiega_arresto, valore
from amplicon16s.steps.base import StepResult

__all__ = [
    "NOME_PUNTO_DI_RIPRESA",
    "Esecutore",
    "EsitoEsecuzione",
    "Conclusione",
    "PuntoDiRipresa",
]

#: Nome dei file del punto di ripresa sotto 99_logs, senza estensione.
NOME_PUNTO_DI_RIPRESA: Final = "punto_di_ripresa"


class Conclusione(StrEnum):
    """Come si è conclusa un'esecuzione."""

    #: Ogni fase richiesta è conclusa.
    COMPLETATA = "completata"
    #: Ferma su un errore, con il punto di ripresa dichiarato.
    ARRESTATA = "arrestata"
    #: Tutte le fasi realizzate sono concluse, ma la prossima non esiste ancora
    #: come codice.
    FASE_NON_REALIZZATA = "fase_non_realizzata"


@dataclass(frozen=True)
class PuntoDiRipresa:
    """Dove e perché l'esecuzione si è fermata, e come ripartire."""

    #: La fase fallita; ``None`` se l'arresto viene dai controlli di avvio e
    #: non c'e' una fase da eseguire.
    passo: Passo | None
    codice: str
    categoria: str
    sintesi: str
    azione: str
    dettaglio: str
    #: Perché non si è ritentato, o perché si è smesso.
    motivo: str
    tentativi: int
    tentativi_massimi: int
    comando: str
    #: ``"fase"`` o ``"controlli di avvio"``.
    origine: str = "fase"

    def come_documento(self) -> dict[str, Any]:
        return {
            "passo": str(self.passo) if self.passo else None,
            "origine": self.origine,
            "codice": self.codice,
            "categoria": self.categoria,
            "sintesi": self.sintesi,
            "azione": self.azione,
            "dettaglio": self.dettaglio,
            "motivo": self.motivo,
            "tentativi": self.tentativi,
            "tentativi_massimi": self.tentativi_massimi,
            "comando": self.comando,
        }

    def testo(self) -> str:
        """La dichiarazione per l'operatore."""
        dove = (
            f"nella fase {self.passo}" if self.origine == "fase" else "ai controlli di avvio"
            + (f", prima della fase {self.passo}" if self.passo else "")
        )
        righe = [
            f"ESECUZIONE FERMATA {dove}",
            f"  codice:     {self.codice} ({self.categoria})",
            f"  problema:   {self.sintesi}",
        ]
        if self.dettaglio:
            righe.append(f"  dettaglio:  {self.dettaglio}")
        righe += [
            f"  tentativi:  {self.tentativi} su {self.tentativi_massimi}",
            f"  perche':    {self.motivo}",
            f"  cosa fare:  {self.azione}",
            f"  ripartire:  {self.comando}",
        ]
        return "\n".join(righe)


@dataclass(frozen=True)
class EsitoEsecuzione:
    conclusione: Conclusione
    #: Le fasi eseguite in questa esecuzione, nell'ordine.
    eseguite: tuple[StepResult, ...] = ()
    punto: PuntoDiRipresa | None = None
    #: La fase non ancora realizzata a cui ci si è fermati.
    non_realizzata: Passo | None = None
    valutazione: Valutazione | None = field(default=None, repr=False)
    #: La configurazione registrata con cui si è eseguito; ``None`` se ci si è
    #: fermati ai controlli di avvio, prima di eseguire alcunche'.
    configurazione: Registrazione | None = None


class _Arresto(Exception):
    def __init__(self, punto: PuntoDiRipresa) -> None:
        super().__init__(punto.testo())
        self.punto = punto


class Esecutore:
    """Esegue le fasi da fare di un :class:`ProjectRun`."""

    def __init__(
        self,
        run: ProjectRun,
        *,
        comando_ripresa: str = "amplicon16s resume --config <configurazione>",
        fino_a: Passo | None = None,
        logger: logging.Logger | None = None,
    ) -> None:
        self.run = run
        self.comando_ripresa = comando_ripresa
        self.fino_a = fino_a
        self.log = logger or run.logger
        self.politica = PoliticaRetry.da_config(run.config)

        # Un aggiustamento dichiarato per un codice non ripetibile e' un
        # difetto della fase: non potrebbe mai essere applicato.
        for passo, fase in run.passi.items():
            for codice in fase.aggiustamenti:
                if not voce(codice).ammette_retry:
                    raise ValueError(
                        f"{passo} dichiara un aggiustamento per {codice}, che non "
                        "ammette il retry"
                    )

    @property
    def config(self) -> Config:
        return self.run.config

    # ----------------------------------------------------------------- #
    # Controlli di avvio                                                 #
    # ----------------------------------------------------------------- #

    def controlli_di_avvio(self, valutazione: Valutazione) -> None:
        """G15 e G14, prima di qualunque fase.

        Solleva :class:`ErroreGate` per G15, che è un errore di
        configurazione, e :class:`_Arresto` per G14, che dichiara un punto di
        ripresa: rimossa la causa, si riprende.
        """
        violazioni = _controlla_coerenza(risolvi(self.config))
        if violazioni:
            raise ErroreGate("G15", violazioni)

        violazioni, _ = _g14_risorse_disponibili(Contesto(self.config))
        if violazioni:
            prima = violazioni[0]
            dettaglio = "; ".join(v.dettaglio for v in violazioni)
            v = voce(prima.codice)
            punto = PuntoDiRipresa(
                passo=valutazione.prossima(),
                codice=prima.codice,
                categoria=v.categoria.value,
                sintesi=v.sintesi,
                azione=v.azione,
                dettaglio=dettaglio,
                motivo=(
                    "Le risorse della macchina si verificano a ogni avvio, anche "
                    "con S0 gia' conclusa."
                ),
                tentativi=0,
                tentativi_massimi=self.politica.tentativi_massimi,
                comando=self.comando_ripresa,
                origine="controlli di avvio",
            )
            raise _Arresto(punto)

    # ----------------------------------------------------------------- #
    # Esecuzione                                                         #
    # ----------------------------------------------------------------- #

    def esegui(self) -> EsitoEsecuzione:
        """Esegue le fasi da fare, fino alla fine, a un arresto o a ``fino_a``."""
        eseguite: list[StepResult] = []
        valutazione = self.run.valuta()
        configurazione: Registrazione | None = None
        try:
            self.controlli_di_avvio(valutazione)
            configurazione = self._registra_configurazione(valutazione)
            ultima: Passo | None = None
            while True:
                valutazione = self.run.valuta()
                passo = self._prossima(valutazione)
                if passo is None:
                    self._rimuovi_punto()
                    self.log.info("esecuzione completata", extra={"eseguite": len(eseguite)})
                    return EsitoEsecuzione(
                        Conclusione.COMPLETATA,
                        tuple(eseguite),
                        valutazione=valutazione,
                        configurazione=configurazione,
                    )
                if valutazione.situazioni[passo].stato is StatoPasso.NON_REALIZZATA:
                    self._rimuovi_punto()
                    self.log.warning(
                        f"{passo} non e' ancora realizzata: l'esecuzione si ferma qui",
                        extra={"passo": str(passo)},
                    )
                    return EsitoEsecuzione(
                        Conclusione.FASE_NON_REALIZZATA,
                        tuple(eseguite),
                        non_realizzata=passo,
                        valutazione=valutazione,
                        configurazione=configurazione,
                    )
                if passo is ultima:
                    # Appena conclusa e di nuovo da eseguire: la fase non ha
                    # scritto cio' che la rende conclusa. Ripeterla non servirebbe.
                    raise RuntimeError(
                        f"{passo} risulta ancora da eseguire subito dopo l'esecuzione: "
                        f"{valutazione.situazioni[passo].motivo}"
                    )
                eseguite.append(self._esegui_con_tentativi(passo, valutazione))
                ultima = passo
        except _Arresto as arresto:
            self._scrivi_punto(arresto.punto)
            return EsitoEsecuzione(
                Conclusione.ARRESTATA,
                tuple(eseguite),
                punto=arresto.punto,
                configurazione=configurazione,
            )

    def _registra_configurazione(self, valutazione: Valutazione) -> Registrazione:
        registrazione = registra_risolta(valutazione.risolta, self.config.io.out_root)
        extra = {
            "file": registrazione.percorso.name,
            "digest": registrazione.digest,
            "registrata_ora": registrazione.nuova,
            "differenze": list(registrazione.differenze),
        }
        if registrazione.nuova and registrazione.differenze:
            self.log.warning(
                "configurazione diversa dall'ultima registrata: registrata accanto",
                extra=extra,
            )
        else:
            self.log.info("configurazione in uso", extra=extra)
        return registrazione

    def _prossima(self, valutazione: Valutazione) -> Passo | None:
        """La prossima fase da eseguire, entro ``fino_a``."""
        ordine = self.run.grafo.ordine()
        for passo in valutazione.da_eseguire:
            if self.fino_a is None or ordine.index(passo) <= ordine.index(self.fino_a):
                return passo
        return None

    def _esegui_con_tentativi(self, passo: Passo, valutazione: Valutazione) -> StepResult:
        fase = self.run.fase(passo)
        config = self.config
        aggiustamenti: list[dict[str, Any]] = []
        tentativo = 1

        while True:
            contesto = self.run.contesto(
                passo,
                valutazione,
                config_effettiva=config if aggiustamenti else None,
                aggiustamenti=tuple(aggiustamenti),
            )
            try:
                risultato = fase.esegui(contesto)
            except ErrorePipeline as e:
                aggiustamento = fase.aggiustamenti.get(e.codice)
                decisione = self.politica.decidi(e.codice, tentativo, aggiustamento, config)
                registra_errore(
                    self.log,
                    e,
                    logging.WARNING if decisione.ritenta else logging.ERROR,
                )
                if not decisione.ritenta:
                    raise _Arresto(self._punto(passo, e, decisione.motivo, tentativo)) from e

                assert aggiustamento is not None
                prima = valore(config, aggiustamento.parametro)
                config = applica(config, aggiustamento.parametro, decisione.nuovo_valore)
                voce_aggiustamento = {
                    "codice": e.codice,
                    "tentativo_fallito": tentativo,
                    "parametro": aggiustamento.parametro,
                    "dichiarato": valore(self.config, aggiustamento.parametro),
                    "precedente": prima,
                    "usato": decisione.nuovo_valore,
                    "azione": aggiustamento.descrizione,
                }
                aggiustamenti.append(voce_aggiustamento)
                self.log.warning(
                    f"{passo}: nuovo tentativo con {aggiustamento.descrizione}",
                    extra={"passo": str(passo), **voce_aggiustamento},
                )
                tentativo += 1
                continue

            if not risultato.completata:
                raise RuntimeError(
                    f"{passo} si e' conclusa senza superare la verifica e senza un "
                    "codice d'errore: la fase deve sollevare il codice che la ferma"
                )
            return risultato

    def _punto(
        self, passo: Passo, e: ErrorePipeline, motivo: Motivo, tentativi: int
    ) -> PuntoDiRipresa:
        return PuntoDiRipresa(
            passo=passo,
            codice=e.codice,
            categoria=e.categoria.value,
            sintesi=e.voce.sintesi,
            azione=e.voce.azione,
            # Un gate respinto porta le violazioni in elenco, non nel dettaglio.
            dettaglio=(
                "; ".join(str(v) for v in e.violazioni)
                if isinstance(e, ErroreGate)
                else e.dettaglio
            ),
            motivo=spiega_arresto(
                e.codice, motivo, tentativi, self.politica.tentativi_massimi
            ),
            tentativi=tentativi,
            tentativi_massimi=self.politica.tentativi_massimi,
            comando=self.comando_ripresa,
        )

    # ----------------------------------------------------------------- #
    # Punto di ripresa su disco                                          #
    # ----------------------------------------------------------------- #

    def _percorsi_punto(self) -> tuple[Path, Path]:
        cartella = self.run.albero.cartella(Fase.LOGS)
        return (
            cartella / f"{NOME_PUNTO_DI_RIPRESA}.json",
            cartella / f"{NOME_PUNTO_DI_RIPRESA}.txt",
        )

    def _scrivi_punto(self, punto: PuntoDiRipresa) -> None:
        self.run.albero.prepara(Fase.LOGS)
        documento, testo = self._percorsi_punto()
        documento.write_text(
            json.dumps(punto.come_documento(), indent=2, ensure_ascii=False) + "\n",
            encoding="utf-8",
        )
        testo.write_text(punto.testo() + "\n", encoding="utf-8")
        self.log.error("punto di ripresa dichiarato", extra=punto.come_documento())

    def _rimuovi_punto(self) -> None:
        for percorso in self._percorsi_punto():
            percorso.unlink(missing_ok=True)
