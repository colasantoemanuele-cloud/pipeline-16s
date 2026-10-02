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
non scrive nulla; una ripresa con un digest diverso (anche solo per
``run.threads`` o un parametro di retry, fuori dall'impronta dei risultati ma
non dal digest) scrive la nuova versione accanto alle precedenti. Il log
dice a ogni avvio con quale versione si esegue.

Poi esegue, una alla volta, le fasi che la :class:`Valutazione` dà da
eseguire. L'esito di una fase fallita dipende dalla categoria del suo codice:

* **revisione umana**: ci si ferma subito;
* **retry automatico** e **retry poi revisione umana**: si ritenta con
  l'azione correttiva dichiarata dalla fase, entro ``retry.max_attempts``, se
  il codice è in ``retry.whitelist`` e ``retry.enabled`` è vero, e se la fase
  non ha dichiarato nell'errore che un nuovo tentativo darebbe lo stesso
  esito; poi ci si ferma, con un messaggio che distingue le due categorie;
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
    _g12_riferimento_verificato,
    _g14_risorse_disponibili,
)
from amplicon16s.io_layer.artifacts import Fase
from amplicon16s.logging.logger import registra_errore
from amplicon16s.runner.graph import Passo
from amplicon16s.runner.project import ProjectRun, StatoPasso, Valutazione
from amplicon16s.runner.retry import (
    RITENTARE_INUTILE,
    Motivo,
    PoliticaRetry,
    applica,
    spiega_arresto,
    valore,
)
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
        """Il punto di ripresa come documento JSON, scritto in ``99_logs``."""
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
    """L'esito di un'esecuzione: come si è conclusa, le fasi eseguite e l'eventuale
    punto di ripresa.
    """
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
    """Arresto ai controlli di avvio, che porta con sé il punto di ripresa dichiarato.
    """
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
        """La configurazione dichiarata dell'esecuzione."""
        return self.run.config

    # ----------------------------------------------------------------- #
    # Controlli di avvio                                                 #
    # ----------------------------------------------------------------- #

    def controlli_di_avvio(self, valutazione: Valutazione) -> None:
        """G15, G14 e G12, prima di qualunque fase.

        Solleva :class:`ErroreGate` per G15, che è un errore di
        configurazione, e :class:`_Arresto` per G14 e G12, che dichiarano un
        punto di ripresa: rimossa la causa, si riprende. G12 verifica il
        riferimento tassonomico contro tax.ref_md5 a ogni avvio, anche con S0
        conclusa: un file alterato dopo S0 non arriva a S8.
        """
        violazioni = _controlla_coerenza(risolvi(self.config))
        if violazioni:
            raise ErroreGate("G15", violazioni)

        completa = Contesto(self.config)
        for controllo, motivo in (
            (_g14_risorse_disponibili,
             "Le risorse della macchina si verificano a ogni avvio, anche con S0 gia' conclusa."),
            (_g12_riferimento_verificato,
             "Il riferimento tassonomico si verifica contro tax.ref_md5 a ogni avvio, "
             "anche con S0 gia' conclusa."),
        ):
            violazioni, _ = controllo(completa)
            if violazioni:
                prima = violazioni[0]
                v = voce(prima.codice)
                raise _Arresto(PuntoDiRipresa(
                    passo=valutazione.prossima(),
                    codice=prima.codice,
                    categoria=v.categoria.value,
                    sintesi=v.sintesi,
                    azione=v.azione,
                    dettaglio="; ".join(x.dettaglio for x in violazioni),
                    motivo=motivo,
                    tentativi=0,
                    tentativi_massimi=self.politica.tentativi_massimi,
                    comando=self.comando_ripresa,
                    origine="controlli di avvio",
                ))

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
            # Una provenienza diversa a parita' di versione non rifa' la fase:
            # si segnala una volta per esecuzione, e resta nel resoconto.
            for passo, situazione in valutazione.situazioni.items():
                if situazione.avviso is not None:
                    self.log.warning(
                        f"{passo} conclusa con una provenienza diversa: {situazione.avviso}",
                        extra={"passo": str(passo), "avviso": situazione.avviso},
                    )
            ultima: Passo | None = None
            while True:
                valutazione = self.run.valuta()
                passo = self._prossima(valutazione)
                if passo is None:
                    self._rimuovi_punto()
                    if self.fino_a is None:
                        self._rimuovi_temporanei(valutazione)
                        valutazione = self.run.valuta()
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

    def _rimuovi_temporanei(self, valutazione: Valutazione) -> None:
        """A esecuzione conclusa, gli artefatti che la configurazione non vuole tenere.

        Solo qui, quando ogni fase e' conclusa: prima, chi li consuma potrebbe
        ancora averne bisogno. La rimozione e' registrata accanto al manifesto
        della fase, che resta conclusa.
        """
        for passo, situazione in valutazione.situazioni.items():
            if situazione.stato is not StatoPasso.COMPLETATA or passo not in self.run.passi:
                continue
            fase = self.run.passi[passo]
            manifesto = self.run.albero.manifesto_passo(passo, fase.cartella)
            if manifesto is None:
                continue
            gia_rimossi = self.run.albero.rimossi_del_passo(manifesto)
            nomi = [n for n in fase.artefatti_temporanei(manifesto, self.config) if n not in gia_rimossi]
            if not nomi:
                continue
            liberati = sum(
                (self.run.albero.cartella(fase.cartella) / n).stat().st_size
                for n in nomi if (self.run.albero.cartella(fase.cartella) / n).exists()
            )
            self.run.albero.rimuovi_artefatti(
                manifesto, nomi,
                "rimossi a esecuzione conclusa, come chiede la configurazione",
            )
            self.log.info(
                f"{passo}: artefatti temporanei rimossi",
                extra={"passo": str(passo), "rimossi": len(nomi), "byte_liberati": liberati},
            )

    def _registra_configurazione(self, valutazione: Valutazione) -> Registrazione:
        """Registra la configurazione risolta nella cartella di output e segnala se
        differisce dall'ultima.
        """
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
        """La prossima fase da eseguire, entro ``fino_a``.

        Con ``fino_a`` si eseguono la fase indicata e i suoi antenati secondo
        le dipendenze del grafo, non tutte le fasi che la precedono
        nell'ordine: S2 non dipende da S1, e chiedere S2 non esegue S1.
        """
        if self.fino_a is None:
            richieste = None
        else:
            richieste = {self.fino_a, *self.run.grafo.antenati(self.fino_a, self.config)}
        for passo in valutazione.da_eseguire:
            if richieste is None or passo in richieste:
                return passo
        return None

    def _esegui_con_tentativi(self, passo: Passo, valutazione: Valutazione) -> StepResult:
        """Esegue una fase applicando la politica di retry, e registra ogni
        aggiustamento.
        """
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
                inutile = e.contesto.get(RITENTARE_INUTILE)
                decisione = self.politica.decidi(
                    e.codice, tentativo, aggiustamento, config, inutile
                )
                registra_errore(
                    self.log,
                    e,
                    logging.WARNING if decisione.ritenta else logging.ERROR,
                )
                if not decisione.ritenta:
                    raise _Arresto(self._punto(passo, e, decisione.motivo, tentativo)) from e

                assert aggiustamento is not None
                if aggiustamento.invariato:
                    # Dichiarato transitorio: si ritenta con la configurazione
                    # di prima, e lo si registra come ogni altro aggiustamento.
                    voce_aggiustamento = {
                        "codice": e.codice,
                        "tentativo_fallito": tentativo,
                        "parametro": None,
                        "dichiarato": None,
                        "precedente": None,
                        "usato": None,
                        "azione": aggiustamento.descrizione,
                    }
                else:
                    assert aggiustamento.parametro is not None
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
        """Il punto di ripresa che corrisponde a un errore della fase indicata."""
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
                e.codice, motivo, tentativi, self.politica.tentativi_massimi,
                e.contesto.get(RITENTARE_INUTILE),
            ),
            tentativi=tentativi,
            tentativi_massimi=self.politica.tentativi_massimi,
            comando=self.comando_ripresa,
        )

    # ----------------------------------------------------------------- #
    # Punto di ripresa su disco                                          #
    # ----------------------------------------------------------------- #

    def _percorsi_punto(self) -> tuple[Path, Path]:
        """I percorsi JSON e testuale del punto di ripresa in ``99_logs``."""
        cartella = self.run.albero.cartella(Fase.LOGS)
        return (
            cartella / f"{NOME_PUNTO_DI_RIPRESA}.json",
            cartella / f"{NOME_PUNTO_DI_RIPRESA}.txt",
        )

    def _scrivi_punto(self, punto: PuntoDiRipresa) -> None:
        """Scrive il punto di ripresa nei due formati e lo registra nel log."""
        self.run.albero.prepara(Fase.LOGS)
        documento, testo = self._percorsi_punto()
        documento.write_text(
            json.dumps(punto.come_documento(), indent=2, ensure_ascii=False) + "\n",
            encoding="utf-8",
        )
        testo.write_text(punto.testo() + "\n", encoding="utf-8")
        self.log.error("punto di ripresa dichiarato", extra=punto.come_documento())

    def _rimuovi_punto(self) -> None:
        """Rimuove il punto di ripresa di un'esecuzione precedente, se presente."""
        for percorso in self._percorsi_punto():
            percorso.unlink(missing_ok=True)
