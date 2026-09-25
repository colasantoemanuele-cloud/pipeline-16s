"""L'esecuzione della pipeline e il punto a cui è arrivata.

:class:`ProjectRun` tiene insieme configurazione risolta, inventario e albero
di output, e sa dire per ogni fase se è conclusa, disattivata o da eseguire, e
quale viene dopo.

**Lo stato è il filesystem.** Non c'è un database né un file di avanzamento:
lo stato si rilegge dal disco. Un terzo canale andrebbe tenuto allineato con
la realtà su disco, e il giorno in cui divergesse direbbe il falso proprio
quando serve riprendere.

**Una valutazione, molte domande.** :meth:`ProjectRun.valuta` legge lo stato
una volta e restituisce una :class:`Valutazione`, a cui si pongono tutte le
domande — concluse, disattivate, prossima — senza ricalcolare nulla. Dentro
una valutazione il checksum di ogni artefatto è calcolato una volta sola. Le
proprieta' omonime di :class:`ProjectRun` sono scorciatoie per una domanda
isolata, e ciascuna compie una valutazione nuova: chi ne pone piu' di una
sullo stesso stato usa la valutazione.

Una fase è **conclusa** quando, in quest'ordine:

1. le fasi da cui dipende sono concluse;
2. il suo manifesto esiste e non è stato alterato;
3. i suoi artefatti ci sono e corrispondono ai checksum;
4. è stata calcolata su ciò che le verrebbe dato adesso: stessi parametri,
   stesse impronte delle fasi a monte, stessi dati esterni.

La condizione 1 propaga l'invalidità a valle: se una fase va rifatta, tutte
quelle che ne dipendono vanno rifatte con lei. La 4 intercetta ciò che la 3
non vede, cioè artefatti integri ma calcolati su ingressi che non esistono
piu'.

L'esecuzione in sequenza delle fasi da eseguire, con i tentativi ripetuti, e
il collegamento con il comando ``resume`` non sono qui.
"""

from __future__ import annotations

import logging
from collections.abc import Mapping
from dataclasses import dataclass
from enum import StrEnum
from pathlib import Path
from typing import Any

from amplicon16s.config.resolve import ConfigRisolta, risolvi
from amplicon16s.config.schema import Config
from amplicon16s.io_layer.artifacts import AlberoOutput, ManifestoPassoNonValido
from amplicon16s.io_layer.checksums import checksum_file
from amplicon16s.logging.logger import ottieni
from amplicon16s.metadata.models import Inventario
from amplicon16s.runner.graph import GRAFO, Grafo, Passo
from amplicon16s.steps.base import PipelineStep, StepContext
from amplicon16s.steps.s00_validate import ValidazioneIngressi, leggi_inventario

__all__ = [
    "ProjectRun",
    "Situazione",
    "StatoPasso",
    "Valutazione",
    "passi_realizzati",
]


def passi_realizzati() -> dict[Passo, PipelineStep]:
    """Le fasi che esistono come codice. Crescono man mano che vengono scritte."""
    return {Passo.S0: ValidazioneIngressi()}


class StatoPasso(StrEnum):
    """Dove si trova una fase rispetto a questa esecuzione."""

    COMPLETATA = "completata"
    DA_ESEGUIRE = "da_eseguire"
    #: Esclusa dalla configurazione: non è lavoro mancante.
    DISATTIVATA = "disattivata"
    #: Prevista dal grafo ma non ancora scritta come codice.
    NON_REALIZZATA = "non_realizzata"


@dataclass(frozen=True)
class Situazione:
    """Lo stato di una fase, con il motivo."""

    passo: Passo
    stato: StatoPasso
    motivo: str
    #: Impronta del manifesto, per le sole fasi concluse.
    impronta: str | None = None


def _differenze(registrata: Mapping[str, Any], attesa: Mapping[str, Any]) -> str:
    """In che cosa gli ingressi registrati differiscono da quelli di adesso."""
    motivi = []
    if registrata.get("configurazione") != attesa["configurazione"]:
        motivi.append("configurazione cambiata")
    prima = registrata.get("a_monte", {})
    ora = attesa["a_monte"]
    if set(prima) != set(ora):
        motivi.append(f"dipendenze cambiate: {sorted(prima)} -> {sorted(ora)}")
    else:
        ricalcolate = sorted(p for p in ora if prima[p] != ora[p])
        if ricalcolate:
            motivi.append(f"ricalcolate a monte: {', '.join(ricalcolate)}")
    if registrata.get("dati_esterni") != attesa["dati_esterni"]:
        motivi.append("dati di ingresso cambiati")
    return "; ".join(motivi) or "calcolata su ingressi diversi"


class _Checksum:
    """Checksum calcolati durante una valutazione, ciascuno una volta sola."""

    def __init__(self) -> None:
        self._noti: dict[Path, str | None] = {}

    def __call__(self, percorso: Path) -> str | None:
        if percorso not in self._noti:
            self._noti[percorso] = (
                checksum_file(percorso) if percorso.is_file() else None
            )
        return self._noti[percorso]


@dataclass(frozen=True)
class Valutazione:
    """Lo stato di un'esecuzione, letto dal disco in un momento dato.

    Rispondere alle domande non rilegge il disco: se nel frattempo un file
    cambia, la valutazione non lo sa, e per saperlo se ne chiede una nuova.
    """

    situazioni: Mapping[Passo, Situazione]
    #: L'inventario prodotto da S0, se S0 risulta conclusa.
    inventario: Inventario | None
    #: La configurazione risolta, con i derivati dai dati quando S0 è conclusa.
    risolta: ConfigRisolta

    def _con_stato(self, *stati: StatoPasso) -> tuple[Passo, ...]:
        return tuple(p for p, s in self.situazioni.items() if s.stato in stati)

    @property
    def completate(self) -> tuple[Passo, ...]:
        return self._con_stato(StatoPasso.COMPLETATA)

    @property
    def disattivate(self) -> tuple[Passo, ...]:
        return self._con_stato(StatoPasso.DISATTIVATA)

    @property
    def da_eseguire(self) -> tuple[Passo, ...]:
        """Il lavoro che manca, in ordine, comprese le fasi non ancora realizzate."""
        return self._con_stato(StatoPasso.DA_ESEGUIRE, StatoPasso.NON_REALIZZATA)

    def prossima(self) -> Passo | None:
        """La prossima fase da eseguire, o ``None`` se l'esecuzione è completa."""
        mancanti = self.da_eseguire
        return mancanti[0] if mancanti else None

    @property
    def completa(self) -> bool:
        """Ogni fase attiva è conclusa; le disattivate non sono lavoro mancante."""
        return not self.da_eseguire


class ProjectRun:
    """Un'esecuzione della pipeline sotto una cartella di output."""

    def __init__(
        self,
        config: Config,
        *,
        passi: Mapping[Passo, PipelineStep] | None = None,
        grafo: Grafo = GRAFO,
        logger: logging.Logger | None = None,
    ) -> None:
        self.config = config
        self.grafo = grafo
        self.passi = dict(passi_realizzati() if passi is None else passi)
        self.logger = logger or ottieni("run")
        self.albero = AlberoOutput(config.io.out_root)
        for passo, fase in self.passi.items():
            if fase.passo is not passo:
                raise ValueError(f"{type(fase).__name__} realizza {fase.passo}, non {passo}")

    def _risolta(
        self, inventario: Inventario | None, config: Config | None = None
    ) -> ConfigRisolta:
        risolta = risolvi(config or self.config)
        if inventario is not None:
            risolta = risolta.con_campioni_biologici(inventario.denominatore_prevalenza())
        return risolta

    # ----------------------------------------------------------------- #
    # Valutazione                                                        #
    # ----------------------------------------------------------------- #

    def valuta(self) -> Valutazione:
        """Legge dal disco lo stato di ogni fase, nell'ordine del grafo."""
        checksum = _Checksum()
        situazioni: dict[Passo, Situazione] = {}
        impronte: dict[Passo, str | None] = {}
        inventario: Inventario | None = None
        risolta = risolvi(self.config)

        for passo in self.grafo.ordine():
            nodo = self.grafo.nodo(passo)
            impronte[passo] = None

            if not nodo.attiva(self.config):
                situazioni[passo] = Situazione(
                    passo, StatoPasso.DISATTIVATA, f"{nodo.parametro_attivazione} falso"
                )
                continue

            fase = self.passi.get(passo)
            if fase is None:
                situazioni[passo] = Situazione(
                    passo, StatoPasso.NON_REALIZZATA, "fase non ancora realizzata"
                )
                continue

            a_monte = {
                d: impronte[d] for d in self.grafo.dipendenze_attive(passo, self.config)
            }
            situazioni[passo] = self._valuta_passo(fase, risolta, a_monte, checksum)
            impronte[passo] = situazioni[passo].impronta

            # Le fasi dopo S0 si calcolano con i derivati dai dati.
            if passo is Passo.S0 and situazioni[passo].stato is StatoPasso.COMPLETATA:
                inventario = leggi_inventario(self.config)
                risolta = self._risolta(inventario)

        return Valutazione(situazioni, inventario, risolta)

    def _valuta_passo(
        self,
        fase: PipelineStep,
        risolta: ConfigRisolta,
        a_monte: Mapping[Passo, str | None],
        checksum: _Checksum,
    ) -> Situazione:
        passo = fase.passo
        da_eseguire = StatoPasso.DA_ESEGUIRE

        incomplete = [str(p) for p, impronta in a_monte.items() if impronta is None]
        if incomplete:
            return Situazione(passo, da_eseguire, f"a monte da eseguire: {', '.join(incomplete)}")

        try:
            manifesto = self.albero.manifesto_passo(passo, fase.cartella)
        except ManifestoPassoNonValido as e:
            return Situazione(passo, da_eseguire, str(e))
        if manifesto is None:
            return Situazione(passo, da_eseguire, "nessun manifesto: mai conclusa")

        non_integri = self.albero.non_integri_del_passo(manifesto, checksum)
        if non_integri:
            return Situazione(
                passo, da_eseguire, f"artefatti mancanti o alterati: {', '.join(non_integri)}"
            )

        attesa = fase.calcolata_su(risolta, a_monte)
        if manifesto.calcolata_su != attesa:
            return Situazione(passo, da_eseguire, _differenze(manifesto.calcolata_su, attesa))

        return Situazione(passo, StatoPasso.COMPLETATA, "conclusa", manifesto.impronta)

    # ----------------------------------------------------------------- #
    # Scorciatoie: una valutazione nuova per ogni domanda                #
    # ----------------------------------------------------------------- #

    def situazione(self) -> dict[Passo, Situazione]:
        """Lo stato di ogni fase, nell'ordine del grafo."""
        return dict(self.valuta().situazioni)

    @property
    def completate(self) -> tuple[Passo, ...]:
        return self.valuta().completate

    @property
    def disattivate(self) -> tuple[Passo, ...]:
        return self.valuta().disattivate

    @property
    def da_eseguire(self) -> tuple[Passo, ...]:
        return self.valuta().da_eseguire

    def prossima(self) -> Passo | None:
        return self.valuta().prossima()

    @property
    def completa(self) -> bool:
        return self.valuta().completa

    @property
    def inventario(self) -> Inventario | None:
        return self.valuta().inventario

    @property
    def risolta(self) -> ConfigRisolta:
        return self.valuta().risolta

    # ----------------------------------------------------------------- #
    # Contesto di una fase                                               #
    # ----------------------------------------------------------------- #

    def contesto(
        self,
        passo: Passo,
        valutazione: Valutazione | None = None,
        *,
        config_effettiva: Config | None = None,
        aggiustamenti: tuple[Mapping[str, Any], ...] = (),
    ) -> StepContext:
        """Il contesto con cui eseguire una fase, secondo la valutazione data.

        Senza valutazione se ne compie una nuova. Le dipendenze non concluse
        compaiono con impronta ``None``: la fase lo rileva nella verifica dei
        prerequisiti e si rifiuta di girare.

        ``config_effettiva`` è la configurazione con i parametri cambiati da
        un'azione correttiva, descritti in ``aggiustamenti``: la fase calcola
        con quella, ma resta giudicata sulla configurazione dichiarata.
        """
        valutazione = valutazione or self.valuta()
        a_monte = {
            d: valutazione.situazioni[d].impronta
            for d in self.grafo.dipendenze_attive(passo, self.config)
        }
        effettiva = (
            self._risolta(valutazione.inventario, config_effettiva)
            if config_effettiva is not None
            else valutazione.risolta
        )
        return StepContext(
            risolta=effettiva,
            albero=self.albero,
            logger=self.logger,
            inventario=valutazione.inventario,
            a_monte=a_monte,
            dichiarata=valutazione.risolta if config_effettiva is not None else None,
            aggiustamenti=aggiustamenti,
        )

    def fase(self, passo: Passo) -> PipelineStep:
        """La fase registrata per ``passo``."""
        try:
            return self.passi[passo]
        except KeyError:
            raise LookupError(f"{passo} non e' ancora realizzata") from None
