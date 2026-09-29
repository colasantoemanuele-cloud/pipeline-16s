"""Fase S5: la tabella delle sequenze.

Riunisce le varianti di tutti i campioni in una tabella campioni x varianti
(``R/05_seqtab.R``, ``dada2::makeSequenceTable``, attraverso il ponte), in
``06_seqtab/``. E' la tabella prima della rimozione delle chimere: resta qui,
e S6 la legge senza ricopiarla.

**Il retry di E-S5-01 non serve, e la fase lo dichiara.** La tabella e' una
matrice densa di interi, campioni x varianti, che ``makeSequenceTable``
costruisce e poi riordina con una copia: la memoria che chiede dipende dal
numero di campioni e di varianti, non da ``run.batch_size`` (sul dataset di
riferimento, 960 campioni per 13.130 varianti, la matrice occupa 53 MB e il
picco ne aggiunge 157 alla lettura delle varianti), e costruirla a
lotti non la ridurrebbe, perche' la tabella finale e' comunque una. Il
numero di varianti distinte si controlla prima di allocarla, contro
``qc.max_asv_count``; sia la soglia superata sia la memoria esaurita durante
la costruzione sono E-S5-01, e in entrambi i casi un nuovo tentativo darebbe
lo stesso esito. La fase lo dichiara nell'errore (``RITENTARE_INUTILE``) e
l'esecuzione si ferma per la revisione.
"""

from __future__ import annotations

import json
from typing import ClassVar, Final

from amplicon16s.errors.exceptions import ErrorePipeline, errore
from amplicon16s.io_layer.artifacts import Fase
from amplicon16s.rbridge.runner import cartella_r, esegui_script
from amplicon16s.runner.graph import Passo
from amplicon16s.runner.retry import RITENTARE_INUTILE
from amplicon16s.steps.base import PipelineStep, Produzione, StepContext
from amplicon16s.steps.s02_filter import NOME_FILTRATE, leggi_conteggi
from amplicon16s.steps.s04_dada import NOME_VARIANTI

__all__ = ["NOME_TABELLA", "TabellaSequenze"]

NOME_SCRIPT: Final = "05_seqtab.R"
#: La tabella prima della rimozione delle chimere.
NOME_TABELLA: Final = "tabella.rds"
NOME_RIEPILOGO: Final = "tabella.json"

_INUTILE: Final = (
    "la tabella e' una matrice densa campioni x varianti, la cui dimensione non "
    "dipende da run.batch_size ne' da altri parametri che un nuovo tentativo "
    "possa cambiare"
)


class TabellaSequenze(PipelineStep):
    """S5: la tabella campioni x varianti."""

    passo: ClassVar[Passo] = Passo.S5
    passi_tracciamento: ClassVar[tuple[str, ...]] = ("tabella",)
    #: La soglia sul numero di varianti. Le varianti arrivano da S4.
    parametri: ClassVar[tuple[str, ...]] = ("qc.max_asv_count",)

    def calcola(self, contesto: StepContext) -> Produzione:
        """Esegue ``R/05_seqtab.R``, che costruisce la tabella campioni per varianti
        dalle varianti di S4.
        """
        if contesto.inventario is None:
            raise RuntimeError("S5 richiede l'inventario prodotto da S0")
        config = contesto.config
        letture = leggi_conteggi(contesto.albero.cartella(Fase.FILTERED) / NOME_FILTRATE)
        senza_letture = sorted(
            c.accession for c in contesto.inventario if letture.get(c.accession, 0) == 0
        )
        try:
            esito = esegui_script(
                cartella_r() / NOME_SCRIPT,
                {
                    "varianti": str(contesto.albero.cartella(Fase.ASV_INFERENCE) / NOME_VARIANTI),
                    "senza_letture": senza_letture,
                    "max_varianti": config.qc.max_asv_count,
                },
                contesto.albero,
                self.cartella,
                codice_memoria="E-S5-01",
                logger=contesto.logger,
            )
        except ErrorePipeline as e:
            if e.codice != "E-S5-01":
                raise
            raise errore(
                "E-S5-01", e.dettaglio, **{**e.contesto, RITENTARE_INUTILE: _INUTILE}
            ) from e

        riepilogo = json.loads(
            (contesto.albero.cartella(self.cartella) / NOME_RIEPILOGO).read_text(encoding="utf-8")
        )
        metriche = {k: riepilogo[k] for k in ("campioni", "varianti", "letture")}
        return Produzione(esito.artefatti, metriche)
