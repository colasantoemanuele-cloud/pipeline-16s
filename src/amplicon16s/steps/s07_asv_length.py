"""Fase S7: il filtro di lunghezza delle varianti.

Tiene della tabella senza chimere le sole varianti di lunghezza fra
``asv.len_min`` e ``asv.len_max`` (``R/07_asv_length.R``), derivati da
``filter.truncLen - filter.trimLeft`` con la tolleranza ``asv.len_tol``, e
scrive in ``07_chimera/`` accanto a S6.

Le letture sono troncate a lunghezza fissa e, in single-end, nessuna fusione
di coppie ne cambia la lunghezza: le varianti hanno tutte la lunghezza del
troncamento, e con la tolleranza zero il filtro non toglie nulla. Sul dataset
di riferimento e' cosi'. La fase resta per i dati in cui la lunghezza varia,
e ``lunghezze.tsv`` mostra comunque la distribuzione delle lunghezze. Che il
filtro tolga le varianti fuori intervallo e' verificato con varianti
sintetiche di lunghezze diverse.

Con S7 il tracciamento delle letture e' completo, un passo per fase da S1:
grezze, prefiltro e filtrate (S2), denoised (S4), tabella (S5), senza
chimere (S6), lunghezza (S7). Ogni fase scrive il proprio file, e la tabella
si ricompone con :func:`amplicon16s.runner.tracciamento.ricomponi`.
"""

from __future__ import annotations

from typing import ClassVar, Final

from amplicon16s.io_layer.artifacts import Fase
from amplicon16s.rbridge.runner import cartella_r, esegui_script
from amplicon16s.runner.graph import Passo
from amplicon16s.steps.base import PipelineStep, Produzione, StepContext
from amplicon16s.steps.s02_filter import NOME_FILTRATE, leggi_conteggi
from amplicon16s.steps.s06_chimera import NOME_SENZA_CHIMERE

__all__ = ["FiltroLunghezza", "NOME_TABELLA_ASV"]

NOME_SCRIPT: Final = "07_asv_length.R"
#: La tabella delle varianti che le fasi successive leggono.
NOME_TABELLA_ASV: Final = "tabella_asv.rds"


class FiltroLunghezza(PipelineStep):
    """S7: le varianti di lunghezza ammessa."""

    passo: ClassVar[Passo] = Passo.S7
    passi_tracciamento: ClassVar[tuple[str, ...]] = ("lunghezza",)
    #: I tre parametri da cui discendono asv.len_min e asv.len_max.
    parametri: ClassVar[tuple[str, ...]] = ("asv", "filter.truncLen", "filter.trimLeft")

    def calcola(self, contesto: StepContext) -> Produzione:
        if contesto.inventario is None:
            raise RuntimeError("S7 richiede l'inventario prodotto da S0")
        derivati = contesto.risolta.derivati
        letture = leggi_conteggi(contesto.albero.cartella(Fase.FILTERED) / NOME_FILTRATE)
        esito = esegui_script(
            cartella_r() / NOME_SCRIPT,
            {
                "tabella": str(contesto.albero.cartella(Fase.CHIMERA) / NOME_SENZA_CHIMERE),
                "senza_letture": sorted(
                    c.accession for c in contesto.inventario if letture.get(c.accession, 0) == 0
                ),
                "len_min": derivati.asv_len_min,
                "len_max": derivati.asv_len_max,
            },
            contesto.albero,
            self.cartella,
            logger=contesto.logger,
        )
        righe = (
            (contesto.albero.cartella(self.cartella) / "lunghezze.tsv")
            .read_text(encoding="utf-8").splitlines()[1:]
        )
        per_lunghezza = [r.split("\t") for r in righe]
        metriche = {
            "len_min": derivati.asv_len_min,
            "len_max": derivati.asv_len_max,
            "varianti_ammesse": sum(int(v) for _, v, _, a in per_lunghezza if a == "si"),
            "varianti_escluse": sum(int(v) for _, v, _, a in per_lunghezza if a == "no"),
            "letture_escluse": sum(int(float(n)) for _, _, n, a in per_lunghezza if a == "no"),
        }
        return Produzione(esito.artefatti, metriche)
