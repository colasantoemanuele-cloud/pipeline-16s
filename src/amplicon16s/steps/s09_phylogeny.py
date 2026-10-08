"""Fase S9: la filogenesi delle varianti finali.

Facoltativa: con ``phylo.enabled`` falso, il predefinito, è disattivata. Quando
è attiva costruisce l'albero filogenetico delle varianti dell'oggetto filtrato
di S13, quelle che restano dopo decontaminazione e filtri e che S14 consegna, e
lo scrive in ``09_phylogeny/`` in formato Newick, con le foglie che portano gli
identificativi delle varianti. S14 lo aggiunge all'oggetto finale e lo esporta.

**Dopo S13, non dopo S7.** Sulle varianti di S7 l'albero conterrebbe anche
contaminanti, mitocondri, cloroplasti e varianti rare che i filtri tolgono:
molte più sequenze delle varianti finali, anche oltre ``phylo.max_seqs``, per
un albero che nessuna analisi userebbe. Il numero della fase e della cartella resta quello del piano; l'ordine
di esecuzione è quello del grafo.

**Le due guardie** sono qui, prima di avviare qualunque calcolo, sul numero di
varianti finali letto dal riepilogo di S13. ``E-S9-01``: oltre
``phylo.max_seqs`` il costo della massima verosimiglianza, che cresce
rapidamente con il numero di sequenze, non è più sostenibile. ``E-S9-02``: con
meno di :data:`VARIANTI_MINIME` varianti un albero radicato non ha una
topologia da stimare. In entrambi i casi un albero richiesto esplicitamente
non si salta in silenzio: l'esecuzione si ferma, a revisione umana.

**Allineamento, modello, ricerca e radicamento** sono in ``R/09_phylogeny.R``:
allineamento multiplo con DECIPHER (``phylo.aligner``), massima verosimiglianza
con phangorn sul modello ``phylo.model`` (GTR+G+I), ricerca deterministica per
scambi fra rami vicini da un albero neighbor-joining, radicamento al punto
medio. Il seme viene da ``run.seed``. L'albero è identico byte per byte fra due
esecuzioni e con un numero di thread diverso: ``run.threads`` accelera solo
l'allineamento, e non entra nell'impronta.

**Un limite da conoscere.** Su letture corte, di un solo tratto del gene, un
albero costruito da zero è debolmente risolto: le relazioni profonde non sono
sostenute dal dato. Per questo la fase è disattivata per difetto.
"""

from __future__ import annotations

import json
from typing import ClassVar, Final

from amplicon16s.config.schema import thread_effettivi
from amplicon16s.errors.exceptions import errore
from amplicon16s.io_layer.artifacts import Fase
from amplicon16s.rbridge.runner import cartella_r, esegui_script
from amplicon16s.runner.graph import Passo
from amplicon16s.steps.base import PipelineStep, Produzione, StepContext
from amplicon16s.steps.s13_filtri import NOME_FILTRATO, NOME_RIEPILOGO

__all__ = [
    "Filogenesi",
    "NOME_ALBERO",
    "NOME_ALLINEAMENTO",
    "NOME_FILOGENESI",
    "VARIANTI_MINIME",
]

NOME_SCRIPT: Final = "09_phylogeny.R"
#: L'albero radicato delle varianti finali, in formato Newick: S14 lo legge.
NOME_ALBERO: Final = "albero.nwk"
NOME_ALLINEAMENTO: Final = "allineamento.fasta"
NOME_FILOGENESI: Final = "filogenesi.json"
#: Sotto questo numero di varianti un albero radicato non ha topologia da
#: stimare: con tre foglie l'albero senza radice e' uno solo. E' una proprieta'
#: del problema, non una soglia da configurare.
VARIANTI_MINIME: Final = 4


class Filogenesi(PipelineStep):
    """S9: l'albero filogenetico radicato delle varianti finali."""

    passo: ClassVar[Passo] = Passo.S9
    versione: ClassVar[int] = 1
    script_r: ClassVar[str | None] = NOME_SCRIPT
    #: Il gruppo phylo per intero (attivazione, limite, allineatore, modello) e
    #: il seme, fissato prima di ogni calcolo.
    parametri: ClassVar[tuple[str, ...]] = ("phylo", "run.seed")

    def calcola(self, contesto: StepContext) -> Produzione:
        """Verifica che le varianti finali siano abbastanza per un albero e non
        piu' di ``phylo.max_seqs``, poi esegue ``R/09_phylogeny.R``.
        """
        config = contesto.config
        albero = contesto.albero
        finale = albero.cartella(Fase.FINAL_INTERMEDI)
        varianti = json.loads(
            (finale / NOME_RIEPILOGO).read_text(encoding="utf-8")
        )["varianti"]["finali"]
        # Prima di allocare il calcolo: le due guardie non devono costare nulla.
        if varianti < VARIANTI_MINIME:
            raise errore(
                "E-S9-02",
                f"{varianti} varianti finali, meno delle {VARIANTI_MINIME} che servono a un albero",
                varianti=varianti,
                minimo=VARIANTI_MINIME,
            )
        if varianti > config.phylo.max_seqs:
            raise errore(
                "E-S9-01",
                f"{varianti} varianti finali, piu' di phylo.max_seqs ({config.phylo.max_seqs})",
                varianti=varianti,
                massimo=config.phylo.max_seqs,
            )

        esito = esegui_script(
            cartella_r() / NOME_SCRIPT,
            {
                "filtrato": str(finale / NOME_FILTRATO),
                "allineatore": config.phylo.aligner,
                "modello": config.phylo.model,
                "seme": config.run.seed,
                "processi": thread_effettivi(config),
            },
            albero,
            self.cartella,
            passo=self.passo,
            tempo_massimo_s=contesto.config.run.r_timeout_s,
            logger=contesto.logger,
        )
        riepilogo = json.loads(
            (albero.cartella(self.cartella) / NOME_FILOGENESI).read_text(encoding="utf-8")
        )
        metriche = {
            k: riepilogo[k] for k in (
                "varianti", "colonne_allineamento", "modello", "radicato",
                "log_verosimiglianza",
            )
        }
        return Produzione(esito.artefatti, metriche)
