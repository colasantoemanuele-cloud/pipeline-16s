"""Fase S14: la serializzazione, gli export e la validazione dell'oggetto finale.

Da S14 esce il risultato dell'intera pipeline: ``ps_final.rds``
(``out.serialization`` = ``rds``) in ``12_final/``, condivisa con S13, con
file del ponte e manifesto propri. Con ``out.export_flat`` anche gli export
piatti, per chi non usa lo stesso ambiente di calcolo: ``conteggi.tsv``,
``tassonomia.tsv``, ``metadati.tsv`` e ``sequenze.fasta``. Sono identici byte
per byte fra due esecuzioni e indipendenti dalle impostazioni locali
(``R/14_finale.R``). ``checksum.sha256`` elenca l'impronta di ogni file
consegnato, nella forma di ``sha256sum``.

**La validazione** (``E-S14-01``) e' strutturale: componenti allineati, solo
campioni biologici, nessun campione e nessuna variante senza letture, l'oggetto
riletto identico a quello scritto, l'oggetto ricostruito dai soli export
identico a quello serializzato, ``ps_filtrato.rds`` di S13 integro rispetto al
suo manifesto, i checksum dei file consegnati uguali a quelli registrati.

**Le due soglie del piano.** ``qc.min_reads_final`` vale per campione e si
applica in S13, dove un campione sotto soglia esce dall'oggetto finale con il
motivo: un singolo campione povero non ferma l'esecuzione, come non la ferma
un bianco azzerato. ``qc.min_frac_reads_retained`` vale per l'insieme dei
campioni finali: le loro letture nell'oggetto finale divise per le loro
letture senza chimere (S7). Una frazione bassa dice che decontaminazione e
filtri hanno tolto la maggior parte del segnale dei campioni tenuti, un
difetto dell'insieme e non di un campione: sotto la soglia, ``E-S14-01``.
"""

from __future__ import annotations

import hashlib
import json
from typing import ClassVar, Final

from amplicon16s.errors.exceptions import errore
from amplicon16s.io_layer.checksums import checksum_file
from amplicon16s.rbridge.runner import cartella_r, esegui_script
from amplicon16s.runner.graph import Passo
from amplicon16s.steps.base import PipelineStep, Produzione, StepContext
from amplicon16s.steps.s13_filtri import NOME_FILTRATO, NOME_RIEPILOGO

__all__ = ["NOME_CHECKSUM", "NOME_FINALE", "Serializzazione"]

NOME_SCRIPT: Final = "14_finale.R"
#: L'oggetto finale della pipeline.
NOME_FINALE: Final = "ps_final.rds"
#: Le impronte dei file consegnati, nella forma di sha256sum.
NOME_CHECKSUM: Final = "checksum.sha256"


class Serializzazione(PipelineStep):
    """S14: l'oggetto finale serializzato, gli export e la validazione."""

    passo: ClassVar[Passo] = Passo.S14
    versione: ClassVar[int] = 1
    script_r: ClassVar[str | None] = NOME_SCRIPT
    #: La serializzazione, l'orientamento verificato, gli export e la frazione
    #: minima delle letture trattenute dall'insieme dei campioni finali.
    parametri: ClassVar[tuple[str, ...]] = (
        "out.serialization", "out.taxa_are_rows", "out.export_flat",
        "qc.min_frac_reads_retained",
    )

    def calcola(self, contesto: StepContext) -> Produzione:
        """Verifica l'ingresso, esegue ``R/14_finale.R``, controlla la frazione
        trattenuta e scrive i checksum dei file consegnati.
        """
        config = contesto.config
        albero = contesto.albero
        cartella = albero.cartella(self.cartella)
        if config.out.serialization != "rds":
            raise RuntimeError(f"out.serialization {config.out.serialization!r} non e' realizzata")

        s13 = albero.manifesto_passo(Passo.S13, self.cartella)
        voce = None if s13 is None else next(
            (v for v in s13.artefatti if v["nome"] == NOME_FILTRATO), None)
        if voce is None or checksum_file(cartella / NOME_FILTRATO) != voce["checksum"]:
            raise errore("E-S14-01", f"{NOME_FILTRATO} di S13 manca o non corrisponde al suo manifesto")

        esito = esegui_script(
            cartella_r() / NOME_SCRIPT,
            {
                "filtrato": str(cartella / NOME_FILTRATO),
                "taxa_are_rows": config.out.taxa_are_rows,
                "export": config.out.export_flat,
            },
            albero,
            self.cartella,
            passo=self.passo,
            logger=contesto.logger,
        )

        riepilogo = json.loads((cartella / NOME_RIEPILOGO).read_text(encoding="utf-8"))
        letture = riepilogo["letture"]
        frazione = (
            letture["finali"] / letture["nonchimeric_finali"]
            if letture["nonchimeric_finali"] else 0.0
        )
        minima = config.qc.min_frac_reads_retained
        if frazione < minima:
            raise errore(
                "E-S14-01",
                f"i campioni finali trattengono {frazione:.4f} delle loro letture senza "
                f"chimere, meno di qc.min_frac_reads_retained ({minima})",
                frazione=round(frazione, 6),
            )

        righe = []
        for artefatto in esito.artefatti:
            ricalcolato = hashlib.sha256(artefatto.percorso.read_bytes()).hexdigest()
            if f"sha256:{ricalcolato}" != artefatto.checksum:
                raise errore("E-S14-01", f"il checksum di {artefatto.nome} non corrisponde")
            righe.append(f"{ricalcolato}  {artefatto.nome}\n")
        somme = albero.scrivi_testo(self.cartella, NOME_CHECKSUM, "".join(righe))

        metriche = {
            "file": [a.nome for a in esito.artefatti],
            "frazione_letture_trattenute": round(frazione, 6),
            "campioni": riepilogo["campioni"]["finali"],
            "varianti": riepilogo["varianti"]["finali"],
        }
        return Produzione(esito.artefatti + (somme,), metriche)
