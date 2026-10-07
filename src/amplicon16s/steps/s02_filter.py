"""Fase S2: filtro e troncamento delle letture.

È la prima trasformazione irreversibile del dato: tronca le letture a
``filter.truncLen`` e scarta quelle che non superano il filtro. Tutte le fasi
successive lavorano sulle sue uscite, in ``03_filtered/``
(``R/02_filter.R``, ``dada2::filterAndTrim``, attraverso il ponte).

**Gli archivi si verificano per intero prima del filtro.** Un archivio gzip
troncato non fa fallire ``filterAndTrim``: viene letto fino a dove arriva, e
le letture perse spariscono senza traccia (verificato su un file della
versione ridotta tagliato a meta'). S2 quindi decomprime ogni file per intero, in parallelo, e un archivio incompleto e'
E-S2-03, errore di lettura: se era transitorio la rilettura riesce, se il file
e' davvero corrotto la fase si ferma dopo i tentativi ammessi.

**I controlli sul risultato valgono per i biologici e i positivi.** Come per
G10, un controllo si applica alle classi per cui ha senso:

* E-S2-01, campioni azzerati oltre ``qc.max_zeroed_samples``: un campione
  biologico azzerato e' un campione perso, un controllo positivo azzerato mette
  in dubbio il lotto. Un controllo negativo azzerato no: e' un bianco pulito, e
  fermare la pipeline per quello sarebbe sbagliato. Si registra e si prosegue.
* E-S2-02, perdita media oltre ``qc.max_frac_lost_filter``: domanda se i
  parametri del filtro sono adatti ai campioni analizzati. La perdita di un
  bianco non risponde a quella domanda: i bianchi hanno poca biomassa, e la
  loro qualita' puo' essere peggiore senza che i parametri siano sbagliati.

L'esclusione dei negativi e' di principio, non una misura: vale anche quando i
negativi conservano quanto i biologici. Le classi sono in
:data:`~amplicon16s.metadata.models.CLASSI_CONTROLLATE`. La media e' per campione, come la chiede il controllo: un campione povero
pesa quanto uno profondo.
"""

from __future__ import annotations

import gzip
import zlib
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Any, ClassVar, Final

from amplicon16s.config.schema import Config, Qc, thread_effettivi
from amplicon16s.errors.exceptions import errore
from amplicon16s.io_layer.artifacts import ManifestoPasso
from amplicon16s.io_layer.conteggi import leggi_conteggi
from amplicon16s.metadata.models import CLASSI_CONTROLLATE, ClasseCampione
from amplicon16s.rbridge.runner import cartella_r, esegui_script
from amplicon16s.runner.graph import Passo
from amplicon16s.runner.retry import Aggiustamento, senza_modifiche
from amplicon16s.steps.base import PipelineStep, Produzione, StepContext

__all__ = [
    "FiltroLetture",
    "archivio_incompleto",
    "controlla_filtro",
    
]

NOME_SCRIPT: Final = "02_filter.R"

NOME_PREFILTRO: Final = "letture_prefiltro.tsv"
NOME_FILTRATE: Final = "letture_filtrate.tsv"
#: Le letture filtrate di un campione: ``<accession>_filt.fastq.gz``.
SUFFISSO_FILTRATI: Final = "_filt.fastq.gz"


def archivio_incompleto(percorso: Path) -> str | None:
    """Il motivo per cui un archivio non si decomprime per intero, o ``None``.

    Un FASTQ non compresso, riconosciuto dai primi byte come in S0, non ha
    nulla da decomprimere: basta che si legga fino in fondo.
    """
    try:
        with open(percorso, "rb") as file:
            compresso = file.read(2) == b"\x1f\x8b"
        with (gzip.open if compresso else open)(percorso, "rb") as file:
            while file.read(1 << 22):
                pass
    except (OSError, EOFError, zlib.error) as e:
        return f"{type(e).__name__}: {e}"
    return None


def controlla_filtro(
    ingresso: dict[str, int],
    uscita: dict[str, int],
    classi: dict[str, ClasseCampione],
    qc: Qc,
) -> dict[str, Any]:
    """Le metriche del filtro per classe, e i controlli E-S2-01 ed E-S2-02.

    I controlli guardano solo le classi di :data:`CLASSI_CONTROLLATE`; un
    controllo negativo azzerato compare fra le metriche e non ferma nulla.
    """
    azzerati: dict[str, list[str]] = {classe.value: [] for classe in ClasseCampione}
    frazioni: dict[str, list[float]] = {classe.value: [] for classe in ClasseCampione}
    for campione, letture in ingresso.items():
        classe = classi[campione].value
        if uscita.get(campione, 0) == 0:
            azzerati[classe].append(campione)
        if letture > 0:
            frazioni[classe].append(uscita.get(campione, 0) / letture)
    metriche: dict[str, Any] = {
        "campioni": len(ingresso),
        "letture_ingresso": sum(ingresso.values()),
        "letture_uscita": sum(uscita.values()),
        "azzerati": {k: sorted(v) for k, v in azzerati.items()},
        "frazione_conservata_media": {
            k: round(sum(v) / len(v), 4) for k, v in frazioni.items() if v
        },
    }

    controllati = sorted(c for classe in CLASSI_CONTROLLATE for c in azzerati[classe.value])
    if len(controllati) > qc.max_zeroed_samples:
        raise errore(
            "E-S2-01",
            f"{len(controllati)} campioni biologici o controlli positivi senza "
            f"letture dopo il filtro, oltre i {qc.max_zeroed_samples} ammessi da "
            f"qc.max_zeroed_samples: {', '.join(controllati)}",
            campioni=controllati,
        )

    considerate = [f for classe in CLASSI_CONTROLLATE for f in frazioni[classe.value]]
    if considerate:
        perdita = 1 - sum(considerate) / len(considerate)
        metriche["perdita_media_controllata"] = round(perdita, 4)
        if perdita > qc.max_frac_lost_filter:
            raise errore(
                "E-S2-02",
                f"il filtro scarta in media il {perdita:.1%} delle letture dei "
                f"campioni biologici e dei controlli positivi, oltre il "
                f"{qc.max_frac_lost_filter:.0%} di qc.max_frac_lost_filter",
                perdita_media=round(perdita, 4),
            )
    return metriche


class FiltroLetture(PipelineStep):
    """S2: filtro e troncamento."""

    passo: ClassVar[Passo] = Passo.S2
    #: 2: riconosce dai primi byte un FASTQ non compresso, che non ha un
    #: archivio da verificare.
    versione: ClassVar[int] = 2
    script_r: ClassVar[str | None] = NOME_SCRIPT
    passi_tracciamento: ClassVar[tuple[str, ...]] = ("prefiltro", "filtrate")
    #: Tutto il gruppo filter e le due soglie dei controlli sul risultato.
    #: run.batch_size, i file passati insieme al filtro, non vi entra: ogni file
    #: si filtra da solo. Variare filter.maxEE o qc.max_frac_lost_filter
    #: ricalcola S2 e le fasi a valle (S3+), ma preserva interamente S0 e S1 gia'
    #: concluse.
    parametri: ClassVar[tuple[str, ...]] = (
        "filter", "qc.max_zeroed_samples", "qc.max_frac_lost_filter",
    )
    aggiustamenti: ClassVar[dict[str, Aggiustamento]] = {
        "E-S2-03": senza_modifiche(
            "un errore di lettura transitorio si corregge rileggendo lo stesso file"
        ),
    }

    def artefatti_temporanei(self, manifesto: ManifestoPasso, config: Config) -> tuple[str, ...]:
        """Le letture filtrate, se ``run.keep_filtered_fastq`` e' falso."""
        if config.run.keep_filtered_fastq:
            return ()
        return tuple(n for n in manifesto.nomi if n.endswith(SUFFISSO_FILTRATI))

    def calcola(self, contesto: StepContext) -> Produzione:
        """Verifica che gli archivi si decomprimano per intero, esegue ``R/02_filter.R``
        e controlla le letture rimaste per classe.
        """
        if contesto.inventario is None:
            raise RuntimeError("S2 richiede l'inventario prodotto da S0")
        config = contesto.config
        filtro = config.filter
        campioni = {
            c.accession: str(c.file) for c in contesto.inventario if c.file is not None
        }
        with ThreadPoolExecutor(max_workers=thread_effettivi(config)) as esecutore:
            motivi = dict(zip(campioni, esecutore.map(archivio_incompleto, map(Path, campioni.values()))))
        incompleti = {a: m for a, m in motivi.items() if m is not None}
        if incompleti:
            raise errore(
                "E-S2-03",
                "archivi che non si decomprimono per intero: " + "; ".join(
                    f"{a} ({Path(campioni[a]).name}: {m})" for a, m in sorted(incompleti.items())
                ),
                campioni=sorted(incompleti),
            )
        esito = esegui_script(
            cartella_r() / NOME_SCRIPT,
            {
                "campioni": campioni,
                "truncLen": filtro.truncLen,
                "trimLeft": filtro.trimLeft,
                "maxEE": filtro.maxEE,
                "truncQ": filtro.truncQ,
                "maxN": filtro.maxN,
                "rm_phix": filtro.rm_phix,
                "processi": thread_effettivi(config),
                "lotto": config.run.batch_size,
            },
            contesto.albero,
            self.cartella,
            passo=self.passo,
            logger=contesto.logger,
        )
        cartella = contesto.albero.cartella(self.cartella)
        metriche = controlla_filtro(
            leggi_conteggi(cartella / NOME_PREFILTRO),
            leggi_conteggi(cartella / NOME_FILTRATE),
            {c.accession: c.classe for c in contesto.inventario},
            config.qc,
        )
        negativi = metriche["azzerati"][ClasseCampione.CONTROLLO_NEGATIVO.value]
        if negativi:
            contesto.logger.info(
                f"{len(negativi)} controlli negativi senza letture dopo il filtro: bianchi puliti",
                extra={"passo": "S2", "campioni": negativi},
            )
        return Produzione(esito.artefatti, metriche)
