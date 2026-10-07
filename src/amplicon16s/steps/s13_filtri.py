"""Fase S13: i filtri finali, per profondita', tassonomici e di prevalenza.

Dall'oggetto decontaminato di S12 ricava l'oggetto finale: **solo i campioni
biologici**, perche' e' il risultato destinato alle analisi ecologiche, e un
controllo al suo interno rischierebbe di essere trattato come un campione
ambientale. I controlli positivi e negativi non si perdono: S14 li consegna
in ``ps_controlli.rds``, dall'oggetto integrato di S10. Scrive in
``12_final/intermedi/``: l'oggetto filtrato e le tabelle dei filtri sono
intermedi di calcolo, e in ``12_final/`` stanno solo i file consegnati da
S14. Il calcolo dei filtri sulle varianti e' in ``R/13_filtri.R``, che
dichiara l'ordine dei filtri e il denominatore della prevalenza.

**Il filtro per profondita'**, che il piano di S13 non nomina, viene per primo
ed e' deciso qui. S11 scrive in ``soglia.json`` la soglia di ogni piastra con
lo stadio delle letture a cui si applica: ``nonchimeric``, le letture senza
chimere dopo il filtro di lunghezza (``letture_lunghezza.tsv`` di S7), o
``raw``, le letture grezze (``letture_prefiltro.tsv`` di S2). Ogni campione
biologico si confronta con la soglia della sua piastra usando la sua
profondita' a quello stadio, letta dal tracciamento: non la profondita'
dell'oggetto, che dopo la decontaminazione e' piu' bassa e non e' la
grandezza con cui la soglia e' stata stimata.

**Gli identificativi delle varianti** non si rinumerano: i vuoti lasciati dai
filtri sono attesi.

**I campioni esclusi**, con il filtro e il motivo di ciascuno, sono in
``esclusioni.tsv``; le varianti rimosse in ``varianti_rimosse.tsv``. Un
campione svuotato dal filtro tassonomico non ha segnale batterico: esce, e la
fase registra ``E-S13-03`` e prosegue. Anche uno svuotato dal filtro di
prevalenza esce con il motivo, e la fase registra ``E-S13-02`` come
degradazione dichiarata: aveva letture batteriche, solo rare, e chi legge il
manifesto deve sapere che la soglia di prevalenza lo ha tolto. Se nessun
campione biologico supera i filtri la fase si ferma con ``E-S13-04``: un
oggetto finale vuoto non e' un risultato.

**I nomi dei taxa** si confrontano senza il prefisso di rango (``p__``,
``o__``) che alcuni riferimenti portano. Se la tassonomia non ha il rango
Phylum, ``filt.remove_na_phylum`` non ha su che cosa operare: la fase lo
dichiara con ``E-S13-05`` e applica i soli ``filt.exclude_taxa``.
"""

from __future__ import annotations

import json
from typing import Any, ClassVar, Final

from amplicon16s.io_layer.artifacts import Fase
from amplicon16s.io_layer.conteggi import leggi_conteggi
from amplicon16s.metadata.models import ClasseCampione
from amplicon16s.rbridge.runner import cartella_r, esegui_script
from amplicon16s.runner.graph import Passo
from amplicon16s.steps.base import PipelineStep, Produzione, StepContext
from amplicon16s.steps.s02_filter import NOME_PREFILTRO
from amplicon16s.steps.s11_controls import NOME_SOGLIA
from amplicon16s.steps.s12_decontam import NOME_OGGETTO_DECONTAMINATO

__all__ = ["FiltriFinali", "NOME_FILTRATO", "NOME_RIEPILOGO", "esclusi_per_profondita"]

NOME_SCRIPT: Final = "13_filtri.R"
#: L'oggetto filtrato che S14 serializza.
NOME_FILTRATO: Final = "ps_filtrato.rds"
NOME_RIEPILOGO: Final = "filtri_riepilogo.json"
#: Le letture senza chimere dopo il filtro di lunghezza, scritte da S7.
NOME_LUNGHEZZA: Final = "letture_lunghezza.tsv"
#: Il riepilogo di S8, con i ranghi della tassonomia assegnata.
NOME_RIEPILOGO_TASSONOMIA: Final = "riepilogo.json"


def esclusi_per_profondita(
    campioni: list[tuple[str, str | None]],
    soglia: dict[str, Any],
    letture: dict[str, dict[str, int]],
) -> tuple[list[str], list[dict[str, str]]]:
    """I campioni sopra e sotto la soglia di profondita' della loro piastra.

    ``campioni`` sono coppie (accession, piastra); ``soglia`` e' il documento
    di S11; ``letture`` le letture per campione di ciascuno stadio. Un campione
    senza piastra usa la soglia ``senza_piastra``. Restituisce i tenuti, in
    ordine, e le esclusioni con il motivo.
    """
    per_piastra = soglia["per_piastra"] or {}
    tenuti: list[str] = []
    esclusi: list[dict[str, str]] = []
    for accession, piastra in campioni:
        voce = per_piastra.get(piastra) if piastra is not None else None
        if voce is None:
            voce = soglia["senza_piastra"]
        stadio = voce["stadio"]
        if stadio not in letture:
            raise RuntimeError(f"stadio della soglia sconosciuto: {stadio!r}")
        profondita = letture[stadio].get(accession, 0)
        if profondita < voce["valore"]:
            esclusi.append({
                "accession": accession,
                "motivo": (
                    f"{profondita} letture allo stadio {stadio}, meno della soglia "
                    f"{voce['valore']} della piastra {piastra or 'non nota'} "
                    f"({voce['origine']})"
                ),
            })
        else:
            tenuti.append(accession)
    return tenuti, esclusi


class FiltriFinali(PipelineStep):
    """S13: l'oggetto dei soli biologici, filtrato per profondita', taxa e prevalenza."""

    passo: ClassVar[Passo] = Passo.S13
    #: 2: scrive in 12_final/intermedi/ e non produce piu' ps_controlli.rds
    #: (lo consegna S14); confronta i taxa senza il prefisso di rango; un
    #: campione svuotato dalla prevalenza esce invece di fermare la fase.
    #: 3: una tassonomia di un solo rango non ferma il filtro sui taxa.
    versione: ClassVar[int] = 3
    script_r: ClassVar[str | None] = NOME_SCRIPT
    passi_tracciamento: ClassVar[tuple[str, ...]] = ("finali",)
    #: I filtri tassonomici (filt), quello di prevalenza (prev) e le letture
    #: finali minime per campione.
    parametri: ClassVar[tuple[str, ...]] = ("filt", "prev", "qc.min_reads_final")

    def calcola(self, contesto: StepContext) -> Produzione:
        """Applica il filtro per profondita' sul tracciamento, poi esegue
        ``R/13_filtri.R`` per gli altri filtri e l'oggetto finale.
        """
        if contesto.inventario is None:
            raise RuntimeError("S13 richiede l'inventario prodotto da S0")
        config = contesto.config
        albero = contesto.albero
        soglia = json.loads(
            (albero.cartella(Fase.CONTROLS) / NOME_SOGLIA).read_text(encoding="utf-8")
        )
        letture = {
            "nonchimeric": leggi_conteggi(albero.cartella(Fase.CHIMERA) / NOME_LUNGHEZZA),
            "raw": leggi_conteggi(albero.cartella(Fase.FILTERED) / NOME_PREFILTRO),
        }
        biologici = [(c.accession, c.piastra) for c in contesto.inventario
                     if c.classe is ClasseCampione.BIOLOGICO]
        tenuti, esclusi = esclusi_per_profondita(biologici, soglia, letture)

        ranghi = json.loads(
            (albero.cartella(Fase.TAXONOMY) / NOME_RIEPILOGO_TASSONOMIA).read_text(encoding="utf-8")
        )["ranghi"]
        if config.filt.remove_na_phylum and "Phylum" not in ranghi:
            contesto.degrada(
                "E-S13-05",
                f"la tassonomia ha i ranghi {', '.join(ranghi)}, senza Phylum: "
                "filt.remove_na_phylum non applicato",
            )

        esito = esegui_script(
            cartella_r() / NOME_SCRIPT,
            {
                "decontaminato": str(albero.cartella(Fase.CONTROLS) / NOME_OGGETTO_DECONTAMINATO),
                "tenuti": tenuti,
                "esclusi_profondita": esclusi,
                "senza_phylum": config.filt.remove_na_phylum,
                "taxa_esclusi": list(config.filt.exclude_taxa),
                "prevalenza": config.prev.apply,
                "frazione": config.prev.min_fraction,
                "minimo_conteggio": config.prev.min_count,
                "letture_minime": config.qc.min_reads_final,
                "letture_nonchimeric": {a: letture["nonchimeric"].get(a, 0) for a in tenuti},
            },
            albero,
            self.cartella,
            passo=self.passo,
            logger=contesto.logger,
        )
        riepilogo = json.loads(
            (albero.cartella(self.cartella) / NOME_RIEPILOGO).read_text(encoding="utf-8")
        )
        svuotati = riepilogo["campioni_svuotati"]
        if svuotati["tassonomico"]:
            contesto.degrada(
                "E-S13-03",
                f"{len(svuotati['tassonomico'])} campioni senza letture dopo il filtro "
                f"tassonomico: {', '.join(svuotati['tassonomico'])}",
                campioni=svuotati["tassonomico"],
            )
        if svuotati["prevalenza"]:
            contesto.degrada(
                "E-S13-02",
                f"{len(svuotati['prevalenza'])} campioni senza letture dopo il filtro di "
                f"prevalenza: {', '.join(svuotati['prevalenza'])}",
                campioni=svuotati["prevalenza"],
            )
        metriche = {
            "campioni": riepilogo["campioni"],
            "varianti": {k: riepilogo["varianti"][k] for k in ("iniziali", "rimosse", "finali")},
            "prevalenza": riepilogo["prevalenza"],
        }
        return Produzione(esito.artefatti, metriche)
