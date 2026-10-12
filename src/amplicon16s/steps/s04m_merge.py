"""Fase SM: la fusione delle coppie, solo in modalita' appaiata.

Dopo l'inferenza ogni campione ha le varianti della lettura forward e quelle
della reverse. La fusione (``R/04m_merge.R``, ``dada2::mergePairs``, attraverso
il ponte) ricompone la regione amplificata: per ogni coppia di letture
attribuite a due varianti allinea la forward con il complemento inverso della
reverse, e la tiene se le due si sovrappongono per almeno ``pair.min_overlap``
basi con non piu' di ``pair.max_mismatch`` discordanze. Scrive in
``05_pair_merge/`` le sequenze fuse di ogni campione, nella forma delle
varianti di S4: da li' la catena e' quella della modalita' single-end.

**Gli ingressi sono i soli file che S4 conserva per la fusione.** Di ogni
campione e di ogni lettura S4 scrive le mappe dalle letture alle sequenze
uniche e da queste alle varianti, e la tabella delle varianti: e' cio' che
``mergePairs`` legge. Le letture filtrate non si rileggono (dal filtro la
fase prende i soli conteggi per campione): cambiare ``pair.min_overlap``
rifa' questa fase e le successive, non il filtro ne' l'inferenza. Se pero' le
letture filtrate sono state tolte a esecuzione conclusa
(``run.keep_filtered_fastq`` falso) vale la regola degli intermedi rimossi:
il filtro si rifa', e con esso l'inferenza.

**La quota di coppie fuse si controlla dove ha senso.** Come la frazione
chimerica, si giudica la quota complessiva dei campioni biologici e dei
controlli positivi (:data:`~amplicon16s.metadata.models.CLASSI_CONTROLLATE`):
le coppie fuse su quelle entrate nella fusione, cioe' con entrambe le letture
attribuite a una variante. Sotto ``pair.min_merged_frac_stop`` la fase si
ferma (``E-SM-01``): troncamenti troppo corti per la regione, o file che non
corrispondono, non lasciano una sovrapposizione, e proseguire consegnerebbe i
soli campioni fortunati. Sotto ``pair.min_merged_frac_warn`` prosegue e lo
registra (``E-SM-02``). I confronti sono stretti.

**La lunghezza delle sequenze fuse si confronta con quella attesa.** La
distribuzione e' in ``lunghezze_fuse.tsv``. Se ``pair.expected_merged_len``
non e' nullo e la mediana, pesata sulle coppie, se ne scosta di piu' di
:data:`SCARTO_DELLA_LUNGHEZZA_ATTESA`, la fase lo registra (``E-SM-03``).
"""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any, ClassVar, Final

from amplicon16s.config.schema import Pair
from amplicon16s.errors.exceptions import errore
from amplicon16s.io_layer.artifacts import Fase
from amplicon16s.io_layer.conteggi import leggi_conteggi
from amplicon16s.metadata.models import CLASSI_CONTROLLATE
from amplicon16s.rbridge.runner import cartella_r, esegui_script
from amplicon16s.runner.graph import Passo
from amplicon16s.steps.base import PipelineStep, Produzione, StepContext
from amplicon16s.steps.s02_filter import NOME_FILTRATE
from amplicon16s.steps.s04_dada import MODELLO_DEL_NOME_PER_LA_FUSIONE

__all__ = [
    "FusioneCoppie",
    "NOME_FUSE",
    "NOME_FUSIONE",
    "NOME_LUNGHEZZE_FUSE",
    "NOME_RIEPILOGO",
    "NOME_VARIANTI_FUSE",
    "SCARTO_DELLA_LUNGHEZZA_ATTESA",
    "controlla_fusione",
    "mediana_pesata",
    "misura_fusione",
]

NOME_SCRIPT: Final = "04m_merge.R"
#: Le sequenze fuse di ogni campione, ingresso di S5 in modalita' appaiata.
NOME_VARIANTI_FUSE: Final = "varianti_per_campione.rds"
NOME_FUSIONE: Final = "fusione.tsv"
NOME_LUNGHEZZE_FUSE: Final = "lunghezze_fuse.tsv"
NOME_RIEPILOGO: Final = "fusione.json"
#: Le coppie fuse per campione, nel formato del tracciamento.
NOME_FUSE: Final = "letture_fuse.tsv"

#: Di quanto la mediana delle lunghezze fuse puo' scostarsi da quella attesa,
#: come frazione dell'attesa, prima dell'avviso. Le regioni ipervariabili del
#: 16S variano di lunghezza fra i taxa di pochi punti percentuali: uno scarto
#: maggiore indica un altro amplicone, o fusioni spurie.
SCARTO_DELLA_LUNGHEZZA_ATTESA: Final = 0.05


def _leggi_fusione(percorso: Path) -> dict[str, dict[str, int]]:
    """``fusione.tsv``: campione -> coppie entrate e coppie fuse."""
    with open(percorso, encoding="utf-8", newline="") as file:
        return {
            r["campione"]: {"coppie": int(float(r["coppie"])), "fuse": int(float(r["fuse"]))}
            for r in csv.DictReader(file, delimiter="\t")
        }


def mediana_pesata(per_lunghezza: dict[int, int]) -> int | None:
    """La lunghezza mediana, pesata sulle coppie: la piu' corta entro cui cade
    almeno meta' delle coppie. ``None`` senza coppie.
    """
    totale = sum(per_lunghezza.values())
    if totale == 0:
        return None
    cumulata = 0
    for lunghezza in sorted(per_lunghezza):
        cumulata += per_lunghezza[lunghezza]
        if 2 * cumulata >= totale:
            return lunghezza
    return None


def _frazione(parte: int, totale: int) -> float | None:
    """Il rapporto arrotondato a quattro decimali, o ``None`` se il totale e' zero."""
    return round(parte / totale, 4) if totale else None


def misura_fusione(
    per_campione: dict[str, dict[str, int]], classi: dict[str, str]
) -> dict[str, Any]:
    """Le quote di coppie fuse per classe, e quella delle classi controllate."""
    per_classe: dict[str, dict[str, int]] = {}
    for campione, conteggi in per_campione.items():
        voce = per_classe.setdefault(classi[campione], {"coppie": 0, "fuse": 0, "campioni": 0})
        voce["coppie"] += conteggi["coppie"]
        voce["fuse"] += conteggi["fuse"]
        voce["campioni"] += 1
    controllate = [per_classe[c.value] for c in CLASSI_CONTROLLATE if c.value in per_classe]
    coppie = sum(c["coppie"] for c in controllate)
    fuse = sum(c["fuse"] for c in controllate)
    return {
        "per_classe": {
            classe: {**voce, "quota": _frazione(voce["fuse"], voce["coppie"])}
            for classe, voce in sorted(per_classe.items())
        },
        "classi_controllate": [c.value for c in CLASSI_CONTROLLATE],
        "coppie_controllate": coppie,
        "fuse_controllate": fuse,
        "quota_controllata": _frazione(fuse, coppie),
    }


def controlla_fusione(misure: dict[str, Any], coppie: Pair) -> str | None:
    """I controlli E-SM-01 ed E-SM-02 sulla quota controllata.

    Solleva E-SM-01 sotto ``pair.min_merged_frac_stop``; sotto
    ``pair.min_merged_frac_warn`` restituisce il dettaglio della degradazione
    da registrare, altrimenti ``None``. Senza coppie da giudicare non c'e'
    una quota, e nessun controllo scatta: lo dira' la tabella vuota.
    """
    if misure["quota_controllata"] is None:
        return None
    # Il confronto e' sul rapporto esatto: quello arrotondato del riepilogo
    # potrebbe raggiungere la soglia da sotto.
    quota = misure["fuse_controllate"] / misure["coppie_controllate"]
    testo = (
        f"il {quota:.1%} delle coppie dei campioni biologici e dei controlli positivi si "
        f"fonde ({misure['fuse_controllate']} su {misure['coppie_controllate']})"
    )
    if quota < coppie.min_merged_frac_stop:
        raise errore(
            "E-SM-01",
            f"{testo}, sotto il {coppie.min_merged_frac_stop:.0%} di "
            "pair.min_merged_frac_stop",
            quota=quota,
        )
    if quota < coppie.min_merged_frac_warn:
        return (f"{testo}, sotto il {coppie.min_merged_frac_warn:.0%} di "
                "pair.min_merged_frac_warn")
    return None


def campioni_sotto_la_quota(
    per_campione: dict[str, dict[str, int]], classi: dict[str, str], minima: float | None
) -> list[tuple[str, int, int]]:
    """I campioni biologici e i controlli positivi la cui quota di coppie fuse
    sta sotto ``minima``: (campione, coppie fuse, coppie), in ordine di nome.

    La quota dell'insieme non vede un campione solo che fonde pochissime
    coppie. Senza soglia (``pair.min_merged_frac_sample_warn`` omessa) non c'e'
    alcun controllo; un campione senza coppie non ha una quota. Il confronto e'
    stretto, sul rapporto esatto.
    """
    if minima is None:
        return []
    controllate = {c.value for c in CLASSI_CONTROLLATE}
    return [
        (campione, conteggi["fuse"], conteggi["coppie"])
        for campione, conteggi in sorted(per_campione.items())
        if classi[campione] in controllate and conteggi["coppie"]
        and conteggi["fuse"] / conteggi["coppie"] < minima
    ]


class FusioneCoppie(PipelineStep):
    """SM: le sequenze fuse di ogni campione."""

    passo: ClassVar[Passo] = Passo.SM
    versione: ClassVar[int] = 1
    script_r: ClassVar[str | None] = NOME_SCRIPT
    #: La fase esiste solo in modalita' appaiata, e i suoi parametri sono tutti
    #: del gruppo pair: sovrapposizione minima, discordanze ammesse, code oltre
    #: l'altra lettura, le due soglie sulla quota e la lunghezza attesa.
    parametri: ClassVar[tuple[str, ...]] = (
        "pair.min_overlap", "pair.max_mismatch", "pair.trim_overhang",
        "pair.min_merged_frac_warn", "pair.min_merged_frac_stop",
        "pair.expected_merged_len",
        # Facoltativa: omessa non entra nell'impronta, e nessun campione e' giudicato.
        "pair.min_merged_frac_sample_warn",
    )
    passi_tracciamento: ClassVar[tuple[str, ...]] = ("fuse",)

    def calcola(self, contesto: StepContext) -> Produzione:
        """Esegue ``R/04m_merge.R`` sui file che S4 ha conservato, scrive il
        riepilogo e applica i controlli sulla quota e sulla lunghezza.
        """
        if contesto.inventario is None:
            raise RuntimeError("SM richiede l'inventario prodotto da S0")
        coppie = contesto.config.pair
        inferenza = contesto.albero.cartella(Fase.ASV_INFERENCE)
        filtrate = leggi_conteggi(contesto.albero.cartella(Fase.FILTERED) / NOME_FILTRATE)
        accession = sorted(c.accession for c in contesto.inventario)
        con_coppie = [a for a in accession if filtrate.get(a, 0) > 0]
        esito = esegui_script(
            cartella_r() / NOME_SCRIPT,
            {
                lettura: {
                    a: str(inferenza / MODELLO_DEL_NOME_PER_LA_FUSIONE.format(
                        lettura=lettura, accession=a))
                    for a in con_coppie
                } for lettura in ("forward", "reverse")
            } | {
                "senza_coppie": [a for a in accession if a not in con_coppie],
                "min_overlap": coppie.min_overlap,
                "max_mismatch": coppie.max_mismatch,
                "trim_overhang": coppie.trim_overhang,
            },
            contesto.albero,
            self.cartella,
            passo=self.passo,
            tempo_massimo_s=contesto.config.run.r_timeout_s,
            logger=contesto.logger,
        )
        cartella = contesto.albero.cartella(self.cartella)
        per_campione = _leggi_fusione(cartella / NOME_FUSIONE)
        classi = {c.accession: c.classe.value for c in contesto.inventario}
        misure = misura_fusione(per_campione, classi)
        with open(cartella / NOME_LUNGHEZZE_FUSE, encoding="utf-8", newline="") as file:
            per_lunghezza = {
                int(r["lunghezza"]): int(float(r["coppie"]))
                for r in csv.DictReader(file, delimiter="\t")
            }
        mediana = mediana_pesata(per_lunghezza)
        attesa = coppie.expected_merged_len
        # Il riepilogo si scrive prima dei controlli: un arresto lascia le
        # misure che l'hanno causato.
        riepilogo = {
            **misure,
            "min_overlap": coppie.min_overlap,
            "max_mismatch": coppie.max_mismatch,
            "trim_overhang": coppie.trim_overhang,
            "lunghezza_minima": min(per_lunghezza) if per_lunghezza else None,
            "lunghezza_massima": max(per_lunghezza) if per_lunghezza else None,
            "lunghezza_mediana": mediana,
            "lunghezza_attesa": attesa,
        }
        artefatto = contesto.albero.scrivi_testo(
            self.cartella, NOME_RIEPILOGO,
            json.dumps(riepilogo, indent=2, ensure_ascii=False, sort_keys=True) + "\n",
        )
        avviso = controlla_fusione(misure, coppie)
        if avviso is not None:
            contesto.degrada("E-SM-02", avviso, quota=misure["quota_controllata"])
        sotto = campioni_sotto_la_quota(per_campione, classi, coppie.min_merged_frac_sample_warn)
        if sotto:
            contesto.degrada(
                "E-SM-04",
                "{} campioni fondono meno del {:.0%} delle coppie "
                "(pair.min_merged_frac_sample_warn): {}".format(
                    len(sotto), coppie.min_merged_frac_sample_warn,
                    ", ".join(f"{campione} ({fuse / totali:.1%}, {fuse} su {totali})"
                              for campione, fuse, totali in sotto)),
                campioni=[campione for campione, _, _ in sotto],
            )
        if (attesa is not None and mediana is not None
                and abs(mediana - attesa) > SCARTO_DELLA_LUNGHEZZA_ATTESA * attesa):
            contesto.degrada(
                "E-SM-03",
                f"la lunghezza mediana delle sequenze fuse e' {mediana} basi, e "
                f"pair.expected_merged_len ne dichiara {attesa}: uno scarto oltre il "
                f"{SCARTO_DELLA_LUNGHEZZA_ATTESA:.0%} dell'attesa",
                mediana=mediana, attesa=attesa,
            )
        metriche = {
            "campioni": len(con_coppie),
            "coppie": sum(v["coppie"] for v in per_campione.values()),
            "fuse": sum(v["fuse"] for v in per_campione.values()),
            "quota_controllata": misure["quota_controllata"],
            "lunghezza_mediana": mediana,
        }
        return Produzione(esito.artefatti + (artefatto,), metriche)
