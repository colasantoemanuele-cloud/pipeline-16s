"""Fase S8: l'assegnazione tassonomica delle varianti.

Assegna a ogni variante della tabella di S7 una linea tassonomica con il
classificatore bayesiano naive di dada2 (``assignTaxonomy``, ``R/08_taxonomy.R``)
sul riferimento ``tax.ref_fasta``, e scrive in ``08_taxonomy/`` la tabella
tassonomica, il bootstrap di ciascun rango (servira' ai filtri di S13 e
all'interpretazione) e la copertura del phylum per classe di campioni.

**Il classificatore.** ``tax.classifier`` vale ``naive_bayes``: il training set
di SILVA per IdTaxa non e' distribuito da alcuna fonte. ``tax.min_boot`` e' 50
e non 80, e ``tax.assign_species`` e' falso: le letture coprono 137 delle circa
253 basi dell'amplicone V4, e su sequenze piu' corte il bootstrap e'
sistematicamente piu' basso a parita' di correttezza.

**Riproducibilita'.** Il seme viene da ``run.seed``. ``run.threads`` resta fuori
dall'impronta: il bootstrap estrae i k-meri con ``runif`` di R prima del calcolo
parallelo, e il numero di thread non cambia i byte degli artefatti (verificato
sul dataset completo e nei test).

**Il difetto noto del riferimento.** Se ``tax.ref_bad_taxa`` indica l'elenco dei
taxa con un difetto noto (per SILVA 138 versione 2, un rango mancante nel
percorso di 10 famiglie e 114 generi, per cui il nome compare nella colonna del
rango superiore), il riferimento non si modifica: ``difetto_riferimento.tsv``
elenca le assegnazioni che ricadono su quei taxa, in qualunque colonna, e il
riepilogo riporta quante varianti e quante letture ne sono interessate.

**La copertura, E-S8-02.** La frazione di varianti con il phylum assegnato si
misura per classe di campioni, sulle varianti presenti nei campioni della
classe; ``qc.min_frac_phylum`` si applica alle classi di :data:`CLASSI_CONTROLLATE`.
"""

from __future__ import annotations

import csv
import hashlib
import json
from collections.abc import Mapping
from pathlib import Path
from typing import Any, ClassVar, Final

from amplicon16s.config.schema import Config, Qc
from amplicon16s.errors.exceptions import errore
from amplicon16s.io_layer.artifacts import Fase
from amplicon16s.metadata.models import ClasseCampione
from amplicon16s.rbridge.runner import cartella_r, esegui_script
from amplicon16s.runner.graph import Passo
from amplicon16s.steps.base import PipelineStep, Produzione, StepContext
from amplicon16s.steps.s07_asv_length import NOME_TABELLA_ASV

__all__ = [
    "AssegnazioneTassonomica",
    "CLASSI_CONTROLLATE",
    "NOME_TASSONOMIA",
    "controlla_copertura",
    "leggi_taxa_difettosi",
    "marca_difetti",
]

NOME_SCRIPT: Final = "08_taxonomy.R"
#: La tabella tassonomica che le fasi successive leggono.
NOME_TASSONOMIA: Final = "tassonomia.tsv"
NOME_COPERTURA: Final = "copertura_per_gruppo.tsv"
NOME_DIFETTI: Final = "difetto_riferimento.tsv"
NOME_RIEPILOGO: Final = "riepilogo.json"
#: Le classi a cui si applica qc.min_frac_phylum (E-S8-02).
CLASSI_CONTROLLATE: Final = (ClasseCampione.BIOLOGICO, ClasseCampione.CONTROLLO_POSITIVO)


def leggi_taxa_difettosi(percorso: Path) -> dict[str, str]:
    """I taxa con un difetto noto del riferimento: nome -> rango atteso."""
    with open(percorso, encoding="utf-8", newline="") as file:
        return {r["name"]: r["rank"] for r in csv.DictReader(file)}


def marca_difetti(
    righe: list[dict[str, str]], ranghi: list[str], difettosi: Mapping[str, str]
) -> list[dict[str, Any]]:
    """Le assegnazioni che ricadono su un taxon difettoso, una per colonna.

    Il nome si cerca in ogni colonna, non in quella del suo rango: il difetto
    e' proprio un rango mancante, che fa comparire il nome in una colonna piu'
    alta di quella attesa.
    """
    marcate = []
    for riga in righe:
        for colonna in ranghi:
            nome = riga[colonna]
            if nome and nome in difettosi:
                marcate.append({
                    "sequenza": riga["sequenza"], "colonna": colonna, "taxon": nome,
                    "rango": difettosi[nome], "letture": int(float(riga["letture"])),
                })
    return marcate


def controlla_copertura(frazioni: Mapping[str, Mapping[str, float | None]], qc: Qc) -> None:
    """E-S8-02 se una classe controllata ha troppe varianti senza phylum."""
    sotto = {
        c.value: frazioni[c.value]["varianti"] for c in CLASSI_CONTROLLATE
        if c.value in frazioni and frazioni[c.value]["varianti"] is not None
        and frazioni[c.value]["varianti"] < qc.min_frac_phylum
    }
    if sotto:
        raise errore(
            "E-S8-02",
            "frazione di varianti con il phylum sotto qc.min_frac_phylum "
            f"({qc.min_frac_phylum}): " + ", ".join(f"{c} {f:.4f}" for c, f in sorted(sotto.items())),
            frazioni=dict(sotto),
        )


class AssegnazioneTassonomica(PipelineStep):
    """S8: la tassonomia delle varianti, con il bootstrap di ogni rango."""

    passo: ClassVar[Passo] = Passo.S8
    versione: ClassVar[int] = 1
    script_r: ClassVar[str | None] = NOME_SCRIPT
    #: Il riferimento e il classificatore (il contenuto del riferimento entra
    #: con tax.ref_md5, che G12 verifica a ogni avvio); il seme del bootstrap;
    #: la soglia di copertura. run.threads no: non cambia i risultati.
    parametri: ClassVar[tuple[str, ...]] = ("tax", "run.seed", "qc.min_frac_phylum")

    def impronta_dati_esterni(self, config: Config) -> str | None:
        """L'impronta del contenuto dell'elenco dei taxa difettosi, se indicato."""
        if config.tax.ref_bad_taxa is None:
            return None
        contenuto = Path(config.tax.ref_bad_taxa).read_bytes()
        return "sha256:" + hashlib.sha256(contenuto).hexdigest()

    def calcola(self, contesto: StepContext) -> Produzione:
        """Esegue ``R/08_taxonomy.R`` sulla tabella di S7, marca i taxa difettosi e
        controlla la copertura del phylum.
        """
        if contesto.inventario is None:
            raise RuntimeError("S8 richiede l'inventario prodotto da S0")
        config = contesto.config
        tax = config.tax
        if tax.assign_species:
            raise RuntimeError(
                "tax.assign_species vero non e' realizzato: S8 assegna fino al genere"
            )
        esito = esegui_script(
            cartella_r() / NOME_SCRIPT,
            {
                "tabella": str(contesto.albero.cartella(Fase.CHIMERA) / NOME_TABELLA_ASV),
                "riferimento": str(tax.ref_fasta),
                "min_boot": tax.min_boot,
                "try_rc": tax.try_rc,
                "seme": config.run.seed,
                "processi": config.run.threads,
                "gruppi": {c.accession: c.classe.value for c in contesto.inventario},
            },
            contesto.albero,
            self.cartella,
            passo=self.passo,
            logger=contesto.logger,
        )
        cartella = contesto.albero.cartella(self.cartella)
        with open(cartella / NOME_TASSONOMIA, encoding="utf-8", newline="") as file:
            lettore = csv.DictReader(file, delimiter="\t")
            ranghi = [c for c in lettore.fieldnames or [] if c not in ("sequenza", "letture")]
            righe = list(lettore)

        artefatti = list(esito.artefatti)
        difetto: dict[str, Any] = {"elenco": None, "varianti": 0, "letture": 0, "taxa": {}}
        if tax.ref_bad_taxa is not None:
            marcate = marca_difetti(righe, ranghi, leggi_taxa_difettosi(Path(tax.ref_bad_taxa)))
            artefatti.append(contesto.albero.scrivi_testo(
                self.cartella, NOME_DIFETTI,
                "sequenza\tcolonna\ttaxon\trango\tletture\n" + "".join(
                    f"{m['sequenza']}\t{m['colonna']}\t{m['taxon']}\t{m['rango']}\t{m['letture']}\n"
                    for m in marcate
                ),
            ))
            per_variante = {m["sequenza"]: m["letture"] for m in marcate}
            per_taxon: dict[str, dict[str, int]] = {}
            for m in marcate:
                voce = per_taxon.setdefault(m["taxon"], {"varianti": 0, "letture": 0})
                voce["varianti"] += 1
                voce["letture"] += m["letture"]
            difetto = {
                "elenco": Path(tax.ref_bad_taxa).name,
                "varianti": len(per_variante),
                "letture": sum(per_variante.values()),
                "taxa": dict(sorted(per_taxon.items())),
            }

        frazioni: dict[str, dict[str, float | None]] = {}
        per_gruppo: dict[str, dict[str, int]] = {}
        with open(cartella / NOME_COPERTURA, encoding="utf-8", newline="") as file:
            for r in csv.DictReader(file, delimiter="\t"):
                conti = {k: int(float(v)) for k, v in r.items() if k != "gruppo"}
                per_gruppo[r["gruppo"]] = conti
                frazioni[r["gruppo"]] = {
                    "varianti": _frazione(conti["varianti_con_phylum"], conti["varianti"]),
                    "letture": _frazione(conti["letture_con_phylum"], conti["letture"]),
                }

        # Il riepilogo si scrive prima del controllo: un arresto lascia le
        # misure che l'hanno causato.
        riepilogo = {
            "riferimento": {"nome": tax.ref_name, "versione": tax.ref_version, "md5": tax.ref_md5},
            "classificatore": tax.classifier, "min_boot": tax.min_boot, "try_rc": tax.try_rc,
            "seme": config.run.seed,
            "ranghi": ranghi,
            "varianti": len(righe),
            "varianti_con_phylum": sum(1 for r in righe if r["Phylum"]),
            "per_gruppo": per_gruppo,
            "frazione_con_phylum": frazioni,
            "classi_controllate": [c.value for c in CLASSI_CONTROLLATE],
            "difetto_riferimento": difetto,
        }
        artefatti.append(contesto.albero.scrivi_testo(
            self.cartella, NOME_RIEPILOGO,
            json.dumps(riepilogo, indent=2, ensure_ascii=False, sort_keys=True) + "\n",
        ))
        controlla_copertura(frazioni, config.qc)

        metriche = {
            "varianti": riepilogo["varianti"],
            "varianti_con_phylum": riepilogo["varianti_con_phylum"],
            "frazione_con_phylum": frazioni,
            "difetto_riferimento": {k: difetto[k] for k in ("varianti", "letture")},
        }
        return Produzione(tuple(artefatti), metriche)


def _frazione(parte: int, totale: int) -> float | None:
    """Il rapporto arrotondato a quattro decimali, o ``None`` se il totale e' zero."""
    return round(parte / totale, 4) if totale else None
