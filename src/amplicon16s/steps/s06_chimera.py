"""Fase S6: la rimozione delle chimere.

Toglie dalla tabella di S5 le varianti che risultano chimere di due varianti
piu' abbondanti (``R/06_chimera.R``, ``dada2::removeBimeraDenovo``, con i
parametri del gruppo ``chimera``), in ``07_chimera/``.

**La tabella prima delle chimere non si duplica.** E' ``tabella.rds`` di S5,
in ``06_seqtab/``, e resta li': ricopiarla in ``07_chimera/`` raddoppierebbe
su disco l'artefatto piu' grande della pipeline senza aggiungere nulla. S6
registra invece in ``chimere.json`` il riferimento a quel file con il suo
checksum, e in ``chimere.tsv`` l'elenco esplicito delle varianti chimeriche
con le loro letture: con il metodo ``consensus`` la tabella senza chimere e'
quella di S5 meno quelle colonne. La tabella di S5 non puo' cambiare sotto
S6 senza che se ne accorga la valutazione: la sua impronta e' fra quelle da
cui S6 dipende, e il checksum sta nel manifesto di S5.

**La frazione chimerica si misura sulle letture.** Una chimera abbondante
pesa sui risultati piu' di cento chimere da una lettura, e contare le varianti
le metterebbe alla pari. La fase riporta per ogni classe entrambe le misure,
letture e varianti; i controlli guardano le letture:

* E-S6-01, frazione oltre ``qc.stop_frac_chimeric``: arresto;
* E-S6-02, frazione oltre ``qc.warn_frac_chimeric``: si prosegue, e la
  degradazione si registra nel manifesto.

La frazione e' quella complessiva delle letture delle classi controllate,
come la misura ``sum(senza chimere) / sum(tabella)`` del flusso di dada2:
ogni lettura pesa uguale, e un campione con cento letture non sposta la
misura quanto uno con centomila. Come per G10 e per S2, il controllo si
applica alle classi per cui ha senso
(:data:`~amplicon16s.metadata.models.CLASSI_CONTROLLATE`): la domanda
e' se la PCR dei campioni analizzati ha prodotto troppe chimere. I biologici
sono i campioni analizzati; i positivi, una comunita' nota, misurano la
formazione di chimere nella stessa PCR. I negativi ne restano fuori: hanno
poca biomassa, spesso poche letture, e la loro frazione dice poco sulla PCR
dei campioni.

L'esclusione e' di principio, non una misura. Sulle varianti la frazione
chimerica e' di norma molto maggiore che sulle letture: le chimere sono tante
varianti rare, ed e' il motivo per cui la misura che conta e' sulle letture.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any, ClassVar, Final

from amplicon16s.config.schema import Qc, thread_effettivi
from amplicon16s.errors.exceptions import errore
from amplicon16s.io_layer.artifacts import Fase
from amplicon16s.io_layer.conteggi import leggi_conteggi
from amplicon16s.metadata.models import CLASSI_CONTROLLATE
from amplicon16s.rbridge.runner import cartella_r, esegui_script
from amplicon16s.runner.graph import Passo
from amplicon16s.steps.base import PipelineStep, Produzione, StepContext
from amplicon16s.steps.s02_filter import NOME_FILTRATE
from amplicon16s.steps.s05_seqtab import NOME_TABELLA

__all__ = [
    "NOME_SENZA_CHIMERE",
    "RimozioneChimere",
    "controlla_chimere",
    "leggi_per_gruppo",
    "misura_chimere",
]

NOME_SCRIPT: Final = "06_chimera.R"
NOME_SENZA_CHIMERE: Final = "tabella_senza_chimere.rds"
NOME_PER_GRUPPO: Final = "chimere_per_gruppo.tsv"
NOME_RIEPILOGO: Final = "chimere.json"



def leggi_per_gruppo(percorso: Path) -> dict[str, dict[str, int]]:
    """``chimere_per_gruppo.tsv``: gruppo -> conteggi."""
    with open(percorso, encoding="utf-8", newline="") as file:
        return {
            r.pop("gruppo"): {k: int(float(v)) for k, v in r.items()}
            for r in csv.DictReader(file, delimiter="\t")
        }


def _frazione(parte: int, totale: int) -> float | None:
    """Il rapporto arrotondato a quattro decimali, o ``None`` se il totale è zero."""
    return round(parte / totale, 4) if totale else None


def misura_chimere(per_classe: dict[str, dict[str, int]]) -> dict[str, Any]:
    """Le frazioni chimeriche per classe, su letture e varianti, e quella controllata."""
    controllate = [per_classe[c.value] for c in CLASSI_CONTROLLATE if c.value in per_classe]
    chimeriche = sum(c["letture_chimeriche"] for c in controllate)
    letture = sum(c["letture"] for c in controllate)
    return {
        "frazione_chimerica": {
            classe: {
                "letture": _frazione(c["letture_chimeriche"], c["letture"]),
                "varianti": _frazione(c["varianti_chimeriche"], c["varianti"]),
            }
            for classe, c in sorted(per_classe.items())
        },
        "classi_controllate": [c.value for c in CLASSI_CONTROLLATE],
        "letture_controllate": letture,
        "letture_chimeriche_controllate": chimeriche,
        "frazione_controllata": _frazione(chimeriche, letture),
    }


def controlla_chimere(misure: dict[str, Any], qc: Qc) -> str | None:
    """I controlli E-S6-01 ed E-S6-02 sulla frazione controllata.

    Solleva E-S6-01 oltre ``qc.stop_frac_chimeric``; oltre
    ``qc.warn_frac_chimeric`` restituisce il dettaglio della degradazione da
    registrare, altrimenti ``None``.
    """
    frazione = misure["frazione_controllata"]
    if frazione is None:
        return None
    testo = (
        f"il {frazione:.1%} delle letture dei campioni biologici e dei controlli "
        f"positivi e' chimerico ({misure['letture_chimeriche_controllate']} su "
        f"{misure['letture_controllate']})"
    )
    if frazione > qc.stop_frac_chimeric:
        raise errore(
            "E-S6-01",
            f"{testo}, oltre il {qc.stop_frac_chimeric:.0%} di qc.stop_frac_chimeric",
            frazione=frazione,
        )
    if frazione > qc.warn_frac_chimeric:
        return f"{testo}, oltre il {qc.warn_frac_chimeric:.0%} di qc.warn_frac_chimeric"
    return None


class RimozioneChimere(PipelineStep):
    """S6: la tabella senza chimere."""

    passo: ClassVar[Passo] = Passo.S6
    versione: ClassVar[int] = 1
    script_r: ClassVar[str | None] = NOME_SCRIPT
    passi_tracciamento: ClassVar[tuple[str, ...]] = ("senza_chimere",)
    #: Il gruppo chimera e le due soglie sulla frazione chimerica.
    parametri: ClassVar[tuple[str, ...]] = (
        "chimera", "qc.warn_frac_chimeric", "qc.stop_frac_chimeric",
    )

    def calcola(self, contesto: StepContext) -> Produzione:
        """Esegue ``R/06_chimera.R`` sulla tabella di S5, scrive il riepilogo per classe
        e applica i controlli sulle chimere.
        """
        if contesto.inventario is None:
            raise RuntimeError("S6 richiede l'inventario prodotto da S0")
        config = contesto.config
        chimera = config.chimera
        letture = leggi_conteggi(contesto.albero.cartella(Fase.FILTERED) / NOME_FILTRATE)
        senza_letture = sorted(
            c.accession for c in contesto.inventario if letture.get(c.accession, 0) == 0
        )
        s5 = contesto.albero.manifesto_passo(Passo.S5.value, Fase.SEQTAB)
        if s5 is None or NOME_TABELLA not in s5.nomi:
            raise RuntimeError("S6 richiede la tabella delle sequenze di S5")
        voce_tabella = next(v for v in s5.artefatti if v["nome"] == NOME_TABELLA)

        esito = esegui_script(
            cartella_r() / NOME_SCRIPT,
            {
                "tabella": str(contesto.albero.cartella(Fase.SEQTAB) / NOME_TABELLA),
                "senza_letture": senza_letture,
                "gruppi": {c.accession: c.classe.value for c in contesto.inventario},
                "method": chimera.method,
                "min_fold_parent_over_abundance": chimera.min_fold_parent_over_abundance,
                "min_parent_abundance": chimera.min_parent_abundance,
                "min_sample_fraction": chimera.min_sample_fraction,
                "allow_one_off": chimera.allow_one_off,
                "processi": thread_effettivi(config),
            },
            contesto.albero,
            self.cartella,
            passo=self.passo,
            tempo_massimo_s=contesto.config.run.r_timeout_s,
            logger=contesto.logger,
        )
        cartella = contesto.albero.cartella(self.cartella)
        per_classe = leggi_per_gruppo(cartella / NOME_PER_GRUPPO)
        misure = misura_chimere(per_classe)
        # Il riepilogo si scrive prima dei controlli: un arresto lascia le
        # misure che l'hanno causato.
        riepilogo = {
            "tabella_prima_delle_chimere": {
                "fase": Passo.S5.value,
                "cartella": Fase.SEQTAB.value,
                "nome": NOME_TABELLA,
                "checksum": voce_tabella["checksum"],
            },
            "per_classe": per_classe,
            **misure,
        }
        artefatto = contesto.albero.scrivi_testo(
            self.cartella, NOME_RIEPILOGO,
            json.dumps(riepilogo, indent=2, ensure_ascii=False, sort_keys=True) + "\n",
        )
        avviso = controlla_chimere(misure, config.qc)
        if avviso is not None:
            contesto.degrada("E-S6-02", avviso, frazione=misure["frazione_controllata"])
        metriche = {
            "varianti_chimeriche": len(
                (cartella / "chimere.tsv").read_text(encoding="utf-8").splitlines()
            ) - 1,
            "frazione_controllata": misure["frazione_controllata"],
            "frazione_chimerica": misure["frazione_chimerica"],
        }
        return Produzione(esito.artefatti + (artefatto,), metriche)
