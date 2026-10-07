"""Fase S11: la validazione della corsa dai controlli positivi (KatharoSeq).

I controlli positivi contengono un solo organismo, ``katharoseq.target_taxon``,
in diluizione seriale. S11 misura per ciascuno la **fedelta'**, la frazione
delle sue letture assegnate al taxon atteso al rango
``katharoseq.collapse_rank`` (sommando tutte le varianti di quel taxon), e
adatta una curva che la lega alla profondita' di lettura: la sigmoide
allosterica ``f = x^h / (k' + x^h)``, con ``x = log10(profondita')``. Dalla
curva ricava la profondita' minima a cui la fedelta' raggiunge
``katharoseq.target_sensitivity``. Il calcolo e' in ``R/11_controls.R`` e
``R/lib/katharoseq.R``, che dichiarano l'equazione, la bonta' (R^2 sulla scala
della fedelta', ``1 - SS_res / SS_tot``) e la regola di scelta fra curva
aggregata e curve per piastra. Scrive in ``11_controls/``, accanto alla futura
S12.

**Profondita' e fedelta' sono misurate allo stadio** ``katharoseq.read_stage``
(``nonchimeric``): le letture dell'oggetto integrato di S10, senza chimere e
dopo il filtro di lunghezza.

**Il livello di diluizione viene dai dati.** E' il numero di cellule seminate,
nella colonna ``katharoseq.cell_count_column`` del file di arricchimento o
della tabella di studio, che S10 porta nell'oggetto; il nome che ha
nell'oggetto si legge da ``colonne_metadati.tsv`` di S10.

**Senza controlli positivi**, o senza la colonna dei livelli, non c'e' nulla
da stimare: nessuna curva, la soglia e' il ripiego per tutti i campioni, e la
fase lo dichiara con ``E-S11-05``, distinto da ``E-S11-02`` (controlli
presenti, curva non attendibile). Un controllo senza piastra entra nella sola
curva aggregata; uno con una sola lettura in nessuna, perche' la curva e'
definita sul logaritmo della profondita'.

**La soglia e' un risultato, non un parametro.** Sta in ``soglia.json``, per
piastra, e ogni valore porta lo stadio a cui si applica: la soglia derivata
vale sulle letture ``nonchimeric``; il ripiego ``qc.min_reads_raw`` sulle
letture grezze (``raw``, ``letture_prefiltro.tsv`` di S2). Con
``qc.min_reads_mode`` = ``katharoseq_if_available`` una piastra senza curva
valida ripiega, e la fase lo registra con ``E-S11-02`` e il motivo (colonna
dei livelli assente, controlli meno di ``ctrl.min_positives``, bonta' sotto
``katharoseq.min_r2``, soglia fuori dalle profondita' osservate); con
``fixed`` si usa sempre il ripiego, senza degradazione.

**La conformita' dei controlli** tiene conto della concentrazione: un
controllo e' non conforme se la sua fedelta' e' anomalamente bassa, o la sua
profondita' anomala, rispetto ai controlli dello stesso livello nelle altre
piastre (z modificato oltre 3,5). A bassa concentrazione la fedelta' attesa e'
bassa: un controllo dominato dai contaminanti come gli altri del suo livello
e' conforme. Se i conformi sono meno di ``ctrl.min_positive_pass_frac`` dei
controlli valutabili, con ``ctrl.positive_gate`` vero la fase si ferma con
``E-S11-03``, a revisione umana; con ``ctrl.positive_gate`` falso prosegue e
registra l'avviso ``E-S11-04``: un avviso e un arresto non condividono il
codice. I controlli non conformi non entrano nella curva.

``profondita_campioni.tsv`` riporta per ogni campione la soglia che gli si
applica e se vi cade sotto: e' una misura, il filtro e' di S13.
"""

from __future__ import annotations

import csv
import json
from typing import Any, ClassVar, Final

from amplicon16s.errors.exceptions import errore
from amplicon16s.io_layer.artifacts import Fase
from amplicon16s.rbridge.runner import cartella_r, esegui_script
from amplicon16s.runner.graph import Passo
from amplicon16s.steps.base import PipelineStep, Produzione, StepContext
from amplicon16s.metadata.tabelle import COLONNE_INVENTARIO
from amplicon16s.steps.s10_phyloseq import NOME_COLONNE, NOME_OGGETTO

__all__ = ["NOME_SOGLIA", "ValidazioneControlli", "colonna_dei_livelli"]

NOME_SCRIPT: Final = "11_controls.R"
#: La soglia di profondita' che le fasi successive leggono.
NOME_SOGLIA: Final = "soglia.json"
NOME_RIEPILOGO: Final = "riepilogo.json"
#: Le letture grezze per campione, scritte da S2.
NOME_GREZZE: Final = "letture_prefiltro.tsv"


def colonna_dei_livelli(
    corrispondenza: list[dict[str, str]], originale: str | None
) -> str | None:
    """Il nome nell'oggetto integrato della colonna dei livelli, se S10 l'ha portata.

    La colonna si cerca per nome originale fra quelle aggiunte all'oggetto dal
    file di arricchimento o dalla tabella di studio (``out.batch_columns``,
    ``out.study_columns``), non fra quelle dell'inventario.
    """
    if originale is None:
        return None
    for voce in corrispondenza:
        if voce["colonna_originale"] == originale and voce["colonna"] not in COLONNE_INVENTARIO:
            return voce["colonna"]
    return None


class ValidazioneControlli(PipelineStep):
    """S11: la curva KatharoSeq, la soglia di profondita' e la conformita' dei positivi."""

    passo: ClassVar[Passo] = Passo.S11
    #: 2: la colonna dei livelli si cerca fra tutte le colonne portate
    #: nell'oggetto, dal file di arricchimento o dalla tabella di studio.
    #: 3: senza controlli positivi o senza la colonna dei livelli dichiara
    #: E-S11-05 e non adatta curve; un controllo senza piastra o con una sola
    #: lettura non ferma la fase.
    versione: ClassVar[int] = 3
    script_r: ClassVar[str | None] = NOME_SCRIPT
    #: La curva (katharoseq), i controlli minimi, la frazione di conformi e il
    #: comportamento sotto di essa (ctrl), la regola e il valore del ripiego (qc).
    parametri: ClassVar[tuple[str, ...]] = (
        "katharoseq",
        "ctrl.min_positives", "ctrl.min_positive_pass_frac", "ctrl.positive_gate",
        "qc.min_reads_mode", "qc.min_reads_raw",
    )

    def calcola(self, contesto: StepContext) -> Produzione:
        """Esegue ``R/11_controls.R`` sull'oggetto di S10, poi registra il ripiego
        (``E-S11-02``) e applica il controllo sulla conformita' dei positivi.
        """
        config = contesto.config
        katharoseq = config.katharoseq
        albero = contesto.albero
        with open(albero.cartella(Fase.PHYLOSEQ) / NOME_COLONNE, encoding="utf-8",
                  newline="") as file:
            corrispondenza = list(csv.DictReader(file, delimiter="\t"))
        colonna = colonna_dei_livelli(corrispondenza, katharoseq.cell_count_column)

        esito = esegui_script(
            cartella_r() / NOME_SCRIPT,
            {
                "oggetto": str(albero.cartella(Fase.PHYLOSEQ) / NOME_OGGETTO),
                "colonna_cellule": colonna,
                "motivo_colonna": (
                    "katharoseq.cell_count_column non dichiarata"
                    if katharoseq.cell_count_column is None else
                    f"la colonna {katharoseq.cell_count_column!r} "
                    "(katharoseq.cell_count_column) non e' fra le colonne dell'oggetto "
                    "integrato"
                ),
                "rango": katharoseq.collapse_rank,
                "taxon": katharoseq.target_taxon,
                "sensibilita": katharoseq.target_sensitivity,
                "min_r2": katharoseq.min_r2,
                "stadio": katharoseq.read_stage,
                "min_positivi": config.ctrl.min_positives,
                "modo": config.qc.min_reads_mode,
                "ripiego": config.qc.min_reads_raw,
                "letture_grezze": str(albero.cartella(Fase.FILTERED) / NOME_GREZZE),
            },
            albero,
            self.cartella,
            passo=self.passo,
            logger=contesto.logger,
        )
        cartella = albero.cartella(self.cartella)
        riepilogo = json.loads((cartella / NOME_RIEPILOGO).read_text(encoding="utf-8"))
        soglia = json.loads((cartella / NOME_SOGLIA).read_text(encoding="utf-8"))

        # Senza controlli positivi, o senza la colonna dei livelli, la curva
        # manca per costruzione e non per un difetto dei controlli: e' un'altra
        # dichiarazione, perche' chi legge deve sapere che non c'era nulla da
        # stimare. Il valore applicato e' lo stesso ripiego.
        if riepilogo["positivi"] == 0 or colonna is None:
            if config.qc.min_reads_mode != "fixed":
                contesto.degrada(
                    "E-S11-05",
                    f"{riepilogo['motivo_scelta']}: qc.min_reads_raw "
                    f"({config.qc.min_reads_raw}) sulle letture grezze per tutti i campioni",
                    scelta=soglia["scelta"],
                )
        elif soglia["ripiego"]["degradazione"]:
            motivi = "; ".join(
                f"{m['piastra']}: {m['motivo']}" for m in soglia["ripiego"]["motivi"]
            )
            contesto.degrada(
                "E-S11-02",
                f"qc.min_reads_raw ({config.qc.min_reads_raw}) sulle letture grezze dove "
                f"la curva manca, per piastra: {motivi}",
                scelta=soglia["scelta"],
            )

        frazione = riepilogo["frazione_conformi"]
        minima = config.ctrl.min_positive_pass_frac
        if frazione is not None and frazione < minima:
            dettaglio = (
                f"{riepilogo['conformi']} controlli positivi conformi su "
                f"{riepilogo['conformi'] + riepilogo['non_conformi']} valutabili "
                f"({frazione:.4f}, minimo {minima}); non conformi: "
                f"{', '.join(riepilogo['non_conformi_elenco'])}"
            )
            if config.ctrl.positive_gate:
                raise errore("E-S11-03", dettaglio, frazione=frazione)
            contesto.degrada("E-S11-04", dettaglio, frazione=frazione)

        metriche: dict[str, Any] = {
            k: riepilogo[k] for k in (
                "positivi", "conformi", "non_conformi", "non_valutabili",
                "frazione_conformi", "scelta", "biologici_sotto_soglia",
            )
        }
        # Senza piastre R scrive un elenco vuoto invece di un oggetto.
        per_piastra = soglia["per_piastra"] or {}
        metriche["soglie"] = {
            p: {"valore": s["valore"], "stadio": s["stadio"]} for p, s in per_piastra.items()
        }
        return Produzione(esito.artefatti, metriche)
