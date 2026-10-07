"""Fase S12: la decontaminazione dai controlli negativi.

Identifica le varianti che vengono dai reagenti e non dal campione con
``decontam::isContaminant`` per prevalenza: una variante piu' prevalente nei
controlli negativi che nei campioni e' un contaminante (``R/12_decontam.R``,
``R/lib/decontaminazione.R``). Lavora sull'oggetto integrato di S10 e scrive in
``11_controls/``, accanto a S11, con i propri file del ponte e il proprio
manifesto.

**Il metodo.** ``decontam.method`` e' ``prevalence``, l'unico praticabile: il
metodo per frequenza richiede una misura di concentrazione del DNA che i
metadati non hanno. ``decontam.threshold`` e' 0,5 e non 0,1: in bassa biomassa
la contaminazione da reagente puo' essere una quota rilevante del segnale, e la
soglia e' orientata alla sensibilita'.

**Chi entra nel confronto.** I controlli negativi da una parte, i campioni
biologici dall'altra, secondo la classe dei metadati. I controlli positivi ne
restano fuori: contengono un organismo aggiunto e, ai livelli diluiti,
contaminanti di reagente amplificati.

**Per piastra e aggregata.** Si calcolano entrambe: per piastra (la colonna
``decontam.batch_column``, nell'oggetto ``piastra``) con le probabilita'
combinate secondo ``decontam.batch_combine``, il ``batch.combine`` di decontam,
il cui predefinito e' ``minimum``; aggregata, tutti i biologici contro tutti i
negativi. Una piastra con meno di ``decontam.min_blanks`` negativi confronta i
propri biologici con i negativi di tutte le piastre. Quale delle due decide i
contaminanti lo dichiara ``decontam.mode``, non l'esito: l'altra resta una
diagnostica, con il confronto numerico in ``decontam_riepilogo.json``. Con
pochi negativi per piastra la modalita' per piastra e' permissiva per
costruzione: il minimo di molte probabilita', ciascuna stimata su pochi
negativi. Se la modalita' dichiarata rimuove dai biologici una frazione
delle letture oltre ``qc.max_frac_contaminant``, la fase si ferma con
``E-S12-02``.

**Senza abbastanza negativi.** Con meno di ``decontam.min_blanks`` controlli
negativi con letture in tutto il dataset la prevalenza nei negativi non e'
stimabile: nessuna variante viene rimossa, e la fase lo dichiara con
``E-S12-03``. L'oggetto prosegue con i suoi contaminanti: e' un risultato non
decontaminato, e il manifesto lo dice.

**La rimozione.** I contaminanti della modalita' dichiarata escono dall'oggetto in
tutti i campioni; le altre varianti tengono il loro identificativo, senza
rinumerazione. L'oggetto ripulito e' ``ps_decontaminato.rds``.

**L'ordine.** La decontaminazione precede il filtro di prevalenza (S13): i
contaminanti da reagente sono prevalenti per costruzione, e filtrare prima li
renderebbe indistinguibili dal segnale. Il grafo lo impone (``E-S13-01``).
"""

from __future__ import annotations

import json
from typing import Any, ClassVar, Final

from amplicon16s.errors.exceptions import errore
from amplicon16s.io_layer.artifacts import Fase
from amplicon16s.rbridge.runner import cartella_r, esegui_script
from amplicon16s.runner.graph import Passo
from amplicon16s.steps.base import PipelineStep, Produzione, StepContext
from amplicon16s.steps.s10_phyloseq import NOME_OGGETTO

__all__ = ["Decontaminazione", "NOME_OGGETTO_DECONTAMINATO"]

NOME_SCRIPT: Final = "12_decontam.R"
#: L'oggetto senza contaminanti che le fasi successive leggono.
NOME_OGGETTO_DECONTAMINATO: Final = "ps_decontaminato.rds"
NOME_RIEPILOGO: Final = "decontam_riepilogo.json"


class Decontaminazione(PipelineStep):
    """S12: i contaminanti da reagente, identificati e rimossi."""

    passo: ClassVar[Passo] = Passo.S12
    #: 3: con meno di decontam.min_blanks negativi con letture non toglie
    #: nulla e lo dichiara (E-S12-03).
    versione: ClassVar[int] = 3
    script_r: ClassVar[str | None] = NOME_SCRIPT
    passi_tracciamento: ClassVar[tuple[str, ...]] = ("decontaminate",)
    #: Il gruppo decontam (metodo, modalita', soglia, negativi minimi, colonna e
    #: combinazione delle piastre) e la frazione massima rimossa.
    parametri: ClassVar[tuple[str, ...]] = ("decontam", "qc.max_frac_contaminant")

    def calcola(self, contesto: StepContext) -> Produzione:
        """Esegue ``R/12_decontam.R`` sull'oggetto di S10 e si ferma con
        ``E-S12-02`` se la modalita' dichiarata rimuove troppo.
        """
        config = contesto.config
        decontam = config.decontam
        if decontam.method != "prevalence":
            raise RuntimeError(
                f"decontam.method {decontam.method!r} non e' realizzato: S12 decontamina "
                "per prevalenza, e i metadati non hanno la concentrazione del DNA che "
                "il metodo per frequenza richiede"
            )
        albero = contesto.albero
        esito = esegui_script(
            cartella_r() / NOME_SCRIPT,
            {
                "oggetto": str(albero.cartella(Fase.PHYLOSEQ) / NOME_OGGETTO),
                "soglia": decontam.threshold,
                "combinazione": decontam.batch_combine,
                "min_negativi": decontam.min_blanks,
                "modalita": decontam.mode,
                "max_frazione": config.qc.max_frac_contaminant,
            },
            albero,
            self.cartella,
            passo=self.passo,
            logger=contesto.logger,
        )
        riepilogo = json.loads(
            (albero.cartella(self.cartella) / NOME_RIEPILOGO).read_text(encoding="utf-8")
        )
        # Senza abbastanza negativi non c'e' stata decontaminazione: l'oggetto
        # passa a S13 con i suoi contaminanti, e va detto.
        if riepilogo["confronto"]["negativi"] < decontam.min_blanks:
            contesto.degrada("E-S12-03", riepilogo["esito"],
                             negativi=riepilogo["confronto"]["negativi"])
        if not riepilogo["entro_max_frazione"]:
            raise errore("E-S12-02", riepilogo["esito"], modalita=riepilogo["modalita"])
        metriche: dict[str, Any] = {
            "modalita_dichiarata": riepilogo["modalita_dichiarata"],
            "contaminanti_rimossi": riepilogo["contaminanti_rimossi"],
            "letture_rimosse": riepilogo["letture_rimosse"],
            "modalita": {
                nome: None if m is None else {
                    "contaminanti": m["contaminanti"],
                    "letture_rimosse_biologico": m["letture_rimosse"]["biologico"],
                }
                for nome, m in riepilogo["modalita"].items()
            },
        }
        return Produzione(esito.artefatti, metriche)
