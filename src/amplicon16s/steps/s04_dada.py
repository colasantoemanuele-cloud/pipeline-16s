"""Fase S4: l'inferenza delle varianti di sequenza.

Separa, campione per campione, le sequenze vere dagli errori di
sequenziamento (``R/04_dada.R``, ``dada2::dada``, attraverso il ponte). Ogni
campione usa il modello d'errore della propria corsa, letto da
``corrispondenza.tsv`` di S3 e non ricalcolato qui.

**Il pseudo-pooling e' realizzato a lotti, in due passate esplicite.**
``dada(pool = "pseudo")`` in una chiamata sola tiene in memoria i risultati di
tutti i campioni, e su centinaia di campioni non e' praticabile. Lo script fa le due
passate di dada, entrambe sui campioni a lotti di ``run.batch_size``: la prima
elabora ogni campione da solo e ne conserva soltanto, per ogni sequenza, in
quanti campioni compare e quante letture ha, e per ciascun modello la somma
delle transizioni; la seconda rielabora ogni campione con le informazioni a
priori. La regola che le sceglie e' quella del sorgente di dada2 1.36.0: le
varianti della prima passata presenti in almeno ``PSEUDO_PREVALENCE`` campioni
(2) oppure con almeno ``PSEUDO_ABUNDANCE`` letture in totale (infinito, quindi
mai). Le soglie si leggono dalla libreria.

Il sorgente mostra anche un fatto che la documentazione non dice: la seconda
passata non usa il modello d'errore fornito, ma quello che ``dada`` ricalcola
alla fine della prima (``loessErrfun`` sulle transizioni di tutti i campioni
della chiamata), anche con ``selfConsist = FALSE``. Misurato sulla versione
ridotta: la seconda passata con il modello di S3 da' varianti diverse da
``dada(pool = "pseudo")`` in 6 campioni su 11 e in 17 su 17 delle due corse;
con il modello ricalcolato coincide in tutti. Lo script lo riproduce, per
ciascun modello di S3 dalle transizioni dei soli suoi campioni, e lo registra
come artefatto. La corsa di un campione continua quindi a decidere il suo
modello in entrambe le passate: nella prima e' quello di S3, nella seconda la
sua stima aggiornata sulla stessa corsa.

**Le informazioni a priori si raccolgono su tutti i campioni, di tutte le
corse.** Il modello d'errore descrive la macchina, e va separato per corsa;
una variante presente in due campioni descrive il campione biologico, e che
i due campioni siano stati sequenziati nella stessa corsa non la rende piu' o
meno vera. L'informazione a priori abbassa soltanto la soglia con cui una
sequenza gia' presente in un campione viene riconosciuta come variante: non
aggiunge sequenze che il campione non ha.

**Il lotto non cambia il risultato.** Ogni campione e' elaborato da solo in
entrambe le passate, e cio' che attraversa i lotti sono somme di conteggi
interi, esatte in qualunque ordine. Il test lo verifica confrontando byte per
byte gli artefatti con due valori di ``run.batch_size``, e confronta il
risultato con ``dada(pool = "pseudo")`` in una sola chiamata. Per questo
``run.batch_size`` e' fra i parametri senza effetto, e il retry di E-S4-02
che lo dimezza e' legittimo: la memoria dipende dal lotto, il risultato no.
Con ``dada.pool`` vero dada riunisce tutti i campioni di un modello in
un'unica inferenza, il lotto non ha effetto, e la fase dichiara inutile il
retry.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import ClassVar, Final

from amplicon16s.config.schema import thread_effettivi
from amplicon16s.errors.exceptions import ErrorePipeline, errore
from amplicon16s.io_layer.artifacts import Fase
from amplicon16s.io_layer.conteggi import leggi_conteggi
from amplicon16s.rbridge.runner import cartella_r, esegui_script
from amplicon16s.runner.graph import Passo
from amplicon16s.runner.retry import RITENTARE_INUTILE, Aggiustamento, dimezza
from amplicon16s.steps.base import PipelineStep, Produzione, StepContext
from amplicon16s.steps.s02_filter import (
    NOME_FILTRATE,
    SUFFISSO_FILTRATI,
    SUFFISSO_FILTRATI_REVERSE,
)
from amplicon16s.steps.s03_learn_errors import NOME_CORRISPONDENZA

__all__ = [
    "InferenzaVarianti", "MODELLO_DEL_NOME_PER_LA_FUSIONE", "NOME_VARIANTI",
    "leggi_corrispondenza",
]

NOME_SCRIPT: Final = "04_dada.R"
#: Le varianti di ogni campione, ingresso di S5.
NOME_VARIANTI: Final = "varianti_per_campione.rds"
NOME_RIEPILOGO: Final = "inferenza.json"


#: Il nome del file in cui S4, in modalita' appaiata, conserva di un campione
#: cio' che la fusione legge per una lettura (``forward`` o ``reverse``): lo
#: stesso modello dello script R.
MODELLO_DEL_NOME_PER_LA_FUSIONE: Final = "fusione_{lettura}_{accession}.rds"


def leggi_corrispondenza(percorso: Path, colonna: str = "modello") -> dict[str, str]:
    """``corrispondenza.tsv`` di S3: campione -> file del modello d'errore. In
    modalita' appaiata la colonna ``modello_reverse`` da' quello della lettura
    reverse.
    """
    with open(percorso, encoding="utf-8", newline="") as file:
        return {r["campione"]: r[colonna] for r in csv.DictReader(file, delimiter="\t")}


class InferenzaVarianti(PipelineStep):
    """S4: le varianti di sequenza di ogni campione."""

    passo: ClassVar[Passo] = Passo.S4
    #: 2: la seconda passata usa la stessa loess_monotono di S3, cambiata
    #: (primo grado, span 2, pesi logaritmici). 3: se il modello della seconda
    #: passata non si stima, l'errore riporta la causa e le transizioni della
    #: prima; il risultato non cambia.
    versione: ClassVar[int] = 3
    script_r: ClassVar[str | None] = NOME_SCRIPT
    passi_tracciamento: ClassVar[tuple[str, ...]] = ("denoised",)
    #: Il gruppo dada: la modalita' di pooling e la soglia omega_a. Le letture
    #: filtrate e i modelli d'errore arrivano da S2 e S3, la cui impronta entra
    #: in quella di S4; run.batch_size non vi entra, perche' non cambia il
    #: risultato (vedi sopra).
    parametri: ClassVar[tuple[str, ...]] = ("dada", "err.error_function")
    #: In modalita' appaiata la fase non usa parametri del gruppo pair: inferisce
    #: anche le letture reverse, che S2 ha filtrato e S3 modellato, e dipende
    #: dalla sola modalita', che ogni fase dichiara.
    passi_tracciamento_appaiata: ClassVar[tuple[str, ...]] = ("denoised_reverse",)
    aggiustamenti: ClassVar[dict[str, Aggiustamento]] = {
        "E-S4-02": dimezza("run.batch_size"),
    }

    def calcola(self, contesto: StepContext) -> Produzione:
        """Esegue ``R/04_dada.R`` sui campioni con letture filtrate, ciascuno con il
        modello d'errore della sua corsa.
        """
        if contesto.inventario is None:
            raise RuntimeError("S4 richiede l'inventario prodotto da S0")
        config = contesto.config
        filtrati = contesto.albero.cartella(Fase.FILTERED)
        modelli_dir = contesto.albero.cartella(Fase.ERROR_MODELS)
        letture = leggi_conteggi(filtrati / NOME_FILTRATE)
        corrispondenza = leggi_corrispondenza(modelli_dir / NOME_CORRISPONDENZA)

        accession = sorted(c.accession for c in contesto.inventario)
        con_letture = [a for a in accession if letture.get(a, 0) > 0]
        senza_modello = [a for a in con_letture if a not in corrispondenza]
        if senza_modello:
            raise RuntimeError(
                f"campioni senza modello d'errore in {NOME_CORRISPONDENZA}: "
                + ", ".join(senza_modello)
            )

        try:
            esito = esegui_script(
                cartella_r() / NOME_SCRIPT,
                {
                    "campioni": {a: str(filtrati / f"{a}{SUFFISSO_FILTRATI}") for a in con_letture},
                    "modelli": {a: str(modelli_dir / corrispondenza[a]) for a in con_letture},
                    "senza_letture": [a for a in accession if a not in con_letture],
                    **self._della_reverse(config, filtrati, modelli_dir, con_letture),
                    "pool": config.dada.pool,
                    "omega_a": config.dada.omega_a,
                    "funzione_errore": config.err.error_function,
                    "lotto": config.run.batch_size,
                    "processi": thread_effettivi(config),
                },
                contesto.albero,
                self.cartella,
                passo=self.passo,
                codice_memoria="E-S4-02",
                tempo_massimo_s=contesto.config.run.r_timeout_s,
                logger=contesto.logger,
            )
        except ErrorePipeline as e:
            if e.codice == "E-S4-02" and config.dada.pool is True:
                inutile = (
                    "con dada.pool vero i campioni di ciascun modello si inferiscono "
                    "insieme in un'unica chiamata, e run.batch_size non cambia la "
                    "memoria richiesta"
                )
                raise errore(
                    "E-S4-02", e.dettaglio, **{**e.contesto, RITENTARE_INUTILE: inutile}
                ) from e
            raise

        cartella = contesto.albero.cartella(self.cartella)
        riepilogo = json.loads((cartella / NOME_RIEPILOGO).read_text(encoding="utf-8"))
        metriche = {
            "campioni": len(con_letture),
            "varianti_distinte": riepilogo["varianti_distinte"],
            "letture_denoised": sum(leggi_conteggi(cartella / "letture_denoised.tsv").values()),
        }
        if "priori" in riepilogo:
            metriche["priori"] = riepilogo["priori"]["priori"]
        if config.appaiata:
            reverse = json.loads(
                (cartella / "inferenza_reverse.json").read_text(encoding="utf-8"))
            metriche["reverse"] = {
                "varianti_distinte": reverse["varianti_distinte"],
                "letture_denoised": sum(
                    leggi_conteggi(cartella / "letture_denoised_reverse.tsv").values()),
            }
        return Produzione(esito.artefatti, metriche)

    @staticmethod
    def _della_reverse(
        config: object, filtrati: Path, modelli_dir: Path, con_letture: list[str]
    ) -> dict[str, dict[str, str]]:
        """In modalita' appaiata, le letture reverse filtrate di ogni campione e
        il modello d'errore della sua corsa per quella lettura; vuoto altrimenti.
        """
        if not config.appaiata:  # type: ignore[attr-defined]
            return {}
        del_reverse = leggi_corrispondenza(modelli_dir / NOME_CORRISPONDENZA, "modello_reverse")
        return {
            "reverse": {
                a: str(filtrati / f"{a}{SUFFISSO_FILTRATI_REVERSE}") for a in con_letture},
            "modelli_reverse": {a: str(modelli_dir / del_reverse[a]) for a in con_letture},
        }
