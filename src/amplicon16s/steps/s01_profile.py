"""Fase S1: profilo di qualità e di lunghezza delle letture.

È la prima fase di calcolo: legge ogni file di letture **per intero**, a
differenza dei gate di S0 che si fermano alle prime ``qc.head_reads``, e
scrive in ``02_qc_profiles/`` per ciascun campione la distribuzione delle
lunghezze, il profilo di qualità per posizione e il conteggio delle letture
(``R/01_profile.R``, attraverso il ponte verso R).

**Chiude il limite noto di G09.** G09 stima la lunghezza minima dalle prime
letture di ogni file: una stima per eccesso, con cui il gate può passare
quando non dovrebbe. S1 conosce il minimo vero e ricontrolla la stessa
condizione: se ``filter.truncLen`` lo supera, le letture più corte verrebbero
scartate, e la fase si ferma con ``E-S1-02``. Il codice è di S1 e non è
``E-S0-09``: la condizione si verifica qui, e il log deve dire dove.

**Registra lo scarto in difetto, E-S1-01.** Se ``filter.truncLen`` è sotto il
minimo osservato di più di ``filter.truncLen_shortfall_warn``, ogni lettura
cede basi che si potrebbero conservare: è una degradazione, che non ferma e
finisce nel manifesto di S1. La registra S1 e non S0 perché la condizione si
valuta sul minimo vero misurato qui, non sulla stima dalle prime letture.

Sul dataset di riferimento il profilo di qualità è piatto (la qualità
mediana non scende sotto 25 in nessuna posizione), e il vincolo sul
troncamento non è la qualità ma l'eterogeneità delle lunghezze.
"""

from __future__ import annotations

import csv
import json
from collections import defaultdict
from typing import Any, ClassVar, Final

from amplicon16s.errors.exceptions import errore
from amplicon16s.rbridge.runner import cartella_r, esegui_script
from amplicon16s.runner.graph import Passo
from amplicon16s.steps.base import PipelineStep, Produzione, StepContext

__all__ = ["ProfiloLetture", "NOME_LUNGHEZZE", "NOME_QUALITA", "NOME_RIEPILOGO"]

NOME_SCRIPT: Final = "01_profile.R"
NOME_LUNGHEZZE: Final = "lunghezze.tsv"
NOME_QUALITA: Final = "qualita.tsv"
NOME_RIEPILOGO: Final = "riepilogo.json"


class ProfiloLetture(PipelineStep):
    """S1: profili di qualità e lunghezza, sul contenuto intero dei file."""

    passo: ClassVar[Passo] = Passo.S1
    versione: ClassVar[int] = 1
    script_r: ClassVar[str | None] = NOME_SCRIPT
    passi_tracciamento: ClassVar[tuple[str, ...]] = ("grezze",)
    #: I profili dipendono solo dalle letture, cioe' da S0; il troncamento e la
    #: sua tolleranza servono ai controlli E-S1-01 ed E-S1-02. Poiche' S1 non
    #: dichiara filter.maxEE, filter.truncQ ne' il gruppo err, la modifica di
    #: questi parametri a valle lascia valido il manifesto manifest_S1.json.
    parametri: ClassVar[tuple[str, ...]] = ("filter.truncLen", "filter.truncLen_shortfall_warn")

    def calcola(self, contesto: StepContext) -> Produzione:
        """Esegue ``R/01_profile.R`` sui file dell'inventario e confronta
        ``filter.truncLen`` con la lunghezza minima vera.
        """
        if contesto.inventario is None:
            raise RuntimeError("S1 richiede l'inventario prodotto da S0")

        campioni = {
            c.accession: str(c.file) for c in contesto.inventario if c.file is not None
        }
        esito = esegui_script(
            cartella_r() / NOME_SCRIPT,
            {"campioni": campioni, "processi": contesto.config.run.threads},
            contesto.albero,
            self.cartella,
            passo=self.passo,
            logger=contesto.logger,
        )

        cartella = contesto.albero.cartella(self.cartella)
        riepilogo = json.loads((cartella / NOME_RIEPILOGO).read_text(encoding="utf-8"))
        minimo = int(riepilogo["lunghezza_minima"])
        metriche = {
            "campioni": riepilogo["campioni"],
            "letture": riepilogo["letture"],
            "lunghezza_minima": minimo,
            "lunghezza_massima": riepilogo["lunghezza_massima"],
            "moda": riepilogo["moda"],
            "qualita_mediana_minima": riepilogo["qualita_mediana_minima"],
        }

        filtro = contesto.config.filter
        troncamento = filtro.truncLen
        if troncamento > minimo:
            raise errore(
                "E-S1-02",
                self._dettaglio_corte(cartella, troncamento, minimo),
                truncLen=troncamento,
                lunghezza_minima=minimo,
            )

        scarto = minimo - troncamento
        if scarto > filtro.truncLen_shortfall_warn:
            contesto.degrada(
                "E-S1-01",
                f"la lettura piu' corta, su tutte le letture, e' di {minimo} bp e "
                f"filter.truncLen vale {troncamento}: uno scarto di {scarto} bp, oltre "
                f"i {filtro.truncLen_shortfall_warn} di filter.truncLen_shortfall_warn",
                lunghezza_minima=minimo,
                truncLen=troncamento,
                scarto=scarto,
            )

        return Produzione(esito.artefatti, metriche)

    @staticmethod
    def _dettaglio_corte(cartella: Any, troncamento: int, minimo: int) -> str:
        """Quante letture perderebbe ciascun campione, dai peggiori."""
        corte: dict[str, int] = defaultdict(int)
        totali: dict[str, int] = defaultdict(int)
        with open(cartella / NOME_LUNGHEZZE, encoding="utf-8", newline="") as file:
            for riga in csv.DictReader(file, delimiter="\t"):
                letture = int(float(riga["letture"]))
                totali[riga["campione"]] += letture
                if int(riga["lunghezza"]) < troncamento:
                    corte[riga["campione"]] += letture
        peggiori = sorted(corte, key=lambda c: corte[c] / totali[c], reverse=True)
        elenco = ", ".join(
            f"{c} ({corte[c]} su {totali[c]})" for c in peggiori[:5]
        )
        return (
            f"filter.truncLen vale {troncamento} ma la lettura piu' corta, su tutte "
            f"le letture, e' di {minimo} bp: {len(corte)} campioni hanno letture "
            f"piu' corte, che verrebbero scartate. I piu' colpiti: {elenco}"
        )
