"""Fase S1: profilo di qualità e di lunghezza delle letture.

È la prima fase di calcolo: legge ogni file di letture **per intero**, a
differenza dei gate di S0 che si fermano alle prime ``qc.head_reads``, e
scrive in ``02_qc_profiles/`` per ciascun campione la distribuzione delle
lunghezze, il profilo di qualità per posizione e il conteggio delle letture
(``R/01_profile.R``, attraverso il ponte verso R).

**Chiude il limite noto di G09.** G09 guarda le prime letture di ogni file,
e può passare quando non dovrebbe. S1 conosce tutte le lunghezze e ripete la
stessa verifica: le letture più corte di ``filter.truncLen`` verrebbero
scartate dal filtro, non accorciate. Una sola lettura corta non ferma nulla:
la fase misura **per classe** la frazione di letture più corte, e si ferma
con ``E-S1-02`` se quella dei campioni biologici, o quella dei controlli
positivi, supera ``qc.max_frac_short_reads``: ogni classe si giudica da sola,
perché in una frazione unica pesata sulle letture una classe poco numerosa
potrebbe perdere tutte le sue letture senza farsi notare. I controlli
negativi non contano: amplificano poco, e le loro letture corte non dicono
nulla del troncamento. Il codice è di S1 e non è ``E-S0-09``: la condizione si verifica
qui, e il log deve dire dove. Le frazioni per classe restano nelle metriche
del manifesto.

**Chiude il limite noto di G07, E-S1-04.** G07 riconosce un file con le due
letture di ogni coppia dalle intestazioni delle prime ``qc.head_reads``
letture: un file con tutte le prime letture seguite da tutte le seconde gli
sfugge se le ispezionate non arrivano al secondo blocco. S1 conta gli stessi
segni (``/1`` e ``/2``, ``1:N:`` e ``2:N:``, lo stesso nome ripetuto) su tutte
le letture, con la stessa funzione e la stessa tolleranza
(``io_layer/reads.py``), registra i conteggi per campione in ``coppie.tsv`` e
si ferma con ``E-S1-04`` se un file contiene le due letture di ogni coppia, o
le sole seconde: proseguire darebbe il doppio delle letture, forward e
inverse mescolate. Il conteggio precede i profili: è una lettura dei soli
nomi, e un dataset da respingere si respinge prima del calcolo.

**Segnala le qualità raggruppate, E-S1-03.** Con quattro valori di qualità
distinti o meno (``valori_qualita.tsv``) le letture vengono da un
sequenziatore che raggruppa le qualità: è una degradazione dichiarata, che
rimanda a ``err.error_function``. Con un solo valore non c'è una funzione da
scegliere (S3 non ha una curva da stimare, ``E-S3-04``), e il messaggio lo dice.

**Registra lo scarto in difetto, E-S1-01.** Se ``filter.truncLen`` è sotto il
minimo osservato di più di ``filter.truncLen_shortfall_warn``, ogni lettura
cede basi che si potrebbero conservare: è una degradazione, che non ferma e
finisce nel manifesto di S1. La registra S1 e non S0 perché la condizione si
valuta sul minimo vero misurato qui, non sulla stima dalle prime letture.
"""

from __future__ import annotations

import csv
import json
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from typing import Any, ClassVar, Final

from amplicon16s.config.schema import thread_effettivi
from amplicon16s.errors.exceptions import errore
from amplicon16s.io_layer.reads import StatisticheFile, conta_coppie
from amplicon16s.metadata.models import CLASSI_CONTROLLATE, ClasseCampione
from amplicon16s.rbridge.runner import cartella_r, esegui_script
from amplicon16s.runner.graph import Passo
from amplicon16s.steps.base import PipelineStep, Produzione, StepContext

__all__ = [
    "ProfiloLetture", "classi_oltre_soglia", "letture_corte", "NOME_COPPIE", "NOME_LUNGHEZZE",
    "NOME_QUALITA", "NOME_RIEPILOGO",
]

NOME_SCRIPT: Final = "01_profile.R"
NOME_LUNGHEZZE: Final = "lunghezze.tsv"
NOME_QUALITA: Final = "qualita.tsv"
NOME_RIEPILOGO: Final = "riepilogo.json"
NOME_VALORI_QUALITA: Final = "valori_qualita.tsv"
#: I segni di coppia contati su tutte le letture di ogni file.
NOME_COPPIE: Final = "coppie.tsv"


#: Con questo numero di valori di qualita' distinti, o meno, le qualita' sono
#: raggruppate (NovaSeq, NextSeq) e S1 lo segnala con E-S1-03.
QUALITA_RAGGRUPPATE: Final = 4


def letture_corte(
    lunghezze: dict[str, dict[int, int]],
    classi: dict[str, ClasseCampione],
    troncamento: int,
) -> dict[str, dict[str, Any]]:
    """Per classe, le letture e quante sono piu' corte del troncamento.

    ``lunghezze`` da' per campione le letture di ogni lunghezza. La voce
    ``controllate`` riunisce le classi da cui ci si attende il segnale
    (biologici e controlli positivi) ed e' un riepilogo: il giudizio e' su
    ciascuna delle due classi, separatamente (:func:`classi_oltre_soglia`).
    """
    per_classe: dict[str, dict[str, Any]] = {}
    for campione, distribuzione in lunghezze.items():
        classe = classi[campione]
        gruppi = [classe.value] + (["controllate"] if classe in CLASSI_CONTROLLATE else [])
        totale = sum(distribuzione.values())
        corte = sum(n for lunghezza, n in distribuzione.items() if lunghezza < troncamento)
        for gruppo in gruppi:
            voce = per_classe.setdefault(gruppo, {"letture": 0, "corte": 0})
            voce["letture"] += totale
            voce["corte"] += corte
    for voce in per_classe.values():
        voce["frazione"] = round(voce["corte"] / voce["letture"], 6) if voce["letture"] else 0.0
    return dict(sorted(per_classe.items()))


def classi_oltre_soglia(corte: dict[str, dict[str, Any]], massima: float) -> list[str]:
    """Le classi controllate la cui frazione di letture corte supera la soglia.

    Ogni classe da sola: tre controlli positivi con tutte le letture corte fra
    cento biologici senza letture corte sono una classe persa per intero, che
    una frazione unica (il 3%) non mostrerebbe.
    """
    return [
        classe.value for classe in CLASSI_CONTROLLATE
        if corte.get(classe.value, {"frazione": 0.0})["frazione"] > massima
    ]


class ProfiloLetture(PipelineStep):
    """S1: profili di qualità e lunghezza, sul contenuto intero dei file."""

    passo: ClassVar[Passo] = Passo.S1
    #: 2: il troncamento si giudica sulla frazione di letture piu' corte, per
    #: classe, e non piu' sulla sola lettura piu' corta; si contano i valori di
    #: qualita' distinti. 3: la frazione si giudica su ciascuna classe
    #: controllata, non sulle due riunite. 4: conta su tutte le letture i segni
    #: di coppia che G07 cerca nelle prime, li scrive in coppie.tsv e ferma un
    #: file con le due letture di ogni coppia (E-S1-04).
    versione: ClassVar[int] = 4
    script_r: ClassVar[str | None] = NOME_SCRIPT
    passi_tracciamento: ClassVar[tuple[str, ...]] = ("grezze",)
    #: I profili dipendono solo dalle letture, cioe' da S0; il troncamento, la
    #: sua tolleranza e la frazione ammessa di letture piu' corte servono ai
    #: controlli E-S1-01 ed E-S1-02. Poiche' S1 non dichiara filter.maxEE,
    #: filter.truncQ ne' il gruppo err, la modifica di questi parametri a valle
    #: lascia valido il manifesto manifest_S1.json.
    parametri: ClassVar[tuple[str, ...]] = (
        "filter.truncLen", "filter.truncLen_shortfall_warn", "qc.max_frac_short_reads",
    )

    def calcola(self, contesto: StepContext) -> Produzione:
        """Esegue ``R/01_profile.R`` sui file dell'inventario e confronta
        ``filter.truncLen`` con le lunghezze vere.
        """
        if contesto.inventario is None:
            raise RuntimeError("S1 richiede l'inventario prodotto da S0")

        campioni = {
            c.accession: str(c.file) for c in contesto.inventario if c.file is not None
        }
        processi = thread_effettivi(contesto.config)
        # Un processo per file: la decompressione e la lettura dei nomi non
        # condividono nulla, e i nomi di un file vivono solo nel suo processo.
        with ProcessPoolExecutor(max_workers=processi) as gruppo:
            coppie = dict(zip(campioni, gruppo.map(conta_coppie, campioni.values())))
        tabella = contesto.albero.scrivi_testo(self.cartella, NOME_COPPIE, self._tsv_coppie(coppie))
        respinti = {c: s.coppie_nello_stesso_file for c, s in coppie.items()
                    if s.coppie_nello_stesso_file}
        if respinti:
            elenco = "; ".join(f"{c} ({coppie[c].nome}): {motivo}"
                               for c, motivo in list(respinti.items())[:5])
            raise errore(
                "E-S1-04",
                f"{len(respinti)} file su {len(coppie)}, su tutte le loro letture, non "
                f"contengono le sole prime letture: {elenco}"
                + ("" if len(respinti) <= 5 else f"; e altri {len(respinti) - 5}")
                + f". I conteggi di ogni file sono in {NOME_COPPIE}",
                campioni=sorted(respinti),
            )
        esito = esegui_script(
            cartella_r() / NOME_SCRIPT,
            {"campioni": campioni, "processi": processi},
            contesto.albero,
            self.cartella,
            passo=self.passo,
            logger=contesto.logger,
        )

        cartella = contesto.albero.cartella(self.cartella)
        riepilogo = json.loads((cartella / NOME_RIEPILOGO).read_text(encoding="utf-8"))
        minimo = int(riepilogo["lunghezza_minima"])
        filtro = contesto.config.filter
        troncamento = filtro.truncLen

        # Le letture piu' corte del troncamento, che il filtro scartera': una
        # sola non e' un motivo per fermarsi, una quota rilevante dei campioni
        # da cui ci si attende il segnale si'.
        lunghezze = self._lunghezze(cartella)
        classi = {c.accession: c.classe for c in contesto.inventario}
        corte = letture_corte(lunghezze, classi, troncamento)
        massima = contesto.config.qc.max_frac_short_reads
        oltre = classi_oltre_soglia(corte, massima)

        with open(cartella / NOME_VALORI_QUALITA, encoding="utf-8", newline="") as file:
            valori = [int(r["qualita"]) for r in csv.DictReader(file, delimiter="\t")]
        metriche = {
            "campioni": riepilogo["campioni"],
            "letture": riepilogo["letture"],
            "lunghezza_minima": minimo,
            "lunghezza_massima": riepilogo["lunghezza_massima"],
            "moda": riepilogo["moda"],
            "qualita_mediana_minima": riepilogo["qualita_mediana_minima"],
            "letture_piu_corte_del_troncamento": corte,
            "valori_di_qualita_distinti": len(valori),
        }

        if oltre:
            raise errore(
                "E-S1-02",
                self._dettaglio_corte(lunghezze, classi, corte, oltre, troncamento, minimo, massima),
                truncLen=troncamento,
                lunghezza_minima=minimo,
                frazione={classe: corte[classe]["frazione"] for classe in oltre},
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
        if len(valori) <= QUALITA_RAGGRUPPATE:
            # Con un solo valore non c'e' una curva da adattare, con nessuna
            # funzione: suggerire loess_monotono manderebbe a un secondo arresto.
            seguito = (
                "con un solo valore di qualita' il modello d'errore non e' stimabile "
                "con nessuna funzione di err.error_function, e S3 si fermera' "
                "(E-S3-04): servono le letture con le qualita' originali"
                if len(valori) <= 1 else
                "con err.error_function loess_monotono la stima del modello di errore "
                "e' vincolata a non crescere con la qualita'"
            )
            contesto.degrada(
                "E-S1-03",
                (f"le letture hanno un solo valore di qualita' ({valori[0]}): "
                 if len(valori) == 1 else
                 f"le letture hanno {len(valori)} valori di qualita' distinti "
                 f"({', '.join(str(v) for v in sorted(valori))}): qualita' raggruppate. ")
                + seguito,
                valori=sorted(valori),
            )

        return Produzione((tabella, *esito.artefatti), metriche)

    @staticmethod
    def _tsv_coppie(coppie: dict[str, StatisticheFile]) -> str:
        """``coppie.tsv``: per campione, le letture e i segni di coppia contati."""
        righe = ["campione\tletture\tprime_di_coppia\tseconde_di_coppia\tnomi_ripetuti\t"
                 "nomi_oltre_due\tcoppie_nello_stesso_file"]
        righe += [
            f"{campione}\t{s.letture_esaminate}\t{s.prime_di_coppia}\t{s.seconde_di_coppia}\t"
            f"{s.nomi_ripetuti}\t{s.nomi_oltre_due}\t"
            f"{'si' if s.coppie_nello_stesso_file else 'no'}"
            for campione, s in coppie.items()
        ]
        return "\n".join(righe) + "\n"

    @staticmethod
    def _lunghezze(cartella: Any) -> dict[str, dict[int, int]]:
        """Per campione, le letture di ogni lunghezza (da ``lunghezze.tsv``)."""
        lunghezze: dict[str, dict[int, int]] = defaultdict(dict)
        with open(cartella / NOME_LUNGHEZZE, encoding="utf-8", newline="") as file:
            for riga in csv.DictReader(file, delimiter="\t"):
                lunghezze[riga["campione"]][int(riga["lunghezza"])] = int(float(riga["letture"]))
        return dict(lunghezze)

    @staticmethod
    def _dettaglio_corte(
        lunghezze: dict[str, dict[int, int]],
        classi: dict[str, ClasseCampione],
        corte: dict[str, dict[str, Any]],
        oltre: list[str],
        troncamento: int,
        minimo: int,
        massima: float,
    ) -> str:
        """Le classi oltre la soglia, la frazione di ogni classe e i campioni
        piu' colpiti fra quelli delle classi oltre la soglia.
        """
        per_campione = {
            c: (sum(n for l, n in d.items() if l < troncamento), sum(d.values()))
            for c, d in lunghezze.items() if classi[c].value in oltre
        }
        colpiti = {c: v for c, v in per_campione.items() if v[0]}
        peggiori = sorted(colpiti, key=lambda c: colpiti[c][0] / colpiti[c][1], reverse=True)
        elenco = ", ".join(f"{c} ({colpiti[c][0]} su {colpiti[c][1]})" for c in peggiori[:5])
        per_classe = "; ".join(
            f"{classe} {voce['corte']} su {voce['letture']} ({voce['frazione']:.1%})"
            for classe, voce in corte.items() if classe != "controllate"
        )
        superano = "; ".join(
            f"il {corte[classe]['frazione']:.1%} delle letture della classe {classe} "
            f"({corte[classe]['corte']} su {corte[classe]['letture']})"
            for classe in oltre
        )
        return (
            f"filter.truncLen vale {troncamento} e {superano} e' piu' corto, oltre il "
            f"{massima:.1%} di qc.max_frac_short_reads; la lettura piu' corta e' di {minimo} "
            f"bp. Per classe: {per_classe}. {len(colpiti)} campioni di "
            f"{'questa classe' if len(oltre) == 1 else 'queste classi'} hanno letture piu' "
            f"corte, che verrebbero scartate; i piu' colpiti: {elenco}"
        )
