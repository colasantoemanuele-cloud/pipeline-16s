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

**Riconosce le letture inverse prese per forward, E-S1-07.** Nessun segno di
coppia distingue un file con le sole seconde letture, o con le due letture a
blocchi, quando le intestazioni numerano le letture come in un file
single-end: lo distingue il contenuto. Nella stessa passata S1 conta, in ogni
tratto di ogni file e nelle letture di posto dispari e pari, quelle con il
segnale di una lettura forward (il motivo conservato dopo il taglio, oppure il
primer nelle basi tagliate: ``io_layer/reads.py``), e si ferma con ``E-S1-07``
se un file di un campione biologico o di un controllo positivo quasi non ne
ha, o se in un file di qualunque classe il segnale sta in una parte del file e
manca in un'altra. La misura serve al solo
giudizio: non entra in ``coppie.tsv`` ne' nelle metriche. Senza motivo e senza
un primer nelle basi tagliate il controllo non si applica, e la fase lo
dichiara nel log.

**Un file che non si legge fino in fondo ferma la fase, E-S1-05.** La stessa
lettura completa trova un record troncato, una riga vuota in mezzo al file o
un archivio che non si decomprime fino in fondo, oltre le letture ispezionate
da S0; le sole righe vuote in fondo al file sono tollerate. Il profilo in R non lo direbbe: il lettore
a blocchi scarta in silenzio un record incompleto, e il conteggio dei segni di
coppia, fermo al punto del guasto, non varrebbe per il file intero.

**Memoria esaurita, E-S1-06.** Il profilo legge un file per processo
(``run.threads``), ognuno con un blocco di letture in memoria: se il sistema
ne uccide alcuni, il ponte lo riconosce e la fase si ferma con ``E-S1-06``, a
revisione umana, che suggerisce più memoria o meno processi.

**Segnala le qualità raggruppate, E-S1-03.** Con quattro valori di qualità
distinti o meno (``valori_qualita.tsv``) le letture vengono da un
sequenziatore che raggruppa le qualità: è una degradazione dichiarata, che
rimanda a ``err.error_function``. Con un solo valore non c'è una funzione da
scegliere (S3 non ha una curva da stimare, ``E-S3-04``), e il messaggio lo dice.

**In modalita' appaiata** la fase legge anche il file reverse di ogni campione
e ripete su di esso gli stessi controlli, con i parametri della lettura
reverse. In piu':

* confronta i nomi delle letture dei due file, una coppia per volta, e si
  ferma con ``E-S1-08`` alla prima coppia fuori posto o se i due file hanno un
  numero diverso di letture: la fusione di coppie non corrispondenti darebbe
  sequenze inesistenti, e con questa verifica il filtro non deve piu'
  confrontare gli identificativi;
* cerca in ogni file il segnale di entrambi i versi: un file forward senza il
  suo segnale il cui compagno lo porta, e viceversa, e' una corsa con i due
  file scambiati (``E-S1-09``, come il marcatore dell'altra lettura nelle
  intestazioni); un file reverse senza il suo segnale e' ``E-S1-09``; un file
  con una quota rilevante dei due versi e' una libreria mescolata
  (``E-S1-10``), che darebbe ogni variante due volte;
* scrive il profilo della lettura reverse accanto a quello della forward, con
  il suffisso ``_reverse``.

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
from concurrent.futures.process import BrokenProcessPool
from functools import partial
from pathlib import Path
from typing import Any, ClassVar, Final

from amplicon16s.config.schema import thread_effettivi
from amplicon16s.errors.exceptions import errore
from amplicon16s.io_layer.reads import (
    LETTURE_MINIME_PER_LA_QUOTA,
    LETTURE_MINIME_PER_UN_NOME_SOSPETTO,
    StatisticheFile,
    conta_coppie,
    corsa_invertita,
    nome_da_seconda_lettura,
    ricerca_del_primer,
    segnale_di_orientamento,
    senza_righe_vuote_finali,
    verifica_sincronia,
    versi_mescolati,
)
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
#: In modalita' appaiata i file della lettura reverse portano questo suffisso
#: prima dell'estensione: lunghezze_reverse.tsv, qualita_reverse.tsv e cosi' via.
SUFFISSO_REVERSE: Final = "_reverse"


def nome_reverse(nome: str) -> str:
    """Il nome dell'artefatto della lettura reverse che corrisponde a ``nome``."""
    base, punto, estensione = nome.rpartition(".")
    return f"{base}{SUFFISSO_REVERSE}{punto}{estensione}"


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


def _descrivi_il_segnale(config: Any) -> str:
    """Il segnale di una lettura forward con questa configurazione, a parole."""
    taglio = config.filter.trimLeft
    parti = []
    if config.taglio_cercato:
        parti.append(
            f"il primer {config.qc.primer_sequence} entro le prime "
            f"{config.filter.primer_max_offset + len(config.qc.primer_sequence)} basi")
        if config.qc.conserved_motif is not None:
            parti.append(f"il motivo {config.qc.conserved_motif!r} all'inizio della lettura")
        return " oppure ".join(parti)
    if config.qc.conserved_motif is not None:
        parti.append(
            f"il motivo {config.qc.conserved_motif!r} "
            + ("all'inizio della lettura" if taglio == 0 else f"alla posizione {taglio}"))
    if config.qc.primer_sequence is not None and taglio > 0:
        parti.append(f"il primer {config.qc.primer_sequence} in testa")
    return " oppure ".join(parti)


class ProfiloLetture(PipelineStep):
    """S1: profili di qualità e lunghezza, sul contenuto intero dei file."""

    passo: ClassVar[Passo] = Passo.S1
    #: 2: il troncamento si giudica sulla frazione di letture piu' corte, per
    #: classe, e non piu' sulla sola lettura piu' corta; si contano i valori di
    #: qualita' distinti. 3: la frazione si giudica su ciascuna classe
    #: controllata, non sulle due riunite. 4: conta su tutte le letture i segni
    #: di coppia che G07 cerca nelle prime, li scrive in coppie.tsv e ferma un
    #: file con le due letture di ogni coppia (E-S1-04). 5: ferma un file che
    #: non si legge fino in fondo (E-S1-05). 6: le righe vuote in fondo a un
    #: file sono tollerate; ogni guasto di decompressione e' E-S1-05; la
    #: memoria esaurita nel profilo e' dichiarata con E-S1-06. 7: fra i nomi
    #: ripetuti contano le letture consecutive che si distinguono per il solo
    #: marcatore di coppia in fondo all'identificativo. 8: il file senza
    #: letture di un controllo negativo non e' un guasto: il campione entra nel
    #: conteggio con zero letture e non ha un profilo. 9: ferma un file che non
    #: contiene le sole letture forward (E-S1-07) e una riga finale di soli
    #: spazi (E-S1-05). 10: un file con il nome di una seconda lettura e'
    #: giudicato sulla quota del segnale da un numero minore di letture; in
    #: modalita' appaiata legge i due file di ogni campione.
    versione: ClassVar[int] = 10
    script_r: ClassVar[str | None] = NOME_SCRIPT
    passi_tracciamento: ClassVar[tuple[str, ...]] = ("grezze",)
    #: I profili dipendono solo dalle letture, cioe' da S0; il troncamento, la
    #: sua tolleranza e la frazione ammessa di letture piu' corte servono ai
    #: controlli E-S1-01 ed E-S1-02. Poiche' S1 non dichiara filter.maxEE,
    #: filter.truncQ ne' il gruppo err, la modifica di questi parametri a valle
    #: lascia valido il manifesto manifest_S1.json. Primer, motivo e taglio
    #: definiscono il segnale di una lettura forward (E-S1-07).
    parametri: ClassVar[tuple[str, ...]] = (
        "filter.truncLen", "filter.truncLen_shortfall_warn", "qc.max_frac_short_reads",
        "qc.primer_sequence", "qc.conserved_motif", "filter.trimLeft",
    )
    #: Il troncamento della reverse e la sua tolleranza (E-S1-01, E-S1-02);
    #: primer, motivo e taglio della reverse, che ne definiscono il segnale
    #: (E-S1-09, E-S1-10).
    parametri_appaiata: ClassVar[tuple[str, ...]] = (
        "pair.truncLen_reverse", "pair.truncLen_shortfall_warn_reverse",
        "pair.primer_sequence_reverse", "pair.conserved_motif_reverse",
        "pair.trimLeft_reverse",
    )
    #: Con il primer cercato: la regola di ricerca, con cui si misurano le
    #: lunghezze dopo il taglio e si riconosce l'orientamento.
    parametri_taglio_cercato: ClassVar[tuple[str, ...]] = (
        "filter.primer_max_mismatch", "filter.primer_max_offset", "filter.primer_untrimmed",
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
        segnale = segnale_di_orientamento(
            contesto.config.qc.primer_sequence, contesto.config.qc.conserved_motif,
            contesto.config.filter.trimLeft, ricerca_del_primer(contesto.config))
        if segnale is None:
            contesto.logger.warning(
                "il controllo sulle letture inverse prese per forward non si applica: "
                "qc.conserved_motif e' nullo e non c'e' un primer da cercare nelle basi "
                "tagliate (filter.trimLeft 0, o qc.primer_sequence nullo). Un file con le "
                "seconde letture di ogni coppia non sarebbe riconosciuto dal contenuto",
                extra={"passo": "S1"},
            )
        appaiata = contesto.config.appaiata
        reverse: dict[str, str] = {}
        segnale_reverse = None
        if appaiata:
            reverse = {
                c.accession: str(c.file_reverse) for c in contesto.inventario
                if c.file is not None and c.file_reverse is not None
            }
            delle_coppie = contesto.config.pair
            segnale_reverse = segnale_di_orientamento(
                delle_coppie.primer_sequence_reverse, delle_coppie.conserved_motif_reverse,
                delle_coppie.trimLeft_reverse, ricerca_del_primer(contesto.config, reverse=True))
            if segnale_reverse is None:
                contesto.logger.warning(
                    "il controllo sull'orientamento delle letture reverse non si applica: "
                    "pair.conserved_motif_reverse e' nullo e non c'e' un primer reverse da "
                    "cercare nelle basi tagliate (pair.trimLeft_reverse 0, o "
                    "pair.primer_sequence_reverse nullo). Una corsa con i due file scambiati "
                    "si riconosce allora solo dal segnale delle letture forward",
                    extra={"passo": "S1"},
                )
        coppie_reverse: dict[str, StatisticheFile] = {}
        sincronie: dict[str, Any] = {}
        # Un processo per file: la decompressione e la lettura dei nomi non
        # condividono nulla, e i nomi di un file vivono solo nel suo processo.
        try:
            with ProcessPoolExecutor(max_workers=processi) as gruppo:
                coppie = dict(zip(campioni, gruppo.map(
                    partial(conta_coppie, orientamento=segnale,
                            altro_orientamento=segnale_reverse if appaiata else None),
                    campioni.values())))
                if appaiata:
                    coppie_reverse = dict(zip(reverse, gruppo.map(
                        partial(conta_coppie, orientamento=segnale_reverse,
                                altro_orientamento=segnale),
                        reverse.values())))
                    sincronie = dict(zip(reverse, gruppo.map(
                        verifica_sincronia, [campioni[c] for c in reverse], reverse.values())))
        except BrokenProcessPool as guasto:
            # Un processo del conteggio e' stato ucciso dal sistema: tiene in
            # memoria i nomi delle letture del suo file, e con molti processi
            # su file grandi la memoria puo' non bastare.
            raise errore(
                "E-S1-06",
                "un processo del conteggio dei segni di coppia e' stato interrotto dal "
                f"sistema operativo, di norma per memoria esaurita ({guasto})",
            ) from guasto
        # Un file che non si legge per intero non ha un conteggio: fermarsi
        # prima di scriverlo, perche' zero segni su una lettura interrotta non
        # dicono che il file e' single-end. Il file senza letture di un
        # controllo negativo non e' un guasto: S0 lo ha ammesso (E-S0-20).
        classi = {c.accession: c.classe for c in contesto.inventario}
        vuoti = [c for c, s in coppie.items()
                 if s.vuoto and classi[c] is ClasseCampione.CONTROLLO_NEGATIVO]
        guasti = {c: s for c, s in coppie.items() if s.errore is not None and c not in vuoti}
        # In modalita' appaiata anche i file reverse: quello senza letture di un
        # controllo negativo segue il suo forward.
        guasti.update({
            f"{c} (reverse)": s for c, s in coppie_reverse.items()
            if s.errore is not None and not (
                s.vuoto and classi[c] is ClasseCampione.CONTROLLO_NEGATIVO)
        })
        if guasti:
            elenco = "; ".join(f"{c} ({s.nome}): {s.errore}" for c, s in list(guasti.items())[:5])
            raise errore(
                "E-S1-05",
                f"{len(guasti)} file su {len(coppie) + len(coppie_reverse)} non si leggono "
                f"fino in fondo: {elenco}"
                + ("" if len(guasti) <= 5 else f"; e altri {len(guasti) - 5}")
                + ". Riscarica i file dalla sorgente e verificane il checksum",
                campioni=sorted(guasti),
            )
        tabella = contesto.albero.scrivi_testo(self.cartella, NOME_COPPIE, self._tsv_coppie(coppie))
        if appaiata:
            return self._appaiata(
                contesto, campioni, reverse, coppie, coppie_reverse, sincronie, classi, vuoti,
                tabella, processi)
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
        # Le letture inverse prese per forward: la quota si giudica dove ci si
        # attende il bersaglio, la struttura a blocchi in ogni file.
        # Un file il cui nome finisce come il secondo di una coppia e' giudicato
        # sulla quota da un numero minore di letture: il nome e' gia' un indizio.
        nomi_dei_file = [Path(p).name for p in campioni.values()]
        inattesi = {
            c: motivo for c, s in coppie.items()
            if (motivo := s.orientamento_inatteso(
                classi[c] in CLASSI_CONTROLLATE,
                LETTURE_MINIME_PER_UN_NOME_SOSPETTO
                if nome_da_seconda_lettura(Path(campioni[c]).name, nomi_dei_file)
                else LETTURE_MINIME_PER_LA_QUOTA))
        }
        if inattesi:
            elenco = "; ".join(f"{c} ({coppie[c].nome}): {motivo}"
                               for c, motivo in list(inattesi.items())[:5])
            raise errore(
                "E-S1-07",
                f"{len(inattesi)} file su {len(coppie)}, su tutte le loro letture, non "
                f"contengono le sole letture forward: {elenco}"
                + ("" if len(inattesi) <= 5 else f"; e altri {len(inattesi) - 5}")
                + ". Il segnale di una lettura forward e' "
                + _descrivi_il_segnale(contesto.config),
                campioni=sorted(inattesi),
            )
        esito, metriche = self._profilo_e_lunghezze(
            contesto, campioni, coppie, classi, vuoti, processi)
        return Produzione((tabella, *esito), metriche)

    @staticmethod
    def _taglio(config: Any, reverse: bool) -> dict[str, Any]:
        """Con il primer cercato lungo la lettura, il parametro con cui lo
        script profila le letture dopo il taglio: e' su quelle che il
        troncamento si misura.
        """
        ricerca = ricerca_del_primer(config, reverse)
        if ricerca is None:
            return {}
        return {"taglio_primer": {
            "primer": ricerca.primer,
            "discordanze_massime": ricerca.discordanze_massime,
            "sfasamento_massimo": ricerca.sfasamento_massimo,
            "scarta": ricerca.scarta_senza_primer,
        }}

    def _profilo_e_lunghezze(
        self,
        contesto: StepContext,
        campioni: dict[str, str],
        coppie: dict[str, StatisticheFile],
        classi: dict[str, ClasseCampione],
        vuoti: list[str],
        processi: int,
        *,
        reverse: bool = False,
    ) -> tuple[tuple[Any, ...], dict[str, Any]]:
        """Il profilo in R dei file dati e i controlli sulle lunghezze contro il
        troncamento: gli artefatti e le metriche. Con ``reverse`` i file sono
        quelli della lettura reverse, i parametri quelli del gruppo ``pair``, e
        gli artefatti portano il suffisso della reverse.
        """
        config = contesto.config
        if reverse:
            troncamento = config.pair.truncLen_reverse
            tolleranza = config.pair.truncLen_shortfall_warn_reverse
            parametro, parametro_tolleranza = (
                "pair.truncLen_reverse", "pair.truncLen_shortfall_warn_reverse")
            quali = " reverse"
        else:
            troncamento = config.filter.truncLen
            tolleranza = config.filter.truncLen_shortfall_warn
            parametro, parametro_tolleranza = "filter.truncLen", "filter.truncLen_shortfall_warn"
            quali = ""

        def nome(file: str) -> str:
            """Il nome del file, nella variante della lettura reverse quando si giudica
            quella."""
            return nome_reverse(file) if reverse else file

        # Le righe vuote in fondo a un file sono tollerate, ma il lettore a
        # blocchi di R le rifiuta: quei soli file gli arrivano in copia ripulita.
        vuote = [c for c, s in coppie.items() if s.righe_vuote_finali]
        if vuote:
            contesto.logger.info(
                f"{len(vuote)} file{quali} hanno righe vuote in fondo: tollerate, lette da una "
                "copia senza quelle righe, scritta nella cartella della fase e tolta alla "
                "fine",
                extra={"passo": "S1", "campioni": sorted(vuote)},
            )
        with senza_righe_vuote_finali(
            campioni, vuote, contesto.albero.cartella(self.cartella)
        ) as leggibili:
            esito = esegui_script(
                cartella_r() / NOME_SCRIPT,
                {
                    # I file senza letture non arrivano al lettore di R: il loro
                    # profilo e' vuoto per definizione, e lo script li conta a zero.
                    "campioni": {c: p for c, p in leggibili.items() if c not in vuoti},
                    "vuoti": vuoti,
                    "ordine": list(campioni),
                    "processi": processi,
                    **({"suffisso": SUFFISSO_REVERSE} if reverse else {}),
                    **self._taglio(config, reverse),
                },
                contesto.albero,
                self.cartella,
                passo=f"{self.passo}{SUFFISSO_REVERSE}" if reverse else self.passo,
                # Un processo per file, ognuno con un blocco di letture e le
                # loro qualita' in memoria: con molti processi e poca memoria
                # il sistema ne uccide alcuni.
                codice_memoria="E-S1-06",
                tempo_massimo_s=contesto.config.run.r_timeout_s,
                logger=contesto.logger,
            )

        cartella = contesto.albero.cartella(self.cartella)
        riepilogo = json.loads((cartella / nome(NOME_RIEPILOGO)).read_text(encoding="utf-8"))
        minimo = int(riepilogo["lunghezza_minima"])

        # Le letture piu' corte del troncamento, che il filtro scartera': una
        # sola non e' un motivo per fermarsi, una quota rilevante dei campioni
        # da cui ci si attende il segnale si'.
        lunghezze = self._lunghezze(cartella, nome(NOME_LUNGHEZZE))
        corte = letture_corte(lunghezze, classi, troncamento)
        massima = contesto.config.qc.max_frac_short_reads
        oltre = classi_oltre_soglia(corte, massima)

        with open(cartella / nome(NOME_VALORI_QUALITA), encoding="utf-8", newline="") as file:
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
                self._dettaglio_corte(lunghezze, classi, corte, oltre, troncamento, minimo,
                                      massima, parametro, quali),
                truncLen=troncamento,
                lunghezza_minima=minimo,
                frazione={classe: corte[classe]["frazione"] for classe in oltre},
            )

        scarto = minimo - troncamento
        if scarto > tolleranza:
            contesto.degrada(
                "E-S1-01",
                f"la lettura{quali} piu' corta, su tutte le letture, e' di {minimo} bp e "
                f"{parametro} vale {troncamento}: uno scarto di {scarto} bp, oltre "
                f"i {tolleranza} di {parametro_tolleranza}",
                lunghezza_minima=minimo,
                truncLen=troncamento,
                scarto=scarto,
            )
        if not reverse and len(valori) <= QUALITA_RAGGRUPPATE:
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

        return esito.artefatti, metriche

    def _appaiata(
        self,
        contesto: StepContext,
        campioni: dict[str, str],
        reverse: dict[str, str],
        coppie: dict[str, StatisticheFile],
        coppie_reverse: dict[str, StatisticheFile],
        sincronie: dict[str, Any],
        classi: dict[str, ClasseCampione],
        vuoti: list[str],
        tabella: Any,
        processi: int,
    ) -> Produzione:
        """Il seguito della fase in modalita' appaiata: i giudizi sulle coppie
        di file, poi il profilo delle due letture.
        """
        albero = contesto.albero
        tabella_reverse = albero.scrivi_testo(
            self.cartella, nome_reverse(NOME_COPPIE), self._tsv_coppie(coppie_reverse))

        def nomi(campione: str) -> str:
            """I nomi dei due file della coppia di un campione, in una frase."""
            return f"{coppie[campione].nome} e {coppie_reverse[campione].nome}"

        def elenca(voci: dict[str, str], quanti: int = 5) -> str:
            """Le prime voci come «campione: motivo», con il numero di quelle taciute."""
            testo = "; ".join(f"{c}: {m}" for c, m in list(voci.items())[:quanti])
            return testo + ("" if len(voci) <= quanti else f"; e altri {len(voci) - quanti}")

        # Le intestazioni: ogni file deve riportare una sola lettura di ogni
        # coppia, la prima nel forward e la seconda nel reverse.
        intercalati: dict[str, str] = {}
        invertite: dict[str, str] = {}
        for campione in reverse:
            for statistiche, attesa in ((coppie[campione], "1"), (coppie_reverse[campione], "2")):
                esito = statistiche.segni_di_coppia_inattesi(attesa)
                if esito is None:
                    continue
                specie, descrizione = esito
                (intercalati if specie == "intercalato" else invertite).setdefault(
                    campione, f"{statistiche.nome}: {descrizione}")
        if intercalati:
            raise errore(
                "E-S1-11",
                f"{len(intercalati)} file, su tutte le loro letture, contengono le due letture "
                f"di ogni coppia: {elenca(intercalati)}. In modalita' appaiata ogni file riporta "
                f"una sola lettura di ogni coppia; i conteggi sono in {NOME_COPPIE} e in "
                f"{nome_reverse(NOME_COPPIE)}",
                campioni=sorted(intercalati),
            )

        fuori_posto = {c: s.difetto for c, s in sincronie.items() if s.difetto is not None}
        if fuori_posto:
            raise errore(
                "E-S1-08",
                f"{len(fuori_posto)} campioni su {len(sincronie)} hanno i due file non in "
                "sincronia: " + elenca({f"{c} ({nomi(c)})": m for c, m in fuori_posto.items()}),
                campioni=sorted(fuori_posto),
            )

        # L'orientamento, dal contenuto. Prima la corsa invertita, che spiega
        # insieme il difetto dei suoi due file; poi ciascun file per se'.
        for campione in reverse:
            if campione not in invertite:
                motivo = corsa_invertita(coppie[campione], coppie_reverse[campione])
                if motivo is not None:
                    invertite[campione] = motivo
        senza_orientamento: dict[str, str] = {}
        inattesi: dict[str, str] = {}
        mescolati: dict[str, str] = {}
        for campione in reverse:
            if campione in invertite:
                continue
            controllata = classi[campione] in CLASSI_CONTROLLATE
            if motivo := coppie[campione].orientamento_inatteso(controllata):
                inattesi[campione] = motivo
            if motivo := coppie_reverse[campione].orientamento_inatteso(controllata):
                senza_orientamento[campione] = f"{coppie_reverse[campione].nome}: {motivo}"
            for statistiche in (coppie[campione], coppie_reverse[campione]):
                if motivo := versi_mescolati(statistiche):
                    mescolati.setdefault(campione, f"{statistiche.nome}: {motivo}")
        if invertite or senza_orientamento:
            parti = []
            if invertite:
                parti.append(
                    f"{len(invertite)} corse hanno i due file scambiati: "
                    + elenca({f"{c} ({nomi(c)})": m for c, m in invertite.items()})
                    + ". Scambia i nomi dei due file di ciascuna")
            if senza_orientamento:
                parti.append(
                    f"{len(senza_orientamento)} file reverse non contengono le sole letture "
                    f"reverse: {elenca(senza_orientamento)}")
            raise errore(
                "E-S1-09", ". ".join(parti),
                campioni=sorted({*invertite, *senza_orientamento}),
                invertite=sorted(invertite),
            )
        if inattesi:
            raise errore(
                "E-S1-07",
                f"{len(inattesi)} file forward su {len(coppie)}, su tutte le loro letture, non "
                "contengono le sole letture forward: "
                + elenca({f"{c} ({coppie[c].nome})": m for c, m in inattesi.items()})
                + ". Il segnale di una lettura forward e' "
                + _descrivi_il_segnale(contesto.config),
                campioni=sorted(inattesi),
            )
        if mescolati:
            raise errore(
                "E-S1-10",
                f"{len(mescolati)} campioni hanno un file con le letture dei due versi: "
                + elenca(mescolati),
                campioni=sorted(mescolati),
            )

        artefatti, metriche = self._profilo_e_lunghezze(
            contesto, campioni, coppie, classi, vuoti, processi)
        vuoti_reverse = [c for c, s in coppie_reverse.items() if s.vuoto]
        artefatti_reverse, metriche_reverse = self._profilo_e_lunghezze(
            contesto, reverse, coppie_reverse, classi, vuoti_reverse, processi, reverse=True)
        metriche["reverse"] = metriche_reverse
        metriche["coppie_in_sincronia"] = sum(s.coppie for s in sincronie.values())
        return Produzione(
            (tabella, tabella_reverse, *artefatti, *artefatti_reverse), metriche)

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
    def _lunghezze(cartella: Any, nome: str = NOME_LUNGHEZZE) -> dict[str, dict[int, int]]:
        """Per campione, le letture di ogni lunghezza (da ``lunghezze.tsv``)."""
        lunghezze: dict[str, dict[int, int]] = defaultdict(dict)
        with open(cartella / nome, encoding="utf-8", newline="") as file:
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
        parametro: str = "filter.truncLen",
        quali: str = "",
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
            f"{parametro} vale {troncamento} e {superano} e' piu' corto, oltre il "
            f"{massima:.1%} di qc.max_frac_short_reads; la lettura{quali} piu' corta e' di {minimo} "
            f"bp. Per classe: {per_classe}. {len(colpiti)} campioni di "
            f"{'questa classe' if len(oltre) == 1 else 'queste classi'} hanno letture piu' "
            f"corte, che verrebbero scartate; i piu' colpiti: {elenco}"
        )
