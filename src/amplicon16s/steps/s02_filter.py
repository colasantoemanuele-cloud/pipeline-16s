"""Fase S2: filtro e troncamento delle letture.

È la prima trasformazione irreversibile del dato: tronca le letture a
``filter.truncLen`` e scarta quelle che non superano il filtro. Tutte le fasi
successive lavorano sulle sue uscite, in ``03_filtered/``
(``R/02_filter.R``, ``dada2::filterAndTrim``, attraverso il ponte).

**Gli archivi si verificano per intero prima del filtro.** Un archivio gzip
troncato non fa fallire ``filterAndTrim``: viene letto fino a dove arriva, e
le letture perse spariscono senza traccia (verificato su un file della
versione ridotta tagliato a meta'). S2 quindi decomprime ogni file per intero, in parallelo, e un archivio incompleto e'
E-S2-03, errore di lettura: se era transitorio la rilettura riesce, se il file
e' davvero corrotto la fase si ferma dopo i tentativi ammessi.

**I controlli sul risultato valgono per i biologici e i positivi.** Come per
G10, un controllo si applica alle classi per cui ha senso:

* E-S2-01, campioni azzerati oltre ``qc.max_zeroed_samples``: un campione
  biologico azzerato e' un campione perso, un controllo positivo azzerato mette
  in dubbio il lotto. Un controllo negativo azzerato no: e' un bianco pulito, e
  fermare la pipeline per quello sarebbe sbagliato. Si registra e si prosegue.
* E-S2-02, perdita media oltre ``qc.max_frac_lost_filter``: domanda se i
  parametri del filtro sono adatti ai campioni analizzati. La perdita di un
  bianco non risponde a quella domanda: i bianchi hanno poca biomassa, e la
  loro qualita' puo' essere peggiore senza che i parametri siano sbagliati.

L'esclusione dei negativi e' di principio, non una misura: vale anche quando i
negativi conservano quanto i biologici. Le classi sono in
:data:`~amplicon16s.metadata.models.CLASSI_CONTROLLATE`. La media e' per campione, come la chiede il controllo: un campione povero
pesa quanto uno profondo.
"""

from __future__ import annotations

import csv
import gzip
import zlib
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Any, ClassVar, Final

from amplicon16s.config import defaults
from amplicon16s.config.schema import Config, Qc, thread_effettivi
from amplicon16s.errors.exceptions import ErrorePipeline, errore
from amplicon16s.io_layer.artifacts import Fase, ManifestoPasso
from amplicon16s.io_layer.conteggi import leggi_conteggi
from amplicon16s.io_layer.reads import senza_righe_vuote_finali
from amplicon16s.metadata.models import CLASSI_CONTROLLATE, ClasseCampione
from amplicon16s.rbridge.runner import cartella_r, esegui_script
from amplicon16s.runner.graph import Passo
from amplicon16s.runner.retry import Aggiustamento, senza_modifiche
from amplicon16s.steps.base import PipelineStep, Produzione, StepContext
from amplicon16s.steps.s01_profile import NOME_COPPIE, SUFFISSO_REVERSE

__all__ = [
    "FiltroLetture",
    "controlla_filtro",
    "esamina_archivio",
]

NOME_SCRIPT: Final = "02_filter.R"

NOME_PREFILTRO: Final = "letture_prefiltro.tsv"
NOME_FILTRATE: Final = "letture_filtrate.tsv"
#: Solo con il primer cercato lungo la lettura: per campione, le letture (o le
#: coppie) in ingresso, quelle con il primer, quelle senza, quelle rimaste.
NOME_TAGLIO: Final = "taglio_primer.tsv"
#: La cartella di lavoro delle letture tagliate: porta il prefisso delle
#: scritture in corso, e a fine fase non esiste piu'.
CARTELLA_DEL_TAGLIO: Final = ".scrittura-taglio"
#: Le letture filtrate di un campione: ``<accession>_filt.fastq.gz``.
SUFFISSO_FILTRATI: Final = "_filt.fastq.gz"
#: In modalita' appaiata, le letture reverse filtrate dello stesso campione.
SUFFISSO_FILTRATI_REVERSE: Final = "_filt_reverse.fastq.gz"


def esamina_archivio(percorso: Path) -> tuple[str | None, bool, bool]:
    """Il motivo per cui un archivio non si decomprime per intero (o ``None``),
    se il file si chiude con righe vuote, e se non ha alcun contenuto.

    Un FASTQ non compresso, riconosciuto dai primi byte come in S0, non ha
    nulla da decomprimere: basta che si legga fino in fondo. Le righe vuote in
    fondo si riconoscono dalla coda del contenuto: sono tollerate, ma
    ``filterAndTrim`` non le ammette, e quel file gli arriva in copia ripulita
    (``io_layer/reads.py``).
    """
    coda = b""
    con_contenuto = False
    try:
        with open(percorso, "rb") as file:
            compresso = file.read(2) == b"\x1f\x8b"
        with (gzip.open if compresso else open)(percorso, "rb") as file:
            while blocco := file.read(1 << 22):
                con_contenuto = con_contenuto or bool(blocco.strip())
                coda = (coda + blocco)[-4096:]
    except (OSError, EOFError, zlib.error) as e:
        return f"{type(e).__name__}: {e}", False, False
    if not con_contenuto:
        return None, False, True
    # Dopo l'ultima riga con del contenuto, piu' di un a capo vuol dire almeno
    # una riga vuota.
    return None, coda[len(coda.rstrip(b" \t\r\n")):].count(b"\n") > 1, False


def controlla_filtro(
    ingresso: dict[str, int],
    uscita: dict[str, int],
    classi: dict[str, ClasseCampione],
    qc: Qc,
) -> dict[str, Any]:
    """Le metriche del filtro per classe, e i controlli E-S2-01 ed E-S2-02.

    I controlli guardano solo le classi di :data:`CLASSI_CONTROLLATE`; un
    controllo negativo azzerato compare fra le metriche e non ferma nulla.
    """
    azzerati: dict[str, list[str]] = {classe.value: [] for classe in ClasseCampione}
    frazioni: dict[str, list[float]] = {classe.value: [] for classe in ClasseCampione}
    for campione, letture in ingresso.items():
        classe = classi[campione].value
        if uscita.get(campione, 0) == 0:
            azzerati[classe].append(campione)
        if letture > 0:
            frazioni[classe].append(uscita.get(campione, 0) / letture)
    metriche: dict[str, Any] = {
        "campioni": len(ingresso),
        "letture_ingresso": sum(ingresso.values()),
        "letture_uscita": sum(uscita.values()),
        "azzerati": {k: sorted(v) for k, v in azzerati.items()},
        "frazione_conservata_media": {
            k: round(sum(v) / len(v), 4) for k, v in frazioni.items() if v
        },
    }

    controllati = sorted(c for classe in CLASSI_CONTROLLATE for c in azzerati[classe.value])
    if len(controllati) > qc.max_zeroed_samples:
        raise errore(
            "E-S2-01",
            f"{len(controllati)} campioni biologici o controlli positivi senza "
            f"letture dopo il filtro, oltre i {qc.max_zeroed_samples} ammessi da "
            f"qc.max_zeroed_samples: {', '.join(controllati)}",
            campioni=controllati,
        )

    considerate = [f for classe in CLASSI_CONTROLLATE for f in frazioni[classe.value]]
    if considerate:
        perdita = 1 - sum(considerate) / len(considerate)
        metriche["perdita_media_controllata"] = round(perdita, 4)
        if perdita > qc.max_frac_lost_filter:
            raise errore(
                "E-S2-02",
                f"il filtro scarta in media il {perdita:.1%} delle letture dei "
                f"campioni biologici e dei controlli positivi, oltre il "
                f"{qc.max_frac_lost_filter:.0%} di qc.max_frac_lost_filter",
                perdita_media=round(perdita, 4),
            )
    return metriche


def taglio_del_primer(config: Any, cartella: Path) -> dict[str, Any] | None:
    """I parametri del taglio dei primer cercati per ``R/02_filter.R``;
    ``None`` con il taglio a posizione fissa.
    """
    if not config.taglio_cercato:
        return None
    filtro = config.filter
    taglio = {
        "primer": config.qc.primer_sequence,
        "discordanze_massime": filtro.primer_max_mismatch,
        "sfasamento_massimo": filtro.primer_max_offset,
        "scarta": filtro.primer_untrimmed == defaults.SENZA_PRIMER_SCARTATA,
        "lavoro": str(cartella / CARTELLA_DEL_TAGLIO),
    }
    if config.appaiata:
        taglio["primer_reverse"] = config.pair.primer_sequence_reverse
    return taglio


def controlla_taglio(
    percorso: Path, classi: dict[str, ClasseCampione], massima: float
) -> tuple[dict[str, Any], str | None]:
    """Le metriche del taglio dei primer cercati e, se la quota di letture
    senza primer supera ``qc.max_frac_untrimmed``, il testo dell'avviso
    (``E-S2-04``).

    La quota e' sull'insieme delle letture dei campioni biologici e dei
    controlli positivi (in modalita' appaiata: le coppie in cui manca almeno un
    primer): i controlli negativi amplificano poco, e le loro letture senza
    primer non dicono nulla del protocollo. Il confronto e' stretto.
    """
    with open(percorso, encoding="utf-8", newline="") as file:
        righe = list(csv.DictReader(file, delimiter="\t"))
    controllate = [r for r in righe if classi[r["campione"]] in CLASSI_CONTROLLATE]
    letture = sum(int(r["letture"]) for r in controllate)
    senza = sum(int(r["senza_primer"]) for r in controllate)
    metriche: dict[str, Any] = {
        "letture": sum(int(r["letture"]) for r in righe),
        "senza_primer": sum(int(r["senza_primer"]) for r in righe),
        "rimaste": sum(int(r["rimaste"]) for r in righe),
        "letture_controllate": letture,
        "senza_primer_controllate": senza,
    }
    if not letture:
        return metriche, None
    metriche["quota_senza_primer_controllata"] = round(senza / letture, 6)
    if senza > massima * letture:
        return metriche, (
            f"il primer non si trova nel {senza / letture:.1%} delle letture dei campioni "
            f"biologici e dei controlli positivi ({senza} su {letture}), oltre il "
            f"{massima:.1%} di qc.max_frac_untrimmed")
    return metriche, None


def _letture_verificate(profili: Path) -> dict[str, tuple[int, int]]:
    """Per campione, le letture che S1 ha contato nel file forward e nel
    reverse dopo averne verificato la sincronia.
    """
    conteggi: list[dict[str, int]] = []
    for nome in (NOME_COPPIE, NOME_COPPIE.replace(".tsv", f"{SUFFISSO_REVERSE}.tsv")):
        with open(profili / nome, encoding="utf-8", newline="") as file:
            conteggi.append(
                {r["campione"]: int(r["letture"]) for r in csv.DictReader(file, delimiter="\t")})
    forward, reverse = conteggi
    return {a: (forward[a], reverse[a]) for a in forward if a in reverse}


class FiltroLetture(PipelineStep):
    """S2: filtro e troncamento."""

    passo: ClassVar[Passo] = Passo.S2
    #: 2: riconosce dai primi byte un FASTQ non compresso, che non ha un
    #: archivio da verificare. 3: un file che si chiude con righe vuote si
    #: filtra da una copia senza quelle righe. 4: un file senza letture (il
    #: bianco che S0 ha ammesso) non arriva al filtro: zero letture in ingresso
    #: e in uscita.
    versione: ClassVar[int] = 4
    script_r: ClassVar[str | None] = NOME_SCRIPT
    passi_tracciamento: ClassVar[tuple[str, ...]] = ("prefiltro", "filtrate")
    #: Con il primer cercato, fra l'ingresso e il filtro: le letture che
    #: restano dopo il taglio.
    passi_tracciamento_di_modalita: ClassVar[tuple[tuple[str, str], ...]] = (
        ("taglio_primer", "prefiltro"),)
    #: Tutto il gruppo filter e le due soglie dei controlli sul risultato.
    #: run.batch_size, i file passati insieme al filtro, non vi entra: ogni file
    #: si filtra da solo. Variare filter.maxEE o qc.max_frac_lost_filter
    #: ricalcola S2 e le fasi a valle (S3+), ma preserva interamente S0 e S1 gia'
    #: concluse.
    parametri: ClassVar[tuple[str, ...]] = (
        "filter", "qc.max_zeroed_samples", "qc.max_frac_lost_filter",
    )
    #: Troncamento, taglio ed errori attesi della lettura reverse.
    parametri_appaiata: ClassVar[tuple[str, ...]] = (
        "pair.truncLen_reverse", "pair.trimLeft_reverse", "pair.maxEE_reverse",
    )
    #: Con il primer cercato: il primer (quello della reverse in modalita'
    #: appaiata) e la soglia sulla quota di letture senza primer; modo,
    #: discordanze, sfasamento e destino delle letture senza primer sono nel
    #: gruppo filter.
    parametri_taglio_cercato: ClassVar[tuple[str, ...]] = (
        "qc.primer_sequence", "pair.primer_sequence_reverse", "qc.max_frac_untrimmed",
    )
    aggiustamenti: ClassVar[dict[str, Aggiustamento]] = {
        "E-S2-03": senza_modifiche(
            "un errore di lettura transitorio si corregge rileggendo lo stesso file"
        ),
    }

    def artefatti_temporanei(self, manifesto: ManifestoPasso, config: Config) -> tuple[str, ...]:
        """Le letture filtrate, se ``run.keep_filtered_fastq`` e' falso."""
        if config.run.keep_filtered_fastq:
            return ()
        return tuple(
            n for n in manifesto.nomi
            if n.endswith(SUFFISSO_FILTRATI) or n.endswith(SUFFISSO_FILTRATI_REVERSE))

    def calcola(self, contesto: StepContext) -> Produzione:
        """Verifica che gli archivi si decomprimano per intero, esegue ``R/02_filter.R``
        e controlla le letture rimaste per classe.
        """
        if contesto.inventario is None:
            raise RuntimeError("S2 richiede l'inventario prodotto da S0")
        config = contesto.config
        filtro = config.filter
        campioni = {
            c.accession: str(c.file) for c in contesto.inventario if c.file is not None
        }
        reverse = {
            c.accession: str(c.file_reverse) for c in contesto.inventario
            if c.file is not None and c.file_reverse is not None
        } if config.appaiata else {}
        with ThreadPoolExecutor(max_workers=thread_effettivi(config)) as esecutore:
            esami = dict(zip(campioni, esecutore.map(esamina_archivio, map(Path, campioni.values()))))
            esami_reverse = dict(zip(
                reverse, esecutore.map(esamina_archivio, map(Path, reverse.values()))))
        incompleti = {a: m for a, (m, _, _) in esami.items() if m is not None}
        incompleti_reverse = {a: m for a, (m, _, _) in esami_reverse.items() if m is not None}
        if incompleti or incompleti_reverse:
            raise errore(
                "E-S2-03",
                "archivi che non si decomprimono per intero: " + "; ".join(
                    [f"{a} ({Path(campioni[a]).name}: {m})" for a, m in sorted(incompleti.items())]
                    + [f"{a} ({Path(reverse[a]).name}: {m})"
                       for a, m in sorted(incompleti_reverse.items())]
                ),
                campioni=sorted({*incompleti, *incompleti_reverse}),
            )
        if config.appaiata:
            return self._appaiata(contesto, campioni, reverse, esami, esami_reverse)
        vuote = [a for a, (_, righe_vuote, _) in esami.items() if righe_vuote]
        # Un file senza letture non ha nulla da filtrare: e' il bianco che S0 ha
        # ammesso, e risulta azzerato come ogni bianco pulito.
        senza_letture = [a for a, (_, _, vuoto) in esami.items() if vuoto]
        with senza_righe_vuote_finali(
            campioni, vuote, contesto.albero.cartella(self.cartella)
        ) as leggibili:
            esito = esegui_script(
                cartella_r() / NOME_SCRIPT,
                {
                    "campioni": {a: p for a, p in leggibili.items() if a not in senza_letture},
                    "vuoti": senza_letture,
                    "ordine": list(campioni),
                    "truncLen": filtro.truncLen,
                    "trimLeft": filtro.trimLeft,
                    "maxEE": filtro.maxEE,
                    "truncQ": filtro.truncQ,
                    "maxN": filtro.maxN,
                    "rm_phix": filtro.rm_phix,
                    "processi": thread_effettivi(config),
                    "lotto": config.run.batch_size,
                    **self._taglio(contesto),
                },
                contesto.albero,
                self.cartella,
                passo=self.passo,
                tempo_massimo_s=contesto.config.run.r_timeout_s,
                logger=contesto.logger,
            )
        metriche = self._controlla(contesto)
        self._giudica_il_taglio(contesto, metriche)
        negativi = metriche["azzerati"][ClasseCampione.CONTROLLO_NEGATIVO.value]
        if negativi:
            contesto.logger.info(
                f"{len(negativi)} controlli negativi senza letture dopo il filtro: bianchi puliti",
                extra={"passo": "S2", "campioni": negativi},
            )
        return Produzione(esito.artefatti, metriche)

    def _taglio(self, contesto: StepContext) -> dict[str, Any]:
        """Il parametro ``taglio_primer`` dello script, solo con il primer cercato."""
        taglio = taglio_del_primer(contesto.config, contesto.albero.cartella(self.cartella))
        return {} if taglio is None else {"taglio_primer": taglio}

    def _controlla(self, contesto: StepContext) -> dict[str, Any]:
        """I controlli sul risultato del filtro. Con il primer cercato un
        arresto dice anche quante letture ha tolto il taglio del primer: la
        perdita e' quella dell'intera fase, e la causa puo' essere li'.
        """
        assert contesto.inventario is not None
        cartella = contesto.albero.cartella(self.cartella)
        classi = {c.accession: c.classe for c in contesto.inventario}
        try:
            return controlla_filtro(
                leggi_conteggi(cartella / NOME_PREFILTRO), leggi_conteggi(cartella / NOME_FILTRATE),
                classi, contesto.config.qc)
        except ErrorePipeline as guasto:
            if not contesto.config.taglio_cercato:
                raise
            del_taglio, _ = controlla_taglio(cartella / NOME_TAGLIO, classi, 1.0)
            raise errore(
                guasto.codice,
                f"{guasto.dettaglio}. Con filter.primer_removal search la perdita comprende le "
                f"letture senza il primer: {del_taglio['senza_primer_controllate']} su "
                f"{del_taglio['letture_controllate']} nei campioni biologici e nei controlli "
                f"positivi ({NOME_TAGLIO}, posizioni_primer.tsv). Se sono molte, verifica "
                "qc.primer_sequence e se filter.primer_max_offset copre le basi che "
                "precedono il primer",
                **guasto.contesto) from guasto

    def _giudica_il_taglio(self, contesto: StepContext, metriche: dict[str, Any]) -> None:
        """Con il primer cercato aggiunge alle metriche quelle del taglio, e
        registra l'avviso se troppe letture non portano il primer.
        """
        config = contesto.config
        if not config.taglio_cercato:
            return
        assert contesto.inventario is not None
        del_taglio, avviso = controlla_taglio(
            contesto.albero.cartella(self.cartella) / NOME_TAGLIO,
            {c.accession: c.classe for c in contesto.inventario},
            config.qc.max_frac_untrimmed)
        metriche["taglio_primer"] = del_taglio
        if avviso is not None:
            contesto.degrada(
                "E-S2-04", avviso, quota=del_taglio["quota_senza_primer_controllata"])

    def _appaiata(
        self,
        contesto: StepContext,
        campioni: dict[str, str],
        reverse: dict[str, str],
        esami: dict[str, tuple[str | None, bool, bool]],
        esami_reverse: dict[str, tuple[str | None, bool, bool]],
    ) -> Produzione:
        """Il filtro in modalita' appaiata: le due letture di ogni coppia si
        tengono o si scartano insieme, ciascuna con il suo troncamento, il suo
        taglio e i suoi errori attesi. I conteggi sono coppie.
        """
        config = contesto.config
        filtro, coppie = config.filter, config.pair
        # Il filtro appaiato non confronta gli identificativi delle letture: le
        # prende a coppie per posizione. Che i due file di ogni campione siano
        # in sincronia lo ha verificato S1 su tutte le letture, e qui si legge
        # il suo esito: senza, il filtro non parte.
        verificate = _letture_verificate(contesto.albero.cartella(Fase.QC_PROFILES))
        non_verificati = sorted(a for a in reverse if a not in verificate)
        if non_verificati:
            raise RuntimeError(
                "S2 in modalita' appaiata richiede la verifica della sincronia delle coppie "
                f"fatta da S1, che manca per: {', '.join(non_verificati[:5])}")
        # Un campione senza letture in uno dei due file non ha coppie da filtrare.
        senza_letture = sorted(
            {a for a, (_, _, vuoto) in esami.items() if vuoto}
            | {a for a, (_, _, vuoto) in esami_reverse.items() if vuoto})
        vuote = [a for a, (_, righe_vuote, _) in esami.items() if righe_vuote]
        vuote_reverse = [a for a, (_, righe_vuote, _) in esami_reverse.items() if righe_vuote]
        cartella = contesto.albero.cartella(self.cartella)
        # Le copie ripulite dei due versi stanno in due cartelle: la seconda
        # chiamata non deve togliere quelle della prima.
        with senza_righe_vuote_finali(campioni, vuote, cartella) as leggibili, \
                senza_righe_vuote_finali(
                    reverse, vuote_reverse, cartella, ".scrittura-letture-reverse",
                ) as leggibili_reverse:
            esito = esegui_script(
                cartella_r() / NOME_SCRIPT,
                {
                    "campioni": {a: p for a, p in leggibili.items() if a not in senza_letture},
                    "reverse": {
                        a: p for a, p in leggibili_reverse.items() if a not in senza_letture},
                    "vuoti": senza_letture,
                    "ordine": list(campioni),
                    "truncLen": filtro.truncLen,
                    "trimLeft": filtro.trimLeft,
                    "maxEE": filtro.maxEE,
                    "truncLen_reverse": coppie.truncLen_reverse,
                    "trimLeft_reverse": coppie.trimLeft_reverse,
                    "maxEE_reverse": coppie.maxEE_reverse,
                    "truncQ": filtro.truncQ,
                    "maxN": filtro.maxN,
                    "rm_phix": filtro.rm_phix,
                    "processi": thread_effettivi(config),
                    "lotto": config.run.batch_size,
                    **self._taglio(contesto),
                },
                contesto.albero,
                self.cartella,
                passo=self.passo,
                tempo_massimo_s=contesto.config.run.r_timeout_s,
                logger=contesto.logger,
            )
        metriche = self._controlla(contesto)
        self._giudica_il_taglio(contesto, metriche)
        negativi = metriche["azzerati"][ClasseCampione.CONTROLLO_NEGATIVO.value]
        if negativi:
            contesto.logger.info(
                f"{len(negativi)} controlli negativi senza coppie dopo il filtro: bianchi puliti",
                extra={"passo": "S2", "campioni": negativi},
            )
        return Produzione(esito.artefatti, metriche)
