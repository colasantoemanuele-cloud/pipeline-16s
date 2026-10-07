"""Fase S3: il modello d'errore di sequenziamento, uno per corsa.

Stima dai dati stessi il profilo degli errori di sequenziamento, su cui S4
distinguera' le varianti vere dagli errori (``R/03_learn_errors.R``,
``dada2::learnErrors``, attraverso il ponte). Un modello sbagliato non fa
fallire nulla: produce varianti spurie o ne perde di vere.

**Un modello per corsa.** Corse diverse hanno profili d'errore diversi, e un
modello unico li medierebbe perdendo accuratezza su entrambe. La corsa di un
campione viene da ``err.batch_column`` in ``io.batch_table``; senza l'una o
l'altra si stima un solo modello su tutti i campioni, e la pipeline resta
utilizzabile senza il file di arricchimento del lotto. Con la colonna attiva,
un campione senza corsa non ha un modello: E-S3-02.

**Quali campioni entrano nella stima.** ``learnErrors`` usa i file nell'ordine
dato finche' le basi non superano ``err.nbases``. L'ordine lo decide questa
fase: i campioni della corsa, mescolati con il seme ``run.seed`` se
``err.randomize`` e' vero, e passati a R con ``randomize = FALSE``. Cosi' la
scelta e' riproducibile e si sa esattamente quali campioni, e quante basi di
ciascuna classe, hanno fatto la stima.

**La corrispondenza campione-modello e' un artefatto**, ``corrispondenza.tsv``:
S4 la legge invece di ricalcolarla. I nomi delle corse diventano parti di nomi
di file, e sono ridotti a lettere, cifre, punto, trattino e trattino basso.

**La convergenza e il suo retry.** Un modello che non converge entro
``err.max_consist`` iterazioni e' E-S3-01, ritentato con ``err.nbases``
raddoppiato. Il raddoppio aiuta solo se la stima non usava gia' tutte le basi
della corsa: in quel caso il nuovo tentativo sarebbe identico, e la fase lo
dichiara nell'errore perche' l'esecutore non lo tenti.
"""

from __future__ import annotations

import json
import random
import re
from collections import defaultdict
from dataclasses import dataclass
from typing import Any, ClassVar, Final

from amplicon16s.config.schema import Config, thread_effettivi
from amplicon16s.errors.exceptions import errore
from amplicon16s.io_layer.artifacts import Fase
from amplicon16s.io_layer.conteggi import leggi_conteggi
from amplicon16s.metadata.models import Inventario
from amplicon16s.rbridge.runner import cartella_r, esegui_script
from amplicon16s.runner.graph import Passo
from amplicon16s.runner.retry import RITENTARE_INUTILE, Aggiustamento, raddoppia
from amplicon16s.steps.base import PipelineStep, Produzione, StepContext
from amplicon16s.steps.s02_filter import NOME_FILTRATE, SUFFISSO_FILTRATI

__all__ = [
    "MODELLO_UNICO",
    "Modello",
    "ModelloErrore",
    "nomi_sicuri",
    "pianifica",
]

NOME_SCRIPT: Final = "03_learn_errors.R"
NOME_CORRISPONDENZA: Final = "corrispondenza.tsv"
NOME_MODELLI: Final = "modelli.json"
NOME_CONVERGENZA: Final = "convergenza.json"

#: Il nome del modello quando se ne stima uno solo.
MODELLO_UNICO: Final = "tutti"

_NON_SICURO: Final = re.compile(r"[^A-Za-z0-9._-]+")


def nomi_sicuri(corse: list[str]) -> dict[str, str]:
    """Un nome utilizzabile in un nome di file per ciascuna corsa, senza collisioni."""
    nomi: dict[str, str] = {}
    usati: set[str] = set()
    for corsa in sorted(corse):
        base = _NON_SICURO.sub("_", corsa).strip("._-") or "corsa"
        nome, indice = base, 2
        while nome in usati or nome == MODELLO_UNICO:
            nome, indice = f"{base}_{indice}", indice + 1
        usati.add(nome)
        nomi[corsa] = nome
    return nomi


@dataclass(frozen=True)
class Modello:
    """Il piano di un modello: quali campioni gli appartengono e quali lo stimano."""

    nome: str
    #: Il valore di err.batch_column; ``None`` per il modello unico.
    corsa: str | None
    #: I campioni della corsa, in ordine di accession.
    campioni: tuple[str, ...]
    #: I campioni con letture filtrate, nell'ordine in cui learnErrors li usa.
    ordine: tuple[str, ...]
    #: Quelli che la stima usa davvero: i primi, finche' le basi superano nbases.
    usati: tuple[str, ...]
    basi_usate: int
    basi_disponibili: int

    @property
    def tutte_usate(self) -> bool:
        """Se la stima usa tutte le basi della corsa: raddoppiare nbases non cambierebbe nulla."""
        return len(self.usati) == len(self.ordine)


def pianifica(
    inventario: Inventario,
    letture: dict[str, int],
    lunghezza: int,
    config: Config,
) -> list[Modello]:
    """I modelli da stimare, dalla corsa di ogni campione e dalle sue letture filtrate."""
    # Razionale biologico: ciascuna corsa di sequenziamento (run_prefix) presenta
    # un proprio profilo fisico di rumore ottico e di fasatura dei cicli Illumina.
    # Stimare un modello parametrico separato per corsa evita medie spurie tra
    # corse eterogenee, preservando in S4 la discriminazione statistica tra errori
    # di sostituzione e varianti ecologiche reali (ASV). Se la partizione per
    # corsa e' attiva ma un campione e' privo di corsa associata, si solleva
    # subito E-S3-02 per impedire assegnazioni arbitrarie silenziose.
    per_lotto = config.err.batch_column is not None and config.io.batch_table is not None
    if per_lotto:
        senza = sorted(c.accession for c in inventario if c.corsa is None)
        if senza:
            raise errore(
                "E-S3-02",
                f"{len(senza)} campioni senza {config.err.batch_column}: {', '.join(senza[:10])}"
                + (f" e altri {len(senza) - 10}" if len(senza) > 10 else ""),
                campioni=senza,
            )
        gruppi: dict[str | None, list[str]] = defaultdict(list)
        for campione in inventario:
            gruppi[campione.corsa].append(campione.accession)
        nomi = nomi_sicuri([c for c in gruppi if c is not None])
    else:
        gruppi = {None: [c.accession for c in inventario]}
        nomi = {}

    modelli = []
    for corsa in sorted(gruppi, key=lambda c: nomi.get(c, MODELLO_UNICO)):
        nome = nomi[corsa] if corsa is not None else MODELLO_UNICO
        campioni = tuple(sorted(gruppi[corsa]))
        utilizzabili = [c for c in campioni if letture.get(c, 0) > 0]
        if config.err.randomize:
            random.Random(f"{config.run.seed}:{nome}").shuffle(utilizzabili)
        usati, basi = [], 0
        for campione in utilizzabili:
            usati.append(campione)
            basi += letture[campione] * lunghezza
            # La stessa regola di learnErrors: si ferma dopo il file che
            # porta le basi oltre nbases.
            if basi > config.err.nbases:
                break
        modelli.append(
            Modello(
                nome=nome,
                corsa=corsa,
                campioni=campioni,
                ordine=tuple(utilizzabili),
                usati=tuple(usati),
                basi_usate=basi,
                basi_disponibili=sum(letture[c] * lunghezza for c in utilizzabili),
            )
        )
    return modelli


class ModelloErrore(PipelineStep):
    """S3: un modello d'errore per corsa."""

    passo: ClassVar[Passo] = Passo.S3
    #: 2: loess_monotono adatta un loess di primo grado con span 2 e pesi
    #: logaritmici (prima span 0,95); una stima fallita con piu' di un valore
    #: di qualita' ha il suo codice (E-S3-05), distinto da E-S3-04, e porta il
    #: messaggio originale della funzione di stima.
    versione: ClassVar[int] = 2
    script_r: ClassVar[str | None] = NOME_SCRIPT
    #: Tutto il gruppo err; il seme per l'ordine dei campioni; troncamento e
    #: taglio iniziale per le basi di ciascuna lettura; io.batch_table, che
    #: decide se si stima un modello per corsa. Variare err.nbases o
    #: err.max_consist ricalcola esclusivamente S3 (e S4+), lasciando validi i
    #: manifesti di S0, S1 e S2.
    parametri: ClassVar[tuple[str, ...]] = (
        "err", "run.seed", "filter.truncLen", "filter.trimLeft", "io.batch_table",
    )
    aggiustamenti: ClassVar[dict[str, Aggiustamento]] = {
        "E-S3-01": raddoppia("err.nbases"),
    }

    def calcola(self, contesto: StepContext) -> Produzione:
        """Pianifica e stima un modello d'errore per corsa con ``R/03_learn_errors.R`` e
        ne verifica la convergenza.
        """
        if contesto.inventario is None:
            raise RuntimeError("S3 richiede l'inventario prodotto da S0")
        config = contesto.config
        filtrati = contesto.albero.cartella(Fase.FILTERED)
        letture = leggi_conteggi(filtrati / NOME_FILTRATE)
        # Le letture filtrate hanno tutte questa lunghezza: S2 le tronca.
        lunghezza = config.filter.truncLen - config.filter.trimLeft
        modelli = pianifica(contesto.inventario, letture, lunghezza, config)
        vuoti = [m.nome for m in modelli if not m.ordine]
        if vuoti:
            # Una corsa i cui campioni non hanno letture dopo il filtro non ha
            # dati su cui stimare il proprio modello: e' una condizione dei
            # dati, e ha il suo codice.
            raise errore(
                "E-S3-03",
                f"nessuna lettura filtrata per {'la corsa' if len(vuoti) == 1 else 'le corse'} "
                f"{', '.join(vuoti)}",
                modelli=vuoti,
            )

        esito = esegui_script(
            cartella_r() / NOME_SCRIPT,
            {
                "modelli": {
                    m.nome: [str(filtrati / f"{c}{SUFFISSO_FILTRATI}") for c in m.ordine]
                    for m in modelli
                },
                "nbases": config.err.nbases,
                "max_consist": config.err.max_consist,
                "funzione_errore": config.err.error_function,
                "seme": config.run.seed,
                "processi": thread_effettivi(config),
            },
            contesto.albero,
            self.cartella,
            passo=self.passo,
            logger=contesto.logger,
        )
        cartella = contesto.albero.cartella(self.cartella)
        convergenza = json.loads((cartella / NOME_CONVERGENZA).read_text(encoding="utf-8"))

        # Razionale biologico: traccia la ripartizione delle basi usate nella stima
        # tra campioni biologici, controlli positivi (mock) e controlli negativi,
        # rendendo ispezionabile la composizione del pool di apprendimento.
        classi = {c.accession: c.classe.value for c in contesto.inventario}
        riepilogo: dict[str, Any] = {}
        for m in modelli:
            per_classe: dict[str, int] = defaultdict(int)
            for campione in m.usati:
                per_classe[classi[campione]] += letture[campione] * lunghezza
            riepilogo[m.nome] = {
                "corsa": m.corsa,
                "modello": f"modello_{m.nome}.rds",
                "grafico": f"modello_{m.nome}.png",
                "campioni": len(m.campioni),
                "campioni_usati": list(m.usati),
                "basi_usate": m.basi_usate,
                "basi_disponibili": m.basi_disponibili,
                "tutte_le_basi_usate": m.tutte_usate,
                "basi_usate_per_classe": dict(sorted(per_classe.items())),
                "nbases": config.err.nbases,
                "convergenza": convergenza[m.nome]["convergenza"],
                "iterazioni": convergenza[m.nome]["iterazioni"],
            }

        righe = ["campione\tcorsa\tmodello"] + [
            f"{c}\t{m.corsa or ''}\tmodello_{m.nome}.rds" for m in modelli for c in m.campioni
        ]
        artefatti = esito.artefatti + (
            contesto.albero.scrivi_testo(self.cartella, NOME_CORRISPONDENZA, "\n".join(righe) + "\n"),
            contesto.albero.scrivi_testo(
                self.cartella, NOME_MODELLI,
                json.dumps(riepilogo, indent=2, ensure_ascii=False, sort_keys=True) + "\n",
            ),
        )

        # Razionale sistemistico: se un modello non converge entro err.max_consist
        # iterazioni (E-S3-01), l'azione correttiva raddoppia err.nbases per ampliare
        # il campione di stima. Tuttavia, se il modello usava gia' tutte le basi
        # disponibili nella corsa (m.tutte_usate), raddoppiare err.nbases passerebbe
        # a learnErrors gli stessi identici dati; la fase valorizza allora
        # RITENTARE_INUTILE per sopprimere un secondo tentativo superfluo.
        falliti = [m for m in modelli if not riepilogo[m.nome]["convergenza"]]
        if falliti:
            contesto_errore: dict[str, Any] = {"modelli": [m.nome for m in falliti]}
            if all(m.tutte_usate for m in falliti):
                contesto_errore[RITENTARE_INUTILE] = (
                    "la stima usava gia' tutte le basi disponibili ("
                    + ", ".join(f"{m.nome}: {m.basi_disponibili} basi" for m in falliti)
                    + f", meno di err.nbases {config.err.nbases:g}), e raddoppiare "
                    "err.nbases darebbe un tentativo identico"
                )
            raise errore(
                "E-S3-01",
                f"i modelli {', '.join(m.nome for m in falliti)} non convergono entro "
                f"err.max_consist {config.err.max_consist} iterazioni",
                **contesto_errore,
            )

        metriche = {
            "modelli": len(modelli),
            "basi_usate": {m.nome: m.basi_usate for m in modelli},
            "iterazioni": {m.nome: riepilogo[m.nome]["iterazioni"] for m in modelli},
        }
        return Produzione(artefatti, metriche)
