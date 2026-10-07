"""Fase S10: l'assemblaggio dell'oggetto integrato.

Riunisce in un oggetto phyloseq (``R/10_phyloseq.R``) quattro componenti
allineati: la tabella dei conteggi di S7, la tabella tassonomica di S8, i
metadati dei campioni e le sequenze di riferimento delle varianti. L'oggetto
integrato non porta mai l'albero filogenetico: S9 lo costruisce sulle varianti
finali, dopo i filtri di S13, ed e' S14 ad aggiungerlo all'oggetto finale.
Scrive in ``10_phyloseq/``.

**Tutti i campioni dell'inventario.** I campioni rimasti senza letture filtrate
non hanno una riga nelle tabelle di S5-S7; l'oggetto integrato li contiene
comunque, con conteggi a zero, e ``riepilogo.json`` li elenca. Un campione non
sparisce in silenzio fra i metadati e i conteggi: il controllo e' che le
colonne della tabella siano esattamente gli accession dell'inventario,
nell'ordine dell'inventario.

**L'orientamento, verificato sugli identificativi.** Con ``out.taxa_are_rows``
vero le varianti stanno sulle righe e i campioni sulle colonne. Lo script lo
verifica a ogni esecuzione confrontando i nomi delle righe con gli
identificativi delle varianti e quelli delle colonne con gli accession, non le
dimensioni: una matrice con tanti campioni quante varianti, trasposta,
supererebbe un controllo sulle dimensioni e darebbe risultati plausibili e
sbagliati. Se non corrisponde, ``E-S10-01``.

**Gli identificativi delle varianti** (``out.asv_id_scheme`` =
``abundance_rank``): ``ASV1``, ``ASV2``, ... in ordine di letture decrescenti
nei campioni **biologici**; a parita', letture in tutti i campioni; a parita'
ancora, la sequenza in ordine lessicografico, che rende la numerazione
deterministica. Il denominatore sono i biologici perche' a questo punto
l'oggetto contiene ancora i controlli: i positivi portano un solo organismo a
concentrazioni alte, e contarli metterebbe in testa alla numerazione il ceppo
di calibrazione invece delle varianti dell'ambiente studiato. Gli
identificativi si assegnano qui e non si rinumerano nelle fasi successive:
dopo i filtri la numerazione avra' dei vuoti, ed e' cio' che permette di
seguire una variante lungo la pipeline. L'identita' stabile di una variante
resta la sua sequenza, conservata fra le sequenze di riferimento e in
``varianti.tsv``.

**Gli identificativi dei campioni** sono gli accession
(``out.sample_id_source``); il nome del campione resta fra i metadati.

**I metadati senza rinomine silenziose.** ``data.frame`` in R rende sintattici
i nomi di colonna senza avvisare: ``Characteristics[Material Type]``
diventerebbe ``Characteristics.Material.Type.``. I nomi nell'oggetto si
scelgono quindi qui, esplicitamente: le colonne dell'inventario hanno i nomi
del crosswalk di S0, quelle chieste con ``out.study_columns`` e
``out.batch_columns`` un nome sintattico derivato dall'originale
(:func:`nome_nell_oggetto`). La corrispondenza con le colonne originali sta in
``colonne_metadati.tsv``, e lo script verifica che l'oggetto porti esattamente
quei nomi. Tutti i valori sono testo, come nelle tabelle d'origine.

**Le informazioni accessorie.** Il bootstrap di ciascun rango (S8) e la
marcatura dei taxa col difetto noto del riferimento non hanno uno slot
nell'oggetto: stanno in ``varianti_accessorie.tsv``, indicizzato per
identificativo di variante.
"""

from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Any, ClassVar, Final

from amplicon16s.errors.exceptions import errore
from amplicon16s.io_layer.artifacts import Artefatto, Fase
from amplicon16s.metadata.crosswalk import chiave_dalla_tabella
from amplicon16s.metadata.tabelle import (
    COLONNE_INVENTARIO,
    intestazione,
    leggi_tsv,
    nome_nell_oggetto,
    tabella_di_studio,
)
from amplicon16s.metadata.models import Inventario
from amplicon16s.rbridge.runner import cartella_r, esegui_script
from amplicon16s.runner.graph import Passo
from amplicon16s.steps.base import PipelineStep, Produzione, StepContext
from amplicon16s.steps.s07_asv_length import NOME_TABELLA_ASV
from amplicon16s.steps.s08_taxonomy import NOME_DIFETTI

__all__ = [
    "AssemblaggioOggetto",
    "COLONNE_INVENTARIO",
    "NOME_OGGETTO",
    "colonne_metadati",
    "nome_nell_oggetto",
]

NOME_SCRIPT: Final = "10_phyloseq.R"
#: L'oggetto integrato che le fasi successive leggono.
NOME_OGGETTO: Final = "ps_integrato.rds"
NOME_METADATI: Final = "metadati_campioni.tsv"
NOME_COLONNE: Final = "colonne_metadati.tsv"
NOME_RIEPILOGO: Final = "riepilogo.json"

#: Le origini registrate in colonne_metadati.tsv.
ORIGINE_LOTTO: Final = "file di arricchimento (io.batch_table)"
ORIGINE_ASSAY: Final = "tabella di assay (io.assay_table)"

def colonne_metadati(config: Any) -> list[dict[str, str]]:
    """Le colonne dei metadati dell'oggetto, con l'origine di ciascuna.

    Ogni voce ha il nome nell'oggetto, la tabella d'origine, la colonna
    originale e come se ne ricava il valore. Una colonna richiesta che non si
    puo' portare nell'oggetto senza ambiguita' solleva ``E-S10-02``.
    """
    io, meta = config.io, config.meta
    lotto = ORIGINE_LOTTO
    # Senza io.study_table classe e variabili vengono dalla tabella di assay.
    studio = (
        "tabella campioni di studio (io.study_table)" if io.study_table is not None
        else ORIGINE_ASSAY
    )
    colonne = [
        {"colonna": "accession", "origine": "tabella di assay (io.assay_table)",
         "colonna_originale": meta.accession_column,
         "valore": "accession estratto con io.accession_regex; identificativo del campione"},
        {"colonna": "sample_name", "origine": "tabella di assay (io.assay_table)",
         "colonna_originale": meta.sample_id_column, "valore": "valore originale"},
        {"colonna": "classe", "origine": studio,
         "colonna_originale": config.ctrl.column,
         "valore": "classe del campione, da ctrl.blank_values, positive_values e "
                   "biological_values"},
        {"colonna": "materiale", "origine": studio,
         "colonna_originale": config.ctrl.column, "valore": "valore originale"},
        {"colonna": "posizione", "origine": studio,
         "colonna_originale": meta.module_column or "",
         "valore": "valore originale; vuoto se non applicabile" if meta.module_column else
                   "vuoto: meta.module_column nullo"},
        {"colonna": "modulo",
         "origine": f"{lotto}, oppure tabella campioni di studio"
                    if io.batch_table is not None or io.study_table is not None else studio,
         "colonna_originale": ", oppure ".join(
             c for c in (meta.batch_module_column, meta.module_column) if c),
         "valore": "dal file di arricchimento se c'e', altrimenti derivato dalla posizione "
                   "con meta.module_regex; vuoto per le posizioni che non sono superfici"},
        {"colonna": "piastra", "origine": lotto,
         "colonna_originale": config.decontam.batch_column or "",
         "valore": "valore originale" if config.decontam.batch_column else
                   "vuoto: decontam.batch_column nullo"},
        {"colonna": "corsa", "origine": lotto,
         "colonna_originale": config.err.batch_column or "",
         "valore": "valore originale" if config.err.batch_column else
                   "vuoto: err.batch_column nullo"},
    ]
    if io.batch_table is None:
        for voce in colonne:
            if voce["colonna"] in ("piastra", "corsa"):
                voce["valore"] = "vuoto: io.batch_table non indicato"

    usati = set(COLONNE_INVENTARIO)
    richieste = [("studio", c) for c in config.out.study_columns]
    richieste += [("lotto", c) for c in config.out.batch_columns]
    intestazioni = {"studio": intestazione(tabella_di_studio(config)[0])}
    if io.batch_table is not None:
        intestazioni["lotto"] = intestazione(Path(io.batch_table))
    origini = {"studio": studio, "lotto": lotto}
    for tabella, originale in richieste:
        if tabella not in intestazioni:
            raise errore(
                "E-S10-02",
                f"{originale!r} e' in out.batch_columns, ma io.batch_table non e' indicato",
                colonna=originale,
            )
        occorrenze = intestazioni[tabella].count(originale)
        if occorrenze != 1:
            motivo = "assente" if occorrenze == 0 else f"ripetuta {occorrenze} volte"
            raise errore(
                "E-S10-02",
                f"{originale!r}: {motivo} nell'intestazione, {origini[tabella]}",
                colonna=originale,
            )
        nome = nome_nell_oggetto(originale)
        if nome in usati:
            raise errore(
                "E-S10-02",
                f"{originale!r} diventerebbe {nome!r}, gia' usato da un'altra colonna",
                colonna=originale,
            )
        usati.add(nome)
        colonne.append({
            "colonna": nome, "origine": origini[tabella],
            "colonna_originale": originale, "valore": "valore originale",
        })
    return colonne


def _valori_metadati(
    config: Any, inventario: Inventario, colonne: list[dict[str, str]]
) -> list[list[str]]:
    """Una riga per campione dell'inventario, nell'ordine di ``colonne``.

    Le colonne aggiuntive si agganciano come in S0: la tabella di studio per
    nome del campione (S0 ha verificato che il nome vi compaia una sola
    volta), il file di arricchimento per accession.
    """
    io, meta = config.io, config.meta
    studio: dict[str, dict[str, str]] = {}
    tabella, colonna_id = tabella_di_studio(config)
    for riga in leggi_tsv(tabella):
        studio.setdefault(riga.get(colonna_id, ""), riga)
    lotto: dict[str, dict[str, str]] = {}
    if io.batch_table is not None:
        espressione = re.compile(io.accession_regex)
        for riga in leggi_tsv(Path(io.batch_table)):
            try:
                accession = chiave_dalla_tabella(riga.get(meta.batch_key_column, ""), espressione)
            except ValueError:
                continue
            lotto.setdefault(accession, riga)

    righe = []
    for c in inventario:
        valori = {
            "accession": c.accession, "sample_name": c.nome, "classe": c.classe.value,
            "materiale": c.materiale, "posizione": c.posizione or "",
            "modulo": c.modulo or "", "piastra": c.piastra or "", "corsa": c.corsa or "",
        }
        riga = []
        for voce in colonne:
            nome = voce["colonna"]
            if nome in valori:
                riga.append(valori[nome])
            elif voce["origine"] != ORIGINE_LOTTO:
                riga.append(studio.get(c.nome, {}).get(voce["colonna_originale"], ""))
            else:
                riga.append(lotto.get(c.accession, {}).get(voce["colonna_originale"], ""))
        righe.append(riga)
    return righe


def _tsv(intestazione: list[str], righe: list[list[str]]) -> str:
    """Una tabella separata da tabulazioni, senza virgolette.

    Un valore con una tabulazione o un a capo spezzerebbe la tabella: non e'
    ammesso, e lo si dice invece di alterarlo.
    """
    for riga in [intestazione, *righe]:
        for valore in riga:
            if "\t" in valore or "\n" in valore or "\r" in valore:
                raise ValueError(f"valore con tabulazione o a capo: {valore!r}")
    return "".join("\t".join(riga) + "\n" for riga in [intestazione, *righe])


class AssemblaggioOggetto(PipelineStep):
    """S10: l'oggetto integrato, con tutti i campioni dell'inventario."""

    passo: ClassVar[Passo] = Passo.S10
    #: 2: classe e variabili dalla tabella di assay quando io.study_table non
    #: e' indicata; chiave del file di arricchimento estratta con io.accession_regex.
    versione: ClassVar[int] = 2
    script_r: ClassVar[str | None] = NOME_SCRIPT
    #: La forma dell'oggetto (orientamento, identificativi, colonne portate;
    #: serializzazione ed export sono di S14); le tabelle e le colonne da cui
    #: vengono i metadati, e i nomi delle colonne d'origine registrati nella
    #: corrispondenza. Il gruppo phylo non c'e': l'oggetto integrato non ha
    #: l'albero, e attivare la filogenesi non lo cambia.
    parametri: ClassVar[tuple[str, ...]] = (
        "out.taxa_are_rows", "out.asv_id_scheme", "out.sample_id_source",
        "out.study_columns", "out.batch_columns", "meta", "ctrl.column",
        "io.assay_table", "io.study_table", "io.batch_table", "io.accession_regex",
        "decontam.batch_column", "err.batch_column",
    )

    def calcola(self, contesto: StepContext) -> Produzione:
        """Scrive i metadati e la corrispondenza delle colonne, poi esegue
        ``R/10_phyloseq.R``, che assembla e verifica l'oggetto.
        """
        if contesto.inventario is None:
            raise RuntimeError("S10 richiede l'inventario prodotto da S0")
        config = contesto.config
        albero = contesto.albero
        inventario = contesto.inventario

        colonne = colonne_metadati(config)
        nomi = [v["colonna"] for v in colonne]
        artefatti: list[Artefatto] = [
            albero.scrivi_testo(self.cartella, NOME_COLONNE, _tsv(
                ["colonna", "origine", "colonna_originale", "valore"],
                [[v["colonna"], v["origine"], v["colonna_originale"], v["valore"]]
                 for v in colonne],
            )),
            albero.scrivi_testo(self.cartella, NOME_METADATI, _tsv(
                nomi, _valori_metadati(config, inventario, colonne)
            )),
        ]

        # La marcatura dei taxa difettosi c'e' solo se S8 l'ha prodotta: lo
        # dice il suo manifesto, non la presenza di un file nella cartella.
        manifesto_s8 = albero.manifesto_passo(Passo.S8, Fase.TAXONOMY)
        difetti = (
            str(albero.cartella(Fase.TAXONOMY) / NOME_DIFETTI)
            if manifesto_s8 is not None and NOME_DIFETTI in manifesto_s8.nomi
            else None
        )
        esito = esegui_script(
            cartella_r() / NOME_SCRIPT,
            {
                "tabella": str(albero.cartella(Fase.CHIMERA) / NOME_TABELLA_ASV),
                "tassonomia": str(albero.cartella(Fase.TAXONOMY) / "tassonomia.rds"),
                "difetti": difetti,
                "metadati": str(albero.cartella(self.cartella) / NOME_METADATI),
                "colonne": nomi,
                "campioni": list(inventario.accessioni),
                "biologici": [c.accession for c in inventario.biologici],
                "taxa_are_rows": config.out.taxa_are_rows,
                "prefisso": "ASV",
            },
            albero,
            self.cartella,
            passo=self.passo,
            logger=contesto.logger,
        )
        artefatti += esito.artefatti
        riepilogo = json.loads(
            (albero.cartella(self.cartella) / NOME_RIEPILOGO).read_text(encoding="utf-8")
        )
        metriche = {
            k: riepilogo[k] for k in (
                "campioni", "varianti", "letture", "campioni_senza_letture", "prima_variante",
            )
        }
        return Produzione(tuple(artefatti), metriche)
