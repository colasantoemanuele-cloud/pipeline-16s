"""Valori predefiniti dei parametri di configurazione.

I valori stanno qui e non dentro lo schema perché la domanda «quali valori
andrebbero rivisti passando a un altro dataset?» deve potersi rispondere
leggendo un solo file.

I parametri sono di tre specie, e la distinzione è sostanziale:

* **Predefiniti generici**: scelte ragionevoli indipendenti dal dataset.
  Restano validi finché non c'è una ragione specifica per cambiarli.

* **Predefiniti tarati sul dataset di riferimento** (:data:`FATTI_OSD734`):
  scelte di metodo il cui valore è stato fissato su un fatto accertato di
  OSD-734. Hanno un predefinito, perché su un dataset simile sono un punto di
  partenza sensato, ma su un altro dataset vanno riesaminati: il report elenca
  quelli che un'esecuzione ha preso senza dichiararli.

* **Parametri obbligatori** (:data:`OBBLIGATORI`): descrivono il formato dei
  metadati e l'esperimento, cioè il dataset e non la pipeline. Non hanno alcun
  predefinito: vanno dichiarati tutti nella configurazione, anche vuoti o nulli
  quando non sono pertinenti, e G15 non parte finché ne manca uno. I valori che
  hanno per OSD-734 sono in :data:`ESEMPIO_OSD734`, a solo titolo di esempio.
"""

from __future__ import annotations

from collections.abc import Mapping
from types import MappingProxyType
from typing import Any, Final, NamedTuple

# --------------------------------------------------------------------------- #
# Predefiniti generici                                                         #
# --------------------------------------------------------------------------- #

# io: individuazione dei file di lettura
IO_FASTQ_GLOB: Final = "*.fastq.gz"

# meta: lettura della tabella dei metadati
META_DERIVE_MODULE: Final = True

# filter: filtraggio e troncamento delle letture
FILTER_TRUNCLEN_SHORTFALL_WARN: Final = 10
FILTER_MAXEE: Final = 2.0
FILTER_TRUNCQ: Final = 2
FILTER_MAXN: Final = 0
FILTER_RM_PHIX: Final = True

# err: apprendimento del modello di errore
ERR_NBASES: Final = 1.0e8
ERR_MAX_CONSIST: Final = 10
ERR_RANDOMIZE: Final = True

# dada: inferenza delle varianti
DADA_POOL: Final = "pseudo"
DADA_OMEGA_A: Final = 1.0e-40

# chimera: rimozione delle chimere
CHIMERA_METHOD: Final = "consensus"
CHIMERA_MIN_FOLD_PARENT_OVER_ABUNDANCE: Final = 2.0
CHIMERA_MIN_PARENT_ABUNDANCE: Final = 2
CHIMERA_MIN_SAMPLE_FRACTION: Final = 0.9
CHIMERA_ALLOW_ONE_OFF: Final = False

# asv: selezione delle sequenze
ASV_LEN_TOL: Final = 0

# tax: assegnazione tassonomica
TAX_CLASSIFIER: Final = "naive_bayes"
TAX_TRY_RC: Final = True
TAX_REF_BAD_TAXA: Final = None

# filt: filtraggio tassonomico
FILT_REMOVE_NA_PHYLUM: Final = True
FILT_EXCLUDE_TAXA: Final = ("Chloroplast", "Mitochondria", "Eukaryota")

# phylo: albero filogenetico
# Disattivata per difetto: su sequenze corte un albero costruito da zero e'
# debolmente risolto, e rischia di suggerire relazioni che il dato non sostiene.
PHYLO_ENABLED: Final = False
# Varianti finali oltre le quali S9 si ferma prima del calcolo (E-S9-01).
PHYLO_MAX_SEQS: Final = 5000
# L'allineatore multiplo e il modello evolutivo della massima verosimiglianza.
PHYLO_ALIGNER: Final = "decipher"
PHYLO_MODEL: Final = "GTR+G+I"

# decontam: rimozione dei contaminanti
DECONTAM_MIN_BLANKS: Final = 5
# Come si combinano le probabilita' delle piastre nella decontaminazione per
# piastra: e' il batch.combine di decontam::isContaminant, il cui predefinito
# (decontam 1.28.0) e' "minimum". Dichiararlo evita di dipendere da un
# predefinito che una versione nuova potrebbe cambiare.
DECONTAM_BATCH_COMBINE: Final = "minimum"

# prev: filtro di prevalenza
PREV_MIN_COUNT: Final = 2
PREV_APPLY: Final = True

# qc: soglie dei controlli di qualità
QC_MIN_READS_RAW: Final = 1000
# Frazione massima delle letture dei campioni biologici che la decontaminazione
# puo' rimuovere (S12, E-S12-02).
QC_MAX_FRAC_CONTAMINANT: Final = 0.40
# Da dove viene la profondita' minima dei campioni: la soglia derivata dai
# controlli positivi (S11) quando la curva e' attendibile, altrimenti
# qc.min_reads_raw sulle letture grezze; "fixed" usa sempre il ripiego.
QC_MIN_READS_MODE: Final = "katharoseq_if_available"
# Letture minime di un campione nell'oggetto finale, dopo tutti i filtri di S13:
# per campione; sotto, il campione esce dall'oggetto finale con il motivo.
QC_MIN_READS_FINAL: Final = 1000
# Frazione minima delle letture trattenute dall'insieme dei campioni finali, dalle
# letture senza chimere (S7) a quelle dell'oggetto finale (S14, E-S14-01): misura
# l'insieme, non il singolo campione.
QC_MIN_FRAC_READS_RETAINED: Final = 0.40
QC_MAX_ASV_COUNT: Final = 300000
QC_MAX_ZEROED_SAMPLES: Final = 0
QC_MAX_FRAC_LOST_FILTER: Final = 0.30
QC_WARN_FRAC_CHIMERIC: Final = 0.25
QC_STOP_FRAC_CHIMERIC: Final = 0.50
QC_MIN_FRAC_PHYLUM: Final = 0.80
# Frazione di letture che iniziano col primer oltre la quale G10 lo considera
# presente e chiede di tagliarlo.
QC_MAX_PRIMER_HIT_FRAC: Final = 0.05
# Letture ispezionate per file dai gate che leggono le sequenze. La
# validazione deve costare minuti, non ore: nessun gate legge un file intero.
QC_HEAD_READS: Final = 10000

# retry: nuovi tentativi sulle fasi fallite
RETRY_ENABLED: Final = True
RETRY_MAX_ATTEMPTS: Final = 2
# Elenco chiuso dei codici di errore per cui e' ammesso un nuovo tentativo
# automatico. Il criterio non dipende dal dataset: si ritenta solo dove
# l'azione correttiva non modifica alcuna assunzione dell'analisi, e vale
# quindi su qualunque dataset.
RETRY_WHITELIST: Final = ("E-S2-03", "E-S3-01", "E-S4-02", "E-S5-01")

# ctrl: validazione della corsa dai controlli positivi (S11)
# Controlli valutabili minimi: per costruire una curva, e per giudicare la
# conformita' dei controlli di un livello di concentrazione.
CTRL_MIN_POSITIVES: Final = 3
# Frazione minima di controlli positivi conformi alla loro concentrazione.
CTRL_MIN_POSITIVE_PASS_FRAC: Final = 0.75
# Se vero, una frazione di conformi sotto il minimo ferma l'esecuzione
# (E-S11-03); se falso si prosegue con un avviso (E-S11-04).
CTRL_POSITIVE_GATE: Final = False

# katharoseq: la curva fedelta'-profondita' dei controlli positivi (S11)
# Rango a cui si sommano le varianti del taxon atteso.
KATHAROSEQ_COLLAPSE_RANK: Final = "Genus"
KATHAROSEQ_CURVE_MODEL: Final = "allosteric_sigmoid"
# Fedelta' a cui si legge la profondita' minima sulla curva.
KATHAROSEQ_TARGET_SENSITIVITY: Final = 0.90
# Bonta' minima dell'adattamento (R^2 sulla scala della fedelta').
KATHAROSEQ_MIN_R2: Final = 0.80
# Stadio delle letture su cui si misurano profondita' e fedelta', e a cui si
# applica la soglia derivata.
KATHAROSEQ_READ_STAGE: Final = "nonchimeric"

# out: forma degli artefatti prodotti
OUT_SERIALIZATION: Final = "rds"
OUT_TAXA_ARE_ROWS: Final = True
OUT_ASV_ID_SCHEME: Final = "abundance_rank"
OUT_EXPORT_FLAT: Final = True
# Identificativi dei campioni nell'oggetto integrato: l'accession, la chiave
# del join in tutta la pipeline. Il nome del campione resta fra i metadati.
OUT_SAMPLE_ID_SOURCE: Final = "accession"

# run: esecuzione
# File di blocco delle versioni dei pacchetti R. Dichiararlo nella
# configurazione lo fa entrare nel digest: senza, due esecuzioni non sarebbero
# confrontabili rispetto alle versioni R impiegate.
RUN_LOCKFILE: Final = "renv.lock"
RUN_SEED: Final = 100
# run.threads non ha un valore fisso: per difetto vale i processori utilizzabili
# dal processo al momento della validazione (config/schema.py).
RUN_BATCH_SIZE: Final = 24
RUN_KEEP_FILTERED_FASTQ: Final = True
# La regola rigorosa sulla provenienza: disattivata per difetto, perche' lo
# sviluppo e i test girano su alberi di lavoro modificati; attiva nelle
# configurazioni congelate.
RUN_STRICT_PROVENANCE: Final = False


# --------------------------------------------------------------------------- #
# Predefiniti tarati sul dataset di riferimento (OSD-734)                      #
# --------------------------------------------------------------------------- #
# Scelte di metodo fissate su un fatto accertato di OSD-734. Su un altro dataset
# vanno riesaminate. Ogni aggiunta qui va riportata anche in FATTI_OSD734, più
# in basso, con il fatto che la giustifica.

# Basi tolte in testa: nessuna, perche' le letture non contengono il primer e
# iniziano nella regione conservata a valle del 515F (G10).
FILTER_TRIMLEFT: Final = 0

# Bootstrap minimo per assegnare un rango, e niente specie: le letture coprono
# 137 delle circa 253 basi dell'amplicone V4, e su sequenze corte il bootstrap
# e' piu' basso a parita' di correttezza.
TAX_MIN_BOOT: Final = 50
TAX_ASSIGN_SPECIES: Final = False

# Decontaminazione per prevalenza: il metodo per frequenza richiede una
# concentrazione del DNA per campione, che i metadati di OSD-734 non hanno.
DECONTAM_METHOD: Final = "prevalence"
# Soglia 0,5 e non 0,1 (il predefinito di decontam): i campioni sono a bassa
# biomassa (i controlli negativi hanno profondita' paragonabili ai biologici),
# e la contaminazione da reagente puo' essere una quota rilevante del segnale.
DECONTAM_THRESHOLD: Final = 0.5

# Frazione minima di letture con il motivo conservato, valutata sulla mediana
# dei soli campioni attesi portatori di segnale: su OSD-734 sei biologici
# legittimi stanno sotto il 25%, uno al 5,8%, e un controllo file per file li
# farebbe fallire.
QC_MIN_MOTIF_FRAC: Final = 0.25

# Frazione minima di campioni in cui una variante deve comparire (S13): 1% e non
# il 5% abituale, perche' i campioni vengono da oltre cento posizioni distinte
# della stazione, e una variante propria di poche posizioni e' segnale, non
# rumore.
PREV_MIN_FRACTION: Final = 0.01

# Modalita' della decontaminazione (S12). Aggregata: le piastre di OSD-734 hanno
# da 5 a 8 controlli negativi, e il minimo di dieci probabilita' stimate
# ciascuna su cosi' pochi negativi, con la soglia 0,5, e' permissivo per
# costruzione (misurato: per piastra si toglierebbe la meta' delle letture dei
# biologici). La modalita' per piastra resta calcolata come diagnostica.
DECONTAM_MODE: Final = "aggregate"


class RiferimentoDataset(NamedTuple):
    """Il valore di un parametro scelto su OSD-734 e il fatto che lo giustifica."""

    #: Il valore adottato per il dataset di riferimento.
    valore: Any
    #: Il fatto accertato su OSD-734 da cui quel valore discende.
    fatto: str


#: I parametri con un predefinito tarato su OSD-734, ciascuno con quel valore e
#: il fatto accertato che lo giustifica. E' l'unica fonte della motivazione: il
#: report la riporta accanto a ogni parametro che un'esecuzione ha preso per
#: difetto senza dichiararlo, e un test la tiene allineata ai marcatori
#: [OSD-734] di config.example.yaml. Un valore ereditato in silenzio su un altro
#: dataset e' il rischio che questo elenco serve a rendere visibile.
FATTI_OSD734: Final[Mapping[str, RiferimentoDataset]] = MappingProxyType({
    "filter.trimLeft": RiferimentoDataset(
        FILTER_TRIMLEFT,
        "le letture non contengono il primer: iniziano nella regione conservata a "
        "valle del 515F (G10)",
    ),
    "tax.min_boot": RiferimentoDataset(
        TAX_MIN_BOOT,
        "le letture coprono 137 delle circa 253 basi dell'amplicone V4: su "
        "sequenze corte il bootstrap è più basso a parità di correttezza",
    ),
    "tax.assign_species": RiferimentoDataset(
        TAX_ASSIGN_SPECIES,
        "137 basi su circa 253 dell'amplicone V4 non bastano a un'assegnazione "
        "attendibile al livello di specie",
    ),
    "decontam.method": RiferimentoDataset(
        DECONTAM_METHOD,
        "i metadati non riportano una concentrazione del DNA per campione, che il "
        "metodo per frequenza richiede",
    ),
    "decontam.threshold": RiferimentoDataset(
        DECONTAM_THRESHOLD,
        "campioni a bassa biomassa: i controlli negativi hanno profondità "
        "paragonabili ai biologici, e la contaminazione da reagente può essere una "
        "quota rilevante del segnale",
    ),
    "decontam.mode": RiferimentoDataset(
        DECONTAM_MODE,
        "le piastre hanno da 5 a 8 controlli negativi dichiarati: per piastra, con "
        "la soglia 0,5, si toglierebbe circa la metà delle letture dei biologici "
        "(misurato)",
    ),
    "qc.min_motif_frac": RiferimentoDataset(
        QC_MIN_MOTIF_FRAC,
        "sei campioni biologici legittimi hanno il motivo conservato sotto il 25% "
        "delle letture, uno al 5,8%: la soglia si valuta sulla mediana",
    ),
    "prev.min_fraction": RiferimentoDataset(
        PREV_MIN_FRACTION,
        "i campioni vengono da oltre cento posizioni distinte della stazione: una "
        "variante propria di poche posizioni è segnale, non rumore",
    ),
})

#: Chiavi il cui valore predefinito e' tarato sul dataset di riferimento.
#: Passando a un altro dataset vanno rivalutate una per una.
DERIVATI_DAL_DATASET: Final[tuple[str, ...]] = tuple(FATTI_OSD734)


# --------------------------------------------------------------------------- #
# Parametri obbligatori: nessun predefinito                                    #
# --------------------------------------------------------------------------- #

#: I parametri che descrivono il dataset e non hanno predefinito. Vanno
#: dichiarati tutti nel file di configurazione; quelli non pertinenti al
#: dataset si dichiarano nulli o vuoti (per esempio le etichette dei controlli
#: positivi in un dataset che non ne ha). G15 respinge la configurazione
#: elencando quelli che mancano (E-G15-10).
OBBLIGATORI: Final[tuple[str, ...]] = (
    # come si riconoscono file e campioni
    "io.accession_regex",
    "meta.sample_id_column",
    "meta.accession_column",
    # il formato dei metadati
    "meta.module_regex",
    "meta.module_column",
    "meta.non_surface_positions",
    "meta.batch_key_column",
    "meta.batch_module_column",
    "err.batch_column",
    "decontam.batch_column",
    "ctrl.column",
    "ctrl.blank_values",
    "ctrl.positive_values",
    "ctrl.biological_values",
    "ctrl.blank_override_column",
    "ctrl.blank_override_values",
    "katharoseq.cell_count_column",
    "out.study_columns",
    "out.batch_columns",
    # l'esperimento
    "filter.truncLen",
    "qc.primer_sequence",
    "qc.conserved_motif",
    "katharoseq.target_taxon",
    "tax.ref_name",
    "tax.ref_version",
)

#: I valori che i parametri obbligatori hanno per OSD-734, con la ragione di
#: ciascuno. Non sono predefiniti e la pipeline non li usa: documentano il
#: dataset di riferimento, sono l'esempio riportato in config.example.yaml e i
#: valori con cui i test costruiscono i loro scenari.
ESEMPIO_OSD734: Final[Mapping[str, RiferimentoDataset]] = MappingProxyType({
    "io.accession_regex": RiferimentoDataset(
        r"[ESD]RX[0-9]{4,}",
        "i file e la tabella di assay riportano l'accession dell'esperimento "
        "(ERX seguito da cifre)",
    ),
    "meta.sample_id_column": RiferimentoDataset(
        "Sample Name",
        "nelle tabelle ISA il nome del campione sta nella colonna Sample Name",
    ),
    "meta.accession_column": RiferimentoDataset(
        "Raw Data File",
        "la tabella di assay riporta il file delle letture, con l'accession nel "
        "nome, nella colonna Raw Data File",
    ),
    "meta.module_regex": RiferimentoDataset(
        r"^([A-Z]{3}[0-9])",
        "i nomi delle posizioni iniziano con il modulo della stazione, tre lettere "
        "e una cifra (LAB1, NOD2, JLP1)",
    ),
    "meta.module_column": RiferimentoDataset(
        "Factor Value[Sample Location]",
        "la posizione di campionamento è dichiarata nella colonna Factor "
        "Value[Sample Location] della tabella campioni di studio",
    ),
    "meta.non_surface_positions": RiferimentoDataset(
        ("Air Sample", "Unopened 3DMM Swab Tube"),
        "i campioni d'aria e i tamponi mai aperti riportano una posizione che non "
        "è una superficie, e non appartengono a un modulo",
    ),
    "meta.batch_key_column": RiferimentoDataset(
        "experiment_accession",
        "il file del lotto identifica i campioni per accession, nella colonna "
        "experiment_accession: i due campioni risequenziati vi si distinguono solo "
        "così",
    ),
    "meta.batch_module_column": RiferimentoDataset(
        "module",
        "il file del lotto dichiara il modulo nella colonna module",
    ),
    "err.batch_column": RiferimentoDataset(
        "run_prefix",
        "nel file del lotto la corsa di sequenziamento sta nella colonna "
        "run_prefix (due corse)",
    ),
    "decontam.batch_column": RiferimentoDataset(
        "extraction_plate_num",
        "nel file del lotto la piastra di estrazione sta nella colonna "
        "extraction_plate_num (dieci piastre da 96)",
    ),
    "ctrl.column": RiferimentoDataset(
        "Characteristics[Material Type]",
        "la classe del campione è dichiarata nella colonna Characteristics[Material "
        "Type] della tabella campioni di studio",
    ),
    "ctrl.blank_values": RiferimentoDataset(
        ("blank control",),
        "i controlli negativi sono dichiarati con il materiale blank control",
    ),
    "ctrl.positive_values": RiferimentoDataset(
        ("Positive Control",),
        "i controlli positivi sono dichiarati con il materiale Positive Control",
    ),
    "ctrl.biological_values": RiferimentoDataset(
        ("Surface swab",),
        "i campioni biologici sono dichiarati con il materiale Surface swab",
    ),
    "ctrl.blank_override_column": RiferimentoDataset(
        "Factor Value[Sample Location]",
        "i tamponi mai aperti si riconoscono dalla posizione dichiarata, nella "
        "colonna Factor Value[Sample Location]",
    ),
    "ctrl.blank_override_values": RiferimentoDataset(
        ("Unopened 3DMM Swab Tube",),
        "33 tamponi dichiarati Surface swab hanno posizione Unopened 3DMM Swab Tube: "
        "non hanno campionato alcuna superficie, e per funzione sono controlli "
        "negativi",
    ),
    "katharoseq.cell_count_column": RiferimentoDataset(
        "katharoseq_cell_count",
        "il file del lotto riporta le cellule seminate in ciascun controllo "
        "positivo nella colonna katharoseq_cell_count",
    ),
    "out.study_columns": RiferimentoDataset(
        ("Factor Value[Spaceflight]", "Parameter Value[Sample Preservation Method]"),
        "la tabella campioni di studio riporta se il campione è stato in volo e "
        "come è stato conservato",
    ),
    "out.batch_columns": RiferimentoDataset(
        ("extraction_plate_id", "well_id", "extraction_date", "primer_plate",
         "run_date", "katharoseq_cell_count"),
        "il file del lotto riporta piastra, pozzetto, date di estrazione e di "
        "sequenziamento, piastra dei primer e cellule seminate nei controlli "
        "positivi",
    ),
    "filter.truncLen": RiferimentoDataset(
        137,
        "le letture del dataset completo sono lunghe 137-151 basi (S1): 137 è la "
        "più corta, e un troncamento maggiore la scarterebbe",
    ),
    "qc.primer_sequence": RiferimentoDataset(
        "GTGYCAGCMGCCGCGGTAA",
        "la regione V4 è stata amplificata con il primer 515F",
    ),
    "qc.conserved_motif": RiferimentoDataset(
        r"TAC[AG].AGG..GC.AGCGTT",
        "le letture iniziano con il motivo conservato a valle del 515F",
    ),
    "katharoseq.target_taxon": RiferimentoDataset(
        "Variovorax",
        "i controlli positivi contengono un solo ceppo, Variovorax sp. OAS795",
    ),
    "tax.ref_name": RiferimentoDataset(
        "SILVA",
        "il riferimento tassonomico adottato per lo studio è SILVA",
    ),
    "tax.ref_version": RiferimentoDataset(
        "138",
        "la versione del riferimento adottata per lo studio è la 138",
    ),
})
assert tuple(ESEMPIO_OSD734) == OBBLIGATORI
