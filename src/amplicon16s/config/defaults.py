"""Valori predefiniti dei parametri di configurazione.

I valori stanno qui e non dentro lo schema perché la domanda «quali valori
andrebbero rivisti passando a un altro dataset?» deve potersi rispondere
leggendo un solo file.

Sono divisi in due categorie, e la distinzione è sostanziale:

* **Predefiniti generici** — scelte ragionevoli indipendenti dal dataset.
  Restano validi finché non c'è una ragione specifica per cambiarli.

* **Predefiniti derivati dal dataset di riferimento** — valori ricavati da
  caratteristiche accertate di OSD-734: nomi di colonne dei metadati, etichette
  usate per distinguere i controlli, taxon atteso nei controlli positivi. Su
  dataset sono quasi certamente sbagliati e vanno rivisti uno per uno. Il
  codice resta generale: è la configurazione a cambiare.

L'elenco delle chiavi della seconda categoria è esposto in
:data:`DERIVATI_DAL_DATASET`, così può essere stampato in un report o
controllato da un test invece di restare una convenzione a voce.
"""

from __future__ import annotations

from typing import Final

# --------------------------------------------------------------------------- #
# Predefiniti generici                                                         #
# --------------------------------------------------------------------------- #

# io — individuazione dei file di lettura e degli accession
IO_FASTQ_GLOB: Final = "*.fastq.gz"
IO_ACCESSION_REGEX: Final = r"(E|S|D)RX[0-9]{4,}"

# meta — lettura della tabella dei metadati
META_SAMPLE_ID_COLUMN: Final = "Sample Name"
META_ACCESSION_COLUMN: Final = "Raw Data File"
META_DERIVE_MODULE: Final = True
META_MODULE_REGEX: Final = r"^([A-Z]{3}[0-9])"

# filter — filtraggio e troncamento delle letture
FILTER_TRUNCLEN: Final = 137
FILTER_TRUNCLEN_SHORTFALL_WARN: Final = 10
FILTER_TRIMLEFT: Final = 0
FILTER_MAXEE: Final = 2.0
FILTER_TRUNCQ: Final = 2
FILTER_MAXN: Final = 0
FILTER_RM_PHIX: Final = True

# err — apprendimento del modello di errore
ERR_NBASES: Final = 1.0e8
ERR_MAX_CONSIST: Final = 10
ERR_RANDOMIZE: Final = True

# dada — inferenza delle varianti
DADA_POOL: Final = "pseudo"
DADA_OMEGA_A: Final = 1.0e-40

# chimera — rimozione delle chimere
CHIMERA_METHOD: Final = "consensus"
CHIMERA_MIN_FOLD_PARENT_OVER_ABUNDANCE: Final = 2.0
CHIMERA_MIN_PARENT_ABUNDANCE: Final = 2
CHIMERA_MIN_SAMPLE_FRACTION: Final = 0.9
CHIMERA_ALLOW_ONE_OFF: Final = False

# asv — selezione delle sequenze
ASV_LEN_TOL: Final = 0

# tax — assegnazione tassonomica
TAX_CLASSIFIER: Final = "naive_bayes"
TAX_MIN_BOOT: Final = 50
TAX_TRY_RC: Final = True
TAX_ASSIGN_SPECIES: Final = False

# filt — filtraggio tassonomico
FILT_REMOVE_NA_PHYLUM: Final = True
FILT_EXCLUDE_TAXA: Final = ("Chloroplast", "Mitochondria", "Eukaryota")

# phylo — albero filogenetico
PHYLO_ENABLED: Final = False
PHYLO_MAX_SEQS: Final = 5000

# decontam — rimozione dei contaminanti
DECONTAM_METHOD: Final = "prevalence"
DECONTAM_THRESHOLD: Final = 0.5
DECONTAM_MIN_BLANKS: Final = 5

# prev — filtro di prevalenza
PREV_MIN_FRACTION: Final = 0.01
PREV_MIN_COUNT: Final = 2
PREV_APPLY: Final = True

# qc — soglie dei controlli di qualità
QC_MIN_READS_RAW: Final = 1000
QC_MIN_READS_FILTERED: Final = 500
QC_MIN_READS_FINAL: Final = 1000
QC_MAX_ASV_COUNT: Final = 300000
QC_WARN_FRAC_CHIMERIC: Final = 0.25
QC_STOP_FRAC_CHIMERIC: Final = 0.50
# Frazione di letture che iniziano col primer oltre la quale G10 lo considera
# presente e chiede di tagliarlo.
QC_MAX_PRIMER_HIT_FRAC: Final = 0.05
# Letture ispezionate per file dai gate che leggono le sequenze. La
# validazione deve costare minuti, non ore: nessun gate legge un file intero.
QC_HEAD_READS: Final = 10000
# Frazione minima di letture che devono presentare il motivo conservato,
# valutata sulla mediana dei soli campioni attesi portatori di segnale.
QC_MIN_MOTIF_FRAC: Final = 0.25

# retry — nuovi tentativi sulle fasi fallite
RETRY_ENABLED: Final = True
RETRY_MAX_ATTEMPTS: Final = 2
# Elenco chiuso dei codici di errore per cui e' ammesso un nuovo tentativo
# automatico. Il criterio non dipende dal dataset: si ritenta solo dove
# l'azione correttiva non modifica alcuna assunzione dell'analisi, e vale
# quindi su qualunque dataset.
RETRY_WHITELIST: Final = ("E-S2-03", "E-S3-01", "E-S4-02", "E-S5-01")

# out — forma degli artefatti prodotti
OUT_SERIALIZATION: Final = "rds"
OUT_TAXA_ARE_ROWS: Final = True
OUT_ASV_ID_SCHEME: Final = "abundance_rank"
OUT_EXPORT_FLAT: Final = True

# run — esecuzione
# File di blocco delle versioni dei pacchetti R. Dichiararlo nella
# configurazione lo fa entrare nel digest: senza, due esecuzioni non sarebbero
# confrontabili rispetto alle versioni R impiegate.
RUN_LOCKFILE: Final = "renv.lock"
RUN_SEED: Final = 100
RUN_THREADS: Final = 16
RUN_BATCH_SIZE: Final = 24


# --------------------------------------------------------------------------- #
# Predefiniti derivati dal dataset di riferimento (OSD-734)                     #
# --------------------------------------------------------------------------- #
# Attenzione: questi valori descrivono OSD-734, non la pipeline. Su un altro
# dataset vanno rivisti tutti. Ogni aggiunta qui va riportata anche in
# DERIVATI_DAL_DATASET, più in basso.

# Valori della colonna della posizione che denotano una posizione che non e'
# una superficie. I campioni che li riportano non appartengono a un modulo, e
# non vi vanno forzati: raggrupparli per modulo mescolerebbe l'aria di un
# locale con le sue superfici. I valori vuoti e i marcatori convenzionali di
# "non applicabile" sono gia' trattati come posizione assente, e non serve
# elencarli qui.
META_NON_SURFACE_POSITIONS: Final = ("Air Sample", "Unopened 3DMM Swab Tube")

# Colonna della tabella campioni di studio da cui si deriva il modulo. Il
# modulo si ricava dalla posizione dichiarata, non dal nome del campione: un
# campione d'aria puo' chiamarsi come una superficie, e derivarlo dal nome lo
# attribuirebbe a un modulo che non gli compete.
META_MODULE_COLUMN: Final = "Factor Value[Sample Location]"

# Colonna con cui il file facoltativo di arricchimento identifica il campione.
# Deve contenere l'accession: il nome del campione si ripete fra repliche.
META_BATCH_KEY_COLUMN: Final = "experiment_accession"

# Colonna del file facoltativo che dichiara il modulo. Quando c'e', il modulo
# viene da li': e' il dato originale, con il suo nome. La derivazione per
# espressione regolare resta il ripiego per i dataset privi di quel file.
META_BATCH_MODULE_COLUMN: Final = "module"

# Colonna dei metadati che identifica il lotto di sequenziamento.
ERR_BATCH_COLUMN: Final = "run_prefix"

# Colonna dei metadati che identifica la piastra di estrazione.
DECONTAM_BATCH_COLUMN: Final = "extraction_plate_num"

# Colonna e etichette che distinguono controlli e campioni biologici.
CTRL_COLUMN: Final = "Characteristics[Material Type]"
CTRL_BLANK_VALUES: Final = ("blank control",)
CTRL_POSITIVE_VALUES: Final = ("Positive Control",)
CTRL_BIOLOGICAL_VALUES: Final = ("Surface swab",)

# Riferimento tassonomico impiegato nello studio. Questi due parametri sono
# obbligatori nello schema e non hanno valore predefinito: i valori qui sotto
# documentano lo studio di riferimento e sono quelli usati in
# config.example.yaml, ma vanno sempre dichiarati esplicitamente.
TAX_REF_NAME: Final = "SILVA"
TAX_REF_VERSION: Final = "138"

# Primer di amplificazione atteso a inizio lettura, in codici IUPAC. Dipende
# dalla coppia di primer dello studio: con una regione amplificata diversa
# sarebbe un'altra sequenza. Qui 515F, con Y = C o T e M = A o C.
QC_PRIMER_SEQUENCE: Final = "GTGYCAGCMGCCGCGGTAA"

# Motivo conservato atteso subito a valle del primer, come espressione
# regolare. Serve a G10 da controllo positivo: cercare la sola assenza del
# primer non distinguerebbe un file corretto da un file privo di segnale.
# Dipende anch'esso dalla regione amplificata.
QC_CONSERVED_MOTIF: Final = r"TAC[AG].AGG..GC.AGCGTT"

# Taxon atteso nei controlli positivi, usato per la calibrazione KatharoSeq.
KATHAROSEQ_TARGET_TAXON: Final = "Variovorax"


#: Chiavi il cui valore predefinito descrive il dataset di riferimento e non la
#: pipeline. Passando a un altro dataset vanno rivalutate una per una.
DERIVATI_DAL_DATASET: Final[tuple[str, ...]] = (
    "err.batch_column",
    "decontam.batch_column",
    "meta.module_column",
    "meta.non_surface_positions",
    "meta.batch_key_column",
    "meta.batch_module_column",
    "ctrl.column",
    "ctrl.blank_values",
    "ctrl.positive_values",
    "ctrl.biological_values",
    "tax.ref_name",
    "tax.ref_version",
    "katharoseq.target_taxon",
    "qc.primer_sequence",
    "qc.conserved_motif",
)
