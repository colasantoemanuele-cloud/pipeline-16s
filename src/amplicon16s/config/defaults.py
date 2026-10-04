"""Valori predefiniti dei parametri di configurazione.

I valori stanno qui e non dentro lo schema perché la domanda «quali valori
andrebbero rivisti passando a un altro dataset?» deve potersi rispondere
leggendo un solo file.

Sono divisi in due categorie, e la distinzione è sostanziale:

* **Predefiniti generici**: scelte ragionevoli indipendenti dal dataset.
  Restano validi finché non c'è una ragione specifica per cambiarli.

* **Predefiniti derivati dal dataset di riferimento**: valori ricavati da
  caratteristiche accertate di OSD-734: nomi di colonne dei metadati, etichette
  usate per distinguere i controlli, taxon atteso nei controlli positivi. Su
  dataset sono quasi certamente sbagliati e vanno rivisti uno per uno. Il
  codice resta generale: è la configurazione a cambiare.

L'elenco delle chiavi della seconda categoria è esposto in
:data:`DERIVATI_DAL_DATASET`, e :data:`FATTI_OSD734` porta per ciascuna il
valore di OSD-734 e il fatto accertato che lo giustifica: il report li riporta
e un test li controlla, invece di restare una convenzione a voce.
"""

from __future__ import annotations

from collections.abc import Mapping
from types import MappingProxyType
from typing import Any, Final, NamedTuple

# --------------------------------------------------------------------------- #
# Predefiniti generici                                                         #
# --------------------------------------------------------------------------- #

# io: individuazione dei file di lettura e degli accession
IO_FASTQ_GLOB: Final = "*.fastq.gz"
IO_ACCESSION_REGEX: Final = r"(E|S|D)RX[0-9]{4,}"

# meta: lettura della tabella dei metadati
META_SAMPLE_ID_COLUMN: Final = "Sample Name"
META_ACCESSION_COLUMN: Final = "Raw Data File"
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
PHYLO_ENABLED: Final = False
PHYLO_MAX_SEQS: Final = 5000

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
RUN_THREADS: Final = 16
RUN_BATCH_SIZE: Final = 24
RUN_KEEP_FILTERED_FASTQ: Final = True


# --------------------------------------------------------------------------- #
# Predefiniti derivati dal dataset di riferimento (OSD-734)                     #
# --------------------------------------------------------------------------- #
# Attenzione: questi valori descrivono OSD-734, non la pipeline. Su un altro
# dataset vanno rivisti tutti. Ogni aggiunta qui va riportata anche in
# FATTI_OSD734, più in basso, con il fatto che la giustifica.

# Lunghezza di troncamento: la lettura piu' corta del dataset completo, 137
# basi su letture di 137-151 (S1). Troncare oltre scarterebbe quelle letture.
FILTER_TRUNCLEN: Final = 137

# Basi tolte in testa: nessuna, perche' le letture non contengono il primer e
# iniziano nella regione conservata a valle del 515F (G10).
FILTER_TRIMLEFT: Final = 0

# Bootstrap minimo per assegnare un rango, e niente specie: le letture coprono
# 137 delle circa 253 basi dell'amplicone V4, e su sequenze corte il bootstrap
# e' piu' basso a parita' di correttezza.
TAX_MIN_BOOT: Final = 50
TAX_ASSIGN_SPECIES: Final = False

# Prefisso della posizione da cui si deriva il modulo quando manca il file di
# arricchimento: tre lettere e una cifra, come i moduli della stazione nei nomi
# delle posizioni di OSD-734 (LAB1, NOD2, JLP1).
META_MODULE_REGEX: Final = r"^([A-Z]{3}[0-9])"

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

# Modalita' della decontaminazione (S12). Aggregata: le piastre di OSD-734 hanno
# da 5 a 8 controlli negativi, e il minimo di dieci probabilita' stimate
# ciascuna su cosi' pochi negativi, con la soglia 0,5, e' permissivo per
# costruzione (misurato: per piastra si toglierebbe la meta' delle letture dei
# biologici). La modalita' per piastra resta calcolata come diagnostica.
DECONTAM_MODE: Final = "aggregate"

# Colonna e etichette che distinguono controlli e campioni biologici.
CTRL_COLUMN: Final = "Characteristics[Material Type]"
CTRL_BLANK_VALUES: Final = ("blank control",)
CTRL_POSITIVE_VALUES: Final = ("Positive Control",)
CTRL_BIOLOGICAL_VALUES: Final = ("Surface swab",)
# Campioni da trattare come controlli negativi qualunque sia il materiale
# dichiarato, riconosciuti dal valore di una colonna della tabella campioni di
# studio. In OSD-734 i tamponi mai aperti ("Unopened 3DMM Swab Tube") sono
# dichiarati "Surface swab" ma non hanno campionato alcuna superficie: per
# funzione sono controlli di campo (bassa biomassa, composizione da reagente).
CTRL_BLANK_OVERRIDE_COLUMN: Final = "Factor Value[Sample Location]"
CTRL_BLANK_OVERRIDE_VALUES: Final = ("Unopened 3DMM Swab Tube",)

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

# Colonna del file di arricchimento con le cellule seminate in ciascun
# controllo positivo della serie di diluizione: e' il livello di diluizione.
KATHAROSEQ_CELL_COUNT_COLUMN: Final = "katharoseq_cell_count"

# Colonne dei metadati portate nell'oggetto integrato (S10) oltre a quelle
# dell'inventario. Dalla tabella campioni di studio: se il campione e' stato in
# volo, e come e' stato conservato. Dal file di arricchimento: identificativo
# della piastra, pozzetto, date di estrazione e di sequenziamento, piastra dei
# primer, e le cellule seminate nei controlli positivi della serie KatharoSeq,
# che servono alla calibrazione di S11.
OUT_STUDY_COLUMNS: Final = (
    "Factor Value[Spaceflight]",
    "Parameter Value[Sample Preservation Method]",
)
OUT_BATCH_COLUMNS: Final = (
    "extraction_plate_id",
    "well_id",
    "extraction_date",
    "primer_plate",
    "run_date",
    "katharoseq_cell_count",
)


class RiferimentoDataset(NamedTuple):
    """Il valore di un parametro scelto su OSD-734 e il fatto che lo giustifica."""

    #: Il valore adottato per il dataset di riferimento.
    valore: Any
    #: Il fatto accertato su OSD-734 da cui quel valore discende.
    fatto: str


#: I parametri il cui valore descrive il dataset di riferimento e non la
#: pipeline, ciascuno con il valore di OSD-734 e il fatto accertato che lo
#: giustifica. E' l'unica fonte della motivazione: il report la riporta accanto
#: a ogni parametro che in un'esecuzione coincide con il valore di OSD-734, e un
#: test la tiene allineata ai marcatori [OSD-734] di config.example.yaml. Un
#: valore ereditato in silenzio su un altro dataset e' il rischio che questo
#: elenco serve a rendere visibile.
FATTI_OSD734: Final[Mapping[str, RiferimentoDataset]] = MappingProxyType({
    "filter.truncLen": RiferimentoDataset(
        FILTER_TRUNCLEN,
        "le letture del dataset completo sono lunghe 137-151 basi (S1): 137 è la "
        "più corta, e un troncamento maggiore la scarterebbe",
    ),
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
    "meta.module_regex": RiferimentoDataset(
        META_MODULE_REGEX,
        "i nomi delle posizioni iniziano con il modulo della stazione, tre lettere "
        "e una cifra (LAB1, NOD2, JLP1)",
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
    "err.batch_column": RiferimentoDataset(
        ERR_BATCH_COLUMN,
        "nel file del lotto la corsa di sequenziamento sta nella colonna "
        "run_prefix (due corse)",
    ),
    "decontam.batch_column": RiferimentoDataset(
        DECONTAM_BATCH_COLUMN,
        "nel file del lotto la piastra di estrazione sta nella colonna "
        "extraction_plate_num (dieci piastre da 96)",
    ),
    "meta.module_column": RiferimentoDataset(
        META_MODULE_COLUMN,
        "la posizione di campionamento è dichiarata nella colonna Factor "
        "Value[Sample Location] della tabella campioni di studio",
    ),
    "meta.non_surface_positions": RiferimentoDataset(
        META_NON_SURFACE_POSITIONS,
        "i campioni d'aria e i tamponi mai aperti riportano una posizione che non "
        "è una superficie, e non appartengono a un modulo",
    ),
    "meta.batch_key_column": RiferimentoDataset(
        META_BATCH_KEY_COLUMN,
        "il file del lotto identifica i campioni per accession, nella colonna "
        "experiment_accession: i due campioni risequenziati vi si distinguono solo "
        "così",
    ),
    "meta.batch_module_column": RiferimentoDataset(
        META_BATCH_MODULE_COLUMN,
        "il file del lotto dichiara il modulo nella colonna module",
    ),
    "ctrl.column": RiferimentoDataset(
        CTRL_COLUMN,
        "la classe del campione è dichiarata nella colonna Characteristics[Material "
        "Type] della tabella campioni di studio",
    ),
    "ctrl.blank_values": RiferimentoDataset(
        CTRL_BLANK_VALUES,
        "i controlli negativi sono dichiarati con il materiale blank control",
    ),
    "ctrl.positive_values": RiferimentoDataset(
        CTRL_POSITIVE_VALUES,
        "i controlli positivi sono dichiarati con il materiale Positive Control",
    ),
    "ctrl.biological_values": RiferimentoDataset(
        CTRL_BIOLOGICAL_VALUES,
        "i campioni biologici sono dichiarati con il materiale Surface swab",
    ),
    "ctrl.blank_override_column": RiferimentoDataset(
        CTRL_BLANK_OVERRIDE_COLUMN,
        "i tamponi mai aperti si riconoscono dalla posizione dichiarata, nella "
        "colonna Factor Value[Sample Location]",
    ),
    "ctrl.blank_override_values": RiferimentoDataset(
        CTRL_BLANK_OVERRIDE_VALUES,
        "33 tamponi dichiarati Surface swab hanno posizione Unopened 3DMM Swab Tube: "
        "non hanno campionato alcuna superficie, e per funzione sono controlli "
        "negativi",
    ),
    "decontam.mode": RiferimentoDataset(
        DECONTAM_MODE,
        "le piastre hanno da 5 a 8 controlli negativi dichiarati: per piastra, con "
        "la soglia 0,5, si toglierebbe circa la metà delle letture dei biologici "
        "(misurato)",
    ),
    "tax.ref_name": RiferimentoDataset(
        TAX_REF_NAME,
        "il riferimento tassonomico adottato per lo studio è SILVA",
    ),
    "tax.ref_version": RiferimentoDataset(
        TAX_REF_VERSION,
        "la versione del riferimento adottata per lo studio è la 138",
    ),
    "katharoseq.target_taxon": RiferimentoDataset(
        KATHAROSEQ_TARGET_TAXON,
        "i controlli positivi contengono un solo ceppo, Variovorax sp. OAS795",
    ),
    "katharoseq.cell_count_column": RiferimentoDataset(
        KATHAROSEQ_CELL_COUNT_COLUMN,
        "il file del lotto riporta le cellule seminate in ciascun controllo "
        "positivo nella colonna katharoseq_cell_count",
    ),
    "qc.primer_sequence": RiferimentoDataset(
        QC_PRIMER_SEQUENCE,
        "la regione V4 è stata amplificata con il primer 515F",
    ),
    "qc.conserved_motif": RiferimentoDataset(
        QC_CONSERVED_MOTIF,
        "le letture iniziano con il motivo conservato a valle del 515F",
    ),
    "out.study_columns": RiferimentoDataset(
        OUT_STUDY_COLUMNS,
        "la tabella campioni di studio riporta se il campione è stato in volo e "
        "come è stato conservato",
    ),
    "out.batch_columns": RiferimentoDataset(
        OUT_BATCH_COLUMNS,
        "il file del lotto riporta piastra, pozzetto, date di estrazione e di "
        "sequenziamento, piastra dei primer e cellule seminate nei controlli "
        "positivi",
    ),
})

#: Chiavi il cui valore predefinito descrive il dataset di riferimento e non la
#: pipeline. Passando a un altro dataset vanno rivalutate una per una.
DERIVATI_DAL_DATASET: Final[tuple[str, ...]] = tuple(FATTI_OSD734)
