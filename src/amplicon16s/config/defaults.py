"""Valori predefiniti dei parametri di configurazione.

I valori stanno qui e non dentro lo schema perché la domanda «da dove viene
questo valore?» deve potersi rispondere leggendo un solo file.

**Nessun predefinito viene da un dataset.** Un parametro ha un predefinito solo
in uno di questi casi, e ogni costante dice in quale ricade:

* **Valore standard del metodo** (:data:`STANDARD_DEL_METODO`): il predefinito
  del pacchetto che realizza il calcolo, o una raccomandazione pubblicata dai
  suoi autori, indipendente da qualunque dataset. Il commento ne cita la fonte.

* **Una sola forma realizzata**: il parametro ammette un solo valore, perché la
  pipeline non realizza le alternative. Non c'è una scelta da ereditare.

* **Parametro dello strumento**: governa come si esegue o in che forma si
  scrive, non che cosa si conclude dai dati (risorse, ripresa, tentativi, seme).

Tutti gli altri parametri dipendono dai dati o dallo studio (regione
amplificata e primer, lunghezza delle letture, biomassa, disegno dei controlli,
soglie di qualità) e sono **obbligatori** (:data:`OBBLIGATORI`): non hanno
alcun predefinito, vanno dichiarati tutti nella configurazione, anche vuoti o
nulli quando non sono pertinenti, e G15 non parte finché ne manca uno. Nel
dubbio un parametro è obbligatorio: ereditare in silenzio una scelta fatta per
un altro dataset è il rischio che questa regola toglie.
"""

from __future__ import annotations

from collections.abc import Mapping
from types import MappingProxyType
from typing import Final

# --------------------------------------------------------------------------- #
# Valori standard del metodo                                                   #
# --------------------------------------------------------------------------- #

# filter: filtraggio delle letture (dada2::filterAndTrim)
# Errori attesi massimi per lettura: il predefinito del pacchetto e' Inf (nessun
# filtro); 2 e' il valore dei "standard filtering parameters" del tutorial di
# dada2 (Callahan et al., https://benjjneb.github.io/dada2/tutorial.html).
FILTER_MAXEE: Final = 2.0
# Predefiniti di dada2::filterAndTrim (dada2 1.36): truncQ = 2, maxN = 0
# (l'inferenza non ammette basi indeterminate), rm.phix = TRUE.
FILTER_TRUNCQ: Final = 2
FILTER_MAXN: Final = 0
FILTER_RM_PHIX: Final = True

# err: apprendimento del modello di errore (dada2::learnErrors)
# Predefiniti di dada2::learnErrors (dada2 1.36): errorEstimationFunction =
# loessErrfun, nbases = 1e8, MAX_CONSIST = 10. La variante monotona della
# funzione serve alle qualita' raggruppate in pochi valori e va scelta.
ERR_ERROR_FUNCTION: Final = "loess"
ERR_NBASES: Final = 1.0e8
ERR_MAX_CONSIST: Final = 10

# dada: inferenza delle varianti (dada2::dada)
# Predefinito dell'opzione OMEGA_A di dada2 (dada2::getDadaOpt, dada2 1.36).
DADA_OMEGA_A: Final = 1.0e-40

# chimera: rimozione delle chimere (dada2::removeBimeraDenovo)
# Predefiniti di dada2 1.36: method = "consensus" (removeBimeraDenovo) e, per
# il metodo di consenso (isBimeraDenovoTable), minParentAbundance = 2,
# minSampleFraction = 0.9, allowOneOff = FALSE.
CHIMERA_METHOD: Final = "consensus"
CHIMERA_MIN_PARENT_ABUNDANCE: Final = 2
CHIMERA_MIN_SAMPLE_FRACTION: Final = 0.9
CHIMERA_ALLOW_ONE_OFF: Final = False

# tax: assegnazione tassonomica (dada2::assignTaxonomy)
# Predefinito di dada2::assignTaxonomy (dada2 1.36): minBoot = 50.
TAX_MIN_BOOT: Final = 50

# decontam: rimozione dei contaminanti (decontam::isContaminant)
# Come si combinano le probabilita' dei lotti: il predefinito di batch.combine
# in decontam::isContaminant (decontam 1.28) e' "minimum". Dichiararlo evita di
# dipendere da un predefinito che una versione nuova potrebbe cambiare.
DECONTAM_BATCH_COMBINE: Final = "minimum"

#: I parametri il cui predefinito e' il valore standard del metodo, con la
#: fonte. E' l'elenco che un test confronta con lo schema: un predefinito che
#: non e' qui, non e' una sola forma realizzata e non e' un parametro dello
#: strumento non e' ammesso.
STANDARD_DEL_METODO: Final[Mapping[str, str]] = MappingProxyType({
    "filter.maxEE": "tutorial di dada2, standard filtering parameters (maxEE = 2)",
    "filter.truncQ": "dada2::filterAndTrim, truncQ = 2",
    "filter.maxN": "dada2::filterAndTrim, maxN = 0",
    "filter.rm_phix": "dada2::filterAndTrim, rm.phix = TRUE",
    "err.error_function": "dada2::learnErrors, errorEstimationFunction = loessErrfun",
    "err.nbases": "dada2::learnErrors, nbases = 1e8",
    "err.max_consist": "dada2::learnErrors, MAX_CONSIST = 10",
    "dada.omega_a": "dada2::setDadaOpt, OMEGA_A = 1e-40",
    "chimera.method": "dada2::removeBimeraDenovo, method = \"consensus\"",
    "chimera.min_parent_abundance": "dada2::isBimeraDenovoTable, minParentAbundance = 2",
    "chimera.min_sample_fraction": "dada2::isBimeraDenovoTable, minSampleFraction = 0.9",
    "chimera.allow_one_off": "dada2::isBimeraDenovoTable, allowOneOff = FALSE",
    "tax.min_boot": "dada2::assignTaxonomy, minBoot = 50",
    "decontam.batch_combine": "decontam::isContaminant, batch.combine = \"minimum\"",
})

# --------------------------------------------------------------------------- #
# Una sola forma realizzata                                                    #
# --------------------------------------------------------------------------- #
# Insiemi chiusi con un solo valore ammesso: vanno estesi, e resi una scelta,
# quando la pipeline realizza un'alternativa.

TAX_CLASSIFIER: Final = "naive_bayes"
# L'assegnazione della specie non e' realizzata: il solo valore ammesso e' falso.
TAX_ASSIGN_SPECIES: Final = False
# L'allineatore multiplo e il modello evolutivo della massima verosimiglianza.
PHYLO_ALIGNER: Final = "decipher"
PHYLO_MODEL: Final = "GTR+G+I"
KATHAROSEQ_CURVE_MODEL: Final = "allosteric_sigmoid"
# Stadio delle letture su cui si misurano profondita' e fedelta', e a cui si
# applica la soglia derivata.
KATHAROSEQ_READ_STAGE: Final = "nonchimeric"
# La sola decontaminazione realizzata e' quella per prevalenza.
DECONTAM_METHOD: Final = "prevalence"
OUT_SERIALIZATION: Final = "rds"
OUT_ASV_ID_SCHEME: Final = "abundance_rank"
# Le fasi dopo l'oggetto integrato lavorano con le varianti sulle righe.
OUT_TAXA_ARE_ROWS: Final = True
# Identificativi dei campioni nell'oggetto integrato: l'accession, la chiave
# del join in tutta la pipeline. Il nome del campione resta fra i metadati.
OUT_SAMPLE_ID_SOURCE: Final = "accession"

# --------------------------------------------------------------------------- #
# Parametri dello strumento                                                    #
# --------------------------------------------------------------------------- #

# asv: selezione delle sequenze
# Tolleranza sulla lunghezza delle varianti: nessuna. Le letture single-end
# sono troncate tutte alla stessa lunghezza e non vengono unite, quindi ogni
# variante ha per costruzione quella lunghezza meno le basi tolte in testa.
ASV_LEN_TOL: Final = 0

# tax
# Elenco facoltativo dei taxa con un difetto noto del riferimento: assente.
TAX_REF_BAD_TAXA: Final = None

# retry: nuovi tentativi sulle fasi fallite
RETRY_ENABLED: Final = True
RETRY_MAX_ATTEMPTS: Final = 2
# Elenco chiuso dei codici di errore per cui e' ammesso un nuovo tentativo
# automatico. Il criterio non dipende dal dataset: si ritenta solo dove
# l'azione correttiva non modifica alcuna assunzione dell'analisi, e vale
# quindi su qualunque dataset.
RETRY_WHITELIST: Final = ("E-S2-03", "E-S4-02")

# out: forma degli artefatti prodotti
# Gli export in testo accanto all'oggetto: stessi dati, altra forma.
OUT_EXPORT_FLAT: Final = True

# run: esecuzione
# File di blocco delle versioni dei pacchetti R. Dichiararlo nella
# configurazione lo fa entrare nel digest: senza, due esecuzioni non sarebbero
# confrontabili rispetto alle versioni R impiegate.
RUN_LOCKFILE: Final = "renv.lock"
# Il seme dei generatori di numeri casuali: un valore arbitrario, fisso perche'
# due esecuzioni diano gli stessi risultati. Non descrive i dati.
RUN_SEED: Final = 100
# run.threads non ha un valore fisso: per difetto e' nullo, cioe' automatico, e i
# processori utilizzabili si contano all'uso (config/schema.py, thread_effettivi).
RUN_BATCH_SIZE: Final = 24
RUN_KEEP_FILTERED_FASTQ: Final = True
# La regola rigorosa sulla provenienza: disattivata per difetto, perche' lo
# sviluppo e i test girano su alberi di lavoro modificati; attiva nelle
# configurazioni di un'esecuzione da consegnare.
RUN_STRICT_PROVENANCE: Final = False


# --------------------------------------------------------------------------- #
# Parametri obbligatori: nessun predefinito                                    #
# --------------------------------------------------------------------------- #

#: I parametri che dipendono dai dati o dallo studio e non hanno predefinito.
#: Vanno dichiarati tutti nel file di configurazione; quelli non pertinenti al
#: dataset si dichiarano nulli o vuoti dove lo schema lo ammette (per esempio
#: le etichette dei controlli positivi in un dataset che non ne ha). G15
#: respinge la configurazione elencando quelli che mancano (E-G15-10). Che
#: cosa significa ciascuno e come si sceglie e' in config/config.example.yaml.
OBBLIGATORI: Final[tuple[str, ...]] = (
    # come si riconoscono file e campioni
    "io.fastq_glob",
    "io.accession_regex",
    "meta.sample_id_column",
    "meta.accession_column",
    # il formato dei metadati
    "meta.derive_module",
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
    # la regione amplificata, le letture e il riferimento
    "filter.truncLen",
    "filter.truncLen_shortfall_warn",
    "filter.trimLeft",
    "qc.primer_sequence",
    "qc.conserved_motif",
    "tax.ref_name",
    "tax.ref_version",
    "tax.try_rc",
    # le scelte di analisi che dipendono dallo studio
    "err.randomize",
    "dada.pool",
    "chimera.min_fold_parent_over_abundance",
    "filt.remove_na_phylum",
    "filt.exclude_taxa",
    "phylo.enabled",
    "phylo.max_seqs",
    "prev.min_fraction",
    "prev.min_count",
    "prev.apply",
    # i controlli sperimentali
    "ctrl.min_positives",
    "ctrl.min_positive_pass_frac",
    "ctrl.positive_gate",
    "katharoseq.target_taxon",
    "katharoseq.collapse_rank",
    "katharoseq.target_sensitivity",
    "katharoseq.min_r2",
    "decontam.threshold",
    "decontam.min_blanks",
    "decontam.mode",
    # le soglie di qualita'
    "qc.head_reads",
    "qc.max_primer_hit_frac",
    "qc.min_motif_frac",
    "qc.max_frac_short_reads",
    "qc.max_zeroed_samples",
    "qc.max_frac_lost_filter",
    "qc.max_asv_count",
    "qc.warn_frac_chimeric",
    "qc.stop_frac_chimeric",
    "qc.min_frac_phylum",
    "qc.min_reads_mode",
    "qc.max_frac_contaminant",
    "qc.min_reads_final",
    "qc.min_frac_reads_retained",
)
