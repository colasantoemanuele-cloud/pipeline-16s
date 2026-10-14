#!/usr/bin/env Rscript
# Installa i pacchetti R richiesti dalla pipeline.
#
# Eseguito durante la costruzione dell'immagine (container/Dockerfile), con le
# variabili BIOCONDUCTOR_VERSION e CRAN_SNAPSHOT; non scrive file oltre alle
# librerie installate. Lo stesso elenco di pacchetti e' in
# container/snapshot_renv_lock.R, che ne registra le versioni.
#
# Usa il repository di binari
# precompilati per questa immagine Bioconductor: evita di ricompilare da
# sorgente e rende l'esito determinato dalla release Bioconductor anziché dallo
# stato del compilatore locale.

options(warn = 1)

bioc_atteso <- Sys.getenv("BIOCONDUCTOR_VERSION")
if (!nzchar(bioc_atteso)) {
  stop("variabile BIOCONDUCTOR_VERSION non impostata")
}

bioc_installato <- as.character(BiocManager::version())
if (!identical(bioc_installato, bioc_atteso)) {
  stop("release Bioconductor attesa ", bioc_atteso, ", trovata ", bioc_installato)
}

# Le dipendenze da CRAN vengono da un'istantanea datata e non da "latest", che
# cambia nel tempo: una build da zero darebbe altrimenti versioni diverse da
# renv.lock. I pacchetti gia' presenti nell'immagine di partenza, ancorata per
# digest, non vengono aggiornati (update = FALSE).
istantanea <- Sys.getenv("CRAN_SNAPSHOT")
if (!nzchar(istantanea)) {
  stop("variabile CRAN_SNAPSHOT non impostata")
}
options(repos = c(CRAN = istantanea))

binari <- sprintf(
  "https://bioconductor.org/packages/%s/container-binaries/bioconductor_docker",
  bioc_atteso
)

# Pacchetti dei calcoli scientifici della pipeline.
pipeline <- c(
  "dada2", "phyloseq", "DECIPHER", "decontam",
  "phangorn", "Biostrings", "ShortRead"
)

# Pacchetti delle analisi a valle. vegan (ordinamento e diversita') e' usato
# dal codice delle analisi ecologiche dello strumento (amplicon16s_eco): e'
# gia' una dipendenza di phyloseq, ed e' elencato perche' quel codice lo usa
# direttamente. DESeq2 (abbondanza differenziale) e randomForest
# (classificazione) non sono usati dal codice dello strumento: restano
# nell'immagine per le analisi a valle dell'utente, che girano cosi' nello
# stesso ambiente bloccato della pipeline. Il codice R delle analisi
# ecologiche si esegue dal repository montato, non dall'immagine.
ecologia <- c("DESeq2", "vegan", "randomForest")

# Strumenti di servizio: renv produce il file di lock, jsonlite lo rilegge in
# fase di verifica.
servizio <- c("renv", "jsonlite")

BiocManager::install(
  c(pipeline, ecologia, servizio),
  site_repository = binari,
  ask = FALSE,
  update = FALSE
)

# Senza questo controllo un'installazione fallita resterebbe un semplice avviso
# e l'immagine risulterebbe costruita ma inutilizzabile.
tutti <- c(pipeline, ecologia, servizio)
mancanti <- tutti[!vapply(tutti, requireNamespace, logical(1), quietly = TRUE)]
if (length(mancanti)) {
  stop("pacchetti non installati: ", paste(mancanti, collapse = ", "))
}

cat("\nVersioni installate:\n")
for (p in tutti) {
  cat(sprintf("  %-12s %s\n", p, as.character(packageVersion(p))))
}
