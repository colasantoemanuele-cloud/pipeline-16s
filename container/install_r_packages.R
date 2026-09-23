#!/usr/bin/env Rscript
# Installa i pacchetti R richiesti dalla pipeline.
#
# Eseguito durante la costruzione dell'immagine. Usa il repository di binari
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

binari <- sprintf(
  "https://bioconductor.org/packages/%s/container-binaries/bioconductor_docker",
  bioc_atteso
)

# Pacchetti dei calcoli scientifici della pipeline.
pipeline <- c(
  "dada2", "phyloseq", "DECIPHER", "decontam",
  "phangorn", "Biostrings", "ShortRead"
)

# Strumenti di servizio: renv produce il file di lock, jsonlite lo rilegge in
# fase di verifica.
servizio <- c("renv", "jsonlite")

BiocManager::install(
  c(pipeline, servizio),
  site_repository = binari,
  ask = FALSE,
  update = FALSE
)

# Senza questo controllo un'installazione fallita resterebbe un semplice avviso
# e l'immagine risulterebbe costruita ma inutilizzabile.
tutti <- c(pipeline, servizio)
mancanti <- tutti[!vapply(tutti, requireNamespace, logical(1), quietly = TRUE)]
if (length(mancanti)) {
  stop("pacchetti non installati: ", paste(mancanti, collapse = ", "))
}

cat("\nVersioni installate:\n")
for (p in tutti) {
  cat(sprintf("  %-12s %s\n", p, as.character(packageVersion(p))))
}
