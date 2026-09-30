#!/usr/bin/env Rscript
# Genera renv.lock dalle versioni effettivamente installate nell'immagine.
#
# Il file di lock non elenca versioni desiderate: registra ciò che
# l'installazione ha realmente prodotto. Per questo viene prodotto dentro il
# container, a installazione conclusa, e non scritto a mano.
#
# Oltre ai pacchetti della pipeline, renv registra l'intera chiusura delle loro
# dipendenze: è quella, non i sette nomi, a determinare il risultato di un
# calcolo.

destinazione <- Sys.getenv("RENV_LOCK_OUT", "/out/renv.lock")

pipeline <- c(
  "dada2", "phyloseq", "DECIPHER", "decontam",
  "phangorn", "Biostrings", "ShortRead"
)

mancanti <- pipeline[!vapply(pipeline, requireNamespace, logical(1), quietly = TRUE)]
if (length(mancanti)) {
  stop("pacchetti non installati nell'immagine: ", paste(mancanti, collapse = ", "))
}

progetto <- tempfile("renv-snapshot-")
dir.create(progetto)

renv::snapshot(
  project  = progetto,
  library  = .libPaths(),
  lockfile = destinazione,
  packages = pipeline,
  prompt   = FALSE,
  force    = TRUE
)

# I pacchetti compilati con una correzione del progetto (per ora dada2, vedi
# container/dada2/) portano nel DESCRIPTION la versione di partenza, il nome
# della correzione e il suo SHA-256: il lock li registra in un blocco Patch, che
# container/verify_renv_lock.R confronta con l'immagine.
lock <- renv::lockfile_read(destinazione)
for (nome in names(lock$Packages)) {
  d <- utils::packageDescription(nome)
  if (!is.null(d[["Amplicon16sPatch"]])) {
    # renv lo attribuirebbe a Bioconductor dai campi git_* del DESCRIPTION: e'
    # invece compilato qui, dai sorgenti ufficiali con la correzione.
    lock$Packages[[nome]]$Source <- "Local"
    lock$Packages[[nome]]$Patch <- list(
      Base = d[["Amplicon16sBase"]],
      File = d[["Amplicon16sPatch"]],
      SHA256 = d[["Amplicon16sPatchSHA256"]]
    )
  }
}
renv::lockfile_write(lock, destinazione)

registrati <- jsonlite::fromJSON(destinazione, simplifyVector = FALSE)$Packages
cat(sprintf("\nrenv.lock scritto in %s: %d pacchetti registrati\n",
            destinazione, length(registrati)))
for (p in pipeline) {
  cat(sprintf("  %-12s %s\n", p, registrati[[p]]$Version))
}
