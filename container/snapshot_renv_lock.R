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

registrati <- jsonlite::fromJSON(destinazione, simplifyVector = FALSE)$Packages
cat(sprintf("\nrenv.lock scritto in %s: %d pacchetti registrati\n",
            destinazione, length(registrati)))
for (p in pipeline) {
  cat(sprintf("  %-12s %s\n", p, registrati[[p]]$Version))
}
