# Addestra con DECIPHER::LearnTaxa i classificatori sintetici di prova per
# IdTaxa, dai riferimenti sintetici delle fixture.
#
#   Rscript --vanilla tests/fixtures/idtaxa/costruisci.R      (dalla radice del clone)
#
# Scrive, in questa cartella, un file .RData per riferimento, con l'oggetto
# trainingSet nella forma dei classificatori pubblicati per DECIPHER: i lignaggi
# cominciano con Root (rango rootrank) e ogni taxon ha il nome del suo rango,
# che qui si ricava dal prefisso del nome sintetico. Servono a provare la
# corrispondenza dei ranghi per nome:
#
# - i generi il cui lignaggio non ha l'ordine (la famiglia discende dalla
#   classe) restano tali, come i lignaggi di SILVA che saltano un rango;
# - nei generi di RANGO_IN_PIU il lignaggio ha, sotto la classe, un rango che
#   la tassonomia della pipeline non ha (subclass);
# - un nome che nel riferimento compare sotto due genitori diversi (LearnTaxa
#   non lo ammette) prende, sotto il secondo, il suffisso del genitore.
#
# I nomi dei taxa sono sintetici: il classificatore non dice nulla della
# composizione reale di alcun dataset.

suppressPackageStartupMessages({
  library(DECIPHER)
  library(Biostrings)
})

RANGO_DEL_PREFISSO <- c(Root = "rootrank", Bacteria = "domain", Phylum = "phylum",
                        Classe = "class", Sottoclasse = "subclass", Ordine = "order",
                        Famiglia = "family", Genere = "genus")
RANGO_IN_PIU <- c("Genere_08", "Genere_09")

addestra <- function(riferimento, uscita) {
  sequenze <- readDNAStringSet(riferimento)
  lignaggi <- strsplit(sub(";$", "", names(sequenze)), ";", fixed = TRUE)

  voci <- list()      # lignaggio del taxon -> c(nome, genitore, livello, rango)
  genitore_di <- list()
  tassonomia <- character(length(lignaggi))
  for (i in seq_along(lignaggi)) {
    nomi <- c("Root", lignaggi[[i]])
    if (nomi[length(nomi)] %in% RANGO_IN_PIU) {
      dopo <- grep("^Classe_", nomi)
      nomi <- append(nomi, sub("^Classe", "Sottoclasse", nomi[dopo]), after = dopo)
    }
    genitore <- ""
    for (livello in seq_along(nomi)) {
      # Lo stesso nome sotto un altro genitore e' un altro taxon.
      if (is.null(genitore_di[[nomi[livello]]])) genitore_di[[nomi[livello]]] <- genitore
      if (genitore_di[[nomi[livello]]] != genitore) {
        nomi[livello] <- paste0(nomi[livello], "_di_", nomi[livello - 1L])
      }
      chiave <- paste(nomi[seq_len(livello)], collapse = "; ")
      if (is.null(voci[[chiave]])) {
        voci[[chiave]] <- c(nomi[livello], genitore, livello - 1L,
                            RANGO_DEL_PREFISSO[[sub("_.*$", "", nomi[livello])]])
      }
      genitore <- chiave
    }
    tassonomia[i] <- paste(nomi, collapse = "; ")
  }
  chiavi <- names(voci)
  rango <- data.frame(
    Index = seq_along(chiavi),
    Name = vapply(voci, `[`, character(1), 1L),
    Parent = match(vapply(voci, `[`, character(1), 2L), chiavi, nomatch = 0L),
    Level = as.integer(vapply(voci, `[`, character(1), 3L)),
    Rank = vapply(voci, `[`, character(1), 4L),
    stringsAsFactors = FALSE, row.names = NULL
  )
  set.seed(100L)
  trainingSet <- LearnTaxa(sequenze, tassonomia, rank = rango, verbose = FALSE)
  save(trainingSet, file = uscita, compress = "xz")
  cat(basename(uscita), ":", length(sequenze), "sequenze,", nrow(rango), "taxa, ranghi",
      paste(sort(unique(trainingSet$ranks)), collapse = ", "), "\n")
}

qui <- "tests/fixtures/idtaxa"
addestra("tests/fixtures/riferimento_sintetico/riferimento.fa.gz",
         file.path(qui, "addestrato_sintetico.RData"))
addestra("tests/fixtures/osd276/appaiato/riferimento/riferimento.fa.gz",
         file.path(qui, "addestrato_appaiato.RData"))
