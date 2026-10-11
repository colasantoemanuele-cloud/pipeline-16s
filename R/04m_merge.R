# SM - fusione delle coppie (dada2::mergePairs), solo in modalita' appaiata.
#
# Legge, per ogni campione, cio' che S4 ha conservato delle due letture
# (fusione_forward_<accession>.rds e fusione_reverse_<accession>.rds: la mappa
# dalle letture alle sequenze uniche, quella dalle sequenze uniche alle
# varianti e la tabella delle varianti) e scrive in 05_pair_merge/:
#
#   varianti_per_campione.rds  le sequenze fuse di ogni campione: lista con
#                              nome (accession) di vettori interi con nome
#                              (sequenza -> coppie), nella forma di quello di
#                              S4; e' l'ingresso della tabella delle sequenze.
#                              I campioni senza coppie fuse non vi compaiono
#   fusione.tsv                per campione: coppie entrate nella fusione
#                              (entrambe le letture attribuite a una variante),
#                              coppie fuse, quota, sequenze fuse distinte
#   lunghezze_fuse.tsv         per lunghezza: sequenze fuse distinte e coppie
#   letture_fuse.tsv           coppie fuse per campione (tracciamento)
#
# Parametri: forward e reverse (accession -> file di S4), senza_coppie
# (accession senza coppie filtrate, zero ovunque), min_overlap, max_mismatch,
# trim_overhang.
#
# Le letture filtrate non si rileggono: lo script legge i soli file di S4.
# mergePairs riceve le parti degli oggetti dada e derep che usa; il test di
# equivalenza con gli oggetti interi lo sorveglia.
#
# trimOverhang toglie dalla sequenza fusa le code di una lettura oltre
# l'inizio dell'altra: con letture piu' lunghe dell'amplicone sono il
# complemento del primer opposto e dell'adattatore.

for (f in c("io_json.R", "errors.R", "letture.R")) {
  source(file.path(Sys.getenv("AMPLICON16S_R_LIB"), f))
}

esegui_fase(function(parametri, cartella) {
  richiedi_pacchetti("dada2")

  forward <- unlist(parametri$forward)
  reverse <- unlist(parametri$reverse)[names(forward)]
  senza_coppie <- as.character(unlist(parametri$senza_coppie))
  if (anyNA(reverse)) {
    stop("campioni senza le varianti della lettura reverse: ",
         paste(names(forward)[is.na(reverse)], collapse = ", "))
  }

  varianti <- list()
  ingresso <- numeric(0)
  fuse <- numeric(0)
  distinte <- integer(0)
  for (a in sort(names(forward), method = "radix")) {
    f <- readRDS(forward[[a]])
    r <- readRDS(reverse[[a]])
    # Il messaggio di mergePairs sulle sequenze fuse ripetute si tace: coppie
    # diverse di varianti possono dare la stessa sequenza, e qui si sommano.
    coppie <- suppressMessages(dada2::mergePairs(
      f$dada, f$derep, r$dada, r$derep,
      minOverlap = as.integer(parametri$min_overlap),
      maxMismatch = as.integer(parametri$max_mismatch),
      trimOverhang = isTRUE(parametri$trim_overhang),
      returnRejects = TRUE,
      verbose = FALSE
    ))
    ingresso[[a]] <- sum(as.numeric(coppie$abundance))
    accettate <- coppie[coppie$accept, , drop = FALSE]
    if (nrow(accettate)) {
      somme <- rowsum(as.numeric(accettate$abundance), accettate$sequence, reorder = FALSE)
      per_sequenza <- stats::setNames(as.integer(somme[, 1L]), rownames(somme))
      ordine <- order(-per_sequenza, names(per_sequenza), method = "radix")
      varianti[[a]] <- per_sequenza[ordine]
    }
    fuse[[a]] <- sum(as.numeric(accettate$abundance))
    distinte[[a]] <- length(unique(accettate$sequence))
  }
  salva_rds(varianti, file.path(cartella, "varianti_per_campione.rds"))

  campioni <- names(ingresso)
  scrivi_atomico(
    c("campione\tcoppie\tfuse\tquota\tsequenze",
      sprintf("%s\t%.0f\t%.0f\t%s\t%d", campioni, ingresso[campioni], fuse[campioni],
              ifelse(ingresso[campioni] > 0,
                     sprintf("%.6f", fuse[campioni] / ingresso[campioni]), ""),
              distinte[campioni])),
    file.path(cartella, "fusione.tsv")
  )

  # La distribuzione delle lunghezze, sulle sequenze distinte e sulle coppie.
  sequenze <- unlist(lapply(varianti, names), use.names = FALSE)
  conteggi <- unlist(lapply(varianti, unname), use.names = FALSE)
  lunghezze <- nchar(sequenze)
  per_lunghezza <- sort(unique(lunghezze))
  scrivi_atomico(
    c("lunghezza\tsequenze\tcoppie",
      vapply(per_lunghezza, function(l) {
        sprintf("%d\t%d\t%.0f", l, length(unique(sequenze[lunghezze == l])),
                sum(as.numeric(conteggi[lunghezze == l])))
      }, character(1))),
    file.path(cartella, "lunghezze_fuse.tsv")
  )

  tracciate <- c(fuse, stats::setNames(numeric(length(senza_coppie)), senza_coppie))
  tracciate <- tracciate[sort(names(tracciate), method = "radix")]
  c("varianti_per_campione.rds", "fusione.tsv", "lunghezze_fuse.tsv",
    traccia_letture(tracciate, "fuse", cartella))
})
