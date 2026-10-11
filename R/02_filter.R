# S2 - filtro e troncamento delle letture (dada2::filterAndTrim).
#
# Scrive in 03_filtered/:
#
#   <accession>_filt.fastq.gz  le letture conservate, troncate a truncLen; un
#                              campione che non ne conserva nessuna non ha file
#   letture_prefiltro.tsv      letture in ingresso al filtro, per campione
#   letture_filtrate.tsv       letture conservate, per campione
#
# Le due tabelle hanno il formato del tracciamento delle letture
# (R/lib/letture.R): ogni fase scrive le proprie, nella propria cartella, e la
# tabella completa si ricompone leggendole tutte.
#
# Parametri: campioni (accession -> file), truncLen, trimLeft, maxEE, truncQ,
# maxN, rm_phix, processi, lotto (quanti file per chiamata al filtro); vuoti
# (accession dei file senza letture, che non si filtrano: zero in ingresso e
# in uscita) e ordine (tutti gli accession, nell'ordine delle tabelle).
#
# In modalita' appaiata arrivano anche reverse (accession -> file della lettura
# reverse), truncLen_reverse, trimLeft_reverse e maxEE_reverse: il filtro e'
# quello appaiato di filterAndTrim, che tiene o scarta le due letture di una
# coppia insieme e scrive anche <accession>_filt_reverse.fastq.gz. I conteggi
# sono coppie. La corrispondenza delle coppie fra i due file e' verificata
# prima, lettura per lettura, dalla fase del profilo: matchIDs resta FALSE.
#
# I file si filtrano a lotti. Un errore di lettura o di decompressione su un
# lotto si dichiara E-S2-03, ammesso al retry: se l'errore era transitorio il
# lotto si rilegge; se il file e' corrotto, fallira' di nuovo.

for (f in c("io_json.R", "errors.R", "letture.R")) {
  source(file.path(Sys.getenv("AMPLICON16S_R_LIB"), f))
}

# Messaggi con cui la lettura di un archivio fallisce: dati compressi non
# validi, file troncato, record FASTQ incompleti o malformati.
ERRORE_DI_LETTURA <- paste(
  c("gzip", "compress", "unexpected end", "end of file", "truncat",
    "incomplete", "does not start with", "cannot open"),
  collapse = "|"
)

esegui_fase(function(parametri, cartella) {
  richiedi_pacchetti("dada2")

  campioni <- unlist(parametri$campioni)
  lotto <- as.integer(parametri$lotto)
  uscite <- file.path(cartella, paste0(names(campioni), "_filt.fastq.gz"))
  names(uscite) <- names(campioni)
  appaiata <- !is.null(facoltativo(parametri, "reverse"))
  if (appaiata) {
    reverse <- unlist(parametri$reverse)[names(campioni)]
    if (anyNA(reverse)) {
      stop("campioni senza il file reverse: ",
           paste(names(campioni)[is.na(reverse)], collapse = ", "))
    }
    uscite_reverse <- file.path(cartella, paste0(names(campioni), "_filt_reverse.fastq.gz"))
    names(uscite_reverse) <- names(campioni)
  }
  indici <- split(seq_along(campioni), ceiling(seq_along(campioni) / lotto))

  conteggi <- NULL
  for (numero in seq_along(indici)) {
    i <- indici[[numero]]
    esito <- tryCatch(
      if (appaiata) {
        dada2::filterAndTrim(
          fwd = unname(campioni[i]),
          filt = unname(uscite[i]),
          rev = unname(reverse[i]),
          filt.rev = unname(uscite_reverse[i]),
          truncLen = c(parametri$truncLen, parametri$truncLen_reverse),
          trimLeft = c(parametri$trimLeft, parametri$trimLeft_reverse),
          maxEE = c(parametri$maxEE, parametri$maxEE_reverse),
          truncQ = parametri$truncQ,
          maxN = parametri$maxN,
          rm.phix = isTRUE(parametri$rm_phix),
          matchIDs = FALSE,
          compress = TRUE,
          multithread = as.integer(parametri$processi),
          verbose = FALSE
        )
      } else {
      dada2::filterAndTrim(
        fwd = unname(campioni[i]),
        filt = unname(uscite[i]),
        truncLen = parametri$truncLen,
        trimLeft = parametri$trimLeft,
        maxEE = parametri$maxEE,
        truncQ = parametri$truncQ,
        maxN = parametri$maxN,
        rm.phix = isTRUE(parametri$rm_phix),
        compress = TRUE,
        multithread = as.integer(parametri$processi),
        verbose = FALSE
      )
      },
      error = function(e) {
        if (grepl(ERRORE_DI_LETTURA, conditionMessage(e), ignore.case = TRUE)) {
          errore_catalogo(
            "E-S2-03",
            sprintf("lotto %d di %d (%s): %s", numero, length(indici),
                    paste(names(campioni)[i], collapse = ", "), conditionMessage(e))
          )
        }
        stop(e)
      }
    )
    rownames(esito) <- names(campioni)[i]
    conteggi <- rbind(conteggi, esito)
  }

  # Con un solo campione l'estrazione di una colonna perde il nome della riga.
  ingresso <- stats::setNames(conteggi[, "reads.in"], rownames(conteggi))
  uscita <- stats::setNames(conteggi[, "reads.out"], rownames(conteggi))
  # I file senza letture: il bianco che non ha amplificato nulla resta nelle
  # tabelle con zero letture, come ogni campione azzerato dal filtro.
  vuoti <- as.character(unlist(facoltativo(parametri, "vuoti", list())))
  if (length(vuoti)) {
    ordine <- as.character(unlist(parametri$ordine))
    zeri <- stats::setNames(numeric(length(vuoti)), vuoti)
    ingresso <- c(ingresso, zeri)[ordine]
    uscita <- c(uscita, zeri)[ordine]
  }
  artefatti <- c(
    traccia_letture(ingresso, "prefiltro", cartella),
    traccia_letture(uscita, "filtrate", cartella)
  )
  # filterAndTrim non scrive il file di un campione azzerato.
  scritti <- names(uscite)[uscita[names(uscite)] > 0 & file.exists(uscite)]
  if (appaiata) {
    senza_compagno <- scritti[!file.exists(uscite_reverse[scritti])]
    if (length(senza_compagno)) {
      stop("letture forward filtrate senza le reverse: ", paste(senza_compagno, collapse = ", "))
    }
    return(c(artefatti, basename(uscite[scritti]), basename(uscite_reverse[scritti])))
  }
  c(artefatti, basename(uscite[scritti]))
})
