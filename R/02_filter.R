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
# maxN, rm_phix, processi, lotto (quanti file per chiamata al filtro).
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
  if (!requireNamespace("dada2", quietly = TRUE)) {
    stop("pacchetto dada2 non disponibile: la fase va eseguita nell'ambiente ",
         "della pipeline, il container")
  }

  campioni <- unlist(parametri$campioni)
  lotto <- as.integer(parametri$lotto)
  uscite <- file.path(cartella, paste0(names(campioni), "_filt.fastq.gz"))
  names(uscite) <- names(campioni)
  indici <- split(seq_along(campioni), ceiling(seq_along(campioni) / lotto))

  conteggi <- NULL
  for (numero in seq_along(indici)) {
    i <- indici[[numero]]
    esito <- tryCatch(
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
      ),
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

  ingresso <- conteggi[, "reads.in"]
  uscita <- conteggi[, "reads.out"]
  artefatti <- c(
    traccia_letture(ingresso, "prefiltro", cartella),
    traccia_letture(uscita, "filtrate", cartella)
  )
  # filterAndTrim non scrive il file di un campione azzerato.
  scritti <- names(uscite)[uscita > 0 & file.exists(uscite)]
  c(artefatti, basename(uscite[scritti]))
})
