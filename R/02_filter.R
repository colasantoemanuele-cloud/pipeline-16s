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
# Con il primer cercato lungo la lettura arriva anche taglio_primer: primer,
# discordanze_massime, sfasamento_massimo, scarta (e primer_reverse in modalita'
# appaiata), e lavoro, una cartella della fase per le letture tagliate. Prima
# del filtro ogni file (o coppia di file) si taglia dove il primer finisce
# (R/lib/primer.R), e il filtro legge le letture tagliate: truncLen e' allora
# la lunghezza dopo il taglio. Si scrivono in piu':
#
#   taglio_primer.tsv          per campione: letture (o coppie) in ingresso,
#                              quante portano il primer (per la reverse:
#                              con_primer_reverse), in quante manca almeno un
#                              primer, quante restano dopo il taglio
#   posizioni_primer.tsv       per campione e lettura: a quale sfasamento e con
#                              quante discordanze il primer e' stato trovato
#   letture_taglio_primer.tsv  le letture che restano dopo il taglio, nel
#                              formato del tracciamento
#
# letture_prefiltro.tsv resta l'ingresso della fase: le letture grezze.
#
# I file si filtrano a lotti. Un errore di lettura o di decompressione su un
# lotto si dichiara E-S2-03, ammesso al retry: se l'errore era transitorio il
# lotto si rilegge; se il file e' corrotto, fallira' di nuovo.

for (f in c("io_json.R", "errors.R", "letture.R", "primer.R")) {
  source(file.path(Sys.getenv("AMPLICON16S_R_LIB"), f))
}

# Il taglio dei primer cercati: scrive nella cartella di lavoro le letture
# tagliate di ogni campione e le tre tabelle, e restituisce i file da filtrare
# al posto degli originali, con i conteggi in ingresso.
taglia_i_primer <- function(parametri, cartella, campioni, reverse) {
  richiedi_pacchetti(c("ShortRead", "Biostrings"))
  taglio <- parametri$taglio_primer
  lavoro <- taglio$lavoro
  dir.create(lavoro, showWarnings = FALSE, recursive = TRUE)
  ricerca <- function(primer) list(
    primer = primer,
    discordanze_massime = as.integer(taglio$discordanze_massime),
    sfasamento_massimo = as.integer(taglio$sfasamento_massimo)
  )
  ricerche <- list(ricerca(taglio$primer))
  if (!is.null(reverse)) ricerche <- c(ricerche, list(ricerca(taglio$primer_reverse)))
  tagliati <- file.path(lavoro, paste0(names(campioni), "_tagliate.fastq.gz"))
  tagliati_reverse <- file.path(lavoro, paste0(names(campioni), "_tagliate_reverse.fastq.gz"))
  esiti <- parallel::mclapply(seq_along(campioni), function(i) {
    ingressi <- c(campioni[[i]], if (!is.null(reverse)) reverse[[i]])
    uscite <- c(tagliati[i], if (!is.null(reverse)) tagliati_reverse[i])
    taglia_file(ingressi, uscite, ricerche, isTRUE(taglio$scarta))
  }, mc.cores = as.integer(parametri$processi), mc.preschedule = FALSE)
  verifica_figli(esiti, names(campioni))
  guasti <- vapply(esiti, inherits, logical(1), what = "try-error")
  if (any(guasti)) {
    messaggio <- as.character(esiti[guasti][[1L]])
    if (grepl(ERRORE_DI_LETTURA, messaggio, ignore.case = TRUE)) {
      errore_catalogo("E-S2-03", sprintf(
        "taglio del primer (%s): %s",
        paste(names(campioni)[guasti], collapse = ", "), messaggio))
    }
    stop("taglio del primer non riuscito per ",
         paste(names(campioni)[guasti], collapse = ", "), ": ", messaggio)
  }
  names(esiti) <- names(campioni)
  list(esiti = esiti, tagliati = stats::setNames(tagliati, names(campioni)),
       tagliati_reverse = stats::setNames(tagliati_reverse, names(campioni)))
}

# Le tabelle del taglio, nell'ordine dei campioni dell'esecuzione: un campione
# senza letture ha zero in ogni colonna.
scrivi_tabelle_del_taglio <- function(esiti, ordine, appaiata, cartella) {
  valore <- function(campione, f) if (is.null(esiti[[campione]])) 0 else f(esiti[[campione]])
  letture <- vapply(ordine, valore, numeric(1), f = function(e) e$letture)
  rimaste <- vapply(ordine, valore, numeric(1), f = function(e) e$rimaste)
  senza <- vapply(ordine, valore, numeric(1), f = function(e) e$senza_primer)
  con <- vapply(ordine, valore, numeric(1), f = function(e) e$con_primer[1L])
  intestazione <- "campione\tletture\tcon_primer"
  righe <- sprintf("%s\t%.0f\t%.0f", ordine, letture, con)
  if (appaiata) {
    con_reverse <- vapply(ordine, valore, numeric(1), f = function(e) e$con_primer[2L])
    intestazione <- paste0(intestazione, "\tcon_primer_reverse")
    righe <- sprintf("%s\t%.0f", righe, con_reverse)
  }
  scrivi_atomico(
    c(paste0(intestazione, "\tsenza_primer\trimaste"),
      sprintf("%s\t%.0f\t%.0f", righe, senza, rimaste)),
    file.path(cartella, "taglio_primer.tsv")
  )
  righe <- "campione\tlettura\tsfasamento\tdiscordanze\tletture"
  for (campione in ordine) {
    if (is.null(esiti[[campione]])) next
    for (i in seq_along(esiti[[campione]]$posizioni)) {
      posizioni <- esiti[[campione]]$posizioni[[i]]
      if (!length(posizioni)) next
      parti <- do.call(rbind, strsplit(names(posizioni), ":", fixed = TRUE))
      sfasamento <- as.integer(parti[, 1L])
      discordanze <- as.integer(parti[, 2L])
      ordinati <- order(sfasamento, discordanze)
      righe <- c(righe, sprintf(
        "%s\t%s\t%d\t%d\t%.0f", campione, c("forward", "reverse")[i],
        sfasamento[ordinati], discordanze[ordinati], posizioni[ordinati]))
    }
  }
  scrivi_atomico(righe, file.path(cartella, "posizioni_primer.tsv"))
  names(rimaste) <- ordine
  c("taglio_primer.tsv", "posizioni_primer.tsv",
    traccia_letture(rimaste, "taglio_primer", cartella))
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
  reverse <- NULL
  if (appaiata) {
    reverse <- unlist(parametri$reverse)[names(campioni)]
    if (anyNA(reverse)) {
      stop("campioni senza il file reverse: ",
           paste(names(campioni)[is.na(reverse)], collapse = ", "))
    }
    uscite_reverse <- file.path(cartella, paste0(names(campioni), "_filt_reverse.fastq.gz"))
    names(uscite_reverse) <- names(campioni)
  }
  # Con il primer cercato il filtro legge le letture tagliate; un campione a
  # cui il taglio non lascia letture non arriva al filtro.
  taglio <- facoltativo(parametri, "taglio_primer")
  artefatti_del_taglio <- character()
  grezze <- NULL
  if (!is.null(taglio)) {
    on.exit(unlink(taglio$lavoro, recursive = TRUE), add = TRUE)
    tagliati <- taglia_i_primer(parametri, cartella, campioni, if (appaiata) reverse else NULL)
    ordine <- as.character(unlist(facoltativo(parametri, "ordine", as.list(names(campioni)))))
    artefatti_del_taglio <- scrivi_tabelle_del_taglio(tagliati$esiti, ordine, appaiata, cartella)
    grezze <- vapply(tagliati$esiti, function(e) e$letture, numeric(1))
    rimaste <- vapply(tagliati$esiti, function(e) e$rimaste, numeric(1))
    senza_resto <- names(campioni)[rimaste == 0]
    campioni <- tagliati$tagliati[setdiff(names(campioni), senza_resto)]
    if (appaiata) reverse <- tagliati$tagliati_reverse[names(campioni)]
    uscite <- uscite[names(campioni)]
    if (appaiata) uscite_reverse <- uscite_reverse[names(campioni)]
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

  if (is.null(conteggi)) {
    # Il taglio del primer non ha lasciato letture ad alcun campione.
    conteggi <- matrix(numeric(), nrow = 0L, ncol = 2L,
                       dimnames = list(NULL, c("reads.in", "reads.out")))
  }
  # Con un solo campione l'estrazione di una colonna perde il nome della riga.
  ingresso <- stats::setNames(conteggi[, "reads.in"], rownames(conteggi))
  uscita <- stats::setNames(conteggi[, "reads.out"], rownames(conteggi))
  if (!is.null(grezze)) {
    # L'ingresso della fase sono le letture grezze, non quelle rimaste dopo il
    # taglio; un campione a cui il taglio non lascia nulla esce con zero.
    uscita <- c(uscita, stats::setNames(
      numeric(length(grezze) - length(uscita)), setdiff(names(grezze), names(uscita))))
    uscita <- uscita[names(grezze)]
    ingresso <- grezze
  }
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
    return(c(artefatti, artefatti_del_taglio,
             basename(uscite[scritti]), basename(uscite_reverse[scritti])))
  }
  c(artefatti, artefatti_del_taglio, basename(uscite[scritti]))
})
