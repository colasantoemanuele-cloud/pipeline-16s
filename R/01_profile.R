# S1 - profilo di qualita' e di lunghezza delle letture.
#
# Legge ogni file per intero, a blocchi, e scrive in 02_qc_profiles/:
#
#   lunghezze.tsv      campione, lunghezza, letture: la distribuzione delle
#                      lunghezze di ogni campione
#   qualita.tsv        campione, posizione, letture, media, p10, p25, mediana,
#                      p75, p90: il profilo di qualita' per posizione
#   letture_grezze.tsv il conteggio delle letture di ogni campione, nel formato
#                      del tracciamento delle letture (R/lib/letture.R)
#   riepilogo.json     lunghezza minima, massima e moda sull'intero insieme, e
#                      per campione; profilo di qualita' aggregato
#
# Parametri: campioni (accession -> percorso del file), processi (quanti file
# leggere in parallelo), blocco (letture per blocco di lettura).
#
# Le tabelle sono in formato lungo, una riga per campione e per valore: si
# leggono con qualunque programma, e 960 campioni restano tre file invece di
# tremila.

for (f in c("io_json.R", "errors.R", "letture.R")) {
  source(file.path(Sys.getenv("AMPLICON16S_R_LIB"), f))
}

# Punteggi Phred possibili: da 0 a 93 nella codifica Sanger.
PUNTEGGI <- 94L

profila <- function(percorso, blocco) {
  lettore <- ShortRead::FastqStreamer(percorso, n = blocco)
  on.exit(close(lettore))

  lunghezze <- integer()
  qualita <- matrix(0, nrow = 0L, ncol = PUNTEGGI)
  letture <- 0

  while (length(parte <- ShortRead::yield(lettore)) > 0L) {
    letture <- letture + length(parte)
    larghezze <- Biostrings::width(ShortRead::sread(parte))
    conteggio <- tabulate(larghezze, nbins = max(larghezze))
    if (length(conteggio) > length(lunghezze)) {
      lunghezze <- c(lunghezze, integer(length(conteggio) - length(lunghezze)))
    }
    lunghezze[seq_along(conteggio)] <- lunghezze[seq_along(conteggio)] + conteggio

    # Una riga per lettura, una colonna per posizione; NA oltre la fine delle
    # letture piu' corte.
    punteggi <- methods::as(Biostrings::quality(parte), "matrix")
    if (ncol(punteggi) > nrow(qualita)) {
      qualita <- rbind(
        qualita,
        matrix(0, nrow = ncol(punteggi) - nrow(qualita), ncol = PUNTEGGI)
      )
    }
    for (posizione in seq_len(ncol(punteggi))) {
      colonna <- punteggi[, posizione]
      colonna <- colonna[!is.na(colonna)]
      qualita[posizione, ] <- qualita[posizione, ] +
        tabulate(colonna + 1L, nbins = PUNTEGGI)
    }
  }
  list(letture = letture, lunghezze = lunghezze, qualita = qualita)
}

# Quantile di un istogramma dei punteggi: il primo punteggio la cui frequenza
# cumulata raggiunge la frazione richiesta.
quantile_istogramma <- function(conteggi, frazione) {
  totale <- sum(conteggi)
  if (totale == 0) return(NA_integer_)
  which(cumsum(conteggi) >= frazione * totale)[1L] - 1L
}

righe_qualita <- function(campione, qualita) {
  valori <- 0:(PUNTEGGI - 1L)
  vapply(seq_len(nrow(qualita)), function(posizione) {
    conteggi <- qualita[posizione, ]
    totale <- sum(conteggi)
    sprintf(
      "%s\t%d\t%.0f\t%.3f\t%d\t%d\t%d\t%d\t%d",
      campione, posizione, totale,
      if (totale > 0) sum(conteggi * valori) / totale else NA_real_,
      quantile_istogramma(conteggi, 0.10), quantile_istogramma(conteggi, 0.25),
      quantile_istogramma(conteggi, 0.50), quantile_istogramma(conteggi, 0.75),
      quantile_istogramma(conteggi, 0.90)
    )
  }, character(1))
}

moda <- function(lunghezze) which.max(lunghezze)

esegui_fase(function(parametri, cartella) {
  for (pacchetto in c("ShortRead", "Biostrings")) {
    if (!requireNamespace(pacchetto, quietly = TRUE)) {
      stop("pacchetto ", pacchetto, " non disponibile: la fase va eseguita ",
           "nell'ambiente della pipeline, il container")
    }
  }

  campioni <- unlist(parametri$campioni)
  processi <- if (is.null(parametri$processi)) 1L else as.integer(parametri$processi)
  blocco <- if (is.null(parametri$blocco)) 100000L else as.integer(parametri$blocco)
  mancanti <- campioni[!file.exists(campioni)]
  if (length(mancanti)) {
    stop("file di letture non trovati: ", paste(mancanti, collapse = ", "))
  }

  profili <- parallel::mclapply(
    campioni, profila, blocco = blocco,
    mc.cores = processi, mc.preschedule = FALSE
  )
  guasti <- vapply(profili, inherits, logical(1), what = "try-error")
  if (any(guasti)) {
    stop("profilo non calcolabile per ", paste(names(campioni)[guasti], collapse = ", "),
         ": ", as.character(profili[guasti][[1L]]))
  }
  names(profili) <- names(campioni)

  # Distribuzione delle lunghezze.
  righe <- "campione\tlunghezza\tletture"
  for (campione in names(profili)) {
    l <- profili[[campione]]$lunghezze
    presenti <- which(l > 0)
    righe <- c(righe, sprintf("%s\t%d\t%.0f", campione, presenti, l[presenti]))
  }
  scrivi_atomico(righe, file.path(cartella, "lunghezze.tsv"))

  # Profilo di qualita' per posizione.
  righe <- "campione\tposizione\tletture\tmedia\tp10\tp25\tmediana\tp75\tp90"
  for (campione in names(profili)) {
    righe <- c(righe, righe_qualita(campione, profili[[campione]]$qualita))
  }
  scrivi_atomico(righe, file.path(cartella, "qualita.tsv"))

  # Conteggio delle letture, nel formato del tracciamento.
  conteggi <- vapply(profili, function(p) p$letture, numeric(1))
  tracciamento <- traccia_letture(conteggi, "grezze", cartella)

  # Riepilogo sull'insieme.
  massimo <- max(vapply(profili, function(p) length(p$lunghezze), integer(1)))
  totale_l <- numeric(massimo)
  larghezza_q <- max(vapply(profili, function(p) nrow(p$qualita), integer(1)))
  totale_q <- matrix(0, nrow = larghezza_q, ncol = PUNTEGGI)
  per_campione <- list()
  for (campione in names(profili)) {
    p <- profili[[campione]]
    totale_l[seq_along(p$lunghezze)] <- totale_l[seq_along(p$lunghezze)] + p$lunghezze
    totale_q[seq_len(nrow(p$qualita)), ] <- totale_q[seq_len(nrow(p$qualita)), ] + p$qualita
    presenti <- which(p$lunghezze > 0)
    per_campione[[campione]] <- list(
      letture = p$letture,
      lunghezza_minima = if (length(presenti)) min(presenti) else NULL,
      lunghezza_massima = if (length(presenti)) max(presenti) else NULL,
      moda = if (length(presenti)) moda(p$lunghezze) else NULL
    )
  }
  presenti <- which(totale_l > 0)
  mediane <- vapply(seq_len(nrow(totale_q)), function(i) {
    quantile_istogramma(totale_q[i, ], 0.5)
  }, integer(1))
  scrivi_json(
    list(
      campioni = length(profili),
      letture = sum(conteggi),
      lunghezza_minima = min(presenti),
      lunghezza_massima = max(presenti),
      moda = moda(totale_l),
      qualita_mediana_minima = min(mediane, na.rm = TRUE),
      qualita_mediana_per_posizione = as.list(mediane),
      per_campione = per_campione
    ),
    file.path(cartella, "riepilogo.json")
  )

  c("lunghezze.tsv", "qualita.tsv", tracciamento, "riepilogo.json")
})
