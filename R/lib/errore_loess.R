# Le funzioni che stimano i tassi d'errore dalle transizioni osservate
# (err.error_function), usate da S3 (dada2::learnErrors) e dalla seconda
# passata di S4.
#
#   loess            dada2::loessErrfun, la stima standard di dada2;
#   loess_monotono   la stessa stima con due modifiche, per le qualita'
#                    raggruppate in pochi valori (NovaSeq, NextSeq).
#
# PERCHE' UNA VARIANTE MONOTONA. Con le qualita' raggruppate ogni valore
# riunisce basi di affidabilita' diversa, i punti su cui il loess si adatta
# sono pochi, e la curva stimata puo' salire con la qualita': un tasso d'errore
# piu' alto a qualita' piu' alta, che non ha senso fisico e fa scambiare a dada
# varianti vere per errori (o il contrario). La variante:
#   1. pesa ogni punto con il logaritmo delle basi osservate e usa un loess di
#      primo grado, meno incline a oscillare fra pochi punti;
#   2. impone che il tasso d'errore non cresca con la qualita': a ogni qualita'
#      il tasso e' almeno quello di tutte le qualita' superiori.
# E' la correzione in uso nella comunita' di dada2 per questi dati; qui e'
# scritta con il solo R di base, senza dipendenze in piu'.

# I limiti dei tassi stimati, gli stessi di dada2::loessErrfun.
TASSO_MASSIMO <- 0.25
TASSO_MINIMO <- 1e-7

loess_monotono <- function(trans) {
  qualita <- as.numeric(colnames(trans))
  basi <- c("A", "C", "G", "T")
  stime <- matrix(0, nrow = 0, ncol = length(qualita))
  for (da in basi) {
    for (a in basi) {
      if (da == a) next
      errori <- trans[paste0(da, "2", a), ]
      totale <- colSums(trans[paste0(da, "2", basi), ])
      # Lo stesso pseudoconteggio di dada2; una qualita' mai osservata non
      # entra nell'adattamento.
      logp <- log10((errori + 1) / totale)
      logp[is.infinite(logp)] <- NA
      dati <- data.frame(q = qualita, logp = logp, peso = log10(pmax(totale, 1)))
      modello <- suppressWarnings(
        stats::loess(logp ~ q, dati, weights = dati$peso, degree = 1, span = 0.95)
      )
      previsto <- stats::predict(modello, qualita)
      # Fuori dalle qualita' osservate vale il valore del bordo.
      validi <- which(!is.na(previsto))
      previsto[seq_along(previsto) > max(validi)] <- previsto[[max(validi)]]
      previsto[seq_along(previsto) < min(validi)] <- previsto[[min(validi)]]
      stime <- rbind(stime, 10^previsto)
    }
  }
  stime[stime > TASSO_MASSIMO] <- TASSO_MASSIMO
  stime[stime < TASSO_MINIMO] <- TASSO_MINIMO
  # Monotonia: il tasso a una qualita' e' il massimo fra il suo e quelli delle
  # qualita' superiori (le colonne sono in ordine di qualita' crescente).
  stime <- t(apply(stime, 1, function(riga) rev(cummax(rev(riga)))))

  # La matrice 16 x qualita' che dada si aspetta: per ogni base di partenza le
  # tre sostituzioni stimate e, per differenza, la probabilita' di non errore.
  err <- rbind(
    1 - colSums(stime[1:3, , drop = FALSE]), stime[1:3, , drop = FALSE],
    stime[4, ], 1 - colSums(stime[4:6, , drop = FALSE]), stime[5:6, , drop = FALSE],
    stime[7:8, , drop = FALSE], 1 - colSums(stime[7:9, , drop = FALSE]), stime[9, ],
    stime[10:12, , drop = FALSE], 1 - colSums(stime[10:12, , drop = FALSE])
  )
  rownames(err) <- paste0(rep(basi, each = 4), "2", basi)
  colnames(err) <- colnames(trans)
  err
}

# La funzione corrispondente a err.error_function.
funzione_errore <- function(nome) {
  switch(
    nome,
    loess = dada2::loessErrfun,
    loess_monotono = loess_monotono,
    stop("err.error_function non realizzata: ", nome)
  )
}
