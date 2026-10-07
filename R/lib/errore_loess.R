# Le funzioni che stimano i tassi d'errore dalle transizioni osservate
# (err.error_function), usate da S3 (dada2::learnErrors) e dalla seconda
# passata di S4.
#
#   loess            dada2::loessErrfun, la stima standard di dada2;
#   loess_monotono   la variante per le qualita' raggruppate in pochi valori
#                    (NovaSeq, NextSeq).
#
# PERCHE' UNA VARIANTE MONOTONA. Con le qualita' raggruppate ogni valore
# riunisce basi di affidabilita' diversa, i punti su cui il loess si adatta
# sono pochi, e la curva stimata puo' salire con la qualita': un tasso d'errore
# piu' alto a qualita' piu' alta, che non ha senso fisico e fa scambiare a dada
# varianti vere per errori (o il contrario).
#
# CHE COSA FA loess_monotono, ESATTAMENTE. Per ognuna delle dodici sostituzioni
# (A>C, A>G, ...):
#   1. calcola a ogni qualita' il logaritmo in base 10 del tasso osservato,
#      (errori + 1) / totale, come dada2::loessErrfun; una qualita' mai
#      osservata, o con una sola base (peso zero), non entra nell'adattamento;
#   2. adatta a quei punti un loess di PRIMO grado (degree = 1), con span = 2
#      e peso log10(totale) per ogni punto. Con span maggiore di 1 ogni
#      adattamento locale usa tutti i punti, con pesi che decrescono piano con
#      la distanza: la curva e' quasi una retta pesata, e resta definita anche
#      con tre o quattro valori di qualita', dove lo span 0,75 del loess
#      standard (secondo grado, pesi pari ai totali) interpola i punti o si
#      ferma. Lo span 2 e i pesi logaritmici sono le modifiche discusse dalla
#      comunita' di dada2 per i dati NovaSeq (benjjneb/dada2, issue 1307, che
#      riprende la 938); il primo grado limita le oscillazioni fra pochi punti;
#   3. con due soli valori di qualita' osservati un loess non e' adattabile: si
#      usa la retta per i due punti (stats::lm). Con uno solo non c'e' una
#      curva, e la funzione si ferma: S3 lo dichiara con E-S3-04;
#   4. fuori dalle qualita' osservate vale la stima del bordo; i tassi sono
#      limitati fra 1e-7 e 0,25, come in dada2;
#   5. impone la monotonia: a ogni qualita' il tasso e' almeno quello di tutte
#      le qualita' superiori (massimo cumulato da destra). E' un vincolo piu'
#      forte di quello della issue, che alza i tassi al solo valore della
#      qualita' massima.
# E' scritta con il solo R di base, senza dipendenze in piu'.

# I limiti dei tassi stimati, gli stessi di dada2::loessErrfun.
TASSO_MASSIMO <- 0.25
TASSO_MINIMO <- 1e-7
# L'ampiezza del loess della variante monotona, e il numero di qualita'
# osservate sotto il quale un loess di primo grado non e' adattabile.
SPAN_MONOTONO <- 2
GRADI_MINIMI_LOESS <- 3L

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
      osservati <- which(!is.na(logp) & dati$peso > 0)
      if (length(osservati) < 2L) {
        # Un solo valore di qualita' (o nessuno): non c'e' una curva da
        # adattare. learnErrors lo riporta come matrice d'errore nulla, e S3
        # lo dichiara con il suo codice.
        stop("loess_monotono: ", length(osservati), " valori di qualita' osservati per ",
             da, "2", a, ", ne servono almeno due")
      }
      punti <- dati[osservati, ]
      if (length(osservati) >= GRADI_MINIMI_LOESS) {
        modello <- suppressWarnings(
          stats::loess(logp ~ q, punti, weights = punti$peso, degree = 1, span = SPAN_MONOTONO)
        )
      } else {
        modello <- stats::lm(logp ~ q, punti, weights = punti$peso)
      }
      previsto <- as.numeric(stats::predict(modello, data.frame(q = qualita)))
      previsto[qualita < min(punti$q) | qualita > max(punti$q)] <- NA
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
