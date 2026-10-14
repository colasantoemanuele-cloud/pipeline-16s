# La decontaminazione per prevalenza (S12), per piastra e aggregata.
#
# Richiede il pacchetto decontam. Interfaccia:
#   probabilita_prevalenza(conteggi, negativo)
#       la probabilita' di decontam::isContaminant (method = "prevalence") per
#       ogni variante, con i nomi delle varianti; NA dove decontam non la calcola
#   combina_probabilita(p, regola)
#       le probabilita' di piu' confronti combinate come batch.combine di
#       decontam ("minimum", "product", "fisher"), variante per variante
#
# Il confronto e' fra i controlli negativi (negativo vero) e i campioni
# (negativo falso). Un campione senza letture non entra: decontam lo
# toglierebbe da se', con un avviso.
#
# La combinazione riproduce quella di decontam 1.28.0 (.is_contaminant):
# "minimum" e' il minimo delle probabilita' non mancanti, NA se mancano tutte;
# "product" il loro prodotto; "fisher" il metodo di Fisher con le mancanti
# sostituite da 0,5. Si riproduce invece di passare batch a decontam perche'
# una piastra con meno di decontam.min_blanks negativi ha un confronto diverso
# dalle altre (i suoi campioni contro i negativi di tutte le piastre), che il
# parametro batch non sa esprimere; che la combinazione coincida con quella di
# decontam quando ogni piastra ha il suo confronto e' verificato nei test.

# La probabilita' di decontam per ogni variante, nel confronto fra i controlli
# negativi e gli altri campioni. conteggi e' la matrice campioni x varianti,
# negativo un vettore logico per campione. I campioni senza letture si tolgono
# prima della chiamata. La soglia passata a decontam non conta: di decontam si
# usa la sola probabilita' (detailed = TRUE), e la decisione sul contaminante
# la prende chi chiama, con la soglia della configurazione. Restituisce un
# vettore con i nomi delle varianti, NA dove decontam non calcola.
probabilita_prevalenza <- function(conteggi, negativo) {
  presenti <- rowSums(conteggi) > 0
  esito <- decontam::isContaminant(
    conteggi[presenti, , drop = FALSE], neg = negativo[presenti],
    method = "prevalence", threshold = 0.5, normalize = TRUE, detailed = TRUE
  )
  stats::setNames(esito$p, colnames(conteggi))
}

# Le probabilita' di piu' confronti combinate in una per variante. p e' una
# matrice varianti x confronti (un vettore vale come un solo confronto), regola
# e' minimum, product o fisher, con il trattamento dei valori mancanti
# descritto in testa al file. Restituisce un vettore con i nomi delle varianti;
# una regola sconosciuta ferma il calcolo.
combina_probabilita <- function(p, regola) {
  if (is.null(dim(p))) p <- matrix(p, ncol = 1L, dimnames = list(names(p), NULL))
  combinate <- switch(
    regola,
    minimum = apply(p, 1, function(v) {
      v <- v[!is.na(v)]
      if (length(v)) min(v) else NA_real_
    }),
    product = apply(p, 1, function(v) prod(v, na.rm = TRUE)),
    fisher = apply(p, 1, function(v) {
      v[is.na(v)] <- 0.5
      stats::pchisq(-2 * log(prod(v)), df = 2 * length(v), lower.tail = FALSE)
    }),
    stop("regola di combinazione sconosciuta: ", regola)
  )
  stats::setNames(combinate, rownames(p))
}
