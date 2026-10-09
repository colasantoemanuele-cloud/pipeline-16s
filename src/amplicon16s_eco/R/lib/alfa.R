# Alfa diversita' su conteggi rarefatti.
#
# Richiede uscite.R. Interfaccia:
#   fissa_seme(seme)                 generatore e seme, prima di ogni operazione
#                                    casuale
#   rarefai(ps, campioni, profondita, seme)   matrice varianti x campioni
#   indici_alfa(rarefatti)           ricchezza, Shannon, Gini-Simpson
#   scrivi_alfa(...)                 le due tabelle
#
# Razionale biologico: la rarefazione e' senza reinserimento, cioe' un
# sottocampione delle letture osservate: con il reinserimento (il predefinito
# di phyloseq::rarefy_even_depth) una variante puo' ricevere piu' letture di
# quante ne abbia, e una variante rara puo' essere estratta piu' volte gonfiando
# gli indici di equiripartizione. Nessuna variante viene tolta dopo la
# rarefazione (trimOTUs = FALSE): la tabella rarefatta ha le stesse righe di
# quella di partenza.

fissa_seme <- function(seme) {
  # Il generatore e' dichiarato, non ereditato: il predefinito di R e' cambiato
  # fra versioni (sample.kind "Rounding" prima della 3.6).
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  set.seed(as.integer(seme))
}

rarefai <- function(ps, campioni, profondita, seme) {
  sotto <- phyloseq::prune_samples(campioni, ps)
  fissa_seme(seme)
  rarefatto <- suppressMessages(phyloseq::rarefy_even_depth(
    sotto, sample.size = profondita, rngseed = as.integer(seme),
    replace = FALSE, trimOTUs = FALSE, verbose = FALSE
  ))
  m <- conteggi_di(rarefatto)
  m[, campioni, drop = FALSE]
}

indici_alfa <- function(rarefatti) {
  per_campione <- t(rarefatti)
  data.frame(
    campione = colnames(rarefatti),
    ricchezza_osservata = as.integer(colSums(rarefatti > 0)),
    shannon = as.numeric(vegan::diversity(per_campione, index = "shannon",
                                         MARGIN = 1, base = exp(1))),
    gini_simpson = as.numeric(vegan::diversity(per_campione, index = "simpson",
                                              MARGIN = 1)),
    stringsAsFactors = FALSE
  )
}

scrivi_alfa <- function(indici, rarefatti, letture, profondita, origine, seme, cartella) {
  dir.create(cartella, showWarnings = FALSE)
  scrivi_tsv(
    data.frame(
      campione = indici$campione,
      letture = as.integer(letture[indici$campione]),
      ricchezza_osservata = indici$ricchezza_osservata,
      shannon = reale(indici$shannon),
      gini_simpson = reale(indici$gini_simpson),
      stringsAsFactors = FALSE
    ),
    file.path(cartella, "alfa_diversita.tsv"),
    intestazione = c(
      "Alfa diversita' per campione, su conteggi rarefatti.",
      sprintf("rarefazione: senza reinserimento, profondita' %d (%s), seme %d",
              profondita, origine, as.integer(seme)),
      "letture: letture del campione prima della rarefazione",
      "ricchezza_osservata: numero di varianti con almeno una lettura nel campione rarefatto",
      "shannon: -sum(p_i * ln(p_i)), logaritmo naturale, p_i proporzione della variante i",
      "gini_simpson: 1 - sum(p_i^2) (indice di Simpson nella forma di Gini-Simpson)"
    )
  )
  scrivi_tsv(
    affianca(list(variante = rownames(rarefatti)), rarefatti, as.integer), file.path(cartella, "conteggi_rarefatti.tsv"),
    intestazione = c(
      "Conteggi rarefatti: una riga per variante, una colonna per campione.",
      sprintf("rarefazione: senza reinserimento, profondita' %d (%s), seme %d",
              profondita, origine, as.integer(seme))
    )
  )
  c("alfa/alfa_diversita.tsv", "alfa/conteggi_rarefatti.tsv")
}
