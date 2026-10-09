# Beta diversita': matrici di distanza fra campioni.
#
# Richiede uscite.R e campioni.R. Interfaccia:
#   clr(conteggi, pseudoconteggio)   trasformata log-rapporto centrata
#   radica(ps, modo)                 l'oggetto con l'albero radicato per UniFrac
#   distanza(metodo, ps, conteggi, pseudoconteggio)   matrice quadrata
#   scrivi_distanze(...)             una tabella per distanza
#
# Razionale biologico: ogni distanza e' calcolata sui conteggi non rarefatti,
# cosi' non dipende dal seme. Jaccard e' in presenza e assenza: il predefinito
# di vegan::vegdist (binary = FALSE) calcola la forma quantitativa, che e'
# un'altra misura. Bray-Curtis sulle proporzioni toglie la differenza di
# profondita' fra campioni; Aitchison lavora sui log-rapporti, e richiede uno
# pseudoconteggio dichiarato perche' il logaritmo di zero non esiste. UniFrac
# misura i rami che portano dalla radice ai taxa di un campione, quindi dipende
# dalla radice: un albero stimato per massima verosimiglianza con un modello
# reversibile non ne ha una, e il modo di ottenerla e' dichiarato.

DESCRIZIONE_DISTANZA <- c(
  bray = "Bray-Curtis sulle abbondanze relative dei conteggi non rarefatti (vegan::vegdist, method = \"bray\")",
  jaccard = "Jaccard in presenza e assenza (vegan::vegdist, method = \"jaccard\", binary = TRUE)",
  aitchison = "Aitchison: distanza euclidea sul CLR di conteggi + pseudoconteggio",
  unifrac_weighted = "UniFrac pesato normalizzato (phyloseq::UniFrac, weighted = TRUE, normalized = TRUE)",
  unifrac_unweighted = "UniFrac non pesato (phyloseq::UniFrac, weighted = FALSE)"
)

DESCRIZIONE_RADICE <- c(
  midpoint = "albero radicato al punto medio del cammino piu' lungo fra due foglie (phangorn::midpoint)",
  existing = "radice dell'albero dell'oggetto"
)

# L'oggetto con l'albero da usare per UniFrac. La radicazione al punto medio
# non usa numeri casuali e si applica anche a un albero gia' radicato: il modo
# dichiarato e' quello eseguito. Cambia il solo oggetto in memoria, letto dalla
# copia: quello di partenza non viene mai scritto.
radica <- function(ps, modo) {
  if (modo == "midpoint") {
    phyloseq::phy_tree(ps) <- phangorn::midpoint(phyloseq::phy_tree(ps))
  }
  ps
}

# Le varianti che entrano nel CLR: quelle con almeno una lettura fra i campioni
# analizzati. Una variante presente solo in campioni esclusi non deve cambiare
# la distanza fra quelli analizzati.
varianti_presenti <- function(conteggi) rowSums(conteggi) > 0

clr <- function(conteggi, pseudoconteggio) {
  # Una riga per campione: logaritmo meno la media dei logaritmi del campione.
  logaritmi <- log(t(conteggi[varianti_presenti(conteggi), , drop = FALSE]) + pseudoconteggio)
  logaritmi - rowMeans(logaritmi)
}

distanza <- function(metodo, ps, conteggi, pseudoconteggio = NULL) {
  campioni <- colnames(conteggi)
  d <- switch(
    metodo,
    bray = vegan::vegdist(t(abbondanze_relative(conteggi)), method = "bray", binary = FALSE),
    jaccard = vegan::vegdist(t(conteggi), method = "jaccard", binary = TRUE),
    aitchison = stats::dist(clr(conteggi, pseudoconteggio), method = "euclidean"),
    unifrac_weighted = phyloseq::UniFrac(
      phyloseq::prune_samples(campioni, ps), weighted = TRUE, normalized = TRUE,
      parallel = FALSE, fast = TRUE),
    unifrac_unweighted = phyloseq::UniFrac(
      phyloseq::prune_samples(campioni, ps), weighted = FALSE, parallel = FALSE,
      fast = TRUE),
    stop("distanza non realizzata: ", metodo)
  )
  # Una matrice senza attributi del metodo: chi la usa non ne eredita scelte.
  m <- as.matrix(d)
  m <- m[campioni, campioni, drop = FALSE]
  dimnames(m) <- list(campioni, campioni)
  m
}

scrivi_distanze <- function(matrici, pseudoconteggio, varianti_nel_clr, cartella,
                            radice = NULL, radicato = NA) {
  dir.create(cartella, showWarnings = FALSE)
  prodotti <- character()
  for (metodo in names(matrici)) {
    m <- matrici[[metodo]]
    scrivi_tsv(
      affianca(list(campione = rownames(m)), m, reale), file.path(cartella, paste0("distanza_", metodo, ".tsv")),
      intestazione = c(
        DESCRIZIONE_DISTANZA[[metodo]],
        if (metodo == "aitchison") sprintf(
          "pseudoconteggio: %s; varianti nel CLR (almeno una lettura fra i campioni analizzati): %d",
          format(pseudoconteggio, digits = 15), varianti_nel_clr),
        if (grepl("^unifrac", metodo)) sprintf(
          "radice (beta.unifrac_root = %s): %s; l'albero dell'oggetto %s",
          radice, DESCRIZIONE_RADICE[[radice]],
          if (radicato) "aveva gia' una radice" else "non aveva una radice"),
        "matrice quadrata simmetrica, campioni in ordine alfabetico (radix)"
      )
    )
    prodotti <- c(prodotti, paste0("beta/distanza_", metodo, ".tsv"))
  }
  prodotti
}
