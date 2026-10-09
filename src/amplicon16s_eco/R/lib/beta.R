# Beta diversita': matrici di distanza fra campioni.
#
# Richiede uscite.R e campioni.R. Interfaccia:
#   clr(conteggi, pseudoconteggio)   trasformata log-rapporto centrata
#   distanza(metodo, ps, conteggi, pseudoconteggio)   matrice quadrata
#   scrivi_distanze(...)             una tabella per distanza
#
# Razionale biologico: ogni distanza e' calcolata sui conteggi non rarefatti,
# cosi' non dipende dal seme. Jaccard e' in presenza e assenza: il predefinito
# di vegan::vegdist (binary = FALSE) calcola la forma quantitativa, che e'
# un'altra misura. Bray-Curtis sulle proporzioni toglie la differenza di
# profondita' fra campioni; Aitchison lavora sui log-rapporti, e richiede uno
# pseudoconteggio dichiarato perche' il logaritmo di zero non esiste.

DESCRIZIONE_DISTANZA <- c(
  bray = "Bray-Curtis sulle abbondanze relative dei conteggi non rarefatti (vegan::vegdist, method = \"bray\")",
  jaccard = "Jaccard in presenza e assenza (vegan::vegdist, method = \"jaccard\", binary = TRUE)",
  aitchison = "Aitchison: distanza euclidea sul CLR di conteggi + pseudoconteggio",
  unifrac_weighted = "UniFrac pesato normalizzato (phyloseq::UniFrac, weighted = TRUE, normalized = TRUE)",
  unifrac_unweighted = "UniFrac non pesato (phyloseq::UniFrac, weighted = FALSE)"
)

clr <- function(conteggi, pseudoconteggio) {
  # Una riga per campione: logaritmo meno la media dei logaritmi del campione.
  logaritmi <- log(t(conteggi) + pseudoconteggio)
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

scrivi_distanze <- function(matrici, pseudoconteggio, cartella) {
  dir.create(cartella, showWarnings = FALSE)
  prodotti <- character()
  for (metodo in names(matrici)) {
    m <- matrici[[metodo]]
    tabella <- data.frame(campione = rownames(m), stringsAsFactors = FALSE,
                          check.names = FALSE)
    for (campione in colnames(m)) tabella[[campione]] <- reale(m[, campione])
    scrivi_tsv(
      tabella, file.path(cartella, paste0("distanza_", metodo, ".tsv")),
      intestazione = c(
        DESCRIZIONE_DISTANZA[[metodo]],
        if (metodo == "aitchison") sprintf("pseudoconteggio: %s", format(pseudoconteggio, digits = 15)),
        "matrice quadrata simmetrica, campioni in ordine alfabetico (radix)"
      )
    )
    prodotti <- c(prodotti, paste0("beta/distanza_", metodo, ".tsv"))
  }
  prodotti
}
