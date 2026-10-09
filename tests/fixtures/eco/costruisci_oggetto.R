# Costruisce un oggetto phyloseq di prova da una descrizione JSON.
#
# Uso: Rscript --vanilla costruisci_oggetto.R descrizione.json oggetto.rds
#
# La descrizione ha: conteggi (una riga per variante, una colonna per
# campione), varianti e campioni (i nomi), tassonomia (rango -> un valore per
# variante, null se non assegnata), metadati (colonna -> un valore per
# campione, null se mancante) e, facoltativo, albero (testo Newick).

argomenti <- commandArgs(trailingOnly = TRUE)
d <- jsonlite::read_json(argomenti[[1L]], simplifyVector = FALSE)

colonna <- function(x) vapply(x, function(v) if (is.null(v)) NA_character_ else as.character(v), "")

conteggi <- do.call(rbind, lapply(d$conteggi, function(riga) as.integer(unlist(riga))))
dimnames(conteggi) <- list(unlist(d$varianti), unlist(d$campioni))
tassonomia <- do.call(cbind, lapply(d$tassonomia, colonna))
dimnames(tassonomia) <- list(unlist(d$varianti), names(d$tassonomia))
metadati <- as.data.frame(lapply(d$metadati, colonna), stringsAsFactors = FALSE,
                          check.names = FALSE)
names(metadati) <- names(d$metadati)
rownames(metadati) <- unlist(d$campioni)

componenti <- list(
  phyloseq::otu_table(conteggi, taxa_are_rows = TRUE),
  phyloseq::tax_table(tassonomia),
  phyloseq::sample_data(metadati)
)
if (!is.null(d$albero)) {
  componenti[[length(componenti) + 1L]] <- phyloseq::phy_tree(ape::read.tree(text = d$albero))
}
saveRDS(do.call(phyloseq::phyloseq, componenti), argomenti[[2L]])
