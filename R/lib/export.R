# La lettura degli export piatti dell'oggetto finale (S14).
#
# Interfaccia:
#   ricostruisci(cartella)  l'oggetto phyloseq ricostruito da conteggi.tsv,
#                           tassonomia.tsv, metadati.tsv e sequenze.fasta, con
#                           l'albero di albero.nwk se il file c'e' (richiede
#                           allora albero.R)
#
# Serve a S14 per verificare che gli export bastino a ricostruire l'oggetto
# serializzato. I file si leggono come testo, con i tipi dichiarati: interi i
# conteggi, testo tutto il resto, il campo vuoto come valore mancante, nessun
# nome di colonna alterato (check.names = FALSE).

# Una tabella degli export letta come testo con i tipi dichiarati in classi:
# nessuna interpretazione delle virgolette ne' dei commenti, campo vuoto come
# valore mancante, nomi di colonna lasciati come sono.
leggi_tabella <- function(percorso, classi) {
  utils::read.delim(percorso, colClasses = classi, quote = "", comment.char = "",
                    check.names = FALSE, na.strings = "", encoding = "UTF-8")
}

# L'oggetto phyloseq ricostruito dai soli export piatti di una cartella. Passi:
# i conteggi, letti come interi dopo aver contato le colonne dall'intestazione;
# la tassonomia e i metadati, letti come testo; le sequenze dal FASTA; l'unione
# dei quattro componenti; infine l'albero, se il file Newick c'e', con le
# foglie numerate nell'ordine delle varianti come fa la fase che scrive
# l'oggetto. Il confronto con l'oggetto serializzato lo fa chi chiama.
ricostruisci <- function(cartella) {
  intestazione <- strsplit(readLines(file.path(cartella, "conteggi.tsv"), n = 1L,
                                     encoding = "UTF-8"), "\t", fixed = TRUE)[[1]]
  conteggi <- leggi_tabella(file.path(cartella, "conteggi.tsv"),
                            c("character", rep("integer", length(intestazione) - 1L)))
  m <- as.matrix(conteggi[, -1L, drop = FALSE])
  rownames(m) <- conteggi$asv_id
  tassonomia <- leggi_tabella(file.path(cartella, "tassonomia.tsv"), "character")
  t <- as.matrix(tassonomia[, -1L, drop = FALSE])
  rownames(t) <- tassonomia$asv_id
  metadati <- leggi_tabella(file.path(cartella, "metadati.tsv"), "character")
  rownames(metadati) <- metadati$accession
  sequenze <- Biostrings::readDNAStringSet(file.path(cartella, "sequenze.fasta"))
  ps <- phyloseq::phyloseq(
    phyloseq::otu_table(m, taxa_are_rows = TRUE),
    phyloseq::tax_table(t),
    phyloseq::sample_data(metadati),
    sequenze
  )
  albero <- file.path(cartella, "albero.nwk")
  if (file.exists(albero)) {
    phyloseq::phy_tree(ps) <- foglie_in_ordine(leggi_newick(albero), phyloseq::taxa_names(ps))
  }
  ps
}
