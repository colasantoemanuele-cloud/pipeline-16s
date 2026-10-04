# La lettura degli export piatti dell'oggetto finale (S14).
#
# Interfaccia:
#   ricostruisci(cartella)  l'oggetto phyloseq ricostruito da conteggi.tsv,
#                           tassonomia.tsv, metadati.tsv e sequenze.fasta
#
# Serve a S14 per verificare che gli export bastino a ricostruire l'oggetto
# serializzato. I file si leggono come testo, con i tipi dichiarati: interi i
# conteggi, testo tutto il resto, il campo vuoto come valore mancante, nessun
# nome di colonna alterato (check.names = FALSE).

leggi_tabella <- function(percorso, classi) {
  utils::read.delim(percorso, colClasses = classi, quote = "", comment.char = "",
                    check.names = FALSE, na.strings = "", encoding = "UTF-8")
}

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
  phyloseq::phyloseq(
    phyloseq::otu_table(m, taxa_are_rows = TRUE),
    phyloseq::tax_table(t),
    phyloseq::sample_data(metadati),
    sequenze
  )
}
