# S14 - serializzazione, export e validazione dell'oggetto finale.
#
# Scrive in 12_final/, condivisa con S13:
#
#   ps_final.rds       l'oggetto finale (out.serialization = rds)
#   conteggi.tsv       con out.export_flat: varianti sulle righe, campioni sulle
#                      colonne, conteggi interi
#   tassonomia.tsv     varianti sulle righe, un rango per colonna, vuoto se non
#                      assegnato
#   metadati.tsv       campioni sulle righe, le colonne dei metadati
#   sequenze.fasta     una voce per variante: identificativo e sequenza
#
# Parametri: filtrato (ps_filtrato.rds di S13), taxa_are_rows, export.
#
# LA VALIDAZIONE (E-S14-01): componenti allineati sugli identificativi
# (conteggi, tassonomia, sequenze, metadati); solo campioni biologici; nessun
# campione e nessuna variante senza letture; l'oggetto riletto da ps_final.rds
# identico a quello scritto, componente per componente; con gli export,
# l'oggetto ricostruito dai soli file piatti identico a quello serializzato,
# componente per componente. I checksum e la
# frazione delle letture trattenute li verifica la fase Python.
#
# GLI EXPORT non dipendono dalle impostazioni locali: interi scritti con
# sprintf("%d"), nessun numero decimale, testo in UTF-8, separatore di
# tabulazione, a capo LF, nessun ordinamento (varianti nell'ordine degli
# identificativi, campioni nell'ordine dell'inventario), valori mancanti come
# campo vuoto. I valori dei metadati non contengono tabulazioni ne' a capo:
# S10 li rifiuta.

for (f in c("io_json.R", "errors.R", "export.R")) {
  source(file.path(Sys.getenv("AMPLICON16S_R_LIB"), f))
}

componenti <- function(ps) {
  list(
    conteggi = methods::as(phyloseq::otu_table(ps), "matrix"),
    tassonomia = methods::as(phyloseq::tax_table(ps), "matrix"),
    metadati = methods::as(phyloseq::sample_data(ps), "data.frame"),
    sequenze = as.character(phyloseq::refseq(ps))
  )
}

vuoto <- function(x) ifelse(is.na(x), "", x)

esegui_fase(function(parametri, cartella) {
  richiedi_pacchetti(c("phyloseq", "Biostrings"))
  ps <- readRDS(parametri$filtrato)
  difetti <- character()
  if (!identical(phyloseq::taxa_are_rows(ps), isTRUE(parametri$taxa_are_rows))) {
    difetti <- c(difetti, "orientamento diverso da out.taxa_are_rows")
  }
  k <- componenti(ps)
  if (!phyloseq::taxa_are_rows(ps)) k$conteggi <- t(k$conteggi)
  varianti <- rownames(k$conteggi)
  campioni <- colnames(k$conteggi)
  if (!identical(rownames(k$tassonomia), varianti)) difetti <- c(difetti, "tassonomia non allineata")
  if (!identical(names(k$sequenze), varianti)) difetti <- c(difetti, "sequenze non allineate")
  if (!identical(rownames(k$metadati), campioni)) difetti <- c(difetti, "metadati non allineati")
  if (any(k$metadati$classe != "biologico")) {
    difetti <- c(difetti, sprintf("%d campioni non biologici", sum(k$metadati$classe != "biologico")))
  }
  if (any(colSums(k$conteggi) == 0)) {
    difetti <- c(difetti, sprintf("%d campioni senza letture", sum(colSums(k$conteggi) == 0)))
  }
  if (any(rowSums(k$conteggi) == 0)) {
    difetti <- c(difetti, sprintf("%d varianti senza letture", sum(rowSums(k$conteggi) == 0)))
  }
  if (length(difetti)) errore_catalogo("E-S14-01", paste(difetti, collapse = "; "))

  # Il confronto e' per componenti: le sequenze (Biostrings) portano
  # riferimenti interni che identical() distingue anche fra copie uguali.
  finale <- file.path(cartella, "ps_final.rds")
  salva_rds(ps, finale)
  riletto <- componenti(readRDS(finale))
  for (nome in names(k)) {
    if (!identical(riletto[[nome]], componenti(ps)[[nome]])) {
      errore_catalogo("E-S14-01", sprintf("ps_final.rds riletto non coincide con l'oggetto: %s", nome))
    }
  }
  scritti <- "ps_final.rds"
  if (!isTRUE(parametri$export)) return(scritti)

  # ---- Gli export piatti ----------------------------------------------------
  m <- k$conteggi
  storage.mode(m) <- "integer"
  scrivi_atomico(
    c(paste(c("asv_id", campioni), collapse = "\t"),
      vapply(seq_along(varianti), function(i) {
        paste(c(varianti[i], sprintf("%d", m[i, ])), collapse = "\t")
      }, character(1))),
    file.path(cartella, "conteggi.tsv")
  )
  scrivi_atomico(
    c(paste(c("asv_id", colnames(k$tassonomia)), collapse = "\t"),
      vapply(seq_along(varianti), function(i) {
        paste(c(varianti[i], vuoto(k$tassonomia[i, ])), collapse = "\t")
      }, character(1))),
    file.path(cartella, "tassonomia.tsv")
  )
  scrivi_atomico(
    c(paste(colnames(k$metadati), collapse = "\t"),
      vapply(seq_along(campioni), function(i) {
        paste(vuoto(vapply(k$metadati[i, ], as.character, character(1))), collapse = "\t")
      }, character(1))),
    file.path(cartella, "metadati.tsv")
  )
  scrivi_atomico(c(rbind(paste0(">", varianti), unname(k$sequenze))),
                 file.path(cartella, "sequenze.fasta"))
  scritti <- c(scritti, "conteggi.tsv", "tassonomia.tsv", "metadati.tsv", "sequenze.fasta")

  # ---- L'oggetto ricostruito dagli export -----------------------------------
  ricostruito <- ricostruisci(cartella)
  r <- componenti(ricostruito)
  for (nome in names(k)) {
    if (!identical(r[[nome]], k[[nome]])) {
      errore_catalogo("E-S14-01", sprintf("gli export non ricostruiscono l'oggetto: %s", nome))
    }
  }
  scritti
})
