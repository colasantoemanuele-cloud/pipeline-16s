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
#   albero.nwk         con la filogenesi attiva: l'albero radicato, le foglie
#                      sono gli identificativi delle varianti
#
# Parametri: filtrato (ps_filtrato.rds di S13), albero (albero.nwk di S9, o
# nullo con la filogenesi disattivata), taxa_are_rows, export.
#
# L'ALBERO. Con la filogenesi attiva l'albero di S9 entra qui nell'oggetto
# finale: S10 assembla l'oggetto integrato sempre senza albero. Prima di unirlo
# si verifica che le foglie siano esattamente gli identificativi delle varianti
# dell'oggetto e che l'albero sia radicato, perche' phyloseq, davanti a foglie e
# varianti diverse, terrebbe in silenzio la sola intersezione. Le foglie si
# numerano nell'ordine delle varianti (R/lib/albero.R), cosi' che l'oggetto
# conservi quell'ordine.
#
# LA VALIDAZIONE (E-S14-01): componenti allineati sugli identificativi
# (conteggi, tassonomia, sequenze, metadati); le foglie dell'albero, se c'e',
# uguali agli identificativi delle varianti; solo campioni biologici; nessun
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
# S10 li rifiuta. L'albero si esporta in formato Newick con le lunghezze dei
# rami a cifre fisse (R/lib/albero.R): e' lo stesso testo scritto da S9.

for (f in c("io_json.R", "errors.R", "albero.R", "export.R")) {
  source(file.path(Sys.getenv("AMPLICON16S_R_LIB"), f))
}

# L'albero come componente confrontabile: il suo testo Newick, che ne fissa
# topologia, foglie e lunghezze; NULL se l'oggetto non ha l'albero.
componenti <- function(ps) {
  albero <- phyloseq::phy_tree(ps, errorIfNULL = FALSE)
  list(
    conteggi = methods::as(phyloseq::otu_table(ps), "matrix"),
    tassonomia = methods::as(phyloseq::tax_table(ps), "matrix"),
    metadati = methods::as(phyloseq::sample_data(ps), "data.frame"),
    sequenze = as.character(phyloseq::refseq(ps)),
    albero = if (is.null(albero)) NULL else testo_newick(albero)
  )
}

vuoto <- function(x) ifelse(is.na(x), "", x)

esegui_fase(function(parametri, cartella) {
  richiedi_pacchetti(c("phyloseq", "Biostrings"))
  ps <- readRDS(parametri$filtrato)
  difetti <- character()
  if (!is.null(phyloseq::phy_tree(ps, errorIfNULL = FALSE))) {
    difetti <- c(difetti, "l'oggetto filtrato di S13 ha gia' un albero")
  }
  if (!is.null(parametri$albero)) {
    richiedi_pacchetti("ape")
    albero <- leggi_newick(parametri$albero)
    difetti_foglie <- difetti_albero(albero, phyloseq::taxa_names(ps))
    if (length(difetti_foglie)) errore_catalogo("E-S14-01", paste(difetti_foglie, collapse = "; "))
    ordine <- phyloseq::taxa_names(ps)
    phyloseq::phy_tree(ps) <- foglie_in_ordine(albero, ordine)
    if (!identical(phyloseq::taxa_names(ps), ordine) ||
        !identical(phyloseq::phy_tree(ps)$tip.label, ordine)) {
      errore_catalogo("E-S14-01", "unire l'albero ha cambiato l'ordine o l'insieme delle varianti")
    }
  }
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
  for (nome in names(riletto)) {
    if (!identical(riletto[[nome]], componenti(ps)[[nome]])) {
      errore_catalogo("E-S14-01", sprintf("ps_final.rds riletto non coincide con l'oggetto: %s", nome))
    }
  }
  if (!identical(names(riletto), names(componenti(ps)))) {
    errore_catalogo("E-S14-01", "ps_final.rds riletto non ha gli stessi componenti dell'oggetto")
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
  if (!is.null(k$albero)) {
    scrivi_atomico(k$albero, file.path(cartella, "albero.nwk"))
    scritti <- c(scritti, "albero.nwk")
  }

  # ---- L'oggetto ricostruito dagli export -----------------------------------
  ricostruito <- ricostruisci(cartella)
  r <- componenti(ricostruito)
  if (!identical(names(r), names(k))) {
    errore_catalogo("E-S14-01", "gli export non ricostruiscono gli stessi componenti dell'oggetto")
  }
  for (nome in names(k)) {
    if (!identical(r[[nome]], k[[nome]])) {
      errore_catalogo("E-S14-01", sprintf("gli export non ricostruiscono l'oggetto: %s", nome))
    }
  }
  scritti
})
