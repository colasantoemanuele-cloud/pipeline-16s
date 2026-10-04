# S10 - assemblaggio dell'oggetto integrato (phyloseq).
#
# Scrive in 10_phyloseq/:
#
#   ps_integrato.rds         l'oggetto phyloseq: conteggi, tassonomia, metadati
#                            dei campioni, sequenze di riferimento; senza albero
#                            (S9 disattivata)
#   varianti.tsv             per variante: identificativo, sequenza, lunghezza,
#                            letture nei biologici e in tutti i campioni
#   varianti_accessorie.tsv  per variante: il bootstrap di ciascun rango (S8) e i
#                            taxa col difetto noto del riferimento
#   riepilogo.json           dimensioni, campioni senza letture, regola degli
#                            identificativi e la variante ASV1
#
# Parametri: tabella (tabella_asv.rds di S7, campioni x varianti), tassonomia
# (tassonomia.rds di S8), difetti (difetto_riferimento.tsv di S8, o null),
# metadati (metadati_campioni.tsv, scritto dalla fase), colonne (i nomi attesi
# delle colonne dei metadati), campioni (gli accession dell'inventario, nel suo
# ordine), biologici (gli accession dei campioni biologici), taxa_are_rows,
# prefisso (degli identificativi delle varianti).
#
# Prima di caricare phyloseq si verificano gli ingressi e si costruiscono le
# matrici, con R di base; poi si assembla l'oggetto e lo si verifica di nuovo,
# slot per slot: phyloseq() tiene in silenzio solo i campioni e le varianti
# comuni a tutti i componenti, e un disallineamento non darebbe alcun errore.

for (f in c("io_json.R", "errors.R", "oggetto.R")) {
  source(file.path(Sys.getenv("AMPLICON16S_R_LIB"), f))
}

esegui_fase(function(parametri, cartella) {
  campioni <- as.character(unlist(parametri$campioni))
  biologici <- as.character(unlist(parametri$biologici))
  attese <- as.character(unlist(parametri$colonne))
  taxa_are_rows <- isTRUE(parametri$taxa_are_rows)

  # ---- Ingressi: la tabella di S7 ha i campioni sulle righe ---------------
  tabella <- readRDS(parametri$tabella)
  tassonomia <- readRDS(parametri$tassonomia)
  tax <- tassonomia$tax
  boot <- tassonomia$boot
  sequenze <- colnames(tabella)
  if (!is.null(rownames(tabella)) && identical(rownames(tabella), rownames(tax)) &&
      all(colnames(tabella) %in% campioni)) {
    errore_catalogo("E-S10-01", paste0(
      "tabella di S7 trasposta: le righe sono le varianti della tassonomia di S8 e ",
      "le colonne gli accession, mentre dada2 mette i campioni sulle righe"
    ))
  }
  verifica_nomi(rownames(tax), sequenze, "tassonomia di S8 - varianti della tabella di S7")
  verifica_nomi(rownames(boot), sequenze, "bootstrap di S8 - varianti della tabella di S7")
  presenti <- rownames(tabella)
  estranei <- setdiff(presenti, campioni)
  if (is.null(presenti) || anyDuplicated(presenti) || length(estranei) > 0L) {
    errore_catalogo("E-S10-01", sprintf(
      "le righe della tabella di S7 non sono accession distinti dell'inventario: %s",
      .descrivi_nomi(if (is.null(presenti)) NULL else c(estranei, presenti[duplicated(presenti)]))
    ))
  }
  if (anyNA(tabella) || any(tabella < 0) || any(tabella != round(tabella))) {
    stop("la tabella di S7 deve contenere conteggi interi non negativi")
  }

  # ---- Tutti i campioni dell'inventario, nell'ordine dell'inventario -------
  completa <- matrix(0L, nrow = length(campioni), ncol = length(sequenze),
                     dimnames = list(campioni, sequenze))
  completa[presenti, ] <- tabella
  storage.mode(completa) <- "integer"
  aggiunti <- setdiff(campioni, presenti)
  senza_letture <- campioni[rowSums(completa) == 0]

  # ---- Identificativi per abbondanza nei biologici -------------------------
  letture_biologici <- colSums(completa[campioni %in% biologici, , drop = FALSE])
  letture_totali <- colSums(completa)
  ordine <- ordine_varianti(letture_biologici, letture_totali, sequenze)
  ids <- paste0(parametri$prefisso, seq_along(ordine))
  ordinate <- sequenze[ordine]

  conteggi <- completa[, ordine, drop = FALSE]
  colnames(conteggi) <- ids
  otu <- if (taxa_are_rows) t(conteggi) else conteggi
  verifica_orientamento(otu, ids, campioni, taxa_are_rows, "tabella dei conteggi")

  tax <- tax[ordinate, , drop = FALSE]
  rownames(tax) <- ids
  boot <- boot[ordinate, , drop = FALSE]
  rownames(boot) <- ids

  # ---- Metadati: i nomi scelti dalla fase, senza rinomine ------------------
  verifica_sintattici(attese, "colonne dei metadati")
  metadati <- utils::read.delim(
    parametri$metadati, colClasses = "character", quote = "", comment.char = "",
    check.names = FALSE, na.strings = "", encoding = "UTF-8"
  )
  verifica_nomi(colnames(metadati), attese, "metadati - colonne")
  verifica_nomi(metadati$accession, campioni, "metadati - campioni")
  rownames(metadati) <- metadati$accession

  # ---- L'oggetto -------------------------------------------------------------
  richiedi_pacchetti(c("phyloseq", "Biostrings"))
  riferimento <- Biostrings::DNAStringSet(stats::setNames(ordinate, ids))
  ps <- phyloseq::phyloseq(
    phyloseq::otu_table(otu, taxa_are_rows = taxa_are_rows),
    phyloseq::tax_table(tax),
    phyloseq::sample_data(metadati),
    riferimento
  )

  # Ogni slot, di nuovo, sugli identificativi.
  if (!identical(phyloseq::taxa_are_rows(ps), taxa_are_rows)) {
    errore_catalogo("E-S10-01", "l'oggetto non ha l'orientamento di out.taxa_are_rows")
  }
  verifica_orientamento(methods::as(phyloseq::otu_table(ps), "matrix"), ids, campioni,
                        taxa_are_rows, "oggetto - tabella dei conteggi")
  verifica_nomi(phyloseq::taxa_names(ps), ids, "oggetto - varianti")
  verifica_nomi(phyloseq::sample_names(ps), campioni, "oggetto - campioni")
  verifica_nomi(rownames(phyloseq::tax_table(ps)), ids, "oggetto - tassonomia")
  verifica_nomi(colnames(phyloseq::tax_table(ps)), colnames(tax), "oggetto - ranghi")
  verifica_nomi(names(phyloseq::refseq(ps)), ids, "oggetto - sequenze di riferimento")
  verifica_nomi(as.character(phyloseq::refseq(ps)), ordinate,
                "oggetto - sequenze di riferimento")
  dati <- methods::as(phyloseq::sample_data(ps), "data.frame")
  verifica_nomi(colnames(dati), attese, "oggetto - colonne dei metadati")
  verifica_nomi(rownames(dati), campioni, "oggetto - righe dei metadati")
  for (colonna in attese) {
    if (!identical(dati[[colonna]], metadati[[colonna]])) {
      errore_catalogo("E-S10-01", sprintf("oggetto - metadati: valori diversi in %s", colonna))
    }
  }
  if (!is.null(phyloseq::phy_tree(ps, errorIfNULL = FALSE))) {
    stop("l'oggetto non deve avere un albero: S9 non e' attiva")
  }
  if (sum(phyloseq::otu_table(ps)) != sum(tabella)) {
    errore_catalogo("E-S10-01", sprintf(
      "oggetto - letture: %.0f nell'oggetto, %.0f nella tabella di S7",
      sum(phyloseq::otu_table(ps)), sum(tabella)
    ))
  }
  salva_rds(ps, file.path(cartella, "ps_integrato.rds"))

  # ---- Le varianti e le informazioni accessorie, per identificativo --------
  scrivi_atomico(
    c("asv_id\tsequenza\tlunghezza\tletture_biologici\tletture_totali",
      sprintf("%s\t%s\t%d\t%.0f\t%.0f", ids, ordinate, nchar(ordinate),
              letture_biologici[ordine], letture_totali[ordine])),
    file.path(cartella, "varianti.tsv")
  )
  marcate <- stats::setNames(character(length(ids)), ordinate)
  if (!is.null(parametri$difetti)) {
    difetti <- utils::read.delim(parametri$difetti, colClasses = "character",
                                 quote = "", comment.char = "", check.names = FALSE)
    estranee <- setdiff(difetti$sequenza, ordinate)
    if (length(estranee) > 0L) {
      errore_catalogo("E-S10-01", sprintf(
        "difetto_riferimento.tsv di S8: %d sequenze non sono varianti della tabella",
        length(estranee)
      ))
    }
    voci <- sprintf("%s:%s(%s)", difetti$colonna, difetti$taxon, difetti$rango)
    for (s in unique(difetti$sequenza)) {
      marcate[[s]] <- paste(voci[difetti$sequenza == s], collapse = ";")
    }
  }
  scrivi_atomico(
    c(paste(c("asv_id", paste0("boot_", colnames(boot)), "difetto_riferimento"),
            collapse = "\t"),
      sprintf("%s\t%s\t%s", ids,
              apply(boot, 1, function(r) paste(sprintf("%d", as.integer(r)), collapse = "\t")),
              unname(marcate[ordinate]))),
    file.path(cartella, "varianti_accessorie.tsv")
  )

  prima <- ordine[1L]
  scrivi_json(
    list(
      componenti = list("otu_table", "tax_table", "sample_data", "refseq"),
      albero = NULL,
      taxa_are_rows = taxa_are_rows,
      campioni = length(campioni),
      varianti = length(ids),
      letture = sum(as.numeric(completa)),
      campioni_senza_letture = as.list(senza_letture),
      campioni_aggiunti_a_zero = as.list(aggiunti),
      identificativi = list(
        schema = "abundance_rank",
        ordine = paste0(
          "letture nei campioni biologici, decrescenti; a parita', letture in tutti ",
          "i campioni, decrescenti; a parita', la sequenza in ordine lessicografico"
        ),
        campioni_biologici = sum(campioni %in% biologici)
      ),
      colonne_metadati = as.list(attese),
      difetto_riferimento_valutato = !is.null(parametri$difetti),
      prima_variante = list(
        asv_id = ids[1L],
        sequenza = sequenze[prima],
        letture_biologici = letture_biologici[[prima]],
        letture_totali = letture_totali[[prima]],
        tassonomia = as.list(tax[1L, ])
      )
    ),
    file.path(cartella, "riepilogo.json")
  )

  c("ps_integrato.rds", "varianti.tsv", "varianti_accessorie.tsv", "riepilogo.json")
})
