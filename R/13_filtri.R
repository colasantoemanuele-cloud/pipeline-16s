# S13 - filtri finali: tassonomico e di prevalenza, dopo quello per profondita'.
#
# Scrive in 12_final/intermedi/: sono gli intermedi dei filtri, non i file
# consegnati, che S14 scrive in 12_final/.
#
#   ps_filtrato.rds          l'oggetto dei soli campioni biologici tenuti, con le
#                            sole varianti tenute; gli identificativi non cambiano
#   esclusioni.tsv           i campioni biologici esclusi, con il filtro e il motivo
#   varianti_rimosse.tsv     le varianti rimosse, con il filtro e il motivo
#   filtri_riepilogo.json    l'ordine dei filtri, il denominatore della prevalenza
#                            e cio' che ciascun filtro ha tolto
#   letture_finali.tsv       letture per campione nell'oggetto finale (tracciamento;
#                            zero per i controlli e per i campioni esclusi)
#
# Parametri: decontaminato (ps_decontaminato.rds di S12), tenuti (accession dei
# biologici sopra la soglia di
# profondita'), esclusi_profondita (le esclusioni per profondita', gia' decise),
# profondita_applicata (falso se S11 non ha dato alcuna soglia),
# senza_phylum (filt.remove_na_phylum), taxa_esclusi (filt.exclude_taxa),
# prevalenza (prev.apply), frazione (prev.min_fraction), minimo_conteggio
# (prev.min_count), letture_minime (qc.min_reads_final), letture_nonchimeric
# (accession -> letture senza chimere, per il riepilogo).
#
# L'ORDINE DEI FILTRI.
#   1. profondita' (campioni): deciso dalla fase Python sul tracciamento, prima
#      di tutto, perche' definisce i campioni dell'oggetto finale. Le soglie di
#      S11 valgono sulle letture senza chimere; se S11 non ne ha (nessuna curva
#      valida, o qc.min_reads_mode "none") il filtro non si applica, e il
#      riepilogo lo dichiara in profondita_applicata;
#   2. tassonomico (varianti): senza phylum, se filt.remove_na_phylum, e i taxa
#      di filt.exclude_taxa cercati in ogni rango; un campione che ne resta
#      vuoto non ha segnale batterico ed esce (E-S13-03). I nomi dei taxa si
#      confrontano senza il prefisso di rango che alcuni riferimenti portano
#      ("p__Proteobacteria", "g__Escherichia"): un nome fatto del solo
#      prefisso e' un rango non assegnato;
#   3. prevalenza (varianti): una variante resta se ha almeno prev.min_count
#      letture in almeno min_campioni dei campioni tenuti; un campione che ne
#      resta vuoto esce con il motivo (E-S13-02);
#   4. letture finali (campioni): sotto qc.min_reads_final il campione esce.
# I due filtri sulle varianti commutano: la prevalenza di una variante non
# dipende dalle altre. Una variante rimasta senza letture nei campioni finali
# esce comunque: dopo il filtro 3 con il motivo della prevalenza, dopo il 4 con
# quello delle letture finali. La prevalenza si conta anche nei campioni che il
# filtro 4 togliera': il 4 dipende dalle letture rimaste dopo il 3, e l'ordine
# dichiarato rompe la circolarita'.
#
# IL DENOMINATORE DELLA PREVALENZA sono i campioni biologici tenuti dopo i
# filtri 1 e 2, non tutti i biologici dell'inventario: la presenza si conta in
# quei campioni, e numeratore e denominatore devono riferirsi allo stesso
# insieme. Un campione escluso per profondita' ha una composizione non
# attendibile (la sua fedelta' e' sotto katharoseq.target_sensitivity): la
# presenza di una variante li' non e' un'osservazione su cui fondarsi.
# min_campioni = ceiling(prev.min_fraction * campioni tenuti): i campioni si
# contano interi, e la soglia si arrotonda per eccesso. Si calcola solo qui,
# perche' i campioni tenuti esistono solo dentro la fase.

for (f in c("io_json.R", "errors.R", "letture.R")) {
  source(file.path(Sys.getenv("AMPLICON16S_R_LIB"), f))
}

esegui_fase(function(parametri, cartella) {
  richiedi_pacchetti("phyloseq")
  ps <- readRDS(parametri$decontaminato)
  if (!phyloseq::taxa_are_rows(ps)) stop("l'oggetto di S12 deve avere le varianti sulle righe")
  dati <- methods::as(phyloseq::sample_data(ps), "data.frame")
  tax <- methods::as(phyloseq::tax_table(ps), "matrix")
  campioni <- phyloseq::sample_names(ps)
  tenuti <- as.character(unlist(parametri$tenuti))
  if (!all(tenuti %in% campioni) || any(dati[tenuti, "classe"] != "biologico")) {
    stop("i campioni tenuti devono essere biologici dell'oggetto di S12")
  }
  tenuti <- campioni[campioni %in% tenuti]
  nonchimeric <- unlist(parametri$letture_nonchimeric)

  esclusioni <- lapply(parametri$esclusi_profondita, function(e) {
    c(accession = e$accession, filtro = "profondita", motivo = e$motivo)
  })
  rimosse <- list()
  conteggi <- methods::as(phyloseq::otu_table(ps), "matrix")[, tenuti, drop = FALSE]
  varianti <- rownames(conteggi)

  # ---- 2. Tassonomico --------------------------------------------------------
  # Senza il prefisso di rango ("p__", "o__"): il confronto con
  # filt.exclude_taxa e il riconoscimento del phylum non assegnato non devono
  # dipendere da come il riferimento scrive i nomi.
  senza_prefisso <- sub("^[A-Za-z]__", "", tax)
  dim(senza_prefisso) <- dim(tax)
  dimnames(senza_prefisso) <- dimnames(tax)
  senza_prefisso[!is.na(senza_prefisso) & !nzchar(senza_prefisso)] <- NA
  motivo_tax <- rep("", length(varianti))
  if (isTRUE(parametri$senza_phylum) && "Phylum" %in% colnames(tax)) {
    motivo_tax[is.na(senza_prefisso[varianti, "Phylum"])] <- "phylum non assegnato"
  }
  esclusi_taxa <- sub("^[A-Za-z]__", "", as.character(unlist(parametri$taxa_esclusi)))
  for (i in seq_along(varianti)) {
    # drop = FALSE e i nomi dalle colonne: con una tassonomia di un solo rango
    # la riga sarebbe un valore senza nome, e il motivo non si potrebbe scrivere.
    ranghi <- senza_prefisso[varianti[i], , drop = FALSE]
    trovato <- which(!is.na(ranghi) & ranghi %in% esclusi_taxa)
    if (length(trovato)) {
      voce <- sprintf("%s = %s, in filt.exclude_taxa", colnames(ranghi)[trovato[1]], ranghi[trovato[1]])
      motivo_tax[i] <- if (nzchar(motivo_tax[i])) paste(motivo_tax[i], voce, sep = "; ") else voce
    }
  }
  via_tax <- nzchar(motivo_tax)
  rimosse$tassonomico <- data.frame(asv_id = varianti[via_tax], motivo = motivo_tax[via_tax],
                                    letture = rowSums(conteggi[via_tax, , drop = FALSE]))
  conteggi <- conteggi[!via_tax, , drop = FALSE]
  vuoti_tax <- colnames(conteggi)[colSums(conteggi) == 0]
  for (a in vuoti_tax) {
    esclusioni[[length(esclusioni) + 1L]] <- c(
      accession = a, filtro = "tassonomico",
      motivo = "nessuna lettura dopo il filtro tassonomico: solo varianti senza phylum o di filt.exclude_taxa")
  }
  conteggi <- conteggi[, setdiff(colnames(conteggi), vuoti_tax), drop = FALSE]

  # ---- 3. Prevalenza ---------------------------------------------------------
  n_tenuti <- ncol(conteggi)
  # Arrotondato prima dell'intero superiore: in virgola mobile 0,07 x 100 vale
  # 7,000000000000001, e l'intero superiore sarebbe 8 invece di 7.
  min_campioni <- ceiling(round(as.numeric(parametri$frazione) * n_tenuti, 9))
  presenze <- rowSums(conteggi >= as.numeric(parametri$minimo_conteggio))
  if (isTRUE(parametri$prevalenza)) {
    via_prev <- presenze < min_campioni
    motivo_prev <- sprintf("presente con almeno %d letture in %d campioni su %d, meno di %d",
                           as.integer(parametri$minimo_conteggio), presenze, n_tenuti, min_campioni)
  } else {
    via_prev <- rep(FALSE, nrow(conteggi))
    motivo_prev <- rep("", nrow(conteggi))
  }
  assenti <- !via_prev & rowSums(conteggi) == 0
  motivo_prev[assenti] <- "nessuna lettura nei campioni tenuti"
  via_prev <- via_prev | assenti
  rimosse$prevalenza <- data.frame(asv_id = rownames(conteggi)[via_prev], motivo = motivo_prev[via_prev],
                                   letture = rowSums(conteggi[via_prev, , drop = FALSE]))
  conteggi <- conteggi[!via_prev, , drop = FALSE]
  vuoti_prev <- colnames(conteggi)[colSums(conteggi) == 0]
  for (a in vuoti_prev) {
    esclusioni[[length(esclusioni) + 1L]] <- c(
      accession = a, filtro = "prevalenza",
      motivo = "nessuna lettura dopo il filtro di prevalenza: solo varianti rare")
  }
  conteggi <- conteggi[, setdiff(colnames(conteggi), vuoti_prev), drop = FALSE]

  # ---- 4. Letture finali ----------------------------------------------------
  minime <- as.numeric(parametri$letture_minime)
  pochi <- colnames(conteggi)[colSums(conteggi) < minime]
  for (a in pochi) {
    esclusioni[[length(esclusioni) + 1L]] <- c(
      accession = a, filtro = "letture_finali",
      motivo = sprintf("%.0f letture finali, meno di qc.min_reads_final (%.0f)",
                       sum(conteggi[, a]), minime))
  }
  conteggi <- conteggi[, setdiff(colnames(conteggi), pochi), drop = FALSE]
  if (ncol(conteggi) == 0L) {
    errore_catalogo("E-S13-04", sprintf(
      "%d campioni biologici, tutti esclusi: %s", sum(dati$classe == "biologico"),
      paste(vapply(c("profondita", "tassonomico", "prevalenza", "letture_finali"), function(f) {
        sprintf("%s %d", f, sum(vapply(esclusioni, function(e) e[["filtro"]] == f, logical(1))))
      }, character(1)), collapse = ", ")))
  }
  orfane <- rowSums(conteggi) == 0
  rimosse$letture_finali <- data.frame(
    asv_id = rownames(conteggi)[orfane], letture = rep(0, sum(orfane)),
    motivo = rep("nessuna lettura nei campioni finali, tolti quelli sotto qc.min_reads_final",
                 sum(orfane)))
  conteggi <- conteggi[!orfane, , drop = FALSE]

  # ---- L'oggetto filtrato ---------------------------------------------------
  finale <- phyloseq::prune_samples(colnames(conteggi), ps)
  finale <- phyloseq::prune_taxa(rownames(conteggi), finale)
  if (!identical(phyloseq::taxa_names(finale), rownames(conteggi)) ||
      !identical(phyloseq::sample_names(finale), colnames(conteggi))) {
    stop("l'oggetto finale non corrisponde alla tabella filtrata")
  }
  salva_rds(finale, file.path(cartella, "ps_filtrato.rds"))

  # ---- Esclusioni, varianti rimosse, tracciamento ---------------------------
  esclusioni <- if (length(esclusioni)) do.call(rbind, esclusioni) else
    matrix(character(), ncol = 3, dimnames = list(NULL, c("accession", "filtro", "motivo")))
  ordine_inv <- match(esclusioni[, "accession"], campioni)
  esclusioni <- esclusioni[order(ordine_inv), , drop = FALSE]
  scrivi_atomico(
    c("accession\tsample_name\tpiastra\tfiltro\tmotivo",
      if (nrow(esclusioni)) sprintf("%s\t%s\t%s\t%s\t%s", esclusioni[, "accession"],
              dati[esclusioni[, "accession"], "sample_name"],
              ifelse(is.na(dati[esclusioni[, "accession"], "piastra"]), "",
                     dati[esclusioni[, "accession"], "piastra"]),
              esclusioni[, "filtro"], esclusioni[, "motivo"])),
    file.path(cartella, "esclusioni.tsv")
  )
  righe_v <- unlist(lapply(names(rimosse), function(f) {
    r <- rimosse[[f]]
    if (!nrow(r)) return(character())
    sprintf("%s\t%s\t%s\t%.0f", r$asv_id, f, r$motivo, r$letture)
  }))
  scrivi_atomico(c("asv_id\tfiltro\tmotivo\tletture_biologici_tenuti", righe_v),
                 file.path(cartella, "varianti_rimosse.tsv"))
  finali <- stats::setNames(numeric(length(campioni)), campioni)
  finali[colnames(conteggi)] <- colSums(conteggi)
  finali <- finali[sort(names(finali), method = "radix")]

  conta <- function(filtro) sum(esclusioni[, "filtro"] == filtro)
  classe_di <- function(accessioni) as.list(table(factor(dati[accessioni, "classe"],
    levels = c("biologico", "controllo_negativo", "controllo_positivo"))))
  tutte <- methods::as(phyloseq::otu_table(ps), "matrix")
  per_classe_varianti <- function(ids) {
    lapply(c(biologico = "biologico", controllo_negativo = "controllo_negativo",
             controllo_positivo = "controllo_positivo"), function(cl) {
      colonne <- dati$classe == cl
      tot <- sum(tutte[, colonne])
      if (tot == 0) NULL else signif(sum(tutte[ids, colonne]) / tot, 6)
    })
  }
  scrivi_json(
    list(
      ordine = list("profondita", "tassonomico", "prevalenza", "letture_finali"),
      profondita_applicata = isTRUE(parametri$profondita_applicata),
      campioni = list(
        biologici = sum(dati$classe == "biologico"),
        controlli_a_parte = sum(dati$classe != "biologico"),
        esclusi = list(profondita = conta("profondita"), tassonomico = conta("tassonomico"),
                       prevalenza = conta("prevalenza"), letture_finali = conta("letture_finali")),
        finali = ncol(conteggi),
        classi_finali = classe_di(colnames(conteggi))),
      prevalenza = list(applicata = isTRUE(parametri$prevalenza), denominatore = n_tenuti,
                        min_campioni = min_campioni,
                        min_conteggio = as.integer(parametri$minimo_conteggio)),
      varianti = list(
        iniziali = length(varianti),
        rimosse = list(tassonomico = nrow(rimosse$tassonomico), prevalenza = nrow(rimosse$prevalenza),
                       letture_finali = nrow(rimosse$letture_finali)),
        finali = nrow(conteggi),
        frazione_letture_rimosse = list(
          tassonomico = per_classe_varianti(rimosse$tassonomico$asv_id),
          prevalenza = per_classe_varianti(rimosse$prevalenza$asv_id),
          letture_finali = per_classe_varianti(rimosse$letture_finali$asv_id))),
      letture = list(
        nonchimeric_finali = sum(as.numeric(nonchimeric[colnames(conteggi)])),
        decontaminate_finali = sum(methods::as(phyloseq::otu_table(ps), "matrix")[, colnames(conteggi)]),
        finali = sum(conteggi)),
      campioni_svuotati = list(tassonomico = as.list(vuoti_tax), prevalenza = as.list(vuoti_prev))
    ),
    file.path(cartella, "filtri_riepilogo.json")
  )
  c("ps_filtrato.rds", "esclusioni.tsv", "varianti_rimosse.tsv",
    "filtri_riepilogo.json", traccia_letture(finali, "finali", cartella))
})
