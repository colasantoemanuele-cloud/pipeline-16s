# S6 - rimozione delle chimere (dada2::removeBimeraDenovo).
#
# Scrive in 07_chimera/:
#
#   tabella_senza_chimere.rds   la tabella di S5 senza le varianti chimeriche
#   chimere.tsv                 le varianti chimeriche: sequenza, letture,
#                               campioni in cui compaiono, e se la variante e'
#                               stata tolta dalla tabella per intero
#   chimere_per_campione.tsv    per campione: letture e varianti prima,
#                               letture e varianti chimeriche
#   chimere_per_gruppo.tsv      lo stesso per gruppo di campioni, con le
#                               varianti contate distinte nel gruppo
#   letture_senza_chimere.tsv   letture per campione dopo la rimozione
#                               (tracciamento)
#
# La tabella prima della rimozione non si ricopia: e' tabella.rds di S5, in
# 06_seqtab/, e la fase ne registra il riferimento con il checksum. Con
# chimere.tsv la differenza fra le due tabelle e' esplicita: con i metodi
# "consensus" e "pooled" una variante e' chimerica ovunque o in nessun
# campione, e la tabella senza chimere e' quella di S5 meno le colonne
# elencate.
#
# Parametri: tabella (il file di S5), senza_letture (accession senza letture
# filtrate, zero nel tracciamento), gruppi (accession -> gruppo), method,
# min_fold_parent_over_abundance, min_parent_abundance, min_sample_fraction
# (solo per "consensus"), allow_one_off, processi. Gli altri argomenti di
# dada2 (ignoreNNegatives, minOneOffParentDistance, maxShift) restano ai
# valori della libreria.

for (f in c("io_json.R", "errors.R", "letture.R")) {
  source(file.path(Sys.getenv("AMPLICON16S_R_LIB"), f))
}

esegui_fase(function(parametri, cartella) {
  richiedi_pacchetti("dada2")

  tabella <- readRDS(parametri$tabella)
  metodo <- parametri$method
  argomenti <- list(
    tabella,
    method = metodo,
    minFoldParentOverAbundance = as.numeric(parametri$min_fold_parent_over_abundance),
    minParentAbundance = as.numeric(parametri$min_parent_abundance),
    allowOneOff = isTRUE(parametri$allow_one_off),
    multithread = as.integer(parametri$processi),
    verbose = FALSE
  )
  # minSampleFraction esiste solo per il metodo consensus: agli altri, che
  # passano gli argomenti a isBimeraDenovo, sarebbe un argomento sconosciuto.
  if (identical(metodo, "consensus")) {
    argomenti$minSampleFraction <- as.numeric(parametri$min_sample_fraction)
  }
  # Razionale biologico: le chimere bimeriche si originano durante la PCR quando
  # un filamento incompleto si appaia come innesco illegittimo su un templato
  # eterologo nel ciclo successivo, generando una sequenza ibrida composta dal
  # segmento 5' di una variante parentale e dal segmento 3' di un'altra. Poiche'
  # i due genitori devono essere piu' abbondanti del ricombinante tardivo
  # (minFoldParentOverAbundance), la rimozione de novo scarta le sequenze spurie
  # preservando le varianti biologiche genuine (ASV).
  senza <- do.call(dada2::removeBimeraDenovo, argomenti)

  if (!identical(rownames(senza), rownames(tabella))) {
    stop("la tabella senza chimere non ha i campioni della tabella di S5")
  }

  # Le letture chimeriche, variante per variante: le colonne tolte per intero
  # e, solo con "per-sample", le letture tolte a una variante in alcuni
  # campioni. Le copie riguardano le sole colonne chimeriche, non la tabella.
  chimeriche <- tabella[, setdiff(colnames(tabella), colnames(senza)), drop = FALSE]
  if (identical(metodo, "per-sample")) {
    parziali <- tabella[, colnames(senza), drop = FALSE] - senza
    chimeriche <- cbind(chimeriche, parziali[, colSums(parziali) > 0, drop = FALSE])
    rm(parziali)
  }

  per_variante <- colSums(chimeriche)
  scrivi_atomico(
    c("sequenza\tletture\tcampioni\trimossa",
      sprintf("%s\t%.0f\t%d\t%s", colnames(chimeriche), per_variante,
              as.integer(colSums(chimeriche > 0)),
              ifelse(colnames(chimeriche) %in% colnames(senza), "no", "si"))),
    file.path(cartella, "chimere.tsv")
  )

  prima <- rowSums(tabella)
  dopo <- rowSums(senza)
  scrivi_atomico(
    c("campione\tletture\tletture_chimeriche\tvarianti\tvarianti_chimeriche",
      sprintf("%s\t%.0f\t%.0f\t%d\t%d", rownames(tabella), prima, prima - dopo,
              rowSums(tabella > 0), rowSums(chimeriche > 0))),
    file.path(cartella, "chimere_per_campione.tsv")
  )

  # Per gruppo di campioni (la classe, per la fase): le varianti si contano
  # distinte, una variante presente in piu' campioni del gruppo conta una volta.
  gruppi <- unlist(parametri$gruppi)[rownames(tabella)]
  righe <- vapply(sort(unique(gruppi), method = "radix"), function(g) {
    r <- rownames(tabella)[gruppi == g]
    sprintf("%s\t%d\t%.0f\t%.0f\t%d\t%d", g, length(r), sum(prima[r]),
            sum(prima[r] - dopo[r]),
            sum(colSums(tabella[r, , drop = FALSE]) > 0),
            sum(colSums(chimeriche[r, , drop = FALSE]) > 0))
  }, character(1))
  scrivi_atomico(
    c("gruppo\tcampioni\tletture\tletture_chimeriche\tvarianti\tvarianti_chimeriche", righe),
    file.path(cartella, "chimere_per_gruppo.tsv")
  )
  rm(chimeriche)

  senza_letture <- as.character(unlist(parametri$senza_letture))
  tracciate <- c(dopo, stats::setNames(numeric(length(senza_letture)), senza_letture))
  tracciate <- tracciate[sort(names(tracciate), method = "radix")]

  salva_rds(senza, file.path(cartella, "tabella_senza_chimere.rds"))
  c("tabella_senza_chimere.rds", "chimere.tsv", "chimere_per_campione.tsv",
    "chimere_per_gruppo.tsv",
    traccia_letture(tracciate, "senza_chimere", cartella))
})
