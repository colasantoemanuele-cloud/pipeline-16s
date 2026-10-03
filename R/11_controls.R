# S11 - validazione della corsa dai controlli positivi (KatharoSeq).
#
# Scrive in 11_controls/, condivisa con S12:
#
#   positivi.tsv             per controllo positivo: piastra, cellule (il livello
#                            di diluizione), profondita', letture del taxon atteso,
#                            fedelta', z per il livello, conformita' e motivo, e se
#                            entra nella curva
#   curve.tsv                per curva (aggregata e per piastra): punti, h, x50,
#                            k', R^2, RMSE, soglia, punto medio, osservazioni fra
#                            punto medio e soglia, se e' valida e perche'
#   soglia.json              la soglia di profondita' scelta, per piastra, con lo
#                            stadio delle letture a cui si applica, l'origine
#                            (curva o ripiego) e il motivo di ogni ripiego
#   profondita_campioni.tsv  per campione: letture grezze e allo stadio della
#                            soglia, la soglia che gli si applica e se vi cade
#                            sotto. E' una misura: il filtro e' di S13
#   riepilogo.json           conformita' dei controlli, scelta della curva e conteggi
#
# Parametri: oggetto (ps_integrato.rds di S10), colonna_cellule (il nome nell'oggetto
# della colonna dei livelli, o null) e motivo_colonna, rango, taxon, sensibilita,
# min_r2, min_positivi, modo (qc.min_reads_mode), ripiego (qc.min_reads_raw),
# letture_grezze (letture_prefiltro.tsv di S2).
#
# La curva, la bonta' e la conformita' sono in R/lib/katharoseq.R.
#
# LA SCELTA FRA CURVA AGGREGATA E CURVE PER PIASTRA. Si adattano entrambe sui
# controlli con un livello noto, con letture, e non giudicati non conformi. Una
# curva e' valida se converge, ha almeno min_positivi punti, ha R^2 non
# inferiore a min_r2, la sua soglia cade dentro l'intervallo di profondita'
# dei suoi punti (una soglia estrapolata non e' sostenuta da alcun controllo)
# ed e' determinata dai dati: almeno un'osservazione fra il punto medio della
# curva e la soglia (R/lib/katharoseq.R).
# La bonta' del modello per piastra e' l'R^2 complessivo delle previsioni di
# ciascuna piastra sui propri punti, calcolato su tutti i punti delle piastre la
# cui curva converge. Si sceglie il modello con la bonta' maggiore fra quelli
# ammissibili (l'aggregato se valido, il modello per piastra se la sua bonta'
# raggiunge min_r2 e almeno una piastra ha una curva valida); a parita',
# l'aggregato. Con il modello per piastra ogni piastra usa la propria soglia,
# e una piastra senza curva valida ripiega su qc.min_reads_raw. Se nessun
# modello e' ammissibile, ripiegano tutte.
#
# GLI STADI. La soglia derivata vale sulle letture dell'oggetto integrato
# (katharoseq.read_stage, "nonchimeric"); il ripiego qc.min_reads_raw vale sulle
# letture grezze ("raw"). Le due grandezze non sono confrontabili, e ogni
# soglia scritta porta il suo stadio.

for (f in c("io_json.R", "errors.R", "letture.R", "katharoseq.R")) {
  source(file.path(Sys.getenv("AMPLICON16S_R_LIB"), f))
}

STADI <- list(
  nonchimeric = paste0("letture dell'oggetto integrato di S10: senza chimere (S6) e ",
                       "dopo il filtro di lunghezza (S7)"),
  raw = "letture grezze, in ingresso al filtro (S2, letture_prefiltro.tsv)"
)

numero <- function(x, cifre = 6L) {
  ifelse(is.na(x), "", trimws(formatC(signif(x, cifre), digits = cifre, format = "g")))
}

esegui_fase(function(parametri, cartella) {
  richiedi_pacchetti("phyloseq")
  sensibilita <- as.numeric(parametri$sensibilita)
  min_r2 <- as.numeric(parametri$min_r2)
  min_positivi <- as.integer(parametri$min_positivi)
  ripiego <- as.numeric(parametri$ripiego)
  stadio <- parametri$stadio

  ps <- readRDS(parametri$oggetto)
  dati <- methods::as(phyloseq::sample_data(ps), "data.frame")
  conteggi <- methods::as(phyloseq::otu_table(ps), "matrix")
  if (!phyloseq::taxa_are_rows(ps)) conteggi <- t(conteggi)
  tax <- methods::as(phyloseq::tax_table(ps), "matrix")
  campioni <- colnames(conteggi)
  if (!identical(rownames(dati), campioni) || !identical(rownames(tax), rownames(conteggi))) {
    stop("l'oggetto di S10 non ha campioni e varianti allineati")
  }
  profondita <- colSums(conteggi)
  grezze <- leggi_letture(parametri$letture_grezze)[campioni]
  grezze[is.na(grezze)] <- 0L
  piastra <- dati$piastra

  # ---- Fedelta' dei controlli positivi -------------------------------------
  motivo_globale <- character()
  rango <- parametri$rango
  if (rango %in% colnames(tax)) {
    bersaglio <- !is.na(tax[, rango]) & tax[, rango] == parametri$taxon
  } else {
    bersaglio <- rep(FALSE, nrow(tax))
    motivo_globale <- c(motivo_globale, sprintf("rango %s assente nella tassonomia", rango))
  }
  positivi <- which(dati$classe == "controllo_positivo")
  n_pos <- profondita[positivi]
  target <- colSums(conteggi[bersaglio, positivi, drop = FALSE])
  fedelta <- ifelse(n_pos > 0, target / pmax(n_pos, 1), NA_real_)

  colonna <- parametri$colonna_cellule
  if (is.null(colonna)) {
    cellule <- rep(NA_real_, length(positivi))
    motivo_globale <- c(motivo_globale, parametri$motivo_colonna)
  } else {
    cellule <- suppressWarnings(as.numeric(dati[[colonna]][positivi]))
    cellule[!is.na(cellule) & cellule <= 0] <- NA_real_
    if (all(is.na(cellule))) {
      motivo_globale <- c(motivo_globale, sprintf(
        "nessun controllo positivo ha un numero di cellule nella colonna %s", colonna))
    }
  }

  conf <- conformita(cellule, fedelta, n_pos, min_positivi)
  nella_curva <- !is.na(cellule) & n_pos >= 1 & !is.na(fedelta) & conf$esito != "non conforme"
  if (length(motivo_globale) == 0L && sum(nella_curva) < min_positivi) {
    motivo_globale <- sprintf(
      "%d controlli positivi utilizzabili per la curva, meno di ctrl.min_positives (%d)",
      sum(nella_curva), min_positivi)
  }

  # ---- Le curve -------------------------------------------------------------
  valuta <- function(nome, indici) {
    curva <- adatta_curva(n_pos[indici], fedelta[indici], sensibilita)
    motivi <- character()
    if (!curva$converge) motivi <- curva$motivo
    if (curva$punti < min_positivi) {
      motivi <- c(motivi, sprintf("%d punti, meno di ctrl.min_positives (%d)",
                                  curva$punti, min_positivi))
    }
    if (curva$converge && curva$r2 < min_r2) {
      motivi <- c(motivi, sprintf("R^2 %.4f sotto katharoseq.min_r2 (%s)", curva$r2,
                                  format(min_r2)))
    }
    if (curva$converge && (curva$soglia < min(n_pos[indici]) ||
                           curva$soglia > max(n_pos[indici]))) {
      motivi <- c(motivi, sprintf(
        "soglia %.0f fuori dalle profondita' dei punti (%.0f-%.0f): sarebbe estrapolata",
        curva$soglia, min(n_pos[indici]), max(n_pos[indici])))
    }
    sostegno <- soglia_determinata(n_pos[indici], curva)
    if (curva$converge && !sostegno$determinata) motivi <- c(motivi, sostegno$motivo)
    c(list(modello = nome, indici = indici, valida = length(motivi) == 0L,
           motivo_validita = paste(motivi, collapse = "; "),
           punto_medio = sostegno$punto_medio, osservazioni_transizione = sostegno$osservazioni),
      curva)
  }
  curve <- list()
  if (length(motivo_globale) == 0L) {
    punti <- which(nella_curva)
    curve[["aggregato"]] <- valuta("aggregato", punti)
    for (p in sort(unique(piastra[positivi][punti]), method = "radix")) {
      curve[[paste0("piastra ", p)]] <- valuta(
        paste0("piastra ", p), punti[piastra[positivi][punti] == p])
    }
  }
  per_piastra <- Filter(function(cv) cv$modello != "aggregato" && cv$converge, curve)
  r2_per_piastra <- NA_real_
  if (length(per_piastra)) {
    osservate <- unlist(lapply(per_piastra, function(cv) fedelta[cv$indici]))
    previste <- unlist(lapply(per_piastra, function(cv) {
      1 / (1 + (cv$x50 / log10(n_pos[cv$indici]))^cv$h)
    }))
    r2_per_piastra <- 1 - sum((osservate - previste)^2) /
      sum((osservate - mean(osservate))^2)
  }
  aggregato <- curve[["aggregato"]]
  r2_aggregato <- if (is.null(aggregato) || !aggregato$converge) NA_real_ else aggregato$r2

  # ---- La scelta ------------------------------------------------------------
  ammesso_agg <- !is.null(aggregato) && aggregato$valida
  ammesso_pp <- !is.na(r2_per_piastra) && r2_per_piastra >= min_r2 &&
    any(vapply(per_piastra, function(cv) cv$valida, logical(1)))
  if (parametri$modo == "fixed") {
    scelta <- "nessuno"
    motivo_scelta <- "qc.min_reads_mode e' fixed: si usa qc.min_reads_raw senza curva"
  } else if (length(motivo_globale)) {
    scelta <- "nessuno"
    motivo_scelta <- paste(motivo_globale, collapse = "; ")
  } else if (!ammesso_agg && !ammesso_pp) {
    scelta <- "nessuno"
    motivo_scelta <- sprintf(
      "nessun modello ammissibile: aggregato %s; per piastra R^2 %s",
      if (is.null(aggregato)) "assente" else if (aggregato$valida) "valido" else aggregato$motivo_validita,
      if (is.na(r2_per_piastra)) "non calcolabile" else sprintf("%.4f", r2_per_piastra))
  } else if (ammesso_pp && (!ammesso_agg || r2_per_piastra > r2_aggregato)) {
    scelta <- "per_piastra"
    motivo_scelta <- sprintf(
      "R^2 per piastra %.4f contro %s dell'aggregato%s", r2_per_piastra,
      if (is.na(r2_aggregato)) "nessuno" else sprintf("%.4f", r2_aggregato),
      if (ammesso_agg) "" else sprintf(" (aggregato non valido: %s)", aggregato$motivo_validita))
  } else {
    scelta <- "aggregato"
    motivo_scelta <- sprintf(
      "R^2 dell'aggregato %.4f contro %s per piastra%s", r2_aggregato,
      if (is.na(r2_per_piastra)) "nessuno" else sprintf("%.4f", r2_per_piastra),
      if (ammesso_pp) "" else " (modello per piastra non ammissibile)")
  }

  # ---- Una soglia per piastra, con il suo stadio ----------------------------
  da_curva <- function(cv, origine) {
    list(valore = ceiling(cv$soglia), stadio = stadio, descrizione_stadio = STADI[[stadio]],
         origine = origine, motivo = "")
  }
  di_ripiego <- function(motivo) {
    list(valore = ripiego, stadio = "raw", descrizione_stadio = STADI$raw,
         origine = "ripiego: qc.min_reads_raw", motivo = motivo)
  }
  piastre <- sort(unique(piastra[!is.na(piastra)]), method = "radix")
  soglie <- list()
  for (p in piastre) {
    cv <- curve[[paste0("piastra ", p)]]
    soglie[[p]] <- switch(
      scelta,
      aggregato = da_curva(aggregato, "curva aggregata"),
      per_piastra = if (!is.null(cv) && cv$valida) da_curva(cv, "curva della piastra") else
        di_ripiego(if (is.null(cv)) "nessun controllo positivo utilizzabile nella piastra" else
          paste("curva della piastra non valida:", cv$motivo_validita)),
      nessuno = di_ripiego(motivo_scelta)
    )
  }
  senza_piastra <- if (scelta == "aggregato") da_curva(aggregato, "curva aggregata") else
    di_ripiego(if (scelta == "nessuno") motivo_scelta else "campione senza piastra")
  degrada <- parametri$modo != "fixed" && any(vapply(
    c(soglie, if (anyNA(piastra)) list(senza_piastra)),
    function(s) s$stadio == "raw", logical(1)))

  # ---- Misura sui campioni: chi cadrebbe sotto la soglia -------------------
  applicata <- lapply(seq_along(campioni), function(i) {
    if (is.na(piastra[i])) senza_piastra else soglie[[piastra[i]]]
  })
  valore <- vapply(applicata, function(s) s$valore, numeric(1))
  stadi <- vapply(applicata, function(s) s$stadio, character(1))
  confrontate <- ifelse(stadi == "raw", grezze, profondita)
  sotto <- confrontate < valore
  scrivi_atomico(
    c("accession\tclasse\tpiastra\tletture_grezze\tletture_nonchimeric\tsoglia\tstadio\tsotto_soglia",
      sprintf("%s\t%s\t%s\t%.0f\t%.0f\t%.0f\t%s\t%s", campioni, dati$classe,
              ifelse(is.na(piastra), "", piastra), grezze, profondita, valore, stadi,
              ifelse(sotto, "si", "no"))),
    file.path(cartella, "profondita_campioni.tsv")
  )

  # ---- Artefatti ------------------------------------------------------------
  scrivi_atomico(
    c(paste("accession", "sample_name", "piastra", "cellule", "profondita",
            "letture_target", "fedelta", "z_fedelta", "z_profondita", "conformita",
            "motivo", "nella_curva", sep = "\t"),
      sprintf("%s\t%s\t%s\t%s\t%.0f\t%.0f\t%s\t%s\t%s\t%s\t%s\t%s",
              campioni[positivi], dati$sample_name[positivi],
              ifelse(is.na(piastra[positivi]), "", piastra[positivi]),
              numero(cellule), n_pos, target, numero(fedelta),
              numero(conf$z_fedelta, 4L), numero(conf$z_profondita, 4L), conf$esito,
              conf$motivo, ifelse(nella_curva, "si", "no"))),
    file.path(cartella, "positivi.tsv")
  )
  righe_curve <- vapply(curve, function(cv) {
    sprintf("%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%d\t%s\t%s", cv$modello, cv$punti,
            numero(cv$h), numero(cv$x50), numero(cv$k), numero(cv$r2, 4L),
            numero(cv$rmse, 4L),
            if (is.na(cv$soglia)) "" else sprintf("%.0f", ceiling(cv$soglia)),
            if (is.na(cv$punto_medio)) "" else sprintf("%.0f", cv$punto_medio),
            cv$osservazioni_transizione, if (cv$valida) "si" else "no", cv$motivo_validita)
  }, character(1))
  scrivi_atomico(
    c(paste0("modello\tpunti\th\tx50\tk\tr2\trmse\tsoglia\tpunto_medio\t",
             "osservazioni_transizione\tvalida\tmotivo"), unname(righe_curve)),
    file.path(cartella, "curve.tsv")
  )

  valutabili <- conf$esito != "non valutabile"
  conformi <- sum(conf$esito == "conforme")
  frazione <- if (any(valutabili)) conformi / sum(valutabili) else NA_real_
  ripieghi <- Filter(function(s) s$stadio == "raw",
                     c(soglie, if (anyNA(piastra)) list(senza_piastra = senza_piastra)))
  scrivi_json(
    list(
      modo = parametri$modo,
      sensibilita = sensibilita,
      modello_curva = "allosteric_sigmoid: f = x^h / (k' + x^h), x = log10(profondita'), k' = x50^h",
      scelta = scelta,
      motivo_scelta = motivo_scelta,
      r2_aggregato = if (is.na(r2_aggregato)) NULL else signif(r2_aggregato, 6),
      r2_per_piastra = if (is.na(r2_per_piastra)) NULL else signif(r2_per_piastra, 6),
      stadi = STADI,
      per_piastra = soglie,
      senza_piastra = senza_piastra,
      ripiego = list(parametro = "qc.min_reads_raw", valore = ripiego, stadio = "raw",
                     degradazione = degrada,
                     motivi = unname(lapply(names(ripieghi), function(n) {
                       list(piastra = n, motivo = ripieghi[[n]]$motivo)
                     })))
    ),
    file.path(cartella, "soglia.json")
  )
  biologici <- dati$classe == "biologico"
  scrivi_json(
    list(
      positivi = length(positivi),
      con_livello = sum(!is.na(cellule)),
      nella_curva = sum(nella_curva),
      conformi = conformi,
      non_conformi = sum(conf$esito == "non conforme"),
      non_valutabili = sum(!valutabili),
      frazione_conformi = if (is.na(frazione)) NULL else signif(frazione, 6),
      non_conformi_elenco = as.list(dati$sample_name[positivi][conf$esito == "non conforme"]),
      scelta = scelta,
      motivo_scelta = motivo_scelta,
      degradazione = degrada,
      biologici = sum(biologici),
      biologici_sotto_soglia = sum(sotto & biologici),
      biologici_sotto_soglia_per_piastra = as.list(vapply(
        piastre, function(p) sum(sotto & biologici & !is.na(piastra) & piastra == p),
        integer(1)))
    ),
    file.path(cartella, "riepilogo.json")
  )
  c("positivi.tsv", "curve.tsv", "soglia.json", "profondita_campioni.tsv", "riepilogo.json")
})
