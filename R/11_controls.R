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
#   soglia.json              la soglia di profondita' di ogni piastra e dei campioni
#                            senza piastra: valore, stadio delle letture, origine
#                            (propria, aggregata, mediana, nessuna) e motivo; il
#                            modello scelto con gli AIC, i punti comuni e l'R^2
#                            dichiarato; la mediana, arrotondata e no
#   profondita_campioni.tsv  per campione: letture allo stadio della soglia, la
#                            soglia che gli si applica con la sua origine e se vi
#                            cade sotto. E' una misura: il filtro e' di S13
#   riepilogo.json           conformita' dei controlli, scelta della curva e conteggi
#
# Parametri: oggetto (ps_integrato.rds di S10), colonna_cellule (il nome nell'oggetto
# della colonna dei livelli, o null) e motivo_colonna, rango, taxon, sensibilita,
# min_r2, min_positivi, stadio (katharoseq.read_stage), modo (qc.min_reads_mode).
#
# La curva, la bonta', l'AIC e la conformita' sono in R/lib/katharoseq.R.
#
# LO STADIO. Ogni soglia vale sulle letture dell'oggetto integrato
# (katharoseq.read_stage, "nonchimeric"), lo stadio a cui profondita' e fedelta'
# dei controlli sono misurate. Nessuna soglia si applica alle letture grezze:
# una soglia su un'altra grandezza non sarebbe confrontabile fra le piastre, e
# i campioni di una piastra passerebbero un filtro piu' largo delle altre (un
# effetto di lotto introdotto dal filtro).
#
# I PUNTI COMUNI. La curva aggregata e le curve per piastra si adattano sugli
# stessi punti: i controlli con un livello noto, con piu' di una lettura e non
# giudicati non conformi, delle piastre che ne hanno almeno min_positivi. Solo a
# parita' di punti le verosimiglianze dei due modelli sono confrontabili. Un
# controllo senza piastra, o di una piastra con meno punti, non entra in alcuna
# curva; uno con una sola lettura nemmeno, perche' la curva e' definita sul
# logaritmo della profondita', che per una lettura e' zero.
#
# LA VALIDITA'. Una curva e' valida se converge, ha almeno min_positivi punti, ha
# R^2 non inferiore a min_r2, la sua soglia cade dentro l'intervallo di
# profondita' dei suoi punti (una soglia estrapolata non e' sostenuta da alcun
# controllo) ed e' determinata dai dati: almeno un'osservazione fra il punto
# medio della curva e la soglia (R/lib/katharoseq.R).
#
# LA SCELTA DEL MODELLO. Si preferisce il modello con l'AIC minore; l'AIC del
# modello per piastra e' la somma degli AIC delle sue curve. L'aggregato e'
# preferito solo se il suo AIC e' strettamente minore: a parita' (una sola
# piastra), o se una delle curve non converge e la somma non e' calcolabile,
# vale il modello per piastra, che usa la curva di ogni piastra dove esiste.
#
# LA SOGLIA DI OGNI PIASTRA, nell'ordine:
#   1. se l'aggregato e' preferito ed e' valido, tutte le piastre usano la
#      soglia aggregata (origine "aggregata");
#   2. altrimenti una piastra con la propria curva valida usa la propria
#      (origine "propria");
#   3. una piastra senza curva valida usa la soglia aggregata, se l'aggregato e'
#      valido (origine "aggregata");
#   4. altrimenti la mediana delle soglie delle piastre con curva propria
#      valida, arrotondata all'intero superiore (origine "mediana"): una piastra
#      i cui controlli non determinano una soglia riceve cosi' un filtro dello
#      stesso ordine di quello delle altre;
#   5. un campione senza piastra segue i passi 3 e 4;
#   6. se nessuna curva e' valida non c'e' alcuna soglia (origine "nessuna"):
#      resta il solo qc.min_reads_final di S13.
# La bonta' dichiarata e' l'R^2 dell'aggregato al passo 1; altrimenti l'R^2
# complessivo delle sole piastre che usano la propria curva, ciascuna con le
# previsioni della propria.
#
# SENZA CONTROLLI POSITIVI, senza la colonna dei livelli o senza una piastra
# con abbastanza punti non si adatta alcuna curva: la scelta e' "nessuno" con il
# motivo. Con qc.min_reads_mode "none" le curve si adattano e si riportano come
# diagnostica, ma nessuna soglia si applica.

for (f in c("io_json.R", "errors.R", "katharoseq.R")) {
  source(file.path(Sys.getenv("AMPLICON16S_R_LIB"), f))
}

STADI <- list(
  nonchimeric = paste0("letture dell'oggetto integrato di S10: senza chimere (S6) e ",
                       "dopo il filtro di lunghezza (S7)")
)

numero <- function(x, cifre = 6L) {
  ifelse(is.na(x), "", trimws(formatC(signif(x, cifre), digits = cifre, format = "g")))
}

esegui_fase(function(parametri, cartella) {
  richiedi_pacchetti("phyloseq")
  sensibilita <- as.numeric(parametri$sensibilita)
  min_r2 <- as.numeric(parametri$min_r2)
  min_positivi <- as.integer(parametri$min_positivi)
  senza_soglia <- identical(parametri$modo, "none")
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
  piastra <- dati$piastra

  # ---- Fedelta' dei controlli positivi -------------------------------------
  motivo_globale <- character()
  rango <- parametri$rango
  positivi <- which(dati$classe == "controllo_positivo")
  if (length(positivi) == 0L) {
    bersaglio <- rep(FALSE, nrow(tax))
    motivo_globale <- "nessun controllo positivo nell'inventario"
  } else if (is.null(parametri$taxon)) {
    bersaglio <- rep(FALSE, nrow(tax))
    motivo_globale <- "katharoseq.target_taxon non dichiarato"
  } else if (rango %in% colnames(tax)) {
    bersaglio <- !is.na(tax[, rango]) & tax[, rango] == parametri$taxon
  } else {
    bersaglio <- rep(FALSE, nrow(tax))
    motivo_globale <- c(motivo_globale, sprintf("rango %s assente nella tassonomia", rango))
  }
  n_pos <- profondita[positivi]
  target <- colSums(conteggi[bersaglio, positivi, drop = FALSE])
  fedelta <- ifelse(n_pos > 0, target / pmax(n_pos, 1), NA_real_)

  colonna <- parametri$colonna_cellule
  if (length(positivi) == 0L) {
    cellule <- numeric()
  } else if (is.null(colonna)) {
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
  nella_curva <- !is.na(cellule) & n_pos > 1 & !is.na(fedelta) & conf$esito != "non conforme"
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
  comuni <- integer()
  piastre_comuni <- character()
  punti_per_piastra <- table(piastra[positivi][nella_curva])
  if (length(motivo_globale) == 0L) {
    # Senza alcuna piastra la tabella e' vuota e non ha nomi.
    piastre_comuni <- sort(as.character(names(punti_per_piastra)[punti_per_piastra >= min_positivi]),
                           method = "radix")
    comuni <- which(nella_curva & piastra[positivi] %in% piastre_comuni)
    if (length(comuni) == 0L) {
      motivo_globale <- sprintf(paste0(
        "nessuna piastra ha almeno ctrl.min_positives (%d) controlli positivi ",
        "utilizzabili per la curva"), min_positivi)
    }
  }
  if (length(motivo_globale) == 0L) {
    curve[["aggregato"]] <- valuta("aggregato", comuni)
    for (p in piastre_comuni) {
      curve[[paste0("piastra ", p)]] <- valuta(
        paste0("piastra ", p), comuni[piastra[positivi][comuni] == p])
    }
  }
  aggregato <- curve[["aggregato"]]
  per_piastra <- curve[names(curve) != "aggregato"]
  r2_complessivo <- function(insieme) {
    if (!length(insieme)) return(NA_real_)
    osservate <- unlist(lapply(insieme, function(cv) fedelta[cv$indici]))
    previste <- unlist(lapply(insieme, function(cv) {
      1 / (1 + (cv$x50 / log10(n_pos[cv$indici]))^cv$h)
    }))
    1 - sum((osservate - previste)^2) / sum((osservate - mean(osservate))^2)
  }
  r2_aggregato <- if (is.null(aggregato) || !aggregato$converge) NA_real_ else aggregato$r2
  r2_per_piastra <- r2_complessivo(Filter(function(cv) cv$converge, per_piastra))

  # ---- Il modello preferito: AIC sugli stessi punti -------------------------
  aic_aggregato <- if (is.null(aggregato)) NA_real_ else aggregato$aic
  aic_per_piastra <- if (length(per_piastra)) {
    sum(vapply(per_piastra, function(cv) cv$aic, numeric(1)))
  } else NA_real_
  preferito <- if (is.null(aggregato)) "nessuno" else if (
    !is.na(aic_aggregato) && !is.na(aic_per_piastra) && aic_aggregato < aic_per_piastra
  ) "aggregato" else "per_piastra"
  aic_testo <- function(v) if (is.na(v)) "non calcolabile" else sprintf("%.2f", v)

  # ---- Una soglia per piastra -----------------------------------------------
  voce <- function(valore, origine, motivo) {
    list(valore = valore, stadio = stadio, descrizione_stadio = STADI[[stadio]],
         origine = origine, motivo = motivo)
  }
  aggregato_valido <- !is.null(aggregato) && aggregato$valida
  proprie <- Filter(function(cv) cv$valida, per_piastra)
  valori_propri <- vapply(proprie, function(cv) ceiling(cv$soglia), numeric(1))
  mediana <- if (length(valori_propri)) stats::median(valori_propri) else NA_real_
  usa_aggregato <- preferito == "aggregato" && aggregato_valido
  # Passi 3, 4 e 6: la soglia di chi non ha una curva propria.
  senza_curva <- function(motivo) {
    if (aggregato_valido) {
      voce(ceiling(aggregato$soglia), "aggregata", motivo)
    } else if (!is.na(mediana)) {
      voce(ceiling(mediana), "mediana", sprintf(
        "%s; curva aggregata non valida: %s", motivo, aggregato$motivo_validita))
    } else {
      voce(NULL, "nessuna", motivo)
    }
  }
  piastre <- sort(unique(piastra[!is.na(piastra)]), method = "radix")
  if (senza_soglia) {
    scelta <- "nessuno"
    motivo_scelta <- "qc.min_reads_mode e' none: nessuna soglia di profondita'"
  } else if (length(motivo_globale)) {
    scelta <- "nessuno"
    motivo_scelta <- paste(motivo_globale, collapse = "; ")
  } else if (usa_aggregato) {
    scelta <- "aggregato"
    motivo_scelta <- sprintf("AIC dell'aggregato %s, minore di %s del modello per piastra",
                             aic_testo(aic_aggregato), aic_testo(aic_per_piastra))
  } else if (aggregato_valido || length(proprie)) {
    scelta <- "per_piastra"
    motivo_scelta <- if (preferito == "aggregato") sprintf(
      "AIC dell'aggregato %s, minore di %s del modello per piastra, ma l'aggregato non e' valido: %s",
      aic_testo(aic_aggregato), aic_testo(aic_per_piastra), aggregato$motivo_validita
    ) else if (is.na(aic_per_piastra)) sprintf(paste0(
      "AIC del modello per piastra non calcolabile (una curva non converge): ",
      "l'aggregato (AIC %s) non e' preferito"), aic_testo(aic_aggregato)
    ) else sprintf("AIC del modello per piastra %s, non maggiore di %s dell'aggregato",
                   aic_testo(aic_per_piastra), aic_testo(aic_aggregato))
  } else {
    scelta <- "nessuno"
    motivo_scelta <- paste0("nessuna curva valida: ", paste(vapply(curve, function(cv) {
      sprintf("%s: %s", cv$modello, cv$motivo_validita)
    }, character(1)), collapse = " | "))
  }
  soglie <- list()
  for (p in piastre) {
    cv <- per_piastra[[paste0("piastra ", p)]]
    soglie[[p]] <- if (scelta == "nessuno") {
      voce(NULL, "nessuna", motivo_scelta)
    } else if (scelta == "aggregato") {
      voce(ceiling(aggregato$soglia), "aggregata", motivo_scelta)
    } else if (!is.null(cv) && cv$valida) {
      voce(ceiling(cv$soglia), "propria", "")
    } else if (is.null(cv)) {
      senza_curva(sprintf(
        "%d controlli positivi utilizzabili nella piastra, meno di ctrl.min_positives (%d)",
        if (p %in% names(punti_per_piastra)) punti_per_piastra[[p]] else 0L, min_positivi))
    } else {
      senza_curva(paste("curva della piastra non valida:", cv$motivo_validita))
    }
  }
  senza_piastra <- if (scelta == "nessuno") {
    voce(NULL, "nessuna", motivo_scelta)
  } else if (scelta == "aggregato") {
    voce(ceiling(aggregato$soglia), "aggregata", motivo_scelta)
  } else {
    senza_curva("campione senza piastra")
  }
  applicate <- c(soglie, if (anyNA(piastra)) list("senza piastra" = senza_piastra))
  non_proprie <- Filter(function(s) s$origine %in% c("aggregata", "mediana"), applicate)
  degrada <- length(non_proprie) > 0L
  usano_la_propria <- per_piastra[paste0("piastra ", names(Filter(
    function(s) s$origine == "propria", soglie)))]
  r2_dichiarato <- switch(scelta, aggregato = r2_aggregato,
                          per_piastra = r2_complessivo(usano_la_propria), NA_real_)

  # ---- Misura sui campioni: chi cadrebbe sotto la soglia -------------------
  applicata <- lapply(seq_along(campioni), function(i) {
    if (is.na(piastra[i])) senza_piastra else soglie[[piastra[i]]]
  })
  valore <- vapply(applicata, function(s) if (is.null(s$valore)) NA_real_ else s$valore,
                   numeric(1))
  sotto <- !is.na(valore) & profondita < valore
  scrivi_atomico(
    c(sprintf("accession\tclasse\tpiastra\tletture_%s\tsoglia\tstadio\torigine\tsotto_soglia",
              stadio),
      sprintf("%s\t%s\t%s\t%.0f\t%s\t%s\t%s\t%s", campioni, dati$classe,
              ifelse(is.na(piastra), "", piastra), profondita,
              ifelse(is.na(valore), "", sprintf("%.0f", valore)), stadio,
              vapply(applicata, function(s) s$origine, character(1)),
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
    sprintf("%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%d\t%s\t%s", cv$modello, cv$punti,
            numero(cv$h), numero(cv$x50), numero(cv$k), numero(cv$r2, 4L),
            numero(cv$rmse, 4L), numero(cv$aic),
            if (is.na(cv$soglia)) "" else sprintf("%.0f", ceiling(cv$soglia)),
            if (is.na(cv$punto_medio)) "" else sprintf("%.0f", cv$punto_medio),
            cv$osservazioni_transizione, if (cv$valida) "si" else "no", cv$motivo_validita)
  }, character(1))
  scrivi_atomico(
    c(paste0("modello\tpunti\th\tx50\tk\tr2\trmse\taic\tsoglia\tpunto_medio\t",
             "osservazioni_transizione\tvalida\tmotivo"), unname(righe_curve)),
    file.path(cartella, "curve.tsv")
  )

  valutabili <- conf$esito != "non valutabile"
  conformi <- sum(conf$esito == "conforme")
  frazione <- if (any(valutabili)) conformi / sum(valutabili) else NA_real_
  o_nullo <- function(v, cifre = 6L) if (is.na(v)) NULL else signif(v, cifre)
  scrivi_json(
    list(
      modo = parametri$modo,
      sensibilita = sensibilita,
      modello_curva = "allosteric_sigmoid: f = x^h / (k' + x^h), x = log10(profondita'), k' = x50^h",
      stadio = stadio,
      descrizione_stadio = STADI[[stadio]],
      scelta = scelta,
      motivo_scelta = motivo_scelta,
      modello = list(
        preferito_aic = preferito,
        aic_aggregato = o_nullo(aic_aggregato, 10L),
        aic_per_piastra = o_nullo(aic_per_piastra, 10L),
        punti_comuni = length(comuni),
        piastre_comuni = as.list(piastre_comuni),
        r2_dichiarato = o_nullo(r2_dichiarato),
        r2_aggregato = o_nullo(r2_aggregato),
        r2_per_piastra = o_nullo(r2_per_piastra),
        aggregato_valido = aggregato_valido
      ),
      # La mediana delle soglie proprie, calcolata sui valori interi applicati
      # alle piastre: e' il ripiego del passo 4, e si registra anche se nessuno
      # la usa.
      mediana = if (is.na(mediana)) NULL else list(
        valore_non_arrotondato = mediana, valore = ceiling(mediana),
        piastre = as.list(sub("^piastra ", "", names(proprie))),
        usata = any(vapply(applicate, function(s) s$origine == "mediana", logical(1)))),
      per_piastra = soglie,
      senza_piastra = senza_piastra,
      degradazione = degrada,
      non_proprie = unname(lapply(names(non_proprie), function(n) {
        list(piastra = n, origine = non_proprie[[n]]$origine,
             valore = non_proprie[[n]]$valore, motivo = non_proprie[[n]]$motivo)
      }))
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
