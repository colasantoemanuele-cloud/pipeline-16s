# Oggetto di partenza, campioni e variabili: cio' che si verifica prima di ogni
# calcolo, e cio' che ogni analisi usa per sapere quali campioni entrano.
#
# Richiede uscite.R. Interfaccia:
#   leggi_oggetto(percorso)          l'oggetto phyloseq, o NULL con il motivo
#   conteggi_di(ps)                  matrice varianti x campioni
#   tabella_campioni(ps)             data.frame di testo, una riga per campione
#   e_mancante(x, marcatori)         valori mancanti di una colonna
#   esamina(ps, parametri)           errori, avvisi e disegno dell'analisi
#   tabella_dei_campioni(disegno)    la tabella campioni.tsv
#
# Razionale biologico: le quantita' per campione (indici, abbondanze relative,
# distanze) non dipendono da una variabile e si calcolano su tutti i campioni
# del sottoinsieme; un valore mancante in una variabile toglie il campione dalle
# sole analisi di quella variabile, e ogni esclusione e' registrata con il
# motivo. Un gruppo sotto la dimensione minima esce dai test, non dalle
# descrizioni: e' una dichiarazione, non un arresto.

leggi_oggetto <- function(percorso) {
  letto <- tryCatch(
    list(riuscito = TRUE, valore = readRDS(percorso)),
    error = function(e) list(riuscito = FALSE, motivo = conditionMessage(e))
  )
  if (!letto$riuscito) {
    return(list(oggetto = NULL,
                motivo = paste0("il file non si legge come RDS: ", letto$motivo)))
  }
  ps <- letto$valore
  if (!methods::is(ps, "phyloseq")) {
    return(list(oggetto = NULL,
                motivo = paste0("il file contiene un oggetto di classe ",
                                paste(class(ps), collapse = "/"), ", non phyloseq")))
  }
  assenti <- c(
    if (is.null(phyloseq::access(ps, "otu_table"))) "tabella dei conteggi",
    if (is.null(phyloseq::access(ps, "sam_data"))) "tabella dei campioni",
    if (is.null(phyloseq::access(ps, "tax_table"))) "tassonomia"
  )
  if (length(assenti)) {
    return(list(oggetto = NULL,
                motivo = paste0("all'oggetto manca: ", paste(assenti, collapse = ", "))))
  }
  # La rarefazione e gli indici sono definiti su conteggi: numeri interi, non
  # negativi, tutti presenti. Proporzioni o conteggi trasformati non lo sono.
  conteggi <- conteggi_di(ps)
  if (!is.numeric(conteggi) || anyNA(conteggi)) {
    return(list(oggetto = NULL, motivo = "la tabella dei conteggi ha valori mancanti o non numerici"))
  }
  if (any(conteggi < 0) || any(conteggi != round(conteggi))) {
    return(list(oggetto = NULL, motivo = paste0(
      "la tabella dei conteggi ha valori negativi o non interi: servono i ",
      "conteggi delle letture, non proporzioni o valori trasformati")))
  }
  if (nrow(conteggi) < 2L) {
    return(list(oggetto = NULL, motivo = sprintf(
      "l'oggetto ha %d varianti: ne servono almeno 2", nrow(conteggi))))
  }
  list(oggetto = ps, motivo = NULL)
}

conteggi_di <- function(ps) {
  m <- methods::as(phyloseq::otu_table(ps), "matrix")
  if (!phyloseq::taxa_are_rows(ps)) m <- t(m)
  m
}

tabella_campioni <- function(ps) {
  # La conversione in data.frame riscrive i nomi di colonna non sintattici
  # (spazi, cifre iniziali): quelli veri sono i nomi delle variabili
  # dell'oggetto, nello stesso ordine.
  grezza <- methods::as(phyloseq::sample_data(ps), "data.frame")
  # I valori si confrontano come testo, senza spazi in testa e in coda.
  testo <- lapply(grezza, function(x) trimws(as.character(x)))
  names(testo) <- phyloseq::sample_variables(ps)
  structure(testo, class = "data.frame", row.names = phyloseq::sample_names(ps))
}

e_mancante <- function(x, marcatori) {
  x <- as.character(x)
  is.na(x) | tolower(trimws(x)) %in% marcatori
}

.avviso <- function(codice, dettaglio) list(codice = codice, dettaglio = dettaglio)

.elenco <- function(x, massimo = 10L) {
  x <- as.character(x)
  if (length(x) > massimo) {
    paste0(paste(x[seq_len(massimo)], collapse = ", "), " e altri ", length(x) - massimo)
  } else {
    paste(x, collapse = ", ")
  }
}

# Gruppi di una variabile fra i campioni indicati: valori, mancanti, dimensioni,
# gruppi che entrano nei test.
gruppi_di <- function(valori, minimo) {
  # `valori` e' un vettore con nome (campione -> valore), NA se mancante.
  presenti <- as.character(valori[!is.na(valori)])
  nomi <- unique(presenti)
  nomi <- nomi[ordine_radix(nomi)]
  # Senza alcun valore i tre elenchi sono vuoti, non nulli.
  dimensioni <- stats::setNames(
    vapply(nomi, function(g) sum(presenti == g), integer(1)), nomi)
  piccoli <- nomi[dimensioni < minimo]
  list(
    dimensioni = dimensioni,
    piccoli = piccoli,
    nei_test = nomi[dimensioni >= minimo]
  )
}

esamina <- function(ps, parametri) {
  errori <- list()
  avvisi <- list()
  errore <- function(codice, dettaglio) {
    errori[[length(errori) + 1L]] <<- .avviso(codice, dettaglio)
  }
  avviso <- function(codice, dettaglio) {
    avvisi[[length(avvisi) + 1L]] <<- .avviso(codice, dettaglio)
  }
  design <- parametri$design
  marcatori <- unlist(parametri$valori_mancanti)
  biologiche <- c(design$variable, unlist(design$other_variables))
  tecniche <- as.character(unlist(design$technical_variables))
  sottoinsieme <- design$subset
  if (is.null(sottoinsieme)) sottoinsieme <- list()

  tabella <- tabella_campioni(ps)
  conteggi <- conteggi_di(ps)
  letture <- colSums(conteggi)

  # --- Cio' che rende impossibile ogni calcolo: ci si ferma qui. -------------
  dichiarate <- unique(c(biologiche, tecniche, names(sottoinsieme),
                         colonne_aggiuntive(parametri)))
  for (colonna in dichiarate[!dichiarate %in% names(tabella)]) {
    errore("E-ECO-05", sprintf(
      "Colonna '%s' non trovata. Colonne dell'oggetto: %s.", colonna,
      paste(names(tabella), collapse = ", ")))
  }
  ranghi <- phyloseq::rank_names(ps)
  if (!parametri$comp$rank %in% ranghi) {
    errore("E-ECO-07", sprintf("comp.rank '%s' non trovato. Ranghi dell'oggetto: %s.",
                              parametri$comp$rank, paste(ranghi, collapse = ", ")))
  }
  distanze <- as.character(unlist(parametri$beta$distances))
  albero <- phyloseq::access(ps, "phy_tree")
  if (any(grepl("^unifrac", distanze))) {
    radice <- parametri$beta$unifrac_root
    if (is.null(albero)) {
      errore("E-ECO-06", "L'oggetto non contiene un albero filogenetico.")
    } else if (length(albero$edge.length) != nrow(albero$edge)) {
      errore("E-ECO-06", paste0(
        "L'albero dell'oggetto non ha le lunghezze dei rami: UniFrac le somma, e ",
        "senza non e' definito nemmeno il punto medio."))
    } else if (radice == "existing" && !ape::is.rooted(albero)) {
      errore("E-ECO-06", paste0(
        "L'albero dell'oggetto non ha la radice e beta.unifrac_root e' existing: ",
        "phyloseq ne sceglierebbe una a caso, e la distanza dipenderebbe da quella ",
        "scelta. Con beta.unifrac_root midpoint l'albero si radica al punto medio."))
    }
  }
  if (length(errori)) return(list(errori = errori, avvisi = avvisi, disegno = NULL))

  # --- Sottoinsieme, prima di tutto il resto. --------------------------------
  tenuti <- rep(TRUE, nrow(tabella))
  for (colonna in names(sottoinsieme)) {
    ammessi <- trimws(as.character(unlist(sottoinsieme[[colonna]])))
    assenti <- ammessi[!ammessi %in% tabella[[colonna]]]
    if (length(assenti)) {
      presenti <- unique(tabella[[colonna]][!is.na(tabella[[colonna]])])
      errore("E-ECO-08", sprintf(
        "design.subset: nella colonna '%s' non compare: %s. Valori presenti: %s.",
        colonna, .elenco(assenti), .elenco(presenti[ordine_radix(presenti)])))
    }
    tenuti <- tenuti & tabella[[colonna]] %in% ammessi
  }
  fuori_sottoinsieme <- rownames(tabella)[!tenuti]
  senza_letture <- rownames(tabella)[tenuti & letture == 0]
  if (length(senza_letture)) {
    avviso("E-ECO-15", sprintf("Campioni senza letture: %d (%s).",
                              length(senza_letture), .elenco(senza_letture)))
  }
  campioni <- rownames(tabella)[tenuti & letture > 0]
  campioni <- campioni[ordine_radix(campioni)]
  if (length(campioni) < 2L) {
    errore("E-ECO-08", sprintf(
      "Campioni con letture dopo il sottoinsieme: %d su %d dell'oggetto.",
      length(campioni), nrow(tabella)))
    return(list(errori = errori, avvisi = avvisi, disegno = NULL))
  }
  if (length(errori)) return(list(errori = errori, avvisi = avvisi, disegno = NULL))

  # --- Profondita' di rarefazione. -------------------------------------------
  dichiarata <- parametri$alpha$rarefy_depth
  if (is.null(dichiarata)) {
    profondita <- as.integer(min(letture[campioni]))
    origine <- "minima fra i campioni analizzati"
    sotto <- character()
  } else {
    profondita <- as.integer(dichiarata)
    origine <- "dichiarata (alpha.rarefy_depth)"
    sotto <- campioni[letture[campioni] < profondita]
    if (length(sotto)) {
      avviso("E-ECO-14", sprintf(
        "Campioni sotto la profondita' %d: %d (%s).", profondita, length(sotto),
        .elenco(sotto)))
    }
    if (length(campioni) - length(sotto) < 2L) {
      errore("E-ECO-10", sprintf(
        "Con alpha.rarefy_depth %d restano %d campioni su %d; la profondita' massima e' %d.",
        profondita, length(campioni) - length(sotto), length(campioni),
        as.integer(max(letture[campioni]))))
    }
  }

  # --- Variabili: valori mancanti e gruppi. ----------------------------------
  minimo <- as.integer(parametri$stat$min_group_size)
  variabili <- list()
  for (nome in biologiche) {
    valori <- stats::setNames(tabella[campioni, nome], campioni)
    valori[e_mancante(valori, marcatori)] <- NA
    mancanti <- names(valori)[is.na(valori)]
    if (length(mancanti)) {
      avviso("E-ECO-13", sprintf("Variabile '%s': %d campioni con valore mancante (%s).",
                                nome, length(mancanti), .elenco(mancanti)))
    }
    g <- gruppi_di(valori, minimo)
    if (length(g$piccoli)) {
      avviso("E-ECO-11", sprintf(
        "Variabile '%s': %d gruppi con meno di %d campioni (%s).", nome,
        length(g$piccoli), minimo,
        .elenco(sprintf("%s: %d", g$piccoli, g$dimensioni[g$piccoli]))))
    }
    if (length(g$nei_test) < 2L) {
      avviso("E-ECO-12", sprintf(
        "Variabile '%s': %d gruppi con almeno %d campioni su %d gruppi.", nome,
        length(g$nei_test), minimo, length(g$dimensioni)))
    }
    variabili[[nome]] <- list(nome = nome, valori = valori, mancanti = mancanti,
                              dimensioni = g$dimensioni, piccoli = g$piccoli,
                              nei_test = g$nei_test)
  }
  tecnici <- list()
  for (nome in tecniche) {
    valori <- stats::setNames(tabella[campioni, nome], campioni)
    valori[e_mancante(valori, marcatori)] <- NA
    tecnici[[nome]] <- valori
  }

  list(
    errori = errori,
    avvisi = avvisi,
    disegno = list(
      tabella = tabella,
      nell_oggetto = rownames(tabella),
      fuori_sottoinsieme = fuori_sottoinsieme,
      senza_letture = senza_letture,
      campioni = campioni,
      letture = letture,
      profondita = profondita,
      origine_profondita = origine,
      sotto_profondita = sotto,
      biologiche = biologiche,
      tecniche = tecniche,
      variabili = variabili,
      tecnici = tecnici,
      ranghi = ranghi,
      albero = !is.null(albero),
      varianti = nrow(conteggi)
    )
  )
}

# Le colonne dichiarate da parametri che non stanno nel gruppo design: gli
# strati della PERMANOVA.
colonne_aggiuntive <- function(parametri) {
  as.character(unlist(parametri$stat$permanova_strata))
}

# Lo stato di ogni campione per una variabile.
stato_per_variabile <- function(v, campioni) {
  stato <- rep("analizzato", length(campioni))
  names(stato) <- campioni
  stato[v$mancanti] <- "escluso: valore mancante"
  piccoli <- campioni[!is.na(v$valori[campioni]) & v$valori[campioni] %in% v$piccoli]
  stato[piccoli] <- "escluso dai test: gruppo sotto la dimensione minima"
  stato
}

tabella_dei_campioni <- function(d, piani = list()) {
  tutti <- d$nell_oggetto[ordine_radix(d$nell_oggetto)]
  generale <- rep("analizzato", length(tutti))
  names(generale) <- tutti
  generale[d$fuori_sottoinsieme] <- "escluso: fuori dal sottoinsieme"
  generale[d$senza_letture] <- "escluso: senza letture"
  alfa <- generale
  alfa[d$sotto_profondita] <- "escluso: sotto la profondita' di rarefazione"
  tabella <- data.frame(
    campione = tutti,
    letture = as.integer(d$letture[tutti]),
    stato = generale,
    stato_alfa = alfa,
    stringsAsFactors = FALSE, check.names = FALSE
  )
  for (v in d$variabili) {
    valori <- rep(NA_character_, length(tutti))
    names(valori) <- tutti
    valori[d$campioni] <- v$valori[d$campioni]
    stato <- generale
    stato[d$campioni] <- stato_per_variabile(v, d$campioni)
    tabella[[paste0("valore:", v$nome)]] <- valori
    tabella[[paste0("stato:", v$nome)]] <- stato
    piano <- piani[[v$nome]]
    if (!is.null(piano)) {
      # Nel confronto dell'alfa diversita' i gruppi si contano sui soli
      # campioni rarefatti.
      alfa_v <- alfa
      alfa_v[d$campioni] <- ifelse(alfa[d$campioni] == "analizzato",
                                   piano$alfa$stato[d$campioni], alfa[d$campioni])
      tabella[[paste0("alfa:", v$nome)]] <- alfa_v
      # Nella PERMANOVA entrano i soli campioni con un valore anche nelle
      # variabili tecniche e negli strati.
      stato[d$campioni] <- piano$permanova$stato[d$campioni]
      tabella[[paste0("permanova:", v$nome)]] <- stato
    }
  }
  for (nome in d$tecniche) {
    valori <- rep(NA_character_, length(tutti))
    names(valori) <- tutti
    valori[d$campioni] <- d$tecnici[[nome]][d$campioni]
    tabella[[paste0("valore:", nome)]] <- valori
  }
  tabella
}

# Il disegno in forma di riepilogo: numeri e nomi, senza i vettori per campione.
riepilogo_del_disegno <- function(d) {
  list(
    campioni = list(
      nell_oggetto = length(d$nell_oggetto),
      analizzati = length(d$campioni),
      esclusi = list(
        fuori_dal_sottoinsieme = length(d$fuori_sottoinsieme),
        senza_letture = length(d$senza_letture)
      )
    ),
    varianti = d$varianti,
    alfa = list(
      profondita_di_rarefazione = d$profondita,
      origine = d$origine_profondita,
      campioni_rarefatti = length(d$campioni) - length(d$sotto_profondita),
      esclusi_sotto_la_profondita = I(d$sotto_profondita)
    ),
    variabili = unname(lapply(d$variabili, function(v) list(
      nome = v$nome,
      campioni_con_valore = sum(!is.na(v$valori)),
      campioni_con_valore_mancante = length(v$mancanti),
      gruppi = as.list(v$dimensioni),
      gruppi_esclusi_dai_test = I(v$piccoli),
      gruppi_nei_test = I(v$nei_test)
    ))),
    variabili_tecniche = I(d$tecniche)
  )
}
