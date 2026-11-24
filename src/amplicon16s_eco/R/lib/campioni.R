# Oggetto di partenza, campioni e variabili: cio' che si verifica prima di ogni
# calcolo, e cio' che ogni analisi usa per sapere quali campioni entrano.
#
# La usano R/valida.R (eco validate, e la validazione che precede eco run),
# R/analisi.R (eco run) e R/descrivi.R (eco init). Legge l'oggetto dalla copia;
# non scrive file: campioni.tsv lo scrive R/analisi.R con la tabella costruita
# qui.
#
# Richiede uscite.R. Interfaccia:
#   leggi_oggetto(percorso)          l'oggetto phyloseq, o NULL con il motivo
#   conteggi_di(ps)                  matrice varianti x campioni
#   tabella_campioni(ps)             data.frame di testo, una riga per campione
#   e_mancante(x, marcatori)         valori mancanti di una colonna
#   deriva_raggruppamenti(...)       i gruppi derivati (design.groupings) come
#                                    nuove colonne della tabella dei campioni
#   filtra_abbondanza(conteggi, soglia, ambito)   il filtro delle varianti rare
#   esamina(ps, parametri)           errori, avvisi e disegno dell'analisi
#   tabella_dei_campioni(disegno)    la tabella campioni.tsv
#   scrivi_filtro(d, cartella)       la tabella delle varianti filtrate
#
# Razionale biologico: le quantita' per campione (indici, abbondanze relative,
# distanze) non dipendono da una variabile e si calcolano su tutti i campioni
# del sottoinsieme; un valore mancante in una variabile toglie il campione dalle
# sole analisi di quella variabile, e ogni esclusione e' registrata con il
# motivo. Un gruppo sotto la dimensione minima esce dai test, non dalle
# descrizioni: e' una dichiarazione, non un arresto.

# Quando una variabile sembra continua: tutti i valori presenti sono numeri, i
# valori distinti sono almeno MINIMO_VALORI_CONTINUA e almeno la quota indicata
# dei campioni con valore.
MINIMO_VALORI_CONTINUA <- 10L
QUOTA_VALORI_CONTINUA <- 0.8

# L'oggetto di partenza letto dal file, o il motivo per cui non si puo' usare.
# Restituisce list(oggetto, motivo), con uno dei due nullo: non si ferma, cosi'
# la validazione puo' riportare il rifiuto insieme agli altri e l'analisi
# dichiararlo con il suo codice. Controlla, in quest'ordine: che il file si
# legga come RDS, che contenga un oggetto phyloseq, che abbia conteggi,
# campioni e tassonomia, che i conteggi siano numeri interi non negativi senza
# mancanti, che le varianti siano almeno due.
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

# La matrice dei conteggi con le varianti sulle righe e i campioni sulle
# colonne, qualunque sia l'orientamento nell'oggetto.
conteggi_di <- function(ps) {
  m <- methods::as(phyloseq::otu_table(ps), "matrix")
  if (!phyloseq::taxa_are_rows(ps)) m <- t(m)
  m
}

# La tabella dei campioni dell'oggetto come data.frame di solo testo, con i
# nomi di colonna veri e i campioni come nomi di riga.
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

# Quali valori di una colonna sono mancanti: NA, oppure uno dei marcatori
# dichiarati, confrontati senza distinguere maiuscole e spazi ai lati.
e_mancante <- function(x, marcatori) {
  x <- as.character(x)
  is.na(x) | tolower(trimws(x)) %in% marcatori
}

# I gruppi derivati (design.groupings) come nuove colonne della tabella dei
# campioni. Per ognuno: il nome non deve essere quello di una colonna
# dell'oggetto; from deve esserlo; con groups ogni valore elencato deve
# comparire nella colonna, e ogni valore va al suo gruppo; con pattern il primo
# gruppo di cattura dell'espressione (PCRE, come re di Python nella
# validazione del file) e' il nome del gruppo. Un valore mancante, che non e'
# elencato o che l'espressione non riconosce resta senza gruppo (NA), e il
# motivo e' registrato per campione. Restituisce list(tabella, derivate,
# errori): la tabella con le nuove colonne, per ogni gruppo derivato la regola
# e i motivi, e i rifiuti E-ECO-40.
deriva_raggruppamenti <- function(tabella, raggruppamenti, marcatori) {
  derivate <- list()
  errori <- list()
  colonne <- names(tabella)
  for (nome in names(raggruppamenti)) {
    g <- raggruppamenti[[nome]]
    da <- g[["from"]]
    if (nome %in% colonne) {
      errori[[length(errori) + 1L]] <- .avviso("E-ECO-40", sprintf(
        "design.groupings.%s: il nome e' quello di una colonna dell'oggetto; scegline un altro.",
        nome))
      next
    }
    if (!da %in% colonne) {
      errori[[length(errori) + 1L]] <- .avviso("E-ECO-40", sprintf(
        "design.groupings.%s: from '%s' non e' una colonna dell'oggetto. Colonne dell'oggetto: %s.",
        nome, da, paste(colonne, collapse = ", ")))
      next
    }
    valori <- tabella[[da]]
    mancante <- e_mancante(valori, marcatori)
    if (!is.null(g$groups)) {
      gruppi <- g$groups
      elencati <- unlist(lapply(gruppi, function(x) trimws(as.character(unlist(x)))))
      mappa <- stats::setNames(rep(names(gruppi), lengths(lapply(gruppi, unlist))), elencati)
      presenti <- unique(valori[!mancante])
      assenti <- setdiff(elencati, presenti)
      if (length(assenti)) {
        errori[[length(errori) + 1L]] <- .avviso("E-ECO-40", sprintf(
          "design.groupings.%s: nella colonna '%s' non compare: %s. Valori presenti: %s.",
          nome, da, .elenco(assenti), .elenco(presenti[ordine_radix(presenti)])))
        next
      }
      derivati <- unname(mappa[valori])
      regola <- sprintf("groups: %s", paste(sprintf("%s = %s", names(gruppi), vapply(
        gruppi, function(x) paste(unlist(x), collapse = ", "), "")), collapse = "; "))
    } else {
      testo <- ifelse(is.na(valori), "", valori)
      catture <- regmatches(testo, regexec(g$pattern, testo, perl = TRUE))
      derivati <- vapply(catture, function(m) {
        if (length(m) >= 2L && !is.na(m[2L]) && nzchar(m[2L])) m[2L] else NA_character_
      }, character(1))
      regola <- sprintf("pattern: %s (primo gruppo di cattura)", g$pattern)
    }
    derivati[mancante] <- NA_character_
    motivi <- ifelse(mancante, sprintf("escluso: valore mancante in '%s'", da), sprintf(
      "escluso: il valore di '%s' non rientra in alcun gruppo di '%s'", da, nome))
    motivi[!is.na(derivati)] <- NA_character_
    tabella[[nome]] <- derivati
    derivate[[nome]] <- list(nome = nome, da = da, regola = regola,
                             motivi = stats::setNames(motivi, rownames(tabella)))
  }
  list(tabella = tabella, derivate = derivate, errori = errori)
}

# Il motivo per cui i campioni senza valore escono dalle analisi di una
# variabile: per un gruppo derivato quello registrato campione per campione,
# altrimenti il valore mancante.
motivo_mancante <- function(v, campioni) {
  if (is.null(v$motivi)) return(rep("escluso: valore mancante", length(campioni)))
  unname(v$motivi[campioni])
}

# Un errore o un avviso come coppia codice e dettaglio, la forma in cui
# finiscono nella validazione e nel riepilogo.
.avviso <- function(codice, dettaglio) list(codice = codice, dettaglio = dettaglio)

# Un elenco di nomi per un messaggio: i primi massimo, e quanti altri restano.
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

# Il filtro delle varianti rare per abbondanza relativa (filter.*), sui
# conteggi dei campioni analizzati (varianti x campioni). soglia e' in
# percentuale; ambito dice su che cosa si misura: in ogni campione (sample: si
# azzerano nel campione le varianti sotto la soglia delle sue letture), sul
# totale (total: si tolgono da tutti i campioni le varianti la cui somma e'
# sotto la soglia delle letture totali, Bokulich et al. 2013) o sulla media
# (mean: si tolgono le varianti la cui abbondanza relativa media fra i campioni
# e' sotto la soglia). Il confronto e' fra interi e prodotti esatti
# (conteggio * 100 < soglia * letture), non fra quozienti arrotondati: una
# variante esattamente alla soglia resta. Restituisce list(conteggi, tolte,
# campioni_toccati, tolte_per_campione): la matrice filtrata con le stesse
# righe e colonne, le letture tolte per variante, in quanti campioni ciascuna e'
# stata tolta e le letture tolte per campione.
filtra_abbondanza <- function(conteggi, soglia, ambito) {
  letture <- colSums(conteggi)
  togli <- switch(
    ambito,
    sample = conteggi * 100 < soglia * rep(letture, each = nrow(conteggi)),
    total = matrix(rowSums(conteggi) * 100 < soglia * sum(letture),
                   nrow(conteggi), ncol(conteggi)),
    mean = matrix(rowMeans(sweep(conteggi, 2L, letture, "/")) * 100 < soglia,
                  nrow(conteggi), ncol(conteggi)),
    stop("ambito del filtro non realizzato: ", ambito)
  )
  togli <- togli & conteggi > 0
  tolti <- conteggi * togli
  filtrati <- conteggi - tolti
  dimnames(filtrati) <- dimnames(conteggi)
  list(conteggi = filtrati, tolte = rowSums(tolti), campioni_toccati = rowSums(togli),
       tolte_per_campione = colSums(tolti))
}

# L'esame dell'oggetto rispetto alla configurazione, prima di ogni calcolo.
# Restituisce list(errori, avvisi, disegno): il disegno e' nullo se c'e' un
# errore che rende impossibile proseguire, altrimenti porta cio' che ogni
# analisi usa (campioni analizzati, letture, profondita' di rarefazione con la
# sua origine, variabili con valori e gruppi, variabili tecniche, continue,
# soggetto; con il filtro, l'oggetto filtrato e cio' che il filtro ha tolto).
# Passi, nell'ordine: colonne dichiarate, rango e albero (gli errori di questo
# passo si raccolgono tutti, poi ci si ferma); sottoinsieme e campioni senza
# letture; filtro di abbondanza e campioni che svuota; singleton per Chao1;
# profondita' di rarefazione (solo con alpha.method rarefy); valori mancanti e gruppi
# di ogni variabile; variabili continue, che devono essere numeri decimali;
# soggetto del disegno a misure ripetute. Non calcola alcuna analisi: serve
# tale e quale alla validazione e all'esecuzione, che cosi' giudicano allo
# stesso modo.
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
  metodo_alfa <- parametri$alpha$method
  # Il gruppo del filtro c'e' solo se il filtro e' dichiarato.
  del_filtro <- facoltativo(parametri, "filter", list())
  marcatori <- unlist(parametri$valori_mancanti)
  biologiche <- c(design$variable, unlist(design$other_variables))
  tecniche <- as.character(unlist(design$technical_variables))
  sottoinsieme <- design$subset
  if (is.null(sottoinsieme)) sottoinsieme <- list()
  # Facoltativi: assenti dai parametri se non dichiarati.
  continue_dichiarate <- as.character(unlist(design$continuous_variables))
  colonna_soggetto <- design$subject

  tabella <- tabella_campioni(ps)
  conteggi <- conteggi_di(ps)
  letture <- colSums(conteggi)
  # Facoltativi: assenti dai parametri se non dichiarati.
  raggruppamenti <- design$groupings
  if (is.null(raggruppamenti)) raggruppamenti <- list()
  colori <- as.character(unlist(parametri$ord$color_by))
  per_composizione <- as.character(unlist(parametri$comp$group_by))
  per_pannelli <- unique(as.character(unlist(c(parametri$comp$facet_wrap,
                                               parametri$comp$facet_grid$rows,
                                               parametri$comp$facet_grid$cols))))
  del_vincolo <- facoltativo(parametri, "constrained", list())
  del_modello <- as.character(unlist(c(del_vincolo$variables, del_vincolo$conditions)))
  del_grafico <- facoltativo(parametri, "plot", list())

  # --- Gruppi derivati: nuove colonne, prima di ogni altro controllo. --------
  derivate <- list()
  if (length(raggruppamenti)) {
    esito <- deriva_raggruppamenti(tabella, raggruppamenti, marcatori)
    tabella <- esito$tabella
    derivate <- esito$derivate
    errori <- c(errori, esito$errori)
  }

  # --- Cio' che rende impossibile ogni calcolo: ci si ferma qui. -------------
  dichiarate <- unique(c(biologiche, tecniche, names(sottoinsieme),
                         colonne_aggiuntive(parametri), continue_dichiarate,
                         colonna_soggetto, colori, per_composizione, del_modello,
                         per_pannelli))
  # Un gruppo derivato respinto non e' una colonna mancante: il suo rifiuto
  # c'e' gia'.
  respinti <- setdiff(names(raggruppamenti), names(derivate))
  for (colonna in dichiarate[!dichiarate %in% c(names(tabella), respinti)]) {
    errore("E-ECO-05", sprintf(
      "Colonna '%s' non trovata. Colonne dell'oggetto: %s.", colonna,
      paste(names(tabella), collapse = ", ")))
  }
  ranghi <- phyloseq::rank_names(ps)
  ranghi_composizione <- as.character(unlist(parametri$comp$rank))
  for (rango in setdiff(ranghi_composizione, ranghi)) {
    errore("E-ECO-07", sprintf("comp.rank '%s' non trovato. Ranghi dell'oggetto: %s.",
                              rango, paste(ranghi, collapse = ", ")))
  }
  livello <- parametri$run$analysis_rank
  if (!is.null(livello)) {
    if (!livello %in% ranghi) {
      errore("E-ECO-07", sprintf("run.analysis_rank '%s' non trovato. Ranghi dell'oggetto: %s.",
                                livello, paste(ranghi, collapse = ", ")))
    } else if (length(ranghi) < 2L) {
      errore("E-ECO-07", sprintf(
        "run.analysis_rank '%s': la tassonomia dell'oggetto ha un solo rango, e le varianti non si agglomerano.",
        livello))
    } else {
      # Dopo l'agglomerazione i ranghi piu' fini non hanno valori: un grafico a
      # uno di quei ranghi sarebbe tutto «non assegnato».
      fini <- ranghi_composizione[ranghi_composizione %in% ranghi &
                                    match(ranghi_composizione, ranghi) > match(livello, ranghi)]
      for (rango in fini) {
        errore("E-ECO-07", sprintf(
          "comp.rank '%s' e' piu' fine di run.analysis_rank '%s': dopo l'agglomerazione i taxa non hanno quel rango. Scegli per i grafici '%s' o un rango piu' generale.",
          rango, livello, livello))
      }
    }
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
  palette <- del_grafico$palette
  if (!is.null(palette) && is.na(origine_palette(palette))) {
    errore("E-ECO-41", sprintf("plot.palette '%s' non e' il nome di una palette conosciuta.",
                               paste(unlist(palette), collapse = ", ")))
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

  # --- Livello dei calcoli: le varianti agglomerate al rango. ----------------
  # Senza rango e senza filtro l'oggetto resta quello letto, e ogni analisi lo
  # usa come sempre: ps_analisi e' nullo. tax_glom con NArm = FALSE tiene le
  # varianti senza nome al rango, separate secondo il loro lignaggio ai ranghi
  # superiori, e pota l'albero alle varianti che rappresentano ogni taxon.
  ps_analisi <- NULL
  filtro <- NULL
  svuotati <- character()
  letture_analisi <- letture
  ps_base <- ps
  agglomerato <- NULL
  if (!is.null(livello)) {
    ps_base <- phyloseq::tax_glom(phyloseq::prune_samples(campioni, ps), taxrank = livello,
                                  NArm = FALSE)
    ps_analisi <- ps_base
    conteggi <- conteggi_di(ps_base)
    letture_analisi <- colSums(conteggi)
    agglomerato <- list(rango = livello, varianti = nrow(conteggi_di(ps)),
                        taxa = phyloseq::ntaxa(ps_base),
                        albero = !is.null(phyloseq::access(ps_base, "phy_tree")))
  }

  # --- Filtro delle varianti rare, prima di ogni analisi. --------------------
  soglia <- del_filtro$min_relative_abundance
  if (!is.null(soglia)) {
    ambito <- del_filtro$scope
    f <- filtra_abbondanza(conteggi[, campioni, drop = FALSE], soglia, ambito)
    svuotati <- campioni[colSums(f$conteggi) == 0]
    if (length(svuotati)) {
      avviso("E-ECO-37", sprintf(
        "Filtro allo %s%% (%s): %d campioni senza letture dopo il filtro (%s).",
        format(soglia), ambito, length(svuotati), .elenco(svuotati)))
    }
    rimasti <- setdiff(campioni, svuotati)
    if (length(rimasti) < 2L) {
      errore("E-ECO-38", sprintf(
        "Con filter.min_relative_abundance %s (scope %s) restano %d campioni con letture su %d.",
        format(soglia), ambito, length(rimasti), length(campioni)))
      return(list(errori = errori, avvisi = avvisi, disegno = NULL))
    }
    # L'oggetto delle analisi: i soli campioni rimasti, i conteggi filtrati, le
    # sole varianti con almeno una lettura (l'albero si pota con loro).
    filtrati <- f$conteggi[, rimasti, drop = FALSE]
    presenti <- rownames(filtrati)[rowSums(filtrati) > 0]
    ps_analisi <- phyloseq::prune_samples(rimasti, ps_base)
    phyloseq::otu_table(ps_analisi) <- phyloseq::otu_table(
      filtrati, taxa_are_rows = TRUE)
    ps_analisi <- phyloseq::prune_taxa(presenti, ps_analisi)
    filtro <- list(soglia = soglia, ambito = ambito, campioni = campioni,
                   prima = rowSums(conteggi[, campioni, drop = FALSE]),
                   tolte = f$tolte, campioni_toccati = f$campioni_toccati,
                   tolte_per_campione = f$tolte_per_campione,
                   varianti_prima = sum(rowSums(conteggi[, campioni, drop = FALSE]) > 0),
                   varianti_dopo = length(presenti))
    campioni <- rimasti
    conteggi <- conteggi_di(ps_analisi)
    letture_analisi <- colSums(conteggi)
  }

  # --- Profondita' di rarefazione. -------------------------------------------
  dichiarata <- parametri$alpha$rarefy_depth
  if (identical(metodo_alfa, "none")) {
    # Nessuna rarefazione: l'alfa diversita' usa tutti i campioni analizzati.
    profondita <- NA_integer_
    origine <- "nessuna rarefazione (alpha.method none)"
    sotto <- character()
  } else if (is.null(dichiarata)) {
    profondita <- as.integer(min(letture_analisi[campioni]))
    origine <- "minima fra i campioni analizzati"
    sotto <- character()
  } else {
    profondita <- as.integer(dichiarata)
    origine <- "dichiarata (alpha.rarefy_depth)"
    sotto <- campioni[letture_analisi[campioni] < profondita]
    if (length(sotto)) {
      avviso("E-ECO-14", sprintf(
        "Campioni sotto la profondita' %d: %d (%s).", profondita, length(sotto),
        .elenco(sotto)))
    }
    if (length(campioni) - length(sotto) < 2L) {
      errore("E-ECO-10", sprintf(
        "Con alpha.rarefy_depth %d restano %d campioni su %d; la profondita' massima e' %d.",
        profondita, length(campioni) - length(sotto), length(campioni),
        as.integer(max(letture_analisi[campioni]))))
    }
  }

  # --- Perplessita' del t-SNE: il vincolo di Rtsne sui campioni analizzati. ---
  # Con meno di quattro campioni nessuna ordinazione si fa (avviso E-ECO-25).
  perplessita <- parametri$ord$tsne_perplexity
  if ("tsne" %in% unlist(parametri$ord$methods) && !is.null(perplessita) &&
      length(campioni) >= 4L && 3 * perplessita > length(campioni) - 1L) {
    errore("E-ECO-42", sprintf(
      "ord.tsne_perplexity %s con %d campioni analizzati: 3 * %s > %d; la perplessita' massima e' %s.",
      format(perplessita), length(campioni), format(perplessita), length(campioni) - 1L,
      format((length(campioni) - 1L) / 3)))
    return(list(errori = errori, avvisi = avvisi, disegno = NULL))
  }

  # --- Singleton dei conteggi di partenza, per Chao1. -------------------------
  # Si contano sui conteggi che entrano nell'alfa diversita', prima della
  # rarefazione: senza singleton Chao1 non ha nulla da cui stimare la ricchezza
  # non osservata (alfa.R).
  singleton <- sum(conteggi[, campioni, drop = FALSE] == 1)
  if (singleton == 0) {
    avviso("E-ECO-36", sprintf(
      "Conteggi di partenza dei %d campioni analizzati: nessuna variante con una sola lettura in un campione; Chao1 coincide con la ricchezza osservata sui conteggi non rarefatti.",
      length(campioni)))
  }

  # --- Variabili: valori mancanti e gruppi. ----------------------------------
  minimo <- as.integer(parametri$stat$min_group_size)
  variabili <- list()
  for (nome in biologiche) {
    valori <- stats::setNames(tabella[campioni, nome], campioni)
    valori[e_mancante(valori, marcatori)] <- NA
    mancanti <- names(valori)[is.na(valori)]
    if (length(mancanti)) {
      avviso("E-ECO-13", sprintf(
        if (nome %in% names(derivate)) "Variabile '%s': %d campioni senza gruppo (%s)." else
          "Variabile '%s': %d campioni con valore mancante (%s).",
        nome, length(mancanti), .elenco(mancanti)))
    }
    g <- gruppi_di(valori, minimo)
    if (length(g$piccoli)) {
      avviso("E-ECO-11", sprintf(
        "Variabile '%s': %d gruppi con meno di %d campioni (%s).", nome,
        length(g$piccoli), minimo,
        .elenco(sprintf("%s: %d", g$piccoli, g$dimensioni[g$piccoli]))))
    }
    # Una variabile numerica con quasi tutti i valori diversi e' con ogni
    # probabilita' una misura continua: resta trattata come categorica, come
    # dichiarato, ma va detto, perche' i suoi "gruppi" sono singoli campioni.
    presenti <- valori[!is.na(valori)]
    distinti <- length(unique(presenti))
    if (length(presenti) > 0L && distinti >= MINIMO_VALORI_CONTINUA &&
        distinti >= QUOTA_VALORI_CONTINUA * length(presenti) &&
        !anyNA(suppressWarnings(as.numeric(presenti)))) {
      avviso("E-ECO-29", sprintf(
        "Variabile '%s': %d valori distinti su %d campioni con valore, tutti numerici.",
        nome, distinti, length(presenti)))
    }
    if (length(g$nei_test) < 2L) {
      avviso("E-ECO-12", sprintf(
        "Variabile '%s': %d gruppi con almeno %d campioni su %d gruppi.", nome,
        length(g$nei_test), minimo, length(g$dimensioni)))
    }
    variabili[[nome]] <- list(nome = nome, valori = valori, mancanti = mancanti,
                              dimensioni = g$dimensioni, piccoli = g$piccoli,
                              nei_test = g$nei_test,
                              motivi = if (nome %in% names(derivate)) derivate[[nome]]$motivi)
  }
  # Le variabili dei soli grafici (ord.color_by): valori e gruppi, nessun test.
  descrittiva <- function(nome) {
    valori <- stats::setNames(tabella[campioni, nome], campioni)
    valori[e_mancante(valori, marcatori)] <- NA
    list(nome = nome, valori = valori, mancanti = names(valori)[is.na(valori)],
         dimensioni = gruppi_di(valori, 1L)$dimensioni)
  }
  colore <- lapply(stats::setNames(nm = colori), descrittiva)
  # Le variabili delle barre della composizione (comp.group_by); senza, quelle
  # del disegno.
  composizione <- if (length(per_composizione)) {
    lapply(stats::setNames(nm = per_composizione), descrittiva)
  } else {
    variabili
  }
  # --- Palette estesa per interpolazione (avviso E-ECO-46). -----------------
  # I grafici con piu' colori di quelli della palette: le variabili del disegno
  # (alfa diversita', ordinazioni, triplot), le colonne di ord.color_by e i
  # taxa della composizione (comp.top_n, o meno se il rango ha meno lignaggi;
  # «altri» e' in grigio, fuori dalla palette).
  palette_scelta <- if (is.null(palette)) PALETTE_GRAFICI else as.character(unlist(palette))
  massimo <- colori_della_palette(palette_scelta)
  palette_estesa <- list()
  if (!is.na(massimo)) {
    oggetto <- if (is.null(ps_analisi)) ps else ps_analisi
    tassonomia <- methods::as(phyloseq::tax_table(oggetto), "matrix")
    lignaggi <- vapply(ranghi_composizione, function(rango) {
      fino_a <- match(rango, colnames(tassonomia))
      if (is.na(fino_a)) return(0L)
      length(unique(apply(tassonomia[, seq_len(fino_a), drop = FALSE], 1L, paste,
                          collapse = ";")))
    }, integer(1))
    taxa_mostrati <- min(as.integer(parametri$comp$top_n), max(lignaggi, 0L))
    richiesti <- c(
      stats::setNames(vapply(variabili, function(v) length(v$dimensioni), integer(1)),
                      sprintf("variabile '%s'", names(variabili))),
      stats::setNames(vapply(colore, function(v) length(v$dimensioni), integer(1)),
                      sprintf("ord.color_by '%s'", names(colore))),
      stats::setNames(taxa_mostrati, "taxa della composizione (comp.top_n)"))
    oltre <- richiesti[richiesti > massimo]
    if (length(oltre)) {
      palette_estesa <- list(palette = nome_palette(palette_scelta), colori = massimo,
                             grafici = unname(Map(function(cosa, n) list(grafici = cosa, colori = n),
                                                  names(oltre), unname(oltre))))
      avviso("E-ECO-46", sprintf(
        "La palette %s ha %d colori; e' estesa per interpolazione nei grafici di: %s.",
        nome_palette(palette_scelta), massimo,
        paste(sprintf("%s (%d colori)", names(oltre), oltre), collapse = "; ")))
    }
  }
  tecnici <- list()
  for (nome in tecniche) {
    valori <- stats::setNames(tabella[campioni, nome], campioni)
    valori[e_mancante(valori, marcatori)] <- NA
    tecnici[[nome]] <- valori
  }

  # --- Variabili continue: i valori devono essere numeri. --------------------
  continue <- list()
  for (nome in continue_dichiarate) {
    grezzi <- stats::setNames(as.character(tabella[campioni, nome]), campioni)
    grezzi[e_mancante(grezzi, marcatori)] <- NA
    numeri <- suppressWarnings(as.numeric(grezzi))
    # Un numero in notazione decimale: as.numeric leggerebbe anche "0x10".
    decimali <- grepl("^[+-]?([0-9]+[.]?[0-9]*|[.][0-9]+)([eE][+-]?[0-9]+)?$", trimws(grezzi))
    non_numerici <- names(grezzi)[!is.na(grezzi) & (!is.finite(numeri) | !decimali)]
    if (length(non_numerici)) {
      errore("E-ECO-30", sprintf(
        "Variabile continua '%s': %d campioni con un valore che non e' un numero (%s).",
        nome, length(non_numerici),
        .elenco(sprintf("%s: %s", non_numerici, grezzi[non_numerici]))))
      next
    }
    mancanti <- names(grezzi)[is.na(grezzi)]
    if (length(mancanti)) {
      avviso("E-ECO-13", sprintf(
        "Variabile continua '%s': %d campioni con valore mancante (%s).",
        nome, length(mancanti), .elenco(mancanti)))
    }
    continue[[nome]] <- list(nome = nome, valori = stats::setNames(numeri, campioni),
                             mancanti = mancanti)
  }
  # --- Soggetto del disegno a misure ripetute. --------------------------------
  soggetto <- NULL
  if (!is.null(colonna_soggetto)) {
    valori <- stats::setNames(as.character(tabella[campioni, colonna_soggetto]), campioni)
    valori[e_mancante(valori, marcatori)] <- NA
    soggetto <- list(nome = colonna_soggetto, valori = valori)
  }

  list(
    errori = errori,
    avvisi = avvisi,
    disegno = list(
      tabella = tabella,
      nell_oggetto = rownames(tabella),
      fuori_sottoinsieme = fuori_sottoinsieme,
      senza_letture = senza_letture,
      svuotati = svuotati,
      campioni = campioni,
      letture = letture,
      letture_analisi = letture_analisi,
      ps = ps_analisi,
      filtro = filtro,
      agglomerato = agglomerato,
      metodo_alfa = if (is.null(metodo_alfa)) "rarefy" else metodo_alfa,
      ripetizioni = if (is.null(parametri$alpha$rarefy_iterations)) 1L else
        as.integer(parametri$alpha$rarefy_iterations),
      profondita = profondita,
      origine_profondita = origine,
      sotto_profondita = sotto,
      singleton = singleton,
      biologiche = biologiche,
      tecniche = tecniche,
      variabili = variabili,
      derivate = derivate,
      colore = colore,
      palette_estesa = palette_estesa,
      composizione = composizione,
      pannelli_wrap = if (!is.null(parametri$comp$facet_wrap)) descrittiva(parametri$comp$facet_wrap),
      pannelli_griglia = if (!is.null(parametri$comp$facet_grid)) list(
        x = variabili[[1L]], righe = descrittiva(parametri$comp$facet_grid$rows),
        colonne = descrittiva(parametri$comp$facet_grid$cols)),
      tecnici = tecnici,
      continue = continue,
      soggetto = soggetto,
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
  stato[v$mancanti] <- motivo_mancante(v, v$mancanti)
  piccoli <- campioni[!is.na(v$valori[campioni]) & v$valori[campioni] %in% v$piccoli]
  stato[piccoli] <- "escluso dai test: gruppo sotto la dimensione minima"
  stato
}

# La tabella campioni.tsv: una riga per campione dell'oggetto, anche se
# escluso, con le letture, lo stato generale, lo stato nell'alfa diversita' e,
# per ogni variabile, il valore e lo stato nelle sue analisi. d e' il disegno
# di esamina; piani, se dati, sono i piani dei test fra gruppi, da cui vengono
# le colonne dello stato nel confronto dell'alfa e nella PERMANOVA; vincolato,
# se dato, e' il piano dell'ordinazione vincolata (colonna vincolata). Le colonne
# delle variabili continue e del soggetto ci sono solo se dichiarate.
tabella_dei_campioni <- function(d, piani = list(), vincolato = NULL) {
  tutti <- d$nell_oggetto[ordine_radix(d$nell_oggetto)]
  generale <- rep("analizzato", length(tutti))
  names(generale) <- tutti
  generale[d$fuori_sottoinsieme] <- "escluso: fuori dal sottoinsieme"
  generale[d$senza_letture] <- "escluso: senza letture"
  generale[d$svuotati] <- "escluso: senza letture dopo il filtro di abbondanza"
  alfa <- generale
  alfa[d$sotto_profondita] <- "escluso: sotto la profondita' di rarefazione"
  tabella <- data.frame(
    campione = tutti,
    letture = as.integer(d$letture[tutti]),
    stato = generale,
    stato_alfa = alfa,
    stringsAsFactors = FALSE, check.names = FALSE
  )
  # Con il filtro, le letture che ha tolto a ogni campione analizzato.
  if (!is.null(d$filtro)) {
    tolte <- rep(NA_integer_, length(tutti))
    names(tolte) <- tutti
    tolte[d$filtro$campioni] <- as.integer(d$filtro$tolte_per_campione[d$filtro$campioni])
    tabella <- cbind(tabella[, 1:2], letture_tolte_dal_filtro = tolte,
                     tabella[, -(1:2)], stringsAsFactors = FALSE)
  }
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
  # I gruppi derivati che non sono variabili del disegno: il gruppo di ogni
  # campione, NA se ne e' senza.
  for (nome in setdiff(names(d$derivate), names(d$variabili))) {
    valori <- rep(NA_character_, length(tutti))
    names(valori) <- tutti
    valori[d$campioni] <- d$tabella[d$campioni, nome]
    tabella[[paste0("valore:", nome)]] <- valori
  }
  # Le colonne delle variabili continue e del soggetto ci sono solo se dichiarate.
  for (v in d$continue) {
    valori <- rep(NA_character_, length(tutti))
    names(valori) <- tutti
    valori[d$campioni] <- as.character(d$tabella[d$campioni, v$nome])
    valori[d$campioni][is.na(v$valori[d$campioni])] <- NA
    tabella[[paste0("valore:", v$nome)]] <- valori
  }
  if (!is.null(d$soggetto)) {
    valori <- rep(NA_character_, length(tutti))
    names(valori) <- tutti
    valori[d$campioni] <- d$soggetto$valori[d$campioni]
    tabella[[paste0("soggetto:", d$soggetto$nome)]] <- valori
  }
  if (!is.null(vincolato)) {
    stato <- generale
    stato[d$campioni] <- vincolato$stato[d$campioni]
    tabella[["vincolata"]] <- stato
  }
  tabella
}

# Scrive nella cartella del filtro la tabella delle varianti a cui il filtro di
# abbondanza ha tolto letture: le letture prima e quelle tolte, in quanti
# campioni, e se la variante e' uscita da tutti. Restituisce il nome del file,
# relativo alla cartella di uscita.
scrivi_filtro <- function(d, cartella) {
  dir.create(cartella, showWarnings = FALSE)
  f <- d$filtro
  toccate <- names(f$tolte)[f$tolte > 0]
  toccate <- toccate[ordine_radix(toccate)]
  scrivi_tsv(
    data.frame(
      variante = toccate,
      letture_prima = as.integer(f$prima[toccate]),
      letture_tolte = as.integer(f$tolte[toccate]),
      campioni_con_letture_tolte = as.integer(f$campioni_toccati[toccate]),
      tolta_da_tutti_i_campioni = ifelse(f$tolte[toccate] == f$prima[toccate], "si", "no"),
      stringsAsFactors = FALSE),
    file.path(cartella, "varianti_filtrate.tsv"),
    intestazione = c(
      sprintf("Varianti a cui il filtro di abbondanza relativa ha tolto letture: soglia %s%% (filter.min_relative_abundance), ambito %s (filter.scope).",
              format(f$soglia), f$ambito),
      switch(f$ambito,
             sample = "sample: in ogni campione si azzerano le varianti con conteggio * 100 < soglia * letture del campione",
             total = "total: si tolgono da tutti i campioni le varianti con letture totali * 100 < soglia * letture totali dei campioni analizzati (Bokulich et al. 2013)",
             mean = "mean: si tolgono le varianti la cui abbondanza relativa media fra i campioni, in percentuale, e' sotto la soglia"),
      sprintf("campioni su cui si misura: %d (dopo design.subset, con letture); varianti con letture prima del filtro: %d, dopo: %d",
              length(f$campioni), f$varianti_prima, f$varianti_dopo),
      sprintf("letture tolte: %d su %d (%s%%); le letture tolte per campione sono in campioni.tsv (letture_tolte_dal_filtro)",
              as.integer(sum(f$tolte)), as.integer(sum(f$prima)),
              formatC(100 * sum(f$tolte) / sum(f$prima), format = "f", digits = 4)),
      "letture_prima: letture della variante nei campioni analizzati prima del filtro; tolta_da_tutti_i_campioni: si se il filtro le ha tolto ogni lettura",
      sprintf("campioni rimasti senza letture, esclusi da tutte le analisi (avviso E-ECO-37): %d", length(d$svuotati))
    ))
  "filtro/varianti_filtrate.tsv"
}

# Il disegno in forma di riepilogo: numeri e nomi, senza i vettori per campione.
riepilogo_del_disegno <- function(d) {
  riepilogo <- list(
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
  # Le voci del filtro, del metodo dell'alfa diversita', delle variabili
  # continue e del soggetto ci sono solo se dichiarate: senza, il riepilogo e'
  # quello di sempre.
  if (!is.null(d$filtro)) {
    riepilogo$campioni$esclusi$svuotati_dal_filtro <- length(d$svuotati)
    riepilogo$filtro <- list(
      soglia_percentuale = d$filtro$soglia, ambito = d$filtro$ambito,
      varianti_prima = d$filtro$varianti_prima, varianti_dopo = d$filtro$varianti_dopo,
      varianti_toccate = sum(d$filtro$tolte > 0),
      letture_tolte = sum(d$filtro$tolte),
      letture_prima = sum(d$filtro$prima))
  }
  if (!is.null(d$agglomerato)) {
    riepilogo$livello_dei_calcoli <- list(
      rango = d$agglomerato$rango, varianti = d$agglomerato$varianti,
      taxa = d$agglomerato$taxa)
  }
  if (length(d$derivate)) {
    # I gruppi ottenuti da ogni regola, fra i campioni analizzati: validate li
    # stampa perche' la regola si controlli prima del calcolo.
    riepilogo$gruppi_derivati <- unname(lapply(d$derivate, function(g) {
      valori <- d$tabella[d$campioni, g$nome]
      list(nome = g$nome, da = g$da, regola = g$regola,
           gruppi = as.list(gruppi_di(stats::setNames(valori, d$campioni), 1L)$dimensioni),
           campioni_senza_gruppo = sum(is.na(valori)))
    }))
  }
  if (length(d$palette_estesa)) {
    riepilogo$palette <- list(
      nome = d$palette_estesa$palette, colori = d$palette_estesa$colori,
      estesa_per_interpolazione = d$palette_estesa$grafici,
      nota = AVVISO_PALETTE)
  }
  if (identical(d$metodo_alfa, "none")) {
    riepilogo$alfa <- list(metodo = "none", campioni = length(d$campioni))
  } else if (d$ripetizioni > 1L) {
    riepilogo$alfa$rarefazioni <- d$ripetizioni
  }
  if (length(d$continue)) {
    riepilogo$variabili_continue <- unname(lapply(d$continue, function(v) list(
      nome = v$nome,
      campioni_con_valore = sum(!is.na(v$valori)),
      campioni_con_valore_mancante = length(v$mancanti),
      valori_distinti = length(unique(v$valori[!is.na(v$valori)]))
    )))
  }
  if (!is.null(d$soggetto)) {
    valori <- d$soggetto$valori
    riepilogo$soggetto <- list(
      colonna = d$soggetto$nome,
      soggetti = length(unique(valori[!is.na(valori)])),
      campioni_senza_soggetto = sum(is.na(valori))
    )
  }
  riepilogo
}
