# Ordinazione vincolata (supervisionata): RDA e db-RDA, con i test per
# permutazione del modello e dei termini.
#
# La usano R/valida.R (il solo piano, con i suoi avvisi) e R/analisi.R (eco
# run), quando constrained.method e' dichiarato. Scrive, nella cartella
# vincolata/ delle uscite, varianza.tsv (quote di varianza e R2 corretto),
# test.tsv (modello e termini), <modello>_coordinate.tsv e un triplot
# <modello>_NN_<variabile>.png per ogni variabile del disegno con un valore fra
# i campioni del modello.
#
# Richiede uscite.R, campioni.R, alfa.R (fissa_seme), beta.R (clr) e test.R
# (.rango, .fattori). Interfaccia:
#   piano_vincolato(d, parametri)    campioni, dati del modello, gradi di
#                                    liberta', avvisi: non calcola nulla, e
#                                    serve anche alla validazione
#   ordinazione_vincolata(...)       il modello, i test e i punteggi
#   scrivi_vincolata(...)            tabelle e triplot
#
# Razionale biologico: la PCoA cerca gli assi di massima variabilita' dei
# campioni senza guardare le variabili (non vincolata); la RDA e la db-RDA
# (Legendre e Anderson 1999) cercano gli assi della sola variabilita' che le
# variabili dichiarate spiegano, e ne dicono la quota. La PERMANOVA verifica se
# i gruppi differiscono; il metodo vincolato misura quanta variabilita' le
# variabili spiegano insieme, al netto delle condizioni (Condition()), e ne
# mostra la direzione. Con molti livelli (una posizione per campione o quasi)
# il modello si avvicina a quello saturo: spiega molto per costruzione, e la
# quota vincolata non dice piu' nulla. Per questo i gradi di liberta' del
# modello si confrontano con una frazione dichiarata dei campioni.

# Le trasformazioni dei conteggi per la RDA, in parole.
DESCRIZIONE_TRASFORMAZIONE <- c(
  hellinger = "trasformazione di Hellinger (vegan::decostand, method = \"hellinger\"; Legendre e Gallagher 2001): radice delle proporzioni",
  clr = "trasformata log-rapporto centrata dei conteggi piu' beta.clr_pseudocount"
)

# Il piano dell'ordinazione vincolata: quali campioni entrano (un valore in
# ogni variabile e condizione del modello), di che tipo e' ogni termine
# (continuo se e' fra le variabili continue del disegno, altrimenti
# categorico), quali termini escono perche' hanno un solo valore fra i
# campioni, i gradi di liberta' del modello a confronto con
# constrained.max_df_fraction. Restituisce list(eseguibile, campioni, dati,
# variabili, condizioni, nomi, gl, stato, avvisi). Non calcola nulla.
piano_vincolato <- function(d, parametri) {
  vincolata <- facoltativo(parametri, "constrained", list())
  avvisi <- list()
  avviso <- function(codice, dettaglio) {
    avvisi[[length(avvisi) + 1L]] <<- .avviso(codice, dettaglio)
  }
  if (is.null(vincolata$method)) return(NULL)
  marcatori <- unlist(parametri$valori_mancanti)
  variabili <- as.character(unlist(vincolata$variables))
  condizioni <- as.character(unlist(vincolata$conditions))
  stato <- stats::setNames(rep("analizzato", length(d$campioni)), d$campioni)
  valori_di <- function(nome) {
    if (nome %in% names(d$continue)) return(d$continue[[nome]]$valori[d$campioni])
    valori <- stats::setNames(d$tabella[d$campioni, nome], d$campioni)
    valori[e_mancante(valori, marcatori)] <- NA
    valori
  }
  colonne <- lapply(stats::setNames(nm = c(variabili, condizioni)), valori_di)
  # Con gli strati della PERMANOVA le permutazioni sono entro gli strati: un
  # campione senza strato non entra.
  strati <- parametri$stat$permanova_strata
  da_controllare <- colonne
  if (!is.null(strati)) da_controllare[[strati]] <- valori_di(strati)
  for (nome in names(da_controllare)) {
    senza <- names(da_controllare[[nome]])[is.na(da_controllare[[nome]]) & stato == "analizzato"]
    if (length(senza)) {
      stato[senza] <- sprintf("escluso: valore mancante in '%s'", nome)
      avviso("E-ECO-13", sprintf(
        "Ordinazione vincolata: %d campioni senza valore in '%s' (%s).", length(senza), nome,
        .elenco(senza)))
    }
  }
  campioni <- names(stato)[stato == "analizzato"]
  # Un termine con un solo valore fra i campioni del modello non spiega nulla.
  for (nome in names(colonne)) {
    if (length(unique(colonne[[nome]][campioni])) < 2L) {
      avviso("E-ECO-45", sprintf(
        "Ordinazione vincolata: '%s' ha un solo valore fra i %d campioni del modello ed esce dal modello.",
        nome, length(campioni)))
      variabili <- setdiff(variabili, nome)
      condizioni <- setdiff(condizioni, nome)
    }
  }
  eseguibile <- length(variabili) > 0L && length(campioni) >= 3L
  if (!eseguibile) {
    avviso("E-ECO-45", sprintf(
      "Ordinazione vincolata: %s.", if (!length(variabili))
        "nessuna variabile del modello ha almeno due valori fra i campioni" else
        sprintf("i campioni con un valore in ogni termine sono %d, ne servono almeno 3",
                length(campioni))))
  }
  # Dati del modello con nomi di colonna neutri: quelli dell'oggetto possono
  # contenere spazi.
  dati <- data.frame(row.names = campioni)
  nomi <- character()
  for (i in seq_along(variabili)) {
    x <- colonne[[variabili[i]]][campioni]
    dati[[paste0("v", i)]] <- if (variabili[i] %in% names(d$continue)) as.numeric(x) else factor(x)
    nomi[paste0("v", i)] <- variabili[i]
  }
  for (i in seq_along(condizioni)) {
    x <- colonne[[condizioni[i]]][campioni]
    dati[[paste0("c", i)]] <- if (condizioni[i] %in% names(d$continue)) as.numeric(x) else factor(x)
    nomi[paste0("c", i)] <- condizioni[i]
  }
  gl <- NA_integer_
  if (eseguibile) {
    condizionati <- dati[, startsWith(names(dati), "c"), drop = FALSE]
    gl <- .rango(dati) - .rango(condizionati)
    frazione <- vincolata$max_df_fraction
    if (gl > frazione * length(campioni)) {
      avviso("E-ECO-44", sprintf(
        "Ordinazione vincolata: il modello usa %d gradi di liberta' con %d campioni, oltre constrained.max_df_fraction (%s, cioe' %s).",
        gl, length(campioni), format(frazione), format(frazione * length(campioni))))
    }
    if (gl >= length(campioni) - .rango(condizionati)) {
      eseguibile <- FALSE
      avviso("E-ECO-45", sprintf(
        "Ordinazione vincolata: il modello e' saturo (%d gradi di liberta' con %d campioni): non resta variabilita' non vincolata.",
        gl, length(campioni)))
    }
  }
  list(eseguibile = eseguibile, campioni = campioni, dati = dati, variabili = variabili,
       condizioni = condizioni, nomi = nomi, gl = gl, stato = stato, avvisi = avvisi,
       metodo = vincolata$method, trasformazione = vincolata$transform,
       distanze = as.character(unlist(vincolata$distances)))
}

# Il segno di ogni asse fissato sui punteggi dei campioni (l'elemento di modulo
# massimo positivo) e applicato allo stesso modo a centroidi e frecce.
.segni <- function(campioni) {
  vapply(seq_len(ncol(campioni)), function(j) {
    if (campioni[which.max(abs(campioni[, j])), j] < 0) -1 else 1
  }, numeric(1))
}

# Il modello vincolato su una risposta (la matrice trasformata per la RDA, la
# distanza per la db-RDA), con i test per permutazione del modello e dei
# termini (vegan::anova.cca, by = "margin") e i punteggi sui primi due assi.
# Restituisce list(varianza, test, coordinate, centroidi, frecce, assi), o
# nulla se il modello ha meno di due assi con inerzia positiva.
ordinazione_vincolata <- function(piano, risposta, metodo, etichetta, valori_strati,
                                  permutazioni, seme) {
  variabili <- names(piano$nomi)[startsWith(names(piano$nomi), "v")]
  condizioni <- names(piano$nomi)[startsWith(names(piano$nomi), "c")]
  destra <- paste(c(variabili, if (length(condizioni))
    sprintf("Condition(%s)", paste(condizioni, collapse = " + "))), collapse = " + ")
  ambiente <- new.env(parent = globalenv())
  assign("risposta", risposta, envir = ambiente)
  formula <- stats::as.formula(paste("risposta ~", destra), env = ambiente)
  modello <- if (metodo == "rda") {
    vegan::rda(formula, data = piano$dati)
  } else {
    vegan::dbrda(formula, data = piano$dati, add = FALSE, sqrt.dist = FALSE)
  }
  controllo <- if (is.null(valori_strati)) {
    permute::how(nperm = permutazioni)
  } else {
    permute::how(nperm = permutazioni, blocks = factor(valori_strati[piano$campioni]))
  }
  fissa_seme(seme)
  globale <- vegan::anova.cca(modello, permutations = controllo, parallel = 1)
  fissa_seme(seme)
  termini <- vegan::anova.cca(modello, by = "margin", permutations = controllo, parallel = 1)
  usate <- nrow(attr(globale, "F.perm"))
  test <- rbind(
    data.frame(termine = "modello", ruolo = "vincolo", gl = globale$Df[1L],
               varianza = globale[[2L]][1L], f = globale$F[1L], p = globale[["Pr(>F)"]][1L],
               stringsAsFactors = FALSE),
    data.frame(termine = unname(piano$nomi[rownames(termini)[-nrow(termini)]]),
               ruolo = "variabile", gl = termini$Df[-nrow(termini)],
               varianza = termini[[2L]][-nrow(termini)], f = termini$F[-nrow(termini)],
               p = termini[["Pr(>F)"]][-nrow(termini)], stringsAsFactors = FALSE),
    data.frame(termine = "residuo", ruolo = "residuo", gl = globale$Df[2L],
               varianza = globale[[2L]][2L], f = NA_real_, p = NA_real_,
               stringsAsFactors = FALSE))
  test$permutazioni <- usate
  r2 <- vegan::RsquareAdj(modello)
  totale <- modello$tot.chi
  condizionata <- if (is.null(modello$pCCA)) 0 else modello$pCCA$tot.chi
  vincolata <- modello$CCA$tot.chi
  varianza <- data.frame(
    campioni = length(piano$campioni), gl_modello = modello$CCA$qrank,
    inerzia_totale = totale, condizionata = condizionata, vincolata = vincolata,
    non_vincolata = modello$CA$tot.chi,
    immaginaria = if (is.null(modello$CA$imaginary.chi)) 0 else modello$CA$imaginary.chi,
    quota_vincolata = vincolata / totale, r2 = r2$r.squared, r2_corretto = r2$adj.r.squared,
    p_modello = globale[["Pr(>F)"]][1L], stringsAsFactors = FALSE)

  # Punteggi sui primi due assi (scaling 2): i campioni, i centroidi dei
  # livelli delle variabili categoriche e le frecce di quelle continue. Con un
  # solo asse vincolato il secondo e' il primo non vincolato.
  # Servono due assi (vincolati o non vincolati) con inerzia positiva.
  if (sum(vegan::eigenvals(modello) > 0) < 2L) return(NULL)
  coordinate <- vegan::scores(modello, display = "sites", choices = 1:2, scaling = 2)
  segni <- .segni(coordinate)
  coordinate <- sweep(coordinate, 2L, segni, "*")
  dimnames(coordinate) <- list(piano$campioni, colnames(coordinate))
  autovalori <- vegan::eigenvals(modello)
  assi <- sprintf("%s (%.1f%% dell'inerzia)", colnames(coordinate),
                  100 * autovalori[colnames(coordinate)] / totale)
  centroidi <- tryCatch(vegan::scores(modello, display = "cn", choices = 1:2, scaling = 2),
                        error = function(e) NULL)
  if (!is.null(centroidi)) centroidi <- sweep(centroidi, 2L, segni, "*")
  continue <- variabili[vapply(variabili, function(v) is.numeric(piano$dati[[v]]), logical(1))]
  frecce <- NULL
  if (length(continue)) {
    bp <- vegan::scores(modello, display = "bp", choices = 1:2, scaling = 2)
    frecce <- sweep(bp[intersect(continue, rownames(bp)), , drop = FALSE], 2L, segni, "*")
  }
  list(varianza = varianza, test = test, coordinate = coordinate, centroidi = centroidi,
       frecce = frecce, assi = assi, etichetta = etichetta, metodo = metodo)
}

# Il triplot di un modello vincolato colorato per una variabile del disegno: i
# campioni (colore del gruppo), i centroidi dei livelli delle variabili
# categoriche del modello e le frecce di quelle continue. Un centroide e' un
# punto grande bordato di nero, pieno del colore del suo gruppo se il termine
# e' la variabile del grafico (bianco altrimenti); la sua etichetta sta su un
# fondo bianco semitrasparente, accanto al centroide dove non copre le altre
# etichette ne' gli altri centroidi (.posa_etichette), con una linea verso il
# centroide se ha dovuto allontanarsi.
.triplot <- function(esito, piano, v, titolo, percorso) {
  coordinate <- esito$coordinate
  campioni <- intersect(rownames(coordinate), names(v$valori)[!is.na(v$valori)])
  gruppi <- names(v$dimensioni)
  colori <- stats::setNames(colori_gruppi(length(gruppi)), gruppi)
  presenti <- table(factor(v$valori[campioni], levels = gruppi))
  gruppi <- gruppi[presenti > 0]
  voci <- sprintf("%s (n = %d)", gruppi, as.integer(presenti[gruppi]))
  colonne <- if (length(gruppi) > SOGLIA_LEGENDA_COLONNE) 2L else 1L
  apri_png(percorso, larghezza = LARGHEZZA_GRAFICO + (colonne - 1L) * PIXEL_COLONNA_LEGENDA)
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::par(mar = c(4.5, 4.5, 3, 4 + 11 * colonne), xpd = NA, cex = 0.7)
  punti <- rbind(coordinate, esito$centroidi)
  if (!is.null(esito$frecce)) {
    # Le frecce in scala con i campioni, come vegan::ordiArrowMul.
    scala <- 0.8 * max(abs(coordinate)) / max(abs(esito$frecce))
    frecce <- esito$frecce * scala
    punti <- rbind(punti, frecce)
  }
  graphics::plot(punti, type = "n", asp = 1, xlab = esito$assi[1L], ylab = esito$assi[2L],
                 main = sprintf("%s (%d campioni)", titolo, length(campioni)))
  graphics::abline(h = 0, v = 0, col = "grey85", xpd = FALSE)
  graphics::points(coordinate[campioni, 1L], coordinate[campioni, 2L], pch = 19,
                   col = colori[v$valori[campioni]])
  if (!is.null(esito$centroidi)) {
    centroidi <- esito$centroidi
    riempimento <- .colori_centroidi(rownames(centroidi), piano, v, colori)
    graphics::points(centroidi, pch = 21, cex = CEX_CENTROIDE, lwd = 1.5, col = "black",
                     bg = riempimento)
    testi <- .etichette_centroidi(rownames(centroidi), piano)
    posa <- .posa_etichette(centroidi[, 1L], centroidi[, 2L], testi, CEX_ETICHETTA,
                            CEX_CENTROIDE)
    lontane <- posa$spostata
    if (any(lontane)) {
      # La linea va dal centroide al punto piu' vicino del riquadro.
      vx <- pmin(pmax(centroidi[lontane, 1L], posa$x[lontane] - posa$larghezza[lontane] / 2),
                 posa$x[lontane] + posa$larghezza[lontane] / 2)
      vy <- pmin(pmax(centroidi[lontane, 2L], posa$y[lontane] - posa$altezza[lontane] / 2),
                 posa$y[lontane] + posa$altezza[lontane] / 2)
      graphics::segments(centroidi[lontane, 1L], centroidi[lontane, 2L], vx, vy,
                         col = "grey30", lwd = 0.8)
    }
    graphics::rect(posa$x - posa$larghezza / 2, posa$y - posa$altezza / 2,
                   posa$x + posa$larghezza / 2, posa$y + posa$altezza / 2,
                   col = grDevices::adjustcolor("white", alpha.f = 0.8), border = "grey55",
                   lwd = 0.6)
    graphics::text(posa$x, posa$y, labels = testi, cex = CEX_ETICHETTA)
  }
  if (!is.null(esito$frecce)) {
    graphics::arrows(0, 0, frecce[, 1L], frecce[, 2L], length = 0.08, lwd = 1.5)
    graphics::text(frecce * 1.08, labels = unname(piano$nomi[rownames(frecce)]), cex = 0.85)
  }
  limiti <- graphics::par("usr")
  graphics::legend(x = limiti[2L] + 0.03 * (limiti[2L] - limiti[1L]), y = limiti[4L],
                   legend = voci, col = colori[gruppi], pch = 19, bty = "n", ncol = colonne,
                   title = v$nome)
}

# Le dimensioni dei centroidi e delle loro etichette nel triplot, rispetto al
# testo del grafico.
CEX_CENTROIDE <- 2.2
CEX_ETICHETTA <- 0.85

# Il colore di riempimento dei centroidi: quello del gruppo se il termine e' la
# variabile v del grafico, bianco per i livelli degli altri termini.
.colori_centroidi <- function(nomi, piano, v, colori) {
  vapply(nomi, function(n) {
    colonna <- names(piano$nomi)[startsWith(n, names(piano$nomi))]
    colonna <- colonna[which.max(nchar(colonna))]
    if (!length(colonna) || piano$nomi[[colonna]] != v$nome) return("white")
    livello <- substring(n, nchar(colonna) + 1L)
    if (livello %in% names(colori)) colori[[livello]] else "white"
  }, character(1), USE.NAMES = FALSE)
}

# Le posizioni delle etichette dei centroidi, in coordinate del grafico (che
# deve essere aperto: le misure del testo dipendono dal dispositivo). Le
# etichette si posano una alla volta, nell'ordine dei centroidi: per ciascuna
# si provano otto posizioni attorno al suo centroide (sopra, sotto, a destra,
# a sinistra, poi in diagonale), prima accostate e poi via via piu' lontane,
# e si tiene la prima il cui riquadro non copre ne' le etichette gia' posate
# ne' alcun centroide. Se nessuna e' libera resta sopra il centroide. Nulla e'
# casuale: gli stessi punti danno le stesse posizioni. Restituisce
# data.frame(x, y, larghezza, altezza, spostata): spostata e' vero se
# l'etichetta non e' accostata al suo centroide (serve una linea).
.posa_etichette <- function(x, y, testi, cex, cex_punto) {
  n <- length(testi)
  larghezza <- graphics::strwidth(testi, cex = cex) + graphics::strwidth("M", cex = cex) * 0.8
  altezza <- rep(graphics::strheight("Mg", cex = cex) * 1.7, n)
  # Il raggio di un centroide, in unita' dell'asse y (con asp = 1 uguale in x).
  raggio <- graphics::strheight("M", cex = cex_punto) * 0.55
  margine <- altezza[1L] * 0.15
  copre <- function(cx, cy, w, h, posate) {
    sui_centroidi <- any(abs(cx - x) < w / 2 + raggio & abs(cy - y) < h / 2 + raggio)
    sulle_etichette <- any(abs(cx - posate$x) < (w + posate$larghezza) / 2 + margine &
                             abs(cy - posate$y) < (h + posate$altezza) / 2 + margine)
    sui_centroidi || sulle_etichette
  }
  # Le direzioni, nell'ordine in cui si provano: sopra, sotto, destra,
  # sinistra e le quattro diagonali.
  dx <- c(0, 0, 1, -1, 1, -1, 1, -1)
  dy <- c(1, -1, 0, 0, 1, 1, -1, -1)
  posate <- data.frame(x = numeric(), y = numeric(), larghezza = numeric(), altezza = numeric(),
                       spostata = logical())
  for (i in seq_len(n)) {
    w <- larghezza[i]
    h <- altezza[i]
    scelta <- NULL
    for (anello in 0:6) {
      distanza <- anello * h
      for (k in seq_along(dx)) {
        cx <- x[i] + dx[k] * (raggio + w / 2 + margine + distanza)
        cy <- y[i] + dy[k] * (raggio + h / 2 + margine + distanza)
        if (!copre(cx, cy, w, h, posate)) {
          scelta <- c(cx, cy, anello > 0)
          break
        }
      }
      if (!is.null(scelta)) break
    }
    if (is.null(scelta)) scelta <- c(x[i], y[i] + raggio + h / 2 + margine, 0)
    posate[i, ] <- list(scelta[1L], scelta[2L], w, h, scelta[3L] > 0)
  }
  posate
}

# I nomi dei centroidi: il livello, preceduto dalla variabile («variabile:
# livello») solo se il modello ha piu' di un termine categorico. vegan li chiama
# con il nome neutro della colonna seguito dal livello (v1NOD1).
.etichette_centroidi <- function(nomi, piano) {
  categorici <- sum(vapply(piano$dati, is.factor, logical(1)))
  vapply(nomi, function(n) {
    colonna <- names(piano$nomi)[startsWith(n, names(piano$nomi))]
    colonna <- colonna[which.max(nchar(colonna))]
    if (!length(colonna)) return(n)
    livello <- substring(n, nchar(colonna) + 1L)
    if (categorici > 1L) paste0(piano$nomi[[colonna]], ": ", livello) else livello
  }, character(1))
}

# Calcola e scrive l'ordinazione vincolata: per la RDA un modello sui conteggi
# trasformati, per la db-RDA uno per distanza; un modello senza variabilita'
# (distanze tutte nulle, meno di due assi) si salta con l'avviso E-ECO-45.
# Restituisce list(file, riepilogo, avvisi).
scrivi_vincolata <- function(piano, ps, matrici, d, parametri, seme, cartella) {
  dir.create(cartella, showWarnings = FALSE)
  prodotti <- character()
  riepilogo <- list()
  permutazioni <- as.integer(parametri$stat$permanova_permutations)
  strati <- parametri$stat$permanova_strata
  valori_strati <- NULL
  if (!is.null(strati)) {
    valori_strati <- stats::setNames(d$tabella[d$campioni, strati], d$campioni)
  }
  modelli <- if (piano$metodo == "rda") {
    conteggi <- conteggi_di(ps)[, piano$campioni, drop = FALSE]
    conteggi <- conteggi[rowSums(conteggi) > 0, , drop = FALSE]
    risposta <- if (piano$trasformazione == "hellinger") {
      vegan::decostand(t(conteggi), method = "hellinger")
    } else {
      clr(conteggi, parametri$beta$clr_pseudocount)
    }
    list(list(base = paste0("rda_", piano$trasformazione), risposta = risposta,
              dati = DESCRIZIONE_TRASFORMAZIONE[[piano$trasformazione]]))
  } else {
    lapply(piano$distanze, function(distanza) {
      list(base = paste0("dbrda_", distanza),
           risposta = stats::as.dist(matrici[[distanza]][piano$campioni, piano$campioni]),
           dati = DESCRIZIONE_DISTANZA[[distanza]], distanza = distanza)
    })
  }
  varianze <- list()
  tutti_i_test <- list()
  avvisi <- list()
  modello_testo <- sprintf("modello: %s%s", paste(unname(piano$nomi[startsWith(names(piano$nomi), "v")]),
                                                  collapse = " + "),
                           if (length(piano$condizioni))
                             sprintf(" + Condition(%s)", paste(piano$condizioni, collapse = " + ")) else "")
  for (m in modelli) {
    # Con distanze tutte nulle non c'e' variabilita' da scomporre.
    if (max(as.matrix(m$risposta)) <= 0 && !is.null(m$distanza)) {
      avvisi[[length(avvisi) + 1L]] <- .avviso("E-ECO-45", sprintf(
        "Ordinazione vincolata, distanza %s: tutte le distanze fra i campioni del modello sono nulle.",
        m$distanza))
      next
    }
    esito <- ordinazione_vincolata(piano, m$risposta, piano$metodo, m$base, valori_strati,
                                   permutazioni, seme)
    if (is.null(esito)) {
      avvisi[[length(avvisi) + 1L]] <- .avviso("E-ECO-45", sprintf(
        "Ordinazione vincolata, %s: meno di due assi con inerzia positiva, nessun grafico in due dimensioni.",
        m$base))
      next
    }
    varianze[[m$base]] <- cbind(data.frame(modello = m$base, stringsAsFactors = FALSE), esito$varianza)
    tutti_i_test[[m$base]] <- cbind(data.frame(modello = m$base, stringsAsFactors = FALSE), esito$test)
    # La palette estesa per interpolazione nei triplot di queste coordinate: i
    # gruppi di ogni variabile presenti fra i campioni del modello.
    note <- unlist(lapply(seq_along(d$variabili), function(i) {
      v <- d$variabili[[i]]
      presenti <- unique(v$valori[intersect(names(v$valori)[!is.na(v$valori)], piano$campioni)])
      nota <- nota_palette(length(presenti))
      if (length(nota)) sprintf("grafico %s_%s.png: %s", m$base, nome_file(i, v$nome), nota)
    }))
    scrivi_tsv(affianca(list(campione = rownames(esito$coordinate)), esito$coordinate, reale),
               file.path(cartella, paste0(m$base, "_coordinate.tsv")), intestazione = c(
                 sprintf("%s: punteggi dei campioni sui primi due assi (vegan::scores, display = \"sites\", scaling = 2).",
                         if (piano$metodo == "rda") "RDA (vegan::rda)" else "db-RDA (vegan::dbrda)"),
                 sprintf("dati: %s", m$dati), modello_testo,
                 "con un solo asse vincolato, il secondo e' il primo non vincolato; il segno di ogni asse e' fissato (l'elemento di modulo massimo e' positivo)",
                 note))
    prodotti <- c(prodotti, paste0("vincolata/", m$base, "_coordinate.tsv"))
    indice <- 0L
    for (v in d$variabili) {
      indice <- indice + 1L
      # Una variabile senza valore fra i campioni del modello non colora nulla:
      # nessun triplot.
      con_valore <- names(v$valori)[!is.na(v$valori)]
      if (!length(v$dimensioni) || !length(intersect(con_valore, piano$campioni))) next
      nome <- paste0(m$base, "_", nome_file(indice, v$nome), ".png")
      .triplot(esito, piano, v, sprintf("%s, %s", if (piano$metodo == "rda") "RDA" else "db-RDA",
                                        if (is.null(m$distanza)) piano$trasformazione else
                                          nome_distanza(m$distanza)),
               file.path(cartella, nome))
      prodotti <- c(prodotti, paste0("vincolata/", nome))
    }
    riepilogo[[length(riepilogo) + 1L]] <- list(
      modello = m$base, campioni = length(piano$campioni),
      quota_vincolata = esito$varianza$quota_vincolata, r2_corretto = esito$varianza$r2_corretto,
      p_modello = esito$varianza$p_modello)
  }
  if (!length(varianze)) return(list(file = character(), riepilogo = riepilogo, avvisi = avvisi))
  varianza <- do.call(rbind, varianze)
  for (colonna in c("inerzia_totale", "condizionata", "vincolata", "non_vincolata", "immaginaria",
                    "quota_vincolata", "r2", "r2_corretto", "p_modello")) {
    varianza[[colonna]] <- reale(varianza[[colonna]])
  }
  scrivi_tsv(varianza, file.path(cartella, "varianza.tsv"), intestazione = c(
    "Ordinazione vincolata: scomposizione della variabilita' (inerzia) e R2.",
    modello_testo,
    "condizionata: tolta dalle condizioni (Condition()); vincolata: spiegata dalle variabili; non_vincolata: il resto; immaginaria: la parte negativa di una distanza non euclidea nella db-RDA",
    "quota_vincolata: vincolata / inerzia_totale; r2 e r2_corretto: vegan::RsquareAdj (R2 corretto per il numero di termini, Peres-Neto et al. 2006)",
    "gl_modello: gradi di liberta' delle variabili; con molti livelli il modello si avvicina a quello saturo e la quota vincolata cresce per costruzione (constrained.max_df_fraction)",
    "p_modello: test per permutazione del modello (vegan::anova.cca); i termini sono in test.tsv"))
  test <- do.call(rbind, tutti_i_test)
  for (colonna in c("varianza", "f", "p")) test[[colonna]] <- reale(test[[colonna]])
  scrivi_tsv(test, file.path(cartella, "test.tsv"), intestazione = c(
    "Ordinazione vincolata: test per permutazione (vegan::anova.cca).",
    modello_testo,
    "modello: il test di tutte le variabili insieme; variabile: ogni variabile al netto delle altre e delle condizioni (by = \"margin\"); residuo: la variabilita' non vincolata",
    sprintf("permutazioni richieste: %d (stat.permanova_permutations), dal seme %d; strati: %s",
            permutazioni, as.integer(seme), if (is.null(strati)) "nessuno" else strati),
    "la PERMANOVA (test/permanova.tsv) resta il test fra gruppi; questo misura quanta variabilita' le variabili spiegano insieme"))
  list(file = c("vincolata/varianza.tsv", "vincolata/test.tsv", prodotti),
       riepilogo = riepilogo, avvisi = avvisi)
}
