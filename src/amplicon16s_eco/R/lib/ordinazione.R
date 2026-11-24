# Ordinazione: PCoA, NMDS e t-SNE sulle matrici di distanza.
#
# La usa R/analisi.R (comando eco run). Scrive, nella cartella ordinazione/
# delle uscite, per ogni metodo e distanza le coordinate, per la PCoA gli
# autovalori, un grafico per variabile e un grafico descrittivo, senza test,
# per ogni colonna di ord.color_by.
#
# Richiede uscite.R e alfa.R (fissa_seme). Interfaccia:
#   pcoa_di(matrice, correzione)     coordinate, autovalori, frazioni
#   nmds_di(matrice, trymax, seme)   coordinate, stress, convergenza
#   tsne_di(matrice, perplessita, iterazioni, theta, seme)   coordinate, costo
#   sottotitolo_tsne(perplessita, theta, iterazioni, seme)   il sottotitolo
#                                    dei grafici del t-SNE, con i parametri
#   scrivi_ordinazioni(...)          tabelle e grafici per variabile
#
# Razionale sistemistico: le distanze non euclidee (Bray-Curtis, Jaccard) danno
# alla PCoA autovalori negativi; la correzione e' una scelta dichiarata, e la
# frazione di varianza si riporta sugli autovalori che ne risultano. La NMDS
# riceve la matrice senza attributi: metaMDS non applica trasformazioni ne' la
# riscalatura in unita' di semicambiamento, che scatterebbe da sola per una
# matrice prodotta da vegdist. Il segno di ogni asse e' arbitrario e viene
# fissato, perche' le coordinate non cambino fra macchine.
#
# Razionale biologico del t-SNE (van der Maaten e Hinton 2008): conserva i
# vicini di ogni campione, non le distanze. Le distanze fra gruppi lontani e
# la dimensione dei gruppi nella figura non hanno significato: e' una figura
# descrittiva, e la differenza fra gruppi si giudica con la PERMANOVA.

SOGLIA_STRESS <- 0.2

# I nomi dei metodi nei titoli dei grafici.
NOMI_ORDINAZIONE <- c(pcoa = "PCoA", nmds = "NMDS", tsne = "t-SNE")

# Le coordinate con il segno di ogni asse fissato: l'elemento di modulo massimo
# della colonna e' positivo.
.fissa_segno <- function(coordinate) {
  for (j in seq_len(ncol(coordinate))) {
    massimo <- which.max(abs(coordinate[, j]))
    if (coordinate[massimo, j] < 0) coordinate[, j] <- -coordinate[, j]
  }
  coordinate
}

# L'analisi delle coordinate principali di una matrice di distanza. correzione
# e' quella dichiarata per gli autovalori negativi. Restituisce
# list(coordinate, autovalori, correzione, applicata, negativi, frazione):
# applicata dice se la correzione e' servita (ape la applica solo con
# autovalori negativi), e la frazione di varianza e' sugli autovalori corretti
# quando lo e', altrimenti su tutti.
pcoa_di <- function(matrice, correzione) {
  esito <- ape::pcoa(stats::as.dist(matrice), correction = correzione)
  valori <- esito$values
  negativi <- sum(valori$Eigenvalues < 0)
  applicata <- correzione != "none" && !is.null(esito$vectors.cor)
  if (applicata) {
    coordinate <- esito$vectors.cor
    corretti <- valori$Corr_eig
    frazione <- valori$Rel_corr_eig
  } else {
    coordinate <- esito$vectors
    corretti <- rep(NA_real_, nrow(valori))
    frazione <- valori$Relative_eig
  }
  coordinate <- .fissa_segno(coordinate)
  rownames(coordinate) <- rownames(matrice)
  list(
    coordinate = coordinate,
    autovalori = data.frame(
      asse = seq_len(nrow(valori)), autovalore = valori$Eigenvalues,
      autovalore_corretto = corretti, frazione_di_varianza = frazione,
      stringsAsFactors = FALSE),
    correzione = correzione, applicata = applicata, negativi = negativi,
    frazione = frazione
  )
}

# Lo scaling multidimensionale non metrico in due dimensioni di una matrice di
# distanza, con trymax avvii casuali dal seme indicato. Restituisce
# list(coordinate, stress, quasi_nullo, convergente, avvii): quasi_nullo
# raccoglie l'avviso di stress quasi nullo, che non deve arrivare a chi chiama
# come avviso imprevisto.
nmds_di <- function(matrice, trymax, seme) {
  fissa_seme(seme)
  # Con pochi campioni, o con gruppi compatti e ben separati che si riducono
  # a pochi punti, metaMDS trova una configurazione senza stress e lo segnala:
  # e' un esito atteso, dichiarato con un codice proprio.
  quasi_nullo <- FALSE
  esito <- withCallingHandlers(vegan::metaMDS(
    stats::as.dist(matrice), k = 2, try = trymax, trymax = trymax,
    engine = "monoMDS", autotransform = FALSE, wascores = FALSE, expand = FALSE,
    trace = 0, plot = FALSE,
    # A monoMDS: regressione monotona globale, stress di tipo 1 (quello a cui
    # si riferisce la soglia 0,2), legami deboli. A postMDS: nessuna
    # riscalatura in unita' di semicambiamento.
    model = "global", stress = 1, weakties = TRUE, scaling = TRUE, halfchange = FALSE
  ), warning = function(w) {
    if (grepl("stress is (nearly) zero", conditionMessage(w), fixed = TRUE)) {
      quasi_nullo <<- TRUE
      invokeRestart("muffleWarning")
    }
  })
  coordinate <- .fissa_segno(unclass(esito$points)[, 1:2, drop = FALSE])
  attributes(coordinate) <- list(dim = dim(coordinate),
                                 dimnames = list(rownames(matrice), c("nmds_1", "nmds_2")))
  list(coordinate = coordinate, stress = esito$stress, quasi_nullo = quasi_nullo,
       # `converged` e' il numero di avvii che hanno ritrovato la soluzione
       # migliore (un valore logico nelle versioni precedenti di vegan).
       convergente = as.numeric(esito$converged) > 0, avvii = as.integer(trymax))
}

# Il t-SNE in due dimensioni di una matrice di distanza (Rtsne::Rtsne, con
# is_distance = TRUE e pca = FALSE: le distanze sono quelle dichiarate, senza
# riduzione preliminare). La disposizione iniziale (Y_init) si estrae dal seme
# con la stessa distribuzione che Rtsne usa da se' (normale, deviazione 1e-4):
# il risultato non dipende dal generatore interno. Con Y_init Rtsne
# toglierebbe da se' la fase di esagerazione iniziale (stop_lying_iter e
# mom_switch_iter a 0, pensati per una disposizione gia' ottimizzata): qui
# restano ai valori standard (250), come senza Y_init. Un solo processo
# (num_threads = 1): l'ordine delle somme non cambia fra esecuzioni.
# Restituisce list(coordinate, costo, iterazioni).
tsne_di <- function(matrice, perplessita, iterazioni, theta, seme) {
  n <- nrow(matrice)
  fissa_seme(seme)
  iniziale <- matrix(stats::rnorm(n * 2L, sd = 1e-4), nrow = n, ncol = 2L)
  esito <- Rtsne::Rtsne(
    unname(matrice), dims = 2L, is_distance = TRUE, pca = FALSE,
    perplexity = perplessita, theta = theta, max_iter = as.integer(iterazioni),
    Y_init = iniziale, stop_lying_iter = 250L, mom_switch_iter = 250L,
    num_threads = 1L, verbose = FALSE)
  coordinate <- esito$Y
  dimnames(coordinate) <- list(rownames(matrice), c("tsne_1", "tsne_2"))
  list(coordinate = coordinate, costo = utils::tail(esito$itercosts, 1L),
       iterazioni = as.integer(iterazioni))
}

# Il sottotitolo dei grafici del t-SNE: i parametri con cui rifare la figura.
sottotitolo_tsne <- function(perplessita, theta, iterazioni, seme) {
  sprintf("perplessità %s, theta %s, %d iterazioni, seme %d (Rtsne)",
          format(perplessita), format(theta), as.integer(iterazioni), as.integer(seme))
}

# Il grafico dei campioni sui primi due assi di un'ordinazione, colorati per
# gruppo della variabile v; i campioni senza valore non compaiono. La legenda
# sta a destra, in una o due colonne, oppure in un pannello sotto il grafico,
# secondo il numero dei gruppi: l'area dei punti resta la stessa. Un
# sottotitolo, se c'e', sta sotto il titolo: per il t-SNE i parametri con cui
# rifare la figura.
.grafico_ordinazione <- function(coordinate, v, titolo, etichette, percorso,
                                 sottotitolo = NULL) {
  con_valore <- names(v$valori)[!is.na(v$valori)]
  gruppi <- names(v$dimensioni)
  colori <- stats::setNames(colori_gruppi(length(gruppi)), gruppi)
  voci <- sprintf("%s (n = %d)", gruppi, v$dimensioni[gruppi])
  # Con il sottotitolo il margine alto ha una riga in piu' e il titolo sale.
  alto <- if (is.null(sottotitolo)) 3 else 4
  punti <- function() {
    principale <- sprintf("%s (%d campioni)", titolo, length(con_valore))
    graphics::plot(
      coordinate[con_valore, 1L], coordinate[con_valore, 2L],
      col = colori[v$valori[con_valore]], pch = 19, asp = 1,
      xlab = etichette[1L], ylab = etichette[2L],
      main = if (is.null(sottotitolo)) principale else NULL
    )
    if (!is.null(sottotitolo)) {
      graphics::title(main = principale, line = 2)
      graphics::mtext(sottotitolo, side = 3, line = 0.6, cex = 0.7)
    }
  }
  if (length(gruppi) <= SOGLIA_LEGENDA_PANNELLO) {
    # Legenda a destra, in una o due colonne: l'immagine si allarga della
    # colonna in piu', e il grafico resta largo uguale.
    colonne <- if (length(gruppi) > SOGLIA_LEGENDA_COLONNE) 2L else 1L
    apri_png(percorso, larghezza = LARGHEZZA_GRAFICO + (colonne - 1L) * PIXEL_COLONNA_LEGENDA)
    on.exit(grDevices::dev.off(), add = TRUE)
    graphics::par(mar = c(4.5, 4.5, alto, 4 + 11 * colonne), xpd = NA, cex = 0.7)
    punti()
    limiti <- graphics::par("usr")
    graphics::legend(
      x = limiti[2L] + 0.03 * (limiti[2L] - limiti[1L]), y = limiti[4L],
      legend = voci, col = colori, pch = 19, bty = "n", ncol = colonne, title = v$nome
    )
  } else {
    # Legenda in un pannello sotto il grafico: tante colonne quante ne stanno
    # nella larghezza, e l'immagine si allunga con il numero di righe. La
    # larghezza di una colonna e' quella delle voci tipiche (nove su dieci): i
    # pochi nomi piu' lunghi si accorciano, tenendo il numero dei campioni,
    # invece di ridurre tutta la legenda a una colonna sola.
    larghezza <- min(max(nchar(voci)), max(20L, as.integer(ceiling(
      stats::quantile(nchar(voci), 0.9, names = FALSE)))))
    lunghe <- nchar(voci) > larghezza
    if (any(lunghe)) {
      coda <- sprintf(" (n = %d)", v$dimensioni[gruppi])
      voci[lunghe] <- paste0(substr(gruppi[lunghe], 1L, pmax(1L, larghezza - nchar(coda[lunghe]) - 3L)),
                             "...", coda[lunghe])
    }
    colonne <- max(1L, min(8L, 150L %/% (max(nchar(voci)) + 6L)))
    righe <- ceiling(length(gruppi) / colonne)
    altezza_legenda <- 140L + righe * PIXEL_RIGA_LEGENDA
    apri_png(percorso, altezza = ALTEZZA_GRAFICO + altezza_legenda)
    on.exit(grDevices::dev.off(), add = TRUE)
    graphics::layout(matrix(1:2, ncol = 1L), heights = c(ALTEZZA_GRAFICO, altezza_legenda))
    graphics::par(mar = c(4.5, 4.5, alto, 2), cex = 0.7)
    punti()
    graphics::par(mar = c(0.2, 1, 0.2, 1), cex = 0.7)
    graphics::plot.new()
    graphics::legend("top", legend = voci, col = colori, pch = 19, bty = "n",
                     ncol = colonne, title = v$nome)
  }
}

# Calcola e scrive le ordinazioni dichiarate: per ogni distanza e per ogni
# metodo le coordinate, per la PCoA anche gli autovalori, e un grafico per
# variabile. Restituisce list(file, riepilogo, avvisi): i nomi dei file
# relativi alla cartella di uscita, una voce di riepilogo per ordinazione e gli
# avvisi. Un'ordinazione che non si puo' fare (distanze tutte nulle, un solo
# asse positivo) si salta con un avviso; stress alto, stress quasi nullo e
# mancata convergenza della NMDS sono avvisi e il risultato si scrive.
scrivi_ordinazioni <- function(matrici, disegno, parametri, seme, cartella) {
  dir.create(cartella, showWarnings = FALSE)
  metodi <- as.character(unlist(parametri$ord$methods))
  prodotti <- character()
  riepilogo <- list()
  avvisi <- list()
  for (distanza in as.character(unlist(parametri$ord$distances))) {
    matrice <- matrici[[distanza]]
    if (max(matrice) <= 0) {
      # Tutti i campioni coincidono per questa distanza: non c'e' nulla da ordinare.
      avvisi[[length(avvisi) + 1L]] <- .avviso("E-ECO-25", sprintf(
        "Distanza %s: tutte le distanze fra i campioni sono nulle.", distanza))
      next
    }
    for (metodo in metodi) {
      base <- paste(metodo, distanza, sep = "_")
      sottotitolo <- NULL
      if (metodo == "pcoa") {
        o <- pcoa_di(matrice, parametri$ord$pcoa_correction)
        if (ncol(o$coordinate) < 2L) {
          avvisi[[length(avvisi) + 1L]] <- .avviso("E-ECO-25", sprintf(
            "PCoA sulla distanza %s: un solo asse con autovalore positivo.", distanza))
          next
        }
        coordinate <- o$coordinate
        colnames(coordinate) <- paste0("asse_", seq_len(ncol(coordinate)))
        correzione <- if (o$applicata) {
          sprintf("correzione %s applicata agli autovalori negativi (%d)", o$correzione, o$negativi)
        } else if (o$negativi > 0L) {
          sprintf("nessuna correzione (ord.pcoa_correction = none); autovalori negativi: %d", o$negativi)
        } else {
          sprintf("nessun autovalore negativo: la correzione dichiarata (%s) non serve e non e' applicata", o$correzione)
        }
        frazioni <- if (o$applicata) {
          "frazione_di_varianza: autovalore corretto / somma degli autovalori corretti"
        } else {
          "frazione_di_varianza: autovalore / somma di tutti gli autovalori, negativi compresi"
        }
        autovalori <- o$autovalori
        for (colonna in c("autovalore", "autovalore_corretto", "frazione_di_varianza")) {
          autovalori[[colonna]] <- reale(autovalori[[colonna]])
        }
        scrivi_tsv(autovalori, file.path(cartella, paste0(base, "_autovalori.tsv")),
                   intestazione = c(
                     sprintf("PCoA (ape::pcoa) sulla distanza %s.", distanza),
                     correzione, frazioni))
        prodotti <- c(prodotti, paste0("ordinazione/", base, "_autovalori.tsv"))
        intestazione <- c(sprintf("PCoA (ape::pcoa) sulla distanza %s: coordinate dei campioni.", distanza),
                          correzione,
                          "il segno di ogni asse e' fissato: l'elemento di modulo massimo e' positivo")
        etichette <- sprintf("asse %d (%.1f%% della varianza)", 1:2, 100 * o$frazione[1:2])
        riepilogo[[length(riepilogo) + 1L]] <- list(
          metodo = "pcoa", distanza = distanza, correzione = o$correzione,
          correzione_applicata = o$applicata, autovalori_negativi = o$negativi,
          frazione_asse_1 = o$frazione[1L], frazione_asse_2 = o$frazione[2L])
      } else if (metodo == "tsne") {
        perplessita <- parametri$ord$tsne_perplexity
        theta <- parametri$ord$tsne_theta
        o <- tsne_di(matrice, perplessita, parametri$ord$tsne_max_iter, theta, seme)
        coordinate <- o$coordinate
        intestazione <- c(
          sprintf("t-SNE (Rtsne::Rtsne) sulla matrice di distanza %s: is_distance = TRUE, pca = FALSE, dims = 2, num_threads = 1.", distanza),
          sprintf("perplessita' %s; iterazioni %d; theta %s%s; seme %d (disposizione iniziale normale con deviazione 1e-4 estratta dal seme; esagerazione iniziale per 250 iterazioni); costo finale (divergenza di Kullback-Leibler): %s",
                  format(perplessita), o$iterazioni, format(theta),
                  if (theta == 0) " (t-SNE esatto)" else " (Barnes-Hut)",
                  as.integer(seme), reale(o$costo)),
          "il t-SNE conserva i vicini di ogni campione, non le distanze: le distanze fra gruppi lontani e le dimensioni dei gruppi non hanno significato (van der Maaten e Hinton 2008); e' una figura descrittiva, i test restano la PERMANOVA")
        etichette <- c("t-SNE 1", "t-SNE 2")
        sottotitolo <- sottotitolo_tsne(perplessita, theta, o$iterazioni, seme)
        riepilogo[[length(riepilogo) + 1L]] <- list(
          metodo = "tsne", distanza = distanza, perplessita = perplessita,
          iterazioni = o$iterazioni, theta = theta, costo_finale = o$costo)
      } else {
        o <- nmds_di(matrice, as.integer(parametri$ord$nmds_trymax), seme)
        coordinate <- o$coordinate
        intestazione <- c(
          sprintf("NMDS (vegan::metaMDS, monoMDS) sulla matrice di distanza %s, k = 2, regressione monotona globale, stress di tipo 1.", distanza),
          sprintf("avvii casuali: %d; seme %d; stress: %s; soluzione convergente: %s",
                  o$avvii, as.integer(seme), reale(o$stress),
                  if (o$convergente) "si" else "no"),
          "nessuna trasformazione dei dati e nessuna riscalatura in unita' di semicambiamento; assi ruotati sulle componenti principali, segno fissato")
        etichette <- c("NMDS 1", "NMDS 2")
        if (o$stress > SOGLIA_STRESS) {
          avvisi[[length(avvisi) + 1L]] <- .avviso("E-ECO-17", sprintf(
            "NMDS sulla distanza %s: stress %s, sopra %s.", distanza, reale(o$stress),
            format(SOGLIA_STRESS)))
        }
        if (o$quasi_nullo) {
          avvisi[[length(avvisi) + 1L]] <- .avviso("E-ECO-28", sprintf(
            "NMDS sulla distanza %s: stress %s con %d campioni.", distanza,
            reale(o$stress), nrow(matrice)))
        }
        if (!o$convergente) {
          avvisi[[length(avvisi) + 1L]] <- .avviso("E-ECO-24", sprintf(
            "NMDS sulla distanza %s: nessuna coppia di soluzioni simili in %d avvii; si riporta la soluzione di stress minimo.",
            distanza, o$avvii))
        }
        riepilogo[[length(riepilogo) + 1L]] <- list(
          metodo = "nmds", distanza = distanza, stress = o$stress,
          stress_quasi_nullo = o$quasi_nullo,
          convergente = o$convergente, avvii = o$avvii)
      }
      scrivi_tsv(affianca(list(campione = rownames(coordinate)), coordinate, reale), file.path(cartella, paste0(base, "_coordinate.tsv")),
                 intestazione = intestazione)
      prodotti <- c(prodotti, paste0("ordinazione/", base, "_coordinate.tsv"))
      indice <- 0L
      for (v in disegno$variabili) {
        indice <- indice + 1L
        if (!length(v$dimensioni)) next
        nome <- paste0(base, "_", nome_file(indice, v$nome), ".png")
        .grafico_ordinazione(coordinate, v,
                             sprintf("%s, distanza %s", NOMI_ORDINAZIONE[[metodo]],
                                     nome_distanza(distanza)),
                             etichette, file.path(cartella, nome), sottotitolo)
        prodotti <- c(prodotti, paste0("ordinazione/", nome))
      }
      # I grafici dei soli colori (ord.color_by): nessun test li accompagna, e
      # il titolo lo dice.
      indice <- 0L
      for (v in disegno$colore) {
        indice <- indice + 1L
        if (!length(v$dimensioni)) next
        nome <- paste0(base, "_colore_", nome_file(indice, v$nome), ".png")
        .grafico_ordinazione(coordinate, v,
                             sprintf("%s, distanza %s, descrittivo, nessun test",
                                     NOMI_ORDINAZIONE[[metodo]], nome_distanza(distanza)),
                             etichette, file.path(cartella, nome), sottotitolo)
        prodotti <- c(prodotti, paste0("ordinazione/", nome))
      }
    }
  }
  list(file = prodotti, riepilogo = riepilogo, avvisi = avvisi)
}
