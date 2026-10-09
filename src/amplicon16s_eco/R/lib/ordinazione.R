# Ordinazione: PCoA e NMDS sulle matrici di distanza.
#
# Richiede uscite.R e alfa.R (fissa_seme). Interfaccia:
#   pcoa_di(matrice, correzione)     coordinate, autovalori, frazioni
#   nmds_di(matrice, trymax, seme)   coordinate, stress, convergenza
#   scrivi_ordinazioni(...)          tabelle e grafici per variabile
#
# Razionale sistemistico: le distanze non euclidee (Bray-Curtis, Jaccard) danno
# alla PCoA autovalori negativi; la correzione e' una scelta dichiarata, e la
# frazione di varianza si riporta sugli autovalori che ne risultano. La NMDS
# riceve la matrice senza attributi: metaMDS non applica trasformazioni ne' la
# riscalatura in unita' di semicambiamento, che scatterebbe da sola per una
# matrice prodotta da vegdist. Il segno di ogni asse e' arbitrario e viene
# fissato, perche' le coordinate non cambino fra macchine.

SOGLIA_STRESS <- 0.2

.fissa_segno <- function(coordinate) {
  for (j in seq_len(ncol(coordinate))) {
    massimo <- which.max(abs(coordinate[, j]))
    if (coordinate[massimo, j] < 0) coordinate[, j] <- -coordinate[, j]
  }
  coordinate
}

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

.grafico_ordinazione <- function(coordinate, v, titolo, etichette, percorso) {
  con_valore <- names(v$valori)[!is.na(v$valori)]
  gruppi <- names(v$dimensioni)
  colori <- stats::setNames(colori_gruppi(length(gruppi)), gruppi)
  apri_png(percorso)
  on.exit(grDevices::dev.off(), add = TRUE)
  colonne <- max(1L, ceiling(length(gruppi) / 30L))
  graphics::par(mar = c(4.5, 4.5, 3, 4 + 11 * colonne), xpd = NA, cex = 0.7)
  graphics::plot(
    coordinate[con_valore, 1L], coordinate[con_valore, 2L],
    col = colori[v$valori[con_valore]], pch = 19, asp = 1,
    xlab = etichette[1L], ylab = etichette[2L],
    main = sprintf("%s (%d campioni)", titolo, length(con_valore))
  )
  limiti <- graphics::par("usr")
  graphics::legend(
    x = limiti[2L] + 0.03 * (limiti[2L] - limiti[1L]), y = limiti[4L],
    legend = sprintf("%s (n = %d)", gruppi, v$dimensioni[gruppi]),
    col = colori, pch = 19, bty = "n", ncol = colonne, title = v$nome
  )
}

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
        .grafico_ordinazione(coordinate, v, sprintf("%s, distanza %s", toupper(metodo), distanza),
                             etichette, file.path(cartella, nome))
        prodotti <- c(prodotti, paste0("ordinazione/", nome))
      }
    }
  }
  list(file = prodotti, riepilogo = riepilogo, avvisi = avvisi)
}
