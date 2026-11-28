# Composizione tassonomica: abbondanze relative agglomerate a un rango.
#
# La usa R/analisi.R (comando eco run). Scrive, nella cartella composizione/
# delle uscite, abbondanze_relative.tsv e, per ogni variabile della
# composizione (comp.group_by; senza, le variabili del disegno), una tabella e
# un grafico gruppi_NN_<variabile>, con NN la posizione nell'elenco; con
# comp.facet_wrap e comp.facet_grid i grafici a pannelli (ggplot2)
# pannelli_<colonna>.png e griglia_<righe>_<colonne>.png, ciascuno con la
# tabella dei valori disegnati accanto. Il principio dei grafici a pannelli e'
# quello di plot_bars del pacchetto AMBER di Aldo Dale, senza dplyr.
#
# Richiede uscite.R e campioni.R. Interfaccia:
#   agglomera(ps, campioni, rango)   conteggi per taxon del rango, con il
#                                    lignaggio e la categoria dei non assegnati
#   abbondanze_relative(conteggi)    proporzioni per campione (colonne a somma 1)
#   scrivi_composizione(...)         tabella completa, tabelle e grafici per
#                                    variabile, grafici a pannelli (ggplot2)
#
# Razionale biologico: con letture corte molte varianti non hanno
# un'assegnazione ai ranghi piu' fini. Con il
# predefinito di phyloseq::tax_glom (NArm = TRUE) sparirebbero in silenzio e le
# abbondanze relative degli altri taxa risulterebbero gonfiate. Qui restano, in
# una categoria dichiarata, e le proporzioni di ogni campione sommano a 1.

NON_ASSEGNATO <- "non assegnato"
ALTRI <- "altri"

# La sovrapposizione fra barre vicine dei grafici per campione, in unita'
# dell'asse x (una barra e' larga 1). Il dispositivo disegna un rettangolo dal
# bordo sinistro e dalla larghezza: con barre larghe esattamente 1 il bordo
# destro di una barra e il sinistro della seguente possono differire
# nell'ultima cifra e cadere ai due lati del centro di un pixel, che nessuna
# delle due copre: resta una riga bianca verticale, che sembra un campione
# senza barra. Un millesimo di sovrapposizione (meno di un decimo di pixel) la
# chiude e non si vede. Per la stessa ragione l'asse x non ha margini: con il
# margine predefinito fra la prima e l'ultima barra e il bordo del pannello
# resterebbe una striscia bianca. Le barre sono rettangoli impilati qui
# (geom_rect), non geom_col: position_stack avvisa, a ragione, se gli
# intervalli si sovrappongono.
SOVRAPPOSIZIONE_BARRE <- 0.001

# I conteggi dei campioni indicati sommati per taxon del rango. Restituisce
# list(conteggi, lignaggio, taxon): la matrice taxa x campioni, il lignaggio di
# ogni riga fino al rango e il nome al rango. Passi: fusione delle varianti con
# lo stesso lignaggio; riconoscimento delle varianti senza assegnazione al
# rango; somma per lignaggio, con una sola riga per tutti i non assegnati.
agglomera <- function(ps, campioni, rango) {
  sotto <- phyloseq::prune_samples(campioni, ps)
  # tax_glom fonde le varianti con lo stesso lignaggio fino al rango; con
  # NArm = FALSE tiene anche quelle senza assegnazione al rango. Su una
  # tassonomia a un solo rango tax_glom fallisce: la somma per nome che segue
  # fa da sola la stessa cosa.
  fuso <- if (length(phyloseq::rank_names(sotto)) > 1L) {
    phyloseq::tax_glom(sotto, taxrank = rango, NArm = FALSE)
  } else {
    sotto
  }
  conteggi <- conteggi_di(fuso)[, campioni, drop = FALSE]
  tassonomia <- methods::as(phyloseq::tax_table(fuso), "matrix")
  fino_a <- match(rango, colnames(tassonomia))
  tassonomia <- tassonomia[rownames(conteggi), seq_len(fino_a), drop = FALSE]
  nome <- tassonomia[, fino_a]
  # Senza nome, o col solo prefisso di rango ("g__"): un rango non assegnato.
  senza <- is.na(nome) | !nzchar(trimws(nome)) | grepl("^[A-Za-z]__$", nome)
  lignaggio <- apply(tassonomia, 1L, function(riga) paste(riga, collapse = ";"))
  lignaggio[senza] <- NON_ASSEGNATO
  nome[senza] <- NON_ASSEGNATO
  # Una sola categoria per i non assegnati, qualunque sia il loro lignaggio ai
  # ranghi superiori; i taxa assegnati restano distinti per lignaggio.
  somma <- rowsum(conteggi, group = lignaggio, reorder = FALSE)
  primo <- match(rownames(somma), lignaggio)
  list(conteggi = somma, lignaggio = rownames(somma), taxon = unname(nome[primo]))
}

# Le proporzioni per campione: ogni colonna divisa per la sua somma.
abbondanze_relative <- function(conteggi) {
  sweep(conteggi, 2L, colSums(conteggi), "/")
}

# Scrive le tabelle e i grafici della composizione. agglomerato e' l'esito di
# agglomera, disegno quello di esamina, top_n il numero di taxa mostrati per
# variabile. Passi: tabella completa delle abbondanze relative, con i taxa in
# ordine decrescente di abbondanza media; poi, per ogni variabile della
# composizione (comp.group_by, o le variabili del disegno) con almeno un
# gruppo, la media per gruppo dei primi taxa, con i rimanenti riuniti in una
# sola voce, in una tabella e in un grafico a barre. Un taxon il cui nome si
# ripete in lignaggi diversi e' etichettato con il lignaggio. Restituisce i
# list(file, pannelli): i nomi dei file, relativi alla cartella di uscita
# (relativa e' la cartella della composizione: composizione/, o
# composizione/<rango>/ con piu' ranghi), e, con comp.facet_wrap, i numeri del
# grafico per campione per il riepilogo (NULL senza).
scrivi_composizione <- function(agglomerato, disegno, rango, top_n, cartella,
                                relativa = "composizione") {
  dir.create(cartella, showWarnings = FALSE, recursive = TRUE)
  relative <- abbondanze_relative(agglomerato$conteggi)
  media <- rowMeans(relative)
  ordine <- order(-media, agglomerato$lignaggio, method = "radix")
  relative <- relative[ordine, , drop = FALSE]
  media <- media[ordine]
  taxon <- agglomerato$taxon[ordine]
  lignaggio <- agglomerato$lignaggio[ordine]

  tabella <- affianca(list(taxon = taxon, lignaggio = lignaggio,
                           abbondanza_media = reale(media)), relative, reale)
  prodotti <- file.path(relativa, "abbondanze_relative.tsv")
  scrivi_tsv(
    tabella, file.path(cartella, "abbondanze_relative.tsv"),
    intestazione = c(
      sprintf("Abbondanze relative per campione, agglomerate al rango %s.", rango),
      "conteggi non rarefatti; ogni colonna di campione somma a 1",
      sprintf("'%s': tutte le varianti senza assegnazione al rango %s", NON_ASSEGNATO, rango),
      "lignaggio: i ranghi fino a quello dichiarato; taxa con lo stesso nome e lignaggio diverso restano distinti",
      "abbondanza_media: media sui campioni analizzati; le righe sono in ordine decrescente di questa media"
    )
  )

  # Etichette dei primi taxa: il nome al rango, o il lignaggio se il nome si ripete.
  etichette <- taxon
  ripetuti <- etichette %in% etichette[duplicated(etichette)]
  etichette[ripetuti] <- lignaggio[ripetuti]
  indice <- 0L
  for (v in disegno$composizione) {
    indice <- indice + 1L
    gruppi <- names(v$dimensioni)
    if (!length(gruppi)) next
    per_gruppo <- vapply(gruppi, function(g) {
      membri <- names(v$valori)[!is.na(v$valori) & v$valori == g]
      rowMeans(relative[, membri, drop = FALSE])
    }, numeric(nrow(relative)))
    per_gruppo <- matrix(per_gruppo, nrow = nrow(relative),
                         dimnames = list(NULL, gruppi))
    primi <- seq_len(min(top_n, nrow(per_gruppo)))
    mostrati <- rbind(per_gruppo[primi, , drop = FALSE], 1 - colSums(per_gruppo[primi, , drop = FALSE]))
    nomi <- c(etichette[primi], ALTRI)
    base <- nome_file(indice, v$nome)
    righe <- affianca(list(taxon = nomi), mostrati, reale)
    scrivi_tsv(
      righe, file.path(cartella, paste0("gruppi_", base, ".tsv")),
      intestazione = c(
        sprintf("Abbondanza relativa media per gruppo della variabile '%s', rango %s.",
                v$nome, rango),
        sprintf("primi %d taxa per abbondanza relativa media sui campioni analizzati; '%s': tutti i rimanenti",
                length(primi), ALTRI),
        paste0("campioni per gruppo: ",
               paste(sprintf("%s = %d", gruppi, v$dimensioni[gruppi]), collapse = "; ")),
        nota_palette(length(primi))
      )
    )
    grafico_composizione(mostrati, nomi, gruppi, v, rango,
                         file.path(cartella, paste0("gruppi_", base, ".png")))
    prodotti <- c(prodotti, file.path(relativa, paste0("gruppi_", base, c(".tsv", ".png"))))
  }
  # I grafici a pannelli: stessi taxa, stessa voce «altri», stessi colori.
  primi <- seq_len(min(top_n, nrow(relative)))
  mostrate <- rbind(relative[primi, , drop = FALSE],
                    1 - colSums(relative[primi, , drop = FALSE]))
  rownames(mostrate) <- c(etichette[primi], ALTRI)
  wrap <- disegno$pannelli_wrap
  pannelli <- NULL
  # Una colonna senza valori fra i campioni analizzati non ha pannelli.
  if (!is.null(wrap) && length(wrap$dimensioni)) {
    pannelli <- .pannelli_per_campione(mostrate, wrap, rango, cartella, relativa)
    prodotti <- c(prodotti, pannelli$file)
    pannelli$file <- NULL
  }
  griglia <- disegno$pannelli_griglia
  if (!is.null(griglia) && length(Reduce(intersect, lapply(griglia, function(v) {
    names(v$valori)[!is.na(v$valori)]
  })))) {
    prodotti <- c(prodotti, .griglia_per_gruppo(mostrate, griglia, rango, cartella, relativa))
  }
  list(file = prodotti, pannelli = pannelli)
}

# I colori dei taxa nei grafici a pannelli, per nome: quelli delle barre per
# gruppo nello stesso ordine (il taxon piu' abbondante per primo), il grigio
# per «altri», l'ultimo.
.colori_dei_taxa <- function(taxa) {
  stats::setNames(c(colori_gruppi(length(taxa) - 1L), "grey70"), taxa)
}

# Le dimensioni di un grafico a pannelli: ogni pannello almeno 500 pixel di
# larghezza (o 8 per barra, se le barre sono di piu') e 900 di altezza, e
# 700 pixel per la legenda: l'immagine cresce con i pannelli.
.dimensioni_pannelli <- function(colonne, righe, barre) {
  c(larghezza = as.integer(700L + colonne * max(500L, 8L * barre)),
    altezza = as.integer(350L + righe * 900L))
}

# Il tema dei grafici a pannelli: sfondo bianco, testo leggibile a 200 punti
# per pollice, legenda a destra.
.tema_pannelli <- function() {
  ggplot2::theme_bw(base_size = 11) + ggplot2::theme(
    panel.grid.major.x = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(),
    strip.text = ggplot2::element_text(size = 10), legend.position = "right",
    plot.title = ggplot2::element_text(face = "bold"))
}

# Barre per campione, un pannello per gruppo della colonna (facet_wrap, scale
# libere sull'asse x): v e' la variabile descrittiva della colonna, mostrate
# la matrice taxa x campioni dei primi taxa e di «altri», nell'ordine
# dell'abbondanza media. Dentro ogni pannello i campioni sono in ordine
# decrescente di abbondanza del primo taxon (il piu' abbondante nel
# complesso), a parita' per nome: chi guarda vede subito come varia il taxon
# dominante. La legenda e' nell'ordine della pila, come nelle barre per
# gruppo: in alto «altri», che sta in cima, in basso il primo taxon. I
# campioni senza valore non compaiono. Un campione senza abbondanze (nessuna lettura al
# rango, quindi proporzioni non calcolabili) non e' disegnato: al suo posto non
# resta uno spazio vuoto, e il grafico, la tabella e il riepilogo lo dicono.
# Scrive il grafico e la tabella dei valori disegnati; restituisce
# list(file, colonna, campioni_disegnati, campioni_senza_abbondanze).
.pannelli_per_campione <- function(mostrate, v, rango, cartella, relativa) {
  campioni <- names(v$valori)[!is.na(v$valori)]
  senza_abbondanze <- campioni[colSums(!is.finite(mostrate[, campioni, drop = FALSE])) > 0L]
  campioni <- setdiff(campioni, senza_abbondanze)
  primo_taxon <- rownames(mostrate)[1L]
  campioni <- campioni[order(v$valori[campioni], -mostrate[1L, campioni], campioni,
                             method = "radix")]
  dimensioni <- table(factor(unname(v$valori[campioni]), levels = names(v$dimensioni)))
  dimensioni <- stats::setNames(as.integer(dimensioni), names(dimensioni))
  gruppi <- names(dimensioni)[dimensioni > 0L]
  etichetta <- stats::setNames(sprintf("%s (n = %d)", gruppi, dimensioni[gruppi]), gruppi)
  nota_senza <- if (length(senza_abbondanze)) sprintf(
    "campioni senza abbondanze al rango %s (nessuna lettura), non disegnati: %d (%s)",
    rango, length(senza_abbondanze), paste(sort(senza_abbondanze), collapse = ", "))
  taxa <- rownames(mostrate)
  lunga <- data.frame(
    campione = rep(campioni, each = length(taxa)),
    gruppo = rep(unname(v$valori[campioni]), each = length(taxa)),
    taxon = rep(taxa, times = length(campioni)),
    abbondanza = as.vector(mostrate[, campioni, drop = FALSE]),
    stringsAsFactors = FALSE)
  base <- paste0("pannelli_", gsub("[^A-Za-z0-9_]", "_", v$nome))
  scrivi_tsv(
    data.frame(campione = lunga$campione, gruppo = lunga$gruppo, taxon = lunga$taxon,
               abbondanza_relativa = reale(lunga$abbondanza), stringsAsFactors = FALSE),
    file.path(cartella, paste0(base, ".tsv")), intestazione = c(
      sprintf("Composizione per campione al rango %s, un pannello per gruppo di '%s' (comp.facet_wrap): i valori del grafico %s.png.",
              rango, v$nome, base),
      sprintf("primi %d taxa per abbondanza relativa media sui campioni analizzati, come nelle barre per gruppo; '%s': tutti i rimanenti",
              length(taxa) - 1L, ALTRI),
      sprintf("righe nell'ordine del grafico: per gruppo, e nel gruppo per abbondanza decrescente di '%s', il taxon piu' abbondante nel complesso",
              primo_taxon),
      sprintf("campioni senza valore in '%s', non disegnati: %d", v$nome, length(v$mancanti)),
      sprintf("campioni senza abbondanze al rango %s (nessuna lettura), non disegnati: %d%s",
              rango, length(senza_abbondanze),
              if (length(senza_abbondanze)) paste0(" (", paste(sort(senza_abbondanze), collapse = ", "), ")") else ""),
      paste0("campioni per gruppo: ", paste(sprintf("%s = %d", gruppi, dimensioni[gruppi]),
                                             collapse = "; ")),
      nota_palette(length(taxa) - 1L)))
  # Le barre come rettangoli: la posizione del campione nel suo pannello
  # (1, 2, ...), la pila dal primo taxon in basso a «altri» in alto.
  posizione <- stats::ave(seq_along(campioni), unname(v$valori[campioni]), FUN = seq_along)
  lunga$x <- rep(posizione, each = length(taxa))
  lunga$alto <- as.vector(apply(mostrate[, campioni, drop = FALSE], 2L, cumsum))
  lunga$basso <- lunga$alto - lunga$abbondanza
  lunga$gruppo <- factor(etichetta[lunga$gruppo], levels = unname(etichetta))
  lunga$taxon <- factor(lunga$taxon, levels = rev(taxa))
  colonne <- ceiling(sqrt(length(gruppi)))
  righe <- ceiling(length(gruppi) / colonne)
  misure <- .dimensioni_pannelli(colonne, righe, max(dimensioni))
  mezza <- 0.5 + SOVRAPPOSIZIONE_BARRE / 2
  grafico <- ggplot2::ggplot(lunga, ggplot2::aes(xmin = .data$x - mezza, xmax = .data$x + mezza,
                                                  ymin = .data$basso, ymax = .data$alto,
                                                  fill = .data$taxon)) +
    ggplot2::geom_rect() +
    ggplot2::facet_wrap(ggplot2::vars(.data$gruppo), scales = "free_x", ncol = colonne) +
    ggplot2::scale_fill_manual(values = .colori_dei_taxa(taxa),
                               breaks = rev(taxa), name = sprintf("taxon (%s)", rango)) +
    ggplot2::scale_x_continuous(expand = c(0, 0)) +
    ggplot2::scale_y_continuous(expand = c(0, 0), limits = c(0, 1.0001)) +
    ggplot2::labs(title = sprintf("Composizione per campione al rango %s, per %s (%d campioni)",
                                  rango, v$nome, length(campioni)),
                  x = sprintf("campioni (una barra per campione), in ordine decrescente di %s, il taxon più abbondante",
                              primo_taxon),
                  y = "abbondanza relativa",
                  caption = nota_senza) +
    .tema_pannelli() +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                   axis.ticks.x = ggplot2::element_blank())
  apri_png(file.path(cartella, paste0(base, ".png")), larghezza = misure[["larghezza"]],
           altezza = misure[["altezza"]])
  print(grafico)
  grDevices::dev.off()
  list(file = file.path(relativa, paste0(base, c(".tsv", ".png"))), colonna = v$nome,
       campioni_disegnati = length(campioni),
       campioni_senza_abbondanze = length(senza_abbondanze))
}

# Barre medie per gruppo della variabile del disegno in una griglia di
# pannelli (facet_grid): righe e colonne sono le due colonne di
# comp.facet_grid. Ogni pannello porta il numero dei suoi campioni; i campioni
# senza valore in una delle tre colonne non entrano. La legenda e'
# nell'ordine della pila, come nelle barre per gruppo. Scrive il grafico e la
# tabella dei valori disegnati; restituisce i nomi dei file.
.griglia_per_gruppo <- function(mostrate, griglia, rango, cartella, relativa) {
  x <- griglia$x
  r <- griglia$righe
  k <- griglia$colonne
  campioni <- Reduce(intersect, lapply(list(x, r, k), function(v) names(v$valori)[!is.na(v$valori)]))
  chiave <- paste(r$valori[campioni], k$valori[campioni], x$valori[campioni], sep = "
")
  celle <- unique(chiave)
  celle <- celle[order(r$valori[campioni][match(celle, chiave)],
                       k$valori[campioni][match(celle, chiave)],
                       x$valori[campioni][match(celle, chiave)], method = "radix")]
  taxa <- rownames(mostrate)
  righe <- lapply(celle, function(cella) {
    membri <- campioni[chiave == cella]
    primo <- membri[1L]
    data.frame(riga = r$valori[[primo]], colonna = k$valori[[primo]], gruppo = x$valori[[primo]],
               campioni = length(membri), taxon = taxa,
               abbondanza = rowMeans(mostrate[, membri, drop = FALSE]), stringsAsFactors = FALSE)
  })
  lunga <- do.call(rbind, righe)
  rownames(lunga) <- NULL
  base <- paste0("griglia_", gsub("[^A-Za-z0-9_]", "_", r$nome), "_",
                 gsub("[^A-Za-z0-9_]", "_", k$nome))
  scrivi_tsv(
    data.frame(riga = lunga$riga, colonna = lunga$colonna, gruppo = lunga$gruppo,
               campioni = as.integer(lunga$campioni), taxon = lunga$taxon,
               abbondanza_relativa_media = reale(lunga$abbondanza), stringsAsFactors = FALSE),
    file.path(cartella, paste0(base, ".tsv")), intestazione = c(
      sprintf("Composizione media al rango %s per gruppo di '%s', in una griglia: righe '%s', colonne '%s' (comp.facet_grid); i valori del grafico %s.png.",
              rango, x$nome, r$nome, k$nome, base),
      sprintf("primi %d taxa per abbondanza relativa media sui campioni analizzati, come nelle barre per gruppo; '%s': tutti i rimanenti",
              length(taxa) - 1L, ALTRI),
      "campioni: i campioni della barra; abbondanza_relativa_media: media delle abbondanze relative dei campioni",
      sprintf("campioni senza valore in una delle tre colonne, non disegnati: %d",
              length(unique(c(x$mancanti, r$mancanti, k$mancanti)))),
      nota_palette(length(taxa) - 1L)))
  # Il numero dei campioni di ogni pannello, scritto in alto a sinistra.
  pannello <- paste(lunga$riga, lunga$colonna, sep = "
")
  per_pannello <- tapply(campioni, paste(r$valori[campioni], k$valori[campioni], sep = "
"),
                         length)
  testi <- unique(lunga[, c("riga", "colonna")])
  testi$etichetta <- sprintf("n = %d", as.integer(per_pannello[paste(testi$riga, testi$colonna,
                                                                       sep = "
")]))
  lunga$taxon <- factor(lunga$taxon, levels = rev(taxa))
  livelli_x <- names(x$dimensioni)
  lunga$gruppo <- factor(lunga$gruppo, levels = livelli_x)
  n_righe <- length(unique(lunga$riga))
  n_colonne <- length(unique(lunga$colonna))
  dimensioni <- .dimensioni_pannelli(n_colonne, n_righe, length(livelli_x) * 6L)
  grafico <- ggplot2::ggplot(lunga, ggplot2::aes(x = .data$gruppo, y = .data$abbondanza,
                                                  fill = .data$taxon)) +
    ggplot2::geom_col(width = 0.85) +
    ggplot2::geom_text(data = testi, ggplot2::aes(x = -Inf, y = Inf, label = .data$etichetta),
                       inherit.aes = FALSE, hjust = -0.1, vjust = 1.4, size = 3.2) +
    ggplot2::facet_grid(rows = ggplot2::vars(.data$riga), cols = ggplot2::vars(.data$colonna),
                        drop = TRUE, scales = "free_x", space = "free_x") +
    ggplot2::scale_fill_manual(values = .colori_dei_taxa(taxa),
                               breaks = rev(taxa), name = sprintf("taxon (%s)", rango)) +
    ggplot2::scale_x_discrete(drop = TRUE) +
    ggplot2::scale_y_continuous(expand = c(0, 0), limits = c(0, 1.12),
                                breaks = seq(0, 1, 0.25)) +
    ggplot2::labs(title = sprintf("Composizione media per %s al rango %s: righe %s, colonne %s (%d campioni)",
                                  x$nome, rango, r$nome, k$nome, length(campioni)),
                  x = x$nome, y = "abbondanza relativa media") +
    .tema_pannelli() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  apri_png(file.path(cartella, paste0(base, ".png")), larghezza = dimensioni[["larghezza"]],
           altezza = dimensioni[["altezza"]])
  print(grafico)
  grDevices::dev.off()
  file.path(relativa, paste0(base, c(".tsv", ".png")))
}

# Il grafico a barre impilate della composizione media per gruppo di una
# variabile: una barra per gruppo, un colore per taxon, grigio per la voce dei
# rimanenti; sotto ogni barra il gruppo con il numero dei suoi campioni.
grafico_composizione <- function(mostrati, nomi, gruppi, v, rango, percorso) {
  colori <- c(colori_gruppi(length(nomi) - 1L), "grey70")
  # Il margine destro ospita la legenda: con nomi di taxa lunghi cresce, e
  # l'immagine con lui, perche' nessun nome resti tagliato; fino a circa sedici
  # caratteri il margine resta al minimo, 16.
  margine <- max(16, ceiling(0.8 * max(nchar(nomi))) + 4)
  # Con molti gruppi l'immagine si allarga: le barre non si stringono.
  apri_png(percorso, larghezza = larghezza_per_gruppi(length(gruppi)) +
             as.integer(round((margine - 16) * 28)))
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::par(mar = c(9, 4.5, 3, margine), xpd = NA, cex = 0.7)
  posizioni <- graphics::barplot(
    mostrati, col = colori, border = NA, names.arg = rep("", length(gruppi)),
    ylim = c(0, 1), ylab = "abbondanza relativa media",
    main = sprintf("%s: composizione al rango %s", v$nome, rango)
  )
  graphics::text(posizioni, -0.03, labels = sprintf("%s (n = %d)", gruppi, v$dimensioni[gruppi]),
                 srt = 45, adj = 1)
  graphics::legend(
    x = max(posizioni) + (posizioni[1L] * 1.4), y = 1, legend = rev(nomi),
    fill = rev(colori), border = NA, bty = "n"
  )
}
