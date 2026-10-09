# Composizione tassonomica: abbondanze relative agglomerate a un rango.
#
# Richiede uscite.R e campioni.R. Interfaccia:
#   agglomera(ps, campioni, rango)   conteggi per taxon del rango, con il
#                                    lignaggio e la categoria dei non assegnati
#   abbondanze_relative(conteggi)    proporzioni per campione (colonne a somma 1)
#   scrivi_composizione(...)         tabella completa, tabelle e grafici per
#                                    variabile
#
# Razionale biologico: a 137 basi molte varianti non hanno un genere. Con il
# predefinito di phyloseq::tax_glom (NArm = TRUE) sparirebbero in silenzio e le
# abbondanze relative degli altri taxa risulterebbero gonfiate. Qui restano, in
# una categoria dichiarata, e le proporzioni di ogni campione sommano a 1.

NON_ASSEGNATO <- "non assegnato"
ALTRI <- "altri"

agglomera <- function(ps, campioni, rango) {
  sotto <- phyloseq::prune_samples(campioni, ps)
  # tax_glom fonde le varianti con lo stesso lignaggio fino al rango; con
  # NArm = FALSE tiene anche quelle senza assegnazione al rango.
  fuso <- phyloseq::tax_glom(sotto, taxrank = rango, NArm = FALSE)
  conteggi <- conteggi_di(fuso)[, campioni, drop = FALSE]
  tassonomia <- methods::as(phyloseq::tax_table(fuso), "matrix")
  fino_a <- match(rango, colnames(tassonomia))
  tassonomia <- tassonomia[rownames(conteggi), seq_len(fino_a), drop = FALSE]
  nome <- tassonomia[, fino_a]
  senza <- is.na(nome) | !nzchar(trimws(nome))
  lignaggio <- apply(tassonomia, 1L, function(riga) paste(riga, collapse = ";"))
  lignaggio[senza] <- NON_ASSEGNATO
  nome[senza] <- NON_ASSEGNATO
  # Una sola categoria per i non assegnati, qualunque sia il loro lignaggio ai
  # ranghi superiori; i taxa assegnati restano distinti per lignaggio.
  somma <- rowsum(conteggi, group = lignaggio, reorder = FALSE)
  primo <- match(rownames(somma), lignaggio)
  list(conteggi = somma, lignaggio = rownames(somma), taxon = unname(nome[primo]))
}

abbondanze_relative <- function(conteggi) {
  sweep(conteggi, 2L, colSums(conteggi), "/")
}

scrivi_composizione <- function(agglomerato, disegno, rango, top_n, cartella) {
  dir.create(cartella, showWarnings = FALSE)
  relative <- abbondanze_relative(agglomerato$conteggi)
  media <- rowMeans(relative)
  ordine <- order(-media, agglomerato$lignaggio, method = "radix")
  relative <- relative[ordine, , drop = FALSE]
  media <- media[ordine]
  taxon <- agglomerato$taxon[ordine]
  lignaggio <- agglomerato$lignaggio[ordine]

  tabella <- data.frame(taxon = taxon, lignaggio = lignaggio,
                        abbondanza_media = reale(media),
                        stringsAsFactors = FALSE, check.names = FALSE)
  for (campione in colnames(relative)) tabella[[campione]] <- reale(relative[, campione])
  prodotti <- "composizione/abbondanze_relative.tsv"
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
  for (v in disegno$variabili) {
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
    righe <- data.frame(taxon = nomi, stringsAsFactors = FALSE, check.names = FALSE)
    for (g in gruppi) righe[[g]] <- reale(mostrati[, g])
    scrivi_tsv(
      righe, file.path(cartella, paste0("gruppi_", base, ".tsv")),
      intestazione = c(
        sprintf("Abbondanza relativa media per gruppo della variabile '%s', rango %s.",
                v$nome, rango),
        sprintf("primi %d taxa per abbondanza relativa media sui campioni analizzati; '%s': tutti i rimanenti",
                length(primi), ALTRI),
        paste0("campioni per gruppo: ",
               paste(sprintf("%s = %d", gruppi, v$dimensioni[gruppi]), collapse = "; "))
      )
    )
    grafico_composizione(mostrati, nomi, gruppi, v, rango,
                         file.path(cartella, paste0("gruppi_", base, ".png")))
    prodotti <- c(prodotti, paste0("composizione/gruppi_", base, c(".tsv", ".png")))
  }
  prodotti
}

grafico_composizione <- function(mostrati, nomi, gruppi, v, rango, percorso) {
  colori <- c(colori_gruppi(length(nomi) - 1L), "grey70")
  apri_png(percorso)
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::par(mar = c(9, 4.5, 3, 16), xpd = NA, cex = 0.7)
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
