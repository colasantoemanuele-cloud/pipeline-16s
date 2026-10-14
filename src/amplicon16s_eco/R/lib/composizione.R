# Composizione tassonomica: abbondanze relative agglomerate a un rango.
#
# La usa R/analisi.R (comando eco run). Scrive, nella cartella composizione/
# delle uscite, abbondanze_relative.tsv e, per ogni variabile, una tabella e
# un grafico gruppi_NN_<variabile>.
#
# Richiede uscite.R e campioni.R. Interfaccia:
#   agglomera(ps, campioni, rango)   conteggi per taxon del rango, con il
#                                    lignaggio e la categoria dei non assegnati
#   abbondanze_relative(conteggi)    proporzioni per campione (colonne a somma 1)
#   scrivi_composizione(...)         tabella completa, tabelle e grafici per
#                                    variabile
#
# Razionale biologico: con letture corte molte varianti non hanno
# un'assegnazione ai ranghi piu' fini. Con il
# predefinito di phyloseq::tax_glom (NArm = TRUE) sparirebbero in silenzio e le
# abbondanze relative degli altri taxa risulterebbero gonfiate. Qui restano, in
# una categoria dichiarata, e le proporzioni di ogni campione sommano a 1.

NON_ASSEGNATO <- "non assegnato"
ALTRI <- "altri"

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
# ordine decrescente di abbondanza media; poi, per ogni variabile con almeno un
# gruppo, la media per gruppo dei primi taxa, con i rimanenti riuniti in una
# sola voce, in una tabella e in un grafico a barre. Un taxon il cui nome si
# ripete in lignaggi diversi e' etichettato con il lignaggio. Restituisce i
# nomi dei file, relativi alla cartella di uscita.
scrivi_composizione <- function(agglomerato, disegno, rango, top_n, cartella) {
  dir.create(cartella, showWarnings = FALSE)
  relative <- abbondanze_relative(agglomerato$conteggi)
  media <- rowMeans(relative)
  ordine <- order(-media, agglomerato$lignaggio, method = "radix")
  relative <- relative[ordine, , drop = FALSE]
  media <- media[ordine]
  taxon <- agglomerato$taxon[ordine]
  lignaggio <- agglomerato$lignaggio[ordine]

  tabella <- affianca(list(taxon = taxon, lignaggio = lignaggio,
                           abbondanza_media = reale(media)), relative, reale)
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
    righe <- affianca(list(taxon = nomi), mostrati, reale)
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

# Il grafico a barre impilate della composizione media per gruppo di una
# variabile: una barra per gruppo, un colore per taxon, grigio per la voce dei
# rimanenti; sotto ogni barra il gruppo con il numero dei suoi campioni.
grafico_composizione <- function(mostrati, nomi, gruppi, v, rango, percorso) {
  colori <- c(colori_gruppi(length(nomi) - 1L), "grey70")
  # Con molti gruppi l'immagine si allarga: le barre non si stringono.
  apri_png(percorso, larghezza = larghezza_per_gruppi(length(gruppi)))
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
