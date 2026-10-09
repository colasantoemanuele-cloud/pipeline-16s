# Genera un oggetto phyloseq sintetico con una struttura nota.
#
# Uso: Rscript --vanilla genera_scenario.R descrizione.json oggetto.rds
#
# La descrizione ha: seme, taxa (pari), profondita e blocchi. Ogni blocco
# genera n campioni con gli stessi parametri:
#   gruppo, lotto        i valori delle due colonne dei metadati
#   scambia              se vero, la composizione di base ha le coppie di taxa
#                        (1,2), (3,4), ... scambiate: stessa forma, altra
#                        posizione
#   spostamento          somma al logaritmo delle proporzioni +s sulla prima
#                        meta' dei taxa e -s sulla seconda (effetto di lotto)
#   sigma                deviazione standard del rumore sul logaritmo delle
#                        proporzioni: la dispersione fra campioni del blocco
#   presenti             numero di taxa presenti (gli altri a zero); null: tutti
#
# Composizione di base: le proporzioni alternano 4 e 1, cosi' lo scambio delle
# coppie cambia la composizione media lasciando identica, per simmetria, la
# distribuzione delle distanze dentro il gruppo. Le letture di ogni campione
# sono un'estrazione multinomiale alla profondita' indicata.

argomenti <- commandArgs(trailingOnly = TRUE)
d <- jsonlite::read_json(argomenti[[1L]], simplifyVector = FALSE)
RNGkind("Mersenne-Twister", "Inversion", "Rejection")
set.seed(as.integer(d$seme))

taxa <- as.integer(d$taxa)
base <- rep(c(4, 1), taxa / 2L)
scambiata <- rep(c(1, 4), taxa / 2L)
direzione <- rep(c(1, -1), each = taxa / 2L)

colonne <- list()
gruppo <- character()
lotto <- character()
for (b in d$blocchi) {
  for (i in seq_len(as.integer(b$n))) {
    logaritmi <- log(if (isTRUE(b$scambia)) scambiata else base)
    if (!is.null(b$spostamento)) logaritmi <- logaritmi + as.numeric(b$spostamento) * direzione
    logaritmi <- logaritmi + stats::rnorm(taxa, sd = as.numeric(b$sigma))
    p <- exp(logaritmi)
    if (!is.null(b$presenti)) p[seq_len(taxa) > as.integer(b$presenti)] <- 0
    colonne[[length(colonne) + 1L]] <- as.integer(stats::rmultinom(1L, as.integer(d$profondita), p))
    gruppo <- c(gruppo, b$gruppo)
    lotto <- c(lotto, if (is.null(b$lotto)) "unico" else b$lotto)
  }
}
conteggi <- do.call(cbind, colonne)
campioni <- sprintf("c%03d", seq_len(ncol(conteggi)))
varianti <- sprintf("v%03d", seq_len(taxa))
dimnames(conteggi) <- list(varianti, campioni)
tassonomia <- cbind(Regno = "Batteri", Famiglia = sprintf("Fam%03d", seq_len(taxa)),
                    Genere = sprintf("Gen%03d", seq_len(taxa)))
rownames(tassonomia) <- varianti
metadati <- data.frame(gruppo = gruppo, lotto = lotto, row.names = campioni,
                       stringsAsFactors = FALSE)
saveRDS(phyloseq::phyloseq(
  phyloseq::otu_table(conteggi, taxa_are_rows = TRUE),
  phyloseq::tax_table(tassonomia),
  phyloseq::sample_data(metadati)
), argomenti[[2L]])
