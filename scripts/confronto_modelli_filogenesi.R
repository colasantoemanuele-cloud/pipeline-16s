#!/usr/bin/env Rscript
# Confronto fra i modelli evolutivi GTR+G+I e GTR+G sulle varianti finali di
# un'esecuzione, con i criteri d'informazione (docs/decision_log.md, 1.5).
#
# Uso, nell'ambiente della pipeline:
#   Rscript scripts/confronto_modelli_filogenesi.R <ps_filtrato.rds> <uscita.json> [processi]
#
# I due modelli si adattano con lo stesso allineamento, lo stesso albero di
# partenza e la stessa ricerca della fase S9 (R/09_phylogeny.R): allineamento
# con DECIPHER::AlignSeqs, albero neighbor-joining su distanze JC69, massima
# verosimiglianza con phangorn, gamma a quattro categorie, scambi fra rami
# vicini (NNI). Cambia solo la quota di siti invarianti: stimata in GTR+G+I,
# assente in GTR+G. Il numero di osservazioni dei criteri e' il numero di
# colonne dell'allineamento.

suppressPackageStartupMessages({
  library(phyloseq); library(DECIPHER); library(phangorn)
})
argomenti <- commandArgs(trailingOnly = TRUE)
if (length(argomenti) < 2) stop("uso: confronto_modelli_filogenesi.R <ps_filtrato.rds> <uscita.json> [processi]")
processi <- if (length(argomenti) >= 3) as.integer(argomenti[3]) else 1L

ps <- readRDS(argomenti[1])
sequenze <- refseq(ps)
set.seed(100)
allineamento <- AlignSeqs(sequenze, processors = processi, verbose = FALSE)
dati <- phyDat(as(allineamento, "matrix"), type = "DNA")
partenza <- NJ(dist.ml(dati))

adatta <- function(invarianti) {
  inizio <- proc.time()[["elapsed"]]
  modello <- pml(partenza, data = dati, model = "GTR", k = 4L, inv = if (invarianti) 0.2 else 0)
  modello <- optim.pml(modello, model = "GTR", optGamma = TRUE, optInv = invarianti, optNni = TRUE,
                       rearrangement = "NNI", control = pml.control(trace = 0L))
  list(
    log_verosimiglianza = signif(modello$logLik, 10),
    parametri = modello$df,
    AIC = signif(AIC(modello), 10),
    BIC = signif(BIC(modello), 10),
    forma_gamma = signif(modello$shape, 8),
    siti_invarianti = signif(modello$inv, 8),
    lunghezza_totale = signif(sum(modello$tree$edge.length), 8),
    secondi = round(proc.time()[["elapsed"]] - inizio)
  )
}

esito <- list(
  varianti = length(sequenze),
  colonne_allineamento = unique(Biostrings::width(allineamento)),
  colonne_costanti = sum(apply(as(allineamento, "matrix"), 2, function(x) length(unique(x[x != "-"])) <= 1)),
  modelli = list("GTR+G+I" = adatta(TRUE), "GTR+G" = adatta(FALSE))
)
esito$differenza_BIC <- signif(esito$modelli[["GTR+G+I"]]$BIC - esito$modelli[["GTR+G"]]$BIC, 6)
esito$differenza_AIC <- signif(esito$modelli[["GTR+G+I"]]$AIC - esito$modelli[["GTR+G"]]$AIC, 6)
writeLines(jsonlite::toJSON(esito, auto_unbox = TRUE, pretty = TRUE, digits = NA), argomenti[2])
cat(readLines(argomenti[2]), sep = "\n")
