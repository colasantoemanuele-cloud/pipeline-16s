# S7 - filtro di lunghezza delle varianti.
#
# Scrive in 07_chimera/, accanto a S6:
#
#   tabella_asv.rds          la tabella senza chimere di S6, con le sole
#                            varianti di lunghezza fra len_min e len_max
#   lunghezze.tsv            per ogni lunghezza presente nella tabella di S6:
#                            varianti, letture, e se e' ammessa
#   letture_lunghezza.tsv    letture per campione dopo il filtro (tracciamento)
#
# Parametri: tabella (il file di S6), senza_letture (accession senza letture
# filtrate, zero nel tracciamento), len_min, len_max.
#
# Con letture troncate a lunghezza fissa e senza fusione di coppie, le
# varianti hanno tutte la lunghezza del troncamento e il filtro non toglie
# nulla. Resta per i dati in cui la lunghezza varia, e per una variante di
# lunghezza inattesa, che qui si vedrebbe in lunghezze.tsv.

for (f in c("io_json.R", "errors.R", "letture.R")) {
  source(file.path(Sys.getenv("AMPLICON16S_R_LIB"), f))
}

esegui_fase(function(parametri, cartella) {
  tabella <- readRDS(parametri$tabella)
  minimo <- as.integer(parametri$len_min)
  massimo <- as.integer(parametri$len_max)

  # Razionale biologico: la regione ipervariabile V4 del gene 16S rRNA
  # (amplificata con la coppia di primer 515F/806R) presenta in natura una
  # lunghezza fisiologica fortemente conservata (circa 250-256 bp sull'amplicone
  # completo, o esattamente la lunghezza di troncamento in letture single-end).
  # Il filtro rigido tra len_min e len_max scarta amplificati non specifici
  # dell'ospite, dimeri residui o inserzioni/delezioni aberranti.
  lunghezza <- nchar(colnames(tabella))
  ammesse <- lunghezza >= minimo & lunghezza <= massimo
  filtrata <- tabella[, ammesse, drop = FALSE]

  per_lunghezza <- sort(unique(lunghezza))
  scrivi_atomico(
    c("lunghezza\tvarianti\tletture\tammessa",
      vapply(per_lunghezza, function(l) {
        colonne <- lunghezza == l
        sprintf("%d\t%d\t%.0f\t%s", l, sum(colonne), sum(tabella[, colonne]),
                if (l >= minimo && l <= massimo) "si" else "no")
      }, character(1))),
    file.path(cartella, "lunghezze.tsv")
  )

  senza_letture <- as.character(unlist(parametri$senza_letture))
  tracciate <- c(rowSums(filtrata),
                 stats::setNames(numeric(length(senza_letture)), senza_letture))
  tracciate <- tracciate[sort(names(tracciate), method = "radix")]

  saveRDS(filtrata, file.path(cartella, "tabella_asv.rds"))
  c("tabella_asv.rds", "lunghezze.tsv",
    traccia_letture(tracciate, "lunghezza", cartella))
})
