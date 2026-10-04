# S5 - la tabella delle sequenze (dada2::makeSequenceTable).
#
# Scrive in 06_seqtab/:
#
#   tabella.rds          la tabella campioni x varianti, conteggi interi, con
#                        le varianti in ordine di abbondanza decrescente; e'
#                        la tabella prima della rimozione delle chimere, e
#                        resta qui: S6 la legge e non la ricopia
#   letture_tabella.tsv  letture per campione nella tabella (tracciamento)
#   tabella.json         campioni, varianti, letture, e quante varianti per
#                        ciascuna lunghezza
#
# Parametri: varianti (il file delle varianti per campione di S4),
# senza_letture (accession senza letture filtrate, zero nel tracciamento),
# max_varianti (qc.max_asv_count).
#
# La memoria. makeSequenceTable costruisce una matrice densa di interi,
# campioni x varianti, e la riordina per abbondanza con una copia. Sul dataset
# di riferimento (960 x 13130, 53 MB) il picco sale di 157 MB oltre la lettura
# delle varianti: circa tre volte la matrice, qualunque sia il lotto con cui
# S4 ha elaborato i campioni. Costruirla a lotti non cambierebbe nulla,
# perche' la tabella finale e' comunque una. Per questo il numero di varianti
# distinte si controlla prima di allocarla: oltre max_varianti la fase si
# ferma con E-S5-01, senza aver chiesto la memoria.

for (f in c("io_json.R", "errors.R", "letture.R", "risorse.R")) {
  source(file.path(Sys.getenv("AMPLICON16S_R_LIB"), f))
}

esegui_fase(function(parametri, cartella) {
  richiedi_pacchetti("dada2")

  varianti <- readRDS(parametri$varianti)
  senza_letture <- as.character(unlist(parametri$senza_letture))
  massimo <- as.numeric(parametri$max_varianti)

  distinte <- length(unique(unlist(lapply(varianti, names), use.names = FALSE)))
  if (distinte > massimo) {
    errore_catalogo(
      "E-S5-01",
      sprintf(paste0(
        "%d varianti distinte su %d campioni, oltre le %.0f di qc.max_asv_count: ",
        "la tabella chiederebbe circa %.1f GB"
      ), distinte, length(varianti), massimo,
      3 * 4 * distinte * length(varianti) / 1e9)
    )
  }

  tabella <- dada2::makeSequenceTable(varianti, orderBy = "abundance")
  salva_rds(tabella, file.path(cartella, "tabella.rds"))

  letture <- c(rowSums(tabella),
               stats::setNames(numeric(length(senza_letture)), senza_letture))
  letture <- letture[sort(names(letture), method = "radix")]

  lunghezze <- table(nchar(colnames(tabella)))
  scrivi_json(
    list(
      campioni = nrow(tabella),
      varianti = ncol(tabella),
      letture = sum(letture),
      varianti_per_lunghezza = as.list(stats::setNames(as.integer(lunghezze), names(lunghezze)))
    ),
    file.path(cartella, "tabella.json")
  )
  message(sprintf("tabella %d x %d; memoria di picco: %s",
                  nrow(tabella), ncol(tabella), memoria_di_picco()))

  c("tabella.rds", traccia_letture(letture, "tabella", cartella), "tabella.json")
})
