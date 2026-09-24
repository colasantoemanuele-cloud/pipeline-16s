# Doppione: termina con successo e scrive un artefatto.
for (f in c("io_json.R", "errors.R", "letture.R")) {
  source(file.path(Sys.getenv("AMPLICON16S_R_LIB"), f))
}

esegui_fase(function(parametri, cartella) {
  cat("uscita standard del doppione\n")
  if (!is.null(parametri$avviso)) message(parametri$avviso)

  # L'artefatto riporta i parametri ricevuti: il test verifica cosi' che il
  # viaggio attraverso il JSON li abbia lasciati intatti.
  scrivi_json(list(ricevuti = parametri), file.path(cartella, "artefatto.json"))

  artefatti <- "artefatto.json"
  if (!is.null(parametri$letture)) {
    conteggi <- unlist(parametri$letture)
    artefatti <- c(artefatti, traccia_letture(conteggi, "doppione", cartella))
  }
  if (isTRUE(parametri$dichiara_inesistente)) {
    artefatti <- c(artefatti, "mai_scritto.rds")
  }
  artefatti
})
