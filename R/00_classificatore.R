# Verifica del classificatore addestrato di IdTaxa (tax.idtaxa_trained).
#
# Non e' una fase. Lo lancia la validazione degli ingressi, con il
# classificatore IdTaxa, dopo aver verificato il checksum del file: lo carica e
# verifica che contenga un oggetto di DECIPHER::LearnTaxa con la tabella dei
# ranghi (R/lib/idtaxa.R), altrimenti si ferma con E-S0-28. Scrive nella
# cartella della validazione:
#
#   classificatore_addestrato.json   i ranghi che il classificatore nomina, e
#                                    quelli che non hanno una colonna nella
#                                    tassonomia; il numero dei taxa
#
# Parametri: addestrato (il percorso del file .RData); facoltativo
# elenco_dei_taxa, il percorso di un file temporaneo in cui scrivere i nomi dei
# taxa del classificatore: la validazione vi cerca i taxa che la configurazione
# nomina (E-S0-32). Non e' un artefatto: lo legge e lo toglie chi lo chiede.
#
# Scoprire un file sbagliato qui costa qualche secondo; scoprirlo
# all'assegnazione tassonomica costerebbe tutta l'inferenza delle varianti.

for (f in c("io_json.R", "errors.R", "idtaxa.R")) {
  source(file.path(Sys.getenv("AMPLICON16S_R_LIB"), f))
}

esegui_fase(function(parametri, cartella) {
  addestrato <- tryCatch(
    carica_addestrato(parametri$addestrato),
    error = function(e) errore_catalogo("E-S0-28", conditionMessage(e))
  )
  ranghi <- ranghi_dell_addestrato(addestrato)
  elenco <- facoltativo(parametri, "elenco_dei_taxa")
  if (!is.null(elenco)) {
    scrivi_json(
      list(taxa = as.list(sort(unique(as.character(addestrato$taxa)), method = "radix"))),
      elenco
    )
  }
  noti <- c(RANGO_RADICE, unlist(COLONNE_DEI_RANGHI))
  scrivi_json(
    list(ranghi = as.list(ranghi),
         ranghi_scartati = as.list(setdiff(ranghi, noti)),
         taxa = length(addestrato$taxonomy)),
    file.path(cartella, "classificatore_addestrato.json")
  )
  "classificatore_addestrato.json"
})
