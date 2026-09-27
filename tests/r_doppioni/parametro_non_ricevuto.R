# Doppione: legge un parametro che la fase non gli ha passato.
#
# In R leggere un nome assente da una lista darebbe NULL, e lo script
# proseguirebbe con il valore predefinito della libreria. Con i parametri
# dichiarati la lettura fallisce: il test verifica che l'errore arrivi.
for (f in c("io_json.R", "errors.R", "letture.R")) {
  source(file.path(Sys.getenv("AMPLICON16S_R_LIB"), f))
}

esegui_fase(function(parametri, cartella) {
  facoltativo(parametri, "assente")  # dichiaratamente facoltativo: nessun errore
  soglia <- parametri$soglia_mai_passata
  scrivi_json(list(soglia = soglia), file.path(cartella, "non_scritto.json"))
  "non_scritto.json"
})
