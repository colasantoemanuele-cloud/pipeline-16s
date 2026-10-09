# Validazione dell'analisi sull'oggetto, prima di ogni calcolo.
#
# Parametri: oggetto (percorso della copia), la configurazione risolta
# (design, run, alpha, comp, beta, stat) e valori_mancanti. Scrive
# validazione.json con gli errori, gli avvisi e il riepilogo del disegno: non
# dichiara un errore del catalogo, perche' i rifiuti possono essere piu' d'uno
# e l'orchestratore li riporta tutti.

for (f in c("io_json.R", "errors.R")) source(file.path(Sys.getenv("AMPLICON16S_R_LIB"), f))
for (f in c("uscite.R", "campioni.R")) source(file.path(Sys.getenv("AMPLICON16S_ECO_R_LIB"), f))

esegui_fase(function(parametri, cartella) {
  richiedi_pacchetti(c("phyloseq", "ape"))
  letto <- leggi_oggetto(parametri$oggetto)
  if (is.null(letto$oggetto)) {
    esito <- list(errori = list(list(codice = "E-ECO-04", dettaglio = paste0(
      toupper(substr(letto$motivo, 1L, 1L)), substring(letto$motivo, 2L), "."))),
      avvisi = list(), disegno = NULL)
  } else {
    esame <- esamina(letto$oggetto, parametri)
    esito <- list(
      errori = esame$errori, avvisi = esame$avvisi,
      disegno = if (!is.null(esame$disegno)) riepilogo_del_disegno(esame$disegno)
    )
  }
  scrivi_json(esito, file.path(cartella, "validazione.json"))
  "validazione.json"
})
