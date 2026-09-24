# Doppione: fallisce dichiarando un codice del catalogo.
#
# Il parametro `modo` produce anche i fallimenti che non sono dichiarazioni:
#   "dichiarato"      (predefinito) errore_catalogo(codice, messaggio)
#   "non_catalogato"  stop() qualunque, senza codice
#   "uscita"          il processo esce con 1 senza dichiarare nulla
#   "segmentazione"   il processo muore per SIGSEGV, attraverso il percorso di
#                     arresto di R per un guasto grave
#   "stallo"          il processo resta fermo senza mai terminare
for (f in c("io_json.R", "errors.R", "letture.R")) {
  source(file.path(Sys.getenv("AMPLICON16S_R_LIB"), f))
}

esegui_fase(function(parametri, cartella) {
  message("il doppione sta per fallire: ", parametri$messaggio)
  modo <- if (is.null(parametri$modo)) "dichiarato" else parametri$modo
  switch(
    modo,
    dichiarato = errore_catalogo(parametri$codice, parametri$messaggio),
    non_catalogato = stop(parametri$messaggio),
    uscita = quit(save = "no", status = 1L),
    # 11 e' SIGSEGV su Linux; tools non ne esporta il nome.
    segmentazione = tools::pskill(Sys.getpid(), 11L),
    stallo = Sys.sleep(3600),
    stop("modo sconosciuto: ", modo)
  )
  "non_raggiunto.json"
})
