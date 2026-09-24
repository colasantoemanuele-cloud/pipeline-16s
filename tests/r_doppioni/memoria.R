# Doppione: esaurisce la memoria, in modo confinato.
#
# Nessuna delle due forme consuma davvero la memoria della macchina:
#   "allocazione"      alloca `elementi` numeri in doppia precisione. Sotto il
#                      limite che il ponte impone al processo figlio
#                      l'allocazione fallisce davvero, dentro R, ed e' R a
#                      intercettarla.
#   "fuori_contratto"  la stessa allocazione prima di esegui_fase: R la
#                      intercetta ma lo script non puo' dichiarare nulla, e
#                      resta solo il messaggio sull'uscita di errore.
#   "terminazione"     il processo riceve SIGKILL, come quando e' il sistema
#                      operativo a ucciderlo per memoria: nessun messaggio,
#                      nessuna dichiarazione.
for (f in c("io_json.R", "errors.R", "letture.R")) {
  source(file.path(Sys.getenv("AMPLICON16S_R_LIB"), f))
}

richiesta <- leggi_json(commandArgs(trailingOnly = TRUE)[[1L]])
if (identical(richiesta$parametri$modo, "fuori_contratto")) {
  x <- numeric(richiesta$parametri$elementi)
}

esegui_fase(function(parametri, cartella) {
  switch(
    parametri$modo,
    allocazione = {
      x <- numeric(parametri$elementi)
      cat("allocati", length(x), "elementi\n")
    },
    terminazione = tools::pskill(Sys.getpid(), tools::SIGKILL),
    stop("modo sconosciuto: ", parametri$modo)
  )
  character()
})
