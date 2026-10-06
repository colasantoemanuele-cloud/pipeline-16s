# Controllo di avvio - l'ambiente R in cui la pipeline sta per calcolare.
#
# Non e' una fase: lo lancia l'esecutore, prima di qualunque fase, quando la
# regola rigorosa sulla provenienza e' attiva (run.strict_provenance). Scrive
# nella cartella indicata:
#
#   ambiente_r.json   la versione di R e, per ogni pacchetto richiesto, la
#                     versione installata e la correzione dichiarata nel suo
#                     DESCRIPTION, se c'e'; null per un pacchetto assente
#
# Parametri: pacchetti (i nomi dei pacchetti del file di blocco).
#
# Lo script non giudica: riporta cio' che e' installato. Il confronto con il
# file di blocco lo fa la parte Python (runner/provenienza.py), che conosce il
# file. E' la verifica piu' forte che la pipeline possa fare da dentro il
# container, dove il digest dell'immagine in esecuzione non e' conoscibile: le
# versioni lette qui sono quelle delle librerie che i calcoli caricheranno.

for (f in c("io_json.R", "errors.R")) {
  source(file.path(Sys.getenv("AMPLICON16S_R_LIB"), f))
}

esegui_fase(function(parametri, cartella) {
  richiesti <- as.character(unlist(parametri$pacchetti))
  installati <- rownames(utils::installed.packages())
  voci <- lapply(richiesti, function(nome) {
    if (!nome %in% installati) return(NULL)
    d <- utils::packageDescription(nome)
    correzione <- if (is.null(d[["Amplicon16sPatch"]]) || is.na(d[["Amplicon16sPatch"]])) NULL else list(
      Base = d[["Amplicon16sBase"]], File = d[["Amplicon16sPatch"]],
      SHA256 = d[["Amplicon16sPatchSHA256"]]
    )
    list(versione = as.character(utils::packageVersion(nome)), correzione = correzione)
  })
  names(voci) <- richiesti
  scrivi_json(
    list(r = paste(R.version$major, R.version$minor, sep = "."), pacchetti = voci),
    file.path(cartella, "ambiente_r.json")
  )
  "ambiente_r.json"
})
