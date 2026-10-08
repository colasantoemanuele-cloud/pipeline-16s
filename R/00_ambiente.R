# L'ambiente R in cui la pipeline sta per calcolare.
#
# Non e' una fase. Lo lancia l'esecutore, prima di qualunque fase, quando la
# regola rigorosa sulla provenienza e' attiva (run.strict_provenance), su tutti
# i pacchetti del file di blocco; e lo lancia S8 a ogni esecuzione, sul solo
# dada2, per sapere quale versione R carichera' e se porta la correzione dei
# pareggi di assignTaxonomy. Scrive nella cartella indicata:
#
#   ambiente_r.json   la versione di R e, per ogni pacchetto richiesto, la
#                     versione installata e la correzione dichiarata nel suo
#                     DESCRIPTION, se c'e'; null per un pacchetto assente
#
# Parametri: pacchetti (i nomi dei pacchetti da leggere), nome (facoltativo:
# il nome del file da scrivere, se diverso da ambiente_r.json).
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
  voci <- lapply(richiesti, function(nome) {
    # Il pacchetto che R caricherebbe: il primo con quel nome nelle librerie.
    if (length(find.package(nome, quiet = TRUE)) == 0L) return(NULL)
    d <- utils::packageDescription(nome)
    correzione <- if (is.null(d[["Amplicon16sPatch"]]) || is.na(d[["Amplicon16sPatch"]])) NULL else list(
      Base = d[["Amplicon16sBase"]], File = d[["Amplicon16sPatch"]],
      SHA256 = d[["Amplicon16sPatchSHA256"]]
    )
    list(versione = as.character(utils::packageVersion(nome)), correzione = correzione)
  })
  names(voci) <- richiesti
  nome <- facoltativo(parametri, "nome", "ambiente_r.json")
  scrivi_json(
    list(r = paste(R.version$major, R.version$minor, sep = "."), pacchetti = voci),
    file.path(cartella, nome)
  )
  nome
})
