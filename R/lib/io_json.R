# Lettura e scrittura dei file JSON del contratto con l'orchestratore.
#
# Interfaccia:
#   leggi_json(percorso)            -> lista
#   scrivi_json(oggetto, percorso)  scrittura atomica, restituisce il percorso
#   scrivi_atomico(righe, percorso) scrittura atomica di testo
#
# La scrittura e' atomica: il contenuto va in un file temporaneo nella stessa
# cartella, che poi viene rinominato. Un processo ucciso a meta' lascia il
# temporaneo, mai un file troncato con il nome definitivo: per la
# dichiarazione d'esito e' cio' che distingue un processo che ha dichiarato da
# uno che non ha potuto farlo.

leggi_json <- function(percorso) {
  # I vettori JSON diventano vettori R, gli oggetti liste con nome; le tabelle
  # restano liste, perche' una conversione implicita in data.frame
  # cambierebbe la forma del dato a seconda del contenuto.
  jsonlite::read_json(
    percorso,
    simplifyVector = TRUE,
    simplifyDataFrame = FALSE,
    simplifyMatrix = FALSE
  )
}

scrivi_atomico <- function(righe, percorso) {
  temporaneo <- tempfile(pattern = ".scrittura-", tmpdir = dirname(percorso))
  on.exit(unlink(temporaneo), add = TRUE)
  con <- file(temporaneo, open = "w", encoding = "UTF-8")
  writeLines(righe, con)
  close(con)
  if (!file.rename(temporaneo, percorso)) {
    stop("impossibile scrivere ", percorso)
  }
  invisible(percorso)
}

scrivi_json <- function(oggetto, percorso) {
  # auto_unbox: un valore singolo diventa uno scalare JSON. Un elenco che deve
  # restare un vettore JSON anche con un solo elemento va passato come lista.
  testo <- jsonlite::toJSON(
    oggetto,
    auto_unbox = TRUE,
    null = "null",
    na = "null",
    digits = NA,
    pretty = TRUE
  )
  scrivi_atomico(testo, percorso)
}
