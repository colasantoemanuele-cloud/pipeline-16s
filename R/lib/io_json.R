# Lettura e scrittura dei file JSON del contratto con l'orchestratore.
#
# Interfaccia:
#   leggi_json(percorso)            -> lista
#   scrivi_json(oggetto, percorso)  scrittura atomica, restituisce il percorso
#   scrivi_atomico(righe, percorso) scrittura atomica di testo
#   salva_rds(oggetto, percorso)    scrittura atomica di un .rds, con
#                                   l'intestazione indipendente dalle
#                                   impostazioni locali
#
# La scrittura e' atomica: il contenuto va in un file temporaneo nella stessa
# cartella, che poi viene rinominato. Un processo ucciso a meta' lascia il
# temporaneo, mai un file troncato con il nome definitivo: per la
# dichiarazione d'esito e' cio' che distingue un processo che ha dichiarato da
# uno che non ha potuto farlo.
#
# Ogni file .rds delle fasi si scrive con salva_rds, mai con saveRDS:
# l'intestazione del formato RDS 3 registra la codifica nativa della sessione,
# "UTF-8" con le impostazioni del container (LANG=C.UTF-8) e "ANSI_X3.4-1968"
# con LC_ALL=C, con il resto del file identico. salva_rds serializza sempre con
# una codifica nativa UTF-8, impostando LC_CTYPE per la sola scrittura se la
# sessione non lo e', perche' gli stessi dati diano gli stessi byte su
# qualunque macchina.

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

salva_rds <- function(oggetto, percorso) {
  if (!isTRUE(l10n_info()[["UTF-8"]])) {
    prima <- Sys.getlocale("LC_CTYPE")
    on.exit(Sys.setlocale("LC_CTYPE", prima), add = TRUE)
    impostata <- FALSE
    for (nome in c("C.UTF-8", "C.utf8", "en_US.UTF-8")) {
      if (nzchar(suppressWarnings(Sys.setlocale("LC_CTYPE", nome)))) {
        impostata <- TRUE
        break
      }
    }
    if (!impostata) stop("nessuna impostazione locale UTF-8 disponibile per serializzare")
  }
  temporaneo <- tempfile(pattern = ".scrittura-", tmpdir = dirname(percorso))
  on.exit(unlink(temporaneo), add = TRUE)
  saveRDS(oggetto, temporaneo)
  if (!file.rename(temporaneo, percorso)) {
    stop("impossibile scrivere ", percorso)
  }
  invisible(percorso)
}
