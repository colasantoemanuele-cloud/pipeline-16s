# Dichiarazione dell'esito verso l'orchestratore, con i codici del catalogo.
#
# Richiede io_json.R. Interfaccia:
#   esegui_fase(principale)          punto d'ingresso di ogni script di fase
#   errore_catalogo(codice, messaggio)  interrompe dichiarando un codice
#   facoltativo(parametri, nome, predefinito)  legge un parametro facoltativo;
#                                    ogni altro accesso a un parametro non
#                                    ricevuto e' un errore
#   richiedi_pacchetti(pacchetti)    carica i pacchetti di calcolo, o si ferma
#                                    con il motivo per cui non si caricano
#
# Uno script di fase ha questa forma:
#
#   for (f in c("io_json.R", "errors.R", "letture.R")) {
#     source(file.path(Sys.getenv("AMPLICON16S_R_LIB"), f))
#   }
#   esegui_fase(function(parametri, cartella) {
#     ...
#     if (condizione) errore_catalogo("E-S2-01", "campione X azzerato")
#     ...
#     c("artefatto1.rds", "artefatto2.tsv")  # nomi relativi a `cartella`
#   })
#
# Il codice di uscita di un processo non puo' portare un codice come E-S2-03,
# quindi l'esito viaggia in un file JSON, scritto per ultimo e in modo
# atomico: la sua assenza e' cio' che l'orchestratore legge come processo
# morto senza aver dichiarato nulla. Il catalogo dei codici sta in Python;
# qui se ne verifica solo la forma, e l'orchestratore rifiuta un codice che
# non vi compare.

PROTOCOLLO_PONTE <- 1L

# I parametri ricevuti dalla fase. La fase Python passa soltanto parametri che
# ha dichiarato, e li legge attraverso una vista che rifiuta gli altri. Qui la
# stessa regola vale per lo script: in R leggere un nome assente da una lista
# restituisce NULL, e per molte funzioni NULL significa "usa il valore
# predefinito della libreria", cioe' un parametro dimenticato che funziona.
# Leggerlo solleva invece un errore; un parametro davvero facoltativo si legge
# con facoltativo(), che lo dice apertamente.
# Razionale sistemistico: il metodo S3 '$.parametri_dichiarati' e
# '[[.parametri_dichiarati' intercetta ogni accesso per nome prima che R possa
# restituire NULL, garantendo la simmetria col contratto VistaConfig di Python.
`$.parametri_dichiarati` <- function(x, name) {
  valori <- unclass(x)
  if (!name %in% names(valori)) {
    stop("parametro non ricevuto: ", name, ". La fase passa allo script solo ",
         "i parametri che dichiara")
  }
  valori[[name]]
}

`[[.parametri_dichiarati` <- function(x, i, ...) {
  if (is.character(i)) {
    return(`$.parametri_dichiarati`(x, i))
  }
  unclass(x)[[i]]
}

# Consente di interrogare esplicitamente un parametro opzionale restituendo il
# valore predefinito indicato senza sollevare l'errore di parametro non ricevuto.
facoltativo <- function(parametri, nome, predefinito = NULL) {
  valori <- unclass(parametri)
  if (nome %in% names(valori)) valori[[nome]] else predefinito
}

# Carica i pacchetti di cui la fase ha bisogno. requireNamespace(quietly =
# TRUE) restituirebbe solo FALSE, e un caricamento fallito per memoria
# esaurita ("cannot allocate vector", "failed to map segment from shared
# object") sembrerebbe un pacchetto mancante: qui il messaggio originale resta
# nell'errore, e il ponte lo riconosce come memoria esaurita.
richiedi_pacchetti <- function(pacchetti) {
  for (pacchetto in pacchetti) {
    tryCatch(
      loadNamespace(pacchetto),
      error = function(e) {
        stop("pacchetto ", pacchetto, " non caricabile: ", conditionMessage(e),
             ". La fase va eseguita nell'ambiente della pipeline, il container",
             call. = FALSE)
      }
    )
  }
  invisible(TRUE)
}

# Codice di uscita di un errore dichiarato. E' solo una conferma: il codice
# del catalogo e' nella dichiarazione.
USCITA_ERRORE_DICHIARATO <- 3L

errore_catalogo <- function(codice, messaggio = "") {
  if (!is.character(codice) || length(codice) != 1L ||
      !grepl("^E-[A-Z0-9]+-[0-9]{2}$", codice)) {
    stop("codice di errore malformato: ", paste(codice, collapse = ", "))
  }
  condizione <- structure(
    class = c("errore_catalogo", "error", "condition"),
    list(message = messaggio, call = sys.call(-1L), codice = codice)
  )
  stop(condizione)
}

esegui_fase <- function(principale) {
  # Gli avvisi vanno sull'uscita di errore appena si presentano, cosi' il log
  # strutturato li raccoglie nell'ordine in cui sono avvenuti.
  options(warn = 1L)

  argomenti <- commandArgs(trailingOnly = TRUE)
  if (length(argomenti) != 1L) {
    stop("atteso un solo argomento, il file della richiesta")
  }
  richiesta <- leggi_json(argomenti[[1L]])
  if (!identical(as.integer(richiesta$protocollo), PROTOCOLLO_PONTE)) {
    stop("protocollo ", richiesta$protocollo, ", atteso ", PROTOCOLLO_PONTE)
  }

  parametri <- richiesta$parametri
  if (is.null(parametri)) parametri <- list()
  class(parametri) <- "parametri_dichiarati"

  esito <- tryCatch(
    {
      artefatti <- principale(parametri, richiesta$cartella_fase)
      if (is.null(artefatti)) artefatti <- character()
      if (!is.character(artefatti)) {
        stop("la funzione principale deve restituire i nomi degli artefatti")
      }
      list(stato = "riuscito", codice = NULL, messaggio = "",
           artefatti = as.list(artefatti))
    },
    errore_catalogo = function(e) {
      message("[", e$codice, "] ", conditionMessage(e))
      list(stato = "errore_catalogo", codice = e$codice,
           messaggio = conditionMessage(e), artefatti = list())
    },
    error = function(e) {
      message("errore non catalogato: ", conditionMessage(e))
      list(stato = "errore_non_catalogato", codice = NULL,
           messaggio = conditionMessage(e), artefatti = list())
    }
  )

  scrivi_json(
    c(list(protocollo = PROTOCOLLO_PONTE, invocazione = richiesta$invocazione),
      esito),
    richiesta$esito
  )
  quit(
    save = "no",
    status = if (identical(esito$stato, "riuscito")) 0L else USCITA_ERRORE_DICHIARATO
  )
}
