# S3 - modello d'errore di sequenziamento, uno per corsa (dada2::learnErrors).
#
# Scrive in 04_error_models/, per ciascun modello:
#
#   modello_<nome>.rds   il modello stimato, letto da S4
#   modello_<nome>.png   il grafico diagnostico: errori osservati e stimati
#                        in funzione della qualita' (dada2::plotErrors)
#
# e convergenza.json, con l'esito dell'auto-consistenza di ogni modello.
#
# Parametri: modelli (nome -> file filtrati, gia' nell'ordine in cui usarli),
# nbases, max_consist, seme, processi. L'ordine e la scelta dei campioni li
# decide la fase Python: learnErrors riceve l'elenco gia' mescolato con il seme
# di run.seed e randomize = FALSE, e usa i file nell'ordine dato finche' le
# basi non superano nbases. Cosi' si sa esattamente su quali campioni e' stata
# fatta la stima.
#
# La convergenza si legge dall'oggetto, non dal testo di un avviso: il ciclo di
# auto-consistenza di dada2 si ferma quando la matrice prodotta e' identica a
# una di quelle usate in ingresso (err_in), oppure dopo MAX_CONSIST giri. Il
# modello e' convergente se err_out e' identica a una delle err_in.
#
# Il grafico e' un PNG del dispositivo cairo, che non registra data ne' ora:
# a parita' di ingresso e' identico byte per byte, come il modello. Un PDF non
# lo sarebbe, perche' ne registra la data di creazione.

for (f in c("io_json.R", "errors.R", "letture.R")) {
  source(file.path(Sys.getenv("AMPLICON16S_R_LIB"), f))
}

# Razionale sistemistico: DADA2 interrompe il ciclo di auto-consistenza quando
# la matrice dei tassi di transizione stimata (err_out) coincide esattamente
# con una delle matrici di ingresso del ciclo (err_in), oppure al raggiungimento
# di MAX_CONSIST. Verificare l'identita' matriciale con identical() evita di
# dipendere dal parsing dei messaggi testuali di console o dalla lingua locale.
convergente <- function(modello) {
  any(vapply(modello$err_in, identical, logical(1), modello$err_out))
}

esegui_fase(function(parametri, cartella) {
  richiedi_pacchetti(c("dada2", "ggplot2"))

  esiti <- list()
  artefatti <- character()
  for (nome in sort(names(parametri$modelli))) {
    file <- unlist(parametri$modelli[[nome]])
    set.seed(as.integer(parametri$seme))
    # Razionale biologico e sistemistico: l'ordine dei campioni e' gia' stato
    # permutato in modo deterministico da Python usando run.seed; randomize = FALSE
    # impone a learnErrors di leggere i file nella sequenza ricevuta fino al
    # superamento di nbases, rendendo tracciabile al singolo campione il pool
    # di basi usato per discriminare errori fisici di lettura da varianti ASV reali.
    modello <- dada2::learnErrors(
      file,
      nbases = as.numeric(parametri$nbases),
      randomize = FALSE,
      MAX_CONSIST = as.integer(parametri$max_consist),
      multithread = as.integer(parametri$processi),
      verbose = 0
    )

    rds <- sprintf("modello_%s.rds", nome)
    saveRDS(modello, file.path(cartella, rds))

    # Razionale sistemistico: l'uso esplicito del backend grafico "cairo" per
    # l'esportazione PNG omette metadati variabili (come i timestamp di creazione
    # presenti nei PDF), garantendo che due esecuzioni sugli stessi dati producano
    # file diagnostici identici byte per byte (stesso checksum SHA-256 e MD5).
    png <- sprintf("modello_%s.png", nome)
    grafico <- dada2::plotErrors(modello, nominalQ = TRUE) +
      ggplot2::ggtitle(sprintf("Modello d'errore: %s", nome))
    grDevices::png(file.path(cartella, png), width = 1200, height = 1200,
                   res = 120, type = "cairo")
    print(grafico)
    invisible(grDevices::dev.off())

    esiti[[nome]] <- list(
      convergenza = convergente(modello),
      iterazioni = length(modello$err_in)
    )
    artefatti <- c(artefatti, rds, png)
  }

  scrivi_json(esiti, file.path(cartella, "convergenza.json"))
  c(artefatti, "convergenza.json")
})
