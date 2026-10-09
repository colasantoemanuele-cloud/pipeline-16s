# Analisi ecologiche di base sull'oggetto finale della pipeline.
#
# Parametri: oggetto (percorso della copia), uscita (cartella in cui scrivere),
# la configurazione risolta e valori_mancanti. L'oggetto si legge e non si
# scrive mai. Ogni file prodotto e' elencato in riepilogo.json insieme al
# disegno, alle profondita' e agli avvisi.

for (f in c("io_json.R", "errors.R")) source(file.path(Sys.getenv("AMPLICON16S_R_LIB"), f))
for (f in c("uscite.R", "campioni.R", "alfa.R", "composizione.R", "beta.R")) {
  source(file.path(Sys.getenv("AMPLICON16S_ECO_R_LIB"), f))
}

esegui_fase(function(parametri, cartella) {
  richiedi_pacchetti(c("phyloseq", "vegan", "ape"))
  uscita <- parametri$uscita
  letto <- leggi_oggetto(parametri$oggetto)
  if (is.null(letto$oggetto)) errore_catalogo("E-ECO-04", letto$motivo)
  ps <- letto$oggetto
  esame <- esamina(ps, parametri)
  if (length(esame$errori)) {
    primo <- esame$errori[[1L]]
    errore_catalogo(primo$codice, primo$dettaglio)
  }
  d <- esame$disegno
  avvisi <- esame$avvisi
  seme <- parametri$run$seed

  # Un avviso di R che nessun passo ha previsto non si perde sull'uscita di
  # errore: finisce nel riepilogo con il suo codice.
  imprevisti <- character()
  prodotti <- withCallingHandlers({
    file <- character()

    scrivi_tsv(tabella_dei_campioni(d), file.path(uscita, "campioni.tsv"), intestazione = c(
      "Campioni dell'oggetto di partenza e loro stato in ogni analisi.",
      "stato: stato generale; stato_alfa: nell'alfa diversita'; stato:<variabile>: nelle analisi di quella variabile",
      "valore:<variabile>: il valore del campione, NA se mancante o se il campione e' escluso"
    ))
    file <- c(file, "campioni.tsv")

    # Alfa diversita'.
    rarefatti_campioni <- setdiff(d$campioni, d$sotto_profondita)
    rarefatti <- rarefai(ps, rarefatti_campioni, d$profondita, seme)
    indici <- indici_alfa(rarefatti)
    file <- c(file, scrivi_alfa(indici, rarefatti, d$letture, d$profondita,
                                d$origine_profondita, seme, file.path(uscita, "alfa")))

    # Composizione.
    agglomerato <- agglomera(ps, d$campioni, parametri$comp$rank)
    file <- c(file, scrivi_composizione(agglomerato, d, parametri$comp$rank,
                                        as.integer(parametri$comp$top_n),
                                        file.path(uscita, "composizione")))

    # Beta diversita'.
    conteggi <- conteggi_di(ps)[, d$campioni, drop = FALSE]
    metodi <- as.character(unlist(parametri$beta$distances))
    matrici <- stats::setNames(lapply(metodi, function(metodo) {
      distanza(metodo, ps, conteggi, parametri$beta$clr_pseudocount)
    }), metodi)
    file <- c(file, scrivi_distanze(matrici, parametri$beta$clr_pseudocount,
                                    file.path(uscita, "beta")))
    file
  }, warning = function(w) {
    imprevisti <<- c(imprevisti, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  for (testo in unique(imprevisti)) {
    avvisi[[length(avvisi) + 1L]] <- list(codice = "E-ECO-16", dettaglio = testo)
  }

  riepilogo <- riepilogo_del_disegno(d)
  riepilogo$composizione <- list(
    rango = parametri$comp$rank,
    taxa = nrow(agglomerato$conteggi),
    categoria_dei_non_assegnati = NON_ASSEGNATO
  )
  riepilogo$beta <- list(distanze = I(metodi))
  riepilogo$avvisi <- avvisi
  riepilogo$ambiente <- list(
    R = as.character(getRversion()),
    pacchetti = lapply(
      stats::setNames(nm = c("ape", "phyloseq", "vegan")),
      function(p) as.character(utils::packageVersion(p))
    )
  )
  riepilogo$file <- I(prodotti)
  scrivi_json(riepilogo, file.path(uscita, "riepilogo.json"))
  character()
})
