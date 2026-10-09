# Analisi ecologiche di base sull'oggetto finale della pipeline.
#
# Parametri: oggetto (percorso della copia), uscita (cartella in cui scrivere),
# la configurazione risolta e valori_mancanti. L'oggetto si legge e non si
# scrive mai. Ogni file prodotto e' elencato in riepilogo.json insieme al
# disegno, alle profondita' e agli avvisi.

for (f in c("io_json.R", "errors.R")) source(file.path(Sys.getenv("AMPLICON16S_R_LIB"), f))
for (f in c("uscite.R", "campioni.R", "alfa.R", "composizione.R", "beta.R",
            "ordinazione.R", "test.R")) {
  source(file.path(Sys.getenv("AMPLICON16S_ECO_R_LIB"), f))
}

esegui_fase(function(parametri, cartella) {
  richiedi_pacchetti(c("phyloseq", "vegan", "ape", "permute"))
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
  piano <- disegno_dei_test(d, parametri)
  avvisi <- c(esame$avvisi, piano$avvisi)
  seme <- parametri$run$seed
  # Radice dell'albero per UniFrac: riguarda le sole distanze.
  radice <- parametri$beta$unifrac_root
  radicato <- NA
  ps_distanze <- ps
  if (!is.null(radice)) {
    if (radice == "midpoint") richiedi_pacchetti("phangorn")
    radicato <- ape::is.rooted(phyloseq::phy_tree(ps))
    ps_distanze <- radica(ps, radice)
  }
  ordinazioni <- list()
  esiti_dei_test <- list()

  # Un avviso di R che nessun passo ha previsto non si perde sull'uscita di
  # errore: finisce nel riepilogo con il suo codice.
  imprevisti <- character()
  prodotti <- withCallingHandlers({
    file <- character()

    scrivi_tsv(tabella_dei_campioni(d, piano$piani), file.path(uscita, "campioni.tsv"), intestazione = c(
      "Campioni dell'oggetto di partenza e loro stato in ogni analisi.",
      "stato: stato generale; stato_alfa: nell'alfa diversita'; stato:<variabile>: nelle analisi di quella variabile; permanova:<variabile>: nella PERMANOVA e nelle dispersioni",
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
      distanza(metodo, ps_distanze, conteggi, parametri$beta$clr_pseudocount)
    }), metodi)
    file <- c(file, scrivi_distanze(matrici, parametri$beta$clr_pseudocount,
                                    sum(varianti_presenti(conteggi)),
                                    file.path(uscita, "beta"), radice, radicato))

    # Ordinazione: con meno di quattro campioni due dimensioni non hanno senso.
    if (length(unlist(parametri$ord$methods))) {
      if (length(d$campioni) < 4L) {
        avvisi[[length(avvisi) + 1L]] <- list(codice = "E-ECO-25", dettaglio = sprintf(
          "Campioni analizzati: %d; ne servono almeno 4.", length(d$campioni)))
      } else {
        ordinate <- scrivi_ordinazioni(matrici, d, parametri, seme,
                                       file.path(uscita, "ordinazione"))
        # Senza alcuna ordinazione producibile la cartella non resta vuota.
        if (!length(ordinate$file)) unlink(file.path(uscita, "ordinazione"), recursive = TRUE)
        file <- c(file, ordinate$file)
        ordinazioni <- ordinate$riepilogo
        avvisi <- c(avvisi, ordinate$avvisi)
      }
    }

    # Test: alfa diversita' fra gruppi, PERMANOVA e dispersioni.
    tabelle_alfa <- list()
    tabelle_multivariate <- list()
    indice_variabile <- 0L
    for (p in piano$piani) {
      indice_variabile <- indice_variabile + 1L
      esito_alfa <- test_alfa(indici, p)
      tabelle_alfa[[length(tabelle_alfa) + 1L]] <- esito_alfa
      if (!is.null(esito_alfa)) {
        nome <- paste0("alfa/gruppi_", nome_file(indice_variabile, p$nome), ".png")
        grafico_alfa(indici, p, d$profondita, as.integer(parametri$stat$min_group_size),
                     file.path(uscita, nome))
        file <- c(file, nome)
      }
      eseguite <- character()
      for (metodo in metodi) {
        esito <- permanova(matrici[[metodo]], metodo, p, d$tecnici, piano$valori_strati,
                           as.integer(parametri$stat$permanova_permutations),
                           parametri$stat$significance_level, seme)
        if (is.null(esito)) next
        eseguite <- c(eseguite, metodo)
        tabelle_multivariate[[length(tabelle_multivariate) + 1L]] <- esito$tabella
        if (!esito$dispersioni_valutabili) {
          avvisi[[length(avvisi) + 1L]] <- list(codice = "E-ECO-27", dettaglio = sprintf(
            "Variabile '%s', distanza %s: somma dei quadrati residua nulla con %d campioni in %d gruppi.",
            p$nome, metodo, length(p$permanova$campioni), length(p$permanova$gruppi)))
        }
        if (esito$dispersioni_diverse) {
          avvisi[[length(avvisi) + 1L]] <- list(codice = "E-ECO-19", dettaglio = sprintf(
            "Variabile '%s', distanza %s: p delle dispersioni %s, sotto %s.", p$nome,
            metodo, reale(esito$p_dispersioni), format(parametri$stat$significance_level)))
        }
      }
      esiti_dei_test[[length(esiti_dei_test) + 1L]] <- list(
        variabile = p$nome,
        alfa = list(eseguito = !is.null(esito_alfa), campioni = length(p$alfa$campioni),
                    gruppi = I(p$alfa$gruppi)),
        permanova = list(eseguita = length(eseguite) > 0L,
                         campioni = length(p$permanova$campioni),
                         gruppi = I(p$permanova$gruppi),
                         termini_tecnici = I(p$permanova$termini),
                         termine_della_variabile_stimabile = p$permanova$stimabile,
                         modello_saturo = p$permanova$saturo,
                         distanze = I(eseguite))
      )
    }
    file <- c(file, scrivi_test(
      if (length(tabelle_alfa)) do.call(rbind, tabelle_alfa),
      if (length(tabelle_multivariate)) do.call(rbind, tabelle_multivariate),
      parametri, d$profondita, file.path(uscita, "test")))
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
  if (!is.null(radice)) {
    riepilogo$beta$albero <- list(radice = radice, radicato_nell_oggetto = radicato)
  }
  riepilogo$ordinazione <- ordinazioni
  riepilogo$test <- esiti_dei_test
  riepilogo$avvisi <- avvisi
  riepilogo$ambiente <- list(
    R = as.character(getRversion()),
    pacchetti = lapply(
      stats::setNames(nm = c("ape", "permute", if (identical(radice, "midpoint")) "phangorn",
                             "phyloseq", "vegan")),
      function(p) as.character(utils::packageVersion(p))
    )
  )
  riepilogo$file <- I(prodotti)
  scrivi_json(riepilogo, file.path(uscita, "riepilogo.json"))
  character()
})
