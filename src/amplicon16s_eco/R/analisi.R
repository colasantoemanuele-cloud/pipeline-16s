# Analisi ecologiche di base sull'oggetto finale della pipeline.
#
# Lo lancia il comando eco run (esecuzione.py, attraverso ponte.py), dopo la
# validazione di R/valida.R, su una copia dell'oggetto.
#
# L'ORDINE DELLE ANALISI e i file che ciascuna scrive nella cartella di uscita
# (la funzione che scrive e' fra parentesi; NN e' la posizione della variabile
# nella configurazione):
#   1. esame dell'oggetto e piani dei test (lib/campioni.R, lib/test.R,
#      lib/continua.R, lib/ripetute.R): nessun calcolo; un rifiuto ferma qui
#        campioni.tsv
#      con il filtro di abbondanza (filter.*), le varianti che ha tolto
#        filtro/varianti_filtrate.tsv
#   2. alfa diversita': ricchezza, Shannon, Gini-Simpson, Chao1 (lib/alfa.R),
#      sui conteggi rarefatti (alpha.method rarefy, una o piu' rarefazioni)
#        alfa/alfa_diversita.tsv, alfa/conteggi_rarefatti.tsv
#        alfa/alfa_per_rarefazione.tsv, con piu' rarefazioni
#      oppure su quelli non rarefatti (alpha.method none), con la correlazione
#      fra letture e indici
#        alfa/alfa_diversita.tsv, alfa/profondita.tsv
#   3. composizione a ogni rango di comp.rank (lib/composizione.R); con piu'
#      ranghi i file stanno in composizione/<rango>/
#        composizione/abbondanze_relative.tsv
#        composizione/gruppi_NN_<variabile>.tsv e .png, per le variabili
#        di comp.group_by o, senza, per quelle del disegno
#        composizione/pannelli_<colonna>.tsv e .png (comp.facet_wrap)
#        composizione/griglia_<righe>_<colonne>.tsv e .png (comp.facet_grid)
#   4. distanze fra campioni (lib/beta.R)
#        beta/distanza_<metodo>.tsv
#   5. ordinazioni PCoA, NMDS, t-SNE, con almeno quattro campioni
#      (lib/ordinazione.R)
#        ordinazione/<metodo>_<distanza>_coordinate.tsv
#        ordinazione/pcoa_<distanza>_autovalori.tsv
#        ordinazione/<metodo>_<distanza>_NN_<variabile>.png
#        ordinazione/<metodo>_<distanza>_colore_NN_<colonna>.png, per
#        ord.color_by (descrittivi, senza test)
#   6. test fra gruppi: alfa diversita', PERMANOVA e dispersioni (lib/test.R)
#        test/alfa_diversita.tsv, test/permanova.tsv
#        alfa/gruppi_NN_<variabile>.png, solo se il confronto si esegue
#   7. variabili continue, se dichiarate (lib/continua.R)
#        test/alfa_continua.tsv, test/permanova_continua.tsv
#        alfa/continua_NN_<variabile>.png, solo se la correlazione si esegue
#   8. ordinazione vincolata, se constrained.method e' dichiarato
#      (lib/vincolata.R): RDA sui conteggi trasformati o db-RDA per distanza
#        vincolata/varianza.tsv, vincolata/test.tsv
#        vincolata/<modello>_coordinate.tsv, vincolata/<modello>_NN_<variabile>.png
#   9. misure ripetute, se il soggetto e' dichiarato (lib/ripetute.R)
#        test/soggetti.tsv, test/alfa_appaiata.tsv, test/permanova_appaiata.tsv
#  10. riepilogo.json: disegno, profondita', esito di ogni analisi, avvisi,
#      versioni di R e dei pacchetti, elenco dei file prodotti
# L'ordine conta: gli indici del passo 2 servono ai passi 6, 7 e 9, le matrici
# del passo 4 ai passi 5, 6, 7, 8 e 9. Gli altri file della cartella di uscita
# (manifesto, configurazione risolta) li scrive la parte Python.
#
# Parametri: oggetto (percorso della copia), uscita (cartella in cui scrivere),
# la configurazione risolta e valori_mancanti. L'oggetto si legge e non si
# scrive mai. Ogni file prodotto e' elencato in riepilogo.json insieme al
# disegno, alle profondita' e agli avvisi.

for (f in c("io_json.R", "errors.R")) source(file.path(Sys.getenv("AMPLICON16S_R_LIB"), f))
for (f in c("uscite.R", "campioni.R", "alfa.R", "composizione.R", "beta.R",
            "ordinazione.R", "test.R", "continua.R", "ripetute.R", "vincolata.R")) {
  source(file.path(Sys.getenv("AMPLICON16S_ECO_R_LIB"), f))
}

esegui_fase(function(parametri, cartella) {
  richiedi_pacchetti(c("phyloseq", "vegan", "ape", "permute"))
  if ("tsne" %in% unlist(parametri$ord$methods)) richiedi_pacchetti("Rtsne")
  if (!is.null(parametri$comp$facet_wrap) || !is.null(parametri$comp$facet_grid)) {
    richiedi_pacchetti("ggplot2")
  }
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
  # L'oggetto delle analisi: quello letto, oppure, con il filtro di abbondanza,
  # quello filtrato (lib/campioni.R).
  if (!is.null(d$ps)) ps <- d$ps
  # Il livello dei calcoli e il filtro si dichiarano nell'intestazione di ogni
  # tabella.
  assign("NOTE_DEI_DATI", c(
    if (!is.null(d$agglomerato)) sprintf(
      "livello dei calcoli: rango %s (run.analysis_rank): le %d varianti agglomerate in %d taxa con phyloseq::tax_glom(taxrank = \"%s\", NArm = FALSE), le varianti senza nome al rango separate secondo il lignaggio ai ranghi superiori%s",
      d$agglomerato$rango, d$agglomerato$varianti, d$agglomerato$taxa, d$agglomerato$rango,
      if (d$agglomerato$albero) "; UniFrac usa l'albero potato da tax_glom alla variante che rappresenta ogni taxon" else ""),
    if (!is.null(d$filtro)) sprintf(
      "filtro di abbondanza relativa: soglia %s%% (filter.min_relative_abundance), ambito %s (filter.scope); %s con letture: %d prima, %d dopo (filtro/varianti_filtrate.tsv)",
      format(d$filtro$soglia), d$filtro$ambito,
      if (is.null(d$agglomerato)) "varianti" else "taxa", d$filtro$varianti_prima,
      d$filtro$varianti_dopo)
  ), envir = globalenv())
  imposta_palette(parametri$plot$palette)
  assign("RIPETIZIONI_ALFA", d$ripetizioni, envir = globalenv())
  piano <- disegno_dei_test(d, parametri)
  # Variabili continue e misure ripetute: piani vuoti se non dichiarate.
  continue <- piano_continue(d, parametri)
  appaiato <- piano_appaiato(d, parametri)
  # Ordinazione vincolata: piano nullo se non dichiarata.
  vincolato <- piano_vincolato(d, parametri)
  avvisi <- c(esame$avvisi, piano$avvisi, continue$avvisi, appaiato$avvisi, vincolato$avvisi)
  esiti_vincolati <- NULL
  esiti_continue <- list()
  esiti_appaiati <- list()
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

    scrivi_tsv(tabella_dei_campioni(d, piano$piani, vincolato), file.path(uscita, "campioni.tsv"), intestazione = c(
      "Campioni dell'oggetto di partenza e loro stato in ogni analisi.",
      "stato: stato generale; stato_alfa: nell'alfa diversita'; stato:<variabile>: nelle analisi di quella variabile; permanova:<variabile>: nella PERMANOVA e nelle dispersioni",
      "valore:<variabile>: il valore del campione, NA se mancante o se il campione e' escluso",
      if (!is.null(vincolato)) "vincolata: nell'ordinazione vincolata (constrained.*)"
    ))
    file <- c(file, "campioni.tsv")
    if (!is.null(d$filtro)) file <- c(file, scrivi_filtro(d, file.path(uscita, "filtro")))

    # Alfa diversita': sui conteggi rarefatti (una o piu' rarefazioni), o su
    # quelli non rarefatti con la correlazione fra letture e indici.
    if (identical(d$metodo_alfa, "none")) {
      indici <- indici_alfa(conteggi_di(ps)[, d$campioni, drop = FALSE])
      file <- c(file, scrivi_alfa(indici, NULL, d$letture_analisi, NA_integer_,
                                  d$origine_profondita, seme, file.path(uscita, "alfa"),
                                  senza_singleton = d$singleton == 0))
      profondita <- profondita_e_indici(indici, d$letture_analisi,
                                        parametri$stat$significance_level,
                                        file.path(uscita, "alfa"))
      file <- c(file, profondita$file)
      avvisi <- c(avvisi, profondita$avvisi)
    } else {
      rarefatti_campioni <- setdiff(d$campioni, d$sotto_profondita)
      rarefatti <- rarefai(ps, rarefatti_campioni, d$profondita, seme, d$ripetizioni)
      per_rarefazione <- lapply(rarefatti, indici_alfa)
      indici <- media_degli_indici(per_rarefazione)
      file <- c(file, scrivi_alfa(indici, rarefatti, d$letture_analisi, d$profondita,
                                  d$origine_profondita, seme, file.path(uscita, "alfa"),
                                  senza_singleton = d$singleton == 0,
                                  per_rarefazione = per_rarefazione))
    }

    # Composizione: a ogni rango di comp.rank; con piu' ranghi, una
    # sottocartella per rango.
    ranghi_composizione <- as.character(unlist(parametri$comp$rank))
    composizioni <- list()
    for (rango in ranghi_composizione) {
      relativa <- if (length(ranghi_composizione) > 1L) {
        file.path("composizione", gsub("[^A-Za-z0-9_]", "_", rango))
      } else {
        "composizione"
      }
      agglomerato <- agglomera(ps, d$campioni, rango)
      composizioni[[rango]] <- list(rango = rango, taxa = nrow(agglomerato$conteggi),
                                    categoria_dei_non_assegnati = NON_ASSEGNATO)
      scritta <- scrivi_composizione(agglomerato, d, rango,
                                     as.integer(parametri$comp$top_n),
                                     file.path(uscita, relativa), relativa)
      file <- c(file, scritta$file)
      # Solo con comp.facet_wrap: i campioni del grafico per campione.
      if (!is.null(scritta$pannelli)) composizioni[[rango]]$pannelli <- scritta$pannelli
    }

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
    note_alfa <- character()
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
        nota <- nota_palette(length(p$alfa$gruppi))
        if (length(nota)) note_alfa <- c(note_alfa, sprintf("grafico %s: %s", nome, nota))
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
      parametri, d$profondita, file.path(uscita, "test"), note_alfa))

    # Variabili continue: correlazione con l'alfa diversita' e PERMANOVA.
    if (length(continue$piani)) {
      correlazioni <- list()
      con_termine <- list()
      indice_variabile <- 0L
      for (p in continue$piani) {
        indice_variabile <- indice_variabile + 1L
        correlazione <- correlazione_alfa(indici, p)
        correlazioni[[length(correlazioni) + 1L]] <- correlazione
        if (!is.null(correlazione)) {
          nome <- paste0("alfa/continua_", nome_file(indice_variabile, p$nome), ".png")
          grafico_continua(indici, p, d$profondita, file.path(uscita, nome))
          file <- c(file, nome)
        }
        eseguite <- character()
        for (metodo in metodi) {
          esito <- permanova_continua(matrici[[metodo]], metodo, p, d$tecnici,
                                      continue$valori_strati,
                                      as.integer(parametri$stat$permanova_permutations), seme)
          if (is.null(esito)) next
          eseguite <- c(eseguite, metodo)
          con_termine[[length(con_termine) + 1L]] <- esito
        }
        esiti_continue[[length(esiti_continue) + 1L]] <- list(
          variabile = p$nome,
          alfa = list(eseguita = !is.null(correlazione), campioni = length(p$alfa$campioni)),
          permanova = list(eseguita = length(eseguite) > 0L,
                           campioni = length(p$permanova$campioni),
                           termini_tecnici = I(p$permanova$termini), distanze = I(eseguite)))
      }
      file <- c(file, scrivi_continue(
        if (length(correlazioni)) do.call(rbind, correlazioni),
        if (length(con_termine)) do.call(rbind, con_termine),
        parametri, d$profondita, file.path(uscita, "test")))
    }

    # Ordinazione vincolata: RDA o db-RDA, con i test per permutazione.
    if (!is.null(vincolato) && vincolato$eseguibile) {
      esiti <- scrivi_vincolata(vincolato, ps, matrici, d, parametri, seme,
                                file.path(uscita, "vincolata"))
      file <- c(file, esiti$file)
      if (!length(esiti$file)) unlink(file.path(uscita, "vincolata"), recursive = TRUE)
      esiti_vincolati <- esiti$riepilogo
      avvisi <- c(avvisi, esiti$avvisi)
    }

    # Misure ripetute: i test entro soggetto di ogni variabile categorica.
    if (length(appaiato$piani)) {
      entro_alfa <- list()
      entro_multivariati <- list()
      for (p in appaiato$piani) {
        esito_alfa <- test_alfa_appaiato(indici, p)
        entro_alfa[[length(entro_alfa) + 1L]] <- esito_alfa
        eseguite <- character()
        for (metodo in metodi) {
          esito <- permanova_appaiata(matrici[[metodo]], metodo, p,
                                      as.integer(parametri$stat$permanova_permutations),
                                      parametri$stat$significance_level, seme)
          if (is.null(esito)) next
          eseguite <- c(eseguite, metodo)
          entro_multivariati[[length(entro_multivariati) + 1L]] <- esito$tabella
        }
        esiti_appaiati[[length(esiti_appaiati) + 1L]] <- list(
          variabile = p$nome,
          alfa = list(eseguito = !is.null(esito_alfa), livelli = I(p$alfa$livelli),
                      soggetti_completi = length(p$alfa$completi),
                      soggetti_esclusi = length(p$alfa$soggetti) - length(p$alfa$completi)),
          permanova = list(eseguita = length(eseguite) > 0L, livelli = I(p$permanova$livelli),
                           soggetti_completi = length(p$permanova$completi),
                           soggetti_esclusi = length(p$permanova$soggetti) -
                             length(p$permanova$completi),
                           distanze = I(eseguite)))
      }
      file <- c(file, scrivi_appaiati(
        if (length(entro_alfa)) do.call(rbind, entro_alfa),
        if (length(entro_multivariati)) do.call(rbind, entro_multivariati),
        appaiato$piani, parametri, d$profondita, file.path(uscita, "test")))
    }
    file
  }, warning = function(w) {
    imprevisti <<- c(imprevisti, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  for (testo in unique(imprevisti)) {
    avvisi[[length(avvisi) + 1L]] <- list(codice = "E-ECO-16", dettaglio = testo)
  }

  riepilogo <- riepilogo_del_disegno(d)
  # Con un rango solo la voce e' quella di sempre; con piu' ranghi, una per rango.
  riepilogo$composizione <- if (length(composizioni) == 1L) composizioni[[1L]] else
    unname(composizioni)
  riepilogo$beta <- list(distanze = I(metodi))
  if (!is.null(radice)) {
    riepilogo$beta$albero <- list(radice = radice, radicato_nell_oggetto = radicato)
  }
  riepilogo$ordinazione <- ordinazioni
  if (!is.null(vincolato)) {
    riepilogo$vincolata <- list(
      metodo = vincolato$metodo, eseguita = vincolato$eseguibile,
      campioni = length(vincolato$campioni), gradi_di_liberta = vincolato$gl,
      variabili = I(vincolato$variabili), condizioni = I(vincolato$condizioni),
      modelli = if (is.null(esiti_vincolati)) list() else esiti_vincolati)
  }
  riepilogo$test <- esiti_dei_test
  # Presenti solo se le analisi corrispondenti sono dichiarate.
  if (length(continue$piani)) riepilogo$test_continue <- esiti_continue
  if (length(appaiato$piani)) {
    riepilogo$misure_ripetute <- list(soggetto = d$soggetto$nome, test = esiti_appaiati)
  }
  riepilogo$avvisi <- avvisi
  riepilogo$ambiente <- list(
    R = as.character(getRversion()),
    pacchetti = lapply(
      stats::setNames(nm = c("ape", "permute", if (identical(radice, "midpoint")) "phangorn",
                             "phyloseq",
                             switch(origine_palette(parametri$plot$palette),
                                    brewer = "RColorBrewer", viridis = "viridisLite"),
                             if ("tsne" %in% unlist(parametri$ord$methods)) "Rtsne",
                             if (!is.null(parametri$comp$facet_wrap) ||
                                 !is.null(parametri$comp$facet_grid)) "ggplot2",
                             "vegan")),
      function(p) as.character(utils::packageVersion(p))
    )
  )
  riepilogo$file <- I(prodotti)
  scrivi_json(riepilogo, file.path(uscita, "riepilogo.json"))
  character()
})
