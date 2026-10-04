# S4 - inferenza delle varianti di sequenza (dada2::dada), campione per
# campione, con il modello d'errore della corsa di ciascuno.
#
# Scrive in 05_asv_inference/:
#
#   varianti_per_campione.rds  le varianti di ogni campione: lista con nome
#                              (accession) di vettori interi con nome
#                              (sequenza -> letture), in ordine di accession;
#                              e' l'ingresso di S5
#   varianti.tsv               per campione: modello, letture, varianti della
#                              prima passata e varianti finali
#   letture_denoised.tsv       letture attribuite a una variante, per campione
#                              (tracciamento, R/lib/letture.R)
#   inferenza.json             il riepilogo: modalita', regola e numero delle
#                              informazioni a priori, conteggi per modello
#
# e, con il pseudo-pooling, anche:
#
#   priori.tsv                 le sequenze usate come informazione a priori,
#                              con prevalenza e abbondanza nella prima passata
#   errori_seconda_passata_<nome>.rds
#                              il modello d'errore della seconda passata, per
#                              ciascun modello di S3 (vedi sotto)
#
# Parametri: campioni (accession -> file filtrato), modelli (accession ->
# file del modello di S3), senza_letture (accession senza letture filtrate),
# pool, omega_a, lotto, processi.
#
# I campioni si elaborano a lotti di `lotto`, e di ogni lotto si conserva solo
# cio' che serve dopo: la memoria dipende dal lotto, non dal numero dei
# campioni. Con pool = TRUE non e' possibile: dada2 riunisce in un'unica
# inferenza tutti i campioni di un modello, e il lotto non ha effetto.
#
# Pseudo-pooling. dada(pool = "pseudo") in una sola chiamata fa due passate
# (sorgente di dada2 1.36.0, funzione dada):
#
#   1. ogni campione e' elaborato da solo, con il modello d'errore fornito;
#   2. le informazioni a priori sono le varianti della prima passata presenti
#      in almeno PSEUDO_PREVALENCE campioni oppure con almeno
#      PSEUDO_ABUNDANCE letture in totale:
#        st <- makeSequenceTable(clustering)
#        colnames(st)[colSums(st > 0) >= PSEUDO_PREVALENCE |
#                     colSums(st) >= PSEUDO_ABUNDANCE]
#      con i valori predefiniti 2 e Inf: conta solo la prevalenza;
#   3. ogni campione e' rielaborato da solo, con le sue sequenze uniche
#      presenti fra le informazioni a priori marcate come tali. Il modello
#      d'errore di questa passata NON e' quello fornito: alla fine di ogni
#      passata dada ricalcola err con errorEstimationFunction (loessErrfun)
#      dalle transizioni di tutti i campioni della chiamata, e la seconda
#      passata usa quello, anche con selfConsist = FALSE.
#
# Qui le due passate sono esplicite, entrambe a lotti. Della prima si
# conservano soltanto, sequenza per sequenza, in quanti campioni compare e
# quante letture ha in totale, e per ciascun modello la somma delle
# transizioni; il resto si scarta alla fine di ogni lotto. Somme di conteggi
# interi sono esatte in qualunque ordine: il risultato non dipende dal lotto.
# Le informazioni a priori si raccolgono su tutti i campioni, di tutte le
# corse; il modello della seconda passata si ricalcola per ciascun modello di
# S3 dalle transizioni dei soli suoi campioni, come farebbe dada in una
# chiamata per corsa. Le soglie si leggono da getDadaOpt(), non si ricopiano.
#
# Un campione solo non attiva il pseudo-pooling (dada: length(derep) <= 1
# implica pool = FALSE), e cosi' qui.

for (f in c("io_json.R", "errors.R", "letture.R", "risorse.R")) {
  source(file.path(Sys.getenv("AMPLICON16S_R_LIB"), f))
}

# Il nome di un file di modello senza prefisso ed estensione: modello_X.rds -> X.
nome_modello <- function(file) {
  sub("^modello_", "", sub("\\.rds$", "", basename(file)))
}

# Somma a `totale` (vettore con nome) i valori di `valori`, per nome.
somma_per_nome <- function(totale, valori) {
  if (length(valori) == 0L) return(totale)
  tutti <- c(totale, valori)
  somme <- rowsum(unname(tutti), names(tutti), reorder = FALSE)
  stats::setNames(somme[, 1L], rownames(somme))
}

esegui_fase(function(parametri, cartella) {
  richiedi_pacchetti("dada2")

  campioni <- unlist(parametri$campioni)
  modelli <- unlist(parametri$modelli)[names(campioni)]
  senza_letture <- as.character(unlist(parametri$senza_letture))
  pool <- parametri$pool
  omega_a <- as.numeric(parametri$omega_a)
  lotto <- as.integer(parametri$lotto)
  processi <- as.integer(parametri$processi)
  if (anyNA(modelli)) {
    stop("campioni senza modello d'errore: ",
         paste(names(campioni)[is.na(modelli)], collapse = ", "))
  }

  errori <- lapply(stats::setNames(nm = sort(unique(modelli), method = "radix")), readRDS)
  pseudo <- identical(pool, "pseudo") && length(campioni) > 1L
  riunito <- isTRUE(pool)

  # dada su un gruppo di campioni, un modello alla volta. Restituisce la
  # lista degli oggetti dada, con il nome dei campioni, nell'ordine dato.
  inferisci <- function(accession, err_per_modello, priori = character(0)) {
    risultati <- list()
    for (file_modello in unique(modelli[accession])) {
      gruppo <- accession[modelli[accession] == file_modello]
      esito <- dada2::dada(
        stats::setNames(campioni[gruppo], gruppo),
        err = err_per_modello[[file_modello]],
        pool = riunito,
        priors = priori,
        multithread = processi,
        verbose = 0,
        OMEGA_A = omega_a
      )
      if (length(gruppo) == 1L) esito <- stats::setNames(list(esito), gruppo)
      risultati[gruppo] <- esito[gruppo]
    }
    risultati[accession]
  }

  lotti <- if (riunito) {
    list(names(campioni))
  } else {
    split(names(campioni), ceiling(seq_along(campioni) / lotto))
  }

  # --- prima passata, o unica ------------------------------------------------
  prevalenza <- numeric(0)
  abbondanza <- numeric(0)
  transizioni <- list()
  varianti <- list()
  for (numero in seq_along(lotti)) {
    message(sprintf("prima passata: lotto %d di %d", numero, length(lotti)))
    esito <- inferisci(lotti[[numero]], errori)
    for (a in names(esito)) {
      unici <- dada2::getUniques(esito[[a]]$clustering)
      if (pseudo) {
        abbondanza <- somma_per_nome(abbondanza, unici)
        prevalenza <- somma_per_nome(prevalenza, (unici > 0) + 0)
        # La funzione con cui dada somma le transizioni dei suoi campioni.
        m <- modelli[[a]]
        transizioni[[m]] <- dada2:::accumulateTrans(
          c(if (!is.null(transizioni[[m]])) list(transizioni[[m]]),
            list(esito[[a]]$trans))
        )
      }
      varianti[[a]] <- unici
    }
    rm(esito)
  }
  varianti_prima <- vapply(varianti, length, integer(1))

  artefatti <- character()
  riepilogo <- list(
    pool = if (is.character(pool)) pool else isTRUE(pool),
    omega_a = omega_a,
    campioni = length(campioni),
    campioni_senza_letture = as.list(sort(senza_letture, method = "radix"))
  )

  # --- seconda passata --------------------------------------------------------
  if (pseudo) {
    soglia_prevalenza <- dada2::getDadaOpt("PSEUDO_PREVALENCE")
    soglia_abbondanza <- dada2::getDadaOpt("PSEUDO_ABUNDANCE")
    abbondanza <- abbondanza[names(prevalenza)]
    scelte <- prevalenza >= soglia_prevalenza | abbondanza >= soglia_abbondanza
    priori <- names(prevalenza)[scelte]

    ordine <- order(-abbondanza[priori], priori, method = "radix")
    priori <- priori[ordine]
    scrivi_atomico(
      c("sequenza\tprevalenza\tabbondanza",
        sprintf("%s\t%.0f\t%.0f", priori, prevalenza[priori], abbondanza[priori])),
      file.path(cartella, "priori.tsv")
    )
    artefatti <- c(artefatti, "priori.tsv")

    # Il modello della seconda passata, come in dada: errorEstimationFunction
    # (loessErrfun) sulle transizioni accumulate della prima.
    seconda <- list()
    for (file_modello in names(errori)) {
      stimato <- tryCatch(
        suppressWarnings(dada2::loessErrfun(transizioni[[file_modello]])),
        error = function(e) NULL
      )
      if (is.null(stimato)) {
        stop("il modello d'errore della seconda passata non si stima dalle ",
             "transizioni della prima per ", basename(file_modello),
             ": dada(pool = \"pseudo\") si fermerebbe allo stesso punto")
      }
      seconda[[file_modello]] <- stimato
      nome <- sprintf("errori_seconda_passata_%s.rds", nome_modello(file_modello))
      salva_rds(stimato, file.path(cartella, nome))
      artefatti <- c(artefatti, nome)
    }
    rm(transizioni)

    varianti <- list()
    for (numero in seq_along(lotti)) {
      message(sprintf("seconda passata: lotto %d di %d", numero, length(lotti)))
      esito <- inferisci(lotti[[numero]], seconda, priori)
      for (a in names(esito)) varianti[[a]] <- dada2::getUniques(esito[[a]]$clustering)
      rm(esito)
    }

    riepilogo$priori <- list(
      regola = "prevalenza >= PSEUDO_PREVALENCE oppure abbondanza >= PSEUDO_ABUNDANCE",
      PSEUDO_PREVALENCE = soglia_prevalenza,
      PSEUDO_ABUNDANCE = format(soglia_abbondanza),
      sequenze_prima_passata = length(prevalenza),
      priori = length(priori)
    )
  }

  varianti <- varianti[sort(names(varianti), method = "radix")]
  salva_rds(varianti, file.path(cartella, "varianti_per_campione.rds"))

  letture <- vapply(varianti, function(v) sum(as.numeric(v)), numeric(1))
  zeri <- stats::setNames(numeric(length(senza_letture)), senza_letture)
  tutte <- c(letture, zeri)
  tutte <- tutte[sort(names(tutte), method = "radix")]
  artefatti <- c(artefatti, traccia_letture(tutte, "denoised", cartella))

  a <- names(varianti)
  scrivi_atomico(
    c("campione\tmodello\tletture\tvarianti_prima_passata\tvarianti",
      sprintf("%s\t%s\t%.0f\t%s\t%d", a, basename(modelli[a]), letture[a],
              if (pseudo) as.character(varianti_prima[a]) else "",
              vapply(varianti, length, integer(1)))),
    file.path(cartella, "varianti.tsv")
  )

  # I livelli espliciti: split() ordinerebbe i modelli secondo la lingua.
  nomi_modelli <- basename(modelli[a])
  per_modello <- lapply(split(a, factor(nomi_modelli, levels = sort(unique(nomi_modelli), method = "radix"))), function(g) {
    list(
      campioni = length(g),
      letture = sum(letture[g]),
      varianti_distinte = length(unique(unlist(lapply(varianti[g], names))))
    )
  })
  riepilogo$modelli <- per_modello
  riepilogo$varianti_distinte <- length(unique(unlist(lapply(varianti, names))))
  scrivi_json(riepilogo, file.path(cartella, "inferenza.json"))
  message(sprintf("memoria di picco (lotti di %d campioni): %s", lotto, memoria_di_picco()))

  # letture_denoised.tsv e' gia' in artefatti, restituito da traccia_letture.
  c(artefatti, "varianti_per_campione.rds", "varianti.tsv", "inferenza.json")
})
