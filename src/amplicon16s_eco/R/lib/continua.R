# Analisi di una variabile continua: correlazione con l'alfa diversita',
# PERMANOVA con un termine numerico, grafico di dispersione.
#
# Richiede uscite.R, campioni.R, alfa.R (fissa_seme) e test.R (.rango,
# .fattori, INDICI, TITOLO_INDICE). Interfaccia:
#   piano_continue(d, parametri)       per ogni variabile continua i campioni di
#                                      ogni analisi e se si esegue, con gli
#                                      avvisi: non calcola nulla, e serve anche
#                                      alla validazione
#   correlazione_alfa(indici, piano)   Spearman per indice, con
#                                      Benjamini-Hochberg sulla famiglia
#   grafico_continua(...)              dispersione per indice
#   permanova_continua(...)            adonis2 per termini, variabile numerica
#   scrivi_continue(...)               le due tabelle
#
# Razionale biologico: una misura (una distanza, un tempo, una concentrazione)
# divisa in gruppi perde l'ordine dei suoi valori, e con valori quasi tutti
# diversi i gruppi sono campioni singoli. La correlazione per ranghi non
# presume una relazione lineare ne' una distribuzione degli indici. Nella
# PERMANOVA la variabile entra come un termine numerico, dopo le variabili
# tecniche. Il test delle dispersioni confronta gruppi, e una variabile continua
# non ne ha: non si esegue, e la tabella lo dichiara.

# Il minimo di campioni per una correlazione: il piu' piccolo n per cui una
# correlazione perfetta scende sotto il livello di significativita'. Il p
# esatto minimo a due code di Spearman e' 2 / n!.
minimo_per_la_correlazione <- function(livello) {
  n <- 2L
  while (2 / factorial(n) >= livello) n <- n + 1L
  n
}

# Con due soli valori distinti la variabile separa due gruppi: e' un confronto
# fra gruppi, non una correlazione.
VALORI_DISTINTI_MINIMI <- 3L

piano_continue <- function(d, parametri) {
  avvisi <- list()
  avviso <- function(codice, dettaglio) {
    avvisi[[length(avvisi) + 1L]] <<- .avviso(codice, dettaglio)
  }
  piani <- list()
  if (!length(d$continue)) return(list(piani = piani, avvisi = avvisi))
  livello <- parametri$stat$significance_level
  minimo <- minimo_per_la_correlazione(livello)
  strati <- parametri$stat$permanova_strata
  marcatori <- unlist(parametri$valori_mancanti)
  valori_strati <- NULL
  if (!is.null(strati)) {
    valori_strati <- stats::setNames(d$tabella[d$campioni, strati], d$campioni)
    valori_strati[e_mancante(valori_strati, marcatori)] <- NA
  }
  rarefatti <- setdiff(d$campioni, d$sotto_profondita)

  giudica <- function(nome, analisi, valori) {
    distinti <- length(unique(valori))
    eseguibile <- length(valori) >= minimo && distinti >= VALORI_DISTINTI_MINIMI
    if (!eseguibile) {
      avviso("E-ECO-31", sprintf(paste0(
        "Variabile continua '%s', %s: %d campioni con un valore (ne servono almeno %d ",
        "perche' una correlazione perfetta scenda sotto %s) e %d valori distinti (ne ",
        "servono almeno %d)."),
        nome, analisi, length(valori), minimo, format(livello), distinti,
        VALORI_DISTINTI_MINIMI))
    }
    eseguibile
  }

  for (v in d$continue) {
    nome <- v$nome
    # --- Alfa diversita': i campioni rarefatti con un valore. -----------------
    x <- v$valori[rarefatti]
    campioni_alfa <- names(x)[!is.na(x)]
    alfa_eseguibile <- giudica(nome, "correlazione con l'alfa diversita'", x[campioni_alfa])

    # --- PERMANOVA: servono anche le variabili tecniche e gli strati. ---------
    completi <- names(v$valori)[!is.na(v$valori)]
    for (tecnica in d$tecniche) {
      senza <- completi[is.na(d$tecnici[[tecnica]][completi])]
      if (length(senza)) {
        avviso("E-ECO-13", sprintf(
          "Variabile continua '%s', PERMANOVA: %d campioni senza valore nella variabile tecnica '%s' (%s).",
          nome, length(senza), tecnica, .elenco(senza)))
        completi <- setdiff(completi, senza)
      }
    }
    if (!is.null(valori_strati)) {
      senza <- completi[is.na(valori_strati[completi])]
      if (length(senza)) {
        avviso("E-ECO-13", sprintf(
          "Variabile continua '%s', PERMANOVA: %d campioni senza valore negli strati '%s' (%s).",
          nome, length(senza), strati, .elenco(senza)))
        completi <- setdiff(completi, senza)
      }
    }
    completi <- completi[ordine_radix(completi)]
    permanova_eseguibile <- giudica(nome, "PERMANOVA", v$valori[completi])
    if (permanova_eseguibile && !is.null(valori_strati)) {
      # Le permutazioni entro gli strati non spostano una variabile che dentro
      # ogni strato ha un solo valore: ogni permutazione da' la stessa
      # statistica, e p varrebbe 1 qualunque sia l'effetto.
      entro <- tapply(v$valori[completi], as.character(valori_strati[completi]),
                      function(x) length(unique(x)))
      if (all(entro == 1L)) {
        permanova_eseguibile <- FALSE
        avviso("E-ECO-31", sprintf(paste0(
          "Variabile continua '%s', PERMANOVA: la variabile ha un solo valore dentro ogni ",
          "strato di '%s', e le permutazioni entro gli strati non la spostano."),
          nome, strati))
      }
    }
    termini <- character()
    stimabile <- permanova_eseguibile
    saturo <- FALSE
    if (permanova_eseguibile) {
      for (tecnica in d$tecniche) {
        livelli <- unique(d$tecnici[[tecnica]][completi])
        if (length(livelli) < 2L) {
          avviso("E-ECO-22", sprintf(
            "Variabile continua '%s': la variabile tecnica '%s' ha il solo livello '%s' fra i %d campioni del test.",
            nome, tecnica, livelli, length(completi)))
          next
        }
        termini <- c(termini, tecnica)
      }
      tecnici <- .fattori(lapply(d$tecnici[termini], function(t) t[completi]))
      modello <- if (length(termini)) cbind(tecnici, .variabile = v$valori[completi]) else
        data.frame(.variabile = v$valori[completi])
      if (length(termini) && .rango(modello) == .rango(tecnici)) {
        stimabile <- FALSE
        avviso("E-ECO-23", sprintf(
          "Variabile continua '%s': e' determinata per intero dalle variabili tecniche (%s) e non le resta alcun grado di liberta'.",
          nome, paste(termini, collapse = ", ")))
      }
      if (.rango(modello) >= length(completi)) {
        saturo <- TRUE
        avviso("E-ECO-26", sprintf(
          "Variabile continua '%s': i termini del modello (%s) assorbono tutti i %d campioni del test.",
          nome, paste(c(termini, nome), collapse = ", "), length(completi)))
      }
    }
    piani[[nome]] <- list(
      nome = nome,
      alfa = list(campioni = campioni_alfa, valori = x[campioni_alfa],
                  eseguibile = alfa_eseguibile),
      permanova = list(campioni = completi, valori = v$valori[completi], termini = termini,
                       stimabile = stimabile, saturo = saturo,
                       eseguibile = permanova_eseguibile)
    )
  }
  list(piani = piani, avvisi = avvisi, valori_strati = valori_strati)
}

# --- Alfa diversita' -------------------------------------------------------- #

# Fino a questo numero di campioni l'algoritmo di stats::cor.test (AS 89)
# enumera le permutazioni; oltre usa una serie di Edgeworth.
CAMPIONI_ENUMERATI_SPEARMAN <- 9L

.spearman <- function(x, y) {
  # La regola di stats::cor.test (exact = NULL), applicata in modo esplicito:
  # senza valori ripetuti (e con meno di 1290 campioni, il limite
  # dell'algoritmo) il p viene dall'algoritmo AS 89, esatto fino a nove
  # campioni e per serie oltre; con valori ripetuti dall'approssimazione t.
  senza_ripetuti <- length(x) < 1290L && !anyDuplicated(x) && !anyDuplicated(y)
  # Un indice costante non ha ranghi da correlare.
  if (length(unique(x)) < 2L || length(unique(y)) < 2L) {
    return(list(test = "Spearman, non calcolabile: valori costanti", rho = NA_real_,
                p = NA_real_))
  }
  esito <- stats::cor.test(x, y, method = "spearman", alternative = "two.sided",
                           exact = senza_ripetuti, continuity = FALSE)
  p <- esito$p.value
  if (senza_ripetuti) {
    # Nessun p esatto a due code sta sotto 2 / n!: la serie, per una
    # correlazione quasi perfetta, puo' dare un valore piu' piccolo, fino a zero.
    p <- max(p, 2 / factorial(length(x)))
  }
  list(
    test = if (!senza_ripetuti) "Spearman, approssimazione t"
           else if (length(x) <= CAMPIONI_ENUMERATI_SPEARMAN) "Spearman, esatto"
           else "Spearman, serie di Edgeworth",
    rho = unname(esito$estimate), p = p)
}

correlazione_alfa <- function(indici, piano) {
  a <- piano$alfa
  if (!a$eseguibile) return(NULL)
  per_campione <- indici[match(a$campioni, indici$campione), , drop = FALSE]
  righe <- lapply(INDICI, function(indice) {
    esito <- .spearman(as.numeric(a$valori), as.numeric(per_campione[[indice]]))
    data.frame(variabile = piano$nome, indice = indice, campioni = length(a$campioni),
               test = esito$test, rho = esito$rho, p = esito$p, stringsAsFactors = FALSE)
  })
  tabella <- do.call(rbind, righe)
  tabella$p[is.nan(tabella$p)] <- NA
  # Una sola famiglia per variabile: i tre indici.
  tabella$p_bh <- stats::p.adjust(tabella$p, method = "BH")
  tabella$dimensione_famiglia <- sum(!is.na(tabella$p))
  tabella
}

# Dispersione di ogni indice contro la variabile, un pannello per indice, sui
# campioni della correlazione. Nessun numero casuale.
grafico_continua <- function(indici, piano, profondita, percorso) {
  a <- piano$alfa
  per_campione <- indici[match(a$campioni, indici$campione), , drop = FALSE]
  x <- as.numeric(a$valori)
  apri_png(percorso)
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::par(mfrow = c(1L, length(INDICI)), oma = c(2.5, 0, 2.5, 0), mar = c(4.5, 4.5, 2, 1))
  graphics::par(cex = 0.7)
  colore <- colori_gruppi(1L)
  for (indice in INDICI) {
    graphics::plot(x, as.numeric(per_campione[[indice]]), pch = 19, cex = 0.8, col = colore,
                   xlab = piano$nome, ylab = TITOLO_INDICE[[indice]],
                   main = TITOLO_INDICE[[indice]])
  }
  graphics::mtext(sprintf("Alfa diversità contro %s (%d campioni)", piano$nome,
                          length(a$campioni)),
                  side = 3, outer = TRUE, line = 0.8, cex = 0.9, font = 2)
  graphics::mtext(sprintf("conteggi rarefatti a %d letture", profondita),
                  side = 1, outer = TRUE, line = 1, cex = 0.7)
}

# --- PERMANOVA -------------------------------------------------------------- #

permanova_continua <- function(matrice, distanza, piano, tecnici, valori_strati,
                               permutazioni, seme) {
  p <- piano$permanova
  if (!p$eseguibile) return(NULL)
  campioni <- p$campioni
  d <- stats::as.dist(matrice[campioni, campioni, drop = FALSE])
  dati <- data.frame(row.names = campioni)
  nomi <- character()
  for (i in seq_along(p$termini)) {
    dati[[paste0("t", i)]] <- factor(tecnici[[p$termini[i]]][campioni])
    nomi[paste0("t", i)] <- p$termini[i]
  }
  ruolo <- stats::setNames(rep("tecnica", length(nomi)), names(nomi))
  if (p$stimabile) {
    dati$v <- as.numeric(p$valori)
    nomi["v"] <- piano$nome
    ruolo["v"] <- "analizzata"
  }
  usate <- NA_integer_
  righe <- list()
  riga <- function(analisi, termine, ruolo, gl, sq, r2, f, pv, nota = "") {
    righe[[length(righe) + 1L]] <<- data.frame(
      variabile = piano$nome, distanza = distanza, analisi = analisi, termine = termine,
      ruolo = ruolo, gl = as.integer(gl), somma_quadrati = sq, r2 = r2, f = f, p = pv,
      permutazioni = usate, campioni = length(campioni), nota = nota,
      stringsAsFactors = FALSE)
  }
  if (p$saturo || !length(nomi)) {
    riga("permanova", piano$nome, "analizzata", NA, NA, NA, NA, NA,
         if (p$saturo) "non eseguita: il modello non lascia gradi di liberta' residui"
         else "non stimabile: determinata per intero dalle variabili tecniche")
  } else {
    controllo <- if (is.null(valori_strati)) {
      permute::how(nperm = permutazioni)
    } else {
      permute::how(nperm = permutazioni, blocks = factor(valori_strati[campioni]))
    }
    formula <- stats::as.formula(paste("d ~", paste(names(nomi), collapse = " + ")))
    fissa_seme(seme)
    grezzo <- vegan::adonis2(
      formula, data = dati, permutations = controllo, by = "terms", parallel = 1,
      sqrt.dist = FALSE, add = FALSE, na.action = stats::na.fail)
    usate <- nrow(attr(grezzo, "F.perm"))
    esito <- as.data.frame(grezzo)
    # Con una distanza non euclidea i termini possono spiegare insieme piu'
    # della variabilita' totale: senza un residuo positivo F non ha significato.
    valutabile <- isTRUE(esito["Residual", "SumOfSqs"] > 0)
    for (termine in names(nomi)) {
      if (termine %in% rownames(esito)) {
        e <- esito[termine, ]
        riga("permanova", nomi[[termine]], ruolo[[termine]], e$Df, e$SumOfSqs, e$R2,
             if (valutabile) e$F else NA, if (valutabile) e[["Pr(>F)"]] else NA,
             if (valutabile) "" else paste0(
               "non valutabile: somma dei quadrati residua non positiva; con questa ",
               "distanza, non euclidea, i termini spiegano piu' della variabilita' totale"))
      } else {
        riga("permanova", nomi[[termine]], ruolo[[termine]], 0L, NA, NA, NA, NA,
             "non stimabile: nessun grado di liberta' dopo i termini precedenti")
      }
    }
    if (!p$stimabile) {
      riga("permanova", piano$nome, "analizzata", 0L, NA, NA, NA, NA,
           "non stimabile: determinata per intero dalle variabili tecniche")
    }
    riga("permanova", "residuo", "residuo", esito["Residual", "Df"],
         esito["Residual", "SumOfSqs"], esito["Residual", "R2"], NA, NA)
    riga("permanova", "totale", "totale", esito["Total", "Df"],
         esito["Total", "SumOfSqs"], esito["Total", "R2"], NA, NA)
  }
  riga("dispersioni", piano$nome, "analizzata", NA, NA, NA, NA, NA,
       "non eseguito: una variabile continua non ha gruppi di cui confrontare la dispersione")
  do.call(rbind, righe)
}

scrivi_continue <- function(alfa, multivariati, parametri, profondita, cartella) {
  dir.create(cartella, showWarnings = FALSE)
  vuota_alfa <- data.frame(
    variabile = character(), indice = character(), campioni = integer(), test = character(),
    rho = character(), p = character(), p_bh = character(), dimensione_famiglia = integer(),
    stringsAsFactors = FALSE)
  if (!is.null(alfa)) {
    alfa$rho <- reale(alfa$rho)
    for (colonna in c("p", "p_bh")) alfa[[colonna]] <- esponenziale(alfa[[colonna]])
  }
  scrivi_tsv(
    if (is.null(alfa)) vuota_alfa else alfa, file.path(cartella, "alfa_continua.tsv"),
    intestazione = c(
      "Correlazione dell'alfa diversita' con ogni variabile continua (design.continuous_variables).",
      sprintf("indici calcolati sui conteggi rarefatti a %d letture; campioni rarefatti con un valore della variabile",
              profondita),
      "correlazione per ranghi di Spearman, a due code; rho: coefficiente di correlazione",
      sprintf("senza valori ripetuti nell'indice e nella variabile il p e' esatto fino a %d campioni e dalla serie di Edgeworth oltre (algoritmo AS 89 di stats::cor.test), mai sotto il minimo esatto 2/n!; con valori ripetuti approssimazione t",
              CAMPIONI_ENUMERATI_SPEARMAN),
      "p_bh: p corretta con Benjamini-Hochberg sulla famiglia dei tre indici della variabile",
      "dimensione_famiglia: numero di test della famiglia; p e p_bh in notazione esponenziale a sei decimali",
      sprintf("una correlazione non si esegue con meno di %d campioni (il minimo perche' una correlazione perfetta scenda sotto %s) o con meno di %d valori distinti",
              minimo_per_la_correlazione(parametri$stat$significance_level),
              format(parametri$stat$significance_level), VALORI_DISTINTI_MINIMI)
    )
  )
  vuota <- data.frame(
    variabile = character(), distanza = character(), analisi = character(),
    termine = character(), ruolo = character(), gl = integer(),
    somma_quadrati = character(), r2 = character(), f = character(), p = character(),
    permutazioni = integer(), campioni = integer(), nota = character(),
    stringsAsFactors = FALSE)
  if (!is.null(multivariati)) {
    for (colonna in c("somma_quadrati", "r2", "f", "p")) {
      multivariati[[colonna]] <- reale(multivariati[[colonna]])
    }
  }
  strati <- parametri$stat$permanova_strata
  scrivi_tsv(
    if (is.null(multivariati)) vuota else multivariati,
    file.path(cartella, "permanova_continua.tsv"),
    intestazione = c(
      "PERMANOVA con ogni variabile continua, per distanza.",
      "permanova: vegan::adonis2, by = \"terms\" (somme dei quadrati sequenziali): le variabili tecniche sono i primi termini, nell'ordine dichiarato, poi la variabile come termine numerico (un grado di liberta')",
      "dispersioni: non eseguito. Il test dell'omogeneita' delle dispersioni confronta gruppi, e una variabile continua non ne ha",
      "con una distanza non euclidea (Bray-Curtis, Jaccard, UniFrac) la somma dei quadrati residua puo' non essere positiva: F e p non sono allora valutabili, e la nota lo dice",
      sprintf("permutazioni richieste: %d; la colonna permutazioni riporta quelle usate; strati entro cui si permuta: %s",
              as.integer(parametri$stat$permanova_permutations),
              if (is.null(strati)) "nessuno" else strati),
      "ruolo: tecnica, analizzata, residuo, totale; gl: gradi di liberta'"
    )
  )
  c("test/alfa_continua.tsv", "test/permanova_continua.tsv")
}
