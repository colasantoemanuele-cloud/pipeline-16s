# Test fra gruppi: alfa diversita', PERMANOVA e omogeneita' delle dispersioni.
#
# Richiede uscite.R, campioni.R e alfa.R (fissa_seme). Interfaccia:
#   disegno_dei_test(d, parametri)   per ogni variabile, i campioni e i termini
#                                    di ogni test, con gli avvisi: non calcola
#                                    nulla, e serve anche alla validazione
#   test_alfa(indici, piano)         Wilcoxon o Kruskal-Wallis, con
#                                    Benjamini-Hochberg sulla famiglia
#   grafico_alfa(...)                scatole per gruppo, un pannello per indice,
#                                    sugli stessi campioni del confronto
#   permanova(matrice, piano, ...)   adonis2 per termini e betadisper
#   scrivi_test(...)                 le due tabelle
#
# Razionale biologico: una PERMANOVA significativa non distingue una differenza
# di posizione (composizione media) da una di dispersione (variabilita' fra
# campioni dello stesso gruppo), quindi ogni PERMANOVA e' accompagnata dal test
# delle dispersioni sugli stessi campioni e sulla stessa matrice. Le variabili
# tecniche entrano come primi termini con somme dei quadrati sequenziali:
# l'effetto della variabile analizzata e' stimato al netto del lotto, e la
# varianza che lotto e variabile condividono va al lotto. Per questo la
# confusione fra le due si dichiara. Con gli strati i due test permutano allo
# stesso modo, entro gli strati: rispondono alla stessa ipotesi di
# scambiabilita' dei campioni.

.rango <- function(tabella) {
  if (!ncol(tabella)) return(1L)
  qr(stats::model.matrix(~ ., data = tabella))$rank
}

.fattori <- function(tabella) {
  as.data.frame(lapply(tabella, factor), stringsAsFactors = FALSE)
}

# Logaritmo del numero di disposizioni distinte dei campioni fra le
# combinazioni dei termini del modello, entro i blocchi: per ogni blocco
# n! / prod(n_c!), con n_c i campioni di ogni combinazione. Permutare due
# campioni con la stessa combinazione non cambia il modello, quindi e' questo
# numero, non n!, a fissare la p minima raggiungibile.
.log_disposizioni <- function(combinazioni, blocchi) {
  sum(vapply(split(combinazioni, blocchi), function(x) {
    lgamma(length(x) + 1) - sum(lgamma(table(x) + 1))
  }, numeric(1)))
}

# Gruppi della variabile che stanno per intero in un solo livello di `altra`.
.annidati <- function(valori, altra) {
  livelli <- tapply(altra, valori, function(x) length(unique(x)))
  nomi <- names(livelli)[livelli == 1L]
  nomi[ordine_radix(nomi)]
}

disegno_dei_test <- function(d, parametri) {
  avvisi <- list()
  avviso <- function(codice, dettaglio) {
    avvisi[[length(avvisi) + 1L]] <<- .avviso(codice, dettaglio)
  }
  minimo <- as.integer(parametri$stat$min_group_size)
  permutazioni <- as.integer(parametri$stat$permanova_permutations)
  strati <- parametri$stat$permanova_strata
  marcatori <- unlist(parametri$valori_mancanti)
  rarefatti <- setdiff(d$campioni, d$sotto_profondita)

  if (!length(d$tecniche)) {
    avviso("E-ECO-20", paste0(
      "design.technical_variables e' vuoto: nella PERMANOVA l'effetto di ogni ",
      "variabile non e' stimato al netto di un lotto."))
  }
  valori_strati <- NULL
  if (!is.null(strati)) {
    valori_strati <- stats::setNames(d$tabella[d$campioni, strati], d$campioni)
    valori_strati[e_mancante(valori_strati, marcatori)] <- NA
  }

  piani <- list()
  for (v in d$variabili) {
    nome <- v$nome
    # --- Alfa diversita': i campioni rarefatti con un valore. -----------------
    valori_alfa <- v$valori[rarefatti]
    g_alfa <- gruppi_di(valori_alfa, minimo)
    if (!setequal(g_alfa$nei_test, v$nei_test)) {
      avviso(if (length(g_alfa$nei_test) < 2L) "E-ECO-12" else "E-ECO-11", sprintf(
        "Variabile '%s', test dell'alfa diversita': fra i campioni rarefatti i gruppi con almeno %d campioni sono %d (%s).",
        nome, minimo, length(g_alfa$nei_test), .elenco(g_alfa$nei_test)))
    }
    campioni_alfa <- names(valori_alfa)[!is.na(valori_alfa) & valori_alfa %in% g_alfa$nei_test]
    stato_alfa <- stats::setNames(rep("analizzato", length(d$campioni)), d$campioni)
    stato_alfa[v$mancanti] <- "escluso: valore mancante"
    stato_alfa[names(valori_alfa)[!is.na(valori_alfa) & valori_alfa %in% g_alfa$piccoli]] <-
      "escluso dai test: gruppo sotto la dimensione minima"

    # --- PERMANOVA: servono anche le variabili tecniche e gli strati. ---------
    stato <- stato_per_variabile(v, d$campioni)
    completi <- names(v$valori)[!is.na(v$valori)]
    for (tecnica in d$tecniche) {
      senza <- completi[is.na(d$tecnici[[tecnica]][completi])]
      if (length(senza)) {
        avviso("E-ECO-13", sprintf(
          "Variabile '%s', PERMANOVA: %d campioni senza valore nella variabile tecnica '%s' (%s).",
          nome, length(senza), tecnica, .elenco(senza)))
        stato[senza] <- sprintf("escluso: valore mancante nella variabile tecnica '%s'", tecnica)
        completi <- setdiff(completi, senza)
      }
    }
    if (!is.null(valori_strati)) {
      senza <- completi[is.na(valori_strati[completi])]
      if (length(senza)) {
        avviso("E-ECO-13", sprintf(
          "Variabile '%s', PERMANOVA: %d campioni senza valore negli strati '%s' (%s).",
          nome, length(senza), strati, .elenco(senza)))
        stato[senza] <- sprintf("escluso: valore mancante negli strati '%s'", strati)
        completi <- setdiff(completi, senza)
      }
    }
    g <- gruppi_di(v$valori[completi], minimo)
    if (!setequal(g$nei_test, v$nei_test)) {
      avviso(if (length(g$nei_test) < 2L) "E-ECO-12" else "E-ECO-11", sprintf(
        "Variabile '%s', PERMANOVA: fra i campioni con tutte le variabili del modello i gruppi con almeno %d campioni sono %d (%s).",
        nome, minimo, length(g$nei_test), .elenco(g$nei_test)))
    }
    piccoli <- completi[v$valori[completi] %in% g$piccoli]
    stato[piccoli] <- "escluso dai test: gruppo sotto la dimensione minima"
    campioni <- completi[v$valori[completi] %in% g$nei_test]
    campioni <- campioni[ordine_radix(campioni)]
    eseguibile <- length(g$nei_test) >= 2L

    termini <- character()
    stimabile <- eseguibile
    saturo <- FALSE
    possibili <- NA_real_
    if (eseguibile) {
      valori <- v$valori[campioni]
      for (tecnica in d$tecniche) {
        livelli <- unique(d$tecnici[[tecnica]][campioni])
        if (length(livelli) < 2L) {
          avviso("E-ECO-22", sprintf(
            "Variabile '%s': la variabile tecnica '%s' ha il solo livello '%s' fra i %d campioni del test.",
            nome, tecnica, livelli, length(campioni)))
          next
        }
        # Un termine tecnico gia' determinato da quelli che lo precedono (una
        # corsa che raccoglie piastre intere, dopo la piastra) resta nel
        # modello senza gradi di liberta': lo si dichiara.
        prima <- .fattori(lapply(d$tecnici[termini], function(x) x[campioni]))
        dopo <- .fattori(lapply(d$tecnici[c(termini, tecnica)], function(x) x[campioni]))
        if (length(termini) && .rango(dopo) == .rango(prima)) {
          avviso("E-ECO-23", sprintf(
            "Variabile '%s': la variabile tecnica '%s' e' determinata per intero dai termini che la precedono (%s) e non le resta alcun grado di liberta'.",
            nome, tecnica, paste(termini, collapse = ", ")))
        }
        termini <- c(termini, tecnica)
        annidati <- .annidati(valori, d$tecnici[[tecnica]][campioni])
        if (length(annidati)) {
          avviso("E-ECO-18", sprintf(
            "Variabile '%s': %d gruppi su %d stanno per intero in un solo livello della variabile tecnica '%s' (%s); le somme sequenziali attribuiscono quella varianza a '%s'.",
            nome, length(annidati), length(g$nei_test), tecnica, .elenco(annidati), tecnica))
        }
      }
      if (length(termini)) {
        tecnici <- .fattori(lapply(d$tecnici[termini], function(x) x[campioni]))
        con <- cbind(tecnici, .variabile = factor(valori))
        if (.rango(con) == .rango(tecnici)) {
          stimabile <- FALSE
          avviso("E-ECO-23", sprintf(
            "Variabile '%s': e' determinata per intero dalle variabili tecniche (%s) e non le resta alcun grado di liberta'.",
            nome, paste(termini, collapse = ", ")))
        }
      }
      modello <- .fattori(c(lapply(d$tecnici[termini], function(x) x[campioni]),
                            list(.variabile = valori)))
      if (.rango(modello) >= length(campioni)) {
        saturo <- TRUE
        avviso("E-ECO-26", sprintf(
          "Variabile '%s': i termini del modello (%s) assorbono tutti i %d campioni del test.",
          nome, paste(c(termini, nome), collapse = ", "), length(campioni)))
      }
      blocchi <- rep("tutti", length(campioni))
      if (!is.null(valori_strati)) {
        blocchi <- valori_strati[campioni]
        # Con un solo strato le permutazioni sono libere: nulla da dichiarare.
        annidati <- if (length(unique(blocchi)) > 1L) .annidati(valori, blocchi) else character()
        if (length(annidati)) {
          avviso("E-ECO-18", sprintf(
            "Variabile '%s': %d gruppi su %d stanno per intero in un solo strato di '%s' (%s); le permutazioni entro gli strati non li spostano.",
            nome, length(annidati), length(g$nei_test), strati, .elenco(annidati)))
        }
      }
      combinazioni <- do.call(paste, c(lapply(modello, as.character), sep = "\r"))
      logaritmo <- .log_disposizioni(combinazioni, blocchi)
      possibili <- if (logaritmo < log(1e15)) round(exp(logaritmo)) else Inf
      if (!saturo && possibili <= permutazioni) {
        avviso("E-ECO-21", sprintf(
          "Variabile '%s': con %d campioni le disposizioni distinte fra le combinazioni del modello sono %.0f, a fronte di %d permutazioni richieste; la p minima raggiungibile e' %s.",
          nome, length(campioni), possibili, permutazioni, reale(1 / possibili)))
      }
    }
    piani[[nome]] <- list(
      nome = nome,
      alfa = list(campioni = campioni_alfa, valori = valori_alfa[campioni_alfa],
                  gruppi = g_alfa$nei_test, stato = stato_alfa,
                  gruppi_esclusi = length(g_alfa$piccoli),
                  campioni_esclusi = sum(valori_alfa %in% g_alfa$piccoli)),
      permanova = list(campioni = campioni, valori = v$valori[campioni],
                       gruppi = g$nei_test, dimensioni = g$dimensioni[g$nei_test],
                       termini = termini, stimabile = stimabile, saturo = saturo,
                       eseguibile = eseguibile, stato = stato,
                       disposizioni = possibili)
    )
  }
  list(piani = piani, avvisi = avvisi, strati = strati, valori_strati = valori_strati)
}

# --- Alfa diversita' ------------------------------------------------------- #

.wilcoxon <- function(x, y) {
  # La regola di stats::wilcox.test (exact = NULL), applicata in modo
  # esplicito: esatto con meno di 50 campioni per gruppo e senza valori
  # ripetuti, altrimenti approssimazione normale con correzione di continuita'.
  esatto <- length(x) < 50L && length(y) < 50L && !anyDuplicated(c(x, y))
  esito <- stats::wilcox.test(x, y, alternative = "two.sided", paired = FALSE,
                              exact = esatto, correct = TRUE, conf.int = FALSE)
  list(
    test = if (esatto) "Wilcoxon della somma dei ranghi, esatto"
           else "Wilcoxon della somma dei ranghi, approssimazione normale con correzione di continuita'",
    statistica = unname(esito$statistic), gl = NA_integer_, p = esito$p.value
  )
}

INDICI <- c("ricchezza_osservata", "shannon", "gini_simpson")

test_alfa <- function(indici, piano) {
  campioni <- piano$alfa$campioni
  gruppi <- piano$alfa$gruppi
  if (length(gruppi) < 2L) return(NULL)
  valori <- piano$alfa$valori
  righe <- list()
  riga <- function(indice, confronto, a, b, n, esito) {
    righe[[length(righe) + 1L]] <<- data.frame(
      variabile = piano$nome, indice = indice, confronto = confronto,
      gruppo_1 = a, gruppo_2 = b, n_1 = n[1L], n_2 = n[2L], test = esito$test,
      statistica = esito$statistica, gl = esito$gl, p = esito$p,
      stringsAsFactors = FALSE)
  }
  per_campione <- indici[match(campioni, indici$campione), , drop = FALSE]
  for (indice in INDICI) {
    x <- as.numeric(per_campione[[indice]])
    if (length(gruppi) == 2L) {
      a <- x[valori == gruppi[1L]]
      b <- x[valori == gruppi[2L]]
      riga(indice, paste(gruppi, collapse = " vs "), gruppi[1L], gruppi[2L],
           c(length(a), length(b)), .wilcoxon(a, b))
    } else {
      globale <- stats::kruskal.test(x, factor(valori, levels = gruppi))
      riga(indice, "globale", NA_character_, NA_character_, c(length(x), NA_integer_),
           list(test = "Kruskal-Wallis", statistica = unname(globale$statistic),
                gl = as.integer(globale$parameter), p = globale$p.value))
      for (i in seq_len(length(gruppi) - 1L)) {
        for (j in seq(i + 1L, length(gruppi))) {
          a <- x[valori == gruppi[i]]
          b <- x[valori == gruppi[j]]
          riga(indice, paste(gruppi[i], gruppi[j], sep = " vs "), gruppi[i], gruppi[j],
               c(length(a), length(b)), .wilcoxon(a, b))
        }
      }
    }
  }
  tabella <- do.call(rbind, righe)
  tabella$p[is.nan(tabella$p)] <- NA
  # Una sola famiglia per variabile: tutti gli indici per tutti i confronti.
  tabella$p_bh <- stats::p.adjust(tabella$p, method = "BH")
  tabella$dimensione_famiglia <- sum(!is.na(tabella$p))
  tabella
}

TITOLO_INDICE <- c(
  ricchezza_osservata = "ricchezza osservata",
  shannon = "Shannon",
  gini_simpson = "Gini-Simpson"
)

# Scatole per gruppo con i punti dei campioni, un pannello per indice. I
# campioni sono quelli del confronto: i gruppi sotto la dimensione minima non
# compaiono, e la didascalia dice quanti sono. La posizione orizzontale dei
# punti entro la scatola segue l'ordine dei campioni: nessun numero casuale,
# quindi gli stessi byte a ogni esecuzione.
grafico_alfa <- function(indici, piano, profondita, minimo, percorso) {
  a <- piano$alfa
  gruppi <- a$gruppi
  per_campione <- indici[match(a$campioni, indici$campione), , drop = FALSE]
  fattore <- factor(a$valori, levels = gruppi)
  posizione <- as.integer(fattore)
  dimensioni <- as.integer(table(fattore))
  scarto <- numeric(length(posizione))
  for (i in seq_along(gruppi)) {
    quali <- which(posizione == i)
    scarto[quali] <- if (length(quali) > 1L) seq(-0.25, 0.25, length.out = length(quali)) else 0
  }
  colori <- colori_gruppi(length(gruppi))
  etichette <- sprintf("%s (n = %d)", gruppi, dimensioni)
  apri_png(percorso)
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::par(mfrow = c(1L, length(INDICI)), oma = c(2.5, 0, 2.5, 0),
                mar = c(min(1.5 + 0.45 * max(nchar(etichette)), 18), 4.5, 2, 1), cex = 0.7)
  for (indice in INDICI) {
    x <- as.numeric(per_campione[[indice]])
    graphics::boxplot(
      x ~ fattore, outline = FALSE, ylim = range(x), names = etichette, las = 2,
      col = grDevices::adjustcolor(colori, alpha.f = 0.25), border = colori,
      xlab = "", ylab = TITOLO_INDICE[[indice]], main = TITOLO_INDICE[[indice]])
    graphics::points(posizione + scarto, x, pch = 19, cex = 0.6, col = colori[posizione])
  }
  graphics::mtext(sprintf("Alfa diversita' per %s (%d campioni in %d gruppi)", piano$nome,
                          length(a$campioni), length(gruppi)),
                  side = 3, outer = TRUE, line = 0.8, cex = 0.9, font = 2)
  graphics::mtext(sprintf(
    "conteggi rarefatti a %d letture; gruppi con almeno %d campioni; esclusi dal confronto: %d gruppi (%d campioni)",
    profondita, minimo, a$gruppi_esclusi, a$campioni_esclusi),
    side = 1, outer = TRUE, line = 1, cex = 0.7)
}

# --- PERMANOVA e dispersioni ----------------------------------------------- #

permanova <- function(matrice, distanza, piano, tecnici, valori_strati, permutazioni,
                      livello, seme) {
  p <- piano$permanova
  if (!p$eseguibile) return(NULL)
  campioni <- p$campioni
  d <- stats::as.dist(matrice[campioni, campioni, drop = FALSE])
  # Nomi di colonna neutri: quelli dell'oggetto possono contenere spazi.
  dati <- data.frame(row.names = campioni)
  nomi <- character()
  for (i in seq_along(p$termini)) {
    dati[[paste0("t", i)]] <- factor(tecnici[[p$termini[i]]][campioni])
    nomi[paste0("t", i)] <- p$termini[i]
  }
  ruolo <- stats::setNames(rep("tecnica", length(nomi)), names(nomi))
  if (p$stimabile) {
    dati$v <- factor(p$valori, levels = p$gruppi)
    nomi["v"] <- piano$nome
    ruolo["v"] <- "analizzata"
  }
  controllo <- if (is.null(valori_strati)) {
    permute::how(nperm = permutazioni)
  } else {
    permute::how(nperm = permutazioni, blocks = factor(valori_strati[campioni]))
  }
  righe <- list()
  riga <- function(analisi, termine, ruolo, gl, sq, r2, f, pv, nota = "") {
    righe[[length(righe) + 1L]] <<- data.frame(
      variabile = piano$nome, distanza = distanza, analisi = analisi, termine = termine,
      ruolo = ruolo, gl = as.integer(gl), somma_quadrati = sq, r2 = r2, f = f, p = pv,
      permutazioni = usate, campioni = length(campioni),
      gruppi = length(p$gruppi), nota = nota, stringsAsFactors = FALSE)
  }
  usate <- NA_integer_
  if (p$saturo) {
    riga("permanova", piano$nome, "analizzata", NA, NA, NA, NA, NA,
         "non eseguita: il modello non lascia gradi di liberta' residui")
  } else {
    formula <- stats::as.formula(paste("d ~", paste(names(nomi), collapse = " + ")))
    fissa_seme(seme)
    grezzo <- vegan::adonis2(
      formula, data = dati, permutations = controllo, by = "terms", parallel = 1,
      sqrt.dist = FALSE, add = FALSE, na.action = stats::na.fail)
    # Le permutazioni effettivamente usate: meno delle richieste quando le
    # disposizioni possibili sono poche e vengono enumerate tutte.
    usate <- nrow(attr(grezzo, "F.perm"))
    esito <- as.data.frame(grezzo)
    una_sola <- is.finite(p$disposizioni) && p$disposizioni <= 1
    for (termine in names(nomi)) {
      if (termine %in% rownames(esito)) {
        e <- esito[termine, ]
        riga("permanova", nomi[[termine]], ruolo[[termine]], e$Df, e$SumOfSqs, e$R2, e$F,
             if (una_sola) NA else e[["Pr(>F)"]],
             if (una_sola) "p non calcolabile: una sola disposizione possibile dei campioni" else "")
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

  # Dispersioni: distanze dalla mediana spaziale del proprio gruppo, con lo
  # stesso schema di permutazione della PERMANOVA (entro gli strati, se
  # dichiarati).
  # Con una distanza non euclidea betadisper azzera le distanze al quadrato
  # negative e lo segnala: e' il comportamento documentato, dichiarato
  # nell'intestazione della tabella. Con residuo nullo anova lo segnala a sua
  # volta, e qui lo si dichiara con un codice proprio.
  attesi <- "squared distances are negative and changed to zero|essentially perfect fit"
  zitto <- function(espressione) withCallingHandlers(espressione, warning = function(w) {
    if (grepl(attesi, conditionMessage(w))) invokeRestart("muffleWarning")
  })
  modello <- zitto(vegan::betadisper(
    d, factor(p$valori, levels = p$gruppi), type = "median",
    bias.adjust = FALSE, sqrt.dist = FALSE, add = FALSE))
  # Con un solo campione per strato non esiste alcuna permutazione: F resta
  # quello dell'analisi della varianza, la p non e' calcolabile.
  senza_permutazioni <- !is.null(valori_strati) && all(table(valori_strati[campioni]) < 2L)
  if (senza_permutazioni) {
    varianza <- zitto(stats::anova(modello))
    prova <- data.frame(
      Df = varianza[["Df"]], "Sum Sq" = varianza[["Sum Sq"]], F = varianza[["F value"]],
      N.Perm = 0L, "Pr(>F)" = NA_real_, check.names = FALSE,
      row.names = c("Groups", "Residuals"))
  } else {
    fissa_seme(seme)
    prova <- zitto(vegan::permutest(
      modello, pairwise = FALSE, permutations = controllo, parallel = 1))$tab
  }
  usate <- as.integer(prova["Groups", "N.Perm"])
  residuo <- prova["Residuals", "Sum Sq"]
  totale <- residuo + prova["Groups", "Sum Sq"]
  # Senza variabilita' residua F e' un rapporto fra errori di arrotondamento.
  valutabile <- is.finite(totale) && totale > 0 && residuo > 1e-12 * totale &&
    prova["Residuals", "Df"] > 0
  p_dispersioni <- if (valutabile) prova["Groups", "Pr(>F)"] else NA_real_
  diverse <- !is.na(p_dispersioni) && p_dispersioni < livello
  riga("dispersioni", piano$nome, "analizzata", prova["Groups", "Df"],
       prova["Groups", "Sum Sq"], NA, if (valutabile) prova["Groups", "F"] else NA,
       p_dispersioni,
       if (!valutabile) {
         "non valutabile: nessuna variabilita' residua delle distanze dal centro dei gruppi"
       } else if (senza_permutazioni) {
         "p non calcolabile: una sola disposizione possibile dei campioni entro gli strati"
       } else if (diverse) sprintf(
         "dispersioni diverse (p < %s): la PERMANOVA non va letta come differenza di posizione",
         format(livello)) else "")
  riga("dispersioni", "residuo", "residuo", prova["Residuals", "Df"], residuo, NA, NA, NA)
  list(tabella = do.call(rbind, righe), dispersioni_diverse = diverse,
       dispersioni_valutabili = valutabile, p_dispersioni = p_dispersioni)
}

scrivi_test <- function(alfa, multivariati, parametri, profondita, cartella) {
  dir.create(cartella, showWarnings = FALSE)
  vuota_alfa <- data.frame(
    variabile = character(), indice = character(), confronto = character(),
    gruppo_1 = character(), gruppo_2 = character(), n_1 = integer(), n_2 = integer(),
    test = character(), statistica = character(), gl = integer(), p = character(),
    p_bh = character(), dimensione_famiglia = integer(), stringsAsFactors = FALSE)
  if (!is.null(alfa)) {
    alfa$statistica <- reale(alfa$statistica)
    for (colonna in c("p", "p_bh")) alfa[[colonna]] <- esponenziale(alfa[[colonna]])
  }
  scrivi_tsv(
    if (is.null(alfa)) vuota_alfa else alfa, file.path(cartella, "alfa_diversita.tsv"),
    intestazione = c(
      "Confronto dell'alfa diversita' fra i gruppi di ogni variabile.",
      sprintf("indici calcolati sui conteggi rarefatti a %d letture; gruppi con almeno %d campioni",
              profondita, as.integer(parametri$stat$min_group_size)),
      "due gruppi: Wilcoxon della somma dei ranghi; piu' di due: Kruskal-Wallis (confronto 'globale') e Wilcoxon per ogni coppia",
      "Wilcoxon esatto con meno di 50 campioni per gruppo e senza valori ripetuti, altrimenti approssimazione normale con correzione di continuita'",
      "p_bh: p corretta con Benjamini-Hochberg sulla famiglia di tutti i test della variabile (ogni indice per ogni confronto)",
      "dimensione_famiglia: numero di test della famiglia; p e p_bh in notazione esponenziale a sei decimali",
      "un test non calcolabile (indice con lo stesso valore in tutti i campioni) ha p NA e non entra nella famiglia"
    )
  )
  vuota <- data.frame(
    variabile = character(), distanza = character(), analisi = character(),
    termine = character(), ruolo = character(), gl = integer(),
    somma_quadrati = character(), r2 = character(), f = character(), p = character(),
    permutazioni = integer(), campioni = integer(), gruppi = integer(),
    nota = character(), stringsAsFactors = FALSE)
  if (!is.null(multivariati)) {
    for (colonna in c("somma_quadrati", "r2", "f", "p")) {
      multivariati[[colonna]] <- reale(multivariati[[colonna]])
    }
  }
  strati <- parametri$stat$permanova_strata
  scrivi_tsv(
    if (is.null(multivariati)) vuota else multivariati, file.path(cartella, "permanova.tsv"),
    intestazione = c(
      "PERMANOVA e omogeneita' delle dispersioni, per variabile e per distanza.",
      "permanova: vegan::adonis2, by = \"terms\" (somme dei quadrati sequenziali): le variabili tecniche sono i primi termini, nell'ordine dichiarato, poi la variabile analizzata",
      "dispersioni: vegan::betadisper (distanze dalla mediana spaziale del gruppo) con permutest, sugli stessi campioni e sulla stessa matrice",
      sprintf("permutazioni richieste: %d; la colonna permutazioni riporta quelle usate (tutte le disposizioni, se sono poche); strati entro cui permutano PERMANOVA e dispersioni: %s",
              as.integer(parametri$stat$permanova_permutations),
              if (is.null(strati)) "nessuno" else strati),
      "con una distanza non euclidea betadisper azzera le distanze al quadrato negative dal centro del gruppo",
      sprintf("livello di significativita' per la nota sulle dispersioni: %s",
              format(parametri$stat$significance_level)),
      "ruolo: tecnica, analizzata, residuo, totale; gl: gradi di liberta'; r2 e f solo per la PERMANOVA"
    )
  )
  c("test/alfa_diversita.tsv", "test/permanova.tsv")
}
