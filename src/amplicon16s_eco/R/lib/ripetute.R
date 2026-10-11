# Disegno a misure ripetute: i test entro soggetto per ogni variabile
# categorica, quando la colonna del soggetto e' dichiarata (design.subject).
#
# Richiede uscite.R, campioni.R, alfa.R (fissa_seme) e test.R (INDICI,
# permanova). Interfaccia:
#   piano_appaiato(d, parametri)       per ogni variabile i soggetti completi di
#                                      ogni analisi e se si esegue, con gli
#                                      avvisi: non calcola nulla, e serve anche
#                                      alla validazione
#   test_alfa_appaiato(indici, piano)  Wilcoxon dei ranghi con segno o Friedman,
#                                      con Benjamini-Hochberg sulla famiglia
#   permanova_appaiata(...)            adonis2 con il soggetto come primo termine
#                                      e le permutazioni entro il soggetto
#   scrivi_appaiati(...)               le tre tabelle
#
# Razionale biologico: quando gli stessi soggetti sono campionati in piu'
# condizioni, la differenza fra soggetti e' di norma maggiore di quella fra
# condizioni. Un test fra gruppi indipendenti la conta come variabilita' entro
# gruppo e non vede l'effetto; un test entro soggetto la toglie, confrontando
# ogni soggetto con se stesso. Vale solo sui soggetti che hanno un campione in
# ogni condizione: gli altri sono esclusi, e lo si dichiara uno per uno.

# Il minimo di soggetti completi: il piu' piccolo n per cui un test dei ranghi
# con segno, con tutte le differenze dello stesso segno, scende sotto il
# livello di significativita' (p esatto minimo a due code 2^(1 - n)). Vale
# anche per il test di Friedman e per la PERMANOVA: una sola regola, prudente.
minimo_di_soggetti <- function(livello) {
  n <- 2L
  while (2^(1 - n) >= livello) n <- n + 1L
  n
}

# I soggetti di un'analisi: chi ha un campione in ogni livello, chi no e
# perche'. `valori` e `soggetti` sono vettori con nome (campione -> valore).
.soggetti_di <- function(valori, soggetti) {
  campioni <- names(valori)[!is.na(valori) & !is.na(soggetti[names(valori)])]
  if (!length(campioni)) {
    # Nessun campione con un valore e un soggetto: non c'e' nulla da appaiare.
    return(list(livelli = character(), soggetti = character(), completi = character(),
                ripetuti = character(), stato = character(),
                campioni = matrix(NA_character_, 0L, 0L), fra_soggetti = FALSE,
                presenti = character(),
                numero = integer()))
  }
  livelli <- unique(as.character(valori[campioni]))
  livelli <- livelli[ordine_radix(livelli)]
  nomi <- unique(as.character(soggetti[campioni]))
  nomi <- nomi[ordine_radix(nomi)]
  conta <- table(factor(soggetti[campioni], levels = nomi),
                 factor(valori[campioni], levels = livelli))
  ripetuti <- nomi[apply(conta > 1L, 1L, any)]
  completi <- nomi[apply(conta == 1L, 1L, all)]
  stato <- vapply(nomi, function(s) {
    if (s %in% ripetuti) {
      return(sprintf("escluso: piu' campioni nello stesso livello (%s)",
                     paste(livelli[conta[s, ] > 1L], collapse = ", ")))
    }
    if (s %in% completi) return("completo")
    sprintf("escluso: manca il livello %s", paste(livelli[conta[s, ] == 0L], collapse = ", "))
  }, character(1))
  # Per ogni soggetto completo il suo campione in ogni livello.
  matrice <- matrix(NA_character_, nrow = length(completi), ncol = length(livelli),
                    dimnames = list(completi, livelli))
  for (campione in campioni) {
    s <- as.character(soggetti[[campione]])
    if (s %in% completi) matrice[s, as.character(valori[[campione]])] <- campione
  }
  list(livelli = livelli, soggetti = nomi, completi = completi, ripetuti = ripetuti,
       stato = stato, campioni = matrice,
       # Vero se nessun soggetto ha campioni in piu' di un livello: la variabile
       # distingue i soggetti, non varia entro soggetto.
       fra_soggetti = all(rowSums(conta > 0L) == 1L),
       presenti = vapply(nomi, function(s) paste(livelli[conta[s, ] > 0L], collapse = ", "),
                         character(1)),
       numero = stats::setNames(as.integer(rowSums(conta)), nomi))
}

piano_appaiato <- function(d, parametri) {
  avvisi <- list()
  avviso <- function(codice, dettaglio) {
    avvisi[[length(avvisi) + 1L]] <<- .avviso(codice, dettaglio)
  }
  piani <- list()
  if (is.null(d$soggetto)) return(list(piani = piani, avvisi = avvisi))
  livello <- parametri$stat$significance_level
  minimo <- minimo_di_soggetti(livello)
  soggetti <- d$soggetto$valori
  senza <- names(soggetti)[is.na(soggetti)]
  if (length(senza)) {
    avviso("E-ECO-13", sprintf(
      "Soggetto '%s': %d campioni con valore mancante (%s); non entrano nelle analisi entro soggetto.",
      d$soggetto$nome, length(senza), .elenco(senza)))
  }
  rarefatti <- setdiff(d$campioni, d$sotto_profondita)
  # I test fra gruppi indipendenti restano, e contano come indipendenti i
  # campioni dello stesso soggetto: chi li legge deve saperlo.
  per_soggetto <- table(soggetti[!is.na(soggetti)])
  if (any(per_soggetto > 1L)) {
    avviso("E-ECO-34", sprintf(paste0(
      "Soggetto '%s': %d soggetti su %d hanno piu' di un campione fra quelli analizzati; ",
      "i confronti fra gruppi indipendenti e le analisi delle variabili continue li ",
      "contano come campioni indipendenti."),
      d$soggetto$nome, sum(per_soggetto > 1L), length(per_soggetto)))
  }

  giudica <- function(nome, analisi, s) {
    motivo <- if (length(s$livelli) < 2L) {
      sprintf("la variabile ha %d livelli fra i campioni con un soggetto (ne servono almeno 2)",
              length(s$livelli))
    } else if (s$fra_soggetti) {
      paste0("ogni soggetto ha campioni in un solo livello: la variabile distingue i ",
             "soggetti e non varia entro soggetto, quindi non c'e' un confronto entro ",
             "soggetto da fare")
    } else if (length(s$ripetuti)) {
      sprintf("%d soggetti hanno piu' di un campione nello stesso livello (%s): non c'e' una sola misura per livello da confrontare",
              length(s$ripetuti), .elenco(s$ripetuti))
    } else if (length(s$completi) < minimo) {
      sprintf("i soggetti con un campione in ogni livello (%s) sono %d su %d, e ne servono almeno %d perche' un test entro soggetto scenda sotto %s",
              paste(s$livelli, collapse = ", "), length(s$completi), length(s$soggetti),
              minimo, format(livello))
    } else {
      NULL
    }
    if (!is.null(motivo)) {
      avviso("E-ECO-32", sprintf("Variabile '%s', soggetto '%s', %s: %s.", nome,
                                 d$soggetto$nome, analisi, motivo))
      return(FALSE)
    }
    esclusi <- setdiff(s$soggetti, s$completi)
    if (length(esclusi)) {
      avviso("E-ECO-33", sprintf(
        "Variabile '%s', soggetto '%s', %s: %d soggetti su %d senza un campione in ogni livello (%s).",
        nome, d$soggetto$nome, analisi, length(esclusi), length(s$soggetti),
        .elenco(esclusi)))
    }
    TRUE
  }

  for (v in d$variabili) {
    alfa <- .soggetti_di(v$valori[rarefatti], soggetti)
    alfa$eseguibile <- giudica(v$nome, "alfa diversita'", alfa)
    multivariata <- .soggetti_di(v$valori, soggetti)
    multivariata$eseguibile <- giudica(v$nome, "PERMANOVA", multivariata)
    piani[[v$nome]] <- list(nome = v$nome, soggetto = d$soggetto$nome, alfa = alfa,
                            permanova = multivariata)
  }
  # Gli strati restringono le permutazioni della PERMANOVA fra gruppi; entro
  # soggetto le permutazioni sono gia' ristrette al soggetto, e gli strati non
  # entrano: chi li ha dichiarati deve saperlo.
  strati <- parametri$stat$permanova_strata
  if (!is.null(strati) &&
      any(vapply(piani, function(p) isTRUE(p$permanova$eseguibile), logical(1)))) {
    avviso("E-ECO-35", sprintf(paste0(
      "Soggetto '%s': stat.permanova_strata ('%s') non si applica alla PERMANOVA entro ",
      "soggetto, le cui permutazioni sono ristrette entro il soggetto."),
      d$soggetto$nome, strati))
  }
  list(piani = piani, avvisi = avvisi)
}

# --- Alfa diversita' -------------------------------------------------------- #

.wilcoxon_con_segno <- function(x, y) {
  # La regola di stats::wilcox.test (exact = NULL) per dati appaiati, applicata
  # in modo esplicito: esatto con meno di 50 coppie, senza differenze nulle e
  # senza valori assoluti ripetuti; altrimenti approssimazione normale con
  # correzione di continuita' (le differenze nulle escono dal test).
  differenze <- x - y
  if (all(differenze == 0)) {
    return(list(test = "Wilcoxon dei ranghi con segno, non calcolabile: nessuna differenza",
                statistica = NA_real_, gl = NA_integer_, p = NA_real_))
  }
  esatto <- length(differenze) < 50L && !any(differenze == 0) &&
    !anyDuplicated(abs(differenze))
  esito <- stats::wilcox.test(x, y, alternative = "two.sided", paired = TRUE,
                              exact = esatto, correct = TRUE, conf.int = FALSE)
  list(
    test = if (esatto) "Wilcoxon dei ranghi con segno, esatto"
           else "Wilcoxon dei ranghi con segno, approssimazione normale con correzione di continuita'",
    statistica = unname(esito$statistic), gl = NA_integer_, p = esito$p.value
  )
}

test_alfa_appaiato <- function(indici, piano) {
  a <- piano$alfa
  if (!a$eseguibile) return(NULL)
  livelli <- a$livelli
  righe <- list()
  riga <- function(indice, confronto, primo, secondo, esito) {
    righe[[length(righe) + 1L]] <<- data.frame(
      variabile = piano$nome, soggetto = piano$soggetto, indice = indice,
      confronto = confronto, livello_1 = primo, livello_2 = secondo,
      soggetti = nrow(a$campioni), test = esito$test, statistica = esito$statistica,
      gl = esito$gl, p = esito$p, stringsAsFactors = FALSE)
  }
  for (indice in INDICI) {
    # Una riga per soggetto completo, una colonna per livello.
    valori <- matrix(as.numeric(indici[[indice]][match(a$campioni, indici$campione)]),
                     nrow = nrow(a$campioni), dimnames = dimnames(a$campioni))
    if (length(livelli) == 2L) {
      riga(indice, paste(livelli, collapse = " vs "), livelli[1L], livelli[2L],
           .wilcoxon_con_segno(valori[, 1L], valori[, 2L]))
    } else {
      # Senza alcuna differenza entro soggetto la statistica non e' definita.
      globale <- if (all(apply(valori, 1L, function(r) length(unique(r)) == 1L))) {
        list(test = "Friedman, non calcolabile: nessuna differenza entro soggetto",
             statistica = NA_real_, gl = NA_integer_, p = NA_real_)
      } else {
        esito <- stats::friedman.test(valori)
        list(test = "Friedman, approssimazione chi quadrato",
             statistica = unname(esito$statistic), gl = as.integer(esito$parameter),
             p = esito$p.value)
      }
      riga(indice, "globale", NA_character_, NA_character_, globale)
      for (i in seq_len(length(livelli) - 1L)) {
        for (j in seq(i + 1L, length(livelli))) {
          riga(indice, paste(livelli[i], livelli[j], sep = " vs "), livelli[i], livelli[j],
               .wilcoxon_con_segno(valori[, i], valori[, j]))
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

# --- PERMANOVA -------------------------------------------------------------- #

# La PERMANOVA entro soggetto e' quella fra gruppi con il soggetto al posto
# delle variabili tecniche e degli strati: stesso calcolo, stesse dispersioni.
permanova_appaiata <- function(matrice, distanza, piano, permutazioni, livello, seme) {
  s <- piano$permanova
  if (!s$eseguibile) return(NULL)
  campioni <- as.vector(s$campioni)
  campioni <- campioni[ordine_radix(campioni)]
  di_chi <- stats::setNames(rep(rownames(s$campioni), times = ncol(s$campioni)),
                            as.vector(s$campioni))[campioni]
  quale <- stats::setNames(rep(colnames(s$campioni), each = nrow(s$campioni)),
                           as.vector(s$campioni))[campioni]
  # Le disposizioni distinte entro soggetto: per ogni soggetto le permutazioni
  # dei suoi livelli.
  logaritmo <- nrow(s$campioni) * lgamma(ncol(s$campioni) + 1)
  finto <- list(nome = piano$nome, permanova = list(
    campioni = campioni, valori = quale, gruppi = s$livelli,
    termini = piano$soggetto, stimabile = TRUE, saturo = FALSE, eseguibile = TRUE,
    disposizioni = if (logaritmo < log(1e15)) round(exp(logaritmo)) else Inf))
  esito <- permanova(matrice, distanza, finto, stats::setNames(list(di_chi), piano$soggetto),
                     di_chi, permutazioni, livello, seme)
  tabella <- esito$tabella
  del_soggetto <- tabella$analisi == "permanova" & tabella$ruolo == "tecnica"
  tabella$ruolo[del_soggetto] <- "soggetto"
  # Le permutazioni entro soggetto non spostano il termine del soggetto.
  tabella$p[del_soggetto] <- NA
  tabella$nota[del_soggetto] <-
    "p non calcolabile: le permutazioni sono entro il soggetto"
  # Con una distanza non euclidea il soggetto e la variabile possono spiegare
  # insieme piu' della variabilita' totale: la somma dei quadrati residua non
  # e' positiva, F e' un rapporto senza significato, e lo si dichiara.
  residua <- tabella$somma_quadrati[tabella$analisi == "permanova" & tabella$ruolo == "residuo"]
  totale <- tabella$somma_quadrati[tabella$analisi == "permanova" & tabella$ruolo == "totale"]
  if (length(residua) && !residuo_positivo(residua, totale)) {
    della_variabile <- tabella$analisi == "permanova" & tabella$ruolo %in% c("analizzata", "soggetto")
    tabella$f[della_variabile] <- NA
    tabella$p[della_variabile] <- NA
    tabella$nota[tabella$analisi == "permanova" & tabella$ruolo == "analizzata"] <-
      nota_residuo_non_positivo(distanza, "il soggetto e la variabile")
  }
  tabella$soggetti <- nrow(s$campioni)
  names(tabella)[names(tabella) == "gruppi"] <- "livelli"
  list(tabella = tabella, dispersioni_diverse = esito$dispersioni_diverse,
       dispersioni_valutabili = esito$dispersioni_valutabili,
       p_dispersioni = esito$p_dispersioni)
}

# I soggetti di ogni variabile e di ogni analisi, uno per riga.
.tabella_dei_soggetti <- function(piani) {
  righe <- list()
  for (p in piani) {
    for (analisi in c("alfa", "permanova")) {
      s <- p[[analisi]]
      if (!length(s$soggetti)) next
      righe[[length(righe) + 1L]] <- data.frame(
        variabile = p$nome, analisi = analisi, soggetto = s$soggetti,
        campioni = as.integer(s$numero[s$soggetti]), livelli_presenti = s$presenti[s$soggetti],
        stato = s$stato[s$soggetti], analisi_eseguita = if (s$eseguibile) "si" else "no",
        stringsAsFactors = FALSE)
    }
  }
  if (length(righe)) do.call(rbind, righe) else data.frame(
    variabile = character(), analisi = character(), soggetto = character(),
    campioni = integer(), livelli_presenti = character(), stato = character(),
    analisi_eseguita = character(), stringsAsFactors = FALSE)
}

scrivi_appaiati <- function(alfa, multivariati, piani, parametri, profondita, cartella) {
  dir.create(cartella, showWarnings = FALSE)
  minimo <- minimo_di_soggetti(parametri$stat$significance_level)
  scrivi_tsv(.tabella_dei_soggetti(piani), file.path(cartella, "soggetti.tsv"), intestazione = c(
    sprintf("Soggetti del disegno a misure ripetute (design.subject: %s), per variabile e per analisi.",
            parametri$design$subject),
    "alfa: fra i campioni rarefatti; permanova: fra i campioni analizzati",
    "completo: un campione in ogni livello della variabile; gli altri sono esclusi, con il motivo",
    sprintf("analisi_eseguita: no se un soggetto ha piu' campioni nello stesso livello, se i livelli sono meno di due o se i soggetti completi sono meno di %d",
            minimo)
  ))
  vuota_alfa <- data.frame(
    variabile = character(), soggetto = character(), indice = character(),
    confronto = character(), livello_1 = character(), livello_2 = character(),
    soggetti = integer(), test = character(), statistica = character(), gl = integer(),
    p = character(), p_bh = character(), dimensione_famiglia = integer(),
    stringsAsFactors = FALSE)
  if (!is.null(alfa)) {
    alfa$statistica <- reale(alfa$statistica)
    for (colonna in c("p", "p_bh")) alfa[[colonna]] <- esponenziale(alfa[[colonna]])
  }
  scrivi_tsv(
    if (is.null(alfa)) vuota_alfa else alfa, file.path(cartella, "alfa_appaiata.tsv"),
    intestazione = c(
      "Confronto dell'alfa diversita' entro soggetto, fra i livelli di ogni variabile.",
      sprintf("indici calcolati sui conteggi rarefatti a %d letture; soli soggetti con un campione rarefatto in ogni livello (test/soggetti.tsv)",
              profondita),
      "due livelli: Wilcoxon dei ranghi con segno; piu' di due: Friedman (confronto 'globale', p dall'approssimazione chi quadrato) e Wilcoxon dei ranghi con segno per ogni coppia di livelli",
      "Wilcoxon esatto con meno di 50 soggetti, senza differenze nulle e senza valori assoluti ripetuti, altrimenti approssimazione normale con correzione di continuita'",
      "p_bh: p corretta con Benjamini-Hochberg sulla famiglia di tutti i test della variabile (ogni indice per ogni confronto)",
      "dimensione_famiglia: numero di test della famiglia; p e p_bh in notazione esponenziale a sei decimali",
      sprintf("il confronto non si esegue con meno di %d soggetti completi: il minimo perche' un test entro soggetto scenda sotto %s",
              minimo, format(parametri$stat$significance_level))
    )
  )
  vuota <- data.frame(
    variabile = character(), distanza = character(), analisi = character(),
    termine = character(), ruolo = character(), gl = integer(),
    somma_quadrati = character(), r2 = character(), f = character(), p = character(),
    permutazioni = integer(), campioni = integer(), livelli = integer(),
    nota = character(), soggetti = integer(), stringsAsFactors = FALSE)
  if (!is.null(multivariati)) {
    for (colonna in c("somma_quadrati", "r2", "f", "p")) {
      multivariati[[colonna]] <- reale(multivariati[[colonna]])
    }
  }
  scrivi_tsv(
    if (is.null(multivariati)) vuota else multivariati,
    file.path(cartella, "permanova_appaiata.tsv"),
    intestazione = c(
      "PERMANOVA e omogeneita' delle dispersioni entro soggetto, per variabile e per distanza.",
      "permanova: vegan::adonis2, by = \"terms\": il soggetto e' il primo termine, poi la variabile; le permutazioni sono ristrette entro il soggetto",
      "soli campioni dei soggetti con un campione in ogni livello (test/soggetti.tsv); le variabili tecniche non entrano nel modello: cio' che e' costante entro soggetto lo assorbe il termine del soggetto",
      "dispersioni: vegan::betadisper sui livelli della variabile, con permutest e le stesse permutazioni entro il soggetto",
      "con una distanza non euclidea (Bray-Curtis, Jaccard, UniFrac) la somma dei quadrati residua puo' non essere positiva: F e p non sono allora valutabili, e la nota lo dice",
      sprintf("permutazioni richieste: %d; la colonna permutazioni riporta quelle usate (tutte le disposizioni, se sono poche)",
              as.integer(parametri$stat$permanova_permutations)),
      if (!is.null(parametri$stat$permanova_strata)) sprintf(
        "stat.permanova_strata (%s) non si applica: le permutazioni sono ristrette entro il soggetto",
        parametri$stat$permanova_strata),
      sprintf("livello di significativita' per la nota sulle dispersioni: %s",
              format(parametri$stat$significance_level)),
      "ruolo: soggetto, analizzata, residuo, totale; gl: gradi di liberta'; r2 e f solo per la PERMANOVA"
    )
  )
  c("test/soggetti.tsv", "test/alfa_appaiata.tsv", "test/permanova_appaiata.tsv")
}
