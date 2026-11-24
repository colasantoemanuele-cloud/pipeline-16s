# Alfa diversita' su conteggi rarefatti, o non rarefatti (alpha.method none).
#
# La usa R/analisi.R (comando eco run). Scrive, nella cartella alfa/ delle
# uscite, alfa_diversita.tsv e, con la rarefazione, conteggi_rarefatti.tsv
# (con piu' rarefazioni anche alfa_per_rarefazione.tsv); senza rarefazione
# profondita.tsv, la correlazione fra letture e indici.
#
# Richiede uscite.R e, per profondita.tsv, continua.R (.spearman). Interfaccia:
#   fissa_seme(seme)                 generatore e seme, prima di ogni operazione
#                                    casuale
#   rarefai(ps, campioni, profondita, seme, ripetizioni)   lista di matrici
#                                    varianti x campioni, una per rarefazione
#   indici_alfa(rarefatti)           ricchezza, Shannon, Gini-Simpson, Chao1
#   media_degli_indici(per_rarefazione)   la media di ogni indice fra le
#                                    rarefazioni
#   testo_alfa(profondita)           i conteggi su cui sono calcolati gli
#                                    indici, per intestazioni e grafici
#   scrivi_alfa(...)                 le tabelle
#   profondita_e_indici(...)         la correlazione fra letture e indici, con
#                                    i suoi avvisi
#
# Razionale biologico: la rarefazione e' senza reinserimento perche' e' un
# sottocampione delle letture osservate. Con il reinserimento (il predefinito
# di phyloseq::rarefy_even_depth) una variante puo' ricevere piu' letture di
# quante ne abbia, e all'estrazione si aggiunge una variabilita' che il dato
# non ha. Nessuna variante viene tolta dopo la
# rarefazione (trimOTUs = FALSE): la tabella rarefatta ha le stesse righe di
# quella di partenza.
#
# Razionale biologico di Chao1: lo stimatore (Chao 1984, nella forma corretta
# per il bias di vegan::estimateR) aggiunge alla ricchezza osservata una stima
# delle varianti non osservate ricavata dai singleton e dai doubleton. dada2 non
# inferisce varianti con una sola lettura: sui conteggi della pipeline i
# singleton mancano e Chao1 coincide con la ricchezza osservata. Dopo una
# rarefazione i singleton nascono dal sottocampionamento, non dal campione. Per
# questo i singleton dei conteggi di partenza si contano prima del calcolo
# (campioni.R) e, se non ce ne sono, l'avviso lo dichiara anche
# nell'intestazione della tabella.

# Il testo dell'avviso E-ECO-36, ripetuto nell'intestazione della tabella.
AVVISO_CHAO1 <- paste0(
  "avviso E-ECO-36: Chao1 non stima la ricchezza non osservata: i conteggi di partenza non ",
  "hanno singleton (dada2 non li inferisce); dopo una rarefazione i singleton nascono dal ",
  "sottocampionamento, non dal campione")

# Dichiara il generatore di numeri casuali e ne fissa il seme. Va chiamata
# prima di ogni operazione casuale: ognuna riparte cosi' dallo stesso stato,
# qualunque cosa sia stata calcolata prima.
fissa_seme <- function(seme) {
  # Il generatore e' dichiarato, non ereditato: il predefinito di R e' cambiato
  # fra versioni (sample.kind "Rounding" prima della 3.6).
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  set.seed(as.integer(seme))
}

# Le rarefazioni ripetute (alpha.rarefy_iterations) e i loro conteggi.
# Restituisce una lista di matrici varianti x campioni, una per rarefazione,
# con le colonne nell'ordine di campioni e tutte le varianti dell'oggetto,
# anche quelle rimaste senza letture. Il seme si fissa una volta, prima della
# prima rarefazione: le successive continuano la stessa sequenza di numeri
# casuali, e la prima e' quella di una rarefazione singola. phyloseq non rifissa
# il seme (rngseed = FALSE): lo stato e' quello dichiarato da fissa_seme, lo
# stesso che phyloseq imposterebbe con rngseed = seme.
rarefai <- function(ps, campioni, profondita, seme, ripetizioni = 1L) {
  sotto <- phyloseq::prune_samples(campioni, ps)
  fissa_seme(seme)
  lapply(seq_len(ripetizioni), function(i) {
    rarefatto <- suppressMessages(phyloseq::rarefy_even_depth(
      sotto, sample.size = profondita, rngseed = FALSE,
      replace = FALSE, trimOTUs = FALSE, verbose = FALSE
    ))
    m <- conteggi_di(rarefatto)
    m[, campioni, drop = FALSE]
  })
}

# La media di ogni indice fra le rarefazioni: una tabella come quella di
# indici_alfa, con la ricchezza come numero reale se le rarefazioni sono piu'
# di una.
media_degli_indici <- function(per_rarefazione) {
  if (length(per_rarefazione) == 1L) return(per_rarefazione[[1L]])
  media <- per_rarefazione[[1L]]
  for (indice in setdiff(names(media), "campione")) {
    media[[indice]] <- rowMeans(do.call(cbind, lapply(per_rarefazione, `[[`, indice)))
  }
  media
}

# I conteggi su cui sono calcolati gli indici, in parole, per le intestazioni
# delle tabelle e le didascalie dei grafici: rarefatti a una profondita' (con
# il numero di rarefazioni, se piu' d'una) o non rarefatti (profondita NA).
# RIPETIZIONI_ALFA lo imposta R/analisi.R.
RIPETIZIONI_ALFA <- 1L
testo_alfa <- function(profondita) {
  if (is.na(profondita)) return("conteggi non rarefatti (alpha.method none)")
  if (RIPETIZIONI_ALFA > 1L) {
    return(sprintf("conteggi rarefatti a %d letture (media di %d rarefazioni)",
                   profondita, RIPETIZIONI_ALFA))
  }
  sprintf("conteggi rarefatti a %d letture", profondita)
}

# Gli indici di alfa diversita' di ogni campione dalla matrice rarefatta
# (varianti x campioni): ricchezza osservata, Shannon con il logaritmo naturale,
# Gini-Simpson e Chao1 corretto per il bias (vegan::estimateR, S.chao1).
# Restituisce una tabella con una riga per campione.
indici_alfa <- function(rarefatti) {
  per_campione <- t(rarefatti)
  # estimateR vuole interi: i conteggi rarefatti lo sono gia'.
  storage.mode(per_campione) <- "integer"
  stimati <- vegan::estimateR(per_campione)
  data.frame(
    campione = colnames(rarefatti),
    ricchezza_osservata = as.integer(colSums(rarefatti > 0)),
    shannon = as.numeric(vegan::diversity(per_campione, index = "shannon",
                                         MARGIN = 1, base = exp(1))),
    gini_simpson = as.numeric(vegan::diversity(per_campione, index = "simpson",
                                              MARGIN = 1)),
    chao1 = as.numeric(stimati["S.chao1", ]),
    stringsAsFactors = FALSE
  )
}

# Scrive nella cartella dell'alfa diversita' la tabella degli indici e, con la
# rarefazione, quella dei conteggi rarefatti (la prima rarefazione) e, con piu'
# rarefazioni, gli indici di ognuna; ogni tabella ha un'intestazione che
# dichiara profondita', sua origine, rarefazioni e seme. rarefatti e' la lista
# delle matrici delle rarefazioni (nulla senza rarefazione), per_rarefazione
# gli indici di ognuna; letture sono le letture dei campioni prima della
# rarefazione; senza_singleton dice se i conteggi di partenza non hanno
# singleton, e allora l'intestazione porta l'avviso su Chao1. Restituisce i
# nomi dei file, relativi alla cartella di uscita.
scrivi_alfa <- function(indici, rarefatti, letture, profondita, origine, seme, cartella,
                        senza_singleton = FALSE, per_rarefazione = NULL) {
  dir.create(cartella, showWarnings = FALSE)
  rarefatta <- !is.na(profondita)
  ripetizioni <- if (rarefatta) length(rarefatti) else 0L
  # Con piu' rarefazioni la ricchezza e' una media, non piu' un intero.
  ricchezza <- if (ripetizioni > 1L) reale(indici$ricchezza_osservata) else
    as.integer(indici$ricchezza_osservata)
  descrizione_rarefazione <- if (!rarefatta) {
    "nessuna rarefazione (alpha.method none): indici sui conteggi non rarefatti, dopo l'eventuale filtro"
  } else if (ripetizioni > 1L) {
    sprintf("rarefazione: senza reinserimento, profondita' %d (%s), %d rarefazioni dal seme %d; ogni indice e' la media sulle rarefazioni",
            profondita, origine, ripetizioni, as.integer(seme))
  } else {
    sprintf("rarefazione: senza reinserimento, profondita' %d (%s), seme %d",
            profondita, origine, as.integer(seme))
  }
  scrivi_tsv(
    data.frame(
      campione = indici$campione,
      letture = as.integer(letture[indici$campione]),
      ricchezza_osservata = ricchezza,
      shannon = reale(indici$shannon),
      gini_simpson = reale(indici$gini_simpson),
      chao1 = reale(indici$chao1),
      stringsAsFactors = FALSE
    ),
    file.path(cartella, "alfa_diversita.tsv"),
    intestazione = c(
      if (rarefatta) "Alfa diversita' per campione, su conteggi rarefatti." else
        "Alfa diversita' per campione, su conteggi non rarefatti.",
      descrizione_rarefazione,
      if (rarefatta) "letture: letture del campione prima della rarefazione" else
        "letture: letture del campione su cui sono calcolati gli indici",
      if (rarefatta) "ricchezza_osservata: numero di varianti con almeno una lettura nel campione rarefatto" else
        "ricchezza_osservata: numero di varianti con almeno una lettura nel campione",
      "shannon: -sum(p_i * ln(p_i)), logaritmo naturale, p_i proporzione della variante i",
      "gini_simpson: 1 - sum(p_i^2) (indice di Simpson nella forma di Gini-Simpson)",
      paste0("chao1: S_oss + F1 (F1 - 1) / (2 (F2 + 1)), Chao1 corretto per il bias ",
             "(vegan::estimateR, S.chao1, senza il fattore (N - 1) / N; Chao 1984); F1 e F2: varianti con una e con due letture"),
      if (senza_singleton) AVVISO_CHAO1
    )
  )
  if (!rarefatta) return("alfa/alfa_diversita.tsv")
  prima <- rarefatti[[1L]]
  scrivi_tsv(
    affianca(list(variante = rownames(prima)), prima, as.integer), file.path(cartella, "conteggi_rarefatti.tsv"),
    intestazione = c(
      "Conteggi rarefatti: una riga per variante, una colonna per campione.",
      descrizione_rarefazione,
      if (ripetizioni > 1L) "i conteggi sono quelli della prima rarefazione; gli indici di ognuna sono in alfa_per_rarefazione.tsv"
    )
  )
  prodotti <- c("alfa/alfa_diversita.tsv", "alfa/conteggi_rarefatti.tsv")
  if (ripetizioni > 1L) {
    lunga <- do.call(rbind, lapply(seq_along(per_rarefazione), function(i) {
      x <- per_rarefazione[[i]]
      data.frame(campione = x$campione, rarefazione = i,
                 ricchezza_osservata = as.integer(x$ricchezza_osservata),
                 shannon = reale(x$shannon), gini_simpson = reale(x$gini_simpson),
                 chao1 = reale(x$chao1), stringsAsFactors = FALSE)
    }))
    lunga <- lunga[order(lunga$campione, lunga$rarefazione, method = "radix"), , drop = FALSE]
    scrivi_tsv(lunga, file.path(cartella, "alfa_per_rarefazione.tsv"), intestazione = c(
      "Indici di alfa diversita' di ogni rarefazione, per campione.",
      descrizione_rarefazione,
      "rarefazione: il numero della rarefazione, nell'ordine in cui sono estratte dal seme",
      "alfa_diversita.tsv riporta, per ogni indice, la media di queste righe"
    ))
    prodotti <- c(prodotti, "alfa/alfa_per_rarefazione.tsv")
  }
  prodotti
}

# Senza rarefazione: la correlazione di Spearman fra le letture di ogni
# campione e ciascun indice, sui campioni dell'alfa diversita'. Con
# p < livello l'indice dipende dalla profondita' di sequenziamento, e i
# confronti fra gruppi possono rifletterla: avviso E-ECO-39. Scrive
# profondita.tsv; restituisce list(file, avvisi).
profondita_e_indici <- function(indici, letture, livello, cartella) {
  x <- as.numeric(letture[indici$campione])
  avvisi <- list()
  righe <- lapply(INDICI, function(indice) {
    esito <- .spearman(x, as.numeric(indici[[indice]]))
    if (!is.na(esito$p) && esito$p < livello) {
      avvisi[[length(avvisi) + 1L]] <<- list(codice = "E-ECO-39", dettaglio = sprintf(
        "Indice %s: rho di Spearman con le letture %s, p %s, sotto %s (%d campioni).",
        indice, reale(esito$rho), esponenziale(esito$p), format(livello), length(x)))
    }
    data.frame(indice = indice, campioni = length(x), test = esito$test,
               rho = reale(esito$rho), p = esponenziale(esito$p),
               dipende_dalla_profondita = if (!is.na(esito$p) && esito$p < livello) "si" else "no",
               stringsAsFactors = FALSE)
  })
  scrivi_tsv(do.call(rbind, righe), file.path(cartella, "profondita.tsv"), intestazione = c(
    "Correlazione fra le letture di ogni campione e ciascun indice di alfa diversita' (alpha.method none).",
    "letture: quelle su cui sono calcolati gli indici, dopo l'eventuale filtro di abbondanza",
    "correlazione per ranghi di Spearman, a due code, come in test/alfa_continua.tsv; rho: coefficiente, p in notazione esponenziale",
    sprintf("dipende_dalla_profondita: si se p < %s (stat.significance_level), con l'avviso E-ECO-39: i confronti fra gruppi dell'indice possono riflettere la profondita' di sequenziamento",
            format(livello)),
    "p non corretta per confronti multipli: e' un controllo, da leggere indice per indice"
  ))
  list(file = "alfa/profondita.tsv", avvisi = avvisi)
}
