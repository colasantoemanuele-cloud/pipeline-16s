# Alfa diversita' su conteggi rarefatti.
#
# La usa R/analisi.R (comando eco run). Scrive, nella cartella alfa/ delle
# uscite, alfa_diversita.tsv e conteggi_rarefatti.tsv.
#
# Richiede uscite.R. Interfaccia:
#   fissa_seme(seme)                 generatore e seme, prima di ogni operazione
#                                    casuale
#   rarefai(ps, campioni, profondita, seme)   matrice varianti x campioni
#   indici_alfa(rarefatti)           ricchezza, Shannon, Gini-Simpson, Chao1
#   scrivi_alfa(...)                 le due tabelle
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

# I conteggi rarefatti, senza reinserimento, dei campioni indicati alla
# profondita' indicata. Restituisce la matrice varianti x campioni, con le
# colonne nell'ordine di campioni e tutte le varianti dell'oggetto, anche
# quelle rimaste senza letture.
rarefai <- function(ps, campioni, profondita, seme) {
  sotto <- phyloseq::prune_samples(campioni, ps)
  fissa_seme(seme)
  rarefatto <- suppressMessages(phyloseq::rarefy_even_depth(
    sotto, sample.size = profondita, rngseed = as.integer(seme),
    replace = FALSE, trimOTUs = FALSE, verbose = FALSE
  ))
  m <- conteggi_di(rarefatto)
  m[, campioni, drop = FALSE]
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

# Scrive nella cartella dell'alfa diversita' la tabella degli indici e quella
# dei conteggi rarefatti, ciascuna con un'intestazione che dichiara
# profondita', sua origine e seme. letture sono le letture dei campioni prima
# della rarefazione; senza_singleton dice se i conteggi di partenza non hanno
# singleton, e allora l'intestazione porta l'avviso su Chao1. Restituisce i
# nomi dei due file, relativi alla cartella di uscita.
scrivi_alfa <- function(indici, rarefatti, letture, profondita, origine, seme, cartella,
                        senza_singleton = FALSE) {
  dir.create(cartella, showWarnings = FALSE)
  scrivi_tsv(
    data.frame(
      campione = indici$campione,
      letture = as.integer(letture[indici$campione]),
      ricchezza_osservata = indici$ricchezza_osservata,
      shannon = reale(indici$shannon),
      gini_simpson = reale(indici$gini_simpson),
      chao1 = reale(indici$chao1),
      stringsAsFactors = FALSE
    ),
    file.path(cartella, "alfa_diversita.tsv"),
    intestazione = c(
      "Alfa diversita' per campione, su conteggi rarefatti.",
      sprintf("rarefazione: senza reinserimento, profondita' %d (%s), seme %d",
              profondita, origine, as.integer(seme)),
      "letture: letture del campione prima della rarefazione",
      "ricchezza_osservata: numero di varianti con almeno una lettura nel campione rarefatto",
      "shannon: -sum(p_i * ln(p_i)), logaritmo naturale, p_i proporzione della variante i",
      "gini_simpson: 1 - sum(p_i^2) (indice di Simpson nella forma di Gini-Simpson)",
      paste0("chao1: S_oss + F1 (F1 - 1) / (2 (F2 + 1)), Chao1 corretto per il bias ",
             "(vegan::estimateR, S.chao1, senza il fattore (N - 1) / N; Chao 1984); F1 e F2: varianti con una e con due letture"),
      if (senza_singleton) AVVISO_CHAO1
    )
  )
  scrivi_tsv(
    affianca(list(variante = rownames(rarefatti)), rarefatti, as.integer), file.path(cartella, "conteggi_rarefatti.tsv"),
    intestazione = c(
      "Conteggi rarefatti: una riga per variante, una colonna per campione.",
      sprintf("rarefazione: senza reinserimento, profondita' %d (%s), seme %d",
              profondita, origine, as.integer(seme))
    )
  )
  c("alfa/alfa_diversita.tsv", "alfa/conteggi_rarefatti.tsv")
}
