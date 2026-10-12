# Ricerca e taglio del primer nelle prime basi delle letture, per i protocolli
# che lo fanno precedere da un numero variabile di basi.
#
# Richiede Biostrings e ShortRead. Interfaccia:
#   trova_primer(sequenze, primer, discordanze_massime, sfasamento_massimo)
#       -> list(sfasamento, discordanze): un intero per lettura, NA dove il
#          primer non si trova
#   taglia_primer(letture, ricerca) -> list(letture, esito)
#       le letture (ShortReadQ) tagliate dopo il primer, o tolte, secondo
#       ricerca$senza_primer, e per ogni lettura in ingresso sfasamento e
#       discordanze
#
# La regola, la stessa di RicercaPrimer in io_layer/reads.py:
# - codici IUPAC nel primer, basi della lettura alla lettera: una N della
#   lettura, o un altro codice di ambiguita', e' una discordanza anche contro
#   una N del primer (nella finestra si sostituiscono con il segno "-", che non
#   concorda con alcun codice);
# - sole sostituzioni, fino a discordanze_massime; nessuna inserzione o
#   delezione, cosi' che la posizione del taglio non sia ambigua;
# - il primer comincia a uno sfasamento fra 0 e sfasamento_massimo e sta per
#   intero nella lettura;
# - fra piu' corrispondenze vale quella con meno discordanze e, a parita', la
#   piu' a sinistra.

trova_primer <- function(sequenze, primer, discordanze_massime, sfasamento_massimo) {
  n <- length(sequenze)
  lunghezza <- nchar(primer)
  sfasamento <- rep(NA_integer_, n)
  discordanze <- rep(NA_integer_, n)
  if (n == 0L) return(list(sfasamento = sfasamento, discordanze = discordanze))

  larghezze <- Biostrings::width(sequenze)
  # Solo le basi in cui il primer puo' stare: una corrispondenza trovata qui
  # comincia entro lo sfasamento massimo.
  ampiezza <- pmin(larghezze, as.integer(sfasamento_massimo) + lunghezza)
  # Ogni lettera della lettura che non e' una base (N e gli altri codici di
  # ambiguita') diventa il segno "-", che non concorda con alcun codice.
  finestra <- Biostrings::chartr(
    "NRYSWKMBDHV", "-----------", Biostrings::subseq(sequenze, 1L, ampiezza))
  modello <- Biostrings::DNAString(primer)

  # Dal numero di discordanze piu' basso: a ogni livello le letture ancora
  # senza corrispondenza prendono la piu' a sinistra fra quelle del livello.
  for (livello in 0:as.integer(discordanze_massime)) {
    restano <- which(is.na(sfasamento))
    if (!length(restano)) break
    inizi <- Biostrings::startIndex(Biostrings::vmatchPattern(
      modello, finestra[restano], max.mismatch = livello,
      fixed = c(pattern = FALSE, subject = TRUE)
    ))
    quanti <- lengths(inizi)
    if (!sum(quanti)) next
    tutti <- unlist(inizi, use.names = FALSE)
    di_chi <- rep(restano, quanti)
    # Una corrispondenza che esce dai bordi della lettura non conta.
    dentro <- tutti >= 1L & tutti + lunghezza - 1L <= ampiezza[di_chi]
    tutti <- tutti[dentro]
    di_chi <- di_chi[dentro]
    # Gli inizi di ogni lettura sono crescenti: assegnandoli a ritroso resta
    # il piu' piccolo.
    sfasamento[rev(di_chi)] <- rev(tutti) - 1L
    discordanze[unique(di_chi)] <- livello
  }
  list(sfasamento = sfasamento, discordanze = discordanze)
}

# Dove comincia cio' che resta di ogni lettura: dopo l'ultima base del primer,
# o NA se il primer non c'e'.
inizio_dopo_il_primer <- function(letture, ricerca) {
  esito <- trova_primer(
    ShortRead::sread(letture), ricerca$primer,
    ricerca$discordanze_massime, ricerca$sfasamento_massimo
  )
  esito$inizio <- esito$sfasamento + nchar(ricerca$primer) + 1L
  esito
}

# Le letture di un blocco dopo il taglio. Con scarta vero restano solo quelle
# con il primer; altrimenti quelle senza restano intere. Una lettura che dopo
# il primer non ha basi non resta in nessun caso: non c'e' nulla da filtrare.
# tieni, se dato, restringe ancora le letture da tenere (in modalita' appaiata:
# la compagna deve restare anch'essa).
dopo_il_taglio <- function(letture, esito, scarta) {
  larghezze <- Biostrings::width(ShortRead::sread(letture))
  trovato <- !is.na(esito$inizio)
  inizio <- ifelse(trovato, esito$inizio, 1L)
  resta <- inizio <= larghezze & (trovato | !scarta)
  list(inizio = as.integer(inizio), resta = resta, trovato = trovato)
}

# Conteggio di sfasamenti e discordanze di un blocco, come tabella
# sfasamento x discordanze sommabile fra blocchi (nomi "s:d").
conta_posizioni <- function(esito) {
  trovato <- !is.na(esito$sfasamento)
  table(paste(esito$sfasamento[trovato], esito$discordanze[trovato], sep = ":"))
}

somma_posizioni <- function(a, b) {
  nomi <- union(names(a), names(b))
  somma <- stats::setNames(numeric(length(nomi)), nomi)
  somma[names(a)] <- somma[names(a)] + as.numeric(a)
  somma[names(b)] <- somma[names(b)] + as.numeric(b)
  somma
}

# Taglia il primer dalle letture di un file (o dei due file di una coppia) e
# scrive cio' che resta. ingressi, uscite e ricerche hanno un elemento per
# lettura (forward, e reverse in modalita' appaiata), nello stesso ordine. I
# due file di una coppia si leggono in parallelo, blocco per blocco: una coppia
# resta solo se restano tutte e due le sue letture, cosi' i file in uscita hanno
# le stesse coppie nello stesso ordine.
#
# Restituisce list(letture, con_primer (uno per lettura), senza_primer (letture
# o coppie in cui manca almeno un primer), rimaste, posizioni (una tabella per
# lettura)).
taglia_file <- function(ingressi, uscite, ricerche, scarta, blocco = 100000L) {
  lettori <- lapply(ingressi, ShortRead::FastqStreamer, n = blocco)
  on.exit(for (lettore in lettori) close(lettore))
  unlink(uscite)
  letture <- 0
  rimaste <- 0
  senza <- 0
  con_primer <- numeric(length(ingressi))
  posizioni <- rep(list(numeric()), length(ingressi))
  repeat {
    parti <- lapply(lettori, ShortRead::yield)
    quante <- vapply(parti, length, integer(1))
    if (length(unique(quante)) != 1L) {
      stop("i due file di una coppia non hanno lo stesso numero di letture: ",
           paste(basename(ingressi), collapse = ", "))
    }
    if (quante[1L] == 0L) break
    letture <- letture + quante[1L]
    esiti <- Map(inizio_dopo_il_primer, parti, ricerche)
    tagli <- Map(dopo_il_taglio, parti, esiti, MoreArgs = list(scarta = scarta))
    resta <- Reduce(`&`, lapply(tagli, `[[`, "resta"))
    senza <- senza + sum(!Reduce(`&`, lapply(tagli, `[[`, "trovato")))
    rimaste <- rimaste + sum(resta)
    for (i in seq_along(parti)) {
      con_primer[i] <- con_primer[i] + sum(tagli[[i]]$trovato)
      posizioni[[i]] <- somma_posizioni(posizioni[[i]], conta_posizioni(esiti[[i]]))
      if (any(resta)) {
        ShortRead::writeFastq(
          ShortRead::narrow(parti[[i]][resta], start = tagli[[i]]$inizio[resta]),
          uscite[i], mode = "a", compress = TRUE
        )
      }
    }
  }
  list(letture = letture, con_primer = con_primer, senza_primer = senza,
       rimaste = rimaste, posizioni = posizioni)
}
