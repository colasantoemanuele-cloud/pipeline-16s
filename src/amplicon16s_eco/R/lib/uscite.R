# Scrittura delle uscite: tabelle e grafici identici byte per byte a parita' di
# oggetto, configurazione e seme.
#
# Richiede io_json.R della pipeline (scrivi_atomico). Interfaccia:
#   reale(x)                         numeri reali in formato fisso
#   esponenziale(x)                  p dei test analitici, in formato fisso
#   scrivi_tsv(tabella, percorso, intestazione)
#   affianca(fisse, matrice, formato)  colonne fisse e una colonna per ogni
#                                    colonna della matrice, qualunque nome abbia
#   ordine_radix(x)                  ordinamento di stringhe indipendente dalla
#                                    lingua della macchina
#   nome_file(indice, nome)          nome di file ricavato da una variabile
#   apri_png(percorso, ...)          dispositivo grafico senza metadati di tempo
#   nome_distanza(metodo)            il nome leggibile di una distanza
#   larghezza_per_gruppi(n)          larghezza di un grafico con n gruppi in ascissa
#
# Razionale sistemistico: le tabelle non passano per write.table, il cui
# formato dei numeri dipende dalle opzioni della sessione (scipen, digits,
# OutDec). Ogni numero reale e' scritto con dieci decimali in formato fisso;
# gli interi come interi. Un valore non calcolabile e' "NA".

DECIMALI <- 10L

reale <- function(x) {
  x <- as.numeric(x)
  # Lo zero negativo e i valori che si arrotondano a zero darebbero "-0.0000000000".
  x[!is.na(x) & abs(x) < 0.5 * 10^(-DECIMALI)] <- 0
  testo <- sprintf(paste0("%.", DECIMALI, "f"), x)
  testo[is.na(x)] <- "NA"
  testo
}

# Per le p dei test analitici, che possono essere molto piccole: notazione
# esponenziale con sei decimali, fissa anch'essa.
esponenziale <- function(x) {
  x <- as.numeric(x)
  testo <- sprintf("%.6e", x)
  testo[is.na(x)] <- "NA"
  testo
}

ordine_radix <- function(x) order(as.character(x), method = "radix")

.cella <- function(x) {
  # Una tabulazione o un a capo dentro un valore romperebbero la tabella.
  x <- gsub("[\t\r\n]", " ", as.character(x))
  x[is.na(x)] <- "NA"
  x
}

scrivi_tsv <- function(tabella, percorso, intestazione = character()) {
  # `tabella` e' un data.frame di colonne gia' in forma di testo o di interi;
  # `intestazione` sono righe di commento, scritte con il prefisso "# " prima
  # della riga dei nomi delle colonne.
  colonne <- lapply(tabella, .cella)
  righe <- if (nrow(tabella) > 0L) do.call(paste, c(colonne, sep = "\t")) else character()
  scrivi_atomico(
    c(if (length(intestazione)) paste0("# ", intestazione),
      paste(.cella(names(tabella)), collapse = "\t"),
      righe),
    percorso
  )
}

affianca <- function(fisse, matrice, formato = as.character) {
  # Le colonne della matrice (campioni o gruppi) si accostano a quelle fisse
  # senza passare per l'assegnazione per nome: un campione o un gruppo che si
  # chiamasse come una colonna fissa la sovrascriverebbe in silenzio.
  colonne <- lapply(seq_len(ncol(matrice)), function(j) formato(matrice[, j]))
  tabella <- c(as.list(fisse), colonne)
  names(tabella) <- c(names(fisse), colnames(matrice))
  structure(tabella, class = "data.frame", row.names = .set_row_names(nrow(matrice)))
}

nome_file <- function(indice, nome) {
  # Il nome di una colonna puo' contenere spazi o segni non ammessi in un nome
  # di file: si tengono lettere, cifre e trattino basso. L'indice, cioe' la
  # posizione della variabile nella configurazione, rende il nome unico anche
  # quando due colonne diventano uguali dopo la sostituzione.
  sprintf("%02d_%s", indice, gsub("[^A-Za-z0-9_]", "_", nome))
}

apri_png <- function(percorso, larghezza = 2400L, altezza = 1500L) {
  # PNG del dispositivo cairo di R: il file non contiene data ne' ora, e due
  # esecuzioni danno gli stessi byte. I formati vettoriali di R (pdf, svg)
  # registrano la data di creazione o dipendono dalla versione della libreria.
  grDevices::png(percorso, width = larghezza, height = altezza, res = 200,
                 type = "cairo", bg = "white")
}

# I nomi con cui le distanze compaiono nei titoli dei grafici.
NOMI_DISTANZA <- c(
  bray = "Bray-Curtis", jaccard = "Jaccard", aitchison = "Aitchison",
  unifrac_weighted = "UniFrac pesata", unifrac_unweighted = "UniFrac non pesata"
)

nome_distanza <- function(metodo) {
  if (metodo %in% names(NOMI_DISTANZA)) NOMI_DISTANZA[[metodo]] else metodo
}

# Legenda dei gruppi. Fino a SOGLIA_LEGENDA_COLONNE gruppi sta a destra del
# grafico in una colonna; oltre, in due colonne, e l'immagine si allarga di
# quanto serve; oltre SOGLIA_LEGENDA_PANNELLO va in un pannello sotto il
# grafico, e l'immagine si allunga con il numero di righe. In tutti i casi
# l'area del grafico resta quella di un grafico con pochi gruppi.
SOGLIA_LEGENDA_COLONNE <- 30L
SOGLIA_LEGENDA_PANNELLO <- 60L
LARGHEZZA_GRAFICO <- 2400L
ALTEZZA_GRAFICO <- 1500L
# Pixel di una colonna di legenda a destra, di una riga di legenda nel
# pannello, e di un gruppo in piu' sull'ascissa (a 200 punti per pollice).
PIXEL_COLONNA_LEGENDA <- 620L
PIXEL_RIGA_LEGENDA <- 34L
PIXEL_PER_GRUPPO <- 45L

# Con molti gruppi in ascissa (scatole, barre) l'immagine si allarga, invece di
# stringere ogni gruppo.
larghezza_per_gruppi <- function(n) {
  LARGHEZZA_GRAFICO + max(0L, n - SOGLIA_LEGENDA_COLONNE) * PIXEL_PER_GRUPPO
}

# Colori dei gruppi: una tavolozza qualitativa fissa per numero di gruppi.
colori_gruppi <- function(n) {
  if (n < 1L) return(character())
  grDevices::hcl.colors(max(n, 2L), palette = "Dark 3")[seq_len(n)]
}
