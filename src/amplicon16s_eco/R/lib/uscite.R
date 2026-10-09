# Scrittura delle uscite: tabelle e grafici identici byte per byte a parita' di
# oggetto, configurazione e seme.
#
# Richiede io_json.R della pipeline (scrivi_atomico). Interfaccia:
#   reale(x)                         numeri reali in formato fisso
#   scrivi_tsv(tabella, percorso, intestazione)
#   ordine_radix(x)                  ordinamento di stringhe indipendente dalla
#                                    lingua della macchina
#   nome_file(indice, nome)          nome di file ricavato da una variabile
#   apri_png(percorso, ...)          dispositivo grafico senza metadati di tempo
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

# Colori dei gruppi: una tavolozza qualitativa fissa per numero di gruppi.
colori_gruppi <- function(n) {
  if (n < 1L) return(character())
  grDevices::hcl.colors(max(n, 2L), palette = "Dark 3")[seq_len(n)]
}
