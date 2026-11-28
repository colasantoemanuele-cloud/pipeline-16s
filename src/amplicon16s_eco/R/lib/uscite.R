# Scrittura delle uscite: tabelle e grafici identici byte per byte a parita' di
# oggetto, configurazione e seme.
#
# La caricano i tre script del pacchetto (R/analisi.R, R/valida.R,
# R/descrivi.R) e, attraverso di loro, ogni altra libreria: nessun file delle
# uscite si scrive in altro modo.
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
#   imposta_palette(palette), colori_gruppi(n), origine_palette(palette)
#                                    la palette dei gruppi (plot.palette)
#   colori_della_palette(), nota_palette(n)
#                                    i colori della palette prima
#                                    dell'interpolazione e la nota E-ECO-46
#
# Razionale sistemistico: le tabelle non passano per write.table, il cui
# formato dei numeri dipende dalle opzioni della sessione (scipen, digits,
# OutDec). Ogni numero reale e' scritto con dieci decimali in formato fisso;
# gli interi come interi. Un valore non calcolabile e' "NA".

DECIMALI <- 10L

# Un numero reale come testo in formato fisso, con DECIMALI decimali; NA per un
# valore mancante.
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

# L'ordine dei valori come stringhe, byte per byte: non dipende dalla lingua
# della macchina.
ordine_radix <- function(x) order(as.character(x), method = "radix")

# Un valore come cella di una tabella: testo, senza tabulazioni ne' a capo, NA
# per un valore mancante.
.cella <- function(x) {
  # Una tabulazione o un a capo dentro un valore romperebbero la tabella.
  x <- gsub("[\t\r\n]", " ", as.character(x))
  x[is.na(x)] <- "NA"
  x
}

# Scrive una tabella separata da tabulazioni in modo atomico, con le righe di
# commento dell'intestazione prima dei nomi delle colonne.
# Le righe che ogni tabella aggiunge alla sua intestazione quando i dati su cui
# si calcola non sono le varianti dell'oggetto: il livello dei calcoli
# (run.analysis_rank) e il filtro di abbondanza. Le imposta R/analisi.R; vuote,
# le intestazioni non hanno quelle righe.
NOTE_DEI_DATI <- character()

scrivi_tsv <- function(tabella, percorso, intestazione = character()) {
  # `tabella` e' un data.frame di colonne gia' in forma di testo o di interi;
  # `intestazione` sono righe di commento, scritte con il prefisso "# " prima
  # della riga dei nomi delle colonne, seguite da NOTE_DEI_DATI.
  intestazione <- c(intestazione, NOTE_DEI_DATI)
  colonne <- lapply(tabella, .cella)
  righe <- if (nrow(tabella) > 0L) do.call(paste, c(colonne, sep = "\t")) else character()
  scrivi_atomico(
    c(if (length(intestazione)) paste0("# ", intestazione),
      paste(.cella(names(tabella)), collapse = "\t"),
      righe),
    percorso
  )
}

# Una tabella con le colonne fisse seguite da una colonna per ogni colonna
# della matrice, formattata con formato; i nomi delle colonne della matrice
# restano quelli, qualunque siano.
affianca <- function(fisse, matrice, formato = as.character) {
  # Le colonne della matrice (campioni o gruppi) si accostano a quelle fisse
  # senza passare per l'assegnazione per nome: un campione o un gruppo che si
  # chiamasse come una colonna fissa la sovrascriverebbe in silenzio.
  colonne <- lapply(seq_len(ncol(matrice)), function(j) formato(matrice[, j]))
  tabella <- c(as.list(fisse), colonne)
  names(tabella) <- c(names(fisse), colnames(matrice))
  structure(tabella, class = "data.frame", row.names = .set_row_names(nrow(matrice)))
}

# La parte di un nome di file ricavata da una variabile: la sua posizione nella
# configurazione, a due cifre, e il nome con i soli caratteri ammessi.
nome_file <- function(indice, nome) {
  # Il nome di una colonna puo' contenere spazi o segni non ammessi in un nome
  # di file: si tengono lettere, cifre e trattino basso. L'indice, cioe' la
  # posizione della variabile nella configurazione, rende il nome unico anche
  # quando due colonne diventano uguali dopo la sostituzione.
  sprintf("%02d_%s", indice, gsub("[^A-Za-z0-9_]", "_", nome))
}

# Apre un dispositivo grafico PNG di dimensioni in pixel date, a 200 punti per
# pollice, su fondo bianco; chi chiama lo chiude.
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

# Il nome leggibile di una distanza per i titoli; il nome del metodo, se non ne
# ha uno.
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

# --- Palette dei gruppi (plot.palette) ------------------------------------ #
#
# Il principio e' quello di get.pal del pacchetto AMBER di Aldo Dale
# (https://github.com/AldoDale/AMBER, R/themes_and_palettes.R), riusato con il
# permesso dell'autore: una palette con nome, estesa per interpolazione quando i
# gruppi sono piu' dei suoi colori. AMBER interpola con ggpubr::get_palette, che
# non e' nell'immagine; qui l'interpolazione e' grDevices::colorRampPalette.
#
# Come si legge il nome, in quest'ordine: "Okabe-Ito" (grDevices::palette.colors,
# adatta ai daltonici, senza il grigio, riservato alla voce «altri»); un nome di RColorBrewer (i colori originali di
# ColorBrewer: per i 26 nomi che hcl.pals() ha anche come approssimazione HCL
# vale RColorBrewer); un nome di grDevices::hcl.pals() (hcl.colors genera tanti
# colori quanti gruppi, come fa il predefinito "Dark 3"); un
# nome di viridisLite. Un elenco di colori esadecimali vale cosi' com'e'.

# La palette in uso: la imposta R/analisi.R con imposta_palette; il predefinito
# e' "Dark 3" (plot.palette).
PALETTE_GRAFICI <- "Dark 3"
OKABE_ITO <- "Okabe-Ito"
OPZIONI_VIRIDIS <- c("magma", "inferno", "plasma", "viridis", "cividis", "rocket",
                     "mako", "turbo")

# Da dove viene una palette: "esadecimali", "okabe_ito", "brewer", "hcl",
# "viridis", oppure NA se il nome non e' riconosciuto.
origine_palette <- function(palette) {
  palette <- as.character(unlist(palette))
  if (length(palette) != 1L || startsWith(palette, "#")) return("esadecimali")
  if (identical(palette, OKABE_ITO)) return("okabe_ito")
  if (requireNamespace("RColorBrewer", quietly = TRUE) &&
      palette %in% rownames(RColorBrewer::brewer.pal.info)) return("brewer")
  if (palette %in% grDevices::hcl.pals()) return("hcl")
  if (palette %in% OPZIONI_VIRIDIS && requireNamespace("viridisLite", quietly = TRUE)) {
    return("viridis")
  }
  NA_character_
}

# I colori di Okabe-Ito senza il grigio, l'ultimo: nei grafici della
# composizione il grigio e' la voce «altri», e un taxon dello stesso colore (o
# di uno interpolato verso il grigio) vi si confonderebbe.
okabe_ito_senza_grigio <- function() {
  colori <- unname(grDevices::palette.colors(NULL, palette = OKABE_ITO))
  colori[toupper(colori) != "#999999"]
}

# L'avviso E-ECO-46, come compare nelle intestazioni delle tabelle.
AVVISO_PALETTE <- "palette estesa per interpolazione: i colori aggiunti non sono garantiti per chi non distingue i colori"

# Quanti colori ha la palette prima di essere estesa per interpolazione; NA per
# le palette che generano tanti colori quanti ne servono (hcl, viridis), che
# non si estendono.
colori_della_palette <- function(palette = PALETTE_GRAFICI) {
  palette <- as.character(unlist(palette))
  origine <- origine_palette(palette)
  if (is.na(origine)) return(NA_integer_)
  switch(
    origine,
    okabe_ito = length(okabe_ito_senza_grigio()),
    brewer = as.integer(RColorBrewer::brewer.pal.info[palette, "maxcolors"]),
    esadecimali = length(palette),
    NA_integer_
  )
}

# Il nome della palette per i messaggi.
nome_palette <- function(palette = PALETTE_GRAFICI) {
  palette <- as.character(unlist(palette))
  if (origine_palette(palette) %in% "esadecimali") {
    sprintf("di %d colori esadecimali", length(palette))
  } else {
    sprintf("'%s'", palette)
  }
}

# La riga per l'intestazione della tabella di un grafico che usa n colori
# della palette: vuota se la palette ne ha abbastanza.
nota_palette <- function(n) {
  massimo <- colori_della_palette()
  if (is.na(massimo) || n <= massimo) return(character())
  sprintf("%s (avviso E-ECO-46): il grafico usa %d colori, la palette %s ne ha %d",
          AVVISO_PALETTE, as.integer(n), nome_palette(), massimo)
}

# Imposta la palette di tutti i grafici dell'analisi.
imposta_palette <- function(palette) {
  assign("PALETTE_GRAFICI", as.character(unlist(palette)), envir = globalenv())
}

# n colori della palette in uso: i primi n, oppure, con piu' gruppi che colori,
# n colori interpolati lungo la palette (grDevices::colorRampPalette). Senza
# numeri casuali: gli stessi colori a ogni esecuzione.
colori_gruppi <- function(n) {
  if (n < 1L) return(character())
  palette <- PALETTE_GRAFICI
  estendi <- function(base) {
    if (n <= length(base)) base[seq_len(n)] else grDevices::colorRampPalette(base)(n)
  }
  switch(
    origine_palette(palette),
    hcl = grDevices::hcl.colors(max(n, 2L), palette = palette)[seq_len(n)],
    viridis = viridisLite::viridis(n, option = palette),
    okabe_ito = estendi(okabe_ito_senza_grigio()),
    brewer = {
      massimo <- RColorBrewer::brewer.pal.info[palette, "maxcolors"]
      if (n <= massimo) RColorBrewer::brewer.pal(max(n, 3L), palette)[seq_len(n)] else
        estendi(RColorBrewer::brewer.pal(massimo, palette))
    },
    esadecimali = estendi(palette),
    stop("palette non riconosciuta: ", paste(palette, collapse = ", "))
  )
}
