# Classificazione tassonomica con DECIPHER::IdTaxa su un classificatore
# addestrato (l'oggetto di DECIPHER::LearnTaxa).
#
# La usano S8 (R/08_taxonomy.R), quando il classificatore dichiarato e' IdTaxa,
# e la verifica del classificatore prima dell'esecuzione
# (R/00_classificatore.R). Legge il file del classificatore; non scrive nulla.
#
# Richiede DECIPHER e Biostrings. Interfaccia:
#   carica_addestrato(percorso)  -> l'oggetto addestrato; si ferma con un
#                                   messaggio che comincia con
#                                   "classificatore non valido:" se il file non
#                                   ne contiene uno con la tabella dei ranghi
#   ranghi_dell_addestrato(a)    -> i ranghi che il classificatore nomina
#   classifica_idtaxa(varianti, addestrato, soglia, filamento, seme, processi)
#                                -> list(tax, confidenza, scartati)
#
# I ranghi si riconoscono per nome, non per posizione: il classificatore di
# SILVA per DECIPHER comincia con Root (rango rootrank) e nomina ogni rango, e
# in alcuni lignaggi ne salta. Le colonne sono sempre le sei della pipeline, da
# Kingdom a Genus; un rango che il classificatore nomina e che non ha una
# colonna (species, subclass, ...) si scarta, e chi chiama lo dichiara. Un rango
# intermedio che il lignaggio non ha resta NA con ranghi assegnati sotto.

COLONNE_DEI_RANGHI <- list(
  Kingdom = c("domain", "kingdom", "superkingdom"),
  Phylum = "phylum",
  Class = "class",
  Order = "order",
  Family = "family",
  Genus = "genus"
)

# Il rango della radice dell'albero dei taxa: non e' un'assegnazione.
RANGO_RADICE <- "rootrank"

# Il classificatore addestrato letto da un file .RData. Il file si carica in un
# ambiente proprio, per non toccare quello globale, e deve contenere un solo
# oggetto di classe Train, con la tabella dei ranghi e almeno un rango fra
# quelli della tassonomia. Ogni rifiuto ha un messaggio che comincia con
# "classificatore non valido:", da cui chi chiama lo riconosce.
carica_addestrato <- function(percorso) {
  ambiente <- new.env()
  nomi <- tryCatch(
    load(percorso, envir = ambiente),
    error = function(e) stop(
      "classificatore non valido: il file non si carica come .RData (",
      conditionMessage(e), ")", call. = FALSE)
  )
  addestrati <- Filter(function(n) inherits(get(n, envir = ambiente), "Train"), nomi)
  if (length(addestrati) != 1L) {
    stop(sprintf(paste0(
      "classificatore non valido: il file contiene %d oggetti di DECIPHER::LearnTaxa ",
      "(classe Train) fra %s: ne serve uno"), length(addestrati),
      if (length(nomi)) paste(nomi, collapse = ", ") else "nessun oggetto"), call. = FALSE)
  }
  addestrato <- get(addestrati, envir = ambiente)
  ranghi <- addestrato$ranks
  if (is.null(ranghi) || !length(ranghi)) {
    stop(paste0(
      "classificatore non valido: l'oggetto addestrato non ha la tabella dei ranghi ",
      "(LearnTaxa senza l'argomento rank): senza i nomi dei ranghi le assegnazioni ",
      "non si possono riportare alle colonne della tassonomia"), call. = FALSE)
  }
  if (!any(tolower(ranghi) %in% unlist(COLONNE_DEI_RANGHI))) {
    stop(sprintf(paste0(
      "classificatore non valido: nessuno dei ranghi del classificatore (%s) e' fra ",
      "quelli della tassonomia (%s)"),
      paste(unique(ranghi), collapse = ", "),
      paste(unlist(COLONNE_DEI_RANGHI), collapse = ", ")), call. = FALSE)
  }
  addestrato
}

# I ranghi che il classificatore nomina, in minuscolo, senza ripetizioni, in
# ordine indipendente dalle impostazioni locali.
ranghi_dell_addestrato <- function(addestrato) {
  sort(unique(tolower(as.character(addestrato$ranks))), method = "radix")
}

# La classificazione delle varianti con DECIPHER::IdTaxa. varianti sono le
# sequenze, soglia la confidenza minima (0-100), filamento l'orientamento da
# provare, seme e processi quelli dell'esecuzione. Restituisce list(tax,
# confidenza, scartati): due matrici varianti x colonne della tassonomia (il
# nome del taxon e la sua confidenza, NA dove non assegnato) e i ranghi del
# classificatore che non hanno una colonna. Passi: seme, classificazione, poi
# per ogni variante e per ogni colonna il primo dei ranghi equivalenti che la
# classificazione porta, cercato per nome.
classifica_idtaxa <- function(varianti, addestrato, soglia, filamento, seme, processi) {
  # Il bootstrap di IdTaxa campiona i k-meri con il generatore di R.
  set.seed(as.integer(seme))
  esito <- DECIPHER::IdTaxa(
    Biostrings::DNAStringSet(varianti), addestrato,
    strand = filamento, threshold = soglia,
    processors = as.integer(processi), verbose = FALSE
  )
  colonne <- names(COLONNE_DEI_RANGHI)
  tax <- matrix(NA_character_, nrow = length(varianti), ncol = length(colonne),
                dimnames = list(varianti, colonne))
  confidenza <- matrix(NA_real_, nrow = length(varianti), ncol = length(colonne),
                       dimnames = list(varianti, colonne))
  for (i in seq_along(esito)) {
    voce <- esito[[i]]
    del_rango <- tolower(voce$rank)
    for (colonna in colonne) {
      posto <- match(COLONNE_DEI_RANGHI[[colonna]], del_rango)
      posto <- posto[!is.na(posto)][1L]
      if (is.na(posto)) next
      nome <- voce$taxon[posto]
      # IdTaxa chiude la classificazione troncata dalla soglia con
      # un'etichetta "unclassified_<ultimo taxon>": non e' un taxon.
      if (is.na(nome) || startsWith(nome, "unclassified_")) next
      tax[i, colonna] <- nome
      confidenza[i, colonna] <- voce$confidence[posto]
    }
  }
  noti <- c(RANGO_RADICE, unlist(COLONNE_DEI_RANGHI))
  list(tax = tax, confidenza = confidenza,
       scartati = setdiff(ranghi_dell_addestrato(addestrato), noti))
}
