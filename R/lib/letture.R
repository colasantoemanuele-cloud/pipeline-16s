# Tracciamento delle letture per campione, passo dopo passo.
#
# La usano gli script delle fasi che registrano un conteggio di letture per
# campione; la tabella complessiva la compone la parte Python
# (runner/tracciamento.py).
#
# Richiede io_json.R. Interfaccia:
#   traccia_letture(conteggi, passo, cartella) -> nome del file scritto
#   leggi_letture(percorso)                    -> vettore intero con nome
#
# Ogni fase registra quante letture restano a ciascun campione dopo il proprio
# passo, in un file della propria cartella: letture_<passo>.tsv, con le
# colonne campione, passo, letture. Un file per passo e non uno condiviso: la
# fase resta proprietaria dei suoi artefatti, il file entra nel suo manifesto
# e una ripresa che riesegue la fase lo sostituisce senza toccare quelli delle
# altre. La tabella complessiva si ottiene accostando i file.
#
# Il nome restituito va incluso negli artefatti dichiarati dalla fase.

# Scrive il file del tracciamento di un passo. conteggi e' un vettore con nome
# (campione -> letture), passo il nome del passo, in minuscole, cifre e
# trattino basso perche' entra nel nome del file. Conteggi non interi, negativi
# o mancanti, e nomi vuoti, ripetuti o con tabulazioni fermano la scrittura:
# una tabella malformata romperebbe il tracciamento di tutta l'esecuzione. I
# conteggi si scrivono senza notazione esponenziale.
traccia_letture <- function(conteggi, passo, cartella) {
  if (!is.character(passo) || length(passo) != 1L ||
      !grepl("^[a-z0-9_]+$", passo)) {
    stop("nome di passo non valido: ", paste(passo, collapse = ", "))
  }
  if (!is.numeric(conteggi) || anyNA(conteggi) || any(conteggi < 0) ||
      any(conteggi != round(conteggi))) {
    stop("i conteggi devono essere interi non negativi")
  }
  campioni <- names(conteggi)
  if (is.null(campioni) || any(!nzchar(campioni)) || anyDuplicated(campioni) ||
      any(grepl("[\t\n]", campioni))) {
    stop("ogni conteggio deve avere il nome univoco del suo campione")
  }

  nome <- sprintf("letture_%s.tsv", passo)
  righe <- c(
    "campione\tpasso\tletture",
    sprintf("%s\t%s\t%.0f", campioni, passo, conteggi)
  )
  scrivi_atomico(righe, file.path(cartella, nome))
  nome
}

# Il file del tracciamento di un passo come vettore intero con nome (campione
# -> letture).
leggi_letture <- function(percorso) {
  tabella <- utils::read.delim(
    percorso,
    colClasses = c("character", "character", "numeric"),
    quote = "",
    comment.char = ""
  )
  stats::setNames(as.integer(tabella$letture), tabella$campione)
}
