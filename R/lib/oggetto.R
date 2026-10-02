# Verifiche dell'oggetto integrato (S10), sugli identificativi.
#
# Richiede errors.R. Interfaccia:
#   verifica_orientamento(matrice, varianti, campioni, taxa_are_rows, cosa)
#       la matrice ha le varianti e i campioni dove out.taxa_are_rows dice,
#       riconosciuti dai nomi di righe e colonne; altrimenti E-S10-01
#   verifica_nomi(trovati, attesi, cosa)
#       i nomi sono esattamente quelli attesi, nello stesso ordine;
#       altrimenti E-S10-01
#   verifica_sintattici(nomi, cosa)
#       nomi che R non altera costruendo un data.frame; altrimenti E-S10-01
#   ordine_varianti(letture_biologici, letture_totali, sequenze)
#       l'ordine delle varianti per gli identificativi (out.asv_id_scheme =
#       abundance_rank): letture nei biologici decrescenti, a parita' letture
#       totali decrescenti, a parita' la sequenza in ordine lessicografico
#
# Il controllo dell'orientamento non guarda le dimensioni. Una matrice con
# tanti campioni quante varianti, trasposta, ha le stesse dimensioni di quella
# giusta: un controllo su nrow e ncol la lascerebbe passare, e ogni calcolo a
# valle darebbe risultati plausibili e sbagliati. I nomi invece non mentono:
# le righe devono essere gli identificativi delle varianti, le colonne gli
# accession, o il contrario con taxa_are_rows falso.

.descrivi_nomi <- function(x, quanti = 3L) {
  if (is.null(x)) return("nessun nome")
  testa <- utils::head(x, quanti)
  paste0(paste(testa, collapse = ", "), if (length(x) > quanti) ", ..." else "")
}

verifica_nomi <- function(trovati, attesi, cosa) {
  if (identical(as.character(trovati), as.character(attesi))) return(invisible(TRUE))
  mancanti <- setdiff(attesi, trovati)
  estranei <- setdiff(trovati, attesi)
  dettaglio <- if (length(mancanti) == 0L && length(estranei) == 0L &&
                   length(trovati) == length(attesi)) {
    "stessi nomi in un ordine diverso"
  } else {
    sprintf("%d attesi, %d trovati; mancanti: %s; estranei: %s",
            length(attesi), length(trovati),
            .descrivi_nomi(mancanti), .descrivi_nomi(estranei))
  }
  errore_catalogo("E-S10-01", sprintf("%s: %s", cosa, dettaglio))
}

verifica_orientamento <- function(matrice, varianti, campioni, taxa_are_rows, cosa) {
  righe <- rownames(matrice)
  colonne <- colnames(matrice)
  attese_righe <- if (isTRUE(taxa_are_rows)) varianti else campioni
  attese_colonne <- if (isTRUE(taxa_are_rows)) campioni else varianti
  if (identical(righe, attese_righe) && identical(colonne, attese_colonne)) {
    return(invisible(TRUE))
  }
  if (identical(righe, attese_colonne) && identical(colonne, attese_righe)) {
    errore_catalogo("E-S10-01", sprintf(
      "%s trasposta: le righe sono %s e le colonne %s, ma out.taxa_are_rows e' %s",
      cosa,
      if (isTRUE(taxa_are_rows)) "i campioni" else "le varianti",
      if (isTRUE(taxa_are_rows)) "le varianti" else "i campioni",
      if (isTRUE(taxa_are_rows)) "vero" else "falso"
    ))
  }
  verifica_nomi(righe, attese_righe, paste(cosa, "- righe"))
  verifica_nomi(colonne, attese_colonne, paste(cosa, "- colonne"))
}

verifica_sintattici <- function(nomi, cosa) {
  alterati <- nomi[make.names(nomi, unique = TRUE) != nomi]
  if (length(alterati) > 0L) {
    errore_catalogo("E-S10-01", sprintf(
      "%s: nomi che R altererebbe: %s", cosa, .descrivi_nomi(alterati, 10L)
    ))
  }
  invisible(TRUE)
}

# L'ultimo criterio rende l'ordine totale: due varianti non hanno mai la stessa
# sequenza. method = "radix" confronta le stringhe byte per byte, come nella
# localizzazione C, qualunque sia la localizzazione del processo: con
# l'ordinamento della localizzazione la numerazione potrebbe cambiare da una
# macchina all'altra.
ordine_varianti <- function(letture_biologici, letture_totali, sequenze) {
  if (anyDuplicated(sequenze)) stop("sequenze ripetute fra le varianti")
  order(-as.numeric(letture_biologici), -as.numeric(letture_totali), sequenze,
        method = "radix")
}
