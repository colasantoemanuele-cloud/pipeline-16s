# L'albero filogenetico delle varianti: scrittura, lettura e verifica.
#
# Richiede ape. Interfaccia:
#   testo_newick(albero)             l'albero in formato Newick, una riga
#   leggi_newick(percorso)           l'albero letto da un file Newick
#   foglie_in_ordine(albero, ordine) lo stesso albero con le foglie numerate
#                                    nell'ordine indicato
#   difetti_albero(albero, varianti) che cosa non va nell'albero rispetto alle
#                                    varianti dell'oggetto; vuoto se nulla
#
# IL FORMATO NEWICK E' L'UNICA FORMA DELL'ALBERO che passa da una fase
# all'altra e che viene consegnata: S9 lo scrive, S14 lo rilegge per metterlo
# nell'oggetto finale e lo esporta. Le lunghezze dei rami hanno CIFRE_NEWICK
# cifre significative, scritte con sprintf("%.8g"): la forma non dipende dalle
# impostazioni locali (il separatore decimale e' sempre il punto) e rileggere e
# riscrivere il file da' lo stesso testo. L'albero nell'oggetto finale e'
# quindi esattamente quello del file, non una versione con piu' cifre che
# l'export non saprebbe ricostruire.
#
# LE FOGLIE portano gli identificativi delle varianti, non le sequenze.
#
# L'ORDINE DELLE FOGLIE. phyloseq, quando un oggetto ha un albero, riordina le
# varianti di ogni componente nell'ordine delle foglie. Gli identificativi
# delle varianti hanno un ordine proprio (out.asv_id_scheme), che gli export
# rispettano: le foglie si numerano quindi in quell'ordine prima di unire
# l'albero all'oggetto, e la topologia non cambia.

CIFRE_NEWICK <- 8L

testo_newick <- function(albero) {
  ape::write.tree(albero, file = "", digits = CIFRE_NEWICK)
}

leggi_newick <- function(percorso) {
  ape::read.tree(file = percorso)
}

foglie_in_ordine <- function(albero, ordine) {
  n <- length(albero$tip.label)
  if (!setequal(albero$tip.label, ordine) || n != length(ordine)) {
    stop("le foglie dell'albero non sono le varianti indicate")
  }
  nuovo <- match(albero$tip.label, ordine)
  foglie <- albero$edge <= n
  albero$edge[foglie] <- nuovo[albero$edge[foglie]]
  albero$tip.label <- as.character(ordine)
  albero
}

difetti_albero <- function(albero, varianti) {
  difetti <- character()
  foglie <- albero$tip.label
  if (anyDuplicated(foglie)) difetti <- c(difetti, "foglie ripetute nell'albero")
  mancanti <- setdiff(varianti, foglie)
  estranee <- setdiff(foglie, varianti)
  if (length(mancanti) > 0L || length(estranee) > 0L) {
    difetti <- c(difetti, sprintf(
      "le foglie dell'albero non coincidono con le varianti: %d varianti senza foglia, %d foglie che non sono varianti",
      length(mancanti), length(estranee)
    ))
  }
  if (!ape::is.rooted(albero)) difetti <- c(difetti, "albero non radicato")
  if (is.null(albero$edge.length) || anyNA(albero$edge.length) || any(albero$edge.length < 0)) {
    difetti <- c(difetti, "lunghezze dei rami assenti o negative")
  }
  difetti
}
