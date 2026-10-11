# Descrizione dell'oggetto di partenza, per la configurazione assistita.
#
# Parametri: oggetto (percorso della copia) e valori_mancanti. Scrive
# descrizione.json: non calcola nulla e non giudica, riporta cio' che l'oggetto
# contiene.
#
#   campioni            numero di campioni e di varianti
#   ranghi              per ogni rango tassonomico, la quota di varianti con
#                       un'assegnazione a quel rango
#   albero              se c'e', se e' radicato, se ha le lunghezze dei rami
#   colonne             per ogni colonna della tabella dei campioni: valori
#                       presenti, valori distinti, gruppi con almeno due
#                       campioni, se i valori sono tutti numerici

for (f in c("io_json.R", "errors.R")) source(file.path(Sys.getenv("AMPLICON16S_R_LIB"), f))
for (f in c("uscite.R", "campioni.R")) source(file.path(Sys.getenv("AMPLICON16S_ECO_R_LIB"), f))

esegui_fase(function(parametri, cartella) {
  richiedi_pacchetti(c("phyloseq", "ape"))
  letto <- leggi_oggetto(parametri$oggetto)
  if (is.null(letto$oggetto)) errore_catalogo("E-ECO-04", letto$motivo)
  ps <- letto$oggetto
  marcatori <- unlist(parametri$valori_mancanti)
  tabella <- tabella_campioni(ps)
  tassonomia <- methods::as(phyloseq::tax_table(ps), "matrix")
  albero <- phyloseq::access(ps, "phy_tree")

  ranghi <- lapply(phyloseq::rank_names(ps), function(rango) {
    valori <- tassonomia[, rango]
    list(rango = rango,
         # Il solo prefisso di rango ("g__") non e' un'assegnazione.
         quota_assegnata = mean(!is.na(valori) & nzchar(trimws(valori)) &
                                  !grepl("^[A-Za-z]__$", valori)))
  })
  colonne <- lapply(names(tabella), function(nome) {
    valori <- tabella[[nome]]
    presenti <- valori[!e_mancante(valori, marcatori)]
    dimensioni <- table(presenti)
    list(
      nome = nome,
      con_valore = length(presenti),
      distinti = length(dimensioni),
      gruppi_ripetuti = sum(dimensioni >= 2L),
      numerica = length(presenti) > 0L &&
        !anyNA(suppressWarnings(as.numeric(presenti)))
    )
  })
  scrivi_json(list(
    campioni = phyloseq::nsamples(ps),
    varianti = phyloseq::ntaxa(ps),
    ranghi = ranghi,
    albero = list(
      presente = !is.null(albero),
      radicato = if (!is.null(albero)) ape::is.rooted(albero) else NULL,
      con_lunghezze = if (!is.null(albero)) length(albero$edge.length) == nrow(albero$edge) else NULL
    ),
    colonne = colonne
  ), file.path(cartella, "descrizione.json"))
  "descrizione.json"
})
