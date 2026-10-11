# S8 - assegnazione tassonomica (dada2::assignTaxonomy).
#
# Scrive in 08_taxonomy/:
#
#   tassonomia.tsv          per variante: sequenza, un rango per colonna
#                           (vuoto dove il bootstrap e' sotto min_boot), letture
#                           nella tabella di S7
#   bootstrap.tsv           per variante: il bootstrap di ciascun rango, 0-100,
#                           anche dove il rango non e' assegnato
#   tassonomia.rds          lista con le matrici tax e boot di assignTaxonomy,
#                           righe nell'ordine delle colonne della tabella
#   copertura_per_gruppo.tsv  per gruppo di campioni: varianti presenti e
#                           quante hanno il phylum, letture e quante sono di
#                           varianti con il phylum
#
# Parametri: tabella (tabella_asv.rds di S7), riferimento, min_boot, try_rc,
# seme, processi, gruppi (accession -> gruppo).
#
# Riproducibilita'. Il bootstrap di assignTaxonomy estrae i k-meri con
# runif di R, prima del calcolo parallelo: fissato il seme con set.seed, il
# sottocampionamento non dipende dal numero di thread. Nella versione
# ufficiale di dada2 resta casuale, con un generatore che il seme non controlla
# (std::random_device), la scelta fra generi a pari probabilita' in
# get_best_genus (src/taxonomy.cpp): la versione corretta dell'immagine usa un
# generatore con seme dal seme di R, dall'indice della sequenza e dal caso.
# Quale delle due R carichi lo accerta la fase prima di chiamare questo script
# (steps/s08_taxonomy.py): senza la correzione l'assegnazione nei pareggi fra
# generi non e' ripetibile, e la fase lo dichiara (E-S8-03).

for (f in c("io_json.R", "errors.R")) {
  source(file.path(Sys.getenv("AMPLICON16S_R_LIB"), f))
}

esegui_fase(function(parametri, cartella) {
  richiedi_pacchetti("dada2")

  tabella <- readRDS(parametri$tabella)
  varianti <- colnames(tabella)
  if (length(varianti) == 0) {
    stop("la tabella di S7 non ha varianti: S7 avrebbe dovuto fermarsi con E-S7-01")
  }

  set.seed(as.integer(parametri$seme))
  esito <- dada2::assignTaxonomy(
    varianti, parametri$riferimento,
    minBoot = as.integer(parametri$min_boot),
    tryRC = isTRUE(parametri$try_rc),
    outputBootstraps = TRUE,
    multithread = as.integer(parametri$processi),
    verbose = FALSE
  )
  tax <- esito$tax
  boot <- esito$boot
  ranghi <- colnames(tax)
  stopifnot(identical(rownames(tax), varianti), identical(colnames(boot), ranghi))
  # assignTaxonomy nomina i ranghi per posizione: il secondo e' il phylum, e
  # oltre il settimo non c'e' un nome. La validazione respinge un riferimento
  # che non li rispetta guardandone le prime intestazioni: qui si verifica il
  # risultato.
  if (anyNA(ranghi) || !("Phylum" %in% ranghi)) {
    stop(sprintf(paste0(
      "la tassonomia assegnata ha i ranghi %s: il riferimento deve avere da due a ",
      "sette ranghi per sequenza (regno, phylum, ... specie), separati da punto e virgola"),
      paste(ifelse(is.na(ranghi), "(senza nome)", ranghi), collapse = ", ")))
  }
  salva_rds(list(tax = tax, boot = boot), file.path(cartella, "tassonomia.rds"))

  letture <- colSums(tabella)
  vuoto <- function(x) ifelse(is.na(x), "", x)
  scrivi_atomico(
    c(paste(c("sequenza", ranghi, "letture"), collapse = "\t"),
      sprintf("%s\t%s\t%.0f", varianti,
              apply(tax, 1, function(r) paste(vuoto(r), collapse = "\t")), letture)),
    file.path(cartella, "tassonomia.tsv")
  )
  scrivi_atomico(
    c(paste(c("sequenza", ranghi), collapse = "\t"),
      sprintf("%s\t%s", varianti, apply(boot, 1, paste, collapse = "\t"))),
    file.path(cartella, "bootstrap.tsv")
  )

  # La copertura per gruppo: le varianti presenti nei campioni del gruppo, e
  # le letture del gruppo, con e senza il phylum assegnato.
  # Un phylum scritto col solo prefisso di rango ("p__") e' un phylum che la
  # sequenza del riferimento non ha: non conta come assegnato.
  con_phylum <- !is.na(tax[, "Phylum"]) & !grepl("^[A-Za-z]__$", tax[, "Phylum"])
  gruppi <- unlist(parametri$gruppi)[rownames(tabella)]
  righe <- vapply(sort(unique(gruppi), method = "radix"), function(g) {
    parziale <- colSums(tabella[gruppi == g, , drop = FALSE])
    presenti <- parziale > 0
    sprintf("%s\t%d\t%d\t%.0f\t%.0f", g, sum(presenti), sum(presenti & con_phylum),
            sum(parziale), sum(parziale[con_phylum]))
  }, character(1))
  scrivi_atomico(
    c("gruppo\tvarianti\tvarianti_con_phylum\tletture\tletture_con_phylum", righe),
    file.path(cartella, "copertura_per_gruppo.tsv")
  )

  c("tassonomia.rds", "tassonomia.tsv", "bootstrap.tsv", "copertura_per_gruppo.tsv")
})
