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
# Con il parametro idtaxa (addestrato, soglia, filamento) il classificatore e'
# DECIPHER::IdTaxa (R/lib/idtaxa.R) invece di assignTaxonomy. tassonomia.tsv e
# copertura_per_gruppo.tsv hanno la stessa forma; al posto di bootstrap.tsv si
# scrive confidenza.tsv (la confidenza di ogni rango assegnato, 0-100, vuota
# dove il rango non e' assegnato), tassonomia.rds porta confidenza al posto di
# boot, e ranghi_scartati.json elenca i ranghi del classificatore che non hanno
# una colonna. Il seme si imposta e DECIPHER si carica solo in questo ramo: lo
# stato del generatore e i pacchetti caricati con assignTaxonomy non dipendono
# da questo ramo.
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

# idtaxa.R definisce solo funzioni: non carica DECIPHER e non tocca il
# generatore di numeri casuali.
for (f in c("io_json.R", "errors.R", "idtaxa.R")) {
  source(file.path(Sys.getenv("AMPLICON16S_R_LIB"), f))
}

# Il ramo IdTaxa: la tassonomia nella stessa forma di assignTaxonomy, con la
# confidenza al posto del bootstrap. Restituisce i nomi degli artefatti propri.
assegna_con_idtaxa <- function(parametri, cartella, varianti) {
  richiedi_pacchetti(c("DECIPHER", "Biostrings"))
  idtaxa <- parametri$idtaxa
  esito <- classifica_idtaxa(
    # Lo stesso controllo della validazione: vale per un classificatore
    # dichiarato dopo che la validazione era gia' conclusa.
    varianti, tryCatch(
      carica_addestrato(idtaxa$addestrato),
      error = function(e) errore_catalogo("E-S0-28", conditionMessage(e))
    ),
    soglia = as.numeric(idtaxa$soglia), filamento = idtaxa$filamento,
    seme = parametri$seme, processi = parametri$processi
  )
  ranghi <- colnames(esito$tax)
  salva_rds(list(tax = esito$tax, confidenza = esito$confidenza),
            file.path(cartella, "tassonomia.rds"))
  # Due decimali: la confidenza e' una media su cento estrazioni.
  scritta <- ifelse(is.na(esito$confidenza), "", sprintf("%.2f", esito$confidenza))
  dim(scritta) <- dim(esito$confidenza)
  scrivi_atomico(
    c(paste(c("sequenza", ranghi), collapse = "\t"),
      sprintf("%s\t%s", varianti, apply(scritta, 1, paste, collapse = "\t"))),
    file.path(cartella, "confidenza.tsv")
  )
  scrivi_json(list(ranghi_scartati = as.list(esito$scartati)),
              file.path(cartella, "ranghi_scartati.json"))
  list(tax = esito$tax, artefatti = c("confidenza.tsv", "ranghi_scartati.json"))
}

esegui_fase(function(parametri, cartella) {
  con_idtaxa <- !is.null(facoltativo(parametri, "idtaxa"))
  if (!con_idtaxa) richiedi_pacchetti("dada2")

  tabella <- readRDS(parametri$tabella)
  varianti <- colnames(tabella)
  if (length(varianti) == 0) {
    stop("la tabella di S7 non ha varianti: S7 avrebbe dovuto fermarsi con E-S7-01")
  }

  if (con_idtaxa) {
    assegnata <- assegna_con_idtaxa(parametri, cartella, varianti)
    tax <- assegnata$tax
    ranghi <- colnames(tax)
  } else {
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
  }

  letture <- colSums(tabella)
  vuoto <- function(x) ifelse(is.na(x), "", x)
  scrivi_atomico(
    c(paste(c("sequenza", ranghi, "letture"), collapse = "\t"),
      sprintf("%s\t%s\t%.0f", varianti,
              apply(tax, 1, function(r) paste(vuoto(r), collapse = "\t")), letture)),
    file.path(cartella, "tassonomia.tsv")
  )
  if (!con_idtaxa) {
    scrivi_atomico(
      c(paste(c("sequenza", ranghi), collapse = "\t"),
        sprintf("%s\t%s", varianti, apply(boot, 1, paste, collapse = "\t"))),
      file.path(cartella, "bootstrap.tsv")
    )
  }

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

  if (con_idtaxa) {
    return(c("tassonomia.rds", "tassonomia.tsv", assegnata$artefatti,
             "copertura_per_gruppo.tsv"))
  }
  c("tassonomia.rds", "tassonomia.tsv", "bootstrap.tsv", "copertura_per_gruppo.tsv")
})
