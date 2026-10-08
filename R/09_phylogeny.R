# S9 - filogenesi delle varianti finali.
#
# Scrive in 09_phylogeny/:
#
#   albero.nwk          l'albero radicato, in formato Newick: le foglie sono gli
#                       identificativi delle varianti
#   allineamento.fasta  l'allineamento multiplo su cui e' stato costruito
#   filogenesi.json     che cosa e' stato fatto: allineatore, modello, ricerca,
#                       radicamento, verosimiglianza e parametri stimati
#
# Parametri: filtrato (ps_filtrato.rds di S13), allineatore, modello, seme,
# processi.
#
# SU CHE COSA. Le sequenze di riferimento delle varianti dell'oggetto filtrato
# di S13: quelle che restano dopo decontaminazione e filtri, e che S14
# consegna. La guardia sul loro numero (phylo.max_seqs, E-S9-01) e' nella fase
# Python, prima che questo processo parta.
#
# ALLINEAMENTO (phylo.aligner = decipher): DECIPHER::AlignSeqs con le
# impostazioni predefinite. Le sequenze sono lette troncate alla stessa
# lunghezza, ma inserzioni e delezioni fra taxa le sfasano: l'allineamento le
# rimette in colonna.
#
# ALBERO (phylo.model = GTR+G+I), con phangorn:
#   1. distanze di massima verosimiglianza fra le sequenze allineate (dist.ml,
#      modello JC69) e albero di partenza neighbor-joining;
#   2. massima verosimiglianza con il modello GTR, eterogeneita' gamma fra i
#      siti in CATEGORIE_GAMMA categorie e una quota di siti invarianti:
#      optim.pml ottimizza frequenze delle basi, tassi di sostituzione, forma
#      della gamma, quota di invarianti, lunghezze dei rami e topologia.
#
# LA RICERCA E' DETERMINISTICA, PER SCELTA. La topologia si migliora con scambi
# fra rami vicini (rearrangement = "NNI") a partire dall'albero
# neighbor-joining, finche' la verosimiglianza non cresce piu'. L'alternativa
# di phangorn, rearrangement = "stochastic" (o "ratchet"), perturba l'albero
# con mosse casuali per uscire dai massimi locali: puo' trovare un albero con
# verosimiglianza un poco piu' alta, ma il risultato dipende dal generatore di
# numeri casuali e costa molte volte tanto. Qui conta che la stessa
# configurazione dia lo stesso albero: nessun passo di questa ricerca estrae
# numeri casuali. Il seme (run.seed) si fissa comunque prima di ogni calcolo,
# cosi' che una versione futura delle librerie che ne introducesse uno resti
# riproducibile. Il numero di processi riguarda solo l'allineamento, e non ne
# cambia il risultato: l'albero e' identico byte per byte fra due esecuzioni e
# con un numero di thread diverso (verificato dai test).
#
# RADICAMENTO: AL PUNTO MEDIO (phangorn::midpoint). La massima verosimiglianza
# con un modello reversibile da' un albero senza radice, e le misure di
# diversita' che usano la filogenesi (UniFrac, la diversita' filogenetica di
# Faith) ne richiedono uno radicato. Fra le varianti non c'e' un gruppo esterno
# su cui radicare: i filtri tolgono cio' che non e' batterico, e scegliere una
# variante a caso come radice, come fa phyloseq quando riceve un albero non
# radicato, renderebbe le distanze dipendenti da quella scelta. Il punto medio
# del cammino piu' lungo fra due foglie e' una regola che dipende solo
# dall'albero. Assume tassi di evoluzione simili nei due rami principali: e'
# un'approssimazione, dichiarata in filogenesi.json.
#
# Le foglie sono numerate nell'ordine degli identificativi delle varianti
# (R/lib/albero.R).

for (f in c("io_json.R", "errors.R", "albero.R")) {
  source(file.path(Sys.getenv("AMPLICON16S_R_LIB"), f))
}

# Categorie della distribuzione gamma discreta dei tassi fra i siti, e quota di
# siti invarianti da cui parte l'ottimizzazione. La ricerca e' una salita
# locale: il punto di arrivo dipende da quello di partenza, che per questo e'
# una costante dichiarata.
CATEGORIE_GAMMA <- 4L
INVARIANTI_INIZIALI <- 0.2
# Sotto questo numero di sequenze un albero radicato non ha topologia da
# stimare. La fase Python lo verifica prima di avviare questo processo
# (E-S9-02): qui resta come invariante.
VARIANTI_MINIME <- 4L

arrotonda <- function(x) signif(as.numeric(x), CIFRE_NEWICK)

esegui_fase(function(parametri, cartella) {
  richiedi_pacchetti(c("phyloseq", "Biostrings", "DECIPHER", "phangorn", "ape"))
  if (!identical(parametri$allineatore, "decipher")) {
    stop("phylo.aligner non realizzato: ", parametri$allineatore)
  }
  if (!identical(parametri$modello, "GTR+G+I")) {
    stop("phylo.model non realizzato: ", parametri$modello)
  }
  ps <- readRDS(parametri$filtrato)
  varianti <- phyloseq::taxa_names(ps)
  sequenze <- phyloseq::refseq(ps)
  if (!identical(names(sequenze), varianti)) {
    stop("le sequenze di riferimento non sono allineate alle varianti dell'oggetto")
  }
  if (length(varianti) < VARIANTI_MINIME) {
    stop(sprintf("servono almeno %d varianti per costruire un albero, trovate %d",
                 VARIANTI_MINIME, length(varianti)))
  }

  set.seed(as.integer(parametri$seme))
  allineamento <- DECIPHER::AlignSeqs(
    sequenze, processors = as.integer(parametri$processi), verbose = FALSE
  )
  dati <- phangorn::phyDat(methods::as(allineamento, "matrix"), type = "DNA")
  partenza <- phangorn::NJ(phangorn::dist.ml(dati))
  adattato <- phangorn::pml(partenza, data = dati, model = "GTR",
                            k = CATEGORIE_GAMMA, inv = INVARIANTI_INIZIALI)
  adattato <- phangorn::optim.pml(
    adattato, model = "GTR", optGamma = TRUE, optInv = TRUE, optNni = TRUE,
    rearrangement = "NNI", control = phangorn::pml.control(trace = 0L)
  )
  albero <- foglie_in_ordine(phangorn::midpoint(adattato$tree), varianti)
  albero$node.label <- NULL

  # Il file e' la forma dell'albero: cio' che si verifica e' cio' che si
  # rilegge da li', con le cifre che il file conserva.
  testo <- testo_newick(albero)
  scrivi_atomico(testo, file.path(cartella, "albero.nwk"))
  riletto <- leggi_newick(file.path(cartella, "albero.nwk"))
  difetti <- difetti_albero(riletto, varianti)
  if (!identical(testo_newick(riletto), testo)) {
    difetti <- c(difetti, "albero.nwk riletto e riscritto non da' lo stesso testo")
  }
  if (length(difetti) > 0L) stop("albero non valido: ", paste(difetti, collapse = "; "))

  scrivi_atomico(c(rbind(paste0(">", varianti), unname(as.character(allineamento)))),
                 file.path(cartella, "allineamento.fasta"))
  tassi <- adattato$Q
  names(tassi) <- c("AC", "AG", "AT", "CG", "CT", "GT")
  frequenze <- adattato$bf
  names(frequenze) <- c("A", "C", "G", "T")
  scrivi_json(
    list(
      varianti = length(varianti),
      allineatore = "decipher: DECIPHER::AlignSeqs, impostazioni predefinite",
      colonne_allineamento = unique(Biostrings::width(allineamento)),
      modello = parametri$modello,
      categorie_gamma = CATEGORIE_GAMMA,
      albero_di_partenza = "neighbor-joining su distanze di massima verosimiglianza (JC69)",
      ricerca = "scambi fra rami vicini (NNI) dall'albero di partenza, senza componenti casuali",
      seme = as.integer(parametri$seme),
      radicamento = paste(
        "punto medio del cammino piu' lungo fra due foglie (phangorn::midpoint):",
        "nessun gruppo esterno fra le varianti; assume tassi simili nei due rami principali"
      ),
      radicato = ape::is.rooted(riletto),
      binario = ape::is.binary(riletto),
      log_verosimiglianza = arrotonda(adattato$logLik),
      forma_gamma = arrotonda(adattato$shape),
      siti_invarianti = arrotonda(adattato$inv),
      frequenze_basi = as.list(arrotonda(frequenze) |> stats::setNames(names(frequenze))),
      tassi_sostituzione = as.list(arrotonda(tassi) |> stats::setNames(names(tassi))),
      lunghezza_totale = arrotonda(sum(riletto$edge.length)),
      cifre_newick = CIFRE_NEWICK
    ),
    file.path(cartella, "filogenesi.json")
  )
  c("albero.nwk", "allineamento.fasta", "filogenesi.json")
})
