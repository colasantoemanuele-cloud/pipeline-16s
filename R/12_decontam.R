# S12 - decontaminazione dai controlli negativi (decontam, per prevalenza).
#
# Scrive in 11_controls/, condivisa con S11:
#
#   decontam_varianti.tsv      per variante: genere, letture e prevalenza nei
#                              negativi, nei biologici e nei positivi, probabilita'
#                              aggregata e per piastra, contaminante in ciascuna
#                              modalita', e se e' rimossa
#   decontam_per_piastra.tsv   per piastra e variante: il confronto usato, la
#                              probabilita' e l'esito
#   decontam_contaminanti.tsv  i contaminanti rimossi, con punteggio e prevalenze
#   decontam_riepilogo.json    le due modalita' a confronto, quella dichiarata e
#                              la frazione rimossa per classe
#   ps_decontaminato.rds       l'oggetto integrato senza i contaminanti: gli
#                              identificativi delle altre varianti non cambiano
#   letture_decontaminate.tsv  letture per campione dopo la rimozione (tracciamento)
#
# Parametri: oggetto (ps_integrato.rds di S10), soglia (decontam.threshold),
# combinazione (decontam.batch_combine), min_negativi (decontam.min_blanks),
# modalita (decontam.mode: aggregate o batch), max_frazione
# (qc.max_frac_contaminant).
#
# CHI ENTRA NEL CONFRONTO. I controlli negativi sono il termine di confronto,
# i campioni biologici l'altro termine. I controlli positivi ne restano fuori:
# contengono un organismo aggiunto e, ai livelli diluiti, contaminanti di
# reagente amplificati, e confonderebbero la prevalenza in entrambi i sensi.
# La classe e' quella dei metadati: i campioni classificati biologici entrano
# come biologici.
#
# LE DUE MODALITA'. Aggregata: tutti i biologici contro tutti i negativi. Per
# piastra: i biologici di ogni piastra contro i negativi della stessa piastra,
# e le probabilita' delle piastre combinate con decontam.batch_combine. Una
# piastra con meno di decontam.min_blanks negativi, o un campione senza
# piastra, non ha un confronto proprio: i suoi biologici si confrontano con i
# negativi di tutte le piastre, e quel confronto entra nella combinazione come
# quello di una piastra. Una variante e' un contaminante se la probabilita' e'
# minore (strettamente) di decontam.threshold.
#
# LA MODALITA'. Decide decontam.mode, dichiarata in configurazione, non
# l'esito: scegliere dopo aver visto quanto si rimuove cambierebbe metodo in
# silenzio da un dataset all'altro. L'altra modalita' si calcola comunque e si
# registra come diagnostica. La frazione delle letture dei biologici che la
# modalita' dichiarata rimuove e' il controllo: oltre qc.max_frac_contaminant
# la fase si ferma con E-S12-02 (lo dichiara la fase Python, dopo che questo
# script ha scritto le misure).

for (f in c("io_json.R", "errors.R", "letture.R", "decontaminazione.R")) {
  source(file.path(Sys.getenv("AMPLICON16S_R_LIB"), f))
}

esegui_fase(function(parametri, cartella) {
  richiedi_pacchetti(c("phyloseq", "decontam"))
  soglia <- as.numeric(parametri$soglia)
  regola <- parametri$combinazione
  min_negativi <- as.integer(parametri$min_negativi)
  max_frazione <- as.numeric(parametri$max_frazione)
  modalita <- parametri$modalita

  ps <- readRDS(parametri$oggetto)
  dati <- methods::as(phyloseq::sample_data(ps), "data.frame")
  conteggi <- methods::as(phyloseq::otu_table(ps), "matrix")
  if (phyloseq::taxa_are_rows(ps)) conteggi <- t(conteggi)
  storage.mode(conteggi) <- "numeric"
  tax <- methods::as(phyloseq::tax_table(ps), "matrix")
  varianti <- colnames(conteggi)
  if (!identical(rownames(dati), rownames(conteggi))) {
    stop("l'oggetto di S10 non ha campioni allineati fra conteggi e metadati")
  }
  classe <- dati$classe
  piastra <- dati$piastra
  negativo <- classe == "controllo_negativo"
  biologico <- classe == "biologico"
  positivo <- classe == "controllo_positivo"
  con_letture <- rowSums(conteggi) > 0
  nel_confronto <- (negativo | biologico) & con_letture

  # ---- Aggregata ------------------------------------------------------------
  righe <- which(nel_confronto)
  p_aggregata <- probabilita_prevalenza(conteggi[righe, , drop = FALSE], negativo[righe])

  # ---- Per piastra ----------------------------------------------------------
  piastre <- sort(unique(piastra[nel_confronto & !is.na(piastra)]), method = "radix")
  confronti <- list()
  descrizione <- list()
  for (p in piastre) {
    della_piastra <- nel_confronto & !is.na(piastra) & piastra == p
    n_neg <- sum(della_piastra & negativo)
    n_bio <- sum(della_piastra & biologico)
    if (n_bio == 0L) next
    if (n_neg >= min_negativi) {
      quali <- which(della_piastra)
      come <- "negativi della piastra"
    } else {
      quali <- which((della_piastra & biologico) | (nel_confronto & negativo))
      come <- sprintf("negativi di tutte le piastre: %d nella piastra, meno di decontam.min_blanks (%d)",
                      n_neg, min_negativi)
    }
    confronti[[p]] <- probabilita_prevalenza(conteggi[quali, , drop = FALSE], negativo[quali])
    descrizione[[p]] <- list(negativi = sum(negativo[quali]), biologici = n_bio, confronto = come)
  }
  senza_piastra <- nel_confronto & biologico & is.na(piastra)
  if (any(senza_piastra)) {
    quali <- which(senza_piastra | (nel_confronto & negativo))
    confronti[["senza piastra"]] <- probabilita_prevalenza(conteggi[quali, , drop = FALSE],
                                                           negativo[quali])
    descrizione[["senza piastra"]] <- list(
      negativi = sum(negativo[quali]), biologici = sum(senza_piastra),
      confronto = "negativi di tutte le piastre: campioni senza piastra")
  }
  per_piastra_disponibile <- length(confronti) > 0L
  if (per_piastra_disponibile) {
    p_matrice <- do.call(cbind, confronti)
    p_per_piastra <- combina_probabilita(p_matrice, regola)
  } else {
    p_per_piastra <- stats::setNames(rep(NA_real_, length(varianti)), varianti)
  }

  contaminante <- function(p) !is.na(p) & p < soglia
  c_aggregata <- contaminante(p_aggregata)
  c_per_piastra <- contaminante(p_per_piastra)

  # ---- Misure e modalita' dichiarata ---------------------------------------
  frazione <- function(c, righe) {
    totale <- sum(conteggi[righe, ])
    if (totale == 0) NA_real_ else sum(conteggi[righe, c]) / totale
  }
  misura <- function(c) {
    list(contaminanti = sum(c), frazione_varianti = signif(mean(c), 6),
         letture_rimosse = list(
           biologico = signif(frazione(c, biologico), 6),
           controllo_negativo = signif(frazione(c, negativo), 6),
           controllo_positivo = signif(frazione(c, positivo), 6)))
  }
  m_aggregata <- misura(c_aggregata)
  m_per_piastra <- if (per_piastra_disponibile) misura(c_per_piastra) else NULL
  dichiarata <- if (modalita == "batch") m_per_piastra else m_aggregata
  if (is.null(dichiarata)) {
    stop("decontam.mode e' batch ma nessun campione ha una piastra")
  }
  rimossa_bio <- dichiarata$letture_rimosse$biologico
  entro <- !is.na(rimossa_bio) && rimossa_bio <= max_frazione
  esito <- sprintf(
    "%s: rimuove %.4f delle letture dei biologici, %s qc.max_frac_contaminant (%s)",
    if (modalita == "batch") "per piastra" else "aggregata", rimossa_bio,
    if (entro) "entro" else "oltre", format(max_frazione))
  rimuovere <- if (!entro) rep(FALSE, length(varianti)) else
    if (modalita == "batch") c_per_piastra else c_aggregata

  # ---- Artefatti per variante -----------------------------------------------
  prev <- function(righe) colSums(conteggi[righe, , drop = FALSE] > 0)
  lett <- function(righe) colSums(conteggi[righe, , drop = FALSE])
  genere <- if ("Genus" %in% colnames(tax)) tax[varianti, "Genus"] else rep(NA, length(varianti))
  numero <- function(x) ifelse(is.na(x), "", trimws(formatC(signif(x, 6), digits = 6, format = "g")))
  si_no <- function(x) ifelse(x, "si", "no")
  tabella <- sprintf(
    "%s\t%s\t%.0f\t%d\t%.0f\t%d\t%.0f\t%d\t%s\t%s\t%s\t%s\t%s",
    varianti, ifelse(is.na(genere), "", genere),
    lett(negativo), prev(negativo), lett(biologico), prev(biologico),
    lett(positivo), prev(positivo), numero(p_aggregata), numero(p_per_piastra),
    si_no(c_aggregata), si_no(c_per_piastra), si_no(rimuovere))
  intestazione <- paste(
    "asv_id", "genere", "letture_negativi", "prevalenza_negativi", "letture_biologici",
    "prevalenza_biologici", "letture_positivi", "prevalenza_positivi", "p_aggregata",
    "p_per_piastra", "contaminante_aggregata", "contaminante_per_piastra", "rimossa",
    sep = "\t")
  scrivi_atomico(c(intestazione, tabella), file.path(cartella, "decontam_varianti.tsv"))
  scrivi_atomico(c(intestazione, tabella[rimuovere]),
                 file.path(cartella, "decontam_contaminanti.tsv"))
  righe_piastra <- unlist(lapply(names(confronti), function(p) {
    v <- confronti[[p]]
    tenute <- !is.na(v)
    sprintf("%s\t%d\t%d\t%s\t%s\t%s", p, descrizione[[p]]$negativi, descrizione[[p]]$biologici,
            varianti[tenute], numero(v[tenute]), si_no(contaminante(v[tenute])))
  }))
  scrivi_atomico(c("piastra\tnegativi\tbiologici\tasv_id\tp\tcontaminante", righe_piastra),
                 file.path(cartella, "decontam_per_piastra.tsv"))

  # ---- L'oggetto ripulito ---------------------------------------------------
  tenute <- varianti[!rimuovere]
  pulito <- phyloseq::prune_taxa(tenute, ps)
  if (!identical(phyloseq::taxa_names(pulito), tenute) ||
      !identical(phyloseq::sample_names(pulito), phyloseq::sample_names(ps))) {
    stop("la rimozione dei contaminanti ha alterato campioni o identificativi")
  }
  salva_rds(pulito, file.path(cartella, "ps_decontaminato.rds"))
  tabella_pulita <- methods::as(phyloseq::otu_table(pulito), "matrix")
  rimaste <- if (phyloseq::taxa_are_rows(pulito)) colSums(tabella_pulita) else rowSums(tabella_pulita)
  rimaste <- rimaste[sort(names(rimaste), method = "radix")]

  # ---- Riepilogo ------------------------------------------------------------
  principali <- order(-lett(biologico) * rimuovere, varianti, method = "radix")[seq_len(min(10L, sum(rimuovere)))]
  scrivi_json(
    list(
      metodo = "prevalence",
      soglia = soglia,
      combinazione = regola,
      min_negativi = min_negativi,
      confronto = list(
        negativi = sum(nel_confronto & negativo),
        biologici = sum(nel_confronto & biologico),
        positivi_esclusi = sum(positivo),
        senza_letture_esclusi = sum((negativo | biologico) & !con_letture)
      ),
      piastre = descrizione,
      modalita = list(aggregata = m_aggregata, per_piastra = m_per_piastra),
      max_frazione = max_frazione,
      modalita_dichiarata = modalita,
      entro_max_frazione = entro,
      esito = esito,
      contaminanti_rimossi = sum(rimuovere),
      letture_rimosse = list(
        biologico = signif(frazione(rimuovere, biologico), 6),
        controllo_negativo = signif(frazione(rimuovere, negativo), 6),
        controllo_positivo = signif(frazione(rimuovere, positivo), 6)),
      principali_per_letture_biologiche = lapply(principali, function(i) list(
        asv_id = varianti[i], genere = if (is.na(genere[i])) NULL else genere[i],
        letture_biologici = lett(biologico)[[i]],
        p = signif(if (modalita == "batch") p_per_piastra[[i]] else p_aggregata[[i]], 6)))
    ),
    file.path(cartella, "decontam_riepilogo.json")
  )
  c("decontam_varianti.tsv", "decontam_per_piastra.tsv", "decontam_contaminanti.tsv",
    "decontam_riepilogo.json", "ps_decontaminato.rds",
    traccia_letture(rimaste, "decontaminate", cartella))
})
