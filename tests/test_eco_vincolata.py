"""Ordinazione vincolata: RDA e db-RDA.

Si verifica che la configurazione chieda i parametri di ciascun metodo e
respinga quelli che non servono; che le quote di varianza, l'R2 corretto e il
test del modello coincidano con una chiamata a vegan scritta qui in modo
indipendente; che un effetto costruito sia trovato al netto della condizione;
che i campioni senza valore escano con il motivo; che il modello con troppi
gradi di liberta' sia segnalato; e che due esecuzioni diano gli stessi byte.
"""

from __future__ import annotations

import json

import pytest

from amplicon16s_eco import config as schema
from amplicon16s_eco.catalogo import CATALOGO, ErroreEco, Specie
from amplicon16s_eco.esecuzione import esegui, valida_sull_oggetto

from eco_aiuti import (  # noqa: F401
    configurazione,
    costruisci,
    eco_r,
    esegui_rscript,
    genera,
    impronte,
    leggi_tsv,
)

DBRDA = {"method": "dbrda", "variables": ["gruppo"], "distances": ["bray"],
         "max_df_fraction": 0.5}


def _codici(dati) -> list[str]:
    with pytest.raises(ErroreEco) as rifiuto:
        schema.valida(dati)
    return sorted(rifiuto.value.codici)


def test_ogni_metodo_chiede_i_suoi_parametri_e_respinge_gli_altri():
    """
    **Obiettivo**: senza metodo il gruppo ``constrained`` non entra nella
    configurazione risolta; con un metodo mancano variabili e frazione dei
    gradi di liberta' (``E-ECO-02``); la RDA chiede la trasformazione, la
    db-RDA le distanze, e ciascuna respinge il parametro dell'altra
    (``E-ECO-01``); con ``clr`` lo pseudoconteggio diventa obbligatorio anche
    senza Aitchison; condizioni senza metodo, variabili fuori dal disegno,
    distanze che ``beta.distances`` non calcola e variabili ripetute sono
    respinte con ``E-ECO-43``.

    **Razionale scientifico e sistemistico**: ogni scelta del metodo vincolato
    va dichiarata, e un parametro senza il suo metodo entrerebbe nel digest
    senza effetto.
    """
    assert "constrained" not in schema.risolta(schema.valida(configurazione()))
    assert _codici(configurazione(constrained={"method": "dbrda"})) == ["E-ECO-02"] * 3
    rda = {"method": "rda", "variables": ["gruppo"], "max_df_fraction": 0.5}
    assert _codici(configurazione(constrained=rda)) == ["E-ECO-02"]
    assert _codici(configurazione(constrained={**rda, "transform": "hellinger",
                                               "distances": ["bray"]})) == ["E-ECO-01"]
    assert _codici(configurazione(constrained={**DBRDA, "transform": "clr"})) == ["E-ECO-01"]
    solo_bray = {"distances": ["bray"], "clr_pseudocount": ...}
    assert _codici(configurazione(beta=solo_bray, constrained={**rda, "transform": "clr"})) == [
        "E-ECO-02"]
    schema.valida(configurazione(beta={**solo_bray, "clr_pseudocount": 1.0},
                                 constrained={**rda, "transform": "clr"}))
    for sbagliata in ({"conditions": ["lotto"]},
                      {**DBRDA, "variables": ["lotto"]},
                      {**DBRDA, "distances": ["unifrac_weighted"]},
                      {**DBRDA, "conditions": ["gruppo"]}):
        assert _codici(configurazione(constrained=sbagliata)) == ["E-ECO-43"], sbagliata
    assert _codici(configurazione(constrained={**DBRDA, "method": "cca"})) == ["E-ECO-03"]
    risolta = schema.risolta(schema.valida(configurazione(constrained=DBRDA)))
    assert risolta["constrained"] == DBRDA
    assert {CATALOGO[c].specie for c in ("E-ECO-44", "E-ECO-45")} == {Specie.AVVISO}


#: Due gruppi che differiscono nella composizione, incrociati con un lotto che
#: sposta i conteggi: la variabile si stima al netto del lotto.
SCENARIO = [
    {"n": 6, "gruppo": "a", "lotto": "x", "sigma": 0.2},
    {"n": 6, "gruppo": "a", "lotto": "y", "sigma": 0.2, "spostamento": 0.6},
    {"n": 6, "gruppo": "b", "lotto": "x", "sigma": 0.2, "scambia": True},
    {"n": 6, "gruppo": "b", "lotto": "y", "sigma": 0.2, "scambia": True, "spostamento": 0.6},
]


def test_la_db_rda_coincide_con_vegan_e_trova_l_effetto_al_netto_della_condizione(
    eco_r, tmp_path
):
    """
    **Obiettivo**: la db-RDA su Bray-Curtis con ``gruppo`` come variabile e
    ``lotto`` come condizione da' la stessa inerzia condizionata, vincolata e
    totale, lo stesso R2 corretto e la stessa statistica F di una chiamata a
    ``vegan::dbrda`` scritta qui; l'effetto costruito e' significativo (p
    minima con 99 permutazioni); coordinate, test e triplot sono scritti; due
    esecuzioni danno gli stessi byte.

    **Razionale scientifico e sistemistico**: la quota vincolata e il suo test
    sono i numeri che si riportano: devono essere quelli del metodo dichiarato,
    al netto della condizione, e non dipendere dall'esecuzione.
    """
    oggetto = genera(tmp_path, SCENARIO, profondita=500)
    config = schema.valida(configurazione(
        design={"technical_variables": ["lotto"]},
        constrained={**DBRDA, "conditions": ["lotto"]},
        stat={"permanova_permutations": 99}))
    esegui(oggetto, config, tmp_path / "a")
    esegui(oggetto, config, tmp_path / "b")
    assert impronte(tmp_path / "a") == impronte(tmp_path / "b")
    cartella = tmp_path / "a" / "vincolata"
    assert (cartella / "dbrda_bray_01_gruppo.png").stat().st_size > 1000
    _, varianza = leggi_tsv(cartella / "varianza.tsv")
    riga = varianza[0]
    assert riga["modello"] == "dbrda_bray" and riga["campioni"] == "24"
    _, test = leggi_tsv(cartella / "test.tsv")
    per = {r["termine"]: r for r in test}
    assert float(per["modello"]["p"]) == pytest.approx(0.01)
    assert per["gruppo"]["ruolo"] == "variabile" and per["gruppo"]["gl"] == "1"
    _, coordinate = leggi_tsv(cartella / "dbrda_bray_coordinate.tsv")
    assert len(coordinate) == 24 and list(coordinate[0]) == ["campione", "dbRDA1", "MDS1"]

    script = tmp_path / "vegan.R"
    atteso = tmp_path / "atteso.json"
    script.write_text(f"""
ps <- readRDS("{oggetto}")
m <- methods::as(phyloseq::otu_table(ps), "matrix")
p <- t(sweep(m, 2, colSums(m), "/"))
p <- p[order(rownames(p), method = "radix"), ]
d <- vegan::vegdist(p, method = "bray")
dati <- methods::as(phyloseq::sample_data(ps), "data.frame")[rownames(p), ]
modello <- vegan::dbrda(d ~ gruppo + Condition(lotto), data = dati)
set.seed(7)
f <- vegan::anova.cca(modello, permutations = 99)$F[1]
jsonlite::write_json(list(totale = modello$tot.chi, condizionata = modello$pCCA$tot.chi,
  vincolata = modello$CCA$tot.chi, r2c = vegan::RsquareAdj(modello)$adj.r.squared, f = f),
  "{atteso}", digits = NA, auto_unbox = TRUE)
""", encoding="utf-8")
    esegui_rscript(script)
    vegan = json.loads(atteso.read_text(encoding="utf-8"))
    assert float(riga["inerzia_totale"]) == pytest.approx(vegan["totale"], abs=1e-9)
    assert float(riga["condizionata"]) == pytest.approx(vegan["condizionata"], abs=1e-9)
    assert float(riga["vincolata"]) == pytest.approx(vegan["vincolata"], abs=1e-9)
    assert float(riga["r2_corretto"]) == pytest.approx(vegan["r2c"], abs=1e-9)
    assert float(per["modello"]["f"]) == pytest.approx(vegan["f"], abs=1e-9)
    assert float(riga["quota_vincolata"]) == pytest.approx(
        vegan["vincolata"] / vegan["totale"], abs=1e-9)


def test_la_rda_di_hellinger_i_campioni_esclusi_e_gli_avvisi_del_modello(eco_r, tmp_path):
    """
    **Obiettivo**: la RDA sui conteggi trasformati con Hellinger scrive le sue
    uscite; un campione senza valore in una variabile del modello esce con il
    motivo in ``campioni.tsv``; un modello che usa piu' gradi di liberta' di
    ``max_df_fraction`` dei campioni da' ``E-ECO-44`` gia' nella validazione;
    una condizione con un solo valore fra i campioni esce dal modello con
    ``E-ECO-45``.

    **Razionale scientifico e sistemistico**: con molti livelli il modello si
    avvicina a quello saturo e la quota vincolata non dice piu' nulla; un
    campione senza valore non puo' entrare in un modello che lo richiede.
    """
    conteggi = [[20 + (3 * i + 5 * j) % 11 + (25 if (j % 2 and i < 2) else 0) for j in range(10)]
                for i in range(5)]
    metadati = {"gruppo": ["a", "b"] * 5, "sito": [f"s{j // 2}" for j in range(10)],
                "costante": ["k"] * 10, "lotto": ["x"] * 9 + [None],
                # Un valore solo nel campione che esce dal modello: nessun triplot.
                "esterno": [None] * 9 + ["e"]}
    oggetto = costruisci(tmp_path, conteggi, metadati)
    rda = {"method": "rda", "variables": ["gruppo"], "transform": "hellinger",
           "max_df_fraction": 0.5, "conditions": ["costante"]}
    validazione = valida_sull_oggetto(oggetto, schema.valida(configurazione(constrained=rda)))
    avvisi = {a.codice: a.dettaglio for a in validazione.avvisi}
    assert "E-ECO-45" in avvisi and "'costante'" in avvisi["E-ECO-45"]
    uscita = tmp_path / "rda"
    esito = esegui(oggetto, schema.valida(configurazione(
        design={"other_variables": ["sito", "esterno"]},
        constrained={**rda, "conditions": ["lotto"], "max_df_fraction": 0.1})), uscita)
    _, campioni = leggi_tsv(uscita / "campioni.tsv")
    assert {r["campione"]: r["vincolata"] for r in campioni}["c10"] == (
        "escluso: valore mancante in 'lotto'")
    assert "E-ECO-44" in {a.codice for a in esito.avvisi}
    for nome in ("varianza.tsv", "test.tsv", "rda_hellinger_coordinate.tsv",
                 "rda_hellinger_01_gruppo.png", "rda_hellinger_02_sito.png"):
        assert (uscita / "vincolata" / nome).is_file(), nome
    assert esito.riepilogo["vincolata"]["campioni"] == 9
    assert not (uscita / "vincolata" / "rda_hellinger_03_esterno.png").exists()
    commenti, _ = leggi_tsv(uscita / "vincolata" / "varianza.tsv")
    assert any("modello: gruppo" in c for c in commenti)


def test_le_due_unifrac_entrano_in_permanova_ordinazioni_tsne_e_db_rda(eco_r, tmp_path):
    """
    **Obiettivo**: con un albero nell'oggetto, UniFrac pesato e non pesato
    dichiarati in ``beta.distances`` producono la loro matrice, la PERMANOVA,
    la PCoA, il t-SNE e un modello di db-RDA ciascuno; ``eco init`` dice in una
    riga che cosa misura ogni distanza.

    **Razionale scientifico e sistemistico**: le UniFrac sono distanze come le
    altre per ogni analisi che parte da una matrice; chi le sceglie deve sapere
    che misurano la parentela filogenetica, con e senza abbondanze.
    """
    from amplicon16s_eco import init

    # Il gruppo b ha meno varianti, e con poche letture la presenza cambia da
    # un campione all'altro: UniFrac non pesato ha piu' di un asse.
    oggetto = genera(tmp_path, [
        {"n": 6, "gruppo": "a", "sigma": 1.0},
        {"n": 6, "gruppo": "b", "sigma": 1.0, "scambia": True, "presenti": 6},
    ], profondita=20, taxa=8)
    # Un albero bilanciato sulle otto varianti del generatore.
    foglie = [f"v{i:03d}" for i in range(1, 9)]
    albero = "(((" + ",".join(f"{f}:1" for f in foglie[:2]) + "):1,(" + ",".join(
        f"{f}:1" for f in foglie[2:4]) + "):1):1,((" + ",".join(
        f"{f}:1" for f in foglie[4:6]) + "):1,(" + ",".join(
        f"{f}:1" for f in foglie[6:]) + "):1):1);"
    script = tmp_path / "albero.R"
    con_albero = tmp_path / "dati" / "con_albero.rds"
    script.write_text(f"""
ps <- readRDS("{oggetto}")
phyloseq::phy_tree(ps) <- phyloseq::phy_tree(ape::read.tree(text = "{albero}"))
saveRDS(ps, "{con_albero}")
""", encoding="utf-8")
    esegui_rscript(script)
    unifrac = ["unifrac_weighted", "unifrac_unweighted"]
    uscita = tmp_path / "uscita"
    esegui(con_albero, schema.valida(configurazione(
        beta={"distances": ["bray", *unifrac], "clr_pseudocount": ..., "unifrac_root": "existing"},
        ord={"methods": ["pcoa", "tsne"], "distances": unifrac, "pcoa_correction": "none",
             "tsne_perplexity": 3.0},
        constrained={**DBRDA, "distances": unifrac},
        stat={"permanova_permutations": 99})), uscita)
    _, permanova = leggi_tsv(uscita / "test" / "permanova.tsv")
    assert {r["distanza"] for r in permanova} == {"bray", *unifrac}
    for distanza in unifrac:
        assert (uscita / "beta" / f"distanza_{distanza}.tsv").is_file()
        for metodo in ("pcoa", "tsne"):
            assert (uscita / "ordinazione" / f"{metodo}_{distanza}_coordinate.tsv").is_file()
        assert (uscita / "vincolata" / f"dbrda_{distanza}_coordinate.tsv").is_file()
    _, varianza = leggi_tsv(uscita / "vincolata" / "varianza.tsv")
    assert [r["modello"] for r in varianza] == [f"dbrda_{d}" for d in unifrac]

    descrizione = {"campioni": 24, "varianti": 8, "ranghi": [{"rango": "Regno",
                                                              "quota_assegnata": 1.0}],
                   "albero": {"presente": True, "radicato": True, "con_lunghezze": True},
                   "colonne": []}
    voci, _ = init.prepara(descrizione)
    note = " ".join(voci["beta.distances"].note)
    assert "parentela filogenetica dei taxa, con le loro abbondanze" in note
    assert "presenti, senza abbondanze" in note and "presenza o assenza" in note
    assert "log-rapporti" in note and "abbondanza relativa" in note


def test_le_etichette_dei_centroidi_non_si_sovrappongono_e_i_centroidi_hanno_il_colore_del_gruppo(
        eco_r, tmp_path):
    """
    **Obiettivo**: nel triplot le etichette di centroidi vicini (qui otto
    centroidi, cinque in un gruppo stretto come NOD1, NOD3, NOD2, JPM1 e COL1
    nella dimostrazione su OSD-734) non si sovrappongono fra loro e non
    coprono alcun centroide; ognuna resta accostata al suo centroide quando
    c'e' posto; le posizioni sono le stesse a ogni esecuzione. Un centroide e'
    pieno del colore del suo gruppo se il termine e' la variabile del grafico,
    bianco altrimenti. Il triplot si disegna senza avvisi.

    **Razionale scientifico e sistemistico**: un'etichetta sopra un'altra non
    si legge, e un'etichetta lontana dal suo centroide si attribuisce al
    centroide sbagliato; il colore lega il centroide ai campioni del gruppo.
    ``ggrepel`` non e' nell'immagine: lo spostamento e' calcolato qui.
    """
    from eco_aiuti import RADICE

    lib = RADICE / "src" / "amplicon16s_eco" / "R" / "lib"
    script = tmp_path / "etichette.R"
    script.write_text("\n".join([
        f'source("{RADICE / "R" / "lib" / "io_json.R"}")',
        f'source("{lib / "uscite.R"}")',
        f'source("{lib / "vincolata.R"}")',
        "x <- c(-0.55, -0.65, 0.18, 0.37, 0.67, 1.03, 1.33, -0.78)",
        "y <- c(0.64, 0.43, -0.17, 0.03, -0.23, -0.30, 0.70, -1.77)",
        "testi <- c('NOD1', 'NOD3', 'NOD2', 'LAB1', 'JPM1', 'COL1', 'JLP1', 'PMM1')",
        f'grDevices::png("{tmp_path / "prova.png"}", width = 2400, height = 1500, res = 200, type = "cairo")',
        "graphics::par(mar = c(4.5, 4.5, 3, 15), cex = 0.7)",
        "graphics::plot(c(-2.7, 3.5), c(-2.6, 6.3), type = 'n', asp = 1)",
        "posa <- .posa_etichette(x, y, testi, CEX_ETICHETTA, CEX_CENTROIDE)",
        "stopifnot(identical(posa, .posa_etichette(x, y, testi, CEX_ETICHETTA, CEX_CENTROIDE)))",
        "raggio <- graphics::strheight('M', cex = CEX_CENTROIDE) * 0.55",
        "for (i in seq_along(x)) {",
        "  for (j in seq_along(x)) if (i < j) {",
        "    stopifnot(abs(posa$x[i] - posa$x[j]) >= (posa$larghezza[i] + posa$larghezza[j]) / 2 ||",
        "              abs(posa$y[i] - posa$y[j]) >= (posa$altezza[i] + posa$altezza[j]) / 2)",
        "  }",
        "  stopifnot(all(abs(posa$x[i] - x) >= posa$larghezza[i] / 2 + raggio |",
        "                abs(posa$y[i] - y) >= posa$altezza[i] / 2 + raggio))",
        "}",
        "stopifnot(!any(posa$spostata))",
        "grDevices::dev.off()",
        "piano <- list(nomi = c(v1 = 'modulo', v2 = 'lotto'),",
        "              dati = data.frame(v1 = factor(c('A', 'B')), v2 = factor(c('x', 'y'))))",
        "colori <- c(A = '#111111', B = '#222222')",
        "stopifnot(identical(.colori_centroidi(c('v1A', 'v1B', 'v2x'), piano,",
        "                                      list(nome = 'modulo'), colori),",
        "                    c('#111111', '#222222', 'white')))",
        "imposta_palette('Okabe-Ito')",
        "set.seed(1)",
        "coordinate <- cbind(rnorm(40), rnorm(40))",
        "gruppi <- rep(c('A', 'B'), each = 20)",
        "rownames(coordinate) <- sprintf('c%02d', 1:40)",
        "centroidi <- rbind(v1A = colMeans(coordinate[1:20, ]), v1B = colMeans(coordinate[21:40, ]))",
        "v <- list(nome = 'modulo', valori = stats::setNames(gruppi, rownames(coordinate)),",
        "          dimensioni = c(A = 20L, B = 20L))",
        "esito <- list(coordinate = coordinate, centroidi = centroidi, frecce = NULL,",
        "              assi = c('dbRDA1', 'dbRDA2'))",
        "avvisi <- NULL",
        "withCallingHandlers(",
        f"  .triplot(esito, list(nomi = c(v1 = 'modulo'), dati = data.frame(v1 = factor(gruppi))), v, 'prova', '{tmp_path / 'triplot.png'}'),",
        "  warning = function(w) { avvisi <<- c(avvisi, conditionMessage(w)); invokeRestart('muffleWarning') })",
        "stopifnot(is.null(avvisi))",
    ]), encoding="utf-8")
    esegui_rscript(script)
    assert (tmp_path / "triplot.png").stat().st_size > 1000
