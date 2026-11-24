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
                "costante": ["k"] * 10, "lotto": ["x"] * 9 + [None]}
    oggetto = costruisci(tmp_path, conteggi, metadati)
    rda = {"method": "rda", "variables": ["gruppo"], "transform": "hellinger",
           "max_df_fraction": 0.5, "conditions": ["costante"]}
    validazione = valida_sull_oggetto(oggetto, schema.valida(configurazione(constrained=rda)))
    avvisi = {a.codice: a.dettaglio for a in validazione.avvisi}
    assert "E-ECO-45" in avvisi and "'costante'" in avvisi["E-ECO-45"]
    uscita = tmp_path / "rda"
    esito = esegui(oggetto, schema.valida(configurazione(
        design={"other_variables": ["sito"]},
        constrained={**rda, "conditions": ["lotto"], "max_df_fraction": 0.1})), uscita)
    _, campioni = leggi_tsv(uscita / "campioni.tsv")
    assert {r["campione"]: r["vincolata"] for r in campioni}["c10"] == (
        "escluso: valore mancante in 'lotto'")
    assert "E-ECO-44" in {a.codice for a in esito.avvisi}
    for nome in ("varianza.tsv", "test.tsv", "rda_hellinger_coordinate.tsv",
                 "rda_hellinger_01_gruppo.png", "rda_hellinger_02_sito.png"):
        assert (uscita / "vincolata" / nome).is_file(), nome
    assert esito.riepilogo["vincolata"]["campioni"] == 9
    commenti, _ = leggi_tsv(uscita / "vincolata" / "varianza.tsv")
    assert any("modello: gruppo" in c for c in commenti)
