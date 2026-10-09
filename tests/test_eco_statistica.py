"""Test fra gruppi delle analisi ecologiche: correttezza con un segnale noto.

I dati sono sintetici e costruiti in R con una struttura decisa in partenza
(differenza di posizione, di dispersione, nessuna differenza, effetto del solo
lotto, differenza di ricchezza): solo cosi' si verifica che un test statistico
dia il risultato giusto, e non soltanto un risultato. Verifica inoltre le
dichiarazioni che accompagnano i test: gruppi esclusi, confusione con il lotto,
termini non stimabili, permutazioni insufficienti.
"""

from __future__ import annotations

import hashlib

import pytest

from amplicon16s_eco import config as schema
from amplicon16s_eco.catalogo import ErroreEco
from amplicon16s_eco.esecuzione import esegui, valida_sull_oggetto

from eco_aiuti import (  # noqa: F401
    configurazione,
    costruisci,
    eco_r,
    genera,
    impronte,
    leggi_tsv,
    righe_permanova,
)

#: Una sola distanza e nessuna ordinazione: i test statistici non ne dipendono.
SOLO_BRAY = {"distances": ["bray"], "clr_pseudocount": ...}


def _config(**modifica):
    base = {"beta": SOLO_BRAY, "stat": {"min_group_size": 5}}
    for gruppo, valori in modifica.items():
        base[gruppo] = {**base.get(gruppo, {}), **valori}
    return schema.valida(configurazione(**base))


def _codici(esito) -> list[str]:
    return [a.codice for a in esito.avvisi]


def test_due_gruppi_che_differiscono_nella_posizione(eco_r, tmp_path):
    """
    **Obiettivo**: due gruppi di venti campioni con composizione media diversa
    e la stessa dispersione danno una PERMANOVA significativa sulla variabile
    (p < 0,01 con 999 permutazioni) e dispersioni non significative, senza
    l'avviso ``E-ECO-19``.

    **Razionale scientifico e sistemistico**: e' il caso in cui la PERMANOVA si
    legge come differenza di posizione; il test delle dispersioni deve
    confermarlo e non segnalare una differenza che non c'e'.
    """
    oggetto = genera(tmp_path, [
        {"n": 20, "gruppo": "a", "sigma": 0.3},
        {"n": 20, "gruppo": "b", "sigma": 0.3, "scambia": True},
    ])
    esito = esegui(oggetto, _config(), tmp_path / "uscita")
    righe = righe_permanova(tmp_path / "uscita")
    variabile = righe[("permanova", "gruppo")]
    assert variabile["ruolo"] == "analizzata" and variabile["gl"] == "1"
    assert float(variabile["p"]) < 0.01
    assert float(variabile["r2"]) > 0.5 and float(variabile["f"]) > 10
    dispersioni = righe[("dispersioni", "gruppo")]
    assert float(dispersioni["p"]) > 0.05 and dispersioni["nota"] == ""
    assert "E-ECO-19" not in _codici(esito)
    assert variabile["permutazioni"] == "999" and variabile["campioni"] == "40"
    totale, residuo = righe[("permanova", "totale")], righe[("permanova", "residuo")]
    assert float(totale["r2"]) == pytest.approx(1.0, abs=1e-9)
    assert float(variabile["r2"]) + float(residuo["r2"]) == pytest.approx(1.0, abs=1e-8)


def test_stessa_composizione_media_e_dispersione_diversa(eco_r, tmp_path):
    """
    **Obiettivo**: due gruppi con la stessa composizione media e dispersione
    diversa danno dispersioni significative; la riga delle dispersioni e
    l'avviso ``E-ECO-19`` dichiarano che la PERMANOVA non va letta come
    differenza di posizione.

    **Razionale scientifico e sistemistico**: la PERMANOVA e' sensibile anche
    alla diversa dispersione; senza il test che la accompagna un risultato
    significativo verrebbe letto come differenza di composizione.
    """
    oggetto = genera(tmp_path, [
        {"n": 20, "gruppo": "a", "sigma": 0.1},
        {"n": 20, "gruppo": "b", "sigma": 0.9},
    ])
    esito = esegui(oggetto, _config(), tmp_path / "uscita")
    dispersioni = righe_permanova(tmp_path / "uscita")[("dispersioni", "gruppo")]
    assert float(dispersioni["p"]) < 0.01
    assert "la PERMANOVA non va letta come differenza di posizione" in dispersioni["nota"]
    avvisi = [a for a in esito.avvisi if a.codice == "E-ECO-19"]
    assert len(avvisi) == 1 and "'gruppo', distanza bray" in avvisi[0].dettaglio


def test_senza_differenza_i_falsi_positivi_restano_al_livello_atteso(eco_r, tmp_path):
    """
    **Obiettivo**: su venti dataset generati con semi diversi, in cui i due
    gruppi hanno la stessa distribuzione, la PERMANOVA da' p sotto 0,05 al
    piu' tre volte, e lo stesso vale per le dispersioni.

    **Razionale scientifico e sistemistico**: un test che rileva sempre il
    segnale costruito potrebbe rilevarlo anche dove non c'e'; la frequenza dei
    falsi positivi sotto l'ipotesi nulla deve restare vicina al livello
    nominale (uno atteso su venti).
    """
    config = _config()
    p_permanova, p_dispersioni = [], []
    for seme in range(1, 21):
        cartella = tmp_path / f"seme_{seme:02d}"
        cartella.mkdir()
        oggetto = genera(cartella, [
            {"n": 10, "gruppo": "a", "sigma": 0.4},
            {"n": 10, "gruppo": "b", "sigma": 0.4},
        ], seme=seme, taxa=30)
        esegui(oggetto, config, cartella / "uscita")
        righe = righe_permanova(cartella / "uscita")
        p_permanova.append(float(righe[("permanova", "gruppo")]["p"]))
        p_dispersioni.append(float(righe[("dispersioni", "gruppo")]["p"]))
    assert sum(p < 0.05 for p in p_permanova) <= 3, p_permanova
    assert sum(p < 0.05 for p in p_dispersioni) <= 3, p_dispersioni
    # Le p non sono tutte uguali ne' tutte al minimo: il test risponde ai dati.
    assert len(set(p_permanova)) > 10 and max(p_permanova) > 0.5


#: La differenza e' tutta del lotto; i gruppi sono sbilanciati fra i lotti.
EFFETTO_LOTTO = [
    {"n": 16, "gruppo": "a", "lotto": "x", "sigma": 0.3, "spostamento": 0.8},
    {"n": 4, "gruppo": "b", "lotto": "x", "sigma": 0.3, "spostamento": 0.8},
    {"n": 4, "gruppo": "a", "lotto": "y", "sigma": 0.3, "spostamento": -0.8},
    {"n": 16, "gruppo": "b", "lotto": "y", "sigma": 0.3, "spostamento": -0.8},
]


def test_un_effetto_del_solo_lotto_non_viene_attribuito_alla_variabile(eco_r, tmp_path):
    """
    **Obiettivo**: quando la differenza e' tutta della variabile tecnica e i
    gruppi sono sbilanciati fra i lotti, con il lotto come primo termine l'R²
    della variabile resta vicino a zero e quello del lotto e' riportato accanto;
    senza la variabile tecnica la stessa variabile risulta significativa con
    un R² molto piu' alto, e l'analisi dichiara ``E-ECO-20``.

    **Razionale scientifico e sistemistico**: le somme dei quadrati sequenziali
    stimano l'effetto della variabile al netto del lotto; senza il termine
    tecnico l'effetto del lotto verrebbe letto come effetto biologico.
    """
    oggetto = genera(tmp_path, EFFETTO_LOTTO)
    con = esegui(oggetto, _config(design={"technical_variables": ["lotto"]}), tmp_path / "con")
    righe = righe_permanova(tmp_path / "con")
    assert list(righe)[:2] == [("permanova", "lotto"), ("permanova", "gruppo")]
    lotto, variabile = righe[("permanova", "lotto")], righe[("permanova", "gruppo")]
    assert lotto["ruolo"] == "tecnica" and variabile["ruolo"] == "analizzata"
    assert float(lotto["r2"]) > 0.4 and float(lotto["p"]) < 0.01
    assert float(variabile["r2"]) < 0.03
    assert "E-ECO-20" not in _codici(con)

    senza = esegui(oggetto, _config(), tmp_path / "senza")
    sola = righe_permanova(tmp_path / "senza")[("permanova", "gruppo")]
    assert float(sola["p"]) < 0.01 and float(sola["r2"]) > 5 * float(variabile["r2"])
    assert "E-ECO-20" in _codici(senza)


def _bh(p: list[float]) -> list[float]:
    """Benjamini-Hochberg calcolato qui, in modo indipendente da R."""
    n = len(p)
    ordine = sorted(range(n), key=lambda i: p[i], reverse=True)
    corretti, minimo = [0.0] * n, 1.0
    for posizione, i in enumerate(ordine):
        rango = n - posizione
        minimo = min(minimo, p[i] * n / rango)
        corretti[i] = minimo
    return corretti


def test_una_differenza_di_ricchezza_e_rilevata_e_corretta_sulla_famiglia(eco_r, tmp_path):
    """
    **Obiettivo**: con tre gruppi, di cui uno con meta' dei taxa presenti,
    Kruskal-Wallis rileva la differenza di ricchezza; i confronti a coppie con
    Wilcoxon la attribuiscono alle coppie giuste; la p corretta e' quella di
    Benjamini-Hochberg sulla famiglia di tutti i test della variabile (tre
    indici per quattro confronti), ricalcolata qui.

    **Razionale scientifico e sistemistico**: la correzione dipende da quali
    test formano la famiglia; dichiararla e ricalcolarla prova che e' quella
    scritta nella tabella.
    """
    oggetto = genera(tmp_path, [
        {"n": 12, "gruppo": "a", "sigma": 0.3},
        {"n": 12, "gruppo": "b", "sigma": 0.3, "presenti": 20},
        {"n": 12, "gruppo": "c", "sigma": 0.3},
    ], profondita=300)
    esegui(oggetto, _config(), tmp_path / "uscita")
    commenti, righe = leggi_tsv(tmp_path / "uscita" / "test" / "alfa_diversita.tsv")
    assert len(righe) == 12 and {r["dimensione_famiglia"] for r in righe} == {"12"}
    assert any("famiglia di tutti i test della variabile" in c for c in commenti)
    assert [r["confronto"] for r in righe[:4]] == ["globale", "a vs b", "a vs c", "b vs c"]
    ricchezza = {r["confronto"]: r for r in righe if r["indice"] == "ricchezza_osservata"}
    assert ricchezza["globale"]["test"] == "Kruskal-Wallis" and ricchezza["globale"]["gl"] == "2"
    assert float(ricchezza["globale"]["p_bh"]) < 0.001
    assert float(ricchezza["a vs b"]["p_bh"]) < 0.001
    assert float(ricchezza["b vs c"]["p_bh"]) < 0.001
    assert float(ricchezza["a vs c"]["p"]) > 0.05
    assert ricchezza["a vs b"]["n_1"] == "12" and ricchezza["a vs b"]["n_2"] == "12"
    assert ricchezza["a vs b"]["test"].startswith("Wilcoxon della somma dei ranghi")

    p = [float(r["p"]) for r in righe]
    for riga, atteso in zip(righe, _bh(p), strict=True):
        assert float(riga["p_bh"]) == pytest.approx(atteso, rel=1e-5)
    assert all(float(r["p_bh"]) >= float(r["p"]) * (1 - 1e-6) for r in righe)


def test_con_due_gruppi_un_wilcoxon_per_indice_esatto_senza_valori_ripetuti(eco_r, tmp_path):
    """
    **Obiettivo**: con due gruppi la famiglia ha tre test, uno per indice; il
    Wilcoxon e' esatto dove i valori non si ripetono (Shannon) e con
    approssimazione normale dove si ripetono, e la tabella lo dice; con cinque
    campioni contro cinque completamente separati la p esatta vale 2/252.

    **Razionale scientifico e sistemistico**: la forma del test cambia la p con
    pochi campioni; la scelta segue una regola dichiarata e resta scritta
    accanto a ogni risultato.
    """
    oggetto = genera(tmp_path, [
        {"n": 5, "gruppo": "a", "sigma": 0.05},
        {"n": 5, "gruppo": "b", "sigma": 0.05, "presenti": 20},
    ])
    esegui(oggetto, _config(), tmp_path / "uscita")
    _, righe = leggi_tsv(tmp_path / "uscita" / "test" / "alfa_diversita.tsv")
    assert [r["indice"] for r in righe] == ["ricchezza_osservata", "shannon", "gini_simpson"]
    assert {r["dimensione_famiglia"] for r in righe} == {"3"}
    shannon = righe[1]
    assert shannon["test"] == "Wilcoxon della somma dei ranghi, esatto"
    assert float(shannon["p"]) == pytest.approx(2 / 252, rel=1e-5)
    assert "approssimazione normale" in righe[0]["test"]
    for riga, atteso in zip(righe, _bh([float(r["p"]) for r in righe]), strict=True):
        assert float(riga["p_bh"]) == pytest.approx(atteso, rel=1e-5)


def test_un_test_non_calcolabile_ha_p_assente_e_non_entra_nella_famiglia(eco_r, tmp_path):
    """
    **Obiettivo**: quando un indice ha lo stesso valore in tutti i campioni
    (ogni variante presente ovunque: ricchezza costante) il suo test ha p
    ``NA``, non entra nella famiglia, e gli altri sono corretti sulla famiglia
    dei soli test calcolati.

    **Razionale scientifico e sistemistico**: un test senza variabilita' non ha
    una p; contarlo nella famiglia renderebbe la correzione piu' severa senza
    che un test sia stato fatto.
    """
    conteggi = [[20 + (3 * i + 5 * j) % 11 + (30 if j >= 6 and i < 2 else 0) for j in range(12)]
                for i in range(4)]
    oggetto = costruisci(tmp_path, conteggi, {"gruppo": ["a"] * 6 + ["b"] * 6})
    esegui(oggetto, _config(), tmp_path / "uscita")
    commenti, righe = leggi_tsv(tmp_path / "uscita" / "test" / "alfa_diversita.tsv")
    assert righe[0]["indice"] == "ricchezza_osservata"
    assert righe[0]["p"] == "NA" and righe[0]["p_bh"] == "NA"
    assert {r["dimensione_famiglia"] for r in righe} == {"2"}
    calcolati = righe[1:]
    for riga, atteso in zip(calcolati, _bh([float(r["p"]) for r in calcolati]), strict=True):
        assert float(riga["p_bh"]) == pytest.approx(atteso, rel=1e-5)
    assert any("non calcolabile" in c for c in commenti)


def test_i_gruppi_piccoli_escono_dai_test_e_sotto_due_gruppi_il_test_non_si_esegue(eco_r, tmp_path):
    """
    **Obiettivo**: un terzo gruppo con due campioni esce dai test con
    ``E-ECO-11`` (la PERMANOVA riporta due gruppi e i campioni dei soli due) e
    resta nelle uscite descrittive; con un solo gruppo sopra il minimo nessun
    test si esegue, con ``E-ECO-12``, e le tabelle dei test hanno la sola
    intestazione.

    **Razionale scientifico e sistemistico**: un gruppo con pochi campioni non
    sostiene un test e non deve entrare nei gradi di liberta'; l'analisi
    prosegue e lo dichiara.
    """
    blocchi = [
        {"n": 8, "gruppo": "a", "sigma": 0.3},
        {"n": 8, "gruppo": "b", "sigma": 0.3, "scambia": True},
        {"n": 2, "gruppo": "c", "sigma": 0.3},
    ]
    oggetto = genera(tmp_path, blocchi)
    esito = esegui(oggetto, _config(), tmp_path / "tre")
    assert _codici(esito).count("E-ECO-11") == 1 and "E-ECO-12" not in _codici(esito)
    variabile = righe_permanova(tmp_path / "tre")[("permanova", "gruppo")]
    assert variabile["gruppi"] == "2" and variabile["campioni"] == "16" and variabile["gl"] == "1"
    _, alfa = leggi_tsv(tmp_path / "tre" / "test" / "alfa_diversita.tsv")
    assert {r["confronto"] for r in alfa} == {"a vs b"}
    _, campioni = leggi_tsv(tmp_path / "tre" / "campioni.tsv")
    stati = {r["permanova:gruppo"] for r in campioni if r["valore:gruppo"] == "c"}
    assert stati == {"escluso dai test: gruppo sotto la dimensione minima"}
    _, gruppi = leggi_tsv(tmp_path / "tre" / "composizione" / "gruppi_01_gruppo.tsv")
    assert "c" in gruppi[0]

    esito = esegui(oggetto, _config(stat={"min_group_size": 9}), tmp_path / "nessuno")
    assert "E-ECO-12" in _codici(esito)
    for nome in ("alfa_diversita.tsv", "permanova.tsv"):
        _, righe = leggi_tsv(tmp_path / "nessuno" / "test" / nome)
        assert righe == []
    (prova,) = esito.riepilogo["test"]
    assert prova["alfa"]["eseguito"] is False and prova["permanova"]["eseguita"] is False


#: Dodici campioni, quattro varianti: i test sul disegno non dipendono dai conteggi.
CONTEGGI_12 = [[5 + (i * j) % 7 for j in range(12)] for i in range(1, 5)]


def _avvisi_del_disegno(tmp_path, metadati, **modifica) -> dict[str, list[str]]:
    oggetto = costruisci(tmp_path, CONTEGGI_12, metadati)
    stat = {"min_group_size": 2, **modifica.pop("stat", {})}
    esito = valida_sull_oggetto(oggetto, _config(stat=stat, **modifica))
    per_codice: dict[str, list[str]] = {}
    for avviso in esito.avvisi:
        per_codice.setdefault(avviso.codice, []).append(avviso.dettaglio)
    return per_codice


def test_la_confusione_fra_variabile_e_lotto_e_dichiarata_prima_del_calcolo(eco_r, tmp_path):
    """
    **Obiettivo**: la sola validazione dichiara ``E-ECO-18`` quando un gruppo
    sta per intero in un lotto, ``E-ECO-23`` quando la variabile e' determinata
    per intero dal lotto, ``E-ECO-22`` quando il lotto ha un solo livello e
    ``E-ECO-20`` quando non e' dichiarata alcuna variabile tecnica; in un
    disegno incrociato non dichiara nulla.

    **Razionale scientifico e sistemistico**: se un gruppo biologico sta in una
    sola piastra, le somme sequenziali attribuiscono quella varianza al lotto;
    chi legge il risultato deve saperlo, e lo si puo' dire dal solo disegno.
    """
    gruppo = ["a"] * 4 + ["b"] * 4 + ["c"] * 4
    tecnica = {"design": {"technical_variables": ["lotto"]}}
    incrociato = _avvisi_del_disegno(tmp_path, {"gruppo": gruppo, "lotto": ["x", "y"] * 6}, **tecnica)
    assert incrociato == {}

    parziale = _avvisi_del_disegno(
        tmp_path, {"gruppo": gruppo, "lotto": ["x"] * 4 + ["x", "y"] * 4}, **tecnica)
    assert list(parziale) == ["E-ECO-18"]
    assert "1 gruppi su 3" in parziale["E-ECO-18"][0] and "(a)" in parziale["E-ECO-18"][0]

    annidato = _avvisi_del_disegno(
        tmp_path, {"gruppo": gruppo, "lotto": ["x"] * 4 + ["y"] * 4 + ["z"] * 4}, **tecnica)
    assert sorted(annidato) == ["E-ECO-18", "E-ECO-23"]

    costante = _avvisi_del_disegno(tmp_path, {"gruppo": gruppo, "lotto": ["x"] * 12}, **tecnica)
    assert list(costante) == ["E-ECO-22"]

    nessuna = _avvisi_del_disegno(tmp_path, {"gruppo": gruppo, "lotto": ["x", "y"] * 6})
    assert list(nessuna) == ["E-ECO-20"]


def test_una_variabile_determinata_dal_lotto_non_ha_un_termine_stimabile(eco_r, tmp_path):
    """
    **Obiettivo**: quando ogni gruppo coincide con un lotto la tabella riporta
    il termine tecnico, il termine della variabile senza gradi di liberta' con
    la nota «non stimabile», e comunque il test delle dispersioni.

    **Razionale scientifico e sistemistico**: ``adonis2`` toglie in silenzio un
    termine senza gradi di liberta'; la tabella deve dire che il termine c'era
    e non era stimabile, non ometterlo.
    """
    oggetto = genera(tmp_path, [
        {"n": 8, "gruppo": "a", "lotto": "x", "sigma": 0.3},
        {"n": 8, "gruppo": "b", "lotto": "y", "sigma": 0.3, "scambia": True},
    ])
    esito = esegui(oggetto, _config(design={"technical_variables": ["lotto"]}), tmp_path / "uscita")
    assert {"E-ECO-18", "E-ECO-23"} <= set(_codici(esito))
    righe = righe_permanova(tmp_path / "uscita")
    assert float(righe[("permanova", "lotto")]["p"]) < 0.01
    variabile = righe[("permanova", "gruppo")]
    assert variabile["gl"] == "0" and variabile["p"] == "NA"
    assert variabile["nota"].startswith("non stimabile")
    assert ("dispersioni", "gruppo") in righe
    (prova,) = esito.riepilogo["test"]
    assert prova["permanova"]["termine_della_variabile_stimabile"] is False


def test_gli_strati_restringono_le_permutazioni_e_sono_dichiarati(eco_r, tmp_path):
    """
    **Obiettivo**: con ``stat.permanova_strata`` la PERMANOVA permuta entro gli
    strati: l'intestazione li nomina, la statistica F e l'R² restano quelli del
    modello senza strati mentre la p cambia; uno strato che contiene per intero
    un gruppo e' dichiarato con ``E-ECO-18``; una colonna inesistente e'
    respinta con ``E-ECO-05``.

    **Razionale scientifico e sistemistico**: gli strati cambiano la
    distribuzione di riferimento, non la statistica; con l'effetto del lotto
    presente la permutazione libera e quella entro i lotti danno p diverse.
    """
    oggetto = genera(tmp_path, EFFETTO_LOTTO)
    esegui(oggetto, _config(), tmp_path / "liberi")
    esegui(oggetto, _config(stat={"permanova_strata": "lotto"}), tmp_path / "strati")
    liberi = righe_permanova(tmp_path / "liberi")[("permanova", "gruppo")]
    strati = righe_permanova(tmp_path / "strati")[("permanova", "gruppo")]
    assert strati["f"] == liberi["f"] and strati["r2"] == liberi["r2"]
    assert float(liberi["p"]) < 0.01 and float(strati["p"]) >= 10 * float(liberi["p"])
    commenti, _ = leggi_tsv(tmp_path / "strati" / "test" / "permanova.tsv")
    assert any("strati entro cui permutano PERMANOVA e dispersioni: lotto" in c
               for c in commenti)
    commenti, _ = leggi_tsv(tmp_path / "liberi" / "test" / "permanova.tsv")
    assert any("strati entro cui permutano PERMANOVA e dispersioni: nessuno" in c
               for c in commenti)

    gruppo = ["a"] * 4 + ["b"] * 4 + ["c"] * 4
    annidati = _avvisi_del_disegno(
        tmp_path, {"gruppo": gruppo, "lotto": ["x"] * 4 + ["x", "y"] * 4},
        stat={"permanova_strata": "lotto"})
    assert "E-ECO-18" in annidati and "strato" in annidati["E-ECO-18"][0]

    with pytest.raises(ErroreEco) as rifiuto:
        valida_sull_oggetto(oggetto, _config(stat={"permanova_strata": "assente"}))
    assert rifiuto.value.codici == ["E-ECO-05"]


def test_con_gli_strati_le_dispersioni_permutano_come_la_permanova(eco_r, tmp_path):
    """
    **Obiettivo**: con sei strati da due campioni, uno per gruppo, le
    permutazioni entro gli strati sono 2^6 = 64: PERMANOVA e test delle
    dispersioni le enumerano entrambe e riportano le stesse 63 permutazioni
    (tutte meno l'identita'), non le 999 richieste; senza strati entrambe ne
    usano 999. La statistica F delle dispersioni non cambia con gli strati. Con
    un campione per strato non esiste alcuna permutazione: nessuno dei due test
    ha una p, e l'analisi si conclude.

    **Razionale scientifico e sistemistico**: i due test si leggono insieme e
    devono rispondere alla stessa ipotesi di scambiabilita': se le dispersioni
    permutassero liberamente mentre la PERMANOVA permuta entro i lotti, la
    verifica non riguarderebbe la stessa distribuzione di riferimento. Il
    numero di permutazioni enumerate rivela lo schema usato.
    """
    blocchi = []
    for strato in "pqrstu":
        blocchi.append({"n": 1, "gruppo": "a", "lotto": strato, "sigma": 0.3})
        blocchi.append({"n": 1, "gruppo": "b", "lotto": strato, "sigma": 0.3, "scambia": True})
    oggetto = genera(tmp_path, blocchi)
    esegui(oggetto, _config(), tmp_path / "liberi")
    esegui(oggetto, _config(stat={"permanova_strata": "lotto"}), tmp_path / "strati")
    liberi, strati = righe_permanova(tmp_path / "liberi"), righe_permanova(tmp_path / "strati")
    assert liberi[("permanova", "gruppo")]["permutazioni"] == "999"
    assert liberi[("dispersioni", "gruppo")]["permutazioni"] == "999"
    assert strati[("permanova", "gruppo")]["permutazioni"] == "63"
    assert strati[("dispersioni", "gruppo")]["permutazioni"] == "63"
    assert strati[("dispersioni", "gruppo")]["f"] == liberi[("dispersioni", "gruppo")]["f"]
    assert strati[("dispersioni", "gruppo")]["campioni"] == "12"

    unici = costruisci(tmp_path, CONTEGGI_12, {
        "gruppo": ["a", "b"] * 6, "lotto": [f"s{i:02d}" for i in range(12)]})
    esito = esegui(unici, _config(stat={"permanova_strata": "lotto"}), tmp_path / "unici")
    assert "E-ECO-16" not in _codici(esito)
    righe = righe_permanova(tmp_path / "unici")
    for analisi in ("permanova", "dispersioni"):
        riga = righe[(analisi, "gruppo")]
        assert riga["p"] == "NA" and riga["permutazioni"] == "0" and riga["f"] != "NA"
        assert "una sola disposizione" in riga["nota"]


def test_il_grafico_dell_alfa_mostra_i_soli_gruppi_del_confronto(eco_r, tmp_path):
    """
    **Obiettivo**: per ogni variabile il cui confronto dell'alfa diversita' si
    esegue l'analisi produce un grafico PNG, elencato nel manifesto e identico
    byte per byte fra due esecuzioni, senza blocchi di data o di testo; il
    grafico cambia se un gruppo sotto la dimensione minima smette di essere
    escluso; se nessun gruppo entra nei test il grafico non si produce,
    l'analisi si conclude e ``E-ECO-12`` lo dichiara.

    **Razionale scientifico e sistemistico**: il grafico accompagna il
    confronto e deve mostrare gli stessi campioni dei test, altrimenti figura e
    tabella direbbero cose diverse; la riproducibilita' si verifica sui byte.
    """
    blocchi = [
        {"n": 8, "gruppo": "a", "sigma": 0.3},
        {"n": 8, "gruppo": "b", "sigma": 0.3, "scambia": True},
        {"n": 2, "gruppo": "c", "sigma": 0.3},
    ]
    oggetto = genera(tmp_path, blocchi)
    config = _config(design={"other_variables": ["lotto"]})
    esito = esegui(oggetto, config, tmp_path / "a")
    esegui(oggetto, config, tmp_path / "b")
    nome = "alfa/gruppi_01_gruppo.png"
    assert nome in esito.file
    # La seconda variabile ha un solo gruppo: nessun confronto, nessun grafico.
    assert not (tmp_path / "a" / "alfa" / "gruppi_02_lotto.png").exists()
    assert "E-ECO-12" in _codici(esito)
    grafico = (tmp_path / "a" / nome).read_bytes()
    assert grafico == (tmp_path / "b" / nome).read_bytes()
    assert grafico.startswith(b"\x89PNG") and len(grafico) > 5000
    for blocco in (b"tIME", b"tEXt", b"iTXt", b"zTXt"):
        assert blocco not in grafico

    # Con il minimo a due il terzo gruppo entra nel confronto e nel grafico.
    esito = esegui(oggetto, _config(stat={"min_group_size": 2}), tmp_path / "tre")
    assert (tmp_path / "tre" / nome).read_bytes() != grafico
    assert "E-ECO-11" not in _codici(esito)

    esito = esegui(oggetto, _config(stat={"min_group_size": 9}), tmp_path / "nessuno")
    assert "E-ECO-12" in _codici(esito) and "E-ECO-16" not in _codici(esito)
    assert not any(f.startswith("alfa/gruppi_") for f in esito.file)
    assert (tmp_path / "nessuno" / "alfa" / "alfa_diversita.tsv").is_file()


def test_le_permutazioni_insufficienti_sono_dichiarate(eco_r, tmp_path):
    """
    **Obiettivo**: con tre campioni contro tre le disposizioni distinte sono
    20, meno delle 999 permutazioni richieste: la validazione lo dichiara con
    ``E-ECO-21`` e la p minima raggiungibile, la tabella riporta le
    permutazioni usate; con 19 permutazioni richieste non lo dichiara.

    **Razionale scientifico e sistemistico**: con pochi campioni la p minima
    raggiungibile e' limitata dal numero di permutazioni che esistono, non da
    quello richiesto.
    """
    oggetto = genera(tmp_path, [
        {"n": 3, "gruppo": "a", "sigma": 0.3},
        {"n": 3, "gruppo": "b", "sigma": 0.3, "scambia": True},
    ])
    piccolo = {"min_group_size": 3}
    esito = esegui(oggetto, _config(stat=piccolo), tmp_path / "uscita")
    avvisi = [a.dettaglio for a in esito.avvisi if a.codice == "E-ECO-21"]
    assert len(avvisi) == 1 and "sono 20," in avvisi[0] and "999" in avvisi[0]
    assert "0.0500000000" in avvisi[0]
    variabile = righe_permanova(tmp_path / "uscita")[("permanova", "gruppo")]
    # Le permutazioni usate sono le 719 dei sei campioni, non le 999 richieste.
    assert variabile["permutazioni"] == "719" and float(variabile["p"]) >= 0.1 - 1e-9
    poche = valida_sull_oggetto(oggetto, _config(stat={**piccolo, "permanova_permutations": 19}))
    assert "E-ECO-21" not in _codici(poche)


def test_un_valore_tecnico_mancante_toglie_il_campione_dalla_sola_permanova(eco_r, tmp_path):
    """
    **Obiettivo**: un campione senza valore nella variabile tecnica esce dalla
    PERMANOVA e dalle dispersioni, con ``E-ECO-13`` e il motivo nella tabella
    dei campioni, e resta nel confronto dell'alfa diversita'.

    **Razionale scientifico e sistemistico**: il modello richiede ogni termine
    per ogni campione; l'esclusione riguarda il solo test che usa quel termine
    e va registrata con il motivo.
    """
    metadati = {
        "gruppo": ["a"] * 6 + ["b"] * 6,
        "lotto": [None, "x", "y", "x", "y", "x"] + ["x", "y"] * 3,
    }
    oggetto = costruisci(tmp_path, CONTEGGI_12, metadati)
    esito = esegui(
        oggetto,
        _config(stat={"min_group_size": 2}, design={"technical_variables": ["lotto"]}),
        tmp_path / "uscita",
    )
    dettagli = [a.dettaglio for a in esito.avvisi if a.codice == "E-ECO-13"]
    assert len(dettagli) == 1 and "variabile tecnica 'lotto'" in dettagli[0] and "c01" in dettagli[0]
    assert righe_permanova(tmp_path / "uscita")[("permanova", "gruppo")]["campioni"] == "11"
    _, campioni = leggi_tsv(tmp_path / "uscita" / "campioni.tsv")
    primo = campioni[0]
    assert primo["stato:gruppo"] == "analizzato"
    assert primo["permanova:gruppo"] == "escluso: valore mancante nella variabile tecnica 'lotto'"
    _, alfa = leggi_tsv(tmp_path / "uscita" / "test" / "alfa_diversita.tsv")
    assert alfa[0]["n_1"] == "6" and alfa[0]["n_2"] == "6"


def test_ogni_variabile_e_una_famiglia_a_se_e_due_esecuzioni_sono_identiche(eco_r, tmp_path):
    """
    **Obiettivo**: una seconda variabile biologica ha le proprie righe e la
    propria famiglia di test; due esecuzioni complete, con ordinazione e test,
    danno tutti i file identici byte per byte e l'oggetto resta immutato.

    **Razionale scientifico e sistemistico**: permutazioni e avvii casuali
    dipendono dal solo seme dichiarato; la correzione per confronti multipli
    non si estende da una variabile all'altra.
    """
    oggetto = genera(tmp_path, [
        {"n": 6, "gruppo": "a", "lotto": "x", "sigma": 0.3},
        {"n": 6, "gruppo": "a", "lotto": "y", "sigma": 0.3},
        {"n": 6, "gruppo": "b", "lotto": "x", "sigma": 0.3, "scambia": True},
        {"n": 6, "gruppo": "b", "lotto": "y", "sigma": 0.3, "scambia": True},
        {"n": 6, "gruppo": "c", "lotto": "x", "sigma": 0.3},
    ], profondita=300)
    prima = hashlib.sha256(oggetto.read_bytes()).hexdigest()
    config = schema.valida(configurazione(
        design={"other_variables": ["lotto"]},
        ord={"methods": ["pcoa", "nmds"], "distances": ["bray", "aitchison"],
             "pcoa_correction": "lingoes", "nmds_trymax": 5},
        stat={"min_group_size": 5, "permanova_permutations": 199},
    ))
    esegui(oggetto, config, tmp_path / "a")
    esegui(oggetto, config, tmp_path / "b")
    assert hashlib.sha256(oggetto.read_bytes()).hexdigest() == prima
    a, b = impronte(tmp_path / "a"), impronte(tmp_path / "b")
    assert a == b
    attesi = {"test/alfa_diversita.tsv", "test/permanova.tsv",
              "ordinazione/pcoa_bray_01_gruppo.png", "ordinazione/nmds_aitchison_02_lotto.png",
              "ordinazione/nmds_bray_coordinate.tsv", "ordinazione/pcoa_aitchison_autovalori.tsv"}
    assert attesi <= set(a)

    _, alfa = leggi_tsv(tmp_path / "a" / "test" / "alfa_diversita.tsv")
    famiglie = {r["variabile"]: r["dimensione_famiglia"] for r in alfa}
    assert famiglie == {"gruppo": "12", "lotto": "3"}
    _, multivariati = leggi_tsv(tmp_path / "a" / "test" / "permanova.tsv")
    assert {(r["variabile"], r["distanza"]) for r in multivariati} == {
        ("gruppo", d) for d in ("bray", "jaccard", "aitchison")} | {
        ("lotto", d) for d in ("bray", "jaccard", "aitchison")}
    assert {r["permutazioni"] for r in multivariati} == {"199"}


def test_con_due_campioni_per_gruppo_le_dispersioni_non_sono_valutabili(eco_r, tmp_path):
    """
    **Obiettivo**: con due campioni per gruppo il test delle dispersioni
    riporta F e p come ``NA`` con la nota «non valutabile» e l'avviso
    ``E-ECO-27``, e non dichiara dispersioni diverse (``E-ECO-19``).

    **Razionale scientifico e sistemistico**: due campioni sono equidistanti
    dal loro centro, la variabilita' residua e' nulla per costruzione e F e'
    un rapporto fra errori di arrotondamento: una p «significativa» sarebbe un
    artefatto.
    """
    oggetto = genera(tmp_path, [
        {"n": 2, "gruppo": "a", "sigma": 0.3},
        {"n": 2, "gruppo": "b", "sigma": 0.3, "scambia": True},
        {"n": 2, "gruppo": "c", "sigma": 0.3},
    ])
    esito = esegui(oggetto, _config(stat={"min_group_size": 2, "permanova_permutations": 99}),
                   tmp_path / "uscita")
    dispersioni = righe_permanova(tmp_path / "uscita")[("dispersioni", "gruppo")]
    assert dispersioni["f"] == "NA" and dispersioni["p"] == "NA"
    assert dispersioni["nota"].startswith("non valutabile")
    assert "E-ECO-27" in _codici(esito) and "E-ECO-19" not in _codici(esito)
    assert "E-ECO-16" not in _codici(esito)


def test_un_modello_saturo_non_viene_stimato_e_lo_si_dichiara(eco_r, tmp_path):
    """
    **Obiettivo**: con una variabile tecnica che ha un valore diverso per ogni
    campione la validazione dichiara ``E-ECO-26``, la PERMANOVA non si esegue
    e la tabella lo dice, senza attribuire zero gradi di liberta' a un termine
    che ne ha; il test delle dispersioni resta.

    **Razionale scientifico e sistemistico**: senza gradi di liberta' residui
    ``adonis2`` restituisce il solo modello complessivo; leggerlo termine per
    termine produrrebbe affermazioni false.
    """
    metadati = {"gruppo": ["a"] * 6 + ["b"] * 6, "lotto": [f"l{i}" for i in range(12)]}
    oggetto = costruisci(tmp_path, CONTEGGI_12, metadati)
    config = _config(stat={"min_group_size": 2}, design={"technical_variables": ["lotto"]})
    assert "E-ECO-26" in _codici(valida_sull_oggetto(oggetto, config))
    esito = esegui(oggetto, config, tmp_path / "uscita")
    righe = righe_permanova(tmp_path / "uscita")
    assert righe[("permanova", "gruppo")]["nota"].startswith("non eseguita")
    assert ("permanova", "lotto") not in righe and ("dispersioni", "gruppo") in righe
    (prova,) = esito.riepilogo["test"]
    assert prova["permanova"]["modello_saturo"] is True
