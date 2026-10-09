"""Alfa diversita', composizione e beta diversita': correttezza su dati noti.

Ogni test costruisce una piccola tabella di cui il risultato giusto si calcola
a mano, o di cui una proprieta' e' nota per costruzione, e verifica che le
uscite del pacchetto la rispettino: solo cosi' si prova che un indice o una
distanza sono quelli dichiarati, e non soltanto che un numero viene prodotto.
Verifica inoltre che due esecuzioni diano gli stessi byte e che l'oggetto di
partenza non cambi.
"""

from __future__ import annotations

import hashlib
import json
import math
import re

import pytest

from amplicon16s_eco import config as schema
from amplicon16s_eco.catalogo import ErroreEco
from amplicon16s_eco.esecuzione import esegui

from eco_aiuti import configurazione, costruisci, eco_r, impronte, leggi_tsv  # noqa: F401

TOLLERANZA = 1e-9

#: Tre campioni con la stessa profondita' (20 letture): la rarefazione alla
#: profondita' minima li lascia come sono, e gli indici si calcolano a mano.
CONTEGGI = [
    [10, 0, 20],
    [5, 10, 0],
    [5, 5, 0],
    [0, 5, 0],
]
METADATI = {"gruppo": ["a", "a", "b"]}


def _colonne(conteggi: list[list[int]]) -> list[list[int]]:
    return [list(c) for c in zip(*conteggi, strict=True)]


def _esegui(tmp_path, conteggi, metadati, nome="uscita", tassonomia=None, **modifica):
    oggetto = costruisci(tmp_path, conteggi, metadati, tassonomia=tassonomia)
    uscita = tmp_path / nome
    esito = esegui(oggetto, schema.valida(configurazione(**modifica)), uscita)
    return oggetto, uscita, esito


def _matrice(percorso) -> dict[tuple[str, str], float]:
    _, righe = leggi_tsv(percorso)
    return {(r["campione"], c): float(v) for r in righe for c, v in r.items() if c != "campione"}


def test_indici_e_distanze_coincidono_con_i_valori_calcolati_a_mano(eco_r, tmp_path):
    """
    **Obiettivo**: su una tabella di quattro varianti per tre campioni,
    ricchezza osservata, Shannon (logaritmo naturale), Gini-Simpson,
    Bray-Curtis sulle proporzioni, Jaccard in presenza e assenza e Aitchison
    coincidono con i valori calcolati qui in modo indipendente.

    **Razionale scientifico e sistemistico**: i predefiniti delle funzioni R
    darebbero altre misure (Jaccard quantitativo, Bray-Curtis sui conteggi); il
    confronto con la formula scritta a mano prova quale misura e' calcolata.
    """
    pseudo = 0.5
    _, uscita, _ = _esegui(tmp_path, CONTEGGI, METADATI, beta={"clr_pseudocount": pseudo})
    campioni = _colonne(CONTEGGI)
    nomi = ["c01", "c02", "c03"]

    _, alfa = leggi_tsv(uscita / "alfa" / "alfa_diversita.tsv")
    for riga, conteggi in zip(alfa, campioni, strict=True):
        p = [x / sum(conteggi) for x in conteggi if x > 0]
        assert int(riga["ricchezza_osservata"]) == len(p)
        assert float(riga["shannon"]) == pytest.approx(-sum(q * math.log(q) for q in p), abs=TOLLERANZA)
        assert float(riga["gini_simpson"]) == pytest.approx(1 - sum(q * q for q in p), abs=TOLLERANZA)
    # Valori scritti per esteso per il primo campione (proporzioni 1/2, 1/4, 1/4).
    assert float(alfa[0]["shannon"]) == pytest.approx(1.5 * math.log(2), abs=TOLLERANZA)
    assert float(alfa[0]["gini_simpson"]) == pytest.approx(0.625, abs=TOLLERANZA)
    assert float(alfa[2]["shannon"]) == 0 and float(alfa[2]["gini_simpson"]) == 0

    bray = _matrice(uscita / "beta" / "distanza_bray.tsv")
    jaccard = _matrice(uscita / "beta" / "distanza_jaccard.tsv")
    aitchison = _matrice(uscita / "beta" / "distanza_aitchison.tsv")
    for i, a in enumerate(campioni):
        for j, b in enumerate(campioni):
            pa, pb = [x / sum(a) for x in a], [x / sum(b) for x in b]
            atteso_bray = sum(abs(x - y) for x, y in zip(pa, pb)) / (sum(pa) + sum(pb))
            presenti_a, presenti_b = {k for k, x in enumerate(a) if x}, {k for k, x in enumerate(b) if x}
            atteso_jaccard = 1 - len(presenti_a & presenti_b) / len(presenti_a | presenti_b)
            la, lb = [math.log(x + pseudo) for x in a], [math.log(x + pseudo) for x in b]
            ca, cb = [x - sum(la) / len(la) for x in la], [x - sum(lb) / len(lb) for x in lb]
            atteso_aitchison = math.sqrt(sum((x - y) ** 2 for x, y in zip(ca, cb)))
            coppia = (nomi[i], nomi[j])
            assert bray[coppia] == pytest.approx(atteso_bray, abs=TOLLERANZA)
            assert jaccard[coppia] == pytest.approx(atteso_jaccard, abs=TOLLERANZA)
            assert aitchison[coppia] == pytest.approx(atteso_aitchison, abs=TOLLERANZA)
    # Valori scritti per esteso: Bray-Curtis 1/2 e Jaccard 1/2 fra i primi due
    # campioni (due varianti condivise su quattro), Jaccard 2/3 fra primo e terzo.
    assert bray[("c01", "c02")] == pytest.approx(0.5, abs=TOLLERANZA)
    assert jaccard[("c01", "c02")] == pytest.approx(0.5, abs=TOLLERANZA)
    assert jaccard[("c01", "c03")] == pytest.approx(2 / 3, abs=TOLLERANZA)
    assert jaccard[("c02", "c03")] == pytest.approx(1.0, abs=TOLLERANZA)


def test_l_intestazione_della_tabella_dichiara_nome_e_formula_di_ogni_indice(eco_r, tmp_path):
    """
    **Obiettivo**: l'intestazione di ``alfa_diversita.tsv`` riporta per ogni
    indice il nome e la formula, e la rarefazione con profondita' e seme.

    **Razionale scientifico e sistemistico**: «Simpson» indica in letteratura
    tre quantita' diverse; una tabella senza la formula non e' interpretabile.
    """
    _, uscita, _ = _esegui(tmp_path, CONTEGGI, METADATI)
    commenti, righe = leggi_tsv(uscita / "alfa" / "alfa_diversita.tsv")
    testo = "\n".join(commenti)
    assert "ricchezza_osservata: numero di varianti con almeno una lettura" in testo
    assert "shannon: -sum(p_i * ln(p_i))" in testo and "logaritmo naturale" in testo
    assert "gini_simpson: 1 - sum(p_i^2)" in testo
    assert "senza reinserimento, profondita' 20" in testo and "seme 7" in testo
    assert list(righe[0]) == ["campione", "letture", "ricchezza_osservata", "shannon", "gini_simpson"]


#: Campioni di profondita' diversa: 30, 60, 300 e 12 letture.
DISUGUALI = [
    [1, 20, 100, 6],
    [2, 20, 100, 3],
    [27, 10, 50, 3],
    [0, 10, 50, 0],
]
METADATI_4 = {"gruppo": ["a", "a", "b", "b"]}


def _rarefatti(uscita) -> dict[str, list[int]]:
    _, righe = leggi_tsv(uscita / "alfa" / "conteggi_rarefatti.tsv")
    return {c: [int(r[c]) for r in righe] for c in righe[0] if c != "variante"}


def test_la_rarefazione_e_senza_reinserimento_alla_profondita_scelta(eco_r, tmp_path):
    """
    **Obiettivo**: con ``alpha.rarefy_depth`` 30 ogni campione rarefatto ha 30
    letture, nessuna variante supera i conteggi di partenza, il campione che ha
    esattamente 30 letture resta identico, e il campione con 12 letture esce
    dall'alfa diversita' con ``E-ECO-14`` restando nelle altre analisi.

    **Razionale scientifico e sistemistico**: senza reinserimento la
    rarefazione e' un sottocampione delle letture osservate; con il
    reinserimento, predefinito di phyloseq, una variante puo' superare i propri
    conteggi e il campione gia' alla profondita' scelta cambierebbe.
    """
    _, uscita, esito = _esegui(tmp_path, DISUGUALI, METADATI_4, alpha={"rarefy_depth": 30})
    rarefatti = _rarefatti(uscita)
    originali = dict(zip(["c01", "c02", "c03", "c04"], _colonne(DISUGUALI), strict=True))
    assert sorted(rarefatti) == ["c01", "c02", "c03"]
    for campione, conteggi in rarefatti.items():
        assert sum(conteggi) == 30
        assert all(r <= o for r, o in zip(conteggi, originali[campione], strict=True))
    assert rarefatti["c01"] == originali["c01"]

    assert [a.codice for a in esito.avvisi].count("E-ECO-14") == 1
    assert esito.riepilogo["alfa"]["profondita_di_rarefazione"] == 30
    assert esito.riepilogo["alfa"]["esclusi_sotto_la_profondita"] == ["c04"]
    _, campioni = leggi_tsv(uscita / "campioni.tsv")
    stato = {r["campione"]: r["stato_alfa"] for r in campioni}
    assert stato["c04"].startswith("escluso: sotto la profondita'") and stato["c01"] == "analizzato"
    _, bray = leggi_tsv(uscita / "beta" / "distanza_bray.tsv")
    assert [r["campione"] for r in bray] == ["c01", "c02", "c03", "c04"]


def test_per_difetto_la_profondita_e_la_minima_dei_campioni_analizzati(eco_r, tmp_path):
    """
    **Obiettivo**: senza ``alpha.rarefy_depth`` la profondita' e' la minima fra
    i campioni analizzati (12), nessun campione e' escluso e il campione minimo
    resta identico; una profondita' che lascia meno di due campioni e' respinta
    con ``E-ECO-10``.

    **Razionale scientifico e sistemistico**: il predefinito e' quello di
    phyloseq, dichiarato; la profondita' usata e' registrata perche' gli indici
    dipendono da essa.
    """
    _, uscita, esito = _esegui(tmp_path, DISUGUALI, METADATI_4)
    rarefatti = _rarefatti(uscita)
    assert {sum(c) for c in rarefatti.values()} == {12} and len(rarefatti) == 4
    assert rarefatti["c04"] == _colonne(DISUGUALI)[3]
    assert esito.riepilogo["alfa"]["origine"] == "minima fra i campioni analizzati"
    assert "E-ECO-14" not in [a.codice for a in esito.avvisi]

    with pytest.raises(ErroreEco) as rifiuto:
        _esegui(tmp_path, DISUGUALI, METADATI_4, nome="troppo", alpha={"rarefy_depth": 100})
    assert rifiuto.value.codici == ["E-ECO-10"]
    assert not (tmp_path / "troppo").exists()


def test_stesso_seme_stessi_conteggi_e_seme_diverso_conteggi_diversi(eco_r, tmp_path):
    """
    **Obiettivo**: due rarefazioni con lo stesso seme danno gli stessi conteggi
    e gli stessi indici; con un seme diverso i conteggi cambiano mentre le
    distanze, calcolate sui conteggi non rarefatti, restano identiche.

    **Razionale scientifico e sistemistico**: la rarefazione e' la sola
    operazione casuale di queste analisi; il seme dichiarato deve bastare a
    riprodurla, e non deve influire su cio' che non e' casuale.
    """
    grandi = [[40 + 3 * i + j for j in range(4)] for i in range(12)]
    grandi[0] = [400, 5, 5, 1]
    alfa = {"rarefy_depth": 50}
    _, prima, _ = _esegui(tmp_path, grandi, METADATI_4, nome="a", alpha=alfa)
    _, seconda, _ = _esegui(tmp_path, grandi, METADATI_4, nome="b", alpha=alfa)
    _, terza, _ = _esegui(tmp_path, grandi, METADATI_4, nome="c", alpha=alfa, run={"seed": 8})
    assert _rarefatti(prima) == _rarefatti(seconda)
    assert _rarefatti(prima) != _rarefatti(terza)
    for nome in ("distanza_bray.tsv", "distanza_jaccard.tsv", "distanza_aitchison.tsv"):
        assert (prima / "beta" / nome).read_bytes() == (terza / "beta" / nome).read_bytes()


def test_i_taxa_non_assegnati_restano_in_una_categoria_e_le_proporzioni_sommano_a_uno(eco_r, tmp_path):
    """
    **Obiettivo**: le varianti senza genere (``NA`` o vuoto), anche di famiglie
    diverse, confluiscono nella sola categoria «non assegnato»; due varianti
    dello stesso genere si sommano; due generi con lo stesso nome in famiglie
    diverse restano distinti; ogni campione somma a 1; l'intestazione dichiara
    la categoria.

    **Razionale scientifico e sistemistico**: con il predefinito di
    ``tax_glom`` le varianti senza genere sparirebbero e le proporzioni degli
    altri taxa risulterebbero gonfiate, senza alcun segnale.
    """
    conteggi = [
        [10, 0, 5, 5],   # GenA, famiglia F1
        [10, 20, 5, 5],  # GenA, famiglia F1: si somma alla precedente
        [20, 20, 10, 0], # senza genere, famiglia F1
        [20, 0, 10, 0],  # senza genere (vuoto), famiglia F2
        [20, 40, 5, 5],  # GenB, famiglia F2
        [20, 20, 5, 5],  # GenB, famiglia F3: omonimo in un'altra famiglia
    ]
    tassonomia = {
        "Regno": ["Batteri"] * 6,
        "Famiglia": ["F1", "F1", "F1", "F2", "F2", "F3"],
        "Genere": ["GenA", "GenA", None, "", "GenB", "GenB"],
    }
    _, uscita, esito = _esegui(tmp_path, conteggi, METADATI_4, tassonomia=tassonomia,
                               comp={"top_n": 2})
    commenti, righe = leggi_tsv(uscita / "composizione" / "abbondanze_relative.tsv")
    assert any("'non assegnato': tutte le varianti senza assegnazione al rango Genere" in c
               for c in commenti)
    per_lignaggio = {r["lignaggio"]: r for r in righe}
    assert sorted(per_lignaggio) == sorted([
        "Batteri;F1;GenA", "Batteri;F2;GenB", "Batteri;F3;GenB", "non assegnato"])
    assert per_lignaggio["non assegnato"]["taxon"] == "non assegnato"
    attese = {
        "Batteri;F1;GenA": [0.2, 0.2, 0.25, 0.5],
        "non assegnato": [0.4, 0.2, 0.5, 0.0],
        "Batteri;F2;GenB": [0.2, 0.4, 0.125, 0.25],
        "Batteri;F3;GenB": [0.2, 0.2, 0.125, 0.25],
    }
    campioni = ["c01", "c02", "c03", "c04"]
    for lignaggio, valori in attese.items():
        for campione, atteso in zip(campioni, valori, strict=True):
            assert float(per_lignaggio[lignaggio][campione]) == pytest.approx(atteso, abs=TOLLERANZA)
    for campione in campioni:
        assert sum(float(r[campione]) for r in righe) == pytest.approx(1.0, abs=1e-8)
    # Righe in ordine decrescente di abbondanza media.
    medie = [float(r["abbondanza_media"]) for r in righe]
    assert medie == sorted(medie, reverse=True)
    assert righe[0]["lignaggio"] == "Batteri;F1;GenA" and medie[0] == pytest.approx(0.2875, abs=TOLLERANZA)
    assert esito.riepilogo["composizione"]["taxa"] == 4

    # Per gruppo: i primi due taxa e «altri», con somma 1 in ogni gruppo.
    commenti, gruppi = leggi_tsv(uscita / "composizione" / "gruppi_01_gruppo.tsv")
    assert [r["taxon"] for r in gruppi] == ["GenA", "non assegnato", "altri"]
    for gruppo in ("a", "b"):
        assert sum(float(r[gruppo]) for r in gruppi) == pytest.approx(1.0, abs=1e-8)
    assert float(gruppi[0]["a"]) == pytest.approx(0.2, abs=TOLLERANZA)
    assert float(gruppi[0]["b"]) == pytest.approx(0.375, abs=TOLLERANZA)
    assert (uscita / "composizione" / "gruppi_01_gruppo.png").read_bytes()[:8] == b"\x89PNG\r\n\x1a\n"


def test_due_esecuzioni_danno_gli_stessi_byte_e_l_oggetto_resta_immutato(eco_r, tmp_path):
    """
    **Obiettivo**: due esecuzioni con lo stesso oggetto, la stessa
    configurazione e lo stesso seme producono gli stessi file, tutti identici
    byte per byte, grafici compresi; lo SHA-256 dell'oggetto e' lo stesso prima
    e dopo, e il manifesto riporta dimensione e SHA-256 di ogni file, il digest
    della configurazione e lo SHA-256 dell'oggetto.

    **Razionale scientifico e sistemistico**: la riproducibilita' si verifica
    sui byte; l'oggetto consegnato dalla pipeline e' un artefatto con un
    checksum pubblicato e non deve cambiare.
    """
    config = schema.valida(configurazione(design={"other_variables": ["lotto"]}))
    oggetto = costruisci(tmp_path, DISUGUALI, {**METADATI_4, "lotto": ["x", "y", "x", "y"]})
    prima = hashlib.sha256(oggetto.read_bytes()).hexdigest()
    esito = esegui(oggetto, config, tmp_path / "a")
    esegui(oggetto, config, tmp_path / "b")
    assert hashlib.sha256(oggetto.read_bytes()).hexdigest() == prima
    assert esito.oggetto_sha256 == f"sha256:{prima}"

    a, b = impronte(tmp_path / "a"), impronte(tmp_path / "b")
    assert a == b and len(a) >= 12
    assert {"composizione/gruppi_01_gruppo.png", "composizione/gruppi_02_lotto.png"} <= set(a)
    assert sorted(a) == sorted(esito.file)

    manifesto = json.loads((tmp_path / "a" / "manifest.json").read_text(encoding="utf-8"))
    assert manifesto["oggetto_di_partenza"]["sha256"] == f"sha256:{prima}"
    assert manifesto["configurazione"]["digest"] == schema.digest(config)
    elencati = {f["nome"]: f for f in manifesto["file"]}
    assert set(elencati) == set(a) - {"manifest.json"}
    for nome, voce in elencati.items():
        assert voce["sha256"] == f"sha256:{a[nome]}"
        assert voce["byte"] == (tmp_path / "a" / nome).stat().st_size
    registrata = json.loads((tmp_path / "a" / "configurazione.json").read_text(encoding="utf-8"))
    assert registrata["digest"] == schema.digest(config)
    assert registrata["parametri"] == schema.risolta(config)
    # Nessun percorso della macchina nelle uscite di testo.
    for nome in a:
        if not nome.endswith(".png"):
            assert str(tmp_path) not in (tmp_path / "a" / nome).read_text(encoding="utf-8")
    # Nessun file di lavoro accanto all'oggetto.
    assert sorted(p.name for p in oggetto.parent.iterdir()) == ["oggetto.rds", "oggetto.rds.json"]


def test_le_tabelle_hanno_formato_fisso_e_ordine_deterministico(eco_r, tmp_path):
    """
    **Obiettivo**: i numeri reali hanno dieci decimali con il punto, senza
    notazione esponenziale; campioni e colonne sono in ordine alfabetico anche
    se l'oggetto li ha in un altro ordine.

    **Razionale scientifico e sistemistico**: un formato che dipende dalle
    opzioni della sessione R o dalla lingua della macchina renderebbe le
    tabelle diverse a parita' di risultato.
    """
    _, uscita, _ = _esegui(tmp_path, DISUGUALI, METADATI_4)
    for nome in ("alfa/alfa_diversita.tsv", "beta/distanza_bray.tsv",
                 "beta/distanza_aitchison.tsv", "composizione/abbondanze_relative.tsv"):
        _, righe = leggi_tsv(uscita / nome)
        for riga in righe:
            for colonna, valore in riga.items():
                if colonna in ("campione", "taxon", "lignaggio", "letture", "ricchezza_osservata"):
                    continue
                assert re.fullmatch(r"-?\d+\.\d{10}", valore), (nome, colonna, valore)
    _, bray = leggi_tsv(uscita / "beta" / "distanza_bray.tsv")
    assert [r["campione"] for r in bray] == sorted(r["campione"] for r in bray)
    assert list(bray[0])[1:] == [r["campione"] for r in bray]
    assert all(float(r[r["campione"]]) == 0 for r in bray)


def test_unifrac_non_pesato_coincide_con_la_frazione_di_rami_non_condivisi(eco_r, tmp_path):
    """
    **Obiettivo**: su un albero radicato di quattro foglie con rami di
    lunghezza 1, UniFrac non pesato vale 1 fra due campioni senza rami in
    comune e 3/5 fra due campioni che condividono due rami su cinque; quello
    pesato e' zero sulla diagonale e simmetrico.

    **Razionale scientifico e sistemistico**: UniFrac dipende dall'albero e
    dalla sua radice; il valore calcolato sui rami prova che l'albero
    dell'oggetto e' quello usato.
    """
    conteggi = [
        [10, 0, 10],
        [10, 0, 0],
        [0, 10, 10],
        [0, 10, 0],
    ]
    oggetto = costruisci(tmp_path, conteggi, {"gruppo": ["a", "a", "b"]},
                         albero="((v01:1,v02:1):1,(v03:1,v04:1):1);")
    config = schema.valida(configurazione(beta={
        "distances": ["unifrac_unweighted", "unifrac_weighted"], "clr_pseudocount": ...,
        "unifrac_root": "existing"}))
    uscita = tmp_path / "uscita"
    esegui(oggetto, config, uscita)
    non_pesato = _matrice(uscita / "beta" / "distanza_unifrac_unweighted.tsv")
    assert non_pesato[("c01", "c02")] == pytest.approx(1.0, abs=TOLLERANZA)
    assert non_pesato[("c01", "c03")] == pytest.approx(0.6, abs=TOLLERANZA)
    assert non_pesato[("c02", "c03")] == pytest.approx(0.6, abs=TOLLERANZA)
    pesato = _matrice(uscita / "beta" / "distanza_unifrac_weighted.tsv")
    assert pesato[("c01", "c01")] == 0
    assert pesato[("c01", "c02")] == pytest.approx(pesato[("c02", "c01")], abs=TOLLERANZA)
    assert pesato[("c01", "c02")] > pesato[("c01", "c03")] > 0


def test_un_albero_senza_radice_si_radica_al_punto_medio_sulla_copia(eco_r, tmp_path):
    """
    **Obiettivo**: su un albero senza radice, con ``beta.unifrac_root:
    midpoint`` le due UniFrac si calcolano e coincidono con quelle ottenute,
    con ``existing``, sullo stesso albero scritto gia' radicato al punto medio
    del cammino piu' lungo; due esecuzioni danno gli stessi byte in ogni file;
    l'oggetto di partenza non cambia; intestazione delle tabelle, riepilogo e
    configurazione registrata riportano il modo di radicare. Lo stesso albero
    con ``existing`` e' respinto con ``E-ECO-06``.

    **Razionale scientifico e sistemistico**: la radice decide quali rami sono
    condivisi fra due campioni; il punto medio e' un procedimento senza numeri
    casuali, verificabile a mano su un albero piccolo, e la radicazione non
    deve toccare l'oggetto consegnato dalla pipeline.
    """
    conteggi = [
        [10, 0, 10, 3],
        [10, 0, 0, 5],
        [0, 10, 10, 2],
        [0, 10, 0, 7],
    ]
    metadati = {"gruppo": ["a", "a", "b", "b"]}
    # Il cammino piu' lungo va da v04 a v01 (o v02) ed e' lungo 5: il punto
    # medio cade sul ramo di v04, a 2,5 dalla foglia.
    senza_radice = costruisci(tmp_path, conteggi, metadati,
                              albero="(v01:1,v02:1,(v03:1,v04:3):1);")
    radicato = costruisci(tmp_path, conteggi, metadati, nome="radicato.rds",
                          albero="(((v01:1,v02:1):1,v03:1):0.5,v04:2.5);")
    distanze = ["unifrac_unweighted", "unifrac_weighted"]

    def config(radice):
        return schema.valida(configurazione(beta={
            "distances": distanze, "clr_pseudocount": ..., "unifrac_root": radice}))

    prima = hashlib.sha256(senza_radice.read_bytes()).hexdigest()
    esito = esegui(senza_radice, config("midpoint"), tmp_path / "a")
    esegui(senza_radice, config("midpoint"), tmp_path / "b")
    assert hashlib.sha256(senza_radice.read_bytes()).hexdigest() == prima
    assert impronte(tmp_path / "a") == impronte(tmp_path / "b")

    esegui(radicato, config("existing"), tmp_path / "atteso")
    for distanza in distanze:
        nome = f"distanza_{distanza}.tsv"
        calcolata = _matrice(tmp_path / "a" / "beta" / nome)
        attesa = _matrice(tmp_path / "atteso" / "beta" / nome)
        assert calcolata == pytest.approx(attesa, abs=TOLLERANZA)
        assert max(calcolata.values()) > 0
        commenti, _ = leggi_tsv(tmp_path / "a" / "beta" / nome)
        assert any("beta.unifrac_root = midpoint" in c and "phangorn::midpoint" in c
                   and "non aveva una radice" in c for c in commenti)
        commenti, _ = leggi_tsv(tmp_path / "atteso" / "beta" / nome)
        assert any("beta.unifrac_root = existing" in c and "aveva gia' una radice" in c
                   for c in commenti)
    assert esito.riepilogo["beta"]["albero"] == {
        "radice": "midpoint", "radicato_nell_oggetto": False}
    assert "phangorn" in esito.riepilogo["ambiente"]["pacchetti"]
    registrata = json.loads((tmp_path / "a" / "configurazione.json").read_text(encoding="utf-8"))
    assert registrata["parametri"]["beta"]["unifrac_root"] == "midpoint"

    with pytest.raises(ErroreEco) as rifiuto:
        esegui(senza_radice, config("existing"), tmp_path / "respinta")
    assert rifiuto.value.codici == ["E-ECO-06"]
    assert not (tmp_path / "respinta").exists()


def test_una_tassonomia_a_un_solo_rango_si_agglomera_per_nome(eco_r, tmp_path):
    """
    **Obiettivo**: con un solo rango tassonomico le varianti con lo stesso
    nome si sommano, quelle senza nome vanno in «non assegnato» e ogni
    campione somma a 1.

    **Razionale scientifico e sistemistico**: ``tax_glom`` non tratta una
    tassonomia a un solo rango; il risultato deve essere lo stesso che darebbe
    su piu' ranghi.
    """
    conteggi = [[10, 0, 5, 5], [10, 20, 5, 5], [20, 20, 10, 0], [20, 0, 10, 0]]
    oggetto = costruisci(tmp_path, conteggi, METADATI_4,
                         tassonomia={"Genere": ["GenA", "GenA", None, "GenB"]})
    esegui(oggetto, schema.valida(configurazione()), tmp_path / "uscita")
    _, righe = leggi_tsv(tmp_path / "uscita" / "composizione" / "abbondanze_relative.tsv")
    per_taxon = {r["taxon"]: r for r in righe}
    assert sorted(per_taxon) == ["GenA", "GenB", "non assegnato"]
    assert float(per_taxon["GenA"]["c01"]) == pytest.approx(1 / 3, abs=TOLLERANZA)
    assert float(per_taxon["non assegnato"]["c02"]) == pytest.approx(0.5, abs=TOLLERANZA)
    for campione in ("c01", "c02", "c03", "c04"):
        assert sum(float(r[campione]) for r in righe) == pytest.approx(1.0, abs=1e-8)


def test_la_distanza_di_aitchison_non_dipende_dai_campioni_esclusi(eco_r, tmp_path):
    """
    **Obiettivo**: la distanza di Aitchison fra i campioni di un sottoinsieme
    e' la stessa che si ottiene da un oggetto con quei soli campioni, anche se
    l'oggetto intero ha varianti presenti solo nei campioni esclusi;
    l'intestazione riporta quante varianti entrano nel CLR.

    **Razionale scientifico e sistemistico**: il CLR divide per la media
    geometrica su tutte le varianti; contare varianti assenti dai campioni
    analizzati farebbe dipendere la distanza da cio' che non si analizza.
    """
    conteggi = [
        [10, 4, 7, 0, 0],
        [5, 9, 3, 0, 0],
        [8, 2, 6, 1, 0],
        [0, 0, 0, 30, 40],
        [0, 0, 0, 25, 10],
    ]
    metadati = {"gruppo": ["a", "a", "b", "b", "b"], "lotto": ["x", "x", "x", "y", "y"]}
    intero = costruisci(tmp_path, conteggi, metadati)
    config = schema.valida(configurazione(design={"subset": {"lotto": ["x"]}}))
    esegui(intero, config, tmp_path / "sottoinsieme")
    ridotto = costruisci(tmp_path, [riga[:3] for riga in conteggi[:3]],
                         {k: v[:3] for k, v in metadati.items()},
                         tassonomia={k: v[:3] for k, v in {
                             "Regno": ["Batteri"] * 5, "Famiglia": ["F"] * 5,
                             "Genere": [f"G{i}" for i in range(5)]}.items()},
                         nome="ridotto.rds")
    esegui(ridotto, schema.valida(configurazione()), tmp_path / "ridotto")
    a = _matrice(tmp_path / "sottoinsieme" / "beta" / "distanza_aitchison.tsv")
    b = _matrice(tmp_path / "ridotto" / "beta" / "distanza_aitchison.tsv")
    assert a.keys() == b.keys()
    for coppia in a:
        assert a[coppia] == pytest.approx(b[coppia], abs=TOLLERANZA)
    commenti, _ = leggi_tsv(tmp_path / "sottoinsieme" / "beta" / "distanza_aitchison.tsv")
    assert any("varianti nel CLR" in c and c.endswith(": 3") for c in commenti)


def test_un_gruppo_che_si_chiama_come_una_colonna_non_sovrascrive_la_tabella(eco_r, tmp_path):
    """
    **Obiettivo**: un gruppo di nome ``taxon`` non sostituisce la colonna dei
    taxa nella tabella per gruppo: le etichette dei taxa restano nella prima
    colonna e i due gruppi hanno ciascuno la propria.

    **Razionale scientifico e sistemistico**: i nomi dei gruppi vengono dai
    metadati dello studio; una coincidenza con un'intestazione non deve
    corrompere una tabella in silenzio.
    """
    oggetto = costruisci(tmp_path, DISUGUALI, {"gruppo": ["taxon", "taxon", "b", "b"]})
    esegui(oggetto, schema.valida(configurazione()), tmp_path / "uscita")
    righe = (tmp_path / "uscita" / "composizione" / "gruppi_01_gruppo.tsv").read_text(
        encoding="utf-8").splitlines()
    corpo = [r.split("\t") for r in righe if not r.startswith("# ")]
    assert corpo[0] == ["taxon", "b", "taxon"]
    assert [r[0] for r in corpo[1:]][-1] == "altri" and corpo[1][0].startswith("Gen")
    assert all(len(r) == 3 for r in corpo)
