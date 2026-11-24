"""Filtro delle varianti rare, alfa diversita' senza rarefazione e rarefazione
ripetuta.

Il filtro per abbondanza relativa (``filter.*``) e' un'alternativa, o un
complemento, alla rarefazione: si verifica sulle tre definizioni della soglia,
calcolate qui a mano su tabelle piccole, che le varianti tolte e le letture
tolte siano quelle attese, e che un campione svuotato esca con il suo avviso.
Con ``alpha.method: none`` gli indici si calcolano sui conteggi non rarefatti e
la correlazione con le letture dice se dipendono dalla profondita'. Con
``alpha.rarefy_iterations`` gli indici sono la media di piu' rarefazioni, dallo
stesso seme: stessi byte a ogni esecuzione.
"""

from __future__ import annotations

import pytest

from amplicon16s_eco import config as schema
from amplicon16s_eco.catalogo import CATALOGO, ErroreEco, Specie
from amplicon16s_eco.esecuzione import esegui, valida_sull_oggetto

from eco_aiuti import configurazione, costruisci, eco_r, impronte, leggi_tsv  # noqa: F401


def _codici(dati) -> list[str]:
    with pytest.raises(ErroreEco) as rifiuto:
        schema.valida(dati)
    return rifiuto.value.codici


# --------------------------------------------------------------------------- #
# 1. Lo schema                                                                 #
# --------------------------------------------------------------------------- #


def test_il_metodo_dell_alfa_e_obbligatorio_e_i_suoi_parametri_lo_seguono():
    """
    **Obiettivo**: ``alpha.method`` e' obbligatorio (``E-ECO-02``) e ammette
    ``rarefy`` e ``none`` (altrimenti ``E-ECO-03``); profondita' e numero di
    rarefazioni sono ammessi solo con ``rarefy`` (``E-ECO-01`` con ``none``);
    con ``rarefy`` il numero di rarefazioni vale 1 se omesso e compare nella
    configurazione risolta, con ``none`` non vi compare.

    **Razionale scientifico e sistemistico**: rarefare o no e' una scelta di
    metodo su cui la letteratura non e' concorde: va dichiarata. Un parametro
    della rarefazione dichiarato senza rarefazione entrerebbe nel digest senza
    effetto.
    """
    assert _codici(configurazione(alpha={"method": ...})) == ["E-ECO-02"]
    assert _codici(configurazione(alpha={"method": "rarefazione"})) == ["E-ECO-03"]
    for parametro, valore in (("rarefy_depth", 100), ("rarefy_iterations", 5)):
        assert _codici(configurazione(alpha={"method": "none", parametro: valore})) == [
            "E-ECO-01"]
    assert _codici(configurazione(alpha={"rarefy_iterations": 0})) == ["E-ECO-01"]
    rarefatta = schema.risolta(schema.valida(configurazione()))
    assert rarefatta["alpha"] == {"method": "rarefy", "rarefy_depth": None,
                                  "rarefy_iterations": 1}
    senza = schema.risolta(schema.valida(configurazione(alpha={"method": "none"})))
    assert senza["alpha"] == {"method": "none"}
    # Omettere il numero di rarefazioni o dichiarare 1 e' la stessa analisi.
    assert schema.digest(schema.valida(configurazione(alpha={"rarefy_iterations": 1}))) == (
        schema.digest(schema.valida(configurazione())))


def test_il_filtro_e_facoltativo_e_il_suo_ambito_lo_segue():
    """
    **Obiettivo**: senza filtro il gruppo ``filter`` non compare nella
    configurazione risolta; con la soglia, l'ambito e' obbligatorio
    (``E-ECO-02``) e ammette ``sample``, ``total``, ``mean`` (altrimenti
    ``E-ECO-03``); senza la soglia l'ambito e' respinto (``E-ECO-01``); la
    soglia sta fra 0 e 100 esclusi.

    **Razionale scientifico e sistemistico**: le tre definizioni della soglia
    tolgono varianti diverse; senza la soglia l'ambito non ha effetto.
    """
    base = schema.risolta(schema.valida(configurazione()))
    assert "filter" not in base
    assert _codici(configurazione(filter={"min_relative_abundance": 0.1})) == ["E-ECO-02"]
    assert _codici(configurazione(filter={"scope": "sample"})) == ["E-ECO-01"]
    assert _codici(configurazione(filter={"min_relative_abundance": 0.1,
                                          "scope": "campione"})) == ["E-ECO-03"]
    for soglia in (0, 100, -1, "0.1"):
        assert _codici(configurazione(filter={"min_relative_abundance": soglia,
                                              "scope": "sample"})) == ["E-ECO-01"]
    dichiarato = schema.risolta(schema.valida(configurazione(
        filter={"min_relative_abundance": 0.1, "scope": "total"})))
    assert dichiarato["filter"] == {"min_relative_abundance": 0.1, "scope": "total"}
    assert CATALOGO["E-ECO-37"].specie is Specie.AVVISO
    assert CATALOGO["E-ECO-38"].specie is Specie.RIFIUTO
    assert CATALOGO["E-ECO-39"].specie is Specie.AVVISO


# --------------------------------------------------------------------------- #
# 2. I tre filtri, calcolati a mano                                            #
# --------------------------------------------------------------------------- #

#: Cinque varianti per quattro campioni (una riga per variante). Letture per
#: campione: 1000, 1000, 200, 100.
CONTEGGI = [
    [900, 500, 150, 60],
    [50, 300, 40, 30],
    [40, 150, 6, 0],
    [9, 45, 4, 9],
    [1, 5, 0, 1],
]
METADATI = {"gruppo": ["a", "a", "b", "b"]}


def _riferimento(conteggi: list[list[int]], soglia: float, ambito: str) -> list[list[int]]:
    """I conteggi filtrati secondo la definizione, scritta qui in modo
    indipendente: una variante resta se la sua quota (in percentuale) non e'
    sotto la soglia."""
    letture = [sum(c) for c in zip(*conteggi)]
    totale = sum(letture)
    filtrati = []
    for riga in conteggi:
        if ambito == "sample":
            filtrati.append([x if 100 * x / n >= soglia else 0 for x, n in zip(riga, letture)])
        elif ambito == "total":
            filtrati.append(riga if 100 * sum(riga) / totale >= soglia else [0] * len(riga))
        else:
            media = sum(100 * x / n for x, n in zip(riga, letture)) / len(letture)
            filtrati.append(riga if media >= soglia else [0] * len(riga))
    return filtrati


def _esegui(tmp_path, nome, conteggi=CONTEGGI, metadati=METADATI, **modifica):
    cartella = tmp_path / nome
    cartella.mkdir()
    oggetto = costruisci(cartella, conteggi, metadati)
    uscita = cartella / "uscita"
    return oggetto, uscita, esegui(oggetto, schema.valida(configurazione(**modifica)), uscita)


@pytest.mark.parametrize("ambito, soglia, tolte_per_variante", [
    # sample, 5%: s1 perde v03 (40 < 50), v04 e v05, e tiene v02 (50, alla
    # soglia); s2 perde v04 (45 < 50) e v05; s3 perde v03 (6 < 10) e v04; s4
    # tiene v04 (9 >= 5) e perde v05.
    ("sample", 5.0, {"v03": 46, "v04": 58, "v05": 7}),
    # total, 2%: letture totali 2300, soglia 46; v04 (67) resta, v05 (7) esce.
    ("total", 2.0, {"v05": 7}),
    # mean, 3%: v04 ha quote 0,9 + 4,5 + 2 + 9 = 16,4 / 4 = 4,1% e resta; v05
    # 0,1 + 0,5 + 0 + 1 = 1,6 / 4 = 0,4% ed esce.
    ("mean", 3.0, {"v05": 7}),
])
def test_ogni_filtro_toglie_le_varianti_calcolate_a_mano(eco_r, tmp_path, ambito, soglia,
                                                         tolte_per_variante):
    """
    **Obiettivo**: per ognuno dei tre ambiti, le letture tolte per variante
    (``filtro/varianti_filtrate.tsv``) e per campione (``campioni.tsv``) e la
    ricchezza sui conteggi non rarefatti coincidono con quelle calcolate qui a
    mano; una variante esattamente alla soglia resta.

    **Razionale scientifico e sistemistico**: le tre definizioni della soglia
    (per campione, sul totale, sulla media) differiscono proprio sulle
    varianti rare in alcuni campioni e abbondanti in altri: solo il confronto
    con il calcolo indipendente prova quale e' applicata.
    """
    _, uscita, esito = _esegui(
        tmp_path, ambito, alpha={"method": "none"},
        filter={"min_relative_abundance": soglia, "scope": ambito})
    atteso = _riferimento(CONTEGGI, soglia, ambito)
    commenti, varianti = leggi_tsv(uscita / "filtro" / "varianti_filtrate.tsv")
    assert {r["variante"]: int(r["letture_tolte"]) for r in varianti} == tolte_per_variante
    for riga in varianti:
        i = int(riga["variante"][1:]) - 1
        assert int(riga["letture_prima"]) == sum(CONTEGGI[i])
        assert int(riga["letture_tolte"]) == sum(CONTEGGI[i]) - sum(atteso[i])
        assert riga["tolta_da_tutti_i_campioni"] == ("si" if sum(atteso[i]) == 0 else "no")
    assert any(ambito in c and "filter.scope" in c for c in commenti)

    _, campioni = leggi_tsv(uscita / "campioni.tsv")
    for j, riga in enumerate(campioni):
        prima = sum(r[j] for r in CONTEGGI)
        dopo = sum(r[j] for r in atteso)
        assert int(riga["letture"]) == prima
        assert int(riga["letture_tolte_dal_filtro"]) == prima - dopo

    _, alfa = leggi_tsv(uscita / "alfa" / "alfa_diversita.tsv")
    for j, riga in enumerate(alfa):
        assert int(riga["ricchezza_osservata"]) == sum(1 for r in atteso if r[j] > 0)
        assert int(riga["letture"]) == sum(r[j] for r in atteso)
    assert esito.riepilogo["filtro"]["letture_tolte"] == sum(tolte_per_variante.values())
    # Le distanze si calcolano sui conteggi filtrati: senza le varianti tolte da
    # tutti i campioni, Jaccard fra s1 e s4 e' quella dei conteggi attesi.
    presenti = [{i for i, r in enumerate(atteso) if r[j] > 0} for j in range(4)]
    _, jaccard = leggi_tsv(uscita / "beta" / "distanza_jaccard.tsv")
    attesa = 1 - len(presenti[0] & presenti[3]) / len(presenti[0] | presenti[3])
    assert float(jaccard[0]["c04"]) == pytest.approx(attesa, abs=1e-9)


def test_un_campione_svuotato_dal_filtro_esce_con_il_suo_avviso(eco_r, tmp_path):
    """
    **Obiettivo**: un campione con le sole varianti che il filtro toglie esce
    da tutte le analisi con l'avviso ``E-ECO-37`` (gia' nella validazione) e il
    motivo in ``campioni.tsv``; se restano meno di due campioni il filtro e'
    respinto con ``E-ECO-38``.

    **Razionale scientifico e sistemistico**: un campione senza letture non ha
    alfa diversita' ne' distanze; uscire in silenzio lo farebbe sparire dai
    gruppi senza che chi legge sappia perche'.
    """
    conteggi = [
        [500, 400, 0, 300],
        [500, 600, 0, 700],
        [0, 0, 3, 0],
    ]
    filtro = {"min_relative_abundance": 1.0, "scope": "total"}
    oggetto, uscita, esito = _esegui(tmp_path, "svuotato", conteggi, METADATI,
                                     filter=filtro, stat={"min_group_size": 2})
    avvisi = [a for a in esito.avvisi if a.codice == "E-ECO-37"]
    assert len(avvisi) == 1 and "c03" in avvisi[0].dettaglio
    validazione = valida_sull_oggetto(oggetto, schema.valida(configurazione(filter=filtro)))
    assert "E-ECO-37" in {a.codice for a in validazione.avvisi}
    _, campioni = leggi_tsv(uscita / "campioni.tsv")
    stato = {r["campione"]: r["stato"] for r in campioni}
    assert stato["c03"] == "escluso: senza letture dopo il filtro di abbondanza"
    assert {stato[c] for c in ("c01", "c02", "c04")} == {"analizzato"}
    _, alfa = leggi_tsv(uscita / "alfa" / "alfa_diversita.tsv")
    assert [r["campione"] for r in alfa] == ["c01", "c02", "c04"]
    assert esito.riepilogo["campioni"]["esclusi"]["svuotati_dal_filtro"] == 1

    with pytest.raises(ErroreEco) as rifiuto:
        valida_sull_oggetto(oggetto, schema.valida(configurazione(
            filter={"min_relative_abundance": 75.0, "scope": "sample"})))
    assert rifiuto.value.codici == ["E-ECO-38"]


# --------------------------------------------------------------------------- #
# 3. Alfa diversita' senza rarefazione                                         #
# --------------------------------------------------------------------------- #


def test_senza_rarefazione_la_dipendenza_dalla_profondita_e_misurata_e_dichiarata(
    eco_r, tmp_path
):
    """
    **Obiettivo**: con ``alpha.method: none`` gli indici sono calcolati sui
    conteggi non rarefatti (non c'e' ``conteggi_rarefatti.tsv``),
    ``alfa/profondita.tsv`` riporta per ogni indice la correlazione di
    Spearman con le letture, e un indice che cresce con le letture e' segnalato
    con ``E-ECO-39``; l'intestazione dei test dichiara i conteggi non
    rarefatti.

    **Razionale scientifico e sistemistico**: senza rarefazione la ricchezza
    osservata cresce con le letture; se la profondita' differisce fra i
    gruppi, un confronto della ricchezza confronta anche le profondita'. Chi
    sceglie questa strada deve vedere se e' cosi'.
    """
    # Otto campioni: il k-esimo ha 100 * k letture su k + 1 varianti.
    colonne = []
    for k in range(1, 9):
        colonna = [0] * 10
        for i in range(k + 1):
            colonna[i] = 100 * k // (k + 1)
        colonna[0] += 100 * k - sum(colonna)
        colonne.append(colonna)
    conteggi = [list(r) for r in zip(*colonne)]
    _, uscita, esito = _esegui(tmp_path, "profondita", conteggi,
                               {"gruppo": ["a"] * 4 + ["b"] * 4},
                               alpha={"method": "none"})
    assert not (uscita / "alfa" / "conteggi_rarefatti.tsv").exists()
    commenti, alfa = leggi_tsv(uscita / "alfa" / "alfa_diversita.tsv")
    assert any("non rarefatti" in c for c in commenti)
    assert [int(r["ricchezza_osservata"]) for r in alfa] == list(range(2, 10))
    _, profondita = leggi_tsv(uscita / "alfa" / "profondita.tsv")
    per = {r["indice"]: r for r in profondita}
    assert set(per) == {"ricchezza_osservata", "shannon", "gini_simpson", "chao1"}
    ricchezza = per["ricchezza_osservata"]
    assert float(ricchezza["rho"]) == pytest.approx(1.0)
    assert ricchezza["dipende_dalla_profondita"] == "si" and ricchezza["campioni"] == "8"
    avvisi = [a.dettaglio for a in esito.avvisi if a.codice == "E-ECO-39"]
    assert any("ricchezza_osservata" in a for a in avvisi)
    commenti_test, _ = leggi_tsv(uscita / "test" / "alfa_diversita.tsv")
    assert any("conteggi non rarefatti (alpha.method none)" in c for c in commenti_test)
    assert esito.riepilogo["alfa"] == {"metodo": "none", "campioni": 8}
    assert "alfa/profondita.tsv" in esito.riepilogo["file"]


# --------------------------------------------------------------------------- #
# 4. Rarefazione ripetuta                                                      #
# --------------------------------------------------------------------------- #

#: Sei campioni di profondita' diverse, tutti con varianti rare.
DISUGUALI = [
    [40, 80, 120, 30, 60, 90],
    [10, 20, 40, 5, 25, 30],
    [5, 7, 20, 3, 9, 11],
    [3, 2, 9, 2, 4, 6],
    [2, 1, 5, 1, 2, 3],
]


def test_la_rarefazione_ripetuta_e_la_media_delle_rarefazioni_e_riproducibile(eco_r, tmp_path):
    """
    **Obiettivo**: con ``alpha.rarefy_iterations: 4`` ogni indice in
    ``alfa_diversita.tsv`` e' la media dei quattro di
    ``alfa_per_rarefazione.tsv``; la prima rarefazione e' quella di una
    rarefazione singola con lo stesso seme (stessi ``conteggi_rarefatti.tsv``);
    due esecuzioni danno gli stessi byte; un seme diverso da' altri valori.

    **Razionale scientifico e sistemistico**: una sola rarefazione aggiunge la
    variabilita' di un'estrazione; la media di piu' estrazioni la riduce
    (Schloss 2024). Tutte vengono dal seme dichiarato, e il risultato non
    dipende dall'esecuzione.
    """
    metadati = {"gruppo": ["a", "a", "a", "b", "b", "b"]}
    _, singola, _ = _esegui(tmp_path, "singola", DISUGUALI, metadati)
    _, ripetuta, esito = _esegui(tmp_path, "ripetuta", DISUGUALI, metadati,
                                 alpha={"rarefy_iterations": 4})
    _, di_nuovo, _ = _esegui(tmp_path, "di_nuovo", DISUGUALI, metadati,
                             alpha={"rarefy_iterations": 4})
    _, altro_seme, _ = _esegui(tmp_path, "altro_seme", DISUGUALI, metadati,
                               alpha={"rarefy_iterations": 4}, run={"seed": 8})
    a, b = impronte(ripetuta), impronte(di_nuovo)
    assert {k: v for k, v in a.items() if k != "manifest.json"} == {
        k: v for k, v in b.items() if k != "manifest.json"}
    assert (singola / "alfa" / "conteggi_rarefatti.tsv").read_bytes().split(b"\n", 2)[2] == (
        ripetuta / "alfa" / "conteggi_rarefatti.tsv").read_bytes().split(b"\n", 3)[3]

    commenti, per_rarefazione = leggi_tsv(ripetuta / "alfa" / "alfa_per_rarefazione.tsv")
    _, alfa = leggi_tsv(ripetuta / "alfa" / "alfa_diversita.tsv")
    assert len(per_rarefazione) == 4 * 6
    for riga in alfa:
        sue = [r for r in per_rarefazione if r["campione"] == riga["campione"]]
        assert [r["rarefazione"] for r in sue] == ["1", "2", "3", "4"]
        for indice in ("ricchezza_osservata", "shannon", "gini_simpson", "chao1"):
            media = sum(float(r[indice]) for r in sue) / 4
            assert float(riga[indice]) == pytest.approx(media, abs=1e-9)
    commenti_alfa, _ = leggi_tsv(ripetuta / "alfa" / "alfa_diversita.tsv")
    assert any("4 rarefazioni dal seme 7" in c for c in commenti_alfa)
    _, test = leggi_tsv(ripetuta / "test" / "alfa_diversita.tsv")
    commenti_test, _ = leggi_tsv(ripetuta / "test" / "alfa_diversita.tsv")
    assert any("(media di 4 rarefazioni)" in c for c in commenti_test)
    assert esito.riepilogo["alfa"]["rarefazioni"] == 4
    _, con_altro = leggi_tsv(altro_seme / "alfa" / "alfa_per_rarefazione.tsv")
    assert [r["shannon"] for r in con_altro] != [r["shannon"] for r in per_rarefazione]
