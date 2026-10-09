"""Ordinazione delle analisi ecologiche: PCoA e NMDS su matrici note.

Verifica le proprieta' che un'ordinazione corretta deve avere per costruzione:
una PCoA su una distanza euclidea riproduce le distanze di partenza; le
frazioni di varianza sono quelle degli autovalori dichiarati; la correzione
degli autovalori negativi e' quella richiesta ed e' scritta nella tabella; la
NMDS separa due gruppi separati, riporta lo stress e avvisa quando e' alto.
"""

from __future__ import annotations

import math

import pytest

from amplicon16s_eco import config as schema
from amplicon16s_eco.esecuzione import esegui

from eco_aiuti import configurazione, eco_r, genera, leggi_tsv  # noqa: F401

DUE_GRUPPI = [
    {"n": 10, "gruppo": "a", "sigma": 0.2},
    {"n": 10, "gruppo": "b", "sigma": 0.2, "scambia": True},
]


def _config(metodi, distanze, **ord_):
    return schema.valida(configurazione(
        ord={"methods": metodi, "distances": distanze, **ord_},
        stat={"min_group_size": 5, "permanova_permutations": 99},
    ))


def _coordinate(percorso) -> dict[str, list[float]]:
    _, righe = leggi_tsv(percorso)
    return {r["campione"]: [float(v) for c, v in r.items() if c != "campione"] for r in righe}


def _distanze(percorso) -> dict[tuple[str, str], float]:
    _, righe = leggi_tsv(percorso)
    return {(r["campione"], c): float(v) for r in righe for c, v in r.items() if c != "campione"}


def test_la_pcoa_di_una_distanza_euclidea_riproduce_le_distanze(eco_r, tmp_path):
    """
    **Obiettivo**: sulla distanza di Aitchison, che e' euclidea, le distanze
    fra le coordinate della PCoA su tutti gli assi coincidono con la matrice di
    partenza; non ci sono autovalori negativi, la correzione dichiarata non e'
    applicata e la tabella lo dice; le frazioni di varianza sommano a 1 e sono
    decrescenti.

    **Razionale scientifico e sistemistico**: e' la proprieta' che definisce la
    PCoA; se le coordinate non riproducono le distanze, l'ordinazione
    rappresenta un'altra matrice.
    """
    oggetto = genera(tmp_path, DUE_GRUPPI, taxa=12)
    esito = esegui(oggetto, _config(["pcoa"], ["aitchison"], pcoa_correction="cailliez"),
                   tmp_path / "uscita")
    cartella = tmp_path / "uscita"
    coordinate = _coordinate(cartella / "ordinazione" / "pcoa_aitchison_coordinate.tsv")
    attese = _distanze(cartella / "beta" / "distanza_aitchison.tsv")
    for (a, b), attesa in attese.items():
        assert math.dist(coordinate[a], coordinate[b]) == pytest.approx(attesa, abs=1e-6)

    commenti, autovalori = leggi_tsv(cartella / "ordinazione" / "pcoa_aitchison_autovalori.tsv")
    assert any("nessun autovalore negativo" in c and "non e' applicata" in c for c in commenti)
    frazioni = [float(r["frazione_di_varianza"]) for r in autovalori]
    assert sum(frazioni) == pytest.approx(1.0, abs=1e-8)
    assert frazioni == sorted(frazioni, reverse=True)
    assert {r["autovalore_corretto"] for r in autovalori} == {"NA"}
    valori = [float(r["autovalore"]) for r in autovalori]
    assert frazioni[0] == pytest.approx(valori[0] / sum(valori), abs=1e-8)
    (riepilogo,) = esito.riepilogo["ordinazione"]
    assert riepilogo["correzione_applicata"] is False and riepilogo["autovalori_negativi"] == 0
    # Il segno di ogni asse e' fissato: l'elemento di modulo massimo e' positivo.
    for asse in range(len(next(iter(coordinate.values())))):
        colonna = [c[asse] for c in coordinate.values()]
        assert max(colonna, key=abs) > 0
    # I due gruppi, separati per costruzione, stanno ai lati opposti del primo asse.
    primo = {nome: c[0] for nome, c in coordinate.items()}
    segni_a = {primo[f"c{i:03d}"] > 0 for i in range(1, 11)}
    segni_b = {primo[f"c{i:03d}"] > 0 for i in range(11, 21)}
    assert len(segni_a) == 1 and len(segni_b) == 1 and segni_a != segni_b
    assert (cartella / "ordinazione" / "pcoa_aitchison_01_gruppo.png").read_bytes()[:4] == b"\x89PNG"


def _rumore(tmp_path, nome="rumore.rds"):
    """Campioni senza struttura e con molto rumore: Bray-Curtis non e'
    euclidea e due dimensioni non bastano a rappresentarla.
    """
    return genera(tmp_path, [{"n": 15, "gruppo": "a", "sigma": 1.5},
                             {"n": 15, "gruppo": "b", "sigma": 1.5}], nome=nome)


def test_la_correzione_degli_autovalori_negativi_e_quella_dichiarata(eco_r, tmp_path):
    """
    **Obiettivo**: su Bray-Curtis, non euclidea, senza correzione la tabella
    riporta gli autovalori negativi e il loro numero, e la frazione e' sul
    totale di tutti gli autovalori; con ``cailliez`` e con ``lingoes`` gli
    autovalori corretti sono non negativi, la frazione e' calcolata su di essi
    e somma a 1, e l'intestazione nomina la correzione.

    **Razionale scientifico e sistemistico**: con autovalori negativi la
    «varianza spiegata» dipende da come li si tratta; la scelta va dichiarata
    accanto al numero.
    """
    oggetto = _rumore(tmp_path)
    risultati = {}
    for correzione in ("none", "cailliez", "lingoes"):
        uscita = tmp_path / correzione
        esito = esegui(oggetto, _config(["pcoa"], ["bray"], pcoa_correction=correzione), uscita)
        commenti, righe = leggi_tsv(uscita / "ordinazione" / "pcoa_bray_autovalori.tsv")
        risultati[correzione] = (commenti, righe, esito.riepilogo["ordinazione"][0])

    commenti, righe, riepilogo = risultati["none"]
    negativi = sum(float(r["autovalore"]) < 0 for r in righe)
    assert negativi > 0 and riepilogo["autovalori_negativi"] == negativi
    assert any(f"autovalori negativi: {negativi}" in c for c in commenti)
    assert any("negativi compresi" in c for c in commenti)
    totale = sum(float(r["autovalore"]) for r in righe)
    assert float(righe[0]["frazione_di_varianza"]) == pytest.approx(
        float(righe[0]["autovalore"]) / totale, abs=1e-8)

    for correzione in ("cailliez", "lingoes"):
        commenti, righe, riepilogo = risultati[correzione]
        assert riepilogo["correzione_applicata"] is True
        assert any(f"correzione {correzione} applicata" in c for c in commenti)
        corretti = [float(r["autovalore_corretto"]) for r in righe if r["autovalore_corretto"] != "NA"]
        assert min(corretti) >= -1e-9
        frazioni = [float(r["frazione_di_varianza"]) for r in righe if r["frazione_di_varianza"] != "NA"]
        assert sum(frazioni) == pytest.approx(1.0, abs=1e-8)
        assert frazioni[0] == pytest.approx(corretti[0] / sum(corretti), abs=1e-8)
        assert [r["autovalore"] for r in righe] == [r["autovalore"] for r in risultati["none"][1]]
    assert (risultati["cailliez"][1][0]["autovalore_corretto"]
            != risultati["lingoes"][1][0]["autovalore_corretto"])


def test_la_nmds_separa_due_gruppi_separati_e_riporta_lo_stress(eco_r, tmp_path):
    """
    **Obiettivo**: su due gruppi ben separati la NMDS in due dimensioni li
    separa lungo il primo asse, lo stress e' basso e senza avviso,
    l'intestazione riporta avvii, seme, stress e convergenza; stesso seme,
    stesse coordinate; il numero di avvii e' quello dichiarato.

    **Razionale scientifico e sistemistico**: la NMDS dipende da avvii casuali;
    il seme dichiarato deve bastare a riprodurla, e lo stress dice quanto la
    figura rispetta le distanze.
    """
    oggetto = genera(tmp_path, DUE_GRUPPI)
    config = _config(["nmds"], ["bray"], nmds_trymax=10)
    esito = esegui(oggetto, config, tmp_path / "a")
    esegui(oggetto, config, tmp_path / "b")
    nome = "ordinazione/nmds_bray_coordinate.tsv"
    assert (tmp_path / "a" / nome).read_bytes() == (tmp_path / "b" / nome).read_bytes()
    assert ((tmp_path / "a" / "ordinazione" / "nmds_bray_01_gruppo.png").read_bytes()
            == (tmp_path / "b" / "ordinazione" / "nmds_bray_01_gruppo.png").read_bytes())

    commenti, righe = leggi_tsv(tmp_path / "a" / nome)
    assert list(righe[0]) == ["campione", "nmds_1", "nmds_2"] and len(righe) == 20
    (riepilogo,) = esito.riepilogo["ordinazione"]
    assert riepilogo["avvii"] == 10 and riepilogo["stress"] < 0.2
    assert any("avvii casuali: 10; seme 7; stress: " in c for c in commenti)
    assert any("nessuna riscalatura in unita' di semicambiamento" in c for c in commenti)
    assert "E-ECO-17" not in [a.codice for a in esito.avvisi]
    primo = [float(r["nmds_1"]) for r in righe]
    assert len({x > 0 for x in primo[:10]}) == 1 and len({x > 0 for x in primo[10:]}) == 1
    assert (primo[0] > 0) != (primo[10] > 0)


def test_uno_stress_alto_e_dichiarato_con_un_avviso(eco_r, tmp_path):
    """
    **Obiettivo**: su campioni senza struttura e con molto rumore lo stress in
    due dimensioni supera 0,2 e l'analisi lo dichiara con ``E-ECO-17``,
    riportando lo stesso valore nell'intestazione e nel riepilogo.

    **Razionale scientifico e sistemistico**: sopra 0,2 la figura in due
    dimensioni deforma le distanze, e chi la guarda deve saperlo.
    """
    oggetto = _rumore(tmp_path)
    esito = esegui(oggetto, _config(["nmds"], ["bray"], nmds_trymax=10), tmp_path / "uscita")
    (riepilogo,) = esito.riepilogo["ordinazione"]
    assert riepilogo["stress"] > 0.2
    avvisi = [a.dettaglio for a in esito.avvisi if a.codice == "E-ECO-17"]
    assert len(avvisi) == 1 and f"{riepilogo['stress']:.10f}" in avvisi[0]
    commenti, _ = leggi_tsv(tmp_path / "uscita" / "ordinazione" / "nmds_bray_coordinate.tsv")
    assert any(f"stress: {riepilogo['stress']:.10f}" in c for c in commenti)


def test_uno_stress_quasi_nullo_e_dichiarato_con_un_avviso_proprio(eco_r, tmp_path):
    """
    **Obiettivo**: con quattro campioni la NMDS in due dimensioni ha stress
    quasi nullo: l'analisi lo dichiara con ``E-ECO-28``, non con l'avviso
    imprevisto ``E-ECO-16``, e riporta comunque coordinate e stress; lo stesso
    con venti campioni in due gruppi compatti e ben separati; con trenta
    campioni senza struttura l'avviso non compare.

    **Razionale scientifico e sistemistico**: con pochi punti, o con gruppi che
    si riducono a pochi punti, due dimensioni riproducono l'ordine delle
    distanze senza sforzo: lo stress nullo non e' un buon adattamento ma
    l'assenza di informazione, ed e' un esito atteso che chi legge la figura
    deve conoscere.
    """
    from eco_aiuti import costruisci
    conteggi = [[30, 2, 11, 7], [4, 25, 9, 1], [8, 6, 20, 13], [1, 9, 3, 28], [12, 5, 7, 2]]
    oggetto = costruisci(tmp_path, conteggi, {"gruppo": ["a", "a", "b", "b"]})
    esito = esegui(oggetto, schema.valida(configurazione(
        ord={"methods": ["nmds"], "distances": ["bray"], "nmds_trymax": 5})),
        tmp_path / "pochi")
    codici = [a.codice for a in esito.avvisi]
    assert codici.count("E-ECO-28") == 1 and "E-ECO-16" not in codici
    (avviso,) = [a.dettaglio for a in esito.avvisi if a.codice == "E-ECO-28"]
    assert "bray" in avviso and "4 campioni" in avviso
    (riepilogo,) = esito.riepilogo["ordinazione"]
    assert riepilogo["stress_quasi_nullo"] is True and riepilogo["stress"] < 1e-3
    _, righe = leggi_tsv(tmp_path / "pochi" / "ordinazione" / "nmds_bray_coordinate.tsv")
    assert len(righe) == 4

    esito = esegui(genera(tmp_path, DUE_GRUPPI), _config(["nmds"], ["bray"], nmds_trymax=10),
                   tmp_path / "venti")
    codici = [a.codice for a in esito.avvisi]
    assert codici.count("E-ECO-28") == 1 and "E-ECO-16" not in codici

    esito = esegui(_rumore(tmp_path), _config(["nmds"], ["bray"], nmds_trymax=10),
                   tmp_path / "rumore")
    assert "E-ECO-28" not in [a.codice for a in esito.avvisi]
    assert esito.riepilogo["ordinazione"][0]["stress_quasi_nullo"] is False


def test_senza_metodi_o_con_troppo_pochi_campioni_non_c_e_ordinazione(eco_r, tmp_path):
    """
    **Obiettivo**: con ``ord.methods`` vuoto non si produce alcuna ordinazione;
    con tre campioni l'ordinazione richiesta non viene prodotta e l'analisi lo
    dichiara con ``E-ECO-25``, concludendo le altre uscite.

    **Razionale scientifico e sistemistico**: tre punti stanno sempre in un
    piano: un'ordinazione in due dimensioni non direbbe nulla, e non deve
    fermare il resto dell'analisi.
    """
    oggetto = genera(tmp_path, DUE_GRUPPI)
    esito = esegui(oggetto, _config([], []), tmp_path / "senza")
    assert not (tmp_path / "senza" / "ordinazione").exists()
    assert esito.riepilogo["ordinazione"] == []

    tre = genera(tmp_path, [{"n": 3, "gruppo": "a", "sigma": 0.2}], nome="tre.rds")
    esito = esegui(tre, _config(["pcoa", "nmds"], ["bray"], pcoa_correction="none", nmds_trymax=5),
                   tmp_path / "tre")
    assert "E-ECO-25" in [a.codice for a in esito.avvisi]
    assert not (tmp_path / "tre" / "ordinazione").exists()
    assert (tmp_path / "tre" / "beta" / "distanza_bray.tsv").is_file()


def test_una_distanza_tutta_nulla_o_con_un_solo_asse_non_ferma_l_analisi(eco_r, tmp_path):
    """
    **Obiettivo**: quando tutte le varianti sono presenti in tutti i campioni
    la distanza di Jaccard e' nulla ovunque: la sua ordinazione si salta con
    ``E-ECO-25`` e quella delle altre distanze si produce; con due soli punti
    distinti la PCoA ha un asse e viene saltata allo stesso modo.

    **Razionale scientifico e sistemistico**: sono ingressi legittimi
    (comunita' piccole, repliche identiche); l'ordinazione che non esiste si
    dichiara, senza fermare cio' che e' calcolabile.
    """
    from eco_aiuti import costruisci
    conteggi = [[20 + (3 * i + 5 * j) % 11 for j in range(8)] for i in range(5)]
    oggetto = costruisci(tmp_path, conteggi, {"gruppo": ["a"] * 4 + ["b"] * 4})
    esito = esegui(oggetto, schema.valida(configurazione(
        ord={"methods": ["pcoa", "nmds"], "distances": ["jaccard", "bray"],
             "pcoa_correction": "none", "nmds_trymax": 5})), tmp_path / "nulla")
    saltate = [a.dettaglio for a in esito.avvisi if a.codice == "E-ECO-25"]
    assert len(saltate) == 1 and "jaccard" in saltate[0]
    prodotti = sorted(p.name for p in (tmp_path / "nulla" / "ordinazione").iterdir())
    assert "pcoa_bray_coordinate.tsv" in prodotti and "nmds_bray_coordinate.tsv" in prodotti
    assert not any("jaccard" in nome for nome in prodotti)

    coppie = [[10, 10, 1, 1], [1, 1, 10, 10], [5, 5, 5, 5]]
    due_punti = costruisci(tmp_path, coppie, {"gruppo": ["a", "a", "b", "b"]}, nome="coppie.rds")
    esito = esegui(due_punti, schema.valida(configurazione(
        ord={"methods": ["pcoa"], "distances": ["bray"], "pcoa_correction": "none"})),
        tmp_path / "coppie")
    assert any(a.codice == "E-ECO-25" and "un solo asse" in a.dettaglio for a in esito.avvisi)
    assert not (tmp_path / "coppie" / "ordinazione").exists()
