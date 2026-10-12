"""Verifica della composizione di una comunita' nota.

Verifica ``dati/comunita_nota/verifica_composizione.py`` su tabelle sintetiche
di cui l'esito si calcola a mano: la soglia di presenza di un genere atteso, la
quota delle letture su generi non attesi con i suoi casi di confine, le
varianti senza genere contate secondo il rango piu' fine che hanno, il
rapporto fra misurato e teorico, e il campione escluso dalla pipeline. Lo
script decide se una comunita' di composizione nota e' stata ricostruita: un
errore a favore non si vedrebbe sui dati reali.
"""

from __future__ import annotations

import csv
import importlib.util
import json
import math
from pathlib import Path

import pytest

RADICE = Path(__file__).resolve().parents[1]
SCRIPT = RADICE / "dati" / "comunita_nota" / "verifica_composizione.py"

pytestmark = pytest.mark.skipif(not SCRIPT.is_file(), reason="dati/comunita_nota non presente")


@pytest.fixture(scope="module")
def vc():
    """Il modulo dello script, caricato dal suo file."""
    specifica = importlib.util.spec_from_file_location("verifica_composizione", SCRIPT)
    modulo = importlib.util.module_from_spec(specifica)
    specifica.loader.exec_module(modulo)
    return modulo


def _tsv(percorso: Path, righe: list[list]) -> None:
    percorso.parent.mkdir(parents=True, exist_ok=True)
    with open(percorso, "w", encoding="utf-8", newline="") as file:
        csv.writer(file, delimiter="\t", lineterminator="\n").writerows(righe)


#: Tre generi attesi, due nella stessa famiglia.
ATTESI = [["genere", "famiglia", "specie", "percentuale", "fonte"],
          ["Aa", "Fa", "Aa prima", "50", "prova"], ["Bb", "Fb", "Bb seconda", "30", "prova"],
          ["Cc-Dd", "Fb", "Cc terza", "20", "prova"]]

#: Le varianti: (genere, famiglia). Le ultime sei non hanno un genere atteso.
VARIANTI = [("Aa", "Fa"), ("g__Bb", "Fb"), ("Cc-Dd", "Fb"), ("Aa", "Fa"),
            ("", "Fb"), ("uncultured", "f__Fa"), ("Zz", "Fb"), ("NA", "Fz"), ("", ""), ("Cc", "Fb")]


def _finale(cartella: Path, colonne: dict[str, list[int]], esclusi: list[list] | None = None) -> Path:
    """Scrive una cartella ``12_final`` con le varianti di :data:`VARIANTI` e,
    per ogni campione, i conteggi dati nello stesso ordine."""
    finale = cartella / "12_final"
    nomi = list(colonne)
    _tsv(finale / "conteggi.tsv", [["asv_id", *[f"ACC{i}" for i in range(len(nomi))]]] + [
        [f"ASV{j + 1}", *[colonne[n][j] for n in nomi]] for j in range(len(VARIANTI))])
    _tsv(finale / "tassonomia.tsv", [["asv_id", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus"]] + [
        [f"ASV{j + 1}", "Bacteria", "P", "C", "O", f, g] for j, (g, f) in enumerate(VARIANTI)])
    _tsv(finale / "metadati.tsv", [["accession", "sample_name", "classe"]] + [
        [f"ACC{i}", n, "biologico"] for i, n in enumerate(nomi)])
    _tsv(finale / "intermedi" / "esclusioni.tsv",
         [["accession", "sample_name", "piastra", "filtro", "motivo"], *(esclusi or [])])
    return finale


def test_la_categoria_di_una_variante_segue_il_rango_piu_fine_assegnato(vc):
    """
    **Obiettivo**: Verificare che una variante conti come attesa solo con il
    nome esatto di un genere atteso (tolto il prefisso di rango); che senza
    genere, o con un'etichetta non informativa, resti «senza genere» se la sua
    famiglia contiene un genere atteso e sia non attesa altrimenti; e che un
    genere diverso, anche della famiglia giusta o con un nome simile, sia non
    atteso.

    **Razionale scientifico e sistemistico**: Con letture corte piu' generi
    della stessa famiglia non si distinguono: fermarsi alla famiglia giusta non
    e' un errore, ma nemmeno una conferma, e va contato a parte. Tutto il resto
    che non e' atteso e' contaminazione o errore di assegnazione.
    """
    attesi = {"Aa": 50.0, "Bb": 30.0, "Cc-Dd": 20.0}
    famiglie = {"Fa", "Fb"}
    c = lambda g, f: vc.categoria(g, f, attesi, famiglie)  # noqa: E731
    assert c("Aa", "Fa") == (vc.ATTESO, "Aa") and c(" g__Bb ", "x") == (vc.ATTESO, "Bb")
    assert c("", "Fb") == (vc.SENZA_GENERE, "Fb") and c("uncultured", "f__Fa") == (vc.SENZA_GENERE, "Fa")
    assert c("NA", "Fz") == (vc.NON_ATTESO, vc.NON_ASSEGNATO)
    assert c("", "") == (vc.NON_ATTESO, vc.NON_ASSEGNATO) and c(None, None)[0] == vc.NON_ATTESO
    assert c("Zz", "Fb") == (vc.NON_ATTESO, "Zz")
    assert c("Cc", "Fb") == (vc.NON_ATTESO, "Cc") and c("aa", "Fa") == (vc.NON_ATTESO, "aa")


def test_i_due_criteri_e_i_rapporti_su_un_campione_calcolato_a_mano(vc, tmp_path, capsys):
    """
    **Obiettivo**: Verificare su tre campioni, dalla riga di comando: le quote
    e i rapporti fra misurato e teorico calcolati a mano; che un genere atteso
    con esattamente lo 0,1% delle letture sia trovato e uno senza letture no;
    che il 2% esatto di letture non attese superi il criterio e il 2,1% no; che
    le varianti senza genere nella famiglia attesa non contino fra le non
    attese; e che la comunita' superi un criterio solo se lo supera ogni
    campione.

    **Razionale scientifico e sistemistico**: Le due soglie decidono l'esito
    della validazione: il loro confine e cio' che entra nel numeratore devono
    essere quelli scritti, non quelli piu' favorevoli.
    """
    _tsv(tmp_path / "attesi.tsv", ATTESI)
    _tsv(tmp_path / "campioni.tsv", [["campione"], ["buono"], ["confine"], ["cattivo"]])
    _tsv(tmp_path / "solo.tsv", [["campione"], ["buono"]])
    colonne = {
        # Aa 300+200, Bb 300, Cc-Dd 190; senza genere in famiglia attesa 3+2;
        # non attese: Zz 3, famiglia ignota 1, non assegnata 1, Cc 0.
        "buono": [300, 300, 190, 200, 3, 2, 3, 1, 1, 0],
        # Cc-Dd con 1 lettura su 1000: lo 0,1% esatto. Non attese 20 su 1000: il 2% esatto.
        "confine": [979, 0, 1, 0, 0, 0, 20, 0, 0, 0],
        # Bb assente; non attese 21 su 1000.
        "cattivo": [500, 0, 479, 0, 0, 0, 0, 0, 0, 21],
    }
    finale = _finale(tmp_path, colonne)
    base = ["--finale", str(finale), "--attesi", str(tmp_path / "attesi.tsv"),
            "--campioni", str(tmp_path / "campioni.tsv")]
    # L'elenco dei campioni attesi e' obbligatorio: senza, una corsa persa non si vedrebbe.
    with pytest.raises(SystemExit):
        vc.main(["--finale", str(finale), "--attesi", str(tmp_path / "attesi.tsv"),
                 "--uscita", str(tmp_path / "z")])
    capsys.readouterr()
    assert vc.main([*base, "--uscita", str(tmp_path / "a")]) == 1
    stampato = capsys.readouterr().out
    assert "cattivo: criterio 1 NON superato (mancano: Bb); criterio 2 NON superato" in stampato
    riepilogo = json.loads((tmp_path / "a" / "riepilogo.json").read_text(encoding="utf-8"))
    assert not riepilogo["criterio_1"]["superato"] and not riepilogo["criterio_2"]["superato"]
    assert riepilogo["famiglie_ammesse"] == ["Fa", "Fb"]
    c = riepilogo["campioni"]

    buono = c["buono"]
    assert buono["letture"] == 1000 and buono["quota_attesi"] == pytest.approx(0.99)
    assert buono["quota_non_attesi"] == pytest.approx(0.005)
    assert buono["quota_senza_genere_in_famiglia_attesa"] == pytest.approx(0.005)
    assert buono["criterio_1"] == {"superato": True, "mancanti": []}
    assert buono["criterio_2"]["superato"]
    assert buono["non_attesi"] == {"Zz": 3, vc.NON_ASSEGNATO: 2}
    generi = {g["genere"]: g for g in buono["generi"]}
    assert generi["Aa"]["misurato"] == pytest.approx(100 * 500 / 990)
    assert generi["Aa"]["teorico"] == pytest.approx(50.0)
    assert generi["Aa"]["rapporto"] == pytest.approx(500 / 990 / 0.5)
    assert generi["Cc-Dd"]["log2"] == pytest.approx(math.log2(190 / 990 / 0.2))
    assert sum(g["misurato"] for g in buono["generi"]) == pytest.approx(100.0)
    assert sum(g["teorico"] for g in buono["generi"]) == pytest.approx(100.0)

    confine = c["confine"]
    assert confine["criterio_1"] == {"superato": False, "mancanti": ["Bb"]}
    assert {g["genere"]: g["trovato"] for g in confine["generi"]} == {"Aa": True, "Bb": False, "Cc-Dd": True}
    assert confine["criterio_2"] == {"superato": True, "valore": pytest.approx(0.02)}
    assert {g["genere"]: g["log2"] for g in confine["generi"]}["Bb"] is None

    cattivo = c["cattivo"]
    assert cattivo["criterio_2"] == {"superato": False, "valore": pytest.approx(0.021)}
    assert cattivo["non_attesi"] == {"Cc": 21}

    # Con i soli campioni che li superano, i criteri della comunita' sono superati.
    finale = _finale(tmp_path / "solo", {"buono": colonne["buono"]})
    assert vc.main(["--finale", str(finale), "--attesi", str(tmp_path / "attesi.tsv"),
                    "--campioni", str(tmp_path / "solo.tsv"), "--uscita", str(tmp_path / "b")]) == 0
    # La stessa cartella, con l'elenco di tutte e tre le corse: le due assenti
    # non superano i criteri, e la comunita' nemmeno.
    assert vc.main(["--finale", str(finale), "--attesi", str(tmp_path / "attesi.tsv"),
                    "--campioni", str(tmp_path / "campioni.tsv"), "--uscita", str(tmp_path / "c")]) == 1
    assert "comunita': criterio 1 superato, criterio 2 superato" in capsys.readouterr().out
    with open(tmp_path / "b" / "composizione.tsv", encoding="utf-8", newline="") as file:
        righe = list(csv.DictReader(file, delimiter="\t"))
    assert [(r["genere"], r["letture"], r["trovato"]) for r in righe] == [
        ("Aa", "500", "si"), ("Bb", "300", "si"), ("Cc-Dd", "190", "si")]


def test_un_campione_escluso_dalla_pipeline_non_supera_i_criteri(vc, tmp_path, capsys):
    """
    **Obiettivo**: Verificare che un campione escluso dai filtri della
    pipeline, o dichiarato nella tabella dei campioni e assente dall'oggetto
    finale, sia riportato con il motivo e faccia fallire i due criteri della
    comunita' anche se gli altri campioni li superano; e che un campione
    senza letture non li superi.

    **Razionale scientifico e sistemistico**: Una replica persa per strada non
    e' una replica riuscita: valutare solo quelle arrivate in fondo darebbe un
    esito migliore del vero.
    """
    _tsv(tmp_path / "attesi.tsv", ATTESI)
    buono = [300, 300, 190, 200, 3, 2, 3, 1, 1, 0]
    finale = _finale(tmp_path, {"buono": buono},
                     esclusi=[["ACC9", "perso", "", "letture_finali", "poche letture"]])
    _tsv(tmp_path / "campioni.tsv", [["campione", "file"], ["buono", "a"], ["perso", "b"], ["mai visto", "c"]])
    assert vc.main(["--finale", str(finale), "--attesi", str(tmp_path / "attesi.tsv"),
                    "--campioni", str(tmp_path / "campioni.tsv"), "--uscita", str(tmp_path / "a")]) == 1
    stampato = capsys.readouterr().out
    assert "perso: escluso dalla pipeline (letture_finali: poche letture)" in stampato
    riepilogo = json.loads((tmp_path / "a" / "riepilogo.json").read_text(encoding="utf-8"))
    assert riepilogo["esclusi_dalla_pipeline"] == {
        "mai visto": "non e' nell'oggetto finale", "perso": "letture_finali: poche letture"}
    assert riepilogo["campioni"]["buono"]["criterio_1"]["superato"]
    assert not riepilogo["criterio_1"]["superato"] and not riepilogo["criterio_2"]["superato"]

    vuoto = vc.valuta({}, {"Aa": 50.0, "Bb": 50.0})
    assert vuoto["letture"] == 0 and not vuoto["criterio_1"]["superato"]
    assert not vuoto["criterio_2"]["superato"] and vuoto["criterio_1"]["mancanti"] == ["Aa", "Bb"]

    # Una tabella degli attesi senza le colonne richieste non si valuta.
    _tsv(tmp_path / "rotta.tsv", [["genere", "percentuale"], ["Aa", "100"]])
    assert vc.main(["--finale", str(finale), "--attesi", str(tmp_path / "rotta.tsv"),
                    "--campioni", str(tmp_path / "campioni.tsv"), "--uscita", str(tmp_path / "b")]) == 2
    assert "verifica non eseguita" in capsys.readouterr().out


def test_gli_altri_nomi_di_un_genere_atteso_si_sommano_al_genere(vc, tmp_path, capsys):
    """
    **Obiettivo**: Verificare che, con la colonna facoltativa ``altri_nomi``,
    le letture delle varianti assegnate a un altro nome dello stesso genere
    atteso si sommino a quelle del genere, nei due criteri e nei rapporti; che
    senza la colonna lo stesso nome conti fra i non attesi e il riepilogo non
    porti la voce degli altri nomi; e che un altro nome ripetuto, o uguale a
    un genere atteso, renda la tabella non valida.

    **Razionale scientifico e sistemistico**: Un riferimento puo' scrivere lo
    stesso organismo con due nomi (un genere diviso in due da una revisione
    della nomenclatura): contarne uno solo darebbe per mancante, o per
    contaminazione, cio' che e' la comunita' stessa. L'equivalenza deve stare
    nella tabella degli attesi, dichiarata, e non nel codice.
    """
    intestazione = ["genere", "famiglia", "specie", "percentuale", "fonte"]
    senza = [intestazione, ["Aa", "Fa", "Aa prima", "50", "prova"], ["Bb", "Fb", "Bb seconda", "30", "prova"],
             ["Cc-Dd", "Fb", "Cc terza", "20", "prova"]]
    con = [[*intestazione, "altri_nomi"], [*senza[1], ""], [*senza[2], " Zz | g__Cc "], [*senza[3], ""]]
    _tsv(tmp_path / "senza.tsv", senza)
    _tsv(tmp_path / "con.tsv", con)
    _tsv(tmp_path / "campioni.tsv", [["campione"], ["uno"]])
    # Aa 500, Bb 0, Cc-Dd 400; Zz 60 e Cc 40, che sono altri nomi di Bb.
    finale = _finale(tmp_path, {"uno": [500, 0, 400, 0, 0, 0, 60, 0, 0, 40]})
    base = ["--finale", str(finale), "--campioni", str(tmp_path / "campioni.tsv")]

    assert vc.main([*base, "--attesi", str(tmp_path / "senza.tsv"), "--uscita", str(tmp_path / "a")]) == 1
    prima = json.loads((tmp_path / "a" / "riepilogo.json").read_text(encoding="utf-8"))
    assert "altri_nomi" not in prima
    assert prima["campioni"]["uno"]["criterio_1"]["mancanti"] == ["Bb"]
    assert prima["campioni"]["uno"]["non_attesi"] == {"Zz": 60, "Cc": 40}

    assert vc.main([*base, "--attesi", str(tmp_path / "con.tsv"), "--uscita", str(tmp_path / "b")]) == 0
    dopo = json.loads((tmp_path / "b" / "riepilogo.json").read_text(encoding="utf-8"))
    assert dopo["altri_nomi"] == {"Zz": "Bb", "Cc": "Bb"}
    uno = dopo["campioni"]["uno"]
    assert uno["criterio_1"] == {"superato": True, "mancanti": []}
    assert uno["quota_non_attesi"] == 0 and uno["non_attesi"] == {}
    generi = {g["genere"]: g for g in uno["generi"]}
    assert generi["Bb"]["letture"] == 100 and generi["Bb"]["misurato"] == pytest.approx(10.0)
    assert dopo["attesi"] == prima["attesi"] and dopo["famiglie_ammesse"] == prima["famiglie_ammesse"]
    capsys.readouterr()

    for altri in ("Zz|Zz", "Aa", "Zz|NA"):
        _tsv(tmp_path / "rotta.tsv", [con[0], con[1], [*senza[2], altri], con[3]])
        assert vc.main([*base, "--attesi", str(tmp_path / "rotta.tsv"), "--uscita", str(tmp_path / "c")]) == 2
        assert "altro nome vuoto o ripetuto" in capsys.readouterr().out
