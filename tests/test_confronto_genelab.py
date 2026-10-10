"""Confronto a livello di genere con le tabelle elaborate da GeneLab.

Verifica ``dati/confronto_genelab/confronta.py`` su tabelle sintetiche piccole,
di cui l'esito si calcola a mano: la normalizzazione dei nomi, i non assegnati
e le etichette non informative, i pareggi nei ranghi e nella classifica, i
doppi zeri, i campioni mancanti da una tabella, il campione con troppo pochi
generi, la correlazione non definita e la quota massima di campioni esclusi.
Lo script misura l'accordo fra due elaborazioni: un errore che lo facesse
apparire migliore del vero non si vedrebbe sui dati reali, dove l'esito atteso
non e' noto.
"""

from __future__ import annotations

import csv
import importlib.util
import json
from pathlib import Path

import pytest

RADICE = Path(__file__).resolve().parents[1]
SCRIPT = RADICE / "dati" / "confronto_genelab" / "confronta.py"

pytestmark = pytest.mark.skipif(not SCRIPT.is_file(), reason="dati/confronto_genelab non presente")


@pytest.fixture(scope="module")
def cg():
    """Il modulo dello script, caricato dal suo file."""
    specifica = importlib.util.spec_from_file_location("confronta_genelab", SCRIPT)
    modulo = importlib.util.module_from_spec(specifica)
    specifica.loader.exec_module(modulo)
    return modulo


#: Dieci generi con abbondanze distinte e uguali nelle due tabelle: il campione
#: di base, con correlazione 1.
BASE = {f"G{i:02d}": 10 * (i + 1) for i in range(10)}


def test_la_normalizzazione_toglie_solo_spazi_e_prefisso_di_rango(cg):
    """
    **Obiettivo**: Verificare che il nome del genere perda gli spazi ai bordi e
    un prefisso di rango in minuscolo, e nient'altro; che vuoto, ``NA`` e le
    etichette non informative (in qualunque combinazione di maiuscole, con lo
    spazio o con il trattino basso) diventino «non assegnato»; e che due
    scritture diverse dello stesso genere restino due nomi diversi.

    **Razionale scientifico e sistemistico**: Ogni equivalenza fra nomi fatta
    dallo script migliorerebbe l'accordo misurato: le sole ammesse sono quelle
    fissate prima di vedere i dati.
    """
    n = cg.genere_normalizzato
    assert n("  g__Bacillus ") == "Bacillus" and n("Bacillus") == "Bacillus"
    assert n("G__Bacillus") == "G__Bacillus" and n("gg__Bacillus") == "gg__Bacillus"
    for vuoto in ("", "  ", None, "NA", " NA ", "g__", "g__NA", "na", "N/A", "NaN", "None",
                  "Unassigned", "g__unassigned"):
        assert n(vuoto) == cg.NON_ASSEGNATO, vuoto
    for etichetta in ("uncultured", "Uncultured bacterium", "g__unidentified", "Unknown Family",
                      "Incertae Sedis", "incertae_sedis", "Clostridiales Incertae Sedis",
                      "gut metagenome", "unclassified_Bacillaceae"):
        assert n(etichetta) == cg.NON_ASSEGNATO, etichetta
    # Solo il valore intero: un genere il cui nome contiene quelle lettere resta.
    assert n("Nana") == "Nana" and n("Nonea") == "Nonea"
    assert n("Escherichia/Shigella") != n("Escherichia-Shigella") and n('"Bacillus"') == '"Bacillus"'
    assert n("Candidatus Nitrosopumilus") == "Candidatus Nitrosopumilus"


def test_spearman_usa_i_ranghi_medi_e_dichiara_la_correlazione_non_definita(cg):
    """
    **Obiettivo**: Verificare la correlazione su una serie con un pareggio,
    calcolata a mano (8 / radice di 95), su due serie in ordine inverso (-1), e
    che con i valori tutti uguali in una serie la correlazione valga 0 e sia
    dichiarata non definita.

    **Razionale scientifico e sistemistico**: Con ranghi minimi invece che
    medi i pareggi, frequenti fra i generi rari, sposterebbero la correlazione;
    e una serie costante non e' accordo, quindi vale 0 e resta nella mediana.
    """
    rho, definita = cg.spearman([1, 2, 3, 4, 5], [5, 6, 7, 8, 7])
    assert definita and rho == pytest.approx(8 / 95 ** 0.5)
    assert cg.spearman([1, 2, 3], [9, 5, 1]) == (pytest.approx(-1.0), True)
    assert cg.spearman([3, 3, 3], [1, 2, 3]) == (0.0, False)
    assert cg.spearman([1, 2, 3], [0.2, 0.2, 0.2]) == (0.0, False)


def test_i_non_assegnati_entrano_nel_totale_ma_non_nella_classifica(cg):
    """
    **Obiettivo**: Verificare che l'abbondanza relativa di un genere sia sul
    totale del campione con i non assegnati, che la quota dei non assegnati sia
    riportata per tabella, e che «non assegnato» non sia fra i generi della
    classifica ne' fra i condivisi.

    **Razionale scientifico e sistemistico**: Togliere i non assegnati dal
    denominatore gonfierebbe le abbondanze della tabella con piu' varianti non
    classificate, e nasconderebbe proprio la differenza fra i classificatori.
    """
    nostra = {"S1": {**BASE, cg.NON_ASSEGNATO: 550}}
    genelab = {"S1": dict(BASE)}
    esito = cg.confronta(nostra, genelab)
    assert esito["quota_non_assegnato"] == {"nostra": pytest.approx(0.5), "genelab": 0.0}
    assert esito["_medie"]["nostra"]["G09"] == pytest.approx(100 / 1100)
    assert esito["_medie"]["genelab"]["G09"] == pytest.approx(100 / 550)
    assert cg.NON_ASSEGNATO not in esito["criterio_1"]["primi_nostra"]
    assert cg.NON_ASSEGNATO not in esito["_condivisi"] and esito["generi"]["condivisi"] == 10
    assert esito["criterio_2"]["quota_letture_nei_generi_condivisi"] == {
        "nostra": pytest.approx(0.5), "genelab": pytest.approx(1.0)}
    assert esito["criterio_2"]["valore"] == pytest.approx(1.0)


def test_la_classifica_dei_primi_dieci_e_il_criterio_degli_otto_in_comune(cg):
    """
    **Obiettivo**: Verificare che i dieci generi piu' abbondanti si scelgano
    per abbondanza relativa media sui campioni comuni e, a pari merito, in
    ordine alfabetico; che con otto in comune il criterio sia superato e con
    sette no.

    **Razionale scientifico e sistemistico**: A pari merito un ordine non
    fissato darebbe un esito diverso da un'esecuzione all'altra, e sceglierlo
    guardando l'altra tabella alzerebbe il numero dei generi in comune.
    """
    # Dodici generi alla pari nella nostra: i primi dieci sono A..J.
    nostra = {"S1": {chr(65 + i): 10 for i in range(12)}}
    # GeneLab: A..H abbondanti (otto in comune), poi due generi solo suoi.
    genelab = {"S1": {**{chr(65 + i): 50 for i in range(8)}, "Y": 40, "Z": 40, "I": 1, "J": 1}}
    esito = cg.confronta(nostra, genelab)
    assert esito["criterio_1"]["primi_nostra"] == list("ABCDEFGHIJ")
    assert esito["criterio_1"]["in_comune"] == list("ABCDEFGH")
    assert esito["criterio_1"]["valore"] == 8 and esito["criterio_1"]["superato"]
    genelab["S1"]["X"] = 45
    del genelab["S1"]["H"]
    esito = cg.confronta(nostra, genelab)
    assert esito["criterio_1"]["valore"] == 7 and not esito["criterio_1"]["superato"]

    # La media e' sui campioni comuni: un genere abbondante solo in un campione
    # che l'altra tabella non ha non entra nella classifica.
    nostra = {"S1": dict(BASE), "S2": {"Raro": 1000}}
    esito = cg.confronta(nostra, {"S1": dict(BASE)})
    assert "Raro" not in esito["criterio_1"]["primi_nostra"]
    assert esito["campioni_comuni"] == ["S1"]


def test_i_doppi_zeri_si_escludono_e_un_genere_mancante_vale_zero(cg):
    """
    **Obiettivo**: Verificare su un campione calcolato a mano che entrino
    nella correlazione i generi condivisi diversi da zero in almeno una delle
    due tabelle (con zero dove mancano), che i generi assenti in entrambe non
    entrino, e che un genere presente in una sola tabella non sia condiviso.

    **Razionale scientifico e sistemistico**: I generi assenti in tutte e due
    le tabelle di un campione sono pareggi a zero che alzano la correlazione
    senza dire nulla sull'accordo; un genere visto da una sola elaborazione va
    invece contato come disaccordo dove l'altra non lo vede.
    """
    # S1: cinque generi; nella nostra Dd ed Ee alla pari, in GeneLab Aa e Bb.
    nostra = {"S1": {"Aa": 50, "Bb": 30, "Cc": 10, "Dd": 5, "Ee": 5},
              "S2": {"Aa": 1, "Bb": 2, "Cc": 3, "Dd": 4, "Ee": 5, "Ff": 6, "Solo": 9}}
    genelab = {"S1": {"Aa": 40, "Bb": 40, "Cc": 10, "Dd": 6, "Ee": 4},
               "S2": {"Aa": 1, "Bb": 2, "Cc": 3, "Dd": 4, "Ee": 5, "Ff": 0, "Altro": 9}}
    esito = cg.confronta(nostra, genelab)
    per_campione = {v["campione"]: v for v in esito["_campioni"]}
    # Ff ha zero in GeneLab ovunque: non e' condiviso. In S1 i generi usati sono cinque.
    assert esito["_condivisi"] == ["Aa", "Bb", "Cc", "Dd", "Ee"]
    assert per_campione["S1"]["generi"] == 5
    assert per_campione["S1"]["correlazione"] == pytest.approx(9 / 9.5)
    # S2: stessi ranghi sui cinque condivisi; Solo, Altro e Ff restano fuori.
    assert per_campione["S2"]["generi"] == 5
    assert per_campione["S2"]["correlazione"] == pytest.approx(1.0)
    assert esito["criterio_2"]["valore"] == pytest.approx((9 / 9.5 + 1) / 2)

    # Un genere condiviso (presente in entrambe in S2) e assente da una tabella
    # in S1 entra in S1 con zero; assente da entrambe in S1 non entra.
    nostra = {"S1": {**BASE, "Kk": 0}, "S2": {**BASE, "Kk": 7, "Ll": 7}}
    genelab = {"S1": {**BASE, "Kk": 500}, "S2": {**BASE, "Kk": 7, "Ll": 7}}
    per_campione = {v["campione"]: v for v in cg.confronta(nostra, genelab)["_campioni"]}
    assert per_campione["S1"]["generi"] == 11 and per_campione["S2"]["generi"] == 12
    # In S1 Kk e' ultimo nella nostra (0) e primo in GeneLab: la correlazione scende.
    assert per_campione["S1"]["correlazione"] == pytest.approx(1 - 6 * (100 + 10) / (11 * 120))


def test_il_campione_con_pochi_generi_e_la_quota_massima_di_esclusi(cg):
    """
    **Obiettivo**: Verificare che un campione con meno di cinque generi usati
    sia escluso dalla mediana ed elencato; che con un escluso su cinque (20%)
    il criterio si giudichi sulla mediana degli altri, e con uno su quattro
    (25%) non sia superato qualunque sia la mediana; e che una correlazione
    non definita resti nella mediana con valore 0.

    **Razionale scientifico e sistemistico**: Una correlazione su tre o
    quattro generi non misura nulla; ma se i campioni scartati sono molti la
    mediana descrive solo i campioni piu' facili, e il criterio non puo' dirsi
    superato.
    """
    pochi = {"G00": 5, "G01": 3, "G02": 2, "G03": 1}
    nostra = {f"S{i}": dict(BASE) for i in range(1, 5)} | {"S5": dict(pochi)}
    genelab = {f"S{i}": dict(BASE) for i in range(1, 5)} | {"S5": dict(pochi)}
    esito = cg.confronta(nostra, genelab)
    c2 = esito["criterio_2"]
    assert c2["esclusi"] == ["S5"] and not c2["esclusi_oltre_il_massimo"]
    assert c2["campioni_nella_mediana"] == 4 and c2["valore"] == pytest.approx(1.0) and c2["superato"]

    del nostra["S4"], genelab["S4"]
    c2 = cg.confronta(nostra, genelab)["criterio_2"]
    assert c2["esclusi"] == ["S5"] and c2["esclusi_oltre_il_massimo"]
    assert c2["valore"] == pytest.approx(1.0) and not c2["superato"]

    # Correlazione non definita: nella nostra i generi usati hanno tutti lo
    # stesso conteggio. Vale 0 e resta nella mediana: con tre campioni a 1 e
    # due a 0 la mediana e' 1, con due e tre e' 0.
    piatto = dict.fromkeys(BASE, 7)
    nostra = {"S1": dict(BASE), "S2": dict(BASE), "S3": dict(BASE), "S4": piatto, "S5": piatto}
    genelab = {f"S{i}": dict(BASE) for i in range(1, 6)}
    esito = cg.confronta(nostra, genelab)
    stati = {v["campione"]: (v["correlazione"], v["stato"]) for v in esito["_campioni"]}
    assert stati["S4"] == (0.0, "nella mediana: non definita, vale 0")
    assert esito["criterio_2"]["campioni_nella_mediana"] == 5
    assert esito["criterio_2"]["valore"] == pytest.approx(1.0)
    nostra["S3"] = piatto
    c2 = cg.confronta(nostra, genelab)["criterio_2"]
    assert c2["valore"] == 0.0 and not c2["superato"]
    assert c2["distribuzione"] == {"minimo": 0.0, "q1": 0.0, "mediana": 0.0, "q3": 1.0, "massimo": 1.0}


def test_la_soglia_della_mediana_e_inclusa(cg, monkeypatch):
    """
    **Obiettivo**: Verificare che una mediana pari alla soglia superi il
    criterio 2 e una appena inferiore no.

    **Razionale scientifico e sistemistico**: «Almeno» comprende il valore
    della soglia: il confine va fissato dal codice come dal testo del criterio.
    """
    nostra = {"S1": {"Aa": 50, "Bb": 30, "Cc": 10, "Dd": 5, "Ee": 5}}
    genelab = {"S1": {"Aa": 40, "Bb": 40, "Cc": 10, "Dd": 6, "Ee": 4}}
    rho = cg.confronta(nostra, genelab)["criterio_2"]["valore"]
    monkeypatch.setattr(cg, "MEDIANA_MINIMA", rho)
    assert cg.confronta(nostra, genelab)["criterio_2"]["superato"]
    monkeypatch.setattr(cg, "MEDIANA_MINIMA", rho + 1e-9)
    assert not cg.confronta(nostra, genelab)["criterio_2"]["superato"]


def _tsv(percorso: Path, righe: list[list]) -> None:
    with open(percorso, "w", encoding="utf-8", newline="") as file:
        csv.writer(file, delimiter="\t", lineterminator="\n").writerows(righe)


def test_dai_file_i_campioni_si_abbinano_per_nome_e_i_controlli_restano_fuori(cg, tmp_path, capsys):
    """
    **Obiettivo**: Verificare dalla riga di comando, su file scritti a mano,
    che le varianti si sommino per genere, che i campioni si abbinino con il
    nome dei metadati e non con la posizione (colonne in ordine diverso nelle
    due tabelle), che le colonne di GeneLab senza un campione nostro (un
    controllo, un campione escluso dai filtri) restino fuori e siano elencate
    con il motivo, che l'elenco delle corse da escludere tolga il campione, e
    che uscite e codice di uscita riportino l'esito.

    **Razionale scientifico e sistemistico**: Un abbinamento per posizione
    confronterebbe campioni diversi senza alcun errore; e un controllo entrato
    nel confronto vi porterebbe una composizione che la pipeline, per scelta,
    non consegna.
    """
    finale = tmp_path / "12_final"
    (finale / "intermedi").mkdir(parents=True)
    generi = sorted(BASE)
    # Due varianti per genere nella nostra tabella: la somma e' quella di BASE.
    conteggi = [["asv_id", "ACC2", "ACC1", "ACC3"]]
    tassonomia = [["asv_id", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus"]]
    for i, g in enumerate(generi):
        meta = BASE[g] // 2
        # Nel terzo campione i generi sono in ordine inverso.
        inverso = BASE[generi[len(generi) - 1 - i]]
        conteggi += [[f"ASV{2 * i + 1}", meta, meta, inverso], [f"ASV{2 * i + 2}", BASE[g] - meta, BASE[g] - meta, 0]]
        tassonomia += [[f"ASV{2 * i + 1}", "Bacteria", "P", "C", "O", "F", g],
                       [f"ASV{2 * i + 2}", "Bacteria", "P", "C", "O", "F", g]]
    conteggi.append(["ASV99", 5, 5, 5])
    tassonomia.append(["ASV99", "Bacteria", "P", "C", "O", "F", ""])
    _tsv(finale / "conteggi.tsv", conteggi)
    _tsv(finale / "tassonomia.tsv", tassonomia)
    _tsv(finale / "metadati.tsv", [["accession", "sample_name", "classe"],
                                   ["ACC1", "uno", "biologico"], ["ACC2", "due", "biologico"],
                                   ["ACC3", "tre", "biologico"]])
    _tsv(finale / "intermedi" / "esclusioni.tsv",
         [["accession", "sample_name", "piastra", "filtro", "motivo"],
          ["ACC4", "quattro", "", "letture_finali", "troppe poche letture"]])
    # GeneLab: colonne in un altro ordine, prefisso di rango, un controllo, il
    # campione escluso dai nostri filtri, e il campione "tre" con i generi nello
    # stesso ordine di BASE (opposto al nostro).
    genelab = [["ASV_ID", "domain", "phylum", "class", "order", "family", "genus", "species",
                "x_tre", "x_bianco", "x_uno", "x_quattro", "x_due"]]
    for i, g in enumerate(generi):
        genelab.append([f"ASV_{i}", "Bacteria", "P", "C", "O", "F", f"g__{g}", "NA",
                        BASE[g], 3, BASE[g], BASE[g], BASE[g]])
    genelab.append(["ASV_x", "Bacteria", "P", "C", "O", "F", "NA", "NA", 0, 9, 0, 0, 0])
    _tsv(tmp_path / "genelab.tsv", genelab)
    _tsv(tmp_path / "abbinamento.tsv", [["Sample Name", "identificativo"]] + [
        [n, f"x_{n}"] for n in ("uno", "due", "tre", "quattro", "bianco")])
    base = ["--finale", str(finale), "--genelab", str(tmp_path / "genelab.tsv"),
            "--abbinamento", str(tmp_path / "abbinamento.tsv")]

    # Tutti e tre: due campioni a 1 e uno a -1, mediana 1; criterio 1 con dieci in comune.
    assert cg.main([*base, "--uscita", str(tmp_path / "a")]) == 0
    stampato = capsys.readouterr().out
    assert "criterio 1" in stampato and "criterio 2" in stampato and "NON superato" not in stampato
    riepilogo = json.loads((tmp_path / "a" / "riepilogo.json").read_text(encoding="utf-8"))
    assert riepilogo["campioni_comuni"] == ["due", "tre", "uno"]
    assert riepilogo["criterio_1"]["valore"] == 10 and riepilogo["criterio_2"]["valore"] == 1.0
    assert riepilogo["solo_nella_tabella_di_genelab"] == {
        "bianco": "non e' fra i campioni biologici dell'oggetto finale",
        "quattro": "letture_finali: troppe poche letture"}
    assert riepilogo["quota_non_assegnato"]["nostra"] == pytest.approx(15 / (3 * 555))
    assert riepilogo["quota_non_assegnato"]["genelab"] == 0
    with open(tmp_path / "a" / "campioni.tsv", encoding="utf-8", newline="") as file:
        campioni = {r["campione"]: r for r in csv.DictReader(file, delimiter="\t")}
    assert float(campioni["uno"]["correlazione"]) == 1.0
    assert float(campioni["tre"]["correlazione"]) == -1.0
    assert campioni["uno"]["letture_nostra"] == "555" and campioni["uno"]["letture_genelab"] == "550"

    # Escludendo "uno" restano 1 e -1: mediana 0, criterio 2 non superato.
    (tmp_path / "escludi.txt").write_text("# corse scambiate\nuno\n", encoding="utf-8")
    assert cg.main([*base, "--escludi", str(tmp_path / "escludi.txt"),
                    "--uscita", str(tmp_path / "b")]) == 1
    assert "NON superato" in capsys.readouterr().out
    riepilogo = json.loads((tmp_path / "b" / "riepilogo.json").read_text(encoding="utf-8"))
    assert riepilogo["campioni_comuni"] == ["due", "tre"] and riepilogo["esclusi_dal_confronto"] == ["uno"]
    assert riepilogo["criterio_2"]["valore"] == 0.0 and not riepilogo["criterio_2"]["superato"]

    # Un identificativo che la tabella non ha, un campione da escludere che
    # l'abbinamento non conosce, una riga con un campo in meno, una virgoletta
    # spaiata: il confronto non si esegue, invece di perdere dati in silenzio.
    _tsv(tmp_path / "refuso.tsv", [["Sample Name", "identificativo"], ["uno", "x_uno"], ["due", "x_duo"]])
    assert cg.main(["--finale", str(finale), "--genelab", str(tmp_path / "genelab.tsv"),
                    "--abbinamento", str(tmp_path / "refuso.tsv"), "--uscita", str(tmp_path / "d")]) == 2
    assert "x_duo" in capsys.readouterr().out
    (tmp_path / "ignoto.txt").write_text("mai visto\n", encoding="utf-8")
    assert cg.main([*base, "--escludi", str(tmp_path / "ignoto.txt"), "--uscita", str(tmp_path / "e")]) == 2
    assert "mai visto" in capsys.readouterr().out
    testo = (tmp_path / "genelab.tsv").read_text(encoding="utf-8").splitlines()
    (tmp_path / "corta.tsv").write_text("\n".join([*testo[:3], testo[3].rsplit("\t", 1)[0], *testo[4:]]) + "\n",
                                        encoding="utf-8")
    assert cg.main(["--finale", str(finale), "--genelab", str(tmp_path / "corta.tsv"),
                    "--abbinamento", str(tmp_path / "abbinamento.tsv"), "--uscita", str(tmp_path / "f")]) == 2
    assert "riga 4" in capsys.readouterr().out
    # Con le virgolette lette come testo la riga successiva non sparisce: il
    # genere fra virgolette e' un altro nome, e i conteggi restano tutti.
    virgolette = [r.replace("g__G03", '"g__G03') if i == 4 else r for i, r in enumerate(testo)]
    (tmp_path / "virgolette.tsv").write_text("\n".join(virgolette) + "\n", encoding="utf-8")
    assert cg.main(["--finale", str(finale), "--genelab", str(tmp_path / "virgolette.tsv"),
                    "--abbinamento", str(tmp_path / "abbinamento.tsv"), "--uscita", str(tmp_path / "g")]) in (0, 1)
    capsys.readouterr()
    riepilogo = json.loads((tmp_path / "g" / "riepilogo.json").read_text(encoding="utf-8"))
    assert riepilogo["letture"]["genelab"] == 3 * 550 and riepilogo["generi"]["genelab"] == 10
    # Un campione nostro senza colonna abbinata e' detto sulla console.
    _tsv(tmp_path / "parziale.tsv", [["Sample Name", "identificativo"], ["uno", "x_uno"], ["due", "x_due"]])
    cg.main(["--finale", str(finale), "--genelab", str(tmp_path / "genelab.tsv"),
             "--abbinamento", str(tmp_path / "parziale.tsv"), "--uscita", str(tmp_path / "h")])
    assert "fuori dal confronto: tre" in capsys.readouterr().out

    # Senza campioni in comune il confronto non si esegue.
    _tsv(tmp_path / "vuoto.tsv", [["Sample Name", "identificativo"], ["altro", "x_bianco"]])
    assert cg.main(["--finale", str(finale), "--genelab", str(tmp_path / "genelab.tsv"),
                    "--abbinamento", str(tmp_path / "vuoto.tsv"), "--uscita", str(tmp_path / "c")]) == 2
    assert "confronto non eseguito" in capsys.readouterr().out
