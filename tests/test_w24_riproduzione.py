r"""Suite di test della settimana 24: la riproduzione dei risultati su OSD-734.

1. Inquadramento nel Piano Operativo
------------------------------------
Settimana 24 (W24), Fase F7, milestone M7 (riproduzione da zero dei risultati
dal repository pubblicato e dalle fonti pubbliche): i checksum attesi
pubblicati e il loro confronto. Anticipa attivita' della fase F9 sulla
riproducibilita'.

2. Moduli sorgente coperti
--------------------------
* ``dati/osd734/confronta_risultati.py``
* ``dati/osd734/checksum_finali.sha256``, ``dati/osd734/checksum_artefatti.tsv``
* ``src/amplicon16s/runner/graph.py`` (l'ordine di esecuzione delle fasi, che
  lo script riporta)
* ``src/amplicon16s/cli.py`` (percorsi relativi della configurazione resi
  assoluti al caricamento)

3. Cosa valuta questo file
--------------------------
- una configurazione con percorsi relativi, come ``config_osd734.yaml``, si
  esegue dalla cartella a cui i percorsi si riferiscono: la configurazione
  registrata porta i percorsi assoluti, i parametri dichiarati restano gli
  stessi, e una fase con un calcolo in R trova i suoi ingressi;
- l'ordine delle fasi dello script e' quello del grafo, con S9 dopo S13;
- ``--aggiorna`` scrive i due file dei checksum attesi da un'esecuzione, e il
  confronto della stessa esecuzione li trova identici, con esito 0;
- un artefatto alterato, uno mancante e uno non atteso danno esito 1, e il
  messaggio indica la prima fase, nell'ordine di esecuzione, in cui la
  differenza compare; manifesti, log e configurazione registrata non contano;
- i checksum pubblicati sono coerenti: i file consegnati di ``12_final/`` sono
  quelli di S14 nell'elenco di tutti gli artefatti, con la stessa impronta;
  ogni fase eseguita con la filogenesi disattivata vi compare; il README della
  cartella riporta il comando del confronto;
- sul dataset completo: la catena eseguita dalla suite da' gli artefatti
  attesi pubblicati.

4. Comandi Bash e scenari di esecuzione
---------------------------------------
    ``<immagine>`` e' l'immagine del container della pipeline; quella corrente
    e' indicata in ``test.txt``, sezione 1.3.

    1. Modalità locale standard (R di base con jsonlite, senza Bioconductor né
       dati reali):
       pytest tests/test_w24_riproduzione.py -v

    2. Modalità container Docker standard (sottoinsieme ridotto con
       R/Bioconductor):
       docker run --rm \
         -e PYTHONPATH=/app/src \
         -e AMPLICON16S_R_DIR=/app/R \
         -v "$(pwd)":/app \
         -w /app \
         <immagine> \
         pytest -o cache_dir=/tmp/.pytest_cache tests/test_w24_riproduzione.py -v

    3. Modalità container Docker completa (con i dati reali OSD-734):
       docker run --rm \
         --memory=24g \
         --memory-swap=24g \
         -e PYTHONPATH=/app/src \
         -e AMPLICON16S_R_DIR=/app/R \
         -e AMPLICON16S_CONFIG_DATI_REALI="$HOME/ASI/config_osd734.yaml" \
         -v "$(pwd)":/app \
         -v "$HOME/ASI":"$HOME/ASI" \
         -w /app \
         <immagine> \
         pytest -o cache_dir=/tmp/.pytest_cache tests/test_w24_riproduzione.py -v

5. Risultato atteso
-------------------
Vedi ``test.txt``, scheda W24. Nell'immagine senza il repository montato, come
nella CI, ``dati/`` manca: i cinque test che usano lo script di confronto e i
checksum pubblicati si saltano, i due sui percorsi relativi girano.

6. Razionale scientifico e sistemistico
---------------------------------------
- La riproducibilita' si prova confrontando byte, non impressioni: chi riproduce
  deve poter verificare con un comando di aver ottenuto gli stessi risultati.
- Quando un risultato finale non coincide, sapere qual e' la prima fase
  divergente restringe la ricerca della causa da quindici fasi a una.
- I checksum attesi sono essi stessi un dato pubblicato: un elenco incoerente
  farebbe fallire una riproduzione corretta.
- La prima riproduzione da un clone si e' fermata in S1 perche' i percorsi
  relativi della configurazione pubblicata non valevano nella cartella in cui
  parte il processo R: il caso e' ora un test.
"""

from __future__ import annotations

import csv
import hashlib
import json
import os
import sys
from pathlib import Path

import pytest
import yaml
from conftest import NEGATIVO, POSITIVO, Campione, crea_scenario, dichiarazione_minima
from sottoinsieme import dati_config, motivo_pacchetti_r_assenti

RADICE = Path(__file__).resolve().parents[1]
DATI = RADICE / "dati" / "osd734"
#: Nell'immagine la cartella dati/ non viene copiata: i test che usano lo script
#: di confronto e i checksum pubblicati si saltano uno per uno, mentre quelli
#: sui percorsi relativi della configurazione girano anche li', e quindi nella CI.
HA_DATI = DATI.is_dir()
richiede_dati = pytest.mark.skipif(
    not HA_DATI, reason="la cartella dati/ non e' presente (nell'immagine non viene copiata)"
)
if HA_DATI:
    sys.path.insert(0, str(DATI))
    import confronta_risultati as confronto  # noqa: E402

import amplicon16s.cli as cli  # noqa: E402
from amplicon16s.config import defaults  # noqa: E402
from amplicon16s.config.schema import carica  # noqa: E402
from amplicon16s.io_layer.artifacts import Fase  # noqa: E402
from amplicon16s.logging.logger import chiudi  # noqa: E402
from amplicon16s.runner.graph import GRAFO, Passo  # noqa: E402
from amplicon16s.runner.project import passi_realizzati  # noqa: E402

#: Tre fasi di prova, con due cartelle: S6 e S7 condividono 07_chimera.
FASI = {
    "S5": ("06_seqtab", {"tabella.rds": b"tabella", "tabella.json": b"{}"}),
    "S6": ("07_chimera", {"chimere.tsv": b"chimere"}),
    "S13": ("12_final", {"ps_filtrato.rds": b"filtrato"}),
    "S14": ("12_final", {"ps_final.rds": b"finale", "conteggi.tsv": b"conteggi"}),
}


def _esecuzione(radice: Path) -> Path:
    """Una cartella di output di prova: artefatti, manifesti di fase, il manifesto
    dei file consegnati, e cio' che il confronto deve ignorare (log e
    configurazione registrata).
    """
    for fase, (cartella, artefatti) in FASI.items():
        (radice / cartella).mkdir(parents=True, exist_ok=True)
        voci = []
        for nome, contenuto in artefatti.items():
            (radice / cartella / nome).write_bytes(contenuto)
            voci.append({"nome": nome, "byte": len(contenuto),
                         "checksum": "sha256:" + hashlib.sha256(contenuto).hexdigest()})
        (radice / cartella / f"manifest_{fase}.json").write_text(json.dumps(
            {"passo": fase, "cartella": cartella, "artefatti": voci, "esecuzione": radice.name}
        ), encoding="utf-8")
    consegnati = FASI["S14"][1]
    (radice / "12_final" / "checksum.sha256").write_text("".join(
        f"{hashlib.sha256(c).hexdigest()}  {n}\n" for n, c in consegnati.items()
    ), encoding="utf-8")
    (radice / "99_logs").mkdir()
    (radice / "99_logs" / "pipeline.jsonl").write_text(f"log di {radice.name}\n", encoding="utf-8")
    (radice / "00_config").mkdir()
    (radice / "00_config" / "resolved.yaml").write_text(f"out_root: {radice}\n", encoding="utf-8")
    return radice


@pytest.fixture
def attesi(tmp_path, monkeypatch) -> Path:
    """I file dei checksum attesi in una cartella temporanea, scritti da
    un'esecuzione di prova: quelli pubblicati non vengono toccati.
    """
    monkeypatch.setattr(confronto, "ATTESI_FINALI", tmp_path / "checksum_finali.sha256")
    monkeypatch.setattr(confronto, "ATTESI_ARTEFATTI", tmp_path / "checksum_artefatti.tsv")
    assert confronto.main(["--uscita", str(_esecuzione(tmp_path / "riferimento")), "--aggiorna"]) == 0
    return tmp_path


_PERCORSI = (
    ("io", "fastq_dir"), ("io", "assay_table"), ("io", "study_table"), ("io", "out_root"),
    ("io", "batch_table"), ("tax", "ref_fasta"), ("tax", "ref_bad_taxa"),
)


def _config_relativa(dati: dict, cartella: Path) -> Path:
    """Scrive in ``cartella`` una configurazione i cui percorsi sono relativi a
    ``cartella`` stessa, come quelli di ``config_osd734.yaml`` lo sono alla
    radice del repository.
    """
    for gruppo, chiave in _PERCORSI:
        if dati[gruppo].get(chiave) is not None:
            dati[gruppo][chiave] = os.path.relpath(dati[gruppo][chiave], cartella)
            assert not Path(dati[gruppo][chiave]).is_absolute()
    percorso = cartella / "config.yaml"
    percorso.write_text(yaml.safe_dump(dati, sort_keys=False), encoding="utf-8")
    return percorso


@pytest.fixture(autouse=True)
def uscite_pulite():
    """Chiude le uscite del log prima e dopo ogni test, perché nessun handler resti
    aperto sulla cartella temporanea.
    """
    chiudi()
    yield
    chiudi()


def test_i_percorsi_relativi_valgono_dalla_cartella_di_lancio(tmp_path, monkeypatch, capsys):
    """
    **Obiettivo**: Verificare che una configurazione con soli percorsi
    relativi, lanciata dalla cartella a cui si riferiscono, concluda S0; che la
    configurazione registrata riporti i percorsi assoluti dei file letti e gli
    stessi parametri dichiarati del file; e che i percorsi gia' assoluti
    restino com'erano.

    **Razionale scientifico e sistemistico**: La configurazione registrata deve
    dire quali file sono stati letti: un percorso relativo lo direbbe solo a
    chi conosce la cartella da cui il comando e' partito.
    """
    scenario = crea_scenario(tmp_path, [
        Campione("ERX3000001", "NOD1D4.L1"), Campione("ERX3000002", "NOD1D4.L2"),
        Campione("ERX3000003", "POS.P1.1", materiale=POSITIVO, posizione="Not Applicable"),
        Campione("ERX3000004", "BLANK.P1.1", materiale=NEGATIVO, posizione="Not Applicable"),
    ], con_letture=True)
    config = scenario.config
    dati = dichiarazione_minima(config)
    percorso = _config_relativa(dati, tmp_path)
    monkeypatch.chdir(tmp_path)
    assert cli.main(["validate", "--config", "config.yaml"]) == 0
    capsys.readouterr()

    registrata = yaml.safe_load((tmp_path / "out" / Fase.CONFIG.value / "resolved.yaml").read_text())
    for gruppo, chiave in _PERCORSI:
        valore = registrata["parametri"][gruppo][chiave]
        if valore is not None:
            assert valore == str(getattr(getattr(config, gruppo), chiave)), chiave
    assert set(registrata["dichiarati"]) == {
        f"{gruppo}.{nome}" for gruppo, valori in dati.items() for nome in valori
    }
    assert set(defaults.OBBLIGATORI) <= set(registrata["dichiarati"])
    assoluta = cli._percorsi_assoluti(config)
    assert assoluta is config and carica(percorso).io.out_root == Path("out")


def test_una_fase_in_r_trova_gli_ingressi_con_percorsi_relativi(tmp_path, monkeypatch, capsys):
    """
    **Obiettivo**: Verificare, sul sottoinsieme di prova, che con una
    configurazione a percorsi relativi lanciata dalla cartella a cui si
    riferiscono si concludano S0 e S1, la prima fase con un calcolo in R.

    **Razionale scientifico e sistemistico**: Il processo R di una fase parte
    nella cartella della fase: con percorsi relativi non trovava la propria
    richiesta, e la riproduzione dal repository pubblicato si fermava in S1 con
    ``E-R-02``.
    """
    motivo = motivo_pacchetti_r_assenti("dada2", "ggplot2", "ShortRead", "Biostrings", "jsonlite")
    if motivo is not None:
        if os.environ.get("AMPLICON16S_RICHIEDI_BIOC") == "1":
            pytest.fail(f"Bioconductor e' richiesto in questo ambiente: {motivo}")
        pytest.skip(motivo)
    _config_relativa(dati_config(tmp_path), tmp_path)
    vere = passi_realizzati()
    monkeypatch.setattr(cli, "_passi", lambda: {p: vere[p] for p in (Passo.S0, Passo.S1)})
    monkeypatch.chdir(tmp_path)
    # Con le sole S0 e S1 registrate l'esecuzione si ferma alla prima fase senza
    # codice, dopo aver concluso le due.
    assert cli.main(["run", "--config", "config.yaml"]) == cli.USCITA_FASE_NON_REALIZZATA
    capsys.readouterr()
    assert (tmp_path / "out" / Fase.QC_PROFILES.value / "manifest_S1.json").is_file()
    assert (tmp_path / "out" / Fase.QC_PROFILES.value / "letture_grezze.tsv").is_file()


def _pubblicati() -> list[dict[str, str]]:
    """Le righe di ``checksum_artefatti.tsv`` pubblicato."""
    with open(DATI / "checksum_artefatti.tsv", encoding="utf-8", newline="") as file:
        return list(csv.DictReader(file, delimiter="\t"))


@richiede_dati
def test_l_ordine_delle_fasi_e_quello_del_grafo():
    """
    **Obiettivo**: Verificare che l'ordine delle fasi dichiarato nello script di
    confronto, che non importa il pacchetto della pipeline, coincida con
    l'ordine di esecuzione del grafo.

    **Razionale scientifico e sistemistico**: La prima fase divergente e' la
    prima nell'ordine di esecuzione: con S9 dopo S13 l'ordine dei numeri
    indicherebbe la fase sbagliata.
    """
    assert confronto.ORDINE == tuple(str(p) for p in GRAFO.ordine())


@richiede_dati
def test_la_stessa_esecuzione_e_identica_agli_attesi(attesi, capsys):
    """
    **Obiettivo**: Verificare che ``--aggiorna`` scriva l'elenco di tutti gli
    artefatti dei manifesti, nell'ordine di esecuzione, e i file consegnati di
    ``12_final/``; e che una seconda esecuzione con gli stessi artefatti, ma
    con manifesti, log e configurazione registrata diversi, risulti identica
    con esito 0.

    **Razionale scientifico e sistemistico**: Manifesti, log e configurazione
    registrata portano date, identificativi e percorsi: differiscono per
    costruzione, e contarli farebbe fallire ogni riproduzione.
    """
    with open(attesi / "checksum_artefatti.tsv", encoding="utf-8", newline="") as file:
        righe = list(csv.DictReader(file, delimiter="\t"))
    assert [(r["fase"], r["nome"]) for r in righe] == [
        ("S5", "tabella.json"), ("S5", "tabella.rds"), ("S6", "chimere.tsv"),
        ("S13", "ps_filtrato.rds"), ("S14", "conteggi.tsv"), ("S14", "ps_final.rds"),
    ]
    assert righe[1]["sha256"] == hashlib.sha256(b"tabella").hexdigest() and righe[1]["byte"] == "7"
    finali = (attesi / "checksum_finali.sha256").read_text(encoding="utf-8")
    assert finali == f"{hashlib.sha256(b'finale').hexdigest()}  ps_final.rds\n" \
                     f"{hashlib.sha256(b'conteggi').hexdigest()}  conteggi.tsv\n"

    capsys.readouterr()
    assert confronto.main(["--uscita", str(_esecuzione(attesi / "nuova"))]) == 0
    uscita = capsys.readouterr().out
    assert "RISULTATI IDENTICI: 6 artefatti di 4 fasi" in uscita
    assert uscita.count("identica") == 4


@richiede_dati
def test_una_differenza_indica_la_prima_fase_divergente(attesi, capsys):
    """
    **Obiettivo**: Verificare che un artefatto alterato in S6 e uno in S14, un
    artefatto mancante in S13 e uno non atteso in S5 diano esito 1, ciascuno
    riportato sotto la sua fase, e che la prima fase divergente sia indicata
    nell'ordine di esecuzione; e che una cartella inesistente dia esito 1.

    **Razionale scientifico e sistemistico**: Una differenza a monte si propaga
    a tutte le fasi che ne dipendono: la causa sta nella prima, e li' va
    cercata.
    """
    nuova = _esecuzione(attesi / "nuova")
    (nuova / "07_chimera" / "chimere.tsv").write_bytes(b"altre chimere")
    (nuova / "12_final" / "ps_final.rds").write_bytes(b"altro finale")
    capsys.readouterr()
    assert confronto.main(["--uscita", str(nuova)]) == 1
    uscita = capsys.readouterr().out
    assert "S5      2 artefatti attesi: identica" in uscita
    assert "chimere.tsv: contenuto diverso" in uscita and "ps_final.rds: contenuto diverso" in uscita
    assert "La prima fase in cui compare una differenza e' S6" in uscita

    (nuova / "07_chimera" / "chimere.tsv").write_bytes(b"chimere")
    (nuova / "12_final" / "ps_filtrato.rds").unlink()
    assert confronto.main(["--uscita", str(nuova)]) == 1
    uscita = capsys.readouterr().out
    assert "ps_filtrato.rds: manca" in uscita and "differenza e' S13" in uscita

    manifesto = nuova / "06_seqtab" / "manifest_S5.json"
    documento = json.loads(manifesto.read_text(encoding="utf-8"))
    documento["artefatti"].append({"nome": "in_piu.tsv", "byte": 1, "checksum": "sha256:0"})
    manifesto.write_text(json.dumps(documento), encoding="utf-8")
    assert confronto.main(["--uscita", str(nuova)]) == 1
    uscita = capsys.readouterr().out
    assert "in_piu.tsv: non atteso" in uscita and "differenza e' S5" in uscita

    assert confronto.main(["--uscita", str(attesi / "assente")]) == 1
    assert "non esiste" in capsys.readouterr().out


@richiede_dati
def test_i_checksum_pubblicati_sono_coerenti():
    """
    **Obiettivo**: Verificare che i checksum pubblicati dei file consegnati
    siano, nell'elenco di tutti gli artefatti, quelli di S14 in ``12_final``
    con la stessa impronta; che l'elenco copra tutte le fasi tranne S9,
    disattivata per difetto, senza artefatti ripetuti e con impronte ben
    formate; che vi compaia l'oggetto finale; e che il README della cartella
    riporti il comando del confronto e i due file.

    **Razionale scientifico e sistemistico**: I due file descrivono la stessa
    esecuzione: se divergessero, il confronto dei risultati finali e quello
    fase per fase darebbero esiti contraddittori.
    """
    righe = _pubblicati()
    assert {r["fase"] for r in righe} == set(confronto.ORDINE) - {"S9"}
    assert len({(r["cartella"], r["nome"]) for r in righe}) == len(righe)
    assert all(len(r["sha256"]) == 64 and int(r["byte"]) >= 0 for r in righe)
    assert [r["fase"] for r in righe] == sorted(
        (r["fase"] for r in righe), key=confronto.ORDINE.index
    )
    per_nome = {r["nome"]: r for r in righe if r["fase"] == "S14"}
    finali = (DATI / "checksum_finali.sha256").read_text(encoding="utf-8").splitlines()
    assert [riga.split()[1] for riga in finali] == [
        "ps_final.rds", "conteggi.tsv", "tassonomia.tsv", "metadati.tsv", "sequenze.fasta",
    ]
    for riga in finali:
        impronta, nome = riga.split()
        assert per_nome[nome]["sha256"] == impronta and per_nome[nome]["cartella"] == "12_final"
    leggimi = (DATI / "README.md").read_text(encoding="utf-8")
    for atteso in ("python3 dati/osd734/confronta_risultati.py", "checksum_finali.sha256",
                   "checksum_artefatti.tsv", "Se i risultati non coincidono"):
        assert atteso in leggimi, atteso


@richiede_dati
@pytest.mark.dati_reali
def test_la_catena_sul_dataset_completo_da_i_checksum_pubblicati(catena_reale, capsys):
    """
    **Obiettivo**: Verificare che gli artefatti della catena S0-S14 eseguita
    dalla suite sul dataset completo coincidano con i checksum attesi
    pubblicati, tranne quelli di S1, che la catena condivisa non esegue, e le
    differenze per costruzione dichiarate in ``dati/osd734/README.md`` finche'
    i checksum pubblicati non si rigenerano: l'ordine dei gate in S0, gli
    intermedi di S13 spostati in ``12_final/intermedi``, e in S14 l'elenco dei
    checksum con ``ps_controlli.rds``. I cinque file consegnati pubblicati
    devono coincidere.

    **Razionale scientifico e sistemistico**: E' la regressione piu' stretta
    che la suite possa esprimere: lo stesso dato, la stessa configurazione e il
    codice corrente devono dare gli stessi byte dell'esecuzione di riferimento.
    """
    run, _ = catena_reale
    capsys.readouterr()
    esito = confronto.main(["--uscita", str(run.config.io.out_root)])
    uscita = capsys.readouterr().out
    diverse = [riga.split()[0] for riga in uscita.splitlines() if "differenze" in riga]
    assert set(diverse) <= {"S0", "S1", "S13", "S14"}, uscita
    assert esito == (1 if diverse else 0)
    # Ogni riga di differenza dev'essere una di quelle dichiarate, per nome:
    # nient'altro puo' cambiare, mancare o comparire.
    pubblicati = _pubblicati()
    di_s1 = {r["nome"] for r in pubblicati if r["fase"] == "S1"}
    di_s13 = {r["nome"]: r["sha256"] for r in pubblicati if r["fase"] == "S13"}
    spostati = set(di_s13) - {"ps_controlli.rds"}
    ammesse = (
        {("S0", "gates.json", "contenuto diverso"), ("S14", "checksum.sha256", "contenuto diverso"),
         ("S14", "ps_controlli.rds", "non atteso"), ("S1", "valori_qualita.tsv", "non atteso")}
        | {("S1", nome, "manca") for nome in di_s1}
        | {("S13", nome, "manca") for nome in spostati}
    )
    fase = None
    for riga in uscita.splitlines():
        if riga[:1] == "S" and "artefatti attesi" in riga:
            fase = riga.split()[0]
        elif riga.startswith(" ") and ": " in riga:
            nome, stato = (parte.strip() for parte in riga.rsplit(": ", 1))
            assert (fase, nome, stato) in ammesse, riga
    # Gli intermedi di S13 hanno gli stessi byte, nella loro cartella nuova; i
    # file consegnati pubblicati coincidono, e sono le prime righe dell'elenco.
    finale = Path(run.config.io.out_root) / "12_final"
    for nome in spostati:
        assert confronto.sha256(finale / "intermedi" / nome) == di_s13[nome], nome
    assert confronto.sha256(finale / "ps_controlli.rds") == di_s13["ps_controlli.rds"]
    attesi_finali = confronto.ATTESI_FINALI.read_text(encoding="utf-8").splitlines()
    for riga in attesi_finali:
        impronta, nome = riga.split(maxsplit=1)
        assert confronto.sha256(finale / nome) == impronta, nome
    elenco = (finale / "checksum.sha256").read_text(encoding="utf-8").splitlines()
    assert elenco[:len(attesi_finali)] == attesi_finali and len(elenco) == len(attesi_finali) + 1
    attesi = sum(r["fase"] != "S1" for r in _pubblicati())
    assert uscita.count("identica") == len(confronto.ORDINE) - 1 - len(diverse) and attesi > 1000
