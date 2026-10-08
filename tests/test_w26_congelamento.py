r"""Suite di test della settimana 26: configurazione congelata e regola rigorosa.

1. Inquadramento nel Piano Operativo
------------------------------------
Settimana 26 (W26), Fase F7 (congelamento della configurazione di OSD-734,
regola rigorosa sulla provenienza).

2. Moduli sorgente coperti
--------------------------
* ``src/amplicon16s/runner/provenienza.py`` (``verifica_git``, ``blocco_r``,
  ``discordanze_ambiente``, ``impronta_rigorosa``)
* ``src/amplicon16s/runner/executor.py`` (controlli di avvio con la regola
  rigorosa), ``src/amplicon16s/runner/project.py``, ``src/amplicon16s/steps/base.py``
* ``R/00_ambiente.R``
* ``src/amplicon16s/config/schema.py``, ``src/amplicon16s/config/defaults.py``
  (``run.strict_provenance``), ``config/config.example.yaml``,
  ``dati/osd734/config_osd734.yaml``
* ``src/amplicon16s/report/builder.py`` (colonna della regola rigorosa)

3. Cosa valuta questo file
--------------------------
- ``run.strict_provenance`` falso per difetto e vero nella configurazione
  congelata; senza la regola l'impronta di fase non cambia, con la regola
  comprende sorgente, file di blocco e immagine dichiarata;
- rifiuto dell'esecuzione, prima di qualunque fase, con git non leggibile
  (``E-PROV-01``), con script R fuori dal repository (``E-PROV-01``), con
  modifiche non committate (``E-PROV-02``), con ambiente R diverso dal file di
  blocco (``E-PROV-03``); esito positivo nell'ambiente dell'immagine;
- con la regola attiva un'immagine dichiarata diversa rende la fase da rifare,
  senza la regola resta conclusa; il report riporta l'esito della regola.

4. Comandi Bash e scenari di esecuzione
---------------------------------------
    ``<immagine>`` e' l'immagine del container della pipeline; quella corrente
    e' indicata in ``README.md``.

    1. Modalità locale standard (R di base con jsonlite, senza Bioconductor né
       dati reali):
       pytest tests/test_w26_congelamento.py -v

    2. Modalità container Docker standard (sottoinsieme ridotto con
       R/Bioconductor):
       docker run --rm \
         -e PYTHONPATH=/app/src \
         -e AMPLICON16S_R_DIR=/app/R \
         -v "$(pwd)":/app \
         -w /app \
         <immagine> \
         pytest -o cache_dir=/tmp/.pytest_cache tests/test_w26_congelamento.py -v

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
         pytest -o cache_dir=/tmp/.pytest_cache tests/test_w26_congelamento.py -v

5. Risultato atteso
-------------------
I conteggi li da' pytest (``pytest --collect-only -q``).

6. Razionale scientifico e sistemistico
---------------------------------------
- Una configurazione congelata certifica i risultati solo se il codice che li
  calcola e' identificato da un commit e l'ambiente R e' quello del file di
  blocco: cio' che non si puo' verificare deve fermare l'esecuzione, non
  passare in silenzio.
"""

from __future__ import annotations

import csv
import json
import re
import subprocess
from pathlib import Path
from typing import Any, ClassVar

import pytest
import yaml
from conftest import NEGATIVO, POSITIVO, Campione, crea_scenario, dichiarazione_minima
from sottoinsieme import motivo_pacchetti_r_assenti

import amplicon16s.cli as cli
import amplicon16s.runner.executor as executor
import amplicon16s.runner.provenienza as provenienza
from amplicon16s.config import defaults
from amplicon16s.config.resolve import risolvi
from amplicon16s.config.schema import Config, carica
from amplicon16s.errors.catalog import CATALOGO, Categoria
from amplicon16s.io_layer.artifacts import Fase
from amplicon16s.logging.logger import NOME_FILE_LOG, chiudi
from amplicon16s.rbridge.runner import cartella_r
from amplicon16s.report.builder import CARTELLA_REPORT, CARTELLA_TABELLE
from amplicon16s.runner.executor import EVENTO_CONTROLLI, NOME_PUNTO_DI_RIPRESA, Esecutore
from amplicon16s.runner.graph import Passo
from amplicon16s.runner.project import ProjectRun, StatoPasso
from amplicon16s.steps.s00_validate import ValidazioneIngressi

RADICE = Path(__file__).resolve().parents[1]
CONGELATA = RADICE / "dati" / "osd734" / "config_osd734.yaml"
senza_dati = pytest.mark.skipif(
    not CONGELATA.is_file(),
    reason="la cartella dati/ non e' presente (nell'immagine non viene copiata)",
)
IMMAGINE_DIVERSA = "amplicon16s@sha256:" + "a" * 64


def _tsv(percorso: Path) -> list[dict[str, str]]:
    """Le righe di una tabella TSV con intestazione."""
    with percorso.open(encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


@pytest.fixture(autouse=True)
def uscite_pulite():
    """Chiude le uscite del log prima e dopo ogni test, perche' nessun handler
    resti aperto sulla cartella temporanea.
    """
    chiudi()
    yield
    chiudi()


@pytest.fixture
def scenario(tmp_path):
    """Lo scenario su disco con tre campioni e le loro letture."""
    campioni = [
        Campione("ERX3000001", "NOD1D4.L1"),
        Campione("ERX3000003", "POS.P1.1", materiale=POSITIVO, posizione="Not Applicable"),
        Campione("ERX3000004", "BLANK.P1.1", materiale=NEGATIVO, posizione="Not Applicable"),
    ]
    return crea_scenario(tmp_path, campioni, con_letture=True)


def _scrivi_config(scenario, percorso: Path, **run: Any) -> Path:
    """Scrive la configurazione minima dello scenario, con i parametri del gruppo
    ``run`` indicati.
    """
    dati = dichiarazione_minima(scenario.config)
    dati["run"].update(run)
    percorso.write_text(yaml.safe_dump(dati, sort_keys=False), encoding="utf-8")
    return percorso


@pytest.fixture
def rigorosa(scenario, tmp_path) -> Path:
    """La configurazione dello scenario con la regola rigorosa attiva."""
    return _scrivi_config(scenario, tmp_path / "rigorosa.yaml", strict_provenance=True)


def _git(cartella: Path, *argomenti: str) -> None:
    """Esegue git nella cartella, con un'identita' di prova."""
    subprocess.run(
        ["git", "-C", str(cartella), "-c", "user.name=prova", "-c", "user.email=prova@example.org",
         *argomenti],
        check=True, capture_output=True,
    )


@pytest.fixture
def clone(tmp_path, monkeypatch) -> Path:
    """Un repository git di prova con ``src/`` e ``R/`` committati, che la
    provenienza considera la radice del repository della pipeline.
    """
    radice = tmp_path / "clone"
    (radice / "src").mkdir(parents=True)
    (radice / "R").mkdir()
    (radice / "src" / "modulo.py").write_text("x = 1\n", encoding="utf-8")
    (radice / "R" / "fase.R").write_text("x <- 1\n", encoding="utf-8")
    _git(radice, "init", "-q")
    _git(radice, "add", ".")
    _git(radice, "commit", "-q", "-m", "stato iniziale")
    monkeypatch.setattr(provenienza, "RADICE_REPOSITORY", radice)
    return radice


def _cli(*argomenti: str, capsys) -> int:
    """Esegue la riga di comando e ne restituisce il codice di uscita."""
    codice = cli.main(list(argomenti))
    capsys.readouterr()
    return codice


def _controlli(radice: Path) -> dict[str, Any]:
    """L'ultimo evento dei controlli di avvio nel log strutturato."""
    righe = (Path(radice) / Fase.LOGS.value / NOME_FILE_LOG).read_text(encoding="utf-8").splitlines()
    return [e for e in map(json.loads, righe) if e.get("evento") == EVENTO_CONTROLLI][-1]


def _punto(radice: Path) -> dict[str, Any]:
    """Il punto di ripresa dichiarato dall'esecuzione."""
    percorso = Path(radice) / Fase.LOGS.value / f"{NOME_PUNTO_DI_RIPRESA}.json"
    return json.loads(percorso.read_text(encoding="utf-8"))


def _rifiutata(scenario, config: Path, codice: str, capsys) -> dict[str, Any]:
    """Avvia la validazione, verifica che sia rifiutata con il codice indicato
    prima di qualunque fase, e restituisce l'esito della regola nel log.
    """
    assert _cli("validate", "--config", str(config), capsys=capsys) == cli.USCITA_ARRESTO
    radice = scenario.config.io.out_root
    punto = _punto(radice)
    assert (punto["codice"], punto["origine"]) == (codice, "controlli di avvio")
    assert not (Path(radice) / Fase.CONFIG.value / "manifest_S0.json").exists()
    assert not list(Path(radice).rglob("manifest_S*.json"))
    return _controlli(radice)["provenienza_rigorosa"]


# --------------------------------------------------------------------------- #
# 1. La regola rigorosa sulla provenienza                                      #
# --------------------------------------------------------------------------- #


def test_la_regola_rigorosa_e_disattivata_per_difetto(scenario):
    """
    **Obiettivo**: Verificare che ``run.strict_provenance`` valga falso per
    difetto e nella configurazione di esempio, che accetti solo un booleano, e
    che senza la regola l'impronta su cui una fase e' calcolata non contenga la
    chiave ``provenienza``.

    **Razionale scientifico e sistemistico**: Durante lo sviluppo il codice ha
    quasi sempre modifiche non committate: la regola deve essere una scelta
    della configurazione congelata, e la sua assenza non deve cambiare le
    impronte delle esecuzioni esistenti.
    """
    assert defaults.RUN_STRICT_PROVENANCE is False
    assert scenario.config.run.strict_provenance is False
    esempio = yaml.safe_load((RADICE / "config" / "config.example.yaml").read_text(encoding="utf-8"))
    assert esempio["run"]["strict_provenance"] is False
    dati = scenario.config.model_dump(mode="python")
    dati["run"]["strict_provenance"] = "true"
    with pytest.raises(ValueError):
        Config.model_validate(dati)
    calcolata = ValidazioneIngressi().calcolata_su(risolvi(scenario.config), {})
    assert "provenienza" not in calcolata


def test_con_la_regola_l_impronta_di_fase_comprende_sorgente_ambiente_e_immagine(scenario):
    """
    **Obiettivo**: Verificare che con la regola attiva l'impronta su cui una
    fase e' calcolata contenga l'impronta del sorgente della fase, quella del
    file di blocco dell'ambiente R e l'immagine dichiarata; che il resto
    dell'impronta sia invariato; che l'impronta cambi con l'immagine dichiarata.

    **Razionale scientifico e sistemistico**: Finche' sorgente e ambiente
    stanno solo nella provenienza, una differenza produce un avviso; nella
    configurazione congelata deve rendere la fase non piu' valida.
    """
    fase = ValidazioneIngressi()
    senza = fase.calcolata_su(risolvi(scenario.config), {})
    dati = scenario.config.model_dump(mode="python")
    dati["run"]["strict_provenance"] = True
    attiva = Config.model_validate(dati)
    con = fase.calcolata_su(risolvi(attiva), {})
    assert {k: v for k, v in con.items() if k != "provenienza"} == senza
    aggiunta = con["provenienza"]
    assert set(aggiunta) == {"sorgente", "blocco_r", "immagine"}
    sorgente, _ = provenienza.impronta_sorgente(provenienza.file_del_sorgente(fase, cartella_r()))
    assert aggiunta["sorgente"] == sorgente and re.fullmatch(r"sha256:[0-9a-f]{64}", sorgente)
    assert re.fullmatch(r"sha256:[0-9a-f]{64}", aggiunta["blocco_r"])
    assert aggiunta["immagine"] == attiva.run.container

    dati["run"]["container"] = IMMAGINE_DIVERSA
    altra = fase.calcolata_su(risolvi(Config.model_validate(dati)), {})
    assert altra["provenienza"]["immagine"] == IMMAGINE_DIVERSA
    assert altra["provenienza"]["sorgente"] == aggiunta["sorgente"]


def test_la_verifica_di_git_distingue_i_quattro_casi(clone, tmp_path, monkeypatch):
    """
    **Obiettivo**: Verificare ``verifica_git``: repository pulito, commit
    restituito; file di ``src/`` o ``R/`` modificato o non tracciato,
    ``E-PROV-02`` con il nome del file; modifica fuori da ``src/`` e ``R/``,
    nessun rifiuto; script R eseguiti da una cartella fuori dal repository,
    ``E-PROV-01``; cartella che non e' un repository, ``E-PROV-01``; git non
    eseguibile, ``E-PROV-01``.

    **Razionale scientifico e sistemistico**: Il commit identifica il codice
    solo se i file eseguiti sono quelli del commit: ogni caso in cui non e'
    dimostrabile deve dare un rifiuto, e git non leggibile non e' un via libera.
    """
    script = clone / "R"
    codice, commit = provenienza.verifica_git(script)
    assert codice is None and re.fullmatch(r"[0-9a-f]{40}", commit)

    (clone / "config.yaml").write_text("run: {}\n", encoding="utf-8")
    assert provenienza.verifica_git(script)[0] is None

    (clone / "src" / "modulo.py").write_text("x = 2\n", encoding="utf-8")
    codice, dettaglio = provenienza.verifica_git(script)
    assert codice == "E-PROV-02" and "src/modulo.py" in dettaglio
    _git(clone, "checkout", "--", "src/modulo.py")
    (clone / "R" / "nuovo.R").write_text("y <- 1\n", encoding="utf-8")
    codice, dettaglio = provenienza.verifica_git(script)
    assert codice == "E-PROV-02" and "R/nuovo.R" in dettaglio
    (clone / "R" / "nuovo.R").unlink()

    altrove = tmp_path / "script_altrove"
    altrove.mkdir()
    codice, dettaglio = provenienza.verifica_git(altrove)
    assert codice == "E-PROV-01" and "fuori dal repository" in dettaglio

    spoglia = tmp_path / "senza_git"
    (spoglia / "R").mkdir(parents=True)
    monkeypatch.setattr(provenienza, "RADICE_REPOSITORY", spoglia)
    assert provenienza.verifica_git(spoglia / "R")[0] == "E-PROV-01"

    monkeypatch.setattr(provenienza, "RADICE_REPOSITORY", clone)

    def assente(*argomenti, **opzioni):
        raise FileNotFoundError("git")

    monkeypatch.setattr(provenienza.subprocess, "run", assente)
    assert provenienza.verifica_git(script)[0] == "E-PROV-01"


def test_con_git_non_leggibile_l_esecuzione_e_rifiutata(scenario, rigorosa, tmp_path, monkeypatch, capsys):
    """
    **Obiettivo**: Verificare che, con la regola attiva e il codice fuori da un
    repository git, l'avvio si fermi con codice di uscita 4 e ``E-PROV-01``
    prima di qualunque fase, che il log riporti la regola attiva con git non
    verificato, e che la stessa configurazione senza la regola venga eseguita.

    **Razionale scientifico e sistemistico**: Senza commit nulla identifica il
    codice che ha prodotto i risultati: la regola deve rifiutare, non lasciar
    passare, ed e' il comportamento che distingue una configurazione congelata.
    """
    spoglia = tmp_path / "senza_git"
    (spoglia / "R").mkdir(parents=True)
    monkeypatch.setattr(provenienza, "RADICE_REPOSITORY", spoglia)
    esito = _rifiutata(scenario, rigorosa, "E-PROV-01", capsys)
    assert esito["attiva"] is True and esito["git"] == {"verificato": False, "commit": None}
    assert "ambiente_r" not in esito

    libera = _scrivi_config(scenario, tmp_path / "libera.yaml")
    assert _cli("validate", "--config", str(libera), capsys=capsys) == 0
    assert _controlli(scenario.config.io.out_root)["provenienza_rigorosa"] == {"attiva": False}


def test_con_modifiche_non_committate_l_esecuzione_e_rifiutata(
    scenario, rigorosa, clone, monkeypatch, capsys
):
    """
    **Obiettivo**: Verificare che, con la regola attiva e un file di ``src/``
    modificato rispetto al commit, l'avvio si fermi con ``E-PROV-02`` prima di
    qualunque fase, con il file nel dettaglio del punto di ripresa; e che il
    report dell'esecuzione riporti il rifiuto nella colonna della regola.

    **Razionale scientifico e sistemistico**: Un risultato prodotto da codice
    modificato porterebbe il commit di un codice diverso: la provenienza
    registrata sarebbe falsa proprio dove deve certificare.
    """
    monkeypatch.setenv("AMPLICON16S_R_DIR", str(clone / "R"))
    (clone / "src" / "modulo.py").write_text("x = 2\n", encoding="utf-8")
    esito = _rifiutata(scenario, rigorosa, "E-PROV-02", capsys)
    assert esito["git"]["verificato"] is False
    radice = scenario.config.io.out_root
    assert "src/modulo.py" in _punto(radice)["dettaglio"]

    monkeypatch.delenv("AMPLICON16S_R_DIR")
    assert _cli("report", "--config", str(rigorosa), capsys=capsys) == 0
    righe = _tsv(Path(radice) / CARTELLA_REPORT / CARTELLA_TABELLE / "esecuzioni.tsv")
    assert righe[-1]["regola rigorosa sulla provenienza"].startswith("rifiutata: codice non identificato")


def test_le_discordanze_fra_ambiente_e_file_di_blocco():
    """
    **Obiettivo**: Verificare ``discordanze_ambiente``: nessuna discordanza
    con versioni uguali, anche scritte con separatori diversi; una per la
    versione di R, per un pacchetto assente, per una versione diversa, per una
    correzione assente o diversa da quella del file di blocco.

    **Razionale scientifico e sistemistico**: Versioni diverse di R e dei
    pacchetti di calcolo possono dare risultati numerici diversi; la correzione
    di dada2 cambia la tassonomia nei pareggi e ha lo stesso numero di versione
    in ogni ambiente che la dichiara: va confrontata per contenuto.
    """
    correzione = {"Base": "1.36.0", "File": "c.patch", "SHA256": "ab" * 32}
    blocco = {
        "R": {"Version": "4.5.2"},
        "Packages": {
            "dada2": {"Version": "1.36.0.1", "Patch": correzione},
            "randomForest": {"Version": "4.7-1.2"},
            "vegan": {"Version": "2.7-1"},
        },
    }
    ambiente = {
        "r": "4.5.2",
        "pacchetti": {
            "dada2": {"versione": "1.36.0.1", "correzione": dict(correzione)},
            "randomForest": {"versione": "4.7.1.2", "correzione": None},
            "vegan": {"versione": "2.7.1"},
        },
    }
    assert provenienza.discordanze_ambiente(blocco, ambiente) == []

    diverso = json.loads(json.dumps(ambiente))
    diverso["r"] = "4.3.3"
    diverso["pacchetti"]["dada2"]["correzione"] = None
    diverso["pacchetti"]["vegan"]["versione"] = "2.6-4"
    diverso["pacchetti"]["randomForest"] = None
    trovate = provenienza.discordanze_ambiente(blocco, diverso)
    assert trovate == [
        "R: atteso 4.5.2, trovato 4.3.3",
        "dada2: correzione diversa da quella del file di blocco",
        "randomForest: non installato",
        "vegan: atteso 2.7-1, trovato 2.6-4",
    ]
    altra = json.loads(json.dumps(ambiente))
    altra["pacchetti"]["dada2"]["correzione"]["SHA256"] = "cd" * 32
    assert len(provenienza.discordanze_ambiente(blocco, altra)) == 1


def test_un_ambiente_r_diverso_dal_file_di_blocco_e_rifiutato(
    scenario, tmp_path, monkeypatch, capsys
):
    """
    **Obiettivo**: Verificare che, con git in regola e un file di blocco che
    dichiara una versione di R e pacchetti diversi da quelli installati, l'avvio
    si fermi con ``E-PROV-03`` prima di qualunque fase, con le discordanze lette
    dalle librerie installate; e che un file di blocco introvabile dia lo
    stesso codice.

    **Razionale scientifico e sistemistico**: ``run.container`` e' una
    dichiarazione: cio' che la pipeline puo' verificare da dentro il container
    sono le versioni che i calcoli caricheranno, ed e' su quelle che il rifiuto
    deve basarsi.
    """
    motivo = motivo_pacchetti_r_assenti("jsonlite")
    if motivo is not None:
        pytest.skip(motivo)
    monkeypatch.setattr(executor, "verifica_git", lambda script: (None, "0" * 40))
    blocco = tmp_path / "blocco.lock"
    blocco.write_text(json.dumps({
        "R": {"Version": "0.0.1"},
        "Packages": {
            "dada2": {"Version": "0.0.1"},
            "jsonlite": {"Version": "0.0.1"},
            "pacchettoinesistente": {"Version": "1.0"},
        },
    }), encoding="utf-8")
    config = _scrivi_config(scenario, tmp_path / "c.yaml", strict_provenance=True, lockfile=str(blocco))
    esito = _rifiutata(scenario, config, "E-PROV-03", capsys)
    assert esito["git"] == {"verificato": True, "commit": "0" * 40}
    ambiente = esito["ambiente_r"]
    assert ambiente["verificato"] is False and ambiente["pacchetti_confrontati"] == 3
    assert ambiente["r"] != "0.0.1" and re.fullmatch(r"\d+\.\d+\.\d+", ambiente["r"])
    assert "R: atteso 0.0.1, trovato " + ambiente["r"] in ambiente["discordanze"]
    assert "pacchettoinesistente: non installato" in ambiente["discordanze"]
    assert any(d.startswith("jsonlite: atteso 0.0.1, trovato ") for d in ambiente["discordanze"])

    assente = _scrivi_config(
        scenario, tmp_path / "d.yaml", strict_provenance=True, lockfile=str(tmp_path / "manca.lock")
    )
    chiudi()
    esito = _rifiutata(scenario, assente, "E-PROV-03", capsys)
    assert esito["ambiente_r"] == {"verificato": False}


def test_nell_ambiente_dell_immagine_la_regola_e_soddisfatta(scenario, rigorosa, monkeypatch, capsys):
    """
    **Obiettivo**: Verificare che nell'ambiente R dell'immagine, con git in
    regola, la regola confronti tutti i pacchetti di ``renv.lock`` con quelli
    installati senza discordanze, compresa la correzione di dada2, che la
    validazione si concluda e che il manifesto di S0 registri l'impronta
    rigorosa; e che il report dichiari che cosa e' verificato e che cosa e'
    dichiarato.

    **Razionale scientifico e sistemistico**: E' la condizione in cui si
    produce il riferimento: se la verifica fallisse nell'immagine pubblicata,
    il file di blocco non descriverebbe l'ambiente che dichiara.
    """
    motivo = motivo_pacchetti_r_assenti("dada2", "phyloseq", "DESeq2")
    if motivo is not None:
        pytest.skip(motivo)
    monkeypatch.setattr(executor, "verifica_git", lambda script: (None, "0" * 40))
    assert _cli("validate", "--config", str(rigorosa), capsys=capsys) == 0
    radice = scenario.config.io.out_root
    esito = _controlli(radice)["provenienza_rigorosa"]
    attesi = json.loads((RADICE / "renv.lock").read_text(encoding="utf-8"))
    assert esito["ambiente_r"] == {
        "verificato": True, "r": attesi["R"]["Version"],
        "pacchetti_confrontati": len(attesi["Packages"]), "discordanze": [],
    }
    assert esito["immagine_dichiarata"] == scenario.config.run.container
    manifesto = json.loads(next(Path(radice).rglob("manifest_S0.json")).read_text(encoding="utf-8"))
    assert set(manifesto["calcolata_su"]["provenienza"]) == {"sorgente", "blocco_r", "immagine"}

    assert _cli("report", "--config", str(rigorosa), capsys=capsys) == 0
    righe = _tsv(Path(radice) / CARTELLA_REPORT / CARTELLA_TABELLE / "esecuzioni.tsv")
    testo = righe[-1]["regola rigorosa sulla provenienza"]
    assert testo.startswith("verificati il codice (commit 000000000000") and "dichiarata" in testo


def test_con_la_regola_un_immagine_dichiarata_diversa_rifa_la_fase(
    scenario, rigorosa, tmp_path, monkeypatch, capsys
):
    """
    **Obiettivo**: Verificare che, conclusa S0 con la regola attiva, dichiarare
    un'immagine diversa renda S0 da eseguire con il motivo della regola
    rigorosa; e che senza la regola la stessa modifica lasci S0 conclusa.

    **Razionale scientifico e sistemistico**: E' la differenza fra avviso e
    vincolo: nella configurazione congelata un risultato non puo' restare
    valido sotto una dichiarazione d'ambiente diversa da quella con cui e'
    stato calcolato.
    """
    monkeypatch.setattr(Esecutore, "_provenienza_rigorosa", lambda self: ({"attiva": True}, None))
    assert _cli("validate", "--config", str(rigorosa), capsys=capsys) == 0
    fasi = {Passo.S0: ValidazioneIngressi()}
    stessa = ProjectRun(carica(rigorosa), passi=fasi).valuta().situazioni[Passo.S0]
    assert stessa.stato is StatoPasso.COMPLETATA

    altra = _scrivi_config(scenario, tmp_path / "altra.yaml", strict_provenance=True)
    dati = yaml.safe_load(altra.read_text(encoding="utf-8"))
    dati["run"]["container"] = IMMAGINE_DIVERSA
    altra.write_text(yaml.safe_dump(dati, sort_keys=False), encoding="utf-8")
    situazione = ProjectRun(carica(altra), passi=fasi).valuta().situazioni[Passo.S0]
    assert situazione.stato is StatoPasso.DA_ESEGUIRE
    assert "regola rigorosa sulla provenienza" in situazione.motivo


def test_senza_la_regola_un_immagine_dichiarata_diversa_non_rifa_la_fase(
    scenario, tmp_path, capsys
):
    """
    **Obiettivo**: Verificare che, conclusa S0 senza la regola, dichiarare
    un'immagine diversa lasci S0 conclusa, e che attivare la regola su
    un'esecuzione calcolata senza renda S0 da eseguire.

    **Razionale scientifico e sistemistico**: Un risultato calcolato senza le
    verifiche della regola non puo' essere promosso a risultato certificato
    attivandola dopo: va ricalcolato sotto la regola.
    """
    libera = _scrivi_config(scenario, tmp_path / "libera.yaml")
    assert _cli("validate", "--config", str(libera), capsys=capsys) == 0
    fasi = {Passo.S0: ValidazioneIngressi()}
    dati = yaml.safe_load(libera.read_text(encoding="utf-8"))
    dati["run"]["container"] = IMMAGINE_DIVERSA
    libera.write_text(yaml.safe_dump(dati, sort_keys=False), encoding="utf-8")
    assert ProjectRun(carica(libera), passi=fasi).valuta().situazioni[Passo.S0].stato is StatoPasso.COMPLETATA

    dati["run"]["strict_provenance"] = True
    libera.write_text(yaml.safe_dump(dati, sort_keys=False), encoding="utf-8")
    situazione = ProjectRun(carica(libera), passi=fasi).valuta().situazioni[Passo.S0]
    assert situazione.stato is StatoPasso.DA_ESEGUIRE
    assert "regola rigorosa sulla provenienza" in situazione.motivo


def test_i_codici_della_regola_dicono_che_cosa_fare():
    """
    **Obiettivo**: Verificare che i tre codici ``E-PROV-*`` siano nel catalogo
    con fase ``PROV`` e revisione umana, e che l'azione di ciascuno nomini il
    parametro ``run.strict_provenance`` e il rimedio (clone del repository,
    commit delle modifiche, immagine indicata).

    **Razionale scientifico e sistemistico**: Un rifiuto all'avvio e' utile
    solo se chi lo riceve capisce che cosa manca; nessuno dei tre casi si
    risolve ritentando.
    """
    attesi = {"E-PROV-01": "clone", "E-PROV-02": "git status", "E-PROV-03": "immagine"}
    for codice, rimedio in attesi.items():
        voce = CATALOGO[codice]
        assert voce.fase == "PROV" and voce.categoria is Categoria.REVISIONE_UMANA
        assert "run.strict_provenance" in voce.azione and rimedio in voce.azione, codice
