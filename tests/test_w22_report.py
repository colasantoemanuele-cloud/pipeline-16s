r"""Suite di test della settimana 22: il report di esecuzione.

1. Inquadramento nel Piano Operativo
------------------------------------
Settimana 22 (W22), Fase F6 (modulo del report: il documento che accompagna
ogni esecuzione, ricavato dalla cartella di output).

2. Moduli sorgente coperti
--------------------------
* ``src/amplicon16s/report/builder.py``, ``src/amplicon16s/report/templates/``
* ``src/amplicon16s/cli.py`` (sottocomando ``report``)
* ``src/amplicon16s/config/resolve.py`` (``parametri_dichiarati``, elenco
  ``dichiarati`` nella configurazione registrata)
* ``src/amplicon16s/config/defaults.py`` (``OBBLIGATORI``, ``STANDARD_DEL_METODO``)
* ``src/amplicon16s/runner/executor.py`` (controlli di avvio come evento
  strutturato, dichiarazione del file di blocco dei pacchetti R)
* ``src/amplicon16s/logging/logger.py`` (registro degli avvii
  ``99_logs/avvii.jsonl``, in aggiunta e senza rotazione)
* ``src/amplicon16s/runner/tracciamento.py`` (``componi``)

3. Cosa valuta questo file
--------------------------
- la configurazione registrata dice quali parametri il file di ingresso ha
  impostato, anche quando il valore coincide con il predefinito, senza che il
  digest ne dipenda;
- l'esecutore registra nel log, a ogni avvio e a ogni ripresa, l'esito di G15,
  G14 e G12 come evento strutturato, anche quando un controllo non e' superato,
  e che cosa il file di blocco dichiara su dada2;
- gli eventi che descrivono un'esecuzione stanno anche nel registro degli
  avvii, che non ruota: dopo molte riprese, con il log ruotato, il report le
  mostra tutte; senza il registro lo dichiara e non ricostruisce dal log;
- il report riporta i quindici gate di S0 e i controlli di ogni avvio, lo stato
  delle fasi, le versioni della configurazione con le differenze e le fasi
  eseguite con ciascuna, l'origine di ogni parametro (dichiarato, predefinito,
  derivato, aggiustato), i tentativi ripetuti e le degradazioni dai manifesti;
- la segnalazione in apertura dei parametri di ``DERIVATI_DAL_DATASET`` che
  coincidono con il valore di OSD-734, comunque siano stati impostati;
- gli avvisi di provenienza in evidenza quando il sorgente registrato differisce
  da quello della pipeline che genera il documento;
- il report e' identico byte per byte fra due generazioni, non altera lo stato
  ne' alcun file fuori dalla sua cartella, non avvia processi e non legge i dati
  grezzi; un artefatto alterato o un manifesto non valido sono segnalati e non
  riportati;
- il documento e' autoconsistente (nessuna risorsa esterna) e rispetta la
  tipografia dei criteri di conformita';
- sul sottoinsieme di prova, dopo la catena completa: tracciamento di tutti i
  passi per ogni campione, decisioni di S11, S12 e S13, dimensioni dell'oggetto
  finale, grafico del modello d'errore incorporato;
- sul dataset completo: il report della catena riporta i 960 campioni e
  l'oggetto finale.

4. Comandi Bash e scenari di esecuzione
---------------------------------------
    ``<immagine>`` e' l'immagine del container della pipeline; quella corrente
    e' indicata in ``README.md``.

    1. Modalità locale standard (R di base con jsonlite, senza Bioconductor né
       dati reali):
       pytest tests/test_w22_report.py -v

    2. Modalità container Docker standard (sottoinsieme ridotto con
       R/Bioconductor):
       docker run --rm \
         -e PYTHONPATH=/app/src \
         -e AMPLICON16S_R_DIR=/app/R \
         -v "$(pwd)":/app \
         -w /app \
         <immagine> \
         pytest -o cache_dir=/tmp/.pytest_cache tests/test_w22_report.py -v

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
         pytest -o cache_dir=/tmp/.pytest_cache tests/test_w22_report.py -v

5. Risultato atteso
-------------------
I conteggi li da' pytest (``pytest --collect-only -q``).

6. Razionale scientifico e sistemistico
---------------------------------------
- Il report e' il primo oggetto che il committente apre: cio' che afferma deve
  venire dagli artefatti, e non da un ricalcolo che potrebbe divergere da essi.
- Un valore predefinito indistinguibile da una scelta deliberata, o un valore
  di OSD-734 ereditato in silenzio su un altro dataset, produce un risultato
  formalmente valido e scientificamente infondato: va reso visibile.
- Un documento che cambia a ogni generazione, o che altera la cartella da cui e'
  ricavato, non puo' accompagnare un'esecuzione riproducibile.
"""

from __future__ import annotations

import functools
import csv
import hashlib
import json
import os
import re
import shutil
import subprocess
from pathlib import Path
from typing import Any, ClassVar

import pytest
import yaml
from conftest import (
    Campione,
    copia_esecuzione,
    crea_scenario,
    dichiarazione_minima,
    FORMATO,
    NEGATIVO,
    nomi_dei_gate,
    POSITIVO,
)
from sottoinsieme import attesi_dataset, motivo_pacchetti_r_assenti

import amplicon16s.cli as cli
import amplicon16s.report.builder as builder
from amplicon16s.config import defaults
from amplicon16s.config.resolve import parametri_dichiarati, risolvi
from amplicon16s.config.schema import Config, carica
from amplicon16s.errors.exceptions import errore
from amplicon16s.io_layer.artifacts import Fase
from amplicon16s.logging.logger import NOME_FILE_AVVII, NOME_FILE_LOG, chiudi, configura
from amplicon16s.report.builder import CARTELLA_REPORT, CARTELLA_TABELLE, NOME_REPORT, costruisci, genera
from amplicon16s.runner.executor import EVENTO_CONTROLLI, _blocco_r
from amplicon16s.runner.graph import Passo
from amplicon16s.runner.project import ProjectRun
from amplicon16s.runner.retry import dimezza
from amplicon16s.steps.base import PipelineStep, Produzione, StepContext
from amplicon16s.steps.s00_validate import ValidazioneIngressi

SEZIONI = (
    "evidenza", "stato", "gate", "configurazione", "provenienza", "decisioni",
    "tracciamento", "risultato",
)


# --------------------------------------------------------------------------- #
# Scenario e fasi di prova                                                     #
# --------------------------------------------------------------------------- #


class Doppia(PipelineStep):
    """Una fase di prova: scrive un artefatto, e secondo lo scenario fallisce al
    primo tentativo con un codice ripetibile o registra una degradazione.
    """

    parametri: ClassVar[tuple[str, ...]] = tuple(Config.model_fields)
    codice: ClassVar[str | None] = None
    degradazione: ClassVar[str | None] = None

    def __init__(self) -> None:
        self.tentativi = 0

    def calcola(self, contesto: StepContext) -> Produzione:
        """Fallisce, degrada o conclude secondo lo scenario della classe."""
        self.tentativi += 1
        if self.codice and self.tentativi == 1:
            raise errore(self.codice, f"fase di prova, lotto {contesto.config.run.batch_size}")
        if self.degradazione:
            contesto.degrada(self.degradazione, "frazione chimerica 0,31")
        artefatto = contesto.albero.scrivi_testo(
            self.cartella, f"{str(self.passo).lower()}.txt", str(self.passo)
        )
        return Produzione((artefatto,))


def _fasi() -> dict[Passo, PipelineStep]:
    """S0 vera e una fase di prova per ogni altra: S4 ritenta dimezzando il
    lotto dopo ``E-S4-02``, S6 registra la degradazione ``E-S6-02``.
    """
    speciali: dict[Passo, dict[str, Any]] = {
        Passo.S4: {"codice": "E-S4-02", "aggiustamenti": {"E-S4-02": dimezza("run.batch_size")}},
        Passo.S6: {"degradazione": "E-S6-02"},
    }
    fasi: dict[Passo, PipelineStep] = {Passo.S0: ValidazioneIngressi()}
    for passo in tuple(Passo)[1:]:
        fasi[passo] = type(f"Doppia{passo}", (Doppia,), {"passo": passo, **speciali.get(passo, {})})()
    return fasi


def _campioni() -> list[Campione]:
    """Due biologici, un positivo, un negativo e un tampone mai aperto, che la
    regola di configurazione riclassifica in controllo negativo.
    """
    return [
        Campione("ERX3000001", "NOD1D4.L1"),
        Campione("ERX3000002", "NOD1D4.L2"),
        Campione("ERX3000003", "POS.P1.1", materiale=POSITIVO, posizione="Not Applicable"),
        Campione("ERX3000004", "BLANK.P1.1", materiale=NEGATIVO, posizione="Not Applicable"),
        Campione("ERX3000005", "TUBO.N1", posizione=FORMATO.riclassificati[0]),
    ]


@pytest.fixture(autouse=True)
def uscite_pulite():
    """Chiude le uscite del log prima e dopo ogni test, perché nessun handler resti
    aperto sulla cartella temporanea.
    """
    chiudi()
    yield
    chiudi()


@pytest.fixture
def scenario(tmp_path):
    """Lo scenario su disco con i cinque campioni e le loro letture."""
    return crea_scenario(tmp_path, _campioni(), con_letture=True)


def _scrivi_config(scenario, percorso: Path, **gruppi: dict[str, Any]) -> Path:
    """Scrive un file di configurazione minimo: i percorsi e i parametri
    obbligatori dello scenario, piu' quelli indicati per gruppo. Tutto il
    resto vale il predefinito.
    """
    dati = dichiarazione_minima(scenario.config)
    for gruppo, valori in gruppi.items():
        dati.setdefault(gruppo, {}).update(valori)
    percorso.write_text(yaml.safe_dump(dati, sort_keys=False), encoding="utf-8")
    return percorso


@pytest.fixture
def file_config(scenario, tmp_path) -> Path:
    """La configurazione minima dello scenario, con ``filter.truncLen`` dichiarato a
    un valore diverso da quello di OSD-734.
    """
    return _scrivi_config(scenario, tmp_path / "config.yaml", filter={"truncLen": 120})


@pytest.fixture
def fasi(monkeypatch) -> dict[Passo, PipelineStep]:
    """Sostituisce le fasi della riga di comando con S0 vera e le fasi di prova."""
    fasi = _fasi()
    monkeypatch.setattr(cli, "_passi", lambda: fasi)
    return fasi


def _cli(*argomenti: str, capsys) -> tuple[int, str]:
    """Esegue la riga di comando e ne restituisce codice di uscita e output."""
    codice = cli.main(list(argomenti))
    return codice, capsys.readouterr().out


def _tabella(radice: Path, nome: str) -> list[dict[str, str]]:
    """Le righe di una tabella del report."""
    percorso = Path(radice) / CARTELLA_REPORT / CARTELLA_TABELLE / f"{nome}.tsv"
    with open(percorso, encoding="utf-8", newline="") as file:
        return list(csv.DictReader(file, delimiter="\t"))


def _documento(radice: Path) -> str:
    """Il testo del documento HTML del report."""
    return (Path(radice) / CARTELLA_REPORT / NOME_REPORT).read_text(encoding="utf-8")


def _impronte(radice: Path, *, report: bool) -> dict[str, str]:
    """L'impronta di ogni file sotto la cartella di output: quelli del report, o
    tutti gli altri.
    """
    return {
        str(p.relative_to(radice)): hashlib.sha256(p.read_bytes()).hexdigest()
        for p in sorted(Path(radice).rglob("*"))
        if p.is_file() and (p.relative_to(radice).parts[0] == CARTELLA_REPORT) == report
    }


def _stato(config: Config, fasi) -> list[tuple[str, str, str, str | None]]:
    """L'esito della valutazione dello stato, fase per fase."""
    valutazione = ProjectRun(config, passi=fasi).valuta()
    return [(str(p), s.stato.value, s.motivo, s.impronta) for p, s in valutazione.situazioni.items()]


def _eventi(radice: Path) -> list[dict[str, Any]]:
    """Gli eventi del log strutturato dell'esecuzione."""
    righe = (Path(radice) / Fase.LOGS.value / NOME_FILE_LOG).read_text(encoding="utf-8").splitlines()
    return [json.loads(r) for r in righe]


# --------------------------------------------------------------------------- #
# 1. Cio' che l'esecuzione registra perche' il report possa leggerlo           #
# --------------------------------------------------------------------------- #


def test_la_configurazione_registrata_elenca_i_parametri_dichiarati(scenario, tmp_path, capsys):
    """
    **Obiettivo**: Verificare che ``parametri_dichiarati`` restituisca i soli
    parametri scritti nel file, compreso uno il cui valore coincide con il
    predefinito, che ``resolved.yaml`` li elenchi in ``dichiarati``, e che il
    digest non dipenda da come un valore e' stato dato.

    **Razionale scientifico e sistemistico**: Confrontare i valori con i
    predefiniti non distingue un predefinito ereditato dallo stesso valore
    scritto di proposito: lo sa solo lo schema, al momento del caricamento.
    """
    minimo = _scrivi_config(scenario, tmp_path / "minimo.yaml")
    esplicito = _scrivi_config(
        scenario, tmp_path / "esplicito.yaml",
        filter={"maxEE": defaults.FILTER_MAXEE},
    )
    dichiarati = parametri_dichiarati(carica(minimo))
    assert "filter.maxEE" not in dichiarati
    assert {"io.out_root", "tax.ref_name", "run.threads"} <= set(dichiarati)
    assert "filter.maxEE" in parametri_dichiarati(carica(esplicito))
    assert risolvi(carica(minimo)).digest == risolvi(carica(esplicito)).digest

    assert _cli("validate", "--config", str(esplicito), capsys=capsys)[0] == 0
    registrata = yaml.safe_load(
        (Path(scenario.config.io.out_root) / Fase.CONFIG.value / "resolved.yaml").read_text()
    )
    assert registrata["dichiarati"] == list(parametri_dichiarati(carica(esplicito)))
    assert not set(registrata["dichiarati"]) & set(registrata["derivati"])


def test_i_controlli_di_avvio_sono_nel_log_a_ogni_avvio(scenario, file_config, capsys):
    """
    **Obiettivo**: Verificare che ogni avvio dell'esecutore, compreso quello che
    trova S0 gia' conclusa, scriva nel log un evento ``controlli_di_avvio`` con
    l'esito di G15, G14 e G12 in forma strutturata.

    **Razionale scientifico e sistemistico**: ``gates.json`` riporta i gate
    come li ha visti S0: una ripresa su un'altra macchina o con un riferimento
    sostituito non lo riscrive, e senza questo evento il suo esito resterebbe
    solo nel testo dei messaggi.
    """
    for _ in range(2):
        assert _cli("validate", "--config", str(file_config), capsys=capsys)[0] == 0
    controlli = [e for e in _eventi(scenario.config.io.out_root) if e.get("evento") == EVENTO_CONTROLLI]
    assert len(controlli) == 2
    for evento in controlli:
        assert [c["gate"] for c in evento["esiti"]] == ["G15", "G14", "G12"]
        assert all(c["eseguito"] and c["superato"] and c["violazioni"] == [] for c in evento["esiti"])
        assert evento["fino_a"] == "S0"


def test_il_report_mostra_tutte_le_riprese_anche_con_il_log_ruotato(
    scenario, tmp_path, monkeypatch, capsys
):
    """
    **Obiettivo**: Verificare che, dopo dodici avvii con un log che ruota
    conservando una sola rotazione, il log abbia perso i controlli dei primi
    avvii mentre il registro ``99_logs/avvii.jsonl`` li ha tutti, che il report
    mostri i dodici avvii con i loro controlli e la loro configurazione, e che
    senza il registro il report dichiari di non averlo invece di ricostruire la
    storia dal log.

    **Razionale scientifico e sistemistico**: Una storia degli avvii incompleta
    presentata come completa nasconde proprio le riprese piu' vecchie, quelle
    che hanno prodotto le fasi a monte: il registro e' solo in aggiunta e non
    ruota, e il report legge soltanto quello.
    """
    monkeypatch.setattr(
        cli, "configura", lambda out_root: configura(out_root, max_byte=1500, rotazioni=1)
    )
    config = tmp_path / "config.yaml"
    avvii = 12
    for numero in range(avvii):
        # Un parametro diverso a ogni ripresa: ogni avvio rifa' S0 e registra
        # una versione della configurazione.
        _scrivi_config(scenario, config, qc={"head_reads": 20 + numero})
        assert _cli("validate", "--config", str(config), capsys=capsys)[0] == 0
    chiudi()

    radice = Path(scenario.config.io.out_root)
    cartella = radice / Fase.LOGS.value
    assert sorted(p.name for p in cartella.iterdir()) == [
        NOME_FILE_AVVII, NOME_FILE_LOG, f"{NOME_FILE_LOG}.1",
    ]
    nel_log = [
        json.loads(riga)
        for file in (cartella / f"{NOME_FILE_LOG}.1", cartella / NOME_FILE_LOG)
        for riga in file.read_text(encoding="utf-8").splitlines()
    ]
    assert sum(e.get("evento") == EVENTO_CONTROLLI for e in nel_log) < avvii
    registro = [json.loads(r) for r in (cartella / NOME_FILE_AVVII).read_text(encoding="utf-8").splitlines()]
    assert all("evento" in e for e in registro)
    assert sum(e["evento"] == EVENTO_CONTROLLI for e in registro) == avvii

    assert _cli("report", "--config", str(config), capsys=capsys)[0] == 0
    esecuzioni = _tabella(radice, "esecuzioni")
    assert [r["avvio"] for r in esecuzioni] == [str(n) for n in range(1, avvii + 1)]
    assert [r["configurazione"] for r in esecuzioni] == (
        ["resolved.yaml"] + [f"resolved_{n}.yaml" for n in range(2, avvii + 1)]
    )
    assert {(r["fasi eseguite"], r["conclusione"]) for r in esecuzioni} == {("S0", "completata")}
    controlli = _tabella(radice, "controlli_di_avvio")
    assert len(controlli) == 3 * avvii and {r["esito"] for r in controlli} == {"superato"}
    assert len(_tabella(radice, "configurazioni")) == avvii

    (cartella / NOME_FILE_AVVII).unlink()
    assert _cli("report", "--config", str(config), capsys=capsys)[0] == 0
    tabelle = radice / CARTELLA_REPORT / CARTELLA_TABELLE
    assert not (tabelle / "esecuzioni.tsv").exists()
    assert not (tabelle / "controlli_di_avvio.tsv").exists()
    documento = _documento(radice)
    assert "Questa esecuzione non ha il registro degli avvii" in documento
    assert "<dt>avvii registrati</dt><dd>registro degli avvii assente</dd>" in documento
    assert {r["fasi eseguite (registro degli avvii)"] for r in _tabella(radice, "configurazioni")} == {
        "non determinabili"
    }


def test_un_controllo_di_avvio_non_superato_resta_nel_log_e_nel_report(
    scenario, file_config, capsys
):
    """
    **Obiettivo**: Verificare che, sostituito il riferimento tassonomico dopo
    S0, la ripresa si fermi su G12 con codice di uscita 4, che l'evento del log
    riporti G12 non superato con la violazione, e che il report mostri il
    controllo del secondo avvio come non superato accanto a quello superato del
    primo, con il punto di ripresa.

    **Razionale scientifico e sistemistico**: Il controllo di avvio e' il
    presidio contro un riferimento alterato fra S0 e S8: il suo esito negativo
    deve restare leggibile, non solo fermare l'esecuzione.
    """
    assert _cli("validate", "--config", str(file_config), capsys=capsys)[0] == 0
    Path(scenario.config.tax.ref_fasta).write_bytes(b">altro\nTTTT\n")
    assert _cli("validate", "--config", str(file_config), capsys=capsys)[0] == cli.USCITA_ARRESTO

    ultimo = [e for e in _eventi(scenario.config.io.out_root) if e.get("evento") == EVENTO_CONTROLLI][-1]
    g12 = ultimo["esiti"][-1]
    assert (g12["gate"], g12["superato"]) == ("G12", False)
    assert g12["violazioni"][0]["codice"].startswith("E-S0-")

    assert _cli("report", "--config", str(file_config), capsys=capsys)[0] == 0
    righe = _tabella(scenario.config.io.out_root, "controlli_di_avvio")
    assert [(r["avvio"], r["gate"], r["esito"]) for r in righe if r["gate"] == "G12"] == [
        ("1", "G12", "superato"), ("2", "G12", "NON SUPERATO"),
    ]
    assert g12["violazioni"][0]["codice"] in righe[-1]["violazioni"]
    assert "Punto di ripresa dichiarato" in _documento(scenario.config.io.out_root)


def test_il_file_di_blocco_dichiara_dada2_con_la_correzione(scenario, tmp_path, fasi, capsys):
    """
    **Obiettivo**: Verificare che ``_blocco_r`` legga dal file di blocco la
    versione di dada2 e la correzione dichiarata, che restituisca ``None`` per
    un file assente, e che il report riporti versione e correzione registrate
    all'avvio, oppure dichiari che non sono determinabili.

    **Razionale scientifico e sistemistico**: La correzione dei pareggi di
    assignTaxonomy rende la tassonomia riproducibile a parita' di seme: chi
    legge un risultato deve sapere se l'esecuzione la dichiarava.
    """
    blocco = tmp_path / "blocco.lock"
    blocco.write_text(json.dumps({"Packages": {"dada2": {
        "Version": "1.36.0.1",
        "Patch": {"Base": "1.36.0", "File": "dada2-1.36.0-pareggi.patch", "SHA256": "ab" * 32},
    }}}), encoding="utf-8")
    letto = _blocco_r(str(blocco))
    assert letto["dada2"] == {
        "versione": "1.36.0.1",
        "correzione": {"Base": "1.36.0", "File": "dada2-1.36.0-pareggi.patch", "SHA256": "ab" * 32},
    }
    assert letto["sha256"] == hashlib.sha256(blocco.read_bytes()).hexdigest()
    assert _blocco_r(str(tmp_path / "assente.lock")) is None

    config = _scrivi_config(scenario, tmp_path / "config.yaml", run={"lockfile": str(blocco)})
    assert _cli("validate", "--config", str(config), capsys=capsys)[0] == 0
    genera(scenario.config.io.out_root, fasi)
    documento = _documento(scenario.config.io.out_root)
    assert "<dd>1.36.0.1</dd>" in documento and "dada2-1.36.0-pareggi.patch" in documento

    shutil.rmtree(Path(scenario.config.io.out_root) / Fase.LOGS.value)
    genera(scenario.config.io.out_root, fasi)
    assert "La versione dichiarata di dada2 non è determinabile" in _documento(scenario.config.io.out_root)


# --------------------------------------------------------------------------- #
# 2. Contenuto del report                                                      #
# --------------------------------------------------------------------------- #


def test_il_report_dopo_la_sola_validazione(scenario, file_config, capsys):
    """
    **Obiettivo**: Verificare che dopo ``validate`` il report esista nella sua
    cartella con tutte le sezioni, riporti i quindici gate superati nell'ordine
    di esecuzione, S0 conclusa, S1 non conclusa e S9 disattivata, e il campione
    riclassificato dalla regola di configurazione.

    **Razionale scientifico e sistemistico**: Il report deve essere leggibile
    in qualunque stato dell'esecuzione, anche parziale: e' dopo un arresto che
    serve di piu'.
    """
    assert _cli("validate", "--config", str(file_config), capsys=capsys)[0] == 0
    codice, uscita = _cli("report", "--config", str(file_config), capsys=capsys)
    radice = Path(scenario.config.io.out_root)
    assert codice == 0 and str(radice / CARTELLA_REPORT / NOME_REPORT) in uscita

    documento = _documento(radice)
    assert re.findall(r'<section id="([a-z]+)">', documento) == list(SEZIONI)
    gate = _tabella(radice, "gate")
    assert [g["gate"] for g in gate] == list(nomi_dei_gate())
    assert {g["esito"] for g in gate} == {"superato"}
    stato = {r["fase"]: r["stato"] for r in _tabella(radice, "fasi")}
    assert stato["S0"] == "conclusa" and stato["S1"] == "non conclusa"
    assert stato["S9"] == "disattivata (phylo.enabled falso)"
    assert "non completa" in documento
    assert [r["campione"] for r in _tabella(radice, "riclassificati")] == ["TUBO.N1"]


def test_il_report_dice_l_origine_di_ogni_parametro(tmp_path, fasi, capsys):
    """
    **Obiettivo**: Verificare che la tabella dei parametri del report dica,
    per ogni parametro, se e' dichiarato nel file, preso dal predefinito,
    derivato o aggiustato da un tentativo ripetuto; che un parametro con un
    predefinito dichiarato al suo stesso valore risulti dichiarato; e che il
    report non segnali piu' alcun predefinito legato a un dataset.

    **Razionale scientifico e sistemistico**: Nessun predefinito viene da un
    dataset, quindi non c'e' un valore ereditato da segnalare in apertura:
    resta la tabella, da cui chi legge vede che cosa e' stato scelto e che
    cosa e' il valore standard del metodo.
    """
    scenario = crea_scenario(tmp_path / "origine", _campioni(), con_arricchimento=True,
                             con_letture=True)
    file_config = _scrivi_config(scenario, tmp_path / "origine.yaml", filter={"truncLen": 120})
    assert _cli("run", "--config", str(file_config), capsys=capsys)[0] == 0
    codice, uscita = _cli("report", "--config", str(file_config), capsys=capsys)
    assert codice == 0
    radice = scenario.config.io.out_root

    parametri = {r["parametro"]: r for r in _tabella(radice, "parametri")}
    assert parametri["filter.truncLen"]["origine"] == "dichiarato"
    for chiave in defaults.OBBLIGATORI:
        assert parametri[chiave]["origine"] == "dichiarato", chiave
    for chiave in defaults.STANDARD_DEL_METODO:
        assert parametri[chiave]["origine"] == "predefinito", chiave
    assert parametri["asv.len_min"]["origine"] == "derivato"
    assert parametri["asv.len_min"]["valore"] == str(
        120 - scenario.config.filter.trimLeft - defaults.ASV_LEN_TOL)
    assert parametri["run.batch_size"]["origine"] == "predefinito"
    assert parametri["run.batch_size"]["valore"] == "24"
    assert parametri["run.batch_size"]["aggiustamento"] == "S4: 24 -> 12 dopo E-S4-02"

    documento = _documento(radice)
    assert "tarat" not in uscita and "tarat" not in documento
    assert not (Path(radice) / CARTELLA_REPORT / CARTELLA_TABELLE / "valori_osd734.tsv").exists()

    # Un predefinito dichiarato, anche con lo stesso valore, diventa una scelta.
    _scrivi_config(scenario, file_config, filter={"truncLen": 120},
                   tax={"min_boot": defaults.TAX_MIN_BOOT})
    assert _cli("resume", "--config", str(file_config), capsys=capsys)[0] == 0
    assert _cli("report", "--config", str(file_config), capsys=capsys)[0] == 0
    parametri = {r["parametro"]: r for r in _tabella(radice, "parametri")}
    assert parametri["tax.min_boot"]["origine"] == "dichiarato"
    assert parametri["filter.maxEE"]["origine"] == "predefinito"


def test_tentativi_ripetuti_e_degradazioni_dai_manifesti(scenario, file_config, fasi, capsys):
    """
    **Obiettivo**: Verificare che il report riporti il tentativo ripetuto di S4
    con il codice, il parametro e i valori dichiarato e usato, e la
    degradazione di S6 con il suo codice, e che la catena risulti completa.

    **Razionale scientifico e sistemistico**: Sono le decisioni che la pipeline
    prende senza intervento umano: chi riceve il risultato deve poterle leggere
    con il loro motivo, dagli stessi manifesti che le registrano.
    """
    assert _cli("run", "--config", str(file_config), capsys=capsys)[0] == 0
    report = costruisci(scenario.config.io.out_root, fasi)
    tabelle = {t.nome: t for t in report.tabelle}
    assert tabelle["tentativi_ripetuti"].righe == (
        (Passo.S4, "E-S4-02", 1, "run.batch_size", "24", "24", "12", "run.batch_size dimezzato"),
    )
    # S0 dichiara per conto suo i controlli che lo scenario non ha (E-S0-17).
    (degradazione,) = [r for r in tabelle["degradazioni"].righe if r[0] is Passo.S6]
    assert degradazione[:2] == (Passo.S6, "E-S6-02") and "0,31" in degradazione[3]
    di_s0 = [r for r in tabelle["degradazioni"].righe if r[0] is Passo.S0]
    assert {r[1] for r in di_s0} == {"E-S0-17"}
    assert f"tentativi ripetuti: 1; degradazioni: {1 + len(di_s0)}" in report.segnalazioni
    assert "<dt>catena</dt><dd>completa</dd>" in report.html


def test_una_ripresa_con_configurazione_cambiata_mostra_entrambe_le_versioni(
    scenario, tmp_path, capsys
):
    """
    **Obiettivo**: Verificare che, ripresa l'esecuzione con un parametro
    cambiato, il report riporti le due versioni della configurazione
    registrata, la differenza fra loro, le fasi eseguite con ciascuna e, come
    parametri in uso, quelli dell'ultima.

    **Razionale scientifico e sistemistico**: Un'esecuzione ripresa puo' avere
    fasi calcolate con configurazioni diverse: senza la storia delle versioni
    non si sa con quali parametri e' stato prodotto ciascun artefatto.
    """
    config = _scrivi_config(scenario, tmp_path / "config.yaml")
    assert _cli("validate", "--config", str(config), capsys=capsys)[0] == 0
    _scrivi_config(scenario, config, qc={"head_reads": 20})
    assert _cli("validate", "--config", str(config), capsys=capsys)[0] == 0
    assert _cli("report", "--config", str(config), capsys=capsys)[0] == 0

    radice = scenario.config.io.out_root
    prima, seconda = _tabella(radice, "configurazioni")
    assert (prima["file"], seconda["file"], seconda["precedente"]) == (
        "resolved.yaml", "resolved_2.yaml", "resolved.yaml",
    )
    assert prima["digest"] != seconda["digest"]
    assert seconda["differenze dalla precedente"] == "qc.head_reads: 10000 -> 20"
    assert (prima["fasi eseguite (registro degli avvii)"], seconda["fasi eseguite (registro degli avvii)"]) == ("avvio 1: S0", "avvio 2: S0")
    parametri = {r["parametro"]: r for r in _tabella(radice, "parametri")}
    assert (parametri["qc.head_reads"]["valore"], parametri["qc.head_reads"]["origine"]) == ("20", "dichiarato")
    assert [r["configurazione"] for r in _tabella(radice, "esecuzioni")] == ["resolved.yaml", "resolved_2.yaml"]
    fase = next(r for r in _tabella(radice, "fasi") if r["fase"] == "S0")
    assert fase["configurazione dell'ultima esecuzione (registro degli avvii)"] == "resolved_2.yaml"


def test_gli_avvisi_di_provenienza_sono_in_evidenza(scenario, file_config, fasi, monkeypatch, capsys):
    """
    **Obiettivo**: Verificare che senza differenze il report dichiari l'assenza
    di avvisi, e che con un sorgente diverso da quello registrato nei manifesti
    ogni fase conclusa compaia fra gli avvisi in apertura, prima dello stato
    della catena, e nella colonna della sezione di provenienza.

    **Razionale scientifico e sistemistico**: A parita' di versione una
    provenienza diversa non rifa' la fase: l'unica difesa contro un risultato
    calcolato da un codice diverso da quello che si crede e' vederlo scritto.
    """
    assert _cli("run", "--config", str(file_config), capsys=capsys)[0] == 0
    radice = scenario.config.io.out_root
    invariato = costruisci(radice, fasi)
    assert "Nessun avviso" in invariato.html
    assert "avvisi_provenienza" not in {t.nome for t in invariato.tabelle}
    assert not any(r[-1] for t in invariato.tabelle if t.nome == "provenienza" for r in t.righe)

    monkeypatch.setattr(builder, "impronta_sorgente", lambda file: ("sha256:" + "0" * 64, {}))
    report = costruisci(radice, fasi)
    avvisi = next(t for t in report.tabelle if t.nome == "avvisi_provenienza")
    assert [r[0] for r in avvisi.righe] == [str(p) for p in Passo if p is not Passo.S9]
    assert all("sorgente diverso" in r[1] for r in avvisi.righe)
    assert any(s.startswith("avvisi di provenienza su S0, S1, ") for s in report.segnalazioni)
    assert report.html.index('<caption>Avvisi di provenienza</caption>') < report.html.index('id="stato"')
    provenienza = next(t for t in report.tabelle if t.nome == "provenienza")
    assert all(r[-1] for r in provenienza.righe)


# --------------------------------------------------------------------------- #
# 3. Riproducibilita' e stato                                                  #
# --------------------------------------------------------------------------- #


def test_il_report_e_identico_fra_due_generazioni_e_non_altera_lo_stato(
    scenario, file_config, fasi, capsys
):
    """
    **Obiettivo**: Verificare che due generazioni dalla stessa esecuzione diano
    gli stessi byte per il documento e per ogni tabella, che nessun file fuori
    dalla cartella del report cambi (log compreso), e che la valutazione dello
    stato dia lo stesso esito prima e dopo.

    **Razionale scientifico e sistemistico**: Il report e' funzione del solo
    contenuto della cartella: una data di generazione o una riga di log in piu'
    lo renderebbero diverso a ogni lettura, e scrivere fra gli artefatti
    invaliderebbe le fasi che descrive.
    """
    assert _cli("run", "--config", str(file_config), capsys=capsys)[0] == 0
    radice = Path(scenario.config.io.out_root)
    config = carica(file_config)
    stato, fuori = _stato(config, fasi), _impronte(radice, report=False)

    assert _cli("report", "--config", str(file_config), capsys=capsys)[0] == 0
    primo = _impronte(radice, report=True)
    assert _cli("report", "--config", str(file_config), capsys=capsys)[0] == 0
    assert _impronte(radice, report=True) == primo
    assert f"{CARTELLA_REPORT}/{NOME_REPORT}" in primo and len(primo) > 8

    assert _impronte(radice, report=False) == fuori
    assert _stato(config, fasi) == stato
    assert all(s[1] in ("completata", "disattivata") for s in stato)


def test_il_report_non_usa_r_ne_i_dati_grezzi(scenario, file_config, fasi, monkeypatch, capsys):
    """
    **Obiettivo**: Verificare che, rimossi le letture, le tabelle dei metadati e
    il riferimento tassonomico, e reso impossibile avviare qualunque processo,
    il sottocomando ``report`` produca lo stesso documento di prima.

    **Razionale scientifico e sistemistico**: Il report deve potersi generare
    da un'esecuzione conclusa in precedenza, su una macchina che ha solo la
    cartella di output: ne' R ne' i dati grezzi vi sono necessari.
    """
    assert _cli("run", "--config", str(file_config), capsys=capsys)[0] == 0
    assert _cli("report", "--config", str(file_config), capsys=capsys)[0] == 0
    radice = Path(scenario.config.io.out_root)
    prima = _impronte(radice, report=True)

    shutil.rmtree(scenario.config.io.fastq_dir)
    for ingresso in (scenario.config.io.assay_table, scenario.config.io.study_table,
                     scenario.config.tax.ref_fasta):
        Path(ingresso).unlink()

    def vietato(*argomenti, **opzioni):
        raise AssertionError(f"il report ha avviato un processo: {argomenti}")

    monkeypatch.setattr(subprocess, "Popen", vietato)
    monkeypatch.setattr(os, "system", vietato)
    shutil.rmtree(radice / CARTELLA_REPORT)
    assert _cli("report", "--config", str(file_config), capsys=capsys)[0] == 0
    assert _impronte(radice, report=True) == prima


def test_un_artefatto_alterato_e_un_manifesto_non_valido_sono_segnalati(
    scenario, file_config, fasi, capsys
):
    """
    **Obiettivo**: Verificare che un artefatto che non corrisponde al checksum
    del suo manifesto sia segnalato in apertura e non riportato, e che un
    manifesto modificato dopo la scrittura renda la fase non valida e la catena
    non completa.

    **Razionale scientifico e sistemistico**: Un report che riportasse il
    contenuto di un file alterato gli darebbe l'autorevolezza dell'esecuzione
    che non lo ha prodotto.
    """
    assert _cli("run", "--config", str(file_config), capsys=capsys)[0] == 0
    radice = Path(scenario.config.io.out_root)
    gate = radice / Fase.INPUT_VALIDATION.value / "gates.json"
    gate.write_text(gate.read_text().replace('"superato": true', '"superato": false', 1))
    manifesto = radice / Fase.QC_PROFILES.value / "manifest_S1.json"
    manifesto.write_text(manifesto.read_text().replace('"metriche": {}', '"metriche": {"x": 1}'))

    report = costruisci(radice, fasi)
    assert "S0: gates.json non corrisponde al checksum del manifesto" in " ".join(report.segnalazioni)
    assert any(s.startswith("manifesto non valido: S1") for s in report.segnalazioni)
    assert "gate" not in {t.nome for t in report.tabelle}
    assert "non completa" in report.html
    assert "non è riportato nel documento" in report.html


def test_le_tabelle_di_una_generazione_precedente_non_restano(scenario, file_config, capsys):
    """
    **Obiettivo**: Verificare che un file presente nella cartella delle tabelle
    e non prodotto dalla generazione corrente venga rimosso, e che i file TSV
    riportino i conteggi senza separatore delle migliaia.

    **Razionale scientifico e sistemistico**: Una tabella rimasta da uno stato
    precedente dell'esecuzione sarebbe letta come parte del report corrente.
    """
    assert _cli("validate", "--config", str(file_config), capsys=capsys)[0] == 0
    assert _cli("report", "--config", str(file_config), capsys=capsys)[0] == 0
    tabelle = Path(scenario.config.io.out_root) / CARTELLA_REPORT / CARTELLA_TABELLE
    (tabelle / "superata.tsv").write_text("vecchia\n", encoding="utf-8")
    assert _cli("report", "--config", str(file_config), capsys=capsys)[0] == 0
    assert not (tabelle / "superata.tsv").exists()

    tabella = builder.Tabella("prova", "Prova", ("campione", "letture"), (("A", 25730344), ("B", 21123.5)))
    assert tabella.tsv() == "campione\tletture\nA\t25730344\nB\t21123.5\n"
    assert "<td>25.730.344</td>" in tabella.html() and "<td>21.123,5</td>" in tabella.html()


def test_il_documento_e_autoconsistente_e_rispetta_la_tipografia(scenario, file_config, fasi, capsys):
    """
    **Obiettivo**: Verificare che il documento non richiami alcuna risorsa
    esterna (stili incorporati, nessuno script, nessun collegamento di rete,
    collegamenti solo interni o alle proprie tabelle) e che ne' il documento
    ne' le sue tabelle contengano trattini lunghi, medi o doppi trattini.

    **Razionale scientifico e sistemistico**: Il documento viene aperto senza
    connessione, anche a distanza di anni; e segue lo stesso registro
    tipografico del resto del progetto.
    """
    assert _cli("run", "--config", str(file_config), capsys=capsys)[0] == 0
    genera(scenario.config.io.out_root, fasi)
    documento = _documento(scenario.config.io.out_root)
    assert documento.startswith("<!DOCTYPE html>") and "<style>" in documento
    assert not re.search(r"<(script|link|iframe)\b", documento)
    assert not re.search(r"(https?:)?//[a-z]", documento)
    collegamenti = re.findall(r'href="([^"]+)"', documento)
    assert collegamenti and all(c.startswith(("#", f"{CARTELLA_TABELLE}/")) for c in collegamenti)
    # I due trattini lunghi per codice: il carattere non compare nemmeno qui.
    tabelle = Path(scenario.config.io.out_root) / CARTELLA_REPORT / CARTELLA_TABELLE
    for testo in (documento, *(p.read_text(encoding="utf-8") for p in sorted(tabelle.iterdir()))):
        assert chr(0x2014) not in testo and chr(0x2013) not in testo and "--" not in testo
    modelli = Path(builder.__file__).parent / "templates"
    assert sorted(p.name for p in modelli.iterdir()) == ["report.html", "stile.css"]


def test_il_report_sostituisce_il_resoconto_provvisorio(file_config, capsys):
    """
    **Obiettivo**: Verificare che la riga di comando non esponga piu' il
    resoconto provvisorio, che l'aiuto del sottocomando ``report`` descriva il
    report, e che ``report`` su una cartella senza esecuzione esca con il
    codice 3 senza crearla.

    **Razionale scientifico e sistemistico**: Due documenti che descrivono la
    stessa esecuzione possono divergere: ne resta uno solo.
    """
    assert not hasattr(cli, "resoconto") and not hasattr(cli, "NOME_RESOCONTO")
    aiuto = cli.build_parser().format_help()
    assert "report HTML" in aiuto and "provvisorio" not in aiuto
    codice, uscita = _cli("report", "--config", str(file_config), capsys=capsys)
    assert codice == cli.USCITA_CONFIGURAZIONE and "Nessuna esecuzione" in uscita
    assert not Path(carica(file_config).io.out_root).exists()


# --------------------------------------------------------------------------- #
# 4. Sul sottoinsieme di prova e sul dataset completo                          #
# --------------------------------------------------------------------------- #


@functools.cache
def _sonda_bioc() -> str | None:
    """Perche' l'ambiente R richiesto non c'e', o ``None``: la sonda parte al
    primo uso, non all'importazione del modulo, e una volta sola.
    """
    return motivo_pacchetti_r_assenti(
        "dada2", "ggplot2", "ShortRead", "jsonlite", "phyloseq", "Biostrings", "decontam"
    )


@pytest.fixture
def bioc():
    """Richiede R con i pacchetti dell'intera catena: salta senza, ma in CI fallisce."""
    if _sonda_bioc() is not None:
        if os.environ.get("AMPLICON16S_RICHIEDI_BIOC") == "1":
            pytest.fail(f"Bioconductor e' richiesto in questo ambiente: {_sonda_bioc()}")
        pytest.skip(_sonda_bioc())


@pytest.fixture
def catena(bioc, finale_calcolata, tmp_path):
    """Una copia della catena completa sul sottoinsieme di prova, ripresa dalla
    riga di comando senza nulla da eseguire: la copia ha cosi' il log di un avvio.
    """
    run = copia_esecuzione(finale_calcolata, tmp_path)
    config = tmp_path / "config.yaml"
    config.write_text(yaml.safe_dump(run.config.model_dump(mode="json"), sort_keys=False), encoding="utf-8")
    assert cli.main(["resume", "--config", str(config)]) == 0
    chiudi()
    return run, config


def test_il_report_della_catena_completa_sul_sottoinsieme(catena, capsys):
    """
    **Obiettivo**: Verificare che, dopo la catena completa sul sottoinsieme di
    prova, il report riporti i quindici gate e i controlli dell'avvio, il
    tracciamento di tutti i passi per ogni campione dell'inventario, le soglie
    di S11, la modalita' di S12, le esclusioni di S13, le dimensioni
    dell'oggetto finale del manifesto di S14, i file consegnati e il grafico
    del modello d'errore incorporato, senza avvisi di provenienza.

    **Razionale scientifico e sistemistico**: E' il report come lo riceve il
    committente: ogni numero deve coincidere con l'artefatto da cui viene.
    """
    run, config = catena
    codice, uscita = _cli("report", "--config", str(config), capsys=capsys)
    assert codice == 0 and "avvisi di provenienza" not in uscita
    radice = Path(run.config.io.out_root)
    documento = _documento(radice)

    assert len(_tabella(radice, "gate")) == 15
    assert [(r["gate"], r["esito"]) for r in _tabella(radice, "controlli_di_avvio")] == [
        ("G15", "superato"), ("G14", "superato"), ("G12", "superato"),
    ]
    assert "<dt>catena</dt><dd>completa</dd>" in documento

    campioni = _tabella(radice, "tracciamento_per_campione")
    passi = ["grezze", "prefiltro", "filtrate", "denoised", "tabella", "senza_chimere",
             "lunghezza", "decontaminate", "finali"]
    assert list(campioni[0])[4:] == passi
    inventario = json.loads((radice / Fase.INPUT_VALIDATION.value / "inventario.json").read_text())
    assert len(campioni) == inventario["campioni"]
    assert all(r[p] != "" for r in campioni for p in passi)
    per_classe = _tabella(radice, "tracciamento_per_classe")
    assert [r["passo"] for r in per_classe] == passi
    assert int(per_classe[0]["biologico: campioni"]) == inventario["conteggi"]["biologico"]
    totali = {r["passo"]: int(r["totale"]) for r in _tabella(radice, "letture_per_classe")}
    assert totali["grezze"] == sum(int(r["grezze"]) for r in campioni)

    finale = run.albero.manifesto_passo(Passo.S14, Fase.FINAL)
    assert f"{finale.metriche['campioni']} campioni x {finale.metriche['varianti']} varianti" in documento
    assert [r["file"] for r in _tabella(radice, "file_finali")] == finale.metriche["file"]
    filtri = json.loads((radice / Fase.FINAL_INTERMEDI.value / "filtri_riepilogo.json").read_text())
    assert [r["filtro"] for r in _tabella(radice, "filtri_finali")] == filtri["ordine"]
    assert int(totali["finali"]) == filtri["letture"]["finali"]
    soglia = json.loads((radice / Fase.CONTROLS.value / "soglia.json").read_text())
    assert len(_tabella(radice, "soglie_profondita")) == len(soglia["per_piastra"]) + ("senza_piastra" in soglia)
    assert {r["modalità"] for r in _tabella(radice, "decontaminazione")} == {"aggregata", "per_piastra"}
    assert documento.count('src="data:image/png;base64,') >= 1


def test_il_report_della_catena_completa_e_riproducibile(catena, capsys):
    """
    **Obiettivo**: Verificare, sulla catena completa del sottoinsieme di prova,
    che due generazioni diano gli stessi byte, che nessun file fuori dalla
    cartella del report cambi e che la valutazione dello stato resti la stessa.

    **Razionale scientifico e sistemistico**: La proprieta' va provata sugli
    artefatti veri delle fasi, non solo su quelli delle fasi di prova.
    """
    run, config = catena
    radice = Path(run.config.io.out_root)
    stato, fuori = _stato(run.config, None), _impronte(radice, report=False)
    assert _cli("report", "--config", str(config), capsys=capsys)[0] == 0
    primo = _impronte(radice, report=True)
    assert _cli("report", "--config", str(config), capsys=capsys)[0] == 0
    assert _impronte(radice, report=True) == primo
    assert _impronte(radice, report=False) == fuori
    assert _stato(run.config, None) == stato
    assert all(s[1] in ("completata", "disattivata") for s in stato)


@pytest.mark.dati_reali
def test_il_report_sul_dataset_completo(catena_reale):
    """
    **Obiettivo**: Verificare che il report costruito dalla catena sul dataset
    completo, senza scrivere nella cartella condivisa, riporti i 960 campioni
    nel tracciamento, i 33 tamponi riclassificati, l'oggetto finale del
    manifesto di S14 e nessun avviso di provenienza.

    **Razionale scientifico e sistemistico**: Sul dataset di riferimento il
    report deve restituire i fatti accertati, a partire dagli artefatti di
    un'esecuzione appena prodotta dallo stesso codice.
    """
    run, _ = catena_reale
    report = costruisci(run.config.io.out_root)
    tabelle = {t.nome: t for t in report.tabelle}
    attesi = attesi_dataset()
    assert len(tabelle["tracciamento_per_campione"].righe) == attesi["campioni"]
    assert len(tabelle["riclassificati"].righe) == attesi["riclassificati"]
    assert "avvisi_provenienza" not in tabelle
    finale = run.albero.manifesto_passo(Passo.S14, Fase.FINAL)
    assert f"{finale.metriche['campioni']} campioni x" in report.html
    assert not any("non integro" in s or "non valido" in s for s in report.segnalazioni)
