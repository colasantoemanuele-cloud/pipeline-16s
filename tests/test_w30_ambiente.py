r"""Suite di test della settimana 30: versione di dada2, ambiente di esecuzione, catalogo.

1. Inquadramento nel Piano Operativo
------------------------------------
Settimana 30 (W30), Fase F8 (generalita'): la versione di dada2 accertata a
ogni esecuzione di S8, la regola rigorosa sulla provenienza eseguibile dal
codice di un clone senza variabili d'ambiente, il tempo massimo dei processi
R, il retry ridotto e il catalogo degli errori allineato al codice.

2. Moduli sorgente coperti
--------------------------
* ``src/amplicon16s/steps/s08_taxonomy.py``, ``R/00_ambiente.R``, ``R/08_taxonomy.R``
* ``src/amplicon16s/report/builder.py`` (versione di dada2 dichiarata e caricata)
* ``src/amplicon16s/runner/provenienza.py`` (``_git``, ``verifica_git``)
* ``scripts/esegui.py``, ``src/amplicon16s/cli.py``
* ``src/amplicon16s/config/schema.py``, ``src/amplicon16s/config/resolve.py``
  (``run.r_timeout_s``)
* ``src/amplicon16s/rbridge/runner.py`` (tempo massimo, ``E-R-05``)
* ``src/amplicon16s/errors/catalog.py`` (``RISERVATI``, ``E-S8-03``, ``E-R-05``)
* ``src/amplicon16s/gates/g01_g15.py`` (codici di G15)

3. Cosa valuta questo file
--------------------------
- S8 legge a ogni esecuzione la versione di dada2 che R carica: con la
  correzione dei pareggi la registra, senza la dichiara con ``E-S8-03``, e il
  report distingue la versione caricata da quella dichiarata nel file di blocco;
- git si legge senza dipendere dall'utente ne' dalla sua configurazione, e il
  rifiuto ``E-PROV-01`` dice la causa vera;
- ``scripts/esegui.py`` esegue il codice e gli script R del clone qualunque
  cosa dica l'ambiente, con una cartella personale non scrivibile;
- ``run.r_timeout_s``: facoltativo, intero positivo, fuori dall'impronta dei
  risultati, passato da ogni fase al ponte, che allo scadere solleva ``E-R-05``;
- il catalogo: ogni codice e' sollevato dal codice o riservato con il motivo,
  ogni codice sollevato e' nel catalogo, ``E-S0-10`` ha un solo significato;
- ogni rifiuto di G15 esce dalla riga di comando con il suo codice e la sua
  azione, con codice di uscita 3;
- il pacchetto vuoto dei controlli e il limite di memoria del ponte non
  esistono piu'.

4. Comandi Bash e scenari di esecuzione
---------------------------------------
    ``<immagine>`` e' l'immagine del container della pipeline; quella corrente
    e' indicata in ``README.md``.

    1. Modalità locale standard (R di base con jsonlite, senza Bioconductor né
       dati reali):
       pytest tests/test_w30_ambiente.py -v

    2. Modalità container Docker standard (sottoinsieme ridotto con
       R/Bioconductor):
       docker run --rm \
         -e PYTHONPATH=/app/src \
         -e AMPLICON16S_R_DIR=/app/R \
         -v "$(pwd)":/app \
         -w /app \
         <immagine> \
         pytest -o cache_dir=/tmp/.pytest_cache tests/test_w30_ambiente.py -v

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
         pytest -o cache_dir=/tmp/.pytest_cache tests/test_w30_ambiente.py -v

5. Risultato atteso
-------------------
I conteggi li da' pytest (``pytest --collect-only -q``).

6. Razionale scientifico e sistemistico
---------------------------------------
- La versione ufficiale di dada2 assegna i pareggi fra generi con un
  generatore che il seme non controlla: una tassonomia non ripetibile va
  dichiarata con i risultati, non scoperta confrontando due esecuzioni.
- Una regola sulla provenienza che funziona solo con le variabili giuste e
  l'utente giusto viene disattivata alla prima difficolta': deve funzionare
  con il comando documentato e dire la causa quando rifiuta.
- Un catalogo con codici che nessuno solleva, o con un codice per tre
  condizioni, non dice piu' che cosa e' successo ne' che cosa fare.
"""

from __future__ import annotations

import importlib.util
import inspect
import json
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

import pytest
import yaml
from conftest import Campione, copia_esecuzione, crea_scenario, dichiarazione_minima
from sottoinsieme import motivo_pacchetti_r_assenti

from amplicon16s.cli import main
from amplicon16s.config import defaults
from amplicon16s.config.resolve import PARAMETRI_SENZA_EFFETTO, risolvi
from amplicon16s.config.schema import ErroreConfigurazione, valida
from amplicon16s.errors.catalog import CATALOGO, RISERVATI, Categoria
from amplicon16s.io_layer.artifacts import AlberoOutput, Fase
from amplicon16s.logging.logger import chiudi
from amplicon16s.rbridge import runner as ponte
from amplicon16s.report.builder import CARTELLA_REPORT, NOME_REPORT, genera
from amplicon16s.runner import provenienza
from amplicon16s.runner.executor import Conclusione, Esecutore
from amplicon16s.runner.graph import Passo
from amplicon16s.runner.project import ProjectRun, passi_realizzati
from amplicon16s.steps.s08_taxonomy import NOME_AMBIENTE, NOME_TASSONOMIA, dada2_caricato

RADICE = Path(__file__).resolve().parents[1]
CORREZIONE = {"Base": "1.36.0", "File": "dada2-1.36.0-pareggi.patch", "SHA256": "ab" * 32}


@pytest.fixture(autouse=True)
def uscite_pulite():
    """Chiude le uscite del log prima e dopo ogni test, perché nessun handler resti
    aperto sulla cartella temporanea.
    """
    chiudi()
    yield
    chiudi()


@pytest.fixture
def bioc():
    """Richiede R con i pacchetti delle fasi: salta senza, ma in CI fallisce."""
    motivo = motivo_pacchetti_r_assenti(
        "dada2", "ggplot2", "ShortRead", "jsonlite", "phyloseq", "Biostrings"
    )
    if motivo is not None:
        if os.environ.get("AMPLICON16S_RICHIEDI_BIOC") == "1":
            pytest.fail(f"Bioconductor e' richiesto in questo ambiente: {motivo}")
        pytest.skip(motivo)


@pytest.fixture
def r_di_base():
    """Richiede Rscript con jsonlite: salta senza."""
    motivo = motivo_pacchetti_r_assenti("jsonlite")
    if motivo is not None:
        pytest.skip(motivo)


def _campioni() -> list[Campione]:
    """Quattro campioni biologici: quanto basta a una configurazione valida."""
    return [Campione(f"ERX300000{i}", f"NOD1D4.L{i}") for i in range(1, 5)]


# --------------------------------------------------------------------------- #
# 1. La versione di dada2 caricata da R                                        #
# --------------------------------------------------------------------------- #


def _libreria_finta(cartella: Path, versione: str, corretta: bool) -> Path:
    """Una libreria R con un pacchetto ``dada2`` che ha il solo ``DESCRIPTION``:
    basta a ``R/00_ambiente.R``, che legge versione e correzione dichiarata.
    """
    pacchetto = cartella / "dada2"
    pacchetto.mkdir(parents=True)
    righe = ["Package: dada2", f"Version: {versione}"]
    if corretta:
        righe += [f"Amplicon16sBase: {CORREZIONE['Base']}",
                  f"Amplicon16sPatch: {CORREZIONE['File']}",
                  f"Amplicon16sPatchSHA256: {CORREZIONE['SHA256']}"]
    (pacchetto / "DESCRIPTION").write_text("\n".join(righe) + "\n", encoding="utf-8")
    return cartella


def test_la_versione_caricata_si_ricava_da_cio_che_r_ha_letto():
    """
    **Obiettivo**: Verificare ``dada2_caricato``: con la correzione dichiarata
    dal pacchetto installato la versione risulta corretta; senza, o con il
    pacchetto assente, non lo e', e la versione e' quella letta (o nulla).

    **Razionale scientifico e sistemistico**: La correzione dei pareggi ha lo
    stesso codice della versione ufficiale salvo una funzione: la si riconosce
    da cio' che la costruzione dell'immagine scrive nel pacchetto, non dal
    numero di versione che chiunque puo' avere.
    """
    con = dada2_caricato({"pacchetti": {"dada2": {"versione": "1.36.0.1", "correzione": CORREZIONE}}})
    assert con == {"versione": "1.36.0.1", "correzione": CORREZIONE, "corretta": True}
    senza = dada2_caricato({"pacchetti": {"dada2": {"versione": "1.36.0", "correzione": None}}})
    assert senza == {"versione": "1.36.0", "correzione": None, "corretta": False}
    assert dada2_caricato({"pacchetti": {"dada2": None}})["corretta"] is False
    assert dada2_caricato({})["versione"] is None


@pytest.mark.parametrize(("versione", "corretta"), [("1.36.0.1", True), ("1.36.0", False)])
def test_s8_accerta_dada2_a_ogni_esecuzione_e_senza_correzione_degrada(
    r_di_base, tmp_path, monkeypatch, versione, corretta
):
    """
    **Obiettivo**: Verificare, con una libreria R finta messa davanti a quelle
    installate, che la lettura che S8 fa prima della tassonomia riporti la
    versione di dada2 che R caricherebbe; che con la correzione non registri
    alcuna degradazione; che con la versione ufficiale registri ``E-S8-03``
    con la versione nel dettaglio; e che la regola rigorosa non sia attiva
    nella configurazione, cioe' che il controllo non dipenda da essa.

    **Razionale scientifico e sistemistico**: Fuori dall'immagine la
    tassonomia nei pareggi fra generi non e' ripetibile: chi esegue senza la
    regola rigorosa deve saperlo dal manifesto e dal report, non scoprirlo
    confrontando due esecuzioni.
    """
    monkeypatch.setenv("R_LIBS", str(_libreria_finta(tmp_path / "lib", versione, corretta)))
    scenario = crea_scenario(tmp_path / "s", _campioni())
    assert scenario.config.run.strict_provenance is False
    run = ProjectRun(scenario.config)
    fase = run.fase(Passo.S8)
    contesto = run.contesto(Passo.S8, run.valuta()).ristretto(fase.parametri)
    esito = fase._accerta_dada2(contesto)
    assert (esito["versione"], esito["corretta"]) == (versione, corretta)
    codici = [d.codice for d in contesto.degradazioni]
    if corretta:
        assert codici == [] and esito["correzione"] == CORREZIONE
    else:
        assert codici == ["E-S8-03"]
        assert f"dada2 {versione}" in contesto.degradazioni[0].dettaglio
        assert "non e' ripetibile" in contesto.degradazioni[0].dettaglio
    letto = json.loads(
        (run.albero.cartella(Fase.LOGS) / "ambiente_r_S8.json").read_text(encoding="utf-8"))
    assert letto["pacchetti"]["dada2"]["versione"] == versione
    voce = CATALOGO["E-S8-03"]
    assert voce.categoria is Categoria.DEGRADAZIONE_AUTOMATICA and voce.fase == "S8"
    assert (ponte.cartella_r() / NOME_AMBIENTE).is_file()


def test_senza_la_correzione_s8_degrada_e_il_report_lo_mostra(
    bioc, catena_calcolata, tassonomia_calcolata, tmp_path, monkeypatch
):
    """
    **Obiettivo**: Verificare sulla catena ridotta che, con un dada2 che si
    dichiara versione ufficiale (una copia del pacchetto installato senza i
    campi della correzione, davanti alle altre librerie), S8 si concluda con
    la degradazione ``E-S8-03`` nel manifesto e la versione caricata nelle
    metriche; che il report mostri la versione caricata, l'assenza della
    correzione e il codice, distinti da cio' che dichiara il file di blocco; e
    che con il dada2 installato la stessa fase non degradi e registri la
    correzione presente.

    **Razionale scientifico e sistemistico**: Il report mostrava la versione
    dichiarata nel file di blocco anche quando R ne caricava un'altra: la
    dichiarazione e la misura vanno riportate entrambe, e separate.
    """
    libreria = tmp_path / "lib"
    libreria.mkdir()
    script = tmp_path / "ufficiale.R"
    script.write_text(
        'origine <- find.package("dada2")\n'
        f"destinazione <- {json.dumps(str(libreria))}\n"
        "file.copy(origine, destinazione, recursive = TRUE)\n"
        'meta <- file.path(destinazione, "dada2", "Meta", "package.rds")\n'
        "d <- readRDS(meta)\n"
        'tieni <- !grepl("^Amplicon16s", names(d$DESCRIPTION))\n'
        "d$DESCRIPTION <- d$DESCRIPTION[tieni]\n"
        'd$DESCRIPTION[["Version"]] <- "1.36.0"\n'
        "saveRDS(d, meta)\n"
        'f <- file.path(destinazione, "dada2", "DESCRIPTION")\n'
        "righe <- readLines(f)\n"
        'righe <- righe[!grepl("^Amplicon16s", righe)]\n'
        'righe[grepl("^Version:", righe)] <- "Version: 1.36.0"\n'
        "writeLines(righe, f)\n",
        encoding="utf-8",
    )
    subprocess.run([str(ponte.trova_rscript()), "--vanilla", str(script)], check=True,
                   capture_output=True)

    corretta, _ = tassonomia_calcolata
    manifesto = corretta.albero.manifesto_passo(Passo.S8, Fase.TAXONOMY)
    assert manifesto.metriche["dada2"]["corretta"] is True
    assert "E-S8-03" not in [d["codice"] for d in manifesto.degradazioni]

    monkeypatch.setenv("R_LIBS", str(libreria))
    run = copia_esecuzione(catena_calcolata, tmp_path / "esecuzione")
    esito = Esecutore(run, fino_a=Passo.S8).esegui()
    assert esito.conclusione is Conclusione.COMPLETATA
    manifesto = run.albero.manifesto_passo(Passo.S8, Fase.TAXONOMY)
    assert [d["codice"] for d in manifesto.degradazioni] == ["E-S8-03"]
    assert manifesto.metriche["dada2"] == {
        "versione": "1.36.0", "correzione": None, "corretta": False}
    # Il codice che calcola e' lo stesso: cambia solo cio' che il pacchetto
    # dichiara di se', e la tassonomia della versione ridotta non ha pareggi.
    assert (run.albero.cartella(Fase.TAXONOMY) / NOME_TASSONOMIA).read_bytes() == (
        corretta.albero.cartella(Fase.TAXONOMY) / NOME_TASSONOMIA).read_bytes()

    radice = Path(run.config.io.out_root)
    genera(radice)
    documento = (radice / CARTELLA_REPORT / NOME_REPORT).read_text(encoding="utf-8")
    assert "versione di dada2 caricata da R in S8" in documento
    assert "assente: versione ufficiale" in documento and "E-S8-03" in documento
    assert "versione di dada2 dichiarata nel file di blocco" in documento or (
        "La versione dichiarata di dada2 non è determinabile" in documento)

    monkeypatch.delenv("R_LIBS")
    copia = copia_esecuzione(tassonomia_calcolata, tmp_path / "corretta")
    genera(Path(copia.config.io.out_root))
    documento = (Path(copia.config.io.out_root) / CARTELLA_REPORT / NOME_REPORT).read_text(
        encoding="utf-8")
    assert "correzione dei pareggi nella versione caricata" in documento
    assert "presente (dada2-" in documento and "assente: versione ufficiale" not in documento


# --------------------------------------------------------------------------- #
# 2. La regola rigorosa senza variabili d'ambiente                             #
# --------------------------------------------------------------------------- #


def _git(cartella: Path, *argomenti: str) -> str:
    """Esegue git nella cartella con un'identita' fissa, e ne restituisce l'uscita."""
    ambiente = {**os.environ, "GIT_AUTHOR_NAME": "prova", "GIT_AUTHOR_EMAIL": "prova@esempio.it",
                "GIT_COMMITTER_NAME": "prova", "GIT_COMMITTER_EMAIL": "prova@esempio.it",
                "GIT_CONFIG_GLOBAL": os.devnull}
    return subprocess.run(["git", "-C", str(cartella), *argomenti], check=True,
                          capture_output=True, text=True, env=ambiente).stdout


@pytest.fixture
def clone(tmp_path) -> Path:
    """Un repository git con il codice corrente: ``src/``, ``R/``, lo script che
    esegue il clone e il file di blocco, in un solo commit.
    """
    if shutil.which("git") is None:
        pytest.skip("git non disponibile")
    radice = tmp_path / "clone"
    for nome in ("src", "R"):
        shutil.copytree(RADICE / nome, radice / nome,
                        ignore=shutil.ignore_patterns("__pycache__", "*.egg-info"))
    (radice / "scripts").mkdir()
    shutil.copy(RADICE / "scripts" / "esegui.py", radice / "scripts" / "esegui.py")
    if (RADICE / "renv.lock").is_file():
        shutil.copy(RADICE / "renv.lock", radice / "renv.lock")
    _git(radice, "init", "-q")
    _git(radice, "add", "-A")
    _git(radice, "commit", "-q", "-m", "codice di prova")
    return radice


def test_git_si_legge_senza_dipendere_dall_utente(clone, tmp_path, monkeypatch):
    """
    **Obiettivo**: Verificare che la lettura di git della regola rigorosa
    riesca con una cartella personale inesistente e con una configurazione
    globale illeggibile; che dichiari sicura la cartella del repository per il
    solo comando e non legga la configurazione globale ne' di sistema; e che
    un guasto riporti il messaggio di git.

    **Razionale scientifico e sistemistico**: In un container il repository
    montato appartiene a un utente diverso da quello del processo, e la
    cartella personale dell'immagine non e' di chi lancia: git rifiuterebbe il
    repository per una ragione che non riguarda il codice.
    """
    monkeypatch.setattr(provenienza, "RADICE_REPOSITORY", clone)
    monkeypatch.setenv("HOME", str(tmp_path / "non_esiste"))
    illeggibile = tmp_path / "gitconfig"
    illeggibile.write_text("[sezione rotta\n", encoding="utf-8")
    monkeypatch.setenv("GIT_CONFIG_GLOBAL", str(illeggibile))
    assert re.fullmatch(r"[0-9a-f]{40}", provenienza._git("rev-parse", "HEAD").strip())
    assert provenienza.verifica_git(clone / "R")[0] is None

    visti = {}
    vero = subprocess.run

    def spia(comando, **opzioni):
        visti["comando"], visti["ambiente"] = comando, opzioni["env"]
        return vero(comando, **opzioni)

    monkeypatch.setattr(provenienza.subprocess, "run", spia)
    provenienza._git("rev-parse", "HEAD")
    assert f"safe.directory={clone}" in visti["comando"]
    assert visti["ambiente"]["GIT_CONFIG_GLOBAL"] == os.devnull
    assert visti["ambiente"]["GIT_CONFIG_SYSTEM"] == os.devnull
    with pytest.raises(provenienza.GitNonLeggibile, match="unknown revision|bad revision|ambiguous"):
        provenienza._git("rev-parse", "ramo-che-non-esiste")


def test_il_rifiuto_della_regola_rigorosa_dice_la_causa_vera(clone, tmp_path, monkeypatch):
    """
    **Obiettivo**: Verificare che ``E-PROV-01`` distingua nel dettaglio le sue
    tre cause: il codice in esecuzione non sta in un clone (pacchetto
    installato), con il comando da usare; git non legge il repository, con il
    messaggio di git; gli script R eseguiti stanno altrove, con il nome della
    variabile che lo causa.

    **Razionale scientifico e sistemistico**: «git non legge un repository»
    per tre cause diverse mandava a cercare un problema di permessi dove il
    problema era il comando lanciato: la causa va detta, con il rimedio.
    """
    installato = tmp_path / "site-packages"
    (installato / "R").mkdir(parents=True)
    monkeypatch.setattr(provenienza, "RADICE_REPOSITORY", installato)
    codice, dettaglio = provenienza.verifica_git(installato / "R")
    assert codice == "E-PROV-01"
    assert "pacchetto installato" in dettaglio and "scripts/esegui.py" in dettaglio

    monkeypatch.setattr(provenienza, "RADICE_REPOSITORY", clone)
    codice, dettaglio = provenienza.verifica_git(tmp_path / "altrove")
    assert codice == "E-PROV-01"
    assert "AMPLICON16S_R_DIR" in dettaglio and "scripts/esegui.py" in dettaglio

    rotto = tmp_path / "rotto"
    (rotto / ".git").mkdir(parents=True)
    (rotto / "R").mkdir()
    monkeypatch.setattr(provenienza, "RADICE_REPOSITORY", rotto)
    codice, dettaglio = provenienza.verifica_git(rotto / "R")
    assert codice == "E-PROV-01" and "git non legge il repository" in dettaglio
    assert "not a git repository" in dettaglio
    assert "scripts/esegui.py" in CATALOGO["E-PROV-01"].azione
    assert "PYTHONPATH" not in CATALOGO["E-PROV-01"].azione


def test_lo_script_del_clone_esegue_la_regola_rigorosa_senza_variabili(
    r_di_base, clone, tmp_path
):
    """
    **Obiettivo**: Verificare che ``scripts/esegui.py``, lanciato dalla radice
    di un clone senza ``PYTHONPATH``, con ``AMPLICON16S_R_DIR`` che punta a
    un'altra cartella e con una cartella personale non scrivibile, superi con
    la regola rigorosa attiva la verifica di git e degli script R (nessun
    ``E-PROV-01`` ne' ``E-PROV-02``): l'esecuzione prosegue oppure si ferma
    solo sull'ambiente R diverso dal file di blocco (``E-PROV-03``), e in
    quel caso il comando di ripresa e' lo stesso script; e che non lasci file
    compilati nel clone.

    **Razionale scientifico e sistemistico**: Dentro l'immagine il comando
    installato esegue il codice copiato alla costruzione, con i suoi script R:
    la regola rigorosa lo rifiuta. Il comando documentato deve eseguire il
    clone senza che chi lancia debba conoscere le variabili che lo rendono
    possibile.
    """
    scenario = crea_scenario(tmp_path / "s", _campioni(), con_letture=True)
    dati = dichiarazione_minima(scenario.config)
    dati["run"]["strict_provenance"] = True
    config = tmp_path / "rigorosa.yaml"
    config.write_text(yaml.safe_dump(dati), encoding="utf-8")
    casa = tmp_path / "casa_di_un_altro"
    casa.mkdir()
    casa.chmod(0o500)
    ambiente = {k: v for k, v in os.environ.items() if k != "PYTHONPATH"}
    ambiente.update(AMPLICON16S_R_DIR=str(tmp_path / "script_dell_immagine"), HOME=str(casa))
    esito = subprocess.run(
        [sys.executable, "scripts/esegui.py", "validate", "--config", str(config)],
        cwd=clone, env=ambiente, capture_output=True, text=True, check=False,
    )
    uscita = esito.stdout + esito.stderr
    assert "E-PROV-01" not in uscita and "E-PROV-02" not in uscita, uscita
    assert esito.returncode in (0, 4), uscita
    if esito.returncode == 4:
        assert "E-PROV-03" in uscita, uscita
        assert "python3 scripts/esegui.py resume --config" in uscita
    else:
        assert "Validazione superata" in uscita
    assert not list((clone / "src").rglob("__pycache__"))
    assert _git(clone, "status", "--porcelain", "--", "src", "R").strip() == ""


# --------------------------------------------------------------------------- #
# 3. Il tempo massimo dei processi R                                           #
# --------------------------------------------------------------------------- #


def test_run_r_timeout_s_e_facoltativo_intero_e_senza_effetto(tmp_path):
    """
    **Obiettivo**: Verificare che ``run.r_timeout_s`` sia nullo per difetto,
    accetti un intero positivo, respinga zero, un negativo, un booleano e un
    numero con la virgola; che sia fra i parametri senza effetto; e che due
    configurazioni che differiscono solo per esso abbiano la stessa impronta
    dei risultati.

    **Razionale scientifico e sistemistico**: Un limite di tempo decide se
    una fase si conclude, non che cosa calcola: cambiarlo non deve far rifare
    ore di calcolo che darebbero gli stessi risultati.
    """
    scenario = crea_scenario(tmp_path, _campioni())
    assert scenario.config.run.r_timeout_s is None
    dati = scenario.config.model_dump(mode="python")
    assert "run.r_timeout_s" in PARAMETRI_SENZA_EFFETTO
    dati["run"]["r_timeout_s"] = 3600
    con = valida(dati)
    assert con.run.r_timeout_s == 3600
    assert risolvi(con).impronta_risultati == risolvi(scenario.config).impronta_risultati
    for valore in (0, -5, True, 1.5, "60"):
        dati["run"]["r_timeout_s"] = valore
        with pytest.raises(ErroreConfigurazione, match="run.r_timeout_s"):
            valida(dati)


def test_ogni_fase_passa_il_tempo_massimo_al_ponte(tmp_path, monkeypatch):
    """
    **Obiettivo**: Verificare che ogni modulo di fase che esegue uno script R
    passi al ponte ``run.r_timeout_s`` come tempo massimo; che il valore
    arrivi davvero al ponte in un'esecuzione di S1; che allo scadere il ponte
    sollevi ``E-R-05``, a revisione umana, la cui azione nomina il parametro;
    e che il ponte non abbia piu' un limite di memoria.

    **Razionale scientifico e sistemistico**: Un processo R bloccato non deve
    poter bloccare l'orchestratore, ma il limite ha senso solo se chi esegue
    lo puo' dichiarare: prima esisteva nel ponte e nessuna fase lo passava.
    """
    passa = "tempo_massimo_s=contesto.config.run.r_timeout_s"
    for percorso in sorted((RADICE / "src" / "amplicon16s" / "steps").glob("s[0-9][0-9]_*.py")):
        testo = percorso.read_text(encoding="utf-8")
        assert testo.count("esegui_script(") == testo.count(passa), percorso.name

    import amplicon16s.steps.s01_profile as s01

    visti = {}

    def spia(script, parametri, *argomenti, **opzioni):
        visti.update(opzioni)
        raise RuntimeError("fermato dal test")

    monkeypatch.setattr(s01, "esegui_script", spia)
    scenario = crea_scenario(tmp_path, _campioni(), con_letture=True,
                             sovrascrivi={"run": {"r_timeout_s": 7, "threads": 1}})
    run = ProjectRun(scenario.config)
    assert Esecutore(run, fino_a=Passo.S0).esegui().conclusione is Conclusione.COMPLETATA
    fase = run.fase(Passo.S1)
    contesto = run.contesto(Passo.S1, run.valuta()).ristretto(fase.parametri)
    with pytest.raises(RuntimeError, match="fermato dal test"):
        fase.calcola(contesto)
    assert visti["tempo_massimo_s"] == 7

    voce = CATALOGO["E-R-05"]
    assert voce.categoria is Categoria.REVISIONE_UMANA and "run.r_timeout_s" in voce.azione
    parametri = inspect.signature(ponte.esegui_script).parameters
    assert "tempo_massimo_s" in parametri and "limite_memoria_byte" not in parametri
    assert not hasattr(ponte, "_limita_memoria")


# --------------------------------------------------------------------------- #
# 4. Il catalogo allineato al codice                                           #
# --------------------------------------------------------------------------- #

_LETTERALE = re.compile(r'"(E-[A-Z0-9]+-\d{2})"')


def _sollevati() -> dict[str, set[str]]:
    """I codici che il codice nomina come letterali, con i file che li nominano.

    Si leggono i sorgenti Python (tranne il catalogo stesso) e gli script R,
    senza le righe di commento; dal modulo dei gate si toglie la tabella
    dichiarativa ``CONTROLLI``, che elenca i controlli di G15 ma non solleva.
    """
    trovati: dict[str, set[str]] = {}
    sorgenti = [p for p in (RADICE / "src" / "amplicon16s").rglob("*.py") if p.name != "catalog.py"]
    sorgenti += list(ponte.cartella_r().rglob("*.R"))
    for percorso in sorgenti:
        testo = percorso.read_text(encoding="utf-8")
        if percorso.name == "g01_g15.py":
            testo = testo[:testo.index("CONTROLLI: Final")] + testo[testo.index("_PER_CODICE: Final"):]
        testo = "\n".join(r for r in testo.split("\n") if not r.lstrip().startswith("#"))
        for codice in _LETTERALE.findall(testo):
            trovati.setdefault(codice, set()).add(percorso.name)
    return trovati


def test_ogni_codice_del_catalogo_e_sollevato_o_riservato_e_viceversa():
    """
    **Obiettivo**: Verificare che ogni codice del catalogo sia nominato, come
    letterale, da un sorgente Python o R (fuori dai commenti e dalla tabella
    dichiarativa dei controlli), oppure sia in ``RISERVATI``
    con un motivo non vuoto; che nessun codice riservato sia anche sollevato;
    che ogni codice nominato dal codice esista nel catalogo; e che la
    whitelist predefinita del retry sia l'insieme dei codici che il catalogo
    ammette al retry.

    **Razionale scientifico e sistemistico**: Il catalogo e' la sola fonte di
    verita' dei codici: una voce che nessuno solleva documenta un
    comportamento che non esiste, e un codice sollevato fuori catalogo ferma
    l'esecuzione senza dire che cosa fare. Il test vede i nomi, non
    l'esecuzione: che ogni codice sia davvero raggiungibile lo provano i test
    che lo provocano, elencati in ``test_w15_consolidamento.py``.
    """
    sollevati = _sollevati()
    senza_origine = set(CATALOGO) - set(sollevati) - set(RISERVATI)
    assert not senza_origine, f"codici del catalogo mai sollevati ne' riservati: {sorted(senza_origine)}"
    fuori_catalogo = set(sollevati) - set(CATALOGO)
    assert not fuori_catalogo, {c: sorted(sollevati[c]) for c in fuori_catalogo}
    assert set(RISERVATI) <= set(CATALOGO)
    assert not set(RISERVATI) & set(sollevati)
    assert all(motivo.strip() for motivo in RISERVATI.values())
    ammessi = {c for c, v in CATALOGO.items() if v.ammette_retry}
    assert set(defaults.RETRY_WHITELIST) == ammessi == {"E-S2-03", "E-S4-02"}
    for codice in ("E-S3-01", "E-S5-01"):
        assert CATALOGO[codice].categoria is Categoria.REVISIONE_UMANA
        assert "Nessun nuovo tentativo automatico" in CATALOGO[codice].azione


def test_e_s0_10_ha_un_solo_significato():
    """
    **Obiettivo**: Verificare che ``E-S0-10`` sia sollevato in un solo punto
    del codice, per il primer ancora presente nelle letture, e che le altre
    due condizioni di G10 abbiano codici propri (``E-S0-16`` il motivo
    conservato assente, ``E-S0-18`` il primer cercato con un taglio iniziale).

    **Razionale scientifico e sistemistico**: Un codice per tre condizioni
    porta un'azione giusta per una sola: chi legge «togli il primer» quando
    manca il segnale della regione amplificata corregge la cosa sbagliata.
    """
    gate = (RADICE / "src" / "amplicon16s" / "gates" / "g01_g15.py").read_text(encoding="utf-8")
    corpo = gate[gate.index("_PER_CODICE: Final"):]
    assert corpo.count('"E-S0-10"') == 1
    assert corpo.count('"E-S0-16"') >= 1 and corpo.count('"E-S0-18"') >= 1
    assert "primer" in CATALOGO["E-S0-10"].sintesi
    assert "motivo" not in CATALOGO["E-S0-10"].sintesi
    sintesi = {CATALOGO[c].sintesi for c in ("E-S0-10", "E-S0-16", "E-S0-18")}
    assert len(sintesi) == 3


@pytest.mark.parametrize(("codice", "gruppo", "modifica"), [
    ("E-G15-02", "decontam", {"threshold": 0.0}),
    ("E-G15-03", "asv", {"len_tol": 500}),
    ("E-G15-04", "prev", {"min_fraction": 2}),
    ("E-G15-05", "qc", {"warn_frac_chimeric": 0.9, "stop_frac_chimeric": 0.5}),
    ("E-G15-06", "retry", {"max_attempts": 0}),
    ("E-G15-08", "asv", {"len_min": 100}),
    ("E-G15-09", "retry", {"whitelist": ["E-S3-01"]}),
    ("E-G15-10", "filter", {"truncLen": ...}),
    ("E-G15-11", "decontam", {"batch_column": "piastra"}),
    ("E-G15-12", "katharoseq", {"target_taxon": None}),
    ("E-G15-13", "out", {"study_columns": ["Sample Name", "sample-name"]}),
    ("E-G15-14", "meta", {"study_sample_id_column": "id"}),
    ("E-G15-99", "qc", {"parametro_che_non_esiste": 1}),
])
def test_ogni_rifiuto_di_g15_esce_dalla_riga_di_comando_con_il_codice(
    tmp_path, capsys, codice, gruppo, modifica
):
    """
    **Obiettivo**: Verificare, per ogni controllo di G15 (quelli dello schema
    e quelli di coerenza), che una configurazione che lo viola faccia uscire
    la riga di comando con il codice 3, e che il messaggio riporti il codice
    del catalogo della violazione e la riga «cosa fare» con la sua azione.

    **Razionale scientifico e sistemistico**: Chi lancia la pipeline da uno
    script vede l'uscita 3 e il testo: senza il codice non puo' cercare la
    voce del catalogo, e senza l'azione deve indovinare la correzione.
    """
    scenario = crea_scenario(tmp_path, _campioni())
    dati = dichiarazione_minima(scenario.config)
    if codice == "E-G15-14":
        dati["io"]["study_table"] = None
        dati["out"]["study_columns"] = []
    for chiave, valore in modifica.items():
        if valore is ...:
            dati[gruppo].pop(chiave)
        else:
            dati.setdefault(gruppo, {})[chiave] = valore
    percorso = tmp_path / "config.yaml"
    percorso.write_text(yaml.safe_dump(dati), encoding="utf-8")
    assert main(["validate", "--config", str(percorso)]) == 3
    uscita = capsys.readouterr().out
    assert "G15 ha respinto la configurazione" in uscita
    assert f"[{codice}]" in uscita, uscita
    assert f"cosa fare [{codice}]: {CATALOGO[codice].azione}" in uscita


def test_il_codice_morto_e_stato_rimosso():
    """
    **Obiettivo**: Verificare che il pacchetto vuoto ``amplicon16s.controls``
    non esista piu', che la categoria «retry poi revisione» sia stata tolta
    con i suoi due codici ora a revisione umana, e che fra i parametri che
    un'azione correttiva puo' cambiare resti il solo ``run.batch_size``.

    **Razionale scientifico e sistemistico**: Codice che nessuno esegue va
    letto, mantenuto e spiegato come quello che calcola, senza calcolare
    nulla: toglierlo riduce cio' che un revisore deve verificare.
    """
    from amplicon16s.runner.retry import PARAMETRI_AGGIUSTABILI

    assert importlib.util.find_spec("amplicon16s.controls") is None
    assert [c.value for c in Categoria] == [
        "revisione_umana", "retry_automatico", "degradazione_automatica"]
    assert PARAMETRI_AGGIUSTABILI == ("run.batch_size",)
    assert not {p for p in passi_realizzati().values() if "E-S3-01" in p.aggiustamenti}
