r"""Suite di test della settimana 25: l'immagine pubblicata e l'ambiente dei test.

1. Inquadramento nel Piano Operativo
------------------------------------
Settimana 25 (W25), Fase F7 (immagine del container pubblicata e identificata
per digest di registro; pacchetti R delle analisi ecologiche nell'immagine;
suite eseguibile su una macchina con un solo processore).

2. Moduli sorgente coperti
--------------------------
* ``container/install_r_packages.R``, ``container/snapshot_renv_lock.R``,
  ``renv.lock``
* ``tests/sottoinsieme.py`` (``processori_disponibili``, ``THREAD_DI_PROVA``)
* ``src/amplicon16s/gates/g01_g15.py`` (G14, ``E-S0-14``)
* ``dati/osd734/config_osd734.yaml`` (``run.container``), ``README.md``,
  ``dati/osd734/README.md``, ``test.txt`` (sezione 1.3)

3. Cosa valuta questo file
--------------------------
- i thread chiesti dalle esecuzioni di prova si adattano ai processori
  utilizzabili dal processo: uno su una macchina con un solo processore, mai
  piu' di due; con un solo processore G14 accetta un thread e ne respinge due
  con ``E-S0-14``;
- ``renv.lock`` contiene i pacchetti delle analisi ecologiche (DESeq2, vegan,
  randomForest) accanto a quelli della pipeline, con Bioconductor 3.21 e dada2
  1.36.0.1 con la correzione; i due script del container elencano gli stessi
  pacchetti;
- nell'immagine i pacchetti delle analisi ecologiche si caricano, alle versioni
  di ``renv.lock``;
- la configurazione di OSD-734 dichiara in ``run.container`` un digest di
  registro reale, e i due README e ``test.txt`` indicano la stessa immagine.

4. Comandi Bash e scenari di esecuzione
---------------------------------------
    ``<immagine>`` e' l'immagine del container della pipeline; quella corrente
    e' indicata in ``test.txt``, sezione 1.3.

    1. Modalità locale standard (R di base con jsonlite, senza Bioconductor né
       dati reali):
       pytest tests/test_w25_immagine.py -v

    2. Modalità container Docker standard (sottoinsieme ridotto con
       R/Bioconductor):
       docker run --rm \
         -e PYTHONPATH=/app/src \
         -e AMPLICON16S_R_DIR=/app/R \
         -v "$(pwd)":/app \
         -w /app \
         <immagine> \
         pytest -o cache_dir=/tmp/.pytest_cache tests/test_w25_immagine.py -v

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
         pytest -o cache_dir=/tmp/.pytest_cache tests/test_w25_immagine.py -v

5. Risultato atteso
-------------------
Vedi ``test.txt``, scheda W25.

6. Razionale scientifico e sistemistico
---------------------------------------
- Un'immagine identificata per digest di registro e' la stessa per chiunque la
  scarichi; un identificativo locale vale solo sulla macchina che l'ha costruita.
- Aggiungere pacchetti a un ambiente bloccato puo' trascinare versioni nuove di
  dipendenze gia' presenti: il file di blocco deve mostrare sole aggiunte.
- Una suite che fallisce su una macchina con un solo processore per i thread
  chiesti dai test segnala un difetto che la pipeline non ha.
"""

from __future__ import annotations

import json
import os
import re
import subprocess
from pathlib import Path

import pytest
import yaml
from conftest import NEGATIVO, POSITIVO, Campione, crea_scenario
from sottoinsieme import (
    THREAD_DI_PROVA,
    dati_config,
    motivo_pacchetti_r_assenti,
    processori_disponibili,
)

from amplicon16s.config.schema import valida
from amplicon16s.gates.g01_g15 import Contesto
from amplicon16s.gates.registry import esegui_gate
from amplicon16s.rbridge.runner import trova_rscript

RADICE = Path(__file__).resolve().parents[1]
ECOLOGIA = ("DESeq2", "vegan", "randomForest")
#: Un riferimento a un'immagine in un registro, ancorato per digest.
DIGEST_DI_REGISTRO = re.compile(r"^ghcr\.io/[a-z0-9-]+/[a-z0-9._-]+@sha256:[0-9a-f]{64}$")


def test_i_thread_di_prova_si_adattano_ai_processori(tmp_path, monkeypatch):
    """
    **Obiettivo**: Verificare che ``processori_disponibili`` conti i processori
    dell'insieme di affinita' del processo, che i thread delle esecuzioni di
    prova siano il minimo fra due e quei processori, e che la configurazione
    del sottoinsieme di prova li usi.

    **Razionale scientifico e sistemistico**: In un container limitato o sotto
    ``taskset`` i processori utilizzabili sono meno di quelli della macchina:
    un test che ne chiede di piu' viene fermato da G14, correttamente.
    """
    assert processori_disponibili() == len(os.sched_getaffinity(0))
    assert THREAD_DI_PROVA == min(2, processori_disponibili()) >= 1
    assert dati_config(tmp_path)["run"]["threads"] == THREAD_DI_PROVA
    monkeypatch.setattr(os, "sched_getaffinity", lambda pid: {0})
    assert processori_disponibili() == 1
    monkeypatch.setattr(os, "sched_getaffinity", lambda pid: set(range(64)))
    assert min(2, processori_disponibili()) == 2


def test_con_un_solo_processore_g14_accetta_un_thread_e_ne_respinge_due(tmp_path, monkeypatch):
    """
    **Obiettivo**: Verificare che con un solo processore utilizzabile G14 sia
    superato con ``run.threads`` 1 e respinga ``run.threads`` 2 con
    ``E-S0-14``.

    **Razionale scientifico e sistemistico**: E' il comportamento corretto
    della pipeline, quello a cui i test si adattano: chiedere piu' thread dei
    processori utilizzabili sovraccaricherebbe la macchina senza accelerare
    nulla.
    """
    monkeypatch.setattr(os, "sched_getaffinity", lambda pid: {0})
    campioni = [
        Campione("ERX3000001", "NOD1D4.L1"),
        Campione("ERX3000003", "POS.P1.1", materiale=POSITIVO, posizione="Not Applicable"),
        Campione("ERX3000004", "BLANK.P1.1", materiale=NEGATIVO, posizione="Not Applicable"),
    ]
    scenario = crea_scenario(tmp_path, campioni, con_letture=True)
    assert scenario.config.run.threads == 1
    assert esegui_gate("G14", Contesto(scenario.config)).superato
    dati = scenario.config.model_dump(mode="python")
    dati["run"]["threads"] = 2
    esito = esegui_gate("G14", Contesto(valida(dati)))
    assert not esito.superato and {v.codice for v in esito.violazioni} == {"E-S0-14"}


def test_il_file_di_blocco_contiene_i_pacchetti_delle_analisi_ecologiche():
    """
    **Obiettivo**: Verificare che ``renv.lock`` registri DESeq2, vegan e
    randomForest con la loro chiusura (locfit), accanto ai pacchetti della
    pipeline, con R 4.5.2, Bioconductor 3.21 e dada2 1.36.0.1 con la
    correzione; e che i due script del container elenchino gli stessi
    pacchetti delle analisi ecologiche.

    **Razionale scientifico e sistemistico**: Le analisi ecologiche girano
    nello stesso ambiente bloccato della pipeline: un pacchetto installato e
    non registrato nel file di blocco non sarebbe verificato alla costruzione.
    """
    blocco = json.loads((RADICE / "renv.lock").read_text(encoding="utf-8"))
    pacchetti = blocco["Packages"]
    assert blocco["R"]["Version"] == "4.5.2" and blocco["Bioconductor"]["Version"] == "3.21"
    for nome in (*ECOLOGIA, "locfit", "dada2", "phyloseq", "DECIPHER", "phangorn", "decontam"):
        assert nome in pacchetti, nome
    assert pacchetti["DESeq2"]["Source"] == "Bioconductor"
    assert pacchetti["dada2"]["Version"] == "1.36.0.1" and pacchetti["dada2"]["Patch"]["Base"] == "1.36.0"
    elenco = 'ecologia <- c("DESeq2", "vegan", "randomForest")'
    for script in ("install_r_packages.R", "snapshot_renv_lock.R"):
        assert elenco in (RADICE / "container" / script).read_text(encoding="utf-8"), script


def test_i_pacchetti_delle_analisi_ecologiche_si_caricano_alle_versioni_del_blocco(tmp_path):
    """
    **Obiettivo**: Verificare che nell'ambiente R dell'immagine DESeq2, vegan e
    randomForest si carichino, e che le loro versioni siano quelle di
    ``renv.lock``.

    **Razionale scientifico e sistemistico**: Le settimane delle analisi
    ecologiche non devono richiedere un'immagine nuova: i pacchetti ci sono
    gia', alle versioni registrate.
    """
    motivo = motivo_pacchetti_r_assenti(*ECOLOGIA)
    if motivo is not None:
        if os.environ.get("AMPLICON16S_RICHIEDI_BIOC") == "1":
            pytest.fail(f"i pacchetti delle analisi ecologiche sono richiesti in questo ambiente: {motivo}")
        pytest.skip(motivo)
    uscita = tmp_path / "versioni.json"
    codice = (
        f'p <- c({", ".join(repr(n).replace(chr(39), chr(34)) for n in ECOLOGIA)}); '
        "for (n in p) loadNamespace(n); "
        "v <- vapply(p, function(n) as.character(packageVersion(n)), character(1)); "
        f"jsonlite::write_json(as.list(v), {json.dumps(str(uscita))}, auto_unbox = TRUE)"
    )
    subprocess.run([str(trova_rscript()), "--vanilla", "-e", codice], check=True, capture_output=True)
    installate = json.loads(uscita.read_text(encoding="utf-8"))
    blocco = json.loads((RADICE / "renv.lock").read_text(encoding="utf-8"))["Packages"]
    for nome in ECOLOGIA:
        # R scrive 4.7-1.2 come 4.7.1.2: si confrontano i numeri, non i separatori.
        assert re.split(r"[.-]", installate[nome]) == re.split(r"[.-]", blocco[nome]["Version"]), nome


@pytest.mark.skipif(
    not (RADICE / "dati" / "osd734").is_dir(),
    reason="la cartella dati/ non e' presente (nell'immagine non viene copiata)",
)
def test_la_configurazione_di_osd734_dichiara_il_digest_di_registro():
    """
    **Obiettivo**: Verificare che ``run.container`` di ``config_osd734.yaml``
    sia un riferimento a un'immagine in un registro ancorato per digest, non il
    segnaposto ne' un identificativo locale, e che il README principale, quello
    della cartella dei dati e ``test.txt`` indichino quella stessa immagine.

    **Razionale scientifico e sistemistico**: Il digest di registro identifica
    gli stessi byte per chiunque scarichi l'immagine; il nome puo' cambiare, il
    digest no. Un riferimento diverso fra configurazione e documentazione
    farebbe eseguire un'immagine e dichiararne un'altra.
    """
    config = yaml.safe_load(
        (RADICE / "dati" / "osd734" / "config_osd734.yaml").read_text(encoding="utf-8")
    )
    immagine = config["run"]["container"]
    assert DIGEST_DI_REGISTRO.match(immagine), immagine
    assert not immagine.endswith("0" * 64)
    valida({**config, "run": {**config["run"]}})
    for documento in ("README.md", "dati/osd734/README.md", "test.txt"):
        assert immagine in (RADICE / documento).read_text(encoding="utf-8"), documento
