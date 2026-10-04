r"""Suite di verifica del sistema di logging strutturato JSONL, della rotazione e dell'esecuzione integrata.

1. Inquadramento nel Piano Operativo
------------------------------------
Settimana 5 (W5), Fase F1 (logging strutturato, tracciabilità e gestione
degli artefatti).

2. Moduli sorgente coperti
--------------------------
* ``src/amplicon16s/logging/logger.py``
* ``src/amplicon16s/io_layer/artifacts.py``
* ``src/amplicon16s/config/resolve.py``

3. Cosa valuta questo file
--------------------------
- doppio canale: console sintetica per l'operatore e file strutturato JSON
  Lines in ``99_logs/``, con la sola console quando manca la radice di output
  e nessun accumulo di uscite a configurazioni ripetute;
- rotazione del file di log (``RotatingFileHandler``) con dimensione massima e
  numero finito di copie;
- ogni riga del file è un oggetto JSON; i campi aggiuntivi dell'evento
  (fase, codice, metriche) diventano campi interrogabili, anche per codice;
- serializzazione resiliente (``default=str``): un valore non serializzabile
  non fa perdere l'evento; le eccezioni finiscono nel log con la categoria di
  gestione dell'errore;
- livello del file più verboso di quello della console;
- esecuzione fittizia che combina ``AlberoOutput``, ``scrivi_risolta``
  (``00_config/resolved.yaml``) e il logger sulla stessa radice.

4. Comandi Bash e scenari di esecuzione
---------------------------------------
    ``<immagine>`` e' l'immagine del container della pipeline; quella corrente
    e' indicata in ``test.txt``, sezione 1.3.

    1. Modalità locale standard (R di base con jsonlite, senza Bioconductor né
       dati reali):
       pytest tests/test_w05_logging.py -v

    2. Modalità container Docker standard (sottoinsieme ridotto con
       R/Bioconductor):
       docker run --rm \
         -e PYTHONPATH=/app/src \
         -e AMPLICON16S_R_DIR=/app/R \
         -v "$(pwd)":/app \
         -w /app \
         <immagine> \
         pytest -o cache_dir=/tmp/.pytest_cache tests/test_w05_logging.py -v

    3. Modalità container Docker completa (con i 2.4 GB di dati reali OSD-734;
       la configurazione e i percorsi che contiene devono stare nella cartella
       montata):
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
         pytest -o cache_dir=/tmp/.pytest_cache tests/test_w05_logging.py -v

5. Risultato atteso
-------------------
14 test totali:
- 14 passed in ambiente locale standard (~0.1s);
- 14 passed nel container Docker standard sul sottoinsieme ridotto (~0.1s);
- 14 passed nel container Docker con i dati reali OSD-734 (~0.1s).

6. Razionale scientifico e sistemistico
---------------------------------------
- In un'elaborazione di 960 campioni attraverso più fasi, il log strutturato è
  lo strumento con cui si ricostruisce che cosa è successo e perché: deve
  poter essere interrogato per fase e per codice.
- Un errore del formattatore farebbe perdere il log proprio nel momento in cui
  si registra un guasto; la rotazione impedisce che un'esecuzione lunga o
  verbosa saturi il disco.
"""

from __future__ import annotations

import json
import logging
import logging.handlers
from pathlib import Path

import pytest
import yaml

from amplicon16s.config.resolve import scrivi_risolta
from amplicon16s.errors.exceptions import errore
from amplicon16s.gates.g01_g15 import esegui_g15
from amplicon16s.io_layer.artifacts import AlberoOutput, Fase
from amplicon16s.logging.logger import (
    MAX_BYTE,
    NOME_FILE_LOG,
    chiudi,
    configura,
    ottieni,
    registra_errore,
)


ESEMPIO = Path(__file__).resolve().parents[1] / "config" / "config.example.yaml"


@pytest.fixture(autouse=True)
def uscite_pulite():
    """Ogni test parte senza uscite di log residue e le chiude alla fine."""
    chiudi()
    yield
    chiudi()


def _righe(percorso) -> list[dict]:
    """Gli eventi di un file JSONL, uno per riga non vuota."""
    testo = percorso.read_text(encoding="utf-8")
    return [json.loads(r) for r in testo.splitlines() if r.strip()]


# --------------------------------------------------------------------------- #
# Le due uscite                                                                #
# --------------------------------------------------------------------------- #


def test_il_log_strutturato_nasce_sotto_99_logs(tmp_path):
    """
    **Obiettivo**: Verificare che ``configura(tmp_path)`` crei la sottocartella
    ``99_logs/`` e restituisca il percorso ``99_logs/pipeline.jsonl``.

    **Razionale scientifico e sistemistico**: Garantisce che la traccia di audit
    della pipeline risieda sempre nella collocazione canonica ``99_logs/``
    dell'albero degli artefatti.
    """
    percorso = configura(tmp_path)
    assert percorso == tmp_path / Fase.LOGS.value / NOME_FILE_LOG
    assert percorso.parent.is_dir()


def test_senza_radice_resta_la_sola_console(tmp_path):
    """
    **Obiettivo**: Verificare che ``configura(None)`` attivi unicamente lo
    ``StreamHandler`` di console senza creare file o cartelle su disco.

    **Razionale scientifico e sistemistico**: Permette ai comandi preliminari o
    alle validazioni in memoria (es. prima che ``out_root`` sia creata) di
    emettere messaggi diagnostici su console senza sporcare il filesystem.
    """
    assert configura(None) is None
    tipi = [type(h) for h in ottieni().handlers]
    assert tipi == [logging.StreamHandler]
    assert not list(tmp_path.iterdir())


def test_le_uscite_sono_due(tmp_path):
    """
    **Obiettivo**: Verificare che dopo ``configura(tmp_path)`` il logger radice
    di ``amplicon16s`` possieda esattamente 2 handler (console + file JSONL).

    **Razionale scientifico e sistemistico**: Accerta l'attivazione simultanea del
    canale leggibile dall'uomo e del canale strutturato per le macchine.
    """
    configura(tmp_path)
    assert len(ottieni().handlers) == 2


def test_configurare_piu_volte_non_accumula_uscite(tmp_path):
    """
    **Obiettivo**: Verificare che chiamate ripetute a ``configura(tmp_path)``
    chiudano gli handler precedenti mantenendo il conteggio fisso a 2.

    **Razionale scientifico e sistemistico**: Previene il classico bug di
    duplicazione dei messaggi nel modulo ``logging`` di Python e la perdita di
    file descriptor aperti durante invocazioni multiple nello stesso processo.
    """
    configura(tmp_path)
    configura(tmp_path)
    configura(tmp_path)
    assert len(ottieni().handlers) == 2


def test_il_file_ruota(tmp_path):
    """
    **Obiettivo**: Verificare che l'handler su file sia un'istanza di
    ``RotatingFileHandler`` configurata con ``maxBytes == MAX_BYTE`` e ``backupCount > 0``.

    **Razionale scientifico e sistemistico**: Protegge il filesystem dall'esaurimento
    dello spazio disco qualora la pipeline produca un volume elevato di log
    diagnostici sui 960 campioni.
    """
    configura(tmp_path)
    rotanti = [
        h for h in ottieni().handlers
        if isinstance(h, logging.handlers.RotatingFileHandler)
    ]
    assert len(rotanti) == 1
    assert rotanti[0].maxBytes == MAX_BYTE
    assert rotanti[0].backupCount > 0


def test_la_rotazione_avviene_davvero(tmp_path):
    """
    **Obiettivo**: Verificare sperimentalmente che superando ``max_byte=500``
    vengano generati su disco sia ``pipeline.jsonl`` sia l'archivio ruotato ``pipeline.jsonl.1``.

    **Razionale scientifico e sistemistico**: Certifica che il meccanismo di
    rollover dei file di log funzioni effettivamente su filesystem reale senza
    interrompere il flusso di scrittura.
    """
    percorso = configura(tmp_path, max_byte=500, rotazioni=2)
    log = ottieni("prova")
    for i in range(200):
        log.info("evento di prova ripetuto per far crescere il file", extra={"i": i})
    chiudi()

    ruotati = sorted(p.name for p in percorso.parent.iterdir())
    assert NOME_FILE_LOG in ruotati
    assert f"{NOME_FILE_LOG}.1" in ruotati


# --------------------------------------------------------------------------- #
# Il file è interrogabile                                                      #
# --------------------------------------------------------------------------- #


def test_ogni_riga_e_un_oggetto_json(tmp_path):
    """
    **Obiettivo**: Verificare che ogni riga di ``pipeline.jsonl`` sia un oggetto
    JSON autonomo e valido dotato dei campi ``messaggio`` e ``livello``.

    **Razionale scientifico e sistemistico**: Il formato JSON Lines (un documento
    JSON per riga) consente la lettura in streaming e il parsing riga per riga
    anche se il processo dovesse interrompersi bruscamente a metà esecuzione.
    """
    percorso = configura(tmp_path)
    log = ottieni("fase")
    log.info("primo")
    log.warning("secondo")
    chiudi()

    righe = _righe(percorso)
    assert len(righe) == 2
    assert [r["messaggio"] for r in righe] == ["primo", "secondo"]
    assert [r["livello"] for r in righe] == ["INFO", "WARNING"]


def test_i_campi_aggiunti_diventano_campi_dell_evento(tmp_path):
    """
    **Obiettivo**: Verificare che gli attributi passati nel dizionario ``extra``
    (``fase="S2"``, ``conservate=0.83``) vengano promossi a chiavi di primo
    livello nell'oggetto JSON insieme a ``origine`` e ``istante``.

    **Razionale scientifico e sistemistico**: Permette di registrare metriche
    quantitative di fase (es. frazione di letture conservate dopo ``filterAndTrim``)
    direttamente interrogabili senza espressioni regolari sul testo.
    """
    percorso = configura(tmp_path)
    ottieni("s2").info("filtro concluso", extra={"fase": "S2", "conservate": 0.83})
    chiudi()

    evento = _righe(percorso)[0]
    assert evento["fase"] == "S2"
    assert evento["conservate"] == 0.83
    assert evento["origine"] == "amplicon16s.s2"
    assert "istante" in evento


def test_gli_eventi_si_filtrano_per_codice(tmp_path):
    """
    **Obiettivo**: Verificare che gli errori registrati tramite ``registra_errore``
    siano filtrabili esaminando ``r.get("codice") == "E-S2-03"``.

    **Razionale scientifico e sistemistico**: Consente al generatore di report e
    agli script di audit di estrarre istantaneamente tutti gli episodi di
    memoria esaurita o degradazione avvenuti su specifici lotti.
    """
    percorso = configura(tmp_path)
    log = ottieni("s2")
    registra_errore(log, errore("E-S2-03", "lotto 1"))
    log.info("proseguo")
    registra_errore(log, errore("E-S2-03", "lotto 4"))
    chiudi()

    del_codice = [r for r in _righe(percorso) if r.get("codice") == "E-S2-03"]
    assert len(del_codice) == 2
    assert [r["dettaglio"] for r in del_codice] == ["lotto 1", "lotto 4"]


def test_l_errore_registrato_porta_con_se_la_gestione(tmp_path):
    """
    **Obiettivo**: Verificare che ``registra_errore`` includa nell'evento JSON
    ``codice``, ``fase``, ``categoria``, i parametri di contesto (``sequenze=12000``)
    e l'``azione`` prescrittiva del catalogo.

    **Razionale scientifico e sistemistico**: Garantisce che ogni evento d'errore
    nel file ``pipeline.jsonl`` sia autosufficiente e contenga già la politica di
    gestione e l'istruzione operativa per il ricercatore.
    """
    percorso = configura(tmp_path)
    registra_errore(ottieni("s9"), errore("E-S9-01", "12000 sequenze", sequenze=12000))
    chiudi()

    evento = _righe(percorso)[0]
    assert evento["codice"] == "E-S9-01"
    assert evento["fase"] == "S9"
    assert evento["categoria"] == "revisione_umana"
    assert evento["sequenze"] == 12000
    assert evento["azione"]


def test_un_valore_non_serializzabile_non_fa_perdere_l_evento(tmp_path):
    """
    **Obiettivo**: Verificare che passare oggetti ``Path`` o ``Enum`` dentro
    ``extra`` non sollevi ``TypeError`` ma li converta in stringa nel JSON finale.

    **Razionale scientifico e sistemistico**: Se ``json.dumps`` fallisse su un
    oggetto ``Path`` passato nel contesto di un errore, l'evento diagnostico
    andrebbe perso proprio nel momento critico del guasto; il fallback
    ``default=str`` rende il logger immune a tipi non nativamente JSON.
    """
    percorso = configura(tmp_path)
    ottieni("s0").info("percorso", extra={"dove": tmp_path, "fase_enum": Fase.CONFIG})
    chiudi()

    evento = _righe(percorso)[0]
    assert evento["dove"] == str(tmp_path)
    assert evento["fase_enum"] == Fase.CONFIG.value


def test_l_eccezione_finisce_nel_log(tmp_path):
    """
    **Obiettivo**: Verificare che ``log.exception(...)`` catturi lo stacktrace
    corrente nel campo ``"eccezione"`` dell'oggetto JSON.

    **Razionale scientifico e sistemistico**: Conserva nel file JSONL la traccia
    completa delle chiamate Python in caso di errore imprevisto senza corrompere
    la struttura a riga singola del file ``pipeline.jsonl``.
    """
    percorso = configura(tmp_path)
    try:
        raise errore("E-S14-01")
    except Exception:
        ottieni("s14").exception("validazione finale fallita")
    chiudi()

    assert "E-S14-01" in _righe(percorso)[0]["eccezione"]


def test_il_livello_del_file_e_piu_verboso_della_console(tmp_path):
    """
    **Obiettivo**: Verificare che anche impostando ``livello_console=logging.WARNING``
    i messaggi ``DEBUG`` vengano comunque scritti in ``pipeline.jsonl``.

    **Razionale scientifico e sistemistico**: Mantiene pulito il terminale
    dell'operatore mostrando solo avvisi ed errori, ma preserva sul disco tutti
    i dettagli diagnostici di livello ``DEBUG`` indispensabili per il post-mortem.
    """
    percorso = configura(tmp_path, livello_console=logging.WARNING)
    ottieni("x").debug("dettaglio diagnostico")
    chiudi()
    assert [r["messaggio"] for r in _righe(percorso)] == ["dettaglio diagnostico"]


# --------------------------------------------------------------------------- #
# Esecuzione fittizia: albero completo, configurazione risolta, log            #
# --------------------------------------------------------------------------- #


def test_esecuzione_fittizia_produce_albero_configurazione_e_log(tmp_path):
    """
    **Obiettivo**: Verificare l'integrazione end-to-end dei tre servizi della
    Fase F1: creazione delle 14 cartelle di ``AlberoOutput``, scrittura e firma
    di ``00_config/resolved.yaml``, e registrazione
    degli eventi in ``99_logs/pipeline.jsonl``.

    **Razionale scientifico e sistemistico**: Simula il ciclo di vita completo
    dell'infrastruttura di base (Settimane W4-W5), dimostrando che
    configurazione risolta, manifesti SHA-256 degli artefatti e log JSONL
    coesistono senza conflitti nella medesima directory di esecuzione.
    """
    radice = tmp_path / "esecuzione"

    albero = AlberoOutput(radice)
    albero.crea()
    percorso_log = configura(radice)
    log = ottieni("runner")

    dati = yaml.safe_load(ESEMPIO.read_text(encoding="utf-8"))
    risolta = esegui_g15(dati)
    scrivi_risolta(risolta, radice)

    log.info(
        "configurazione risolta",
        extra={"digest": risolta.digest, "campioni_biologici": 803},
    )
    registra_errore(log, errore("E-S1-01", "scarto di 14 basi", scarto=14))
    albero.scrivi_testo(Fase.INPUT_VALIDATION, "riepilogo.tsv", "campioni\t803\n")
    chiudi()

    # L'albero completo
    cartelle = sorted(p.name for p in radice.iterdir() if p.is_dir())
    assert cartelle == sorted(f.value for f in Fase)
    assert len(cartelle) == 14

    # La configurazione risolta, registrata come artefatto
    assert albero.fase_completa(Fase.CONFIG)
    assert "resolved.yaml" in albero.manifesto(Fase.CONFIG)

    # Il log strutturato
    eventi = _righe(percorso_log)
    assert any(e.get("digest") == risolta.digest for e in eventi)
    degradazioni = [e for e in eventi if e.get("categoria") == "degradazione_automatica"]
    assert [e["codice"] for e in degradazioni] == ["E-S1-01"]

    # Gli artefatti sono verificabili
    assert albero.non_integri(Fase.INPUT_VALIDATION) == ()
