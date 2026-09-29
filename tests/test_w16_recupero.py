r"""Suite di test del recupero della settimana 16: riproducibilità degli artefatti, file del ponte per fase, tabella vuota.

1. Inquadramento nel Piano Operativo
------------------------------------
Settimana 16 (W16), Fase F4 (recupero sulle fasi S0-S7 prima della fase
tassonomica: artefatti identici fra esecuzioni, traccia del ponte in una
cartella condivisa, arresto su una tabella senza varianti).

2. Moduli sorgente coperti
--------------------------
* ``src/amplicon16s/gates/registry.py`` (``EsitoGate.come_voce``)
* ``src/amplicon16s/steps/s00_validate.py`` (``gates.json`` e log delle durate)
* ``src/amplicon16s/rbridge/payload.py`` (``nomi_del_contratto``)
* ``src/amplicon16s/rbridge/runner.py`` (``esegui_script``, argomento ``passo``)
* ``src/amplicon16s/steps/s07_asv_length.py`` (``E-S7-01``)
* ``src/amplicon16s/errors/catalog.py`` (voce ``E-S7-01``)
* ``R/07_asv_length.R``

3. Cosa valuta questo file
--------------------------
- ``gates.json`` e le metriche del manifesto di S0 non registrano durate; le
  durate di ciascun gate finiscono nel log strutturato; due esecuzioni di S0
  sugli stessi ingressi danno gli stessi byte;
- una seconda esecuzione da zero di S0-S7 sulla versione ridotta dà, fase per
  fase, gli stessi byte della prima (sostituisce le due verifiche separate di
  S3 e di S4-S7 delle settimane 13 e 14, e vi aggiunge S0, S1 e S2);
- i file del contratto del ponte portano il nome della fase: S6 e S7, che
  condividono ``07_chimera/``, lasciano ciascuna la propria richiesta e la
  propria dichiarazione d'esito; nessuna coppia di fasi che condivide una
  cartella condivide un file del ponte o un manifesto di fase;
- S7 si ferma con ``E-S7-01``, di revisione umana, se la tabella non ha più
  varianti: sia quando la tabella di S6 è già vuota, sia quando tutte le
  varianti cadono fuori dall'intervallo di lunghezza.

4. Comandi Bash e scenari di esecuzione
---------------------------------------
    1. Modalità locale standard (R di base con jsonlite, senza Bioconductor né
       dati reali):
       pytest tests/test_w16_recupero.py -v

    2. Modalità container Docker standard (sottoinsieme ridotto con
       R/Bioconductor):
       docker run --rm \
         -e PYTHONPATH=/app/src \
         -e AMPLICON16S_R_DIR=/app/R \
         -v "$(pwd)":/app \
         -w /app \
         amplicon16s:dev \
         pytest -o cache_dir=/tmp/.pytest_cache tests/test_w16_recupero.py -v

    3. Modalità container Docker completa (con i 2.4 GB di dati reali OSD-734;
       la configurazione e i percorsi che contiene devono stare nella cartella
       montata):
       docker run --rm \
         -e PYTHONPATH=/app/src \
         -e AMPLICON16S_R_DIR=/app/R \
         -e AMPLICON16S_CONFIG_DATI_REALI="$HOME/ASI/config_osd734.yaml" \
         -v "$(pwd)":/app \
         -v "$HOME/ASI":"$HOME/ASI" \
         -w /app \
         amplicon16s:dev \
         pytest -o cache_dir=/tmp/.pytest_cache tests/test_w16_recupero.py -v

5. Risultato atteso
-------------------
9 test totali:
- 5 passed, 4 skipped in ambiente locale standard (~0.1s): i 4 test che eseguono
  le fasi richiedono dada2, ShortRead e ggplot2;
- 9 passed nel container Docker standard sul sottoinsieme ridotto (~67s);
- 9 passed nel container Docker con i dati reali OSD-734 (~67s): il modulo non ha
  test sui dati reali.

6. Razionale scientifico e sistemistico
---------------------------------------
- La verifica finale di riproducibilità confronta i checksum degli artefatti
  di due esecuzioni: un artefatto che registra una durata cambia a ogni
  esecuzione e la rende impossibile. Le durate descrivono l'esecuzione, non il
  risultato, e il loro posto è il log.
- I file del ponte non sono artefatti, ma sono l'unica traccia di come è
  andata l'invocazione di R: in una cartella condivisa l'ultima fase non deve
  cancellare quella della precedente.
- Una tabella senza varianti farebbe fallire le fasi successive con un errore
  che non ne indica la causa; fermarsi in S7 con un codice del catalogo dice
  dove guardare (``lunghezze.tsv``) e che un nuovo tentativo non servirebbe.
"""

from __future__ import annotations

import hashlib
import json
import os
import shutil
import subprocess
from pathlib import Path

import pytest
from conftest import NEGATIVO, POSITIVO, Campione, copia_esecuzione, crea_scenario
from sottoinsieme import config_ridotta, motivo_pacchetti_r_assenti

from amplicon16s.errors.catalog import Categoria
from amplicon16s.errors.exceptions import ErrorePipeline
from amplicon16s.io_layer.artifacts import Fase, nome_manifesto_passo
from amplicon16s.logging.logger import NOME_FILE_LOG, chiudi, configura
from amplicon16s.rbridge.payload import PREFISSO, nomi_del_contratto
from amplicon16s.rbridge.runner import trova_rscript
from amplicon16s.runner.executor import Conclusione, Esecutore
from amplicon16s.runner.graph import GRAFO, Passo
from amplicon16s.runner.project import ProjectRun
from amplicon16s.steps.s00_validate import (
    NOME_CROSSWALK,
    NOME_ESITI,
    NOME_INVENTARIO,
    NOME_LETTURE,
    esegui_s0,
)


@pytest.fixture(autouse=True)
def uscite_pulite():
    """Chiude le uscite del log prima e dopo ogni test, perché nessun handler resti
    aperto sulla cartella temporanea.
    """
    chiudi()
    yield
    chiudi()


_MOTIVO_ASSENTI = motivo_pacchetti_r_assenti("dada2", "ggplot2", "ShortRead", "jsonlite")


@pytest.fixture
def dada2():
    """Richiede R con dada2: salta senza, ma in CI fallisce."""
    if _MOTIVO_ASSENTI is not None:
        if os.environ.get("AMPLICON16S_RICHIEDI_BIOC") == "1":
            pytest.fail(f"Bioconductor e' richiesto in questo ambiente: {_MOTIVO_ASSENTI}")
        pytest.skip(_MOTIVO_ASSENTI)


def _campioni() -> list[Campione]:
    """Quattro campioni sulla piastra 1: due biologici, un positivo, un negativo."""
    return [
        Campione("ERX3000001", "NOD1D4.L1", piastra="1"),
        Campione("ERX3000002", "NOD1D4.L2", piastra="1"),
        Campione("ERX3000003", "POS.P1.1", materiale=POSITIVO,
                 posizione="Not Applicable", piastra="1"),
        Campione("ERX3000004", "BLANK.P1.1", materiale=NEGATIVO,
                 posizione="Not Applicable", piastra="1"),
    ]


def _chiavi(valore) -> set[str]:
    """Tutte le chiavi di un documento JSON, a qualunque profondità."""
    if isinstance(valore, dict):
        return set(valore) | {k for v in valore.values() for k in _chiavi(v)}
    if isinstance(valore, list):
        return {k for v in valore for k in _chiavi(v)}
    return set()


def _impronte(cartella: Path) -> dict[str, str]:
    """L'MD5 di ogni file della cartella, esclusi i file del ponte e i manifesti."""
    return {
        p.name: hashlib.md5(p.read_bytes()).hexdigest()
        for p in sorted(cartella.iterdir())
        if not p.name.startswith((PREFISSO, "manifest"))
    }


#: Le cartelle delle fasi realizzate, da S0 a S7.
_CARTELLE_S0_S7 = (
    Fase.INPUT_VALIDATION, Fase.QC_PROFILES, Fase.FILTERED, Fase.ERROR_MODELS,
    Fase.ASV_INFERENCE, Fase.SEQTAB, Fase.CHIMERA,
)


# --------------------------------------------------------------------------- #
# 1a. Le durate fuori dagli artefatti                                          #
# --------------------------------------------------------------------------- #


def test_gates_json_e_il_manifesto_di_s0_non_registrano_durate(tmp_path):
    """
    **Obiettivo**: Verificare che dopo S0 né ``gates.json`` né le metriche del
    manifesto di S0 contengano il campo ``secondi``, a nessuna profondità.

    **Razionale scientifico e sistemistico**: ``gates.json`` era l'unico
    artefatto su 996 a differire dall'esecuzione di riferimento, per le durate
    dei gate: un artefatto che cambia a ogni esecuzione rende impossibile il
    confronto dei checksum da un clone pulito.
    """
    scenario = crea_scenario(tmp_path, _campioni(), con_letture=True)
    esegui_s0(scenario.config)
    cartella = Path(scenario.config.io.out_root) / Fase.INPUT_VALIDATION.value
    esiti = json.loads((cartella / NOME_ESITI).read_text(encoding="utf-8"))
    assert esiti["superata"] is True
    assert "secondi" not in _chiavi(esiti)
    manifesto = json.loads((cartella / nome_manifesto_passo("S0")).read_text(encoding="utf-8"))
    assert "secondi" not in _chiavi(manifesto["metriche"])


def test_due_esecuzioni_di_s0_danno_gli_stessi_byte(tmp_path):
    """
    **Obiettivo**: Verificare che due esecuzioni di S0 sugli stessi ingressi
    scrivano ``gates.json``, ``crosswalk.tsv``, ``inventario.json`` e
    ``letture_ispezionate.tsv`` identici byte per byte.

    **Razionale scientifico e sistemistico**: E' la proprieta' che la verifica
    finale di riproducibilita' richiede a ogni artefatto; il test gira senza R,
    quindi anche nel job leggero della CI.
    """
    scenario = crea_scenario(tmp_path, _campioni(), con_letture=True)
    cartella = Path(scenario.config.io.out_root) / Fase.INPUT_VALIDATION.value
    nomi = (NOME_ESITI, NOME_CROSSWALK, NOME_INVENTARIO, NOME_LETTURE)
    esegui_s0(scenario.config)
    prima = {n: (cartella / n).read_bytes() for n in nomi}
    shutil.rmtree(scenario.config.io.out_root)
    esegui_s0(scenario.config)
    assert {n: (cartella / n).read_bytes() for n in nomi} == prima


def test_le_durate_dei_gate_finiscono_nel_log(tmp_path):
    """
    **Obiettivo**: Verificare che per ogni gate eseguito il log strutturato
    registri un evento di S0 con il nome del gate e la durata in secondi.

    **Razionale scientifico e sistemistico**: Le durate restano disponibili per
    capire dove si spende il tempo, ma nel posto che descrive l'esecuzione e
    non nel risultato.
    """
    scenario = crea_scenario(tmp_path, _campioni(), con_letture=True)
    percorso = configura(scenario.config.io.out_root)
    esegui_s0(scenario.config)
    chiudi()
    assert percorso is not None and percorso.name == NOME_FILE_LOG
    eventi = [json.loads(r) for r in percorso.read_text(encoding="utf-8").splitlines() if r]
    durate = {e["gate"]: e["secondi"] for e in eventi if e.get("fase") == "S0" and "gate" in e
              and "secondi" in e}
    assert set(durate) == {f"G{i:02d}" for i in range(1, 16)}
    assert all(isinstance(s, float) and s >= 0 for s in durate.values())


def test_una_seconda_esecuzione_da_zero_di_s0_s7_da_gli_stessi_byte(dada2, catena_calcolata, tmp_path):
    """
    **Obiettivo**: Verificare che una seconda esecuzione da zero di S0-S7 sulla
    versione ridotta, in un'altra cartella di output, produca in ciascuna delle
    sette cartelle di fase gli stessi byte della prima, ``gates.json``, modelli
    d'errore e grafici compresi.

    **Razionale scientifico e sistemistico**: E' la verifica di riproducibilita'
    che si fara' sul dataset completo da un clone pulito, qui sulla versione
    ridotta. Una sola esecuzione condivisa sostituisce le due riesecuzioni
    parziali di S3 e di S4-S7 e copre anche S0, S1 e S2.
    """
    run, _ = catena_calcolata
    seconda = Esecutore(ProjectRun(config_ridotta(tmp_path)), fino_a=Passo.S7)
    esito = seconda.esegui()
    assert esito.conclusione is Conclusione.COMPLETATA
    for fase in _CARTELLE_S0_S7:
        assert _impronte(seconda.run.albero.cartella(fase)) == _impronte(run.albero.cartella(fase)), fase


# --------------------------------------------------------------------------- #
# 1b. I file del ponte distinti per fase                                       #
# --------------------------------------------------------------------------- #


def test_i_file_del_ponte_portano_il_nome_della_fase():
    """
    **Obiettivo**: Verificare che ``nomi_del_contratto`` dia a fasi diverse
    nomi diversi, con il prefisso comune ``rbridge_``, e rifiuti una fase vuota
    o con una barra.

    **Razionale scientifico e sistemistico**: Il nome e' l'unica cosa che
    separa i file di due fasi nella stessa cartella; un nome che uscisse dalla
    cartella di fase non sarebbe accettabile.
    """
    s6, s7 = nomi_del_contratto("S6"), nomi_del_contratto("S7")
    assert s6 == ("rbridge_richiesta_S6.json", "rbridge_esito_S6.json")
    assert not set(s6) & set(s7)
    for nome in (*s6, *s7):
        assert nome.startswith(PREFISSO)
    for non_valida in ("", "S6/../S7"):
        with pytest.raises(ValueError):
            nomi_del_contratto(non_valida)


def test_le_fasi_che_condividono_una_cartella_non_condividono_file_di_servizio():
    """
    **Obiettivo**: Verificare, per ogni cartella del grafo in cui scrivono piu'
    fasi (oggi ``07_chimera``, ``11_controls`` e ``12_final``), che i file del
    ponte e i manifesti di fase abbiano nomi distinti per ciascuna fase.

    **Razionale scientifico e sistemistico**: Lo stesso problema di S6 e S7 si
    ripresenterebbe con S11 e S12, e con S13 e S14: il controllo lo esclude
    anche per le fasi non ancora realizzate.
    """
    per_cartella: dict[Fase, list[Passo]] = {}
    for passo in GRAFO.ordine():
        per_cartella.setdefault(GRAFO.nodo(passo).cartella, []).append(passo)
    condivise = {c: p for c, p in per_cartella.items() if len(p) > 1}
    assert {c.value for c in condivise} == {"07_chimera", "11_controls", "12_final"}
    for cartella, passi in condivise.items():
        nomi = [n for p in passi for n in (*nomi_del_contratto(str(p)), nome_manifesto_passo(str(p)))]
        assert len(nomi) == len(set(nomi)), cartella


def test_s6_e_s7_lasciano_ciascuna_la_propria_traccia_del_ponte(dada2, catena_calcolata):
    """
    **Obiettivo**: Verificare che dopo S6 e S7 la cartella ``07_chimera``
    contenga la richiesta e la dichiarazione d'esito di entrambe le fasi, con
    lo stesso identificativo d'invocazione fra richiesta ed esito di ciascuna,
    e nessun file del ponte senza il nome della fase.

    **Razionale scientifico e sistemistico**: Prima S7 sovrascriveva i file di
    S6, e con loro si perdeva la traccia di come era andata S6.
    """
    run, _ = catena_calcolata
    cartella = run.albero.cartella(Fase.CHIMERA)
    invocazioni = set()
    for passo, script in (("S6", "06_chimera.R"), ("S7", "07_asv_length.R")):
        nome_richiesta, nome_esito = nomi_del_contratto(passo)
        richiesta = json.loads((cartella / nome_richiesta).read_text(encoding="utf-8"))
        esito = json.loads((cartella / nome_esito).read_text(encoding="utf-8"))
        assert esito["stato"] == "riuscito", passo
        assert esito["invocazione"] == richiesta["invocazione"], passo
        assert richiesta["esito"] == str(cartella / nome_esito), passo
        invocazioni.add(richiesta["invocazione"])
    assert len(invocazioni) == 2
    del_ponte = sorted(p.name for p in cartella.iterdir() if p.name.startswith(PREFISSO))
    assert del_ponte == sorted([*nomi_del_contratto("S6"), *nomi_del_contratto("S7")])


# --------------------------------------------------------------------------- #
# 1c. Una tabella senza varianti                                               #
# --------------------------------------------------------------------------- #


_RENDI_VUOTA = {
    # Una tabella di S6 gia' vuota: tutte le varianti rimosse a monte.
    "tabella_di_s6_vuota": "t <- t[, 0, drop = FALSE]",
    # Varianti tutte di 100 basi, fuori dall'intervallo 137 +/- asv.len_tol.
    "varianti_fuori_intervallo": "colnames(t) <- substr(colnames(t), 1, 100)",
}


@pytest.mark.parametrize("caso", sorted(_RENDI_VUOTA))
def test_s7_senza_varianti_si_ferma_con_e_s7_01(dada2, catena_calcolata, tmp_path, caso):
    """
    **Obiettivo**: Verificare che S7, eseguita su una tabella senza chimere
    senza varianti o con varianti tutte fuori intervallo, si fermi con
    ``E-S7-01`` di revisione umana, con un dettaglio che distingue i due casi,
    senza concludersi e lasciando ``lunghezze.tsv`` per la diagnosi.

    **Razionale scientifico e sistemistico**: Nessuna configurazione rende
    vuota la tabella della versione ridotta, le cui varianti hanno tutte la
    lunghezza del troncamento; la tabella di S6 si altera quindi dopo aver
    preparato il contesto di S7, e la fase vera la legge.
    """
    run = copia_esecuzione(catena_calcolata, tmp_path)
    contesto = run.contesto(Passo.S7)
    (run.albero.cartella(Fase.CHIMERA) / nome_manifesto_passo("S7")).unlink()
    tabella = run.albero.cartella(Fase.CHIMERA) / "tabella_senza_chimere.rds"
    subprocess.run(
        [str(trova_rscript()), "--vanilla", "-e",
         f"t <- readRDS('{tabella}'); {_RENDI_VUOTA[caso]}; saveRDS(t, '{tabella}')"],
        check=True, capture_output=True,
    )
    with pytest.raises(ErrorePipeline) as info:
        run.fase(Passo.S7).esegui(contesto)
    errore = info.value
    assert errore.codice == "E-S7-01"
    assert errore.categoria is Categoria.REVISIONE_UMANA
    assert not errore.ammette_retry
    if caso == "tabella_di_s6_vuota":
        assert "non contiene varianti" in errore.dettaglio
        assert errore.contesto["varianti_escluse"] == 0
    else:
        assert "fuori da" in errore.dettaglio
        assert errore.contesto["varianti_escluse"] > 0
    assert run.albero.manifesto_passo(Passo.S7, Fase.CHIMERA) is None
    assert (run.albero.cartella(Fase.CHIMERA) / "lunghezze.tsv").exists()
