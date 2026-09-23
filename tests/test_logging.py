"""Test della registrazione degli eventi e dell'esecuzione fittizia completa."""

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
    testo = percorso.read_text(encoding="utf-8")
    return [json.loads(r) for r in testo.splitlines() if r.strip()]


# --------------------------------------------------------------------------- #
# Le due uscite                                                                #
# --------------------------------------------------------------------------- #


def test_il_log_strutturato_nasce_sotto_99_logs(tmp_path):
    percorso = configura(tmp_path)
    assert percorso == tmp_path / Fase.LOGS.value / NOME_FILE_LOG
    assert percorso.parent.is_dir()


def test_senza_radice_resta_la_sola_console(tmp_path):
    assert configura(None) is None
    tipi = [type(h) for h in ottieni().handlers]
    assert tipi == [logging.StreamHandler]
    assert not list(tmp_path.iterdir())


def test_le_uscite_sono_due(tmp_path):
    configura(tmp_path)
    assert len(ottieni().handlers) == 2


def test_configurare_piu_volte_non_accumula_uscite(tmp_path):
    configura(tmp_path)
    configura(tmp_path)
    configura(tmp_path)
    assert len(ottieni().handlers) == 2


def test_il_file_ruota(tmp_path):
    configura(tmp_path)
    rotanti = [
        h for h in ottieni().handlers
        if isinstance(h, logging.handlers.RotatingFileHandler)
    ]
    assert len(rotanti) == 1
    assert rotanti[0].maxBytes == MAX_BYTE
    assert rotanti[0].backupCount > 0


def test_la_rotazione_avviene_davvero(tmp_path):
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
    percorso = configura(tmp_path)
    ottieni("s2").info("filtro concluso", extra={"fase": "S2", "conservate": 0.83})
    chiudi()

    evento = _righe(percorso)[0]
    assert evento["fase"] == "S2"
    assert evento["conservate"] == 0.83
    assert evento["origine"] == "amplicon16s.s2"
    assert "istante" in evento


def test_gli_eventi_si_filtrano_per_codice(tmp_path):
    """È il motivo del formato: si interroga, non si legge."""
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
    """Perdere il log proprio quando serve sarebbe il guasto peggiore."""
    percorso = configura(tmp_path)
    ottieni("s0").info("percorso", extra={"dove": tmp_path, "fase_enum": Fase.CONFIG})
    chiudi()

    evento = _righe(percorso)[0]
    assert evento["dove"] == str(tmp_path)
    assert evento["fase_enum"] == Fase.CONFIG.value


def test_l_eccezione_finisce_nel_log(tmp_path):
    percorso = configura(tmp_path)
    try:
        raise errore("E-S14-01")
    except Exception:
        ottieni("s14").exception("validazione finale fallita")
    chiudi()

    assert "E-S14-01" in _righe(percorso)[0]["eccezione"]


def test_il_livello_del_file_e_piu_verboso_della_console(tmp_path):
    percorso = configura(tmp_path, livello_console=logging.WARNING)
    ottieni("x").debug("dettaglio diagnostico")
    chiudi()
    assert [r["messaggio"] for r in _righe(percorso)] == ["dettaglio diagnostico"]


# --------------------------------------------------------------------------- #
# Esecuzione fittizia: albero completo, configurazione risolta, log            #
# --------------------------------------------------------------------------- #


def test_esecuzione_fittizia_produce_albero_configurazione_e_log(tmp_path):
    """Mette insieme i tre servizi come farà una fase vera."""
    radice = tmp_path / "esecuzione"

    albero = AlberoOutput(radice)
    albero.crea()
    percorso_log = configura(radice)
    log = ottieni("runner")

    dati = yaml.safe_load(ESEMPIO.read_text(encoding="utf-8"))
    risolta = esegui_g15(dati).con_campioni_biologici(803)
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
