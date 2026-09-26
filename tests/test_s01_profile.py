"""Test del sottoinsieme di prova e della fase S1.

La prima parte verifica la versione ridotta del sottoinsieme, che vive nel
repository: integrità e proprietà che ne giustificano la selezione. Non
richiede R.

La seconda esegue S1 davvero, sulla versione ridotta: R, ShortRead e
Biostrings attraverso il ponte. Senza Bioconductor quei test si saltano; nella
catena di integrazione continua girano dentro l'immagine della pipeline, dove
``AMPLICON16S_RICHIEDI_BIOC`` trasforma il salto in un fallimento.

L'ultima parte, sul dataset completo, richiede i dati reali e Bioconductor:
gira nel container, che è il riferimento.
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import itertools
import json
import logging
import os
import subprocess
import sys
from collections import Counter
from pathlib import Path

import pytest
from sottoinsieme import RIDOTTO, config_ridotta, manifesto, selezione

from amplicon16s.config.schema import carica, valida
from amplicon16s.io_layer.artifacts import Fase
from amplicon16s.logging.logger import NOME_FILE_LOG, chiudi, configura
from amplicon16s.rbridge.runner import trova_rscript
from amplicon16s.runner.executor import Conclusione, Esecutore
from amplicon16s.runner.graph import Passo
from amplicon16s.runner.project import ProjectRun

RADICE = Path(__file__).resolve().parents[1]


@pytest.fixture(autouse=True)
def uscite_pulite():
    chiudi()
    yield
    chiudi()


# --------------------------------------------------------------------------- #
# La selezione e la versione ridotta                                           #
# --------------------------------------------------------------------------- #


def test_la_selezione_e_motivata_campione_per_campione():
    righe = selezione()
    assert len(righe) == 28
    assert len({r["accession"] for r in righe}) == 28
    for riga in righe:
        assert riga["motivo"].strip(), riga["accession"]
        assert len(riga["md5"]) == 32
        assert riga["run_ena"].startswith("ERR")


def test_la_selezione_esercita_decontaminazione_e_calibrazione():
    righe = selezione()
    # Le tre classi e le due corse.
    assert Counter(r["classe"] for r in righe) == {
        "biologico": 9, "controllo_positivo": 9, "controllo_negativo": 10,
    }
    assert len({r["corsa"] for r in righe}) == 2
    # Due piastre con i cinque negativi che decontam.min_blanks richiede, e
    # biologici della stessa piastra con cui confrontarli.
    for piastra in ("4", "10"):
        della = [r for r in righe if r["piastra"] == piastra]
        assert sum(r["classe"] == "controllo_negativo" for r in della) == 5
        assert any(r["classe"] == "biologico" for r in della)
    # Una serie di calibrazione completa, otto livelli nella stessa piastra.
    serie = [r for r in righe if r["classe"] == "controllo_positivo" and r["piastra"] == "10"]
    assert len(serie) == 8


def test_la_selezione_contiene_i_casi_noti():
    per_accession = {r["accession"]: r for r in selezione()}
    assert per_accession["ERX12083125"]["lunghezza_minima"] == "137"  # BLANK.3DMM.P4.8C
    assert per_accession["ERX12083170"]["campione"] == "JLP1A1.L4"
    letture = sorted(int(r["letture"]) for r in per_accession.values())
    assert letture[0] == 44 and letture[-1] == 228431  # il piu' povero e il piu' profondo


def test_la_versione_ridotta_e_integra_e_piccola():
    voci = manifesto()
    cartella = RIDOTTO / "fastq"
    assert {v["file"] for v in voci} == {p.name for p in cartella.iterdir()}
    for voce in voci:
        percorso = cartella / voce["file"]
        assert hashlib.md5(percorso.read_bytes()).hexdigest() == voce["md5"], voce["file"]
    totale = sum(p.stat().st_size for p in cartella.iterdir())
    assert totale < 3 * 1024 * 1024


def _lunghezze(percorso: Path) -> list[int]:
    with gzip.open(percorso, "rt") as file:
        return [len(r.rstrip("\n")) for i, r in enumerate(file) if i % 4 == 1]


def test_la_versione_ridotta_conserva_le_proprieta_della_selezione():
    per_accession = {r["accession"]: r for r in selezione()}
    for voce in manifesto():
        accession = voce["file"].split("_")[0]
        lunghezze = _lunghezze(RIDOTTO / "fastq" / voce["file"])
        assert len(lunghezze) == int(voce["letture"])
        # Ogni file conserva la propria lunghezza minima, 137 compreso.
        assert min(lunghezze) == int(per_accession[accession]["lunghezza_minima"])
    letture = {v["file"].split("_")[0]: int(v["letture"]) for v in manifesto()}
    # Il contrasto fra il campione piu' profondo e il piu' povero resta.
    assert letture["ERX12084006"] > 50 * letture["ERX12083285"]
    assert letture["ERX12083285"] == 44  # sotto il minimo: tutte le letture


def test_i_metadati_ridotti_coprono_esattamente_la_selezione():
    accession = {r["accession"] for r in selezione()}
    with open(RIDOTTO / "metadati" / "lotti.tsv", encoding="utf-8") as file:
        assert {r["experiment_accession"] for r in csv.DictReader(file, delimiter="\t")} == accession


def test_i_fastq_di_prova_non_sono_ignorati_da_git():
    dentro = subprocess.run(
        ["git", "rev-parse", "--is-inside-work-tree"],
        cwd=RADICE, capture_output=True, check=False,
    )
    if dentro.returncode != 0:
        pytest.skip("non e' un repository git: nell'immagine .git non viene copiato")
    esempio = next((RIDOTTO / "fastq").iterdir())
    esito = subprocess.run(
        ["git", "check-ignore", "-q", str(esempio)], cwd=RADICE, check=False
    )
    assert esito.returncode == 1  # 1: il file non e' ignorato


# --------------------------------------------------------------------------- #
# S1 sulla versione ridotta                                                    #
# --------------------------------------------------------------------------- #


def _motivo_bioconductor_assente() -> str | None:
    rscript = trova_rscript()
    if rscript is None:
        return "Rscript non disponibile"
    prova = subprocess.run(
        [str(rscript), "--vanilla", "-e",
         "quit(status = !all(vapply(c('ShortRead', 'Biostrings', 'jsonlite'), "
         "requireNamespace, logical(1), quietly = TRUE)))"],
        capture_output=True, check=False,
    )
    if prova.returncode != 0:
        return "ShortRead, Biostrings o jsonlite non installati"
    return None


_MOTIVO_BIOC_ASSENTE = _motivo_bioconductor_assente()


@pytest.fixture
def bioconductor():
    """Richiede R con Bioconductor: salta senza, ma in CI fallisce."""
    if _MOTIVO_BIOC_ASSENTE is not None:
        if os.environ.get("AMPLICON16S_RICHIEDI_BIOC") == "1":
            pytest.fail(f"Bioconductor e' richiesto in questo ambiente: {_MOTIVO_BIOC_ASSENTE}")
        pytest.skip(_MOTIVO_BIOC_ASSENTE)


def _fino_a_s1(config):
    run = ProjectRun(config)
    return run, Esecutore(run, fino_a=Passo.S1).esegui()


def _profili(run) -> Path:
    return run.albero.cartella(Fase.QC_PROFILES)


@pytest.fixture(scope="module")
def eseguita(tmp_path_factory):
    """S0 e S1 sulla versione ridotta, una volta sola per i test che leggono.

    Caricare ShortRead costa da solo alcuni secondi: ripeterlo per ogni test
    non aggiungerebbe nulla alla verifica.
    """
    if _MOTIVO_BIOC_ASSENTE is not None:
        return None
    return _fino_a_s1(config_ridotta(tmp_path_factory.mktemp("s1")))


def test_s1_scrive_i_profili_con_il_proprio_manifesto(bioconductor, eseguita):
    run, esito = eseguita
    assert esito.conclusione is Conclusione.COMPLETATA
    assert [r.passo for r in esito.eseguite] == [Passo.S0, Passo.S1]

    manifesto_s1 = run.albero.manifesto_passo(Passo.S1, Fase.QC_PROFILES)
    assert set(manifesto_s1.nomi) == {
        "lunghezze.tsv", "qualita.tsv", "letture_grezze.tsv", "riepilogo.json"
    }
    assert run.valuta().situazioni[Passo.S1].stato.value == "completata"


def test_s1_legge_tutte_le_letture(bioconductor, eseguita):
    run, _ = eseguita
    attese = {v["file"].split("_")[0]: int(v["letture"]) for v in manifesto()}
    with open(_profili(run) / "letture_grezze.tsv", encoding="utf-8") as file:
        contate = {r["campione"]: int(r["letture"]) for r in csv.DictReader(file, delimiter="\t")}
    assert contate == attese

    riepilogo = json.loads((_profili(run) / "riepilogo.json").read_text())
    assert riepilogo["campioni"] == 28
    assert riepilogo["letture"] == sum(attese.values())
    assert riepilogo["lunghezza_minima"] == 137
    assert riepilogo["moda"] == 151


def test_s1_scrive_distribuzione_e_profilo_per_campione(bioconductor, eseguita):
    run, _ = eseguita
    with open(_profili(run) / "lunghezze.tsv", encoding="utf-8") as file:
        lunghezze = list(csv.DictReader(file, delimiter="\t"))
    blank = {int(r["lunghezza"]): int(r["letture"]) for r in lunghezze if r["campione"] == "ERX12083125"}
    attese = Counter(_lunghezze(next((RIDOTTO / "fastq").glob("ERX12083125_*"))))
    assert blank == dict(attese)

    with open(_profili(run) / "qualita.tsv", encoding="utf-8") as file:
        qualita = [r for r in csv.DictReader(file, delimiter="\t") if r["campione"] == "ERX12083125"]
    assert [int(r["posizione"]) for r in qualita] == list(range(1, 152))
    # Alla posizione 137 ci sono tutte le letture, alla 151 solo le complete.
    assert int(qualita[136]["letture"]) == sum(attese.values())
    assert int(qualita[150]["letture"]) == attese[151]
    for riga in qualita:
        assert int(riga["p10"]) <= int(riga["mediana"]) <= int(riga["p90"])


def test_con_truncLen_137_non_c_e_scarto(bioconductor, eseguita):
    run, _ = eseguita
    assert run.albero.manifesto_passo(Passo.S1, Fase.QC_PROFILES).degradazioni == ()


def test_con_truncLen_120_lo_scarto_e_registrato_da_s1_non_da_s0(bioconductor, tmp_path):
    configura(tmp_path / "out")
    run, esito = _fino_a_s1(config_ridotta(tmp_path, filter={"truncLen": 120}))
    assert esito.conclusione is Conclusione.COMPLETATA

    (degradazione,) = run.albero.manifesto_passo(Passo.S1, Fase.QC_PROFILES).degradazioni
    assert degradazione["codice"] == "E-S1-01"
    assert degradazione["scarto"] == 17
    s0 = run.albero.manifesto_passo(Passo.S0, Fase.INPUT_VALIDATION).degradazioni
    assert "E-S1-01" not in {d["codice"] for d in s0}

    righe = [json.loads(r) for r in (tmp_path / "out" / "99_logs" / NOME_FILE_LOG).read_text().splitlines()]
    eventi = [r for r in righe if r.get("codice") == "E-S1-01"]
    assert len(eventi) == 1 and eventi[0]["fase"] == "S1"


def test_s1_ferma_se_il_troncamento_supera_il_minimo_vero(bioconductor, tmp_path):
    """G09 vede solo le prime letture e passa; S1 le legge tutte e si ferma."""
    def prime(percorso, k):
        with gzip.open(percorso, "rt") as file:
            return [len(r.rstrip()) for i, r in enumerate(itertools.islice(file, 4 * k)) if i % 4 == 1]

    stima = min(min(prime(p, 20)) for p in (RIDOTTO / "fastq").iterdir())
    assert stima > 137
    config = config_ridotta(tmp_path, qc={"head_reads": 20}, filter={"truncLen": stima})
    run, esito = _fino_a_s1(config)

    assert esito.conclusione is Conclusione.ARRESTATA
    assert [r.passo for r in esito.eseguite] == [Passo.S0]  # G09 e' passato
    assert esito.punto.passo is Passo.S1
    assert esito.punto.codice == "E-S1-02"
    assert "137 bp" in esito.punto.dettaglio
    # I profili restano per decidere il troncamento, ma S1 non e' conclusa.
    assert (_profili(run) / "riepilogo.json").exists()
    assert run.albero.manifesto_passo(Passo.S1, Fase.QC_PROFILES) is None


# --------------------------------------------------------------------------- #
# Dataset completo, nel container                                              #
# --------------------------------------------------------------------------- #


@pytest.mark.dati_reali
def test_s1_sul_dataset_completo(bioconductor, tmp_path):
    percorso = os.environ.get("AMPLICON16S_CONFIG_DATI_REALI")
    if not percorso or not Path(percorso).expanduser().is_file():
        pytest.skip("AMPLICON16S_CONFIG_DATI_REALI non impostata o file assente")
    dati = carica(Path(percorso).expanduser()).model_dump(mode="python")
    dati["io"]["out_root"] = str(tmp_path / "out")
    logging.getLogger("amplicon16s").setLevel(logging.WARNING)

    run, esito = _fino_a_s1(valida(dati))
    assert esito.conclusione is Conclusione.COMPLETATA
    riepilogo = json.loads((_profili(run) / "riepilogo.json").read_text())
    assert riepilogo["campioni"] == 960
    assert riepilogo["lunghezza_minima"] == 137
    assert riepilogo["moda"] == 151
    assert riepilogo["qualita_mediana_minima"] >= 25


@pytest.mark.dati_reali
def test_lo_script_ricostruisce_la_versione_ridotta_identica(tmp_path, monkeypatch, capsys):
    """Dai dati locali, la versione ridotta rigenerata coincide con quella versionata."""
    percorso = os.environ.get("AMPLICON16S_CONFIG_DATI_REALI")
    if not percorso or not Path(percorso).expanduser().is_file():
        pytest.skip("AMPLICON16S_CONFIG_DATI_REALI non impostata o file assente")
    import importlib.util

    spec = importlib.util.spec_from_file_location(
        "build_test_subset", RADICE / "scripts" / "build_test_subset.py"
    )
    modulo = importlib.util.module_from_spec(spec)
    # Le dataclass del modulo cercano il modulo in sys.modules.
    monkeypatch.setitem(sys.modules, "build_test_subset", modulo)
    spec.loader.exec_module(modulo)

    modulo.main(["ridotto", "--config", str(Path(percorso).expanduser()),
                 "--destinazione", str(tmp_path)])
    rigenerato = (tmp_path / "manifesto.tsv").read_text(encoding="utf-8")
    assert rigenerato == (RIDOTTO / "manifesto.tsv").read_text(encoding="utf-8")
    for nome in ("assay.txt", "studio.txt", "lotti.tsv"):
        assert (tmp_path / "metadati" / nome).read_bytes() == (RIDOTTO / "metadati" / nome).read_bytes()
