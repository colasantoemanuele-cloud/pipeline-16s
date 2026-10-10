"""Catena completa e analisi ecologiche sul sottoinsieme ridotto di OSD-170.

Il sottoinsieme (``tests/fixtures/osd170/ridotto/``: sei campioni, i due
controlli negativi, poche migliaia di letture per file e un riferimento
sintetico) viene da un dataset che non ha guidato lo sviluppo: letture piu'
lunghe, controlli negativi senza controlli positivi, nessun file del lotto. Il
modulo esegue su di esso la catena dalla riga di comando, con la configurazione
della fixture, e ne verifica invarianti strutturali: nessun numero atteso, che
dipenderebbe dalle versioni dei pacchetti di calcolo.

Richiede R con i pacchetti di calcolo:

    docker run --rm -e PYTHONPATH=/app/src -e AMPLICON16S_R_DIR=/app/R \\
      -v "$(pwd)":/app -w /app <immagine> pytest tests/test_catena_osd170_ridotto.py
"""

from __future__ import annotations

import csv
import functools
import hashlib
import json
import os
import subprocess
from pathlib import Path

import pytest
import yaml
from eco_aiuti import impronte
from sottoinsieme import motivo_pacchetti_r_assenti

from amplicon16s import cli
from amplicon16s.logging.logger import chiudi
from amplicon16s.rbridge.runner import trova_rscript
from amplicon16s_eco.cli import main as eco_main

RADICE = Path(__file__).resolve().parents[1]
FIXTURE = Path(__file__).resolve().parent / "fixtures" / "osd170" / "ridotto"
FINALI = ("conteggi.tsv", "tassonomia.tsv", "metadati.tsv", "sequenze.fasta", "ps_final.rds",
          "ps_controlli.rds", "checksum.sha256")

# Tutti i test del modulo sullo stesso processo: le due esecuzioni della catena
# (fixture di modulo) si preparano una volta sola, invece che una per ogni
# processo a cui pytest-xdist distribuirebbe i test.
pytestmark = pytest.mark.xdist_group("catena_osd170_ridotto")


@functools.cache
def _sonda() -> str | None:
    """Perche' l'ambiente R richiesto non c'e', o ``None``."""
    return motivo_pacchetti_r_assenti(
        "dada2", "ShortRead", "jsonlite", "phyloseq", "Biostrings", "decontam", "vegan", "ape")


@pytest.fixture(scope="module")
def esecuzioni(tmp_path_factory) -> dict[str, Path]:
    """Due esecuzioni della catena sul sottoinsieme, con la configurazione della
    fixture: cambiano solo i percorsi, resi assoluti, e la cartella di uscita.
    """
    if _sonda() is not None:
        if os.environ.get("AMPLICON16S_RICHIEDI_BIOC") == "1":
            pytest.fail(f"Bioconductor e' richiesto in questo ambiente: {_sonda()}")
        pytest.skip(_sonda())
    radice = tmp_path_factory.mktemp("osd170_ridotto")
    dati = yaml.safe_load((FIXTURE / "config.yaml").read_text(encoding="utf-8"))
    for chiave in ("fastq_dir", "assay_table", "study_table"):
        dati["io"][chiave] = str(RADICE / dati["io"][chiave])
    dati["tax"]["ref_fasta"] = str(RADICE / dati["tax"]["ref_fasta"])
    uscite = {}
    for nome in ("prima", "seconda"):
        dati["io"]["out_root"] = str(radice / nome)
        config = radice / f"{nome}.yaml"
        config.write_text(yaml.safe_dump(dati, sort_keys=False), encoding="utf-8")
        chiudi()
        assert cli.main(["validate", "--config", str(config)]) == 0
        chiudi()
        assert cli.main(["resume", "--config", str(config)]) == 0
        chiudi()
        assert cli.main(["report", "--config", str(config)]) == 0
        chiudi()
        uscite[nome] = radice / nome
    return uscite


def _manifesti(uscita: Path) -> list[Path]:
    return sorted(uscita.rglob("manifest_S*.json"))


def _colonna(percorso: Path, nome: str) -> list[str]:
    with open(percorso, encoding="utf-8", newline="") as file:
        return [riga[nome] for riga in csv.DictReader(file, delimiter="\t")]


def _campioni_di(oggetto: Path) -> list[str]:
    """I nomi dei campioni di un oggetto phyloseq salvato."""
    processo = subprocess.run(
        [str(trova_rscript()), "--vanilla", "-e",
         f"cat(phyloseq::sample_names(readRDS('{oggetto}')), sep = '\\n')"],
        capture_output=True, text=True, check=False)
    assert processo.returncode == 0, processo.stderr[-2000:]
    return processo.stdout.split()


def test_la_fixture_basta_a_se_stessa_e_resta_piccola():
    """
    **Obiettivo**: Verificare che la configurazione della fixture indichi solo
    file della fixture stessa, tutti presenti, che l'MD5 dichiarato sia quello
    del riferimento sintetico, che ogni file del manifesto abbia le letture
    dichiarate nei metadati e che l'insieme resti sotto i 5 MB.

    **Razionale scientifico e sistemistico**: Il test deve poter girare dove
    non c'e' nient'altro che il repository: un percorso verso un dato da
    scaricare lo farebbe saltare proprio nell'ambiente che deve presidiare.
    """
    dati = yaml.safe_load((FIXTURE / "config.yaml").read_text(encoding="utf-8"))
    relativa = FIXTURE.relative_to(RADICE).as_posix()
    for chiave in ("fastq_dir", "assay_table", "study_table"):
        assert dati["io"][chiave].startswith(relativa) and (RADICE / dati["io"][chiave]).exists()
    riferimento = RADICE / dati["tax"]["ref_fasta"]
    assert dati["tax"]["ref_fasta"].startswith(relativa)
    assert hashlib.md5(riferimento.read_bytes()).hexdigest() == dati["tax"]["ref_md5"]
    file = sorted(p.name for p in (FIXTURE / "fastq").iterdir())
    assert file == sorted(_colonna(FIXTURE / "manifesto.tsv", "file"))
    campioni = _colonna(FIXTURE / "metadati" / "assay.txt", "Sample Name")
    assert sorted(f"{c}.fastq.gz" for c in campioni) == file
    assert sum(p.stat().st_size for p in FIXTURE.rglob("*") if p.is_file()) < 5 * 1024 * 1024


def test_ogni_fase_ha_il_manifesto_e_ogni_artefatto_corrisponde(esecuzioni):
    """
    **Obiettivo**: Verificare che dopo ``validate``, ``resume`` e ``report``
    (usciti con 0) ogni fase eseguita abbia il proprio manifesto concluso, che
    lo sha256 e la dimensione di ogni artefatto dichiarato corrispondano al
    file, e che i file consegnati siano tutti in ``12_final``.

    **Razionale scientifico e sistemistico**: Il manifesto e' cio' con cui una
    ripresa decide se una fase e' valida e con cui un risultato si riconosce:
    su un dataset nuovo deve valere come su quello di riferimento.
    """
    uscita = esecuzioni["prima"]
    manifesti = _manifesti(uscita)
    passi = {json.loads(m.read_text(encoding="utf-8"))["passo"] for m in manifesti}
    assert passi == {f"S{i}" for i in range(15)} - {"S9"}
    artefatti = 0
    for manifesto in manifesti:
        contenuto = json.loads(manifesto.read_text(encoding="utf-8"))
        assert contenuto["conclusa"], manifesto
        for artefatto in contenuto["artefatti"]:
            file = manifesto.parent / artefatto["nome"]
            assert file.is_file(), file
            assert file.stat().st_size == artefatto["byte"], file
            assert "sha256:" + hashlib.sha256(file.read_bytes()).hexdigest() == artefatto["checksum"], file
            artefatti += 1
    assert artefatti > 50
    for nome in FINALI:
        assert (uscita / "12_final" / nome).is_file(), nome
    assert (uscita / "report" / "report.html").is_file()


def test_campioni_e_varianti_sono_gli_stessi_in_tutti_i_file_finali(esecuzioni):
    """
    **Obiettivo**: Verificare che conteggi e metadati abbiano gli stessi
    campioni, tutti biologici; che conteggi, tassonomia e sequenze abbiano le
    stesse varianti nello stesso ordine; che i due controlli negativi siano in
    ``ps_controlli.rds`` e non in ``ps_final.rds``, e i campioni dell'oggetto
    finale siano quelli delle tabelle.

    **Razionale scientifico e sistemistico**: I file piatti sono lo stesso
    risultato in un'altra forma: uno scarto fra loro, o un controllo rimasto
    fra i campioni, darebbe ad analisi diverse dati diversi.
    """
    finale = esecuzioni["prima"] / "12_final"
    with open(finale / "conteggi.tsv", encoding="utf-8") as file:
        campioni = file.readline().rstrip("\n").split("\t")[1:]
    assert campioni and campioni == _colonna(finale / "metadati.tsv", "accession")
    assert set(_colonna(finale / "metadati.tsv", "classe")) == {"biologico"}
    varianti = _colonna(finale / "conteggi.tsv", "asv_id")
    assert varianti and varianti == _colonna(finale / "tassonomia.tsv", "asv_id")
    intestazioni = [r[1:].split()[0] for r in (finale / "sequenze.fasta").read_text(
        encoding="ascii").splitlines() if r.startswith(">")]
    assert intestazioni == varianti

    negativi = {c for c, t in zip(_colonna(FIXTURE / "metadati" / "studio.txt", "Sample Name"),
                                  _colonna(FIXTURE / "metadati" / "studio.txt", "Factor Value[Time]"))
                if t == "Negative control"}
    assert len(negativi) == 2
    assert sorted(_campioni_di(finale / "ps_final.rds")) == sorted(campioni)
    assert set(_campioni_di(finale / "ps_controlli.rds")) == negativi
    assert not negativi & set(campioni)


def test_due_esecuzioni_danno_gli_stessi_byte(esecuzioni):
    """
    **Obiettivo**: Verificare che due esecuzioni della stessa configurazione
    in cartelle diverse diano, per ogni artefatto di ogni fase, lo stesso
    sha256, salvo la configurazione registrata (che riporta la cartella di
    uscita), e che i file consegnati siano identici byte per byte.

    **Razionale scientifico e sistemistico**: La riproducibilita' e' un
    requisito dello strumento, non una proprieta' del dataset su cui e' nato:
    va dimostrata anche su dati con un'altra struttura.
    """
    sha = {}
    for nome, uscita in esecuzioni.items():
        sha[nome] = {
            f"{m.parent.name}/{a['nome']}": a["checksum"]
            for m in _manifesti(uscita)
            for a in json.loads(m.read_text(encoding="utf-8"))["artefatti"]}
    assert sha["prima"] == sha["seconda"]
    for nome in FINALI:
        assert (esecuzioni["prima"] / "12_final" / nome).read_bytes() == (
            esecuzioni["seconda"] / "12_final" / nome).read_bytes(), nome


def test_le_analisi_ecologiche_si_eseguono_sul_risultato(esecuzioni, tmp_path):
    """
    **Obiettivo**: Verificare che ``eco run`` sull'oggetto finale, con la
    configurazione della fixture, esca con 0; che il manifesto dell'analisi
    elenchi tutti e soli i file scritti, con sha256 e dimensione di ciascuno;
    che riporti l'impronta dell'oggetto di partenza; e che una seconda
    esecuzione dia gli stessi file.

    **Razionale scientifico e sistemistico**: Le analisi ecologiche partono
    dall'oggetto consegnato dalla pipeline: la catena e' completa solo se
    quell'oggetto, prodotto da un dataset qualunque, e' davvero utilizzabile.
    """
    oggetto = esecuzioni["prima"] / "12_final" / "ps_final.rds"
    config = FIXTURE / "eco.yaml"
    for nome in ("a", "b"):
        assert eco_main(["run", "--object", str(oggetto), "--config", str(config),
                         "--out", str(tmp_path / nome)]) == 0
    uscita = tmp_path / "a"
    manifesto = json.loads((uscita / "manifest.json").read_text(encoding="utf-8"))
    scritti = impronte(uscita)
    assert set(scritti) - {"manifest.json"} == {f["nome"] for f in manifesto["file"]} | {
        manifesto["configurazione"]["file"]} - {"manifest.json"}
    for voce in manifesto["file"]:
        assert "sha256:" + scritti[voce["nome"]] == voce["sha256"], voce["nome"]
        assert (uscita / voce["nome"]).stat().st_size == voce["byte"]
    assert manifesto["oggetto_di_partenza"]["sha256"] == (
        "sha256:" + hashlib.sha256(oggetto.read_bytes()).hexdigest())
    assert impronte(tmp_path / "a") == impronte(tmp_path / "b")
