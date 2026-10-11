"""Catena completa in modalita' appaiata sul sottoinsieme ridotto di OSD-276.

Il sottoinsieme (``tests/fixtures/osd276/appaiato/``: quattro campioni, le
prime coppie di ogni corsa, un riferimento tassonomico sintetico) e' un dataset
con le due letture di ogni frammento. Il modulo esegue su di esso la catena
dalla riga di comando, con la configurazione della fixture, e ne verifica
invarianti strutturali: nessun numero atteso, che dipenderebbe dalle versioni
dei pacchetti di calcolo.

Verifica inoltre cio' che la modalita' appaiata aggiunge al calcolo: che le
varianti di ciascuna lettura siano quelle di ``dada(pool = "pseudo")``, che la
fusione dalle parti conservate degli oggetti di dada2 sia quella di
``mergePairs`` sugli oggetti interi, e che la fusione si possa rifare senza
cambiare un byte, anche dopo che le letture filtrate sono state tolte.

Richiede R con i pacchetti di calcolo:

    docker run --rm -e PYTHONPATH=/app/src -e AMPLICON16S_R_DIR=/app/R \\
      -v "$(pwd)":/app -w /app <immagine> pytest tests/test_catena_appaiata.py
"""

from __future__ import annotations

import csv
import functools
import hashlib
import json
import os
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path

import pytest
import yaml
from sottoinsieme import motivo_pacchetti_r_assenti

from amplicon16s import cli, interfaccia
from amplicon16s.config.schema import carica, valida
from amplicon16s.io_layer import reads
from amplicon16s.logging.logger import chiudi
from amplicon16s.rbridge.runner import trova_rscript
from amplicon16s.runner.project import ProjectRun, StatoPasso

RADICE = Path(__file__).resolve().parents[1]
FIXTURE = Path(__file__).resolve().parent / "fixtures" / "osd276" / "appaiato"
FINALI = ("conteggi.tsv", "tassonomia.tsv", "metadati.tsv", "sequenze.fasta", "ps_final.rds",
          "checksum.sha256")
# Le fasi che non leggono ne' le letture filtrate ne' cio' che ne deriva.
PRIMA_DEL_FILTRO = ("S0", "S1")
FINO_ALL_INFERENZA = (*PRIMA_DEL_FILTRO, "S2", "S3", "S4")

# Tutti i test del modulo sullo stesso processo: le esecuzioni della catena
# (fixture di modulo) si preparano una volta sola, invece che una per ogni
# processo a cui pytest-xdist distribuirebbe i test.
pytestmark = pytest.mark.xdist_group("catena_appaiata")


@functools.cache
def _sonda() -> str | None:
    """Perche' l'ambiente R richiesto non c'e', o ``None``."""
    return motivo_pacchetti_r_assenti(
        "dada2", "ShortRead", "jsonlite", "phyloseq", "Biostrings", "decontam")


def _dati() -> dict:
    """La configurazione della fixture, con i percorsi resi assoluti."""
    dati = yaml.safe_load((FIXTURE / "config.yaml").read_text(encoding="utf-8"))
    for chiave in ("fastq_dir", "assay_table", "study_table"):
        dati["io"][chiave] = str(RADICE / dati["io"][chiave])
    dati["tax"]["ref_fasta"] = str(RADICE / dati["tax"]["ref_fasta"])
    return dati


def _manifesti(uscita: Path) -> dict[str, Path]:
    """I manifesti di fase di un'esecuzione, per passo."""
    return {m.stem.removeprefix("manifest_"): m for m in uscita.rglob("manifest_S*.json")}


def _impronte(uscita: Path) -> dict[str, str]:
    """Lo sha256 dichiarato di ogni artefatto di ogni fase."""
    return {
        f"{m.parent.name}/{a['nome']}": a["checksum"]
        for m in _manifesti(uscita).values()
        for a in json.loads(m.read_text(encoding="utf-8"))["artefatti"]}


def _colonna(percorso: Path, nome: str) -> list[str]:
    with open(percorso, encoding="utf-8", newline="") as file:
        return [riga[nome] for riga in csv.DictReader(file, delimiter="\t")]


def _righe(percorso: Path) -> list[dict[str, str]]:
    with open(percorso, encoding="utf-8", newline="") as file:
        return list(csv.DictReader(file, delimiter="\t"))


def _esegui(dati: dict, radice: Path, nome: str, *, report: bool = True) -> Path:
    """Scrive la configurazione ed esegue la catena nella cartella indicata."""
    dati["io"]["out_root"] = str(radice / nome)
    config = radice / f"{nome}.yaml"
    config.write_text(yaml.safe_dump(dati, sort_keys=False), encoding="utf-8")
    comandi = ("validate", "resume", "report") if report else ("validate", "resume")
    for comando in comandi:
        chiudi()
        assert cli.main([comando, "--config", str(config)]) == 0, (nome, comando)
    chiudi()
    return radice / nome


@dataclass(frozen=True)
class Ripresa:
    """Una ripresa dopo aver tolto il manifesto della fusione: lo stato prima,
    le impronte e i manifesti prima e dopo.
    """

    uscita: Path
    da_eseguire: tuple[str, ...]
    impronte_prima: dict[str, str]
    impronte_dopo: dict[str, str]
    manifesti_prima: dict[str, bytes]
    manifesti_dopo: dict[str, bytes]

    @property
    def rifatte(self) -> set[str]:
        """Le fasi il cui manifesto e' stato riscritto."""
        return {p for p, testo in self.manifesti_dopo.items() if testo != self.manifesti_prima[p]}


def _rifai_la_fusione(uscita: Path, config: Path) -> Ripresa:
    """Toglie il manifesto della fusione e riprende l'esecuzione."""
    impronte = _impronte(uscita)
    manifesti = {p: m.read_bytes() for p, m in _manifesti(uscita).items()}
    _manifesti(uscita)["SM"].unlink()
    chiudi()
    situazioni = ProjectRun(carica(config)).valuta().situazioni
    da_eseguire = tuple(
        str(p) for p, s in situazioni.items() if s.stato is StatoPasso.DA_ESEGUIRE)
    assert "SM" in [
        f["passo"] for f in interfaccia.stato(uscita)["fasi"] if f["stato"] == "da_eseguire"]
    assert cli.main(["resume", "--config", str(config)]) == 0
    chiudi()
    return Ripresa(
        uscita, da_eseguire, impronte, _impronte(uscita), manifesti,
        {p: m.read_bytes() for p, m in _manifesti(uscita).items()})


@dataclass(frozen=True)
class Esecuzioni:
    """Le esecuzioni della catena sulla fixture, con le impronte di ciascuna
    come erano alla fine, prima di ogni ripresa.
    """

    prima: Path
    seconda: Path
    senza_filtrate: Path
    impronte: dict[str, dict[str, str]]
    ripresa: Ripresa
    ripresa_senza_filtrate: Ripresa


@pytest.fixture(scope="module")
def esecuzioni(tmp_path_factory) -> Esecuzioni:
    """Tre esecuzioni della catena con la configurazione della fixture: due
    uguali, e una che toglie le letture filtrate a catena conclusa; poi, sulla
    seconda e sulla terza, una ripresa che rifa' la fusione.
    """
    if _sonda() is not None:
        if os.environ.get("AMPLICON16S_RICHIEDI_BIOC") == "1":
            pytest.fail(f"Bioconductor e' richiesto in questo ambiente: {_sonda()}")
        pytest.skip(_sonda())
    radice = tmp_path_factory.mktemp("osd276_appaiato")
    prima = _esegui(_dati(), radice, "prima")
    seconda = _esegui(_dati(), radice, "seconda", report=False)
    dati = _dati()
    dati["run"]["keep_filtered_fastq"] = False
    senza = _esegui(dati, radice, "senza_filtrate", report=False)
    impronte = {u.name: _impronte(u) for u in (prima, seconda, senza)}
    return Esecuzioni(
        prima, seconda, senza, impronte,
        _rifai_la_fusione(seconda, radice / "seconda.yaml"),
        _rifai_la_fusione(senza, radice / "senza_filtrate.yaml"))


def _r(codice: str) -> str:
    """Esegue codice R e ne restituisce l'uscita."""
    processo = subprocess.run(
        [str(trova_rscript()), "--vanilla", "-e", codice],
        capture_output=True, text=True, check=False)
    assert processo.returncode == 0, processo.stderr[-3000:]
    return processo.stdout


def test_la_fixture_basta_a_se_stessa_e_resta_piccola():
    """
    **Obiettivo**: Verificare che la configurazione della fixture sia valida e
    in modalita' appaiata, indichi solo file della fixture stessa, tutti
    presenti; che l'MD5 dichiarato sia quello del riferimento sintetico; che
    ogni campione abbia i suoi due file, con lo stesso numero di letture nello
    stesso ordine, come dichiara il manifesto; e che l'insieme resti sotto i
    3 MB.

    **Razionale scientifico e sistemistico**: Il test deve poter girare dove
    non c'e' nient'altro che il repository, e le coppie della fixture devono
    essere coppie vere: due file tagliati in punti diversi darebbero una
    catena che si ferma, o peggio che fonde letture di frammenti diversi.
    """
    dati = yaml.safe_load((FIXTURE / "config.yaml").read_text(encoding="utf-8"))
    relativa = FIXTURE.relative_to(RADICE).as_posix()
    for chiave in ("fastq_dir", "assay_table", "study_table"):
        assert dati["io"][chiave].startswith(relativa) and (RADICE / dati["io"][chiave]).exists()
    assert dati["tax"]["ref_fasta"].startswith(relativa)
    riferimento = RADICE / dati["tax"]["ref_fasta"]
    assert hashlib.md5(riferimento.read_bytes()).hexdigest() == dati["tax"]["ref_md5"]
    config = valida(_dati())
    assert config.appaiata

    file = sorted(p.name for p in (FIXTURE / "fastq").iterdir())
    assert file == sorted(_colonna(FIXTURE / "manifesto.tsv", "file"))
    abbinamento = reads.abbina_file(file, config.pair.file_markers)
    assert abbinamento.completo and len(abbinamento.coppie) == len(file) // 2 == 4
    dichiarate = {r["file"]: int(r["coppie"]) for r in _righe(FIXTURE / "manifesto.tsv")}
    for forward, reverse in abbinamento.coppie.items():
        sincronia = reads.verifica_sincronia(
            FIXTURE / "fastq" / forward, FIXTURE / "fastq" / reverse)
        assert sincronia.difetto is None, sincronia.difetto
        assert sincronia.coppie == dichiarate[forward] == dichiarate[reverse]
    assert sum(p.stat().st_size for p in FIXTURE.rglob("*") if p.is_file()) < 3 * 1024 * 1024


def test_ogni_fase_ha_il_manifesto_e_ogni_artefatto_corrisponde(esecuzioni):
    """
    **Obiettivo**: Verificare che dopo ``validate``, ``resume`` e ``report``
    (usciti con 0) ogni fase eseguita, fusione compresa, abbia il proprio
    manifesto concluso; che lo sha256 e la dimensione di ogni artefatto
    dichiarato corrispondano al file; che ogni manifesto dichiari la modalita'
    fra i parametri da cui la fase dipende; e che ``status`` dia l'esecuzione
    per conclusa, con la fusione fra l'inferenza e la tabella.

    **Razionale scientifico e sistemistico**: Il manifesto e' cio' con cui una
    ripresa decide se una fase e' valida: in modalita' appaiata ogni fase
    dipende dalla modalita', perche' le stesse letture forward in un dataset
    single-end darebbero un altro risultato.
    """
    uscita = esecuzioni.prima
    manifesti = _manifesti(uscita)
    assert list(manifesti) and "SM" in manifesti
    assert set(manifesti) == {"S0", "S1", "S2", "S3", "S4", "SM", "S5", "S6", "S7", "S8",
                              "S10", "S11", "S12", "S13", "S14"}
    for passo, manifesto in manifesti.items():
        contenuto = json.loads(manifesto.read_text(encoding="utf-8"))
        assert "io.layout" in contenuto["calcolata_su"]["configurazione"]["parametri"], passo
        for artefatto in contenuto["artefatti"]:
            file = manifesto.parent / artefatto["nome"]
            assert file.stat().st_size == artefatto["byte"], file
            assert "sha256:" + hashlib.sha256(file.read_bytes()).hexdigest() == (
                artefatto["checksum"]), file
    for nome in FINALI:
        assert (uscita / "12_final" / nome).is_file(), nome
    report = (uscita / "report" / "report.html").read_text(encoding="utf-8")
    assert "Fusione delle coppie" in report and "quota di coppie fuse" in report
    assert len(_righe(uscita / "report" / "tabelle" / "fusione_per_campione.tsv")) == 4
    stato = interfaccia.stato(uscita)
    passi = [f["passo"] for f in stato["fasi"]]
    assert stato["stato"] == "conclusa"
    assert passi.index("SM") == passi.index("S4") + 1 == passi.index("S5") - 1


def test_le_due_letture_si_filtrano_e_si_modellano_separatamente(esecuzioni):
    """
    **Obiettivo**: Verificare che ogni campione abbia i due file filtrati con
    lo stesso numero di letture; che il modello d'errore delle letture reverse
    sia un file distinto da quello delle forward, con contenuto diverso; e che
    l'inferenza abbia, per ciascuna lettura, le varianti di ogni campione e
    cio' che la fusione legge.

    **Razionale scientifico e sistemistico**: Il filtro appaiato scarta la
    coppia se una delle due letture non passa: i due file devono restare
    allineati. Le letture reverse hanno un profilo d'errore proprio, e un solo
    modello per entrambe attribuirebbe alle une gli errori delle altre.
    """
    uscita = esecuzioni.prima
    filtrate = dict(zip(_colonna(uscita / "03_filtered" / "letture_filtrate.tsv", "campione"),
                        _colonna(uscita / "03_filtered" / "letture_filtrate.tsv", "letture")))
    assert len(filtrate) == 4
    for accession, letture in filtrate.items():
        coppia = (uscita / "03_filtered" / f"{accession}_filt.fastq.gz",
                  uscita / "03_filtered" / f"{accession}_filt_reverse.fastq.gz")
        sincronia = reads.verifica_sincronia(*coppia)
        assert sincronia.difetto is None and sincronia.coppie == int(letture) > 0, accession

    corrispondenza = _righe(uscita / "04_error_models" / "corrispondenza.tsv")
    assert {r["campione"] for r in corrispondenza} == set(filtrate)
    for riga in corrispondenza:
        assert riga["modello"] != riga["modello_reverse"]
        due = [(uscita / "04_error_models" / riga[c]).read_bytes()
               for c in ("modello", "modello_reverse")]
        assert due[0] != due[1]

    inferenza = uscita / "05_asv_inference"
    for nome in ("varianti_per_campione.rds", "varianti_per_campione_reverse.rds",
                 "letture_denoised.tsv", "letture_denoised_reverse.tsv"):
        assert (inferenza / nome).is_file(), nome
    for accession in filtrate:
        for lettura in ("forward", "reverse"):
            assert (inferenza / f"fusione_{lettura}_{accession}.rds").is_file()


def test_la_fusione_rende_conto_di_ogni_coppia(esecuzioni):
    """
    **Obiettivo**: Verificare che per ogni campione le coppie entrate nella
    fusione non superino le letture attribuite a una variante di ciascuna
    lettura, e le coppie fuse non superino quelle entrate; che la quota fusa
    del riepilogo sia quella delle somme e stia sopra la soglia d'arresto
    dichiarata; che la distribuzione delle lunghezze sommi alle coppie fuse; e
    che il tracciamento porti, fra l'inferenza e la tabella, le letture
    reverse e le coppie fuse, mai crescenti da un passo al successivo.

    **Razionale scientifico e sistemistico**: La quota di coppie fuse e' la
    misura con cui si giudica se i troncamenti lasciano sovrapposizione: deve
    essere una frazione delle coppie che potevano fondersi, e ogni lettura
    persa deve risultare nel passo in cui e' stata persa.
    """
    uscita = esecuzioni.prima
    fusione = _righe(uscita / "05_pair_merge" / "fusione.tsv")
    riepilogo = json.loads((uscita / "05_pair_merge" / "fusione.json").read_text())
    config = valida(_dati())
    denoised = {
        lato: dict(zip(_colonna(uscita / "05_asv_inference" / nome, "campione"),
                       map(int, _colonna(uscita / "05_asv_inference" / nome, "letture"))))
        for lato, nome in (("forward", "letture_denoised.tsv"),
                           ("reverse", "letture_denoised_reverse.tsv"))}
    assert len(fusione) == 4
    for riga in fusione:
        coppie, fuse = int(riga["coppie"]), int(riga["fuse"])
        assert 0 < fuse <= coppie <= min(
            denoised["forward"][riga["campione"]], denoised["reverse"][riga["campione"]])
        assert abs(float(riga["quota"]) - fuse / coppie) < 1e-6
        assert int(riga["sequenze"]) > 0
    coppie = sum(int(r["coppie"]) for r in fusione)
    fuse = sum(int(r["fuse"]) for r in fusione)
    assert (riepilogo["coppie_controllate"], riepilogo["fuse_controllate"]) == (coppie, fuse)
    assert abs(riepilogo["quota_controllata"] - fuse / coppie) < 1e-4
    assert riepilogo["quota_controllata"] >= config.pair.min_merged_frac_stop
    assert riepilogo["min_overlap"] == config.pair.min_overlap
    assert riepilogo["lunghezza_attesa"] == config.pair.expected_merged_len
    lunghezze = _righe(uscita / "05_pair_merge" / "lunghezze_fuse.tsv")
    assert sum(int(r["coppie"]) for r in lunghezze) == fuse
    assert riepilogo["lunghezza_minima"] == int(lunghezze[0]["lunghezza"])
    assert riepilogo["lunghezza_massima"] == int(lunghezze[-1]["lunghezza"])
    assert riepilogo["lunghezza_minima"] <= riepilogo["lunghezza_mediana"] <= (
        riepilogo["lunghezza_massima"])

    tracciamento = _righe(uscita / "report" / "tabelle" / "tracciamento_per_campione.tsv")
    colonne = list(tracciamento[0])
    passi = colonne[colonne.index("grezze"):]
    assert passi[passi.index("denoised"):passi.index("tabella") + 1] == [
        "denoised", "denoised_reverse", "fuse", "tabella"]
    per_campione = {r["campione"]: int(r["fuse"]) for r in fusione}
    for riga in tracciamento:
        assert int(riga["fuse"]) == int(riga["tabella"]) == per_campione[riga["accession"]]
        valori = [int(riga[p]) for p in passi if p != "denoised_reverse"]
        assert valori == sorted(valori, reverse=True), riga["accession"]
        assert int(riga["denoised_reverse"]) <= int(riga["filtrate"])


def test_le_varianti_sono_sequenze_fuse_della_lunghezza_ammessa(esecuzioni):
    """
    **Obiettivo**: Verificare che il filtro di lunghezza ammetta le sole
    varianti fra ``pair.len_min`` e ``pair.len_max``, che ne scarti almeno una
    piu' corta, che le sequenze consegnate siano tutte in quell'intervallo
    (piu' lunghe di una lettura: sono fusioni, non letture), e che conteggi,
    tassonomia, sequenze e metadati descrivano le stesse varianti e gli stessi
    quattro campioni.

    **Razionale scientifico e sistemistico**: In modalita' appaiata la
    lunghezza delle varianti e' quella della regione amplificata, e
    l'intervallo ammesso si dichiara: le fusioni fuori intervallo sono
    amplificazioni aspecifiche, e una variante lunga quanto una lettura
    direbbe che la fusione non e' avvenuta.
    """
    uscita = esecuzioni.prima
    config = valida(_dati())
    minima, massima = config.pair.len_min, config.pair.len_max
    lunghezze = _righe(uscita / "07_chimera" / "lunghezze.tsv")
    for riga in lunghezze:
        assert (riga["ammessa"] == "si") == (minima <= int(riga["lunghezza"]) <= massima), riga
    assert any(r["ammessa"] == "no" and int(r["lunghezza"]) < minima for r in lunghezze)

    finale = uscita / "12_final"
    testo = (finale / "sequenze.fasta").read_text(encoding="ascii").splitlines()
    sequenze = [r for r in testo if not r.startswith(">")]
    assert sequenze and all(minima <= len(s) <= massima for s in sequenze)
    assert min(map(len, sequenze)) > config.filter.truncLen
    varianti = _colonna(finale / "conteggi.tsv", "asv_id")
    assert varianti == _colonna(finale / "tassonomia.tsv", "asv_id")
    assert varianti == [r[1:].split()[0] for r in testo if r.startswith(">")]
    with open(finale / "conteggi.tsv", encoding="utf-8") as file:
        campioni = file.readline().rstrip("\n").split("\t")[1:]
    assert len(campioni) == 4 and campioni == _colonna(finale / "metadati.tsv", "accession")


def test_le_varianti_e_la_fusione_sono_quelle_di_dada2(esecuzioni):
    """
    **Obiettivo**: Verificare che le varianti di ciascuna lettura siano
    esattamente quelle di ``dada(pool = "pseudo")`` in una sola chiamata sulle
    letture filtrate con il modello della lettura, e che le sequenze fuse di
    ogni campione siano esattamente quelle di ``mergePairs`` sugli oggetti
    interi di dada2 (gli oggetti ``dada`` e le dereplicazioni dai file), con
    le stesse coppie per sequenza.

    **Razionale scientifico e sistemistico**: La fase di fusione non rilegge
    le letture: riceve dall'inferenza le sole parti degli oggetti di dada2 che
    ``mergePairs`` usa. Il risparmio e' legittimo solo se il risultato e'
    quello della funzione sugli oggetti interi, a ogni aggiornamento di dada2.
    """
    uscita = esecuzioni.prima
    config = valida(_dati())
    righe = _righe(uscita / "04_error_models" / "corrispondenza.tsv")
    campioni = sorted(r["campione"] for r in righe)
    assert len({(r["modello"], r["modello_reverse"]) for r in righe}) == 1
    filtrate, modelli, inferenza = (
        uscita / "03_filtered", uscita / "04_error_models", uscita / "05_asv_inference")
    elenco = ", ".join(f'"{a}"' for a in campioni)
    codice = f"""
        suppressPackageStartupMessages(library(dada2))
        campioni <- c({elenco})
        file <- function(suffisso) stats::setNames(
          file.path("{filtrate}", paste0(campioni, suffisso)), campioni)
        inferisci <- function(suffisso, modello) dada(
          file(suffisso), err = readRDS(file.path("{modelli}", modello)), pool = "pseudo",
          OMEGA_A = {config.dada.omega_a!r}, multithread = 2, verbose = 0)
        f <- inferisci("_filt.fastq.gz", "{righe[0]["modello"]}")
        r <- inferisci("_filt_reverse.fastq.gz", "{righe[0]["modello_reverse"]}")
        uguali <- function(oggetti, nome) {{
          attese <- lapply(oggetti, getUniques)
          scritte <- readRDS(file.path("{inferenza}", nome))
          identical(attese[names(scritte)], scritte)
        }}
        coppie <- mergePairs(
          f, file("_filt.fastq.gz"), r, file("_filt_reverse.fastq.gz"),
          minOverlap = {config.pair.min_overlap}L, maxMismatch = {config.pair.max_mismatch}L,
          trimOverhang = TRUE, verbose = FALSE)
        attese <- lapply(coppie, function(x) {{
          somme <- rowsum(as.numeric(x$abundance), x$sequence)
          v <- stats::setNames(as.integer(somme[, 1L]), rownames(somme))
          v[order(names(v), method = "radix")]
        }})
        fuse <- readRDS("{uscita / "05_pair_merge" / "varianti_per_campione.rds"}")
        fuse <- lapply(fuse, function(v) v[order(names(v), method = "radix")])
        cat(uguali(f, "varianti_per_campione.rds"),
            uguali(r, "varianti_per_campione_reverse.rds"),
            identical(attese[names(fuse)], fuse), length(fuse), sum(lengths(fuse)))
    """
    forward, reverse, fusione, quanti, sequenze = _r(codice).split()
    assert (forward, reverse, fusione) == ("TRUE", "TRUE", "TRUE")
    assert int(quanti) == 4 and int(sequenze) > 0


def test_due_esecuzioni_danno_gli_stessi_byte(esecuzioni):
    """
    **Obiettivo**: Verificare che due esecuzioni della stessa configurazione
    in cartelle diverse diano, per ogni artefatto di ogni fase, lo stesso
    sha256, salvo la configurazione registrata (che riporta la cartella di
    uscita), e che l'esecuzione che toglie le letture filtrate dia gli stessi
    artefatti di quella che le tiene.

    **Razionale scientifico e sistemistico**: La riproducibilita' vale per la
    modalita' appaiata come per l'altra; e ``run.keep_filtered_fastq`` governa
    lo spazio su disco, non il risultato.
    """
    prima, seconda, senza = (esecuzioni.impronte[n] for n in ("prima", "seconda", "senza_filtrate"))
    assert len(prima) > 80
    configurazione = {k for k in prima if k.startswith("00_config/")}
    assert {k: v for k, v in prima.items() if k not in configurazione} == {
        k: v for k, v in seconda.items() if k not in configurazione}
    assert {k: v for k, v in prima.items() if k not in configurazione} == {
        k: v for k, v in senza.items() if k not in configurazione}
    rimossi = json.loads(
        (esecuzioni.senza_filtrate / "03_filtered" / "rimossi_S2.json").read_text())
    assert len(rimossi["rimossi"]) == 8
    assert not list((esecuzioni.senza_filtrate / "03_filtered").glob("*_filt*.fastq.gz"))


def test_la_fusione_si_rifa_senza_rifare_l_inferenza(esecuzioni):
    """
    **Obiettivo**: Verificare che, tolto il manifesto della fusione da
    un'esecuzione conclusa che ha tenuto le letture filtrate, risultino da
    eseguire la fusione e le fasi che ne dipendono e nessuna di quelle che la
    precedono; che la ripresa non riscriva i manifesti fino
    all'inferenza; e che ogni artefatto abbia lo stesso sha256 di prima.

    **Razionale scientifico e sistemistico**: L'inferenza delle due letture e'
    la parte lunga della catena: cambiare la sovrapposizione minima o
    riprendere dopo un'interruzione nella fusione non deve ripeterla.
    """
    ripresa = esecuzioni.ripresa
    assert ripresa.da_eseguire[0] == "SM" and "S5" in ripresa.da_eseguire
    assert not set(ripresa.da_eseguire) & set(FINO_ALL_INFERENZA)
    assert not ripresa.rifatte & set(FINO_ALL_INFERENZA)
    assert ripresa.rifatte == set(ripresa.manifesti_prima) - set(FINO_ALL_INFERENZA)
    assert ripresa.impronte_dopo == ripresa.impronte_prima


def test_la_fusione_si_rifa_dopo_che_le_letture_filtrate_sono_state_tolte(esecuzioni):
    """
    **Obiettivo**: Verificare che, tolto il manifesto della fusione da
    un'esecuzione conclusa con ``run.keep_filtered_fastq: false``, la ripresa
    si concluda, non rifaccia le fasi che precedono il filtro, e dia per ogni
    artefatto lo stesso sha256 di prima; e che le letture filtrate, rigenerate
    per le fasi che le richiedono, siano di nuovo tolte alla fine.

    **Razionale scientifico e sistemistico**: Una ripresa non deve dipendere
    da cio' che la pulizia ha tolto: cio' che manca si rigenera dalle letture
    grezze, e il risultato e' lo stesso byte per byte.
    """
    ripresa = esecuzioni.ripresa_senza_filtrate
    assert "SM" in ripresa.da_eseguire
    assert not ripresa.rifatte & set(PRIMA_DEL_FILTRO)
    assert "SM" in ripresa.rifatte and "S14" in ripresa.rifatte
    assert ripresa.impronte_dopo == ripresa.impronte_prima
    assert ripresa.impronte_dopo == {
        k: v for k, v in esecuzioni.impronte["senza_filtrate"].items()}
    assert not list((ripresa.uscita / "03_filtered").glob("*_filt*.fastq.gz"))
    assert (ripresa.uscita / "03_filtered" / "rimossi_S2.json").is_file()



def _copia_con_altre_soglie(esecuzioni: Esecuzioni, tmp_path: Path, **coppie) -> tuple[Path, Path]:
    """Una copia dell'esecuzione conclusa, e la sua configurazione con altri
    valori del gruppo ``pair``.
    """
    uscita = tmp_path / "copia"
    shutil.copytree(esecuzioni.prima, uscita)
    dati = _dati()
    dati["io"]["out_root"] = str(uscita)
    dati["pair"].update(coppie)
    config = tmp_path / "copia.yaml"
    config.write_text(yaml.safe_dump(dati, sort_keys=False), encoding="utf-8")
    return uscita, config


def test_una_quota_di_coppie_fuse_sotto_la_soglia_ferma_la_fusione(esecuzioni, tmp_path):
    """
    **Obiettivo**: Verificare che, alzata la soglia d'arresto sopra la quota
    di coppie fuse misurata, la ripresa si fermi alla fusione con ``E-SM-01``
    (uscita 4) senza rifare l'inferenza, che il punto di ripresa riporti la
    quota, e che il riepilogo della fusione resti scritto con le misure.

    **Razionale scientifico e sistemistico**: Una quota bassa dice che le
    letture troncate non si sovrappongono: proseguire darebbe una tabella con
    i soli taxa dalla regione corta. Chi legge l'arresto deve trovare le
    misure che lo hanno causato.
    """
    uscita, config = _copia_con_altre_soglie(
        esecuzioni, tmp_path, min_merged_frac_warn=0.999, min_merged_frac_stop=0.999)
    inferenza = _manifesti(uscita)["S4"].read_bytes()
    chiudi()
    assert cli.main(["resume", "--config", str(config)]) == 4
    chiudi()
    punto = json.loads((uscita / "99_logs" / "punto_di_ripresa.json").read_text())
    assert (punto["passo"], punto["codice"]) == ("SM", "E-SM-01")
    assert "pair.min_merged_frac_stop" in punto["dettaglio"]
    riepilogo = json.loads((uscita / "05_pair_merge" / "fusione.json").read_text())
    assert 0 < riepilogo["quota_controllata"] < 0.999
    assert _manifesti(uscita)["S4"].read_bytes() == inferenza
    assert "SM" not in _manifesti(uscita)


def test_gli_avvisi_della_fusione_restano_nel_manifesto(esecuzioni, tmp_path):
    """
    **Obiettivo**: Verificare che, con la soglia d'avviso sopra la quota
    misurata e una lunghezza attesa lontana dalla mediana, la ripresa si
    concluda (uscita 0) rifacendo la fusione e le fasi successive ma non
    l'inferenza, che il manifesto della fusione registri ``E-SM-02`` ed
    ``E-SM-03``, e che i file consegnati siano gli stessi byte.

    **Razionale scientifico e sistemistico**: I due avvisi non cambiano il
    risultato: dicono che i troncamenti o la regione dichiarata vanno
    riguardati. Devono restare negli artefatti, dove li trova chi legge il
    report, e non soltanto nel log.
    """
    uscita, config = _copia_con_altre_soglie(
        esecuzioni, tmp_path, min_merged_frac_warn=0.999, expected_merged_len=400)
    inferenza = _manifesti(uscita)["S4"].read_bytes()
    chiudi()
    assert cli.main(["resume", "--config", str(config)]) == 0
    chiudi()
    manifesto = json.loads(_manifesti(uscita)["SM"].read_text(encoding="utf-8"))
    assert sorted(d["codice"] for d in manifesto["degradazioni"]) == ["E-SM-02", "E-SM-03"]
    assert _manifesti(uscita)["S4"].read_bytes() == inferenza
    for nome in FINALI:
        assert (uscita / "12_final" / nome).read_bytes() == (
            esecuzioni.prima / "12_final" / nome).read_bytes(), nome


def test_un_intervallo_di_lunghezze_senza_sequenze_fuse_ferma_il_filtro_di_lunghezza(
    esecuzioni, tmp_path
):
    """
    **Obiettivo**: Verificare che, dichiarato un intervallo di lunghezze in cui
    non cade alcuna sequenza fusa, la ripresa rifaccia il solo filtro di
    lunghezza e si fermi con ``E-S7-02`` (uscita 4), il codice della modalita'
    appaiata, il cui testo rimanda a ``pair.len_min`` e ``pair.len_max``.

    **Razionale scientifico e sistemistico**: In modalita' appaiata
    l'intervallo ammesso si dichiara, e sbagliarlo svuota la tabella: chi
    legge l'arresto deve trovare i parametri che lo governano, non quelli
    della modalita' single-end.
    """
    from amplicon16s.errors.catalog import CATALOGO

    uscita, config = _copia_con_altre_soglie(esecuzioni, tmp_path, len_min=400, len_max=410)
    prima = {p: m.read_bytes() for p, m in _manifesti(uscita).items()}
    chiudi()
    assert cli.main(["resume", "--config", str(config)]) == 4
    chiudi()
    punto = json.loads((uscita / "99_logs" / "punto_di_ripresa.json").read_text())
    assert (punto["passo"], punto["codice"]) == ("S7", "E-S7-02")
    assert "400-410" in punto["dettaglio"]
    for passo in ("S4", "SM", "S5", "S6"):
        assert _manifesti(uscita)[passo].read_bytes() == prima[passo], passo
    assert "pair.len_min" in CATALOGO["E-S7-02"].azione
    assert "asv.len_tol" not in CATALOGO["E-S7-02"].azione


def test_soglie_propone_le_soglie_della_fusione_dalla_misura(esecuzioni, tmp_path, capsys):
    """
    **Obiettivo**: Verificare che ``soglie --accetta-proposte`` su
    un'esecuzione appaiata conclusa scriva un file valido in cui le due soglie
    sulla quota di coppie fuse sono proposte dalla misura, non piu' strette di
    essa e con l'avviso non sotto l'arresto, senza modificare l'esecuzione.

    **Razionale scientifico e sistemistico**: Le soglie della modalita'
    appaiata seguono lo stesso percorso delle altre: esplorative alla prima
    esecuzione, poi fissate sulla misura dei propri dati.
    """
    radice = esecuzioni.prima.parent
    nuovo = tmp_path / "con_soglie.yaml"
    prima = _impronte(esecuzioni.prima)
    chiudi()
    assert cli.main(["soglie", "--config", str(radice / "prima.yaml"), "--output", str(nuovo),
                     "--accetta-proposte"]) == 0
    chiudi()
    assert "pair.min_merged_frac_stop" in capsys.readouterr().out
    config = carica(nuovo)
    quota = json.loads(
        (esecuzioni.prima / "05_pair_merge" / "fusione.json").read_text())["quota_controllata"]
    assert 0 < config.pair.min_merged_frac_stop <= config.pair.min_merged_frac_warn <= quota
    assert _impronte(esecuzioni.prima) == prima
