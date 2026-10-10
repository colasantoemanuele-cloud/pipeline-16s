"""Dataset senza controlli, e controlli negativi senza letture.

Verifica due proprieta' legittime di un dataset che la pipeline deve dichiarare
invece di respingere: l'assenza di una colonna che distingua i campioni
(``ctrl.column`` nullo: ogni campione e' biologico) e un controllo negativo il
cui file non contiene alcuna lettura (un bianco che non ha amplificato nulla).
I due casi si provano su un dataset sintetico (``tests/sintetico.py``) che
attraversa tutte le fasi, con i codici con cui ciascuna dichiara cio' che non
ha potuto fare.

I test che eseguono la catena richiedono R con i pacchetti di calcolo:

    docker run --rm -e PYTHONPATH=/app/src -e AMPLICON16S_R_DIR=/app/R \\
      -v "$(pwd)":/app -w /app <immagine> pytest tests/test_controlli_assenti.py
"""

from __future__ import annotations

import csv
import functools
import os
from pathlib import Path

import pytest
import sintetico
from sottoinsieme import motivo_pacchetti_r_assenti

from amplicon16s.config.schema import ErroreConfigurazione, valida
from amplicon16s.errors.catalog import CATALOGO
from amplicon16s.io_layer.artifacts import Fase
from amplicon16s.logging.logger import chiudi
from amplicon16s.metadata.models import ClasseCampione
from amplicon16s.report.builder import genera
from amplicon16s.runner.executor import Conclusione, Esecutore
from amplicon16s.runner.graph import Passo
from amplicon16s.runner.project import ProjectRun
from amplicon16s.steps.s00_validate import esegui_s0
from amplicon16s.steps.s14_finale import NOME_CONTROLLI


@pytest.fixture(autouse=True)
def uscite_pulite():
    """Chiude le uscite del log prima e dopo ogni test."""
    chiudi()
    yield
    chiudi()


@functools.cache
def _sonda_bioc() -> str | None:
    """Perche' l'ambiente R richiesto non c'e', o ``None``."""
    return motivo_pacchetti_r_assenti(
        "dada2", "ShortRead", "jsonlite", "phyloseq", "Biostrings", "decontam"
    )


@pytest.fixture
def bioc():
    """Richiede R con i pacchetti di tutte le fasi: salta senza, ma in CI fallisce."""
    if _sonda_bioc() is not None:
        if os.environ.get("AMPLICON16S_RICHIEDI_BIOC") == "1":
            pytest.fail(f"Bioconductor e' richiesto in questo ambiente: {_sonda_bioc()}")
        pytest.skip(_sonda_bioc())


def _tsv(percorso: Path) -> dict[str, dict[str, str]]:
    """Una tabella del tracciamento, per campione."""
    with open(percorso, encoding="utf-8", newline="") as file:
        return {r["campione"]: r for r in csv.DictReader(file, delimiter="\t")}


def _degradazioni(run: ProjectRun) -> dict[str, list[str]]:
    """Per fase, i codici delle degradazioni registrate nel manifesto."""
    return {
        str(passo): [d["codice"] for d in
                     run.albero.manifesto_passo(passo, run.fase(passo).cartella).degradazioni]
        for passo in Passo if passo is not Passo.S9
    }


def test_senza_colonna_del_tipo_le_etichette_si_dichiarano_vuote(tmp_path):
    """
    **Obiettivo**: Verificare che ``ctrl.column`` nullo sia accettato con le tre
    etichette vuote e respinto se una di esse ha valori o se e' dichiarata una
    riclassificazione; che con la colonna dichiarata le etichette dei
    biologici restino obbligatorie; e che S0 costruisca allora un inventario di
    soli biologici, dichiarando l'assenza dei controlli con ``E-S0-17``.

    **Razionale scientifico e sistemistico**: Un dataset senza controlli non
    ha una colonna che li distingua: obbligarlo a indicarne una fittizia
    metterebbe nella configurazione un valore senza significato. L'assenza si
    dichiara, e le fasi che dipendono dai controlli la riportano con il loro
    codice.
    """
    dataset = sintetico.genera(tmp_path / "dati", biologici=3, letture=40, con_tipo=False)
    dati = sintetico.configurazione(dataset, tmp_path / "out")
    config = valida(dati)
    assert config.ctrl.column is None and config.ctrl.biological_values == []

    for chiave, valore in (("biological_values", ["x"]), ("blank_values", ["x"]),
                           ("positive_values", ["x"])):
        sbagliata = sintetico.configurazione(dataset, tmp_path / "out", ctrl={chiave: valore})
        with pytest.raises(ErroreConfigurazione, match="column e' nullo"):
            valida(sbagliata)
    con_colonna = sintetico.configurazione(
        dataset, tmp_path / "out", ctrl={"column": sintetico.COLONNA_SITO})
    with pytest.raises(ErroreConfigurazione, match="ctrl.biological_values"):
        valida(con_colonna)

    risultato = esegui_s0(config)
    assert risultato.superata
    assert {c.classe for c in risultato.inventario} == {ClasseCampione.BIOLOGICO}
    assert {c.materiale for c in risultato.inventario} == {""}
    avvisi = [a for a in risultato.avvisi if a.codice == "E-S0-17"]
    assert len(avvisi) == 2 and all("ctrl.column e' nullo" in a.dettaglio for a in avvisi)


def test_il_file_vuoto_di_un_controllo_negativo_supera_s0_e_si_dichiara(tmp_path):
    """
    **Obiettivo**: Verificare che un controllo negativo il cui file non
    contiene alcuna lettura superi G13 con l'avviso ``E-S0-20``, che resti
    nell'inventario come controllo negativo, e che lo stesso file attribuito a
    un campione biologico fermi S0 con ``E-S0-13``.

    **Razionale scientifico e sistemistico**: Nei campioni a bassa biomassa un
    bianco di estrazione senza letture attesta l'assenza di contaminazione:
    respingerlo obbligherebbe a togliere dal dataset proprio il controllo
    riuscito. Per un campione da cui ci si attende il segnale, invece, zero
    letture restano un difetto del file.
    """
    dataset = sintetico.genera(tmp_path / "dati", biologici=3, negativi=1,
                               negativi_vuoti=1, letture=40, letture_negativi=10)
    (vuoto,) = dataset.vuoti
    risultato = esegui_s0(valida(sintetico.configurazione(dataset, tmp_path / "out")))
    assert risultato.superata
    (avviso,) = [a for a in risultato.avvisi if a.codice == "E-S0-20"]
    assert vuoto in avviso.dettaglio
    assert risultato.inventario[vuoto].classe is ClasseCampione.CONTROLLO_NEGATIVO
    assert CATALOGO["E-S0-20"].categoria.value == "degradazione_automatica"

    # Lo stesso file, dichiarato biologico: un difetto, non un bianco pulito.
    tabella = dataset.tabella.read_text(encoding="utf-8").replace(
        sintetico.ETICHETTA_NEGATIVO, sintetico.ETICHETTA_BIOLOGICO)
    dataset.tabella.write_text(tabella, encoding="utf-8")
    chiudi()
    config = valida(sintetico.configurazione(
        dataset, tmp_path / "out2", ctrl={"blank_values": []}))
    fallita = esegui_s0(config, solleva=False)
    assert not fallita.superata
    assert [v.codice for e in fallita.falliti for v in e.violazioni] == ["E-S0-13"]


def test_un_controllo_negativo_senza_letture_attraversa_tutte_le_fasi(bioc, tmp_path):
    """
    **Obiettivo**: Verificare su un dataset sintetico che la catena S0-S14 si
    concluda con un controllo negativo il cui file e' vuoto: il campione ha
    zero letture in ogni tabella del tracciamento, non ha un profilo ne' un
    file filtrato, e' consegnato fra i controlli; nei manifesti ci sono solo
    codici del catalogo, fra cui ``E-S0-20``, e nessuno del ponte verso R; la
    decontaminazione si esegue sui controlli negativi con letture; il report
    si genera.

    **Razionale scientifico e sistemistico**: Accettare il file in S0 senza
    provarlo fino in fondo sposterebbe l'arresto in una fase successiva, dopo
    il calcolo. Il campione senza letture e' un caso che ogni fase deve
    trattare come tale, come gia' fa per un bianco azzerato dal filtro.
    """
    dataset = sintetico.genera(tmp_path / "dati", biologici=8, negativi=3, negativi_vuoti=1)
    (vuoto,) = dataset.vuoti
    run = ProjectRun(valida(sintetico.configurazione(dataset, tmp_path / "out")))
    esito = Esecutore(run).esegui()
    assert esito.conclusione is Conclusione.COMPLETATA, esito.punto

    degradazioni = _degradazioni(run)
    assert "E-S0-20" in degradazioni["S0"]
    tutti = [c for codici in degradazioni.values() for c in codici]
    assert not [c for c in tutti if c not in CATALOGO or c.startswith("E-R-")], degradazioni
    assert "E-S12-03" not in degradazioni["S12"]

    profili = run.albero.cartella(Fase.QC_PROFILES)
    filtrati = run.albero.cartella(Fase.FILTERED)
    assert _tsv(profili / "letture_grezze.tsv")[vuoto]["letture"] == "0"
    assert _tsv(profili / "coppie.tsv")[vuoto]["letture"] == "0"
    assert _tsv(filtrati / "letture_prefiltro.tsv")[vuoto]["letture"] == "0"
    assert _tsv(filtrati / "letture_filtrate.tsv")[vuoto]["letture"] == "0"
    assert not (filtrati / f"{vuoto}_filt.fastq.gz").exists()
    # Gli altri campioni hanno il loro ordine e le loro letture.
    assert list(_tsv(profili / "letture_grezze.tsv")) == sorted(dataset.campioni)
    assert int(_tsv(profili / "letture_grezze.tsv")["LIB001"]["letture"]) == 1500

    finale = run.albero.cartella(Fase.FINAL)
    assert (finale / NOME_CONTROLLI).is_file() and (finale / "ps_final.rds").is_file()
    with open(finale / "conteggi.tsv", encoding="utf-8") as file:
        intestazione = file.readline().rstrip("\n").split("\t")
    assert vuoto not in intestazione and "LIB001" in intestazione
    chiudi()
    assert "Risultato finale" in genera(run.config.io.out_root).read_text(encoding="utf-8")


def test_un_dataset_senza_colonna_del_tipo_attraversa_tutte_le_fasi(bioc, tmp_path):
    """
    **Obiettivo**: Verificare su un dataset sintetico senza colonna del tipo
    di campione (``ctrl.column`` nullo) che la catena S0-S14 si concluda, che
    l'assenza dei controlli sia dichiarata con i codici esistenti (``E-S0-17``
    in S0, ``E-S11-05`` in S11, ``E-S12-03`` in S12), che i controlli non siano
    fra i file consegnati e che tutti i campioni siano nell'oggetto finale.

    **Razionale scientifico e sistemistico**: Senza controlli la pipeline non
    puo' stimare una soglia di profondita' ne' i contaminanti: prosegue, e lo
    dichiara con i codici con cui lo dichiarerebbe per un dataset che ha la
    colonna ma nessun controllo. Nessun codice nuovo e nessun percorso a parte.
    """
    dataset = sintetico.genera(tmp_path / "dati", biologici=8, con_tipo=False)
    # Con "none" l'assenza della soglia sarebbe una scelta dichiarata, non una
    # degradazione: qui la si chiede ai controlli, che non ci sono.
    run = ProjectRun(valida(sintetico.configurazione(
        dataset, tmp_path / "out", qc={"min_reads_mode": "katharoseq_if_available"})))
    esito = Esecutore(run).esegui()
    assert esito.conclusione is Conclusione.COMPLETATA, esito.punto

    degradazioni = _degradazioni(run)
    assert degradazioni["S0"].count("E-S0-17") == 2
    assert "E-S11-05" in degradazioni["S11"] and "E-S12-03" in degradazioni["S12"]
    tutti = [c for codici in degradazioni.values() for c in codici]
    assert not [c for c in tutti if c not in CATALOGO or c.startswith("E-R-")], degradazioni

    finale = run.albero.cartella(Fase.FINAL)
    assert not (finale / NOME_CONTROLLI).exists()
    with open(finale / "conteggi.tsv", encoding="utf-8") as file:
        intestazione = file.readline().rstrip("\n").split("\t")
    assert sorted(set(intestazione) & set(dataset.campioni)) == sorted(dataset.campioni)
