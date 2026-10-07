r"""Suite di test della settimana 27: generalita' degli ingressi, configurazione e S0.

1. Inquadramento nel Piano Operativo
------------------------------------
Settimana 27 (W27), Fase F8 (un dataset 16S single-end diverso da OSD-734
supera la fase S0, o vi fallisce per la ragione giusta, prima di qualunque
calcolo).

2. Moduli sorgente coperti
--------------------------
* ``src/amplicon16s/config/defaults.py`` (``OBBLIGATORI``, ``ESEMPIO_OSD734``,
  ``FATTI_OSD734``), ``src/amplicon16s/config/schema.py``,
  ``src/amplicon16s/config/resolve.py`` (``registra_risolta``)
* ``src/amplicon16s/gates/g01_g15.py``, ``src/amplicon16s/gates/registry.py``
* ``src/amplicon16s/metadata/tabelle.py``, ``src/amplicon16s/metadata/crosswalk.py``
* ``src/amplicon16s/io_layer/reads.py``
* ``src/amplicon16s/steps/s00_validate.py``, ``src/amplicon16s/cli.py``,
  ``src/amplicon16s/steps/s02_filter.py`` (``archivio_incompleto``),
  ``src/amplicon16s/steps/s11_controls.py`` (``colonna_dei_livelli``)
* ``src/amplicon16s/report/builder.py`` (parametri tarati e non dichiarati)
* ``config/config.example.yaml``, ``dati/osd734/config_osd734.yaml``,
  ``dati/osd276/`` (``scarica_letture.py``)
* ``src/amplicon16s/errors/catalog.py`` (``E-G15-14``, ``E-S0-19``)

3. Cosa valuta questo file
--------------------------
- i parametri che descrivono il dataset sono obbligatori: una configurazione
  che non li dichiara e' respinta da G15 con ``E-G15-10`` e l'elenco completo
  dei mancanti; un parametro dichiarato nullo o vuoto e' dichiarato; G15
  respinge le colonne del lotto senza il file (``E-G15-11``) e i controlli
  positivi senza taxon o colonna delle cellule (``E-G15-12``);
- l'ordine dei gate e' quello delle dipendenze: un'etichetta non dichiarata e'
  diagnosticata da G11 e non da G10;
- G10 con ``filter.trimLeft``: primer cercato solo con il taglio a zero, motivo
  cercato alla posizione del taglio, codici distinti per primer presente
  (``E-S0-10``) e segnale assente (``E-S0-16``), motivo nullo dichiarato
  (``E-S0-18``);
- S0 verifica ogni colonna dichiarata (G02 sulle tabelle di assay e di studio,
  G08 sul file del lotto) e respinge un file del lotto incompleto o ambiguo;
- un dataset con nomi di colonne, etichette e chiavi dei campioni diversi da
  OSD-734 supera S0: tabella di studio facoltativa, colonna identificativa
  distinta, colonne del modulo e della piastra nulle, file nominati per
  campione;
- lettura dei metadati con BOM, intestazioni fra virgolette e fine riga CRLF;
  etichette senza distinzione di maiuscole, anche nella disgiunzione dello
  schema; marcatori delle letture inverse ancorati prima dell'estensione;
  FASTQ compressi o no riconosciuti dai primi byte;
- S0 dichiara i controlli che il dataset non ha (``E-S0-17``);
  ``run.threads`` e' nullo per difetto, i thread effettivi sono i processori
  utilizzabili contati all'uso, e il digest della configurazione e' lo stesso
  con 4 e con 16 processori; ``run.container`` e' facoltativo;
- senza ``meta.module_column`` il modulo viene dal file del lotto; G15 respinge
  un parametro dichiarato senza quello da cui dipende (``E-G15-14``); un
  gruppo mancante elenca tutti i suoi parametri senza predefinito;
- G10 cerca il primer nelle sole classi controllate; G08 dichiara le righe del
  lotto senza campione (``E-S0-19``);
- lo script che ricava le letture forward del secondo dataset conta i record
  ``/1``, ``/2`` e gli altri, mette da parte un file non conforme e raccoglie i
  guasti nel riepilogo;
- il report elenca solo i parametri tarati su OSD-734 presi per difetto, e
  dichiararli al loro valore li toglie dall'elenco anche a una ripresa;
- le configurazioni di OSD-734 e del secondo dataset dichiarano tutto, e nessun
  valore del secondo dataset sta nel codice.

4. Comandi Bash e scenari di esecuzione
---------------------------------------
    ``<immagine>`` e' l'immagine del container della pipeline; quella corrente
    e' indicata in ``test.txt``, sezione 1.3.

    1. Modalità locale standard (R di base con jsonlite, senza Bioconductor né
       dati reali):
       pytest tests/test_w27_ingressi.py -v

    2. Modalità container Docker standard (sottoinsieme ridotto con
       R/Bioconductor):
       docker run --rm \
         -e PYTHONPATH=/app/src \
         -e AMPLICON16S_R_DIR=/app/R \
         -v "$(pwd)":/app \
         -w /app \
         <immagine> \
         pytest -o cache_dir=/tmp/.pytest_cache tests/test_w27_ingressi.py -v

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
         pytest -o cache_dir=/tmp/.pytest_cache tests/test_w27_ingressi.py -v

5. Risultato atteso
-------------------
Vedi ``test.txt``, scheda W27.

6. Razionale scientifico e sistemistico
---------------------------------------
- Un valore ereditato da un altro dataset (una colonna, un'etichetta, un
  troncamento) produce un risultato plausibile e sbagliato: i parametri che
  descrivono il dataset non devono avere un predefinito.
- Un errore nella descrizione del dataset deve fermare l'esecuzione in S0, in
  secondi e con la diagnosi giusta, non dopo ore di calcolo o con il messaggio
  di un altro controllo.
"""

from __future__ import annotations

import copy
import gzip
import json
import os
import re
from pathlib import Path
from typing import Any

import pytest
import yaml
from conftest import (
    INIZIO_CON_MOTIVO,
    INIZIO_CON_PRIMER,
    NEGATIVO,
    POSITIVO,
    Campione,
    crea_scenario,
    dichiarazione_minima,
    lettura,
    parametri_osd734,
    scrivi_fastq,
)
from sottoinsieme import dati_esempio

import amplicon16s.cli as cli
from amplicon16s.config import defaults
from amplicon16s.config.resolve import parametri_dichiarati, risolvi, scrivi_risolta
from amplicon16s.config.schema import (
    INTESTAZIONE_MANCANTI,
    Config,
    ErroreConfigurazione,
    carica,
    obbligatori_mancanti,
    processori_disponibili,
    thread_effettivi,
    valida,
)
from amplicon16s.errors.catalog import CATALOGO, Categoria
from amplicon16s.gates.g01_g15 import Contesto, ErroreGate, esegui_g15
from amplicon16s.gates.registry import esegui_gate, esegui_tutti, nomi_dei_gate
from amplicon16s.io_layer.artifacts import Fase
from amplicon16s.io_layer.reads import scansiona_file
from amplicon16s.logging.logger import chiudi
from amplicon16s.metadata.crosswalk import (
    AccessionNonEstraibile,
    chiave_dalla_tabella,
)
from amplicon16s.metadata.models import ClasseCampione
from amplicon16s.metadata.tabelle import intestazione, leggi_tsv, tabella_di_studio
from amplicon16s.report.builder import CARTELLA_REPORT, NOME_REPORT
from amplicon16s.runner.graph import Passo
from amplicon16s.steps.s00_validate import ValidazioneIngressi, esegui_s0
from amplicon16s.steps.s02_filter import archivio_incompleto

RADICE = Path(__file__).resolve().parents[1]
PRIMER = "GTGYCAGCMGCCGCGGTAA"
ORDINE_DELLE_DIPENDENZE = (
    "G15", "G01", "G02", "G04", "G05", "G06", "G03", "G11", "G13", "G07", "G08", "G09",
    "G10", "G12", "G14",
)


@pytest.fixture(autouse=True)
def uscite_pulite():
    """Chiude le uscite del log prima e dopo ogni test, perche' nessun handler
    resti aperto sulla cartella temporanea.
    """
    chiudi()
    yield
    chiudi()


def _campioni() -> list[Campione]:
    """Due biologici, un positivo e un negativo, con la forma di OSD-734."""
    return [
        Campione("ERX3000001", "NOD1D4.L1"),
        Campione("ERX3000002", "NOD1D4.L2"),
        Campione("ERX3000003", "POS.P1.1", materiale=POSITIVO, posizione="Not Applicable"),
        Campione("ERX3000004", "BLANK.P1.1", materiale=NEGATIVO, posizione="Not Applicable"),
    ]


def _esiti(config: Config) -> dict[str, Any]:
    """Gli esiti di tutti i gate, per nome."""
    return {e.gate: e for e in esegui_tutti(Contesto(config))}


def _primo_fallito(esiti: dict[str, Any]) -> str | None:
    """Il primo gate eseguito e non superato."""
    return next((n for n, e in esiti.items() if e.eseguito and not e.superato), None)


def _con(config: Config, **gruppi: dict[str, Any]) -> Config:
    """La configurazione con i parametri indicati sostituiti, rivalidata."""
    dati = config.model_dump(mode="python")
    for gruppo, valori in gruppi.items():
        dati[gruppo].update(valori)
    return valida(dati)


def _tsv(percorso: Path, righe: list[list[str]], *, bom: bool = False, fine: str = "\n",
         virgolette: bool = False) -> Path:
    """Scrive una tabella separata da tabulazioni, con le varianti di formato
    che un foglio di calcolo puo' produrre.
    """
    testo = fine.join(
        "\t".join(f'"{c}"' if virgolette and i == 0 else c for c in riga)
        for i, riga in enumerate(righe)
    ) + fine
    percorso.write_bytes((b"\xef\xbb\xbf" if bom else b"") + testo.encode("utf-8"))
    return percorso


def _dataset_diverso(radice: Path, *, con_studio: bool = True, **gruppi: dict[str, Any]) -> Config:
    """Un dataset che non ha nulla della forma di OSD-734.

    Sei campioni: file nominati per campione (``S01_R1.fastq.gz``), tabella di
    assay con le colonne ``campione`` e ``chiave``, classi ``tampone``,
    ``bianco`` e ``mock`` in una colonna ``tipo``, nessun file del lotto,
    nessuna posizione. Con ``con_studio`` la classe sta in una tabella di
    studio la cui colonna identificativa ha un altro nome (``id``), scritta con
    BOM e fine riga CRLF; senza, sta nella tabella di assay.
    """
    fastq = radice / "letture"
    fastq.mkdir(parents=True)
    campioni = [("S01", "tampone"), ("S02", "Tampone"), ("S03", "tampone"),
                ("S04", "bianco"), ("S05", "mock"), ("S06", "tampone")]
    for nome, _ in campioni:
        scrivi_fastq(fastq / f"{nome}_R1.fastq.gz", [lettura(INIZIO_CON_MOTIVO, 120)] * 30)
    assay = [["campione", "chiave"] + ([] if con_studio else ["tipo", "nota"])]
    assay += [[f"camp-{n}", n] + ([] if con_studio else [t, "x"]) for n, t in campioni]
    _tsv(radice / "assay.tsv", assay)
    if con_studio:
        _tsv(radice / "studio.tsv",
             [["id", "tipo", "nota"], *[[f"camp-{n}", t, "x"] for n, t in campioni]],
             bom=True, fine="\r\n", virgolette=True)
    riferimento = radice / "rif.fa.gz"
    riferimento.write_bytes(b">s\nACGT\n")
    import hashlib

    dati: dict[str, Any] = {
        "io": {"fastq_dir": str(fastq), "assay_table": str(radice / "assay.tsv"),
               "study_table": str(radice / "studio.tsv") if con_studio else None,
               "out_root": str(radice / "out"), "batch_table": None,
               "accession_regex": r"^(S[0-9]+)_R1"},
        "meta": {"sample_id_column": "campione", "accession_column": "chiave",
                 "study_sample_id_column": "id" if con_studio else None,
                 "module_regex": None, "module_column": None, "non_surface_positions": [],
                 "batch_key_column": None, "batch_module_column": None},
        "filter": {"truncLen": 100},
        "err": {"batch_column": None},
        "tax": {"ref_fasta": str(riferimento),
                "ref_md5": hashlib.md5(riferimento.read_bytes()).hexdigest(),
                "ref_name": "RIF", "ref_version": "1"},
        "ctrl": {"column": "tipo", "blank_values": ["Bianco"], "positive_values": ["MOCK"],
                 "biological_values": ["tampone"], "blank_override_column": None,
                 "blank_override_values": []},
        "katharoseq": {"target_taxon": "Escherichia", "cell_count_column": "nota"},
        "decontam": {"batch_column": None, "min_blanks": 1},
        "qc": {"primer_sequence": PRIMER, "conserved_motif": "TAC[AG].AGG..GC.AGCGTT"},
        "out": {"study_columns": ["nota"], "batch_columns": []},
        "run": {"threads": 1},
    }
    for gruppo, valori in gruppi.items():
        dati.setdefault(gruppo, {}).update(valori)
    return valida(dati)


# --------------------------------------------------------------------------- #
# 1. Parametri obbligatori (criterio C)                                        #
# --------------------------------------------------------------------------- #


def test_senza_i_parametri_obbligatori_g15_respinge_con_l_elenco_dei_mancanti(tmp_path, capsys):
    """
    **Obiettivo**: Verificare che una configurazione con i soli percorsi sia
    respinta da G15 con ``E-G15-10`` e una sola violazione che elenca tutti i
    25 parametri obbligatori; che togliendone tre da una configurazione
    completa siano elencati quei tre e solo quelli; che la riga di comando
    esca con il codice della configurazione non valida, stampando il codice,
    l'elenco e l'azione del catalogo, senza scrivere nulla.

    **Razionale scientifico e sistemistico**: Chi descrive un dataset nuovo
    deve vedere in una volta tutto cio' che manca: un parametro alla volta
    costerebbe venticinque avvii, e un predefinito al suo posto darebbe un
    risultato plausibile su una colonna sbagliata.
    """
    soli_percorsi = {
        "io": {"fastq_dir": "a", "assay_table": "b", "out_root": str(tmp_path / "out")},
        "tax": {"ref_fasta": "c", "ref_md5": "0" * 32},
    }
    assert obbligatori_mancanti(soli_percorsi) == list(defaults.OBBLIGATORI)
    assert len(defaults.OBBLIGATORI) == 25
    with pytest.raises(ErroreGate) as rifiuto:
        esegui_g15(soli_percorsi)
    assert rifiuto.value.gate == "G15" and rifiuto.value.codice == "E-G15-10"
    (violazione,) = rifiuto.value.violazioni
    assert violazione.dettaglio.startswith(f"{INTESTAZIONE_MANCANTI} (25): ")
    for chiave in defaults.OBBLIGATORI:
        assert chiave in violazione.dettaglio, chiave
    assert CATALOGO["E-G15-10"].categoria is Categoria.REVISIONE_UMANA

    completa = dati_esempio()
    tolti = ("ctrl.column", "filter.truncLen", "out.batch_columns")
    for chiave in tolti:
        gruppo, nome = chiave.split(".")
        del completa[gruppo][nome]
    with pytest.raises(ErroreGate) as rifiuto:
        esegui_g15(completa)
    assert rifiuto.value.violazioni[0].dettaglio.startswith(
        f"{INTESTAZIONE_MANCANTI} (3): ctrl.column, out.batch_columns, filter.truncLen."
    )
    assert len(rifiuto.value.violazioni) == 1

    percorso = tmp_path / "config.yaml"
    percorso.write_text(yaml.safe_dump(soli_percorsi), encoding="utf-8")
    assert cli.main(["validate", "--config", str(percorso)]) == cli.USCITA_CONFIGURAZIONE
    uscita = capsys.readouterr().out
    assert "G15 ha respinto la configurazione" in uscita and "[E-G15-10]" in uscita
    assert "meta.accession_column" in uscita and "qc.primer_sequence" in uscita
    assert CATALOGO["E-G15-10"].azione in uscita
    assert not (tmp_path / "out").exists()


def test_un_parametro_obbligatorio_dichiarato_nullo_o_vuoto_e_dichiarato(tmp_path):
    """
    **Obiettivo**: Verificare che i parametri obbligatori non pertinenti si
    dichiarino nulli o vuoti (etichette dei controlli, taxon atteso, colonna
    delle cellule, motivo conservato, colonne del lotto e del modulo), che
    ``ctrl.biological_values`` non possa essere vuoto, e che ogni parametro
    di ``OBBLIGATORI`` sia senza predefinito nello schema.

    **Razionale scientifico e sistemistico**: Obbligatorio significa dichiarato
    in modo esplicito, non provvisto di un valore: un dataset senza controlli
    positivi lo dice scrivendo un elenco vuoto, e non puo' ometterlo.
    """
    config = _dataset_diverso(
        tmp_path,
        ctrl={"blank_values": [], "positive_values": []},
        katharoseq={"target_taxon": None, "cell_count_column": None},
        qc={"conserved_motif": None},
    )
    assert config.ctrl.positive_values == [] and config.katharoseq.target_taxon is None
    assert config.meta.module_column is None and config.decontam.batch_column is None
    assert obbligatori_mancanti(config.model_dump(mode="json")) == []
    esegui_g15(config.model_dump(mode="python"))

    dati = config.model_dump(mode="python")
    dati["ctrl"]["biological_values"] = []
    with pytest.raises(ErroreConfigurazione, match="ctrl.biological_values"):
        valida(dati)
    for chiave in defaults.OBBLIGATORI:
        gruppo, nome = chiave.split(".")
        campo = Config.model_fields[gruppo].annotation.model_fields[nome]
        assert campo.is_required(), chiave
    assert not set(defaults.OBBLIGATORI) & set(defaults.FATTI_OSD734)


def test_g15_respinge_le_colonne_del_lotto_senza_il_file_in_s0_non_in_s10(tmp_path, capsys):
    """
    **Obiettivo**: Verificare che senza ``io.batch_table`` e con
    ``out.batch_columns`` non vuoto G15 respinga la configurazione con
    ``E-G15-11`` nominando il parametro; che la riga di comando si fermi ai
    controlli di avvio, prima di qualunque fase, senza alcun manifesto; e che
    lo stesso valga per le altre colonne del lotto.

    **Razionale scientifico e sistemistico**: Una colonna chiesta a un file che
    non c'e' si scopriva in S10, dopo tutte le fasi di calcolo (dodici minuti
    sul secondo dataset, piu' di un'ora su OSD-734): e' un errore della
    configurazione, e si vede leggendo la configurazione.
    """
    scenario = crea_scenario(tmp_path, _campioni(), con_letture=True)
    incoerente = _con(scenario.config, out={"batch_columns": ["well_id"]})
    esito = esegui_gate("G15", Contesto(incoerente))
    assert not esito.superato and {v.codice for v in esito.violazioni} == {"E-G15-11"}
    assert "out.batch_columns" in esito.violazioni[0].dettaglio

    percorso = tmp_path / "c.yaml"
    percorso.write_text(yaml.safe_dump(dichiarazione_minima(incoerente)), encoding="utf-8")
    assert cli.main(["run", "--config", str(percorso)]) != 0
    uscita = capsys.readouterr().out
    assert "E-G15-11" in uscita
    radice = Path(scenario.config.io.out_root)
    assert not list(radice.rglob("manifest_S*.json"))

    for gruppo, nome in (("err", "batch_column"), ("decontam", "batch_column"),
                         ("meta", "batch_key_column"), ("meta", "batch_module_column")):
        altra = _con(scenario.config, **{gruppo: {nome: "colonna"}})
        violazioni = esegui_gate("G15", Contesto(altra)).violazioni
        assert [v.codice for v in violazioni] == ["E-G15-11"], nome
        assert f"{gruppo}.{nome}" in violazioni[0].dettaglio


def test_g15_respinge_i_controlli_positivi_senza_taxon_o_colonna_delle_cellule(tmp_path):
    """
    **Obiettivo**: Verificare che con ``ctrl.positive_values`` non vuoto e
    ``katharoseq.target_taxon`` o ``katharoseq.cell_count_column`` nulli G15
    respinga con ``E-G15-12``, e che con ``ctrl.positive_values`` vuoto gli
    stessi nulli siano accettati.

    **Razionale scientifico e sistemistico**: I controlli positivi dichiarati
    ma non valutabili sarebbero ignorati in silenzio da S11: la coerenza fra i
    parametri si verifica dove costa meno.
    """
    scenario = crea_scenario(tmp_path, _campioni(), con_letture=True)
    for nome in ("target_taxon", "cell_count_column"):
        esito = esegui_gate("G15", Contesto(_con(scenario.config, katharoseq={nome: None})))
        assert [v.codice for v in esito.violazioni] == ["E-G15-12"]
        assert f"katharoseq.{nome}" in esito.violazioni[0].dettaglio
    senza = _con(scenario.config, ctrl={"positive_values": []},
                 katharoseq={"target_taxon": None, "cell_count_column": None})
    assert esegui_gate("G15", Contesto(senza)).superato


# --------------------------------------------------------------------------- #
# 2. Ordine dei gate                                                           #
# --------------------------------------------------------------------------- #


def test_i_gate_seguono_l_ordine_delle_dipendenze(tmp_path):
    """
    **Obiettivo**: Verificare che il registro esegua i gate nell'ordine G15,
    G01, G02, G04, G05, G06, G03, G11, G13, G07, G08, G09, G10, G12, G14, e che
    un'etichetta di classe non dichiarata sia diagnosticata da G11 con
    ``E-S0-11``, con G10 non eseguito.

    **Razionale scientifico e sistemistico**: G10 giudica il segnale sui
    campioni di una classe: se nessun campione ha una classe, il suo messaggio
    parla di primer e motivo mentre l'errore e' un'etichetta. Sul secondo
    dataset la diagnosi sbagliata e' costata un tentativo.
    """
    assert nomi_dei_gate() == ORDINE_DELLE_DIPENDENZE
    campioni = _campioni()
    for campione in campioni:
        campione.materiale = "Cells"
    scenario = crea_scenario(tmp_path, campioni, con_letture=True)
    esiti = _esiti(scenario.config)
    assert _primo_fallito(esiti) == "G11"
    assert {v.codice for v in esiti["G11"].violazioni} == {"E-S0-11"}
    assert "'Cells'" in esiti["G11"].violazioni[0].dettaglio
    assert not esiti["G10"].eseguito and not esiti["G13"].eseguito


# --------------------------------------------------------------------------- #
# 3. G10 e filter.trimLeft (criterio E)                                        #
# --------------------------------------------------------------------------- #


def _con_primer(tmp_path: Path, **gruppi: dict[str, Any]) -> Config:
    """Lo scenario con il primer in testa a ogni lettura, seguito dal motivo."""
    campioni = _campioni()
    for campione in campioni:
        campione.inizio_letture = INIZIO_CON_PRIMER + INIZIO_CON_MOTIVO
    scenario = crea_scenario(tmp_path, campioni, con_letture=True)
    return _con(scenario.config, **gruppi) if gruppi else scenario.config


def test_con_il_primer_nelle_letture_e_trimleft_pari_alla_sua_lunghezza_g10_passa(tmp_path):
    """
    **Obiettivo**: Verificare che con il primer in testa alle letture e
    ``filter.trimLeft`` a zero G10 fallisca con il solo ``E-S0-10``, indicando
    di impostare ``filter.trimLeft`` a 19; che applicata quella correzione G10
    passi, con il motivo conservato trovato alla posizione 19, anche quando
    il motivo e' ancorato all'inizio con ``^``; e che l'intera S0 sia superata.

    **Razionale scientifico e sistemistico**: La correzione che il gate indica
    deve farlo superare: prima cercava primer e motivo sempre all'inizio della
    lettura, e un dataset con il primer nelle letture non poteva passare
    nemmeno dichiarando il taglio.
    """
    senza_taglio = _con_primer(tmp_path / "a")
    esito = _esiti(senza_taglio)["G10"]
    assert [v.codice for v in esito.violazioni] == ["E-S0-10"]
    assert f"filter.trimLeft a {len(PRIMER)}" in esito.violazioni[0].dettaglio

    con_taglio = _con_primer(tmp_path / "b", filter={"trimLeft": len(PRIMER)})
    contesto = Contesto(con_taglio)
    esito = esegui_gate("G10", contesto)
    assert esito.superato and not esito.avvisi
    assert all(s.frazione_motivo == 1.0 and s.con_primer == 0 for s in contesto.scansione.values())
    assert esegui_s0(con_taglio).superata
    ancorato = _con(con_taglio, qc={"conserved_motif": "^" + con_taglio.qc.conserved_motif})
    assert esegui_gate("G10", Contesto(ancorato)).superato


def test_primer_presente_e_segnale_assente_hanno_codici_distinti(tmp_path):
    """
    **Obiettivo**: Verificare che un taglio che non corrisponde a cio' che
    precede il motivo faccia fallire G10 con ``E-S0-16`` (segnale assente),
    nominando la posizione cercata, e non con ``E-S0-10`` (primer presente);
    che i due codici abbiano sintesi e azioni diverse nel catalogo.

    **Razionale scientifico e sistemistico**: Sono due difetti con due rimedi:
    il primer si toglie con il taglio, il segnale assente dice che i file o la
    regione dichiarata sono sbagliati. Sotto lo stesso codice l'azione
    suggerita era giusta per uno solo.
    """
    sbagliato = _con_primer(tmp_path, filter={"trimLeft": 5})
    esito = _esiti(sbagliato)["G10"]
    assert [v.codice for v in esito.violazioni] == ["E-S0-16"]
    assert "alla posizione 5" in esito.violazioni[0].dettaglio
    primer, segnale = CATALOGO["E-S0-10"], CATALOGO["E-S0-16"]
    assert primer.sintesi != segnale.sintesi and primer.azione != segnale.azione
    assert "filter.trimLeft" in primer.azione and "qc.conserved_motif" in segnale.azione
    assert primer.categoria is segnale.categoria is Categoria.REVISIONE_UMANA


def test_senza_motivo_conservato_g10_verifica_il_solo_primer(tmp_path):
    """
    **Obiettivo**: Verificare che con ``qc.conserved_motif`` nullo G10 passi su
    letture senza primer dichiarando ``E-S0-18`` (degradazione, registrata nel
    manifesto di S0), e fallisca comunque con ``E-S0-10`` se il primer e' in
    testa.

    **Razionale scientifico e sistemistico**: Una regione amplificata senza un
    motivo noto non deve fermare la pipeline, ma il controllo che non si e'
    potuto fare va dichiarato: senza, un file privo di segnale passerebbe per
    verificato.
    """
    scenario = crea_scenario(tmp_path / "a", _campioni(), con_letture=True)
    senza = _con(scenario.config, qc={"conserved_motif": None})
    esito = esegui_gate("G10", Contesto(senza))
    assert esito.superato and [a.codice for a in esito.avvisi] == ["E-S0-18"]
    assert CATALOGO["E-S0-18"].categoria is Categoria.DEGRADAZIONE_AUTOMATICA
    risultato = esegui_s0(senza)
    assert risultato.superata
    manifesto = json.loads(next(Path(senza.io.out_root).rglob("manifest_S0.json")).read_text())
    assert "E-S0-18" in {d["codice"] for d in manifesto["degradazioni"]}

    col_primer = _con_primer(tmp_path / "b", qc={"conserved_motif": None})
    esito = esegui_gate("G10", Contesto(col_primer))
    assert [v.codice for v in esito.violazioni] == ["E-S0-10"] and not esito.avvisi


# --------------------------------------------------------------------------- #
# 4. Colonne dichiarate e file del lotto                                       #
# --------------------------------------------------------------------------- #


def test_g02_verifica_ogni_colonna_dichiarata_delle_tabelle(tmp_path):
    """
    **Obiettivo**: Verificare che G02 respinga con ``E-S0-02``, nominando
    parametro e colonna, una colonna di ``out.study_columns`` e la colonna di
    ``ctrl.blank_override_column`` assenti dalla tabella di studio; che la
    colonna di ``katharoseq.cell_count_column``, anche se non e' fra quelle
    richieste per l'oggetto, sia cercata nelle tabelle e respinta se manca da
    tutte; che una colonna ripetuta nell'intestazione sia respinta allo stesso
    modo; e che due colonne richieste che nell'oggetto prenderebbero lo stesso
    nome, o quello di una colonna dell'inventario, siano respinte da G15 con
    ``E-G15-13`` (prima le scopriva S10, a calcolo concluso).

    **Razionale scientifico e sistemistico**: Ogni colonna che la
    configurazione nomina va cercata prima del calcolo: scoperta assente da
    S10 o da S11 costa l'intera catena.
    """
    scenario = crea_scenario(tmp_path, _campioni(), con_letture=True)
    casi = {
        "out.study_columns": {"out": {"study_columns": ["katharoseq_cell_count",
                                                        "Colonna che non c'e'"]}},
        "ctrl.blank_override_column": {"ctrl": {"blank_override_column": "Altra colonna"}},
    }
    non_portata = _con(scenario.config, katharoseq={"cell_count_column": "cellule"})
    assert esegui_gate("G15", Contesto(non_portata)).superato
    esito = _esiti(non_portata)["G02"]
    assert {v.codice for v in esito.violazioni} == {"E-S0-02"}
    assert "katharoseq.cell_count_column" in esito.violazioni[0].dettaglio
    portata = _con(scenario.config, katharoseq={"cell_count_column": "cellule"},
                   out={"study_columns": ["cellule"]})
    assert esegui_gate("G15", Contesto(portata)).superato
    esito = _esiti(portata)["G02"]
    assert {v.codice for v in esito.violazioni} == {"E-S0-02"}
    assert "'cellule'" in esito.violazioni[0].dettaglio
    for parametro, modifica in casi.items():
        esito = _esiti(_con(scenario.config, **modifica))["G02"]
        assert not esito.superato and {v.codice for v in esito.violazioni} == {"E-S0-02"}
        assert parametro in esito.violazioni[0].dettaglio, parametro
        assert "non compare" in esito.violazioni[0].dettaglio

    studio = Path(scenario.config.io.study_table)
    righe = studio.read_text(encoding="utf-8").split("\n")
    righe[0] += "\tnota\tnota"
    studio.write_text("\n".join(r + "\tx\ty" if i and r else r for i, r in enumerate(righe)),
                      encoding="utf-8")
    esito = _esiti(_con(scenario.config,
                        out={"study_columns": ["katharoseq_cell_count", "nota"]}))["G02"]
    assert "compare 2 volte" in esito.violazioni[0].dettaglio

    for colonne in ({"study_columns": ["katharoseq_cell_count", "Classe"]},
                    {"study_columns": ["katharoseq_cell_count", "Var Uno", "Var-Uno"]},
                    {"study_columns": ["katharoseq_cell_count"],
                     "batch_columns": ["katharoseq_cell_count"]}):
        esito = esegui_gate("G15", Contesto(_con(scenario.config, out=colonne)))
        assert "E-G15-13" in [v.codice for v in esito.violazioni], colonne


def test_g08_verifica_le_colonne_e_la_completezza_del_file_del_lotto(tmp_path):
    """
    **Obiettivo**: Verificare che G08 respinga con ``E-S0-08`` una colonna di
    ``out.batch_columns`` o di ``meta.batch_module_column`` assente dal file
    del lotto; un file in cui un campione non ha riga (incompleto); un file in
    cui un campione ha due righe (ambiguo); un campione con la corsa vuota o
    con la piastra vuota. Il file completo passa.

    **Razionale scientifico e sistemistico**: Un file del lotto incompleto o
    ambiguo era un errore di S3, dopo il filtro di tutte le letture: riguarda
    gli ingressi, e gli ingressi si giudicano in S0.
    """
    completo = crea_scenario(tmp_path / "a", _campioni(), con_letture=True, con_arricchimento=True)
    assert _esiti(completo.config)["G08"].superato

    for parametro, modifica in (
        ("out.batch_columns", {"out": {"batch_columns": ["well_id"]}}),
        ("meta.batch_module_column", {"meta": {"batch_module_column": "modulo_assente"}}),
    ):
        esito = _esiti(_con(completo.config, **modifica))["G08"]
        assert {v.codice for v in esito.violazioni} == {"E-S0-08"}
        assert parametro in esito.violazioni[0].dettaglio

    lotto = Path(completo.config.io.batch_table)
    righe = lotto.read_text(encoding="utf-8").rstrip("\n").split("\n")

    lotto.write_text("\n".join(righe[:-1]) + "\n", encoding="utf-8")
    esito = _esiti(completo.config)["G08"]
    assert not esito.superato and "non hanno una riga" in esito.violazioni[0].dettaglio
    assert "ERX3000004" in esito.violazioni[0].dettaglio

    lotto.write_text("\n".join([*righe, righe[1]]) + "\n", encoding="utf-8")
    esito = _esiti(completo.config)["G08"]
    assert "piu' di una riga" in esito.violazioni[0].dettaglio
    assert "ERX3000001 (2 righe)" in esito.violazioni[0].dettaglio

    campioni = _campioni()
    campioni[1].corsa = ""
    senza_corsa = crea_scenario(tmp_path / "b", campioni, con_letture=True, con_arricchimento=True)
    esito = _esiti(senza_corsa.config)["G08"]
    assert {v.codice for v in esito.violazioni} == {"E-S0-08"}
    assert "err.batch_column" in esito.violazioni[0].dettaglio
    assert "ERX3000002" in esito.violazioni[0].dettaglio

    campioni = _campioni()
    campioni[2].piastra = "Not Applicable"
    senza_piastra = crea_scenario(tmp_path / "c", campioni, con_letture=True, con_arricchimento=True)
    esito = _esiti(senza_piastra.config)["G08"]
    assert {v.codice for v in esito.violazioni} == {"E-S0-08"}
    assert "decontam.batch_column" in esito.violazioni[0].dettaglio
    assert "ERX3000003" in esito.violazioni[0].dettaglio


# --------------------------------------------------------------------------- #
# 5. Un dataset che non ha la forma di OSD-734 (criterio E)                    #
# --------------------------------------------------------------------------- #


def test_un_dataset_con_colonne_etichette_e_chiavi_diverse_supera_s0(tmp_path):
    """
    **Obiettivo**: Verificare che un dataset con file nominati per campione,
    chiave estratta da un gruppo di cattura, tabella di studio con una colonna
    identificativa di nome diverso, scritta con BOM, CRLF e intestazione fra
    virgolette, etichette scritte con maiuscole diverse da quelle dichiarate,
    nessun file del lotto e nessuna colonna del modulo superi tutti i quindici
    gate; e che l'inventario abbia le sei chiavi, le classi giuste, e piastra,
    corsa, modulo e posizione vuoti.

    **Razionale scientifico e sistemistico**: E' cio' che la fase S0 deve
    lasciare disponibile: accettare qualunque dataset 16S single-end descritto
    dalla configurazione, senza che nessun nome di colonna, etichetta o formato
    di accession di OSD-734 sia presupposto dal codice.
    """
    config = _dataset_diverso(tmp_path)
    risultato = esegui_s0(config)
    assert risultato.superata, [str(v) for e in risultato.falliti for v in e.violazioni]
    assert [e.gate for e in risultato.esiti] == list(ORDINE_DELLE_DIPENDENZE)
    inventario = risultato.inventario
    assert inventario.accessioni == ("S01", "S02", "S03", "S04", "S05", "S06")
    conteggi = inventario.conteggi()
    assert (conteggi[ClasseCampione.BIOLOGICO], conteggi[ClasseCampione.CONTROLLO_NEGATIVO],
            conteggi[ClasseCampione.CONTROLLO_POSITIVO]) == (4, 1, 1)
    assert inventario["S02"].nome == "camp-S02" and inventario["S02"].materiale == "Tampone"
    assert all(c.piastra is None and c.corsa is None and c.modulo is None and c.posizione is None
               for c in inventario)
    assert inventario["S01"].file.name == "S01_R1.fastq.gz"


def test_senza_tabella_di_studio_la_classe_si_legge_dalla_tabella_di_assay(tmp_path):
    """
    **Obiettivo**: Verificare che con ``io.study_table`` nullo classe e colonne
    da portare nell'oggetto si leggano dalla tabella di assay, che S0 passi, e
    che una colonna dichiarata e assente dalla tabella di assay sia respinta da
    G02 nominando ``io.assay_table``.

    **Razionale scientifico e sistemistico**: Non ogni dataset ha due tabelle:
    quando i metadati stanno in una sola, pretendere la seconda obbligherebbe a
    fabbricarla.
    """
    config = _dataset_diverso(tmp_path / "a", con_studio=False)
    assert tabella_di_studio(config) == (Path(config.io.assay_table), "campione")
    risultato = esegui_s0(config)
    assert risultato.superata
    assert risultato.inventario.conteggi()[ClasseCampione.BIOLOGICO] == 4

    manca = _dataset_diverso(tmp_path / "b", con_studio=False,
                             out={"study_columns": ["nota", "assente"]})
    esito = _esiti(manca)["G02"]
    assert {v.codice for v in esito.violazioni} == {"E-S0-02"}
    assert "io.assay_table" in esito.violazioni[0].dettaglio


def test_la_chiave_del_campione_si_ricava_dal_nome_del_file_e_dalla_colonna():
    """
    **Obiettivo**: Verificare ``chiave_dalla_tabella``: da un valore che e' il
    nome di un file la chiave si estrae con la stessa espressione dei file; un
    valore che e' gia' la chiave vale cosi' com'e' quando l'espressione non vi
    trova nulla; un valore vuoto o con due chiavi e' un errore.

    **Razionale scientifico e sistemistico**: In OSD-734 la colonna riporta il
    nome del file, in altri dataset il nome del campione o la corsa: la stessa
    regola deve servire entrambi i casi senza un parametro in piu'.
    """
    accession = re.compile(r"[ESD]RX[0-9]{4,}")
    assert chiave_dalla_tabella("GLDS-653_Amplicon_ERX12083297_raw.fastq.gz", accession) == "ERX12083297"
    per_campione = re.compile(r"^(S[0-9]+)_R1")
    assert chiave_dalla_tabella("S01", per_campione) == "S01"
    assert chiave_dalla_tabella("S01_R1.fastq.gz", per_campione) == "S01"
    with pytest.raises(AccessionNonEstraibile):
        chiave_dalla_tabella("", per_campione)
    with pytest.raises(AccessionNonEstraibile, match="ambiguo"):
        chiave_dalla_tabella("ERX1000001_ERX1000002", accession)


# --------------------------------------------------------------------------- #
# 6. Lettura dei metadati, etichette, layout e formato dei file                #
# --------------------------------------------------------------------------- #


def test_le_tabelle_si_leggono_con_bom_virgolette_e_crlf(tmp_path):
    """
    **Obiettivo**: Verificare che ``intestazione`` e ``leggi_tsv`` diano gli
    stessi nomi di colonna, ripuliti, per una tabella con BOM, fine riga CRLF,
    virgolette e spazi attorno ai nomi, e per la stessa tabella senza; che una
    colonna senza nome sia ignorata.

    **Razionale scientifico e sistemistico**: Una tabella salvata da un foglio
    di calcolo ha un BOM in testa: senza toglierlo la prima colonna, di solito
    quella del nome del campione, non verrebbe trovata, e l'errore parlerebbe
    di una colonna che a occhio c'e'.
    """
    righe = [["Sample Name", " tipo ", ""], ["A", ' "x" ', "y"]]
    pulita = _tsv(tmp_path / "pulita.tsv", [["Sample Name", "tipo", ""], ["A", "x", "y"]])
    sporca = _tsv(tmp_path / "sporca.tsv", righe, bom=True, fine="\r\n", virgolette=True)
    assert intestazione(sporca) == intestazione(pulita) == ["Sample Name", "tipo", ""]
    assert leggi_tsv(sporca) == leggi_tsv(pulita) == [{"Sample Name": "A", "tipo": "x"}]
    assert intestazione(_tsv(tmp_path / "vuota.tsv", [])) in ([], [""])


def test_le_etichette_si_confrontano_senza_distinzione_di_maiuscole(tmp_path):
    """
    **Obiettivo**: Verificare che lo schema respinga la stessa etichetta
    dichiarata in due classi anche se scritta con maiuscole diverse, e che la
    classificazione riconosca un materiale scritto con maiuscole diverse da
    quelle dichiarate.

    **Razionale scientifico e sistemistico**: La classificazione ignora le
    maiuscole: il controllo di disgiunzione deve usare lo stesso confronto,
    altrimenti "Blank" e "blank" in due elenchi passerebbero lo schema e uno
    dei due vincerebbe in silenzio.
    """
    dati = parametri_osd734()
    with pytest.raises(ValueError, match="due categorie"):
        Config.model_fields["ctrl"].annotation.model_validate(
            {**dati["ctrl"], "blank_values": ["Blank"], "biological_values": ["blank ", "swab"]}
        )
    config = _dataset_diverso(tmp_path)
    inventario = esegui_s0(config).inventario
    assert inventario["S04"].classe is ClasseCampione.CONTROLLO_NEGATIVO
    assert inventario["S05"].classe is ClasseCampione.CONTROLLO_POSITIVO


@pytest.mark.parametrize(("nome", "inversa"), [
    ("ERX3000009_x_R2.fastq.gz", True),
    ("ERX3000009_x.R2.fastq.gz", True),
    ("ERX3000009_x_2.fastq.gz", True),
    ("ERX3000009_x_R2_001.fastq.gz", True),
    ("ERX3000009_LAB_R2D2.fastq.gz", False),
    ("ERX3000009_R2_coda.fastq.gz", False),
    ("ERX3000009_x_R1.fastq.gz", False),
    ("ERX3000009_plate_2_x.fastq.gz", False),
])
def test_g07_riconosce_le_letture_inverse_dal_marcatore_prima_dell_estensione(tmp_path, nome, inversa):
    """
    **Obiettivo**: Verificare che G07 consideri lettura inversa un file il cui
    nome termina, prima dell'estensione, con ``_R2``, ``.R2`` o ``_2`` (anche
    seguito dal numero di blocco), e non un file che contiene ``_R2`` o ``_2``
    in un altro punto del nome.

    **Razionale scientifico e sistemistico**: Il marcatore cercato ovunque nel
    nome respingeva come paired-end un campione chiamato ``LAB_R2D2``: il nome
    del campione e' scelto da chi ha fatto l'esperimento, non dalla pipeline.
    """
    campioni = [*_campioni(), Campione("ERX3000009", "x", file=nome)]
    scenario = crea_scenario(tmp_path, campioni, con_letture=True)
    esito = esegui_gate("G07", Contesto(scenario.config))
    assert esito.superato is (not inversa)
    if inversa:
        assert {v.codice for v in esito.violazioni} == {"E-S0-07"} and nome in str(esito.violazioni[0])


def test_un_fastq_si_riconosce_dai_primi_byte_non_dall_estensione(tmp_path):
    """
    **Obiettivo**: Verificare che ``scansiona_file`` legga allo stesso modo un
    FASTQ compresso, lo stesso file non compresso con l'estensione ``.gz``, e
    lo stesso file compresso con l'estensione ``.fastq``; anche con fine riga
    CRLF.

    **Razionale scientifico e sistemistico**: L'estensione e' un'etichetta, il
    contenuto e' il dato: un archivio decompresso senza rinominarlo non deve
    essere respinto come corrotto.
    """
    sequenze = [lettura(INIZIO_CON_MOTIVO, 100)] * 5
    compresso = tmp_path / "a.fastq.gz"
    scrivi_fastq(compresso, sequenze)
    testo = gzip.decompress(compresso.read_bytes())
    (tmp_path / "non_compresso.fastq.gz").write_bytes(testo)
    (tmp_path / "compresso.fastq").write_bytes(compresso.read_bytes())
    (tmp_path / "crlf.fastq").write_bytes(testo.replace(b"\n", b"\r\n"))
    attese = None
    for nome in ("a.fastq.gz", "non_compresso.fastq.gz", "compresso.fastq", "crlf.fastq"):
        s = scansiona_file(tmp_path / nome, 100, "GTG[CT]CAGC", "TAC[AG].AGG..GC.AGCGTT")
        misure = (s.valido, s.letture_esaminate, s.lunghezza_minima, s.lunghezza_massima,
                  s.con_primer, s.con_motivo)
        attese = attese or misure
        assert misure == attese == (True, 5, 100, 100, 0, 5), nome
        # Il controllo di S2 prima del filtro legge il file allo stesso modo:
        # un FASTQ non compresso non e' un archivio corrotto.
        assert archivio_incompleto(tmp_path / nome) is None, nome
    troncato = tmp_path / "troncato.fastq.gz"
    troncato.write_bytes(compresso.read_bytes()[:-20])
    assert archivio_incompleto(troncato) is not None


# --------------------------------------------------------------------------- #
# 7. Controlli presenti, thread e immagine                                     #
# --------------------------------------------------------------------------- #


def test_s0_dichiara_i_controlli_che_il_dataset_non_ha(tmp_path):
    """
    **Obiettivo**: Verificare che S0 dichiari con ``E-S0-17`` un dataset senza
    controlli positivi e uno con meno controlli negativi di
    ``decontam.min_blanks``, distinguendo le etichette vuote da quelle che non
    riconoscono campioni; che le dichiarazioni siano degradazioni nel manifesto
    di S0 e in ``gates.json`` sotto G11; e che con i controlli presenti non ce
    ne siano.

    **Razionale scientifico e sistemistico**: Senza positivi non c'e' curva di
    calibrazione, senza negativi non c'e' decontaminazione: chi esegue deve
    saperlo dalla validazione, non scoprirlo leggendo l'oggetto finale. Sul
    secondo dataset S0 passava senza dire nulla.
    """
    solo_biologici = [Campione(f"ERX300000{i}", f"NOD1D4.L{i}") for i in range(1, 4)]
    senza = crea_scenario(
        tmp_path / "a", solo_biologici, con_letture=True,
        sovrascrivi={"ctrl": {"positive_values": [], "blank_values": [],
                              "blank_override_values": []},
                     "katharoseq": {"target_taxon": None, "cell_count_column": None}},
    )
    risultato = esegui_s0(senza.config)
    assert risultato.superata
    g11 = next(e for e in risultato.esiti if e.gate == "G11")
    assert [a.codice for a in g11.avvisi] == ["E-S0-17", "E-S0-17"]
    testi = " | ".join(a.dettaglio for a in g11.avvisi)
    minimo = defaults.DECONTAM_MIN_BLANKS
    assert f"0 controlli negativi, meno dei {minimo}" in testi
    assert "ctrl.blank_values e' vuoto" in testi
    assert "non ha controlli positivi (ctrl.positive_values e' vuoto)" in testi
    radice = Path(senza.config.io.out_root)
    manifesto = json.loads(next(radice.rglob("manifest_S0.json")).read_text())
    assert [d["codice"] for d in manifesto["degradazioni"]] == ["E-S0-17", "E-S0-17"]
    gates = json.loads((radice / Fase.INPUT_VALIDATION.value / "gates.json").read_text())
    assert len(next(g for g in gates["gate"] if g["gate"] == "G11")["avvisi"]) == 2

    # Etichette dichiarate, ma nessun campione le porta.
    non_riconosciuti = crea_scenario(tmp_path / "b", solo_biologici, con_letture=True)
    avvisi = esegui_gate("G11", Contesto(non_riconosciuti.config)).avvisi
    assert "nessun campione porta le etichette di ctrl.positive_values" in avvisi[1].dettaglio

    negativi = [Campione(f"ERX31000{i:02d}", f"B{i}", materiale=NEGATIVO) for i in range(5)]
    completo = crea_scenario(tmp_path / "c", [*_campioni(), *negativi], con_letture=True)
    assert not esegui_gate("G11", Contesto(completo.config)).avvisi


def test_i_thread_valgono_i_processori_disponibili_e_l_immagine_e_facoltativa(tmp_path, monkeypatch):
    """
    **Obiettivo**: Verificare che senza ``run.threads`` il parametro resti
    nullo (automatico) e non risulti dichiarato, che i thread effettivi siano i
    processori utilizzabili dal processo, che G14 lo accetti, e che senza
    ``run.container`` la configurazione sia valida, S0 passi e la provenienza
    registri un'immagine nulla.

    **Razionale scientifico e sistemistico**: Un numero fisso di thread fermava
    G14 su ogni macchina piu' piccola di quella di sviluppo senza cambiare
    alcun risultato; l'immagine e' una dichiarazione, e chi esegue fuori da un
    container non ne ha una da dichiarare.
    """
    scenario = crea_scenario(tmp_path, _campioni(), con_letture=True)
    dati = scenario.config.model_dump(mode="python")
    del dati["run"]
    monkeypatch.setattr(os, "sched_getaffinity", lambda pid: {0, 1, 2})
    config = valida(dati)
    assert processori_disponibili() == 3 and config.run.threads is None
    assert thread_effettivi(config) == 3
    assert config.run.container is None
    assert not {"run.threads", "run.container"} & set(parametri_dichiarati(config))
    assert esegui_gate("G14", Contesto(config)).superato
    assert esegui_s0(config).superata
    manifesto = json.loads(next(Path(config.io.out_root).rglob("manifest_S0.json")).read_text())
    assert manifesto["provenienza"]["immagine"] is None
    assert defaults.RUN_SEED == 100 and not hasattr(defaults, "RUN_THREADS")


# --------------------------------------------------------------------------- #
# 8. Report: i parametri tarati e non dichiarati                               #
# --------------------------------------------------------------------------- #


def test_il_report_elenca_solo_i_parametri_tarati_presi_per_difetto(tmp_path, capsys):
    """
    **Obiettivo**: Verificare che con una configurazione che dichiara solo i
    parametri obbligatori il report elenchi in apertura gli otto parametri con
    un predefinito tarato su OSD-734, ciascuno con il suo fatto, e nessun
    parametro obbligatorio; che dichiarandone uno al suo stesso valore e
    riprendendo l'esecuzione ne restino sette; che dichiarandoli tutti la
    sezione dica in una riga che non ce ne sono.

    **Razionale scientifico e sistemistico**: La segnalazione serve a far
    riesaminare le scelte di metodo ereditate: deve contenere solo quelle, e
    sparire quando chi esegue le ha fatte proprie dichiarandole, altrimenti
    resterebbe sempre accesa e nessuno la leggerebbe.
    """
    scenario = crea_scenario(tmp_path, _campioni(), con_letture=True)
    percorso = tmp_path / "c.yaml"
    dati = dichiarazione_minima(scenario.config)
    radice = Path(scenario.config.io.out_root)

    def report(dichiarazione: dict[str, Any]) -> tuple[str, str]:
        percorso.write_text(yaml.safe_dump(dichiarazione, sort_keys=False), encoding="utf-8")
        assert cli.main(["validate", "--config", str(percorso)]) == 0
        capsys.readouterr()
        assert cli.main(["report", "--config", str(percorso)]) == 0
        return capsys.readouterr().out, (radice / CARTELLA_REPORT / NOME_REPORT).read_text()

    uscita, documento = report(dati)
    assert len(defaults.FATTI_OSD734) == 8
    assert "8 parametri non dichiarati valgono il predefinito tarato su OSD-734" in uscita
    for chiave, riferimento in defaults.FATTI_OSD734.items():
        assert f">{chiave}<" in documento, chiave
    assert defaults.FATTI_OSD734["decontam.threshold"].fatto.split(":")[0] in documento
    tabella = documento[documento.index("Parametri non dichiarati che valgono"):]
    tabella = tabella[: tabella.index("</table>")]
    for chiave in defaults.OBBLIGATORI:
        assert f">{chiave}<" not in tabella, chiave

    dati.setdefault("prev", {})["min_fraction"] = defaults.PREV_MIN_FRACTION
    uscita, documento = report(dati)
    assert "7 parametri non dichiarati" in uscita

    for chiave, riferimento in defaults.FATTI_OSD734.items():
        gruppo, nome = chiave.split(".")
        dati.setdefault(gruppo, {})[nome] = riferimento.valore
    uscita, documento = report(dati)
    assert "tarato su OSD-734" not in uscita
    assert "Nessun parametro vale un predefinito tarato su OSD-734" in documento
    assert "Parametri non dichiarati che valgono" not in documento


# --------------------------------------------------------------------------- #
# 9. Le configurazioni dei due dataset                                         #
# --------------------------------------------------------------------------- #


@pytest.mark.skipif(not (RADICE / "dati" / "osd734").is_dir(),
                    reason="la cartella dati/ non e' presente (nell'immagine non viene copiata)")
def test_le_configurazioni_dei_due_dataset_dichiarano_tutto():
    """
    **Obiettivo**: Verificare che ``dati/osd734/config_osd734.yaml`` dichiari
    tutti i parametri obbligatori e tutti quelli tarati su OSD-734, con i
    valori di ``ESEMPIO_OSD734`` e di ``FATTI_OSD734``; che
    ``dati/osd276/config_osd276.yaml`` sia valida, dichiari tutti gli
    obbligatori e descriva un dataset senza controlli e senza file del lotto;
    e che nessun valore proprio del secondo dataset compaia nel codice.

    **Razionale scientifico e sistemistico**: Con tutto dichiarato il report
    di OSD-734 non ha valori ereditati da segnalare, e i suoi risultati non
    dipendono piu' dai predefiniti; il secondo dataset deve stare tutto nella
    sua configurazione, perche' il codice resti generale.
    """
    osd734 = carica(RADICE / "dati" / "osd734" / "config_osd734.yaml")
    dichiarati = set(parametri_dichiarati(osd734))
    assert set(defaults.OBBLIGATORI) | set(defaults.FATTI_OSD734) <= dichiarati
    for chiave, riferimento in {**defaults.ESEMPIO_OSD734, **defaults.FATTI_OSD734}.items():
        gruppo, nome = chiave.split(".")
        atteso = list(riferimento.valore) if isinstance(riferimento.valore, tuple) else riferimento.valore
        assert getattr(getattr(osd734, gruppo), nome) == atteso, chiave
    esegui_g15(yaml.safe_load((RADICE / "dati/osd734/config_osd734.yaml").read_text(encoding="utf-8")))

    percorso = RADICE / "dati" / "osd276" / "config_osd276.yaml"
    grezzo = yaml.safe_load(percorso.read_text(encoding="utf-8"))
    assert obbligatori_mancanti(grezzo) == []
    osd276 = esegui_g15(grezzo).config
    assert osd276.io.batch_table is None and osd276.out.batch_columns == []
    assert osd276.ctrl.positive_values == [] and osd276.ctrl.blank_values == []
    assert osd276.katharoseq.target_taxon is None and osd276.err.batch_column is None

    propri = ("SRR5336", "OSD-276", "osd276", "GLDS-276", "PRJNA376404")
    for sorgente in [*(RADICE / "src").rglob("*.py"), *(RADICE / "R").rglob("*.R")]:
        testo = sorgente.read_text(encoding="utf-8")
        assert not [p for p in propri if p in testo], sorgente


def test_s0_dichiara_i_parametri_che_i_suoi_gate_leggono(tmp_path):
    """
    **Obiettivo**: Verificare che S0 dichiari fra i propri parametri quelli che
    i gate leggono dopo questa settimana (``filter.trimLeft``,
    ``out.study_columns``, ``out.batch_columns``,
    ``katharoseq.cell_count_column``), che la sua versione sia almeno 3, e che
    cambiare ``filter.trimLeft`` cambi l'impronta su cui S0 e' calcolata.

    **Razionale scientifico e sistemistico**: Un parametro letto e non
    dichiarato non entra nell'impronta della fase: cambiarlo lascerebbe S0
    conclusa con gli esiti dei gate calcolati sul valore vecchio.
    """
    fase = ValidazioneIngressi()
    assert fase.versione >= 3 and fase.passo is Passo.S0
    assert {"filter.trimLeft", "out.study_columns", "out.batch_columns",
            "katharoseq.cell_count_column"} <= set(fase.parametri)
    scenario = crea_scenario(tmp_path, _campioni(), con_letture=True)
    prima = fase.calcolata_su(risolvi(scenario.config), {})
    dopo = fase.calcolata_su(risolvi(_con(scenario.config, filter={"trimLeft": 3})), {})
    assert prima["configurazione"]["impronta"] != dopo["configurazione"]["impronta"]
    assert copy.deepcopy(prima) == fase.calcolata_su(risolvi(scenario.config), {})


# --------------------------------------------------------------------------- #
# 9. Correzioni dopo la revisione esterna                                      #
# --------------------------------------------------------------------------- #


def test_il_digest_non_dipende_dai_processori_della_macchina(tmp_path, monkeypatch):
    """
    **Obiettivo**: Verificare che la stessa configurazione, senza
    ``run.threads``, abbia lo stesso digest e la stessa impronta dei risultati
    su una macchina con 4 processori e su una con 16 (simulati), che
    ``resolved.yaml`` riporti ``threads`` nullo, e che il valore effettivo si
    calcoli all'uso: 4 e 16. Un valore dichiarato resta quello dichiarato,
    entra nel digest e G14 lo confronta con i processori utilizzabili.

    **Razionale scientifico e sistemistico**: Il digest identifica la
    configurazione: se dipendesse dalla macchina, due esecuzioni della stessa
    configurazione congelata risulterebbero diverse, e ogni ripresa su un'altra
    macchina registrerebbe una configurazione nuova.
    """
    scenario = crea_scenario(tmp_path, _campioni(), con_letture=True)
    dati = scenario.config.model_dump(mode="python")
    del dati["run"]
    visti = {}
    for processori in (4, 16):
        monkeypatch.setattr(os, "sched_getaffinity", lambda pid, n=processori: set(range(n)))
        config = valida(copy.deepcopy(dati))
        risolta = risolvi(config)
        assert config.run.threads is None and thread_effettivi(config) == processori
        assert risolta.come_mappa()["run"]["threads"] is None
        assert esegui_gate("G14", Contesto(config)).superato
        registrata = yaml.safe_load(
            scrivi_risolta(risolta, tmp_path / f"out_{processori}").read_text(encoding="utf-8"))
        assert registrata["parametri"]["run"]["threads"] is None
        assert registrata["digest"] == risolta.digest
        visti[processori] = (risolta.digest, risolta.impronta_risultati)
    assert visti[4] == visti[16]

    dichiarata = _con(config, run={"threads": 2})
    assert thread_effettivi(dichiarata) == 2
    assert risolvi(dichiarata).digest != visti[16][0]
    assert risolvi(dichiarata).impronta_risultati == visti[16][1]
    troppi = esegui_gate("G14", Contesto(_con(config, run={"threads": 17})))
    assert [v.codice for v in troppi.violazioni] == ["E-S0-14"]


def test_senza_colonna_della_posizione_il_modulo_viene_dal_file_del_lotto(tmp_path):
    """
    **Obiettivo**: Verificare che con ``meta.module_column`` nullo e
    ``meta.batch_module_column`` dichiarata ogni campione prenda il modulo dal
    file del lotto e resti senza posizione; che con la colonna della posizione
    dichiarata la regola delle non superfici continui a valere; e che G15
    respinga con ``E-G15-14`` ``meta.module_regex`` e
    ``meta.non_surface_positions`` valorizzati con ``meta.module_column``
    nullo.

    **Razionale scientifico e sistemistico**: Senza una posizione dichiarata la
    regola delle non superfici non ha nulla da giudicare, e scartare il modulo
    che il file del lotto dichiara lascerebbe i campioni senza raggruppamento
    senza alcun messaggio; un'espressione senza la colonna a cui applicarla non
    deriverebbe mai nulla.
    """
    campioni = _campioni()
    senza = crea_scenario(
        tmp_path / "a", campioni, con_arricchimento=True, con_letture=True,
        sovrascrivi={"meta": {"module_column": None, "module_regex": None,
                              "non_surface_positions": []},
                     "ctrl": {"blank_override_column": None, "blank_override_values": []}},
    )
    assert esegui_gate("G15", Contesto(senza.config)).superato
    inventario = Contesto(senza.config).inventario
    assert [c.modulo for c in inventario] == ["Modulo Uno"] * len(campioni)
    assert [c.posizione for c in inventario] == [None] * len(campioni)

    con = crea_scenario(tmp_path / "b", campioni, con_arricchimento=True, con_letture=True)
    moduli = {c.nome: c.modulo for c in Contesto(con.config).inventario}
    assert moduli["NOD1D4.L1"] == "Modulo Uno" and moduli["POS.P1.1"] is None

    incoerente = crea_scenario(
        tmp_path / "c", campioni, con_arricchimento=True,
        sovrascrivi={"meta": {"module_column": None},
                     "ctrl": {"blank_override_column": None, "blank_override_values": []}},
    )
    # L'espressione e le posizioni non di superficie, entrambe senza la colonna:
    # due violazioni, dette insieme.
    esito = esegui_gate("G15", Contesto(incoerente.config))
    assert [v.codice for v in esito.violazioni] == ["E-G15-14", "E-G15-14"]
    testi = " | ".join(v.dettaglio for v in esito.violazioni)
    assert "meta.module_regex" in testi and "meta.non_surface_positions" in testi
    assert CATALOGO["E-G15-14"].categoria is Categoria.REVISIONE_UMANA


def test_g10_non_cerca_il_primer_nei_controlli_negativi(tmp_path):
    """
    **Obiettivo**: Verificare che un controllo negativo con tre letture, di cui
    una comincia con il primer, non fermi G10; e che lo stesso primer in testa
    alle letture di un campione biologico lo fermi con ``E-S0-10``, nominando
    quel solo file.

    **Razionale scientifico e sistemistico**: In un bianco con una manciata di
    letture una frazione non ha significato: una lettura su tre supera
    qualunque soglia senza dire nulla di come le letture sono state prodotte.
    Il controllo del primer riguarda le stesse classi del controllo del motivo.
    """
    scenario = crea_scenario(tmp_path, _campioni(), con_letture=True)
    bianco = next(Path(scenario.config.io.fastq_dir).glob("*ERX3000004*"))
    scrivi_fastq(bianco, [lettura(INIZIO_CON_PRIMER)] + [lettura()] * 2)
    contesto = Contesto(scenario.config)
    assert contesto.scansione["ERX3000004"].frazione_primer > scenario.config.qc.max_primer_hit_frac
    assert esegui_gate("G10", contesto).superato

    biologico = next(Path(scenario.config.io.fastq_dir).glob("*ERX3000001*"))
    scrivi_fastq(biologico, [lettura(INIZIO_CON_PRIMER)] * 40)
    esito = esegui_gate("G10", Contesto(scenario.config))
    assert [v.codice for v in esito.violazioni] == ["E-S0-10"]
    assert biologico.name in esito.violazioni[0].dettaglio
    assert bianco.name not in esito.violazioni[0].dettaglio
    assert esito.violazioni[0].dettaglio.startswith("1 file")


@pytest.mark.parametrize("gruppo, attesi", [
    ("tax", {"tax.ref_fasta", "tax.ref_md5"}),
    ("io", {"io.fastq_dir", "io.assay_table", "io.out_root"}),
])
def test_un_gruppo_mancante_elenca_tutti_i_suoi_parametri_senza_predefinito(tmp_path, gruppo, attesi):
    """
    **Obiettivo**: Verificare che, tolto dalla configurazione un intero gruppo,
    l'elenco dei problemi nomini ogni suo parametro senza predefinito: quelli
    di ``OBBLIGATORI`` nella voce che li riunisce, gli altri (percorsi,
    checksum del riferimento) ciascuno con la propria riga; e che nessuno sia
    detto due volte.

    **Razionale scientifico e sistemistico**: Chi compila la configurazione di
    un dataset nuovo deve vedere l'elenco intero in una sola esecuzione: un
    parametro taciuto perche' il suo gruppo manca verrebbe scoperto solo dopo
    aver corretto gli altri.
    """
    scenario = crea_scenario(tmp_path, _campioni())
    dati = scenario.config.model_dump(mode="json")
    del dati[gruppo]
    with pytest.raises(ErroreConfigurazione) as info:
        valida(dati)
    problemi = info.value.problemi
    obbligatori = {c for c in defaults.OBBLIGATORI if c.startswith(f"{gruppo}.")}
    assert problemi[0].startswith(INTESTAZIONE_MANCANTI)
    assert all(c in problemi[0] for c in obbligatori)
    singoli = {p.split(":")[0] for p in problemi[1:]}
    assert singoli == attesi and not singoli & obbligatori
    assert all("parametro obbligatorio mancante" in p for p in problemi[1:])
    assert len(problemi) == 1 + len(attesi)
    with pytest.raises(ErroreGate) as rifiuto:
        esegui_g15(dati)
    assert [v.codice for v in rifiuto.value.violazioni] == ["E-G15-10"] + ["E-G15-99"] * len(attesi)


def test_g15_respinge_un_parametro_dichiarato_senza_quello_da_cui_dipende(tmp_path):
    """
    **Obiettivo**: Verificare che G15 respinga con ``E-G15-14``
    ``ctrl.blank_override_values`` non vuoto con ``ctrl.blank_override_column``
    nullo, e ``meta.study_sample_id_column`` dichiarato con ``io.study_table``
    nullo; che
    ``run.threads: true`` sia respinto; e che le stesse configurazioni,
    corrette, passino.

    **Razionale scientifico e sistemistico**: Un parametro che ha effetto solo
    insieme a un altro, dichiarato da solo, verrebbe ignorato in silenzio: chi
    lo ha scritto crederebbe di avere riclassificato dei campioni, o di leggere
    i nomi da un'altra colonna.
    """
    scenario = crea_scenario(tmp_path, _campioni())
    dati = scenario.config.model_dump(mode="json")
    assert esegui_g15(dati).config.ctrl.blank_override_values

    senza_colonna = copy.deepcopy(dati)
    senza_colonna["ctrl"]["blank_override_column"] = None
    with pytest.raises(ErroreGate) as info:
        esegui_g15(senza_colonna)
    assert [v.codice for v in info.value.violazioni] == ["E-G15-14"]
    assert "blank_override_values" in info.value.violazioni[0].dettaglio

    senza_studio = copy.deepcopy(dati)
    senza_studio["io"]["study_table"] = None
    senza_studio["meta"]["study_sample_id_column"] = "Source Name"
    with pytest.raises(ErroreGate) as info:
        esegui_g15(senza_studio)
    assert [v.codice for v in info.value.violazioni] == ["E-G15-14"]
    assert "meta.study_sample_id_column" in info.value.violazioni[0].dettaglio
    assert "io.study_table" in info.value.violazioni[0].dettaglio

    senza_studio["meta"]["study_sample_id_column"] = None
    assert esegui_g15(senza_studio).config.io.study_table is None

    booleano = copy.deepcopy(dati)
    booleano["run"]["threads"] = True
    with pytest.raises(ErroreGate) as info:
        esegui_g15(booleano)
    assert "run.threads" in str(info.value)


def test_g08_dichiara_le_righe_del_lotto_senza_campione(tmp_path):
    """
    **Obiettivo**: Verificare che, con due righe del file del lotto che non
    corrispondono ad alcun campione, G08 passi con un avviso ``E-S0-19`` che ne
    riporta il numero e le chiavi, che S0 lo registri fra le degradazioni del
    manifesto, e che senza righe in piu' l'avviso non ci sia.

    **Razionale scientifico e sistemistico**: Una riga senza campione puo'
    essere legittima (un file che copre piu' assay), ma e' anche il sintomo di
    una chiave scritta in modo diverso da quella dei campioni: ignorarla in
    silenzio nasconderebbe il secondo caso.
    """
    scenario = crea_scenario(tmp_path, _campioni(), con_arricchimento=True, con_letture=True)
    assert not [a for a in esegui_gate("G08", Contesto(scenario.config)).avvisi
                if a.codice == "E-S0-19"]

    lotto = Path(scenario.config.io.batch_table)
    with open(lotto, "a", encoding="utf-8") as file:
        file.write("ERX3999998\t1\tcorsa_A\tModulo Uno\n")
        file.write("ERX3999999\t1\tcorsa_A\tModulo Uno\n")
    esito = esegui_gate("G08", Contesto(scenario.config))
    assert esito.superato
    (avviso,) = [a for a in esito.avvisi if a.codice == "E-S0-19"]
    assert avviso.dettaglio.startswith("2 righe")
    assert "ERX3999998" in avviso.dettaglio and "ERX3999999" in avviso.dettaglio
    assert CATALOGO["E-S0-19"].categoria is Categoria.DEGRADAZIONE_AUTOMATICA

    risultato = esegui_s0(scenario.config)
    assert risultato.superata
    manifesto = json.loads(
        next(Path(scenario.config.io.out_root).rglob("manifest_S0.json")).read_text())
    assert "E-S0-19" in [d["codice"] for d in manifesto["degradazioni"]]


def _modulo_delle_letture_del_secondo_dataset():
    """Lo script ``dati/osd276/scarica_letture.py``, caricato con un nome suo
    (quello di ``dati/osd734`` ha lo stesso nome di file).
    """
    import importlib.util

    percorso = RADICE / "dati" / "osd276" / "scarica_letture.py"
    specifica = importlib.util.spec_from_file_location("scarica_letture_secondo", percorso)
    modulo = importlib.util.module_from_spec(specifica)
    specifica.loader.exec_module(modulo)
    return modulo


def _deposito(percorso: Path, intestazioni: list[str]) -> None:
    """Un file come quelli del deposito: un record per intestazione."""
    percorso.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(percorso, "wb") as file:
        for intestazione in intestazioni:
            file.write(f"@{intestazione}\nACGT\n+\nIIII\n".encode())


@pytest.mark.skipif(not (RADICE / "dati" / "osd276").is_dir(),
                    reason="la cartella dati/ non e' presente (nell'immagine non viene copiata)")
def test_le_forward_si_ricavano_solo_da_un_deposito_con_le_coppie_complete(tmp_path, monkeypatch, capsys):
    """
    **Obiettivo**: Verificare che lo script che ricava le letture forward del
    secondo dataset conti i record ``/1``, ``/2`` e gli altri, e rifiuti un
    deposito con record che non sono di una coppia, con ``/1`` e ``/2`` in
    numero diverso, o con una somma diversa dalle letture dichiarate; che un
    file forward con un'impronta diversa dall'attesa venga messo da parte e non
    resti in ``fastq/``; e che il guasto di una corsa (un archivio illeggibile)
    finisca nel riepilogo senza fermare le altre, con esito 1.

    **Razionale scientifico e sistemistico**: Un deposito che non ha la forma
    attesa darebbe letture forward plausibili e sbagliate, e un file non
    conforme lasciato nella cartella delle letture verrebbe analizzato: i
    controlli devono fermare prima, e dire per quale corsa.
    """
    import hashlib

    script = _modulo_delle_letture_del_secondo_dataset()
    coppie = [f"c.{i} {i}/1" for i in range(1, 4)] + [f"c.{i} {i}/2" for i in range(1, 4)]
    buono = tmp_path / "prove" / "buono.fastq.gz"
    _deposito(buono, coppie)
    assert script.conta_letture(buono) == (3, 3, 0)
    assert script.verifica_conteggi(buono, 6) == 3
    for nome, intestazioni, dichiarate, atteso in (
        ("altri", coppie + ["c.9 senza marcatore"], 7, "ne' la prima ne' la seconda"),
        ("spaiate", coppie[:5], 5, "non sono in numero uguale"),
        ("somma", coppie, 8, "letture dichiarato dal deposito (8)"),
    ):
        percorso = tmp_path / "prove" / f"{nome}.fastq.gz"
        _deposito(percorso, intestazioni)
        with pytest.raises(script.DepositoInatteso, match=re.escape(atteso)):
            script.verifica_conteggi(percorso, dichiarate)

    # Tre corse gia' scaricate: una conforme, una le cui forward non hanno
    # l'impronta attesa, una con un archivio illeggibile.
    cartella = tmp_path / "dati"
    forward = "".join(f"@{i}\nACGT\n+\nIIII\n" for i in coppie[:3]).encode()
    righe, impronte = [], []
    for corsa, contenuto_md5 in (("CORSA1", hashlib.md5(forward).hexdigest()),
                                 ("CORSA2", "0" * 32), ("CORSA3", "0" * 32)):
        deposito = cartella / "ena" / f"{corsa}.fastq.gz"
        if corsa == "CORSA3":
            deposito.parent.mkdir(parents=True, exist_ok=True)
            deposito.write_bytes(b"non e' un archivio")
        else:
            _deposito(deposito, coppie)
        nome = f"{corsa}_campione_R1.fastq.gz"
        righe.append([nome, corsa, "ESP1", "campione", "https://esempio.invalid/" + corsa,
                      hashlib.md5(deposito.read_bytes()).hexdigest(),
                      str(deposito.stat().st_size), "6", "3"])
        impronte.append(f"{contenuto_md5}  {nome}\n")
    elenco = tmp_path / "letture_ena.tsv"
    elenco.write_text("\t".join(
        ["file", "run_accession", "experiment_accession", "sample_name", "fastq_url",
         "fastq_md5", "fastq_bytes", "letture_deposito", "letture_forward"]) + "\n"
        + "".join("\t".join(r) + "\n" for r in righe), encoding="utf-8")
    (tmp_path / "letture_forward.md5").write_text("".join(impronte), encoding="utf-8")
    monkeypatch.setattr(script, "ELENCO", elenco)
    monkeypatch.setattr(script, "IMPRONTE", tmp_path / "letture_forward.md5")

    assert script.main(["--cartella", str(cartella), "--pausa", "0"]) == 1
    uscita = capsys.readouterr().out
    assert "file forward conformi 1, mancanti o non conformi 2" in uscita
    presenti = sorted(p.name for p in (cartella / "fastq").iterdir())
    assert presenti == ["CORSA1_campione_R1.fastq.gz", "CORSA2_campione_R1.fastq.gz.md5_errato"]
    assert "CORSA2_campione_R1.fastq.gz: DepositoInatteso" in uscita and "messo da parte" in uscita
    assert "CORSA3_campione_R1.fastq.gz: BadGzipFile" in uscita
    assert script.main(["--cartella", str(cartella), "--solo-verifica"]) == 1

    # Un file forward gia' presente ma illeggibile e' non conforme, non un arresto.
    (cartella / "fastq" / "CORSA1_campione_R1.fastq.gz").write_bytes(b"\x1f\x8b troncato")
    capsys.readouterr()
    assert script.main(["--cartella", str(cartella), "--solo-verifica"]) == 1
    uscita = capsys.readouterr().out
    assert "CORSA1_campione_R1.fastq.gz: illeggibile" in uscita
    assert "file forward conformi 0, mancanti o non conformi 3" in uscita
