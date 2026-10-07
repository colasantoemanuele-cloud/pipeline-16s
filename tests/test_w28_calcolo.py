r"""Suite di test della settimana 28: generalita' delle fasi di calcolo e del report.

1. Inquadramento nel Piano Operativo
------------------------------------
Settimana 28 (W28), Fase F8 (un dataset 16S single-end diverso da OSD-734
attraversa le fasi di calcolo S1-S14 e il report, o vi si ferma con un codice
del catalogo e mai con un errore generico di R).

2. Moduli sorgente coperti
--------------------------
* ``src/amplicon16s/gates/g01_g15.py`` (G09, G02, G08, G15),
  ``src/amplicon16s/io_layer/reads.py``, ``src/amplicon16s/config/schema.py``,
  ``src/amplicon16s/config/resolve.py``, ``src/amplicon16s/metadata/tabelle.py``
* ``src/amplicon16s/steps/s01_profile.py``, ``R/01_profile.R``
* ``src/amplicon16s/steps/s03_learn_errors.py``, ``R/03_learn_errors.R``,
  ``R/04_dada.R``, ``R/lib/errore_loess.R``
* ``src/amplicon16s/steps/s10_phyloseq.py`` (colonna dei livelli)
* ``src/amplicon16s/steps/s11_controls.py``, ``R/11_controls.R``
* ``src/amplicon16s/steps/s12_decontam.py``, ``R/12_decontam.R``
* ``src/amplicon16s/steps/s13_filtri.py``, ``R/13_filtri.R``
* ``src/amplicon16s/steps/s14_finale.py``, ``R/14_finale.R``
* ``src/amplicon16s/report/builder.py`` (troncamento suggerito, classi vuote)
* ``src/amplicon16s/errors/catalog.py`` (``E-S1-03``, ``E-S3-03``,
  ``E-S3-04``, ``E-S11-05``, ``E-S12-03``, ``E-S13-04``, ``E-S13-05``)

3. Cosa valuta questo file
--------------------------
- il troncamento si giudica sulla frazione di letture piu' corte dei campioni
  biologici e dei controlli positivi (``qc.max_frac_short_reads``): poche
  letture corte non fermano ne' G09 ne' S1, una quota rilevante si', e i
  controlli negativi non contano; ``filter.minLen`` non esiste piu';
- il modello di errore: ``err.error_function`` accetta ``loess`` e
  ``loess_monotono``, la variante monotona da' tassi non crescenti con la
  qualita'; S1 avvisa con ``E-S1-03`` se le qualita' distinte sono poche; una
  corsa senza letture filtrate ferma S3 con ``E-S3-03``, e letture con un solo
  valore di qualita' con ``E-S3-04``;
- un dataset senza controlli: S11 senza positivi dichiara ``E-S11-05`` e non
  adatta curve, non si ferma per un controllo senza piastra o con una sola
  lettura; S12 con pochi negativi dichiara ``E-S12-03`` e non toglie nulla;
- la colonna dei livelli dei controlli positivi si cerca nel file del lotto e
  nella tabella di studio, ed entra nell'oggetto da sola;
- i filtri finali: i nomi dei taxa si confrontano senza il prefisso di rango;
  una tassonomia senza Phylum e' dichiarata con ``E-S13-05``; se nessun
  campione supera i filtri la fase si ferma con ``E-S13-04``;
  ``tax.assign_species`` vero e' respinto dallo schema;
- ``12_final``: i file consegnati di un'esecuzione precedente si tolgono
  prima del calcolo, e senza controlli ``ps_controlli.rds`` non esiste;
- il report: troncamento suggerito, e classi di controllo vuote;
- la catena intera su quattro insiemi ridotti (senza positivi, senza negativi,
  senza file del lotto, con una sola piastra) si conclude, o si ferma con un
  codice del catalogo che non e' del ponte verso R, e il report si genera.

4. Comandi Bash e scenari di esecuzione
---------------------------------------
    ``<immagine>`` e' l'immagine del container della pipeline; quella corrente
    e' indicata in ``test.txt``, sezione 1.3.

    1. Modalità locale standard (R di base con jsonlite, senza Bioconductor né
       dati reali):
       pytest tests/test_w28_calcolo.py -v

    2. Modalità container Docker standard (sottoinsieme ridotto con
       R/Bioconductor):
       docker run --rm \
         -e PYTHONPATH=/app/src \
         -e AMPLICON16S_R_DIR=/app/R \
         -v "$(pwd)":/app \
         -w /app \
         <immagine> \
         pytest -o cache_dir=/tmp/.pytest_cache tests/test_w28_calcolo.py -v

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
         pytest -o cache_dir=/tmp/.pytest_cache tests/test_w28_calcolo.py -v

5. Risultato atteso
-------------------
Vedi ``test.txt``, scheda W28.

6. Razionale scientifico e sistemistico
---------------------------------------
- Una fase che presuppone la struttura di OSD-734 (controlli positivi e
  negativi, piastre, letture tutte piu' lunghe del troncamento) si ferma su un
  altro dataset con un errore di R che non dice che cosa manca: ogni assenza
  dev'essere una condizione dichiarata, con il suo codice.
- Un ripiego silenzioso (nessuna decontaminazione, nessuna curva) produce un
  risultato plausibile e non confrontabile: resta nel manifesto e nel report.
"""

from __future__ import annotations

import csv
import json
import os
import subprocess
from pathlib import Path
from typing import Any

import pytest
from conftest import NEGATIVO, POSITIVO, Campione, copia_esecuzione, crea_scenario, lettura, scrivi_fastq
from sottoinsieme import RIDOTTO, dati_config, motivo_pacchetti_r_assenti, selezione

from amplicon16s.config.resolve import risolvi
from amplicon16s.config.schema import Config, ErroreConfigurazione, valida
from amplicon16s.errors.catalog import CATALOGO
from amplicon16s.errors.exceptions import ErrorePipeline
from amplicon16s.gates.g01_g15 import Contesto
from amplicon16s.gates.registry import esegui_gate
from amplicon16s.io_layer.artifacts import Fase
from amplicon16s.logging.logger import chiudi
from amplicon16s.metadata.models import ClasseCampione
from amplicon16s.rbridge.runner import cartella_r, trova_rscript
from amplicon16s.report.builder import _pct, genera, troncamento_suggerito
from amplicon16s.runner.executor import Conclusione, Esecutore
from amplicon16s.runner.graph import Passo
from amplicon16s.runner.project import ProjectRun
from amplicon16s.steps.s01_profile import letture_corte
from amplicon16s.steps.s10_phyloseq import NOME_OGGETTO, colonne_metadati
from amplicon16s.steps.s14_finale import NOME_CONTROLLI


@pytest.fixture(autouse=True)
def uscite_pulite():
    """Chiude le uscite del log prima e dopo ogni test, perche' nessun handler
    resti aperto sulla cartella temporanea.
    """
    chiudi()
    yield
    chiudi()


_MOTIVO_R = motivo_pacchetti_r_assenti("jsonlite")
_MOTIVO_BIOC = motivo_pacchetti_r_assenti(
    "dada2", "ggplot2", "ShortRead", "jsonlite", "phyloseq", "Biostrings", "decontam"
)


@pytest.fixture
def r():
    """Richiede R con jsonlite: salta senza, ma in CI fallisce."""
    if _MOTIVO_R is not None:
        if os.environ.get("AMPLICON16S_RICHIEDI_R") == "1":
            pytest.fail(f"R e' richiesto in questo ambiente: {_MOTIVO_R}")
        pytest.skip(_MOTIVO_R)


@pytest.fixture
def bioc():
    """Richiede R con i pacchetti di tutte le fasi: salta senza, ma in CI fallisce."""
    if _MOTIVO_BIOC is not None:
        if os.environ.get("AMPLICON16S_RICHIEDI_BIOC") == "1":
            pytest.fail(f"Bioconductor e' richiesto in questo ambiente: {_MOTIVO_BIOC}")
        pytest.skip(_MOTIVO_BIOC)


def _tsv(percorso: Path) -> list[dict[str, str]]:
    """Le righe di una tabella separata da tabulazioni."""
    with open(percorso, encoding="utf-8", newline="") as file:
        return list(csv.DictReader(file, delimiter="\t"))


def _r(codice: str, cartella: Path) -> Any:
    """Esegue codice R che scrive un oggetto JSON in ``uscita``."""
    cartella.mkdir(parents=True, exist_ok=True)
    uscita = cartella / "uscita.json"
    script = cartella / "codice.R"
    script.write_text(f"uscita <- {json.dumps(str(uscita))}\n{codice}\n", encoding="utf-8")
    subprocess.run([str(trova_rscript()), "--vanilla", str(script)], check=True,
                   capture_output=True)
    return json.loads(uscita.read_text(encoding="utf-8"))


def _con(config: Config, **gruppi: dict[str, Any]) -> Config:
    """La configurazione con i parametri indicati sostituiti, rivalidata."""
    dati = config.model_dump(mode="python")
    for gruppo, valori in gruppi.items():
        dati[gruppo].update(valori)
    return valida(dati)


def _campioni() -> list[Campione]:
    """Due biologici, un positivo e un negativo, con la forma di OSD-734."""
    return [
        Campione("ERX3000001", "NOD1D4.L1"),
        Campione("ERX3000002", "NOD1D4.L2"),
        Campione("ERX3000003", "POS.P1.1", materiale=POSITIVO, posizione="Not Applicable"),
        Campione("ERX3000004", "BLANK.P1.1", materiale=NEGATIVO, posizione="Not Applicable"),
    ]


def _file_di(scenario, accession: str) -> Path:
    """Il file FASTQ di un campione dello scenario."""
    (percorso,) = Path(scenario.config.io.fastq_dir).glob(f"*{accession}*")
    return percorso


def _degradazioni(run, passo: Passo, fase: Fase) -> list[str]:
    """I codici delle degradazioni registrate nel manifesto di una fase."""
    return [d["codice"] for d in run.albero.manifesto_passo(passo, fase).degradazioni]


# --------------------------------------------------------------------------- #
# 1. Troncamento: la frazione di letture piu' corte (criterio H)               #
# --------------------------------------------------------------------------- #


def test_poche_letture_corte_non_fermano_g09(tmp_path):
    """
    **Obiettivo**: Verificare che, con letture a lunghezza variabile in cui il
    2,5% dei campioni biologici e' piu' corto di ``filter.truncLen``, G09
    passi; che con ``qc.max_frac_short_reads`` a 0,01 la stessa quota lo
    fermi con ``E-S0-09``, nominando frazione e parametro.

    **Razionale scientifico e sistemistico**: In un dataset a lunghezza
    variabile qualche lettura corta c'e' sempre, e il filtro la toglie senza
    danno: fermarsi per una sola lettura rende la pipeline inutilizzabile
    fuori da OSD-734, non fermarsi mai azzera campioni in silenzio.
    """
    scenario = crea_scenario(tmp_path, _campioni(), con_letture=True)
    for campione in scenario.campioni[:2]:
        scrivi_fastq(_file_di(scenario, campione.accession),
                     [lettura(lunghezza=120)] + [lettura()] * 39)
    assert esegui_gate("G09", Contesto(scenario.config)).superato

    severa = _con(scenario.config, qc={"max_frac_short_reads": 0.01})
    esito = esegui_gate("G09", Contesto(severa))
    assert [v.codice for v in esito.violazioni] == ["E-S0-09"]
    assert "qc.max_frac_short_reads" in esito.violazioni[0].dettaglio
    assert "2 su 120" in esito.violazioni[0].dettaglio


def test_le_letture_corte_dei_controlli_negativi_non_contano(tmp_path):
    """
    **Obiettivo**: Verificare che un controllo negativo fatto di sole letture
    piu' corte del troncamento non fermi G09, e che ``letture_corte`` riporti
    la sua frazione a parte senza includerla in quella dei campioni
    controllati.

    **Razionale scientifico e sistemistico**: Un bianco amplifica poco e male,
    spesso solo dimeri corti: e' l'esito atteso, e non dice nulla sul
    troncamento adatto ai campioni.
    """
    scenario = crea_scenario(tmp_path, _campioni(), con_letture=True)
    scrivi_fastq(_file_di(scenario, "ERX3000004"), [lettura(lunghezza=60)] * 40)
    assert esegui_gate("G09", Contesto(scenario.config)).superato

    corte = letture_corte(
        {"a": {151: 95, 120: 5}, "b": {60: 40}, "c": {151: 100}},
        {"a": ClasseCampione.BIOLOGICO, "b": ClasseCampione.CONTROLLO_NEGATIVO,
         "c": ClasseCampione.CONTROLLO_POSITIVO},
        137,
    )
    assert corte["controllate"] == {"letture": 200, "corte": 5, "frazione": 0.025}
    assert corte["controllo_negativo"]["frazione"] == 1.0
    assert corte["biologico"]["corte"] == 5 and corte["controllo_positivo"]["corte"] == 0


def test_filter_minlen_non_esiste_e_la_specie_e_respinta(tmp_path):
    """
    **Obiettivo**: Verificare che ``filter.minLen`` non sia piu' fra i
    derivati e sia respinto come chiave sconosciuta; che ``E-G15-01`` non sia
    nel catalogo; che ``tax.assign_species: true`` sia respinto dallo schema;
    e che ``err.error_function`` accetti solo ``loess`` e ``loess_monotono``.

    **Razionale scientifico e sistemistico**: Un derivato che coincide per
    costruzione con un altro parametro e il controllo che lo sorveglia non
    verificano nulla; un parametro accettato e non realizzato promette un
    risultato che la pipeline non da'.
    """
    scenario = crea_scenario(tmp_path, _campioni())
    assert "filter.minLen" not in risolvi(scenario.config).derivati.come_chiavi()
    assert "E-G15-01" not in CATALOGO
    dati = scenario.config.model_dump(mode="python")
    for gruppo, chiave, valore in (("filter", "minLen", 137), ("tax", "assign_species", True),
                                   ("err", "error_function", "lineare")):
        modificati = json.loads(json.dumps(dati, default=str))
        modificati[gruppo][chiave] = valore
        with pytest.raises(ErroreConfigurazione):
            valida(modificati)
    assert _con(scenario.config, err={"error_function": "loess_monotono"}).err.error_function \
        == "loess_monotono"


def test_s1_non_si_ferma_per_poche_letture_corte(bioc, tmp_path):
    """
    **Obiettivo**: Verificare sul sottoinsieme ridotto, le cui letture vanno da
    137 a 151 bp, che con ``filter.truncLen`` a 151 (il 3,7% delle letture dei
    campioni controllati e' piu' corto) S1 si concluda e registri nel
    manifesto la frazione per classe; e che con ``qc.max_frac_short_reads`` a
    0,01 la stessa quota fermi l'esecuzione, gia' in S0 con ``E-S0-09`` perche'
    G09 qui legge tutte le letture; e che, quando G09 guarda la sola prima
    lettura di ogni file (``qc.head_reads`` 1, troncamento a 149), sia S1 a
    fermarsi con ``E-S1-02``, con la frazione per classe nel dettaglio.

    **Razionale scientifico e sistemistico**: S1 legge tutte le letture e
    chiude il limite di G09: la sua regola dev'essere la stessa, la frazione e
    non la singola lettura piu' corta.
    """
    def configura(cartella: Path, **qc: Any) -> Config:
        return valida(dati_config(cartella, filter={"truncLen": 151}, qc=qc))

    run = ProjectRun(configura(tmp_path / "a"))
    esito = Esecutore(run, fino_a=Passo.S1).esegui()
    assert esito.conclusione is Conclusione.COMPLETATA
    corte = run.albero.manifesto_passo(Passo.S1, Fase.QC_PROFILES).metriche[
        "letture_piu_corte_del_troncamento"]
    assert 0 < corte["controllate"]["frazione"] < 0.05
    assert corte["controllate"]["corte"] > 600

    severa = ProjectRun(configura(tmp_path / "b", max_frac_short_reads=0.01))
    esito = Esecutore(severa, fino_a=Passo.S1).esegui()
    assert esito.conclusione is Conclusione.ARRESTATA
    assert esito.punto.codice == "E-S0-09" and "qc.max_frac_short_reads" in esito.punto.dettaglio

    oltre_la_testa = ProjectRun(valida(dati_config(
        tmp_path / "c", filter={"truncLen": 149},
        qc={"max_frac_short_reads": 0.01, "head_reads": 1})))
    esito = Esecutore(oltre_la_testa, fino_a=Passo.S1).esegui()
    assert esito.conclusione is Conclusione.ARRESTATA
    assert [r.passo for r in esito.eseguite] == [Passo.S0]
    assert esito.punto.codice == "E-S1-02" and "Per classe: biologico" in esito.punto.dettaglio


# --------------------------------------------------------------------------- #
# 2. Modello di errore                                                         #
# --------------------------------------------------------------------------- #


def test_la_funzione_di_errore_monotona_non_cresce_con_la_qualita(r, tmp_path):
    """
    **Obiettivo**: Verificare che ``loess_monotono`` (``R/lib/errore_loess.R``),
    su conteggi di transizione in cui il tasso osservato risale alle qualita'
    alte, dia per ogni sostituzione un tasso non crescente con la qualita',
    entro i limiti dichiarati, e che le probabilita' di ogni base di partenza
    sommino a uno.

    **Razionale scientifico e sistemistico**: Con qualita' raggruppate in pochi
    valori la stima loess standard puo' dare un tasso di errore che cresce con
    la qualita', e dada tratterebbe come piu' affidabili le basi peggiori.
    """
    esito = _r(f"""
source({json.dumps(str(cartella_r() / "lib" / "errore_loess.R"))})
basi <- c("A", "C", "G", "T")
qualita <- c(2, 12, 23, 37)
nomi <- paste0(rep(basi, each = 4), "2", basi)
trans <- matrix(0, nrow = 16, ncol = 4, dimnames = list(nomi, qualita))
tasso <- c(0.05, 0.002, 0.0005, 0.004)
for (n in nomi) trans[n, ] <- if (substr(n, 1, 1) == substr(n, 3, 3)) 1e6 else round(1e6 * tasso)
err <- loess_monotono(trans)
errori <- err[substr(nomi, 1, 1) != substr(nomi, 3, 3), ]
jsonlite::write_json(list(
  monotona = all(apply(errori, 1, function(r) all(diff(r) <= 1e-15))),
  entro = all(errori >= 1e-7 & errori <= 0.25),
  somme = max(abs(vapply(basi, function(b) max(abs(colSums(err[paste0(b, "2", basi), ]) - 1)), 1))),
  forma = dim(err), nomi = identical(rownames(err), nomi),
  grezzo_risale = tasso[4] > tasso[3]
), uscita, auto_unbox = TRUE)
""", tmp_path)
    assert esito["grezzo_risale"] and esito["monotona"] and esito["entro"]
    assert esito["somme"] < 1e-12 and esito["forma"] == [16, 4] and esito["nomi"]


def test_con_poche_qualita_distinte_s1_avvisa_con_e_s1_03(bioc, tmp_path):
    """
    **Obiettivo**: Verificare che su letture con un solo valore di qualita' S1
    si concluda registrando la degradazione ``E-S1-03``, con i valori trovati
    e l'indicazione di ``err.error_function``, e scriva ``valori_qualita.tsv``;
    e che proseguendo S3 si fermi con ``E-S3-04``, di revisione umana, con
    entrambe le funzioni di errore, e non con un errore generico del ponte.

    **Razionale scientifico e sistemistico**: Le piattaforme a qualita'
    raggruppate danno pochi valori distinti, su cui la stima standard del
    modello di errore non e' affidabile: chi conduce l'analisi deve saperlo
    prima della stima, per scegliere ``err.error_function``.
    """
    scenario = crea_scenario(tmp_path, _campioni(), con_arricchimento=True, con_letture=True)
    run = ProjectRun(scenario.config)
    esito = Esecutore(run, fino_a=Passo.S1).esegui()
    assert esito.conclusione is Conclusione.COMPLETATA
    manifesto = run.albero.manifesto_passo(Passo.S1, Fase.QC_PROFILES)
    (degradazione,) = [d for d in manifesto.degradazioni if d["codice"] == "E-S1-03"]
    assert "1 valori" in degradazione["dettaglio"] and "loess" in degradazione["dettaglio"]
    valori = _tsv(run.albero.cartella(Fase.QC_PROFILES) / "valori_qualita.tsv")
    assert [v["qualita"] for v in valori] == ["40"]

    # Con un solo valore di qualita' nessuna funzione ha una curva da adattare.
    for funzione in ("loess", "loess_monotono"):
        copia = copia_esecuzione((run, esito), tmp_path / funzione, err={"error_function": funzione})
        fermo = Esecutore(copia, fino_a=Passo.S3).esegui()
        assert fermo.conclusione is Conclusione.ARRESTATA, funzione
        assert fermo.punto.passo is Passo.S3 and fermo.punto.codice == "E-S3-04", fermo.punto
        assert fermo.punto.categoria == "revisione_umana"


def test_una_corsa_senza_letture_filtrate_ferma_con_e_s3_03(bioc, ridotta_calcolata, tmp_path):
    """
    **Obiettivo**: Verificare che, se nessun campione di una corsa ha letture
    dopo il filtro, il calcolo di S3 si fermi con ``E-S3-03``, di revisione
    umana, nominando la corsa.

    **Razionale scientifico e sistemistico**: Senza letture filtrate una corsa
    non ha dati su cui stimare il proprio modello di errore: e' una condizione
    dei dati, non un difetto del programma, e ha il suo codice.
    """
    run = copia_esecuzione(ridotta_calcolata, tmp_path)
    corsa = sorted({c.corsa for c in run.valuta().inventario})[0]
    vuoti = {c.accession for c in run.valuta().inventario if c.corsa == corsa}
    tabella = run.albero.cartella(Fase.FILTERED) / "letture_filtrate.tsv"
    righe = tabella.read_text(encoding="utf-8").splitlines()
    tabella.write_text("\n".join(
        "\t".join([*r.split("\t")[:2], "0"]) if r.split("\t")[0] in vuoti else r for r in righe
    ) + "\n", encoding="utf-8")
    fase = run.fase(Passo.S3)
    contesto = run.contesto(Passo.S3, run.valuta()).ristretto(fase.parametri)
    with pytest.raises(ErrorePipeline) as info:
        fase.calcola(contesto)
    assert info.value.codice == "E-S3-03" and info.value.categoria.value == "revisione_umana"
    assert corsa in info.value.dettaglio


def test_la_funzione_di_errore_monotona_si_usa_in_s3(bioc, ridotta_calcolata, tmp_path):
    """
    **Obiettivo**: Verificare che cambiare ``err.error_function`` in
    ``loess_monotono`` rifaccia S3, che la fase si concluda, che la funzione
    arrivi allo script e che i modelli stimati siano diversi da quelli della
    funzione standard.

    **Razionale scientifico e sistemistico**: Il parametro deve arrivare a
    ``learnErrors``: dichiarato e non usato, darebbe il modello standard con
    un'impronta che dice il contrario.
    """
    run = copia_esecuzione(ridotta_calcolata, tmp_path, err={"error_function": "loess_monotono"})
    esito = Esecutore(run, fino_a=Passo.S3).esegui()
    assert esito.conclusione is Conclusione.COMPLETATA
    assert [r.passo for r in esito.eseguite] == [Passo.S3]
    richiesta = json.loads(
        (run.albero.cartella(Fase.ERROR_MODELS) / "rbridge_richiesta_S3.json").read_text())
    assert "loess_monotono" in json.dumps(richiesta)

    def modelli(esecuzione) -> dict[str, str]:
        return {v["nome"]: v["checksum"]
                for v in esecuzione.albero.manifesto_passo(Passo.S3, Fase.ERROR_MODELS).artefatti
                if v["nome"].endswith(".rds")}

    standard, monotona = modelli(ridotta_calcolata[0]), modelli(run)
    assert standard and set(standard) == set(monotona)
    assert all(standard[n] != monotona[n] for n in standard)


# --------------------------------------------------------------------------- #
# 3. Dataset senza controlli: S10, S11, S12                                    #
# --------------------------------------------------------------------------- #


def test_la_colonna_dei_livelli_entra_nell_oggetto_da_sola(tmp_path):
    """
    **Obiettivo**: Verificare che una ``katharoseq.cell_count_column`` non
    elencata in ``out.study_columns`` ne' in ``out.batch_columns`` superi G02
    e G08 se una delle due tabelle la contiene, e che ``colonne_metadati`` la
    porti nell'oggetto dalla tabella in cui si trova; che elencata non venga
    aggiunta due volte.

    **Razionale scientifico e sistemistico**: La calibrazione legge i livelli
    dall'oggetto integrato: farla dipendere da un secondo elenco in cui
    ripetere la colonna e' un modo sicuro di restare senza curva a calcolo
    concluso.
    """
    scenario = crea_scenario(tmp_path, _campioni(), con_arricchimento=True, con_letture=True)
    colonna = scenario.config.katharoseq.cell_count_column
    assert colonna in scenario.config.out.study_columns

    sola = _con(scenario.config, out={"study_columns": []})
    for gate in ("G15", "G02", "G08"):
        assert esegui_gate(gate, Contesto(sola)).superato, gate
    portate = [c for c in colonne_metadati(sola) if c["colonna_originale"] == colonna]
    assert len(portate) == 1 and "studio" in portate[0]["origine"]
    elencata = [c for c in colonne_metadati(scenario.config) if c["colonna_originale"] == colonna]
    assert elencata == portate

    assente = _con(sola, katharoseq={"cell_count_column": "colonna_che_non_esiste"})
    esito = esegui_gate("G02", Contesto(assente))
    assert [v.codice for v in esito.violazioni] == ["E-S0-02"]


def test_la_decontaminazione_per_piastra_richiede_la_colonna_della_piastra(tmp_path):
    """
    **Obiettivo**: Verificare che ``decontam.mode: batch`` con
    ``decontam.batch_column`` nullo sia respinto da G15 con ``E-G15-11``.

    **Razionale scientifico e sistemistico**: Senza la colonna della piastra la
    modalita' per piastra non ha gruppi: prima lo si scopriva in S12, dopo
    tutto il calcolo, con un errore di R.
    """
    scenario = crea_scenario(tmp_path, _campioni(), con_letture=True)
    config = _con(scenario.config, decontam={"mode": "batch", "batch_column": None})
    esito = esegui_gate("G15", Contesto(config))
    assert "E-G15-11" in [v.codice for v in esito.violazioni]
    assert any("decontam.mode" in v.dettaglio for v in esito.violazioni)


#: Un oggetto integrato costruito sui campioni dell'oggetto vero: ``classe`` e
#: ``piastra`` si riscrivono con le espressioni R date, i conteggi restano.
RISCRIVI_OGGETTO = r"""
suppressPackageStartupMessages(library(phyloseq))
ps <- readRDS(oggetto)
dati <- as(sample_data(ps), "data.frame")
m <- as(otu_table(ps), "matrix")
MODIFICA
sample_data(ps) <- sample_data(dati)
otu_table(ps) <- otu_table(m, taxa_are_rows = TRUE)
saveRDS(ps, oggetto)
jsonlite::write_json(list(classi = as.list(table(dati$classe))), uscita, auto_unbox = TRUE)
"""


def _con_oggetto_riscritto(base, cartella: Path, modifica: str, **sovrascrivi):
    """Una copia della catena fino a S10 con l'oggetto integrato riscritto da
    ``modifica`` (codice R su ``dati`` e ``m``); restituisce l'esecuzione.
    """
    run = copia_esecuzione(base, cartella / "copia", **sovrascrivi)
    oggetto = run.albero.cartella(Fase.PHYLOSEQ) / NOME_OGGETTO
    _r(f"oggetto <- {json.dumps(str(oggetto))}\n" + RISCRIVI_OGGETTO.replace("MODIFICA", modifica),
       cartella)
    return run


def _calcola(run, passo: Passo):
    """Il calcolo di una fase sul contesto dell'esecuzione: il contesto, con le
    degradazioni registrate.
    """
    fase = run.fase(passo)
    contesto = run.contesto(passo, run.valuta()).ristretto(fase.parametri)
    fase.calcola(contesto)
    return contesto


def test_senza_controlli_positivi_s11_dichiara_e_s11_05(bioc, oggetto_calcolato, tmp_path):
    """
    **Obiettivo**: Verificare che su un oggetto senza controlli positivi S11
    si concluda senza adattare curve, con ``E-S11-05`` fra le degradazioni e
    non ``E-S11-02``, con ``positivi.tsv`` vuoto e la soglia di ripiego per
    ogni piastra; e che con ``qc.min_reads_mode: fixed`` non dichiari nulla.

    **Razionale scientifico e sistemistico**: Senza controlli positivi non
    c'e' una curva che sia riuscita male: non c'e' nulla da stimare. Sono due
    condizioni diverse, e chi legge il manifesto deve poterle distinguere.
    """
    modifica = 'dati$classe[dati$classe == "controllo_positivo"] <- "biologico"'
    run = _con_oggetto_riscritto(oggetto_calcolato, tmp_path / "a", modifica)
    contesto = _calcola(run, Passo.S11)
    assert [d.codice for d in contesto.degradazioni] == ["E-S11-05"]
    assert "nessun controllo positivo" in contesto.degradazioni[0].dettaglio
    cartella = run.albero.cartella(Fase.CONTROLS)
    assert _tsv(cartella / "positivi.tsv") == [] and _tsv(cartella / "curve.tsv") == []
    soglia = json.loads((cartella / "soglia.json").read_text())
    assert soglia["scelta"] == "nessuno"
    assert {v["stadio"] for v in soglia["per_piastra"].values()} == {"raw"}

    fisso = _con_oggetto_riscritto(oggetto_calcolato, tmp_path / "b", modifica,
                                   qc={"min_reads_mode": "fixed"})
    assert _calcola(fisso, Passo.S11).degradazioni == []


def test_un_controllo_senza_piastra_o_con_una_lettura_non_ferma_s11(
    bioc, oggetto_calcolato, tmp_path
):
    """
    **Obiettivo**: Verificare che S11 si concluda con un codice del catalogo,
    e non con un errore di R, quando un controllo positivo non ha la piastra e
    un altro ha una sola lettura; e quando nessun campione ha la piastra.

    **Razionale scientifico e sistemistico**: La curva e' definita sul
    logaritmo della profondita', che per una lettura vale zero, e le curve per
    piastra si indicizzano sulla piastra: due casi limite di dati reali che
    non devono diventare un arresto senza diagnosi.
    """
    modifica = """
pos <- which(dati$classe == "controllo_positivo")
dati$piastra[pos[1]] <- NA
m[, pos[2]] <- 0L; m[1, pos[2]] <- 1L
"""
    run = _con_oggetto_riscritto(oggetto_calcolato, tmp_path / "a", modifica)
    contesto = _calcola(run, Passo.S11)
    assert {d.codice for d in contesto.degradazioni} <= {"E-S11-02", "E-S11-04"}
    positivi = _tsv(run.albero.cartella(Fase.CONTROLS) / "positivi.tsv")
    assert [p["nella_curva"] for p in positivi if p["profondita"] == "1"] == ["no"]
    assert sum(1 for p in positivi if p["piastra"] == "") == 1

    senza = _con_oggetto_riscritto(oggetto_calcolato, tmp_path / "b", "dati$piastra <- NA_character_")
    contesto = _calcola(senza, Passo.S11)
    soglia = json.loads((senza.albero.cartella(Fase.CONTROLS) / "soglia.json").read_text())
    assert not soglia["per_piastra"] and soglia["senza_piastra"]["valore"] > 0


def test_con_pochi_negativi_s12_dichiara_e_s12_03(bioc, oggetto_calcolato, tmp_path):
    """
    **Obiettivo**: Verificare che su un oggetto senza controlli negativi, e su
    uno con meno negativi di ``decontam.min_blanks``, S12 si concluda senza
    togliere alcuna variante, con ``E-S12-03`` fra le degradazioni e l'esito
    scritto in ``decontam_riepilogo.json``; anche con ``decontam.mode: batch``.

    **Razionale scientifico e sistemistico**: Senza abbastanza negativi la
    prevalenza nei negativi non e' stimabile: l'oggetto passa ai filtri con i
    suoi contaminanti, e un risultato non decontaminato va dichiarato tale.
    """
    casi = {
        "nessuno": 'dati$classe[dati$classe == "controllo_negativo"] <- "biologico"',
        "pochi": 'neg <- which(dati$classe == "controllo_negativo"); '
                 'dati$classe[neg[-(1:2)]] <- "biologico"',
    }
    for nome, modifica in casi.items():
        for modalita in ("aggregate", "batch"):
            run = _con_oggetto_riscritto(oggetto_calcolato, tmp_path / f"{nome}_{modalita}",
                                         modifica, decontam={"mode": modalita})
            contesto = _calcola(run, Passo.S12)
            assert [d.codice for d in contesto.degradazioni] == ["E-S12-03"], (nome, modalita)
            riepilogo = json.loads(
                (run.albero.cartella(Fase.CONTROLS) / "decontam_riepilogo.json").read_text())
            assert riepilogo["contaminanti_rimossi"] == 0
            assert "nessuna decontaminazione" in riepilogo["esito"]
            assert "decontam.min_blanks" in contesto.degradazioni[0].dettaglio


# --------------------------------------------------------------------------- #
# 4. Filtri finali e cartella consegnata                                       #
# --------------------------------------------------------------------------- #


def _s13_su_decontaminato(base, cartella: Path, modifica_r: str, **sovrascrivi):
    """Una copia della catena fino a S12 con l'oggetto decontaminato modificato
    da ``modifica_r`` (che lavora su ``ps``); restituisce l'esecuzione.
    """
    run = copia_esecuzione(base, cartella / "copia", **sovrascrivi)
    oggetto = run.albero.cartella(Fase.CONTROLS) / "ps_decontaminato.rds"
    _r(f"""
suppressPackageStartupMessages(library(phyloseq))
ps <- readRDS({json.dumps(str(oggetto))})
{modifica_r}
saveRDS(ps, {json.dumps(str(oggetto))})
jsonlite::write_json(list(), uscita)
""", cartella)
    return run


def test_i_taxa_si_confrontano_senza_il_prefisso_di_rango(bioc, decontam_calcolata, tmp_path):
    """
    **Obiettivo**: Verificare che, con i nomi della tassonomia riscritti nella
    forma ``p__Nome`` (e i ranghi non assegnati nella forma del solo
    prefisso), S13 rimuova le stesse varianti che rimuove con i nomi senza
    prefisso: per ``filt.exclude_taxa`` e per il phylum non assegnato.

    **Razionale scientifico e sistemistico**: Alcuni riferimenti (Greengenes,
    GTDB) scrivono il rango nel nome: senza togliere il prefisso il filtro
    tassonomico non troverebbe mai cloroplasti e mitocondri, e non lo direbbe.
    """
    genere = _r(f"""
suppressPackageStartupMessages(library(phyloseq))
ps <- readRDS({json.dumps(str(decontam_calcolata[0].albero.cartella(Fase.CONTROLS) / "ps_decontaminato.rds"))})
t <- as(tax_table(ps), "matrix")
jsonlite::write_json(list(genere = names(sort(table(t[, "Genus"]), decreasing = TRUE))[1]),
                     uscita, auto_unbox = TRUE)
""", tmp_path / "genere")["genere"]
    filt = {"exclude_taxa": ["Chloroplast", "Mitochondria", "Eukaryota", genere]}
    prefissi = """
t <- as(tax_table(ps), "matrix")
for (j in seq_len(ncol(t))) {
  p <- paste0(tolower(substr(colnames(t)[j], 1, 1)), "__")
  t[, j] <- ifelse(is.na(t[, j]), p, paste0(p, t[, j]))
}
tax_table(ps) <- tax_table(t)
"""
    rimosse = {}
    for nome, modifica in (("senza", ""), ("con", prefissi)):
        run = _s13_su_decontaminato(decontam_calcolata, tmp_path / nome, modifica, filt=filt,
                                    prev={"apply": False}, qc={"min_reads_final": 1})
        _calcola(run, Passo.S13)
        righe = _tsv(run.albero.cartella(Fase.FINAL_INTERMEDI) / "varianti_rimosse.tsv")
        rimosse[nome] = {r["asv_id"] for r in righe if r["filtro"] == "tassonomico"}
    assert rimosse["senza"] and rimosse["con"] == rimosse["senza"]


def test_senza_il_rango_phylum_s13_dichiara_e_s13_05(bioc, decontam_calcolata, tmp_path):
    """
    **Obiettivo**: Verificare che, con una tassonomia priva del rango Phylum,
    S13 si concluda registrando ``E-S13-05`` e applichi i soli
    ``filt.exclude_taxa``, senza l'errore di R di una colonna inesistente.

    **Razionale scientifico e sistemistico**: ``filt.remove_na_phylum`` non ha
    su che cosa operare: applicarlo a meta' o fermarsi con un indice fuori dai
    limiti nasconderebbe che il filtro non e' stato applicato.
    """
    run = _s13_su_decontaminato(
        decontam_calcolata, tmp_path,
        't <- as(tax_table(ps), "matrix"); '
        'tax_table(ps) <- tax_table(t[, colnames(t) != "Phylum", drop = FALSE])',
    )
    riepilogo = run.albero.cartella(Fase.TAXONOMY) / "riepilogo.json"
    dati = json.loads(riepilogo.read_text())
    dati["ranghi"] = [r for r in dati["ranghi"] if r != "Phylum"]
    riepilogo.write_text(json.dumps(dati), encoding="utf-8")
    contesto = _calcola(run, Passo.S13)
    assert "E-S13-05" in [d.codice for d in contesto.degradazioni]
    righe = _tsv(run.albero.cartella(Fase.FINAL_INTERMEDI) / "varianti_rimosse.tsv")
    assert not [r for r in righe if "phylum non assegnato" in r["motivo"]]


def test_se_nessun_campione_supera_i_filtri_s13_ferma_con_e_s13_04(
    bioc, decontam_calcolata, tmp_path
):
    """
    **Obiettivo**: Verificare che, con ``qc.min_reads_final`` sopra le letture
    di ogni campione, S13 si fermi con ``E-S13-04``, di revisione umana, con
    il numero di esclusi per filtro nel dettaglio.

    **Razionale scientifico e sistemistico**: Un oggetto finale senza campioni
    non e' un risultato: prima era un ``stop`` di R, riportato come errore
    generico del ponte, senza indicare quale filtro avesse tolto tutto.
    """
    run = copia_esecuzione(decontam_calcolata, tmp_path, qc={"min_reads_final": 10**9})
    fase = run.fase(Passo.S13)
    contesto = run.contesto(Passo.S13, run.valuta()).ristretto(fase.parametri)
    with pytest.raises(ErrorePipeline) as info:
        fase.calcola(contesto)
    assert info.value.codice == "E-S13-04" and info.value.categoria.value == "revisione_umana"
    assert "letture_finali" in info.value.dettaglio


def test_s14_toglie_i_file_consegnati_di_un_esecuzione_precedente(
    bioc, finale_calcolata, tmp_path
):
    """
    **Obiettivo**: Verificare che, rifatta S14 con ``out.export_flat`` falso
    su un'esecuzione che aveva gli export e un ``albero.nwk`` rimasto da una
    configurazione precedente e i file che S13 scriveva nella stessa cartella
    prima di averne una propria, in ``12_final`` restino solo i file elencati
    in ``checksum.sha256`` (l'oggetto finale e i controlli), il manifesto e i
    file del ponte, e la sottocartella degli intermedi.

    **Razionale scientifico e sistemistico**: Un export o un albero di una
    configurazione precedente, lasciato nella cartella consegnata, si
    scambierebbe per un risultato di questa esecuzione senza che alcun
    checksum lo copra.
    """
    run = copia_esecuzione(finale_calcolata, tmp_path, out={"export_flat": False})
    finale = run.albero.cartella(Fase.FINAL)
    (finale / "albero.nwk").write_text("(a,b);\n", encoding="utf-8")
    # Una cartella prodotta quando S13 scriveva qui: i suoi file non sono di
    # questa esecuzione.
    for vecchio in ("esclusioni.tsv", "ps_filtrato.rds", "manifest_S13.json"):
        (finale / vecchio).write_text("vecchio\n", encoding="utf-8")
    assert (finale / "conteggi.tsv").is_file()
    esito = Esecutore(run, fino_a=Passo.S14).esegui()
    assert esito.conclusione is Conclusione.COMPLETATA
    assert [r.passo for r in esito.eseguite] == [Passo.S14]
    elencati = {r.split("  ")[1] for r in (finale / "checksum.sha256").read_text().splitlines()}
    assert elencati == {"ps_final.rds", NOME_CONTROLLI}
    consegnati = {p.name for p in finale.iterdir()
                  if p.is_file() and not p.name.startswith(("rbridge_", "manifest"))}
    assert consegnati == elencati | {"checksum.sha256"}
    assert {p.name for p in finale.iterdir() if p.is_dir()} == {"intermedi"}
    assert not (finale / "manifest_S13.json").exists()


# --------------------------------------------------------------------------- #
# 5. Report                                                                    #
# --------------------------------------------------------------------------- #


def test_il_troncamento_suggerito_e_il_minore_fra_lunghezza_e_qualita():
    """
    **Obiettivo**: Verificare che ``troncamento_suggerito`` dia per lunghezza
    il troncamento piu' lungo che scarta non oltre il 5% delle letture di
    biologici e positivi (i negativi non contano), per qualita' l'ultima
    posizione con mediana dei biologici almeno 30, e come suggerito il minore
    dei due; e che senza letture dia valori nulli, e ``_pct`` un trattino per
    una frazione che non esiste.

    **Razionale scientifico e sistemistico**: Chi imposta ``filter.truncLen``
    su un dataset nuovo ha bisogno di un'indicazione ricavata dalle letture,
    con la regola dichiarata: resta un'indicazione, e il report non la applica.
    """
    classi = {"b1": "biologico", "b2": "biologico", "p": "controllo_positivo",
              "n": "controllo_negativo"}
    lunghezze = [
        {"campione": "b1", "lunghezza": "150", "letture": "90"},
        {"campione": "b1", "lunghezza": "140", "letture": "6"},
        {"campione": "b1", "lunghezza": "100", "letture": "4"},
        {"campione": "p", "lunghezza": "150", "letture": "100"},
        {"campione": "n", "lunghezza": "50", "letture": "1000"},
    ]
    qualita = [
        {"campione": c, "posizione": str(p), "mediana": str(q)}
        for c, valori in (("b1", [38, 36, 31, 28]), ("b2", [38, 34, 30, 20]), ("n", [10, 10, 10, 10]))
        for p, q in enumerate(valori, start=1)
    ]
    esito = troncamento_suggerito(lunghezze, qualita, classi)
    # A 140 si scartano 4 letture su 200 (2%); a 150, 10 su 200 (5%): ammesso.
    assert esito == {"per_lunghezza": 150, "per_qualita": 3, "suggerito": 3}
    assert troncamento_suggerito(lunghezze, qualita[:4], classi)["suggerito"] == 3
    vuoto = troncamento_suggerito([], [], classi)
    assert vuoto == {"per_lunghezza": None, "per_qualita": None, "suggerito": None}
    assert _pct(None) == "-" and _pct(0.125) == "12,5%"


# --------------------------------------------------------------------------- #
# 6. La catena intera su insiemi senza controlli, lotto o piastre (criterio G) #
# --------------------------------------------------------------------------- #


def _insieme_ridotto(cartella: Path, tieni, *, senza_lotto: bool = False) -> dict[str, Any]:
    """Un sottoinsieme della versione ridotta: metadati filtrati sui campioni
    scelti da ``tieni`` (una funzione sulla riga di ``selezione.tsv``) e
    collegamenti alle loro letture. Restituisce i percorsi per ``io``.

    Senza il file del lotto la colonna dei livelli dei controlli positivi si
    aggiunge alla tabella di studio, come in un dataset che la porta li'.
    """
    scelti = {r["accession"]: r["campione"] for r in selezione() if tieni(r)}
    metadati = RIDOTTO / "metadati"
    cartella.mkdir(parents=True)
    (cartella / "fastq").mkdir()
    for file in sorted((RIDOTTO / "fastq").iterdir()):
        if file.name.split("_")[0] in scelti:
            (cartella / "fastq" / file.name).symlink_to(file)

    def filtra(nome: str, colonna: str, valori: set[str], aggiunta: dict[str, str] | None = None):
        with open(metadati / nome, encoding="utf-8", newline="") as file:
            righe = list(csv.reader(file, delimiter="\t"))
        indice = righe[0].index(colonna)
        tenute = [righe[0]] + [r for r in righe[1:] if r[indice] in valori]
        if aggiunta is not None:
            tenute = [tenute[0] + ["katharoseq_cell_count"]] + [
                r + [aggiunta.get(r[indice], "")] for r in tenute[1:]]
        with open(cartella / nome, "w", encoding="utf-8", newline="") as file:
            csv.writer(file, delimiter="\t", lineterminator="\n").writerows(tenute)

    lotto = _tsv(metadati / "lotti.tsv")
    cellule = {r["sample_name_ena"]: r["katharoseq_cell_count"] for r in lotto}
    filtra("assay.txt", "Sample Name", set(scelti.values()))
    filtra("studio.txt", "Sample Name", set(scelti.values()), cellule if senza_lotto else None)
    filtra("lotti.tsv", "experiment_accession", set(scelti))
    return {
        "fastq_dir": str(cartella / "fastq"), "assay_table": str(cartella / "assay.txt"),
        "study_table": str(cartella / "studio.txt"),
        "batch_table": None if senza_lotto else str(cartella / "lotti.tsv"),
    }


def _senza_curve(run) -> None:
    """Senza controlli positivi: nessun controllo valutato, nessuna curva."""
    cartella = run.albero.cartella(Fase.CONTROLS)
    assert _tsv(cartella / "positivi.tsv") == [] and _tsv(cartella / "curve.tsv") == []


def _nulla_rimosso(run) -> None:
    """Senza controlli negativi: nessuna variante tolta dalla decontaminazione."""
    riepilogo = json.loads(
        (run.albero.cartella(Fase.CONTROLS) / "decontam_riepilogo.json").read_text())
    assert riepilogo["contaminanti_rimossi"] == 0 and riepilogo["confronto"]["negativi"] == 0


def _livelli_dalla_tabella_di_studio(run) -> None:
    """Senza file del lotto: i livelli dei positivi arrivano a S11 dalla tabella
    di studio, e nessun campione ha piastra o corsa.
    """
    colonne = _tsv(run.albero.cartella(Fase.PHYLOSEQ) / "colonne_metadati.tsv")
    (livelli,) = [c for c in colonne if c["colonna_originale"] == "katharoseq_cell_count"]
    assert "studio" in livelli["origine"]
    positivi = _tsv(run.albero.cartella(Fase.CONTROLS) / "positivi.tsv")
    assert len(positivi) == 9 and all(p["cellule"] and p["piastra"] == "" for p in positivi)
    soglia = json.loads((run.albero.cartella(Fase.CONTROLS) / "soglia.json").read_text())
    assert not soglia["per_piastra"]


def _una_sola_piastra(run) -> None:
    """Una sola piastra: una sola soglia per piastra, e i controlli consegnati."""
    soglia = json.loads((run.albero.cartella(Fase.CONTROLS) / "soglia.json").read_text())
    assert list(soglia["per_piastra"]) == ["10"]
    assert (run.albero.cartella(Fase.FINAL) / NOME_CONTROLLI).is_file()


#: Per ogni insieme: quali campioni tenere, se togliere il file del lotto, i
#: parametri che lo descrivono, le degradazioni attese per fase e la verifica
#: di cio' che l'insieme deve dimostrare.
INSIEMI: dict[str, dict[str, Any]] = {
    "senza_positivi": {
        "tieni": lambda r: r["classe"] != "controllo_positivo",
        "config": {"ctrl": {"positive_values": []},
                   "katharoseq": {"target_taxon": None, "cell_count_column": None}},
        "attese": {Passo.S11: "E-S11-05"},
        "verifica": _senza_curve,
    },
    "senza_negativi": {
        "tieni": lambda r: r["classe"] != "controllo_negativo",
        "config": {"ctrl": {"blank_values": [], "blank_override_column": None,
                            "blank_override_values": []}},
        "attese": {Passo.S12: "E-S12-03"},
        "verifica": _nulla_rimosso,
    },
    "senza_lotto": {
        "tieni": lambda r: True,
        "senza_lotto": True,
        "config": {"out": {"batch_columns": []}, "err": {"batch_column": None},
                   "decontam": {"batch_column": None},
                   "meta": {"batch_key_column": None, "batch_module_column": None}},
        "attese": {},
        "verifica": _livelli_dalla_tabella_di_studio,
    },
    "una_piastra": {
        "tieni": lambda r: r["piastra"] == "10",
        "config": {},
        "attese": {},
        "verifica": _una_sola_piastra,
    },
}


@pytest.mark.parametrize("nome", sorted(INSIEMI))
def test_la_catena_intera_non_da_mai_un_errore_generico_di_r(bioc, tmp_path, nome):
    """
    **Obiettivo**: Verificare che la catena S0-S14, su quattro insiemi ricavati
    dal sottoinsieme ridotto (senza controlli positivi, senza controlli
    negativi, senza file del lotto, con una sola piastra), si concluda con
    tutte le fasi; che nei manifesti ci siano solo codici del catalogo, mai del
    ponte verso R, e fra essi le degradazioni attese; che ogni insieme mostri
    cio' che deve (nessuna curva, nessuna variante rimossa, i livelli letti
    dalla tabella di studio, una sola soglia per piastra); che senza controlli
    ``ps_controlli.rds`` non compaia fra i file consegnati; e che il report si
    generi in ogni caso.

    **Razionale scientifico e sistemistico**: E' il criterio di generalita'
    delle fasi di calcolo: l'assenza di una classe di controlli, del file del
    lotto o di piu' piastre e' una proprieta' legittima di un dataset, e deve
    dare un risultato dichiarato o un arresto diagnosticato, mai un errore di
    R senza codice.
    """
    insieme = INSIEMI[nome]
    percorsi = _insieme_ridotto(tmp_path / "dati", insieme["tieni"],
                                senza_lotto=insieme.get("senza_lotto", False))
    sovrascrivi = {g: dict(v) for g, v in insieme["config"].items()}
    sovrascrivi.setdefault("io", {}).update(percorsi)
    run = ProjectRun(valida(dati_config(tmp_path, **sovrascrivi)))
    esito = Esecutore(run).esegui()

    # Sul sottoinsieme ridotto le quattro catene arrivano in fondo: un arresto,
    # anche con un codice del catalogo, qui sarebbe una regressione.
    assert esito.conclusione is Conclusione.COMPLETATA, esito.punto
    assert [r.passo for r in esito.eseguite] == [p for p in Passo if p is not Passo.S9]
    for passo in Passo:
        if passo is not Passo.S9:
            codici = _degradazioni(run, passo, run.fase(passo).cartella)
            assert not [c for c in codici if c not in CATALOGO or c.startswith("E-R-")], (passo, codici)
    for passo, codice in insieme["attese"].items():
        assert codice in _degradazioni(run, passo, run.fase(passo).cartella), (passo, codice)
    finale = run.albero.cartella(Fase.FINAL)
    elencati = {r.split("  ")[1] for r in (finale / "checksum.sha256").read_text().splitlines()}
    assert (NOME_CONTROLLI in elencati) == (finale / NOME_CONTROLLI).is_file()
    assert (finale / NOME_CONTROLLI).is_file() == any(
        c.classe.e_controllo for c in run.valuta().inventario)
    insieme["verifica"](run)
    chiudi()
    report = genera(run.config.io.out_root)
    assert report.is_file() and "Decisioni prese automaticamente" in report.read_text(encoding="utf-8")
    print(f"\n{nome}: {esito.conclusione.value}, fasi {[str(r.passo) for r in esito.eseguite]}")
