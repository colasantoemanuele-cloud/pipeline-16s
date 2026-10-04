r"""Suite di test della settimana 18: fino_a sulle dipendenze del grafo, fase S10.

1. Inquadramento nel Piano Operativo
------------------------------------
Settimana 18 (W18), Fase F5 (assemblaggio dell'oggetto integrato S10 in
``10_phyloseq/``), preceduta dalla correzione dell'esecutore: con ``fino_a``
si eseguono la fase richiesta e i suoi antenati nel grafo, non tutte le fasi
che la precedono.

2. Moduli sorgente coperti
--------------------------
* ``src/amplicon16s/runner/graph.py`` (``Grafo.antenati``),
  ``src/amplicon16s/runner/executor.py`` (``fino_a``)
* ``src/amplicon16s/steps/s10_phyloseq.py``, ``R/10_phyloseq.R``, ``R/lib/oggetto.R``
* ``src/amplicon16s/config/schema.py`` (``out.sample_id_source``,
  ``out.study_columns``, ``out.batch_columns``),
  ``src/amplicon16s/errors/catalog.py`` (``E-S10-01``, ``E-S10-02``)

3. Cosa valuta questo file
--------------------------
- ``fino_a`` esegue la fase indicata e i soli antenati secondo le dipendenze
  attive del grafo (S2 senza S1); una ripresa completa esegue ancora tutto;
- i nomi delle colonne dei metadati nell'oggetto sono scelti dalla fase e
  sintattici, con la corrispondenza agli originali; una colonna che non si
  puo' portare senza ambiguita' ferma con ``E-S10-02``;
- l'orientamento si verifica sugli identificativi: una matrice quadrata
  trasposta e' intercettata con ``E-S10-01``, sia nella funzione di verifica sia
  come tabella di S7 data allo script di S10;
- l'ordine degli identificativi delle varianti e' deterministico a parita' di
  abbondanza;
- S10 sul sottoinsieme di prova: quattro componenti allineati, campioni
  identificati dall'accession nell'ordine dell'inventario, varianti ``ASV1``,
  ``ASV2``, ... per abbondanza nei biologici, nomi delle colonne inalterati,
  informazioni accessorie per variante, stessi byte in due esecuzioni;
- un campione dell'inventario senza riga nella tabella di S7 entra
  nell'oggetto con conteggi a zero;
- sul dataset completo: dimensioni dell'oggetto e variante ASV1.

4. Comandi Bash e scenari di esecuzione
---------------------------------------
    ``<immagine>`` e' l'immagine del container della pipeline; quella corrente
    e' indicata in ``test.txt``, sezione 1.3.

    1. Modalità locale standard (R di base con jsonlite, senza Bioconductor né
       dati reali):
       pytest tests/test_w18_oggetto.py -v

    2. Modalità container Docker standard (sottoinsieme ridotto con
       R/Bioconductor):
       docker run --rm \
         -e PYTHONPATH=/app/src \
         -e AMPLICON16S_R_DIR=/app/R \
         -v "$(pwd)":/app \
         -w /app \
         <immagine> \
         pytest -o cache_dir=/tmp/.pytest_cache tests/test_w18_oggetto.py -v

    3. Modalità container Docker completa (con i dati reali OSD-734; la
       configurazione e i percorsi che contiene devono stare nella cartella
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
         pytest -o cache_dir=/tmp/.pytest_cache tests/test_w18_oggetto.py -v

5. Risultato atteso
-------------------
30 test totali:
- 24 passed, 6 skipped in ambiente locale standard (~3s): i 6 test che eseguono
  S10 richiedono phyloseq e le fasi a monte, uno anche i dati reali;
- 29 passed, 1 skipped nel container Docker standard sul sottoinsieme ridotto:
  resta saltato il test sui dati reali;
- 30 passed nel container Docker con i dati reali OSD-734.

6. Razionale scientifico e sistemistico
---------------------------------------
- Con ``fino_a`` eseguire fasi che la richiesta non chiede e' un difetto di
  comportamento prima che di tempo: S1 non serve a S2 e alle fasi a valle.
- Una matrice trasposta produce risultati plausibili e completamente errati:
  le dimensioni non la distinguono quando campioni e varianti sono tanti
  quanti, gli identificativi si'.
- Un campione senza letture e' un dato, non un'assenza: farlo sparire fra i
  metadati e i conteggi cambierebbe in silenzio i denominatori a valle.
- ``data.frame`` in R rinomina in silenzio i nomi non sintattici: i nomi
  scelti esplicitamente e verificati, con la corrispondenza agli originali,
  mantengono il legame con la tabella d'origine.
"""

from __future__ import annotations

import csv
import hashlib
import json
import os
import subprocess
from pathlib import Path

import pytest
from conftest import copia_esecuzione, crea_scenario
from sottoinsieme import config_ridotta, motivo_pacchetti_r_assenti
from test_w08_w09_graph_resume import _campioni, _passi

from amplicon16s.errors.exceptions import ErrorePipeline
from amplicon16s.io_layer.artifacts import AlberoOutput, Fase
from amplicon16s.logging.logger import chiudi
from amplicon16s.metadata.models import ClasseCampione
from amplicon16s.rbridge.payload import PREFISSO
from amplicon16s.rbridge.runner import cartella_r, esegui_script, trova_rscript
from amplicon16s.runner.executor import Conclusione, Esecutore
from amplicon16s.runner.graph import GRAFO, Passo
from amplicon16s.runner.project import ProjectRun, StatoPasso
from amplicon16s.steps.s10_phyloseq import (
    COLONNE_INVENTARIO,
    NOME_OGGETTO,
    colonne_metadati,
    nome_nell_oggetto,
)


@pytest.fixture(autouse=True)
def uscite_pulite():
    """Chiude le uscite del log prima e dopo ogni test, perché nessun handler resti
    aperto sulla cartella temporanea.
    """
    chiudi()
    yield
    chiudi()


_MOTIVO_R = motivo_pacchetti_r_assenti("jsonlite")
_MOTIVO_BIOC = motivo_pacchetti_r_assenti(
    "dada2", "ggplot2", "ShortRead", "jsonlite", "phyloseq", "Biostrings"
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
    """Richiede R con phyloseq e i pacchetti delle fasi a monte: salta senza, ma in
    CI fallisce.
    """
    if _MOTIVO_BIOC is not None:
        if os.environ.get("AMPLICON16S_RICHIEDI_BIOC") == "1":
            pytest.fail(f"Bioconductor e' richiesto in questo ambiente: {_MOTIVO_BIOC}")
        pytest.skip(_MOTIVO_BIOC)


def _tsv(percorso: Path) -> list[dict[str, str]]:
    """Le righe di una tabella separata da tabulazioni."""
    with open(percorso, encoding="utf-8", newline="") as file:
        return list(csv.DictReader(file, delimiter="\t"))


def _impronte(cartella: Path) -> dict[str, str]:
    """L'MD5 di ogni file della cartella, esclusi i file del ponte e i manifesti."""
    return {
        p.name: hashlib.md5(p.read_bytes()).hexdigest()
        for p in sorted(cartella.iterdir())
        if not p.name.startswith((PREFISSO, "manifest"))
    }


def _r(codice: str, cartella: Path) -> object:
    """Esegue codice R che scrive un oggetto JSON in ``uscita``, e lo restituisce."""
    uscita = cartella / "uscita.json"
    script = cartella / "lettura.R"
    script.write_text(
        f"uscita <- {json.dumps(str(uscita))}\n{codice}\n", encoding="utf-8"
    )
    subprocess.run(
        [str(trova_rscript()), "--vanilla", str(script)], check=True, capture_output=True,
    )
    return json.loads(uscita.read_text(encoding="utf-8"))


#: Legge l'oggetto integrato e ne descrive ogni slot.
LEGGI_OGGETTO = r"""
suppressPackageStartupMessages(library(phyloseq))
ps <- readRDS(oggetto)
otu <- methods::as(otu_table(ps), "matrix")
dati <- methods::as(sample_data(ps), "data.frame")
jsonlite::write_json(list(
  classe = class(ps)[1],
  taxa_are_rows = taxa_are_rows(ps),
  righe = rownames(otu), colonne = colnames(otu),
  varianti = taxa_names(ps), campioni = sample_names(ps),
  tassonomia_righe = rownames(tax_table(ps)), ranghi = colnames(tax_table(ps)),
  sequenze_nomi = names(refseq(ps)), sequenze = as.character(refseq(ps)),
  albero = !is.null(phy_tree(ps, errorIfNULL = FALSE)),
  colonne_metadati = colnames(dati), righe_metadati = rownames(dati),
  tipi_metadati = as.list(vapply(dati, function(x) class(x)[1], "")),
  per_campione = as.list(colSums(otu)), per_variante = as.list(rowSums(otu)),
  intero = is.integer(otu),
  metadati = dati
), uscita, auto_unbox = TRUE, digits = NA, na = "null")
"""


def _leggi_oggetto(run: ProjectRun, cartella: Path) -> dict:
    """La descrizione dell'oggetto integrato di un'esecuzione."""
    percorso = run.albero.cartella(Fase.PHYLOSEQ) / NOME_OGGETTO
    return _r(f"oggetto <- {json.dumps(str(percorso))}\n{LEGGI_OGGETTO}", cartella)


# --------------------------------------------------------------------------- #
# 1. fino_a: la fase richiesta e i suoi antenati                               #
# --------------------------------------------------------------------------- #


def test_gli_antenati_seguono_le_dipendenze_attive_del_grafo(tmp_path):
    """
    **Obiettivo**: Verificare che ``Grafo.antenati`` restituisca, in ordine, le
    fasi da cui una fase dipende anche indirettamente: S2 ha il solo S0; S1
    non e' antenato di S8 ne' di S10; con la filogenesi disattivata S9 non e'
    antenato di S10, attivata si'.

    **Razionale scientifico e sistemistico**: Gli antenati sono cio' che
    ``fino_a`` esegue: devono seguire le dipendenze di dato, le stesse con cui
    la ripresa giudica cosa rifare.
    """
    config = crea_scenario(tmp_path, _campioni()).config
    assert GRAFO.antenati(Passo.S0, config) == ()
    assert GRAFO.antenati(Passo.S1, config) == (Passo.S0,)
    assert GRAFO.antenati(Passo.S2, config) == (Passo.S0,)
    assert GRAFO.antenati(Passo.S8, config) == (
        Passo.S0, Passo.S2, Passo.S3, Passo.S4, Passo.S5, Passo.S6, Passo.S7,
    )
    assert GRAFO.antenati(Passo.S10, config) == (
        Passo.S0, Passo.S2, Passo.S3, Passo.S4, Passo.S5, Passo.S6, Passo.S7, Passo.S8,
    )
    dati = config.model_dump(mode="python")
    dati["phylo"]["enabled"] = True
    con_albero = type(config).model_validate(dati)
    assert Passo.S9 in GRAFO.antenati(Passo.S10, con_albero)
    assert Passo.S1 not in GRAFO.antenati(Passo.S14, con_albero)


@pytest.mark.parametrize(
    ("fino_a", "attese"),
    [
        (Passo.S0, [Passo.S0]),
        (Passo.S1, [Passo.S0, Passo.S1]),
        (Passo.S2, [Passo.S0, Passo.S2]),
        (Passo.S8, [Passo.S0, Passo.S2, Passo.S3, Passo.S4, Passo.S5, Passo.S6,
                    Passo.S7, Passo.S8]),
    ],
)
def test_fino_a_esegue_soltanto_gli_antenati(tmp_path, fino_a, attese):
    """
    **Obiettivo**: Verificare che l'esecutore con ``fino_a`` esegua, nell'ordine
    del grafo, la fase indicata e i soli suoi antenati: S2 senza S1, S8 senza
    S1; e che l'esecuzione risulti completata.

    **Razionale scientifico e sistemistico**: Prima l'esecutore eseguiva tutte
    le fasi che precedono ``fino_a`` nell'ordine, anche quelle da cui la fase
    non dipende: lavoro non richiesto, e nei test S1 girava a vuoto.
    """
    registro: list[Passo] = []
    run = ProjectRun(crea_scenario(tmp_path, _campioni(), con_letture=True).config,
                     passi=_passi(registro))
    esito = Esecutore(run, fino_a=fino_a).esegui()
    assert esito.conclusione is Conclusione.COMPLETATA
    assert registro == attese
    assert [r.passo for r in esito.eseguite] == attese
    if fino_a is not Passo.S1:
        assert run.valuta().situazioni[Passo.S1].stato is StatoPasso.DA_ESEGUIRE


def test_una_ripresa_completa_esegue_ancora_tutto(tmp_path):
    """
    **Obiettivo**: Verificare che, dopo un'esecuzione con ``fino_a=S8`` che
    lascia S1 da eseguire, una ripresa senza ``fino_a`` esegua S1 e tutte le
    fasi attive rimaste, S10-S14, senza ripetere S0 e S2-S8; e che un'ultima
    ripresa non esegua nulla.

    **Razionale scientifico e sistemistico**: Restringere ``fino_a`` agli
    antenati non deve togliere nulla alla ripresa completa: le fasi non
    richieste restano da eseguire, e si eseguono quando si chiede tutto.
    """
    registro: list[Passo] = []
    run = ProjectRun(crea_scenario(tmp_path, _campioni(), con_letture=True).config,
                     passi=_passi(registro))
    Esecutore(run, fino_a=Passo.S8).esegui()
    registro.clear()
    esito = Esecutore(run).esegui()
    assert esito.conclusione is Conclusione.COMPLETATA
    assert registro == [Passo.S1, Passo.S10, Passo.S11, Passo.S12, Passo.S13, Passo.S14]
    assert run.valuta().completa
    registro.clear()
    assert Esecutore(run).esegui().eseguite == ()
    assert registro == []


# --------------------------------------------------------------------------- #
# 2. Le colonne dei metadati                                                   #
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize(
    ("originale", "atteso"),
    [
        ("Characteristics[Material Type]", "characteristics_material_type"),
        ("Factor Value[Spaceflight]", "factor_value_spaceflight"),
        ("Parameter Value[Sample Preservation Method]",
         "parameter_value_sample_preservation_method"),
        ("katharoseq_cell_count", "katharoseq_cell_count"),
        ("Well-ID", "well_id"),
        ("16S copies", "x_16s_copies"),
        ("NA", "na_"),
    ],
)
def test_il_nome_nell_oggetto_e_sintattico(originale, atteso):
    """
    **Obiettivo**: Verificare che ``nome_nell_oggetto`` produca nomi minuscoli
    con trattini bassi, prefissi un nome che non comincia con una lettera ed
    eviti le parole riservate di R.

    **Razionale scientifico e sistemistico**: Un nome che R non altera e'
    l'unico modo di sapere quale colonna si legge a valle: con
    ``Characteristics[Material Type]`` ``data.frame`` produrrebbe in silenzio
    ``Characteristics.Material.Type.``.
    """
    assert nome_nell_oggetto(originale) == atteso


def test_la_corrispondenza_dei_nomi_conserva_gli_originali(tmp_path):
    """
    **Obiettivo**: Verificare che le colonne dei metadati per la versione
    ridotta siano quelle dell'inventario seguite dalle richieste, che ognuna
    riporti la tabella e la colonna originale, e che ``materiale`` e
    ``classe`` vengano da ``Characteristics[Material Type]``.

    **Razionale scientifico e sistemistico**: Il nome nell'oggetto e' una
    scelta della pipeline: la corrispondenza scritta e' cio' che lega la
    colonna al dato d'origine.
    """
    config = config_ridotta(tmp_path)
    colonne = colonne_metadati(config)
    nomi = [v["colonna"] for v in colonne]
    assert tuple(nomi[: len(COLONNE_INVENTARIO)]) == COLONNE_INVENTARIO
    assert nomi[len(COLONNE_INVENTARIO):] == [
        nome_nell_oggetto(c) for c in (*config.out.study_columns, *config.out.batch_columns)
    ]
    per_nome = {v["colonna"]: v for v in colonne}
    assert per_nome["materiale"]["colonna_originale"] == "Characteristics[Material Type]"
    assert per_nome["classe"]["colonna_originale"] == "Characteristics[Material Type]"
    assert per_nome["factor_value_spaceflight"]["colonna_originale"] == "Factor Value[Spaceflight]"
    assert per_nome["katharoseq_cell_count"]["origine"].startswith("file di arricchimento")


@pytest.mark.parametrize(
    ("sovrascrivi", "motivo"),
    [
        ({"study_columns": ["Factor Value[Inesistente]"]}, "assente"),
        ({"study_columns": ["Term Source REF"]}, "ripetuta 5 volte"),
        ({"study_columns": ["Sample Name", "sample-name"]}, "gia' usato"),
        ({"study_columns": [], "batch_columns": ["run_date", "Run Date"]}, "assente"),
        ({"batch_columns": ["module"], "study_columns": ["Classe"]}, "assente"),
    ],
)
def test_una_colonna_ambigua_ferma_con_e_s10_02(tmp_path, sovrascrivi, motivo):
    """
    **Obiettivo**: Verificare che una colonna richiesta assente, ripetuta
    nell'intestazione (le colonne ``Term Source REF`` delle tabelle ISA) o con
    un nome nell'oggetto gia' usato sollevi ``E-S10-02``.

    **Razionale scientifico e sistemistico**: Scegliere una delle colonne
    ripetute, o lasciarne una sovrascrivere un'altra, attaccherebbe ai campioni
    un dato diverso da quello richiesto.
    """
    config = config_ridotta(tmp_path, out=sovrascrivi)
    with pytest.raises(ErrorePipeline) as info:
        colonne_metadati(config)
    assert info.value.codice == "E-S10-02"
    assert motivo in info.value.dettaglio


def test_senza_file_di_arricchimento_le_sue_colonne_fermano(tmp_path):
    """
    **Obiettivo**: Verificare che, senza ``io.batch_table``, una colonna in
    ``out.batch_columns`` sollevi ``E-S10-02``, e che con l'elenco vuoto le
    colonne siano quelle dell'inventario e della tabella di studio.

    **Razionale scientifico e sistemistico**: Un dato richiesto e mancante non
    diventa una colonna vuota in silenzio.
    """
    with pytest.raises(ErrorePipeline) as info:
        colonne_metadati(config_ridotta(tmp_path, io={"batch_table": None}))
    assert info.value.codice == "E-S10-02"
    config = config_ridotta(tmp_path, io={"batch_table": None}, out={"batch_columns": []})
    assert len(colonne_metadati(config)) == len(COLONNE_INVENTARIO) + len(config.out.study_columns)


# --------------------------------------------------------------------------- #
# 3. L'orientamento e l'ordine, in R                                           #
# --------------------------------------------------------------------------- #


#: Uno script che verifica l'orientamento di una matrice 3 x 3, giusta o trasposta.
VERIFICA_QUADRATA = """
for (f in c("io_json.R", "errors.R", "oggetto.R")) {
  source(file.path(Sys.getenv("AMPLICON16S_R_LIB"), f))
}
esegui_fase(function(parametri, cartella) {
  varianti <- c("ASV1", "ASV2", "ASV3")
  campioni <- c("ERX1", "ERX2", "ERX3")
  nomi <- if (isTRUE(parametri$trasposta)) list(campioni, varianti) else list(varianti, campioni)
  m <- matrix(1:9, nrow = 3, dimnames = nomi)
  verifica_orientamento(m, varianti, campioni, parametri$taxa_are_rows, "tabella di prova")
  character()
})
"""


@pytest.mark.parametrize("taxa_are_rows", [True, False])
def test_una_matrice_quadrata_trasposta_e_intercettata(r, tmp_path, taxa_are_rows):
    """
    **Obiettivo**: Verificare che ``verifica_orientamento`` (``R/lib/oggetto.R``),
    su una matrice 3 x 3 con tanti campioni quante varianti, accetti
    l'orientamento di ``out.taxa_are_rows`` e intercetti quello opposto con
    ``E-S10-01``, dicendo che la matrice e' trasposta.

    **Razionale scientifico e sistemistico**: Con una matrice quadrata le
    dimensioni sono le stesse nei due orientamenti: solo gli identificativi
    distinguono la matrice giusta da quella trasposta.
    """
    script = tmp_path / "quadrata.R"
    script.write_text(VERIFICA_QUADRATA, encoding="utf-8")
    albero = AlberoOutput(tmp_path / "out")
    giusta = esegui_script(
        script, {"trasposta": not taxa_are_rows, "taxa_are_rows": taxa_are_rows},
        albero, Fase.PHYLOSEQ, passo="S10",
    )
    assert giusta.riuscito
    with pytest.raises(ErrorePipeline) as info:
        esegui_script(
            script, {"trasposta": taxa_are_rows, "taxa_are_rows": taxa_are_rows},
            albero, Fase.PHYLOSEQ, passo="S10",
        )
    assert info.value.codice == "E-S10-01"
    assert "trasposta" in info.value.dettaglio


def _ingressi_s10(cartella: Path, trasposta: bool) -> dict:
    """Ingressi sintetici di S10: tre campioni e tre varianti, con la tabella di S7
    nell'orientamento di dada2 o trasposta.
    """
    cartella.mkdir(parents=True, exist_ok=True)
    tabella = cartella / "tabella_asv.rds"
    tassonomia = cartella / "tassonomia.rds"
    _r(f"""
sequenze <- c("ACGT", "CCGT", "GCGT")
campioni <- c("ERX1", "ERX2", "ERX3")
t <- matrix(c(5L, 0L, 1L, 2L, 7L, 0L, 0L, 3L, 9L), nrow = 3,
            dimnames = list(campioni, sequenze))
if ({'TRUE' if trasposta else 'FALSE'}) t <- t(t)
saveRDS(t, {json.dumps(str(tabella))})
tax <- matrix(c("Bacteria", "Bacteria", "Bacteria", "A", "B", NA), nrow = 3,
              dimnames = list(sequenze, c("Kingdom", "Phylum")))
boot <- matrix(c(100L, 100L, 100L, 90L, 80L, 20L), nrow = 3, dimnames = dimnames(tax))
saveRDS(list(tax = tax, boot = boot), {json.dumps(str(tassonomia))})
jsonlite::write_json(list(), uscita)
""", cartella)
    metadati = cartella / "metadati.tsv"
    metadati.write_text("accession\tclasse\nERX1\tbiologico\nERX2\tbiologico\nERX3\tcontrollo_negativo\n",
                        encoding="utf-8")
    return {
        "tabella": str(tabella), "tassonomia": str(tassonomia), "difetti": None,
        "metadati": str(metadati), "colonne": ["accession", "classe"],
        "campioni": ["ERX1", "ERX2", "ERX3"], "biologici": ["ERX1", "ERX2"],
        "taxa_are_rows": True, "prefisso": "ASV",
    }


def test_s10_intercetta_una_tabella_di_s7_quadrata_trasposta(r, tmp_path):
    """
    **Obiettivo**: Verificare che ``R/10_phyloseq.R``, data una tabella di S7
    quadrata (tre campioni, tre varianti) con le varianti sulle righe invece
    dei campioni, si fermi con ``E-S10-01`` senza produrre l'oggetto.

    **Razionale scientifico e sistemistico**: E' lo scenario che un controllo
    sulle dimensioni non vedrebbe: lo script riconosce righe e colonne dagli
    accession dell'inventario e dalle varianti della tassonomia di S8.
    """
    parametri = _ingressi_s10(tmp_path / "ingressi", trasposta=True)
    albero = AlberoOutput(tmp_path / "out")
    with pytest.raises(ErrorePipeline) as info:
        esegui_script(cartella_r() / "10_phyloseq.R", parametri, albero, Fase.PHYLOSEQ, passo="S10")
    assert info.value.codice == "E-S10-01"
    assert "trasposta" in info.value.dettaglio
    assert not (albero.cartella(Fase.PHYLOSEQ) / NOME_OGGETTO).exists()


def test_l_ordine_degli_identificativi_e_deterministico_a_parita(r, tmp_path):
    """
    **Obiettivo**: Verificare che ``ordine_varianti`` ordini per letture nei
    biologici, poi per letture totali, poi per sequenza in ordine
    lessicografico byte per byte, anche in una localizzazione che ordinerebbe
    diversamente le minuscole.

    **Razionale scientifico e sistemistico**: Senza un criterio totale la
    numerazione di due varianti a pari abbondanza dipenderebbe dall'ordine in
    cui arrivano, e cambierebbe fra due esecuzioni.
    """
    libreria = cartella_r() / "lib"
    risultato = _r(f"""
for (f in c("io_json.R", "errors.R", "oggetto.R")) source(file.path({json.dumps(str(libreria))}, f))
Sys.setlocale("LC_COLLATE", "en_US.UTF-8")
sequenze <- c("TTT", "GGG", "aCC", "CCC", "AAA")
bio <- c(5, 10, 10, 10, 0)
tot <- c(50, 12, 10, 10, 99)
jsonlite::write_json(sequenze[ordine_varianti(bio, tot, sequenze)], uscita)
""", tmp_path)
    # GGG ha piu' letture totali di CCC e aCC; fra questi due decide la
    # sequenza, in ordine di byte (la maiuscola prima della minuscola).
    assert risultato == ["GGG", "CCC", "aCC", "TTT", "AAA"]


# --------------------------------------------------------------------------- #
# 4. S10 sul sottoinsieme di prova, nel container                               #
# --------------------------------------------------------------------------- #


def test_s10_assembla_quattro_componenti_allineati(bioc, oggetto_calcolato, tmp_path):
    """
    **Obiettivo**: Verificare che S10 sulla versione ridotta si concluda e
    produca un oggetto phyloseq con conteggi, tassonomia, metadati e
    sequenze allineati: varianti sulle righe con identificativi ``ASV1``...
    ``ASVn``, campioni sulle colonne con gli accession nell'ordine
    dell'inventario, tutti i 28 campioni, nessun albero; sequenze e letture
    coerenti con ``varianti.tsv`` e con la tabella di S7.

    **Razionale scientifico e sistemistico**: ``phyloseq()`` tiene in silenzio
    solo campioni e varianti comuni a tutti i componenti: l'allineamento si
    verifica sugli identificativi, slot per slot.
    """
    run, esito = oggetto_calcolato
    assert esito.conclusione is Conclusione.COMPLETATA
    assert [r.passo for r in esito.eseguite] == [Passo.S10]
    inventario = run.valuta().inventario
    oggetto = _leggi_oggetto(run, tmp_path)
    varianti = _tsv(run.albero.cartella(Fase.PHYLOSEQ) / "varianti.tsv")
    ids = [f"ASV{i}" for i in range(1, len(varianti) + 1)]

    assert oggetto["classe"] == "phyloseq"
    assert oggetto["taxa_are_rows"] is True
    assert oggetto["righe"] == oggetto["varianti"] == oggetto["tassonomia_righe"] == ids
    assert oggetto["sequenze_nomi"] == ids
    assert oggetto["colonne"] == oggetto["campioni"] == oggetto["righe_metadati"]
    assert oggetto["campioni"] == list(inventario.accessioni)
    assert len(oggetto["campioni"]) == 28
    assert oggetto["albero"] is False
    assert oggetto["intero"] is True
    assert oggetto["ranghi"] == ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus"]
    assert [v["asv_id"] for v in varianti] == ids
    assert oggetto["sequenze"] == [v["sequenza"] for v in varianti]
    assert [oggetto["per_variante"][i] for i in ids] == [int(v["letture_totali"]) for v in varianti]
    biologici = [int(v["letture_biologici"]) for v in varianti]
    assert biologici == sorted(biologici, reverse=True)

    riepilogo = json.loads((run.albero.cartella(Fase.PHYLOSEQ) / "riepilogo.json").read_text())
    tracciamento = {
        r["campione"]: int(r["letture"])
        for r in _tsv(run.albero.cartella(Fase.CHIMERA) / "letture_lunghezza.tsv")
    }
    assert oggetto["per_campione"] == tracciamento
    assert riepilogo["letture"] == sum(tracciamento.values())
    assert riepilogo["campioni_aggiunti_a_zero"] == []
    assert riepilogo["prima_variante"]["asv_id"] == "ASV1"
    assert riepilogo["prima_variante"]["sequenza"] == varianti[0]["sequenza"]
    print(f"\nS10 sulla versione ridotta: {dict(esito.eseguite[0].metriche)}")


def test_i_nomi_delle_colonne_dei_metadati_non_sono_alterati(bioc, oggetto_calcolato, tmp_path):
    """
    **Obiettivo**: Verificare che le colonne dei metadati dell'oggetto siano
    esattamente quelle di ``colonne_metadati.tsv``, nello stesso ordine, che
    nessuna contenga un punto o una parentesi, che siano tutte testo, e che i
    valori corrispondano alle tabelle d'origine per ogni campione.

    **Razionale scientifico e sistemistico**: Una colonna rinominata in
    silenzio da R spezzerebbe il legame con la colonna originale; una colonna
    convertita in numero o in fattore cambierebbe il dato.
    """
    run, _ = oggetto_calcolato
    cartella = run.albero.cartella(Fase.PHYLOSEQ)
    corrispondenza = _tsv(cartella / "colonne_metadati.tsv")
    oggetto = _leggi_oggetto(run, tmp_path)
    nomi = [v["colonna"] for v in corrispondenza]
    assert oggetto["colonne_metadati"] == nomi
    assert all(c not in n for n in nomi for c in ".[] ")
    assert set(oggetto["tipi_metadati"].values()) == {"character"}

    studio = {r["Sample Name"]: r for r in _tsv(Path(run.config.io.study_table))}
    lotto = {r["experiment_accession"]: r for r in _tsv(Path(run.config.io.batch_table))}
    for campione, riga in zip(run.valuta().inventario, oggetto["metadati"]):
        assert riga["accession"] == campione.accession
        assert riga["sample_name"] == campione.nome
        assert riga["classe"] == campione.classe.value
        assert riga["materiale"] == studio[campione.nome]["Characteristics[Material Type]"]
        assert riga["piastra"] == lotto[campione.accession]["extraction_plate_num"]
        for voce in corrispondenza[len(COLONNE_INVENTARIO):]:
            origine = studio[campione.nome] if voce["origine"].startswith("tabella") else lotto[campione.accession]
            atteso = origine[voce["colonna_originale"]].strip().strip('"') or None
            assert riga.get(voce["colonna"]) == atteso, (campione.accession, voce["colonna"])


def test_le_informazioni_accessorie_sono_indicizzate_per_variante(bioc, oggetto_calcolato):
    """
    **Obiettivo**: Verificare che ``varianti_accessorie.tsv`` abbia una riga per
    variante, nell'ordine degli identificativi, con il bootstrap di ogni rango
    di S8 e le marcature dei taxa difettosi di ``difetto_riferimento.tsv``.

    **Razionale scientifico e sistemistico**: Il bootstrap servira' ai filtri
    di S13 e all'interpretazione, la marcatura a leggere le assegnazioni sui
    taxa col difetto noto: senza uno slot nell'oggetto, restano accanto,
    raggiungibili dall'identificativo.
    """
    run, _ = oggetto_calcolato
    cartella = run.albero.cartella(Fase.PHYLOSEQ)
    varianti = {v["sequenza"]: v["asv_id"] for v in _tsv(cartella / "varianti.tsv")}
    accessorie = _tsv(cartella / "varianti_accessorie.tsv")
    assert [a["asv_id"] for a in accessorie] == [f"ASV{i}" for i in range(1, len(varianti) + 1)]
    per_id = {a["asv_id"]: a for a in accessorie}
    for riga in _tsv(run.albero.cartella(Fase.TAXONOMY) / "bootstrap.tsv"):
        accessoria = per_id[varianti[riga["sequenza"]]]
        for rango, valore in riga.items():
            if rango != "sequenza":
                assert accessoria[f"boot_{rango}"] == valore
    difetti = _tsv(run.albero.cartella(Fase.TAXONOMY) / "difetto_riferimento.tsv")
    assert difetti, "il riferimento sintetico ha taxa difettosi"
    for d in difetti:
        voce = f"{d['colonna']}:{d['taxon']}({d['rango']})"
        assert voce in per_id[varianti[d["sequenza"]]]["difetto_riferimento"].split(";")
    marcate = {varianti[d["sequenza"]] for d in difetti}
    assert {i for i, a in per_id.items() if a["difetto_riferimento"]} == marcate


def test_un_campione_senza_letture_resta_con_conteggi_a_zero(bioc, tassonomia_calcolata, tmp_path):
    """
    **Obiettivo**: Verificare che, togliendo dalla tabella di S7 la riga di un
    controllo negativo (come accade a un campione senza letture filtrate),
    S10 produca un oggetto che contiene comunque tutti i 28 campioni
    dell'inventario, con quel campione a conteggi nulli e i suoi metadati, e
    che il riepilogo lo elenchi fra i campioni aggiunti a zero e senza
    letture.

    **Razionale scientifico e sistemistico**: Sul dataset di riferimento nessun
    campione resta senza letture, quindi il caso si costruisce: si esegue il
    calcolo di S10 sugli artefatti della versione ridotta con la tabella
    modificata. Un campione non puo' sparire in silenzio fra i metadati e i
    conteggi.
    """
    run = copia_esecuzione(tassonomia_calcolata, tmp_path / "copia")
    inventario = run.valuta().inventario
    bianco = inventario.di_classe(ClasseCampione.CONTROLLO_NEGATIVO)[0].accession
    tabella = run.albero.cartella(Fase.CHIMERA) / "tabella_asv.rds"
    _r(f"""
t <- readRDS({json.dumps(str(tabella))})
saveRDS(t[rownames(t) != {json.dumps(bianco)}, , drop = FALSE], {json.dumps(str(tabella))})
jsonlite::write_json(list(), uscita)
""", tmp_path)

    fase = run.fase(Passo.S10)
    contesto = run.contesto(Passo.S10, run.valuta())
    produzione = fase.calcola(contesto.ristretto(fase.parametri))
    assert produzione.metriche["campioni"] == 28
    assert produzione.metriche["campioni_senza_letture"] == [bianco]

    oggetto = _leggi_oggetto(run, tmp_path)
    assert oggetto["campioni"] == list(inventario.accessioni)
    assert oggetto["per_campione"][bianco] == 0
    assert all(n > 0 for a, n in oggetto["per_campione"].items() if a != bianco)
    riga = next(r for r in oggetto["metadati"] if r["accession"] == bianco)
    assert riga["classe"] == "controllo_negativo"
    riepilogo = json.loads((run.albero.cartella(Fase.PHYLOSEQ) / "riepilogo.json").read_text())
    assert riepilogo["campioni_aggiunti_a_zero"] == [bianco]
    assert riepilogo["campioni_senza_letture"] == [bianco]


def test_s10_da_gli_stessi_byte_in_due_esecuzioni(bioc, oggetto_calcolato, tmp_path):
    """
    **Obiettivo**: Verificare che S10, rieseguita sugli stessi artefatti di
    S0-S8, produca in ``10_phyloseq`` gli stessi byte, oggetto serializzato
    compreso.

    **Razionale scientifico e sistemistico**: L'oggetto integrato e' il punto
    di partenza di tutte le fasi successive: due esecuzioni sugli stessi
    ingressi devono dare lo stesso file, non solo un oggetto equivalente.
    """
    run, _ = oggetto_calcolato
    attese = _impronte(run.albero.cartella(Fase.PHYLOSEQ))
    assert NOME_OGGETTO in attese
    copia = copia_esecuzione(oggetto_calcolato, tmp_path)
    assert copia.valuta().situazioni[Passo.S10].stato is StatoPasso.COMPLETATA
    copia.albero.rimuovi_manifesto_passo(Passo.S10, Fase.PHYLOSEQ)
    esito = Esecutore(copia, fino_a=Passo.S10).esegui()
    assert [r.passo for r in esito.eseguite] == [Passo.S10]
    assert _impronte(copia.albero.cartella(Fase.PHYLOSEQ)) == attese


# --------------------------------------------------------------------------- #
# Dataset completo, nel container                                              #
# --------------------------------------------------------------------------- #


@pytest.mark.dati_reali
def test_s10_sul_dataset_completo(bioc, catena_reale, tmp_path):
    """
    **Obiettivo**: Verificare che sul dataset completo S10 si concluda dopo
    S0-S8 con tutti i 960 campioni e tutte le varianti di S7, e riportarne
    durata, dimensioni e la variante ASV1.

    **Razionale scientifico e sistemistico**: Le dimensioni dell'oggetto
    devono essere quelle dell'inventario e della tabella di S7: nessun
    campione e nessuna variante si perdono nell'assemblaggio.
    """
    run, esito = catena_reale
    assert esito.conclusione is Conclusione.COMPLETATA
    assert Passo.S1 not in [r.passo for r in esito.eseguite]
    s10 = next(r for r in esito.eseguite if r.passo is Passo.S10)
    print(f"\nS10: {s10.secondi} s, {dict(s10.metriche)}")
    assert s10.metriche["campioni"] == 960
    oggetto = _leggi_oggetto(run, tmp_path)
    assert len(oggetto["campioni"]) == 960
    s7 = json.loads((run.albero.cartella(Fase.TAXONOMY) / "riepilogo.json").read_text())
    assert len(oggetto["varianti"]) == s7["varianti"]
