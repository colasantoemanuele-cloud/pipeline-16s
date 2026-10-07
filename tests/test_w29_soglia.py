r"""Suite di test della settimana 29: soglia di profondita' omogenea e letture di coppia.

1. Inquadramento nel Piano Operativo
------------------------------------
Settimana 29 (W29), Fase F8 (generalita'): soglia di profondita' confrontabile
fra le piastre, conteggio dei segni di coppia su tutte le letture, strumento di
sensibilita' senza valori di un dataset, riferimento tassonomico comune.

2. Moduli sorgente coperti
--------------------------
* ``src/amplicon16s/steps/s11_controls.py``, ``R/11_controls.R``, ``R/lib/katharoseq.R``
* ``src/amplicon16s/steps/s13_filtri.py``, ``R/13_filtri.R``
* ``src/amplicon16s/steps/s01_profile.py``, ``src/amplicon16s/io_layer/reads.py``
* ``src/amplicon16s/config/schema.py`` (``qc.min_reads_mode``, ``qc.min_reads_raw`` rimosso)
* ``src/amplicon16s/errors/catalog.py`` (``E-S1-04``, ``E-S11-02``, ``E-S11-05``)
* ``src/amplicon16s/report/builder.py`` (soglia per piastra, troncamento suggerito)
* ``src/amplicon16s/runner/graph.py`` (dipendenze di S11 e S13)
* ``scripts/sensitivity.py``, ``docs/sensibilita/griglie.yaml``
* ``dati/riferimento/``, ``dati/osd734/config_osd734.yaml``, ``dati/osd276/config_osd276.yaml``

3. Cosa valuta questo file
--------------------------
- lo schema respinge ``qc.min_reads_raw`` e ``qc.min_reads_mode: fixed`` con
  un messaggio che spiega il cambiamento, e accetta ``none``;
- la regola della soglia, su oggetti costruiti: piastre con e senza curva
  valida e aggregato non valido (mediana); aggregato valido non preferito
  (le piastre senza curva propria usano l'aggregato); aggregato preferito
  dall'AIC e valido (tutte lo usano); nessuna curva valida (nessuna soglia);
  campione senza piastra; una curva che non converge;
- l'AIC dei due modelli e' calcolato sugli stessi punti, e quello per piastra
  e' la somma degli AIC delle curve;
- in ogni caso la soglia attesa si applica alle letture senza chimere, e
  nessun ramo usa le letture grezze;
- S13 non applica il filtro per profondita' quando S11 non da' una soglia, e
  lo dichiara;
- il report: una riga per piastra e una per i campioni senza piastra, con
  biologici in ingresso e conservati; troncamento suggerito per classe, con
  una classe minoritaria di sole letture corte;
- G07: un nome in tre copie non spegne il riconoscimento delle coppie per
  nome; sotto le venti letture le sole seconde letture non cedono a un record;
- S1 conta i segni di coppia su tutte le letture, li registra in
  ``coppie.tsv`` e ferma con ``E-S1-04`` un file che G07 non vede;
- lo strumento di sensibilita': griglie lette dal file, valore corrente dalla
  configurazione, un arresto non previsto ferma lo strumento;
- le due configurazioni pubblicate puntano al riferimento comune;
- sui dati reali: i 15 file ENA del secondo dataset si fermano tutti in S0 o
  in S1, le forward ricavate e i file del dataset di riferimento no.

4. Comandi Bash e scenari di esecuzione
---------------------------------------
    ``<immagine>`` e' l'immagine del container della pipeline; quella corrente
    e' indicata in ``test.txt``, sezione 1.3.

    1. Modalità locale standard (R di base con jsonlite, senza Bioconductor né
       dati reali):
       pytest tests/test_w29_soglia.py -v

    2. Modalità container Docker standard (sottoinsieme ridotto con
       R/Bioconductor):
       docker run --rm \
         -e PYTHONPATH=/app/src \
         -e AMPLICON16S_R_DIR=/app/R \
         -v "$(pwd)":/app \
         -w /app \
         <immagine> \
         pytest -o cache_dir=/tmp/.pytest_cache tests/test_w29_soglia.py -v

    3. Modalità container Docker completa (con i dati reali OSD-734; i file
       del secondo dataset si indicano con ``AMPLICON16S_DATI_OSD276``, la
       cartella che contiene ``ena/`` e ``fastq/``):
       docker run --rm \
         --memory=24g \
         --memory-swap=24g \
         -e PYTHONPATH=/app/src \
         -e AMPLICON16S_R_DIR=/app/R \
         -e AMPLICON16S_CONFIG_DATI_REALI="$HOME/ASI/config_osd734.yaml" \
         -e AMPLICON16S_DATI_OSD276="$HOME/ASI/w27/dati/osd276" \
         -v "$(pwd)":/app \
         -v "$HOME/ASI":"$HOME/ASI" \
         -w /app \
         <immagine> \
         pytest -o cache_dir=/tmp/.pytest_cache tests/test_w29_soglia.py -v

5. Risultato atteso
-------------------
Vedi ``test.txt``, scheda W29.

6. Razionale scientifico e sistemistico
---------------------------------------
- Una soglia di profondita' che in alcune piastre vale sulle letture senza
  chimere e in altre sulle letture grezze seleziona i campioni in modo diverso
  da piastra a piastra: e' un effetto di lotto introdotto dal filtro, che le
  analisi ecologiche leggerebbero come una differenza fra i campioni.
- La scelta fra curva aggregata e curve per piastra e' un confronto fra
  modelli: vale solo se i due sono adattati sugli stessi punti.
- Un file con le due letture di ogni coppia trattato come single-end da' il
  doppio delle letture, forward e inverse mescolate, senza alcun errore: va
  riconosciuto su tutte le letture, non solo sulle prime.
- Uno strumento di analisi che porta nel codice i valori di un dataset non e'
  riusabile su un altro, e i suoi guasti non devono confondersi con le misure.
"""

from __future__ import annotations

import csv
import gzip
import importlib.util
import json
import math
import os
import statistics
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import pytest
import yaml
from conftest import Campione, copia_esecuzione, crea_scenario, lettura
from sottoinsieme import config_ridotta, dati_config, motivo_pacchetti_r_assenti

from amplicon16s.config.schema import ErroreConfigurazione, carica
from amplicon16s.errors.catalog import CATALOGO, Categoria
from amplicon16s.gates.g01_g15 import Contesto
from amplicon16s.gates.registry import esegui_gate
from amplicon16s.io_layer.artifacts import Fase
from amplicon16s.io_layer.reads import conta_coppie, scansiona_file
from amplicon16s.logging.logger import chiudi
from amplicon16s.metadata.models import ClasseCampione
from amplicon16s.rbridge.runner import cartella_r, trova_rscript
from amplicon16s.report.builder import (
    CARTELLA_REPORT,
    CARTELLA_TABELLE,
    NOME_REPORT,
    genera,
    troncamento_suggerito,
)
from amplicon16s.runner.executor import Conclusione, Esecutore
from amplicon16s.runner.graph import GRAFO, Passo
from amplicon16s.runner.project import ProjectRun
from amplicon16s.steps.s01_profile import NOME_COPPIE, classi_oltre_soglia, letture_corte
from amplicon16s.steps.s10_phyloseq import NOME_OGGETTO
from amplicon16s.steps.s13_filtri import esclusi_per_profondita

RADICE = Path(__file__).resolve().parents[1]
DATI = RADICE / "dati"
senza_dati = pytest.mark.skipif(
    not (DATI / "riferimento").is_dir(),
    reason="la cartella dati/ non e' presente (nell'immagine non viene copiata)",
)


@pytest.fixture(autouse=True)
def uscite_pulite():
    """Chiude le uscite del log prima e dopo ogni test, perché nessun handler resti
    aperto sulla cartella temporanea.
    """
    chiudi()
    yield
    chiudi()


_MOTIVO_BIOC = motivo_pacchetti_r_assenti(
    "dada2", "ggplot2", "ShortRead", "jsonlite", "phyloseq", "Biostrings"
)


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


def _scrivi_record(percorso: Path, record: list[tuple[str, str]]) -> None:
    """Scrive un file FASTQ compresso con le intestazioni e le sequenze indicate."""
    with gzip.open(percorso, "wt", encoding="utf-8") as file:
        for intestazione, sequenza in record:
            file.write(f"@{intestazione}\n{sequenza}\n+\n{'I' * len(sequenza)}\n")


# --------------------------------------------------------------------------- #
# 1. Schema, catalogo e grafo                                                  #
# --------------------------------------------------------------------------- #


def test_min_reads_raw_e_il_modo_fixed_sono_respinti_con_un_messaggio(tmp_path):
    """
    **Obiettivo**: Verificare che una configurazione che dichiara
    ``qc.min_reads_raw``, o ``qc.min_reads_mode: fixed``, sia respinta dallo
    schema con un messaggio che dice che la soglia sulle letture grezze non
    esiste piu' e che cosa la sostituisce; che ``none`` e
    ``katharoseq_if_available`` siano accettati, e il secondo sia il
    predefinito.

    **Razionale scientifico e sistemistico**: Una configurazione scritta per la
    soglia fissa chiedeva un filtro che non c'e' piu': accettarla in silenzio,
    o respingerla come chiave sconosciuta, lascerebbe credere che il filtro sia
    solo cambiato di nome.
    """
    with pytest.raises(ErroreConfigurazione) as raw:
        config_ridotta(tmp_path, qc={"min_reads_raw": 1000})
    testo = "\n".join(raw.value.problemi)
    assert "min_reads_raw non esiste piu'" in testo and "letture grezze" in testo
    assert "min_reads_final" in testo and "mediana" in testo

    with pytest.raises(ErroreConfigurazione) as fisso:
        config_ridotta(tmp_path, qc={"min_reads_mode": "fixed"})
    testo = "\n".join(fisso.value.problemi)
    assert "fixed non esiste piu'" in testo
    assert "katharoseq_if_available" in testo and "none" in testo

    assert config_ridotta(tmp_path).qc.min_reads_mode == "katharoseq_if_available"
    assert config_ridotta(tmp_path, qc={"min_reads_mode": "none"}).qc.min_reads_mode == "none"
    assert not hasattr(config_ridotta(tmp_path).qc, "min_reads_raw")


def test_i_codici_della_soglia_e_delle_coppie_sono_nel_catalogo():
    """
    **Obiettivo**: Verificare che ``E-S11-02`` ed ``E-S11-05`` siano
    degradazioni e descrivano la regola nuova (soglia aggregata o mediana;
    nessuna soglia) senza nominare le letture grezze come ripiego, e che
    ``E-S1-04`` sia un arresto a revisione umana con l'indicazione di G07,
    estrarre le sole letture forward.

    **Razionale scientifico e sistemistico**: Il catalogo e' l'unica fonte dei
    messaggi: un'azione che promettesse ancora la soglia fissa manderebbe chi
    legge a cercare un parametro rimosso.
    """
    ripiego, nessuna, coppie = CATALOGO["E-S11-02"], CATALOGO["E-S11-05"], CATALOGO["E-S1-04"]
    assert ripiego.categoria is Categoria.DEGRADAZIONE_AUTOMATICA
    assert nessuna.categoria is Categoria.DEGRADAZIONE_AUTOMATICA
    assert "aggregata" in ripiego.azione and "mediana" in ripiego.azione
    assert "qc.min_reads_final" in nessuna.azione
    for voce in (ripiego, nessuna):
        assert "min_reads_raw" not in voce.sintesi + voce.azione
    assert coppie.categoria is Categoria.REVISIONE_UMANA and coppie.fase == "S1"
    assert "estrai le sole letture forward" in coppie.azione
    assert "estrai le sole letture forward" in CATALOGO["E-S0-07"].azione


def test_s11_e_s13_non_leggono_piu_le_letture_grezze():
    """
    **Obiettivo**: Verificare che nel grafo S11 dipenda dalla sola S10 e S13
    non dipenda da S2, che S11 non dichiari ``qc.min_reads_raw``, e che il
    sorgente di S11 e di S13 (Python e R) non nomini ne' il parametro rimosso
    ne' la tabella delle letture grezze.

    **Razionale scientifico e sistemistico**: Le dipendenze del grafo sono di
    dato: una fase dipende dalle fasi di cui legge gli artefatti. Tolto il
    ramo sulle letture grezze, una dipendenza rimasta rifarebbe S11 e S13 a
    ogni modifica del filtro senza che nulla di cio' che leggono sia cambiato.
    """
    assert GRAFO.nodo(Passo.S11).dipendenze == (Passo.S10,)
    assert Passo.S2 not in GRAFO.nodo(Passo.S13).dipendenze
    from amplicon16s.runner.project import passi_realizzati

    assert "qc.min_reads_raw" not in passi_realizzati()[Passo.S11].parametri
    import amplicon16s.steps.s11_controls as s11
    import amplicon16s.steps.s13_filtri as s13

    for percorso in (Path(s11.__file__), Path(s13.__file__),
                     cartella_r() / "11_controls.R", cartella_r() / "13_filtri.R"):
        testo = percorso.read_text(encoding="utf-8")
        assert "min_reads_raw" not in testo and "letture_prefiltro" not in testo, percorso.name


# --------------------------------------------------------------------------- #
# 2. La regola della soglia, su oggetti costruiti                              #
# --------------------------------------------------------------------------- #

#: Profondita' di una serie di otto controlli positivi e scostamenti dalla
#: curva, fissi: le curve adattate non sono mai perfette, e l'AIC e' finito.
PROFONDITA = (120000, 64000, 35000, 9000, 4500, 3000, 2200, 2500)
SCOSTAMENTI = (0.004, -0.006, 0.008, -0.005, 0.007, -0.004, 0.003, 0.002)
#: Fedelta' senza relazione con la profondita': una curva converge, con una
#: bonta' molto sotto il minimo.
DISORDINE = (0.7, 0.3, 0.8, 0.1, 0.9, 0.2)
PROFONDITA_DISORDINE = (200000, 100000, 40000, 18000, 8000, 4000)


def _serie(x50: float, profondita=PROFONDITA, pendenza: float = 20.0) -> list[tuple[int, float]]:
    """Una serie di controlli positivi che segue la sigmoide con punto medio
    ``x50`` (in log10 della profondita'), con gli scostamenti fissi.
    """
    punti = []
    for n, scarto in zip(profondita, SCOSTAMENTI, strict=False):
        fedelta = 1.0 / (1.0 + (x50 / math.log10(n)) ** pendenza) + scarto
        punti.append((n, min(max(fedelta, 0.0), 1.0)))
    return punti


#: Costruisce l'oggetto integrato dalla descrizione JSON: due varianti, quella
#: del taxon atteso e un'altra; ogni controllo ha un livello di concentrazione
#: tutto suo, cosi' la conformita' per livello non ne esclude nessuno.
COSTRUISCI = r"""
suppressPackageStartupMessages(library(phyloseq))
d <- jsonlite::read_json(descrizione, simplifyVector = TRUE)
campioni <- d$accession
bersaglio <- as.integer(d$bersaglio); totale <- as.integer(d$letture)
otu <- rbind(Var = bersaglio, Altro = totale - bersaglio)
colnames(otu) <- campioni
tax <- matrix(c("Bacteria", "Bacteria", "Proteobacteria", "Cyanobacteria", "C", "C", "O", "O",
                "Comamonadaceae", "F", d$taxon, "G"), nrow = 2,
              dimnames = list(c("Var", "Altro"), c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")))
dati <- data.frame(accession = campioni, sample_name = campioni, classe = d$classe,
                   piastra = d$piastra, katharoseq_cell_count = d$cellule,
                   row.names = campioni, stringsAsFactors = FALSE)
ps <- phyloseq(otu_table(otu, taxa_are_rows = TRUE), tax_table(tax), sample_data(dati))
saveRDS(ps, d$oggetto)
"""


def _s11_costruito(base, cartella: Path, positivi: dict[str | None, list[tuple[int, float]]],
                   biologici: dict[str | None, list[int]], **sovrascrivi):
    """Sostituisce l'oggetto di S10 con uno costruito ed esegue il calcolo di S11.

    ``positivi`` da' per piastra (``None``: senza piastra) le coppie
    (profondita', fedelta') dei controlli; ``biologici`` per piastra le
    profondita' dei campioni. Restituisce l'esecuzione e il contesto, con le
    degradazioni registrate.
    """
    run = copia_esecuzione(base, cartella / "copia", **sovrascrivi)
    righe: dict[str, list[Any]] = {k: [] for k in (
        "accession", "classe", "piastra", "cellule", "letture", "bersaglio")}
    for piastra, punti in positivi.items():
        for n, fedelta in punti:
            i = len(righe["accession"])
            righe["accession"].append(f"ERXP{i:03d}")
            righe["classe"].append("controllo_positivo")
            righe["piastra"].append(piastra)
            righe["cellule"].append(str(1000 + i))
            righe["letture"].append(n)
            righe["bersaglio"].append(round(n * fedelta))
    for piastra, profondita in biologici.items():
        for n in profondita:
            i = len(righe["accession"])
            righe["accession"].append(f"ERXB{i:03d}")
            righe["classe"].append("biologico")
            righe["piastra"].append(piastra)
            righe["cellule"].append(None)
            righe["letture"].append(n)
            righe["bersaglio"].append(0)
    descrizione = cartella / "descrizione.json"
    descrizione.write_text(json.dumps({
        **righe, "taxon": run.config.katharoseq.target_taxon,
        "oggetto": str(run.albero.cartella(Fase.PHYLOSEQ) / NOME_OGGETTO),
    }), encoding="utf-8")
    script = cartella / "costruisci.R"
    script.write_text(f"descrizione <- {json.dumps(str(descrizione))}\n" + COSTRUISCI,
                      encoding="utf-8")
    subprocess.run([str(trova_rscript()), "--vanilla", str(script)], check=True,
                   capture_output=True)
    fase = run.fase(Passo.S11)
    contesto = run.contesto(Passo.S11, run.valuta()).ristretto(fase.parametri)
    fase.calcola(contesto)
    return run, contesto


def _esiti(run) -> tuple[dict, dict[str, dict[str, str]], list[dict[str, str]]]:
    """``soglia.json``, le curve per modello e la misura per campione di S11."""
    cartella = run.albero.cartella(Fase.CONTROLS)
    soglia = json.loads((cartella / "soglia.json").read_text(encoding="utf-8"))
    curve = {c["modello"]: c for c in _tsv(cartella / "curve.tsv")}
    return soglia, curve, _tsv(cartella / "profondita_campioni.tsv")


def _applicata(soglia: dict, campioni: list[dict[str, str]]) -> dict[str, bool]:
    """Per campione biologico, se il filtro di S13 lo tiene.

    Applica la funzione di S13 alle letture senza chimere misurate da S11 e
    verifica che il suo esito coincida con la misura ``sotto_soglia`` e con il
    confronto diretto fra profondita' e soglia attesa: e' la soglia di
    ``soglia.json``, sullo stadio senza chimere, quella che decide.
    """
    biologici = [r for r in campioni if r["classe"] == ClasseCampione.BIOLOGICO.value]
    letture = {r["accession"]: int(r["letture_nonchimeric"]) for r in biologici}
    tenuti, esclusi = esclusi_per_profondita(
        [(r["accession"], r["piastra"] or None) for r in biologici], soglia, letture)
    assert {e["accession"] for e in esclusi} == {
        r["accession"] for r in biologici if r["sotto_soglia"] == "si"}
    assert all("stadio nonchimeric" in e["motivo"] for e in esclusi)
    esito = {}
    for r in biologici:
        voce = soglia["per_piastra"].get(r["piastra"]) if r["piastra"] else soglia["senza_piastra"]
        assert voce["stadio"] == "nonchimeric" and r["stadio"] == "nonchimeric"
        attesa = voce["valore"] is None or letture[r["accession"]] >= voce["valore"]
        esito[r["accession"]] = r["accession"] in tenuti
        assert esito[r["accession"]] is attesa, r
    return esito


#: Profondita' dei biologici di ogni piastra nei casi costruiti: attorno a
#: tutte le soglie in gioco, cosi' ogni soglia tiene e toglie qualcosa.
BIOLOGICI = [1000, 12000, 14000, 16000, 60000, 140000, 500000]


def test_piastre_senza_curva_valida_e_aggregato_non_valido_usano_la_mediana(
    bioc, oggetto_calcolato, tmp_path
):
    """
    **Obiettivo**: Verificare che, con due piastre dalla curva valida e soglie
    lontane (aggregato non valido), una piastra dalla curva non valida, una
    con due soli controlli, una senza controlli e un campione senza piastra,
    le prime due usino la propria soglia e tutti gli altri la mediana delle
    due, arrotondata all'intero superiore, con il valore non arrotondato
    registrato; che ``E-S11-02`` riporti per ognuno origine e motivo; e che
    il filtro applichi quella soglia alle letture senza chimere.

    **Razionale scientifico e sistemistico**: Una piastra i cui controlli non
    determinano una soglia deve ricevere un filtro dello stesso ordine di
    quello delle altre, sulla stessa grandezza: e' cio' che evita l'effetto di
    lotto del ripiego sulle letture grezze.
    """
    positivi = {
        "A": _serie(3.7),
        "B": _serie(4.6, tuple(10 * n for n in PROFONDITA)),
        "C": list(zip(PROFONDITA_DISORDINE, DISORDINE, strict=True)),
        "D": _serie(3.7)[:2],
    }
    biologici = {p: BIOLOGICI for p in ("A", "B", "C", "D", "E", None)}
    run, contesto = _s11_costruito(oggetto_calcolato, tmp_path, positivi, biologici)
    soglia, curve, campioni = _esiti(run)

    assert set(curve) == {"aggregato", "piastra A", "piastra B", "piastra C"}
    assert curve["piastra A"]["valida"] == curve["piastra B"]["valida"] == "si"
    assert curve["piastra C"]["valida"] == curve["aggregato"]["valida"] == "no"
    assert soglia["scelta"] == "per_piastra" and soglia["modello"]["aggregato_valido"] is False
    proprie = {p: int(curve[f"piastra {p}"]["soglia"]) for p in ("A", "B")}
    mediana = statistics.median(proprie.values())
    assert soglia["mediana"] == {
        "valore_non_arrotondato": mediana, "valore": math.ceil(mediana),
        "piastre": ["A", "B"], "usata": True,
    }
    per_piastra = soglia["per_piastra"]
    for piastra, valore in proprie.items():
        assert (per_piastra[piastra]["valore"], per_piastra[piastra]["origine"]) == (valore, "propria")
    for voce in (per_piastra["C"], per_piastra["D"], per_piastra["E"], soglia["senza_piastra"]):
        assert (voce["valore"], voce["origine"]) == (math.ceil(mediana), "mediana")
        assert "curva aggregata non valida" in voce["motivo"]
    assert "curva della piastra non valida" in per_piastra["C"]["motivo"]
    assert "2 controlli positivi utilizzabili nella piastra" in per_piastra["D"]["motivo"]
    assert "0 controlli positivi utilizzabili nella piastra" in per_piastra["E"]["motivo"]
    assert soglia["senza_piastra"]["motivo"].startswith("campione senza piastra")

    (degradazione,) = [d for d in contesto.degradazioni if d.codice == "E-S11-02"]
    for piastra in ("C", "D", "E", "senza piastra"):
        assert f"{piastra}: soglia mediana ({math.ceil(mediana)} letture nonchimeric)" \
            in degradazione.dettaglio
    assert "A: soglia" not in degradazione.dettaglio

    tenuti = _applicata(soglia, campioni)
    # In ogni piastra la soglia tiene qualcuno ed esclude qualcuno.
    for piastra in ("A", "B", "C", "D", "E", ""):
        esiti = {tenuti[r["accession"]] for r in campioni
                 if r["classe"] == "biologico" and r["piastra"] == piastra}
        assert esiti == {True, False}, piastra


def test_la_mediana_di_un_numero_dispari_di_soglie_e_quella_centrale(
    bioc, oggetto_calcolato, tmp_path
):
    """
    **Obiettivo**: Verificare, con tre piastre dalla curva valida e punti medi
    scelti a distanza (soglie attese dall'equazione attorno a 13.500, 43.000 e
    136.000 letture), che la piastra senza controlli riceva la soglia
    centrale, senza arrotondamenti, e che le tre soglie proprie coincidano
    entro l'1% con quelle calcolate qui dall'equazione della curva.

    **Razionale scientifico e sistemistico**: Le attese degli altri casi sono
    lette da ``curve.tsv``, cioe' dal codice sotto prova: qui la soglia viene
    dall'equazione dichiarata, n* = 10^(x50 (s / (1 - s))^(1 / h)), e la
    mediana da un conto indipendente.
    """
    punti_medi = {"A": 3.7, "B": 4.15, "C": 4.6}
    scala = {"A": 1, "B": 3, "C": 10}
    positivi = {p: _serie(x50, tuple(scala[p] * n for n in PROFONDITA))
                for p, x50 in punti_medi.items()}
    run, _ = _s11_costruito(oggetto_calcolato, tmp_path, positivi,
                            {p: BIOLOGICI for p in ("A", "B", "C", "D")})
    soglia, curve, campioni = _esiti(run)
    attese = {p: 10 ** (x50 * 9 ** (1 / 20)) for p, x50 in punti_medi.items()}
    for piastra, attesa in attese.items():
        voce = soglia["per_piastra"][piastra]
        assert voce["origine"] == "propria"
        assert voce["valore"] == pytest.approx(attesa, rel=0.01), piastra
    centrale = soglia["per_piastra"]["B"]["valore"]
    assert soglia["mediana"]["valore_non_arrotondato"] == soglia["mediana"]["valore"] == centrale
    assert soglia["mediana"]["piastre"] == ["A", "B", "C"]
    assert (soglia["per_piastra"]["D"]["valore"], soglia["per_piastra"]["D"]["origine"]) == (
        centrale, "mediana")
    assert curve["aggregato"]["valida"] == "no"
    _applicata(soglia, campioni)


def test_con_l_aggregato_valido_le_piastre_senza_curva_propria_lo_usano(
    bioc, oggetto_calcolato, tmp_path
):
    """
    **Obiettivo**: Verificare che, con due piastre dalla curva valida e
    vicina (aggregato valido, ma con AIC maggiore del modello per piastra),
    una piastra la cui curva converge e non e' determinata dai dati, una senza
    controlli e un campione senza piastra, le prime due usino la propria
    soglia e gli altri quella aggregata; che la mediana sia registrata e non
    usata; che l'R^2 dichiarato sia quello delle sole piastre con la curva
    propria; e che ``E-S11-02`` nomini chi usa l'aggregata.

    **Razionale scientifico e sistemistico**: Dove i controlli di tutte le
    piastre descrivono bene una sola curva, quella curva e' un'informazione
    migliore della mediana di soglie altrui: viene da osservazioni, non da un
    riassunto.
    """
    positivi = {
        "A": _serie(3.70),
        "B": _serie(3.78, tuple(round(1.2 * n) for n in PROFONDITA)),
        "C": [_serie(3.74)[i] for i in (1, 2, 5, 6)],
    }
    biologici = {p: BIOLOGICI for p in ("A", "B", "C", "D", None)}
    run, contesto = _s11_costruito(oggetto_calcolato, tmp_path, positivi, biologici)
    soglia, curve, campioni = _esiti(run)

    assert curve["aggregato"]["valida"] == "si" and curve["piastra C"]["valida"] == "no"
    assert "non determinata dai dati" in curve["piastra C"]["motivo"]
    modello = soglia["modello"]
    assert modello["aggregato_valido"] is True and modello["preferito_aic"] == "per_piastra"
    assert modello["aic_per_piastra"] < modello["aic_aggregato"]
    assert soglia["scelta"] == "per_piastra"
    aggregata = int(curve["aggregato"]["soglia"])
    per_piastra = soglia["per_piastra"]
    for piastra in ("A", "B"):
        assert per_piastra[piastra]["origine"] == "propria"
        assert per_piastra[piastra]["valore"] == int(curve[f"piastra {piastra}"]["soglia"])
    for voce in (per_piastra["C"], per_piastra["D"], soglia["senza_piastra"]):
        assert (voce["valore"], voce["origine"]) == (aggregata, "aggregata")
    assert soglia["mediana"]["usata"] is False and soglia["mediana"]["piastre"] == ["A", "B"]
    # La bonta' dichiarata e' delle sole piastre A e B: non quella di tutte le
    # curve convergenti, che comprende la C.
    assert modello["r2_dichiarato"] != modello["r2_per_piastra"]
    assert modello["r2_dichiarato"] >= 0.99

    (degradazione,) = [d for d in contesto.degradazioni if d.codice == "E-S11-02"]
    for piastra in ("C", "D", "senza piastra"):
        assert f"{piastra}: soglia aggregata ({aggregata} letture nonchimeric)" in degradazione.dettaglio
    _applicata(soglia, campioni)


def test_l_aic_che_preferisce_l_aggregato_valido_lo_da_a_tutte_le_piastre(
    bioc, oggetto_calcolato, tmp_path
):
    """
    **Obiettivo**: Verificare che, con tre piastre dai controlli identici
    (ognuna con la propria curva valida), l'AIC preferisca l'aggregato, che
    ha gli stessi residui con un terzo dei parametri; che tutte le piastre e
    il campione senza piastra usino la soglia aggregata, anche se hanno una
    curva propria valida; e che l'R^2 dichiarato sia quello dell'aggregato.

    **Razionale scientifico e sistemistico**: Se una sola curva descrive i
    controlli di tutte le piastre quanto le curve separate, le differenze fra
    le soglie per piastra sono rumore di stima: applicarle darebbe filtri
    diversi a piastre che i controlli non distinguono.
    """
    fedelta = (0.98, 0.96, 0.91, 0.78, 0.40, 0.08, 0.02, 0.003)
    serie = list(zip(PROFONDITA, fedelta, strict=True))
    positivi = {"A": serie, "B": serie, "C": serie}
    biologici = {p: BIOLOGICI for p in ("A", "B", "C", None)}
    run, contesto = _s11_costruito(oggetto_calcolato, tmp_path, positivi, biologici)
    soglia, curve, campioni = _esiti(run)

    assert {curve[m]["valida"] for m in curve} == {"si"}
    modello = soglia["modello"]
    assert modello["preferito_aic"] == "aggregato" and soglia["scelta"] == "aggregato"
    # Stessi residui per punto: la differenza e' la sola penalita' dei sei
    # parametri in piu' del modello per piastra (tre per curva).
    assert modello["aic_per_piastra"] - modello["aic_aggregato"] == pytest.approx(12.0, abs=1e-6)
    aggregata = int(curve["aggregato"]["soglia"])
    for voce in (*soglia["per_piastra"].values(), soglia["senza_piastra"]):
        assert (voce["valore"], voce["origine"], voce["stadio"]) == (aggregata, "aggregata", "nonchimeric")
        assert "AIC dell'aggregato" in voce["motivo"]
    assert modello["r2_dichiarato"] == modello["r2_aggregato"]
    assert [d.codice for d in contesto.degradazioni] == ["E-S11-02"]
    _applicata(soglia, campioni)


def test_senza_alcuna_curva_valida_non_c_e_soglia(bioc, oggetto_calcolato, tmp_path):
    """
    **Obiettivo**: Verificare che, con controlli positivi la cui fedelta' non
    dipende dalla profondita' (nessuna curva valida, nemmeno l'aggregata),
    nessuna piastra e nessun campione senza piastra abbia una soglia (origine
    ``nessuna``, valore nullo), che la fase dichiari ``E-S11-05`` e non
    ``E-S11-02``, e che il filtro per profondita' non escluda alcun campione,
    nemmeno il meno profondo.

    **Razionale scientifico e sistemistico**: Senza una curva non c'e' un
    valore da cui ripiegare: una soglia fissa sarebbe una scelta non derivata
    dai dati, su una grandezza diversa. Resta il minimo sulle letture finali.
    """
    positivi = {
        "A": list(zip(PROFONDITA_DISORDINE, DISORDINE, strict=True)),
        "B": list(zip((2 * n for n in PROFONDITA_DISORDINE), reversed(DISORDINE), strict=True)),
    }
    biologici = {p: BIOLOGICI for p in ("A", "B", None)}
    run, contesto = _s11_costruito(oggetto_calcolato, tmp_path, positivi, biologici)
    soglia, curve, campioni = _esiti(run)

    assert {c["valida"] for c in curve.values()} == {"no"}
    assert soglia["scelta"] == "nessuno" and soglia["motivo_scelta"].startswith("nessuna curva valida")
    assert soglia["mediana"] is None and soglia["degradazione"] is False
    for voce in (*soglia["per_piastra"].values(), soglia["senza_piastra"]):
        assert (voce["valore"], voce["origine"]) == (None, "nessuna")
    assert [d.codice for d in contesto.degradazioni] == ["E-S11-05"]
    assert "qc.min_reads_final" in contesto.degradazioni[0].dettaglio
    assert set(_applicata(soglia, campioni).values()) == {True}
    assert {r["soglia"] for r in campioni} == {""}


def test_una_curva_che_non_converge_non_fa_preferire_l_aggregato(bioc, oggetto_calcolato, tmp_path):
    """
    **Obiettivo**: Verificare che, se la curva di una piastra con abbastanza
    punti non converge, l'AIC del modello per piastra non sia calcolabile,
    l'aggregato non sia preferito, la piastra con la curva valida usi la
    propria soglia e quella senza una soglia presa altrove, dichiarata.

    **Razionale scientifico e sistemistico**: Un confronto fra modelli con un
    termine mancante non e' un confronto: preferire l'aggregato per difetto
    toglierebbe la curva propria a una piastra che ne ha una valida.
    """
    positivi = {
        "A": _serie(3.7),
        "B": list(zip((100000, 50000, 20000, 9000, 4000, 2000),
                      (0.2, 0.9, 0.1, 0.8, 0.3, 0.7), strict=True)),
    }
    run, contesto = _s11_costruito(oggetto_calcolato, tmp_path, positivi, {"A": BIOLOGICI, "B": BIOLOGICI})
    soglia, curve, campioni = _esiti(run)
    assert curve["piastra B"]["aic"] == "" and "non convergente" in curve["piastra B"]["motivo"]
    modello = soglia["modello"]
    assert modello["aic_per_piastra"] is None and modello["preferito_aic"] == "per_piastra"
    assert "non calcolabile" in soglia["motivo_scelta"]
    propria = int(curve["piastra A"]["soglia"])
    assert (soglia["per_piastra"]["A"]["valore"], soglia["per_piastra"]["A"]["origine"]) == (
        propria, "propria")
    # L'aggregato, con i punti disordinati della B, non e' valido: la B riceve
    # la mediana delle soglie proprie, che con una sola piastra e' quella di A.
    assert curve["aggregato"]["valida"] == "no"
    assert (soglia["per_piastra"]["B"]["valore"], soglia["per_piastra"]["B"]["origine"]) == (
        propria, "mediana")
    assert soglia["mediana"]["valore_non_arrotondato"] == propria
    assert "E-S11-02" in [d.codice for d in contesto.degradazioni]
    _applicata(soglia, campioni)


def test_l_aic_dei_due_modelli_e_calcolato_sugli_stessi_punti(bioc, oggetto_calcolato, tmp_path):
    """
    **Obiettivo**: Verificare che la curva aggregata e le curve per piastra
    siano adattate sugli stessi punti: i controlli delle sole piastre con
    almeno ``ctrl.min_positives`` punti, senza quelli di una piastra con meno
    punti ne' quelli senza piastra; che l'AIC del modello per piastra sia la
    somma degli AIC delle curve; e che l'AIC di ogni curva sia quello della
    verosimiglianza gaussiana con tre parametri (due della curva e la
    varianza), ricalcolato dai residui.

    **Razionale scientifico e sistemistico**: L'AIC confronta modelli sugli
    stessi dati: una curva aggregata adattata anche sui controlli di una
    piastra che non ha una curva propria avrebbe una verosimiglianza su piu'
    punti, e la differenza fra i due AIC non direbbe nulla.
    """
    positivi = {
        "A": _serie(3.70),
        "B": _serie(3.78, tuple(round(1.2 * n) for n in PROFONDITA)),
        "D": _serie(3.7)[:2],
        None: _serie(3.7)[2:6],
    }
    run, _ = _s11_costruito(oggetto_calcolato, tmp_path, positivi, {"A": BIOLOGICI})
    soglia, curve, _ = _esiti(run)
    positivi_scritti = _tsv(run.albero.cartella(Fase.CONTROLS) / "positivi.tsv")
    # Tutti e 22 i controlli sono utilizzabili, ma i punti comuni sono i 16
    # delle piastre A e B.
    assert sum(p["nella_curva"] == "si" for p in positivi_scritti) == 22
    modello = soglia["modello"]
    assert modello["piastre_comuni"] == ["A", "B"] and modello["punti_comuni"] == 16
    assert set(curve) == {"aggregato", "piastra A", "piastra B"}
    assert int(curve["aggregato"]["punti"]) == 16
    assert int(curve["piastra A"]["punti"]) + int(curve["piastra B"]["punti"]) == 16
    somma = float(curve["piastra A"]["aic"]) + float(curve["piastra B"]["aic"])
    assert modello["aic_per_piastra"] == pytest.approx(somma, rel=1e-5)
    assert modello["aic_aggregato"] == pytest.approx(float(curve["aggregato"]["aic"]), rel=1e-5)
    for nome in curve:
        n, rmse = int(curve[nome]["punti"]), float(curve[nome]["rmse"])
        atteso = n * math.log(rmse ** 2) + n * (1 + math.log(2 * math.pi)) + 2 * 3
        # rmse e' scritto con quattro cifre significative.
        assert float(curve[nome]["aic"]) == pytest.approx(atteso, abs=n * 2e-3), nome


def test_con_il_modo_none_le_curve_restano_diagnostica_e_s13_non_filtra(
    bioc, decontam_calcolata, tmp_path
):
    """
    **Obiettivo**: Verificare sulla catena ridotta che con
    ``qc.min_reads_mode: none`` S11 riporti le curve ma nessuna soglia e
    nessuna degradazione sulla soglia; che S13 non escluda alcun campione per
    profondita' e lo dichiari in ``filtri_riepilogo.json``
    (``profondita_applicata`` falso); e che con il modo predefinito, sulla
    stessa catena, il filtro sia applicato e ogni esclusione per profondita'
    sia motivata sulle letture senza chimere.

    **Razionale scientifico e sistemistico**: Un filtro non applicato deve
    risultare tale negli artefatti: un riepilogo con zero esclusi non
    distingue una soglia che nessuno viola da una soglia che non c'e'.
    """
    senza = copia_esecuzione(decontam_calcolata, tmp_path / "none", qc={"min_reads_mode": "none"})
    esito = Esecutore(senza, fino_a=Passo.S13).esegui()
    assert esito.conclusione is Conclusione.COMPLETATA
    # S12 legge il solo oggetto integrato: il modo della soglia non la rifa'.
    assert [r.passo for r in esito.eseguite] == [Passo.S11, Passo.S13]
    soglia, curve, _ = _esiti(senza)
    assert "piastra 10" in curve and soglia["scelta"] == "nessuno"
    assert "none" in soglia["motivo_scelta"]
    assert {v["origine"] for v in soglia["per_piastra"].values()} == {"nessuna"}
    s11 = senza.albero.manifesto_passo(Passo.S11, Fase.CONTROLS)
    assert not {d["codice"] for d in s11.degradazioni} & {"E-S11-02", "E-S11-05"}
    riepilogo = json.loads((senza.albero.cartella(Fase.FINAL_INTERMEDI)
                            / "filtri_riepilogo.json").read_text(encoding="utf-8"))
    assert riepilogo["profondita_applicata"] is False
    assert riepilogo["campioni"]["esclusi"]["profondita"] == 0

    con = copia_esecuzione(decontam_calcolata, tmp_path / "soglia")
    assert Esecutore(con, fino_a=Passo.S13).esegui().conclusione is Conclusione.COMPLETATA
    riepilogo = json.loads((con.albero.cartella(Fase.FINAL_INTERMEDI)
                            / "filtri_riepilogo.json").read_text(encoding="utf-8"))
    assert riepilogo["profondita_applicata"] is True
    esclusi = [e for e in _tsv(con.albero.cartella(Fase.FINAL_INTERMEDI) / "esclusioni.tsv")
               if e["filtro"] == "profondita"]
    assert len(esclusi) == riepilogo["campioni"]["esclusi"]["profondita"]
    assert all("stadio nonchimeric" in e["motivo"] for e in esclusi)


# --------------------------------------------------------------------------- #
# 3. Report                                                                    #
# --------------------------------------------------------------------------- #


def test_il_report_ha_una_riga_per_piastra_e_una_per_i_campioni_senza_piastra(
    bioc, finale_calcolata, tmp_path
):
    """
    **Obiettivo**: Verificare che la sezione della soglia del report abbia una
    riga per piastra e una per i campioni senza piastra, con soglia, stadio,
    origine, biologici in ingresso, conservati e quota; che i conteggi
    coincidano con l'inventario e con le esclusioni per profondita' di S13; e
    che sopra la tabella compaiano il modello scelto, i due AIC, l'R^2
    dichiarato e la mediana.

    **Razionale scientifico e sistemistico**: Chi valuta il risultato deve
    vedere, piastra per piastra, da dove viene la soglia e quanti campioni
    conserva: e' il confronto che rivela un filtro non omogeneo.
    """
    run = copia_esecuzione(finale_calcolata, tmp_path)
    radice = Path(run.config.io.out_root)
    genera(radice)
    with open(radice / CARTELLA_REPORT / CARTELLA_TABELLE / "soglie_profondita.tsv",
              encoding="utf-8", newline="") as file:
        righe = list(csv.DictReader(file, delimiter="\t"))
    soglia, _, _ = _esiti(run)
    assert [r["piastra"] for r in righe] == ["2", "4", "8", "10", "senza piastra"]
    assert set(righe[0]) >= {"soglia (letture)", "stadio delle letture", "origine",
                             "biologici in ingresso", "biologici conservati", "quota conservata"}
    inventario = run.valuta().inventario
    esclusi = {e["accession"] for e in _tsv(run.albero.cartella(Fase.FINAL_INTERMEDI) / "esclusioni.tsv")
               if e["filtro"] == "profondita"}
    for riga in righe[:-1]:
        voce = soglia["per_piastra"][riga["piastra"]]
        assert (riga["soglia (letture)"], riga["stadio delle letture"], riga["origine"]) == (
            str(voce["valore"]), "nonchimeric", voce["origine"])
        biologici = [c.accession for c in inventario
                     if c.classe is ClasseCampione.BIOLOGICO and c.piastra == riga["piastra"]]
        assert int(riga["biologici in ingresso"]) == len(biologici)
        assert int(riga["biologici conservati"]) == len(set(biologici) - esclusi)
    assert (righe[-1]["biologici in ingresso"], righe[-1]["quota conservata"]) == ("0", "-")
    documento = (radice / CARTELLA_REPORT / NOME_REPORT).read_text(encoding="utf-8")
    for voce in ("AIC della curva aggregata", "AIC del modello per piastra", "R² dichiarato",
                 "mediana delle soglie proprie", "punti comuni ai due modelli"):
        assert voce in documento, voce


def test_il_troncamento_suggerito_giudica_ogni_classe_da_sola():
    """
    **Obiettivo**: Verificare, con 1.000 letture lunghe dei biologici e 30
    letture corte dei controlli positivi, che il troncamento suggerito per
    lunghezza non sia quello lungo (che sulle due classi riunite scarterebbe
    il 2,9%, ma tutte le letture dei positivi) e sia la lunghezza dei
    positivi; che il valore suggerito sia accettato dalla regola di S1 e
    quello lungo no; che la lunghezza senza perdite sia la lettura piu' corta
    delle classi controllate; e che la perdita peggiore per campione sia
    riportata con il nome del campione.

    **Razionale scientifico e sistemistico**: Un'indicazione che la fase
    successiva rifiuterebbe non e' un'indicazione: il report e S1 devono usare
    la stessa regola, ogni classe controllata giudicata da sola.
    """
    classi = {"b1": "biologico", "b2": "biologico", "p": "controllo_positivo",
              "n": "controllo_negativo"}
    lunghezze = [
        {"campione": "b1", "lunghezza": "150", "letture": "480"},
        {"campione": "b1", "lunghezza": "120", "letture": "20"},
        {"campione": "b2", "lunghezza": "150", "letture": "500"},
        {"campione": "p", "lunghezza": "100", "letture": "30"},
        {"campione": "n", "lunghezza": "40", "letture": "5000"},
    ]
    qualita = [{"campione": c, "posizione": str(p), "mediana": "38"}
               for c in ("b1", "b2") for p in range(1, 151)]
    esito = troncamento_suggerito(lunghezze, qualita, classi, 0.05)
    # Riunite: a 150 sono corte 50 letture su 1.030 (4,9%), sotto la quota.
    assert (20 + 30) / 1030 < 0.05
    assert esito["senza_perdite"] == 100
    # A 120 nessuna lettura dei biologici e' piu' corta, ma lo sono tutte
    # quelle dei positivi: l'unica lunghezza ammessa e' 100.
    assert esito["per_lunghezza"] == 100 and esito["suggerito"] == 100
    assert (esito["perdita_peggiore"], esito["campione_peggiore"]) == (0.0, "b1")

    per_campione = {"b1": {150: 480, 120: 20}, "b2": {150: 500}, "p": {100: 30}}
    classi_s1 = {"b1": ClasseCampione.BIOLOGICO, "b2": ClasseCampione.BIOLOGICO,
                 "p": ClasseCampione.CONTROLLO_POSITIVO}
    assert classi_oltre_soglia(letture_corte(per_campione, classi_s1, 150), 0.05) == [
        "controllo_positivo"]
    assert classi_oltre_soglia(letture_corte(per_campione, classi_s1, esito["suggerito"]), 0.05) == []

    # Senza i positivi la lunghezza ammessa e' 150, e il campione che perde di
    # piu' e' b1, con il 4% delle sue letture: la quota vale per la classe.
    senza = troncamento_suggerito(lunghezze[:3], qualita, classi, 0.05)
    assert (senza["senza_perdite"], senza["per_lunghezza"]) == (120, 150)
    assert (senza["perdita_peggiore"], senza["campione_peggiore"]) == (0.04, "b1")
    # Il minore fra lunghezza ammessa e qualita'.
    assert troncamento_suggerito(lunghezze[:3], qualita[:130], classi, 0.05)["suggerito"] == 130


def test_il_report_mostra_le_tre_lunghezze_del_troncamento(bioc, finale_calcolata, tmp_path):
    """
    **Obiettivo**: Verificare che il report della catena ridotta riporti il
    troncamento dichiarato, la lunghezza senza perdite, la lunghezza ammessa
    per classe con ``qc.max_frac_short_reads``, la perdita peggiore per
    singolo campione con il suo nome, e il testo che lo dichiara
    un'indicazione e avverte della perdita diversa da campione a campione.

    **Razionale scientifico e sistemistico**: Una quota ammessa per classe non
    dice quanto perde il singolo campione: troncare oltre la lunghezza senza
    perdite cambia la profondita' in modo diverso da campione a campione.
    """
    run = copia_esecuzione(finale_calcolata, tmp_path)
    radice = Path(run.config.io.out_root)
    genera(radice)
    documento = (radice / CARTELLA_REPORT / NOME_REPORT).read_text(encoding="utf-8")
    for voce in ("lunghezza senza perdite", "lunghezza ammessa", "qc.max_frac_short_reads",
                 "perdita peggiore per singolo campione", "è un'indicazione",
                 "in modo diverso da campione a campione"):
        assert voce in documento, voce


# --------------------------------------------------------------------------- #
# 4. G07 e S1: le due letture di ogni coppia                                   #
# --------------------------------------------------------------------------- #


def _campioni() -> list[Campione]:
    """Due biologici, un positivo e un negativo."""
    return [
        Campione("ERX3000001", "NOD1D4.L1"),
        Campione("ERX3000002", "NOD1D4.L2"),
        Campione("ERX3000003", "POS.P1.1", materiale="Positive Control", posizione="Not Applicable"),
        Campione("ERX3000004", "BLANK.P1.1", materiale="blank control", posizione="Not Applicable"),
    ]


def _file_di(scenario, accession: str) -> Path:
    """Il file FASTQ di un campione dello scenario."""
    (percorso,) = Path(scenario.config.io.fastq_dir).glob(f"*{accession}*")
    return percorso


def test_un_nome_in_tre_copie_non_spegne_il_riconoscimento_per_nome(tmp_path):
    """
    **Obiettivo**: Verificare che un file di coppie riconoscibili dal solo
    nome (ogni nome due volte) resti respinto anche se un nome compare tre
    volte; che un'intestazione uguale per tutte le letture, o vuota, continui
    a passare; e che un file di pochi nomi ripetuti ciascuno molte volte passi.

    **Razionale scientifico e sistemistico**: Un record anomalo isolato non
    deve fermare un dataset single-end, ma nemmeno assolvere un file di
    coppie: la tolleranza vale nei due versi, sulla stessa frazione.
    """
    percorso = tmp_path / "letture.fastq.gz"
    coppie = [(f"r{i}", lettura()) for i in range(20) for _ in range(2)]
    _scrivi_record(percorso, coppie + [("r0", lettura())])
    statistiche = scansiona_file(percorso, 1000, None, None)
    assert (statistiche.nomi_ripetuti, statistiche.nomi_oltre_due) == (19, 1)
    assert "compaiono due volte" in statistiche.coppie_nello_stesso_file

    for intestazione in ("campione", ""):
        _scrivi_record(percorso, [(intestazione, lettura())] * 40)
        assert scansiona_file(percorso, 1000, None, None).coppie_nello_stesso_file is None
    _scrivi_record(percorso, [(f"n{i % 4}", lettura()) for i in range(40)])
    assert scansiona_file(percorso, 1000, None, None).coppie_nello_stesso_file is None


def test_sotto_le_venti_letture_le_seconde_non_cedono_a_un_record(tmp_path):
    """
    **Obiettivo**: Verificare che un file di quindici seconde letture resti
    respinto anche con una sola prima lettura, o con una sola lettura senza
    marcatore; che con due letture diverse passi; e che un file di quindici
    forward con una sola seconda lettura passi, come prima.

    **Razionale scientifico e sistemistico**: Sotto le venti letture il 5%
    vale meno di una lettura: senza un minimo assoluto la tolleranza al
    record isolato varrebbe per assolvere le forward ma non per riconoscere
    le seconde.
    """
    percorso = tmp_path / "letture.fastq.gz"
    seconde = [(f"r{i}/2", lettura()) for i in range(15)]
    for anomalo in ("r99/1", "r99"):
        _scrivi_record(percorso, seconde + [(anomalo, lettura())])
        esito = scansiona_file(percorso, 1000, None, None).coppie_nello_stesso_file
        assert esito is not None and "la seconda lettura di una coppia" in esito, anomalo
    _scrivi_record(percorso, seconde + [("r98", lettura()), ("r99", lettura())])
    assert scansiona_file(percorso, 1000, None, None).coppie_nello_stesso_file is None
    _scrivi_record(percorso, [(f"r{i}/1", lettura()) for i in range(15)] + [("r99/2", lettura())])
    assert scansiona_file(percorso, 1000, None, None).coppie_nello_stesso_file is None


@pytest.mark.parametrize("marcatori", [("/1", "/2"), (" 1:N:0:ACGT", " 2:N:0:ACGT")],
                         ids=["barra", "commento"])
def test_un_file_a_blocchi_sfugge_a_g07_e_non_al_conteggio_su_tutte_le_letture(tmp_path, marcatori):
    """
    **Obiettivo**: Verificare che un file con tutte le prime letture seguite
    da tutte le seconde superi G07 quando le letture ispezionate non arrivano
    al secondo blocco, e che il conteggio su tutte le letture
    (``conta_coppie``) lo riconosca con lo stesso giudizio di G07, nelle due
    convenzioni dei marcatori; e che un file di sole forward non sia toccato.

    **Razionale scientifico e sistemistico**: E' il limite dichiarato di G07:
    il costo di S0 e' costante perche' legge le prime letture, e un deposito
    che riporta le coppie a blocchi ha in testa solo prime letture.
    """
    prima, seconda = marcatori
    scenario = crea_scenario(tmp_path, _campioni(), con_letture=True,
                             sovrascrivi={"qc": {"head_reads": 20}})
    percorso = _file_di(scenario, "ERX3000001")
    _scrivi_record(percorso, [(f"r{i}{prima}", lettura()) for i in range(30)]
                   + [(f"r{i}{seconda}", lettura()) for i in range(30)])
    assert esegui_gate("G07", Contesto(scenario.config)).superato
    tutte = conta_coppie(percorso)
    assert (tutte.letture_esaminate, tutte.prime_di_coppia, tutte.seconde_di_coppia) == (60, 30, 30)
    assert "la prima lettura di una coppia e 30 la seconda" in tutte.coppie_nello_stesso_file
    forward = conta_coppie(_file_di(scenario, "ERX3000002"))
    assert forward.letture_esaminate == 40 and forward.coppie_nello_stesso_file is None


def test_un_file_con_le_due_letture_di_ogni_coppia_ferma_s1(bioc, tmp_path):
    """
    **Obiettivo**: Verificare che, con un file a blocchi che G07 non vede, S0
    si concluda e S1 si fermi con ``E-S1-04``, a revisione umana, nominando il
    file, i conteggi e ``coppie.tsv``, con l'azione di G07 (estrarre le sole
    letture forward), prima di calcolare i profili; e che ``coppie.tsv``
    riporti i conteggi di ogni campione.

    **Razionale scientifico e sistemistico**: S1 legge ogni file per intero:
    e' il primo punto in cui il layout si puo' giudicare su tutte le letture,
    e fermarsi li' costa secondi invece dell'intera catena su un dataset con
    il doppio delle letture.
    """
    scenario = crea_scenario(tmp_path, _campioni(), con_arricchimento=True, con_letture=True,
                             sovrascrivi={"qc": {"head_reads": 20}})
    percorso = _file_di(scenario, "ERX3000002")
    _scrivi_record(percorso, [(f"r{i}/1", lettura()) for i in range(30)]
                   + [(f"r{i}/2", lettura()) for i in range(30)])
    run = ProjectRun(scenario.config)
    esito = Esecutore(run, fino_a=Passo.S1).esegui()
    assert [r.passo for r in esito.eseguite] == [Passo.S0]
    assert esito.conclusione is Conclusione.ARRESTATA
    punto = esito.punto
    assert (punto.passo, punto.codice, punto.categoria) == (Passo.S1, "E-S1-04", "revisione_umana")
    assert "1 file su 4" in punto.dettaglio and percorso.name in punto.dettaglio
    assert "30 intestazioni su 60" in punto.dettaglio and NOME_COPPIE in punto.dettaglio
    assert "estrai le sole letture forward" in punto.azione
    cartella = run.albero.cartella(Fase.QC_PROFILES)
    coppie = {r["campione"]: r for r in _tsv(cartella / NOME_COPPIE)}
    assert coppie["ERX3000002"] == {
        "campione": "ERX3000002", "letture": "60", "prime_di_coppia": "30",
        "seconde_di_coppia": "30", "nomi_ripetuti": "30", "nomi_oltre_due": "0",
        "coppie_nello_stesso_file": "si"}
    assert coppie["ERX3000001"]["coppie_nello_stesso_file"] == "no"
    # Il conteggio precede i profili: la fase si e' fermata prima di R.
    assert not (cartella / "riepilogo.json").exists()
    assert run.albero.manifesto_passo(Passo.S1, Fase.QC_PROFILES) is None


def test_un_file_troncato_oltre_le_letture_ispezionate_ferma_s1(bioc, tmp_path):
    """
    **Obiettivo**: Verificare che un file con l'ultimo record incompleto,
    oltre le letture ispezionate da S0, superi la validazione e fermi S1 con
    ``E-S1-05``, a revisione umana, nominando il file e il punto del guasto,
    sia quando il file contiene le due letture di ogni coppia sia quando e'
    di sole forward; e che ``coppie.tsv`` non venga scritto.

    **Razionale scientifico e sistemistico**: Il lettore a blocchi del profilo
    scarta in silenzio un record incompleto, e il filtro verifica la sola
    decompressione: senza questo arresto un file troncato attraverserebbe la
    catena, e un file di coppie troncato risulterebbe privo di segni di coppia.
    """
    for nome, record in (
        ("coppie", [(f"r{i}/1", lettura()) for i in range(30)]
                   + [(f"r{i}/2", lettura()) for i in range(30)]),
        ("forward", [(f"r{i}", lettura()) for i in range(60)]),
    ):
        scenario = crea_scenario(tmp_path / nome, _campioni(), con_arricchimento=True,
                                 con_letture=True, sovrascrivi={"qc": {"head_reads": 20}})
        percorso = _file_di(scenario, "ERX3000002")
        _scrivi_record(percorso, record)
        with gzip.open(percorso, "at", encoding="utf-8") as file:
            file.write("@ultima\nACGT\n")
        run = ProjectRun(scenario.config)
        esito = Esecutore(run, fino_a=Passo.S1).esegui()
        assert [r.passo for r in esito.eseguite] == [Passo.S0], nome
        punto = esito.punto
        assert (punto.passo, punto.codice, punto.categoria) == (
            Passo.S1, "E-S1-05", "revisione_umana"), nome
        assert percorso.name in punto.dettaglio and "record troncato dopo 60 letture" in punto.dettaglio
        assert not (run.albero.cartella(Fase.QC_PROFILES) / NOME_COPPIE).exists()
        chiudi()


def test_s1_registra_i_conteggi_di_coppia_fra_i_suoi_artefatti(bioc, ridotta_calcolata):
    """
    **Obiettivo**: Verificare che sulla versione ridotta S1 si concluda con
    ``coppie.tsv`` fra gli artefatti del suo manifesto, una riga per campione
    con le letture contate uguali a quelle del profilo, e nessun segno di
    coppia.

    **Razionale scientifico e sistemistico**: Il conteggio e' un artefatto
    come gli altri: ha un checksum, e chi rilegge l'esecuzione puo' verificare
    che il layout single-end sia stato controllato su tutte le letture.
    """
    run, _ = ridotta_calcolata
    manifesto = run.albero.manifesto_passo(Passo.S1, Fase.QC_PROFILES)
    assert NOME_COPPIE in {a["nome"] for a in manifesto.artefatti}
    cartella = run.albero.cartella(Fase.QC_PROFILES)
    coppie = _tsv(cartella / NOME_COPPIE)
    grezze = {r["campione"]: r["letture"] for r in _tsv(cartella / "letture_grezze.tsv")}
    assert {r["campione"]: r["letture"] for r in coppie} == grezze
    assert {r["coppie_nello_stesso_file"] for r in coppie} == {"no"}
    assert {(r["prime_di_coppia"], r["seconde_di_coppia"], r["nomi_ripetuti"]) for r in coppie} == {
        ("0", "0", "0")}


# --------------------------------------------------------------------------- #
# 5. Lo strumento di sensibilita'                                              #
# --------------------------------------------------------------------------- #


def _strumento():
    """Il modulo ``scripts/sensitivity.py``, caricato dal file: la cartella degli
    script non e' un pacchetto.
    """
    specifica = importlib.util.spec_from_file_location(
        "sensitivity_w29", RADICE / "scripts" / "sensitivity.py")
    modulo = importlib.util.module_from_spec(specifica)
    specifica.loader.exec_module(modulo)
    return modulo


senza_script = pytest.mark.skipif(
    not (RADICE / "scripts" / "sensitivity.py").is_file(),
    reason="la cartella scripts/ non e' presente (nell'immagine non viene copiata)",
)


@senza_script
def test_lo_strumento_di_sensibilita_legge_griglie_e_valori_correnti(tmp_path):
    """
    **Obiettivo**: Verificare che lo strumento non porti nel codice griglie,
    valori correnti ne' modalita'; che legga le griglie dal file indicato, il
    valore corrente dalla configurazione e riconosca la modalita' corrente dai
    parametri; e che una configurazione che corrisponde a nessuna o a piu'
    modalita' sia un errore dello strumento.

    **Razionale scientifico e sistemistico**: Uno strumento con i valori di un
    dataset nel codice, applicato a un altro, confronterebbe la configurazione
    con una griglia che non la riguarda e fallirebbe su un'asserzione.
    """
    S = _strumento()
    for nome in ("GRIGLIE", "FUORI_GRIGLIA", "MODALITA", "MODALITA_CORRENTE"):
        assert not hasattr(S, nome), nome
    griglie = tmp_path / "griglie.yaml"
    griglie.write_text(yaml.safe_dump({
        "numerici": {"prev.min_fraction": {"griglia": [0.1, 0.2, 0.3], "fuori_griglia": [0.5]},
                     "filt.remove_na_phylum": {"griglia": [False, True]}},
        "modalita": {"nome": "confronto", "valori": {
            "uno": {"decontam.mode": "aggregate"},
            "due": {"decontam.mode": "batch", "decontam.batch_combine": "fisher"}}},
    }, sort_keys=False), encoding="utf-8")
    numerici, nome, modalita = S.leggi_griglie(griglie)
    assert numerici["prev.min_fraction"] == {"griglia": (0.1, 0.2, 0.3), "fuori_griglia": (0.5,)}
    assert numerici["filt.remove_na_phylum"]["fuori_griglia"] == ()
    assert nome == "confronto" and list(modalita) == ["uno", "due"]

    config = config_ridotta(tmp_path, prev={"min_fraction": 0.2})
    assert S.valore_corrente(config, "prev.min_fraction") == 0.2
    assert S.modalita_corrente(config, modalita) == "uno"
    with pytest.raises(S.ErroreSensibilita, match="0 delle modalita'"):
        S.modalita_corrente(config, {"due": modalita["due"]})
    with pytest.raises(S.ErroreSensibilita, match="2 delle modalita'"):
        S.modalita_corrente(config, {"uno": modalita["uno"], "anche": {"decontam.mode": "aggregate"}})
    assert S._per_gruppo({"decontam.mode": "batch", "decontam.batch_combine": "fisher",
                          "qc.max_frac_contaminant": 1.0}) == {
        "decontam": {"mode": "batch", "batch_combine": "fisher"}, "qc": {"max_frac_contaminant": 1.0}}


def _esecuzione_finta(cartella: Path, degradazioni: list[str]) -> Path:
    """Una cartella con i soli artefatti che lo strumento di sensibilita' legge."""
    (cartella / "11_controls").mkdir(parents=True)
    (cartella / "12_final" / "intermedi").mkdir(parents=True)
    (cartella / "11_controls" / "soglia.json").write_text(json.dumps({
        "scelta": "per_piastra", "per_piastra": {
            "2": {"valore": 900, "origine": "mediana"}, "10": {"valore": 800, "origine": "propria"}},
    }), encoding="utf-8")
    (cartella / "11_controls" / "decontam_riepilogo.json").write_text(json.dumps({
        "contaminanti_rimossi": 3, "letture_rimosse": {"biologico": 0.05}}), encoding="utf-8")
    (cartella / "12_final" / "tassonomia.tsv").write_text("asv_id\nASV1\nASV2\n", encoding="utf-8")
    (cartella / "12_final" / "conteggi.tsv").write_text("asv_id\tA\tB\n", encoding="utf-8")
    (cartella / "12_final" / "intermedi" / "filtri_riepilogo.json").write_text(json.dumps({
        "letture": {"finali": 100}, "varianti": {"finali": 2, "rimosse": {"prevalenza": 1}},
        "campioni": {"finali": 2, "biologici": 4, "esclusi": {"profondita": 2}},
    }), encoding="utf-8")
    (cartella / "12_final" / "intermedi" / "manifest_S13.json").write_text(json.dumps({
        "degradazioni": [{"codice": c} for c in degradazioni]}), encoding="utf-8")
    return cartella


@senza_script
def test_solo_un_arresto_previsto_dalla_regola_rende_un_valore_non_ammissibile(
    tmp_path, monkeypatch
):
    """
    **Obiettivo**: Verificare che lo strumento giudichi non ammissibile un
    valore solo per le condizioni che la regola prevede (arresto ``E-S12-02``
    o ``E-S14-01``, campioni svuotati dal filtro di prevalenza dichiarati con
    ``E-S13-02``), e che un arresto con qualunque altro codice lo fermi con
    un errore che nomina la variante e il codice; che le origini non proprie
    delle soglie siano riportate per piastra, in ordine numerico.

    **Razionale scientifico e sistemistico**: Un guasto contato come "valore
    non ammissibile" finirebbe nella tabella come una proprieta' del
    parametro: la regola sceglierebbe un valore in base a un errore.
    """
    S = _strumento()
    assert set(S.ARRESTI_PREVISTI) == {"E-S12-02", "E-S14-01"}
    assert S.SVUOTATI_DALLA_PREVALENZA == "E-S13-02"
    assert {S.ARRESTI_PREVISTI.keys() <= set(CATALOGO), S.SVUOTATI_DALLA_PREVALENZA in CATALOGO} == {True}

    completa = S.misura(_esecuzione_finta(tmp_path / "a", []))
    assert completa["ammissibile"] and completa["esito"] == "completata"
    assert completa["origini"] == {"2": "mediana", "10": "propria"}
    svuotati = S.misura(_esecuzione_finta(tmp_path / "b", ["E-S13-03", "E-S13-02"]))
    assert not svuotati["ammissibile"] and "E-S13-02" in svuotati["esito"]
    fermata = S.misura(_esecuzione_finta(tmp_path / "c", []), "E-S12-02")
    assert not fermata["ammissibile"] and "campioni" not in fermata
    riga = S._riga("p", 1, completa, completa)
    assert riga[-3:-1] == ["2: mediana", "2: 900; 10: 800"]
    assert S._riga("p", 1, fermata, completa)[3:13] == [""] * 10

    origine = _esecuzione_finta(tmp_path / "origine", [])
    base = dati_config(tmp_path / "base")

    def esecutore_con(codice: str | None, dettaglio: str = ""):
        punto = None if codice is None else SimpleNamespace(
            passo=Passo.S12, codice=codice, dettaglio=dettaglio)
        esito = SimpleNamespace(
            conclusione=Conclusione.COMPLETATA if codice is None else Conclusione.ARRESTATA,
            punto=punto, eseguite=[])
        return lambda run: SimpleNamespace(esegui=lambda: esito)

    import amplicon16s.runner.executor as modulo_esecutore

    monkeypatch.setattr(modulo_esecutore, "Esecutore", esecutore_con("E-S12-02"))
    prevista = S.esegui_variante("prevista", base, origine, tmp_path / "lavoro", {}, False)
    assert not prevista["ammissibile"] and prevista["codice_arresto"] == "E-S12-02"
    monkeypatch.setattr(modulo_esecutore, "Esecutore", esecutore_con("E-R-03"))
    with pytest.raises(S.ErroreSensibilita, match="imprevista.*E-R-03.*non e' un arresto previsto"):
        S.esegui_variante("imprevista", base, origine, tmp_path / "lavoro", {}, False)
    # E-S14-01 e' anche un oggetto finale non valido: e' previsto dalla regola
    # solo quando riguarda la frazione di letture trattenute.
    monkeypatch.setattr(modulo_esecutore, "Esecutore", esecutore_con(
        "E-S14-01", "frazione 0.31 sotto qc.min_frac_reads_retained (0.4)"))
    trattenute = S.esegui_variante("trattenute", base, origine, tmp_path / "lavoro", {}, False)
    assert trattenute["codice_arresto"] == "E-S14-01" and not trattenute["ammissibile"]
    monkeypatch.setattr(modulo_esecutore, "Esecutore", esecutore_con(
        "E-S14-01", "ps_filtrato.rds di S13 manca o non corrisponde al suo manifesto"))
    with pytest.raises(S.ErroreSensibilita, match="strutturale.*E-S14-01.*ps_filtrato"):
        S.esegui_variante("strutturale", base, origine, tmp_path / "lavoro", {}, False)
    monkeypatch.setattr(modulo_esecutore, "Esecutore", esecutore_con(None))
    assert S.esegui_variante("conclusa", base, origine, tmp_path / "lavoro", {}, False)["ammissibile"]


# --------------------------------------------------------------------------- #
# 6. Il riferimento tassonomico comune                                         #
# --------------------------------------------------------------------------- #


@senza_dati
def test_le_configurazioni_pubblicate_puntano_al_riferimento_comune():
    """
    **Obiettivo**: Verificare che le configurazioni dei due dataset indichino
    il riferimento tassonomico in ``dati/riferimento/``, con lo stesso MD5
    dello script di scarico, di ``impronte.md5`` e di ``FONTI.tsv``; che la
    configurazione del secondo dataset non indichi alcun percorso nella
    cartella del primo; e che nessuna delle due dichiari i parametri rimossi.

    **Razionale scientifico e sistemistico**: Il riferimento non appartiene a
    un dataset: una configurazione che lo cerca nella cartella di un altro
    dataset non si puo' usare senza scaricare quel dataset.
    """
    sys.path.insert(0, str(DATI / "riferimento"))
    try:
        import scarica_riferimento
    finally:
        sys.path.pop(0)
    assert scarica_riferimento.CARTELLA == DATI / "riferimento"
    impronte = dict(reversed(r.split("  ")) for r in
                    (DATI / "riferimento" / "impronte.md5").read_text(encoding="utf-8").splitlines())
    assert impronte == {nome: voce.md5 for nome, voce in scarica_riferimento.FILE.items()}
    fonti = {r["percorso"]: r["controllo"] for r in _tsv(DATI / "riferimento" / "FONTI.tsv")}
    assert fonti == {nome: f"md5:{md5}" for nome, md5 in impronte.items()}
    assert not [r for r in _tsv(DATI / "osd734" / "FONTI.tsv") if "riferimento" in r["percorso"]]

    for dataset in ("osd734", "osd276"):
        percorso = DATI / dataset / f"config_{dataset}.yaml"
        config = carica(percorso)
        assert config.tax.ref_fasta.parent == Path("dati/riferimento"), dataset
        assert config.tax.ref_md5 == impronte[config.tax.ref_fasta.name]
        assert "min_reads_raw" not in percorso.read_text(encoding="utf-8")
        assert config.qc.min_reads_mode == "katharoseq_if_available"
    assert "osd734" not in (DATI / "osd276" / "config_osd276.yaml").read_text(encoding="utf-8")


# --------------------------------------------------------------------------- #
# 7. Dati reali                                                                #
# --------------------------------------------------------------------------- #

_CONFIG_REALE = os.environ.get("AMPLICON16S_CONFIG_DATI_REALI")
_OSD276 = os.environ.get("AMPLICON16S_DATI_OSD276")


def _giudizi(cartella: Path, head_reads: int) -> tuple[list[str], list[str], list[str]]:
    """I file che G07 respinge, quelli che respinge solo S1 e quelli che passano."""
    g07, s1, passano = [], [], []
    for percorso in sorted(cartella.glob("*.fastq.gz")):
        if scansiona_file(percorso, head_reads, None, None).coppie_nello_stesso_file:
            g07.append(percorso.name)
        elif conta_coppie(percorso).coppie_nello_stesso_file:
            s1.append(percorso.name)
        else:
            passano.append(percorso.name)
    return g07, s1, passano


@pytest.mark.dati_reali
def test_i_file_di_ena_del_secondo_dataset_si_fermano_tutti_in_s0_o_in_s1(tmp_path):
    """
    **Obiettivo**: Verificare che i 15 file che ENA distribuisce per il
    secondo dataset, ciascuno con le due letture di ogni coppia, si fermino
    tutti: in S0 quelli che G07 riconosce dalle prime letture, in S1 gli
    altri, sul conteggio di tutte le letture; e che dei 15 file di letture
    forward ricavate non se ne fermi nessuno.

    **Razionale scientifico e sistemistico**: E' il caso reale che ha
    mostrato il limite di G07: otto dei quindici file cominciano con le prime
    letture, e presi da soli superavano la validazione.
    """
    if not _OSD276 or not (Path(_OSD276) / "ena").is_dir():
        pytest.skip("AMPLICON16S_DATI_OSD276 non impostata: file del secondo dataset non disponibili")
    head_reads = config_ridotta(tmp_path).qc.head_reads
    g07, s1, passano = _giudizi(Path(_OSD276) / "ena", head_reads)
    assert len(g07) + len(s1) == 15 and passano == []
    assert s1, "nessun file a blocchi: il conteggio di S1 non sarebbe messo alla prova"
    g07, s1, passano = _giudizi(Path(_OSD276) / "fastq", head_reads)
    assert (g07, s1) == ([], []) and len(passano) == 15
    print(f"\nfile ENA fermati da G07 e da S1: {15 - len(s1)} e {len(s1)}")


@pytest.mark.dati_reali
def test_nessun_file_del_dataset_di_riferimento_ha_segni_di_coppia():
    """
    **Obiettivo**: Verificare che nessuno dei file di letture del dataset di
    riferimento sia respinto dal conteggio su tutte le letture, e che nessuno
    porti un marcatore di coppia o un nome ripetuto.

    **Razionale scientifico e sistemistico**: Il controllo nuovo di S1 non
    deve fermare un dataset single-end legittimo: e' la verifica di assenza di
    falsi positivi sul dataset su cui la pipeline e' collaudata.
    """
    if not _CONFIG_REALE or not Path(_CONFIG_REALE).is_file():
        pytest.skip("AMPLICON16S_CONFIG_DATI_REALI non impostata o file assente")
    from concurrent.futures import ProcessPoolExecutor

    file = sorted(Path(carica(_CONFIG_REALE).io.fastq_dir).glob("*.fastq.gz"))
    assert len(file) >= 900
    with ProcessPoolExecutor(max_workers=len(os.sched_getaffinity(0))) as gruppo:
        esiti = list(gruppo.map(conta_coppie, file))
    assert [s.nome for s in esiti if s.coppie_nello_stesso_file] == []
    assert sum(s.prime_di_coppia + s.seconde_di_coppia + s.nomi_ripetuti + s.nomi_oltre_due
               for s in esiti) == 0
    print(f"\n{len(file)} file, {sum(s.letture_esaminate for s in esiti)} letture, nessun segno di coppia")
