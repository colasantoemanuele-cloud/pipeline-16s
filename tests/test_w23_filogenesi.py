r"""Suite di test della settimana 23: la filogenesi opzionale, fase S9.

1. Inquadramento nel Piano Operativo
------------------------------------
Settimana 23 (W23), Fase F6 (filogenesi opzionale S9 sulle varianti finali,
dopo S13 e prima di S14; albero nell'oggetto finale e negli export).

2. Moduli sorgente coperti
--------------------------
* ``src/amplicon16s/steps/s09_phylogeny.py``, ``R/09_phylogeny.R``,
  ``R/lib/albero.R``
* ``src/amplicon16s/runner/graph.py`` (S9 dopo S13, S14 dipende da S9 solo se
  attiva, S10 non ne dipende)
* ``src/amplicon16s/steps/s10_phyloseq.py``, ``src/amplicon16s/steps/s14_finale.py``,
  ``R/14_finale.R``, ``R/lib/export.R``
* ``src/amplicon16s/config/schema.py``, ``src/amplicon16s/config/defaults.py``
  (``phylo.aligner``, ``phylo.model``)

3. Cosa valuta questo file
--------------------------
- il grafo: S9 dipende da S13 e si esegue fra S13 e S14; S14 dipende da S9 solo
  con ``phylo.enabled`` vero; S10 non dipende da S9; la precedenza S12-S13
  resta;
- la valutazione dello stato: con la filogenesi disattivata S9 e' disattivata,
  non incompleta, e la catena e' completa; attivarla dopo rifa' solo S9 e S14,
  disattivarla di nuovo rifa' solo S14;
- i parametri: ``phylo.aligner`` e ``phylo.model`` a insieme chiuso; il gruppo
  ``phylo`` e il seme cambiano l'impronta di S9, ``phylo.enabled`` quella di
  S14, nessuno quella di S10 e delle fasi a monte;
- ``E-S9-01``: con ``phylo.max_seqs`` sotto il numero di varianti finali la fase
  si ferma prima di avviare il calcolo, a revisione umana, senza nuovi tentativi;
- ``E-S9-02``: con meno di quattro varianti finali la fase si ferma allo stesso
  modo, prima di avviare il calcolo;
- sul sottoinsieme di prova, con la filogenesi attiva: l'albero e' radicato, le
  sue foglie sono gli identificativi delle varianti finali, l'oggetto finale lo
  contiene senza cambiare l'ordine delle varianti, l'export Newick coincide con
  quello di S9 ed e' nel manifesto dei checksum;
- l'albero, l'allineamento, il riepilogo e l'oggetto finale sono identici byte
  per byte fra due esecuzioni e con un numero di thread diverso;
- il testo Newick ha lunghezze a cifre fisse, e riletto e riscritto resta
  uguale; un albero con foglie diverse dalle varianti, o non radicato, e'
  respinto;
- con la filogenesi disattivata l'oggetto finale non ha l'albero;
- sul dataset completo: S9 e S14 sulla catena, con l'albero delle varianti
  finali.

4. Comandi Bash e scenari di esecuzione
---------------------------------------
    ``<immagine>`` e' l'immagine del container della pipeline; quella corrente
    e' indicata in ``test.txt``, sezione 1.3.

    1. Modalità locale standard (R di base con jsonlite, senza Bioconductor né
       dati reali):
       pytest tests/test_w23_filogenesi.py -v

    2. Modalità container Docker standard (sottoinsieme ridotto con
       R/Bioconductor):
       docker run --rm \
         -e PYTHONPATH=/app/src \
         -e AMPLICON16S_R_DIR=/app/R \
         -v "$(pwd)":/app \
         -w /app \
         <immagine> \
         pytest -o cache_dir=/tmp/.pytest_cache tests/test_w23_filogenesi.py -v

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
         pytest -o cache_dir=/tmp/.pytest_cache tests/test_w23_filogenesi.py -v

5. Risultato atteso
-------------------
Vedi ``test.txt``, scheda W23.

6. Razionale scientifico e sistemistico
---------------------------------------
- L'albero utile alle analisi e' quello delle varianti consegnate: costruito
  prima dei filtri conterrebbe contaminanti e varianti rare, e sul dataset di
  riferimento supererebbe il limite di fattibilita'.
- Le misure di diversita' filogenetica richiedono un albero radicato con le
  stesse varianti dell'oggetto: una foglia in piu' o in meno verrebbe scartata
  in silenzio dalla libreria.
- Una ricerca dell'albero con componenti casuali darebbe, a parita' di dati e
  di configurazione, alberi diversi: qui e' deterministica, e i test lo provano.
"""

from __future__ import annotations

import hashlib
import json
import os
import re
import shutil
import subprocess
from pathlib import Path
from typing import Any, ClassVar

import pytest
from conftest import NEGATIVO, POSITIVO, Campione, copia_esecuzione, crea_scenario
from sottoinsieme import config_ridotta, motivo_pacchetti_r_assenti, processori_disponibili

import amplicon16s.steps.s09_phylogeny as s09
from amplicon16s.config.resolve import risolvi
from amplicon16s.config.schema import Config, ErroreConfigurazione, valida
from amplicon16s.errors.catalog import Categoria, voce
from amplicon16s.errors.exceptions import ErrorePipeline
from amplicon16s.io_layer.artifacts import AlberoOutput, Fase
from amplicon16s.logging.logger import chiudi, ottieni
from amplicon16s.rbridge.payload import PREFISSO
from amplicon16s.rbridge.runner import cartella_r, trova_rscript
from amplicon16s.runner.executor import Conclusione, Esecutore
from amplicon16s.runner.graph import GRAFO, PRECEDENZE_OBBLIGATORIE, Passo
from amplicon16s.runner.project import ProjectRun, StatoPasso, passi_realizzati
from amplicon16s.runner.provenienza import leggi_registro
from amplicon16s.steps.base import PipelineStep, Produzione, StepContext
from amplicon16s.steps.s00_validate import ValidazioneIngressi
from amplicon16s.steps.s09_phylogeny import NOME_ALBERO, NOME_ALLINEAMENTO, NOME_FILOGENESI, Filogenesi
from amplicon16s.steps.s13_filtri import NOME_RIEPILOGO

ORDINE = (
    Passo.S0, Passo.S1, Passo.S2, Passo.S3, Passo.S4, Passo.S5, Passo.S6, Passo.S7,
    Passo.S8, Passo.S10, Passo.S11, Passo.S12, Passo.S13, Passo.S9, Passo.S14,
)


@pytest.fixture(autouse=True)
def uscite_pulite():
    """Chiude le uscite del log prima e dopo ogni test, perché nessun handler resti
    aperto sulla cartella temporanea.
    """
    chiudi()
    yield
    chiudi()


# --------------------------------------------------------------------------- #
# 1. Il grafo e la valutazione dello stato                                     #
# --------------------------------------------------------------------------- #


class Doppia(PipelineStep):
    """Una fase di prova: annota di essere stata eseguita e scrive un artefatto."""

    parametri: ClassVar[tuple[str, ...]] = ("run.seed",)

    def __init__(self, eseguite: list[Passo]) -> None:
        self.eseguite = eseguite

    def calcola(self, contesto: StepContext) -> Produzione:
        """Registra l'esecuzione e scrive l'artefatto della fase."""
        self.eseguite.append(self.passo)
        artefatto = contesto.albero.scrivi_testo(
            self.cartella, f"{str(self.passo).lower()}.txt", str(self.passo)
        )
        return Produzione((artefatto,))


def _fasi(eseguite: list[Passo]) -> dict[Passo, PipelineStep]:
    """S0 vera e una fase di prova per ogni altra, con i parametri della fase vera
    per S9 e S14: sono quelli che decidono che cosa rifa' un cambio di ``phylo``.
    """
    vere = passi_realizzati()
    fasi: dict[Passo, PipelineStep] = {Passo.S0: ValidazioneIngressi()}
    for passo in ORDINE[1:]:
        parametri = vere[passo].parametri if passo in (Passo.S9, Passo.S10, Passo.S14) else ("run.seed",)
        fasi[passo] = type(f"Doppia{passo}", (Doppia,), {"passo": passo, "parametri": parametri})(eseguite)
    return fasi


@pytest.fixture
def scenario(tmp_path):
    """Uno scenario su disco con quattro campioni e le loro letture."""
    return crea_scenario(tmp_path, [
        Campione("ERX3000001", "NOD1D4.L1"),
        Campione("ERX3000002", "NOD1D4.L2"),
        Campione("ERX3000003", "POS.P1.1", materiale=POSITIVO, posizione="Not Applicable"),
        Campione("ERX3000004", "BLANK.P1.1", materiale=NEGATIVO, posizione="Not Applicable"),
    ], con_letture=True)


def _con(config: Config, **gruppi: dict[str, Any]) -> Config:
    """La configurazione con i parametri indicati per gruppo sostituiti."""
    dati = config.model_dump(mode="python")
    for gruppo, valori in gruppi.items():
        dati[gruppo].update(valori)
    return valida(dati)


def _esegui(config: Config, eseguite: list[Passo]) -> tuple[ProjectRun, Any]:
    """Esegue la catena con le fasi di prova; ``eseguite`` riparte da vuoto."""
    eseguite.clear()
    run = ProjectRun(config, passi=_fasi(eseguite))
    return run, Esecutore(run).esegui()


def test_s9_sta_fra_s13_e_s14_e_s10_non_ne_dipende():
    """
    **Obiettivo**: Verificare che nel grafo S9 dipenda dalla sola S13, che
    l'ordine di esecuzione la collochi fra S13 e S14, che S10 dipenda da S0, S7
    e S8 e non da S9, che S9 resti l'unica fase facoltativa e che la precedenza
    obbligatoria sia ancora la sola S12-S13.

    **Razionale scientifico e sistemistico**: L'albero si costruisce sulle
    varianti che restano dopo i filtri, quelle consegnate; l'oggetto integrato
    di S10, con tutte le varianti, non lo porta.
    """
    assert GRAFO.ordine() == ORDINE
    assert GRAFO.nodo(Passo.S9).dipendenze == (Passo.S13,)
    assert GRAFO.nodo(Passo.S10).dipendenze == (Passo.S0, Passo.S7, Passo.S8)
    assert GRAFO.nodo(Passo.S14).dipendenze == (Passo.S10, Passo.S13, Passo.S9)
    assert [n.passo for n in GRAFO if n.facoltativa] == [Passo.S9]
    assert GRAFO.nodo(Passo.S9).parametro_attivazione == "phylo.enabled"
    assert GRAFO.nodo(Passo.S9).cartella is Fase.PHYLOGENY
    assert GRAFO.discendenti(Passo.S9) == (Passo.S14,)
    assert PRECEDENZE_OBBLIGATORIE == ((Passo.S12, Passo.S13, "E-S13-01"),)
    assert set(passi_realizzati()) == set(Passo)


def test_s14_dipende_da_s9_solo_con_la_filogenesi_attiva(scenario):
    """
    **Obiettivo**: Verificare che le dipendenze attive di S14 siano la sola S13
    con ``phylo.enabled`` falso e S13 e S9 con ``phylo.enabled`` vero, e che
    gli antenati di S14 comprendano S9 solo nel secondo caso.

    **Razionale scientifico e sistemistico**: Una fase disattivata non e'
    lavoro mancante: S14 non deve restare in attesa di un albero che la
    configurazione non chiede.
    """
    spenta, accesa = scenario.config, _con(scenario.config, phylo={"enabled": True})
    assert GRAFO.dipendenze_attive(Passo.S14, spenta) == (Passo.S10, Passo.S13)
    assert GRAFO.dipendenze_attive(Passo.S14, accesa) == (Passo.S10, Passo.S13, Passo.S9)
    assert Passo.S9 not in GRAFO.antenati(Passo.S14, spenta)
    assert Passo.S9 in GRAFO.antenati(Passo.S14, accesa)
    assert Passo.S9 not in GRAFO.attive(spenta) and Passo.S9 in GRAFO.attive(accesa)


def test_con_la_filogenesi_disattivata_s9_e_disattivata_e_la_catena_completa(scenario):
    """
    **Obiettivo**: Verificare che con ``phylo.enabled`` falso la catena si
    concluda senza eseguire S9, che la valutazione dia S9 disattivata con il
    parametro nel motivo e non fra le fasi da eseguire, e che il manifesto di
    S14 registri a monte la sola S13.

    **Razionale scientifico e sistemistico**: Nella nuova posizione S9 sta fra
    due fasi sempre attive: disattivata non deve fermare S14 ne' far risultare
    incompleta l'esecuzione.
    """
    eseguite: list[Passo] = []
    run, esito = _esegui(scenario.config, eseguite)
    assert esito.conclusione is Conclusione.COMPLETATA
    assert eseguite == [p for p in ORDINE[1:] if p is not Passo.S9]
    valutazione = run.valuta()
    assert valutazione.completa and valutazione.disattivate == (Passo.S9,)
    assert valutazione.situazioni[Passo.S9].stato is StatoPasso.DISATTIVATA
    assert "phylo.enabled" in valutazione.situazioni[Passo.S9].motivo
    assert Passo.S9 not in valutazione.da_eseguire
    manifesto = run.albero.manifesto_passo(Passo.S14, Fase.FINAL)
    assert list(manifesto.calcolata_su["a_monte"]) == ["S10", "S13"]
    assert not run.albero.cartella(Fase.PHYLOGENY).exists()


def test_attivare_la_filogenesi_rifa_solo_s9_e_s14(scenario):
    """
    **Obiettivo**: Verificare che, a catena conclusa, attivare la filogenesi
    faccia eseguire S9 e poi S14 lasciando concluse S10-S13; che cambiare
    ``phylo.max_seqs`` rifaccia di nuovo solo S9 e S14; e che disattivarla
    rifaccia la sola S14, con S9 disattivata.

    **Razionale scientifico e sistemistico**: Chiedere l'albero a esecuzione
    conclusa non deve ricalcolare l'oggetto integrato, i controlli, la
    decontaminazione e i filtri: l'albero non li tocca.
    """
    eseguite: list[Passo] = []
    _esegui(scenario.config, eseguite)

    accesa = _con(scenario.config, phylo={"enabled": True})
    run, esito = _esegui(accesa, eseguite)
    assert esito.conclusione is Conclusione.COMPLETATA
    assert eseguite == [Passo.S9, Passo.S14]
    assert set(run.albero.manifesto_passo(Passo.S14, Fase.FINAL).calcolata_su["a_monte"]) == {"S9", "S10", "S13"}

    _, esito = _esegui(_con(accesa, phylo={"max_seqs": 4000}), eseguite)
    assert eseguite == [Passo.S9, Passo.S14]

    run, esito = _esegui(scenario.config, eseguite)
    assert eseguite == [Passo.S14]
    assert run.valuta().situazioni[Passo.S9].stato is StatoPasso.DISATTIVATA
    assert run.valuta().completa


def test_i_parametri_della_filogenesi(tmp_path):
    """
    **Obiettivo**: Verificare i valori predefiniti del gruppo ``phylo``
    (disattivata, 5.000 varianti, ``decipher``, ``GTR+G+I``), che un allineatore
    o un modello non realizzati siano respinti dallo schema, e che il gruppo
    ``phylo`` e ``run.seed`` cambino l'impronta di S9, ``phylo.enabled`` anche
    quella di S14, e nessuno di questi quella di S10.

    **Razionale scientifico e sistemistico**: Un nome di modello accettato e
    non realizzato darebbe un albero calcolato con un altro modello; e un
    parametro dimenticato nella dichiarazione non rifarebbe la fase.
    """
    config = config_ridotta(tmp_path)
    assert (config.phylo.enabled, config.phylo.max_seqs, config.phylo.aligner, config.phylo.model) == (
        False, 5000, "decipher", "GTR+G+I",
    )
    for chiave, valore in (("aligner", "mafft"), ("model", "JC")):
        with pytest.raises(ErroreConfigurazione, match=f"phylo.{chiave}"):
            config_ridotta(tmp_path, phylo={chiave: valore})

    passi = passi_realizzati()
    assert passi[Passo.S9].parametri == ("phylo", "run.seed")
    assert "phylo.enabled" in passi[Passo.S14].parametri
    assert not any(p.startswith("phylo") for p in passi[Passo.S10].parametri)
    prima = risolvi(config)
    for variazione, cambiate in (
        ({"phylo": {"enabled": True}}, {Passo.S9, Passo.S14}),
        ({"phylo": {"max_seqs": 100}}, {Passo.S9}),
    ):
        dopo = risolvi(config_ridotta(tmp_path, **variazione))
        diverse = {p for p, f in passi.items() if f.calcolata_su(prima, {}) != f.calcolata_su(dopo, {})}
        assert diverse == cambiate, variazione
    seme = risolvi(config_ridotta(tmp_path, run={"seed": 7}))
    assert passi[Passo.S9].calcolata_su(prima, {}) != passi[Passo.S9].calcolata_su(seme, {})


def test_e_s9_01_ferma_prima_di_avviare_il_calcolo(tmp_path, monkeypatch):
    """
    **Obiettivo**: Verificare che con piu' varianti finali di ``phylo.max_seqs``
    S9 sollevi ``E-S9-01`` senza avviare il processo R, che con un numero di
    varianti pari al limite lo avvii, e che il codice sia a revisione umana e
    fuori dai codici ammessi al nuovo tentativo.

    **Razionale scientifico e sistemistico**: Il costo della massima
    verosimiglianza cresce rapidamente con il numero di sequenze: la guardia
    deve scattare prima di allocare calcolo, e un albero richiesto non si
    salta in silenzio.
    """
    config = config_ridotta(tmp_path, phylo={"enabled": True, "max_seqs": 40})
    albero = AlberoOutput(config.io.out_root)
    avviati: list[dict[str, Any]] = []

    def finto(script, parametri, *argomenti, **opzioni):
        avviati.append(parametri)
        raise RuntimeError("processo R avviato")

    monkeypatch.setattr(s09, "esegui_script", finto)
    contesto = StepContext(
        risolta=risolvi(config), albero=albero, logger=ottieni("prova"), inventario=None,
    )

    def con_varianti(numero: int) -> None:
        albero.prepara(Fase.FINAL_INTERMEDI)
        (albero.cartella(Fase.FINAL_INTERMEDI) / NOME_RIEPILOGO).write_text(
            json.dumps({"varianti": {"finali": numero}}), encoding="utf-8"
        )

    con_varianti(41)
    with pytest.raises(ErrorePipeline) as fermata:
        Filogenesi().calcola(contesto)
    assert fermata.value.codice == "E-S9-01"
    assert "41 varianti finali" in fermata.value.dettaglio and "(40)" in fermata.value.dettaglio
    assert avviati == []

    con_varianti(40)
    with pytest.raises(RuntimeError, match="processo R avviato"):
        Filogenesi().calcola(contesto)
    assert avviati[0]["modello"] == "GTR+G+I" and avviati[0]["seme"] == config.run.seed

    assert voce("E-S9-01").categoria is Categoria.REVISIONE_UMANA
    assert not voce("E-S9-01").ammette_retry
    assert "E-S9-01" not in config.retry.whitelist and Filogenesi.aggiustamenti == {}


def test_e_s9_02_ferma_con_troppo_poche_varianti(tmp_path, monkeypatch):
    """
    **Obiettivo**: Verificare che con meno di quattro varianti finali S9
    sollevi ``E-S9-02`` senza avviare il processo R, anche con zero varianti;
    che con quattro lo avvii; e che il codice sia a revisione umana, fuori dai
    codici ammessi al nuovo tentativo, con un messaggio che indica dove
    guardare e le due vie d'uscita.

    **Razionale scientifico e sistemistico**: Con tre foglie l'albero senza
    radice e' uno solo: non c'e' topologia da stimare. E' un esito prevedibile
    di filtri severi su un dataset piccolo, e va detto all'operatore con un
    codice del catalogo, non con un errore non catalogato del processo R.
    """
    config = config_ridotta(tmp_path, phylo={"enabled": True})
    albero = AlberoOutput(config.io.out_root)
    avviati: list[dict[str, Any]] = []

    def finto(script, parametri, *argomenti, **opzioni):
        avviati.append(parametri)
        raise RuntimeError("processo R avviato")

    monkeypatch.setattr(s09, "esegui_script", finto)
    contesto = StepContext(
        risolta=risolvi(config), albero=albero, logger=ottieni("prova"), inventario=None,
    )
    albero.prepara(Fase.FINAL_INTERMEDI)
    for numero in (0, 3):
        (albero.cartella(Fase.FINAL_INTERMEDI) / NOME_RIEPILOGO).write_text(
            json.dumps({"varianti": {"finali": numero}}), encoding="utf-8"
        )
        with pytest.raises(ErrorePipeline) as fermata:
            Filogenesi().calcola(contesto)
        assert fermata.value.codice == "E-S9-02"
        assert f"{numero} varianti finali, meno delle 4" in fermata.value.dettaglio
    assert avviati == []

    (albero.cartella(Fase.FINAL_INTERMEDI) / NOME_RIEPILOGO).write_text(
        json.dumps({"varianti": {"finali": s09.VARIANTI_MINIME}}), encoding="utf-8"
    )
    with pytest.raises(RuntimeError, match="processo R avviato"):
        Filogenesi().calcola(contesto)

    v = voce("E-S9-02")
    assert v.fase == "S9" and v.categoria is Categoria.REVISIONE_UMANA and not v.ammette_retry
    assert "E-S9-02" not in config.retry.whitelist
    assert "filtri_riepilogo.json" in v.azione and "phylo.enabled: false" in v.azione


# --------------------------------------------------------------------------- #
# 2. S9 e S14 sul sottoinsieme di prova                                        #
# --------------------------------------------------------------------------- #


_MOTIVO_BIOC = motivo_pacchetti_r_assenti(
    "dada2", "ggplot2", "ShortRead", "jsonlite", "phyloseq", "Biostrings", "decontam",
    "DECIPHER", "phangorn", "ape",
)


@pytest.fixture
def bioc():
    """Richiede R con i pacchetti della catena, DECIPHER e phangorn: salta senza,
    ma in CI fallisce.
    """
    if _MOTIVO_BIOC is not None:
        if os.environ.get("AMPLICON16S_RICHIEDI_BIOC") == "1":
            pytest.fail(f"Bioconductor e' richiesto in questo ambiente: {_MOTIVO_BIOC}")
        pytest.skip(_MOTIVO_BIOC)


def _con_filogenesi(base, cartella: Path, **sovrascrivi: dict[str, Any]) -> tuple[ProjectRun, Any]:
    """Una copia della catena conclusa, ripresa con la filogenesi attiva."""
    phylo = {"enabled": True, **sovrascrivi.pop("phylo", {})}
    run = copia_esecuzione(base, cartella, phylo=phylo, **sovrascrivi)
    return run, Esecutore(run).esegui()


@pytest.fixture
def con_albero(bioc, finale_calcolata, tmp_path):
    """La catena completa sul sottoinsieme di prova, con S9 e S14 eseguite a
    filogenesi attiva su una copia.
    """
    return _con_filogenesi(finale_calcolata, tmp_path / "attiva")


def _r(codice: str, cartella: Path) -> Any:
    """Esegue codice R che scrive un oggetto JSON in ``uscita``; le funzioni di
    ``R/lib/albero.R`` sono disponibili.
    """
    cartella.mkdir(parents=True, exist_ok=True)
    uscita = cartella / "uscita.json"
    script = cartella / "codice.R"
    script.write_text(
        f"uscita <- {json.dumps(str(uscita))}\n"
        f"source({json.dumps(str(cartella_r() / 'lib' / 'albero.R'))})\n{codice}\n",
        encoding="utf-8",
    )
    subprocess.run([str(trova_rscript()), "--vanilla", str(script)], check=True, capture_output=True)
    return json.loads(uscita.read_text(encoding="utf-8"))


def _impronte(cartella: Path) -> dict[str, str]:
    """L'impronta di ogni file della cartella, esclusi i file del ponte e i manifesti."""
    return {
        p.name: hashlib.sha256(p.read_bytes()).hexdigest()
        for p in sorted(cartella.iterdir())
        if p.is_file() and not p.name.startswith((PREFISSO, "manifest"))
    }


def _foglie(newick: str) -> list[str]:
    """Le foglie di un albero in formato Newick, nell'ordine in cui compaiono."""
    return re.findall(r"[(,]([A-Za-z0-9_.]+):", newick)


def test_l_albero_delle_varianti_finali_entra_nell_oggetto_finale(
    con_albero, finale_calcolata, tmp_path
):
    """
    **Obiettivo**: Verificare che attivando la filogenesi sulla catena conclusa
    si eseguano S9 e S14; che le foglie dell'albero di S9 siano esattamente gli
    identificativi delle varianti finali; che l'oggetto finale contenga
    l'albero, radicato, con le foglie nell'ordine delle varianti; che la
    provenienza di S9 e S14 sia quella del registro del sorgente; che l'ordine
    delle varianti e i conteggi siano quelli dell'oggetto senza albero; e che
    l'export Newick sia lo stesso testo scritto da S9, elencato nel manifesto
    dei checksum.

    **Razionale scientifico e sistemistico**: Le misure di diversita'
    filogenetica usano l'albero dell'oggetto: deve avere le stesse varianti, ne'
    una di piu' ne' una di meno, e non deve spostare le altre componenti.
    """
    run, esito = con_albero
    assert esito.conclusione is Conclusione.COMPLETATA
    assert [r.passo for r in esito.eseguite] == [Passo.S9, Passo.S14]
    filogenesi, finale = run.albero.cartella(Fase.PHYLOGENY), run.albero.cartella(Fase.FINAL)
    manifesto = run.albero.manifesto_passo(Passo.S9, Fase.PHYLOGENY)
    assert set(manifesto.nomi) == {NOME_ALBERO, NOME_ALLINEAMENTO, NOME_FILOGENESI}
    # Il sorgente eseguito e' quello del registro, come per le altre fasi.
    registro = leggi_registro()
    for passo, fase in ((Passo.S9, Fase.PHYLOGENY), (Passo.S14, Fase.FINAL)):
        provenienza = run.albero.manifesto_passo(passo, fase).provenienza
        assert provenienza["sorgente"] == registro[str(passo)]["sorgente"], passo
        assert provenienza["versione"] == registro[str(passo)]["versione"], passo

    letti = _r(f"""
suppressPackageStartupMessages(library(phyloseq))
f <- readRDS({json.dumps(str(finale / 'ps_final.rds'))})
a <- phy_tree(f)
jsonlite::write_json(list(
  varianti = taxa_names(f), foglie = a$tip.label, radicato = ape::is.rooted(a),
  newick = testo_newick(a), righe = rownames(as(otu_table(f), "matrix")),
  rami_negativi = sum(a$edge.length < 0)
), uscita, auto_unbox = TRUE)
""", tmp_path / "r")
    newick = (filogenesi / NOME_ALBERO).read_text(encoding="utf-8")
    assert sorted(_foglie(newick)) == sorted(letti["varianti"])
    assert len(set(letti["varianti"])) == len(letti["varianti"]) >= 4
    assert letti["foglie"] == letti["varianti"] == letti["righe"]
    assert letti["radicato"] is True and letti["rami_negativi"] == 0
    assert letti["newick"] + "\n" == newick
    assert (finale / NOME_ALBERO).read_bytes() == (filogenesi / NOME_ALBERO).read_bytes()

    senza = finale_calcolata[0].albero.cartella(Fase.FINAL)
    for nome in ("conteggi.tsv", "tassonomia.tsv", "metadati.tsv", "sequenze.fasta"):
        assert (finale / nome).read_bytes() == (senza / nome).read_bytes(), nome
    somme = (finale / "checksum.sha256").read_text(encoding="utf-8")
    assert f"{hashlib.sha256((finale / NOME_ALBERO).read_bytes()).hexdigest()}  {NOME_ALBERO}\n" in somme
    assert run.albero.manifesto_passo(Passo.S14, Fase.FINAL).metriche["albero"] is True

    riepilogo = json.loads((filogenesi / NOME_FILOGENESI).read_text(encoding="utf-8"))
    assert riepilogo["varianti"] == len(letti["varianti"]) == manifesto.metriche["varianti"]
    assert riepilogo["modello"] == "GTR+G+I" and riepilogo["radicato"] is True
    assert riepilogo["seme"] == run.config.run.seed
    assert "punto medio" in riepilogo["radicamento"] and "NNI" in riepilogo["ricerca"]
    allineate = (filogenesi / NOME_ALLINEAMENTO).read_text(encoding="utf-8").splitlines()
    assert [r[1:] for r in allineate[0::2]] == letti["varianti"]
    assert {len(r) for r in allineate[1::2]} == {riepilogo["colonne_allineamento"]}
    print(f"\nS9: {dict(esito.eseguite[0].metriche)}")


def test_l_albero_e_identico_fra_due_esecuzioni_e_con_thread_diversi(
    con_albero, finale_calcolata, tmp_path
):
    """
    **Obiettivo**: Verificare che l'albero, l'allineamento e il riepilogo di S9,
    e l'oggetto finale e gli export di S14, siano identici byte per byte in una
    seconda esecuzione e in una terza con un numero di thread diverso (con lo
    stesso numero su una macchina con un solo processore, dove un altro valore
    non e' ammesso).

    **Razionale scientifico e sistemistico**: La ricerca dell'albero e'
    configurata senza componenti casuali, e il numero di processi riguarda il
    solo allineamento: lo stesso dato e la stessa configurazione danno lo
    stesso albero su qualunque macchina.
    """
    run, _ = con_albero
    attese = {f: _impronte(run.albero.cartella(f)) for f in (Fase.PHYLOGENY, Fase.FINAL)}
    assert NOME_ALBERO in attese[Fase.PHYLOGENY] and "ps_final.rds" in attese[Fase.FINAL]
    thread = run.config.run.threads
    altri = 1 if thread > 1 else min(2, processori_disponibili())
    for nome, processi in (("seconda", thread), ("thread", altri)):
        altra, esito = _con_filogenesi(finale_calcolata, tmp_path / nome, run={"threads": processi})
        assert esito.conclusione is Conclusione.COMPLETATA
        for fase, impronte in attese.items():
            assert _impronte(altra.albero.cartella(fase)) == impronte, (nome, fase)


def test_e_s9_01_ferma_la_catena_prima_del_calcolo(bioc, finale_calcolata, tmp_path):
    """
    **Obiettivo**: Verificare che con ``phylo.max_seqs`` sotto il numero di
    varianti finali l'esecuzione si arresti in S9 con ``E-S9-01`` al primo
    tentativo, senza che la cartella della filogenesi contenga la richiesta di
    un processo R, e che S14 non venga eseguita.

    **Razionale scientifico e sistemistico**: Un albero richiesto e non
    costruibile non si salta: l'oggetto finale senza la componente domandata
    sarebbe consegnato come completo.
    """
    run, esito = _con_filogenesi(finale_calcolata, tmp_path / "limite", phylo={"max_seqs": 1})
    assert esito.conclusione is Conclusione.ARRESTATA
    assert (esito.punto.passo, esito.punto.codice, esito.punto.tentativi) == (Passo.S9, "E-S9-01", 1)
    assert esito.eseguite == ()
    cartella = run.albero.cartella(Fase.PHYLOGENY)
    assert not cartella.exists() or not any(cartella.iterdir())
    situazioni = run.valuta().situazioni
    assert situazioni[Passo.S9].stato is StatoPasso.DA_ESEGUIRE
    assert situazioni[Passo.S14].motivo == "a monte da eseguire: S9"


def test_senza_filogenesi_l_oggetto_finale_non_ha_l_albero(bioc, finale_calcolata, tmp_path):
    """
    **Obiettivo**: Verificare che con ``phylo.enabled`` falso la catena completa
    dia S9 disattivata, nessuna cartella della filogenesi, nessun ``albero.nwk``
    fra i file consegnati e un oggetto finale privo dell'albero.

    **Razionale scientifico e sistemistico**: E' la configurazione predefinita:
    lo spostamento di S9 non deve aver cambiato l'oggetto che si consegna senza
    filogenesi.
    """
    run, _ = finale_calcolata
    valutazione = run.valuta()
    assert valutazione.completa and valutazione.situazioni[Passo.S9].stato is StatoPasso.DISATTIVATA
    finale = run.albero.cartella(Fase.FINAL)
    assert not (finale / NOME_ALBERO).exists()
    assert NOME_ALBERO not in (finale / "checksum.sha256").read_text(encoding="utf-8")
    assert run.albero.manifesto_passo(Passo.S14, Fase.FINAL).metriche["albero"] is False
    letto = _r(f"""
suppressPackageStartupMessages(library(phyloseq))
f <- readRDS({json.dumps(str(finale / 'ps_final.rds'))})
jsonlite::write_json(list(albero = !is.null(phy_tree(f, errorIfNULL = FALSE))), uscita, auto_unbox = TRUE)
""", tmp_path / "r")
    assert letto["albero"] is False


def test_il_newick_ha_cifre_fisse_e_le_foglie_si_verificano(bioc, tmp_path):
    """
    **Obiettivo**: Verificare che ``testo_newick`` scriva le lunghezze dei rami
    con otto cifre significative e il punto come separatore anche con
    impostazioni locali che usano la virgola, che rileggere e riscrivere dia lo
    stesso testo, che numerare le foglie in un altro ordine non cambi il testo,
    e che ``difetti_albero`` segnali una foglia mancante, una estranea e un
    albero non radicato.

    **Razionale scientifico e sistemistico**: Il file Newick e' la forma
    consegnata dell'albero: deve essere lo stesso su ogni macchina, e un albero
    con varianti diverse da quelle dell'oggetto verrebbe sfoltito in silenzio.
    """
    letto = _r("""
options(OutDec = ",")
a <- ape::read.tree(text = "((ASV1:0.123456789012,ASV2:1e-09):0.25,(ASV3:2,ASV7:0.333333333333):0.5);")
testo <- testo_newick(a)
riordinato <- foglie_in_ordine(a, c("ASV7", "ASV3", "ASV2", "ASV1"))
senza_radice <- ape::unroot(a)
jsonlite::write_json(list(
  testo = testo,
  riletto = testo_newick(ape::read.tree(text = testo)),
  riordinato = testo_newick(riordinato), foglie = riordinato$tip.label,
  giusto = difetti_albero(a, c("ASV1", "ASV2", "ASV3", "ASV7")),
  mancante = difetti_albero(a, c("ASV1", "ASV2", "ASV3", "ASV7", "ASV9")),
  estranea = difetti_albero(a, c("ASV1", "ASV2", "ASV3")),
  non_radicato = difetti_albero(senza_radice, c("ASV1", "ASV2", "ASV3", "ASV7"))
), uscita, auto_unbox = FALSE)
""", tmp_path)
    assert letto["testo"] == ["((ASV1:0.12345679,ASV2:1e-09):0.25,(ASV3:2,ASV7:0.33333333):0.5);"]
    assert letto["riletto"] == letto["testo"] == letto["riordinato"]
    assert letto["foglie"] == ["ASV7", "ASV3", "ASV2", "ASV1"]
    assert letto["giusto"] == []
    assert "1 varianti senza foglia, 0 foglie che non sono varianti" in letto["mancante"][0]
    assert "0 varianti senza foglia, 1 foglie che non sono varianti" in letto["estranea"][0]
    assert letto["non_radicato"] == ["albero non radicato"]


def test_un_albero_alterato_ferma_s14_con_e_s14_01(con_albero, tmp_path):
    """
    **Obiettivo**: Verificare che, alterato ``albero.nwk`` di S9 dopo la sua
    conclusione, S14 eseguita da sola si fermi con ``E-S14-01`` perche' il file
    non corrisponde al manifesto di S9.

    **Razionale scientifico e sistemistico**: L'albero che entra nell'oggetto
    finale e' quello che S9 ha calcolato e registrato, non un file trovato nella
    cartella.
    """
    run, _ = con_albero
    percorso = run.albero.cartella(Fase.PHYLOGENY) / NOME_ALBERO
    percorso.write_text(percorso.read_text(encoding="utf-8").replace("ASV", "VAR"), encoding="utf-8")
    with pytest.raises(ErrorePipeline) as fermata:
        run.fase(Passo.S14).calcola(run.contesto(Passo.S14).ristretto(run.fase(Passo.S14).parametri))
    assert fermata.value.codice == "E-S14-01" and NOME_ALBERO in fermata.value.dettaglio


# --------------------------------------------------------------------------- #
# 3. Sul dataset completo                                                      #
# --------------------------------------------------------------------------- #


@pytest.mark.dati_reali
def test_la_filogenesi_sul_dataset_completo(bioc, catena_reale, tmp_path):
    """
    **Obiettivo**: Verificare, su una copia per collegamenti fisici della catena
    del dataset completo, che attivando la filogenesi si eseguano le sole S9 e
    S14, che le foglie dell'albero siano le varianti finali di S13, e che la
    catena condivisa non venga modificata.

    **Razionale scientifico e sistemistico**: Sul dataset di riferimento le
    varianti finali sono 1.753, sotto ``phylo.max_seqs``: e' il caso d'uso per
    cui S9 e' stata collocata dopo i filtri. La copia per collegamenti fisici
    non occupa spazio, perche' ogni scrittura della pipeline sostituisce il file
    per rinomina.
    """
    base, _ = catena_reale
    origine = Path(base.config.io.out_root)
    prima = {str(p.relative_to(origine)): p.stat().st_mtime_ns for p in origine.rglob("*") if p.is_file()}
    shutil.copytree(origine, tmp_path / "out", copy_function=os.link)
    dati = base.config.model_dump(mode="python")
    dati["io"]["out_root"] = str(tmp_path / "out")
    dati["phylo"]["enabled"] = True
    run = ProjectRun(valida(dati))
    # Fino a S14, come la catena condivisa: S1 non e' fra i suoi antenati, e la
    # catena condivisa non la esegue.
    esito = Esecutore(run, fino_a=Passo.S14).esegui()
    assert esito.conclusione is Conclusione.COMPLETATA
    assert [r.passo for r in esito.eseguite] == [Passo.S9, Passo.S14]
    finali = json.loads((run.albero.cartella(Fase.FINAL_INTERMEDI) / NOME_RIEPILOGO).read_text())["varianti"]["finali"]
    newick = (run.albero.cartella(Fase.FINAL) / NOME_ALBERO).read_text(encoding="utf-8")
    assert len(set(_foglie(newick))) == finali == esito.eseguite[0].metriche["varianti"]
    dopo = {str(p.relative_to(origine)): p.stat().st_mtime_ns for p in origine.rglob("*") if p.is_file()}
    assert dopo == prima
    print(f"\nS9 sul dataset completo: {esito.eseguite[0].secondi} s, {finali} varianti")
