r"""Suite di test della settimana 17: codice nell'impronta, G12 come precondizione, fase S8.

1. Inquadramento nel Piano Operativo
------------------------------------
Settimana 17 (W17), Fase F5 (assegnazione tassonomica S8 in ``08_taxonomy/``,
preceduta dalla versione del calcolo nell'impronta di fase con il registro del
sorgente, e da G12 spostato fra le precondizioni).

2. Moduli sorgente coperti
--------------------------
* ``src/amplicon16s/runner/provenienza.py``, ``scripts/registro_sorgente.py``,
  ``src/amplicon16s/steps/registro_sorgente.json``
* ``src/amplicon16s/steps/base.py`` (``versione``, ``script_r``, provenienza)
* ``src/amplicon16s/io_layer/artifacts.py`` (provenienza nel manifesto)
* ``src/amplicon16s/runner/project.py``, ``src/amplicon16s/runner/executor.py``
* ``src/amplicon16s/steps/s00_validate.py``, ``src/amplicon16s/gates/g01_g15.py`` (G12)
* ``src/amplicon16s/steps/s08_taxonomy.py``, ``R/08_taxonomy.R``
* ``scripts/costruisci_riferimento_sintetico.py``, ``tests/fixtures/riferimento_sintetico/``
* ``container/dada2/`` (correzione dei pareggi di dada2), ``container/verify_renv_lock.R``,
  ``renv.lock``

3. Cosa valuta questo file
--------------------------
- il registro del sorgente corrisponde al codice, e ogni fase realizzata
  dichiara la propria versione; una modifica al sorgente non registrata fa
  fallire il controllo, e l'aggiornamento chiede una scelta esplicita fra
  incrementare la versione e dichiarare la modifica senza effetto;
- la provenienza registra l'impronta degli script R effettivamente eseguiti
  (la cartella indicata da ``AMPLICON16S_R_DIR``), non delle copie del
  repository: uno script diverso risulta nel manifesto e produce un avviso
  alla ripresa, senza ricalcolo; una versione incrementata invece rifa' la
  fase;
- S0 non dichiara il riferimento tassonomico: cambiarlo lascia concluse S0-S7 e
  rifa' S8; G12 respinge all'avvio un riferimento alterato, anche con S0
  conclusa;
- S8 sul sottoinsieme di prova, con il riferimento sintetico: tabella
  tassonomica e bootstrap per rango, marcatura dei taxa difettosi, copertura del
  phylum per classe, arresto con ``E-S8-02``, stessi byte fra due esecuzioni e
  con un numero di thread diverso;
- ``renv.lock`` dichiara dada2 1.36.0.1 con la correzione dei pareggi e il suo
  SHA-256, e nell'immagine l'R della pipeline carica quella versione;
- sul dataset completo, S8 dopo S0-S7 con SILVA 138: ai due livelli di
  diluizione piu' alti la variante dominante dei controlli positivi e'
  Variovorax in ogni piastra.

4. Comandi Bash e scenari di esecuzione
---------------------------------------
    ``<immagine>`` e' l'immagine del container della pipeline; quella corrente
    e' indicata in ``test.txt``, sezione 1.3.

    1. Modalità locale standard (R di base con jsonlite, senza Bioconductor né
       dati reali):
       pytest tests/test_w17_tassonomia.py -v

    2. Modalità container Docker standard (sottoinsieme ridotto con
       R/Bioconductor):
       docker run --rm \
         -e PYTHONPATH=/app/src \
         -e AMPLICON16S_R_DIR=/app/R \
         -v "$(pwd)":/app \
         -w /app \
         <immagine> \
         pytest -o cache_dir=/tmp/.pytest_cache tests/test_w17_tassonomia.py -v

    3. Modalità container Docker completa (con i 2.4 GB di dati reali OSD-734;
       la configurazione e i percorsi che contiene devono stare nella cartella
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
         pytest -o cache_dir=/tmp/.pytest_cache tests/test_w17_tassonomia.py -v

5. Risultato atteso
-------------------
21 test totali:
- 12 passed, 9 skipped in ambiente locale standard (~0.4s): i 9 test che eseguono
  le fasi richiedono dada2, ShortRead e ggplot2, uno anche i dati reali;
- 20 passed, 1 skipped nel container Docker standard sul sottoinsieme ridotto
  (~41s): resta saltato il test sui dati reali;
- 21 passed nel container Docker con i dati reali OSD-734 (~39s oltre la catena
  S0-S14 condivisa ``catena_reale``, calcolata una volta per sessione).

6. Razionale scientifico e sistemistico
---------------------------------------
- Senza il codice nell'impronta, correggere un difetto di una fase lascerebbe
  validi gli artefatti prodotti dal codice difettoso. La versione decide la
  validita' senza ricalcolare a ogni commento; il registro obbliga a
  scegliere; la provenienza sugli script eseguiti rende visibile un'immagine
  che esegue codice diverso da quello del repository.
- Il riferimento tassonomico e' un ingresso di S8, non di S0: cambiarlo non
  deve rifare tre quarti d'ora di inferenza.
- Il bootstrap di assignTaxonomy e' casuale: il seme viene da ``run.seed``, e
  gli stessi byte con un numero di thread diverso confermano che
  ``run.threads`` resta legittimamente fuori dall'impronta.
- La variante dominante dei controlli positivi, un solo ceppo noto, e' un
  riscontro con una verita' nota sull'intera catena, dal FASTQ alla tassonomia.
"""

from __future__ import annotations

import csv
import hashlib
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path
from typing import ClassVar

import pytest
from conftest import NEGATIVO, POSITIVO, Campione, copia_esecuzione, crea_scenario
from sottoinsieme import RIFERIMENTO, config_ridotta, motivo_pacchetti_r_assenti

from amplicon16s.config.schema import valida
from amplicon16s.errors.catalog import Categoria
from amplicon16s.errors.exceptions import ErrorePipeline
from amplicon16s.io_layer.artifacts import Fase
from amplicon16s.logging.logger import NOME_FILE_LOG, chiudi, configura
from amplicon16s.metadata.models import ClasseCampione
from amplicon16s.rbridge.payload import PREFISSO
from amplicon16s.rbridge.runner import VARIABILE_CARTELLA_R, trova_rscript
from amplicon16s.runner.executor import Conclusione, Esecutore
from amplicon16s.runner.graph import Passo
from amplicon16s.runner.project import ProjectRun, StatoPasso, passi_realizzati
from amplicon16s.runner.provenienza import (
    RegistroNonAggiornabile,
    aggiorna_registro,
    cartella_r_repository,
    differenze_registro,
    leggi_registro,
    provenienza,
    voce_registro,
)
from amplicon16s.steps.base import PipelineStep, Produzione, StepContext
from amplicon16s.steps.s00_validate import ValidazioneIngressi, impronta_dati_grezzi
from amplicon16s.steps.s04_dada import InferenzaVarianti
from amplicon16s.steps.s08_taxonomy import controlla_copertura, marca_difetti


@pytest.fixture(autouse=True)
def uscite_pulite():
    """Chiude le uscite del log prima e dopo ogni test, perché nessun handler resti
    aperto sulla cartella temporanea.
    """
    chiudi()
    yield
    chiudi()


_MOTIVO_ASSENTI = motivo_pacchetti_r_assenti("dada2", "ggplot2", "ShortRead", "jsonlite")
#: La radice del repository, o della sua copia nell'immagine, dalla posizione dei test.
RADICE = Path(__file__).resolve().parents[1]


@pytest.fixture
def dada2():
    """Richiede R con dada2: salta senza, ma in CI fallisce."""
    if _MOTIVO_ASSENTI is not None:
        if os.environ.get("AMPLICON16S_RICHIEDI_BIOC") == "1":
            pytest.fail(f"Bioconductor e' richiesto in questo ambiente: {_MOTIVO_ASSENTI}")
        pytest.skip(_MOTIVO_ASSENTI)


def _fasi() -> dict[str, PipelineStep]:
    """Le fasi realizzate, per nome."""
    return {str(p): f for p, f in passi_realizzati().items()}


def _copia_r(destinazione: Path) -> Path:
    """Una copia degli script R del repository, da poter modificare."""
    return Path(shutil.copytree(cartella_r_repository(), destinazione))


def _campioni() -> list[Campione]:
    """Quattro campioni sulla piastra 1: due biologici, un positivo, un negativo."""
    return [
        Campione("ERX3000001", "NOD1D4.L1", piastra="1"),
        Campione("ERX3000002", "NOD1D4.L2", piastra="1"),
        Campione("ERX3000003", "POS.P1.1", materiale=POSITIVO,
                 posizione="Not Applicable", piastra="1"),
        Campione("ERX3000004", "BLANK.P1.1", materiale=NEGATIVO,
                 posizione="Not Applicable", piastra="1"),
    ]


def _impronte(cartella: Path) -> dict[str, str]:
    """L'MD5 di ogni file della cartella, esclusi i file del ponte e i manifesti."""
    return {
        p.name: hashlib.md5(p.read_bytes()).hexdigest()
        for p in sorted(cartella.iterdir())
        if not p.name.startswith((PREFISSO, "manifest"))
    }


def _tsv(percorso: Path) -> list[dict[str, str]]:
    """Le righe di una tabella separata da tabulazioni."""
    with open(percorso, encoding="utf-8", newline="") as file:
        return list(csv.DictReader(file, delimiter="\t"))


# --------------------------------------------------------------------------- #
# 1. La versione e il registro del sorgente                                     #
# --------------------------------------------------------------------------- #


def test_il_registro_corrisponde_al_codice_e_ogni_fase_dichiara_la_versione():
    """
    **Obiettivo**: Verificare che il registro del sorgente corrisponda alle
    copie del repository per ogni fase realizzata, e che ogni fase dichiari la
    propria versione nella sua classe, non ereditandola.

    **Razionale scientifico e sistemistico**: E' il controllo che fallisce
    quando si modifica una fase senza aggiornare il registro: chi la modifica
    deve scegliere fra incrementare la versione e dichiarare la modifica senza
    effetto (``scripts/registro_sorgente.py``).
    """
    assert differenze_registro(leggi_registro(), _fasi()) == {}
    for nome, fase in _fasi().items():
        assert "versione" in type(fase).__dict__, nome
        assert leggi_registro()[nome]["versione"] == fase.versione


def test_da_un_altra_fase_si_importano_solo_nomi_e_i_moduli_comuni_entrano_nel_sorgente():
    """
    **Obiettivo**: Verificare che nessun modulo di fase importi da un'altra fase
    qualcosa che non sia una costante (un nome in maiuscolo, come il nome di un
    artefatto), e che i moduli del calcolo ricavati dagli import comprendano
    quelli condivisi: il lettore del tracciamento per S3-S7 e S13, i gate e il
    crosswalk per S0, il lettore delle tabelle dei metadati per S0 e S10.

    **Razionale scientifico e sistemistico**: Il sorgente di una fase comprende
    i moduli che importa, esclusi quelli delle altre fasi, da cui importa solo
    nomi di artefatti. Una funzione presa da un'altra fase sarebbe calcolo fuori
    dall'impronta: modificarla non cambierebbe il sorgente registrato della fase
    che la usa. Prima ``leggi_conteggi`` stava nel modulo di S2 e partecipava al
    calcolo di S5, S6 e S13 senza entrarne nel sorgente.
    """
    import ast

    from amplicon16s.runner.provenienza import MODULO_DI_FASE, moduli_del_calcolo

    for nome, fase in _fasi().items():
        modulo = type(fase).__module__
        sorgente = Path(sys.modules[modulo].__file__).read_text(encoding="utf-8")
        for nodo in ast.walk(ast.parse(sorgente)):
            if isinstance(nodo, ast.ImportFrom) and nodo.module and MODULO_DI_FASE.match(nodo.module):
                if nodo.module != modulo:
                    assert all(a.name.isupper() for a in nodo.names), (nome, nodo.module)
    conteggi = "amplicon16s.io_layer.conteggi"
    for passo in ("S3", "S4", "S5", "S6", "S7", "S13"):
        assert conteggi in moduli_del_calcolo(type(_fasi()[passo]).__module__), passo
    assert {"amplicon16s.gates.g01_g15", "amplicon16s.metadata.crosswalk",
            "amplicon16s.metadata.tabelle"} <= set(moduli_del_calcolo("amplicon16s.steps.s00_validate"))
    assert "amplicon16s.metadata.tabelle" in moduli_del_calcolo("amplicon16s.steps.s10_phyloseq")


def test_le_esclusioni_dal_sorgente_sono_esplicite_e_non_toccano_gli_artefatti():
    """
    **Obiettivo**: Verificare che ogni modulo escluso dal sorgente delle fasi
    sia un modulo esistente con la sua ragione scritta; che il lettore
    dell'inventario, che il contesto passa alle fasi senza import, entri nel
    sorgente di ogni fase che dipende da S0 (e di S0, che ne importa il nome
    del file) e solo di quelle; che il catalogo, escluso, entri soltanto in S0,
    che ne scrive la sintesi in gates.json; e che i moduli che
    scrivono, formattano o serializzano il contenuto di un artefatto non siano
    esclusi: per ogni fase con uno script R entrano nel sorgente il ponte
    (richiesta e ambiente di R), la scrittura degli artefatti e i checksum, e
    per S2-S5 il calcolo dei valori corretti di un nuovo tentativo.

    **Razionale scientifico e sistemistico**: Un modulo escluso puo' cambiare
    senza che cambi l'impronta del sorgente di alcuna fase: se cambiasse i byte
    di un artefatto, il registro non lo vedrebbe. L'elenco esplicito obbliga a
    dichiarare la ragione di ogni esclusione.
    """
    import importlib

    from amplicon16s.runner.provenienza import (
        DAL_CONTESTO,
        ECCEZIONI_PER_FASE,
        ESCLUSI,
        moduli_del_calcolo,
        moduli_della_fase,
    )

    for modulo, ragione in ESCLUSI.items():
        importlib.import_module(modulo)
        assert ragione.strip(), modulo
    for modulo, (produttore, ragione) in DAL_CONTESTO.items():
        importlib.import_module(modulo)
        assert ragione.strip() and modulo not in ESCLUSI, modulo
        for nome, fase in _fasi().items():
            dipende = produttore in {str(d) for d in fase.nodo.dipendenze}
            assert (modulo in moduli_della_fase(fase)) is (dipende or nome == produttore), (nome, modulo)
    for nome, eccezioni in ECCEZIONI_PER_FASE.items():
        for modulo, ragione in eccezioni.items():
            assert ragione.strip() and modulo in ESCLUSI
            assert modulo in moduli_della_fase(_fasi()[nome])
    # Il catalogo, escluso, entra solo in S0 (gates.json).
    for nome, fase in _fasi().items():
        assert ("amplicon16s.errors.catalog" in moduli_della_fase(fase)) is (nome == "S0"), nome
    calcolo = {"amplicon16s.io_layer.artifacts", "amplicon16s.io_layer.checksums",
               "amplicon16s.rbridge.payload", "amplicon16s.rbridge.runner",
               "amplicon16s.runner.retry"}
    assert not calcolo & set(ESCLUSI)
    for nome, fase in _fasi().items():
        moduli = set(moduli_del_calcolo(type(fase).__module__))
        assert {"amplicon16s.io_layer.artifacts", "amplicon16s.io_layer.checksums"} <= moduli, nome
        if fase.script_r is not None:
            assert {"amplicon16s.rbridge.payload", "amplicon16s.rbridge.runner"} <= moduli, nome
    for nome in ("S2", "S3", "S4", "S5"):
        assert "amplicon16s.runner.retry" in moduli_del_calcolo(type(_fasi()[nome]).__module__)


def test_una_modifica_al_sorgente_non_registrata_fa_fallire_il_controllo(tmp_path):
    """
    **Obiettivo**: Verificare che, modificando in una copia lo script R di S4,
    il confronto con il registro segnali S4 con il sorgente cambiato a parita'
    di versione, che l'aggiornamento sia rifiutato senza una scelta, e che sia
    accettato dichiarando la modifica senza effetto.

    **Razionale scientifico e sistemistico**: Un aggiornamento automatico del
    registro renderebbe il controllo inutile: la scelta deve restare scritta.
    """
    cartella = _copia_r(tmp_path / "R")
    with open(cartella / "04_dada.R", "a", encoding="utf-8") as file:
        file.write("# una modifica\n")
    registro = leggi_registro()
    assert differenze_registro(registro, _fasi(), cartella) == {
        "S4": "sorgente cambiato a parita' di versione"
    }
    with pytest.raises(RegistroNonAggiornabile, match="S4"):
        aggiorna_registro(registro, _fasi(), cartella_r=cartella)
    with pytest.raises(RegistroNonAggiornabile, match="non ne hanno bisogno: S5"):
        aggiorna_registro(registro, _fasi(), frozenset({"S4", "S5"}), cartella)
    nuovo = aggiorna_registro(registro, _fasi(), frozenset({"S4"}), cartella)
    assert nuovo["S4"]["sorgente"] != registro["S4"]["sorgente"]
    assert nuovo["S4"]["versione"] == registro["S4"]["versione"]


def test_una_versione_incrementata_si_registra_senza_dichiarazioni(tmp_path):
    """
    **Obiettivo**: Verificare che, con il sorgente cambiato e la versione
    incrementata, il confronto riporti la versione incrementata e
    l'aggiornamento non chieda altro; e che una versione scesa sia un errore.

    **Razionale scientifico e sistemistico**: E' la strada per le modifiche che
    cambiano cio' che la fase calcola: la ripresa la rifara'.
    """
    cartella = _copia_r(tmp_path / "R")
    with open(cartella / "04_dada.R", "a", encoding="utf-8") as file:
        file.write("# il calcolo cambia\n")

    class Successiva(InferenzaVarianti):
        versione: ClassVar[int] = InferenzaVarianti.versione + 1

    fasi = {**_fasi(), "S4": Successiva()}
    # La sottoclasse vive in questo modulo: il sorgente cambia comunque.
    assert differenze_registro(leggi_registro(), fasi, cartella) == {"S4": "versione incrementata"}
    assert (aggiorna_registro(leggi_registro(), fasi, cartella_r=cartella)["S4"]["versione"]
            == InferenzaVarianti.versione + 1)

    class Precedente(InferenzaVarianti):
        versione: ClassVar[int] = 0

    assert "errore" in differenze_registro(leggi_registro(), {**_fasi(), "S4": Precedente()})["S4"]


def test_lo_strumento_verifica_il_registro():
    """
    **Obiettivo**: Verificare che ``scripts/registro_sorgente.py`` in modalita'
    di verifica esca con 0 sul repository allineato e rifiuti
    ``--senza-effetto`` senza ``--aggiorna``.

    **Razionale scientifico e sistemistico**: Lo strumento e' il modo previsto
    di aggiornare il registro, invece di farlo a mano.
    """
    import importlib.util

    specifica = importlib.util.spec_from_file_location(
        "registro_sorgente", RADICE / "scripts" / "registro_sorgente.py"
    )
    modulo = importlib.util.module_from_spec(specifica)
    specifica.loader.exec_module(modulo)
    assert modulo.main([]) == 0
    with pytest.raises(SystemExit):
        modulo.main(["--senza-effetto", "S4"])


def test_la_provenienza_usa_gli_script_r_eseguiti(tmp_path, monkeypatch):
    """
    **Obiettivo**: Verificare che l'impronta del sorgente nella provenienza di
    S4 si calcoli nella cartella indicata da ``AMPLICON16S_R_DIR``: una copia
    identica degli script ha l'impronta del registro, una copia con
    ``04_dada.R`` modificato ha un'impronta diversa, e il file diverso e' quello.

    **Razionale scientifico e sistemistico**: Un'immagine che esegue script R
    diversi da quelli del repository (per esempio quelli copiati quando e'
    stata costruita) deve risultare tale nel manifesto, invece di passare
    inosservata.
    """
    config = config_ridotta(tmp_path)
    cartella = _copia_r(tmp_path / "R")
    monkeypatch.setenv(VARIABILE_CARTELLA_R, str(cartella))
    uguale = provenienza(InferenzaVarianti(), config)
    assert uguale["sorgente"] == leggi_registro()["S4"]["sorgente"]

    with open(cartella / "04_dada.R", "a", encoding="utf-8") as file:
        file.write("# uno script diverso da quello del repository\n")
    diversa = provenienza(InferenzaVarianti(), config)
    assert diversa["sorgente"] != leggi_registro()["S4"]["sorgente"]
    assert [n for n in diversa["file"] if diversa["file"][n] != uguale["file"][n]] == ["R/04_dada.R"]
    assert diversa["versione"] == InferenzaVarianti.versione
    assert diversa["immagine"] == config.run.container


class _ConScript(PipelineStep):
    """Fase doppione di S1 con uno script R dichiarato, che non esegue."""

    passo: ClassVar[Passo] = Passo.S1
    parametri: ClassVar[tuple[str, ...]] = ("filter.truncLen",)
    script_r: ClassVar[str | None] = "doppione.R"

    def calcola(self, contesto: StepContext) -> Produzione:
        """Scrive un artefatto fisso."""
        return Produzione((contesto.albero.scrivi_testo(self.cartella, "profilo.txt", "x"),))


def test_una_provenienza_diversa_avvisa_e_una_versione_nuova_rifa(tmp_path, monkeypatch):
    """
    **Obiettivo**: Verificare che, conclusa una fase con uno script R nella
    cartella di ``AMPLICON16S_R_DIR``, modificare lo script lasci la fase
    conclusa con un avviso che nomina lo script, segnalato nel log alla ripresa
    senza eseguire nulla; e che incrementare la versione la renda da rifare.

    **Razionale scientifico e sistemistico**: A parita' di versione una
    provenienza diversa e' un avviso, non un ricalcolo; la versione e'
    nell'impronta e decide la validita'.
    """
    cartella = tmp_path / "R"
    cartella.mkdir()
    (cartella / "doppione.R").write_text("# versione 1\n", encoding="utf-8")
    monkeypatch.setenv(VARIABILE_CARTELLA_R, str(cartella))
    scenario = crea_scenario(tmp_path / "s", _campioni(), con_letture=True)
    run = ProjectRun(scenario.config, passi={Passo.S0: ValidazioneIngressi(), Passo.S1: _ConScript()})
    assert Esecutore(run, fino_a=Passo.S1).esegui().conclusione is Conclusione.COMPLETATA
    manifesto = run.albero.manifesto_passo(Passo.S1, Fase.QC_PROFILES)
    assert manifesto.calcolata_su["versione"] == 1
    # Il sorgente della fase doppione: il suo script, il suo modulo (questo) e i
    # moduli del calcolo, dagli import e, poiche' e' S1 e dipende da S0, dal
    # lettore dell'inventario.
    from amplicon16s.runner.provenienza import moduli_della_fase

    assert set(manifesto.provenienza["file"]) == {"R/doppione.R", f"python/{__name__}"} | {
        f"python/{m}" for m in moduli_della_fase(_ConScript())}

    (cartella / "doppione.R").write_text("# versione 1, modificata\n", encoding="utf-8")
    situazione = run.valuta().situazioni[Passo.S1]
    assert situazione.stato is StatoPasso.COMPLETATA
    assert "R/doppione.R" in situazione.avviso
    percorso = configura(scenario.config.io.out_root)
    assert Esecutore(run, fino_a=Passo.S1).esegui().eseguite == ()
    chiudi()
    eventi = [json.loads(r) for r in percorso.read_text(encoding="utf-8").splitlines() if r]
    assert any(e.get("avviso") and "R/doppione.R" in e["avviso"] for e in eventi)

    class Seconda(_ConScript):
        versione: ClassVar[int] = 2

    nuova = ProjectRun(scenario.config, passi={Passo.S0: ValidazioneIngressi(), Passo.S1: Seconda()})
    situazione = nuova.valuta().situazioni[Passo.S1]
    assert situazione.stato is StatoPasso.DA_ESEGUIRE
    assert "versione del calcolo cambiata: 1 -> 2" in situazione.motivo


def test_i_manifesti_della_catena_registrano_il_sorgente_del_registro(dada2, finale_calcolata):
    """
    **Obiettivo**: Verificare che nei manifesti di ogni fase eseguita sulla
    versione ridotta, S0-S8 e S10-S14 (S9 e' disattivata per difetto, e la
    verifica ``tests/test_w23_filogenesi.py``), la provenienza riporti la
    versione e l'impronta del sorgente del registro.

    **Razionale scientifico e sistemistico**: Nel job del container, con
    l'immagine costruita dal repository, gli script eseguiti sono quelli del
    repository; con un'immagine che esegue script diversi questo test fallisce,
    ed e' il modo in cui la differenza si vede.
    """
    run, esito = finale_calcolata
    assert esito.conclusione is Conclusione.COMPLETATA
    registro = leggi_registro()
    for passo, fase in passi_realizzati().items():
        if passo is Passo.S9:
            continue
        manifesto = run.albero.manifesto_passo(passo, fase.cartella)
        assert manifesto.provenienza["sorgente"] == registro[str(passo)]["sorgente"], passo
        assert manifesto.provenienza["versione"] == registro[str(passo)]["versione"], passo


# --------------------------------------------------------------------------- #
# 2. G12 precondizione, il riferimento fuori da S0                              #
# --------------------------------------------------------------------------- #


def test_s0_non_dipende_dal_riferimento_tassonomico(tmp_path):
    """
    **Obiettivo**: Verificare che S0 non dichiari parametri del gruppo tax e
    che l'impronta dei suoi dati esterni non cambi cambiando il riferimento.

    **Razionale scientifico e sistemistico**: Cambiare il database deve
    invalidare S8 e le successive, non S0 e con lei l'intera catena.
    """
    assert not any(p.startswith("tax") for p in ValidazioneIngressi.parametri)
    config = config_ridotta(tmp_path)
    altro = tmp_path / "altro.fa.gz"
    altro.write_bytes(b">Bacteria;\nACGT\n")
    dati = config.model_dump(mode="python")
    dati["tax"].update(ref_fasta=str(altro), ref_md5=hashlib.md5(altro.read_bytes()).hexdigest())
    assert impronta_dati_grezzi(valida(dati)) == impronta_dati_grezzi(config)


def test_cambiare_il_riferimento_lascia_concluse_s0(tmp_path):
    """
    **Obiettivo**: Verificare che, conclusa S0, una configurazione con un altro
    riferimento valido lasci S0 conclusa.

    **Razionale scientifico e sistemistico**: E' il caso di G15 con il gruppo
    qc: G12 e' una precondizione, non un risultato di S0.
    """
    scenario = crea_scenario(tmp_path, _campioni(), con_letture=True)
    run = ProjectRun(scenario.config, passi={Passo.S0: ValidazioneIngressi()})
    assert Esecutore(run, fino_a=Passo.S0).esegui().conclusione is Conclusione.COMPLETATA
    altro = tmp_path / "altro.fa.gz"
    altro.write_bytes(b">Bacteria;\nTTTT\n")
    dati = scenario.config.model_dump(mode="python")
    dati["tax"].update(ref_fasta=str(altro), ref_md5=hashlib.md5(altro.read_bytes()).hexdigest())
    cambiata = ProjectRun(valida(dati), passi={Passo.S0: ValidazioneIngressi()})
    assert cambiata.valuta().situazioni[Passo.S0].stato is StatoPasso.COMPLETATA


def test_g12_respinge_all_avvio_un_riferimento_alterato(tmp_path):
    """
    **Obiettivo**: Verificare che, conclusa S0, alterare il file del
    riferimento fermi l'esecuzione successiva ai controlli di avvio con
    E-S0-12, a revisione umana, senza eseguire alcuna fase.

    **Razionale scientifico e sistemistico**: Tolto il riferimento
    dall'impronta di S0, G12 deve verificarlo a ogni avvio: un file alterato
    dopo S0 non deve arrivare a S8.
    """
    scenario = crea_scenario(tmp_path, _campioni(), con_letture=True)
    run = ProjectRun(scenario.config, passi={Passo.S0: ValidazioneIngressi()})
    assert Esecutore(run, fino_a=Passo.S0).esegui().conclusione is Conclusione.COMPLETATA
    Path(scenario.config.tax.ref_fasta).write_bytes(b">alterato\nTTTT\n")
    esito = Esecutore(run).esegui()
    assert esito.conclusione is Conclusione.ARRESTATA
    assert esito.eseguite == ()
    assert esito.punto.codice == "E-S0-12"
    assert esito.punto.categoria == Categoria.REVISIONE_UMANA.value
    assert esito.punto.origine == "controlli di avvio"


def test_cambiare_il_riferimento_rifa_solo_s8(dada2, tassonomia_calcolata, tmp_path):
    """
    **Obiettivo**: Verificare che, con S0-S8 concluse sulla versione ridotta,
    un altro riferimento valido lasci concluse S0-S7 e che la ripresa rifaccia
    soltanto S8.

    **Razionale scientifico e sistemistico**: Rifare S0-S7 costerebbe tre
    quarti d'ora sul dataset completo senza cambiare nulla.
    """
    altro = tmp_path / "altro.fa.gz"
    contenuto = subprocess.run(["gzip", "-dc", str(RIFERIMENTO)], capture_output=True, check=True).stdout
    altro.write_bytes(contenuto.replace(b"Genere_00", b"Genere_zero"))
    run = copia_esecuzione(
        tassonomia_calcolata, tmp_path,
        tax={"ref_fasta": str(altro), "ref_md5": hashlib.md5(altro.read_bytes()).hexdigest()},
    )
    situazioni = run.valuta().situazioni
    for passo in (Passo.S0, Passo.S1, Passo.S2, Passo.S3, Passo.S4, Passo.S5, Passo.S6, Passo.S7):
        assert situazioni[passo].stato is StatoPasso.COMPLETATA, passo
    assert situazioni[Passo.S8].motivo == "configurazione cambiata"
    esito = Esecutore(run, fino_a=Passo.S8).esegui()
    assert [r.passo for r in esito.eseguite] == [Passo.S8]


# --------------------------------------------------------------------------- #
# dada2 con i pareggi riproducibili                                            #
# --------------------------------------------------------------------------- #


def test_renv_lock_dichiara_la_correzione_di_dada2():
    """
    **Obiettivo**: Verificare che ``renv.lock`` registri dada2 come versione
    1.36.0.1, con un blocco ``Patch`` che nomina la correzione del repository,
    la versione di partenza 1.36.0 e lo SHA-256 del file della correzione.

    **Razionale scientifico e sistemistico**: Una correzione cambiata senza
    rigenerare il lock farebbe fallire la costruzione dell'immagine; qui lo si
    vede gia' nel job leggero, senza costruirla.
    """
    lock = json.loads((RADICE / "renv.lock").read_text(encoding="utf-8"))
    dada2 = lock["Packages"]["dada2"]
    assert dada2["Version"] == "1.36.0.1"
    correzione = dada2["Patch"]
    assert correzione["Base"] == "1.36.0"
    file = RADICE / "container" / "dada2" / correzione["File"]
    assert correzione["SHA256"] == hashlib.sha256(file.read_bytes()).hexdigest()


def test_il_dada2_caricato_e_quello_corretto(dada2):
    """
    **Obiettivo**: Verificare che l'R della pipeline carichi dada2 1.36.0.1,
    con la correzione dei pareggi dichiarata nel DESCRIPTION.

    **Razionale scientifico e sistemistico**: Con la versione ufficiale la
    tassonomia di S8 non sarebbe riproducibile: i pareggi fra generi si
    risolverebbero con un seme preso dal sistema.
    """
    uscita = subprocess.run(
        [str(trova_rscript()), "--vanilla", "-e",
         "d <- packageDescription('dada2'); cat(d$Version, d$Amplicon16sPatch)"],
        capture_output=True, text=True, check=True,
    ).stdout.split()
    assert uscita == ["1.36.0.1", "dada2-1.36.0-pareggi.patch"]


# --------------------------------------------------------------------------- #
# 3. La fase S8                                                                #
# --------------------------------------------------------------------------- #


def test_s8_produce_la_tabella_tassonomica(dada2, tassonomia_calcolata):
    """
    **Obiettivo**: Verificare che S8 sulla versione ridotta si concluda con
    una riga per variante della tabella di S7 in ``tassonomia.tsv`` e in
    ``bootstrap.tsv``, sei ranghi dal regno al genere, bootstrap fra 0 e 100,
    regno assegnato a ogni variante e la copertura per ciascuna classe.

    **Razionale scientifico e sistemistico**: E' la tabella che le fasi
    successive leggono; il bootstrap per rango servira' ai filtri e
    all'interpretazione.
    """
    run, esito = tassonomia_calcolata
    assert esito.conclusione is Conclusione.COMPLETATA
    cartella = run.albero.cartella(Fase.TAXONOMY)
    tassonomia = _tsv(cartella / "tassonomia.tsv")
    bootstrap = _tsv(cartella / "bootstrap.tsv")
    ranghi = ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus"]
    assert list(tassonomia[0]) == ["sequenza", *ranghi, "letture"]
    assert list(bootstrap[0]) == ["sequenza", *ranghi]
    s7 = run.albero.manifesto_passo(Passo.S7, Fase.CHIMERA).metriche
    assert len(tassonomia) == len(bootstrap) == s7["varianti_ammesse"]
    assert [r["sequenza"] for r in tassonomia] == [r["sequenza"] for r in bootstrap]
    assert all(r["Kingdom"] == "Bacteria" for r in tassonomia)
    assert all(0 <= int(r[g]) <= 100 for r in bootstrap for g in ranghi)
    riepilogo = json.loads((cartella / "riepilogo.json").read_text(encoding="utf-8"))
    assert set(riepilogo["frazione_con_phylum"]) == {c.value for c in ClasseCampione}
    assert riepilogo["varianti"] == len(tassonomia)


def test_le_varianti_del_riferimento_hanno_la_loro_linea(dada2, tassonomia_calcolata):
    """
    **Obiettivo**: Verificare che le varianti che compongono il riferimento
    sintetico ricevano ciascuna la propria linea tassonomica, fino al rango che
    il bootstrap consente, senza contraddirla in alcun rango.

    **Razionale scientifico e sistemistico**: E' il riscontro con una verita'
    nota che il riferimento sintetico permette in CI: la classificazione non
    assegna a una sequenza un taxon diverso dal proprio.
    """
    run, _ = tassonomia_calcolata
    contenuto = subprocess.run(["gzip", "-dc", str(RIFERIMENTO)], capture_output=True, check=True).stdout
    righe = contenuto.decode("ascii").split()
    attese = {righe[i + 1]: righe[i][1:].rstrip(";").split(";") for i in range(0, len(righe), 2)}
    tassonomia = {r["sequenza"]: r for r in _tsv(run.albero.cartella(Fase.TAXONOMY) / "tassonomia.tsv")}
    ranghi = ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus"]
    assegnate = 0
    for sequenza, linea in attese.items():
        assert sequenza in tassonomia
        valori = [tassonomia[sequenza][g] for g in ranghi[: len(linea)]]
        assert all(v in ("", atteso) for v, atteso in zip(valori, linea)), sequenza
        assegnate += valori[1] != ""
    assert assegnate == len(attese)


def test_le_assegnazioni_sui_taxa_difettosi_sono_marcate(dada2, tassonomia_calcolata):
    """
    **Obiettivo**: Verificare che ``difetto_riferimento.tsv`` elenchi le
    assegnazioni alla famiglia difettosa del riferimento sintetico nella
    colonna dell'ordine, dove il difetto la fa comparire, e che il riepilogo
    riporti quante varianti e quante letture ne sono interessate.

    **Razionale scientifico e sistemistico**: Il riferimento non si modifica:
    le assegnazioni che ricadono sui taxa col difetto noto si marcano, perche'
    un rango mancante li' e' un artefatto del database e non un dato.
    """
    run, _ = tassonomia_calcolata
    cartella = run.albero.cartella(Fase.TAXONOMY)
    marcate = _tsv(cartella / "difetto_riferimento.tsv")
    famiglia = [m for m in marcate if m["taxon"] == "Famiglia_10"]
    assert famiglia and all(m["colonna"] == "Order" and m["rango"] == "family" for m in famiglia)
    riepilogo = json.loads((cartella / "riepilogo.json").read_text(encoding="utf-8"))
    difetto = riepilogo["difetto_riferimento"]
    assert difetto["varianti"] == len({m["sequenza"] for m in marcate})
    assert difetto["letture"] == sum({m["sequenza"]: int(m["letture"]) for m in marcate}.values())
    assert difetto["letture"] > 0


def test_la_marcatura_cerca_il_nome_in_ogni_colonna():
    """
    **Obiettivo**: Verificare che ``marca_difetti`` trovi un taxon difettoso
    in qualunque colonna e ignori le assegnazioni vuote e i taxa non elencati.

    **Razionale scientifico e sistemistico**: Il difetto e' un rango mancante:
    il nome compare nella colonna del rango superiore, non nella propria.
    """
    righe = [
        {"sequenza": "A", "Order": "Famiglia_X", "Family": "Genere_Y", "letture": "10"},
        {"sequenza": "B", "Order": "Ordine_1", "Family": "", "letture": "5"},
    ]
    marcate = marca_difetti(righe, ["Order", "Family"], {"Famiglia_X": "family", "Genere_Y": "genus"})
    assert [(m["sequenza"], m["colonna"], m["taxon"]) for m in marcate] == [
        ("A", "Order", "Famiglia_X"), ("A", "Family", "Genere_Y"),
    ]


def test_s8_da_gli_stessi_byte_in_due_esecuzioni_e_con_thread_diversi(
    dada2, tassonomia_calcolata, tmp_path
):
    """
    **Obiettivo**: Verificare che S8, rieseguita sulla stessa tabella con lo
    stesso numero di thread e con uno diverso, produca in ``08_taxonomy`` gli
    stessi byte, e che cambiare ``run.threads`` non la renda da rifare.

    **Razionale scientifico e sistemistico**: Il bootstrap e' casuale e il seme
    viene da ``run.seed``; ``run.threads`` e' escluso dall'impronta perche'
    ritenuto senza effetto: se l'esito dipendesse dai thread l'esclusione
    sarebbe sbagliata.
    """
    run, _ = tassonomia_calcolata
    attese = _impronte(run.albero.cartella(Fase.TAXONOMY))
    # Su una macchina con un solo processore i due valori coincidono: le due
    # copie stanno comunque in cartelle distinte.
    for numero, thread in enumerate((run.config.run.threads, 1)):
        copia = copia_esecuzione(
            tassonomia_calcolata, tmp_path / f"copia{numero}", run={"threads": thread}
        )
        assert copia.valuta().situazioni[Passo.S8].stato is StatoPasso.COMPLETATA
        copia.albero.rimuovi_manifesto_passo(Passo.S8, Fase.TAXONOMY)
        esito = Esecutore(copia, fino_a=Passo.S8).esegui()
        assert [r.passo for r in esito.eseguite] == [Passo.S8]
        assert _impronte(copia.albero.cartella(Fase.TAXONOMY)) == attese, thread


def test_e_s8_02_ferma_con_una_copertura_insufficiente(dada2, tassonomia_calcolata, tmp_path):
    """
    **Obiettivo**: Verificare che S8, con ``qc.min_frac_phylum`` sopra la
    copertura della versione ridotta, si fermi con E-S8-02 di revisione
    umana, senza tentativi e senza manifesto, lasciando il riepilogo con le
    misure.

    **Razionale scientifico e sistemistico**: Una copertura insufficiente
    indica un riferimento che non copre la regione amplificata o un bootstrap
    troppo severo: si capisce, non si ritenta.
    """
    run = copia_esecuzione(tassonomia_calcolata, tmp_path, qc={"min_frac_phylum": 0.9999})
    esito = Esecutore(run, fino_a=Passo.S8).esegui()
    assert esito.conclusione is Conclusione.ARRESTATA
    assert esito.punto.passo is Passo.S8
    assert esito.punto.codice == "E-S8-02"
    assert esito.punto.categoria == Categoria.REVISIONE_UMANA.value
    assert esito.punto.tentativi == 1
    assert run.albero.manifesto_passo(Passo.S8, Fase.TAXONOMY) is None
    assert (run.albero.cartella(Fase.TAXONOMY) / "riepilogo.json").exists()


def test_la_copertura_si_controlla_sulle_classi_scelte():
    """
    **Obiettivo**: Verificare che ``controlla_copertura`` sollevi E-S8-02
    per una classe controllata sotto soglia e ignori i controlli negativi,
    qualunque sia la loro copertura.

    **Razionale scientifico e sistemistico**: I bianchi a bassa biomassa
    raccolgono contaminanti e rumore: la loro copertura non misura la
    qualita' della classificazione dei campioni.
    """
    config = config_ridotta(Path("/nonusata"))
    buone = {c.value: {"varianti": 0.95, "letture": 0.99} for c in ClasseCampione}
    controlla_copertura({**buone, "controllo_negativo": {"varianti": 0.1, "letture": 0.1}}, config.qc)
    with pytest.raises(ErrorePipeline) as info:
        controlla_copertura({**buone, "biologico": {"varianti": 0.5, "letture": 0.9}}, config.qc)
    assert info.value.codice == "E-S8-02"
    assert "biologico" in info.value.dettaglio


# --------------------------------------------------------------------------- #
# Dataset completo, nel container                                              #
# --------------------------------------------------------------------------- #


@pytest.mark.dati_reali
def test_s8_sul_dataset_completo_e_variovorax_nei_controlli_positivi(dada2, catena_reale):
    """
    **Obiettivo**: Verificare che sul dataset completo S8 si concluda dopo
    S0-S7 con SILVA 138, e che nei controlli positivi ai due livelli di
    diluizione piu' alti (2.000.000 e 400.000 cellule) la variante dominante
    sia assegnata al genere Variovorax in ogni piastra; riportare durata,
    copertura per classe, taxa difettosi e la dominante di ogni controllo.

    **Razionale scientifico e sistemistico**: I controlli positivi contengono
    un solo ceppo, Variovorax sp. OAS795: e' un riscontro con una verita' nota
    sull'intera catena, dal FASTQ alla tassonomia. Vale dove il ceppo domina
    la biomassa: nella serie KatharoSeq, scendendo a poche centinaia di
    cellule, prevalgono i contaminanti di reagente (misurato: livelli 6-8), e
    nelle piastre 1 e 2 un contaminante cloroplastico domina gia' dai livelli
    3-5. Sono misure, non un difetto della classificazione.
    """
    run, esito = catena_reale
    assert esito.conclusione is Conclusione.COMPLETATA
    s8 = next(r for r in esito.eseguite if r.passo is Passo.S8)
    print(f"\nS8: {s8.secondi} s, {dict(s8.metriche)}")
    cartella = run.albero.cartella(Fase.TAXONOMY)
    genere = {r["sequenza"]: r["Genus"] for r in _tsv(cartella / "tassonomia.tsv")}
    positivi = {c.accession: c.nome for c in run.valuta().inventario.controlli_positivi}
    codice = (
        f"t <- readRDS('{run.albero.cartella(Fase.CHIMERA) / 'tabella_asv.rds'}'); "
        f"p <- intersect(c({', '.join(repr(a) for a in positivi)}), rownames(t)); "
        "cat(paste(p, colnames(t)[apply(t[p, , drop = FALSE], 1, which.max)], sep = '\\t'), sep = '\\n')"
    )
    uscita = subprocess.run(
        [str(trova_rscript()), "--vanilla", "-e", codice], capture_output=True, text=True, check=True,
    ).stdout
    dominanti = {positivi[a]: genere[s] or "-" for a, s in (r.split("\t") for r in uscita.splitlines() if r)}
    assert len(dominanti) == len(positivi) == 80
    print(" ".join(f"{n}:{g}" for n, g in sorted(dominanti.items())))
    alti = {n: g for n, g in dominanti.items() if n.rsplit(".", 1)[1] in ("1", "2")}
    assert len(alti) == 20
    assert {n: g for n, g in alti.items() if g != "Variovorax"} == {}
