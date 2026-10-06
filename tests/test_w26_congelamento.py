r"""Suite di test della settimana 26: analisi di sensibilita' e congelamento.

1. Inquadramento nel Piano Operativo
------------------------------------
Settimana 26 (W26), Fase F7 (regola di decisione dell'analisi di sensibilita',
congelamento della configurazione di OSD-734, regola rigorosa sulla
provenienza).

2. Moduli sorgente coperti
--------------------------
* ``scripts/sensitivity.py`` (la regola: ``jaccard``, ``cambiamento``,
  ``applica_regola``, ``scegli_modalita``; le griglie e le soglie)
* ``docs/decision_log.md``, ``docs/sensibilita/`` (misure, passi, decisioni)
* ``src/amplicon16s/runner/provenienza.py`` (``verifica_git``, ``blocco_r``,
  ``discordanze_ambiente``, ``impronta_rigorosa``)
* ``src/amplicon16s/runner/executor.py`` (controlli di avvio con la regola
  rigorosa), ``src/amplicon16s/runner/project.py``, ``src/amplicon16s/steps/base.py``
* ``R/00_ambiente.R``
* ``src/amplicon16s/config/schema.py``, ``src/amplicon16s/config/defaults.py``
  (``run.strict_provenance``), ``config/config.example.yaml``,
  ``dati/osd734/config_osd734.yaml``
* ``src/amplicon16s/report/builder.py`` (colonna della regola rigorosa)

3. Cosa valuta questo file
--------------------------
- la regola di decisione, come funzione pura: valore mantenuto se non
  instabile; instabile solo oltre la soglia assoluta e quella di sproporzione
  insieme; un vicino non ammissibile rende instabile; sostituzione con
  l'alternativa in zona stabile piu' vicina; regola delle modalita';
- griglie e soglie dello strumento uguali a quelle del registro delle
  decisioni; decisioni pubblicate riottenibili dai passi pubblicati; valori
  della configurazione congelata uguali a quelli scelti;
- ``run.strict_provenance`` falso per difetto e vero nella configurazione
  congelata; senza la regola l'impronta di fase non cambia, con la regola
  comprende sorgente, file di blocco e immagine dichiarata;
- rifiuto dell'esecuzione, prima di qualunque fase, con git non leggibile
  (``E-PROV-01``), con script R fuori dal repository (``E-PROV-01``), con
  modifiche non committate (``E-PROV-02``), con ambiente R diverso dal file di
  blocco (``E-PROV-03``); esito positivo nell'ambiente dell'immagine;
- con la regola attiva un'immagine dichiarata diversa rende la fase da rifare,
  senza la regola resta conclusa; il report riporta l'esito della regola.

4. Comandi Bash e scenari di esecuzione
---------------------------------------
    ``<immagine>`` e' l'immagine del container della pipeline; quella corrente
    e' indicata in ``test.txt``, sezione 1.3.

    1. Modalità locale standard (R di base con jsonlite, senza Bioconductor né
       dati reali):
       pytest tests/test_w26_congelamento.py -v

    2. Modalità container Docker standard (sottoinsieme ridotto con
       R/Bioconductor):
       docker run --rm \
         -e PYTHONPATH=/app/src \
         -e AMPLICON16S_R_DIR=/app/R \
         -v "$(pwd)":/app \
         -w /app \
         <immagine> \
         pytest -o cache_dir=/tmp/.pytest_cache tests/test_w26_congelamento.py -v

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
         pytest -o cache_dir=/tmp/.pytest_cache tests/test_w26_congelamento.py -v

5. Risultato atteso
-------------------
Vedi ``test.txt``, scheda W26.

6. Razionale scientifico e sistemistico
---------------------------------------
- Una regola di decisione fissata prima dei risultati vale solo se lo strumento
  che la applica e' quello descritto: griglie, soglie e decisioni pubblicate
  devono derivare l'una dall'altra.
- Una configurazione congelata certifica i risultati solo se il codice che li
  calcola e' identificato da un commit e l'ambiente R e' quello del file di
  blocco: cio' che non si puo' verificare deve fermare l'esecuzione, non
  passare in silenzio.
"""

from __future__ import annotations

import csv
import importlib.util
import json
import re
import subprocess
from pathlib import Path
from typing import Any, ClassVar

import pytest
import yaml
from conftest import NEGATIVO, POSITIVO, Campione, crea_scenario
from sottoinsieme import motivo_pacchetti_r_assenti

import amplicon16s.cli as cli
import amplicon16s.runner.executor as executor
import amplicon16s.runner.provenienza as provenienza
from amplicon16s.config import defaults
from amplicon16s.config.resolve import risolvi
from amplicon16s.config.schema import Config, carica
from amplicon16s.errors.catalog import CATALOGO, Categoria
from amplicon16s.io_layer.artifacts import Fase
from amplicon16s.logging.logger import NOME_FILE_LOG, chiudi
from amplicon16s.rbridge.runner import cartella_r
from amplicon16s.report.builder import CARTELLA_REPORT, CARTELLA_TABELLE
from amplicon16s.runner.executor import EVENTO_CONTROLLI, NOME_PUNTO_DI_RIPRESA, Esecutore
from amplicon16s.runner.graph import Passo
from amplicon16s.runner.project import ProjectRun, StatoPasso
from amplicon16s.steps.s00_validate import ValidazioneIngressi

RADICE = Path(__file__).resolve().parents[1]
SENSIBILITA = RADICE / "docs" / "sensibilita"
CONGELATA = RADICE / "dati" / "osd734" / "config_osd734.yaml"
senza_dati = pytest.mark.skipif(
    not CONGELATA.is_file(),
    reason="la cartella dati/ non e' presente (nell'immagine non viene copiata)",
)
senza_docs = pytest.mark.skipif(
    not SENSIBILITA.is_dir(),
    reason="la cartella docs/ non e' presente (nell'immagine non viene copiata)",
)
IMMAGINE_DIVERSA = "amplicon16s@sha256:" + "a" * 64


def _strumento():
    """Il modulo ``scripts/sensitivity.py``, caricato dal file: la cartella degli
    script non e' un pacchetto.
    """
    specifica = importlib.util.spec_from_file_location("sensitivity", RADICE / "scripts" / "sensitivity.py")
    modulo = importlib.util.module_from_spec(specifica)
    specifica.loader.exec_module(modulo)
    return modulo


S = _strumento()


def _passo(varianti: float = 0.0, campioni: float = 0.0, letture: float = 0.0) -> dict[str, float]:
    """Un cambiamento di passo con le tre misure della regola."""
    return {"varianti": varianti, "campioni": campioni, "letture": letture}


def _tsv(percorso: Path) -> list[dict[str, str]]:
    """Le righe di una tabella TSV con intestazione."""
    with percorso.open(encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


@pytest.fixture(autouse=True)
def uscite_pulite():
    """Chiude le uscite del log prima e dopo ogni test, perche' nessun handler
    resti aperto sulla cartella temporanea.
    """
    chiudi()
    yield
    chiudi()


@pytest.fixture
def scenario(tmp_path):
    """Lo scenario su disco con tre campioni e le loro letture."""
    campioni = [
        Campione("ERX3000001", "NOD1D4.L1"),
        Campione("ERX3000003", "POS.P1.1", materiale=POSITIVO, posizione="Not Applicable"),
        Campione("ERX3000004", "BLANK.P1.1", materiale=NEGATIVO, posizione="Not Applicable"),
    ]
    return crea_scenario(tmp_path, campioni, con_letture=True)


def _scrivi_config(scenario, percorso: Path, **run: Any) -> Path:
    """Scrive la configurazione minima dello scenario, con i parametri del gruppo
    ``run`` indicati.
    """
    config = scenario.config
    dati = {
        "io": {n: str(getattr(config.io, n)) for n in ("fastq_dir", "assay_table", "study_table", "out_root")},
        "tax": {n: str(getattr(config.tax, n)) for n in ("ref_fasta", "ref_md5", "ref_name", "ref_version")},
        "run": {"container": config.run.container, "threads": 1, **run},
    }
    percorso.write_text(yaml.safe_dump(dati, sort_keys=False), encoding="utf-8")
    return percorso


@pytest.fixture
def rigorosa(scenario, tmp_path) -> Path:
    """La configurazione dello scenario con la regola rigorosa attiva."""
    return _scrivi_config(scenario, tmp_path / "rigorosa.yaml", strict_provenance=True)


def _git(cartella: Path, *argomenti: str) -> None:
    """Esegue git nella cartella, con un'identita' di prova."""
    subprocess.run(
        ["git", "-C", str(cartella), "-c", "user.name=prova", "-c", "user.email=prova@example.org",
         *argomenti],
        check=True, capture_output=True,
    )


@pytest.fixture
def clone(tmp_path, monkeypatch) -> Path:
    """Un repository git di prova con ``src/`` e ``R/`` committati, che la
    provenienza considera la radice del repository della pipeline.
    """
    radice = tmp_path / "clone"
    (radice / "src").mkdir(parents=True)
    (radice / "R").mkdir()
    (radice / "src" / "modulo.py").write_text("x = 1\n", encoding="utf-8")
    (radice / "R" / "fase.R").write_text("x <- 1\n", encoding="utf-8")
    _git(radice, "init", "-q")
    _git(radice, "add", ".")
    _git(radice, "commit", "-q", "-m", "stato iniziale")
    monkeypatch.setattr(provenienza, "RADICE_REPOSITORY", radice)
    return radice


def _cli(*argomenti: str, capsys) -> int:
    """Esegue la riga di comando e ne restituisce il codice di uscita."""
    codice = cli.main(list(argomenti))
    capsys.readouterr()
    return codice


def _controlli(radice: Path) -> dict[str, Any]:
    """L'ultimo evento dei controlli di avvio nel log strutturato."""
    righe = (Path(radice) / Fase.LOGS.value / NOME_FILE_LOG).read_text(encoding="utf-8").splitlines()
    return [e for e in map(json.loads, righe) if e.get("evento") == EVENTO_CONTROLLI][-1]


def _punto(radice: Path) -> dict[str, Any]:
    """Il punto di ripresa dichiarato dall'esecuzione."""
    percorso = Path(radice) / Fase.LOGS.value / f"{NOME_PUNTO_DI_RIPRESA}.json"
    return json.loads(percorso.read_text(encoding="utf-8"))


def _rifiutata(scenario, config: Path, codice: str, capsys) -> dict[str, Any]:
    """Avvia la validazione, verifica che sia rifiutata con il codice indicato
    prima di qualunque fase, e restituisce l'esito della regola nel log.
    """
    assert _cli("validate", "--config", str(config), capsys=capsys) == cli.USCITA_ARRESTO
    radice = scenario.config.io.out_root
    punto = _punto(radice)
    assert (punto["codice"], punto["origine"]) == (codice, "controlli di avvio")
    assert not (Path(radice) / Fase.CONFIG.value / "manifest_S0.json").exists()
    assert not list(Path(radice).rglob("manifest_S*.json"))
    return _controlli(radice)["provenienza_rigorosa"]


# --------------------------------------------------------------------------- #
# 1. La regola di decisione                                                    #
# --------------------------------------------------------------------------- #


def test_i_cambiamenti_di_passo_sono_quelli_della_regola():
    """
    **Obiettivo**: Verificare che ``jaccard`` e ``cambiamento`` calcolino le tre
    misure come definite dal registro: 1 meno Jaccard sulle varianti, campioni
    che cambiano stato sui biologici del dataset, differenza assoluta di
    letture sulle letture correnti; nessun cambiamento se una delle due
    configurazioni non e' ammissibile.

    **Razionale scientifico e sistemistico**: La decisione dipende solo da
    queste tre grandezze: un denominatore diverso da quello dichiarato
    sposterebbe i valori rispetto alle soglie senza che nulla lo segnali.
    """
    assert S.jaccard({"a", "b", "c"}, {"b", "c", "d"}) == 0.5
    assert S.jaccard(set(), set()) == 1.0
    a = {"ammissibile": True, "varianti": ["v1", "v2", "v3"], "campioni": ["c1", "c2"], "letture": 1000}
    b = {"ammissibile": True, "varianti": ["v2", "v3", "v4"], "campioni": ["c2", "c3"], "letture": 900}
    assert S.cambiamento(a, b, biologici=10, letture_correnti=2000) == {
        "varianti": 0.5, "campioni": 0.2, "letture": 0.05,
    }
    assert S.cambiamento(a, b, 10, 2000) == S.cambiamento(b, a, 10, 2000)
    assert S.cambiamento(a, {"ammissibile": False}, 10, 2000) is None


def test_un_valore_non_instabile_si_mantiene():
    """
    **Obiettivo**: Verificare che con cambiamenti di passo sotto le soglie
    assolute il valore corrente sia mantenuto, anche quando esistono
    alternative in zona stabile.

    **Razionale scientifico e sistemistico**: Il principio della regola e' la
    conservazione: un'alternativa altrettanto stabile non e' una ragione per
    cambiare un valore, altrimenti la scelta tornerebbe arbitraria.
    """
    passi = [_passo(0.02), _passo(0.03), _passo(0.04), _passo(0.03)]
    d = S.applica_regola([1.0, 2.0, 3.0, 4.0, 5.0], 3.0, passi)
    assert (d["instabile"], d["scelto"]) == (False, 3.0)
    assert d["alternative_in_zona_stabile"] == [2.0, 4.0]
    assert d["esito"].startswith("mantenuto")


def test_l_instabilita_richiede_entrambe_le_soglie():
    """
    **Obiettivo**: Verificare che un cambiamento oltre la soglia assoluta ma
    non oltre tre volte la mediana dei passi non renda instabile il valore
    corrente, e che lo renda instabile quando supera entrambe.

    **Razionale scientifico e sistemistico**: Un parametro che cambia il
    risultato a ogni passo in misura simile non ha una zona piu' stabile di
    un'altra: spostare il valore non lo renderebbe piu' affidabile. E' il caso
    misurato per la soglia di prevalenza.
    """
    uniformi = [_passo(0.20), _passo(0.18), _passo(0.11), _passo(0.10)]
    d = S.applica_regola([1.0, 2.0, 3.0, 4.0, 5.0], 2.0, uniformi)
    assert d["instabile"] is False and d["scelto"] == 2.0
    assert d["soglie_di_sproporzione"]["varianti"] == pytest.approx(3 * 0.145)

    scarto = [_passo(0.01), _passo(0.02), _passo(campioni=0.30), _passo(0.01), _passo(0.02)]
    d = S.applica_regola([1.0, 2.0, 3.0, 4.0, 5.0, 6.0], 3.0, scarto)
    assert d["instabile"] is True
    assert len(d["motivi_di_instabilita"]) == 1 and "campioni" in d["motivi_di_instabilita"][0]


def test_un_valore_instabile_si_sostituisce_con_l_alternativa_stabile_piu_vicina():
    """
    **Obiettivo**: Verificare che un valore instabile sia sostituito
    dall'alternativa in zona stabile piu' vicina, che a parita' di distanza i
    candidati siano entrambi riportati, e che senza alternative stabili il
    valore resti con l'instabilita' dichiarata come limite noto.

    **Razionale scientifico e sistemistico**: La sostituzione e' ammessa solo
    verso una zona in cui il risultato non dipende dal valore esatto; in
    mancanza, cambiare sposterebbe l'instabilita' senza toglierla.
    """
    valori = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0]
    passi = [_passo(0.01), _passo(0.01), _passo(0.40), _passo(0.01), _passo(0.01), _passo(0.01)]
    d = S.applica_regola(valori, 3.0, passi)
    assert d["instabile"] and d["alternative_in_zona_stabile"] == [2.0, 5.0, 6.0]
    assert d["candidati"] == [2.0] and d["esito"].startswith("sostituito")

    d = S.applica_regola(valori, 4.0, passi)
    assert d["candidati"] == [5.0]

    simmetrici = [_passo(0.01), _passo(0.01), _passo(0.40), _passo(0.40), _passo(0.01), _passo(0.01)]
    d = S.applica_regola(valori, 4.0, simmetrici)
    assert d["candidati"] == [2.0, 6.0]

    senza = [_passo(0.20), _passo(0.90), _passo(0.20)]
    d = S.applica_regola([1.0, 2.0, 3.0, 4.0], 2.0, senza)
    assert d["instabile"] and d["scelto"] == 2.0 and "limite noto" in d["esito"]


def test_un_vicino_non_ammissibile_rende_instabile_il_valore():
    """
    **Obiettivo**: Verificare che, se un vicino del valore corrente non e'
    ammissibile (passo non definito), il valore sia instabile da quel lato, che
    il passo non definito non entri nella mediana, e che un'alternativa con un
    vicino non ammissibile non sia in zona stabile.

    **Razionale scientifico e sistemistico**: Un valore adiacente a una
    configurazione che la pipeline ferma e' al bordo della regione in cui la
    catena si conclude: e' la forma piu' netta di instabilita'.
    """
    passi = [_passo(0.01), None, _passo(0.02), _passo(0.03)]
    d = S.applica_regola([1.0, 2.0, 3.0, 4.0, 5.0], 3.0, passi)
    assert d["instabile"] and d["motivi_di_instabilita"] == ["vicino inferiore non ammissibile"]
    assert d["mediane_dei_passi"]["varianti"] == 0.02
    assert d["alternative_in_zona_stabile"] == [4.0]
    assert d["candidati"] == [4.0]


def test_la_modalita_corrente_si_mantiene_se_ammissibile():
    """
    **Obiettivo**: Verificare la regola delle modalita' di decontaminazione:
    la corrente resta se ammissibile, qualunque sia la somiglianza delle altre;
    se non lo e', si adotta l'ammissibile piu' simile; se nessuna lo e', resta.

    **Razionale scientifico e sistemistico**: Le modalita' rispondono a domande
    statistiche diverse e non hanno vicini: la scelta e' di metodo e cambia
    solo se viola un vincolo che la pipeline gia' dichiara.
    """
    somiglianza = {"a": 1.0, "b": 0.6, "c": 0.9}
    d = S.scegli_modalita("a", {"a": True, "b": False, "c": True}, somiglianza)
    assert d["scelto"] == "a" and d["esito"].startswith("mantenuta")
    d = S.scegli_modalita("b", {"a": False, "b": False, "c": True, "d": True}, {**somiglianza, "d": 0.95})
    assert d["scelto"] == "d" and d["esito"].startswith("sostituita")
    d = S.scegli_modalita("a", {"a": False, "b": False}, somiglianza)
    assert d["scelto"] == "a" and "limite noto" in d["esito"]


@senza_docs
def test_griglie_e_soglie_dello_strumento_sono_quelle_del_registro():
    """
    **Obiettivo**: Verificare che le griglie, i valori correnti, le soglie
    assolute e il fattore di sproporzione di ``scripts/sensitivity.py``
    coincidano con quelli scritti nel registro delle decisioni, e che in ogni
    griglia il valore corrente abbia due vicini a passo uniforme.

    **Razionale scientifico e sistemistico**: La regola e' stata fissata nel
    registro prima dei risultati: lo strumento che la applica non puo'
    discostarsene senza che la garanzia cada.
    """
    registro = (RADICE / "docs" / "decision_log.md").read_text(encoding="utf-8")

    def numeri(testo: str) -> tuple[float, ...]:
        return tuple(float(n.replace(",", ".")) for n in re.findall(r"\d+,\d+", testo))

    for parametro, (corrente, griglia) in S.GRIGLIE.items():
        riga = next(r for r in registro.splitlines() if r.startswith(f"| `{parametro}` |"))
        _, _, scritto, valori, _ = (c.strip() for c in riga.split("|"))
        assert numeri(scritto) == (corrente,) and numeri(valori) == griglia, parametro
        i = griglia.index(corrente)
        assert 0 < i < len(griglia) - 1
        passi = {round(b - a, 9) for a, b in zip(griglia, griglia[1:], strict=False)}
        assert len(passi) == 1, parametro
    assert S.SOGLIE_ASSOLUTE == {"varianti": 0.10, "campioni": 0.05, "letture": 0.05}
    assert "0,10 per le varianti, 0,05 per i campioni, 0,05 per le" in registro
    assert S.FATTORE_SPROPORZIONE == 3.0 and "tre volte la mediana" in registro
    assert S.MODALITA[S.MODALITA_CORRENTE] == {"mode": "aggregate", "batch_combine": "minimum"}
    assert set(S.MODALITA) == {"aggregata", "per piastra, minimum", "per piastra, fisher"}


@senza_docs
def test_le_decisioni_pubblicate_derivano_dai_passi_pubblicati():
    """
    **Obiettivo**: Verificare che, applicando la regola ai cambiamenti di passo
    pubblicati in ``docs/sensibilita/passi.tsv``, si riottengano per ogni
    parametro l'instabilita' e il valore scelto di ``decisioni.json``; che la
    tabella delle misure contenga ogni valore di ogni griglia, con 0,80 e 0,90
    per la sensibilita' della curva, e le tre modalita' di decontaminazione.

    **Razionale scientifico e sistemistico**: Il rapporto di sensibilita' e'
    credibile se chi legge puo' rifare la decisione dalle tabelle, senza
    rieseguire la pipeline.
    """
    decisioni = json.loads((SENSIBILITA / "decisioni.json").read_text(encoding="utf-8"))
    passi = _tsv(SENSIBILITA / "passi.tsv")
    misure = _tsv(SENSIBILITA / "misure.tsv")
    for parametro, (corrente, griglia) in S.GRIGLIE.items():
        del_parametro = [
            {m: float(r[m]) for m in S.MISURE} for r in passi if r["parametro"] == parametro
        ]
        assert len(del_parametro) == len(griglia) - 1
        rifatta = S.applica_regola(list(griglia), corrente, del_parametro)
        pubblicata = decisioni["decisioni"][parametro]
        assert rifatta["instabile"] == pubblicata["instabile"], parametro
        assert rifatta["scelto"] == pubblicata["scelto"], parametro
        misurati = {float(r["valore"]) for r in misure if r["parametro"] == parametro}
        assert set(griglia) | set(S.FUORI_GRIGLIA.get(parametro, ())) == misurati
    assert {0.80, 0.90} <= set(S.GRIGLIE["katharoseq.target_sensitivity"][1])
    modalita = {r["valore"] for r in misure if r["parametro"].startswith("modalita")}
    assert modalita == set(S.MODALITA)
    assert decisioni["soglie_assolute"] == S.SOGLIE_ASSOLUTE


@senza_dati
@senza_docs
def test_la_configurazione_congelata_ha_i_valori_scelti_e_la_regola_attiva():
    """
    **Obiettivo**: Verificare che ``dati/osd734/config_osd734.yaml`` dichiari
    in modo esplicito i valori scelti dall'analisi di sensibilita' per i
    quattro parametri, attivi ``run.strict_provenance`` e dichiari un'immagine
    ancorata per digest di registro; che la configurazione sia valida.

    **Razionale scientifico e sistemistico**: Congelare significa che i valori
    non dipendono piu' dai predefiniti del codice, che possono cambiare: sono
    scritti, e l'esecuzione che li usa e' vincolata a codice e ambiente
    verificati.
    """
    dati = yaml.safe_load(CONGELATA.read_text(encoding="utf-8"))
    scelte = json.loads((SENSIBILITA / "decisioni.json").read_text(encoding="utf-8"))["decisioni"]
    assert dati["katharoseq"]["target_sensitivity"] == scelte["katharoseq.target_sensitivity"]["scelto"]
    assert dati["decontam"]["threshold"] == scelte["decontam.threshold"]["scelto"]
    assert dati["prev"]["min_fraction"] == scelte["prev.min_fraction"]["scelto"]
    modalita = S.MODALITA[scelte["modalita' di decontaminazione"]["scelto"]]
    assert {k: dati["decontam"][k] for k in ("mode", "batch_combine")} == modalita
    assert dati["run"]["strict_provenance"] is True
    assert re.fullmatch(r"ghcr\.io/[a-z0-9-]+/[a-z0-9._-]+@sha256:[0-9a-f]{64}", dati["run"]["container"])
    Config.model_validate(dati)


# --------------------------------------------------------------------------- #
# 2. La regola rigorosa sulla provenienza                                      #
# --------------------------------------------------------------------------- #


def test_la_regola_rigorosa_e_disattivata_per_difetto(scenario):
    """
    **Obiettivo**: Verificare che ``run.strict_provenance`` valga falso per
    difetto e nella configurazione di esempio, che accetti solo un booleano, e
    che senza la regola l'impronta su cui una fase e' calcolata non contenga la
    chiave ``provenienza``.

    **Razionale scientifico e sistemistico**: Durante lo sviluppo il codice ha
    quasi sempre modifiche non committate: la regola deve essere una scelta
    della configurazione congelata, e la sua assenza non deve cambiare le
    impronte delle esecuzioni esistenti.
    """
    assert defaults.RUN_STRICT_PROVENANCE is False
    assert scenario.config.run.strict_provenance is False
    esempio = yaml.safe_load((RADICE / "config" / "config.example.yaml").read_text(encoding="utf-8"))
    assert esempio["run"]["strict_provenance"] is False
    dati = scenario.config.model_dump(mode="python")
    dati["run"]["strict_provenance"] = "true"
    with pytest.raises(ValueError):
        Config.model_validate(dati)
    calcolata = ValidazioneIngressi().calcolata_su(risolvi(scenario.config), {})
    assert "provenienza" not in calcolata


def test_con_la_regola_l_impronta_di_fase_comprende_sorgente_ambiente_e_immagine(scenario):
    """
    **Obiettivo**: Verificare che con la regola attiva l'impronta su cui una
    fase e' calcolata contenga l'impronta del sorgente della fase, quella del
    file di blocco dell'ambiente R e l'immagine dichiarata; che il resto
    dell'impronta sia invariato; che l'impronta cambi con l'immagine dichiarata.

    **Razionale scientifico e sistemistico**: Finche' sorgente e ambiente
    stanno solo nella provenienza, una differenza produce un avviso; nella
    configurazione congelata deve rendere la fase non piu' valida.
    """
    fase = ValidazioneIngressi()
    senza = fase.calcolata_su(risolvi(scenario.config), {})
    dati = scenario.config.model_dump(mode="python")
    dati["run"]["strict_provenance"] = True
    attiva = Config.model_validate(dati)
    con = fase.calcolata_su(risolvi(attiva), {})
    assert {k: v for k, v in con.items() if k != "provenienza"} == senza
    aggiunta = con["provenienza"]
    assert set(aggiunta) == {"sorgente", "blocco_r", "immagine"}
    sorgente, _ = provenienza.impronta_sorgente(provenienza.file_del_sorgente(fase, cartella_r()))
    assert aggiunta["sorgente"] == sorgente and re.fullmatch(r"sha256:[0-9a-f]{64}", sorgente)
    assert re.fullmatch(r"sha256:[0-9a-f]{64}", aggiunta["blocco_r"])
    assert aggiunta["immagine"] == attiva.run.container

    dati["run"]["container"] = IMMAGINE_DIVERSA
    altra = fase.calcolata_su(risolvi(Config.model_validate(dati)), {})
    assert altra["provenienza"]["immagine"] == IMMAGINE_DIVERSA
    assert altra["provenienza"]["sorgente"] == aggiunta["sorgente"]


def test_la_verifica_di_git_distingue_i_quattro_casi(clone, tmp_path, monkeypatch):
    """
    **Obiettivo**: Verificare ``verifica_git``: repository pulito, commit
    restituito; file di ``src/`` o ``R/`` modificato o non tracciato,
    ``E-PROV-02`` con il nome del file; modifica fuori da ``src/`` e ``R/``,
    nessun rifiuto; script R eseguiti da una cartella fuori dal repository,
    ``E-PROV-01``; cartella che non e' un repository, ``E-PROV-01``; git non
    eseguibile, ``E-PROV-01``.

    **Razionale scientifico e sistemistico**: Il commit identifica il codice
    solo se i file eseguiti sono quelli del commit: ogni caso in cui non e'
    dimostrabile deve dare un rifiuto, e git non leggibile non e' un via libera.
    """
    script = clone / "R"
    codice, commit = provenienza.verifica_git(script)
    assert codice is None and re.fullmatch(r"[0-9a-f]{40}", commit)

    (clone / "config.yaml").write_text("run: {}\n", encoding="utf-8")
    assert provenienza.verifica_git(script)[0] is None

    (clone / "src" / "modulo.py").write_text("x = 2\n", encoding="utf-8")
    codice, dettaglio = provenienza.verifica_git(script)
    assert codice == "E-PROV-02" and "src/modulo.py" in dettaglio
    _git(clone, "checkout", "--", "src/modulo.py")
    (clone / "R" / "nuovo.R").write_text("y <- 1\n", encoding="utf-8")
    codice, dettaglio = provenienza.verifica_git(script)
    assert codice == "E-PROV-02" and "R/nuovo.R" in dettaglio
    (clone / "R" / "nuovo.R").unlink()

    altrove = tmp_path / "script_altrove"
    altrove.mkdir()
    codice, dettaglio = provenienza.verifica_git(altrove)
    assert codice == "E-PROV-01" and "fuori dal repository" in dettaglio

    spoglia = tmp_path / "senza_git"
    (spoglia / "R").mkdir(parents=True)
    monkeypatch.setattr(provenienza, "RADICE_REPOSITORY", spoglia)
    assert provenienza.verifica_git(spoglia / "R")[0] == "E-PROV-01"

    monkeypatch.setattr(provenienza, "RADICE_REPOSITORY", clone)

    def assente(*argomenti, **opzioni):
        raise FileNotFoundError("git")

    monkeypatch.setattr(provenienza.subprocess, "run", assente)
    assert provenienza.verifica_git(script)[0] == "E-PROV-01"


def test_con_git_non_leggibile_l_esecuzione_e_rifiutata(scenario, rigorosa, tmp_path, monkeypatch, capsys):
    """
    **Obiettivo**: Verificare che, con la regola attiva e il codice fuori da un
    repository git, l'avvio si fermi con codice di uscita 4 e ``E-PROV-01``
    prima di qualunque fase, che il log riporti la regola attiva con git non
    verificato, e che la stessa configurazione senza la regola venga eseguita.

    **Razionale scientifico e sistemistico**: Senza commit nulla identifica il
    codice che ha prodotto i risultati: la regola deve rifiutare, non lasciar
    passare, ed e' il comportamento che distingue una configurazione congelata.
    """
    spoglia = tmp_path / "senza_git"
    (spoglia / "R").mkdir(parents=True)
    monkeypatch.setattr(provenienza, "RADICE_REPOSITORY", spoglia)
    esito = _rifiutata(scenario, rigorosa, "E-PROV-01", capsys)
    assert esito["attiva"] is True and esito["git"] == {"verificato": False, "commit": None}
    assert "ambiente_r" not in esito

    libera = _scrivi_config(scenario, tmp_path / "libera.yaml")
    assert _cli("validate", "--config", str(libera), capsys=capsys) == 0
    assert _controlli(scenario.config.io.out_root)["provenienza_rigorosa"] == {"attiva": False}


def test_con_modifiche_non_committate_l_esecuzione_e_rifiutata(
    scenario, rigorosa, clone, monkeypatch, capsys
):
    """
    **Obiettivo**: Verificare che, con la regola attiva e un file di ``src/``
    modificato rispetto al commit, l'avvio si fermi con ``E-PROV-02`` prima di
    qualunque fase, con il file nel dettaglio del punto di ripresa; e che il
    report dell'esecuzione riporti il rifiuto nella colonna della regola.

    **Razionale scientifico e sistemistico**: Un risultato prodotto da codice
    modificato porterebbe il commit di un codice diverso: la provenienza
    registrata sarebbe falsa proprio dove deve certificare.
    """
    monkeypatch.setenv("AMPLICON16S_R_DIR", str(clone / "R"))
    (clone / "src" / "modulo.py").write_text("x = 2\n", encoding="utf-8")
    esito = _rifiutata(scenario, rigorosa, "E-PROV-02", capsys)
    assert esito["git"]["verificato"] is False
    radice = scenario.config.io.out_root
    assert "src/modulo.py" in _punto(radice)["dettaglio"]

    monkeypatch.delenv("AMPLICON16S_R_DIR")
    assert _cli("report", "--config", str(rigorosa), capsys=capsys) == 0
    righe = _tsv(Path(radice) / CARTELLA_REPORT / CARTELLA_TABELLE / "esecuzioni.tsv")
    assert righe[-1]["regola rigorosa sulla provenienza"].startswith("rifiutata: codice non identificato")


def test_le_discordanze_fra_ambiente_e_file_di_blocco():
    """
    **Obiettivo**: Verificare ``discordanze_ambiente``: nessuna discordanza
    con versioni uguali, anche scritte con separatori diversi; una per la
    versione di R, per un pacchetto assente, per una versione diversa, per una
    correzione assente o diversa da quella del file di blocco.

    **Razionale scientifico e sistemistico**: Versioni diverse di R e dei
    pacchetti di calcolo possono dare risultati numerici diversi; la correzione
    di dada2 cambia la tassonomia nei pareggi e ha lo stesso numero di versione
    in ogni ambiente che la dichiara: va confrontata per contenuto.
    """
    correzione = {"Base": "1.36.0", "File": "c.patch", "SHA256": "ab" * 32}
    blocco = {
        "R": {"Version": "4.5.2"},
        "Packages": {
            "dada2": {"Version": "1.36.0.1", "Patch": correzione},
            "randomForest": {"Version": "4.7-1.2"},
            "vegan": {"Version": "2.7-1"},
        },
    }
    ambiente = {
        "r": "4.5.2",
        "pacchetti": {
            "dada2": {"versione": "1.36.0.1", "correzione": dict(correzione)},
            "randomForest": {"versione": "4.7.1.2", "correzione": None},
            "vegan": {"versione": "2.7.1"},
        },
    }
    assert provenienza.discordanze_ambiente(blocco, ambiente) == []

    diverso = json.loads(json.dumps(ambiente))
    diverso["r"] = "4.3.3"
    diverso["pacchetti"]["dada2"]["correzione"] = None
    diverso["pacchetti"]["vegan"]["versione"] = "2.6-4"
    diverso["pacchetti"]["randomForest"] = None
    trovate = provenienza.discordanze_ambiente(blocco, diverso)
    assert trovate == [
        "R: atteso 4.5.2, trovato 4.3.3",
        "dada2: correzione diversa da quella del file di blocco",
        "randomForest: non installato",
        "vegan: atteso 2.7-1, trovato 2.6-4",
    ]
    altra = json.loads(json.dumps(ambiente))
    altra["pacchetti"]["dada2"]["correzione"]["SHA256"] = "cd" * 32
    assert len(provenienza.discordanze_ambiente(blocco, altra)) == 1


def test_un_ambiente_r_diverso_dal_file_di_blocco_e_rifiutato(
    scenario, tmp_path, monkeypatch, capsys
):
    """
    **Obiettivo**: Verificare che, con git in regola e un file di blocco che
    dichiara una versione di R e pacchetti diversi da quelli installati, l'avvio
    si fermi con ``E-PROV-03`` prima di qualunque fase, con le discordanze lette
    dalle librerie installate; e che un file di blocco introvabile dia lo
    stesso codice.

    **Razionale scientifico e sistemistico**: ``run.container`` e' una
    dichiarazione: cio' che la pipeline puo' verificare da dentro il container
    sono le versioni che i calcoli caricheranno, ed e' su quelle che il rifiuto
    deve basarsi.
    """
    motivo = motivo_pacchetti_r_assenti("jsonlite")
    if motivo is not None:
        pytest.skip(motivo)
    monkeypatch.setattr(executor, "verifica_git", lambda script: (None, "0" * 40))
    blocco = tmp_path / "blocco.lock"
    blocco.write_text(json.dumps({
        "R": {"Version": "0.0.1"},
        "Packages": {
            "dada2": {"Version": "0.0.1"},
            "jsonlite": {"Version": "0.0.1"},
            "pacchettoinesistente": {"Version": "1.0"},
        },
    }), encoding="utf-8")
    config = _scrivi_config(scenario, tmp_path / "c.yaml", strict_provenance=True, lockfile=str(blocco))
    esito = _rifiutata(scenario, config, "E-PROV-03", capsys)
    assert esito["git"] == {"verificato": True, "commit": "0" * 40}
    ambiente = esito["ambiente_r"]
    assert ambiente["verificato"] is False and ambiente["pacchetti_confrontati"] == 3
    assert ambiente["r"] != "0.0.1" and re.fullmatch(r"\d+\.\d+\.\d+", ambiente["r"])
    assert "R: atteso 0.0.1, trovato " + ambiente["r"] in ambiente["discordanze"]
    assert "pacchettoinesistente: non installato" in ambiente["discordanze"]
    assert any(d.startswith("jsonlite: atteso 0.0.1, trovato ") for d in ambiente["discordanze"])

    assente = _scrivi_config(
        scenario, tmp_path / "d.yaml", strict_provenance=True, lockfile=str(tmp_path / "manca.lock")
    )
    chiudi()
    esito = _rifiutata(scenario, assente, "E-PROV-03", capsys)
    assert esito["ambiente_r"] == {"verificato": False}


def test_nell_ambiente_dell_immagine_la_regola_e_soddisfatta(scenario, rigorosa, monkeypatch, capsys):
    """
    **Obiettivo**: Verificare che nell'ambiente R dell'immagine, con git in
    regola, la regola confronti tutti i pacchetti di ``renv.lock`` con quelli
    installati senza discordanze, compresa la correzione di dada2, che la
    validazione si concluda e che il manifesto di S0 registri l'impronta
    rigorosa; e che il report dichiari che cosa e' verificato e che cosa e'
    dichiarato.

    **Razionale scientifico e sistemistico**: E' la condizione in cui si
    produce il riferimento: se la verifica fallisse nell'immagine pubblicata,
    il file di blocco non descriverebbe l'ambiente che dichiara.
    """
    motivo = motivo_pacchetti_r_assenti("dada2", "phyloseq", "DESeq2")
    if motivo is not None:
        pytest.skip(motivo)
    monkeypatch.setattr(executor, "verifica_git", lambda script: (None, "0" * 40))
    assert _cli("validate", "--config", str(rigorosa), capsys=capsys) == 0
    radice = scenario.config.io.out_root
    esito = _controlli(radice)["provenienza_rigorosa"]
    attesi = json.loads((RADICE / "renv.lock").read_text(encoding="utf-8"))
    assert esito["ambiente_r"] == {
        "verificato": True, "r": attesi["R"]["Version"],
        "pacchetti_confrontati": len(attesi["Packages"]), "discordanze": [],
    }
    assert esito["immagine_dichiarata"] == scenario.config.run.container
    manifesto = json.loads(next(Path(radice).rglob("manifest_S0.json")).read_text(encoding="utf-8"))
    assert set(manifesto["calcolata_su"]["provenienza"]) == {"sorgente", "blocco_r", "immagine"}

    assert _cli("report", "--config", str(rigorosa), capsys=capsys) == 0
    righe = _tsv(Path(radice) / CARTELLA_REPORT / CARTELLA_TABELLE / "esecuzioni.tsv")
    testo = righe[-1]["regola rigorosa sulla provenienza"]
    assert testo.startswith("verificati il codice (commit 000000000000") and "dichiarata" in testo


def test_con_la_regola_un_immagine_dichiarata_diversa_rifa_la_fase(
    scenario, rigorosa, tmp_path, monkeypatch, capsys
):
    """
    **Obiettivo**: Verificare che, conclusa S0 con la regola attiva, dichiarare
    un'immagine diversa renda S0 da eseguire con il motivo della regola
    rigorosa; e che senza la regola la stessa modifica lasci S0 conclusa.

    **Razionale scientifico e sistemistico**: E' la differenza fra avviso e
    vincolo: nella configurazione congelata un risultato non puo' restare
    valido sotto una dichiarazione d'ambiente diversa da quella con cui e'
    stato calcolato.
    """
    monkeypatch.setattr(Esecutore, "_provenienza_rigorosa", lambda self: ({"attiva": True}, None))
    assert _cli("validate", "--config", str(rigorosa), capsys=capsys) == 0
    fasi = {Passo.S0: ValidazioneIngressi()}
    stessa = ProjectRun(carica(rigorosa), passi=fasi).valuta().situazioni[Passo.S0]
    assert stessa.stato is StatoPasso.COMPLETATA

    altra = _scrivi_config(scenario, tmp_path / "altra.yaml", strict_provenance=True)
    dati = yaml.safe_load(altra.read_text(encoding="utf-8"))
    dati["run"]["container"] = IMMAGINE_DIVERSA
    altra.write_text(yaml.safe_dump(dati, sort_keys=False), encoding="utf-8")
    situazione = ProjectRun(carica(altra), passi=fasi).valuta().situazioni[Passo.S0]
    assert situazione.stato is StatoPasso.DA_ESEGUIRE
    assert "regola rigorosa sulla provenienza" in situazione.motivo


def test_senza_la_regola_un_immagine_dichiarata_diversa_non_rifa_la_fase(
    scenario, tmp_path, capsys
):
    """
    **Obiettivo**: Verificare che, conclusa S0 senza la regola, dichiarare
    un'immagine diversa lasci S0 conclusa, e che attivare la regola su
    un'esecuzione calcolata senza renda S0 da eseguire.

    **Razionale scientifico e sistemistico**: Un risultato calcolato senza le
    verifiche della regola non puo' essere promosso a risultato certificato
    attivandola dopo: va ricalcolato sotto la regola.
    """
    libera = _scrivi_config(scenario, tmp_path / "libera.yaml")
    assert _cli("validate", "--config", str(libera), capsys=capsys) == 0
    fasi = {Passo.S0: ValidazioneIngressi()}
    dati = yaml.safe_load(libera.read_text(encoding="utf-8"))
    dati["run"]["container"] = IMMAGINE_DIVERSA
    libera.write_text(yaml.safe_dump(dati, sort_keys=False), encoding="utf-8")
    assert ProjectRun(carica(libera), passi=fasi).valuta().situazioni[Passo.S0].stato is StatoPasso.COMPLETATA

    dati["run"]["strict_provenance"] = True
    libera.write_text(yaml.safe_dump(dati, sort_keys=False), encoding="utf-8")
    situazione = ProjectRun(carica(libera), passi=fasi).valuta().situazioni[Passo.S0]
    assert situazione.stato is StatoPasso.DA_ESEGUIRE
    assert "regola rigorosa sulla provenienza" in situazione.motivo


def test_i_codici_della_regola_dicono_che_cosa_fare():
    """
    **Obiettivo**: Verificare che i tre codici ``E-PROV-*`` siano nel catalogo
    con fase ``PROV`` e revisione umana, e che l'azione di ciascuno nomini il
    parametro ``run.strict_provenance`` e il rimedio (clone del repository,
    commit delle modifiche, immagine indicata).

    **Razionale scientifico e sistemistico**: Un rifiuto all'avvio e' utile
    solo se chi lo riceve capisce che cosa manca; nessuno dei tre casi si
    risolve ritentando.
    """
    attesi = {"E-PROV-01": "clone", "E-PROV-02": "git status", "E-PROV-03": "immagine"}
    for codice, rimedio in attesi.items():
        voce = CATALOGO[codice]
        assert voce.fase == "PROV" and voce.categoria is Categoria.REVISIONE_UMANA
        assert "run.strict_provenance" in voce.azione and rimedio in voce.azione, codice
