r"""Suite di test della settimana 30: indipendenza dei test dal dataset di riferimento.

1. Inquadramento nel Piano Operativo
------------------------------------
Settimana 30 (W30), Fase F8 (generalita'): i test non dipendono dai valori del
dataset di riferimento, che si verificano in un solo punto.

2. Moduli sorgente coperti
--------------------------
* ``src/amplicon16s/config/defaults.py`` (``OBBLIGATORI``, ``STANDARD_DEL_METODO``)
* ``dati/osd734/config_osd734.yaml``, ``config/config.example.yaml``
* ``tests/conftest.py`` (``Formato``, ``crea_scenario``, catene condivise)
* ``tests/sottoinsieme.py``, ``tests/fixtures/osd734/`` (configurazione della
  versione ridotta, selezione, valori attesi del dataset completo)
* ``pyproject.toml`` (marcatori)

3. Cosa valuta questo file
--------------------------
- **il solo test sui valori del dataset di riferimento**: la configurazione
  pubblicata di OSD-734 riporta, per ogni
  parametro derivato dal dataset, il valore registrato nell'elenco;
- il generatore di scenari non contiene etichette, colonne ne' nomi del dataset
  di riferimento, e ogni scenario dichiara tutti i parametri che descrivono il
  dataset;
- un formato diverso da quello predefinito attraversa la validazione;
- la versione ridotta ha la propria configurazione, che non dipende dai
  predefiniti ne' dall'esempio;
- nessun modulo di test avvia R all'importazione, e l'attesa di una catena
  condivisa ha un limite;
- i test sui dati reali non stampano misure, e i numeri del dataset completo
  stanno nel solo file dei valori attesi;
- il confronto con i checksum pubblicati ha un marcatore proprio.

4. Comandi Bash e scenari di esecuzione
---------------------------------------
    ``<immagine>`` e' l'immagine del container della pipeline; quella corrente
    e' indicata in ``README.md``.

    1. Modalità locale standard (R di base con jsonlite, senza Bioconductor né
       dati reali):
       pytest tests/test_w30_indipendenza.py -v

    2. Modalità container Docker standard (sottoinsieme ridotto con
       R/Bioconductor):
       docker run --rm \
         -e PYTHONPATH=/app/src \
         -e AMPLICON16S_R_DIR=/app/R \
         -v "$(pwd)":/app \
         -w /app \
         <immagine> \
         pytest -o cache_dir=/tmp/.pytest_cache tests/test_w30_indipendenza.py -v

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
         pytest -o cache_dir=/tmp/.pytest_cache tests/test_w30_indipendenza.py -v

5. Risultato atteso
-------------------
I conteggi li da' pytest (``pytest --collect-only -q``).

6. Razionale scientifico e sistemistico
---------------------------------------
- Un test che verifica un valore del dataset di riferimento invece di un
  comportamento cade appena quel valore viene corretto, e non dice nulla su un
  altro dataset: i valori si verificano una volta, dove sono dichiarati.
- Un generatore di scenari con le etichette del dataset di riferimento prova la
  pipeline solo su metadati fatti come i suoi.
- Una suite che avvia R per decidere che cosa saltare paga quel costo anche
  per la sola raccolta dei test, e un'attesa senza limite trasforma un
  processo morto in una suite ferma.
"""

from __future__ import annotations

import ast
import dataclasses
import re
from pathlib import Path

import pytest
import yaml
from conftest import (
    ATTESA_MASSIMA_S,
    FORMATO,
    MARCATORI_SUL_DATASET_COMPLETO,
    Campione,
    Formato,
    crea_scenario,
    esegui_gate_metadati,
    parametri_del_formato,
)
from sottoinsieme import (
    PARAMETRI,
    VALORI_ATTESI,
    attesi_dataset,
    config_ridotta,
)

from amplicon16s.config import defaults
from amplicon16s.config.resolve import parametri_dichiarati
from amplicon16s.config.schema import carica
from amplicon16s.logging.logger import chiudi
from amplicon16s.metadata.models import ClasseCampione
from amplicon16s.steps.s00_validate import esegui_s0

RADICE = Path(__file__).resolve().parents[1]
TEST = RADICE / "tests"
CONFIG_OSD734 = RADICE / "dati" / "osd734" / "config_osd734.yaml"


@pytest.fixture(autouse=True)
def uscite_pulite():
    """Chiude le uscite del log prima e dopo ogni test, perché nessun handler resti
    aperto sulla cartella temporanea.
    """
    chiudi()
    yield
    chiudi()


def _parametri_di_osd734() -> dict[str, object]:
    """I valori dei parametri del dataset di riferimento, dal file dei valori attesi."""
    return attesi_dataset()["parametri"]


# --------------------------------------------------------------------------- #
# 1. Il solo test sui valori del dataset di riferimento                        #
# --------------------------------------------------------------------------- #


@pytest.mark.skipif(not CONFIG_OSD734.is_file(),
                    reason="la cartella dati/ non e' presente (nell'immagine non viene copiata)")
def test_i_valori_di_osd734_sono_quelli_dell_elenco_dei_parametri_derivati_dal_dataset():
    """
    **Obiettivo**: Verificare che ``dati/osd734/config_osd734.yaml`` dichiari,
    per ogni parametro dell'elenco dei valori del dataset di riferimento
    (``tests/fixtures/osd734/valori_attesi.yaml``, voce ``parametri``), proprio
    quel valore; e che l'elenco comprenda tutti i parametri obbligatori.

    **Razionale scientifico e sistemistico**: E' l'unico test che conosce i
    valori del dataset di riferimento. Se un valore dell'elenco cambia senza
    che cambi la configurazione pubblicata, o viceversa, fallisce questo test e
    nessun altro: gli altri verificano comportamenti, su dati che dichiarano i
    propri parametri, e non devono accorgersene.
    """
    pubblicata = carica(CONFIG_OSD734)
    dichiarati = set(parametri_dichiarati(pubblicata))
    elenco = _parametri_di_osd734()
    assert set(defaults.OBBLIGATORI) <= set(elenco)
    for chiave, atteso in elenco.items():
        gruppo, nome = chiave.split(".")
        assert chiave in dichiarati, f"{chiave} non e' dichiarato nella configurazione pubblicata"
        assert getattr(getattr(pubblicata, gruppo), nome) == atteso, chiave


def test_nessun_valore_del_dataset_di_riferimento_fa_da_predefinito():
    """
    **Obiettivo**: Verificare che nessun parametro dell'elenco dei valori del
    dataset di riferimento abbia un predefinito nello schema, salvo quelli il
    cui predefinito e' il valore standard del metodo, con la fonte; e che
    nessuna delle stringhe proprie del dataset compaia nei predefiniti, nello
    schema o nel modello di configurazione.

    **Razionale scientifico e sistemistico**: Un predefinito preso dal dataset
    su cui lo strumento e' nato verrebbe ereditato in silenzio da ogni altro
    dataset: cio' che dipende dai dati si dichiara, e resta scritto solo nella
    configurazione del dataset.
    """
    from amplicon16s.config.schema import Config

    for chiave in _parametri_di_osd734():
        gruppo, nome = chiave.split(".")
        campo = Config.model_fields[gruppo].annotation.model_fields[nome]
        assert campo.is_required() or chiave in defaults.STANDARD_DEL_METODO, chiave
    radice = TEST.parent
    sorgenti = [radice / "src" / "amplicon16s" / "config" / nome
                for nome in ("defaults.py", "schema.py")]
    sorgenti.append(radice / "config" / "config.example.yaml")
    for percorso in sorgenti:
        testo = percorso.read_text(encoding="utf-8")
        # Come valore: un letterale fra virgolette nel codice, il valore di una
        # chiave o la voce di un elenco nel modello. Una parola comune dentro
        # un commento non e' un valore del dataset.
        presenti = sorted(
            s for s in _stringhe_del_dataset_di_riferimento()
            if f'"{s}"' in testo or f"'{s}'" in testo
            or re.search(rf"^\s*(#\s*)?(- |[A-Za-z_]+:\s*){re.escape(s)}\s*(#.*)?$", testo, flags=re.M)
        )
        assert presenti == [], percorso.name
        assert not re.search(r"OSD-?[0-9]{3}", testo), percorso.name


# --------------------------------------------------------------------------- #
# 2. Il generatore di scenari                                                  #
# --------------------------------------------------------------------------- #


def _stringhe_del_dataset_di_riferimento() -> set[str]:
    """Etichette, colonne e nomi propri del dataset di riferimento: i valori
    testuali dei parametri che lo descrivono, tranne quelli che sono
    convenzioni generali (un'espressione regolare, il primer, il motivo, il
    modello dei nomi dei file) o voci di un vocabolario chiuso dello schema.
    """
    generali = {"io.accession_regex", "io.fastq_glob", "meta.module_regex",
                "qc.primer_sequence", "qc.conserved_motif",
                "dada.pool", "katharoseq.collapse_rank", "decontam.mode",
                "decontam.batch_combine", "qc.min_reads_mode"}
    stringhe: set[str] = set()
    for chiave, valore in _parametri_di_osd734().items():
        if chiave in generali:
            continue
        valori = valore if isinstance(valore, list) else [valore]
        stringhe |= {v for v in valori if isinstance(v, str) and len(v) > 3}
    return stringhe


def test_il_generatore_di_scenari_non_contiene_nomi_del_dataset_di_riferimento():
    """
    **Obiettivo**: Verificare che nel sorgente del generatore di scenari
    (``tests/conftest.py``) non compaia alcuna etichetta, colonna o nome del
    dataset di riferimento, ne' il suo identificativo; che il generatore non
    legga i valori d'esempio del dataset; e che il formato predefinito non
    condivida con il dataset di riferimento alcun nome di colonna, etichetta
    di classe o taxon atteso.

    **Razionale scientifico e sistemistico**: Uno scenario fatto con i nomi
    del dataset di riferimento prova la pipeline solo su metadati come i suoi,
    e nasconde ogni punto del codice che quei nomi li presuppone.
    """
    sorgente = (TEST / "conftest.py").read_text(encoding="utf-8")
    # Come letterali fra virgolette: una parola comune dentro un commento non
    # e' un nome del dataset.
    presenti = sorted(s for s in _stringhe_del_dataset_di_riferimento()
                      if f'"{s}"' in sorgente or f"'{s}'" in sorgente)
    assert presenti == []
    assert "valori_attesi" not in sorgente
    assert not re.search(r"OSD-?734|osd734", sorgente)
    del_formato: set[str] = set()
    for campo in dataclasses.fields(FORMATO):
        valore = getattr(FORMATO, campo.name)
        if isinstance(valore, str):
            del_formato.add(valore)
        elif isinstance(valore, tuple):
            del_formato |= {v for v in valore if isinstance(v, str)}
    assert not del_formato & _stringhe_del_dataset_di_riferimento()


def test_ogni_scenario_dichiara_i_parametri_che_descrivono_il_dataset(tmp_path):
    """
    **Obiettivo**: Verificare che la configurazione di uno scenario dichiari
    tutti i parametri obbligatori, comprese le scelte di analisi e le soglie,
    e che i valori dichiarati siano quelli del formato dello scenario.

    **Razionale scientifico e sistemistico**: Un parametro ereditato cambia
    quando cambia il predefinito: il test che lo usa cadrebbe, o peggio
    cambierebbe significato, per una modifica che non lo riguarda.
    """
    dati = parametri_del_formato()
    dichiarati = {f"{g}.{n}" for g, valori in dati.items() for n in valori}
    assert set(defaults.OBBLIGATORI) <= dichiarati
    scenario = crea_scenario(tmp_path, [Campione("ERX1000001", "MOD1A2.S1")])
    config = scenario.config
    assert config.meta.sample_id_column == FORMATO.colonna_campione
    assert config.ctrl.biological_values == [FORMATO.biologico]
    assert config.filter.truncLen == FORMATO.troncamento
    for chiave, valore in FORMATO.scelte:
        gruppo, nome = chiave.split(".")
        atteso = list(valore) if isinstance(valore, tuple) else valore
        assert getattr(getattr(config, gruppo), nome) == atteso, chiave
    assert {c for c, _ in FORMATO.scelte} <= set(defaults.OBBLIGATORI) | {"tax.min_boot"}


def test_un_formato_diverso_attraversa_la_validazione(tmp_path):
    """
    **Obiettivo**: Verificare che uno scenario costruito con un formato
    diverso da quello predefinito in ogni nome (colonne, etichette delle
    classi, file, colonne del lotto, troncamento) superi tutti i gate di S0, e
    che l'inventario abbia le classi, la piastra e la corsa attese per ogni
    campione, compreso quello che la regola del formato riclassifica.

    **Razionale scientifico e sistemistico**: Il formato e' un parametro del
    generatore solo se cambiarlo per intero non rompe nulla: e' la prova che
    ne' il generatore ne' la validazione presuppongono i nomi di un dataset.
    """
    formato = dataclasses.replace(
        FORMATO,
        colonna_campione="ID", colonna_file="FASTQ", colonna_classe="Kind",
        colonna_posizione="Site", colonna_cellule="Cells",
        biologico="sample", positivo="mock", negativo="blank",
        non_superfici=("air",), riclassificati=("sealed swab",),
        modello_file="{accession}.fastq.gz", modello_file_assay="{accession}.fq.gz",
        colonna_chiave_lotto="run", colonna_piastra="plate", colonna_corsa="lane",
        colonna_modulo_lotto="area", troncamento=120, taxon_atteso="Altro_genere",
    )
    campioni = [
        Campione("ERX7000001", "a1", materiale="sample", posizione="LAB2B1", piastra="P1"),
        Campione("ERX7000002", "a2", materiale="sample", posizione="sealed swab", piastra="P1"),
        Campione("ERX7000003", "p1", materiale="mock", posizione="n/a", piastra="P1"),
        Campione("ERX7000004", "n1", materiale="blank", posizione="n/a", piastra="P2",
                 corsa="corsa_B"),
    ]
    scenario = crea_scenario(tmp_path, campioni, formato=formato, con_arricchimento=True,
                             con_letture=True)
    assert scenario.formato is formato
    intestazione = Path(scenario.config.io.study_table).read_text(encoding="utf-8").splitlines()[0]
    assert intestazione.split("\t") == ["ID", "Kind", "Site", "Cells"]
    risultato = esegui_s0(scenario.config, solleva=False)
    assert risultato.superata, [(e.gate, [str(v) for v in e.violazioni]) for e in risultato.falliti]
    inventario = esegui_gate_metadati(scenario.config)
    assert {c.nome: c.classe for c in inventario} == {
        "a1": ClasseCampione.BIOLOGICO, "a2": ClasseCampione.CONTROLLO_NEGATIVO,
        "p1": ClasseCampione.CONTROLLO_POSITIVO, "n1": ClasseCampione.CONTROLLO_NEGATIVO,
    }
    assert inventario["ERX7000004"].piastra == "P2" and inventario["ERX7000004"].corsa == "corsa_B"
    assert inventario["ERX7000001"].file.name == "ERX7000001.fastq.gz"


# --------------------------------------------------------------------------- #
# 3. La versione ridotta                                                       #
# --------------------------------------------------------------------------- #


def test_la_versione_ridotta_ha_la_propria_configurazione(tmp_path):
    """
    **Obiettivo**: Verificare che il file dei parametri accanto alla versione
    ridotta dichiari tutti i parametri obbligatori; che la configurazione
    della versione ridotta abbia quei valori; e che il modulo che la
    costruisce non legga il modello di configurazione ne' la configurazione
    pubblicata del dataset di riferimento.

    **Razionale scientifico e sistemistico**: La versione ridotta e' un
    dataset con la sua configurazione: i test che la attraversano verificano
    il comportamento della pipeline su quei dati, e non devono cambiare esito
    quando cambia un altro file.
    """
    propri = yaml.safe_load(PARAMETRI.read_text(encoding="utf-8"))
    dichiarati = {f"{g}.{n}" for g, valori in propri.items() for n in valori}
    assert set(defaults.OBBLIGATORI) <= dichiarati
    config = config_ridotta(tmp_path)
    for gruppo, valori in propri.items():
        for nome, valore in valori.items():
            if (gruppo, nome) == ("katharoseq", "target_taxon"):
                continue  # nel riferimento sintetico il taxon ha un nome inventato
            assert getattr(getattr(config, gruppo), nome) == valore, f"{gruppo}.{nome}"
    sorgente = (TEST / "sottoinsieme.py").read_text(encoding="utf-8")
    assert "config.example" not in sorgente and "config_osd734" not in sorgente


# --------------------------------------------------------------------------- #
# 4. L'infrastruttura dei test                                                 #
# --------------------------------------------------------------------------- #


def _moduli_di_test() -> list[Path]:
    """I moduli di test e i due di infrastruttura."""
    return sorted(TEST.glob("test_*.py")) + [TEST / "conftest.py", TEST / "sottoinsieme.py"]


def test_nessun_modulo_di_test_avvia_r_all_importazione():
    """
    **Obiettivo**: Verificare, sull'albero sintattico di ogni modulo di test,
    che nessuna istruzione eseguita all'importazione (assegnazioni, decoratori
    e argomenti predefiniti a livello di modulo) chiami la sonda dei pacchetti
    R, la ricerca dell'interprete o un sottoprocesso; e che la sonda tenga in
    memoria il suo esito.

    **Razionale scientifico e sistemistico**: La sola raccolta dei test, che
    serve anche a contarli e a generare la documentazione, non deve avviare
    decine di processi R: le sonde partono dalle fixture, al primo test che ne
    ha bisogno.
    """
    vietate = {"motivo_pacchetti_r_assenti", "trova_rscript", "run", "check_output", "Popen"}

    def chiamate_all_importazione(albero: ast.Module):
        for nodo in albero.body:
            if isinstance(nodo, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
                eseguiti = list(nodo.decorator_list)
                if isinstance(nodo, ast.FunctionDef):
                    eseguiti += nodo.args.defaults + [d for d in nodo.args.kw_defaults if d]
            else:
                eseguiti = [nodo]
            for radice in eseguiti:
                for interno in ast.walk(radice):
                    if isinstance(interno, ast.Call):
                        funzione = interno.func
                        nome = funzione.id if isinstance(funzione, ast.Name) else getattr(
                            funzione, "attr", "")
                        if nome in vietate or nome.startswith(("_sonda_", "_motivo_")):
                            yield nome, interno.lineno

    trovate = {
        percorso.name: sorted(chiamate_all_importazione(ast.parse(percorso.read_text("utf-8"))))
        for percorso in _moduli_di_test()
    }
    assert {nome: voci for nome, voci in trovate.items() if voci} == {}
    from sottoinsieme import motivo_pacchetti_r_assenti

    assert hasattr(motivo_pacchetti_r_assenti, "cache_info")


def test_l_attesa_della_catena_condivisa_ha_un_limite(tmp_path, monkeypatch):
    """
    **Obiettivo**: Verificare che un processo in attesa di una catena condivisa
    che un altro sta calcolando, e che non arriva mai, si fermi allo scadere
    del limite con un errore che nomina la catena e la variabile che cambia il
    limite; e che il limite predefinito sia finito.

    **Razionale scientifico e sistemistico**: Se il processo che calcola la
    catena viene ucciso dal sistema non lascia traccia: senza un limite gli
    altri resterebbero in attesa per sempre, e la suite non finirebbe mai
    invece di fallire.
    """
    import conftest

    assert 0 < ATTESA_MASSIMA_S < 24 * 3600

    class _Fabbrica:
        def getbasetemp(self):
            return tmp_path / "base"

    (tmp_path / "condivise").mkdir()
    (tmp_path / "condivise" / "mai.blocco").touch()
    monkeypatch.setenv("PYTEST_XDIST_WORKER", "gw9")
    monkeypatch.setattr(conftest, "ATTESA_MASSIMA_S", 1)
    with pytest.raises(RuntimeError) as info:
        conftest._condivisa(_Fabbrica(), "mai", lambda cartella: None)
    assert "mai" in str(info.value) and conftest.VARIABILE_ATTESA in str(info.value)


def test_i_test_non_stampano_misure_e_i_numeri_del_dataset_stanno_in_un_file():
    """
    **Obiettivo**: Verificare che nessun modulo di test chiami ``print``; che
    i fatti del dataset completo si leggano dal solo file dei valori attesi,
    con le tre classi che sommano ai campioni e le piastre che li contengono
    tutti; e che i numeri di campioni e di classi del dataset completo non
    compaiano come letterali in alcun altro modulo di test.

    **Razionale scientifico e sistemistico**: Una misura stampata non
    verifica nulla: passa qualunque valore assuma. E un numero del dataset
    scritto in dieci moduli va corretto in dieci posti, e in nessuno si sa da
    dove venga.
    """
    con_stampe = {}
    for percorso in _moduli_di_test():
        albero = ast.parse(percorso.read_text(encoding="utf-8"))
        righe = [n.lineno for n in ast.walk(albero)
                 if isinstance(n, ast.Call) and isinstance(n.func, ast.Name) and n.func.id == "print"]
        if righe:
            con_stampe[percorso.name] = righe
    assert con_stampe == {}

    attesi = attesi_dataset()
    assert sum(attesi["classi"].values()) == attesi["campioni"]
    assert set(attesi["classi"]) == {c.value for c in ClasseCampione}
    lotto = attesi["lotto"]
    assert lotto["piastre"] * lotto["campioni_per_piastra"] == attesi["campioni"]
    assert lotto["corse"] * lotto["campioni_per_corsa"] == attesi["campioni"]
    assert VALORI_ATTESI.parent == PARAMETRI.parent.parent

    # I conteggi propri del dataset completo: campioni, classi, riclassificati.
    numeri = {attesi["campioni"], *attesi["classi"].values()} - {80}
    altrove = {}
    for percorso in sorted(TEST.glob("test_*.py")):
        albero = ast.parse(percorso.read_text(encoding="utf-8"))
        trovati = sorted({n.value for n in ast.walk(albero)
                          if isinstance(n, ast.Constant) and type(n.value) is int
                          and n.value in numeri})
        if trovati:
            altrove[percorso.name] = trovati
    assert altrove == {}


def test_il_confronto_con_i_checksum_ha_un_marcatore_proprio():
    """
    **Obiettivo**: Verificare che ``pyproject.toml`` registri il marcatore
    ``riferimento`` accanto a ``dati_reali``; che il test che confronta la
    catena sul dataset completo con i checksum pubblicati porti il primo e non
    il secondo; che nessun altro test lo porti; e che i due marcatori
    condividano il gruppo che li tiene su un solo processo.

    **Razionale scientifico e sistemistico**: Il confronto con i checksum non
    dice se un comportamento e' corretto, dice se i risultati sono quelli
    pubblicati: fallisce a ragione dopo ogni modifica voluta dei risultati.
    Mescolato ai test di correttezza, un suo fallimento atteso li
    nasconderebbe, e viceversa.
    """
    progetto = (RADICE / "pyproject.toml").read_text(encoding="utf-8")
    assert re.search(r'^\s+"riferimento: ', progetto, flags=re.M)
    assert re.search(r'^\s+"dati_reali: ', progetto, flags=re.M)
    assert MARCATORI_SUL_DATASET_COMPLETO == ("dati_reali", "riferimento")
    marcati = {}
    for percorso in sorted(TEST.glob("test_*.py")):
        albero = ast.parse(percorso.read_text(encoding="utf-8"))
        for nodo in albero.body:
            if isinstance(nodo, ast.FunctionDef) and nodo.name.startswith("test_"):
                nomi = {d.attr for d in nodo.decorator_list if isinstance(d, ast.Attribute)}
                if "riferimento" in nomi:
                    marcati[nodo.name] = (percorso.name, "dati_reali" in nomi)
    assert marcati == {
        "test_la_catena_sul_dataset_completo_da_i_checksum_pubblicati":
            ("test_w24_riproduzione.py", False),
    }


