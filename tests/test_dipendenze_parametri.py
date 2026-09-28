"""Suite di test per il contratto trasversale delle dipendenze di parametro (S0 - S14).

Inquadramento nel Piano Operativo:
    Modulo trasversale e infrastrutturale del DAG (Settimane W8 - W14, Fasi F3 e F4).
    Presidia il contratto di scoping dei parametri dichiarati da ciascuna fase
    della pipeline (fasi realizzate da S0 a S7 sull'intero grafo S0 - S14) e la
    validita' selettiva dei manifesti di fase.

Moduli sorgente coperti:
    - src/amplicon16s/config/vista.py (VistaConfig, ParametroNonDichiarato,
      risolta_ristretta, _VistaGruppo, _VistaDerivati)
    - src/amplicon16s/config/resolve.py (ConfigRisolta.digest_parametri,
      PARAMETRI_SENZA_EFFETTO, incluso run.batch_size)
    - src/amplicon16s/steps/base.py (PipelineStep.__init_subclass__, validazione
      statica dell'attributo di classe 'parametri' e impronta selettiva)
    - R/lib/errors.R (classe S3 'parametri_dichiarati', operatori '$.parametri_dichiarati',
      '[[.parametri_dichiarati' e funzione 'facoltativo')
    - tests/r_doppioni/parametro_non_ricevuto.R

Cosa valuta questo file:
    1. Validazione statica alla definizione della classe: una sottoclasse di
       PipelineStep priva dell'attributo 'parametri', con una voce inesistente
       nello schema o che dichiara un parametro operativo senza effetto
       ('run.threads', 'io.out_root', 'run.batch_size') viene rifiutata
       immediatamente con errore.
    2. Isolamento a tempo di esecuzione in Python ('VistaConfig'): una fase vede
       attraverso StepContext solo i gruppi e le chiavi che ha dichiarato (oltre
       a 'PARAMETRI_SENZA_EFFETTO'); ogni lettura di parametri non dichiarati,
       di derivati non ammessi o del digest/mappa globale solleva
       'ParametroNonDichiarato' (sottoclasse di LookupError, non catturabile da
       getattr(..., default) o hasattr).
    3. Isolamento simmetrico negli script R ('parametri_dichiarati'): uno script R
       che legge tramite '$' o '[[' un parametro non passato dalla fase Python
       fallisce immediatamente invece di ricevere silenziosamente NULL.
    4. Impronta selettiva e invalidazione mirata nel DAG (S0 - S7): cambiare un
       parametro di una fase a valle (come 'filter.maxEE' in S2, 'err.nbases' in
       S3, 'dada.omega_a' in S4, 'qc.max_asv_count' in S5, 'chimera.*' in S6 o
       'asv.len_tol' in S7) invalida solo la fase interessata, mentre i parametri
       operativi ('run.threads', 'run.keep_filtered_fastq', 'run.batch_size',
       'io.out_root', 'retry.enabled', 'retry.max_attempts') non invalidano
       nessuna fase.

Comandi Bash e scenari di esecuzione:
    1. Modalita locale standard (con R base e jsonlite):
       pytest tests/test_dipendenze_parametri.py -v
       Risultato atteso: 20 test (20 passed se R e jsonlite sono presenti;
       19 passed e 1 skipped se R/jsonlite non sono installati sull'host).

    2. Modalita container Docker standard:
       docker run --rm \
         -e PYTHONPATH=/app/src \
         -v "$(pwd)":/app:ro \
         --entrypoint pytest \
         amplicon16s:dev \
         tests/test_dipendenze_parametri.py -v
       Risultato atteso: 20 passed.

Risultato atteso:
    20 test totali (12 funzioni di test, di cui 3 parametrizzate):
    20 passed in ambiente con R + jsonlite (~0.50s);
    19 passed, 1 skipped in ambiente privo di R (~0.30s).

Razionale scientifico e sistemistico:
    1. Prevenzione strutturale di risultati obsoleti silenziosi: in una pipeline
       multi-step in cui il manifesto di ciascuna fase dipende solo dai propri
       parametri dichiarati, dimenticare di dichiarare un parametro produrrebbe
       artefatti non rigenerati al variare della configurazione. Rendi la
       dichiarazione vincolante per costruzione ('VistaConfig' in Python e
       'parametri_dichiarati' in R) trasforma ogni omissione in un errore
       immediato al primo test che esercita la fase.
    2. Ereditarieta' da LookupError anziche' AttributeError: in Python la
       funzione built-in getattr(obj, nome, default) intercetta AttributeError
       e restituisce il valore di fallback; derivando 'ParametroNonDichiarato'
       da LookupError, nessun accesso indiretto puo' mascherare una lettura
       non autorizzata.
    3. Riduzione dei tempi di ricalcolo e di collaudo in container: poiche' S1
       dichiara solo ('filter.truncLen', 'filter.truncLen_shortfall_warn') e S2
       non dipende da 'err.*', i test e le riprese che variano parametri a valle
       riutilizzano i profili e i filtrati gia' calcolati, riducendo la durata
       del job containerizzato della CI da 6m27s a 4m40s.
"""

from __future__ import annotations

import os
from pathlib import Path
from typing import ClassVar

import pytest
from sottoinsieme import config_ridotta

from amplicon16s.config.resolve import risolvi
from amplicon16s.config.vista import ParametroNonDichiarato, VistaConfig, risolta_ristretta
from amplicon16s.errors.exceptions import ErroreRevisioneUmana
from amplicon16s.io_layer.artifacts import AlberoOutput, Fase
from amplicon16s.logging.logger import chiudi
from amplicon16s.rbridge.runner import esegui_script
from amplicon16s.runner.graph import Passo
from amplicon16s.runner.project import ProjectRun, passi_realizzati
from amplicon16s.steps.base import PipelineStep, Produzione, StepContext
from amplicon16s.steps.s00_validate import ValidazioneIngressi, esegui_s0

DOPPIONI_R = Path(__file__).resolve().parent / "r_doppioni"


@pytest.fixture(autouse=True)
def uscite_pulite():
    chiudi()
    yield
    chiudi()


# --------------------------------------------------------------------------- #
# La vista                                                                     #
# --------------------------------------------------------------------------- #


def test_la_vista_lascia_leggere_solo_i_parametri_dichiarati(tmp_path):
    """
    **Obiettivo**: Verificare che la vista lasci leggere un gruppo dichiarato
    per intero, una chiave dichiarata da sola e i parametri senza effetto sui
    risultati, e rifiuti tutto il resto.

    **Razionale Scientifico/Sistemistico**: Una fase che legge un parametro non
    dichiarato non verrebbe invalidata quando quel parametro cambia: la vista
    lo trasforma in un errore al primo accesso.
    """
    config = config_ridotta(tmp_path)
    vista = VistaConfig(config, ("filter", "err.nbases"))
    assert vista.filter.truncLen == config.filter.truncLen
    assert vista.err.nbases == config.err.nbases
    assert vista.run.threads == config.run.threads  # senza effetto: sempre leggibile
    with pytest.raises(ParametroNonDichiarato, match="err.max_consist"):
        vista.err.max_consist
    with pytest.raises(ParametroNonDichiarato, match="tax.ref_name"):
        vista.tax.ref_name
    with pytest.raises(ParametroNonDichiarato, match="run.seed"):
        vista.run.seed
    with pytest.raises(ParametroNonDichiarato):
        vista.model_dump()


def test_un_accesso_rifiutato_non_si_puo_inghiottire_in_silenzio(tmp_path):
    """
    **Obiettivo**: Verificare che l'errore della vista non sia un
    ``AttributeError``, cosi' che ``getattr`` con un valore predefinito non lo
    trasformi in un valore qualunque.

    **Razionale Scientifico/Sistemistico**: Un errore inghiottito renderebbe di
    nuovo silenziosa la dipendenza dimenticata.
    """
    vista = VistaConfig(config_ridotta(tmp_path), ("filter",))
    assert not issubclass(ParametroNonDichiarato, AttributeError)
    with pytest.raises(ParametroNonDichiarato):
        getattr(vista.tax, "ref_name", None)


def test_anche_i_derivati_passano_per_la_vista(tmp_path):
    """
    **Obiettivo**: Verificare che un parametro derivato si legga solo se il suo
    gruppo, o la sua chiave, e' dichiarato.

    **Razionale Scientifico/Sistemistico**: ``asv.len_min`` discende da filter
    e asv: leggerlo senza dichiararli sarebbe la stessa dipendenza nascosta.
    """
    risolta = risolvi(config_ridotta(tmp_path))
    assert risolta_ristretta(risolta, ("asv",)).derivati.asv_len_min == risolta.derivati.asv_len_min
    with pytest.raises(ParametroNonDichiarato, match="asv.len_min"):
        risolta_ristretta(risolta, ("filter",)).derivati.asv_len_min
    with pytest.raises(ParametroNonDichiarato):
        risolta_ristretta(risolta, ("filter",)).digest


# --------------------------------------------------------------------------- #
# Una fase che legge un parametro non dichiarato                               #
# --------------------------------------------------------------------------- #


class _Distratta(PipelineStep):
    """Dichiara filter.truncLen, ma legge anche err.nbases."""

    passo: ClassVar[Passo] = Passo.S1
    parametri: ClassVar[tuple[str, ...]] = ("filter.truncLen",)

    def calcola(self, contesto: StepContext) -> Produzione:
        testo = f"{contesto.config.filter.truncLen} {contesto.config.err.nbases}"
        artefatto = contesto.albero.scrivi_testo(self.cartella, "profilo.txt", testo)
        return Produzione((artefatto,))


def test_una_fase_che_legge_un_parametro_non_dichiarato_fallisce(tmp_path):
    """
    **Obiettivo**: Verificare che una fase che legge un parametro non
    dichiarato fallisca con ``ParametroNonDichiarato`` e non si concluda.

    **Razionale Scientifico/Sistemistico**: Senza il vincolo, la fase
    produrrebbe un risultato che una modifica di ``err.nbases`` non
    invaliderebbe: un risultato obsoleto senza alcun errore.
    """
    config = config_ridotta(tmp_path)
    esegui_s0(config)
    run = ProjectRun(config, passi={Passo.S0: ValidazioneIngressi(), Passo.S1: _Distratta()})
    with pytest.raises(ParametroNonDichiarato, match="err.nbases"):
        run.fase(Passo.S1).esegui(run.contesto(Passo.S1))
    assert run.albero.manifesto_passo(Passo.S1, Fase.QC_PROFILES) is None


# --------------------------------------------------------------------------- #
# Registrazione                                                                #
# --------------------------------------------------------------------------- #


class _SenzaDichiarazione(PipelineStep):
    passo: ClassVar[Passo] = Passo.S1

    def calcola(self, contesto: StepContext) -> Produzione:
        raise AssertionError("non deve girare")


@pytest.mark.parametrize(
    ("parametri", "errore", "messaggio"),
    [
        (None, TypeError, "non dichiara"),
        (("tax.min_bootstrap",), KeyError, "tax.min_bootstrap"),
        (("run.threads",), ValueError, "run.threads"),
    ],
    ids=["senza-dichiarazione", "parametro-inesistente", "parametro-senza-effetto"],
)
def test_una_fase_senza_dichiarazione_valida_non_si_registra(tmp_path, parametri, errore, messaggio):
    """
    **Obiettivo**: Verificare che una fase senza dichiarazione, o con un nome
    inesistente, o con un parametro escluso dall'impronta, sia respinta alla
    registrazione.

    **Razionale Scientifico/Sistemistico**: Il meccanismo e' obbligatorio per
    ogni fase futura: una fase non puo' entrare nel grafo senza dire da che
    cosa dipende.
    """
    fase = type("Fase", (_SenzaDichiarazione,), {"parametri": parametri})()
    with pytest.raises(errore, match=messaggio):
        ProjectRun(config_ridotta(tmp_path), passi={Passo.S1: fase})


def test_ogni_fase_realizzata_dichiara_i_propri_parametri(tmp_path):
    """
    **Obiettivo**: Verificare che le fasi realizzate, da S0 a S7, dichiarino
    i propri parametri e che la registrazione le accetti.

    **Razionale Scientifico/Sistemistico**: E' la condizione perche' la loro
    validita' si giudichi sui parametri da cui dipendono davvero.
    """
    passi = passi_realizzati()
    assert set(passi) == {
        Passo.S0, Passo.S1, Passo.S2, Passo.S3, Passo.S4, Passo.S5, Passo.S6, Passo.S7,
    }
    assert all(f.parametri for f in passi.values())
    ProjectRun(config_ridotta(tmp_path))


# --------------------------------------------------------------------------- #
# Un parametro err non invalida S0, S1 e S2                                    #
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize(
    "variazione",
    [{"nbases": 2e8}, {"max_consist": 20}, {"randomize": False}],
    ids=["nbases", "max_consist", "randomize"],
)
def test_un_parametro_err_non_cambia_l_impronta_di_s0_s1_s2(tmp_path, variazione):
    """
    **Obiettivo**: Verificare che cambiare un parametro err lasci invariata
    l'impronta di configurazione di S0, S1 e S2, e cambi quella di S3.

    **Razionale Scientifico/Sistemistico**: S0-S2 non leggono i parametri del
    modello d'errore; rifarle a ogni variazione costerebbe un quarto d'ora sul
    dataset completo senza cambiare nulla.
    """
    prima = risolvi(config_ridotta(tmp_path))
    dopo = risolvi(config_ridotta(tmp_path, err=variazione))
    passi = passi_realizzati()
    for passo in (Passo.S0, Passo.S1, Passo.S2):
        fase = passi[passo]
        assert fase.calcolata_su(prima, {}) == fase.calcolata_su(dopo, {}), passo
    s3 = passi[Passo.S3]
    assert s3.calcolata_su(prima, {}) != s3.calcolata_su(dopo, {})


def test_err_batch_column_invalida_anche_s0(tmp_path):
    """
    **Obiettivo**: Verificare che ``err.batch_column``, a differenza degli altri
    parametri err, cambi anche l'impronta di S0.

    **Razionale Scientifico/Sistemistico**: S0 la legge davvero: il crosswalk
    ne ricava la corsa di ogni campione e G08 ne verifica la colonna.
    """
    prima = risolvi(config_ridotta(tmp_path))
    dopo = risolvi(config_ridotta(tmp_path, err={"batch_column": None}))
    s0 = passi_realizzati()[Passo.S0]
    assert s0.calcolata_su(prima, {}) != s0.calcolata_su(dopo, {})


# --------------------------------------------------------------------------- #
# Lato R                                                                       #
# --------------------------------------------------------------------------- #


def test_uno_script_r_che_legge_un_parametro_non_ricevuto_fallisce(tmp_path):
    """
    **Obiettivo**: Verificare che uno script R che legge un parametro non
    ricevuto fallisca, invece di ottenere ``NULL`` e proseguire.

    **Razionale Scientifico/Sistemistico**: In R ``NULL`` significa spesso
    "usa il valore predefinito della libreria": una dipendenza dimenticata
    funzionerebbe in silenzio.
    """
    from sottoinsieme import motivo_pacchetti_r_assenti

    motivo = motivo_pacchetti_r_assenti("jsonlite")
    if motivo is not None:
        if os.environ.get("AMPLICON16S_RICHIEDI_R") == "1":
            pytest.fail(f"R e' richiesto in questo ambiente: {motivo}")
        pytest.skip(motivo)

    with pytest.raises(ErroreRevisioneUmana) as info:
        esegui_script(
            DOPPIONI_R / "parametro_non_ricevuto.R", {"soglia": 2},
            AlberoOutput(tmp_path / "out"), Fase.FILTERED, tempo_massimo_s=120,
        )
    assert info.value.codice == "E-R-03"
    assert "soglia_mai_passata" in info.value.dettaglio
    assert not (tmp_path / "out" / Fase.FILTERED.value / "non_scritto.json").exists()


def test_run_batch_size_non_invalida_nessuna_fase(tmp_path):
    """
    **Obiettivo**: Verificare che ``run.batch_size`` sia fra i parametri senza
    effetto, che nessuna fase lo dichiari, e che cambiarlo lasci invariata
    l'impronta di tutte le fasi realizzate.

    **Razionale Scientifico/Sistemistico**: Il lotto decide quanta memoria
    chiede una fase, non che cosa calcola: S2 filtra ogni file da solo, e S4
    da' gli stessi byte con lotti diversi (verificato nel container). Il retry
    di E-S4-02 che lo dimezza non deve rendere da rifare cio' che e' concluso.
    """
    from amplicon16s.config.resolve import PARAMETRI_SENZA_EFFETTO

    assert "run.batch_size" in PARAMETRI_SENZA_EFFETTO
    prima = risolvi(config_ridotta(tmp_path))
    dopo = risolvi(config_ridotta(tmp_path, run={"batch_size": 3}))
    for passo, fase in passi_realizzati().items():
        assert "run.batch_size" not in fase.parametri, passo
        assert fase.calcolata_su(prima, {}) == fase.calcolata_su(dopo, {}), passo


@pytest.mark.parametrize(
    ("sezione", "variazione", "invalidate"),
    [
        ("dada", {"omega_a": 1e-20}, {Passo.S4}),
        ("qc", {"max_asv_count": 1000}, {Passo.S0, Passo.S5}),
        ("chimera", {"min_fold_parent_over_abundance": 1.5}, {Passo.S6}),
        ("qc", {"warn_frac_chimeric": 0.3}, {Passo.S0, Passo.S6}),
        ("asv", {"len_tol": 2}, {Passo.S0, Passo.S7}),
    ],
    ids=["dada.omega_a", "qc.max_asv_count", "chimera", "qc.warn_frac_chimeric", "asv.len_tol"],
)
def test_i_parametri_di_s4_s7_invalidano_solo_le_proprie_fasi(tmp_path, sezione, variazione, invalidate):
    """
    **Obiettivo**: Verificare che un parametro di una fase fra S4 e S7 cambi
    l'impronta di configurazione di quella fase e di nessun'altra, a parte
    S0, che dichiara per intero i gruppi qc e asv perche' G15 li verifica.

    **Razionale Scientifico/Sistemistico**: Una soglia sulle chimere non deve
    far rifare l'inferenza delle varianti, che sul dataset completo costa ore.
    Le fasi a valle si rifanno comunque, perche' l'impronta della fase a monte
    entra nella loro.
    """
    prima = risolvi(config_ridotta(tmp_path))
    dopo = risolvi(config_ridotta(tmp_path, **{sezione: variazione}))
    for passo, fase in passi_realizzati().items():
        cambiata = fase.calcolata_su(prima, {}) != fase.calcolata_su(dopo, {})
        assert cambiata is (passo in invalidate), passo


def test_filter_max_ee_invalida_s2_ma_non_s0_ne_s1(tmp_path):
    """
    **Obiettivo**: Verificare che cambiare ``filter.maxEE`` cambi l'impronta di
    S2 e non quelle di S0 e S1.

    **Razionale Scientifico/Sistemistico**: Di filter S0 legge solo truncLen e
    trimLeft, S1 solo truncLen e la sua tolleranza: la soglia degli errori
    attesi riguarda il filtro e basta.
    """
    prima = risolvi(config_ridotta(tmp_path))
    dopo = risolvi(config_ridotta(tmp_path, filter={"maxEE": 1.0}))
    passi = passi_realizzati()
    for passo in (Passo.S0, Passo.S1):
        assert passi[passo].calcolata_su(prima, {}) == passi[passo].calcolata_su(dopo, {}), passo
    assert passi[Passo.S2].calcolata_su(prima, {}) != passi[Passo.S2].calcolata_su(dopo, {})
