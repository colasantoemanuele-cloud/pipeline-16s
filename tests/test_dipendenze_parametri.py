"""Test delle dipendenze di parametro dichiarate da ogni fase.

Ogni fase dichiara i parametri da cui dipende, e la sua validita' si giudica
su quelli. La dichiarazione e' vincolante: il codice Python di una fase vede
la configurazione attraverso una vista ristretta, gli script R ricevono solo
cio' che la fase passa e non possono leggere altro, e una fase senza
dichiarazione non si registra. Questi test mostrano che una dipendenza
dimenticata e' un errore immediato, non un risultato obsoleto.
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
    **Obiettivo**: Verificare che S0, S1, S2 e S3 dichiarino i propri
    parametri e che la registrazione le accetti.

    **Razionale Scientifico/Sistemistico**: E' la condizione perche' la loro
    validita' si giudichi sui parametri da cui dipendono davvero.
    """
    passi = passi_realizzati()
    assert set(passi) == {Passo.S0, Passo.S1, Passo.S2, Passo.S3}
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
