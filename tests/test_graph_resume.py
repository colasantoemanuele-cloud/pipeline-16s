"""Test del grafo delle fasi, del loro completamento e della ripresa.

S0 è quella vera, eseguita su uno scenario sintetico scritto su disco; le
fasi da S1 a S14 non esistono ancora e sono sostituite da doppioni. Un
doppione scrive un artefatto che dipende dai propri ingressi, come farebbe
una fase vera, e annota ogni propria esecuzione: così i test misurano quale
lavoro la ripresa rifà davvero, non quale dichiara di rifare.

La ripresa stessa — eseguire in ordine le fasi da eseguire — è qui una
funzione dei test: l'esecutore della pipeline non è ancora realizzato.
"""

from __future__ import annotations

import json
import os
import shutil
from typing import Any, ClassVar

import pytest
from conftest import NEGATIVO, POSITIVO, Campione, crea_scenario

from amplicon16s.config.resolve import risolvi
from amplicon16s.config.schema import Config, valida
from amplicon16s.errors.exceptions import ErroreRevisioneUmana
from amplicon16s.gates.g01_g15 import Contesto
from amplicon16s.io_layer.artifacts import (
    AlberoOutput,
    Fase,
    nome_manifesto_passo,
)
from amplicon16s.runner.graph import GRAFO, Grafo, Nodo, Passo
from amplicon16s.runner.project import (
    ProjectRun,
    StatoPasso,
    passi_realizzati,
)
from amplicon16s.steps.base import (
    Esito,
    PipelineStep,
    Produzione,
    StepContext,
    impronta_parametri,
)
from amplicon16s.steps.s00_validate import ValidazioneIngressi, esegui_s0

TUTTE = tuple(Passo)


# --------------------------------------------------------------------------- #
# Scenario e doppioni                                                          #
# --------------------------------------------------------------------------- #


def _campioni() -> list[Campione]:
    return [
        Campione("ERX3000001", "NOD1D4.L1", piastra="1"),
        Campione("ERX3000002", "NOD1D4.L2", piastra="1"),
        Campione("ERX3000003", "POS.P1.1", materiale=POSITIVO,
                 posizione="Not Applicable", piastra="1"),
        Campione("ERX3000004", "BLANK.P1.1", materiale=NEGATIVO,
                 posizione="Not Applicable", piastra="1"),
    ]


def _variante(config: Config, **valori: Any) -> Config:
    """La stessa configurazione con alcuni parametri cambiati: ``gruppo__chiave``."""
    dati = config.model_dump(mode="python")
    for nome, valore in valori.items():
        gruppo, chiave = nome.split("__")
        dati[gruppo][chiave] = valore
    return valida(dati)


class Doppione(PipelineStep):
    """Una fase finta che scrive un artefatto derivato dai propri ingressi."""

    def __init__(self, registro: list[Passo]) -> None:
        self.registro = registro

    @property
    def nome_artefatto(self) -> str:
        return f"{str(self.passo).lower()}.json"

    def calcola(self, contesto: StepContext) -> Produzione:
        self.registro.append(self.passo)
        contenuto = {
            "passo": str(self.passo),
            "a_monte": {str(p): i for p, i in contesto.a_monte.items()},
            "configurazione": impronta_parametri(contesto.risolta, self.parametri),
        }
        artefatto = contesto.albero.scrivi_testo(
            self.cartella, self.nome_artefatto, json.dumps(contenuto, sort_keys=True)
        )
        return Produzione((artefatto,), {"righe": 1})


class S0Contata(ValidazioneIngressi):
    """La S0 vera, che annota le proprie esecuzioni."""

    def __init__(self, registro: list[Passo]) -> None:
        super().__init__()
        self.registro = registro

    def calcola(self, contesto: StepContext) -> Produzione:
        self.registro.append(self.passo)
        return super().calcola(contesto)


def _passi(
    registro: list[Passo], parametri: dict[Passo, tuple[str, ...] | None] | None = None
) -> dict[Passo, PipelineStep]:
    """S0 vera piu' un doppione per ogni altra fase, con i parametri indicati."""
    parametri = parametri or {}
    classe_s0 = type(
        "S0Parametrica", (S0Contata,), {"parametri": parametri.get(Passo.S0)}
    )
    passi: dict[Passo, PipelineStep] = {Passo.S0: classe_s0(registro)}
    for passo in TUTTE[1:]:
        classe = type(
            f"Doppione{passo}",
            (Doppione,),
            {"passo": passo, "parametri": parametri.get(passo)},
        )
        passi[passo] = classe(registro)
    return passi


#: Parametri ristretti: ogni fase dipende solo dal proprio gruppo, cosi' si
#: vede quali fasi invalida un singolo cambiamento.
RISTRETTI: dict[Passo, tuple[str, ...]] = {
    Passo.S0: ("io", "meta", "ctrl", "qc"),
    Passo.S1: ("qc",),
    Passo.S2: ("filter",),
    Passo.S3: ("err",),
    Passo.S4: ("asv",),
    Passo.S5: ("asv",),
    Passo.S6: ("chimera",),
    Passo.S7: ("asv",),
    Passo.S8: ("tax",),
    Passo.S9: ("phylo",),
    Passo.S10: ("out",),
    Passo.S11: ("katharoseq",),
    Passo.S12: ("decontam",),
    Passo.S13: ("prev",),
    Passo.S14: ("out",),
}


def riprendi(run: ProjectRun, registro: list[Passo]) -> list[Passo]:
    """Esegue in ordine le fasi da eseguire e restituisce quelle eseguite."""
    registro.clear()
    for _ in range(2 * len(TUTTE)):
        valutazione = run.valuta()
        passo = valutazione.prossima()
        if passo is None:
            return list(registro)
        run.fase(passo).esegui(run.contesto(passo, valutazione))
    raise AssertionError("la ripresa non converge")


@pytest.fixture
def scenario(tmp_path):
    return crea_scenario(tmp_path, _campioni(), con_letture=True)


@pytest.fixture
def registro() -> list[Passo]:
    return []


@pytest.fixture
def eseguita(scenario, registro):
    """Un'esecuzione completa, con la filogenesi disattivata come per difetto."""
    run = ProjectRun(scenario.config, passi=_passi(registro))
    riprendi(run, registro)
    return run


def _attive_da(passo: Passo, config: Config) -> list[Passo]:
    """La fase data e tutte quelle attive che ne dipendono, in ordine."""
    coinvolte = {passo, *GRAFO.discendenti(passo)}
    return [p for p in GRAFO.attive(config) if p in coinvolte]


def _file(run: ProjectRun, passo: Passo):
    fase = run.fase(passo)
    return run.albero.cartella(fase.cartella) / fase.nome_artefatto


# --------------------------------------------------------------------------- #
# Il grafo                                                                     #
# --------------------------------------------------------------------------- #


def test_il_grafo_ordina_le_quindici_fasi():
    assert GRAFO.ordine() == TUTTE
    assert len(GRAFO) == 15


def test_la_decontaminazione_precede_il_filtro_di_prevalenza():
    ordine = GRAFO.ordine()
    assert ordine.index(Passo.S12) < ordine.index(Passo.S13)
    assert Passo.S12 in GRAFO.nodo(Passo.S13).dipendenze
    assert (Passo.S12, Passo.S13, "E-S13-01") in GRAFO.precedenze


def test_ogni_dipendenza_precede_la_fase_che_la_usa():
    ordine = GRAFO.ordine()
    for nodo in GRAFO:
        for dipendenza in nodo.dipendenze:
            assert ordine.index(dipendenza) < ordine.index(nodo.passo)


def test_la_sola_fase_facoltativa_e_la_filogenesi():
    assert [n.passo for n in GRAFO if n.facoltativa] == [Passo.S9]
    assert GRAFO.nodo(Passo.S9).parametro_attivazione == "phylo.enabled"


def test_fasi_che_condividono_una_cartella():
    """Le cartelle sono quattordici, le fasi quindici."""
    per_cartella: dict[Fase, list[Passo]] = {}
    for nodo in GRAFO:
        per_cartella.setdefault(nodo.cartella, []).append(nodo.passo)
    condivise = {f: p for f, p in per_cartella.items() if len(p) > 1}
    assert condivise == {
        Fase.CHIMERA: [Passo.S6, Passo.S7],
        Fase.CONTROLS: [Passo.S11, Passo.S12],
        Fase.FINAL: [Passo.S13, Passo.S14],
    }


def _nodi_base() -> list[Nodo]:
    return [
        Nodo(Passo.S12, "a", Fase.CONTROLS),
        Nodo(Passo.S13, "b", Fase.FINAL, (Passo.S12,)),
    ]


def test_il_grafo_rifiuta_una_precedenza_obbligatoria_non_dichiarata():
    with pytest.raises(ValueError, match="S12 deve precedere S13"):
        Grafo([Nodo(Passo.S12, "a", Fase.CONTROLS), Nodo(Passo.S13, "b", Fase.FINAL)])


def test_il_grafo_rifiuta_una_precedenza_obbligatoria_facoltativa():
    nodi = _nodi_base()
    nodi[0] = Nodo(Passo.S12, "a", Fase.CONTROLS, attiva_se=lambda c: True)
    with pytest.raises(ValueError, match="facoltativa"):
        Grafo(nodi)


def test_il_grafo_rifiuta_una_dipendenza_che_non_precede():
    with pytest.raises(ValueError, match="non la precede"):
        Grafo(
            [Nodo(Passo.S13, "b", Fase.FINAL, (Passo.S12,)),
             Nodo(Passo.S12, "a", Fase.CONTROLS)],
            precedenze=(),
        )


def test_con_la_filogenesi_disattivata_s10_non_ne_dipende(scenario):
    config = scenario.config
    assert not config.phylo.enabled
    assert Passo.S9 not in GRAFO.attive(config)
    assert GRAFO.dipendenze_attive(Passo.S10, config) == (Passo.S0, Passo.S7, Passo.S8)

    attiva = _variante(config, phylo__enabled=True)
    assert Passo.S9 in GRAFO.dipendenze_attive(Passo.S10, attiva)


# --------------------------------------------------------------------------- #
# Prima esecuzione e ripresa di un'esecuzione completa                         #
# --------------------------------------------------------------------------- #


def test_la_prima_esecuzione_esegue_ogni_fase_attiva_una_volta(scenario, registro):
    run = ProjectRun(scenario.config, passi=_passi(registro))
    assert run.prossima() is Passo.S0

    eseguite = riprendi(run, registro)

    assert eseguite == [p for p in TUTTE if p is not Passo.S9]
    assert run.completa
    assert run.prossima() is None


def test_un_esecuzione_completa_non_ripete_alcun_lavoro(eseguita, scenario, registro):
    assert riprendi(eseguita, registro) == []

    # Anche un oggetto nuovo, che sa solo cio' che trova su disco.
    nuova = ProjectRun(scenario.config, passi=_passi(registro))
    assert nuova.completa
    assert riprendi(nuova, registro) == []


def test_la_filogenesi_disattivata_non_e_lavoro_mancante(eseguita):
    situazione = eseguita.situazione()
    assert situazione[Passo.S9].stato is StatoPasso.DISATTIVATA
    assert "phylo.enabled" in situazione[Passo.S9].motivo
    assert eseguita.disattivate == (Passo.S9,)
    assert Passo.S9 not in eseguita.da_eseguire
    assert eseguita.completa
    assert set(eseguita.completate) == set(TUTTE) - {Passo.S9}


def test_senza_filogenesi_s10_e_calcolata_senza_s9(eseguita):
    manifesto = eseguita.albero.manifesto_passo(Passo.S10, Fase.PHYLOSEQ)
    assert set(manifesto.calcolata_su["a_monte"]) == {"S0", "S7", "S8"}


# --------------------------------------------------------------------------- #
# Artefatti cancellati e alterati                                              #
# --------------------------------------------------------------------------- #


def test_cancellando_un_artefatto_si_riparte_da_quella_fase(eseguita, scenario, registro):
    _file(eseguita, Passo.S4).unlink()

    situazione = eseguita.situazione()
    assert situazione[Passo.S4].stato is StatoPasso.DA_ESEGUIRE
    assert "s4.json" in situazione[Passo.S4].motivo
    assert situazione[Passo.S5].motivo == "a monte da eseguire: S4"
    assert eseguita.prossima() is Passo.S4

    eseguite = riprendi(eseguita, registro)
    assert eseguite == _attive_da(Passo.S4, scenario.config)
    # A monte e di lato non si rifa' nulla.
    assert not {Passo.S0, Passo.S1, Passo.S2, Passo.S3} & set(eseguite)
    assert eseguita.completa


def test_cancellando_un_artefatto_di_s0_si_rifa_tutto(eseguita, scenario, registro):
    cartella = eseguita.albero.cartella(Fase.INPUT_VALIDATION)
    (cartella / "crosswalk.tsv").unlink()
    assert riprendi(eseguita, registro) == list(GRAFO.attive(scenario.config))


def test_ricalcolare_una_fase_senza_dipendenti_non_tocca_le_altre(eseguita, registro):
    """S1 produce profili che nessuna fase legge: rifarla non invalida il resto."""
    _file(eseguita, Passo.S1).unlink()
    assert riprendi(eseguita, registro) == [Passo.S1]


def test_un_artefatto_alterato_fa_rieseguire_la_fase(eseguita, scenario, registro):
    percorso = _file(eseguita, Passo.S8)
    percorso.write_text(percorso.read_text() + " ", encoding="utf-8")

    situazione = eseguita.situazione()
    assert situazione[Passo.S8].stato is StatoPasso.DA_ESEGUIRE
    assert "alterati: s8.json" in situazione[Passo.S8].motivo

    assert riprendi(eseguita, registro) == _attive_da(Passo.S8, scenario.config)


def test_un_manifesto_alterato_non_vale_come_conclusione(eseguita, registro):
    percorso = eseguita.albero.percorso_manifesto_passo(Passo.S5, Fase.SEQTAB)
    documento = json.loads(percorso.read_text())
    documento["metriche"]["righe"] = 999
    percorso.write_text(json.dumps(documento))

    situazione = eseguita.situazione()
    assert situazione[Passo.S5].stato is StatoPasso.DA_ESEGUIRE
    assert "alterato" in situazione[Passo.S5].motivo


def test_una_fase_interrotta_non_risulta_conclusa(eseguita, scenario, registro):
    """Il manifesto precedente sparisce all'avvio: un calcolo a meta' non conta."""

    class Interrotta(Doppione):
        passo: ClassVar[Passo] = Passo.S6

        def calcola(self, contesto):
            contesto.albero.scrivi_testo(self.cartella, "parziale.tsv", "meta'")
            raise KeyboardInterrupt

    with pytest.raises(KeyboardInterrupt):
        Interrotta(registro).esegui(eseguita.contesto(Passo.S6))

    assert eseguita.albero.manifesto_passo(Passo.S6, Fase.CHIMERA) is None
    assert eseguita.prossima() is Passo.S6
    assert riprendi(eseguita, registro) == _attive_da(Passo.S6, scenario.config)


# --------------------------------------------------------------------------- #
# Su che cosa e' stata calcolata una fase                                      #
# --------------------------------------------------------------------------- #


def test_cambiando_un_parametro_nessuna_fase_resta_conclusa(eseguita, scenario, registro):
    """Per difetto una fase dipende dall'intera configurazione."""
    cambiata = _variante(scenario.config, decontam__threshold=0.4)
    run = ProjectRun(cambiata, passi=_passi(registro))

    situazione = run.situazione()
    assert run.completate == ()
    assert situazione[Passo.S0].motivo == "configurazione cambiata"
    assert situazione[Passo.S9].stato is StatoPasso.DISATTIVATA

    assert riprendi(run, registro) == list(GRAFO.attive(cambiata))

    # Tornare al valore precedente e' a sua volta un cambiamento.
    tornata = ProjectRun(scenario.config, passi=_passi(registro))
    assert tornata.completate == ()


def test_con_parametri_ristretti_si_rifa_solo_cio_che_ne_dipende(scenario, registro):
    config = _variante(scenario.config, phylo__enabled=True)
    run = ProjectRun(config, passi=_passi(registro, RISTRETTI))
    assert riprendi(run, registro) == list(TUTTE)

    cambiata = _variante(config, tax__min_boot=80)
    run = ProjectRun(cambiata, passi=_passi(registro, RISTRETTI))
    situazione = run.situazione()
    assert situazione[Passo.S8].motivo == "configurazione cambiata"
    assert situazione[Passo.S10].motivo == "a monte da eseguire: S8"

    eseguite = riprendi(run, registro)
    assert eseguite == [Passo.S8, Passo.S10, Passo.S11, Passo.S12, Passo.S13, Passo.S14]
    assert Passo.S9 not in eseguite  # la filogenesi non legge la tassonomia


def test_attivare_la_filogenesi_rifa_s9_e_cio_che_ne_dipende(scenario, registro):
    run = ProjectRun(scenario.config, passi=_passi(registro, RISTRETTI))
    riprendi(run, registro)

    attiva = _variante(scenario.config, phylo__enabled=True)
    run = ProjectRun(attiva, passi=_passi(registro, RISTRETTI))
    situazione = run.situazione()
    assert situazione[Passo.S9].motivo == "nessun manifesto: mai conclusa"
    assert situazione[Passo.S8].stato is StatoPasso.COMPLETATA

    assert riprendi(run, registro) == [
        Passo.S9, Passo.S10, Passo.S11, Passo.S12, Passo.S13, Passo.S14
    ]

    # E disattivarla di nuovo rifa' S10, che non deve piu' contenere l'albero.
    run = ProjectRun(scenario.config, passi=_passi(registro, RISTRETTI))
    assert "dipendenze cambiate" in run.situazione()[Passo.S10].motivo
    assert riprendi(run, registro) == [
        Passo.S10, Passo.S11, Passo.S12, Passo.S13, Passo.S14
    ]


def test_dati_di_ingresso_sostituiti_invalidano_s0_e_il_resto(eseguita, scenario, registro):
    """Un file di letture sostituito sotto lo stesso nome, a configurazione invariata."""
    fastq = sorted(scenario.config.io.fastq_dir.glob(scenario.config.io.fastq_glob))[0]
    stato = fastq.stat()
    os.utime(fastq, ns=(stato.st_atime_ns, stato.st_mtime_ns + 1_000_000_000))

    situazione = eseguita.situazione()
    assert situazione[Passo.S0].motivo == "dati di ingresso cambiati"
    assert riprendi(eseguita, registro) == list(GRAFO.attive(scenario.config))


def test_un_parametro_dichiarato_male_e_un_errore(scenario):
    risolta = risolvi(scenario.config)
    with pytest.raises(KeyError, match="tax.min_bootstrap"):
        impronta_parametri(risolta, ("tax.min_bootstrap",))
    with pytest.raises(KeyError, match="tassonomia"):
        impronta_parametri(risolta, ("tassonomia",))
    assert impronta_parametri(risolta, None) == risolta.impronta_risultati
    assert risolta.impronta_risultati != risolta.digest
    with pytest.raises(ValueError, match="run.threads"):
        impronta_parametri(risolta, ("run.threads",))


# --------------------------------------------------------------------------- #
# Parametri che non incidono sui risultati                                     #
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize(
    "variazione",
    [
        {"run__threads": 2},
        {"retry__enabled": False},
        {"retry__max_attempts": 7},
        {"run__threads": 2, "retry__enabled": False, "retry__max_attempts": 7},
    ],
    ids=["threads", "retry.enabled", "retry.max_attempts", "tutti"],
)
def test_un_parametro_escluso_non_invalida_alcuna_fase(
    eseguita, scenario, registro, variazione
):
    cambiata = _variante(scenario.config, **variazione)
    run = ProjectRun(cambiata, passi=_passi(registro))
    assert run.completa
    assert riprendi(run, registro) == []

    # Il digest scritto in 00_config identifica comunque la configurazione.
    assert risolvi(cambiata).digest != risolvi(scenario.config).digest


def test_spostare_la_cartella_di_output_non_rende_incompleta_l_esecuzione(
    eseguita, scenario, registro, tmp_path
):
    nuova_radice = tmp_path / "spostata"
    shutil.copytree(scenario.config.io.out_root, nuova_radice)
    spostata = _variante(scenario.config, io__out_root=nuova_radice)

    run = ProjectRun(spostata, passi=_passi(registro))
    assert run.completa
    assert riprendi(run, registro) == []


@pytest.mark.parametrize(
    "variazione",
    [{"run__batch_size": 7}, {"run__lockfile": "altro.lock"}, {"decontam__threshold": 0.4}],
    ids=["batch_size", "lockfile", "decontam.threshold"],
)
def test_un_parametro_incluso_invalida_le_fasi(eseguita, scenario, registro, variazione):
    cambiata = _variante(scenario.config, **variazione)
    run = ProjectRun(cambiata, passi=_passi(registro))
    assert run.completate == ()
    assert run.situazione()[Passo.S0].motivo == "configurazione cambiata"


def test_da_un_gruppo_i_parametri_esclusi_restano_fuori(scenario):
    """Una fase che dichiara il gruppo run non si invalida per run.threads."""
    prima = impronta_parametri(risolvi(scenario.config), ("run",))
    dopo = impronta_parametri(risolvi(_variante(scenario.config, run__threads=2)), ("run",))
    assert prima == dopo
    lotti = impronta_parametri(
        risolvi(_variante(scenario.config, run__batch_size=7)), ("run",)
    )
    assert lotti != prima


# --------------------------------------------------------------------------- #
# Costo della valutazione                                                      #
# --------------------------------------------------------------------------- #


@pytest.fixture
def conteggio_checksum(monkeypatch):
    """Conta i checksum calcolati durante la valutazione dello stato."""
    import amplicon16s.runner.project as progetto

    calcolati: list = []
    originale = progetto.checksum_file

    def contato(percorso):
        calcolati.append(percorso)
        return originale(percorso)

    monkeypatch.setattr(progetto, "checksum_file", contato)
    return calcolati


def test_una_valutazione_calcola_ogni_checksum_una_volta(eseguita, conteggio_checksum):
    valutazione = eseguita.valuta()
    artefatti = {
        eseguita.albero.cartella(eseguita.fase(p).cartella) / nome
        for p in valutazione.completate
        for nome in eseguita.albero.manifesto_passo(p, eseguita.fase(p).cartella).nomi
    }
    assert len(artefatti) == 17  # 4 di S0 e uno per ciascuna delle 13 fasi doppione
    assert len(conteggio_checksum) == len(artefatti)
    assert set(conteggio_checksum) == artefatti

    # Tutte le domande poste alla stessa valutazione non ne calcolano altri.
    for _ in range(3):
        valutazione.prossima()
        valutazione.completate
        valutazione.disattivate
        valutazione.da_eseguire
        valutazione.completa
        eseguita.contesto(Passo.S14, valutazione)
    assert len(conteggio_checksum) == len(artefatti)


def test_le_scorciatoie_compiono_una_valutazione_nuova(eseguita, conteggio_checksum):
    """Una domanda isolata rilegge il disco: e' cio' che la rende attuale."""
    eseguita.prossima()
    primo = len(conteggio_checksum)
    _file(eseguita, Passo.S3).unlink()
    assert eseguita.prossima() is Passo.S3
    assert len(conteggio_checksum) > primo


# --------------------------------------------------------------------------- #
# Cartelle condivise                                                           #
# --------------------------------------------------------------------------- #


def test_la_cartella_completa_non_rende_completa_la_fase_accanto(scenario, registro):
    """Il difetto del manifesto per cartella, e la sua soluzione."""
    run = ProjectRun(scenario.config, passi=_passi(registro))
    for passo in GRAFO.attive(scenario.config):
        if passo is Passo.S12:
            break
        run.fase(passo).esegui(run.contesto(passo))

    # Il manifesto della cartella direbbe che 11_controls e' completa...
    assert run.albero.fase_completa(Fase.CONTROLS)
    # ...ma S12 non ha mai girato, e il suo criterio e' il suo manifesto.
    situazione = run.situazione()
    assert situazione[Passo.S11].stato is StatoPasso.COMPLETATA
    assert situazione[Passo.S12].stato is StatoPasso.DA_ESEGUIRE
    assert situazione[Passo.S12].motivo == "nessun manifesto: mai conclusa"


@pytest.mark.parametrize(
    ("prima", "dopo", "cartella"),
    [
        (Passo.S6, Passo.S7, Fase.CHIMERA),
        (Passo.S11, Passo.S12, Fase.CONTROLS),
        (Passo.S13, Passo.S14, Fase.FINAL),
    ],
)
def test_due_fasi_nella_stessa_cartella_si_concludono_separatamente(
    eseguita, registro, prima, dopo, cartella
):
    albero = eseguita.albero
    assert albero.manifesto_passo(prima, cartella).nomi == (f"{prima.lower()}.json",)
    assert albero.manifesto_passo(dopo, cartella).nomi == (f"{dopo.lower()}.json",)

    # Rifare la seconda non tocca la conclusione della prima.
    impronta_prima = albero.manifesto_passo(prima, cartella).impronta
    _file(eseguita, dopo).unlink()
    situazione = eseguita.situazione()
    assert situazione[prima].stato is StatoPasso.COMPLETATA
    assert situazione[dopo].stato is StatoPasso.DA_ESEGUIRE

    eseguite = riprendi(eseguita, registro)
    assert prima not in eseguite and dopo in eseguite
    assert albero.manifesto_passo(prima, cartella).impronta == impronta_prima


def test_una_fase_che_riscrive_il_file_dell_altra_la_invalida(eseguita):
    """Nella stessa cartella: il checksum registrato dalla prima non torna piu'."""
    albero = eseguita.albero
    albero.scrivi_testo(Fase.CONTROLS, "s11.json", "sovrascritto da S12")
    assert eseguita.situazione()[Passo.S11].stato is StatoPasso.DA_ESEGUIRE


def test_i_manifesti_di_fase_hanno_nomi_distinti():
    assert nome_manifesto_passo("S11") != nome_manifesto_passo("S12")
    assert nome_manifesto_passo("S11") != "manifest.json"


# --------------------------------------------------------------------------- #
# Prerequisiti                                                                 #
# --------------------------------------------------------------------------- #


def test_s13_prima_di_s12_e_rifiutata(eseguita):
    eseguita.albero.rimuovi_manifesto_passo(Passo.S12, Fase.CONTROLS)

    with pytest.raises(ErroreRevisioneUmana) as info:
        eseguita.fase(Passo.S13).esegui(eseguita.contesto(Passo.S13))
    assert info.value.codice == "E-S13-01"
    assert info.value.contesto["mancanti"] == ["S12"]
    # Il rifiuto avviene prima di toccare alcunche'.
    assert eseguita.albero.manifesto_passo(Passo.S13, Fase.FINAL) is not None


def test_nessuna_fase_gira_prima_delle_sue_dipendenze(scenario, registro):
    """La violazione generica ha il proprio codice, non quello di S13."""
    run = ProjectRun(scenario.config, passi=_passi(registro))
    with pytest.raises(ErroreRevisioneUmana, match="S0") as info:
        run.fase(Passo.S2).esegui(run.contesto(Passo.S2))
    assert info.value.codice == "E-GRAFO-01"
    assert info.value.fase == "GRAFO"
    assert registro == []


def test_s13_senza_s11_ne_s12_porta_il_codice_della_precedenza(eseguita):
    """Se manca anche S12, conta la precedenza di metodo."""
    for passo in (Passo.S11, Passo.S12):
        eseguita.albero.rimuovi_manifesto_passo(passo, Fase.CONTROLS)
    contesto = eseguita.contesto(Passo.S13)
    with pytest.raises(ErroreRevisioneUmana) as info:
        eseguita.fase(Passo.S13).esegui(contesto)
    assert info.value.codice == "E-S13-01"


def test_s12_senza_s11_porta_il_codice_generico(eseguita):
    eseguita.albero.rimuovi_manifesto_passo(Passo.S11, Fase.CONTROLS)
    with pytest.raises(ErroreRevisioneUmana) as info:
        eseguita.fase(Passo.S12).esegui(eseguita.contesto(Passo.S12))
    assert info.value.codice == "E-GRAFO-01"


def test_una_precedenza_con_un_codice_inesistente_e_rifiutata():
    with pytest.raises(KeyError, match="E-S99-01"):
        Grafo(_nodi_base(), precedenze=((Passo.S12, Passo.S13, "E-S99-01"),))


# --------------------------------------------------------------------------- #
# S0 sulla classe base, e l'esecuzione reale di oggi                           #
# --------------------------------------------------------------------------- #


def test_s0_invocata_da_sola_risulta_conclusa_per_l_esecuzione(scenario):
    esegui_s0(scenario.config)
    run = ProjectRun(scenario.config)
    assert run.situazione()[Passo.S0].stato is StatoPasso.COMPLETATA


def test_s0_non_superata_non_viene_registrata_come_conclusa(tmp_path):
    scenario = crea_scenario(
        tmp_path, _campioni(), con_letture=True, sovrascrivi={"run": {"threads": 100000}}
    )
    fase = ValidazioneIngressi(solleva=False)
    run = ProjectRun(scenario.config)
    risultato = fase.esegui(run.contesto(Passo.S0))

    assert risultato.esito is Esito.NON_SUPERATA
    assert risultato.impronta is None
    assert run.albero.manifesto_passo(Passo.S0, Fase.INPUT_VALIDATION) is None
    assert run.prossima() is Passo.S0


def test_il_risultato_porta_artefatti_con_checksum_e_metriche(scenario):
    run = ProjectRun(scenario.config)
    risultato = run.fase(Passo.S0).esegui(run.contesto(Passo.S0))

    assert risultato.completata
    assert {a.nome for a in risultato.artefatti} == {
        "gates.json", "crosswalk.tsv", "inventario.json", "letture_ispezionate.tsv"
    }
    assert all(a.checksum.startswith("sha256:") for a in risultato.artefatti)
    assert risultato.metriche["campioni"] == 4
    manifesto = run.albero.manifesto_passo(Passo.S0, Fase.INPUT_VALIDATION)
    assert manifesto.impronta == risultato.impronta


def test_l_inventario_e_riletto_dall_artefatto_di_s0(scenario):
    run = ProjectRun(scenario.config)
    assert run.inventario is None
    esegui_s0(scenario.config)

    assert run.inventario == Contesto(scenario.config).inventario
    assert run.risolta.derivati.prev_min_samples is not None


def test_oggi_esiste_solo_s0(scenario):
    run = ProjectRun(scenario.config)
    assert set(passi_realizzati()) == {Passo.S0}
    esegui_s0(scenario.config)

    situazione = run.situazione()
    assert situazione[Passo.S0].stato is StatoPasso.COMPLETATA
    assert situazione[Passo.S1].stato is StatoPasso.NON_REALIZZATA
    assert run.prossima() is Passo.S1
    assert not run.completa
    with pytest.raises(LookupError, match="S1"):
        run.fase(Passo.S1)


def test_l_albero_e_quello_della_configurazione(scenario):
    run = ProjectRun(scenario.config)
    assert run.albero.radice == AlberoOutput(scenario.config.io.out_root).radice
