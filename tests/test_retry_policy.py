"""Test dell'esecutore, della politica dei tentativi e della riga di comando.

S0 è quella vera, su uno scenario sintetico; le fasi da S1 a S14 sono doppioni.
Un doppione annota ogni esecuzione con i parametri che ha visto, così i test
verificano l'aggiustamento effettivamente usato, non quello dichiarato. Un
doppione *fragile* fallisce con un codice dato per i primi tentativi e poi
riesce, come farebbe una fase vera dopo l'azione correttiva.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, ClassVar

import pytest
import yaml
from conftest import NEGATIVO, POSITIVO, Campione, crea_scenario

import amplicon16s.cli as cli
from amplicon16s.config.resolve import risolvi, versioni_registrate
from amplicon16s.config.schema import Config, valida
from amplicon16s.errors.exceptions import errore
from amplicon16s.gates.g01_g15 import ErroreGate
from amplicon16s.io_layer.artifacts import Fase
from amplicon16s.logging.logger import NOME_FILE_LOG, chiudi, configura
from amplicon16s.runner.executor import Conclusione, Esecutore
from amplicon16s.runner.graph import GRAFO, Passo
from amplicon16s.runner.project import ProjectRun, StatoPasso
from amplicon16s.runner.retry import (
    Aggiustamento,
    PoliticaRetry,
    dimezza,
    raddoppia,
)
from amplicon16s.steps.base import PipelineStep, Produzione, StepContext
from amplicon16s.steps.s00_validate import ValidazioneIngressi

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
    dati = config.model_dump(mode="python")
    for nome, v in valori.items():
        gruppo, chiave = nome.split("__")
        dati[gruppo][chiave] = v
    return valida(dati)


class Doppione(PipelineStep):
    """Scrive un artefatto e annota i parametri che ha visto."""

    #: Per i fragili: codice sollevato e quanti tentativi falliscono.
    codice: ClassVar[str | None] = None
    fallimenti: ClassVar[int] = 0
    #: "degrada": registra una degradazione e prosegue; "solleva": solleva
    #: una DegradazioneRichiesta; con ``ripiego`` la fase sa proseguire.
    degradazione: ClassVar[str | None] = None
    codice_degradazione: ClassVar[str] = "E-S6-02"
    ripiego: ClassVar[bool] = False

    def __init__(self, registro: list[tuple[Passo, int, float]]) -> None:
        self.registro = registro
        self.tentativi = 0

    def _scrivi(self, contesto: StepContext, testo: str = "") -> Produzione:
        artefatto = contesto.albero.scrivi_testo(
            self.cartella, f"{str(self.passo).lower()}.txt", f"{self.passo} {testo}"
        )
        return Produzione((artefatto,))

    def calcola(self, contesto: StepContext) -> Produzione:
        config = contesto.config
        self.registro.append((self.passo, config.run.batch_size, config.err.nbases))
        self.tentativi += 1
        if self.codice and self.tentativi <= self.fallimenti:
            raise errore(self.codice, f"doppione, tentativo {self.tentativi}")
        if self.degradazione == "degrada":
            contesto.degrada(self.codice_degradazione, "frazione chimerica 0,31")
        elif self.degradazione == "solleva":
            raise errore(self.codice_degradazione, "soglia non attendibile")
        return self._scrivi(contesto)

    def ripiega(self, contesto, degradazione):
        if not self.ripiego:
            return super().ripiega(contesto, degradazione)
        return self._scrivi(contesto, "ripiego")


def _passi(registro, **speciali: dict[str, Any]) -> dict[Passo, PipelineStep]:
    """S0 vera e un doppione per ogni altra fase; ``S4={...}`` ne cambia uno."""
    passi: dict[Passo, PipelineStep] = {Passo.S0: ValidazioneIngressi()}
    for passo in TUTTE[1:]:
        attributi = {"passo": passo, **speciali.get(str(passo), {})}
        passi[passo] = type(f"Doppione{passo}", (Doppione,), attributi)(registro)
    return passi


FRAGILE_S4 = {"codice": "E-S4-02", "fallimenti": 1,
              "aggiustamenti": {"E-S4-02": dimezza("run.batch_size")}}


@pytest.fixture(autouse=True)
def uscite_pulite():
    chiudi()
    yield
    chiudi()


@pytest.fixture
def scenario(tmp_path):
    return crea_scenario(tmp_path, _campioni(), con_letture=True)


@pytest.fixture
def registro() -> list:
    return []


def _esegui(config, registro, **speciali):
    run = ProjectRun(config, passi=_passi(registro, **speciali))
    return run, Esecutore(run, comando_ripresa="amplicon16s resume --config prova.yaml").esegui()


def _visti(registro, passo) -> list[int]:
    return [batch for p, batch, _ in registro if p is passo]


# --------------------------------------------------------------------------- #
# La politica                                                                  #
# --------------------------------------------------------------------------- #


def test_max_attempts_conta_i_tentativi_totali_compreso_il_primo(scenario):
    politica = PoliticaRetry.da_config(scenario.config)
    assert politica.tentativi_massimi == 2
    config = scenario.config
    adatta = dimezza("run.batch_size")
    assert politica.decidi("E-S4-02", 1, adatta, config).ritenta
    assert not politica.decidi("E-S4-02", 2, adatta, config).ritenta


def test_un_aggiustamento_su_un_parametro_non_ammesso_e_rifiutato():
    with pytest.raises(ValueError, match="filter.maxEE"):
        Aggiustamento("filter.maxEE", lambda v: v * 2, "soglia alzata")


def test_una_fase_non_puo_aggiustare_un_codice_non_ripetibile(scenario, registro):
    passi = _passi(registro, S2={"aggiustamenti": {"E-S2-01": dimezza("run.batch_size")}})
    with pytest.raises(ValueError, match="E-S2-01"):
        Esecutore(ProjectRun(scenario.config, passi=passi))


# --------------------------------------------------------------------------- #
# Ritentare con l'aggiustamento                                                #
# --------------------------------------------------------------------------- #


def test_un_codice_ammesso_viene_ritentato_con_l_aggiustamento(scenario, registro):
    run, esito = _esegui(scenario.config, registro, S4=FRAGILE_S4)

    assert esito.conclusione is Conclusione.COMPLETATA
    assert _visti(registro, Passo.S4) == [24, 12]
    # L'esecuzione prosegue dopo il tentativo riuscito, con il valore dichiarato.
    assert _visti(registro, Passo.S5) == [24]
    assert run.valuta().completa


def test_l_aggiustamento_e_nel_manifesto_della_fase(scenario, registro):
    run, _ = _esegui(scenario.config, registro, S4=FRAGILE_S4)

    manifesto = run.albero.manifesto_passo(Passo.S4, Fase.ASV_INFERENCE)
    (aggiustamento,) = manifesto.aggiustamenti
    assert aggiustamento["codice"] == "E-S4-02"
    assert aggiustamento["parametro"] == "run.batch_size"
    assert aggiustamento["dichiarato"] == 24
    assert aggiustamento["usato"] == 12
    assert aggiustamento["tentativo_fallito"] == 1
    # Le altre fasi non ne hanno.
    assert run.albero.manifesto_passo(Passo.S5, Fase.SEQTAB).aggiustamenti == ()


def test_l_aggiustamento_e_nel_log(scenario, registro):
    configura(scenario.config.io.out_root)
    _esegui(scenario.config, registro, S4=FRAGILE_S4)
    righe = [
        json.loads(r)
        for r in (Path(scenario.config.io.out_root) / Fase.LOGS.value / NOME_FILE_LOG)
        .read_text().splitlines()
    ]
    ritento = [r for r in righe if r.get("parametro") == "run.batch_size"]
    assert len(ritento) == 1
    assert ritento[0]["codice"] == "E-S4-02"
    assert ritento[0]["usato"] == 12
    assert any(r.get("codice") == "E-S4-02" and r["livello"] == "WARNING" for r in righe)


def test_la_ripresa_non_rifa_la_fase_riuscita_con_l_aggiustamento(scenario, registro):
    _esegui(scenario.config, registro, S4=FRAGILE_S4)
    registro.clear()

    run, esito = _esegui(scenario.config, registro)
    situazione = run.valuta().situazioni[Passo.S4]
    assert situazione.stato is StatoPasso.COMPLETATA
    assert esito.conclusione is Conclusione.COMPLETATA
    assert esito.eseguite == ()
    assert registro == []


def test_tre_tentativi_totali_con_due_fallimenti(scenario, registro):
    config = _variante(scenario.config, retry__max_attempts=3)
    fragile = {**FRAGILE_S4, "fallimenti": 2}
    run, esito = _esegui(config, registro, S4=fragile)

    assert esito.conclusione is Conclusione.COMPLETATA
    assert _visti(registro, Passo.S4) == [24, 12, 6]
    usati = [a["usato"] for a in run.albero.manifesto_passo(Passo.S4, Fase.ASV_INFERENCE).aggiustamenti]
    assert usati == [12, 6]


def test_il_modello_d_errore_si_ritenta_con_piu_basi(scenario, registro):
    fragile = {"codice": "E-S3-01", "fallimenti": 1,
               "aggiustamenti": {"E-S3-01": raddoppia("err.nbases")}}
    _, esito = _esegui(scenario.config, registro, S3=fragile)
    assert esito.conclusione is Conclusione.COMPLETATA
    assert [n for p, _, n in registro if p is Passo.S3] == [1e8, 2e8]


# --------------------------------------------------------------------------- #
# Fermarsi                                                                     #
# --------------------------------------------------------------------------- #


def test_fallendo_due_volte_si_ferma_con_il_punto_di_ripresa(scenario, registro):
    fragile = {**FRAGILE_S4, "fallimenti": 99}
    run, esito = _esegui(scenario.config, registro, S4=fragile)

    assert esito.conclusione is Conclusione.ARRESTATA
    punto = esito.punto
    assert punto.passo is Passo.S4
    assert punto.codice == "E-S4-02"
    assert punto.tentativi == 2 and punto.tentativi_massimi == 2
    assert "run.batch_size ridotto" in punto.azione  # il messaggio del catalogo
    assert punto.comando == "amplicon16s resume --config prova.yaml"
    assert _visti(registro, Passo.S4) == [24, 12]
    assert _visti(registro, Passo.S5) == []

    # Anche su disco, per dopo la chiusura del terminale.
    cartella = run.albero.cartella(Fase.LOGS)
    documento = json.loads((cartella / "punto_di_ripresa.json").read_text())
    assert documento["passo"] == "S4"
    assert documento["codice"] == "E-S4-02"
    assert documento["tentativi"] == 2
    assert documento["comando"] == punto.comando
    testo = (cartella / "punto_di_ripresa.txt").read_text()
    assert "ESECUZIONE FERMATA nella fase S4" in testo and punto.comando in testo


def test_un_esecuzione_completata_rimuove_il_vecchio_punto_di_ripresa(scenario, registro):
    _esegui(scenario.config, registro, S4={**FRAGILE_S4, "fallimenti": 99})
    cartella = Path(scenario.config.io.out_root) / Fase.LOGS.value
    assert (cartella / "punto_di_ripresa.json").exists()

    _, esito = _esegui(scenario.config, registro)
    assert esito.conclusione is Conclusione.COMPLETATA
    assert not (cartella / "punto_di_ripresa.json").exists()


def test_un_codice_a_revisione_umana_si_ferma_al_primo_tentativo(scenario, registro):
    fragile = {"codice": "E-S2-01", "fallimenti": 99}
    _, esito = _esegui(scenario.config, registro, S2=fragile)
    assert esito.punto.tentativi == 1
    assert esito.punto.categoria == "revisione_umana"
    assert "revisione umana" in esito.punto.motivo
    assert _visti(registro, Passo.S2) == [24]


def test_un_codice_ripetibile_fuori_whitelist_si_ferma_al_primo(scenario, registro):
    config = _variante(scenario.config, retry__whitelist=["E-S2-03"])
    _, esito = _esegui(config, registro, S4={**FRAGILE_S4, "fallimenti": 99})
    assert esito.punto.tentativi == 1
    assert "retry.whitelist" in esito.punto.motivo
    assert _visti(registro, Passo.S4) == [24]


def test_con_retry_disattivato_nessun_codice_viene_ritentato(scenario, registro):
    config = _variante(scenario.config, retry__enabled=False)
    fragile_s3 = {"codice": "E-S3-01", "fallimenti": 1,
                  "aggiustamenti": {"E-S3-01": raddoppia("err.nbases")}}
    _, esito = _esegui(config, registro, S3=fragile_s3, S4=FRAGILE_S4)
    assert esito.conclusione is Conclusione.ARRESTATA
    assert esito.punto.passo is Passo.S3 and esito.punto.tentativi == 1
    assert "retry.enabled" in esito.punto.motivo


def test_senza_azione_correttiva_non_si_ritenta_identico(scenario, registro):
    fragile = {"codice": "E-S4-02", "fallimenti": 1}
    _, esito = _esegui(scenario.config, registro, S4=fragile)
    assert esito.punto.tentativi == 1
    assert "azione correttiva" in esito.punto.motivo


def test_un_azione_correttiva_al_limite_non_si_ripete(scenario, registro):
    config = _variante(scenario.config, run__batch_size=1)
    _, esito = _esegui(config, registro, S4=FRAGILE_S4)
    assert esito.punto.tentativi == 1
    assert "limite" in esito.punto.motivo


def test_il_messaggio_distingue_retry_e_retry_poi_revisione(scenario, registro):
    s3 = {"codice": "E-S3-01", "fallimenti": 99,
          "aggiustamenti": {"E-S3-01": raddoppia("err.nbases")}}
    _, dopo_s3 = _esegui(scenario.config, registro, S3=s3)
    assert dopo_s3.punto.categoria == "retry_poi_revisione"
    assert "serve la revisione umana" in dopo_s3.punto.motivo

    _, dopo_s4 = _esegui(
        _variante(scenario.config, io__out_root=Path(scenario.radice) / "altra"),
        registro, S4={**FRAGILE_S4, "fallimenti": 99},
    )
    assert dopo_s4.punto.categoria == "retry_automatico"
    assert "serve la revisione umana" not in dopo_s4.punto.motivo
    assert "ammette il retry" in dopo_s4.punto.motivo


# --------------------------------------------------------------------------- #
# Whitelist incoerente con il catalogo                                         #
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize(
    "whitelist",
    [["E-S4-02", "E-S2-01"], ["E-S99-01"]],
    ids=["codice-non-ripetibile", "codice-inesistente"],
)
def test_una_whitelist_incoerente_e_respinta_prima_di_ogni_fase(
    scenario, registro, whitelist
):
    config = _variante(scenario.config, retry__whitelist=whitelist)
    with pytest.raises(ErroreGate) as info:
        _esegui(config, registro)
    assert info.value.gate == "G15"
    assert "E-G15-09" in {v.codice for v in info.value.violazioni}
    assert registro == []
    assert not (Path(config.io.out_root) / Fase.INPUT_VALIDATION.value).exists()


def test_una_whitelist_ristretta_e_ammessa(scenario, registro):
    config = _variante(scenario.config, retry__whitelist=["E-S4-02"])
    _, esito = _esegui(config, registro)
    assert esito.conclusione is Conclusione.COMPLETATA


# --------------------------------------------------------------------------- #
# Degradazioni                                                                 #
# --------------------------------------------------------------------------- #


def test_una_degradazione_non_interrompe_ed_e_registrata(scenario, registro):
    run, esito = _esegui(scenario.config, registro, S6={"degradazione": "degrada"})

    assert esito.conclusione is Conclusione.COMPLETATA
    (degradazione,) = run.albero.manifesto_passo(Passo.S6, Fase.CHIMERA).degradazioni
    assert degradazione["codice"] == "E-S6-02"
    assert degradazione["categoria"] == "degradazione_automatica"
    assert degradazione["dettaglio"] == "frazione chimerica 0,31"


def test_una_degradazione_sollevata_con_ripiego_non_interrompe(scenario, registro):
    s11 = {"degradazione": "solleva", "codice_degradazione": "E-S11-02", "ripiego": True}
    run, esito = _esegui(scenario.config, registro, S11=s11)

    assert esito.conclusione is Conclusione.COMPLETATA
    manifesto = run.albero.manifesto_passo(Passo.S11, Fase.CONTROLS)
    assert [d["codice"] for d in manifesto.degradazioni] == ["E-S11-02"]
    assert "ripiego" in (run.albero.cartella(Fase.CONTROLS) / "s11.txt").read_text()


def test_una_degradazione_sollevata_senza_ripiego_e_un_difetto(scenario, registro):
    with pytest.raises(RuntimeError, match="senza dichiarare un ripiego"):
        _esegui(scenario.config, registro, S11={"degradazione": "solleva",
                                                "codice_degradazione": "E-S11-02"})


def test_degradare_con_un_codice_che_ferma_e_rifiutato(scenario):
    run = ProjectRun(scenario.config)
    contesto = run.contesto(Passo.S0)
    with pytest.raises(ValueError, match="revisione_umana"):
        contesto.degrada("E-S2-01", "non e' una degradazione")


def test_le_degradazioni_di_s0_finiscono_nel_suo_manifesto(tmp_path, registro):
    """Oggi nascono come avvisi dei gate e non fermano: E-S0-15 da G08, e
    E-S1-01 da G09, perche' le letture sintetiche sono di 151 bp contro un
    troncamento di 137."""
    scenario = crea_scenario(
        tmp_path, _campioni(), con_letture=True, con_arricchimento=True,
        sovrascrivi={"decontam": {"min_blanks": 3}},
    )
    run, esito = _esegui(scenario.config, registro)
    assert esito.conclusione is Conclusione.COMPLETATA
    codici = [d["codice"] for d in run.albero.manifesto_passo(Passo.S0, Fase.INPUT_VALIDATION).degradazioni]
    assert sorted(codici) == ["E-S0-15", "E-S1-01"]


# --------------------------------------------------------------------------- #
# Controlli di risorse a ogni avvio                                            #
# --------------------------------------------------------------------------- #


def test_i_controlli_di_risorse_girano_anche_con_s0_conclusa(scenario, registro):
    _esegui(scenario.config, registro)
    registro.clear()

    # run.threads non incide sui risultati: S0 resta conclusa.
    troppi = _variante(scenario.config, run__threads=100000)
    run = ProjectRun(troppi, passi=_passi(registro))
    assert run.valuta().situazioni[Passo.S0].stato is StatoPasso.COMPLETATA
    impronta_s0 = run.valuta().situazioni[Passo.S0].impronta

    esito = Esecutore(run).esegui()
    assert esito.conclusione is Conclusione.ARRESTATA
    assert esito.punto.codice == "E-S0-14"
    assert esito.punto.origine == "controlli di avvio"
    assert "100000" in esito.punto.dettaglio
    # Nessuna fase eseguita, S0 compresa.
    assert registro == [] and esito.eseguite == ()
    assert run.valuta().situazioni[Passo.S0].impronta == impronta_s0


# --------------------------------------------------------------------------- #
# Riga di comando                                                              #
# --------------------------------------------------------------------------- #


@pytest.fixture
def file_config(scenario, tmp_path) -> Path:
    percorso = tmp_path / "config.yaml"
    percorso.write_text(
        yaml.safe_dump(scenario.config.model_dump(mode="json"), sort_keys=False),
        encoding="utf-8",
    )
    return percorso


@pytest.fixture
def con_doppioni(monkeypatch, registro):
    """Sostituisce le fasi realizzate con S0 vera piu' i doppioni."""
    speciali: dict[str, Any] = {}
    monkeypatch.setattr(cli, "_passi", lambda: _passi(registro, **speciali))
    return speciali


def _cli(*argomenti, capsys) -> tuple[int, str]:
    codice = cli.main(list(argomenti))
    return codice, capsys.readouterr().out


def test_cli_run_completa_con_uscita_zero(file_config, con_doppioni, capsys):
    codice, uscita = _cli("run", "--config", str(file_config), capsys=capsys)
    assert codice == cli.USCITA_SUCCESSO == 0
    assert "Esecuzione completata" in uscita


def test_cli_run_non_sovrascrive_un_esecuzione_esistente(
    file_config, scenario, con_doppioni, capsys
):
    assert _cli("run", "--config", str(file_config), capsys=capsys)[0] == 0
    radice = Path(scenario.config.io.out_root)
    prima = {p: p.stat().st_mtime_ns for p in radice.rglob("*") if p.is_file()}

    codice, uscita = _cli("run", "--config", str(file_config), capsys=capsys)
    assert codice == cli.USCITA_CONFIGURAZIONE == 3
    assert "non sovrascrive" in uscita
    assert f"amplicon16s resume --config {file_config.resolve()}" in uscita
    dopo = {p: p.stat().st_mtime_ns for p in radice.rglob("*") if p.is_file()}
    assert dopo == prima  # nessun file toccato, nemmeno i log


def test_cli_run_rifiuta_anche_una_cartella_con_altri_file(file_config, scenario, capsys):
    radice = Path(scenario.config.io.out_root)
    radice.mkdir(parents=True)
    (radice / "appunti.txt").write_text("miei")
    assert _cli("run", "--config", str(file_config), capsys=capsys)[0] == 3
    assert (radice / "appunti.txt").read_text() == "miei"


def test_cli_arresto_e_ripresa(file_config, con_doppioni, registro, capsys):
    con_doppioni["S4"] = {**FRAGILE_S4, "fallimenti": 99}
    codice, uscita = _cli("run", "--config", str(file_config), capsys=capsys)
    assert codice == cli.USCITA_ARRESTO == 4
    assert "ESECUZIONE FERMATA nella fase S4" in uscita
    assert f"amplicon16s resume --config {file_config.resolve()}" in uscita

    # Rimossa la causa, resume riparte da S4 senza rifare S0-S3.
    con_doppioni.clear()
    registro.clear()
    codice, _ = _cli("resume", "--config", str(file_config), capsys=capsys)
    assert codice == 0
    assert [p for p, _, _ in registro][0] is Passo.S4


def test_cli_validate_esegue_solo_s0(file_config, scenario, con_doppioni, registro, capsys):
    codice, uscita = _cli("validate", "--config", str(file_config), capsys=capsys)
    assert codice == 0 and "S0 conclusa" in uscita
    assert registro == []  # nessun doppione eseguito
    radice = Path(scenario.config.io.out_root)
    assert (radice / Fase.INPUT_VALIDATION.value / "manifest_S0.json").exists()

    codice, uscita = _cli("validate", "--config", str(file_config), capsys=capsys)
    assert codice == 0 and "gia' conclusa" in uscita


def test_cli_report_e_un_resoconto_provvisorio(file_config, con_doppioni, capsys):
    con_doppioni["S4"] = FRAGILE_S4
    con_doppioni["S6"] = {"degradazione": "degrada"}
    assert _cli("run", "--config", str(file_config), capsys=capsys)[0] == 0

    codice, uscita = _cli("report", "--config", str(file_config), capsys=capsys)
    assert codice == 0
    assert "RESOCONTO PROVVISORIO" in uscita
    assert "S9   disattivata" in uscita
    assert "aggiustamento [E-S4-02]: run.batch_size 24 -> 12" in uscita
    assert "degradazione [E-S6-02]" in uscita
    assert "Esecuzione completa." in uscita


def test_cli_report_senza_esecuzione(file_config, capsys):
    assert _cli("report", "--config", str(file_config), capsys=capsys)[0] == 3


def test_cli_oggi_si_ferma_alla_prima_fase_non_realizzata(file_config, capsys):
    codice, uscita = _cli("run", "--config", str(file_config), capsys=capsys)
    assert codice == cli.USCITA_FASE_NON_REALIZZATA == 5
    assert "S1, non e' ancora realizzata" in uscita


def test_cli_errore_di_configurazione(tmp_path, capsys):
    rotta = tmp_path / "rotta.yaml"
    rotta.write_text("filter: [\n")
    assert _cli("run", "--config", str(rotta), capsys=capsys)[0] == 3
    assert _cli("resume", "--config", str(tmp_path / "assente.yaml"), capsys=capsys)[0] == 3


def test_cli_whitelist_incoerente_e_errore_di_configurazione(
    scenario, tmp_path, con_doppioni, registro, capsys
):
    dati = scenario.config.model_dump(mode="json")
    dati["retry"]["whitelist"] = ["E-S4-02", "E-S14-01"]
    percorso = tmp_path / "incoerente.yaml"
    percorso.write_text(yaml.safe_dump(dati))
    codice, uscita = _cli("resume", "--config", str(percorso), capsys=capsys)
    assert codice == 3 and "E-G15-09" in uscita
    assert registro == []


def test_cli_controlli_di_risorse_nella_ripresa(scenario, tmp_path, con_doppioni, registro, capsys):
    percorso = tmp_path / "c.yaml"
    dati = scenario.config.model_dump(mode="json")
    percorso.write_text(yaml.safe_dump(dati))
    assert _cli("run", "--config", str(percorso), capsys=capsys)[0] == 0

    dati["run"]["threads"] = 100000
    percorso.write_text(yaml.safe_dump(dati))
    registro.clear()
    codice, uscita = _cli("resume", "--config", str(percorso), capsys=capsys)
    assert codice == 4 and "E-S0-14" in uscita and "controlli di avvio" in uscita
    assert registro == []


def test_il_file_di_configurazione_non_viene_mai_modificato(
    file_config, con_doppioni, capsys
):
    contenuto = file_config.read_bytes()
    istante = file_config.stat().st_mtime_ns
    con_doppioni["S4"] = FRAGILE_S4
    for comando in ("validate", "run", "resume", "report"):
        _cli(comando, "--config", str(file_config), capsys=capsys)
    assert file_config.read_bytes() == contenuto
    assert file_config.stat().st_mtime_ns == istante


def test_l_aggiustamento_non_modifica_la_configurazione_in_memoria(scenario, registro):
    run, _ = _esegui(scenario.config, registro, S4=FRAGILE_S4)
    assert run.config.run.batch_size == 24
    assert scenario.config.run.batch_size == 24


def test_un_codice_ripetibile_senza_aggiustamento_e_un_arresto_non_un_errore(
    scenario, registro
):
    """Una fase sconosciuta all'esecutore non lo fa uscire con un'eccezione."""
    _, esito = _esegui(scenario.config, registro,
                       S5={"codice": "E-S5-01", "fallimenti": 99})
    assert esito.conclusione is Conclusione.ARRESTATA
    assert esito.punto.passo is Passo.S5


def test_le_fasi_eseguite_sono_quelle_attive(scenario, registro):
    _, esito = _esegui(scenario.config, registro)
    assert [r.passo for r in esito.eseguite] == list(GRAFO.attive(scenario.config))


# --------------------------------------------------------------------------- #
# La configurazione registrata in 00_config                                   #
# --------------------------------------------------------------------------- #


def _scrivi_config(percorso: Path, config: Config, **valori: Any) -> Path:
    dati = config.model_dump(mode="json")
    for nome, v in valori.items():
        gruppo, chiave = nome.split("__")
        dati[gruppo][chiave] = v
    percorso.write_text(yaml.safe_dump(dati, sort_keys=False), encoding="utf-8")
    return percorso


def _registrate(config: Config) -> list[dict]:
    return [
        yaml.safe_load(p.read_text()) | {"_nome": p.name}
        for p in versioni_registrate(config.io.out_root)
    ]


def test_run_registra_la_configurazione_con_il_digest(
    file_config, scenario, con_doppioni, capsys
):
    assert _cli("run", "--config", str(file_config), capsys=capsys)[0] == 0
    (registrata,) = _registrate(scenario.config)
    assert registrata["_nome"] == "resolved.yaml"
    assert registrata["digest"] == risolvi(scenario.config).digest


def test_resume_con_lo_stesso_digest_non_scrive_nulla(
    file_config, scenario, con_doppioni, capsys
):
    _cli("run", "--config", str(file_config), capsys=capsys)
    cartella = Path(scenario.config.io.out_root) / Fase.CONFIG.value
    prima = {p.name: p.stat().st_mtime_ns for p in cartella.iterdir()}

    assert _cli("resume", "--config", str(file_config), capsys=capsys)[0] == 0
    assert {p.name: p.stat().st_mtime_ns for p in cartella.iterdir()} == prima


def test_resume_con_un_digest_diverso_registra_accanto_senza_sovrascrivere(
    file_config, scenario, con_doppioni, registro, tmp_path, capsys
):
    _cli("run", "--config", str(file_config), capsys=capsys)
    originale = versioni_registrate(scenario.config.io.out_root)[0].read_bytes()

    # run.threads e' fuori dall'impronta dei risultati ma non dal digest.
    registro.clear()
    diversa = _scrivi_config(tmp_path / "c2.yaml", scenario.config, run__threads=2)
    codice, uscita = _cli("resume", "--config", str(diversa), capsys=capsys)
    assert codice == 0 and registro == []  # nessuna fase rifatta
    assert "registrata in resolved_2.yaml" in uscita

    prima, seconda = _registrate(scenario.config)
    assert versioni_registrate(scenario.config.io.out_root)[0].read_bytes() == originale
    assert seconda["_nome"] == "resolved_2.yaml"
    assert seconda["precedente"] == "resolved.yaml"
    assert seconda["differenze_dalla_precedente"] == ["run.threads: 1 -> 2"]
    assert seconda["digest"] != prima["digest"]

    # Un parametro di retry e' un'altra configurazione ancora.
    terza = _scrivi_config(tmp_path / "c3.yaml", scenario.config, run__threads=2,
                           retry__max_attempts=3)
    _cli("resume", "--config", str(terza), capsys=capsys)
    assert _registrate(scenario.config)[-1]["differenze_dalla_precedente"] == [
        "retry.max_attempts: 2 -> 3"
    ]

    # Tornare alla prima e' a sua volta un cambiamento, e resta ricostruibile.
    _cli("resume", "--config", str(file_config), capsys=capsys)
    versioni = _registrate(scenario.config)
    assert [v["_nome"] for v in versioni] == [
        "resolved.yaml", "resolved_2.yaml", "resolved_3.yaml", "resolved_4.yaml"
    ]
    assert versioni[-1]["digest"] == versioni[0]["digest"]


def test_validate_su_cartella_nuova_segue_la_regola_di_run(
    file_config, scenario, con_doppioni, capsys
):
    assert _cli("validate", "--config", str(file_config), capsys=capsys)[0] == 0
    assert [v["_nome"] for v in _registrate(scenario.config)] == ["resolved.yaml"]


def test_validate_con_s0_conclusa_segue_la_regola_di_resume(
    file_config, scenario, con_doppioni, tmp_path, capsys
):
    _cli("validate", "--config", str(file_config), capsys=capsys)
    _cli("validate", "--config", str(file_config), capsys=capsys)
    assert len(_registrate(scenario.config)) == 1

    diversa = _scrivi_config(tmp_path / "c2.yaml", scenario.config, retry__enabled=False)
    codice, uscita = _cli("validate", "--config", str(diversa), capsys=capsys)
    assert codice == 0 and "gia' conclusa" in uscita
    assert [v["_nome"] for v in _registrate(scenario.config)] == [
        "resolved.yaml", "resolved_2.yaml"
    ]


def test_un_arresto_ai_controlli_di_avvio_non_registra_nulla(
    file_config, scenario, con_doppioni, tmp_path, capsys
):
    _cli("run", "--config", str(file_config), capsys=capsys)
    troppi = _scrivi_config(tmp_path / "c2.yaml", scenario.config, run__threads=100000)
    assert _cli("resume", "--config", str(troppi), capsys=capsys)[0] == 4
    assert len(_registrate(scenario.config)) == 1


def test_il_log_dice_con_quale_configurazione_si_esegue(
    file_config, scenario, con_doppioni, capsys
):
    _cli("run", "--config", str(file_config), capsys=capsys)
    _cli("resume", "--config", str(file_config), capsys=capsys)
    righe = [
        json.loads(r)
        for r in (Path(scenario.config.io.out_root) / Fase.LOGS.value / NOME_FILE_LOG)
        .read_text().splitlines()
    ]
    in_uso = [r for r in righe if r["messaggio"] == "configurazione in uso"]
    assert [(r["file"], r["registrata_ora"]) for r in in_uso] == [
        ("resolved.yaml", True), ("resolved.yaml", False)
    ]


def test_le_versioni_si_ordinano_per_numero(tmp_path):
    cartella = tmp_path / Fase.CONFIG.value
    cartella.mkdir()
    for nome in ("resolved_10.yaml", "resolved.yaml", "resolved_2.yaml", "resolved_9.yaml",
                 "manifest.json", "resolved_x.yaml"):
        (cartella / nome).write_text("digest: x\n")
    assert [p.name for p in versioni_registrate(tmp_path)] == [
        "resolved.yaml", "resolved_2.yaml", "resolved_9.yaml", "resolved_10.yaml"
    ]
