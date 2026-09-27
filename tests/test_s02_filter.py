"""Test della fase S2, del retry senza modifiche, della rimozione dei file
filtrati e del tracciamento delle letture.

La prima parte non richiede R: prova i meccanismi con S0 vera e fasi
doppione. La seconda esegue S2 davvero, sulla versione ridotta del
sottoinsieme di prova, con dada2: si salta senza Bioconductor e nella catena
di integrazione continua gira dentro l'immagine, dove non puo' saltarsi.
"""

from __future__ import annotations

import gzip
import json
import os
import shutil
from pathlib import Path
from typing import Any, ClassVar

import pytest
import yaml
from conftest import NEGATIVO, POSITIVO, Campione, crea_scenario
from sottoinsieme import RIDOTTO, config_ridotta, motivo_pacchetti_r_assenti

from amplicon16s.config.resolve import PARAMETRI_SENZA_EFFETTO, risolvi
from amplicon16s.config.schema import Config, valida
from amplicon16s.errors.exceptions import ErroreRevisioneUmana, errore
from amplicon16s.io_layer.artifacts import Fase, ManifestoPasso, nome_registro_rimozioni
from amplicon16s.logging.logger import chiudi
from amplicon16s.metadata.models import ClasseCampione
from amplicon16s.runner.executor import Conclusione, Esecutore
from amplicon16s.runner.graph import GRAFO, Passo
from amplicon16s.runner.project import ProjectRun, StatoPasso
from amplicon16s.runner.retry import Aggiustamento, Motivo, PoliticaRetry, senza_modifiche
from amplicon16s.runner.tracciamento import ricomponi
from amplicon16s.steps.base import PipelineStep, Produzione, StepContext
from amplicon16s.steps.s00_validate import ValidazioneIngressi
from amplicon16s.steps.s02_filter import (
    CLASSI_CONTROLLATE,
    archivio_incompleto,
    controlla_filtro,
)

ESEMPIO = Path(__file__).resolve().parents[1] / "config" / "config.example.yaml"
TUTTE = tuple(Passo)


@pytest.fixture(autouse=True)
def uscite_pulite():
    chiudi()
    yield
    chiudi()


# --------------------------------------------------------------------------- #
# Parametri                                                                    #
# --------------------------------------------------------------------------- #


def test_i_parametri_nuovi_hanno_i_valori_del_piano():
    esempio = yaml.safe_load(ESEMPIO.read_text(encoding="utf-8"))
    assert esempio["qc"]["max_zeroed_samples"] == 0
    assert esempio["qc"]["max_frac_lost_filter"] == 0.30
    assert esempio["run"]["keep_filtered_fastq"] is True


def test_conservare_i_file_filtrati_non_incide_sui_risultati(tmp_path):
    assert "run.keep_filtered_fastq" in PARAMETRI_SENZA_EFFETTO
    dati = yaml.safe_load(ESEMPIO.read_text(encoding="utf-8"))
    prima = risolvi(valida(dati)).impronta_risultati
    dati["run"]["keep_filtered_fastq"] = False
    assert risolvi(valida(dati)).impronta_risultati == prima


def test_i_negativi_sono_esclusi_dai_controlli_del_filtro():
    assert {c.value for c in CLASSI_CONTROLLATE} == {"biologico", "controllo_positivo"}


# --------------------------------------------------------------------------- #
# Controlli sul risultato, per classe                                          #
# --------------------------------------------------------------------------- #

B, P, N = ClasseCampione.BIOLOGICO, ClasseCampione.CONTROLLO_POSITIVO, ClasseCampione.CONTROLLO_NEGATIVO
CLASSI = {"b1": B, "b2": B, "p1": P, "n1": N, "n2": N}
QC = valida(yaml.safe_load(ESEMPIO.read_text(encoding="utf-8"))).qc


def _controlla(uscita: dict[str, int], ingresso: int = 1000):
    return controlla_filtro({c: ingresso for c in CLASSI}, uscita, CLASSI, QC)


def test_un_bianco_azzerato_non_ferma():
    metriche = _controlla({"b1": 990, "b2": 985, "p1": 980, "n1": 0, "n2": 970})
    assert metriche["azzerati"]["controllo_negativo"] == ["n1"]


@pytest.mark.parametrize("azzerato", ["b1", "p1"])
def test_un_biologico_o_un_positivo_azzerato_ferma(azzerato):
    uscita = {"b1": 990, "b2": 985, "p1": 980, "n1": 970, "n2": 970, azzerato: 0}
    with pytest.raises(ErroreRevisioneUmana) as info:
        _controlla(uscita)
    assert info.value.codice == "E-S2-01"
    assert info.value.contesto["campioni"] == [azzerato]


def test_la_perdita_dei_bianchi_non_conta_per_e_s2_02():
    metriche = _controlla({"b1": 990, "b2": 985, "p1": 980, "n1": 100, "n2": 50})
    assert metriche["perdita_media_controllata"] < QC.max_frac_lost_filter


def test_una_perdita_media_oltre_il_30_per_cento_ferma():
    with pytest.raises(ErroreRevisioneUmana) as info:
        _controlla({"b1": 600, "b2": 650, "p1": 700, "n1": 990, "n2": 990})
    assert info.value.codice == "E-S2-02"


# --------------------------------------------------------------------------- #
# Retry senza modifiche                                                        #
# --------------------------------------------------------------------------- #


def test_un_nuovo_tentativo_senza_modifiche_va_dichiarato_e_motivato():
    with pytest.raises(ValueError, match="motivato"):
        senza_modifiche("  ")
    with pytest.raises(ValueError, match="regola"):
        Aggiustamento(None, lambda v: v, "incoerente")
    assert senza_modifiche("errore transitorio").invariato


def test_la_politica_ritenta_senza_modifiche_solo_se_dichiarato(tmp_path):
    config = _scenario(tmp_path).config
    politica = PoliticaRetry.da_config(config)
    assert politica.decidi("E-S2-03", 1, senza_modifiche("transitorio"), config).ritenta
    # Senza dichiarazione la regola generale resta: non si ritenta identico.
    assert politica.decidi("E-S2-03", 1, None, config).motivo is Motivo.NESSUNA_AZIONE
    # E i tentativi restano quelli ammessi.
    assert not politica.decidi("E-S2-03", 2, senza_modifiche("transitorio"), config).ritenta


# --------------------------------------------------------------------------- #
# Doppioni                                                                     #
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


def _scenario(tmp_path, **sovrascrivi):
    return crea_scenario(tmp_path, _campioni(), con_letture=True, sovrascrivi=sovrascrivi or None)


class Doppione(PipelineStep):
    """Scrive un artefatto; fallisce per i primi ``fallimenti`` tentativi."""

    codice: ClassVar[str | None] = None
    fallimenti: ClassVar[int] = 0
    #: Nome dell'artefatto; per il doppione di S2, un file "filtrato".
    nome: ClassVar[str | None] = None

    def __init__(self, registro: list) -> None:
        self.registro = registro
        self.tentativi = 0

    def calcola(self, contesto: StepContext) -> Produzione:
        self.registro.append(self.passo)
        self.tentativi += 1
        if self.codice and self.tentativi <= self.fallimenti:
            raise errore(self.codice, f"doppione, tentativo {self.tentativi}")
        nome = self.nome or f"{str(self.passo).lower()}.txt"
        artefatto = contesto.albero.scrivi_testo(self.cartella, nome, f"{self.passo}")
        return Produzione((artefatto,))


class FiltroDoppione(Doppione):
    """Un S2 finto con le letture filtrate temporanee, come FiltroLetture."""

    passo: ClassVar[Passo] = Passo.S2
    nome: ClassVar[str] = "ERX3000001_filt.fastq.gz"

    def artefatti_temporanei(self, manifesto: ManifestoPasso, config: Config) -> tuple[str, ...]:
        if config.run.keep_filtered_fastq:
            return ()
        return tuple(n for n in manifesto.nomi if n.endswith("_filt.fastq.gz"))


def _passi(registro, **speciali: dict[str, Any]) -> dict[Passo, PipelineStep]:
    passi: dict[Passo, PipelineStep] = {Passo.S0: ValidazioneIngressi()}
    for passo in TUTTE[1:]:
        base = FiltroDoppione if passo is Passo.S2 else Doppione
        attributi = {"passo": passo, **speciali.get(str(passo), {})}
        passi[passo] = type(f"Doppione{passo}", (base,), attributi)(registro)
    return passi


def _esegui(config, registro, **speciali):
    run = ProjectRun(config, passi=_passi(registro, **speciali))
    return run, Esecutore(run).esegui()


# --------------------------------------------------------------------------- #
# E-S2-03 ritentato senza modifiche                                            #
# --------------------------------------------------------------------------- #

TRANSITORIO = {
    "codice": "E-S2-03", "fallimenti": 1,
    "aggiustamenti": {"E-S2-03": senza_modifiche("errore di lettura transitorio")},
}


def test_e_s2_03_si_ritenta_senza_modifiche_e_si_registra(tmp_path):
    registro: list = []
    run, esito = _esegui(_scenario(tmp_path).config, registro, S2=TRANSITORIO)

    assert esito.conclusione is Conclusione.COMPLETATA
    assert registro.count(Passo.S2) == 2
    (voce,) = run.albero.manifesto_passo(Passo.S2, Fase.FILTERED).aggiustamenti
    assert voce["codice"] == "E-S2-03"
    assert voce["parametro"] is None and voce["usato"] is None
    assert "senza modifiche" in voce["azione"]


def test_un_errore_di_lettura_che_persiste_si_ferma_dopo_i_tentativi(tmp_path):
    registro: list = []
    persistente = {**TRANSITORIO, "fallimenti": 99}
    _, esito = _esegui(_scenario(tmp_path).config, registro, S2=persistente)

    assert esito.conclusione is Conclusione.ARRESTATA
    assert esito.punto.codice == "E-S2-03"
    assert esito.punto.tentativi == 2
    assert registro.count(Passo.S2) == 2


# --------------------------------------------------------------------------- #
# run.keep_filtered_fastq                                                      #
# --------------------------------------------------------------------------- #


def _filtrato(run) -> Path:
    return run.albero.cartella(Fase.FILTERED) / "ERX3000001_filt.fastq.gz"


def test_per_difetto_i_file_filtrati_restano(tmp_path):
    registro: list = []
    run, esito = _esegui(_scenario(tmp_path).config, registro)
    assert esito.conclusione is Conclusione.COMPLETATA
    assert _filtrato(run).exists()


def test_rimossi_a_fine_esecuzione_senza_rifare_nulla(tmp_path):
    registro: list = []
    config = _scenario(tmp_path, run={"keep_filtered_fastq": False}).config
    run, esito = _esegui(config, registro)

    assert esito.conclusione is Conclusione.COMPLETATA
    assert not _filtrato(run).exists()
    registro_rimozioni = run.albero.cartella(Fase.FILTERED) / nome_registro_rimozioni("S2")
    documento = json.loads(registro_rimozioni.read_text())
    assert [v["nome"] for v in documento["rimossi"]] == ["ERX3000001_filt.fastq.gz"]

    # La rimozione e' intenzionale: S2 resta conclusa e la ripresa non fa nulla.
    registro.clear()
    ripresa = ProjectRun(config, passi=_passi(registro))
    valutazione = ripresa.valuta()
    assert valutazione.completa
    assert "rimossi di proposito" in valutazione.situazioni[Passo.S2].motivo
    assert Esecutore(ripresa).esegui().eseguite == ()
    assert registro == []


def test_un_file_filtrato_perso_senza_registro_fa_rifare_s2(tmp_path):
    registro: list = []
    run, _ = _esegui(_scenario(tmp_path).config, registro)
    _filtrato(run).unlink()
    situazione = run.valuta().situazioni[Passo.S2]
    assert situazione.stato is StatoPasso.DA_ESEGUIRE
    assert "mancanti" in situazione.motivo


def test_se_una_fase_che_li_legge_va_rifatta_s2_torna_da_eseguire(tmp_path):
    registro: list = []
    config = _scenario(tmp_path, run={"keep_filtered_fastq": False}).config
    run, _ = _esegui(config, registro)
    run.albero.rimuovi_manifesto_passo(Passo.S3, Fase.ERROR_MODELS)

    situazione = run.valuta().situazioni[Passo.S2]
    assert situazione.stato is StatoPasso.DA_ESEGUIRE
    assert "richiesti da S3" in situazione.motivo

    registro.clear()
    esito = Esecutore(ProjectRun(config, passi=_passi(registro))).esegui()
    assert esito.conclusione is Conclusione.COMPLETATA
    rifatte = [r.passo for r in esito.eseguite]
    assert rifatte[:2] == [Passo.S2, Passo.S3]
    assert Passo.S0 not in rifatte and Passo.S1 not in rifatte


def test_nessuna_rimozione_prima_della_fine(tmp_path):
    registro: list = []
    config = _scenario(tmp_path, run={"keep_filtered_fastq": False}).config
    fragile = {"codice": "E-S6-01", "fallimenti": 99}
    run, esito = _esegui(config, registro, S6=fragile)
    assert esito.conclusione is Conclusione.ARRESTATA
    assert _filtrato(run).exists()

    # Ne' con un'esecuzione limitata a una parte del grafo.
    registro.clear()
    run = ProjectRun(_scenario(tmp_path / "b", run={"keep_filtered_fastq": False}).config,
                     passi=_passi(registro))
    Esecutore(run, fino_a=Passo.S2).esegui()
    assert _filtrato(run).exists()


# --------------------------------------------------------------------------- #
# Integrita' degli archivi                                                     #
# --------------------------------------------------------------------------- #


def test_un_archivio_troncato_non_si_decomprime_per_intero(tmp_path):
    originale = next((RIDOTTO / "fastq").glob("ERX12084006_*"))
    troncato = tmp_path / originale.name
    dati = originale.read_bytes()
    troncato.write_bytes(dati[: len(dati) * 7 // 10])
    assert archivio_incompleto(originale) is None
    assert archivio_incompleto(troncato) is not None


# --------------------------------------------------------------------------- #
# S2 vera, sulla versione ridotta                                              #
# --------------------------------------------------------------------------- #

_MOTIVO_ASSENTI = motivo_pacchetti_r_assenti("dada2", "ShortRead", "Biostrings", "jsonlite")


@pytest.fixture
def dada2():
    """Richiede R con dada2: salta senza, ma in CI fallisce."""
    if _MOTIVO_ASSENTI is not None:
        if os.environ.get("AMPLICON16S_RICHIEDI_BIOC") == "1":
            pytest.fail(f"Bioconductor e' richiesto in questo ambiente: {_MOTIVO_ASSENTI}")
        pytest.skip(_MOTIVO_ASSENTI)


@pytest.fixture(scope="module")
def filtrata(tmp_path_factory):
    """S0, S1 e S2 sulla versione ridotta, una volta sola per i test che leggono."""
    if _MOTIVO_ASSENTI is not None:
        return None
    run = ProjectRun(config_ridotta(tmp_path_factory.mktemp("s2")))
    return run, Esecutore(run, fino_a=Passo.S2).esegui()


def _conteggi(percorso: Path) -> dict[str, int]:
    import csv

    with open(percorso, encoding="utf-8") as file:
        return {r["campione"]: int(float(r["letture"])) for r in csv.DictReader(file, delimiter="\t")}


def test_s2_produce_file_filtrati_e_tracciamento_con_il_suo_manifesto(dada2, filtrata):
    run, esito = filtrata
    assert esito.conclusione is Conclusione.COMPLETATA
    manifesto = run.albero.manifesto_passo(Passo.S2, Fase.FILTERED)
    nomi = set(manifesto.nomi)
    assert {"letture_prefiltro.tsv", "letture_filtrate.tsv"} <= nomi
    filtrati = {n for n in nomi if n.endswith("_filt.fastq.gz")}
    assert len(filtrati) == 28  # nessun campione azzerato sulla versione ridotta


def test_tutte_le_letture_filtrate_hanno_137_basi(dada2, filtrata):
    run, _ = filtrata
    cartella = run.albero.cartella(Fase.FILTERED)
    uscita = _conteggi(cartella / "letture_filtrate.tsv")
    for percorso in cartella.glob("*_filt.fastq.gz"):
        with gzip.open(percorso, "rt") as file:
            lunghezze = {len(r.rstrip("\n")) for i, r in enumerate(file) if i % 4 == 1}
        assert lunghezze == {137}, percorso.name
    contate = {}
    for percorso in cartella.glob("*_filt.fastq.gz"):
        with gzip.open(percorso, "rt") as file:
            contate[percorso.name.split("_")[0]] = sum(1 for _ in file) // 4
    assert contate == {c: n for c, n in uscita.items() if n > 0}


def test_le_letture_in_ingresso_sono_quelle_contate_da_s1(dada2, filtrata):
    run, _ = filtrata
    grezze = _conteggi(run.albero.cartella(Fase.QC_PROFILES) / "letture_grezze.tsv")
    prefiltro = _conteggi(run.albero.cartella(Fase.FILTERED) / "letture_prefiltro.tsv")
    assert prefiltro == grezze


def test_il_tracciamento_si_ricompone_senza_toccare_i_file_di_s2(dada2, filtrata):
    run, _ = filtrata
    manifesto_s2 = run.albero.manifesto_passo(Passo.S2, Fase.FILTERED)

    # Una fase successiva registra il proprio passo, solo nei propri artefatti.
    class DenoiseDoppione(PipelineStep):
        passo: ClassVar[Passo] = Passo.S3
        passi_tracciamento: ClassVar[tuple[str, ...]] = ("modello",)

        def calcola(self, contesto: StepContext) -> Produzione:
            filtrate = _conteggi(contesto.albero.cartella(Fase.FILTERED) / "letture_filtrate.tsv")
            righe = ["campione\tpasso\tletture"] + [f"{c}\tmodello\t{n}" for c, n in filtrate.items()]
            artefatto = contesto.albero.scrivi_testo(
                self.cartella, "letture_modello.tsv", "\n".join(righe) + "\n"
            )
            return Produzione((artefatto,))

    passi = {**run.passi, Passo.S3: DenoiseDoppione()}
    successiva = ProjectRun(run.config, passi=passi)
    successiva.fase(Passo.S3).esegui(successiva.contesto(Passo.S3))

    tracciamento = ricomponi(successiva)
    assert tracciamento.passi == ("grezze", "prefiltro", "filtrate", "modello")
    assert tracciamento.origine == {
        "grezze": "S1", "prefiltro": "S2", "filtrate": "S2", "modello": "S3",
    }
    # S2 e' intatta: stesso manifesto, artefatti integri, ancora conclusa.
    assert successiva.albero.manifesto_passo(Passo.S2, Fase.FILTERED) == manifesto_s2
    assert successiva.albero.non_integri_del_passo(manifesto_s2) == ()
    assert successiva.valuta().situazioni[Passo.S2].stato is StatoPasso.COMPLETATA


def test_con_max_ee_molto_restrittivo_il_filtro_si_ferma(dada2, tmp_path):
    run = ProjectRun(config_ridotta(tmp_path, filter={"maxEE": 0.01}))
    esito = Esecutore(run, fino_a=Passo.S2).esegui()
    assert esito.conclusione is Conclusione.ARRESTATA
    assert esito.punto.passo is Passo.S2
    assert esito.punto.codice in {"E-S2-01", "E-S2-02"}
    assert run.albero.manifesto_passo(Passo.S2, Fase.FILTERED) is None


def test_un_archivio_corrotto_si_ferma_dopo_i_tentativi(dada2, tmp_path):
    """Troncato oltre le letture che S0 ispeziona: S0 passa, S2 se ne accorge."""
    fastq = tmp_path / "fastq"
    shutil.copytree(RIDOTTO / "fastq", fastq)
    bersaglio = next(fastq.glob("ERX12084006_*"))
    dati = bersaglio.read_bytes()
    bersaglio.write_bytes(dati[: len(dati) * 7 // 10])

    config = config_ridotta(tmp_path, io={"fastq_dir": str(fastq)}, qc={"head_reads": 20})
    run = ProjectRun(config)
    esito = Esecutore(run, fino_a=Passo.S2).esegui()

    assert esito.conclusione is Conclusione.ARRESTATA
    assert esito.punto.passo is Passo.S2
    assert esito.punto.codice == "E-S2-03"
    assert esito.punto.tentativi == 2
    assert "ERX12084006" in esito.punto.dettaglio
    assert Passo.S0 in {r.passo for r in esito.eseguite}
    assert GRAFO.nodo(Passo.S2).cartella is Fase.FILTERED


# --------------------------------------------------------------------------- #
# Dataset completo, nel container                                              #
# --------------------------------------------------------------------------- #


@pytest.mark.dati_reali
def test_s2_sul_dataset_completo(dada2, tmp_path):
    """Nessun campione azzerato, perdita per classe misurata, tutte a 137 basi."""
    from concurrent.futures import ThreadPoolExecutor
    import logging

    from amplicon16s.config.schema import carica

    percorso = os.environ.get("AMPLICON16S_CONFIG_DATI_REALI")
    if not percorso or not Path(percorso).expanduser().is_file():
        pytest.skip("AMPLICON16S_CONFIG_DATI_REALI non impostata o file assente")
    dati = carica(Path(percorso).expanduser()).model_dump(mode="python")
    dati["io"]["out_root"] = str(tmp_path / "out")
    logging.getLogger("amplicon16s").setLevel(logging.WARNING)

    run = ProjectRun(valida(dati))
    esito = Esecutore(run, fino_a=Passo.S2).esegui()
    assert esito.conclusione is Conclusione.COMPLETATA, esito.punto and esito.punto.testo()
    s2 = next(r for r in esito.eseguite if r.passo is Passo.S2)
    print("\nS2 sul dataset completo:", json.dumps(dict(s2.metriche), indent=1), f"{s2.secondi} s")

    assert s2.metriche["campioni"] == 960
    assert all(not v for v in s2.metriche["azzerati"].values())
    assert s2.metriche["perdita_media_controllata"] < 0.05

    def lunghezze(percorso: Path) -> set[int]:
        with gzip.open(percorso, "rt") as file:
            return {len(r.rstrip("\n")) for i, r in enumerate(file) if i % 4 == 1}

    file = list(run.albero.cartella(Fase.FILTERED).glob("*_filt.fastq.gz"))
    assert len(file) == 960
    with ThreadPoolExecutor(max_workers=8) as esecutore:
        assert set().union(*esecutore.map(lunghezze, file)) == {137}
