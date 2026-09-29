r"""Suite di test per la Fase S3 (modello parametrico d'errore per corsa di sequenziamento).

1. Inquadramento nel Piano Operativo
------------------------------------
Settimana 13 (W13), Fase F4 (Fasi di calcolo R/Bioconductor: Fase S3, Modello
d'errore per corsa di sequenziamento in ``04_error_models/``).

2. Moduli sorgente coperti
--------------------------
* ``src/amplicon16s/steps/s03_learn_errors.py``
* ``R/03_learn_errors.R``
* ``src/amplicon16s/runner/retry.py``
* ``src/amplicon16s/runner/executor.py``

3. Cosa valuta questo file
--------------------------
Verifica l'intero contratto scientifico e sistemistico della Fase S3:
- stima dei modelli separati per corsa di sequenziamento (``err.batch_column``)
  e fallback a modello singolo (``modello_tutti.rds``) con colonna nulla o
  senza tabella di arricchimento del lotto (``io.batch_table``);
- determinismo dell'ordine di campionamento (``run.seed``) e riproducibilita'
  crittografica byte per byte degli artefatti (inclusi i grafici PNG generati
  con dispositivo grafico ``cairo`` privo di metadati temporali);
- determinazione della convergenza di ``dada2::learnErrors`` tramite confronto
  matriciale ``identical(err_in, err_out)``;
- gestione di ``E-S3-01`` con retry automatico a ``err.nbases`` raddoppiato
  oppure soppressione motivata del retry (``RITENTARE_INUTILE``) quando le basi
  disponibili nella corsa sono gia' tutte utilizzate;
- gestione del codice ``E-S3-02`` quando ``err.batch_column`` e' attivo ma un
  campione risulta privo di corsa associata.

4. Comandi Bash e scenari di esecuzione:
    1. Modalita locale standard (senza Bioconductor R):
       pytest tests/test_w13_s03_batch.py -v
       Risultato atteso: 20 test (12 passed in Python, 8 skipped per assenza
       di dada2/ggplot2/ShortRead in R locale e dei dati reali).

    2. Modalita container Docker standard (subset ridotto con Bioconductor):
       docker run --rm \
         -e PYTHONPATH=/app/src \
         -v "$(pwd)":/app \
         -w /app \
         amplicon16s:dev \
         pytest -o cache_dir=/tmp/.pytest_cache tests/test_w13_s03_batch.py -v
       Risultato atteso: 19 passed, 1 skipped in ~90s (resta saltato solo il
       test sui 960 file FASTQ reali).

    3. Modalita container Docker completa (100% verde con dati reali OSD-734):
       docker run --rm \
         -e PYTHONPATH=/app/src \
         -e AMPLICON16S_CONFIG_DATI_REALI=/home/nemo/ASI/config_osd734.yaml \
         -v "$(pwd)":/app \
         -v /home/nemo/ASI:/home/nemo/ASI \
         -w /app \
         amplicon16s:dev \
         pytest -o cache_dir=/tmp/.pytest_cache tests/test_w13_s03_batch.py -v
       Risultato atteso: 20 passed in ~15 minuti.

    Accorgimenti operativi per il container Docker:
    - Impostare '-e PYTHONPATH=/app/src' per caricare i moduli aggiornati da /app/src.
    - Usare '-o cache_dir=/tmp/.pytest_cache' per proteggere i permessi della cartella locale.
    - Montare '-v /home/nemo/ASI:/home/nemo/ASI' per rendere accessibili i 2.4 GB di dati reali.

5. Risultato atteso
-------------------
20 test totali (12 passed, 8 skipped in ~0.80s in ambiente locale privo di
Bioconductor e dei dati reali; 19 passed, 1 skipped nel container CI;
20 passed nel container con ``AMPLICON16S_CONFIG_DATI_REALI``).

6. Razionale scientifico e sistemistico
---------------------------------------
- Corse di sequenziamento Illumina MiSeq distinte presentano profili fisici di
  errore differenti: stimare un modello parametrico separato per ciascuna corsa
  (``3DMM_Plates_1-5_S1_L001`` e ``3DMM_rerun_Plates_6-10_S1_L001``) evita medie
  spurie che altererebbero la chiamata delle ASV in S4.
- L'uso del dispositivo ``cairo`` per i PNG diagnostici e la pre-permutazione
  deterministica dei campioni in Python con ``run.seed`` garantiscono che due
  esecuzioni da zero producano in ``04_error_models/`` artefatti identici byte
  per byte, verificabili tramite checksum.
"""

from __future__ import annotations

import dataclasses
import hashlib
import json
import os
import shutil
import subprocess
from pathlib import Path
from typing import ClassVar

import pytest
from sottoinsieme import config_ridotta, motivo_pacchetti_r_assenti, selezione

from amplicon16s.config.schema import valida
from amplicon16s.errors.exceptions import (
    ErroreRevisioneUmana,
    ErroreRitentabileConRevisione,
    errore,
)
from amplicon16s.io_layer.artifacts import Fase
from amplicon16s.logging.logger import chiudi
from amplicon16s.metadata.models import Inventario
from amplicon16s.rbridge.runner import trova_rscript
from amplicon16s.runner.executor import Conclusione, Esecutore
from amplicon16s.runner.graph import Passo
from amplicon16s.runner.project import ProjectRun, StatoPasso
from amplicon16s.runner.retry import RITENTARE_INUTILE, Motivo, PoliticaRetry, raddoppia
from amplicon16s.steps.base import PipelineStep, Produzione, StepContext
from amplicon16s.steps.s00_validate import esegui_s0
from amplicon16s.steps.s03_learn_errors import MODELLO_UNICO, nomi_sicuri, pianifica

CORSA_A = "3DMM_Plates_1-5_S1_L001"
CORSA_B = "3DMM_rerun_Plates_6-10_S1_L001"


@pytest.fixture(autouse=True)
def uscite_pulite():
    """Chiude le uscite del log prima e dopo ogni test, perché nessun handler resti
    aperto sulla cartella temporanea.
    """
    chiudi()
    yield
    chiudi()


# --------------------------------------------------------------------------- #
# Nomi sicuri                                                                  #
# --------------------------------------------------------------------------- #


def test_i_nomi_delle_corse_diventano_nomi_di_file_sicuri():
    """
    **Obiettivo**: Verificare che i nomi delle corse, usati nei nomi dei file
    dei modelli, siano ridotti a caratteri sicuri, senza collisioni e senza
    occupare il nome riservato al modello unico.

    **Razionale scientifico e sistemistico**: Un nome di corsa con barre o spazi
    produrrebbe percorsi sbagliati; due corse che si riducono allo stesso nome
    si sovrascriverebbero il modello a vicenda.
    """
    nomi = nomi_sicuri([CORSA_A, "corsa/1 con spazi", "a b", "a_b", "tutti", ".."])
    assert nomi[CORSA_A] == CORSA_A
    assert nomi["corsa/1 con spazi"] == "corsa_1_con_spazi"
    assert {nomi["a b"], nomi["a_b"]} == {"a_b", "a_b_2"}
    assert nomi["tutti"] != MODELLO_UNICO
    assert nomi[".."] == "corsa"
    assert len(set(nomi.values())) == len(nomi)


# --------------------------------------------------------------------------- #
# Pianificazione dei modelli                                                   #
# --------------------------------------------------------------------------- #


@pytest.fixture(scope="module")
def inventario_ridotto(tmp_path_factory) -> Inventario:
    """L'inventario della versione ridotta, da S0 vera (che non richiede R)."""
    cartella = tmp_path_factory.mktemp("s0")
    config = config_ridotta(cartella)
    esegui_s0(config)
    inventario = ProjectRun(config).inventario
    assert inventario is not None
    return inventario


def _letture(inventario: Inventario, quante: int = 1000) -> dict[str, int]:
    """Lo stesso numero di letture filtrate per ogni campione dell'inventario."""
    return {c.accession: quante for c in inventario}


def test_un_modello_per_corsa(inventario_ridotto, tmp_path):
    """
    **Obiettivo**: Verificare che con ``err.batch_column`` attivo si
    pianifichino due modelli, uno per ciascuna corsa, e che ogni campione
    appartenga al modello della propria corsa.

    **Razionale scientifico e sistemistico**: Corse diverse hanno profili
    d'errore diversi; un modello unico li medierebbe.
    """
    modelli = pianifica(inventario_ridotto, _letture(inventario_ridotto), 137, config_ridotta(tmp_path))
    assert [m.nome for m in modelli] == [CORSA_A, CORSA_B]
    corse = {r["accession"]: r["corsa"] for r in selezione()}
    for modello in modelli:
        assert {corse[c] for c in modello.campioni} == {modello.corsa}
    assert sum(len(m.campioni) for m in modelli) == 28


@pytest.mark.parametrize(
    "sovrascrivi",
    [{"err": {"batch_column": None}}, {"io": {"batch_table": None}}],
    ids=["colonna-nulla", "senza-file-del-lotto"],
)
def test_senza_corsa_un_solo_modello(inventario_ridotto, tmp_path, sovrascrivi):
    """
    **Obiettivo**: Verificare che con ``err.batch_column`` nullo, o senza il
    file di arricchimento del lotto, si pianifichi un solo modello su tutti i
    campioni.

    **Razionale scientifico e sistemistico**: La pipeline deve restare
    utilizzabile senza il file di arricchimento del lotto.
    """
    config = config_ridotta(tmp_path, **sovrascrivi)
    (modello,) = pianifica(inventario_ridotto, _letture(inventario_ridotto), 137, config)
    assert modello.nome == MODELLO_UNICO
    assert modello.corsa is None
    assert len(modello.campioni) == 28


def test_un_campione_senza_corsa_con_la_colonna_attiva_ferma(inventario_ridotto, tmp_path):
    """
    **Obiettivo**: Verificare che un campione senza corsa, con
    ``err.batch_column`` attivo, sollevi ``E-S3-02``.

    **Razionale scientifico e sistemistico**: Un campione senza corsa non ha un
    modello d'errore da cui farsi correggere; assegnarlo a caso sarebbe una
    scelta silenziosa.
    """
    campioni = list(inventario_ridotto)
    campioni[0] = dataclasses.replace(campioni[0], corsa=None)
    with pytest.raises(ErroreRevisioneUmana) as info:
        pianifica(Inventario(tuple(campioni)), _letture(inventario_ridotto), 137, config_ridotta(tmp_path))
    assert info.value.codice == "E-S3-02"


def test_l_ordine_dei_campioni_viene_dal_seme(inventario_ridotto, tmp_path):
    """
    **Obiettivo**: Verificare che l'ordine in cui i campioni entrano nella stima
    dipenda da ``run.seed`` quando ``err.randomize`` e' vero, e sia l'ordine
    delle accession quando e' falso.

    **Razionale scientifico e sistemistico**: Due esecuzioni identiche devono
    stimare il modello sugli stessi campioni.
    """
    letture = _letture(inventario_ridotto)
    primo = pianifica(inventario_ridotto, letture, 137, config_ridotta(tmp_path))
    stesso = pianifica(inventario_ridotto, letture, 137, config_ridotta(tmp_path))
    altro = pianifica(inventario_ridotto, letture, 137, config_ridotta(tmp_path, run={"seed": 7}))
    fisso = pianifica(inventario_ridotto, letture, 137, config_ridotta(tmp_path, err={"randomize": False}))
    assert [m.ordine for m in primo] == [m.ordine for m in stesso]
    assert [m.ordine for m in primo] != [m.ordine for m in altro]
    assert [m.ordine for m in fisso] == [tuple(sorted(m.campioni)) for m in fisso]


def test_la_stima_si_ferma_oltre_err_nbases(inventario_ridotto, tmp_path):
    """
    **Obiettivo**: Verificare che i campioni usati siano i primi dell'ordine,
    finche' le basi superano ``err.nbases``, e che con ``err.nbases`` oltre le
    basi della corsa si usino tutte.

    **Razionale scientifico e sistemistico**: E' la stessa regola di
    ``learnErrors``; saperla applicare fuori da R dice su quali campioni e'
    stata fatta la stima, e se raddoppiare ``err.nbases`` servirebbe.
    """
    letture = _letture(inventario_ridotto)  # 1000 letture da 137 basi ciascuna
    pochi = pianifica(inventario_ridotto, letture, 137, config_ridotta(tmp_path, err={"nbases": 300000}))
    for modello in pochi:
        assert len(modello.usati) == 3  # 137000, 274000, 411000 > 300000
        assert modello.usati == modello.ordine[:3]
        assert not modello.tutte_usate
    tanti = pianifica(inventario_ridotto, letture, 137, config_ridotta(tmp_path))
    assert all(m.tutte_usate for m in tanti)
    assert all(m.basi_usate == m.basi_disponibili for m in tanti)


def test_un_campione_senza_letture_filtrate_non_entra_nella_stima(inventario_ridotto, tmp_path):
    """
    **Obiettivo**: Verificare che un campione azzerato dal filtro resti nel
    modello della sua corsa ma non entri nella stima.

    **Razionale scientifico e sistemistico**: Non ha un file filtrato da cui
    stimare; la corrispondenza gli assegna comunque il modello della corsa.
    """
    letture = _letture(inventario_ridotto)
    letture["ERX12083285"] = 0
    modelli = pianifica(inventario_ridotto, letture, 137, config_ridotta(tmp_path))
    (modello,) = [m for m in modelli if "ERX12083285" in m.campioni]
    assert "ERX12083285" not in modello.ordine


# --------------------------------------------------------------------------- #
# Retry dichiarato inutile                                                     #
# --------------------------------------------------------------------------- #


def test_la_politica_non_ritenta_se_la_fase_lo_dichiara_inutile(tmp_path):
    """
    **Obiettivo**: Verificare che un errore che dichiara inutile l'azione
    correttiva non venga ritentato, con il motivo nell'arresto.

    **Razionale scientifico e sistemistico**: Raddoppiare ``err.nbases`` non
    cambia nulla se la stima usava gia' tutte le basi disponibili.
    """
    config = config_ridotta(tmp_path)
    politica = PoliticaRetry.da_config(config)
    decisione = politica.decidi("E-S3-01", 1, raddoppia("err.nbases"), config, "tutte le basi gia' usate")
    assert decisione.motivo is Motivo.AZIONE_INUTILE
    assert politica.decidi("E-S3-01", 1, raddoppia("err.nbases"), config).ritenta


class _S3Doppione(PipelineStep):
    """Un S3 finto che non converge mai, dichiarando o no il retry inutile."""

    passo: ClassVar[Passo] = Passo.S3
    parametri: ClassVar[tuple[str, ...]] = ("err",)
    aggiustamenti: ClassVar = {"E-S3-01": raddoppia("err.nbases")}
    inutile: ClassVar[bool] = False

    def __init__(self) -> None:
        self.nbases: list[float] = []

    def calcola(self, contesto: StepContext) -> Produzione:
        """Annota ``err.nbases`` e fallisce sempre con ``E-S3-01``, con
        ``RITENTARE_INUTILE`` se configurato.
        """
        self.nbases.append(contesto.config.err.nbases)
        extra = {RITENTARE_INUTILE: "tutte le basi gia' usate"} if self.inutile else {}
        raise errore("E-S3-01", "doppione che non converge", **extra)


@pytest.mark.parametrize("inutile", [True, False], ids=["dichiarato-inutile", "utile"])
def test_l_esecutore_ritenta_solo_se_serve(inventario_ridotto, tmp_path, inutile):
    """
    **Obiettivo**: Verificare che l'esecutore non ritenti E-S3-01 quando la
    fase lo dichiara inutile, e lo ritenti con ``err.nbases`` raddoppiato
    quando non lo dichiara.

    **Razionale scientifico e sistemistico**: Un nuovo tentativo identico al
    primo non serve; uno con piu' basi puo' servire.
    """
    from amplicon16s.runner.project import passi_realizzati

    config = config_ridotta(tmp_path)
    esegui_s0(config)
    doppione = type("S3Doppione", (_S3Doppione,), {"inutile": inutile})()
    passi = {Passo.S0: passi_realizzati()[Passo.S0], Passo.S3: doppione}

    class _Fatta(PipelineStep):
        parametri: ClassVar[tuple[str, ...]] = ()

        def calcola(self, contesto):
            a = contesto.albero.scrivi_testo(self.cartella, f"{self.passo}.txt", "x")
            return Produzione((a,))

    for passo in (Passo.S1, Passo.S2):
        passi[passo] = type(f"Fatta{passo}", (_Fatta,), {"passo": passo})()
    run = ProjectRun(config, passi=passi)
    esito = Esecutore(run, fino_a=Passo.S3).esegui()

    assert esito.conclusione is Conclusione.ARRESTATA
    assert esito.punto.codice == "E-S3-01"
    if inutile:
        assert doppione.nbases == [1e8]
        assert esito.punto.tentativi == 1
        assert "tutte le basi gia' usate" in esito.punto.motivo
    else:
        assert doppione.nbases == [1e8, 2e8]
        assert esito.punto.tentativi == 2
        assert "revisione umana" in esito.punto.motivo


# --------------------------------------------------------------------------- #
# S3 vera, sulla versione ridotta                                              #
# --------------------------------------------------------------------------- #

_MOTIVO_ASSENTI = motivo_pacchetti_r_assenti("dada2", "ggplot2", "ShortRead", "jsonlite")


@pytest.fixture
def dada2():
    """Richiede R con dada2: salta senza, ma in CI fallisce."""
    if _MOTIVO_ASSENTI is not None:
        if os.environ.get("AMPLICON16S_RICHIEDI_BIOC") == "1":
            pytest.fail(f"Bioconductor e' richiesto in questo ambiente: {_MOTIVO_ASSENTI}")
        pytest.skip(_MOTIVO_ASSENTI)


def _fino_a_s3(config):
    """Esegue il grafo fino a S3 e restituisce l'esecuzione e il suo esito."""
    run = ProjectRun(config)
    return run, Esecutore(run, fino_a=Passo.S3).esegui()


@pytest.fixture(scope="module")
def stimata(ridotta_calcolata):
    """L'esecuzione di base condivisa, S0-S3, in sola lettura."""
    return ridotta_calcolata


def _con_err(stimata, tmp_path: Path, **err):
    """L'albero della stima di base, copiato, con i soli parametri err cambiati.

    S0, S1 e S2 non dichiarano i parametri err (salvo err.batch_column, letto
    da S0): restano concluse, e l'esecuzione rifa' soltanto S3. Il riferimento
    tassonomico resta quello della stima di base, perche' S0 lo dichiara.
    """
    run, _ = stimata
    dati = run.config.model_dump(mode="python")
    dati["io"]["out_root"] = str(tmp_path / "out")
    dati["err"].update(err)
    shutil.copytree(run.config.io.out_root, tmp_path / "out")
    return _fino_a_s3(valida(dati))


def _modelli(run) -> Path:
    """La cartella dei modelli d'errore dell'esecuzione."""
    return run.albero.cartella(Fase.ERROR_MODELS)


def _impronte(cartella: Path) -> dict[str, str]:
    """L'MD5 di ogni file della cartella, esclusi i file del ponte e i manifesti."""
    return {
        p.name: hashlib.md5(p.read_bytes()).hexdigest()
        for p in sorted(cartella.iterdir())
        if not p.name.startswith(("rbridge_", "manifest"))
    }


def test_s3_produce_due_modelli_e_la_corrispondenza(dada2, stimata):
    """
    **Obiettivo**: Verificare che S3 sulla versione ridotta produca un modello
    e un grafico per ciascuna corsa, la corrispondenza campione-modello e il
    proprio manifesto.

    **Razionale scientifico e sistemistico**: S4 leggera' la corrispondenza
    invece di ricalcolarla.
    """
    run, esito = stimata
    assert esito.conclusione is Conclusione.COMPLETATA
    nomi = set(run.albero.manifesto_passo(Passo.S3, Fase.ERROR_MODELS).nomi)
    for corsa in (CORSA_A, CORSA_B):
        assert {f"modello_{corsa}.rds", f"modello_{corsa}.png"} <= nomi
    assert {"corrispondenza.tsv", "modelli.json", "convergenza.json"} <= nomi

    righe = (_modelli(run) / "corrispondenza.tsv").read_text().splitlines()
    corse = {r["accession"]: r["corsa"] for r in selezione()}
    assert len(righe) == 29
    for riga in righe[1:]:
        campione, corsa, modello = riga.split("\t")
        assert corsa == corse[campione]
        assert modello == f"modello_{corsa}.rds"

    modelli = json.loads((_modelli(run) / "modelli.json").read_text())
    assert all(m["convergenza"] for m in modelli.values())
    # Sulla versione ridotta ogni corsa ha meno basi di err.nbases: si usano tutte.
    assert all(m["tutte_le_basi_usate"] for m in modelli.values())


def test_i_grafici_mostrano_errori_decrescenti_con_la_qualita(dada2, stimata):
    """
    **Obiettivo**: Verificare che i grafici siano PNG validi e che il modello
    che rappresentano abbia, per ogni sostituzione, un tasso d'errore stimato
    piu' basso a qualita' 40 che a qualita' 15.

    **Razionale scientifico e sistemistico**: E' la forma attesa del modello: a
    qualita' piu' alta, meno errori. Un modello piatto o crescente sarebbe il
    segno di una stima sbagliata.
    """
    run, _ = stimata
    for corsa in (CORSA_A, CORSA_B):
        assert (_modelli(run) / f"modello_{corsa}.png").read_bytes()[:8] == b"\x89PNG\r\n\x1a\n"
        verifica = subprocess.run(
            [str(trova_rscript()), "--vanilla", "-e",
             f"e <- readRDS('{_modelli(run) / f'modello_{corsa}.rds'}')$err_out; "
             "s <- !rownames(e) %in% c('A2A','C2C','G2G','T2T'); "
             "quit(status = !all(e[s, '40'] < e[s, '15']))"],
            capture_output=True, check=False,
        )
        assert verifica.returncode == 0, corsa


def test_due_esecuzioni_identiche_danno_artefatti_identici(dada2, stimata, tmp_path):
    """
    **Obiettivo**: Verificare che una seconda esecuzione da zero produca in
    ``04_error_models`` gli stessi byte della prima, modelli e grafici compresi.

    **Razionale scientifico e sistemistico**: La verifica finale di
    riproducibilita' confronta i checksum da un clone pulito; un grafico con
    la data di creazione la renderebbe impossibile.
    """
    run, _ = stimata
    secondo, esito = _fino_a_s3(config_ridotta(tmp_path))
    assert esito.conclusione is Conclusione.COMPLETATA
    assert _impronte(_modelli(secondo)) == _impronte(_modelli(run))


def test_cambiare_err_nbases_rifa_solo_s3(dada2, stimata, tmp_path):
    """
    **Obiettivo**: Verificare, con le fasi vere, che dopo S0-S3 un cambiamento
    di ``err.nbases`` lasci concluse S0, S1 e S2 e che la ripresa rifaccia
    soltanto S3.

    **Razionale scientifico e sistemistico**: Le fasi dichiarano i parametri da
    cui dipendono: S0-S2 non leggono err.nbases, e rifarle costerebbe un quarto
    d'ora sul dataset completo senza cambiare nulla.
    """
    run, _ = stimata
    dati = run.config.model_dump(mode="python")
    dati["io"]["out_root"] = str(tmp_path / "out")
    dati["err"]["nbases"] = 2e8
    shutil.copytree(run.config.io.out_root, tmp_path / "out")
    cambiata = ProjectRun(valida(dati))

    situazione = cambiata.valuta().situazioni
    for passo in (Passo.S0, Passo.S1, Passo.S2):
        assert situazione[passo].stato is StatoPasso.COMPLETATA, passo
    assert situazione[Passo.S3].motivo == "configurazione cambiata"

    esito = Esecutore(cambiata, fino_a=Passo.S3).esegui()
    assert esito.conclusione is Conclusione.COMPLETATA
    assert [r.passo for r in esito.eseguite] == [Passo.S3]


def test_con_err_batch_column_nullo_un_solo_modello(dada2, tmp_path):
    """
    **Obiettivo**: Verificare che con ``err.batch_column`` nullo S3 produca un
    solo modello su tutti i campioni e che l'esecuzione resti valida.

    **Razionale scientifico e sistemistico**: E' il ramo che rende la pipeline
    utilizzabile senza il file di arricchimento del lotto.
    """
    run, esito = _fino_a_s3(config_ridotta(tmp_path, err={"batch_column": None}))
    assert esito.conclusione is Conclusione.COMPLETATA
    nomi = set(run.albero.manifesto_passo(Passo.S3, Fase.ERROR_MODELS).nomi)
    assert {f"modello_{MODELLO_UNICO}.rds", f"modello_{MODELLO_UNICO}.png"} <= nomi
    assert not any(CORSA_A in n or CORSA_B in n for n in nomi)
    righe = (_modelli(run) / "corrispondenza.tsv").read_text().splitlines()[1:]
    assert {r.split("\t")[2] for r in righe} == {f"modello_{MODELLO_UNICO}.rds"}
    assert run.valuta().situazioni[Passo.S3].stato is StatoPasso.COMPLETATA


def test_la_mancata_convergenza_senza_basi_in_piu_non_si_ritenta(dada2, stimata, tmp_path):
    """
    **Obiettivo**: Verificare che con ``err.max_consist = 1`` il modello non
    converga, che S3 sollevi ``E-S3-01``, e che, avendo gia' usato tutte le
    basi disponibili, il retry non venga tentato e il motivo sia dichiarato.

    **Razionale scientifico e sistemistico**: Sulla versione ridotta le basi
    disponibili sono meno di 1e8: raddoppiare ``err.nbases`` darebbe un
    tentativo identico.
    """
    run, esito = _con_err(stimata, tmp_path, max_consist=1)
    assert esito.eseguite == ()  # S0-S2 restano concluse; S3 non si conclude
    assert esito.conclusione is Conclusione.ARRESTATA
    assert esito.punto.passo is Passo.S3
    assert esito.punto.codice == "E-S3-01"
    assert esito.punto.tentativi == 1
    assert "tutte le basi disponibili" in esito.punto.motivo
    assert run.albero.manifesto_passo(Passo.S3, Fase.ERROR_MODELS) is None


def test_la_mancata_convergenza_con_basi_in_piu_si_ritenta(dada2, stimata, tmp_path):
    """
    **Obiettivo**: Verificare che con poche basi per la stima, meno di quelle
    disponibili, la mancata convergenza venga ritentata con ``err.nbases``
    raddoppiato, e si fermi dopo i tentativi ammessi.

    **Razionale scientifico e sistemistico**: Qui il raddoppio cambia davvero i
    dati della stima: il nuovo tentativo ha senso.
    """
    run, esito = _con_err(stimata, tmp_path, max_consist=1, nbases=100000)
    assert esito.conclusione is Conclusione.ARRESTATA
    assert esito.punto.codice == "E-S3-01"
    assert esito.punto.tentativi == 2
    assert "revisione umana" in esito.punto.motivo
    modelli = json.loads((_modelli(run) / "modelli.json").read_text())
    assert {m["nbases"] for m in modelli.values()} == {200000}


def test_l_errore_e_ritentabile_con_revisione():
    """
    **Obiettivo**: Verificare che ``E-S3-01`` sia della categoria retry poi
    revisione umana.

    **Razionale scientifico e sistemistico**: Se non converge nemmeno con piu'
    basi, la causa va capita: il lotto potrebbe raccogliere dati eterogenei.
    """
    assert isinstance(errore("E-S3-01"), ErroreRitentabileConRevisione)


# --------------------------------------------------------------------------- #
# Dataset completo, nel container                                              #
# --------------------------------------------------------------------------- #


@pytest.mark.dati_reali
def test_s3_sul_dataset_completo(dada2, tmp_path):
    """
    **Obiettivo**: Verificare che sul dataset completo S3 produca i due modelli,
    uno per corsa, entrambi convergenti, stimati su una parte delle basi
    disponibili e non su tutte.

    **Razionale scientifico e sistemistico**: E' la situazione in cui
    ``err.nbases`` conta davvero: ogni corsa ha circa due miliardi di basi, e la
    stima ne usa poco piu' di 1e8.
    """
    import logging

    from amplicon16s.config.schema import carica, valida

    percorso = os.environ.get("AMPLICON16S_CONFIG_DATI_REALI")
    if not percorso or not Path(percorso).expanduser().is_file():
        pytest.skip("AMPLICON16S_CONFIG_DATI_REALI non impostata o file assente")
    dati = carica(Path(percorso).expanduser()).model_dump(mode="python")
    dati["io"]["out_root"] = str(tmp_path / "out")
    logging.getLogger("amplicon16s").setLevel(logging.WARNING)

    run, esito = _fino_a_s3(valida(dati))
    assert esito.conclusione is Conclusione.COMPLETATA
    s3 = next(r for r in esito.eseguite if r.passo is Passo.S3)
    print(f"\nS3 sul dataset completo: {s3.secondi} s, {dict(s3.metriche)}")
    modelli = json.loads((_modelli(run) / "modelli.json").read_text())
    assert sorted(modelli) == [CORSA_A, CORSA_B]
    for modello in modelli.values():
        assert modello["campioni"] == 480
        assert modello["convergenza"]
        assert not modello["tutte_le_basi_usate"]
        assert 1e8 < modello["basi_usate"] < modello["basi_disponibili"]
