"""Test del ponte verso R.

La prima parte prova il contratto e la classificazione dell'esito senza
lanciare R. La seconda lancia davvero gli script doppioni di
``tests/r_doppioni/``: il ponte e' il componente da cui dipendono tutte le
fasi di calcolo, e simularne l'esecuzione proverebbe la simulazione.

Senza R quei test si saltano. In CI la variabile ``AMPLICON16S_RICHIEDI_R``
trasforma il salto in un fallimento: un ambiente di integrazione che ha perso
R non deve poter passare per verde.
"""

from __future__ import annotations

import functools
import json
import logging
import os
import signal
import subprocess
import sys
from pathlib import Path

import pytest

from amplicon16s.errors.exceptions import (
    ErrorePipeline,
    ErroreRevisioneUmana,
    ErroreRitentabile,
    ErroreRitentabileConRevisione,
)
from amplicon16s.io_layer.artifacts import AlberoOutput, Fase
from amplicon16s.logging.logger import NOME_FILE_LOG, chiudi, configura
from amplicon16s.rbridge.payload import (
    NOME_ESITO,
    NOME_RICHIESTA,
    PROTOCOLLO,
    Dichiarazione,
    DichiarazioneNonValida,
    Stato,
    leggi_dichiarazione,
    scrivi_richiesta,
)
from amplicon16s.rbridge.runner import (
    VARIABILE_CARTELLA_R,
    VARIABILE_RSCRIPT,
    Condizione,
    cartella_r,
    classifica,
    esegui_script as _esegui_script,
    trova_rscript,
)

DOPPIONI = Path(__file__).resolve().parent / "r_doppioni"
FASE = Fase.FILTERED

#: Limite di memoria virtuale imposto al figlio nei test di esaurimento: basta
#: ad avviare R e jsonlite, non a quanto il doppione cerca di allocare.
LIMITE = 1024**3
#: Circa 1,9 GB di numeri in doppia precisione: oltre il limite.
TROPPI = 250_000_000

#: Ogni lancio nei test ha un tetto di tempo: un processo R bloccato deve far
#: fallire il test, non lasciare appesa la suite.
esegui_script = functools.partial(_esegui_script, tempo_massimo_s=120)


@pytest.fixture(autouse=True)
def uscite_pulite():
    chiudi()
    yield
    chiudi()


@pytest.fixture
def albero(tmp_path) -> AlberoOutput:
    return AlberoOutput(tmp_path / "out")


# --------------------------------------------------------------------------- #
# Disponibilita' di R                                                          #
# --------------------------------------------------------------------------- #


def _motivo_r_assente() -> str | None:
    rscript = trova_rscript()
    if rscript is None:
        return "Rscript non disponibile"
    prova = subprocess.run(
        [str(rscript), "--vanilla", "-e",
         "quit(status = !requireNamespace('jsonlite', quietly = TRUE))"],
        capture_output=True,
        check=False,
    )
    if prova.returncode != 0:
        return "il pacchetto R jsonlite non e' installato"
    return None


_MOTIVO_R_ASSENTE = _motivo_r_assente()


@pytest.fixture
def r():
    """Richiede R: salta senza, ma in CI fallisce."""
    if _MOTIVO_R_ASSENTE is not None:
        if os.environ.get("AMPLICON16S_RICHIEDI_R") == "1":
            pytest.fail(f"R e' richiesto in questo ambiente: {_MOTIVO_R_ASSENTE}")
        pytest.skip(_MOTIVO_R_ASSENTE)


#: Ambiente di una macchina in italiano, come lo vede il processo chiamante.
AMBIENTE_ITALIANO = {"LANGUAGE": "it", "LC_ALL": "it_IT.UTF-8"}

#: Un'allocazione impossibile: fallisce subito, senza consumare memoria, con
#: lo stesso messaggio di un esaurimento vero.
ALLOCAZIONE_IMPOSSIBILE = (
    "cat(tryCatch({numeric(1e15); 'allocato'}, "
    "error = function(e) conditionMessage(e)))"
)


def _motivo_r_non_traduce() -> str | None:
    """Perche' R, in un ambiente italiano e senza forzatura, parlerebbe inglese.

    La sonda non cerca la traduzione di un testo preciso, che cambia fra le
    versioni di R: provoca un'allocazione impossibile e guarda il messaggio
    che R produce davvero. Serve per non far passare il test senza aver
    provato nulla. Con R 4.5 il test si salta sempre: il messaggio e'
    diventato "cannot allocate vector of size %0.1f %s" e i cataloghi di
    traduzione non ne contengono alcuna versione, in nessuna lingua. Il
    difetto si presenta quindi solo con versioni precedenti, come la 4.3.3
    su cui e' stato trovato.
    """
    if _MOTIVO_R_ASSENTE is not None:
        return _MOTIVO_R_ASSENTE
    sonda = subprocess.run(
        [str(trova_rscript()), "--vanilla", "-e", ALLOCAZIONE_IMPOSSIBILE],
        env={**os.environ, **AMBIENTE_ITALIANO},
        capture_output=True,
        text=True,
        encoding="utf-8",
        check=False,
    )
    messaggio = sonda.stdout.strip()
    if sonda.returncode != 0 or not messaggio:
        return f"sonda della lingua fallita: {sonda.stderr.strip()[-200:]}"
    if "cannot allocate" in messaggio:
        return (
            "in un ambiente italiano R riporta l'allocazione fallita in inglese "
            f"anche senza forzatura ({messaggio!r}): il test non proverebbe nulla"
        )
    return None


_MOTIVO_R_NON_TRADUCE = _motivo_r_non_traduce()


@pytest.fixture
def r_che_traduce(r):
    """Richiede un R che traduca il messaggio di allocazione fallita."""
    if _MOTIVO_R_NON_TRADUCE is not None:
        pytest.skip(_MOTIVO_R_NON_TRADUCE)


solo_linux = pytest.mark.skipif(
    not sys.platform.startswith("linux"),
    reason="limite di memoria e numeri dei segnali verificati solo su Linux",
)


def _righe_log(radice: Path) -> list[dict]:
    percorso = radice / Fase.LOGS.value / NOME_FILE_LOG
    return [json.loads(r) for r in percorso.read_text(encoding="utf-8").splitlines()]


# --------------------------------------------------------------------------- #
# Il contratto, senza R                                                        #
# --------------------------------------------------------------------------- #


def _dichiara(percorso: Path, invocazione: str = "abc", **campi) -> None:
    documento = {
        "protocollo": PROTOCOLLO,
        "invocazione": invocazione,
        "stato": "riuscito",
        "codice": None,
        "messaggio": "",
        "artefatti": [],
        **campi,
    }
    percorso.write_text(json.dumps(documento), encoding="utf-8")


def test_la_richiesta_porta_parametri_e_riferimenti(tmp_path):
    richiesta, esito = scrivi_richiesta(
        tmp_path, "abc", {"soglia": 2.5, "file": tmp_path / "x.fastq.gz"}
    )
    documento = json.loads(richiesta.read_text(encoding="utf-8"))

    assert richiesta.name == NOME_RICHIESTA
    assert documento["protocollo"] == PROTOCOLLO
    assert documento["invocazione"] == "abc"
    assert documento["esito"] == str(esito)
    assert documento["cartella_fase"] == str(tmp_path)
    assert documento["parametri"] == {
        "soglia": 2.5, "file": str(tmp_path / "x.fastq.gz")
    }


def test_la_richiesta_rimuove_l_esito_di_un_tentativo_precedente(tmp_path):
    _dichiara(tmp_path / NOME_ESITO)
    _, esito = scrivi_richiesta(tmp_path, "nuova", {})
    assert not esito.exists()


@pytest.mark.parametrize("valore", [float("nan"), float("inf"), object()])
def test_un_parametro_non_json_e_respinto(tmp_path, valore):
    with pytest.raises((ValueError, TypeError)):
        scrivi_richiesta(tmp_path, "abc", {"x": valore})
    assert not (tmp_path / NOME_RICHIESTA).exists()


def test_dichiarazione_assente_e_nessuna_dichiarazione(tmp_path):
    assert leggi_dichiarazione(tmp_path / NOME_ESITO, "abc") is None


def test_dichiarazione_di_un_errore(tmp_path):
    percorso = tmp_path / NOME_ESITO
    _dichiara(percorso, stato="errore_catalogo", codice="E-S2-03", messaggio="lotto 3")
    assert leggi_dichiarazione(percorso, "abc") == Dichiarazione(
        Stato.ERRORE_CATALOGO, "E-S2-03", "lotto 3", ()
    )


@pytest.mark.parametrize(
    "campi",
    [
        {"invocazione": "un'altra"},
        {"protocollo": PROTOCOLLO + 1},
        {"stato": "forse"},
        {"stato": "errore_catalogo", "codice": None},
        {"stato": "riuscito", "codice": "E-S2-03"},
        {"artefatti": ["../fuori.rds"]},
        {"artefatti": ["/assoluto.rds"]},
        {"artefatti": "uno.rds"},
    ],
    ids=lambda c: next(iter(c)),
)
def test_dichiarazione_che_viola_il_contratto(tmp_path, campi):
    percorso = tmp_path / NOME_ESITO
    _dichiara(percorso, **campi)
    with pytest.raises(DichiarazioneNonValida):
        leggi_dichiarazione(percorso, "abc")


def test_dichiarazione_troncata(tmp_path):
    percorso = tmp_path / NOME_ESITO
    percorso.write_text('{"protocollo": 1, "stato": "riusc', encoding="utf-8")
    with pytest.raises(DichiarazioneNonValida, match="illeggibile"):
        leggi_dichiarazione(percorso, "abc")


# --------------------------------------------------------------------------- #
# Classificazione, senza R                                                     #
# --------------------------------------------------------------------------- #


def _d(stato: Stato, codice: str | None = None, messaggio: str = "") -> Dichiarazione:
    return Dichiarazione(stato, codice, messaggio, ())


@pytest.mark.parametrize(
    ("uscita", "dichiarazione", "stderr", "attesa"),
    [
        (0, _d(Stato.RIUSCITO), "", Condizione.RIUSCITO),
        (3, _d(Stato.ERRORE_CATALOGO, "E-S2-03"), "", Condizione.ERRORE_DICHIARATO),
        (3, _d(Stato.ERRORE_NON_CATALOGATO, messaggio="subscript out of bounds"),
         "", Condizione.ERRORE_NON_CATALOGATO),
        (3, _d(Stato.ERRORE_NON_CATALOGATO,
               messaggio="cannot allocate vector of size 1.9 Gb"),
         "", Condizione.MEMORIA_ESAURITA),
        (3, _d(Stato.ERRORE_NON_CATALOGATO, messaggio="std::bad_alloc"),
         "", Condizione.MEMORIA_ESAURITA),
        (1, None, "Error: cannot allocate vector of size 381.5 Mb",
         Condizione.MEMORIA_ESAURITA),
        (-signal.SIGKILL, None, "", Condizione.MEMORIA_ESAURITA),
        (128 + signal.SIGKILL, None, "", Condizione.MEMORIA_ESAURITA),
        (-signal.SIGSEGV, None, "*** caught segfault ***",
         Condizione.NESSUNA_DICHIARAZIONE),
        (1, None, "Error in source: cannot open file",
         Condizione.NESSUNA_DICHIARAZIONE),
        (0, None, "", Condizione.NESSUNA_DICHIARAZIONE),
        # Successo dichiarato, poi il processo e' morto: non e' un successo.
        (-signal.SIGSEGV, _d(Stato.RIUSCITO), "", Condizione.NESSUNA_DICHIARAZIONE),
        (-signal.SIGKILL, _d(Stato.RIUSCITO), "", Condizione.MEMORIA_ESAURITA),
    ],
)
def test_classificazione(uscita, dichiarazione, stderr, attesa):
    assert classifica(uscita, dichiarazione, stderr) is attesa


# --------------------------------------------------------------------------- #
# Interprete e cartella degli script                                           #
# --------------------------------------------------------------------------- #


def test_rscript_indicato_che_non_esiste(tmp_path):
    assert trova_rscript(tmp_path / "Rscript") is None


def test_rscript_dalla_variabile_d_ambiente(tmp_path, monkeypatch):
    finto = tmp_path / "Rscript"
    finto.write_text("#!/bin/sh\n", encoding="utf-8")
    finto.chmod(0o755)
    monkeypatch.setenv(VARIABILE_RSCRIPT, str(finto))
    assert trova_rscript() == finto


def test_la_cartella_r_predefinita_e_quella_del_repository():
    assert (cartella_r() / "lib" / "errors.R").is_file()


def test_la_cartella_r_si_indica(tmp_path, monkeypatch):
    monkeypatch.setenv(VARIABILE_CARTELLA_R, str(tmp_path))
    assert cartella_r() == tmp_path


def test_senza_interprete_l_errore_e_del_catalogo(albero, tmp_path):
    with pytest.raises(ErroreRevisioneUmana) as info:
        esegui_script(
            DOPPIONI / "successo.R", {}, albero, FASE, rscript=tmp_path / "nessuno"
        )
    assert info.value.codice == "E-R-01"
    assert info.value.contesto["condizione_r"] == "non_avviato"


def test_un_codice_di_memoria_inesistente_e_un_difetto_del_chiamante(albero):
    with pytest.raises(KeyError):
        esegui_script(DOPPIONI / "successo.R", {}, albero, FASE, codice_memoria="E-S99-01")


# --------------------------------------------------------------------------- #
# Successo                                                                     #
# --------------------------------------------------------------------------- #


def test_successo(r, albero):
    parametri = {
        "soglia": 2.5,
        "nomi": ["a", "b"],
        "annidati": {"intero": 137, "vero": True},
        "letture": {"C1": 1200, "C2": 0},
    }
    esito = esegui_script(DOPPIONI / "successo.R", parametri, albero, FASE)

    assert esito.riuscito
    assert esito.codice_uscita == 0
    cartella = albero.cartella(FASE)
    assert (cartella / "artefatto.json").is_file()
    assert "uscita standard del doppione" in esito.stdout

    # I parametri attraversano il confine intatti, in andata e ritorno.
    ricevuti = json.loads((cartella / "artefatto.json").read_text(encoding="utf-8"))
    assert ricevuti["ricevuti"] == parametri

    # Gli artefatti dichiarati sono nel manifesto, quelli del contratto no.
    assert [a.nome for a in esito.artefatti] == ["artefatto.json", "letture_doppione.tsv"]
    assert set(albero.manifesto(FASE)) == {"artefatto.json", "letture_doppione.tsv"}
    assert albero.fase_completa(FASE)


def test_tracciamento_delle_letture(r, albero):
    esegui_script(
        DOPPIONI / "successo.R", {"letture": {"C1": 1200, "C2": 0}}, albero, FASE
    )
    righe = (albero.cartella(FASE) / "letture_doppione.tsv").read_text().splitlines()
    assert righe == [
        "campione\tpasso\tletture",
        "C1\tdoppione\t1200",
        "C2\tdoppione\t0",
    ]


def test_un_secondo_tentativo_non_legge_l_esito_del_primo(r, albero):
    esegui_script(DOPPIONI / "successo.R", {}, albero, FASE)
    with pytest.raises(ErrorePipeline) as info:
        esegui_script(
            DOPPIONI / "fallimento_dichiarato.R",
            {"modo": "uscita", "messaggio": "x"},
            albero,
            FASE,
        )
    assert info.value.codice == "E-R-02"


def test_un_artefatto_dichiarato_ma_assente_non_e_un_successo(r, albero):
    with pytest.raises(ErroreRevisioneUmana) as info:
        esegui_script(
            DOPPIONI / "successo.R", {"dichiara_inesistente": True}, albero, FASE
        )
    assert info.value.codice == "E-R-02"
    assert "mai_scritto.rds" in info.value.dettaglio
    assert not albero.manifesto(FASE)


# --------------------------------------------------------------------------- #
# Fallimenti                                                                   #
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize(
    ("codice", "classe"),
    [
        ("E-S2-01", ErroreRevisioneUmana),
        ("E-S2-03", ErroreRitentabile),
        ("E-S3-01", ErroreRitentabileConRevisione),
    ],
)
def test_fallimento_dichiarato(r, albero, codice, classe):
    with pytest.raises(classe) as info:
        esegui_script(
            DOPPIONI / "fallimento_dichiarato.R",
            {"codice": codice, "messaggio": "campione C7 azzerato"},
            albero,
            FASE,
        )
    e = info.value
    assert e.codice == codice
    assert e.dettaglio == "campione C7 azzerato"
    assert e.contesto["condizione_r"] == "errore_dichiarato"
    assert e.contesto["codice_uscita"] == 3


def test_un_codice_dichiarato_fuori_catalogo(r, albero):
    with pytest.raises(ErroreRevisioneUmana) as info:
        esegui_script(
            DOPPIONI / "fallimento_dichiarato.R",
            {"codice": "E-S99-01", "messaggio": "x"},
            albero,
            FASE,
        )
    assert info.value.codice == "E-R-03"
    assert info.value.contesto["codice_dichiarato"] == "E-S99-01"


def test_errore_r_non_catalogato(r, albero):
    with pytest.raises(ErroreRevisioneUmana) as info:
        esegui_script(
            DOPPIONI / "fallimento_dichiarato.R",
            {"modo": "non_catalogato", "messaggio": "indice fuori dai limiti"},
            albero,
            FASE,
        )
    assert info.value.codice == "E-R-03"
    assert info.value.dettaglio == "indice fuori dai limiti"


@pytest.mark.parametrize("modo", ["uscita", "segmentazione"])
def test_processo_morto_senza_dichiarare(r, albero, modo):
    """Distinto da un fallimento dichiarato: altro codice, altra condizione."""
    with pytest.raises(ErroreRevisioneUmana) as info:
        esegui_script(
            DOPPIONI / "fallimento_dichiarato.R",
            {"modo": modo, "messaggio": "x"},
            albero,
            FASE,
        )
    e = info.value
    assert e.codice == "E-R-02"
    assert e.contesto["condizione_r"] == "nessuna_dichiarazione"
    assert not (albero.cartella(FASE) / NOME_ESITO).exists()


def test_processo_bloccato_ucciso_allo_scadere_del_tempo(r, albero):
    """Un figlio che non termina non blocca l'orchestratore."""
    with pytest.raises(ErroreRevisioneUmana) as info:
        _esegui_script(
            DOPPIONI / "fallimento_dichiarato.R",
            {"modo": "stallo", "messaggio": "x"},
            albero,
            FASE,
            tempo_massimo_s=3,
        )
    e = info.value
    assert e.codice == "E-R-02"
    assert e.contesto["condizione_r"] == "tempo_scaduto"
    assert "tempo massimo" in e.dettaglio


@solo_linux
def test_un_guasto_grave_di_r_non_ferma_il_chiamante(r, albero):
    """R muore per SIGSEGV; il processo Python prosegue e puo' rilanciarlo."""
    with pytest.raises(ErrorePipeline) as info:
        esegui_script(
            DOPPIONI / "fallimento_dichiarato.R",
            {"modo": "segmentazione", "messaggio": "x"},
            albero,
            FASE,
        )
    assert info.value.contesto["segnale"] == "SIGSEGV"

    # Il chiamante e' ancora vivo e il ponte ancora utilizzabile.
    assert esegui_script(DOPPIONI / "successo.R", {}, albero, FASE).riuscito


# --------------------------------------------------------------------------- #
# Memoria esaurita                                                             #
# --------------------------------------------------------------------------- #


@solo_linux
def test_memoria_intercettata_da_r(r, albero):
    """L'allocazione fallisce davvero, sotto il limite imposto al figlio."""
    with pytest.raises(ErroreRitentabile) as info:
        esegui_script(
            DOPPIONI / "memoria.R",
            {"modo": "allocazione", "elementi": TROPPI},
            albero,
            FASE,
            codice_memoria="E-S4-02",
            limite_memoria_byte=LIMITE,
        )
    e = info.value
    assert e.codice == "E-S4-02"
    assert e.contesto["condizione_r"] == "memoria_esaurita"
    assert "cannot allocate" in e.dettaglio


@solo_linux
def test_sotto_lo_stesso_limite_un_allocazione_piccola_riesce(r, albero):
    """Controllo: il fallimento dipende dalla dimensione, non dal limite in se'."""
    esito = esegui_script(
        DOPPIONI / "memoria.R",
        {"modo": "allocazione", "elementi": 1_000_000},
        albero,
        FASE,
        codice_memoria="E-S4-02",
        limite_memoria_byte=LIMITE,
    )
    assert esito.riuscito


@solo_linux
def test_memoria_esaurita_prima_di_poter_dichiarare(r, albero):
    """R intercetta l'allocazione, ma resta solo l'uscita di errore."""
    with pytest.raises(ErroreRitentabileConRevisione) as info:
        esegui_script(
            DOPPIONI / "memoria.R",
            {"modo": "fuori_contratto", "elementi": TROPPI},
            albero,
            FASE,
            codice_memoria="E-S5-01",
            limite_memoria_byte=LIMITE,
        )
    assert info.value.codice == "E-S5-01"
    assert info.value.contesto["condizione_r"] == "memoria_esaurita"


@solo_linux
def test_memoria_processo_ucciso_dal_sistema(r, albero):
    """Nessun messaggio, nessuna dichiarazione: solo SIGKILL."""
    with pytest.raises(ErroreRitentabile) as info:
        esegui_script(
            DOPPIONI / "memoria.R",
            {"modo": "terminazione"},
            albero,
            FASE,
            codice_memoria="E-S4-02",
        )
    e = info.value
    assert e.codice == "E-S4-02"
    assert e.contesto["condizione_r"] == "memoria_esaurita"
    assert e.contesto["segnale"] == "SIGKILL"


@solo_linux
def test_memoria_in_una_fase_senza_codice(r, albero):
    """Il ponte riconosce la condizione; senza un codice della fase, E-R-04."""
    with pytest.raises(ErroreRevisioneUmana) as info:
        esegui_script(DOPPIONI / "memoria.R", {"modo": "terminazione"}, albero, FASE)
    assert info.value.codice == "E-R-04"


@solo_linux
def test_il_riconoscimento_non_dipende_dalla_lingua(r_che_traduce, albero, monkeypatch):
    """Con la macchina in italiano R tradurrebbe il messaggio di allocazione.

    La fixture ha gia' verificato che in questo ambiente R, senza forzatura,
    parla italiano: se il test passa, e' merito della forzatura del ponte.
    L'invariante e' provato in ogni ambiente dal test deterministico
    sull'ambiente del figlio; questo ne mostra l'effetto su un R reale, dove
    la versione lo consente.
    """
    for nome, valore in AMBIENTE_ITALIANO.items():
        monkeypatch.setenv(nome, valore)
    with pytest.raises(ErroreRitentabile) as info:
        esegui_script(
            DOPPIONI / "memoria.R",
            {"modo": "allocazione", "elementi": TROPPI},
            albero,
            FASE,
            codice_memoria="E-S4-02",
            limite_memoria_byte=LIMITE,
        )
    assert info.value.codice == "E-S4-02"


@pytest.mark.skipif(os.name != "posix", reason="l'interprete finto e' uno script sh")
@pytest.mark.parametrize(
    "ambiente_chiamante",
    [
        {"LANGUAGE": "it", "LC_ALL": "it_IT.UTF-8"},
        {"LANGUAGE": "de:fr", "LC_MESSAGES": "de_DE.UTF-8"},
        {"LANGUAGE": None, "LC_ALL": None, "LANG": "it_IT.UTF-8"},
    ],
    ids=["italiano", "tedesco-francese", "senza-LANGUAGE"],
)
def test_il_figlio_riceve_sempre_language_en(
    albero, tmp_path, monkeypatch, ambiente_chiamante
):
    """L'invariante che protegge dal difetto, verificato senza R.

    Se il figlio riceve LANGUAGE=en, i messaggi di R escono in inglese con
    qualunque versione e il riconoscimento della memoria esaurita non dipende
    dalla lingua della macchina. Un interprete finto registra l'ambiente con
    cui il ponte lo lancia davvero; non dichiara nulla, quindi il ponte
    solleva E-R-02, che qui non interessa.
    """
    for nome, valore in ambiente_chiamante.items():
        if valore is None:
            monkeypatch.delenv(nome, raising=False)
        else:
            monkeypatch.setenv(nome, valore)

    registro = tmp_path / "ambiente_del_figlio.txt"
    finto = tmp_path / "Rscript"
    finto.write_text(f"#!/bin/sh\nenv > '{registro}'\n", encoding="utf-8")
    finto.chmod(0o755)

    with pytest.raises(ErrorePipeline):
        esegui_script(DOPPIONI / "successo.R", {}, albero, FASE, rscript=finto)

    ambiente = dict(
        riga.split("=", 1) for riga in registro.read_text().splitlines() if "=" in riga
    )
    assert ambiente["LANGUAGE"] == "en"


# --------------------------------------------------------------------------- #
# Log strutturato                                                              #
# --------------------------------------------------------------------------- #


def test_l_uscita_di_errore_finisce_nel_log_strutturato(r, tmp_path):
    radice = tmp_path / "out"
    configura(radice, livello_console=logging.CRITICAL)
    with pytest.raises(ErrorePipeline):
        esegui_script(
            DOPPIONI / "fallimento_dichiarato.R",
            {"codice": "E-S2-01", "messaggio": "traccia riconoscibile 7Q"},
            AlberoOutput(radice),
            FASE,
        )
    righe = _righe_log(radice)

    stderr = [r for r in righe if r.get("flusso") == "stderr"]
    assert len(stderr) == 1
    assert "traccia riconoscibile 7Q" in stderr[0]["testo"]
    assert stderr[0]["livello"] == "WARNING"
    assert stderr[0]["script"] == "fallimento_dichiarato.R"

    conclusione = next(r for r in righe if r["messaggio"] == "processo R concluso")
    assert conclusione["condizione_r"] == "errore_dichiarato"
    assert conclusione["codice_uscita"] == 3


def test_anche_un_successo_lascia_le_sue_uscite_nel_log(r, tmp_path):
    radice = tmp_path / "out"
    configura(radice, livello_console=logging.CRITICAL)
    esegui_script(
        DOPPIONI / "successo.R", {"avviso": "messaggio di servizio"},
        AlberoOutput(radice), FASE,
    )
    righe = _righe_log(radice)
    flussi = {r["flusso"]: r["testo"] for r in righe if "flusso" in r}
    assert "uscita standard del doppione" in flussi["stdout"]
    assert "messaggio di servizio" in flussi["stderr"]
