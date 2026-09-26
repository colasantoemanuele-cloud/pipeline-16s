"""Suite di verifica del ponte di esecuzione subprocess verso R/Rscript (``rbridge``).

Inquadramento nel Piano Operativo
---------------------------------
* **Settimana di riferimento**: **Settimana 5 / Settimana 8 (W5/W8 — Fase F1/F3:
  Ponte di comunicazione ed esecuzione subprocess verso R/Rscript)**.
* **Moduli sorgente coperti**:
  - ``src/amplicon16s/rbridge/payload.py``
  - ``src/amplicon16s/rbridge/runner.py``
  - ``R/lib/io_json.R``
  - ``R/lib/errors.R``

Scopo sperimentale e razionale scientifico/sistemistico
-------------------------------------------------------
Il ponte ``rbridge`` è il meccanismo attraverso cui l'orchestratore Python
invoca i pacchetti R/Bioconductor (``dada2``, ``DECIPHER``, ``phangorn``,
``phyloseq``, ``decontam``). I test provano sia il contratto formale in memoria
sia l'esecuzione reale degli script in ``tests/r_doppioni/`` su un vero
interprete ``Rscript``, presidiando cinque proprietà critiche:

1. **Isolamento dei processi (``Rscript --vanilla`` in nuova sessione POSIX)**:
   R non viene mai caricato *in-process* (nessun ``rpy2`` o memoria condivisa C).
   In questo modo un crash catastrofico nelle estensioni C++/Rcpp di DADA2 o
   DECIPHER (``SIGSEGV``, ``std::bad_alloc``, ``malloc`` failure) o un ``SIGKILL``
   dell'OOM killer Linux uccide soltanto il figlio ``Rscript``, lasciando
   intatto l'orchestratore Python per registrare il guasto o applicare il retry.
2. **Contratto atomico su filesystem con UUID di invocazione**: lo scambio avviene
   tramite ``rbridge_richiesta.json`` e ``rbridge_esito.json`` legati da un
   identificativo univoco ``invocazione``, impedendo che un file ``rbridge_esito.json``
   residuo di un tentativo precedente venga scambiato per il risultato attuale.
3. **Riconoscimento multi-segnale della memoria esaurita (OOM)**: intercetta sia
   l'errore R ``cannot allocate vector of size ...`` / ``std::bad_alloc`` (con o
   senza dichiarazione JSON), sia l'uccisione diretta da parte del kernel Linux
   tramite ``SIGKILL`` (codice ``-9`` o ``137``), mappandoli sul ``codice_memoria``
   della fase (es. ``E-S2-03``, ``E-S4-02``, ``E-S5-01``) o su ``E-R-04``.
4. **Normalizzazione linguistica (``LANGUAGE=en``)**: forza l'ambiente del
   sottoprocesso R a emettere diagnostica in inglese anche su sistemi operativi
   configurati in italiano/tedesco/francese, evitando che la localizzazione dei
   messaggi di errore rompa il riconoscimento regex dell'esaurimento di memoria.
5. **Governo dei timeout con ``os.killpg``**: uccide l'intero process group del
   figlio allo scadere di ``tempo_massimo_s``, impedendo che uno script R in
   stallo blocchi indefinitamente la pipeline.
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
    """
    **Obiettivo**: Verificare che ``scrivi_richiesta`` produca ``rbridge_richiesta.json``
    con versione del ``protocollo``, ID di ``invocazione``, percorsi e parametri
    (convertendo ``Path`` in stringa).

    **Razionale Scientifico/Sistemistico**: Garantisce che lo script R riceva
    tutti i parametri e i riferimenti alle directory in un unico documento JSON
    auto-contenuto e tipizzato.
    """
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
    """
    **Obiettivo**: Verificare che ``scrivi_richiesta`` cancelli preventivamente
    un eventuale ``rbridge_esito.json`` già presente nella cartella della fase.

    **Razionale Scientifico/Sistemistico**: Impedisce che, qualora il nuovo
    processo R muoia prima di poter scrivere il proprio esito, il runner legga
    per errore l'esito residuo dell'esecuzione precedente.
    """
    _dichiara(tmp_path / NOME_ESITO)
    _, esito = scrivi_richiesta(tmp_path, "nuova", {})
    assert not esito.exists()


@pytest.mark.parametrize("valore", [float("nan"), float("inf"), object()])
def test_un_parametro_non_json_e_respinto(tmp_path, valore):
    """
    **Obiettivo**: Verificare che ``scrivi_richiesta`` rifiuti con ``ValueError``
    o ``TypeError`` valori ``NaN``, ``Inf`` o oggetti Python arbitrari senza
    lasciare ``rbridge_richiesta.json`` su disco.

    **Razionale Scientifico/Sistemistico**: Lo standard JSON (RFC 8259) e
    ``jsonlite`` in R non ammettono letterali ``NaN``/``Infinity`` non quotati;
    bloccarli sul lato Python evita errori di parsing opachi dentro R.
    """
    with pytest.raises((ValueError, TypeError)):
        scrivi_richiesta(tmp_path, "abc", {"x": valore})
    assert not (tmp_path / NOME_RICHIESTA).exists()


def test_dichiarazione_assente_e_nessuna_dichiarazione(tmp_path):
    """
    **Obiettivo**: Verificare che ``leggi_dichiarazione`` restituisca ``None``
    quando il file ``rbridge_esito.json`` non esiste su disco.

    **Razionale Scientifico/Sistemistico**: Consente al classificatore di
    distinguere un processo R che ha chiuso ordinatamente il contratto da uno
    terminato prematuramente (crash, segfault, ``quit()``).
    """
    assert leggi_dichiarazione(tmp_path / NOME_ESITO, "abc") is None


def test_dichiarazione_di_un_errore(tmp_path):
    """
    **Obiettivo**: Verificare che ``leggi_dichiarazione`` decodifichi
    correttamente uno stato ``errore_catalogo`` con codice ``E-S2-03`` e
    messaggio di dettaglio.

    **Razionale Scientifico/Sistemistico**: Permette alle funzioni R di
    segnalare condizioni di errore note al catalogo Python preservando il codice
    esatto e il dettaglio contestuale (es. quale lotto ha fallito).
    """
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
    """
    **Obiettivo**: Verificare che ``leggi_dichiarazione`` sollevi
    ``DichiarazioneNonValida`` se ``invocazione`` o ``protocollo`` non
    coincidono, se ``stato``/``codice`` sono incoerenti o se gli ``artefatti``
    contengono path assoluti o *path traversal* (``../``).

    **Razionale Scientifico/Sistemistico**: Blindatura del confine tra R e
    Python: rifiuta risposte appartenenti ad altre invocazioni (anti-replay) e
    impedisce a uno script R di dichiarare artefatti fuori dalla cartella della
    fase corrente.
    """
    percorso = tmp_path / NOME_ESITO
    _dichiara(percorso, **campi)
    with pytest.raises(DichiarazioneNonValida):
        leggi_dichiarazione(percorso, "abc")


def test_dichiarazione_troncata(tmp_path):
    """
    **Obiettivo**: Verificare che un file ``rbridge_esito.json`` con JSON
    troncato a metà sollevi ``DichiarazioneNonValida("illeggibile")``.

    **Razionale Scientifico/Sistemistico**: Intercetta scritture interrotte da
    disco pieno o uccisione del processo durante il flush di ``jsonlite``.
    """
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
    """
    **Obiettivo**: Verificare la tabella decisionale di ``classifica(uscita,
    dichiarazione, stderr)`` su tutte le combinazioni di exit code, segnali
    POSIX (``SIGKILL``, ``SIGSEGV``), dichiarazioni JSON e diagnostica ``stderr``.

    **Razionale Scientifico/Sistemistico**: Un processo R può esaurire la RAM in
    tre modi diversi (errore intercettato da ``tryCatch`` in R, ``std::bad_alloc``
    in C++ su ``stderr``, oppure ``SIGKILL`` del kernel Linux) e può persino
    morire di ``SIGSEGV`` *dopo* aver scritto ``stato="riuscito"`` durante la
    deallocazione finale. ``classifica`` unifica tutti i casi di OOM in
    ``Condizione.MEMORIA_ESAURITA`` e nega il successo se il codice di uscita
    non è zero.
    """
    assert classifica(uscita, dichiarazione, stderr) is attesa


# --------------------------------------------------------------------------- #
# Interprete e cartella degli script                                           #
# --------------------------------------------------------------------------- #


def test_rscript_indicato_che_non_esiste(tmp_path):
    """
    **Obiettivo**: Verificare che ``trova_rscript`` restituisca ``None`` se il
    percorso esplicito fornito non esiste.

    **Razionale Scientifico/Sistemistico**: Consente al chiamante di rilevare
    immediatamente un eseguibile ``Rscript`` mancante prima di tentare ``Popen``.
    """
    assert trova_rscript(tmp_path / "Rscript") is None


def test_rscript_dalla_variabile_d_ambiente(tmp_path, monkeypatch):
    """
    **Obiettivo**: Verificare che ``trova_rscript`` rispetti la variabile
    d'ambiente ``AMPLICON16S_RSCRIPT`` quando punta a un file eseguibile.

    **Razionale Scientifico/Sistemistico**: Permette di selezionare una
    specifica installazione di R (es. dentro un modulo HPC o un ambiente Conda/renv
    dedicato) senza modificare il ``PATH`` di sistema.
    """
    finto = tmp_path / "Rscript"
    finto.write_text("#!/bin/sh\n", encoding="utf-8")
    finto.chmod(0o755)
    monkeypatch.setenv(VARIABILE_RSCRIPT, str(finto))
    assert trova_rscript() == finto


def test_la_cartella_r_predefinita_e_quella_del_repository():
    """
    **Obiettivo**: Verificare che ``cartella_r()`` individui la directory ``R/``
    del repository contenente ``lib/errors.R``.

    **Razionale Scientifico/Sistemistico**: Garantisce che gli script R trovino
    sempre le librerie condivise ``io_json.R`` ed ``errors.R``.
    """
    assert (cartella_r() / "lib" / "errors.R").is_file()


def test_la_cartella_r_si_indica(tmp_path, monkeypatch):
    """
    **Obiettivo**: Verificare che la variabile d'ambiente ``AMPLICON16S_R_DIR``
    sovrascriva il percorso restituito da ``cartella_r()``.

    **Razionale Scientifico/Sistemistico**: Consente il disaccoppiamento del
    percorso degli script R quando il pacchetto Python è installato in
    ``site-packages`` dentro un container.
    """
    monkeypatch.setenv(VARIABILE_CARTELLA_R, str(tmp_path))
    assert cartella_r() == tmp_path


def test_senza_interprete_l_errore_e_del_catalogo(albero, tmp_path):
    """
    **Obiettivo**: Verificare che se ``rscript`` non esiste, ``esegui_script``
    sollevi ``ErroreRevisioneUmana`` con codice ``E-R-01`` e
    ``condizione_r == "non_avviato"``.

    **Razionale Scientifico/Sistemistico**: Traduce l'assenza dell'interprete R
    in un errore strutturato del catalogo (``E-R-01``) con istruzioni chiare per
    l'utente anziché in un ``FileNotFoundError`` di sistema.
    """
    with pytest.raises(ErroreRevisioneUmana) as info:
        esegui_script(
            DOPPIONI / "successo.R", {}, albero, FASE, rscript=tmp_path / "nessuno"
        )
    assert info.value.codice == "E-R-01"
    assert info.value.contesto["condizione_r"] == "non_avviato"


def test_un_codice_di_memoria_inesistente_e_un_difetto_del_chiamante(albero):
    """
    **Obiettivo**: Verificare che passare a ``esegui_script`` un
    ``codice_memoria="E-S99-01"`` non presente nel catalogo sollevi subito
    ``KeyError`` prima ancora di avviare R.

    **Razionale Scientifico/Sistemistico**: Valida preventivamente il codice di
    fallback per OOM fornito dallo step Python, evitando che un codice errato
    resti latente finché non si verifica davvero un esaurimento di memoria.
    """
    with pytest.raises(KeyError):
        esegui_script(DOPPIONI / "successo.R", {}, albero, FASE, codice_memoria="E-S99-01")


# --------------------------------------------------------------------------- #
# Successo                                                                     #
# --------------------------------------------------------------------------- #


def test_successo(r, albero):
    """
    **Obiettivo**: Verificare l'esecuzione reale di ``successo.R`` tramite
    ``Rscript``: passaggio intatto dei parametri annidati, scrittura degli
    artefatti e registrazione nel manifesto della fase (escludendo i file di
    servizio del ponte).

    **Razionale Scientifico/Sistemistico**: Dimostra su un vero processo R che
    numeri, liste e dizionari attraversano il confine Python → JSON → R → JSON → Python
    senza perdita di tipo o precisione e che gli artefatti prodotti da R vengono
    immediatamente sigillati con SHA-256 in ``AlberoOutput``.
    """
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
    """
    **Obiettivo**: Verificare che la funzione helper R ``scrivi_tracciamento_letture``
    produca il file TSV ``letture_doppione.tsv`` con le colonne ``campione``,
    ``passo`` e ``letture``.

    **Razionale Scientifico/Sistemistico**: Standardizza il tracciamento della
    perdita di letture campione per campione attraverso i passaggi R (filtraggio,
    denoising, fusione, rimozione chimere) per il report QC finale.
    """
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
    """
    **Obiettivo**: Verificare che se dopo un'esecuzione riuscita di ``successo.R``
    viene lanciato ``fallimento_dichiarato.R`` (che esce senza scrivere esito),
    ``esegui_script`` sollevi ``E-R-02`` e non rilegga il successo precedente.

    **Razionale Scientifico/Sistemistico**: Impedisce il bug critico di *stale
    result*, in cui il crash silenzioso di un ricalcolo verrebbe mascherato dal
    file ``rbridge_esito.json`` lasciato da una corsa precedente.
    """
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
    """
    **Obiettivo**: Verificare che se lo script R dichiara ``stato="riuscito"``
    elencando un artefatto ``mai_scritto.rds`` che non esiste su disco, il ponte
    sollevi ``ErroreRevisioneUmana`` con codice ``E-R-02`` senza aggiornare il manifesto.

    **Razionale Scientifico/Sistemistico**: Non si fida ciecamente della
    dichiarazione verbale dello script R ma verifica l'esistenza fisica su
    filesystem di ogni singolo file promesso prima di sigillare il manifesto.
    """
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
    """
    **Obiettivo**: Verificare che quando lo script R invoca ``ferma_con_codice(codice, ...)``
    il ponte sollevi in Python esattamente la sottoclasse di ``ErrorePipeline``
    prevista dalla categoria del codice (``ErroreRevisioneUmana``, ``ErroreRitentabile``
    o ``ErroreRitentabileConRevisione``).

    **Razionale Scientifico/Sistemistico**: Consente agli script R di partecipare
    direttamente alla macchina a stati di gestione degli errori e di attivare il
    retry automatico con aggiustamento dei parametri sul lato Python.
    """
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
    """
    **Obiettivo**: Verificare che se uno script R dichiara un codice ``E-S99-01``
    non esistente nel catalogo Python, il ponte sollevi ``ErroreRevisioneUmana``
    con codice ``E-R-03``.

    **Razionale Scientifico/Sistemistico**: Intercetta eventuali disallineamenti
    tra i codici usati negli script R e il catalogo centrale Python, evitando
    ``KeyError`` non gestiti.
    """
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
    """
    **Obiettivo**: Verificare che un ``stop("indice fuori dai limiti")`` generico
    in R venga intercettato dal gestore globale di ``errors.R`` e tradotto in
    ``ErroreRevisioneUmana`` con codice ``E-R-03``.

    **Razionale Scientifico/Sistemistico**: Cattura qualsiasi eccezione R
    imprevista (es. bug in un pacchetto Bioconductor o dato malformato)
    riportando il messaggio originale di R dentro il contesto strutturato ``E-R-03``.
    """
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
    """
    **Obiettivo**: Verificare che se R termina con ``q(status=1)`` (``uscita``)
    o muore per ``SIGSEGV`` (``segmentazione``) senza scrivere ``rbridge_esito.json``,
    il ponte sollevi ``ErroreRevisioneUmana`` con codice ``E-R-02`` e
    ``condizione_r == "nessuna_dichiarazione"``.

    **Razionale Scientifico/Sistemistico**: Distingue un errore applicativo R
    catturato da ``tryCatch`` (``E-R-03``) da un aborto brusco del processo o
    della libreria C sottostante (``E-R-02``).
    """
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
    """
    **Obiettivo**: Verificare che uno script R in ciclo infinito (``modo="stallo"``)
    venga terminato tramite ``os.killpg`` allo scadere di ``tempo_massimo_s=3``
    sollevando ``E-R-02`` con ``condizione_r == "tempo_scaduto"``.

    **Razionale Scientifico/Sistemistico**: Impedisce che un deadlock nei thread
    C++ o nell'I/O di R lasci appesa indefinitamente la pipeline o il job HPC,
    uccidendo l'intero gruppo di processi del figlio.
    """
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
    """
    **Obiettivo**: Verificare che dopo l'uccisione di un processo R per
    ``SIGSEGV`` il processo Python chiamante resti vivo e possa immediatamente
    eseguire con successo un nuovo script ``successo.R``.

    **Razionale Scientifico/Sistemistico**: È la dimostrazione sperimentale del
    perché l'architettura a sottoprocessi separati (``Rscript``) è superiore a
    ``rpy2``: un *segmentation fault* in codice C/Fortran di R non abbatte
    l'interprete Python.
    """
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
    """
    **Obiettivo**: Verificare che imponendo ``limite_memoria_byte=1 GiB``
    (``RLIMIT_AS``) e tentando di allocare 250 milioni di ``numeric`` (~1,9 GB),
    il ponte intercetti l'OOM di R e sollevi ``ErroreRitentabile("E-S4-02")`` con
    ``condizione_r == "memoria_esaurita"``.

    **Razionale Scientifico/Sistemistico**: Consente a ``PoliticaRetry`` di
    catturare il fallimento di allocazione di ``dada()`` in S4 su piastre ad
    alta profondità e di rilanciare automaticamente il passo dimezzando
    ``run.batch_size``.
    """
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
    """
    **Obiettivo**: Verificare che sotto il medesimo tetto ``limite_memoria_byte=1 GiB``
    un'allocazione di 1 milione di elementi (~8 MB) completi con ``esito.riuscito is True``.

    **Razionale Scientifico/Sistemistico**: Funge da controllo sperimentale per
    il test precedente, dimostrando che il fallimento era effettivamente causato
    dalla dimensione del vettore allocato e non dall'impossibilità di avviare R
    entro 1 GiB di spazio d'indirizzamento.
    """
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
    """
    **Obiettivo**: Verificare che se R esaurisce la memoria fuori dal blocco
    ``con_contratto`` (così che non può nemmeno scrivere ``rbridge_esito.json``),
    il ponte riconosca comunque ``cannot allocate`` su ``stderr`` e sollevi
    ``ErroreRitentabileConRevisione("E-S5-01")``.

    **Razionale Scientifico/Sistemistico**: Quando la RAM è completamente
    saturata, persino la serializzazione JSON dell'errore dentro R può fallire
    per mancanza di memoria; l'ispezione di ``stderr`` garantisce che l'OOM
    venga riconosciuto anche in assenza di dichiarazione JSON.
    """
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
    """
    **Obiettivo**: Verificare che un processo R terminato da ``SIGKILL``
    (senza alcun messaggio su ``stderr`` né dichiarazione JSON) venga classificato
    come ``memoria_esaurita`` e sollevi ``ErroreRitentabile("E-S4-02")``.

    **Razionale Scientifico/Sistemistico**: Su Linux l'OOM killer del kernel
    invia ``SIGKILL`` (9) al processo che eccede la RAM fisica/cgroup senza
    dargli modo di eseguire handler o stampare messaggi; mappare ``SIGKILL`` su
    ``MEMORIA_ESAURITA`` abilita il retry con riduzione del ``batch_size``.
    """
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
    """
    **Obiettivo**: Verificare che se si verifica un esaurimento di memoria in
    una fase che non ha specificato un ``codice_memoria`` ritentabile, il ponte
    sollevi ``ErroreRevisioneUmana("E-R-04")``.

    **Razionale Scientifico/Sistemistico**: Nelle fasi in cui non esiste un
    parametro di partizionamento riducibile automaticamente, l'OOM viene
    comunque diagnosticato con precisione (``E-R-04``) ma richiede la revisione
    umana.
    """
    with pytest.raises(ErroreRevisioneUmana) as info:
        esegui_script(DOPPIONI / "memoria.R", {"modo": "terminazione"}, albero, FASE)
    assert info.value.codice == "E-R-04"


@solo_linux
def test_il_riconoscimento_non_dipende_dalla_lingua(r_che_traduce, albero, monkeypatch):
    """
    **Obiettivo**: Verificare che anche impostando ``LANGUAGE=it`` e
    ``LC_ALL=it_IT.UTF-8`` nell'ambiente del processo Python chiamante, un R che
    normalmente tradurrebbe l'errore in italiano emetta invece il messaggio in
    inglese permettendo il riconoscimento di ``E-S4-02``.

    **Razionale Scientifico/Sistemistico**: Su workstation Linux localizzate in
    italiano (es. ``"impossibile allocare un vettore di dimensione..."``), le
    espressioni regolari basate su ``"cannot allocate vector"`` fallirebbero se
    il ponte non forzasse ``LANGUAGE=en`` nell'ambiente del sottoprocesso R.
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
    """
    **Obiettivo**: Verificare tramite uno script POSIX spia che, qualunque sia
    la combinazione di ``LANGUAGE``, ``LC_ALL``, ``LC_MESSAGES`` e ``LANG`` nel
    chiamante, il sottoprocesso figlio riceva sempre ``LANGUAGE=en``.

    **Razionale Scientifico/Sistemistico**: Prova l'invariante di normalizzazione
    della lingua in modo deterministico e indipendente dalla versione di R
    installata sulla macchina di test.
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
    """
    **Obiettivo**: Verificare che lo ``stderr`` di uno script R fallito venga
    registrato in ``99_logs/pipeline.jsonl`` con ``flusso="stderr"`` e ``livello="WARNING"``
    insieme all'evento di conclusione del processo R.

    **Razionale Scientifico/Sistemistico**: Conserva nel log JSONL tutti i
    messaggi diagnostici emessi da R/Bioconductor su ``stderr`` prima dell'uscita.
    """
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
    """
    **Obiettivo**: Verificare che anche quando lo script R termina con successo
    i flussi ``stdout`` e ``stderr`` vengano integralmente registrati in ``pipeline.jsonl``.

    **Razionale Scientifico/Sistemistico**: Molte funzioni di ``dada2`` e
    ``phyloseq`` emettono su ``stdout``/``stderr`` statistiche di convergenza e
    avvisi non bloccanti durante corse riuscite; catturarli nel log JSONL
    garantisce la completa ispezionabilità a posteriori.
    """
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
