"""Il ponte verso R: l'unico punto del sistema che esegue codice R.

Ogni script R gira come **processo separato**, mai come libreria caricata nel
processo Python. Un guasto grave dentro una routine di calcolo (un errore di
segmentazione in codice compilato, l'intervento del sistema operativo che
uccide il processo per memoria) abbatterebbe l'intero processo, orchestratore
compreso. Con un processo separato muore il figlio: il padre legge l'esito e
decide cosa fare.

Il contratto passa per il filesystem (:mod:`amplicon16s.rbridge.payload`), e
il ponte traduce ciò che trova in un'eccezione della gerarchia di
:class:`~amplicon16s.errors.exceptions.ErrorePipeline`, così un errore R si
cattura e si registra come qualunque altro errore della pipeline. I casi sono
sei:

* **riuscito**: lo script ha dichiarato il successo, è uscito con 0 e gli
  artefatti dichiarati esistono; vengono registrati nel manifesto della fase;
* **errore dichiarato**: lo script ha dichiarato un codice del catalogo, e
  l'eccezione porta esattamente quel codice;
* **memoria esaurita**: riconosciuta dal ponte, ma il codice lo sceglie la
  fase con ``codice_memoria``, perché lo stesso guasto ha codici diversi in
  fasi diverse (``E-S1-06`` nel profilo, ``E-S4-02`` durante l'inferenza,
  ``E-S5-01`` sulla tabella delle varianti). Una fase che non ne indica uno riceve ``E-R-04``;
* **errore non catalogato**: un errore R che lo script non ha ricondotto a un
  codice, ``E-R-03``;
* **nessuna dichiarazione**: il processo è morto senza poter scrivere
  l'esito, ``E-R-02``;
* **tempo scaduto**: il processo ha superato ``run.r_timeout_s`` ed è stato
  ucciso con il suo gruppo, ``E-R-05``, perché un figlio bloccato non deve
  poter bloccare l'orchestratore. Senza il parametro non c'è limite.

La memoria esaurita si manifesta in tre modi, e il ponte li riconosce
tutti. Se il sistema uccide un processo figlio del calcolo parallelo di R, il
genitore sopravvive con un risultato mancante e un avviso
(``did not deliver results``): l'errore che ne segue, qualunque sia, e' letto
come memoria esaurita. Se l'interprete intercetta l'allocazione fallita, l'errore arriva
con il messaggio di R (dichiarato come errore non catalogato, oppure solo
sull'uscita di errore se il processo non è riuscito a dichiarare). Se invece è
il sistema operativo a uccidere il processo, non resta alcun messaggio, solo
il segnale ``SIGKILL``: nessuna dichiarazione e ``SIGKILL`` sono letti come
memoria esaurita. È una presunzione (anche un operatore può inviare quel
segnale), ma è quella giusta: il solo retry ammesso per la memoria riduce la
dimensione del lotto, e non cambia alcuna assunzione metodologica.

I messaggi di R sono tradotti nella lingua della macchina, quindi il ponte
impone al figlio ``LANGUAGE=en``: il riconoscimento della memoria esaurita
non può dipendere dalla localizzazione. Il difetto dipende dalla versione di
R. Con la 4.3.3 in italiano l'allocazione fallita diventa «non è possibile
allocare un vettore di dimensione…». In R 4.5 il messaggio è cambiato in
``cannot allocate vector of size %0.1f %s`` e i cataloghi di traduzione non
ne hanno alcuna versione, in nessuna lingua: con la 4.5.2 di riferimento
esce in inglese comunque. La forzatura protegge chi esegue la pipeline con
versioni precedenti.

Il ponte impone anche ``LC_ALL=C.UTF-8``: ordinamenti e confronti di stringhe
negli script R non dipendono cosi' dalla macchina, e la codifica nativa e'
sempre UTF-8. Gli script ordinano comunque con ``method = "radix"``, che non
dipende dalla lingua, e serializzano con ``salva_rds`` (``R/lib/io_json.R``).
"""

from __future__ import annotations

import logging
import os
import re
import shutil
import signal
import subprocess
import time
import uuid
from collections.abc import Mapping
from dataclasses import dataclass
from enum import StrEnum
from pathlib import Path
from typing import Any, Final

from amplicon16s.errors.catalog import voce
from amplicon16s.errors.exceptions import ErrorePipeline, errore
from amplicon16s.io_layer.artifacts import AlberoOutput, Artefatto, Fase
from amplicon16s.logging.logger import ottieni
from amplicon16s.rbridge.payload import (
    Dichiarazione,
    DichiarazioneNonValida,
    Stato,
    leggi_dichiarazione,
    scrivi_richiesta,
)

__all__ = [
    "VARIABILE_CARTELLA_R",
    "VARIABILE_LIB_R",
    "VARIABILE_RSCRIPT",
    "Condizione",
    "EsitoR",
    "cartella_r",
    "classifica",
    "esegui_script",
    "trova_rscript",
]

#: Variabile d'ambiente con cui indicare l'interprete, se non e' nel PATH.
VARIABILE_RSCRIPT: Final = "AMPLICON16S_RSCRIPT"

#: Variabile d'ambiente con cui indicare la cartella degli script R. Serve
#: quando il pacchetto e' installato in modo non modificabile e non sta piu'
#: accanto al repository, come nel container.
VARIABILE_CARTELLA_R: Final = "AMPLICON16S_R_DIR"

#: Variabile che il ponte imposta nel figlio: la cartella delle funzioni R
#: condivise, da cui ogni script carica il contratto.
VARIABILE_LIB_R: Final = "AMPLICON16S_R_LIB"

#: Quanto delle due uscite del processo finisce nel log. Si conserva la coda,
#: dove sta l'errore; il resto e' progresso dei calcoli.
MAX_CARATTERI_USCITA: Final = 256 * 1024

#: Messaggi con cui un'allocazione fallita arriva dall'interprete: quelli di R,
#: quello di R su macOS, quello del codice C++ dei pacchetti di calcolo,
#: quello di R_Calloc e R_Realloc, con cui alloca il codice C di dada2
#: (``'R_Calloc' could not allocate memory (1000000 of 16 bytes)``, osservato
#: in S4; ``'Calloc'`` prima di R 4.3), e quello del caricatore dinamico
#: quando non c'e' spazio per mappare la libreria di un pacchetto.
_MEMORIA_ESAURITA: Final = re.compile(
    r"cannot allocate (?:vector of size|memory block)"
    r"|vector memory (?:exhausted|limit)"
    r"|memory exhausted"
    r"|std::bad_alloc"
    r"|'(?:R_)?Calloc' could not allocate memory"
    r"|'(?:R_)?Realloc' could not re-allocate memory"
    r"|failed to map segment from shared object"
    # Un figlio del calcolo parallelo che non restituisce il risultato: lo
    # dichiarano gli script (verifica_figli, R/lib/errors.R) e, per il codice
    # parallelo dei pacchetti, l'avviso di parallel::mclapply.
    r"|processi paralleli su \d+ non hanno restituito un risultato",
    re.IGNORECASE,
)

#: L'avviso con cui ``parallel::mclapply`` dice che un processo figlio e' morto
#: senza restituire il risultato: il sistema lo ha ucciso, di norma per memoria
#: esaurita, e l'errore che segue nel genitore e' una conseguenza qualunque.
_FIGLI_SENZA_RISULTATO: Final = re.compile(
    r"did not deliver (?:a )?results?", re.IGNORECASE
)

#: Codici di uscita con cui si presenta un processo ucciso con SIGKILL: il
#: negativo del segnale se lanciato direttamente, 128 + segnale se passato per
#: una shell.
_UCCISO: Final = frozenset({-signal.SIGKILL, 128 + signal.SIGKILL})


class Condizione(StrEnum):
    """Come si è concluso il processo R, secondo il ponte."""

    RIUSCITO = "riuscito"
    ERRORE_DICHIARATO = "errore_dichiarato"
    MEMORIA_ESAURITA = "memoria_esaurita"
    ERRORE_NON_CATALOGATO = "errore_non_catalogato"
    NESSUNA_DICHIARAZIONE = "nessuna_dichiarazione"
    TEMPO_SCADUTO = "tempo_scaduto"
    NON_AVVIATO = "non_avviato"


@dataclass(frozen=True)
class EsitoR:
    """Esito di uno script R concluso con successo."""

    script: Path
    fase: Fase
    invocazione: str
    codice_uscita: int
    durata_s: float
    #: Artefatti dichiarati dallo script, gia' registrati nel manifesto.
    artefatti: tuple[Artefatto, ...]
    stdout: str
    stderr: str

    @property
    def riuscito(self) -> bool:
        """Vero se il processo R è uscito con codice 0."""
        return self.codice_uscita == 0


# --------------------------------------------------------------------------- #
# Dove stanno l'interprete e gli script                                        #
# --------------------------------------------------------------------------- #


def trova_rscript(esplicito: Path | str | None = None) -> Path | None:
    """L'interprete R: quello indicato, poi la variabile d'ambiente, poi il PATH.

    Nessun percorso assoluto è scritto qui: nel container ``Rscript`` è nel
    PATH, su una macchina di sviluppo pure, e dove non lo è si indica.
    """
    for percorso_r in (esplicito, os.environ.get(VARIABILE_RSCRIPT)):
        if percorso_r:
            trovato = shutil.which(os.fspath(percorso_r))
            return Path(trovato) if trovato else None
    trovato = shutil.which("Rscript")
    return Path(trovato) if trovato else None


def cartella_r() -> Path:
    """La cartella ``R/`` con gli script di fase e le funzioni condivise.

    Per difetto è quella del repository, che con l'installazione in modalità
    sviluppo sta accanto a ``src/``; altrove si indica con
    ``AMPLICON16S_R_DIR``.
    """
    indicata = os.environ.get(VARIABILE_CARTELLA_R)
    if indicata:
        return Path(indicata)
    return Path(__file__).resolve().parents[3] / "R"


# --------------------------------------------------------------------------- #
# Classificazione dell'esito                                                   #
# --------------------------------------------------------------------------- #


def _memoria_esaurita(testo: str) -> bool:
    """Vero se il testo contiene uno dei messaggi con cui R segnala la memoria esaurita.
    """
    return bool(_MEMORIA_ESAURITA.search(testo))


def classifica(
    codice_uscita: int,
    dichiarazione: Dichiarazione | None,
    stderr: str,
) -> Condizione:
    """Stabilisce la condizione dall'esito del processo.

    È una funzione pura: tutto ciò che il ponte sa del processo concluso sta
    nei tre argomenti. Una dichiarazione non valida si passa come ``None``.
    """
    if dichiarazione is not None:
        if dichiarazione.stato is Stato.ERRORE_CATALOGO:
            return Condizione.ERRORE_DICHIARATO
        if dichiarazione.stato is Stato.ERRORE_NON_CATALOGATO:
            # Il messaggio dell'errore, oppure l'avviso dei figli uccisi che lo
            # precede sull'uscita di errore: l'errore dichiarato e' allora solo
            # la conseguenza di un risultato mancante.
            if _memoria_esaurita(dichiarazione.messaggio) or _FIGLI_SENZA_RISULTATO.search(stderr):
                return Condizione.MEMORIA_ESAURITA
            return Condizione.ERRORE_NON_CATALOGATO
        if codice_uscita == 0:
            # Un figlio del calcolo parallelo morto senza risultato non fa
            # fallire tutti gli script: alcune funzioni dei pacchetti proseguono
            # con i risultati che restano, e l'esito sarebbe un artefatto
            # sbagliato dichiarato riuscito. Il successo non e' credibile.
            if _FIGLI_SENZA_RISULTATO.search(stderr):
                return Condizione.MEMORIA_ESAURITA
            return Condizione.RIUSCITO
        # Successo dichiarato ma uscita non nulla: il processo e' morto dopo
        # aver scritto l'esito, e il successo non e' piu' credibile. Si decide
        # come se la dichiarazione non ci fosse.

    if codice_uscita in _UCCISO or _memoria_esaurita(stderr):
        return Condizione.MEMORIA_ESAURITA
    return Condizione.NESSUNA_DICHIARAZIONE


def _descrivi_uscita(codice_uscita: int) -> dict[str, Any]:
    """Il codice di uscita e, se il processo è stato terminato da un segnale, il nome
    del segnale.
    """
    descrizione: dict[str, Any] = {"codice_uscita": codice_uscita}
    if codice_uscita < 0:
        try:
            descrizione["segnale"] = signal.Signals(-codice_uscita).name
        except ValueError:
            descrizione["segnale"] = -codice_uscita
    return descrizione


def _coda(testo: str, limite: int = MAX_CARATTERI_USCITA) -> tuple[str, bool]:
    """Gli ultimi ``limite`` caratteri del testo, e se è stato troncato."""
    if len(testo) <= limite:
        return testo, False
    return testo[-limite:], True


# --------------------------------------------------------------------------- #
# Esecuzione                                                                   #
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class _Concluso:
    """Ciò che resta di un processo concluso: codice di uscita, uscite e scadenza del
    tempo.
    """
    codice_uscita: int
    stdout: str
    stderr: str
    scaduto: bool


def _lancia(
    comando: list[str],
    cartella: Path,
    ambiente: Mapping[str, str],
    tempo_massimo_s: float | None,
) -> _Concluso:
    """Avvia il processo in una sessione propria e ne attende la fine, uccidendo il
    gruppo se scade il tempo.
    """
    # Il figlio apre una propria sessione: allo scadere del tempo si uccide
    # l'intero gruppo, compresi i processi che R avesse generato a sua volta.
    processo = subprocess.Popen(
        comando,
        cwd=cartella,
        env=dict(ambiente),
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        encoding="utf-8",
        errors="replace",
        start_new_session=True,
    )
    try:
        stdout, stderr = processo.communicate(timeout=tempo_massimo_s)
        scaduto = False
    except subprocess.TimeoutExpired:
        _uccidi_gruppo(processo)
        stdout, stderr = processo.communicate()
        scaduto = True
    except BaseException:
        # L'orchestratore viene interrotto (Ctrl-C, un segnale, un errore): il
        # figlio ha una sessione propria e gli sopravvivrebbe, continuando a
        # scrivere nella cartella di fase mentre una ripresa ne avvia un altro.
        _uccidi_gruppo(processo)
        processo.wait()
        raise
    return _Concluso(processo.returncode, stdout, stderr, scaduto)


def _uccidi_gruppo(processo: subprocess.Popen[str]) -> None:
    """Uccide il processo R e tutti quelli della sua sessione."""
    try:
        os.killpg(processo.pid, signal.SIGKILL)
    except ProcessLookupError:
        pass


def _registra_uscita(
    logger: logging.Logger,
    flusso: str,
    testo: str,
    livello: int,
    contesto: Mapping[str, Any],
) -> None:
    """Registra nel log la coda di un flusso di uscita del processo R, se non è vuoto.
    """
    if not testo:
        return
    coda, troncato = _coda(testo)
    logger.log(
        livello,
        f"{flusso} del processo R",
        extra={**contesto, "flusso": flusso, "testo": coda, "troncato": troncato},
    )


def esegui_script(
    script: Path | str,
    parametri: Mapping[str, Any],
    albero: AlberoOutput,
    fase: Fase,
    *,
    passo: str,
    codice_memoria: str | None = None,
    tempo_massimo_s: float | None = None,
    rscript: Path | str | None = None,
    logger: logging.Logger | None = None,
) -> EsitoR:
    """Esegue uno script R e ne traduce l'esito.

    ``parametri`` arrivano allo script come oggetto JSON. Lo script scrive i
    propri artefatti nella cartella di ``fase``; quelli che dichiara vengono
    registrati nel manifesto. ``passo`` è la fase che invoca lo script: dà il
    nome ai file del contratto, distinti anche in una cartella condivisa. ``codice_memoria`` è il codice del catalogo che
    la fase associa alla memoria esaurita; ``tempo_massimo_s``
    (``run.r_timeout_s``) è la durata oltre la quale il processo viene ucciso
    con tutto il suo gruppo (``E-R-05``), nullo per nessun limite.

    Restituisce l'esito se lo script è riuscito, altrimenti solleva un
    :class:`~amplicon16s.errors.exceptions.ErrorePipeline`. Qualunque sia il
    modo in cui il processo R termina, il processo chiamante prosegue.
    """
    if codice_memoria is not None:
        voce(codice_memoria)  # un codice inesistente e' un difetto del chiamante

    logger = logger or ottieni("rbridge")
    script = Path(script)
    invocazione = uuid.uuid4().hex
    contesto: dict[str, Any] = {
        "script": script.name,
        # Non "fase": nell'evento di un errore quel campo e' la fase del
        # catalogo, e questa e' la cartella in cui lo script ha lavorato.
        "cartella": fase.value,
        "passo": str(passo),
        "invocazione": invocazione,
    }

    interprete = trova_rscript(rscript)
    lib = cartella_r() / "lib"
    if interprete is None:
        raise errore(
            "E-R-01",
            "Rscript non trovato nel PATH ne' indicato",
            **contesto,
            condizione_r=Condizione.NON_AVVIATO.value,
        )
    for richiesto in (script, lib):
        if not richiesto.exists():
            raise errore(
                "E-R-01",
                f"non trovato: {richiesto}",
                **contesto,
                condizione_r=Condizione.NON_AVVIATO.value,
            )

    cartella = albero.prepara(fase)
    richiesta, percorso_esito = scrivi_richiesta(cartella, invocazione, parametri, str(passo))

    ambiente = {**os.environ, "LANGUAGE": "en", "LC_ALL": "C.UTF-8", VARIABILE_LIB_R: str(lib)}
    comando = [str(interprete), "--vanilla", str(script.resolve()), str(richiesta)]

    logger.info(
        "avvio del processo R",
        extra={**contesto, "tempo_massimo_s": tempo_massimo_s},
    )
    inizio = time.monotonic()
    try:
        processo = _lancia(comando, cartella, ambiente, tempo_massimo_s)
    except OSError as e:
        raise errore(
            "E-R-01",
            f"l'interprete {interprete} non si avvia: {e}",
            **contesto,
            condizione_r=Condizione.NON_AVVIATO.value,
        ) from e
    durata = round(time.monotonic() - inizio, 3)

    uscita = _descrivi_uscita(processo.codice_uscita)
    motivo_non_valida: str | None = None
    try:
        dichiarazione = leggi_dichiarazione(percorso_esito, invocazione)
    except DichiarazioneNonValida as e:
        dichiarazione, motivo_non_valida = None, str(e)

    if processo.scaduto:
        # Ucciso dal ponte: il SIGKILL non e' del sistema operativo, e
        # qualunque cosa lo script abbia dichiarato prima non e' un esito.
        condizione = Condizione.TEMPO_SCADUTO
        motivo_non_valida = f"tempo massimo di {tempo_massimo_s} s (run.r_timeout_s) superato"
    else:
        condizione = classifica(
            processo.codice_uscita, dichiarazione, processo.stderr
        )

    # Un successo e' credibile solo se cio' che dichiara esiste davvero.
    artefatti: tuple[Artefatto, ...] = ()
    if condizione is Condizione.RIUSCITO:
        assert dichiarazione is not None
        mancanti = [n for n in dichiarazione.artefatti if not (cartella / n).is_file()]
        if mancanti:
            condizione = Condizione.NESSUNA_DICHIARAZIONE
            motivo_non_valida = f"artefatti dichiarati ma assenti: {mancanti}"
        else:
            artefatti = tuple(albero.registra(fase, n) for n in dichiarazione.artefatti)

    riuscito = condizione is Condizione.RIUSCITO
    livello = logging.DEBUG if riuscito else logging.WARNING
    _registra_uscita(logger, "stdout", processo.stdout, logging.DEBUG, contesto)
    _registra_uscita(logger, "stderr", processo.stderr, livello, contesto)
    logger.log(
        logging.INFO if riuscito else logging.WARNING,
        "processo R concluso",
        extra={
            **contesto,
            **uscita,
            "durata_s": durata,
            "condizione_r": condizione.value,
            "artefatti": [a.nome for a in artefatti],
        },
    )

    if riuscito:
        return EsitoR(
            script=script,
            fase=fase,
            invocazione=invocazione,
            codice_uscita=processo.codice_uscita,
            durata_s=durata,
            artefatti=artefatti,
            stdout=processo.stdout,
            stderr=processo.stderr,
        )

    raise _eccezione(
        condizione,
        dichiarazione,
        motivo_non_valida,
        codice_memoria,
        {**contesto, **uscita, "condizione_r": condizione.value},
    )


def _eccezione(
    condizione: Condizione,
    dichiarazione: Dichiarazione | None,
    motivo_non_valida: str | None,
    codice_memoria: str | None,
    contesto: dict[str, Any],
) -> ErrorePipeline:
    """L'eccezione del catalogo che corrisponde alla condizione di uscita del processo
    R.
    """
    messaggio = dichiarazione.messaggio if dichiarazione else ""

    if condizione is Condizione.ERRORE_DICHIARATO:
        assert dichiarazione is not None and dichiarazione.codice is not None
        try:
            return errore(dichiarazione.codice, messaggio, **contesto)
        except KeyError:
            return errore(
                "E-R-03",
                f"lo script ha dichiarato {dichiarazione.codice}, che non e' nel "
                f"catalogo: {messaggio}",
                **contesto,
                codice_dichiarato=dichiarazione.codice,
            )

    if condizione is Condizione.MEMORIA_ESAURITA:
        return errore(
            codice_memoria or "E-R-04",
            messaggio or "processo R ucciso dal sistema operativo (SIGKILL)",
            **contesto,
        )

    if condizione is Condizione.ERRORE_NON_CATALOGATO:
        return errore("E-R-03", messaggio, **contesto)

    if condizione is Condizione.TEMPO_SCADUTO:
        return errore("E-R-05", f"processo ucciso: {motivo_non_valida}", **contesto)

    dettaglio = "nessuna dichiarazione d'esito"
    if motivo_non_valida:
        dettaglio = f"dichiarazione d'esito non valida: {motivo_non_valida}"
    return errore("E-R-02", dettaglio, **contesto)
