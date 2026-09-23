"""Registrazione degli eventi: console per l'operatore, file per la tracciabilità.

Le due uscite servono a due lettori diversi e hanno perciò forme diverse.

* La **console** è per chi guarda l'esecuzione mentre avviene: una riga per
  evento, leggibile, con il minimo indispensabile.
* Il **file** sotto ``99_logs`` è per chi ricostruisce a posteriori cosa è
  successo in un'esecuzione durata ore. Il formato è JSON Lines — un oggetto
  JSON per riga — perché quel lettore non legge, interroga: «tutti gli eventi
  del codice E-S2-03», «tutte le degradazioni», «quanto è durata la fase S4».
  Con righe di testo libero ognuna di quelle domande diventa un'espressione
  regolare da indovinare; con JSON Lines è un filtro.

Il file ruota: un'esecuzione su 960 campioni produce molte righe, e un log che
cresce senza limite finisce per essere il motivo per cui il disco si riempie.

I campi passati in ``extra`` finiscono nel JSON accanto ai campi standard,
quindi un evento si arricchisce senza cambiare formato::

    log.info("fase conclusa", extra={"fase": "S2", "campioni": 803})
"""

from __future__ import annotations

import json
import logging
import logging.handlers
from pathlib import Path
from typing import Any, Final

from amplicon16s.io_layer.artifacts import AlberoOutput, Fase

__all__ = [
    "NOME_FILE_LOG",
    "FormattatoreJsonl",
    "chiudi",
    "configura",
    "ottieni",
    "registra_errore",
]

#: Nome del logger radice del pacchetto. I logger delle fasi sono suoi figli.
NOME_RADICE: Final = "amplicon16s"

#: Nome del file di log strutturato sotto 99_logs.
NOME_FILE_LOG: Final = "pipeline.jsonl"

#: Dimensione oltre la quale il file ruota, e quante rotazioni conservare.
MAX_BYTE: Final = 16 * 1024 * 1024
ROTAZIONI: Final = 5

#: Attributi che ogni LogRecord possiede: tutto il resto è stato aggiunto dal
#: chiamante con ``extra`` ed è quindi un campo dell'evento.
_CAMPI_DI_SERVIZIO: Final[frozenset[str]] = frozenset(
    {
        "args", "asctime", "created", "exc_info", "exc_text", "filename",
        "funcName", "levelname", "levelno", "lineno", "module", "msecs",
        "message", "msg", "name", "pathname", "process", "processName",
        "relativeCreated", "stack_info", "taskName", "thread", "threadName",
    }
)


class FormattatoreJsonl(logging.Formatter):
    """Un oggetto JSON per riga, con i campi aggiunti via ``extra``."""

    def format(self, record: logging.LogRecord) -> str:
        evento: dict[str, Any] = {
            "istante": self.formatTime(record, "%Y-%m-%dT%H:%M:%S%z"),
            "livello": record.levelname,
            "origine": record.name,
            "messaggio": record.getMessage(),
        }

        for chiave, valore in record.__dict__.items():
            if chiave in _CAMPI_DI_SERVIZIO or chiave.startswith("_"):
                continue
            evento[chiave] = valore

        if record.exc_info:
            evento["eccezione"] = self.formatException(record.exc_info)

        # default=str: un Path o un Enum finiscono nel log come testo invece
        # di far fallire la scrittura del log, che sarebbe il modo peggiore di
        # perdere la tracciabilita' proprio quando serve.
        return json.dumps(evento, ensure_ascii=False, default=str)


class _FormattatoreConsole(logging.Formatter):
    """Una riga leggibile, con il codice di errore in evidenza se presente."""

    def format(self, record: logging.LogRecord) -> str:
        codice = getattr(record, "codice", None)
        prefisso = f"[{codice}] " if codice else ""
        return "{:<8} {}{}".format(record.levelname, prefisso, record.getMessage())


def ottieni(nome: str | None = None) -> logging.Logger:
    """Restituisce il logger del pacchetto, o un suo figlio."""
    return logging.getLogger(NOME_RADICE if nome is None else f"{NOME_RADICE}.{nome}")


def configura(
    out_root: Path | str | None = None,
    *,
    livello_console: int = logging.INFO,
    livello_file: int = logging.DEBUG,
    max_byte: int = MAX_BYTE,
    rotazioni: int = ROTAZIONI,
) -> Path | None:
    """Predispone le due uscite e restituisce il percorso del file di log.

    Senza ``out_root`` resta la sola console: è il caso dei comandi che non
    producono artefatti, per i quali non esiste una cartella dove scrivere.

    Chiamarla più volte non accumula uscite: quelle predisposte da una
    chiamata precedente vengono rimosse.
    """
    logger = ottieni()
    chiudi()

    # Il livello del logger deve lasciar passare l'uscita piu' verbosa: sono
    # poi le singole uscite a filtrare.
    livelli = [livello_console] + ([livello_file] if out_root is not None else [])
    logger.setLevel(min(livelli))
    logger.propagate = False

    console = logging.StreamHandler()
    console.setLevel(livello_console)
    console.setFormatter(_FormattatoreConsole())
    logger.addHandler(console)

    if out_root is None:
        return None

    albero = AlberoOutput(out_root)
    percorso = albero.prepara(Fase.LOGS) / NOME_FILE_LOG

    file = logging.handlers.RotatingFileHandler(
        percorso,
        maxBytes=max_byte,
        backupCount=rotazioni,
        encoding="utf-8",
    )
    file.setLevel(livello_file)
    file.setFormatter(FormattatoreJsonl())
    logger.addHandler(file)

    return percorso


def chiudi() -> None:
    """Rimuove e chiude le uscite predisposte.

    Serve soprattutto ai test e alle esecuzioni ripetute nello stesso
    processo: un handler di file lasciato aperto tiene il descrittore e, su
    alcune piattaforme, impedisce di rimuovere la cartella di lavoro.
    """
    logger = ottieni()
    for handler in list(logger.handlers):
        logger.removeHandler(handler)
        handler.close()


def registra_errore(
    logger: logging.Logger,
    errore: Any,
    livello: int = logging.ERROR,
) -> None:
    """Registra un :class:`~amplicon16s.errors.exceptions.ErrorePipeline`.

    I campi del catalogo — codice, fase, categoria, azione — finiscono nel log
    strutturato come campi distinti, così una ricerca per codice o per
    categoria non deve passare dal testo del messaggio.
    """
    logger.log(livello, errore.voce.sintesi, extra=errore.come_evento())
