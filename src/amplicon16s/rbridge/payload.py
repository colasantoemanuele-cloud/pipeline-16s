"""Il contratto fra Python e R: i due file che attraversano il confine.

Nessun oggetto passa in memoria fra i due processi. Python scrive la
**richiesta** in un file JSON e lancia lo script; lo script scrive i propri
artefatti nella cartella di fase e, per ultima cosa, la **dichiarazione
d'esito** in un secondo file JSON.

Il secondo file esiste perché il codice di uscita di un processo è un numero
fra 0 e 255 e non può portare un codice del catalogo come ``E-S2-03``. Il
codice viaggia quindi per il filesystem, come tutto il resto, e il codice di
uscita resta solo una conferma. La distinzione che conta davvero è fra un
processo che *ha dichiarato* qualcosa e uno che *non ha potuto farlo*, e
passa per l'esistenza della dichiarazione:

* la dichiarazione viene scritta su un file temporaneo e poi rinominata, che
  sullo stesso filesystem è un'operazione atomica: un processo ucciso a metà
  scrittura lascia un file temporaneo, mai una dichiarazione troncata;
* il ponte cancella la dichiarazione di un'esecuzione precedente prima di
  lanciare lo script, e la richiesta porta un identificativo di invocazione
  che la dichiarazione deve riportare: un file rimasto da un tentativo
  precedente non può essere scambiato per l'esito di quello attuale.

Una dichiarazione assente, illeggibile o di un'altra invocazione equivale
quindi a un processo che non ha dichiarato nulla.

I due file portano nel nome la fase che li scrive (``rbridge_richiesta_S6.json``,
``rbridge_esito_S6.json``): più fasi possono condividere una cartella, come S6
e S7 in ``07_chimera/``, e con un nome unico l'ultima cancellerebbe la traccia
di come è andata la precedente. Due tentativi della stessa fase usano invece
lo stesso nome, ed è l'identificativo di invocazione a distinguerli.

Forma della richiesta::

    {"protocollo": 1, "invocazione": "…", "cartella_fase": "…",
     "esito": "…/rbridge_esito_S6.json", "parametri": {…}}

Forma della dichiarazione::

    {"protocollo": 1, "invocazione": "…", "stato": "riuscito",
     "codice": null, "messaggio": "", "artefatti": ["nome", …]}

``stato`` vale ``riuscito``, ``errore_catalogo`` (``codice`` è allora il
codice dichiarato) oppure ``errore_non_catalogato`` (un errore R che lo
script non ha ricondotto a un codice; ``messaggio`` è quello di R). Le
funzioni R che scrivono la dichiarazione sono in ``R/lib/errors.R``.
"""

from __future__ import annotations

import json
import os
import tempfile
from collections.abc import Mapping
from dataclasses import dataclass
from enum import StrEnum
from pathlib import Path, PurePosixPath
from typing import Any, Final

__all__ = [
    "PROTOCOLLO",
    "nomi_del_contratto",
    "Dichiarazione",
    "DichiarazioneNonValida",
    "Stato",
    "leggi_dichiarazione",
    "scrivi_richiesta",
]

#: Versione del contratto. Va incrementata a ogni modifica incompatibile della
#: forma dei due file, e ``R/lib/errors.R`` la verifica.
PROTOCOLLO: Final = 1

#: Prefisso comune dei file del contratto: non sono artefatti.
PREFISSO: Final = "rbridge_"


def nomi_del_contratto(passo: str) -> tuple[str, str]:
    """I nomi della richiesta e della dichiarazione d'esito della fase indicata."""
    if not passo or "/" in passo:
        raise ValueError(f"fase non valida per i file del contratto: {passo!r}")
    return f"{PREFISSO}richiesta_{passo}.json", f"{PREFISSO}esito_{passo}.json"


class Stato(StrEnum):
    """Ciò che uno script R può dichiarare di sé."""

    RIUSCITO = "riuscito"
    ERRORE_CATALOGO = "errore_catalogo"
    ERRORE_NON_CATALOGATO = "errore_non_catalogato"


@dataclass(frozen=True)
class Dichiarazione:
    """L'esito come lo script R lo ha dichiarato."""

    stato: Stato
    codice: str | None
    messaggio: str
    #: Nomi degli artefatti, relativi alla cartella di fase.
    artefatti: tuple[str, ...]


class DichiarazioneNonValida(ValueError):
    """La dichiarazione esiste ma non rispetta il contratto."""


def _json_nativo(valore: Any) -> Any:
    """Converte in testo i percorsi; ogni altro tipo non JSON è un errore."""
    # I percorsi sono l'unico tipo non nativo che ha senso passare a R: tutto
    # il resto deve gia' essere JSON, e un oggetto qualunque trasformato in
    # testo arriverebbe dall'altra parte come una stringa senza significato.
    if isinstance(valore, os.PathLike):
        return os.fspath(valore)
    raise TypeError(
        f"parametro non rappresentabile in JSON: {type(valore).__name__}"
    )


def _scrivi_atomico(percorso: Path, documento: Mapping[str, Any]) -> None:
    """Scrive il documento JSON su un file temporaneo nella stessa cartella e lo
    sostituisce al percorso finale.
    """
    # allow_nan=False: NaN e Infinity non sono JSON, e il lettore R li
    # rifiuterebbe o, peggio, li leggerebbe come stringhe.
    testo = json.dumps(
        documento,
        ensure_ascii=False,
        allow_nan=False,
        indent=2,
        default=_json_nativo,
    )
    descrittore, temporaneo = tempfile.mkstemp(
        prefix=".scrittura-", suffix=".json", dir=percorso.parent
    )
    try:
        with os.fdopen(descrittore, "w", encoding="utf-8") as file:
            file.write(testo + "\n")
        os.replace(temporaneo, percorso)
    except BaseException:
        Path(temporaneo).unlink(missing_ok=True)
        raise


def scrivi_richiesta(
    cartella_fase: Path,
    invocazione: str,
    parametri: Mapping[str, Any],
    passo: str,
) -> tuple[Path, Path]:
    """Scrive la richiesta e rimuove la dichiarazione di un'esecuzione precedente.

    Restituisce i percorsi della richiesta e della dichiarazione attesa, con i
    nomi della fase ``passo``.
    """
    nome_richiesta, nome_esito = nomi_del_contratto(passo)
    richiesta = cartella_fase / nome_richiesta
    esito = cartella_fase / nome_esito
    esito.unlink(missing_ok=True)

    _scrivi_atomico(
        richiesta,
        {
            "protocollo": PROTOCOLLO,
            "invocazione": invocazione,
            "cartella_fase": str(cartella_fase),
            "esito": str(esito),
            "parametri": dict(parametri),
        },
    )
    return richiesta, esito


def _nome_relativo(nome: Any) -> str:
    """Il nome dell'artefatto, se è un percorso relativo che non esce dalla cartella di
    fase.
    """
    if not isinstance(nome, str) or not nome:
        raise DichiarazioneNonValida(f"nome di artefatto non valido: {nome!r}")
    forma = PurePosixPath(nome)
    if forma.is_absolute() or ".." in forma.parts:
        raise DichiarazioneNonValida(
            f"l'artefatto {nome!r} non e' relativo alla cartella di fase"
        )
    return nome


def leggi_dichiarazione(percorso: Path, invocazione: str) -> Dichiarazione | None:
    """Legge la dichiarazione d'esito.

    Restituisce ``None`` se il file non esiste: lo script non ha dichiarato
    nulla. Solleva :class:`DichiarazioneNonValida` se il file esiste ma non
    rispetta il contratto, compreso il caso in cui appartenga a un'altra
    invocazione.
    """
    if not percorso.is_file():
        return None

    try:
        documento = json.loads(percorso.read_text(encoding="utf-8"))
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as e:
        raise DichiarazioneNonValida(f"dichiarazione illeggibile: {e}") from e

    if not isinstance(documento, dict):
        raise DichiarazioneNonValida("la dichiarazione non e' un oggetto JSON")
    if documento.get("protocollo") != PROTOCOLLO:
        raise DichiarazioneNonValida(
            f"protocollo {documento.get('protocollo')!r}, atteso {PROTOCOLLO}"
        )
    if documento.get("invocazione") != invocazione:
        raise DichiarazioneNonValida(
            "la dichiarazione appartiene a un'altra invocazione"
        )

    try:
        stato = Stato(documento.get("stato"))
    except ValueError:
        raise DichiarazioneNonValida(
            f"stato sconosciuto: {documento.get('stato')!r}"
        ) from None

    codice = documento.get("codice")
    if stato is Stato.ERRORE_CATALOGO:
        if not isinstance(codice, str) or not codice:
            raise DichiarazioneNonValida("errore dichiarato senza codice")
    elif codice is not None:
        raise DichiarazioneNonValida(
            f"lo stato {stato.value} non ammette un codice, trovato {codice!r}"
        )

    messaggio = documento.get("messaggio", "")
    if not isinstance(messaggio, str):
        raise DichiarazioneNonValida("il messaggio non e' un testo")

    artefatti = documento.get("artefatti", [])
    if not isinstance(artefatti, list):
        raise DichiarazioneNonValida("gli artefatti non sono un elenco")
    nomi = tuple(_nome_relativo(n) for n in artefatti)
    # Un nome ripetuto finirebbe due volte nel manifesto della fase: e' un
    # difetto dello script, che deve dichiarare ogni artefatto una volta sola.
    ripetuti = sorted({n for n in nomi if nomi.count(n) > 1})
    if ripetuti:
        raise DichiarazioneNonValida(f"artefatti dichiarati piu' volte: {ripetuti}")

    return Dichiarazione(
        stato=stato,
        codice=codice,
        messaggio=messaggio,
        artefatti=nomi,
    )
