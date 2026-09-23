"""Calcolo e verifica dei checksum degli artefatti.

Il checksum non è un ornamento: serve a stabilire se una fase è già stata
completata, cioè se i suoi artefatti esistono *e* sono ancora quelli che quella
fase ha prodotto. Un file presente ma troncato — un'esecuzione interrotta a
metà scrittura — è indistinguibile da un file completo se ci si limita a
guardare se esiste.

La lettura avviene a blocchi: gli artefatti della pipeline arrivano ai
gigabyte e caricarli in memoria per calcolarne l'impronta sarebbe un modo di
esaurirla.
"""

from __future__ import annotations

import hashlib
from pathlib import Path
from typing import Final

__all__ = [
    "ALGORITMO",
    "checksum_bytes",
    "checksum_file",
    "corrisponde",
]

#: Algoritmo usato per tutti i checksum degli artefatti.
ALGORITMO: Final = "sha256"

#: Dimensione dei blocchi di lettura, 1 MiB.
_BLOCCO: Final = 1024 * 1024


def _formatta(digest: str) -> str:
    return f"{ALGORITMO}:{digest}"


def checksum_file(percorso: Path | str) -> str:
    """Checksum di un file, nella forma ``sha256:<esadecimale>``."""
    impronta = hashlib.new(ALGORITMO)
    with open(percorso, "rb") as file:
        while blocco := file.read(_BLOCCO):
            impronta.update(blocco)
    return _formatta(impronta.hexdigest())


def checksum_bytes(dati: bytes) -> str:
    """Checksum di un contenuto già in memoria."""
    return _formatta(hashlib.new(ALGORITMO, dati).hexdigest())


def corrisponde(percorso: Path | str, atteso: str) -> bool:
    """Se il file esiste e il suo checksum è quello atteso.

    Un file assente non è un errore da sollevare ma una risposta: la fase non
    è completa.
    """
    percorso = Path(percorso)
    if not percorso.is_file():
        return False
    return checksum_file(percorso) == atteso
