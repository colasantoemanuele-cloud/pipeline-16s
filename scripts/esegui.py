#!/usr/bin/env python3
"""Esegue la pipeline dal codice di questo repository.

Uso, dalla radice di un clone (gli argomenti sono quelli di ``amplicon16s``)::

    python3 scripts/esegui.py run --config dati/osd734/config_osd734.yaml

Il comando ``amplicon16s`` esegue il pacchetto installato nell'ambiente: dentro
l'immagine e' il codice copiato quando l'immagine e' stata costruita, con gli
script R indicati dalla sua variabile ``AMPLICON16S_R_DIR``. Questo script
esegue invece il codice del clone in cui si trova, ``src/`` e ``R/``, senza che
chi lancia debba impostare alcuna variabile d'ambiente: e' cio' che la regola
rigorosa sulla provenienza (``run.strict_provenance``) richiede, perche' il
commit del clone identifichi davvero il codice che calcola.

Sistema anche la cartella personale: in un container avviato con l'utente di
chi lancia (``-u "$(id -u):$(id -g)"``) ``HOME`` punta spesso a una cartella di
un altro utente, non scrivibile, e i programmi che vi tengono una cache
fallirebbero o cambierebbero comportamento. In quel caso ``HOME`` diventa una
cartella temporanea, eliminata alla fine.
"""

from __future__ import annotations

import os
import sys
import tempfile
from pathlib import Path

RADICE = Path(__file__).resolve().parents[1]


def main() -> int:
    """Prepara l'ambiente del processo ed esegue la riga di comando del clone."""
    # Il codice del clone prima di ogni copia installata, e i suoi script R al
    # posto di quelli che l'ambiente indicasse.
    sys.path.insert(0, str(RADICE / "src"))
    os.environ["AMPLICON16S_R_DIR"] = str(RADICE / "R")
    # Nessun file compilato accanto ai sorgenti: il clone resta com'e'.
    sys.dont_write_bytecode = True

    casa = os.environ.get("HOME")
    temporanea = None
    if not casa or not os.access(casa, os.W_OK):
        temporanea = tempfile.TemporaryDirectory(prefix="amplicon16s-home-")
        os.environ["HOME"] = temporanea.name
    try:
        from amplicon16s.cli import main as riga_di_comando

        return riga_di_comando(sys.argv[1:], programma=f"python3 {sys.argv[0]}")
    finally:
        if temporanea is not None:
            temporanea.cleanup()


if __name__ == "__main__":
    sys.exit(main())
