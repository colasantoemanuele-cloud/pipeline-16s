#!/usr/bin/env python3
"""Scarica da OSDR le letture forward di OSD-170.

OSD-170 e' un sequenziamento a letture appaiate: OSDR pubblica per ognuno dei
30 campioni due file, R1 e R2. La pipeline tratta dati single-end e usa il solo
file forward di ogni campione. La regola e' in ``dati/letture_forward.py``: si
tiene il file le cui letture iniziano con il primer forward dichiarato nei
metadati; se nessuno dei due lo porta si tiene l'R1 e lo si dichiara. Le
letture di questo dataset sono depositate senza primer: lo script lo riporta
per tutte le corse, e l'inversione non e' verificabile con questa regola.

Il file tenuto si salva in ``<cartella>/fastq/<campione>.fastq.gz``, senza
marcatori di coppia nel nome. OSDR dichiara di ogni file la sola dimensione:
l'MD5 del file tenuto e' registrato in ``letture_osdr.tsv`` dal primo scarico, e
ogni scarico successivo deve dare gli stessi byte.

    python3 dati/osd170/scarica_letture.py                  # in dati/osd170/
    python3 dati/osd170/scarica_letture.py --solo-verifica  # controlla e basta
    python3 dati/osd170/scarica_letture.py --cartella ALTRA

L'esito e' 0 se alla fine tutti i file forward sono presenti e conformi.
"""

from __future__ import annotations

import sys
from pathlib import Path

QUI = Path(__file__).resolve().parent
sys.path.insert(0, str(QUI.parent))
import letture_forward  # noqa: E402

#: Il primer forward dichiarato nella tabella di assay (``Primer Info``).
PRIMER_FORWARD = "TCCTACGGGAGGCAGCAGT"

if __name__ == "__main__":
    sys.exit(letture_forward.main(QUI / "letture_osdr.tsv", PRIMER_FORWARD, QUI,
                                  __doc__.splitlines()[0]))
