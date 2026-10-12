#!/usr/bin/env python3
"""Scarica da OSDR le letture di OSD-191.

OSD-191 e' un sequenziamento a letture appaiate della regione V3-V4 (primer
341F e 805R, Illumina MiSeq 2 x 300): OSDR pubblica per ognuno dei 159 campioni
due file, R1 e R2. La regione e' piu' lunga di una lettura: la pipeline la
copre in modalita' appaiata, con i due file di ogni campione.

    python3 dati/osd191/scarica_letture.py --coppie                  # i due file, in dati/osd191/fastq_coppie/
    python3 dati/osd191/scarica_letture.py --coppie --solo-verifica  # controlla e basta
    python3 dati/osd191/scarica_letture.py --coppie --cartella ALTRA

Con ``--coppie`` i due file di ogni campione si salvano in
``<cartella>/fastq_coppie/`` come ``<campione>_R1.fastq.gz`` e
``<campione>_R2.fastq.gz``, nell'ordine del deposito e senza scambi. Senza
l'opzione lo script scarica in ``<cartella>/fastq/`` il solo file forward di
ogni campione (``<campione>.fastq.gz``), riconosciuto dal primer forward in
testa alle letture (la regola e' in ``dati/letture_forward.py``).

OSDR dichiara di ogni file la sola dimensione: l'MD5 di ciascun file e'
registrato in ``letture_osdr.tsv`` dal primo scarico, e ogni scarico successivo
deve dare gli stessi byte.

L'esito e' 0 se alla fine tutti i file richiesti sono presenti e conformi.
"""

from __future__ import annotations

import sys
from pathlib import Path

QUI = Path(__file__).resolve().parent
sys.path.insert(0, str(QUI.parent))
import letture_forward  # noqa: E402

#: Il primer forward dichiarato nella tabella di assay (``Primer Info``): il 341F.
PRIMER_FORWARD = "CCTACGGGNGGCWGCAG"

if __name__ == "__main__":
    sys.exit(letture_forward.main(QUI / "letture_osdr.tsv", PRIMER_FORWARD, QUI,
                                  __doc__.splitlines()[0]))
