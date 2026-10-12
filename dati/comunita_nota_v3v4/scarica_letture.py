#!/usr/bin/env python3
"""Scarica da ENA le letture della comunita' nota sequenziata sulla regione V3-V4.

Le sei corse dello studio PRJEB30924 sono a letture appaiate (Illumina MiSeq,
2 x 300, primer 341F e 806R): ENA pubblica per ognuna due file, ``_1`` e ``_2``,
con dimensione e MD5 (``letture_ena.tsv``). La regione e' piu' lunga di una
lettura: la pipeline la copre in modalita' appaiata, con i due file di ogni
corsa.

    python3 dati/comunita_nota_v3v4/scarica_letture.py --coppie                  # in dati/comunita_nota_v3v4/fastq_coppie/
    python3 dati/comunita_nota_v3v4/scarica_letture.py --coppie --solo-verifica  # controlla e basta
    python3 dati/comunita_nota_v3v4/scarica_letture.py --coppie --cartella ALTRA

Con ``--coppie`` i due file di ogni corsa si salvano in
``<cartella>/fastq_coppie/`` come ``<corsa>_R1.fastq.gz`` e
``<corsa>_R2.fastq.gz``, nell'ordine del deposito e senza scambi. Senza
l'opzione lo script scarica in ``<cartella>/fastq/`` il solo file forward di
ogni corsa (``<corsa>.fastq.gz``), riconosciuto dal primer forward in testa
alle letture (la regola e' in ``dati/letture_forward.py``).

L'esito e' 0 se alla fine tutti i file richiesti sono presenti e conformi.
"""

from __future__ import annotations

import sys
from pathlib import Path

QUI = Path(__file__).resolve().parent
sys.path.insert(0, str(QUI.parent))
import letture_forward  # noqa: E402

#: Il primer forward 341F (Klindworth et al. 2013).
PRIMER_FORWARD = "CCTACGGGNGGCWGCAG"

if __name__ == "__main__":
    sys.exit(letture_forward.main(QUI / "letture_ena.tsv", PRIMER_FORWARD, QUI,
                                  __doc__.splitlines()[0]))
