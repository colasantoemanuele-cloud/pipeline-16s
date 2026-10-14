#!/usr/bin/env python3
"""Scarica da ENA le letture forward della comunita' nota.

Le cinque corse dello studio PRJNA758928 sono a letture appaiate: ENA pubblica
per ognuna due file, ``_1`` e ``_2``, con dimensione e MD5. La configurazione di
base del dataset usa le sole letture forward, in modalita' single-end (ce n'e'
anche una appaiata, per cui lo script scarica i due file con ``--coppie``). La
regola e' in
``dati/letture_forward.py``: si tiene il file le cui letture iniziano con il
primer forward 515F, e una corsa in cui lo porta il secondo file e' dichiarata
invertita.

Il file tenuto si salva in ``<cartella>/fastq/<corsa>.fastq.gz``, senza
marcatori di coppia nel nome, dopo la verifica di dimensione e MD5 dichiarati
da ENA (``letture_ena.tsv``).

    python3 dati/comunita_nota/scarica_letture.py                  # in dati/comunita_nota/
    python3 dati/comunita_nota/scarica_letture.py --solo-verifica  # controlla e basta
    python3 dati/comunita_nota/scarica_letture.py --cartella ALTRA

L'esito e' 0 se alla fine tutti i file forward sono presenti e conformi.
"""

from __future__ import annotations

import sys
from pathlib import Path

QUI = Path(__file__).resolve().parent
sys.path.insert(0, str(QUI.parent))
import letture_forward  # noqa: E402

#: Il primer forward 515F nella forma di Parada et al. 2016, che comprende
#: quella di Caporaso et al. 2011 usata dallo studio.
PRIMER_FORWARD = "GTGYCAGCMGCCGCGGTAA"

if __name__ == "__main__":
    sys.exit(letture_forward.main(QUI / "letture_ena.tsv", PRIMER_FORWARD, QUI,
                                  __doc__.splitlines()[0]))
