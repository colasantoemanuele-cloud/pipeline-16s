#!/usr/bin/env python3
"""Scarica il riferimento tassonomico usato con OSD-734: SILVA 138 per dada2.

Tre file, ciascuno con il checksum della sua fonte:

* ``silva_nr99_v138_train_set.fa.gz``, il training set di ``assignTaxonomy``
  (``tax.ref_fasta``; l'MD5 e' ``tax.ref_md5``), dal record Zenodo 3986799,
  versione 2 (DOI 10.5281/zenodo.3986799), licenza CC BY 4.0;
* ``SILVA_LICENSE.txt``, la licenza dichiarata dallo stesso record;
* ``silva_138_v2_bad-taxa.csv``, l'elenco dei taxa con un rango mancante nella
  versione 2 (``tax.ref_bad_taxa``), dal repository degli autori dei file
  (mikemc/dada2-reference-databases, licenza MIT), fissato a un commit.

Come ``scarica_letture.py``, lo script e' ripetibile e riprende uno scarico
interrotto; un file gia' presente e integro non viene scaricato di nuovo.

    python3 dati/osd734/scarica_riferimento.py                  # in dati/osd734/riferimento/
    python3 dati/osd734/scarica_riferimento.py --cartella ALTRA
    python3 dati/osd734/scarica_riferimento.py --solo-verifica
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Final

sys.path.insert(0, str(Path(__file__).resolve().parent))
from scaricamento import Atteso, ErroreScarico, metti_da_parte, scarica, valido  # noqa: E402

QUI = Path(__file__).resolve().parent
CARTELLA = QUI / "riferimento"

ZENODO: Final = "https://zenodo.org/records/3986799/files"
#: Il commit del repository mikemc/dada2-reference-databases da cui viene
#: l'elenco dei taxa difettosi: fissato, perche' il ramo principale puo' cambiare.
COMMIT_BAD_TAXA: Final = "73d8ddf62ed322c5d47cbcf0bf529d64ba31700d"

FILE: Final[dict[str, Atteso]] = {
    "silva_nr99_v138_train_set.fa.gz": Atteso(
        f"{ZENODO}/silva_nr99_v138_train_set.fa.gz?download=1",
        "56821dd6f365f64c8518a21e049c3232", 137973851,
    ),
    "SILVA_LICENSE.txt": Atteso(
        f"{ZENODO}/SILVA_LICENSE.txt?download=1",
        "5311c80d208d3ae43d5dc0136df2fcab", 449,
    ),
    "silva_138_v2_bad-taxa.csv": Atteso(
        "https://raw.githubusercontent.com/mikemc/dada2-reference-databases/"
        f"{COMMIT_BAD_TAXA}/silva-138/v2/bad-taxa.csv",
        "5d919588f2c291367330ad515883e3c3", 12309,
    ),
}


def main(argomenti: list[str] | None = None) -> int:
    """Controlla, e se serve scarica, i file del riferimento."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--cartella", type=Path, default=CARTELLA,
                        help="dove stanno, o andranno, i file (predefinita: dati/osd734/riferimento)")
    parser.add_argument("--solo-verifica", action="store_true",
                        help="controlla soltanto, senza scaricare ne' spostare nulla")
    opzioni = parser.parse_args(argomenti)

    mancanti = []
    for nome, voce in FILE.items():
        percorso = opzioni.cartella / nome
        if valido(percorso, voce):
            print(f"{nome}: presente e integro")
            continue
        if opzioni.solo_verifica:
            mancanti.append(f"{nome}: {'non corrispondente' if percorso.exists() else 'assente'}")
            continue
        if percorso.exists():
            print(f"{nome}: MD5 o dimensione diversi, messo da parte in "
                  f"{metti_da_parte(percorso).name}")
        try:
            scarica(voce, percorso)
        except ErroreScarico as guasto:
            mancanti.append(str(guasto))
        else:
            print(f"{nome}: scaricato e verificato")
    for voce in mancanti:
        print(f"  {voce}")
    return 0 if not mancanti else 1


if __name__ == "__main__":
    sys.exit(main())
