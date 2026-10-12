#!/usr/bin/env python3
"""Scarica un classificatore addestrato per IdTaxa (``tax.idtaxa_trained``).

Con ``tax.classifier: idtaxa`` la pipeline non usa un riferimento FASTA ma un
classificatore gia' addestrato nel formato di DECIPHER: un file ``.RData`` con
l'oggetto ``trainingSet`` di ``LearnTaxa``. Lo script ne scarica due, entrambi
di SILVA, addestrati dagli autori di DECIPHER:

* ``2019``: ``SILVA_SSU_r138_2019.RData`` (SILVA SSU r138). E' il riferimento
  delle versioni A e B della pipeline di NASA GeneLab (GL-DPPD-7104): serve a
  chi confronta i propri risultati con tabelle elaborate con quelle versioni;
* ``138.2``: ``SILVA_SSU_r138_2_2024.RData`` (SILVA SSU r138.2): quello da
  usare per un'analisi nuova.

**La fonte.** Il sito di DECIPHER non serve il file del 2019 e non pubblica
checksum dei suoi file. Le copie scaricate qui sono quelle archiviate su
figshare dal gruppo di elaborazione dati di GeneLab (licenza CC BY 4.0), con la
dimensione e l'MD5 che figshare pubblica per ciascuna: che la copia del 2019
sia il file allora distribuito da DECIPHER lo dichiara chi l'ha archiviata.

Come ``scarica_riferimento.py``, lo script e' ripetibile e riprende uno
scarico interrotto; un file gia' presente e integro non viene scaricato di
nuovo, e uno con il contenuto sbagliato e' messo da parte.

    python3 dati/riferimento/scarica_addestrato.py 138.2            # in dati/riferimento/
    python3 dati/riferimento/scarica_addestrato.py 2019 --cartella ALTRA
    python3 dati/riferimento/scarica_addestrato.py 2019 138.2 --solo-verifica
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Final

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from scaricamento import Atteso, ErroreScarico, metti_da_parte, scarica, valido  # noqa: E402

CARTELLA = Path(__file__).resolve().parent

#: Per ogni classificatore: il nome del file e, dalla scheda figshare della sua
#: copia archiviata, indirizzo, MD5 e dimensione. I DOI sono in FONTI.tsv.
FILE: Final[dict[str, tuple[str, Atteso]]] = {
    "2019": ("SILVA_SSU_r138_2019.RData", Atteso(
        "https://ndownloader.figshare.com/files/46245217",
        "cb983b6a5e8cdb46f8c88b5afae21f66", 151849504,
    )),
    "138.2": ("SILVA_SSU_r138_2_2024.RData", Atteso(
        "https://ndownloader.figshare.com/files/52846199",
        "4e272e39c2d71f5d3e7a31b00dbb1df4", 298719756,
    )),
}


def main(argomenti: list[str] | None = None) -> int:
    """Controlla, e se serve scarica, i classificatori addestrati richiesti."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("quali", nargs="+", choices=sorted(FILE),
                        help="quale classificatore: 2019 (SILVA r138) o 138.2 (SILVA r138.2)")
    parser.add_argument("--cartella", type=Path, default=CARTELLA,
                        help="dove stanno, o andranno, i file (predefinita: dati/riferimento)")
    parser.add_argument("--solo-verifica", action="store_true",
                        help="controlla soltanto, senza scaricare ne' spostare nulla")
    opzioni = parser.parse_args(argomenti)

    mancanti = []
    for quale in dict.fromkeys(opzioni.quali):
        nome, voce = FILE[quale]
        percorso = opzioni.cartella / nome
        if valido(percorso, voce):
            print(f"{nome}: presente e integro")
        elif opzioni.solo_verifica:
            mancanti.append(f"{nome}: {'non corrispondente' if percorso.exists() else 'assente'}")
            continue
        else:
            if percorso.exists():
                print(f"{nome}: MD5 o dimensione diversi, messo da parte in "
                      f"{metti_da_parte(percorso).name}")
            try:
                scarica(voce, percorso)
            except ErroreScarico as guasto:
                mancanti.append(str(guasto))
                continue
            print(f"{nome}: scaricato e verificato")
        print(f"  tax.idtaxa_trained: {percorso}\n  tax.idtaxa_trained_md5: \"{voce.md5}\"")
    for voce_mancante in mancanti:
        print(f"  {voce_mancante}")
    return 0 if not mancanti else 1


if __name__ == "__main__":
    sys.exit(main())
