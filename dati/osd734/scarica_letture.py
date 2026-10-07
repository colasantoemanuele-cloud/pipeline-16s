#!/usr/bin/env python3
"""Scarica da ENA le 960 letture del sequenziamento 16S di OSD-734.

L'elenco, con l'indirizzo, l'MD5 e la dimensione di ogni file dichiarati da
ENA, e' ``letture_ena.tsv`` accanto a questo script. I file si chiamano
``<ERX>_<nome del campione>_<ERR>.fastq.gz``: un solo accession dell'esperimento
nel nome, quello che la pipeline estrae (``io.accession_regex``).

Lo script e' **ripetibile**: a ogni avvio controlla tutti i file e scarica solo
cio' che manca o non corrisponde all'MD5; uno scarico interrotto riprende dal
punto a cui era arrivato (``dati/scaricamento.py``). Una cartella che contiene gia'
tutti i file integri non provoca alcuno scarico.

    python3 dati/osd734/scarica_letture.py                  # in dati/osd734/fastq/
    python3 dati/osd734/scarica_letture.py --solo-verifica  # controlla e basta
    python3 dati/osd734/scarica_letture.py --cartella ALTRA --accession ERX12083091

L'esito e' 0 se alla fine tutti i file richiesti sono presenti e integri.
"""

from __future__ import annotations

import argparse
import csv
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from scaricamento import Atteso, ErroreScarico, metti_da_parte, scarica, valido  # noqa: E402

QUI = Path(__file__).resolve().parent
ELENCO = QUI / "letture_ena.tsv"
CARTELLA = QUI / "fastq"


def leggi_elenco(percorso: Path | None = None) -> list[dict[str, str]]:
    """Le righe dell'elenco delle letture (predefinito :data:`ELENCO`), nell'ordine
    del file.
    """
    with open(percorso or ELENCO, encoding="utf-8", newline="") as file:
        return list(csv.DictReader(file, delimiter="\t"))


def atteso(riga: dict[str, str]) -> Atteso:
    """Il file atteso per una riga dell'elenco."""
    return Atteso(riga["fastq_url"], riga["fastq_md5"], int(riga["fastq_bytes"]))


def main(argomenti: list[str] | None = None) -> int:
    """Controlla, e se serve scarica, le letture richieste."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--cartella", type=Path, default=CARTELLA,
                        help="dove stanno, o andranno, i file (predefinita: dati/osd734/fastq)")
    parser.add_argument("--solo-verifica", action="store_true",
                        help="controlla soltanto, senza scaricare ne' spostare nulla")
    parser.add_argument("--accession", action="append", default=[],
                        help="limita a questi accession ERX o ERR (ripetibile)")
    parser.add_argument("--pausa", type=float, default=0.5,
                        help="secondi fra una richiesta e la successiva")
    parser.add_argument("--tentativi", type=int, default=4)
    opzioni = parser.parse_args(argomenti)

    righe = leggi_elenco()
    if opzioni.accession:
        scelti = set(opzioni.accession)
        righe = [r for r in righe if {r["experiment_accession"], r["run_accession"]} & scelti]
        if not righe:
            print("nessuna riga dell'elenco corrisponde agli accession indicati", file=sys.stderr)
            return 2

    presenti = scaricati = 0
    mancanti: list[str] = []
    for indice, riga in enumerate(righe, start=1):
        percorso = opzioni.cartella / riga["file"]
        voce = atteso(riga)
        if valido(percorso, voce):
            presenti += 1
            continue
        if opzioni.solo_verifica:
            stato = "non corrispondente" if percorso.exists() else "assente"
            mancanti.append(f"{riga['file']}: {stato}")
            continue
        if percorso.exists():
            print(f"{riga['file']}: MD5 o dimensione diversi, messo da parte in "
                  f"{metti_da_parte(percorso).name}", flush=True)
        try:
            scarica(voce, percorso, tentativi=opzioni.tentativi)
        except ErroreScarico as guasto:
            mancanti.append(str(guasto))
        else:
            scaricati += 1
            print(f"[{indice}/{len(righe)}] {riga['file']}: scaricato e verificato", flush=True)
        time.sleep(opzioni.pausa)

    print(f"file richiesti {len(righe)}: gia' presenti e integri {presenti}, "
          f"scaricati {scaricati}, mancanti o non integri {len(mancanti)}")
    for voce in mancanti:
        print(f"  {voce}")
    return 0 if not mancanti else 1


if __name__ == "__main__":
    sys.exit(main())
