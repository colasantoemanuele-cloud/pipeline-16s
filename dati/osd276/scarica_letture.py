#!/usr/bin/env python3
"""Scarica da ENA le letture di OSD-276 e ne ricava le letture forward.

OSD-276 e' un sequenziamento paired-end. ENA distribuisce, per ciascuna delle
15 corse, **un solo file** che contiene entrambe le letture di ogni coppia, in
due blocchi: i record la cui intestazione termina con ``/1`` (letture forward)
e quelli che terminano con ``/2`` (letture inverse). La pipeline tratta dati
single-end: qui si usano le **sole letture forward**, ed e' questo script a
separarle.

Per ogni corsa:

1. il file di ENA si scarica in ``<cartella>/ena/`` e si verifica con l'MD5 e la
   dimensione dichiarati da ENA (``letture_ena.tsv``, ``dati/scaricamento.py``): e' la sola fonte, ed e' verificata;
2. si contano i record ``/1``, i record ``/2`` e gli altri. Il file si rifiuta
   se un record non e' ne' ``/1`` ne' ``/2``, se le ``/1`` e le ``/2`` non sono
   in numero uguale, o se la loro somma non e' il numero di letture dichiarato
   dal deposito (``letture_deposito``): in ciascuno dei tre casi il file non e'
   quello che questo script sa separare, e scriverne una parte darebbe letture
   forward plausibili e sbagliate;
3. i record ``/1`` si scrivono, nell'ordine e senza modifiche, in
   ``<cartella>/fastq/<SRR>_<campione>_R1.fastq.gz``; l'MD5 del contenuto non
   compresso deve essere quello di ``letture_forward.md5`` accanto a questo
   script. Un file con un'impronta diversa non resta in ``fastq/``, dove la
   pipeline lo leggerebbe: si mette da parte con il suffisso ``.md5_errato``.

Lo script e' **ripetibile**: un file di ENA gia' integro non si scarica di
nuovo, e un file forward il cui contenuto corrisponde non si riscrive. Un guasto
su una corsa (deposito illeggibile, conteggi incoerenti) non ferma le altre: e'
raccolto nel riepilogo finale, e l'esito e' 1.

    python3 dati/osd276/scarica_letture.py                  # in dati/osd276/
    python3 dati/osd276/scarica_letture.py --solo-verifica  # controlla e basta
    python3 dati/osd276/scarica_letture.py --cartella ALTRA

L'esito e' 0 se alla fine tutti i file forward sono presenti e conformi.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import os
import sys
import tempfile
import time
import zlib
from pathlib import Path

QUI = Path(__file__).resolve().parent
sys.path.insert(0, str(QUI.parent))
from scaricamento import Atteso, ErroreScarico, metti_da_parte, scarica, valido  # noqa: E402

ELENCO = QUI / "letture_ena.tsv"
IMPRONTE = QUI / "letture_forward.md5"


def leggi_elenco() -> list[dict[str, str]]:
    """Le righe dell'elenco delle corse, nell'ordine del file."""
    with open(ELENCO, encoding="utf-8", newline="") as file:
        return list(csv.DictReader(file, delimiter="\t"))


def leggi_impronte() -> dict[str, str]:
    """L'MD5 atteso del contenuto non compresso di ogni file forward."""
    impronte = {}
    for riga in IMPRONTE.read_text(encoding="utf-8").splitlines():
        md5, nome = riga.split()
        impronte[nome] = md5
    return impronte


def impronta_contenuto(percorso: Path) -> tuple[str, int]:
    """MD5 del contenuto non compresso di un FASTQ e numero dei suoi record."""
    impronta = hashlib.md5()
    righe = 0
    with gzip.open(percorso, "rb") as file:
        for riga in file:
            impronta.update(riga)
            righe += 1
    return impronta.hexdigest(), righe // 4


class DepositoInatteso(ValueError):
    """Il file del deposito non ha la forma che lo script sa separare."""


#: I guasti di una corsa che finiscono nel riepilogo senza fermare le altre:
#: scarico fallito, deposito inatteso, archivio illeggibile o troncato.
GUASTI = (ErroreScarico, DepositoInatteso, OSError, EOFError, zlib.error)


def conta_letture(deposito: Path) -> tuple[int, int, int]:
    """Quanti record di ``deposito`` terminano con ``/1``, quanti con ``/2`` e
    quanti con altro.
    """
    prime = seconde = altre = 0
    with gzip.open(deposito, "rb") as ingresso:
        while intestazione := ingresso.readline():
            if not (ingresso.readline() and ingresso.readline() and ingresso.readline()):
                raise DepositoInatteso(f"{deposito.name}: record incompleto")
            fine = intestazione.rstrip()[-2:]
            prime += fine == b"/1"
            seconde += fine == b"/2"
            altre += fine not in (b"/1", b"/2")
    return prime, seconde, altre


def verifica_conteggi(deposito: Path, letture_deposito: int) -> int:
    """Verifica che ``deposito`` contenga le due letture di ogni coppia e nulla
    d'altro, e restituisce il numero delle forward attese.

    Solleva :class:`DepositoInatteso` se un record non e' ne' ``/1`` ne' ``/2``,
    se le ``/1`` e le ``/2`` sono in numero diverso, o se la loro somma non e'
    il numero di letture dichiarato dal deposito.
    """
    prime, seconde, altre = conta_letture(deposito)
    conteggi = f"{prime} record /1, {seconde} record /2, {altre} altri"
    if altre:
        raise DepositoInatteso(f"{deposito.name}: {conteggi}: ci sono record che non sono "
                               "ne' la prima ne' la seconda lettura di una coppia")
    if prime != seconde:
        raise DepositoInatteso(f"{deposito.name}: {conteggi}: le due letture delle coppie "
                               "non sono in numero uguale")
    if prime + seconde != letture_deposito:
        raise DepositoInatteso(f"{deposito.name}: {conteggi}: la somma non e' il numero di "
                               f"letture dichiarato dal deposito ({letture_deposito})")
    return prime


def ricava_forward(deposito: Path, destinazione: Path) -> int:
    """Scrive in ``destinazione`` i record ``/1`` di ``deposito``; ne restituisce
    il numero. La scrittura e' atomica, e l'archivio non porta ne' data ne' nome.
    """
    destinazione.parent.mkdir(parents=True, exist_ok=True)
    descrittore, temporaneo = tempfile.mkstemp(dir=destinazione.parent, prefix=".scrittura-")
    scritti = 0
    try:
        with (
            gzip.open(deposito, "rb") as ingresso,
            os.fdopen(descrittore, "wb") as grezzo,
            gzip.GzipFile(filename="", mode="wb", fileobj=grezzo, mtime=0) as uscita,
        ):
            while intestazione := ingresso.readline():
                record = [intestazione, ingresso.readline(), ingresso.readline(), ingresso.readline()]
                if not record[3]:
                    raise DepositoInatteso(f"{deposito.name}: record incompleto")
                if intestazione.rstrip().endswith(b"/1"):
                    uscita.writelines(record)
                    scritti += 1
        os.replace(temporaneo, destinazione)
    except BaseException:
        Path(temporaneo).unlink(missing_ok=True)
        raise
    return scritti


def prepara_corsa(riga: dict[str, str], cartella: Path, attese: tuple[str, int],
                  tentativi: int, pausa: float) -> None:
    """Scarica se serve il file di una corsa, lo verifica e ne ricava le letture
    forward. Solleva un'eccezione con il motivo se la corsa non va a buon fine.
    """
    forward = cartella / "fastq" / riga["file"]
    deposito = cartella / "ena" / f"{riga['run_accession']}.fastq.gz"
    voce = Atteso(riga["fastq_url"], riga["fastq_md5"], int(riga["fastq_bytes"]))
    if not valido(deposito, voce):
        if deposito.exists():
            print(f"{deposito.name}: MD5 o dimensione diversi, messo da parte in "
                  f"{metti_da_parte(deposito).name}", flush=True)
        scarica(voce, deposito, tentativi=tentativi)
        time.sleep(pausa)
    previste = verifica_conteggi(deposito, int(riga["letture_deposito"]))
    if previste != attese[1]:
        raise DepositoInatteso(f"{deposito.name}: {previste} letture forward, attese {attese[1]}")
    ricava_forward(deposito, forward)
    if impronta_contenuto(forward) != attese:
        # Non resta dove la pipeline lo leggerebbe.
        raise DepositoInatteso(
            f"{riga['file']}: le letture forward ricavate non sono quelle attese; file messo "
            f"da parte in {metti_da_parte(forward).name}")


def main(argomenti: list[str] | None = None) -> int:
    """Controlla, e se serve scarica e ricava, le letture forward."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--cartella", type=Path, default=QUI,
                        help="dove stanno, o andranno, ena/ e fastq/ (predefinita: dati/osd276)")
    parser.add_argument("--solo-verifica", action="store_true",
                        help="controlla soltanto, senza scaricare ne' scrivere nulla")
    parser.add_argument("--pausa", type=float, default=0.5,
                        help="secondi fra una richiesta e la successiva")
    parser.add_argument("--tentativi", type=int, default=4)
    opzioni = parser.parse_args(argomenti)

    righe = leggi_elenco()
    impronte = leggi_impronte()
    conformi = 0
    mancanti: list[str] = []
    for indice, riga in enumerate(righe, start=1):
        forward = opzioni.cartella / "fastq" / riga["file"]
        attese = (impronte[riga["file"]], int(riga["letture_forward"]))
        # Un file forward illeggibile (troncato, non compresso) e' non conforme
        # come uno dal contenuto diverso: non deve fermare le altre corse.
        try:
            conforme = forward.is_file() and impronta_contenuto(forward) == attese
        except GUASTI as guasto:
            conforme = False
            print(f"{riga['file']}: illeggibile ({type(guasto).__name__}: {guasto})", flush=True)
        if conforme:
            conformi += 1
            continue
        if opzioni.solo_verifica:
            mancanti.append(f"{riga['file']}: {'non conforme' if forward.exists() else 'assente'}")
            continue

        # Ogni guasto di una corsa finisce nel riepilogo: le altre proseguono.
        try:
            prepara_corsa(riga, opzioni.cartella, attese, opzioni.tentativi, opzioni.pausa)
        except GUASTI as guasto:
            mancanti.append(f"{riga['file']}: {type(guasto).__name__}: {guasto}")
            continue
        conformi += 1
        print(f"[{indice}/{len(righe)}] {riga['file']}: {attese[1]} letture forward", flush=True)

    print(f"corse {len(righe)}: file forward conformi {conformi}, "
          f"mancanti o non conformi {len(mancanti)}")
    for voce in mancanti:
        print(f"  {voce}")
    return 0 if not mancanti else 1


if __name__ == "__main__":
    sys.exit(main())
