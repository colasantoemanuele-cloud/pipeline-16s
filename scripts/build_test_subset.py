#!/usr/bin/env python3
"""Ricostruisce il sottoinsieme di prova del dataset di riferimento.

La selezione vive nel repository (``tests/fixtures/osd734/selezione.tsv``):
accession, MD5 del file originale e motivo della scelta di ciascun campione.
I FASTQ no: si ricostruiscono dai dati locali, su due livelli.

* ``completo``: copia i file selezionati, verificandone l'MD5, insieme alle
  righe dei metadati che li riguardano. È il sottoinsieme su cui si sviluppa
  in locale, con le letture intere.
* ``ridotto``: sottocampiona le letture in modo deterministico e scrive la
  versione piccola che vive in ``tests/fixtures/osd734/ridotto/`` e su cui
  gira la catena di integrazione continua, dove i dati locali non ci sono.

**Il sottocampionamento conserva le proprietà che giustificano la
selezione.** Tiene di ogni campione una frazione fissa delle letture, con un
minimo, così il contrasto fra campioni profondi e poveri resta; tiene sempre
le letture alla lunghezza minima del file, altrimenti la lunghezza minima
globale di 137 bp (una lettura sola su decine di migliaia) sparirebbe; non
tocca classi, corse e piastre, che vengono dai metadati. Il seme è fisso, e
da esso e dall'accession deriva quello di ciascun file: aggiungere o togliere
un campione non cambia le letture scelte per gli altri. Anche il file
compresso è deterministico: nessun nome né istante nell'intestazione gzip.

I dati di origine sono pubblici, nell'archivio ENA: la colonna ``run_ena``
della selezione indica il run da cui riscaricare ciascun file.

Uso::

    python scripts/build_test_subset.py completo --config config.yaml --destinazione DIR
    python scripts/build_test_subset.py ridotto  --config config.yaml
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import math
import random
import re
import shutil
import sys
from collections.abc import Iterator
from dataclasses import dataclass
from pathlib import Path
from typing import Final

from amplicon16s.config.schema import Config, carica

RADICE: Final = Path(__file__).resolve().parents[1]
SELEZIONE: Final = RADICE / "tests" / "fixtures" / "osd734" / "selezione.tsv"
RIDOTTO: Final = RADICE / "tests" / "fixtures" / "osd734" / "ridotto"

#: Seme del sottocampionamento. Cambiarlo cambia la versione ridotta.
SEME: Final = 734
#: Frazione delle letture conservata, e minimo per campione: un campione con
#: meno letture del minimo le conserva tutte.
FRAZIONE: Final = 0.02
MINIMO: Final = 100

#: Nomi dei file dei metadati ristretti, dentro ``metadati/``.
NOME_ASSAY: Final = "assay.txt"
NOME_STUDIO: Final = "studio.txt"
NOME_LOTTI: Final = "lotti.tsv"


@dataclass(frozen=True)
class Scelto:
    accession: str
    campione: str
    letture: int
    lunghezza_minima: int
    md5: str


def leggi_selezione(percorso: Path = SELEZIONE) -> list[Scelto]:
    with open(percorso, encoding="utf-8", newline="") as file:
        return [
            Scelto(
                r["accession"], r["campione"], int(r["letture"]),
                int(r["lunghezza_minima"]), r["md5"],
            )
            for r in csv.DictReader(file, delimiter="\t")
        ]


def md5(percorso: Path) -> str:
    impronta = hashlib.md5()
    with open(percorso, "rb") as file:
        while blocco := file.read(1 << 20):
            impronta.update(blocco)
    return impronta.hexdigest()


def _file_di(config: Config, accession: str) -> Path:
    trovati = [
        p for p in Path(config.io.fastq_dir).glob(config.io.fastq_glob)
        if p.name.startswith(f"{accession}_")
    ]
    if len(trovati) != 1:
        raise SystemExit(f"{accession}: attesi un file, trovati {len(trovati)}")
    return trovati[0]


# --------------------------------------------------------------------------- #
# Metadati ristretti                                                           #
# --------------------------------------------------------------------------- #


def _filtra_tabella(sorgente: Path, destinazione: Path, tieni) -> int:
    """Copia l'intestazione e le righe per cui ``tieni(riga)`` è vero, intatte."""
    tenute = 0
    with open(sorgente, encoding="utf-8", newline="") as ingresso, \
            open(destinazione, "w", encoding="utf-8", newline="") as uscita:
        intestazione = ingresso.readline()
        uscita.write(intestazione)
        colonne = next(csv.reader([intestazione], delimiter="\t"))
        for riga in ingresso:
            valori = dict(zip(colonne, next(csv.reader([riga], delimiter="\t"))))
            if tieni(valori):
                uscita.write(riga)
                tenute += 1
    return tenute


def scrivi_metadati(config: Config, scelti: list[Scelto], cartella: Path) -> None:
    """Le righe di assay, studio e lotti che riguardano i campioni scelti."""
    cartella.mkdir(parents=True, exist_ok=True)
    accession = {s.accession for s in scelti}
    nomi = {s.campione for s in scelti}
    regex = re.compile(config.io.accession_regex)

    def dell_assay(riga: dict[str, str]) -> bool:
        trovato = regex.search(riga.get(config.meta.accession_column, ""))
        return bool(trovato) and trovato.group(0) in accession

    conteggi = {
        NOME_ASSAY: _filtra_tabella(Path(config.io.assay_table), cartella / NOME_ASSAY, dell_assay),
        NOME_STUDIO: _filtra_tabella(
            Path(config.io.study_table), cartella / NOME_STUDIO,
            lambda r: r.get(config.meta.sample_id_column) in nomi,
        ),
    }
    if config.io.batch_table is not None:
        conteggi[NOME_LOTTI] = _filtra_tabella(
            Path(config.io.batch_table), cartella / NOME_LOTTI,
            lambda r: r.get(config.meta.batch_key_column) in accession,
        )
    for nome, tenute in conteggi.items():
        if tenute != len(scelti):
            raise SystemExit(f"{nome}: {tenute} righe per {len(scelti)} campioni")


# --------------------------------------------------------------------------- #
# Sottocampionamento                                                           #
# --------------------------------------------------------------------------- #


def _record(percorso: Path) -> Iterator[tuple[str, str, str, str]]:
    with gzip.open(percorso, "rt", encoding="ascii") as file:
        while intestazione := file.readline():
            yield intestazione, file.readline(), file.readline(), file.readline()


def quante(letture: int) -> int:
    """Quante letture conservare di un campione che ne ha ``letture``."""
    return min(letture, max(MINIMO, math.ceil(letture * FRAZIONE)))


def seme_del_file(accession: str) -> int:
    """Seme proprio di un file, derivato dal seme fisso e dall'accession."""
    impronta = hashlib.sha256(f"{SEME}:{accession}".encode()).digest()
    return int.from_bytes(impronta[:8], "big")


def sottocampiona(sorgente: Path, destinazione: Path, accession: str) -> tuple[int, int]:
    """Scrive le letture scelte, nell'ordine originale. Restituisce (tenute, minimo)."""
    lunghezze = [len(r[1].rstrip("\n")) for r in _record(sorgente)]
    minimo = min(lunghezze)
    obbligate = {i for i, lunghezza in enumerate(lunghezze) if lunghezza == minimo}
    restanti = [i for i in range(len(lunghezze)) if i not in obbligate]
    da_estrarre = max(0, quante(len(lunghezze)) - len(obbligate))
    scelte = obbligate | set(random.Random(seme_del_file(accession)).sample(restanti, da_estrarre))

    with open(destinazione, "wb") as grezzo, \
            gzip.GzipFile(filename="", mode="wb", fileobj=grezzo, mtime=0, compresslevel=9) as uscita:
        for indice, record in enumerate(_record(sorgente)):
            if indice in scelte:
                uscita.write("".join(record).encode("ascii"))
    return len(scelte), minimo


# --------------------------------------------------------------------------- #
# Comandi                                                                      #
# --------------------------------------------------------------------------- #


def _verifica(config: Config, scelti: list[Scelto]) -> dict[str, Path]:
    percorsi = {}
    for scelto in scelti:
        percorso = _file_di(config, scelto.accession)
        if md5(percorso) != scelto.md5:
            raise SystemExit(f"{percorso.name}: MD5 diverso da quello della selezione")
        percorsi[scelto.accession] = percorso
    return percorsi


def completo(config: Config, destinazione: Path) -> None:
    scelti = leggi_selezione()
    percorsi = _verifica(config, scelti)
    (destinazione / "fastq").mkdir(parents=True, exist_ok=True)
    for scelto in scelti:
        shutil.copyfile(percorsi[scelto.accession], destinazione / "fastq" / percorsi[scelto.accession].name)
    scrivi_metadati(config, scelti, destinazione / "metadati")
    print(f"{len(scelti)} campioni in {destinazione}")


def ridotto(config: Config, destinazione: Path) -> None:
    scelti = leggi_selezione()
    percorsi = _verifica(config, scelti)
    cartella = destinazione / "fastq"
    if cartella.exists():
        shutil.rmtree(cartella)
    cartella.mkdir(parents=True)

    righe = []
    for scelto in scelti:
        uscita = cartella / percorsi[scelto.accession].name
        tenute, minimo = sottocampiona(percorsi[scelto.accession], uscita, scelto.accession)
        if minimo != scelto.lunghezza_minima:
            raise SystemExit(f"{scelto.accession}: lunghezza minima {minimo}, attesa {scelto.lunghezza_minima}")
        righe.append([uscita.name, scelto.letture, tenute, minimo, md5(uscita)])

    scrivi_metadati(config, scelti, destinazione / "metadati")
    with open(destinazione / "manifesto.tsv", "w", encoding="utf-8", newline="") as file:
        scrittore = csv.writer(file, delimiter="\t", lineterminator="\n")
        scrittore.writerow(["file", "letture_originali", "letture", "lunghezza_minima", "md5"])
        scrittore.writerows(righe)

    totale = sum(p.stat().st_size for p in cartella.iterdir())
    print(f"{len(righe)} campioni, {sum(r[2] for r in righe)} letture, "
          f"{totale / 1e6:.2f} MB in {cartella}")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    comandi = parser.add_subparsers(dest="comando", required=True)
    for nome in ("completo", "ridotto"):
        sotto = comandi.add_parser(nome)
        sotto.add_argument("--config", required=True, type=Path,
                           help="configurazione che indica dati e metadati locali")
        sotto.add_argument("--destinazione", type=Path,
                           default=None if nome == "completo" else RIDOTTO)
    args = parser.parse_args(argv)
    if args.destinazione is None:
        parser.error("--destinazione e' obbligatoria per 'completo'")

    config = carica(args.config)
    (completo if args.comando == "completo" else ridotto)(config, args.destinazione)
    return 0


if __name__ == "__main__":
    sys.exit(main())
