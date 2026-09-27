"""Infrastruttura di caricamento e configurazione per il sottoinsieme di prova OSD-734.

Inquadramento nel Piano Operativo:
    - **Settimana di riferimento**: **Settimana 11 (W11: Fase F4, Fase S1 e
      sottoinsieme di prova OSD-734)**.
    - **Moduli sorgente supportati**: ``src/amplicon16s/steps/s01_profile.py``,
      ``R/01_profile.R``, ``scripts/build_test_subset.py``.

La versione ridotta del sottoinsieme di prova vive in
``tests/fixtures/osd734/ridotto/`` ed è prodotta in modo deterministico da
``scripts/build_test_subset.py``. Questo modulo espone le funzioni di accesso
alla tabella di selezione, al manifesto delle letture sottocampionate e al
costruttore della configurazione validata che punta alla versione ridotta,
ereditando tutti i parametri scientifici da ``config/config.example.yaml``.

Il database tassonomico non fa parte della versione ridotta: pesa centinaia
di megabyte, e fino all'assegnazione tassonomica nessuna fase lo legge. Al
suo posto viene generato un file segnaposto con il proprio MD5, affinché la
Fase S0 (Gate G12) ne verifichi regolarmente l'integrità.
"""

from __future__ import annotations

import csv
import hashlib
from pathlib import Path
from typing import Any, Final

import yaml

from amplicon16s.config.schema import Config, valida

RADICE: Final = Path(__file__).resolve().parent / "fixtures" / "osd734"
SELEZIONE: Final = RADICE / "selezione.tsv"
RIDOTTO: Final = RADICE / "ridotto"
ESEMPIO: Final = Path(__file__).resolve().parents[1] / "config" / "config.example.yaml"


def selezione() -> list[dict[str, str]]:
    """Carica e restituisce le righe di ``tests/fixtures/osd734/selezione.tsv``.

    Ogni dizionario descrive uno dei 28 campioni selezionati dal dataset di
    riferimento OSD-734 (accession, run ENA, classe, piastra, corsa, numero di
    letture, lunghezza minima originale, checksum MD5 e motivazione).
    """
    with open(SELEZIONE, encoding="utf-8", newline="") as file:
        return list(csv.DictReader(file, delimiter="\t"))


def manifesto() -> list[dict[str, str]]:
    """Carica e restituisce le righe di ``tests/fixtures/osd734/ridotto/manifesto.tsv``.

    Ogni voce riporta il nome del file FASTQ sottocampionato, il numero di
    letture conservate, la lunghezza minima osservata e l'impronta MD5 del file
    compresso deterministico.
    """
    with open(RIDOTTO / "manifesto.tsv", encoding="utf-8", newline="") as file:
        return list(csv.DictReader(file, delimiter="\t"))


def dati_config(cartella: Path, **sovrascrivi: dict[str, Any]) -> dict[str, Any]:
    """I parametri d'esempio con i percorsi della versione ridotta."""
    dati = yaml.safe_load(ESEMPIO.read_text(encoding="utf-8"))
    riferimento = cartella / "riferimento.fa.gz"
    if not riferimento.exists():
        riferimento.write_bytes(b">segnaposto\nACGT\n")
    metadati = RIDOTTO / "metadati"
    dati["io"].update(
        fastq_dir=str(RIDOTTO / "fastq"),
        assay_table=str(metadati / "assay.txt"),
        study_table=str(metadati / "studio.txt"),
        batch_table=str(metadati / "lotti.tsv"),
        out_root=str(cartella / "out"),
    )
    dati["tax"].update(
        ref_fasta=str(riferimento),
        ref_md5=hashlib.md5(riferimento.read_bytes()).hexdigest(),
    )
    dati["run"]["threads"] = 2
    for gruppo, valori in sovrascrivi.items():
        dati.setdefault(gruppo, {}).update(valori)
    return dati


def config_ridotta(cartella: Path, **sovrascrivi: dict[str, Any]) -> Config:
    """Costruisce e valida un'istanza immutabile ``Config`` rivolta al subset ridotto.

    Crea un file FASTA segnaposto con relativo checksum MD5 nella directory
    temporanea indicata da ``cartella``, indirizza ``io.*`` ai 28 campioni
    sottocampionati di ``tests/fixtures/osd734/ridotto/``, applica eventuali
    sovrascritture per sezione e restituisce il modello Pydantic validato.
    """
    return valida(dati_config(cartella, **sovrascrivi))


def motivo_pacchetti_r_assenti(*pacchetti: str) -> str | None:
    """Perche' i pacchetti R indicati non sono utilizzabili, o ``None`` se lo sono."""
    import subprocess

    from amplicon16s.rbridge.runner import trova_rscript

    rscript = trova_rscript()
    if rscript is None:
        return "Rscript non disponibile"
    elenco = ", ".join(f"'{p}'" for p in pacchetti)
    prova = subprocess.run(
        [str(rscript), "--vanilla", "-e",
         f"quit(status = !all(vapply(c({elenco}), requireNamespace, logical(1), quietly = TRUE)))"],
        capture_output=True, check=False,
    )
    if prova.returncode != 0:
        return f"pacchetti R non installati fra: {', '.join(pacchetti)}"
    return None
