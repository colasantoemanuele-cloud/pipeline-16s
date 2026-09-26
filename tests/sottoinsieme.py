"""La versione ridotta del sottoinsieme di prova, pronta per i test.

Vive in ``tests/fixtures/osd734/ridotto/`` ed è prodotta da
``scripts/build_test_subset.py``. Qui si costruisce la configurazione che la
indica: i parametri sono quelli di ``config/config.example.yaml``, cioè del
dataset di riferimento, con i percorsi rivolti alla versione ridotta.

Il database tassonomico non fa parte della versione ridotta: pesa centinaia
di megabyte, e fino all'assegnazione tassonomica nessuna fase lo legge. Al
suo posto c'e' un file segnaposto con il proprio MD5, perche' S0 ne verifica
l'integrita'.
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
    with open(SELEZIONE, encoding="utf-8", newline="") as file:
        return list(csv.DictReader(file, delimiter="\t"))


def manifesto() -> list[dict[str, str]]:
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
    return valida(dati_config(cartella, **sovrascrivi))
