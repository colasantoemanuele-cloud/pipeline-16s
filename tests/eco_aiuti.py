"""Aiuti condivisi dei test delle analisi ecologiche.

Costruiscono oggetti phyloseq di prova da tabelle scritte nel test, una
configurazione valida da cui ogni test parte, e leggono le tabelle prodotte.
Gli oggetti nascono in R (``tests/fixtures/eco/costruisci_oggetto.R``) perche'
il formato e' quello che la pipeline consegna, e Python non lo scrive.
"""

from __future__ import annotations

import copy
import functools
import hashlib
import json
import os
import subprocess
from pathlib import Path
from typing import Any, Final

import pytest

from amplicon16s.rbridge.runner import trova_rscript

from sottoinsieme import motivo_pacchetti_r_assenti

RADICE: Final = Path(__file__).resolve().parents[1]
FIXTURE: Final = Path(__file__).resolve().parent / "fixtures" / "eco"
COSTRUTTORE: Final = FIXTURE / "costruisci_oggetto.R"

#: I ranghi degli oggetti di prova.
RANGHI: Final = ("Regno", "Famiglia", "Genere")


@functools.cache
def _sonda() -> str | None:
    """Perche' l'ambiente R delle analisi non c'e', o ``None``: la sonda parte
    al primo uso e una volta sola.
    """
    return motivo_pacchetti_r_assenti("phyloseq", "vegan", "ape", "permute", "jsonlite")


@pytest.fixture
def eco_r():
    """Richiede R con phyloseq, vegan e ape: salta senza, ma fallisce dove
    l'ambiente della pipeline e' richiesto.
    """
    if _sonda() is not None:
        if os.environ.get("AMPLICON16S_RICHIEDI_BIOC") == "1":
            pytest.fail(f"Bioconductor e' richiesto in questo ambiente: {_sonda()}")
        pytest.skip(_sonda())


def configurazione(**sovrascrivi: dict[str, Any]) -> dict[str, Any]:
    """Una configurazione valida per gli oggetti di prova: la variabile
    biologica e' la colonna ``gruppo``. Ogni gruppo passato sostituisce o
    aggiunge le chiavi indicate; un valore ``...`` toglie la chiave.
    """
    dati: dict[str, Any] = {
        "design": {"variable": "gruppo"},
        "run": {"seed": 7},
        "comp": {"rank": "Genere", "top_n": 3},
        "beta": {"distances": ["bray", "jaccard", "aitchison"], "clr_pseudocount": 0.5},
        "stat": {"min_group_size": 2},
    }
    dati = copy.deepcopy(dati)
    for gruppo, valori in sovrascrivi.items():
        dati.setdefault(gruppo, {})
        for chiave, valore in valori.items():
            if valore is ...:
                dati[gruppo].pop(chiave, None)
            else:
                dati[gruppo][chiave] = valore
    return dati


def tassonomia_semplice(n: int) -> dict[str, list[str | None]]:
    """Una tassonomia in cui ogni variante ha il proprio genere."""
    return {
        "Regno": ["Batteri"] * n,
        "Famiglia": [f"Fam{i // 2 + 1}" for i in range(n)],
        "Genere": [f"Gen{i + 1}" for i in range(n)],
    }


def costruisci(
    cartella: Path,
    conteggi: list[list[int]],
    metadati: dict[str, list[str | None]],
    tassonomia: dict[str, list[str | None]] | None = None,
    albero: str | None = None,
    nome: str = "oggetto.rds",
) -> Path:
    """Scrive in ``cartella/dati`` un oggetto phyloseq con i conteggi dati (una
    riga per variante, una colonna per campione) e ne restituisce il percorso.
    I campioni si chiamano ``c01``..., le varianti ``v01``...
    """
    n_varianti, n_campioni = len(conteggi), len(conteggi[0])
    descrizione = {
        "conteggi": conteggi,
        "varianti": [f"v{i + 1:02d}" for i in range(n_varianti)],
        "campioni": [f"c{i + 1:02d}" for i in range(n_campioni)],
        "tassonomia": tassonomia or tassonomia_semplice(n_varianti),
        "metadati": metadati,
        "albero": albero,
    }
    dati = cartella / "dati"
    dati.mkdir(exist_ok=True)
    richiesta = dati / (nome + ".json")
    richiesta.write_text(json.dumps(descrizione), encoding="utf-8")
    oggetto = dati / nome
    esegui_rscript(COSTRUTTORE, richiesta, oggetto)
    return oggetto


def esegui_rscript(script: Path, *argomenti: Any) -> None:
    """Esegue uno script R di prova e fallisce con la sua uscita di errore."""
    processo = subprocess.run(
        [str(trova_rscript()), "--vanilla", str(script), *map(str, argomenti)],
        capture_output=True, text=True, check=False,
    )
    assert processo.returncode == 0, processo.stderr[-3000:]


def leggi_tsv(percorso: Path) -> tuple[list[str], list[dict[str, str]]]:
    """Le righe di commento (senza il prefisso) e le righe di una tabella."""
    righe = percorso.read_text(encoding="utf-8").splitlines()
    commenti = [r[2:] for r in righe if r.startswith("# ")]
    corpo = [r.split("\t") for r in righe if not r.startswith("# ")]
    nomi = corpo[0]
    return commenti, [dict(zip(nomi, r, strict=True)) for r in corpo[1:]]


def impronte(cartella: Path) -> dict[str, str]:
    """SHA-256 di ogni file della cartella, per nome relativo."""
    return {
        p.relative_to(cartella).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest()
        for p in sorted(cartella.rglob("*")) if p.is_file()
    }
