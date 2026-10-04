"""Lettura delle tabelle dei metadati separate da tabulazioni.

Le tabelle ISA-Tab e il file di arricchimento del lotto sono compilati a mano:
spazi e virgolette ai margini dei valori compaiono in modo incostante, e una
colonna senza nome nasce da una tabulazione in piu' in fondo alla riga. Il
crosswalk di S0 e l'oggetto integrato di S10 le leggono allo stesso modo.
"""

from __future__ import annotations

import csv
from pathlib import Path

__all__ = ["leggi_tsv", "pulisci"]


def pulisci(valore: str | None) -> str:
    """Toglie spazi e virgolette: le tabelle ISA le usano in modo incostante."""
    return (valore or "").strip().strip('"').strip()


def leggi_tsv(percorso: Path) -> list[dict[str, str]]:
    """Le righe di una tabella separata da tabulazioni, con i valori ripuliti."""
    with open(percorso, encoding="utf-8", newline="") as file:
        return [
            {chiave: pulisci(valore) for chiave, valore in riga.items() if chiave}
            for riga in csv.DictReader(file, delimiter="\t")
        ]
