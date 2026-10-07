"""Lettura delle tabelle dei metadati separate da tabulazioni.

Le tabelle ISA-Tab e il file di arricchimento del lotto sono compilati a mano:
spazi e virgolette ai margini dei valori compaiono in modo incostante, una
colonna senza nome nasce da una tabulazione in piu' in fondo alla riga, e un
foglio di calcolo puo' salvare il file con un segno d'ordine dei byte (BOM) in
testa. Tutta la pipeline le legge da qui (i gate di S0, il crosswalk, l'oggetto
integrato di S10): intestazioni e valori sono ripuliti in un solo modo, e una
colonna dichiarata in configurazione si confronta sempre con lo stesso nome.
"""

from __future__ import annotations

import csv
from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from amplicon16s.config.schema import Config

__all__ = ["intestazione", "leggi_tsv", "pulisci", "tabella_di_studio"]

#: UTF-8 con o senza BOM: con "utf-8-sig" il segno, se c'e', non finisce nel
#: nome della prima colonna.
_CODIFICA = "utf-8-sig"


def pulisci(valore: str | None) -> str:
    """Toglie spazi e virgolette: le tabelle ISA le usano in modo incostante."""
    return (valore or "").strip().strip('"').strip()


def intestazione(percorso: Path) -> list[str]:
    """I nomi delle colonne di una tabella, ripuliti come quelli che
    :func:`leggi_tsv` usa per chiave; elenco vuoto se il file non ha righe.
    """
    with open(percorso, encoding=_CODIFICA, newline="") as file:
        return [pulisci(c) for c in next(csv.reader(file, delimiter="\t"), [])]


def leggi_tsv(percorso: Path) -> list[dict[str, str]]:
    """Le righe di una tabella separata da tabulazioni, con nomi delle colonne e
    valori ripuliti. Le colonne senza nome si ignorano.
    """
    with open(percorso, encoding=_CODIFICA, newline="") as file:
        return [
            {nome: pulisci(valore) for chiave, valore in riga.items()
             if chiave and (nome := pulisci(chiave))}
            for riga in csv.DictReader(file, delimiter="\t")
        ]


def tabella_di_studio(config: Config) -> tuple[Path, str]:
    """La tabella da cui si leggono classe e variabili dei campioni, e la sua
    colonna con il nome del campione.

    E' ``io.study_table`` se indicata; altrimenti la tabella di assay, che
    allora porta anche quelle colonne.
    """
    if config.io.study_table is not None:
        return Path(config.io.study_table), config.meta.colonna_id_studio
    return Path(config.io.assay_table), config.meta.sample_id_column
