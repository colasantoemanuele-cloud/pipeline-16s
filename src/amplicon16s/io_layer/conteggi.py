"""Lettura delle tabelle del tracciamento delle letture.

Ogni fase che conta le letture per campione scrive una tabella
``letture_<passo>.tsv`` con le colonne ``campione``, ``passo`` e ``letture``
(``R/lib/letture.R``). Piu' fasi leggono le tabelle di altre fasi, per
riconoscere i campioni rimasti senza letture o per confrontare una profondita'
con una soglia: il lettore sta qui, fuori dai moduli delle fasi, perche' fa
parte del calcolo di ciascuna e ne entra nel sorgente
(:mod:`amplicon16s.runner.provenienza`).
"""

from __future__ import annotations

import csv
from pathlib import Path

from amplicon16s.io_layer.artifacts import AlberoOutput, Fase

__all__ = ["leggi_conteggi", "letture_in_ingresso_alla_tabella"]


def leggi_conteggi(percorso: Path) -> dict[str, int]:
    """Una tabella del tracciamento: campione -> letture."""
    with open(percorso, encoding="utf-8", newline="") as file:
        return {r["campione"]: int(float(r["letture"])) for r in csv.DictReader(file, delimiter="\t")}


def letture_in_ingresso_alla_tabella(albero: AlberoOutput, appaiata: bool) -> dict[str, int]:
    """Le letture con cui ogni campione arriva alla tabella delle sequenze: in
    modalita' single-end le letture filtrate, in modalita' appaiata le coppie
    fuse. Chi costruisce e filtra la tabella vi riconosce i campioni rimasti
    senza letture, che nella tabella non compaiono.
    """
    if appaiata:
        return leggi_conteggi(albero.cartella(Fase.PAIR_MERGE) / "letture_fuse.tsv")
    return leggi_conteggi(albero.cartella(Fase.FILTERED) / "letture_filtrate.tsv")
