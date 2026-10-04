"""Lettura dell'inventario dei campioni scritto da S0.

S0 scrive l'inventario in ``01_input_validation/crosswalk.tsv``; le fasi a
valle non lo ricostruiscono dai metadati ma lo rileggono da li', attraverso il
contesto che l'esecuzione passa loro (:mod:`amplicon16s.runner.project`). Il
lettore sta in un modulo proprio perche' decide su quali campioni, con quali
classi, piastre e corse, le fasi calcolano: senza essere importato da quelle
fasi, entra comunque nel sorgente di ogni fase che dipende da S0
(:mod:`amplicon16s.runner.provenienza`).
"""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Final

from amplicon16s.config.schema import Config
from amplicon16s.io_layer.artifacts import AlberoOutput, Fase
from amplicon16s.metadata.models import Campione, ClasseCampione, Inventario

__all__ = ["NOME_CROSSWALK", "leggi_inventario"]

#: L'inventario scritto da S0, una riga per campione.
NOME_CROSSWALK: Final = "crosswalk.tsv"


def leggi_inventario(config: Config) -> Inventario:
    """Rilegge l'inventario dal crosswalk scritto da S0.

    È l'inventario su cui le fasi a valle sono state calcolate: rileggerlo
    dall'artefatto, invece di ricostruirlo dai metadati, fa sì che una ripresa
    usi esattamente quello.
    """
    percorso = AlberoOutput(config.io.out_root).cartella(Fase.INPUT_VALIDATION) / NOME_CROSSWALK
    campioni = []
    with open(percorso, encoding="utf-8", newline="") as file:
        for riga in csv.DictReader(file, delimiter="\t"):
            campioni.append(
                Campione(
                    accession=riga["accession"],
                    nome=riga["sample_name"],
                    classe=ClasseCampione(riga["classe"]),
                    materiale=riga["materiale"],
                    file=Path(config.io.fastq_dir) / riga["file"] if riga["file"] else None,
                    posizione=riga["posizione"] or None,
                    modulo=riga["modulo"] or None,
                    piastra=riga["piastra"] or None,
                    corsa=riga["corsa"] or None,
                )
            )
    return Inventario(tuple(campioni))
