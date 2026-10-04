"""Il tracciamento delle letture, ricomposto dalle fasi che lo producono.

Ogni fase registra quante letture restano a ciascun campione dopo il proprio
passo, in un file della propria cartella (``letture_<passo>.tsv``, scritto con
``traccia_letture`` di ``R/lib/letture.R``). Nessuna fase aggiunge colonne a
un file di un'altra: lo modificherebbe, ne romperebbe il checksum, e la fase
che l'ha prodotto risulterebbe da rifare.

La tabella completa non esiste come file: si ricompone qui leggendo, nell'ordine
del grafo, i file di tracciamento elencati nei manifesti delle fasi concluse.
Un file di una fase non conclusa non entra, perche' potrebbe appartenere a un
calcolo interrotto o superato.
"""

from __future__ import annotations

import csv
from collections.abc import Iterable
from dataclasses import dataclass
from typing import Final

from amplicon16s.io_layer.artifacts import AlberoOutput
from amplicon16s.runner.project import ProjectRun, StatoPasso, Valutazione
from amplicon16s.steps.base import PipelineStep

__all__ = ["PREFISSO", "Tracciamento", "componi", "ricomponi"]

#: I file di tracciamento si riconoscono dal nome.
PREFISSO: Final = "letture_"


@dataclass(frozen=True)
class Tracciamento:
    """Letture per campione, un passo dopo l'altro."""

    #: I passi del tracciamento, nell'ordine in cui avvengono.
    passi: tuple[str, ...]
    #: campione -> passo -> letture.
    letture: dict[str, dict[str, int]]
    #: Da quale fase viene ciascun passo.
    origine: dict[str, str]

    def totali(self) -> dict[str, int]:
        """Il totale delle letture di tutti i campioni per ciascun passo."""
        return {
            passo: sum(v.get(passo, 0) for v in self.letture.values()) for passo in self.passi
        }


def ricomponi(run: ProjectRun, valutazione: Valutazione | None = None) -> Tracciamento:
    """La tabella completa, dai file delle fasi concluse. Non scrive nulla."""
    valutazione = valutazione or run.valuta()
    concluse = [
        run.passi[fase]
        for fase, situazione in valutazione.situazioni.items()
        if situazione.stato is StatoPasso.COMPLETATA and fase in run.passi
    ]
    return componi(run.albero, concluse)


def componi(albero: AlberoOutput, concluse: Iterable[PipelineStep]) -> Tracciamento:
    """La tabella completa dalle fasi indicate, concluse e nell'ordine del grafo.

    Si leggono solo i passi che ciascuna fase dichiara in
    ``passi_tracciamento``, nell'ordine dichiarato: il nome di un file non
    basta a dire che e' un tracciamento, ne' in che ordine vengono i passi.
    Chi decide quali fasi sono concluse e' il chiamante: la valutazione dello
    stato per :func:`ricomponi`, i manifesti su disco per il report.
    """
    passi: list[str] = []
    origine: dict[str, str] = {}
    letture: dict[str, dict[str, int]] = {}

    for fase in concluse:
        manifesto = albero.manifesto_passo(fase.passo, fase.cartella)
        assert manifesto is not None
        for passo in fase.passi_tracciamento:
            nome = f"{PREFISSO}{passo}.tsv"
            if nome not in manifesto.nomi:
                raise ValueError(f"{fase.passo} dichiara il passo {passo} ma non ha {nome}")
            if passo in origine:
                raise ValueError(
                    f"il passo {passo} e' registrato sia da {origine[passo]} sia da {fase.passo}"
                )
            passi.append(passo)
            origine[passo] = str(fase.passo)
            with open(albero.cartella(fase.cartella) / nome, encoding="utf-8", newline="") as file:
                for riga in csv.DictReader(file, delimiter="\t"):
                    if riga["passo"] != passo:
                        raise ValueError(f"{nome} contiene il passo {riga['passo']}")
                    letture.setdefault(riga["campione"], {})[passo] = int(float(riga["letture"]))

    return Tracciamento(tuple(passi), letture, origine)
