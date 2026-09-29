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
from dataclasses import dataclass
from typing import Final

from amplicon16s.runner.project import ProjectRun, StatoPasso, Valutazione

__all__ = ["PREFISSO", "Tracciamento", "ricomponi"]

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
    """La tabella completa, dai file delle fasi concluse. Non scrive nulla.

    Si leggono solo i passi che ciascuna fase dichiara in
    ``passi_tracciamento``, nell'ordine dichiarato: il nome di un file non
    basta a dire che e' un tracciamento, ne' in che ordine vengono i passi.
    """
    valutazione = valutazione or run.valuta()
    passi: list[str] = []
    origine: dict[str, str] = {}
    letture: dict[str, dict[str, int]] = {}

    for fase, situazione in valutazione.situazioni.items():
        if situazione.stato is not StatoPasso.COMPLETATA or fase not in run.passi:
            continue
        cartella = run.grafo.nodo(fase).cartella
        manifesto = run.albero.manifesto_passo(fase, cartella)
        assert manifesto is not None
        for passo in run.passi[fase].passi_tracciamento:
            nome = f"{PREFISSO}{passo}.tsv"
            if nome not in manifesto.nomi:
                raise ValueError(f"{fase} dichiara il passo {passo} ma non ha {nome}")
            if passo in origine:
                raise ValueError(f"il passo {passo} e' registrato sia da {origine[passo]} sia da {fase}")
            passi.append(passo)
            origine[passo] = str(fase)
            with open(run.albero.cartella(cartella) / nome, encoding="utf-8", newline="") as file:
                for riga in csv.DictReader(file, delimiter="\t"):
                    if riga["passo"] != passo:
                        raise ValueError(f"{nome} contiene il passo {riga['passo']}")
                    letture.setdefault(riga["campione"], {})[passo] = int(float(riga["letture"]))

    return Tracciamento(tuple(passi), letture, origine)
