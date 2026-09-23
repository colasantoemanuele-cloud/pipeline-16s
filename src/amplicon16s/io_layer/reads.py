"""Ispezione delle prime letture di ciascun file, per i gate di validazione.

I gate che guardano le sequenze — struttura del file, lunghezze, presenza del
primer — devono costare minuti, non ore. Nessuno di essi legge un file intero:
si fermano alle prime ``qc.head_reads`` letture. L'integrità completa dei file
è già garantita altrove, dai checksum del manifesto.

La lettura avviene **una volta sola per file**: tutte le statistiche che
servono ai vari gate si ricavano dalla stessa passata, perché riaprire e
decomprimere lo stesso file tre volte triplicherebbe il costo senza aggiungere
nulla.

La scansione è **sequenziale**. Distribuirla su più processi è stato provato e
scartato: sull'intero dataset di riferimento, 960 file, la passata sequenziale
costa una ventina di secondi, e il vincolo da rispettare è che la validazione
duri minuti e non ore. Il parallelismo faceva risparmiare pochi secondi in
cambio di una fragilità vera — con ``fork`` il pool può bloccarsi quando il
processo genitore ha sostituito i flussi standard, e con ``spawn`` ogni
chiamante sarebbe costretto a proteggere il proprio modulo principale con
``if __name__ == "__main__"``. Non è un prezzo che una libreria debba far
pagare a chi la usa.

**Limite dichiarato.** Le lunghezze osservate sono quelle delle prime letture,
non di tutto il file. Un file le cui letture più corte stiano oltre quelle
esaminate darebbe un minimo sovrastimato. È il prezzo del vincolo di costo, e
va tenuto presente da chi usa :attr:`StatisticheFile.lunghezza_minima`.
"""

from __future__ import annotations

import gzip
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Final

__all__ = [
    "StatisticheFile",
    "espandi_iupac",
    "scansiona",
    "scansiona_file",
]

#: Corrispondenza fra codici IUPAC e classi di caratteri.
_IUPAC: Final[dict[str, str]] = {
    "A": "A", "C": "C", "G": "G", "T": "T",
    "R": "[AG]", "Y": "[CT]", "S": "[GC]", "W": "[AT]",
    "K": "[GT]", "M": "[AC]", "B": "[CGT]", "D": "[AGT]",
    "H": "[ACT]", "V": "[ACG]", "N": ".",
}


def espandi_iupac(sequenza: str) -> str:
    """Traduce una sequenza con codici degenerati in un'espressione regolare.

    ``GTGYCAGC`` diventa ``GTG[CT]CAGC``: cercare la sequenza alla lettera
    mancherebbe le letture che nelle posizioni degenerate portano l'altra
    base ammessa, e il primer risulterebbe assente anche quando c'è.
    """
    try:
        return "".join(_IUPAC[carattere] for carattere in sequenza.upper())
    except KeyError as guasto:
        raise ValueError(f"codice IUPAC sconosciuto: {guasto.args[0]!r}") from None


@dataclass(frozen=True)
class StatisticheFile:
    """Ciò che si è potuto misurare sulle prime letture di un file."""

    nome: str
    letture_esaminate: int
    lunghezza_minima: int | None
    lunghezza_massima: int | None
    con_primer: int
    con_motivo: int
    #: Descrizione del guasto strutturale, se il file non è leggibile come
    #: FASTQ compresso. ``None`` se la struttura è valida.
    errore: str | None = None
    #: Vero se il file è finito prima di raggiungere il numero richiesto: non
    #: è un problema, ma dice che le statistiche coprono tutto il file.
    esaurito: bool = False

    @property
    def valido(self) -> bool:
        return self.errore is None

    @property
    def frazione_primer(self) -> float:
        return self.con_primer / self.letture_esaminate if self.letture_esaminate else 0.0

    @property
    def frazione_motivo(self) -> float:
        return self.con_motivo / self.letture_esaminate if self.letture_esaminate else 0.0


def scansiona_file(
    percorso: Path | str,
    head_reads: int,
    primer: str,
    motivo: str,
) -> StatisticheFile:
    """Legge le prime ``head_reads`` letture e ne ricava tutte le statistiche.

    ``primer`` e ``motivo`` sono espressioni regolari già pronte, passate come
    stringhe gia' compilate a ogni chiamata.
    Entrambe si cercano **ancorate a inizio lettura**: il primer starebbe in
    testa se non fosse stato tolto, e il motivo conservato è quello che apre la
    regione amplificata.
    """
    percorso = Path(percorso)
    espressione_primer = re.compile(primer)
    espressione_motivo = re.compile(motivo)

    letture = con_primer = con_motivo = 0
    minima: int | None = None
    massima: int | None = None
    esaurito = True

    try:
        with gzip.open(percorso, "rt", encoding="utf-8", errors="replace") as file:
            while letture < head_reads:
                intestazione = file.readline()
                if not intestazione:
                    break  # fine del file su un confine di record: corretto
                sequenza = file.readline()
                separatore = file.readline()
                qualita = file.readline()

                if not (sequenza and separatore and qualita):
                    return StatisticheFile(
                        percorso.name, letture, minima, massima, con_primer, con_motivo,
                        errore=f"record troncato dopo {letture} letture: "
                               f"il file non ha quattro righe per record",
                    )
                if not intestazione.startswith("@"):
                    return StatisticheFile(
                        percorso.name, letture, minima, massima, con_primer, con_motivo,
                        errore=f"alla lettura {letture + 1} l'intestazione non comincia "
                               f"con '@': {intestazione[:40]!r}",
                    )
                if not separatore.startswith("+"):
                    return StatisticheFile(
                        percorso.name, letture, minima, massima, con_primer, con_motivo,
                        errore=f"alla lettura {letture + 1} la terza riga non comincia "
                               f"con '+': {separatore[:40]!r}",
                    )

                sequenza = sequenza.rstrip("\n")
                if len(sequenza) != len(qualita.rstrip("\n")):
                    return StatisticheFile(
                        percorso.name, letture, minima, massima, con_primer, con_motivo,
                        errore=f"alla lettura {letture + 1} sequenza e qualita' hanno "
                               f"lunghezze diverse",
                    )

                letture += 1
                lunghezza = len(sequenza)
                minima = lunghezza if minima is None else min(minima, lunghezza)
                massima = lunghezza if massima is None else max(massima, lunghezza)
                if espressione_primer.match(sequenza):
                    con_primer += 1
                if espressione_motivo.match(sequenza):
                    con_motivo += 1
            else:
                esaurito = False
    except (OSError, EOFError, gzip.BadGzipFile) as guasto:
        return StatisticheFile(
            percorso.name, letture, minima, massima, con_primer, con_motivo,
            errore=f"archivio non leggibile: {guasto}",
        )

    if letture == 0:
        return StatisticheFile(
            percorso.name, 0, None, None, 0, 0,
            errore="il file non contiene alcuna lettura",
        )

    return StatisticheFile(
        percorso.name, letture, minima, massima, con_primer, con_motivo,
        esaurito=esaurito,
    )


def scansiona(
    percorsi: dict[str, Path],
    head_reads: int,
    primer: str,
    motivo: str,
) -> dict[str, StatisticheFile]:
    """Scansiona i file indicati e restituisce le statistiche di ciascuno."""
    return {
        chiave: scansiona_file(percorso, head_reads, primer, motivo)
        for chiave, percorso in percorsi.items()
    }
