"""Ispezione delle prime letture di ciascun file, per i gate di validazione.

I gate che guardano le sequenze (struttura del file, lunghezze, presenza del
primer) devono costare minuti, non ore. Nessuno di essi legge un file intero:
si fermano alle prime ``qc.head_reads`` letture. L'integrità completa dei file
è già garantita altrove, dai checksum del manifesto.

La lettura avviene **una volta sola per file**: tutte le statistiche che
servono ai vari gate si ricavano dalla stessa passata, perché riaprire e
decomprimere lo stesso file tre volte triplicherebbe il costo senza aggiungere
nulla.

La scansione è **sequenziale**. Distribuirla su più processi è stato provato e
scartato: la passata sequenziale su un dataset di centinaia di file costa
secondi, e il vincolo da rispettare è che la validazione duri minuti e non ore. Il parallelismo faceva risparmiare pochi secondi in
cambio di una fragilità vera: con ``fork`` il pool può bloccarsi quando il
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
from typing import IO, Final

__all__ = [
    "StatisticheFile",
    "apri_fastq",
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


#: I primi due byte di un archivio gzip.
_MAGIA_GZIP: Final = b"\x1f\x8b"


def apri_fastq(percorso: Path | str) -> IO[str]:
    """Apre in lettura un FASTQ, compresso con gzip o no.

    Il formato si riconosce dai primi due byte, non dall'estensione: un file
    ``.fastq.gz`` non compresso, o un ``.fastq`` compresso, si leggono per
    quello che sono.
    """
    with open(percorso, "rb") as file:
        compresso = file.read(2) == _MAGIA_GZIP
    if compresso:
        return gzip.open(percorso, "rt", encoding="utf-8", errors="replace")
    return open(percorso, encoding="utf-8", errors="replace", newline="\n")


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
    #: FASTQ. ``None`` se la struttura è valida.
    errore: str | None = None
    #: Vero se il file è finito prima di raggiungere il numero richiesto: non
    #: è un problema, ma dice che le statistiche coprono tutto il file.
    esaurito: bool = False
    #: Letture esaminate piu' corte della lunghezza richiesta alla scansione
    #: (``filter.truncLen``): il filtro le scarterebbe.
    piu_corte: int = 0

    @property
    def frazione_corte(self) -> float:
        """La frazione delle letture esaminate piu' corte della lunghezza richiesta."""
        return self.piu_corte / self.letture_esaminate if self.letture_esaminate else 0.0

    @property
    def valido(self) -> bool:
        """Vero se il file si è letto senza errori."""
        return self.errore is None

    @property
    def frazione_primer(self) -> float:
        """La frazione delle letture esaminate che iniziano con il primer; 0 se nessuna
        è stata esaminata.
        """
        return self.con_primer / self.letture_esaminate if self.letture_esaminate else 0.0

    @property
    def frazione_motivo(self) -> float:
        """La frazione delle letture esaminate che contengono il motivo conservato; 0 se
        nessuna è stata esaminata.
        """
        return self.con_motivo / self.letture_esaminate if self.letture_esaminate else 0.0


def scansiona_file(
    percorso: Path | str,
    head_reads: int,
    primer: str | None,
    motivo: str | None,
    inizio_motivo: int = 0,
    lunghezza_richiesta: int = 0,
) -> StatisticheFile:
    """Legge le prime ``head_reads`` letture e ne ricava tutte le statistiche.

    ``lunghezza_richiesta`` e' la lunghezza sotto la quale una lettura si conta
    fra le piu' corte (``filter.truncLen``); zero non ne conta nessuna.

    ``primer`` e ``motivo`` sono espressioni regolari. Il primer si cerca
    **ancorato a inizio lettura**: starebbe in testa se non fosse stato tolto.
    Il motivo conservato, che apre la regione amplificata, si cerca ancorato
    alla posizione ``inizio_motivo``: zero se le letture iniziano dalla
    regione, la lunghezza di cio' che la precede (``filter.trimLeft``) se
    portano ancora il primer che il filtro togliera'. Con ``None`` la ricerca
    corrispondente non si fa e il conteggio resta zero.
    """
    percorso = Path(percorso)
    espressione_primer = re.compile(primer) if primer is not None else None
    espressione_motivo = re.compile(motivo) if motivo is not None else None

    letture = con_primer = con_motivo = corte = 0
    minima: int | None = None
    massima: int | None = None
    esaurito = True

    try:
        with apri_fastq(percorso) as file:
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

                sequenza = sequenza.rstrip("\r\n")
                if len(sequenza) != len(qualita.rstrip("\r\n")):
                    return StatisticheFile(
                        percorso.name, letture, minima, massima, con_primer, con_motivo,
                        errore=f"alla lettura {letture + 1} sequenza e qualita' hanno "
                               f"lunghezze diverse",
                    )

                letture += 1
                lunghezza = len(sequenza)
                minima = lunghezza if minima is None else min(minima, lunghezza)
                massima = lunghezza if massima is None else max(massima, lunghezza)
                if lunghezza < lunghezza_richiesta:
                    corte += 1
                if espressione_primer is not None and espressione_primer.match(sequenza):
                    con_primer += 1
                # Sulla lettura tagliata, non dalla posizione del taglio: un
                # motivo ancorato con ^ deve valere all'inizio di cio' che resta.
                if espressione_motivo is not None and espressione_motivo.match(
                    sequenza[inizio_motivo:]
                ):
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
        esaurito=esaurito, piu_corte=corte,
    )


def scansiona(
    percorsi: dict[str, Path],
    head_reads: int,
    primer: str | None,
    motivo: str | None,
    inizio_motivo: int = 0,
    lunghezza_richiesta: int = 0,
) -> dict[str, StatisticheFile]:
    """Scansiona i file indicati e restituisce le statistiche di ciascuno."""
    return {
        chiave: scansiona_file(percorso, head_reads, primer, motivo, inizio_motivo,
                               lunghezza_richiesta)
        for chiave, percorso in percorsi.items()
    }
