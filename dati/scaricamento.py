"""Scarico verificato e riprendibile dei file pubblici dei dataset e del riferimento.

Usato dagli script ``scarica_letture.py`` dei dataset e da
``riferimento/scarica_riferimento.py``; solo libreria
standard, perche' gli script devono girare prima che l'ambiente della pipeline
esista.

**Ogni file ha un checksum dichiarato dalla fonte** (MD5 di ENA o di Zenodo):
un file e' presente solo se la dimensione e l'MD5 corrispondono. Un file con il
nome giusto ma il contenuto sbagliato non viene sovrascritto in silenzio: e'
messo da parte con il suffisso ``.md5_errato``, perche' va guardato, e il file
viene scaricato di nuovo.

**Lo scarico riprende da dove si e' interrotto.** Il trasferimento avviene su
``<nome>.parziale``; a un nuovo avvio si chiede al server soltanto il resto
(intestazione HTTP ``Range``). Se il server non lo concede, il file riparte da
capo. Il nome definitivo si assegna solo dopo la verifica dell'MD5: un
trasferimento interrotto non lascia mai un file incompleto con il nome giusto.

**Una richiesta alla volta, con una pausa fra l'una e l'altra**: ENA e Zenodo
sono archivi pubblici, e aprire molte connessioni in parallelo e' un modo di
abusarne. Un rifiuto temporaneo (per esempio un 403 o un 429 di limitazione
della frequenza) si ritenta dopo una pausa crescente.
"""

from __future__ import annotations

import hashlib
import time
import urllib.error
import urllib.request
from dataclasses import dataclass
from pathlib import Path
from typing import Final

__all__ = ["Atteso", "ErroreScarico", "md5_file", "metti_da_parte", "scarica", "valido"]

#: Dimensione dei blocchi letti e scritti: memoria costante anche sui file grandi.
BLOCCO: Final = 1024 * 1024
#: Il tempo massimo di attesa di una risposta, in secondi.
TEMPO_MASSIMO: Final = 60


class ErroreScarico(RuntimeError):
    """Un file non e' stato ottenuto integro dopo i tentativi ammessi."""


@dataclass(frozen=True)
class Atteso:
    """Un file da ottenere: da dove, con quale MD5 e quanti byte."""

    url: str
    md5: str
    byte: int


def md5_file(percorso: Path) -> str:
    """L'MD5 di un file, letto a blocchi."""
    impronta = hashlib.md5()
    with open(percorso, "rb") as file:
        while blocco := file.read(BLOCCO):
            impronta.update(blocco)
    return impronta.hexdigest()


def valido(percorso: Path, atteso: Atteso) -> bool:
    """Vero se il file esiste con la dimensione e l'MD5 attesi.

    La dimensione si controlla per prima perche' non costa nulla: un file
    interrotto ha quasi sempre una lunghezza diversa.
    """
    if not percorso.is_file() or percorso.stat().st_size != atteso.byte:
        return False
    return md5_file(percorso) == atteso.md5


def metti_da_parte(percorso: Path) -> Path:
    """Rinomina un file non conforme con il suffisso ``.md5_errato``."""
    da_parte = percorso.with_name(percorso.name + ".md5_errato")
    percorso.replace(da_parte)
    return da_parte


def _trasferisci(url: str, parziale: Path) -> None:
    """Scrive in ``parziale`` il contenuto di ``url``, riprendendo dal byte a cui
    il file parziale e' arrivato se il server lo concede.
    """
    gia = parziale.stat().st_size if parziale.is_file() else 0
    richiesta = urllib.request.Request(url, headers={"User-Agent": "amplicon16s-dati"})
    if gia:
        richiesta.add_header("Range", f"bytes={gia}-")
    with urllib.request.urlopen(richiesta, timeout=TEMPO_MASSIMO) as risposta:
        # 206: il server ha mandato soltanto il resto; altrimenti il file intero.
        ripresa = gia and getattr(risposta, "status", None) == 206
        with open(parziale, "ab" if ripresa else "wb") as file:
            while blocco := risposta.read(BLOCCO):
                file.write(blocco)


def scarica(
    atteso: Atteso,
    destinazione: Path,
    *,
    tentativi: int = 4,
    pausa: float = 5.0,
) -> None:
    """Porta ``destinazione`` allo stato atteso, scaricando e verificando.

    Solleva :class:`ErroreScarico` se dopo ``tentativi`` il file non e' integro.
    Un MD5 sbagliato a trasferimento concluso non si ritenta riprendendo, ma da
    capo: il file parziale e' messo da parte.
    """
    destinazione.parent.mkdir(parents=True, exist_ok=True)
    parziale = destinazione.with_name(destinazione.name + ".parziale")
    ultimo = "nessun tentativo"
    for tentativo in range(1, tentativi + 1):
        try:
            _trasferisci(atteso.url, parziale)
        except (urllib.error.URLError, OSError) as guasto:
            ultimo = f"trasferimento interrotto: {guasto}"
            time.sleep(pausa * tentativo)
            continue
        if parziale.stat().st_size < atteso.byte:
            ultimo = f"{parziale.stat().st_size} byte su {atteso.byte}"
            continue
        if valido(parziale, atteso):
            parziale.replace(destinazione)
            return
        ultimo = f"MD5 diverso da {atteso.md5}"
        metti_da_parte(parziale)
    raise ErroreScarico(f"{destinazione.name}: {ultimo}")
