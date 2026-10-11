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

**Poche connessioni alla volta, con una pausa fra una richiesta e l'altra**:
ENA e Zenodo sono archivi pubblici, e aprire molte connessioni in parallelo e'
un modo di abusarne. :func:`scarica_in_parallelo` ne tiene aperte al piu'
:data:`CONNESSIONI` (modificabile da chi lancia lo script). Un rifiuto
temporaneo (per esempio un 403 o un 429 di limitazione della frequenza) si
ritenta dopo una pausa crescente.

**Un guasto di rete e' un esito, non un'eccezione**: connessione caduta,
risposta interrotta a meta', nome non risolto finiscono in
:class:`ErroreScarico`, con il file e che cosa fare.
"""

from __future__ import annotations

import hashlib
import http.client
import time
import urllib.error
import urllib.request
from collections.abc import Iterator, Sequence
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Final

__all__ = [
    "CONNESSIONI", "GUASTI_DI_RETE", "Atteso", "ErroreScarico", "md5_file", "metti_da_parte",
    "scarica", "scarica_in_parallelo", "valido",
]

#: Dimensione dei blocchi letti e scritti: memoria costante anche sui file grandi.
BLOCCO: Final = 1024 * 1024
#: Il tempo massimo di attesa di una risposta, in secondi.
TEMPO_MASSIMO: Final = 60
#: Le connessioni tenute aperte insieme verso un archivio. Con centinaia di
#: file di pochi megabyte il tempo e' quello di apertura di ogni connessione,
#: non la banda: poche connessioni lo dividono per altrettanto. Di piu' non
#: aiuta, e verso un archivio pubblico e' il modo di farsi limitare (429).
CONNESSIONI: Final = 4
#: I guasti di rete che si ritentano e che, esauriti i tentativi, diventano un
#: :class:`ErroreScarico`: errori di urllib e del sistema (connessione
#: rifiutata o caduta, tempo scaduto, nome non risolto) e risposte HTTP
#: interrotte o malformate, che non sono errori del sistema.
GUASTI_DI_RETE: Final = (urllib.error.URLError, OSError, http.client.HTTPException)


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
    try:
        risposta = urllib.request.urlopen(richiesta, timeout=TEMPO_MASSIMO)
    except urllib.error.HTTPError as guasto:
        # 416: il file parziale e' gia' lungo quanto il file, o di piu' (uno
        # scarico interrotto fra la fine del trasferimento e la verifica): non
        # c'e' un resto da chiedere. Si riparte da capo invece di fallire a ogni
        # rilancio.
        if not (gia and guasto.code == 416):
            raise
        parziale.unlink()
        gia = 0
        richiesta.remove_header("Range")
        risposta = urllib.request.urlopen(richiesta, timeout=TEMPO_MASSIMO)
    with risposta:
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
        except GUASTI_DI_RETE as guasto:
            ultimo = f"trasferimento interrotto ({type(guasto).__name__}: {guasto})"
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
    raise ErroreScarico(
        f"{destinazione.name}: non ottenuto dopo {tentativi} tentativi, {ultimo}. Controlla "
        "la connessione e rilancia lo script: riprende dai file mancanti")


def scarica_in_parallelo(
    voci: Sequence[tuple[Atteso, Path]],
    *,
    connessioni: int = CONNESSIONI,
    tentativi: int = 4,
    pausa: float = 0.5,
) -> Iterator[tuple[int, ErroreScarico | None]]:
    """Scarica piu' file tenendo aperte al piu' ``connessioni`` connessioni.

    Restituisce, man mano che i file si concludono, la posizione di ciascuno
    in ``voci`` e l'errore, o ``None`` se il file e' stato ottenuto integro.
    ``pausa`` separa una richiesta dalla successiva sulla stessa connessione.
    Un guasto su un file non ferma gli altri.
    """
    def uno(posizione: int) -> tuple[int, ErroreScarico | None]:
        atteso, destinazione = voci[posizione]
        try:
            scarica(atteso, destinazione, tentativi=tentativi)
        except ErroreScarico as guasto:
            return posizione, guasto
        finally:
            time.sleep(pausa)
        return posizione, None

    gruppo = ThreadPoolExecutor(max_workers=max(1, connessioni))
    try:
        for concluso in as_completed([gruppo.submit(uno, i) for i in range(len(voci))]):
            yield concluso.result()
    finally:
        # Un'interruzione (Ctrl-C) non deve aspettare i file ancora in coda: si
        # annullano, e si attendono solo i trasferimenti gia' cominciati.
        gruppo.shutdown(wait=True, cancel_futures=True)
