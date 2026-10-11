"""Scarico delle letture forward di un dataset depositato a letture appaiate.

La pipeline tratta dati single-end: di ogni corsa appaiata usa il solo file le
cui letture partono dal primer forward. Il nome del file nel deposito (R1,
``_1``) lo dice quasi sempre, ma non sempre: in alcuni depositi i due file di
una corsa sono scambiati. Questo modulo, comune agli script dei dataset,

1. legge le prime letture dei due file di ogni corsa (i soli primi byte, con
   una richiesta parziale) e conta quante iniziano con il primer forward;
2. tiene il file che lo porta. Se e' il secondo, la corsa e' **invertita** e lo
   dichiara. Se nessuno dei due lo porta (letture depositate senza primer)
   l'inversione non si puo' stabilire: tiene il file dichiarato forward e lo
   dichiara;
3. scarica per intero il solo file tenuto, lo verifica (dimensione dichiarata
   dalla fonte e MD5) e lo salva con un nome senza marcatori di coppia
   (``R1``, ``R2``, ``_1``, ``_2``): con una corsa invertita il file forward e'
   il secondo del deposito, e con quel nome la pipeline lo respingerebbe.

L'elenco delle corse e' una tabella TSV accanto allo script del dataset, con le
colonne ``file`` (il nome da dare al file tenuto), ``url_1``, ``byte_1``,
``md5_1``, ``url_2``, ``byte_2``, ``md5_2`` (come li dichiara la fonte; un MD5
che la fonte non dichiara e' vuoto e si registra al primo scarico con
``--registra``) e ``lato`` (``1`` o ``2``: il file tenuto, registrato allo
stesso modo). Lo scarico e' ripetibile: un file gia' presente e conforme non si
scarica di nuovo.

**Entrambe le letture, per la modalita' appaiata** (``--coppie``). I due file
di ogni corsa si scaricano in una cartella separata, ``fastq_coppie/``, con il
nome del file forward piu' ``_R1`` per il primo file del deposito e ``_R2`` per
il secondo: **nell'ordine del deposito, senza scambi**. Il riconoscimento del
lato forward qui non si applica di proposito: se in una corsa il deposito ha i
due file invertiti deve accorgersene la pipeline, che in modalita' appaiata lo
verifica sul contenuto. Ogni file si verifica con la dimensione e, dove c'e',
con l'MD5 dell'elenco; con ``--registra`` gli MD5 che la fonte non dichiara si
registrano dal primo scarico. La cartella ``fastq/`` non viene toccata.
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
import time
import urllib.request
import zlib
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from scaricamento import (  # noqa: E402
    GUASTI_DI_RETE,
    Atteso,
    ErroreScarico,
    md5_file,
    scarica,
    valido,
)

#: Quanti byte di ogni file si leggono per riconoscere il primer, e quante
#: letture al massimo.
BYTE_DI_TESTA = 400_000
LETTURE_DI_TESTA = 2000
#: Frazione di letture con il primer in testa oltre la quale il file lo porta.
FRAZIONE_MINIMA = 0.5

_IUPAC = {"A": "A", "C": "C", "G": "G", "T": "T", "R": "[AG]", "Y": "[CT]", "S": "[CG]",
          "W": "[AT]", "K": "[GT]", "M": "[AC]", "B": "[CGT]", "D": "[AGT]", "H": "[ACT]",
          "V": "[ACG]", "N": "[ACGT]"}


class CoppiaAmbigua(ValueError):
    """Il primer forward apre le letture di tutti e due i file di una corsa."""


def espressione(primer: str) -> re.Pattern[str]:
    """L'espressione regolare di un primer in codici IUPAC."""
    return re.compile("".join(_IUPAC[base] for base in primer.upper()))


def letture_di_testa(dati: bytes) -> list[str]:
    """Le sequenze delle prime letture di un FASTQ compresso, dai suoi primi
    byte: l'ultima riga, che il taglio puo' avere interrotto, si scarta."""
    decompressore = zlib.decompressobj(16 + zlib.MAX_WBITS)
    try:
        testo = decompressore.decompress(dati)
    except zlib.error:
        testo = b""
    righe = testo.decode("ascii", "replace").split("\n")[:-1]
    return [righe[i] for i in range(1, len(righe), 4)][:LETTURE_DI_TESTA]


def frazione_con_primer(letture: list[str], primer: str) -> float:
    """La frazione di letture che iniziano con il primer."""
    if not letture:
        return 0.0
    cerca = espressione(primer)
    return sum(bool(cerca.match(lettura)) for lettura in letture) / len(letture)


def lato_forward(frazione_1: float, frazione_2: float) -> tuple[str, str]:
    """Quale file di una corsa tenere, e che cosa dichiarare.

    Restituisce il lato (``"1"`` o ``"2"``) e l'esito: ``forward`` (il primer
    apre il primo file), ``invertita`` (apre il secondo) o ``senza primer``
    (nessuno dei due: si tiene il primo, senza poterlo verificare).
    """
    porta_1, porta_2 = frazione_1 >= FRAZIONE_MINIMA, frazione_2 >= FRAZIONE_MINIMA
    if porta_1 and porta_2:
        raise CoppiaAmbigua("il primer forward apre le letture di entrambi i file")
    if porta_1:
        return "1", "forward"
    if porta_2:
        return "2", "invertita"
    return "1", "senza primer"


def _testa(url: str, tentativi: int = 3) -> bytes:
    """I primi byte di un file remoto; una connessione caduta si ritenta."""
    richiesta = urllib.request.Request(
        url, headers={"User-Agent": "amplicon16s-dati", "Range": f"bytes=0-{BYTE_DI_TESTA - 1}"})
    for tentativo in range(1, tentativi + 1):
        try:
            with urllib.request.urlopen(richiesta, timeout=120) as risposta:
                return risposta.read(BYTE_DI_TESTA)
        except GUASTI_DI_RETE as guasto:
            if tentativo == tentativi:
                raise ErroreScarico(
                    f"{url}: prime letture non ottenute dopo {tentativi} tentativi "
                    f"({type(guasto).__name__}: {guasto}). Controlla la connessione e "
                    "rilancia lo script: riprende dai file mancanti") from guasto
            time.sleep(5.0 * tentativo)
    raise AssertionError("non raggiungibile")


def leggi_elenco(percorso: Path) -> list[dict[str, str]]:
    """Le corse dell'elenco."""
    with open(percorso, encoding="utf-8", newline="") as file:
        return list(csv.DictReader(file, delimiter="\t"))


def scrivi_elenco(percorso: Path, righe: list[dict[str, str]]) -> None:
    """Riscrive l'elenco, con gli MD5 e i lati registrati."""
    with open(percorso, "w", encoding="utf-8", newline="") as file:
        scrittore = csv.DictWriter(file, fieldnames=list(righe[0]), delimiter="\t",
                                   lineterminator="\n")
        scrittore.writeheader()
        scrittore.writerows(righe)


def prepara_corsa(riga: dict[str, str], cartella: Path, primer: str, *, registra: bool,
                  tentativi: int) -> str:
    """Porta in ``cartella/fastq`` il file forward di una corsa e restituisce
    l'esito (``forward``, ``invertita`` o ``senza primer``). Con ``registra``
    scrive nella riga il lato tenuto e, se la fonte non lo dichiara, il suo MD5.
    """
    frazioni = [frazione_con_primer(letture_di_testa(_testa(riga[f"url_{lato}"], tentativi)), primer)
                for lato in ("1", "2")]
    lato, esito = lato_forward(*frazioni)
    if riga["lato"] and riga["lato"] != lato:
        raise ErroreScarico(
            f"{riga['file']}: il primer indica il file {lato}, l'elenco registra il {riga['lato']}")
    _scarica_un_lato(riga, lato, cartella / "fastq" / riga["file"], registra=registra,
                     tentativi=tentativi)
    if registra:
        riga["lato"] = lato
    return esito


def _scarica_un_lato(riga: dict[str, str], lato: str, destinazione: Path, *, registra: bool,
                     tentativi: int) -> None:
    """Scarica in ``destinazione`` il file ``lato`` (``"1"`` o ``"2"``) di una
    corsa e lo verifica: con l'MD5 dell'elenco se c'e'; altrimenti, con
    ``registra``, con la sola dimensione, e l'MD5 del file ottenuto entra nella
    riga come riferimento per gli scarichi successivi.
    """
    md5 = riga[f"md5_{lato}"]
    if md5:
        scarica(Atteso(riga[f"url_{lato}"], md5, int(riga[f"byte_{lato}"])), destinazione,
                tentativi=tentativi)
        return
    if not registra:
        raise ErroreScarico(
            f"{riga['file']}: l'elenco non ha l'MD5 del file {lato}; registralo con --registra")
    # La fonte dichiara la sola dimensione: si verifica quella, e l'MD5 del
    # file ottenuto diventa il riferimento per gli scarichi successivi.
    parziale = destinazione.with_name(destinazione.name + ".parziale")
    destinazione.parent.mkdir(parents=True, exist_ok=True)
    ultimo = "nessun tentativo"
    for tentativo in range(1, tentativi + 1):
        try:
            urllib.request.urlretrieve(riga[f"url_{lato}"], parziale)
        except GUASTI_DI_RETE as guasto:
            ultimo = f"{type(guasto).__name__}: {guasto}"
            time.sleep(5.0 * tentativo)
            continue
        if parziale.stat().st_size == int(riga[f"byte_{lato}"]):
            parziale.replace(destinazione)
            break
        ultimo = f"{parziale.stat().st_size} byte su {riga[f'byte_{lato}']}"
    else:
        raise ErroreScarico(f"{riga['file']}: {ultimo}")
    riga[f"md5_{lato}"] = md5_file(destinazione)


#: La cartella delle due letture di ogni corsa (``--coppie``), separata da
#: quella delle sole forward, e i marcatori che ne distinguono i due file.
CARTELLA_COPPIE = "fastq_coppie"
MARCATORI_DI_COPPIA = {"1": "_R1", "2": "_R2"}


def nome_di_coppia(nome: str, lato: str) -> str:
    """Il nome, in :data:`CARTELLA_COPPIE`, del file ``lato`` del deposito per
    la corsa il cui file forward si chiama ``nome``: il marcatore di coppia
    prima dell'estensione. Il primo file del deposito e' sempre ``_R1``.
    """
    for estensione in (".fastq.gz", ".fq.gz", ".fastq", ".fq"):
        if nome.endswith(estensione):
            return nome[:-len(estensione)] + MARCATORI_DI_COPPIA[lato] + estensione
    raise ValueError(f"{nome}: estensione non riconosciuta")


def coppie(elenco: Path, opzioni: argparse.Namespace) -> int:
    """Scarica e verifica i due file di ogni corsa, nell'ordine del deposito."""
    righe = leggi_elenco(elenco)
    conformi, guasti = 0, []
    for indice, riga in enumerate(righe, start=1):
        percorsi = {lato: opzioni.cartella / CARTELLA_COPPIE / nome_di_coppia(riga["file"], lato)
                    for lato in ("1", "2")}
        presenti = all(
            riga[f"md5_{lato}"] and valido(percorso, Atteso(
                riga[f"url_{lato}"], riga[f"md5_{lato}"], int(riga[f"byte_{lato}"])))
            for lato, percorso in percorsi.items())
        if presenti:
            conformi += 1
            continue
        if opzioni.solo_verifica:
            guasti.append(f"{riga['file']}: i due file assenti o non conformi")
            continue
        try:
            for lato, percorso in percorsi.items():
                _scarica_un_lato(riga, lato, percorso, registra=opzioni.registra,
                                 tentativi=opzioni.tentativi)
        except (ErroreScarico, OSError) as guasto:
            guasti.append(f"{riga['file']}: {type(guasto).__name__}: {guasto}")
            continue
        conformi += 1
        print(f"[{indice}/{len(righe)}] {percorsi['1'].name} e {percorsi['2'].name}", flush=True)
        time.sleep(opzioni.pausa)
    if opzioni.registra and not opzioni.solo_verifica:
        scrivi_elenco(elenco, righe)
    print(f"corse {len(righe)}: coppie di file conformi {conformi}, mancanti o non conformi "
          f"{len(guasti)}; i file sono nell'ordine del deposito, senza scambi")
    for voce in guasti:
        print(f"  {voce}")
    return 0 if not guasti else 1


def main(elenco: Path, primer: str, predefinita: Path, descrizione: str,
         argomenti: list[str] | None = None) -> int:
    """La riga di comando comune agli script dei dataset."""
    parser = argparse.ArgumentParser(description=descrizione)
    parser.add_argument("--cartella", type=Path, default=predefinita,
                        help="dove sta, o andra', la cartella fastq/")
    parser.add_argument("--solo-verifica", action="store_true",
                        help="controlla i file presenti, senza scaricare nulla")
    parser.add_argument("--registra", action="store_true",
                        help="registra nell'elenco il lato tenuto e gli MD5 che la fonte non dichiara")
    parser.add_argument("--pausa", type=float, default=0.5)
    parser.add_argument("--tentativi", type=int, default=3)
    parser.add_argument("--coppie", action="store_true",
                        help="scarica i due file di ogni corsa in fastq_coppie/, nell'ordine "
                             "del deposito e senza scambi, per la modalita' appaiata")
    opzioni = parser.parse_args(argomenti)
    if opzioni.coppie:
        return coppie(elenco, opzioni)

    righe = leggi_elenco(elenco)
    esiti: dict[str, list[str]] = {"forward": [], "invertita": [], "senza primer": []}
    conformi, guasti = 0, []
    for indice, riga in enumerate(righe, start=1):
        destinazione = opzioni.cartella / "fastq" / riga["file"]
        lato = riga["lato"]
        if lato and riga[f"md5_{lato}"] and valido(destinazione, Atteso(
                riga[f"url_{lato}"], riga[f"md5_{lato}"], int(riga[f"byte_{lato}"]))):
            conformi += 1
            continue
        if opzioni.solo_verifica:
            guasti.append(f"{riga['file']}: {'non conforme' if destinazione.exists() else 'assente'}")
            continue
        try:
            esito = prepara_corsa(riga, opzioni.cartella, primer, registra=opzioni.registra,
                                  tentativi=opzioni.tentativi)
        except (ErroreScarico, CoppiaAmbigua, OSError) as guasto:
            guasti.append(f"{riga['file']}: {type(guasto).__name__}: {guasto}")
            continue
        esiti[esito].append(riga["file"])
        conformi += 1
        print(f"[{indice}/{len(righe)}] {riga['file']}: {esito}", flush=True)
        time.sleep(opzioni.pausa)
    if opzioni.registra and not opzioni.solo_verifica:
        scrivi_elenco(elenco, righe)
    print(f"corse {len(righe)}: file forward conformi {conformi}, mancanti o non conformi {len(guasti)}")
    if esiti["invertita"]:
        print(f"corse invertite (tenuto il secondo file): {', '.join(esiti['invertita'])}")
    if esiti["senza primer"]:
        print(f"letture senza il primer forward in testa in {len(esiti['senza primer'])} corse: "
              "tenuto il file dichiarato forward, l'inversione non e' verificabile")
    for voce in guasti:
        print(f"  {voce}")
    return 0 if not guasti else 1
