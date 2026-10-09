"""Orchestrazione di un'analisi: copia dell'oggetto, validazione, calcolo,
manifesto.

L'oggetto consegnato dalla pipeline non si tocca: se ne calcola lo SHA-256, lo
si copia in una cartella temporanea, R lavora sulla copia, e alla fine lo
SHA-256 dell'originale si ricalcola. Le uscite si scrivono in una cartella
temporanea e passano in quella indicata solo ad analisi conclusa: un'analisi
interrotta non lascia una cartella di uscita a meta'.

Le uscite non contengono date, durate, percorsi ne' identificativi di
esecuzione: a parita' di oggetto, configurazione e seme i file sono identici
byte per byte. I file del contratto con R, che portano l'identificativo
dell'invocazione, restano nella cartella temporanea.
"""

from __future__ import annotations

import json
import shutil
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Final

from amplicon16s.io_layer.artifacts import scrivi_atomico
from amplicon16s.io_layer.checksums import checksum_file
from amplicon16s.metadata.crosswalk import _VALORI_ASSENTI

from amplicon16s_eco import __version__
from amplicon16s_eco.catalogo import ErroreEco
from amplicon16s_eco.config import ConfigEco, digest, risolta
from amplicon16s_eco.ponte import esegui_r

__all__ = ["Avviso", "Esito", "NOME_MANIFESTO", "esegui", "valida_sull_oggetto"]

NOME_MANIFESTO: Final = "manifest.json"
NOME_CONFIGURAZIONE: Final = "configurazione.json"
NOME_RIEPILOGO: Final = "riepilogo.json"


@dataclass(frozen=True)
class Avviso:
    """Un avviso del catalogo con il suo dettaglio."""

    codice: str
    dettaglio: str


@dataclass(frozen=True)
class Esito:
    """Cio' che un'analisi, o la sua sola validazione, ha stabilito."""

    #: SHA-256 dell'oggetto di partenza, nella forma ``sha256:<esadecimale>``.
    oggetto_sha256: str
    #: Digest della configurazione risolta.
    digest_configurazione: str
    avvisi: tuple[Avviso, ...]
    #: Il riepilogo del disegno (validazione) o dell'analisi (esecuzione).
    riepilogo: dict[str, Any]
    #: I file prodotti, relativi alla cartella di uscita; vuoto per la validazione.
    file: tuple[str, ...] = ()


def _json(documento: Any) -> str:
    """JSON in forma stabile: chiavi nell'ordine dato, due spazi, a capo finale."""
    return json.dumps(documento, ensure_ascii=False, indent=2) + "\n"


def _parametri_r(config: ConfigEco, copia: Path) -> dict[str, Any]:
    """I parametri per gli script R: la configurazione risolta, la copia
    dell'oggetto e le convenzioni di valore mancante della pipeline.
    """
    return {
        **risolta(config),
        "oggetto": str(copia),
        "valori_mancanti": sorted(_VALORI_ASSENTI),
    }


def _verifica_oggetto(oggetto: Path) -> None:
    if not oggetto.is_file():
        raise ErroreEco([("E-ECO-04", f"File non trovato: {oggetto}.")])


def _verifica_uscita(oggetto: Path, uscita: Path) -> None:
    """La cartella di uscita deve essere assente o vuota e fuori dalla cartella
    dell'oggetto.
    """
    casa_oggetto = oggetto.resolve().parent
    assoluta = uscita.resolve()
    if assoluta == casa_oggetto or casa_oggetto in assoluta.parents:
        raise ErroreEco([(
            "E-ECO-09",
            f"{uscita} sta nella cartella dell'oggetto di partenza ({casa_oggetto}).",
        )])
    if assoluta.exists() and not assoluta.is_dir():
        raise ErroreEco([("E-ECO-09", f"{uscita} esiste e non e' una cartella.")])
    if assoluta.is_dir() and any(assoluta.iterdir()):
        raise ErroreEco([("E-ECO-09", f"{uscita} non e' vuota.")])


def _valida(config: ConfigEco, copia: Path, lavoro: Path) -> tuple[tuple[Avviso, ...], dict[str, Any]]:
    """Esegue la validazione sull'oggetto e ne traduce l'esito."""
    esegui_r("valida.R", _parametri_r(config, copia), lavoro, "validazione")
    esito = json.loads((lavoro / "validazione.json").read_text(encoding="utf-8"))
    errori = [(e["codice"], e["dettaglio"]) for e in esito.get("errori") or []]
    if errori:
        raise ErroreEco(errori)
    avvisi = tuple(Avviso(a["codice"], a["dettaglio"]) for a in esito.get("avvisi") or [])
    return avvisi, esito.get("disegno") or {}


def _copia(oggetto: Path, lavoro: Path, atteso: str) -> Path:
    """Copia l'oggetto nella cartella di lavoro e verifica che sia identica."""
    copia = lavoro / "oggetto.rds"
    shutil.copyfile(oggetto, copia)
    if checksum_file(copia) != atteso:
        raise ErroreEco([("E-ECO-91", f"La copia di {oggetto} non coincide con l'originale.")])
    return copia


def valida_sull_oggetto(oggetto: Path | str, config: ConfigEco) -> Esito:
    """Valida la configurazione contro l'oggetto, senza calcolare nulla.

    Solleva :class:`~amplicon16s_eco.catalogo.ErroreEco` con tutti i rifiuti;
    altrimenti restituisce gli avvisi e il disegno dell'analisi.
    """
    oggetto = Path(oggetto)
    _verifica_oggetto(oggetto)
    impronta = checksum_file(oggetto)
    with tempfile.TemporaryDirectory(prefix="amplicon16s-eco-") as temporanea:
        lavoro = Path(temporanea)
        avvisi, disegno = _valida(config, _copia(oggetto, lavoro, impronta), lavoro)
    return Esito(impronta, digest(config), avvisi, disegno)


def esegui(oggetto: Path | str, config: ConfigEco, uscita: Path | str) -> Esito:
    """Esegue le analisi sull'oggetto e scrive le uscite nella cartella indicata.

    L'ordine e': verifiche sui percorsi, impronta e copia dell'oggetto,
    validazione sull'oggetto, calcolo, manifesto, verifica finale dell'impronta
    dell'originale. Nulla viene scritto nella cartella di uscita finche' il
    calcolo non e' concluso.
    """
    oggetto, uscita = Path(oggetto), Path(uscita)
    _verifica_oggetto(oggetto)
    _verifica_uscita(oggetto, uscita)
    impronta = checksum_file(oggetto)

    with tempfile.TemporaryDirectory(prefix="amplicon16s-eco-") as temporanea:
        lavoro = Path(temporanea)
        copia = _copia(oggetto, lavoro, impronta)
        _valida(config, copia, lavoro)

        prodotta = lavoro / "uscita"
        prodotta.mkdir()
        esegui_r("analisi.R", {**_parametri_r(config, copia), "uscita": str(prodotta)},
                 lavoro, "analisi")
        riepilogo = json.loads((prodotta / NOME_RIEPILOGO).read_text(encoding="utf-8"))

        if checksum_file(oggetto) != impronta:
            raise ErroreEco([("E-ECO-91", f"L'impronta di {oggetto} non e' piu' {impronta}.")])

        scrivi_atomico(prodotta / NOME_CONFIGURAZIONE, _json({
            "digest": digest(config),
            "parametri": risolta(config),
        }))
        elencati = sorted(
            p.relative_to(prodotta).as_posix() for p in prodotta.rglob("*") if p.is_file()
        )
        scrivi_atomico(prodotta / NOME_MANIFESTO, _json({
            "pacchetto": {"nome": "amplicon16s_eco", "versione": __version__},
            "oggetto_di_partenza": {"sha256": impronta, "byte": oggetto.stat().st_size},
            "configurazione": {"digest": digest(config), "file": NOME_CONFIGURAZIONE},
            "file": [
                {"nome": nome, "byte": (prodotta / nome).stat().st_size,
                 "sha256": checksum_file(prodotta / nome)}
                for nome in elencati
            ],
        }))

        uscita.mkdir(parents=True, exist_ok=True)
        for voce in sorted(prodotta.iterdir()):
            shutil.move(str(voce), str(uscita / voce.name))

    avvisi = tuple(Avviso(a["codice"], a["dettaglio"]) for a in riepilogo.get("avvisi") or [])
    return Esito(impronta, digest(config), avvisi, riepilogo,
                 tuple([*elencati, NOME_MANIFESTO]))
