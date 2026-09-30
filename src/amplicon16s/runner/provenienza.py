"""La provenienza di una fase: da quale codice è stato calcolato un risultato.

L'impronta di una fase (``calcolata_su``) dice su quali parametri, fasi a
monte e dati è stata calcolata. Del codice registra solo la **versione** che la
fase dichiara (:attr:`~amplicon16s.steps.base.PipelineStep.versione`), da
incrementare quando cambia ciò che la fase calcola: una modifica senza effetto
sui risultati, un commento o una riorganizzazione, non deve far ricalcolare la
catena.

Accanto alla versione il manifesto registra la **provenienza**: l'impronta del
sorgente della fase, il commit del repository e l'immagine dichiarata. Non
decide la validità: una provenienza diversa a parità di versione produce un
avviso alla ripresa, non un ricalcolo (:mod:`amplicon16s.runner.project`).

**Il sorgente di una fase** è il suo modulo Python, i moduli in cui vive il suo
calcolo (:attr:`~amplicon16s.steps.base.PipelineStep.moduli_sorgente`), il suo
script R e i file di ``R/lib`` che lo script carica. L'impronta si calcola sui
file **effettivamente usati**: i moduli importati e gli script nella cartella
che il ponte esegue (:func:`~amplicon16s.rbridge.runner.cartella_r`), non le
copie del repository. Un'immagine che esegue script R diversi da quelli del
repository (per esempio perché ``AMPLICON16S_R_DIR`` punta agli script copiati
quando l'immagine è stata costruita) lascia quindi nel manifesto un'impronta
diversa da quella registrata nel repository, invece di passare inosservata.

**Il registro** (``steps/registro_sorgente.json``) tiene, per ogni fase
realizzata, la versione e l'impronta del sorgente del repository. Un test
fallisce se il sorgente cambia senza aggiornarlo, e lo strumento
``scripts/registro_sorgente.py`` lo aggiorna obbligando a scegliere fra
incrementare la versione e dichiarare la modifica senza effetto.
"""

from __future__ import annotations

import functools
import hashlib
import importlib
import json
import re
import subprocess
from collections.abc import Mapping
from pathlib import Path
from typing import TYPE_CHECKING, Any, Final

if TYPE_CHECKING:
    from amplicon16s.config.schema import Config
    from amplicon16s.steps.base import PipelineStep

__all__ = [
    "REGISTRO",
    "RADICE_REPOSITORY",
    "cartella_r_repository",
    "RegistroNonAggiornabile",
    "aggiorna_registro",
    "differenze_registro",
    "file_del_sorgente",
    "impronta_sorgente",
    "leggi_registro",
    "provenienza",
    "scrivi_registro",
    "voce_registro",
]

#: La radice del repository quando il pacchetto si usa dai sorgenti; installato
#: (come nell'immagine) e' invece una cartella dell'ambiente Python.
RADICE_REPOSITORY: Final = Path(__file__).resolve().parents[3]
#: Il registro delle versioni e delle impronte del sorgente di ogni fase: e' un
#: dato del pacchetto, installato accanto ai moduli delle fasi.
REGISTRO: Final = Path(__file__).resolve().parents[1] / "steps" / "registro_sorgente.json"


def cartella_r_repository() -> Path:
    """Gli script R del repository, su cui si calcola il registro.

    Dai sorgenti sono quelli accanto a ``src/``; nel pacchetto installato
    nell'immagine sono quelli copiati dal repository quando e' stata costruita,
    che il ponte usa (:func:`~amplicon16s.rbridge.runner.cartella_r`).
    """
    candidata = RADICE_REPOSITORY / "R"
    if candidata.is_dir():
        return candidata
    from amplicon16s.rbridge.runner import cartella_r

    return cartella_r()

#: Il caricamento dei file condivisi negli script R:
#: ``for (f in c("io_json.R", "errors.R")) { source(...) }``.
_CARICAMENTO: Final = re.compile(r"for\s*\(\s*f\s+in\s+c\(([^)]*)\)\s*\)")
_NOME_R: Final = re.compile(r'"([^"]+\.R)"')


def _senza_commenti(testo: str) -> str:
    """Il testo di uno script R senza le righe di commento."""
    return "\n".join(r for r in testo.splitlines() if not r.lstrip().startswith("#"))


def _librerie_caricate(script: Path) -> list[str]:
    """I file di ``R/lib`` che lo script carica, nell'ordine in cui li carica.

    Si fermano i casi che non si sanno leggere: un ``source`` fuori dal ciclo
    di caricamento renderebbe l'impronta incompleta senza avviso.
    """
    testo = _senza_commenti(script.read_text(encoding="utf-8"))
    cicli = _CARICAMENTO.findall(testo)
    if testo.count("source(") != len(cicli):
        raise ValueError(
            f"{script.name}: non riconosco tutti i file caricati con source(); "
            "usa il ciclo for (f in c(...)) sui file di R/lib"
        )
    return [n for ciclo in cicli for n in _NOME_R.findall(ciclo)]


def file_del_sorgente(fase: PipelineStep, cartella_r: Path) -> dict[str, Path]:
    """I file del sorgente di una fase, per nome logico.

    Il nome logico non dipende da dove stanno i file (``python/<modulo>``,
    ``R/<script>``, ``R/lib/<file>``): la stessa copia in due posti ha la
    stessa impronta.
    """
    file: dict[str, Path] = {}
    for modulo in (type(fase).__module__, *fase.moduli_sorgente):
        percorso = importlib.import_module(modulo).__file__
        assert percorso is not None, modulo
        file[f"python/{modulo}"] = Path(percorso)
    if fase.script_r is not None:
        script = cartella_r / fase.script_r
        file[f"R/{fase.script_r}"] = script
        for libreria in _librerie_caricate(script):
            file[f"R/lib/{libreria}"] = cartella_r / "lib" / libreria
    return file


def impronta_sorgente(file: Mapping[str, Path]) -> tuple[str, dict[str, str]]:
    """L'impronta complessiva del sorgente e quella di ciascun file."""
    per_file = {
        nome: "sha256:" + hashlib.sha256(Path(p).read_bytes()).hexdigest()
        for nome, p in sorted(file.items())
    }
    canonico = json.dumps(per_file, sort_keys=True, separators=(",", ":"))
    return "sha256:" + hashlib.sha256(canonico.encode("utf-8")).hexdigest(), per_file


@functools.cache
def _stato_git() -> tuple[str | None, bool | None]:
    """Commit del repository e presenza di modifiche non committate.

    ``(None, None)`` se il codice non sta in un repository git leggibile, come
    nell'immagine, dove ``.git`` non viene copiato: la provenienza lo dice
    invece di inventare un commit.
    """
    try:
        commit = subprocess.run(
            ["git", "-C", str(RADICE_REPOSITORY), "rev-parse", "HEAD"],
            capture_output=True, text=True, timeout=10, check=True,
        ).stdout.strip()
        modifiche = subprocess.run(
            ["git", "-C", str(RADICE_REPOSITORY), "status", "--porcelain", "--", "src", "R"],
            capture_output=True, text=True, timeout=10, check=True,
        ).stdout.strip()
    except (OSError, subprocess.SubprocessError):
        return None, None
    return commit or None, bool(modifiche)


def provenienza(fase: PipelineStep, config: Config) -> dict[str, Any]:
    """La provenienza da registrare nel manifesto di una fase conclusa.

    ``sorgente`` è calcolata sui file effettivamente usati: i moduli importati
    e gli script nella cartella che il ponte esegue.
    """
    from amplicon16s.rbridge.runner import cartella_r

    impronta, per_file = impronta_sorgente(file_del_sorgente(fase, cartella_r()))
    commit, modifiche = _stato_git()
    return {
        "versione": fase.versione,
        "sorgente": impronta,
        "file": per_file,
        "commit": commit,
        "modifiche_non_committate": modifiche,
        "immagine": config.run.container,
    }


def voce_registro(fase: PipelineStep, cartella_r: Path | None = None) -> dict[str, Any]:
    """La voce del registro per una fase: versione e impronta delle copie del repository."""
    cartella = cartella_r if cartella_r is not None else cartella_r_repository()
    impronta, per_file = impronta_sorgente(file_del_sorgente(fase, cartella))
    return {"versione": fase.versione, "sorgente": impronta, "file": sorted(per_file)}


def leggi_registro(percorso: Path = REGISTRO) -> dict[str, dict[str, Any]]:
    """Il registro delle versioni e delle impronte, per fase."""
    return json.loads(percorso.read_text(encoding="utf-8"))


# --------------------------------------------------------------------------- #
# Il registro e il suo aggiornamento                                           #
# --------------------------------------------------------------------------- #


class RegistroNonAggiornabile(ValueError):
    """L'aggiornamento chiesto lascerebbe una modifica senza una scelta esplicita."""


def differenze_registro(
    registro: Mapping[str, Mapping[str, Any]],
    fasi: Mapping[str, PipelineStep],
    cartella_r: Path | None = None,
) -> dict[str, str]:
    """Per ogni fase che non corrisponde al registro, in che cosa differisce.

    I casi sono quattro: ``nuova`` (la fase non e' registrata), ``rimossa``
    (registrata ma non realizzata), ``versione incrementata`` (il sorgente e'
    cambiato e la versione e' salita: la ripresa rifara' la fase) e
    ``sorgente cambiato a parita' di versione``, il caso che chiede una scelta.
    Una versione scesa o cambiata senza modifiche al sorgente e' un errore.
    """
    differenze: dict[str, str] = {}
    for nome in sorted(set(registro) | set(fasi)):
        if nome not in fasi:
            differenze[nome] = "rimossa"
            continue
        attuale = voce_registro(fasi[nome], cartella_r)
        if nome not in registro:
            differenze[nome] = "nuova"
            continue
        registrata = registro[nome]
        if attuale["versione"] < registrata["versione"]:
            differenze[nome] = "versione scesa: errore"
        elif attuale["sorgente"] == registrata["sorgente"]:
            if attuale["versione"] != registrata["versione"] or attuale["file"] != registrata["file"]:
                differenze[nome] = "versione o file cambiati senza modifiche al sorgente: errore"
        elif attuale["versione"] > registrata["versione"]:
            differenze[nome] = "versione incrementata"
        else:
            differenze[nome] = "sorgente cambiato a parita' di versione"
    return differenze


def aggiorna_registro(
    registro: Mapping[str, Mapping[str, Any]],
    fasi: Mapping[str, PipelineStep],
    senza_effetto: frozenset[str] = frozenset(),
    cartella_r: Path | None = None,
) -> dict[str, dict[str, Any]]:
    """Il registro aggiornato, se ogni modifica ha una scelta esplicita.

    Una fase il cui sorgente e' cambiato a parita' di versione si aggiorna solo
    se compare in ``senza_effetto``: chi l'ha modificata dichiara che la
    modifica non cambia cio' che la fase calcola. Altrimenti deve incrementarne
    la versione. ``senza_effetto`` non puo' nominare fasi che non ne hanno
    bisogno: una dichiarazione superflua nasconderebbe un errore di battitura.
    """
    differenze = differenze_registro(registro, fasi, cartella_r)
    errori = [f"{n}: {d}" for n, d in differenze.items() if d.endswith("errore")]
    da_scegliere = {n for n, d in differenze.items() if d.startswith("sorgente cambiato")}
    senza_scelta = sorted(da_scegliere - senza_effetto)
    if senza_scelta:
        errori.append(
            f"{', '.join(senza_scelta)}: sorgente cambiato a parita' di versione. Se "
            "cambia cio' che la fase calcola, incrementa la sua versione; se la "
            "modifica e' senza effetto sui risultati, dichiaralo con --senza-effetto"
        )
    superflue = sorted(senza_effetto - da_scegliere)
    if superflue:
        errori.append(f"--senza-effetto nomina fasi che non ne hanno bisogno: {', '.join(superflue)}")
    if errori:
        raise RegistroNonAggiornabile("; ".join(errori))
    return {nome: voce_registro(fasi[nome], cartella_r) for nome in sorted(fasi)}


def scrivi_registro(registro: Mapping[str, Any], percorso: Path = REGISTRO) -> None:
    """Scrive il registro in forma stabile: fasi e chiavi ordinate."""
    percorso.write_text(
        json.dumps(registro, indent=2, ensure_ascii=False, sort_keys=True) + "\n",
        encoding="utf-8",
    )
