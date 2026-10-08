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
calcolo (:func:`moduli_del_calcolo`), il suo script R e i file di ``R/lib`` che
lo script carica. I moduli del calcolo si ricavano dagli import, non da un
elenco scritto a mano che potrebbe restare indietro: sono i moduli del pacchetto
che il modulo della fase importa, direttamente o attraverso altri, esclusi
l'infrastruttura elencata con la sua ragione in :data:`ESCLUSI` e i moduli
delle altre fasi. Da questi ultimi una fase importa solo nomi di artefatti, e
un test lo verifica: una funzione condivisa fra fasi sta in un modulo comune,
e cosi' entra nel sorgente di ciascuna. Il criterio e' chi puo' cambiare cio'
che una fase calcola: per questo entrano anche i moduli che producono un dato
che il contesto passa alle fasi senza import (:data:`DAL_CONTESTO`, il lettore
dell'inventario per le fasi che dipendono da S0), e per una fase sola un modulo
altrimenti escluso (:data:`ECCEZIONI_PER_FASE`). L'impronta si calcola sui
file **effettivamente usati**: i moduli importati e gli script nella cartella
che il ponte esegue (:func:`~amplicon16s.rbridge.runner.cartella_r`), non le
copie del repository. Un'immagine che esegue script R diversi da quelli del
repository (per esempio perché ``AMPLICON16S_R_DIR`` punta agli script copiati
quando l'immagine è stata costruita) lascia quindi nel manifesto un'impronta
diversa da quella registrata nel repository, invece di passare inosservata.

**La regola rigorosa** (``run.strict_provenance``, disattivata per difetto e
attiva nelle configurazioni congelate) rende la provenienza vincolante. Va
letta distinguendo ciò che la pipeline **verifica** da ciò che può solo
**registrare come dichiarato**:

* *verificato, il codice*: il repository git è leggibile e ``src/`` e ``R/``
  non hanno modifiche non committate, e gli script R eseguiti sono quelli del
  repository (:func:`verifica_git`); l'impronta del
  sorgente di ogni fase è calcolata sui file effettivamente eseguiti ed entra
  nell'impronta della fase (:func:`impronta_rigorosa`). Il codice eseguito è
  quindi quello del commit registrato;
* *verificato, l'ambiente R*: la versione di R e quelle di tutti i pacchetti
  del file di blocco, compresa la correzione di dada2, sono lette dalle
  librerie installate (``R/00_ambiente.R``) e confrontate con ``run.lockfile``
  (:func:`discordanze_ambiente`); l'impronta del file di blocco entra
  nell'impronta di ogni fase;
* *dichiarato, l'immagine*: ``run.container`` entra nell'impronta di ogni fase,
  ma dall'interno del container il digest dell'immagine in esecuzione non è
  conoscibile. È una dichiarazione di chi lancia l'esecuzione: ciò che la
  sostiene è la verifica dell'ambiente R, non il digest. Non sono verificate le
  librerie di sistema dell'immagine, che il file di blocco non descrive.

Se una verifica fallisce l'esecuzione non parte (``E-PROV-01``, ``E-PROV-02``,
``E-PROV-03``): in particolare, se git non è leggibile la regola rifiuta,
invece di lasciar passare un codice che nulla identifica.

**Il registro** (``steps/registro_sorgente.json``) tiene, per ogni fase
realizzata, la versione e l'impronta del sorgente del repository. Un test
fallisce se il sorgente cambia senza aggiornarlo, e lo strumento
``scripts/registro_sorgente.py`` lo aggiorna obbligando a scegliere fra
incrementare la versione e dichiarare la modifica senza effetto.
"""

from __future__ import annotations

import ast
import functools
import hashlib
import importlib
import json
import os
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
    "DAL_CONTESTO",
    "blocco_r",
    "discordanze_ambiente",
    "impronta_rigorosa",
    "verifica_git",
    "ECCEZIONI_PER_FASE",
    "ESCLUSI",
    "file_del_sorgente",
    "impronta_sorgente",
    "moduli_del_calcolo",
    "moduli_della_fase",
    "leggi_registro",
    "provenienza",
    "scrivi_registro",
    "voce_registro",
]

#: La radice del repository quando il pacchetto si usa dai sorgenti; installato
#: (come nell'immagine) e' invece una cartella dell'ambiente Python.
RADICE_REPOSITORY: Final = Path(__file__).resolve().parents[3]
#: Lo script che esegue la pipeline dal codice di un clone, relativo alla sua
#: radice: lo nominano i messaggi della regola rigorosa.
LANCIATORE: Final = "scripts/esegui.py"
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


#: I moduli del pacchetto esclusi dal sorgente delle fasi, ciascuno con la
#: ragione. Il criterio: un modulo resta fuori solo se non puo' cambiare cio'
#: che una fase calcola, o se cio' che determina del risultato entra
#: nell'impronta per un'altra via (i valori dei parametri, le fasi a monte).
#: Restano invece fra i moduli del calcolo di chi li importa
#: ``io_layer.artifacts`` (scrive gli artefatti di Python e ne codifica il testo;
#: calcola i checksum che S6 riporta in un artefatto), ``io_layer.checksums``
#: (la forma di quei checksum), ``rbridge.payload`` (serializza i parametri su
#: cui calcola R), ``rbridge.runner`` (l'ambiente di R: lingua, ``LC_ALL``,
#: librerie, che cambiano per esempio l'intestazione degli .rds e gli
#: ordinamenti) e ``runner.retry`` (calcola i valori corretti di un nuovo
#: tentativo). L'elenco e' esplicito: un modulo nuovo dell'infrastruttura entra
#: nel sorgente finche' non vi si aggiunge con la sua ragione.
ESCLUSI: Final[Mapping[str, str]] = {
    "amplicon16s.config": "pacchetto, vuoto",
    "amplicon16s.config.schema": (
        "valida e normalizza i parametri: i valori, cosi' come risultano, entrano "
        "nell'impronta dei parametri della fase"
    ),
    "amplicon16s.config.defaults": (
        "i valori predefiniti entrano nell'impronta come valori dei parametri"
    ),
    "amplicon16s.config.resolve": (
        "calcola i parametri derivati, che il contesto passa alle fasi: entrano "
        "nell'impronta come valori dei gruppi dichiarati; resolved.yaml non e' un "
        "artefatto di fase"
    ),
    "amplicon16s.config.vista": (
        "restituisce i valori della configurazione senza trasformarli; un accesso "
        "non dichiarato solleva un errore, non cambia un valore"
    ),
    "amplicon16s.errors": "pacchetto, vuoto",
    "amplicon16s.errors.catalog": (
        "la categoria di un codice decide se una fase si ferma, non che cosa "
        "calcola (nessuna fase dichiara un ripiego); i messaggi vanno nel log e nei "
        "manifesti. Fa eccezione S0, in ECCEZIONI_PER_FASE"
    ),
    "amplicon16s.errors.exceptions": (
        "costruisce le eccezioni dai codici del catalogo, senza cambiare un valore"
    ),
    "amplicon16s.logging": "pacchetto, vuoto",
    "amplicon16s.logging.logger": "il log in 99_logs, che non e' un artefatto di fase",
    "amplicon16s.rbridge": "pacchetto, vuoto",
    "amplicon16s.runner": "pacchetto, vuoto",
    "amplicon16s.runner.executor": (
        "ordine di esecuzione, tentativi e punto di ripresa in 99_logs; i valori "
        "corretti di un tentativo li calcola runner.retry, che resta nel sorgente"
    ),
    "amplicon16s.runner.graph": "le dipendenze: entrano nell'impronta come fasi a monte",
    "amplicon16s.runner.project": (
        "costruisce il contesto delle fasi: la configurazione risolta, i cui valori "
        "entrano nell'impronta, e l'inventario, letto da un modulo in DAL_CONTESTO"
    ),
    "amplicon16s.runner.provenienza": "il registro e la provenienza, nei manifesti",
    "amplicon16s.runner.tracciamento": (
        "ricompone il tracciamento per il report; nessuna fase ne legge il risultato"
    ),
    "amplicon16s.steps.base": (
        "lo scheletro delle fasi: prerequisiti, vista dei parametri, manifesto di "
        "fase; il calcolo e' della fase"
    ),
}

#: I moduli che producono dati che il contesto passa alle fasi senza che le
#: fasi li importino, con la fase che produce il dato e la ragione: entrano nel
#: sorgente di ogni fase che dipende da quella nel grafo.
DAL_CONTESTO: Final[Mapping[str, tuple[str, str]]] = {
    "amplicon16s.metadata.lettura_inventario": (
        "S0",
        "rilegge l'inventario di S0 e decide su quali campioni, con quali classi, "
        "piastre e corse, le fasi calcolano",
    ),
}

#: Moduli esclusi che per una fase determinano comunque il contenuto di un suo
#: artefatto, con la ragione.
ECCEZIONI_PER_FASE: Final[Mapping[str, Mapping[str, str]]] = {
    "S0": {
        "amplicon16s.errors.catalog": (
            "il dettaglio di una violazione di G15 riporta la sintesi del catalogo, e "
            "S0 lo scrive in gates.json"
        ),
    },
}

#: I moduli delle fasi: ``amplicon16s.steps.s00_validate`` e cosi' via.
MODULO_DI_FASE: Final = re.compile(r"^amplicon16s\.steps\.s\d{2}_")


def _importati(modulo: str) -> set[str]:
    """I moduli del pacchetto che ``modulo`` importa, dal suo sorgente."""
    percorso = importlib.import_module(modulo).__file__
    assert percorso is not None, modulo
    importati = set()
    for nodo in ast.walk(ast.parse(Path(percorso).read_text(encoding="utf-8"))):
        if isinstance(nodo, ast.ImportFrom) and nodo.module and nodo.level == 0:
            importati.add(nodo.module)
        elif isinstance(nodo, ast.Import):
            importati.update(alias.name for alias in nodo.names)
    return {m for m in importati if m == "amplicon16s" or m.startswith("amplicon16s.")}


def moduli_del_calcolo(modulo: str) -> tuple[str, ...]:
    """I moduli in cui vive il calcolo della fase di ``modulo``, in ordine.

    Quelli che il modulo importa, direttamente o attraverso altri, esclusi
    quelli di :data:`ESCLUSI` e i moduli delle altre fasi.
    """
    trovati: set[str] = set()
    da_visitare = [modulo]
    while da_visitare:
        for importato in _importati(da_visitare.pop()):
            if (
                importato in trovati
                or importato == modulo
                or importato in ESCLUSI
                or MODULO_DI_FASE.match(importato)
            ):
                continue
            trovati.add(importato)
            da_visitare.append(importato)
    return tuple(sorted(trovati))


def moduli_della_fase(fase: PipelineStep) -> tuple[str, ...]:
    """I moduli del calcolo di una fase: quelli ricavati dagli import
    (:func:`moduli_del_calcolo`), quelli di :data:`DAL_CONTESTO` se la fase
    dipende dalla fase che produce il dato, con i loro import, e quelli di
    :data:`ECCEZIONI_PER_FASE`.
    """
    modulo_fase = type(fase).__module__
    moduli = set(moduli_del_calcolo(modulo_fase))
    dipendenze = {str(d) for d in fase.nodo.dipendenze}
    for modulo, (produttore, _) in DAL_CONTESTO.items():
        if produttore in dipendenze:
            moduli |= {modulo, *moduli_del_calcolo(modulo)}
    moduli |= set(ECCEZIONI_PER_FASE.get(str(fase.passo), {}))
    moduli.discard(modulo_fase)
    return tuple(sorted(moduli))


def file_del_sorgente(fase: PipelineStep, cartella_r: Path) -> dict[str, Path]:
    """I file del sorgente di una fase, per nome logico.

    Il nome logico non dipende da dove stanno i file (``python/<modulo>``,
    ``R/<script>``, ``R/lib/<file>``): la stessa copia in due posti ha la
    stessa impronta.
    """
    file: dict[str, Path] = {}
    for modulo in (type(fase).__module__, *moduli_della_fase(fase)):
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


class GitNonLeggibile(Exception):
    """git non ha risposto sul repository: il messaggio dice perche'."""


def _git(*argomenti: str) -> str:
    """L'uscita di un comando git di sola lettura sul repository del codice.

    Il comando non dipende dall'utente che lo lancia ne' dalla sua
    configurazione. In un container il repository montato appartiene di norma
    a un utente diverso da quello del processo, e git lo rifiuterebbe come
    proprieta' dubbia: la cartella e' dichiarata sicura per questo solo
    comando, che legge e non esegue nulla del repository. La configurazione
    globale e di sistema non si legge: una cartella personale assente o non
    leggibile (l'utente di chi lancia, in un'immagine che ne prevede un altro)
    non deve fermare la lettura, e nessuna impostazione locale deve cambiarne
    l'esito.

    Solleva :class:`GitNonLeggibile` con la causa: git assente, la cartella
    non e' un repository, o l'errore che git ha scritto.
    """
    comando = [
        "git", "-c", f"safe.directory={RADICE_REPOSITORY}", "-C", str(RADICE_REPOSITORY),
        *argomenti,
    ]
    ambiente = {**os.environ, "GIT_CONFIG_GLOBAL": os.devnull, "GIT_CONFIG_SYSTEM": os.devnull,
                "GIT_OPTIONAL_LOCKS": "0", "LC_ALL": "C"}
    try:
        esito = subprocess.run(
            comando, capture_output=True, text=True, timeout=30, check=False, env=ambiente,
        )
    except FileNotFoundError as e:
        raise GitNonLeggibile("il comando git non e' installato") from e
    except (OSError, subprocess.SubprocessError) as e:
        raise GitNonLeggibile(f"git non si avvia: {e}") from e
    if esito.returncode != 0:
        messaggio = " ".join(esito.stderr.split()) or f"uscita {esito.returncode}"
        raise GitNonLeggibile(messaggio)
    return esito.stdout


@functools.cache
def _stato_git() -> tuple[str | None, bool | None]:
    """Commit del repository e presenza di modifiche non committate.

    ``(None, None)`` se il codice non sta in un repository git leggibile, come
    nell'immagine, dove ``.git`` non viene copiato: la provenienza lo dice
    invece di inventare un commit.
    """
    try:
        commit = _git("rev-parse", "HEAD").strip()
        modifiche = _git("status", "--porcelain", "--", "src", "R").strip()
    except GitNonLeggibile:
        return None, None
    return commit or None, bool(modifiche)


def modifiche_non_committate() -> list[str] | None:
    """I file di ``src/`` e ``R/`` con modifiche non committate, o ``None`` se
    git non è leggibile.
    """
    try:
        uscita = _git("status", "--porcelain", "--", "src", "R")
    except GitNonLeggibile:
        return None
    return [riga[3:] for riga in uscita.splitlines() if riga.strip()]


def verifica_git(script_r: Path) -> tuple[str | None, str]:
    """Il codice del catalogo che la regola rigorosa solleva per lo stato di
    git, con il dettaglio; ``(None, commit)`` se il repository è leggibile e
    il codice non ha modifiche non committate.

    Un repository non leggibile è un rifiuto (``E-PROV-01``), non un via
    libera: senza commit nulla identifica il codice eseguito. Lo stesso vale
    se gli script R eseguiti (``script_r``) non sono quelli del repository:
    il commit descriverebbe file diversi da quelli che calcolano.
    """
    # Lo stato si legge adesso, non dalla memoria di :func:`_stato_git`: la
    # regola decide su cio' che c'e' al momento dell'avvio.
    if not (RADICE_REPOSITORY / ".git").exists():
        # Il caso piu' comune, e va detto per quello che e': il codice in
        # esecuzione e' una copia installata, non un clone.
        return "E-PROV-01", (
            f"il codice in esecuzione sta in {Path(__file__).resolve().parents[1]}, che non "
            "fa parte di un clone del repository (nessuna cartella .git in "
            f"{RADICE_REPOSITORY}): e' il pacchetto installato, per esempio quello "
            "dell'immagine. Dalla radice di un clone lancia la pipeline con "
            f"python3 {LANCIATORE}, che esegue il codice del clone"
        )
    try:
        commit = _git("rev-parse", "HEAD").strip()
        modificati = [
            riga[3:] for riga in _git("status", "--porcelain", "--", "src", "R").splitlines()
            if riga.strip()
        ]
    except GitNonLeggibile as guasto:
        return "E-PROV-01", f"git non legge il repository in {RADICE_REPOSITORY}: {guasto}"
    if not commit:
        return "E-PROV-01", f"il repository in {RADICE_REPOSITORY} non ha alcun commit"
    if script_r.resolve() != (RADICE_REPOSITORY / "R").resolve():
        return "E-PROV-01", (
            f"gli script R eseguiti stanno in {script_r}, fuori dal repository "
            f"{RADICE_REPOSITORY} a cui il commit si riferisce: la variabile "
            f"AMPLICON16S_R_DIR indica un'altra cartella (nell'immagine, gli script "
            f"copiati quando e' stata costruita). Lancia la pipeline con python3 "
            f"{LANCIATORE}, che usa gli script del clone"
        )
    if modificati:
        return "E-PROV-02", f"{len(modificati)} file modificati: {', '.join(modificati[:10])}"
    return None, commit


def blocco_r(lockfile: str) -> tuple[Path, dict[str, Any]] | None:
    """Il file di blocco dei pacchetti R e il suo contenuto, se si trova.

    ``run.lockfile`` è cercato com'è (dalla cartella di lavoro, come
    nell'immagine) e poi nella radice del repository.
    """
    for percorso in (Path(lockfile), RADICE_REPOSITORY / lockfile):
        try:
            contenuto = json.loads(percorso.read_bytes())
            contenuto["Packages"]["dada2"]
        except (OSError, ValueError, KeyError, TypeError):
            continue
        return percorso, contenuto
    return None


def _numeri(versione: str) -> list[str]:
    """Le componenti di una versione: R tratta ``-`` e ``.`` come separatori
    equivalenti (``4.7-1.2`` e ``4.7.1.2`` sono la stessa versione).
    """
    return re.split(r"[.-]", versione)


def discordanze_ambiente(blocco: Mapping[str, Any], ambiente: Mapping[str, Any]) -> list[str]:
    """In che cosa l'ambiente R installato differisce dal file di blocco.

    ``ambiente`` è ciò che ``R/00_ambiente.R`` ha letto dalle librerie
    installate: la versione di R e, per pacchetto, versione e correzione.
    Elenco vuoto se corrispondono.
    """
    discordanze = []
    if ambiente.get("r") != blocco["R"]["Version"]:
        discordanze.append(f"R: atteso {blocco['R']['Version']}, trovato {ambiente.get('r')}")
    installati = ambiente.get("pacchetti", {})
    for nome, voce in sorted(blocco["Packages"].items()):
        trovato = installati.get(nome)
        if trovato is None:
            discordanze.append(f"{nome}: non installato")
            continue
        if _numeri(trovato["versione"]) != _numeri(voce["Version"]):
            discordanze.append(f"{nome}: atteso {voce['Version']}, trovato {trovato['versione']}")
        if (voce.get("Patch") or None) != (trovato.get("correzione") or None):
            discordanze.append(f"{nome}: correzione diversa da quella del file di blocco")
    return discordanze


def impronta_rigorosa(fase: PipelineStep, config: Config) -> dict[str, Any]:
    """Ciò che la regola rigorosa aggiunge all'impronta di una fase: il sorgente
    eseguito, il file di blocco dell'ambiente R e l'immagine dichiarata.

    Con la regola attiva, cambiare il codice di una fase, l'ambiente o
    l'immagine dichiarata rende la fase da rifare, e non solo da segnalare.
    """
    from amplicon16s.rbridge.runner import cartella_r

    sorgente, _ = impronta_sorgente(file_del_sorgente(fase, cartella_r()))
    trovato = blocco_r(config.run.lockfile)
    return {
        "sorgente": sorgente,
        "blocco_r": None if trovato is None else (
            "sha256:" + hashlib.sha256(trovato[0].read_bytes()).hexdigest()
        ),
        "immagine": config.run.container,
    }


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
