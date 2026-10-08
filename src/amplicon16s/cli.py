"""Riga di comando della pipeline 16S.

Quattro sottocomandi, tutti con ``--config`` che indica il file di
configurazione. Il file non viene mai modificato: le azioni correttive dei
tentativi ripetuti cambiano una copia in memoria. I percorsi relativi del file
valgono rispetto alla cartella da cui il comando e' lanciato, e sono resi
assoluti al caricamento.

* ``validate`` esegue solo S0. Se S0 è già conclusa e valida per questa
  configurazione non la riesegue: la dichiara conclusa.
* ``run`` esegue dall'inizio, in una cartella di output assente o vuota. **Non
  sovrascrive**: ogni esecuzione deve restare ispezionabile, quindi se
  ``io.out_root`` contiene già qualcosa ``run`` si rifiuta e indica le due
  strade: ``resume`` per continuare quell'esecuzione, oppure un'altra
  ``io.out_root`` per cominciarne una nuova. Non esiste un'opzione per
  sovrascrivere: contraddirebbe proprio quella decisione.
* ``resume`` riprende dagli artefatti esistenti: esegue le sole fasi che la
  valutazione dello stato dà da eseguire.
* ``report`` genera il report dell'esecuzione: un documento HTML e le sue
  tabelle, nella cartella ``report`` sotto ``io.out_root``, ricavati da ciò che
  l'esecuzione ha lasciato su disco, senza rieseguire nulla
  (:mod:`amplicon16s.report.builder`).

**Codici di uscita**, per chi lancia la pipeline da uno script o da uno
scheduler:

* ``0``: successo;
* ``1``: errore imprevisto, un difetto del programma: il log in ``99_logs``
  ne riporta la traccia;
* ``2``: riga di comando non valida (argomenti mancanti o sconosciuti);
* ``3``: errore di configurazione: il file non è valido, G15 lo respinge,
  oppure ``run`` trova la cartella di output già usata. Nessuna fase è
  partita. Il rifiuto di G15 riporta il codice del catalogo di ogni
  violazione (``E-G15-*``) e la sua azione;
* ``4``: arresto con punto di ripresa dichiarato, stampato e scritto in
  ``99_logs/punto_di_ripresa.json`` e ``.txt``;
* ``5``: una fase prevista dal grafo non esiste come codice. Con le quindici
  fasi realizzate non può presentarsi in un'esecuzione normale: resta come
  guardia, perché un grafo esteso prima della fase che lo realizza fermi
  l'esecuzione con un codice proprio invece di un errore imprevisto.
"""

from __future__ import annotations

import argparse
import logging
import os
import shlex
import sys
from collections.abc import Callable, Sequence
from pathlib import Path
from typing import Final

from amplicon16s import __version__
from amplicon16s.config.schema import Config, ErroreConfigurazione, carica
from amplicon16s.errors.catalog import voce
from amplicon16s.gates.g01_g15 import ErroreGate, rifiuto_di_g15
from amplicon16s.logging.logger import chiudi, configura, ottieni
from amplicon16s.runner.executor import Conclusione, Esecutore, EsitoEsecuzione
from amplicon16s.runner.graph import Passo
from amplicon16s.report.builder import CARTELLA_TABELLE, costruisci, scrivi
from amplicon16s.runner.project import ProjectRun, passi_realizzati
from amplicon16s.steps.base import PipelineStep

__all__ = [
    "USCITA_ARRESTO",
    "USCITA_CONFIGURAZIONE",
    "USCITA_ERRORE_IMPREVISTO",
    "USCITA_FASE_NON_REALIZZATA",
    "USCITA_SUCCESSO",
    "build_parser",
    "main",
]

USCITA_SUCCESSO: Final = 0
USCITA_ERRORE_IMPREVISTO: Final = 1
# 2 e' il codice con cui argparse segnala una riga di comando non valida.
USCITA_CONFIGURAZIONE: Final = 3
USCITA_ARRESTO: Final = 4
USCITA_FASE_NON_REALIZZATA: Final = 5

#: I parametri che sono percorsi di file o cartelle, per gruppo.
_PERCORSI: Final = {
    "io": ("fastq_dir", "assay_table", "study_table", "out_root", "batch_table"),
    "tax": ("ref_fasta", "ref_bad_taxa"),
}


def _percorsi_assoluti(config: Config) -> Config:
    """La configurazione con i percorsi relativi resi assoluti rispetto alla
    cartella da cui il comando e' lanciato.

    Un file di configurazione puo' indicare percorsi relativi, come quello di
    ``dati/osd734/``, pensato per essere usato dalla radice del repository. I
    processi R delle fasi partono pero' nella cartella della propria fase: un
    percorso relativo vi indicherebbe un altro file, e la prima fase con un
    calcolo in R si fermerebbe senza trovare i suoi ingressi. Il percorso si
    fissa quindi qui, una volta, e la configurazione registrata in
    ``00_config`` riporta quello assoluto, cioe' i file davvero letti. I
    parametri impostati nel file restano gli stessi: si aggiorna la copia, non
    si ricostruisce la configurazione.
    """
    gruppi = {}
    for gruppo, chiavi in _PERCORSI.items():
        modello = getattr(config, gruppo)
        # La tilde si espande prima: senza, "~/dati" sarebbe un percorso
        # relativo, e nascerebbe una cartella di nome "~" sotto quella corrente.
        relativi = {
            chiave: assoluto
            for chiave in chiavi
            if (valore := getattr(modello, chiave)) is not None
            and (assoluto := (
                valore.expanduser() if valore.expanduser().is_absolute()
                else Path.cwd() / valore
            )) != valore
        }
        if relativi:
            gruppi[gruppo] = modello.model_copy(update=relativi)
    return config.model_copy(update=gruppi) if gruppi else config


def _uscita_non_utilizzabile(radice: Path) -> str | None:
    """Perche' la cartella di output non si puo' usare, o ``None``.

    Si verifica prima di aprire il log, che e' la prima scrittura: una cartella
    che non e' una cartella, o che l'utente non puo' creare o scrivere (un
    volume di un altro utente montato nel container), e' un errore di chi
    configura, con il suo codice di uscita, non un errore imprevisto.
    """
    if radice.exists() and not radice.is_dir():
        return f"{radice} esiste e non e' una cartella"
    esistente = radice
    while not esistente.exists():
        esistente = esistente.parent
    if not esistente.is_dir():
        return f"{esistente} non e' una cartella: {radice} non si puo' creare"
    if not os.access(esistente, os.W_OK | os.X_OK):
        return (f"{esistente} non e' scrivibile dall'utente che esegue: "
                f"{radice} non si puo' creare ne' scrivere")
    return None


def _stampa_rifiuto(rifiuto: ErroreGate) -> None:
    """Il rifiuto di G15 sull'uscita: ogni violazione con il suo codice del
    catalogo, poi l'azione di ogni codice, una volta sola.

    Vale per i due modi in cui G15 respinge: una configurazione che lo schema
    non accetta, e una che lo schema accetta ma i controlli di coerenza no.
    """
    _stampa(str(rifiuto))
    for codice in dict.fromkeys(v.codice for v in rifiuto.violazioni):
        _stampa(f"  cosa fare [{codice}]: {voce(codice).azione}")


def _passi() -> dict[Passo, PipelineStep]:
    """Le fasi realizzate. Un punto solo, cosi' i test possono sostituirle."""
    return passi_realizzati()


def _comando_ripresa(percorso: Path, programma: str = "amplicon16s") -> str:
    """Il comando da stampare per riprendere l'esecuzione con questa
    configurazione, con il programma con cui e' stata lanciata.
    """
    return f"{programma} resume --config {shlex.quote(str(percorso.resolve()))}"


def _stampa(testo: str = "") -> None:
    """Scrive una riga sullo standard output senza bufferizzarla."""
    print(testo, flush=True)


# --------------------------------------------------------------------------- #
# Esecuzione                                                                   #
# --------------------------------------------------------------------------- #


def _esegui(
    config: Config, percorso: Path, *, fino_a: Passo | None = None,
    programma: str = "amplicon16s",
) -> EsitoEsecuzione:
    """Esegue il grafo fino alla fase indicata, o fino all'ultima realizzata."""
    configura(config.io.out_root)
    run = ProjectRun(config, passi=_passi(), logger=ottieni("run"))
    esecutore = Esecutore(
        run, comando_ripresa=_comando_ripresa(percorso, programma), fino_a=fino_a)
    return esecutore.esegui()


def _uscita(esito: EsitoEsecuzione) -> int:
    """Riporta l'esito dell'esecuzione e lo traduce nel codice di uscita del comando."""
    registrazione = esito.configurazione
    if registrazione is not None and registrazione.nuova and registrazione.differenze:
        _stampa(
            f"Configurazione diversa dall'ultima registrata: registrata in "
            f"{registrazione.percorso.name}, accanto alle precedenti "
            f"({'; '.join(registrazione.differenze)})."
        )
    eseguite = ", ".join(str(r.passo) for r in esito.eseguite) or "nessuna"
    if esito.conclusione is Conclusione.COMPLETATA:
        _stampa(f"Esecuzione completata. Fasi eseguite ora: {eseguite}.")
        return USCITA_SUCCESSO
    if esito.conclusione is Conclusione.FASE_NON_REALIZZATA:
        _stampa(
            f"Fasi eseguite ora: {eseguite}. Tutte le fasi realizzate sono concluse; "
            f"la prossima, {esito.non_realizzata}, non e' ancora realizzata."
        )
        return USCITA_FASE_NON_REALIZZATA
    assert esito.punto is not None
    _stampa(esito.punto.testo())
    return USCITA_ARRESTO


def _cmd_validate(config: Config, percorso: Path, args: argparse.Namespace) -> int:
    """Comando ``validate``: esegue la sola S0 e ne riporta l'esito."""
    esito = _esegui(config, percorso, fino_a=Passo.S0, programma=args.programma)
    if esito.conclusione is Conclusione.COMPLETATA:
        registrazione = esito.configurazione
        if registrazione is not None and registrazione.nuova and registrazione.differenze:
            _stampa(
                f"Configurazione diversa dall'ultima registrata: registrata in "
                f"{registrazione.percorso.name}."
            )
        if esito.eseguite:
            _stampa("Validazione superata: S0 conclusa.")
        else:
            _stampa("S0 era gia' conclusa e valida per questa configurazione.")
        return USCITA_SUCCESSO
    return _uscita(esito)


def _cartella_usata(radice: Path) -> bool:
    """Vero se la cartella di output esiste e contiene già almeno un elemento."""
    return radice.exists() and any(radice.iterdir())


def _cmd_run(config: Config, percorso: Path, args: argparse.Namespace) -> int:
    """Comando ``run``: avvia un'esecuzione nuova, solo in una cartella di output vuota.
    """
    radice = Path(config.io.out_root)
    if _cartella_usata(radice):
        _stampa(
            f"La cartella di output {radice} non e' vuota: contiene gia' "
            "un'esecuzione, o altri file. 'run' non sovrascrive, perche' ogni "
            "esecuzione deve restare ispezionabile.\n"
            f"  per continuare quell'esecuzione:  "
            f"{_comando_ripresa(percorso, args.programma)}\n"
            "  per cominciarne una nuova:        indica in io.out_root una cartella "
            "nuova o vuota"
        )
        return USCITA_CONFIGURAZIONE
    return _uscita(_esegui(config, percorso, programma=args.programma))


def _cmd_resume(config: Config, percorso: Path, args: argparse.Namespace) -> int:
    """Comando ``resume``: riprende l'esecuzione dalla prima fase non valida."""
    return _uscita(_esegui(config, percorso, programma=args.programma))


def _cmd_report(config: Config, percorso: Path, args: argparse.Namespace) -> int:
    """Comando ``report``: genera il report dell'esecuzione in ``io.out_root``.

    Della configurazione usa solo ``io.out_root``: il report si ricava da ciò
    che l'esecuzione ha lasciato nella cartella, e non scrive nel log.
    """
    radice = Path(config.io.out_root)
    if not _cartella_usata(radice):
        _stampa(f"Nessuna esecuzione in {radice}: non c'e' nulla da riportare.")
        return USCITA_CONFIGURAZIONE
    report = costruisci(radice, _passi())
    destinazione = scrivi(report, radice)
    _stampa(f"Report scritto in {destinazione}")
    _stampa(f"Tabelle: {len(report.tabelle)} file in {destinazione.parent / CARTELLA_TABELLE}")
    _stampa("Segnalazioni in apertura:" if report.segnalazioni else "Nessuna segnalazione in apertura.")
    for segnalazione in report.segnalazioni:
        _stampa(f"  - {segnalazione}")
    return USCITA_SUCCESSO


# --------------------------------------------------------------------------- #
# Parser e punto di ingresso                                                   #
# --------------------------------------------------------------------------- #

_Comando = Callable[[Config, Path, argparse.Namespace], int]


def build_parser() -> argparse.ArgumentParser:
    """Costruisce il parser della riga di comando."""
    parser = argparse.ArgumentParser(
        prog="amplicon16s",
        description="Pipeline per l'analisi di dati 16S rRNA da sequenziamento Illumina single-end.",
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")

    subparsers = parser.add_subparsers(dest="command", metavar="comando", required=True)

    comandi: list[tuple[str, str, _Comando]] = [
        ("run", "esegue la pipeline dall'inizio, in una cartella di output nuova", _cmd_run),
        ("resume", "riprende un'esecuzione dagli artefatti gia' prodotti", _cmd_resume),
        ("validate", "esegue solo la fase di validazione degli input (S0)", _cmd_validate),
        ("report", "genera il report HTML di un'esecuzione, dalla sua cartella di output", _cmd_report),
    ]
    for nome, aiuto, funzione in comandi:
        sotto = subparsers.add_parser(nome, help=aiuto)
        sotto.add_argument(
            "--config", "-c", required=True, type=Path,
            help="file di configurazione YAML (non viene mai modificato)",
        )
        sotto.set_defaults(func=funzione)

    return parser


def main(argv: Sequence[str] | None = None, programma: str = "amplicon16s") -> int:
    """Punto di ingresso della riga di comando.

    ``programma`` e' il comando con cui la pipeline e' stata lanciata: compare
    nel comando di ripresa che un arresto dichiara, perche' chi ha lanciato il
    codice di un clone (``scripts/esegui.py``) riprenda con lo stesso.
    """
    parser = build_parser()
    args = parser.parse_args(argv)
    args.programma = programma

    try:
        config = _percorsi_assoluti(carica(args.config))
    except ErroreConfigurazione as e:
        # Una configurazione che lo schema non accetta e' respinta da G15, il
        # gate della configurazione: ogni problema porta il codice del suo
        # controllo, e i parametri obbligatori non dichiarati sono elencati
        # tutti insieme (E-G15-10), con l'azione del catalogo.
        _stampa(f"configurazione non valida ({e.origine})" if e.origine else
                "configurazione non valida")
        _stampa_rifiuto(rifiuto_di_g15(e))
        return USCITA_CONFIGURAZIONE
    except OSError as e:
        _stampa(f"configurazione non leggibile: {e}")
        return USCITA_CONFIGURAZIONE

    problema = _uscita_non_utilizzabile(Path(config.io.out_root))
    if problema is not None:
        _stampa(f"configurazione non valida: io.out_root: {problema}")
        return USCITA_CONFIGURAZIONE

    try:
        return args.func(config, args.config, args)
    except ErroreGate as e:
        # Solo G15 arriva fin qui: gli altri gate sono dentro S0, e un loro
        # fallimento e' un arresto con punto di ripresa.
        _stampa_rifiuto(e)
        return USCITA_CONFIGURAZIONE
    except Exception as e:  # noqa: BLE001 - l'ultimo argine prima dell'uscita
        ottieni().exception("errore imprevisto", extra={"comando": args.command})
        _stampa(f"errore imprevisto: {type(e).__name__}: {e}")
        return USCITA_ERRORE_IMPREVISTO
    finally:
        chiudi()


if __name__ == "__main__":
    logging.captureWarnings(True)
    sys.exit(main())
