"""Riga di comando della pipeline 16S.

Quattro sottocomandi, tutti con ``--config`` che indica il file di
configurazione. Il file non viene mai modificato: le azioni correttive dei
tentativi ripetuti cambiano una copia in memoria.

* ``validate`` esegue solo S0. Se S0 è già conclusa e valida per questa
  configurazione non la riesegue: la dichiara conclusa.
* ``run`` esegue dall'inizio, in una cartella di output assente o vuota. **Non
  sovrascrive**: ogni esecuzione deve restare ispezionabile, quindi se
  ``io.out_root`` contiene già qualcosa ``run`` si rifiuta e indica le due
  strade — ``resume`` per continuare quell'esecuzione, oppure un'altra
  ``io.out_root`` per cominciarne una nuova. Non esiste un'opzione per
  sovrascrivere: contraddirebbe proprio quella decisione.
* ``resume`` riprende dagli artefatti esistenti: esegue le sole fasi che la
  valutazione dello stato dà da eseguire.
* ``report`` scrive un resoconto provvisorio dello stato, ricavato dai
  manifesti delle fasi; il report completo non è ancora realizzato.

**Codici di uscita**, per chi lancia la pipeline da uno script o da uno
scheduler:

* ``0`` — successo;
* ``1`` — errore imprevisto, un difetto del programma: il log in ``99_logs``
  ne riporta la traccia;
* ``2`` — riga di comando non valida (argomenti mancanti o sconosciuti);
* ``3`` — errore di configurazione: il file non è valido, G15 lo respinge,
  oppure ``run`` trova la cartella di output già usata. Nessuna fase è
  partita;
* ``4`` — arresto con punto di ripresa dichiarato, stampato e scritto in
  ``99_logs/punto_di_ripresa.json`` e ``.txt``;
* ``5`` — tutte le fasi realizzate sono concluse, ma la prossima non esiste
  ancora come codice. È uno stato transitorio dello sviluppo: oggi, dopo S0.
"""

from __future__ import annotations

import argparse
import json
import logging
import shlex
import sys
from collections.abc import Callable, Sequence
from pathlib import Path
from typing import Any, Final

from amplicon16s import __version__
from amplicon16s.config.schema import Config, ErroreConfigurazione, carica
from amplicon16s.gates.g01_g15 import ErroreGate
from amplicon16s.io_layer.artifacts import Fase
from amplicon16s.logging.logger import chiudi, configura, ottieni
from amplicon16s.runner.executor import Conclusione, Esecutore, EsitoEsecuzione
from amplicon16s.runner.graph import Passo
from amplicon16s.runner.project import ProjectRun, StatoPasso, passi_realizzati
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

#: Resoconto provvisorio dello stato, sotto 99_logs.
NOME_RESOCONTO: Final = "resoconto_stato.json"


def _passi() -> dict[Passo, PipelineStep]:
    """Le fasi realizzate. Un punto solo, cosi' i test possono sostituirle."""
    return passi_realizzati()


def _comando_ripresa(percorso: Path) -> str:
    return f"amplicon16s resume --config {shlex.quote(str(percorso.resolve()))}"


def _stampa(testo: str = "") -> None:
    print(testo, flush=True)


# --------------------------------------------------------------------------- #
# Esecuzione                                                                   #
# --------------------------------------------------------------------------- #


def _esegui(
    config: Config, percorso: Path, *, fino_a: Passo | None = None
) -> EsitoEsecuzione:
    configura(config.io.out_root)
    run = ProjectRun(config, passi=_passi(), logger=ottieni("run"))
    esecutore = Esecutore(run, comando_ripresa=_comando_ripresa(percorso), fino_a=fino_a)
    return esecutore.esegui()


def _uscita(esito: EsitoEsecuzione) -> int:
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
    esito = _esegui(config, percorso, fino_a=Passo.S0)
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
    return radice.exists() and any(radice.iterdir())


def _cmd_run(config: Config, percorso: Path, args: argparse.Namespace) -> int:
    radice = Path(config.io.out_root)
    if _cartella_usata(radice):
        _stampa(
            f"La cartella di output {radice} non e' vuota: contiene gia' "
            "un'esecuzione, o altri file. 'run' non sovrascrive, perche' ogni "
            "esecuzione deve restare ispezionabile.\n"
            f"  per continuare quell'esecuzione:  {_comando_ripresa(percorso)}\n"
            "  per cominciarne una nuova:        indica in io.out_root una cartella "
            "nuova o vuota"
        )
        return USCITA_CONFIGURAZIONE
    return _uscita(_esegui(config, percorso))


def _cmd_resume(config: Config, percorso: Path, args: argparse.Namespace) -> int:
    return _uscita(_esegui(config, percorso))


def resoconto(run: ProjectRun) -> dict[str, Any]:
    """Il resoconto provvisorio dello stato, dai manifesti delle fasi."""
    valutazione = run.valuta()
    fasi = []
    for passo, situazione in valutazione.situazioni.items():
        voce: dict[str, Any] = {
            "passo": str(passo),
            "descrizione": run.grafo.nodo(passo).descrizione,
            "stato": situazione.stato.value,
            "motivo": situazione.motivo,
        }
        if situazione.stato is StatoPasso.COMPLETATA:
            manifesto = run.albero.manifesto_passo(passo, run.grafo.nodo(passo).cartella)
            assert manifesto is not None
            voce["conclusa"] = manifesto.conclusa
            voce["artefatti"] = list(manifesto.nomi)
            voce["aggiustamenti"] = list(manifesto.aggiustamenti)
            voce["degradazioni"] = [
                {"codice": d["codice"], "dettaglio": d.get("dettaglio")}
                for d in manifesto.degradazioni
            ]
        fasi.append(voce)
    return {
        "provvisorio": True,
        "completa": valutazione.completa,
        "prossima": str(valutazione.prossima()) if valutazione.prossima() else None,
        "fasi": fasi,
    }


def _testo_resoconto(documento: dict[str, Any]) -> str:
    righe = [
        "RESOCONTO PROVVISORIO DELLO STATO (non e' il report definitivo)",
        "",
    ]
    for fase in documento["fasi"]:
        righe.append(f"  {fase['passo']:<4} {fase['stato']:<15} {fase['descrizione']}")
        for a in fase.get("aggiustamenti", []):
            righe.append(
                f"         aggiustamento [{a['codice']}]: {a['parametro']} "
                f"{a['dichiarato']} -> {a['usato']}"
            )
        for d in fase.get("degradazioni", []):
            righe.append(f"         degradazione [{d['codice']}]: {d['dettaglio']}")
    righe += [
        "",
        "Esecuzione completa." if documento["completa"]
        else f"Esecuzione non completa: prossima fase {documento['prossima']}.",
    ]
    return "\n".join(righe)


def _cmd_report(config: Config, percorso: Path, args: argparse.Namespace) -> int:
    radice = Path(config.io.out_root)
    if not _cartella_usata(radice):
        _stampa(f"Nessuna esecuzione in {radice}: non c'e' nulla da riportare.")
        return USCITA_CONFIGURAZIONE
    run = ProjectRun(config, passi=_passi())
    documento = resoconto(run)
    destinazione = run.albero.prepara(Fase.LOGS) / NOME_RESOCONTO
    destinazione.write_text(
        json.dumps(documento, indent=2, ensure_ascii=False, default=str) + "\n",
        encoding="utf-8",
    )
    _stampa(_testo_resoconto(documento))
    _stampa(f"\nResoconto scritto in {destinazione}")
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
        ("report", "resoconto provvisorio dello stato dell'esecuzione", _cmd_report),
    ]
    for nome, aiuto, funzione in comandi:
        sotto = subparsers.add_parser(nome, help=aiuto)
        sotto.add_argument(
            "--config", "-c", required=True, type=Path,
            help="file di configurazione YAML (non viene mai modificato)",
        )
        sotto.set_defaults(func=funzione)

    return parser


def main(argv: Sequence[str] | None = None) -> int:
    """Punto di ingresso della riga di comando."""
    parser = build_parser()
    args = parser.parse_args(argv)

    try:
        config = carica(args.config)
    except ErroreConfigurazione as e:
        _stampa(str(e))
        return USCITA_CONFIGURAZIONE
    except OSError as e:
        _stampa(f"configurazione non leggibile: {e}")
        return USCITA_CONFIGURAZIONE

    try:
        return args.func(config, args.config, args)
    except ErroreGate as e:
        # Solo G15 arriva fin qui: gli altri gate sono dentro S0, e un loro
        # fallimento e' un arresto con punto di ripresa.
        _stampa(str(e))
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
