"""Riga di comando della pipeline 16S.

Al momento i quattro sottocomandi sono segnaposto: dichiarano l'interfaccia e
stampano quale comando è stato invocato, senza alcuna logica applicativa.
"""

from __future__ import annotations

import argparse
import sys
from collections.abc import Sequence

from amplicon16s import __version__

_NOT_IMPLEMENTED = "non ancora implementato"


def _cmd_run(args: argparse.Namespace) -> int:
    print(f"[amplicon16s] comando 'run': esegue la pipeline completa — {_NOT_IMPLEMENTED}")
    return 0


def _cmd_resume(args: argparse.Namespace) -> int:
    print(
        f"[amplicon16s] comando 'resume': riprende un'esecuzione interrotta dagli "
        f"artefatti già prodotti — {_NOT_IMPLEMENTED}"
    )
    return 0


def _cmd_validate(args: argparse.Namespace) -> int:
    print(
        f"[amplicon16s] comando 'validate': esegue solo la fase di validazione degli "
        f"input — {_NOT_IMPLEMENTED}"
    )
    return 0


def _cmd_report(args: argparse.Namespace) -> int:
    print(
        f"[amplicon16s] comando 'report': rigenera il report da un'esecuzione già "
        f"conclusa — {_NOT_IMPLEMENTED}"
    )
    return 0


def build_parser() -> argparse.ArgumentParser:
    """Costruisce il parser della riga di comando."""
    parser = argparse.ArgumentParser(
        prog="amplicon16s",
        description="Pipeline per l'analisi di dati 16S rRNA da sequenziamento Illumina single-end.",
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")

    subparsers = parser.add_subparsers(dest="command", metavar="comando", required=True)

    p_run = subparsers.add_parser("run", help="esegue la pipeline completa")
    p_run.set_defaults(func=_cmd_run)

    p_resume = subparsers.add_parser(
        "resume", help="riprende un'esecuzione interrotta dagli artefatti già prodotti"
    )
    p_resume.set_defaults(func=_cmd_resume)

    p_validate = subparsers.add_parser(
        "validate", help="esegue solo la fase di validazione degli input"
    )
    p_validate.set_defaults(func=_cmd_validate)

    p_report = subparsers.add_parser(
        "report", help="rigenera il report da un'esecuzione già conclusa"
    )
    p_report.set_defaults(func=_cmd_report)

    return parser


def main(argv: Sequence[str] | None = None) -> int:
    """Punto di ingresso della riga di comando."""
    parser = build_parser()
    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
