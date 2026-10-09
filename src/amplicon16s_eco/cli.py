"""Riga di comando delle analisi ecologiche di base.

Due sottocomandi::

    python3 -m amplicon16s_eco validate --object ps_final.rds --config eco.yaml
    python3 -m amplicon16s_eco run --object ps_final.rds --config eco.yaml --out CARTELLA

* ``validate`` verifica la configurazione, da sola e contro l'oggetto (colonne,
  gruppi, rango, albero, disegno dei test), senza calcolare nulla.
* ``run`` ripete la validazione ed esegue le analisi di base: alfa diversita',
  composizione, distanze, ordinazione e test fra gruppi. La cartella di uscita deve
  essere assente o vuota e fuori dalla cartella dell'oggetto: l'oggetto non
  viene mai modificato.

**Codici di uscita**:

* ``0``: successo (gli avvisi sono stampati con il loro codice);
* ``1``: errore imprevisto, un difetto del programma;
* ``2``: riga di comando non valida;
* ``3``: rifiuto con uno o piu' codici ``E-ECO-``: ogni rifiuto dice la causa e
  che cosa correggere. Nessuna uscita e' stata scritta.
"""

from __future__ import annotations

import argparse
import sys
from collections.abc import Sequence
from pathlib import Path
from typing import Final

from amplicon16s_eco import __version__
from amplicon16s_eco.catalogo import ErroreEco, descrivi
from amplicon16s_eco.config import carica
from amplicon16s_eco.esecuzione import Esito, esegui, valida_sull_oggetto

__all__ = ["USCITA_RIFIUTO", "USCITA_SUCCESSO", "build_parser", "main"]

USCITA_SUCCESSO: Final = 0
USCITA_ERRORE_IMPREVISTO: Final = 1
# 2 e' il codice con cui argparse segnala una riga di comando non valida.
USCITA_RIFIUTO: Final = 3


def build_parser(programma: str = "python3 -m amplicon16s_eco") -> argparse.ArgumentParser:
    """Il parser dei due sottocomandi."""
    parser = argparse.ArgumentParser(
        prog=programma,
        description="Analisi ecologiche di base sull'oggetto finale della pipeline 16S.",
    )
    parser.add_argument("--version", action="version", version=f"amplicon16s_eco {__version__}")
    sotto = parser.add_subparsers(dest="comando", required=True, metavar="COMANDO")
    for nome, aiuto in (
        ("validate", "verifica configurazione e oggetto, senza calcolare"),
        ("run", "valida ed esegue le analisi"),
    ):
        comando = sotto.add_parser(nome, help=aiuto)
        comando.add_argument("--object", dest="oggetto", required=True, type=Path,
                             metavar="FILE", help="l'oggetto phyloseq finale (ps_final.rds)")
        comando.add_argument("--config", required=True, type=Path, metavar="FILE",
                             help="la configurazione dell'analisi")
        if nome == "run":
            comando.add_argument("--out", dest="uscita", required=True, type=Path,
                                 metavar="CARTELLA",
                                 help="cartella di uscita, assente o vuota")
    return parser


def _stampa_esito(esito: Esito) -> None:
    campioni = esito.riepilogo.get("campioni", {})
    alfa = esito.riepilogo.get("alfa", {})
    print(f"oggetto di partenza: {esito.oggetto_sha256}")
    print(f"configurazione: digest {esito.digest_configurazione}")
    print(f"campioni: {campioni.get('analizzati')} analizzati su "
          f"{campioni.get('nell_oggetto')} dell'oggetto")
    print(f"profondita' di rarefazione: {alfa.get('profondita_di_rarefazione')} "
          f"({alfa.get('origine')}), {alfa.get('campioni_rarefatti')} campioni")
    for avviso in esito.avvisi:
        print(f"avviso: {descrivi(avviso.codice, avviso.dettaglio)}")


def main(argv: Sequence[str] | None = None, programma: str = "python3 -m amplicon16s_eco") -> int:
    """Punto d'ingresso: restituisce il codice di uscita."""
    args = build_parser(programma).parse_args(argv)
    try:
        config = carica(args.config)
        if args.comando == "validate":
            esito = valida_sull_oggetto(args.oggetto, config)
            _stampa_esito(esito)
            print("validazione superata")
        else:
            esito = esegui(args.oggetto, config, args.uscita)
            _stampa_esito(esito)
            print(f"file prodotti: {len(esito.file)} in {args.uscita}")
    except ErroreEco as e:
        print(str(e), file=sys.stderr)
        return USCITA_RIFIUTO
    return USCITA_SUCCESSO
