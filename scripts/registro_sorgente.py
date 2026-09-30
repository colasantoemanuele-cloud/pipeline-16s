"""Verifica e aggiorna il registro delle versioni e del sorgente delle fasi.

Il registro, ``src/amplicon16s/steps/registro_sorgente.json``, riporta per ogni
fase realizzata la versione dichiarata e l'impronta del suo sorgente (modulo,
script R, file di ``R/lib`` caricati). Un test fallisce se non corrisponde al
codice: chi modifica una fase sceglie qui fra le due strade.

    python scripts/registro_sorgente.py
        verifica: elenca le fasi che non corrispondono, esce con 1 se ce ne sono

    python scripts/registro_sorgente.py --aggiorna
        aggiorna il registro; rifiuta se una fase ha il sorgente cambiato a
        parita' di versione, finche' non si sceglie:
          - se cambia cio' che la fase calcola, si incrementa la sua versione
            (la ripresa rifara' la fase e le successive);
          - se la modifica e' senza effetto sui risultati (commenti,
            riorganizzazione), la si dichiara:

    python scripts/registro_sorgente.py --aggiorna --senza-effetto S4 S6
"""

from __future__ import annotations

import argparse
import sys

from amplicon16s.runner.project import passi_realizzati
from amplicon16s.runner.provenienza import (
    REGISTRO,
    RegistroNonAggiornabile,
    aggiorna_registro,
    differenze_registro,
    leggi_registro,
    scrivi_registro,
)


def main(argv: list[str] | None = None) -> int:
    """Punto d'ingresso da riga di comando: verifica, o aggiorna con ``--aggiorna``."""
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--aggiorna", action="store_true", help="scrive il registro aggiornato")
    parser.add_argument(
        "--senza-effetto", nargs="+", default=[], metavar="FASE",
        help="fasi modificate senza effetto sui risultati, a parita' di versione",
    )
    args = parser.parse_args(argv)
    if args.senza_effetto and not args.aggiorna:
        parser.error("--senza-effetto vale solo con --aggiorna")

    fasi = {str(p): f for p, f in passi_realizzati().items()}
    registro = leggi_registro() if REGISTRO.exists() else {}
    differenze = differenze_registro(registro, fasi)

    if not args.aggiorna:
        for nome, differenza in differenze.items():
            print(f"{nome}: {differenza}")
        print("registro allineato al codice" if not differenze else
              f"{len(differenze)} fasi non corrispondono al registro")
        return 1 if differenze else 0

    try:
        nuovo = aggiorna_registro(registro, fasi, frozenset(args.senza_effetto))
    except RegistroNonAggiornabile as e:
        print(f"registro non aggiornato: {e}", file=sys.stderr)
        return 1
    scrivi_registro(nuovo)
    for nome, differenza in differenze.items():
        print(f"{nome}: {differenza}, registrato")
    print(f"registro scritto in {REGISTRO}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
