"""Riga di comando delle analisi ecologiche di base.

I sottocomandi::

    python3 -m amplicon16s_eco init --object ps_final.rds --output eco.yaml
    python3 -m amplicon16s_eco validate --object ps_final.rds --config eco.yaml [--json]
    python3 -m amplicon16s_eco run --object ps_final.rds --config eco.yaml --out CARTELLA
    python3 -m amplicon16s_eco status --out CARTELLA [--json]
    python3 -m amplicon16s_eco schema

* ``init`` esamina l'oggetto e scrive una configurazione di partenza: cio' che
  rileva, cio' che propone e cio' che resta da scegliere
  (:mod:`amplicon16s_eco.init`). Non scrive su un file esistente.
* ``schema`` esporta in JSON lo schema della configurazione; ``status`` dice a
  che punto e' l'analisi di una cartella di uscita; con ``--json`` ``validate``
  e ``status`` scrivono sullo standard output un solo documento JSON
  (:mod:`amplicon16s_eco.interfaccia`).

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
import json
import sys
from collections.abc import Sequence
from pathlib import Path
from typing import Final

from amplicon16s_eco import __version__
from amplicon16s_eco.catalogo import ErroreEco, descrivi
from amplicon16s_eco import interfaccia
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
        else:
            comando.add_argument("--json", action="store_true",
                                 help="scrive sullo standard output l'esito come documento JSON")
    inizio = sotto.add_parser("init", help="scrive una configurazione di partenza dall'oggetto")
    inizio.add_argument("--object", dest="oggetto", required=True, type=Path, metavar="FILE",
                        help="l'oggetto phyloseq finale (ps_final.rds)")
    inizio.add_argument("--output", required=True, type=Path, metavar="FILE",
                        help="il file di configurazione da scrivere: non deve esistere")
    inizio.add_argument("--accetta-proposte", action="store_true",
                        help="scrive come righe attive le proposte e i valori di partenza "
                             "(non i condizionali): restano da scrivere le scelte")
    sotto.add_parser("schema", help="esporta in JSON lo schema della configurazione")
    stato = sotto.add_parser("status", help="dice a che punto e' l'analisi di una cartella di uscita")
    stato.add_argument("--out", dest="uscita", required=True, type=Path, metavar="CARTELLA",
                       help="la cartella di uscita dell'analisi")
    stato.add_argument("--json", action="store_true",
                       help="scrive sullo standard output lo stato come documento JSON")
    return parser


def _json(documento: dict) -> None:
    print(json.dumps(documento, ensure_ascii=False, indent=2), flush=True)


_STATI: Final = {
    "non_avviata": "nessuna analisi in questa cartella",
    "in_corso": "analisi in corso",
    "interrotta": "analisi interrotta senza dichiarazione (il processo non risponde)",
    "conclusa": "analisi conclusa",
}


def _cmd_init(args: argparse.Namespace) -> int:
    """Comando ``init``: scrive la configurazione di partenza e il riepilogo."""
    from amplicon16s.init.scrittura import SENZA_VALORE, Specie, valore_yaml

    from amplicon16s_eco.init import scrivi

    if args.output.exists():
        print(f"init: {args.output} esiste gia': init non scrive su una configurazione "
              "esistente. Indica con --output un file nuovo", file=sys.stderr)
        return USCITA_RIFIUTO
    accettate = bool(args.accetta_proposte)
    voci, fatti = scrivi(args.oggetto, args.output, accettate)
    print(f"Configurazione di partenza scritta in {args.output}")
    print("Rilevato sull'oggetto:")
    for fatto in fatti:
        print(f"  {fatto}")
    print("Proposti (righe attive, accettati con --accetta-proposte: controllali):"
          if accettate else "Proposti (righe in commento: togli il commento per accettare):")
    for voce in voci.values():
        if voce.specie is Specie.PROPOSTO:
            print(f"  {voce.chiave}: {valore_yaml(voce.valore)}"
                  + ("" if voce.attiva or not accettate else " (condizionale: in commento)"))
    partenze = [voce for voce in voci.values() if voce.specie is Specie.PARTENZA]
    if partenze:
        print("Valori di partenza accettati (righe attive):")
        for voce in partenze:
            print(f"  {voce.chiave}: {valore_yaml(voce.valore)}")
    print("Restano da scegliere (righe in commento):" if accettate else "Da scegliere:")
    for voce in voci.values():
        if voce.specie is Specie.DA_SCEGLIERE:
            partenza = "" if voce.valore is SENZA_VALORE else f" (partenza: {valore_yaml(voce.valore)})"
            print(f"  {voce.chiave}{partenza}")
    print(f"Poi: rivedi {args.output} e lancia eco validate, che nomina i parametri mancanti.")
    return USCITA_SUCCESSO


def _cmd_status(args: argparse.Namespace) -> int:
    """Comando ``status``: lo stato dell'analisi, letto dalla cartella."""
    documento = interfaccia.stato(args.uscita)
    if args.json:
        _json(documento)
        return USCITA_SUCCESSO
    print(f"{args.uscita}: {_STATI[documento['stato']]}")
    if documento["esecuzione"] is not None:
        print(f"  ultimo segnale {documento['esecuzione']['secondi_dall_ultimo_battito']} s fa")
    if documento["stato"] == "conclusa":
        print(f"  file prodotti: {documento['file']}; avvisi: "
              f"{', '.join(documento['avvisi']) or 'nessuno'}")
    return USCITA_SUCCESSO


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
    if args.comando == "schema":
        _json(interfaccia.schema())
        return USCITA_SUCCESSO
    if args.comando == "status":
        return _cmd_status(args)
    in_json = bool(getattr(args, "json", False))
    try:
        if args.comando == "init":
            return _cmd_init(args)
        config = carica(args.config)
        if args.comando == "validate":
            esito = valida_sull_oggetto(args.oggetto, config)
            if in_json:
                _json(interfaccia.esito_validazione(
                    configurazione=args.config, oggetto=args.oggetto,
                    codice_uscita=USCITA_SUCCESSO, problemi=[],
                    avvisi=[interfaccia.problema(a.codice, a.dettaglio) for a in esito.avvisi]))
                return USCITA_SUCCESSO
            _stampa_esito(esito)
            print("validazione superata")
        else:
            esito = esegui(args.oggetto, config, args.uscita)
            _stampa_esito(esito)
            print(f"file prodotti: {len(esito.file)} in {args.uscita}")
    except ErroreEco as e:
        print(str(e), file=sys.stderr)
        if in_json:
            _json(interfaccia.esito_validazione(
                configurazione=args.config, oggetto=args.oggetto,
                codice_uscita=USCITA_RIFIUTO,
                problemi=[interfaccia.problema(c, d) for c, d in e.rifiuti]))
        return USCITA_RIFIUTO
    return USCITA_SUCCESSO
