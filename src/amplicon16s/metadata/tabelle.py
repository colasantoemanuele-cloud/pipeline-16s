"""Lettura delle tabelle dei metadati separate da tabulazioni.

Le tabelle ISA-Tab e il file di arricchimento del lotto sono compilati a mano:
spazi e virgolette ai margini dei valori compaiono in modo incostante, una
colonna senza nome nasce da una tabulazione in piu' in fondo alla riga, e un
foglio di calcolo puo' salvare il file con un segno d'ordine dei byte (BOM) in
testa. Tutta la pipeline le legge da qui (i gate di S0, il crosswalk, l'oggetto
integrato di S10): intestazioni e valori sono ripuliti in un solo modo, e una
colonna dichiarata in configurazione si confronta sempre con lo stesso nome.
"""

from __future__ import annotations

import csv
import re
from collections.abc import Iterable
from pathlib import Path
from typing import TYPE_CHECKING, Final

if TYPE_CHECKING:
    from amplicon16s.config.schema import Config

__all__ = [
    "COLONNE_INVENTARIO",
    "intestazione",
    "leggi_tsv",
    "nome_nell_oggetto",
    "nomi_in_collisione",
    "pulisci",
    "valori_non_tabellari",
    "tabella_delle_cellule",
    "tabella_di_studio",
]

#: UTF-8 con o senza BOM: con "utf-8-sig" il segno, se c'e', non finisce nel
#: nome della prima colonna.
_CODIFICA = "utf-8-sig"


def pulisci(valore: str | None) -> str:
    """Toglie spazi e virgolette: le tabelle ISA le usano in modo incostante."""
    return (valore or "").strip().strip('"').strip()


def intestazione(percorso: Path) -> list[str]:
    """I nomi delle colonne di una tabella, ripuliti come quelli che
    :func:`leggi_tsv` usa per chiave; elenco vuoto se il file non ha righe.
    """
    with open(percorso, encoding=_CODIFICA, newline="") as file:
        return [pulisci(c) for c in next(csv.reader(file, delimiter="\t"), [])]


def leggi_tsv(percorso: Path) -> list[dict[str, str]]:
    """Le righe di una tabella separata da tabulazioni, con nomi delle colonne e
    valori ripuliti. Le colonne senza nome si ignorano.
    """
    with open(percorso, encoding=_CODIFICA, newline="") as file:
        return [
            {nome: pulisci(valore) for chiave, valore in riga.items()
             if chiave and (nome := pulisci(chiave))}
            for riga in csv.DictReader(file, delimiter="\t")
        ]


def valori_non_tabellari(percorso: Path, colonne: Iterable[str]) -> list[tuple[int, str]]:
    """I valori delle colonne indicate che contengono una tabulazione o un a
    capo, come coppie (riga del file, colonna); la prima riga di dati e' la 2.

    Un campo fra virgolette puo' contenerli, e il lettore lo restituisce per
    intero: ma le tabelle che la pipeline scrive (il crosswalk, i metadati
    dell'oggetto) sono senza virgolette, e un valore cosi' le spezzerebbe. Si
    guardano le sole colonne che la pipeline legge: un campo di testo libero
    che non entra in alcun artefatto non e' un difetto.
    """
    attese = set(colonne)
    trovati = []
    for numero, riga in enumerate(leggi_tsv(percorso), start=2):
        for nome, valore in riga.items():
            if nome in attese and any(c in valore for c in "\t\n\r"):
                trovati.append((numero, nome))
    return trovati


def tabella_di_studio(config: Config) -> tuple[Path, str]:
    """La tabella da cui si leggono classe e variabili dei campioni, e la sua
    colonna con il nome del campione.

    E' ``io.study_table`` se indicata; altrimenti la tabella di assay, che
    allora porta anche quelle colonne.
    """
    if config.io.study_table is not None:
        return Path(config.io.study_table), config.meta.colonna_id_studio
    return Path(config.io.assay_table), config.meta.sample_id_column


#: Le colonne dell'inventario nell'oggetto integrato, con i nomi del crosswalk
#: di S0: una colonna richiesta non puo' prenderne il nome.
COLONNE_INVENTARIO: Final = (
    "accession", "sample_name", "classe", "materiale",
    "posizione", "modulo", "piastra", "corsa",
)

#: Un nome che R conserva cosi' com'e': make.names lo lascia invariato, e non
#: e' una parola riservata (nessuna delle riservate ha questa forma minuscola
#: con trattini bassi, salvo quelle escluse sotto).
_SINTATTICO: Final = re.compile(r"^[a-z][a-z0-9_]*$")
_RISERVATE: Final = frozenset({
    "if", "else", "repeat", "while", "function", "for", "next", "break",
    "in", "true", "false", "null", "inf", "nan", "na",
})


def nome_nell_oggetto(originale: str) -> str:
    """Il nome sintattico di una colonna dei metadati nell'oggetto.

    Minuscole, e ogni sequenza di caratteri diversi da lettere e cifre
    diventa un trattino basso: ``Factor Value[Spaceflight]`` diventa
    ``factor_value_spaceflight``. Un nome che non comincia con una lettera
    prende il prefisso ``x_``. Il risultato e' un nome che R non altera.
    """
    nome = re.sub(r"[^a-z0-9]+", "_", originale.casefold()).strip("_")
    if not nome or not nome[0].isalpha():
        nome = f"x_{nome}"
    if nome in _RISERVATE:
        nome = f"{nome}_"
    if not _SINTATTICO.match(nome):
        raise ValueError(f"nome non sintattico per la colonna {originale!r}: {nome!r}")
    return nome


def nomi_in_collisione(originali: list[str]) -> list[tuple[str, str]]:
    """Le colonne richieste per l'oggetto il cui nome sintattico e' gia' preso.

    Per ciascuna, la coppia (nome originale, nome nell'oggetto): gia' usato da
    una colonna dell'inventario o da una colonna che la precede in
    ``originali``. Due colonne diverse nelle tabelle possono dare lo stesso
    nome (``Var Uno`` e ``Var-Uno``), e la stessa colonna puo' essere chiesta
    a due tabelle.
    """
    usati = set(COLONNE_INVENTARIO)
    collisioni = []
    for originale in originali:
        nome = nome_nell_oggetto(originale)
        if nome in usati:
            collisioni.append((originale, nome))
        usati.add(nome)
    return collisioni


def tabella_delle_cellule(config: Config) -> str | None:
    """Da quale tabella si porta nell'oggetto la colonna dei livelli dei
    controlli positivi (``katharoseq.cell_count_column``) quando non e' gia'
    fra le colonne richieste: ``"lotto"``, ``"studio"`` o ``None``.

    ``None`` se la colonna non e' dichiarata, se ``out.batch_columns`` o
    ``out.study_columns`` la portano gia', o se nessuna delle due tabelle la
    contiene (lo segnala G02). Si cerca prima nel file di arricchimento e poi
    nella tabella di studio: la calibrazione non deve dipendere dal fatto che
    la colonna sia stata anche elencata fra quelle da portare nell'oggetto.
    """
    colonna = config.katharoseq.cell_count_column
    if colonna is None or colonna in (*config.out.batch_columns, *config.out.study_columns):
        return None
    if config.io.batch_table is not None and colonna in intestazione(Path(config.io.batch_table)):
        return "lotto"
    if colonna in intestazione(tabella_di_studio(config)[0]):
        return "studio"
    return None
