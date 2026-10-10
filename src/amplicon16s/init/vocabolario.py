"""Il vocabolario con cui ``init`` legge un dataset.

Qui stanno le sole conoscenze generali che il rilevamento usa: nessuna viene da
un dataset. Ogni elenco dice la sua fonte o la sua ragione.
"""

from __future__ import annotations

from typing import Final, NamedTuple

__all__ = [
    "ESTENSIONI_FASTQ",
    "PAROLE_CELLULE",
    "PAROLE_CORSA",
    "PAROLE_ESTRAZIONE",
    "PAROLE_NEGATIVO",
    "PAROLE_PIASTRA",
    "PAROLE_POSITIVO",
    "PAROLE_TIPO",
    "PAROLA_CONTROLLO",
    "PRIMER_16S",
    "Primer",
    "RIFERIMENTI",
]

#: Le estensioni con cui si riconosce un file di letture, dalla piu' lunga.
ESTENSIONI_FASTQ: Final = (".fastq.gz", ".fq.gz", ".fastq", ".fq")


class Primer(NamedTuple):
    """Un primer forward del gene 16S, con la pubblicazione che lo definisce."""

    nome: str
    sequenza: str
    fonte: str


#: Primer forward 16S di uso comune, in codici IUPAC. L'ordine conta solo a
#: parita' di frazione di letture: prima la forma piu' recente di ogni primer.
PRIMER_16S: Final[tuple[Primer, ...]] = (
    Primer("515F", "GTGYCAGCMGCCGCGGTAA",
           "Parada et al. 2016, Environ Microbiol 18:1403 (Earth Microbiome Project, V4)"),
    Primer("515F (originale)", "GTGCCAGCMGCCGCGGTAA",
           "Caporaso et al. 2011, PNAS 108:4516 (V4)"),
    Primer("341F", "CCTACGGGNGGCWGCAG",
           "Klindworth et al. 2013, Nucleic Acids Res 41:e1 (V3-V4)"),
    Primer("357F", "CCTACGGGAGGCAGCAG",
           "Muyzer et al. 1993, Appl Environ Microbiol 59:695 (V3)"),
    Primer("27F", "AGAGTTTGATCMTGGCTCAG",
           "Lane 1991, in Nucleic Acid Techniques in Bacterial Systematics (V1-V2)"),
    Primer("8F", "AGAGTTTGATCCTGGCTCAG",
           "Edwards et al. 1989, Nucleic Acids Res 17:7843 (V1-V2)"),
    Primer("799F", "AACMGGATTAGATACCCKG",
           "Chelius e Triplett 2001, Microb Ecol 41:252 (V5-V7)"),
)

#: Nomi dei riferimenti tassonomici distribuiti nel formato di dada2, come
#: compaiono nei nomi dei file, con il nome da riportare.
RIFERIMENTI: Final[dict[str, str]] = {
    "silva": "SILVA", "gtdb": "GTDB", "rdp": "RDP", "greengenes2": "Greengenes2",
    "greengenes": "Greengenes", "gg2": "Greengenes2", "gg": "Greengenes",
    "unite": "UNITE", "pr2": "PR2", "homd": "HOMD",
}

#: Parole con cui i metadati chiamano i controlli negativi e i positivi. La
#: sola parola «control» non basta a dire di quale specie sia.
PAROLE_NEGATIVO: Final = ("blank", "negative", "ntc", "no template")
PAROLE_POSITIVO: Final = ("positive", "mock")
PAROLA_CONTROLLO: Final = "control"

#: Parole nel nome di una colonna che dice di che tipo e' il campione.
PAROLE_TIPO: Final = ("type", "tipo", "material")

#: Parole nel nome delle colonne di un file del lotto.
PAROLE_PIASTRA: Final = ("plate", "piastra")
PAROLE_ESTRAZIONE: Final = ("extraction", "estrazione")
PAROLE_CORSA: Final = ("run", "corsa", "lane", "flowcell")
PAROLE_CELLULE: Final = ("cell", "cells", "cellule")
