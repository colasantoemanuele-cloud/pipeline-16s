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
    "PRIMER_16S_REVERSE",
    "Primer",
    "RIFERIMENTI",
]

#: Le estensioni con cui si riconosce un file di letture, dalla piu' lunga.
ESTENSIONI_FASTQ: Final = (".fastq.gz", ".fq.gz", ".fastq", ".fq")


class Primer(NamedTuple):
    """Un primer del gene 16S, con la pubblicazione che lo definisce.

    ``posizione`` e' la posizione, sul 16S di *Escherichia coli*, della base
    all'estremita' 5' del primer (la prima base letta): la numerazione e'
    quella di Brosius et al. 1978, PNAS 75:4801, da cui i primer prendono il
    nome, e l'estensione di ciascuno e' quella della sua pubblicazione. Per un
    primer forward e' la posizione piu' bassa che copre, per un reverse la piu'
    alta: la regione amplificata da una coppia, primer compresi, e' lunga su
    *E. coli* posizione del reverse meno posizione del forward piu' uno.
    """

    nome: str
    sequenza: str
    fonte: str
    posizione: int | None = None


#: Primer forward 16S di uso comune, in codici IUPAC. L'ordine conta solo a
#: parita' di frazione di letture: prima la forma piu' recente di ogni primer.
PRIMER_16S: Final[tuple[Primer, ...]] = (
    Primer("515F", "GTGYCAGCMGCCGCGGTAA",
           "Parada et al. 2016, Environ Microbiol 18:1403 (Earth Microbiome Project, V4)", 515),
    Primer("515F (originale)", "GTGCCAGCMGCCGCGGTAA",
           "Caporaso et al. 2011, PNAS 108:4516 (V4)", 515),
    Primer("341F", "CCTACGGGNGGCWGCAG",
           "Klindworth et al. 2013, Nucleic Acids Res 41:e1 (V3-V4)", 341),
    Primer("357F", "CCTACGGGAGGCAGCAG",
           "Muyzer et al. 1993, Appl Environ Microbiol 59:695 (V3)", 341),
    Primer("27F", "AGAGTTTGATCMTGGCTCAG",
           "Lane 1991, in Nucleic Acid Techniques in Bacterial Systematics (V1-V2)", 8),
    Primer("8F", "AGAGTTTGATCCTGGCTCAG",
           "Edwards et al. 1989, Nucleic Acids Res 17:7843 (V1-V2)", 8),
    Primer("799F", "AACMGGATTAGATACCCKG",
           "Chelius e Triplett 2001, Microb Ecol 41:252 (V5-V7)", 781),
)

#: Primer reverse 16S di uso comune, in codici IUPAC, scritti da 5' a 3' come
#: compaiono in testa alle letture reverse. Stesso criterio d'ordine dei forward.
PRIMER_16S_REVERSE: Final[tuple[Primer, ...]] = (
    Primer("806R", "GGACTACNVGGGTWTCTAAT",
           "Apprill et al. 2015, Aquat Microb Ecol 75:129 (Earth Microbiome Project, V4)", 806),
    Primer("806R (originale)", "GGACTACHVGGGTWTCTAAT",
           "Caporaso et al. 2011, PNAS 108:4516 (V4)", 806),
    Primer("805R", "GACTACHVGGGTATCTAATCC",
           "Herlemann et al. 2011, ISME J 5:1571; lo stesso oligonucleotide e' il 785R di "
           "Klindworth et al. 2013, Nucleic Acids Res 41:e1 (V3-V4)", 805),
    Primer("926R", "CCGYCAATTYMTTTRAGTTT",
           "Parada et al. 2016, Environ Microbiol 18:1403 (V4-V5)", 926),
    Primer("518R", "ATTACCGCGGCTGCTGG",
           "Muyzer et al. 1993, Appl Environ Microbiol 59:695 (V3)", 534),
    Primer("338R", "GCTGCCTCCCGTAGGAGT",
           "Amann et al. 1990, Appl Environ Microbiol 56:1919 (V1-V2)", 355),
    Primer("1492R", "TACGGYTACCTTGTTACGACTT",
           "Lane 1991, in Nucleic Acid Techniques in Bacterial Systematics (V7-V9)", 1513),
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
