"""Coerenza fra cio' che i metadati dichiarano di una riga e il file collegato.

Un deposito puo' descrivere nella stessa tabella campioni di ampliconi diversi
(16S e ITS dello stesso studio), e collegare una riga al file sbagliato: il
campione sarebbe analizzato con le letture di un altro, senza che alcun
controllo sui nomi se ne accorga. Dove i metadati dichiarano per riga
l'amplicone o il primer, qui lo si confronta con il primer che le letture del
file collegato portano davvero in testa.

Il confronto si fa solo dove e' senza ambiguita':

* la dichiarazione deve essere interpretabile: un nome di amplicone noto
  (:data:`AMPLICONI`), una sequenza nucleotidica in codici IUPAC, o il nome di
  un primer del vocabolario;
* il file deve portare in testa un primer riconosciuto, fra quelli che la
  tabella dichiara e quelli del vocabolario: letture a cui il primer e' gia'
  stato tolto non dicono nulla, e la riga non si giudica;
* la riga e' incoerente solo se il primer trovato e' di un altro amplicone, o
  e' quello che la tabella dichiara per altre righe e non per questa.

Uno scambio fra due campioni dello stesso amplicone non lascia traccia nelle
letture: non si riconosce dal contenuto, e resta un limite dichiarato.

Il modulo legge e confronta, non decide: ``init`` ne riporta l'esito nel
riepilogo, la validazione ne fa un arresto (``E-S0-30``).
"""

from __future__ import annotations

import re
import zlib
from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Final

# Importato per nome di modulo: cosi' il vocabolario dei primer, che decide
# l'esito del confronto, entra nel sorgente della validazione.
from amplicon16s.init.vocabolario import ESTENSIONI_FASTQ, PRIMER_16S, PRIMER_16S_REVERSE
from amplicon16s.io_layer.reads import apri_fastq, espandi_iupac
from amplicon16s.metadata.crosswalk import elementi_di_una_cella

__all__ = [
    "AMPLICONI",
    "COLONNA_AMPLICONE",
    "COLONNA_PRIMER",
    "LETTURE_PER_IL_PRIMER",
    "QUOTA_DEL_PRIMER_IN_TESTA",
    "Dichiarazione",
    "Incoerenza",
    "RigaCollegata",
    "amplicone_dichiarato",
    "collega",
    "confronta",
    "dichiarazione",
    "file_della_cella",
    "primer_in_testa",
]

#: Le colonne con cui una tabella di assay ISA-Tab dichiara per riga l'amplicone
#: e i primer: i nomi sono quelli dello standard dei depositi, non di un dataset.
COLONNA_AMPLICONE: Final = "Parameter Value[Target Molecule]"
COLONNA_PRIMER: Final = "Parameter Value[Primer Info]"

#: I nomi con cui i metadati chiamano i marcatori del sequenziamento di
#: ampliconi, con il nome a cui si riportano.
AMPLICONI: Final[Mapping[str, str]] = {
    "16s": "16S", "18s": "18S", "23s": "23S", "28s": "28S", "12s": "12S",
    "its": "ITS", "its1": "ITS", "its2": "ITS", "coi": "COI", "co1": "COI",
}

#: Il solo amplicone di cui il vocabolario conosce i primer.
_AMPLICONE_DEL_VOCABOLARIO: Final = "16S"

#: Quante letture in testa a ogni file si guardano per riconoscere il primer.
#: E' la dimensione del campione di un riconoscimento, non un parametro
#: dell'analisi: un primer in testa lo e' in gran parte delle letture.
LETTURE_PER_IL_PRIMER: Final = 2000

#: La quota di letture che devono cominciare con un primer perche' il file lo
#: porti in testa: la stessa con cui la validazione distingue il primer in testa
#: da uno preceduto da un numero variabile di basi.
QUOTA_DEL_PRIMER_IN_TESTA: Final = 0.5

_SEQUENZA: Final = re.compile(r"(?<![A-Za-z])[ACGTRYSWKMBDHVN]{15,}(?![A-Za-z])")


def amplicone_dichiarato(valore: str) -> str | None:
    """L'amplicone che un valore dichiara, o ``None`` se non lo dice senza
    ambiguita': nessun nome noto, o piu' di uno.
    """
    trovati = {
        AMPLICONI[elemento] for elemento in re.findall(r"[a-z0-9]+", valore.casefold())
        if elemento in AMPLICONI}
    return trovati.pop() if len(trovati) == 1 else None


@dataclass(frozen=True)
class Dichiarazione:
    """Che cosa i metadati dichiarano di una riga."""

    #: L'amplicone, con il nome di :data:`AMPLICONI`; ``None`` se non e' dichiarato.
    amplicone: str | None
    #: I primer, come sequenze in codici IUPAC: quelle scritte nella cella e
    #: quelle dei primer del vocabolario che la cella nomina.
    primer: tuple[str, ...]
    #: Il testo dichiarato, per i messaggi.
    testo: str

    @property
    def interpretabile(self) -> bool:
        """Vero se la riga dichiara qualcosa con cui un primer si puo' confrontare."""
        return self.amplicone is not None or bool(self.primer)


def dichiarazione(
    riga: Mapping[str, str], colonna_amplicone: str | None, colonna_primer: str | None
) -> Dichiarazione:
    """La dichiarazione di una riga, dalle due colonne (una puo' mancare)."""
    amplicone = riga.get(colonna_amplicone, "") if colonna_amplicone else ""
    primer = riga.get(colonna_primer, "") if colonna_primer else ""
    sequenze = list(dict.fromkeys(_SEQUENZA.findall(primer.upper())))
    for noto in (*PRIMER_16S, *PRIMER_16S_REVERSE):
        nome = noto.nome.split()[0]
        if re.search(rf"(?<![A-Za-z0-9]){re.escape(nome)}(?![A-Za-z0-9])", primer,
                     re.IGNORECASE) and noto.sequenza not in sequenze:
            sequenze.append(noto.sequenza)
    # Per i messaggi basta l'amplicone, se la riga lo dichiara: i primer che lo
    # accompagnano allungherebbero ogni voce di un elenco senza aggiungere nulla.
    return Dichiarazione(
        amplicone_dichiarato(amplicone), tuple(sequenze), amplicone.strip() or primer.strip())


def _e_un_file_di_letture(nome: str) -> bool:
    """Vero se il nome ha un'estensione di file FASTQ, senza distinguere le maiuscole."""
    return nome.casefold().endswith(ESTENSIONI_FASTQ)


def file_della_cella(
    valore: str, presenti: Mapping[str, Path]
) -> tuple[str | None, tuple[str, ...]]:
    """Il file di letture a cui una cella collega la riga, e i file che nomina
    senza che siano nella cartella.

    La cella puo' elencare piu' file (i due di una coppia): vale il primo che
    e' fra i ``presenti``. Se la cella nomina solo file di letture e nessuno e'
    presente, il secondo valore li elenca. Una cella che non nomina file (una
    chiave, un codice) da' ``(None, ())``: la riga si collega per chiave.
    """
    elementi = elementi_di_una_cella(valore)
    for elemento in elementi:
        if elemento in presenti:
            return elemento, ()
    if elementi and all(_e_un_file_di_letture(e) for e in elementi):
        return None, tuple(elementi)
    return None, ()


def primer_in_testa(
    percorso: Path, candidati: Sequence[str], letture: int = LETTURE_PER_IL_PRIMER
) -> dict[str, float]:
    """Per ogni primer candidato (in codici IUPAC), la frazione delle prime
    ``letture`` letture del file che cominciano con esso. Vuoto se il file non
    si legge o non ha letture: lo dichiarano i controlli sui file.
    """
    espressioni = [re.compile(espandi_iupac(c)) for c in candidati]
    conteggi = [0] * len(candidati)
    esaminate = 0
    try:
        with apri_fastq(percorso) as file:
            for numero, riga in enumerate(file):
                if numero % 4 != 1:
                    continue
                if esaminate >= letture:
                    break
                esaminate += 1
                sequenza = riga.upper()
                for indice, espressione in enumerate(espressioni):
                    if espressione.match(sequenza):
                        conteggi[indice] += 1
    except (OSError, EOFError, zlib.error, UnicodeDecodeError):
        return {}
    if not esaminate:
        return {}
    return {c: n / esaminate for c, n in zip(candidati, conteggi)}


@dataclass(frozen=True)
class RigaCollegata:
    """Una riga dei metadati con il file a cui e' collegata."""

    nome: str
    file: Path
    dichiarazione: Dichiarazione


def collega(
    righe: Sequence[Mapping[str, str]], colonna_nome: str, colonna_file: str,
    presenti: Mapping[str, Path], colonna_amplicone: str | None, colonna_primer: str | None,
    file_per_valore: Callable[[str], Path | None] | None = None,
) -> tuple[list[RigaCollegata], list[tuple[str, str, Dichiarazione]]]:
    """Le righe con il file di letture a cui sono collegate, e quelle che
    nominano un file assente dalla cartella: (nome del campione, file, cio' che
    la riga dichiara).

    La riga si collega al file che la cella di ``colonna_file`` nomina, se e'
    fra i ``presenti``; se la cella non nomina file, con ``file_per_valore``
    (la chiave del campione estratta dal valore), quando e' data.
    """
    collegate: list[RigaCollegata] = []
    assenti: list[tuple[str, str, Dichiarazione]] = []
    for riga in righe:
        nome = riga.get(colonna_nome, "")
        valore = riga.get(colonna_file, "")
        dichiarata = dichiarazione(riga, colonna_amplicone, colonna_primer)
        file, mancanti = file_della_cella(valore, presenti)
        percorso = presenti[file] if file is not None else None
        if percorso is None and file_per_valore is not None:
            # Anche una cella che nomina un file con un altro nome (un percorso
            # davanti, il nome del deposito) si collega per chiave, se la
            # chiave porta a un file della cartella.
            percorso = file_per_valore(valore)
        if percorso is None and mancanti:
            assenti.append((nome, mancanti[0], dichiarata))
        if percorso is not None:
            collegate.append(RigaCollegata(nome, percorso, dichiarata))
    return collegate, assenti


@dataclass(frozen=True)
class Incoerenza:
    """Una riga collegata a un file le cui letture portano il primer di un altro."""

    nome: str
    file: str
    dichiarato: str
    #: Il primer trovato in testa, con la frazione di letture e che cosa e'.
    trovato: str

    def __str__(self) -> str:
        return (f"{self.nome} (dichiara {self.dichiarato}; file {self.file}: {self.trovato})")


def _nome_del_primer(sequenza: str) -> str:
    """Il nome con cui il vocabolario conosce una sequenza, se la conosce."""
    for noto in (*PRIMER_16S, *PRIMER_16S_REVERSE):
        if noto.sequenza == sequenza:
            return f"{noto.nome.split()[0]} {sequenza}"
    return sequenza


def confronta(
    righe: Sequence[RigaCollegata], letture: int = LETTURE_PER_IL_PRIMER
) -> tuple[list[Incoerenza], int]:
    """Le righe collegate a un file di un altro amplicone, e quante righe si
    sono potute giudicare.

    I primer candidati sono quelli che la tabella dichiara, in qualunque riga,
    e i primer forward del vocabolario; di ciascuno si sa a quali ampliconi
    appartiene dalle righe che lo dichiarano (i primer del vocabolario valgono
    per le righe 16S, ma non provano che una riga di un altro amplicone sia
    sbagliata: alcuni amplificano anche il 18S). Di ogni file si guardano le prime ``letture`` letture: i
    primer con cui ne comincia almeno :data:`QUOTA_DEL_PRIMER_IN_TESTA` sono
    quelli che il file porta. Una riga non si giudica se non dichiara nulla di
    interpretabile o se il suo file non porta alcun primer candidato; e'
    coerente se il file porta un primer che la riga dichiara o un primer del
    suo amplicone; e' incoerente se il file porta soltanto primer di altri
    ampliconi, o dichiarati da altre righe e non da questa.
    """
    del_vocabolario = {p.sequenza for p in PRIMER_16S}
    # A quali ampliconi appartiene un primer lo dicono le sole righe della
    # tabella: i primer del vocabolario sono universali, e alcuni amplificano
    # anche il 18S. Servono a riconoscere il file di una riga 16S, non a
    # smentire una riga di un altro amplicone.
    ampliconi_di: dict[str, set[str]] = {}
    dichiarati_da: dict[str, int] = {}
    for riga in righe:
        for sequenza in riga.dichiarazione.primer:
            dichiarati_da[sequenza] = dichiarati_da.get(sequenza, 0) + 1
            if riga.dichiarazione.amplicone is not None:
                ampliconi_di.setdefault(sequenza, set()).add(riga.dichiarazione.amplicone)
    candidati = sorted({*del_vocabolario, *dichiarati_da})

    in_testa: dict[Path, dict[str, float]] = {}
    incoerenze: list[Incoerenza] = []
    giudicate = 0
    for riga in righe:
        dichiarata = riga.dichiarazione
        if not dichiarata.interpretabile:
            continue
        if riga.file not in in_testa:
            in_testa[riga.file] = {
                c: f for c, f in primer_in_testa(riga.file, candidati, letture).items()
                if f >= QUOTA_DEL_PRIMER_IN_TESTA}
        trovati = in_testa[riga.file]
        if not trovati:
            continue
        propri = set(dichiarata.primer)
        if dichiarata.amplicone == _AMPLICONE_DEL_VOCABOLARIO:
            propri |= del_vocabolario
        dello_stesso_amplicone = dichiarata.amplicone is not None and any(
            dichiarata.amplicone in ampliconi_di.get(c, ()) for c in trovati)
        if propri & set(trovati) or dello_stesso_amplicone:
            giudicate += 1
            continue
        # Nessun primer della riga: e' di un altro solo se si sa di chi e'.
        estranei = [
            c for c in trovati
            if (dichiarata.amplicone is not None and ampliconi_di.get(c)
                and dichiarata.amplicone not in ampliconi_di[c])
            or (dichiarata.primer and c in dichiarati_da)
        ]
        if not estranei:
            continue
        giudicate += 1
        migliore = max(estranei, key=lambda c: (trovati[c], c))
        di_chi = ", ".join(sorted(ampliconi_di.get(migliore, ())))
        incoerenze.append(Incoerenza(
            riga.nome, riga.file.name, dichiarata.testo,
            f"il {trovati[migliore]:.0%} delle letture comincia con il primer "
            f"{_nome_del_primer(migliore)}" + (f", di {di_chi}" if di_chi else "")))
    return incoerenze, giudicate
