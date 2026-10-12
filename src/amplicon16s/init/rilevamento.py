"""Che cosa ``init`` ricava dai file: fatti verificati ed euristiche.

Le funzioni di questo modulo leggono e non decidono. Un **fatto** e' qualcosa
che si verifica su tutti i file con le stesse funzioni che la pipeline usera'
(l'estrazione della chiave, il riconoscimento delle coppie, la regola del
troncamento): se la verifica riesce il valore e' rilevato. Un'**euristica**
guarda nomi e parole, e puo' sbagliare: il suo esito e' sempre una proposta, e
quando piu' risposte sono plausibili le riporta tutte.

Le regole sono riusate, non riscritte: la chiave del campione si estrae con
:func:`~amplicon16s.metadata.crosswalk.estrai_accession` e
:func:`~amplicon16s.metadata.crosswalk.chiave_dalla_tabella`, le coppie si
riconoscono con la scansione di S0 e S1 e con la regola sui nomi di G07, il
profilo e' quello dello script di S1, il troncamento quello del report.
"""

from __future__ import annotations

import csv
import hashlib
import re
import sys
import zlib
from collections import Counter, defaultdict
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass, field, replace
from pathlib import Path
from typing import Any

from amplicon16s.gates.g01_g15 import PRIMER_IN_TESTA, _e_lettura_inversa
from amplicon16s.init import vocabolario as v
from amplicon16s.io_layer.reads import (
    FRAZIONE_DI_ERRORI,
    QUOTA_FUORI_POSIZIONE,
    SFASAMENTO_CERCATO,
    RicercaPrimer,
    discordanze_di_partenza,
    MARCATORI_DI_FILE,
    StatisticheFile,
    abbina_file,
    apri_fastq,
    espandi_iupac,
    scansiona_file,
)
from amplicon16s.metadata.crosswalk import (
    AccessionNonEstraibile,
    chiave_dalla_tabella,
    estrai_accession,
)
from amplicon16s.metadata.tabelle import intestazione, leggi_tsv

__all__ = [
    "Chiave",
    "Ispezione",
    "PrimerOsservato",
    "PrimerSfasato",
    "LETTURE_PER_IL_MOTIVO",
    "SOGLIA_BASE_DEL_MOTIVO",
    "colonna_del_campione",
    "colonne_del_lotto",
    "colonne_del_tipo",
    "chiave_del_campione",
    "chiave_del_lotto",
    "chiavi_dei_file",
    "estensione_di",
    "file_di_letture",
    "ispeziona",
    "md5_di",
    "motivo_dal_consenso",
    "nome_e_versione_del_riferimento",
    "primer_dai_metadati",
    "primer_nominati",
    "primer_nelle_letture",
    "secondi_di_coppia",
]

#: Quante letture per file entrano nel consenso da cui si propone il motivo
#: conservato. E' la dimensione del campione di un'euristica, non un parametro
#: dell'analisi: la proposta si verifica poi con ``validate``.
LETTURE_PER_IL_MOTIVO = 2000

#: La quota di letture che una base (o una coppia di basi) deve coprire perche'
#: una posizione del motivo proposto sia quella base (o quella classe).
SOGLIA_BASE_DEL_MOTIVO = 0.9

#: Quante basi in testa a ogni lettura si contano per il consenso.
_BASI_IN_TESTA = 96


# --------------------------------------------------------------------------- #
# File di letture                                                              #
# --------------------------------------------------------------------------- #


def estensione_di(nome: str) -> str | None:
    """L'estensione FASTQ del nome, o ``None`` se non e' un file di letture."""
    minuscolo = nome.casefold()
    for estensione in v.ESTENSIONI_FASTQ:
        if minuscolo.endswith(estensione):
            return nome[len(nome) - len(estensione):]
    return None


def file_di_letture(cartella: Path) -> dict[str, list[Path]]:
    """I file di letture della cartella, per estensione, senza sottocartelle."""
    per_estensione: dict[str, list[Path]] = defaultdict(list)
    for percorso in sorted(cartella.iterdir()):
        estensione = estensione_di(percorso.name) if percorso.is_file() else None
        if estensione is not None:
            per_estensione[estensione].append(percorso)
    return dict(per_estensione)


def coppie_di_file(nomi: Sequence[str]) -> tuple[tuple[str, str], dict[str, str]] | None:
    """I marcatori con cui ogni file della cartella ha il suo compagno, e le
    coppie (file forward -> file reverse); ``None`` se con nessuna coppia di
    marcatori i file si abbinano tutti.

    Si provano le convenzioni di :data:`MARCATORI_DI_FILE`, dalla piu'
    specifica. Basta un file senza compagno perche' la cartella non sia quella
    di un dataset appaiato: sara' chi lancia il comando a dire che cosa sono i
    file rimasti.
    """
    for marcatori in MARCATORI_DI_FILE:
        abbinamento = abbina_file(nomi, marcatori)
        if abbinamento.completo:
            return marcatori, dict(abbinamento.coppie)
    return None


# --------------------------------------------------------------------------- #
# Lunghezza dell'inserto, dalle coppie                                         #
# --------------------------------------------------------------------------- #

#: Quante coppie per file entrano nella stima della lunghezza dell'inserto. E'
#: la dimensione del campione di un'euristica: la stima si verifica poi sulla
#: quota di coppie fuse della prima esecuzione.
COPPIE_PER_L_INSERTO = 500

#: I semi con cui si cerca dove le due letture si sovrappongono: tratti di
#: questa lunghezza della lettura forward, uno ogni ``_PASSO_DEI_SEMI`` basi.
_SEME = 12
_PASSO_DEI_SEMI = 4

#: Una sovrapposizione conta se e' lunga almeno queste basi e se non piu' di
#: questa quota discorda: le code delle letture hanno errori, e un
#: allineamento esatto ne perderebbe la maggior parte.
SOVRAPPOSIZIONE_MINIMA_STIMATA = 20
DISCORDANZE_NELLA_SOVRAPPOSIZIONE = 0.1

_COMPLEMENTO = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def complemento_inverso(sequenza: str) -> str:
    """Il complemento inverso di una sequenza di basi."""
    return sequenza.translate(_COMPLEMENTO)[::-1]


def lunghezza_dell_inserto(forward: str, reverse: str) -> int | None:
    """La lunghezza del frammento di cui ``forward`` e ``reverse`` sono le due
    estremita', o ``None`` se le due letture non si sovrappongono.

    La lettura reverse si porta sullo stesso filamento della forward
    (complemento inverso): se le due letture coprono lo stesso frammento, un
    tratto della forward compare nella reverse, a una distanza che dice la
    lunghezza del frammento. I semi propongono la distanza, la sovrapposizione
    intera la conferma. Comprende i primer, se le letture li portano; vale
    anche quando il frammento e' piu' corto delle letture.
    """
    altra = complemento_inverso(reverse)
    voti: Counter[int] = Counter()
    for inizio in range(0, len(forward) - _SEME + 1, _PASSO_DEI_SEMI):
        seme = forward[inizio:inizio + _SEME]
        posto = altra.find(seme)
        while posto >= 0:
            voti[inizio - posto] += 1
            posto = altra.find(seme, posto + 1)
    for scarto, _ in voti.most_common(3):
        da, a = max(0, scarto), min(len(forward), scarto + len(altra))
        if a - da < SOVRAPPOSIZIONE_MINIMA_STIMATA:
            continue
        diverse = sum(x != y for x, y in zip(forward[da:a], altra[da - scarto:a - scarto]))
        if diverse <= DISCORDANZE_NELLA_SOVRAPPOSIZIONE * (a - da):
            return scarto + len(altra)
    return None


def inserti_di_una_coppia_di_file(argomenti: tuple[Any, ...]) -> list[int | None]:
    """La lunghezza dell'inserto di ciascuna delle prime
    :data:`COPPIE_PER_L_INSERTO` coppie dei due file (``None`` dove le letture
    non si sovrappongono). Una sola tupla di argomenti: la funzione gira in un
    processo per campione.

    Con due ricerche del primer in coda alla tupla (quella della forward e
    quella della reverse) la lunghezza e' quella dell'inserto senza i primer e
    senza le basi che li precedono: dal frammento misurato si toglie, per ogni
    coppia, cio' che il taglio toglie alle due letture; una coppia in cui manca
    un primer non si misura.
    """
    forward, reverse, *ricerche = argomenti
    lunghezze: list[int | None] = []
    try:
        with apri_fastq(forward) as primo, apri_fastq(reverse) as secondo:
            for numero, (riga_1, riga_2) in enumerate(zip(primo, secondo)):
                if numero % 4 != 1:
                    continue
                if len(lunghezze) >= COPPIE_PER_L_INSERTO:
                    break
                prima, seconda = riga_1.strip().upper(), riga_2.strip().upper()
                misurata = lunghezza_dell_inserto(prima, seconda)
                if ricerche and misurata is not None:
                    tagli = (ricerche[0].fine(prima), ricerche[1].fine(seconda))
                    misurata = (
                        None if None in tagli or misurata <= sum(tagli)  # type: ignore[arg-type]
                        else misurata - sum(tagli))  # type: ignore[arg-type]
                lunghezze.append(misurata)
    except (OSError, EOFError, zlib.error):
        pass
    return lunghezze


@dataclass(frozen=True)
class Inserto:
    """La lunghezza dell'inserto stimata dalle coppie campionate."""

    #: Coppie campionate e coppie le cui letture si sovrappongono.
    coppie: int
    sovrapposte: int
    #: Mediana e percentili 0,5 e 99,5 delle lunghezze misurate.
    mediana: int
    minima: int
    massima: int


def stima_dell_inserto(per_file: Iterable[Sequence[int | None]]) -> Inserto | None:
    """La distribuzione della lunghezza dell'inserto, o ``None`` se non si
    stima: nessuna coppia campionata, o meno di meta' con le letture
    sovrapposte (le letture non coprono il frammento, e la lunghezza non si
    ricava da loro).

    I percentili 0,5 e 99,5, e non il minimo e il massimo: fra migliaia di
    coppie ce n'e' sempre qualcuna fuori bersaglio o chimerica.
    """
    tutte = [x for lunghezze in per_file for x in lunghezze]
    misurate = sorted(x for x in tutte if x is not None)
    if not tutte or 2 * len(misurate) < len(tutte):
        return None

    def percentile(quota: float) -> int:
        return misurate[min(len(misurate) - 1, max(0, round(quota * (len(misurate) - 1))))]

    return Inserto(len(tutte), len(misurate), percentile(0.5), percentile(0.005),
                   percentile(0.995))


def secondi_di_coppia(nomi: Sequence[str]) -> list[str]:
    """I file che il nome indica come seconda lettura di una coppia, con la
    regola di G07.
    """
    presenti = frozenset(n.casefold() for n in nomi)
    return sorted(n for n in nomi if _e_lettura_inversa(n, presenti))


@dataclass(frozen=True)
class Ispezione:
    """Cio' che si e' misurato sulle letture di un file."""

    statistiche: StatisticheFile
    #: Letture che iniziano con ciascun primer dell'elenco, nell'ordine dell'elenco.
    con_primer: tuple[int, ...]
    #: Per posizione, quante delle prime letture portano ciascuna base.
    basi: tuple[dict[str, int], ...]
    #: Quante letture sono entrate in ``basi``.
    letture_nel_consenso: int
    #: Per ciascun primer cercato lungo la lettura, fra le letture entrate nel
    #: consenso: a quale sfasamento si trova, e in quante letture.
    sfasamenti: tuple[dict[int, int], ...] = ()
    #: Per ciascuno di quei primer e per posizione dopo la sua fine, quante
    #: letture portano ciascuna base.
    basi_dopo: tuple[tuple[dict[str, int], ...], ...] = ()

    def dopo_il_primer(self, indice: int) -> Ispezione:
        """La stessa ispezione con, al posto delle basi in testa, quelle che
        seguono il primer trovato: per il consenso del motivo.
        """
        return replace(self, basi=self.basi_dopo[indice])


#: La quota di letture con il primer sotto la quale uno sfasamento non conta
#: per lo sfasamento massimo proposto: qualche lettura isolata non allarga la
#: finestra.
QUOTA_DI_UNO_SFASAMENTO = 0.01


def ispeziona(argomenti: tuple[Any, ...]) -> Ispezione:
    """Legge un file di letture: i segni di coppia e la struttura con la
    scansione di S0 e S1, e le letture che iniziano con ciascun primer.

    ``argomenti`` e' la terna (percorso, letture da ispezionare o ``None`` per
    tutte, espressioni dei primer): una sola tupla, perche' la funzione gira in
    un processo per file. Un quarto elemento, le sequenze IUPAC degli stessi
    primer, li fa cercare anche lungo la lettura
    (:class:`~amplicon16s.io_layer.reads.RicercaPrimer`, fino a
    :data:`SFASAMENTO_CERCATO` basi di sfasamento) nelle letture del consenso.
    """
    percorso, quante, primer, *resto = argomenti
    limite = sys.maxsize if quante is None else quante
    statistiche = scansiona_file(percorso, limite, None, None)
    espressioni = [re.compile(p) for p in primer]
    conteggi = [0] * len(espressioni)
    basi: list[Counter[str]] = [Counter() for _ in range(_BASI_IN_TESTA)]
    nel_consenso = 0
    ricerche = [
        RicercaPrimer(sequenza, discordanze_di_partenza(sequenza), SFASAMENTO_CERCATO)
        for sequenza in (resto[0] if resto else ())]
    sfasamenti: list[Counter[int]] = [Counter() for _ in ricerche]
    dopo: list[list[Counter[str]]] = [
        [Counter() for _ in range(_BASI_IN_TESTA)] for _ in ricerche]
    if statistiche.valido:
        with apri_fastq(percorso) as file:
            for numero, riga in enumerate(file):
                if numero % 4 != 1:
                    continue
                lettura = numero // 4
                if lettura >= limite:
                    break
                for indice, espressione in enumerate(espressioni):
                    if espressione.match(riga):
                        conteggi[indice] += 1
                if lettura < LETTURE_PER_IL_MOTIVO:
                    nel_consenso += 1
                    for posizione, base in enumerate(riga[:_BASI_IN_TESTA].rstrip("\r\n")):
                        basi[posizione][base] += 1
                    sequenza = riga.rstrip("\r\n")
                    for indice, ricerca in enumerate(ricerche):
                        trovata = ricerca.trova(sequenza)
                        if trovata is None:
                            continue
                        sfasamenti[indice][trovata[0]] += 1
                        fine = trovata[0] + len(ricerca.primer)
                        for posizione, base in enumerate(
                                sequenza[fine:fine + _BASI_IN_TESTA]):
                            dopo[indice][posizione][base] += 1
    return Ispezione(
        statistiche, tuple(conteggi), tuple(dict(b) for b in basi), nel_consenso,
        tuple(dict(s) for s in sfasamenti),
        tuple(tuple(dict(b) for b in del_primer) for del_primer in dopo))


@dataclass(frozen=True)
class PrimerSfasato:
    """Un primer dell'elenco trovato lungo le letture, con gli sfasamenti."""

    primer: v.Primer
    #: Il posto del primer nell'elenco: l'indice di ``Ispezione.basi_dopo``.
    indice: int
    #: Mediana, minimo e massimo fra i file della frazione di letture con il
    #: primer entro :data:`SFASAMENTO_CERCATO` basi.
    mediana: float
    minima: float
    massima: float
    #: Sull'insieme dei file: sfasamento -> letture.
    sfasamenti: dict[int, int]

    @property
    def fuori_posizione(self) -> float:
        """La quota di letture con il primer a uno sfasamento diverso dal piu' frequente."""
        totale = sum(self.sfasamenti.values())
        return 1 - max(self.sfasamenti.values()) / totale if totale else 0.0

    @property
    def sfasamento_massimo(self) -> int:
        """Il piu' grande sfasamento che porta almeno
        :data:`QUOTA_DI_UNO_SFASAMENTO` delle letture con il primer.
        """
        totale = sum(self.sfasamenti.values())
        return max(
            (s for s, n in self.sfasamenti.items() if n >= QUOTA_DI_UNO_SFASAMENTO * totale),
            default=0)

    @property
    def da_cercare(self) -> bool:
        """Vero se il primer non sta a una posizione fissa: va cercato."""
        return self.fuori_posizione > QUOTA_FUORI_POSIZIONE


def primer_sfasato(
    ispezioni: Mapping[str, Ispezione], elenco: Sequence[v.Primer],
) -> PrimerSfasato | None:
    """Il primer dell'elenco che le letture portano entro
    :data:`SFASAMENTO_CERCATO` basi dall'inizio, se ce n'e' uno.

    La regola e' quella di :func:`primer_nelle_letture`, con il primer cercato
    lungo la lettura invece che in testa: la mediana fra i file della frazione
    di letture che lo portano deve raggiungere ``PRIMER_IN_TESTA``; fra piu'
    primer vale la mediana piu' alta, poi la forma meno degenerata, poi
    l'ordine dell'elenco. Chi chiama guarda :attr:`PrimerSfasato.da_cercare`:
    un primer che sta sempre allo stesso posto non richiede la ricerca.
    """
    import statistics

    valide = {n: i for n, i in ispezioni.items()
              if i.statistiche.valido and i.letture_nel_consenso and i.sfasamenti}
    if not valide:
        return None
    frazioni = [
        {n: sum(i.sfasamenti[k].values()) / i.letture_nel_consenso for n, i in valide.items()}
        for k in range(len(elenco))
    ]
    mediane = [statistics.median(f.values()) for f in frazioni]
    migliore = max(
        range(len(mediane)),
        key=lambda k: (mediane[k], -_degenerazione(elenco[k].sequenza), -k))
    if mediane[migliore] < PRIMER_IN_TESTA:
        return None
    totali: Counter[int] = Counter()
    for ispezione in valide.values():
        totali.update(ispezione.sfasamenti[migliore])
    return PrimerSfasato(
        elenco[migliore], migliore, mediane[migliore], min(frazioni[migliore].values()),
        max(frazioni[migliore].values()), dict(sorted(totali.items())))


@dataclass(frozen=True)
class PrimerOsservato:
    """Che cosa si e' visto di un primer dell'elenco in testa alle letture."""

    primer: v.Primer
    #: Mediana, minimo e massimo fra i file della frazione di letture che
    #: iniziano con il primer.
    mediana: float
    minima: float
    massima: float
    #: I file con letture in cui la frazione non raggiunge la soglia.
    sotto_soglia: tuple[str, ...]

    @property
    def in_ogni_file(self) -> bool:
        """Vero se ogni file con letture porta il primer in testa."""
        return not self.sotto_soglia


def _degenerazione(sequenza: str) -> int:
    """Quante posizioni della sequenza non sono una base sola."""
    return sum(1 for base in sequenza if base not in "ACGT")


def primer_nelle_letture(
    ispezioni: dict[str, Ispezione], elenco: Sequence[v.Primer] = v.PRIMER_16S,
) -> tuple[PrimerOsservato | None, list[tuple[v.Primer, float]]]:
    """Il primer dell'elenco che sta in testa alle letture, se ce n'e' uno.

    Per ogni primer si misura, in ogni file con letture, la frazione di letture
    che iniziano con esso. Il primer e' osservato se la mediana fra i file
    raggiunge la soglia con cui G10 distingue il primer in testa
    (``PRIMER_IN_TESTA``); fra piu' primer vale la mediana piu' alta, poi la
    forma meno degenerata (una forma degenerata corrisponde anche alle letture
    della forma originale), poi l'ordine dell'elenco. Chi chiama decide che
    cosa farne: e' un fatto solo se vale per ogni file
    (:attr:`PrimerOsservato.in_ogni_file`). Il secondo valore sono i primer
    visti sotto la soglia, con la frazione sul totale delle letture.
    """
    import statistics

    valide = {n: i for n, i in ispezioni.items()
              if i.statistiche.valido and i.statistiche.letture_esaminate}
    if not valide:
        return None, []
    frazioni = [
        {n: i.con_primer[k] / i.statistiche.letture_esaminate for n, i in valide.items()}
        for k in range(len(elenco))
    ]
    mediane = [statistics.median(f.values()) for f in frazioni]
    migliore = max(
        range(len(mediane)),
        key=lambda k: (mediane[k], -_degenerazione(elenco[k].sequenza), -k))
    if mediane[migliore] >= PRIMER_IN_TESTA:
        del_primer = frazioni[migliore]
        return PrimerOsservato(
            elenco[migliore], mediane[migliore], min(del_primer.values()),
            max(del_primer.values()),
            tuple(sorted(n for n, f in del_primer.items() if f < PRIMER_IN_TESTA)),
        ), []
    # Sotto la soglia: quanto se ne vede comunque, su tutte le letture.
    totale = sum(i.statistiche.letture_esaminate for i in valide.values())
    parziali = [
        (elenco[k], sum(i.con_primer[k] for i in valide.values()) / totale)
        for k in range(len(elenco))
    ]
    return None, [(p, f) for p, f in parziali if f > 0.01]


#: La quota di posizioni in cui la base piu' frequente di un file deve
#: differire da quella della maggioranza dei file perche' il file dissenta: due
#: regioni diverse dello stesso gene concordano per caso in circa un quarto
#: delle posizioni, un file della stessa regione nella maggior parte.
QUOTA_DI_DISSENSO = 0.5


def file_che_dissentono(
    ispezioni: Mapping[str, Ispezione], inizio: int, lunghezza: int
) -> list[str]:
    """I file le cui letture non iniziano come quelle della maggioranza dei file.

    Per ogni file la base piu' frequente a ciascuna delle ``lunghezza``
    posizioni da ``inizio``; il riferimento e', posizione per posizione, la
    base scelta da piu' file (un voto per file, non per lettura: un file
    profondo non decide per gli altri). Dissente il file che ne differisce in
    piu' di :data:`QUOTA_DI_DISSENSO` delle posizioni. Un file di letture
    inverse fra file di letture forward dissente: lasciato nel consenso
    allargherebbe il motivo fino a comprenderlo. Se dissente la meta' dei file
    o piu' non c'e' una maggioranza, e l'elenco e' vuoto.
    """
    prime: dict[str, list[str]] = {}
    for nome, ispezione in ispezioni.items():
        if ispezione.letture_nel_consenso and inizio + lunghezza <= len(ispezione.basi):
            prime[nome] = [
                max(sorted(conteggi), key=conteggi.get) if conteggi else "."
                for conteggi in ispezione.basi[inizio:inizio + lunghezza]
            ]
    if len(prime) < 3:
        return []
    riferimento = [
        Counter(basi[posizione] for basi in prime.values()).most_common(1)[0][0]
        for posizione in range(lunghezza)
    ]
    dissentono = sorted(
        nome for nome, basi in prime.items()
        if sum(a != b for a, b in zip(basi, riferimento)) > QUOTA_DI_DISSENSO * lunghezza
    )
    return dissentono if 2 * len(dissentono) < len(prime) else []


def motivo_dal_consenso(
    ispezioni: Iterable[Ispezione], inizio: int, lunghezza: int
) -> str | None:
    """Il motivo conservato proposto: il consenso delle basi da ``inizio``.

    Una posizione e' una base se questa copre almeno ``SOGLIA_BASE_DEL_MOTIVO``
    delle letture contate, una classe di due basi se insieme la coprono,
    altrimenti ``.``. I ``.`` finali si tolgono. ``None`` se non resta nulla o
    se le letture non arrivano fin li'.
    """
    somme: list[Counter[str]] = [Counter() for _ in range(lunghezza)]
    for ispezione in ispezioni:
        for posizione in range(lunghezza):
            if inizio + posizione < len(ispezione.basi):
                somme[posizione].update(ispezione.basi[inizio + posizione])
    simboli = []
    for conteggi in somme:
        totale = sum(conteggi.values())
        if not totale:
            return None
        frequenti = [(b, n) for b, n in conteggi.most_common() if b in "ACGT"]
        if frequenti and frequenti[0][1] / totale >= SOGLIA_BASE_DEL_MOTIVO:
            simboli.append(frequenti[0][0])
        elif len(frequenti) > 1 and (frequenti[0][1] + frequenti[1][1]) / totale >= SOGLIA_BASE_DEL_MOTIVO:
            simboli.append("[" + "".join(sorted(b for b, _ in frequenti[:2])) + "]")
        else:
            simboli.append(".")
    while simboli and simboli[-1] == ".":
        simboli.pop()
    # Un motivo quasi tutto libero non verifica nulla.
    if sum(1 for s in simboli if s != ".") < max(4, lunghezza // 2):
        return None
    return "".join(simboli)


# --------------------------------------------------------------------------- #
# Chiave del campione                                                          #
# --------------------------------------------------------------------------- #

_SEPARATORI = "_.-"


def _senza_estensione(nome: str) -> str:
    estensione = estensione_di(nome)
    return nome[: -len(estensione)] if estensione else nome


def _forma(chiavi: Sequence[str]) -> str | None:
    """Un'espressione che descrive la forma comune delle chiavi: i tratti di
    lettere uguali in tutte restano alla lettera, le cifre diventano
    ``[0-9]+``. ``None`` se le chiavi non hanno la stessa struttura.
    """
    def tratti(chiave: str) -> list[tuple[str, str]]:
        return [
            ("L" if t[0].isalpha() else "C" if t[0].isdigit() else "A", t)
            for t in re.findall(r"[A-Za-z]+|[0-9]+|[^A-Za-z0-9]+", chiave)
        ]

    spezzate = [tratti(c) for c in chiavi]
    struttura = [tipo for tipo, _ in spezzate[0]]
    if not struttura or any([tipo for tipo, _ in s] != struttura for s in spezzate):
        return None
    parti = []
    for posizione, tipo in enumerate(struttura):
        valori = {s[posizione][1] for s in spezzate}
        if len(valori) == 1:
            parti.append(re.escape(next(iter(valori))))
        elif tipo == "C":
            parti.append("[0-9]+")
        elif tipo == "L":
            parti.append(
                "[A-Z]+" if all(x.isupper() for x in valori)
                else "[a-z]+" if all(x.islower() for x in valori) else "[A-Za-z]+"
            )
        else:
            return None
    return "".join(parti)


def _espressioni_candidate(nomi: Sequence[str], chiavi: Sequence[str]) -> list[str]:
    """Le espressioni che potrebbero estrarre ``chiavi[i]`` da ``nomi[i]``:
    la forma delle chiavi da sola, poi ancorata a cio' che la precede o la
    segue nel nome quando e' costante, poi la parte fino al primo separatore.
    """
    forme = [f for f in (_forma(chiavi), "[A-Za-z0-9]+") if f]
    prima = {n[: n.find(c)] if c in n else None for n, c in zip(nomi, chiavi)}
    dopo = {n[n.find(c) + len(c):] if c in n else None for n, c in zip(nomi, chiavi)}
    candidate = [forme[0]] if _forma(chiavi) else []
    for forma in forme:
        if len(prima) == 1 and None not in prima:
            candidate.append(f"^{re.escape(next(iter(prima)))}({forma})")
        if len(dopo) == 1 and None not in dopo:
            candidate.append(f"({forma}){re.escape(next(iter(dopo)))}$")
    if len(prima) == 1 and len(dopo) == 1 and None not in prima | dopo:
        candidate.append(f"^{re.escape(next(iter(prima)))}(.+){re.escape(next(iter(dopo)))}$")
    if prima == {""}:
        for separatore in _SEPARATORI:
            if all(separatore not in c for c in chiavi):
                candidate.append(f"^([^{re.escape(separatore)}]+)")
    return list(dict.fromkeys(candidate))


def _verificata(espressione: str, nomi: Sequence[str], valori: Sequence[str]) -> list[str] | None:
    """Le chiavi che l'espressione estrae dai nomi dei file, se sono tutte
    distinte e sono le stesse che la pipeline ricava dai valori della colonna,
    una riga per file. Altrimenti ``None``.
    """
    try:
        compilata = re.compile(espressione)
    except re.error:
        return None
    try:
        dai_file = [estrai_accession(n, compilata) for n in nomi]
        dalla_tabella = [chiave_dalla_tabella(x, compilata) for x in valori]
    except AccessionNonEstraibile:
        return None
    if len(set(dai_file)) != len(dai_file) or len(set(dalla_tabella)) != len(dalla_tabella):
        return None
    return dai_file if set(dai_file) == set(dalla_tabella) else None


def _chiavi_possibili(nomi: Sequence[str], valori: Sequence[str]) -> list[list[tuple[str, str]]]:
    """I modi di mettere in corrispondenza file e valori di una colonna: per
    ciascuno, le coppie (nome del file, chiave).

    Due modi: il valore intero (senza estensione) compare nel nome del file,
    delimitato; oppure valore e nome condividono un elemento alfanumerico che
    nessun altro file ha.
    """
    radici = {n: _senza_estensione(n) for n in nomi}
    modi: list[list[tuple[str, str]]] = []

    # 1. Il valore intero nel nome del file, delimitato da cio' che non e' una
    #    lettera o una cifra: "S1" non corrisponde a "S10".
    coppie = []
    for valore in valori:
        chiave = _senza_estensione(valore)
        delimitata = re.compile(rf"(?<![A-Za-z0-9]){re.escape(chiave)}(?![A-Za-z0-9])")
        trovati = [n for n, radice in radici.items() if delimitata.search(radice)]
        if len(trovati) != 1:
            coppie = []
            break
        coppie.append((trovati[0], chiave))
    if coppie and len({n for n, _ in coppie}) == len(nomi):
        modi.append(coppie)

    # 2. Un elemento alfanumerico in comune, proprio di un solo file.
    elementi = {n: re.findall(r"[A-Za-z0-9]+", radice) for n, radice in radici.items()}
    quanti = Counter(e for lista in elementi.values() for e in set(lista))
    di_chi = {e: n for n, lista in elementi.items() for e in lista if quanti[e] == 1}
    condivisi: list[tuple[str, list[str]]] = []
    for valore in valori:
        del_valore = set(re.findall(r"[A-Za-z0-9]+", _senza_estensione(valore)))
        file = {di_chi[e] for e in del_valore if e in di_chi}
        if len(file) != 1:
            condivisi = []
            break
        (nome,) = file
        condivisi.append((nome, [e for e in elementi[nome] if e in del_valore and di_chi.get(e) == nome]))
    if condivisi and len({n for n, _ in condivisi}) == len(nomi):
        quante = {len(e) for _, e in condivisi}
        for indice in range(min(quante)):
            modi.append([(nome, comuni[indice]) for nome, comuni in condivisi])
    return modi


@dataclass(frozen=True)
class Chiave:
    """Una colonna da cui si ricava la chiave del campione, con l'espressione."""

    colonna: str
    espressione: str
    #: Vero se i valori della colonna sono nomi di file di letture.
    di_file: bool
    #: Chiave -> nome del file.
    file_per_chiave: dict[str, str] = field(compare=False, default_factory=dict)
    #: Le altre espressioni verificate per la stessa colonna.
    alternative: tuple[str, ...] = ()

    @property
    def debole(self) -> bool:
        """Vero se la corrispondenza puo' essere una coincidenza: una chiave di
        sole cifre, o piu' corta di tre caratteri, la da' anche una colonna di
        numeri d'ordine o di misure.
        """
        return any(c.isdigit() or len(c) < 3 for c in self.file_per_chiave)


def chiave_del_campione(
    nomi: Sequence[str], righe: Sequence[dict[str, str]], colonne: Sequence[str]
) -> list[Chiave]:
    """Le colonne della tabella di assay che corrispondono uno a uno ai file.

    Per ogni colonna con valori tutti presenti e distinti, tanti quanti i file,
    si cercano le espressioni che estraggono da ogni file una chiave unica e
    dalla colonna la stessa chiave: la verifica e' fatta su tutti i file e su
    tutte le righe, con le funzioni della pipeline.
    """
    candidate = []
    for colonna in dict.fromkeys(colonne):
        valori = [r.get(colonna, "") for r in righe]
        if len(valori) != len(nomi) or not all(valori) or len(set(valori)) != len(valori):
            continue
        verificate: dict[str, list[str]] = {}
        for coppie in _chiavi_possibili(nomi, valori):
            ordinati = [n for n, _ in coppie]
            for espressione in _espressioni_candidate(ordinati, [c for _, c in coppie]):
                if espressione in verificate:
                    continue
                chiavi = _verificata(espressione, ordinati, valori)
                if chiavi is not None:
                    verificate[espressione] = [f"{c}\t{n}" for c, n in zip(chiavi, ordinati)]
        if verificate:
            # Prima le espressioni che descrivono la sola forma della chiave,
            # senza ancorarla a cio' che la circonda: sono le piu' leggibili.
            verificate = dict(sorted(verificate.items(), key=lambda voce: "(" in voce[0]))
            prima = next(iter(verificate))
            di_file = all(
                estensione_di(x) is not None and not re.search(r"[,;\s]", x) for x in valori
            )
            candidate.append(Chiave(
                colonna, prima, di_file,
                dict(riga.split("\t") for riga in verificate[prima]),
                tuple(list(verificate)[1:]),
            ))
    return candidate


def chiavi_dei_file(espressione: str, nomi: Sequence[str]) -> dict[str, str]:
    """Nome del file -> chiave estratta dall'espressione (gia' verificata)."""
    compilata = re.compile(espressione)
    return {n: estrai_accession(n, compilata) for n in nomi}


def chiave_del_lotto(
    chiavi: Iterable[str], espressione: str, righe: Sequence[dict[str, str]],
    colonne: Sequence[str],
) -> list[str]:
    """Le colonne del file del lotto che danno, per ogni campione, la sua
    chiave in una e una sola riga.
    """
    attese = set(chiavi)
    compilata = re.compile(espressione)
    trovate = []
    for colonna in colonne:
        conteggio: Counter[str] = Counter()
        for riga in righe:
            try:
                conteggio[chiave_dalla_tabella(riga.get(colonna, ""), compilata)] += 1
            except AccessionNonEstraibile:
                continue
        if attese and all(conteggio.get(c) == 1 for c in attese):
            trovate.append(colonna)
    return trovate


# --------------------------------------------------------------------------- #
# Tabelle                                                                      #
# --------------------------------------------------------------------------- #


def colonna_del_campione(
    assay: Sequence[dict[str, str]], colonne_assay: Sequence[str],
    studio: Sequence[dict[str, str]] | None, colonne_studio: Sequence[str],
) -> list[tuple[str, str]]:
    """Le coppie (colonna dell'assay, colonna dello studio) che agganciano ogni
    riga dell'assay a una e una sola riga dello studio.

    Prima le colonne con lo stesso nome nelle due tabelle; se nessuna aggancia,
    le coppie di nome diverso. Senza tabella di studio l'elenco e' vuoto.
    """
    if studio is None:
        return []

    def aggancia(colonna_assay: str, colonna_studio: str) -> bool:
        nello_studio = Counter(r.get(colonna_studio, "") for r in studio)
        valori = [r.get(colonna_assay, "") for r in assay]
        return bool(valori) and all(x and nello_studio.get(x) == 1 for x in valori)

    omonime = [(c, c) for c in colonne_assay if c in colonne_studio and aggancia(c, c)]
    if omonime:
        return omonime
    return [(a, s) for a in colonne_assay for s in colonne_studio if a != s and aggancia(a, s)]


def _contiene(testo: str, parole: Iterable[str]) -> bool:
    """Vero se il testo contiene una delle parole, come parola intera."""
    minuscolo = testo.casefold()
    return any(re.search(rf"(?<![a-z]){re.escape(p)}(?![a-z])", minuscolo) for p in parole)


@dataclass(frozen=True)
class TipoDiCampione:
    """Una colonna candidata a dire il tipo dei campioni, con i suoi valori."""

    colonna: str
    negativi: tuple[str, ...]
    positivi: tuple[str, ...]
    biologici: tuple[str, ...]
    #: Valori con la sola parola «control»: non attribuiti.
    non_attribuiti: tuple[str, ...]
    #: Come la colonna e' stata riconosciuta: "parole" o "nome".
    da: str


def colonne_del_tipo(righe: Sequence[dict[str, str]], colonne: Sequence[str]) -> list[TipoDiCampione]:
    """Le colonne candidate a ``ctrl.column``, la piu' plausibile per prima.

    Candidate le colonne con almeno un valore che contiene una parola dei
    controlli; fra piu' candidate vengono prima quelle che lasciano almeno un
    valore ai biologici, poi quella con meno valori distinti (una
    classificazione e' piu' grossolana di una posizione), poi l'ordine della
    tabella. Senza candidate, le colonne il cui nome contiene una parola
    del tipo, con tutti i valori come biologici.
    """
    per_parole, per_nome = [], []
    for posizione, colonna in enumerate(colonne):
        valori = sorted({r.get(colonna, "") for r in righe} - {""})
        if not valori:
            continue
        negativi = [x for x in valori if _contiene(x, v.PAROLE_NEGATIVO)]
        positivi = [x for x in valori if x not in negativi and _contiene(x, v.PAROLE_POSITIVO)]
        dubbi = [x for x in valori if x not in negativi + positivi
                 and _contiene(x, (v.PAROLA_CONTROLLO,))]
        biologici = [x for x in valori if x not in negativi + positivi + dubbi]
        if negativi or positivi:
            # Una colonna che non lascia alcun valore ai biologici non descrive
            # il tipo dei campioni (un esito positivo o negativo, per esempio):
            # viene dopo quelle che ne lasciano.
            per_parole.append((not biologici, len(valori), posizione, TipoDiCampione(
                colonna, tuple(negativi), tuple(positivi), tuple(biologici), tuple(dubbi),
                "parole")))
        elif _contiene(colonna, v.PAROLE_TIPO) and len(valori) < max(2, len(righe)):
            per_nome.append((len(valori), posizione, TipoDiCampione(
                colonna, (), (), tuple(biologici), tuple(dubbi), "nome")))
    return [t for *_, t in sorted(per_parole, key=lambda x: x[:3])] or [
        t for *_, t in sorted(per_nome, key=lambda x: x[:2])]


def colonne_del_lotto(
    colonne: Sequence[str], righe: Sequence[dict[str, str]] = (), chiave: str | None = None
) -> dict[str, list[str]]:
    """Le colonne del file del lotto candidate a piastra, corsa e cellule, per
    le parole del loro nome; fra le piastre prima quelle di estrazione.

    Una piastra o una corsa raggruppa i campioni: non puo' essere la colonna
    della chiave ne' una colonna con un valore diverso per ogni riga. La
    colonna delle cellule ha la parola intera nel nome e, dove ha un valore, un
    numero.
    """
    def raggruppa(colonna: str) -> bool:
        valori = [r.get(colonna, "") for r in righe if r.get(colonna, "")]
        return colonna != chiave and (not righe or len(set(valori)) < len(valori))

    def numerica(colonna: str) -> bool:
        valori = [r.get(colonna, "") for r in righe if r.get(colonna, "")]
        return all(re.fullmatch(r"[0-9]+(\.[0-9]+)?([eE][+-]?[0-9]+)?", x) for x in valori)

    def nome(colonna: str) -> str:
        return colonna.replace("_", " ")

    piastre = [c for c in colonne if _contiene(nome(c), v.PAROLE_PIASTRA) and raggruppa(c)]
    piastre.sort(key=lambda c: not _contiene(nome(c), v.PAROLE_ESTRAZIONE))
    return {
        "piastra": piastre,
        "corsa": [c for c in colonne if _contiene(nome(c), v.PAROLE_CORSA) and raggruppa(c)],
        "cellule": [c for c in colonne
                    if c != chiave and _contiene(nome(c), v.PAROLE_CELLULE) and numerica(c)],
    }


_SEQUENZA = re.compile(r"(?<![A-Za-z])[ACGTRYSWKMBDHVN]{15,}(?![A-Za-z])")
_ETICHETTA_FORWARD = re.compile(
    r"(?<![A-Za-z])(?:F|FW|FWD|FORWARD|[0-9]+F)(?![A-Za-z])", re.IGNORECASE)
_ETICHETTA_REVERSE = re.compile(
    r"(?<![A-Za-z])(?:R|RV|REV|REVERSE|[0-9]+R)(?![A-Za-z])", re.IGNORECASE)


def primer_dai_metadati(
    tabelle: Iterable[Sequence[dict[str, str]]],
    primer_reverse: list[str] | None = None,
) -> tuple[str | None, list[str], bool]:
    """Il primer forward che i metadati riportano, tutte le sequenze trovate, e
    se il forward e' dichiarato tale.

    Si cercano nei valori le sequenze nucleotidiche di almeno 15 basi.
    L'etichetta di una sequenza e' nel testo che la precede nello stesso
    valore, oppure nel nome della colonna (``forward_primer``, ``515F``...). Il
    primer e' proposto se una sola sequenza ha un'etichetta di forward; oppure,
    senza alcuna etichetta, se la sequenza e' una sola, e allora il terzo
    valore e' falso: nessuno dice che sia il forward. Una sequenza con
    un'etichetta di reverse non si propone. Se ``primer_reverse`` e' un elenco
    vi si aggiungono le sequenze con un'etichetta di reverse: servono alla
    modalita' appaiata.
    """
    sequenze: list[str] = []
    forward: list[str] = []
    reverse: list[str] = []
    viste: set[tuple[str, str]] = set()
    for tabella in tabelle:
        for riga in tabella:
            for colonna, testo in riga.items():
                if not testo or (colonna, testo) in viste or not _SEQUENZA.search(testo):
                    continue
                viste.add((colonna, testo))
                nome = colonna.replace("_", " ")
                fine_precedente = 0
                for trovata in _SEQUENZA.finditer(testo):
                    sequenza = trovata.group(0)
                    prima = testo[fine_precedente:trovata.start()]
                    fine_precedente = trovata.end()
                    if sequenza not in sequenze:
                        sequenze.append(sequenza)
                    # L'etichetta accanto alla sequenza vale piu' del nome della colonna.
                    dove = prima if (_ETICHETTA_FORWARD.search(prima)
                                     or _ETICHETTA_REVERSE.search(prima)) else nome
                    if _ETICHETTA_FORWARD.search(dove) and not _ETICHETTA_REVERSE.search(dove):
                        if sequenza not in forward:
                            forward.append(sequenza)
                    elif _ETICHETTA_REVERSE.search(dove) and sequenza not in reverse:
                        reverse.append(sequenza)
    if primer_reverse is not None:
        primer_reverse.extend(reverse)
    if len(forward) == 1:
        return forward[0], sequenze, True
    if not forward and len(sequenze) == 1 and not reverse:
        return sequenze[0], sequenze, False
    return None, sequenze, False


def primer_nominati(
    tabelle: Iterable[Sequence[dict[str, str]]],
    elenco: Sequence[v.Primer] = v.PRIMER_16S,
) -> list[v.Primer]:
    """I primer dell'elenco il cui nome compare nei metadati (``515F``,
    ``341F``...): un nome non dice la sequenza, perche' lo stesso nome ha piu'
    forme, ma dice fra quali cercarla.
    """
    testi = {x for tabella in tabelle for riga in tabella for x in riga.values() if x}
    trovati = []
    for primer in elenco:
        nome = primer.nome.split()[0]
        cerca = re.compile(rf"(?<![A-Za-z0-9]){re.escape(nome)}(?![A-Za-z0-9])", re.IGNORECASE)
        if any(cerca.search(testo) for testo in testi):
            trovati.append(primer)
    return trovati


# --------------------------------------------------------------------------- #
# Riferimento tassonomico                                                      #
# --------------------------------------------------------------------------- #


def md5_di(percorso: Path) -> str:
    """L'MD5 del file, letto a blocchi."""
    impronta = hashlib.md5()
    with open(percorso, "rb") as file:
        while blocco := file.read(1 << 20):
            impronta.update(blocco)
    return impronta.hexdigest()


def nome_e_versione_del_riferimento(nome_file: str) -> tuple[str | None, str | None]:
    """Nome e versione del riferimento letti dal nome del file, se vi compaiono.

    Il nome e' uno di quelli noti; la versione il primo elemento numerico che
    lo segue, con o senza ``v`` davanti (``v138.1``, ``18``).
    """
    radice = re.sub(r"\.(fa|fasta|fna|fas)(\.gz)?$", "", nome_file, flags=re.IGNORECASE)
    elementi = [e for e in re.split(r"[_\-\s]+", radice) if e]
    for posizione, elemento in enumerate(elementi):
        nome = v.RIFERIMENTI.get(elemento.casefold())
        if nome is None:
            continue
        for seguente in elementi[posizione + 1:]:
            trovata = re.fullmatch(r"[vV]?([0-9]+(?:\.[0-9]+)*)", seguente)
            if trovata:
                return nome, trovata.group(1)
        return nome, None
    return None, None


def leggi_tabella(percorso: Path) -> tuple[list[str], list[dict[str, str]]]:
    """Intestazione e righe di una tabella, con il lettore della pipeline."""
    try:
        colonne, righe = intestazione(percorso), leggi_tsv(percorso)
    except (OSError, UnicodeDecodeError, csv.Error) as guasto:
        raise ValueError(f"{percorso} non si legge come tabella UTF-8 separata da tabulazioni: {guasto}") from guasto
    if not colonne or not righe:
        raise ValueError(
            f"{percorso} non ha " + ("un'intestazione" if not colonne else "alcuna riga di dati")
            + ": serve una tabella con una riga per campione")
    return colonne, righe


def valori_per_colonna(righe: Sequence[dict[str, str]], colonna: str) -> list[Any]:
    """I valori distinti di una colonna, in ordine."""
    return sorted({r.get(colonna, "") for r in righe} - {""})
