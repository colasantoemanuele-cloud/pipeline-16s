"""Ispezione delle prime letture di ciascun file, per i gate di validazione.

I gate che guardano le sequenze (struttura del file, lunghezze, presenza del
primer) devono costare minuti, non ore. Nessuno di essi legge un file intero:
si fermano alle prime ``qc.head_reads`` letture. L'integrità completa dei file
è già garantita altrove, dai checksum del manifesto.

La lettura avviene **una volta sola per file**: tutte le statistiche che
servono ai vari gate si ricavano dalla stessa passata, perché riaprire e
decomprimere lo stesso file tre volte triplicherebbe il costo senza aggiungere
nulla.

La scansione è **sequenziale**. Distribuirla su più processi è stato provato e
scartato: la passata sequenziale su un dataset di centinaia di file costa
secondi, e il vincolo da rispettare è che la validazione duri minuti e non ore. Il parallelismo faceva risparmiare pochi secondi in
cambio di una fragilità vera: con ``fork`` il pool può bloccarsi quando il
processo genitore ha sostituito i flussi standard, e con ``spawn`` ogni
chiamante sarebbe costretto a proteggere il proprio modulo principale con
``if __name__ == "__main__"``. Non è un prezzo che una libreria debba far
pagare a chi la usa.

**Limite dichiarato.** Le lunghezze osservate sono quelle delle prime letture,
non di tutto il file. Un file le cui letture più corte stiano oltre quelle
esaminate darebbe un minimo sovrastimato. È il prezzo del vincolo di costo, e
va tenuto presente da chi usa :attr:`StatisticheFile.lunghezza_minima`.
"""

from __future__ import annotations

import gzip
import re
import sys
import tempfile
import zlib
from collections import Counter
from collections.abc import Iterable, Iterator
from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path
from typing import IO, Final

__all__ = [
    "StatisticheFile",
    "apri_fastq",
    "conta_coppie",
    "espandi_iupac",
    "scansiona",
    "scansiona_file",
    "senza_righe_vuote_finali",
]

#: Corrispondenza fra codici IUPAC e classi di caratteri.
_IUPAC: Final[dict[str, str]] = {
    "A": "A", "C": "C", "G": "G", "T": "T",
    "R": "[AG]", "Y": "[CT]", "S": "[GC]", "W": "[AT]",
    "K": "[GT]", "M": "[AC]", "B": "[CGT]", "D": "[AGT]",
    "H": "[ACT]", "V": "[ACG]", "N": ".",
}


def espandi_iupac(sequenza: str) -> str:
    """Traduce una sequenza con codici degenerati in un'espressione regolare.

    ``GTGYCAGC`` diventa ``GTG[CT]CAGC``: cercare la sequenza alla lettera
    mancherebbe le letture che nelle posizioni degenerate portano l'altra
    base ammessa, e il primer risulterebbe assente anche quando c'è.
    """
    try:
        return "".join(_IUPAC[carattere] for carattere in sequenza.upper())
    except KeyError as guasto:
        raise ValueError(f"codice IUPAC sconosciuto: {guasto.args[0]!r}") from None


#: I primi due byte di un archivio gzip.
_MAGIA_GZIP: Final = b"\x1f\x8b"


def apri_fastq(percorso: Path | str) -> IO[str]:
    """Apre in lettura un FASTQ, compresso con gzip o no.

    Il formato si riconosce dai primi due byte, non dall'estensione: un file
    ``.fastq.gz`` non compresso, o un ``.fastq`` compresso, si leggono per
    quello che sono.
    """
    with open(percorso, "rb") as file:
        compresso = file.read(2) == _MAGIA_GZIP
    if compresso:
        return gzip.open(percorso, "rt", encoding="utf-8", errors="replace")
    return open(percorso, encoding="utf-8", errors="replace", newline="\n")


@dataclass(frozen=True)
class StatisticheFile:
    """Ciò che si è potuto misurare sulle prime letture di un file."""

    nome: str
    letture_esaminate: int
    lunghezza_minima: int | None
    lunghezza_massima: int | None
    con_primer: int
    con_motivo: int
    #: Descrizione del guasto strutturale, se il file non è leggibile come
    #: FASTQ. ``None`` se la struttura è valida.
    errore: str | None = None
    #: Vero se il file è finito prima di raggiungere il numero richiesto: non
    #: è un problema, ma dice che le statistiche coprono tutto il file.
    esaurito: bool = False
    #: Letture esaminate piu' corte della lunghezza richiesta alla scansione
    #: (``filter.truncLen``): il filtro le scarterebbe.
    piu_corte: int = 0
    #: Letture esaminate la cui intestazione le marca come prima o seconda
    #: lettura di una coppia (``/1``, ``/2``, o ``1:`` e ``2:`` nel commento), e
    #: nomi di lettura comparsi esattamente due volte nello stesso file, come
    #: le due letture di una coppia; e nomi comparsi piu' di due volte, che di
    #: una coppia non possono essere (un'intestazione vuota o costante).
    prime_di_coppia: int = 0
    seconde_di_coppia: int = 0
    nomi_ripetuti: int = 0
    nomi_oltre_due: int = 0
    #: Vero se dopo l'ultima lettura il file ha solo righe vuote: sono
    #: tollerate, ma i lettori di R non le ammettono
    #: (:func:`senza_righe_vuote_finali`).
    righe_vuote_finali: bool = False
    #: Quante delle coppie contate in ``nomi_ripetuti`` sono letture vicine con
    #: lo stesso nome a meno del marcatore di coppia (:class:`SerieDiCompagne`).
    compagne_per_marcatore: int = 0

    @property
    def coppie_nello_stesso_file(self) -> str | None:
        """Perche' il file non contiene le sole prime letture, o ``None``.

        Tre segni, giudicati sulla frazione delle letture esaminate e non su un
        record isolato (:data:`FRAZIONE_DI_COPPIA`): intestazioni che marcano
        sia la prima sia la seconda lettura di una coppia; nomi di lettura che
        compaiono due volte (con i nomi comparsi piu' di due volte sotto la
        stessa frazione: un'intestazione vuota o uguale per tutte le letture
        non e' un segno di coppia, ma un solo nome in tre copie non assolve un
        file di coppie); oppure quasi sole seconde letture, come nella testa di
        un file che riporta prima tutte le seconde e poi tutte le prime. Vale
        per le letture esaminate: G07 guarda le prime ``qc.head_reads``, e un
        file che riporta prima tutte le prime letture e poi le seconde gli
        sfugge se le esaminate non arrivano al secondo blocco; S1 ripete il
        giudizio su tutte le letture (:func:`conta_coppie`).
        """
        if not self.letture_esaminate:
            return None
        # Almeno due letture, oltre la frazione: in un file di poche letture
        # un solo record anomalo non deve bastare.
        tollerate = FRAZIONE_DI_COPPIA * self.letture_esaminate
        minimo = max(2.0, tollerate)
        if self.prime_di_coppia >= minimo and self.seconde_di_coppia >= minimo:
            return (
                f"{self.prime_di_coppia} intestazioni su {self.letture_esaminate} marcano "
                f"la prima lettura di una coppia e {self.seconde_di_coppia} la seconda"
            )
        # Lo stesso nome due volte e' il segno di una coppia; piu' di due volte
        # e' un'intestazione che non distingue le letture, e non dice nulla. I
        # nomi oltre le due copie si tollerano come ogni altro record anomalo:
        # sotto la frazione non spengono il riconoscimento.
        if self.nomi_ripetuti >= minimo and self.nomi_oltre_due < minimo:
            if self.compagne_per_marcatore:
                return (
                    f"{self.nomi_ripetuti} nomi di lettura su {self.letture_esaminate} letture "
                    f"compaiono due volte, {self.compagne_per_marcatore} dei quali in letture "
                    "consecutive che si distinguono per il solo marcatore di coppia in fondo "
                    "all'identificativo (1 e 2, oppure R1 e R2)"
                )
            return (
                f"{self.nomi_ripetuti} nomi di lettura su {self.letture_esaminate} letture "
                "compaiono due volte"
            )
        if (
            self.seconde_di_coppia >= minimo
            and self.prime_di_coppia < minimo
            # Almeno un record tollerato anche nei file piccoli, dove la
            # frazione vale meno di una lettura: come per le altre soglie, un
            # record isolato non decide in nessuno dei due versi.
            and self.seconde_di_coppia >= self.letture_esaminate - max(1.0, tollerate)
        ):
            return (
                f"{self.seconde_di_coppia} intestazioni su {self.letture_esaminate} marcano "
                f"la seconda lettura di una coppia, e {self.prime_di_coppia} la prima"
            )
        return None

    @property
    def frazione_corte(self) -> float:
        """La frazione delle letture esaminate piu' corte della lunghezza richiesta."""
        return self.piu_corte / self.letture_esaminate if self.letture_esaminate else 0.0

    @property
    def valido(self) -> bool:
        """Vero se il file si è letto senza errori."""
        return self.errore is None

    @property
    def frazione_primer(self) -> float:
        """La frazione delle letture esaminate che iniziano con il primer; 0 se nessuna
        è stata esaminata.
        """
        return self.con_primer / self.letture_esaminate if self.letture_esaminate else 0.0

    @property
    def frazione_motivo(self) -> float:
        """La frazione delle letture esaminate che contengono il motivo conservato; 0 se
        nessuna è stata esaminata.
        """
        return self.con_motivo / self.letture_esaminate if self.letture_esaminate else 0.0


#: La frazione delle letture esaminate oltre la quale un segno di coppia conta:
#: un record isolato con un marcatore o un nome ripetuto e' un difetto del
#: file, non la prova che contenga le due letture di ogni coppia. E' una
#: tolleranza al rumore, non una proprieta' di un dataset: in un file con le
#: coppie i segni riguardano circa meta' delle letture.
FRAZIONE_DI_COPPIA: Final = 0.05

#: La lettura di una coppia dichiarata dall'intestazione, nelle due convenzioni
#: diffuse: il suffisso ``/1`` o ``/2`` in fondo al nome o all'intestazione, e
#: il commento ``1:N:0:...`` o ``2:N:0:...`` dei sequenziatori Illumina.
_SUFFISSO_COPPIA: Final = re.compile(r"/([12])$")
_COMMENTO_COPPIA: Final = re.compile(r"^([12]):[YN]:")


def _nome_e_coppia(intestazione: str) -> tuple[str, str | None]:
    """L'intestazione della lettura senza il marcatore di coppia, e il marcatore
    (``"1"``, ``"2"`` o ``None``).

    Il nome con cui si riconosce una ripetizione e' l'intera intestazione, non
    il suo primo campo: alcuni file ripetono in testa il nome del campione e
    numerano la lettura nel campo seguente. Tolto il marcatore, le due letture
    di una coppia hanno la stessa intestazione.
    """
    riga = intestazione[1:].rstrip("\r\n")
    nome, separatore, commento = riga.partition(" ")
    for testo, resto in ((nome, separatore + commento), (riga, "")):
        trovato = _SUFFISSO_COPPIA.search(testo)
        if trovato:
            return _SUFFISSO_COPPIA.sub("", testo) + resto, trovato.group(1)
    trovato = _COMMENTO_COPPIA.match(commento)
    if trovato:
        return f"{nome} {commento[1:]}", trovato.group(1)
    return riga, None


#: Il marcatore di coppia in fondo all'identificativo della lettura: l'ultimo
#: elemento, quando vale esattamente 1, 2, R1 o R2 ed e' preceduto da un
#: separatore. Il separatore e' obbligatorio: in ``x.11`` l'ultimo elemento e'
#: 11, e non c'e' marcatore.
_MARCATORE_IN_FONDO: Final = re.compile(r"^(.+)[/._-](R?[12])$")
_PRIMA_E_SECONDA: Final = frozenset({("1", "2"), ("R1", "R2")})


def marcatore_di_coppia(intestazione: str) -> tuple[str, str] | None:
    """Il prefisso e il marcatore di coppia dell'identificativo di una lettura,
    o ``None`` se l'identificativo non finisce con un marcatore.

    L'identificativo e' il primo campo dell'intestazione, senza ``@``; il
    prefisso e' l'identificativo senza il separatore e il marcatore.
    """
    campi = intestazione[1:].split(None, 1)
    trovato = _MARCATORE_IN_FONDO.match(campi[0]) if campi else None
    return (trovato.group(1), trovato.group(2)) if trovato else None


class SerieDiCompagne:
    """Conta le coppie di letture vicine che si distinguono per il solo
    marcatore di coppia: lo stesso nome ripetuto, a meno del marcatore.

    Il suffisso da solo non basta. Gli archivi numerano le letture di un file
    single-end in fondo all'identificativo (``x.1``, ``x.2``, ``x.3``, ...), e
    le prime due avrebbero lo stesso prefisso e i marcatori 1 e 2 senza essere
    una coppia. Si guarda quindi la struttura:

    * due letture consecutive in posizione dispari e pari (1-2, 3-4, ...) sono
      compagne se hanno lo stesso prefisso e i marcatori 1 poi 2, oppure R1 poi
      R2 (:func:`marcatore_di_coppia`);
    * contano solo le compagne che stanno in una serie di almeno due coppie di
      compagne consecutive con prefissi diversi. In un file single-end numerato
      la sola coppia 1-2 e' isolata e non conta; in un file con le due letture
      di ogni coppia intercalate contano quasi tutte.

    Si confrontano solo letture vicine: la memoria non cresce con il file. Le
    due letture il cui nome e' gia' uguale senza il marcatore (``/1`` e ``/2``)
    non si contano qui: sono gia' fra i nomi ripetuti.
    """

    def __init__(self) -> None:
        #: Le coppie di compagne contate.
        self.contate = 0
        self._dispari: tuple[tuple[str, str] | None, str] | None = None
        self._in_serie = 0
        self._prefisso_precedente: str | None = None

    def aggiungi(self, intestazione: str, nome: str) -> None:
        """Esamina la lettura successiva; ``nome`` e' quello con cui si
        riconoscono i nomi ripetuti (:func:`_nome_e_coppia`).
        """
        marcatore = marcatore_di_coppia(intestazione)
        if self._dispari is None:
            self._dispari = (marcatore, nome)
            return
        (prima, nome_della_prima), self._dispari = self._dispari, None
        compagne = (
            prima is not None and marcatore is not None
            and prima[0] == marcatore[0]
            and (prima[1], marcatore[1]) in _PRIMA_E_SECONDA
            and nome_della_prima != nome
        )
        if not compagne:
            self._in_serie, self._prefisso_precedente = 0, None
            return
        if self._in_serie and prima[0] != self._prefisso_precedente:
            self._in_serie += 1
            # Alla seconda coppia della serie conta anche la prima.
            self.contate += 2 if self._in_serie == 2 else 1
        else:
            self._in_serie = 1
        self._prefisso_precedente = prima[0]


def _solo_righe_vuote(file: IO[str]) -> bool:
    """Vero se fino alla fine del file restano solo righe vuote."""
    return all(not riga.strip() for riga in file)


def scansiona_file(
    percorso: Path | str,
    head_reads: int,
    primer: str | None,
    motivo: str | None,
    inizio_motivo: int = 0,
    lunghezza_richiesta: int = 0,
) -> StatisticheFile:
    """Legge le prime ``head_reads`` letture e ne ricava tutte le statistiche.

    ``lunghezza_richiesta`` e' la lunghezza sotto la quale una lettura si conta
    fra le piu' corte (``filter.truncLen``); zero non ne conta nessuna.

    ``primer`` e ``motivo`` sono espressioni regolari. Il primer si cerca
    **ancorato a inizio lettura**: starebbe in testa se non fosse stato tolto.
    Il motivo conservato, che apre la regione amplificata, si cerca ancorato
    alla posizione ``inizio_motivo``: zero se le letture iniziano dalla
    regione, la lunghezza di cio' che la precede (``filter.trimLeft``) se
    portano ancora il primer che il filtro togliera'. Con ``None`` la ricerca
    corrispondente non si fa e il conteggio resta zero.
    """
    percorso = Path(percorso)
    espressione_primer = re.compile(primer) if primer is not None else None
    espressione_motivo = re.compile(motivo) if motivo is not None else None

    letture = con_primer = con_motivo = corte = 0
    prime = seconde = 0
    nomi: Counter[str] = Counter()
    compagne = SerieDiCompagne()
    minima: int | None = None
    massima: int | None = None
    esaurito = True
    vuote_finali = False

    try:
        with apri_fastq(percorso) as file:
            while letture < head_reads:
                intestazione = file.readline()
                if not intestazione:
                    break  # fine del file su un confine di record: corretto
                if not intestazione.strip():
                    # Una riga vuota al posto di un'intestazione: in fondo al
                    # file e' solo spaziatura (alcuni strumenti chiudono cosi'
                    # i file), in mezzo spezza i record e chi legge a blocchi
                    # di quattro righe li sfaserebbe tutti.
                    if _solo_righe_vuote(file):
                        vuote_finali = True
                        break
                    return StatisticheFile(
                        percorso.name, letture, minima, massima, con_primer, con_motivo,
                        errore=f"riga vuota in mezzo al file dopo {letture} letture: "
                               f"seguono altre righe",
                    )
                sequenza = file.readline()
                separatore = file.readline()
                qualita = file.readline()

                if not (sequenza and separatore and qualita):
                    return StatisticheFile(
                        percorso.name, letture, minima, massima, con_primer, con_motivo,
                        errore=f"record troncato dopo {letture} letture: "
                               f"il file non ha quattro righe per record",
                    )
                if not intestazione.startswith("@"):
                    return StatisticheFile(
                        percorso.name, letture, minima, massima, con_primer, con_motivo,
                        errore=f"alla lettura {letture + 1} l'intestazione non comincia "
                               f"con '@': {intestazione[:40]!r}",
                    )
                if not separatore.startswith("+"):
                    return StatisticheFile(
                        percorso.name, letture, minima, massima, con_primer, con_motivo,
                        errore=f"alla lettura {letture + 1} la terza riga non comincia "
                               f"con '+': {separatore[:40]!r}",
                    )

                sequenza = sequenza.rstrip("\r\n")
                if len(sequenza) != len(qualita.rstrip("\r\n")):
                    return StatisticheFile(
                        percorso.name, letture, minima, massima, con_primer, con_motivo,
                        errore=f"alla lettura {letture + 1} sequenza e qualita' hanno "
                               f"lunghezze diverse",
                    )

                letture += 1
                nome, coppia = _nome_e_coppia(intestazione)
                prime += coppia == "1"
                seconde += coppia == "2"
                nomi[nome] += 1
                compagne.aggiungi(intestazione, nome)
                lunghezza = len(sequenza)
                minima = lunghezza if minima is None else min(minima, lunghezza)
                massima = lunghezza if massima is None else max(massima, lunghezza)
                if lunghezza < lunghezza_richiesta:
                    corte += 1
                if espressione_primer is not None and espressione_primer.match(sequenza):
                    con_primer += 1
                # Sulla lettura tagliata, non dalla posizione del taglio: un
                # motivo ancorato con ^ deve valere all'inizio di cio' che resta.
                if espressione_motivo is not None and espressione_motivo.match(
                    sequenza[inizio_motivo:]
                ):
                    con_motivo += 1
            else:
                esaurito = False
    except (OSError, EOFError, zlib.error) as guasto:
        # Ogni guasto di decompressione: archivio interrotto (EOFError), flusso
        # o intestazione gzip corrotti (gzip.BadGzipFile, che e' un OSError, e
        # zlib.error, che non lo e'), errore di lettura dal disco.
        return StatisticheFile(
            percorso.name, letture, minima, massima, con_primer, con_motivo,
            errore=f"archivio non leggibile dopo {letture} letture "
                   f"({type(guasto).__name__}: {guasto})",
        )

    if letture == 0:
        return StatisticheFile(
            percorso.name, 0, None, None, 0, 0,
            errore="il file non contiene alcuna lettura",
        )

    return StatisticheFile(
        percorso.name, letture, minima, massima, con_primer, con_motivo,
        esaurito=esaurito, piu_corte=corte,
        prime_di_coppia=prime, seconde_di_coppia=seconde,
        # Lo stesso nome a meno del marcatore di coppia e' un nome ripetuto: le
        # compagne riconosciute dal marcatore si sommano ai nomi comparsi due
        # volte, e ne seguono soglia e giudizio.
        nomi_ripetuti=sum(1 for volte in nomi.values() if volte == 2) + compagne.contate,
        nomi_oltre_due=sum(1 for volte in nomi.values() if volte > 2),
        righe_vuote_finali=vuote_finali,
        compagne_per_marcatore=compagne.contate,
    )


@contextmanager
def senza_righe_vuote_finali(
    percorsi: dict[str, str], con_righe_vuote: Iterable[str]
) -> Iterator[dict[str, str]]:
    """I percorsi da dare ai lettori di R, con una copia ripulita dei file indicati.

    Le righe vuote in fondo a un FASTQ sono tollerate, ma i lettori a blocchi
    di R (``ShortRead``, quindi anche ``dada2::filterAndTrim``) le rifiutano
    come un record che non comincia con ``@``. Dei soli file in
    ``con_righe_vuote`` (chiavi di ``percorsi``) si scrive una copia senza
    quelle righe in una cartella temporanea, che vive quanto il blocco
    ``with``: i file di ingresso non si toccano mai, e senza file da ripulire
    non si scrive nulla.
    """
    chiavi = sorted(con_righe_vuote)
    if not chiavi:
        yield percorsi
        return
    with tempfile.TemporaryDirectory(prefix="amplicon16s-letture-") as cartella:
        ripuliti = dict(percorsi)
        for numero, chiave in enumerate(chiavi):
            # Una sottocartella per file: due file di cartelle diverse possono
            # avere lo stesso nome.
            destinazione = Path(cartella) / str(numero) / Path(percorsi[chiave]).name
            destinazione.parent.mkdir()
            _copia_senza_righe_vuote_finali(percorsi[chiave], destinazione)
            ripuliti[chiave] = str(destinazione)
        yield ripuliti


def _copia_senza_righe_vuote_finali(origine: Path | str, destinazione: Path) -> None:
    """Copia un FASTQ, compresso, senza le righe vuote che lo chiudono."""
    with apri_fastq(origine) as ingresso, gzip.open(
        destinazione, "wt", encoding="utf-8", newline="\n", compresslevel=1
    ) as uscita:
        for riga in ingresso:
            if not riga.strip():
                # La scansione ha gia' stabilito che da qui in poi sono vuote.
                break
            uscita.write(riga)


def conta_coppie(percorso: Path | str) -> StatisticheFile:
    """I segni di coppia su **tutte** le letture del file, con le regole di G07.

    E' la stessa scansione dei gate senza il limite delle prime letture e senza
    le ricerche di primer e motivo: il giudizio
    (:attr:`StatisticheFile.coppie_nello_stesso_file`) e' lo stesso, e cio' che
    G07 non vede oltre ``qc.head_reads`` qui si vede. I nomi di lettura di un
    file restano in memoria solo per la durata della sua scansione. Dice anche
    se il file si chiude con righe vuote (``righe_vuote_finali``).
    """
    return scansiona_file(percorso, sys.maxsize, None, None)


def scansiona(
    percorsi: dict[str, Path],
    head_reads: int,
    primer: str | None,
    motivo: str | None,
    inizio_motivo: int = 0,
    lunghezza_richiesta: int = 0,
) -> dict[str, StatisticheFile]:
    """Scansiona i file indicati e restituisce le statistiche di ciascuno."""
    return {
        chiave: scansiona_file(percorso, head_reads, primer, motivo, inizio_motivo,
                               lunghezza_richiesta)
        for chiave, percorso in percorsi.items()
    }
