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
import shutil
import sys
import tempfile
import zlib
from collections import Counter
from collections.abc import Iterable, Iterator
from contextlib import contextmanager
from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path
from typing import IO, Any, Final

__all__ = [
    "Abbinamento",
    "LETTURE_MINIME_PER_LA_QUOTA",
    "LETTURE_MINIME_PER_UN_NOME_SOSPETTO",
    "MARCATORI_DI_FILE",
    "QUOTA_MINIMA_DEL_SEGNALE",
    "RAPPORTO_MASSIMO_FRA_LE_PARTI",
    "SEGNALE_MINIMO_DI_UNA_PARTE",
    "TRATTI",
    "FRAZIONE_DI_ERRORI",
    "QUOTA_FUORI_POSIZIONE",
    "SFASAMENTO_CERCATO",
    "RicercaPrimer",
    "discordanze_di_partenza",
    "sonda_del_taglio_fisso",
    "ricerca_del_primer",
    "SegnaleDiOrientamento",
    "Sincronia",
    "StatisticheFile",
    "abbina_file",
    "apri_fastq",
    "compagni_per_nome",
    "corsa_invertita",
    "conta_coppie",
    "espandi_iupac",
    "nome_da_seconda_lettura",
    "scansiona",
    "scansiona_file",
    "segnale_di_orientamento",
    "senza_righe_vuote_finali",
    "verifica_sincronia",
    "versi_mescolati",
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


#: Le basi che ciascun codice IUPAC del primer ammette nella lettura. Le basi
#: della lettura si prendono alla lettera: una N della lettura non e' fra
#: quelle ammesse da alcun codice, nemmeno dalla N del primer.
_BASI_AMMESSE: Final[dict[str, str]] = {
    codice: (classe.strip("[]") if codice != "N" else "ACGT")
    for codice, classe in _IUPAC.items()
}


@dataclass(frozen=True)
class RicercaPrimer:
    """La ricerca del primer nelle prime basi di una lettura, per i protocolli
    che lo fanno precedere da un numero variabile di basi.

    Il primer (codici IUPAC) deve cominciare a uno sfasamento fra zero e
    ``sfasamento_massimo`` e stare per intero nella lettura, con al piu'
    ``discordanze_massime`` basi diverse: sole sostituzioni, nessuna inserzione
    o delezione, e una N della lettura e' una discordanza. Fra piu'
    corrispondenze vale quella con meno discordanze e, a parita', la piu' a
    sinistra. E' la stessa regola di ``R/lib/primer.R``, che taglia i file
    interi: qui serve a validazione, segnale di orientamento e ``init``.

    Razionale sistemistico: senza inserzioni ne' delezioni la posizione del
    taglio e' determinata senza ambiguita', e ogni lettura tagliata comincia
    alla stessa base della regione amplificata.
    """

    primer: str
    discordanze_massime: int
    sfasamento_massimo: int
    #: Vero se una lettura senza il primer viene scartata dal filtro.
    scarta_senza_primer: bool = True

    def __post_init__(self) -> None:
        primer = self.primer.upper()
        try:
            ammesse = tuple(frozenset(_BASI_AMMESSE[c]) for c in primer)
        except KeyError as guasto:
            raise ValueError(f"codice IUPAC sconosciuto: {guasto.args[0]!r}") from None
        object.__setattr__(self, "_ammesse", ammesse)
        # La corrispondenza senza discordanze piu' a sinistra, in un passo solo:
        # e' il caso di quasi tutte le letture.
        object.__setattr__(self, "_esatta", re.compile(
            ".{0,%d}?(%s)" % (
                self.sfasamento_massimo,
                "".join(f"[{_BASI_AMMESSE[c]}]" for c in primer))))

    def trova(self, sequenza: str) -> tuple[int, int] | None:
        """Sfasamento e discordanze della corrispondenza scelta, o ``None``.

        Le basi si confrontano in maiuscolo, come le legge R: una lettura
        scritta in minuscolo porta il primer quanto una in maiuscolo.
        """
        sequenza = sequenza.upper()
        esatta = self._esatta.match(sequenza)
        if esatta is not None:
            return esatta.start(1), 0
        if self.discordanze_massime == 0:
            return None
        ammesse = self._ammesse
        lunghezza = len(ammesse)
        migliore: tuple[int, int] | None = None
        tetto = self.discordanze_massime
        for sfasamento in range(min(self.sfasamento_massimo, len(sequenza) - lunghezza) + 1):
            discordanze = 0
            for base, classe in zip(sequenza[sfasamento:sfasamento + lunghezza], ammesse):
                if base not in classe:
                    discordanze += 1
                    if discordanze > tetto:
                        break
            else:
                # Solo una corrispondenza con meno discordanze sostituisce una
                # piu' a sinistra.
                migliore = (sfasamento, discordanze)
                tetto = discordanze - 1
                if tetto < 1:
                    # Una senza discordanze l'avrebbe trovata la ricerca esatta.
                    break
        return migliore

    def fine(self, sequenza: str) -> int | None:
        """La posizione dopo l'ultima base del primer: dove la lettura si
        taglia; ``None`` se il primer non c'e'.
        """
        trovata = self.trova(sequenza)
        return None if trovata is None else trovata[0] + len(self.primer)


#: Fino a quante basi di sfasamento si cerca un primer di cui non e' dichiarata
#: la finestra (``init``, e la validazione con il taglio a posizione fissa): il doppio delle
#: otto lunghezze, da zero a sette basi, dei distanziatori di Fadrosh et al.
#: 2014 (Microbiome 2:6), cosi' che un protocollo con distanziatori piu' lunghi
#: di quelli pubblicati si riconosca comunque.
SFASAMENTO_CERCATO = 16

#: La frazione di errori ammessa fra primer e lettura, da cui il numero di
#: discordanze di partenza: il predefinito di cutadapt (``-e 0.1``,
#: documentazione di cutadapt, «Error tolerance»), arrotondato per difetto
#: sulla lunghezza del primer.
FRAZIONE_DI_ERRORI = 0.1

#: La quota di letture con il primer che deve portarlo a uno sfasamento diverso
#: dal piu' frequente perche' il primer vada cercato: ``init`` propone allora la ricerca, e la
#: validazione avvisa se e' dichiarato il taglio a posizione fissa. Un errore di
#: inserzione o delezione sposta il primer in ben meno di una lettura su cento;
#: un protocollo con distanziatori distribuisce le letture su piu' sfasamenti
#: con quote di decine di punti.
QUOTA_FUORI_POSIZIONE = 0.05

def discordanze_di_partenza(*primer: str) -> int:
    """Le discordanze ammesse di partenza per i primer dati: la frazione di
    errori di cutadapt sulla lunghezza del piu' corto.
    """
    return int(FRAZIONE_DI_ERRORI * min(len(p) for p in primer))



def sonda_del_taglio_fisso(primer: str | None, taglio: int) -> RicercaPrimer | None:
    """La ricerca con cui la validazione verifica un taglio a posizione fissa:
    il primer dichiarato, cercato lungo la lettura con i valori di partenza
    della ricerca. ``None`` se la verifica non si applica: primer nullo, nessun
    taglio, o un taglio piu' corto del primer (che non lo toglierebbe intero).

    La posizione attesa del primer e' ``taglio`` meno la sua lunghezza; la
    finestra arriva :data:`SFASAMENTO_CERCATO` basi oltre, cosi' che si veda un
    primer sia prima sia dopo quella posizione.
    """
    if primer is None or taglio <= 0 or taglio < len(primer):
        return None
    return RicercaPrimer(
        primer, discordanze_di_partenza(primer), taglio - len(primer) + SFASAMENTO_CERCATO)


def ricerca_del_primer(config: Any, reverse: bool = False) -> RicercaPrimer | None:
    """La ricerca del primer che la configurazione dichiara per le letture
    forward, o per le reverse; ``None`` con il taglio a posizione fissa.

    ``config`` e' la configurazione, o la vista che ne ha una fase: si leggono
    il primer, i tre parametri del taglio cercato e, per le reverse, il primer
    reverse.
    """
    if not config.taglio_cercato:
        return None
    primer = config.pair.primer_sequence_reverse if reverse else config.qc.primer_sequence
    return RicercaPrimer(
        primer, config.filter.primer_max_mismatch, config.filter.primer_max_offset,
        scarta_senza_primer=config.filter.primer_untrimmed == "discard")


@dataclass(frozen=True)
class SegnaleDiOrientamento:
    """Che cosa distingue le letture dell'orientamento atteso in un file.

    Una lettura porta il segnale se il motivo corrisponde alla posizione
    ``taglio``, oppure se il primer e' in testa: le due condizioni si uniscono
    per lettura, perche' una lettura con un errore nel primer e il motivo
    integro (o viceversa) e' dell'orientamento atteso quanto una che li ha
    entrambi, e una lettura dell'altro orientamento non ha ne' l'uno ne'
    l'altro. ``primer`` e ``motivo`` sono espressioni regolari; una delle due
    puo' mancare.

    Vale per qualunque estremita' dell'amplicone: con il primer forward e il
    motivo con cui inizia la regione descrive un file di letture forward, con
    il primer inverso (e il motivo dell'altra estremita', se noto) un file di
    letture inverse.
    """

    primer: str | None
    motivo: str | None
    taglio: int = 0
    #: Con il primer cercato lungo la lettura: il segnale e' il primer trovato
    #: con la regola di ricerca, in unione con il motivo dopo di esso.
    ricerca: RicercaPrimer | None = None

    def lo_porta(self, sequenza: str, in_testa: Any, dopo_il_taglio: Any) -> bool:
        """Se la lettura porta il segnale; ``in_testa`` e ``dopo_il_taglio``
        sono le espressioni compilate di ``primer`` e ``motivo``.
        """
        if self.ricerca is not None:
            fine = self.ricerca.fine(sequenza)
            return fine is not None or bool(
                dopo_il_taglio is not None and dopo_il_taglio.match(sequenza))
        return bool(
            (dopo_il_taglio is not None and dopo_il_taglio.match(sequenza[self.taglio:]))
            or (in_testa is not None and in_testa.match(sequenza)))


def segnale_di_orientamento(
    primer: str | None, motivo: str | None, taglio: int,
    ricerca: RicercaPrimer | None = None,
) -> SegnaleDiOrientamento | None:
    """Il segnale da cercare nelle letture, dati il primer (in codici IUPAC),
    il motivo conservato e le basi tolte in testa; ``None`` se non ce n'e'.

    Il primer si cerca solo se il taglio lo toglie: con il taglio a zero il
    primer deve essere assente dalle letture, e il segnale e' il solo motivo.
    Senza motivo e senza un primer da cercare nelle basi tagliate non c'e'
    nulla che distingua un orientamento dall'altro.
    """
    if ricerca is not None:
        # Il primer cercato lungo la lettura e' esso stesso il segnale; il
        # motivo, che si trova dopo il primer, vale per le letture in cui il
        # primer non c'e' e la regione comincia alla prima base.
        return SegnaleDiOrientamento(None, motivo, 0, ricerca)
    in_testa = espandi_iupac(primer) if primer is not None and taglio > 0 else None
    if in_testa is None and motivo is None:
        return None
    return SegnaleDiOrientamento(in_testa, motivo, taglio)


#: Le letture con cui la quota del segnale di un file si giudica: con meno,
#: zero letture con il segnale non escludono una quota vera dell'uno per cento
#: (0,99 alla 500 e' sotto l'uno per cento).
LETTURE_MINIME_PER_LA_QUOTA: Final = 500

#: La quota del segnale sotto la quale un file non e' dell'orientamento atteso.
#: In una lettura dell'altro orientamento il motivo compare solo per caso (un
#: motivo con quindici basi determinate corrisponde a una sequenza casuale una
#: volta su 4 alla 15), quindi la quota attesa e' zero: uno per mille tollera
#: le letture attribuite al file sbagliato. Un file dell'orientamento atteso,
#: anche dominato da sequenze fuori bersaglio, conserva una quota di ordini di
#: grandezza superiore.
QUOTA_MINIMA_DEL_SEGNALE: Final = Fraction(1, 1000)

#: Le letture con cui la quota si giudica in un file che il nome rende sospetto
#: (finisce come il secondo file di una coppia): con 50 letture e nessuna con il
#: segnale si esclude una quota vera del dieci per cento (0,9 alla 50 e' sotto
#: l'uno per cento). A un file gia' sospetto non serve la prova piu' forte
#: chiesta a un file qualunque.
LETTURE_MINIME_PER_UN_NOME_SOSPETTO: Final = 50

#: Il contrasto fra le parti di un file che ne rivela la struttura a blocchi
#: (tutte le prime letture, poi tutte le seconde, anche in piu' blocchi o di
#: lunghezza diversa) o intercalata (una lettura per verso): una parte ha
#: almeno SEGNALE_MINIMO_DI_UNA_PARTE letture con il segnale, e la quota di
#: un'altra e' sotto RAPPORTO_MASSIMO_FRA_LE_PARTI della sua. Le parti sono i
#: TRATTI tratti consecutivi di uguale lunghezza in cui si divide il file, e le
#: letture di posto dispari e pari. In un file di un solo orientamento, letto
#: nell'ordine del sequenziatore, le parti hanno la stessa composizione: con 30
#: letture in una e meno di 3 attese nell'altra, la probabilita' che si dividano
#: cosi' per caso e' sotto una su un milione. Un rapporto, e non una soglia di
#: significativita': su un file profondo anche una differenza piccola fra le
#: parti e' significativa, e non indica nulla. Dieci tratti: abbastanza da
#: vedere blocchi di lunghezza diversa o ripetuti, pochi abbastanza perche' un
#: tratto di un file di poche migliaia di letture ne abbia a sufficienza.
SEGNALE_MINIMO_DI_UNA_PARTE: Final = 30
RAPPORTO_MASSIMO_FRA_LE_PARTI: Final = Fraction(1, 10)
TRATTI: Final = 10


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
    #: Vero se il file si legge senza guasti e non contiene alcuna lettura.
    #: Resta non valido (``errore`` lo dice), ma chi conosce la classe del
    #: campione lo distingue da un file guasto: per un controllo negativo
    #: l'assenza di letture e' un esito legittimo.
    vuoto: bool = False
    #: Le letture con il segnale dell'orientamento atteso in ogni tratto
    #: consecutivo del file, ciascuno come (letture con il segnale, letture);
    #: ``None`` se il segnale non e' stato cercato.
    segnale_per_tratto: tuple[tuple[int, int], ...] | None = None
    #: Lo stesso per le letture di posto dispari e per quelle di posto pari.
    segnale_alterno: tuple[tuple[int, int], tuple[int, int]] | None = None
    #: Le letture con il segnale dell'altro orientamento, quando la scansione lo
    #: ha cercato (modalita' appaiata); ``None`` altrimenti.
    con_altro_segnale: int | None = None
    #: Con una sonda del primer: a quale sfasamento il primer e' stato trovato,
    #: e in quante letture. Serve al giudizio della validazione sul taglio a
    #: posizione fissa, e non compare in alcuna tabella.
    sfasamenti_del_primer: tuple[tuple[int, int], ...] | None = None

    @property
    def con_segnale(self) -> int | None:
        """Le letture con il segnale dell'orientamento atteso, o ``None`` se
        non e' stato cercato."""
        if self.segnale_per_tratto is None:
            return None
        return sum(k for k, _ in self.segnale_per_tratto)

    def segni_di_coppia_inattesi(self, attesa: str) -> tuple[str, str] | None:
        """In modalita' appaiata, che cosa dicono le intestazioni di un file che
        deve contenere le sole letture ``attesa`` (``"1"`` o ``"2"``) di ogni
        coppia: ``None`` se nulla di inatteso, altrimenti la specie del difetto
        e la sua descrizione.

        ``"intercalato"``: il file porta i segni di entrambe le letture (i due
        marcatori, o lo stesso nome due volte): contiene le coppie intere.
        ``"invertito"``: il file porta quasi solo il marcatore dell'altra
        lettura. Le soglie sono quelle di
        :attr:`coppie_nello_stesso_file`: un record isolato non decide.
        """
        if not self.letture_esaminate:
            return None
        tollerate = FRAZIONE_DI_COPPIA * self.letture_esaminate
        minimo = max(2.0, tollerate)
        proprie, altre = (
            (self.prime_di_coppia, self.seconde_di_coppia) if attesa == "1"
            else (self.seconde_di_coppia, self.prime_di_coppia))
        if proprie >= minimo and altre >= minimo:
            return "intercalato", (
                f"{self.prime_di_coppia} intestazioni su {self.letture_esaminate} marcano "
                f"la prima lettura di una coppia e {self.seconde_di_coppia} la seconda")
        if self.nomi_ripetuti >= minimo and self.nomi_oltre_due < minimo:
            return "intercalato", (
                f"{self.nomi_ripetuti} nomi di lettura su {self.letture_esaminate} letture "
                "compaiono due volte")
        if (altre >= minimo and proprie < minimo
                and altre >= self.letture_esaminate - max(1.0, tollerate)):
            quale = "seconda" if attesa == "1" else "prima"
            return "invertito", (
                f"{altre} intestazioni su {self.letture_esaminate} marcano la {quale} "
                "lettura di una coppia")
        return None

    def orientamento_inatteso(
        self, giudica_la_quota: bool = True,
        letture_minime: int = LETTURE_MINIME_PER_LA_QUOTA,
    ) -> str | None:
        """Perche' il file non contiene le sole letture dell'orientamento
        atteso, o ``None``.

        Due giudizi, sulle letture esaminate (tutte, per essere validi:
        :func:`conta_coppie`). Il primo riconosce un file di sole letture
        dell'altro orientamento: la quota delle letture con il segnale e' sotto
        :data:`QUOTA_MINIMA_DEL_SEGNALE`, su almeno
        :data:`LETTURE_MINIME_PER_LA_QUOTA` letture. Si chiede con
        ``giudica_la_quota``: da un controllo negativo non ci si attende il
        bersaglio, e un bianco con soli prodotti aspecifici ha legittimamente
        quota zero. Il secondo riconosce, in un file di qualunque classe, le
        letture dei due orientamenti a blocchi o intercalate: una parte del
        file ha il segnale e un'altra quasi nulla
        (:data:`SEGNALE_MINIMO_DI_UNA_PARTE`,
        :data:`RAPPORTO_MASSIMO_FRA_LE_PARTI`). ``letture_minime`` abbassa il
        minimo del primo giudizio per un file che il nome rende sospetto
        (:data:`LETTURE_MINIME_PER_UN_NOME_SOSPETTO`).
        """
        if self.segnale_per_tratto is None:
            return None
        letture = sum(n for _, n in self.segnale_per_tratto)
        con_segnale = sum(k for k, _ in self.segnale_per_tratto)
        if (giudica_la_quota and letture >= letture_minime
                and con_segnale < QUOTA_MINIMA_DEL_SEGNALE * letture):
            return (
                f"{con_segnale} letture su {letture} portano il segnale dell'orientamento "
                f"atteso (quota {con_segnale / letture:.4f}, sotto "
                f"{float(QUOTA_MINIMA_DEL_SEGNALE)})"
            )
        tratti = len(self.segnale_per_tratto)
        suddivisioni = [(
            self.segnale_per_tratto,
            [f"tratto {i + 1} di {tratti}" for i in range(tratti)],
            "il file riporta a blocchi letture dei due orientamenti")]
        if self.segnale_alterno is not None:
            suddivisioni.append((
                self.segnale_alterno, ["letture di posto dispari", "letture di posto pari"],
                "il file alterna letture dei due orientamenti"))
        for parti, nomi, conclusione in suddivisioni:
            # Frazioni esatte: al confine del rapporto un errore di virgola
            # mobile deciderebbe il giudizio.
            quote = [Fraction(k, n) if n else None for k, n in parti]
            candidate = [i for i, (k, _) in enumerate(parti) if k >= SEGNALE_MINIMO_DI_UNA_PARTE]
            presenti = [i for i, q in enumerate(quote) if q is not None]
            if not candidate or len(presenti) < 2:
                continue
            forte = max(candidate, key=lambda i: quote[i])
            debole = min(presenti, key=lambda i: quote[i])
            if quote[debole] < RAPPORTO_MASSIMO_FRA_LE_PARTI * quote[forte]:
                return (
                    f"il segnale dell'orientamento atteso e' nel {nomi[forte]} "
                    f"({parti[forte][0]} letture su {parti[forte][1]}, quota "
                    f"{float(quote[forte]):.4f}) e quasi assente nel {nomi[debole]} "
                    f"({parti[debole][0]} su {parti[debole][1]}, quota "
                    f"{float(quote[debole]):.4f}): {conclusione}"
                ).replace("nel letture", "nelle letture")
        return None

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


def _resto_del_file(file: IO[str]) -> str | None:
    """Che cosa c'e' fino alla fine del file dopo una riga vuota: ``None`` se
    solo righe vuote, altrimenti la descrizione del difetto.

    Una riga di soli spazi o tabulazioni non e' una riga vuota: i lettori di R
    la leggono come l'inizio di un record.
    """
    for riga in file:
        if riga.strip():
            return "seguono altre righe"
        if riga.strip("\r\n"):
            return "segue una riga di soli spazi"
    return None


def scansiona_file(
    percorso: Path | str,
    head_reads: int,
    primer: str | None,
    motivo: str | None,
    inizio_motivo: int = 0,
    lunghezza_richiesta: int = 0,
    orientamento: SegnaleDiOrientamento | None = None,
    altro_orientamento: SegnaleDiOrientamento | None = None,
    ricerca: RicercaPrimer | None = None,
    sonda: RicercaPrimer | None = None,
) -> StatisticheFile:
    """Legge le prime ``head_reads`` letture e ne ricava tutte le statistiche.

    Con ``sonda`` si registra soltanto a quale sfasamento il primer si trova
    (:attr:`StatisticheFile.sfasamenti_del_primer`): nessun'altra statistica
    cambia.

    Con ``ricerca`` il primer si cerca lungo la lettura (``primer`` e
    ``inizio_motivo`` non si usano): si contano le letture in cui si trova, il
    motivo si cerca dopo il primer trovato (dalla prima base se non c'e'), e le
    lunghezze sono quelle della lettura dopo il taglio; una lettura senza
    primer che il filtro scartera' non entra nelle lunghezze.

    ``lunghezza_richiesta`` e' la lunghezza sotto la quale una lettura si conta
    fra le piu' corte (``filter.truncLen``); zero non ne conta nessuna.

    ``primer`` e ``motivo`` sono espressioni regolari. Il primer si cerca
    **ancorato a inizio lettura**: starebbe in testa se non fosse stato tolto.
    Il motivo conservato, che apre la regione amplificata, si cerca ancorato
    alla posizione ``inizio_motivo``: zero se le letture iniziano dalla
    regione, la lunghezza di cio' che la precede (``filter.trimLeft``) se
    portano ancora il primer che il filtro togliera'. Con ``None`` la ricerca
    corrispondente non si fa e il conteggio resta zero.

    Con ``orientamento`` si registra per ogni lettura se porta il segnale
    dell'orientamento atteso, e le statistiche ne riportano il conteggio per
    tratto del file e per posto dispari o pari
    (:attr:`StatisticheFile.segnale_per_tratto`,
    :attr:`StatisticheFile.segnale_alterno`): un byte per lettura, per la
    durata della scansione. Con ``altro_orientamento`` si contano anche le
    letture che portano il segnale dell'altro verso
    (:attr:`StatisticheFile.con_altro_segnale`): serve, in modalita' appaiata,
    a riconoscere una corsa invertita e una libreria con i due versi mescolati.
    """
    percorso = Path(percorso)
    espressione_primer = re.compile(primer) if primer is not None else None
    espressione_motivo = re.compile(motivo) if motivo is not None else None
    in_testa = dopo_il_taglio = None
    if orientamento is not None:
        in_testa = re.compile(orientamento.primer) if orientamento.primer is not None else None
        dopo_il_taglio = (
            re.compile(orientamento.motivo) if orientamento.motivo is not None else None)
    con_segnale = bytearray()
    altro_in_testa = altro_dopo_il_taglio = None
    if altro_orientamento is not None:
        altro_in_testa = (
            re.compile(altro_orientamento.primer) if altro_orientamento.primer is not None
            else None)
        altro_dopo_il_taglio = (
            re.compile(altro_orientamento.motivo) if altro_orientamento.motivo is not None
            else None)
    con_altro = 0
    sfasamenti: Counter[int] = Counter()

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
                    # di quattro righe li sfaserebbe tutti. Una riga di soli
                    # spazi non e' spaziatura: il lettore di R la prende per un
                    # record, e va respinta qui.
                    difetto = ("e' di soli spazi" if intestazione.strip("\r\n")
                               else _resto_del_file(file))
                    if difetto is None:
                        vuote_finali = True
                        break
                    return StatisticheFile(
                        percorso.name, letture, minima, massima, con_primer, con_motivo,
                        errore=(f"riga di soli spazi dopo {letture} letture: non e' una riga "
                                "vuota, e i lettori la prendono per un record"
                                if "soli spazi" in difetto else
                                f"riga vuota in mezzo al file dopo {letture} letture: "
                                f"seguono altre righe"),
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
                if sonda is not None:
                    trovata = sonda.trova(sequenza)
                    if trovata is not None:
                        sfasamenti[trovata[0]] += 1
                lunghezza = len(sequenza)
                misurata = True
                if ricerca is not None:
                    fine = ricerca.fine(sequenza)
                    if fine is not None:
                        con_primer += 1
                        lunghezza -= fine
                    else:
                        misurata = not ricerca.scarta_senza_primer
                    inizio_motivo = fine or 0
                elif espressione_primer is not None and espressione_primer.match(sequenza):
                    con_primer += 1
                if misurata:
                    minima = lunghezza if minima is None else min(minima, lunghezza)
                    massima = lunghezza if massima is None else max(massima, lunghezza)
                    if lunghezza < lunghezza_richiesta:
                        corte += 1
                # Sulla lettura tagliata, non dalla posizione del taglio: un
                # motivo ancorato con ^ deve valere all'inizio di cio' che resta.
                if espressione_motivo is not None and espressione_motivo.match(
                    sequenza[inizio_motivo:]
                ):
                    con_motivo += 1
                if orientamento is not None:
                    con_segnale.append(orientamento.lo_porta(sequenza, in_testa, dopo_il_taglio))
                if altro_orientamento is not None and altro_orientamento.lo_porta(
                        sequenza, altro_in_testa, altro_dopo_il_taglio):
                    con_altro += 1
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
            righe_vuote_finali=vuote_finali, vuoto=True,
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
        segnale_per_tratto=_per_tratto(con_segnale) if orientamento is not None else None,
        segnale_alterno=(
            ((sum(con_segnale[0::2]), len(con_segnale[0::2])),
             (sum(con_segnale[1::2]), len(con_segnale[1::2])))
            if orientamento is not None else None),
        con_altro_segnale=con_altro if altro_orientamento is not None else None,
        sfasamenti_del_primer=tuple(sorted(sfasamenti.items())) if sonda is not None else None,
    )


def _per_tratto(con_segnale: bytearray) -> tuple[tuple[int, int], ...]:
    """Le letture con il segnale in ciascuno dei :data:`TRATTI` tratti
    consecutivi di uguale lunghezza del file, con le letture del tratto.
    """
    confini = [i * len(con_segnale) // TRATTI for i in range(TRATTI + 1)]
    return tuple((sum(con_segnale[a:b]), b - a) for a, b in zip(confini, confini[1:]))


#: La cartella delle copie ripulite, dentro la cartella della fase che le usa.
#: Il nome ha il prefisso delle scritture in corso (``io_layer/artifacts.py``):
#: se la fase e' interrotta, l'avvio successivo la toglie.
CARTELLA_DELLE_COPIE: Final = ".scrittura-letture"


@contextmanager
def senza_righe_vuote_finali(
    percorsi: dict[str, str], con_righe_vuote: Iterable[str], in_cartella: Path | str | None = None,
    nome_cartella: str = CARTELLA_DELLE_COPIE,
) -> Iterator[dict[str, str]]:
    """I percorsi da dare ai lettori di R, con una copia ripulita dei file indicati.

    Le righe vuote in fondo a un FASTQ sono tollerate, ma i lettori a blocchi
    di R (``ShortRead``, quindi anche ``dada2::filterAndTrim``) le rifiutano
    come un record che non comincia con ``@``. Dei soli file in
    ``con_righe_vuote`` (chiavi di ``percorsi``) si scrive una copia senza
    quelle righe, che vive quanto il blocco ``with``: i file di ingresso non
    si toccano mai, e senza file da ripulire non si scrive nulla.

    Le copie stanno in :data:`CARTELLA_DELLE_COPIE` sotto ``in_cartella`` (la
    cartella della fase): occupano fino a una volta i file ripuliti, e lo
    spazio e' quello della cartella di uscita, che la validazione ha
    verificato, non quello della cartella temporanea del sistema. Senza
    ``in_cartella`` si usa quest'ultima. ``nome_cartella`` distingue due insiemi
    di copie aperti insieme nella stessa cartella (le due letture di una
    coppia): deve cominciare con il prefisso delle scritture in corso.
    """
    chiavi = sorted(con_righe_vuote)
    if not chiavi:
        yield percorsi
        return
    if in_cartella is None:
        cartella = Path(tempfile.mkdtemp(prefix="amplicon16s-letture-"))
    else:
        cartella = Path(in_cartella) / nome_cartella
        shutil.rmtree(cartella, ignore_errors=True)
        cartella.mkdir(parents=True)
    try:
        ripuliti = dict(percorsi)
        for numero, chiave in enumerate(chiavi):
            # Una sottocartella per file: due file di cartelle diverse possono
            # avere lo stesso nome.
            destinazione = cartella / str(numero) / Path(percorsi[chiave]).name
            destinazione.parent.mkdir()
            _copia_senza_righe_vuote_finali(percorsi[chiave], destinazione)
            ripuliti[chiave] = str(destinazione)
        yield ripuliti
    finally:
        shutil.rmtree(cartella, ignore_errors=True)


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


def conta_coppie(
    percorso: Path | str, orientamento: SegnaleDiOrientamento | None = None,
    altro_orientamento: SegnaleDiOrientamento | None = None,
) -> StatisticheFile:
    """I segni di coppia su **tutte** le letture del file, con le regole di G07.

    E' la stessa scansione dei gate senza il limite delle prime letture e senza
    le ricerche di primer e motivo: il giudizio
    (:attr:`StatisticheFile.coppie_nello_stesso_file`) e' lo stesso, e cio' che
    G07 non vede oltre ``qc.head_reads`` qui si vede. I nomi di lettura di un
    file restano in memoria solo per la durata della sua scansione. Dice anche
    se il file si chiude con righe vuote (``righe_vuote_finali``) e, con
    ``orientamento``, quante letture portano il segnale dell'orientamento
    atteso in ogni tratto del file (``segnale_per_tratto``).
    """
    return scansiona_file(percorso, sys.maxsize, None, None, orientamento=orientamento,
                          altro_orientamento=altro_orientamento)


def scansiona(
    percorsi: dict[str, Path],
    head_reads: int,
    primer: str | None,
    motivo: str | None,
    inizio_motivo: int = 0,
    lunghezza_richiesta: int = 0,
    ricerca: RicercaPrimer | None = None,
    sonda: RicercaPrimer | None = None,
) -> dict[str, StatisticheFile]:
    """Scansiona i file indicati e restituisce le statistiche di ciascuno."""
    return {
        chiave: scansiona_file(percorso, head_reads, primer, motivo, inizio_motivo,
                               lunghezza_richiesta, ricerca=ricerca, sonda=sonda)
        for chiave, percorso in percorsi.items()
    }


# --------------------------------------------------------------------------- #
# Modalita' appaiata: i due file di ogni campione                              #
# --------------------------------------------------------------------------- #

#: Le coppie di marcatori con cui i depositi e i sequenziatori distinguono nel
#: nome i due file di una corsa, dalla piu' specifica: ``R1`` e ``R2``, oppure
#: il solo numero, dopo un trattino basso o un punto. Sono convenzioni diffuse,
#: non una proprieta' di un dataset.
MARCATORI_DI_FILE: Final[tuple[tuple[str, str], ...]] = (
    ("_R1", "_R2"), (".R1", ".R2"), ("_1", "_2"), (".1", ".2"),
)


#: Il secondo file di una coppia, riconosciuto dal marcatore che precede
#: l'estensione: ``_R2``, ``.R2``, ``_2`` o ``.2``, con un eventuale numero di
#: blocco (``_001``) o il suffisso ``_raw``. Ancorato alla fine del nome: un campione che si chiama
#: ``LAB_R2D2`` non e' una lettura inversa. Sono convenzioni diffuse, non una
#: proprieta' di un dataset.
NOME_DI_LETTURA_INVERSA: Final = re.compile(
    r"(?P<prima>.*)(?P<segno>[._]R2|[._]2)(?P<dopo>(?:_[0-9]+|_raw)?\.(?:fastq|fq)(?:\.gz)?)$",
    re.IGNORECASE,
)


def nome_da_seconda_lettura(nome: str, nomi: Iterable[str]) -> bool:
    """Vero se il nome del file indica la seconda lettura di una coppia.

    ``R2`` basta da solo. Il solo numero (``_2``, ``.2``) e' anche il modo in
    cui si numerano campioni e repliche di un dataset single-end: conta come
    seconda lettura solo se fra ``nomi`` (i nomi dei file della cartella, senza
    distinzione fra maiuscole e minuscole) c'e' il compagno con ``1`` allo
    stesso posto.

    E' la sola definizione di nome sospetto: la usano la validazione, per
    respingere il file dove il contenuto non si puo' giudicare, e il controllo
    di orientamento, per giudicarlo da meno letture
    (:data:`LETTURE_MINIME_PER_UN_NOME_SOSPETTO`).
    """
    trovato = NOME_DI_LETTURA_INVERSA.match(nome)
    if trovato is None:
        return False
    segno = trovato.group("segno")
    if segno[1:].upper() == "R2":
        return True
    compagno = f"{trovato.group('prima')}{segno[0]}1{trovato.group('dopo')}"
    return compagno.casefold() in {n.casefold() for n in nomi}


#: Cio' che puo' seguire il marcatore di coppia nel nome di un file: un numero
#: di blocco di tre cifre (``_001``, la forma dei file Illumina) oppure il
#: suffisso ``_raw`` (la forma dei file grezzi di alcuni archivi), e
#: l'estensione. Un numero piu' corto e' parte del nome: ``X_1_2`` e' il
#: reverse del campione ``X_1``, non il forward ``X`` del blocco 2.
_DOPO_IL_MARCATORE: Final = r"(?P<dopo>(?:_[0-9]{3}|_raw)?\.(?:fastq|fq)(?:\.gz)?)$"


def _con_l_altro_marcatore(nome: str, da: str, a: str) -> str | None:
    """Il nome con il marcatore ``da`` sostituito da ``a``, o ``None`` se il
    nome non finisce con ``da`` (seguito al piu' dal numero di blocco o da
    ``_raw``, e dall'estensione).

    Il marcatore conta solo in fondo al nome: un campione che lo contiene nel
    proprio nome (``X_1`` con i marcatori ``_1`` e ``_2``) non e' per questo
    il forward di un altro file.
    """
    trovato = re.match(rf"(?P<prima>.*){re.escape(da)}{_DOPO_IL_MARCATORE}", nome, re.IGNORECASE)
    if trovato is None:
        return None
    return f"{trovato.group('prima')}{a}{trovato.group('dopo')}"


@dataclass(frozen=True)
class Abbinamento:
    """I file di una cartella abbinati per coppia di marcatori."""

    #: File forward -> il suo file reverse.
    coppie: dict[str, str]
    #: I file che non sono ne' il forward ne' il reverse di un altro file.
    senza_compagno: tuple[str, ...]
    #: I file che risultano insieme forward e reverse: il nome porta i due
    #: marcatori in modo che la regola non li distingua.
    ambigui: tuple[str, ...]

    @property
    def completo(self) -> bool:
        """Vero se ogni file ha il suo compagno, senza ambiguita'."""
        return bool(self.coppie) and not self.senza_compagno and not self.ambigui


def abbina_file(nomi: Iterable[str], marcatori: tuple[str, str] | list[str]) -> Abbinamento:
    """Abbina i file per nome: un file e' forward se il suo nome finisce con
    il primo marcatore (prima dell'estensione, del numero di blocco o di
    ``_raw``) e lo
    stesso nome con il secondo marcatore al suo posto e' un altro file
    dell'insieme; reverse nel caso simmetrico.

    Il marcatore conta solo in fondo al nome: e' la sua ultima occorrenza, e
    un campione che lo contiene nel proprio nome non si confonde.
    """
    del_forward, del_reverse = marcatori
    presenti = set(nomi)
    coppie: dict[str, str] = {}
    reverse: set[str] = set()
    for nome in sorted(presenti):
        compagno = _con_l_altro_marcatore(nome, del_forward, del_reverse)
        if compagno is not None and compagno != nome and compagno in presenti:
            coppie[nome] = compagno
            reverse.add(compagno)
    ambigui = sorted(set(coppie) & reverse)
    for nome in ambigui:
        del coppie[nome]
    abbinati = set(coppie) | reverse
    return Abbinamento(
        coppie={f: r for f, r in coppie.items() if r not in ambigui},
        senza_compagno=tuple(sorted(presenti - abbinati)),
        ambigui=tuple(ambigui),
    )


def compagni_per_nome(nomi: Iterable[str]) -> tuple[str, str] | None:
    """I marcatori di :data:`MARCATORI_DI_FILE` per cui i nomi dati (i file di
    uno stesso campione) sono esattamente una coppia forward e reverse, o
    ``None``.
    """
    nomi = sorted(nomi)
    if len(nomi) != 2:
        return None
    for marcatori in MARCATORI_DI_FILE:
        abbinamento = abbina_file(nomi, marcatori)
        if abbinamento.completo and len(abbinamento.coppie) == 1:
            return marcatori
    return None


def _identificativo(intestazione: str) -> tuple[str, str | None]:
    """L'identificativo di una lettura (il primo campo dell'intestazione) senza
    il marcatore di coppia, e il marcatore: ``/1`` e ``/2`` in fondo, oppure 1
    e 2, R1 e R2 dopo un separatore.

    Solo l'identificativo: il commento che lo segue puo' differire fra le due
    letture di una coppia (la lunghezza, il numero della lettura).
    """
    campi = intestazione[1:].split(None, 1)
    nome = campi[0] if campi else ""
    trovato = _SUFFISSO_COPPIA.search(nome)
    if trovato:
        return nome[:trovato.start()], trovato.group(1)
    trovato = _MARCATORE_IN_FONDO.match(nome)
    if trovato:
        return trovato.group(1), trovato.group(2)
    return nome, None


#: L'identificativo che un archivio assegna numerando le letture: un prefisso,
#: un punto e il numero progressivo della lettura.
_NUMERATA: Final = re.compile(r"^(.+)\.([0-9]+)$")


def _relazione(prima: str, seconda: str) -> str | int | None:
    """In che modo due intestazioni sono le due letture della stessa coppia, o
    ``None`` se non lo sono.

    * ``"uguali"``: lo stesso identificativo;
    * ``"marcatori"``: lo stesso identificativo a meno del marcatore di coppia
      (``/1`` e ``/2``, oppure 1 e 2, R1 e R2 dopo un separatore);
    * un intero: gli identificativi sono numeri progressivi dello stesso
      archivio (``prefisso.numero``) e questo e' il loro scarto. Un archivio
      che deposita le due letture in un solo file le numera di seguito: la
      seconda lettura della coppia k e' la lettura N + k, e lo scarto N e' lo
      stesso per ogni coppia finche' i due file sono in sincronia. Che lo
      scarto non sia minore delle letture del file lo verifica chi arriva in
      fondo (:func:`verifica_sincronia`).
    """
    campi_1, campi_2 = prima[1:].split(None, 1), seconda[1:].split(None, 1)
    nome_1 = campi_1[0] if campi_1 else ""
    nome_2 = campi_2[0] if campi_2 else ""
    if nome_1 == nome_2:
        return "uguali"
    (senza_1, marcatore_1), (senza_2, marcatore_2) = _identificativo(prima), _identificativo(seconda)
    if (senza_1 == senza_2 and marcatore_1 is not None and marcatore_2 is not None
            and ((marcatore_1, marcatore_2) in _PRIMA_E_SECONDA
                 # Un file gia' scambiato porta i marcatori nell'altro ordine: la
                 # sincronia c'e', e l'inversione la giudica chi guarda i marcatori.
                 or (marcatore_2, marcatore_1) in _PRIMA_E_SECONDA)):
        return "marcatori"
    numerata_1, numerata_2 = _NUMERATA.match(nome_1), _NUMERATA.match(nome_2)
    if numerata_1 and numerata_2 and numerata_1.group(1) == numerata_2.group(1):
        return int(numerata_2.group(2)) - int(numerata_1.group(2))
    return None


@dataclass(frozen=True)
class Sincronia:
    """L'esito del confronto, lettura per lettura, dei due file di un campione."""

    #: Le coppie lette prima della fine o del primo difetto.
    coppie: int
    #: Perche' i due file non sono in sincronia, o ``None``.
    difetto: str | None = None


def verifica_sincronia(forward: Path | str, reverse: Path | str) -> Sincronia:
    """Confronta i nomi delle letture dei due file di un campione, una coppia
    per volta: la lettura n dell'uno e la lettura n dell'altro devono essere la
    stessa coppia, e i due file devono finire insieme.

    Che cosa vuol dire "la stessa coppia" lo dice la prima (:func:`_relazione`):
    nomi uguali, nomi che differiscono per il solo marcatore di coppia, oppure
    numeri progressivi di un archivio a distanza costante. Ogni coppia
    successiva deve stare nella stessa relazione: una lettura in piu' o in meno
    in uno dei due file la rompe dal punto in cui manca.

    Si leggono le sole intestazioni, in flusso: la memoria non cresce con i
    file. Il confronto si ferma al primo difetto e ne dice il posto. Un guasto
    di lettura non e' un difetto di sincronia: lo dichiara la scansione di
    ciascun file, e qui il confronto si ferma dove arriva.
    """
    coppie = 0
    attesa: str | int | None = None
    try:
        with apri_fastq(forward) as primo, apri_fastq(reverse) as secondo:
            while True:
                intestazione_1, intestazione_2 = primo.readline(), secondo.readline()
                # Le righe vuote in fondo a un file sono tollerate.
                fine_1 = not intestazione_1.strip()
                fine_2 = not intestazione_2.strip()
                if fine_1 or fine_2:
                    if fine_1 and fine_2:
                        # Numeri progressivi: le seconde letture sono numerate
                        # dopo tutte le prime, quindi lo scarto non puo' essere
                        # minore delle letture del file. Uno scarto minore e'
                        # quello di due file sfasati di altrettante letture.
                        if isinstance(attesa, int) and abs(attesa) < coppie:
                            return Sincronia(coppie, (
                                f"i nomi delle due letture sono numeri progressivi a distanza "
                                f"{attesa}, minore delle {coppie} letture di ciascun file: in "
                                "un archivio che numera le seconde letture dopo le prime la "
                                "distanza e' almeno il numero delle letture, e una distanza "
                                "minore dice che i due file sono sfasati di altrettante letture"))
                        return Sincronia(coppie)
                    corto = "forward" if fine_1 else "reverse"
                    return Sincronia(coppie, (
                        f"il file {corto} finisce dopo {coppie} letture e l'altro ne ha "
                        "ancora: i due file non hanno lo stesso numero di letture"))
                relazione = _relazione(intestazione_1, intestazione_2)
                if coppie == 0:
                    # La prima coppia dice come i due file nominano le letture;
                    # le successive devono farlo allo stesso modo.
                    attesa = relazione
                if relazione is None or relazione != attesa:
                    return Sincronia(coppie, (
                        f"alla lettura {coppie + 1} i nomi non corrispondono: "
                        f"{intestazione_1.strip()[:60]!r} nel file forward e "
                        f"{intestazione_2.strip()[:60]!r} nel file reverse"
                        + ("" if relazione is None or coppie == 0 else
                           " (fino a qui i nomi delle due letture "
                           + ("coincidevano" if attesa == "uguali" else
                              "si distinguevano per il solo marcatore di coppia"
                              if attesa == "marcatori" else
                              f"erano numeri progressivi a distanza {attesa}") + ")")))
                for _ in range(3):
                    primo.readline()
                    secondo.readline()
                coppie += 1
    except (OSError, EOFError, zlib.error):
        return Sincronia(coppie)


def corsa_invertita(forward: StatisticheFile, reverse: StatisticheFile) -> str | None:
    """Perche' i due file di una corsa risultano scambiati, dal contenuto, o
    ``None``.

    Per ogni segnale cercato, il file che dovrebbe portarlo ne e' privo (almeno
    :data:`LETTURE_MINIME_PER_LA_QUOTA` letture, quota sotto
    :data:`QUOTA_MINIMA_DEL_SEGNALE`) e il compagno lo porta (almeno
    :data:`SEGNALE_MINIMO_DI_UNA_PARTE` letture con quel segnale). La prova e'
    positiva, quindi vale per ogni classe di campione. Serve almeno un segnale
    cercato in entrambi i file.
    """
    prove = []
    # Il segnale forward: atteso nel file forward, cercato come "altro" nel reverse.
    for nome, atteso_in, cercato_in in (
        ("forward", forward, reverse), ("reverse", reverse, forward),
    ):
        propri, altrui = atteso_in.con_segnale, cercato_in.con_altro_segnale
        if propri is None or altrui is None:
            continue
        letture = atteso_in.letture_esaminate
        privo = (letture >= LETTURE_MINIME_PER_LA_QUOTA
                 and propri < QUOTA_MINIMA_DEL_SEGNALE * letture)
        if not (privo and altrui >= SEGNALE_MINIMO_DI_UNA_PARTE):
            return None
        dove_manca = "forward" if nome == "forward" else "reverse"
        dove_sta = "reverse" if nome == "forward" else "forward"
        prove.append(
            f"il segnale delle letture {nome} e' in {propri} letture su {letture} del file "
            f"{dove_manca} e in {altrui} su {cercato_in.letture_esaminate} del file {dove_sta}")
    return "; ".join(prove) if prove else None


def versi_mescolati(statistiche: StatisticheFile) -> str | None:
    """Perche' un file riporta letture dei due versi mescolate, o ``None``.

    Almeno :data:`SEGNALE_MINIMO_DI_UNA_PARTE` letture portano il segnale
    atteso e almeno altrettante quello dell'altro verso, e la quota minore non
    e' sotto :data:`RAPPORTO_MASSIMO_FRA_LE_PARTI` della maggiore. Servono
    entrambi i segnali.
    """
    propri, altrui = statistiche.con_segnale, statistiche.con_altro_segnale
    if propri is None or altrui is None:
        return None
    minore, maggiore = sorted((propri, altrui))
    if (minore >= SEGNALE_MINIMO_DI_UNA_PARTE
            and minore >= RAPPORTO_MASSIMO_FRA_LE_PARTI * maggiore):
        letture = statistiche.letture_esaminate
        return (
            f"{propri} letture su {letture} portano il segnale del verso atteso (quota "
            f"{propri / letture:.4f}) e {altrui} quello dell'altro (quota {altrui / letture:.4f})")
    return None
