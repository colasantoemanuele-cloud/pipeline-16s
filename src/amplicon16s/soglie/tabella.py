"""Le soglie il cui solo effetto e' fermare l'esecuzione o dare un avviso.

Nessuna di queste soglie cambia conteggi, tassonomia, campioni o varianti:
ognuna si confronta con una misura dell'esecuzione, e il controllo scatta se la
misura la supera (soglia massima) o le sta sotto (soglia minima). Tutti i
confronti sono stretti: una misura uguale alla soglia non fa scattare nulla.

Per ogni soglia la tabella dichiara cio' che il codice del controllo e lo
schema della configurazione stabiliscono:

* la **direzione** e l'unita' (frazione, conteggio, basi);
* il **limite**: il valore oltre il quale il controllo non puo' piu' scattare.
  Una frazione non supera 1 e non scende sotto 0; i campioni azzerati non sono
  piu' dei campioni controllati; per lo scarto del troncamento lo schema ammette
  al piu' ``filter.truncLen - 1``, che vale come limite; le varianti distinte
  non hanno un limite che la configurazione conosca;
* quando il controllo **non si esegue**, qualunque sia la soglia (il primer non
  si cerca con ``filter.trimLeft`` maggiore di zero, il motivo con
  ``qc.conserved_motif`` nullo, la frazione di contaminanti senza controlli
  negativi dichiarati, omessa o dichiarata che sia);
* il **margine** con cui ``soglie`` propone un valore dalla misura osservata.

Una soglia uguale al limite, o oltre nella sua direzione, e' un **controllo
spento**: la configurazione e' valida, ma quel controllo non protegge. Chi
esegue deve saperlo: :func:`controlli_spenti` li elenca, per la riga di comando
e per il report.

Razionale sistemistico: qui non c'e' alcun valore di un dataset. I limiti
vengono dagli intervalli dello schema e dalla forma dei confronti; i margini
sono costanti con un nome e un motivo.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from decimal import Decimal
from enum import StrEnum
from typing import Any, Final

from amplicon16s.config import defaults
from amplicon16s.config.schema import Config

__all__ = [
    "MARGINE_ARRESTO_CHIMERE",
    "MARGINE_AZZERATI",
    "MARGINE_BASI",
    "MARGINE_FRAZIONE",
    "MARGINE_VARIANTI",
    "SOGLIA_SENZA_LIMITE_AMMESSO",
    "SOGLIE",
    "ControlloSpento",
    "Direzione",
    "Soglia",
    "Unita",
    "controlli_spenti",
    "limite",
    "non_eseguito",
    "oltre_il_limite",
    "righe_di_avviso",
    "soglia",
    "valore_esplorativo",
    "valori_piatti",
]

#: Margine delle soglie in frazione: cinque punti percentuali. Assorbe lo
#: spostamento di una misura d'insieme quando si aggiungono campioni dello
#: stesso studio o si ripete l'esecuzione dopo un aggiornamento che non cambia
#: il metodo, senza lasciar passare un peggioramento netto.
MARGINE_FRAZIONE: Final = Decimal("0.05")

#: Margine dell'arresto sulle chimere: il doppio di quello dell'avviso, cosi'
#: che fra avviso e arresto resti una fascia in cui l'esecuzione prosegue e
#: segnala.
MARGINE_ARRESTO_CHIMERE: Final = Decimal("0.10")

#: Margine sulle varianti distinte, in proporzione alla misura: un conteggio
#: senza scala propria, che cresce con campioni e profondita'. La meta' in piu'
#: regge l'aggiunta di campioni; un primer non tolto o letture di bassa qualita'
#: moltiplicano le varianti di piu' volte.
MARGINE_VARIANTI: Final = Decimal("0.5")

#: Margine sui campioni azzerati: nessuno. Ogni campione azzerato in piu' e' un
#: campione perso, un evento discreto e non una fluttuazione.
MARGINE_AZZERATI: Final = 0

#: Margine sullo scarto del troncamento, in basi: nessuno. Lo scarto cresce
#: solo se il troncamento si abbassa o esce il campione con la lettura piu'
#: corta; e' la stessa regola del valore di partenza di ``init``.
MARGINE_BASI: Final = 0


class Direzione(StrEnum):
    """Da che parte la misura fa scattare il controllo."""

    MASSIMA = "massima"
    MINIMA = "minima"


class Unita(StrEnum):
    """In che cosa si esprime la soglia."""

    FRAZIONE = "frazione"
    CONTEGGIO = "conteggio"
    BASI = "basi"


@dataclass(frozen=True)
class Soglia:
    """Una soglia che ferma o avvisa."""

    chiave: str
    direzione: Direzione
    unita: Unita
    #: I codici del catalogo con cui il controllo scatta.
    codici: tuple[str, ...]
    #: Vero se il controllo ferma l'esecuzione, falso se da' un avviso.
    ferma: bool
    #: Il margine del suggerimento: assoluto, oppure una quota della misura.
    margine: Decimal | int
    margine_proporzionale: bool = False

    @property
    def testo_margine(self) -> str:
        """Il margine a parole, per chi legge la regola applicata."""
        if self.margine_proporzionale:
            return f"{self.margine} della misura"
        if self.unita is Unita.FRAZIONE:
            return str(self.margine)
        return f"{self.margine} {'basi' if self.unita is Unita.BASI else 'campioni'}"


_MASSIMA, _MINIMA = Direzione.MASSIMA, Direzione.MINIMA
_FRAZIONE = Unita.FRAZIONE

#: Le dodici soglie, nell'ordine delle fasi che le applicano.
SOGLIE: Final[tuple[Soglia, ...]] = (
    Soglia("qc.max_primer_hit_frac", _MASSIMA, _FRAZIONE, ("E-S0-10",), True, MARGINE_FRAZIONE),
    Soglia("qc.min_motif_frac", _MINIMA, _FRAZIONE, ("E-S0-16",), True, MARGINE_FRAZIONE),
    Soglia("qc.max_frac_short_reads", _MASSIMA, _FRAZIONE, ("E-S0-09", "E-S1-02"), True,
           MARGINE_FRAZIONE),
    Soglia("filter.truncLen_shortfall_warn", _MASSIMA, Unita.BASI, ("E-S1-01",), False,
           MARGINE_BASI),
    Soglia("qc.max_zeroed_samples", _MASSIMA, Unita.CONTEGGIO, ("E-S2-01",), True,
           MARGINE_AZZERATI),
    Soglia("qc.max_frac_lost_filter", _MASSIMA, _FRAZIONE, ("E-S2-02",), True, MARGINE_FRAZIONE),
    Soglia("qc.max_asv_count", _MASSIMA, Unita.CONTEGGIO, ("E-S5-01",), True, MARGINE_VARIANTI,
           margine_proporzionale=True),
    Soglia("qc.warn_frac_chimeric", _MASSIMA, _FRAZIONE, ("E-S6-02",), False, MARGINE_FRAZIONE),
    Soglia("qc.stop_frac_chimeric", _MASSIMA, _FRAZIONE, ("E-S6-01",), True,
           MARGINE_ARRESTO_CHIMERE),
    Soglia("qc.min_frac_phylum", _MINIMA, _FRAZIONE, ("E-S8-02",), True, MARGINE_FRAZIONE),
    Soglia("qc.max_frac_contaminant", _MASSIMA, _FRAZIONE, ("E-S12-02",), True, MARGINE_FRAZIONE),
    Soglia("qc.min_frac_reads_retained", _MINIMA, _FRAZIONE, ("E-S14-01",), True,
           MARGINE_FRAZIONE),
)

_PER_CHIAVE: Final = {s.chiave: s for s in SOGLIE}

#: La soglia per cui nessun valore ammesso dallo schema spegne di sicuro il
#: controllo: lo scarto non ha un massimo, e lo schema respinge un valore pari o
#: superiore a ``filter.truncLen``. Il massimo ammesso vale come limite.
SOGLIA_SENZA_LIMITE_AMMESSO: Final = "filter.truncLen_shortfall_warn"


def soglia(chiave: str) -> Soglia:
    """La voce della tabella per un parametro."""
    return _PER_CHIAVE[chiave]


def valori_piatti(config: Config) -> dict[str, Any]:
    """I valori della configurazione nella forma ``gruppo.parametro``."""
    return {
        f"{gruppo}.{nome}": valore
        for gruppo, valori in config.model_dump(mode="python").items()
        for nome, valore in valori.items()
    }


def non_eseguito(voce: Soglia, valori: Mapping[str, Any]) -> str | None:
    """Perche' il controllo della soglia non si esegue con questa
    configurazione, qualunque valore abbia la soglia; ``None`` se si esegue.
    """
    if valori.get(voce.chiave) is None:
        return "il parametro e' omesso: il suo controllo non riguarda questo dataset"
    if voce.chiave in defaults.CONDIZIONALI and not defaults.condizione_vera(voce.chiave, valori):
        # Dichiarato comunque, ma la condizione che lo rende pertinente e' falsa.
        return (f"serve solo se {defaults.testo_condizione(voce.chiave)}: dichiarato senza "
                "che la condizione sia vera, non ha effetto")
    if voce.chiave == "qc.max_primer_hit_frac" and (valori.get("filter.trimLeft") or 0) > 0:
        return ("con filter.trimLeft maggiore di zero la validazione non cerca il primer in "
                "testa alle letture: la soglia non ha effetto")
    if voce.chiave == "qc.min_motif_frac" and valori.get("qc.conserved_motif") is None:
        return ("qc.conserved_motif e' nullo: la validazione non cerca il motivo, e la soglia "
                "non ha effetto")
    return None


def limite(
    voce: Soglia, valori: Mapping[str, Any], campioni_controllati: int | None = None
) -> int | None:
    """Il valore oltre il quale il controllo non puo' piu' scattare, o
    ``None`` se non e' definito o non si puo' ancora stabilire.

    ``campioni_controllati`` e' il numero dei campioni biologici e dei
    controlli positivi dell'inventario: serve alla sola soglia sui campioni
    azzerati, e senza di esso quel limite non si stabilisce.
    """
    if voce.unita is Unita.FRAZIONE:
        return 1 if voce.direzione is Direzione.MASSIMA else 0
    if voce.chiave == "qc.max_zeroed_samples":
        return campioni_controllati
    if voce.chiave == "filter.truncLen_shortfall_warn":
        troncamento = valori.get("filter.truncLen")
        return None if troncamento is None else int(troncamento) - 1
    return None


def oltre_il_limite(voce: Soglia, valore: Any, al_limite: int | None) -> bool:
    """Vero se una soglia con questo valore non puo' far scattare il controllo."""
    if al_limite is None or valore is None:
        return False
    if voce.direzione is Direzione.MASSIMA:
        return valore >= al_limite
    return valore <= al_limite


def valore_esplorativo(voce: Soglia) -> float | None:
    """Il valore che non ferma, per una prima esecuzione esplorativa.

    Lo hanno le sole soglie in frazione, il cui limite e' un estremo fisso
    dell'intervallo ammesso; le altre hanno un valore di partenza che ``init``
    ricava dai dati o dalla macchina, o restano da scegliere.
    """
    if voce.unita is not Unita.FRAZIONE:
        return None
    return 1.0 if voce.direzione is Direzione.MASSIMA else 0.0


@dataclass(frozen=True)
class ControlloSpento:
    """Una soglia al limite o oltre: il suo controllo non puo' scattare."""

    chiave: str
    valore: Any
    limite: int
    codici: tuple[str, ...]
    ferma: bool

    @property
    def testo(self) -> str:
        """Una riga per l'operatore."""
        verso = "pari o superiore a" if soglia(self.chiave).direzione is Direzione.MASSIMA else (
            "pari o inferiore a")
        che_cosa = "l'arresto" if self.ferma else "l'avviso"
        if self.chiave == SOGLIA_SENZA_LIMITE_AMMESSO:
            # Lo schema non ammette un valore piu' largo: il controllo e' al
            # suo minimo, non impossibile.
            return (f"{self.chiave}: {self.valore} (il massimo che lo schema ammette, "
                    f"filter.truncLen meno 1: {che_cosa} {', '.join(self.codici)} scatta solo "
                    f"se la lettura piu' corta supera il troncamento di oltre {self.limite} basi)")
        return (f"{self.chiave}: {self.valore} ({verso} {self.limite}: {che_cosa} "
                f"{', '.join(self.codici)} non puo' scattare)")


def controlli_spenti(
    valori: Mapping[str, Any], campioni_controllati: int | None = None
) -> list[ControlloSpento]:
    """Le soglie della configurazione che non possono far scattare il loro
    controllo, nell'ordine della tabella.

    Non comprende le soglie il cui controllo non si esegue comunque
    (:func:`non_eseguito`): non sono spente dal loro valore, e nessun valore le
    accenderebbe.
    """
    spenti = []
    for voce in SOGLIE:
        if voce.chiave not in valori or non_eseguito(voce, valori) is not None:
            continue
        al_limite = limite(voce, valori, campioni_controllati)
        if oltre_il_limite(voce, valori[voce.chiave], al_limite):
            assert al_limite is not None
            spenti.append(ControlloSpento(
                voce.chiave, valori[voce.chiave], al_limite, voce.codici, voce.ferma))
    return spenti


def righe_di_avviso(spenti: list[ControlloSpento]) -> list[str]:
    """L'avviso sui controlli spenti, una riga per elemento; vuoto se non ce ne sono."""
    if not spenti:
        return []
    return [
        f"Avviso: {len(spenti)} "
        + ("controllo di qualita' e' spento" if len(spenti) == 1 else
           "controlli di qualita' sono spenti")
        + ": la soglia e' al valore che non ferma, o oltre.",
        *[f"  - {s.testo}" for s in spenti],
        "  Vanno bene per una prima esecuzione esplorativa. Poi il sottocomando soglie "
        "propone i valori dalle misure dell'esecuzione.",
    ]
