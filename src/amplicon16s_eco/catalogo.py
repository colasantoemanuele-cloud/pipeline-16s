"""Catalogo dei codici del pacchetto delle analisi ecologiche.

Ogni codice ``E-ECO-NN`` porta una sintesi (che cosa e' successo) e un'azione
(che cosa correggere, oppure che cosa il pacchetto ha fatto da se'). I codici
sono di due specie:

* **rifiuto**: l'analisi non parte, o si ferma. La causa va corretta nella
  configurazione, nell'oggetto o nella riga di comando;
* **avviso**: l'analisi prosegue, e il riepilogo riporta il codice con il
  dettaglio. Un avviso non cambia un metodo: dichiara una condizione che chi
  legge i risultati deve conoscere.

Un rifiuto e un avviso non condividono mai lo stesso codice. Il catalogo e'
separato da quello della pipeline: il pacchetto lavora a valle dell'oggetto
consegnato e non partecipa alle fasi che lo producono.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import StrEnum
from typing import Final

__all__ = ["CATALOGO", "ErroreEco", "Specie", "VoceEco", "descrivi", "voce"]


class Specie(StrEnum):
    """Che cosa comporta un codice: l'arresto oppure una dichiarazione."""

    RIFIUTO = "rifiuto"
    AVVISO = "avviso"


@dataclass(frozen=True)
class VoceEco:
    """Un codice del pacchetto con il suo significato."""

    codice: str
    specie: Specie
    #: Che cosa e' successo, in una riga.
    sintesi: str
    #: Che cosa correggere, o che cosa il pacchetto ha fatto da se'.
    azione: str


_R = Specie.RIFIUTO
_A = Specie.AVVISO

_VOCI: Final[tuple[VoceEco, ...]] = (
    VoceEco(
        "E-ECO-01", _R,
        "La configurazione dell'analisi non e' valida.",
        "Correggi il parametro indicato: una chiave sconosciuta e' un refuso o un "
        "parametro non previsto, un valore fuori dominio va riportato entro i "
        "limiti descritti nel modello config/eco.example.yaml.",
    ),
    VoceEco(
        "E-ECO-02", _R,
        "Manca un parametro obbligatorio.",
        "Dichiara il parametro indicato nella configurazione dell'analisi: cio' "
        "che dipende dallo studio non ha un valore predefinito, a cominciare "
        "dalla variabile biologica (design.variable).",
    ),
    VoceEco(
        "E-ECO-03", _R,
        "E' richiesto un metodo che il pacchetto non realizza.",
        "Scegli uno dei metodi ammessi, elencati nel messaggio.",
    ),
    VoceEco(
        "E-ECO-04", _R,
        "L'oggetto di partenza non e' utilizzabile.",
        "Indica con --object l'oggetto phyloseq finale consegnato dalla pipeline "
        "(ps_final.rds), con la tabella dei conteggi, quella dei campioni e la "
        "tassonomia.",
    ),
    VoceEco(
        "E-ECO-05", _R,
        "Una colonna dichiarata non esiste nella tabella dei campioni dell'oggetto.",
        "Correggi il nome della colonna: quelle disponibili sono elencate nel "
        "messaggio.",
    ),
    VoceEco(
        "E-ECO-06", _R,
        "E' richiesta una distanza UniFrac ma l'oggetto non ha un albero radicato.",
        "Togli le distanze UniFrac da beta.distances, oppure usa un oggetto che "
        "contenga l'albero filogenetico con la radice.",
    ),
    VoceEco(
        "E-ECO-07", _R,
        "Il rango tassonomico dichiarato non e' fra quelli dell'oggetto.",
        "Porta comp.rank a uno dei ranghi elencati nel messaggio.",
    ),
    VoceEco(
        "E-ECO-08", _R,
        "Il sottoinsieme dichiarato non seleziona campioni, o nomina un valore "
        "che non compare.",
        "Correggi design.subset: ogni valore ammesso deve comparire nella sua "
        "colonna, e almeno due campioni devono restare.",
    ),
    VoceEco(
        "E-ECO-09", _R,
        "La cartella di uscita non e' utilizzabile.",
        "Indica con --out una cartella assente o vuota, fuori dalla cartella che "
        "contiene l'oggetto di partenza: le uscite dell'analisi non si mescolano "
        "con gli artefatti della pipeline ne' con un'analisi precedente.",
    ),
    VoceEco(
        "E-ECO-10", _R,
        "La profondita' di rarefazione dichiarata lascia meno di due campioni.",
        "Abbassa alpha.rarefy_depth, oppure toglilo per usare la profondita' "
        "minima fra i campioni analizzati.",
    ),
    VoceEco(
        "E-ECO-11", _A,
        "Un gruppo ha meno campioni di stat.min_group_size.",
        "Il gruppo esce dai test di quella variabile; resta nelle tabelle e nei "
        "grafici descrittivi.",
    ),
    VoceEco(
        "E-ECO-12", _A,
        "Restano meno di due gruppi per una variabile.",
        "I test di quella variabile non si eseguono; le uscite descrittive "
        "restano valide.",
    ),
    VoceEco(
        "E-ECO-13", _A,
        "Alcuni campioni hanno un valore mancante nella variabile analizzata.",
        "Quei campioni escono dalle analisi di quella variabile e sono elencati "
        "nella tabella dei campioni con il motivo.",
    ),
    VoceEco(
        "E-ECO-14", _A,
        "Alcuni campioni hanno meno letture della profondita' di rarefazione "
        "dichiarata.",
        "Quei campioni escono dall'alfa diversita' e sono elencati nella tabella "
        "dei campioni con il motivo.",
    ),
    VoceEco(
        "E-ECO-15", _A,
        "Alcuni campioni non hanno letture.",
        "Quei campioni escono da tutte le analisi e sono elencati nella tabella "
        "dei campioni con il motivo.",
    ),
    VoceEco(
        "E-ECO-16", _A,
        "Un calcolo in R ha emesso un avviso non previsto.",
        "Il testo dell'avviso e' nel dettaglio: valuta se riguarda i risultati "
        "che usi.",
    ),
    VoceEco(
        "E-ECO-90", _R,
        "Il calcolo in R non si e' concluso.",
        "Leggi il messaggio di R riportato nel dettaglio. Le analisi vanno "
        "eseguite nell'ambiente della pipeline (il container), che contiene i "
        "pacchetti richiesti.",
    ),
    VoceEco(
        "E-ECO-91", _R,
        "L'oggetto di partenza e' cambiato durante l'analisi.",
        "Le uscite non sono affidabili: verifica che nessun altro processo "
        "scriva l'oggetto e ripeti l'analisi in una cartella di uscita nuova.",
    ),
)

#: Il catalogo, per codice.
CATALOGO: Final[dict[str, VoceEco]] = {v.codice: v for v in _VOCI}
assert len(CATALOGO) == len(_VOCI), "codice ripetuto nel catalogo"


def voce(codice: str) -> VoceEco:
    """La voce del codice; un codice assente e' un difetto di chi lo solleva."""
    return CATALOGO[codice]


class ErroreEco(Exception):
    """Uno o piu' rifiuti del catalogo: l'analisi non parte, o si ferma.

    ``rifiuti`` e' l'elenco delle coppie (codice, dettaglio). Il messaggio
    riporta per ciascuna la sintesi, il dettaglio e l'azione.
    """

    def __init__(self, rifiuti: list[tuple[str, str]]) -> None:
        for codice, _ in rifiuti:
            if voce(codice).specie is not Specie.RIFIUTO:
                raise TypeError(f"{codice} e' un avviso e non puo' fermare l'analisi")
        self.rifiuti = list(rifiuti)
        super().__init__("\n".join(descrivi(c, d) for c, d in rifiuti))

    @property
    def codici(self) -> list[str]:
        """I codici dei rifiuti, nell'ordine in cui sono stati rilevati."""
        return [codice for codice, _ in self.rifiuti]


def descrivi(codice: str, dettaglio: str = "") -> str:
    """Il testo completo di un codice: sintesi, dettaglio se presente, azione."""
    v = voce(codice)
    parti = [f"[{codice}] {v.sintesi}"]
    if dettaglio:
        parti.append(dettaglio)
    parti.append(v.azione)
    return " ".join(parti)
