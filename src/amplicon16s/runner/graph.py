"""Il grafo delle fasi: ordine, dipendenze, fasi facoltative.

Le fasi sono quindici, da S0 a S14, e sono dichiarate qui in un solo posto,
nell'ordine di esecuzione. Il grafo non esegue nulla e non conosce il
filesystem: dice quali fasi esistono, in quale cartella scrivono, da quali
fasi dipendono e quando sono attive. Stabilire quali sono già completate è
compito di :mod:`amplicon16s.runner.project`, che lo legge dal disco.

**Le dipendenze sono dipendenze di dato.** Una fase dipende da quelle di cui
consuma gli artefatti, non da tutte quelle che la precedono: S2 filtra le
letture indicate dall'inventario di S0 e non legge i profili di qualità di S1,
quindi ricalcolare S1 non invalida S2 e tutto ciò che segue. La ripresa rifà
così solo il lavoro che dipende davvero da ciò che è cambiato. Il rovescio è
che una dipendenza dimenticata produce un risultato stantio senza alcun
errore: quando una fase viene realizzata, le sue dipendenze qui vanno
verificate contro gli artefatti che legge davvero.

**Due vincoli sono rappresentati esplicitamente.**

* S9, la filogenesi, è **facoltativa**: con ``phylo.enabled`` falso è
  disattivata. Una fase disattivata non è una fase incompleta, e chi ne
  dipende la ignora: S10 senza albero dipende solo dalle altre fasi.
* S12 **precede obbligatoriamente** S13: la decontaminazione va fatta prima del
  filtro di prevalenza. Il vincolo è una dipendenza diretta che il grafo
  verifica alla costruzione, e non può essere resa facoltativa. Porta un
  codice proprio, ``E-S13-01``: violarlo è un errore di metodo, distinto dalla
  violazione di una dipendenza qualunque, che è ``E-GRAFO-01``.

**Limite noto: il costo della valutazione.** Stabilire se una fase è
conclusa ricalcola il checksum di ogni suo artefatto, perché solo il
checksum garantisce che un artefatto alterato faccia rieseguire la fase: una
scorciatoia su dimensione e data di modifica non vedrebbe un file alterato
che le conserva. Dentro una valutazione ogni checksum è calcolato una volta
sola, ma una valutazione completa li calcola tutti: con le letture filtrate
di S2 e i risultati dell'inferenza di S4, che arrivano ai gigabyte, costerà
secondi.
"""

from __future__ import annotations

from collections.abc import Callable, Iterable, Iterator
from dataclasses import dataclass
from enum import StrEnum
from typing import Final

from amplicon16s.config.schema import Config
from amplicon16s.errors.catalog import voce
from amplicon16s.io_layer.artifacts import Fase

__all__ = [
    "GRAFO",
    "PRECEDENZE_OBBLIGATORIE",
    "Grafo",
    "Nodo",
    "Passo",
]


class Passo(StrEnum):
    """Le quindici fasi della pipeline, nell'ordine di esecuzione."""

    S0 = "S0"
    S1 = "S1"
    S2 = "S2"
    S3 = "S3"
    S4 = "S4"
    S5 = "S5"
    S6 = "S6"
    S7 = "S7"
    S8 = "S8"
    S9 = "S9"
    S10 = "S10"
    S11 = "S11"
    S12 = "S12"
    S13 = "S13"
    S14 = "S14"


@dataclass(frozen=True)
class Nodo:
    """Una fase nel grafo."""

    passo: Passo
    descrizione: str
    #: Cartella di output in cui la fase scrive. Piu' fasi possono
    #: condividerla: il completamento resta comunque per fase.
    cartella: Fase
    #: Fasi di cui la fase consuma gli artefatti.
    dipendenze: tuple[Passo, ...] = ()
    #: Condizione di attivazione; ``None`` per le fasi sempre attive.
    attiva_se: Callable[[Config], bool] | None = None
    #: Il parametro che governa la condizione, per i messaggi.
    parametro_attivazione: str | None = None

    @property
    def facoltativa(self) -> bool:
        return self.attiva_se is not None

    def attiva(self, config: Config) -> bool:
        return self.attiva_se is None or self.attiva_se(config)


#: Terne (prima, dopo, codice): la prima fase deve essere completata prima che
#: la seconda possa girare, per una ragione di metodo e non solo di dato, e
#: violarlo solleva il codice indicato.
PRECEDENZE_OBBLIGATORIE: Final[tuple[tuple[Passo, Passo, str], ...]] = (
    # La prevalenza va calcolata su dati gia' decontaminati: un contaminante
    # diffuso supererebbe il filtro proprio perche' e' ovunque.
    (Passo.S12, Passo.S13, "E-S13-01"),
)


class Grafo:
    """L'insieme ordinato delle fasi, verificato alla costruzione."""

    def __init__(
        self,
        nodi: Iterable[Nodo],
        precedenze: Iterable[tuple[Passo, Passo, str]] = PRECEDENZE_OBBLIGATORIE,
    ) -> None:
        self._nodi: dict[Passo, Nodo] = {}
        for nodo in nodi:
            if nodo.passo in self._nodi:
                raise ValueError(f"fase dichiarata due volte: {nodo.passo}")
            # Ogni dipendenza deve precedere la fase che la usa: l'ordine di
            # dichiarazione e' quindi un ordine di esecuzione valido, e un
            # ciclo e' impossibile per costruzione.
            for dipendenza in nodo.dipendenze:
                if dipendenza not in self._nodi:
                    raise ValueError(
                        f"{nodo.passo} dipende da {dipendenza}, che non la precede"
                    )
            self._nodi[nodo.passo] = nodo

        self.precedenze = tuple(precedenze)
        for prima, dopo, codice in self.precedenze:
            voce(codice)  # un codice inesistente e' un difetto del grafo
            if prima not in self._nodi or dopo not in self._nodi:
                raise ValueError(f"precedenza su una fase assente: {prima} -> {dopo}")
            if prima not in self._nodi[dopo].dipendenze:
                raise ValueError(
                    f"{prima} deve precedere {dopo}: dev'essere una sua dipendenza diretta"
                )
            if self._nodi[prima].facoltativa or self._nodi[dopo].facoltativa:
                raise ValueError(
                    f"la precedenza {prima} -> {dopo} non puo' coinvolgere una fase "
                    "facoltativa: disattivarla la aggirerebbe"
                )

    def __iter__(self) -> Iterator[Nodo]:
        return iter(self._nodi.values())

    def __len__(self) -> int:
        return len(self._nodi)

    def nodo(self, passo: Passo) -> Nodo:
        return self._nodi[passo]

    def ordine(self) -> tuple[Passo, ...]:
        """Tutte le fasi, nell'ordine di esecuzione."""
        return tuple(self._nodi)

    def attive(self, config: Config) -> tuple[Passo, ...]:
        """Le fasi che questa configurazione esegue, in ordine."""
        return tuple(p for p, n in self._nodi.items() if n.attiva(config))

    def dipendenze_attive(self, passo: Passo, config: Config) -> tuple[Passo, ...]:
        """Le dipendenze della fase, meno quelle disattivate dalla configurazione."""
        return tuple(
            d for d in self._nodi[passo].dipendenze if self._nodi[d].attiva(config)
        )

    def discendenti(self, passo: Passo) -> tuple[Passo, ...]:
        """Le fasi che dipendono, anche indirettamente, da quella data."""
        raggiunti = {passo}
        for nodo in self._nodi.values():
            if raggiunti & set(nodo.dipendenze):
                raggiunti.add(nodo.passo)
        raggiunti.discard(passo)
        return tuple(p for p in self._nodi if p in raggiunti)


def _filogenesi_attiva(config: Config) -> bool:
    return config.phylo.enabled


# S7 non ha una cartella propria: filtra per lunghezza la tabella senza
# chimere di S6 e ne raccoglie il tracciamento delle letture, quindi scrive
# accanto a S6 in 07_chimera.
GRAFO: Final = Grafo(
    [
        Nodo(Passo.S0, "validazione degli ingressi", Fase.INPUT_VALIDATION),
        Nodo(Passo.S1, "profilo di qualita'", Fase.QC_PROFILES, (Passo.S0,)),
        Nodo(Passo.S2, "filtro e troncamento", Fase.FILTERED, (Passo.S0,)),
        Nodo(Passo.S3, "modello d'errore", Fase.ERROR_MODELS, (Passo.S2,)),
        Nodo(
            Passo.S4, "inferenza delle varianti", Fase.ASV_INFERENCE,
            (Passo.S2, Passo.S3),
        ),
        Nodo(Passo.S5, "tabella delle sequenze", Fase.SEQTAB, (Passo.S4,)),
        Nodo(Passo.S6, "rimozione delle chimere", Fase.CHIMERA, (Passo.S5,)),
        Nodo(
            Passo.S7, "filtro di lunghezza e tracciamento delle letture",
            Fase.CHIMERA,
            # Il tracciamento raccoglie i conteggi di ogni passo precedente.
            (Passo.S2, Passo.S4, Passo.S5, Passo.S6),
        ),
        Nodo(Passo.S8, "assegnazione tassonomica", Fase.TAXONOMY, (Passo.S7,)),
        Nodo(
            Passo.S9, "filogenesi", Fase.PHYLOGENY, (Passo.S7,),
            attiva_se=_filogenesi_attiva,
            parametro_attivazione="phylo.enabled",
        ),
        Nodo(
            Passo.S10, "assemblaggio dell'oggetto integrato", Fase.PHYLOSEQ,
            # S0 per i metadati dei campioni, S9 solo se attiva.
            (Passo.S0, Passo.S7, Passo.S8, Passo.S9),
        ),
        Nodo(
            Passo.S11, "validazione dai controlli positivi", Fase.CONTROLS,
            (Passo.S10,),
        ),
        Nodo(
            Passo.S12, "decontaminazione", Fase.CONTROLS,
            (Passo.S10, Passo.S11),
        ),
        Nodo(
            Passo.S13, "filtri tassonomici e di prevalenza", Fase.FINAL,
            (Passo.S12,),
        ),
        Nodo(
            Passo.S14, "serializzazione e validazione finale", Fase.FINAL,
            (Passo.S13,),
        ),
    ]
)
