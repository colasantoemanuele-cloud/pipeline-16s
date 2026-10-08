"""Modelli dell'inventario dei campioni.

L'inventario è ciò che la pipeline sa dei campioni prima di toccare una
lettura: quale file corrisponde a quale campione, a quale classe appartiene,
in quale lotto è stato trattato. Su questo si reggono la decontaminazione, il
modello d'errore e il denominatore dei filtri di prevalenza.

**La chiave è l'accession, non il nome del campione.** Il nome è un'etichetta
umana: può ripetersi fra repliche ed essere riusato per i controlli. Usarlo
come chiave produce appaiamenti fra file e metadati che nessun controllo a
valle è in grado di rilevare, perché il risultato resta plausibile.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from enum import StrEnum
from typing import Final
from pathlib import Path

__all__ = ["CLASSI_CONTROLLATE", "Campione", "ClasseCampione", "Inventario"]


class ClasseCampione(StrEnum):
    """Classe di un campione, ricavata dalla colonna indicata da ``ctrl.column``."""

    BIOLOGICO = "biologico"
    CONTROLLO_POSITIVO = "controllo_positivo"
    CONTROLLO_NEGATIVO = "controllo_negativo"


#: Le classi su cui si applicano i controlli di qualita' sul segnale: G10 (il
#: motivo conservato), S2 (campioni azzerati e letture perse), S6 (frazione
#: chimerica), S8 (copertura del phylum). Sono quelle da cui ci si attende il
#: segnale del bersaglio. I controlli negativi ne sono esclusi per principio, non
#: per misura: un bianco non contiene il materiale studiato, e un bianco
#: azzerato dal filtro o con poco segnale e' un bianco pulito, non un difetto dei
#: parametri. Definite una volta sola: e' una decisione di metodo.
CLASSI_CONTROLLATE: Final[tuple[ClasseCampione, ...]] = (
    ClasseCampione.BIOLOGICO,
    ClasseCampione.CONTROLLO_POSITIVO,
)


@dataclass(frozen=True)
class Campione:
    """Un campione dell'assay, con la sua chiave, la sua classe e il suo lotto."""

    #: Chiave del join: l'accession estratto dal nome del file grezzo.
    accession: str
    #: Nome biologico. Serve a leggere i risultati, mai come chiave.
    nome: str
    classe: ClasseCampione
    #: Valore grezzo della colonna di classificazione, conservato perché un
    #: risultato deve restare riconducibile al dato che l'ha prodotto.
    materiale: str
    #: File di letture corrispondente, quando l'inventario è costruito con i file.
    file: Path | None = None
    #: Posizione dichiarata nei metadati, da cui si deriva il modulo.
    posizione: str | None = None
    #: Modulo del campione: il valore dichiarato dal file di arricchimento se
    #: c'è, altrimenti quello derivato dalla posizione. ``None`` per le
    #: posizioni che non sono superfici, in entrambe le modalità.
    modulo: str | None = None
    #: Lotto di estrazione e corsa di sequenziamento. Restano ``None`` quando
    #: il file di arricchimento non c'è o non copre questo campione.
    piastra: str | None = None
    corsa: str | None = None

    @property
    def ha_lotto(self) -> bool:
        """Vero se il campione ha almeno una fra piastra e corsa."""
        return self.piastra is not None or self.corsa is not None


@dataclass(frozen=True)
class Inventario:
    """L'insieme dei campioni dell'assay, indicizzato per accession."""

    campioni: tuple[Campione, ...]

    def __len__(self) -> int:
        return len(self.campioni)

    def __iter__(self):
        return iter(self.campioni)

    def __getitem__(self, accession: str) -> Campione:
        for campione in self.campioni:
            if campione.accession == accession:
                return campione
        raise KeyError(f"accession non presente nell'inventario: {accession}")

    @property
    def accessioni(self) -> tuple[str, ...]:
        """Gli accession dei campioni, nell'ordine dell'inventario."""
        return tuple(c.accession for c in self.campioni)

    def di_classe(self, classe: ClasseCampione) -> tuple[Campione, ...]:
        """I campioni della classe indicata, nell'ordine dell'inventario."""
        return tuple(c for c in self.campioni if c.classe is classe)

    @property
    def biologici(self) -> tuple[Campione, ...]:
        """I campioni biologici."""
        return self.di_classe(ClasseCampione.BIOLOGICO)

    def denominatore_prevalenza(self) -> int:
        """Numero di campioni su cui si calcolano le prevalenze.

        Sono i soli campioni biologici, non il totale. I controlli non sono
        campioni dell'ambiente studiato: includerli gonfierebbe il
        denominatore e abbasserebbe ogni prevalenza, rendendo piu' facile
        superare la soglia di un filtro che serve proprio a non farlo.

        Esiste come metodo con un nome proprio perche' a valle non si possa
        usare per sbaglio ``len(inventario)``: i due numeri sono entrambi
        plausibili, e sbagliarli non produce alcun errore visibile.
        """
        return len(self.biologici)

    def conteggi(self) -> dict[ClasseCampione, int]:
        """Il numero di campioni per classe, con zero per le classi assenti."""
        conteggio = Counter(c.classe for c in self.campioni)
        return {classe: conteggio.get(classe, 0) for classe in ClasseCampione}

    @property
    def piastre(self) -> tuple[str, ...]:
        """Le piastre distinte, in ordine alfabetico."""
        return tuple(sorted({c.piastra for c in self.campioni if c.piastra is not None}))

    @property
    def corse(self) -> tuple[str, ...]:
        """Le corse di sequenziamento distinte, in ordine alfabetico."""
        return tuple(sorted({c.corsa for c in self.campioni if c.corsa is not None}))

    @property
    def moduli(self) -> tuple[str, ...]:
        """I moduli distinti, in ordine alfabetico."""
        return tuple(sorted({c.modulo for c in self.campioni if c.modulo is not None}))

    @property
    def senza_lotto(self) -> tuple[Campione, ...]:
        """Campioni per cui il lotto non è noto.

        Senza il file di arricchimento sono tutti: è il funzionamento a corsa
        singola, non un errore.
        """
        return tuple(c for c in self.campioni if not c.ha_lotto)

    def campioni_per_piastra(self) -> dict[str, int]:
        """Il numero di campioni di ciascuna piastra."""
        return dict(Counter(c.piastra for c in self.campioni if c.piastra is not None))
