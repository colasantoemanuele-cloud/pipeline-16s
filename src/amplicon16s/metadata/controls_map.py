"""Classificazione dei campioni in biologici e controlli.

La classe non si deduce dal nome del campione: si legge dalla colonna
dichiarata in ``ctrl.column``, confrontandone il valore con gli elenchi di
etichette dichiarati in configurazione. Dedurla dal nome funzionerebbe su
questo dataset e fallirebbe sul successivo, che userà altre convenzioni.

Un valore che non compare in nessuno dei tre elenchi non viene attribuito a
una categoria di ripiego: resta **non mappato**, e il gate G11 ferma
l'esecuzione. È la differenza fra non sapere e credere di sapere.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

from amplicon16s.metadata.models import ClasseCampione

if TYPE_CHECKING:
    from amplicon16s.config.schema import Ctrl

__all__ = ["MappaControlli"]


@dataclass(frozen=True)
class MappaControlli:
    """Corrispondenza fra le etichette dichiarate e le classi dei campioni."""

    colonna: str
    per_etichetta: dict[str, ClasseCampione]
    #: La regola di riclassificazione in controllo negativo: la colonna e i
    #: valori normalizzati che la attivano; nessuna regola con l'insieme vuoto.
    colonna_negativi: str | None = None
    valori_negativi: frozenset[str] = frozenset()

    @classmethod
    def da_configurazione(cls, ctrl: Ctrl) -> MappaControlli:
        """Costruisce la mappa dagli elenchi del gruppo ``ctrl``.

        Che le tre categorie siano disgiunte è già garantito dallo schema, che
        rifiuta una configurazione in cui la stessa etichetta compaia in due
        elenchi: qui non serve ricontrollarlo.
        """
        mappa: dict[str, ClasseCampione] = {}
        for etichette, classe in (
            (ctrl.biological_values, ClasseCampione.BIOLOGICO),
            (ctrl.positive_values, ClasseCampione.CONTROLLO_POSITIVO),
            (ctrl.blank_values, ClasseCampione.CONTROLLO_NEGATIVO),
        ):
            for etichetta in etichette:
                mappa[cls._normalizza(etichetta)] = classe
        return cls(
            colonna=ctrl.column,
            per_etichetta=mappa,
            colonna_negativi=ctrl.blank_override_column,
            valori_negativi=frozenset(cls._normalizza(v) for v in ctrl.blank_override_values),
        )

    @staticmethod
    def _normalizza(valore: str) -> str:
        """Confronto insensibile a maiuscole e spazi ai bordi.

        Le tabelle di metadati sono compilate a mano: "Positive Control" e
        "positive control" sono la stessa cosa, e far fallire l'esecuzione su
        una maiuscola sarebbe pedanteria, non rigore. La distinzione fra
        etichette diverse resta intatta.
        """
        return valore.strip().casefold()

    def classifica(self, valore: str | None) -> ClasseCampione | None:
        """Classe corrispondente al valore, o ``None`` se non è mappato."""
        if valore is None:
            return None
        return self.per_etichetta.get(self._normalizza(valore))

    def riclassifica(
        self, classe: ClasseCampione, riga: dict[str, str]
    ) -> ClasseCampione:
        """La classe dopo la regola di riclassificazione in controllo negativo.

        ``riga`` e' la riga della tabella campioni di studio del campione. Se il
        valore della colonna della regola e' fra quelli dichiarati, il campione
        e' un controllo negativo qualunque sia il materiale: un tampone mai
        aperto dichiarato come superficie non ha campionato alcuna superficie.
        """
        if self.colonna_negativi is None or not self.valori_negativi:
            return classe
        valore = riga.get(self.colonna_negativi)
        if valore is not None and self._normalizza(valore) in self.valori_negativi:
            return ClasseCampione.CONTROLLO_NEGATIVO
        return classe

    @property
    def etichette(self) -> tuple[str, ...]:
        """Le etichette mappate, normalizzate e in ordine alfabetico."""
        return tuple(sorted(self.per_etichetta))

    def etichette_di(self, classe: ClasseCampione) -> tuple[str, ...]:
        """Le etichette normalizzate assegnate alla classe indicata, in ordine
        alfabetico.
        """
        return tuple(sorted(e for e, c in self.per_etichetta.items() if c is classe))
