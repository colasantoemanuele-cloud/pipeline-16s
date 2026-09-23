"""Eccezioni della pipeline, legate al catalogo dei codici.

Ogni errore della pipeline porta un codice, e dal codice discendono il
messaggio operativo e la categoria di gestione: chi solleva l'errore dichiara
*cosa* è successo, non *come* va gestito. La gestione è una proprietà del
codice, decisa una volta nel catalogo.

Le sottoclassi corrispondono alle quattro categorie, così un chiamante può
distinguere i casi con ``except`` invece di ispezionare una stringa. La
scelta della sottoclasse non è però lasciata a chi solleva l'errore:
:func:`errore` la deriva dal catalogo, e costruire direttamente la sottoclasse
sbagliata viene respinto.
"""

from __future__ import annotations

from typing import Any

from amplicon16s.errors.catalog import Categoria, VoceCatalogo, voce

__all__ = [
    "DegradazioneRichiesta",
    "ErrorePipeline",
    "ErroreRevisioneUmana",
    "ErroreRitentabile",
    "ErroreRitentabileConRevisione",
    "errore",
]


class ErrorePipeline(Exception):
    """Errore identificato da un codice del catalogo.

    ``dettaglio`` aggiunge ciò che il catalogo non può sapere: quale file,
    quale campione, quale valore. ``contesto`` raccoglie dati strutturati utili
    al log, che non appartengono al messaggio per un lettore umano.
    """

    #: Categorie che questa sottoclasse rappresenta. La base le accetta tutte.
    categorie_ammesse: tuple[Categoria, ...] = tuple(Categoria)

    def __init__(
        self,
        codice: str,
        dettaglio: str = "",
        **contesto: Any,
    ) -> None:
        self.voce: VoceCatalogo = voce(codice)
        self.codice = codice
        self.dettaglio = dettaglio
        self.contesto = contesto

        if self.voce.categoria not in type(self).categorie_ammesse:
            raise TypeError(
                f"{codice} e' classificato come {self.voce.categoria.value} e non "
                f"puo' essere sollevato come {type(self).__name__}; usa "
                f"errore({codice!r}) per ottenere la classe corretta"
            )

        super().__init__(self._testo())

    def _testo(self) -> str:
        parti = [f"[{self.codice}] {self.voce.sintesi}"]
        if self.dettaglio:
            parti.append(self.dettaglio)
        parti.append(self.voce.azione)
        return " ".join(parti)

    @property
    def categoria(self) -> Categoria:
        return self.voce.categoria

    @property
    def fase(self) -> str:
        return self.voce.fase

    @property
    def ammette_retry(self) -> bool:
        return self.voce.ammette_retry

    def come_evento(self) -> dict[str, Any]:
        """Forma strutturata dell'errore, per il log interrogabile."""
        return {
            "codice": self.codice,
            "fase": self.fase,
            "categoria": self.categoria.value,
            "sintesi": self.voce.sintesi,
            "azione": self.voce.azione,
            "dettaglio": self.dettaglio or None,
            **self.contesto,
        }


class ErroreRevisioneUmana(ErrorePipeline):
    """L'esecuzione si ferma: la decisione spetta all'operatore."""

    categorie_ammesse = (Categoria.REVISIONE_UMANA,)


class ErroreRitentabile(ErrorePipeline):
    """Si ritenta entro il limite di ``retry.max_attempts``."""

    categorie_ammesse = (Categoria.RETRY_AUTOMATICO,)


class ErroreRitentabileConRevisione(ErrorePipeline):
    """Si ritenta e, se fallisce ancora, l'esecuzione si ferma."""

    categorie_ammesse = (Categoria.RETRY_POI_REVISIONE,)


class DegradazioneRichiesta(ErrorePipeline):
    """Si prosegue con un comportamento di ripiego, registrandolo.

    È un'eccezione perché il punto in cui la condizione viene rilevata non è
    quello che sa come ripiegare: viene sollevata e raccolta da chi governa la
    fase, che registra la degradazione e prosegue.
    """

    categorie_ammesse = (Categoria.DEGRADAZIONE_AUTOMATICA,)


_PER_CATEGORIA: dict[Categoria, type[ErrorePipeline]] = {
    Categoria.REVISIONE_UMANA: ErroreRevisioneUmana,
    Categoria.RETRY_AUTOMATICO: ErroreRitentabile,
    Categoria.RETRY_POI_REVISIONE: ErroreRitentabileConRevisione,
    Categoria.DEGRADAZIONE_AUTOMATICA: DegradazioneRichiesta,
}


def errore(codice: str, dettaglio: str = "", **contesto: Any) -> ErrorePipeline:
    """Costruisce l'eccezione della categoria giusta per il codice dato.

    È il modo previsto per sollevare un errore della pipeline::

        raise errore("E-S2-01", f"campione {nome} azzerato", campione=nome)
    """
    classe = _PER_CATEGORIA[voce(codice).categoria]
    return classe(codice, dettaglio, **contesto)
