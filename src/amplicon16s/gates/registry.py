"""Registro dei gate di validazione e loro esecuzione in sequenza.

I quindici gate sono dichiarati qui in un solo posto, con l'ordine in cui
vanno eseguiti. L'ordine non è quello numerico per una ragione e una
eccezione:

* **G15 apre la sequenza.** Verifica la coerenza interna della configurazione,
  e un errore lì va scoperto prima di aprire un solo file: sarebbe assurdo
  ispezionare 960 archivi per poi fermarsi su una soglia incoerente.
* gli altri seguono l'ordine numerico, che è anche l'ordine delle dipendenze:
  senza ingressi leggibili non si aprono le tabelle, senza tabelle non si
  estraggono gli accession, senza accession non si confrontano gli insiemi.

**L'esecuzione si ferma al primo gate fallito.** I gate successivi
dipendono dai precedenti, e proseguire produrrebbe errori derivati che
confonderebbero la diagnosi invece di arricchirla. I gate non eseguiti sono
riportati come tali, non come superati.
"""

from __future__ import annotations

import time
from dataclasses import dataclass, field
from typing import Any, Final

from amplicon16s.gates.g01_g15 import (
    Avviso,
    Contesto,
    Violazione,
    _adatta,
    _g01_ingressi_leggibili,
    _g02_tabelle_apribili,
    _g03_join_ristretto,
    _g04_accession_estraibile,
    _g05_accession_univoci,
    _g06_insiemi_simmetrici,
    _g07_layout_single_end,
    _g08_lotto_coerente,
    _g09_troncamento_compatibile,
    _g10_primer_assente,
    _g11_classi_complete,
    _g12_riferimento_verificato,
    _g13_archivi_validi,
    _g14_risorse_disponibili,
    _g15_coerenza_configurazione,
)

__all__ = ["REGISTRO", "EsitoGate", "esegui_tutti", "nomi_dei_gate"]


@dataclass(frozen=True)
class EsitoGate:
    """Esito di un singolo gate."""

    gate: str
    descrizione: str
    eseguito: bool
    superato: bool
    violazioni: tuple[Violazione, ...] = ()
    avvisi: tuple[Avviso, ...] = ()
    secondi: float = 0.0

    def come_voce(self) -> dict[str, Any]:
        """Forma registrabile su disco e interrogabile."""
        return {
            "gate": self.gate,
            "descrizione": self.descrizione,
            "eseguito": self.eseguito,
            "superato": self.superato,
            "secondi": round(self.secondi, 3),
            "violazioni": [
                {"codice": v.codice, "dettaglio": v.dettaglio} for v in self.violazioni
            ],
            "avvisi": [
                {"codice": a.codice, "dettaglio": a.dettaglio} for a in self.avvisi
            ],
        }


@dataclass(frozen=True)
class _Voce:
    nome: str
    descrizione: str
    controllo: Any


REGISTRO: Final[tuple[_Voce, ...]] = (
    _Voce("G15", "coerenza interna della configurazione", _g15_coerenza_configurazione),
    _Voce("G01", "gli ingressi dichiarati esistono e sono leggibili", _g01_ingressi_leggibili),
    _Voce("G02", "le tabelle di metadati si aprono e hanno le colonne attese", _g02_tabelle_apribili),
    _Voce("G03", "il join resta ristretto alla tabella di assay", _adatta(_g03_join_ristretto)),
    _Voce("G04", "l'accession e' estraibile da ogni nome e non e' ambiguo", _adatta(_g04_accession_estraibile)),
    _Voce("G05", "gli accession sono univoci", _adatta(_g05_accession_univoci)),
    _Voce("G06", "gli insiemi dei file e dei metadati coincidono", _adatta(_g06_insiemi_simmetrici)),
    _Voce("G07", "il layout e' single-end", _g07_layout_single_end),
    _Voce("G08", "l'informazione di lotto e' coerente e le piastre sono plausibili", _g08_lotto_coerente),
    _Voce("G09", "filter.truncLen e' compatibile con le lunghezze osservate", _g09_troncamento_compatibile),
    _Voce("G10", "il primer non e' presente e il segnale atteso c'e'", _g10_primer_assente),
    _Voce("G11", "ogni campione ricade in una classe dichiarata", _adatta(_g11_classi_complete)),
    _Voce("G12", "il database tassonomico e' quello dichiarato", _g12_riferimento_verificato),
    _Voce("G13", "ogni file di letture e' un archivio valido", _g13_archivi_validi),
    _Voce("G14", "le risorse richieste sono disponibili", _g14_risorse_disponibili),
)


def nomi_dei_gate() -> tuple[str, ...]:
    """I nomi dei gate, nell'ordine di esecuzione."""
    return tuple(voce.nome for voce in REGISTRO)


def esegui_tutti(contesto: Contesto) -> list[EsitoGate]:
    """Esegue i gate in sequenza, fermandosi al primo che fallisce."""
    esiti: list[EsitoGate] = []
    fermato = False

    for voce in REGISTRO:
        if fermato:
            esiti.append(
                EsitoGate(voce.nome, voce.descrizione, eseguito=False, superato=False)
            )
            continue

        inizio = time.perf_counter()
        violazioni, avvisi = voce.controllo(contesto)
        durata = time.perf_counter() - inizio

        esiti.append(
            EsitoGate(
                gate=voce.nome,
                descrizione=voce.descrizione,
                eseguito=True,
                superato=not violazioni,
                violazioni=tuple(violazioni),
                avvisi=tuple(avvisi),
                secondi=durata,
            )
        )
        if violazioni:
            fermato = True

    return esiti
