"""Fase S0 — validazione iniziale.

È la barriera che gira prima di allocare qualunque calcolo costoso. Il suo
scopo è far emergere in pochi minuti un problema che altrimenti si scoprirebbe
dopo ore di elaborazione: un file corrotto, una colonna sbagliata, un
troncamento incompatibile con le lunghezze reali.

**Nessun gate legge un file intero.** Le verifiche sulle sequenze si fermano
alle prime ``qc.head_reads`` letture; l'integrità completa dei file è già
garantita dai checksum. È il vincolo che tiene la validazione nell'ordine dei
minuti.

La fase scrive in ``01_input_validation/`` l'esito di ogni gate, il crosswalk
e l'inventario, passando dal servizio di scrittura degli artefatti: ogni file
prodotto entra nel manifesto con il proprio checksum, come qualunque altro
artefatto della pipeline.
"""

from __future__ import annotations

import csv
import json
import time
from dataclasses import dataclass
from io import StringIO
from pathlib import Path
from typing import Any, Final

from amplicon16s.config.schema import Config
from amplicon16s.gates.g01_g15 import Contesto, ErroreGate, Violazione
from amplicon16s.gates.registry import EsitoGate, esegui_tutti
from amplicon16s.io_layer.artifacts import AlberoOutput, Fase
from amplicon16s.logging.logger import ottieni
from amplicon16s.metadata.models import ClasseCampione, Inventario

__all__ = ["RisultatoS0", "esegui_s0"]

NOME_ESITI: Final = "gates.json"
NOME_CROSSWALK: Final = "crosswalk.tsv"
NOME_INVENTARIO: Final = "inventario.json"
NOME_LETTURE: Final = "letture_ispezionate.tsv"


@dataclass(frozen=True)
class RisultatoS0:
    """Esito della fase, con gli artefatti prodotti."""

    esiti: tuple[EsitoGate, ...]
    inventario: Inventario | None
    artefatti: tuple[Path, ...]
    secondi: float

    @property
    def superata(self) -> bool:
        return all(e.superato for e in self.esiti)

    @property
    def falliti(self) -> tuple[EsitoGate, ...]:
        return tuple(e for e in self.esiti if e.eseguito and not e.superato)

    @property
    def avvisi(self) -> tuple[Any, ...]:
        return tuple(a for e in self.esiti for a in e.avvisi)


def _crosswalk_tsv(inventario: Inventario) -> str:
    """Il crosswalk in forma tabellare: una riga per campione."""
    buffer = StringIO()
    scrittore = csv.writer(buffer, delimiter="\t", lineterminator="\n")
    scrittore.writerow(
        ["accession", "sample_name", "classe", "materiale", "file",
         "posizione", "modulo", "piastra", "corsa"]
    )
    for campione in inventario:
        scrittore.writerow(
            [campione.accession, campione.nome, campione.classe.value,
             campione.materiale, campione.file.name if campione.file else "",
             campione.posizione or "", campione.modulo or "",
             campione.piastra or "", campione.corsa or ""]
        )
    return buffer.getvalue()


def _inventario_json(inventario: Inventario) -> str:
    conteggi = inventario.conteggi()
    documento = {
        "campioni": len(inventario),
        "conteggi": {classe.value: numero for classe, numero in conteggi.items()},
        # Il denominatore e' esplicito perche' a valle non lo si confonda con
        # il totale dei campioni: sono entrambi plausibili.
        "denominatore_prevalenza": inventario.denominatore_prevalenza(),
        "moduli": list(inventario.moduli),
        "piastre": list(inventario.piastre),
        "corse": list(inventario.corse),
        "campioni_per_piastra": inventario.campioni_per_piastra(),
        "senza_lotto": len(inventario.senza_lotto),
    }
    return json.dumps(documento, indent=2, ensure_ascii=False) + "\n"


def _letture_tsv(scansione: dict[str, Any], head_reads: int) -> str:
    """Le statistiche per file, con il limite della stima messo agli atti.

    Le lunghezze sono quelle delle prime ``head_reads`` letture, non di tutto
    il file: il minimo puo' quindi essere piu' alto del vero. La colonna
    ``file_esaurito`` dice per quali file la statistica copre tutto, e per
    quali e' una stima.
    """
    buffer = StringIO()
    scrittore = csv.writer(buffer, delimiter="\t", lineterminator="\n")
    scrittore.writerow(
        ["accession", "file", "letture_esaminate", "file_esaurito",
         "lunghezza_minima", "lunghezza_massima", "frazione_primer",
         "frazione_motivo", "errore"]
    )
    for accession in sorted(scansione):
        s = scansione[accession]
        scrittore.writerow(
            [accession, s.nome, s.letture_esaminate, "si" if s.esaurito else "no",
             s.lunghezza_minima if s.lunghezza_minima is not None else "",
             s.lunghezza_massima if s.lunghezza_massima is not None else "",
             f"{s.frazione_primer:.6f}", f"{s.frazione_motivo:.6f}",
             s.errore or ""]
        )
    return buffer.getvalue()


def _esiti_json(esiti: tuple[EsitoGate, ...], secondi: float) -> str:
    documento = {
        "superata": all(e.superato for e in esiti),
        "secondi": round(secondi, 3),
        "gate": [e.come_voce() for e in esiti],
    }
    return json.dumps(documento, indent=2, ensure_ascii=False) + "\n"


def esegui_s0(config: Config, *, solleva: bool = True) -> RisultatoS0:
    """Esegue i quindici gate e registra l'esito in ``01_input_validation/``.

    Con ``solleva`` falso restituisce il risultato anche quando un gate
    fallisce, invece di sollevare: serve a chi vuole ispezionare l'esito
    completo, per esempio il comando ``validate``.
    """
    log = ottieni("s0")
    inizio = time.perf_counter()

    contesto = Contesto(config)
    esiti = tuple(esegui_tutti(contesto))
    secondi = time.perf_counter() - inizio

    # L'inventario esiste solo se i gate che lo costruiscono sono passati.
    inventario: Inventario | None = None
    if all(e.superato for e in esiti if e.gate in ("G03", "G04", "G05", "G06", "G11")):
        try:
            inventario = contesto.inventario
        except Exception:  # noqa: BLE001 - l'inventario e' accessorio al referto
            inventario = None

    albero = AlberoOutput(config.io.out_root)
    artefatti = [
        albero.scrivi_testo(
            Fase.INPUT_VALIDATION, NOME_ESITI, _esiti_json(esiti, secondi)
        ).percorso
    ]
    if inventario is not None:
        artefatti.append(
            albero.scrivi_testo(
                Fase.INPUT_VALIDATION, NOME_CROSSWALK, _crosswalk_tsv(inventario)
            ).percorso
        )
        artefatti.append(
            albero.scrivi_testo(
                Fase.INPUT_VALIDATION, NOME_INVENTARIO, _inventario_json(inventario)
            ).percorso
        )
    if contesto.scansione_gia_fatta:
        artefatti.append(
            albero.scrivi_testo(
                Fase.INPUT_VALIDATION,
                NOME_LETTURE,
                _letture_tsv(contesto.scansione, config.qc.head_reads),
            ).percorso
        )

    for esito in esiti:
        for avviso in esito.avvisi:
            log.warning(
                avviso.dettaglio,
                extra={"gate": esito.gate, "codice": avviso.codice, "fase": "S0"},
            )

    risultato = RisultatoS0(
        esiti=esiti,
        inventario=inventario,
        artefatti=tuple(artefatti),
        secondi=secondi,
    )

    if risultato.superata:
        log.info(
            "S0 superata",
            extra={
                "fase": "S0",
                "secondi": round(secondi, 3),
                "campioni": len(inventario) if inventario else 0,
                "avvisi": len(risultato.avvisi),
            },
        )
        return risultato

    falliti = risultato.falliti
    for esito in falliti:
        for violazione in esito.violazioni:
            log.error(
                violazione.dettaglio,
                extra={"gate": esito.gate, "codice": violazione.codice, "fase": "S0"},
            )

    if solleva:
        primo = falliti[0]
        raise ErroreGate(primo.gate, primo.violazioni)
    return risultato
