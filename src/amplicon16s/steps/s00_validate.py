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

È realizzata come :class:`ValidazioneIngressi`, sulla classe base comune a
tutte le fasi; :func:`esegui_s0` resta il modo di invocarla da sola, con la
stessa interfaccia di prima. S0 è l'unica fase che legge dati fuori
dall'albero di output, e ne registra un'impronta: letture, tabelle e
riferimento sostituiti sotto lo stesso percorso la rendono non piu' conclusa,
anche a configurazione invariata.
"""

from __future__ import annotations

import csv
import hashlib
import json
import time
from dataclasses import dataclass
from io import StringIO
from pathlib import Path
from typing import Any, ClassVar, Final

from amplicon16s.config.resolve import risolvi
from amplicon16s.config.schema import Config
from amplicon16s.gates.g01_g15 import Contesto, ErroreGate, Violazione
from amplicon16s.gates.registry import EsitoGate, esegui_tutti
from amplicon16s.io_layer.artifacts import AlberoOutput, Fase
from amplicon16s.logging.logger import ottieni
from amplicon16s.metadata.models import Campione, ClasseCampione, Inventario
from amplicon16s.runner.graph import Passo
from amplicon16s.steps.base import PipelineStep, Produzione, StepContext

__all__ = [
    "RisultatoS0",
    "ValidazioneIngressi",
    "esegui_s0",
    "impronta_dati_grezzi",
    "leggi_inventario",
]

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


def leggi_inventario(config: Config) -> Inventario:
    """Rilegge l'inventario dal crosswalk scritto da S0.

    È l'inventario su cui le fasi a valle sono state calcolate: rileggerlo
    dall'artefatto, invece di ricostruirlo dai metadati, fa sì che una ripresa
    usi esattamente quello.
    """
    percorso = AlberoOutput(config.io.out_root).cartella(Fase.INPUT_VALIDATION) / NOME_CROSSWALK
    campioni = []
    with open(percorso, encoding="utf-8", newline="") as file:
        for riga in csv.DictReader(file, delimiter="\t"):
            campioni.append(
                Campione(
                    accession=riga["accession"],
                    nome=riga["sample_name"],
                    classe=ClasseCampione(riga["classe"]),
                    materiale=riga["materiale"],
                    file=Path(config.io.fastq_dir) / riga["file"] if riga["file"] else None,
                    posizione=riga["posizione"] or None,
                    modulo=riga["modulo"] or None,
                    piastra=riga["piastra"] or None,
                    corsa=riga["corsa"] or None,
                )
            )
    return Inventario(tuple(campioni))


def impronta_dati_grezzi(config: Config) -> str:
    """Impronta dei dati letti da S0: letture, tabelle, riferimento.

    Si basa su nome, dimensione e istante di modifica di ogni file, non sul
    contenuto: rileggere 2,5 GB a ogni ripresa costerebbe piu' della
    validazione stessa. Un file sostituito sotto lo stesso nome cambia quasi
    sempre dimensione o istante di modifica; il contenuto del riferimento e'
    comunque verificato da G12 contro tax.ref_md5.
    """

    def descrivi(percorso: Path) -> list[Any]:
        try:
            stato = percorso.stat()
        except OSError:
            return [str(percorso), "assente"]
        return [str(percorso), stato.st_size, stato.st_mtime_ns]

    io = config.io
    file: list[Path] = sorted(Path(io.fastq_dir).glob(io.fastq_glob))
    file += [Path(io.assay_table), Path(io.study_table), Path(config.tax.ref_fasta)]
    if io.batch_table is not None:
        file.append(Path(io.batch_table))

    canonico = json.dumps([descrivi(f) for f in file], separators=(",", ":"))
    return "sha256:" + hashlib.sha256(canonico.encode("utf-8")).hexdigest()


class ValidazioneIngressi(PipelineStep):
    """S0: i quindici gate, e l'esito in ``01_input_validation/``.

    Con ``solleva`` falso restituisce il risultato anche quando un gate
    fallisce, invece di sollevare: la fase risulta allora non superata e non
    viene registrata come conclusa.
    """

    passo: ClassVar[Passo] = Passo.S0

    def __init__(self, *, solleva: bool = True) -> None:
        self.solleva = solleva

    def impronta_dati_esterni(self, config: Config) -> str | None:
        return impronta_dati_grezzi(config)

    def calcola(self, contesto: StepContext) -> Produzione:
        config = contesto.config
        log = contesto.logger
        inizio = time.perf_counter()

        gate = Contesto(config)
        esiti = tuple(esegui_tutti(gate))
        secondi = time.perf_counter() - inizio

        # L'inventario esiste solo se i gate che lo costruiscono sono passati.
        inventario: Inventario | None = None
        if all(e.superato for e in esiti if e.gate in ("G03", "G04", "G05", "G06", "G11")):
            try:
                inventario = gate.inventario
            except Exception:  # noqa: BLE001 - l'inventario e' accessorio al referto
                inventario = None

        albero = contesto.albero
        artefatti = [
            albero.scrivi_testo(
                Fase.INPUT_VALIDATION, NOME_ESITI, _esiti_json(esiti, secondi)
            )
        ]
        if inventario is not None:
            artefatti.append(
                albero.scrivi_testo(
                    Fase.INPUT_VALIDATION, NOME_CROSSWALK, _crosswalk_tsv(inventario)
                )
            )
            artefatti.append(
                albero.scrivi_testo(
                    Fase.INPUT_VALIDATION, NOME_INVENTARIO, _inventario_json(inventario)
                )
            )
        if gate.scansione_gia_fatta:
            artefatti.append(
                albero.scrivi_testo(
                    Fase.INPUT_VALIDATION,
                    NOME_LETTURE,
                    _letture_tsv(gate.scansione, config.qc.head_reads),
                )
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
            artefatti=tuple(a.percorso for a in artefatti),
            secondi=secondi,
        )
        metriche = {
            "secondi": round(secondi, 3),
            "campioni": len(inventario) if inventario else 0,
            "avvisi": len(risultato.avvisi),
        }

        if risultato.superata:
            log.info("S0 superata", extra={"fase": "S0", **metriche})
            return Produzione(tuple(artefatti), metriche, True, risultato)

        falliti = risultato.falliti
        for esito in falliti:
            for violazione in esito.violazioni:
                log.error(
                    violazione.dettaglio,
                    extra={"gate": esito.gate, "codice": violazione.codice, "fase": "S0"},
                )

        if self.solleva:
            primo = falliti[0]
            raise ErroreGate(primo.gate, primo.violazioni)
        return Produzione(tuple(artefatti), metriche, False, risultato)


def esegui_s0(config: Config, *, solleva: bool = True) -> RisultatoS0:
    """Esegue i quindici gate e registra l'esito in ``01_input_validation/``.

    Con ``solleva`` falso restituisce il risultato anche quando un gate
    fallisce, invece di sollevare: serve a chi vuole ispezionare l'esito
    completo, per esempio il comando ``validate``.
    """
    contesto = StepContext(
        risolta=risolvi(config),
        albero=AlberoOutput(config.io.out_root),
        logger=ottieni("s0"),
        inventario=None,
    )
    return ValidazioneIngressi(solleva=solleva).esegui(contesto).dettaglio
