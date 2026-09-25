"""Generatore di scenari sintetici su filesystem per i test della pipeline 16S.

Inquadramento nel Piano Operativo:
    - **Settimane di riferimento**: Trasversale a **W6–W10** (Fasi **F2** e **F3**).
    - **Scopo del modulo**: Fornisce le primitive e la factory ``crea_scenario()``
      per materializzare su disco (tramite la fixture ``tmp_path`` di ``pytest``)
      mini-dataset sintetici conformi o deliberatamente corrotti rispetto allo
      standard ISA-Tab di NASA GeneLab (modellato sul dataset di riferimento
      **OSD-734**). Il dataset reale OSD-734 è integro e supera tutti i 15 gate:
      questo modulo permette di iniettare su file reali (archivi ``.fastq.gz``,
      Assay Table, Study Table, Batch Table e FASTA tassonomico) tutte le
      patologie bioinformatiche e strutturali che i gate ``G01–G14`` devono
      intercettare senza ricorrere a mock in memoria.
    - **Moduli sorgente coperti**:
        * ``src/amplicon16s/metadata/crosswalk.py``
        * ``src/amplicon16s/metadata/controls_map.py``
        * ``src/amplicon16s/io_layer/reads.py``
        * ``src/amplicon16s/gates/g01_g15.py``
        * ``src/amplicon16s/steps/s00_validate.py``
"""

from __future__ import annotations

import csv
import gzip
import hashlib
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from amplicon16s.config.schema import Config, valida

#: Etichette di classificazione dei campioni nella colonna ``ctrl.column``
#: (``Characteristics[Material Type]``), allineate ai valori predefiniti per OSD-734.
BIOLOGICO = "Surface swab"
POSITIVO = "Positive Control"
NEGATIVO = "blank control"

#: Digest SHA-256 formalmente valido per il campo obbligatorio ``run.container``.
CONTAINER = "registro.esempio/amplicon16s@sha256:" + "0" * 64

#: Prefisso nucleotidico che contiene il motivo conservato della regione V4
#: (``TAC[AG].AGG..GC.AGCGTT``) ed è privo del primer forward 515F in testa.
INIZIO_CON_MOTIVO = "TACGGAGGGTGCAAGCGTT"

#: Prefisso nucleotidico che inizia col primer forward 515F (``GTGYCAGCMGCCGCGGTAA``)
#: non rimosso, condizione che il Gate G10 deve bloccare con ``E-S0-10``.
INIZIO_CON_PRIMER = "GTGCCAGCAGCCGCGGTAA"

#: Sequenza omopolimerica priva sia del primer 515F sia del motivo conservato V4,
#: usata per simulare letture prive di segnale biologico 16S o controlli negativi.
INIZIO_MUTO = "CCCCCCCCCCCCCCCCCCC"

#: Lunghezza di lettura standard (151 nt) delle corse Illumina MiSeq di OSD-734.
LUNGHEZZA = 151


def lettura(inizio: str = INIZIO_CON_MOTIVO, lunghezza: int = LUNGHEZZA) -> str:
    """Costruisce una sequenza nucleotidica sintetica di lunghezza prefissata.

    **Obiettivo**: Generare una stringa di basi azotate che inizia con ``inizio``
    ed è completata con code di adenina fino a ``lunghezza`` nucleotidi.

    **Razionale Scientifico/Sistemistico**: Consente ai test dei gate G09 e G10
    di controllare indipendentemente la presenza del primer in 5', la presenza
    del motivo V4 e la lunghezza esatta delle letture rispetto a ``filter.truncLen``.
    """
    return (inizio + "A" * lunghezza)[:lunghezza]


def scrivi_fastq(percorso: Path, sequenze: list[str]) -> None:
    """Materializza su disco un archivio FASTQ compresso GZIP conforme allo standard.

    **Obiettivo**: Scrivere ciascuna sequenza come record FASTQ canonico a 4 righe
    (header ``@``, sequenza, separatore ``+``, qualità Phred ``I`` = Q40) in ``.fastq.gz``.

    **Razionale Scientifico/Sistemistico**: Garantisce che lo scanner in streaming
    ``reads.scansiona_file()`` e il Gate G13 eseguano la vera decompressione ``gzip``
    e il parsing a 4 righe esattamente come avviene sui 960 file di OSD-734.
    """
    with gzip.open(percorso, "wt", encoding="utf-8") as file:
        for indice, sequenza in enumerate(sequenze, start=1):
            file.write(f"@lettura{indice}\n{sequenza}\n+\n{'I' * len(sequenza)}\n")


@dataclass
class Campione:
    """Descrittore dichiarativo di un campione sintetico per la costruzione dello scenario.

    Raccoglie gli attributi bioinformatici (accession ENA, ``Sample Name``, tipo
    di materiale biologico o di controllo, piastra di estrazione, prefisso corsa)
    e le proprietà fisiche del file FASTQ associato.
    """

    accession: str
    nome: str
    materiale: str = BIOLOGICO
    posizione: str = "NOD1D4"
    #: Nome esplicito del file FASTQ; se ``None`` viene generato dall'accession,
    #: mentre se ``""`` omette la creazione del file su disco per testare G06.
    file: str | None = None
    #: Identificativi di piastra e corsa scritti nella ``batch_table`` opzionale (G08).
    piastra: str = "1"
    corsa: str = "corsa_A"
    #: Valore della colonna modulo nella tabella di arricchimento dei lotti.
    modulo_arricchimento: str = "Modulo Uno"
    #: Chiave di join nella ``batch_table``; per default coincide con ``accession``.
    chiave_arricchimento: str | None = None
    #: Sequenza in 5' iniettata nelle letture sintetiche del file FASTQ.
    inizio_letture: str = INIZIO_CON_MOTIVO
    #: Lunghezza in nucleotidi e numerosità dei record scritti nel file FASTQ.
    lunghezza_letture: int = LUNGHEZZA
    numero_letture: int = 40
    #: Payload binario arbitrario per simulare archivi GZIP corrotti o troncati (G13).
    contenuto_grezzo: bytes | None = None


@dataclass
class Scenario:
    """Contenitore immutabile dello scenario materializzato su disco e della sua ``Config``."""

    radice: Path
    config: Config
    campioni: list[Campione] = field(default_factory=list)


def _scrivi_tsv(percorso: Path, intestazione: list[str], righe: list[list[Any]]) -> None:
    """Scrive una tabella TSV con terminatori di riga POSIX (LF) per i metadati ISA-Tab."""
    with open(percorso, "w", encoding="utf-8", newline="") as file:
        scrittore = csv.writer(file, delimiter="\t", lineterminator="\n")
        scrittore.writerow(intestazione)
        scrittore.writerows(righe)


def crea_scenario(
    radice: Path,
    campioni: list[Campione],
    *,
    righe_studio_extra: list[tuple[str, str, str]] | None = None,
    righe_studio_ripetute: list[str] | None = None,
    con_arricchimento: bool = False,
    con_letture: bool = False,
    colonna_chiave_arricchimento: str = "experiment_accession",
    colonna_modulo_arricchimento: str | None = "module",
    file_in_piu: list[str] | None = None,
    sovrascrivi: dict[str, Any] | None = None,
) -> Scenario:
    """Costruisce su filesystem un ambiente sperimentale completo pronto per la validazione.

    **Obiettivo**: Generare nella directory temporanea ``radice`` gli archivi FASTQ,
    l'Assay Table (``assay.txt``), la Study Sample Table (``studio.txt``), l'eventuale
    Batch Table (``lotti.tsv``) e il database FASTA di riferimento con checksum MD5,
    restituendo l'istanza ``Scenario`` con l'oggetto Pydantic ``Config`` già validato.

    **Razionale Scientifico/Sistemistico**: In studi multi-omics NASA GeneLab come
    OSD-734, la Study Table contiene più righe dell'Assay Table 16S (1.072 contro 868)
    perché elenca anche campioni destinati ad altri assay (es. metagenomica shotgun).
    I parametri ``righe_studio_extra`` e ``righe_studio_ripetute`` permettono di
    riprodurre esattamente questa struttura relazionale per collaudare il join
    ristretto (G03) e l'integrità del crosswalk su file fisici reali.
    """
    fastq = radice / "fastq"
    fastq.mkdir(parents=True, exist_ok=True)

    for campione in campioni:
        nome_file = (
            campione.file
            if campione.file is not None
            else f"{campione.accession}_{campione.nome}.fastq.gz"
        )
        if not nome_file:
            continue
        percorso = fastq / nome_file
        if campione.contenuto_grezzo is not None:
            percorso.write_bytes(campione.contenuto_grezzo)
        elif con_letture:
            scrivi_fastq(
                percorso,
                [lettura(campione.inizio_letture, campione.lunghezza_letture)]
                * campione.numero_letture,
            )
        else:
            percorso.write_bytes(b"")
    for nome_file in file_in_piu or []:
        (fastq / nome_file).write_bytes(b"")

    assay = radice / "assay.txt"
    _scrivi_tsv(
        assay,
        ["Sample Name", "Raw Data File"],
        [[c.nome, f"GLDS-000_Amplicon_{c.accession}_raw.fastq.gz"] for c in campioni],
    )

    righe_studio = [
        [c.nome, c.materiale, c.posizione] for c in campioni
    ]
    for nome in righe_studio_ripetute or []:
        originale = next(c for c in campioni if c.nome == nome)
        righe_studio.append([nome, originale.materiale, originale.posizione])
    for nome, materiale, posizione in righe_studio_extra or []:
        righe_studio.append([nome, materiale, posizione])

    studio = radice / "studio.txt"
    _scrivi_tsv(
        studio,
        ["Sample Name", "Characteristics[Material Type]", "Factor Value[Sample Location]"],
        righe_studio,
    )

    arricchimento = None
    if con_arricchimento:
        arricchimento = radice / "lotti.tsv"
        intestazione = [
            colonna_chiave_arricchimento, "extraction_plate_num", "run_prefix"
        ]
        righe = [
            [c.chiave_arricchimento or c.accession, c.piastra, c.corsa]
            for c in campioni
        ]
        if colonna_modulo_arricchimento:
            intestazione.append(colonna_modulo_arricchimento)
            for riga, campione in zip(righe, campioni):
                riga.append(campione.modulo_arricchimento)
        _scrivi_tsv(arricchimento, intestazione, righe)

    # Genera un file FASTA di riferimento minimale e ne calcola il vero digest MD5
    # affinché il Gate G12 passi salvo esplicita manomissione nei test.
    riferimento = radice / "riferimento.fa.gz"
    riferimento.write_bytes(b">seq1\nACGT\n")
    md5_riferimento = hashlib.md5(riferimento.read_bytes()).hexdigest()

    dati: dict[str, Any] = {
        "io": {
            "fastq_dir": str(fastq),
            "assay_table": str(assay),
            "study_table": str(studio),
            "out_root": str(radice / "out"),
            "batch_table": str(arricchimento) if arricchimento else None,
        },
        "tax": {
            "ref_fasta": str(riferimento),
            "ref_md5": md5_riferimento,
            "ref_name": "SILVA",
            "ref_version": "138",
        },
        "run": {"container": CONTAINER, "threads": 1},
    }
    for gruppo, valori in (sovrascrivi or {}).items():
        dati.setdefault(gruppo, {}).update(valori)

    return Scenario(radice=radice, config=valida(dati), campioni=campioni)
