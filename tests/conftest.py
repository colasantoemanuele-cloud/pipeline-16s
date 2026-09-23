"""Costruzione di scenari sintetici per i test del crosswalk.

Il dato reale è corretto, quindi non contiene i casi che i gate devono
intercettare. Questi scenari li costruiscono in piccolo: poche righe, scritte
su file veri in una cartella temporanea, così il percorso di lettura provato
dai test è lo stesso che userà la pipeline e non una simulazione.
"""

from __future__ import annotations

import csv
import gzip
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from amplicon16s.config.schema import Config, valida

#: Etichette di classificazione usate negli scenari, coerenti con i valori
#: predefiniti della configurazione.
BIOLOGICO = "Surface swab"
POSITIVO = "Positive Control"
NEGATIVO = "blank control"

#: Un digest sintatticamente valido per run.container.
CONTAINER = "registro.esempio/amplicon16s@sha256:" + "0" * 64

#: Inizio di lettura che presenta il motivo conservato predefinito
#: ``TAC[AG].AGG..GC.AGCGTT`` e non comincia col primer 515F.
INIZIO_CON_MOTIVO = "TACGGAGGGTGCAAGCGTT"

#: Inizio di lettura che comincia col primer 515F, senza codici degenerati.
INIZIO_CON_PRIMER = "GTGCCAGCAGCCGCGGTAA"

#: Inizio di lettura che non presenta ne' il primer ne' il motivo.
INIZIO_MUTO = "CCCCCCCCCCCCCCCCCCC"

#: Lunghezza predefinita delle letture sintetiche.
LUNGHEZZA = 151


def lettura(inizio: str = INIZIO_CON_MOTIVO, lunghezza: int = LUNGHEZZA) -> str:
    """Una sequenza sintetica che comincia con l'inizio dato."""
    return (inizio + "A" * lunghezza)[:lunghezza]


def scrivi_fastq(percorso: Path, sequenze: list[str]) -> None:
    """Scrive un FASTQ compresso con quattro righe per record."""
    with gzip.open(percorso, "wt", encoding="utf-8") as file:
        for indice, sequenza in enumerate(sequenze, start=1):
            file.write(f"@lettura{indice}\n{sequenza}\n+\n{'I' * len(sequenza)}\n")


@dataclass
class Campione:
    """Un campione dello scenario, con tutto ciò che lo descrive."""

    accession: str
    nome: str
    materiale: str = BIOLOGICO
    posizione: str = "NOD1D4"
    #: Nome del file di letture. Quando è ``None`` viene costruito
    #: dall'accession; quando è ``""`` il file non viene creato.
    file: str | None = None
    #: Piastra e corsa, usate solo se lo scenario ha l'arricchimento.
    piastra: str = "1"
    corsa: str = "corsa_A"
    #: Modulo dichiarato dal file di arricchimento.
    modulo_arricchimento: str = "Modulo Uno"
    #: Chiave con cui il campione compare nel file di arricchimento. Per
    #: difetto e' il suo accession, come pretende la pipeline.
    chiave_arricchimento: str | None = None
    #: Inizio delle letture sintetiche scritte nel file.
    inizio_letture: str = INIZIO_CON_MOTIVO
    #: Lunghezza delle letture e quante scriverne.
    lunghezza_letture: int = LUNGHEZZA
    numero_letture: int = 40
    #: Contenuto grezzo che sostituisce il FASTQ, per provare i file guasti.
    contenuto_grezzo: bytes | None = None


@dataclass
class Scenario:
    """Uno scenario scritto su disco, con la configurazione che lo indica."""

    radice: Path
    config: Config
    campioni: list[Campione] = field(default_factory=list)


def _scrivi_tsv(percorso: Path, intestazione: list[str], righe: list[list[Any]]) -> None:
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
    """Scrive su disco le sorgenti di uno scenario e restituisce la configurazione.

    ``righe_studio_extra`` aggiunge alla tabella campioni di studio righe che
    non appartengono a questo assay: servono a verificare che il join resti
    ristretto. ``righe_studio_ripetute`` duplica il nome di un campione, che è
    il modo in cui la restrizione può rompersi davvero.
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

    riferimento = radice / "riferimento.fa.gz"
    riferimento.write_bytes(b">seq1\nACGT\n")
    import hashlib

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
