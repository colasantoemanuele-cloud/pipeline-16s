#!/usr/bin/env python3
"""Ricostruisce il file di arricchimento del lotto di OSD-734 dalla fonte pubblica.

Il file (``io.batch_table``) associa a ogni campione la piastra di estrazione, il
pozzetto, la corsa di sequenziamento, il modulo della stazione e, per i
controlli positivi, le cellule seminate: senza, la pipeline stima un solo
modello d'errore e nessuna soglia di profondita' per piastra.

**La fonte** sono i metadati Qiita (studio 14542) pubblicati dagli autori dello
studio nel loro repository, ``RodolfoSalido/3DMM``, file
``metadata/2023_03_21_3DMM_metadata_to_update_qiita.txt``, fissato a un commit
e verificato con l'impronta git del contenuto. Il repository non dichiara una
licenza: il file ricostruito non e' nel repository della pipeline, e questo
script ne documenta la provenienza e lo ricostruisce dove serve. Qiita non
espone i metadati senza autenticazione.

**La ricostruzione**, colonna per colonna dalla fonte:

* ``sample_name_ena``: il nome del campione senza il prefisso di Qiita
  (``<artefatto>.14542.``, uno per corsa), cioe' il nome visibile su ENA e nei
  metadati ISA;
* ``qiita_sample_id``: l'identificativo di Qiita (``#SampleID``);
* ``class``: ``SWAB``, ``KATHARO`` o ``BLANK`` da ``sample_type``;
* ``extraction_plate_num``, ``extraction_date``: il numero della piastra e la
  data di estrazione, ricavati da ``sample_plate``
  (``3DMM_<data>_14542_Plate_<n>``), che resta in ``extraction_plate_id``;
* ``well_id``, ``primer_plate``, ``run_prefix``, ``run_date``, ``module``,
  ``katharoseq_cell_count``: copiati come sono;
* ``experiment_accession``, in prima colonna: l'accession ENA dell'esperimento,
  dalla tabella di assay ISA. E' la chiave del join (``meta.batch_key_column``),
  perche' il nome si ripete fra repliche: due campioni sono stati sequenziati
  due volte, una per corsa. Per un nome ripetuto si sceglie l'accession della
  replica che cade nell'intervallo di accession della sua corsa, ricavato dalle
  righe non ambigue; se i candidati non sono esattamente uno lo script si ferma.

Le righe sono ordinate per piastra e pozzetto (come testo).

    python3 dati/osd734/ricostruisci_lotto.py                 # scarica la fonte e scrive
    python3 dati/osd734/ricostruisci_lotto.py --fonte FILE    # da una copia locale
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import re
import sys
import urllib.request
from collections import defaultdict
from pathlib import Path
from typing import Final

QUI = Path(__file__).resolve().parent
ASSAY = QUI / "metadati" / "a_OSD-734_amplicon-sequencing_16s_Illumina MiSeq.txt"
USCITA = QUI / "lotto" / "plate_well_map_960.tsv"
FONTE_LOCALE = QUI / "lotto" / "fonte" / "2023_03_21_3DMM_metadata_to_update_qiita.txt"

#: Il commit del repository degli autori a cui la fonte e' fissata.
COMMIT: Final = "10b69ca2a90c0ed0b6e494fa77fd1a42ea3fa194"
URL_FONTE: Final = (
    f"https://raw.githubusercontent.com/RodolfoSalido/3DMM/{COMMIT}/"
    "metadata/2023_03_21_3DMM_metadata_to_update_qiita.txt"
)
#: L'impronta git (blob SHA-1) del file a quel commit, come la riporta GitHub.
IMPRONTA_GIT: Final = "07f6c07c1b27c94b16d785cef86a8627bec5787a"

PREFISSO_QIITA: Final = re.compile(r"^\d+\.14542\.")
CLASSI: Final = {"surface swab": "SWAB", "control positive": "KATHARO", "control blank": "BLANK"}
PIASTRA: Final = re.compile(r"^3DMM_(\d{8})_14542_Plate_(\d+)$")
ACCESSION: Final = re.compile(r"(?:E|S|D)RX[0-9]{4,}")
COLONNE: Final = (
    "experiment_accession", "sample_name_ena", "qiita_sample_id", "class",
    "extraction_plate_num", "extraction_plate_id", "extraction_date", "well_id",
    "primer_plate", "run_prefix", "run_date", "module", "katharoseq_cell_count",
)


def impronta_git(contenuto: bytes) -> str:
    """L'impronta con cui git identifica il contenuto di un file (blob SHA-1)."""
    return hashlib.sha1(b"blob %d\0" % len(contenuto) + contenuto).hexdigest()


def leggi_fonte(contenuto: bytes) -> list[dict[str, str]]:
    """Le righe della fonte, verificata contro :data:`IMPRONTA_GIT`."""
    if impronta_git(contenuto) != IMPRONTA_GIT:
        raise SystemExit("la fonte non coincide con il file del commit fissato: nulla e' stato scritto")
    testo = contenuto.decode("utf-8")
    return list(csv.DictReader(testo.splitlines(), delimiter="\t"))


def mappa_lotto(righe: list[dict[str, str]]) -> list[dict[str, str]]:
    """Le colonne del file di arricchimento, senza l'accession, nell'ordine finale."""
    mappate = []
    for riga in righe:
        identificativo = riga["#SampleID"]
        piastra = PIASTRA.match(riga["sample_plate"])
        prefisso = PREFISSO_QIITA.match(identificativo)
        if prefisso is None or piastra is None:
            raise SystemExit(f"riga della fonte inattesa: {identificativo!r}, {riga['sample_plate']!r}")
        mappate.append({
            "sample_name_ena": identificativo[prefisso.end():],
            "qiita_sample_id": identificativo,
            "class": CLASSI[riga["sample_type"]],
            "extraction_plate_num": piastra.group(2),
            "extraction_plate_id": riga["sample_plate"],
            "extraction_date": piastra.group(1),
            "well_id": riga["well_id"],
            "primer_plate": riga["primer_plate"],
            "run_prefix": riga["run_prefix"],
            "run_date": riga["run_date"],
            "module": riga["module"],
            "katharoseq_cell_count": riga["katharoseq_cell_count"],
        })
    return sorted(mappate, key=lambda r: (int(r["extraction_plate_num"]), r["well_id"]))


def _numero(accession: str) -> int:
    return int(re.sub(r"^[A-Z]+", "", accession))


def accession_per_nome(assay: Path = ASSAY) -> dict[str, str]:
    """Nome del campione -> accession dell'esperimento, dalla tabella di assay."""
    per_nome = {}
    with open(assay, encoding="utf-8", newline="") as file:
        for riga in csv.DictReader(file, delimiter="\t"):
            trovati = ACCESSION.findall((riga["Raw Data File"] or "").strip('"'))
            if len(trovati) != 1:
                raise SystemExit(f"accession non estraibile dalla tabella di assay: {riga['Raw Data File']!r}")
            per_nome[riga["Sample Name"].strip('"')] = trovati[0]
    return per_nome


def aggiungi_accession(righe: list[dict[str, str]], per_nome: dict[str, str]) -> None:
    """Scrive in ogni riga ``experiment_accession``; per i nomi ripetuti sceglie la
    replica nell'intervallo di accession della corsa.
    """
    quante: dict[str, int] = defaultdict(int)
    for riga in righe:
        quante[riga["sample_name_ena"]] += 1
    per_corsa: dict[str, list[int]] = defaultdict(list)
    for riga in righe:
        nome = riga["sample_name_ena"]
        if quante[nome] == 1 and nome in per_nome:
            per_corsa[riga["run_prefix"]].append(_numero(per_nome[nome]))
    limiti = {corsa: (min(v), max(v)) for corsa, v in per_corsa.items()}
    ordinati = sorted(limiti.values())
    if any(fine >= inizio for (_, fine), (inizio, _) in zip(ordinati, ordinati[1:])):
        raise SystemExit("gli intervalli di accession delle corse si sovrappongono")
    for riga in righe:
        nome = riga["sample_name_ena"]
        if quante[nome] == 1:
            if nome not in per_nome:
                raise SystemExit(f"il campione {nome!r} non compare nella tabella di assay")
            riga["experiment_accession"] = per_nome[nome]
            continue
        inizio, fine = limiti[riga["run_prefix"]]
        candidati = [a for n, a in per_nome.items()
                     if (n == nome or n.startswith(nome + "_")) and inizio <= _numero(a) <= fine]
        if len(candidati) != 1:
            raise SystemExit(f"il campione replicato {nome!r} ha {len(candidati)} candidati nella sua corsa")
        riga["experiment_accession"] = candidati[0]
    if len({r["experiment_accession"] for r in righe}) != len(righe):
        raise SystemExit("gli accession ricavati non sono univoci: nulla e' stato scritto")


def scrivi(righe: list[dict[str, str]], uscita: Path) -> None:
    """Scrive il file di arricchimento, separato da tabulazioni, a capo LF."""
    uscita.parent.mkdir(parents=True, exist_ok=True)
    temporaneo = uscita.with_name(uscita.name + ".parziale")
    with open(temporaneo, "w", encoding="utf-8", newline="") as file:
        scrittore = csv.DictWriter(file, fieldnames=COLONNE, delimiter="\t", lineterminator="\n")
        scrittore.writeheader()
        scrittore.writerows(righe)
    temporaneo.replace(uscita)


def main(argomenti: list[str] | None = None) -> int:
    """Ottiene la fonte, la verifica e scrive il file di arricchimento."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--fonte", type=Path,
                        help="una copia locale della fonte, invece di scaricarla")
    parser.add_argument("--uscita", type=Path, default=USCITA)
    opzioni = parser.parse_args(argomenti)

    if opzioni.fonte is not None:
        contenuto = opzioni.fonte.read_bytes()
    else:
        with urllib.request.urlopen(URL_FONTE, timeout=60) as risposta:
            contenuto = risposta.read()
        FONTE_LOCALE.parent.mkdir(parents=True, exist_ok=True)
        FONTE_LOCALE.write_bytes(contenuto)
    righe = mappa_lotto(leggi_fonte(contenuto))
    aggiungi_accession(righe, accession_per_nome())
    scrivi(righe, opzioni.uscita)
    print(f"{opzioni.uscita}: {len(righe)} campioni")
    return 0


if __name__ == "__main__":
    sys.exit(main())
