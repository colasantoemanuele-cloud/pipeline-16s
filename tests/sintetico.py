"""Un dataset sintetico completo, che attraversa tutte le fasi della pipeline.

Le letture si ricavano dalle sequenze del riferimento tassonomico sintetico
(``tests/fixtures/riferimento_sintetico``): ogni campione e' una miscela di
alcune di quelle sequenze, con errori di sostituzione rari e qualita' variabili,
cosi' che il modello d'errore si stimi, l'inferenza ritrovi le sequenze e la
tassonomia le assegni. Nulla viene da un dataset reale: nomi dei file, delle
colonne e delle etichette sono scelti qui, e sono volutamente diversi da quelli
dei dataset distribuiti con il repository.

I metadati sono una sola tabella dei campioni separata da tabulazioni (non
ISA-Tab), con accanto, a richiesta, un file del lotto. Il generatore e'
deterministico: lo stesso seme da' gli stessi byte.
"""

from __future__ import annotations

import gzip
import hashlib
import random
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

RIFERIMENTO = Path(__file__).parent / "fixtures" / "riferimento_sintetico" / "riferimento.fa.gz"

#: Il primer forward messo in testa alle letture quando richiesto: il 341F di
#: Klindworth et al. 2013, diverso da quello dei dataset del repository.
PRIMER_341F = "CCTACGGGNGGCWGCAG"

#: Basi aggiunte in coda a ogni lettura, oltre la regione del riferimento.
CODA = "GTAGTCCACGCCGTAAACGATG"

#: Nomi delle colonne della tabella dei campioni e del file del lotto.
COLONNA_CAMPIONE = "campione"
COLONNA_FILE = "file_letture"
COLONNA_TIPO = "tipo"
COLONNA_SITO = "sito"
COLONNA_GRUPPO = "trattamento"
COLONNA_CHIAVE_LOTTO = "libreria"
COLONNA_PIASTRA = "piastra_estrazione"
COLONNA_CORSA = "corsa_sequenziamento"

ETICHETTA_BIOLOGICO = "tampone"
ETICHETTA_NEGATIVO = "bianco di estrazione (blank)"


def sequenze_di_riferimento() -> list[tuple[str, str]]:
    """Le coppie (linea tassonomica, sequenza) del riferimento sintetico."""
    righe = gzip.decompress(RIFERIMENTO.read_bytes()).decode().splitlines()
    return [(righe[i][1:], righe[i + 1]) for i in range(0, len(righe), 2)]


@dataclass(frozen=True)
class DatasetSintetico:
    """Dove sta un dataset generato, e come e' fatto."""

    cartella: Path
    letture: Path
    tabella: Path
    lotto: Path | None
    #: Chiave del campione (quella del nome del file) -> tipo dichiarato.
    campioni: dict[str, str]
    #: Chiavi dei campioni il cui file non contiene alcuna lettura.
    vuoti: tuple[str, ...]
    primer: str | None
    lunghezza: int
    con_tipo: bool = True
    nomi: dict[str, str] = field(default_factory=dict)

    @property
    def troncamento(self) -> int:
        """La lunghezza a cui troncare: primer, regione del riferimento, nessuna coda."""
        return len(self.primer or "") + 137


def _lettura(
    sequenza: str, caso: random.Random, raggruppate: bool = False
) -> tuple[str, str]:
    """Una lettura con qualche errore di sostituzione e le sue qualita'."""
    basi = list(sequenza)
    qualita = []
    for i in range(len(basi)):
        # Qualita' fra 28 e 40, piu' basse verso la coda: abbastanza valori
        # distinti perche' il modello d'errore abbia una curva da stimare. Con
        # le qualita' raggruppate, tre soli valori, come nei sequenziatori che
        # le riducono a poche classi.
        if raggruppate:
            q = caso.choice((37, 37, 37, 25)) if i < len(basi) * 0.8 else caso.choice((37, 25, 11))
        else:
            q = caso.choice((40, 40, 40, 38, 37, 36, 34)) if i < len(basi) * 0.8 else caso.choice((36, 34, 32, 30, 28))
        if caso.random() < 10 ** (-q / 10) * 4:
            basi[i] = caso.choice([b for b in "ACGT" if b != basi[i]])
        qualita.append(chr(33 + q))
    return "".join(basi), "".join(qualita)


def _scrivi_fastq(percorso: Path, record: list[tuple[str, str, str]]) -> None:
    """Scrive un FASTQ compresso senza nome ne' istante nell'intestazione gzip."""
    testo = "".join(f"@{nome}\n{sequenza}\n+\n{qualita}\n" for nome, sequenza, qualita in record)
    with open(percorso, "wb") as file:
        with gzip.GzipFile(fileobj=file, mode="wb", mtime=0, filename="") as archivio:
            archivio.write(testo.encode())


def genera(
    cartella: Path,
    *,
    biologici: int = 8,
    negativi: int = 1,
    negativi_vuoti: int = 0,
    letture: int = 1500,
    letture_negativi: int = 60,
    primer: str | None = None,
    con_tipo: bool = True,
    con_lotto: bool = False,
    qualita_raggruppate: bool = False,
    seme: int = 20,
) -> DatasetSintetico:
    """Genera il dataset nella cartella indicata e ne restituisce la descrizione.

    ``negativi_vuoti`` controlli negativi hanno il file senza alcuna lettura.
    Con ``con_tipo`` falso la tabella non ha la colonna del tipo di campione, e
    non si generano controlli. Con ``primer`` ogni lettura comincia con quel
    primer (le basi degenerate si risolvono a caso).
    """
    caso = random.Random(seme)
    riferimento = [sequenza for _, sequenza in sequenze_di_riferimento()]
    cartella.mkdir(parents=True, exist_ok=True)
    cartella_letture = cartella / "letture"
    cartella_letture.mkdir()
    iupac = {"N": "ACGT", "W": "AT", "R": "AG", "Y": "CT", "M": "AC", "K": "GT", "S": "CG"}

    if not con_tipo:
        negativi = negativi_vuoti = 0
    tipi = (
        [ETICHETTA_BIOLOGICO] * biologici
        + [ETICHETTA_NEGATIVO] * (negativi + negativi_vuoti)
    )
    campioni: dict[str, str] = {}
    vuoti: list[str] = []
    righe = []
    righe_lotto = []
    nomi = {}
    for indice, tipo in enumerate(tipi, start=1):
        chiave = f"LIB{indice:03d}"
        nome = f"{'T' if tipo == ETICHETTA_BIOLOGICO else 'B'}{indice:02d}"
        file = f"{chiave}_{nome}.fq.gz"
        campioni[chiave] = tipo
        nomi[chiave] = nome
        e_vuoto = tipo == ETICHETTA_NEGATIVO and indice > biologici + negativi
        record = []
        if not e_vuoto:
            # Sei sequenze per campione, scelte in modo che ogni sequenza
            # compaia in piu' campioni; proporzioni disuguali.
            scelte = [riferimento[(indice * 3 + k * 5) % len(riferimento)] for k in range(6)]
            pesi = [caso.uniform(1, 6) for _ in scelte]
            quante = letture if tipo == ETICHETTA_BIOLOGICO else letture_negativi
            for numero in range(quante):
                testa = "".join(caso.choice(iupac.get(b, b)) for b in (primer or ""))
                sequenza, qualita = _lettura(
                    testa + caso.choices(scelte, pesi)[0] + CODA, caso, qualita_raggruppate)
                record.append((f"{chiave}.{numero + 1}", sequenza, qualita))
        else:
            vuoti.append(chiave)
        _scrivi_fastq(cartella_letture / file, record)
        riga = {
            COLONNA_CAMPIONE: nome, COLONNA_FILE: file,
            COLONNA_SITO: f"sito {1 + indice % 3}",
            COLONNA_GRUPPO: "trattato" if indice % 2 else "non trattato",
        }
        if con_tipo:
            riga[COLONNA_TIPO] = tipo
        righe.append(riga)
        righe_lotto.append({
            COLONNA_CHIAVE_LOTTO: chiave, COLONNA_PIASTRA: f"P{1 + indice % 2}",
            COLONNA_CORSA: "corsa_A",
        })

    def scrivi(percorso: Path, tabella: list[dict[str, str]]) -> None:
        colonne = list(tabella[0])
        percorso.write_text(
            "\t".join(colonne) + "\n"
            + "".join("\t".join(r[c] for c in colonne) + "\n" for r in tabella),
            encoding="utf-8",
        )

    tabella = cartella / "campioni.tsv"
    scrivi(tabella, righe)
    lotto = None
    if con_lotto:
        lotto = cartella / "lotto.tsv"
        scrivi(lotto, righe_lotto)
    return DatasetSintetico(
        cartella, cartella_letture, tabella, lotto, campioni, tuple(vuoti), primer,
        len(primer or "") + 137 + len(CODA), con_tipo, nomi,
    )


def configurazione(dataset: DatasetSintetico, uscita: Path, **sovrascrivi: dict[str, Any]) -> dict[str, Any]:
    """Una configurazione completa per il dataset generato: ogni parametro
    obbligatorio dichiarato, con il riferimento sintetico.
    """
    con_negativi = dataset.con_tipo and ETICHETTA_NEGATIVO in dataset.campioni.values()
    dati: dict[str, dict[str, Any]] = {
        "io": {
            "fastq_dir": str(dataset.letture), "assay_table": str(dataset.tabella),
            "study_table": None, "out_root": str(uscita),
            "batch_table": str(dataset.lotto) if dataset.lotto else None,
            "fastq_glob": "*.fq.gz", "accession_regex": "LIB[0-9]+",
        },
        "meta": {
            "sample_id_column": COLONNA_CAMPIONE, "accession_column": COLONNA_FILE,
            "derive_module": False, "module_regex": None, "module_column": COLONNA_SITO,
            "non_surface_positions": [],
            "batch_key_column": COLONNA_CHIAVE_LOTTO if dataset.lotto else None,
            "batch_module_column": None,
        },
        "filter": {"truncLen": dataset.troncamento, "truncLen_shortfall_warn": 30,
                   "trimLeft": len(dataset.primer or "")},
        "err": {"randomize": False, "batch_column": COLONNA_CORSA if dataset.lotto else None},
        "dada": {"pool": "pseudo"},
        "chimera": {"min_fold_parent_over_abundance": 1.5},
        "tax": {
            "ref_fasta": str(RIFERIMENTO),
            "ref_md5": hashlib.md5(RIFERIMENTO.read_bytes()).hexdigest(),
            "ref_name": "sintetico", "ref_version": "1", "try_rc": False,
        },
        "filt": {"remove_na_phylum": True, "exclude_taxa": []},
        "phylo": {"enabled": False, "max_seqs": 1000},
        "ctrl": {
            "column": COLONNA_TIPO if dataset.con_tipo else None,
            "blank_values": [ETICHETTA_NEGATIVO] if con_negativi else [],
            "positive_values": [],
            "biological_values": [ETICHETTA_BIOLOGICO] if dataset.con_tipo else [],
            "blank_override_column": None, "blank_override_values": [],
            "min_positives": 3, "min_positive_pass_frac": 0.5, "positive_gate": False,
        },
        "katharoseq": {"target_taxon": None, "cell_count_column": None,
                       "collapse_rank": "Genus", "target_sensitivity": 0.8, "min_r2": 0.5},
        "decontam": {"threshold": 0.1, "min_blanks": 2, "mode": "aggregate",
                     "batch_column": COLONNA_PIASTRA if dataset.lotto else None},
        "prev": {"apply": True, "min_fraction": 0.2, "min_count": 1},
        "qc": {
            "head_reads": 500, "primer_sequence": dataset.primer or "GGGGGGGGGGGGGGGGGGGG",
            "max_primer_hit_frac": 0.1, "conserved_motif": None, "min_motif_frac": 0.2,
            "max_frac_short_reads": 0.02, "max_zeroed_samples": 0,
            "max_frac_lost_filter": 0.5, "max_asv_count": 10000,
            "warn_frac_chimeric": 0.2, "stop_frac_chimeric": 0.6, "min_frac_phylum": 0.5,
            "min_reads_mode": "none", "max_frac_contaminant": 0.5,
            "min_reads_final": 100, "min_frac_reads_retained": 0.2,
        },
        "out": {"study_columns": [COLONNA_GRUPPO], "batch_columns": []},
        "run": {"threads": 2},
    }
    for gruppo, valori in sovrascrivi.items():
        dati.setdefault(gruppo, {}).update(valori)
    return dati
