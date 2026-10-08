"""Generatore di scenari sintetici su filesystem per i test della pipeline 16S.

Inquadramento nel Piano Operativo:
    - **Settimane di riferimento**: Trasversale a **W6-W10** (Fasi **F2** e **F3**),
      reso indipendente dal dataset di riferimento nella **W30**.
    - **Scopo del modulo**: Fornisce le primitive e la factory ``crea_scenario()``
      per materializzare su disco (tramite la fixture ``tmp_path`` di ``pytest``)
      mini-dataset sintetici, conformi o deliberatamente corrotti: archivi
      ``.fastq.gz``, tabella di assay, tabella campioni di studio, file del lotto
      e riferimento tassonomico. Permette di iniettare su file veri le patologie
      che i gate ``G01-G14`` devono intercettare, senza mock in memoria.
    - **Il formato e' un parametro**: nomi delle colonne, etichette delle classi,
      nomi dei file e valori dei parametri che descrivono il dataset stanno in
      :class:`Formato`, con valori sintetici propri. Nessuna etichetta, colonna
      o nome del dataset di riferimento e' scritto qui: uno scenario dichiara
      nella propria configurazione ogni parametro che descrive il dataset, e
      non ne eredita alcuno dai predefiniti.
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

from amplicon16s.config import defaults
from amplicon16s.config.schema import Config, valida


@dataclass(frozen=True)
class Formato:
    """Il formato di un dataset sintetico: come si chiamano colonne, classi e
    file, e i valori dei parametri che descrivono il dataset.

    I valori predefiniti sono inventati per i test. Un test che esercita un
    formato diverso ne costruisce un altro con ``dataclasses.replace`` e lo
    passa a :func:`crea_scenario`.
    """

    # --- tabella di assay e tabella campioni di studio ----------------------
    colonna_campione: str = "campione"
    colonna_file: str = "file delle letture"
    colonna_classe: str = "tipo di materiale"
    colonna_posizione: str = "posizione di prelievo"
    colonna_cellule: str = "cellule seminate"
    # --- etichette delle classi dei campioni ---------------------------------
    biologico: str = "tampone di superficie"
    positivo: str = "controllo positivo"
    negativo: str = "bianco di estrazione"
    #: Posizioni che non sono superfici, e posizione che riclassifica un
    #: campione dichiarato biologico come controllo negativo.
    non_superfici: tuple[str, ...] = ("campione d'aria", "tampone mai aperto")
    riclassificati: tuple[str, ...] = ("tampone mai aperto",)
    senza_posizione: str = "non applicabile"
    # --- nomi dei file --------------------------------------------------------
    #: Il file delle letture di un campione e il nome riportato nell'assay.
    modello_file: str = "{accession}_{nome}.fastq.gz"
    modello_file_assay: str = "LETTURE_{accession}_grezze.fastq.gz"
    #: La chiave del campione: un accession di esperimento degli archivi INSDC.
    regex_accession: str = r"[ESD]RX[0-9]{4,}"
    #: Il modulo: le prime tre lettere e la cifra della posizione.
    regex_modulo: str = r"^([A-Z]{3}[0-9])"
    # --- file del lotto -------------------------------------------------------
    colonna_chiave_lotto: str = "accession"
    colonna_piastra: str = "piastra"
    colonna_corsa: str = "corsa"
    colonna_modulo_lotto: str = "modulo"
    # --- l'esperimento --------------------------------------------------------
    troncamento: int = 140
    primer: str = "GTGYCAGCMGCCGCGGTAA"
    motivo: str = r"TAC[AG].AGG..GC.AGCGTT"
    taxon_atteso: str = "Genere_bersaglio"
    riferimento: tuple[str, str] = ("riferimento di prova", "1")
    # --- scelte di analisi e soglie di qualita' ------------------------------
    #: Dichiarate dallo scenario come ogni altro parametro obbligatorio: la
    #: pipeline non ha per esse alcun predefinito. Sono valori di prova, scelti
    #: perche' le letture sintetiche attraversino i controlli.
    scelte: tuple[tuple[str, Any], ...] = (
        ("io.fastq_glob", "*.fastq.gz"),
        ("meta.derive_module", True),
        ("filter.truncLen_shortfall_warn", 10),
        ("filter.trimLeft", 0),
        ("err.randomize", True),
        ("dada.pool", "pseudo"),
        ("chimera.min_fold_parent_over_abundance", 2.0),
        ("tax.try_rc", True),
        ("tax.min_boot", 60),
        ("filt.remove_na_phylum", True),
        ("filt.exclude_taxa", ("Organello_di_prova", "Dominio_escluso")),
        ("phylo.enabled", False),
        ("phylo.max_seqs", 5000),
        ("prev.apply", True),
        ("prev.min_fraction", 0.02),
        ("prev.min_count", 2),
        ("ctrl.min_positives", 3),
        ("ctrl.min_positive_pass_frac", 0.75),
        ("ctrl.positive_gate", False),
        ("katharoseq.collapse_rank", "Genus"),
        ("katharoseq.target_sensitivity", 0.90),
        ("katharoseq.min_r2", 0.80),
        ("decontam.threshold", 0.4),
        ("decontam.min_blanks", 5),
        ("decontam.mode", "aggregate"),
        ("qc.head_reads", 10000),
        ("qc.max_primer_hit_frac", 0.05),
        ("qc.min_motif_frac", 0.3),
        ("qc.max_frac_short_reads", 0.05),
        ("qc.max_zeroed_samples", 0),
        ("qc.max_frac_lost_filter", 0.30),
        ("qc.max_asv_count", 300000),
        ("qc.warn_frac_chimeric", 0.25),
        ("qc.stop_frac_chimeric", 0.50),
        ("qc.min_frac_phylum", 0.80),
        ("qc.min_reads_mode", "katharoseq_if_available"),
        ("qc.max_frac_contaminant", 0.40),
        ("qc.min_reads_final", 1000),
        ("qc.min_frac_reads_retained", 0.40),
    )


#: Il formato degli scenari che non ne indicano uno.
FORMATO = Formato()

#: Etichette di classificazione dei campioni nella colonna ``ctrl.column`` degli
#: scenari con il formato predefinito.
BIOLOGICO = FORMATO.biologico
POSITIVO = FORMATO.positivo
NEGATIVO = FORMATO.negativo

#: Digest SHA-256 formalmente valido per il campo ``run.container``.
CONTAINER = "registro.esempio/amplicon16s@sha256:" + "0" * 64

#: Prefisso nucleotidico che contiene il motivo conservato del formato
#: predefinito ed è privo del primer in testa.
INIZIO_CON_MOTIVO = "TACGGAGGGTGCAAGCGTT"

#: Prefisso nucleotidico che inizia col primer del formato predefinito, non
#: rimosso: la condizione che il Gate G10 deve bloccare con ``E-S0-10``.
INIZIO_CON_PRIMER = "GTGCCAGCAGCCGCGGTAA"

#: Sequenza omopolimerica priva sia del primer sia del motivo conservato,
#: usata per simulare letture prive del segnale atteso o controlli negativi.
INIZIO_MUTO = "CCCCCCCCCCCCCCCCCCC"

#: Lunghezza delle letture sintetiche, maggiore del troncamento del formato.
LUNGHEZZA = 151


def lettura(inizio: str = INIZIO_CON_MOTIVO, lunghezza: int = LUNGHEZZA) -> str:
    """Costruisce una sequenza nucleotidica sintetica di lunghezza prefissata.

    **Obiettivo**: Generare una stringa di basi azotate che inizia con ``inizio``
    ed è completata con code di adenina fino a ``lunghezza`` nucleotidi.

    **Razionale scientifico e sistemistico**: Consente ai test dei gate G09 e G10
    di controllare indipendentemente la presenza del primer in 5', la presenza
    del motivo conservato e la lunghezza esatta delle letture rispetto a
    ``filter.truncLen``.
    """
    return (inizio + "A" * lunghezza)[:lunghezza]


def scrivi_fastq(percorso: Path, sequenze: list[str]) -> None:
    """Materializza su disco un archivio FASTQ compresso GZIP conforme allo standard.

    **Obiettivo**: Scrivere ciascuna sequenza come record FASTQ canonico a 4 righe
    (header ``@``, sequenza, separatore ``+``, qualità Phred ``I`` = Q40) in ``.fastq.gz``.

    **Razionale scientifico e sistemistico**: Garantisce che lo scanner in streaming
    ``reads.scansiona_file()`` e il Gate G13 eseguano la vera decompressione ``gzip``
    e il parsing a 4 righe, come su file reali.
    """
    with gzip.open(percorso, "wt", encoding="utf-8") as file:
        for indice, sequenza in enumerate(sequenze, start=1):
            file.write(f"@lettura{indice}\n{sequenza}\n+\n{'I' * len(sequenza)}\n")


@dataclass
class Campione:
    """Descrittore dichiarativo di un campione sintetico per la costruzione dello scenario.

    Raccoglie gli attributi del campione (accession, nome, classe, posizione,
    piastra e corsa) e le proprietà fisiche del file FASTQ associato.
    """

    accession: str
    nome: str
    materiale: str = BIOLOGICO
    #: La posizione di prelievo: le prime tre lettere e la cifra sono il modulo.
    posizione: str = "MOD1A2"
    #: Nome esplicito del file FASTQ; se ``None`` viene generato dal modello del
    #: formato, mentre se ``""`` omette la creazione del file su disco (G06).
    file: str | None = None
    #: Identificativi di piastra e corsa scritti nel file del lotto (G08).
    piastra: str = "1"
    corsa: str = "corsa_A"
    #: Valore della colonna del modulo nel file del lotto.
    modulo_arricchimento: str = "Modulo Uno"
    #: Chiave di join nel file del lotto; per default coincide con ``accession``.
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
    """Lo scenario materializzato su disco, la sua ``Config`` e il suo formato."""

    radice: Path
    config: Config
    campioni: list[Campione] = field(default_factory=list)
    formato: Formato = FORMATO


def _scrivi_tsv(percorso: Path, intestazione: list[str], righe: list[list[Any]]) -> None:
    """Scrive una tabella TSV con terminatori di riga POSIX (LF)."""
    with open(percorso, "w", encoding="utf-8", newline="") as file:
        scrittore = csv.writer(file, delimiter="\t", lineterminator="\n")
        scrittore.writerow(intestazione)
        scrittore.writerows(righe)


def parametri_del_formato(formato: Formato = FORMATO, *, con_lotto: bool = True) -> dict[str, Any]:
    """I parametri che descrivono il dataset, per gruppo, come li dichiara uno
    scenario con il formato dato: tutti gli obbligatori, comprese le scelte di
    analisi e le soglie. Senza file del lotto le sue colonne sono nulle.
    """
    dati: dict[str, Any] = {
        "io": {"accession_regex": formato.regex_accession},
        "meta": {
            "sample_id_column": formato.colonna_campione,
            "accession_column": formato.colonna_file,
            "module_regex": formato.regex_modulo,
            "module_column": formato.colonna_posizione,
            "non_surface_positions": list(formato.non_superfici),
            "batch_key_column": formato.colonna_chiave_lotto if con_lotto else None,
            "batch_module_column": formato.colonna_modulo_lotto if con_lotto else None,
        },
        "err": {"batch_column": formato.colonna_corsa if con_lotto else None},
        "decontam": {"batch_column": formato.colonna_piastra if con_lotto else None},
        "ctrl": {
            "column": formato.colonna_classe,
            "blank_values": [formato.negativo],
            "positive_values": [formato.positivo],
            "biological_values": [formato.biologico],
            "blank_override_column": formato.colonna_posizione,
            "blank_override_values": list(formato.riclassificati),
        },
        "katharoseq": {
            "cell_count_column": formato.colonna_cellule,
            "target_taxon": formato.taxon_atteso,
        },
        # Delle colonne da portare nell'oggetto resta la sola colonna delle
        # cellule, che la calibrazione legge dall'oggetto.
        "out": {"study_columns": [formato.colonna_cellule], "batch_columns": []},
        "filter": {"truncLen": formato.troncamento},
        "qc": {"primer_sequence": formato.primer, "conserved_motif": formato.motivo},
        "tax": {"ref_name": formato.riferimento[0], "ref_version": formato.riferimento[1]},
    }
    for gruppo, valori in scelte_del_formato(formato).items():
        dati.setdefault(gruppo, {}).update(valori)
    dichiarati = {f"{g}.{n}" for g, valori in dati.items() for n in valori}
    assert set(defaults.OBBLIGATORI) <= dichiarati
    return dati


def scelte_del_formato(formato: Formato = FORMATO) -> dict[str, dict[str, Any]]:
    """Le scelte di analisi e le soglie di qualita' del formato, per gruppo: per
    i test che scrivono a mano la descrizione di un dataset e prendono dal
    formato il resto dei parametri obbligatori.
    """
    dati: dict[str, dict[str, Any]] = {}
    for chiave, valore in formato.scelte:
        gruppo, nome = chiave.split(".")
        dati.setdefault(gruppo, {})[nome] = list(valore) if isinstance(valore, tuple) else valore
    return dati


def dichiarazione_minima(config: Config) -> dict[str, Any]:
    """La configurazione piu' breve che la pipeline accetta: i percorsi, il
    riferimento e i parametri obbligatori, con i valori che hanno in
    ``config``. Tutto il resto resta al predefinito.

    Serve ai test che scrivono un file di configurazione e lo passano alla riga
    di comando: senza i parametri obbligatori G15 lo respingerebbe.
    """
    completa = config.model_dump(mode="json")
    dati: dict[str, Any] = {
        "io": {n: completa["io"][n] for n in ("fastq_dir", "assay_table", "study_table",
                                              "out_root", "batch_table")},
        "tax": {n: completa["tax"][n] for n in ("ref_fasta", "ref_md5")},
        "run": {"container": completa["run"]["container"], "threads": 1},
    }
    for chiave in defaults.OBBLIGATORI:
        gruppo, nome = chiave.split(".")
        dati.setdefault(gruppo, {})[nome] = completa[gruppo][nome]
    return dati


def configurazione_di_prova(formato: Formato = FORMATO) -> dict[str, Any]:
    """Una configurazione completa e valida, come dizionario modificabile: i
    parametri del formato, percorsi che non esistono, e ogni altro parametro
    dello schema scritto con il suo valore predefinito.

    Serve ai test dello schema e della risoluzione, che cambiano un parametro
    alla volta: nessun file viene letto.
    """
    dati = parametri_del_formato(formato)
    dati["io"].update(
        fastq_dir="/dati/letture", assay_table="/dati/assay.txt",
        study_table="/dati/studio.txt", batch_table="/dati/lotto.tsv", out_root="/uscita",
    )
    dati["tax"].update(ref_fasta="/dati/riferimento.fa.gz", ref_md5="0" * 32)
    dati["out"]["batch_columns"] = ["data di estrazione"]
    completa = valida(dati).model_dump(mode="json")
    # I parametri senza valore restano fuori, come in un file di configurazione.
    for gruppo in ("run", "meta"):
        completa[gruppo] = {k: v for k, v in completa[gruppo].items() if v is not None}
    return completa


def crea_scenario(
    radice: Path,
    campioni: list[Campione],
    *,
    formato: Formato = FORMATO,
    righe_studio_extra: list[tuple[str, str, str]] | None = None,
    righe_studio_ripetute: list[str] | None = None,
    con_arricchimento: bool = False,
    con_letture: bool = False,
    colonna_chiave_arricchimento: str | None = None,
    colonna_modulo_arricchimento: str | None = ...,  # type: ignore[assignment]
    file_in_piu: list[str] | None = None,
    sovrascrivi: dict[str, Any] | None = None,
) -> Scenario:
    """Costruisce su filesystem un ambiente sperimentale completo pronto per la validazione.

    **Obiettivo**: Generare nella directory temporanea ``radice`` gli archivi FASTQ,
    la tabella di assay (``assay.txt``), la tabella campioni di studio
    (``studio.txt``), l'eventuale file del lotto (``lotti.tsv``) e il riferimento
    tassonomico con il suo MD5, con i nomi di colonne, etichette e file di
    ``formato``, restituendo lo ``Scenario`` con la ``Config`` già validata.

    **Razionale scientifico e sistemistico**: In uno studio con piu' assay la
    tabella di studio ha piu' righe di quella di assay, perche' elenca anche i
    campioni destinati agli altri. ``righe_studio_extra`` e
    ``righe_studio_ripetute`` riproducono questa struttura per collaudare il
    join ristretto (G03) e l'integrità del crosswalk su file fisici.

    ``colonna_chiave_arricchimento`` e ``colonna_modulo_arricchimento`` cambiano
    le sole intestazioni del file del lotto (``None`` per il modulo lo omette),
    lasciando alla configurazione i nomi del formato: servono ai test che
    provano un file del lotto diverso da quello dichiarato.
    """
    fastq = radice / "fastq"
    fastq.mkdir(parents=True, exist_ok=True)

    for campione in campioni:
        nome_file = (
            campione.file
            if campione.file is not None
            else formato.modello_file.format(accession=campione.accession, nome=campione.nome)
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
        [formato.colonna_campione, formato.colonna_file],
        [[c.nome, formato.modello_file_assay.format(accession=c.accession, nome=c.nome)]
         for c in campioni],
    )

    righe_studio = [
        [c.nome, c.materiale, c.posizione] for c in campioni
    ]
    for nome in righe_studio_ripetute or []:
        originale = next(c for c in campioni if c.nome == nome)
        righe_studio.append([nome, originale.materiale, originale.posizione])
    for nome, materiale, posizione in righe_studio_extra or []:
        righe_studio.append([nome, materiale, posizione])

    # La colonna delle cellule dei controlli positivi sta nella tabella di
    # studio, vuota: lo scenario non ha un file del lotto che la porti, e la
    # configurazione che dichiara controlli positivi deve dichiararla.
    studio = radice / "studio.txt"
    _scrivi_tsv(
        studio,
        [formato.colonna_campione, formato.colonna_classe, formato.colonna_posizione,
         formato.colonna_cellule],
        [[*riga, ""] for riga in righe_studio],
    )

    arricchimento = None
    if con_arricchimento:
        arricchimento = radice / "lotti.tsv"
        intestazione = [
            colonna_chiave_arricchimento or formato.colonna_chiave_lotto,
            formato.colonna_piastra, formato.colonna_corsa,
        ]
        righe = [
            [c.chiave_arricchimento or c.accession, c.piastra, c.corsa]
            for c in campioni
        ]
        colonna_modulo = (
            formato.colonna_modulo_lotto if colonna_modulo_arricchimento is ...
            else colonna_modulo_arricchimento
        )
        if colonna_modulo:
            intestazione.append(colonna_modulo)
            for riga, campione in zip(righe, campioni):
                riga.append(campione.modulo_arricchimento)
        _scrivi_tsv(arricchimento, intestazione, righe)

    # Genera un file FASTA di riferimento minimale e ne calcola il vero digest MD5
    # affinché il Gate G12 passi salvo esplicita manomissione nei test.
    riferimento = radice / "riferimento.fa.gz"
    riferimento.write_bytes(b">seq1\nACGT\n")
    md5_riferimento = hashlib.md5(riferimento.read_bytes()).hexdigest()

    dati = parametri_del_formato(formato, con_lotto=arricchimento is not None)
    dati["io"].update(
        fastq_dir=str(fastq), assay_table=str(assay), study_table=str(studio),
        out_root=str(radice / "out"),
        batch_table=str(arricchimento) if arricchimento else None,
    )
    dati["tax"].update(ref_fasta=str(riferimento), ref_md5=md5_riferimento)
    dati["run"] = {"container": CONTAINER, "threads": 1}
    for gruppo, valori in (sovrascrivi or {}).items():
        dati.setdefault(gruppo, {}).update(valori)

    return Scenario(radice=radice, config=valida(dati), campioni=campioni, formato=formato)


# --------------------------------------------------------------------------- #
# Aiuti che compongono funzioni di produzione                                  #
# --------------------------------------------------------------------------- #
# La pipeline non ha bisogno di queste composizioni (la riga di comando e
# l'esecutore chiamano le stesse funzioni una alla volta): servono ai test, e
# per questo stanno qui e non nel pacchetto.


def esegui_g15(dati: dict[str, Any]):
    """G15 su una configurazione non ancora validata: lo schema e poi i
    controlli di coerenza, come li eseguono la riga di comando e l'esecutore.

    Restituisce la configurazione risolta; solleva ``ErroreGate`` con il codice
    del controllo fallito.
    """
    from amplicon16s.config.resolve import risolvi
    from amplicon16s.config.schema import ErroreConfigurazione
    from amplicon16s.gates.g01_g15 import ErroreGate, controlla_coerenza, rifiuto_di_g15

    try:
        config = valida(dict(dati))
    except ErroreConfigurazione as errore:
        raise rifiuto_di_g15(errore) from errore
    risolta = risolvi(config)
    violazioni = controlla_coerenza(risolta)
    if violazioni:
        raise ErroreGate("G15", violazioni)
    return risolta


#: I gate sui metadati, nell'ordine del registro.
GATE_METADATI = ("G04", "G05", "G06", "G03", "G11")


def esegui_gate_metadati(config: Config):
    """L'inventario dei campioni, dopo i cinque gate sui metadati eseguiti dal
    registro; solleva ``ErroreGate`` al primo che trova violazioni.
    """
    from amplicon16s.gates.g01_g15 import Contesto, ErroreGate
    from amplicon16s.gates.registry import esegui_gate

    contesto = Contesto(config)
    for nome in GATE_METADATI:
        esito = esegui_gate(nome, contesto)
        if esito.violazioni:
            raise ErroreGate(nome, esito.violazioni)
    return contesto.inventario


def nomi_dei_gate() -> tuple[str, ...]:
    """I nomi dei gate, nell'ordine di esecuzione del registro."""
    from amplicon16s.gates.registry import REGISTRO

    return tuple(voce.nome for voce in REGISTRO)


def chiavi_schema() -> frozenset[str]:
    """Tutte le chiavi dello schema nella forma ``gruppo.parametro``."""
    return frozenset(
        f"{gruppo}.{campo}" for gruppo, modello in Config.model_fields.items()
        for campo in modello.annotation.model_fields
    )


def codici_con_retry() -> frozenset[str]:
    """I codici che il catalogo ammette al retry automatico."""
    from amplicon16s.errors.catalog import CATALOGO

    return frozenset(codice for codice, voce in CATALOGO.items() if voce.ammette_retry)


def fase_completa(albero, fase, attesi=None) -> bool:
    """Se la cartella di una fase ha artefatti registrati, tutti integri, e fra
    essi quelli attesi.
    """
    registrati = albero.manifesto(fase)
    if not registrati:
        return False
    if attesi is not None and not set(attesi) <= set(registrati):
        return False
    return not albero.non_integri(fase)


def discendenti(passo):
    """Le fasi del grafo che dipendono, anche indirettamente, da quella data."""
    from amplicon16s.runner.graph import GRAFO, Passo

    # Fino al punto fisso: l'ordine dei numeri non e' quello delle dipendenze
    # (la filogenesi, S9, dipende da S13).
    raggiunti = {passo}
    while True:
        nuovi = {p for p in Passo if raggiunti & set(GRAFO.nodo(p).dipendenze)} - raggiunti
        if not nuovi:
            break
        raggiunti |= nuovi
    raggiunti.discard(passo)
    return tuple(p for p in Passo if p in raggiunti)


def rscript_con_limite(cartella: Path, byte: int) -> Path:
    """Un interprete R con la memoria virtuale limitata, per provare la memoria
    esaurita senza toccare quella della macchina.

    E' un involucro di shell da passare al ponte come interprete: applica il
    limite al solo processo R e riduce a uno i thread dell'algebra lineare,
    perche' con OpenBLAS multithread R intercetta l'allocazione fallita ma poi
    resta bloccato in uscita.
    """
    import shutil
    import stat

    involucro = cartella / "Rscript-limitato"
    involucro.write_text(
        "#!/bin/sh\n"
        f"ulimit -v {byte // 1024}\n"
        "export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1\n"
        f'exec {shutil.which("Rscript")} "$@"\n',
        encoding="utf-8",
    )
    involucro.chmod(involucro.stat().st_mode | stat.S_IXUSR)
    return involucro


# --------------------------------------------------------------------------- #
# Esecuzione di base sulla versione ridotta del sottoinsieme di prova          #
# --------------------------------------------------------------------------- #

import pytest  # noqa: E402

#: Il gruppo pytest-xdist dei test sui dati reali.
GRUPPO_DATI_REALI = "dati_reali"
#: I marcatori dei test che eseguono sul dataset completo: quelli di
#: correttezza (``dati_reali``) e il confronto con i checksum pubblicati
#: (``riferimento``), che condivide con i primi la catena e quindi il gruppo.
MARCATORI_SUL_DATASET_COMPLETO = ("dati_reali", "riferimento")

#: Quanto un processo attende, al massimo, una catena condivisa che un altro
#: sta calcolando: due ore, piu' della catena completa sul dataset di
#: riferimento su una macchina occupata. Si cambia con la variabile d'ambiente.
VARIABILE_ATTESA = "AMPLICON16S_ATTESA_CONDIVISA_S"
ATTESA_MASSIMA_S = int(__import__("os").environ.get(VARIABILE_ATTESA, "7200"))


@pytest.hookimpl(tryfirst=True)
def pytest_collection_modifyitems(config, items):
    """Mette tutti i test sul dataset completo (``dati_reali`` e ``riferimento``)
    in un solo gruppo di pytest-xdist.

    Con ``--dist loadgroup`` (addopts in ``pyproject.toml``) i test di un gruppo
    girano tutti sullo stesso processo, uno dopo l'altro, e gli altri restano
    distribuiti uno per uno sui quattro processi, come con ``load``. Sul
    dataset completo i test reali lanciano catene con decine di processi R:
    distribuiti, giravano in parallelo fra loro e saturavano la memoria. Il
    gruppo si aggiunge prima che pytest-xdist legga i gruppi (tryfirst).
    """
    for item in items:
        if any(item.get_closest_marker(m) is not None for m in MARCATORI_SUL_DATASET_COMPLETO):
            item.add_marker(pytest.mark.xdist_group(GRUPPO_DATI_REALI))


def _condivisa(tmp_path_factory, nome, calcola):
    """La coppia (esecuzione, esito) di una catena condivisa, calcolata una volta
    sola per tutta la suite, anche con piu' processi (pytest-xdist).

    ``calcola`` riceve la cartella in cui calcolare e restituisce la coppia, o
    ``None``. In sequenza la si calcola qui. Con piu' processi la calcola il
    primo che la chiede, in una cartella comune a tutti; gli altri ne attendono
    la conclusione e ricostruiscono l'esecuzione dalla configurazione, con
    l'esito scritto dal primo: senza, ogni processo rifarebbe la stessa catena
    per i moduli che gli toccano (sul dataset completo, tre quarti d'ora per
    processo). La catena resta in sola lettura, come in sequenza: i test che
    devono modificarla ne usano una copia, e restano indipendenti fra loro.
    """
    import os
    import pickle
    import time

    from amplicon16s.runner.project import ProjectRun

    if os.environ.get("PYTEST_XDIST_WORKER") is None:
        return calcola(tmp_path_factory.mktemp(nome))
    comune = tmp_path_factory.getbasetemp().parent / "condivise"
    comune.mkdir(exist_ok=True)
    pronta = comune / f"{nome}.pickle"
    guasta = comune / f"{nome}.guasta"
    try:
        os.close(os.open(comune / f"{nome}.blocco", os.O_CREAT | os.O_EXCL | os.O_WRONLY))
    except FileExistsError:
        # L'attesa ha un limite: se il processo che calcola muore senza lasciare
        # traccia (ucciso dal sistema, per esempio per memoria), gli altri non
        # devono restare fermi per sempre.
        scadenza = time.monotonic() + ATTESA_MASSIMA_S
        while not pronta.exists():
            if guasta.exists():
                raise RuntimeError(
                    f"la catena condivisa {nome} e' fallita in un altro processo: "
                    f"{guasta.read_text(encoding='utf-8')}"
                ) from None
            if time.monotonic() > scadenza:
                raise RuntimeError(
                    f"la catena condivisa {nome} non e' pronta dopo {ATTESA_MASSIMA_S} s: il "
                    "processo che la calcola e' probabilmente morto senza dichiararlo. "
                    f"Il limite si cambia con la variabile {VARIABILE_ATTESA}"
                ) from None
            time.sleep(0.5)
        valore = pickle.loads(pronta.read_bytes())
        return None if valore is None else (ProjectRun(valore[0]), valore[1])
    try:
        risultato = calcola(comune / nome)
    except BaseException as guasto:
        guasta.write_text(repr(guasto), encoding="utf-8")
        raise
    provvisoria = comune / f"{nome}.provvisoria"
    provvisoria.write_bytes(
        pickle.dumps(None if risultato is None else (risultato[0].config, risultato[1]))
    )
    provvisoria.replace(pronta)
    return risultato


@pytest.fixture(scope="session")
def ridotta_calcolata(tmp_path_factory):
    """S0-S3 sulla versione ridotta, con S1, una volta sola per l'intera sessione.

    Serve ai test che leggono gli artefatti delle fasi: rifarla in ogni modulo
    costerebbe mezzo minuto a modulo senza verificare nulla di piu'. S1 non e'
    fra gli antenati di S3, e ``fino_a=S3`` non la esegue: la si chiede a parte,
    perche' i test di S1, del filtro e del tracciamento ne leggono gli
    artefatti. L'esito restituito e' quello dell'esecuzione fino a S3. E' in
    sola lettura: un test che deve modificare l'albero ne usa una copia.
    ``None`` dove mancano R e Bioconductor; le fixture dei moduli decidono se
    saltare.
    """
    from sottoinsieme import config_ridotta, motivo_pacchetti_r_assenti

    from amplicon16s.runner.executor import Esecutore
    from amplicon16s.runner.graph import Passo
    from amplicon16s.runner.project import ProjectRun

    if motivo_pacchetti_r_assenti("dada2", "ggplot2", "ShortRead", "Biostrings", "jsonlite"):
        return None

    def calcola(cartella):
        """S0-S3 e poi S1 sulla versione ridotta in ``cartella``: la coppia
        (esecuzione, esito fino a S3) che :func:`_condivisa` mette in comune."""
        run = ProjectRun(config_ridotta(cartella))
        esito = Esecutore(run, fino_a=Passo.S3).esegui()
        Esecutore(run, fino_a=Passo.S1).esegui()
        return run, esito

    return _condivisa(tmp_path_factory, "ridotta", calcola)


def _sopra(base, tmp_path_factory, nome, fino_a):
    """Una catena condivisa che prosegue fino a ``fino_a`` su una copia di ``base``."""
    from amplicon16s.runner.executor import Esecutore

    if base is None:
        return None

    def calcola(cartella):
        """Copia ``base`` in ``cartella`` e prosegue fino a ``fino_a``: la coppia
        (esecuzione, esito) che :func:`_condivisa` mette in comune."""
        run = copia_esecuzione(base, cartella)
        return run, Esecutore(run, fino_a=fino_a).esegui()

    return _condivisa(tmp_path_factory, nome, calcola)


def copia_esecuzione(base, cartella, **sovrascrivi):
    """L'albero di un'esecuzione, copiato in ``cartella``, con parametri cambiati.

    ``base`` e' la coppia (esecuzione, esito) di una fixture condivisa;
    ``sovrascrivi`` cambia chiavi per gruppo, come ``run={"batch_size": 5}``.
    """
    import shutil

    from amplicon16s.config.schema import valida
    from amplicon16s.runner.project import ProjectRun

    run, _ = base
    dati = run.config.model_dump(mode="python")
    dati["io"]["out_root"] = str(cartella / "out")
    for gruppo, valori in sovrascrivi.items():
        dati[gruppo].update(valori)
    shutil.copytree(run.config.io.out_root, cartella / "out")
    return ProjectRun(valida(dati))


@pytest.fixture(scope="session")
def catena_calcolata(ridotta_calcolata, tmp_path_factory):
    """S4-S7 sulla versione ridotta, sopra S0-S3, una volta per la sessione.

    In sola lettura, come ``ridotta_calcolata``; ``None`` dove mancano R e
    Bioconductor.
    """
    from amplicon16s.runner.graph import Passo

    return _sopra(ridotta_calcolata, tmp_path_factory, "catena", Passo.S7)


@pytest.fixture(scope="session")
def tassonomia_calcolata(catena_calcolata, tmp_path_factory):
    """S8 sulla versione ridotta, sopra S0-S7, con il riferimento sintetico.

    In sola lettura, come ``catena_calcolata``; ``None`` dove mancano R e
    Bioconductor.
    """
    from amplicon16s.runner.graph import Passo

    return _sopra(catena_calcolata, tmp_path_factory, "tassonomia", Passo.S8)


@pytest.fixture(scope="session")
def oggetto_calcolato(tassonomia_calcolata, tmp_path_factory):
    """S10 sulla versione ridotta, sopra S0-S8, una volta per la sessione.

    In sola lettura, come ``tassonomia_calcolata``; ``None`` dove mancano R e
    Bioconductor.
    """
    from amplicon16s.runner.graph import Passo

    return _sopra(tassonomia_calcolata, tmp_path_factory, "oggetto", Passo.S10)


@pytest.fixture(scope="session")
def controlli_calcolati(oggetto_calcolato, tmp_path_factory):
    """S11 sulla versione ridotta, sopra S0-S10, una volta per la sessione.

    In sola lettura, come ``oggetto_calcolato``; ``None`` dove mancano R e
    Bioconductor.
    """
    from amplicon16s.runner.graph import Passo

    return _sopra(oggetto_calcolato, tmp_path_factory, "controlli", Passo.S11)


@pytest.fixture(scope="session")
def decontam_calcolata(controlli_calcolati, tmp_path_factory):
    """S12 sulla versione ridotta, sopra S0-S11, una volta per la sessione.

    In sola lettura, come ``controlli_calcolati``; ``None`` dove mancano R e
    Bioconductor.
    """
    from amplicon16s.runner.graph import Passo

    return _sopra(controlli_calcolati, tmp_path_factory, "decontam", Passo.S12)


@pytest.fixture(scope="session")
def finale_calcolata(decontam_calcolata, tmp_path_factory):
    """S13 e S14 sulla versione ridotta, sopra S0-S12: la catena intera.

    In sola lettura, come ``decontam_calcolata``; ``None`` dove mancano R e
    Bioconductor.
    """
    from amplicon16s.runner.graph import Passo

    return _sopra(decontam_calcolata, tmp_path_factory, "finale", Passo.S14)


@pytest.fixture(scope="session")
def catena_reale(tmp_path_factory):
    """S0-S14 sul dataset completo, una volta sola per la sessione.

    Serve ai test sui dati reali delle fasi da S4 in poi: ciascuna catena
    completa costa tre quarti d'ora, e condividerla fra i moduli e fra i
    processi evita di ripeterla. Con fino_a si eseguono solo gli antenati di S14: S1 no.
    Salta senza ``AMPLICON16S_CONFIG_DATI_REALI``; in sola lettura.
    """
    import logging
    import os

    from amplicon16s.config.schema import carica, valida
    from amplicon16s.runner.executor import Esecutore
    from amplicon16s.runner.graph import Passo
    from amplicon16s.runner.project import ProjectRun

    percorso = os.environ.get("AMPLICON16S_CONFIG_DATI_REALI")
    if not percorso or not Path(percorso).expanduser().is_file():
        pytest.skip("AMPLICON16S_CONFIG_DATI_REALI non impostata o file assente")
    logging.getLogger("amplicon16s").setLevel(logging.WARNING)

    def calcola(cartella):
        """La catena fino a S14 sul dataset completo, con gli artefatti in
        ``cartella/out``: la coppia (esecuzione, esito) che :func:`_condivisa`
        mette in comune."""
        dati = carica(Path(percorso).expanduser()).model_dump(mode="python")
        dati["io"]["out_root"] = str(cartella / "out")
        run = ProjectRun(valida(dati))
        return run, Esecutore(run, fino_a=Passo.S14).esegui()

    return _condivisa(tmp_path_factory, "reale", calcola)
