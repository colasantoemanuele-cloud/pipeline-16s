"""Il sottocomando ``init``: dai dati a una configurazione di partenza.

Uso, dalla radice di un clone e nel container::

    python3 scripts/esegui.py init --reads CARTELLA --metadata TABELLA_O_CARTELLA \\
        --output configurazione.yaml [--batch-table FILE] [--reference FILE] \\
        [--out-root CARTELLA]

``--metadata`` e' una tabella dei campioni separata da tabulazioni oppure una
cartella ISA-Tab (una tabella di assay ``a_*.txt`` e una di studio
``s_*.txt``). Il comando legge e non modifica nulla: non scarica, non tocca i
file di ingresso, e non scrive su un file di configurazione che esiste gia'.

Che cosa fa, in ordine: elenca i file di letture; se risultano appaiati (dal
nome o dal contenuto, con le regole di G07 e S1) si ferma senza scrivere e
dice di usare le sole forward; mette in corrispondenza file e tabella; cerca
nelle letture un primer di un elenco di primer 16S comuni; calcola il profilo
di qualita' con lo script di S1 e ne ricava il troncamento che il report
suggerirebbe; propone colonne ed etichette; scrive il file e un riepilogo.

**Codici di uscita**: ``0`` file scritto; ``3`` rifiuto (file gia' esistente,
letture appaiate, ingressi non trovati o non leggibili), nulla e' stato scritto.
"""

from __future__ import annotations

import argparse
import csv
import os
import statistics
import tempfile
from collections import Counter
from collections.abc import Callable
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from typing import Any, Final

from amplicon16s.config import defaults
from amplicon16s.config.schema import Config, processori_disponibili
from amplicon16s.errors.exceptions import ErrorePipeline
from amplicon16s.init import rilevamento as r
from amplicon16s.init import vocabolario as v
from amplicon16s.init.scrittura import SENZA_VALORE, Specie, Voce, testo_yaml
from amplicon16s.io_layer.artifacts import AlberoOutput, Fase, scrivi_atomico
from amplicon16s.io_layer.reads import espandi_iupac
from amplicon16s.metadata.models import ClasseCampione
from amplicon16s.rbridge.runner import cartella_r, esegui_script
from amplicon16s.report.builder import QUALITA_MINIMA_SUGGERITA, troncamento_suggerito
from amplicon16s.steps.s01_profile import (
    NOME_LUNGHEZZE,
    NOME_QUALITA,
    NOME_SCRIPT,
    NOME_VALORI_QUALITA,
    QUALITA_RAGGRUPPATE,
)

__all__ = ["USCITA_RIFIUTO", "Rifiuto", "aggiungi_argomenti", "esegui", "prepara"]

USCITA_RIFIUTO: Final = 3

#: La quota della memoria della macchina che la tabella delle sequenze puo'
#: occupare nel valore di partenza di ``qc.max_asv_count``, e i byte di una
#: cella (un intero di R).
_QUOTA_DI_MEMORIA: Final = 0.1
_BYTE_PER_CELLA: Final = 4

RILEVATO, INDICATO = Specie.RILEVATO, Specie.INDICATO
PROPOSTO, SCELTA = Specie.PROPOSTO, Specie.DA_SCEGLIERE


class Rifiuto(Exception):
    """``init`` non scrive alcuna configurazione, e dice perche'."""


def _intero_positivo(testo: str) -> int:
    """Un argomento intero di almeno uno."""
    try:
        valore = int(testo)
    except ValueError:
        raise argparse.ArgumentTypeError(f"{testo!r} non e' un numero intero") from None
    if valore < 1:
        raise argparse.ArgumentTypeError(f"{valore} non e' ammesso: serve almeno 1")
    return valore


def _sequenza_iupac(testo: str) -> str:
    """Una sequenza nucleotidica in codici IUPAC, in maiuscolo."""
    sequenza = testo.strip().upper()
    try:
        espandi_iupac(sequenza)
    except ValueError as guasto:
        raise argparse.ArgumentTypeError(str(guasto)) from None
    if not sequenza:
        raise argparse.ArgumentTypeError("sequenza vuota")
    return sequenza


def aggiungi_argomenti(parser: argparse.ArgumentParser) -> None:
    """Gli argomenti del sottocomando."""
    parser.add_argument("--reads", required=True, type=Path, metavar="CARTELLA",
                        help="la cartella con i file delle letture, uno per campione")
    parser.add_argument("--metadata", required=True, type=Path, metavar="FILE_O_CARTELLA",
                        help="la tabella dei campioni (TSV) oppure una cartella ISA-Tab")
    parser.add_argument("--assay", metavar="NOME",
                        help="in una cartella ISA-Tab con piu' tabelle di assay, il nome "
                             "del file da usare")
    parser.add_argument("--study-table", type=Path, metavar="FILE",
                        help="la tabella di studio, se --metadata e' una tabella di assay")
    parser.add_argument("--batch-table", type=Path, metavar="FILE",
                        help="il file del lotto (piastra e corsa per campione), se c'e'")
    parser.add_argument("--primer", type=_sequenza_iupac, metavar="SEQUENZA",
                        help="il primer forward usato per l'amplificazione, in codici IUPAC, "
                             "se non e' fra quelli che init riconosce da se'")
    parser.add_argument("--reference", type=Path, metavar="FILE",
                        help="il riferimento tassonomico, nel formato di dada2")
    parser.add_argument("--out-root", type=Path, metavar="CARTELLA",
                        help="la cartella di uscita prevista per l'esecuzione")
    parser.add_argument("--output", required=True, type=Path, metavar="FILE",
                        help="il file di configurazione da scrivere: non deve esistere")
    parser.add_argument("--head-reads", type=_intero_positivo, metavar="N",
                        help="ispeziona le sole prime N letture di ogni file per le coppie "
                             "e il primer (senza: tutte)")
    parser.add_argument("--skip-profile", action="store_true",
                        help="non calcola il profilo di qualita': nessuna proposta di "
                             "troncamento ne' di funzione d'errore")
    parser.add_argument("--threads", type=_intero_positivo, metavar="N",
                        help="processi per l'ispezione (senza: quelli utilizzabili)")


# --------------------------------------------------------------------------- #
# Ingressi                                                                     #
# --------------------------------------------------------------------------- #


def _tabelle(args: argparse.Namespace) -> tuple[Path, Path | None]:
    """La tabella di assay e quella di studio, dagli argomenti."""
    metadati: Path = args.metadata
    if not metadati.exists():
        raise Rifiuto(f"--metadata: {metadati} non esiste")
    if metadati.is_file():
        if args.study_table is not None and not args.study_table.is_file():
            raise Rifiuto(f"--study-table: {args.study_table} non esiste")
        return metadati, args.study_table
    assay = sorted(metadati.glob("a_*.txt"))
    if args.assay is not None:
        assay = [p for p in assay if p.name == args.assay]
        if not assay:
            raise Rifiuto(f"--assay: {args.assay} non e' fra le tabelle di assay di {metadati}")
    if not assay:
        raise Rifiuto(
            f"{metadati} non contiene una tabella di assay ISA-Tab (a_*.txt): indica con "
            "--metadata la tabella dei campioni")
    if len(assay) > 1:
        raise Rifiuto(
            f"{metadati} contiene {len(assay)} tabelle di assay: indica con --assay quella "
            "del sequenziamento 16S, fra " + ", ".join(repr(p.name) for p in assay))
    studio = args.study_table
    if studio is None:
        trovate = sorted(metadati.glob("s_*.txt"))
        if len(trovate) > 1:
            raise Rifiuto(
                f"{metadati} contiene {len(trovate)} tabelle di studio: indica quella da "
                "usare con --study-table")
        studio = trovate[0] if trovate else None
    return assay[0], studio


def _memoria_della_macchina() -> int | None:
    """I byte di memoria della macchina, o ``None`` se non si possono leggere."""
    try:
        return os.sysconf("SC_PHYS_PAGES") * os.sysconf("SC_PAGE_SIZE")
    except (ValueError, OSError, AttributeError):
        return None


def _due_cifre(valore: float) -> int:
    """Il valore arrotondato per difetto a due cifre significative."""
    intero = int(valore)
    if intero < 100:
        return max(intero, 1)
    scala = 10 ** (len(str(intero)) - 2)
    return intero // scala * scala


# --------------------------------------------------------------------------- #
# Profilo                                                                      #
# --------------------------------------------------------------------------- #


def _profilo(
    file_per_chiave: dict[str, Path], negativi: set[str], processi: int
) -> dict[str, Any]:
    """Il profilo delle letture con lo script di S1, e cio' che se ne ricava.

    Lo script e' quello della fase, eseguito in una cartella temporanea: le
    tabelle sono le stesse che S1 scrivera', e il troncamento suggerito e'
    quello che il report calcolerebbe da esse. La frazione ammessa di letture
    corte e' zero: nessuna tolleranza e' ancora stata scelta.
    """
    with tempfile.TemporaryDirectory(prefix="amplicon16s-init-") as temporanea:
        albero = AlberoOutput(temporanea)
        esegui_script(
            cartella_r() / NOME_SCRIPT,
            {"campioni": {c: str(p) for c, p in file_per_chiave.items()}, "processi": processi},
            albero, Fase.QC_PROFILES, passo="S1",
        )
        cartella = albero.cartella(Fase.QC_PROFILES)

        def righe(nome: str) -> list[dict[str, str]]:
            with open(cartella / nome, encoding="utf-8", newline="") as file:
                return list(csv.DictReader(file, delimiter="\t"))

        lunghezze, qualita = righe(NOME_LUNGHEZZE), righe(NOME_QUALITA)
        valori = [int(x["qualita"]) for x in righe(NOME_VALORI_QUALITA)]
    classi = {
        c: (ClasseCampione.CONTROLLO_NEGATIVO if c in negativi else ClasseCampione.BIOLOGICO).value
        for c in file_per_chiave
    }
    suggerito = troncamento_suggerito(lunghezze, qualita, classi, 0.0)
    distribuzione: Counter[int] = Counter()
    for riga in lunghezze:
        if riga["campione"] not in negativi:
            distribuzione[int(riga["lunghezza"])] += int(float(riga["letture"]))
    totale = sum(distribuzione.values())
    moda = max(distribuzione, key=lambda x: (distribuzione[x], x)) if distribuzione else None
    return {
        "suggerito": suggerito,
        "valori_di_qualita": sorted(valori),
        "minima": min(distribuzione) if distribuzione else None,
        "massima": max(distribuzione) if distribuzione else None,
        "moda": moda,
        "perdita_alla_moda": (
            sum(n for lunghezza, n in distribuzione.items() if lunghezza < moda) / totale
            if moda is not None and totale else None
        ),
    }


# --------------------------------------------------------------------------- #
# Dalle osservazioni alle voci                                                 #
# --------------------------------------------------------------------------- #


def prepara(args: argparse.Namespace, avanzamento: Callable[[str], None] = lambda _: None) -> tuple[
    dict[str, Voce], list[str], list[str], list[str]
]:
    """Esamina i dati e restituisce le voci della configurazione, le righe
    d'intestazione, cio' che ``init`` non stabilisce e le note del riepilogo.

    Solleva :class:`Rifiuto` se non c'e' una configurazione da proporre.
    """
    voci: dict[str, Voce] = {}
    note: list[str] = []
    non_indovinati: list[str] = []

    def metti(chiave: str, specie: Specie, valore: Any, origine: str, *altre: str) -> None:
        voci[chiave] = Voce(chiave, specie, valore, origine, tuple(a for a in altre if a))

    # ---- Letture ---------------------------------------------------------
    cartella: Path = args.reads
    if not cartella.is_dir():
        raise Rifiuto(f"--reads: {cartella} non e' una cartella")
    per_estensione = r.file_di_letture(cartella)
    if not per_estensione:
        raise Rifiuto(
            f"{cartella} non contiene file di letture ({', '.join(v.ESTENSIONI_FASTQ)})")
    file = [p for lista in per_estensione.values() for p in lista]
    nomi = [p.name for p in file]
    inverse = r.secondi_di_coppia(nomi)
    if inverse:
        raise Rifiuto(_letture_appaiate(
            f"{len(inverse)} file hanno il nome della seconda lettura di una coppia "
            f"({', '.join(inverse[:3])}{'...' if len(inverse) > 3 else ''})"))

    avanzamento(f"letture: {len(file)} file in {cartella}")
    # Il primer indicato da chi lancia il comando si cerca nelle letture come
    # quelli dell'elenco, e viene per primo.
    cercati = v.PRIMER_16S
    if args.primer is not None:
        cercati = (v.Primer("indicato", args.primer, "--primer"),
                  *(p for p in v.PRIMER_16S if p.sequenza != args.primer))
    primer = tuple(espandi_iupac(p.sequenza) for p in cercati)
    processi = args.threads or processori_disponibili()
    with ProcessPoolExecutor(max_workers=processi) as gruppo:
        ispezioni = dict(zip(nomi, gruppo.map(
            r.ispeziona, [(str(p), args.head_reads, primer) for p in file])))
    intercalati = sorted(
        f"{n} ({i.statistiche.coppie_nello_stesso_file})" for n, i in ispezioni.items()
        if i.statistiche.valido and i.statistiche.coppie_nello_stesso_file)
    if intercalati:
        raise Rifiuto(_letture_appaiate(
            f"{len(intercalati)} file non contengono le sole prime letture delle coppie: "
            + "; ".join(intercalati[:3]) + ("..." if len(intercalati) > 3 else "")))
    guasti = sorted(n for n, i in ispezioni.items()
                    if not i.statistiche.valido and not i.statistiche.vuoto)
    vuoti = sorted(n for n, i in ispezioni.items() if i.statistiche.vuoto)
    if guasti:
        note.append(
            f"{len(guasti)} file non si leggono come FASTQ ({', '.join(guasti[:3])}"
            f"{'...' if len(guasti) > 3 else ''}): validate li segnalera' (E-S0-13).")
    if vuoti:
        note.append(
            f"{len(vuoti)} file non contengono letture ({', '.join(vuoti[:3])}"
            f"{'...' if len(vuoti) > 3 else ''}): ammesso solo per i controlli negativi.")
    quante = "tutte le letture" if args.head_reads is None else f"le prime {args.head_reads} letture"

    metti("io.fastq_dir", INDICATO, str(cartella), "--reads")
    if len(per_estensione) == 1:
        (estensione,) = per_estensione
        metti("io.fastq_glob", RILEVATO, f"*{estensione}",
              f"i {len(file)} file di letture della cartella hanno tutti questa estensione")
    else:
        metti("io.fastq_glob", SCELTA, SENZA_VALORE,
              "la cartella contiene file di letture con estensioni diverse ("
              + ", ".join(f"{len(l)} {e}" for e, l in per_estensione.items())
              + "): indica il modello che seleziona quelli da usare")

    # ---- Tabelle ----------------------------------------------------------
    assay, studio = _tabelle(args)
    try:
        colonne_assay, righe_assay = r.leggi_tabella(assay)
        colonne_studio, righe_studio = r.leggi_tabella(studio) if studio else ([], None)
    except ValueError as guasto:
        raise Rifiuto(str(guasto)) from guasto
    metti("io.assay_table", INDICATO, str(assay), "--metadata")
    metti("io.study_table", INDICATO, str(studio) if studio else None,
          "--study-table, o la tabella s_*.txt della cartella ISA-Tab" if studio else
          "nessuna tabella di studio: tipo e variabili si leggono dalla tabella dei campioni")
    if args.out_root is not None:
        metti("io.out_root", INDICATO, str(args.out_root), "--out-root")

    # ---- Chiave del campione ---------------------------------------------
    candidate = r.chiave_del_campione(nomi, righe_assay, colonne_assay)
    di_file = [c for c in candidate if c.di_file]
    migliori = di_file or candidate
    chiave: r.Chiave | None = migliori[0] if len(migliori) == 1 else None

    # ---- File del lotto ----------------------------------------------------
    lotto: Path | None = args.batch_table
    colonne_lotto: list[str] = []
    righe_lotto: list[dict[str, str]] = []
    if lotto is not None:
        if not lotto.is_file():
            raise Rifiuto(f"--batch-table: {lotto} non esiste")
        try:
            colonne_lotto, righe_lotto = r.leggi_tabella(lotto)
        except ValueError as guasto:
            raise Rifiuto(str(guasto)) from guasto
    metti("io.batch_table", INDICATO, str(lotto) if lotto else None,
          "--batch-table" if lotto else "nessun file del lotto indicato: corsa singola, nessun lotto")

    colonne_chiave_lotto: list[str] = []
    if chiave is not None and lotto is not None:
        # Fra le espressioni verificate si preferisce quella con cui anche il
        # file del lotto si aggancia: la chiave deve essere la stessa ovunque.
        tutte = (chiave.espressione, *chiave.alternative)
        for espressione in tutte:
            chiavi = r.chiavi_dei_file(espressione, nomi)
            trovate = r.chiave_del_lotto(chiavi.values(), espressione, righe_lotto, colonne_lotto)
            if trovate:
                chiave = r.Chiave(
                    chiave.colonna, espressione, chiave.di_file,
                    {c: n for n, c in chiavi.items()},
                    tuple(e for e in tutte if e != espressione))
                colonne_chiave_lotto = trovate
                break

    if chiave is not None:
        # Una chiave di sole cifre o cortissima puo' corrispondere ai file per
        # coincidenza: l'abbinamento e' allora una proposta, da guardare.
        specie_chiave = PROPOSTO if chiave.debole else RILEVATO
        esempi = "abbinamento, per esempio: " + "; ".join(
            f"{nome} -> {c}" for c, nome in list(chiave.file_per_chiave.items())[:3])
        avvertenza = (
            "attenzione: la chiave e' di sole cifre o molto corta, e la corrispondenza "
            "con i file puo' essere una coincidenza: controlla l'abbinamento"
        ) if chiave.debole else ""
        metti("meta.accession_column", specie_chiave, chiave.colonna,
              "la sola colonna "
              + ("con nomi di file di letture " if chiave.di_file else "")
              + f"i cui valori corrispondono uno a uno ai file ({len(nomi)})",
              esempi, avvertenza)
        metti("io.accession_regex", specie_chiave, chiave.espressione,
              f"estrae da ogni file una chiave unica, la stessa che si ricava dalla "
              f"colonna {chiave.colonna!r}: verificato su tutti i file e le righe ({len(nomi)})",
              ("altre espressioni verificate: " + ", ".join(chiave.alternative[:4]))
              if chiave.alternative else "", avvertenza)
        if chiave.debole:
            note.append(
                "La chiave del campione e' proposta e non rilevata: e' di sole cifre o molto "
                "corta, e la corrispondenza con i file puo' essere una coincidenza.")
    else:
        if migliori:
            elenco = "; ".join(f"{c.colonna!r} con {c.espressione!r}" for c in migliori)
            origine = (f"{len(migliori)} colonne corrispondono uno a uno ai file, e init non "
                       f"sceglie fra loro: {elenco}")
        else:
            origine = (
                "nessuna colonna della tabella corrisponde uno a uno ai file: i valori "
                f"devono contenere una parte del nome del file unica per campione ("
                f"{len(nomi)} file, {len(righe_assay)} righe)")
        metti("meta.accession_column", SCELTA, SENZA_VALORE, origine)
        metti("io.accession_regex", SCELTA, SENZA_VALORE,
              "si sceglie insieme a meta.accession_column: la parte del nome del file "
              "che e' unica per campione")
        note.append("Chiave del campione non rilevata: " + origine + ".")

    # ---- Nome del campione --------------------------------------------------
    agganci = r.colonna_del_campione(righe_assay, colonne_assay, righe_studio, colonne_studio)
    colonna_campione: str | None = None
    colonna_studio: str | None = None
    if len(agganci) == 1:
        colonna_campione, colonna_studio = agganci[0]
        metti("meta.sample_id_column", RILEVATO, colonna_campione,
              "la sola colonna che aggancia ogni riga della tabella di assay a una e una "
              "sola riga della tabella di studio")
        if colonna_studio != colonna_campione:
            metti("meta.study_sample_id_column", RILEVATO, colonna_studio,
                  "la colonna della tabella di studio che corrisponde a meta.sample_id_column")
    elif righe_studio is not None:
        metti("meta.sample_id_column", SCELTA, SENZA_VALORE,
              (f"{len(agganci)} coppie di colonne agganciano le due tabelle: "
               + "; ".join(f"{a!r} con {s!r}" for a, s in agganci[:6])) if agganci else
              "nessuna colonna aggancia ogni riga della tabella di assay a una sola riga "
              "della tabella di studio")
    else:
        uniche = [c for c in colonne_assay
                  if c != (chiave.colonna if chiave else None)
                  and len({x.get(c, "") for x in righe_assay} - {""}) == len(righe_assay)]
        if uniche:
            colonna_campione = uniche[0]
            metti("meta.sample_id_column", PROPOSTO, colonna_campione,
                  "la prima colonna con un valore diverso per ogni riga",
                  ("altre colonne con valori tutti distinti: " + ", ".join(map(repr, uniche[1:5])))
                  if uniche[1:] else "")
        else:
            metti("meta.sample_id_column", SCELTA, SENZA_VALORE,
                  "nessuna colonna, oltre a quella della chiave, ha un valore diverso per "
                  "ogni riga")

    # ---- Tipo di campione ----------------------------------------------------
    righe_classe = righe_assay
    colonne_classe = colonne_assay
    if righe_studio is not None:
        colonne_classe = colonne_studio
        if colonna_campione is not None:
            nomi_assay = {x.get(colonna_campione, "") for x in righe_assay}
            righe_classe = [x for x in righe_studio if x.get(colonna_studio or "", "") in nomi_assay]
        else:
            righe_classe = righe_studio
    escluse = {chiave.colonna if chiave else None, colonna_campione, colonna_studio}
    tipi = r.colonne_del_tipo(righe_classe, [c for c in colonne_classe if c not in escluse])
    negativi_proposti: set[str] = set()
    con_positivi = False
    if tipi:
        tipo = tipi[0]
        altre = ("altre colonne candidate: " + ", ".join(repr(t.colonna) for t in tipi[1:5])
                 ) if tipi[1:] else ""
        come = ("ha valori con le parole dei controlli" if tipo.da == "parole" else
                "il nome dice il tipo del campione; nessun valore ha le parole dei controlli")
        metti("ctrl.column", PROPOSTO, tipo.colonna, f"colonna proposta perche' {come}", altre)
        metti("ctrl.blank_values", PROPOSTO, list(tipo.negativi),
              "i valori con una parola dei controlli negativi ("
              + ", ".join(v.PAROLE_NEGATIVO) + ")" if tipo.negativi else
              "nessun valore ha una parola dei controlli negativi")
        metti("ctrl.positive_values", PROPOSTO, list(tipo.positivi),
              "i valori con una parola dei controlli positivi ("
              + ", ".join(v.PAROLE_POSITIVO) + ")" if tipo.positivi else
              "nessun valore ha una parola dei controlli positivi")
        metti("ctrl.biological_values", PROPOSTO, list(tipo.biologici),
              "i valori senza parole dei controlli",
              ("valori con la sola parola «control», non attribuiti: da assegnare a mano a "
               "una delle tre etichette: " + ", ".join(map(repr, tipo.non_attribuiti)))
              if tipo.non_attribuiti else "")
        con_positivi = bool(tipo.positivi)
        if chiave is not None and colonna_campione is not None:
            tipo_per_nome = {x.get(colonna_studio or colonna_campione, ""): x.get(tipo.colonna, "")
                             for x in righe_classe}
            for riga in righe_assay:
                if tipo_per_nome.get(riga.get(colonna_campione, "")) in tipo.negativi:
                    della_riga = _chiave_della_riga(chiave, riga)
                    if della_riga is not None:
                        negativi_proposti.add(della_riga)
    else:
        metti("ctrl.column", PROPOSTO, None,
              "nessuna colonna ha valori con le parole dei controlli ne' un nome che dica "
              "il tipo del campione: nullo vuol dire nessun controllo, ogni campione e' "
              "biologico")
        for nome in ("blank_values", "positive_values", "biological_values"):
            metti(f"ctrl.{nome}", PROPOSTO, [], "senza la colonna del tipo le etichette sono vuote")
    metti("ctrl.blank_override_column", SCELTA, None,
          "partenza: nessuna riclassificazione. Campioni dichiarati di un tipo ma da "
          "trattare come controlli negativi (per esempio materiale mai usato) si "
          "riconoscono solo conoscendo lo studio: init non li puo' rilevare")
    metti("ctrl.blank_override_values", SCELTA, [], "partenza: nessuna riclassificazione")

    # ---- Modulo, colonne da portare nell'oggetto ---------------------------------
    metti("meta.derive_module", SCELTA, False,
          "partenza: nessun raggruppamento di posizioni da derivare")
    metti("meta.module_regex", SCELTA, None,
          "partenza: nessuna; serve solo se le posizioni hanno un prefisso che dice il gruppo")
    metti("meta.module_column", SCELTA, None,
          "partenza: nessuna colonna della posizione di prelievo",
          "colonne della tabella: " + ", ".join(map(repr, colonne_classe[:20])))
    metti("meta.non_surface_positions", SCELTA, [],
          "partenza: nessuna; ha senso solo con meta.module_column")
    metti("out.study_columns", SCELTA, [],
          "partenza: nessuna colonna oltre a quelle dell'inventario; elenca qui le "
          "variabili dello studio che le analisi a valle useranno",
          "colonne disponibili: " + ", ".join(map(repr, colonne_classe[:20])))

    # ---- File del lotto: colonne ---------------------------------------------------
    if lotto is None:
        for nome in ("meta.batch_key_column", "meta.batch_module_column", "err.batch_column",
                     "decontam.batch_column"):
            metti(nome, RILEVATO, None,
                  "nessun file del lotto indicato: e' il solo valore ammesso (E-G15-11)")
        metti("out.batch_columns", RILEVATO, [],
              "nessun file del lotto indicato: e' il solo valore ammesso (E-G15-11)")
        metti("katharoseq.cell_count_column", SCELTA, None,
              "partenza: nessuna; senza file del lotto la colonna delle cellule dei "
              "controlli positivi andrebbe cercata nella tabella dei campioni")
    else:
        if len(colonne_chiave_lotto) == 1:
            metti("meta.batch_key_column", RILEVATO, colonne_chiave_lotto[0],
                  "la sola colonna del file del lotto che da' per ogni campione la sua "
                  "chiave, in una e una sola riga")
        else:
            metti("meta.batch_key_column", SCELTA, SENZA_VALORE,
                  (f"{len(colonne_chiave_lotto)} colonne del file del lotto danno la chiave "
                   "di ogni campione: " + ", ".join(map(repr, colonne_chiave_lotto)))
                  if colonne_chiave_lotto else
                  "nessuna colonna del file del lotto da' per ogni campione la chiave "
                  "estratta dai file" + ("" if chiave else " (la chiave non e' stata rilevata)"))
        del_lotto = r.colonne_del_lotto(
            colonne_lotto, righe_lotto,
            colonne_chiave_lotto[0] if len(colonne_chiave_lotto) == 1 else None)
        for nome, tipo_colonna, che_cosa in (
            ("decontam.batch_column", "piastra", "la piastra di estrazione"),
            ("err.batch_column", "corsa", "la corsa di sequenziamento"),
            ("katharoseq.cell_count_column", "cellule", "le cellule dei controlli positivi"),
        ):
            trovate = del_lotto[tipo_colonna]
            if trovate:
                metti(nome, PROPOSTO, trovate[0],
                      f"il nome della colonna dice {che_cosa}",
                      ("altre colonne candidate: " + ", ".join(map(repr, trovate[1:])))
                      if trovate[1:] else "")
            else:
                metti(nome, SCELTA, None,
                      f"partenza: nessuna; nessuna colonna del file del lotto ha un nome "
                      f"che dica {che_cosa}",
                      "colonne del file: " + ", ".join(map(repr, colonne_lotto)))
        metti("meta.batch_module_column", SCELTA, None,
              "partenza: nessuna colonna del modulo nel file del lotto")
        metti("out.batch_columns", SCELTA, [],
              "partenza: nessuna colonna del file del lotto oltre a piastra e corsa",
              "colonne disponibili: " + ", ".join(map(repr, colonne_lotto)))
    con_cellule = con_positivi and voci["katharoseq.cell_count_column"].valore not in (None, SENZA_VALORE)

    # ---- Primer e motivo -----------------------------------------------------------
    osservato, parziali = r.primer_nelle_letture(ispezioni, cercati)
    da_considerare = [
        i for n, i in ispezioni.items()
        if chiave is None or _chiave_del_file(chiave, n) not in negativi_proposti
    ] or list(ispezioni.values())
    esclusi = len(ispezioni) - len(da_considerare)
    senza_negativi = (
        f"esclusi {esclusi} file dei campioni proposti come controlli negativi" if esclusi else "")
    dichiarato, sequenze, dichiarato_forward = r.primer_dai_metadati(
        [righe_assay] + ([righe_studio] if righe_studio is not None else []))
    primer_noto: str | None = None
    taglio = 0
    if osservato is not None:
        trovato = osservato.primer
        primer_noto, taglio = trovato.sequenza, len(trovato.sequenza)
        dove = (f"in testa alle letture di ogni file: dal {osservato.minima:.0%} al "
                f"{osservato.massima:.0%} delle letture, su {quante}")
        # E' un fatto solo se vale per ogni file e su tutte le letture: un
        # taglio fisso toglierebbe basi vere ai file che non portano il primer.
        if osservato.in_ogni_file and args.head_reads is None:
            specie_primer, limite = RILEVATO, ""
        elif not osservato.in_ogni_file:
            sotto = osservato.sotto_soglia
            specie_primer = PROPOSTO
            dove = (f"in testa alle letture nella mediana dei file ({osservato.mediana:.0%}), "
                    f"ma non in tutti: dal {osservato.minima:.0%} al {osservato.massima:.0%}")
            limite = (
                f"attenzione: {len(sotto)} file su {len(ispezioni)} non lo portano in testa "
                f"({', '.join(sotto[:5])}{'...' if len(sotto) > 5 else ''}): un taglio fisso "
                "toglierebbe loro basi vere. Vanno capiti prima di accettare")
            note.append(
                f"Il primer {trovato.nome} e' in testa alle letture della maggior parte dei "
                f"file ma non di {len(sotto)}: primer e taglio sono proposti, non rilevati.")
        else:
            specie_primer = PROPOSTO
            limite = ("verificato sulle sole prime letture (--head-reads): senza l'opzione "
                      "init lo verifica su tutte")
        metti("qc.primer_sequence", specie_primer, trovato.sequenza,
              f"primer {trovato.nome} ({trovato.fonte}), {dove}", limite)
        metti("filter.trimLeft", specie_primer, taglio,
              f"la lunghezza del primer {trovato.nome}, che le letture portano in testa", limite)
    else:
        variabile = "; ".join(f"{p.nome} nel {f:.1%}" for p, f in parziali)
        metti("filter.trimLeft", PROPOSTO, 0,
              f"nessun primer dell'elenco e' in testa alle letture (su {quante}); un primer "
              "fuori elenco va indicato a mano, con la sua lunghezza qui",
              "primer cercati: " + ", ".join(p.nome for p in cercati),
              (f"attenzione: una parte delle letture inizia con un primer ({variabile}): e' "
               "in posizione variabile, e un taglio fisso non lo toglie") if variabile else "")
        if variabile:
            note.append(
                f"Una parte delle letture inizia con un primer dell'elenco ({variabile}), "
                "sotto la soglia per dirlo in testa: e' probabilmente preceduto da basi di "
                "lunghezza variabile e va rimosso per sequenza prima della pipeline.")
        nominati = r.primer_nominati(
            [righe_assay] + ([righe_studio] if righe_studio is not None else []))
        if args.primer is not None:
            primer_noto = args.primer
            metti("qc.primer_sequence", INDICATO, args.primer,
                  "--primer; non e' in testa alle letture")
        elif dichiarato is not None:
            primer_noto = dichiarato
            metti("qc.primer_sequence", PROPOSTO, dichiarato,
                  ("sequenza che i metadati dichiarano come primer forward; non e' nelle letture"
                   if dichiarato_forward else
                   "la sola sequenza nucleotidica trovata nei metadati, senza indicazione che "
                   "sia il primer forward: da verificare; non e' nelle letture"),
                  ("sequenze trovate nei metadati: " + ", ".join(sequenze)) if len(sequenze) > 1 else "")
        else:
            metti("qc.primer_sequence", SCELTA, SENZA_VALORE,
                  "il primer forward non e' nelle letture e i metadati non ne dichiarano "
                  "la sequenza: va indicato quello usato per l'amplificazione",
                  ("i metadati riportano sequenze senza dire quale sia il forward: "
                   + ", ".join(sequenze)) if sequenze else "",
                  ("i metadati nominano un primer senza riportarne la sequenza; quelle che "
                   "init conosce con quel nome: "
                   + "; ".join(f"{p.nome} {p.sequenza} ({p.fonte})" for p in nominati)
                   + ". Se e' una di queste, scrivila qui o rilancia init con --primer")
                  if nominati else "")
            non_indovinati.append(
                "qc.primer_sequence: il primer forward non e' nelle letture ne' nei metadati.")
    motivo = None
    if primer_noto is not None:
        motivo = r.motivo_dal_consenso(da_considerare, taglio, len(primer_noto))
    if motivo is not None:
        metti("qc.conserved_motif", PROPOSTO, motivo,
              f"il consenso delle {len(primer_noto or '')} basi "
              + ("che seguono il primer" if taglio else "con cui iniziano le letture")
              + f", su non piu' di {r.LETTURE_PER_IL_MOTIVO} letture per file: una base o "
              f"una coppia di basi dove copre almeno il {r.SOGLIA_BASE_DEL_MOTIVO:.0%} delle letture",
              "validate ne misura la frequenza con qc.min_motif_frac", senza_negativi)
    else:
        metti("qc.conserved_motif", SCELTA, None,
              "partenza: nullo (si verifica la sola assenza del primer). "
              + ("Il consenso delle letture non da' un motivo abbastanza definito"
                 if primer_noto else
                 "Senza conoscere il primer init non propone il motivo della regione amplificata"))
        if primer_noto is None:
            non_indovinati.append(
                "qc.conserved_motif: senza il primer la regione amplificata non e' nota "
                "(nullo e' ammesso).")

    # ---- Profilo: troncamento e funzione d'errore --------------------------------
    profilo: dict[str, Any] | None = None
    if args.skip_profile:
        motivo_senza_profilo = "il profilo di qualita' non e' stato calcolato (--skip-profile)"
    else:
        avanzamento("profilo di qualita' con lo script della fase S1...")
        leggibili = {
            (_chiave_del_file(chiave, n) if chiave else n): p
            for n, p in zip(nomi, file) if ispezioni[n].statistiche.valido
        }
        try:
            profilo = _profilo(leggibili, negativi_proposti, processi)
            motivo_senza_profilo = ""
        except ErrorePipeline as guasto:
            motivo_senza_profilo = (
                f"il profilo di qualita' non e' stato calcolato ({guasto.codice}: "
                f"{guasto.dettaglio or guasto.voce.sintesi})")
            note.append(
                "Profilo di qualita' non calcolato: " + motivo_senza_profilo
                + ". init va eseguito nel container, come gli altri comandi.")
    suggerito = profilo["suggerito"]["suggerito"] if profilo else None
    if profilo and suggerito is not None:
        s = profilo["suggerito"]
        metti("filter.truncLen", PROPOSTO, int(suggerito),
              f"il minore fra la lettura piu' corta dei campioni non di controllo negativo "
              f"({s['senza_perdite']}) e l'ultima posizione con qualita' mediana almeno "
              f"{QUALITA_MINIMA_SUGGERITA} ({s['per_qualita']}): non scarta alcuna lettura",
              f"lunghezze osservate: da {profilo['minima']} a {profilo['massima']}, moda "
              f"{profilo['moda']}; troncando alla moda si scarterebbe il "
              f"{profilo['perdita_alla_moda']:.2%} delle letture",
              "comprende le basi di filter.trimLeft", senza_negativi)
        scarto = max(int(s["senza_perdite"]) - int(suggerito), 0)
        metti("filter.truncLen_shortfall_warn", SCELTA, min(scarto, int(suggerito) - 1),
              "partenza: lo scarto fra la lettura piu' corta e il troncamento proposto; con "
              "un valore minore il profilo segnalerebbe il troncamento proposto (E-S1-01)")
        metti("qc.max_frac_short_reads", SCELTA, 0.0,
              "partenza: zero, perche' il troncamento proposto non scarta letture; va "
              "alzata solo se si sceglie un troncamento piu' lungo")
    else:
        metti("filter.truncLen", SCELTA, SENZA_VALORE,
              (motivo_senza_profilo or "le letture non bastano a suggerire un troncamento")
              + ": si sceglie dal profilo delle lunghezze e della qualita'")
        metti("filter.truncLen_shortfall_warn", SCELTA, SENZA_VALORE,
              "si sceglie insieme a filter.truncLen")
        metti("qc.max_frac_short_reads", SCELTA, SENZA_VALORE,
              "si sceglie insieme a filter.truncLen")
    if profilo and len(profilo["valori_di_qualita"]) <= QUALITA_RAGGRUPPATE:
        metti("err.error_function", PROPOSTO, "loess_monotono",
              f"le letture hanno {len(profilo['valori_di_qualita'])} valori di qualita' "
              f"distinti ({', '.join(map(str, profilo['valori_di_qualita']))}): qualita' "
              "raggruppate; il predefinito e' loess (dada2::learnErrors)")

    # ---- Riferimento ---------------------------------------------------------------
    riferimento: Path | None = args.reference
    if riferimento is not None:
        if not riferimento.is_file():
            raise Rifiuto(f"--reference: {riferimento} non esiste")
        metti("tax.ref_fasta", INDICATO, str(riferimento), "--reference")
        metti("tax.ref_md5", RILEVATO, r.md5_di(riferimento), "calcolato sul file indicato")
        nome, versione = r.nome_e_versione_del_riferimento(riferimento.name)
        for parametro, valore, che_cosa in (("tax.ref_name", nome, "il nome"),
                                            ("tax.ref_version", versione, "la versione")):
            if valore is not None:
                metti(parametro, PROPOSTO, valore, f"letto dal nome del file ({riferimento.name})")
            else:
                metti(parametro, SCELTA, SENZA_VALORE,
                      f"{che_cosa} del riferimento non si legge dal nome del file "
                      f"({riferimento.name})")
    else:
        for parametro in ("tax.ref_name", "tax.ref_version"):
            metti(parametro, SCELTA, SENZA_VALORE, "nessun riferimento indicato (--reference)")

    # ---- Scelte di analisi: valori di partenza -----------------------------------
    letture = [i.statistiche.letture_esaminate for i in ispezioni.values() if i.statistiche.valido]
    mediana = int(statistics.median(letture)) if letture else 0
    metti("qc.head_reads", SCELTA, SENZA_VALORE,
          "letture ispezionate per file dalla validazione: abbastanza da rappresentare il "
          "file, poche abbastanza da validare in minuti",
          (f"il file mediano ha {mediana} letture: un valore pari o superiore legge quel "
           "file per intero") if args.head_reads is None else
          f"init ha ispezionato le prime {args.head_reads} letture di ogni file")
    metti("err.randomize", SCELTA, False,
          "partenza: il predefinito di dada2::learnErrors (randomize = FALSE); true se "
          "l'ordine dei file segue una struttura dello studio")
    metti("dada.pool", SCELTA, False,
          "partenza: il predefinito di dada2::dada (pool = FALSE); \"pseudo\" recupera le "
          "varianti rare presenti in piu' campioni")
    metti("chimera.min_fold_parent_over_abundance", SCELTA, 1.5,
          "partenza: il predefinito del metodo di consenso (dada2::isBimeraDenovoTable, "
          "minFoldParentOverAbundance = 1.5)")
    metti("tax.try_rc", SCELTA, False,
          "partenza: il predefinito di dada2::assignTaxonomy (tryRC = FALSE); true se le "
          "letture possono essere nell'orientamento opposto al riferimento")
    metti("filt.remove_na_phylum", SCELTA, SENZA_VALORE,
          "dipende da che cosa lo studio considera fuori bersaglio")
    metti("filt.exclude_taxa", SCELTA, SENZA_VALORE,
          "i nomi dipendono dallo studio e da come li scrive il riferimento; [] per non "
          "togliere nulla")
    metti("phylo.enabled", SCELTA, False,
          "partenza: disattivata; la filogenesi e' una fase facoltativa, che serve alle "
          "distanze filogenetiche fra campioni")
    metti("phylo.max_seqs", SCELTA, SENZA_VALORE,
          "con phylo.enabled false non ha effetto; altrimenti si sceglie dopo aver visto "
          "nel report quante varianti restano")
    for parametro, origine in (
        ("prev.apply", "se filtrare le varianti presenti in pochi campioni"),
        ("prev.min_fraction", "dipende da quanto sono eterogenei i campioni"),
        ("prev.min_count", "letture minime perche' una variante conti come presente"),
    ):
        metti(parametro, SCELTA, SENZA_VALORE, origine)

    # Controlli positivi e calibrazione.
    metti("ctrl.min_positives", SCELTA, 3,
          "partenza: un punto piu' dei due parametri della curva di calibrazione, il "
          "minimo perche' l'adattamento abbia un residuo"
          + ("" if con_positivi else "; senza controlli positivi non ha effetto"))
    metti("ctrl.min_positive_pass_frac", SCELTA, SENZA_VALORE,
          "quanti controlli positivi devono essere conformi al loro livello"
          + ("" if con_positivi else "; senza controlli positivi non ha effetto"))
    metti("ctrl.positive_gate", SCELTA, False,
          "partenza: un avviso e non un arresto, finche' la conformita' dei controlli "
          "non e' nota dal report della prima esecuzione")
    if con_positivi:
        metti("katharoseq.target_taxon", SCELTA, SENZA_VALORE,
              "il taxon atteso nei controlli positivi: lo sa solo chi li ha allestiti")
        non_indovinati.append(
            "katharoseq.target_taxon: il taxon dei controlli positivi.")
    else:
        metti("katharoseq.target_taxon", SCELTA, None,
              "partenza: nullo, coerente con l'assenza di controlli positivi nella proposta")
    metti("katharoseq.collapse_rank", SCELTA, "Genus",
          "partenza: il rango piu' fine che la pipeline assegna")
    metti("katharoseq.target_sensitivity", SCELTA, SENZA_VALORE,
          "la fedelta' a cui si legge la profondita' minima: una scelta su quanta "
          "contaminazione si accetta nei campioni a bassa biomassa"
          + ("" if con_positivi else "; senza controlli positivi non ha effetto"))
    metti("katharoseq.min_r2", SCELTA, SENZA_VALORE,
          "la bonta' minima di una curva valida"
          + ("" if con_positivi else "; senza controlli positivi non ha effetto"))
    metti("qc.min_reads_mode", SCELTA, "katharoseq_if_available" if con_cellule else "none",
          "partenza: le soglie dei controlli positivi, che la proposta riconosce con la "
          "colonna delle cellule" if con_cellule else
          "partenza: nessuna soglia di profondita', perche' la proposta non ha controlli "
          "positivi con la colonna delle cellule")

    # Decontaminazione.
    metti("decontam.threshold", SCELTA, SENZA_VALORE,
          "la soglia di decontaminazione dipende dalla biomassa: init non la propone",
          "il predefinito di decontam::isContaminant e' 0.1; un valore piu' alto toglie di piu'")
    metti("decontam.min_blanks", SCELTA, SENZA_VALORE,
          "quanti controlli negativi con letture servono per stimare i contaminanti: "
          "dipende dal disegno dei controlli")
    metti("decontam.mode", SCELTA, "aggregate",
          "partenza: tutti i biologici contro tutti i negativi; l'altra modalita' si "
          "calcola comunque come diagnostica")
    non_indovinati += [
        "decontam.threshold: la soglia di decontaminazione.",
        "Ogni scelta legata alla biomassa o al disegno dello studio: le soglie di qualita' "
        "(qc.*), la prevalenza (prev.*), i taxa da escludere (filt.*), la fedelta' della "
        "calibrazione, la riclassificazione di campioni in controlli negativi.",
    ]

    # Soglie di qualita'.
    memoria = _memoria_della_macchina()
    if memoria is not None:
        massimo = memoria * _QUOTA_DI_MEMORIA / (_BYTE_PER_CELLA * max(len(file), 1))
        if args.head_reads is None and letture:
            massimo = min(massimo, sum(letture))
        metti("qc.max_asv_count", SCELTA, _due_cifre(massimo),
              f"partenza: la tabella delle sequenze piu' grande che occupa il "
              f"{_QUOTA_DI_MEMORIA:.0%} della memoria di questa macchina "
              f"({memoria / 2**30:.0f} GB), con {len(file)} campioni e {_BYTE_PER_CELLA} byte "
              "per cella" + (", non oltre il numero di letture" if args.head_reads is None
                             else "; le letture non sono state contate (--head-reads)"))
    else:
        metti("qc.max_asv_count", SCELTA, SENZA_VALORE, "secondo la memoria disponibile")
    metti("qc.max_zeroed_samples", SCELTA, 0,
          "partenza: nessun campione biologico perso nel filtro senza che lo si veda")
    for parametro, misura in (
        ("qc.max_primer_hit_frac", "frazione di letture con il primer in testa, in "
                                   "01_input_validation/letture_ispezionate.tsv"),
        ("qc.min_motif_frac", "frazione di letture con il motivo, nello stesso file"),
        ("qc.max_frac_lost_filter", "perdita media nel filtro, nel tracciamento delle letture del report"),
        ("qc.warn_frac_chimeric", "frazione di letture chimeriche, nel report"),
        ("qc.stop_frac_chimeric", "frazione di letture chimeriche, nel report"),
        ("qc.min_frac_phylum", "frazione di varianti con il phylum assegnato, nel report"),
        ("qc.max_frac_contaminant", "frazione di letture rimosse dalla decontaminazione, nel report"),
        ("qc.min_reads_final", "letture per campione nell'oggetto finale, nel report"),
        ("qc.min_frac_reads_retained", "letture trattenute dall'insieme, nel report"),
    ):
        metti(parametro, SCELTA, SENZA_VALORE,
              f"soglia di qualita' da fissare sui propri dati; misura da guardare: {misura}")

    mancanti = [c for c in defaults.OBBLIGATORI if c not in voci]
    assert not mancanti, mancanti

    intestazione = [
        f"Letture: {len(file)} file in {cartella}, ispezionate {quante} di ogni file.",
        f"Metadati: {assay.name}" + (f" e {studio.name}" if studio else "")
        + (f"; file del lotto {lotto.name}" if lotto else "; nessun file del lotto") + ".",
    ]
    return voci, intestazione, non_indovinati, note


def _chiave_del_file(chiave: r.Chiave, nome: str) -> str | None:
    """La chiave del campione a cui appartiene un file."""
    for c, n in chiave.file_per_chiave.items():
        if n == nome:
            return c
    return None


def _chiave_della_riga(chiave: r.Chiave, riga: dict[str, str]) -> str | None:
    """La chiave del campione di una riga della tabella di assay."""
    import re

    from amplicon16s.metadata.crosswalk import AccessionNonEstraibile, chiave_dalla_tabella

    try:
        return chiave_dalla_tabella(riga.get(chiave.colonna, ""), re.compile(chiave.espressione))
    except AccessionNonEstraibile:
        return None


def _letture_appaiate(dettaglio: str) -> str:
    """Il messaggio con cui ``init`` rifiuta un dataset di letture appaiate."""
    return (
        f"le letture risultano appaiate: {dettaglio}. La pipeline tratta dati single-end e "
        "init non propone una configurazione. Usa le sole letture forward: un file per "
        "campione con la prima lettura di ogni coppia, in una cartella a parte, e indica "
        "quella con --reads (GUIDA.md, «Letture appaiate»)"
    )


def esegui(args: argparse.Namespace, stampa: Callable[[str], None] = print) -> int:
    """Esegue ``init`` e restituisce il codice di uscita."""
    uscita: Path = args.output
    if uscita.exists():
        stampa(f"init: {uscita} esiste gia': init non scrive su una configurazione "
               "esistente. Indica con --output un file nuovo")
        return USCITA_RIFIUTO
    if not uscita.parent.is_dir():
        stampa(f"init: la cartella {uscita.parent} non esiste")
        return USCITA_RIFIUTO
    if not os.access(uscita.parent, os.W_OK | os.X_OK):
        # Si verifica prima di leggere le letture: una cartella montata nel
        # container e creata da Docker appartiene a un altro utente.
        stampa(f"init: la cartella {uscita.parent} non e' scrivibile dall'utente che esegue: "
               "creala prima di lanciare il comando, oppure indica un altro --output. "
               "Nulla e' stato scritto")
        return USCITA_RIFIUTO
    try:
        voci, intestazione, non_indovinati, note = prepara(args, lambda t: stampa(f"init: {t}"))
    except Rifiuto as rifiuto:
        stampa(f"init: {rifiuto}. Nulla e' stato scritto")
        return USCITA_RIFIUTO
    scrivi_atomico(uscita, testo_yaml(voci, intestazione, non_indovinati) + "\n")

    da_schema = [
        f"{gruppo}.{nome}" for gruppo, campo in Config.model_fields.items()
        for nome in campo.annotation.model_fields
    ]
    voci = {chiave: voci[chiave] for chiave in da_schema if chiave in voci}
    per_specie = {s: [x for x in voci.values() if x.specie is s] for s in Specie}
    stampa(f"Configurazione di partenza scritta in {uscita}")
    for riga in intestazione:
        stampa(f"  {riga}")
    stampa(f"  rilevati: {len(per_specie[Specie.RILEVATO])}, indicati: "
           f"{len(per_specie[Specie.INDICATO])}, proposti: {len(per_specie[Specie.PROPOSTO])}, "
           f"da scegliere: {len(per_specie[Specie.DA_SCEGLIERE])}")
    stampa("Rilevati sui file (righe attive):")
    for voce in per_specie[Specie.RILEVATO]:
        stampa(f"  {voce.chiave}: {_breve(voce.valore)}")
    stampa("Proposti (righe in commento: togli il commento per accettare):")
    for voce in per_specie[Specie.PROPOSTO]:
        stampa(f"  {voce.chiave}: {_breve(voce.valore)}")
    stampa("Da scegliere:")
    for voce in per_specie[Specie.DA_SCEGLIERE]:
        partenza = "" if voce.valore is SENZA_VALORE else f" (partenza: {_breve(voce.valore)})"
        stampa(f"  {voce.chiave}{partenza}")
    if non_indovinati:
        stampa("Che cosa init non stabilisce:")
        for testo in non_indovinati:
            stampa(f"  - {testo}")
    for testo in note:
        stampa(f"Nota: {testo}")
    stampa(f"Poi: rivedi {uscita}, e lancia validate: elenca i parametri ancora da "
           "dichiarare (E-G15-10).")
    return 0


def _breve(valore: Any) -> str:
    """Un valore su una riga, accorciato se lungo."""
    from amplicon16s.init.scrittura import valore_yaml

    testo = valore_yaml(valore)
    return testo if len(testo) <= 90 else testo[:87] + "..."
