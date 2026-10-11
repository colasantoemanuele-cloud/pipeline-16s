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

**Letture appaiate.** Se ogni file della cartella ha il suo compagno, con una
delle convenzioni dei nomi (``_R1`` e ``_R2``, ``_1`` e ``_2``), ``init``
rileva la modalita' appaiata e propone anche i parametri del gruppo ``pair``:
il primer e il motivo della lettura reverse come per la forward, la lunghezza
dell'inserto misurata sulle coppie, i due troncamenti dal profilo di qualita'
entro il vincolo che la fusione richiede, l'intervallo delle lunghezze fuse.

**Codici di uscita**: ``0`` file scritto; ``3`` rifiuto (file gia' esistente,
file con le due letture di ogni coppia, ingressi non trovati o non leggibili),
nulla e' stato scritto.
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import statistics
import tempfile
from collections import Counter
from collections.abc import Callable, Mapping
from concurrent.futures import ProcessPoolExecutor
from dataclasses import replace
from pathlib import Path
from typing import Any, Final

from amplicon16s.config import defaults
from amplicon16s.config.schema import Config, processori_disponibili
from amplicon16s.errors.exceptions import ErrorePipeline
from amplicon16s.init import rilevamento as r
from amplicon16s.init import vocabolario as v
from amplicon16s.init.scrittura import SENZA_VALORE, Specie, Voce, testo_yaml
from amplicon16s.io_layer.artifacts import AlberoOutput, Fase, scrivi_atomico
from amplicon16s.io_layer.reads import MARCATORI_DI_FILE, abbina_file, espandi_iupac
from amplicon16s.metadata.models import ClasseCampione
from amplicon16s.rbridge.runner import cartella_r, esegui_script
from amplicon16s.soglie.tabella import SOGLIE, non_eseguito, valore_esplorativo
from amplicon16s.report.builder import QUALITA_MINIMA_SUGGERITA, troncamento_suggerito
from amplicon16s.steps.s01_profile import (
    NOME_LUNGHEZZE,
    NOME_QUALITA,
    NOME_SCRIPT,
    NOME_VALORI_QUALITA,
    QUALITA_RAGGRUPPATE,
)

__all__ = ["USCITA_RIFIUTO", "Rifiuto", "accetta", "aggiungi_argomenti", "esegui", "prepara"]

USCITA_RIFIUTO: Final = 3

#: La quota della memoria della macchina che la tabella delle sequenze puo'
#: occupare nel valore di partenza di ``qc.max_asv_count``, e i byte di una
#: cella (un intero di R).
_QUOTA_DI_MEMORIA: Final = 0.1
_BYTE_PER_CELLA: Final = 4

#: Il valore di partenza di qc.head_reads: le letture con cui la frazione
#: misurata su un file (primer, motivo, letture corte) ha, nel caso peggiore di
#: una frazione di 0,5, un intervallo di confidenza al 95% di semiampiezza non
#: superiore a PRECISIONE_DELLE_FRAZIONI. Un punto percentuale e' un quinto del
#: margine con cui il sottocomando soglie propone le soglie in frazione:
#: l'errore di stima resta piccolo rispetto alla tolleranza.
QUANTILE_NORMALE_95: Final = 1.96
#: Le basi di sovrapposizione che i troncamenti proposti lasciano oltre la
#: sovrapposizione minima della fusione: con il predefinito di 12 fanno le 20
#: basi che il tutorial di dada2 chiede oltre la variabilita' di lunghezza
#: della regione ("truncLen must be large enough to maintain 20 +
#: biological.length.variation nucleotides of overlap").
MARGINE_DI_SOVRAPPOSIZIONE: Final = 8

#: Di quanto la lunghezza misurata dell'inserto puo' scostarsi da quella attesa
#: dai primer prima che init lo dica: la stessa quota con cui la fusione
#: confronta la mediana delle lunghezze fuse con l'attesa.
SCARTO_DALL_ATTESA_DEI_PRIMER: Final = 0.05

PRECISIONE_DELLE_FRAZIONI: Final = 0.01
LETTURE_PER_UNA_FRAZIONE: Final = math.ceil(
    (QUANTILE_NORMALE_95 * 0.5 / PRECISIONE_DELLE_FRAZIONI) ** 2)

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
    parser.add_argument("--primer-reverse", type=_sequenza_iupac, metavar="SEQUENZA",
                        help="in modalita' appaiata, il primer reverse usato per "
                             "l'amplificazione, in codici IUPAC e scritto da 5' a 3', se non e' "
                             "fra quelli che init riconosce da se'")
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
    parser.add_argument("--accetta-proposte", action="store_true",
                        help="scrive come righe attive le proposte, i valori di partenza e, "
                             "per le soglie di qualita' senza partenza, il valore che non "
                             "ferma: restano da scrivere le sole scelte senza valore")


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
            # Percorsi assoluti: lo script lavora in una cartella temporanea, e
            # una cartella delle letture indicata in forma relativa non vi si
            # troverebbe.
            {"campioni": {c: str(Path(p).resolve()) for c, p in file_per_chiave.items()},
             "processi": processi},
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
    # L'ultima posizione in cui la qualita' mediana dei campioni non di
    # controllo negativo (la mediana, fra i campioni, della mediana di ciascuno)
    # raggiunge la soglia del troncamento suggerito. Serve alla modalita'
    # appaiata: una lettura reverse puo' avere la prima base sotto la soglia e
    # il resto sopra, e la regola del primo calo non le darebbe un troncamento.
    mediane: dict[int, list[float]] = {}
    for riga in qualita:
        if riga["campione"] not in negativi:
            mediane.setdefault(int(riga["posizione"]), []).append(float(riga["mediana"]))
    buone = [posizione for posizione, valori_q in mediane.items()
             if statistics.median(valori_q) >= QUALITA_MINIMA_SUGGERITA]
    return {
        "ultima_di_qualita": max(buone) if buone else None,
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
# Modalita' appaiata                                                           #
# --------------------------------------------------------------------------- #


def _voci_delle_coppie(
    args: argparse.Namespace,
    metti: Callable[..., None],
    voci: dict[str, Voce],
    note: list[str],
    non_indovinati: list[str],
    avanzamento: Callable[[str], None],
    *,
    file: list[Path],
    nomi: list[str],
    reverse_di: dict[str, Path],
    chiave: r.Chiave | None,
    negativi: set[str],
    ispezioni: dict[str, r.Ispezione],
    cercati: tuple[v.Primer, ...],
    inserti: list[list[int | None]],
    quante: str,
    processi: int,
    tabelle: list[Any],
    taglio_forward: int,
    primer_forward: v.Primer | None,
    profilo_forward: dict[str, Any] | None,
) -> None:
    """Le voci del gruppo ``pair`` per un dataset appaiato, e i troncamenti
    delle due letture entro il vincolo della fusione.

    Segue per la lettura reverse le stesse regole della forward (primer
    dall'elenco, motivo dal consenso, troncamento dal profilo), e in piu'
    misura sulle coppie la lunghezza dell'inserto: da quella vengono il vincolo
    sui troncamenti e l'intervallo delle lunghezze fuse.
    """
    PROPOSTO, RILEVATO, SCELTA, INDICATO = (
        Specie.PROPOSTO, Specie.RILEVATO, Specie.DA_SCEGLIERE, Specie.INDICATO)
    da_considerare = {
        n: i for n, i in ispezioni.items()
        if chiave is None or _chiave_del_file(chiave, n) not in negativi
    } or dict(ispezioni)

    # ---- Primer reverse e taglio ---------------------------------------------
    osservato, parziali = r.primer_nelle_letture(ispezioni, cercati)
    dai_metadati: list[str] = []
    r.primer_dai_metadati(tabelle, dai_metadati)
    primer_noto: str | None = None
    taglio = 0
    if osservato is not None:
        trovato = osservato.primer
        primer_noto, taglio = trovato.sequenza, len(trovato.sequenza)
        if osservato.in_ogni_file and args.head_reads is None:
            specie, limite = RILEVATO, ""
            dove = (f"in testa alle letture reverse di ogni file: dal {osservato.minima:.0%} al "
                    f"{osservato.massima:.0%} delle letture, su {quante}")
        elif not osservato.in_ogni_file:
            sotto = osservato.sotto_soglia
            specie = PROPOSTO
            dove = (f"in testa alle letture reverse nella mediana dei file "
                    f"({osservato.mediana:.0%}), ma non in tutti: dal {osservato.minima:.0%} al "
                    f"{osservato.massima:.0%}")
            limite = (
                f"attenzione: {len(sotto)} file su {len(ispezioni)} non lo portano in testa "
                f"({', '.join(sotto[:5])}{'...' if len(sotto) > 5 else ''}): un taglio fisso "
                "toglierebbe loro basi vere. Vanno capiti prima di accettare")
            note.append(
                f"Il primer reverse {trovato.nome} e' in testa alle letture reverse della "
                f"maggior parte dei file ma non di {len(sotto)}: primer e taglio sono proposti.")
        else:
            specie = PROPOSTO
            dove = f"in testa alle letture reverse di ogni file, su {quante}"
            limite = ("verificato sulle sole prime letture (--head-reads): senza l'opzione "
                      "init lo verifica su tutte")
        metti("pair.primer_sequence_reverse", specie, trovato.sequenza,
              f"primer {trovato.nome} ({trovato.fonte}), {dove}", limite)
        metti("pair.trimLeft_reverse", specie, taglio,
              f"la lunghezza del primer {trovato.nome}, che le letture reverse portano in testa",
              limite)
    else:
        variabile = "; ".join(f"{p.nome} nel {f:.1%}" for p, f in parziali)
        metti("pair.trimLeft_reverse", PROPOSTO, 0,
              f"nessun primer reverse dell'elenco e' in testa alle letture reverse (su "
              f"{quante}); un primer fuori elenco va indicato a mano, con la sua lunghezza qui",
              "primer cercati: " + ", ".join(p.nome for p in cercati),
              (f"attenzione: una parte delle letture reverse inizia con un primer "
               f"({variabile}): e' in posizione variabile, e un taglio fisso non lo toglie")
              if variabile else "")
        if args.primer_reverse is not None:
            primer_noto = args.primer_reverse
            metti("pair.primer_sequence_reverse", INDICATO, args.primer_reverse,
                  "--primer-reverse; non e' in testa alle letture reverse")
        elif len(dai_metadati) == 1:
            primer_noto = dai_metadati[0]
            metti("pair.primer_sequence_reverse", PROPOSTO, primer_noto,
                  "sequenza che i metadati dichiarano come primer reverse; non e' nelle "
                  "letture reverse")
        else:
            nominati = r.primer_nominati(tabelle, v.PRIMER_16S_REVERSE)
            metti("pair.primer_sequence_reverse", SCELTA, None,
                  "partenza: nullo, perche' il primer reverse non e' nelle letture reverse e "
                  "i metadati non ne riportano una sola sequenza. Con il primer nullo la "
                  "validazione non lo cerca. Se e' noto, rilancia init con --primer-reverse "
                  "SEQUENZA: init propone allora anche il motivo delle letture reverse",
                  ("i metadati riportano piu' sequenze con un'etichetta di reverse: "
                   + ", ".join(dai_metadati)) if dai_metadati else "",
                  ("i metadati nominano un primer reverse senza riportarne la sequenza; le "
                   "forme che init conosce con quel nome: "
                   + "; ".join(f"{p.nome} {p.sequenza} ({p.fonte})" for p in nominati))
                  if nominati else "")
            non_indovinati.append(
                "pair.primer_sequence_reverse: il primer reverse non e' nelle letture reverse, "
                "e i metadati non ne riportano una sola sequenza.")

    # ---- Motivo della reverse --------------------------------------------------
    motivo = None
    dissentono: list[str] = []
    if primer_noto is not None:
        dissentono = r.file_che_dissentono(da_considerare, taglio, len(primer_noto))
        motivo = r.motivo_dal_consenso(
            [i for n, i in da_considerare.items() if n not in dissentono],
            taglio, len(primer_noto))
    if dissentono:
        elenco = ", ".join(reverse_di[n].name for n in dissentono[:5]) + (
            "..." if len(dissentono) > 5 else "")
        note.append(
            f"{len(dissentono)} file reverse non iniziano come gli altri e sono esclusi dal "
            f"consenso del motivo delle letture reverse ({elenco}): da guardare (una corsa con "
            "i due file scambiati, un'altra regione, sequenze fuori bersaglio).")
    if motivo is not None:
        metti("pair.conserved_motif_reverse", PROPOSTO, motivo,
              f"il consenso delle {len(primer_noto or '')} basi "
              + ("che seguono il primer reverse" if taglio else
                 "con cui iniziano le letture reverse")
              + f", su non piu' di {r.LETTURE_PER_IL_MOTIVO} letture per file: una base o "
              f"una coppia di basi dove copre almeno il {r.SOGLIA_BASE_DEL_MOTIVO:.0%} delle "
              "letture",
              "serve a riconoscere una corsa con i due file scambiati e una libreria con i "
              "due versi mescolati",
              (f"attenzione: {len(dissentono)} file reverse non iniziano come gli altri e "
               "sono esclusi dal consenso: da guardare prima di eseguire") if dissentono else "")
    else:
        metti("pair.conserved_motif_reverse", SCELTA, None,
              "partenza: nullo. "
              + ("Il consenso delle letture reverse non da' un motivo abbastanza definito"
                 if primer_noto else
                 "Senza conoscere il primer reverse init non propone il motivo con cui "
                 "iniziano le letture reverse")
              + ". Con motivo e primer reverse nulli l'orientamento delle letture reverse "
              "non si verifica")

    # ---- Lunghezza dell'inserto ---------------------------------------------------
    inserto = r.stima_dell_inserto(inserti)
    dei_primer = None
    if primer_forward is not None and osservato is not None:
        if primer_forward.posizione is not None and osservato.primer.posizione is not None:
            dei_primer = osservato.primer.posizione - primer_forward.posizione + 1
    tagli = taglio_forward + taglio
    if inserto is not None:
        confronto = ""
        if dei_primer is not None:
            concordano = abs(inserto.mediana - dei_primer) <= (
                SCARTO_DALL_ATTESA_DEI_PRIMER * dei_primer)
            confronto = (
                f"dai primer riconosciuti nelle letture ({primer_forward.nome} e "
                f"{osservato.primer.nome}) su Escherichia coli se ne attendono {dei_primer}: "
                + ("concordano" if concordano else
                   "NON concordano, e vale la lunghezza misurata sulle letture"))
            if not concordano:
                note.append(
                    f"La lunghezza dell'inserto misurata sulle letture (mediana "
                    f"{inserto.mediana}) non concorda con quella attesa dai primer riconosciuti "
                    f"({dei_primer}): vale la misurata. Verifica la regione amplificata.")
        misura = (
            f"misurata su {inserto.sovrapposte} coppie con le letture sovrapposte, su "
            f"{inserto.coppie} campionate (le prime {r.COPPIE_PER_L_INSERTO} di ogni "
            f"campione): mediana {inserto.mediana}, dal {inserto.minima} al {inserto.massima} "
            "(percentili 0,5 e 99,5)"
            + ("; comprende i primer, che le letture portano" if tagli else ""))
        metti("pair.expected_merged_len", PROPOSTO, inserto.mediana - tagli,
              "la mediana della lunghezza dell'inserto, tolte le basi dei due tagli", misura,
              confronto)
        metti("pair.len_min", PROPOSTO, max(1, inserto.minima - tagli),
              "il percentile 0,5 della lunghezza dell'inserto, tolte le basi dei due tagli",
              "una sequenza fusa piu' corta esce dalla tabella: allarga l'intervallo se la "
              "regione varia di piu' fra i taxa attesi",
              "l'intervallo descrive le coppie misurate, non la regione: se molto sotto la "
              "lunghezza attesa c'e' una popolazione a se' (amplificazioni fuori bersaglio), "
              "il percentile basso la comprende. Dopo la prima esecuzione guarda 05_pair_merge/lunghezze_fuse.tsv "
              "e stringi l'intervallo attorno alla regione")
        metti("pair.len_max", PROPOSTO, max(1, inserto.massima - tagli),
              "il percentile 99,5 della lunghezza dell'inserto, tolte le basi dei due tagli")
    else:
        campionate = sum(len(x) for x in inserti)
        sovrapposte = sum(1 for x in inserti for y in x if y is not None)
        perche = (
            f"la lunghezza dell'inserto non si stima dalle letture: solo {sovrapposte} coppie "
            f"su {campionate} campionate hanno le due letture sovrapposte")
        metti("pair.expected_merged_len", SCELTA, None,
              f"partenza: nullo, perche' {perche}. Senza una lunghezza attesa il confronto "
              "della mediana delle sequenze fuse non si fa")
        metti("pair.len_min", SCELTA, SENZA_VALORE,
              f"{perche}: l'intervallo delle lunghezze fuse va dichiarato conoscendo la "
              "regione amplificata")
        metti("pair.len_max", SCELTA, SENZA_VALORE, "si sceglie insieme a pair.len_min")
        note.append(
            perche[0].upper() + perche[1:] + ". Le letture non coprono l'intera regione "
            "amplificata, e le coppie non si potrebbero fondere: come alternativa, la "
            "modalita' single-end con le sole letture forward (io.layout: \"single-end\", "
            "senza il gruppo pair, e i soli file forward nella cartella).")
        non_indovinati.append(
            "pair.len_min, pair.len_max e i due troncamenti: la lunghezza dell'inserto non si "
            "stima dalle letture.")

    # ---- Profilo della reverse e troncamenti ------------------------------------
    profilo: dict[str, Any] | None = None
    if not args.skip_profile and profilo_forward is not None:
        avanzamento("profilo di qualita' delle letture reverse...")
        leggibili = {
            (_chiave_del_file(chiave, n) if chiave else n): reverse_di[n]
            for n in nomi if ispezioni[n].statistiche.valido
        }
        try:
            profilo = _profilo(leggibili, negativi, processi)
        except ErrorePipeline as guasto:
            note.append(
                f"Profilo di qualita' delle letture reverse non calcolato ({guasto.codice}: "
                f"{guasto.dettaglio or guasto.voce.sintesi}).")
    sovrapposizione = defaults.PAIR_MIN_OVERLAP + MARGINE_DI_SOVRAPPOSIZIONE
    vincolo = (
        "vincolo della fusione: (filter.truncLen - filter.trimLeft) + (pair.truncLen_reverse "
        "- pair.trimLeft_reverse) deve superare la lunghezza della regione senza i primer di "
        f"almeno pair.min_overlap piu' {MARGINE_DI_SOVRAPPOSIZIONE} basi di margine")
    forward = profilo_forward["suggerito"] if profilo_forward else None
    della_reverse = profilo["suggerito"] if profilo else None
    proposti: tuple[int, int] | None = None
    if (inserto is not None and forward and della_reverse
            and forward["senza_perdite"] is not None
            and della_reverse["senza_perdite"] is not None):
        necessarie = inserto.massima + sovrapposizione
        massimo_forward = int(forward["senza_perdite"])
        massimo_reverse = int(della_reverse["senza_perdite"])

        def dalla_qualita(del_profilo: dict[str, Any], massimo: int) -> int:
            # Il troncamento suggerito dal profilo; dove la regola del primo
            # calo non ne da' uno (la prima base sotto la soglia), l'ultima
            # posizione con la qualita' mediana alla soglia; mai oltre la
            # lunghezza che non scarta letture.
            suggerito = del_profilo["suggerito"]["suggerito"]
            if suggerito is None:
                suggerito = del_profilo["ultima_di_qualita"] or 0
            return min(int(suggerito), massimo)

        t_forward = dalla_qualita(profilo_forward, massimo_forward)
        t_reverse = dalla_qualita(profilo, massimo_reverse)
        allungati = False
        if t_forward + t_reverse < necessarie:
            # Prima la forward, che di norma ha la qualita' migliore in coda;
            # poi la reverse; mai oltre la lunghezza che non scarta letture.
            allungati = True
            t_forward = min(massimo_forward, max(t_forward, necessarie - t_reverse))
            t_reverse = min(massimo_reverse, max(t_reverse, necessarie - t_forward))
        if t_forward + t_reverse >= necessarie:
            proposti = (t_forward, t_reverse)
            conto = (
                f"{t_forward} + {t_reverse} = {t_forward + t_reverse} basi, e ne servono "
                f"{necessarie}: l'inserto piu' lungo misurato ({inserto.massima}) piu' "
                f"{sovrapposizione} di sovrapposizione")
            come = (
                "allungato oltre il calo della qualita' quanto serve a fondere le coppie, "
                "senza scartare letture" if allungati else
                "dal profilo di qualita' di questa lettura; non scarta alcuna lettura")
            metti("filter.truncLen", PROPOSTO, t_forward,
                  f"troncamento della lettura forward, {come}", conto, vincolo,
                  "comprende le basi di filter.trimLeft")
            metti("pair.truncLen_reverse", PROPOSTO, t_reverse,
                  f"troncamento della lettura reverse, {come}", conto, vincolo,
                  "comprende le basi di pair.trimLeft_reverse")
            for parametro, troncamento, massimo in (
                ("filter.truncLen_shortfall_warn", t_forward, massimo_forward),
                ("pair.truncLen_shortfall_warn_reverse", t_reverse, massimo_reverse),
            ):
                metti(parametro, SCELTA, min(max(massimo - troncamento, 0), troncamento - 1),
                      "partenza: lo scarto fra la lettura piu' corta e il troncamento "
                      "proposto; con un valore minore il profilo segnalerebbe il troncamento "
                      "proposto (E-S1-01)")
        else:
            note.append(
                f"Nessuna coppia di troncamenti che non scarti letture permette la fusione: "
                f"le due letture intere danno {massimo_forward} + {massimo_reverse} = "
                f"{massimo_forward + massimo_reverse} basi, e ne servono {necessarie} "
                f"(l'inserto piu' lungo misurato, {inserto.massima}, piu' {sovrapposizione} di "
                "sovrapposizione). I due troncamenti restano da scegliere; come alternativa, "
                "la modalita' single-end con le sole letture forward (io.layout: "
                "\"single-end\", senza il gruppo pair, e i soli file forward nella cartella).")
    if proposti is None:
        perche = (
            "la lunghezza dell'inserto non e' stimata" if inserto is None else
            "il profilo di qualita' non e' stato calcolato" if not (forward and della_reverse)
            else "nessuna coppia di troncamenti che non scarti letture rispetta il vincolo")
        for parametro in ("filter.truncLen", "pair.truncLen_reverse"):
            metti(parametro, SCELTA, SENZA_VALORE,
                  f"da scegliere insieme all'altro troncamento: {perche}", vincolo)
        for parametro in ("filter.truncLen_shortfall_warn",
                          "pair.truncLen_shortfall_warn_reverse"):
            metti(parametro, SCELTA, SENZA_VALORE, "si sceglie insieme al suo troncamento")
        metti("qc.max_frac_short_reads", SCELTA, SENZA_VALORE,
              "si sceglie insieme ai troncamenti")

    # ---- Soglie della modalita' appaiata -------------------------------------------
    for parametro, misura in (
        ("pair.max_primer_hit_frac_reverse",
         "frazione di letture reverse con il primer reverse in testa, in "
         "01_input_validation/letture_ispezionate_reverse.tsv"),
        ("pair.min_merged_frac_warn", "quota di coppie fuse, in 05_pair_merge/fusione.json"),
        ("pair.min_merged_frac_stop", "quota di coppie fuse, in 05_pair_merge/fusione.json"),
    ):
        metti(parametro, SCELTA, SENZA_VALORE,
              f"soglia di qualita' da fissare sui propri dati; misura da guardare: {misura}")
    mancanti = [c for c in defaults.PARAMETRI_DELLE_COPPIE if c not in voci]
    assert not mancanti, mancanti


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
    # Se ogni file ha il suo compagno il dataset e' appaiato: da qui in poi
    # ``file`` e ``nomi`` sono i soli file forward, e i reverse li seguono.
    appaiamento = r.coppie_di_file(nomi)
    reverse_di: dict[str, Path] = {}
    if appaiamento is not None:
        marcatori, coppie_di_file = appaiamento
        per_nome = {p.name: p for p in file}
        reverse_di = {f: per_nome[x] for f, x in coppie_di_file.items()}
        file = [per_nome[f] for f in sorted(coppie_di_file)]
        nomi = [p.name for p in file]
        for estensione in list(per_estensione):
            per_estensione[estensione] = [
                p for p in per_estensione[estensione] if p.name in coppie_di_file]
    else:
        inverse = r.secondi_di_coppia(nomi)
        if inverse:
            raise Rifiuto(_senza_compagno(nomi, inverse))

    avanzamento(f"letture: {len(file)} file in {cartella}"
                + (f", ciascuno con il suo file reverse ({marcatori[0]} e {marcatori[1]})"
                   if reverse_di else ""))
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
        # In modalita' appaiata: le letture reverse, con i primer reverse, e la
        # lunghezza dell'inserto su un campione di coppie di ogni campione.
        cercati_reverse = v.PRIMER_16S_REVERSE
        if args.primer_reverse is not None:
            cercati_reverse = (
                v.Primer("indicato", args.primer_reverse, "--primer-reverse"),
                *(p for p in v.PRIMER_16S_REVERSE if p.sequenza != args.primer_reverse))
        ispezioni_reverse: dict[str, r.Ispezione] = {}
        inserti: list[list[int | None]] = []
        if reverse_di:
            primer_reverse = tuple(espandi_iupac(p.sequenza) for p in cercati_reverse)
            ispezioni_reverse = dict(zip(nomi, gruppo.map(
                r.ispeziona,
                [(str(reverse_di[n]), args.head_reads, primer_reverse) for n in nomi])))
            inserti = list(gruppo.map(
                r.inserti_di_una_coppia_di_file,
                [(str(p), str(reverse_di[p.name])) for p in file]))
    if reverse_di:
        # Un file con i segni di entrambe le letture contiene le coppie intere;
        # il marcatore dell'altra lettura dice una corsa con i file scambiati,
        # che init segnala e la validazione ferma.
        intercalati, invertite = [], []
        for n in nomi:
            for ispezione, attesa in ((ispezioni[n], "1"), (ispezioni_reverse[n], "2")):
                esito = (ispezione.statistiche.segni_di_coppia_inattesi(attesa)
                         if ispezione.statistiche.valido else None)
                if esito is not None:
                    (intercalati if esito[0] == "intercalato" else invertite).append(
                        f"{ispezione.statistiche.nome} ({esito[1]})")
        intercalati.sort()
        if invertite:
            note.append(
                f"{len(invertite)} file portano nelle intestazioni il marcatore dell'altra "
                f"lettura ({'; '.join(sorted(invertite)[:3])}"
                f"{'...' if len(invertite) > 3 else ''}): i due file di quelle corse sono "
                "scambiati, e validate si fermera' (E-S0-26) finche' non se ne scambiano i nomi.")
    else:
        intercalati = sorted(
            f"{n} ({i.statistiche.coppie_nello_stesso_file})" for n, i in ispezioni.items()
            if i.statistiche.valido and i.statistiche.coppie_nello_stesso_file)
    if intercalati:
        raise Rifiuto(_letture_appaiate(
            f"{len(intercalati)} file non contengono le sole prime letture delle coppie: "
            + "; ".join(intercalati[:3]) + ("..." if len(intercalati) > 3 else "")
        ) if not reverse_di else (
            f"{len(intercalati)} file contengono le due letture di ogni coppia nello stesso "
            "file: " + "; ".join(intercalati[:3]) + ("..." if len(intercalati) > 3 else "")
            + ". In modalita' appaiata ogni file riporta una sola lettura di ogni coppia: "
            "separa le due letture in due file per campione"))
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
    if reverse_di:
        metti("io.layout", RILEVATO, defaults.LAYOUT_APPAIATO,
              f"ognuno dei {2 * len(file)} file della cartella ha il suo compagno: "
              f"{len(file)} coppie di file che si distinguono per {marcatori[0]} e {marcatori[1]}",
              "per usare le sole letture forward: \"single-end\", togliendo il gruppo pair e "
              "lasciando nella cartella (o in io.fastq_glob) i soli file forward")
        metti("pair.file_markers", RILEVATO, list(marcatori),
              "i due marcatori con cui i nomi dei file di ogni coppia si distinguono, prima "
              "quello del file forward",
              "abbinamento, per esempio: " + "; ".join(
                  f"{n} e {reverse_di[n].name}" for n in nomi[:2]))
    else:
        metti("io.layout", RILEVATO, defaults.LAYOUT_SINGOLO,
              "nessun file della cartella ha un compagno che se ne distingua per un "
              "marcatore di coppia (R1 e R2, 1 e 2): un file per campione")
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
            metti("qc.primer_sequence", SCELTA, None,
                  "partenza: nullo, perche' il primer forward non e' nelle letture e i "
                  "metadati non ne riportano la sequenza. Con il primer nullo la validazione "
                  "non cerca alcun primer in testa alle letture. Se il primer e' noto, "
                  "rilancia init con --primer SEQUENZA: init propone allora anche il motivo "
                  "conservato",
                  ("i metadati riportano sequenze senza dire quale sia il forward: "
                   + ", ".join(sequenze)) if sequenze else "",
                  ("i metadati nominano un primer senza riportarne la sequenza; le forme che "
                   "init conosce con quel nome: "
                   + "; ".join(f"{p.nome} {p.sequenza} ({p.fonte})" for p in nominati)
                   + ". Come scegliere fra piu' forme: GUIDA.md, sezione 6.3")
                  if nominati else "")
            non_indovinati.append(
                "qc.primer_sequence: il primer forward non e' nelle letture, e i metadati "
                + ("lo nominano senza riportarne la sequenza." if nominati else
                   "non ne riportano la sequenza."))
    motivo = None
    dissentono: list[str] = []
    if primer_noto is not None:
        # I file che non iniziano come la maggioranza escono dal consenso: un
        # file di letture inverse lo allargherebbe fino a comprendere se stesso.
        dissentono = r.file_che_dissentono(
            {n: i for n, i in ispezioni.items() if any(i is c for c in da_considerare)},
            taglio, len(primer_noto))
        motivo = r.motivo_dal_consenso(
            [i for n, i in ispezioni.items()
             if any(i is c for c in da_considerare) and n not in dissentono],
            taglio, len(primer_noto))
    elenco_diversi = ", ".join(dissentono[:5]) + ("..." if len(dissentono) > 5 else "")
    quanti_diversi = (
        "1 file non inizia come gli altri ed e' escluso" if len(dissentono) == 1 else
        f"{len(dissentono)} file non iniziano come gli altri e sono esclusi")
    diversi = (
        f"attenzione: {quanti_diversi} dal consenso ({elenco_diversi}): da guardare prima "
        "di eseguire. Puo' contenere le letture inverse di una corsa appaiata, un'altra "
        "regione, o un campione dominato da sequenze fuori bersaglio"
    ) if dissentono else ""
    if dissentono:
        note.append(
            f"{quanti_diversi} dal consenso del motivo ({elenco_diversi}): da guardare "
            "(letture inverse, un'altra regione, sequenze fuori bersaglio).")
    if motivo is not None:
        metti("qc.conserved_motif", PROPOSTO, motivo,
              f"il consenso delle {len(primer_noto or '')} basi "
              + ("che seguono il primer" if taglio else "con cui iniziano le letture")
              + f", su non piu' di {r.LETTURE_PER_IL_MOTIVO} letture per file: una base o "
              f"una coppia di basi dove copre almeno il {r.SOGLIA_BASE_DEL_MOTIVO:.0%} delle letture",
              "validate ne misura la frequenza con qc.min_motif_frac", senza_negativi, diversi)
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

    # ---- Modalita' appaiata: la lettura reverse e la fusione ----------------------
    if reverse_di:
        _voci_delle_coppie(
            args, metti, voci, note, non_indovinati, avanzamento,
            file=file, nomi=nomi, reverse_di=reverse_di, chiave=chiave,
            negativi=negativi_proposti, ispezioni=ispezioni_reverse, cercati=cercati_reverse,
            inserti=inserti, quante=quante, processi=processi,
            tabelle=[righe_assay] + ([righe_studio] if righe_studio is not None else []),
            taglio_forward=taglio, primer_forward=osservato.primer if osservato else None,
            profilo_forward=profilo,
        )

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
    # Con ogni file letto per intero il file piu' grande limita il valore: uno
    # maggiore leggerebbe comunque tutte le letture di ogni file.
    intere = args.head_reads is None and bool(letture)
    ispezionate = min(LETTURE_PER_UNA_FRAZIONE, max(letture)) if intere else LETTURE_PER_UNA_FRAZIONE
    metti("qc.head_reads", SCELTA, ispezionate,
          "partenza: "
          + (f"le letture del file piu' grande ({max(letture)}), meno delle "
             f"{LETTURE_PER_UNA_FRAZIONE} con cui" if intere and ispezionate < LETTURE_PER_UNA_FRAZIONE
             else "le letture con cui")
          + " una frazione misurata su un file ha un intervallo di confidenza al "
          f"95% di non oltre {PRECISIONE_DELLE_FRAZIONI:.0%} in piu' o in meno; non e' una scelta "
          "scientifica: con piu' letture la validazione dura di piu' e le frazioni sono "
          "piu' precise",
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
    if voci["ctrl.blank_values"].valore:
        non_indovinati.append("decontam.threshold: la soglia di decontaminazione.")
    non_indovinati += [
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

    mancanti = [c for c in (*defaults.OBBLIGATORI, *defaults.CONDIZIONALI) if c not in voci]
    assert not mancanti, mancanti

    # I condizionali la cui condizione e' falsa con i valori rilevati, proposti
    # o di partenza: non si chiede di sceglierli. Un parametro senza valore
    # (da scegliere, o gia' non richiesto) vale come la sua assenza, e la
    # condizione che ne dipende e' falsa: validate chiedera' il condizionale se
    # l'utente la rende vera. L'ordine dell'elenco mette prima le etichette, da
    # cui dipendono i parametri dei controlli positivi.
    for chiave in defaults.CONDIZIONALI:
        valori = {c: v.valore for c, v in voci.items() if v.valore is not SENZA_VALORE}
        if not defaults.condizione_vera(chiave, valori):
            voci[chiave] = Voce(
                chiave, Specie.CONDIZIONALE, SENZA_VALORE,
                f"serve solo se {defaults.testo_condizione(chiave)}: con i valori di "
                "questo file non e' richiesto e si puo' lasciare in commento")

    intestazione = [
        (f"Letture appaiate: {len(file)} coppie di file in {cartella}, ispezionate {quante} "
         "di ogni file." if reverse_di else
         f"Letture: {len(file)} file in {cartella}, ispezionate {quante} di ogni file."),
        f"Metadati: {assay.name}" + (f" e {studio.name}" if studio else "")
        + (f"; file del lotto {lotto.name}" if lotto else "; nessun file del lotto") + ".",
    ]
    return voci, intestazione, non_indovinati, note


def accetta(voci: Mapping[str, Voce]) -> dict[str, Voce]:
    """Le voci con proposte e valori di partenza accettati in blocco.

    E' il gesto di ``--accetta-proposte``: ogni proposta diventa una riga
    attiva, ogni scelta con un valore di partenza lo prende (``[PARTENZA]``), e
    ogni soglia che ferma o avvisa rimasta senza valore prende quello che non
    ferma (``[ESPLORATIVA]``), se ne ha uno. Nessun valore e' inventato: una
    scelta senza valore di partenza resta da scegliere, un condizionale non
    pertinente resta omesso.
    """
    accettate = {}
    valori = {c: x.valore for c, x in voci.items() if x.valore is not SENZA_VALORE}
    esplorative = {s.chiave: s for s in SOGLIE}
    for chiave, voce in voci.items():
        if voce.specie is PROPOSTO:
            voce = replace(voce, accettata=True,
                           origine=voce.origine + "; accettato con --accetta-proposte")
        elif voce.specie is SCELTA and voce.valore is not SENZA_VALORE:
            voce = replace(voce, specie=Specie.PARTENZA,
                           origine=voce.origine.removeprefix("partenza: ")
                           + "; accettato con --accetta-proposte")
        elif voce.specie is SCELTA and chiave in esplorative:
            soglia = esplorative[chiave]
            valore = valore_esplorativo(soglia)
            if valore is not None:
                senza_effetto = non_eseguito(soglia, {**valori, chiave: valore})
                voce = replace(
                    voce, specie=Specie.ESPLORATIVA, valore=valore,
                    origine=(
                        "il valore che non ferma: il controllo "
                        f"({', '.join(soglia.codici)}) e' spento. Serve alla prima "
                        "esecuzione; poi il sottocomando soglie propone il valore dalla "
                        "misura" if senza_effetto is None else
                        f"il valore che non ferma; {senza_effetto}"),
                    note=(voce.origine, *voce.note))
        accettate[chiave] = voce
    return accettate


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


def _senza_compagno(nomi: list[str], inverse: list[str]) -> str:
    """Il messaggio con cui ``init`` rifiuta una cartella in cui alcuni file
    hanno il nome di una seconda lettura ma non ogni file ha il suo compagno:
    nomina i file rimasti soli e le due strade possibili.
    """
    soli: list[str] = []
    for marcatori in MARCATORI_DI_FILE:
        abbinamento = abbina_file(nomi, marcatori)
        if abbinamento.coppie:
            soli = [*abbinamento.senza_compagno, *abbinamento.ambigui]
            break
    elenco = soli or inverse
    return (
        f"{len(inverse)} file hanno il nome della seconda lettura di una coppia, ma non ogni "
        f"file della cartella ha il suo compagno: {', '.join(elenco[:5])}"
        f"{'...' if len(elenco) > 5 else ''}"
        + (" (senza compagno)" if soli else "")
        + ". Per la modalita' appaiata aggiungi i file che mancano o togli dalla cartella "
        "quelli che non appartengono a una coppia; per usare le sole letture forward tieni "
        "nella cartella indicata con --reads un file per campione, con la prima lettura di "
        "ogni coppia (GUIDA.md, «Letture appaiate»)"
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
    accettate = bool(getattr(args, "accetta_proposte", False))
    if accettate:
        voci = accetta(voci)
    scrivi_atomico(uscita, testo_yaml(voci, intestazione, non_indovinati, accettate) + "\n")

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
           f"da scegliere: {len(per_specie[Specie.DA_SCEGLIERE])}, "
           f"non richiesti: {len(per_specie[Specie.CONDIZIONALE])}")
    stampa("Rilevati sui file (righe attive):")
    for voce in per_specie[Specie.RILEVATO]:
        stampa(f"  {voce.chiave}: {_breve(voce.valore)}")
    stampa("Proposti (righe attive, accettati con --accetta-proposte: controllali):"
           if accettate else
           "Proposti (righe in commento: togli il commento per accettare):")
    for voce in per_specie[Specie.PROPOSTO]:
        stampa(f"  {voce.chiave}: {_breve(voce.valore)}")
    if accettate:
        stampa("Valori di partenza accettati (righe attive):")
        for voce in per_specie[Specie.PARTENZA]:
            stampa(f"  {voce.chiave}: {_breve(voce.valore)}")
        assenze = [x.chiave for x in per_specie[Specie.PARTENZA] if x.valore in (None, [])]
        if assenze:
            stampa("  Fra questi dichiarano un'assenza (nullo o elenco vuoto), da rivedere se "
                   "lo studio ha quell'informazione: " + ", ".join(assenze))
        stampa("Soglie esplorative, al valore che non ferma (righe attive):")
        for voce in per_specie[Specie.ESPLORATIVA]:
            stampa(f"  {voce.chiave}: {_breve(voce.valore)}")
    stampa("Restano da scegliere (righe in commento, da scrivere):" if accettate else
           "Da scegliere:")
    for voce in per_specie[Specie.DA_SCEGLIERE]:
        partenza = "" if voce.valore is SENZA_VALORE else f" (partenza: {_breve(voce.valore)})"
        stampa(f"  {voce.chiave}{partenza}")
    if per_specie[Specie.CONDIZIONALE]:
        stampa("Non richiesti con le proposte attuali (condizionali):")
        for voce in per_specie[Specie.CONDIZIONALE]:
            stampa(f"  {voce.chiave}: {voce.origine.split(':')[0]}")
    if non_indovinati:
        stampa("Che cosa init non stabilisce:")
        for testo in non_indovinati:
            stampa(f"  - {testo}")
    for testo in note:
        stampa(f"Nota: {testo}")
    stampa(f"Poi: rivedi {uscita}, e lancia validate: elenca i parametri ancora da "
           "dichiarare (E-G15-10).")
    if accettate:
        stampa("Dopo la prima esecuzione, il sottocomando soglie propone dalle misure le "
               "soglie ora esplorative.")
    return 0


def _breve(valore: Any) -> str:
    """Un valore su una riga, accorciato se lungo."""
    from amplicon16s.init.scrittura import valore_yaml

    testo = valore_yaml(valore)
    return testo if len(testo) <= 90 else testo[:87] + "..."
