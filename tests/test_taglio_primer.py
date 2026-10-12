"""Il taglio dei primer cercati lungo la lettura (``filter.primer_removal: search``).

Verifica la regola di ricerca (codici IUPAC nel primer, basi della lettura alla
lettera, sole sostituzioni, finestra, scelta fra piu' corrispondenze), che il
taglio a posizione fissa non veda quello cercato (stesso digest, stessa
configurazione registrata, parametri respinti fuori dalla loro modalita'), che
validazione, profilo e filtro misurino ogni lunghezza dopo il primer trovato,
che le due letture di una coppia restino abbinate, e che ``init`` proponga la
ricerca solo quando il primer non sta a una posizione fissa.

I test che eseguono il filtro richiedono R con dada2 e ShortRead: senza si
saltano, e nell'immagine dell'integrazione continua devono girare.

    docker run --rm -u "$(id -u):$(id -g)" -e HOME=/tmp -e PYTHONPATH=/app/src \\
      -e AMPLICON16S_R_DIR=/app/R -v "$(pwd)":/app -w /app <immagine> \\
      pytest tests/test_taglio_primer.py
"""

from __future__ import annotations

import copy
import csv
import gzip
import json
import os
import random
from collections import Counter
from pathlib import Path

import pytest
import sintetico
import yaml
from init_aiuti import righe_dei_parametri
from sottoinsieme import motivo_pacchetti_r_assenti

from amplicon16s.cli import main
from amplicon16s.config import defaults
from amplicon16s.config.resolve import parametri_dichiarati, risolvi, scrivi_risolta
from amplicon16s.config.schema import ErroreConfigurazione, valida
from amplicon16s.errors.catalog import CATALOGO, Categoria
from amplicon16s.gates.g01_g15 import rifiuto_di_g15
from amplicon16s.init import rilevamento
from amplicon16s.io_layer.reads import RicercaPrimer, apri_fastq, scansiona_file
from amplicon16s.logging.logger import chiudi
from amplicon16s.runner.executor import Conclusione, Esecutore
from amplicon16s.runner.graph import Passo
from amplicon16s.runner.project import ProjectRun
from amplicon16s.runner.tracciamento import componi
from amplicon16s.soglie import tabella

#: Il primer forward dei dataset di questo modulo: il 515F di Parada et al.
#: 2016, con due posizioni degenerate.
PRIMER = "GTGYCAGCMGCCGCGGTAA"

#: La regola di ricerca dei test: una discordanza, fino a sette basi prima del
#: primer.
CERCA = {"primer_removal": "search", "primer_max_mismatch": 1, "primer_max_offset": 7,
         "primer_untrimmed": "discard", "trimLeft": 0}


@pytest.fixture(autouse=True)
def uscite_pulite():
    """Chiude le uscite del log prima e dopo ogni test."""
    chiudi()
    yield
    chiudi()


@pytest.fixture
def con_r():
    """Richiede R con i pacchetti del filtro: salta senza, ma in CI fallisce."""
    motivo = motivo_pacchetti_r_assenti("dada2", "ShortRead", "Biostrings", "jsonlite")
    if motivo is not None:
        if os.environ.get("AMPLICON16S_RICHIEDI_BIOC") == "1":
            pytest.fail(f"Bioconductor e' richiesto in questo ambiente: {motivo}")
        pytest.skip(motivo)


def _letture(percorso: Path) -> dict[str, str]:
    """Identificativo -> sequenza di ogni lettura di un FASTQ."""
    with apri_fastq(percorso) as file:
        righe = file.read().splitlines()
    return {righe[i][1:].split()[0]: righe[i + 1] for i in range(0, len(righe), 4)}


def _ordine(percorso: Path) -> list[str]:
    """Gli identificativi delle letture di un FASTQ, nell'ordine del file."""
    with apri_fastq(percorso) as file:
        return [riga[1:].split()[0] for numero, riga in enumerate(file) if numero % 4 == 0]


def _tabella(percorso: Path) -> list[dict[str, str]]:
    with open(percorso, encoding="utf-8", newline="") as file:
        return list(csv.DictReader(file, delimiter="\t"))


def _cercato(dataset, uscita: Path, **gruppi) -> dict:
    """La configurazione del dataset sintetico con il primer cercato: il
    troncamento e' la regione amplificata, senza il primer.
    """
    dati = sintetico.configurazione(dataset, uscita)
    dati["filter"].update(CERCA, truncLen=137)
    dati["qc"]["max_frac_untrimmed"] = 0.5
    for gruppo, valori in gruppi.items():
        dati.setdefault(gruppo, {}).update(valori)
    return dati


def _esegui(dati: dict, fino_a: Passo):
    """Esegue la catena fino alla fase indicata: esito ed esecuzione."""
    run = ProjectRun(valida(dati))
    esito = Esecutore(run, fino_a=fino_a).esegui()
    chiudi()
    return esito, run


# --------------------------------------------------------------------------- #
# La regola di ricerca                                                         #
# --------------------------------------------------------------------------- #


def test_la_regola_di_ricerca_taglia_alla_base_giusta():
    """
    **Obiettivo**: Verificare su letture costruite che il primer preceduto da
    zero a sette basi casuali sia trovato alla sua posizione, con e senza una
    discordanza; che una N della lettura conti come discordanza anche contro
    una N del primer; che il primer oltre la finestra, a cavallo della fine
    della lettura o con una discordanza di troppo non sia trovato; che fra due
    corrispondenze valga quella con meno discordanze e, a parita', la piu' a
    sinistra; e che non si ammettano inserzioni ne' delezioni.

    **Razionale scientifico e sistemistico**: Ogni lettura deve cominciare,
    dopo il taglio, alla stessa base della regione amplificata: una lettura
    tagliata una base prima o dopo e' una variante diversa dello stesso taxon.
    Una corrispondenza accettata con troppe discordanze taglierebbe basi vere
    a una lettura che il primer non lo porta.
    """
    caso = random.Random(7)
    ricerca = RicercaPrimer(PRIMER, 1, 7)
    regione = "TACGGAGGGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGCACGCAGG"
    concreto = "GTGCCAGCAGCCGCGGTAA"
    assert len(concreto) == len(PRIMER)
    for sfasamento in range(8):
        distanziatore = "".join(caso.choice("ACGT") for _ in range(sfasamento))
        # Le due posizioni degenerate ammettono entrambe le basi.
        for primer in (concreto, "GTGTCAGCCGCCGCGGTAA"):
            lettura = distanziatore + primer + regione
            trovata = ricerca.trova(lettura)
            assert trovata is not None and trovata[1] == 0
            assert lettura[ricerca.fine(lettura):] == regione
        # Una discordanza in una posizione qualunque del primer.
        posizione = caso.randrange(len(concreto))
        sbagliato = (concreto[:posizione] + ("A" if concreto[posizione] != "A" else "C")
                     + concreto[posizione + 1:])
        if RicercaPrimer(PRIMER, 0, 0).trova(sbagliato) is None:
            lettura = distanziatore + sbagliato + regione
            assert lettura[ricerca.fine(lettura):] == regione
            assert ricerca.trova(lettura)[1] == 1
        # Una N nella lettura e' una discordanza: una si ammette, due no.
        con_n = distanziatore + "N" + concreto[1:] + regione
        assert ricerca.trova(con_n) == (sfasamento, 1)
        assert ricerca.trova(distanziatore + "NN" + concreto[2:] + regione) is None
        assert RicercaPrimer(PRIMER, 0, 7).trova(con_n) is None

    # Una N della lettura non concorda nemmeno con una N del primer.
    con_la_n = RicercaPrimer("ACGNACGT", 0, 0)
    assert con_la_n.trova("ACGTACGTAA") == (0, 0)
    assert con_la_n.trova("ACGNACGTAA") is None
    assert RicercaPrimer("ACGNACGT", 1, 0).trova("ACGNACGTAA") == (0, 1)

    # Oltre la finestra, a cavallo della fine, due discordanze.
    assert ricerca.trova("ACGTACGT" + concreto + regione) is None
    assert RicercaPrimer(PRIMER, 1, 8).trova("ACGTACGT" + concreto + regione) == (8, 0)
    assert ricerca.trova("ACG" + concreto[:-1]) is None
    assert ricerca.trova("ACG" + concreto) == (3, 0)
    assert ricerca.trova("ATACCAGCAGCCGCGGTAA" + regione) is None
    # Nessuna inserzione ne' delezione: una base in piu' o in meno dentro il
    # primer sposta tutto il resto, e le discordanze sono troppe.
    assert ricerca.trova(concreto[:9] + "T" + concreto[9:] + regione) is None
    assert ricerca.trova(concreto[:9] + concreto[10:] + regione) is None

    # Due corrispondenze nella finestra: meno discordanze, poi piu' a sinistra.
    breve = RicercaPrimer("ACGTAC", 1, 8)
    assert breve.trova("ACGAACACGTACTT") == (6, 0)      # una a 0 con una discordanza
    assert breve.trova("ACGTACACGTACTT") == (0, 0)      # due esatte: la prima
    assert breve.trova("ACGAACTACGTTCTT") == (0, 1)     # due con una discordanza: la prima
    # Le basi si confrontano in maiuscolo; un codice di ambiguita' nella
    # lettura e' una discordanza come la N.
    assert ricerca.trova("ac" + concreto.lower() + regione.lower()) == (2, 0)
    assert ricerca.trova("GTGRCAGCAGCCGCGGTAA" + regione) == (0, 1)
    assert RicercaPrimer("ACGNACGT", 0, 0).trova("ACGRACGTAA") is None
    # Un codice sconosciuto nel primer e' un errore, non una corrispondenza.
    with pytest.raises(ValueError, match="IUPAC"):
        RicercaPrimer("ACGX", 0, 0)


def test_la_scansione_misura_le_lunghezze_dopo_il_primer(tmp_path):
    """
    **Obiettivo**: Verificare che con la ricerca del primer la scansione di un
    file conti le letture in cui il primer si trova, cerchi il motivo dopo il
    primer trovato, e misuri lunghezza minima, massima e letture corte sulla
    lettura tagliata; e che una lettura senza primer entri nelle lunghezze solo
    se il filtro la terra'.

    **Razionale scientifico e sistemistico**: Con il primer cercato il
    troncamento e' la lunghezza dopo il taglio: confrontarlo con la lettura
    intera farebbe passare letture che il filtro scartera', e fermare un
    troncamento corretto.
    """
    regione = "TACGGAGGGTGCAAGCGTTAATCGGAATTAC"
    primer = "GTGCCAGCAGCCGCGGTAA"
    letture = [
        primer + regione,                    # sfasamento 0: restano 31 basi
        "ACGT" + primer + regione,           # sfasamento 4: restano 31
        "AC" + primer + regione[:20],        # restano 20: corta
        "T" * 60,                            # senza primer
    ]
    percorso = tmp_path / "a.fastq"
    percorso.write_text("".join(f"@r{i}\n{s}\n+\n{'I' * len(s)}\n" for i, s in enumerate(letture)))
    scarta = scansiona_file(percorso, 100, None, "^TACGGAGG", lunghezza_richiesta=31,
                            ricerca=RicercaPrimer(PRIMER, 1, 7))
    assert (scarta.letture_esaminate, scarta.frazione_primer, scarta.frazione_motivo) == (
        4, 0.75, 0.75)
    assert (scarta.lunghezza_minima, scarta.lunghezza_massima, scarta.piu_corte) == (20, 31, 1)
    tiene = scansiona_file(percorso, 100, None, "^TACGGAGG", lunghezza_richiesta=31,
                           ricerca=RicercaPrimer(PRIMER, 1, 7, scarta_senza_primer=False))
    assert (tiene.lunghezza_minima, tiene.lunghezza_massima, tiene.piu_corte) == (20, 60, 1)
    # Senza la ricerca nulla cambia: il motivo si cerca alla posizione data.
    fisso = scansiona_file(percorso, 100, None, "^TACGGAGG", inizio_motivo=19,
                           lunghezza_richiesta=31)
    assert (fisso.frazione_motivo, fisso.lunghezza_minima, fisso.piu_corte) == (0.25, 41, 0)


# --------------------------------------------------------------------------- #
# La configurazione                                                            #
# --------------------------------------------------------------------------- #


def test_il_taglio_a_posizione_fissa_dichiarato_vale_come_omesso(tmp_path):
    """
    **Obiettivo**: Verificare che una configurazione che dichiara
    ``filter.primer_removal: fixed`` abbia lo stesso digest, la stessa impronta
    dei risultati, gli stessi parametri dichiarati e lo stesso
    ``resolved.yaml``, byte per byte, di quella che non lo dichiara; che la
    configurazione registrata non nomini il modo ne' alcun parametro del taglio
    cercato; e che con ``search`` li riporti tutti.

    **Razionale scientifico e sistemistico**: Il taglio cercato non deve
    cambiare alcun artefatto di un'esecuzione che non lo usa: una configurazione
    esistente conserva digest e impronte delle fasi, e nessuna fase si rifa'.
    """
    dataset = sintetico.genera(tmp_path / "dati", biologici=2, letture=30, primer=PRIMER)
    omesso = sintetico.configurazione(dataset, tmp_path / "out")
    dichiarato = copy.deepcopy(omesso)
    dichiarato["filter"]["primer_removal"] = "fixed"
    a, b = risolvi(valida(omesso)), risolvi(valida(dichiarato))
    assert not a.config.taglio_cercato and not b.config.taglio_cercato
    assert a.digest == b.digest and a.impronta_risultati == b.impronta_risultati
    assert parametri_dichiarati(a.config) == parametri_dichiarati(b.config)
    testi = [scrivi_risolta(r, tmp_path / nome).read_text(encoding="utf-8")
             for nome, r in (("a", a), ("b", b))]
    assert testi[0] == testi[1]
    assenti = ("primer_removal", "primer_max_mismatch", "primer_max_offset",
               "primer_untrimmed", "max_frac_untrimmed")
    assert not any(chiave in testi[0] for chiave in assenti)
    mappa = a.come_mappa()
    assert not set(assenti) & (set(mappa["filter"]) | set(mappa["qc"]))

    cercato = risolvi(valida(_cercato(dataset, tmp_path / "out")))
    assert cercato.config.taglio_cercato and cercato.digest != a.digest
    assert cercato.come_mappa()["filter"]["primer_removal"] == "search"
    assert cercato.come_mappa()["qc"]["max_frac_untrimmed"] == 0.5
    assert {"filter.primer_removal", *defaults.PARAMETRI_DEL_TAGLIO_CERCATO} <= set(
        parametri_dichiarati(cercato.config))
    # Le lunghezze attese delle varianti: il troncamento, che e' gia' senza primer.
    assert (cercato.derivati.asv_len_min, cercato.derivati.asv_len_max) == (137, 137)


def test_i_parametri_del_taglio_cercato_valgono_solo_con_la_ricerca(tmp_path):
    """
    **Obiettivo**: Verificare che un parametro del taglio cercato dichiarato
    senza ``filter.primer_removal: search`` sia respinto con ``E-G15-16``, che
    lo nomina; che con ``search`` un taglio fisso diverso da zero e un primer
    nullo siano respinti con lo stesso codice, ciascuno con la sua condizione;
    che in modalita' appaiata valgano le stesse condizioni per la lettura
    reverse; e che i parametri mancanti siano elencati da ``E-G15-10`` con la
    loro condizione.

    **Razionale scientifico e sistemistico**: Un parametro accettato senza
    effetto farebbe credere che il primer sia cercato; un taglio fisso insieme
    alla ricerca toglierebbe basi due volte; senza la sequenza del primer non
    c'e' nulla da cercare.
    """
    dataset = sintetico.genera(tmp_path / "dati", biologici=2, letture=30, primer=PRIMER)

    def rifiuto(dati: dict) -> list[tuple[str, str]]:
        with pytest.raises(ErroreConfigurazione) as respinta:
            valida(dati)
        return [(v.codice, v.dettaglio) for v in rifiuto_di_g15(respinta.value).violazioni]

    for chiave, valore, gruppo in (
        ("primer_max_mismatch", 1, "filter"), ("primer_max_offset", 7, "filter"),
        ("primer_untrimmed", "keep", "filter"), ("max_frac_untrimmed", 0.1, "qc"),
    ):
        fisso = sintetico.configurazione(dataset, tmp_path / "out")
        fisso[gruppo][chiave] = valore
        ((codice, dettaglio),) = rifiuto(fisso)
        assert codice == "E-G15-16" and f"{gruppo}.{chiave}" in dettaglio
        assert "filter.primer_removal" in dettaglio

    ((codice, dettaglio),) = rifiuto(_cercato(dataset, tmp_path / "out", filter={"trimLeft": 19}))
    assert codice == "E-G15-16" and "filter.trimLeft deve essere 0" in dettaglio
    assert "qc.primer_sequence" not in dettaglio
    ((codice, dettaglio),) = rifiuto(
        _cercato(dataset, tmp_path / "out", qc={"primer_sequence": None}))
    assert codice == "E-G15-16" and "qc.primer_sequence non puo' essere nullo" in dettaglio
    assert "filter.trimLeft" not in dettaglio
    ((codice, dettaglio),) = rifiuto(_cercato(
        dataset, tmp_path / "out", filter={"trimLeft": 3}, qc={"primer_sequence": None}))
    assert "filter.trimLeft deve essere 0" in dettaglio and "qc.primer_sequence" in dettaglio

    # Tante discordanze quante le basi del primer: ogni lettura lo porterebbe.
    ((codice, dettaglio),) = rifiuto(_cercato(
        dataset, tmp_path / "out", filter={"primer_max_mismatch": len(PRIMER)}))
    assert codice == "E-G15-16" and "filter.primer_max_mismatch" in dettaglio

    mancanti = _cercato(dataset, tmp_path / "out")
    del mancanti["filter"]["primer_max_offset"], mancanti["qc"]["max_frac_untrimmed"]
    ((codice, dettaglio),) = rifiuto(mancanti)
    assert codice == "E-G15-10"
    assert "filter.primer_max_offset" in dettaglio and "qc.max_frac_untrimmed" in dettaglio
    assert "perche' filter.primer_removal e' search" in dettaglio
    # Un valore non ammesso resta un problema di schema, con il suo parametro.
    with pytest.raises(ErroreConfigurazione, match="filter.primer_untrimmed"):
        valida(_cercato(dataset, tmp_path / "out", filter={"primer_untrimmed": "tieni"}))

    coppie = sintetico.genera_coppie(
        tmp_path / "coppie", biologici=2, coppie=20, primer=PRIMER,
        primer_reverse=sintetico.PRIMER_805R)
    appaiata = sintetico.configurazione_appaiata(coppie, tmp_path / "out2")
    appaiata["filter"].update(CERCA)
    appaiata["qc"]["max_frac_untrimmed"] = 0.5
    ((codice, dettaglio),) = rifiuto(appaiata)
    assert codice == "E-G15-16" and "pair.trimLeft_reverse deve essere 0" in dettaglio
    appaiata["pair"].update(trimLeft_reverse=0, primer_sequence_reverse=None)
    ((codice, dettaglio),) = rifiuto(appaiata)
    assert codice == "E-G15-16" and "pair.primer_sequence_reverse non puo' essere nullo" in dettaglio
    appaiata["pair"]["primer_sequence_reverse"] = sintetico.PRIMER_805R
    assert valida(appaiata).taglio_cercato
    assert CATALOGO["E-G15-16"].categoria is Categoria.REVISIONE_UMANA


def test_la_soglia_delle_letture_senza_primer_sta_fra_le_soglie(tmp_path):
    """
    **Obiettivo**: Verificare che ``qc.max_frac_untrimmed`` sia una soglia del
    meccanismo comune solo con il primer cercato: una frazione massima che da'
    l'avviso ``E-S2-04``, con il margine delle frazioni, il valore che non
    ferma a uno e il controllo spento al limite; che con il taglio a posizione
    fissa non sia fra le soglie di cui si parla; e che con il primer cercato la
    soglia del primer in testa risulti senza effetto.

    **Razionale scientifico e sistemistico**: Una soglia fuori dal meccanismo
    comune non avrebbe valore esplorativo, non sarebbe proposta da ``soglie`` e
    il suo controllo spento non sarebbe dichiarato. Con il primer cercato il
    primer deve esserci: giudicarne la presenza come un difetto fermerebbe
    un'esecuzione corretta.
    """
    dataset = sintetico.genera(tmp_path / "dati", biologici=2, letture=30, primer=PRIMER)
    voce = tabella.soglia("qc.max_frac_untrimmed")
    assert (voce.direzione, voce.unita, voce.codici, voce.ferma, voce.margine) == (
        tabella.Direzione.MASSIMA, tabella.Unita.FRAZIONE, ("E-S2-04",), False,
        tabella.MARGINE_FRAZIONE)
    assert CATALOGO["E-S2-04"].categoria is Categoria.DEGRADAZIONE_AUTOMATICA
    assert tabella.valore_esplorativo(voce) == 1.0

    fisso = tabella.valori_piatti(valida(sintetico.configurazione(dataset, tmp_path / "out")))
    assert "qc.max_frac_untrimmed" not in {s.chiave for s in tabella.soglie_di(fisso)}
    assert not [s for s in tabella.controlli_spenti(fisso, 2)
                if s.chiave == "qc.max_frac_untrimmed"]

    cercato = tabella.valori_piatti(valida(_cercato(dataset, tmp_path / "out")))
    assert "qc.max_frac_untrimmed" in {s.chiave for s in tabella.soglie_di(cercato)}
    assert tabella.non_eseguito(voce, cercato) is None
    assert "filter.primer_removal" in tabella.non_eseguito(
        tabella.soglia("qc.max_primer_hit_frac"), cercato)
    spento = tabella.valori_piatti(valida(_cercato(
        dataset, tmp_path / "out", qc={"max_frac_untrimmed": 1.0})))
    (controllo,) = [s for s in tabella.controlli_spenti(spento, 2)
                    if s.chiave == "qc.max_frac_untrimmed"]
    assert controllo.codici == ("E-S2-04",) and not controllo.ferma


# --------------------------------------------------------------------------- #
# La validazione                                                               #
# --------------------------------------------------------------------------- #


def test_la_validazione_misura_dopo_il_primer_trovato(tmp_path):
    """
    **Obiettivo**: Verificare su un dataset sintetico con il primer preceduto
    da zero a sette basi che la validazione accetti il primer presente (la
    soglia del primer in testa non si applica), riporti in
    ``letture_ispezionate.tsv`` la quota di letture in cui lo trova e la
    lunghezza delle letture dopo il taglio, giudichi il motivo dopo il primer
    trovato, e si fermi con ``E-S0-09`` se il troncamento supera cio' che resta
    dopo il taglio pur stando sotto la lettura intera.

    **Razionale scientifico e sistemistico**: Con il primer cercato
    ``filter.truncLen`` e' la lunghezza dopo il taglio. Un troncamento valido
    per la lettura intera e non per quella tagliata azzererebbe i campioni nel
    filtro, ore dopo.
    """
    dataset = sintetico.genera(
        tmp_path / "dati", biologici=3, negativi=0, letture=200, primer=PRIMER, sfasamento=7)
    # Il motivo con cui comincia ogni sequenza del riferimento sintetico dopo il primer.
    prime = {s[:4] for s in _letture(sorted(dataset.letture.iterdir())[0]).values()}
    assert len(prime) > 4   # le letture intere non cominciano allo stesso modo
    dati = _cercato(dataset, tmp_path / "out", qc={"max_primer_hit_frac": 0.0, "head_reads": 200})
    esito, run = _esegui(dati, Passo.S0)
    assert esito.conclusione is Conclusione.COMPLETATA
    righe = _tabella(tmp_path / "out" / "01_input_validation" / "letture_ispezionate.tsv")
    assert len(righe) == 3
    for riga in righe:
        assert float(riga["frazione_primer"]) > 0.95
        # Regione e coda: la stessa lunghezza qualunque sia lo sfasamento.
        assert riga["lunghezza_minima"] == riga["lunghezza_massima"] == str(
            137 + len(sintetico.CODA))

    # Una finestra troppo stretta: il primer non si trova in gran parte delle
    # letture, e l'arresto sul motivo dice dove guardare.
    stretta = _cercato(dataset, tmp_path / "out3", filter={"primer_max_offset": 1},
                       qc={"conserved_motif": "^TAC", "min_motif_frac": 0.8})
    esito, _ = _esegui(stretta, Passo.S0)
    assert esito.conclusione is Conclusione.ARRESTATA and esito.punto.codice == "E-S0-16"
    assert "filter.primer_max_offset" in esito.punto.dettaglio
    assert "filter.trimLeft" not in esito.punto.dettaglio

    lungo = _cercato(dataset, tmp_path / "out2", filter={"truncLen": 137 + len(sintetico.CODA) + 5,
                                                         "truncLen_shortfall_warn": 0})
    esito, _ = _esegui(lungo, Passo.S0)
    assert esito.conclusione is Conclusione.ARRESTATA and esito.punto.codice == "E-S0-09"
    # Lo stesso troncamento sta sotto la lettura intera piu' corta: a posizione
    # fissa, senza taglio, non sarebbe un difetto di lunghezza.
    intere = min(len(s) for p in dataset.letture.iterdir() for s in _letture(p).values())
    assert lungo["filter"]["truncLen"] <= intere


# --------------------------------------------------------------------------- #
# Il filtro                                                                    #
# --------------------------------------------------------------------------- #


def _attese(percorso: Path, ricerca: RicercaPrimer) -> tuple[dict[str, str], Counter]:
    """Con la regola in Python: per ogni lettura con il primer cio' che resta
    dopo il taglio, e quante letture per (sfasamento, discordanze).
    """
    resto: dict[str, str] = {}
    posizioni: Counter = Counter()
    for nome, sequenza in _letture(percorso).items():
        trovata = ricerca.trova(sequenza)
        if trovata is not None:
            posizioni[trovata] += 1
            resto[nome] = sequenza[trovata[0] + len(ricerca.primer):]
    return resto, posizioni


def test_il_filtro_taglia_ogni_lettura_dove_il_primer_finisce(con_r, tmp_path):
    """
    **Obiettivo**: Verificare, su un dataset sintetico con il primer preceduto
    da zero a sette basi e una quota di letture senza primer, che il filtro
    trovi il primer in ogni lettura agli stessi sfasamenti e con le stesse
    discordanze della regola in Python; che ogni lettura filtrata sia la
    lettura grezza tagliata dopo il primer e troncata; che le letture senza
    primer siano scartate e contate; che il tracciamento riporti il passo del
    taglio fra l'ingresso e il filtro; e che l'avviso ``E-S2-04`` resti nel
    manifesto quando la quota supera la soglia.

    **Razionale scientifico e sistemistico**: La validazione e ``init``
    cercano il primer in Python, il filtro in R sui file interi: se le due
    realizzazioni della regola divergessero, il troncamento verificato non
    sarebbe quello applicato. La prova lettura per lettura e' la sola che
    esclude un taglio alla base sbagliata.
    """
    dataset = sintetico.genera(
        tmp_path / "dati", biologici=3, negativi=1, letture=400, primer=PRIMER,
        sfasamento=7, senza_primer=0.05)
    dati = _cercato(dataset, tmp_path / "out", qc={"max_frac_untrimmed": 0.01})
    # Il filtro non dipende dal profilo: le due fasi si chiedono una per una.
    for passo in (Passo.S2, Passo.S1):
        esito, run = _esegui(dati, passo)
        assert esito.conclusione is Conclusione.COMPLETATA
    filtrati = tmp_path / "out" / "03_filtered"
    ricerca = RicercaPrimer(PRIMER, 1, 7)
    taglio = {r["campione"]: r for r in _tabella(filtrati / "taglio_primer.tsv")}
    posizioni = _tabella(filtrati / "posizioni_primer.tsv")
    senza_totali = 0
    for grezzo in sorted(dataset.letture.iterdir()):
        campione = grezzo.name.split("_")[0]
        resto, attese = _attese(grezzo, ricerca)
        grezze = _letture(grezzo)
        riga = taglio[campione]
        assert int(riga["letture"]) == len(grezze)
        assert int(riga["con_primer"]) == int(riga["rimaste"]) == len(resto)
        assert int(riga["senza_primer"]) == len(grezze) - len(resto) > 0
        senza_totali += len(grezze) - len(resto)
        del_campione = {
            (int(r["sfasamento"]), int(r["discordanze"])): int(r["letture"])
            for r in posizioni if r["campione"] == campione}
        assert del_campione == dict(attese)
        assert {s for s, _ in del_campione} == set(range(8))
        filtrate = _letture(filtrati / f"{campione}_filt.fastq.gz")
        assert filtrate and set(filtrate) <= set(resto)
        for nome, sequenza in filtrate.items():
            assert sequenza == resto[nome][:137]
    assert not (filtrati / ".scrittura-taglio").exists()
    # Il tracciamento: ingresso, dopo il taglio, dopo il filtro.
    tracciamento = componi(run.albero, [run.fase(p) for p in (Passo.S1, Passo.S2)])
    assert tracciamento.passi == ("grezze", "prefiltro", "taglio_primer", "filtrate")
    for campione, letture in tracciamento.letture.items():
        assert letture["grezze"] == letture["prefiltro"] == int(taglio[campione]["letture"])
        assert letture["taglio_primer"] == int(taglio[campione]["rimaste"])
        assert letture["filtrate"] <= letture["taglio_primer"]
    manifesto = json.loads((filtrati / "manifest_S2.json").read_text(encoding="utf-8"))
    assert manifesto["metriche"]["taglio_primer"]["senza_primer"] == senza_totali
    assert [d["codice"] for d in manifesto["degradazioni"]] == ["E-S2-04"]
    assert "qc.max_frac_untrimmed" in manifesto["degradazioni"][0]["dettaglio"]
    # Il profilo di S1 e' quello delle letture tagliate.
    lunghezze = {int(r["lunghezza"]) for r in _tabella(
        tmp_path / "out" / "02_qc_profiles" / "lunghezze.tsv")}
    assert lunghezze == {137 + len(sintetico.CODA)}


def test_il_taglio_cercato_equivale_al_fisso_dove_il_primer_e_in_testa(con_r, tmp_path):
    """
    **Obiettivo**: Verificare, su letture con il primer a posizione fissa, che
    il taglio cercato con il troncamento del taglio fisso meno la lunghezza del
    primer dia le stesse letture filtrate, base per base e qualita' per
    qualita', limitando il confronto alle letture in cui la ricerca trova il
    primer alla prima base; e che le letture filtrate dal solo taglio fisso
    siano quelle in cui la ricerca non lo trova.

    **Razionale scientifico e sistemistico**: Sui protocolli senza
    distanziatori i due modi devono coincidere: una differenza direbbe che il
    taglio cercato sposta il punto in cui le letture cominciano o cambia cio'
    che il filtro vede.
    """
    dataset = sintetico.genera(
        tmp_path / "dati", biologici=3, negativi=0, letture=400, primer=PRIMER,
        senza_primer=0.04)
    fisso = sintetico.configurazione(dataset, tmp_path / "fisso")
    fisso["qc"]["max_primer_hit_frac"] = 1.0
    assert (fisso["filter"]["trimLeft"], fisso["filter"]["truncLen"]) == (19, 156)
    cercato = _cercato(dataset, tmp_path / "cercato", filter={"primer_max_offset": 0})
    assert cercato["filter"]["truncLen"] == fisso["filter"]["truncLen"] - len(PRIMER)
    for dati in (fisso, cercato):
        esito, _ = _esegui(dati, Passo.S2)
        assert esito.conclusione is Conclusione.COMPLETATA

    def record(percorso: Path) -> dict[str, tuple[str, str]]:
        with gzip.open(percorso, "rt") as file:
            righe = file.read().splitlines()
        return {righe[i][1:].split()[0]: (righe[i + 1], righe[i + 3])
                for i in range(0, len(righe), 4)}

    ricerca = RicercaPrimer(PRIMER, 1, 0)
    for grezzo in sorted(dataset.letture.iterdir()):
        campione = grezzo.name.split("_")[0]
        con_primer = {n for n, s in _letture(grezzo).items() if ricerca.trova(s) is not None}
        del_fisso = record(tmp_path / "fisso" / "03_filtered" / f"{campione}_filt.fastq.gz")
        del_cercato = record(tmp_path / "cercato" / "03_filtered" / f"{campione}_filt.fastq.gz")
        assert del_cercato and set(del_cercato) <= con_primer
        assert del_cercato == {n: r for n, r in del_fisso.items() if n in con_primer}
        assert set(del_fisso) - set(del_cercato) == set(del_fisso) - con_primer
        assert set(del_fisso) - con_primer   # il taglio fisso tiene letture senza primer


def test_le_letture_senza_primer_si_possono_tenere_intere(con_r, tmp_path):
    """
    **Obiettivo**: Verificare che con ``filter.primer_untrimmed: keep`` le
    letture senza primer arrivino al filtro intere (le loro prime basi sono
    quelle della lettura grezza), quelle con il primer tagliate, e che la
    tabella del taglio le conti come rimaste.

    **Razionale scientifico e sistemistico**: Chi sceglie di tenere le letture
    senza primer deve ritrovarle nel risultato: scartarle comunque, o tagliarle
    a una posizione qualunque, tradirebbe il parametro dichiarato.
    """
    dataset = sintetico.genera(
        tmp_path / "dati", biologici=2, negativi=0, letture=300, primer=PRIMER,
        sfasamento=5, senza_primer=0.1)
    dati = _cercato(dataset, tmp_path / "out", filter={"primer_untrimmed": "keep"},
                    qc={"max_frac_untrimmed": 1.0})
    esito, _ = _esegui(dati, Passo.S2)
    assert esito.conclusione is Conclusione.COMPLETATA
    ricerca = RicercaPrimer(PRIMER, 1, 7)
    filtrati = tmp_path / "out" / "03_filtered"
    taglio = {r["campione"]: r for r in _tabella(filtrati / "taglio_primer.tsv")}
    for grezzo in sorted(dataset.letture.iterdir()):
        campione = grezzo.name.split("_")[0]
        grezze = _letture(grezzo)
        resto, _ = _attese(grezzo, ricerca)
        assert int(taglio[campione]["rimaste"]) == int(taglio[campione]["letture"]) == len(grezze)
        assert int(taglio[campione]["senza_primer"]) == len(grezze) - len(resto) > 0
        filtrate = _letture(filtrati / f"{campione}_filt.fastq.gz")
        intere = [n for n in filtrate if n not in resto]
        assert intere
        for nome, sequenza in filtrate.items():
            assert sequenza == (resto[nome] if nome in resto else grezze[nome])[:137]


def test_una_coppia_resta_solo_se_restano_le_sue_due_letture(con_r, tmp_path):
    """
    **Obiettivo**: Verificare, in modalita' appaiata con i due primer preceduti
    da basi di sfasamento, che il primer forward si cerchi nelle forward e il
    reverse nelle reverse; che una coppia in cui manca uno dei due primer sia
    scartata per intero; che i due file filtrati abbiano le stesse coppie nello
    stesso ordine; e che ogni lettura filtrata sia la sua lettura grezza
    tagliata dopo il proprio primer.

    **Razionale scientifico e sistemistico**: La fusione prende le coppie per
    posizione: una lettura rimasta senza compagna sfaserebbe tutte le
    successive, e ogni coppia fusa unirebbe due frammenti diversi.
    """
    dataset = sintetico.genera_coppie(
        tmp_path / "dati", biologici=2, coppie=300, lunghezza=130, primer=PRIMER,
        primer_reverse=sintetico.PRIMER_805R, sfasamento=7, senza_primer=0.05)
    dati = sintetico.configurazione_appaiata(dataset, tmp_path / "out")
    dati["filter"].update(CERCA, truncLen=100)
    dati["pair"].update(trimLeft_reverse=0, truncLen_reverse=100)
    dati["qc"]["max_frac_untrimmed"] = 1.0
    esito, _ = _esegui(dati, Passo.S2)
    assert esito.conclusione is Conclusione.COMPLETATA
    filtrati = tmp_path / "out" / "03_filtered"
    taglio = {r["campione"]: r for r in _tabella(filtrati / "taglio_primer.tsv")}
    assert {r["lettura"] for r in _tabella(filtrati / "posizioni_primer.tsv")} == {
        "forward", "reverse"}
    del_forward, del_reverse = RicercaPrimer(PRIMER, 1, 7), RicercaPrimer(sintetico.PRIMER_805R, 1, 7)
    for campione, (forward, reverse) in dataset.coppie.items():
        resto_f, _ = _attese(forward, del_forward)
        resto_r, _ = _attese(reverse, del_reverse)
        riga = taglio[campione]
        entrambe = set(resto_f) & set(resto_r)
        assert int(riga["con_primer"]) == len(resto_f)
        assert int(riga["con_primer_reverse"]) == len(resto_r)
        assert int(riga["rimaste"]) == len(entrambe) < int(riga["letture"]) == 300
        assert int(riga["senza_primer"]) == 300 - len(entrambe)
        assert set(resto_f) != set(resto_r)   # in qualche coppia manca un solo primer
        ordine_f = _ordine(filtrati / f"{campione}_filt.fastq.gz")
        ordine_r = _ordine(filtrati / f"{campione}_filt_reverse.fastq.gz")
        assert ordine_f == ordine_r and set(ordine_f) <= entrambe
        filtrate_f = _letture(filtrati / f"{campione}_filt.fastq.gz")
        filtrate_r = _letture(filtrati / f"{campione}_filt_reverse.fastq.gz")
        for nome in ordine_f:
            assert filtrate_f[nome] == resto_f[nome][:100]
            assert filtrate_r[nome] == resto_r[nome][:100]


def test_la_perdita_del_filtro_dice_quanto_viene_dal_taglio(con_r, tmp_path):
    """
    **Obiettivo**: Verificare che, con una finestra troppo stretta per i
    distanziatori del dataset, l'arresto del filtro per la perdita di letture
    (``E-S2-02``) dica quante letture sono senza primer e nomini
    ``filter.primer_max_offset`` e la tabella del taglio.

    **Razionale scientifico e sistemistico**: Con il primer cercato la perdita
    del filtro comprende le letture senza primer: un messaggio che mandasse a
    rivedere gli errori attesi e il troncamento porterebbe lontano dalla causa.
    """
    dataset = sintetico.genera(
        tmp_path / "dati", biologici=2, negativi=0, letture=300, primer=PRIMER, sfasamento=7)
    dati = _cercato(dataset, tmp_path / "out", filter={"primer_max_offset": 1},
                    qc={"max_frac_lost_filter": 0.5, "max_frac_untrimmed": 1.0})
    esito, _ = _esegui(dati, Passo.S2)
    assert esito.conclusione is Conclusione.ARRESTATA and esito.punto.codice == "E-S2-02"
    dettaglio = esito.punto.dettaglio
    assert "filter.primer_max_offset" in dettaglio and "taglio_primer.tsv" in dettaglio
    senza = sum(int(r["senza_primer"]) for r in _tabella(
        tmp_path / "out" / "03_filtered" / "taglio_primer.tsv"))
    assert f"{senza} su 600" in dettaglio and senza > 300


def test_nessuna_lettura_con_il_primer_ferma_il_profilo(con_r, tmp_path):
    """
    **Obiettivo**: Verificare che con il primer cercato un dataset le cui
    letture non portano il primer si fermi nel profilo con ``E-S1-12``, che
    nomina il primer e la finestra, senza arrivare al filtro.

    **Razionale scientifico e sistemistico**: Con il primer cercato ogni
    lettura senza primer e' scartata: se nessuna lo porta il filtro azzererebbe
    ogni campione, e l'arresto direbbe «campioni senza letture» invece della
    causa.
    """
    dataset = sintetico.genera(tmp_path / "dati", biologici=2, negativi=0, letture=120)
    dati = _cercato(dataset, tmp_path / "out", qc={"primer_sequence": PRIMER,
                                                   "max_frac_short_reads": 1.0})
    esito, _ = _esegui(dati, Passo.S1)
    assert esito.conclusione is Conclusione.ARRESTATA
    assert (str(esito.punto.passo), esito.punto.codice) == ("S1", "E-S1-12")
    assert PRIMER in esito.punto.dettaglio and "26 basi" in esito.punto.dettaglio
    # L'esecuzione intera si ferma li': il profilo viene prima del filtro.
    config = tmp_path / "config.yaml"
    config.write_text(yaml.safe_dump(dati), encoding="utf-8")
    assert main(["resume", "--config", str(config)]) == 4
    assert not (tmp_path / "out" / "03_filtered" / "manifest_S2.json").exists()
    assert CATALOGO["E-S1-12"].categoria is Categoria.REVISIONE_UMANA


# --------------------------------------------------------------------------- #
# init                                                                         #
# --------------------------------------------------------------------------- #


def _init(dataset, uscita: Path, *altri: str) -> dict[str, tuple[bool, str]]:
    assert main(["init", "--reads", str(dataset.letture), "--metadata", str(dataset.tabella),
                 "--output", str(uscita), "--skip-profile", *altri]) == 0
    return righe_dei_parametri(uscita.read_text("utf-8"))


def test_init_propone_la_ricerca_solo_se_il_primer_non_ha_una_posizione_fissa(tmp_path, capsys):
    """
    **Obiettivo**: Verificare che ``init`` proponga ``filter.primer_removal:
    search`` (con il primer, il taglio fisso a zero, le discordanze di partenza
    e lo sfasamento massimo osservato) su letture con il primer preceduto da
    zero a sette basi; che su letture con il primer in testa proponga il taglio
    fisso e non scriva alcuna riga del taglio cercato, anche se in poche
    letture una base in piu' sposta il primer; e che il motivo proposto con la
    ricerca sia quello delle basi che seguono il primer.

    **Razionale scientifico e sistemistico**: La ricerca serve ai protocolli
    con i distanziatori. Proporla per qualche errore di sequenza cambierebbe il
    metodo di un dataset a posizione fissa senza motivo; non proporla dove
    serve lascerebbe le letture sfasate, e lo stesso taxon in piu' varianti.
    """
    sfasato = sintetico.genera(
        tmp_path / "sfasato", biologici=4, negativi=0, letture=300, primer=PRIMER, sfasamento=7)
    righe = _init(sfasato, tmp_path / "sfasato.yaml")
    assert righe["filter.primer_removal"] == (False, '"search"')
    assert righe["filter.trimLeft"] == (False, "0")
    assert righe["filter.primer_max_mismatch"] == (False, "1")
    assert righe["filter.primer_max_offset"] == (False, "7")
    assert righe["filter.primer_untrimmed"] == (False, '"discard"')
    assert righe["qc.max_frac_untrimmed"] == (False, "")
    assert righe["qc.primer_sequence"] == (False, f'"{PRIMER}"')
    testo = (tmp_path / "sfasato.yaml").read_text("utf-8")
    assert "cutadapt" in testo and "sfasamento" in testo
    # Il motivo e' quello delle basi dopo il primer: lo portano le letture tagliate.
    motivo = yaml.safe_load(righe["qc.conserved_motif"][1])
    ricerca = RicercaPrimer(PRIMER, 1, 7)
    import re
    letture = [s for p in sfasato.letture.iterdir() for s in _letture(p).values()]
    con_motivo = sum(bool(re.match(motivo, s[ricerca.fine(s) or 0:])) for s in letture)
    assert con_motivo > 0.9 * len(letture)
    # Con l'opzione: proposte attive, soglia al valore che non ferma, e il file
    # completato delle sole scelte supera la validazione dello schema.
    attive = _init(sfasato, tmp_path / "accettate.yaml", "--accetta-proposte")
    assert attive["filter.primer_removal"] == (True, '"search"')
    assert attive["filter.primer_max_offset"] == (True, "7")
    assert attive["qc.max_frac_untrimmed"] == (True, "1.0")

    fisso = sintetico.genera(
        tmp_path / "fisso", biologici=4, negativi=0, letture=300, primer=PRIMER)
    # In tre letture su cento una base in piu' sposta il primer: errori di
    # sequenza isolati, sotto la quota che distingue un protocollo sfasato.
    for percorso in fisso.letture.iterdir():
        with gzip.open(percorso, "rt") as file:
            record = file.read().splitlines()
        for i in range(0, len(record), 4 * 33):
            record[i + 1], record[i + 3] = "A" + record[i + 1], "I" + record[i + 3]
        with gzip.open(percorso, "wt") as file:
            file.write("\n".join(record) + "\n")
    righe = _init(fisso, tmp_path / "fisso.yaml")
    assert righe["filter.trimLeft"][1] == str(len(PRIMER))
    assert not {"filter.primer_removal", *defaults.PARAMETRI_DEL_TAGLIO_CERCATO} & set(righe)
    assert "primer_removal" not in (tmp_path / "fisso.yaml").read_text("utf-8")
    spostate = sum(
        1 for p in fisso.letture.iterdir() for s in _letture(p).values()
        if RicercaPrimer(PRIMER, 1, 7).trova(s) == (1, 0))
    assert 0 < spostate < rilevamento.QUOTA_FUORI_POSIZIONE * 4 * 300
    capsys.readouterr()


def test_init_misura_l_inserto_senza_i_primer(tmp_path, capsys):
    """
    **Obiettivo**: Verificare che in modalita' appaiata, con i due primer
    preceduti da basi di sfasamento, ``init`` proponga la ricerca per entrambe
    le letture (taglio reverse a zero, primer reverse) e misuri la lunghezza
    dell'inserto senza i primer e senza le basi che li precedono: la lunghezza
    attesa e l'intervallo delle sequenze fuse sono quelli della regione.

    **Razionale scientifico e sistemistico**: Con il primer cercato i
    troncamenti sono lunghezze dopo il taglio, e il vincolo della fusione si
    scrive sull'inserto senza i primer: una lunghezza che li comprendesse
    porterebbe a troncamenti piu' lunghi del necessario e a un intervallo che
    nessuna sequenza fusa rispetta.
    """
    dataset = sintetico.genera_coppie(
        tmp_path / "dati", biologici=3, coppie=300, lunghezza=130, primer=PRIMER,
        primer_reverse=sintetico.PRIMER_805R, sfasamento=7)
    righe = _init(dataset, tmp_path / "c.yaml")
    assert righe["io.layout"] == (True, '"paired-end"')
    assert righe["filter.primer_removal"] == (False, '"search"')
    assert righe["filter.trimLeft"] == (False, "0")
    assert righe["pair.trimLeft_reverse"] == (False, "0")
    assert righe["pair.primer_sequence_reverse"] == (False, f'"{sintetico.PRIMER_805R}"')
    assert righe["pair.expected_merged_len"] == (False, str(sintetico.INSERTO))
    assert righe["pair.len_min"][1] == righe["pair.len_max"][1] == str(sintetico.INSERTO)
    capsys.readouterr()


# --------------------------------------------------------------------------- #
# La catena sulla fixture derivata                                             #
# --------------------------------------------------------------------------- #

#: La fixture appaiata del repository: le sue letture non portano primer.
FIXTURE_APPAIATA = Path(__file__).resolve().parent / "fixtures" / "osd276" / "appaiato"
#: Il classificatore addestrato sul riferimento sintetico di quella fixture.
ADDESTRATO_APPAIATO = (
    Path(__file__).resolve().parent / "fixtures" / "idtaxa" / "addestrato_appaiato.RData")


def _fixture_derivata(cartella: Path) -> dict:
    """La fixture derivata e la sua configurazione.

    Le letture sono quelle della fixture appaiata con, in testa a ciascuna, il
    suo primer preceduto da zero a sette basi casuali: aggiunte in modo
    sintetico e dichiarato, per provare il taglio cercato su letture vere. Dopo
    il taglio le letture sono di nuovo quelle della fixture. La tassonomia e'
    quella del riferimento sintetico (nomi inventati): la catena non dice nulla
    della composizione reale dei campioni.
    """
    dati = yaml.safe_load((FIXTURE_APPAIATA / "config.yaml").read_text(encoding="utf-8"))
    radice = Path(__file__).resolve().parents[1]
    for chiave in ("assay_table", "study_table"):
        dati["io"][chiave] = str(radice / dati["io"][chiave])
    letture = cartella / "letture"
    letture.mkdir(parents=True)
    for numero, origine in enumerate(sorted((FIXTURE_APPAIATA / "fastq").iterdir())):
        primer = (dati["qc"]["primer_sequence"] if "_R1" in origine.name
                  else dati["pair"]["primer_sequence_reverse"])
        sintetico.anteponi_primer(origine, letture / origine.name, primer, 7, seme=numero)
    dati["io"]["fastq_dir"] = str(letture)
    dati["filter"].update(CERCA)
    dati["qc"]["max_frac_untrimmed"] = 0.1
    return dati


def _con_idtaxa(dati: dict) -> dict:
    """La stessa configurazione con IdTaxa sul classificatore sintetico."""
    import hashlib

    dati = copy.deepcopy(dati)
    for chiave in ("ref_fasta", "ref_md5", "try_rc"):
        del dati["tax"][chiave]
    dati["tax"].update(
        classifier="idtaxa", idtaxa_trained=str(ADDESTRATO_APPAIATO),
        idtaxa_trained_md5=hashlib.md5(ADDESTRATO_APPAIATO.read_bytes()).hexdigest(),
        idtaxa_threshold=60, idtaxa_strand="both")
    return dati


def _catena(dati: dict, cartella: Path, nome: str) -> Path:
    dati = copy.deepcopy(dati)
    dati["io"]["out_root"] = str(cartella / nome)
    config = cartella / f"{nome}.yaml"
    config.write_text(yaml.safe_dump(dati, sort_keys=False), encoding="utf-8")
    assert main(["run", "--config", str(config)]) == 0, nome
    chiudi()
    return cartella / nome


def _impronte(uscita: Path, *cartelle: str) -> dict[str, str]:
    """Lo sha256 dichiarato degli artefatti delle cartelle di fase indicate."""
    return {
        f"{m.parent.name}/{a['nome']}": a["checksum"]
        for m in uscita.rglob("manifest_S*.json") if m.parent.name in cartelle
        for a in json.loads(m.read_text(encoding="utf-8"))["artefatti"]}


@pytest.mark.xdist_group("catena_appaiata")
def test_la_catena_appaiata_con_il_primer_cercato_e_idtaxa_arriva_in_fondo(con_r, tmp_path):
    """
    **Obiettivo**: Verificare che la catena completa in modalita' appaiata, sulla
    fixture derivata (le letture della fixture appaiata con i due primer e le
    basi di sfasamento anteposti), con il primer cercato lungo la lettura e il
    classificatore IdTaxa, arrivi all'oggetto finale: ogni coppia porta i due
    primer e resta; le letture filtrate sono quelle della fixture di partenza,
    base per base; le coppie si fondono nella lunghezza attesa della regione; e
    la tassonomia ha le sei colonne con la confidenza.

    **Razionale scientifico e sistemistico**: Le due funzioni si usano insieme
    sui protocolli con i distanziatori. Una lettura tagliata una base fuori
    posto si vedrebbe qui: le letture filtrate non coinciderebbero con quelle
    di partenza, e le sequenze fuse uscirebbero dall'intervallo delle
    lunghezze.
    """
    dati = _con_idtaxa(_fixture_derivata(tmp_path))
    uscita = _catena(dati, tmp_path, "out")
    taglio = _tabella(uscita / "03_filtered" / "taglio_primer.tsv")
    assert len(taglio) == 4
    for riga in taglio:
        assert riga["letture"] == riga["con_primer"] == riga["con_primer_reverse"] == riga["rimaste"]
        assert riga["senza_primer"] == "0"
    sfasamenti = {(r["lettura"], int(r["sfasamento"]))
                  for r in _tabella(uscita / "03_filtered" / "posizioni_primer.tsv")}
    assert sfasamenti == {(lettura, s) for lettura in ("forward", "reverse") for s in range(8)}
    # Le letture filtrate sono quelle della fixture senza primer.
    for filtrato in sorted((uscita / "03_filtered").glob("*_filt*.fastq.gz")):
        campione = filtrato.name.split("_filt")[0]
        verso = "_R2" if "reverse" in filtrato.name else "_R1"
        (origine,) = [p for p in (FIXTURE_APPAIATA / "fastq").iterdir()
                      if p.name.startswith(campione) and verso in p.name]
        di_partenza = _letture(origine)
        filtrate = _letture(filtrato)
        assert len(filtrate) > 1000
        assert all(sequenza == di_partenza[nome] for nome, sequenza in filtrate.items())
    fusione = json.loads((uscita / "05_pair_merge" / "fusione.json").read_text(encoding="utf-8"))
    assert fusione["quota_controllata"] > 0.5 and fusione["lunghezza_mediana"] == 253
    righe = _tabella(uscita / "12_final" / "tassonomia.tsv")
    assert list(righe[0]) == ["asv_id", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus"]
    assert len(righe) > 20 and any(r["Genus"] for r in righe)
    assert (uscita / "08_taxonomy" / "confidenza.tsv").is_file()
    assert not (uscita / "08_taxonomy" / "bootstrap.tsv").exists()
    manifesto = json.loads(
        (uscita / "03_filtered" / "manifest_S2.json").read_text(encoding="utf-8"))
    assert not manifesto.get("degradazioni")


@pytest.mark.lenta
@pytest.mark.xdist_group("catena_appaiata")
def test_la_fixture_derivata_da_gli_stessi_artefatti_della_fixture_senza_primer(con_r, tmp_path):
    """
    **Obiettivo**: Verificare che la catena appaiata sulla fixture derivata,
    con il primer cercato e il classificatore della fixture, dia dall'inferenza
    in poi gli stessi artefatti, byte per byte, della catena sulla fixture di
    partenza con il taglio a posizione fissa; che due esecuzioni sulla fixture
    derivata diano gli stessi byte in ogni fase; e che in modalita' single-end,
    sulle sole letture forward della fixture derivata, la catena arrivi in
    fondo con varianti della lunghezza del troncamento.

    **Razionale scientifico e sistemistico**: Dopo il taglio del primer le
    letture derivate sono quelle di partenza: ogni differenza a valle direbbe
    che il taglio cercato cambia qualcosa oltre le basi che toglie. E il taglio
    non deve introdurre nulla di non ripetibile.
    """
    derivata = _fixture_derivata(tmp_path)
    radice = Path(__file__).resolve().parents[1]
    derivata["tax"]["ref_fasta"] = str(radice / derivata["tax"]["ref_fasta"])
    di_partenza = copy.deepcopy(derivata)
    di_partenza["io"]["fastq_dir"] = str(FIXTURE_APPAIATA / "fastq")
    for chiave in CERCA:
        del di_partenza["filter"][chiave]
    di_partenza["filter"]["trimLeft"] = 0
    del di_partenza["qc"]["max_frac_untrimmed"]
    a_valle = ("04_error_models", "05_asv_inference", "05_pair_merge", "06_seqtab",
               "07_chimera", "08_taxonomy", "10_phyloseq", "11_controls", "12_final")
    cercata = _catena(derivata, tmp_path, "cercata")
    fissa = _catena(di_partenza, tmp_path, "fissa")
    assert len(_impronte(fissa, *a_valle)) > 50
    assert _impronte(cercata, *a_valle) == _impronte(fissa, *a_valle)
    ripetuta = _catena(derivata, tmp_path, "ripetuta")
    tutte = (*a_valle, "01_input_validation", "02_qc_profiles", "03_filtered")
    assert _impronte(cercata, *tutte) == _impronte(ripetuta, *tutte)

    sole_forward = _fixture_derivata(tmp_path / "forward")
    for file in Path(sole_forward["io"]["fastq_dir"]).glob("*_R2*"):
        file.unlink()
    sole_forward["tax"]["ref_fasta"] = str(radice / sole_forward["tax"]["ref_fasta"])
    del sole_forward["io"]["layout"], sole_forward["pair"]
    sole_forward["filter"]["truncLen"] = 140
    singola = _catena(sole_forward, tmp_path, "singola")
    with open(singola / "12_final" / "sequenze.fasta", encoding="utf-8") as file:
        lunghezze = {len(riga.strip()) for riga in file if not riga.startswith(">")}
    assert lunghezze == {140}


# --------------------------------------------------------------------------- #
# Il taglio a posizione fissa su letture sfasate                               #
# --------------------------------------------------------------------------- #


def _avvisi_del_taglio_fisso(dati: dict, uscita: Path) -> tuple[list[dict], dict]:
    """Esegue la validazione: gli avvisi ``E-S0-29`` di ``gates.json`` e il
    manifesto della fase.
    """
    esito, _ = _esegui(dati, Passo.S0)
    assert esito.conclusione is Conclusione.COMPLETATA
    cartella = uscita / "01_input_validation"
    gates = json.loads((cartella / "gates.json").read_text(encoding="utf-8"))
    avvisi = [a for g in gates["gate"] for a in g["avvisi"] if a["codice"] == "E-S0-29"]
    return avvisi, json.loads((cartella / "manifest_S0.json").read_text(encoding="utf-8"))


def test_il_taglio_fisso_su_letture_sfasate_da_un_avviso(tmp_path):
    """
    **Obiettivo**: Verificare che, con il taglio a posizione fissa dichiarato su
    letture con il primer preceduto da zero a sette basi, la validazione
    prosegua e dia l'avviso ``E-S0-29``, in ``gates.json`` e nel manifesto, con
    la quota delle letture fuori posizione, le posizioni trovate e
    l'indicazione di passare a ``filter.primer_removal: search``; che su
    letture con il primer a posizione fissa non dia alcun avviso, e
    ``gates.json`` non porti nulla in piu'; e che in modalita' appaiata l'avviso
    riguardi anche la lettura reverse.

    **Razionale scientifico e sistemistico**: Un taglio fisso su letture
    sfasate lascia residui del primer in testa a gran parte delle letture: la
    catena arriverebbe in fondo, con lo stesso taxon diviso in piu' varianti,
    senza che nulla lo dica.
    """
    sfasato = sintetico.genera(
        tmp_path / "sfasato", biologici=3, negativi=0, letture=200, primer=PRIMER, sfasamento=7)
    dati = sintetico.configurazione(sfasato, tmp_path / "out")
    dati["qc"].update(max_primer_hit_frac=1.0, head_reads=200, max_frac_short_reads=1.0)
    assert dati["filter"]["trimLeft"] == len(PRIMER)
    (avviso,), manifesto = _avvisi_del_taglio_fisso(dati, tmp_path / "out")
    assert "filter.trimLeft 19" in avviso["dettaglio"]
    assert 'filter.primer_removal: "search"' in avviso["dettaglio"]
    # Un ottavo delle letture ha il primer alla prima base, le altre no.
    assert "0: 1" in avviso["dettaglio"] and "7: 1" in avviso["dettaglio"]
    assert [d["codice"] for d in manifesto["degradazioni"] if d["codice"] == "E-S0-29"] == [
        "E-S0-29"]
    assert CATALOGO["E-S0-29"].categoria is Categoria.DEGRADAZIONE_AUTOMATICA

    fisso = sintetico.genera(
        tmp_path / "fisso", biologici=3, negativi=0, letture=200, primer=PRIMER)
    dati = sintetico.configurazione(fisso, tmp_path / "out2")
    dati["qc"].update(max_primer_hit_frac=1.0, head_reads=200)
    avvisi, manifesto = _avvisi_del_taglio_fisso(dati, tmp_path / "out2")
    assert avvisi == [] and "E-S0-29" not in json.dumps(manifesto)
    righe = _tabella(tmp_path / "out2" / "01_input_validation" / "letture_ispezionate.tsv")
    assert list(righe[0]) == [
        "accession", "file", "letture_esaminate", "file_esaurito", "lunghezza_minima",
        "lunghezza_massima", "frazione_primer", "frazione_motivo", "errore"]

    coppie = sintetico.genera_coppie(
        tmp_path / "coppie", biologici=2, coppie=200, lunghezza=130, primer=PRIMER,
        primer_reverse=sintetico.PRIMER_805R, sfasamento=7)
    dati = sintetico.configurazione_appaiata(coppie, tmp_path / "out3")
    dati["qc"].update(max_frac_short_reads=1.0)
    dati["pair"]["conserved_motif_reverse"] = None
    avvisi, _ = _avvisi_del_taglio_fisso(dati, tmp_path / "out3")
    assert len(avvisi) == 2
    assert "pair.trimLeft_reverse 21" in avvisi[1]["dettaglio"]
    assert sintetico.PRIMER_805R in avvisi[1]["dettaglio"]


def test_il_controllo_del_taglio_fisso_si_applica_solo_se_il_taglio_toglie_il_primer(tmp_path):
    """
    **Obiettivo**: Verificare che il controllo sulla posizione del primer non
    si applichi, e non dia alcun avviso, con il primer nullo, con
    ``filter.trimLeft`` zero e con un taglio piu' corto del primer; che con
    meno di trenta letture con il primer non giudichi la quota; e che una quota
    di letture fuori posizione sotto la soglia non dia l'avviso.

    **Razionale scientifico e sistemistico**: Il controllo giudica un taglio
    pensato per togliere il primer: dove il taglio toglie altro (un codice a
    barre, poche basi di bassa qualita') la posizione del primer non dice
    nulla, e un avviso sarebbe falso. Su poche letture una quota non ha
    significato.
    """
    from amplicon16s.io_layer.reads import SEGNALE_MINIMO_DI_UNA_PARTE, sonda_del_taglio_fisso

    assert sonda_del_taglio_fisso(None, 19) is None
    assert sonda_del_taglio_fisso(PRIMER, 0) is None
    assert sonda_del_taglio_fisso(PRIMER, len(PRIMER) - 1) is None
    sonda = sonda_del_taglio_fisso(PRIMER, len(PRIMER) + 4)
    assert (sonda.discordanze_massime, sonda.sfasamento_massimo) == (1, 4 + 16)

    sfasato = sintetico.genera(
        tmp_path / "sfasato", biologici=3, negativi=0, letture=200, primer=PRIMER, sfasamento=7)
    for numero, modifiche in enumerate((
        {"filter": {"trimLeft": 10}},
        {"filter": {"trimLeft": 0}, "qc": {"max_primer_hit_frac": 1.0}},
        {"qc": {"primer_sequence": None}},
        # Nove letture per file: ventisette in tutto, meno del minimo.
        {"qc": {"head_reads": 9}},
    )):
        dati = sintetico.configurazione(sfasato, tmp_path / f"out{numero}")
        dati["qc"].update(max_primer_hit_frac=1.0, head_reads=200, max_frac_short_reads=1.0)
        for gruppo, valori in modifiche.items():
            dati[gruppo].update(valori)
        avvisi, _ = _avvisi_del_taglio_fisso(dati, tmp_path / f"out{numero}")
        assert avvisi == [], modifiche
    assert 3 * 9 < SEGNALE_MINIMO_DI_UNA_PARTE

    # Tre letture su cento con una base in piu': sotto la quota, nessun avviso.
    quasi_fisso = sintetico.genera(
        tmp_path / "quasi", biologici=3, negativi=0, letture=200, primer=PRIMER)
    for percorso in quasi_fisso.letture.iterdir():
        with gzip.open(percorso, "rt") as file:
            record = file.read().splitlines()
        for i in range(0, len(record), 4 * 33):
            record[i + 1], record[i + 3] = "A" + record[i + 1], "I" + record[i + 3]
        with gzip.open(percorso, "wt") as file:
            file.write("\n".join(record) + "\n")
    dati = sintetico.configurazione(quasi_fisso, tmp_path / "quasi_out")
    dati["qc"].update(max_primer_hit_frac=1.0, head_reads=200)
    avvisi, _ = _avvisi_del_taglio_fisso(dati, tmp_path / "quasi_out")
    assert avvisi == []
