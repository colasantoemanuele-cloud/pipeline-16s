"""I file dei dataset sotto ``dati/``: configurazioni, elenchi delle letture,
tabelle di abbinamento e degli attesi, fonti.

Verifica cio' che si puo' verificare senza scaricare nulla: che ogni
configurazione della pipeline superi il controllo dei parametri e indichi
tabelle che esistono; che le configurazioni che dichiarano di differire da
un'altra solo per la modalita' o per il classificatore differiscano davvero
solo in quello; che gli elenchi delle letture abbiano la forma che lo scarico
si attende; che le tabelle di abbinamento e dei campioni nominino i campioni
dei metadati; che le tabelle degli attesi di una stessa comunita' riportino la
stessa composizione. Chi esegue una di queste configurazioni si fida del suo
nome: un parametro in piu' cambiato per errore non si vedrebbe nei risultati.
"""

from __future__ import annotations

import csv
import re
from pathlib import Path

import pytest
import yaml

from tests.conftest import esegui_g15

RADICE = Path(__file__).resolve().parents[1]
DATI = RADICE / "dati"

pytestmark = pytest.mark.skipif(not (DATI / "comunita_nota").is_dir(),
                                reason="la cartella dati/ non e' presente")

#: Le configurazioni della pipeline (non quelle delle analisi ecologiche).
CONFIGURAZIONI = sorted(DATI.glob("*/config_*.yaml"))

#: Cio' che puo' cambiare fra due configurazioni dello stesso dataset senza
#: essere una scelta di analisi: la cartella di uscita e i parametri
#: dell'esecuzione.
DELL_ESECUZIONE = ("io.out_root", "run.")
#: Cio' che cambia passando a un altro classificatore.
DEL_CLASSIFICATORE = (*DELL_ESECUZIONE, "tax.", "filt.exclude_taxa")
#: Cio' che cambia passando dalle sole forward alla modalita' appaiata.
DELLA_MODALITA = (*DELL_ESECUZIONE, "io.layout", "io.fastq_dir", "io.fastq_glob", "pair.")
#: In piu', fra una configurazione appaiata e la sua versione con le sole
#: forward scritta a parte: troncamento e taglio della lettura forward.
DELLA_FORWARD = (*DELLA_MODALITA, "filter.truncLen", "filter.trimLeft", "filter.primer_",
                 "qc.max_frac_short_reads", "qc.max_frac_untrimmed", "qc.max_primer_hit_frac",
                 "asv.len_tol")


def _righe(percorso: Path) -> list[dict[str, str]]:
    with open(percorso, encoding="utf-8", newline="") as file:
        return list(csv.DictReader(file, delimiter="\t", quoting=csv.QUOTE_NONE))


def _piatta(percorso: Path) -> dict[str, object]:
    dati = yaml.safe_load(percorso.read_text(encoding="utf-8"))
    return {f"{gruppo}.{nome}": valore for gruppo, voci in dati.items()
            for nome, valore in (voci or {}).items()}


def _diverse(a: Path, b: Path) -> set[str]:
    """I parametri dichiarati con valori diversi, o da una sola delle due."""
    pa, pb = _piatta(a), _piatta(b)
    assente = object()
    return {chiave for chiave in set(pa) | set(pb) if pa.get(chiave, assente) != pb.get(chiave, assente)}


def _coppie_derivate() -> list[tuple[Path, Path, tuple[str, ...]]]:
    """Per ogni configurazione il cui nome dice modalita' e classificatore, la
    configurazione da cui dichiara di differire e i parametri ammessi."""
    coppie = []
    for percorso in CONFIGURAZIONI:
        trovato = re.fullmatch(r"(config_.+?)_(appaiata|forward)_(naive_bayes|idtaxa)", percorso.stem)
        if not trovato:
            continue
        base, modalita, classificatore = trovato.groups()
        cartella = percorso.parent
        semplice = cartella / f"{base}.yaml"
        appaiata = cartella / f"{base}_appaiata_naive_bayes.yaml"
        if classificatore == "idtaxa":
            partenza = appaiata if modalita == "appaiata" else (
                cartella / f"{base}_forward_naive_bayes.yaml")
            partenza = partenza if partenza.is_file() else semplice
            coppie.append((percorso, partenza, DEL_CLASSIFICATORE))
        elif modalita == "appaiata" and semplice.is_file():
            coppie.append((percorso, semplice, DELLA_MODALITA))
        elif modalita == "forward":
            coppie.append((percorso, appaiata, DELLA_FORWARD))
    return coppie


@pytest.mark.parametrize("percorso", CONFIGURAZIONI, ids=lambda p: p.name)
def test_ogni_configurazione_supera_il_controllo_dei_parametri(percorso):
    """
    **Obiettivo**: Verificare che ogni configurazione della pipeline sotto
    ``dati/`` superi lo schema e i controlli di coerenza fra parametri, che le
    tabelle dei metadati che indica esistano nel repository e che dichiari
    l'immagine con cui si esegue.

    **Razionale scientifico e sistemistico**: Una configurazione pubblicata
    accanto ai dati e' un'istruzione per chi la esegue: se non supera il primo
    controllo dell'esecuzione, o indica una tabella che non c'e', l'errore si
    scopre solo dopo aver scaricato le letture.
    """
    dati = yaml.safe_load(percorso.read_text(encoding="utf-8"))
    risolta = esegui_g15(dati)
    for tabella in (risolta.config.io.assay_table, risolta.config.io.study_table):
        if tabella is not None:
            assert (RADICE / tabella).is_file(), tabella
    assert risolta.config.run.container and "@sha256:" in risolta.config.run.container
    cartella = percorso.parent.relative_to(RADICE).as_posix()
    assert str(risolta.config.io.fastq_dir).startswith(cartella + "/")


def test_ci_sono_configurazioni_che_dichiarano_di_derivare_da_un_altra():
    """
    **Obiettivo**: Verificare che il riconoscimento delle configurazioni
    derivate dal nome ne trovi, e che per ciascuna la configurazione di
    partenza esista.

    **Razionale scientifico e sistemistico**: Il controllo sulle differenze
    ammesse e' fatto per coppie ricavate dai nomi dei file: se una convenzione
    di nome cambiasse, il controllo passerebbe senza confrontare nulla.
    """
    coppie = _coppie_derivate()
    assert coppie
    for derivata, partenza, _ in coppie:
        assert partenza.is_file(), f"{derivata.name}: manca {partenza.name}"
        assert partenza != derivata


@pytest.mark.parametrize("derivata,partenza,ammessi", _coppie_derivate(),
                         ids=lambda v: v.name if isinstance(v, Path) else None)
def test_una_configurazione_derivata_cambia_solo_cio_che_il_nome_dice(derivata, partenza, ammessi):
    """
    **Obiettivo**: Verificare che una configurazione il cui nome dichiara la
    modalita' e il classificatore differisca da quella di partenza dello stesso
    dataset solo nei parametri della modalita' (gruppo ``pair``, cartella e
    selezione dei file; per le sole forward anche troncamento e taglio della
    lettura forward), del classificatore (gruppo ``tax`` e i nomi dei taxa da
    escludere) e dell'esecuzione; che con ``idtaxa`` il classificatore sia
    davvero cambiato; e che la modalita' dichiarata sia quella del nome.

    **Razionale scientifico e sistemistico**: Due esecuzioni si possono
    confrontare per dire che cosa cambia con la fusione delle coppie, o con un
    altro classificatore, solo se tutto il resto e' uguale: una soglia, un
    filtro o un parametro dell'inferenza cambiati insieme renderebbero la
    differenza non attribuibile.
    """
    diverse = _diverse(derivata, partenza)
    fuori = sorted(chiave for chiave in diverse if not chiave.startswith(ammessi))
    assert fuori == [], f"{derivata.name} contro {partenza.name}"
    piatta = _piatta(derivata)
    if derivata.stem.endswith("_idtaxa"):
        assert piatta["tax.classifier"] == "idtaxa" and "tax.classifier" in diverse
    else:
        assert piatta.get("tax.classifier", "naive_bayes") == "naive_bayes"
    appaiata = "_appaiata_" in derivata.stem
    assert (piatta.get("io.layout", "single-end") == "paired-end") is appaiata
    assert any(chiave.startswith("pair.") for chiave in piatta) is appaiata


def test_gli_elenchi_delle_letture_hanno_la_forma_che_lo_scarico_si_attende():
    """
    **Obiettivo**: Verificare, per ogni elenco delle letture nel formato dello
    scarico comune (``letture_osdr.tsv``, ``letture_ena.tsv`` con le colonne
    dei due file), che i nomi dei file siano distinti, che i due indirizzi di
    ogni corsa siano diversi, che le dimensioni siano interi positivi e che un
    MD5, dove c'e', abbia 32 cifre esadecimali.

    **Razionale scientifico e sistemistico**: L'elenco e' l'unica descrizione
    dei file da scaricare: una riga ripetuta o due indirizzi uguali darebbero
    un campione in meno, o le stesse letture come forward e come reverse.
    """
    elenchi = [p for p in sorted(DATI.glob("*/letture_*.tsv")) if "url_1" in p.read_text(encoding="utf-8")[:200]]
    assert elenchi
    for elenco in elenchi:
        righe = _righe(elenco)
        assert righe, elenco
        assert set(righe[0]) == {"file", "url_1", "byte_1", "md5_1", "url_2", "byte_2", "md5_2", "lato"}, elenco
        assert len({r["file"] for r in righe}) == len(righe), elenco
        indirizzi = [r[c] for r in righe for c in ("url_1", "url_2")]
        assert len(set(indirizzi)) == len(indirizzi), elenco
        for riga in righe:
            assert riga["file"].endswith(".fastq.gz")
            assert int(riga["byte_1"]) > 0 and int(riga["byte_2"]) > 0
            assert riga["lato"] in ("", "1", "2")
            for colonna in ("md5_1", "md5_2"):
                assert riga[colonna] == "" or re.fullmatch(r"[0-9a-f]{32}", riga[colonna]), elenco


def test_le_tabelle_di_abbinamento_nominano_i_campioni_dei_metadati():
    """
    **Obiettivo**: Verificare che ogni tabella di abbinamento con la tabella
    elaborata da un'altra pipeline abbia due colonne, nomi e identificativi
    distinti, e che i nomi siano fra i valori di una colonna della tabella di
    assay dello stesso dataset.

    **Razionale scientifico e sistemistico**: Il confronto abbina i campioni
    per nome: un nome che i metadati non hanno toglie un campione dal
    confronto, e il risultato si leggerebbe su meno campioni di quanti si
    crede.
    """
    tabelle = sorted(DATI.glob("*/abbinamento_*.tsv"))
    assert tabelle
    for tabella in tabelle:
        righe = _righe(tabella)
        nome, identificativo = list(righe[0])
        nomi = [r[nome] for r in righe]
        assert len(set(nomi)) == len(nomi) and len({r[identificativo] for r in righe}) == len(righe)
        assay = _righe(next(tabella.parent.glob("metadati/a_*.txt")))
        assert nome in assay[0], tabella
        assert set(nomi) <= {r[nome].strip('"') for r in assay}, tabella


def test_le_tabelle_dei_campioni_e_degli_elenchi_delle_letture_concordano():
    """
    **Obiettivo**: Verificare, per i dataset descritti da una tabella dei
    campioni semplice, che ogni campione indichi un file dell'elenco delle
    letture e che ogni file dell'elenco abbia il suo campione.

    **Razionale scientifico e sistemistico**: La pipeline abbina letture e
    metadati dal nome del file: una corsa nell'elenco senza una riga nella
    tabella fermerebbe la validazione dopo lo scarico.
    """
    tabelle = sorted(DATI.glob("*/campioni.tsv"))
    assert tabelle
    for tabella in tabelle:
        campioni = _righe(tabella)
        elenco = _righe(tabella.with_name("letture_ena.tsv"))
        assert {r["file_letture"] for r in campioni} == {r["file"] for r in elenco}, tabella
        assert len({r["campione"] for r in campioni}) == len(campioni), tabella


def test_le_tabelle_degli_attesi_di_una_comunita_riportano_la_stessa_composizione():
    """
    **Obiettivo**: Verificare che le tabelle degli attesi di una stessa
    cartella, una per riferimento tassonomico, abbiano le stesse specie con le
    stesse percentuali e la stessa fonte, e differiscano solo nei nomi di
    genere e famiglia e negli altri nomi; che nessun altro nome coincida con un
    genere atteso.

    **Razionale scientifico e sistemistico**: La composizione teorica e' una
    sola, quella dichiarata dal produttore; cio' che cambia da un riferimento
    all'altro e' come ciascuno scrive i taxa. Una percentuale diversa fra due
    tabelle farebbe dipendere l'esito della verifica dal classificatore per
    una ragione che non e' il classificatore.
    """
    for cartella in sorted({p.parent for p in DATI.glob("*/attesi*.tsv")}):
        tabelle = sorted(cartella.glob("attesi*.tsv"))
        base = _righe(tabelle[0])
        for tabella in tabelle:
            righe = _righe(tabella)
            assert [(r["specie"], r["percentuale"], r["fonte"]) for r in righe] == [
                (r["specie"], r["percentuale"], r["fonte"]) for r in base], tabella
            generi = [r["genere"] for r in righe]
            assert len(set(generi)) == len(generi)
            altri = [n for r in righe for n in (r.get("altri_nomi") or "").split("|") if n]
            assert not set(altri) & set(generi) and len(set(altri)) == len(altri), tabella


def test_le_fonti_elencano_file_che_esistono():
    """
    **Obiettivo**: Verificare che ogni file dichiarato «nel repository» nella
    tabella delle fonti di un dataset esista nella sua cartella, e che un
    controllo MD5 dichiarato per un file del repository corrisponda al file.

    **Razionale scientifico e sistemistico**: La tabella delle fonti dice da
    dove viene ogni file e come verificarlo: una voce che non corrisponde a un
    file, o un'impronta che non e' la sua, toglie valore a tutte le altre.
    """
    import hashlib

    for fonti in sorted(DATI.glob("*/FONTI.tsv")):
        for riga in _righe(fonti):
            if riga["nel_repository"] != "si":
                continue
            file = fonti.parent / riga["percorso"]
            assert file.is_file(), f"{fonti}: {riga['percorso']}"
            if riga["controllo"].startswith("md5:"):
                assert hashlib.md5(file.read_bytes()).hexdigest() == riga["controllo"][4:], file
