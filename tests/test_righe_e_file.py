"""Righe dei metadati e file di letture: celle con piu' file, righe di un altro
amplicone, righe collegate al file sbagliato.

Verifica cio' che la validazione e ``init`` fanno di una tabella dei campioni
che non corrisponde riga per riga ai file della cartella: una cella che elenca
i due file di una coppia, righe che descrivono campioni di un altro amplicone
dello stesso studio, righe collegate al file di un altro. I casi si costruiscono
su un dataset sintetico (``tests/sintetico.py``), a cui si aggiungono righe e
file di un secondo amplicone con un altro primer in testa alle letture.

Nessun test richiede R: ``init`` si lancia senza il profilo di qualita'.
"""

from __future__ import annotations

import csv
import gzip
import random
import re
from pathlib import Path

import pytest
import sintetico
from init_aiuti import righe_dei_parametri

from amplicon16s import cli
from amplicon16s.config.resolve import risolvi
from amplicon16s.config.schema import valida
from amplicon16s.gates import g01_g15
from amplicon16s.gates.registry import esegui_gate
from amplicon16s.logging.logger import chiudi
from amplicon16s.metadata import ampliconi
from amplicon16s.metadata.crosswalk import (
    AccessionNonEstraibile,
    analizza,
    chiave_dalla_tabella,
)
from amplicon16s.steps.s00_validate import esegui_s0
from conftest import esegui_g15

#: Il primer forward dell'altro amplicone dei casi di prova: l'ITS1F di Gardes
#: e Bruns 1993, che apre la regione ITS dei funghi.
PRIMER_ITS = "CTTGGTCATTTAGAGGAAGTAA"

COLONNA_AMPLICONE = ampliconi.COLONNA_AMPLICONE
COLONNA_PRIMER = ampliconi.COLONNA_PRIMER


@pytest.fixture(autouse=True)
def uscite_pulite():
    """Chiude le uscite del log prima e dopo ogni test."""
    chiudi()
    yield
    chiudi()


def _leggi(tabella: Path) -> list[dict[str, str]]:
    with open(tabella, encoding="utf-8", newline="") as file:
        return list(csv.DictReader(file, delimiter="\t"))


def _scrivi(tabella: Path, righe: list[dict[str, str]]) -> None:
    colonne = list(righe[0])
    tabella.write_text(
        "\t".join(colonne) + "\n"
        + "".join("\t".join(r.get(c, "") for c in colonne) + "\n" for r in righe),
        encoding="utf-8")


def _file_dell_altro_amplicone(cartella: Path, nome: str, letture: int = 60) -> None:
    """Un file di letture che cominciano con il primer dell'altro amplicone."""
    caso = random.Random(nome)
    with gzip.open(cartella / nome, "wt", encoding="ascii") as file:
        for numero in range(letture):
            sequenza = PRIMER_ITS + "".join(caso.choice("ACGT") for _ in range(120))
            file.write(f"@{nome}.{numero + 1}\n{sequenza}\n+\n{'I' * len(sequenza)}\n")


def _due_ampliconi(tmp_path: Path, scambiate: bool = False, con_file_its: bool = False):
    """Un dataset con quattro campioni 16S (primer 341F in testa alle letture)
    e, nella stessa tabella, due righe di un altro amplicone; ogni riga dichiara
    amplicone e primer nelle colonne dei depositi ISA-Tab.

    Senza altro, le righe dell'altro amplicone nominano file che non sono nella
    cartella. Con ``con_file_its`` quei file ci sono. Con ``scambiate`` due
    righe sono collegate ciascuna al file dell'altra: una riga 16S a un file
    dell'altro amplicone e viceversa.
    """
    dataset = sintetico.genera(
        tmp_path / "dati", biologici=4, negativi=0, letture=80, primer=sintetico.PRIMER_341F)
    righe = _leggi(dataset.tabella)
    for riga in righe:
        riga[COLONNA_AMPLICONE] = "16S"
        riga[COLONNA_PRIMER] = f"341F, 5-{sintetico.PRIMER_341F}-3"
    for indice in (5, 6):
        nome_file = f"LIB{indice:03d}_F{indice:02d}.fq.gz"
        if con_file_its or scambiate:
            _file_dell_altro_amplicone(dataset.letture, nome_file)
        righe.append({
            **{c: "" for c in righe[0]},
            sintetico.COLONNA_CAMPIONE: f"F{indice:02d}", sintetico.COLONNA_FILE: nome_file,
            sintetico.COLONNA_TIPO: sintetico.ETICHETTA_BIOLOGICO,
            sintetico.COLONNA_SITO: "sito 1", sintetico.COLONNA_GRUPPO: "trattato",
            COLONNA_AMPLICONE: "ITS", COLONNA_PRIMER: f"ITS1F, 5-{PRIMER_ITS}-3",
        })
    if scambiate:
        righe[0][sintetico.COLONNA_FILE], righe[4][sintetico.COLONNA_FILE] = (
            righe[4][sintetico.COLONNA_FILE], righe[0][sintetico.COLONNA_FILE])
    _scrivi(dataset.tabella, righe)
    return dataset


def _init(dataset, uscita: Path, *altri: str) -> int:
    return cli.main(["init", "--reads", str(dataset.letture), "--metadata", str(dataset.tabella),
                     "--output", str(uscita), "--skip-profile", "--threads", "1", *altri])


# --------------------------------------------------------------------------- #
# Una cella con i due file della coppia                                        #
# --------------------------------------------------------------------------- #


def test_una_cella_con_i_due_file_della_coppia_da_una_chiave_sola(tmp_path):
    """
    **Obiettivo**: Verificare che ``chiave_dalla_tabella`` legga come elenco una
    cella con piu' file separati da virgola: la stessa chiave in ogni elemento
    e' la chiave della riga; chiavi diverse, o un elemento senza chiave,
    restano un errore di ambiguita'; un valore solo con due chiavi resta
    ambiguo; e che l'analisi dei metadati di un dataset appaiato la cui tabella
    riporta i due file di ogni campione non segnali alcuna riga.

    **Razionale scientifico e sistemistico**: Nei depositi la cella del file di
    un campione appaiato riporta i due file della coppia: la stessa chiave
    trovata due volte non e' un'ambiguita', e respingerla obbligherebbe a
    riscrivere i metadati di ogni dataset appaiato.
    """
    espressione = re.compile(r"LIB[0-9]+")
    assert chiave_dalla_tabella("LIB001_a_R1.fq.gz, LIB001_a_R2.fq.gz", espressione) == "LIB001"
    assert chiave_dalla_tabella("LIB001_a_R1.fq.gz", espressione) == "LIB001"
    for valore in ("LIB001_R1.fq.gz, LIB002_R2.fq.gz", "LIB001_LIB002.fq.gz",
                   "LIB001_R1.fq.gz, LIB001_R2.fq.gz, note", "LIB001_LIB001.fq.gz"):
        with pytest.raises(AccessionNonEstraibile, match="ambiguo"):
            chiave_dalla_tabella(valore, espressione)

    dataset = sintetico.genera_coppie(tmp_path / "dati", biologici=3, coppie=40)
    righe = _leggi(dataset.tabella)
    for riga in righe:
        forward, reverse = next(
            coppia for coppia in dataset.coppie.values()
            if coppia[0].name == riga[sintetico.COLONNA_FILE])
        riga[sintetico.COLONNA_FILE] = f"{forward.name}, {reverse.name}"
    _scrivi(dataset.tabella, righe)
    analisi = analizza(valida(sintetico.configurazione_appaiata(dataset, tmp_path / "out")))
    assert not analisi.assay_senza_accession and not analisi.accession_ripetuti_nell_assay
    assert len(analisi.campioni) == 3 and set(analisi.reverse_per_accession) == set(
        dataset.coppie)


def test_init_trova_la_chiave_con_la_cella_dei_due_file(tmp_path):
    """
    **Obiettivo**: Verificare che ``init``, su un dataset appaiato la cui
    tabella riporta in ogni cella i due file della coppia, rilevi la colonna
    della chiave e un'espressione che ``validate`` accetta su quelle celle.

    **Razionale scientifico e sistemistico**: Chi parte dai metadati di un
    deposito non deve scrivere a mano la chiave perche' la cella elenca due
    file: la corrispondenza uno a uno fra righe e campioni c'e', e va trovata.
    """
    dataset = sintetico.genera_coppie(tmp_path / "dati", biologici=3, coppie=40)
    righe = _leggi(dataset.tabella)
    for riga in righe:
        forward, reverse = next(
            coppia for coppia in dataset.coppie.values()
            if coppia[0].name == riga[sintetico.COLONNA_FILE])
        riga[sintetico.COLONNA_FILE] = f"{forward.name}, {reverse.name}"
    _scrivi(dataset.tabella, righe)
    uscita = tmp_path / "c.yaml"
    assert _init(dataset, uscita) == 0
    parametri = righe_dei_parametri(uscita.read_text(encoding="utf-8"))
    assert parametri["meta.accession_column"] == (True, f'"{sintetico.COLONNA_FILE}"')
    attiva, espressione = parametri["io.accession_regex"]
    assert attiva
    compilata = re.compile(espressione.strip('"').encode().decode("unicode_escape"))
    chiavi = {chiave_dalla_tabella(r[sintetico.COLONNA_FILE], compilata) for r in righe}
    assert len(chiavi) == 3


# --------------------------------------------------------------------------- #
# Che cosa una riga dichiara, e a quale file e' collegata                      #
# --------------------------------------------------------------------------- #


def test_la_dichiarazione_di_una_riga_si_interpreta_solo_senza_ambiguita(tmp_path):
    """
    **Obiettivo**: Verificare che l'amplicone dichiarato si riconosca dai nomi
    noti e solo se e' uno; che i primer di una riga siano le sequenze IUPAC
    scritte nella cella e quelle dei primer del vocabolario che la cella
    nomina; che una riga senza l'una ne' gli altri non sia interpretabile; e
    che una cella si colleghi al primo file presente che nomina, elenchi i file
    assenti quando nomina solo file, e non colleghi nulla quando non ne nomina.

    **Razionale scientifico e sistemistico**: Il confronto con le letture puo'
    fermare un'analisi: si fa solo dove la dichiarazione dei metadati non
    lascia dubbi su che cosa dice.
    """
    assert ampliconi.QUOTA_DEL_PRIMER_IN_TESTA == g01_g15.PRIMER_IN_TESTA
    assert ampliconi.amplicone_dichiarato("16S") == "16S"
    assert ampliconi.amplicone_dichiarato("16S rRNA gene") == "16S"
    assert ampliconi.amplicone_dichiarato("ITS2") == "ITS"
    assert ampliconi.amplicone_dichiarato("16S and ITS") is None
    assert ampliconi.amplicone_dichiarato("DNA") is None

    riga = {COLONNA_AMPLICONE: "16S", COLONNA_PRIMER: "515F e 806R"}
    dichiarata = ampliconi.dichiarazione(riga, COLONNA_AMPLICONE, COLONNA_PRIMER)
    assert dichiarata.amplicone == "16S" and "GTGYCAGCMGCCGCGGTAA" in dichiarata.primer
    con_sequenza = ampliconi.dichiarazione(
        {COLONNA_PRIMER: f"forward 5-{PRIMER_ITS}-3"}, None, COLONNA_PRIMER)
    assert con_sequenza.primer == (PRIMER_ITS,) and con_sequenza.amplicone is None
    assert con_sequenza.interpretabile
    assert not ampliconi.dichiarazione({COLONNA_PRIMER: "vedi protocollo"}, None,
                                       COLONNA_PRIMER).interpretabile

    presenti = {"a_R1.fq.gz": tmp_path / "a_R1.fq.gz", "a_R2.fq.gz": tmp_path / "a_R2.fq.gz"}
    assert ampliconi.file_della_cella("a_R1.fq.gz, a_R2.fq.gz", presenti) == ("a_R1.fq.gz", ())
    assert ampliconi.file_della_cella("b_R1.fq.gz, b_R2.fq.gz", presenti) == (
        None, ("b_R1.fq.gz", "b_R2.fq.gz"))
    assert ampliconi.file_della_cella("campione 7", presenti) == (None, ())


def test_una_riga_collegata_al_file_di_un_altro_amplicone_ferma_la_validazione(tmp_path):
    """
    **Obiettivo**: Verificare che, con una riga 16S collegata a un file le cui
    letture cominciano con il primer dell'altro amplicone e una riga dell'altro
    amplicone collegata a un file 16S, la validazione si fermi a G04 con
    ``E-S0-30`` nominando i due campioni, il file collegato e il primer
    trovato; che il gate successivo non si esegua; e che la stessa tabella con
    ogni riga collegata al proprio file superi G04 senza alcun avviso.

    **Razionale scientifico e sistemistico**: Una riga collegata al file di un
    altro campione non si vede dai nomi: i conteggi sarebbero attribuiti al
    campione sbagliato senza alcun errore. Dove i metadati dichiarano che cosa
    il file deve contenere, il contenuto lo puo' smentire prima di ogni calcolo.
    """
    dataset = _due_ampliconi(tmp_path, scambiate=True)
    config = valida(sintetico.configurazione(dataset, tmp_path / "out"))
    esito = esegui_gate("G04", g01_g15.Contesto(config))
    assert not esito.superato and [v.codice for v in esito.violazioni] == ["E-S0-30"]
    dettaglio = esito.violazioni[0].dettaglio
    assert "2 righe" in dettaglio and "T01" in dettaglio and "F05" in dettaglio
    assert "LIB005_F05.fq.gz" in dettaglio and "LIB001_T01.fq.gz" in dettaglio
    assert PRIMER_ITS in dettaglio and "341F" in dettaglio
    assert "copia della tabella" in dettaglio
    risultato = esegui_s0(config, solleva=False)
    eseguiti = {e.gate: e.eseguito for e in risultato.esiti}
    assert eseguiti["G04"] and not eseguiti["G05"]

    coerente = _due_ampliconi(tmp_path / "coerente", con_file_its=True)
    dati = sintetico.configurazione(coerente, tmp_path / "out2")
    dati["io"]["accession_regex"] = "LIB[0-9]+"
    giusto = esegui_gate("G04", g01_g15.Contesto(valida(dati)))
    assert giusto.superato and not giusto.avvisi and not giusto.violazioni


def test_senza_un_primer_in_testa_o_senza_dichiarazioni_il_confronto_non_si_fa(tmp_path):
    """
    **Obiettivo**: Verificare che il confronto fra dichiarazione e letture non
    giudichi alcuna riga quando le letture non portano un primer riconosciuto
    (gia' tolto) anche se la tabella dichiara un altro amplicone, e che la
    validazione non costruisca il confronto quando la tabella non ha le colonne
    dell'amplicone e del primer.

    **Razionale scientifico e sistemistico**: Un arresto fondato su un indizio
    che manca sarebbe un falso allarme su ogni dataset a cui i primer sono
    stati tolti prima del deposito: senza un primer in testa non c'e' nulla da
    confrontare, e il controllo deve tacere.
    """
    dataset = sintetico.genera(tmp_path / "dati", biologici=3, negativi=0, letture=60)
    config = valida(sintetico.configurazione(dataset, tmp_path / "out"))
    assert g01_g15.righe_e_file(config, analizza(config)) is None

    righe = _leggi(dataset.tabella)
    for riga in righe:
        riga[COLONNA_AMPLICONE] = "ITS"
        riga[COLONNA_PRIMER] = f"ITS1F, 5-{PRIMER_ITS}-3"
    _scrivi(dataset.tabella, righe)
    collegate, assenti = g01_g15.righe_e_file(config, analizza(config))
    assert len(collegate) == 3 and not assenti
    assert ampliconi.confronta(collegate) == ([], 0)
    assert esegui_gate("G04", g01_g15.Contesto(config)).superato


# --------------------------------------------------------------------------- #
# Righe di un altro amplicone                                                  #
# --------------------------------------------------------------------------- #


def test_init_dice_quante_righe_non_ha_potuto_confrontare_con_le_letture(tmp_path, capsys):
    """
    **Obiettivo**: Verificare che, quando una parte delle righe collegate a un
    file non si puo' confrontare con le letture (la riga dichiara l'amplicone
    in un modo che non si interpreta, e non dichiara un primer) mentre le
    altre risultano coerenti, il riepilogo di ``init`` dica quante sono e che
    il loro collegamento resta da verificare; e che con tutte le righe
    giudicate la nota non compaia.

    **Razionale scientifico e sistemistico**: Una riga non giudicata puo'
    essere collegata al file di un altro amplicone senza che nulla lo segnali:
    il silenzio, quando le altre righe sono coerenti, verrebbe letto come una
    verifica riuscita su tutta la tabella.
    """
    tutte = _due_ampliconi(tmp_path / "tutte", con_file_its=True)
    assert _init(tutte, tmp_path / "tutte.yaml") == 0
    assert "non si sono potute confrontare" not in capsys.readouterr().out

    dataset = _due_ampliconi(tmp_path / "parte", con_file_its=True)
    righe = _leggi(dataset.tabella)
    righe[1][COLONNA_AMPLICONE], righe[1][COLONNA_PRIMER] = "amplicone da definire", ""
    _scrivi(dataset.tabella, righe)
    assert _init(dataset, tmp_path / "parte.yaml") == 0
    console = " ".join(capsys.readouterr().out.split())
    assert "1 righe su 6 collegate a un file non si sono potute confrontare con le letture" \
        in console
    assert "resta da verificare" in console


def test_le_righe_da_tenere_si_dichiarano_e_le_altre_sono_escluse_con_un_avviso(tmp_path):
    """
    **Obiettivo**: Verificare che una tabella con righe di un altro amplicone,
    i cui file non sono nella cartella, sia respinta da G06 (``E-S0-06``) senza
    la regola; che con ``meta.assay_keep_column`` e ``meta.assay_keep_values``
    la validazione passi, G04 dichiari le righe escluse con ``E-S0-31``
    nominandole e l'inventario abbia i soli campioni tenuti; che dichiarare la
    colonna senza i valori sia respinto con ``E-G15-14`` e una colonna assente
    con ``E-S0-02``; e che senza la regola i due parametri non compaiano nella
    configurazione risolta, il cui digest non dipende dalla loro esistenza.

    **Razionale scientifico e sistemistico**: La tabella di assay di uno studio
    con piu' ampliconi descrive campioni che un'analisi 16S non tratta: vanno
    esclusi per dichiarazione, con un elenco che si possa leggere, non
    cancellando righe dai metadati originali ne' in silenzio.
    """
    dataset = _due_ampliconi(tmp_path)
    senza = valida(sintetico.configurazione(dataset, tmp_path / "out"))
    risultato = esegui_s0(senza, solleva=False)
    respinti = [v.codice for e in risultato.esiti for v in e.violazioni]
    assert respinti == ["E-S0-06"]
    mappa = risolvi(senza).come_mappa()
    assert "assay_keep_column" not in mappa["meta"] and "assay_keep_values" not in mappa["meta"]
    assert {"meta.assay_keep_column", "meta.assay_keep_values"} <= senza.chiavi_assenti()

    regola = {"assay_keep_column": COLONNA_AMPLICONE, "assay_keep_values": ["16s"]}
    con = valida(sintetico.configurazione(dataset, tmp_path / "out", meta=regola))
    assert risolvi(con).come_mappa()["meta"]["assay_keep_values"] == ["16s"]
    assert risolvi(con).digest != risolvi(senza).digest
    risultato = esegui_s0(con, solleva=False)
    assert risultato.superata
    g04 = next(e for e in risultato.esiti if e.gate == "G04")
    assert [a.codice for a in g04.avvisi] == ["E-S0-31"]
    assert "2 righe" in g04.avvisi[0].dettaglio
    assert "F05 (ITS)" in g04.avvisi[0].dettaglio and "F06 (ITS)" in g04.avvisi[0].dettaglio
    assert sorted(c.nome for c in risultato.inventario) == ["T01", "T02", "T03", "T04"]

    with pytest.raises(g01_g15.ErroreGate) as rifiuto:
        esegui_g15(sintetico.configurazione(
            dataset, tmp_path / "out", meta={"assay_keep_column": COLONNA_AMPLICONE}))
    assert rifiuto.value.codice == "E-G15-14"
    assente = valida(sintetico.configurazione(
        dataset, tmp_path / "out",
        meta={"assay_keep_column": "colonna che non c'e'", "assay_keep_values": ["16S"]}))
    codici = [v.codice for e in esegui_s0(assente, solleva=False).esiti for v in e.violazioni]
    assert codici == ["E-S0-02"]


def test_init_cerca_la_chiave_fra_le_righe_con_un_file_e_propone_le_righe_16s(tmp_path, capsys):
    """
    **Obiettivo**: Verificare che ``init``, su una tabella con piu' righe dei
    file perche' descrive anche un altro amplicone, trovi la chiave fra le sole
    righe che nominano un file della cartella, proponga di tenere le righe 16S
    con la colonna che dichiara l'amplicone, lo scriva nel riepilogo insieme
    alle righe che nominano un file assente, e che la configurazione con le
    proposte accettate superi i gate sui metadati escludendo le altre righe.

    **Razionale scientifico e sistemistico**: Senza cercare la chiave fra le
    righe che hanno un file, ``init`` non la troverebbe in nessuna tabella di
    uno studio con piu' ampliconi, e chi comincia dovrebbe scriverla a mano.
    """
    dataset = _due_ampliconi(tmp_path)
    uscita = tmp_path / "c.yaml"
    assert _init(dataset, uscita, "--accetta-proposte") == 0
    console = " ".join(capsys.readouterr().out.split())
    testo = uscita.read_text(encoding="utf-8")
    parametri = righe_dei_parametri(testo)
    assert parametri["meta.accession_column"] == (True, f'"{sintetico.COLONNA_FILE}"')
    assert parametri["meta.assay_keep_column"] == (True, f'"{COLONNA_AMPLICONE}"')
    assert parametri["meta.assay_keep_values"] == (True, '["16S"]')
    piatto = " ".join(testo.replace("#", " ").split())
    assert "fra le 4 righe su 6 che nominano un file della cartella" in piatto
    assert "proposto di tenere le sole 4 righe 16S" in console
    assert "2 righe della tabella nominano un file che non e' nella cartella" in console
    assert "F05 (LIB005_F05.fq.gz, ITS)" in console

    # Un dataset senza altri ampliconi non riceve i due parametri.
    semplice = sintetico.genera(tmp_path / "semplice", biologici=3, negativi=0, letture=40)
    altra = tmp_path / "s.yaml"
    assert _init(semplice, altra) == 0
    assert "assay_keep" not in altra.read_text(encoding="utf-8")


def test_init_non_propone_una_chiave_con_righe_collegate_al_file_sbagliato(tmp_path, capsys):
    """
    **Obiettivo**: Verificare che ``init``, su una tabella in cui una riga 16S
    e una dell'altro amplicone sono collegate ciascuna al file dell'altra, non
    proponga la chiave del campione, nomini nel riepilogo le righe incoerenti
    con cio' che dichiarano e cio' che le letture portano, e dica di correggere
    una copia della tabella.

    **Razionale scientifico e sistemistico**: Una chiave proposta su una
    tabella che collega i campioni ai file sbagliati darebbe una configurazione
    che la validazione ferma, o peggio che esegue: chi comincia deve sapere
    subito che il difetto e' nei metadati.
    """
    dataset = _due_ampliconi(tmp_path, scambiate=True)
    uscita = tmp_path / "c.yaml"
    assert _init(dataset, uscita, "--accetta-proposte") == 0
    console = " ".join(capsys.readouterr().out.split())
    parametri = righe_dei_parametri(uscita.read_text(encoding="utf-8"))
    assert parametri["meta.accession_column"] == (False, "")
    assert parametri["io.accession_regex"] == (False, "")
    assert "2 righe della tabella sono collegate a un file" in console
    assert "T01" in console and "F05" in console and PRIMER_ITS in console
    assert "correggi una copia della tabella" in console and "E-S0-30" in console


def test_fra_due_colonne_di_file_init_sceglie_quella_con_un_file_per_cella(tmp_path):
    """
    **Obiettivo**: Verificare che, quando una colonna riporta il file di ogni
    campione e un'altra elenca in ogni cella due file con altri nomi che
    contengono la stessa chiave, ``init`` rilevi la chiave nella prima e non
    resti senza scelta fra le due.

    **Razionale scientifico e sistemistico**: Le tabelle dei depositi hanno
    spesso piu' colonne di file (grezzi, filtrati, elaborati): leggere come
    elenco le celle con piu' file non deve togliere la precedenza alla colonna
    che nomina un file per campione.
    """
    dataset = sintetico.genera(tmp_path / "dati", biologici=3, negativi=0, letture=40)
    righe = _leggi(dataset.tabella)
    for riga in righe:
        radice = riga[sintetico.COLONNA_FILE].removesuffix(".fq.gz")
        riga["file_elaborati"] = f"{radice}_a.fastq.gz, {radice}_b.fastq.gz"
    _scrivi(dataset.tabella, righe)
    uscita = tmp_path / "c.yaml"
    assert _init(dataset, uscita) == 0
    parametri = righe_dei_parametri(uscita.read_text(encoding="utf-8"))
    assert parametri["meta.accession_column"] == (True, f'"{sintetico.COLONNA_FILE}"')


def test_una_cella_con_un_altro_nome_del_file_si_collega_per_chiave(tmp_path):
    """
    **Obiettivo**: Verificare che una riga la cui cella nomina il file con un
    nome diverso da quello nella cartella (un percorso davanti) si colleghi per
    chiave e non risulti con il file assente; che una riga 18S che dichiara un
    primer universale dell'elenco, collegata al proprio file che lo porta in
    testa, non sia incoerente; e che un primer del solo elenco non basti a
    dire incoerente una riga di un altro amplicone.

    **Razionale scientifico e sistemistico**: Un arresto su una tabella giusta
    costa piu' di un'incoerenza non vista: i primer universali amplificano
    anche il 18S, e il nome di un file nei metadati non e' sempre quello sul
    disco. Il confronto ferma solo dove la tabella stessa dice di chi e' il
    primer trovato.
    """
    dataset = sintetico.genera(
        tmp_path / "dati", biologici=3, negativi=0, letture=60, primer=sintetico.PRIMER_341F)
    righe = _leggi(dataset.tabella)
    presenti = {p.name: p for p in dataset.letture.iterdir()}
    for riga in righe:
        riga[sintetico.COLONNA_FILE] = "deposito/grezzi/" + riga[sintetico.COLONNA_FILE]
        riga[COLONNA_AMPLICONE] = "18S rRNA"
        riga[COLONNA_PRIMER] = "341F"
    per_chiave = {re.search(r"LIB[0-9]+", nome).group(0): p for nome, p in presenti.items()}
    collegate, assenti = ampliconi.collega(
        righe, sintetico.COLONNA_CAMPIONE, sintetico.COLONNA_FILE, presenti,
        COLONNA_AMPLICONE, COLONNA_PRIMER,
        lambda valore: per_chiave.get(re.search(r"LIB[0-9]+", valore).group(0)))
    assert len(collegate) == 3 and not assenti
    assert ampliconi.confronta(collegate) == ([], 3)
    # Senza la dichiarazione del primer: il file porta un primer che solo
    # l'elenco conosce, e la riga 18S non si giudica.
    for riga in righe:
        riga[COLONNA_PRIMER] = ""
    collegate, _ = ampliconi.collega(
        righe, sintetico.COLONNA_CAMPIONE, sintetico.COLONNA_FILE, presenti,
        COLONNA_AMPLICONE, COLONNA_PRIMER,
        lambda valore: per_chiave.get(re.search(r"LIB[0-9]+", valore).group(0)))
    assert ampliconi.confronta(collegate) == ([], 0)
    # Senza modo di collegare per chiave, il file nominato risulta assente.
    _, assenti = ampliconi.collega(
        righe, sintetico.COLONNA_CAMPIONE, sintetico.COLONNA_FILE, presenti,
        COLONNA_AMPLICONE, COLONNA_PRIMER)
    assert len(assenti) == 3


def test_init_avverte_dei_file_delle_righe_che_la_regola_esclude(tmp_path, capsys):
    """
    **Obiettivo**: Verificare che, con i file dei due ampliconi nella cartella
    e ogni riga collegata al proprio, ``init`` proponga le righe 16S e scriva
    nel riepilogo che i file delle righe escluse vanno tolti dalla cartella,
    nominandoli, perche' la validazione si fermerebbe su un file senza riga.

    **Razionale scientifico e sistemistico**: La regola toglie le righe, non i
    file: chi accetta la proposta deve sapere prima di ``validate`` che cosa
    resta da fare, non scoprirlo da un arresto.
    """
    dataset = _due_ampliconi(tmp_path, con_file_its=True)
    assert _init(dataset, tmp_path / "c.yaml", "--accetta-proposte") == 0
    console = " ".join(capsys.readouterr().out.split())
    assert "2 file della cartella appartengono alle righe che la regola esclude" in console
    assert "LIB005_F05.fq.gz" in console and "E-S0-06" in console
