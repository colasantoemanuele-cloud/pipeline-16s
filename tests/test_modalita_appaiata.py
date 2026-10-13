"""La modalita' appaiata: configurazione, ingressi e controlli che precedono il calcolo.

Verifica cio' che distingue un dataset con le due letture di ogni frammento
(``io.layout: paired-end``) da uno single-end, prima di qualunque fase in R:

* la configurazione (``amplicon16s/config``): la modalita' single-end dichiarata
  vale come omessa, il gruppo ``pair`` esiste solo in modalita' appaiata;
* l'abbinamento dei file e la sincronia delle coppie
  (``amplicon16s/io_layer/reads.py``, i gate di S0, la fase S1);
* l'orientamento delle due letture: una corsa con i due file scambiati, un
  file con i due versi mescolati, il primer reverse lasciato in testa;
* in modalita' single-end, un file con il nome di una seconda lettura;
* le soglie della modalita' appaiata nel meccanismo delle soglie
  (``amplicon16s/soglie``), il grafo delle fasi e ``init``.

I dataset sono sintetici (``tests/sintetico.py``, :func:`genera_coppie`): ogni
frammento ha lunghezza nota, e le due letture ne sono le estremita'. Gli arresti
di S1 avvengono prima del profilo in R: nessun test di questo modulo richiede i
pacchetti di calcolo. La catena completa e' in ``tests/test_catena_appaiata.py``.
"""

from __future__ import annotations

import math
import copy
import gzip
import json
import random
from collections import Counter
from pathlib import Path

import pytest
import sintetico
import yaml
from init_aiuti import righe_dei_parametri

from amplicon16s import interfaccia
from amplicon16s.cli import main
from amplicon16s.config import defaults
from amplicon16s.config.resolve import parametri_dichiarati, risolvi, scrivi_risolta
from amplicon16s.config.schema import ErroreConfigurazione, valida
from amplicon16s.errors.catalog import CATALOGO, Categoria
from amplicon16s.gates.g01_g15 import rifiuto_di_g15
from amplicon16s.init import rilevamento
from amplicon16s.init import vocabolario
from amplicon16s.io_layer import reads
from amplicon16s.io_layer.artifacts import Fase
from amplicon16s.logging.logger import chiudi
from amplicon16s.runner.executor import Conclusione, Esecutore
from amplicon16s.runner.graph import GRAFO, Passo
from amplicon16s.runner.project import ProjectRun
from amplicon16s.soglie import tabella


@pytest.fixture(autouse=True)
def uscite_pulite():
    """Chiude le uscite del log prima e dopo ogni test."""
    chiudi()
    yield
    chiudi()


def _config(dataset: sintetico.DatasetAppaiato, uscita: Path, **gruppi):
    """La configurazione appaiata del dataset sintetico, validata."""
    return valida(sintetico.configurazione_appaiata(dataset, uscita, **gruppi))


def _arresto(config, fino_a: Passo = Passo.S1) -> tuple[Conclusione, str | None, str]:
    """Esegue la catena fino alla fase indicata e restituisce la conclusione,
    il codice dell'arresto e il suo dettaglio.
    """
    esito = Esecutore(ProjectRun(config), fino_a=fino_a).esegui()
    chiudi()
    punto = esito.punto
    return esito.conclusione, punto.codice if punto else None, (punto.dettaglio or "") if punto else ""


def _record(percorso: Path) -> list[tuple[str, str, str]]:
    """I record di un FASTQ: nome, sequenza, qualita'."""
    with gzip.open(percorso, "rt", encoding="ascii") as file:
        righe = file.read().splitlines()
    return [(righe[i][1:], righe[i + 1], righe[i + 3]) for i in range(0, len(righe), 4)]


def _motivo_forward(lunghezza: int = 19) -> str:
    """Le prime basi delle sequenze del riferimento sintetico come espressione:
    la base dove quasi tutte concordano, un punto dove variano.
    """
    sequenze = [s for _, s in sintetico.sequenze_di_riferimento()]
    simboli = []
    for posizione in range(lunghezza):
        base, quante = Counter(s[posizione] for s in sequenze).most_common(1)[0]
        simboli.append(base if quante >= 0.9 * len(sequenze) else ".")
    return "".join(simboli)


# --------------------------------------------------------------------------- #
# Configurazione                                                               #
# --------------------------------------------------------------------------- #


def test_la_modalita_single_end_dichiarata_vale_come_omessa(tmp_path):
    """
    **Obiettivo**: una configurazione che dichiara ``io.layout: single-end`` ha
    lo stesso digest, la stessa impronta dei risultati, gli stessi parametri
    dichiarati e lo stesso ``resolved.yaml``, byte per byte, di quella che non
    dichiara la modalita'; in nessuna delle due la configurazione registrata
    nomina la modalita' ne' il gruppo ``pair``.

    **Razionale scientifico e sistemistico**: la modalita' appaiata non deve
    cambiare alcun artefatto di un'esecuzione single-end: una configurazione
    che non nomina la modalita' conserva il proprio digest, e resta la stessa
    quando ``init`` scrive la modalita' per esteso.
    """
    dataset = sintetico.genera(tmp_path / "dati", biologici=2, letture=30)
    omessa = sintetico.configurazione(dataset, tmp_path / "out")
    dichiarata = copy.deepcopy(omessa)
    dichiarata["io"]["layout"] = "single-end"
    a, b = risolvi(valida(omessa)), risolvi(valida(dichiarata))
    assert not a.config.appaiata and not b.config.appaiata
    assert a.digest == b.digest and a.impronta_risultati == b.impronta_risultati
    assert parametri_dichiarati(a.config) == parametri_dichiarati(b.config)
    testi = [scrivi_risolta(r, tmp_path / nome).read_text(encoding="utf-8")
             for nome, r in (("a", a), ("b", b))]
    assert testi[0] == testi[1]
    assert "layout" not in testi[0] and "pair" not in yaml.safe_load(testi[0])["parametri"]
    assert "pair" not in a.come_mappa() and "layout" not in a.come_mappa()["io"]


def test_il_gruppo_pair_vale_solo_in_modalita_appaiata(tmp_path):
    """
    **Obiettivo**: un parametro del gruppo ``pair`` dichiarato in modalita'
    single-end e' respinto con ``E-G15-15``, che lo nomina e dice come
    rimediare; in modalita' appaiata i parametri senza predefinito che mancano
    sono elencati da ``E-G15-10`` con la loro condizione; e una configurazione
    appaiata completa registra la modalita' e il gruppo.

    **Razionale scientifico e sistemistico**: un parametro accettato senza
    effetto farebbe credere che le letture reverse siano usate; un parametro
    mancante va chiesto prima di aprire un file.
    """
    dataset = sintetico.genera_coppie(tmp_path / "dati", biologici=2, coppie=20)
    completa = sintetico.configurazione_appaiata(dataset, tmp_path / "out")
    risolta = risolvi(valida(completa))
    assert risolta.config.appaiata
    assert risolta.come_mappa()["io"]["layout"] == "paired-end"
    assert risolta.come_mappa()["pair"]["min_overlap"] == defaults.PAIR_MIN_OVERLAP
    assert risolta.derivati.asv_len_min == completa["pair"]["len_min"]
    assert risolta.derivati.asv_len_max == completa["pair"]["len_max"]
    assert "io.layout" in parametri_dichiarati(risolta.config)

    senza_modalita = copy.deepcopy(completa)
    del senza_modalita["io"]["layout"]
    with pytest.raises(ErroreConfigurazione) as rifiuto:
        valida(senza_modalita)
    assert [v.codice for v in rifiuto_di_g15(rifiuto.value).violazioni] == ["E-G15-15"]
    assert "pair.file_markers" in str(rifiuto.value) and "paired-end" in str(rifiuto.value)

    incompleta = copy.deepcopy(completa)
    for nome in ("truncLen_reverse", "len_min"):
        del incompleta["pair"][nome]
    incompleta["pair"]["min_merged_frac_stop"] = None
    with pytest.raises(ErroreConfigurazione) as rifiuto:
        valida(incompleta)
    testo = str(rifiuto.value)
    assert [v.codice for v in rifiuto_di_g15(rifiuto.value).violazioni] == ["E-G15-10"]
    for chiave in ("pair.truncLen_reverse", "pair.len_min", "pair.min_merged_frac_stop"):
        assert f"{chiave} perche' io.layout e' paired-end" in testo, chiave
    # Il primer e il motivo della reverse, e la lunghezza attesa, si dichiarano
    # anche nulli: nullo non e' mancante.
    for nome in ("primer_sequence_reverse", "conserved_motif_reverse", "expected_merged_len"):
        completa["pair"][nome] = None
    assert valida(completa).pair.expected_merged_len is None
    assert CATALOGO["E-G15-15"].categoria is Categoria.REVISIONE_UMANA


@pytest.mark.parametrize(("modifica", "atteso"), [
    ({"file_markers": ["_R1", "_R1"]}, "i due marcatori sono uguali"),
    ({"truncLen_reverse": 5, "trimLeft_reverse": 5}, "deve essere minore di truncLen_reverse"),
    ({"truncLen_shortfall_warn_reverse": 100}, "truncLen_shortfall_warn_reverse"),
    ({"min_merged_frac_warn": 0.1, "min_merged_frac_stop": 0.5}, "non puo' stare sotto"),
    ({"len_min": 200, "len_max": 100}, "l'intervallo delle lunghezze ammesse sarebbe vuoto"),
    ({"trim_overhang": False}, "trim_overhang: false non e' realizzato"),
])
def test_le_combinazioni_impossibili_dei_parametri_delle_coppie_sono_respinte(
    tmp_path, modifica, atteso
):
    """
    **Obiettivo**: ogni combinazione dei parametri di ``pair`` che non puo'
    valere (marcatori uguali, taglio non minore del troncamento, tolleranza non
    minore del troncamento, avviso sotto l'arresto, intervallo vuoto, code non
    tolte) e' respinta con ``E-G15-15`` e un messaggio che dice quale.

    **Razionale scientifico e sistemistico**: sono errori che si scoprirebbero
    dopo il filtro e l'inferenza di entrambe le letture, o mai: una soglia
    d'avviso sotto quella d'arresto non avviserebbe mai.
    """
    dataset = sintetico.genera_coppie(tmp_path / "dati", biologici=2, coppie=20)
    dati = sintetico.configurazione_appaiata(dataset, tmp_path / "out", pair=modifica)
    with pytest.raises(ErroreConfigurazione) as rifiuto:
        valida(dati)
    assert atteso in str(rifiuto.value)
    assert [v.codice for v in rifiuto_di_g15(rifiuto.value).violazioni] == ["E-G15-15"]


# --------------------------------------------------------------------------- #
# Abbinamento dei file                                                         #
# --------------------------------------------------------------------------- #


def test_i_file_si_abbinano_per_l_ultima_occorrenza_del_marcatore():
    """
    **Obiettivo**: :func:`abbina_file` riconosce il file forward e il suo
    reverse sostituendo l'ultima occorrenza del primo marcatore, anche quando
    il nome del campione contiene il marcatore e con un numero di blocco dopo;
    elenca i file senza compagno; e :func:`compagni_per_nome` riconosce i due
    file di uno stesso campione con le convenzioni note.

    **Razionale scientifico e sistemistico**: il marcatore di coppia sta in
    fondo al nome; un campione che lo contiene nel proprio nome non deve essere
    scambiato per un file reverse, ne' restare senza compagno.
    """
    esito = reads.abbina_file(
        ["x_2_1.fq", "x_2_2.fq", "y_1_001.fq", "y_2_001.fq"], ("_1", "_2"))
    assert esito.coppie == {"x_2_1.fq": "x_2_2.fq", "y_1_001.fq": "y_2_001.fq"}
    assert esito.completo
    # Due campioni che si chiamano come i due marcatori: quattro file, due coppie.
    for marcatori, nomi in ((("_1", "_2"), ["X_1_1.fq.gz", "X_1_2.fq.gz", "X_2_1.fq.gz",
                                             "X_2_2.fq.gz"]),
                            (("_R1", "_R2"), ["A_R1_R1_001.fastq", "A_R1_R2_001.fastq",
                                              "A_R2_R1_001.fastq", "A_R2_R2_001.fastq"])):
        esito = reads.abbina_file(nomi, marcatori)
        assert esito.completo and esito.coppie == {nomi[0]: nomi[1], nomi[2]: nomi[3]}, nomi
    # Il marcatore conta solo in fondo al nome.
    assert not reads.abbina_file(["a_R1_x.fq", "a_R2_x.fq"], ("_R1", "_R2")).coppie
    esito = reads.abbina_file(["a_R1.fq", "a_R2.fq", "b_R1.fq", "note.fq"], ("_R1", "_R2"))
    assert esito.coppie == {"a_R1.fq": "a_R2.fq"}
    assert esito.senza_compagno == ("b_R1.fq", "note.fq") and not esito.completo
    assert reads.compagni_per_nome(["s.1.fq", "s.2.fq"]) == (".1", ".2")
    assert reads.compagni_per_nome(["s_R1_001.fq", "s_R2_001.fq"]) == ("_R1", "_R2")
    assert reads.compagni_per_nome(["s_rep1.fq", "s_rep2.fq"]) is None
    assert reads.compagni_per_nome(["s_R1.fq", "s_R2.fq", "s_R3.fq"]) is None


def test_i_file_con_il_suffisso_raw_dopo_il_marcatore_si_abbinano():
    """
    **Obiettivo**: :func:`abbina_file` e :func:`compagni_per_nome` riconoscono
    le due letture di un campione quando il marcatore e' seguito da ``_raw``
    prima dell'estensione, senza distinzione fra maiuscole e minuscole;
    :func:`nome_da_seconda_lettura` riconosce il file reverse con quel nome; un
    altro suffisso dopo il marcatore (``_trimmed``, ``_L001``) resta non
    riconosciuto, e il solo ``_raw`` senza marcatore non fa di un file una
    lettura di coppia.

    **Razionale scientifico e sistemistico**: alcuni archivi pubblici
    distribuiscono le letture grezze con ``_raw`` dopo il marcatore della
    lettura. E' una forma in piu' della stessa regola, ancorata alla fine del
    nome: i nomi gia' riconosciuti si abbinano come prima, e le forme che la
    regola non nomina non diventano coppie per somiglianza.
    """
    nomi = ["a_R1_raw.fastq.gz", "a_R2_raw.fastq.gz", "b_R1_raw.fastq.gz", "b_R2_raw.fastq.gz"]
    esito = reads.abbina_file(nomi, ("_R1", "_R2"))
    assert esito.completo
    assert esito.coppie == {nomi[0]: nomi[1], nomi[2]: nomi[3]}
    assert reads.compagni_per_nome(["s_R1_raw.fq", "s_R2_raw.fq"]) == ("_R1", "_R2")
    assert reads.compagni_per_nome(["s_1_raw.fq.gz", "s_2_raw.fq.gz"]) == ("_1", "_2")
    assert reads.compagni_per_nome(["s_R1_RAW.fastq", "s_R2_RAW.fastq"]) == ("_R1", "_R2")
    assert reads.nome_da_seconda_lettura("a_R2_raw.fastq.gz", [])
    assert reads.nome_da_seconda_lettura("s_2_raw.fq.gz", ["s_1_raw.fq.gz", "s_2_raw.fq.gz"])
    assert not reads.nome_da_seconda_lettura("s_2_raw.fq.gz", ["s_2_raw.fq.gz"])
    assert not reads.nome_da_seconda_lettura("corsa2_raw.fastq.gz", [])
    for dopo in ("_trimmed", "_filtered", "_L001", "_raw_001"):
        coppia = [f"a_R1{dopo}.fastq.gz", f"a_R2{dopo}.fastq.gz"]
        assert not reads.abbina_file(coppia, ("_R1", "_R2")).coppie, dopo
        assert reads.compagni_per_nome(coppia) is None, dopo
        assert not reads.nome_da_seconda_lettura(coppia[1], coppia), dopo
    # I nomi senza suffisso e con il numero di blocco si abbinano come prima.
    assert reads.compagni_per_nome(["s_R1.fq", "s_R2.fq"]) == ("_R1", "_R2")
    assert reads.compagni_per_nome(["s_R1_001.fq", "s_R2_001.fq"]) == ("_R1", "_R2")


def test_un_file_senza_compagno_ferma_la_validazione(tmp_path):
    """
    **Obiettivo**: in modalita' appaiata, tolto il file reverse di un campione,
    la validazione si ferma con ``E-S0-23``, nomina il file rimasto solo e i
    marcatori con cui lo ha cercato; con tutti i file al loro posto passa, i gate
    dichiarano in ``gates.json`` che il layout verificato e' quello appaiato, e
    l'inventario porta per ogni campione il file reverse.

    **Razionale scientifico e sistemistico**: un campione con una sola lettura
    non si puo' fondere; scoprirlo dopo il filtro costerebbe il calcolo di
    tutti gli altri.
    """
    dataset = sintetico.genera_coppie(tmp_path / "dati", biologici=3, coppie=40)
    conclusione, codice, _ = _arresto(_config(dataset, tmp_path / "buono"), Passo.S0)
    assert (conclusione, codice) == (Conclusione.COMPLETATA, None)
    validazione = tmp_path / "buono" / Fase.INPUT_VALIDATION.value
    gate = {g["gate"]: g for g in json.loads((validazione / "gates.json").read_text())["gate"]}
    assert gate["G05"]["descrizione"].endswith("ogni file ha il suo compagno")
    assert gate["G07"]["descrizione"].startswith("il layout e' appaiato")
    inventario = (validazione / "crosswalk.tsv").read_text(encoding="utf-8").splitlines()
    assert inventario[0].split("\t")[-1] == "file_reverse"
    assert all(riga.split("\t")[-1].endswith("_R2.fq.gz") for riga in inventario[1:])
    assert (validazione / "letture_ispezionate_reverse.tsv").is_file()

    solo = dataset.coppie["LIB002"][0]
    dataset.coppie["LIB002"][1].unlink()
    conclusione, codice, dettaglio = _arresto(_config(dataset, tmp_path / "out"), Passo.S0)
    assert (conclusione, codice) == (Conclusione.ARRESTATA, "E-S0-23")
    assert solo.name in dettaglio and "senza compagno" in dettaglio
    assert "'_R1'" in dettaglio and "'_R2'" in dettaglio
    assert "pair.file_markers" in CATALOGO["E-S0-23"].azione


def test_le_due_letture_in_modalita_single_end_propongono_la_modalita_appaiata(tmp_path):
    """
    **Obiettivo**: una cartella con ``_R1`` e ``_R2`` di ogni campione, data a
    una configurazione single-end, ferma la validazione con ``E-S0-24`` (e non
    con il codice della chiave ripetuta), e il messaggio nomina i due file e
    scrive la modalita' e i marcatori da dichiarare.

    **Razionale scientifico e sistemistico**: chi ha un dataset appaiato e
    legge «accession ripetuto» cerca un errore nei metadati che non c'e': il
    messaggio deve dire che cosa sono quei file e le due strade possibili.
    """
    dataset = sintetico.genera_coppie(tmp_path / "dati", biologici=2, coppie=40)
    config = valida(sintetico.configurazione(dataset, tmp_path / "out"))
    conclusione, codice, dettaglio = _arresto(config, Passo.S0)
    assert (conclusione, codice) == (Conclusione.ARRESTATA, "E-S0-24")
    assert "LIB001" in dettaglio and "_R1 e _R2" in dettaglio
    assert 'io.layout: "paired-end"' in dettaglio
    assert 'pair.file_markers: ["_R1", "_R2"]' in dettaglio
    assert "E-S0-05" not in dettaglio
    assert "paired-end" in CATALOGO["E-S0-24"].azione


# --------------------------------------------------------------------------- #
# Sincronia delle coppie                                                       #
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize("intestazioni", ["illumina", "barra", "uguali", "numerate"])
def test_la_sincronia_si_verifica_comunque_le_coppie_siano_nominate(tmp_path, intestazioni):
    """
    **Obiettivo**: :func:`verifica_sincronia` accetta due file in cui la
    lettura n dell'uno e dell'altro sono la stessa coppia, con ognuno dei modi
    di nominarle (stesso identificativo con ``1:N`` e ``2:N`` nel commento,
    ``/1`` e ``/2``, lo stesso nome, numeri progressivi di un archivio a
    distanza costante) e ne conta le coppie; tolta una lettura a meta' del file
    reverse, dichiara il posto esatto in cui i nomi smettono di corrispondere.

    **Razionale scientifico e sistemistico**: la fusione appaia le letture per
    posizione: un file filtrato o troncato senza l'altro darebbe, da li' in
    poi, coppie di frammenti diversi, senza alcun errore. I depositi nominano
    le coppie in modi diversi, e nessuno di essi deve far respingere due file
    in ordine.
    """
    dataset = sintetico.genera_coppie(
        tmp_path / "dati", biologici=1, coppie=60, intestazioni=intestazioni)
    forward, reverse = dataset.coppie["LIB001"]
    assert reads.verifica_sincronia(forward, reverse) == reads.Sincronia(60)
    record = _record(reverse)
    sintetico._scrivi_fastq(reverse, record[:30] + record[31:])
    esito = reads.verifica_sincronia(forward, reverse)
    assert esito.coppie == 30 and esito.difetto is not None
    assert "alla lettura 31 i nomi non corrispondono" in esito.difetto


def test_due_file_numerati_sfasati_di_poche_letture_non_sono_in_sincronia(tmp_path):
    """
    **Obiettivo**: con letture nominate da numeri progressivi,
    :func:`verifica_sincronia` accetta due file in cui le seconde letture sono
    numerate dopo tutte le prime (lo scarto e' almeno il numero delle letture,
    anche quando i file sono le prime letture di un archivio piu' grande), e
    respinge due file con lo stesso numero di letture sfasati di qualche
    lettura (il forward senza le ultime, il reverse senza le prime), in cui lo
    scarto fra i numeri e' costante ma minore delle letture.

    **Razionale scientifico e sistemistico**: uno scarto costante da solo non
    prova la sincronia: e' anche quello di due file tagliati in punti diversi,
    che darebbero coppie di frammenti diversi a ogni riga, quasi nessuna
    fusione in quel campione e nessun arresto.
    """
    def scrivi(percorso: Path, numeri) -> Path:
        sintetico._scrivi_fastq(percorso, [(f"ARCH.{n} {n}", "ACGT" * 5, "I" * 20) for n in numeri])
        return percorso

    archivio = (scrivi(tmp_path / "a_R1.fq.gz", range(1, 201)),
                scrivi(tmp_path / "a_R2.fq.gz", range(5001, 5201)))
    assert reads.verifica_sincronia(*archivio) == reads.Sincronia(200)
    sfasati = (scrivi(tmp_path / "b_R1.fq.gz", range(1, 191)),
               scrivi(tmp_path / "b_R2.fq.gz", range(11, 201)))
    esito = reads.verifica_sincronia(*sfasati)
    assert esito.coppie == 190 and esito.difetto is not None
    assert "a distanza 10, minore delle 190 letture" in esito.difetto
    assert "sfasati" in esito.difetto


def test_una_coppia_fuori_posto_e_un_numero_diverso_di_letture_fermano_s1(tmp_path):
    """
    **Obiettivo**: in modalita' appaiata, scambiate due letture a meta' del
    file reverse di un campione, S1 si ferma con ``E-S1-08`` nominando il
    campione, i due file e la lettura a cui i nomi non corrispondono; tolta
    invece l'ultima lettura, si ferma con lo stesso codice dicendo che i due
    file non hanno lo stesso numero di letture.

    **Razionale scientifico e sistemistico**: con questa verifica, fatta su
    tutte le letture prima di ogni calcolo, il filtro appaiato puo' scrivere i
    due file senza confrontarne gli identificativi (``matchIDs = FALSE``).
    """
    dataset = sintetico.genera_coppie(tmp_path / "dati", biologici=3, coppie=200)
    forward, reverse = dataset.coppie["LIB002"]
    record = _record(reverse)
    scambiati = list(record)
    scambiati[99], scambiati[100] = scambiati[100], scambiati[99]
    sintetico._scrivi_fastq(reverse, scambiati)
    conclusione, codice, dettaglio = _arresto(_config(dataset, tmp_path / "scambio"))
    assert (conclusione, codice) == (Conclusione.ARRESTATA, "E-S1-08")
    assert "1 campioni su 3" in dettaglio and "LIB002" in dettaglio
    assert forward.name in dettaglio and reverse.name in dettaglio
    assert "alla lettura 100 i nomi non corrispondono" in dettaglio

    sintetico._scrivi_fastq(reverse, record[:-1])
    conclusione, codice, dettaglio = _arresto(_config(dataset, tmp_path / "corto"))
    assert (conclusione, codice) == (Conclusione.ARRESTATA, "E-S1-08")
    assert "il file reverse finisce dopo 199 letture" in dettaglio
    assert "non hanno lo stesso numero di letture" in dettaglio
    assert CATALOGO["E-S1-08"].categoria is Categoria.REVISIONE_UMANA


# --------------------------------------------------------------------------- #
# Orientamento                                                                 #
# --------------------------------------------------------------------------- #


def test_una_corsa_invertita_si_riconosce_dal_contenuto(tmp_path):
    """
    **Obiettivo**: scambiati i due file di una corsa in un dataset le cui
    intestazioni non distinguono le due letture, la validazione passa e S1 si
    ferma con ``E-S1-09``: nomina la corsa e i suoi due file, riporta le
    letture con il segnale delle reverse nei due file, e dice di scambiarne i
    nomi. Le altre corse, in ordine, non sono nominate.

    **Razionale scientifico e sistemistico**: con i file scambiati la lettura
    forward verrebbe troncata, filtrata e modellata con i parametri della
    reverse: la fusione riuscirebbe, e le varianti sarebbero il complemento
    inverso di quelle degli altri campioni. Nei depositi succede, e nessun
    segno nei nomi lo dice.
    """
    dataset = sintetico.genera_coppie(
        tmp_path / "dati", biologici=4, coppie=600, intestazioni="uguali")
    forward, reverse = dataset.coppie["LIB003"]
    temporaneo = forward.with_name("scambio.tmp")
    forward.rename(temporaneo)
    reverse.rename(forward)
    temporaneo.rename(reverse)
    conclusione, codice, _ = _arresto(_config(dataset, tmp_path / "s0"), Passo.S0)
    assert (conclusione, codice) == (Conclusione.COMPLETATA, None)
    conclusione, codice, dettaglio = _arresto(_config(dataset, tmp_path / "out"))
    assert (conclusione, codice) == (Conclusione.ARRESTATA, "E-S1-09")
    assert "1 corse hanno i due file scambiati" in dettaglio
    assert "LIB003" in dettaglio and forward.name in dettaglio and reverse.name in dettaglio
    assert "il segnale delle letture reverse e' in 0 letture su 600 del file reverse" in dettaglio
    assert "Scambia i nomi dei due file" in dettaglio
    assert "LIB001" not in dettaglio and "LIB002" not in dettaglio and "LIB004" not in dettaglio
    assert "scambia i nomi dei due file" in CATALOGO["E-S1-09"].azione


def test_il_giudizio_sulla_corsa_invertita_chiede_una_prova_nei_due_file():
    """
    **Obiettivo**: :func:`corsa_invertita` dichiara lo scambio solo se, per
    ogni segnale cercato, il file che dovrebbe portarlo ne e' privo su almeno
    500 letture e il compagno lo porta in almeno 30; non lo dichiara per un
    file povero del segnale il cui compagno non lo porta, ne' sotto le 500
    letture, ne' senza un segnale cercato in entrambi i file.

    **Razionale scientifico e sistemistico**: un campione senza il bersaglio ha
    poche letture con il segnale in entrambi i file, e non e' una corsa
    invertita: lo scambio si dichiara su una prova positiva, non su
    un'assenza.
    """
    def statistiche(letture: int, propri: int | None, altrui: int | None):
        return reads.StatisticheFile(
            "f.fq.gz", letture, 100, 100, 0, 0, esaurito=True,
            segnale_per_tratto=None if propri is None else ((propri, letture),),
            con_altro_segnale=altrui)

    assert reads.corsa_invertita(statistiche(600, 0, 580), statistiche(600, 0, 590)) is not None
    assert reads.corsa_invertita(statistiche(600, 0, 0), statistiche(600, 590, 0)) is None
    assert reads.corsa_invertita(statistiche(600, 590, 0), statistiche(600, 585, 0)) is None
    assert reads.corsa_invertita(statistiche(499, 0, 480), statistiche(499, 0, 480)) is None
    assert reads.corsa_invertita(statistiche(600, 0, 29), statistiche(600, 0, 590)) is None
    # Un solo segnale cercato (il motivo forward e' nullo): basta quello.
    solo_reverse = reads.corsa_invertita(statistiche(600, None, 590), statistiche(600, 0, None))
    assert solo_reverse is not None and "file forward" in solo_reverse
    assert reads.corsa_invertita(statistiche(600, None, None), statistiche(600, None, None)) is None


def test_una_corsa_invertita_si_riconosce_dalle_intestazioni(tmp_path):
    """
    **Obiettivo**: scambiati i due file di una corsa le cui letture portano il
    marcatore di coppia (``1:N`` e ``2:N`` nel commento), la validazione si
    ferma gia' sulle prime letture con ``E-S0-26``, nominando i due file e la
    lettura che ciascuno dichiara; un file che contiene entrambe le letture di
    ogni coppia e' respinto con il codice della modalita' appaiata
    (``E-S0-27``), che dice di separarlo in due file.

    **Razionale scientifico e sistemistico**: in modalita' appaiata il
    marcatore uniforme e' atteso, 1 nel file forward e 2 nel reverse: quello
    opposto e' il segno di uno scambio, i due insieme di un file che non e'
    stato separato.
    """
    dataset = sintetico.genera_coppie(tmp_path / "dati", biologici=3, coppie=60)
    forward, reverse = dataset.coppie["LIB001"]
    record_forward, record_reverse = _record(forward), _record(reverse)
    sintetico._scrivi_fastq(forward, record_reverse)
    sintetico._scrivi_fastq(reverse, record_forward)
    conclusione, codice, dettaglio = _arresto(_config(dataset, tmp_path / "out"), Passo.S0)
    assert (conclusione, codice) == (Conclusione.ARRESTATA, "E-S0-26")
    assert forward.name in dettaglio and reverse.name in dettaglio
    assert "marcano la seconda lettura" in dettaglio and "marcano la prima lettura" in dettaglio
    assert "scambiane i nomi" in dettaglio

    intercalati = [r for coppia in zip(record_forward, record_reverse) for r in coppia]
    sintetico._scrivi_fastq(forward, intercalati)
    sintetico._scrivi_fastq(reverse, record_reverse)
    conclusione, codice, dettaglio = _arresto(_config(dataset, tmp_path / "due"), Passo.S0)
    assert (conclusione, codice) == (Conclusione.ARRESTATA, "E-S0-27")
    assert "le due letture di ogni coppia nello stesso file" in dettaglio
    assert "single-end" not in CATALOGO["E-S0-27"].azione
    assert "separa" in CATALOGO["E-S0-27"].azione


def test_un_file_con_le_coppie_intere_oltre_le_prime_letture_ferma_s1(tmp_path):
    """
    **Obiettivo**: in modalita' appaiata, un file forward che dopo le prime
    letture (quelle che la validazione ispeziona) contiene le due letture di
    ogni coppia supera la validazione e ferma S1 con ``E-S1-11``, che nomina
    il file e non suggerisce di passare a un'altra pipeline.

    **Razionale scientifico e sistemistico**: la validazione guarda le prime
    letture per durare poco; un file composto male piu' avanti si riconosce
    solo contando i segni di coppia su tutte le letture, e il rimedio in
    modalita' appaiata e' separare le due letture, non scartare le seconde.
    """
    dataset = sintetico.genera_coppie(tmp_path / "dati", biologici=3, coppie=200)
    forward, reverse = dataset.coppie["LIB002"]
    record_forward, record_reverse = _record(forward), _record(reverse)
    misto = record_forward[:40] + [
        r for coppia in zip(record_forward[40:], record_reverse[40:]) for r in coppia]
    sintetico._scrivi_fastq(forward, misto)
    config = _config(dataset, tmp_path / "out", qc={"head_reads": 30})
    conclusione, codice, dettaglio = _arresto(config)
    assert (conclusione, codice) == (Conclusione.ARRESTATA, "E-S1-11")
    assert forward.name in dettaglio and "LIB001" not in dettaglio
    assert "contengono le due letture di ogni coppia" in dettaglio
    assert "single-end" not in CATALOGO["E-S1-11"].azione
    assert "letture singole" not in CATALOGO["E-S1-11"].azione


def test_un_file_con_i_due_versi_mescolati_ferma_s1(tmp_path):
    """
    **Obiettivo**: in un dataset che porta i due primer nelle letture (i due
    segnali sono quindi entrambi definiti), scambiate fra i due file le letture
    di meta' delle coppie di un campione, scelte a caso, S1 si ferma con
    ``E-S1-10``, nomina il file e riporta le due quote; una corsa in ordine non
    e' nominata.

    **Razionale scientifico e sistemistico**: una libreria costruita con gli
    adattatori nei due versi ha in ogni file frammenti di entrambi gli
    orientamenti: le coppie si fondono tutte, e ogni variante esce due volte,
    una come complemento inverso dell'altra. Non e' una corsa invertita (i
    due file sono uguali) ne' un file a blocchi (le letture sono mescolate).
    """
    dataset = sintetico.genera_coppie(
        tmp_path / "dati", biologici=3, coppie=600, intestazioni="uguali",
        primer=sintetico.PRIMER_341F, primer_reverse=sintetico.PRIMER_805R, lunghezza=120)
    forward, reverse = dataset.coppie["LIB002"]
    record_forward, record_reverse = _record(forward), _record(reverse)
    caso = random.Random(7)
    for i in range(len(record_forward)):
        if caso.random() < 0.5:
            record_forward[i], record_reverse[i] = record_reverse[i], record_forward[i]
    sintetico._scrivi_fastq(forward, record_forward)
    sintetico._scrivi_fastq(reverse, record_reverse)
    conclusione, codice, dettaglio = _arresto(_config(dataset, tmp_path / "out"))
    assert (conclusione, codice) == (Conclusione.ARRESTATA, "E-S1-10")
    assert "1 campioni hanno un file con le letture dei due versi" in dettaglio
    assert "LIB002" in dettaglio and "LIB001" not in dettaglio
    assert "portano il segnale del verso atteso" in dettaglio
    assert CATALOGO["E-S1-10"].categoria is Categoria.REVISIONE_UMANA


def test_il_giudizio_sui_versi_mescolati_chiede_una_quota_rilevante_di_entrambi():
    """
    **Obiettivo**: :func:`versi_mescolati` segnala un file in cui almeno 30
    letture portano il segnale atteso e almeno 30 quello dell'altro verso, con
    la quota minore non sotto un decimo della maggiore; non segnala un file con
    poche letture dell'altro verso, ne' uno in cui uno dei due segnali non e'
    stato cercato.

    **Razionale scientifico e sistemistico**: qualche lettura dell'altro verso
    c'e' in ogni corsa (attribuzioni sbagliate fra i campioni): la regola deve
    fermare una libreria mescolata, non un file con l'uno per cento di
    letture estranee.
    """
    def statistiche(letture: int, propri: int | None, altrui: int | None):
        return reads.StatisticheFile(
            "f.fq.gz", letture, 100, 100, 0, 0, esaurito=True,
            segnale_per_tratto=None if propri is None else ((propri, letture),),
            con_altro_segnale=altrui)

    assert reads.versi_mescolati(statistiche(1000, 480, 470)) is not None
    assert reads.versi_mescolati(statistiche(1000, 900, 90)) is not None      # un decimo
    assert reads.versi_mescolati(statistiche(1000, 900, 89)) is None
    assert reads.versi_mescolati(statistiche(1000, 950, 29)) is None          # meno di 30
    assert reads.versi_mescolati(statistiche(1000, 100, 900)) is not None     # quasi tutto inverso
    assert reads.versi_mescolati(statistiche(1000, None, 500)) is None
    assert reads.versi_mescolati(statistiche(1000, 500, None)) is None


def test_il_primer_reverse_in_testa_senza_taglio_ferma_la_validazione(tmp_path):
    """
    **Obiettivo**: con le letture reverse che iniziano con il primer reverse e
    ``pair.trimLeft_reverse`` a zero, la validazione si ferma con ``E-S0-25``,
    nomina il primer e i file, e suggerisce il taglio pari alla sua lunghezza;
    con il taglio dichiarato passa.

    **Razionale scientifico e sistemistico**: un primer degenerato lasciato in
    testa alla lettura reverse finisce in coda alle sequenze fuse, e moltiplica
    le varianti quante sono le sue forme.
    """
    dataset = sintetico.genera_coppie(
        tmp_path / "dati", biologici=2, coppie=60, lunghezza=120,
        primer=sintetico.PRIMER_341F, primer_reverse=sintetico.PRIMER_805R)
    conclusione, codice, _ = _arresto(_config(dataset, tmp_path / "con_taglio"), Passo.S0)
    assert (conclusione, codice) == (Conclusione.COMPLETATA, None)
    config = _config(dataset, tmp_path / "out", pair={
        "trimLeft_reverse": 0, "conserved_motif_reverse": None})
    conclusione, codice, dettaglio = _arresto(config, Passo.S0)
    assert (conclusione, codice) == (Conclusione.ARRESTATA, "E-S0-25")
    assert sintetico.PRIMER_805R in dettaglio and "_R2.fq.gz" in dettaglio
    assert f"Imposta pair.trimLeft_reverse a {len(sintetico.PRIMER_805R)}" in dettaglio


# --------------------------------------------------------------------------- #
# Modalita' single-end: un file con il nome di una seconda lettura             #
# --------------------------------------------------------------------------- #


def _single_end_con_un_nome_sospetto(
    tmp_path: Path, suffisso: str, *, inverse: bool = False, letture: int = 200,
    negativo: bool = False, motivo: str | None = "noto",
):
    """Un dataset single-end in cui il file di un campione finisce con il
    suffisso indicato; a richiesta contiene letture inverse, o e' di un
    controllo negativo.
    """
    dataset = sintetico.genera(
        tmp_path / "dati", biologici=3, negativi=1, letture=letture, letture_negativi=letture)
    chiave = "LIB004" if negativo else "LIB002"
    (file,) = dataset.letture.glob(f"{chiave}_*.fq.gz")
    nuovo = file.with_name(f"{chiave}{suffisso}.fq.gz")
    record = _record(file)
    if inverse:
        record = [(n, sintetico.complemento_inverso(s), q[::-1]) for n, s, q in record]
    sintetico._scrivi_fastq(nuovo, record)
    file.unlink()
    testo = dataset.tabella.read_text(encoding="utf-8").replace(file.name, nuovo.name)
    dataset.tabella.write_text(testo, encoding="utf-8")
    qc = {"conserved_motif": _motivo_forward() if motivo else None, "min_motif_frac": 0.0}
    return valida(sintetico.configurazione(dataset, tmp_path / "out", qc=qc)), nuovo


@pytest.mark.parametrize("suffisso", ["_R2", "_2"])
def test_un_nome_da_seconda_lettura_non_basta_se_il_contenuto_si_puo_giudicare(
    tmp_path, suffisso
):
    """
    **Obiettivo**: in modalita' single-end, un file di letture forward di un
    campione biologico il cui nome finisce con ``_R2`` o con ``_2``, senza un
    compagno nella cartella, non e' respinto per il nome quando c'e' un segnale
    da cercare: la validazione passa.

    **Razionale scientifico e sistemistico**: un campione puo' chiamarsi
    davvero cosi'. Se il contenuto del file si puo' giudicare, e' il contenuto
    a dire se le letture sono forward, non il nome.
    """
    config, _ = _single_end_con_un_nome_sospetto(tmp_path, suffisso)
    conclusione, codice, _ = _arresto(config, Passo.S0)
    assert (conclusione, codice) == (Conclusione.COMPLETATA, None)


def test_un_file_sospetto_di_letture_inverse_e_respinto_anche_sotto_le_cinquecento_letture(
    tmp_path
):
    """
    **Obiettivo**: un file il cui nome finisce con ``_R2``, se contiene
    letture inverse, supera la validazione (che non lo respinge per il nome) e
    ferma S1 con ``E-S1-07`` anche con 200 letture, meno delle 500 chieste a
    un file qualunque; il giudizio vale da
    :data:`LETTURE_MINIME_PER_UN_NOME_SOSPETTO`. Il nome e' sospetto con la
    stessa regola per la validazione e per S1: ``R2`` da solo, il solo ``2``
    soltanto se nella cartella c'e' il compagno con ``1``.

    **Razionale scientifico e sistemistico**: togliere il rifiuto per nome non
    deve lasciar passare le seconde letture: per i file che il nome rende
    sospetti il controllo sul contenuto vale da un numero minore di letture.
    Un campione numerato ``_2`` senza un ``_1`` accanto non e' sospetto, e per
    lui il controllo resta quello di ogni file.
    """
    config, file = _single_end_con_un_nome_sospetto(tmp_path, "_R2", inverse=True)
    conclusione, codice, dettaglio = _arresto(config)
    assert (conclusione, codice) == (Conclusione.ARRESTATA, "E-S1-07")
    assert file.name in dettaglio and "0 letture su 200" in dettaglio
    assert reads.LETTURE_MINIME_PER_UN_NOME_SOSPETTO == 50
    assert reads.nome_da_seconda_lettura(file.name, [file.name])
    assert not reads.nome_da_seconda_lettura("LIB002_T02.fq.gz", [])
    assert not reads.nome_da_seconda_lettura("sito_2.fq.gz", ["sito_2.fq.gz", "altro_1.fq.gz"])
    assert reads.nome_da_seconda_lettura("sito_2.fq.gz", ["sito_2.fq.gz", "SITO_1.fq.gz"])
    assert reads.nome_da_seconda_lettura("u.R2_001.fastq.gz", [])
    # Sotto il minimo dei nomi sospetti il giudizio sulla quota non si da'.
    senza_segnale = reads.StatisticheFile(
        "f_R2.fq.gz", 49, 100, 100, 0, 0, esaurito=True, segnale_per_tratto=((0, 49),))
    assert senza_segnale.orientamento_inatteso(True, 50) is None
    assert senza_segnale.orientamento_inatteso(True, 49) is not None


@pytest.mark.parametrize(("caso", "argomenti"), [
    ("senza segnale", {"motivo": None}),
    ("poche letture", {"letture": 40}),
    ("controllo negativo", {"negativo": True}),
])
def test_il_rifiuto_per_nome_resta_dove_il_contenuto_non_si_puo_giudicare(
    tmp_path, caso, argomenti
):
    """
    **Obiettivo**: un file il cui nome finisce con ``_R2`` resta respinto per
    il nome (``E-S0-07``) quando il controllo di orientamento non lo copre:
    senza un segnale da cercare, con meno letture del minimo dei nomi sospetti,
    o se il campione e' un controllo negativo.

    **Razionale scientifico e sistemistico**: dove il contenuto non si puo'
    giudicare il nome e' l'unico segno disponibile, e la regola sbaglia per
    prudenza.
    """
    config, file = _single_end_con_un_nome_sospetto(tmp_path, "_R2", **argomenti)
    conclusione, codice, dettaglio = _arresto(config, Passo.S0)
    assert (conclusione, codice) == (Conclusione.ARRESTATA, "E-S0-07"), caso
    assert file.name in dettaglio and "sembrano letture inverse" in dettaglio


# --------------------------------------------------------------------------- #
# Soglie, grafo, stato                                                         #
# --------------------------------------------------------------------------- #


def test_la_quota_di_coppie_fuse_per_campione_e_facoltativa(tmp_path):
    """
    **Obiettivo**: Verificare che ``pair.min_merged_frac_sample_warn`` omessa
    (o nulla) non compaia nella configurazione risolta ne' fra i parametri
    assenti dichiarati, e che l'impronta dei parametri della fusione sia la
    stessa di una configurazione che la dichiara nulla; che dichiarata entri
    nella configurazione risolta e cambi quell'impronta; che
    ``campioni_sotto_la_quota`` senza soglia non giudichi nessuno, e con la
    soglia nomini i soli campioni biologici e controlli positivi con coppie la
    cui quota sta strettamente sotto; e che ``soglie`` legga come misura la
    quota piu' bassa fra quei campioni, non quella dell'insieme.

    **Razionale scientifico e sistemistico**: La quota dell'insieme non vede
    un campione solo che fonde pochissime coppie. Il controllo per campione e'
    una scelta in piu': chi non lo dichiara deve ritrovare gli stessi digest,
    le stesse impronte e gli stessi artefatti di prima.
    """
    from amplicon16s.soglie import misure
    from amplicon16s.steps.base import impronta_parametri
    from amplicon16s.steps.s04m_merge import FusioneCoppie, campioni_sotto_la_quota

    dataset = sintetico.genera_coppie(tmp_path / "dati", biologici=2, coppie=20)
    omessa = _config(dataset, tmp_path / "out")
    nulla = _config(dataset, tmp_path / "out", pair={"min_merged_frac_sample_warn": None})
    dichiarata = _config(dataset, tmp_path / "out", pair={"min_merged_frac_sample_warn": 0.5})
    assert "pair.min_merged_frac_sample_warn" in omessa.chiavi_assenti()
    assert "min_merged_frac_sample_warn" not in risolvi(omessa).come_mappa()["pair"]
    assert risolvi(omessa).digest == risolvi(nulla).digest != risolvi(dichiarata).digest
    assert risolvi(dichiarata).come_mappa()["pair"]["min_merged_frac_sample_warn"] == 0.5
    fase = FusioneCoppie()
    impronte = [impronta_parametri(risolvi(c), fase.parametri_per(c))
                for c in (omessa, nulla, dichiarata)]
    assert impronte[0] == impronte[1] != impronte[2]

    per_campione = {
        "a": {"coppie": 1000, "fuse": 950}, "b": {"coppie": 1000, "fuse": 100},
        "c": {"coppie": 0, "fuse": 0}, "n": {"coppie": 50, "fuse": 1},
        "p": {"coppie": 200, "fuse": 100},
    }
    classi = {"a": "biologico", "b": "biologico", "c": "biologico",
              "n": "controllo_negativo", "p": "controllo_positivo"}
    assert campioni_sotto_la_quota(per_campione, classi, None) == []
    assert campioni_sotto_la_quota(per_campione, classi, 0.5) == [("b", 100, 1000)]
    assert campioni_sotto_la_quota(per_campione, classi, 0.6) == [
        ("b", 100, 1000), ("p", 100, 200)]
    assert campioni_sotto_la_quota(per_campione, classi, 0.0) == []

    radice = tmp_path / "esecuzione"
    (radice / "01_input_validation").mkdir(parents=True)
    (radice / "05_pair_merge").mkdir()
    (radice / "01_input_validation" / "crosswalk.tsv").write_text(
        "accession\tsample_name\tclasse\n"
        + "".join(f"{c}\t{c}\t{classe}\n" for c, classe in classi.items()), encoding="utf-8")
    (radice / "05_pair_merge" / "fusione.tsv").write_text(
        "campione\tcoppie\tfuse\tquota\tsequenze\n"
        + "".join(f"{c}\t{v['coppie']}\t{v['fuse']}\t0\t1\n" for c, v in per_campione.items()),
        encoding="utf-8")
    misura = misure._fusione_per_campione(radice)
    assert float(misura.valore) == 0.1 and "minimo" in misura.statistica
    assert misura.origine == "05_pair_merge/fusione.tsv"
    assert CATALOGO["E-SM-04"].categoria is Categoria.DEGRADAZIONE_AUTOMATICA


def test_la_lunghezza_attesa_sta_nell_intervallo_e_la_tolleranza_non_entra_nel_filtro(tmp_path):
    """
    **Obiettivo**: Verificare che una lunghezza attesa delle sequenze fuse
    fuori dall'intervallo fra ``pair.len_min`` e ``pair.len_max`` sia respinta
    come errore di coerenza della modalita' appaiata (``E-G15-15``), e che agli
    estremi e nulla sia accettata; che in modalita' appaiata il filtro di
    lunghezza dichiari il solo intervallo, cosi' che cambiare ``asv.len_tol`` o
    il troncamento non ne cambi l'impronta dei parametri, mentre la cambia
    l'intervallo; e che in modalita' single-end la fase dipenda dalla
    tolleranza come prima.

    **Razionale scientifico e sistemistico**: Una lunghezza attesa fuori
    dall'intervallo ammesso dichiara due cose che non possono valere insieme:
    le sequenze attese uscirebbero dalla tabella. Un parametro che in una
    modalita' non ha effetto non deve far rifare una fase.
    """
    from amplicon16s.steps.base import impronta_parametri
    from amplicon16s.steps.s07_asv_length import FiltroLunghezza

    dataset = sintetico.genera_coppie(tmp_path / "dati", biologici=2, coppie=20)
    base = sintetico.configurazione_appaiata(dataset, tmp_path / "out")
    minima, massima = base["pair"]["len_min"], base["pair"]["len_max"]
    for attesa in (minima, massima, None):
        assert _config(dataset, tmp_path / "out", pair={"expected_merged_len": attesa})
    for attesa in (minima - 1, massima + 1):
        with pytest.raises(ErroreConfigurazione) as rifiuto:
            _config(dataset, tmp_path / "out", pair={"expected_merged_len": attesa})
        violazione = rifiuto_di_g15(rifiuto.value).violazioni[0]
        assert violazione.codice == "E-G15-15" and "expected_merged_len" in violazione.dettaglio

    fase = FiltroLunghezza()
    appaiata = _config(dataset, tmp_path / "out")
    assert set(fase.parametri_per(appaiata)) == {
        "asv.len_min", "asv.len_max", "io.layout", "pair.len_min", "pair.len_max"}

    def impronta(**gruppi) -> str:
        config = _config(dataset, tmp_path / "out", **gruppi)
        return impronta_parametri(risolvi(config), fase.parametri_per(config))

    assert impronta() == impronta(asv={"len_tol": 3})
    assert impronta() == impronta(filter={"truncLen": base["filter"]["truncLen"] - 1})
    assert impronta() != impronta(pair={"len_max": massima + 5})

    singolo = sintetico.genera(tmp_path / "singolo", biologici=2, negativi=0, letture=20)
    config = valida(sintetico.configurazione(singolo, tmp_path / "out"))
    con_tolleranza = valida(sintetico.configurazione(singolo, tmp_path / "out",
                                                     asv={"len_tol": 3}))
    assert fase.parametri_per(config) == fase.parametri
    assert impronta_parametri(risolvi(config), fase.parametri_per(config)) != (
        impronta_parametri(risolvi(con_tolleranza), fase.parametri_per(con_tolleranza)))


def test_le_soglie_della_modalita_appaiata_stanno_nel_meccanismo_delle_soglie(tmp_path):
    """
    **Obiettivo**: la tabella delle soglie ha le cinque della modalita'
    appaiata, con direzione, codici e margine; in modalita' single-end non se
    ne parla (:func:`soglie_di`); in modalita' appaiata i loro valori che non
    fermano sono riconosciuti come controlli spenti, e la soglia del primer
    reverse non ha effetto con il taglio o con il primer nullo.

    **Razionale scientifico e sistemistico**: una soglia nuova fuori dal
    meccanismo sarebbe una soglia che ``init`` non rende esplorativa,
    ``soglie`` non propone e nessuno dichiara spenta.
    """
    attese = {
        "pair.max_primer_hit_frac_reverse": ("massima", ("E-S0-25",), True),
        "pair.truncLen_shortfall_warn_reverse": ("massima", ("E-S1-01",), False),
        "pair.min_merged_frac_warn": ("minima", ("E-SM-02",), False),
        "pair.min_merged_frac_stop": ("minima", ("E-SM-01",), True),
        "pair.min_merged_frac_sample_warn": ("minima", ("E-SM-04",), False),
    }
    assert tabella.SOGLIE_DELLE_COPPIE == set(attese)
    for chiave, (direzione, codici, ferma) in attese.items():
        voce = tabella.soglia(chiave)
        assert (voce.direzione.value, voce.codici, voce.ferma) == (direzione, codici, ferma)
        assert all(c in CATALOGO for c in voce.codici)
    assert not {s.chiave for s in tabella.soglie_di({})} & set(attese)
    assert tabella.soglia("pair.min_merged_frac_stop").margine == tabella.MARGINE_ARRESTO_CHIMERE

    dataset = sintetico.genera_coppie(tmp_path / "dati", biologici=2, coppie=20)
    config = _config(dataset, tmp_path / "out", pair={
        "primer_sequence_reverse": sintetico.PRIMER_805R, "max_primer_hit_frac_reverse": 1.0,
        "min_merged_frac_warn": 0.0, "min_merged_frac_stop": 0.0,
        "min_merged_frac_sample_warn": 0.0,
        "truncLen_shortfall_warn_reverse": dataset.lunghezza - 1})
    valori = tabella.valori_piatti(config)
    assert set(attese) <= {s.chiave for s in tabella.soglie_di(valori)}
    spenti = {s.chiave: s for s in tabella.controlli_spenti(valori, 2)}
    assert set(attese) <= set(spenti)
    # La soglia per campione, omessa, non e' un controllo spento: non esiste.
    senza = tabella.valori_piatti(_config(dataset, tmp_path / "out"))
    assert "pair.min_merged_frac_sample_warn" not in senza
    assert "pair.min_merged_frac_sample_warn" not in {
        s.chiave for s in tabella.controlli_spenti(senza, 2)}
    assert "pair.truncLen_reverse meno 1" in spenti["pair.truncLen_shortfall_warn_reverse"].testo
    assert tabella.valore_esplorativo(tabella.soglia("pair.min_merged_frac_stop")) == 0.0
    assert tabella.valore_esplorativo(tabella.soglia("pair.max_primer_hit_frac_reverse")) == 1.0
    voce = tabella.soglia("pair.max_primer_hit_frac_reverse")
    assert "non cerca il primer" in tabella.non_eseguito(
        voce, {**valori, "pair.trimLeft_reverse": 20})
    assert "nullo" in tabella.non_eseguito(voce, {**valori, "pair.primer_sequence_reverse": None})


def test_la_fusione_e_una_fase_della_sola_modalita_appaiata(tmp_path):
    """
    **Obiettivo**: in modalita' single-end la valutazione dello stato non
    conosce la fase della fusione, S5 non ne dipende e i parametri dichiarati
    da ogni fase sono i soli dichiarati dalla fase; in modalita' appaiata la fusione sta
    fra S4 e S5, S5 ne dipende, e ogni fase dichiara anche la modalita'. Dopo
    la validazione ``status`` elenca la fusione solo in modalita' appaiata.

    **Razionale scientifico e sistemistico**: una fase che esiste solo
    nell'altra modalita' non e' una fase disattivata: elencarla cambierebbe lo
    stato, il report e i manifesti di ogni esecuzione single-end.
    """
    singolo = sintetico.genera(tmp_path / "singolo", biologici=2, letture=30)
    config_singola = valida(sintetico.configurazione(singolo, tmp_path / "out_singolo"))
    coppie = sintetico.genera_coppie(tmp_path / "coppie", biologici=2, coppie=40)
    config_appaiata = _config(coppie, tmp_path / "out_coppie")

    run = ProjectRun(config_singola)
    assert Passo.SM not in run.valuta().situazioni
    assert GRAFO.dipendenze_attive(Passo.S5, config_singola) == (Passo.S0, Passo.S2, Passo.S4)
    # Il filtro appaiato legge da S1 la verifica della sincronia delle coppie.
    assert GRAFO.dipendenze_attive(Passo.S2, config_singola) == (Passo.S0,)
    assert GRAFO.dipendenze_attive(Passo.S2, config_appaiata) == (Passo.S0, Passo.S1)
    assert Passo.S1 in GRAFO.antenati(Passo.S6, config_appaiata)
    assert Passo.S1 not in GRAFO.antenati(Passo.S6, config_singola)
    assert Passo.SM not in GRAFO.elencate(config_singola)
    for fase in run.passi.values():
        assert fase.parametri_per(config_singola) == fase.parametri

    run = ProjectRun(config_appaiata)
    ordine = list(run.valuta().situazioni)
    assert ordine.index(Passo.SM) == ordine.index(Passo.S4) + 1 == ordine.index(Passo.S5) - 1
    assert Passo.SM in GRAFO.dipendenze_attive(Passo.S5, config_appaiata)
    assert GRAFO.nodo(Passo.SM).cartella is Fase.PAIR_MERGE
    for passo, fase in run.passi.items():
        assert "io.layout" in fase.parametri_per(config_appaiata), passo
    # I parametri della fusione li dichiara la sola fusione: cambiarli non
    # rifa' il filtro ne' l'inferenza.
    for parametro in ("pair.min_overlap", "pair.max_mismatch", "pair.min_merged_frac_stop"):
        assert [p for p, f in run.passi.items()
                if parametro in f.parametri_per(config_appaiata)] == [Passo.SM], parametro

    for config, attesa in ((config_singola, False), (config_appaiata, True)):
        assert _arresto(config, Passo.S0)[0] is Conclusione.COMPLETATA
        fasi = {f["passo"]: f["stato"] for f in interfaccia.stato(config.io.out_root)["fasi"]}
        assert ("SM" in fasi) is attesa
        assert fasi["S0"] == "conclusa" and fasi["S9"] == "disattivata"
        if attesa:
            assert fasi["SM"] == fasi["S5"] == "da_eseguire"


def test_la_quota_di_coppie_fuse_si_giudica_sulle_classi_controllate(tmp_path):
    """
    **Obiettivo**: :func:`misura_fusione` somma coppie e coppie fuse per
    classe e calcola la quota sui soli campioni biologici e controlli positivi;
    :func:`controlla_fusione` solleva ``E-SM-01`` sotto
    ``pair.min_merged_frac_stop``, restituisce il testo dell'avviso
    (``E-SM-02``) fra le due soglie e nulla sopra; una quota uguale alla soglia
    non fa scattare il controllo, e senza coppie da giudicare non si giudica.
    :func:`mediana_pesata` pesa le lunghezze sulle coppie.

    **Razionale scientifico e sistemistico**: nei controlli negativi quasi
    nessuna coppia si fonde, perche' quasi nulla e' amplificato: contarli
    abbasserebbe la quota di un dataset in ordine. Il confronto stretto e' cio'
    che rende 0 il valore che non ferma.
    """
    from amplicon16s.errors.exceptions import ErrorePipeline
    from amplicon16s.steps import s04m_merge as fusione

    per_campione = {"a": {"coppie": 1000, "fuse": 800}, "b": {"coppie": 1000, "fuse": 600},
                    "n": {"coppie": 500, "fuse": 0}}
    classi = {"a": "biologico", "b": "controllo_positivo", "n": "controllo_negativo"}
    misure = fusione.misura_fusione(per_campione, classi)
    assert (misure["coppie_controllate"], misure["fuse_controllate"]) == (2000, 1400)
    assert misure["quota_controllata"] == 0.7
    assert misure["per_classe"]["controllo_negativo"] == {
        "coppie": 500, "fuse": 0, "campioni": 1, "quota": 0.0}

    dataset = sintetico.genera_coppie(tmp_path / "dati", biologici=2, coppie=20)

    def soglie(avviso: float, arresto: float):
        return _config(dataset, tmp_path / "out", pair={
            "min_merged_frac_warn": avviso, "min_merged_frac_stop": arresto}).pair

    with pytest.raises(ErrorePipeline) as arresto:
        fusione.controlla_fusione(misure, soglie(0.9, 0.8))
    assert arresto.value.codice == "E-SM-01"
    assert "il 70.0% delle coppie" in str(arresto.value) and "1400 su 2000" in str(arresto.value)
    avviso = fusione.controlla_fusione(misure, soglie(0.9, 0.5))
    assert avviso is not None and "pair.min_merged_frac_warn" in avviso
    assert fusione.controlla_fusione(misure, soglie(0.7, 0.7)) is None
    assert fusione.controlla_fusione(misure, soglie(0.0, 0.0)) is None
    senza = fusione.misura_fusione({"n": per_campione["n"]}, classi)
    assert senza["quota_controllata"] is None
    assert fusione.controlla_fusione(senza, soglie(0.9, 0.8)) is None
    assert (CATALOGO["E-SM-01"].categoria, CATALOGO["E-SM-02"].categoria) == (
        Categoria.REVISIONE_UMANA, Categoria.DEGRADAZIONE_AUTOMATICA)

    assert fusione.mediana_pesata({250: 10, 253: 80, 400: 10}) == 253
    assert fusione.mediana_pesata({250: 50, 253: 50}) == 250
    assert fusione.mediana_pesata({}) is None


# --------------------------------------------------------------------------- #
# init                                                                         #
# --------------------------------------------------------------------------- #


def test_la_lunghezza_dell_inserto_si_misura_sulla_sovrapposizione():
    """
    **Obiettivo**: :func:`lunghezza_dell_inserto` restituisce la lunghezza del
    frammento di cui due letture sono le estremita', anche con qualche errore
    nella sovrapposizione e quando il frammento e' piu' corto delle letture;
    restituisce ``None`` per due letture che non si sovrappongono. La stima su
    piu' coppie non si da' se meno di meta' si sovrappone.

    **Razionale scientifico e sistemistico**: la lunghezza della regione si
    misura sulle letture, non si prende dai primer dichiarati nei metadati, che
    possono essere quelli di un'altra regione.
    """
    caso = random.Random(3)
    frammento = "".join(caso.choice("ACGT") for _ in range(253))
    forward = frammento[:150]
    reverse = rilevamento.complemento_inverso(frammento)[:150]
    assert rilevamento.lunghezza_dell_inserto(forward, reverse) == 253
    con_errori = list(reverse)
    for posto in (5, 40):
        con_errori[posto] = "A" if con_errori[posto] != "A" else "C"
    assert rilevamento.lunghezza_dell_inserto(forward, "".join(con_errori)) == 253
    corto = frammento[:120]
    assert rilevamento.lunghezza_dell_inserto(
        corto, rilevamento.complemento_inverso(corto)) == 120
    estranea = "".join(caso.choice("ACGT") for _ in range(150))
    assert rilevamento.lunghezza_dell_inserto(forward, estranea) is None
    assert rilevamento.stima_dell_inserto([[253, 253, None], [None, None, None]]) is None
    stima = rilevamento.stima_dell_inserto([[253, 252, 254, None], [253, 253, 253]])
    assert stima is not None and (stima.mediana, stima.coppie, stima.sovrapposte) == (253, 7, 6)
    assert rilevamento.stima_dell_inserto([]) is None


def test_i_primer_reverse_dell_elenco_hanno_la_fonte_e_la_posizione():
    """
    **Obiettivo**: ogni primer dell'elenco dei reverse, e dei forward, porta la
    pubblicazione che lo definisce e la posizione della sua estremita' 5' sul
    16S di *E. coli*; e per le coppie di primer piu' usate la regione fra i due
    ha la lunghezza nota (253 basi senza i primer per 515F e 806R).

    **Razionale scientifico e sistemistico**: le posizioni servono solo a
    confrontare la lunghezza misurata sulle letture con quella attesa dai
    primer riconosciuti: sono conoscenza di metodo, e devono avere una fonte.
    """
    per_nome = {p.nome: p for p in (*vocabolario.PRIMER_16S, *vocabolario.PRIMER_16S_REVERSE)}
    assert all(p.fonte and p.posizione for p in per_nome.values())
    assert {"806R", "806R (originale)", "805R", "926R"} <= set(per_nome)
    forward, reverse = per_nome["515F"], per_nome["806R"]
    con_i_primer = reverse.posizione - forward.posizione + 1
    assert con_i_primer - len(forward.sequenza) - len(reverse.sequenza) == 253


def test_init_riconosce_un_dataset_appaiato_e_ne_propone_i_parametri(
        tmp_path, capsys, monkeypatch):
    """
    **Obiettivo**: su una cartella con le due letture di ogni frammento, che
    portano i due primer, ``init`` rileva la modalita' appaiata e i marcatori
    dei file, rileva il primer reverse dell'elenco in testa alle letture
    reverse e il taglio pari alla sua lunghezza, propone il motivo con cui le
    letture reverse proseguono, e propone lunghezza attesa e intervallo delle
    sequenze fuse dalla lunghezza dell'inserto misurata, tolti i due tagli. La
    lunghezza attesa dai primer riconosciuti (quella della regione vera, piu'
    lunga del frammento sintetico) non concorda: ``init`` lo dice, e vale la
    misurata. Un intervallo misurato piu' stretto della variabilita' di una
    regione (qui ogni frammento ha la stessa lunghezza) non e' proposto: i due
    estremi restano da scegliere con un intervallo raccomandato; uno largo
    quanto una regione vera lo e'. Con ``--accetta-proposte`` le soglie della
    modalita' appaiata sono esplorative.

    **Razionale scientifico e sistemistico**: nessun valore dei dati entra
    nella configurazione scelto a mano: cio' che la fusione richiede si misura
    sulle letture, e dove una fonte indiretta lo contraddice si dichiara.
    """
    dataset = sintetico.genera_coppie(
        tmp_path / "dati", biologici=3, coppie=300, lunghezza=120,
        primer=sintetico.PRIMER_341F, primer_reverse=sintetico.PRIMER_805R)
    uscita = tmp_path / "c.yaml"
    assert main(["init", "--reads", str(dataset.letture), "--metadata", str(dataset.tabella),
                 "--output", str(uscita), "--skip-profile", "--accetta-proposte"]) == 0
    stampato = capsys.readouterr().out
    righe = righe_dei_parametri(uscita.read_text(encoding="utf-8"))
    assert righe["io.layout"] == (True, '"paired-end"')
    assert righe["pair.file_markers"] == (True, '["_R1", "_R2"]')
    assert righe["io.accession_regex"][0] is True
    assert righe["qc.primer_sequence"] == (True, f'"{sintetico.PRIMER_341F}"')
    assert righe["pair.primer_sequence_reverse"] == (True, f'"{sintetico.PRIMER_805R}"')
    assert righe["pair.trimLeft_reverse"] == (True, str(len(sintetico.PRIMER_805R)))
    attivo, motivo = righe["pair.conserved_motif_reverse"]
    assert attivo and yaml.safe_load(motivo).startswith(sintetico.MOTIVO_REVERSE[:8])
    assert righe["pair.expected_merged_len"] == (True, str(sintetico.INSERTO))
    # Ogni frammento sintetico ha la stessa lunghezza: un intervallo largo zero
    # non descrive la regione, e i due estremi restano da scegliere, con
    # l'intervallo raccomandato nel commento e nel riepilogo.
    testo = " ".join(uscita.read_text(encoding="utf-8").replace("#", " ").split())
    for chiave in ("pair.len_min", "pair.len_max"):
        assert righe[chiave] == (False, ""), chiave
    largo = (f"da {math.floor(sintetico.INSERTO * 0.95)} a "
             f"{math.ceil(sintetico.INSERTO * 1.05)}")
    assert f"e' largo meno del 2% della mediana" in testo and largo in testo
    assert "Intervallo delle lunghezze fuse non proposto" in " ".join(stampato.split())
    assert "non concorda con quella attesa dai primer riconosciuti" in stampato
    assert "vale la misurata" in stampato
    for chiave in ("pair.max_primer_hit_frac_reverse", "pair.min_merged_frac_warn",
                   "pair.min_merged_frac_stop"):
        assert righe[chiave][0] is True, chiave
    assert righe["pair.min_merged_frac_stop"][1] == "0.0"
    # Senza il profilo di qualita' i troncamenti restano da scegliere, con il vincolo.
    assert righe["filter.truncLen"] == (False, "") and righe["pair.truncLen_reverse"] == (False, "")
    assert "vincolo della fusione" in uscita.read_text(encoding="utf-8")
    # I predefiniti del gruppo restano in commento con il loro valore.
    assert righe["pair.min_overlap"] == (False, str(defaults.PAIR_MIN_OVERLAP))

    # Lunghezze che variano come in una regione vera: l'intervallo e' proposto.
    from amplicon16s.init import rilevamento

    vera = rilevamento.stima_dell_inserto

    def con_variazione(per_file, *altri):
        return vera([[None if x is None else x + (n % 9) - 4 for n, x in enumerate(f)]
                     for f in per_file], *altri)

    monkeypatch.setattr(rilevamento, "stima_dell_inserto", con_variazione)
    larga = tmp_path / "larga.yaml"
    assert main(["init", "--reads", str(dataset.letture), "--metadata", str(dataset.tabella),
                 "--output", str(larga), "--skip-profile", "--accetta-proposte"]) == 0
    righe = righe_dei_parametri(larga.read_text(encoding="utf-8"))
    assert righe["pair.len_min"] == (True, str(sintetico.INSERTO - 4))
    assert righe["pair.len_max"] == (True, str(sintetico.INSERTO + 4))
    assert "Intervallo delle lunghezze fuse non proposto" not in capsys.readouterr().out


def _profilo_finto(
        frequente: int, ammessa: int, per_qualita: int, per_fondere: int | None = None,
        corte: int = 0) -> dict:
    """Il profilo che ``init`` ricava dallo script di S1, costruito a mano: 1.000
    letture di ``frequente`` basi (``corte`` di 20), con la lunghezza ammessa
    data, qualita' mediana 36 fino a ``per_qualita``, 24 fino a ``per_fondere``
    e 12 oltre.
    """
    per_fondere = per_qualita if per_fondere is None else per_fondere
    distribuzione = {frequente: 1000 - corte, **({20: corte} if corte else {})}
    mediana = {p: 36.0 if p <= per_qualita else 24.0 if p <= per_fondere else 12.0
               for p in range(1, frequente + 1)}
    return {
        "suggerito": {
            "senza_perdite": min(distribuzione), "per_lunghezza": ammessa,
            "per_qualita": per_qualita, "suggerito": min(ammessa, per_qualita),
            "perdita_peggiore": 0.0, "campione_peggiore": None, "perdita_per_classe": {}},
        "mediana": mediana, "per_fondere": per_fondere,
        "valori_di_qualita": list(range(2, 41)),
        "minima": min(distribuzione), "massima": frequente, "moda": frequente,
        "distribuzione": distribuzione, "totale": 1000,
        "qualita": {p: (q, q - 2) for p, q in mediana.items()},
    }


def test_init_non_allunga_i_troncamenti_dentro_basi_di_qualita_bassa(
        tmp_path, capsys, monkeypatch):
    """
    **Obiettivo**: con letture reverse accorciate prima del deposito e di
    qualita' bassa quasi dall'inizio (la gran parte lunga 100 basi, qualita'
    mediana che regge a 30 fino alla posizione 15 e a 20 fino alla 25; forward
    che regge a 30 fino alla 60 e a 20 fino alla 70, su 120),
    ``init --accetta-proposte`` non propone i due troncamenti: restano da
    scegliere, senza valore, e il riepilogo dice che le coppie non si fondono
    senza basi di qualita' bassa, riporta per le due letture la lunghezza piu'
    frequente e il limite di qualita', e consiglia le sole forward. Quando per
    fondere basta allungare una lettura dentro posizioni in cui la qualita'
    mediana e' almeno 20, i due troncamenti sono proposti e il commento dice
    di quante basi e fino a quale qualita'; una base oltre, no; e quando la
    qualita' regge senza allungare la proposta e' quella di sempre.

    **Razionale scientifico e sistemistico**: il vincolo della fusione si puo'
    sempre rispettare allungando i troncamenti, ma una sovrapposizione fatta di
    basi di qualita' bassa fa scartare le letture al filtro e fonde le coppie
    su basi inaffidabili. In quel caso la scelta giusta e' un'altra (le sole
    forward), e una proposta attiva la nasconderebbe.
    """
    from amplicon16s.init import comando

    dataset = sintetico.genera_coppie(
        tmp_path / "dati", biologici=3, coppie=300, lunghezza=120,
        primer=sintetico.PRIMER_341F, primer_reverse=sintetico.PRIMER_805R)

    def lancia(nome: str, forward: dict, reverse: dict) -> tuple[dict, str, str]:
        profili = iter((forward, reverse))
        monkeypatch.setattr(comando, "_profilo", lambda *a, **k: next(profili))
        uscita = tmp_path / nome
        assert main(["init", "--reads", str(dataset.letture), "--metadata", str(dataset.tabella),
                     "--output", str(uscita), "--accetta-proposte"]) == 0
        testo = uscita.read_text(encoding="utf-8")
        return (righe_dei_parametri(testo), " ".join(capsys.readouterr().out.split()),
                " ".join(testo.replace("#", " ").split()))

    assert comando.QUALITA_MINIMA_PER_FONDERE == 20
    # Reverse corte e scadenti: nessuna proposta, e la nota che lo spiega.
    righe, stampato, testo = lancia(
        "scadenti.yaml", _profilo_finto(120, 120, 60, 70),
        _profilo_finto(100, 100, 15, 25, corte=50))
    assert righe["filter.truncLen"] == (False, "")
    assert righe["pair.truncLen_reverse"] == (False, "")
    assert "Le coppie non si fondono senza basi di qualita' bassa" in stampato
    assert "fino alla posizione 60 nelle letture forward (lunghezza piu' frequente 120)" in stampato
    assert "fino alla 15 nelle reverse (lunghezza piu' frequente 100)" in stampato
    assert "fin dove la qualita' mediana e' almeno 20 le due letture danno 70 + 25 = 95 basi" \
        in stampato
    assert "Si consiglia di usare le sole letture forward: prepara una cartella con i soli " \
           "file forward (GUIDA.md, passo 2) e rilancia init con --reads su quella cartella " \
           "e un altro --output." in stampato
    assert "le coppie non si fondono senza basi di qualita' bassa" in testo
    necessarie = int(stampato.split("e ne servono ")[1].split()[0])
    # La somma delle lunghezze ammesse basterebbe: e' la qualita' a fermare.
    assert 120 + 100 >= necessarie > 95

    # Dentro basi di qualita' almeno 20: proposti, la forward allungata di quanto serve.
    reverse_buona = necessarie - 100
    righe, stampato, testo = lancia(
        "entro.yaml", _profilo_finto(120, 120, 90, 100),
        _profilo_finto(120, 120, reverse_buona))
    assert righe["filter.truncLen"] == (True, "100")
    assert righe["pair.truncLen_reverse"] == (True, str(reverse_buona))
    assert ("allungato di 10 basi oltre il calo della qualita' quanto serve a fondere le coppie, "
            "dentro posizioni in cui la qualita' mediana e' almeno 20 (la piu' bassa: 24)") in testo
    assert "non si fondono" not in stampato
    # Una base in meno di qualita' almeno 20: di nuovo da scegliere.
    righe, stampato, _ = lancia(
        "oltre.yaml", _profilo_finto(120, 120, 90, 99), _profilo_finto(120, 120, reverse_buona))
    assert righe["filter.truncLen"] == (False, "") and "non si fondono" in stampato

    # La qualita' regge: proposti dove regge, senza allungare, come prima.
    righe, stampato, testo = lancia(
        "buone.yaml", _profilo_finto(120, 120, 115), _profilo_finto(120, 120, 108))
    assert righe["filter.truncLen"] == (True, "115")
    assert righe["pair.truncLen_reverse"] == (True, "108")
    assert "allungato" not in testo and "non si fondono" not in stampato


def test_il_taglio_delle_reverse_resta_da_scegliere_se_un_file_non_porta_il_primer(
    tmp_path, capsys
):
    """
    **Obiettivo**: Verificare che in modalita' appaiata, quando il primer
    reverse e' in testa alle letture reverse di quasi tutti i file ma non di
    uno, ``init`` proponga il primer reverse e lasci da scegliere, senza
    valore, il taglio delle reverse, nominando il file che non lo porta e il
    taglio che sarebbe; che il riepilogo lo dica; che il taglio delle forward,
    che ogni file porta, resti un fatto rilevato; e che con il primer in ogni
    file reverse anche quel taglio sia rilevato.

    **Razionale scientifico e sistemistico**: Un taglio fisso toglie le prime
    basi di ogni lettura: al file che non porta il primer toglierebbe basi
    vere della regione, spostando la sovrapposizione delle coppie di quel solo
    campione. Vale per le reverse come per le forward, con la stessa regola.
    """
    def lancia(dataset: sintetico.DatasetAppaiato, uscita: Path) -> dict:
        assert main(["init", "--reads", str(dataset.letture), "--metadata",
                     str(dataset.tabella), "--output", str(uscita), "--skip-profile",
                     "--threads", "1"]) == 0
        return righe_dei_parametri(uscita.read_text(encoding="utf-8"))

    forward, reverse = sintetico.PRIMER_341F, vocabolario.PRIMER_16S_REVERSE[0].sequenza
    intero = sintetico.genera_coppie(tmp_path / "intero", biologici=4, coppie=80,
                                     primer=forward, primer_reverse=reverse)
    righe = lancia(intero, tmp_path / "intero.yaml")
    assert righe["pair.trimLeft_reverse"] == (True, str(len(reverse)))
    capsys.readouterr()

    misto = sintetico.genera_coppie(tmp_path / "misto", biologici=4, coppie=80,
                                    primer=forward, primer_reverse=reverse)
    senza = sorted(misto.coppie)[0]
    del_reverse = misto.coppie[senza][1]
    sintetico._scrivi_fastq(del_reverse, [
        (nome, sequenza[len(reverse):], qualita[len(reverse):])
        for nome, sequenza, qualita in _record(del_reverse)])
    righe = lancia(misto, tmp_path / "misto.yaml")
    console = " ".join(capsys.readouterr().out.split())
    assert righe["pair.trimLeft_reverse"] == (False, "")
    assert righe["pair.primer_sequence_reverse"][0] is False
    assert yaml.safe_load(righe["pair.primer_sequence_reverse"][1]) == reverse
    assert righe["filter.trimLeft"] == (True, str(len(forward)))
    piatto = " ".join((tmp_path / "misto.yaml").read_text(encoding="utf-8").replace(
        "#", " ").split())
    assert "1 file su 4 non lo portano in testa" in piatto
    assert f"il taglio sarebbe {len(reverse)}" in piatto
    assert "maggior parte dei file ma non di 1: il primer e' proposto, il taglio resta da " \
           "scegliere" in console
