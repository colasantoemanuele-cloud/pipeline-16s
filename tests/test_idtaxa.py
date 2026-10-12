"""Il classificatore IdTaxa (``tax.classifier: idtaxa``) accanto al bayesiano ingenuo.

Verifica che il classificatore di chi non ne dichiara uno resti rappresentato
com'era (stesso digest, stessa configurazione registrata), che i parametri di
un classificatore siano respinti con l'altro, che il classificatore addestrato
sia verificato prima delle fasi di calcolo (checksum e contenuto), che i ranghi
si riconoscano per nome (``Root`` scartato, un rango in piu' scartato e
dichiarato, un rango saltato lasciato vuoto con i ranghi assegnati sotto), che
il rango saltato arrivi all'oggetto finale e alle analisi ecologiche, e che due
esecuzioni diano gli stessi byte.

Il classificatore di prova e' addestrato sul riferimento sintetico
(``tests/fixtures/idtaxa/costruisci.R``): i nomi dei taxa sono inventati.

I test che classificano richiedono R con DECIPHER: senza si saltano, e
nell'immagine dell'integrazione continua devono girare.

    docker run --rm -u "$(id -u):$(id -g)" -e HOME=/tmp -e PYTHONPATH=/app/src \\
      -e AMPLICON16S_R_DIR=/app/R -v "$(pwd)":/app -w /app <immagine> \\
      pytest tests/test_idtaxa.py
"""

from __future__ import annotations

import copy
import csv
import functools
import hashlib
import json
import os
import subprocess
from pathlib import Path

import pytest
import sintetico
import yaml
from eco_aiuti import configurazione as configurazione_eco
from eco_aiuti import leggi_tsv
from sottoinsieme import motivo_pacchetti_r_assenti

from amplicon16s import cli
from amplicon16s.config import defaults
from amplicon16s.config.resolve import parametri_dichiarati, risolvi, scrivi_risolta
from amplicon16s.config.schema import ErroreConfigurazione, valida
from amplicon16s.errors.catalog import CATALOGO, Categoria
from amplicon16s.gates.g01_g15 import Contesto, rifiuto_di_g15
from amplicon16s.gates.registry import esegui_gate
from amplicon16s.logging.logger import chiudi
from amplicon16s.rbridge.runner import trova_rscript
from amplicon16s.runner.executor import Conclusione, Esecutore
from amplicon16s.runner.graph import Passo
from amplicon16s.runner.project import ProjectRun
from amplicon16s_eco.cli import main as main_eco

RADICE = Path(__file__).resolve().parents[1]
#: Il classificatore addestrato sul riferimento sintetico dei dataset generati.
ADDESTRATO = RADICE / "tests" / "fixtures" / "idtaxa" / "addestrato_sintetico.RData"
PRIMER = "GTGYCAGCMGCCGCGGTAA"
COLONNE = ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus"]

pytestmark = pytest.mark.xdist_group("catena_idtaxa")


@pytest.fixture(autouse=True)
def uscite_pulite():
    """Chiude le uscite del log prima e dopo ogni test."""
    chiudi()
    yield
    chiudi()


@functools.cache
def _sonda() -> str | None:
    """Perche' l'ambiente R richiesto non c'e', o ``None``."""
    return motivo_pacchetti_r_assenti(
        "dada2", "ShortRead", "jsonlite", "phyloseq", "Biostrings", "DECIPHER", "decontam")


def _richiedi_r() -> None:
    if _sonda() is not None:
        if os.environ.get("AMPLICON16S_RICHIEDI_BIOC") == "1":
            pytest.fail(f"Bioconductor e' richiesto in questo ambiente: {_sonda()}")
        pytest.skip(_sonda())


@pytest.fixture
def con_r():
    """Richiede R con DECIPHER: salta senza, ma in CI fallisce."""
    _richiedi_r()


def _con_idtaxa(dati: dict, addestrato: Path = ADDESTRATO, **tax) -> dict:
    """La configurazione con IdTaxa al posto del bayesiano ingenuo."""
    dati = copy.deepcopy(dati)
    for chiave in ("ref_fasta", "ref_md5", "try_rc", "min_boot", "ref_bad_taxa"):
        dati["tax"].pop(chiave, None)
    dati["tax"].update(
        classifier="idtaxa", idtaxa_trained=str(addestrato),
        idtaxa_trained_md5=hashlib.md5(Path(addestrato).read_bytes()).hexdigest(),
        idtaxa_threshold=60, idtaxa_strand="top")
    dati["tax"].update(tax)
    return dati


def _tabella(percorso: Path) -> list[dict[str, str]]:
    with open(percorso, encoding="utf-8", newline="") as file:
        return list(csv.DictReader(file, delimiter="\t"))


def _impronte(uscita: Path) -> dict[str, str]:
    """Lo sha256 dichiarato di ogni artefatto di ogni fase."""
    return {
        f"{m.parent.name}/{a['nome']}": a["checksum"]
        for m in uscita.rglob("manifest_S*.json")
        for a in json.loads(m.read_text(encoding="utf-8"))["artefatti"]}


# --------------------------------------------------------------------------- #
# La configurazione                                                            #
# --------------------------------------------------------------------------- #


def test_il_bayesiano_ingenuo_dichiarato_vale_come_omesso(tmp_path):
    """
    **Obiettivo**: Verificare che una configurazione senza ``tax.classifier``
    e una che dichiara ``naive_bayes`` abbiano lo stesso digest, la stessa
    impronta dei risultati e la stessa configurazione registrata, con
    ``classifier: naive_bayes`` e senza alcun parametro di IdTaxa; e che con
    ``idtaxa`` la configurazione registrata porti i quattro parametri propri e
    nessuno di quelli del bayesiano ingenuo.

    **Razionale scientifico e sistemistico**: Il secondo classificatore non
    deve cambiare alcun artefatto di chi usa il primo: stesso digest, stesse
    impronte di fase, nessuna fase da rifare. E una configurazione registrata
    che riportasse il riferimento FASTA con IdTaxa direbbe usato un file che
    non lo e'.
    """
    dataset = sintetico.genera(tmp_path / "dati", biologici=2, letture=30)
    omesso = sintetico.configurazione(dataset, tmp_path / "out")
    assert "classifier" not in omesso["tax"]
    dichiarato = copy.deepcopy(omesso)
    dichiarato["tax"]["classifier"] = "naive_bayes"
    a, b = risolvi(valida(omesso)), risolvi(valida(dichiarato))
    assert not a.config.idtaxa and a.digest == b.digest
    assert a.impronta_risultati == b.impronta_risultati
    testi = [scrivi_risolta(r, tmp_path / nome).read_text(encoding="utf-8")
             for nome, r in (("a", a), ("b", b))]
    registrata = yaml.safe_load(testi[0])["parametri"]["tax"]
    assert registrata["classifier"] == "naive_bayes" and "idtaxa" not in testi[0]
    assert {"ref_fasta", "ref_md5", "min_boot", "try_rc", "ref_bad_taxa"} <= set(registrata)
    # La sola differenza fra i due file: il parametro e' fra i dichiarati.
    assert yaml.safe_load(testi[0])["parametri"] == yaml.safe_load(testi[1])["parametri"]
    assert set(parametri_dichiarati(b.config)) - set(parametri_dichiarati(a.config)) == {
        "tax.classifier"}

    idtaxa = risolvi(valida(_con_idtaxa(omesso)))
    assert idtaxa.config.idtaxa and idtaxa.digest != a.digest
    del_classificatore = idtaxa.come_mappa()["tax"]
    assert set(del_classificatore) == {
        "ref_name", "ref_version", "classifier", "assign_species", "idtaxa_trained",
        "idtaxa_trained_md5", "idtaxa_threshold", "idtaxa_strand"}
    assert (del_classificatore["idtaxa_threshold"], del_classificatore["idtaxa_strand"]) == (
        60.0, "top")


def test_i_parametri_di_un_classificatore_sono_respinti_con_l_altro(tmp_path):
    """
    **Obiettivo**: Verificare che con ``tax.classifier: idtaxa`` ogni
    parametro del bayesiano ingenuo dichiarato (riferimento FASTA, checksum,
    ``try_rc``, ``min_boot``, elenco dei taxa difettosi) sia respinto con
    ``E-G15-17``, che lo nomina; che con il bayesiano ingenuo lo sia ogni
    parametro di IdTaxa; che i parametri mancanti del classificatore scelto
    siano elencati da ``E-G15-10`` con la condizione; e che una soglia fuori da
    0-100 o un filamento sconosciuto siano respinti.

    **Razionale scientifico e sistemistico**: Un parametro accettato senza
    effetto farebbe credere a un bootstrap minimo applicato a IdTaxa, o a un
    classificatore addestrato usato da ``assignTaxonomy``: la tassonomia
    riportata non sarebbe quella dichiarata.
    """
    dataset = sintetico.genera(tmp_path / "dati", biologici=2, letture=30)
    base = sintetico.configurazione(dataset, tmp_path / "out")

    def rifiuto(dati: dict) -> list[tuple[str, str]]:
        with pytest.raises(ErroreConfigurazione) as respinta:
            valida(dati)
        return [(v.codice, v.dettaglio) for v in rifiuto_di_g15(respinta.value).violazioni]

    for chiave, valore in (
        ("ref_fasta", base["tax"]["ref_fasta"]), ("ref_md5", base["tax"]["ref_md5"]),
        ("try_rc", False), ("min_boot", 50), ("ref_bad_taxa", None),
    ):
        ((codice, dettaglio),) = rifiuto(_con_idtaxa(base, **{chiave: valore}))
        assert codice == "E-G15-17" and f"tax.{chiave}" in dettaglio and "idtaxa" in dettaglio
    for chiave, valore in (
        ("idtaxa_trained", str(ADDESTRATO)), ("idtaxa_trained_md5", "0" * 32),
        ("idtaxa_threshold", 60), ("idtaxa_strand", "both"),
    ):
        con = copy.deepcopy(base)
        con["tax"][chiave] = valore
        ((codice, dettaglio),) = rifiuto(con)
        assert codice == "E-G15-17" and f"tax.{chiave}" in dettaglio
        assert "tax.classifier" in dettaglio

    manca = _con_idtaxa(base)
    del manca["tax"]["idtaxa_threshold"], manca["tax"]["idtaxa_trained"]
    ((codice, dettaglio),) = rifiuto(manca)
    assert codice == "E-G15-10"
    assert "tax.idtaxa_trained" in dettaglio and "tax.idtaxa_threshold" in dettaglio
    assert "perche' tax.classifier e' idtaxa" in dettaglio
    # Il nullo scritto a mano vale come l'omissione.
    ((codice, dettaglio),) = rifiuto(_con_idtaxa(base, idtaxa_strand=None))
    assert codice == "E-G15-10" and "tax.idtaxa_strand" in dettaglio
    for chiave, valore in (("idtaxa_threshold", 101), ("idtaxa_threshold", -1),
                           ("idtaxa_threshold", True), ("idtaxa_strand", "forward")):
        with pytest.raises(ErroreConfigurazione, match=f"tax.{chiave}"):
            valida(_con_idtaxa(base, **{chiave: valore}))
    with pytest.raises(ErroreConfigurazione, match="tax.classifier"):
        valida(_con_idtaxa(base, classifier="rdp"))
    # Con il bayesiano ingenuo il bootstrap minimo non puo' essere nullo.
    senza_soglia = copy.deepcopy(base)
    senza_soglia["tax"]["min_boot"] = None
    with pytest.raises(ErroreConfigurazione, match="tax.min_boot"):
        valida(senza_soglia)
    assert CATALOGO["E-G15-17"].categoria is Categoria.REVISIONE_UMANA


def test_il_checksum_del_classificatore_addestrato_si_verifica_a_ogni_avvio(tmp_path):
    """
    **Obiettivo**: Verificare che con IdTaxa il gate del riferimento confronti
    l'MD5 del file di ``tax.idtaxa_trained`` con quello dichiarato e si fermi
    con ``E-S0-12`` se differisce o se il file manca, nominando il parametro; e
    che non applichi al classificatore addestrato i controlli del riferimento
    FASTA.

    **Razionale scientifico e sistemistico**: Il classificatore addestrato
    decide ogni assegnazione: un file sostituito dopo la validazione darebbe
    una tassonomia diversa sotto lo stesso nome e la stessa versione
    dichiarati.
    """
    dataset = sintetico.genera(tmp_path / "dati", biologici=2, letture=30)
    base = sintetico.configurazione(dataset, tmp_path / "out")
    assert esegui_gate("G12", Contesto(valida(_con_idtaxa(base)))).superato

    diverso = esegui_gate("G12", Contesto(valida(_con_idtaxa(base, idtaxa_trained_md5="0" * 32))))
    (violazione,) = diverso.violazioni
    assert violazione.codice == "E-S0-12" and "tax.idtaxa_trained_md5" in violazione.dettaglio
    assente = esegui_gate("G12", Contesto(valida(_con_idtaxa(
        base, idtaxa_trained=str(tmp_path / "manca.RData")))))
    (violazione,) = assente.violazioni
    assert violazione.codice == "E-S0-12" and "tax.idtaxa_trained" in violazione.dettaglio
    assert "non esiste" in violazione.dettaglio


# --------------------------------------------------------------------------- #
# Il contenuto del classificatore addestrato                                   #
# --------------------------------------------------------------------------- #


def _salva(percorso: Path, espressione: str) -> Path:
    """Scrive un file .RData con R: ``espressione`` crea gli oggetti da salvare
    e ne restituisce i nomi.
    """
    subprocess.run(
        [str(trova_rscript()), "--vanilla", "-e",
         f"nomi <- local({{ {espressione} }}); save(list = nomi, file = '{percorso}')"],
        check=True, capture_output=True)
    return percorso


def test_un_file_che_non_e_un_classificatore_addestrato_ferma_la_validazione(con_r, tmp_path):
    """
    **Obiettivo**: Verificare che la validazione carichi il file di
    ``tax.idtaxa_trained`` e si fermi con ``E-S0-28``, in ``gates.json`` e nel
    punto di ripresa, se il file non e' un ``.RData``, se non contiene un
    oggetto di ``LearnTaxa``, o se l'oggetto non ha la tabella dei ranghi; e
    che con un classificatore valido scriva i ranghi che questo nomina e
    quelli senza una colonna.

    **Razionale scientifico e sistemistico**: Un file sbagliato si scoprirebbe
    all'assegnazione tassonomica, dopo l'inferenza delle varianti: ore di
    calcolo. Senza i nomi dei ranghi le assegnazioni non si possono riportare
    alle colonne, e ogni rango finirebbe in quella sbagliata.
    """
    dataset = sintetico.genera(tmp_path / "dati", biologici=2, negativi=0, letture=40)
    base = sintetico.configurazione(dataset, tmp_path / "out")
    testo = tmp_path / "testo.RData"
    testo.write_text("non e' un file di R\n")
    casi = [
        (testo, "non si carica come .RData"),
        (_salva(tmp_path / "altro.RData", "x <- 1:3; assign('x', x, envir = globalenv()); 'x'"),
         "0 oggetti di DECIPHER::LearnTaxa"),
        (_salva(tmp_path / "senza_ranghi.RData",
                "assign('trainingSet', structure(list(taxonomy = 'Root;'), "
                "class = c('Taxa', 'Train')), envir = globalenv()); 'trainingSet'"),
         "non ha la tabella dei ranghi"),
    ]
    for numero, (file, atteso) in enumerate(casi):
        dati = _con_idtaxa(base, addestrato=file)
        dati["io"]["out_root"] = str(tmp_path / f"out{numero}")
        esito = Esecutore(ProjectRun(valida(dati)), fino_a=Passo.S0).esegui()
        chiudi()
        assert esito.conclusione is Conclusione.ARRESTATA and esito.punto.codice == "E-S0-28"
        assert atteso in esito.punto.dettaglio and "tax.idtaxa_trained" in esito.punto.dettaglio
        gates = json.loads(
            (tmp_path / f"out{numero}" / "01_input_validation" / "gates.json").read_text())
        (g12,) = [g for g in gates["gate"] if g["gate"] == "G12"]
        assert not g12["superato"] and g12["violazioni"][0]["codice"] == "E-S0-28"

    esito = Esecutore(ProjectRun(valida(_con_idtaxa(base))), fino_a=Passo.S0).esegui()
    chiudi()
    assert esito.conclusione is Conclusione.COMPLETATA
    scritto = json.loads(
        (tmp_path / "out" / "01_input_validation" / "classificatore_addestrato.json").read_text())
    assert scritto["ranghi_scartati"] == ["subclass"] and "rootrank" in scritto["ranghi"]
    assert CATALOGO["E-S0-28"].categoria is Categoria.REVISIONE_UMANA


# --------------------------------------------------------------------------- #
# La catena: primer cercato e IdTaxa insieme                                   #
# --------------------------------------------------------------------------- #


def _dati_della_catena(cartella: Path) -> dict:
    """Un dataset sintetico con il primer preceduto da basi di sfasamento, e
    la sua configurazione con il primer cercato e IdTaxa.
    """
    dataset = sintetico.genera(
        cartella / "dati", biologici=8, negativi=1, letture=1500, primer=PRIMER, sfasamento=7)
    dati = _con_idtaxa(sintetico.configurazione(dataset, cartella / "out"))
    dati["filter"].update(
        primer_removal="search", primer_max_mismatch=1, primer_max_offset=7,
        primer_untrimmed="discard", trimLeft=0, truncLen=137)
    dati["qc"]["max_frac_untrimmed"] = 0.5
    return dati


def _esegui(dati: dict, cartella: Path, nome: str, *comandi: str) -> Path:
    """Scrive la configurazione ed esegue i comandi nella cartella indicata."""
    dati = copy.deepcopy(dati)
    dati["io"]["out_root"] = str(cartella / nome)
    config = cartella / f"{nome}.yaml"
    config.write_text(yaml.safe_dump(dati, sort_keys=False), encoding="utf-8")
    for comando in comandi or ("run",):
        chiudi()
        assert cli.main([comando, "--config", str(config)]) == 0, (nome, comando)
    chiudi()
    return cartella / nome


@pytest.fixture(scope="module")
def catena(tmp_path_factory) -> tuple[Path, dict]:
    """La catena completa, una volta per il modulo: primer cercato lungo la
    lettura e IdTaxa sul classificatore sintetico.
    """
    _richiedi_r()
    cartella = tmp_path_factory.mktemp("idtaxa")
    dati = _dati_della_catena(cartella)
    return _esegui(dati, cartella, "out", "run", "report"), dati


def test_i_ranghi_si_riconoscono_per_nome(catena):
    """
    **Obiettivo**: Verificare sulla tassonomia assegnata da IdTaxa che le
    colonne siano le sei della pipeline; che ``Root`` non compaia in alcuna
    colonna e la prima sia il dominio; che i nomi del rango in piu' del
    classificatore (``subclass``) non compaiano, e il riepilogo lo dichiari
    scartato; che nessuna etichetta ``unclassified_`` resti; che una variante
    del lignaggio senza l'ordine abbia l'ordine vuoto con famiglia e genere
    assegnati; e che ogni nome stia nella colonna del suo rango.

    **Razionale scientifico e sistemistico**: Il classificatore comincia con
    la radice e, in certi lignaggi, salta o aggiunge un rango: una
    corrispondenza per posizione metterebbe il dominio nella colonna del
    phylum e, dove un rango manca, ogni nome successivo una colonna prima.
    """
    uscita, _ = catena
    cartella = uscita / "08_taxonomy"
    righe = _tabella(cartella / "tassonomia.tsv")
    assert list(righe[0]) == ["sequenza", *COLONNE, "letture"]
    prefisso = {"Phylum": "Phylum_", "Class": "Classe_", "Order": "Ordine_",
                "Family": "Famiglia_", "Genus": "Genere_"}
    for riga in righe:
        assert riga["Kingdom"] == "Bacteria"
        for colonna, atteso in prefisso.items():
            assert riga[colonna] == "" or riga[colonna].startswith(atteso), riga
    tutti = {riga[c] for riga in righe for c in COLONNE}
    assert "Root" not in tutti
    assert not any(n.startswith(("Sottoclasse", "unclassified_")) for n in tutti)
    # Il rango in piu' e' stato attraversato: i suoi generi sono assegnati.
    assert {"Genere_08", "Genere_09"} <= tutti
    saltato = [r for r in righe if r["Family"] == "Famiglia_10"]
    assert saltato and all(r["Order"] == "" and r["Class"] == "Classe_2" for r in saltato)
    assert any(r["Genus"] for r in saltato)

    riepilogo = json.loads((cartella / "riepilogo.json").read_text(encoding="utf-8"))
    assert (riepilogo["classificatore"], riepilogo["soglia"], riepilogo["filamento"]) == (
        "idtaxa", 60.0, "top")
    assert riepilogo["ranghi"] == COLONNE and riepilogo["ranghi_scartati"] == ["subclass"]
    assert "min_boot" not in riepilogo and "try_rc" not in riepilogo
    assert riepilogo["riferimento"]["md5"] == hashlib.md5(ADDESTRATO.read_bytes()).hexdigest()
    assert riepilogo["varianti_con_phylum"] == len(righe)


def test_la_confidenza_ha_un_artefatto_proprio(catena):
    """
    **Obiettivo**: Verificare che con IdTaxa la fase scriva ``confidenza.tsv``
    e non ``bootstrap.tsv``: una riga per variante, la confidenza (fra la
    soglia e 100) dove il rango e' assegnato e vuota dove non lo e'; che la
    tabella accessoria dell'oggetto integrato porti colonne ``confidenza_``; e
    che il manifesto non riporti la versione di dada2, che qui non assegna.

    **Razionale scientifico e sistemistico**: Confidenza di IdTaxa e bootstrap
    di ``assignTaxonomy`` non sono la stessa misura: scritte sotto lo stesso
    nome, una soglia pensata per l'una verrebbe letta sull'altra.
    """
    uscita, _ = catena
    cartella = uscita / "08_taxonomy"
    assert (cartella / "confidenza.tsv").is_file() and not (cartella / "bootstrap.tsv").exists()
    tassonomia = {r["sequenza"]: r for r in _tabella(cartella / "tassonomia.tsv")}
    confidenze = _tabella(cartella / "confidenza.tsv")
    assert [r["sequenza"] for r in confidenze] == list(tassonomia)
    for riga in confidenze:
        for colonna in COLONNE:
            assegnato = tassonomia[riga["sequenza"]][colonna] != ""
            assert (riga[colonna] != "") is assegnato
            if assegnato:
                assert 60 <= float(riga[colonna]) <= 100
    accessorie = _tabella(uscita / "10_phyloseq" / "varianti_accessorie.tsv")
    assert list(accessorie[0]) == [
        "asv_id", *(f"confidenza_{c}" for c in COLONNE), "difetto_riferimento"]
    manifesto = json.loads((cartella / "manifest_S8.json").read_text(encoding="utf-8"))
    assert "dada2" not in manifesto["metriche"]
    assert {"confidenza.tsv", "ranghi_scartati.json", "tassonomia.tsv", "tassonomia.rds"} <= {
        a["nome"] for a in manifesto["artefatti"]}
    # Il report dice quale classificatore ha assegnato, e non da' la fase per
    # non conclusa perche' manca la versione di dada2.
    report = (uscita / "report" / "report.html").read_text(encoding="utf-8")
    assert "DECIPHER::IdTaxa" in report and "la fase non è conclusa" not in report
    assert manifesto["metriche"]["classificatore"] == "idtaxa"
    # La validazione ha verificato il classificatore prima di ogni calcolo.
    assert (uscita / "01_input_validation" / "classificatore_addestrato.json").is_file()


def test_il_rango_saltato_arriva_all_oggetto_finale_e_alle_analisi(catena, tmp_path):
    """
    **Obiettivo**: Verificare che la variante con l'ordine saltato resti
    nell'oggetto finale con l'ordine vuoto e famiglia e genere assegnati; che
    tutte le varianti abbiano la lunghezza del troncamento, senza il primer; e
    che nelle analisi ecologiche, raggruppando all'ordine, le sue letture
    stiano fra i non assegnati, mentre raggruppando al genere resti un taxon
    distinto con il suo lignaggio.

    **Razionale scientifico e sistemistico**: Un rango piu' profondo assegnato
    sotto un rango vuoto con il bayesiano ingenuo non esiste: chi cercasse il
    primo rango vuoto per sapere fin dove la variante e' assegnata la darebbe
    ferma alla classe, e la sua famiglia e il suo genere sparirebbero dalla
    composizione.
    """
    uscita, _ = catena
    finale = uscita / "12_final"
    righe = _tabella(finale / "tassonomia.tsv")
    (saltata,) = [r for r in righe if r["Genus"] == "Genere_21"]
    assert (saltata["Class"], saltata["Order"], saltata["Family"]) == (
        "Classe_2", "", "Famiglia_10")
    with open(finale / "sequenze.fasta", encoding="utf-8") as file:
        sequenze = [riga.strip() for riga in file if not riga.startswith(">")]
    assert sequenze and {len(s) for s in sequenze} == {137}
    conteggi = {r["asv_id"]: sum(int(v) for k, v in r.items() if k != "asv_id")
                for r in _tabella(finale / "conteggi.tsv")}
    letture_saltata = conteggi[saltata["asv_id"]]
    assert letture_saltata > 0

    def composizione(rango: str) -> list[dict[str, str]]:
        config = tmp_path / f"eco_{rango}.yaml"
        config.write_text(yaml.safe_dump(configurazione_eco(
            design={"variable": sintetico.COLONNA_GRUPPO},
            comp={"rank": rango, "top_n": 5},
            beta={"distances": ["bray"], "clr_pseudocount": ...})), encoding="utf-8")
        cartella = tmp_path / f"eco_{rango}"
        assert main_eco(["run", "--object", str(finale / "ps_final.rds"),
                         "--config", str(config), "--out", str(cartella)]) == 0
        return leggi_tsv(cartella / "composizione" / "abbondanze_relative.tsv")[1]

    al_genere = composizione("Genus")
    (taxon,) = [r for r in al_genere if r["taxon"] == "Genere_21"]
    assert taxon["lignaggio"] == "Bacteria;Phylum_B;Classe_2;NA;Famiglia_10;Genere_21"
    all_ordine = composizione("Order")
    assert not any("Famiglia_10" in r["lignaggio"] for r in all_ordine)
    (non_assegnato,) = [r for r in all_ordine if r["taxon"] == "non assegnato"]
    assert float(non_assegnato["abbondanza_media"]) > 0
    # All'ordine i taxa sono gli ordini assegnati piu' la categoria dei non assegnati.
    assert {r["taxon"] for r in all_ordine} - {"non assegnato"} == {
        r["Order"] for r in righe if r["Order"]}


@pytest.mark.lenta
def test_due_esecuzioni_con_idtaxa_danno_gli_stessi_byte(catena, tmp_path):
    """
    **Obiettivo**: Verificare che una seconda esecuzione della stessa
    configurazione, e una terza con un solo processo invece di due, diano per
    ogni artefatto di ogni fase lo stesso sha256 della prima, salvo la
    configurazione registrata; e che un altro seme cambi le confidenze.

    **Razionale scientifico e sistemistico**: IdTaxa campiona i k-meri a caso:
    senza il seme, o se il risultato dipendesse dal numero di processi, la
    stessa configurazione darebbe confidenze diverse a ogni esecuzione, e una
    variante vicina alla soglia cambierebbe assegnazione.
    """
    uscita, dati = catena
    prima = _impronte(uscita)
    assert len(prima) > 60
    seconda = _impronte(_esegui(dati, tmp_path, "seconda"))
    un_processo = copy.deepcopy(dati)
    un_processo["run"]["threads"] = 1
    terza = _impronte(_esegui(un_processo, tmp_path, "terza"))

    def senza_configurazione(impronte: dict[str, str]) -> dict[str, str]:
        return {k: v for k, v in impronte.items() if not k.startswith("00_config/")}

    assert senza_configurazione(prima) == senza_configurazione(seconda)
    assert senza_configurazione(prima) == senza_configurazione(terza)
    altro_seme = copy.deepcopy(dati)
    altro_seme["run"]["seed"] = dati["run"].get("seed", defaults.RUN_SEED) + 1
    quarta = _impronte(_esegui(altro_seme, tmp_path, "quarta"))
    assert quarta["08_taxonomy/confidenza.tsv"] != prima["08_taxonomy/confidenza.tsv"]


@pytest.mark.lenta
def test_passare_a_idtaxa_rifa_la_tassonomia_e_non_l_inferenza(con_r, tmp_path):
    """
    **Obiettivo**: Verificare che, su un'esecuzione conclusa con il bayesiano
    ingenuo, dichiarare IdTaxa e riprendere rifaccia l'assegnazione tassonomica
    e le fasi che ne dipendono, non la validazione ne' l'inferenza; che gli
    artefatti propri del classificatore precedente escano dalla cartella e dal
    suo manifesto; e che un file che non e' un classificatore, dichiarato a
    validazione conclusa, fermi l'assegnazione con ``E-S0-28``.

    **Razionale scientifico e sistemistico**: Confrontare due classificatori
    sulle stesse varianti non deve costare una seconda inferenza; e i file del
    primo, lasciati accanto a quelli del secondo, direbbero un bootstrap che
    questa esecuzione non ha calcolato.
    """
    dataset = sintetico.genera(tmp_path / "dati", biologici=6, negativi=0, letture=800)
    base = sintetico.configurazione(dataset, tmp_path / "out")
    uscita = _esegui(base, tmp_path, "out")
    cartella = uscita / "08_taxonomy"
    assert (cartella / "bootstrap.tsv").is_file()
    prima = {m.name: m.read_bytes() for m in uscita.rglob("manifest_S*.json")}
    _esegui(_con_idtaxa(base), tmp_path, "out", "resume")
    dopo = {m.name: m.read_bytes() for m in uscita.rglob("manifest_S*.json")}
    rifatte = {n.removeprefix("manifest_").removesuffix(".json")
               for n in dopo if dopo[n] != prima.get(n)}
    assert rifatte == {"S8", "S10", "S11", "S12", "S13", "S14"}
    assert (cartella / "confidenza.tsv").is_file() and not (cartella / "bootstrap.tsv").exists()
    registrati = json.loads((cartella / "manifest.json").read_text(encoding="utf-8"))
    assert "bootstrap.tsv" not in {a["nome"] for a in registrati["artefatti"]}

    falso = tmp_path / "falso.RData"
    falso.write_text("non e' un classificatore\n")
    dati = _con_idtaxa(base, addestrato=falso)
    dati["io"]["out_root"] = str(uscita)
    config = tmp_path / "falso.yaml"
    config.write_text(yaml.safe_dump(dati, sort_keys=False), encoding="utf-8")
    assert cli.main(["resume", "--config", str(config)]) == 4
    punto = json.loads((uscita / "99_logs" / "punto_di_ripresa.json").read_text(encoding="utf-8"))
    assert (punto["passo"], punto["codice"]) == ("S8", "E-S0-28")
