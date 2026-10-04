r"""Suite di test per i dati di OSD-734 recuperabili dal repository (``dati/osd734/``).

1. Inquadramento nel Piano Operativo
------------------------------------
Settimana 22 (W22), Fase F9 (riproducibilita' e consegna, anticipata: i dati del
dataset di riferimento recuperabili da un clone del repository).

2. Moduli sorgente coperti
--------------------------
* ``dati/osd734/scaricamento.py``, ``dati/osd734/scarica_letture.py``,
  ``dati/osd734/scarica_riferimento.py``, ``dati/osd734/ricostruisci_lotto.py``
* ``dati/osd734/letture_ena.tsv``, ``dati/osd734/FONTI.tsv``,
  ``dati/osd734/config_osd734.yaml``, ``.gitignore``

3. Cosa valuta questo file
--------------------------
- l'elenco delle letture copre esattamente i campioni della tabella di assay, con
  un MD5 e una dimensione per file, e nomi che contengono un solo accession;
- i file del repository corrispondono ai checksum dichiarati in ``FONTI.tsv``;
- lo scarico (su un server HTTP locale): un file integro non si scarica, uno
  interrotto riprende dal punto a cui era arrivato, uno con l'MD5 sbagliato
  viene messo da parte e non sostituisce quello atteso;
- la ricostruzione del file del lotto: classi, piastre, ordine e accession,
  anche per un campione sequenziato due volte; una fonte diversa da quella
  fissata e' respinta;
- la configurazione di esempio e' valida, con i percorsi in ``dati/osd734/``;
- i file scaricati sono esclusi da git.

4. Comandi Bash e scenari di esecuzione
---------------------------------------
    1. Modalità locale standard:
       pytest tests/test_w22_dati_osd734.py -v

    2. Modalità container Docker standard: la cartella ``dati/`` non e' copiata
       nell'immagine, quindi i test si saltano, salvo montare il repository:
       docker run --rm \
         -e PYTHONPATH=/app/src \
         -e AMPLICON16S_R_DIR=/app/R \
         -v "$(pwd)":/app \
         -w /app \
         amplicon16s:dev \
         pytest -o cache_dir=/tmp/.pytest_cache tests/test_w22_dati_osd734.py -v

    3. Modalità container Docker completa: nessun test del modulo usa i dati
       reali, che si verificano scaricandoli nella prima esecuzione completa.

5. Risultato atteso
-------------------
Vedi ``test.txt``, scheda W22.

6. Razionale scientifico e sistemistico
---------------------------------------
- Chi clona il repository deve poter riprodurre i risultati su OSD-734: senza il
  file del lotto la pipeline degrada in silenzio (un solo modello d'errore,
  nessuna soglia di profondita' per piastra).
- Un file scaricato vale solo se coincide con il checksum della sua fonte: un
  archivio troncato o alterato produrrebbe risultati diversi senza errori.
- I test non interrogano la rete: ENA, Zenodo e GitHub non sono sotto il
  controllo della suite.
"""

from __future__ import annotations

import csv
import hashlib
import re
import subprocess
import sys
import threading
from functools import partial
from http.server import HTTPServer, SimpleHTTPRequestHandler
from pathlib import Path

import pytest

RADICE = Path(__file__).resolve().parents[1]
DATI = RADICE / "dati" / "osd734"
if not DATI.is_dir():
    pytest.skip("la cartella dati/ non e' presente (nell'immagine non viene copiata)",
                allow_module_level=True)
sys.path.insert(0, str(DATI))

import ricostruisci_lotto  # noqa: E402
import scarica_letture  # noqa: E402
import scarica_riferimento  # noqa: E402
from scaricamento import Atteso, ErroreScarico, scarica, valido  # noqa: E402

from amplicon16s.config.schema import carica  # noqa: E402


def _md5(contenuto: bytes) -> str:
    return hashlib.md5(contenuto).hexdigest()


class _ConRipresa(SimpleHTTPRequestHandler):
    """Un server di file che rispetta l'intestazione ``Range`` e registra le
    richieste ricevute.
    """

    richieste: list[str | None] = []

    def log_message(self, *argomenti):  # noqa: D401 - silenzioso nei test
        """Nessun messaggio sulla console."""

    def do_GET(self):  # noqa: N802 - nome imposto da http.server
        """Serve il file intero, o il resto da un byte se richiesto."""
        intervallo = self.headers.get("Range")
        type(self).richieste.append(intervallo)
        percorso = Path(self.translate_path(self.path))
        contenuto = percorso.read_bytes()
        trovato = re.fullmatch(r"bytes=(\d+)-", intervallo or "")
        if trovato:
            inizio = int(trovato.group(1))
            self.send_response(206)
            corpo = contenuto[inizio:]
        else:
            self.send_response(200)
            corpo = contenuto
        self.send_header("Content-Length", str(len(corpo)))
        self.end_headers()
        self.wfile.write(corpo)


@pytest.fixture
def server(tmp_path):
    """Un server HTTP locale sulla cartella ``origine``; restituisce (cartella, url)."""
    origine = tmp_path / "origine"
    origine.mkdir()
    _ConRipresa.richieste = []
    istanza = HTTPServer(("127.0.0.1", 0), partial(_ConRipresa, directory=str(origine)))
    filo = threading.Thread(target=istanza.serve_forever, daemon=True)
    filo.start()
    yield origine, f"http://127.0.0.1:{istanza.server_address[1]}"
    istanza.shutdown()


# --------------------------------------------------------------------------- #
# Elenchi e checksum                                                           #
# --------------------------------------------------------------------------- #


def test_l_elenco_delle_letture_copre_la_tabella_di_assay():
    """
    **Obiettivo**: Verificare che ``letture_ena.tsv`` abbia una riga per ciascuno
    dei 960 campioni della tabella di assay, con accession univoci, un MD5 e
    una dimensione per file, un indirizzo ENA, e un nome di file che contiene
    un solo accession dell'esperimento, quello della riga.

    **Razionale scientifico e sistemistico**: La pipeline estrae l'accession dal
    nome del file (``io.accession_regex``): un nome con due accession, o con uno
    diverso, romperebbe il join con i metadati.
    """
    righe = scarica_letture.leggi_elenco()
    assay = ricostruisci_lotto.accession_per_nome()
    assert len(righe) == 960
    assert {r["experiment_accession"] for r in righe} == set(assay.values())
    assert len({r["run_accession"] for r in righe}) == 960
    for riga in righe:
        assert re.findall(r"(?:E|S|D)RX[0-9]{4,}", riga["file"]) == [riga["experiment_accession"]]
        assert re.fullmatch(r"[0-9a-f]{32}", riga["fastq_md5"])
        assert int(riga["fastq_bytes"]) > 0
        assert riga["fastq_url"].startswith("https://ftp.sra.ebi.ac.uk/")


def test_i_file_del_repository_corrispondono_alle_fonti():
    """
    **Obiettivo**: Verificare che ogni file di ``FONTI.tsv`` presente nel
    repository abbia l'MD5 dichiarato, con fonte e licenza indicate, e che il
    file del lotto risulti non ridistribuito.

    **Razionale scientifico e sistemistico**: Un metadato modificato per errore
    cambierebbe l'inventario dei campioni; la provenienza di ogni file deve
    restare verificabile.
    """
    with open(DATI / "FONTI.tsv", encoding="utf-8", newline="") as file:
        fonti = list(csv.DictReader(file, delimiter="\t"))
    nel_repository = [f for f in fonti if f["nel_repository"] == "si"]
    assert len(nel_repository) == 3
    for fonte in fonti:
        assert fonte["fonte"] and fonte["licenza"], fonte["percorso"]
    for fonte in nel_repository:
        contenuto = (DATI / fonte["percorso"]).read_bytes()
        assert fonte["controllo"] == "md5:" + _md5(contenuto), fonte["percorso"]
    lotto = next(f for f in fonti if f["percorso"].startswith("lotto/"))
    assert lotto["nel_repository"] == "no" and "non ridistribuito" in lotto["licenza"]


# --------------------------------------------------------------------------- #
# Scarico verificato e riprendibile                                            #
# --------------------------------------------------------------------------- #


def test_lo_scarico_riprende_dal_punto_interrotto(server, tmp_path):
    """
    **Obiettivo**: Verificare che, con un file ``.parziale`` gia' presente, lo
    scarico chieda al server solo il resto (``Range``), ricomponga il file e
    lo verifichi; e che un secondo avvio non chieda nulla, perche' il file e'
    integro.

    **Razionale scientifico e sistemistico**: Scaricare 2,4 GB da un archivio
    pubblico deve poter riprendere dopo un'interruzione, senza ricominciare e
    senza lasciare un file incompleto con il nome definitivo.
    """
    origine, url = server
    contenuto = bytes(range(256)) * 400
    (origine / "lettura.fastq.gz").write_bytes(contenuto)
    atteso = Atteso(f"{url}/lettura.fastq.gz", _md5(contenuto), len(contenuto))
    destinazione = tmp_path / "fastq" / "lettura.fastq.gz"
    destinazione.parent.mkdir()
    (destinazione.parent / "lettura.fastq.gz.parziale").write_bytes(contenuto[:30000])

    scarica(atteso, destinazione, pausa=0)
    assert destinazione.read_bytes() == contenuto
    assert _ConRipresa.richieste == ["bytes=30000-"]
    assert not (destinazione.parent / "lettura.fastq.gz.parziale").exists()
    assert valido(destinazione, atteso)


def test_un_file_con_l_md5_sbagliato_non_diventa_il_file_atteso(server, tmp_path):
    """
    **Obiettivo**: Verificare che, se la fonte serve un contenuto diverso dal
    checksum dichiarato, lo scarico fallisca con ``ErroreScarico`` dopo i
    tentativi, che il file definitivo non esista e che il contenuto sbagliato
    resti da parte con il suffisso ``.md5_errato``.

    **Razionale scientifico e sistemistico**: Un file alterato con il nome
    giusto passerebbe inosservato fino ai risultati: deve fermare lo script,
    non sostituire in silenzio quello atteso.
    """
    origine, url = server
    (origine / "lettura.fastq.gz").write_bytes(b"contenuto alterato")
    atteso = Atteso(f"{url}/lettura.fastq.gz", _md5(b"contenuto originale"), 18)
    destinazione = tmp_path / "lettura.fastq.gz"
    with pytest.raises(ErroreScarico, match="MD5"):
        scarica(atteso, destinazione, tentativi=2, pausa=0)
    assert not destinazione.exists()
    assert (tmp_path / "lettura.fastq.gz.parziale.md5_errato").read_bytes() == b"contenuto alterato"


def test_lo_script_non_scarica_i_file_gia_integri(tmp_path, monkeypatch):
    """
    **Obiettivo**: Verificare che ``scarica_letture.py``, su una cartella che
    contiene gia' il file integro, lo riconosca senza alcuna richiesta di rete
    ed esca con 0; e che con ``--solo-verifica`` un file assente faccia uscire
    con 1 senza creare nulla.

    **Razionale scientifico e sistemistico**: Le letture gia' presenti (per
    esempio collegate da un'altra cartella) non si riscaricano: e' un peso
    inutile per un servizio pubblico.
    """
    contenuto = b"@lettura\nACGT\n+\nIIII\n"
    elenco = tmp_path / "elenco.tsv"
    elenco.write_text(
        "file\trun_accession\texperiment_accession\tsample_name\tfastq_url\tfastq_md5\tfastq_bytes\n"
        f"ERX1_A_ERR1.fastq.gz\tERR1\tERX1\tA\thttp://127.0.0.1:9/x\t{_md5(contenuto)}\t{len(contenuto)}\n",
        encoding="utf-8",
    )
    monkeypatch.setattr(scarica_letture, "ELENCO", elenco)
    cartella = tmp_path / "fastq"
    assert scarica_letture.main(["--cartella", str(cartella), "--solo-verifica"]) == 1
    assert not cartella.exists()
    cartella.mkdir()
    (cartella / "ERX1_A_ERR1.fastq.gz").write_bytes(contenuto)
    assert scarica_letture.main(["--cartella", str(cartella), "--pausa", "0"]) == 0


# --------------------------------------------------------------------------- #
# Ricostruzione del file del lotto                                             #
# --------------------------------------------------------------------------- #


def _fonte(*righe: tuple[str, ...]) -> bytes:
    """Una fonte sintetica nel formato dei metadati Qiita degli autori."""
    colonne = ("#SampleID", "sample_type", "sample_plate", "well_id", "primer_plate",
               "run_prefix", "run_date", "module", "katharoseq_cell_count")
    testo = "\t".join(colonne) + "\n" + "".join("\t".join(r) + "\n" for r in righe)
    return testo.encode("utf-8")


def test_il_lotto_si_ricostruisce_dalla_fonte(tmp_path):
    """
    **Obiettivo**: Verificare su una fonte sintetica che la ricostruzione tolga
    il prefisso di Qiita (diverso per corsa), mappi le classi, ricavi numero di
    piastra e data dal nome della piastra, ordini per piastra e pozzetto, e
    attribuisca a un campione sequenziato due volte l'accession della replica
    nell'intervallo della sua corsa.

    **Razionale scientifico e sistemistico**: Il nome del campione si ripete
    fra le repliche: senza l'accession giusto, una replica prenderebbe il
    lotto dell'altra senza alcun errore.
    """
    p1, p6 = "3DMM_20220614_14542_Plate_1", "3DMM_20220623_14542_Plate_6"
    c1, c2 = "3DMM_Plates_1-5_S1_L001", "3DMM_rerun_Plates_6-10_S1_L001"
    fonte = _fonte(
        ("150297.14542.LAB1.F1", "surface swab", p1, "B1", "1.0", c1, "6/8/22", "LAB 1", ""),
        ("150297.14542.LAB1P3.L1", "surface swab", p1, "A10", "1.0", c1, "6/8/22", "LAB 1", ""),
        ("150297.14542.KATHARO.P1.1", "control positive", p1, "A2", "1.0", c1, "6/8/22",
         "not applicable", "2000000"),
        ("151322.14542.LAB1P3.L1", "surface swab", p6, "C3", "4.0", c2, "6/16/22", "LAB 1", ""),
        ("151322.14542.BLANK.3DMM.P6.1A", "control blank", p6, "A1", "4.0", c2, "6/16/22",
         "not applicable", ""),
        ("151322.14542.NOD1.F1", "surface swab", p6, "A2", "4.0", c2, "6/16/22", "Node 1", ""),
    )
    per_nome = {"LAB1.F1": "ERX10", "LAB1P3.L1_rep1": "ERX11", "KATHARO.P1.1": "ERX12",
                "LAB1P3.L1_rep2": "ERX21", "BLANK.3DMM.P6.1A": "ERX20", "NOD1.F1": "ERX22"}
    lotto = ricostruisci_lotto.mappa_lotto(
        list(csv.DictReader(fonte.decode().splitlines(), delimiter="\t")))
    ricostruisci_lotto.aggiungi_accession(lotto, per_nome)
    assert [(r["sample_name_ena"], r["well_id"]) for r in lotto] == [
        ("LAB1P3.L1", "A10"), ("KATHARO.P1.1", "A2"), ("LAB1.F1", "B1"),
        ("BLANK.3DMM.P6.1A", "A1"), ("NOD1.F1", "A2"), ("LAB1P3.L1", "C3")]
    assert [r["class"] for r in lotto] == ["SWAB", "KATHARO", "SWAB", "BLANK", "SWAB", "SWAB"]
    assert [r["experiment_accession"] for r in lotto] == [
        "ERX11", "ERX12", "ERX10", "ERX20", "ERX22", "ERX21"]
    assert (lotto[3]["extraction_plate_num"], lotto[3]["extraction_date"]) == ("6", "20220623")
    uscita = tmp_path / "lotto.tsv"
    ricostruisci_lotto.scrivi(lotto, uscita)
    assert uscita.read_text(encoding="utf-8").splitlines()[0].split("\t") == list(
        ricostruisci_lotto.COLONNE)


def test_una_fonte_diversa_da_quella_fissata_e_respinta():
    """
    **Obiettivo**: Verificare che un contenuto con un'impronta git diversa da
    quella del commit fissato fermi la ricostruzione, e che l'impronta si
    calcoli come la calcola git.

    **Razionale scientifico e sistemistico**: Il file degli autori puo' cambiare
    nel loro repository: il lotto deve venire sempre dalla stessa versione.
    """
    vuoto = subprocess.run(["git", "hash-object", "--stdin"], input=b"x\n",
                           capture_output=True, check=False)
    if vuoto.returncode == 0:
        assert ricostruisci_lotto.impronta_git(b"x\n") == vuoto.stdout.decode().strip()
    with pytest.raises(SystemExit, match="commit fissato"):
        ricostruisci_lotto.leggi_fonte(_fonte())


# --------------------------------------------------------------------------- #
# Configurazione ed esclusioni                                                 #
# --------------------------------------------------------------------------- #


def test_la_configurazione_di_esempio_usa_la_cartella_dei_dati():
    """
    **Obiettivo**: Verificare che ``config_osd734.yaml`` sia valida, che i
    percorsi dei dati stiano in ``dati/osd734/`` e siano relativi, che il lotto
    sia il file ricostruito e che ``tax.ref_md5`` sia l'MD5 che lo script del
    riferimento verifica.

    **Razionale scientifico e sistemistico**: La configurazione e gli script
    devono indicare gli stessi file: un percorso o un checksum diversi farebbero
    girare la pipeline su un riferimento non verificato.
    """
    config = carica(DATI / "config_osd734.yaml")
    for percorso in (config.io.fastq_dir, config.io.assay_table, config.io.study_table,
                     config.io.batch_table, config.tax.ref_fasta, config.tax.ref_bad_taxa):
        assert not percorso.is_absolute() and percorso.parts[:2] == ("dati", "osd734"), percorso
    assert RADICE / config.io.batch_table == ricostruisci_lotto.USCITA
    assert config.tax.ref_md5 == scarica_riferimento.FILE[config.tax.ref_fasta.name].md5
    assert (RADICE / config.io.assay_table).is_file() and (RADICE / config.io.study_table).is_file()


def test_i_dati_scaricati_sono_esclusi_da_git():
    """
    **Obiettivo**: Verificare che letture, riferimento e lotto scaricati in
    ``dati/osd734/`` siano esclusi da git, e che i file del repository no.

    **Razionale scientifico e sistemistico**: Gigabyte di letture o un file
    senza licenza di ridistribuzione non devono poter entrare in git per errore.
    """
    if not (RADICE / ".git").exists():
        pytest.skip("non e' un repository git")
    esclusi = ["dati/osd734/fastq/ERX1_A_ERR1.fastq.gz",
               "dati/osd734/riferimento/silva_nr99_v138_train_set.fa.gz",
               "dati/osd734/lotto/plate_well_map_960.tsv",
               "dati/osd734/fastq/ERX1_A_ERR1.fastq.gz.parziale"]
    tenuti = ["dati/osd734/letture_ena.tsv", "dati/osd734/metadati/s_OSD-734.txt"]
    for percorso in esclusi:
        assert subprocess.run(["git", "check-ignore", "-q", percorso], cwd=RADICE).returncode == 0, percorso
    for percorso in tenuti:
        assert subprocess.run(["git", "check-ignore", "-q", percorso], cwd=RADICE).returncode == 1, percorso
