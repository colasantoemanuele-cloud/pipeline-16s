r"""Suite di test della settimana 30: correzioni a lettura dei file, memoria e validazione.

1. Inquadramento nel Piano Operativo
------------------------------------
Settimana 30 (W30), Fase F8 (generalita'): archivi illeggibili e righe vuote
nei file di letture, memoria esaurita nel profilo, osservazioni minori su
configurazione, gate, console e permessi dei file.

2. Moduli sorgente coperti
--------------------------
* ``src/amplicon16s/io_layer/reads.py`` (``scansiona_file``, ``senza_righe_vuote_finali``)
* ``src/amplicon16s/steps/s01_profile.py``, ``R/01_profile.R``, ``R/lib/errors.R``
* ``src/amplicon16s/steps/s02_filter.py`` (``esamina_archivio``)
* ``src/amplicon16s/steps/s14_finale.py``
* ``src/amplicon16s/rbridge/runner.py`` (``classifica``), ``src/amplicon16s/rbridge/payload.py``
* ``src/amplicon16s/errors/catalog.py`` (``E-S1-05``, ``E-S1-06``)
* ``src/amplicon16s/config/schema.py`` (interi in senso stretto)
* ``src/amplicon16s/gates/g01_g15.py`` (G02, G08), ``src/amplicon16s/metadata/tabelle.py``,
  ``src/amplicon16s/metadata/crosswalk.py``
* ``src/amplicon16s/logging/logger.py`` (console), ``src/amplicon16s/io_layer/artifacts.py``
* ``dati/osd276/scarica_letture.py``

3. Cosa valuta questo file
--------------------------
- ogni guasto di decompressione di un file di letture diventa un esito
  dichiarato, mai un'eccezione, e ferma S1 con ``E-S1-05`` e l'indicazione di
  riscaricare il file;
- le righe vuote in fondo a un FASTQ sono tollerate (i lettori di R ricevono
  una copia senza quelle righe, il file di ingresso non si tocca); una riga
  vuota in mezzo e un record troncato restano ``E-S1-05``, con messaggi distinti;
- un processo figlio del calcolo parallelo di R ucciso dal sistema e' letto
  dal ponte come memoria esaurita, e S1 lo dichiara con ``E-S1-06``;
- gli interi della configurazione respingono ``true`` e ``false``;
- G02 non attribuisce alla tabella di studio un guasto della tabella di assay,
  e nomina la riga con una tabulazione nell'identificativo;
- un file del lotto non UTF-8 ferma G08 con ``E-S0-08``;
- S14 si ferma con ``E-S14-01`` se ``ps_filtrato.rds`` manca;
- due avvisi con lo stesso codice sono righe distinte sulla console;
- i file scritti in modo atomico hanno i permessi della umask dell'utente.

4. Comandi Bash e scenari di esecuzione
---------------------------------------
    ``<immagine>`` e' l'immagine del container della pipeline; quella corrente
    e' indicata in ``README.md``.

    1. Modalità locale standard (R di base con jsonlite, senza Bioconductor né
       dati reali):
       pytest tests/test_w30_correzioni.py -v

    2. Modalità container Docker standard (sottoinsieme ridotto con
       R/Bioconductor):
       docker run --rm \
         -e PYTHONPATH=/app/src \
         -e AMPLICON16S_R_DIR=/app/R \
         -v "$(pwd)":/app \
         -w /app \
         <immagine> \
         pytest -o cache_dir=/tmp/.pytest_cache tests/test_w30_correzioni.py -v

    3. Modalità container Docker completa (con i dati reali OSD-734):
       docker run --rm \
         --memory=24g \
         --memory-swap=24g \
         -e PYTHONPATH=/app/src \
         -e AMPLICON16S_R_DIR=/app/R \
         -e AMPLICON16S_CONFIG_DATI_REALI="$HOME/ASI/config_osd734.yaml" \
         -v "$(pwd)":/app \
         -v "$HOME/ASI":"$HOME/ASI" \
         -w /app \
         <immagine> \
         pytest -o cache_dir=/tmp/.pytest_cache tests/test_w30_correzioni.py -v

5. Risultato atteso
-------------------
I conteggi li da' pytest (``pytest --collect-only -q``).

6. Razionale scientifico e sistemistico
---------------------------------------
- Un file di letture corrotto o troncato usato in parte da' abbondanze sbagliate
  senza alcun errore: ogni guasto di lettura deve fermare con la sua causa.
- Una riga vuota in coda non altera alcuna lettura: fermare l'esecuzione
  sarebbe un rifiuto senza motivo, tollerarla senza ripulirla romperebbe i
  lettori di R con un errore fuori catalogo.
- La memoria esaurita che arriva come errore generico manda chi esegue a
  cercare un difetto dello script che non c'e'.
- Un refuso nella configurazione (``true`` per un intero) non deve diventare
  un valore valido.
"""

from __future__ import annotations

import gzip
import importlib.util
import logging
import os
import random
import stat
import zlib
from pathlib import Path

import pytest
from conftest import NEGATIVO, POSITIVO, Campione, copia_esecuzione, crea_scenario, lettura
from sottoinsieme import motivo_pacchetti_r_assenti

from amplicon16s.config.schema import ErroreConfigurazione, valida
from amplicon16s.errors.catalog import CATALOGO, Categoria
from amplicon16s.errors.exceptions import ErrorePipeline
from amplicon16s.gates.g01_g15 import Contesto
from amplicon16s.gates.registry import esegui_gate, esegui_tutti
from amplicon16s.io_layer import reads
from amplicon16s.io_layer.artifacts import AlberoOutput, Fase, scrivi_atomico
from amplicon16s.io_layer.conteggi import leggi_conteggi
from amplicon16s.io_layer.reads import scansiona_file, senza_righe_vuote_finali
from amplicon16s.logging.logger import _FormattatoreConsole, chiudi
from amplicon16s.rbridge.payload import Dichiarazione, Stato
from amplicon16s.rbridge.runner import Condizione, classifica, esegui_script, trova_rscript
from amplicon16s.runner.executor import Conclusione, Esecutore
from amplicon16s.runner.graph import Passo
from amplicon16s.runner.project import ProjectRun
from amplicon16s.steps.s02_filter import esamina_archivio
from amplicon16s.steps.s13_filtri import NOME_FILTRATO

RADICE = Path(__file__).resolve().parents[1]


@pytest.fixture(autouse=True)
def uscite_pulite():
    """Chiude le uscite del log prima e dopo ogni test, perché nessun handler resti
    aperto sulla cartella temporanea.
    """
    chiudi()
    yield
    chiudi()


@pytest.fixture
def bioc():
    """Richiede R con i pacchetti delle fasi: salta senza, ma in CI fallisce."""
    motivo = motivo_pacchetti_r_assenti(
        "dada2", "ggplot2", "ShortRead", "jsonlite", "phyloseq", "Biostrings"
    )
    if motivo is not None:
        if os.environ.get("AMPLICON16S_RICHIEDI_BIOC") == "1":
            pytest.fail(f"Bioconductor e' richiesto in questo ambiente: {motivo}")
        pytest.skip(motivo)


@pytest.fixture
def r_di_base():
    """Richiede Rscript con jsonlite e parallel: salta senza."""
    motivo = motivo_pacchetti_r_assenti("jsonlite", "parallel")
    if motivo is not None:
        pytest.skip(motivo)


def _record(n: int, lunghezza: int = 151, seme: int = 0) -> str:
    """``n`` record FASTQ con sequenze casuali (seme fisso), poco comprimibili."""
    generatore = random.Random(seme)
    righe = []
    for i in range(n):
        sequenza = "".join(generatore.choice("ACGT") for _ in range(lunghezza))
        righe.append(f"@lettura{i}\n{sequenza}\n+\n{'I' * lunghezza}\n")
    return "".join(righe)


def _scrivi(percorso: Path, testo: str) -> Path:
    """Scrive ``testo`` come FASTQ compresso."""
    with gzip.open(percorso, "wt", encoding="utf-8", newline="\n") as file:
        file.write(testo)
    return percorso


def _corruzione_zlib(percorso: Path) -> bytes:
    """Il contenuto di ``percorso`` con un byte alterato nell'ultimo terzo, scelto
    come il primo che fa fallire la decompressione con ``zlib.error``.
    """
    originale = percorso.read_bytes()
    for posizione in range(2 * len(originale) // 3, len(originale) - 8):
        alterato = bytearray(originale)
        alterato[posizione] ^= 0xFF
        try:
            gzip.decompress(bytes(alterato))
        except zlib.error:
            return bytes(alterato)
        except (OSError, EOFError):
            continue
    raise AssertionError("nessuna alterazione produce zlib.error")


def _campioni() -> list[Campione]:
    """Due biologici, un positivo e un negativo."""
    return [
        Campione("ERX3000001", "NOD1D4.L1"),
        Campione("ERX3000002", "NOD1D4.L2"),
        Campione("ERX3000003", "POS.P1.1", materiale=POSITIVO, posizione="Not Applicable"),
        Campione("ERX3000004", "BLANK.P1.1", materiale=NEGATIVO, posizione="Not Applicable"),
    ]


def _file_di(scenario, accession: str) -> Path:
    """Il file FASTQ di un campione dello scenario."""
    (percorso,) = Path(scenario.config.io.fastq_dir).glob(f"*{accession}*")
    return percorso


# --------------------------------------------------------------------------- #
# 1. Archivi illeggibili e righe vuote                                         #
# --------------------------------------------------------------------------- #


def test_un_errore_di_zlib_e_un_esito_dichiarato_non_un_eccezione(tmp_path, monkeypatch):
    """
    **Obiettivo**: Verificare che la scansione di un archivio il cui flusso
    compresso e' alterato (``zlib.error``, che non e' un ``OSError``)
    restituisca un esito con l'errore «archivio non leggibile» e il tipo del
    guasto, senza sollevare; e che su duecento alterazioni di un byte la
    scansione non sollevi mai.

    **Razionale scientifico e sistemistico**: Un'eccezione fuori catalogo
    ferma l'esecuzione senza dire quale file riscaricare: il guasto di
    decompressione dev'essere un esito del file, qualunque forma prenda.
    """
    percorso = _scrivi(tmp_path / "a.fastq.gz", _record(400))
    originale = percorso.read_bytes()
    percorso.write_bytes(_corruzione_zlib(percorso))
    esito = scansiona_file(percorso, 10**9, None, None)
    assert esito.errore is not None and esito.errore.startswith("archivio non leggibile dopo ")
    assert "error" in esito.errore and not esito.valido

    generatore = random.Random(1)
    for _ in range(200):
        alterato = bytearray(originale)
        alterato[generatore.randrange(10, len(alterato))] ^= generatore.randrange(1, 256)
        percorso.write_bytes(bytes(alterato))
        scansiona_file(percorso, 10**9, None, None)  # non deve sollevare

    def guasto(_percorso):
        raise zlib.error("Error -3 while decompressing data: invalid block type")

    monkeypatch.setattr(reads, "apri_fastq", guasto)
    esito = scansiona_file(percorso, 10, None, None)
    assert "zlib.error" in esito.errore or "error: Error -3" in esito.errore


def test_le_righe_vuote_in_fondo_sono_tollerate_e_quelle_in_mezzo_no(tmp_path):
    """
    **Obiettivo**: Verificare che un file che si chiude con righe vuote (una,
    piu' d'una, con spazi, con ritorno carrello) sia valido, con tutte le sue
    letture contate e ``righe_vuote_finali`` vero; che una riga vuota seguita
    da altri record dia l'errore «riga vuota in mezzo al file»; che un record
    incompleto dia «record troncato»; e che un file senza righe vuote non le
    dichiari.

    **Razionale scientifico e sistemistico**: Una riga vuota in coda non
    tocca alcuna lettura; una in mezzo sfasa i record di chi legge a blocchi
    di quattro righe. Sono due difetti diversi dal record troncato, e il
    messaggio deve dire quale.
    """
    base = _record(30, 20)
    for coda in ("\n", "\n\n\n", "  \n", "\r\n\r\n"):
        esito = scansiona_file(_scrivi(tmp_path / "v.fastq.gz", base + coda), 10**9, None, None)
        assert esito.valido and esito.letture_esaminate == 30 and esito.righe_vuote_finali, repr(coda)
    pulito = scansiona_file(_scrivi(tmp_path / "p.fastq.gz", base), 10**9, None, None)
    assert pulito.valido and not pulito.righe_vuote_finali

    in_mezzo = scansiona_file(
        _scrivi(tmp_path / "m.fastq.gz", base + "\n" + _record(5, 20)), 10**9, None, None)
    assert in_mezzo.errore == "riga vuota in mezzo al file dopo 30 letture: seguono altre righe"
    troncato = scansiona_file(
        _scrivi(tmp_path / "t.fastq.gz", base + "@ultima\nACGT\n"), 10**9, None, None)
    assert troncato.errore.startswith("record troncato dopo 30 letture")
    assert "riga vuota" not in troncato.errore


def test_la_copia_senza_righe_vuote_non_tocca_il_file_di_ingresso(tmp_path):
    """
    **Obiettivo**: Verificare che ``senza_righe_vuote_finali`` restituisca gli
    stessi percorsi quando nessun file e' indicato, e per un file indicato una
    copia temporanea con le stesse letture e senza righe vuote, lasciando il
    file di ingresso identico byte per byte; che la copia sparisca all'uscita
    dal blocco; e che ``esamina_archivio`` di S2 riconosca le righe vuote
    finali, anche con il ritorno carrello, e non le veda dove non ci sono.

    **Razionale scientifico e sistemistico**: Nessuna fase modifica i dati di
    ingresso: la tolleranza si realizza su una copia, che vive quanto il
    calcolo che la legge.
    """
    testo = _record(12, 20)
    con = _scrivi(tmp_path / "con.fastq.gz", testo + "\n\n")
    senza = _scrivi(tmp_path / "senza.fastq.gz", testo)
    prima = con.read_bytes()
    percorsi = {"a": str(con), "b": str(senza)}
    with senza_righe_vuote_finali(percorsi, []) as stessi:
        assert stessi is percorsi
    with senza_righe_vuote_finali(percorsi, ["a"]) as ripuliti:
        assert ripuliti["b"] == str(senza) and ripuliti["a"] != str(con)
        copia = Path(ripuliti["a"])
        assert copia.name == con.name
        assert gzip.decompress(copia.read_bytes()).decode() == testo
    assert not copia.exists() and con.read_bytes() == prima

    assert esamina_archivio(con) == (None, True)
    assert esamina_archivio(senza) == (None, False)
    ritorno = _scrivi(tmp_path / "crlf.fastq.gz", testo.replace("\n", "\r\n") + "\r\n")
    assert esamina_archivio(ritorno) == (None, True)
    assert esamina_archivio(_scrivi(tmp_path / "crlf2.fastq.gz", testo.replace("\n", "\r\n"))) == (
        None, False)


def test_un_archivio_corrotto_oltre_le_letture_ispezionate_ferma_s1(tmp_path):
    """
    **Obiettivo**: Verificare che un archivio con il flusso compresso alterato
    oltre le letture ispezionate da S0 superi la validazione e fermi S1 con
    ``E-S1-05``, a revisione umana, con il nome del file, «archivio non
    leggibile» e l'indicazione di riscaricarlo; e che una riga vuota in mezzo
    al file si fermi allo stesso modo, con il proprio messaggio.

    **Razionale scientifico e sistemistico**: S0 legge solo le prime letture:
    il guasto piu' avanti nel file va trovato dalla prima fase che lo legge
    per intero, con un codice del catalogo e non con un'eccezione.
    """
    record = "".join(
        f"@r{i}\n{lettura()}\n+\n{'I' * 151}\n" for i in range(40)) + _record(400, 151, seme=3)
    for nome, prepara, atteso in (
        ("zlib", lambda p: p.write_bytes(_corruzione_zlib(_scrivi(p, record))),
         "archivio non leggibile"),
        ("vuota", lambda p: _scrivi(p, record + "\n" + record), "riga vuota in mezzo al file"),
    ):
        scenario = crea_scenario(tmp_path / nome, _campioni(), con_arricchimento=True,
                                 con_letture=True, sovrascrivi={"qc": {"head_reads": 20}})
        percorso = _file_di(scenario, "ERX3000002")
        prepara(percorso)
        run = ProjectRun(scenario.config)
        esito = Esecutore(run, fino_a=Passo.S1).esegui()
        assert [r.passo for r in esito.eseguite] == [Passo.S0], nome
        assert esito.conclusione is Conclusione.ARRESTATA
        punto = esito.punto
        assert (punto.passo, punto.codice, punto.categoria) == (
            Passo.S1, "E-S1-05", "revisione_umana"), nome
        assert percorso.name in punto.dettaglio and atteso in punto.dettaglio
        assert "Riscarica i file dalla sorgente" in punto.dettaglio
        chiudi()
    assert "riga vuota in mezzo" in CATALOGO["E-S1-05"].azione
    assert "Riscarica" in CATALOGO["E-S1-05"].azione


def test_s1_e_s2_leggono_un_file_che_si_chiude_con_righe_vuote(bioc, tmp_path):
    """
    **Obiettivo**: Verificare che, con un file di letture che si chiude con
    righe vuote, S1 e S2 si concludano; che S1 conti tutte le letture del
    file, quante ne conta per un file identico senza righe vuote; che S2 le
    filtri; e che il file di ingresso resti identico byte per byte.

    **Razionale scientifico e sistemistico**: I lettori a blocchi di R
    rifiutano una riga vuota in coda come un record malformato: senza la
    copia ripulita la tolleranza dichiarata diventerebbe un errore generico
    di R nella fase successiva.
    """
    scenario = crea_scenario(tmp_path, _campioni(), con_arricchimento=True, con_letture=True)
    percorso = _file_di(scenario, "ERX3000002")
    testo = gzip.decompress(percorso.read_bytes()).decode()
    _scrivi(percorso, testo + "\n\n")
    prima = percorso.read_bytes()
    run = ProjectRun(scenario.config)
    assert Esecutore(run, fino_a=Passo.S1).esegui().conclusione is Conclusione.COMPLETATA
    assert Esecutore(run, fino_a=Passo.S2).esegui().conclusione is Conclusione.COMPLETATA
    assert percorso.read_bytes() == prima
    grezze = leggi_conteggi(run.albero.cartella(Fase.QC_PROFILES) / "letture_grezze.tsv")
    assert grezze["ERX3000002"] == grezze["ERX3000001"] == testo.count("\n") // 4
    assert (run.albero.cartella(Fase.FILTERED) / "ERX3000002_filt.fastq.gz").is_file()


# --------------------------------------------------------------------------- #
# 2. Memoria esaurita nel calcolo parallelo                                    #
# --------------------------------------------------------------------------- #


def _non_catalogato(messaggio: str) -> Dichiarazione:
    """La dichiarazione di un errore che lo script non ha ricondotto a un codice."""
    return Dichiarazione(stato=Stato.ERRORE_NON_CATALOGATO, codice=None, messaggio=messaggio,
                         artefatti=())


def test_i_figli_uccisi_del_calcolo_parallelo_sono_memoria_esaurita():
    """
    **Obiettivo**: Verificare che il ponte classifichi come memoria esaurita
    un errore non catalogato preceduto, sull'uscita di errore, dall'avviso
    con cui ``parallel::mclapply`` dichiara i figli senza risultato (al
    singolare e al plurale, anche nella forma dei calcoli prenotati); che
    riconosca il messaggio di ``verifica_figli``; e che lo stesso errore senza
    quell'avviso resti non catalogato.

    **Razionale scientifico e sistemistico**: Quando il sistema uccide un
    processo figlio per memoria, il genitore prosegue con un risultato nullo e
    fallisce piu' avanti con un errore qualunque: la causa sta nell'avviso, e
    senza leggerlo il guasto risulta un difetto dello script.
    """
    errore = _non_catalogato("argument must be coercible to non-negative integer")
    for avviso in (
        "Warning in parallel::mclapply(campioni, profila) :\n  7 parallel function calls did "
        "not deliver results\n",
        "Warning: 1 parallel function call did not deliver a result\n",
        "scheduled cores 3, 5 did not deliver results, all values of the jobs will be affected\n",
    ):
        assert classifica(3, errore, avviso) is Condizione.MEMORIA_ESAURITA, avviso
    assert classifica(3, errore, "Warning: NAs introduced by coercion\n") is (
        Condizione.ERRORE_NON_CATALOGATO)
    dichiarato = _non_catalogato(
        "6 processi paralleli su 15 non hanno restituito un risultato (a, b): interrotti "
        "dal sistema operativo, di norma per memoria esaurita")
    assert classifica(3, dichiarato, "") is Condizione.MEMORIA_ESAURITA


def test_un_figlio_ucciso_in_r_arriva_con_il_codice_di_memoria_della_fase(r_di_base, tmp_path):
    """
    **Obiettivo**: Verificare, con un vero processo R in cui due dei quattro
    figli di ``mclapply`` si uccidono con ``SIGKILL``, che ``verifica_figli``
    fermi lo script, che il ponte lo legga come memoria esaurita e sollevi il
    codice indicato dalla fase (``E-S1-06``), a revisione umana, con il numero
    dei figli mancanti nel dettaglio; e che ``E-S1-06`` suggerisca piu'
    memoria o meno processi.

    **Razionale scientifico e sistemistico**: E' il caso osservato con il
    container limitato: il profilo legge un file per processo e il sistema ne
    uccide alcuni. Il codice deve dire la causa e il rimedio, non «difetto
    dello script».
    """
    script = tmp_path / "figli.R"
    script.write_text(
        'for (f in c("io_json.R", "errors.R")) source(file.path(Sys.getenv("AMPLICON16S_R_LIB"), f))\n'
        "esegui_fase(function(parametri, cartella) {\n"
        "  esiti <- parallel::mclapply(1:4, function(i) {\n"
        "    if (i %% 2L == 0L) tools::pskill(Sys.getpid(), tools::SIGKILL)\n"
        "    i\n"
        "  }, mc.cores = 4L, mc.preschedule = FALSE)\n"
        "  verifica_figli(esiti, paste0('campione', 1:4))\n"
        "  character()\n"
        "})\n",
        encoding="utf-8",
    )
    with pytest.raises(ErrorePipeline) as info:
        esegui_script(script, {}, AlberoOutput(tmp_path / "out"), Fase.QC_PROFILES,
                      passo="S1", codice_memoria="E-S1-06")
    assert info.value.codice == "E-S1-06"
    assert "2 processi paralleli su 4 non hanno restituito un risultato" in info.value.dettaglio
    assert "campione2, campione4" in info.value.dettaglio
    voce = CATALOGO["E-S1-06"]
    assert voce.categoria is Categoria.REVISIONE_UMANA and voce.fase == "S1"
    assert "run.threads" in voce.azione and "memoria" in voce.azione


# --------------------------------------------------------------------------- #
# 3. Osservazioni minori                                                       #
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize("parametro", [
    "run.batch_size", "run.seed", "run.threads", "qc.head_reads", "decontam.min_blanks",
    "filter.truncLen", "filter.trimLeft", "tax.min_boot", "retry.max_attempts",
    "ctrl.min_positives", "prev.min_count",
])
@pytest.mark.parametrize("valore", [True, False])
def test_un_booleano_non_e_un_intero_nella_configurazione(tmp_path, parametro, valore):
    """
    **Obiettivo**: Verificare che ``true`` e ``false`` dati a un parametro
    intero siano respinti dallo schema, nominando il parametro, e non letti
    come 1 e 0; e che il valore intero corrente dello stesso parametro resti
    accettato.

    **Razionale scientifico e sistemistico**: YAML legge ``true`` come
    booleano, che per Python e' l'intero 1: un refuso diventerebbe un lotto di
    un campione o zero controlli richiesti, senza alcun messaggio.
    """
    scenario = crea_scenario(tmp_path, _campioni())
    dati = scenario.config.model_dump(mode="python")
    gruppo, nome = parametro.split(".")
    assert valida(dati) is not None
    dati[gruppo][nome] = valore
    with pytest.raises(ErroreConfigurazione) as info:
        valida(dati)
    assert parametro in str(info.value)


def test_g02_non_attribuisce_alla_tabella_di_studio_il_guasto_dell_assay(tmp_path):
    """
    **Obiettivo**: Verificare che, con la tabella di assay che non si legge,
    G02 dichiari il guasto per ``io.assay_table`` e nessuna violazione lo
    attribuisca a ``io.study_table``; e che una tabulazione dentro
    l'identificativo di una riga di studio di un campione dell'assay fermi
    G02 con ``E-S0-02``, nominando la riga e la colonna.

    **Razionale scientifico e sistemistico**: Un messaggio che indica la
    tabella sbagliata, o un join che dice soltanto «manca la riga», manda chi
    corregge i metadati a cercare nel posto sbagliato.
    """
    scenario = crea_scenario(tmp_path / "assay", _campioni())
    Path(scenario.config.io.assay_table).write_bytes(b"Sample Name\tRaw Data File\n\xff\xfe\x00\n")
    esito = esegui_gate("G02", Contesto(scenario.config))
    assert not esito.superato
    testi = [str(v) for v in esito.violazioni]
    assert any("io.assay_table" in t for t in testi)
    assert not any("io.study_table" in t for t in testi), testi

    scenario = crea_scenario(tmp_path / "studio", _campioni())
    studio = Path(scenario.config.io.study_table)
    studio.write_text(
        studio.read_text(encoding="utf-8").replace("NOD1D4.L2\t", '"NOD1D4\t.L2"\t'),
        encoding="utf-8")
    esito = esegui_gate("G02", Contesto(scenario.config))
    assert [v.codice for v in esito.violazioni] == ["E-S0-02"]
    assert "riga 3" in str(esito.violazioni[0]) and "tabulazione" in str(esito.violazioni[0])


@pytest.mark.parametrize("dove", ["intestazione", "dati"])
def test_un_file_del_lotto_non_utf8_ferma_g08_con_un_codice_del_catalogo(tmp_path, dove):
    """
    **Obiettivo**: Verificare che un file del lotto con un byte non UTF-8,
    nell'intestazione o in una riga di dati, non produca un'eccezione di
    decodifica: i gate che precedono G08 passano, G08 si ferma con
    ``E-S0-08`` e il messaggio dice che le tabelle devono essere UTF-8.

    **Razionale scientifico e sistemistico**: Un file salvato da un foglio di
    calcolo in un'altra codifica e' un difetto comune dei metadati: va
    dichiarato come tale, con il rimedio, non come errore imprevisto.
    """
    scenario = crea_scenario(tmp_path, _campioni(), con_arricchimento=True, con_letture=True)
    lotto = Path(scenario.config.io.batch_table)
    contenuto = lotto.read_bytes()
    lotto.write_bytes(
        contenuto.replace(scenario.config.err.batch_column.encode(), b"c\xf2rsa", 1)
        if dove == "intestazione"
        else contenuto.replace(b"corsa_A", b"cors\xe0_A", 1))
    esiti = {e.gate: e for e in esegui_tutti(Contesto(scenario.config))}
    falliti = [g for g, e in esiti.items() if e.eseguito and not e.superato]
    assert falliti == ["G08"]
    (violazione,) = esiti["G08"].violazioni
    assert violazione.codice == "E-S0-08" and "UTF-8" in str(violazione)
    assert lotto.name in str(violazione)


def test_s14_si_ferma_se_l_oggetto_filtrato_manca(bioc, finale_calcolata, tmp_path):
    """
    **Obiettivo**: Verificare che, tolto ``ps_filtrato.rds`` dagli intermedi
    di S13, il calcolo di S14 si fermi con ``E-S14-01`` nominando il file,
    senza toccare i file consegnati.

    **Razionale scientifico e sistemistico**: L'oggetto finale si costruisce
    da quello filtrato: senza l'ingresso registrato da S13 non c'e' nulla da
    consegnare, e l'arresto deve avere un codice del catalogo.
    """
    run = copia_esecuzione(finale_calcolata, tmp_path)
    finale = run.albero.cartella(Fase.FINAL)
    prima = (finale / "checksum.sha256").read_bytes()
    (run.albero.cartella(Fase.FINAL_INTERMEDI) / NOME_FILTRATO).unlink()
    fase = run.fase(Passo.S14)
    contesto = run.contesto(Passo.S14, run.valuta()).ristretto(fase.parametri)
    with pytest.raises(ErrorePipeline) as info:
        fase.calcola(contesto)
    assert info.value.codice == "E-S14-01" and NOME_FILTRATO in info.value.dettaglio
    assert (finale / "checksum.sha256").read_bytes() == prima


def test_due_avvisi_con_lo_stesso_codice_sono_righe_distinte_sulla_console():
    """
    **Obiettivo**: Verificare che la riga della console di un avviso con un
    codice del catalogo riporti anche il dettaglio, cosi' che due avvisi con
    lo stesso codice e dettagli diversi diano due righe diverse; che un
    messaggio senza codice resti invariato; e che un errore non ripeta il
    dettaglio, che ha gia' nel riepilogo dell'arresto.

    **Razionale scientifico e sistemistico**: Due righe identiche sulla
    console sembrano un avviso ripetuto: chi esegue deve leggere che mancano
    sia i controlli positivi sia i negativi.
    """
    formattatore = _FormattatoreConsole()

    def riga(livello: int, messaggio: str, **extra) -> str:
        record = logging.LogRecord("amplicon16s", livello, __file__, 0, messaggio, (), None)
        record.__dict__.update(extra)
        return formattatore.format(record)

    sintesi = CATALOGO["E-S0-17"].sintesi
    prima = riga(logging.WARNING, sintesi, codice="E-S0-17", dettaglio="0 controlli negativi")
    seconda = riga(logging.WARNING, sintesi, codice="E-S0-17", dettaglio="nessun controllo positivo")
    assert prima != seconda
    assert prima == f"WARNING  [E-S0-17] {sintesi} 0 controlli negativi"
    assert riga(logging.INFO, "S0 completata") == "INFO     S0 completata"
    assert riga(logging.ERROR, sintesi, codice="E-S0-17", dettaglio="x").endswith(sintesi)


def test_i_file_scritti_in_modo_atomico_hanno_i_permessi_della_umask(tmp_path):
    """
    **Obiettivo**: Verificare che un file scritto con ``scrivi_atomico``, la
    richiesta del ponte verso R e un file forward ricavato dallo script di
    scarico del secondo dataset abbiano i permessi di un file creato
    normalmente dall'utente (666 meno la umask), non 600.

    **Razionale scientifico e sistemistico**: Un artefatto leggibile dal solo
    proprietario non si legge da un container avviato con un altro utente ne'
    da un collega dello stesso gruppo: i permessi li decide la umask di chi
    esegue, non la funzione che crea il file temporaneo.
    """
    maschera = os.umask(0)
    os.umask(maschera)
    attesi = 0o666 & ~maschera
    assert attesi != 0o600 or maschera == 0o177

    def permessi(percorso: Path) -> int:
        return stat.S_IMODE(percorso.stat().st_mode)

    scrivi_atomico(tmp_path / "a.txt", "testo")
    assert permessi(tmp_path / "a.txt") == attesi

    from amplicon16s.rbridge.payload import scrivi_richiesta

    richiesta, _ = scrivi_richiesta(tmp_path, "0" * 32, {"a": 1}, "S1")
    assert permessi(richiesta) == attesi

    percorso = RADICE / "dati" / "osd276" / "scarica_letture.py"
    if not percorso.is_file():
        pytest.skip("la cartella dati/ non e' presente (nell'immagine non viene copiata)")
    specifica = importlib.util.spec_from_file_location("scarica_letture_osd276", percorso)
    modulo = importlib.util.module_from_spec(specifica)
    import sys

    sys.path.insert(0, str(percorso.parents[1]))
    try:
        specifica.loader.exec_module(modulo)
    finally:
        sys.path.pop(0)
    deposito = _scrivi(tmp_path / "deposito.fastq.gz",
                       "@r1/1\nACGT\n+\nIIII\n@r1/2\nTTTT\n+\nIIII\n")
    assert modulo.ricava_forward(deposito, tmp_path / "forward.fastq.gz") == 1
    assert permessi(tmp_path / "forward.fastq.gz") == attesi


# --------------------------------------------------------------------------- #
# 4. Dalla revisione completa                                                  #
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize(("gruppo", "nome", "valore", "nel_messaggio"), [
    ("decontam", "method", "frequency", "decontam.method"),
    ("out", "taxa_are_rows", False, "taxa_are_rows"),
    ("filter", "maxEE", float("inf"), "filter.maxEE"),
    ("err", "nbases", float("inf"), "err.nbases"),
    ("io", "fastq_glob", "*/*.fastq.gz", "separatori di percorso"),
    ("io", "fastq_glob", "/dati/*.fastq.gz", "separatori di percorso"),
    ("run", "seed", 2**40, "run.seed"),
])
def test_lo_schema_respinge_cio_che_le_fasi_non_realizzano(tmp_path, gruppo, nome, valore,
                                                           nel_messaggio):
    """
    **Obiettivo**: Verificare che lo schema respinga, nominandoli, i valori che
    una fase a valle non realizza o non puo' ricevere: un metodo di
    decontaminazione diverso dalla prevalenza, le varianti sulle colonne, un
    numero infinito, un modello dei file con sottocartelle o assoluto, un seme
    fuori dagli interi di R.

    **Razionale scientifico e sistemistico**: Ognuno di questi valori era
    accettato all'avvio e fermava la catena piu' avanti, anche dopo ore di
    calcolo, con un errore fuori catalogo: cio' che non e' realizzato si
    respinge prima di cominciare, con il nome del parametro.
    """
    dati = crea_scenario(tmp_path, _campioni()).config.model_dump(mode="python")
    dati[gruppo][nome] = valore
    with pytest.raises(ErroreConfigurazione) as info:
        valida(dati)
    assert nel_messaggio in str(info.value)


def test_la_riga_di_comando_rifiuta_con_il_codice_3_gli_errori_di_chi_configura(tmp_path, capsys):
    """
    **Obiettivo**: Verificare che escano con il codice 3, e non con l'errore
    imprevisto, una configurazione con un gruppo scritto due volte (con il
    nome della chiave ripetuta), una non UTF-8, e una con ``io.out_root`` che
    e' un file o sta in una cartella non scrivibile; e che un percorso con la
    tilde sia espanso nella cartella personale, senza creare una cartella di
    nome ``~``.

    **Razionale scientifico e sistemistico**: Sono errori di chi configura,
    non difetti del programma: il codice di uscita e il messaggio devono dirlo.
    Una chiave ripetuta, in particolare, farebbe tornare ai predefiniti i
    parametri del primo blocco senza alcun avviso.
    """
    import yaml
    from conftest import dichiarazione_minima

    from amplicon16s.cli import _percorsi_assoluti, main
    from amplicon16s.config.schema import carica

    scenario = crea_scenario(tmp_path / "s", _campioni(), con_letture=True)
    dati = dichiarazione_minima(scenario.config)
    valido = yaml.safe_dump(dati, sort_keys=False)

    ripetuta = tmp_path / "ripetuta.yaml"
    ripetuta.write_text(valido + "qc:\n  head_reads: 50\n", encoding="utf-8")
    assert main(["validate", "--config", str(ripetuta)]) == 3
    assert "'qc' compare due volte" in capsys.readouterr().out

    latina = tmp_path / "latina.yaml"
    latina.write_bytes(valido.encode("utf-8") + "# qualit\xe0\n".encode("latin-1"))
    assert main(["validate", "--config", str(latina)]) == 3
    assert "UTF-8" in capsys.readouterr().out

    occupata = tmp_path / "file_al_posto_della_cartella"
    occupata.write_text("x", encoding="utf-8")
    sola_lettura = tmp_path / "sola_lettura"
    sola_lettura.mkdir()
    sola_lettura.chmod(0o500)
    for nome, radice in (("file", occupata), ("chiusa", sola_lettura / "out")):
        percorso = tmp_path / f"{nome}.yaml"
        percorso.write_text(
            yaml.safe_dump({**dati, "io": {**dati["io"], "out_root": str(radice)}}),
            encoding="utf-8")
        if os.access(sola_lettura, os.W_OK):  # da amministratore i permessi non fermano
            continue
        assert main(["validate", "--config", str(percorso)]) == 3, nome
        assert "io.out_root" in capsys.readouterr().out
    sola_lettura.chmod(0o700)

    tilde = tmp_path / "tilde.yaml"
    tilde.write_text(
        yaml.safe_dump({**dati, "io": {**dati["io"], "out_root": "~/uscita_di_prova"}}),
        encoding="utf-8")
    assoluta = _percorsi_assoluti(carica(tilde))
    assert Path(assoluta.io.out_root) == Path.home() / "uscita_di_prova"


def test_g07_respinge_il_numero_due_solo_se_c_e_il_compagno(tmp_path):
    """
    **Obiettivo**: Verificare che un dataset single-end con i file numerati
    (``..._1``, ``..._2``, ``..._3``: campioni, non letture) superi G07 quando
    il file con il 2 non ha accanto lo stesso nome con l'1; che lo stesso nome
    con ``_1`` e ``_2`` sia respinto; e che ``R2`` sia respinto anche da solo.

    **Razionale scientifico e sistemistico**: Il numero in fondo al nome e'
    anche il modo piu' comune di numerare campioni e repliche: respingere il
    solo file numero 2 fermerebbe un dataset legittimo senza un rimedio che
    non sia rinominare i file.
    """
    def esito(nomi: list[str]):
        campioni = [Campione(f"ERX500000{i}", f"c{i}", file=nome)
                    for i, nome in enumerate(nomi, start=1)]
        scenario = crea_scenario(tmp_path / str(abs(hash(tuple(nomi)))), campioni, con_letture=True)
        return esegui_gate("G07", Contesto(scenario.config))

    assert esito(["ERX5000001_a_3.fastq.gz", "ERX5000002_b_2.fastq.gz",
                  "ERX5000003_c.2.fastq.gz"]).superato
    con_compagno = esito(["ERX5000001_x_1.fastq.gz", "ERX5000002_y_R2.fastq.gz"])
    assert [v.codice for v in con_compagno.violazioni] == ["E-S0-07"]
    assert "ERX5000002_y_R2.fastq.gz" in str(con_compagno.violazioni[0])
    from amplicon16s.gates.g01_g15 import _e_lettura_inversa

    presenti = frozenset({"s_1.fastq.gz", "s_2.fastq.gz", "t.2.fq"})
    assert _e_lettura_inversa("s_2.fastq.gz", presenti)
    assert not _e_lettura_inversa("t.2.fq", presenti)
    assert not _e_lettura_inversa("s_1.fastq.gz", presenti)
    assert _e_lettura_inversa("u.R2_001.fastq.gz", frozenset())


def test_le_righe_di_sole_tabulazioni_non_sono_campioni(tmp_path):
    """
    **Obiettivo**: Verificare che le righe fatte di sole tabulazioni in fondo
    alla tabella di assay e al file del lotto siano ignorate dal lettore delle
    tabelle, che i gate sui metadati passino, e che l'inventario abbia i soli
    campioni veri; e che la numerazione delle righe nei messaggi resti quella
    del file anche dopo una riga vuota.

    **Razionale scientifico e sistemistico**: Un foglio di calcolo esporta
    spesso righe finali vuote: fermarsi con «valore vuoto» senza dire quale
    riga manda a cercare un campione che non esiste.
    """
    from conftest import esegui_gate_metadati

    from amplicon16s.metadata.tabelle import leggi_tsv, valori_non_tabellari

    scenario = crea_scenario(tmp_path, _campioni(), con_arricchimento=True, con_letture=True)
    for tabella in (scenario.config.io.assay_table, scenario.config.io.batch_table):
        with open(tabella, "a", encoding="utf-8") as file:
            file.write("\t\t\n\t\t\n")
    assert len(leggi_tsv(Path(scenario.config.io.assay_table))) == len(_campioni())
    assert len(esegui_gate_metadati(scenario.config)) == len(_campioni())

    prova = tmp_path / "prova.tsv"
    prova.write_text('a\tb\n\t\nx\t"y\tz"\n', encoding="utf-8")
    assert valori_non_tabellari(prova, ["b"]) == [(3, "b")]


def test_un_primer_a_posizione_variabile_non_si_cura_con_un_taglio_fisso(tmp_path):
    """
    **Obiettivo**: Verificare che, con il primer in testa a tutte le letture,
    ``E-S0-10`` indichi ``filter.trimLeft`` pari alla lunghezza del primer; e
    che con il primer in testa solo a una parte delle letture (preceduto
    altrove da basi di lunghezza variabile) il messaggio dica che un taglio
    fisso non lo toglie e che va rimosso per sequenza, senza prescrivere
    ``trimLeft``.

    **Razionale scientifico e sistemistico**: Con i distanziatori di
    eterogeneita' il rimedio indicato portava a un secondo arresto
    (``E-S0-16``) senza uscita: il messaggio deve distinguere i due casi.
    """
    from conftest import INIZIO_CON_PRIMER, scrivi_fastq

    scenario = crea_scenario(tmp_path / "testa", [
        Campione(f"ERX600000{i}", f"c{i}", inizio_letture=INIZIO_CON_PRIMER) for i in range(1, 4)
    ], con_letture=True)
    (violazione,) = esegui_gate("G10", Contesto(scenario.config)).violazioni
    assert violazione.codice == "E-S0-10"
    assert f"filter.trimLeft a {len(scenario.config.qc.primer_sequence)}" in str(violazione)

    scenario = crea_scenario(tmp_path / "variabile", [
        Campione(f"ERX600000{i}", f"c{i}") for i in range(1, 4)], con_letture=True)
    for percorso in Path(scenario.config.io.fastq_dir).glob("*.fastq.gz"):
        scrivi_fastq(percorso, [lettura("ACGTACG"[:n] + INIZIO_CON_PRIMER) for n in range(8)] * 5)
    (violazione,) = esegui_gate("G10", Contesto(scenario.config)).violazioni
    assert violazione.codice == "E-S0-10"
    assert "taglio fisso" in str(violazione) and "Imposta filter.trimLeft" not in str(violazione)


def test_un_successo_con_figli_senza_risultato_non_e_credibile():
    """
    **Obiettivo**: Verificare che il ponte classifichi come memoria esaurita
    anche uno script che dichiara il successo ed esce con 0, se sull'uscita di
    errore c'e' l'avviso dei figli del calcolo parallelo morti senza
    risultato; e che lo stesso successo senza l'avviso resti riuscito.

    **Razionale scientifico e sistemistico**: Alcune funzioni proseguono con i
    risultati dei figli sopravvissuti: lo script si conclude e l'artefatto e'
    sbagliato (una rimozione delle chimere su una parte dei campioni). Un
    risultato parziale dichiarato riuscito e' peggio di un arresto.
    """
    riuscito = Dichiarazione(stato=Stato.RIUSCITO, codice=None, messaggio="", artefatti=("a.rds",))
    avviso = "Warning in mclapply(...) :\n  2 parallel function calls did not deliver results\n"
    assert classifica(0, riuscito, avviso) is Condizione.MEMORIA_ESAURITA
    assert classifica(0, riuscito, "Warning: altro\n") is Condizione.RIUSCITO


def test_s2_filtra_un_dataset_di_un_solo_campione(bioc, tmp_path):
    """
    **Obiettivo**: Verificare che S2 si concluda su un dataset di un solo
    campione, con il suo conteggio delle letture nelle due tabelle del
    tracciamento sotto il nome del campione.

    **Razionale scientifico e sistemistico**: Con una sola riga l'estrazione
    di una colonna in R perde il nome del campione, e la fase si fermava con
    un errore generico: un dataset minimo e' un dataset legittimo.
    """
    scenario = crea_scenario(tmp_path, [Campione("ERX7000001", "solo")], con_letture=True)
    run = ProjectRun(scenario.config)
    esito = Esecutore(run, fino_a=Passo.S2).esegui()
    assert esito.conclusione is Conclusione.COMPLETATA, esito.punto and esito.punto.testo()
    for nome in ("letture_prefiltro.tsv", "letture_filtrate.tsv"):
        assert set(leggi_conteggi(run.albero.cartella(Fase.FILTERED) / nome)) == {"ERX7000001"}
