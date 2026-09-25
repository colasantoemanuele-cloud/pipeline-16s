"""Suite di verifica per il gestore dell'albero di output e dei manifesti SHA-256.

Inquadramento nel Piano Operativo:
    - **Settimana di riferimento**: **Settimana 5 (W5 — Fase F1: Catalogo degli errori,
      gestione degli artefatti e logging)**, con estensioni utilizzate dal meccanismo
      di ripresa su filesystem della **Settimana 9 (W9 — Fase F3)**.
    - **Scopo del modulo**: Verifica che la creazione delle 14 directory canoniche
      di output (``00_config`` .. ``99_logs``), il calcolo streaming dei digest
      crittografici SHA-256 e l'aggiornamento atomico dei manifesti JSON garantiscano
      il rilevamento deterministico di qualsiasi cancellazione, troncamento o
      manomissione anche di un singolo byte sugli artefatti intermedi della pipeline.
    - **Moduli sorgente coperti**:
        * ``src/amplicon16s/io_layer/checksums.py``
        * ``src/amplicon16s/io_layer/artifacts.py``
"""

from __future__ import annotations

import json

import pytest

from amplicon16s.io_layer.artifacts import NOME_MANIFESTO, AlberoOutput, Fase
from amplicon16s.io_layer.checksums import (
    ALGORITMO,
    checksum_bytes,
    checksum_file,
    corrisponde,
)

CARTELLE_ATTESE = (
    "00_config", "01_input_validation", "02_qc_profiles", "03_filtered",
    "04_error_models", "05_asv_inference", "06_seqtab", "07_chimera",
    "08_taxonomy", "09_phylogeny", "10_phyloseq", "11_controls", "12_final",
    "99_logs",
)


@pytest.fixture
def albero(tmp_path) -> AlberoOutput:
    """Fornisce un'istanza isolata di ``AlberoOutput`` ancorata a ``tmp_path``."""
    return AlberoOutput(tmp_path)


# --------------------------------------------------------------------------- #
# Checksum                                                                     #
# --------------------------------------------------------------------------- #


def test_il_checksum_dichiara_il_proprio_algoritmo(tmp_path):
    """
    **Obiettivo**: Verificare che l'impronta restituita da ``checksum_file`` sia
    prefissata con l'identificativo dell'algoritmo (``sha256:<esadecimale>``).

    **Razionale Scientifico/Sistemistico**: Esplicitare il prefisso ``sha256:``
    nei manifesti previene ambiguità crittografiche (es. confusione con MD5 usato
    solo per il database FASTA esterno in G12) e garantisce l'auto-descrittività
    a lungo termine degli artefatti depositati a fine studio.
    """
    file = tmp_path / "x.txt"
    file.write_text("contenuto")
    assert checksum_file(file).startswith(f"{ALGORITMO}:")


def test_checksum_di_file_e_di_byte_coincidono(tmp_path):
    """
    **Obiettivo**: Verificare l'identità matematica tra l'hash calcolato su un
    buffer in memoria (``checksum_bytes``) e quello calcolato leggendo il file da
    disco (``checksum_file``).

    **Razionale Scientifico/Sistemistico**: Quando ``AlberoOutput.scrivi_testo``
    o ``scrivi_bytes`` persiste un artefatto, calcola il checksum dai byte in
    memoria, mentre alla ripresa (`resume`) ``ProjectRun`` ricalcola il checksum
    leggendo il file dal disco: se le due funzioni divergessero, ogni ripresa
    invaliderebbe falsamente tutte le fasi completate.
    """
    dati = b"contenuto qualunque"
    file = tmp_path / "x.bin"
    file.write_bytes(dati)
    assert checksum_file(file) == checksum_bytes(dati)


def test_contenuti_diversi_hanno_checksum_diversi():
    """
    **Obiettivo**: Verificare la sensibilità della funzione di hashing alla
    variazione di un singolo byte del payload.

    **Razionale Scientifico/Sistemistico**: Assicura l'assenza di banalità
    nell'implementazione del wrapper SHA-256, condizione necessaria affinché
    anche una minima modifica ad una tabella di abbondanza ASV o tassonomica
    venga intercettata dal controllo di integrità.
    """
    assert checksum_bytes(b"a") != checksum_bytes(b"b")


def test_un_file_grande_non_viene_caricato_in_memoria(tmp_path):
    """
    **Obiettivo**: Verificare la correttezza del calcolo SHA-256 quando la
    dimensione del file (> 3 MiB) supera i blocchi di lettura (1 MiB) di ``checksum_file``.

    **Razionale Scientifico/Sistemistico**: Gli artefatti delle fasi S2 (FASTQ
    filtrati per 960 campioni) e S4/S10 (oggetti ``.rds`` DADA2 e ``phyloseq``)
    possono raggiungere centinaia di megabyte o gigabyte. La lettura chunked a
    blocchi evita picchi di allocazione RAM (OOM) nell'orchestratore Python e
    deve concatenare correttamente tutti i blocchi parziali nello stato SHA-256.
    """
    dati = b"x" * (3 * 1024 * 1024 + 7)
    file = tmp_path / "grande.bin"
    file.write_bytes(dati)
    assert checksum_file(file) == checksum_bytes(dati)


def test_un_file_assente_non_corrisponde(tmp_path):
    """
    **Obiettivo**: Verificare che ``corrisponde()`` restituisca ``False`` (senza
    sollevare ``FileNotFoundError``) quando il file indicato non esiste su disco.

    **Razionale Scientifico/Sistemistico**: Nel modello di ripresa basato sul
    filesystem, l'assenza di un artefatto non è un crash del programma ma il
    segnale fisiologico che indica che la fase è incompleta e deve essere
    pianificata per l'esecuzione.
    """
    assert corrisponde(tmp_path / "mai_scritto", checksum_bytes(b"")) is False


# --------------------------------------------------------------------------- #
# Albero delle cartelle                                                        #
# --------------------------------------------------------------------------- #


def test_le_fasi_sono_quattordici():
    """
    **Obiettivo**: Verificare che l'enumerazione ``Fase`` contenga esattamente
    le 14 directory di output previste dall'architettura.

    **Razionale Scientifico/Sistemistico**: Blocca regressioni strutturali
    nell'albero di output: mentre i passi logici sono 15 (`S0..S14`), le cartelle
    fisiche sono 14 perché `07_chimera` (S6, S7), `11_controls` (S11, S12) e
    `12_final` (S13, S14) aggregano fasi contigue più `00_config` e `99_logs`.
    """
    assert len(list(Fase)) == 14


def test_i_nomi_delle_cartelle_sono_quelli_previsti():
    """
    **Obiettivo**: Verificare che i nomi delle 14 directory corrispondano
    esattamente alla nomenclatura numerata (``00_config`` .. ``12_final``, ``99_logs``).

    **Razionale Scientifico/Sistemistico**: Gli script R a valle, i moduli di
    analisi ecologica (`amplicon16s_eco`) e il generatore di report si aspettano
    percorsi deterministici e ordinati lessicograficamente secondo l'ordine del
    flusso bioinformatico.
    """
    assert tuple(f.value for f in Fase) == CARTELLE_ATTESE


def test_creare_l_albero_produce_tutte_le_cartelle(albero, tmp_path):
    """
    **Obiettivo**: Verificare che ``albero.crea()`` materializzi su disco tutte
    le 14 sottocartelle sotto ``io.out_root``.

    **Razionale Scientifico/Sistemistico**: Assicura che fin dall'avvio della
    pipeline tutte le destinazioni per i manifesti, gli artefatti R e il file
    di log strutturato ``99_logs/pipeline.jsonl`` siano pronte e scrivibili.
    """
    albero.crea()
    presenti = sorted(p.name for p in tmp_path.iterdir() if p.is_dir())
    assert presenti == sorted(CARTELLE_ATTESE)


def test_creare_l_albero_e_ripetibile(albero, tmp_path):
    """
    **Obiettivo**: Verificare che invocazioni successive di ``albero.crea()``
    siano idempotenti e non sovrascrivano né cancellino i file già presenti.

    **Razionale Scientifico/Sistemistico**: Durante un comando ``amplicon16s resume``,
    l'inizializzazione dell'albero viene rieseguita su una directory di output
    popolata da ore di calcolo DADA2: una creazione distruttiva azzererebbe il
    lavoro già completato.
    """
    albero.crea()
    albero.scrivi_testo(Fase.FINAL, "nota.txt", "contenuto")
    albero.crea()
    assert (tmp_path / Fase.FINAL.value / "nota.txt").is_file()


def test_costruire_l_albero_non_tocca_il_filesystem(tmp_path):
    """
    **Obiettivo**: Verificare che la sola istanziazione di ``AlberoOutput(radice)``
    non crei directory su disco prima di una scrittura esplicita.

    **Razionale Scientifico/Sistemistico**: Evita effetti collaterali indesiderati
    quando si ispeziona o si valida una configurazione in sola lettura senza
    voler ancora sporcare il filesystem.
    """
    radice = tmp_path / "non_creata"
    AlberoOutput(radice)
    assert not radice.exists()


def test_la_cartella_nasce_alla_prima_scrittura(albero, tmp_path):
    """
    **Obiettivo**: Verificare che ``scrivi_testo`` crei automaticamente (lazy
    creation) la directory della fase bersaglio qualora non esista ancora.

    **Razionale Scientifico/Sistemistico**: Rende ogni singola fase autonoma e
    resiliente anche se invocata isolatamente o se una cartella vuota è stata
    rimossa manualmente dall'operatore prima di una ripresa.
    """
    assert not (tmp_path / Fase.SEQTAB.value).exists()
    albero.scrivi_testo(Fase.SEQTAB, "x.txt", "y")
    assert (tmp_path / Fase.SEQTAB.value).is_dir()


# --------------------------------------------------------------------------- #
# Manifesto                                                                    #
# --------------------------------------------------------------------------- #


def test_l_artefatto_scritto_finisce_nel_manifesto(albero):
    """
    **Obiettivo**: Verificare che ogni scrittura tramite ``AlberoOutput`` registri
    immediatamente nel ``manifest.json`` il nome, il checksum SHA-256 e la
    dimensione in byte del file.

    **Razionale Scientifico/Sistemistico**: Garantisce l'accoppiamento atomico
    tra dato e metadato di integrità: nessun artefatto può nascere senza che la
    sua impronta crittografica venga contestualmente blindata nel manifesto.
    """
    artefatto = albero.scrivi_testo(Fase.TAXONOMY, "taxa.tsv", "a\tb\n")

    manifesto = albero.manifesto(Fase.TAXONOMY)
    assert manifesto["taxa.tsv"]["checksum"] == artefatto.checksum
    assert manifesto["taxa.tsv"]["byte"] == len("a\tb\n")


def test_il_manifesto_e_interrogabile(albero):
    """
    **Obiettivo**: Verificare che il file ``manifest.json`` scritto su disco
    rispetti lo schema JSON atteso (chiavi ``fase`` e lista ``artefatti``).

    **Razionale Scientifico/Sistemistico**: Consente sia a ``ProjectRun`` sia
    agli strumenti esterni di audit e reportistica di ispezionare i manifesti
    di fase come documenti JSON standard e auto-contenuti.
    """
    albero.scrivi_testo(Fase.QC_PROFILES, "profili.tsv", "x")
    documento = json.loads(
        (albero.cartella(Fase.QC_PROFILES) / NOME_MANIFESTO).read_text(encoding="utf-8")
    )
    assert documento["fase"] == "02_qc_profiles"
    assert [v["nome"] for v in documento["artefatti"]] == ["profili.tsv"]


def test_il_manifesto_non_dipende_dall_ordine_di_scrittura(tmp_path):
    """
    **Obiettivo**: Verificare che la serializzazione di ``manifest.json`` ordini
    lessicograficamente le voci degli artefatti producendo un file identico bit
    a bit indipendentemente dall'ordine temporale di scrittura dei file.

    **Razionale Scientifico/Sistemistico**: Garantisce che l'hash SHA-256 del
    manifesto stesso sia strettamente deterministico e riproducibile (idempotenza).
    Se l'impronta del manifesto dipendesse dall'ordine non deterministico con cui
    thread o sottoprocessi R chiudono i file, le fasi a valle rileverebbero una
    falsa variazione dell'impronta ``a_monte`` e ricalcolerebbero inutilmente
    l'intera pipeline, sprecando ore di calcolo computazionale.
    """
    uno, due = AlberoOutput(tmp_path / "a"), AlberoOutput(tmp_path / "b")
    for nome in ("z.txt", "a.txt", "m.txt"):
        uno.scrivi_testo(Fase.FILTERED, nome, nome)
    for nome in ("m.txt", "z.txt", "a.txt"):
        due.scrivi_testo(Fase.FILTERED, nome, nome)

    assert (uno.cartella(Fase.FILTERED) / NOME_MANIFESTO).read_text() == (
        due.cartella(Fase.FILTERED) / NOME_MANIFESTO
    ).read_text()


def test_riscrivere_un_artefatto_aggiorna_la_sua_voce(albero):
    """
    **Obiettivo**: Verificare che la sovrascrittura di un artefatto esistente
    aggiorni in-place la voce corrispondente nel manifesto senza creare duplicati.

    **Razionale Scientifico/Sistemistico**: Quando una fase viene rieseguita
    dopo una ripresa o un tentativo di retry automatico, i nuovi artefatti
    sostituiscono quelli del tentativo precedente: il manifesto deve riflettere
    esclusivamente lo stato dell'ultima esecuzione valida.
    """
    albero.scrivi_testo(Fase.SEQTAB, "tab.tsv", "primo")
    secondo = albero.scrivi_testo(Fase.SEQTAB, "tab.tsv", "secondo")

    manifesto = albero.manifesto(Fase.SEQTAB)
    assert len(manifesto) == 1
    assert manifesto["tab.tsv"]["checksum"] == secondo.checksum


def test_il_manifesto_di_una_fase_mai_eseguita_e_vuoto(albero):
    """
    **Obiettivo**: Verificare che ``albero.manifesto()`` restituisca un dizionario
    vuoto ``{}`` per una fase nella cui cartella non è ancora stato scritto nulla.

    **Razionale Scientifico/Sistemistico**: Permette di interrogare in sicurezza
    lo stato di qualsiasi fase (incluse fasi opzionali disattivate come S9
    filogenesi) senza sollevare eccezioni di I/O.
    """
    assert albero.manifesto(Fase.PHYLOGENY) == {}


# --------------------------------------------------------------------------- #
# Integrità                                                                    #
# --------------------------------------------------------------------------- #


def test_un_artefatto_alterato_non_e_integro(albero):
    """
    **Obiettivo**: Verificare che la modifica del contenuto di un artefatto già
    registrato faccia passare ``artefatto.integro`` e ``fase_completa`` a ``False``,
    elencando il file in ``non_integri()``.

    **Razionale Scientifico/Sistemistico**: Previene l'uso silenzioso di file
    corrotti o modificati accidentalmente su disco (es. un file ``.rds`` o TSV
    alterato manualmente dopo il calcolo), obbligando il grafo di ripresa a
    invalidare e ricalcolare la fase compromessa e tutte le sue dipendenti.
    """
    artefatto = albero.scrivi_testo(Fase.FINAL, "ps_final.txt", "contenuto originale")
    assert artefatto.integro

    artefatto.percorso.write_text("contenuto manomesso", encoding="utf-8")

    assert not artefatto.integro
    assert albero.non_integri(Fase.FINAL) == ("ps_final.txt",)
    assert not albero.fase_completa(Fase.FINAL)


def test_un_artefatto_rimosso_non_e_integro(albero):
    """
    **Obiettivo**: Verificare che la cancellazione fisica di un file registrato
    nel manifesto venga rilevata da ``non_integri()`` e invalidi ``fase_completa()``.

    **Razionale Scientifico/Sistemistico**: Garantisce che il manifesto da solo
    non basti a dichiarare conclusa una fase se il file fisico corrispondente è
    stato rimosso dal filesystem.
    """
    albero.scrivi_testo(Fase.CHIMERA, "chimere.tsv", "x")
    (albero.cartella(Fase.CHIMERA) / "chimere.tsv").unlink()

    assert albero.non_integri(Fase.CHIMERA) == ("chimere.tsv",)
    assert not albero.fase_completa(Fase.CHIMERA)


def test_un_troncamento_viene_rilevato(albero):
    """
    **Obiettivo**: Verificare che un artefatto troncato a metà rispetto alla
    sua lunghezza originale venga immediatamente segnalato come non integro.

    **Razionale Scientifico/Sistemistico**: Simula il caso classico di un
    processo interrotto da ``SIGKILL`` (OOM killer) o da saturazione dello
    spazio disco mentre stava scaricando un file su disco: il file esiste e non
    è vuoto, ma il confronto SHA-256 impedisce che una tabella troncata venga
    letta dalla fase successiva producendo risultati scientificamente falsati.
    """
    albero.scrivi_testo(Fase.ERROR_MODELS, "modello.txt", "a" * 1000)
    (albero.cartella(Fase.ERROR_MODELS) / "modello.txt").write_text("a" * 400)
    assert albero.non_integri(Fase.ERROR_MODELS) == ("modello.txt",)


def test_una_fase_con_artefatti_integri_risulta_completa(albero):
    """
    **Obiettivo**: Verificare che ``fase_completa()`` restituisca ``True`` e
    ``non_integri()`` restituisca ``()`` quando tutti i file registrati nel
    manifesto esistono e corrispondono ai rispettivi hash SHA-256.

    **Razionale Scientifico/Sistemistico**: Costituisce la condizione necessaria
    affinché ``ProjectRun`` possa saltare le fasi già concluse durante un
    ``resume`` senza ripetere calcoli ridondanti.
    """
    albero.scrivi_testo(Fase.CONTROLS, "katharoseq.tsv", "x")
    albero.scrivi_testo(Fase.CONTROLS, "positivi.tsv", "y")
    assert albero.fase_completa(Fase.CONTROLS)
    assert albero.non_integri(Fase.CONTROLS) == ()


def test_una_fase_mai_eseguita_non_risulta_completa(albero):
    """
    **Obiettivo**: Verificare che ``fase_completa()`` restituisca ``False`` per
    una fase priva di manifesto su disco.

    **Razionale Scientifico/Sistemistico**: Evita che l'assenza di file non
    integri (insieme vuoto su cartella mai eseguita) venga scambiata per
    completamento con successo.
    """
    assert not albero.fase_completa(Fase.ASV_INFERENCE)


def test_un_manifesto_parziale_non_significa_fase_conclusa(albero):
    """
    **Obiettivo**: Verificare che ``fase_completa(fase, attesi=[...])`` restituisca
    ``False`` se nel manifesto manca anche uno solo degli artefatti richiesti.

    **Razionale Scientifico/Sistemistico**: Se una fase deve produrre due file
    (es. ``seqtab.tsv`` e ``conteggi.tsv``) e si interrompe dopo aver scritto e
    registrato solo il primo, il primo file è integro ma la fase è incompleta:
    il controllo sull'insieme atteso impedisce riprese su output parziali.
    """
    albero.scrivi_testo(Fase.SEQTAB, "seqtab.tsv", "x")
    assert albero.fase_completa(Fase.SEQTAB)
    assert not albero.fase_completa(Fase.SEQTAB, attesi=["seqtab.tsv", "conteggi.tsv"])


# --------------------------------------------------------------------------- #
# Artefatti prodotti fuori da Python                                           #
# --------------------------------------------------------------------------- #


def test_un_file_scritto_da_altri_puo_essere_registrato(albero):
    """
    **Obiettivo**: Verificare che ``albero.registra(fase, nome)`` calcoli il
    checksum SHA-256 di un file già scritto su disco da un processo esterno e
    lo inserisca nel manifesto della fase.

    **Razionale Scientifico/Sistemistico**: Le fasi di calcolo statistico e
    bioinformatico (`S1..S13`) vengono eseguite da sottoprocessi ``Rscript``
    indipendenti che salvano su disco oggetti ``.rds`` e tabelle TSV. Al termine
    del processo R, il ponte Python ``rbridge`` utilizza ``albero.registra`` per
    sottoporre gli artefatti prodotti da R allo stesso regime di verifica SHA-256
    degli artefatti scritti direttamente in Python.
    """
    cartella = albero.prepara(Fase.PHYLOSEQ)
    (cartella / "ps.rds").write_bytes(b"finto contenuto rds")

    artefatto = albero.registra(Fase.PHYLOSEQ, "ps.rds")

    assert artefatto.integro
    assert albero.manifesto(Fase.PHYLOSEQ)["ps.rds"]["checksum"] == artefatto.checksum
    assert albero.fase_completa(Fase.PHYLOSEQ)


def test_registrare_un_file_inesistente_e_un_errore(albero):
    """
    **Obiettivo**: Verificare che ``albero.registra`` sollevi ``FileNotFoundError``
    se il file dichiarato non esiste fisicamente nella cartella della fase.

    **Razionale Scientifico/Sistemistico**: Impedisce che uno script R buggato
    possa dichiarare nel proprio JSON d'esito di aver prodotto un artefatto
    (es. ``ps.rds``) senza averlo realmente scritto su disco.
    """
    albero.prepara(Fase.PHYLOSEQ)
    with pytest.raises(FileNotFoundError, match="non trovato"):
        albero.registra(Fase.PHYLOSEQ, "mai_scritto.rds")


def test_artefatti_binari(albero):
    """
    **Obiettivo**: Verificare che ``scrivi_bytes`` persista flussi binari grezzi
    preservandone i byte esatti e registrandone il relativo checksum SHA-256.

    **Razionale Scientifico/Sistemistico**: Garantisce che file binari non
    testuali (come archivi compressi o strutture serializzate) non subiscano
    alterazioni di codifica UTF-8 o conversioni di fine riga.
    """
    dati = b"\x00\x01\x02 contenuto binario"
    artefatto = albero.scrivi_bytes(Fase.SEQTAB, "seqtab.bin", dati)
    assert artefatto.percorso.read_bytes() == dati
    assert artefatto.checksum == checksum_bytes(dati)
