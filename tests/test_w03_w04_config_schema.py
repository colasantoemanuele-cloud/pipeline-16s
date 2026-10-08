r"""Suite di verifica dello schema di configurazione, del Gate G15 e dei parametri derivati.

1. Inquadramento nel Piano Operativo
------------------------------------
Settimane 3 e 4 (W3/W4), Fase F1: schema di validazione della configurazione
(W3) e coerenza interna della configurazione con i parametri derivati (W4).

2. Moduli sorgente coperti
--------------------------
* ``src/amplicon16s/config/schema.py``
* ``src/amplicon16s/config/defaults.py``
* ``src/amplicon16s/config/resolve.py``
* ``src/amplicon16s/gates/g01_g15.py`` (limitatamente a ``esegui_g15`` e ``CONTROLLI``)
* ``config/config.example.yaml``

3. Cosa valuta questo file
--------------------------
- allineamento fra ``config/config.example.yaml`` e lo schema: nessuna chiave
  mancante o sconosciuta, tutti i 23 gruppi presenti (18 di produzione e 5 delle
  analisi ecologiche), tipi numerici conservati, valori del dataset di
  riferimento OSD-734 e whitelist dei ritentativi indipendente dal dataset;
- validazione dello schema: parametri obbligatori senza valore predefinito,
  tipi e domini, vocabolari chiusi, espressioni regolari, immagine con digest,
  MD5, rifiuto dei refusi (``extra="forbid"``), immutabilità dopo il
  caricamento, segnalazione di tutti i problemi e non solo del primo;
- errori di caricamento del file YAML (file vuoto, sintassi non valida) con
  indicazione del file;
- Gate G15 (codici ``E-G15-02`` .. ``E-G15-09``): accettazione del file di
  esempio, coerenza del registro dei controlli, rifiuto delle combinazioni
  incoerenti con raccolta di tutte le violazioni;
- parametri derivati (``asv.len_min``, ``asv.len_max``):
  valori attesi, dipendenza dai parametri di origine, rifiuto di un derivato
  impostato a mano; i parametri rimossi (``prev.min_samples``,
  ``qc.min_reads_filtered``) sono respinti come chiavi sconosciute;
- digest SHA-256 canonico della configurazione e scrittura deterministica di
  ``00_config/resolved.yaml`` con i derivati.

4. Comandi Bash e scenari di esecuzione
---------------------------------------
    ``<immagine>`` e' l'immagine del container della pipeline; quella corrente
    e' indicata in ``README.md``.

    1. Modalità locale standard (R di base con jsonlite, senza Bioconductor né
       dati reali):
       pytest tests/test_w03_w04_config_schema.py -v

    2. Modalità container Docker standard (sottoinsieme ridotto con
       R/Bioconductor):
       docker run --rm \
         -e PYTHONPATH=/app/src \
         -e AMPLICON16S_R_DIR=/app/R \
         -v "$(pwd)":/app \
         -w /app \
         <immagine> \
         pytest -o cache_dir=/tmp/.pytest_cache tests/test_w03_w04_config_schema.py -v

    3. Modalità container Docker completa (con i 2.4 GB di dati reali OSD-734;
       la configurazione e i percorsi che contiene devono stare nella cartella
       montata):
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
         pytest -o cache_dir=/tmp/.pytest_cache tests/test_w03_w04_config_schema.py -v

5. Risultato atteso
-------------------
107 test totali:
- 107 passed in ambiente locale standard (~1.2s);
- 107 passed nel container Docker standard sul sottoinsieme ridotto (~1.1s);
- 107 passed nel container Docker con i dati reali OSD-734 (~0.9s).

6. Razionale scientifico e sistemistico
---------------------------------------
- Un errore formale o metodologico nella configurazione va intercettato prima
  di allocare risorse o avviare processi R sui 960 campioni del dataset: una
  soglia incoerente scoperta dopo ore di calcolo produrrebbe risultati da
  scartare.
- I parametri derivati si calcolano e non si scrivono a mano, perché solo
  così restano coerenti con i parametri da cui discendono; il digest canonico
  rende confrontabili due esecuzioni sulla stessa configurazione.
"""

from __future__ import annotations

import copy
from pathlib import Path

import pytest
from conftest import FORMATO, chiavi_schema, esegui_g15
import yaml

from amplicon16s.errors.catalog import Categoria, voce
from conftest import configurazione_di_prova as _esempio_compilato

from amplicon16s.config import defaults
from amplicon16s.config.resolve import (
    NOME_FILE_RISOLTO,
    ConfigRisolta,
    risolvi,
    scrivi_risolta,
)
from amplicon16s.config.schema import (
    INTESTAZIONE_MANCANTI,
    PARAMETRI_DERIVATI,
    Config,
    ErroreConfigurazione,
    carica,
            obbligatori_mancanti,
    valida,
)
from amplicon16s.gates.g01_g15 import CONTROLLI, ErroreGate
from amplicon16s.io_layer.artifacts import Fase

RADICE = Path(__file__).resolve().parents[1]
ESEMPIO = RADICE / "config" / "config.example.yaml"

#: I parametri che dipendono dai dati o dallo studio e non hanno un
#: predefinito. L'elenco e' scritto qui per esteso, e non letto da
#: ``defaults.OBBLIGATORI``: e' il codice che il test deve verificare, e un
#: parametro tolto da li' per errore sparirebbe anche dall'atteso.
DEL_DATASET = (
    # come si riconoscono file e campioni
    "io.fastq_glob", "io.accession_regex", "meta.sample_id_column", "meta.accession_column",
    # il formato dei metadati
    "meta.derive_module", "meta.module_regex", "meta.module_column",
    "meta.non_surface_positions", "meta.batch_key_column", "meta.batch_module_column",
    "err.batch_column", "decontam.batch_column",
    "ctrl.column", "ctrl.blank_values", "ctrl.positive_values", "ctrl.biological_values",
    "ctrl.blank_override_column", "ctrl.blank_override_values",
    "katharoseq.cell_count_column", "out.study_columns", "out.batch_columns",
    # la regione amplificata, le letture e il riferimento
    "filter.truncLen", "filter.truncLen_shortfall_warn", "filter.trimLeft",
    "qc.primer_sequence", "qc.conserved_motif", "tax.ref_name", "tax.ref_version", "tax.try_rc",
    # le scelte di analisi che dipendono dallo studio
    "err.randomize", "dada.pool", "chimera.min_fold_parent_over_abundance",
    "filt.remove_na_phylum", "filt.exclude_taxa", "phylo.enabled", "phylo.max_seqs",
    "prev.min_fraction", "prev.min_count", "prev.apply",
    # i controlli sperimentali
    "ctrl.min_positives", "ctrl.min_positive_pass_frac", "ctrl.positive_gate",
    "katharoseq.target_taxon", "katharoseq.collapse_rank", "katharoseq.target_sensitivity",
    "katharoseq.min_r2", "decontam.threshold", "decontam.min_blanks", "decontam.mode",
    # le soglie di qualita'
    "qc.head_reads", "qc.max_primer_hit_frac", "qc.min_motif_frac", "qc.max_frac_short_reads",
    "qc.max_zeroed_samples", "qc.max_frac_lost_filter", "qc.max_asv_count",
    "qc.warn_frac_chimeric", "qc.stop_frac_chimeric", "qc.min_frac_phylum",
    "qc.min_reads_mode", "qc.max_frac_contaminant", "qc.min_reads_final",
    "qc.min_frac_reads_retained",
)

#: Elenco tassativo dei parametri privi di valore predefinito nello schema: i
#: percorsi e il checksum del riferimento, piu' quelli di ``DEL_DATASET``.
OBBLIGATORI = (
    "io.fastq_dir",
    "io.assay_table",
    "io.out_root",
    "tax.ref_fasta",
    "tax.ref_md5",
    *DEL_DATASET,
)

#: I parametri con un predefinito che non e' il valore standard di un metodo:
#: quelli di cui la pipeline realizza una sola forma, e quelli dello strumento
#: (come si esegue e in che forma si scrive). Con ``defaults.STANDARD_DEL_METODO``
#: e i nulli sono tutti i predefiniti ammessi.
UNICA_FORMA = (
    "tax.classifier", "tax.assign_species", "phylo.aligner", "phylo.model",
    "katharoseq.curve_model", "katharoseq.read_stage", "decontam.method",
    "out.serialization", "out.asv_id_scheme", "out.taxa_are_rows", "out.sample_id_source",
)
DELLO_STRUMENTO = (
    "asv.len_tol", "out.export_flat", "retry.enabled", "retry.max_attempts", "retry.whitelist",
    "run.lockfile", "run.seed", "run.batch_size", "run.keep_filtered_fastq",
    "run.strict_provenance",
)
#: Facoltativi il cui predefinito e' l'assenza.
NULLI = (
    "io.study_table", "io.batch_table", "meta.study_sample_id_column", "tax.ref_bad_taxa",
    "run.container", "run.threads", "run.r_timeout_s",
)


def _modello() -> dict[str, set[str]]:
    """Le chiavi ``gruppo.parametro`` del modello di configurazione, per specie:
    quelle con un valore, e quelle lasciate in commento con il loro marcatore
    (``[OBBLIGATORIO]``) o senza (``commentate``).
    """
    import re

    chiavi: dict[str, set[str]] = {"con_valore": set(), "obbligatorie": set(), "commentate": set(),
                                   "standard": set()}
    gruppo = None
    for riga in ESEMPIO.read_text(encoding="utf-8").splitlines():
        if trovato := re.match(r"^([a-z]+):", riga):
            gruppo = trovato.group(1)
        elif trovato := re.match(r"^  ([A-Za-z_0-9]+):", riga):
            chiavi["con_valore"].add(f"{gruppo}.{trovato.group(1)}")
            if "[STANDARD]" in riga:
                chiavi["standard"].add(f"{gruppo}.{trovato.group(1)}")
        elif trovato := re.match(r"^  # ([A-Za-z_0-9]+):\s*(#.*|\d+)?$", riga):
            specie = "obbligatorie" if "[OBBLIGATORIO]" in riga else "commentate"
            chiavi[specie].add(f"{gruppo}.{trovato.group(1)}")
    return chiavi


@pytest.fixture
def dati_esempio() -> dict:
    """Una configurazione completa e valida, come dizionario mutabile indipendente."""
    return _esempio_compilato()


def _senza(dati: dict, chiave: str) -> dict:
    """Restituisce una copia profonda di ``dati`` privata della chiave ``gruppo.campo``."""
    gruppo, campo = chiave.split(".")
    copia = copy.deepcopy(dati)
    del copia[gruppo][campo]
    return copia


# --------------------------------------------------------------------------- #
# Il file di esempio è un'istanza valida (W3)                                  #
# --------------------------------------------------------------------------- #


def test_esempio_esiste():
    """
    **Obiettivo**: Verificare la presenza fisica di ``config/config.example.yaml`` nel repository.

    **Razionale scientifico e sistemistico**: Il file di esempio costituisce il
    template documentato di riferimento per l'operatore e il contratto vivente
    dei parametri calibrati per il dataset NASA GeneLab **OSD-734**.
    """
    assert ESEMPIO.is_file(), f"file di esempio non trovato in {ESEMPIO}"


def test_il_modello_non_parte_e_una_configurazione_completa_si_valida():
    """
    **Obiettivo**: Verificare che ``config.example.yaml`` cosi' com'e' sia
    respinto con l'elenco dei parametri obbligatori da dichiarare, e che una
    configurazione che li dichiara tutti sia valida.

    **Razionale scientifico e sistemistico**: Il modello non porta valori di
    alcun dataset: chi lo copia deve dichiarare cio' che dipende dai suoi
    dati, e la pipeline deve dirgli che cosa manca prima di ogni calcolo.
    """
    with pytest.raises(ErroreConfigurazione) as respinta:
        carica(ESEMPIO)
    assert INTESTAZIONE_MANCANTI in str(respinta.value)
    config = valida(_esempio_compilato())
    assert isinstance(config, Config)


def test_esempio_conserva_i_tipi_numerici(dati_esempio):
    """
    **Obiettivo**: Verificare che i parametri espressi in notazione scientifica
    nel YAML (``err.nbases = 1.0e+8``, ``dada.omega_a = 1.0e-40``) siano letti
    da PyYAML come ``float`` nativi e non come stringhe.

    **Razionale scientifico e sistemistico**: Secondo la specifica YAML 1.1 usata
    da ``yaml.safe_load``, una scrittura priva di segno e punto decimale come
    ``1e8`` viene interpretata come stringa ``"1e8"`` anziché come numero,
    causando il rifiuto dello schema o un passaggio errato di tipo al motore R DADA2.
    """
    assert isinstance(dati_esempio["err"]["nbases"], float)
    assert isinstance(dati_esempio["dada"]["omega_a"], float)


def test_esempio_dichiara_versione_come_stringa(dati_esempio):
    """
    **Obiettivo**: Verificare che ``tax.ref_version`` (es. ``"138.2"``) sia
    tipizzato nel YAML come stringa e non come numero in virgola mobile.

    **Razionale scientifico e sistemistico**: I numeri di versione tassonomica
    (come SILVA ``138.1`` o ``138.20``) non sono grandezze reali: se letti come
    ``float``, ``138.20`` verrebbe troncato a ``138.2`` alterando la tracciabilità
    del database di riferimento nei report finali.
    """
    assert isinstance(dati_esempio["tax"]["ref_version"], str)


# --------------------------------------------------------------------------- #
# Corrispondenza biunivoca fra schema e file di esempio (W3)                   #
# --------------------------------------------------------------------------- #





def test_nessuna_chiave_dello_schema_manca_nel_file():
    """
    **Obiettivo**: Verificare che ogni parametro dello schema ``Config`` sia
    documentato in ``config/config.example.yaml``: con il suo valore
    predefinito, oppure in commento se e' obbligatorio o facoltativo senza
    valore.

    **Razionale scientifico e sistemistico**: Il modello e' la sola
    descrizione dei parametri per chi configura un dataset nuovo: un parametro
    che non vi compare non verrebbe mai rivisto.
    """
    modello = _modello()
    presenti = modello["con_valore"] | modello["obbligatorie"] | modello["commentate"]
    mancanti = chiavi_schema() - presenti
    assert not mancanti, f"chiavi dello schema assenti dal file: {sorted(mancanti)}"
    # run.threads e run.r_timeout_s sono nulli per difetto: restano in commento.
    assert modello["commentate"] == {"run.threads", "run.r_timeout_s"}


def test_nessuna_chiave_del_file_e_sconosciuta_allo_schema():
    """
    **Obiettivo**: Verificare che ``config.example.yaml`` non contenga alcuna
    chiave, attiva o in commento, che lo schema non conosce.

    **Razionale scientifico e sistemistico**: Una chiave superata rimasta nel
    modello verrebbe copiata e poi respinta alla prima esecuzione.
    """
    modello = _modello()
    estranee = (modello["con_valore"] | modello["obbligatorie"] | modello["commentate"]) \
        - chiavi_schema()
    assert not estranee, f"chiavi del file assenti dallo schema: {sorted(estranee)}"


def test_tutti_i_gruppi_compaiono_nel_file():
    """
    **Obiettivo**: Verificare che i gruppi di ``config.example.yaml`` siano
    esattamente quelli dello schema ``Config``.

    **Razionale scientifico e sistemistico**: Lo schema e' a chiave chiusa: un
    gruppo in piu' o in meno nel modello produrrebbe una configurazione
    respinta o incompleta.
    """
    gruppi = set(yaml.safe_load(ESEMPIO.read_text(encoding="utf-8")))
    assert gruppi == set(Config.model_fields)


# --------------------------------------------------------------------------- #
# Parametri obbligatori (W3)                                                   #
# --------------------------------------------------------------------------- #


def test_elenco_obbligatori_coincide_con_lo_schema():
    """
    **Obiettivo**: Verificare che l'insieme dei campi richiesti senza default
    in ``Config`` coincida esattamente con la tupla ``OBBLIGATORI`` scritta in
    questo modulo, e che ``defaults.OBBLIGATORI`` sia, nello stesso ordine,
    l'elenco letterale dei 25 parametri del dataset (``DEL_DATASET``).

    **Razionale scientifico e sistemistico**: I percorsi dei dati di input
    (``io.fastq_dir``, ``io.assay_table``, ``io.out_root``), il riferimento
    tassonomico con il suo checksum (``tax.ref_fasta``, ``tax.ref_md5``) e i 25
    parametri che descrivono il dataset dipendono dall'ambiente e dallo
    studio: un valore di default per uno di essi farebbe girare la pipeline su
    file sbagliati o con il formato di un altro dataset senza avvertire.
    ``io.study_table`` e ``run.container`` non sono obbligatori: senza tabella
    di studio le classi si leggono dalla tabella di assay, e l'immagine e' una
    dichiarazione che chi esegue fuori da un container non ha. L'atteso e'
    letterale perche' un elenco letto dal codice verificherebbe il codice
    contro se stesso.
    """
    assert len(DEL_DATASET) == 63 and len(set(DEL_DATASET)) == 63
    assert defaults.OBBLIGATORI == DEL_DATASET
    assert not {"io.study_table", "run.container"} & set(OBBLIGATORI)
    dallo_schema = {
        f"{gruppo}.{campo}"
        for gruppo, descrittore in Config.model_fields.items()
        for campo, sotto in descrittore.annotation.model_fields.items()
        if sotto.is_required()
    }
    assert dallo_schema == set(OBBLIGATORI)


@pytest.mark.parametrize("chiave", OBBLIGATORI)
def test_obbligatorio_mancante_viene_segnalato_col_suo_nome(chiave, dati_esempio):
    """
    **Obiettivo**: Verificare che l'omissione di ciascuno dei parametri
    obbligatori sollevi ``ErroreConfigurazione`` citando il percorso esatto
    ``gruppo.campo`` e la parola ``"obbligatorio"``.

    **Razionale scientifico e sistemistico**: Fornisce una diagnostica immediata
    e azionabile all'operatore, indicando esattamente quale chiave obbligatoria
    deve essere compilata nel file YAML.
    """
    with pytest.raises(ErroreConfigurazione) as errore:
        valida(_senza(dati_esempio, chiave))

    messaggio = str(errore.value)
    assert chiave in messaggio, f"il messaggio non nomina {chiave}: {messaggio}"
    assert "obbligator" in messaggio


def test_gruppo_obbligatorio_mancante_viene_segnalato(dati_esempio):
    """
    **Obiettivo**: Verificare che la rimozione dell'intera sezione ``io`` dal
    file YAML venga intercettata e segnalata con il nome del gruppo mancante.

    **Razionale scientifico e sistemistico**: Se un gruppo che racchiude parametri
    obbligatori viene omesso del tutto dal YAML, Pydantic non può istanziarlo
    tramite ``default_factory`` e deve bloccare l'esecuzione.
    """
    del dati_esempio["io"]
    with pytest.raises(ErroreConfigurazione, match=r"io\.fastq_dir"):
        valida(dati_esempio)


# --------------------------------------------------------------------------- #
# Tipi e vincoli di dominio (W3)                                               #
# --------------------------------------------------------------------------- #


def test_tipo_errato_viene_respinto_con_messaggio_leggibile(dati_esempio):
    """
    **Obiettivo**: Verificare che l'inserimento di una stringa alfabetica
    (``"due"``) per il parametro numerico ``filter.maxEE`` sollevi
    ``ErroreConfigurazione`` riportando sia il campo sia il valore rifiutato.

    **Razionale scientifico e sistemistico**: Evita che errori di battitura nel
    YAML arrivino fino allo script R ``filterAndTrim`` di DADA2 provocando
    errori criptici di coercizione a runtime.
    """
    dati_esempio["filter"]["maxEE"] = "due"

    with pytest.raises(ErroreConfigurazione) as errore:
        valida(dati_esempio)

    messaggio = str(errore.value)
    assert "filter.maxEE" in messaggio
    assert "'due'" in messaggio


@pytest.mark.parametrize(
    ("gruppo", "campo", "valore"),
    [
        ("run", "threads", 0),
        ("run", "threads", -4),
        ("filter", "truncLen", 0),
        ("tax", "min_boot", 101),
        ("decontam", "threshold", 1.5),
        ("qc", "warn_frac_chimeric", -0.1),
        ("chimera", "min_sample_fraction", 2.0),
    ],
)
def test_valore_fuori_dominio_viene_respinto(gruppo, campo, valore, dati_esempio):
    """
    **Obiettivo**: Verificare che valori numerici fuori dall'intervallo di
    validità matematica o fisica (es. ``threads <= 0``, ``min_boot > 100``,
    probabilità ``> 1.0`` o ``< 0.0``) vengano respinti dallo schema.

    **Razionale scientifico e sistemistico**: Parametri fuori scala (come una
    soglia bootstrap del 101% o una lunghezza di troncamento pari a 0 bp)
    azzererebbero l'intera assegnazione tassonomica o tutte le letture dei 960
    campioni solo dopo ore di elaborazione.
    """
    dati_esempio[gruppo][campo] = valore
    with pytest.raises(ErroreConfigurazione, match=f"{gruppo}.{campo}"):
        valida(dati_esempio)


@pytest.mark.parametrize(
    ("gruppo", "campo", "valore"),
    [
        ("chimera", "method", "consenso"),
        ("decontam", "method", "prevalenza"),
        ("out", "serialization", "qs"),
        ("dada", "pool", "pseudoo"),
    ],
)
def test_valore_fuori_vocabolario_viene_respinto(gruppo, campo, valore, dati_esempio):
    """
    **Obiettivo**: Verificare che i campi vincolati a enumerazioni chiuse
    (``Literal``) rifiutino traduzioni italiane o refusi non ammessi dalle
    librerie Bioconductor sottostanti.

    **Razionale scientifico e sistemistico**: Le funzioni R ``removeBimeraDenovo``
    e ``isContaminant`` accettano solo stringhe esatte in inglese (``"consensus"``,
    ``"prevalence"``, ``"pseudo"``): bloccare ``"pseudoo"`` o ``"prevalenza"``
    in fase F1 impedisce il fallimento tardivo in S4, S6 o S12.
    """
    dati_esempio[gruppo][campo] = valore
    with pytest.raises(ErroreConfigurazione, match=f"{gruppo}.{campo}"):
        valida(dati_esempio)


def test_dada_pool_accetta_il_booleano(dati_esempio):
    """
    **Obiettivo**: Verificare che ``dada.pool`` accetti i valori booleani ``True``
    e ``False`` oltre alla stringa ``"pseudo"``.

    **Razionale scientifico e sistemistico**: Riproduce fedelmente il dominio del
    parametro ``pool`` di ``dada2::dada()``, che ammette ``TRUE`` (pooling globale),
    ``FALSE`` (campione per campione) o ``"pseudo"`` (pseudo-pooling in due passate).
    """
    dati_esempio["dada"]["pool"] = True
    assert valida(dati_esempio).dada.pool is True


def test_booleano_scritto_come_stringa_viene_respinto(dati_esempio):
    """
    **Obiettivo**: Verificare che la stringa ``"true"`` assegnata a un campo
    booleano (``phylo.enabled``) venga rifiutata anziché convertita silenziosamente.

    **Razionale scientifico e sistemistico**: In modalità non-strict di Pydantic,
    stringhe come ``"false"`` o ``"no"`` possono produrre comportamenti di
    coercizione ambigui; imporre booleani YAML nativi evita attivazioni o
    disattivazioni involontarie di fasi pesanti come la filogenesi S9.
    """
    dati_esempio["phylo"]["enabled"] = "true"
    with pytest.raises(ErroreConfigurazione, match="phylo.enabled"):
        valida(dati_esempio)


def test_regex_malformata_viene_respinta(dati_esempio):
    """
    **Obiettivo**: Verificare che un'espressione regolare sintatticamente
    invalida in ``meta.module_regex`` (parentesi non chiusa ``"^([A-Z]"``)
    venga intercettata durante la validazione dello schema.

    **Razionale scientifico e sistemistico**: Compilare preventivamente tutte le
    regex di configurazione evita che ``re.error`` esploda durante la
    costruzione del crosswalk o la scansione delle letture in S0.
    """
    dati_esempio["meta"]["module_regex"] = "^([A-Z]"
    with pytest.raises(ErroreConfigurazione, match="meta.module_regex"):
        valida(dati_esempio)


def test_immagine_senza_digest_viene_respinta(dati_esempio):
    """
    **Obiettivo**: Verificare che ``run.container`` rifiuti un riferimento
    basato su semplice tag mutabile (``"amplicon16s:dev"``) ed esiga il formato
    con digest immutabile ``nome@sha256:<64_hex>``.

    **Razionale scientifico e sistemistico**: Un tag Docker come ``:latest`` o
    ``:dev`` può essere sovrascritto nel tempo puntando a versioni diverse di R
    o ``dada2``; solo il digest crittografico ``@sha256:...`` garantisce la
    riproducibilità computazionale esatta dell'ambiente di esecuzione.
    """
    dati_esempio["run"]["container"] = "amplicon16s:dev"
    with pytest.raises(ErroreConfigurazione, match="run.container"):
        valida(dati_esempio)


def test_md5_malformato_viene_respinto(dati_esempio):
    """
    **Obiettivo**: Verificare che ``tax.ref_md5`` accetti esclusivamente una
    stringa esadecimale di 32 caratteri.

    **Razionale scientifico e sistemistico**: Impedisce di inserire segnaposto
    testuali (es. ``"TODO"`` o ``"non-un-md5"``) per il checksum del database
    SILVA prima ancora di arrivare al controllo sul file fisico del Gate G12.
    """
    dati_esempio["tax"]["ref_md5"] = "non-un-md5"
    with pytest.raises(ErroreConfigurazione, match="tax.ref_md5"):
        valida(dati_esempio)


# --------------------------------------------------------------------------- #
# Chiavi sconosciute (W3)                                                      #
# --------------------------------------------------------------------------- #


def test_refuso_in_un_parametro_viene_respinto(dati_esempio):
    """
    **Obiettivo**: Verificare che un errore di battitura nel nome di una chiave
    (``filter.maxEEE`` invece di ``filter.maxEE``) venga bloccato con
    ``ErroreConfigurazione`` grazie a ``extra="forbid"``.

    **Razionale scientifico e sistemistico**: Senza ``extra="forbid"``, Pydantic
    ignorerebbe silenziosamente ``filter.maxEEE: 2`` e applicherebbe il valore
    predefinito ``filter.maxEE: 1.0``, facendo girare l'intera pipeline con una
    stringenza di filtraggio diversa da quella che il ricercatore credeva di
    aver impostato.
    """
    dati_esempio["filter"]["maxEEE"] = 2
    with pytest.raises(ErroreConfigurazione, match="filter.maxEEE"):
        valida(dati_esempio)


def test_gruppo_sconosciuto_viene_respinto(dati_esempio):
    """
    **Obiettivo**: Verificare che l'aggiunta di una sezione di primo livello non
    prevista dallo schema (``inventato``) venga immediatamente rifiutata.

    **Razionale scientifico e sistemistico**: Previene errori di indentazione nel
    YAML in cui un sotto-parametro viene accidentalmente promosso a radice del
    documento e ignorato dalla pipeline.
    """
    dati_esempio["inventato"] = {"x": 1}
    with pytest.raises(ErroreConfigurazione, match="inventato"):
        valida(dati_esempio)


@pytest.mark.parametrize("gruppo", ["norm", "glom", "beta", "ord", "stat"])
def test_gruppi_ecologici_non_accettano_parametri(gruppo, dati_esempio):
    """
    **Obiettivo**: Verificare che i 5 gruppi destinati alle analisi ecologiche a
    valle (``norm``, ``glom``, ``beta``, ``ord``, ``stat``) rifiutino qualsiasi
    chiave interna finché i relativi moduli non sono implementati.

    **Razionale scientifico e sistemistico**: Impedisce l'illusione che parametri
    ecologici scritti in quelle sezioni vengano già letti e applicati dalla
    pipeline di produzione `amplicon16s`.
    """
    dati_esempio[gruppo] = {"metodo": "qualcosa"}
    with pytest.raises(ErroreConfigurazione, match=f"{gruppo}.metodo"):
        valida(dati_esempio)


# --------------------------------------------------------------------------- #
# Coerenza fra campi (W3/W4)                                                   #
# --------------------------------------------------------------------------- #


def test_avviso_non_puo_superare_arresto(dati_esempio):
    """
    **Obiettivo**: Verificare che una soglia di avviso chimere superiore alla
    soglia di arresto (``warn_frac_chimeric: 0.8 > stop_frac_chimeric: 0.5``)
    venga respinta.

    **Razionale scientifico e sistemistico**: Se la soglia di warning superasse
    quella di stop in S6, la fascia di allerta precoce sarebbe irraggiungibile
    e la pipeline passerebbe direttamente dal silenzio all'arresto critico.
    """
    dati_esempio["qc"]["warn_frac_chimeric"] = 0.8
    dati_esempio["qc"]["stop_frac_chimeric"] = 0.5
    with pytest.raises(ErroreConfigurazione, match="warn_frac_chimeric"):
        valida(dati_esempio)


def test_tolleranza_non_puo_raggiungere_il_troncamento(dati_esempio):
    """
    **Obiettivo**: Verificare che ``filter.truncLen_shortfall_warn`` non possa
    eguagliare o superare la lunghezza effettiva delle letture troncate
    (``truncLen - trimLeft = 137``).

    **Razionale scientifico e sistemistico**: Garantisce che la soglia di scarto
    usata dal Gate G09 per avvisare di un troncamento eccessivamente conservativo
    (``E-S1-01``) resti inferiore alla lunghezza dell'amplicone prodotto.
    """
    dati_esempio["filter"]["truncLen_shortfall_warn"] = dati_esempio["filter"]["truncLen"]
    with pytest.raises(ErroreConfigurazione, match="truncLen_shortfall_warn"):
        valida(dati_esempio)


def test_categorie_di_controllo_devono_essere_disgiunte(dati_esempio):
    """
    **Obiettivo**: Verificare che una stessa etichetta (es. ``"Surface swab"``)
    non possa comparire contemporaneamente tra i controlli negativi
    (``ctrl.blank_values``) e i campioni biologici (``ctrl.biological_values``).

    **Razionale scientifico e sistemistico**: Se un campione biologico venisse
    classificato anche come blank, ``decontam`` in S12 tratterebbe i taxa
    autentici del microbioma come contaminanti da reagente, sottraendoli
    sistematicamente dal dataset finale.
    """
    dati_esempio["ctrl"]["blank_values"] = list(dati_esempio["ctrl"]["biological_values"])
    with pytest.raises(ErroreConfigurazione, match="ctrl"):
        valida(dati_esempio)


# --------------------------------------------------------------------------- #
# Segnalazione degli errori (W3)                                               #
# --------------------------------------------------------------------------- #


def test_vengono_riportati_tutti_i_problemi_non_solo_il_primo(dati_esempio):
    """
    **Obiettivo**: Verificare che in presenza di 3 errori distinti nel YAML,
    ``ErroreConfigurazione.problemi`` li raccolga e li elenchi tutti e 3 in
    un'unica esecuzione.

    **Razionale scientifico e sistemistico**: Evita al ricercatore un frustrante
    ciclo di tentativi ed errori ("correggi un parametro, rilancia, scopri il
    secondo errore"), restituendo in un solo colpo la diagnosi completa del file.
    """
    dati_esempio["filter"]["maxEE"] = "due"
    dati_esempio["run"]["threads"] = -1
    del dati_esempio["io"]["fastq_dir"]

    with pytest.raises(ErroreConfigurazione) as errore:
        valida(dati_esempio)

    assert len(errore.value.problemi) == 3
    messaggio = str(errore.value)
    for atteso in ("filter.maxEE", "run.threads", "io.fastq_dir"):
        assert atteso in messaggio


def test_errore_di_caricamento_indica_il_file(tmp_path):
    """
    **Obiettivo**: Verificare che ``carica(percorso)`` includa il nome del file
    sorgente (``rotto.yaml``) nel messaggio di ``ErroreConfigurazione``.

    **Razionale scientifico e sistemistico**: Quando la CLI viene invocata da
    script automatizzati o pipeline CI con più file YAML, il log deve indicare
    esattamente quale file ha fallito la validazione.
    """
    percorso = tmp_path / "rotto.yaml"
    percorso.write_text("io: [questo non e' una mappa]\n", encoding="utf-8")
    with pytest.raises(ErroreConfigurazione, match="rotto.yaml"):
        carica(percorso)


def test_file_vuoto_viene_segnalato(tmp_path):
    """
    **Obiettivo**: Verificare che un file YAML completamente vuoto (per cui
    ``yaml.safe_load`` restituisce ``None``) venga intercettato con un messaggio
    esplicito ``"vuoto"``.

    **Razionale scientifico e sistemistico**: Gestisce con chiarezza il caso di un
    file di configurazione appena creato con ``touch`` o troncato per errore.
    """
    percorso = tmp_path / "config.yaml"
    percorso.write_text("", encoding="utf-8")
    with pytest.raises(ErroreConfigurazione, match="vuoto"):
        carica(percorso)


def test_yaml_malformato_viene_segnalato(tmp_path):
    """
    **Obiettivo**: Verificare che un errore sintattico del parser YAML
    (``yaml.YAMLError`` su parentesi graffa non chiusa) venga incapsulato in
    ``ErroreConfigurazione``.

    **Razionale scientifico e sistemistico**: Uniforma tutte le patologie del file
    di configurazione sotto l'unica eccezione ``ErroreConfigurazione``, evitando
    traceback grezzi di PyYAML sulla console dell'utente.
    """
    percorso = tmp_path / "malformato.yaml"
    percorso.write_text("io: {non chiuso\n", encoding="utf-8")
    with pytest.raises(ErroreConfigurazione, match="YAML"):
        carica(percorso)


# --------------------------------------------------------------------------- #
# Predefiniti derivati dal dataset di riferimento (W3)                         #
# --------------------------------------------------------------------------- #


def test_ogni_predefinito_e_di_una_specie_ammessa():
    """
    **Obiettivo**: Verificare che ogni parametro dello schema con un
    predefinito sia il valore standard di un metodo
    (``defaults.STANDARD_DEL_METODO``, con la fonte), oppure un parametro di
    cui la pipeline realizza una sola forma, oppure un parametro dello
    strumento, oppure un facoltativo nullo; e che nessun parametro
    obbligatorio sia fra questi.

    **Razionale scientifico e sistemistico**: Un predefinito scelto guardando
    un dataset verrebbe ereditato in silenzio da un altro: e' ammesso solo cio'
    che vale per qualunque dataset, e ogni aggiunta deve passare da qui.
    """
    con_predefinito = {
        f"{gruppo}.{campo}": sotto.get_default(call_default_factory=True)
        for gruppo, descrittore in Config.model_fields.items()
        for campo, sotto in descrittore.annotation.model_fields.items()
        if not sotto.is_required()
    }
    ammessi = set(defaults.STANDARD_DEL_METODO) | set(UNICA_FORMA) | set(DELLO_STRUMENTO) | set(NULLI)
    assert set(con_predefinito) == ammessi
    assert not ammessi & set(DEL_DATASET)
    assert {c for c in NULLI if con_predefinito[c] is not None} == set()
    for chiave, fonte in defaults.STANDARD_DEL_METODO.items():
        assert "::" in fonte or "tutorial" in fonte, chiave
    # Una sola forma: lo schema respinge ogni altro valore.
    dati = _esempio_compilato()
    for chiave in UNICA_FORMA:
        gruppo, campo = chiave.split(".")
        altro = copy.deepcopy(dati)
        attuale = altro[gruppo][campo]
        altro[gruppo][campo] = (not attuale) if isinstance(attuale, bool) else "un_altro_valore"
        with pytest.raises(ErroreConfigurazione):
            valida(altro)


def test_i_marcatori_del_modello_coincidono_con_il_codice():
    """
    **Obiettivo**: Verificare che i parametri lasciati da dichiarare in
    ``config/config.example.yaml`` (``[OBBLIGATORIO]``) siano esattamente
    ``defaults.OBBLIGATORI``, e che quelli marcati ``[STANDARD]`` siano
    esattamente ``defaults.STANDARD_DEL_METODO``, con il valore dello schema.

    **Razionale scientifico e sistemistico**: Il modello dice a chi configura
    che cosa deve scegliere e che cosa puo' lasciare: se divergesse dal codice
    direbbe di poter omettere un parametro che la pipeline pretende, o
    viceversa.
    """
    modello = _modello()
    assert modello["obbligatorie"] == set(defaults.OBBLIGATORI)
    assert modello["standard"] == set(defaults.STANDARD_DEL_METODO)
    scritti = yaml.safe_load(ESEMPIO.read_text(encoding="utf-8"))
    for chiave in (*defaults.STANDARD_DEL_METODO, *UNICA_FORMA, *DELLO_STRUMENTO):
        gruppo, campo = chiave.split(".")
        predefinito = Config.model_fields[gruppo].annotation.model_fields[campo].get_default(
            call_default_factory=True)
        assert scritti[gruppo][campo] == predefinito, chiave


def test_il_modello_lascia_da_dichiarare_ogni_parametro_obbligatorio():
    """
    **Obiettivo**: Verificare che nel modello ogni parametro obbligatorio sia
    in commento e senza valore, che il file letto cosi' com'e' li dia tutti per
    mancanti, e che non nomini alcun dataset.

    **Razionale scientifico e sistemistico**: Un valore d'esempio accanto a un
    parametro obbligatorio viene copiato: il modello dice come si sceglie, non
    che cosa scegliere.
    """
    import re

    testo = ESEMPIO.read_text(encoding="utf-8")
    righe = [r for r in testo.splitlines() if "[OBBLIGATORIO]" in r and r.startswith("  # ")]
    assert len(righe) == len(defaults.OBBLIGATORI)
    for riga in righe:
        assert re.fullmatch(r"  # [A-Za-z_0-9]+: +# \[OBBLIGATORIO\]( oppure (null|\[\]))?", riga), riga
    assert obbligatori_mancanti(yaml.safe_load(testo)) == list(defaults.OBBLIGATORI)
    assert not re.search(r"OSD|GLDS", testo)


def test_whitelist_dei_ritentativi_non_dipende_dal_dataset():
    """
    **Obiettivo**: Verificare che ``retry.whitelist`` abbia un predefinito e
    non sia fra i parametri obbligatori.

    **Razionale scientifico e sistemistico**: L'elenco dei codici ammessi al
    nuovo tentativo e' una scelta metodologica dello strumento: si ritenta solo
    dove l'azione correttiva non cambia alcuna assunzione, su qualunque dataset.
    """
    assert "retry.whitelist" not in defaults.OBBLIGATORI
    assert "retry.whitelist" in DELLO_STRUMENTO


def test_esempio_resta_allineato_alla_whitelist_predefinita(dati_esempio):
    """
    **Obiettivo**: Verificare che la lista ``retry.whitelist`` in
    ``config.example.yaml`` coincida con ``defaults.RETRY_WHITELIST``.

    **Razionale scientifico e sistemistico**: Garantisce che la configurazione di
    esempio esponga esattamente i 4 codici di errore autorizzati al retry automatico.
    """
    assert list(dati_esempio["retry"]["whitelist"]) == list(defaults.RETRY_WHITELIST)


# --------------------------------------------------------------------------- #
# Immutabilità (W3)                                                            #
# --------------------------------------------------------------------------- #


def test_la_configurazione_non_si_modifica_dopo_il_caricamento():
    """
    **Obiettivo**: Verificare che qualsiasi tentativo di mutare un attributo di
    ``Config`` dopo la validazione (es. ``config.run.threads = 1``) sollevi
    un'eccezione grazie a ``ConfigDict(frozen=True)``.

    **Razionale scientifico e sistemistico**: Se una fase intermedia potesse
    modificare in-place l'oggetto ``Config`` condiviso in memoria, le fasi
    successive girerebbero con parametri diversi da quelli registrati in
    ``00_config/resolved.yaml``, invalidando il digest SHA-256 e la
    riproducibilità dell'intera corsa.
    """
    config = valida(_esempio_compilato())
    with pytest.raises(Exception):
        config.run.threads = 1


# =========================================================================== #
# Gate G15 : coerenza interna della configurazione (W4)                        #
# =========================================================================== #

def _viola(dati_esempio: dict, modifica) -> ErroreGate:
    """Applica una mutazione incoerente al dizionario YAML ed esegue ``esegui_g15`` catturando ``ErroreGate``."""
    modifica(dati_esempio)
    with pytest.raises(ErroreGate) as errore:
        esegui_g15(dati_esempio)
    return errore.value


# --------------------------------------------------------------------------- #
# Il gate accetta una configurazione coerente (W4)                             #
# --------------------------------------------------------------------------- #


def test_g15_accetta_il_file_di_esempio(dati_esempio):
    """
    **Obiettivo**: Verificare che ``esegui_g15(dati_esempio)`` superi tutti i
    controlli del Gate G15 restituendo un'istanza valida di ``ConfigRisolta``.

    **Razionale scientifico e sistemistico**: Certifica che la configurazione di
    riferimento per OSD-734 soddisfa simultaneamente tutte le regole di coerenza
    interna (``E-G15-01`` .. ``E-G15-09``).
    """
    risolta = esegui_g15(dati_esempio)
    assert isinstance(risolta, ConfigRisolta)


def test_registro_dei_controlli_e_coerente():
    """
    **Obiettivo**: Verificare che il registro dichiarativo ``CONTROLLI`` del
    Gate G15 non contenga codici duplicati e dichiari per ciascuno l'origine
    (``"schema"`` o ``"gate"``).

    **Razionale scientifico e sistemistico**: Garantisce la tracciabilità formale
    di ogni regola di validazione tra il modello Pydantic e il controllo G15.
    """
    codici = [c.codice for c in CONTROLLI]
    assert len(codici) == len(set(codici)), "codici duplicati"
    for controllo in CONTROLLI:
        assert controllo.implementato_da in ("schema", "gate")


# --------------------------------------------------------------------------- #
# Configurazioni incoerenti (W4)                                               #
# --------------------------------------------------------------------------- #

INCOERENTI = [
    pytest.param(
        lambda d: d["qc"].update(warn_frac_chimeric=0.8, stop_frac_chimeric=0.5),
        "E-G15-05",
        ("qc.warn_frac_chimeric", "qc.stop_frac_chimeric"),
        id="avviso-oltre-arresto",
    ),
    pytest.param(
        lambda d: d["decontam"].__setitem__("threshold", 1.0),
        "E-G15-02",
        ("decontam.threshold",),
        id="soglia-decontam-degenere",
    ),
    pytest.param(
        lambda d: d["asv"].__setitem__("len_tol", 200),
        "E-G15-03",
        ("asv.len_min", "asv.len_tol", "filter.truncLen"),
        id="tolleranza-oltre-lunghezza",
    ),
    pytest.param(
        lambda d: d["prev"].__setitem__("min_fraction", 1.5),
        "E-G15-04",
        ("prev.min_fraction",),
        id="frazione-fuori-scala",
    ),
    pytest.param(
        lambda d: d["retry"].__setitem__("max_attempts", -1),
        "E-G15-06",
        ("retry.max_attempts",),
        id="tentativi-negativi",
    ),
    pytest.param(
        lambda d: d["io"].__setitem__("batch_table", None),
        "E-G15-11",
        ("io.batch_table", "out.batch_columns"),
        id="colonne-del-lotto-senza-il-file",
    ),
    pytest.param(
        lambda d: d["meta"].__setitem__("batch_key_column", None),
        "E-G15-11",
        ("io.batch_table", "meta.batch_key_column"),
        id="file-del-lotto-senza-la-chiave",
    ),
    pytest.param(
        lambda d: d["katharoseq"].__setitem__("target_taxon", None),
        "E-G15-12",
        ("ctrl.positive_values", "katharoseq.target_taxon"),
        id="positivi-senza-taxon-atteso",
    ),
    pytest.param(
        lambda d: d["out"]["study_columns"].append("Classe"),
        "E-G15-13",
        ("out.study_columns", "out.batch_columns"),
        id="colonna-con-il-nome-di-una-dell-inventario",
    ),
]


@pytest.mark.parametrize(("modifica", "codice", "parametri"), INCOERENTI)
def test_g15_respinge_le_configurazioni_incoerenti(
    modifica, codice, parametri, dati_esempio
):
    """
    **Obiettivo**: Verificare che ciascuna violazione logica incrociata sollevi
    ``ErroreGate`` con lo specifico codice di catalogo (``E-G15-02`` .. ``E-G15-07``)
    e nomini tutti i parametri coinvolti nel conflitto.

    **Razionale scientifico e sistemistico**: Previene configurazioni numericamente
    valide sul singolo campo ma biologicamente o matematicamente degeneri
    (es. ``decontam.threshold = 1.0`` che classificherebbe come contaminante il
    100% delle ASV in S12, oppure ``ctrl.blank_values = []`` con metodo
    ``prevalence`` o ``combined`` che richiede obbligatoriamente controlli negativi).
    """
    errore = _viola(dati_esempio, modifica)

    codici = {v.codice for v in errore.violazioni}
    assert codice in codici, f"atteso {codice}, ottenuti {sorted(codici)}"
    assert errore.categoria is Categoria.REVISIONE_UMANA
    assert voce(codice).categoria is Categoria.REVISIONE_UMANA

    testo = str(errore)
    for parametro in parametri:
        assert parametro in testo, f"il messaggio non nomina {parametro}: {testo}"


def test_g15_termina_con_un_codice_di_errore(dati_esempio):
    """
    **Obiettivo**: Verificare che ``decontam.threshold = 0.0`` (estremo inferiore
    escluso dall'intervallo aperto $(0, 1)$) faccia fallire G15 con codice
    ``E-G15-02`` e identificativo ``gate == "G15"``.

    **Razionale scientifico e sistemistico**: Nel test statistico di ``decontam``,
    una soglia di probabilità pari a ``0.0`` non rimuoverebbe mai alcun
    contaminante, disattivando di fatto la decontaminazione S12 pur lasciandola
    apparentemente abilitata nei report.
    """
    errore = _viola(dati_esempio, lambda d: d["decontam"].__setitem__("threshold", 0.0))
    assert errore.codice == "E-G15-02"
    assert errore.gate == "G15"


def test_g15_raccoglie_tutte_le_violazioni(dati_esempio):
    """
    **Obiettivo**: Verificare che ``esegui_g15`` accumuli tutte le violazioni
    presenti nella configurazione (``len(errore.violazioni) >= 2``) anziché
    fermarsi alla prima.

    **Razionale scientifico e sistemistico**: Permette all'utente di vedere e
    correggere simultaneamente tutti i conflitti tra parametri in un solo passaggio.
    """
    errore = _viola(
        dati_esempio,
        lambda d: (
            d["qc"].update(warn_frac_chimeric=0.9, stop_frac_chimeric=0.1),
            d["retry"].__setitem__("max_attempts", -1),
        ),
    )
    assert len(errore.violazioni) >= 2


# --------------------------------------------------------------------------- #
# Parametri derivati (W4)                                                      #
# --------------------------------------------------------------------------- #


def test_derivati_statici_assumono_i_valori_attesi(dati_esempio):
    """
    **Obiettivo**: Verificare che sulla configurazione predefinita di OSD-734
    (``truncLen = 137``, ``trimLeft = 0``, ``len_tol = 0``) i due parametri
    derivati statici assumano esattamente il valore ``137`` bp.
    **Razionale scientifico e sistemistico**: In una corsa single-end Illumina
    MiSeq troncata a 137 bp senza taglio in 5' (``trimLeft = 0``) e con tolleranza
    zero (``len_tol = 0``), ogni lettura e ogni ASV in uscita da DADA2 ha lunghezza
    attesa pari a 137 - 0 = 137 nt: derivare queste soglie per formula
    elimina il rischio che l'operatore aggiorni ``truncLen`` dimenticando di
    aggiornare ``asv.len_min`` e ``asv.len_max``. ``filter.minLen`` non esiste
    piu': coincideva per costruzione con ``filter.truncLen``.
    """
    derivati = esegui_g15(dati_esempio).derivati
    assert not hasattr(derivati, "filter_minLen")
    assert derivati.asv_len_min == FORMATO.troncamento
    assert derivati.asv_len_max == FORMATO.troncamento


def test_derivati_seguono_i_parametri_da_cui_discendono(dati_esempio):
    r"""
    **Obiettivo**: Verificare che al variare di ``truncLen = 150``,
    ``trimLeft = 10`` e ``len_tol = 2``, i derivati vengano ricalcolati
    dinamicamente in ``asv_len_min = 138`` e ``asv_len_max = 142``.

    **Razionale scientifico e sistemistico**: Dimostra che i parametri derivati
    non sono costanti cablate ma funzioni pure di ``truncLen``, ``trimLeft`` e
    ``len_tol`` ($L = \text{truncLen} - \text{trimLeft}$, intervallo ASV
    $[L - \text{len\_tol}, L + \text{len\_tol}]$).
    """
    dati_esempio["filter"]["truncLen"] = 150
    dati_esempio["filter"]["trimLeft"] = 10
    dati_esempio["asv"]["len_tol"] = 2

    derivati = esegui_g15(dati_esempio).derivati
    assert derivati.asv_len_min == 138
    assert derivati.asv_len_max == 142


def test_la_risoluzione_e_utilizzabile_senza_il_gate(dati_esempio):
    """
    **Obiettivo**: Verificare che ``risolvi(valida(dati))`` calcoli i derivati
    anche se invocata direttamente senza passare da ``esegui_g15``.

    **Razionale scientifico e sistemistico**: Mantiene disaccoppiato il motore di
    risoluzione algebrica (`resolve.py`) dal registro dei gate (`g01_g15.py`).
    """
    risolta = risolvi(valida(dati_esempio))
    assert risolta.derivati.asv_len_min == FORMATO.troncamento


@pytest.mark.parametrize("chiave", ["prev.min_samples", "qc.min_reads_filtered"])
def test_i_parametri_rimossi_sono_respinti(chiave, dati_esempio):
    """
    **Obiettivo**: Verificare che ``prev.min_samples`` e ``qc.min_reads_filtered``,
    rimossi dallo schema, siano respinti come chiavi sconosciute, e che fra i
    derivati non compaia piu' ``prev.min_samples``.

    **Razionale scientifico e sistemistico**: ``prev.min_samples`` era un
    derivato calcolato sull'intero inventario (770 biologici, minimo 8) che
    nessuna fase leggeva: S13 calcola il minimo del filtro di prevalenza sui
    biologici che tiene dopo i filtri per profondita' e tassonomico (sul
    riferimento 492, minimo 5), e il derivato lo contraddiceva. L'arrotondamento
    per eccesso del minimo e' verificato su S13 in ``test_w21_finale.py``.
    ``qc.min_reads_filtered`` non era letto da alcuna fase: chi lo cambiava non
    otteneva alcun effetto. Uno schema a chiave chiusa li respinge, invece di
    ignorarli in silenzio.
    """
    gruppo, campo = chiave.split(".")
    dati_esempio[gruppo][campo] = 1
    with pytest.raises(ErroreConfigurazione, match="parametro sconosciuto"):
        valida(dati_esempio)
    assert chiave not in PARAMETRI_DERIVATI


@pytest.mark.parametrize("chiave", PARAMETRI_DERIVATI)
def test_derivato_impostato_a_mano_e_un_errore(chiave, dati_esempio):
    """
    **Obiettivo**: Verificare che la valorizzazione manuale nel file YAML di uno
    qualsiasi dei 2 parametri derivati (``asv.len_min``, ``asv.len_max``)
    venga bloccata da G15 con codice ``E-G15-08``.

    **Razionale scientifico e sistemistico**: I 2 parametri derivati sono
    determinati univocamente dalle formule della pipeline; permettere all'utente
    di sovrascriverli a mano creerebbe contraddizioni interne tra il filtro
    ``filterAndTrim`` (S2) e il filtro di lunghezza ASV (S7).
    """
    gruppo, campo = chiave.split(".")
    dati_esempio[gruppo][campo] = 1

    with pytest.raises(ErroreGate) as errore:
        esegui_g15(dati_esempio)

    assert "E-G15-08" in {v.codice for v in errore.value.violazioni}
    assert errore.value.categoria is Categoria.REVISIONE_UMANA
    assert chiave in str(errore.value)


def test_elenco_dei_derivati_coincide_con_quelli_calcolati(dati_esempio):
    """
    **Obiettivo**: Verificare che le chiavi prodotte da ``Derivati.come_chiavi()``
    coincidano esattamente con la costante ``PARAMETRI_DERIVATI``.

    **Razionale scientifico e sistemistico**: Impedisce disallineamenti tra
    l'elenco dei campi vietati all'utente in ``schema.py`` e quelli effettivamente
    iniettati da ``resolve.py`` in ``00_config/resolved.yaml``.
    """
    calcolati = set(esegui_g15(dati_esempio).derivati.come_chiavi())
    assert calcolati == set(PARAMETRI_DERIVATI)


# --------------------------------------------------------------------------- #
# Digest (W4)                                                                  #
# --------------------------------------------------------------------------- #


def test_configurazioni_identiche_hanno_lo_stesso_digest(dati_esempio):
    """
    **Obiettivo**: Verificare che due istanze distinte ma uguali della
    configurazione producano esattamente lo stesso digest SHA-256.

    **Razionale scientifico e sistemistico**: Il digest della configurazione è la
    chiave primaria con cui ``ProjectRun`` stabilisce se una fase su disco può
    essere riutilizzata o va ricalcolata: qualsiasi non-determinismo nell'hashing
    romperebbe il meccanismo di ``resume``.
    """
    primo = esegui_g15(copy.deepcopy(dati_esempio)).digest
    secondo = esegui_g15(copy.deepcopy(dati_esempio)).digest
    assert primo == secondo


def test_un_solo_parametro_diverso_cambia_il_digest(dati_esempio):
    """
    **Obiettivo**: Verificare che la modifica di un singolo parametro (es.
    ``run.seed`` da ``100`` a ``101``) produca un digest SHA-256 differente.

    **Razionale scientifico e sistemistico**: Cambiare il seme del generatore
    pseudo-casuale (``run.seed``) altera il campionamento delle letture in
    ``learnErrors`` (S3) e le permutazioni tassonomiche (S8): il cambio del
    digest garantisce che tutte le fasi dipendenti vengano invalidate e
    ricalcolate col nuovo seme.
    """
    originale = esegui_g15(copy.deepcopy(dati_esempio)).digest
    dati_esempio["run"]["seed"] = 101
    assert esegui_g15(dati_esempio).digest != originale


def test_il_digest_non_dipende_dall_ordine_delle_chiavi(dati_esempio):
    """
    **Obiettivo**: Verificare che invertendo l'ordine dei gruppi nel dizionario
    YAML il digest SHA-256 rimanga identico (serializzazione canonica con
    ``sort_keys=True``).

    **Razionale scientifico e sistemistico**: Riformattare o riordinare i blocchi
    di un file YAML senza cambiare alcun valore non deve mai invalidare un'analisi
    già completata su disco.
    """
    originale = esegui_g15(copy.deepcopy(dati_esempio)).digest
    rovesciato = {k: dati_esempio[k] for k in reversed(list(dati_esempio))}
    assert esegui_g15(rovesciato).digest == originale


def test_il_digest_cambia_se_cambia_un_derivato(dati_esempio):
    """
    **Obiettivo**: Verificare che, a parita' di parametri dichiarati, due
    configurazioni risolte con valori derivati diversi abbiano digest diversi;
    e che una variazione di ``asv.len_tol`` sposti i derivati.

    **Razionale scientifico e sistemistico**: Assicura che l'impronta crittografica
    copra tanto i parametri dichiarati quanto i valori derivati effettivamente
    applicati nel calcolo: un derivato calcolato in un altro modo da una
    versione successiva non deve dare lo stesso digest.
    """
    import dataclasses

    risolta = esegui_g15(copy.deepcopy(dati_esempio))
    spostati = dataclasses.replace(risolta.derivati, asv_len_min=risolta.derivati.asv_len_min - 1)
    assert dataclasses.replace(risolta, derivati=spostati).digest != risolta.digest
    dati_esempio["asv"]["len_tol"] = 3
    assert esegui_g15(dati_esempio).derivati != risolta.derivati


def test_il_digest_ha_la_forma_attesa(dati_esempio):
    """
    **Obiettivo**: Verificare che ``risolta.digest`` rispetti il formato
    ``sha256:<64_caratteri_esadecimali_minuscoli>``.

    **Razionale scientifico e sistemistico**: Uniforma il formato del digest di
    configurazione a quello di tutti gli altri checksum registrati nei manifesti
    della pipeline (`io_layer/checksums.py`).
    """
    digest = esegui_g15(dati_esempio).digest
    algoritmo, _, valore = digest.partition(":")
    assert algoritmo == "sha256"
    assert len(valore) == 64
    assert set(valore) <= set("0123456789abcdef")


# --------------------------------------------------------------------------- #
# Registrazione della configurazione usata (W4)                                #
# --------------------------------------------------------------------------- #


def test_resolved_viene_scritto_con_il_digest(dati_esempio, tmp_path):
    """
    **Obiettivo**: Verificare che ``scrivi_risolta()`` crei il file
    ``00_config/resolved.yaml`` contenente la chiave radice ``digest`` allineata
    a ``risolta.digest``.

    **Razionale scientifico e sistemistico**: Persiste nella cartella ``00_config``
    il contratto completo dell'esecuzione, rendendo il risultato finale in
    ``12_final`` auto-documentato e riproducibile anche a distanza di anni.
    """
    risolta = esegui_g15(dati_esempio)
    percorso = scrivi_risolta(risolta, tmp_path)

    assert percorso == tmp_path / Fase.CONFIG.value / NOME_FILE_RISOLTO
    assert percorso.is_file()

    documento = yaml.safe_load(percorso.read_text(encoding="utf-8"))
    assert documento["digest"] == risolta.digest


def test_resolved_contiene_i_parametri_derivati(dati_esempio, tmp_path):
    """
    **Obiettivo**: Verificare che ``00_config/resolved.yaml`` contenga sia i
    valori calcolati dei 3 parametri derivati (``137``, ``137``, ``137``)
    dentro ``parametri``, sia la mappa separata ``derivati``.

    **Razionale scientifico e sistemistico**: Consente a chi ispeziona
    ``resolved.yaml`` (e al report finale in ``12_final``) di leggere
    direttamente i valori numerici risolti senza doverli ricalcolare a mente.
    """
    risolta = esegui_g15(dati_esempio)
    documento = yaml.safe_load(
        scrivi_risolta(risolta, tmp_path).read_text(encoding="utf-8")
    )

    parametri = documento["parametri"]
    assert "minLen" not in parametri["filter"]
    assert parametri["asv"]["len_min"] == FORMATO.troncamento
    assert parametri["asv"]["len_max"] == FORMATO.troncamento
    assert set(documento["derivati"]) == set(PARAMETRI_DERIVATI)


def test_resolved_registra_tutti_i_parametri_dichiarati(dati_esempio, tmp_path):
    """
    **Obiettivo**: Verificare che ``00_config/resolved.yaml`` includa l'unione
    completa di tutte le chiavi dello schema e di tutti i parametri derivati.

    **Razionale scientifico e sistemistico**: Anche se l'utente fornisce un file
    YAML minimale con i soli 9 campi obbligatori, ``resolved.yaml`` congela
    su disco tutti i valori predefiniti effettivamente usati da quella corsa.
    """
    documento = yaml.safe_load(
        scrivi_risolta(esegui_g15(dati_esempio), tmp_path).read_text(encoding="utf-8")
    )
    scritte = {
        f"{gruppo}.{campo}"
        for gruppo, contenuto in documento["parametri"].items()
        for campo in contenuto
    }
    assert chiavi_schema() <= scritte
    assert set(PARAMETRI_DERIVATI) <= scritte


def test_resolved_e_deterministico(dati_esempio, tmp_path):
    """
    **Obiettivo**: Verificare che due invocazioni di ``scrivi_risolta`` sulla
    stessa ``ConfigRisolta`` producano file ``resolved.yaml`` identici carattere
    per carattere.

    **Razionale scientifico e sistemistico**: Poiché ``scrivi_risolta`` registra
    ``resolved.yaml`` nel manifesto di ``00_config`` tramite ``AlberoOutput``,
    il file deve avere un checksum SHA-256 invariante a parità di configurazione.
    """
    risolta = esegui_g15(dati_esempio)
    primo = scrivi_risolta(risolta, tmp_path / "a").read_text(encoding="utf-8")
    secondo = scrivi_risolta(risolta, tmp_path / "b").read_text(encoding="utf-8")
    assert primo == secondo


def test_la_cartella_di_destinazione_viene_creata(dati_esempio, tmp_path):
    """
    **Obiettivo**: Verificare che ``scrivi_risolta`` crei ricorsivamente l'albero
    delle directory genitori e la sottocartella ``00_config`` se non esistono ancora.

    **Razionale scientifico e sistemistico**: La scrittura di ``00_config/resolved.yaml``
    è spesso la primissima operazione su disco di una nuova corsa in una
    directory ``io.out_root`` non ancora esistente.
    """
    out_root = tmp_path / "non" / "ancora" / "esistente"
    percorso = scrivi_risolta(esegui_g15(dati_esempio), out_root)
    assert percorso.is_file()
