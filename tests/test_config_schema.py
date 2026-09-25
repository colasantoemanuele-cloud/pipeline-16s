"""Suite di verifica dello schema di configurazione, del Gate G15 e dei parametri derivati.

Inquadramento nel Piano Operativo:
    - **Settimane di riferimento**:
        * **Settimana 3 (W3 — Fase F1: Schema di validazione della configurazione)**:
          modellazione dei 22 gruppi Pydantic (17 di produzione + 5 ecologici),
          parametri obbligatori privi di default, vincoli di dominio, blocco dei
          refusi (``extra="forbid"``), immutabilità (``frozen=True``) e allineamento
          con il dataset di riferimento **OSD-734**.
        * **Settimana 4 (W4 — Fase F1: Coerenza interna della configurazione e
          parametri derivati)**: controlli relazionali incrociati del **Gate G15**
          (codici ``E-G15-01`` .. ``E-G15-09``), derivazione matematica di
          ``filter.minLen``, ``asv.len_min``, ``asv.len_max`` e ``prev.min_samples``
          (con arrotondamento per eccesso ``math.ceil``), calcolo del digest
          canonico SHA-256 e scrittura di ``00_config/resolved.yaml``.
    - **Scopo del modulo**: Garantire che qualsiasi errore formale, semantico o
      metodologico nel file YAML dell'utente venga intercettato istantaneamente
      prima di allocare risorse computazionali o avviare processi R sui 960
      campioni del dataset.
    - **Moduli sorgente coperti**:
        * ``src/amplicon16s/config/schema.py``
        * ``src/amplicon16s/config/defaults.py``
        * ``src/amplicon16s/config/resolve.py``
        * ``src/amplicon16s/gates/g01_g15.py`` (limitatamente a ``esegui_g15`` e ``CONTROLLI``)
        * ``config/config.example.yaml``
"""

from __future__ import annotations

import copy
from pathlib import Path

import pytest
import yaml

from amplicon16s.config import defaults
from amplicon16s.config.resolve import (
    NOME_FILE_RISOLTO,
    ConfigRisolta,
    risolvi,
    scrivi_risolta,
)
from amplicon16s.config.schema import (
    PARAMETRI_DERIVATI,
    Config,
    ErroreConfigurazione,
    carica,
    chiavi_schema,
    gruppi_schema,
    valida,
)
from amplicon16s.gates.g01_g15 import CONTROLLI, ErroreGate, esegui_g15
from amplicon16s.io_layer.artifacts import Fase

RADICE = Path(__file__).resolve().parents[1]
ESEMPIO = RADICE / "config" / "config.example.yaml"

#: Elenco tassativo dei parametri obbligatori privi di valore predefinito nello schema.
OBBLIGATORI = (
    "io.fastq_dir",
    "io.assay_table",
    "io.study_table",
    "io.out_root",
    "tax.ref_fasta",
    "tax.ref_name",
    "tax.ref_version",
    "tax.ref_md5",
    "run.container",
)


@pytest.fixture
def dati_esempio() -> dict:
    """Carica in un dizionario mutabile indipendente il file ``config.example.yaml``."""
    return yaml.safe_load(ESEMPIO.read_text(encoding="utf-8"))


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

    **Razionale Scientifico/Sistemistico**: Il file di esempio costituisce il
    template documentato di riferimento per l'operatore e il contratto vivente
    dei parametri calibrati per il dataset NASA GeneLab **OSD-734**.
    """
    assert ESEMPIO.is_file(), f"file di esempio non trovato in {ESEMPIO}"


def test_esempio_si_carica_e_valida():
    """
    **Obiettivo**: Verificare che ``carica(ESEMPIO)`` superi l'intera validazione
    Pydantic producendo un'istanza valida di ``Config``.

    **Razionale Scientifico/Sistemistico**: Impedisce che modifiche allo schema
    Python lascino disallineato o invalido il file di configurazione di esempio
    distribuito con la pipeline.
    """
    config = carica(ESEMPIO)
    assert isinstance(config, Config)


def test_esempio_conserva_i_tipi_numerici(dati_esempio):
    """
    **Obiettivo**: Verificare che i parametri espressi in notazione scientifica
    nel YAML (``err.nbases = 1.0e+8``, ``dada.omega_a = 1.0e-40``) siano letti
    da PyYAML come ``float`` nativi e non come stringhe.

    **Razionale Scientifico/Sistemistico**: Secondo la specifica YAML 1.1 usata
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

    **Razionale Scientifico/Sistemistico**: I numeri di versione tassonomica
    (come SILVA ``138.1`` o ``138.20``) non sono grandezze reali: se letti come
    ``float``, ``138.20`` verrebbe troncato a ``138.2`` alterando la tracciabilità
    del database di riferimento nei report finali.
    """
    assert isinstance(dati_esempio["tax"]["ref_version"], str)


# --------------------------------------------------------------------------- #
# Corrispondenza biunivoca fra schema e file di esempio (W3)                   #
# --------------------------------------------------------------------------- #


def _chiavi_del_file(dati: dict) -> set[str]:
    """Estrae l'insieme di tutte le chiavi ``gruppo.campo`` presenti nel dizionario YAML."""
    return {
        f"{gruppo}.{campo}"
        for gruppo, contenuto in dati.items()
        if isinstance(contenuto, dict)
        for campo in contenuto
    }


def test_nessuna_chiave_dello_schema_manca_nel_file(dati_esempio):
    """
    **Obiettivo**: Verificare che ogni parametro definito nello schema Pydantic
    ``Config`` sia esplicitamente documentato in ``config/config.example.yaml``.

    **Razionale Scientifico/Sistemistico**: Evita la presenza di "parametri
    nascosti" nel codice Python di cui il biologo computazionale ignora
    l'esistenza consultando il file di configurazione di riferimento.
    """
    mancanti = chiavi_schema() - _chiavi_del_file(dati_esempio)
    assert not mancanti, f"chiavi dello schema assenti dal file: {sorted(mancanti)}"


def test_nessuna_chiave_del_file_e_sconosciuta_allo_schema(dati_esempio):
    """
    **Obiettivo**: Verificare che ``config.example.yaml`` non contenga alcuna
    chiave obsoleta o estranea allo schema ``Config``.

    **Razionale Scientifico/Sistemistico**: Garantisce la sincronia bidirezionale
    1:1 tra documentazione YAML e modello formale di validazione.
    """
    estranee = _chiavi_del_file(dati_esempio) - chiavi_schema()
    assert not estranee, f"chiavi del file assenti dallo schema: {sorted(estranee)}"


def test_tutti_i_gruppi_compaiono_nel_file(dati_esempio):
    """
    **Obiettivo**: Verificare che tutti i 22 gruppi di primo livello (inclusi i
    5 gruppi ecologici ``norm``, ``glom``, ``beta``, ``ord``, ``stat`` privi di
    sotto-chiavi) siano presenti in ``config.example.yaml``.

    **Razionale Scientifico/Sistemistico**: Poiché i gruppi riservati alle
    analisi ecologiche a valle sono mappe vuote ``{}``, un controllo basato solo
    sulle foglie ``gruppo.campo`` non li vedrebbe; questo test assicura che
    l'intera architettura a 22 sezioni sia rappresentata nel file di esempio.
    """
    mancanti = set(gruppi_schema()) - set(dati_esempio)
    assert not mancanti, f"gruppi assenti dal file: {sorted(mancanti)}"
    estranei = set(dati_esempio) - set(gruppi_schema())
    assert not estranei, f"gruppi del file assenti dallo schema: {sorted(estranei)}"


# --------------------------------------------------------------------------- #
# Parametri obbligatori (W3)                                                   #
# --------------------------------------------------------------------------- #


def test_elenco_obbligatori_coincide_con_lo_schema():
    """
    **Obiettivo**: Verificare che l'insieme dei campi richiesti senza default
    in ``Config`` coincida esattamente con la tupla ``OBBLIGATORI``.

    **Razionale Scientifico/Sistemistico**: I percorsi dei dati di input
    (``io.fastq_dir``, ``io.assay_table``, ``io.study_table``, ``io.out_root``),
    il riferimento tassonomico con il suo checksum (``tax.ref_*``) e il digest
    del container (``run.container``) dipendono dall'ambiente e dallo studio:
    introdurre un valore di default per uno di essi rischierebbe di far girare
    la pipeline su file o database sbagliati senza avvertire l'utente.
    """
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
    **Obiettivo**: Verificare che l'omissione di ciascuno dei 9 parametri
    obbligatori sollevi ``ErroreConfigurazione`` citando il percorso esatto
    ``gruppo.campo`` e la parola ``"obbligatorio"``.

    **Razionale Scientifico/Sistemistico**: Fornisce una diagnostica immediata
    e azionabile all'operatore, indicando esattamente quale chiave obbligatoria
    deve essere compilata nel file YAML.
    """
    with pytest.raises(ErroreConfigurazione) as errore:
        valida(_senza(dati_esempio, chiave))

    messaggio = str(errore.value)
    assert chiave in messaggio, f"il messaggio non nomina {chiave}: {messaggio}"
    assert "obbligatorio" in messaggio


def test_gruppo_obbligatorio_mancante_viene_segnalato(dati_esempio):
    """
    **Obiettivo**: Verificare che la rimozione dell'intera sezione ``io`` dal
    file YAML venga intercettata e segnalata con il nome del gruppo mancante.

    **Razionale Scientifico/Sistemistico**: Se un gruppo che racchiude parametri
    obbligatori viene omesso del tutto dal YAML, Pydantic non può istanziarlo
    tramite ``default_factory`` e deve bloccare l'esecuzione.
    """
    del dati_esempio["io"]
    with pytest.raises(ErroreConfigurazione, match="io"):
        valida(dati_esempio)


# --------------------------------------------------------------------------- #
# Tipi e vincoli di dominio (W3)                                               #
# --------------------------------------------------------------------------- #


def test_tipo_errato_viene_respinto_con_messaggio_leggibile(dati_esempio):
    """
    **Obiettivo**: Verificare che l'inserimento di una stringa alfabetica
    (``"due"``) per il parametro numerico ``filter.maxEE`` sollevi
    ``ErroreConfigurazione`` riportando sia il campo sia il valore rifiutato.

    **Razionale Scientifico/Sistemistico**: Evita che errori di battitura nel
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

    **Razionale Scientifico/Sistemistico**: Parametri fuori scala (come una
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

    **Razionale Scientifico/Sistemistico**: Le funzioni R ``removeBimeraDenovo``
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

    **Razionale Scientifico/Sistemistico**: Riproduce fedelmente il dominio del
    parametro ``pool`` di ``dada2::dada()``, che ammette ``TRUE`` (pooling globale),
    ``FALSE`` (campione per campione) o ``"pseudo"`` (pseudo-pooling in due passate).
    """
    dati_esempio["dada"]["pool"] = True
    assert valida(dati_esempio).dada.pool is True


def test_booleano_scritto_come_stringa_viene_respinto(dati_esempio):
    """
    **Obiettivo**: Verificare che la stringa ``"true"`` assegnata a un campo
    booleano (``phylo.enabled``) venga rifiutata anziché convertita silenziosamente.

    **Razionale Scientifico/Sistemistico**: In modalità non-strict di Pydantic,
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

    **Razionale Scientifico/Sistemistico**: Compilare preventivamente tutte le
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

    **Razionale Scientifico/Sistemistico**: Un tag Docker come ``:latest`` o
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

    **Razionale Scientifico/Sistemistico**: Impedisce di inserire segnaposto
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

    **Razionale Scientifico/Sistemistico**: Senza ``extra="forbid"``, Pydantic
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

    **Razionale Scientifico/Sistemistico**: Previene errori di indentazione nel
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

    **Razionale Scientifico/Sistemistico**: Impedisce l'illusione che parametri
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

    **Razionale Scientifico/Sistemistico**: Se la soglia di warning superasse
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

    **Razionale Scientifico/Sistemistico**: Garantisce che la soglia di scarto
    usata dal Gate G09 per avvisare di un troncamento eccessivamente conservativo
    (``E-S1-01``) resti inferiore alla lunghezza dell'amplicone prodotto.
    """
    dati_esempio["filter"]["truncLen_shortfall_warn"] = 137
    with pytest.raises(ErroreConfigurazione, match="truncLen_shortfall_warn"):
        valida(dati_esempio)


def test_categorie_di_controllo_devono_essere_disgiunte(dati_esempio):
    """
    **Obiettivo**: Verificare che una stessa etichetta (es. ``"Surface swab"``)
    non possa comparire contemporaneamente tra i controlli negativi
    (``ctrl.blank_values``) e i campioni biologici (``ctrl.biological_values``).

    **Razionale Scientifico/Sistemistico**: Se un campione biologico venisse
    classificato anche come blank, ``decontam`` in S12 tratterebbe i taxa
    autentici del microbioma come contaminanti da reagente, sottraendoli
    sistematicamente dal dataset finale.
    """
    dati_esempio["ctrl"]["blank_values"] = ["Surface swab"]
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

    **Razionale Scientifico/Sistemistico**: Evita al ricercatore un frustrante
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

    **Razionale Scientifico/Sistemistico**: Quando la CLI viene invocata da
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

    **Razionale Scientifico/Sistemistico**: Gestisce con chiarezza il caso di un
    file di configurazione appena creato con ``touch`` o troncato per errore.
    """
    percorso = tmp_path / "vuoto.yaml"
    percorso.write_text("", encoding="utf-8")
    with pytest.raises(ErroreConfigurazione, match="vuoto"):
        carica(percorso)


def test_yaml_malformato_viene_segnalato(tmp_path):
    """
    **Obiettivo**: Verificare che un errore sintattico del parser YAML
    (``yaml.YAMLError`` su parentesi graffa non chiusa) venga incapsulato in
    ``ErroreConfigurazione``.

    **Razionale Scientifico/Sistemistico**: Uniforma tutte le patologie del file
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


def test_chiavi_derivate_dal_dataset_esistono_nello_schema():
    """
    **Obiettivo**: Verificare che tutte le chiavi elencate nel registro
    ``defaults.DERIVATI_DAL_DATASET`` corrispondano a campi reali dello schema.

    **Razionale Scientifico/Sistemistico**: ``DERIVATI_DAL_DATASET`` censisce i
    19 parametri i cui valori di default sono calibrati specificamente sulla
    chimica 515F/806R e sul disegno sperimentale di **OSD-734**: questo test
    impedisce che rinomine nello schema rendano orfano quel registro.
    """
    ignote = set(defaults.DERIVATI_DAL_DATASET) - chiavi_schema()
    assert not ignote, f"chiavi inesistenti nello schema: {sorted(ignote)}"


def test_esempio_usa_i_valori_del_dataset_di_riferimento(dati_esempio):
    """
    **Obiettivo**: Verificare che ``config/config.example.yaml`` riporti gli
    stessi valori predefiniti di ``defaults.py`` per le colonne batch, controllo,
    database SILVA e taxon bersaglio KatharoSeq di OSD-734.

    **Razionale Scientifico/Sistemistico**: Garantisce che il file di esempio
    resti sempre allineato alle costanti di riferimento di **OSD-734** senza
    derive silenziose tra codice e configurazione.
    """
    assert dati_esempio["err"]["batch_column"] == defaults.ERR_BATCH_COLUMN
    assert dati_esempio["decontam"]["batch_column"] == defaults.DECONTAM_BATCH_COLUMN
    assert dati_esempio["ctrl"]["column"] == defaults.CTRL_COLUMN
    assert dati_esempio["tax"]["ref_name"] == defaults.TAX_REF_NAME
    assert dati_esempio["tax"]["ref_version"] == defaults.TAX_REF_VERSION
    assert dati_esempio["katharoseq"]["target_taxon"] == defaults.KATHAROSEQ_TARGET_TAXON


def test_whitelist_dei_ritentativi_non_dipende_dal_dataset():
    """
    **Obiettivo**: Verificare che ``retry.whitelist`` NON compaia in
    ``defaults.DERIVATI_DAL_DATASET``.

    **Razionale Scientifico/Sistemistico**: La whitelist dei 4 codici ripetibili
    (``E-S2-03``, ``E-S3-01``, ``E-S4-02``, ``E-S5-01``) è un invariante
    architetturale della pipeline — ammette il retry automatico solo dove
    l'azione correttiva non altera alcuna assunzione biologica — e non deve mai
    essere considerata un parametro da ritarare cambiando dataset.
    """
    assert "retry.whitelist" not in defaults.DERIVATI_DAL_DATASET


def test_esempio_resta_allineato_alla_whitelist_predefinita(dati_esempio):
    """
    **Obiettivo**: Verificare che la lista ``retry.whitelist`` in
    ``config.example.yaml`` coincida con ``defaults.RETRY_WHITELIST``.

    **Razionale Scientifico/Sistemistico**: Garantisce che la configurazione di
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

    **Razionale Scientifico/Sistemistico**: Se una fase intermedia potesse
    modificare in-place l'oggetto ``Config`` condiviso in memoria, le fasi
    successive girerebbero con parametri diversi da quelli registrati in
    ``00_config/resolved.yaml``, invalidando il digest SHA-256 e la
    riproducibilità dell'intera corsa.
    """
    config = carica(ESEMPIO)
    with pytest.raises(Exception):
        config.run.threads = 1


# =========================================================================== #
# Gate G15 — coerenza interna della configurazione (W4)                        #
# =========================================================================== #

#: Numero esatto di campioni biologici del dataset di riferimento OSD-734 (su 960 totali).
CAMPIONI_BIOLOGICI_OSD734 = 803


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

    **Razionale Scientifico/Sistemistico**: Certifica che la configurazione di
    riferimento per OSD-734 soddisfa simultaneamente tutte le regole di coerenza
    interna (``E-G15-01`` .. ``E-G15-09``).
    """
    risolta = esegui_g15(dati_esempio)
    assert isinstance(risolta, ConfigRisolta)


def test_registro_dei_controlli_e_coerente():
    """
    **Obiettivo**: Verificare che il registro dichiarativo ``CONTROLLI`` del
    Gate G15 non contenga codici duplicati e dichiari per ciascuno l'origine
    (``"schema"`` o ``"gate"``) e la descrizione.

    **Razionale Scientifico/Sistemistico**: Garantisce la tracciabilità formale
    di ogni regola di validazione tra il modello Pydantic e il controllo G15.
    """
    codici = [c.codice for c in CONTROLLI]
    assert len(codici) == len(set(codici)), "codici duplicati"
    for controllo in CONTROLLI:
        assert controllo.implementato_da in ("schema", "gate")
        assert controllo.descrizione


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
        lambda d: d["ctrl"].__setitem__("blank_values", []),
        "E-G15-07",
        ("decontam.method", "ctrl.blank_values"),
        id="nessun-controllo-negativo",
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

    **Razionale Scientifico/Sistemistico**: Previene configurazioni numericamente
    valide sul singolo campo ma biologicamente o matematicamente degeneri
    (es. ``decontam.threshold = 1.0`` che classificherebbe come contaminante il
    100% delle ASV in S12, oppure ``ctrl.blank_values = []`` con metodo
    ``prevalence`` o ``combined`` che richiede obbligatoriamente controlli negativi).
    """
    errore = _viola(dati_esempio, modifica)

    codici = {v.codice for v in errore.violazioni}
    assert codice in codici, f"atteso {codice}, ottenuti {sorted(codici)}"

    testo = str(errore)
    for parametro in parametri:
        assert parametro in testo, f"il messaggio non nomina {parametro}: {testo}"


def test_g15_termina_con_un_codice_di_errore(dati_esempio):
    """
    **Obiettivo**: Verificare che ``decontam.threshold = 0.0`` (estremo inferiore
    escluso dall'intervallo aperto $(0, 1)$) faccia fallire G15 con codice
    ``E-G15-02`` e identificativo ``gate == "G15"``.

    **Razionale Scientifico/Sistemistico**: Nel test statistico di ``decontam``,
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

    **Razionale Scientifico/Sistemistico**: Permette all'utente di vedere e
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
    (``truncLen = 137``, ``trimLeft = 0``, ``len_tol = 0``) i tre parametri
    derivati statici assumano esattamente il valore ``137`` bp.
    **Razionale Scientifico/Sistemistico**: In una corsa single-end Illumina
    MiSeq troncata a 137 bp senza taglio in 5' (``trimLeft = 0``) e con tolleranza
    zero (``len_tol = 0``), ogni lettura e ogni ASV in uscita da DADA2 ha lunghezza
    attesa pari a 137 - 0 = 137 nt: derivare queste soglie per formula
    elimina il rischio che l'operatore aggiorni ``truncLen`` dimenticando di
    aggiornare ``filter.minLen``, ``asv.len_min`` e ``asv.len_max``.
    """
    derivati = esegui_g15(dati_esempio).derivati
    assert derivati.filter_minLen == 137
    assert derivati.asv_len_min == 137
    assert derivati.asv_len_max == 137


def test_derivati_seguono_i_parametri_da_cui_discendono(dati_esempio):
    r"""
    **Obiettivo**: Verificare che al variare di ``truncLen = 150``,
    ``trimLeft = 10`` e ``len_tol = 2``, i derivati vengano ricalcolati
    dinamicamente in ``filter_minLen = 150``, ``asv_len_min = 138`` e
    ``asv_len_max = 142``.

    **Razionale Scientifico/Sistemistico**: Dimostra che i parametri derivati
    non sono costanti cablate ma funzioni pure di ``truncLen``, ``trimLeft`` e
    ``len_tol`` ($L = \text{truncLen} - \text{trimLeft}$, intervallo ASV
    $[L - \text{len\_tol}, L + \text{len\_tol}]$).
    """
    dati_esempio["filter"]["truncLen"] = 150
    dati_esempio["filter"]["trimLeft"] = 10
    dati_esempio["asv"]["len_tol"] = 2

    derivati = esegui_g15(dati_esempio).derivati
    assert derivati.filter_minLen == 150
    assert derivati.asv_len_min == 138
    assert derivati.asv_len_max == 142


def test_la_prima_fase_e_utilizzabile_senza_il_gate(dati_esempio):
    """
    **Obiettivo**: Verificare che ``risolvi(valida(dati))`` calcoli i derivati
    statici lasciando ``prev_min_samples = None`` anche se invocata
    direttamente senza passare da ``esegui_g15``.

    **Razionale Scientifico/Sistemistico**: Mantiene disaccoppiato il motore di
    risoluzione algebrica (`resolve.py`) dal registro dei gate (`g01_g15.py`).
    """
    risolta = risolvi(valida(dati_esempio))
    assert risolta.derivati.filter_minLen == 137
    assert risolta.derivati.prev_min_samples is None


def test_min_samples_non_e_noto_prima_dei_metadati(dati_esempio):
    """
    **Obiettivo**: Verificare che all'uscita del Gate G15 (prima della lettura
    delle tabelle ISA-Tab in S0) ``prev_min_samples`` valga ``None`` e
    ``risolta.completa`` sia ``False``.

    **Razionale Scientifico/Sistemistico**: Il numero di campioni biologici
    autentici non può essere indovinato dalla sola configurazione YAML ma deve
    emergere dal crosswalk validato in S0 (escludendo i 92 controlli negativi
    e i 65 controlli positivi di OSD-734).
    """
    risolta = esegui_g15(dati_esempio)
    assert risolta.derivati.prev_min_samples is None
    assert not risolta.completa


def test_min_samples_calcolato_quando_il_dato_diventa_disponibile(dati_esempio):
    r"""
    **Obiettivo**: Verificare che fornendo il conteggio dei campioni biologici di
    OSD-734 (``803``) a ``con_campioni_biologici(803)``, ``prev_min_samples``
    diventi ``9`` (``ceil(0.01 * 803)``) e ``completa`` passi a ``True``.

    **Razionale Scientifico/Sistemistico**: Certifica il valore esatto della
    soglia di prevalenza applicata in S13 sul dataset reale OSD-734: una variante
    ASV deve comparire in almeno 9 campioni biologici su 803 ($\ge 1\%$) per
    superare il filtro di prevalenza.
    """
    completa = esegui_g15(dati_esempio).con_campioni_biologici(
        CAMPIONI_BIOLOGICI_OSD734
    )
    assert completa.derivati.prev_min_samples == 9
    assert completa.completa


def test_min_samples_arrotonda_per_eccesso(dati_esempio):
    r"""
    **Obiettivo**: Verificare che ``prev_min_samples`` venga calcolato tramite
    arrotondamento all'intero superiore (``math.ceil``): con ``min_fraction = 0.01``,
    100 campioni danno ``1`` mentre 101 campioni danno ``2``.

    **Razionale Scientifico/Sistemistico**: Arrotondare per difetto (``int()`` o
    ``round()``) su 101 campioni darebbe $101 \times 0.01 = 1.01 \to 1$ campione
    ($0.99\% < 1\%$), ammettendo nel dataset finale batteri e singleton presenti
    in una frazione di campioni inferiore alla soglia biologica minima stabilita
    dal protocollo e inquinando il filtraggio di prevalenza in S13.
    """
    risolta = esegui_g15(dati_esempio)
    assert risolta.con_campioni_biologici(100).derivati.prev_min_samples == 1
    assert risolta.con_campioni_biologici(101).derivati.prev_min_samples == 2
    assert risolta.con_campioni_biologici(0).derivati.prev_min_samples == 0


def test_numero_di_campioni_negativo_viene_respinto(dati_esempio):
    """
    **Obiettivo**: Verificare che ``con_campioni_biologici(-1)`` sollevi
    immediatamente ``ValueError``.

    **Razionale Scientifico/Sistemistico**: Impedisce che un contatore errato
    produca una soglia di prevalenza negativa che disabiliterebbe silenziosamente
    il filtro in S13.
    """
    with pytest.raises(ValueError, match="negativo"):
        esegui_g15(dati_esempio).con_campioni_biologici(-1)


@pytest.mark.parametrize("chiave", PARAMETRI_DERIVATI)
def test_derivato_impostato_a_mano_e_un_errore(chiave, dati_esempio):
    """
    **Obiettivo**: Verificare che la valorizzazione manuale nel file YAML di uno
    qualsiasi dei 4 parametri derivati (``filter.minLen``, ``asv.len_min``,
    ``asv.len_max``, ``prev.min_samples``) venga bloccata da G15 con codice
    ``E-G15-08``.

    **Razionale Scientifico/Sistemistico**: I 4 parametri derivati sono
    determinati univocamente dalle formule della pipeline; permettere all'utente
    di sovrascriverli a mano creerebbe contraddizioni interne tra il filtro
    ``filterAndTrim`` (S2), il filtro di lunghezza ASV (S7) e la soglia di
    prevalenza (S13).
    """
    gruppo, campo = chiave.split(".")
    dati_esempio[gruppo][campo] = 1

    with pytest.raises(ErroreGate) as errore:
        esegui_g15(dati_esempio)

    assert "E-G15-08" in {v.codice for v in errore.value.violazioni}
    assert chiave in str(errore.value)


def test_elenco_dei_derivati_coincide_con_quelli_calcolati(dati_esempio):
    """
    **Obiettivo**: Verificare che le chiavi prodotte da ``Derivati.come_chiavi()``
    coincidano esattamente con la costante ``PARAMETRI_DERIVATI``.

    **Razionale Scientifico/Sistemistico**: Impedisce disallineamenti tra
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

    **Razionale Scientifico/Sistemistico**: Il digest della configurazione è la
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

    **Razionale Scientifico/Sistemistico**: Cambiare il seme del generatore
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

    **Razionale Scientifico/Sistemistico**: Riformattare o riordinare i blocchi
    di un file YAML senza cambiare alcun valore non deve mai invalidare un'analisi
    già completata su disco.
    """
    originale = esegui_g15(copy.deepcopy(dati_esempio)).digest
    rovesciato = {k: dati_esempio[k] for k in reversed(list(dati_esempio))}
    assert esegui_g15(rovesciato).digest == originale


def test_il_digest_cambia_se_cambia_un_derivato(dati_esempio):
    """
    **Obiettivo**: Verificare che una variazione di ``asv.len_tol`` (che sposta
    ``asv.len_min`` e ``asv.len_max``) modifichi il digest finale.

    **Razionale Scientifico/Sistemistico**: Assicura che l'impronta crittografica
    copra tanto i parametri dichiarati quanto i valori derivati effettivamente
    applicati nel calcolo.
    """
    originale = esegui_g15(copy.deepcopy(dati_esempio)).digest
    dati_esempio["asv"]["len_tol"] = 3
    assert esegui_g15(dati_esempio).digest != originale


def test_il_digest_non_cambia_fra_le_due_fasi(dati_esempio):
    """
    **Obiettivo**: Verificare che l'arricchimento di ``ConfigRisolta`` con il
    numero di campioni biologici (``con_campioni_biologici(803)``) mantenga
    invariato ``risolta.digest``.

    **Razionale Scientifico/Sistemistico**: Il digest identifica univocamente la
    configurazione scelta dall'utente (che contiene già ``prev.min_fraction = 0.01``),
    evitando che il passaggio dalla fase pre-S0 alla fase post-S0 alteri
    l'identificativo della corsa registrato nei log iniziali.
    """
    risolta = esegui_g15(dati_esempio)
    completa = risolta.con_campioni_biologici(CAMPIONI_BIOLOGICI_OSD734)
    assert completa.digest == risolta.digest


def test_il_digest_ha_la_forma_attesa(dati_esempio):
    """
    **Obiettivo**: Verificare che ``risolta.digest`` rispetti il formato
    ``sha256:<64_caratteri_esadecimali_minuscoli>``.

    **Razionale Scientifico/Sistemistico**: Uniforma il formato del digest di
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

    **Razionale Scientifico/Sistemistico**: Persiste nella cartella ``00_config``
    il contratto completo dell'esecuzione, rendendo il risultato finale in
    ``12_final`` auto-documentato e riproducibile anche a distanza di anni.
    """
    risolta = esegui_g15(dati_esempio).con_campioni_biologici(
        CAMPIONI_BIOLOGICI_OSD734
    )
    percorso = scrivi_risolta(risolta, tmp_path)

    assert percorso == tmp_path / Fase.CONFIG.value / NOME_FILE_RISOLTO
    assert percorso.is_file()

    documento = yaml.safe_load(percorso.read_text(encoding="utf-8"))
    assert documento["digest"] == risolta.digest


def test_resolved_contiene_i_parametri_derivati(dati_esempio, tmp_path):
    """
    **Obiettivo**: Verificare che ``00_config/resolved.yaml`` contenga sia i
    valori calcolati dei 4 parametri derivati (``137``, ``137``, ``137``, ``9``)
    dentro ``parametri``, sia la mappa separata ``derivati``.

    **Razionale Scientifico/Sistemistico**: Consente a chi ispeziona
    ``resolved.yaml`` (e al report finale in ``12_final``) di leggere
    direttamente i valori numerici risolti senza doverli ricalcolare a mente.
    """
    risolta = esegui_g15(dati_esempio).con_campioni_biologici(
        CAMPIONI_BIOLOGICI_OSD734
    )
    documento = yaml.safe_load(
        scrivi_risolta(risolta, tmp_path).read_text(encoding="utf-8")
    )

    parametri = documento["parametri"]
    assert parametri["filter"]["minLen"] == 137
    assert parametri["asv"]["len_min"] == 137
    assert parametri["asv"]["len_max"] == 137
    assert parametri["prev"]["min_samples"] == 9
    assert set(documento["derivati"]) == set(PARAMETRI_DERIVATI)


def test_resolved_registra_il_derivato_non_ancora_noto(dati_esempio, tmp_path):
    """
    **Obiettivo**: Verificare che se ``scrivi_risolta`` viene chiamata prima del
    completamento di S0, ``parametri.prev.min_samples`` venga serializzato
    esplicitamente come ``null`` (``None``).

    **Razionale Scientifico/Sistemistico**: Documenta trasparentemente su disco
    che la risoluzione statica è avvenuta ma che il denominatore biologico non
    è stato ancora estratto dai metadati.
    """
    documento = yaml.safe_load(
        scrivi_risolta(esegui_g15(dati_esempio), tmp_path).read_text(encoding="utf-8")
    )
    assert documento["parametri"]["prev"]["min_samples"] is None


def test_resolved_registra_tutti_i_parametri_dichiarati(dati_esempio, tmp_path):
    """
    **Obiettivo**: Verificare che ``00_config/resolved.yaml`` includa l'unione
    completa di tutte le chiavi dello schema e di tutti i parametri derivati.

    **Razionale Scientifico/Sistemistico**: Anche se l'utente fornisce un file
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

    **Razionale Scientifico/Sistemistico**: Poiché ``scrivi_risolta`` registra
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

    **Razionale Scientifico/Sistemistico**: La scrittura di ``00_config/resolved.yaml``
    è spesso la primissima operazione su disco di una nuova corsa in una
    directory ``io.out_root`` non ancora esistente.
    """
    out_root = tmp_path / "non" / "ancora" / "esistente"
    percorso = scrivi_risolta(esegui_g15(dati_esempio), out_root)
    assert percorso.is_file()
