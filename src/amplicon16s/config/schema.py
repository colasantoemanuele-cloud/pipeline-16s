"""Schema dei parametri di configurazione.

La pipeline non contiene valori scritti nel codice: ogni scelta che dipende dal
dataset è un parametro. Il file di configurazione è quindi il vero punto di
controllo dell'esecuzione, e un errore al suo interno deve emergere prima che
venga allocato qualunque calcolo — non a metà di un'elaborazione che dura ore.

Lo schema è dichiarato con pydantic. Tre proprietà contano più delle altre:

* **Nessuna chiave sconosciuta.** Ogni gruppo rifiuta i campi che non conosce.
  Senza questo vincolo un refuso come ``maxEEE`` verrebbe ignorato in silenzio
  e l'esecuzione proseguirebbe con il valore predefinito: un errore che non si
  manifesta è peggio di uno che ferma il programma.
* **Vincoli di dominio, non solo di tipo.** Una frazione fuori da [0, 1] o un
  numero di thread negativo sono errori tanto quanto una stringa al posto di un
  numero.
* **Coerenza fra campi.** Alcuni controlli riguardano più parametri insieme e
  non si possono esprimere campo per campo.

I valori predefiniti stanno in :mod:`amplicon16s.config.defaults`.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Annotated, Any, Final, Literal

import yaml
from pydantic import (
    AfterValidator,
    BaseModel,
    ConfigDict,
    Field,
    StrictBool,
    ValidationError,
    model_validator,
)

from amplicon16s.config import defaults as d

__all__ = [
    "Config",
    "ErroreConfigurazione",
    "PARAMETRI_DERIVATI",
    "carica",
    "chiavi_schema",
    "gruppi_schema",
    "valida",
]


# --------------------------------------------------------------------------- #
# Tipi riutilizzabili                                                          #
# --------------------------------------------------------------------------- #


def _iupac_valido(valore: str) -> str:
    """Rifiuta una sequenza che contenga simboli non nucleotidici."""
    ammessi = set("ACGTRYSWKMBDHVN")
    estranei = sorted(set(valore.upper()) - ammessi)
    if estranei:
        raise ValueError(
            "simboli non ammessi nella sequenza: {}; sono ammessi i codici IUPAC "
            "{}".format(", ".join(estranei), "".join(sorted(ammessi)))
        )
    return valore.upper()


def _regex_valida(valore: str) -> str:
    """Rifiuta un'espressione regolare che non compila.

    Un'espressione malformata scoperta a metà esecuzione costerebbe l'intera
    elaborazione già svolta.
    """
    try:
        re.compile(valore)
    except re.error as errore:
        raise ValueError(f"espressione regolare non valida: {errore}") from errore
    return valore


Regex = Annotated[str, Field(min_length=1), AfterValidator(_regex_valida)]
StringaNonVuota = Annotated[str, Field(min_length=1)]
InteroPositivo = Annotated[int, Field(gt=0)]
InteroNonNegativo = Annotated[int, Field(ge=0)]
RealePositivo = Annotated[float, Field(gt=0)]
Frazione = Annotated[float, Field(ge=0.0, le=1.0)]
ElencoNonVuoto = Annotated[list[StringaNonVuota], Field(min_length=1)]

#: Riferimento a un'immagine container ancorata per digest. Un tag può essere
#: riassegnato a un'immagine diversa, un digest no: accettare un tag mobile
#: qui vanificherebbe la riproducibilità che l'immagine deve garantire.
RiferimentoImmagine = Annotated[
    str,
    Field(
        pattern=r"^\S+@sha256:[0-9a-f]{64}$",
        description="immagine ancorata per digest, es. registro/nome@sha256:<64 esadecimali>",
    ),
]

Md5 = Annotated[str, Field(pattern=r"^[0-9a-fA-F]{32}$")]

#: Sequenza nucleotidica in codici IUPAC, normalizzata in maiuscolo.
SequenzaIupac = Annotated[str, Field(min_length=1), AfterValidator(_iupac_valido)]


#: Parametri calcolati dalla pipeline a partire da altri, e per questo non
#: ammessi nel file di configurazione: scriverli a mano significherebbe poterli
#: mettere in contraddizione con i parametri da cui dipendono. Il calcolo vive
#: in :mod:`amplicon16s.config.resolve`; qui servono i soli nomi, per poter
#: spiegare il rifiuto invece di limitarsi a segnalare una chiave sconosciuta.
PARAMETRI_DERIVATI: Final[tuple[str, ...]] = (
    "filter.minLen",
    "asv.len_min",
    "asv.len_max",
    "prev.min_samples",
)


class _Gruppo(BaseModel):
    """Base comune a tutti i gruppi di parametri."""

    model_config = ConfigDict(
        extra="forbid",       # un refuso e' un errore, non un campo in piu'
        frozen=True,          # la configurazione non cambia durante l'esecuzione
        validate_default=True,
        str_strip_whitespace=True,
    )


# --------------------------------------------------------------------------- #
# Gruppi della pipeline di produzione                                          #
# --------------------------------------------------------------------------- #


class Io(_Gruppo):
    """Percorsi di ingresso e di uscita, e individuazione dei file di lettura.

    L'esistenza dei percorsi non viene controllata qui: lo schema è una
    proprietà del file, non della macchina su cui gira. I controlli sul
    filesystem appartengono alla fase di validazione degli input.
    """

    fastq_dir: Path
    assay_table: Path
    # Tabella campioni di studio: porta la colonna della classe del campione,
    # che nella tabella di assay non c'e'. E' condivisa fra piu' assay dello
    # stesso studio, quindi il join verso di essa va tenuto ristretto.
    study_table: Path
    out_root: Path
    # Associa a ciascun campione la piastra di estrazione e la corsa di
    # sequenziamento. E' facoltativo: quando manca, entrambe restano nulle e la
    # pipeline procede in modalita' a corsa singola e senza lotto.
    batch_table: Path | None = None
    fastq_glob: StringaNonVuota = d.IO_FASTQ_GLOB
    accession_regex: Regex = d.IO_ACCESSION_REGEX


class Meta(_Gruppo):
    """Lettura della tabella dei metadati e derivazione dei campi."""

    sample_id_column: StringaNonVuota = d.META_SAMPLE_ID_COLUMN
    accession_column: StringaNonVuota = d.META_ACCESSION_COLUMN
    derive_module: StrictBool = d.META_DERIVE_MODULE
    module_regex: Regex = d.META_MODULE_REGEX
    # Derivati dal dataset di riferimento.
    module_column: StringaNonVuota = d.META_MODULE_COLUMN
    # Posizioni che non sono superfici. Vale in entrambe le vie di
    # attribuzione del modulo: anche un modulo dichiarato dal file di
    # arricchimento non viene attribuito a un campione che non sta su una
    # superficie, altrimenti lo stesso dataset darebbe raggruppamenti diversi
    # a seconda che quel file ci sia o no.
    non_surface_positions: list[StringaNonVuota] = Field(
        default_factory=lambda: list(d.META_NON_SURFACE_POSITIONS)
    )
    # Colonna con cui il file facoltativo identifica il campione: deve
    # contenere l'accession, non il nome.
    batch_key_column: StringaNonVuota = d.META_BATCH_KEY_COLUMN
    # Colonna del file facoltativo che dichiara il modulo. Ha la precedenza
    # sulla derivazione da module_regex, che resta il ripiego quando il file
    # non c'e'.
    batch_module_column: StringaNonVuota = d.META_BATCH_MODULE_COLUMN


class Filter(_Gruppo):
    """Filtraggio e troncamento delle letture grezze."""

    truncLen: InteroPositivo = d.FILTER_TRUNCLEN
    truncLen_shortfall_warn: InteroNonNegativo = d.FILTER_TRUNCLEN_SHORTFALL_WARN
    trimLeft: InteroNonNegativo = d.FILTER_TRIMLEFT
    maxEE: RealePositivo = d.FILTER_MAXEE
    truncQ: InteroNonNegativo = d.FILTER_TRUNCQ
    maxN: InteroNonNegativo = d.FILTER_MAXN
    rm_phix: StrictBool = d.FILTER_RM_PHIX

    @model_validator(mode="after")
    def _tolleranza_minore_del_troncamento(self) -> Filter:
        if self.truncLen_shortfall_warn >= self.truncLen:
            raise ValueError(
                "truncLen_shortfall_warn ({}) deve essere minore di truncLen ({}): "
                "la tolleranza si confronta con lo scarto fra il troncamento e la "
                "lettura piu' corta osservata, e una tolleranza grande quanto il "
                "troncamento non potrebbe mai essere superata".format(
                    self.truncLen_shortfall_warn, self.truncLen
                )
            )
        return self


class Err(_Gruppo):
    """Apprendimento del modello di errore di sequenziamento."""

    nbases: RealePositivo = d.ERR_NBASES
    max_consist: InteroPositivo = d.ERR_MAX_CONSIST
    randomize: StrictBool = d.ERR_RANDOMIZE
    # Derivato dal dataset di riferimento.
    batch_column: StringaNonVuota = d.ERR_BATCH_COLUMN


class Dada(_Gruppo):
    """Inferenza delle varianti di sequenza."""

    # dada2 accetta TRUE, FALSE oppure "pseudo": il parametro e' passato cosi'
    # com'e' ed eredita il vocabolario della libreria.
    pool: StrictBool | Literal["pseudo"] = d.DADA_POOL
    omega_a: RealePositivo = d.DADA_OMEGA_A


class Chimera(_Gruppo):
    """Rimozione delle sequenze chimeriche."""

    # Vocabolario di dada2::removeBimeraDenovo.
    method: Literal["consensus", "pooled", "per-sample"] = d.CHIMERA_METHOD
    min_fold_parent_over_abundance: Annotated[float, Field(ge=1.0)] = (
        d.CHIMERA_MIN_FOLD_PARENT_OVER_ABUNDANCE
    )
    min_parent_abundance: InteroPositivo = d.CHIMERA_MIN_PARENT_ABUNDANCE
    min_sample_fraction: Frazione = d.CHIMERA_MIN_SAMPLE_FRACTION
    allow_one_off: StrictBool = d.CHIMERA_ALLOW_ONE_OFF


class Asv(_Gruppo):
    """Selezione delle sequenze varianti."""

    len_tol: InteroNonNegativo = d.ASV_LEN_TOL


class Tax(_Gruppo):
    """Assegnazione tassonomica e riferimento impiegato.

    Nome, versione e checksum del riferimento sono obbligatori: senza di essi
    un risultato non è riconducibile ai dati che l'hanno prodotto.
    """

    ref_fasta: Path
    ref_md5: Md5
    # Obbligatori e senza valore predefinito: ereditare in silenzio un
    # riferimento tassonomico renderebbe il risultato non riconducibile ai
    # dati che l'hanno prodotto. I valori dello studio di riferimento sono
    # documentati in defaults.py e compaiono in config.example.yaml.
    ref_name: StringaNonVuota
    ref_version: StringaNonVuota
    # Insieme chiuso: va esteso quando un altro classificatore viene realizzato.
    classifier: Literal["naive_bayes"] = d.TAX_CLASSIFIER
    min_boot: Annotated[int, Field(ge=0, le=100)] = d.TAX_MIN_BOOT
    try_rc: StrictBool = d.TAX_TRY_RC
    assign_species: StrictBool = d.TAX_ASSIGN_SPECIES


class Filt(_Gruppo):
    """Filtraggio tassonomico successivo all'assegnazione."""

    remove_na_phylum: StrictBool = d.FILT_REMOVE_NA_PHYLUM
    exclude_taxa: list[StringaNonVuota] = Field(
        default_factory=lambda: list(d.FILT_EXCLUDE_TAXA)
    )


class Phylo(_Gruppo):
    """Costruzione dell'albero filogenetico."""

    enabled: StrictBool = d.PHYLO_ENABLED
    max_seqs: InteroPositivo = d.PHYLO_MAX_SEQS


class Ctrl(_Gruppo):
    """Riconoscimento di controlli negativi, positivi e campioni biologici.

    Tutti i valori sono derivati dal dataset di riferimento: sono le etichette
    usate in quello studio, non una proprietà della pipeline.
    """

    column: StringaNonVuota = d.CTRL_COLUMN
    blank_values: ElencoNonVuoto = Field(
        default_factory=lambda: list(d.CTRL_BLANK_VALUES)
    )
    positive_values: ElencoNonVuoto = Field(
        default_factory=lambda: list(d.CTRL_POSITIVE_VALUES)
    )
    biological_values: ElencoNonVuoto = Field(
        default_factory=lambda: list(d.CTRL_BIOLOGICAL_VALUES)
    )

    @model_validator(mode="after")
    def _categorie_disgiunte(self) -> Ctrl:
        categorie = {
            "blank_values": set(self.blank_values),
            "positive_values": set(self.positive_values),
            "biological_values": set(self.biological_values),
        }
        nomi = sorted(categorie)
        for i, primo in enumerate(nomi):
            for secondo in nomi[i + 1 :]:
                comuni = categorie[primo] & categorie[secondo]
                if comuni:
                    raise ValueError(
                        "le etichette {} compaiono sia in {} sia in {}: un campione "
                        "non puo' appartenere a due categorie".format(
                            sorted(comuni), primo, secondo
                        )
                    )
        return self


class Katharoseq(_Gruppo):
    """Calibrazione KatharoSeq sui controlli positivi."""

    # Derivato dal dataset di riferimento.
    target_taxon: StringaNonVuota = d.KATHAROSEQ_TARGET_TAXON


class Decontam(_Gruppo):
    """Rimozione dei contaminanti a partire dai controlli negativi."""

    # Vocabolario di decontam::isContaminant.
    method: Literal[
        "auto", "frequency", "prevalence", "combined", "minimum", "either", "both"
    ] = d.DECONTAM_METHOD
    threshold: Frazione = d.DECONTAM_THRESHOLD
    min_blanks: InteroPositivo = d.DECONTAM_MIN_BLANKS
    # Derivato dal dataset di riferimento.
    batch_column: StringaNonVuota = d.DECONTAM_BATCH_COLUMN


class Prev(_Gruppo):
    """Filtro di prevalenza sulle varianti."""

    min_fraction: Frazione = d.PREV_MIN_FRACTION
    min_count: InteroNonNegativo = d.PREV_MIN_COUNT
    apply: StrictBool = d.PREV_APPLY


class Qc(_Gruppo):
    """Soglie dei controlli di qualità che governano avvisi e arresti."""

    min_reads_raw: InteroNonNegativo = d.QC_MIN_READS_RAW
    min_reads_filtered: InteroNonNegativo = d.QC_MIN_READS_FILTERED
    min_reads_final: InteroNonNegativo = d.QC_MIN_READS_FINAL
    max_asv_count: InteroPositivo = d.QC_MAX_ASV_COUNT
    warn_frac_chimeric: Frazione = d.QC_WARN_FRAC_CHIMERIC
    stop_frac_chimeric: Frazione = d.QC_STOP_FRAC_CHIMERIC
    # Letture ispezionate per file dai gate che leggono le sequenze.
    head_reads: InteroPositivo = d.QC_HEAD_READS
    max_primer_hit_frac: Frazione = d.QC_MAX_PRIMER_HIT_FRAC
    min_motif_frac: Frazione = d.QC_MIN_MOTIF_FRAC
    # Derivati dal dataset di riferimento: dipendono dalla regione amplificata.
    primer_sequence: SequenzaIupac = d.QC_PRIMER_SEQUENCE
    conserved_motif: Regex = d.QC_CONSERVED_MOTIF

    @model_validator(mode="after")
    def _avviso_prima_dell_arresto(self) -> Qc:
        if self.warn_frac_chimeric > self.stop_frac_chimeric:
            raise ValueError(
                "warn_frac_chimeric ({}) non puo' superare stop_frac_chimeric ({}): "
                "l'esecuzione si fermerebbe prima di emettere l'avviso".format(
                    self.warn_frac_chimeric, self.stop_frac_chimeric
                )
            )
        return self


class Retry(_Gruppo):
    """Nuovi tentativi sulle fasi fallite."""

    enabled: StrictBool = d.RETRY_ENABLED
    max_attempts: InteroPositivo = d.RETRY_MAX_ATTEMPTS
    # Codici di errore per cui e' ammesso un nuovo tentativo automatico.
    # L'elenco e' chiuso per scelta metodologica e non per una proprieta'
    # del dataset: si ritenta solo dove l'azione correttiva non modifica
    # alcuna assunzione dell'analisi.
    whitelist: list[StringaNonVuota] = Field(
        default_factory=lambda: list(d.RETRY_WHITELIST)
    )


class Out(_Gruppo):
    """Forma degli artefatti prodotti."""

    # Insiemi chiusi: vanno estesi quando un'altra forma viene realizzata.
    serialization: Literal["rds"] = d.OUT_SERIALIZATION
    asv_id_scheme: Literal["abundance_rank"] = d.OUT_ASV_ID_SCHEME
    taxa_are_rows: StrictBool = d.OUT_TAXA_ARE_ROWS
    export_flat: StrictBool = d.OUT_EXPORT_FLAT


class Run(_Gruppo):
    """Parametri dell'esecuzione: ambiente, parallelismo, riproducibilità."""

    container: RiferimentoImmagine
    # Il file di blocco delle versioni R. E' referenziato dal Dockerfile e
    # dagli script R, quindi l'ambiente sarebbe riproducibile comunque; entra
    # qui perche' cosi' finisce nella configurazione risolta e nel suo digest,
    # e due esecuzioni diventano confrontabili anche rispetto alle versioni R.
    lockfile: StringaNonVuota = d.RUN_LOCKFILE
    seed: int = d.RUN_SEED
    threads: InteroPositivo = d.RUN_THREADS
    batch_size: InteroPositivo = d.RUN_BATCH_SIZE


# --------------------------------------------------------------------------- #
# Gruppi delle analisi ecologiche a valle                                      #
# --------------------------------------------------------------------------- #
# I nomi dei gruppi sono riservati fin d'ora, ma le analisi ecologiche non sono
# ancora realizzate e nessun parametro e' stato definito. I gruppi restano
# quindi vuoti: essendo a chiave chiusa, qualunque parametro vi venga scritto
# oggi verrebbe respinto, ed e' il comportamento corretto finche' non esiste
# codice che sappia interpretarlo.


class Norm(_Gruppo):
    """Normalizzazione delle abbondanze. Nessun parametro ancora definito."""


class Glom(_Gruppo):
    """Aggregazione a un livello tassonomico. Nessun parametro ancora definito."""


class Beta(_Gruppo):
    """Misure di diversità beta. Nessun parametro ancora definito."""


class Ord(_Gruppo):
    """Ordinamento e riduzione dimensionale. Nessun parametro ancora definito."""


class Stat(_Gruppo):
    """Test statistici. Nessun parametro ancora definito."""


# --------------------------------------------------------------------------- #
# Radice                                                                       #
# --------------------------------------------------------------------------- #


class Config(_Gruppo):
    """Configurazione completa: tutti i gruppi, tutti i parametri."""

    # Pipeline di produzione.
    io: Io
    meta: Meta = Field(default_factory=Meta)
    filter: Filter = Field(default_factory=Filter)
    err: Err = Field(default_factory=Err)
    dada: Dada = Field(default_factory=Dada)
    chimera: Chimera = Field(default_factory=Chimera)
    asv: Asv = Field(default_factory=Asv)
    tax: Tax
    filt: Filt = Field(default_factory=Filt)
    phylo: Phylo = Field(default_factory=Phylo)
    ctrl: Ctrl = Field(default_factory=Ctrl)
    katharoseq: Katharoseq = Field(default_factory=Katharoseq)
    decontam: Decontam = Field(default_factory=Decontam)
    prev: Prev = Field(default_factory=Prev)
    qc: Qc = Field(default_factory=Qc)
    retry: Retry = Field(default_factory=Retry)
    out: Out = Field(default_factory=Out)
    run: Run

    # Analisi ecologiche a valle.
    norm: Norm = Field(default_factory=Norm)
    glom: Glom = Field(default_factory=Glom)
    beta: Beta = Field(default_factory=Beta)
    ord: Ord = Field(default_factory=Ord)
    stat: Stat = Field(default_factory=Stat)


# --------------------------------------------------------------------------- #
# Ispezione dello schema                                                       #
# --------------------------------------------------------------------------- #


def gruppi_schema() -> tuple[str, ...]:
    """Nomi dei gruppi dichiarati, nell'ordine in cui compaiono nello schema."""
    return tuple(Config.model_fields)


def chiavi_schema() -> frozenset[str]:
    """Tutte le chiavi dello schema nella forma ``gruppo.parametro``.

    Serve a confrontare lo schema con un file di configurazione senza
    riscrivere a mano l'elenco dei parametri, che si disallineerebbe al primo
    cambiamento.
    """
    chiavi: set[str] = set()
    for nome_gruppo, campo_gruppo in Config.model_fields.items():
        classe = campo_gruppo.annotation
        for nome_campo in classe.model_fields:
            chiavi.add(f"{nome_gruppo}.{nome_campo}")
    return frozenset(chiavi)


# --------------------------------------------------------------------------- #
# Caricamento e validazione                                                    #
# --------------------------------------------------------------------------- #


class ErroreConfigurazione(Exception):
    """Configurazione non valida.

    Raccoglie tutti i problemi trovati invece di fermarsi al primo: chi corregge
    un file di configurazione vuole l'elenco completo, non un errore alla volta.
    """

    def __init__(self, problemi: list[str], origine: str | None = None) -> None:
        self.problemi = problemi
        self.origine = origine
        intestazione = (
            f"configurazione non valida ({origine})"
            if origine
            else "configurazione non valida"
        )
        super().__init__(
            "{}: {} problem{}\n{}".format(
                intestazione,
                len(problemi),
                "a" if len(problemi) == 1 else "i",
                "\n".join(f"  - {p}" for p in problemi),
            )
        )


#: Messaggi per i tipi di errore che non dipendono dal contesto.
_MESSAGGI_SEMPLICI: Final[dict[str, str]] = {
    "missing": "parametro obbligatorio mancante",
    "extra_forbidden": "parametro sconosciuto (refuso o parametro non previsto)",
    "int_type": "deve essere un numero intero",
    "int_parsing": "deve essere un numero intero",
    "int_from_float": "deve essere un numero intero, non decimale",
    "float_type": "deve essere un numero",
    "float_parsing": "deve essere un numero",
    "string_type": "deve essere una stringa",
    "bool_type": "deve essere true oppure false",
    "bool_parsing": "deve essere true oppure false",
    "list_type": "deve essere un elenco",
    "dict_type": "deve essere una mappa di parametri",
    "model_type": "deve essere una mappa di parametri",
    "path_type": "deve essere un percorso",
}


def _spiega(errore: dict[str, Any]) -> tuple[str, bool]:
    """Restituisce la spiegazione dell'errore e se convenga mostrare l'input.

    I messaggi di pydantic sono in inglese: quelli ricorrenti vengono tradotti,
    per gli altri si ricade sul testo originale invece di inventare una
    traduzione approssimativa.
    """
    tipo = errore["type"]
    contesto = errore.get("ctx") or {}

    if tipo in _MESSAGGI_SEMPLICI:
        # Per un parametro mancante o sconosciuto il valore ricevuto non
        # aggiunge nulla: nel primo caso non esiste, nel secondo e' gia' noto.
        mostra_input = tipo not in ("missing", "extra_forbidden")
        return _MESSAGGI_SEMPLICI[tipo], mostra_input

    if tipo == "value_error":
        # Errore sollevato da un validatore nostro: il messaggio e' gia'
        # scritto per essere letto, e l'input sarebbe l'intero gruppo.
        return str(contesto.get("error", errore["msg"])), False

    secondo_contesto = {
        "greater_than": "deve essere maggiore di {gt}",
        "greater_than_equal": "deve essere maggiore o uguale a {ge}",
        "less_than": "deve essere minore di {lt}",
        "less_than_equal": "deve essere minore o uguale a {le}",
        "string_pattern_mismatch": "non rispetta il formato richiesto ({pattern})",
        "too_short": "deve contenere almeno {min_length} elementi",
        "string_too_short": "deve contenere almeno {min_length} caratteri",
    }
    if tipo in secondo_contesto:
        try:
            return secondo_contesto[tipo].format(**contesto), True
        except KeyError:  # contesto inatteso: meglio il messaggio originale
            pass

    if tipo == "literal_error":
        return f"valore non ammesso; ammessi: {contesto.get('expected', '?')}", True

    return errore["msg"], True


def _descrivi(errore: dict[str, Any]) -> str:
    """Traduce un errore di pydantic in una riga leggibile con chiave puntata."""
    chiave = ".".join(str(p) for p in errore["loc"]) or "<radice>"

    if errore["type"] == "extra_forbidden" and chiave in PARAMETRI_DERIVATI:
        return (
            f"{chiave}: parametro derivato, calcolato dalla pipeline a partire "
            f"dagli altri parametri; non va impostato nella configurazione"
        )

    spiegazione, mostra_input = _spiega(errore)

    if not mostra_input:
        return f"{chiave}: {spiegazione}"

    ricevuto = repr(errore.get("input"))
    if len(ricevuto) > 60:  # un gruppo intero renderebbe illeggibile la riga
        ricevuto = ricevuto[:57] + "..."
    return f"{chiave}: {spiegazione} (ricevuto: {ricevuto})"


def valida(dati: dict[str, Any], origine: str | None = None) -> Config:
    """Valida un dizionario già caricato.

    Solleva :class:`ErroreConfigurazione` con l'elenco completo dei problemi.
    """
    try:
        return Config.model_validate(dati)
    except ValidationError as errore:
        problemi = [_descrivi(e) for e in errore.errors()]
        raise ErroreConfigurazione(problemi, origine) from errore


def carica(percorso: str | Path) -> Config:
    """Carica e valida un file di configurazione YAML."""
    percorso = Path(percorso)
    testo = percorso.read_text(encoding="utf-8")

    try:
        dati = yaml.safe_load(testo)
    except yaml.YAMLError as errore:
        raise ErroreConfigurazione([f"YAML non leggibile: {errore}"], str(percorso)) from errore

    if dati is None:
        raise ErroreConfigurazione(["il file è vuoto"], str(percorso))
    if not isinstance(dati, dict):
        raise ErroreConfigurazione(
            [f"il file deve contenere una mappa di gruppi, trovato {type(dati).__name__}"],
            str(percorso),
        )

    return valida(dati, origine=str(percorso))
