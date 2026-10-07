"""Schema dei parametri di configurazione.

La pipeline non contiene valori scritti nel codice: ogni scelta che dipende dal
dataset è un parametro. Il file di configurazione è quindi il vero punto di
controllo dell'esecuzione, e un errore al suo interno deve emergere prima che
venga allocato qualunque calcolo, non a metà di un'elaborazione che dura ore.

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

import os
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
    "INTESTAZIONE_MANCANTI",
    "PARAMETRI_DERIVATI",
    "obbligatori_mancanti",
    "processori_disponibili",
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
)


def processori_disponibili() -> int:
    """I processori utilizzabili dal processo: quelli della sua affinita', che in
    un container limitato sono meno di quelli della macchina.
    """
    try:
        return len(os.sched_getaffinity(0))
    except AttributeError:  # piattaforme senza affinita' di processo
        return os.cpu_count() or 1


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
    # Facoltativa: senza, classe e variabili dei campioni si leggono dalla
    # tabella di assay, che allora deve portarne le colonne.
    study_table: Path | None = None
    out_root: Path
    # Associa a ciascun campione la piastra di estrazione e la corsa di
    # sequenziamento. E' facoltativo: quando manca, entrambe restano nulle e la
    # pipeline procede in modalita' a corsa singola e senza lotto.
    batch_table: Path | None = None
    fastq_glob: StringaNonVuota = d.IO_FASTQ_GLOB
    # Obbligatorio: estrae dal nome di ogni file la chiave del campione, che e'
    # la corrispondenza intera oppure, se l'espressione ha un gruppo di
    # cattura, il primo gruppo. La stessa chiave si ricava dal valore di
    # meta.accession_column.
    accession_regex: Regex


class Meta(_Gruppo):
    """Lettura della tabella dei metadati e derivazione dei campi."""

    # I parametri senza predefinito sono obbligatori (defaults.OBBLIGATORI):
    # descrivono il formato dei metadati. Quelli che ammettono il valore nullo
    # si dichiarano nulli quando il dataset non ha l'informazione.
    sample_id_column: StringaNonVuota
    accession_column: StringaNonVuota
    # Colonna del nome del campione nella tabella di studio, quando ha un nome
    # diverso da quello della tabella di assay; nullo: lo stesso nome.
    study_sample_id_column: StringaNonVuota | None = None
    derive_module: StrictBool = d.META_DERIVE_MODULE
    module_regex: Regex | None
    module_column: StringaNonVuota | None
    # Posizioni che non sono superfici. Vale in entrambe le vie di
    # attribuzione del modulo: anche un modulo dichiarato dal file di
    # arricchimento non viene attribuito a un campione che non sta su una
    # superficie, altrimenti lo stesso dataset darebbe raggruppamenti diversi
    # a seconda che quel file ci sia o no.
    non_surface_positions: list[StringaNonVuota]
    # Colonna con cui il file facoltativo identifica il campione: deve
    # contenere l'accession, non il nome.
    batch_key_column: StringaNonVuota | None
    # Colonna del file facoltativo che dichiara il modulo. Ha la precedenza
    # sulla derivazione da module_regex, che resta il ripiego quando il file
    # non c'e'.
    batch_module_column: StringaNonVuota | None

    @property
    def colonna_id_studio(self) -> str:
        """La colonna del nome del campione nella tabella di studio."""
        return self.study_sample_id_column or self.sample_id_column


class Filter(_Gruppo):
    """Filtraggio e troncamento delle letture grezze."""

    truncLen: InteroPositivo  # obbligatorio: dipende dalle letture del dataset
    truncLen_shortfall_warn: InteroNonNegativo = d.FILTER_TRUNCLEN_SHORTFALL_WARN
    trimLeft: InteroNonNegativo = d.FILTER_TRIMLEFT
    maxEE: RealePositivo = d.FILTER_MAXEE
    truncQ: InteroNonNegativo = d.FILTER_TRUNCQ
    maxN: InteroNonNegativo = d.FILTER_MAXN
    rm_phix: StrictBool = d.FILTER_RM_PHIX

    @model_validator(mode="after")
    def _tolleranza_minore_del_troncamento(self) -> Filter:
        """Respinge una tolleranza sullo scarto non minore di ``truncLen``, che non
        potrebbe mai essere superata.
        """
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
    # Obbligatorio, anche nullo. Colonna di io.batch_table con la
    # corsa di sequenziamento: un modello d'errore per corsa. Con null, o senza
    # io.batch_table, si stima un solo modello su tutti i campioni.
    batch_column: StringaNonVuota | None


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
    # Elenco dei taxa con un difetto noto del riferimento, nella forma
    # rank,name,path,ranks_present,ranks_expected: S8 marca le assegnazioni che
    # vi ricadono. Facoltativo: non tutti i riferimenti ne hanno uno.
    ref_bad_taxa: Path | None = d.TAX_REF_BAD_TAXA


class Filt(_Gruppo):
    """Filtraggio tassonomico successivo all'assegnazione."""

    remove_na_phylum: StrictBool = d.FILT_REMOVE_NA_PHYLUM
    exclude_taxa: list[StringaNonVuota] = Field(
        default_factory=lambda: list(d.FILT_EXCLUDE_TAXA)
    )


class Phylo(_Gruppo):
    """Costruzione dell'albero filogenetico."""

    enabled: StrictBool = d.PHYLO_ENABLED
    # Varianti dell'oggetto filtrato (S13) oltre le quali S9 si ferma prima del
    # calcolo (E-S9-01).
    max_seqs: InteroPositivo = d.PHYLO_MAX_SEQS
    # Insiemi chiusi: vanno estesi quando un'altra forma viene realizzata.
    # L'allineatore multiplo (DECIPHER::AlignSeqs) e il modello evolutivo
    # dell'albero di massima verosimiglianza (phangorn): GTR con eterogeneita'
    # gamma fra i siti e una quota di siti invarianti.
    aligner: Literal["decipher"] = d.PHYLO_ALIGNER
    model: Literal["GTR+G+I"] = d.PHYLO_MODEL


class Ctrl(_Gruppo):
    """Riconoscimento di controlli negativi, positivi e campioni biologici, e
    validazione della corsa dai controlli positivi.

    La colonna, le etichette e la regola di riclassificazione sono derivate dal
    dataset di riferimento: sono quelle usate in quello studio, non una
    proprietà della pipeline.
    """

    # Obbligatori. Le etichette dei controlli possono essere elenchi vuoti (un
    # dataset senza controlli positivi o negativi lo dichiara cosi'); quelle dei
    # campioni biologici no.
    column: StringaNonVuota
    blank_values: list[StringaNonVuota]
    positive_values: list[StringaNonVuota]
    biological_values: ElencoNonVuoto
    # Riclassificazione in controllo negativo, indipendente dal materiale
    # dichiarato: i campioni il cui valore nella colonna indicata (della tabella
    # campioni di studio) e' fra quelli elencati. Il materiale resta quello
    # originale nell'inventario, cosi' la riclassificazione resta tracciabile.
    # Con l'elenco vuoto non si riclassifica nulla.
    blank_override_column: StringaNonVuota | None
    blank_override_values: list[StringaNonVuota]
    # Validazione della corsa dai controlli positivi (S11).
    min_positives: InteroPositivo = d.CTRL_MIN_POSITIVES
    min_positive_pass_frac: Frazione = d.CTRL_MIN_POSITIVE_PASS_FRAC
    positive_gate: StrictBool = d.CTRL_POSITIVE_GATE

    @model_validator(mode="after")
    def _categorie_disgiunte(self) -> Ctrl:
        """Respinge un'etichetta assegnata a più di una classe di campioni.

        Il confronto ignora maiuscole e spazi ai bordi, come la classificazione
        dei campioni: "Blank" e "blank" sono la stessa etichetta.
        """
        categorie = {
            "blank_values": {v.strip().casefold() for v in self.blank_values},
            "positive_values": {v.strip().casefold() for v in self.positive_values},
            "biological_values": {v.strip().casefold() for v in self.biological_values},
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
        if self.blank_override_values and self.blank_override_column is None:
            raise ValueError(
                "blank_override_values e' indicato ma blank_override_column e' nullo: "
                "senza colonna la regola non riconosce alcun campione"
            )
        return self


class Katharoseq(_Gruppo):
    """Calibrazione KatharoSeq sui controlli positivi."""

    # Obbligatori; nulli se il dataset non ha controlli positivi (G15 verifica
    # che non lo siano quando ctrl.positive_values non e' vuoto).
    target_taxon: StringaNonVuota | None
    # Colonna del file di arricchimento con le cellule di ciascun controllo
    # positivo, il livello di diluizione; nome originale della colonna.
    cell_count_column: StringaNonVuota | None
    collapse_rank: Literal["Phylum", "Class", "Order", "Family", "Genus"] = (
        d.KATHAROSEQ_COLLAPSE_RANK
    )
    # Insiemi chiusi: vanno estesi quando un'altra forma viene realizzata.
    curve_model: Literal["allosteric_sigmoid"] = d.KATHAROSEQ_CURVE_MODEL
    # Strettamente fra 0 e 1: a fedelta' 1 la profondita' richiesta e' infinita.
    target_sensitivity: Annotated[float, Field(gt=0.0, lt=1.0)] = (
        d.KATHAROSEQ_TARGET_SENSITIVITY
    )
    min_r2: Frazione = d.KATHAROSEQ_MIN_R2
    read_stage: Literal["nonchimeric"] = d.KATHAROSEQ_READ_STAGE


class Decontam(_Gruppo):
    """Rimozione dei contaminanti a partire dai controlli negativi."""

    # Vocabolario di decontam::isContaminant.
    method: Literal[
        "auto", "frequency", "prevalence", "combined", "minimum", "either", "both"
    ] = d.DECONTAM_METHOD
    threshold: Frazione = d.DECONTAM_THRESHOLD
    min_blanks: InteroPositivo = d.DECONTAM_MIN_BLANKS
    # batch.combine di decontam::isContaminant, per la decontaminazione per piastra.
    batch_combine: Literal["minimum", "product", "fisher"] = d.DECONTAM_BATCH_COMBINE
    # Modalita' che decide i contaminanti rimossi: aggregate (tutti i biologici
    # contro tutti i negativi) o batch (per piastra, con batch_combine). L'altra
    # si calcola come diagnostica.
    mode: Literal["aggregate", "batch"] = d.DECONTAM_MODE
    # Obbligatorio, anche nullo: la colonna di io.batch_table con la piastra.
    batch_column: StringaNonVuota | None


class Prev(_Gruppo):
    """Filtro di prevalenza sulle varianti."""

    min_fraction: Frazione = d.PREV_MIN_FRACTION
    min_count: InteroNonNegativo = d.PREV_MIN_COUNT
    apply: StrictBool = d.PREV_APPLY


class Qc(_Gruppo):
    """Soglie dei controlli di qualità che governano avvisi e arresti."""

    min_reads_raw: InteroNonNegativo = d.QC_MIN_READS_RAW
    min_reads_mode: Literal["fixed", "katharoseq_if_available"] = d.QC_MIN_READS_MODE
    # Frazione massima delle letture dei biologici rimossa come contaminante (S12).
    max_frac_contaminant: Frazione = d.QC_MAX_FRAC_CONTAMINANT
    # Per campione, sulle letture dell'oggetto finale (S13): sotto, il campione
    # esce dall'oggetto finale.
    min_reads_final: InteroNonNegativo = d.QC_MIN_READS_FINAL
    # Per l'insieme dei campioni finali: letture finali su letture senza chimere
    # (S14, E-S14-01).
    min_frac_reads_retained: Frazione = d.QC_MIN_FRAC_READS_RETAINED
    max_asv_count: InteroPositivo = d.QC_MAX_ASV_COUNT
    # Controlli sul risultato del filtro (S2), applicati ai campioni biologici
    # e ai controlli positivi, non ai negativi: vedi steps/s02_filter.py.
    max_zeroed_samples: InteroNonNegativo = d.QC_MAX_ZEROED_SAMPLES
    max_frac_lost_filter: Frazione = d.QC_MAX_FRAC_LOST_FILTER
    warn_frac_chimeric: Frazione = d.QC_WARN_FRAC_CHIMERIC
    stop_frac_chimeric: Frazione = d.QC_STOP_FRAC_CHIMERIC
    # Frazione minima di varianti con il phylum assegnato (S8, E-S8-02): vedi
    # steps/s08_taxonomy.py per le classi a cui si applica.
    min_frac_phylum: Frazione = d.QC_MIN_FRAC_PHYLUM
    # Letture ispezionate per file dai gate che leggono le sequenze.
    head_reads: InteroPositivo = d.QC_HEAD_READS
    max_primer_hit_frac: Frazione = d.QC_MAX_PRIMER_HIT_FRAC
    min_motif_frac: Frazione = d.QC_MIN_MOTIF_FRAC
    # Obbligatori: dipendono dalla regione amplificata. Il motivo puo' essere
    # nullo, per una regione senza un motivo noto: G10 verifica allora la sola
    # assenza del primer.
    primer_sequence: SequenzaIupac
    conserved_motif: Regex | None

    @model_validator(mode="after")
    def _avviso_prima_dell_arresto(self) -> Qc:
        """Respinge una soglia d'avviso sulle chimere superiore a quella d'arresto."""
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
    # Da dove vengono gli identificativi dei campioni nell'oggetto integrato.
    sample_id_source: Literal["accession"] = d.OUT_SAMPLE_ID_SOURCE
    # Colonne dei metadati portate nell'oggetto integrato oltre a quelle
    # dell'inventario, con il loro nome originale: dalla tabella campioni di
    # studio e dal file di arricchimento. Nell'oggetto prendono un nome
    # sintattico, e la corrispondenza resta in 10_phyloseq/colonne_metadati.tsv.
    # Obbligatori; elenchi vuoti se non si porta alcuna colonna.
    study_columns: list[StringaNonVuota]
    batch_columns: list[StringaNonVuota]


class Run(_Gruppo):
    """Parametri dell'esecuzione: ambiente, parallelismo, riproducibilità."""

    # Facoltativo: l'immagine con cui si dichiara di eseguire. Non sceglie
    # l'immagine, la registra nella provenienza.
    container: RiferimentoImmagine | None = None
    # Il file di blocco delle versioni R. E' referenziato dal Dockerfile e
    # dagli script R, quindi l'ambiente sarebbe riproducibile comunque; entra
    # qui perche' cosi' finisce nella configurazione risolta e nel suo digest,
    # e due esecuzioni diventano confrontabili anche rispetto alle versioni R.
    lockfile: StringaNonVuota = d.RUN_LOCKFILE
    seed: int = d.RUN_SEED
    # Per difetto i processori utilizzabili dal processo: non incide sui
    # risultati, e un valore fisso fermerebbe G14 sulle macchine piu' piccole.
    threads: InteroPositivo = Field(default_factory=lambda: processori_disponibili())
    batch_size: InteroPositivo = d.RUN_BATCH_SIZE
    # Se conservare le letture filtrate da S2 a esecuzione conclusa. Con false
    # si rimuovono solo quando tutte le fasi sono concluse, e la rimozione e'
    # registrata: una ripresa non la scambia per un artefatto perso.
    keep_filtered_fastq: StrictBool = d.RUN_KEEP_FILTERED_FASTQ
    # La regola rigorosa sulla provenienza. Con true l'esecuzione parte solo da
    # un repository git leggibile, senza modifiche non committate al codice, e
    # in un ambiente R che corrisponde a run.lockfile; l'impronta del sorgente,
    # quella del file di blocco e l'immagine dichiarata entrano nell'impronta di
    # ogni fase. Vedi runner/provenienza.py per cio' che e' verificato e cio'
    # che e' solo dichiarato.
    strict_provenance: StrictBool = d.RUN_STRICT_PROVENANCE


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
    meta: Meta
    filter: Filter
    err: Err
    dada: Dada = Field(default_factory=Dada)
    chimera: Chimera = Field(default_factory=Chimera)
    asv: Asv = Field(default_factory=Asv)
    tax: Tax
    filt: Filt = Field(default_factory=Filt)
    phylo: Phylo = Field(default_factory=Phylo)
    ctrl: Ctrl
    katharoseq: Katharoseq
    decontam: Decontam
    prev: Prev = Field(default_factory=Prev)
    qc: Qc
    retry: Retry = Field(default_factory=Retry)
    out: Out
    run: Run = Field(default_factory=Run)

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


#: Come inizia il problema che elenca i parametri obbligatori non dichiarati:
#: G15 lo riconosce da qui per attribuirgli il proprio codice (E-G15-10).
INTESTAZIONE_MANCANTI: Final = "parametri obbligatori non dichiarati"


def obbligatori_mancanti(dati: dict[str, Any]) -> list[str]:
    """I parametri di :data:`defaults.OBBLIGATORI` assenti dal dizionario.

    Conta la presenza della chiave, non il valore: un parametro dichiarato
    nullo o vuoto e' dichiarato.
    """
    mancanti = []
    for chiave in d.OBBLIGATORI:
        gruppo, nome = chiave.split(".")
        contenuto = dati.get(gruppo)
        if not isinstance(contenuto, dict) or nome not in contenuto:
            mancanti.append(chiave)
    return mancanti


def _voce_mancanti(mancanti: list[str]) -> str:
    """La voce che elenca i parametri obbligatori non dichiarati."""
    return (
        f"{INTESTAZIONE_MANCANTI} ({len(mancanti)}): {', '.join(mancanti)}. Descrivono il "
        "dataset e non hanno un valore predefinito: vanno dichiarati tutti nel file di "
        "configurazione, nulli o vuoti quando non sono pertinenti (per esempio le "
        "etichette dei controlli positivi in un dataset che non ne ha)"
    )


def valida(dati: dict[str, Any], origine: str | None = None) -> Config:
    """Valida un dizionario già caricato.

    Solleva :class:`ErroreConfigurazione` con l'elenco completo dei problemi. I
    parametri obbligatori non dichiarati vi compaiono in una sola voce, tutti
    insieme: chi compila la configurazione per un dataset nuovo deve vedere
    l'elenco intero, non scoprirli uno per esecuzione.
    """
    mancanti = obbligatori_mancanti(dati)
    try:
        return Config.model_validate(dati)
    except ValidationError as errore:
        gia_detti = set(mancanti) | {m.split(".")[0] for m in mancanti}
        problemi = [
            _descrivi(e) for e in errore.errors()
            if not (e["type"] == "missing" and ".".join(str(p) for p in e["loc"]) in gia_detti)
        ]
        if mancanti:
            problemi.insert(0, _voce_mancanti(mancanti))
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
