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
    "thread_effettivi",
    "carica",
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


def _glob_di_soli_nomi(valore: str) -> str:
    """Rifiuta un modello con un separatore di percorso: i file di letture si
    cercano nella sola cartella ``io.fastq_dir``.
    """
    if "/" in valore or "\\" in valore:
        raise ValueError(
            "il modello si applica ai nomi dei file dentro io.fastq_dir e non puo' "
            "contenere separatori di percorso: metti la cartella in io.fastq_dir e "
            "raccogli i file in una sola cartella"
        )
    return valore


Regex = Annotated[str, Field(min_length=1), AfterValidator(_regex_valida)]
StringaNonVuota = Annotated[str, Field(min_length=1)]
# Interi in senso stretto: YAML legge "true" e "false" come booleani, che per
# Python sono anche gli interi 1 e 0, e un refuso diventerebbe un valore valido
# (un solo thread, un lotto di un campione, zero controlli richiesti).
InteroPositivo = Annotated[int, Field(gt=0, strict=True)]
InteroNonNegativo = Annotated[int, Field(ge=0, strict=True)]
# Finiti: un infinito non si scrive nel JSON con cui i parametri arrivano a R,
# e la prima fase che lo riceve si fermerebbe con un errore fuori catalogo.
RealePositivo = Annotated[float, Field(gt=0, allow_inf_nan=False)]
Frazione = Annotated[float, Field(ge=0.0, le=1.0, allow_inf_nan=False)]
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


def thread_effettivi(config: Any) -> int:
    """I thread con cui eseguire: ``run.threads`` se dichiarato, altrimenti i
    processori utilizzabili dal processo in questo momento.

    Si calcola all'uso (G14 e le fasi) e non alla validazione: un valore
    ricavato dalla macchina non deve entrare nella configurazione, altrimenti
    il suo digest cambierebbe da una macchina all'altra a parita' di file.
    """
    return config.run.threads or processori_disponibili()


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
    # Un modello sui soli nomi dei file di io.fastq_dir: le sottocartelle non
    # sono ammesse, perche' l'inventario registra il nome del file e le fasi lo
    # cercano in io.fastq_dir.
    fastq_glob: Annotated[StringaNonVuota, AfterValidator(_glob_di_soli_nomi)]
    # Obbligatorio: estrae dal nome di ogni file la chiave del campione, che e'
    # la corrispondenza intera oppure, se l'espressione ha un gruppo di
    # cattura, il primo gruppo. La stessa chiave si ricava dal valore di
    # meta.accession_column.
    accession_regex: Regex


class Meta(_Gruppo):
    """Lettura della tabella dei metadati e derivazione dei campi."""

    # I parametri senza predefinito sono obbligatori (defaults.OBBLIGATORI):
    # in tutto lo schema, un campo senza valore dipende dai dati o dallo
    # studio e va dichiarato. Quelli che ammettono il valore nullo si
    # dichiarano nulli quando il dataset non ha l'informazione.
    sample_id_column: StringaNonVuota
    accession_column: StringaNonVuota
    # Colonna del nome del campione nella tabella di studio, quando ha un nome
    # diverso da quello della tabella di assay; nullo: lo stesso nome.
    study_sample_id_column: StringaNonVuota | None = None
    derive_module: StrictBool
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
    truncLen_shortfall_warn: InteroNonNegativo
    trimLeft: InteroNonNegativo
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

    error_function: Literal["loess", "loess_monotono"] = d.ERR_ERROR_FUNCTION
    nbases: RealePositivo = d.ERR_NBASES
    max_consist: InteroPositivo = d.ERR_MAX_CONSIST
    randomize: StrictBool
    # Obbligatorio, anche nullo. Colonna di io.batch_table con la
    # corsa di sequenziamento: un modello d'errore per corsa. Con null, o senza
    # io.batch_table, si stima un solo modello su tutti i campioni.
    batch_column: StringaNonVuota | None


class Dada(_Gruppo):
    """Inferenza delle varianti di sequenza."""

    # dada2 accetta TRUE, FALSE oppure "pseudo": il parametro e' passato cosi'
    # com'e' ed eredita il vocabolario della libreria.
    pool: StrictBool | Literal["pseudo"]
    omega_a: RealePositivo = d.DADA_OMEGA_A


class Chimera(_Gruppo):
    """Rimozione delle sequenze chimeriche."""

    # Vocabolario di dada2::removeBimeraDenovo.
    method: Literal["consensus", "pooled", "per-sample"] = d.CHIMERA_METHOD
    # Obbligatorio: quanto un genitore deve essere piu' abbondante della
    # chimera. Il predefinito del pacchetto dipende dal metodo e dalla
    # versione, e il valore adatto dalla profondita' dei campioni.
    min_fold_parent_over_abundance: Annotated[float, Field(ge=1.0)]
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
    # dati che l'hanno prodotto.
    ref_name: StringaNonVuota
    ref_version: StringaNonVuota
    # Insieme chiuso: va esteso quando un altro classificatore viene realizzato.
    classifier: Literal["naive_bayes"] = d.TAX_CLASSIFIER
    min_boot: Annotated[int, Field(ge=0, le=100, strict=True)] = d.TAX_MIN_BOOT
    try_rc: StrictBool
    assign_species: StrictBool = d.TAX_ASSIGN_SPECIES

    @model_validator(mode="after")
    def _specie_non_realizzata(self) -> Tax:
        """Respinge ``assign_species`` vero: l'assegnazione della specie non e'
        realizzata, e accettare il parametro farebbe credere il contrario.
        """
        if self.assign_species:
            raise ValueError(
                "assign_species: true non e' realizzato: la pipeline assegna la "
                "tassonomia fino al genere. Imposta false"
            )
        return self
    # Elenco dei taxa con un difetto noto del riferimento, nella forma
    # rank,name,path,ranks_present,ranks_expected: S8 marca le assegnazioni che
    # vi ricadono. Facoltativo: non tutti i riferimenti ne hanno uno.
    ref_bad_taxa: Path | None = d.TAX_REF_BAD_TAXA


class Filt(_Gruppo):
    """Filtraggio tassonomico successivo all'assegnazione."""

    remove_na_phylum: StrictBool
    exclude_taxa: list[StringaNonVuota]


class Phylo(_Gruppo):
    """Costruzione dell'albero filogenetico."""

    enabled: StrictBool
    # Varianti dell'oggetto filtrato (S13) oltre le quali S9 si ferma prima del
    # calcolo (E-S9-01).
    max_seqs: InteroPositivo
    # Insiemi chiusi: vanno estesi quando un'altra forma viene realizzata.
    # L'allineatore multiplo (DECIPHER::AlignSeqs) e il modello evolutivo
    # dell'albero di massima verosimiglianza (phangorn): GTR con eterogeneita'
    # gamma fra i siti e una quota di siti invarianti.
    aligner: Literal["decipher"] = d.PHYLO_ALIGNER
    model: Literal["GTR+G+I"] = d.PHYLO_MODEL


class Ctrl(_Gruppo):
    """Riconoscimento di controlli negativi, positivi e campioni biologici, e
    validazione della corsa dai controlli positivi.

    La colonna, le etichette e la regola di riclassificazione descrivono il
    dataset, non la pipeline: sono tutte da dichiarare.
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
    min_positives: InteroPositivo
    min_positive_pass_frac: Frazione
    positive_gate: StrictBool

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
    # Colonna con le cellule di ciascun controllo positivo, il livello di
    # diluizione; nome originale della colonna, cercata nel file di
    # arricchimento e poi nella tabella di studio (S10 la porta nell'oggetto).
    cell_count_column: StringaNonVuota | None
    collapse_rank: Literal["Phylum", "Class", "Order", "Family", "Genus"]
    # Insiemi chiusi: vanno estesi quando un'altra forma viene realizzata.
    curve_model: Literal["allosteric_sigmoid"] = d.KATHAROSEQ_CURVE_MODEL
    # Strettamente fra 0 e 1: a fedelta' 1 la profondita' richiesta e' infinita.
    target_sensitivity: Annotated[float, Field(gt=0.0, lt=1.0)]
    min_r2: Frazione
    read_stage: Literal["nonchimeric"] = d.KATHAROSEQ_READ_STAGE


class Decontam(_Gruppo):
    """Rimozione dei contaminanti a partire dai controlli negativi."""

    # Vocabolario di decontam::isContaminant.
    # Insieme chiuso: la sola decontaminazione realizzata e' quella per
    # prevalenza, che non richiede la concentrazione del DNA per campione. Gli
    # altri metodi di decontam vanno aggiunti qui quando una fase li realizza:
    # accettarli prima fermerebbe la catena in S12, dopo tutto il calcolo.
    method: Literal["prevalence"] = d.DECONTAM_METHOD
    threshold: Frazione
    min_blanks: InteroPositivo
    # batch.combine di decontam::isContaminant, per la decontaminazione per piastra.
    batch_combine: Literal["minimum", "product", "fisher"] = d.DECONTAM_BATCH_COMBINE
    # Modalita' che decide i contaminanti rimossi: aggregate (tutti i biologici
    # contro tutti i negativi) o batch (per piastra, con batch_combine). L'altra
    # si calcola come diagnostica.
    mode: Literal["aggregate", "batch"]
    # Obbligatorio, anche nullo: la colonna di io.batch_table con la piastra.
    batch_column: StringaNonVuota | None


class Prev(_Gruppo):
    """Filtro di prevalenza sulle varianti."""

    min_fraction: Frazione
    min_count: InteroNonNegativo
    apply: StrictBool


class Qc(_Gruppo):
    """Soglie dei controlli di qualità che governano avvisi e arresti."""

    # katharoseq_if_available: le soglie di S11, sulle letture senza chimere;
    # none: nessuna soglia di profondita', resta min_reads_final.
    min_reads_mode: Literal["katharoseq_if_available", "none"]
    # Frazione massima delle letture dei biologici rimossa come contaminante (S12).
    max_frac_contaminant: Frazione
    # Per campione, sulle letture dell'oggetto finale (S13): sotto, il campione
    # esce dall'oggetto finale.
    min_reads_final: InteroNonNegativo
    # Per l'insieme dei campioni finali: letture finali su letture senza chimere
    # (S14, E-S14-01).
    min_frac_reads_retained: Frazione
    max_asv_count: InteroPositivo
    # Controlli sul risultato del filtro (S2), applicati ai campioni biologici
    # e ai controlli positivi, non ai negativi: vedi steps/s02_filter.py.
    max_zeroed_samples: InteroNonNegativo
    max_frac_lost_filter: Frazione
    max_frac_short_reads: Frazione
    warn_frac_chimeric: Frazione
    stop_frac_chimeric: Frazione
    # Frazione minima di varianti con il phylum assegnato (S8, E-S8-02): vedi
    # steps/s08_taxonomy.py per le classi a cui si applica.
    min_frac_phylum: Frazione
    # Letture ispezionate per file dai gate che leggono le sequenze.
    head_reads: InteroPositivo
    max_primer_hit_frac: Frazione
    min_motif_frac: Frazione
    # Obbligatori: dipendono dalla regione amplificata. Il motivo puo' essere
    # nullo, per una regione senza un motivo noto: G10 verifica allora la sola
    # assenza del primer.
    primer_sequence: SequenzaIupac
    conserved_motif: Regex | None

    @model_validator(mode="before")
    @classmethod
    def _soglia_sulle_letture_grezze_rimossa(cls, valori: Any) -> Any:
        """Respinge ``min_reads_raw`` e il modo ``fixed`` spiegando il cambiamento.

        Una configurazione scritta per la soglia fissa sulle letture grezze
        chiedeva un filtro che non esiste piu': respingerla come chiave
        sconosciuta non direbbe che cosa lo sostituisce.
        """
        if isinstance(valori, dict):
            if "min_reads_raw" in valori:
                raise ValueError(
                    "min_reads_raw non esiste piu': nessuna soglia di profondita' si "
                    "applica alle letture grezze. La soglia viene dai controlli positivi "
                    "(S11), sulle letture senza chimere, e dove una piastra non ha una "
                    "curva valida si usa la soglia aggregata o la mediana delle altre "
                    "piastre; senza curve valide resta min_reads_final. Togli la chiave"
                )
            if valori.get("min_reads_mode") == "fixed":
                raise ValueError(
                    "min_reads_mode: fixed non esiste piu', perche' applicava "
                    "min_reads_raw alle letture grezze. I valori ammessi sono "
                    "katharoseq_if_available (le soglie di S11, sulle letture senza "
                    "chimere) e none (nessuna soglia di profondita', resta "
                    "min_reads_final)"
                )
        return valori

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

    @model_validator(mode="after")
    def _orientamento_realizzato(self) -> Out:
        """Respinge ``taxa_are_rows`` falso: i filtri finali lavorano con le
        varianti sulle righe, e accettare l'altro orientamento fermerebbe la
        catena dopo tutto il calcolo.
        """
        if not self.taxa_are_rows:
            raise ValueError(
                "taxa_are_rows: false non e' realizzato: le fasi dopo l'oggetto "
                "integrato lavorano con le varianti sulle righe. Imposta true"
            )
        return self
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
    # Entro gli interi di R, a cui il seme viene passato.
    seed: Annotated[int, Field(strict=True, ge=-2147483647, le=2147483647)] = d.RUN_SEED
    # Nullo per difetto, cioe' automatico: si usano i processori utilizzabili
    # dal processo, contati all'uso (thread_effettivi). Il valore ricavato dalla
    # macchina non entra nella configurazione ne' nel suo digest, che resta lo
    # stesso su macchine diverse; un valore fisso per difetto fermerebbe G14
    # sulle macchine piu' piccole.
    threads: InteroPositivo | None = None
    # Tempo massimo, in secondi, di ogni processo R di una fase. Nullo per
    # difetto: nessun limite, perche' la durata di una fase dipende dal dataset
    # e dalla macchina. Allo scadere il ponte uccide il processo e il suo
    # gruppo (E-R-05). Non incide sui risultati: o la fase si conclude, o si
    # ferma.
    r_timeout_s: InteroPositivo | None = None
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
    dada: Dada
    chimera: Chimera
    asv: Asv = Field(default_factory=Asv)
    tax: Tax
    filt: Filt
    phylo: Phylo
    ctrl: Ctrl
    katharoseq: Katharoseq
    decontam: Decontam
    prev: Prev
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
        "dataset o le scelte di analisi che ne dipendono, e non hanno un valore "
        "predefinito: vanno dichiarati tutti nel file di configurazione, nulli o vuoti "
        "quando non sono pertinenti e lo schema lo ammette (per esempio le etichette dei "
        "controlli positivi in un dataset che non ne ha). Che cosa significa ciascuno e "
        "come si sceglie e' in config/config.example.yaml"
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
        # Si tacciono le sole chiavi gia' elencate fra gli obbligatori. Un gruppo
        # assente e' per pydantic un solo errore sul gruppo: lo si scioglie nei
        # suoi campi senza predefinito, cosi' quelli che non sono in OBBLIGATORI
        # (un percorso, un checksum) restano nell'elenco invece di sparire con
        # il gruppo.
        gia_detti = set(mancanti)
        problemi = []
        for e in errore.errors():
            chiave = ".".join(str(p) for p in e["loc"])
            if e["type"] != "missing":
                problemi.append(_descrivi(e))
            elif chiave in Config.model_fields:
                problemi += [
                    f"{chiave}.{nome}: {_MESSAGGI_SEMPLICI['missing']}"
                    for nome, campo in Config.model_fields[chiave].annotation.model_fields.items()
                    if campo.is_required() and f"{chiave}.{nome}" not in gia_detti
                ]
            elif chiave not in gia_detti:
                problemi.append(_descrivi(e))
        if mancanti:
            problemi.insert(0, _voce_mancanti(mancanti))
        raise ErroreConfigurazione(problemi, origine) from errore


class _LettoreSenzaRipetizioni(yaml.SafeLoader):
    """Il lettore YAML sicuro, che in piu' rifiuta una chiave ripetuta.

    PyYAML tiene in silenzio l'ultima occorrenza: un gruppo scritto due volte
    perderebbe i parametri del primo blocco, che tornerebbero ai predefiniti
    senza alcun messaggio. E' lo stesso rischio di un refuso in una chiave.
    """

    def construct_mapping(self, node: yaml.MappingNode, deep: bool = False) -> dict[Any, Any]:
        viste: set[Any] = set()
        for chiave_nodo, _ in node.value:
            chiave = self.construct_object(chiave_nodo, deep=deep)
            if chiave in viste:
                raise yaml.constructor.ConstructorError(
                    None, None,
                    f"la chiave {chiave!r} compare due volte nella stessa mappa: i "
                    "parametri di un gruppo vanno scritti in un solo blocco",
                    chiave_nodo.start_mark,
                )
            viste.add(chiave)
        return super().construct_mapping(node, deep=deep)


def carica(percorso: str | Path) -> Config:
    """Carica e valida un file di configurazione YAML."""
    percorso = Path(percorso)
    try:
        testo = percorso.read_text(encoding="utf-8")
    except UnicodeDecodeError as errore:
        raise ErroreConfigurazione(
            [f"il file non e' testo UTF-8: {errore}. Salvalo in UTF-8 (per esempio con iconv)"],
            str(percorso),
        ) from errore

    try:
        dati = yaml.load(testo, Loader=_LettoreSenzaRipetizioni)  # noqa: S506 - deriva da SafeLoader
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
