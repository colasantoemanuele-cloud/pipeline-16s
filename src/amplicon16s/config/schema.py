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
from collections.abc import Iterable
from pathlib import Path
from typing import Annotated, Any, Final, Literal

import yaml
from pydantic import (
    AfterValidator,
    BaseModel,
    BeforeValidator,
    ConfigDict,
    Field,
    StrictBool,
    ValidationError,
    field_validator,
    model_serializer,
    model_validator,
)

from amplicon16s.config import defaults as d

__all__ = [
    "Config",
    "ErroreConfigurazione",
    "INTESTAZIONE_CLASSIFICATORE",
    "INTESTAZIONE_COPPIE",
    "INTESTAZIONE_TAGLIO",
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


#: I costrutti di un'espressione regolare le cui lettere non sono basi: le
#: sequenze con la barra rovesciata (``\d``, ``\x41``), i gruppi con nome, con
#: opzioni o condizionali (``(?P<nome>``, ``(?P=nome)``, ``(?i)``, ``(?(1)``).
_LETTERE_NON_BASI: Final = re.compile(
    r"\\.|\(\?P<\w+>|\(\?P=\w+\)|\(\?[aiLmsux-]+[:)]|\(\?\(\w+\)")
#: Quelli in cui portare in maiuscolo le altre lettere non e' sicuro: i gruppi
#: di cui sopra e le sequenze con argomento (``\x41``, ``\N{...}``, ``\g<nome>``).
_COSTRUTTI_CON_LETTERE: Final = re.compile(r"\(\?[A-Za-z(]|\\[xuUNgk]")


def _motivo_in_maiuscolo(valore: str) -> str:
    r"""Porta in maiuscolo le basi di un motivo scritto in minuscolo.

    Le letture sono in maiuscolo: un motivo in minuscolo non corrisponderebbe
    mai. Una lettera preceduta da una barra rovesciata non e' una base (``\d``,
    ``\s``, ``\w``, ``\b`` in maiuscolo sono i loro opposti) e resta com'e'. Se
    l'espressione ha basi in minuscolo e un costrutto in cui la regola non e'
    sicura, si respinge: le basi vanno scritte in maiuscolo da chi la conosce.
    Un'espressione senza basi in minuscolo resta identica.
    """
    basi = _LETTERE_NON_BASI.sub("", valore)
    if basi == basi.upper():
        return valore
    if _COSTRUTTI_CON_LETTERE.search(valore):
        raise ValueError(
            "l'espressione ha lettere minuscole e un costrutto (un gruppo con nome o con "
            "opzioni, una sequenza con argomento) in cui non si possono portare in "
            "maiuscolo senza cambiarne il significato: scrivi le basi in maiuscolo, perche' "
            "le letture lo sono"
        )
    return _regex_valida(re.sub(
        r"\\.|[a-z]", lambda t: t.group() if t.group()[0] == "\\" else t.group().upper(), valore))


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
def _non_booleano(valore: Any) -> Any:
    """Rifiuta ``true`` e ``false`` dove serve un numero.

    YAML li legge come booleani, che per Python sono 1 e 0: una soglia scritta
    ``true`` varrebbe 1 senza alcun messaggio.
    """
    if isinstance(valore, bool):
        raise ValueError("deve essere un numero, non true o false")
    return valore


# Finiti: un infinito non si scrive nel JSON con cui i parametri arrivano a R,
# e la prima fase che lo riceve si fermerebbe con un errore fuori catalogo.
# Il rifiuto dei booleani viene dopo i vincoli nell'annotazione: cosi' lo
# schema JSON del parametro resta quello di un numero con il suo intervallo.
NonBooleano = BeforeValidator(_non_booleano)
RealePositivo = Annotated[float, Field(gt=0, allow_inf_nan=False), NonBooleano]
Frazione = Annotated[float, Field(ge=0.0, le=1.0, allow_inf_nan=False), NonBooleano]
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

def _md5_non_numerico(valore: Any) -> Any:
    """Rifiuta un MD5 che YAML ha letto come numero.

    Un MD5 di sole cifre scritto senza virgolette e' per YAML un numero intero:
    convertirlo in testo ne perderebbe gli zeri iniziali (e, se comincia per
    zero con sole cifre fino a 7, ne cambierebbe la base), dando un checksum
    diverso da quello scritto. Va scritto fra virgolette.
    """
    if isinstance(valore, (int, float)) and not isinstance(valore, bool):
        raise ValueError(
            "e' stato letto come numero: un MD5 di sole cifre va scritto fra virgolette "
            "(per esempio \"01234567890123456789012345678901\"), altrimenti YAML lo legge "
            "come numero e ne perde gli zeri iniziali"
        )
    return valore


# Il rifiuto del numero viene dopo il vincolo nell'annotazione: cosi' lo schema
# JSON del parametro resta quello di una stringa con il suo formato.
Md5 = Annotated[
    str, Field(pattern=r"^[0-9a-fA-F]{32}$"), BeforeValidator(_md5_non_numerico)]

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

    fastq_dir: Path = Field(description=(
        "Cartella con i file delle letture, uno per campione, senza sottocartelle."))
    assay_table: Path = Field(description=(
        "Tabella di assay, o tabella dei campioni: una riga per campione, con il nome "
        "del campione e il riferimento al suo file di letture. Testo UTF-8 separato da "
        "tabulazioni."))
    # Tabella campioni di studio: porta la colonna della classe del campione,
    # che nella tabella di assay non c'e'. E' condivisa fra piu' assay dello
    # stesso studio, quindi il join verso di essa va tenuto ristretto.
    # Facoltativa: senza, classe e variabili dei campioni si leggono dalla
    # tabella di assay, che allora deve portarne le colonne.
    study_table: Path | None = Field(None, description=(
        "Tabella di studio con il tipo di ogni campione e le sue variabili. Nullo: "
        "tipo e variabili si leggono dalla tabella di assay."))
    out_root: Path = Field(description=(
        "Cartella di uscita dell'esecuzione: una sottocartella per fase."))
    # Associa a ciascun campione la piastra di estrazione e la corsa di
    # sequenziamento. E' facoltativo: quando manca, entrambe restano nulle e la
    # pipeline procede in modalita' a corsa singola e senza lotto.
    batch_table: Path | None = Field(None, description=(
        "File del lotto: per ogni campione la piastra di estrazione e la corsa di "
        "sequenziamento. Nullo: corsa singola e nessun lotto."))
    # Un modello sui soli nomi dei file di io.fastq_dir: le sottocartelle non
    # sono ammesse, perche' l'inventario registra il nome del file e le fasi lo
    # cercano in io.fastq_dir.
    fastq_glob: Annotated[StringaNonVuota, AfterValidator(_glob_di_soli_nomi)] = Field(
        description=(
            "Modello dei nomi dei file delle letture dentro io.fastq_dir, per esempio "
            "l'estensione comune preceduta da un asterisco."))
    # Obbligatorio: estrae dal nome di ogni file la chiave del campione, che e'
    # la corrispondenza intera oppure, se l'espressione ha un gruppo di
    # cattura, il primo gruppo. La stessa chiave si ricava dal valore di
    # meta.accession_column.
    accession_regex: Regex = Field(description=(
        "Espressione regolare che estrae dal nome di ogni file la chiave del campione: "
        "la corrispondenza intera oppure il primo gruppo di cattura. La stessa chiave "
        "deve ricavarsi dal valore di meta.accession_column, ed essere unica per "
        "campione."))
    # La modalita' di lettura. L'assenza e' la modalita' single-end: dichiararla
    # single-end vale come ometterla, e il valore si normalizza a nullo perche'
    # le due forme siano la stessa configurazione, con lo stesso digest. In
    # modalita' single-end la chiave non compare nella configurazione risolta
    # (Config._senza_la_modalita_appaiata).
    layout: Literal["single-end", "paired-end"] | None = Field(None, description=(
        "Modalita' di lettura: single-end (un file per campione, le letture forward) "
        "oppure paired-end (due file per campione, fusi dopo l'inferenza). Omesso: "
        "single-end."))

    @field_validator("layout", mode="after")
    @classmethod
    def _single_end_e_l_assenza(cls, valore: str | None) -> str | None:
        """Porta ``single-end`` al valore nullo: e' la modalita' di chi non dichiara nulla."""
        return None if valore == d.LAYOUT_SINGOLO else valore


class Meta(_Gruppo):
    """Lettura della tabella dei metadati e derivazione dei campi."""

    # I parametri senza predefinito sono obbligatori (defaults.OBBLIGATORI):
    # in tutto lo schema, un campo senza valore dipende dai dati o dallo
    # studio e va dichiarato. Quelli che ammettono il valore nullo si
    # dichiarano nulli quando il dataset non ha l'informazione. I condizionali
    # (defaults.CONDIZIONALI) hanno invece un valore di assenza, elenco vuoto o
    # nullo: si dichiarano solo quando la loro condizione e' vera, e lo
    # verifica la radice (Config._condizionali_richiesti).
    sample_id_column: StringaNonVuota = Field(description=(
        "Colonna della tabella di assay con il nome del campione."))
    accession_column: StringaNonVuota = Field(description=(
        "Colonna della tabella di assay da cui si ricava la chiave del campione con "
        "io.accession_regex: di solito quella che riporta il file delle letture."))
    # Colonna del nome del campione nella tabella di studio, quando ha un nome
    # diverso da quello della tabella di assay; nullo: lo stesso nome.
    study_sample_id_column: StringaNonVuota | None = Field(None, description=(
        "Colonna del nome del campione nella tabella di studio, se ha un nome diverso "
        "da meta.sample_id_column. Nullo: lo stesso nome."))
    derive_module: StrictBool = Field(description=(
        "Se derivare il modulo (il gruppo di posizioni a cui il campione appartiene) "
        "dalla posizione dichiarata, con meta.module_regex, quando il file del lotto "
        "non lo dichiara."))
    module_regex: Regex | None = Field(None, description=(
        "Espressione regolare, con un gruppo di cattura, che estrae il modulo dalla "
        "posizione del campione. Nullo: nessun modulo derivato. Serve solo con "
        "meta.derive_module vero."))
    module_column: StringaNonVuota | None = Field(description=(
        "Colonna della tabella di studio con la posizione di prelievo. Nullo: il "
        "dataset non dichiara posizioni."))
    # Posizioni che non sono superfici. Vale in entrambe le vie di
    # attribuzione del modulo: anche un modulo dichiarato dal file di
    # arricchimento non viene attribuito a un campione che non sta su una
    # superficie, altrimenti lo stesso dataset darebbe raggruppamenti diversi
    # a seconda che quel file ci sia o no.
    non_surface_positions: list[StringaNonVuota] = Field(description=(
        "Valori di meta.module_column che non sono superfici: i campioni che li "
        "riportano non appartengono a un modulo. Elenco vuoto se non ce ne sono."))
    # Colonna con cui il file facoltativo identifica il campione: deve
    # contenere l'accession, non il nome.
    batch_key_column: StringaNonVuota | None = Field(description=(
        "Colonna con cui il file del lotto identifica il campione: deve contenere la "
        "chiave estratta da io.accession_regex. Nullo senza file del lotto."))
    # Colonna del file facoltativo che dichiara il modulo. Ha la precedenza
    # sulla derivazione da module_regex, che resta il ripiego quando il file
    # non c'e'.
    batch_module_column: StringaNonVuota | None = Field(description=(
        "Colonna del file del lotto che dichiara il modulo; ha la precedenza sulla "
        "derivazione da meta.module_regex. Nullo se il file non la porta."))

    @property
    def colonna_id_studio(self) -> str:
        """La colonna del nome del campione nella tabella di studio."""
        return self.study_sample_id_column or self.sample_id_column


class Filter(_Gruppo):
    """Filtraggio e troncamento delle letture grezze."""

    truncLen: InteroPositivo = Field(description=(
        "Lunghezza a cui si troncano le letture; quelle piu' corte sono scartate. Si "
        "sceglie dal profilo delle lunghezze e della qualita' delle proprie letture."))
    truncLen_shortfall_warn: InteroNonNegativo = Field(description=(
        "Di quante basi il troncamento puo' stare sotto la lettura piu' corta "
        "osservata prima che il profilo lo segnali (E-S1-01)."))
    trimLeft: InteroNonNegativo = Field(description=(
        "Basi tolte in testa a ogni lettura: la lunghezza del primer se le letture lo "
        "contengono ancora, zero se e' gia' stato rimosso."))
    maxEE: RealePositivo = Field(d.FILTER_MAXEE, description=(
        "Errori attesi massimi per lettura (dada2::filterAndTrim, maxEE)."))
    truncQ: InteroNonNegativo = Field(d.FILTER_TRUNCQ, description=(
        "Qualita' alla quale una lettura viene troncata (dada2::filterAndTrim, truncQ)."))
    maxN: InteroNonNegativo = Field(d.FILTER_MAXN, description=(
        "Basi indeterminate ammesse per lettura (dada2::filterAndTrim, maxN)."))
    rm_phix: StrictBool = Field(d.FILTER_RM_PHIX, description=(
        "Se scartare le letture del genoma di controllo PhiX (dada2::filterAndTrim, "
        "rm.phix)."))
    # Come si toglie il primer. A posizione fissa e' l'assenza: dichiararlo
    # fixed vale come ometterlo, e il valore si normalizza a nullo perche' le
    # due forme siano la stessa configurazione, con lo stesso digest. Con il
    # taglio a posizione fissa la chiave e i parametri del taglio cercato non
    # compaiono nella configurazione risolta (Config._senza_le_modalita_assenti).
    primer_removal: Literal["fixed", "search"] | None = Field(None, description=(
        "Come si toglie il primer: fixed (le prime filter.trimLeft basi di ogni "
        "lettura) oppure search (il primer si cerca nelle prime basi della lettura, "
        "che si taglia dove finisce: per i protocolli che lo fanno precedere da un "
        "numero variabile di basi). Omesso: fixed."))
    primer_max_mismatch: InteroNonNegativo | None = Field(None, description=(
        "Taglio cercato: discordanze ammesse fra il primer e la lettura (sole "
        "sostituzioni; una N della lettura e' una discordanza)."))
    primer_max_offset: InteroNonNegativo | None = Field(None, description=(
        "Taglio cercato: quante basi possono precedere il primer. Zero: il primer "
        "deve cominciare alla prima base."))
    primer_untrimmed: Literal["discard", "keep"] | None = Field(None, description=(
        "Taglio cercato: che cosa succede a una lettura in cui il primer non si "
        "trova. discard: e' scartata (in modalita' appaiata con la sua compagna); "
        "keep: resta intera."))

    @field_validator("primer_removal", mode="after")
    @classmethod
    def _fisso_e_l_assenza(cls, valore: str | None) -> str | None:
        """Porta ``fixed`` al valore nullo: e' il modo di chi non dichiara nulla."""
        return None if valore == d.PRIMER_FISSO else valore

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

    error_function: Literal["loess", "loess_monotono"] = Field(
        d.ERR_ERROR_FUNCTION, description=(
            "Funzione che stima i tassi d'errore: loess (quella di dada2) oppure "
            "loess_monotono, per le qualita' raggruppate in pochi valori."))
    nbases: RealePositivo = Field(d.ERR_NBASES, description=(
        "Basi minime su cui si stima il modello d'errore (dada2::learnErrors, nbases)."))
    max_consist: InteroPositivo = Field(d.ERR_MAX_CONSIST, description=(
        "Iterazioni massime della stima del modello d'errore (dada2::learnErrors, "
        "MAX_CONSIST)."))
    randomize: StrictBool = Field(description=(
        "Se scegliere a caso, invece che nell'ordine dei file, i campioni da cui si "
        "stima il modello d'errore (dada2::learnErrors, randomize)."))
    # Obbligatorio, anche nullo. Colonna di io.batch_table con la
    # corsa di sequenziamento: un modello d'errore per corsa. Con null, o senza
    # io.batch_table, si stima un solo modello su tutti i campioni.
    batch_column: StringaNonVuota | None = Field(description=(
        "Colonna del file del lotto con la corsa di sequenziamento: un modello "
        "d'errore per corsa. Nullo: un solo modello."))


class Dada(_Gruppo):
    """Inferenza delle varianti di sequenza."""

    # dada2 accetta TRUE, FALSE oppure "pseudo": il parametro e' passato cosi'
    # com'e' ed eredita il vocabolario della libreria.
    pool: StrictBool | Literal["pseudo"] = Field(description=(
        "Come i campioni condividono l'informazione nell'inferenza: false (ogni "
        "campione da solo), true (tutti insieme) o pseudo (due passate) "
        "(dada2::dada, pool)."))
    omega_a: RealePositivo = Field(d.DADA_OMEGA_A, description=(
        "Soglia di significativita' con cui una sequenza e' riconosciuta come "
        "variante (dada2::setDadaOpt, OMEGA_A)."))


class Chimera(_Gruppo):
    """Rimozione delle sequenze chimeriche."""

    # Vocabolario di dada2::removeBimeraDenovo.
    method: Literal["consensus", "pooled", "per-sample"] = Field(
        d.CHIMERA_METHOD, description=(
            "Metodo di riconoscimento delle chimere (dada2::removeBimeraDenovo, method)."))
    # Obbligatorio: quanto un genitore deve essere piu' abbondante della
    # chimera. Il predefinito del pacchetto dipende dal metodo e dalla
    # versione, e il valore adatto dalla profondita' dei campioni.
    min_fold_parent_over_abundance: Annotated[float, Field(ge=1.0), NonBooleano] = Field(description=(
        "Quanto un genitore deve essere piu' abbondante di una variante perche' questa "
        "possa essere giudicata chimerica (dada2, minFoldParentOverAbundance)."))
    min_parent_abundance: InteroPositivo = Field(
        d.CHIMERA_MIN_PARENT_ABUNDANCE, description=(
            "Letture minime di un genitore (dada2::isBimeraDenovoTable, "
            "minParentAbundance)."))
    min_sample_fraction: Frazione = Field(d.CHIMERA_MIN_SAMPLE_FRACTION, description=(
        "Frazione minima dei campioni in cui una variante deve risultare chimerica "
        "(dada2::isBimeraDenovoTable, minSampleFraction)."))
    allow_one_off: StrictBool = Field(d.CHIMERA_ALLOW_ONE_OFF, description=(
        "Se ammettere una differenza dal modello a due genitori "
        "(dada2::isBimeraDenovoTable, allowOneOff)."))


class Asv(_Gruppo):
    """Selezione delle sequenze varianti."""

    len_tol: InteroNonNegativo = Field(d.ASV_LEN_TOL, description=(
        "Tolleranza sulla lunghezza delle varianti rispetto a quella attesa "
        "(filter.truncLen meno filter.trimLeft)."))


class Tax(_Gruppo):
    """Assegnazione tassonomica e riferimento impiegato.

    Nome, versione e checksum del riferimento sono obbligatori: senza di essi
    un risultato non è riconducibile ai dati che l'hanno prodotto.
    """

    # Il riferimento FASTA, il suo MD5, tax.try_rc, tax.min_boot e
    # tax.ref_bad_taxa sono del bayesiano ingenuo; i parametri idtaxa_* di
    # IdTaxa. Quelli del classificatore non scelto sono respinti e non
    # compaiono nella configurazione risolta (defaults.SOLO_DEL_BAYESIANO,
    # defaults.PARAMETRI_DI_IDTAXA): senza valore sono nulli, e che siano
    # dichiarati con il loro classificatore lo verifica la radice.
    ref_fasta: Path | None = Field(None, description=(
        "Riferimento tassonomico nel formato dei training set di dada2."))
    ref_md5: Md5 | None = Field(None, description=(
        "MD5 del file del riferimento, verificato a ogni avvio."))
    # Obbligatori e senza valore predefinito: ereditare in silenzio un
    # riferimento tassonomico renderebbe il risultato non riconducibile ai
    # dati che l'hanno prodotto.
    ref_name: StringaNonVuota = Field(description=(
        "Nome del riferimento, come va riportato con i risultati."))
    ref_version: StringaNonVuota = Field(description=(
        "Versione del riferimento, come va riportata con i risultati."))
    classifier: Literal["naive_bayes", "idtaxa"] = Field(d.TAX_CLASSIFIER, description=(
        "Classificatore: naive_bayes (il bayesiano ingenuo di dada2::assignTaxonomy, "
        "su un riferimento FASTA) oppure idtaxa (DECIPHER::IdTaxa, su un classificatore "
        "addestrato). Omesso: naive_bayes."))
    min_boot: Annotated[int, Field(ge=0, le=100, strict=True)] = Field(
        d.TAX_MIN_BOOT, description=(
            "Bootstrap minimo per assegnare un rango (dada2::assignTaxonomy, minBoot)."))
    try_rc: StrictBool | None = Field(None, description=(
        "Se provare anche il complemento inverso di ogni variante "
        "(dada2::assignTaxonomy, tryRC)."))
    idtaxa_trained: Path | None = Field(None, description=(
        "Classificatore addestrato nel formato di DECIPHER: un file .RData con "
        "l'oggetto trainingSet prodotto da LearnTaxa (classificatore idtaxa)."))
    idtaxa_trained_md5: Md5 | None = Field(None, description=(
        "MD5 del file del classificatore addestrato, verificato a ogni avvio "
        "(classificatore idtaxa)."))
    idtaxa_threshold: Annotated[
        float, Field(ge=0.0, le=100.0, allow_inf_nan=False), NonBooleano] | None = Field(
        None, description=(
            "Confidenza minima, da 0 a 100, per assegnare un rango (DECIPHER::IdTaxa, "
            "threshold; classificatore idtaxa)."))
    idtaxa_strand: Literal["both", "top", "bottom"] | None = Field(None, description=(
        "Filamento su cui si classifica ogni variante: both (anche il complemento "
        "inverso), top, bottom (DECIPHER::IdTaxa, strand; classificatore idtaxa)."))
    assign_species: StrictBool = Field(d.TAX_ASSIGN_SPECIES, description=(
        "Assegnazione della specie: non realizzata, il solo valore ammesso e' false."))

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
    ref_bad_taxa: Path | None = Field(d.TAX_REF_BAD_TAXA, description=(
        "Elenco dei taxa con un difetto noto del riferimento (colonne name e rank): "
        "le assegnazioni che vi ricadono sono marcate. Nullo se non ce n'e' uno."))


class Filt(_Gruppo):
    """Filtraggio tassonomico successivo all'assegnazione."""

    remove_na_phylum: StrictBool = Field(description=(
        "Se togliere le varianti senza phylum assegnato."))
    exclude_taxa: list[StringaNonVuota] = Field(description=(
        "Nomi di taxa, a qualunque rango, le cui varianti si tolgono, scritti come li "
        "scrive il riferimento. Elenco vuoto: nessuno."))


class Phylo(_Gruppo):
    """Costruzione dell'albero filogenetico."""

    enabled: StrictBool = Field(description=(
        "Se costruire l'albero filogenetico delle varianti finali e aggiungerlo "
        "all'oggetto finale."))
    # Varianti dell'oggetto filtrato (S13) oltre le quali S9 si ferma prima del
    # calcolo (E-S9-01).
    max_seqs: InteroPositivo | None = Field(None, description=(
        "Varianti finali oltre le quali la filogenesi si ferma prima del calcolo "
        "(E-S9-01). Serve solo con phylo.enabled vero."))
    # Insiemi chiusi: vanno estesi quando un'altra forma viene realizzata.
    # L'allineatore multiplo (DECIPHER::AlignSeqs) e il modello evolutivo
    # dell'albero di massima verosimiglianza (phangorn): GTR con eterogeneita'
    # gamma fra i siti e una quota di siti invarianti.
    aligner: Literal["decipher"] = Field(d.PHYLO_ALIGNER, description=(
        "Allineatore multiplo: DECIPHER::AlignSeqs."))
    model: Literal["GTR+G+I"] = Field(d.PHYLO_MODEL, description=(
        "Modello evolutivo dell'albero di massima verosimiglianza (phangorn)."))


class Ctrl(_Gruppo):
    """Riconoscimento di controlli negativi, positivi e campioni biologici, e
    validazione della corsa dai controlli positivi.

    La colonna, le etichette e la regola di riclassificazione descrivono il
    dataset, non la pipeline: sono tutte da dichiarare.
    """

    # Obbligatori. Le etichette dei controlli possono essere elenchi vuoti (un
    # dataset senza controlli positivi o negativi lo dichiara cosi'); quelle dei
    # campioni biologici no, salvo che la colonna sia nulla: allora il dataset
    # non ha una colonna che distingua i campioni, sono tutti biologici e i tre
    # elenchi sono vuoti.
    column: StringaNonVuota | None = Field(description=(
        "Colonna (della tabella di studio, o di assay se quella manca) che dice di "
        "che tipo e' il campione. Nullo: il dataset non ha controlli e ogni campione "
        "e' biologico; le tre etichette si dichiarano allora vuote."))
    blank_values: list[StringaNonVuota] = Field([], description=(
        "Valori di ctrl.column per i controlli negativi. Elenco vuoto: nessuno."))
    positive_values: list[StringaNonVuota] = Field([], description=(
        "Valori di ctrl.column per i controlli positivi. Elenco vuoto: nessuno."))
    biological_values: list[StringaNonVuota] = Field([], description=(
        "Valori di ctrl.column per i campioni biologici: almeno uno, salvo con "
        "ctrl.column nullo."))
    # Riclassificazione in controllo negativo, indipendente dal materiale
    # dichiarato: i campioni il cui valore nella colonna indicata (della tabella
    # campioni di studio) e' fra quelli elencati. Il materiale resta quello
    # originale nell'inventario, cosi' la riclassificazione resta tracciabile.
    # Con l'elenco vuoto non si riclassifica nulla.
    blank_override_column: StringaNonVuota | None = Field(description=(
        "Colonna con cui riconoscere campioni da trattare come controlli negativi "
        "qualunque sia il tipo dichiarato. Nullo: nessuna riclassificazione."))
    blank_override_values: list[StringaNonVuota] = Field([], description=(
        "Valori di ctrl.blank_override_column che fanno di un campione un controllo "
        "negativo. Elenco vuoto: nessuno."))
    # Validazione della corsa dai controlli positivi (S11).
    # Condizionali: servono solo con controlli positivi dichiarati, e senza
    # valore sono nulli (defaults.CONDIZIONALI).
    min_positives: InteroPositivo | None = Field(None, description=(
        "Controlli positivi minimi per adattare una curva di calibrazione e per "
        "giudicare i controlli di un livello di diluizione. Serve solo con "
        "ctrl.positive_values non vuoto."))
    min_positive_pass_frac: Frazione | None = Field(None, description=(
        "Frazione minima di controlli positivi conformi al loro livello. Serve solo "
        "con ctrl.positive_values non vuoto."))
    positive_gate: StrictBool | None = Field(None, description=(
        "Che cosa fare sotto ctrl.min_positive_pass_frac: true ferma l'esecuzione "
        "(E-S11-03), false prosegue con un avviso (E-S11-04). Serve solo con "
        "ctrl.positive_values non vuoto."))

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
        if self.column is None:
            # Senza colonna non c'e' nulla con cui confrontare un'etichetta:
            # dichiararne una farebbe credere che qualche campione sia un controllo.
            dichiarate = [nome for nome in nomi if categorie[nome]]
            if dichiarate:
                raise ValueError(
                    "column e' nullo ma {} {} valori: senza la colonna del tipo di "
                    "campione ogni campione e' biologico e le tre etichette vanno "
                    "dichiarate vuote".format(
                        ", ".join(dichiarate), "ha" if len(dichiarate) == 1 else "hanno"
                    )
                )
        elif not self.biological_values and "biological_values" in self.model_fields_set:
            # Omesso del tutto e' un condizionale mancante: lo dice la radice.
            raise ValueError(
                "ctrl.biological_values deve contenere almeno 1 elemento: con ctrl.column "
                "dichiarata servono le etichette dei campioni biologici"
            )
        return self


class Katharoseq(_Gruppo):
    """Calibrazione KatharoSeq sui controlli positivi."""

    # Obbligatori; nulli se il dataset non ha controlli positivi (G15 verifica
    # che non lo siano quando ctrl.positive_values non e' vuoto).
    target_taxon: StringaNonVuota | None = Field(None, description=(
        "Taxon atteso nei controlli positivi, come lo scrive il riferimento al rango "
        "di katharoseq.collapse_rank. Nullo senza controlli positivi."))
    # Colonna con le cellule di ciascun controllo positivo, il livello di
    # diluizione; nome originale della colonna, cercata nel file di
    # arricchimento e poi nella tabella di studio (S10 la porta nell'oggetto).
    cell_count_column: StringaNonVuota | None = Field(description=(
        "Colonna con le cellule di ciascun controllo positivo, cioe' il livello di "
        "diluizione: nel file del lotto o nella tabella di studio. Nullo senza "
        "controlli positivi."))
    collapse_rank: Literal["Phylum", "Class", "Order", "Family", "Genus"] | None = Field(
        None, description=(
            "Rango a cui si sommano le varianti del taxon atteso: il piu' fine a cui il "
            "riferimento lo distingue con le letture disponibili. Serve solo con "
            "ctrl.positive_values non vuoto."))
    # Insiemi chiusi: vanno estesi quando un'altra forma viene realizzata.
    curve_model: Literal["allosteric_sigmoid"] = Field(
        d.KATHAROSEQ_CURVE_MODEL, description=(
            "Forma della curva di calibrazione: la sola realizzata."))
    # Strettamente fra 0 e 1: a fedelta' 1 la profondita' richiesta e' infinita.
    target_sensitivity: Annotated[float, Field(gt=0.0, lt=1.0), NonBooleano] | None = Field(
        None, description=(
            "Fedelta' (frazione delle letture del controllo sul taxon atteso) a cui si "
            "legge sulla curva la profondita' minima dei campioni. Serve solo con "
            "ctrl.positive_values non vuoto."))
    min_r2: Frazione | None = Field(None, description=(
        "Bonta' minima dell'adattamento (R^2 sulla fedelta') perche' una curva sia "
        "valida. Serve solo con ctrl.positive_values non vuoto."))
    read_stage: Literal["nonchimeric"] = Field(d.KATHAROSEQ_READ_STAGE, description=(
        "Stadio delle letture su cui si misurano profondita' e fedelta': senza "
        "chimere."))


class Decontam(_Gruppo):
    """Rimozione dei contaminanti a partire dai controlli negativi."""

    # Vocabolario di decontam::isContaminant.
    # Insieme chiuso: la sola decontaminazione realizzata e' quella per
    # prevalenza, che non richiede la concentrazione del DNA per campione. Gli
    # altri metodi di decontam vanno aggiunti qui quando una fase li realizza:
    # accettarli prima fermerebbe la catena in S12, dopo tutto il calcolo.
    method: Literal["prevalence"] = Field(d.DECONTAM_METHOD, description=(
        "Metodo di decontam::isContaminant: per prevalenza nei controlli negativi, il "
        "solo realizzato."))
    # Condizionali: soglia, negativi minimi e modalita' servono solo con
    # controlli negativi dichiarati (ctrl.blank_values o
    # ctrl.blank_override_values non vuoto), e senza valore sono nulli
    # (defaults.CONDIZIONALI).
    threshold: Frazione | None = Field(None, description=(
        "Soglia della probabilita' sotto la quale una variante e' un contaminante "
        "(decontam::isContaminant, threshold). Serve solo con controlli negativi "
        "dichiarati (ctrl.blank_values o ctrl.blank_override_values non vuoto)."))
    min_blanks: InteroPositivo | None = Field(None, description=(
        "Controlli negativi con letture minimi perche' la decontaminazione si esegua, "
        "e perche' una piastra abbia un confronto proprio. Serve solo con controlli "
        "negativi dichiarati (ctrl.blank_values o ctrl.blank_override_values non "
        "vuoto)."))
    # batch.combine di decontam::isContaminant, per la decontaminazione per piastra.
    batch_combine: Literal["minimum", "product", "fisher"] = Field(
        d.DECONTAM_BATCH_COMBINE, description=(
            "Come si combinano le probabilita' delle piastre (decontam::isContaminant, "
            "batch.combine)."))
    # Modalita' che decide i contaminanti rimossi: aggregate (tutti i biologici
    # contro tutti i negativi) o batch (per piastra, con batch_combine). L'altra
    # si calcola come diagnostica.
    mode: Literal["aggregate", "batch"] | None = Field(None, description=(
        "Modalita' che decide i contaminanti rimossi: aggregate (tutti i biologici "
        "contro tutti i negativi) o batch (per piastra). L'altra si calcola come "
        "diagnostica. Serve solo con controlli negativi dichiarati (ctrl.blank_values "
        "o ctrl.blank_override_values non vuoto)."))
    # Obbligatorio, anche nullo: la colonna di io.batch_table con la piastra.
    batch_column: StringaNonVuota | None = Field(description=(
        "Colonna del file del lotto con la piastra di estrazione. Nullo senza file "
        "del lotto."))


class Prev(_Gruppo):
    """Filtro di prevalenza sulle varianti."""

    # Condizionali: servono solo con prev.apply vero, e senza valore sono
    # nulli (defaults.CONDIZIONALI).
    min_fraction: Frazione | None = Field(None, description=(
        "Frazione minima dei campioni biologici tenuti in cui una variante deve "
        "comparire. Serve solo con prev.apply vero."))
    min_count: InteroNonNegativo | None = Field(None, description=(
        "Letture minime perche' una variante conti come presente in un campione. "
        "Serve solo con prev.apply vero."))
    apply: StrictBool = Field(description="Se applicare il filtro di prevalenza.")


class Qc(_Gruppo):
    """Soglie dei controlli di qualità che governano avvisi e arresti."""

    # katharoseq_if_available: le soglie di S11, sulle letture senza chimere;
    # none: nessuna soglia di profondita', resta min_reads_final.
    min_reads_mode: Literal["katharoseq_if_available", "none"] = Field(description=(
        "Da dove viene la profondita' minima dei campioni: katharoseq_if_available "
        "(le soglie dei controlli positivi) oppure none (nessuna soglia di "
        "profondita')."))
    # Frazione massima delle letture dei biologici rimossa come contaminante (S12).
    # Condizionale: senza valore e' nullo (defaults.CONDIZIONALI).
    max_frac_contaminant: Frazione | None = Field(None, description=(
        "Frazione massima delle letture dei biologici che la decontaminazione puo' "
        "rimuovere (E-S12-02). Serve solo con controlli negativi dichiarati "
        "(ctrl.blank_values o ctrl.blank_override_values non vuoto)."))
    # Per campione, sulle letture dell'oggetto finale (S13): sotto, il campione
    # esce dall'oggetto finale.
    min_reads_final: InteroNonNegativo = Field(description=(
        "Letture minime di un campione nell'oggetto finale: sotto, il campione esce."))
    # Per l'insieme dei campioni finali: letture finali su letture senza chimere
    # (S14, E-S14-01).
    min_frac_reads_retained: Frazione = Field(description=(
        "Frazione minima di letture trattenute dall'insieme dei campioni finali, "
        "rispetto a quelle senza chimere (E-S14-01)."))
    max_asv_count: InteroPositivo = Field(description=(
        "Varianti distinte oltre le quali la tabella delle sequenze non si costruisce "
        "(E-S5-01): secondo la memoria disponibile."))
    # Controlli sul risultato del filtro (S2), applicati ai campioni biologici
    # e ai controlli positivi, non ai negativi: vedi steps/s02_filter.py.
    max_zeroed_samples: InteroNonNegativo = Field(description=(
        "Campioni biologici o controlli positivi che il filtro puo' lasciare senza "
        "letture (E-S2-01)."))
    max_frac_lost_filter: Frazione = Field(description=(
        "Perdita media massima di letture nel filtro, sui biologici e sui controlli "
        "positivi (E-S2-02)."))
    max_frac_short_reads: Frazione = Field(description=(
        "Frazione massima di letture piu' corte di filter.truncLen, per classe "
        "(E-S0-09, E-S1-02)."))
    warn_frac_chimeric: Frazione = Field(description=(
        "Frazione di letture chimeriche che da' un avviso; non puo' superare "
        "qc.stop_frac_chimeric."))
    stop_frac_chimeric: Frazione = Field(description=(
        "Frazione di letture chimeriche che ferma la rimozione delle chimere."))
    # Frazione minima di varianti con il phylum assegnato (S8, E-S8-02): vedi
    # steps/s08_taxonomy.py per le classi a cui si applica.
    min_frac_phylum: Frazione = Field(description=(
        "Frazione minima di varianti con il phylum assegnato (E-S8-02)."))
    # Letture ispezionate per file dai gate che leggono le sequenze.
    head_reads: InteroPositivo = Field(description=(
        "Letture ispezionate per file dalla validazione iniziale."))
    max_primer_hit_frac: Frazione = Field(description=(
        "Frazione di letture che iniziano con il primer oltre la quale lo si "
        "considera presente (E-S0-10)."))
    min_motif_frac: Frazione = Field(description=(
        "Frazione minima di letture che devono portare il motivo conservato, sulla "
        "mediana dei campioni biologici e dei controlli positivi (E-S0-16)."))
    max_frac_untrimmed: Frazione | None = Field(None, description=(
        "Taglio cercato: quota di letture senza il primer, sui campioni biologici e "
        "sui controlli positivi, oltre la quale il filtro da' un avviso (E-S2-04)."))
    # Obbligatori: dipendono dalla regione amplificata. Il motivo puo' essere
    # nullo, per una regione senza un motivo noto: G10 verifica allora la sola
    # assenza del primer. Anche il primer puo' essere nullo, se non e' noto:
    # G10 non cerca allora alcun primer, e lo dichiara.
    primer_sequence: SequenzaIupac | None = Field(description=(
        "Primer di amplificazione forward, in codici IUPAC. Nullo: la validazione "
        "non cerca alcun primer in testa alle letture."))
    conserved_motif: Annotated[Regex, AfterValidator(_motivo_in_maiuscolo)] | None = Field(
        description=(
            "Motivo conservato atteso all'inizio della regione amplificata, come "
            "espressione regolare sulle basi in maiuscolo. Nullo: si verifica la sola "
            "assenza del primer."))

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


class Pair(_Gruppo):
    """Modalita' appaiata: la lettura reverse e la fusione delle coppie.

    Il gruppo esiste solo con ``io.layout: paired-end``: in modalita' single-end
    ogni sua chiave e' respinta e il gruppo non compare nella configurazione
    risolta. I parametri senza predefinito sono condizionali della modalita'
    (defaults.PARAMETRI_DELLE_COPPIE): senza valore sono nulli, e che siano
    dichiarati con la modalita' appaiata lo verifica la radice.
    """

    file_markers: Annotated[list[StringaNonVuota], Field(min_length=2, max_length=2)] | None = (
        Field(None, description=(
            "I due marcatori che distinguono nel nome i file di una coppia, prima quello "
            "del file forward e poi quello del reverse: il file reverse ha il nome del "
            "file forward con l'ultima occorrenza del primo marcatore sostituita dal "
            "secondo.")))
    truncLen_reverse: InteroPositivo | None = Field(None, description=(
        "Lunghezza a cui si troncano le letture reverse; quelle piu' corte sono "
        "scartate, con la loro compagna."))
    trimLeft_reverse: InteroNonNegativo | None = Field(None, description=(
        "Basi tolte in testa a ogni lettura reverse: la lunghezza del primer reverse "
        "se le letture lo contengono ancora, zero se e' gia' stato rimosso."))
    maxEE_reverse: RealePositivo = Field(d.PAIR_MAXEE_REVERSE, description=(
        "Errori attesi massimi per lettura reverse (dada2::filterAndTrim, maxEE)."))
    truncLen_shortfall_warn_reverse: InteroNonNegativo | None = Field(None, description=(
        "Di quante basi pair.truncLen_reverse puo' stare sotto la lettura reverse piu' "
        "corta osservata prima che il profilo lo segnali (E-S1-01)."))
    primer_sequence_reverse: SequenzaIupac | None = Field(None, description=(
        "Primer di amplificazione reverse, in codici IUPAC, scritto da 5' a 3' come "
        "compare in testa alle letture reverse. Nullo: non si cerca."))
    conserved_motif_reverse: Annotated[Regex, AfterValidator(_motivo_in_maiuscolo)] | None = (
        Field(None, description=(
            "Motivo conservato con cui iniziano le letture reverse dopo le basi tolte in "
            "testa, come espressione regolare sulle basi in maiuscolo. Nullo: non si "
            "cerca.")))
    max_primer_hit_frac_reverse: Frazione | None = Field(None, description=(
        "Frazione di letture reverse che iniziano con il primer reverse oltre la quale "
        "lo si considera presente (E-S0-25)."))
    min_overlap: InteroPositivo = Field(d.PAIR_MIN_OVERLAP, description=(
        "Sovrapposizione minima fra le due letture perche' una coppia si fonda "
        "(dada2::mergePairs, minOverlap)."))
    max_mismatch: InteroNonNegativo = Field(d.PAIR_MAX_MISMATCH, description=(
        "Discordanze ammesse nella sovrapposizione (dada2::mergePairs, maxMismatch)."))
    trim_overhang: StrictBool = Field(d.PAIR_TRIM_OVERHANG, description=(
        "Le code di una lettura oltre l'inizio dell'altra si tolgono dalla sequenza "
        "fusa (dada2::mergePairs, trimOverhang): il solo valore ammesso e' true."))
    min_merged_frac_warn: Frazione | None = Field(None, description=(
        "Quota di coppie fuse, sui campioni biologici e sui controlli positivi, sotto "
        "la quale la fusione da' un avviso (E-SM-02); non puo' stare sotto "
        "pair.min_merged_frac_stop."))
    min_merged_frac_stop: Frazione | None = Field(None, description=(
        "Quota di coppie fuse sotto la quale la fusione si ferma (E-SM-01)."))
    len_min: InteroPositivo | None = Field(None, description=(
        "Lunghezza minima ammessa delle sequenze fuse."))
    len_max: InteroPositivo | None = Field(None, description=(
        "Lunghezza massima ammessa delle sequenze fuse."))
    expected_merged_len: InteroPositivo | None = Field(None, description=(
        "Lunghezza attesa delle sequenze fuse: la mediana osservata si confronta con "
        "questa (E-SM-03). Nullo: non e' nota, e il confronto non si fa."))

    @model_validator(mode="after")
    def _coerenza(self) -> Pair:
        """Respinge le combinazioni che non possono valere insieme."""
        if not self.trim_overhang:
            raise ValueError(
                "trim_overhang: false non e' realizzato: le code oltre l'inizio dell'altra "
                "lettura si tolgono sempre. Imposta true")
        if self.file_markers is not None and self.file_markers[0] == self.file_markers[1]:
            raise ValueError(
                f"file_markers: i due marcatori sono uguali ({self.file_markers[0]!r}): "
                "non distinguono il file forward dal reverse")
        if (self.truncLen_reverse is not None and self.trimLeft_reverse is not None
                and self.trimLeft_reverse >= self.truncLen_reverse):
            raise ValueError(
                f"trimLeft_reverse ({self.trimLeft_reverse}) deve essere minore di "
                f"truncLen_reverse ({self.truncLen_reverse}): il troncamento si misura "
                "sulla lettura intera, e non resterebbe alcuna base")
        if (self.truncLen_reverse is not None
                and self.truncLen_shortfall_warn_reverse is not None
                and self.truncLen_shortfall_warn_reverse >= self.truncLen_reverse):
            raise ValueError(
                f"truncLen_shortfall_warn_reverse ({self.truncLen_shortfall_warn_reverse}) "
                f"deve essere minore di truncLen_reverse ({self.truncLen_reverse}): una "
                "tolleranza grande quanto il troncamento non potrebbe mai essere superata")
        if (self.min_merged_frac_warn is not None and self.min_merged_frac_stop is not None
                and self.min_merged_frac_warn < self.min_merged_frac_stop):
            raise ValueError(
                f"min_merged_frac_warn ({self.min_merged_frac_warn}) non puo' stare sotto "
                f"min_merged_frac_stop ({self.min_merged_frac_stop}): l'esecuzione si "
                "fermerebbe prima di dare l'avviso")
        if (self.len_min is not None and self.len_max is not None
                and self.len_min > self.len_max):
            raise ValueError(
                f"len_min ({self.len_min}) supera len_max ({self.len_max}): l'intervallo "
                "delle lunghezze ammesse sarebbe vuoto")
        return self


class Retry(_Gruppo):
    """Nuovi tentativi sulle fasi fallite."""

    enabled: StrictBool = Field(d.RETRY_ENABLED, description=(
        "Se ritentare in automatico gli errori ammessi."))
    max_attempts: InteroPositivo = Field(d.RETRY_MAX_ATTEMPTS, description=(
        "Numero totale di tentativi, compreso il primo."))
    # Codici di errore per cui e' ammesso un nuovo tentativo automatico.
    # L'elenco e' chiuso per scelta metodologica e non per una proprieta'
    # del dataset: si ritenta solo dove l'azione correttiva non modifica
    # alcuna assunzione dell'analisi.
    whitelist: list[StringaNonVuota] = Field(
        default_factory=lambda: list(d.RETRY_WHITELIST), description=(
            "Codici per cui e' ammesso un nuovo tentativo automatico: solo dove "
            "l'azione correttiva non modifica alcuna assunzione dell'analisi."))


class Out(_Gruppo):
    """Forma degli artefatti prodotti."""

    # Insiemi chiusi: vanno estesi quando un'altra forma viene realizzata.
    serialization: Literal["rds"] = Field(d.OUT_SERIALIZATION, description=(
        "Formato dell'oggetto finale: RDS."))
    asv_id_scheme: Literal["abundance_rank"] = Field(d.OUT_ASV_ID_SCHEME, description=(
        "Come si nominano le varianti: per abbondanza decrescente nei campioni "
        "biologici."))
    taxa_are_rows: StrictBool = Field(d.OUT_TAXA_ARE_ROWS, description=(
        "Orientamento della tabella dei conteggi: varianti sulle righe, il solo "
        "realizzato."))

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
    export_flat: StrictBool = Field(d.OUT_EXPORT_FLAT, description=(
        "Se scrivere accanto all'oggetto finale gli stessi dati in testo."))
    # Da dove vengono gli identificativi dei campioni nell'oggetto integrato.
    sample_id_source: Literal["accession"] = Field(d.OUT_SAMPLE_ID_SOURCE, description=(
        "Identificativo dei campioni nell'oggetto: la chiave del campione."))
    # Colonne dei metadati portate nell'oggetto integrato oltre a quelle
    # dell'inventario, con il loro nome originale: dalla tabella campioni di
    # studio e dal file di arricchimento. Nell'oggetto prendono un nome
    # sintattico, e la corrispondenza resta in 10_phyloseq/colonne_metadati.tsv.
    # Obbligatori; elenchi vuoti se non si porta alcuna colonna.
    study_columns: list[StringaNonVuota] = Field(description=(
        "Colonne della tabella di studio da portare nell'oggetto oltre a quelle "
        "dell'inventario. Elenco vuoto: nessuna."))
    batch_columns: list[StringaNonVuota] = Field(description=(
        "Colonne del file del lotto da portare nell'oggetto. Elenco vuoto: nessuna."))


class Run(_Gruppo):
    """Parametri dell'esecuzione: ambiente, parallelismo, riproducibilità."""

    # Facoltativo: l'immagine con cui si dichiara di eseguire. Non sceglie
    # l'immagine, la registra nella provenienza.
    container: RiferimentoImmagine | None = Field(None, description=(
        "Immagine con cui si dichiara di eseguire, ancorata per digest: e' registrata "
        "nella provenienza, non sceglie l'immagine."))
    # Il file di blocco delle versioni R. E' referenziato dal Dockerfile e
    # dagli script R, quindi l'ambiente sarebbe riproducibile comunque; entra
    # qui perche' cosi' finisce nella configurazione risolta e nel suo digest,
    # e due esecuzioni diventano confrontabili anche rispetto alle versioni R.
    lockfile: StringaNonVuota = Field(d.RUN_LOCKFILE, description=(
        "File di blocco delle versioni dei pacchetti R."))
    # Entro gli interi di R, a cui il seme viene passato.
    seed: Annotated[int, Field(strict=True, ge=-2147483647, le=2147483647)] = Field(
        d.RUN_SEED, description="Seme dei generatori di numeri casuali.")
    # Nullo per difetto, cioe' automatico: si usano i processori utilizzabili
    # dal processo, contati all'uso (thread_effettivi). Il valore ricavato dalla
    # macchina non entra nella configurazione ne' nel suo digest, che resta lo
    # stesso su macchine diverse; un valore fisso per difetto fermerebbe G14
    # sulle macchine piu' piccole.
    threads: InteroPositivo | None = Field(None, description=(
        "Processori da usare. Nullo: quelli utilizzabili dal processo."))
    # Tempo massimo, in secondi, di ogni processo R di una fase. Nullo per
    # difetto: nessun limite, perche' la durata di una fase dipende dal dataset
    # e dalla macchina. Allo scadere il ponte uccide il processo e il suo
    # gruppo (E-R-05). Non incide sui risultati: o la fase si conclude, o si
    # ferma.
    r_timeout_s: InteroPositivo | None = Field(None, description=(
        "Tempo massimo in secondi di ogni processo R (E-R-05). Nullo: nessun limite."))
    batch_size: InteroPositivo = Field(d.RUN_BATCH_SIZE, description=(
        "Campioni elaborati insieme: governa la memoria richiesta, non i risultati."))
    # Se conservare le letture filtrate da S2 a esecuzione conclusa. Con false
    # si rimuovono solo quando tutte le fasi sono concluse, e la rimozione e'
    # registrata: una ripresa non la scambia per un artefatto perso.
    keep_filtered_fastq: StrictBool = Field(d.RUN_KEEP_FILTERED_FASTQ, description=(
        "Se conservare le letture filtrate a esecuzione conclusa."))
    # La regola rigorosa sulla provenienza. Con true l'esecuzione parte solo da
    # un repository git leggibile, senza modifiche non committate al codice, e
    # in un ambiente R che corrisponde a run.lockfile; l'impronta del sorgente,
    # quella del file di blocco e l'immagine dichiarata entrano nell'impronta di
    # ogni fase. Vedi runner/provenienza.py per cio' che e' verificato e cio'
    # che e' solo dichiarato.
    strict_provenance: StrictBool = Field(d.RUN_STRICT_PROVENANCE, description=(
        "Regola rigorosa sulla provenienza: con true l'esecuzione parte solo da un "
        "clone senza modifiche non committate al codice e in un ambiente R identico "
        "al file di blocco."))


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
    pair: Pair = Field(default_factory=Pair)
    retry: Retry = Field(default_factory=Retry)
    out: Out
    run: Run = Field(default_factory=Run)

    # Analisi ecologiche a valle.
    norm: Norm = Field(default_factory=Norm)
    glom: Glom = Field(default_factory=Glom)
    beta: Beta = Field(default_factory=Beta)
    ord: Ord = Field(default_factory=Ord)
    stat: Stat = Field(default_factory=Stat)

    @model_validator(mode="after")
    def _condizionali_richiesti(self) -> Config:
        """Respinge un condizionale non dichiarato quando la sua condizione e' vera.

        I campi condizionali hanno un valore di assenza perche' si possano
        omettere: che siano dichiarati quando servono lo garantisce questo
        controllo, sui campi ricevuti in ingresso (``model_fields_set``).
        """
        dichiarati = {
            f"{gruppo}.{nome}": getattr(getattr(self, gruppo), nome)
            for gruppo in type(self).model_fields
            for nome in getattr(self, gruppo).model_fields_set
        }
        mancanti = condizionali_mancanti(dichiarati)
        if self.tax.try_rc is None and not self.idtaxa:
            # Obbligatorio con il bayesiano ingenuo, il classificatore di chi
            # non ne dichiara uno: il campo ammette il nullo solo perche' con
            # IdTaxa non esiste.
            mancanti = ["tax.try_rc", *mancanti]
        if mancanti:
            raise ValueError(_voce_mancanti(mancanti))
        return self

    @property
    def appaiata(self) -> bool:
        """Vero in modalita' appaiata (``io.layout: paired-end``)."""
        return self.io.layout == d.LAYOUT_APPAIATO

    @model_validator(mode="after")
    def _coppie_solo_in_modalita_appaiata(self) -> Config:
        """Respinge un parametro di ``pair`` dichiarato in modalita' single-end.

        Non avrebbe effetto, e non puo' essere accettato in silenzio: in
        modalita' single-end il gruppo non esiste nella configurazione risolta.
        """
        if not self.appaiata and self.pair.model_fields_set:
            raise ValueError(
                "{}: {} {} ma io.layout non e' {}: i parametri del gruppo pair valgono "
                "solo in modalita' appaiata. Dichiara io.layout: \"{}\", oppure togli il "
                "gruppo pair".format(
                    INTESTAZIONE_COPPIE,
                    ", ".join(f"pair.{nome}" for nome in sorted(self.pair.model_fields_set)),
                    "e' dichiarato" if len(self.pair.model_fields_set) == 1 else "sono dichiarati",
                    d.LAYOUT_APPAIATO, d.LAYOUT_APPAIATO))
        if self.appaiata and self.pair.trimLeft_reverse is not None:
            # Con la modalita' appaiata il troncamento forward deve lasciare basi
            # dopo il taglio, come il reverse (Pair._coerenza).
            if self.filter.trimLeft >= self.filter.truncLen:
                raise ValueError(
                    f"{INTESTAZIONE_COPPIE}: filter.trimLeft ({self.filter.trimLeft}) deve "
                    f"essere minore di filter.truncLen ({self.filter.truncLen}): non "
                    "resterebbe alcuna base da fondere")
        return self

    @property
    def taglio_cercato(self) -> bool:
        """Vero se il primer si toglie cercandolo (``filter.primer_removal: search``)."""
        return self.filter.primer_removal == d.PRIMER_CERCATO

    @property
    def idtaxa(self) -> bool:
        """Vero se il classificatore e' IdTaxa (``tax.classifier: idtaxa``)."""
        return self.tax.classifier == d.CLASSIFICATORE_IDTAXA

    def chiavi_assenti(self) -> frozenset[str]:
        """I parametri che con le modalita' di questa configurazione non
        esistono: quelli del taglio cercato con il taglio a posizione fissa, e
        quelli del classificatore non scelto. Non sono nella configurazione
        risolta ne' nell'impronta di alcuna fase.
        """
        assenti: set[str] = set()
        if not self.taglio_cercato:
            assenti.update(("filter.primer_removal", *d.PARAMETRI_DEL_TAGLIO_CERCATO))
        assenti.update(d.SOLO_DEL_BAYESIANO if self.idtaxa else d.PARAMETRI_DI_IDTAXA)
        return frozenset(assenti)

    @model_validator(mode="after")
    def _modalita_coerenti(self) -> Config:
        """Respinge i parametri di una modalita' non dichiarata, e le
        combinazioni che il taglio cercato non ammette.

        Un parametro del taglio cercato senza ``filter.primer_removal: search``,
        o del classificatore non scelto, non avrebbe effetto: non puo' essere
        accettato in silenzio, perche' non esiste nella configurazione risolta.
        """
        dichiarati = {
            f"{gruppo}.{nome}"
            for gruppo in ("filter", "qc", "tax", "pair")
            for nome in getattr(self, gruppo).model_fields_set
        }

        def elenco(chiavi: Iterable[str]) -> str:
            trovate = [c for c in chiavi if c in dichiarati]
            return ", ".join(trovate) + (
                " e' dichiarato" if len(trovate) == 1 else " sono dichiarati")

        if not self.taglio_cercato:
            if dichiarati & set(d.PARAMETRI_DEL_TAGLIO_CERCATO):
                raise ValueError(
                    f"{INTESTAZIONE_TAGLIO}: {elenco(d.PARAMETRI_DEL_TAGLIO_CERCATO)} ma "
                    f"filter.primer_removal non e' {d.PRIMER_CERCATO}: valgono solo con il "
                    "primer cercato lungo la lettura. Dichiara filter.primer_removal: "
                    f"\"{d.PRIMER_CERCATO}\", oppure toglili")
        else:
            mancano = []
            if self.filter.trimLeft != 0:
                mancano.append(
                    f"filter.trimLeft deve essere 0 (e' {self.filter.trimLeft}): il "
                    "taglio lo fa la ricerca del primer")
            if self.qc.primer_sequence is None:
                mancano.append("qc.primer_sequence non puo' essere nullo: e' il primer da cercare")
            if self.appaiata:
                if self.pair.trimLeft_reverse not in (None, 0):
                    mancano.append(
                        f"pair.trimLeft_reverse deve essere 0 (e' {self.pair.trimLeft_reverse})")
                if (self.pair.primer_sequence_reverse is None
                        and "primer_sequence_reverse" in self.pair.model_fields_set):
                    mancano.append(
                        "pair.primer_sequence_reverse non puo' essere nullo: e' il primer "
                        "da cercare nelle letture reverse")
            lunghezze = [len(p) for p in (
                self.qc.primer_sequence,
                self.pair.primer_sequence_reverse if self.appaiata else None) if p]
            if (lunghezze and self.filter.primer_max_mismatch is not None
                    and self.filter.primer_max_mismatch >= min(lunghezze)):
                mancano.append(
                    f"filter.primer_max_mismatch ({self.filter.primer_max_mismatch}) deve "
                    f"essere minore della lunghezza del primer ({min(lunghezze)}): altrimenti "
                    "ogni lettura lo porterebbe")
            if mancano:
                raise ValueError(
                    f"{INTESTAZIONE_TAGLIO}: con filter.primer_removal "
                    f"{d.PRIMER_CERCATO} " + "; ".join(mancano))
        if self.idtaxa:
            # tax.min_boot ha un predefinito: conta se e' stato scritto.
            if dichiarati & set(d.SOLO_DEL_BAYESIANO):
                raise ValueError(
                    f"{INTESTAZIONE_CLASSIFICATORE}: {elenco(d.SOLO_DEL_BAYESIANO)} ma "
                    f"tax.classifier e' {d.CLASSIFICATORE_IDTAXA}: sono parametri di "
                    f"{d.TAX_CLASSIFIER}. Toglili: IdTaxa usa tax.idtaxa_trained, "
                    "tax.idtaxa_trained_md5, tax.idtaxa_threshold e tax.idtaxa_strand")
        elif dichiarati & set(d.PARAMETRI_DI_IDTAXA):
            raise ValueError(
                f"{INTESTAZIONE_CLASSIFICATORE}: {elenco(d.PARAMETRI_DI_IDTAXA)} ma "
                f"tax.classifier non e' {d.CLASSIFICATORE_IDTAXA}: valgono solo con "
                f"IdTaxa. Dichiara tax.classifier: \"{d.CLASSIFICATORE_IDTAXA}\", oppure "
                "toglili")
        return self

    @model_serializer(mode="wrap")
    def _senza_le_modalita_assenti(self, serializza: Any) -> dict[str, Any]:
        """La configurazione come mappa, senza cio' che con le sue modalita' non
        esiste: in modalita' single-end la chiave della modalita' e il gruppo
        ``pair``; con il taglio a posizione fissa la chiave del modo e i
        parametri del taglio cercato; i parametri del classificatore non scelto.

        E' la forma su cui si calcolano il digest e le impronte delle fasi e
        che finisce nella configurazione registrata: una configurazione che non
        dichiara una modalita' ha esattamente le chiavi che aveva senza.
        """
        mappa = serializza(self)
        if not self.appaiata:
            mappa.pop("pair", None)
            if isinstance(mappa.get("io"), dict):
                mappa["io"].pop("layout", None)
        for chiave in self.chiavi_assenti():
            gruppo, nome = chiave.split(".")
            if isinstance(mappa.get(gruppo), dict):
                mappa[gruppo].pop(nome, None)
        return mappa


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
        ammessi = str(contesto.get("expected", "?")).replace(" or ", " oppure ")
        return f"valore non ammesso; ammessi: {ammessi}", True

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


#: Come inizia un problema di coerenza della modalita' appaiata sollevato dalla
#: radice: G15 lo riconosce da qui per attribuirgli il proprio codice (E-G15-15).
INTESTAZIONE_COPPIE: Final = "modalita' appaiata"

#: Come iniziano i problemi del taglio dei primer cercati e del classificatore
#: sollevati dalla radice: G15 li riconosce da qui (E-G15-16, E-G15-17).
INTESTAZIONE_TAGLIO: Final = "taglio dei primer"
INTESTAZIONE_CLASSIFICATORE: Final = "classificatore tassonomico"

#: Come inizia il problema che elenca i parametri obbligatori non dichiarati:
#: G15 lo riconosce da qui per attribuirgli il proprio codice (E-G15-10).
INTESTAZIONE_MANCANTI: Final = "parametri obbligatori non dichiarati"


def _valori_dichiarati(dati: dict[str, Any]) -> dict[str, Any]:
    """I parametri presenti nel dizionario, come ``gruppo.chiave`` -> valore."""
    return {
        f"{gruppo}.{nome}": valore
        for gruppo, contenuto in dati.items() if isinstance(contenuto, dict)
        for nome, valore in contenuto.items()
    }


def condizionali_mancanti(valori: dict[str, Any]) -> list[str]:
    """I condizionali richiesti e non dichiarati, dati i parametri dichiarati.

    Un condizionale manca se la sua condizione e' vera e la chiave e' assente,
    oppure, per quelli che non hanno un valore di assenza, e' nulla. La
    condizione si valuta sui parametri dichiarati: uno assente vale come la
    sua assenza, e cio' che ne dipende non e' richiesto.
    """
    nulli_se_omessi = (
        d.CONDIZIONALI_NULLI_SE_OMESSI | d.COPPIE_NULLE_SE_OMESSE | d.MODALITA_NULLE_SE_OMESSE)
    return [
        chiave for chiave in (*d.CONDIZIONALI, *d.CONDIZIONALI_RESPINTI)
        if d.condizione_vera(chiave, valori) and (
            chiave not in valori
            or (chiave in nulli_se_omessi and valori[chiave] is None)
        )
    ]


def obbligatori_mancanti(dati: dict[str, Any]) -> list[str]:
    """I parametri obbligatori assenti dal dizionario: quelli di
    :data:`defaults.OBBLIGATORI`, e i condizionali la cui condizione e' vera.

    Conta la presenza della chiave, non il valore: un parametro dichiarato
    nullo o vuoto e' dichiarato (salvo i condizionali senza un valore di
    assenza, per i quali il nullo vale come l'omissione).
    """
    valori = _valori_dichiarati(dati)
    # Con IdTaxa i parametri del bayesiano ingenuo non si dichiarano.
    non_chiesti = (
        set(d.SOLO_DEL_BAYESIANO) if valori.get("tax.classifier") == d.CLASSIFICATORE_IDTAXA
        else set())
    # tax.try_rc ammette il nullo solo perche' con IdTaxa non esiste: con il
    # bayesiano ingenuo il nullo vale come l'omissione.
    def manca(chiave: str) -> bool:
        return chiave not in valori or (chiave == "tax.try_rc" and valori[chiave] is None)

    return [c for c in d.OBBLIGATORI if manca(c) and c not in non_chiesti] + (
        condizionali_mancanti(valori))


def _voce_mancanti(mancanti: list[str]) -> str:
    """La voce che elenca i parametri obbligatori non dichiarati."""
    voce = (
        f"{INTESTAZIONE_MANCANTI} ({len(mancanti)}): {', '.join(mancanti)}. Descrivono il "
        "dataset o le scelte di analisi che ne dipendono, e non hanno un valore "
        "predefinito: vanno dichiarati nel file di configurazione. Che cosa significa "
        "ciascuno e come si sceglie e' in GUIDA.md, sezione 6, e in "
        "config/config.example.yaml"
    )
    condizionali = [
        c for c in mancanti if c in d.CONDIZIONALI or c in d.CONDIZIONALI_RESPINTI]
    if condizionali:
        voce += ". Richiesti da una condizione: " + "; ".join(
            f"{c} perche' {d.testo_condizione(c)}" for c in condizionali
        )
    return voce


def valida(dati: dict[str, Any], origine: str | None = None) -> Config:
    """Valida un dizionario già caricato.

    Solleva :class:`ErroreConfigurazione` con l'elenco completo dei problemi. I
    parametri obbligatori non dichiarati vi compaiono in una sola voce, tutti
    insieme: chi compila la configurazione per un dataset nuovo deve vedere
    l'elenco intero, non scoprirli uno per esecuzione.
    """
    # Un gruppo con tutte le chiavi in commento arriva come nullo: vale come
    # un gruppo vuoto, cosi' il rifiuto nomina i parametri che gli mancano.
    dati = {gruppo: {} if valori is None else valori for gruppo, valori in dati.items()}
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
            if not chiave and INTESTAZIONE_MANCANTI in str(e.get("msg", "")):
                # Il rifiuto dei condizionali del modello: e' gia' nell'elenco.
                continue
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
