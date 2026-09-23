"""Catalogo dei codici di errore della pipeline.

Ogni codice porta due cose che non possono stare altrove senza perdersi:

* un **messaggio operativo**, che dice cosa fare e non solo cosa è andato
  storto. Un messaggio che si limita a constatare il guasto lascia all'utente
  il lavoro di capire da dove ripartire, e quel lavoro va fatto una volta qui
  invece che ogni volta davanti al terminale;
* una **categoria di gestione**, che stabilisce cosa fa la pipeline quando il
  codice si presenta.

Il catalogo è la sola fonte di verità dei codici. I gate di validazione
(:mod:`amplicon16s.gates.g01_g15`) usano i propri codici ``E-G15-*``, che sono
registrati qui insieme a tutti gli altri: un solo insieme di codici, un solo
posto dove leggere cosa significano.

**Numerazione dei codici di fase.** Per S0 i codici da ``E-S0-01`` a
``E-S0-14`` corrispondono ai gate da ``G01`` a ``G14``. ``E-S0-15`` non
corrisponde a un gate: è l'avviso che G08 emette sulla composizione delle
piastre, distinto dal proprio codice di errore perché un avviso e un arresto
non possono condividere una categoria di gestione.

**Il retry automatico è un elenco chiuso.** Solo quattro codici lo ammettono, e
sono gli stessi dichiarati in ``retry.whitelist``. Il criterio non è la gravità
ma la natura dell'azione correttiva: ridurre la dimensione di un lotto o
aumentare i dati di una stima non modifica alcuna assunzione metodologica,
mentre spostare una soglia è una decisione scientifica e non può essere presa
da un programma. Un test verifica nei due versi che l'elenco resti quello.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import StrEnum
from typing import Final

__all__ = [
    "CATALOGO",
    "Categoria",
    "VoceCatalogo",
    "codici_con_retry",
    "codici_di_fase",
    "voce",
]


class Categoria(StrEnum):
    """Come la pipeline reagisce a un codice di errore."""

    #: L'esecuzione si ferma, decide l'operatore.
    REVISIONE_UMANA = "revisione_umana"

    #: Si ritenta entro il limite di retry.max_attempts.
    RETRY_AUTOMATICO = "retry_automatico"

    #: Si ritenta e, se fallisce ancora, l'esecuzione si ferma.
    RETRY_POI_REVISIONE = "retry_poi_revisione"

    #: Si prosegue con un comportamento di ripiego, registrandolo.
    DEGRADAZIONE_AUTOMATICA = "degradazione_automatica"

    @property
    def ammette_retry(self) -> bool:
        return self in (Categoria.RETRY_AUTOMATICO, Categoria.RETRY_POI_REVISIONE)

    @property
    def ferma_esecuzione(self) -> bool:
        return self in (Categoria.REVISIONE_UMANA, Categoria.RETRY_POI_REVISIONE)


@dataclass(frozen=True)
class VoceCatalogo:
    """Un codice di errore con il suo significato e la sua gestione."""

    codice: str
    #: Fase a cui il codice appartiene: ``S0``…``S14`` per le fasi della
    #: pipeline, ``G15`` per il gate di coerenza della configurazione.
    fase: str
    #: Che cosa è andato storto, in una riga.
    sintesi: str
    #: Che cosa deve fare l'operatore, o che cosa fa la pipeline da sé.
    azione: str
    categoria: Categoria

    @property
    def ammette_retry(self) -> bool:
        return self.categoria.ammette_retry

    @property
    def messaggio(self) -> str:
        """Messaggio completo: la constatazione seguita dall'azione."""
        return f"{self.sintesi} {self.azione}"

    def __str__(self) -> str:
        return f"[{self.codice}] {self.messaggio}"


def _v(
    codice: str,
    fase: str,
    sintesi: str,
    azione: str,
    categoria: Categoria,
) -> VoceCatalogo:
    return VoceCatalogo(codice, fase, sintesi, azione, categoria)


_UMANA = Categoria.REVISIONE_UMANA
_RETRY = Categoria.RETRY_AUTOMATICO
_RETRY_UMANA = Categoria.RETRY_POI_REVISIONE
_DEGRADA = Categoria.DEGRADAZIONE_AUTOMATICA


_VOCI: Final[tuple[VoceCatalogo, ...]] = (
    # ---------------------------------------------------------------- G15 ---
    _v(
        "E-G15-01", "G15",
        "filter.minLen supera filter.truncLen.",
        "E' un'incoerenza nella derivazione dei parametri, non nella "
        "configurazione: segnalala, perche' filter.minLen e' calcolato da "
        "filter.truncLen e non dovrebbe poterlo superare.",
        _UMANA,
    ),
    _v(
        "E-G15-02", "G15",
        "decontam.threshold non e' strettamente compreso fra 0 e 1.",
        "Porta decontam.threshold entro l'intervallo aperto (0, 1): a 0 nessuna "
        "variante verrebbe mai classificata come contaminante, a 1 lo sarebbero "
        "tutte.",
        _UMANA,
    ),
    _v(
        "E-G15-03", "G15",
        "L'intervallo delle lunghezze ammesse per le ASV non e' valido.",
        "Riduci asv.len_tol: con il valore attuale l'intervallo derivato da "
        "filter.truncLen meno filter.trimLeft risulta vuoto o non positivo.",
        _UMANA,
    ),
    _v(
        "E-G15-04", "G15",
        "prev.min_fraction e' fuori dall'intervallo [0, 1].",
        "E' una frazione di campioni: riportala fra 0 e 1.",
        _UMANA,
    ),
    _v(
        "E-G15-05", "G15",
        "qc.warn_frac_chimeric supera qc.stop_frac_chimeric.",
        "Abbassa qc.warn_frac_chimeric sotto qc.stop_frac_chimeric, altrimenti "
        "l'esecuzione si fermerebbe prima di emettere l'avviso.",
        _UMANA,
    ),
    _v(
        "E-G15-06", "G15",
        "retry.max_attempts non e' un numero di tentativi ammissibile.",
        "Impostalo ad almeno 1; per disattivare i nuovi tentativi usa "
        "retry.enabled: false.",
        _UMANA,
    ),
    _v(
        "E-G15-07", "G15",
        "La decontaminazione e' attiva ma ctrl.blank_values e' vuoto.",
        "Elenca in ctrl.blank_values le etichette dei controlli negativi: senza "
        "di essi non e' possibile individuare i contaminanti.",
        _UMANA,
    ),
    _v(
        "E-G15-08", "G15",
        "Un parametro derivato e' stato impostato a mano.",
        "Rimuovilo dalla configurazione: e' calcolato dalla pipeline a partire "
        "dagli altri parametri, e scriverlo significherebbe poterlo mettere in "
        "contraddizione con quelli da cui dipende.",
        _UMANA,
    ),
    _v(
        "E-G15-99", "G15",
        "La configurazione non e' valida parametro per parametro.",
        "Correggi i parametri segnalati nell'elenco che accompagna l'errore: "
        "ciascuno indica la chiave e il valore rifiutato.",
        _UMANA,
    ),
    # ----------------------------------------------------------------- S0 ---
    _v(
        "E-S0-01", "S0",
        "Un file di ingresso non esiste o non e' leggibile.",
        "Verifica che io.fastq_dir e io.assay_table indichino percorsi "
        "esistenti e leggibili dall'utente che esegue la pipeline; se i dati "
        "stanno su un volume montato, controlla che sia montato nel container.",
        _UMANA,
    ),
    _v(
        "E-S0-02", "S0",
        "Una tabella di metadati non si apre o non ha le colonne attese.",
        "Verifica che i nomi dichiarati in meta.sample_id_column, "
        "meta.accession_column, ctrl.column e meta.module_column corrispondano "
        "alle intestazioni effettive delle tabelle: nelle tabelle ISA i nomi "
        "contengono spazi e parentesi e vanno riportati alla lettera.",
        _UMANA,
    ),
    _v(
        "E-S0-03", "S0",
        "Il join fra file e metadati non e' ristretto alla tabella di assay.",
        "Verifica che io.assay_table sia la tabella dell'assay 16S e non quella "
        "di studio, e che meta.accession_column indichi la colonna giusta: sono "
        "comparse righe estranee all'assay.",
        _UMANA,
    ),
    _v(
        "E-S0-04", "S0",
        "Un accession non e' estraibile o e' ambiguo.",
        "Adatta io.accession_regex al formato degli accession presenti in "
        "meta.accession_column: l'espressione attuale non produce una "
        "corrispondenza unica.",
        _UMANA,
    ),
    _v(
        "E-S0-05", "S0",
        "Due o piu' campioni condividono lo stesso accession.",
        "Correggi la tabella di assay o restringi il join prima di rieseguire: "
        "con accession duplicati non e' possibile associare le letture ai "
        "campioni in modo univoco.",
        _UMANA,
    ),
    _v(
        "E-S0-06", "S0",
        "Gli insiemi dei file e dei metadati non coincidono.",
        "Verifica io.fastq_glob e meta.accession_column, poi decidi se i "
        "campioni spaiati vanno esclusi o recuperati: la pipeline non li scarta "
        "da se', perche' sarebbe una scelta sui dati.",
        _UMANA,
    ),
    _v(
        "E-S0-07", "S0",
        "I dati non hanno un layout single-end.",
        "Seleziona il sottoinsieme single-end oppure usa una pipeline "
        "paired-end: questa tratta solo letture singole.",
        _UMANA,
    ),
    _v(
        "E-S0-08", "S0",
        "L'informazione di lotto non e' coerente.",
        "Verifica che le colonne dichiarate in decontam.batch_column ed "
        "err.batch_column esistano nel file indicato da io.batch_table: se non "
        "ci sono, il lotto resterebbe nullo per ogni campione senza che nulla "
        "lo segnali. In alternativa togli io.batch_table e la pipeline "
        "procedera' dichiaratamente senza lotto.",
        _UMANA,
    ),
    _v(
        "E-S0-15", "S0",
        "Una piastra non ha abbastanza controlli negativi.",
        "L'esecuzione prosegue e la piastra viene registrata. La "
        "decontaminazione su di essa sara' fondata su meno bianchi di quanti "
        "decontam.min_blanks ne chieda, quindi i suoi contaminanti saranno "
        "stimati con piu' incertezza: tienine conto leggendo i risultati di "
        "quei campioni. Non e' un errore da correggere a posteriori, e' una "
        "proprieta' di come la piastra e' stata allestita.",
        _DEGRADA,
    ),
    _v(
        "E-S0-09", "S0",
        "filter.truncLen e' incompatibile con le lunghezze osservate.",
        "Abbassa filter.truncLen dopo aver esaminato i profili di lunghezza in "
        "02_qc_profiles: al valore attuale una parte rilevante delle letture "
        "verrebbe scartata.",
        _UMANA,
    ),
    _v(
        "E-S0-10", "S0",
        "Le letture contengono ancora il primer.",
        "Imposta filter.trimLeft alla lunghezza del primer, oppure rimuovilo a "
        "monte: lasciarlo falsa l'inferenza delle varianti.",
        _UMANA,
    ),
    _v(
        "E-S0-11", "S0",
        "Alcuni campioni non ricadono in nessuna categoria di controllo.",
        "Verifica ctrl.column e le etichette in ctrl.blank_values, "
        "ctrl.positive_values e ctrl.biological_values: ogni campione deve "
        "ricadere in una e una sola categoria.",
        _UMANA,
    ),
    _v(
        "E-S0-12", "S0",
        "Il database tassonomico non supera la verifica di integrita'.",
        "Il checksum di tax.ref_fasta non corrisponde a tax.ref_md5: riscarica "
        "il riferimento, oppure aggiorna tax.ref_md5 dopo aver verificato "
        "l'origine del file.",
        _UMANA,
    ),
    _v(
        "E-S0-13", "S0",
        "Un file di letture non ha la struttura attesa.",
        "L'archivio non si decomprime, oppure i record non hanno quattro righe "
        "ciascuno. Riscarica il file dalla sorgente e verificane il checksum "
        "prima di rieseguire: la validazione ispeziona solo le prime "
        "qc.head_reads letture, quindi un guasto oltre quelle non verrebbe "
        "visto qui.",
        _UMANA,
    ),
    _v(
        "E-S0-14", "S0",
        "Le risorse richieste non sono disponibili.",
        "Riduci run.threads al numero di CPU effettivamente utilizzabili, "
        "oppure libera spazio sul volume che ospita io.out_root: la pipeline "
        "scrive artefatti intermedi di dimensione paragonabile ai dati di "
        "ingresso.",
        _UMANA,
    ),
    # ----------------------------------------------------------------- S1 ---
    _v(
        "E-S1-01", "S1",
        "Lo scarto fra filter.truncLen e la lunghezza minima osservata supera "
        "filter.truncLen_shortfall_warn.",
        "L'esecuzione prosegue e lo scarto viene registrato in 02_qc_profiles. "
        "Se lo scarto e' ampio, valuta di abbassare filter.truncLen: i campioni "
        "piu' corti perderanno letture.",
        _DEGRADA,
    ),
    # ----------------------------------------------------------------- S2 ---
    _v(
        "E-S2-01", "S2",
        "Uno o piu' campioni non conservano alcuna lettura dopo il filtro.",
        "Allenta filter.maxEE o filter.truncLen, oppure escludi quei campioni "
        "consapevolmente: proseguire con campioni vuoti falserebbe le fasi "
        "successive.",
        _UMANA,
    ),
    _v(
        "E-S2-02", "S2",
        "La frazione media di letture conservate e' sotto la soglia attesa.",
        "Rivedi filter.maxEE e filter.truncLen sulla base dei profili di "
        "qualita' in 02_qc_profiles: il filtro sta scartando troppo.",
        _UMANA,
    ),
    _v(
        "E-S2-03", "S2",
        "Errore di lettura o file non decomprimibile su un lotto.",
        "Il lotto viene rielaborato con run.batch_size ridotto. Se l'errore "
        "persiste oltre retry.max_attempts, il file e' probabilmente corrotto: "
        "verificane l'integrita' alla sorgente.",
        _RETRY,
    ),
    # ----------------------------------------------------------------- S3 ---
    _v(
        "E-S3-01", "S3",
        "Il modello d'errore non converge.",
        "L'apprendimento viene ritentato con piu' basi (err.nbases). Se non "
        "converge ancora, esamina il lotto indicato da err.batch_column: "
        "potrebbe raccogliere dati eterogenei che vanno separati.",
        _RETRY_UMANA,
    ),
    # ----------------------------------------------------------------- S4 ---
    _v(
        "E-S4-02", "S4",
        "Memoria esaurita durante l'inferenza delle varianti.",
        "L'inferenza viene ritentata con run.batch_size ridotto. Se l'errore "
        "persiste, aumenta la memoria disponibile al container o riduci "
        "run.threads: ogni thread ne consuma una quota.",
        _RETRY,
    ),
    # ----------------------------------------------------------------- S5 ---
    _v(
        "E-S5-01", "S5",
        "Le varianti inferite eccedono la memoria disponibile.",
        "Si ritenta con run.batch_size ridotto. Se il numero di varianti resta "
        "oltre qc.max_asv_count, la decisione e' scientifica: valuta un filtro "
        "piu' severo a monte, anziche' alzare la soglia.",
        _RETRY_UMANA,
    ),
    # ----------------------------------------------------------------- S6 ---
    _v(
        "E-S6-01", "S6",
        "La frazione chimerica supera qc.stop_frac_chimeric.",
        "Verifica la presenza di primer residui e i parametri del gruppo "
        "chimera prima di proseguire: alzare la soglia non risolve il problema, "
        "lo nasconde.",
        _UMANA,
    ),
    _v(
        "E-S6-02", "S6",
        "La frazione chimerica supera qc.warn_frac_chimeric.",
        "L'esecuzione prosegue e la frazione viene registrata in 07_chimera. Se "
        "l'avviso ricorre su piu' lotti, esamina i parametri del gruppo chimera.",
        _DEGRADA,
    ),
    # ----------------------------------------------------------------- S8 ---
    _v(
        "E-S8-02", "S8",
        "La copertura tassonomica e' insufficiente.",
        "Troppe varianti restano senza assegnazione: verifica che tax.ref_fasta "
        "copra la regione amplificata e che tax.min_boot non sia troppo alto.",
        _UMANA,
    ),
    # ----------------------------------------------------------------- S9 ---
    _v(
        "E-S9-01", "S9",
        "Le sequenze superano phylo.max_seqs con la filogenesi attiva.",
        "L'esecuzione si ferma prima di allocare il calcolo. Il codice si "
        "presenta solo se phylo.enabled e' true, cioe' se l'albero e' stato "
        "richiesto esplicitamente: saltarlo in silenzio consegnerebbe un "
        "oggetto privo proprio della componente domandata. Scegli se alzare "
        "phylo.max_seqs, sapendo che il costo cresce rapidamente con il numero "
        "di sequenze, oppure rinunciare all'albero con phylo.enabled: false, "
        "accettando che le analisi che lo richiedono non saranno disponibili.",
        _UMANA,
    ),
    # ---------------------------------------------------------------- S10 ---
    _v(
        "E-S10-01", "S10",
        "Gli slot dell'oggetto integrato non sono coerenti fra loro.",
        "Campioni o varianti non concordano fra gli slot. Conserva 10_phyloseq "
        "e verifica con il manifesto l'integrita' degli artefatti delle fasi "
        "precedenti: l'incoerenza nasce quasi sempre a monte.",
        _UMANA,
    ),
    # ---------------------------------------------------------------- S11 ---
    _v(
        "E-S11-02", "S11",
        "La soglia di profondita' derivata dai controlli positivi non e' "
        "attendibile.",
        "Si prosegue usando qc.min_reads_final come soglia, e il ripiego viene "
        "registrato in 11_controls. Per una soglia derivata, verifica che i "
        "controlli positivi siano sufficienti e conformi.",
        _DEGRADA,
    ),
    _v(
        "E-S11-03", "S11",
        "I controlli positivi non sono conformi.",
        "Verifica ctrl.positive_values e katharoseq.target_taxon: un controllo "
        "positivo anomalo mette in dubbio l'intero lotto, quindi la decisione "
        "di proseguire non puo' essere automatica.",
        _UMANA,
    ),
    # ---------------------------------------------------------------- S12 ---
    _v(
        "E-S12-02", "S12",
        "I contaminanti individuati superano la quota attesa.",
        "Esamina i controlli negativi e decontam.batch_column prima di "
        "accettare il risultato: una quota anomala indica piu' spesso un "
        "raggruppamento sbagliato che una contaminazione reale.",
        _UMANA,
    ),
    # ---------------------------------------------------------------- S13 ---
    _v(
        "E-S13-01", "S13",
        "Una fase e' stata eseguita prima di una da cui dipende.",
        "Riparti da un'esecuzione pulita, oppure usa 'amplicon16s resume', che "
        "ricostruisce l'ordine dagli artefatti gia' prodotti.",
        _UMANA,
    ),
    _v(
        "E-S13-02", "S13",
        "Il filtro di prevalenza lascia alcuni campioni senza varianti.",
        "Abbassa prev.min_fraction o prev.min_count, oppure escludi quei "
        "campioni consapevolmente: un campione vuoto non e' un campione a "
        "prevalenza bassa.",
        _UMANA,
    ),
    # ---------------------------------------------------------------- S14 ---
    _v(
        "E-S14-01", "S14",
        "La validazione finale dell'oggetto non e' superata.",
        "Non usare ps_final.rds: conserva 12_final e il log in 99_logs, che "
        "insieme indicano quale controllo non e' stato superato.",
        _UMANA,
    ),
)


#: Il catalogo, indicizzato per codice.
CATALOGO: Final[dict[str, VoceCatalogo]] = {v.codice: v for v in _VOCI}


def voce(codice: str) -> VoceCatalogo:
    """Restituisce la voce di catalogo, o solleva ``KeyError`` se non esiste.

    Un codice non catalogato è un difetto di programmazione, non una
    condizione da gestire a runtime: fallire subito lo rende visibile.
    """
    try:
        return CATALOGO[codice]
    except KeyError:
        raise KeyError(f"codice non presente nel catalogo degli errori: {codice}") from None


def codici_con_retry() -> frozenset[str]:
    """Codici ammessi al retry automatico. È un elenco chiuso."""
    return frozenset(v.codice for v in _VOCI if v.ammette_retry)


def codici_di_fase(fase: str) -> tuple[str, ...]:
    """Codici appartenenti a una fase, in ordine."""
    return tuple(v.codice for v in _VOCI if v.fase == fase)
