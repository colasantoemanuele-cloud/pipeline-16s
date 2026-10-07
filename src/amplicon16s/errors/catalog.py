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
ma la natura dell'azione correttiva: rileggere un file (E-S2-03), aumentare i
dati di una stima (E-S3-01) o ridurre la dimensione di un lotto (E-S4-02) non
modifica alcuna assunzione metodologica, mentre spostare una soglia è una
decisione scientifica e non può essere presa da un programma. E-S5-01 resta
ammesso, ma la sua fase dichiara nell'errore che ritentare non servirebbe (la
memoria della tabella non dipende dal lotto), e l'esecuzione si ferma. Un test
verifica nei due versi che l'elenco resti quello.

**I codici del ponte verso R.** I codici ``E-R-*`` non appartengono a una
fase ma al ponte (:mod:`amplicon16s.rbridge`) che esegue gli script R di tutte
le fasi. Descrivono ciò che la fase non può dichiarare da sé: un interprete
assente, un processo morto prima di scrivere l'esito, un errore che lo script
non ha ricondotto a un codice. Sono tutti a revisione umana: nessuno di questi
guasti ha un'azione correttiva che non richieda di capirne prima la causa.

**I codici della provenienza.** I codici ``E-PROV-*`` appartengono alla regola
rigorosa sulla provenienza (``run.strict_provenance``), un controllo di avvio
che l'esecutore compie prima di qualunque fase: repository git non leggibile,
codice con modifiche non committate, ambiente R diverso dal file di blocco.
Sono a revisione umana.

**Il codice del grafo.** ``E-GRAFO-01`` segnala una fase avviata prima che le
fasi da cui dipende fossero concluse, qualunque sia la fase. Non va confuso
con ``E-S13-01``, che riguarda la sola precedenza di metodo fra
decontaminazione e filtro di prevalenza: un codice per ogni significato, cosi'
il log non attribuisce a S13 violazioni che riguardano altre fasi.
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
        """Vero per le categorie che prevedono un nuovo tentativo automatico."""
        return self in (Categoria.RETRY_AUTOMATICO, Categoria.RETRY_POI_REVISIONE)

    @property
    def ferma_esecuzione(self) -> bool:
        """Vero per le categorie che, esauriti o esclusi i tentativi, fermano
        l'esecuzione.
        """
        return self in (Categoria.REVISIONE_UMANA, Categoria.RETRY_POI_REVISIONE)


@dataclass(frozen=True)
class VoceCatalogo:
    """Un codice di errore con il suo significato e la sua gestione."""

    codice: str
    #: Fase a cui il codice appartiene: ``S0``…``S14`` per le fasi della
    #: pipeline, ``G15`` per il gate di coerenza della configurazione, ``R``
    #: per il ponte che esegue gli script R, ``GRAFO`` per l'ordine delle fasi.
    fase: str
    #: Che cosa è andato storto, in una riga.
    sintesi: str
    #: Che cosa deve fare l'operatore, o che cosa fa la pipeline da sé.
    azione: str
    categoria: Categoria

    @property
    def ammette_retry(self) -> bool:
        """Vero se la categoria del codice prevede un nuovo tentativo automatico."""
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
    """Costruisce una voce del catalogo; abbrevia la tabella dei codici."""
    return VoceCatalogo(codice, fase, sintesi, azione, categoria)


_UMANA = Categoria.REVISIONE_UMANA
_RETRY = Categoria.RETRY_AUTOMATICO
_RETRY_UMANA = Categoria.RETRY_POI_REVISIONE
_DEGRADA = Categoria.DEGRADAZIONE_AUTOMATICA


_VOCI: Final[tuple[VoceCatalogo, ...]] = (
    # ---------------------------------------------------------------- G15 ---
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
        "E-G15-08", "G15",
        "Un parametro derivato e' stato impostato a mano.",
        "Rimuovilo dalla configurazione: e' calcolato dalla pipeline a partire "
        "dagli altri parametri, e scriverlo significherebbe poterlo mettere in "
        "contraddizione con quelli da cui dipende.",
        _UMANA,
    ),
    _v(
        "E-G15-09", "G15",
        "retry.whitelist contiene un codice che il catalogo non ammette al "
        "retry.",
        "Togli il codice da retry.whitelist. Il retry automatico e' ammesso "
        "solo per i codici la cui azione correttiva non modifica alcuna "
        "assunzione metodologica, e il catalogo li classifica come tali; un "
        "codice assente dal catalogo e' un refuso. L'elenco puo' essere "
        "ristretto, non allargato.",
        _UMANA,
    ),
    _v(
        "E-G15-99", "G15",
        "La configurazione non e' valida parametro per parametro.",
        "Correggi i parametri segnalati nell'elenco che accompagna l'errore: "
        "ciascuno indica la chiave e il valore rifiutato.",
        _UMANA,
    ),
    _v(
        "E-G15-10", "G15",
        "Uno o piu' parametri obbligatori non sono dichiarati.",
        "Dichiara nel file di configurazione tutti i parametri elencati: "
        "descrivono il dataset (formato dei metadati, etichette, primer, "
        "troncamento) e non hanno un valore predefinito, perche' un valore "
        "ereditato da un altro dataset sarebbe quasi certamente sbagliato. Un "
        "parametro non pertinente si dichiara nullo o con un elenco vuoto, non si "
        "omette. config/config.example.yaml li elenca con un esempio per ciascuno.",
        _UMANA,
    ),
    _v(
        "E-G15-11", "G15",
        "Un parametro del file del lotto e' in contraddizione con io.batch_table.",
        "Senza io.batch_table le colonne del file del lotto vanno dichiarate "
        "vuote o nulle (out.batch_columns, err.batch_column, decontam.batch_column, "
        "meta.batch_key_column, meta.batch_module_column): una colonna chiesta a un "
        "file che non c'e' resterebbe vuota senza che nulla lo segnali. Con "
        "io.batch_table indicato, meta.batch_key_column deve nominare la colonna "
        "con la chiave del campione.",
        _UMANA,
    ),
    _v(
        "E-G15-12", "G15",
        "I parametri dei controlli positivi sono in contraddizione fra loro.",
        "Se ctrl.positive_values elenca delle etichette, katharoseq.target_taxon "
        "e katharoseq.cell_count_column devono essere indicati: senza il taxon "
        "atteso e le cellule seminate i controlli positivi non sono valutabili. "
        "La colonna delle cellule si cerca nel file di arricchimento e nella "
        "tabella di studio, ed entra nell'oggetto anche se non e' fra le "
        "colonne richieste. Se il dataset non ha "
        "controlli positivi, dichiara ctrl.positive_values vuoto e i due "
        "parametri nulli.",
        _UMANA,
    ),
    _v(
        "E-G15-13", "G15",
        "Due colonne dei metadati prenderebbero lo stesso nome nell'oggetto.",
        "Le colonne di out.study_columns e out.batch_columns entrano "
        "nell'oggetto con un nome sintattico (minuscole, trattini bassi): due "
        "colonne diverse possono dare lo stesso nome, o quello di una colonna "
        "dell'inventario (accession, sample_name, classe, materiale, posizione, "
        "modulo, piastra, corsa). Togli una delle due dagli elenchi, o "
        "rinominala nella tabella: scoperta da S10, la collisione costerebbe "
        "tutto il calcolo che la precede.",
        _UMANA,
    ),
    _v(
        "E-G15-14", "G15",
        "Un parametro dei metadati e' dichiarato senza quello da cui dipende.",
        "Quattro parametri hanno effetto solo insieme a un altro, e dichiarati da "
        "soli verrebbero ignorati senza che nulla lo segnali: "
        "ctrl.blank_override_values richiede ctrl.blank_override_column (la "
        "colonna in cui cercare i valori); meta.study_sample_id_column richiede "
        "io.study_table (senza tabella di studio il nome del campione si legge "
        "dalla tabella di assay, in meta.sample_id_column); meta.module_regex "
        "e meta.non_surface_positions richiedono meta.module_column (la posizione "
        "da cui derivare il modulo, e in cui riconoscere le non superfici). "
        "Dichiara il parametro mancante, oppure rendi nullo o vuoto quello che "
        "ne dipende.",
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
        "Adatta io.accession_regex a come sono nominati i file: deve trovare nel "
        "nome di ciascuno una e una sola chiave del campione (la corrispondenza "
        "intera, o il primo gruppo di cattura se l'espressione ne ha uno), e la "
        "stessa chiave deve ricavarsi dal valore di meta.accession_column.",
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
        "paired-end: questa tratta solo letture singole. Se un file contiene le "
        "due letture di ogni coppia (lo dicono le intestazioni, con /1 e /2 o con "
        "lo stesso nome due volte), o le sole seconde letture (intestazioni tutte "
        "con /2), estrai le sole letture forward in un file per campione.",
        _UMANA,
    ),
    _v(
        "E-S0-08", "S0",
        "L'informazione di lotto non e' coerente.",
        "Verifica che le colonne dichiarate per il file del lotto "
        "(decontam.batch_column, err.batch_column, meta.batch_module_column, "
        "out.batch_columns) esistano nel file "
        "indicato da io.batch_table, e che il file abbia una e una sola riga per "
        "ogni campione, con piastra e corsa compilate dove sono dichiarate: un "
        "campione senza riga, o con piu' righe, resterebbe senza lotto. In "
        "alternativa togli io.batch_table e dichiara nulle le sue colonne: la "
        "pipeline procedera' dichiaratamente senza lotto.",
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
        "E-S0-19", "S0",
        "Il file del lotto ha righe che non corrispondono ad alcun campione.",
        "L'esecuzione prosegue e quelle righe sono ignorate: il loro numero e' "
        "registrato nel manifesto di S0. E' normale se il file descrive piu' "
        "campioni di quanti l'assay ne contenga. Se invece le righe dovevano "
        "corrispondere a dei campioni, la chiave in meta.batch_key_column e' "
        "scritta in modo diverso dall'accession dei file e della tabella di "
        "assay: confronta le chiavi elencate con 01_input_validation/crosswalk.tsv.",
        _DEGRADA,
    ),
    _v(
        "E-S0-09", "S0",
        "filter.truncLen e' incompatibile con le lunghezze osservate.",
        "Fra le prime qc.head_reads letture dei campioni biologici, o fra quelle "
        "dei controlli positivi (le due classi si giudicano separatamente), "
        "quelle piu' corte di filter.truncLen superano "
        "qc.max_frac_short_reads: il filtro le scarterebbe, non le accorcerebbe. "
        "Abbassa filter.truncLen dopo aver esaminato le lunghezze, oppure alza "
        "qc.max_frac_short_reads se la perdita e' accettata. Poche letture corte "
        "non fermano l'esecuzione.",
        _UMANA,
    ),
    _v(
        "E-S0-10", "S0",
        "Le letture contengono ancora il primer.",
        "Imposta filter.trimLeft alla lunghezza del primer, oppure rimuovilo a "
        "monte: lasciarlo falsa l'inferenza delle varianti. Con filter.trimLeft "
        "maggiore di zero il gate non cerca piu' il primer in testa e verifica il "
        "motivo conservato a partire da quella posizione.",
        _UMANA,
    ),
    _v(
        "E-S0-16", "S0",
        "Le letture non contengono il segnale atteso della regione amplificata.",
        "Il motivo conservato dichiarato in qc.conserved_motif non compare, alla "
        "posizione filter.trimLeft, in una quota sufficiente delle letture dei "
        "campioni da cui ci si attende il segnale. Verifica che i file siano quelli "
        "dell'amplicone dichiarato, che qc.conserved_motif e qc.primer_sequence "
        "corrispondano alla regione amplificata e che filter.trimLeft sia la "
        "lunghezza di cio' che precede il motivo; per una regione senza un motivo "
        "noto dichiara qc.conserved_motif nullo.",
        _UMANA,
    ),
    _v(
        "E-S0-17", "S0",
        "Il dataset non ha i controlli che la decontaminazione o la calibrazione "
        "richiedono.",
        "L'esecuzione prosegue e la mancanza viene registrata. Senza controlli "
        "positivi la soglia di profondita' non puo' essere derivata da una curva; "
        "con meno controlli negativi di decontam.min_blanks i contaminanti non sono "
        "stimabili. Se i controlli esistono ma non sono stati riconosciuti, "
        "correggi ctrl.positive_values, ctrl.blank_values o ctrl.column.",
        _DEGRADA,
    ),
    _v(
        "E-S0-18", "S0",
        "Il segnale atteso non e' stato verificato: nessun motivo conservato "
        "dichiarato.",
        "L'esecuzione prosegue. Con qc.conserved_motif nullo il gate verifica solo "
        "che il primer non sia in testa alle letture, e non puo' distinguere un "
        "file corretto da uno privo della regione amplificata: se un motivo e' "
        "noto, dichiaralo.",
        _DEGRADA,
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
        "L'esecuzione prosegue e lo scarto viene registrato nel manifesto di S1. "
        "Il troncamento e' molto sotto la lettura piu' corta: ogni lettura cede "
        "basi che si potrebbero conservare. Valuta di alzare filter.truncLen "
        "fino al minimo riportato in 02_qc_profiles/riepilogo.json.",
        _DEGRADA,
    ),
    _v(
        "E-S1-02", "S1",
        "Le letture piu' corte di filter.truncLen, misurate su tutte le letture, "
        "superano qc.max_frac_short_reads.",
        "G09, in S0, guarda le sole prime qc.head_reads letture di ogni file e "
        "non le aveva viste. Le letture piu' corte del troncamento verrebbero "
        "scartate, non accorciate: abbassa filter.truncLen (il report indica un "
        "valore suggerito dal profilo di lunghezza e di qualita'), oppure alza "
        "qc.max_frac_short_reads accettando la perdita, che l'errore quantifica "
        "per classe e per campione. La frazione si misura separatamente sui "
        "campioni biologici e sui controlli positivi, e basta che una delle due "
        "classi superi la soglia; i controlli negativi non contano.",
        _UMANA,
    ),
    _v(
        "E-S1-03", "S1",
        "Le letture hanno pochi valori di qualita' distinti.",
        "L'esecuzione prosegue. Qualita' raggruppate in pochi valori sono proprie "
        "di alcuni sequenziatori (NovaSeq, NextSeq): con esse la stima standard "
        "del modello d'errore puo' dare tassi non monotoni rispetto alla qualita'. "
        "Valuta err.error_function: loess_monotono, e controlla i grafici in "
        "04_error_models. Con un solo valore di qualita' nessuna funzione puo' "
        "stimare il modello, e S3 si ferma con E-S3-04.",
        _DEGRADA,
    ),
    _v(
        "E-S1-04", "S1",
        "Un file di letture, su tutte le sue letture, non ha un layout single-end.",
        "G07, in S0, guarda le sole prime qc.head_reads letture di ogni file: un "
        "file con tutte le prime letture seguite da tutte le seconde gli sfugge. "
        "Qui i segni di coppia sono contati su tutte le letture "
        "(02_qc_profiles/coppie.tsv). Se un file contiene le due letture di ogni "
        "coppia (lo dicono le intestazioni, con /1 e /2 o con lo stesso nome due "
        "volte), o le sole seconde letture (intestazioni tutte con /2), estrai le "
        "sole letture forward in un file per campione: questa pipeline tratta "
        "solo letture singole.",
        _UMANA,
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
        "Il lotto viene riletto senza modifiche: un errore di lettura "
        "transitorio non richiede altro. Se l'errore persiste oltre "
        "retry.max_attempts, il file e' probabilmente corrotto: verificane "
        "l'integrita' alla sorgente.",
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
    _v(
        "E-S3-02", "S3",
        "Alcuni campioni non hanno la corsa di sequenziamento, con "
        "err.batch_column attivo.",
        "Un campione senza corsa non ha un modello d'errore da cui farsi "
        "correggere. Completa io.batch_table per i campioni indicati, oppure "
        "imposta err.batch_column a null per stimare un solo modello su tutti "
        "i campioni.",
        _UMANA,
    ),
    _v(
        "E-S3-03", "S3",
        "Una corsa di sequenziamento non ha letture filtrate.",
        "Nessun campione della corsa indicata ha letture dopo il filtro, quindi "
        "il suo modello d'errore non si puo' stimare. Guarda quante letture il "
        "filtro ha tolto a quei campioni in 03_filtered/letture_filtrate.tsv: se "
        "la corsa e' fatta di soli controlli negativi vuoti, togli i campioni o "
        "imposta err.batch_column a null per stimare un solo modello.",
        _UMANA,
    ),
    _v(
        "E-S3-04", "S3",
        "Il modello d'errore non e' stimabile dalle qualita' delle letture.",
        "La stima lega il tasso d'errore alla qualita' delle basi: con un solo "
        "valore di qualita' in tutte le letture (file a qualita' normalizzata, "
        "come quelli di alcuni archivi pubblici) non c'e' una curva da adattare, "
        "con nessuna delle funzioni di err.error_function. I valori presenti "
        "sono in 02_qc_profiles/valori_qualita.tsv. Servono le letture con le "
        "qualita' originali del sequenziatore.",
        _UMANA,
    ),
    _v(
        "E-S3-05", "S3",
        "La stima del modello d'errore e' fallita.",
        "La funzione di err.error_function non ha potuto adattare la curva dei "
        "tassi d'errore, e le letture hanno piu' di un valore di qualita': non e' "
        "la condizione di E-S3-04. Il messaggio originale accompagna l'errore. "
        "Le cause note sono poche letture filtrate nella corsa (guarda "
        "03_filtered/letture_filtrate.tsv) e qualita' raggruppate in due o tre "
        "valori (02_qc_profiles/valori_qualita.tsv): nel secondo caso prova "
        "l'altra funzione di err.error_function. Se nessuna delle due spiega "
        "l'arresto, segnalalo con il log in 99_logs.",
        _UMANA,
    ),
    # ----------------------------------------------------------------- S4 ---
    _v(
        "E-S4-02", "S4",
        "Memoria esaurita durante l'inferenza delle varianti.",
        "L'inferenza viene ritentata con run.batch_size ridotto: il lotto "
        "cambia la memoria richiesta, non le varianti inferite (con dada.pool "
        "vero il lotto non ha effetto, e non si ritenta). Se l'errore "
        "persiste, aumenta la memoria disponibile al container o riduci "
        "run.threads: ogni thread ne consuma una quota.",
        _RETRY,
    ),
    # ----------------------------------------------------------------- S5 ---
    _v(
        "E-S5-01", "S5",
        "Le varianti distinte sono troppe per la tabella delle sequenze.",
        "Le varianti distinte superano qc.max_asv_count, oppure la memoria si e' "
        "esaurita costruendo la tabella. Un nuovo tentativo non cambierebbe "
        "nulla: la tabella e' una matrice densa campioni x varianti, la cui "
        "dimensione non dipende da run.batch_size, e la fase lo dichiara. La "
        "decisione e' scientifica: valuta un filtro piu' severo a monte, "
        "anziche' alzare la soglia; solo se il numero di varianti e' plausibile, "
        "aumenta la memoria disponibile al container.",
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
    # ----------------------------------------------------------------- S7 ---
    _v(
        "E-S7-01", "S7",
        "Nessuna variante resta nella tabella dopo il filtro di lunghezza.",
        "Guarda lunghezze.tsv in 07_chimera: se le varianti hanno lunghezze fuori "
        "da asv.len_min-asv.len_max, verifica filter.truncLen, filter.trimLeft e "
        "asv.len_tol; se la tabella senza chimere di S6 era gia' vuota, la causa "
        "e' a monte, nel filtro, nell'inferenza o nella rimozione delle chimere. "
        "Un nuovo tentativo darebbe la stessa tabella vuota.",
        _UMANA,
    ),
    # ----------------------------------------------------------------- S8 ---
    _v(
        "E-S8-02", "S8",
        "La copertura tassonomica e' insufficiente.",
        "In una classe controllata la frazione di varianti con il phylum assegnato "
        "e' sotto qc.min_frac_phylum: verifica che tax.ref_fasta copra la regione "
        "amplificata e che tax.min_boot non sia troppo alto; le misure per classe "
        "sono in 08_taxonomy/riepilogo.json.",
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
    _v(
        "E-S9-02", "S9",
        "Le varianti finali sono troppo poche per costruire un albero.",
        "L'esecuzione si ferma prima di avviare il calcolo: un albero radicato ha "
        "una topologia da stimare solo con almeno quattro sequenze. Il codice si "
        "presenta solo se phylo.enabled e' true. Controlla in "
        "12_final/intermedi/filtri_riepilogo.json e in varianti_rimosse.tsv, "
        "nella stessa cartella, quali filtri "
        "hanno tolto le varianti: se cosi' poche varianti finali sono un esito "
        "inatteso, rivedi le soglie di S13 (filt, prev) e la decontaminazione; se "
        "sono attese, rinuncia all'albero con phylo.enabled: false.",
        _UMANA,
    ),
    # ---------------------------------------------------------------- S10 ---
    _v(
        "E-S10-01", "S10",
        "Gli slot dell'oggetto integrato non sono coerenti fra loro.",
        "Campioni o varianti non concordano fra gli slot, oppure la tabella dei "
        "conteggi non ha l'orientamento di out.taxa_are_rows: righe e colonne "
        "si riconoscono dagli identificativi (varianti e accession), non dalle "
        "dimensioni. Il dettaglio dice quale slot non torna. Conserva "
        "10_phyloseq e verifica con il manifesto l'integrita' degli artefatti "
        "delle fasi precedenti: l'incoerenza nasce quasi sempre a monte.",
        _UMANA,
    ),
    _v(
        "E-S10-02", "S10",
        "Una colonna dei metadati richiesta per l'oggetto integrato non si puo' "
        "portare nell'oggetto senza ambiguita'.",
        "Il dettaglio indica la colonna di out.study_columns o out.batch_columns "
        "e il motivo: assente nella tabella, ripetuta nella sua intestazione, "
        "chiesta al file di arricchimento senza io.batch_table, oppure con un "
        "nome nell'oggetto gia' usato da un'altra colonna. Correggi l'elenco "
        "nella configurazione: indovinare quale colonna si intendesse "
        "attaccherebbe ai campioni un dato diverso da quello richiesto.",
        _UMANA,
    ),
    # ---------------------------------------------------------------- S11 ---
    _v(
        "E-S11-02", "S11",
        "Almeno una piastra, o un campione senza piastra, non usa una soglia "
        "di profondita' derivata da una curva propria.",
        "Si prosegue con la soglia della curva aggregata, se e' valida, altrimenti "
        "con la mediana delle soglie delle piastre con una curva propria valida, "
        "sempre sulle letture senza chimere: il dettaglio e 11_controls/soglia.json "
        "riportano per ognuno l'origine della soglia e il motivo. Per una curva "
        "propria servono almeno ctrl.min_positives controlli utilizzabili nella "
        "piastra, una bonta' non inferiore a katharoseq.min_r2 e una soglia dentro "
        "le profondita' osservate e determinata dai dati: le misure di ogni curva "
        "sono in curve.tsv.",
        _DEGRADA,
    ),
    _v(
        "E-S11-03", "S11",
        "I controlli positivi non sono conformi.",
        "Verifica ctrl.positive_values e katharoseq.target_taxon, poi i "
        "controlli non conformi in 11_controls/positivi.tsv: un controllo "
        "positivo anomalo per la sua concentrazione mette in dubbio l'intero "
        "lotto, quindi la decisione di proseguire non puo' essere automatica.",
        _UMANA,
    ),
    _v(
        "E-S11-04", "S11",
        "I controlli positivi conformi alla loro concentrazione sono meno di "
        "ctrl.min_positive_pass_frac.",
        "L'esecuzione prosegue perche' ctrl.positive_gate e' falso; i controlli "
        "non conformi e il motivo sono in 11_controls/positivi.tsv. Esaminali "
        "prima di usare i risultati: se il problema e' sistematico, imposta "
        "ctrl.positive_gate a true perche' fermi l'esecuzione (E-S11-03).",
        _DEGRADA,
    ),
    # ---------------------------------------------------------------- S12 ---
    _v(
        "E-S11-05", "S11",
        "Nessuna curva di calibrazione valida: nessuna soglia di profondita' "
        "dai controlli positivi.",
        "Si prosegue senza il filtro per profondita': in S13 resta il solo "
        "qc.min_reads_final, sulle letture dell'oggetto finale, e "
        "11_controls/soglia.json ne registra il motivo. Succede senza controlli "
        "positivi, senza la colonna delle cellule seminate "
        "(katharoseq.cell_count_column), senza una piastra con almeno "
        "ctrl.min_positives controlli utilizzabili, o se nessuna curva e' valida "
        "(curve.tsv). Va dichiarato con i risultati: i campioni non sono stati "
        "selezionati per profondita'.",
        _DEGRADA,
    ),
    _v(
        "E-S12-02", "S12",
        "I contaminanti individuati superano la quota attesa.",
        "La modalita' dichiarata in decontam.mode rimuoverebbe dai campioni "
        "biologici una frazione delle letture oltre qc.max_frac_contaminant: le "
        "misure di entrambe le modalita' sono in 11_controls/decontam_riepilogo.json. Esamina i controlli negativi e "
        "decontam.batch_column prima di accettare il risultato: una quota anomala "
        "indica piu' spesso un raggruppamento sbagliato che una contaminazione reale.",
        _UMANA,
    ),
    _v(
        "E-S12-03", "S12",
        "I controlli negativi non bastano: nessuna decontaminazione.",
        "Si prosegue senza togliere alcuna variante: i controlli negativi con "
        "letture sono meno di decontam.min_blanks, e la prevalenza nei negativi "
        "non e' stimabile. I contaminanti da reagente restano nell'oggetto "
        "finale: va dichiarato con i risultati. Se i controlli esistono ma non "
        "sono stati riconosciuti, correggi ctrl.blank_values.",
        _DEGRADA,
    ),
    _v(
        "E-S12-04", "S12",
        "Nessun campione biologico ha letture da decontaminare.",
        "Nell'oggetto integrato di S10 i campioni biologici non hanno letture: "
        "non c'e' nulla da confrontare con i controlli negativi, in nessuna "
        "modalita' di decontam.mode. Guarda in 99_logs il tracciamento delle "
        "letture per campione per capire in quale fase i biologici le hanno "
        "perse (filtro, inferenza, chimere, lunghezza delle varianti), e "
        "verifica che ctrl.biological_values riconosca i campioni giusti.",
        _UMANA,
    ),
    _v(
        "E-S12-05", "S12",
        "La decontaminazione toglierebbe tutte le varianti.",
        "Ogni variante risulta piu' prevalente nei controlli negativi che nei "
        "campioni, e qc.max_frac_contaminant e' abbastanza alto da ammetterlo: "
        "l'oggetto decontaminato sarebbe vuoto. Le misure sono nel messaggio. "
        "Verifica ctrl.blank_values e ctrl.biological_values (le classi "
        "potrebbero essere scambiate), decontam.threshold e decontam.mode; un "
        "valore di qc.max_frac_contaminant pari a 1 toglie il controllo che "
        "avrebbe fermato prima, con E-S12-02.",
        _UMANA,
    ),
    # ---------------------------------------------------------------- S13 ---
    _v(
        "E-S13-01", "S13",
        "Il filtro di prevalenza e' stato avviato prima che la decontaminazione "
        "fosse conclusa.",
        "La prevalenza va calcolata su dati gia' decontaminati: un contaminante "
        "diffuso supererebbe il filtro proprio perche' compare ovunque. Completa "
        "S12 prima di S13; se l'errore compare durante una ripresa, e' un "
        "difetto dell'orchestrazione: segnalalo insieme al log in 99_logs.",
        _UMANA,
    ),
    _v(
        "E-S13-02", "S13",
        "Il filtro di prevalenza lascia alcuni campioni senza varianti.",
        "Si prosegue: i campioni contenevano solo varianti che il filtro di "
        "prevalenza toglie, ed escono dall'oggetto finale con il motivo in "
        "12_final/intermedi/esclusioni.tsv. Se sono molti, prev.min_fraction o "
        "prev.min_count sono troppo severi per questo dataset: abbassali.",
        _DEGRADA,
    ),
    _v(
        "E-S13-03", "S13",
        "Il filtro tassonomico lascia alcuni campioni senza letture.",
        "Si prosegue: i campioni contenevano solo letture senza phylum o dei taxa "
        "di filt.exclude_taxa (organelli, eucarioti), cioe' nessun segnale "
        "batterico, ed escono dall'oggetto finale con il motivo in "
        "12_final/intermedi/esclusioni.tsv. Se sono molti, verifica il "
        "campionamento e filt.exclude_taxa.",
        _DEGRADA,
    ),
    _v(
        "E-S13-04", "S13",
        "Nessun campione biologico resta nell'oggetto finale.",
        "I filtri per profondita', tassonomico, di prevalenza e delle letture "
        "finali hanno escluso tutti i campioni: i motivi, campione per campione, "
        "sono in 12_final/intermedi/esclusioni.tsv. Rivedi la soglia di "
        "profondita' (11_controls/soglia.json), prev.min_fraction e "
        "qc.min_reads_final.",
        _UMANA,
    ),
    _v(
        "E-S13-05", "S13",
        "La tassonomia non ha il rango Phylum: il filtro sulle varianti senza "
        "phylum non e' stato applicato.",
        "Si prosegue con il solo filtro sui taxa di filt.exclude_taxa. Con "
        "filt.remove_na_phylum vero ci si attende una colonna Phylum nella "
        "tassonomia: verifica i ranghi del riferimento in 08_taxonomy/riepilogo.json.",
        _DEGRADA,
    ),
    # ---------------------------------------------------------------- S14 ---
    _v(
        "E-S14-01", "S14",
        "La validazione finale dell'oggetto non e' superata.",
        "Non usare ps_final.rds: conserva 12_final e il log in 99_logs, che "
        "insieme indicano quale controllo non e' stato superato (componenti non "
        "allineati, un campione o una variante vuoti, export che non ricostruiscono "
        "l'oggetto, checksum, frazione delle letture trattenute sotto "
        "qc.min_frac_reads_retained).",
        _UMANA,
    ),
    # -------------------------------------------------------------- GRAFO ---
    _v(
        "E-GRAFO-01", "GRAFO",
        "Una fase e' stata avviata prima che le fasi da cui dipende fossero "
        "concluse.",
        "Le fasi mancanti sono indicate nell'errore: eseguile prima, oppure "
        "riprendi l'esecuzione, che le esegue nell'ordine del grafo. Se l'errore "
        "compare durante una ripresa, e' un difetto dell'orchestrazione e non "
        "dei dati: segnalalo insieme al log in 99_logs.",
        _UMANA,
    ),
    # --------------------------------------------------------------- PROV ---
    _v(
        "E-PROV-01", "PROV",
        "La regola rigorosa sulla provenienza e' attiva e nessun commit leggibile "
        "identifica il codice eseguito.",
        "Con run.strict_provenance: true l'esecuzione parte solo da un clone del "
        "repository: il commit e' cio' che identifica il codice eseguito. Lancia la "
        "pipeline dalla radice di un clone, con il repository montato nel container e "
        "il comando del README (PYTHONPATH e AMPLICON16S_R_DIR che puntano al clone), "
        "con lo stesso utente proprietario dei file, altrimenti git rifiuta il "
        "repository. Senza un clone, disattiva la regola: i risultati non saranno "
        "certificati rispetto al codice.",
        _UMANA,
    ),
    _v(
        "E-PROV-02", "PROV",
        "La regola rigorosa sulla provenienza e' attiva e il codice ha modifiche non committate.",
        "Con run.strict_provenance: true il codice eseguito deve coincidere con un "
        "commit: le modifiche non committate in src/ o in R/ non sarebbero "
        "identificate da nulla. Committa le modifiche o ripristina i file (git "
        "status, git stash), poi riprendi; durante lo sviluppo, disattiva la regola.",
        _UMANA,
    ),
    _v(
        "E-PROV-03", "PROV",
        "La regola rigorosa sulla provenienza e' attiva e l'ambiente R non corrisponde al file di blocco.",
        "Con run.strict_provenance: true la versione di R e quelle dei pacchetti "
        "installati devono essere quelle di run.lockfile, compresa la correzione di "
        "dada2. Il dettaglio elenca le discordanze. Esegui nell'immagine indicata "
        "dal README, scaricata con il suo digest, e non con l'R della macchina; se "
        "il file di blocco non si trova, lancia la pipeline dalla radice del "
        "repository.",
        _UMANA,
    ),
    # ------------------------------------------------------------------ R ---
    _v(
        "E-R-01", "R",
        "L'interprete R o gli script R della pipeline non sono disponibili.",
        "Installa R oppure indica l'interprete con la variabile d'ambiente "
        "AMPLICON16S_RSCRIPT; se gli script non vengono trovati, indica la "
        "cartella R/ del repository con AMPLICON16S_R_DIR. Nel container "
        "entrambi sono gia' presenti: se l'errore compare li', l'immagine non "
        "e' quella prevista.",
        _UMANA,
    ),
    _v(
        "E-R-02", "R",
        "Il processo R e' terminato senza dichiarare un esito valido.",
        "Il processo si e' interrotto prima di poter scrivere l'esito: per un "
        "segnale, per un guasto dell'interprete, o per un errore avvenuto "
        "prima che lo script caricasse le funzioni condivise. Il codice di "
        "uscita accompagna l'errore e l'uscita di errore del processo e' nel "
        "log strutturato in 99_logs: parti da li', perche' rieseguire senza "
        "conoscerne la causa ripeterebbe il guasto.",
        _UMANA,
    ),
    _v(
        "E-R-03", "R",
        "Lo script R si e' interrotto con un errore privo di un codice del "
        "catalogo.",
        "E' un difetto dello script, non dei dati: ogni condizione prevista "
        "va dichiarata con un codice del catalogo. Il messaggio originale di R "
        "accompagna l'errore; segnalalo insieme al log in 99_logs.",
        _UMANA,
    ),
    _v(
        "E-R-04", "R",
        "Memoria esaurita in un processo R di una fase che non prevede il "
        "retry.",
        "Aumenta la memoria disponibile al container oppure riduci "
        "run.threads, poi riprendi l'esecuzione: per questa fase non esiste "
        "un'azione correttiva automatica che lasci intatte le sue assunzioni.",
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
