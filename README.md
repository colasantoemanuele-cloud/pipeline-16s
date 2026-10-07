# amplicon16s

Pipeline per l'analisi di dati 16S rRNA da sequenziamento Illumina single-end.

## Contesto

Il progetto è commissionato dall'Agenzia Spaziale Italiana. Il dataset di riferimento è
OSD-734 (NASA GeneLab/OSDR), 960 campioni prelevati dalla Stazione Spaziale
Internazionale.

Il dataset è rappresentativo, non esclusivo: la pipeline è progettata per essere
generalizzabile a qualunque dataset 16S single-end con caratteristiche simili. Ogni
valore che dipende dal dataset è un parametro di configurazione, mai un valore scritto
nel codice. Questo vale anche per i percorsi dei dati, che non risiedono nel repository.

## Impostazione tecnica

Scelta di progetto: due linguaggi con ruoli distinti.

- **Python ≥ 3.11** orchestra l'esecuzione.
- **R ≥ 4.1** esegue i calcoli scientifici (dada2, phyloseq, DECIPHER, phangorn, decontam) come
  processi separati, non come libreria caricata nel processo Python.

Lo scambio fra i due passa per il filesystem: parametri in ingresso, artefatti su disco
in uscita. La ragione è l'isolamento dei guasti: un errore in una routine di calcolo
termina il processo figlio senza abbattere l'orchestratore.

## Struttura del repository

```
src/amplicon16s/      Pacchetto Python della pipeline: orchestrazione e riga di comando
src/amplicon16s_eco/  Pacchetto Python separato per le analisi ecologiche a valle
R/                    Script R dei calcoli scientifici; R/lib/ per le funzioni condivise
config/               File di configurazione (config.example.yaml)
tests/                Suite di test
dati/osd734/          Dati del dataset di riferimento: metadati, elenchi con i checksum, script di scarico
dati/osd276/          Dati del secondo dataset, su cui si verifica la generalità: metadati, script di scarico, configurazione
dati/riferimento/     Riferimento tassonomico comune ai dataset: script di scarico, fonti e impronte
docs/                 Registro delle decisioni di metodo (decision_log.md) e analisi di sensibilità (sensibilita.md, sensibilita/)
scripts/              Script di utilità
container/            Definizione dell'ambiente riproducibile (Dockerfile e script R)
```

## Requisiti e installazione

Richiede Python ≥ 3.11. Installare sempre in un ambiente virtuale dedicato:

```bash
python3.11 -m venv .venv
source .venv/bin/activate
pip install -e ".[dev]"
```

È possibile installare le dipendenze sia tramite `pip install -e ".[dev]"` sia tramite
`pip install -r requirements-dev.txt` (seguito da `pip install --no-deps -e .` per
registrare il pacchetto in modalità modificabile). L'installazione rende disponibile il
comando `amplicon16s`. La suite di test si esegue
con `pytest`, su quattro processi con pytest-xdist (`addopts` in `pyproject.toml`), in
locale come nella CI; `pytest -n 0` la esegue in sequenza. Le esecuzioni condivise
del sottoinsieme di prova si calcolano una volta sola per tutti i processi. I test sui
dati reali (marcatore `dati_reali`, attivi con `AMPLICON16S_CONFIG_DATI_REALI`)
girano tutti sullo stesso processo, uno dopo l'altro, mentre gli altri restano
distribuiti sui quattro: `--dist loadgroup` e il gruppo `dati_reali`, assegnato in
`tests/conftest.py`. Sul dataset completo ciascuno lancia catene con decine di
processi R, e in parallelo saturavano la memoria. I calcoli scientifici girano nell'immagine descritta qui sotto, che porta
con sé il proprio R. I test del ponte verso R lanciano invece script R veri, e chiedono
solo `Rscript` nel PATH (o indicato con `AMPLICON16S_RSCRIPT`) e il pacchetto `jsonlite`;
i test delle fasi di calcolo chiedono anche i pacchetti Bioconductor, e girano
nell'immagine. Dove mancano, quei test si saltano.

### Ambiente containerizzato

L'immagine contiene sia Python sia R con tutte le librerie della pipeline, ed è il modo
previsto per eseguire la pipeline in modo riproducibile. La via principale è scaricare
l'immagine pubblicata, identificata per digest:

```bash
docker pull ghcr.io/colasantoemanuele-cloud/amplicon16s@sha256:930dda008581b24cc3a3fea718c0fed3690b011958bf1716c643f92056970154
docker run --rm ghcr.io/colasantoemanuele-cloud/amplicon16s@sha256:930dda008581b24cc3a3fea718c0fed3690b011958bf1716c643f92056970154
```

Il digest identifica i byte dell'immagine: il nome con cui è pubblicata potrà cambiare,
il digest no. È l'immagine a cui si riferiscono i checksum attesi dei risultati in
`dati/osd734/`, e quella dichiarata in `run.container` da
`dati/osd734/config_osd734.yaml`. Lo scarico è di 1,86 GB.

In alternativa l'immagine si costruisce dal Dockerfile:

```bash
docker build -f container/Dockerfile -t <immagine> .
```

Due costruzioni non sono identiche byte per byte (date dei file, compressione,
metadati): un'immagine costruita in locale ha un altro identificativo e non è l'immagine
pubblicata, anche se contiene gli stessi pacchetti alle stesse versioni, verificati
contro `renv.lock`. Nelle prove fatte, tre immagini costruite in momenti diversi hanno
dato gli stessi 1.031 artefatti sul dataset di riferimento.

Nei comandi dei test `<immagine>` è il riferimento all'immagine in uso: quella corrente
del progetto è indicata in `test.txt` (sezione 1.3). Il comando completo per eseguire
la pipeline nel container, con il repository montato, è in `dati/osd734/README.md`.

L'immagine contiene anche i pacchetti R delle analisi ecologiche a valle (DESeq2, vegan,
randomForest), registrati in `renv.lock` come quelli della pipeline. Il codice di quelle
analisi non è nell'immagine: si eseguirà dal repository montato.

Le versioni sono bloccate su entrambi i fronti: le dipendenze Python in
`pyproject.toml` (e nei file `requirements.txt` / `requirements-dev.txt`), quelle R in
`renv.lock`. `renv.lock` è generato dall'immagine già
costruita, quindi registra le versioni effettivamente ottenute e non quelle attese; la
build lo rilegge e fallisce se anche una sola versione, o una correzione dichiarata, non
coincide. Dopo ogni modifica ai pacchetti R va rigenerato: si costruisce l'immagine da
una copia del repository senza `renv.lock`, che salta la verifica, si rigenera il file
da quell'immagine e si ricostruisce con la verifica attiva:

```bash
scripts/genera_renv_lock.sh <immagine>
```

Le immagini di partenza sono ancorate per digest e non per tag, perché un tag può essere
riassegnato a un'immagine diversa mentre un digest no. Per la stessa ragione le
dipendenze R installate da CRAN vengono da un'istantanea datata del repository P3M
(`CRAN_SNAPSHOT` nel Dockerfile) e non dall'indirizzo `latest` dell'immagine di partenza,
che cambia nel tempo: una costruzione da zero, come quella della CI, otterrebbe
altrimenti versioni diverse da `renv.lock`.

**L'immagine contiene dada2 con una correzione.** In dada2 1.36.0, `assignTaxonomy`
sceglie fra generi a pari probabilità con `std::random_device` (`src/taxonomy.cpp`,
`get_best_genus`), un seme preso dal sistema che `set.seed()` non controlla: due
esecuzioni con lo stesso seme non danno la stessa tassonomia. La correzione
(`container/dada2/dada2-1.36.0-pareggi.patch`) inizializza il generatore dei pareggi
con un seme che dipende dal seme di R, dall'indice della sequenza e dal caso
(classificazione diretta, complementare inversa, replica di bootstrap): ogni pareggio
ha un seme fisso qualunque thread lo esegua. La scelta resta uniforme fra i generi a
pari probabilità, e fuori dai pareggi il risultato non cambia. Il Dockerfile scarica i
sorgenti ufficiali di dada2 1.36.0, ne verifica lo SHA-256, applica la correzione e li
compila come versione 1.36.0.1; il DESCRIPTION del pacchetto e `renv.lock` registrano la
versione di partenza, la correzione e il suo SHA-256, e la verifica di `renv.lock`
accetta la versione modificata solo se questi coincidono.
Misure sul dataset di riferimento (12.045 varianti, SILVA 138, stesso seme). Con la
versione ufficiale due esecuzioni a 12 thread differiscono in 3 varianti per la
tassonomia e in 99 per il bootstrap, e in dieci esecuzioni 160 varianti non restano
identiche. Il bootstrap cambia al più di 8 punti, e le assegnazioni cambiano a cavallo
di `tax.min_boot`. Nessuna di quelle 160 varianti è priva di pareggi: 820 varianti ne
incontrano almeno uno (120 nella classificazione principale, 1.662 repliche di bootstrap
su 1.204.500).
Con la correzione, tre esecuzioni a 12, 12 e 4 thread danno gli stessi byte. Confrontata
con dieci esecuzioni ufficiali, differisce solo sulle 160 varianti che la versione
ufficiale stessa fa variare, tutte con pareggi, e su nessuna senza pareggi.

### Dati del dataset di riferimento

I dati con cui la pipeline si esegue su OSD-734 si recuperano da `dati/osd734/` (vedi il
suo README): nel repository ci sono le tabelle ISA usate dalla pipeline, l'elenco delle
960 corse con l'MD5 e la dimensione dichiarati da ENA, la provenienza e la licenza di
ogni file (`FONTI.tsv`) e una configurazione con i percorsi relativi alla radice del
repository;
le letture, il riferimento tassonomico e il file del lotto si ottengono con tre script,
dalla radice del repository e con la sola libreria standard di Python. Il riferimento
tassonomico è comune ai dataset e sta in `dati/riferimento/`, con le sue fonti e le
sue impronte:

```bash
python3 dati/osd734/scarica_letture.py            # 960 file FASTQ da ENA, 2,52 GB
python3 dati/riferimento/scarica_riferimento.py   # SILVA 138 per dada2 da Zenodo
python3 dati/osd734/ricostruisci_lotto.py         # piastre, pozzetti e corse dalla fonte degli autori
```

La pipeline si esegue poi nel container, con il repository montato: il README della
cartella riporta i comandi completi, i valori da adattare e i requisiti misurati di
spazio, durata e memoria. A esecuzione
conclusa i risultati si verificano con un comando contro i checksum attesi pubblicati
nella stessa cartella:

```bash
python3 dati/osd734/confronta_risultati.py    # 1.031 artefatti, fase per fase
```

Gli script verificano ogni file contro il checksum della sua fonte, riprendono uno scarico
interrotto e non scaricano ciò che è già presente e integro. Il file del lotto non è nel
repository perché il repository degli autori da cui deriva non dichiara una licenza: lo
script lo ricostruisce da quella fonte, fissata a un commit. Senza quel file la pipeline
procede con un solo modello d'errore e senza soglie di profondità per piastra. I file
scaricati sono esclusi da git.

### Identificare l'immagine per digest

L'immagine pubblicata si cita con il digest di registro riportato sopra. Un tag locale
identifica l'immagine solo su quella macchina e può
essere riassegnato. Per riferirsi senza ambiguità a un'immagine costruita in locale, si
usa il suo digest:

```bash
# Digest dell'immagine costruita in locale
docker image inspect <immagine> --format '{{.Id}}'

# Digest con cui l'immagine è pubblicata in un registry, disponibile dopo il push
docker image inspect <immagine> --format '{{index .RepoDigests 0}}'
```

Il primo comando restituisce l'identificatore del contenuto locale; il secondo il
digest con cui l'immagine viene recuperata da un registry, ed è quello da citare quando
si deve riprodurre un'analisi a distanza di tempo.

## Uso della riga di comando

Ogni sottocomando richiede il file di configurazione, che non viene mai modificato; i
suoi percorsi relativi valgono rispetto alla cartella da cui il comando è lanciato:

```bash
amplicon16s validate --config config.yaml   # esegue solo S0
amplicon16s run      --config config.yaml   # esegue dall'inizio, in una cartella nuova
amplicon16s resume   --config config.yaml   # riprende dagli artefatti esistenti
amplicon16s report   --config config.yaml   # genera il report dell'esecuzione
```

`run` non sovrascrive mai un'esecuzione: ogni esecuzione deve restare ispezionabile. Se
la cartella indicata da `io.out_root` non è vuota, `run` si rifiuta e indica le due
strade: `resume` per continuare quell'esecuzione, oppure una `io.out_root` nuova per
cominciarne un'altra. Non esiste un'opzione per sovrascrivere. `validate` non riesegue
S0 se è già conclusa e valida per la configurazione data.

**La configurazione usata resta registrata.** All'avvio, superati i controlli di
coerenza della configurazione e delle risorse, e prima di qualunque fase, la
configurazione risolta viene registrata in `00_config/` con il suo digest, e una
versione registrata non viene mai sovrascritta:

- `run` scrive `resolved.yaml`;
- `resume`, se il digest è identico a quello dell'ultima versione registrata, non
  scrive nulla; se è diverso (anche solo per `run.threads` o un parametro di retry,
  che sono fuori dall'impronta dei risultati ma non dal digest), conserva le versioni
  precedenti e scrive la nuova accanto: `resolved_2.yaml`, `resolved_3.yaml` e così
  via, ciascuna con il nome della precedente e i parametri che ne differiscono.
  Tornare a una configurazione già usata registra a sua volta una nuova versione;
- `validate` segue la regola di `run` su una cartella nuova e quella di `resume` se
  S0 è già conclusa.

Il log in `99_logs` registra a ogni avvio con quale versione si esegue. Un arresto ai
controlli di avvio non registra alcuna configurazione, perché nessuna fase è partita;
l'esito dei controlli resta comunque nel log.

**Il registro degli avvii.** Gli eventi che descrivono un'esecuzione (controlli di
avvio, configurazione in uso, fasi eseguite, conclusione) sono scritti, oltre che nel
log, in `99_logs/avvii.jsonl`, nella stessa forma. Il log ruota (16 MB, cinque
rotazioni) e a lungo andare perde gli avvii più vecchi; il registro è solo in aggiunta
e non ruota, poche righe per avvio, e conserva la storia intera di un'esecuzione
ripresa molte volte.

**Il report di esecuzione.** `report` scrive in `report/`, sotto `io.out_root`, il
documento `report.html` e, in `report/tabelle/`, ogni sua tabella come file TSV. Lo
ricava da ciò che l'esecuzione ha lasciato su disco (manifesti di fase, artefatti,
configurazioni registrate, registro degli avvii): non riesegue alcun calcolo, non usa R e
non legge i dati grezzi, quindi si genera anche su un'esecuzione conclusa in
precedenza e su una macchina che ha la sola cartella di output. Della configurazione
indicata usa soltanto `io.out_root`. Il documento è autoconsistente (stili e grafici
incorporati, nessuna risorsa esterna), non riporta la data in cui è generato e non
scrive nel log: due generazioni dalla stessa cartella danno gli stessi byte, e la
valutazione dello stato dà lo stesso esito prima e dopo. Senza un'esecuzione in
`io.out_root` il comando esce con il codice 3.

Codici di uscita, per chi lancia la pipeline da uno script o da uno scheduler:

| Codice | Significato |
|---|---|
| 0 | successo |
| 1 | errore imprevisto, cioè un difetto del programma; la traccia è nel log in `99_logs` |
| 2 | riga di comando non valida (argomenti mancanti o sconosciuti) |
| 3 | errore di configurazione: il file non è valido, G15 lo respinge, oppure `run` trova la cartella di output già usata; nessuna fase è partita |
| 4 | arresto con punto di ripresa dichiarato, stampato e scritto in `99_logs/punto_di_ripresa.json` e `.txt` |
| 5 | riservato: una fase prevista dal grafo non esiste come codice. Con le quindici fasi realizzate non si presenta in alcuna esecuzione; il codice resta come guardia per un grafo esteso prima della fase che lo realizza |

## Stato dell'implementazione

Sono realizzati:

- la struttura del repository e il packaging Python, installabile in modalità sviluppo;
- l'ambiente di esecuzione containerizzato, che contiene Python 3.11 e R 4.5 con i
  pacchetti Bioconductor della pipeline, con le versioni bloccate e verificate a ogni
  costruzione dell'immagine;
- **l'immagine pubblicata** nel registro dei container, identificata per digest
  (`ghcr.io/colasantoemanuele-cloud/amplicon16s@sha256:930dda008581b24cc3a3fea718c0fed3690b011958bf1716c643f92056970154`),
  1,86 GB da scaricare. Prima della pubblicazione ne è stato verificato il contenuto:
  nessun percorso personale, nessuna credenziale, nessun dato oltre alle 28 letture del
  sottoinsieme di prova, nulla di ciò che `.dockerignore` esclude. Rispetto alle
  immagini precedenti contiene i pacchetti R delle analisi ecologiche (DESeq2, vegan,
  randomForest): in `renv.lock` sono tre voci in più (DESeq2, locfit, randomForest;
  vegan c'era già, come dipendenza di phyloseq) e nessuna versione cambiata. Con questa
  immagine la catena completa sul dataset di riferimento dà i 1.031 artefatti dei
  checksum pubblicati, e sul sottoinsieme di prova gli stessi 99 artefatti
  dell'immagine locale. `dati/osd734/config_osd734.yaml` la dichiara in
  `run.container`. La suite di test passa anche su una macchina con un solo
  processore: le esecuzioni di prova chiedono due thread dove ci sono, uno altrimenti
  (`tests/sottoinsieme.py`), perché G14 respinge, correttamente, più thread dei
  processori utilizzabili;
- la validazione della configurazione: lo schema copre tutti i parametri della
  pipeline e ne verifica tipo e dominio, respingendo le chiavi sconosciute;
  `config/config.example.yaml` ne è un'istanza completa;
- il gate G15, che verifica la coerenza fra parametri (le combinazioni singolarmente
  valide ma insensate messe insieme) e risolve i parametri che discendono da altri.
  La configurazione effettivamente usata viene registrata in `00_config/resolved.yaml`,
  con l'elenco dei parametri impostati nel file di ingresso e
  con un digest che ne identifica la combinazione, all'avvio di ogni esecuzione e
  senza mai sovrascrivere le versioni precedenti, come descritto nella sezione sull'uso
  della riga di comando. Tutto questo avviene prima che venga
  allocato qualunque calcolo: è il gate che apre la sequenza di S0, perché un errore
  di configurazione va scoperto prima di aprire un solo file. I parametri derivati
  (`asv.len_min`, `asv.len_max`) discendono dalla sola
  configurazione; il numero minimo di campioni del filtro di prevalenza lo calcola
  S13 sui biologici che tiene;
- l'inventario dei campioni: il crosswalk fra i file di letture, la tabella di assay
  dell'amplicone e la tabella campioni di studio, con la classe di ciascun campione.
  La chiave del join è l'accession estratto dal nome del file, non il nome del
  campione, che può ripetersi fra repliche: `io.accession_regex` la estrae dal nome di
  ogni file (la corrispondenza intera, o il primo gruppo di cattura se l'espressione
  ne ha uno, per i file nominati per campione o per corsa), e la stessa chiave si
  ricava dal valore di `meta.accession_column`, che può essere il nome di un file o la
  chiave stessa. La tabella campioni di studio è facoltativa: senza `io.study_table`
  classe e variabili si leggono dalla tabella di assay; se c'è e la sua colonna del
  nome del campione ha un altro nome, lo dichiara `meta.study_sample_id_column`.
  Tutte le tabelle si leggono da un solo modulo (`metadata/tabelle.py`): UTF-8 con o
  senza BOM, nomi delle colonne e valori ripuliti da spazi e virgolette. Il join verso la tabella campioni di
  studio resta ristretto alle righe dell'assay, perché quella tabella è condivisa fra
  più assay dello stesso studio. Un file facoltativo associa a ogni campione la
  piastra di estrazione, la corsa di sequenziamento e il modulo; anche quel file si
  aggancia per accession, e quando manca piastra e corsa restano nulle e la pipeline
  procede in modalità a corsa singola. Una regola di configurazione
  (`ctrl.blank_override_column` e `ctrl.blank_override_values`) tratta come controlli
  negativi i campioni con un dato valore in una colonna della tabella campioni di
  studio, qualunque sia il materiale dichiarato, che resta il valore originale
  nell'inventario: in OSD-734 sono i 33 tamponi mai aperti ("Unopened 3DMM Swab Tube"),
  dichiarati superfici ma senza alcuna superficie campionata, e i campioni diventano
  770 biologici, 80 controlli positivi e 110 controlli negativi;
- l'attribuzione del modulo, con una precedenza dichiarata: quando il file facoltativo
  fornisce la colonna del modulo, il modulo viene da lì con il nome originale;
  altrimenti si ricade sull'estrazione del prefisso della posizione tramite
  `meta.module_regex`. In entrambi i casi vale la stessa regola: un campione la cui
  posizione non è una superficie (aria, contenitori non aperti, posizioni non
  dichiarate) resta senza modulo, perché raggrupparlo per modulo mescolerebbe l'aria
  di un locale con le sue superfici. Le posizioni da considerare tali si dichiarano in
  `meta.non_surface_positions`.

  Le due vie attribuiscono quindi un modulo agli stessi campioni, e differiscono solo
  nei nomi e nella copertura. Sul dataset di riferimento la copertura differisce per
  16 campioni: il file riconosce l'Airlock, mentre la derivazione per prefisso no,
  perché `A/L1` non ha la forma che l'espressione predefinita cattura. I moduli sono
  nove con il file e otto senza;
- **i parametri che descrivono il dataset sono obbligatori** (`OBBLIGATORI` in
  `config/defaults.py`, 25 parametri): il formato dei metadati (colonne, etichette
  delle classi, colonne del file del lotto, colonne da portare nell'oggetto),
  l'esperimento (`filter.truncLen`, primer, motivo conservato, taxon atteso nei
  controlli positivi, nome e versione del riferimento) e il modo di riconoscere file e
  campioni (`io.accession_regex`, `meta.sample_id_column`, `meta.accession_column`).
  Non hanno un valore predefinito: una configurazione che non li dichiara è respinta
  da G15 con `E-G15-10` e l'elenco completo dei mancanti, prima di aprire un solo
  file. Obbligatorio significa dichiarato, non provvisto di un valore: un parametro
  non pertinente si dichiara nullo o vuoto (le etichette dei controlli in un dataset
  che non ne ha, il motivo conservato di una regione che non ne ha uno noto, le
  colonne del lotto senza il file). G15 ne verifica la coerenza: senza
  `io.batch_table` le colonne del lotto devono essere vuote o nulle (`E-G15-11`); con
  etichette di controlli positivi, taxon atteso e colonna delle cellule non possono
  essere nulli (`E-G15-12`). `config.example.yaml` li lascia da compilare, con il
  valore di OSD-734 come esempio in commento; `dati/osd734/config_osd734.yaml` e
  `dati/osd276/config_osd276.yaml` li dichiarano tutti. Otto parametri restano con un
  predefinito tarato su OSD-734 (`FATTI_OSD734`: scelte di metodo come la soglia di
  decontaminazione o la prevalenza minima), e il report segnala quelli presi per
  difetto. `run.threads` è nullo per difetto, cioè automatico: le fasi e G14 contano
  i processori utilizzabili dal processo quando servono, e il valore ricavato dalla
  macchina non entra nella configurazione registrata né nel suo digest, che è lo
  stesso su macchine diverse; `run.container` è facoltativo. G15 respinge un
  parametro dichiarato senza quello da cui dipende (`E-G15-14`: la regola di
  riclassificazione senza la sua colonna, `meta.study_sample_id_column` senza
  tabella di studio, `meta.module_regex` o `meta.non_surface_positions` senza
  `meta.module_column`), e se manca un
  intero gruppo l'elenco dei problemi nomina ogni suo parametro senza predefinito;
- **la prova su un secondo dataset** (`dati/osd276/`): OSD-276, 15 tamponi di
  superfici della stessa stazione, di un altro laboratorio; paired-end, di cui si
  usano le sole letture forward come dato single-end; 151 basi, nessun controllo,
  nessun file del lotto, file nominati per corsa. Con i predefiniti di OSD-734 la
  pipeline respingeva i file per il formato dell'accession, diagnosticava l'etichetta
  non dichiarata dei campioni come primer nelle letture, e si fermava in S10 per una
  colonna inesistente dopo dodici minuti di calcolo. Con la configurazione che
  dichiara il dataset la validazione passa in pochi secondi, dichiarando che non ci
  sono controlli positivi né negativi (`E-S0-17`), e la catena arriva in fondo:
  S11 dichiara che non c'è alcuna curva da adattare (`E-S11-05`), S12 che non c'è
  stata decontaminazione (`E-S12-03`), e S14 consegna l'oggetto finale (15 campioni,
  1.005 varianti, 774.790 letture) senza `ps_controlli.rds`; il report si genera, con
  il troncamento suggerito dalle letture. Senza controlli positivi non si applica
  alcuna soglia di profondità: resta il solo `qc.min_reads_final`;
- **un dataset senza una classe di controlli, senza file del lotto o con una sola
  piastra** attraversa le fasi di calcolo con le assenze dichiarate nei manifesti e nel
  report, o si ferma con un codice del catalogo: lo verificano quattro catene intere
  su insiemi ricavati dal sottoinsieme di prova (`tests/test_w28_calcolo.py`);
- **la fase S0, la validazione iniziale**, con tutti e quindici i gate, eseguiti
  nell'ordine delle dipendenze: G15, G01, G02, G04, G05, G06, G03, G11, G13, G07,
  G08, G09, G10, G12, G14. Le classi dei campioni (G11) si verificano prima dei gate
  che leggono le sequenze, così un'etichetta non dichiarata è diagnosticata come tale
  e non come segnale assente. G02 e G08 verificano ogni colonna che la configurazione
  nomina, sulle tabelle di assay e di studio e sul file del lotto, comprese quelle da
  portare nell'oggetto (`out.study_columns`, `out.batch_columns`) e la colonna delle
  cellule dei controlli positivi, cercata nel file del lotto e poi nella tabella di
  studio, che S10 porta nell'oggetto anche se quegli elenchi non la nominano; due
  colonne richieste che nell'oggetto prenderebbero lo stesso nome, o quello di una
  colonna dell'inventario, sono respinte da G15 (`E-G15-13`) invece che da S10 a
  calcolo concluso; G08 respinge un file del lotto in cui un
  campione non ha riga, ne ha più d'una, o non ha la piastra o la corsa dichiarate,
  e dichiara con un avviso (`E-S0-19`) le righe che non corrispondono ad alcun
  campione. G02 e G08 respingono anche un valore con una tabulazione o un a capo
  dentro un campo fra virgolette, nelle colonne che la pipeline legge: le tabelle
  che scrive non hanno virgolette, e prima se ne accorgeva S10. Senza
  `meta.module_column` il modulo viene dal file del lotto, se
  `meta.batch_module_column` è dichiarata.
  G07 riconosce le letture inverse dal marcatore che precede l'estensione e, nelle
  letture ispezionate, i file con le due letture di ogni coppia (intestazioni con
  `/1` e `/2`, o la stessa intestazione ripetuta) e quelli di sole seconde letture,
  giudicando sulla frazione delle letture e non su un record isolato (un nome in più
  di due copie, o un record di altro tipo fra le seconde letture, si tollera entro la
  stessa frazione). Un file con tutte le prime letture seguite da tutte le seconde
  sfugge a G07 se le `qc.head_reads` letture ispezionate non arrivano al secondo
  blocco: lo ferma S1, che conta gli stessi segni su tutte le letture (`E-S1-04`; gli
  altri limiti del riconoscimento sono in `docs/decision_log.md`, sezione 3.15). I FASTQ,
  compressi o no, si riconoscono dai primi byte, in S0 come nel controllo che S2 fa
  sugli archivi prima del filtro. G09 giudica la frazione di letture più corte
  del troncamento separatamente per i campioni biologici e per i controlli positivi:
  riunite in una sola frazione, una classe poco numerosa potrebbe perdere tutte le
  sue letture senza farsi notare. G10 cerca il primer in testa alle
  letture solo con `filter.trimLeft` a zero e solo nei campioni biologici e nei
  controlli positivi (in un bianco con poche letture una frazione non ha
  significato), e il motivo conservato a partire dalla
  posizione `filter.trimLeft`: così la correzione che il gate indica (il taglio pari
  alla lunghezza del primer) lo fa superare. Primer presente (`E-S0-10`) e segnale
  assente (`E-S0-16`) hanno codici e rimedi distinti; con `qc.conserved_motif` nullo
  il gate verifica il solo primer e lo dichiara (`E-S0-18`). S0 conta i controlli e
  dichiara subito, come degradazione registrata nel suo manifesto (`E-S0-17`), un
  dataset senza controlli positivi o con meno controlli negativi di
  `decontam.min_blanks`.
  Il controllo sul motivo verifica la
  presenza della regione amplificata dichiarata, non la qualità del campione:
  sul dataset di riferimento il motivo conservato compare nei controlli
  negativi quanto nei biologici, perché i bianchi amplificano contaminanti. È la barriera
  che precede qualunque calcolo costoso, e per restare tale nessun gate legge un file
  per intero: le verifiche sulle sequenze si fermano alle prime `qc.head_reads`
  letture, mentre l'integrità completa è già garantita dai checksum. Sul dataset di
  riferimento, 960 campioni, la fase impiega una ventina di secondi. L'esecuzione si
  ferma al primo gate fallito, e i gate non eseguiti sono riportati come tali;
- l'esito di S0 in `01_input_validation/`: l'esito di ogni gate, il crosswalk,
  l'inventario e le statistiche delle letture ispezionate, ciascuno registrato nel
  manifesto con il proprio checksum. Le durate dei gate vanno nel log strutturato e non
  in `gates.json`: descrivono l'esecuzione, non il risultato, e un artefatto deve avere
  lo stesso checksum fra due esecuzioni sugli stessi ingressi;
- il catalogo degli errori (76 codici totali): ogni codice porta un messaggio che dice
  cosa fare e una categoria di gestione fra revisione umana, retry automatico, retry
  seguito da revisione, e degradazione automatica. Il retry automatico è un elenco chiuso di
  quattro codici, gli stessi dichiarati in `retry.whitelist`. Sono catalogati i codici
  di tutte le fasi, comprese quelle non ancora realizzate: il catalogo esiste, il
  codice che solleverà quegli errori no. Accanto ai codici di fase ci sono i quattro
  codici `E-R-*` del ponte verso R, tutti a revisione umana, per ciò che una fase non
  può dichiarare da sé: interprete assente, processo morto senza esito, errore R privo
  di codice, memoria esaurita in una fase che non prevede il retry. Il codice
  `E-GRAFO-01`, anch'esso a revisione umana, segnala una fase avviata prima delle sue
  dipendenze. G15 respinge con `E-G15-09` una `retry.whitelist` che contenga un codice
  che il catalogo non ammette al retry: la whitelist può restringere l'elenco del
  catalogo, non allargarlo;
- **il ponte verso R** (`src/amplicon16s/rbridge/`), l'unico punto che esegue codice R.
  Ogni script gira come processo separato, lanciato con l'`Rscript` del PATH o quello
  indicato da `AMPLICON16S_RSCRIPT`: un guasto grave di R, anche un errore di
  segmentazione, termina il figlio e non l'orchestratore. Il contratto passa per due
  file JSON nella cartella di fase: la richiesta con i parametri, scritta da Python, e
  la dichiarazione d'esito, scritta dallo script per ultima e in modo atomico. I due
  file portano nel nome la fase (`rbridge_richiesta_S6.json`, `rbridge_esito_S6.json`):
  in una cartella condivisa, come `07_chimera/` per S6 e S7, ciascuna fase conserva la
  propria traccia. Una dichiarazione che elenca due volte lo stesso artefatto non è
  valida. È la
  dichiarazione a portare il codice del catalogo, che il codice di uscita di un processo
  non può rappresentare; la sua assenza distingue un processo morto senza dichiarare
  nulla da un fallimento dichiarato. Il ponte traduce l'esito in un'eccezione della
  gerarchia della pipeline con il codice dichiarato, registra nel manifesto gli
  artefatti dichiarati e nel log strutturato le uscite standard e di errore dello
  script. Riconosce la memoria esaurita in entrambe le forme (l'allocazione fallita
  intercettata da R, anche quella del codice C dei pacchetti con `R_Calloc` e quella
  del caricamento di una libreria, o il processo ucciso dal sistema con `SIGKILL`) e
  la traduce nel codice che la fase indica; gli script caricano i pacchetti con
  `richiedi_pacchetti`, che conserva il messaggio originale invece di ridurlo a
  «pacchetto non disponibile»; per rendere il riconoscimento indipendente dalla lingua
  della macchina, impone a R i messaggi in inglese (`LANGUAGE=en`), e impone
  `LC_ALL=C.UTF-8`, così che ordinamenti e codifica non dipendano dalla macchina;
  gli script ordinano le stringhe con `method = "radix"`. Può limitare la memoria virtuale
  del processo figlio, riducendo allora a uno i thread dell'algebra lineare: con
  OpenBLAS multithread, come nell'immagine, R sotto quel limite resta bloccato in
  uscita. Accetta un tempo massimo, allo scadere del quale uccide il processo e il suo
  gruppo: un figlio bloccato non blocca l'orchestratore;
- le funzioni R condivise in `R/lib/`: `io_json.R` per leggere e scrivere il JSON in
  modo atomico e per scrivere ogni file `.rds` con `salva_rds`, `errors.R` con il punto d'ingresso degli script di fase e la
  dichiarazione degli errori con un codice del catalogo, `letture.R` per registrare
  per ogni passo quante letture restano a ciascun campione, `risorse.R` per riportare
  nel log la memoria di picco del processo. Le usano gli script delle fasi e gli
  script doppioni dei test, in `tests/r_doppioni/`;
- la gestione degli artefatti: l'albero delle cartelle di output sotto
  `io.out_root` (quattordici, più `12_final/intermedi`) e, in ciascuna, un manifesto che registra ogni file scritto con il suo
  checksum; gli artefatti scritti dai processi R vi si registrano allo stesso modo di
  quelli scritti da Python. Ogni scrittura, da Python come da R, passa per un file
  temporaneo nella stessa cartella (`.scrittura-*`) poi rinominato: un'interruzione
  a metà non lascia mai un file troncato con il nome definitivo. Il completamento di una fase non si legge da quel
  manifesto, che è per cartella, ma da quello proprio della fase, descritto qui sotto;
- **la classe base delle fasi, il grafo e lo stato di un'esecuzione**
  (`steps/base.py`, `runner/graph.py`, `runner/project.py`). Ogni fase eredita da
  `PipelineStep` lo stesso scheletro: verifica dei prerequisiti, calcolo, validazione
  degli artefatti, registrazione dell'esito; S0 è realizzata su questa base, con il
  comportamento di prima. Il grafo dichiara le quindici fasi da S0 a S14 in ordine,
  con la cartella in cui ciascuna scrive e le fasi di cui consuma gli artefatti, e
  solo quelle: S0 per ogni fase che legge l'inventario dei campioni, S2 per S5, S6 e
  S7, che ne leggono le letture filtrate; S12 dipende dalla sola S10, così cambiare
  i parametri di S11 non rifà la decontaminazione. L'ordine di esecuzione non è
  quello dei numeri: S9, la filogenesi, è attiva solo con `phylo.enabled`, dipende da
  S13 e si esegue fra S13 e S14, perché l'albero si costruisce sulle varianti che
  restano dopo i filtri finali; S14 dipende da S13 e, solo con la filogenesi attiva,
  da S9; S10 non dipende da S9. Disattivata, S9 non è lavoro mancante e S14 non la
  richiede; attivarla a catena conclusa rifà le sole S9 e S14. Il numero della fase e
  la sua cartella `09_phylogeny` restano quelli del piano. S12 precede
  obbligatoriamente S13.
  Ogni fase conclusa scrive nella propria cartella un manifesto suo, `manifest_S<n>.json`,
  così due fasi che condividono una cartella si concludono separatamente. Il manifesto
  registra anche su che cosa la fase è stata calcolata: l'impronta della
  configurazione, quella di ogni fase a monte e, per S0, quella dei dati di ingresso.
  Una fase è conclusa se il manifesto c'è, i suoi artefatti sono integri e questi
  ingressi coincidono con quelli di adesso: cambiare un parametro o ricalcolare una
  fase a monte la rende da rifare, insieme a tutte quelle che ne dipendono. Ogni fase
  dichiara i parametri da cui dipende, e solo quelli contano: cambiare `err.nbases`
  rifà S3 e ciò che segue, non S0, S1 e S2. La dichiarazione è obbligatoria, perché
  una fase che non la fa non si registra, ed è vincolante: il codice Python di una
  fase vede la configurazione attraverso una vista ristretta ai parametri dichiarati
  (`config/vista.py`), che solleva un errore su qualunque altro accesso, e gli script
  R ricevono soltanto i parametri che la fase passa e falliscono se ne leggono uno
  non ricevuto. Una dipendenza dimenticata diventa così un errore al primo test che
  esercita la fase, non un risultato obsoleto. Le eccezioni sono G15, che verifica
  la coerenza dell'intera configurazione, e G12, che verifica il riferimento
  tassonomico contro `tax.ref_md5`: sono precondizioni, S0 li esegue sulla
  configurazione completa prima di vederla ristretta, e i parametri che servono solo a
  loro non entrano nella sua impronta, perché si ripetono comunque a ogni avvio.
  Cambiare una soglia di S6 rifà S6 e ciò che segue, non S0; cambiare il riferimento
  tassonomico rifà S8 e ciò che segue, perché il riferimento entra nell'impronta di S8
  con `tax.ref_md5`. Del gruppo `ctrl` S0 dichiara le sole chiavi che legge (la
  colonna, le tre etichette delle classi e la regola di riclassificazione in
  controllo negativo): i parametri dei controlli positivi sono di S11, e cambiarli
  rifà soltanto S11. I parametri che non incidono sui
  risultati, dichiarati in un solo elenco in `config/resolve.py` (`run.threads`,
  `run.keep_filtered_fastq`, `run.batch_size`, `io.out_root`, `retry.enabled`,
  `retry.max_attempts`), sono leggibili da ogni fase e non entrano in nessuna impronta: cambiarli, o spostare
  la cartella di un'esecuzione conclusa, non la rende da rifare. Il digest scritto in
  `00_config/resolved.yaml` resta calcolato sull'intera configurazione.
  Una fase avviata prima delle fasi da cui dipende solleva `E-GRAFO-01`, oppure
  `E-S13-01` se a mancare è la decontaminazione prima del filtro di prevalenza.
  `ProjectRun` legge questo stato dal disco in una valutazione, a cui si chiede quali
  fasi sono concluse, quali disattivate, quali da eseguire e quale è la prossima;
  dentro una valutazione il checksum di ogni artefatto è calcolato una volta sola, ma
  una valutazione completa li calcola tutti, e con gli artefatti di S2 e S4, da
  gigabyte, costa secondi (7,0 sul dataset di riferimento con
  S0-S8 e S10-S14 concluse); le quindici fasi esistono tutte come codice, e S9, la
  filogenesi, è disattivata per difetto;
- **la versione del calcolo e la provenienza** (`runner/provenienza.py`). Ogni fase
  dichiara una versione, che entra nell'impronta: si incrementa quando cambia ciò che
  la fase calcola, e la ripresa rifà allora quella fase e le successive. Il registro
  `src/amplicon16s/steps/registro_sorgente.json` riporta per ogni fase la versione e
  l'impronta del suo sorgente (il modulo, i moduli in cui vive il calcolo, lo script R e
  i file di `R/lib` che carica). I moduli del calcolo si ricavano dagli import del
  modulo della fase, direttamente o attraverso altri moduli, esclusi i moduli delle
  altre fasi, da cui si importano solo nomi di artefatti (un test lo verifica), e i
  moduli dell'infrastruttura elencati uno per uno, con la ragione, in
  `runner/provenienza.py` (`ESCLUSI`). Il criterio è chi può cambiare ciò che una
  fase calcola: restano nel sorgente la scrittura degli artefatti e i checksum, il
  ponte verso R (richiesta e ambiente del processo) e il calcolo dei valori corretti
  di un nuovo tentativo; una funzione condivisa fra fasi, come il lettore del
  tracciamento, sta in un modulo comune ed entra nel sorgente di ciascuna. Entra
  anche ciò che il contesto passa alle fasi senza import (`DAL_CONTESTO`): il
  lettore dell'inventario (`metadata/lettura_inventario.py`), nel sorgente di ogni
  fase che dipende da S0; e, per S0 soltanto, il catalogo degli errori, la cui
  sintesi compare nel dettaglio di una violazione di G15 scritto in `gates.json`
  (`ECCEZIONI_PER_FASE`); un test fallisce se non corrisponde al codice, e
  `scripts/registro_sorgente.py --aggiorna` lo aggiorna solo se, per ogni fase
  modificata a parità di versione, si dichiara la modifica senza effetto con
  `--senza-effetto`. Il manifesto di ogni fase registra la provenienza: versione,
  impronta del sorgente, commit del repository se disponibile, immagine dichiarata in
  `run.container`. L'impronta del sorgente si calcola sui file effettivamente usati:
  i moduli importati e gli script R della cartella che il ponte esegue
  (`AMPLICON16S_R_DIR`), non le copie del repository, così un'immagine che esegue
  script diversi risulta tale nel manifesto. A parità di versione una provenienza
  diversa produce un avviso nel log e nel report, non un ricalcolo;
- **la regola rigorosa sulla provenienza** (`run.strict_provenance`, falso per difetto,
  vero nella configurazione congelata di OSD-734). Con la regola attiva l'esecutore,
  dopo G15, G14 e G12 e prima di qualunque fase, verifica che il repository git sia
  leggibile, che `src/` e `R/` non abbiano modifiche non committate né file non
  tracciati e che gli script R eseguiti siano quelli del repository; poi legge dalle
  librerie installate (`R/00_ambiente.R`) la versione di R e quelle di tutti i
  pacchetti di `run.lockfile`, compresa la correzione di dada2 con la sua impronta, e
  le confronta con il file. Se una verifica fallisce l'esecuzione non parte:
  `E-PROV-01` (nessun commit leggibile identifica il codice: git non leggibile non è
  un via libera), `E-PROV-02` (modifiche non committate), `E-PROV-03` (ambiente R
  diverso dal file di blocco), tutti a revisione umana. Con la regola l'impronta del
  sorgente della fase, quella del file di blocco e l'immagine dichiarata entrano in
  ciò su cui la fase è calcolata: una differenza rende la fase da rifare, non solo da
  segnalare, e un'esecuzione calcolata senza la regola va ricalcolata per intero
  quando la si attiva. L'immagine resta una dichiarazione: da dentro il container il
  digest dell'immagine in esecuzione non è conoscibile, e non sono verificate le
  librerie di sistema né i pacchetti Python di terze parti. L'esito della regola è
  registrato a ogni avvio e riportato dal report, che distingue ciò che è stato
  verificato da ciò che è dichiarato;
- **l'analisi di sensibilità e il registro delle decisioni** (`scripts/sensitivity.py`,
  `scripts/confronto_modelli_filogenesi.R`, `docs/`). `docs/decision_log.md` registra
  la regola con cui si scelgono i valori dei parametri sottoposti ad analisi di
  sensibilità, fissata in una revisione che precede ogni risultato, l'esito
  dell'analisi e le decisioni di metodo del progetto con motivazione e misure.
  `scripts/sensitivity.py` copia un'esecuzione conclusa con collegamenti fisici, la
  riprende cambiando un parametro alla volta, misura campioni conservati, varianti
  finali, letture, soglie per piastra con la loro origine e contaminanti, e applica la
  regola. Non porta nel codice alcun valore di un dataset: le griglie e le modalità
  vengono da un file (`--griglie`; per OSD-734 `docs/sensibilita/griglie.yaml`, con
  `katharoseq.target_sensitivity`, `decontam.threshold`, `prev.min_fraction` e la
  modalità di decontaminazione), il valore corrente dalla configurazione, e i
  percorsi relativi sono resi assoluti come dalla riga di comando. Un valore è non
  ammissibile solo se produce un arresto previsto dalla regola; ogni altro arresto
  ferma lo strumento; `docs/sensibilita.md` riporta parametri, intervalli, misure e
  decisioni, con le tabelle in `docs/sensibilita/`. Su OSD-734 la regola mantiene i
  quattro valori correnti, e il confronto fra GTR+G+I e GTR+G con il criterio
  d'informazione bayesiano mantiene GTR+G+I come predefinito di `phylo.model`. Con la
  configurazione congelata e la regola rigorosa attiva la catena completa, rieseguita
  da zero nell'immagine pubblicata, dà i 1.031 artefatti dei checksum pubblicati in
  `dati/osd734/`;
- **l'esecutore e la politica dei tentativi** (`runner/executor.py`,
  `runner/retry.py`). A ogni avvio, con `run` come con `resume`, l'esecutore ripete
  la verifica di coerenza della configurazione (G15), quella delle risorse della
  macchina (G14) e quella del riferimento tassonomico contro `tax.ref_md5` (G12),
  senza rieseguire S0: processori, spazio su disco e integrità del riferimento si
  controllano anche quando S0 è già conclusa, e un riferimento alterato dopo S0 ferma
  l'esecuzione con E-S0-12 prima di arrivare a S8. Poi esegue in ordine le fasi da
  eseguire. Chiamato con una fase di arrivo (`validate` arriva a S0), esegue quella
  fase e i soli suoi antenati secondo le dipendenze del grafo, non tutte le fasi che
  la precedono: arrivare a S2 non esegue S1, da cui S2 non dipende; le fasi non
  richieste restano da eseguire, e una ripresa completa le esegue. Una fase
  fallita si comporta secondo la categoria del suo codice: a revisione umana ci si
  ferma subito; un codice ripetibile si ritenta solo se è in `retry.whitelist` e
  `retry.enabled` è vero, entro `retry.max_attempts` tentativi totali compreso il
  primo, e solo se la fase dichiara per quel codice un'azione correttiva, perché
  ritentare identico darebbe lo stesso esito. L'azione correttiva cambia un parametro
  in una copia in memoria della configurazione (oggi `run.batch_size` o `err.nbases`)
  e il manifesto della fase registra il codice, il valore dichiarato e quello usato.
  Per un errore transitorio la fase può invece dichiarare esplicitamente, per quel
  codice, un nuovo tentativo senza modifiche, registrato come ogni altro
  aggiustamento: è il caso di E-S2-03, errore di lettura;
  la fase resta giudicata sulla configurazione dichiarata, quindi una ripresa non la
  rifà. Le degradazioni non fermano l'esecuzione: la fase le registra e proseguono nel
  log e nel manifesto. Sono la categoria «degradazione automatica» del catalogo, fra
  cui E-S0-15 (G08, una piastra con pochi controlli negativi), E-S1-01 (scarto del
  troncamento sotto la lettura più corta), E-S6-02 (frazione chimerica oltre
  l'avviso), E-S11-02 (una piastra senza curva propria usa la soglia aggregata o la
  mediana), E-S11-05 (nessuna curva valida, nessuna soglia), E-S11-04 (controlli
  positivi conformi sotto il minimo, con `ctrl.positive_gate` falso), E-S13-03
  (campioni svuotati dal filtro tassonomico).
  Quando l'esecuzione si ferma, l'esecutore dichiara il punto di ripresa: fase,
  codice, messaggio del catalogo, tentativi fatti e comando per ripartire. Delle fasi
  con codici ripetibili esistono oggi S2, il cui E-S2-03 è provato su un archivio
  davvero corrotto, S3, il cui E-S3-01 è provato su una stima che davvero non
  converge, S4, il cui E-S4-02 è provato limitando davvero la memoria del processo R,
  e S5. Una fase può anche dichiarare, nell'errore, che un nuovo tentativo non
  cambierebbe l'esito: allora non si ritenta, e il motivo compare nel punto di
  ripresa. È il caso di E-S5-01, sempre, e di E-S4-02 con `dada.pool` vero;
- la registrazione degli eventi su due uscite: la console per chi segue l'esecuzione e
  un file JSON Lines con rotazione sotto `99_logs`, in un formato che si interroga per
  codice, fase o categoria invece di doversi leggere;
- la catena di integrazione continua, con due job. Il primo installa il pacchetto con
  Python 3.11, R 4.5.2 e jsonlite 2.0.0 ed esegue la suite di test: è il riscontro
  rapido, e vi girano davvero i test del ponte verso R. Il secondo costruisce
  l'immagine della pipeline dallo stesso Dockerfile ed esegue la suite al suo interno,
  così i test delle fasi di calcolo usano i pacchetti Bioconductor alle versioni del
  container, verificate contro `renv.lock` durante la costruzione; costa alcuni minuti
  per push. In entrambi l'assenza di R, e nel secondo quella di Bioconductor, fa
  fallire i test che li richiedono invece di saltarli;
- **i dati di OSD-734 recuperabili dal repository** (`dati/osd734/`): le tabelle ISA
  usate dalla pipeline, l'elenco delle 960 corse con i checksum di ENA, fonte e licenza
  di ogni file, una configurazione con i percorsi relativi e tre script, con la sola
  libreria standard: `scarica_letture.py` (ENA) e
  `dati/riferimento/scarica_riferimento.py` (SILVA 138 da Zenodo e l'elenco dei suoi
  taxa difettosi, in una cartella comune ai dataset, con `FONTI.tsv` e
  `impronte.md5`) scaricano verificando l'MD5 e riprendendo uno scarico interrotto; `ricostruisci_lotto.py` ricostruisce il file del
  lotto dai metadati Qiita pubblicati dagli autori, fissati a un commit. Il README
  della cartella riporta il comando per eseguire la pipeline nel container e i
  requisiti misurati;
- **la riproduzione da zero dei risultati su OSD-734**, eseguita il 5 ottobre 2026 come
  la eseguirebbe un utente esterno: repository clonato da GitHub (commit `59ef39d`),
  immagine costruita dal suo Dockerfile senza cache, 960 letture scaricate da ENA (2,52
  GB in 2 ore e 48 minuti), SILVA 138 da Zenodo e file del lotto ricostruito, ogni file
  verificato contro il checksum della fonte; catena S0-S14 in 53 minuti e 20 secondi,
  con 12,3 GB di memoria di picco del container. Tutti i 1.031 artefatti dei manifesti
  di fase sono identici, byte per byte, a quelli dell'esecuzione di riferimento
  precedente. I checksum attesi sono pubblicati in `dati/osd734/`
  (`checksum_finali.sha256` per i file consegnati, `checksum_artefatti.tsv` per tutti
  gli artefatti), e `confronta_risultati.py` li confronta con un'esecuzione indicando,
  se differiscono, la prima fase in cui la differenza compare; legge i manifesti nelle
  sole cartelle delle fasi e in `12_final/intermedi`, e si ferma se una fase ne ha due. Si riferiscono
  all'immagine pubblicata nel registro, con cui la catena completa è stata rieseguita
  (1.031 artefatti identici, 48 minuti e 11 secondi, 12,6 GB di picco). La prova ha
  fatto emergere un difetto, corretto: i percorsi relativi della
  configurazione fermavano S1, perché il processo R di una fase parte nella cartella
  della fase; ora la riga di comando rende assoluti, rispetto alla cartella di lancio,
  i percorsi relativi della configurazione (`cli.py`), e la configurazione registrata
  riporta i percorsi dei file letti;
- **il sottoinsieme di prova** (`tests/fixtures/osd734/`, `scripts/build_test_subset.py`).
  Ventotto campioni del dataset di riferimento, scelti per le fasi successive: due
  piastre, una per corsa, con i cinque controlli negativi che `decontam.min_blanks`
  richiede e biologici con cui confrontarli; una serie completa degli otto livelli di
  diluizione dei controlli positivi, più un positivo anomalo; il biologico più povero
  e il più profondo; la lunghezza minima globale di 137 bp; JLP1A1.L4, con il motivo
  conservato nel 5,8% delle letture; un tubo non aperto, un campione d'aria e una
  superficie dell'Airlock. La selezione, con il motivo di ogni campione e gli MD5 dei
  file originali, sta nel repository; lo script ne ricostruisce dai dati locali la
  versione completa, verificando gli MD5, e la versione ridotta, sottocampionata con
  un seme fisso, che vive nel repository (circa 2 MB) e su cui girano i test;
- **la fase S1, il profilo delle letture** (`steps/s01_profile.py`,
  `R/01_profile.R`), la prima fase di calcolo e la prima a usare Bioconductor
  attraverso il ponte. Legge ogni file per intero e scrive in `02_qc_profiles/` la
  distribuzione delle lunghezze, il profilo di qualità per posizione e il conteggio
  delle letture di ogni campione, in tabelle a una riga per campione e valore, con un
  riepilogo sull'intero insieme. Sul dataset di riferimento, nel container, la
  lunghezza minima è 137, la moda 151, e la qualità mediana non scende sotto 25 in
  nessuna posizione; la fase impiega alcuni minuti. S1 chiude il limite noto di G09,
  che misura le prime `qc.head_reads` letture: su tutte le letture calcola, per
  classe, la frazione più corta di `filter.truncLen`, e si ferma con `E-S1-02` se
  quella dei campioni biologici, o quella dei controlli positivi, supera
  `qc.max_frac_short_reads` (0,05; ogni classe si giudica da sola, e i controlli
  negativi non contano). Una sola
  lettura corta non ferma più né G09 né S1: in un dataset a lunghezza variabile ce
  n'è quasi sempre qualcuna, e il filtro la toglie. S1 conta anche i valori di
  qualità distinti (`valori_qualita.tsv`) e, se sono quattro o meno (qualità
  raggruppate), lo dichiara con `E-S1-03`, perché la stima standard del modello di
  errore vi è poco affidabile: `err.error_function` sceglie fra la funzione standard
  di dada2 (`loess`) e una variante pesata e monotona (`loess_monotono`,
  `R/lib/errore_loess.R`: loess di primo grado con span 2 e pesi pari al logaritmo
  delle basi osservate, poi tassi non crescenti con la qualità; regge due, tre o
  quattro valori di qualità), usata da S3 e da S4. Una corsa senza letture filtrate
  ferma S3 con `E-S3-03`; letture con un solo valore di qualità, da cui nessuna
  funzione può stimare il modello, con `E-S3-04` (e in quel caso `E-S1-03` non
  suggerisce una funzione); ogni altra stima fallita, con più valori di qualità,
  con `E-S3-05`, che riporta il messaggio originale della funzione di stima. S1 chiude
  anche il limite noto di G07: prima dei profili conta su tutte le letture di ogni
  file i segni di coppia che G07 cerca nelle prime (`/1` e `/2`, `1:N:` e `2:N:`, lo
  stesso nome ripetuto), con la stessa funzione e la stessa tolleranza
  (`io_layer/reads.py`), li registra per campione in `coppie.tsv` e si ferma con
  `E-S1-04` se un file contiene le due letture di ogni coppia, o le sole seconde: un
  file con tutte le prime letture seguite da tutte le seconde supera G07 quando le
  letture ispezionate non arrivano al secondo blocco. Registra inoltre E-S1-01, lo scarto del troncamento sotto il minimo oltre
  `filter.truncLen_shortfall_warn`, che prima registrava S0 dalla stima;
- **la fase S2, filtro e troncamento** (`steps/s02_filter.py`, `R/02_filter.R`), con
  `dada2::filterAndTrim` e i parametri del gruppo `filter`; scrive in `03_filtered/`
  le letture filtrate di ogni campione, tutte di `filter.truncLen` basi, e le letture
  in ingresso e in uscita per campione. Prima del filtro decomprime per intero ogni
  archivio: un gzip troncato non fa fallire `filterAndTrim`, che lo legge fino a dove
  arriva, e un archivio incompleto diventa E-S2-03, ritentato senza modifiche e poi
  fermato. I controlli sul risultato si applicano ai campioni biologici e ai
  controlli positivi, non ai negativi, perché un bianco azzerato è un bianco pulito:
  E-S2-01 ferma quando i campioni azzerati superano `qc.max_zeroed_samples` (0),
  E-S2-02 quando la frazione media di letture scartate supera
  `qc.max_frac_lost_filter` (0,30). Sul dataset di riferimento, nel container, il
  filtro conserva in media il 98,7% delle letture in ogni classe e non azzera alcun
  campione. `run.keep_filtered_fastq` (vero per difetto) conserva le letture filtrate;
  con falso vengono rimosse solo quando tutte le fasi sono concluse, e la rimozione è
  registrata accanto al manifesto di S2, così la valutazione dello stato non la
  scambia per un artefatto perso e una ripresa non rifà nulla; se una fase che le
  legge deve poi essere ripetuta, S2 torna da eseguire prima di lei;
- **la fase S3, il modello d'errore** (`steps/s03_learn_errors.py`,
  `R/03_learn_errors.R`), con `dada2::learnErrors` e i parametri del gruppo `err`: un
  modello per corsa di sequenziamento, dalla colonna `err.batch_column` del file di
  arricchimento del lotto; con la colonna nulla, o senza quel file, un solo modello su
  tutti i campioni. Scrive in `04_error_models/` per ciascun modello il modello stesso e
  il grafico diagnostico, errori osservati e stimati in funzione della qualità, e per
  S4 la corrispondenza fra ciascun campione e il suo modello, `corrispondenza.tsv`.
  L'ordine in cui i campioni entrano nella stima viene dal seme `run.seed`, e S3
  registra quali campioni e quante basi di ciascuna classe l'hanno fatta. Modelli e
  grafici sono identici byte per byte fra due esecuzioni identiche: il grafico è un PNG,
  che non registra la data. La convergenza si riconosce da una proprietà del modello,
  non dal testo di un avviso; un modello che non converge è E-S3-01, ritentato con
  `err.nbases` raddoppiato solo se la stima non usava già tutte le basi della corsa.
  Un campione senza corsa, con la colonna attiva, è E-S3-02;
- **la fase S4, l'inferenza delle varianti** (`steps/s04_dada.py`, `R/04_dada.R`), con
  `dada2::dada` e i parametri del gruppo `dada`: ogni campione usa il modello d'errore
  della sua corsa, letto da `corrispondenza.tsv` di S3. Il pseudo-pooling è realizzato
  in due passate esplicite, entrambe a lotti di `run.batch_size` campioni: la prima
  elabora ogni campione da solo e conserva soltanto, per ogni sequenza, in quanti
  campioni compare e quante letture ha; la seconda rielabora ogni campione con le
  informazioni a priori, scelte con la regola del sorgente di dada2 1.36.0 (presenti
  in almeno `PSEUDO_PREVALENCE` campioni, 2, o con almeno `PSEUDO_ABUNDANCE` letture,
  infinito). Come in `dada(pool = "pseudo")`, la seconda passata usa il modello
  d'errore che dada ricalcola dalle transizioni della prima, uno per modello di S3.
  Le informazioni a priori si raccolgono su tutti i campioni, di entrambe le corse.
  Il risultato è quello di `dada(pool = "pseudo")` in una sola chiamata e non dipende
  dal lotto: i test lo verificano byte per byte, e per questo `run.batch_size` è fuori
  dall'impronta e il retry di E-S4-02 che lo dimezza è legittimo. Scrive in
  `05_asv_inference/` le varianti di ogni campione, le informazioni a priori, i modelli
  della seconda passata e un riepilogo. Sul dataset di riferimento, nel container con
  12 processori, impiega 26 minuti e trova 13.130 varianti distinte, con 4.090
  informazioni a priori su 12.857 sequenze della prima passata; due esecuzioni
  complete danno gli stessi byte in ogni artefatto;
- **la fase S5, la tabella delle sequenze** (`steps/s05_seqtab.py`, `R/05_seqtab.R`),
  in `06_seqtab/`: è la tabella prima della rimozione delle chimere. Il numero di
  varianti si controlla contro `qc.max_asv_count` prima di allocarla; oltre, o con la
  memoria esaurita, è E-S5-01, che la fase dichiara inutile ritentare, perché la
  tabella è una matrice densa campioni x varianti la cui memoria non dipende dal lotto;
- **la fase S6, la rimozione delle chimere** (`steps/s06_chimera.py`,
  `R/06_chimera.R`), con `dada2::removeBimeraDenovo` e i parametri del gruppo `chimera`,
  in `07_chimera/`. Non ricopia la tabella di S5: ne registra il riferimento con il
  checksum ed elenca le varianti chimeriche tolte. Riporta per ogni classe la frazione
  chimerica sulle letture e sulle varianti; i controlli guardano la frazione delle
  letture di biologici e positivi: E-S6-02 oltre `qc.warn_frac_chimeric` (0,25) si
  registra e si prosegue, E-S6-01 oltre `qc.stop_frac_chimeric` (0,50) ferma. Sul
  dataset di riferimento le letture chimeriche sono lo 0,44% nei biologici, lo 0,02%
  nei positivi e lo 0,04% nei negativi, contro l'8,8% delle varianti nei biologici;
- **la fase S7, il filtro di lunghezza** (`steps/s07_asv_length.py`,
  `R/07_asv_length.R`), che tiene le varianti fra `asv.len_min` e `asv.len_max` e
  scrive la tabella delle varianti accanto a S6. Con letture troncate a lunghezza fissa
  le varianti hanno tutte la stessa lunghezza, e sul dataset di riferimento il filtro
  non toglie nulla: 12.045 varianti di 137 basi. Che tolga le varianti fuori intervallo
  è verificato con varianti sintetiche. Una tabella che resta senza varianti, perché
  quella di S6 era già vuota o perché tutte cadono fuori intervallo, ferma la fase con
  E-S7-01, a revisione umana: le fasi successive fallirebbero più avanti con un errore
  che non ne indica la causa;
- **la fase S8, l'assegnazione tassonomica** (`steps/s08_taxonomy.py`,
  `R/08_taxonomy.R`), con `dada2::assignTaxonomy`, il classificatore bayesiano naive
  (`tax.classifier: naive_bayes`: il training set di SILVA per IdTaxa non è distribuito
  da alcuna fonte), sul riferimento `tax.ref_fasta` (SILVA 138), in `08_taxonomy/`.
  `tax.min_boot` vale 50 e `tax.assign_species` è falso, perché le letture coprono 137
  delle circa 253 basi dell'amplicone e su sequenze corte il bootstrap è più basso a
  parità di correttezza. Scrive la tabella tassonomica fino al genere, il bootstrap di
  ogni rango, la copertura del phylum per classe di campioni e un riepilogo. Il seme
  del bootstrap viene da `run.seed`; con la correzione di dada2 descritta sopra due
  esecuzioni danno gli stessi byte, anche con un numero di thread diverso, e
  `run.threads` resta fuori dall'impronta. Il riferimento non si modifica: se
  `tax.ref_bad_taxa` indica l'elenco dei taxa con un difetto noto (per SILVA 138
  versione 2, un rango mancante in 10 famiglie e 114 generi, che fa comparire il nome
  nella colonna del rango superiore), `difetto_riferimento.tsv` marca le assegnazioni
  che vi ricadono, in qualunque colonna, e il riepilogo conta varianti e letture
  interessate. E-S8-02 ferma la fase se nei campioni biologici o nei controlli
  positivi la frazione di varianti con il phylum è sotto `qc.min_frac_phylum` (0,80).
  Sul dataset di riferimento S8 impiega circa 5 minuti con 12 thread. Il phylum è
  assegnato al 99,0% delle varianti dei biologici, al 99,1% di quelle dei controlli
  positivi e al 99,4% di quelle dei negativi, con oltre il 99,9% delle letture in
  ogni classe. Le assegnazioni sui taxa col difetto noto riguardano 320 varianti e
  506.697 letture (1,6%), quasi tutte anaerobi dell'ordine
  Peptostreptococcales-Tissierellales (Anaerococcus, Finegoldia, Peptoniphilus). Nei
  controlli positivi la variante dominante è Variovorax, il ceppo della serie, ai due
  livelli di diluizione più alti in tutte le dieci piastre; ai livelli più bassi
  prevalgono i contaminanti di reagente, come atteso nella serie KatharoSeq;
- **la fase S10, l'oggetto integrato** (`steps/s10_phyloseq.py`, `R/10_phyloseq.R`,
  `R/lib/oggetto.R`), in `10_phyloseq/`: un oggetto phyloseq, `ps_integrato.rds`, con
  quattro componenti allineati, la tabella dei conteggi di S7, la tabella tassonomica
  di S8, i metadati dei campioni e le sequenze di riferimento delle varianti; l'oggetto
  integrato non porta mai l'albero, che S9 costruisce sulle varianti finali e S14
  aggiunge all'oggetto finale. L'oggetto contiene
  tutti i campioni dell'inventario, nell'ordine dell'inventario: un campione rimasto
  senza letture, che non ha una riga nelle tabelle di S5-S7, vi entra con conteggi a
  zero e il riepilogo lo elenca (sul dataset di riferimento non ce ne sono; il caso è
  verificato togliendo un controllo negativo dalla tabella di S7). I campioni sono
  identificati dall'accession (`out.sample_id_source`), le varianti da `ASV1`, `ASV2`,
  ... (`out.asv_id_scheme: abundance_rank`): letture decrescenti nei campioni
  biologici, a parità letture in tutti i campioni, a parità ancora la sequenza in
  ordine di byte, così la numerazione è deterministica. Il denominatore sono i
  biologici perché a questo punto l'oggetto contiene ancora i controlli. Gli
  identificativi non si rinumerano nelle fasi successive, e la sequenza resta
  l'identità stabile della variante. Con `out.taxa_are_rows` vero le varianti stanno
  sulle righe: lo script lo verifica a ogni esecuzione sui nomi di righe e colonne,
  identificativi delle varianti e accession, non sulle dimensioni, che una matrice
  quadrata trasposta non distinguerebbe; e verifica di nuovo ogni componente
  dell'oggetto assemblato, perché `phyloseq()` tiene in silenzio solo campioni e
  varianti comuni. Un disallineamento ferma con `E-S10-01`. I metadati portano le
  colonne dell'inventario (accession, nome, classe, materiale, posizione, modulo,
  piastra, corsa) e quelle indicate in `out.study_columns` e `out.batch_columns`, con
  nomi sintattici scelti dalla fase, che R non altera; la corrispondenza con le
  colonne originali è in `colonne_metadati.tsv`, e una colonna assente, ripetuta
  nell'intestazione o con un nome già usato ferma con `E-S10-02`. Il bootstrap di
  ogni rango e la marcatura dei taxa col difetto noto del riferimento, che non hanno
  uno slot nell'oggetto, stanno in `varianti_accessorie.tsv`, per identificativo di
  variante; `varianti.tsv` dà sequenza e letture di ogni identificativo. Due
  esecuzioni danno gli stessi byte, oggetto compreso. Sul dataset di riferimento S10
  impiega 6 secondi: 960 campioni, 12.045 varianti, 30.912.706 letture, oggetto di
  1,1 MB; ASV1 è un Pseudomonas, con 4.114.563 letture nei 770 biologici;
- **la fase S11, la validazione della corsa dai controlli positivi**
  (`steps/s11_controls.py`, `R/11_controls.R`, `R/lib/katharoseq.R`), in
  `11_controls/`, sull'oggetto di S10. Per ogni controllo positivo misura la
  fedeltà, la frazione delle sue letture assegnate al taxon atteso
  (`katharoseq.target_taxon`, sommando le varianti al rango
  `katharoseq.collapse_rank`), e la lega alla profondità con la sigmoide allosterica
  di KatharoSeq, f = x^h / (k' + x^h) con x = log10(profondità), stimata con `nls`
  nella forma equivalente con k' = x50^h, da valori iniziali ricavati dai dati senza
  numeri casuali. La bontà è R² = 1 − SS_res / SS_tot sulla scala della fedeltà, che
  per un modello non lineare dice solo quanto la curva riduce l'errore rispetto alla
  media. La profondità minima è quella a cui la curva raggiunge
  `katharoseq.target_sensitivity` (0,90). Il livello di diluizione viene dai dati, la
  colonna `katharoseq.cell_count_column` (del file del lotto o della tabella di studio), non dai nomi dei
  campioni. Si adattano la curva aggregata e una per piastra sugli stessi punti (i
  controlli utilizzabili delle piastre con almeno `ctrl.min_positives` punti), e si
  preferisce il modello con il criterio di Akaike minore (quello per piastra è la
  somma delle sue curve); una curva vale se ha almeno
  `ctrl.min_positives` punti, R² non inferiore a `katharoseq.min_r2`, la soglia
  dentro le profondità osservate e determinata dai dati: almeno un'osservazione deve
  cadere fra il punto medio della curva, dove la fedeltà vale 1/2, e la soglia,
  altrimenti la posizione della soglia verrebbe solo dalla forma del modello. La soglia è un artefatto, `soglia.json`, per
  piastra e per i campioni senza piastra, e ogni valore porta lo stadio, l'origine e
  il motivo: ogni soglia vale sulle letture dell'oggetto integrato (senza chimere,
  `katharoseq.read_stage`), e nessuna sulle letture grezze. Se l'aggregato è
  preferito ed è valido, tutte le piastre usano la sua soglia (origine `aggregata`);
  altrimenti una piastra con la curva valida usa la propria (`propria`), e una senza
  usa l'aggregata se valida, o la mediana delle soglie proprie arrotondata all'intero
  superiore (`mediana`); un campione senza piastra segue la stessa strada. Chi non usa
  una curva propria è dichiarato con `E-S11-02`, con origine e motivo (bontà
  insufficiente, soglia non determinata, controlli insufficienti). Se nessuna curva è
  valida (nessun controllo positivo, colonna dei livelli assente, nessuna piastra con
  abbastanza punti) non c'è alcuna soglia (`nessuna`, `E-S11-05`) e in S13 resta il
  solo `qc.min_reads_final`; con `qc.min_reads_mode: none` le curve si riportano come
  diagnostica e nessuna soglia si applica. `soglia.json` registra anche gli AIC dei
  due modelli, i punti comuni, l'R² dichiarato (dell'aggregato, o delle sole piastre
  che usano la propria curva) e la mediana, arrotondata e no. Un controllo è conforme se fedeltà e profondità non
  si discostano da quelle dei controlli dello stesso livello di concentrazione
  nelle altre piastre (z modificato entro 3,5): ai livelli diluiti la fedeltà attesa
  è bassa e i contaminanti non lo rendono non conforme. Sotto
  `ctrl.min_positive_pass_frac` conformi la fase si ferma con `E-S11-03` se
  `ctrl.positive_gate` è vero, altrimenti registra l'avviso `E-S11-04`; i non
  conformi non entrano nella curva. `profondita_campioni.tsv` dice per ogni campione
  se cadrebbe sotto la sua soglia, come misura: il filtro è di S13. Due esecuzioni
  danno gli stessi byte. Sul dataset di riferimento S11 impiega 5 secondi: 74
  controlli conformi su 80; il modello per piastra (R² complessivo 0,925) prevale
  sull'aggregato (0,566), con soglie fra 6.376 e 41.968 letture nelle piastre 3-10.
  Le piastre 1 e 2 ripiegano su 1.000 letture grezze: nella 1, dove un contaminante
  cloroplastico gonfia la profondità dei controlli diluiti, la curva è un gradino in
  un tratto senza osservazioni, fra 49.417 e 99.521 letture, e la soglia non è
  determinata; nella 2 l'R² è 0,532. Sotto la propria soglia cadrebbero 278 campioni
  biologici su 770;
- **la fase S12, la decontaminazione dai controlli negativi** (`steps/s12_decontam.py`,
  `R/12_decontam.R`, `R/lib/decontaminazione.R`), in `11_controls/` accanto a S11, con
  i propri file del ponte e il proprio manifesto. Con `decontam::isContaminant` per
  prevalenza (`decontam.method`; il metodo per frequenza richiede una concentrazione
  del DNA che i metadati non hanno) confronta i controlli negativi con i campioni
  biologici, secondo le classi dell'inventario; i controlli positivi restano fuori dal
  confronto. Una variante è un contaminante se la probabilità è sotto
  `decontam.threshold` (0,5, orientata alla sensibilità in bassa biomassa). Calcola
  due modalità: aggregata, tutti i biologici contro tutti i negativi, e per piastra,
  con le probabilità delle piastre combinate secondo `decontam.batch_combine`, il
  `batch.combine` di decontam, dichiarato al suo predefinito `minimum`; una piastra
  con meno di `decontam.min_blanks` negativi confronta i propri biologici con i
  negativi di tutte le piastre. Quale delle due decide i contaminanti lo dichiara
  `decontam.mode`, non l'esito; l'altra resta una diagnostica. Per OSD-734 è
  l'aggregata: le piastre hanno pochi negativi, e il minimo di dieci probabilità
  stimate ciascuna su così pochi negativi è permissivo per costruzione. Se la
  modalità dichiarata rimuove dai biologici una frazione delle letture oltre
  `qc.max_frac_contaminant` (0,40), `E-S12-02` ferma la fase; se nessun campione
  biologico ha letture la ferma `E-S12-04`, e se tutte le varianti risultano
  contaminanti entro la quota ammessa `E-S12-05`. Il confronto numerico e
  la frazione rimossa per classe sono in `decontam_riepilogo.json`; le probabilità e
  le prevalenze di ogni variante in `decontam_varianti.tsv`, quelle di ogni piastra
  in `decontam_per_piastra.tsv`, l'elenco dei contaminanti in
  `decontam_contaminanti.tsv`. I contaminanti escono dall'oggetto in tutti i campioni,
  `ps_decontaminato.rds`, e le altre varianti tengono il loro identificativo. Due
  esecuzioni danno gli stessi byte. La precedenza su S13 è nel grafo: avviare il
  filtro di prevalenza prima della decontaminazione solleva `E-S13-01`. Sul dataset di
  riferimento, con i 110 controlli negativi (i tamponi mai aperti compresi), S12
  impiega 34 secondi: l'aggregata trova 759 contaminanti, con il 4,9% delle letture
  dei biologici, l'87,1% di quelle dei negativi e il 9,7% di quelle dei positivi; per
  piastra, come diagnostica, sarebbero 1.138 e il 44,2% delle letture dei biologici.
  Il contaminante cloroplastico dei negativi delle piastre 1 e 2 è identificato
  nell'aggregata e, per piastra, nelle piastre 1 e 2; ASV1, un Pseudomonas con 4,1
  milioni di letture nei biologici, non è un contaminante (probabilità 1);
- **la fase S13, i filtri finali** (`steps/s13_filtri.py`, `R/13_filtri.R`), in
  `12_final/intermedi/`: l'oggetto filtrato e le tabelle dei filtri sono intermedi di
  calcolo, e in `12_final/` stanno solo i file consegnati da S14. Dall'oggetto
  decontaminato ricava l'oggetto dei soli campioni biologici: un controllo
  nell'oggetto finale sarebbe trattato come un campione ambientale. I controlli
  positivi e negativi li consegna S14, dall'oggetto integrato di S10, in
  `ps_controlli.rds`.
  I filtri, in quest'ordine: (1) **profondità**, sui campioni: ogni biologico si
  confronta con la soglia della sua piastra in `soglia.json` di S11 sulle letture
  senza chimere del tracciamento (S7), non con la profondità dell'oggetto
  decontaminato, che non è la grandezza su cui la soglia è stata stimata; dove S11 non
  dà una soglia il filtro non si applica, e `filtri_riepilogo.json` lo dichiara
  (`profondita_applicata`); (2) **tassonomico**, sulle varianti:
  senza phylum (`filt.remove_na_phylum`) e i taxa di `filt.exclude_taxa` cercati in
  ogni rango, confrontati senza il prefisso di rango che alcuni riferimenti portano
  (`p__`, `o__`); se la tassonomia non ha il rango Phylum il primo non si applica, e
  la fase lo dichiara con `E-S13-05`; (3) **prevalenza**, sulle varianti: resta una variante con almeno
  `prev.min_count` letture in almeno `ceil(prev.min_fraction × n)` campioni, dove
  `n` sono i biologici tenuti dopo i filtri 1 e 2, perché numeratore e denominatore
  si riferiscono allo stesso insieme e un campione sotto soglia di profondità ha una
  composizione non attendibile; (4) **letture finali**, sui campioni: sotto
  `qc.min_reads_final` il campione esce, senza fermare l'esecuzione. I due filtri
  sulle varianti commutano. Un campione svuotato dal filtro tassonomico non ha
  segnale batterico ed esce con la degradazione `E-S13-03`; uno svuotato dal filtro
  di prevalenza esce anch'esso con il motivo, e la fase lo dichiara con la
  degradazione `E-S13-02`: aveva letture batteriche, solo rare, e chi legge il
  manifesto deve sapere che la soglia lo ha tolto. Se nessun campione supera i filtri
  la fase si ferma con `E-S13-04`. Gli identificativi delle varianti non
  si rinumerano. I campioni esclusi, con il filtro e il motivo, sono in
  `esclusioni.tsv`; le varianti rimosse in `varianti_rimosse.tsv`; l'ordine, il
  denominatore e ciò che ciascun filtro ha tolto, per classe, in
  `filtri_riepilogo.json`;
- **la fase S14, la serializzazione** (`steps/s14_finale.py`, `R/14_finale.R`,
  `R/lib/export.R`), in `12_final/`: `ps_final.rds`
  (`out.serialization`) e, con `out.export_flat`, gli export piatti `conteggi.tsv`,
  `tassonomia.tsv`, `metadati.tsv` e `sequenze.fasta` (e `albero.nwk` con la
  filogenesi attiva), scritti senza numeri decimali,
  in UTF-8, con valori mancanti come campo vuoto; se l'inventario ha controlli,
  `ps_controlli.rds`; `checksum.sha256` riporta l'impronta di ogni file consegnato,
  `ps_controlli.rds` compreso, nella forma di `sha256sum`. Prima del calcolo S14
  toglie dalla cartella i file consegnati di un'esecuzione precedente (un albero o
  degli export che la configurazione non produce più). La validazione,
  `E-S14-01`: componenti allineati, solo biologici, nessun campione e nessuna variante
  senza letture, l'oggetto riletto e quello ricostruito dai soli export identici a
  quello serializzato, `ps_filtrato.rds` di S13 e `ps_integrato.rds` di S10 integri rispetto ai manifesti di chi li ha scritti, la
  frazione delle letture trattenute dall'insieme dei campioni finali (letture finali
  su letture senza chimere) non sotto `qc.min_frac_reads_retained` (0,40). Le due
  soglie del piano agiscono a livelli diversi: `qc.min_reads_final` per campione, in
  S13, `qc.min_frac_reads_retained` sull'insieme, in S14. Gli artefatti di S13 e S14
  sono identici byte per byte fra due esecuzioni e con impostazioni locali diverse
  (`LC_ALL=C` e `en_US.UTF-8`): l'intestazione del formato RDS registra la codifica
  della sessione, e ogni fase scrive i suoi file `.rds` con `salva_rds`
  (`R/lib/io_json.R`), sempre con una codifica UTF-8. Sul
  dataset di riferimento S13 impiega 9 secondi e S14 7: il filtro per profondità
  esclude 278 biologici su 770 (277 nelle piastre 3-10, sulle letture senza chimere;
  uno nella piastra 1, sul ripiego delle letture grezze), gli altri nessuno; il
  tassonomico toglie 624 varianti su 11.286 (il 3,4% delle letture dei biologici,
  l'1,1% di quelle dei negativi), la prevalenza, con denominatore 492 e quindi almeno 5 campioni, 8.909
  varianti (l'1,0%). L'oggetto finale ha 492 campioni e 1.753 varianti, con 21,2
  milioni di letture, il 92,2% delle loro letture senza chimere; il campione più
  povero ne ha 1.085;
- **la fase S9, la filogenesi opzionale** (`steps/s09_phylogeny.py`,
  `R/09_phylogeny.R`, `R/lib/albero.R`), in `09_phylogeny/`. È disattivata per difetto
  (`phylo.enabled` falso): su letture di 137 basi, di un solo tratto del gene, un
  albero costruito da zero è debolmente risolto, e le relazioni profonde non sono
  sostenute dal dato. Attivata, costruisce l'albero delle varianti dell'oggetto
  filtrato di S13, quelle consegnate: sulle varianti di S7 sarebbero 12.045 sul
  dataset di riferimento, oltre `phylo.max_seqs` (5.000), contro le 1.753 finali.
  **Le due guardie**, a revisione umana, sono nella fase Python, prima di avviare il
  calcolo, sul numero di varianti finali: `E-S9-01` se supera `phylo.max_seqs`,
  `E-S9-02` se sono meno di quattro, perché con tre foglie un albero non ha topologia
  da stimare.
  **Il metodo**: allineamento multiplo con `DECIPHER::AlignSeqs` (`phylo.aligner`,
  `decipher`); massima verosimiglianza con phangorn sul modello `phylo.model`
  (`GTR+G+I`: GTR, gamma a quattro categorie, siti invarianti), a partire da un albero
  neighbor-joining su distanze JC69, ottimizzando frequenze, tassi, forma della
  gamma, quota di invarianti, lunghezze dei rami e topologia. **La ricerca è
  deterministica per scelta**: la topologia si migliora con scambi fra rami vicini
  (NNI), senza la perturbazione casuale (`stochastic`, `ratchet`) che phangorn offre
  per uscire dai massimi locali; il seme viene comunque da `run.seed`, e
  `run.threads` riguarda il solo allineamento. È una salita locale: può fermarsi su
  un albero con verosimiglianza un poco più bassa di quello che una ricerca casuale
  troverebbe, e dipende dal punto di partenza, dichiarato nello script. L'albero,
  l'allineamento e il riepilogo sono identici byte per byte fra due esecuzioni e con
  un numero di thread diverso. **Il radicamento è al punto medio**
  (`phangorn::midpoint`): fra le varianti non c'è un gruppo esterno, e le misure di
  diversità che usano la filogenesi richiedono un albero radicato; assume tassi
  simili nei due rami principali. Scrive `albero.nwk`, in formato Newick con le
  foglie che portano gli identificativi delle varianti e le lunghezze dei rami a otto
  cifre significative, `allineamento.fasta` e `filogenesi.json`, con metodo,
  verosimiglianza e parametri stimati. **L'albero nell'oggetto finale**: con la
  filogenesi attiva S14 legge `albero.nwk`, verificato contro il manifesto di S9,
  controlla che le foglie siano esattamente gli identificativi delle varianti
  dell'oggetto e che l'albero sia radicato (`E-S14-01` altrimenti), lo aggiunge a
  `ps_final.rds` senza cambiare l'ordine delle varianti e, con `out.export_flat`, lo
  esporta in `12_final/albero.nwk`, elencato in `checksum.sha256`. Sul dataset di
  riferimento, su una copia dell'esecuzione, S9 impiega 2 minuti e 36 secondi con 12
  thread per 1.753 varianti (allineamento di 145 colonne): il limite di 5.000 lascia
  margine. Sulla stessa copia, rieseguite S10-S14 con la filogenesi disattivata, i 29
  artefatti sono identici a quelli del riferimento;
- il tracciamento delle letture (`runner/tracciamento.py`): ogni fase registra i propri
  passi in file suoi, `letture_<passo>.tsv`, e la tabella completa si ricompone
  leggendo quelli delle fasi concluse nell'ordine del grafo, senza che una fase
  modifichi mai un file di un'altra. I passi sono uno per fase, salvo S2, che ne
  registra due: le letture
  grezze (S1), quelle in ingresso e in uscita dal filtro (S2), quelle attribuite a una
  variante (S4), quelle nella tabella (S5), quelle senza chimere (S6), quelle dopo il
  filtro di lunghezza (S7), quelle dopo la rimozione dei contaminanti (S12) e quelle
  nell'oggetto finale (S13, zero per i controlli e per i campioni esclusi);
- i quattro sottocomandi della riga di comando, descritti sopra, con i codici di
  uscita documentati;
- **il report di esecuzione** (`report/builder.py`, con i modelli del documento in
  `report/templates/`). Non è
  una fase: non entra nel grafo e non ha un manifesto. Si apre con le segnalazioni e
  prosegue con otto sezioni:
  - *segnalazioni in apertura*: i parametri che hanno il valore di OSD-734, gli avvisi
    di provenienza, l'integrità di ciò che è stato letto, il numero di tentativi
    ripetuti e di degradazioni;
  - *stato della catena*: per ogni fase, se è conclusa secondo il suo manifesto, se
    è stata calcolata sugli artefatti a monte oggi su disco, e con quale versione
    della configurazione ha girato l'ultima volta; il punto di ripresa, se dichiarato;
  - *gate e controlli di avvio*: l'esito dei quindici gate come li ha visti S0
    (`gates.json`) e, per ogni avvio, l'esito di G15, G14 e G12, che l'esecutore
    ripete a ogni avvio e ripresa e registra come evento strutturato
    (`controlli_di_avvio`), superati o no, nel log e nel registro degli avvii;
  - *configurazione*: gli avvii del registro, tutte le versioni della
    configurazione registrata con le differenze e le fasi eseguite con ciascuna, e i
    parametri in uso con l'origine di ogni valore;
  - *provenienza*: versione, impronta del sorgente, commit e immagine dichiarata di
    ogni fase, con l'avviso dove il sorgente registrato differisce da quello della
    pipeline che genera il documento, e la versione di dada2 con la correzione;
  - *decisioni prese automaticamente*, lette dagli artefatti che le registrano:
    tentativi ripetuti con l'aggiustamento applicato e degradazioni (manifesti), soglia
    di profondità effettiva per piastra e per i campioni senza piastra, con stadio,
    origine, biologici in ingresso e conservati, e sopra il modello scelto, gli AIC,
    l'R² dichiarato e la mediana (`soglia.json`, `profondita_campioni.tsv`),
    troncamento suggerito con la regola di S1 (lunghezza senza perdite, lunghezza
    ammessa per classe, perdita peggiore per singolo campione), modalità di
    decontaminazione (`decontam_riepilogo.json`), campioni riclassificati dalla regola
    `ctrl.blank_override_*` (`crosswalk.tsv`), campioni e varianti esclusi dai filtri
    finali con il motivo (`esclusioni.tsv`, `varianti_rimosse.tsv`);
  - *tracciamento delle letture*, per campione e riassunto per classe;
  - *risultato finale*: dimensioni dell'oggetto, letture per classe lungo la catena,
    impronte dei file consegnati, grafici del modello d'errore di S3.

  **L'origine dei valori.** La configurazione registrata in `00_config` elenca in
  `dichiarati` i parametri impostati nel file di ingresso, come li riconosce lo schema
  al caricamento; il report dice così, per ogni parametro, se il valore è dichiarato,
  preso dal predefinito, derivato da altri parametri, o aggiustato da un tentativo
  ripetuto. I parametri che descrivono il dataset sono obbligatori e non possono essere
  ereditati. Per gli otto con un predefinito tarato su OSD-734 (`FATTI_OSD734` in
  `config/defaults.py`, unica fonte, tenuta allineata ai marcatori `[OSD-734]`
  dell'esempio da un test) il report elenca in apertura quelli che l'esecuzione ha
  preso per difetto senza dichiararli, con il fatto accertato su OSD-734 che
  giustificava il valore e che cosa fare: controllare se quel fatto vale per il
  dataset in uso, e dichiarare il parametro con lo stesso valore o con quello adatto.
  Un parametro dichiarato non compare più, anche a una ripresa: dichiararlo registra
  una nuova versione della configurazione pur senza cambiare alcun valore. Se non ce
  ne sono, la sezione lo dice in una riga.

  **Che cosa legge e che cosa no.** Lo stato riportato è quello dei manifesti: se una
  fase è ancora valida per la configurazione e i dati di adesso lo giudica `resume`.
  Ogni artefatto letto è verificato contro il checksum del manifesto della sua fase:
  uno alterato è segnalato in apertura e non riportato. L'unico confronto con qualcosa
  fuori dalla cartella è l'avviso di provenienza, che per definizione confronta il
  sorgente registrato con quello della pipeline in uso. La versione di dada2 è quella
  dichiarata dal file di blocco `run.lockfile`, registrata all'avvio se il file
  è raggiungibile (dalla cartella di lavoro o dalla radice del repository): non è una
  misura della libreria caricata da R. Gli avvii si leggono dal solo registro
  `99_logs/avvii.jsonl`, mai dal log, che ruotando ne darebbe una storia incompleta
  senza poterlo dire. Un'esecuzione prodotta prima di queste registrazioni non ha il
  registro degli avvii, e quindi né i controlli di avvio né il file di blocco, e non ha
  l'elenco dei parametri dichiarati: il report lo dichiara, senza ricostruirli. Sul
  dataset di riferimento il report si genera in meno di un secondo.

Le quindici fasi sono tutte realizzate. Con `phylo.enabled` falso, il predefinito, S9
è disattivata e la catena arriva a `ps_final.rds` senza albero.

Questa sezione viene aggiornata a ogni avanzamento del lavoro.
