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
- **R ≥ 4.1** esegue i calcoli scientifici (dada2, phyloseq, DECIPHER, decontam) come
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
docs/                 Documentazione
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
con `pytest`. I calcoli scientifici girano nell'immagine descritta qui sotto, che porta
con sé il proprio R. I test del ponte verso R lanciano invece script R veri, e chiedono
solo `Rscript` nel PATH (o indicato con `AMPLICON16S_RSCRIPT`) e il pacchetto `jsonlite`;
dove mancano, quei test si saltano.

### Ambiente containerizzato

L'immagine contiene sia Python sia R con tutte le librerie della pipeline, ed è il modo
previsto per eseguire la pipeline in modo riproducibile:

```bash
docker build -f container/Dockerfile -t amplicon16s:dev .
docker run --rm amplicon16s:dev
```

Le versioni sono bloccate su entrambi i fronti: le dipendenze Python in
`pyproject.toml` (e nei file `requirements.txt` / `requirements-dev.txt`), quelle R in
`renv.lock`. `renv.lock` è generato dall'immagine già
costruita, quindi registra le versioni effettivamente ottenute e non quelle attese; la
build lo rilegge e fallisce se anche una sola versione non coincide. Dopo ogni modifica
ai pacchetti R va rigenerato:

```bash
scripts/genera_renv_lock.sh amplicon16s:dev
```

Le immagini di partenza sono ancorate per digest e non per tag, perché un tag può essere
riassegnato a un'immagine diversa mentre un digest no.

### Identificare l'immagine per digest

Un tag locale come `amplicon16s:dev` identifica l'immagine solo su quella macchina e può
essere riassegnato. Per riferirsi senza ambiguità a un'immagine costruita, si usa il suo
digest:

```bash
# Digest dell'immagine costruita in locale
docker image inspect amplicon16s:dev --format '{{.Id}}'

# Digest con cui l'immagine è pubblicata in un registry, disponibile dopo il push
docker image inspect amplicon16s:dev --format '{{index .RepoDigests 0}}'
```

Il primo comando restituisce l'identificatore del contenuto locale; il secondo il
digest con cui l'immagine viene recuperata da un registry, ed è quello da citare quando
si deve riprodurre un'analisi a distanza di tempo.

## Uso della riga di comando

Ogni sottocomando richiede il file di configurazione, che non viene mai modificato:

```bash
amplicon16s validate --config config.yaml   # esegue solo S0
amplicon16s run      --config config.yaml   # esegue dall'inizio, in una cartella nuova
amplicon16s resume   --config config.yaml   # riprende dagli artefatti esistenti
amplicon16s report   --config config.yaml   # resoconto provvisorio dello stato
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
controlli di avvio non registra nulla, perché nessuna fase è partita.

Codici di uscita, per chi lancia la pipeline da uno script o da uno scheduler:

| Codice | Significato |
|---|---|
| 0 | successo |
| 1 | errore imprevisto, cioè un difetto del programma; la traccia è nel log in `99_logs` |
| 2 | riga di comando non valida (argomenti mancanti o sconosciuti) |
| 3 | errore di configurazione: il file non è valido, G15 lo respinge, oppure `run` trova la cartella di output già usata; nessuna fase è partita |
| 4 | arresto con punto di ripresa dichiarato, stampato e scritto in `99_logs/punto_di_ripresa.json` e `.txt` |
| 5 | tutte le fasi realizzate sono concluse, ma la prossima non esiste ancora come codice: uno stato transitorio dello sviluppo, oggi dopo S0 |

## Stato dell'implementazione

Sono realizzati:

- la struttura del repository e il packaging Python, installabile in modalità sviluppo;
- l'ambiente di esecuzione containerizzato, che contiene Python 3.11 e R 4.5 con i
  pacchetti Bioconductor della pipeline, con le versioni bloccate e verificate a ogni
  costruzione dell'immagine;
- la validazione della configurazione: lo schema copre tutti i parametri della
  pipeline e ne verifica tipo e dominio, respingendo le chiavi sconosciute;
  `config/config.example.yaml` ne è un'istanza completa;
- il gate G15, che verifica la coerenza fra parametri (le combinazioni singolarmente
  valide ma insensate messe insieme) e risolve i parametri che discendono da altri.
  La configurazione effettivamente usata viene registrata in `00_config/resolved.yaml`
  con un digest che ne identifica la combinazione, all'avvio di ogni esecuzione e
  senza mai sovrascrivere le versioni precedenti, come descritto nella sezione sull'uso
  della riga di comando. Tutto questo avviene prima che venga
  allocato qualunque calcolo: è il gate che apre la sequenza di S0, perché un errore
  di configurazione va scoperto prima di aprire un solo file. Un parametro derivato,
  `prev.min_samples`, dipende dal numero di campioni biologici e viene calcolato
  quando i metadati sono stati letti, non prima;
- l'inventario dei campioni: il crosswalk fra i file di letture, la tabella di assay
  dell'amplicone e la tabella campioni di studio, con la classe di ciascun campione.
  La chiave del join è l'accession estratto dal nome del file, non il nome del
  campione, che può ripetersi fra repliche. Il join verso la tabella campioni di
  studio resta ristretto alle righe dell'assay, perché quella tabella è condivisa fra
  più assay dello stesso studio. Un file facoltativo associa a ogni campione la
  piastra di estrazione, la corsa di sequenziamento e il modulo; anche quel file si
  aggancia per accession, e quando manca piastra e corsa restano nulle e la pipeline
  procede in modalità a corsa singola;
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
- **la fase S0, la validazione iniziale**, con tutti e quindici i gate: dagli ingressi
  leggibili e dalle tabelle apribili fino al troncamento compatibile con le lunghezze
  osservate, all'assenza del primer e alla disponibilità delle risorse.
  L'assenza del primer è accompagnata da un controllo positivo, che verifica la
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
  manifesto con il proprio checksum;
- il catalogo degli errori (47 codici totali): ogni codice porta un messaggio che dice
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
  la dichiarazione d'esito, scritta dallo script per ultima e in modo atomico. È la
  dichiarazione a portare il codice del catalogo, che il codice di uscita di un processo
  non può rappresentare; la sua assenza distingue un processo morto senza dichiarare
  nulla da un fallimento dichiarato. Il ponte traduce l'esito in un'eccezione della
  gerarchia della pipeline con il codice dichiarato, registra nel manifesto gli
  artefatti dichiarati e nel log strutturato le uscite standard e di errore dello
  script. Riconosce la memoria esaurita in entrambe le forme (l'allocazione fallita
  intercettata da R o il processo ucciso dal sistema con `SIGKILL`) e la traduce nel
  codice che la fase indica; per rendere il riconoscimento indipendente dalla lingua
  della macchina, impone a R i messaggi in inglese. Può limitare la memoria virtuale
  del processo figlio, riducendo allora a uno i thread dell'algebra lineare: con
  OpenBLAS multithread, come nell'immagine, R sotto quel limite resta bloccato in
  uscita. Accetta un tempo massimo, allo scadere del quale uccide il processo e il suo
  gruppo: un figlio bloccato non blocca l'orchestratore;
- le funzioni R condivise in `R/lib/`: `io_json.R` per leggere e scrivere il JSON in
  modo atomico, `errors.R` con il punto d'ingresso degli script di fase e la
  dichiarazione degli errori con un codice del catalogo, `letture.R` per registrare
  per ogni passo quante letture restano a ciascun campione. Nessuno script di fase le
  usa ancora: finora le esercitano solo gli script doppioni dei test, in
  `tests/r_doppioni/`;
- la gestione degli artefatti: l'albero delle quattordici cartelle di output sotto
  `io.out_root` e, in ciascuna, un manifesto che registra ogni file scritto con il suo
  checksum; gli artefatti scritti dai processi R vi si registrano allo stesso modo di
  quelli scritti da Python. Il completamento di una fase non si legge da quel
  manifesto, che è per cartella, ma da quello proprio della fase, descritto qui sotto;
- **la classe base delle fasi, il grafo e lo stato di un'esecuzione**
  (`steps/base.py`, `runner/graph.py`, `runner/project.py`). Ogni fase eredita da
  `PipelineStep` lo stesso scheletro: verifica dei prerequisiti, calcolo, validazione
  degli artefatti, registrazione dell'esito; S0 è realizzata su questa base, con il
  comportamento di prima. Il grafo dichiara le quindici fasi da S0 a S14 in ordine,
  con la cartella in cui ciascuna scrive e le fasi di cui consuma gli artefatti; S9, la
  filogenesi, è attiva solo con `phylo.enabled`, e S12 precede obbligatoriamente S13.
  Ogni fase conclusa scrive nella propria cartella un manifesto suo, `manifest_S<n>.json`,
  così due fasi che condividono una cartella si concludono separatamente. Il manifesto
  registra anche su che cosa la fase è stata calcolata: l'impronta della
  configurazione, quella di ogni fase a monte e, per S0, quella dei dati di ingresso.
  Una fase è conclusa se il manifesto c'è, i suoi artefatti sono integri e questi
  ingressi coincidono con quelli di adesso: cambiare un parametro o ricalcolare una
  fase a monte la rende da rifare, insieme a tutte quelle che ne dipendono. Per
  difetto una fase dipende dall'intera configurazione, meno i parametri che non
  incidono sui risultati, dichiarati in un solo elenco in `config/resolve.py`:
  `run.threads`, `io.out_root`, `retry.enabled`, `retry.max_attempts`. Cambiarli, o
  spostare la cartella di un'esecuzione conclusa, non la rende da rifare; il digest
  scritto in `00_config/resolved.yaml` resta calcolato sull'intera configurazione.
  Una fase avviata prima delle fasi da cui dipende solleva `E-GRAFO-01`, oppure
  `E-S13-01` se a mancare è la decontaminazione prima del filtro di prevalenza.
  `ProjectRun` legge questo stato dal disco in una valutazione, a cui si chiede quali
  fasi sono concluse, quali disattivate, quali da eseguire e quale è la prossima;
  dentro una valutazione il checksum di ogni artefatto è calcolato una volta sola, ma
  una valutazione completa li calcola tutti, e con gli artefatti di S2 e S4, da
  gigabyte, costerà secondi; delle quindici fasi oggi esiste come codice solo S0,
  e le altre risultano non realizzate;
- **l'esecutore e la politica dei tentativi** (`runner/executor.py`,
  `runner/retry.py`). A ogni avvio, con `run` come con `resume`, l'esecutore ripete
  la verifica di coerenza della configurazione (G15) e quella delle risorse della
  macchina (G14), senza rieseguire S0: processori e spazio su disco si controllano
  anche quando S0 è già conclusa. Poi esegue in ordine le fasi da eseguire. Una fase
  fallita si comporta secondo la categoria del suo codice: a revisione umana ci si
  ferma subito; un codice ripetibile si ritenta solo se è in `retry.whitelist` e
  `retry.enabled` è vero, entro `retry.max_attempts` tentativi totali compreso il
  primo, e solo se la fase dichiara per quel codice un'azione correttiva, perché
  ritentare identico darebbe lo stesso esito. L'azione correttiva cambia un parametro
  in una copia in memoria della configurazione (oggi `run.batch_size` o `err.nbases`)
  e il manifesto della fase registra il codice, il valore dichiarato e quello usato;
  la fase resta giudicata sulla configurazione dichiarata, quindi una ripresa non la
  rifà. Le degradazioni non fermano l'esecuzione: la fase le registra e proseguono nel
  log e nel manifesto; S0 vi registra oggi E-S0-15 ed E-S1-01, emesse dai suoi gate.
  Quando l'esecuzione si ferma, l'esecutore dichiara il punto di ripresa: fase,
  codice, messaggio del catalogo, tentativi fatti e comando per ripartire. Le fasi da
  S2 a S5, i cui codici sono ripetibili, non esistono ancora: il meccanismo è provato
  con fasi doppione;
- la registrazione degli eventi su due uscite: la console per chi segue l'esecuzione e
  un file JSON Lines con rotazione sotto `99_logs`, in un formato che si interroga per
  codice, fase o categoria invece di doversi leggere;
- la catena di integrazione continua, che a ogni push installa il pacchetto con
  Python 3.11 ed esegue la suite di test. La catena installa anche R 4.5.2 e jsonlite
  2.0.0, le stesse versioni del container, così i test del ponte lanciano davvero gli
  script R; lì l'assenza di R fa fallire quei test invece di saltarli;
- i quattro sottocomandi della riga di comando, descritti sopra, con i codici di
  uscita documentati. `report` produce oggi un **resoconto provvisorio** dello stato,
  ricavato dai manifesti delle fasi: fasi concluse, disattivate e da eseguire,
  aggiustamenti applicati, degradazioni registrate. Non è il report definitivo, che
  non è ancora realizzato.

L'implementazione delle fasi di analisi non è ancora iniziata.

Questa sezione viene aggiornata a ogni avanzamento del lavoro.
