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
i test delle fasi di calcolo chiedono anche i pacchetti Bioconductor, e girano
nell'immagine. Dove mancano, quei test si saltano.

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
| 5 | tutte le fasi realizzate sono concluse, ma la prossima non esiste ancora come codice: uno stato transitorio dello sviluppo, oggi dopo S7 |

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
- il catalogo degli errori (50 codici totali): ogni codice porta un messaggio che dice
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
  intercettata da R, anche quella del codice C dei pacchetti con `R_Calloc` e quella
  del caricamento di una libreria, o il processo ucciso dal sistema con `SIGKILL`) e
  la traduce nel codice che la fase indica; gli script caricano i pacchetti con
  `richiedi_pacchetti`, che conserva il messaggio originale invece di ridurlo a
  «pacchetto non disponibile»; per rendere il riconoscimento indipendente dalla lingua
  della macchina, impone a R i messaggi in inglese. Può limitare la memoria virtuale
  del processo figlio, riducendo allora a uno i thread dell'algebra lineare: con
  OpenBLAS multithread, come nell'immagine, R sotto quel limite resta bloccato in
  uscita. Accetta un tempo massimo, allo scadere del quale uccide il processo e il suo
  gruppo: un figlio bloccato non blocca l'orchestratore;
- le funzioni R condivise in `R/lib/`: `io_json.R` per leggere e scrivere il JSON in
  modo atomico, `errors.R` con il punto d'ingresso degli script di fase e la
  dichiarazione degli errori con un codice del catalogo, `letture.R` per registrare
  per ogni passo quante letture restano a ciascun campione, `risorse.R` per riportare
  nel log la memoria di picco del processo. Le usano gli script delle fasi e gli
  script doppioni dei test, in `tests/r_doppioni/`;
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
  fase a monte la rende da rifare, insieme a tutte quelle che ne dipendono. Ogni fase
  dichiara i parametri da cui dipende, e solo quelli contano: cambiare `err.nbases`
  rifà S3 e ciò che segue, non S0, S1 e S2. La dichiarazione è obbligatoria, perché
  una fase che non la fa non si registra, ed è vincolante: il codice Python di una
  fase vede la configurazione attraverso una vista ristretta ai parametri dichiarati
  (`config/vista.py`), che solleva un errore su qualunque altro accesso, e gli script
  R ricevono soltanto i parametri che la fase passa e falliscono se ne leggono uno
  non ricevuto. Una dipendenza dimenticata diventa così un errore al primo test che
  esercita la fase, non un risultato obsoleto. L'unica eccezione è G15, che verifica
  la coerenza dell'intera configurazione: S0 lo esegue sulla configurazione completa,
  prima di vederla ristretta, e i parametri che servono solo a G15 non entrano nella
  sua impronta, perché G15 si ripete comunque a ogni avvio; cambiare una soglia di S6
  rifà S6 e ciò che segue, non S0. I parametri che non incidono sui
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
  gigabyte, costa secondi (5,6 sul dataset di riferimento con S0-S7 concluse); delle quindici fasi oggi esistono come codice le prime
  otto, da S0 a S7, e le altre risultano non realizzate;
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
  e il manifesto della fase registra il codice, il valore dichiarato e quello usato.
  Per un errore transitorio la fase può invece dichiarare esplicitamente, per quel
  codice, un nuovo tentativo senza modifiche, registrato come ogni altro
  aggiustamento: è il caso di E-S2-03, errore di lettura;
  la fase resta giudicata sulla configurazione dichiarata, quindi una ripresa non la
  rifà. Le degradazioni non fermano l'esecuzione: la fase le registra e proseguono nel
  log e nel manifesto; S0 vi registra E-S0-15, emessa da G08, e S1 vi registra E-S1-01.
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
  che stima la lunghezza minima dalle prime `qc.head_reads` letture: ricontrolla la
  condizione sul minimo vero e, se `filter.truncLen` lo supera, si ferma con
  `E-S1-02`. Registra inoltre E-S1-01, lo scarto del troncamento sotto il minimo oltre
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
  nei positivi e lo 0,03% nei negativi, contro l'8,7% delle varianti nei biologici;
- **la fase S7, il filtro di lunghezza** (`steps/s07_asv_length.py`,
  `R/07_asv_length.R`), che tiene le varianti fra `asv.len_min` e `asv.len_max` e
  scrive la tabella delle varianti accanto a S6. Con letture troncate a lunghezza fissa
  le varianti hanno tutte la stessa lunghezza, e sul dataset di riferimento il filtro
  non toglie nulla: 12.045 varianti di 137 basi. Che tolga le varianti fuori intervallo
  è verificato con varianti sintetiche;
- il tracciamento delle letture (`runner/tracciamento.py`): ogni fase registra i propri
  passi in file suoi, `letture_<passo>.tsv`, e la tabella completa si ricompone
  leggendo quelli delle fasi concluse nell'ordine del grafo, senza che una fase
  modifichi mai un file di un'altra. Oggi i passi sono uno per fase: le letture
  grezze (S1), quelle in ingresso e in uscita dal filtro (S2), quelle attribuite a una
  variante (S4), quelle nella tabella (S5), quelle senza chimere (S6) e quelle dopo il
  filtro di lunghezza (S7);
- i quattro sottocomandi della riga di comando, descritti sopra, con i codici di
  uscita documentati. `report` produce oggi un **resoconto provvisorio** dello stato,
  ricavato dai manifesti delle fasi: fasi concluse, disattivate e da eseguire,
  aggiustamenti applicati, degradazioni registrate. Non è il report definitivo, che
  non è ancora realizzato.

Delle fasi di analisi sono realizzate quelle da S1 a S7; le altre, da S8 a S14, non
sono ancora realizzate.

Questa sezione viene aggiornata a ogni avanzamento del lavoro.
