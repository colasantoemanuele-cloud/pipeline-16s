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

L'installazione rende disponibile il comando `amplicon16s`. La suite di test si esegue
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
`pyproject.toml`, quelle R in `renv.lock`. `renv.lock` è generato dall'immagine già
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

## Stato dell'implementazione

Sono realizzati:

- la struttura del repository e il packaging Python, installabile in modalità sviluppo;
- l'ambiente di esecuzione containerizzato, che contiene Python 3.11 e R 4.5 con i
  pacchetti Bioconductor della pipeline, con le versioni bloccate e verificate a ogni
  costruzione dell'immagine;
- la validazione della configurazione: lo schema copre tutti i parametri della
  pipeline e ne verifica tipo e dominio, respingendo le chiavi sconosciute;
  `config/config.example.yaml` ne è un'istanza completa;
- il gate G15, che verifica la coerenza fra parametri — le combinazioni singolarmente
  valide ma insensate messe insieme — risolve i parametri che discendono da altri e
  registra la configurazione effettivamente usata in `00_config/resolved.yaml` con un
  digest che ne identifica la combinazione. Tutto questo avviene prima che venga
  allocato qualunque calcolo — è il gate che apre la sequenza di S0, perché un errore
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
  posizione non è una superficie — aria, contenitori non aperti, posizioni non
  dichiarate — resta senza modulo, perché raggrupparlo per modulo mescolerebbe l'aria
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
  presenza della regione amplificata dichiarata — non la qualità del campione:
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
- il catalogo degli errori: ogni codice porta un messaggio che dice cosa fare e una
  categoria di gestione fra revisione umana, retry automatico, retry seguito da
  revisione, e degradazione automatica. Il retry automatico è un elenco chiuso di
  quattro codici, gli stessi dichiarati in `retry.whitelist`. Sono catalogati i codici
  di tutte le fasi, comprese quelle non ancora realizzate: il catalogo esiste, il
  codice che solleverà quegli errori no. Accanto ai codici di fase ci sono i quattro
  codici `E-R-*` del ponte verso R, tutti a revisione umana, per ciò che una fase non
  può dichiarare da sé: interprete assente, processo morto senza esito, errore R privo
  di codice, memoria esaurita in una fase che non prevede il retry;
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
  script. Riconosce la memoria esaurita in entrambe le forme — l'allocazione fallita
  intercettata da R, il processo ucciso dal sistema con `SIGKILL` — e la traduce nel
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
- la gestione degli artefatti: l'albero delle cartelle di output sotto `io.out_root`,
  una per fase, e un manifesto per fase che registra ogni artefatto con il suo
  checksum. Il manifesto permette di stabilire se una fase ha prodotto i propri
  artefatti e se sono ancora integri; gli artefatti scritti dai processi R vi si
  registrano allo stesso modo di quelli scritti da Python;
- la registrazione degli eventi su due uscite: la console per chi segue l'esecuzione e
  un file JSON Lines con rotazione sotto `99_logs`, in un formato che si interroga per
  codice, fase o categoria invece di doversi leggere;
- la catena di integrazione continua, che a ogni push installa il pacchetto con
  Python 3.11 ed esegue la suite di test. La catena installa anche R 4.5.2 e jsonlite
  2.0.0, le stesse versioni del container, così i test del ponte lanciano davvero gli
  script R; lì l'assenza di R fa fallire quei test invece di saltarli;
- i quattro sottocomandi della riga di comando — `run`, `resume`, `validate` e
  `report` — come **segnaposto non operativi**: sono invocabili e dichiarano
  l'interfaccia prevista, ma non eseguono alcuna elaborazione.

L'implementazione delle fasi di analisi non è ancora iniziata.

Questa sezione viene aggiornata a ogni avanzamento del lavoro.
