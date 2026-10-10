# Guida: configurare ed eseguire la pipeline su un nuovo dataset

Questa guida accompagna chi ha un dataset 16S single-end proprio, dai file grezzi
all'oggetto finale e alle analisi ecologiche di base. Il README dice che cosa fa lo
strumento; qui c'è come lo si usa, passo per passo, e che cosa significa ogni parametro
da decidere.

Nei comandi `<immagine>` è l'immagine pubblicata, indicata nel README. I comandi si
lanciano dalla radice di un clone del repository.

## 1. Preparare i dati

### 1.1 Le letture

- Un file FASTQ per campione, compresso con gzip o no (`.fastq.gz`, `.fq.gz`, `.fastq`,
  `.fq`), tutti nella stessa cartella, senza sottocartelle. Conviene che abbiano tutti
  la stessa estensione.
- Il nome di ogni file deve contenere una parte unica per campione (un accession, un
  codice di libreria): è la chiave con cui file e metadati si agganciano.
- Un file senza alcuna lettura è ammesso solo per un controllo negativo.

### 1.2 Letture appaiate

La pipeline tratta dati single-end. Di una corsa appaiata si usano le sole letture
forward: un file per campione con la prima lettura di ogni coppia, in una cartella a
parte. I file delle seconde letture (`_R2`, `_2`) non vanno nella cartella. Se un
archivio distribuisce le due letture di ogni coppia nello stesso file, vanno separate
prima: `dati/osd276/scarica_letture.py` mostra come, per un dataset pubblico.

`init` e la validazione riconoscono le letture appaiate dal nome dei file e dal
contenuto, e si fermano senza proporre né eseguire nulla.

### 1.3 I metadati

Tabelle di testo UTF-8 separate da tabulazioni, con una riga di intestazione. Due forme:

- **una tabella dei campioni**: una riga per campione, con almeno una colonna da cui si
  ricava la chiave del file (di solito il nome del file), una con il nome del campione
  e, se ci sono controlli, una con il tipo del campione;
- **ISA-Tab**, come negli archivi OSDR di NASA GeneLab: una tabella di assay (`a_*.txt`,
  un campione per riga con il file delle letture) e una di studio (`s_*.txt`, con tipo
  e variabili dei campioni). La tabella di studio può contenere campioni di altri
  assay: la pipeline usa le sole righe dei campioni dell'assay.

I valori non devono contenere tabulazioni né a capo.

### 1.4 Il file del lotto (facoltativo)

Una tabella con una riga per campione che dichiara la piastra di estrazione e la corsa
di sequenziamento, e per i controlli positivi in diluizione le cellule seminate. Deve
identificare il campione con la stessa chiave dei file. Senza, la pipeline stima un
solo modello d'errore e non distingue le piastre.

### 1.5 Il riferimento tassonomico

Un file nel formato dei training set di dada2 (`assignTaxonomy`). Quello usato per i
dataset del repository si scarica con `python3 dati/riferimento/scarica_riferimento.py`.

## 2. Scrivere la configurazione con `init`

`init` guarda letture e metadati e scrive un file di configurazione di partenza. Non
scarica nulla, non modifica i file di ingresso e non scrive su un file che esiste già.

Con i dati fuori dal repository, montati nel container. La cartella di lavoro deve
esistere ed essere scrivibile da chi lancia il comando (`mkdir -p /percorso/di/lavoro`):
se manca, Docker la crea di un altro utente e `init` si rifiuta di scrivervi.

```bash
docker run --rm -u "$(id -u):$(id -g)" \
  -v "$(pwd)":/app -v /percorso/dei/dati:/dati:ro -v /percorso/di/lavoro:/lavoro -w /app <immagine> \
  python3 scripts/esegui.py init \
    --reads /dati/letture --metadata /dati/campioni.tsv \
    --reference /dati/riferimento.fa.gz \
    --out-root /lavoro/esecuzione --output /lavoro/configurazione.yaml
```

| Argomento | Significato |
|---|---|
| `--reads CARTELLA` | la cartella delle letture (obbligatorio) |
| `--metadata FILE_O_CARTELLA` | la tabella dei campioni, oppure la cartella ISA-Tab (obbligatorio) |
| `--assay NOME` | in una cartella ISA-Tab con più tabelle di assay, il nome di quella del 16S |
| `--study-table FILE` | la tabella di studio, se `--metadata` è una tabella di assay |
| `--batch-table FILE` | il file del lotto |
| `--reference FILE` | il riferimento tassonomico |
| `--primer SEQUENZA` | il primer forward usato, in codici IUPAC, se `init` non lo riconosce da sé |
| `--out-root CARTELLA` | la cartella di uscita prevista per l'esecuzione |
| `--output FILE` | il file di configurazione da scrivere (obbligatorio; non deve esistere) |
| `--head-reads N` | ispeziona le sole prime N letture di ogni file; senza, tutte |
| `--skip-profile` | non calcola il profilo di qualità: nessuna proposta di troncamento |
| `--threads N` | processi per l'ispezione; senza, quelli utilizzabili |

I percorsi scritti nella configurazione sono quelli visti dal container (`/dati/...`,
`/lavoro/...`, `/app/...` per ciò che sta nel clone): ogni comando successivo va
lanciato con gli stessi `-v`. Dati in più cartelle si montano con un `-v` ciascuna.

`init` legge tutte le letture per riconoscere le coppie e cercare il primer, e calcola
il profilo di qualità con lo stesso calcolo della fase S1: su 960 campioni e 21 milioni
di letture impiega alcuni minuti. Esce con 0 se ha scritto il file, con 3 se si rifiuta
(file già esistente, letture appaiate, ingressi non trovati): in quel caso non scrive
nulla.

### 2.1 Le tre specie di valore

Ogni riga del file dice da dove viene il suo valore, nel commento che la precede:

| Marcatore | Che cos'è | Come è scritta |
|---|---|---|
| `[RILEVATO]` | un fatto verificato su tutti i file | riga attiva |
| `[INDICATO]` | un percorso dato sulla riga di comando | riga attiva |
| `[PROPOSTO]` | l'esito di un'euristica, che può sbagliare | riga in commento, con il valore già scritto |
| `[DA SCEGLIERE]` | una scelta di chi conduce l'analisi | riga in commento, con un valore di partenza e il suo motivo, oppure senza valore |

I parametri che hanno un predefinito (`[STANDARD]`, `[UNICA FORMA]`, `[STRUMENTO]`,
`[FACOLTATIVO]`) sono in commento con quello: si possono lasciare così.

Nessun valore entra nella configurazione senza un gesto: una proposta si accetta
togliendo il commento alla sua riga (`  # column: "tipo"` diventa `  column: "tipo"`),
una scelta si scrive. Le spiegazioni cominciano con `##`, le righe dei parametri in
commento con un solo `#`. In testa al file, e sulla console, c'è quanti valori sono
rilevati, proposti e da scegliere, l'elenco di questi ultimi, e ciò che `init` dichiara
di non poter stabilire.

### 2.2 Che cosa viene rilevato

- Il modello dei nomi dei file (`io.fastq_glob`), se hanno tutti la stessa estensione.
- La colonna e l'espressione della chiave del campione (`meta.accession_column`,
  `io.accession_regex`): la colonna i cui valori corrispondono uno a uno ai file, e
  un'espressione che estrae da ogni file una chiave unica e la stessa chiave dalla
  colonna. La verifica è fatta su tutti i file e tutte le righe. L'espressione può
  essere scritta diversamente da come la scriverebbe una persona: conta che estragga la
  stessa chiave. Se nessuna colonna corrisponde, o più d'una, `init` non sceglie: elenca
  le candidate e lascia la riga da scegliere.
- La colonna del nome del campione (`meta.sample_id_column`), con due tabelle: quella
  che aggancia ogni riga della tabella di assay a una sola riga della tabella di studio.
- La colonna della chiave nel file del lotto (`meta.batch_key_column`). Senza file del
  lotto, i parametri che lo riguardano sono rilevati nulli o vuoti.
- L'MD5 del riferimento (`tax.ref_md5`).
- Il primer, se le letture di ogni file iniziano con uno dei primer 16S comuni
  dell'elenco (515F nelle due forme, 341F, 357F, 27F, 8F, 799F): `qc.primer_sequence` e
  `filter.trimLeft` pari alla sua lunghezza. Se una parte dei file non lo porta, o con
  `--head-reads`, primer e taglio sono una proposta, e i file che non lo portano sono
  nominati: un taglio fisso toglierebbe loro basi vere.

Una chiave del campione di sole cifre, o più corta di tre caratteri, può corrispondere
ai file per coincidenza: in quel caso colonna ed espressione sono una proposta, e il
commento mostra alcuni abbinamenti file-chiave da controllare.

### 2.3 Che cosa viene proposto

- La colonna del tipo di campione e le etichette (`ctrl.column`, `ctrl.blank_values`,
  `ctrl.positive_values`, `ctrl.biological_values`), dalle parole dei valori: blank,
  negative, NTC, no template per i controlli negativi; positive, mock per i positivi.
  Un valore con la sola parola «control» non è attribuito: va assegnato a mano. Senza
  controlli e senza una colonna del tipo la proposta è `ctrl.column: null`.
- Nel file del lotto, le colonne di piastra, corsa e cellule dei controlli positivi,
  dal loro nome.
- Il troncamento (`filter.truncLen`), dal profilo di qualità: la lunghezza più lunga
  che non scarta alcuna lettura e che precede il calo della qualità mediana. È lo
  stesso valore che il report suggerisce dopo la prima esecuzione.
- `err.error_function: loess_monotono`, se le qualità sono raggruppate in pochi valori.
- Nome e versione del riferimento, se si leggono dal nome del file.
- Se nessun primer dell'elenco è nelle letture: `filter.trimLeft: 0` (un primer fuori
  elenco va indicato a mano) e `qc.primer_sequence` dai metadati, se riportano una
  sequenza marcata come forward, o una sola sequenza (e allora il commento dice che
  l'indicazione manca).
- Il motivo conservato (`qc.conserved_motif`), solo se il primer è noto: il consenso
  delle basi con cui inizia la regione amplificata.

### 2.4 Che cosa non viene stabilito

`init` lo dice nel riepilogo: il taxon dei controlli positivi, la soglia di
decontaminazione, il motivo conservato se il primer non è noto (nullo è ammesso), la
riclassificazione di campioni in controlli negativi, e ogni scelta legata alla biomassa
o al disegno dello studio (soglie di qualità, prevalenza, taxa da escludere).

### 2.5 Completare il file

1. Rivedere le righe `[RILEVATO]`: sono già attive.
2. Per ogni riga `[PROPOSTO]`: controllare il valore sui propri dati, poi togliere il
   commento, o correggerlo.
3. Per ogni riga `[DA SCEGLIERE]`: scrivere il valore. La sezione 6 dice, parametro per
   parametro, che cosa significa, come si ricava e che cosa guardare per confermarlo.
4. Lanciare `validate` (sezione 3): finché manca qualcosa elenca i parametri ancora da
   dichiarare (`E-G15-10`). È la lista delle cose da fare.

Un parametro non pertinente si dichiara nullo (`null`) o vuoto (`[]`), non si omette.
Gli elenchi si scrivono in linea (`["a", "b"]`), le stringhe fra virgolette.

Una proposta è ciò che i dati suggeriscono, non ciò che uno studio ha scelto: sui due
dataset del repository alcune proposte differiscono dalla configurazione pubblicata.
Per OSD-276 `init` propone `filter.truncLen: 151`, la lunghezza di tutte le letture,
mentre la configurazione pubblicata tronca a 150: i due valori danno oggetti finali
diversi (990 e 1.005 varianti). Chi vuole riprodurre un risultato pubblicato usa i
valori pubblicati.

## 3. Eseguire

Tutti i comandi hanno la stessa forma; cambiano il sottocomando e i suoi argomenti.

```bash
docker run --rm -u "$(id -u):$(id -g)" \
  -v "$(pwd)":/app -v /percorso/dei/dati:/dati:ro -v /percorso/di/lavoro:/lavoro -w /app <immagine> \
  python3 scripts/esegui.py validate --config /lavoro/configurazione.yaml
```

| Passo | Sottocomando | Che cosa fa |
|---|---|---|
| 1 | `validate --config FILE` | controlla configurazione, metadati e file; scrive l'esito nella cartella di uscita (`00_config/`, `01_input_validation/`, `99_logs/`) |
| 2 | `resume --config FILE` | esegue le fasi che mancano: dopo un `validate` riuscito, tutte le successive |
| 3 | `status --out CARTELLA` | dice a che punto è l'esecuzione, anche mentre è in corso |
| 4 | `report --config FILE` | scrive `report/report.html` dall'esecuzione |
| 5 | `resume --config FILE` | dopo un arresto o una modifica alla configurazione, rifà solo ciò che serve |

`run --config FILE` fa in un solo comando i passi 1 e 2, ma solo in una cartella di
uscita nuova o vuota: dopo un `validate` nella stessa cartella si rifiuta e indica
`resume`.

- `validate` esce con 0 se tutto è a posto; con 3 se la configurazione non è valida
  (nessun file è stato aperto); con 4 se un controllo sui dati fallisce, e dice quale,
  con il codice e che cosa fare.
- `run` non sovrascrive: se la cartella di uscita non è vuota si rifiuta.
- `status` legge soltanto: si può lanciare da un altro terminale durante `run`.
- Un arresto (codice di uscita 4) dichiara la fase, il codice, che cosa fare e il
  comando per ripartire; corretta la causa, `resume` riparte da lì.
- La prima esecuzione serve anche a fissare le soglie: il report riporta le misure
  accanto a ciascuna (sezione 6). Cambiare una soglia e lanciare `resume` rifà le sole
  fasi che ne dipendono.

## 4. Analisi ecologiche di base

Lavorano sull'oggetto finale (`12_final/ps_final.rds`), che non modificano, con una
configurazione propria.

```bash
docker run --rm -u "$(id -u):$(id -g)" \
  -v "$(pwd)":/app -v /percorso/di/lavoro:/lavoro -w /app <immagine> \
  python3 scripts/esegui.py eco init \
    --object /lavoro/esecuzione/12_final/ps_final.rds --output /lavoro/eco.yaml
```

`eco init` riporta ciò che l'oggetto contiene (ranghi tassonomici con la quota di
varianti assegnate, albero e sua radice) e scrive una configurazione di partenza in cui
ogni riga è in commento:

- **proposti**: `design.variable`, solo se una sola colonna ha almeno due gruppi di
  almeno due campioni (con più candidate la riga resta da scegliere e le elenca: quale
  variabile lo studio confronta lo sa chi lo conduce; senza candidate i test fra gruppi
  non si potranno eseguire); `design.technical_variables`, piastra e corsa se hanno più
  di un valore, dalla più grossolana alla più annidata; `comp.rank`, il rango più fine con almeno
  metà delle varianti assegnate; `beta.unifrac_root: existing` se l'albero è radicato;
- **da scegliere**: distanze (`beta.distances`), metodi e distanze dell'ordinazione
  (`ord.methods`, `ord.distances`), seme (`run.seed`), numero di taxa nei grafici
  (`comp.top_n`), dimensione minima dei gruppi (`stat.min_group_size`), livello di
  significatività (`stat.significance_level`), pseudoconteggio della distanza di
  Aitchison (`beta.clr_pseudocount`, solo con quella distanza). Hanno un valore di
  partenza solo `ord.nmds_trymax` (20, `vegan::metaMDS`) e `ord.pcoa_correction`
  (`none`, `ape::pcoa`).

I nomi delle colonne sono quelli della tabella dei campioni dell'oggetto, non quelli
dei metadati di partenza: la posizione di prelievo è `posizione`, e una colonna portata
con `out.study_columns` prende un nome in minuscolo con i trattini bassi (`Factor
Value[Spaceflight]` diventa `factor_value_spaceflight`; la corrispondenza è in
`10_phyloseq/colonne_metadati.tsv`). `eco init` elenca i nomi che trova.

Il significato di ogni parametro è in `config/eco.example.yaml`. Completate le scelte:

```bash
docker run --rm -u "$(id -u):$(id -g)" \
  -v "$(pwd)":/app -v /percorso/di/lavoro:/lavoro -w /app <immagine> \
  python3 scripts/esegui.py eco validate \
    --object /lavoro/esecuzione/12_final/ps_final.rds --config /lavoro/eco.yaml

docker run --rm -u "$(id -u):$(id -g)" \
  -v "$(pwd)":/app -v /percorso/di/lavoro:/lavoro -w /app <immagine> \
  python3 scripts/esegui.py eco run \
    --object /lavoro/esecuzione/12_final/ps_final.rds --config /lavoro/eco.yaml \
    --out /lavoro/analisi
```

La cartella di uscita deve essere assente o vuota e fuori dalla cartella dell'oggetto.
Contiene tabelle (alfa diversità, composizione, distanze, ordinazione, test), grafici
PNG, `riepilogo.json` con gli avvisi e `manifest.json` con dimensione e SHA-256 di ogni
altro file. A parità di oggetto, configurazione e seme i file sono identici byte per byte.
Gli avvisi (`E-ECO-`) non fermano l'analisi: dichiarano una condizione che chi legge i
risultati deve conoscere, per esempio un gruppo troppo piccolo per i test o una
variabile che sembra continua.

## 5. Integrazione con altri programmi

Tre comandi scrivono sullo standard output un solo documento JSON, per la pipeline e,
con `eco` davanti, per le analisi ecologiche. Ogni documento ha `versione_formato`
(oggi 1), `comando` e `strumento`. La versione cambia solo se un campo esistente cambia
significato o sparisce: chi legge deve ignorare i campi che non conosce. Con `--json`
il testo per l'operatore va sullo standard error.

### 5.1 `schema`

```bash
docker run --rm -v "$(pwd)":/app -w /app <immagine> python3 scripts/esegui.py schema
docker run --rm -v "$(pwd)":/app -w /app <immagine> python3 scripts/esegui.py eco schema
```

| Campo | Contenuto |
|---|---|
| `gruppi` | per ogni gruppo di parametri, una riga di descrizione |
| `parametri` | un elemento per parametro (vedi sotto) |
| `json_schema` | lo schema JSON del modello di configurazione |

Ogni elemento di `parametri`: `nome` (`gruppo.chiave`), `gruppo`, `chiave`, `tipo`
(`stringa`, `intero`, `reale`, `logico`, `percorso`, `scelta`, `mappa`, `elenco di ...`,
oppure due tipi uniti da «oppure»), `obbligatorio` (senza predefinito: va scritto nel
file), `ammette_nullo`, `valori_ammessi` (per un insieme chiuso, o per gli elementi di
un elenco; altrimenti nullo), `specie` e `descrizione`; `predefinito` c'è solo se il
parametro non è obbligatorio, `vincoli` (minimo, massimo, lunghezza minima, formato)
solo se ne ha. La specie è `obbligatorio`, `percorso`,
`standard` (con `fonte`), `unica_forma`, `strumento` o `facoltativo`; per le analisi
ecologiche anche `condizionale`, con `richiesto_da` (il parametro e i valori che lo
rendono obbligatorio). Le descrizioni non entrano nel digest della configurazione.

### 5.2 `validate --json`

Nel container, con la stessa forma di `docker run` della sezione 3:

```bash
python3 scripts/esegui.py validate --config configurazione.yaml --json
python3 scripts/esegui.py eco validate --object ps_final.rds --config eco.yaml --json
```

Il controllo è lo stesso di `validate`, e così il codice di uscita.

| Campo | Contenuto |
|---|---|
| `valida` | vero se non ci sono problemi |
| `codice_uscita` | il codice di uscita del comando |
| `problemi` | un elemento per problema |
| `avvisi` | stessa forma: condizioni dichiarate che non fermano |

Ogni problema: `parametro` (il primo parametro coinvolto, o nullo se il problema non
riguarda un parametro), `parametri` (tutti), `codice` (del catalogo), `categoria`,
`messaggio`, `azione` (che cosa fare) e, per la pipeline, `gate` (il controllo che l'ha
trovato). Un parametro obbligatorio non dichiarato è un problema a sé, con `E-G15-10`.

### 5.3 `status --json`

Nel container, come sopra:

```bash
python3 scripts/esegui.py status --out CARTELLA --json
python3 scripts/esegui.py eco status --out CARTELLA --json
```

Legge la cartella di uscita senza scrivervi e senza ricalcolare i checksum: si può
lanciare durante un'esecuzione, anche da un altro container.

| Campo | Contenuto |
|---|---|
| `stato` | `non_avviata`, `in_corso`, `interrotta`, `ferma`, `conclusa`, `incompleta` |
| `fasi` | per ogni fase: `passo`, `descrizione`, `cartella`, `stato` (`conclusa`, `in_corso`, `da_eseguire`, `disattivata`, `sconosciuta`) e, se conclusa, `artefatti`, `versione`, `degradazioni` |
| `fasi_concluse` | i passi conclusi, nell'ordine di esecuzione |
| `fase_in_corso` | il passo in esecuzione, o nullo |
| `esecuzione` | se un'esecuzione ha lasciato il suo stato: `passo`, `avviata`, `ultimo_battito`, `secondi_dall_ultimo_battito`, `battito_recente` |
| `arresto` | se c'è un punto di ripresa: `passo`, `codice`, `sintesi`, `dettaglio`, `azione`, `comando` |
| `configurazione` | file e digest dell'ultima configurazione registrata |

- `in_corso`: l'esecuzione aggiorna il suo stato ogni 20 secondi. Se l'ultimo
  aggiornamento è più vecchio di un minuto lo stato è `interrotta`: il processo è stato
  ucciso senza poter dichiarare nulla, e si riprende con `resume`.
- `ferma`: l'esecuzione si è fermata con un codice; `arresto` dice che cosa fare.
- `incompleta`: alcune fasi sono concluse e nessun processo è attivo, per esempio dopo
  `validate`, che esegue la sola prima fase.
- Una fase risulta conclusa se il suo manifesto esiste: se i suoi artefatti siano
  ancora integri e validi per una configurazione lo stabilisce `resume`.

Per le analisi ecologiche lo stato è `non_avviata`, `in_corso`, `interrotta` o
`conclusa` (con `file`, il numero dei file prodotti compreso il manifesto, `avvisi`,
`configurazione`, `oggetto_di_partenza`). Un rifiuto
non lascia nulla nella cartella: lo riporta il comando stesso, e `eco validate --json`
lo anticipa.

## 6. I parametri obbligatori, uno per uno

Per ogni parametro: che cosa significa; come si ricava dai propri dati; il valore di
partenza che `init` scrive, con la sua fonte; che cosa guardare nel report della prima
esecuzione per confermarlo. «Nessuno» nella colonna della partenza vuol dire che `init`
lascia la riga senza valore: non c'è un valore difendibile senza conoscere lo studio.
Le sezioni del report sono quelle descritte nel README.

Le configurazioni dei due dataset del repository sono esempi commentati di ogni scelta:
`dati/osd734/config_osd734.yaml` (960 campioni, controlli negativi e positivi, file del
lotto) e `dati/osd276/config_osd276.yaml` (15 campioni, nessun controllo).

### 6.1 File e campioni

| Parametro | Significato | Come si ricava | Partenza (fonte) | Conferma |
|---|---|---|---|---|
| `io.fastq_glob` | modello dei nomi dei file delle letture | dall'estensione dei file | rilevato | S0: tutti i file attesi nell'inventario |
| `io.accession_regex` | espressione che estrae la chiave del campione dal nome del file | la parte del nome unica per campione | rilevata, se una sola colonna corrisponde ai file | `01_input_validation/crosswalk.tsv`: una riga per campione |
| `meta.accession_column` | colonna da cui si ricava la stessa chiave | di solito quella con il nome del file | rilevata | come sopra |
| `meta.sample_id_column` | colonna con il nome del campione | la colonna che aggancia le due tabelle | rilevata con due tabelle, proposta con una | `crosswalk.tsv`, colonna `sample_name` |
| `ctrl.column` | colonna con il tipo del campione; nullo: nessun controllo | la colonna che distingue campioni e controlli | proposta (parole dei valori, o nome della colonna) | report, «Gate di validazione»: conteggi per classe |
| `ctrl.blank_values` | valori dei controlli negativi; `[]` se non ce ne sono | copiati dalla colonna | proposti (blank, negative, NTC, no template) | conteggi per classe |
| `ctrl.positive_values` | valori dei controlli positivi; `[]` se non ce ne sono | copiati dalla colonna | proposti (positive, mock) | conteggi per classe |
| `ctrl.biological_values` | valori dei campioni biologici | copiati dalla colonna | proposti: i valori senza parole di controllo | conteggi per classe |
| `ctrl.blank_override_column` | colonna con cui trattare come controlli negativi campioni dichiarati di un altro tipo; nullo: nessuna | solo conoscendo lo studio (materiale mai usato, tamponi non aperti) | nullo (assenza) | `crosswalk.tsv`: classe e materiale |
| `ctrl.blank_override_values` | i valori di quella colonna che riclassificano; `[]` | come sopra | `[]` (assenza) | come sopra |
| `meta.module_column` | colonna con la posizione di prelievo; nullo | la colonna della posizione, se lo studio ne ha una | nullo (assenza) | `crosswalk.tsv`, colonna `posizione` |
| `meta.derive_module` | se derivare il gruppo di posizioni dalla posizione | vero se le posizioni hanno un prefisso comune per gruppo | falso (assenza) | `crosswalk.tsv`, colonna `modulo` |
| `meta.module_regex` | espressione con un gruppo di cattura che estrae il gruppo dalla posizione; nullo | guardando i valori della posizione | nullo (assenza) | come sopra |
| `meta.non_surface_positions` | posizioni che non appartengono a un gruppo (aria, materiale non usato); `[]` | dai valori della posizione | `[]` (assenza) | come sopra |
| `out.study_columns` | colonne della tabella dei campioni da portare nell'oggetto; `[]` | le variabili che le analisi a valle useranno | `[]`; il commento elenca le colonne | `10_phyloseq/colonne_metadati.tsv` |

### 6.2 File del lotto

Senza file del lotto questi parametri sono rilevati nulli o vuoti.

| Parametro | Significato | Come si ricava | Partenza (fonte) | Conferma |
|---|---|---|---|---|
| `meta.batch_key_column` | colonna del file del lotto con la chiave del campione | la colonna con gli stessi valori estratti dai file | rilevata | S0, gate G08: ogni campione ha una riga |
| `decontam.batch_column` | colonna della piastra di estrazione | dal file del lotto | proposta (plate, piastra; prima extraction) | report: campioni per piastra |
| `err.batch_column` | colonna della corsa di sequenziamento: un modello d'errore per corsa | dal file del lotto | proposta (run, corsa, lane, flowcell) | `04_error_models/`: un modello per corsa |
| `katharoseq.cell_count_column` | colonna con le cellule dei controlli positivi; nullo | dal file del lotto o dalla tabella dei campioni | proposta (cell); nullo senza | report, «Decisioni»: curve di calibrazione |
| `meta.batch_module_column` | colonna del file del lotto che dichiara il gruppo di posizioni; nullo | se il file la porta | nullo (assenza) | `crosswalk.tsv`, colonna `modulo` |
| `out.batch_columns` | colonne del file del lotto da portare nell'oggetto; `[]` | le variabili tecniche utili a valle | `[]`; il commento elenca le colonne | `10_phyloseq/colonne_metadati.tsv` |

### 6.3 Letture, primer e riferimento

| Parametro | Significato | Come si ricava | Partenza (fonte) | Conferma |
|---|---|---|---|---|
| `filter.truncLen` | lunghezza a cui si troncano le letture; le più corte sono scartate | dal profilo di lunghezze e qualità | proposta: la più lunga che non scarta letture e precede il calo di qualità (profilo di S1) | report, «Troncamento suggerito dalle letture» |
| `filter.truncLen_shortfall_warn` | di quante basi il troncamento può stare sotto la lettura più corta prima di un avviso | dalla variabilità delle lunghezze | lo scarto fra lettura più corta e troncamento proposto (dai dati) | segnalazione `E-S1-01` |
| `filter.trimLeft` | basi tolte in testa: la lunghezza del primer se è nelle letture, altrimenti 0 | dal primer | rilevato se un primer dell'elenco è nelle letture; altrimenti proposto 0 | S0, gate G10 |
| `qc.primer_sequence` | primer forward, in codici IUPAC | dal protocollo di amplificazione (vedi sotto) | rilevato nelle letture; proposto dai metadati; indicato con `--primer`; altrimenti nessuno | `01_input_validation/letture_ispezionate.tsv`, `frazione_primer` |
| `qc.conserved_motif` | espressione del motivo con cui inizia la regione amplificata; nullo | dal consenso delle prime basi (vedi sotto) | proposto se il primer è noto; altrimenti nullo | `letture_ispezionate.tsv`, `frazione_motivo` |
| `qc.max_primer_hit_frac` | frazione di letture con il primer in testa oltre la quale è considerato presente | bassa: il primer o c'è o non c'è | nessuno | `frazione_primer` per file |
| `qc.min_motif_frac` | frazione minima di letture con il motivo, sulla mediana dei campioni | sotto il valore osservato nei campioni buoni | nessuno | `frazione_motivo` per file |
| `qc.head_reads` | letture ispezionate per file dalla validazione | abbastanza da rappresentare il file | nessuno; il commento riporta le letture del file mediano | durata della validazione |
| `qc.max_frac_short_reads` | frazione massima di letture più corte del troncamento, per classe | zero se il troncamento non scarta letture | 0 con il troncamento proposto (dai dati) | report, «Troncamento suggerito» |
| `tax.ref_name`, `tax.ref_version` | nome e versione del riferimento, riportati con i risultati | dalla fonte del riferimento | proposti dal nome del file, se vi compaiono | report, «Configurazione» |
| `tax.try_rc` | se provare anche il complemento inverso | vero se le letture possono essere nell'altro orientamento | falso (`dada2::assignTaxonomy`) | report: frazione di varianti con il phylum |

**Il primer.** Il nome di un primer nei metadati non basta: lo stesso nome ha più
forme. Le sequenze che `init` riconosce sono: 515F `GTGYCAGCMGCCGCGGTAA` (Parada 2016) e
`GTGCCAGCMGCCGCGGTAA` (Caporaso 2011); 341F `CCTACGGGNGGCWGCAG`; 357F
`CCTACGGGAGGCAGCAG`; 27F `AGAGTTTGATCMTGGCTCAG`; 8F `AGAGTTTGATCCTGGCTCAG`; 799F
`AACMGGATTAGATACCCKG`. Se i metadati nominano uno di questi `init` riporta le sue forme
nel commento. Per un primer diverso, o per farsi proporre il motivo quando il primer
non è nelle letture, si rilancia `init` con `--primer SEQUENZA` e un altro `--output`.

**Il motivo conservato.** È un'espressione regolare sulle prime basi della regione
amplificata, quelle che seguono il primer: una base dove quasi tutte le letture
concordano, una classe (`[AG]`) dove ne compaiono due, un punto dove variano. `init` lo
propone dal consenso delle letture quando conosce il primer. Senza una proposta la
scelta ragionevole è `null`: la validazione verifica allora la sola assenza del primer.

### 6.4 Scelte di analisi

| Parametro | Significato | Come si ricava | Partenza (fonte) | Conferma |
|---|---|---|---|---|
| `err.randomize` | se scegliere a caso i campioni per il modello d'errore | vero se l'ordine dei file segue piastre o gruppi | falso (`dada2::learnErrors`) | `04_error_models/`: campioni usati |
| `dada.pool` | come i campioni condividono l'informazione: `false`, `true`, `pseudo` | `pseudo` se contano le varianti rare presenti in più campioni | falso (`dada2::dada`) | report: varianti e letture dopo l'inferenza |
| `chimera.min_fold_parent_over_abundance` | quanto un genitore deve superare una variante perché sia chimerica | il valore del metodo, salvo ragioni | 1,5 (`dada2::isBimeraDenovoTable`) | report: frazione di letture chimeriche |
| `filt.remove_na_phylum` | se togliere le varianti senza phylum | vero se ciò che il riferimento non riconosce è fuori bersaglio | nessuno | report, «Risultato finale»: varianti rimosse |
| `filt.exclude_taxa` | taxa da togliere, come li scrive il riferimento; `[]` | ciò che il gene cattura senza essere l'oggetto dello studio | nessuno | `12_final/intermedi/varianti_rimosse.tsv` |
| `phylo.enabled` | se costruire l'albero delle varianti finali | vero se servono distanze filogenetiche | falso (fase facoltativa) | `09_phylogeny/` |
| `phylo.max_seqs` | varianti oltre le quali la filogenesi si ferma | dal tempo e dalla memoria disponibili; senza effetto con `phylo.enabled` falso | nessuno | report: varianti finali |
| `prev.apply` | se filtrare le varianti presenti in pochi campioni | secondo lo studio | nessuno | report: varianti rimosse per prevalenza |
| `prev.min_fraction` | frazione minima di campioni in cui una variante deve comparire | più bassa con campioni eterogenei | nessuno | come sopra |
| `prev.min_count` | letture minime perché una variante conti come presente | secondo la profondità | nessuno | come sopra |

### 6.5 Controlli positivi e negativi

Senza controlli positivi i parametri della calibrazione non hanno effetto, ma vanno
dichiarati: qualunque valore valido. Per questo la configurazione di OSD-276, che non
ha controlli, dichiara `qc.min_reads_mode: katharoseq_if_available` e dà lo stesso
oggetto finale che darebbe con `none`: nel primo caso l'assenza della soglia è
riportata come avviso (`E-S11-05`), nel secondo è una scelta dichiarata. Allo stesso
modo `meta.derive_module: true` senza un'espressione non deriva nulla.

| Parametro | Significato | Come si ricava | Partenza (fonte) | Conferma |
|---|---|---|---|---|
| `ctrl.min_positives` | controlli positivi minimi per una curva e per giudicare un livello | dal numero di controlli per piastra e per livello | 3: un punto più dei due parametri della curva | report, «Decisioni»: curve per piastra |
| `ctrl.min_positive_pass_frac` | frazione minima di controlli positivi conformi | secondo quanto ci si fida dei controlli | nessuno | report: controlli non conformi |
| `ctrl.positive_gate` | vero: sotto la frazione l'esecuzione si ferma; falso: avvisa | falso finché la conformità non è nota | falso | segnalazioni `E-S11-03`, `E-S11-04` |
| `katharoseq.target_taxon` | taxon atteso nei controlli positivi, come lo scrive il riferimento; nullo | da chi ha allestito i controlli | nessuno con controlli positivi; nullo senza | `11_controls/`: fedeltà dei controlli |
| `katharoseq.collapse_rank` | rango a cui si sommano le varianti del taxon atteso | il più fine a cui il riferimento lo distingue | `Genus`: il più fine che la pipeline assegna | come sopra |
| `katharoseq.target_sensitivity` | fedeltà a cui si legge la profondità minima dei campioni | quanta contaminazione si accetta nei campioni poveri | nessuno | report: soglia per piastra, campioni esclusi |
| `katharoseq.min_r2` | bontà minima di una curva valida | secondo la qualità dei controlli | nessuno | report: R² delle curve |
| `qc.min_reads_mode` | da dove viene la profondità minima: `katharoseq_if_available` o `none` | `none` senza controlli positivi in diluizione | secondo i controlli della proposta | report: campioni esclusi per profondità |
| `decontam.threshold` | soglia sotto la quale una variante è un contaminante | più alta con bassa biomassa; il predefinito di `decontam` è 0,1 | nessuno | `11_controls/decontam_riepilogo.json` |
| `decontam.min_blanks` | controlli negativi con letture minimi per decontaminare | dal numero di controlli negativi | nessuno | segnalazione `E-S12-03` |
| `decontam.mode` | `aggregate` (tutti contro tutti) o `batch` (per piastra) | `batch` solo con abbastanza negativi per piastra | `aggregate`: l'altra si calcola comunque | `decontam_riepilogo.json`: le due misure |

### 6.6 Soglie di qualità

Non hanno un valore di partenza: si fissano sui propri dati. Per una prima esecuzione
esplorativa si può dichiarare la soglia che non ferma nulla, indicata nell'ultima
colonna, leggere nel report la misura, e poi fissare la soglia e lanciare `resume`.
Una soglia lasciata al valore che non ferma è un controllo spento.

| Parametro | Significato | Misura nel report | Non ferma |
|---|---|---|---|
| `qc.max_zeroed_samples` | campioni che il filtro può lasciare senza letture (partenza di `init`: 0) | tracciamento delle letture | un numero alto |
| `qc.max_frac_lost_filter` | perdita media massima di letture nel filtro | tracciamento delle letture | 1 |
| `qc.max_asv_count` | varianti oltre le quali la tabella non si costruisce (partenza di `init`: dalla memoria della macchina) | varianti dopo l'inferenza | un numero alto |
| `qc.warn_frac_chimeric` | frazione di letture chimeriche che dà un avviso | frazione chimerica | 1 |
| `qc.stop_frac_chimeric` | frazione di letture chimeriche che ferma | frazione chimerica | 1 |
| `qc.min_frac_phylum` | frazione minima di varianti con il phylum assegnato | varianti con il phylum | 0 |
| `qc.max_frac_contaminant` | frazione massima di letture rimosse come contaminanti | `decontam_riepilogo.json` | sotto 1 |
| `qc.min_reads_final` | letture minime di un campione nell'oggetto finale: sotto, il campione esce | letture per campione | 0 |
| `qc.min_frac_reads_retained` | frazione minima di letture trattenute dall'insieme | letture finali su letture senza chimere | 0 |

`qc.min_reads_final` cambia il risultato, non solo un arresto: decide quali campioni
restano.
