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
- Un file può chiudersi con righe vuote. Una riga finale di soli spazi è un difetto del
  file: la validazione la respinge (`E-S0-13`) se cade fra le letture che ispeziona,
  altrimenti la respinge S1 (`E-S1-05`).

### 1.2 Letture appaiate

La pipeline tratta dati single-end. Di una corsa appaiata si usano le sole letture
forward: un file per campione con la prima lettura di ogni coppia, in una cartella a
parte. I file delle seconde letture (`_R2`, `_2`) non vanno nella cartella. Se un
archivio distribuisce le due letture di ogni coppia nello stesso file, vanno separate
prima: `dati/osd276/scarica_letture.py` mostra come, per un dataset pubblico.

Il file forward è quello le cui letture partono dal primer forward, non
necessariamente quello che il deposito chiama R1 o `_1`: in alcuni depositi i due file
di una corsa sono scambiati, e con il file sbagliato la pipeline leggerebbe l'altra
estremità dell'amplicone senza alcun errore. `dati/letture_forward.py` (usato da
`dati/osd170/scarica_letture.py` e `dati/comunita_nota/scarica_letture.py`) lo verifica
sulle prime letture dei due file, tiene quello che porta il primer e dichiara le corse
invertite. Se le letture sono depositate senza primer il controllo non è possibile: si
tiene il file dichiarato forward, e conviene guardare il motivo che `init` propone
(sezione 6.3), che dice con quale regione iniziano le letture. Il file tenuto va salvato
con un nome senza `R1`, `R2`, `_1` o `_2`: quei marcatori dicono alla validazione che il
file fa parte di una coppia.

`init` e la validazione riconoscono le letture appaiate dal nome dei file e dalle
intestazioni delle letture, e si fermano senza proporre né eseguire nulla.

Un file con le sole seconde letture, o con tutte le prime seguite da tutte le seconde,
può non avere alcun segno nelle intestazioni. Lo riconosce il contenuto: la fase S1
conta, su tutte le letture di ogni file, quelle che portano il segnale di una lettura
forward (il motivo conservato, `qc.conserved_motif`, dopo le basi di
`filter.trimLeft`; con il taglio, anche il primer in testa) e si ferma con `E-S1-07`,
nominando il file e la quota misurata, in due casi:

- un file di un campione biologico o di un controllo positivo con almeno 500 letture in
  cui meno di una lettura su mille porta il segnale: sono le letture inverse. Vale
  anche se ogni file del dataset è così, e qualunque siano le soglie di qualità;
- un file, di qualunque tipo, in cui una parte ha il segnale (almeno 30 letture) e
  un'altra meno di un decimo della sua quota: riporta le letture dei due versi, a
  blocchi (il file si confronta in dieci tratti consecutivi) o alternate (le letture
  di posto dispari con quelle di posto pari).

Due casi in cui il controllo ferma un file di sole letture forward: un campione che non
ha quasi nulla del bersaglio (va tolto dal dataset: non avrebbe dato varianti utili), e
un file riordinato per sequenza, in cui le letture col motivo stanno tutte insieme (va
riportato a un ordine non legato alla sequenza).

Il controllo ha bisogno del segnale: con `qc.conserved_motif` nullo e nessun primer
nelle basi tagliate non si applica, e l'esecuzione lo dichiara nel log. Il motivo deve
venire da chi conosce la regione amplificata, o dalla proposta di `init`. `init` ricava
il motivo dal consenso dei file che iniziano come la maggioranza: un file che inizia in
un altro modo (un file di letture inverse fra file forward, per esempio) esce dal
consenso ed è nominato nel riepilogo e nel commento del motivo, e va guardato prima di
eseguire. Se però la maggior parte dei file è di letture inverse, il motivo proposto
descrive quelle, e il controllo passa. Non
sono riconosciuti un file di letture inverse con meno di 500 letture, uno di un
controllo negativo, e una libreria in cui i due versi sono mescolati in ogni file.

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
Se il riferimento ha un elenco di taxa con un difetto noto (`tax.ref_bad_taxa`), il
file è una tabella CSV con le colonne `name` e `rank`: la validazione lo verifica a
ogni avvio (`E-S0-12`). L'MD5 da dichiarare in `tax.ref_md5` è quello di `md5sum FILE`,
fra virgolette; `init --reference` lo calcola da sé.

L'intestazione di ogni sequenza porta i ranghi separati da punto e virgola, dal regno
in giù. L'assegnazione dà loro un nome per posizione (`Kingdom`, `Phylum`, `Class`,
`Order`, `Family`, `Genus`, `Species`), qualunque nome abbiano nel riferimento: il
phylum è il secondo. Servono da due a sette ranghi; la validazione respinge
(`E-S0-22`) un riferimento in cui ogni sequenza ne ha uno solo, o in cui una ne ha più
di sette. Vanno bene i training set di SILVA, RDP, GTDB e Greengenes per dada2; come
scrivere i filtri tassonomici per ciascuno è nella sezione 6.4.


## 2. Scrivere la configurazione con `init`

`init` guarda letture e metadati e scrive un file di configurazione di partenza. Non
scarica nulla, non modifica i file di ingresso e non scrive su un file che esiste già.

### Il percorso consigliato, in cinque passi

Lo strumento propone, chi conduce l'analisi conferma. La conferma può essere un solo
gesto per tutto il file, l'opzione `--accetta-proposte`; restano da scrivere le scelte
che nessuno strumento può fare al posto di chi conosce lo studio.

1. **`init --accetta-proposte`** scrive come righe attive ciò che ha rilevato, ciò che
   propone e i valori di partenza; le soglie di qualità senza un valore di partenza
   prendono il valore che non ferma (`[ESPLORATIVA]`, sezione 6.6).
   Se `init` non trova il primer nelle letture né la sua sequenza nei metadati, lo
   scrive nullo: la validazione non cerca allora alcun primer, e il motivo conservato
   resta nullo. Se il primer è noto conviene rilanciare `init` con `--primer SEQUENZA`
   e un altro `--output` prima di proseguire: solo così propone anche il motivo
   (sezione 6.3), che serve a verificare la regione amplificata e a riconoscere un
   file di letture inverse (sezione 1.2). Scritto a mano nel file, il primer lascia il
   motivo nullo (`E-S0-18` a ogni esecuzione).
2. **Le scelte.** `validate` elenca i parametri che mancano (`E-G15-10`): sono le
   scelte senza valore di partenza (di solito i filtri tassonomici, la prevalenza, le
   letture minime, con i controlli negativi la soglia di decontaminazione). Si
   scrivono nel file, con la sezione 6 accanto; le righe attive si controllano e,
   dove serve, si correggono. Vanno guardati in particolare i valori di partenza che
   dichiarano un'assenza (nullo, o un elenco vuoto), che `init` elenca in testa al
   file: sono validi, ma se lo studio ha quell'informazione va scritta. Senza
   `meta.module_column` l'oggetto finale non ha la posizione di prelievo, e senza
   `out.study_columns` non ha le variabili dello studio che le analisi a valle
   confrontano.
3. **Una prima esecuzione esplorativa**, con `run`. Alla fine il comando elenca i
   controlli spenti: è atteso.
4. **`soglie --accetta-proposte`** legge le misure di quell'esecuzione e scrive un
   nuovo file con le soglie proposte, ciascuna accanto alla sua misura (sezione 6.6).
5. **`resume`** con il nuovo file: le soglie diventano attive.

Con i comandi (la forma di `docker run` è quella mostrata sotto; qui il solo
sottocomando):

```bash
python3 scripts/esegui.py init --reads /dati/letture --metadata /dati/campioni.tsv \
    --reference /dati/riferimento.fa.gz --out-root /lavoro/esecuzione \
    --output /lavoro/configurazione.yaml --accetta-proposte
python3 scripts/esegui.py validate --config /lavoro/configurazione.yaml   # elenca le scelte
# ... si scrivono le scelte nel file, poi, in una cartella di uscita vuota:
python3 scripts/esegui.py run --config /lavoro/configurazione.yaml
python3 scripts/esegui.py soglie --config /lavoro/configurazione.yaml \
    --output /lavoro/con_soglie.yaml --accetta-proposte
python3 scripts/esegui.py resume --config /lavoro/con_soglie.yaml
python3 scripts/esegui.py report --config /lavoro/con_soglie.yaml
```

`validate` scrive il suo esito nella cartella di uscita, e `run` parte solo da una
cartella vuota: dopo `validate` si cancella la cartella di uscita, oppure si prosegue
con `resume` al posto di `run`. La ripresa del passo 5 non cambia conteggi, tassonomia,
campioni né varianti, ma ripete i calcoli: ogni fase dichiara le soglie che applica, e
la prima a cambiare è la validazione, da cui tutte dipendono. Dura quanto la prima
esecuzione.

L'alternativa è accettare le proposte una per una: `init` senza l'opzione scrive
proposte e valori di partenza in commento, e nulla entra nella configurazione finché
non si toglie il commento alla sua riga (sezioni 2.1 e 2.5). Le due strade portano
allo stesso file; cambia quando si guarda ogni valore.

### Il comando

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
| `--accetta-proposte` | scrive come righe attive le proposte, i valori di partenza e le soglie esplorative; senza, restano in commento |

I percorsi scritti nella configurazione sono quelli visti dal container (`/dati/...`,
`/lavoro/...`, `/app/...` per ciò che sta nel clone): ogni comando successivo va
lanciato con gli stessi `-v`. Dati in più cartelle si montano con un `-v` ciascuna.

`init` legge tutte le letture per riconoscere le coppie e cercare il primer, e calcola
il profilo di qualità con lo stesso calcolo della fase S1: su 960 campioni e 21 milioni
di letture impiega alcuni minuti. Esce con 0 se ha scritto il file, con 3 se si rifiuta
(file già esistente, letture appaiate, ingressi non trovati): in quel caso non scrive
nulla.

### 2.1 Le specie di valore

Ogni riga del file dice da dove viene il suo valore, nel commento che la precede:

| Marcatore | Che cos'è | Come è scritta |
|---|---|---|
| `[RILEVATO]` | un fatto verificato su tutti i file | riga attiva |
| `[INDICATO]` | un valore dato sulla riga di comando: un percorso, il primer di `--primer`, o l'assenza di un file non indicato | riga attiva |
| `[PROPOSTO]` | l'esito di un'euristica, che può sbagliare | riga in commento, con il valore già scritto |
| `[DA SCEGLIERE]` | una scelta di chi conduce l'analisi | riga in commento, con un valore di partenza e il suo motivo, oppure senza valore |
| `[CONDIZIONALE]` | un parametro che serve solo se una condizione è vera, e con i valori del file non lo è | riga in commento, senza valore, con la condizione: non va scelto |

Con `--accetta-proposte` le righe `[PROPOSTO]` sono attive, e compaiono due marcatori:

| Marcatore | Che cos'è | Come è scritta |
|---|---|---|
| `[PARTENZA]` | il valore di partenza di una scelta (il predefinito del metodo, o un valore ricavato dai dati), accettato | riga attiva |
| `[ESPLORATIVA]` | una soglia di qualità al valore che non ferma: il suo controllo è spento | riga attiva; `soglie` la sostituisce dopo la prima esecuzione |

Restano `[DA SCEGLIERE]`, in commento, le sole scelte senza un valore di partenza. Ogni
riga attiva porta nel commento la sua origine e la fonte, e va comunque letta: una
proposta può sbagliare (sezione 2.3).

I parametri che hanno un predefinito (`[STANDARD]`, `[UNICA FORMA]`, `[STRUMENTO]`,
`[FACOLTATIVO]`) sono in commento con quello: si possono lasciare così. Un parametro
`[UNICA FORMA]` ammette un solo valore, perché la pipeline realizza quella sola forma
del metodo: `katharoseq.curve_model`, per esempio, nomina la curva di calibrazione,
che è una, e scriverlo non cambia nulla.

Nessun valore entra nella configurazione senza un gesto: una proposta si accetta
togliendo il commento alla sua riga (`  # column: "tipo"` diventa `  column: "tipo"`),
oppure tutte insieme con `--accetta-proposte`; una scelta si scrive. Le spiegazioni cominciano con `##`, le righe dei parametri in
commento con un solo `#`. In testa al file, e sulla console, c'è quanti valori sono
rilevati, proposti e da scegliere, l'elenco di questi ultimi, l'elenco dei condizionali
non richiesti con le proposte attuali, e ciò che `init` dichiara di non poter stabilire.

Un parametro **condizionale** è obbligatorio solo quando la sua condizione è vera: i
parametri della calibrazione servono solo se `ctrl.positive_values` non è vuoto, quelli
della decontaminazione solo se il dataset dichiara controlli negativi
(`ctrl.blank_values` o `ctrl.blank_override_values` non vuoto), la frazione e il
conteggio minimi della prevalenza solo con `prev.apply` vero, il limite delle varianti
solo con `phylo.enabled` vero (l'elenco è nella sezione 6). `init` valuta la condizione
sui valori che ha rilevato, proposto o scritto come partenza; se il parametro che la
decide non ha ancora un valore, la considera falsa: `prev.apply` è una scelta senza
valore di partenza, quindi `prev.min_fraction` e `prev.min_count` compaiono fra i non
richiesti finché non lo si dichiara vero. Un condizionale non richiesto si lascia in
commento. Se correggendo una proposta o facendo una scelta la condizione diventa vera
(per esempio si aggiungono le etichette dei controlli positivi), `validate` chiede il
condizionale, e allora si scrive come ogni altra scelta.

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
  controlli e senza una colonna del tipo la proposta è `ctrl.column: null`. La proposta
  guarda le parole, non il disegno dello studio: se i campioni da analizzare sono essi
  stessi i controlli di un altro studio (le repliche di uno standard, per esempio) li
  propone tutti come controlli positivi, e nessuno resterebbe nell'oggetto finale. In
  quel caso si dichiara `ctrl.column: null`, oppure si mettono quei valori fra le
  etichette dei campioni biologici.
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
  l'indicazione manca). Senza una sequenza nei metadati il primer è una scelta con
  partenza nulla, e il commento elenca le forme che `init` conosce del primer che i
  metadati nominano.
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

Un parametro condizionale non pertinente si omette: la sua riga resta in commento.
Dichiarato comunque non cambia i risultati, e il report lo segnala; il valore però è
controllato come ogni altro, anche se non serve: per esempio `decontam.threshold` a 0
o a 1 è respinto (`E-G15-02`), `decontam.mode: batch` senza `decontam.batch_column`
pure (`E-G15-11`), e così `meta.module_regex` senza una colonna da cui derivare il
modulo (`E-G15-14`). Per quattro elenchi, con la condizione falsa, è accettato solo
l'elenco vuoto: le tre etichette (`ctrl.blank_values`, `ctrl.positive_values`,
`ctrl.biological_values`) con `ctrl.column` nullo, e `ctrl.blank_override_values` senza
`ctrl.blank_override_column`; un valore non vuoto è respinto, perché farebbe credere
che qualche campione sia un controllo. Gli altri obbligatori si dichiarano sempre; dove
il dataset non ha l'informazione il valore è nullo (`null`) o vuoto (`[]`), come dice
la sezione 6. Gli elenchi si scrivono in linea (`["a", "b"]`), le stringhe fra
virgolette.

Omettere un condizionale non pertinente e dichiararlo con un valore non sono la stessa
configurazione: il digest è diverso (è lo stesso solo se si dichiara il valore che
vale l'omissione: l'elenco vuoto, o `null`), e passare dall'una all'altra forma su
un'esecuzione già fatta rifà le fasi che dichiarano quel parametro e quelle che ne
dipendono. Conteggi, tassonomia, campioni e varianti sono
identici nei due casi; che cosa cambia negli altri artefatti dipende dal parametro:

- per le etichette, `ctrl.blank_override_values`, `meta.module_regex`,
  `phylo.max_seqs` e i parametri dei controlli positivi (`ctrl.min_positives`,
  `ctrl.min_positive_pass_frac`, `ctrl.positive_gate`, `katharoseq.target_taxon`,
  `katharoseq.collapse_rank`, `katharoseq.min_r2`) cambia solo la configurazione
  registrata, `00_config/resolved.yaml`;
- per `katharoseq.target_sensitivity`, `decontam.threshold`, `decontam.min_blanks`,
  `decontam.mode`, `qc.max_frac_contaminant`, `prev.min_fraction` e `prev.min_count`
  cambiano anche i campi che li riportano nei riepiloghi delle fasi
  (`11_controls/soglia.json`, `11_controls/decontam_riepilogo.json`,
  `12_final/intermedi/filtri_riepilogo.json`): dichiarati ne portano il valore, omessi
  valgono `null`, e il dettaglio dei messaggi che riportavano il valore di
  `decontam.min_blanks` (gli avvisi sui controlli negativi in
  `01_input_validation/gates.json`, l'esito della decontaminazione) non lo riporta. Il
  report mostra quei campi come «non pertinente».

Una proposta è ciò che i dati suggeriscono, non ciò che uno studio ha scelto: sui due
dataset OSD-734 e OSD-276 alcune proposte differiscono dalla configurazione del repository.
Per OSD-276 `init` propone `filter.truncLen: 151`, la lunghezza di tutte le letture,
mentre la configurazione del repository tronca a 150 (è il valore dell'esempio, non
una scelta degli autori dello studio): i due valori danno oggetti finali
diversi (990 e 1.005 varianti). Chi vuole riprodurre un risultato usa i valori della
configurazione che lo ha prodotto.

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
| 6 | `soglie --config FILE --output NUOVO` | da un'esecuzione conclusa, scrive un nuovo file con le soglie di qualità proposte dalle sue misure (sezione 6.6) |

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
- La prima esecuzione serve anche a fissare le soglie: `soglie` le propone dalle
  misure, e il report riporta le misure accanto a ciascuna (sezione 6.6). Cambiare una
  soglia e lanciare `resume` rifà la fase che la applica e quelle che ne dipendono.
- `validate` guarda le prime letture di ogni file: un file di letture inverse, o a
  blocchi, lo riconosce la fase S1, che le legge tutte (sezione 1.2).
- Una soglia al valore che non ferma è una configurazione valida con un controllo
  spento: `validate`, `run` e `resume` lo dicono in coda, con l'elenco delle soglie, e
  il report lo riporta fra le segnalazioni in apertura.

### 3.1 Eseguire con Apptainer

Dove Docker non è disponibile (i cluster di calcolo) la stessa immagine si esegue con
Apptainer, senza privilegi. L'immagine si converte una volta, **per digest**, in un file
`.sif`:

```bash
apptainer build amplicon16s.sif docker://<immagine>
```

`<immagine>` è quella indicata nel README, con il suo `@sha256:...`: il file `.sif`
contiene esattamente quell'immagine, e `run.container` nella configurazione resta lo
stesso. Ogni comando di questa guida ha poi la stessa forma, con `apptainer exec` al
posto di `docker run`:

```bash
apptainer exec --cleanenv --containall \
  --bind "$(pwd)":/app --bind /percorso/dei/dati:/dati:ro --bind /percorso/di/lavoro:/lavoro \
  --pwd /app amplicon16s.sif \
  python3 scripts/esegui.py validate --config /lavoro/configurazione.yaml
```

| Docker | Apptainer | Perché |
|---|---|---|
| `-v A:B`, `-v A:B:ro` | `--bind A:B`, `--bind A:B:ro` | i montaggi; i percorsi nel container sono gli stessi |
| `-w /app` | `--pwd /app` | la cartella di lavoro |
| `-u "$(id -u):$(id -g)"` | nulla | Apptainer esegue già con l'utente di chi lancia |
| isolamento per difetto | `--cleanenv --containall` | Apptainer per difetto porta nel container le variabili d'ambiente e la cartella personale di chi lancia (`PYTHONPATH`, `R_LIBS_USER`, `.Rprofile`): con le due opzioni il container vede solo ciò che è montato, come con Docker |

Con `--containall` la cartella temporanea del container è piccola e in memoria: se un
comando segnala spazio esaurito in `/tmp`, le si monta una cartella del disco
(`--bind /percorso/tmp:/tmp`). Le righe `INFO:` che Apptainer stampa prima dell'uscita
dei comandi sono sue e non riguardano la pipeline. I risultati non dipendono da quale dei due programmi
esegue l'immagine: gli artefatti di fase hanno gli stessi checksum, e differiscono
solo i campi che descrivono l'esecuzione (istante, identificativo, e le impronte che
ne derivano).

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
  `design.continuous_variables`, le colonne che sembrano misure (valori tutti numerici,
  almeno dieci distinti, distinti in almeno quattro quinti dei campioni);
- **da scegliere**: distanze (`beta.distances`), metodi e distanze dell'ordinazione
  (`ord.methods`, `ord.distances`), seme (`run.seed`), numero di taxa nei grafici
  (`comp.top_n`), dimensione minima dei gruppi (`stat.min_group_size`), livello di
  significatività (`stat.significance_level`), pseudoconteggio della distanza di
  Aitchison (`beta.clr_pseudocount`, solo con quella distanza). Hanno un valore di
  partenza solo `ord.nmds_trymax` (20, `vegan::metaMDS`) e `ord.pcoa_correction`
  (`none`, `ape::pcoa`).

Con `--accetta-proposte` le proposte sono righe attive (la variabile biologica, le
variabili tecniche, il rango, le variabili continue). Fanno eccezione i parametri che servono solo con un
metodo (`beta.unifrac_root`, `ord.pcoa_correction`, `ord.nmds_trymax`), i soli ad
avere un valore proposto o di partenza oltre a quelli: restano in commento con il loro
valore, perché dichiarati senza quel metodo sono respinti, e quali distanze e quali
ordinazioni usare è una scelta. Restano da scrivere le scelte; `eco validate` nomina
quelle che mancano, e dopo la scelta dei metodi i parametri che questi richiedono.

I nomi delle colonne sono quelli della tabella dei campioni dell'oggetto, non quelli
dei metadati di partenza: la posizione di prelievo è `posizione`, e una colonna portata
con `out.study_columns` prende un nome in minuscolo con i trattini bassi (`Factor
Value[Spaceflight]` diventa `factor_value_spaceflight`; la corrispondenza è in
`10_phyloseq/colonne_metadati.tsv`). `eco init` elenca i nomi che trova.

Il significato di ogni parametro è in `config/eco.example.yaml`; esempi completi e
commentati sono `dati/osd734/eco_osd734.yaml`, `dati/osd276/eco_osd276.yaml` e
`dati/osd170/eco_osd170.yaml`. Per le scelte senza un valore di partenza:

- `design.variable` va dichiarata anche quando nessuna colonna ha gruppi ripetuti: si
  indica la colonna che descrive i campioni (per esempio la posizione). I test fra
  gruppi non si eseguono, e l'analisi lo dichiara con un avviso; tabelle e grafici
  descrittivi si producono comunque;
- `comp.top_n`: quanti taxa distinguere nei grafici di composizione prima di riunire
  gli altri; conta la leggibilità (una decina di colori si distinguono, trenta no);
- `stat.min_group_size`: quanti campioni deve avere un gruppo per entrare nei test, non
  meno di 2; con gruppi più piccoli un test fra ranghi non può dare un valore piccolo
  di p, per quanto i gruppi differiscano;
- `stat.significance_level`: il livello sotto il quale le dispersioni dei gruppi si
  dichiarano diverse; è la soglia di significatività che lo studio adotta per i suoi
  test (di consueto 0,05);
- `beta.clr_pseudocount`, solo con la distanza di Aitchison: il valore sommato a ogni
  conteggio prima del logaritmo, perché gli zeri non hanno logaritmo; 1 è la scelta più
  diffusa, e un valore più piccolo dà più peso alle varianti rare;
- `run.seed`: un intero qualunque, che fissa rarefazione e permutazioni.

Completate le scelte:

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

### 4.1 Variabili continue

Una misura (una distanza, un tempo, una concentrazione) non si dichiara fra le
variabili categoriche: divisa in gruppi diventa un gruppo per campione, nessun test si
esegue, e l'analisi avvisa (`E-ECO-29`). Si dichiara in `design.continuous_variables`:

```yaml
design:
  variable: "trattamento"            # resta obbligatoria: una colonna categorica
  continuous_variables: ["distanza_m"]
```

I valori devono essere numeri, con il punto decimale (altrimenti `E-ECO-30`). Per ogni
variabile continua l'analisi scrive:

- `test/alfa_continua.tsv`: la correlazione per ranghi di Spearman fra ogni indice di
  alfa diversità e la variabile, con il valore p e la correzione di
  Benjamini-Hochberg sui tre indici;
- `alfa/continua_NN_<variabile>.png`: l'alfa diversità contro la variabile;
- `test/permanova_continua.tsv`: la PERMANOVA con la variabile come termine numerico,
  dopo le variabili tecniche. Il test delle dispersioni non si esegue (confronta
  gruppi, e una variabile continua non ne ha): la tabella lo dichiara. Con
  `stat.permanova_strata`, una variabile che ha un solo valore dentro ogni strato non
  ha la PERMANOVA (`E-ECO-31`): le permutazioni entro gli strati non la spostano.

La colonna deve essere nella tabella dei campioni dell'oggetto: una variabile dello
studio vi arriva dichiarandola in `out.study_columns` (o `out.batch_columns`, se sta
nel file del lotto) nella configurazione della pipeline, prima dell'esecuzione o con un
`resume`; `eco init` elenca in testa al file tutte le colonne che l'oggetto ha.

La correlazione non si esegue (avviso `E-ECO-31`, che dice che cosa manca) con meno
campioni di quanti ne servono perché una correlazione perfetta scenda sotto
`stat.significance_level` (5 campioni al 5%), o con meno di tre valori distinti: con
due soli valori la variabile separa due gruppi, e va dichiarata categorica.

### 4.2 Misure ripetute

Se gli stessi soggetti (individui, siti, colture) sono campionati in più condizioni, la
colonna che li identifica si dichiara in `design.subject`:

```yaml
design:
  variable: "tempo"
  subject: "paziente"
```

Ogni variabile categorica è allora analizzata anche **entro soggetto**, sui soli
soggetti che hanno un campione in ogni livello della variabile:

- `test/alfa_appaiata.tsv`: con due livelli il test di Wilcoxon dei ranghi con segno;
  con più di due il test di Friedman e un Wilcoxon dei ranghi con segno per ogni
  coppia di livelli;
- `test/permanova_appaiata.tsv`: la PERMANOVA con il soggetto come primo termine e le
  permutazioni ristrette entro il soggetto, con il test delle dispersioni. Con una
  distanza non euclidea (Bray-Curtis, Jaccard) il soggetto e la variabile possono
  spiegare insieme più della variabilità totale: il test non è allora valutabile, e la
  tabella lo dice; la distanza di Aitchison non ha questo limite;
- `test/soggetti.tsv`: ogni soggetto, i livelli che ha e se è completo o perché è
  escluso.

I confronti fra gruppi indipendenti (`test/alfa_diversita.tsv`, `test/permanova.tsv`)
restano, e non tengono conto del soggetto: con differenze fra soggetti grandi possono
non vedere un effetto che il test entro soggetto trova, e per una variabile che
distingue i soggetti i loro valori p sono troppo piccoli, perché contano come
indipendenti i campioni dello stesso soggetto (avviso `E-ECO-34`; lo stesso vale per
le variabili continue).

L'analisi entro soggetto di una variabile **non si esegue** (avviso `E-ECO-32`, che dice
che cosa manca) se un soggetto ha più di un campione nello stesso livello (se sono
repliche, se ne tiene una con `design.subset`, che seleziona i campioni per i valori
di una colonna, anche `accession`, o meglio con una colonna dei metadati che le
distingua), o se i soggetti completi sono meno del
minimo perché un test entro soggetto scenda sotto `stat.significance_level` (6 soggetti
al 5%). I soggetti senza un campione in ogni livello sono esclusi e dichiarati
(`E-ECO-33`).

Quando un'analisi di una variabile continua o entro soggetto non si esegue, la sua
tabella ha la sola intestazione: il motivo è nell'avviso, sulla console e in
`riepilogo.json`.

I due parametri sono facoltativi: una configurazione che non li dichiara dà gli stessi
file, byte per byte, e lo stesso digest.

### 4.3 Uscite

La cartella di uscita deve essere assente o vuota e fuori dalla cartella dell'oggetto.
Contiene tabelle (alfa diversità, composizione, distanze, ordinazione, test), grafici
PNG, `riepilogo.json` con gli avvisi e `manifest.json` con dimensione e SHA-256 di ogni
altro file. A parità di oggetto, configurazione e seme i file sono identici byte per byte.
Gli avvisi (`E-ECO-`) non fermano l'analisi: dichiarano una condizione che chi legge i
risultati deve conoscere, per esempio un gruppo troppo piccolo per i test, una
variabile che sembra continua, un'analisi entro soggetto che non si è potuta eseguire.

## 5. Integrazione con altri programmi

Quattro comandi scrivono sullo standard output un solo documento JSON: `schema`,
`validate --json` e `status --json` per la pipeline e, con `eco` davanti, per le
analisi ecologiche; `soglie --json` per la sola pipeline. Ogni documento ha `versione_formato`
(un intero: vale 1), `comando` e `strumento`. La versione cambia solo se un campo esistente
cambia significato o sparisce: chi legge deve ignorare i campi che non conosce. Con `--json`
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
solo se ne ha. La specie è `obbligatorio`, `condizionale`, `percorso`,
`standard` (con `fonte`), `unica_forma`, `strumento` o `facoltativo`. Un `condizionale`
porta `richiesto_da`, che dice quando è obbligatorio. Per la pipeline: `parametro` e
`condizione` (`non_nullo`, `non_vuoto` o `vero`, riferita al valore di quel parametro;
con `vero` c'è anche `valori: [true]`). Quando le vie che lo rendono obbligatorio sono
due c'è anche `oppure`, un elenco di altre condizioni della stessa forma: ne basta una
vera, fra `richiesto_da` e quelle di `oppure` (lo hanno i parametri della
decontaminazione, che servono se i controlli negativi sono dichiarati con le etichette
oppure con la riclassificazione). Se nessuna è vera il parametro si omette; dichiarato
comunque è controllato come ogni valore e non cambia i risultati (sezione 2.5). Per le
analisi ecologiche: `parametro` e `valori` (i metodi che lo rendono obbligatorio); se
nessuno è richiesto, dichiararlo è un errore. Le descrizioni non entrano nel digest
della configurazione.

Un condizionale ha `obbligatorio: false`, perché lo schema gli dà un valore di assenza:
un programma che costruisce un modulo di configurazione non deve fermarsi a quel campo.
Deve leggere `specie` e, per la specie `condizionale`, `richiesto_da`: il parametro va
chiesto quando la sua condizione, o una di quelle di `oppure`, è vera con i valori già
raccolti. Chi guardasse
solo `obbligatorio` non chiederebbe mai un condizionale, e `validate` respingerebbe la
configurazione con `E-G15-10`.

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
| `controlli_spenti` | per la pipeline: le soglie al valore oltre il quale il loro controllo non può scattare (sezione 6.6); per ciascuna `parametro`, `valore`, `limite`, `codici`, `ferma` (falso se il controllo è un avviso), `messaggio`. Nullo se la configurazione è stata respinta prima di valutarli |

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
| `fasi` | per ogni fase: `passo`, `descrizione`, `cartella`, `stato` (`conclusa`, `in_corso`, `da_eseguire`, `disattivata`, `sconosciuta`) e, se conclusa, `artefatti`, `versione`, `degradazioni`; `da_rifare: true` se il manifesto è di un'esecuzione precedente e quella in corso deve ricalcolare la fase |
| `fasi_concluse` | i passi conclusi, nell'ordine di esecuzione |
| `fasi_da_rifare` | i passi con `da_rifare`: hanno un manifesto, ma l'esecuzione in corso (o interrotta) li ricalcola |
| `fase_in_corso` | il passo in esecuzione, o nullo |
| `esecuzione` | se un'esecuzione ha lasciato il suo stato: `passo`, `avviata`, `ultimo_battito`, `secondi_dall_ultimo_battito`, `battito_recente`, `da_rifare` (le fasi che deve ancora calcolare) |
| `arresto` | se c'è un punto di ripresa: `passo`, `codice`, `sintesi`, `dettaglio`, `azione`, `comando` |
| `configurazione` | file e digest dell'ultima configurazione registrata |

- `in_corso`: l'esecuzione aggiorna il suo stato ogni 20 secondi. Se l'ultimo
  aggiornamento è più vecchio di un minuto lo stato è `interrotta`: il processo è stato
  ucciso senza poter dichiarare nulla, e si riprende con `resume`.
- `ferma`: l'esecuzione si è fermata con un codice; `arresto` dice che cosa fare.
- `incompleta`: alcune fasi sono concluse e nessun processo è attivo, per esempio dopo
  `validate`, che esegue la sola prima fase.
- Una fase risulta conclusa se il suo manifesto esiste: se i suoi artefatti siano
  ancora integri e validi per una configurazione lo stabilisce `resume`. Durante una
  ripresa che ricalcola le fasi (dopo `soglie`, o dopo una modifica alla
  configurazione) i manifesti dell'esecuzione precedente restano sul disco finché la
  fase non è rifatta: quelle fasi hanno `stato` `conclusa` e in più `da_rifare`, e
  `status` senza `--json` le mostra come `da_rifare`.

Per le analisi ecologiche lo stato è `non_avviata`, `in_corso`, `interrotta` o
`conclusa` (con `file`, il numero dei file prodotti compreso il manifesto, `avvisi`,
`configurazione`, `oggetto_di_partenza`). Un rifiuto
non lascia nulla nella cartella: lo riporta il comando stesso, e `eco validate --json`
lo anticipa.

### 5.4 `soglie --json`

Nel container, come sopra:

```bash
python3 scripts/esegui.py soglie --config configurazione.yaml --output nuova.yaml --json
```

Il comando è lo stesso di `soglie` (sezione 6.6): scrive il file e, sullo standard
output, il documento.

| Campo | Contenuto |
|---|---|
| `configurazione`, `cartella` | il file letto e la cartella dell'esecuzione |
| `file_scritto` | il file scritto, o nullo se il comando si è rifiutato |
| `proposte_accettate` | vero con `--accetta-proposte`: le soglie proposte sono righe attive |
| `codice_uscita`, `rifiuto` | 0 e nullo, oppure 3 e il motivo del rifiuto |
| `soglie` | un elemento per soglia (vedi sotto) |
| `controlli_spenti` | i parametri il cui controllo resta spento anche con il valore del file scritto |
| `avvertenza` | che cosa una soglia proposta protegge e che cosa no |

Ogni elemento di `soglie`: `parametro`; `esito` (`proposta`, `tenuta`, `senza_misura`,
`esclusa`, `omessa`, `non_proposta`); `direzione` (`massima` o `minima`), `unita` (`frazione`,
`conteggio`, `basi`), `ferma` e `codici` del controllo; `misura`, `statistica` (che
cosa è stato misurato) e `origine_della_misura` (l'artefatto); `valore_configurazione`
(nel file letto), `valore_suggerito` (misura e margine), `valore` (quello che il file
scritto dichiara); `margine`, `regola` (a parole), `controllo_spento`, `note`.

## 6. I parametri obbligatori e condizionali, uno per uno

Per ogni parametro: che cosa significa; come si ricava dai propri dati; il valore di
partenza che `init` scrive, con la sua fonte; che cosa guardare nel report della prima
esecuzione per confermarlo. «Nessuno» nella colonna della partenza vuol dire che `init`
lascia la riga senza valore: non c'è un valore difendibile senza conoscere lo studio.
Le sezioni del report sono quelle descritte nel README.

I parametri **condizionali** sono indicati nella colonna del significato con la loro
condizione («solo se ...»): con la condizione falsa si omettono. Sono le tre etichette
dei campioni (solo se `ctrl.column` non è nullo), `ctrl.blank_override_values` (solo se
`ctrl.blank_override_column` non è nullo), `meta.module_regex` (solo se
`meta.derive_module` è vero), `phylo.max_seqs` (solo se `phylo.enabled` è vero),
`prev.min_fraction` e `prev.min_count` (solo se `prev.apply` è vero); solo se
`ctrl.positive_values` non è vuoto, `ctrl.min_positives`,
`ctrl.min_positive_pass_frac`, `ctrl.positive_gate`, `katharoseq.target_taxon`,
`katharoseq.collapse_rank`, `katharoseq.min_r2` e `katharoseq.target_sensitivity`; e,
solo se il dataset dichiara controlli negativi (`ctrl.blank_values` non vuoto oppure
`ctrl.blank_override_values` non vuoto), `decontam.threshold`, `decontam.min_blanks`,
`decontam.mode` e `qc.max_frac_contaminant`.

Le configurazioni dei dataset del repository sono esempi commentati di ogni scelta:
`dati/osd734/config_osd734.yaml` (960 campioni, controlli negativi e positivi, file del
lotto), `dati/osd276/config_osd276.yaml` (15 campioni, nessun controllo),
`dati/osd170/config_osd170.yaml` (30 campioni, due controlli negativi, condizionali
omessi) e `dati/comunita_nota/config_comunita_nota.yaml` (uno standard di composizione
nota, tabella dei campioni in TSV semplice).

### 6.1 File e campioni

| Parametro | Significato | Come si ricava | Partenza (fonte) | Conferma |
|---|---|---|---|---|
| `io.fastq_glob` | modello dei nomi dei file delle letture | dall'estensione dei file | rilevato | S0: tutti i file attesi nell'inventario |
| `io.accession_regex` | espressione che estrae la chiave del campione dal nome del file | la parte del nome unica per campione | rilevata, se una sola colonna corrisponde ai file | `01_input_validation/crosswalk.tsv`: una riga per campione |
| `meta.accession_column` | colonna da cui si ricava la stessa chiave | di solito quella con il nome del file | rilevata | come sopra |
| `meta.sample_id_column` | colonna con il nome del campione | la colonna che aggancia le due tabelle | rilevata con due tabelle, proposta con una | `crosswalk.tsv`, colonna `sample_name` |
| `ctrl.column` | colonna con il tipo del campione; nullo: nessun controllo | la colonna che distingue campioni e controlli | proposta (parole dei valori, o nome della colonna) | report, «Gate di validazione»: conteggi per classe |
| `ctrl.blank_values` | valori dei controlli negativi; `[]` se non ce ne sono; solo se `ctrl.column` non è nullo | copiati dalla colonna | proposti (blank, negative, NTC, no template) | conteggi per classe |
| `ctrl.positive_values` | valori dei controlli positivi; `[]` se non ce ne sono; solo se `ctrl.column` non è nullo | copiati dalla colonna | proposti (positive, mock) | conteggi per classe |
| `ctrl.biological_values` | valori dei campioni biologici; solo se `ctrl.column` non è nullo | copiati dalla colonna | proposti: i valori senza parole di controllo | conteggi per classe |
| `ctrl.blank_override_column` | colonna con cui trattare come controlli negativi campioni dichiarati di un altro tipo; nullo: nessuna | solo conoscendo lo studio (materiale mai usato, tamponi non aperti) | nullo (assenza) | `crosswalk.tsv`: classe e materiale |
| `ctrl.blank_override_values` | i valori di quella colonna che riclassificano; solo se `ctrl.blank_override_column` non è nullo | come sopra | non richiesto senza la colonna | come sopra |
| `meta.module_column` | colonna con la posizione di prelievo; nullo | la colonna della posizione, se lo studio ne ha una | nullo (assenza) | `crosswalk.tsv`, colonna `posizione` |
| `meta.derive_module` | se derivare il gruppo di posizioni dalla posizione | vero se le posizioni hanno un prefisso comune per gruppo | falso (assenza) | `crosswalk.tsv`, colonna `modulo` |
| `meta.module_regex` | espressione con un gruppo di cattura che estrae il gruppo dalla posizione; nullo; solo se `meta.derive_module` è vero | guardando i valori della posizione | non richiesto con `meta.derive_module` falso | come sopra |
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
| `qc.primer_sequence` | primer forward, in codici IUPAC; nullo se non è noto | dal protocollo di amplificazione (vedi sotto) | rilevato nelle letture; proposto dai metadati; indicato con `--primer`; altrimenti nullo | `01_input_validation/letture_ispezionate.tsv`, `frazione_primer` |
| `qc.conserved_motif` | espressione del motivo con cui inizia la regione amplificata; nullo | dal consenso delle prime basi (vedi sotto) | proposto se il primer è noto; altrimenti nullo | `letture_ispezionate.tsv`, `frazione_motivo` |
| `qc.max_primer_hit_frac` | frazione di letture con il primer in testa oltre la quale è considerato presente | bassa: il primer o c'è o non c'è | nessuno | `frazione_primer` per file |
| `qc.min_motif_frac` | frazione minima di letture con il motivo, sulla mediana dei campioni | sotto il valore osservato nei campioni buoni | nessuno | `frazione_motivo` per file |
| `qc.head_reads` | letture ispezionate per file dalla validazione | un compromesso fra durata della validazione e precisione delle frazioni misurate | 9.604: le letture con cui una frazione ha un intervallo di confidenza al 95% di un punto percentuale; le letture del file più grande, se sono meno | durata della validazione |
| `qc.max_frac_short_reads` | frazione massima di letture più corte del troncamento, per classe | zero se il troncamento non scarta letture | 0 con il troncamento proposto (dai dati) | report, «Troncamento suggerito» |
| `tax.ref_name`, `tax.ref_version` | nome e versione del riferimento, riportati con i risultati | dalla fonte del riferimento | proposti dal nome del file, se vi compaiono | report, «Configurazione» |
| `tax.try_rc` | se provare anche il complemento inverso | vero se le letture possono essere nell'altro orientamento | falso (`dada2::assignTaxonomy`) | report: frazione di varianti con il phylum |

**Il primer.** Il nome di un primer nei metadati non basta: lo stesso nome ha più
forme. Del 515F ce ne sono due, che differiscono per la quarta base: la forma
modificata (Parada 2016) porta `Y`, cioè C oppure T, dove la forma originale (Caporaso
2011) porta C. Si dichiara quella del protocollo di amplificazione usato; se il
protocollo dice solo «515F», si dichiara la forma modificata: comprende l'originale
(ogni lettura che inizia con il primer originale inizia anche con quello modificato),
quindi la verifica che il primer sia assente dalle letture vale per entrambe, e il
motivo conservato, che descrive le basi dopo il primer, è lo stesso. Le sequenze che `init` riconosce sono: 515F `GTGYCAGCMGCCGCGGTAA` (Parada 2016) e
`GTGCCAGCMGCCGCGGTAA` (Caporaso 2011); 341F `CCTACGGGNGGCWGCAG`; 357F
`CCTACGGGAGGCAGCAG`; 27F `AGAGTTTGATCMTGGCTCAG`; 8F `AGAGTTTGATCCTGGCTCAG`; 799F
`AACMGGATTAGATACCCKG`. Se i metadati nominano uno di questi `init` riporta le sue forme
nel commento. Per un primer diverso, o per farsi proporre il motivo quando il primer
non è nelle letture, si rilancia `init` con `--primer SEQUENZA` e un altro `--output`.

Con il primer nullo la validazione non cerca alcun primer in testa alle letture e lo
dichiara (`E-S0-21`); `qc.max_primer_hit_frac` non ha effetto. È la scelta di chi non
conosce il primer, non di chi lo conosce: se le letture lo portano ancora e nessuno lo
dichiara, entra nelle varianti. La sequenza si può scrivere in minuscolo.

I metadati possono dichiarare un primer che le letture non portano (tolto prima del
deposito), o persino il primer di un'altra regione: `init` propone allora la sequenza dei
metadati dicendo che non è nelle letture, e `filter.trimLeft: 0`. La validazione ne
verifica la sola assenza; quale regione le letture coprano davvero lo dice il motivo.

**Il motivo conservato.** È un'espressione regolare sulle prime basi della regione
amplificata, quelle che seguono il primer: una base dove quasi tutte le letture
concordano, una classe (`[AG]`) dove ne compaiono due, un punto dove variano. `init` lo
propone dal consenso delle letture quando conosce il primer. Senza una proposta la
scelta ragionevole è `null`: la validazione verifica allora la sola assenza del primer,
e S1 non può riconoscere un file di letture inverse (sezione 1.2). Le basi si scrivono
in maiuscolo, come le letture; scritte in minuscolo sono portate in maiuscolo, salvo in
un'espressione con gruppi con nome o con opzioni, che è respinta con la richiesta di
scriverle in maiuscolo.

### 6.4 Scelte di analisi

| Parametro | Significato | Come si ricava | Partenza (fonte) | Conferma |
|---|---|---|---|---|
| `err.randomize` | se scegliere a caso i campioni per il modello d'errore | vero se l'ordine dei file segue piastre o gruppi | falso (`dada2::learnErrors`) | `04_error_models/`: campioni usati |
| `dada.pool` | come i campioni condividono l'informazione: `false`, `true`, `pseudo` | `pseudo` se contano le varianti rare presenti in più campioni | falso (`dada2::dada`) | report: varianti e letture dopo l'inferenza |
| `chimera.min_fold_parent_over_abundance` | quanto un genitore deve superare una variante perché sia chimerica | il valore del metodo, salvo ragioni | 1,5 (`dada2::isBimeraDenovoTable`) | report: frazione di letture chimeriche |
| `filt.remove_na_phylum` | se togliere le varianti senza phylum | vero se ciò che il riferimento non riconosce è fuori bersaglio | nessuno | report, «Risultato finale»: varianti rimosse |
| `filt.exclude_taxa` | taxa da togliere, a qualunque rango, come li scrive il riferimento; `[]` | ciò che il gene cattura senza essere l'oggetto dello studio: di solito cloroplasti, mitocondri ed eucarioti (vedi sotto) | nessuno | `12_final/intermedi/varianti_rimosse.tsv` |
| `phylo.enabled` | se costruire l'albero delle varianti finali | vero se servono distanze filogenetiche | falso (fase facoltativa) | `09_phylogeny/` |
| `phylo.max_seqs` | varianti oltre le quali la filogenesi si ferma; solo se `phylo.enabled` è vero | dal tempo e dalla memoria disponibili | non richiesto con `phylo.enabled` falso; altrimenti nessuno | report: varianti finali |
| `prev.apply` | se filtrare le varianti presenti in pochi campioni | vero se servono le varianti condivise fra campioni; falso se conta ogni variante, anche di un solo campione (vedi sotto) | nessuno | report: varianti rimosse per prevalenza |
| `prev.min_fraction` | frazione minima di campioni in cui una variante deve comparire; solo se `prev.apply` è vero | più bassa con campioni eterogenei o numerosi (vedi sotto) | non richiesto finché `prev.apply` non è vero; poi nessuno | come sopra |
| `prev.min_count` | letture minime perché una variante conti come presente; solo se `prev.apply` è vero | 1 conta ogni lettura; 2 non conta le varianti con una sola lettura in un campione | non richiesto finché `prev.apply` non è vero; poi nessuno | come sopra |
| `qc.min_reads_final` | letture minime di un campione nell'oggetto finale: sotto, il campione esce | la profondità sotto la quale un campione non è confrontabile con gli altri; 0 non toglie nessuno | nessuno | report, «Risultato finale»: letture per campione, campioni esclusi |

**La prevalenza.** Il filtro tiene una variante se compare, con almeno `prev.min_count`
letture, in almeno `prev.min_fraction` dei campioni rimasti dopo i filtri precedenti,
arrotondato per eccesso: con 20 campioni e 0,2 servono 4 campioni, con 200 campioni e
0,02 ne servono ancora 4. Toglie le varianti viste in un solo campione, che sono spesso errori
o contaminazioni sporadiche, ma anche i taxa davvero rari: con pochi campioni, o con
campioni di ambienti diversi fra loro, una frazione alta toglie gran parte delle
varianti. Se la domanda riguarda la composizione di ogni singolo campione (per esempio
la verifica di una comunità nota) il filtro si lascia spento.

**Le letture minime.** `qc.min_reads_final` cambia il risultato, quindi si dichiara
prima dell'esecuzione. Se non c'è un valore da cui partire si può dichiarare 0 (nessun
campione escluso), guardare nel report della prima esecuzione le letture per campione,
e poi alzarlo: è una nuova scelta di analisi, da dichiarare come tale, e la ripresa
rifà i filtri finali.

**Che cosa togliere con i filtri tassonomici.** I filtri tolgono soltanto ciò che si
dichiara. I primer del 16S amplificano anche il 16S dei cloroplasti e dei mitocondri e,
per molti eucarioti, il gene 18S: in un campione con cellule vegetali, animali o
fungine quelle letture arrivano fino all'oggetto finale se nessun filtro le nomina. Di
solito si tolgono tutti e tre, scritti come li scrive il riferimento; per SILVA 138:

```yaml
filt:
  remove_na_phylum: true
  exclude_taxa: ["Chloroplast", "Mitochondria", "Eukaryota"]
```

Il nome si confronta come stringa esatta con ogni rango della tassonomia (senza un
eventuale prefisso di rango come `o__`), quindi va verificato nel file del riferimento,
dove le intestazioni portano i ranghi separati da punto e virgola:

```bash
zcat RIFERIMENTO.fa.gz | grep '^>' | sed 's/^>//' | tr ';' '\n' | sort -u \
  | grep -i -E 'chloroplast|mitochondria|eukaryot'
```

**Con un altro riferimento** i nomi cambiano, e con loro ciò che i filtri possono
togliere. Si guarda come il riferimento scrive i taxa, con il comando qui sopra, e si
scrivono così:

- un training set con i **prefissi di rango** (Greengenes: `p__Proteobacteria`,
  `c__Chloroplast`) si confronta senza il prefisso: `Chloroplast` e `c__Chloroplast`
  valgono lo stesso. Un rango scritto col solo prefisso (`g__`) è un rango che la
  sequenza non ha: non conta come assegnato, né per `filt.remove_na_phylum` né per
  `qc.min_frac_phylum`. In Greengenes i cloroplasti sono la classe `Chloroplast` e i
  mitocondri la famiglia `mitochondria`, in minuscolo, dentro i Rickettsiales: si
  scrivono così (`["Chloroplast", "mitochondria"]`);
- un training set **GTDB** ha sei ranghi senza prefissi, e i nomi portano i suffissi
  di GTDB (`Firmicutes_A`, `Blautia_A`): `Firmicutes` non toglie `Firmicutes_A`, vanno
  scritti entrambi. GTDB contiene solo genomi di batteri e archei: **non ha cloroplasti,
  mitocondri né eucarioti**, quindi nessun nome li toglie, e le loro letture vengono
  assegnate al taxon batterico più simile, anche con un bootstrap alto, o restano senza
  phylum. Se il campione può
  contenerne (tessuti, piante, superfici toccate) e vanno tolte, conviene un
  riferimento che le nomini; con GTDB `filt.remove_na_phylum: true` toglie solo quelle
  rimaste senza phylum.

Un nome di `filt.exclude_taxa` che il riferimento non ha non è un errore e non toglie
nulla: che cosa è stato tolto si controlla in `varianti_rimosse.tsv`.

Un riferimento può dare a una variante il solo regno: in SILVA 138 nel formato di dada2
le sequenze eucariotiche hanno `Eukaryota` e nessun altro rango. Una variante così non
ha phylum, quindi `filt.remove_na_phylum: true` la toglie anche senza `Eukaryota`
nell'elenco; toglie però con lei ogni variante batterica o archeale a cui il phylum non
è stato assegnato. Chi vuole tenere queste ultime dichiara `filt.remove_na_phylum:
false` e scrive `Eukaryota` in `filt.exclude_taxa`. Se la tassonomia non ha il rango
`Phylum`, `filt.remove_na_phylum` non si applica e l'esecuzione lo dichiara
(`E-S13-05`). Che cosa è stato tolto, e perché, è in
`12_final/intermedi/varianti_rimosse.tsv`.

### 6.5 Controlli positivi e negativi

Senza controlli positivi (`ctrl.positive_values` vuoto, o `ctrl.column` nullo) i
parametri della tabella indicati «solo con controlli positivi» non servono e si
omettono. Senza controlli negativi (`ctrl.blank_values` e `ctrl.blank_override_values`
vuoti) non servono e si omettono i tre parametri di `decontam` indicati «solo con
controlli negativi» e `qc.max_frac_contaminant` (sezione 6.6): un campione è un
controllo negativo per la sua etichetta oppure per la regola di riclassificazione, e
basta una delle due perché quei parametri siano richiesti. Si dichiarano sempre
`katharoseq.cell_count_column` (nullo senza colonna), `decontam.batch_column` (nullo
senza file del lotto) e `qc.min_reads_mode`.

Senza controlli negativi la decontaminazione non toglie nulla e lo dichiara
(`E-S12-03`); senza controlli positivi nessuna soglia di profondità si applica, e con
`qc.min_reads_mode: katharoseq_if_available` l'assenza è riportata come avviso
(`E-S11-05`), con `none` è una scelta dichiarata: l'oggetto finale è lo stesso. Un
condizionale dichiarato senza che serva è accettato e non cambia conteggi, tassonomia,
campioni né varianti (la configurazione di OSD-276, che non ha controlli, li dichiara
tutti); che cosa cambia fra ometterlo e dichiararlo è nella sezione 2.5.

| Parametro | Significato | Come si ricava | Partenza (fonte) | Conferma |
|---|---|---|---|---|
| `ctrl.min_positives` | controlli positivi minimi per una curva e per giudicare un livello; solo con controlli positivi | dal numero di controlli per piastra e per livello | 3: un punto più dei due parametri della curva | report, «Decisioni»: curve per piastra |
| `ctrl.min_positive_pass_frac` | frazione minima di controlli positivi conformi; solo con controlli positivi | secondo quanto ci si fida dei controlli; a 0 il controllo è spento (sezione 6.6) | nessuno | report: controlli non conformi |
| `ctrl.positive_gate` | vero: sotto la frazione l'esecuzione si ferma; falso: avvisa; solo con controlli positivi | falso finché la conformità non è nota | falso | segnalazioni `E-S11-03`, `E-S11-04` |
| `katharoseq.target_taxon` | taxon atteso nei controlli positivi, come lo scrive il riferimento (con o senza il prefisso di rango); solo con controlli positivi | da chi ha allestito i controlli | nessuno | `11_controls/`: fedeltà dei controlli |
| `katharoseq.collapse_rank` | rango a cui si sommano le varianti del taxon atteso; solo con controlli positivi | il più fine a cui il riferimento lo distingue | `Genus`: il più fine che la pipeline assegna | come sopra |
| `katharoseq.target_sensitivity` | fedeltà a cui si legge la profondità minima dei campioni; solo con controlli positivi | quanta contaminazione si accetta nei campioni poveri | nessuno | report: soglia per piastra, campioni esclusi |
| `katharoseq.min_r2` | bontà minima di una curva valida; solo con controlli positivi | secondo la qualità dei controlli | nessuno | report: R² delle curve |
| `qc.min_reads_mode` | da dove viene la profondità minima: `katharoseq_if_available` o `none` | `none` senza controlli positivi in diluizione | secondo i controlli della proposta | report: campioni esclusi per profondità |
| `decontam.threshold` | soglia sotto la quale una variante è un contaminante; solo con controlli negativi | più alta con bassa biomassa; il predefinito di `decontam` è 0,1 | nessuno | `11_controls/decontam_riepilogo.json` |
| `decontam.min_blanks` | controlli negativi con letture minimi per decontaminare; solo con controlli negativi | dal numero di controlli negativi | nessuno | segnalazione `E-S12-03` |
| `decontam.mode` | `aggregate` (tutti contro tutti) o `batch` (per piastra); solo con controlli negativi | `batch` solo con abbastanza negativi per piastra | `aggregate`: l'altra si calcola comunque | `decontam_riepilogo.json`: le due misure |

### 6.6 Soglie di qualità

I parametri di questa tabella hanno un solo effetto: fermare l'esecuzione o dare un
avviso quando una misura supera la soglia (o le sta sotto). Non cambiano conteggi,
tassonomia, campioni né varianti. Si fissano sui propri dati, dopo una prima
esecuzione; fanno eccezione gli ultimi due, che `soglie` non propone: uno protegge la
memoria, l'altro è un giudizio sui controlli positivi.

**Il valore che non ferma.** Ogni confronto è stretto: una misura uguale alla soglia
non fa scattare il controllo. Una soglia al limite della tabella, o oltre, è quindi un
controllo spento: la configurazione è valida, ma `validate`, `run`, `resume` e il
report lo segnalano. `init --accetta-proposte` scrive il valore che non ferma
(`[ESPLORATIVA]`) per le soglie in frazione che non hanno un valore di partenza e che
`soglie` propone poi dalla misura.

| Parametro | Significato | Scatta se | Non ferma (limite) | Partenza di `init` |
|---|---|---|---|---|
| `qc.max_primer_hit_frac` | frazione di letture con il primer in testa oltre la quale è considerato presente | la frazione di un file di un campione biologico o di un controllo positivo la supera (`E-S0-10`); con `filter.trimLeft` maggiore di zero, o con `qc.primer_sequence` nullo, il primer non si cerca e la soglia non ha effetto | 1 | nessuna |
| `qc.min_motif_frac` | frazione minima di letture con il motivo conservato | la mediana fra i file le sta sotto (`E-S0-16`); con `qc.conserved_motif` nullo non ha effetto | 0 | nessuna |
| `qc.max_frac_short_reads` | frazione massima di letture più corte del troncamento, per classe | la frazione di una classe la supera, sulle prime letture o su tutte (`E-S0-09`, `E-S1-02`) | 1 | 0 con il troncamento proposto |
| `filter.truncLen_shortfall_warn` | di quante basi il troncamento può stare sotto la lettura più corta | lo scarto la supera (avviso `E-S1-01`) | `filter.truncLen` meno 1, il massimo ammesso | lo scarto del troncamento proposto |
| `qc.max_zeroed_samples` | campioni che il filtro può lasciare senza letture | i campioni biologici e i controlli positivi azzerati sono di più (`E-S2-01`) | il numero dei campioni biologici e dei controlli positivi | 0 |
| `qc.max_frac_lost_filter` | perdita media massima di letture nel filtro | la perdita media la supera (`E-S2-02`) | 1 | nessuna |
| `qc.warn_frac_chimeric` | frazione di letture chimeriche che dà un avviso | la frazione chimerica la supera (avviso `E-S6-02`) | 1 | nessuna |
| `qc.stop_frac_chimeric` | frazione di letture chimeriche che ferma | la frazione chimerica la supera (`E-S6-01`) | 1 | nessuna |
| `qc.min_frac_phylum` | frazione minima di varianti con il phylum assegnato | la frazione di una classe le sta sotto (`E-S8-02`) | 0 | nessuna |
| `qc.max_frac_contaminant` | frazione massima di letture rimosse come contaminanti; solo con controlli negativi | la frazione rimossa dai biologici la supera (`E-S12-02`) | 1 | nessuna |
| `qc.min_frac_reads_retained` | frazione minima di letture trattenute dall'insieme | letture finali su letture senza chimere le sta sotto (`E-S14-01`) | 0 | nessuna |
| `qc.max_asv_count` | varianti oltre le quali la tabella non si costruisce: protegge la memoria | le varianti distinte sono di più (`E-S5-01`) | nessun limite: non è mai riportata come spenta | dalla memoria della macchina; `soglie` non la propone |
| `ctrl.min_positive_pass_frac` | frazione minima di controlli positivi conformi; solo con controlli positivi | i conformi sono meno: arresto `E-S11-03` con `ctrl.positive_gate` vero, altrimenti avviso `E-S11-04` | 0 | nessuna, mai: è una scelta; `soglie` non la propone |

**`soglie`.** Da un'esecuzione conclusa e dalla configurazione con cui è stata fatta
scrive un nuovo file, uguale al primo salvo le soglie; non ricalcola nulla e non
modifica la cartella.

```bash
python3 scripts/esegui.py soglie --config /lavoro/configurazione.yaml \
    --output /lavoro/con_soglie.yaml [--out CARTELLA] [--accetta-proposte] [--json]
```

Senza `--out` la cartella è `io.out_root` della configurazione. Senza
`--accetta-proposte` ogni soglia proposta è una riga in commento, e il parametro manca
finché non si toglie il commento; con l'opzione è attiva. Il comando si rifiuta
(uscita 3, senza scrivere) se il file di uscita esiste; se l'esecuzione non è conclusa;
se la configurazione non è quella con cui ogni fase è stata calcolata (dopo aver
cambiato un parametro va lanciato `resume`, non il solo `validate`); se un artefatto
manca o non si legge; se il file non ha un gruppo per blocco e un parametro per riga,
la sola forma che il comando sa riscrivere.

Per ogni soglia legge **la stessa misura che il controllo confronta con la soglia**
(il valore peggiore fra i file o fra le classi dove il controllo giudica ciascuno, la
mediana o la misura d'insieme dove giudica quella) e propone:

- per una soglia massima la misura più un margine, per una minima la misura meno un
  margine: 0,05 per le frazioni; 0,10 per l'arresto delle chimere, così che fra avviso
  e arresto resti una fascia; nessun margine per i campioni azzerati (ogni campione
  perso deve essere visto) e per le basi dello scarto;
- **non propone** `qc.max_asv_count` né `ctrl.min_positive_pass_frac`: restano al
  valore del file letto, e il comando lo dichiara con il motivo. Le varianti crescono
  in modo naturale con i campioni: una soglia ricavata dalla misura fermerebbe
  un'analisi con più campioni senza un peggioramento, e la soglia protegge la memoria.
  Quanti controlli positivi devono essere conformi è un giudizio, non una misura da
  proteggere;
- mai un valore più stretto della misura: la ripresa sugli stessi dati non avvisa né
  si ferma per una soglia appena scritta. Le frazioni hanno quattro decimali,
  arrotondate verso l'esterno. Una sola eccezione, dichiarata dal comando: se la
  lettura più corta supera il troncamento di più di `filter.truncLen` meno 1 basi, lo
  schema non ammette uno scarto più largo, e l'avviso `E-S1-01` continua a scattare;
- se la soglia del file letto è attiva, è rispettata dalla misura e non è più larga
  del suggerimento, **tiene quella** e lo dice;
- l'avviso delle chimere non supera mai l'arresto: se l'arresto tenuto è sotto
  l'avviso proposto, l'avviso è riportato all'arresto;
- se il valore proposto coincide con il limite (la misura è già lì, o misura e margine
  escono dall'intervallo), lo dichiara come **controllo spento**;
- una soglia **senza misura** in quell'esecuzione (motivo nullo, primer con il taglio,
  decontaminazione senza abbastanza controlli negativi) o la cui misura non si legge
  dagli artefatti (letture corte fra le prime letture di un file non letto per intero)
  resta com'è, e il comando lo dice; una condizionale omessa resta omessa.

Nel file scritto ogni soglia è preceduta da un commento: `[PROPOSTO] da soglie` con
misura, artefatto, regola e valore che aveva, per una soglia proposta; `[SOGLIE]` con
l'esito, per una soglia tenuta, senza misura, esclusa o non proposta, che conserva la
sua riga. In
testa al file c'è il riepilogo; le intestazioni che seguono, se il file veniva da
`init`, descrivono il file com'era prima.

**Che cosa protegge una soglia proposta, e che cosa no.** Non giudica se la misura è
buona: la protegge dai peggioramenti, per esempio dopo aver aggiunto campioni, cambiato
una scelta o aggiornato lo strumento. Il giudizio resta a chi conduce l'analisi: il
file, la console e il documento JSON riportano ogni misura accanto al valore proposto.
Una frazione chimerica alta, proposta con il suo margine, resta una frazione chimerica
alta: va guardata prima di accettarla.

`qc.min_reads_final` non è fra queste soglie: decide quali campioni restano nell'oggetto
finale (un campione con meno letture esce), quindi cambia il risultato e si sceglie
prima dell'esecuzione.
