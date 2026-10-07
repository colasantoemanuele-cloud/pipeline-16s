# Registro delle decisioni di metodo

Questo documento registra le decisioni di metodo della pipeline 16S rRNA single-end e
la regola con cui si scelgono i valori dei parametri sottoposti ad analisi di
sensibilità. È scritto per chi deve valutare o riprodurre i risultati senza aver
partecipato allo sviluppo.

## 1. Regola di decisione per l'analisi di sensibilità

La regola è fissata prima di eseguire qualunque calcolo di sensibilità: la storia del
repository lo attesta, perché questa sezione è registrata in una revisione che precede
ogni risultato. Serve a impedire che la scelta di un valore sia guidata dall'esito che
produce.

### 1.1 Che cosa si decide

Per ciascuno dei parametri elencati sotto si decide se mantenere il valore corrente
della configurazione di OSD-734 o sostituirlo con un altro valore della griglia. Si
varia un parametro alla volta, tenendo gli altri al valore corrente.

| Parametro | Valore corrente | Griglia |
|---|---|---|
| `katharoseq.target_sensitivity` | 0,90 | 0,70; 0,75; 0,80; 0,85; 0,90; 0,95 |
| `decontam.threshold` | 0,5 | 0,1; 0,2; 0,3; 0,4; 0,5; 0,6; 0,7 |
| `prev.min_fraction` | 0,01 | 0,005; 0,010; 0,015; 0,020; 0,025; 0,030 |
| modalità di decontaminazione (`decontam.mode`, `decontam.batch_combine`) | aggregata | aggregata; per piastra con `minimum`; per piastra con `fisher` |

Le griglie numeriche hanno passo uniforme, perché la regola confronta il cambiamento
fra valori adiacenti: con passi diversi i cambiamenti non sarebbero confrontabili. Il
valore corrente non sta mai all'estremo della griglia, così ha sempre due vicini. Per
`prev.min_fraction` si riporta in più il valore 0,05, consueto in letteratura: è fuori
dalla griglia di decisione perché dista dal valore corrente più passi di quanti la
regola ne confronti.

### 1.2 Il criterio

Il criterio è tecnico. Per ogni configurazione si misurano, sull'oggetto finale:

- l'insieme delle varianti finali, e la sua similarità di Jaccard con quello di
  un'altra configurazione;
- l'insieme dei campioni biologici conservati, e quanti campioni cambiano stato
  (conservati in una configurazione ed esclusi nell'altra, o viceversa);
- le letture conservate nell'oggetto finale.

Da queste misure si ricavano, per ogni coppia di valori adiacenti della griglia, tre
**cambiamenti di passo**:

- varianti: 1 meno la similarità di Jaccard fra i due insiemi di varianti finali;
- campioni: campioni biologici che cambiano stato, divisi per i campioni biologici del
  dataset;
- letture: differenza assoluta fra le letture conservate, divisa per le letture
  conservate nella configurazione corrente.

**Due esclusioni.**

- Non si usa la separazione fra gruppi sperimentali, né alcuna misura che dipenda dal
  disegno dello studio (moduli, superfici, voli): un parametro scelto perché separa
  meglio i gruppi renderebbe circolare ogni analisi successiva su quei gruppi.
- Il numero di campioni, di varianti o di letture conservati si riporta ma non si
  massimizza, e non entra nella decisione con il suo segno: conta quanto un risultato
  cambia fra valori vicini, non in quale direzione. Un criterio che premiasse chi
  conserva di più sceglierebbe sempre il valore più permissivo.

### 1.3 La regola per i parametri numerici

Principio: si mantiene il valore corrente, a meno che i risultati siano instabili
attorno a esso mentre un'alternativa si trova in una zona stabile.

1. **Ammissibilità.** Un valore è ammissibile se con esso la catena si conclude senza
   che scatti un arresto già previsto dalla pipeline (frazione di letture rimosse come
   contaminanti oltre `qc.max_frac_contaminant`, campioni svuotati dal filtro di
   prevalenza, frazione di letture trattenute sotto `qc.min_frac_reads_retained`). Un
   valore non ammissibile non può essere scelto, e il cambiamento di passo verso di
   esso non è definito.
2. **Instabilità del valore corrente.** Il valore corrente è instabile se, per almeno
   una delle tre misure, il cambiamento di passo verso almeno uno dei suoi due vicini
   supera **entrambe** le soglie seguenti:
   - una soglia assoluta: 0,10 per le varianti, 0,05 per i campioni, 0,05 per le
     letture;
   - una soglia di sproporzione: tre volte la mediana dei cambiamenti di passo di quella
     misura su tutte le coppie adiacenti ammissibili della griglia.
   Se uno dei due vicini non è ammissibile, il valore corrente è instabile da quel lato.
3. **Zona stabile.** Un valore alternativo è in una zona stabile se ha due vicini
   ammissibili nella griglia e, per tutte e tre le misure, i cambiamenti di passo verso
   entrambi restano sotto le soglie assolute.
4. **Decisione.** Se il valore corrente non è instabile, si mantiene. Se è instabile e
   nessun'alternativa è in una zona stabile, si mantiene e l'instabilità si dichiara
   come limite noto. Se è instabile ed esiste almeno un'alternativa in una zona
   stabile, si adotta quella più vicina al valore corrente; a parità di distanza,
   quella il cui insieme di varianti è più simile a quello corrente.

**Motivazione delle soglie.**

- Varianti, 0,10: passando al valore adiacente, più di una variante su dieci entra o
  esce dall'insieme finale. Sotto, due valori vicini danno insiemi che per le analisi a
  valle sono lo stesso insieme con un margine di varianti rare.
- Campioni, 0,05: più di un campione biologico su venti cambia stato. È l'ordine di
  grandezza oltre il quale la composizione del dataset analizzato, e non solo la sua
  profondità, dipende dal valore scelto.
- Letture, 0,05: la stessa soglia, sulle letture conservate.
- Sproporzione, tre volte la mediana: un cambiamento è sproporzionato se è nettamente
  maggiore di quello tipico fra valori vicini dello stesso parametro. Un parametro che
  cambia i risultati a ogni passo in misura simile non ha una zona più stabile di
  un'altra, e spostare il valore non lo renderebbe più affidabile: per questo la soglia
  assoluta da sola non basta. Il fattore tre separa uno scarto netto dalla variabilità
  ordinaria fra i passi di una griglia di sei o sette valori, in cui la mediana è
  calcolata su cinque o sei coppie.

Le soglie assolute sono volutamente larghe rispetto alla precisione delle misure, che
sono conteggi esatti e non stime: non c'è rumore di misura da cui proteggersi, c'è solo
da stabilire quanto cambiamento sia rilevante.

### 1.4 La regola per la modalità di decontaminazione

Le tre modalità non sono valori ordinati di una stessa grandezza: non hanno vicini, e il
criterio di stabilità fra valori adiacenti non si applica. Vale la regola seguente.

- Una modalità è ammissibile se la catena si conclude con essa senza arresti, nel senso
  del punto 1 sopra, con le soglie di controllo della configurazione corrente.
- Si mantiene la modalità corrente se è ammissibile. La si sostituisce solo se non lo
  è, con la modalità ammissibile il cui insieme di varianti finali è più simile a
  quello corrente.
- Le misure delle modalità non scelte si riportano come diagnostica.

Motivazione: fra modalità che rispondono a domande statistiche diverse (un confronto su
tutti i controlli negativi, o uno per piastra con pochi negativi ciascuno) non esiste un
criterio tecnico di stabilità che ne indichi una; la scelta è di metodo, e la si cambia
solo se quella adottata viola un vincolo che la pipeline già dichiara.

### 1.5 Il modello evolutivo della filogenesi

La fase di filogenesi è disattivata nella configurazione di OSD-734; si decide il valore
predefinito di `phylo.model` fra GTR+G+I e GTR+G.

- I due modelli si adattano alle varianti finali dell'esecuzione di riferimento, con lo
  stesso allineamento, lo stesso albero di partenza e la stessa ricerca.
- Si confrontano con il criterio d'informazione bayesiano (BIC), calcolato sul numero
  di colonne dell'allineamento; si riporta anche il criterio di Akaike (AIC).
- Si adotta il modello con il BIC più basso. Se la differenza di BIC è inferiore a 2,
  i due modelli non sono distinguibili e si adotta il più semplice, GTR+G.

Motivazione: il modello con i siti invarianti ha un parametro in più, e la quota di
siti invarianti e la forma della distribuzione gamma descrivono in parte lo stesso
fenomeno, la presenza di siti che variano poco; un criterio d'informazione dice se il
parametro in più è giustificato dai dati. Il BIC è preferito all'AIC come criterio di
decisione perché penalizza di più i parametri aggiuntivi, ed è la scelta prudente
quando l'allineamento ha poche colonne.

### 1.6 Che cosa era già noto quando la regola è stata fissata

Per trasparenza: alcune misure sull'esecuzione di riferimento esistevano già, perché
la pipeline le produce come diagnostica a ogni esecuzione.

- La modalità di decontaminazione per piastra con `minimum` rimuove circa il 44% delle
  letture dei campioni biologici, contro il 4,9% della modalità aggregata; è la ragione
  per cui la modalità aggregata era stata adottata, ed è superiore a
  `qc.max_frac_contaminant` (0,40).
- Nel modello GTR+G+I adattato alle varianti finali la quota di siti invarianti resta
  al valore di partenza dell'ottimizzazione.
- Con `katharoseq.target_sensitivity` 0,90 il filtro per profondità esclude 278
  campioni biologici su 770.

Nessuna misura esisteva per gli altri valori delle griglie.

## 2. Esito dell'analisi di sensibilità e configurazione congelata

L'analisi è descritta per intero in [sensibilita.md](sensibilita.md), con le tabelle in
[sensibilita/](sensibilita/). Qui se ne registra l'esito.

| Parametro | Valore prima | Esito della regola | Valore congelato |
|---|---|---|---|
| `katharoseq.target_sensitivity` | 0,90 | non instabile; nessuna alternativa in zona stabile | 0,90 |
| `decontam.threshold` | 0,5 | non instabile | 0,5 |
| `prev.min_fraction` | 0,01 | non instabile; nessuna alternativa in zona stabile | 0,01 |
| `decontam.mode`, `decontam.batch_combine` | `aggregate`, `minimum` | modalità corrente ammissibile | `aggregate`, `minimum` |
| `phylo.model` (predefinito; la fase è disattivata) | GTR+G+I | BIC più basso di 270,2 rispetto a GTR+G | GTR+G+I |

Nessun valore è cambiato, e il risultato finale su OSD-734 è lo stesso di prima
dell'analisi: 492 campioni, 1.753 varianti, 21.202.825 letture.

**La configurazione congelata** è `dati/osd734/config_osd734.yaml`. Dichiara in modo
esplicito i quattro valori della tabella, anche se coincidono con i predefiniti, perché
non dipendano da essi, e con loro tutti i parametri che descrivono il dataset e tutti
quelli il cui predefinito è tarato su OSD-734 (sezione 3.15); dichiara l'immagine del
container con il suo digest di registro; attiva la regola rigorosa sulla provenienza
(sezione 3.14).

**Limiti che l'analisi ha messo in evidenza**, da dichiarare con i risultati:

- L'insieme dei campioni analizzati dipende da `katharoseq.target_sensitivity` a ogni
  passo della griglia: fra valori adiacenti cambia stato dal 3% al 12% dei campioni
  biologici, e nessun valore è in una zona stabile. Sotto 0,90 una parte dell'effetto
  non viene dalla curva ma dal suo venir meno: a 0,85 e a 0,80 le piastre 3 e 4
  ripiegano sulla soglia fissa.
- L'insieme delle varianti finali dipende da `prev.min_fraction` in modo sostanziale e
  continuo (dal 7% al 20% delle varianti per ogni mezzo punto percentuale); le letture
  conservate no (meno dello 0,2% per passo). Il numero di varianti è il numero di
  quelle presenti in almeno 5 campioni su 492, non una proprietà del dataset.
- La regola non confronta fra loro le modalità di decontaminazione ammissibili: quella
  per piastra con `fisher` toglie 167 contaminanti invece di 759 e resta un'alternativa
  di metodo difendibile.
- Nel modello GTR+G+I la quota di siti invarianti resta al valore di partenza
  dell'ottimizzazione: il confronto fra i modelli vale per quella quota, non per la sua
  stima di massima verosimiglianza.

## 3. Decisioni di metodo

Ogni voce riporta la decisione, la motivazione e le misure che la sostengono. Dove non
è indicato altro, le misure sono sull'esecuzione completa di OSD-734: 960 campioni da
due corse di sequenziamento e dieci piastre da 96, letture single-end della regione V4
del gene 16S rRNA.

### 3.1 Dati di ingresso e classi dei campioni

**Nessuna rimozione dei primer.**
- Motivazione: le letture depositate non contengono il primer; iniziano nella regione
  conservata a valle del 515F. Tagliare una lunghezza fissa all'inizio toglierebbe
  sequenza biologica.
- Misure: il motivo conservato atteso dopo il primer è all'inizio delle letture nella
  grande maggioranza dei campioni, compresi i controlli negativi (dal 79% al 95%).
  Misura la presenza della regione amplificata, non la qualità del campione: non è
  usato come criterio di esclusione.

**Tre classi di campioni, riconosciute dai metadati.**
- Decisione: campioni biologici, controlli positivi e controlli negativi si
  riconoscono dal tipo di materiale dichiarato (`ctrl.column` e le tre liste di
  etichette); le tre classi sono disgiunte e ogni controllo di qualità le distingue.
- Motivazione: in un ambiente a bassa biomassa un controllo negativo con zero letture è
  il risultato atteso, non un guasto; una soglia calcolata su tutti i campioni insieme
  scambierebbe l'uno per l'altro.

**I tamponi mai aperti sono controlli negativi.**
- Decisione: i 33 campioni la cui posizione è "Unopened 3DMM Swab Tube", dichiarati
  "Surface swab" nei metadati, sono trattati come controlli negativi
  (`ctrl.blank_override_column`, `ctrl.blank_override_values`). Il materiale dichiarato
  resta quello originale nei metadati. Le classi diventano 770 biologici, 80 positivi,
  110 negativi (nei metadati 803, 80, 77).
- Motivazione: non hanno campionato alcuna superficie; la loro composizione è quella
  dei reagenti.
- Misure: mediana di 2.427 letture, contro 28.829 dei controlli negativi dichiarati e
  20.832 dei biologici. Lasciati fra i biologici, la decontaminazione toglieva il 27,7%
  delle loro letture contro il 5,1% degli altri biologici, e Variovorax, il ceppo dei
  controlli positivi, risultava contaminante per contaminazione crociata nei 77
  negativi dichiarati; con la riclassificazione non lo è più (p = 0,88).

**La chiave dei campioni è l'accession dell'esperimento.**
- Motivazione: due campioni risequenziati (LAB1P3.L1, NOD2S4.R6) hanno lo stesso nome
  biologico nel file che associa i campioni alle piastre, e si distinguono solo per
  accession.

**L'appartenenza alle piastre viene da un file esterno ai metadati depositati.**
- Decisione: piastra e pozzetto di ogni campione vengono dal file degli autori dello
  studio, ricostruito da una revisione fissata del loro repository; il file non è
  ridistribuito perché il repository di origine non dichiara una licenza.
- Misure: la piastra 10 ha 83 campioni biologici e 5 controlli negativi dichiarati,
  invece di 80 e 8: è la composizione reale della piastra.

### 3.2 Controlli di qualità per classe

**Le soglie di arresto si calcolano sui campioni biologici e sui controlli positivi,
mai sui negativi.**
- Decisione: la profondità minima all'ingresso (sulla mediana della classe), i
  campioni azzerati dal filtro e la perdita media di letture, la frazione di letture
  chimeriche e la quota di varianti senza phylum si valutano sulle due classi in cui
  ci si aspetta DNA amplificabile. I controlli negativi sono rendicontati a parte.
- Motivazione: un controllo negativo che perde tutte le letture nel filtro attesta
  l'assenza di contaminazione; contarlo fra i campioni azzerati fermerebbe
  un'esecuzione corretta.

**Un avviso e un arresto non condividono mai lo stesso codice.**
- Motivazione: chi legge il registro di un'esecuzione deve poter distinguere dal solo
  codice se la catena si è fermata o ha proseguito segnalando.

### 3.3 Filtro di qualità e modello d'errore

**Troncamento a 137 basi (`filter.truncLen`).**
- Motivazione: il modello d'errore di DADA2 richiede letture allineate per posizione e
  della stessa lunghezza; le letture del dataset sono lunghe da 137 a 151 basi, e 137 è
  la lunghezza che non ne scarta nessuna per lunghezza insufficiente.

**Un modello d'errore per corsa di sequenziamento (`err.batch_column`).**
- Motivazione: i tassi d'errore dipendono dalla corsa; un modello unico per due corse
  ne descriverebbe la media.

**Le basi usate per la stima (`err.nbases`) possono essere raddoppiate da un nuovo
tentativo automatico se la stima non converge.**
- Motivazione: raddoppiarle cambia la quantità di dati della stima, non il metodo.
- Misure: da 1e8 a 2e8 basi i tassi d'errore cambiano di un fattore mediano di 10^0,014
  e 10^0,037 nelle due corse; le varianti passano da 13.130 a 13.135, con 13.129 in
  comune, e quelle non comuni raccolgono 25 e 94 letture su 31 milioni; la dissimilarità
  di Bray-Curtis per campione fra le due tabelle ha mediana 0 e massimo 0,004. L'effetto
  non è nullo: il parametro resta fra quelli che determinano il risultato, e il valore
  effettivamente usato è registrato nel manifesto della fase.

### 3.4 Inferenza delle varianti

**Pseudo-pooling, eseguito in due passate su lotti di campioni.**
- Decisione: le varianti si inferiscono campione per campione, poi una seconda passata
  usa come informazione a priori le varianti trovate in almeno due campioni; i campioni
  si elaborano a lotti (`run.batch_size`) per contenere la memoria.
- Motivazione: il pseudo-pooling recupera le varianti rare condivise fra campioni, che
  in un ambiente a bassa biomassa sono una parte rilevante del segnale, senza il costo
  di memoria dell'inferenza congiunta su 960 campioni.
- Misure: gli artefatti sono identici byte per byte con lotti di 24 e di 5 campioni, e
  l'elaborazione a lotti dà esattamente le varianti di `dada(pool = "pseudo")` in una
  sola chiamata. Per questo la dimensione del lotto non entra fra i parametri che
  determinano il risultato, e dimezzarla dopo un esaurimento di memoria è un nuovo
  tentativo legittimo.

### 3.5 Chimere

**La frazione chimerica si misura sulle letture, non sulle varianti.**
- Motivazione: le chimere sono molte come sequenze distinte e poche come letture; la
  frazione sulle varianti sovrastima il problema e dipende dalla profondità.
- Misure: letture chimeriche 0,44% nei biologici (massimo per campione 4,4%), 0,02% nei
  positivi, 0,03% nei negativi; sulle varianti sarebbero 8,7%, 2,3% e 2,8%. Le soglie
  di avviso (0,25) e di arresto (0,50) sono lontane in ogni classe.

### 3.6 Tassonomia

**Classificatore bayesiano ingenuo (`assignTaxonomy` di DADA2) su SILVA 138, con
bootstrap minimo 50, senza assegnazione della specie, con prova del filamento
complementare.**
- Motivazione: l'insieme di addestramento di SILVA 138 per il classificatore IdTaxa non
  è più distribuito da alcuna fonte; la specie non si assegna perché le letture coprono
  137 delle circa 253 basi dell'amplicone.

**DADA2 è usato nella versione 1.36.0 con una correzione che rende riproducibile la
risoluzione dei pareggi fra generi.**
- Motivazione: nella versione ufficiale, quando due generi hanno la stessa probabilità
  la scelta usa un generatore di numeri casuali che `set.seed` non controlla: due
  esecuzioni identiche danno tassonomie diverse. La correzione deriva il seme dei
  pareggi dal seme dell'esecuzione (`run.seed`), dall'indice della sequenza e dal caso;
  la scelta resta uniforme fra i generi a pari probabilità.
- Misure: sulle 12.045 varianti del dataset, dieci esecuzioni della versione ufficiale
  differiscono fra loro in 160 varianti; la versione corretta dà lo stesso risultato in
  tre esecuzioni con numeri diversi di thread, e differisce dalle ufficiali solo in
  quelle 160 varianti, tutte con almeno un pareggio. I pareggi riguardano 820 varianti
  su 12.045 (143.775 letture), 120 nella classificazione principale.
- La versione modificata è registrata nel file di blocco dei pacchetti R con il file
  della correzione e la sua impronta, e la costruzione dell'immagine la accetta solo se
  coincidono.

**Il riferimento tassonomico è verificato per impronta a ogni avvio.**
- Motivazione: un riferimento sostituito fra una ripresa e l'altra cambierebbe la
  tassonomia senza lasciare traccia.

### 3.7 Oggetto integrato

**Tutti i campioni dell'inventario entrano nell'oggetto, anche quelli senza letture
dopo il filtro, con conteggi a zero.**
- Motivazione: un campione che scompare non è distinguibile da un campione mai
  sequenziato.

**Le varianti sono numerate (ASV1, ASV2, ...) per letture decrescenti nei soli campioni
biologici, poi per letture totali, poi per sequenza.**
- Motivazione: la numerazione deve riflettere l'abbondanza nei campioni che si
  studiano, ed essere determinata dai dati senza ambiguità.
- Misure: numerando su tutti i campioni, la terza e la quarta variante sarebbero il
  cloroplasto dei controlli negativi e Variovorax dei controlli positivi (159esimo fra
  i biologici); 828 varianti non hanno letture nei biologici.

### 3.8 Controlli positivi e soglia di profondità

**La profondità minima dei campioni si legge su una curva KatharoSeq stimata per
piastra.**
- Decisione: per ogni piastra, la fedeltà dei controlli positivi (frazione delle
  letture assegnate al ceppo atteso, su otto livelli di diluizione) si descrive in
  funzione del logaritmo della profondità con la curva f = 1 / (1 + (x50 / x)^h); la
  soglia è la profondità a cui la curva raggiunge `katharoseq.target_sensitivity`.
  La stima è deterministica (nessun numero casuale). Si usa il modello per piastra o
  quello aggregato secondo la bontà di adattamento (R^2 sulla fedeltà).
- Misure: R^2 per piastra 0,92 contro 0,57 dell'aggregato, che resta sotto il minimo
  richiesto (`katharoseq.min_r2`, 0,8). Soglie fra 6.376 e 41.968 letture nelle otto
  piastre con curva valida.

**Una curva è valida solo se la soglia è determinata dai dati.**
- Decisione: oltre al numero minimo di controlli e all'R^2, si richiede almeno un
  controllo osservato fra il punto medio della curva e la soglia. Una piastra senza
  curva valida ripiega su una soglia fissa (`qc.min_reads_raw`, 1.000 letture grezze) e
  l'esecuzione lo segnala.
- Motivazione: se fra il punto medio e la soglia non c'è alcuna osservazione, la soglia
  è un'estrapolazione della forma della curva, non una misura.
- Misure: ripiegano le piastre 1 e 2. Nella piastra 1 la curva è un gradino fra 20.257
  e 98.611 letture senza osservazioni intermedie; nella 2 l'R^2 è 0,53. In entrambe un
  contaminante cloroplastico fa il 69,5% delle letture dei controlli negativi e domina
  i controlli positivi più diluiti. Un intervallo di confidenza della soglia è stato
  valutato come criterio alternativo e scartato: risultava stretto proprio dove i dati
  non dicono nulla (piastra 1) e avrebbe respinto piastre sostenute da osservazioni.

**Sensibilità richiesta 0,90 (`katharoseq.target_sensitivity`).**
- Motivazione: è il valore del piano del progetto, più severo dello 0,80 usato dagli
  autori del metodo. L'analisi di sensibilità (sezione 2) lo mantiene e ne dichiara il
  limite.
- Misure: con 0,90 il filtro per profondità esclude 278 campioni biologici su 770; con
  0,80 ne esclude 182, ma due piastre in più perdono la curva.

**La curva si stima sulle letture senza chimere, compresi i taxa che i filtri finali
escluderanno.**
- Misure: escludendo cloroplasti, mitocondri ed eucarioti dalla stima, le piastre 1 e 2
  ripiegano comunque e le altre soglie cambiano di poco (da 5.885 a 41.288 contro da
  6.376 a 41.968): miglioramento non netto, a fronte di una dipendenza in più.

**La conformità dei controlli positivi si giudica per livello di diluizione, e non
ferma l'esecuzione.**
- Decisione: un controllo è non conforme se la sua fedeltà o la sua profondità si
  discostano da quelle degli altri controlli dello stesso livello (punteggio z
  modificato oltre 3,5); i non conformi escono dalla stima della curva. La non
  conformità produce un avviso (`ctrl.positive_gate: false`).
- Motivazione: ai livelli più diluiti i contaminanti di reagente prevalgono sul ceppo
  atteso per costruzione: giudicare ogni controllo sulla dominanza del ceppo
  respingerebbe controlli che si comportano come previsto.
- Misure: conformi 74 controlli su 80; per dominanza del ceppo sarebbero 47 su 80.

### 3.9 Decontaminazione

**Metodo della prevalenza di decontam, soglia 0,5, controlli positivi esclusi dal
confronto.**
- Decisione: una variante è contaminante se la sua prevalenza nei controlli negativi,
  rispetto ai campioni biologici, dà una probabilità sotto `decontam.threshold`. La
  soglia 0,5 classifica contaminante ciò che è più prevalente nei negativi che nei
  campioni.
- Motivazione: il metodo della frequenza richiede la concentrazione del DNA, che i
  metadati non riportano; i controlli positivi non sono né campioni né bianchi.
- Misure: 759 contaminanti, che raccolgono il 4,9% delle letture dei biologici, l'87,1%
  di quelle dei negativi e il 9,7% di quelle dei positivi. Il contaminante cloroplastico
  delle piastre 1 e 2 è riconosciuto. La variante più abbondante del dataset
  (Pseudomonas), presente in 735 biologici su 770 e in 81 negativi su 110, non è
  classificata contaminante.

**Confronto aggregato su tutti i controlli negativi (`decontam.mode: aggregate`), non
per piastra.**
- Motivazione: per piastra ogni confronto dispone di 6-14 controlli negativi, e la
  combinazione predefinita (`minimum`) prende la più bassa di dieci probabilità stimate
  ciascuna su pochi controlli: è permissiva per costruzione.
- Misure: per piastra con `minimum`, 1.138 contaminanti e il 44,2% delle letture dei
  biologici, comprese varianti di generi attesi sulle superfici abitate
  (Staphylococcus, Corynebacterium, Streptococcus); con `fisher`, 167 contaminanti e il
  2,1%. La modalità non dichiarata si calcola comunque a ogni esecuzione, come
  diagnostica.

**Un limite alla frazione di letture rimosse (`qc.max_frac_contaminant`, 0,40) ferma
l'esecuzione.**
- Motivazione: una decontaminazione che toglie quasi metà delle letture dei campioni
  indica un confronto mal posto, non un dataset molto contaminato; va esaminata prima
  di proseguire.

**Almeno 5 controlli negativi per un confronto (`decontam.min_blanks`).**
- Motivazione: sotto questo numero la prevalenza nei negativi non è stimabile; una
  piastra che non li raggiunge si confronta con i negativi di tutte le piastre.

### 3.10 Filtri finali

**Ordine dei filtri: profondità, tassonomico, prevalenza, letture finali. La
decontaminazione li precede.**
- Motivazione: la prevalenza va calcolata sui campioni che si terranno e dopo aver
  tolto contaminanti e taxa non batterici, altrimenti una variante sarebbe giudicata
  su campioni e conteggi che non fanno parte del risultato.

**Il filtro per profondità confronta la soglia della piastra con le letture del
campione allo stadio a cui la soglia è stata stimata.**
- Motivazione: la soglia viene dalla curva dei controlli positivi, stimata sulle
  letture senza chimere; confrontarla con le letture dopo la decontaminazione
  mescolerebbe due grandezze.

**Filtro tassonomico: via le varianti senza phylum e quelle di cloroplasti, mitocondri
ed eucarioti.**
- Misure: in quattro campioni della stessa posizione la variante dominante è 16S
  mitocondriale (dal 52% all'87% delle letture). Un campione svuotato da questo filtro
  è escluso con una segnalazione, senza fermare l'esecuzione: è un fatto biologico.

**Prevalenza minima 1% dei campioni conservati (`prev.min_fraction`), con numeratore e
denominatore sullo stesso insieme.**
- Decisione: il minimo è `ceil(prev.min_fraction x campioni conservati)`, calcolato sui
  campioni biologici rimasti dopo i filtri per profondità e tassonomico.
- Motivazione: l'1% e non il 5% consueto perché i campioni vengono da oltre cento
  posizioni distinte della stazione: una variante caratteristica di una posizione
  compare in pochi campioni per disegno. Il denominatore sono i campioni conservati
  perché quelli esclusi per profondità hanno una composizione non attendibile.
- Misure: 492 campioni conservati, minimo 5 campioni, 1.753 varianti. Sul totale dei
  770 biologici il minimo sarebbe 8 e le varianti 1.439. Con il 5% le varianti sarebbero
  792 (sezione 2).

**Due soglie sulle letture finali, con ruoli diversi.**
- Decisione: un campione sotto `qc.min_reads_final` letture dopo i filtri è escluso,
  senza arresto; se l'insieme dei campioni finali trattiene meno di
  `qc.min_frac_reads_retained` (0,40) delle letture senza chimere, l'esecuzione si
  ferma.
- Misure: nessun campione sotto la soglia per campione (minimo 1.085 letture);
  frazione trattenuta dall'insieme 0,922. La frazione non si applica per campione:
  otto campioni sarebbero sotto 0,40 per il 16S mitocondriale tolto dal filtro
  tassonomico, che non è un difetto.

### 3.11 Oggetto finale

**L'oggetto finale contiene i soli campioni biologici; i controlli sono consegnati in
un oggetto separato.**
- Motivazione: le analisi ecologiche non devono poter includere un controllo per
  errore; i controlli restano disponibili per la verifica.

**Accanto all'oggetto R si consegnano tabelle di testo da cui l'oggetto è
ricostruibile, con le impronte di ogni file.**
- Motivazione: il risultato deve restare leggibile senza R e verificabile senza la
  pipeline. La ricostruzione dalle tabelle è verificata a ogni esecuzione.

### 3.12 Filogenesi

**Fase opzionale, disattivata per OSD-734, eseguita sulle varianti finali.**
- Motivazione: su letture di 137 basi di una sola regione del gene un albero costruito
  da zero è debolmente risolto; costruirlo sulle varianti finali evita di calcolarlo
  su varianti che i filtri toglieranno, e lo rende coerente con l'oggetto consegnato.

**Allineamento con DECIPHER; massima verosimiglianza con phangorn, modello GTR+G+I;
ricerca deterministica; radice al punto medio.**
- Motivazione: la ricerca per scambi fra rami vicini a partire dall'albero
  neighbor-joining non usa numeri casuali, e dà lo stesso albero a ogni esecuzione e
  con qualunque numero di thread; le ricerche stocastiche possono trovare una
  verosimiglianza poco più alta ma non sono riproducibili allo stesso modo. Fra le
  varianti non c'è un gruppo esterno su cui radicare: il punto medio dipende solo
  dall'albero.
- Misure: confronto dei modelli nella sezione 2 e in sensibilita.md, sezione 4.

### 3.13 Errori e nuovi tentativi

**Ogni condizione di arresto o di avviso ha un codice in un catalogo unico, con la
descrizione e l'azione da compiere.**

**Un nuovo tentativo automatico è ammesso solo dove l'azione correttiva non cambia
alcuna assunzione di metodo.**
- Decisione: quattro casi. Errore transitorio di lettura degli archivi (si ritenta
  senza modifiche); stima del modello d'errore non convergente (si raddoppiano le basi
  usate, sezione 3.3); memoria esaurita nell'inferenza (si dimezza il lotto, sezione
  3.4); memoria esaurita nella costruzione della tabella, per il quale nessuna azione
  correttiva esiste e il tentativo è dichiarato inutile. Ogni altro arresto richiede
  una persona.
- Motivazione: un tentativo che cambiasse una soglia o un metodo per far concludere la
  catena consegnerebbe un risultato diverso da quello dichiarato. Ogni aggiustamento
  è annotato nel manifesto della fase.

### 3.14 Riproducibilità e provenienza

**Gli artefatti non contengono date, durate né percorsi.**
- Motivazione: due esecuzioni della stessa configurazione devono dare file identici
  byte per byte, perché il confronto delle impronte sia il criterio di riproduzione.
- Misure: due esecuzioni complete indipendenti, una delle quali a partire dai dati
  scaricati da zero dagli archivi pubblici, danno le stesse impronte per tutti i 1.031
  artefatti. Le impronte attese sono pubblicate in `dati/osd734/`.

**Ogni fase dichiara i parametri da cui dipende, e si rifà solo quando cambiano quelli,
il suo codice o ciò che riceve dalle fasi a monte.**
- Decisione: i parametri che governano le risorse e non il calcolo (numero di thread,
  dimensione del lotto, cartella di uscita, conservazione degli intermedi, nuovi
  tentativi) non rendono invalida alcuna fase. Nel dubbio un parametro si considera
  influente.
- Misure: la tassonomia è identica con 4 e con 12 thread, l'albero filogenetico con
  numeri diversi di thread, gli artefatti dell'inferenza con lotti di 24 e di 5 campioni.

**Ambiente di calcolo fissato: immagine del container identificata per digest, versioni
dei pacchetti R in un file di blocco verificato alla costruzione.**
- Motivazione: versioni diverse di R e di Bioconductor possono dare risultati numerici
  diversi; un'etichetta di immagine può essere riassegnata, un digest no.

**Scritture atomiche; ordinamenti e codifiche indipendenti dalle impostazioni locali.**
- Misure: con le impostazioni locali del sistema cambiate, gli artefatti restano
  identici; prima che la scrittura degli oggetti R fosse resa indipendente dalla
  codifica della sessione, ne cambiavano 12 su 95.

**Regola rigorosa sulla provenienza (`run.strict_provenance`), attiva nella
configurazione congelata.**
- Decisione: con la regola attiva l'esecuzione parte solo se il codice è identificato
  da un commit e l'ambiente R corrisponde al file di blocco; l'impronta del codice di
  ogni fase, quella del file di blocco e l'immagine dichiarata entrano fra le
  condizioni di validità della fase, così che cambiarle la rende da rifare e non solo
  da segnalare. La regola è disattivata per difetto, perché durante lo sviluppo il
  codice ha modifiche non committate.
- Che cosa è **verificato** a ogni avvio:
  - il repository git è leggibile, e i file di `src/` e `R/` non hanno modifiche
    rispetto al commit né file non tracciati; gli script R eseguiti sono quelli del
    repository. Se git non è leggibile l'esecuzione è rifiutata (`E-PROV-01`): senza
    commit nulla identifica il codice. Con modifiche non committate è rifiutata
    (`E-PROV-02`);
  - la versione di R e quelle di tutti i pacchetti del file di blocco, compresa la
    correzione di DADA2 con la sua impronta, sono lette dalle librerie installate e
    confrontate con il file (`E-PROV-03` se differiscono). Sono le librerie che i
    calcoli caricheranno.
- Che cosa è **dichiarato** e non verificabile: l'immagine (`run.container`). Da dentro
  un container il digest dell'immagine in esecuzione non è conoscibile: il valore è una
  dichiarazione di chi lancia l'esecuzione, registrata e vincolante per la validità
  delle fasi, ma sostenuta dalla verifica dell'ambiente R e non da una misura propria.
  Non sono verificate le librerie di sistema dell'immagine, che il file di blocco non
  descrive, né il codice Python di terze parti, le cui versioni sono fissate nei file
  dei requisiti e installate alla costruzione dell'immagine.
- Motivazione: una configurazione congelata certifica i risultati solo se ciò che li
  ha calcolati è identificato; ciò che non si può verificare deve fermare l'esecuzione
  o essere dichiarato come tale, non passare in silenzio.

### 3.15 Generalità: che cosa descrive il dataset e che cosa la pipeline

**I parametri che descrivono il dataset non hanno un valore predefinito.**
- Decisione: sono obbligatori i parametri sul formato dei metadati (colonne del nome
  del campione e del file, colonna e etichette delle classi, colonne della posizione,
  del modulo, della piastra e della corsa, colonne da portare nell'oggetto), quelli
  sull'esperimento (lunghezza di troncamento, primer, motivo conservato, taxon atteso
  nei controlli positivi, nome e versione del riferimento) e l'espressione che estrae
  dal nome dei file la chiave del campione: 25 in tutto. Una configurazione che ne
  omette uno non parte, e l'errore li elenca tutti. Un parametro che per il dataset
  non ha contenuto si dichiara nullo o vuoto: un dataset senza controlli positivi lo
  dice con un elenco vuoto.
- Motivazione: un valore ereditato da un altro dataset (una colonna che per caso
  esiste, un'etichetta, un troncamento) non produce un errore ma un risultato
  plausibile e sbagliato. Dichiarare è l'unico modo di sapere che chi esegue ha
  guardato i propri dati.
- Misure: su un secondo dataset (OSD-276, 15 tamponi della stessa stazione, altro
  laboratorio) la pipeline con i predefiniti di OSD-734 respingeva i file per il
  formato dell'accession, diagnosticava un'etichetta non dichiarata come primer nelle
  letture, e si fermava per una colonna inesistente dopo dodici minuti di calcolo; il
  troncamento ereditato toglieva 14 basi a ogni lettura senza alcun arresto.

**Otto parametri conservano un predefinito tarato sul dataset di riferimento, e il
report di ogni esecuzione elenca quelli presi per difetto.**
- Decisione: taglio iniziale, bootstrap minimo e assegnazione della specie, metodo,
  soglia e modalità della decontaminazione, frazione minima del motivo conservato,
  prevalenza minima. Sono scelte di metodo, non descrizioni del dataset: hanno un
  valore di partenza sensato, ma il fatto che lo giustifica va riesaminato.
- Motivazione: renderli tutti obbligatori chiederebbe a chi esegue di scegliere prima
  di aver visto un risultato; lasciarli impliciti nasconderebbe che sono stati scelti
  su altri dati. La segnalazione scompare quando il parametro viene dichiarato.
- Misure: con 15 campioni la prevalenza minima dell'1% equivale a un solo campione,
  cioè a nessun filtro.

**La validazione iniziale verifica tutto ciò che la configurazione dichiara, nell'ordine
in cui le verifiche dipendono l'una dall'altra.**
- Decisione: ogni colonna nominata si cerca nelle intestazioni prima del calcolo; il
  file del lotto deve avere una e una sola riga per campione; le classi dei campioni
  si verificano prima dei controlli sulle letture; il primer si cerca in testa alle
  letture solo se il filtro non lo toglierà, e il motivo conservato dove inizierà la
  lettura filtrata; un dataset senza controlli positivi, o con pochi controlli
  negativi, viene dichiarato tale subito.
- Motivazione: la diagnosi deve nominare l'errore vero. Un controllo eseguito prima di
  ciò da cui dipende fallisce con il messaggio di un altro problema.

**Ciò che la validazione iniziale lasciava passare in silenzio ha un esito dichiarato.**
- Decisione: un parametro che ha effetto solo insieme a un altro, dichiarato da solo,
  è respinto (`E-G15-14`: i valori della riclassificazione in controllo negativo
  senza la colonna, la colonna del nome nella tabella di studio senza la tabella,
  l'espressione del modulo o le posizioni non di superficie senza la colonna della
  posizione). Senza colonna della
  posizione il modulo viene dal file del lotto, se lo dichiara. Un valore dei
  metadati con una tabulazione o un a capo dentro un campo fra virgolette, in una
  colonna che la pipeline legge e in una riga di un campione dell'assay, ferma la
  validazione (`E-S0-02` per le tabelle di
  assay e di studio, `E-S0-08` per il file del lotto). Le righe del file del lotto
  che non corrispondono ad alcun campione non fermano, ma sono dichiarate con il
  loro numero (`E-S0-19`). Un file che nelle letture ispezionate porta le due
  letture di ogni coppia (intestazioni con `/1` e `/2`, o la stessa intestazione
  ripetuta), o le sole seconde letture, è respinto come layout non single-end
  (`E-S0-07`). Un record anomalo isolato non ferma un dataset single-end: i
  marcatori di ciascuna delle due letture devono comparire in almeno il 5% delle
  letture ispezionate (e in almeno due), e così i nomi che compaiono due volte, che
  si contano per nome; un file è di sole seconde
  letture se le marcate come seconde sono almeno il 95% e le prime meno del 5%. Un
  nome che compare più di due volte (un'intestazione vuota o uguale per tutte le
  letture) non è un segno di coppia. Il primer in testa alle letture
  si cerca nei soli campioni biologici e controlli positivi, come il motivo
  conservato.
- Motivazione: sono tutti casi in cui la pipeline proseguiva con un risultato
  plausibile (un modulo mancante, una riclassificazione non avvenuta, il doppio delle
  letture) o si fermava molto dopo con un errore che non nominava la causa. Per il
  primer: in un controllo negativo con tre letture una sola che comincia come il
  primer supera qualunque soglia in frazione.
- Limite dichiarato: i file con le due letture di ogni coppia si riconoscono nelle
  sole letture ispezionate (`qc.head_reads`). Un file con tutte le prime letture
  seguite da tutte le seconde non si riconosce se le letture ispezionate non
  arrivano al secondo blocco. Dei 15 file che ENA distribuisce per il secondo
  dataset (da 35.489 a 105.423 coppie per file) i 7 che cominciano con le seconde
  letture sono respinti con il valore predefinito di 10.000, e bastano a fermare la
  validazione; gli 8 che cominciano con le prime, presi da soli, passerebbero.
  Chiudere il limite richiede di contare i marcatori su tutte le letture, nel
  profilo di qualità. Non si riconoscono nemmeno le coppie marcate in altri modi
  (suffissi `.1` e `.2` o `_1` e `_2`, commento separato da una tabulazione o da
  più spazi, letture rinumerate senza marcatore), né, fra quelle riconoscibili dal
  solo nome, i file in cui anche un solo nome compare più di due volte. E un dataset single-end fatto di seconde
  letture, se le intestazioni lo dichiarano, è respinto: non c'è un parametro per
  dichiararlo voluto.

**Il numero di thread non fa parte dell'identità della configurazione.**
- Decisione: `run.threads` è nullo per difetto, cioè automatico: i processori
  utilizzabili si contano quando servono, e il numero non entra nella configurazione
  registrata né nel suo digest. Un valore dichiarato resta dichiarato, e la
  validazione lo confronta con i processori utilizzabili.
- Motivazione: il digest identifica la configurazione; con il numero dei processori
  al suo interno la stessa configurazione avrebbe avuto un digest diverso su ogni
  macchina, e una ripresa altrove avrebbe registrato una configurazione nuova.
- Misure: la stessa configurazione dà lo stesso digest con 4 e con 16 processori.

### 3.16 Generalità delle fasi di calcolo: ciò che un dataset può non avere

**Il troncamento si giudica sulla frazione di letture più corte, non sulla più corta.**
- Decisione: G09 (sulle prime letture) e S1 (su tutte) si fermano solo se le letture
  più corte di `filter.truncLen` superano `qc.max_frac_short_reads` (0,05) delle
  letture dei campioni biologici, o di quelle dei controlli positivi: ogni classe si
  giudica da sola, perché in una frazione unica pesata sulle letture la classe meno
  numerosa (tre controlli positivi fra cento campioni) potrebbe perdere tutte le sue
  letture restando sotto la soglia. I controlli negativi non contano. `filter.minLen`, che coincideva per costruzione con `filter.truncLen`, e il
  controllo che lo sorvegliava (`E-G15-01`) sono rimossi. Il report indica un
  troncamento suggerito: il minore fra il più lungo che scarta non oltre il 5% delle
  letture e l'ultima posizione con qualità mediana dei biologici almeno 30. È
  un'indicazione: il valore applicato resta quello dichiarato.
- Motivazione: in un dataset a lunghezza variabile una lettura corta c'è sempre, e il
  filtro la toglie senza danno; fermarsi per quella rendeva la pipeline inutilizzabile
  fuori dal dataset di riferimento, le cui letture sono tutte più lunghe del
  troncamento. Un bianco amplifica poco e male, spesso solo dimeri: non dice nulla sul
  troncamento adatto ai campioni.
- Misure: sul dataset di riferimento nessuna lettura è più corta di 137 basi, e la
  frazione è zero.

**Il modello di errore ha una variante per le qualità raggruppate, e S1 dice quando
serve.**
- Decisione: `err.error_function` vale `loess` (la funzione standard di dada2,
  predefinita) o `loess_monotono`: un loess di primo grado con span 2 e pesi pari al
  logaritmo in base 10 delle basi osservate a ogni qualità (span e pesi sono quelli
  discussi dalla comunità di dada2 per i dati a qualità raggruppate), poi reso non
  crescente con la qualità (a ogni qualità il tasso è almeno quello di tutte le
  qualità superiori); con due soli valori di qualità la retta per i due punti. S1
  conta i valori di qualità distinti e con quattro o meno lo dichiara (`E-S1-03`).
  Una corsa senza letture filtrate ferma S3 con un codice proprio (`E-S3-03`);
  letture con un solo valore di qualità, da cui nessuna funzione può stimare il
  modello, la fermano con `E-S3-04`, e solo quelle: ogni altra stima fallita ha un
  codice diverso (`E-S3-05`), deciso contando i valori di qualità delle letture su
  cui si stimava.
- Motivazione: con pochi valori di qualità la stima standard può dare un tasso di
  errore che cresce con la qualità, e l'inferenza tratterebbe come più affidabili le
  basi peggiori. La scelta resta di chi conduce l'analisi: la pipeline segnala, non
  sostituisce il modello da sé. Lo span largo e il primo grado servono proprio ai
  pochi punti: con lo span stretto la curva interpolava esattamente tre punti o non
  era adattabile.
- Misure: su matrici di transizione con due, tre e quattro qualità osservate la
  variante dà tassi finiti e non crescenti; la funzione standard si ferma con due.

**L'assenza di una classe di controlli è una condizione dichiarata, non un errore.**
- Decisione: senza controlli positivi, o senza la colonna dei livelli, S11 non adatta
  alcuna curva e lo dichiara (`E-S11-05`), distinto dal ripiego di una curva non
  attendibile (`E-S11-02`); un controllo senza piastra entra nella sola curva
  aggregata e uno con una sola lettura in nessuna. Con meno di `decontam.min_blanks`
  controlli negativi con letture S12 non toglie nulla e lo dichiara (`E-S12-03`).
  Se nessun campione biologico ha letture la fase si ferma (`E-S12-04`), e così se
  tutte le varianti risultano contaminanti entro la quota ammessa (`E-S12-05`).
  La colonna dei livelli si cerca nel file del lotto e poi nella tabella di studio, ed
  entra nell'oggetto anche se non è fra le colonne richieste. Il valore del ripiego
  della soglia di profondità non cambia.
- Motivazione: un ripiego silenzioso (nessuna decontaminazione, nessuna curva) dà un
  risultato plausibile e non confrontabile con uno decontaminato e calibrato; un
  errore di R senza codice non dice che cosa manca.

**Nei filtri finali un campione esce con il motivo; la cartella consegnata contiene
solo ciò che si consegna.**
- Decisione: un campione svuotato dal filtro di prevalenza esce, e la fase lo
  dichiara (`E-S13-02`, da arresto a degradazione); se nessun campione supera i
  filtri la fase si ferma (`E-S13-04`). I nomi dei taxa si confrontano senza il
  prefisso di rango (`p__`, `o__`), e una tassonomia senza il rango Phylum è
  dichiarata (`E-S13-05`). `tax.assign_species` vero è respinto dallo schema, perché
  non è realizzato. S13 scrive in `12_final/intermedi/`; S14 consegna in `12_final/`
  l'oggetto finale, gli export e i controlli (`ps_controlli.rds`, solo se ci sono
  controlli), toglie prima i file di un'esecuzione precedente, e `checksum.sha256`
  li elenca tutti.
- Motivazione: un solo campione di sole varianti rare non deve fermare un dataset; un
  intermedio o un file di una configurazione precedente nella cartella consegnata si
  scambia per un risultato che nessun checksum copre.
- Misure: sul dataset di riferimento l'oggetto finale, gli export e le tabelle di
  calcolo restano identici byte per byte; cambiano la collocazione degli intermedi di
  S13, la fase che scrive `ps_controlli.rds` e l'elenco di `checksum.sha256`.
