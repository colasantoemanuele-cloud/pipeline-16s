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
