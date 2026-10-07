# Analisi di sensibilità dei parametri sul dataset OSD-734

Questo rapporto riporta come cambia il risultato finale della pipeline al variare di
quattro scelte di metodo, e applica ai numeri la regola di decisione fissata in
[decision_log.md](decision_log.md), sezione 1, prima che i numeri esistessero. Le
tabelle complete sono in [sensibilita/](sensibilita/): `misure.tsv` (una riga per
configurazione), `passi.tsv` (i cambiamenti fra valori adiacenti), `decisioni.json`
(l'applicazione della regola) e `confronto_modelli_filogenesi.json`.

## 1. Che cosa è stato fatto

**Punto di partenza.** L'esecuzione completa della pipeline su OSD-734 (960 campioni:
770 biologici, 80 controlli positivi, 110 controlli negativi) con la configurazione
`dati/osd734/config_osd734.yaml`: oggetto finale di 492 campioni per 1.753 varianti,
21.202.825 letture.

**Parametri e intervalli.** Si varia un parametro alla volta, tenendo gli altri al
valore corrente.

| Parametro | Che cosa governa | Valore corrente | Valori esaminati |
|---|---|---|---|
| `katharoseq.target_sensitivity` | la fedeltà attesa dei controlli positivi a cui si legge, sulla curva KatharoSeq di ogni piastra, la profondità minima dei campioni | 0,90 | 0,70; 0,75; 0,80; 0,85; 0,90; 0,95 |
| `decontam.threshold` | la soglia di probabilità sotto la quale una variante è classificata contaminante | 0,5 | 0,1; 0,2; 0,3; 0,4; 0,5; 0,6; 0,7 |
| `prev.min_fraction` | la frazione minima di campioni conservati in cui una variante deve comparire | 0,01 | 0,005; 0,010; 0,015; 0,020; 0,025; 0,030; in più 0,05, fuori dalla griglia di decisione |
| modalità di decontaminazione (`decontam.mode`, `decontam.batch_combine`) | confronto con tutti i controlli negativi, o per piastra con combinazione delle probabilità | aggregata | aggregata; per piastra con `minimum`; per piastra con `fisher` |

Il valore 0,80 di `katharoseq.target_sensitivity` è quello usato dagli autori del
metodo; 0,90 è quello adottato dal progetto.

**Esecuzione.** Per ogni configurazione lo strumento `scripts/sensitivity.py` copia
l'esecuzione di partenza con collegamenti fisici, cambia il parametro e riprende la
catena: la pipeline riesegue solo le fasi che dipendono dal parametro cambiato (dalla
curva dei controlli positivi, dalla decontaminazione o dai filtri finali in poi) e
riusa le altre. L'esecuzione di partenza non viene modificata: dopo le venti riprese i
suoi 1.106 file hanno gli stessi inode, le stesse dimensioni e gli stessi istanti di
modifica di prima. Tutte le riprese sono avvenute nell'immagine pubblicata della
pipeline.

**Misure.** Per ogni configurazione, sull'oggetto finale: campioni biologici
conservati e quanti cambiano stato rispetto alla configurazione corrente; varianti
finali e similarità di Jaccard del loro insieme con quello corrente; letture
conservate; soglie di profondità per piastra e piastre in cui la curva non è valida e
si ripiega sulla soglia fissa (`qc.min_reads_raw`, 1.000 letture grezze); contaminanti
rimossi e frazione delle letture dei campioni biologici che essi raccolgono.

**Come rieseguire.** Dalla radice del repository, con l'esecuzione di partenza in
`output/osd734` (si veda `dati/osd734/README.md`) e `<immagine>` l'immagine indicata
nel README:

```bash
docker run --rm \
  --memory=24g \
  --memory-swap=24g \
  -u "$(id -u):$(id -g)" \
  -e HOME=/tmp \
  -e PYTHONPATH=/app/src \
  -e AMPLICON16S_R_DIR=/app/R \
  -v "$(pwd)":/app \
  -w /app \
  <immagine> \
  python scripts/sensitivity.py \
    --config dati/osd734/config_osd734.yaml \
    --lavoro output/sensibilita \
    --uscita docs/sensibilita
```

Le varianti già calcolate si riusano; `--rifai` le ricalcola. Le riprese partono da
una copia e cambiano un parametro: con la regola rigorosa sulla provenienza attiva
nella configurazione, lo strumento va eseguito da un clone senza modifiche non
committate, come ogni altra esecuzione. Ogni ripresa dura circa tre minuti con dodici
processori; le copie occupano in tutto 2,3 GB.

## 2. Risultati

Nelle tabelle, "campioni che cambiano" conta i campioni biologici conservati in una
configurazione ed esclusi nell'altra, rispetto alla configurazione corrente; "rispetto
alla corrente" è il rapporto fra le letture conservate.

### 2.1 `katharoseq.target_sensitivity`

| Valore | Campioni conservati | Campioni che cambiano | Varianti finali | Jaccard con la corrente | Letture conservate | Rispetto alla corrente | Piastre con ripiego |
|---|---|---|---|---|---|---|---|
| 0,70 | 697 | 205 | 1.540 | 0,865 | 22.663.111 | 1,069 | 1, 2, 3, 4, 6, 7, 8, 9 |
| 0,75 | 671 | 179 | 1.619 | 0,900 | 22.576.646 | 1,065 | 1, 2, 3, 4, 6, 7, 8 |
| 0,80 | 580 | 88 | 1.676 | 0,938 | 22.048.971 | 1,040 | 1, 2, 3, 4 |
| 0,85 | 560 | 68 | 1.660 | 0,937 | 21.840.454 | 1,030 | 1, 2, 3, 4 |
| 0,90 (corrente) | 492 | 0 | 1.753 | 1,000 | 21.202.825 | 1,000 | 1, 2 |
| 0,95 | 450 | 42 | 1.709 | 0,975 | 20.116.186 | 0,949 | 1, 2 |

Soglie di profondità per piastra (letture senza chimere; 1.000 indica il ripiego, che
si applica alle letture grezze):

| Valore | P1 | P2 | P3 | P4 | P5 | P6 | P7 | P8 | P9 | P10 |
|---|---|---|---|---|---|---|---|---|---|---|
| 0,70 | 1.000 | 1.000 | 1.000 | 1.000 | 6.284 | 1.000 | 1.000 | 1.000 | 1.000 | 14.061 |
| 0,75 | 1.000 | 1.000 | 1.000 | 1.000 | 6.663 | 1.000 | 1.000 | 1.000 | 8.531 | 14.962 |
| 0,80 | 1.000 | 1.000 | 1.000 | 1.000 | 7.127 | 11.562 | 12.884 | 27.878 | 9.323 | 16.074 |
| 0,85 | 1.000 | 1.000 | 1.000 | 1.000 | 7.738 | 13.969 | 15.828 | 33.168 | 10.394 | 17.542 |
| 0,90 | 1.000 | 1.000 | 9.984 | 6.376 | 8.642 | 18.063 | 20.948 | 41.968 | 12.034 | 19.727 |
| 0,95 | 1.000 | 1.000 | 12.480 | 8.105 | 10.360 | 27.763 | 33.513 | 62.071 | 15.320 | 23.916 |

Che cosa si osserva.

- Il parametro agisce in due modi distinti. Dove la curva della piastra è valida, la
  soglia cresce con continuità con la sensibilità richiesta. Ma abbassando la
  sensibilità la soglia scende verso il punto medio della curva, e quando fra il punto
  medio e la soglia non resta alcun controllo positivo osservato la curva non è più
  ritenuta valida: la soglia non sarebbe determinata dai dati, e la piastra ripiega
  sulla soglia fissa di 1.000 letture grezze. A 0,90 ripiegano le piastre 1 e 2; a 0,85
  e a 0,80 anche la 3 e la 4; a 0,75 sette piastre; a 0,70 otto su dieci.
- Per questo i campioni conservati in più ai valori bassi non sono tutti campioni
  ammessi da una soglia più permissiva letta sulla curva: una parte consistente viene
  da piastre la cui curva ha smesso di determinare la soglia. A 0,85, 43 dei 68
  campioni in più rispetto a 0,90 sono delle piastre 3 e 4, passate dalla soglia della
  curva (9.984 e 6.376 letture) a quella fissa; a 0,80, il valore degli autori del
  metodo, sono 43 su 88.
- Il numero di varianti finali non cresce con il numero di campioni: a 0,90 le varianti
  sono 1.753 con 492 campioni, a 0,80 sono 1.676 con 580. Il filtro di prevalenza chiede
  una frazione dei campioni conservati: con più campioni il minimo sale, e i campioni
  aggiunti, poco profondi, portano poche varianti condivise.
- Le letture cambiano poco (dal 5% in meno al 7% in più), perché i campioni che entrano
  o escono sono i meno profondi.
- I contaminanti non dipendono da questo parametro: 759 in ogni configurazione.

### 2.2 `decontam.threshold`

| Valore | Campioni conservati | Varianti finali | Jaccard con la corrente | Letture conservate | Rispetto alla corrente | Contaminanti | Letture dei biologici rimosse come contaminanti |
|---|---|---|---|---|---|---|---|
| 0,1 | 492 | 1.991 | 0,880 | 21.680.062 | 1,022 | 136 | 0,2% |
| 0,2 | 492 | 1.944 | 0,902 | 21.439.265 | 1,011 | 365 | 1,9% |
| 0,3 | 492 | 1.893 | 0,926 | 21.415.877 | 1,010 | 491 | 2,1% |
| 0,4 | 492 | 1.829 | 0,958 | 21.255.547 | 1,002 | 645 | 4,5% |
| 0,5 (corrente) | 492 | 1.753 | 1,000 | 21.202.825 | 1,000 | 759 | 4,9% |
| 0,6 | 492 | 1.632 | 0,931 | 21.122.378 | 0,996 | 898 | 5,4% |
| 0,7 | 492 | 1.540 | 0,878 | 21.084.498 | 0,994 | 1.777 | 5,8% |

Che cosa si osserva.

- I campioni conservati non cambiano: la profondità si giudica sulle letture prima
  della decontaminazione.
- L'insieme delle varianti cambia con gradualità: fra valori adiacenti entra o esce dal
  2% al 7% delle varianti, e le letture conservate variano di meno del 3% sull'intero
  intervallo.
- Fra 0,6 e 0,7 i contaminanti raddoppiano (da 898 a 1.777), ma le varianti finali
  scendono solo da 1.632 a 1.540 e le letture dello 0,2%: le varianti classificate in
  più sono quasi tutte varianti rare, che il filtro di prevalenza toglierebbe comunque.

### 2.3 `prev.min_fraction`

| Valore | Campioni conservati | Varianti finali | Jaccard con la corrente | Letture conservate | Rispetto alla corrente |
|---|---|---|---|---|---|
| 0,005 | 492 | 2.179 | 0,804 | 21.230.741 | 1,001 |
| 0,010 (corrente) | 492 | 1.753 | 1,000 | 21.202.825 | 1,000 |
| 0,015 | 492 | 1.439 | 0,821 | 21.171.310 | 0,999 |
| 0,020 | 492 | 1.291 | 0,737 | 21.149.651 | 0,998 |
| 0,025 | 492 | 1.151 | 0,657 | 21.105.891 | 0,995 |
| 0,030 | 492 | 1.070 | 0,610 | 21.082.665 | 0,994 |
| 0,050 (fuori dalla griglia di decisione) | 492 | 792 | 0,452 | 20.972.942 | 0,989 |

Che cosa si osserva.

- È il parametro che più cambia l'insieme delle varianti: ogni mezzo punto percentuale
  toglie o aggiunge dal 7% al 20% delle varianti. Con il 5% abituale in letteratura le
  varianti sarebbero 792, meno della metà.
- Le letture quasi non cambiano: fra 0,005 e 0,05 l'1,2%. Le varianti che il filtro
  decide sono rare per definizione; ciò che cambia è la coda dell'insieme, non la massa
  dei conteggi.
- I campioni conservati non cambiano in nessuna configurazione: nessun campione resta
  senza letture.

### 2.4 Modalità di decontaminazione

| Modalità | Esito | Campioni conservati | Varianti finali | Jaccard con la corrente | Letture conservate | Contaminanti | Letture dei biologici rimosse come contaminanti |
|---|---|---|---|---|---|---|---|
| aggregata (corrente) | completata | 492 | 1.753 | 1,000 | 21.202.825 | 759 | 4,9% |
| per piastra, `minimum` | arresto: letture rimosse oltre `qc.max_frac_contaminant` (`E-S12-02`) | | | | | | |
| per piastra, `minimum`, con il vincolo sospeso (diagnostica) | completata | 489 | 1.134 | 0,632 | 12.757.037 | 1.138 | 44,2% |
| per piastra, `fisher` | completata | 492 | 1.895 | 0,916 | 21.587.534 | 167 | 2,1% |

Che cosa si osserva.

- Per piastra con `minimum` la catena si ferma: le letture dei campioni biologici
  classificate contaminanti sarebbero il 44,2%, oltre il limite dichiarato di 0,40. La
  riga di diagnostica, ottenuta sospendendo il limite, mostra che cosa accadrebbe: il
  40% delle letture finali in meno. Non è una configurazione candidata.
- Per piastra con `fisher` la catena si conclude e classifica 167 contaminanti contro
  759: è la più prudente nel togliere, e il suo insieme di varianti coincide per il 92%
  con quello corrente.

### 2.5 Cambiamenti di passo

I cambiamenti fra valori adiacenti, come definiti dalla regola: varianti (1 meno la
similarità di Jaccard), campioni (campioni che cambiano stato sui 770 biologici),
letture (differenza sulle letture della configurazione corrente).

| Parametro | Da | A | Varianti | Campioni | Letture |
|---|---|---|---|---|---|
| `katharoseq.target_sensitivity` | 0,70 | 0,75 | 0,050 | 0,034 | 0,004 |
| `katharoseq.target_sensitivity` | 0,75 | 0,80 | 0,059 | 0,118 | 0,025 |
| `katharoseq.target_sensitivity` | 0,80 | 0,85 | 0,009 | 0,026 | 0,010 |
| `katharoseq.target_sensitivity` | 0,85 | 0,90 | 0,063 | 0,088 | 0,030 |
| `katharoseq.target_sensitivity` | 0,90 | 0,95 | 0,025 | 0,054 | 0,051 |
| `decontam.threshold` | 0,1 | 0,2 | 0,024 | 0,000 | 0,011 |
| `decontam.threshold` | 0,2 | 0,3 | 0,026 | 0,000 | 0,001 |
| `decontam.threshold` | 0,3 | 0,4 | 0,034 | 0,000 | 0,008 |
| `decontam.threshold` | 0,4 | 0,5 | 0,042 | 0,000 | 0,003 |
| `decontam.threshold` | 0,5 | 0,6 | 0,069 | 0,000 | 0,004 |
| `decontam.threshold` | 0,6 | 0,7 | 0,056 | 0,000 | 0,002 |
| `prev.min_fraction` | 0,005 | 0,010 | 0,196 | 0,000 | 0,001 |
| `prev.min_fraction` | 0,010 | 0,015 | 0,179 | 0,000 | 0,002 |
| `prev.min_fraction` | 0,015 | 0,020 | 0,103 | 0,000 | 0,001 |
| `prev.min_fraction` | 0,020 | 0,025 | 0,108 | 0,000 | 0,002 |
| `prev.min_fraction` | 0,025 | 0,030 | 0,070 | 0,000 | 0,001 |

## 3. Applicazione della regola

La regola (decision_log.md, 1.3): il valore corrente è instabile se, per almeno una
misura, il cambiamento verso un vicino supera insieme la soglia assoluta (0,10 per le
varianti, 0,05 per i campioni, 0,05 per le letture) e tre volte la mediana dei
cambiamenti di passo di quella misura sulla griglia. Si sostituisce solo un valore
instabile, e solo con un'alternativa in zona stabile.

| Parametro | Misura | Verso il vicino inferiore | Verso il vicino superiore | Soglia assoluta | Tre volte la mediana | Supera entrambe |
|---|---|---|---|---|---|---|
| `katharoseq.target_sensitivity` | varianti | 0,063 | 0,025 | 0,10 | 0,150 | no |
| | campioni | 0,088 | 0,054 | 0,05 | 0,164 | no (supera la sola soglia assoluta, da entrambi i lati) |
| | letture | 0,030 | 0,051 | 0,05 | 0,075 | no (supera la sola soglia assoluta, verso 0,95) |
| `decontam.threshold` | varianti | 0,042 | 0,069 | 0,10 | 0,113 | no |
| | campioni | 0,000 | 0,000 | 0,05 | 0,000 | no |
| | letture | 0,003 | 0,004 | 0,05 | 0,009 | no |
| `prev.min_fraction` | varianti | 0,196 | 0,179 | 0,10 | 0,325 | no (supera la sola soglia assoluta, da entrambi i lati) |
| | campioni | 0,000 | 0,000 | 0,05 | 0,000 | no |
| | letture | 0,001 | 0,002 | 0,05 | 0,004 | no |

**Decisioni.**

| Parametro | Valore corrente | Instabile secondo la regola | Alternative in zona stabile | Decisione |
|---|---|---|---|---|
| `katharoseq.target_sensitivity` | 0,90 | no | nessuna | mantenuto: 0,90 |
| `decontam.threshold` | 0,5 | no | 0,2; 0,3; 0,4; 0,6 | mantenuto: 0,5 |
| `prev.min_fraction` | 0,01 | no | nessuna | mantenuto: 0,01 |
| modalità di decontaminazione | aggregata | ammissibile | | mantenuta: aggregata |

Nessun valore cambia: la configurazione congelata coincide, per questi quattro
parametri, con quella di partenza, e il risultato finale non cambia.

**Che cosa la regola non dice, e va letto insieme alla decisione.**

- `katharoseq.target_sensitivity`: il valore è mantenuto, ma il risultato non è
  insensibile a esso. Passando a un valore adiacente cambia stato dal 5% al 9% dei
  campioni biologici, oltre la soglia assoluta; la regola non lo giudica instabile
  perché cambiamenti di questa entità avvengono a ogni passo della griglia (mediana
  0,055), e nessun altro valore della griglia è in una zona stabile. Non esiste, in
  altre parole, un valore attorno al quale l'insieme dei campioni non dipenda dalla
  scelta: è un limite del dato (poche osservazioni per curva, otto livelli di
  diluizione per piastra) e va dichiarato a chi usa il risultato. Il numero dei
  campioni analizzati (492 su 770) dipende da questa scelta più che da ogni altra.
- `katharoseq.target_sensitivity`, valori bassi: il confronto fra 0,90 e 0,80 non è un
  confronto fra una soglia severa e una permissiva sulla stessa curva. A 0,80 le
  piastre 3 e 4 perdono la curva e ricevono la soglia fissa, la più bassa possibile:
  metà degli 88 campioni in più viene da lì. Scegliere 0,80 per conservare più campioni
  significherebbe, per quelle piastre, rinunciare al criterio dei controlli positivi.
- `prev.min_fraction`: è mantenuto perché il cambiamento, pur ampio (dal 18% al 20%
  delle varianti verso i vicini), è dello stesso ordine a ogni passo: nessun valore
  della griglia è più stabile di un altro. L'insieme delle varianti finali dipende
  quindi da questa soglia in modo sostanziale e continuo; le letture no (meno dello
  0,2% per passo). Le analisi a valle che pesano le varianti per abbondanza ne
  risentono poco; quelle che contano le varianti (ricchezza, presenza o assenza) ne
  risentono molto, e vanno interpretate sapendo che 1.753 è il numero di varianti
  presenti in almeno 5 dei 492 campioni, non una proprietà del dataset.
- `decontam.threshold`: è il parametro attorno al cui valore il risultato è più
  stabile. Anche quattro alternative sono in zona stabile; la regola non le sceglie
  perché il valore corrente non è instabile.
- Modalità di decontaminazione: la regola mantiene l'aggregata perché è ammissibile, e
  non confronta fra loro le modalità ammissibili. La modalità per piastra con `fisher`
  resta una scelta di metodo difendibile, che toglie meno (167 contaminanti contro
  759); la differenza sull'insieme finale è di 142 varianti in più e dell'1,8% delle
  letture.

## 4. Il modello evolutivo della filogenesi

La fase di filogenesi è disattivata nella configurazione di OSD-734; il confronto
decide il valore predefinito di `phylo.model`. I due modelli sono stati adattati alle
1.753 varianti finali dell'esecuzione di partenza con `scripts/confronto_modelli_filogenesi.R`:
stesso allineamento (145 colonne, 20 costanti), stesso albero di partenza, stessa
ricerca deterministica della fase.

| Modello | Log-verosimiglianza | Parametri | AIC | BIC | Forma della gamma | Quota di siti invarianti |
|---|---|---|---|---|---|---|
| GTR+G+I | -41.734,9 | 3.513 | 90.495,8 | 100.953,0 | 0,619 | 0,2 |
| GTR+G | -41.872,5 | 3.512 | 90.768,9 | 101.223,2 | 0,481 | 0 |

La differenza di BIC è 270,2 a favore di GTR+G+I (di AIC 273,2): ben oltre la soglia
di 2 fissata dalla regola. **Il valore predefinito resta GTR+G+I.**

Una riserva, da leggere con il risultato. Nel modello GTR+G+I la quota di siti
invarianti resta al valore da cui parte l'ottimizzazione (0,2): la ricerca non trova
un valore vicino che migliori la verosimiglianza, e partendo da valori diversi (0 e
0,1) arriva ad alberi con verosimiglianza più bassa. Il confronto dice quindi che il
modello con una quota di siti invarianti di 0,2 descrive i dati meglio di quello senza,
non che 0,2 sia la stima di massima verosimiglianza di quella quota. Con 145 colonne e
oltre 3.500 parametri (quasi tutti lunghezze dei rami) entrambi i modelli sono
fortemente parametrizzati rispetto ai dati: è la ragione per cui la fase resta
disattivata, e l'albero, quando richiesto, va considerato debolmente risolto.

Per rieseguire il confronto, nell'immagine della pipeline:

```bash
Rscript --vanilla scripts/confronto_modelli_filogenesi.R \
  output/osd734/12_final/intermedi/ps_filtrato.rds confronto_modelli_filogenesi.json 12
```
