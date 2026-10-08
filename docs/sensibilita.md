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
`dati/osd734/config_osd734.yaml` e la soglia di profondità omogenea
([decision_log.md](decision_log.md), sezione 4): oggetto finale di 451 campioni per
1.737 varianti, 20.922.721 letture.

**Perché l'analisi è stata ripetuta.** La prima analisi è stata fatta quando una
piastra senza curva valida riceveva una soglia fissa di 1.000 letture grezze. Con la
soglia omogenea quella piastra riceve la mediana delle soglie delle altre, sulle
letture senza chimere: cambiano i campioni conservati (da 492 a 451) e con loro il
denominatore della prevalenza. La regola di decisione è la stessa, non modificata; le
griglie sono le stesse (`sensibilita/griglie.yaml`).

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
suoi file hanno gli stessi inode, le stesse dimensioni e gli stessi istanti di
modifica di prima. Tutte le riprese sono avvenute nell'immagine pubblicata della
pipeline, da un clone senza modifiche, con la regola rigorosa sulla provenienza
attiva. Lo strumento non contiene valori del dataset: legge le griglie dal file e il
valore corrente di ogni parametro dalla configurazione.

**Ammissibilità.** Un valore non è ammissibile solo se produce una delle condizioni
previste dalla regola (letture rimosse come contaminanti oltre
`qc.max_frac_contaminant`, campioni svuotati dal filtro di prevalenza, frazione di
letture trattenute sotto `qc.min_frac_reads_retained`). Ogni altro arresto ferma lo
strumento: in questa analisi non ce ne sono stati.

**Misure.** Per ogni configurazione, sull'oggetto finale: campioni biologici
conservati e quanti cambiano stato rispetto alla configurazione corrente; varianti
finali e similarità di Jaccard del loro insieme con quello corrente; letture
conservate; soglie di profondità per piastra e piastre che non usano una curva
propria, con l'origine della soglia; contaminanti rimossi e frazione delle letture dei
campioni biologici che essi raccolgono.

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
    --griglie docs/sensibilita/griglie.yaml \
    --lavoro output/sensibilita \
    --uscita docs/sensibilita
```

Le varianti già calcolate si riusano; `--rifai` le ricalcola. Le riprese partono da
una copia e cambiano un parametro: con la regola rigorosa sulla provenienza attiva
nella configurazione, lo strumento va eseguito da un clone senza modifiche non
committate, come ogni altra esecuzione. Le venti riprese sono durate in tutto 17
minuti con dodici processori; le copie occupano 2,3 GB.

## 2. Risultati

Nelle tabelle, "campioni che cambiano" conta i campioni biologici conservati in una
configurazione ed esclusi nell'altra, rispetto alla configurazione corrente; "rispetto
alla corrente" è il rapporto fra le letture conservate.

### 2.1 `katharoseq.target_sensitivity`

| Valore | Campioni conservati | Campioni che cambiano | Varianti finali | Jaccard con la corrente | Letture conservate | Rispetto alla corrente | Piastre senza curva propria |
|---|---|---|---|---|---|---|---|
| 0,70 | 525 | 88 | 1.657 | 0,945 | 21.984.405 | 1,051 | 1, 2, 3, 4, 6, 7, 8, 9 |
| 0,75 | 545 | 98 | 1.673 | 0,945 | 22.141.599 | 1,058 | 1, 2, 3, 4, 6, 7, 8 |
| 0,80 | 489 | 66 | 1.769 | 0,975 | 21.582.124 | 1,032 | 1, 2, 3, 4 |
| 0,85 | 459 | 42 | 1.743 | 0,986 | 21.253.582 | 1,016 | 1, 2, 3, 4 |
| 0,90 (corrente) | 451 | 0 | 1.737 | 1,000 | 20.922.721 | 1,000 | 1, 2 |
| 0,95 | 394 | 57 | 1.829 | 0,937 | 19.611.097 | 0,937 | 1, 2 |

Soglie di profondità per piastra, tutte sulle letture senza chimere; (m) indica la
mediana delle soglie delle piastre con curva propria valida:

| Valore | P1 | P2 | P3 | P4 | P5 | P6 | P7 | P8 | P9 | P10 |
|---|---|---|---|---|---|---|---|---|---|---|
| 0,70 | 10.173 (m) | 10.173 (m) | 10.173 (m) | 10.173 (m) | 6.284 | 10.173 (m) | 10.173 (m) | 10.173 (m) | 10.173 (m) | 14.061 |
| 0,75 | 8.531 (m) | 8.531 (m) | 8.531 (m) | 8.531 (m) | 6.663 | 8.531 (m) | 8.531 (m) | 8.531 (m) | 8.531 | 14.962 |
| 0,80 | 12.223 (m) | 12.223 (m) | 12.223 (m) | 12.223 (m) | 7.127 | 11.562 | 12.884 | 27.878 | 9.323 | 16.074 |
| 0,85 | 14.899 (m) | 14.899 (m) | 14.899 (m) | 14.899 (m) | 7.738 | 13.969 | 15.828 | 33.168 | 10.394 | 17.542 |
| 0,90 | 15.049 (m) | 15.049 (m) | 9.984 | 6.376 | 8.642 | 18.063 | 20.948 | 41.968 | 12.034 | 19.727 |
| 0,95 | 19.618 (m) | 19.618 (m) | 12.480 | 8.105 | 10.360 | 27.763 | 33.513 | 62.071 | 15.320 | 23.916 |

Che cosa si osserva.

- Il parametro agisce in due modi distinti. Dove la curva della piastra è valida, la
  soglia cresce con continuità con la sensibilità richiesta. Ma abbassando la
  sensibilità la soglia scende verso il punto medio della curva, e quando fra il punto
  medio e la soglia non resta alcun controllo positivo osservato la curva non è più
  ritenuta valida: la piastra riceve allora la mediana delle soglie delle piastre
  rimaste valide. A 0,90 non hanno una curva propria le piastre 1 e 2; a 0,85 e a 0,80
  anche la 3 e la 4; a 0,75 sette piastre; a 0,70 otto su dieci.
- Con la soglia omogenea il venir meno di una curva non apre più la piastra a tutti i
  suoi campioni: la piastra riceve una soglia dello stesso ordine delle altre. A 0,85
  le piastre 3 e 4 passano dalla propria soglia (9.984 e 6.376 letture) alla mediana
  (14.899), più alta, e perdono 17 campioni invece di guadagnarne; i 25 campioni in
  più vengono dalle piastre 6-10, la cui soglia scende lungo la curva. Con il ripiego
  precedente, a 0,85 le piastre 3 e 4 guadagnavano 43 campioni per il solo passaggio
  alla soglia fissa.
- Ai valori bassi la mediana si calcola su poche piastre (tre a 0,75, due a 0,70) e
  diventa essa stessa instabile: scende da 12.223 a 8.531 fra 0,80 e 0,75 e risale a
  10.173 a 0,70. I campioni conservati non sono quindi monotoni nel parametro (489,
  545, 525). Sotto 0,80 la soglia della maggior parte delle piastre non viene più dai
  loro controlli.
- Il numero di varianti finali non cresce con il numero di campioni: il filtro di
  prevalenza chiede una frazione dei campioni conservati, e i campioni aggiunti, poco
  profondi, portano poche varianti condivise.
- Verso 0,95 escono 57 campioni, da tutte le piastre, e il 6,3% delle letture: è il
  cambiamento più ampio attorno al valore corrente.
- I contaminanti non dipendono da questo parametro: 759 in ogni configurazione.

### 2.2 `decontam.threshold`

| Valore | Campioni conservati | Varianti finali | Jaccard con la corrente | Letture conservate | Rispetto alla corrente | Contaminanti | Letture dei biologici rimosse come contaminanti |
|---|---|---|---|---|---|---|---|
| 0,1 | 451 | 1.966 | 0,883 | 21.367.488 | 1,021 | 136 | 0,2% |
| 0,2 | 451 | 1.920 | 0,905 | 21.146.069 | 1,011 | 365 | 1,9% |
| 0,3 | 451 | 1.872 | 0,928 | 21.123.748 | 1,010 | 491 | 2,1% |
| 0,4 | 451 | 1.809 | 0,960 | 20.974.592 | 1,002 | 645 | 4,5% |
| 0,5 (corrente) | 451 | 1.737 | 1,000 | 20.922.721 | 1,000 | 759 | 4,9% |
| 0,6 | 451 | 1.617 | 0,931 | 20.844.186 | 0,996 | 898 | 5,4% |
| 0,7 | 451 | 1.525 | 0,878 | 20.807.209 | 0,995 | 1.777 | 5,8% |

Che cosa si osserva.

- I campioni conservati non cambiano: la profondità si giudica sulle letture prima
  della decontaminazione.
- L'insieme delle varianti cambia con gradualità: fra valori adiacenti entra o esce dal
  2% al 7% delle varianti, e le letture conservate variano di meno del 3% sull'intero
  intervallo.
- Fra 0,6 e 0,7 i contaminanti raddoppiano (da 898 a 1.777), ma le varianti finali
  scendono solo da 1.617 a 1.525 e le letture dello 0,2%: le varianti classificate in
  più sono quasi tutte varianti rare, che il filtro di prevalenza toglierebbe comunque.

### 2.3 `prev.min_fraction`

| Valore | Campioni conservati | Varianti finali | Jaccard con la corrente | Letture conservate | Rispetto alla corrente |
|---|---|---|---|---|---|
| 0,005 | 451 | 2.161 | 0,804 | 20.950.337 | 1,001 |
| 0,010 (corrente) | 451 | 1.737 | 1,000 | 20.922.721 | 1,000 |
| 0,015 | 451 | 1.504 | 0,866 | 20.901.461 | 0,999 |
| 0,020 | 451 | 1.286 | 0,740 | 20.871.507 | 0,998 |
| 0,025 | 451 | 1.178 | 0,678 | 20.833.210 | 0,996 |
| 0,030 | 451 | 1.095 | 0,630 | 20.812.245 | 0,995 |
| 0,050 (fuori dalla griglia di decisione) | 451 | 814 | 0,469 | 20.710.029 | 0,990 |

Che cosa si osserva.

- È il parametro che più cambia l'insieme delle varianti: ogni mezzo punto percentuale
  toglie o aggiunge dal 7% al 20% delle varianti. Con il 5% abituale in letteratura le
  varianti sarebbero 814, meno della metà.
- Con 451 campioni conservati il minimo è 5 campioni all'1%, come con i 492 di prima
  (l'intero superiore di 4,51 e di 4,92): il minimo non è cambiato, è cambiato
  l'insieme dei campioni su cui si conta. A 0,005 è 3 campioni, a 0,015 è 7.
- Le letture quasi non cambiano: fra 0,005 e 0,05 l'1,2%. Le varianti che il filtro
  decide sono rare per definizione; ciò che cambia è la coda dell'insieme, non la massa
  dei conteggi.
- I campioni conservati non cambiano in nessuna configurazione: nessun campione resta
  senza letture.

### 2.4 Modalità di decontaminazione

| Modalità | Esito | Campioni conservati | Varianti finali | Jaccard con la corrente | Letture conservate | Contaminanti | Letture dei biologici rimosse come contaminanti |
|---|---|---|---|---|---|---|---|
| aggregata (corrente) | completata | 451 | 1.737 | 1,000 | 20.922.721 | 759 | 4,9% |
| per piastra, `minimum` | arresto: letture rimosse oltre `qc.max_frac_contaminant` (`E-S12-02`) | | | | | | |
| per piastra, `minimum`, con il vincolo sospeso (diagnostica) | completata | 451 | 1.115 | 0,631 | 12.593.728 | 1.138 | 44,2% |
| per piastra, `fisher` | completata | 451 | 1.870 | 0,920 | 21.276.987 | 167 | 2,1% |

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
| `katharoseq.target_sensitivity` | 0,70 | 0,75 | 0,010 | 0,031 | 0,007 |
| `katharoseq.target_sensitivity` | 0,75 | 0,80 | 0,063 | 0,073 | 0,027 |
| `katharoseq.target_sensitivity` | 0,80 | 0,85 | 0,015 | 0,039 | 0,016 |
| `katharoseq.target_sensitivity` | 0,85 | 0,90 | 0,014 | 0,054 | 0,016 |
| `katharoseq.target_sensitivity` | 0,90 | 0,95 | 0,063 | 0,074 | 0,063 |
| `decontam.threshold` | 0,1 | 0,2 | 0,023 | 0,000 | 0,011 |
| `decontam.threshold` | 0,2 | 0,3 | 0,025 | 0,000 | 0,001 |
| `decontam.threshold` | 0,3 | 0,4 | 0,034 | 0,000 | 0,007 |
| `decontam.threshold` | 0,4 | 0,5 | 0,040 | 0,000 | 0,003 |
| `decontam.threshold` | 0,5 | 0,6 | 0,069 | 0,000 | 0,004 |
| `decontam.threshold` | 0,6 | 0,7 | 0,057 | 0,000 | 0,002 |
| `prev.min_fraction` | 0,005 | 0,010 | 0,196 | 0,000 | 0,001 |
| `prev.min_fraction` | 0,010 | 0,015 | 0,134 | 0,000 | 0,001 |
| `prev.min_fraction` | 0,015 | 0,020 | 0,145 | 0,000 | 0,001 |
| `prev.min_fraction` | 0,020 | 0,025 | 0,084 | 0,000 | 0,002 |
| `prev.min_fraction` | 0,025 | 0,030 | 0,070 | 0,000 | 0,001 |

## 3. Applicazione della regola

La regola (decision_log.md, 1.3): il valore corrente è instabile se, per almeno una
misura, il cambiamento verso un vicino supera insieme la soglia assoluta (0,10 per le
varianti, 0,05 per i campioni, 0,05 per le letture) e tre volte la mediana dei
cambiamenti di passo di quella misura sulla griglia. Si sostituisce solo un valore
instabile, e solo con un'alternativa in zona stabile.

| Parametro | Misura | Verso il vicino inferiore | Verso il vicino superiore | Soglia assoluta | Tre volte la mediana | Supera entrambe |
|---|---|---|---|---|---|---|
| `katharoseq.target_sensitivity` | varianti | 0,014 | 0,063 | 0,10 | 0,044 | no (verso 0,95 supera la sola sproporzione) |
| | campioni | 0,054 | 0,074 | 0,05 | 0,164 | no (supera la sola soglia assoluta, da entrambi i lati) |
| | letture | 0,016 | 0,063 | 0,05 | 0,047 | **sì, verso 0,95** |
| `decontam.threshold` | varianti | 0,040 | 0,069 | 0,10 | 0,110 | no |
| | campioni | 0,000 | 0,000 | 0,05 | 0,000 | no |
| | letture | 0,003 | 0,004 | 0,05 | 0,009 | no |
| `prev.min_fraction` | varianti | 0,196 | 0,134 | 0,10 | 0,402 | no (supera la sola soglia assoluta, da entrambi i lati) |
| | campioni | 0,000 | 0,000 | 0,05 | 0,000 | no |
| | letture | 0,001 | 0,001 | 0,05 | 0,004 | no |

**Decisioni.**

| Parametro | Valore corrente | Instabile secondo la regola | Alternative in zona stabile | Decisione |
|---|---|---|---|---|
| `katharoseq.target_sensitivity` | 0,90 | sì (letture, verso 0,95) | nessuna | mantenuto: 0,90, con l'instabilità dichiarata come limite |
| `decontam.threshold` | 0,5 | no | 0,2; 0,3; 0,4; 0,6 | mantenuto: 0,5 |
| `prev.min_fraction` | 0,01 | no | 0,025 | mantenuto: 0,01 |
| modalità di decontaminazione | aggregata | ammissibile | | mantenuta: aggregata |

Nessun valore cambia: la regola indica per ogni parametro il valore già congelato, e
la configurazione congelata resta quella.

**Che cosa è cambiato rispetto alla prima analisi.** Allora
`katharoseq.target_sensitivity` non risultava instabile; ora lo è, per le letture
verso 0,95 (6,3% contro una soglia di sproporzione del 4,7%). Non è il passo verso
0,95 a essere cresciuto (era il 5,1%): è la mediana dei passi a essere scesa, da 0,025
a 0,016, perché con la soglia omogenea i passi fra i valori bassi non sono più
gonfiati dalle piastre che ricevevano la soglia fissa. La regola chiede in questo caso
un'alternativa in zona stabile, cioè con due vicini ammissibili e tutti i cambiamenti
sotto le soglie assolute: non ce n'è nessuna (0,85 cambia il 5,4% dei campioni verso
0,90; 0,80 il 7,3% verso 0,75; 0,75 il 7,3% verso 0,80). Il valore si mantiene, e
l'instabilità si dichiara come limite noto.

**Che cosa la regola non dice, e va letto insieme alla decisione.**

- `katharoseq.target_sensitivity`: il valore è mantenuto, ma il risultato non è
  insensibile a esso. Verso 0,95 escono 57 campioni e il 6,3% delle letture; verso
  0,85 cambiano stato 42 campioni. Nessun valore della griglia è in una zona stabile:
  non esiste un valore attorno al quale l'insieme dei campioni non dipenda dalla
  scelta. È un limite del dato (poche osservazioni per curva, otto livelli di
  diluizione per piastra) e va dichiarato a chi usa il risultato. Il numero dei
  campioni analizzati (451 su 770) dipende da questa scelta più che da ogni altra.
- `katharoseq.target_sensitivity`, valori bassi: sotto 0,90 le piastre 3 e 4 perdono
  la curva, sotto 0,80 la perdono altre tre o quattro piastre, e la soglia di tutte
  viene dalla mediana di due o tre curve. Scegliere un valore basso per conservare più
  campioni significherebbe, per la maggior parte delle piastre, rinunciare al
  criterio dei propri controlli positivi.
- `prev.min_fraction`: è mantenuto perché il cambiamento, pur ampio (dal 13% al 20%
  delle varianti verso i vicini), è dello stesso ordine a ogni passo. Un'alternativa
  risulta in zona stabile (0,025), ma la regola non la sceglie perché il valore
  corrente non è instabile; con essa le varianti sarebbero 1.178. L'insieme delle
  varianti finali dipende da questa soglia in modo sostanziale e continuo; le letture
  no (meno dello 0,2% per passo). Le analisi a valle che pesano le varianti per
  abbondanza ne risentono poco; quelle che contano le varianti (ricchezza, presenza o
  assenza) ne risentono molto, e vanno interpretate sapendo che 1.737 è il numero di
  varianti presenti in almeno 5 dei 451 campioni, non una proprietà del dataset.
- `decontam.threshold`: è il parametro attorno al cui valore il risultato è più
  stabile. Anche quattro alternative sono in zona stabile; la regola non le sceglie
  perché il valore corrente non è instabile.
- Modalità di decontaminazione: la regola mantiene l'aggregata perché è ammissibile, e
  non confronta fra loro le modalità ammissibili. La modalità per piastra con `fisher`
  resta una scelta di metodo difendibile, che toglie meno (167 contaminanti contro
  759); la differenza sull'insieme finale è di 133 varianti in più e dell'1,7% delle
  letture.

## 4. Il modello evolutivo della filogenesi

La fase di filogenesi è disattivata nella configurazione di OSD-734; il confronto
decide il valore predefinito di `phylo.model`. Non è stato ripetuto con la soglia
omogenea: le misure che seguono sono della prima analisi, sulle 1.753 varianti finali
dell'esecuzione di riferimento di allora (oggi le varianti finali sono 1.737, di cui
tutte fra quelle 1.753). I due modelli sono stati adattati con `scripts/confronto_modelli_filogenesi.R`:
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
