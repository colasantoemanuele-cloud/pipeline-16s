# Guida: da un dataset 16S al report e alle analisi ecologiche

Per chi sa usare un terminale e ha un dataset 16S Illumina, proprio o di OSDR (l'archivio aperto
della NASA). Si seguono i passi nell'ordine: ogni passo dà il comando, che cosa si vede e che cosa
si scrive. Per un'analisi normale non serve altro: [RIFERIMENTO.md](RIFERIMENTO.md) spiega ogni
parametro e i casi rari.

Nei messaggi dello strumento «GUIDA.md, passo N» è un passo di questa guida; una sezione numerata
(«RIFERIMENTO.md, sezione 6», «GUIDA.md, 6.3») o un titolo fra virgolette è in RIFERIMENTO.md.

## Passo 1. Requisiti

- **Docker** (su un cluster senza Docker: Apptainer, RIFERIMENTO.md 3.1) e **Python 3** per gli
  script di scarico.
- **L'immagine pubblicata**, scaricata per digest (1,86 GB da scaricare, 7,6 GB su disco):

  ```bash
  docker pull ghcr.io/colasantoemanuele-cloud/amplicon16s@sha256:930dda008581b24cc3a3fea718c0fed3690b011958bf1716c643f92056970154
  ```

- **Un clone di questo repository**: ogni comando si lancia dalla sua radice.

  ```bash
  git clone https://github.com/colasantoemanuele-cloud/pipeline-16s.git
  cd pipeline-16s
  ```

- **Memoria e tempo**: per una decina di campioni pochi GB e 10-20 minuti; per 960, 13 GB e 50.

## Passo 2. Preparare i dati

**Da OSDR**, con l'identificativo dello studio (`OSD-` e il numero) e una cartella dei dati fuori dal clone:

```bash
python3 scripts/scarica_osdr.py OSD-NNN /percorso/dei/dati
```

Lo script crea `metadati/` (le tabelle ISA-Tab dello studio) e `letture/` (i file elencati nella
tabella di assay del 16S). Prima di scaricare dice quanti file e quanti GB; alla fine stampa il
comando di `init` del passo 4, da copiare. Si può rilanciare dopo un'interruzione: riprende. Se
lo studio ha più tabelle del 16S le nomina e si ferma: si rilancia con `--assay NOME`.

**Coppie o sole forward.** Se lo studio ha due file per campione (`_R1` e `_R2`) si parte con le
coppie. Si passa alle sole forward (i file `_R1`) in due casi: lo script avverte che nei nomi il
marcatore non è in fondo (`_R1_raw`: lo strumento non li riconosce come coppie); oppure `init`, al
passo 4, scrive la nota «Le coppie non si fondono senza basi di qualità bassa» o lascia da
scegliere i troncamenti consigliando le forward (le reverse sono troppo corte o di qualità bassa).
Si preparano così, e si rilancia il comando del passo 4 con `--reads /dati/letture_forward`:

```bash
mkdir /percorso/dei/dati/letture_forward
cp /percorso/dei/dati/letture/*_R1*.fastq.gz /percorso/dei/dati/letture_forward/
```

**Con dati propri** la cartella dei dati, fuori dal clone, contiene:

- le letture, in una cartella senza sottocartelle: file `.fastq.gz` (o `.fq.gz`, `.fastq`, `.fq`),
  uno per campione, oppure due con lo stesso nome salvo `_R1`/`_R2` (o `_1`/`_2`) in fondo. Il
  nome di ogni file contiene una parte unica per campione, che compare anche nei metadati. Non
  vanno messi file di campioni assenti dai metadati né di un altro amplicone (ITS, 18S);
- i metadati: una tabella separata da tabulazioni, una riga per campione (colonne con il nome del
  file, il nome del campione, il tipo se ci sono controlli, le variabili), o le tabelle ISA-Tab.

**Il riferimento tassonomico**, una volta sola, dal clone:

```bash
python3 dati/riferimento/scarica_riferimento.py     # SILVA 138, 138 MB, in dati/riferimento/
```

## Passo 3. Un modo solo di lanciare i comandi

Nel terminale, dalla radice del clone, si definisce una volta il comando `p16s`:

```bash
IMMAGINE=ghcr.io/colasantoemanuele-cloud/amplicon16s@sha256:930dda008581b24cc3a3fea718c0fed3690b011958bf1716c643f92056970154
DATI=/percorso/dei/dati        # la cartella dei dati del passo 2
LAVORO=/percorso/di/lavoro     # una cartella nuova, per configurazione e uscite
mkdir -p "$LAVORO"
p16s() { docker run --rm -u "$(id -u):$(id -g)" -v "$(pwd)":/app -v "$DATI":/dati:ro \
           -v "$LAVORO":/lavoro -w /app "$IMMAGINE" python3 scripts/esegui.py "$@"; }
```

Da qui ogni passo è una riga che comincia con `p16s`, lanciata dalla radice del clone. Dentro il
container la cartella dei dati si chiama `/dati`, quella di lavoro `/lavoro`, il clone `/app`: nei
comandi e nella configurazione si scrivono questi nomi, non i percorsi del proprio disco. In un
terminale nuovo le righe qui sopra si ripetono. La cartella di lavoro va creata prima (`mkdir`):
se la crea Docker appartiene a root. Dove un messaggio dice `python3 scripts/esegui.py ...`, si
scrive `p16s ...`.

## Passo 4. `init`: la configurazione di partenza

```bash
p16s init --reads /dati/letture --metadata /dati/metadati \
  --reference /app/dati/riferimento/silva_nr99_v138_train_set.fa.gz \
  --out-root /lavoro/esecuzione --output /lavoro/configurazione.yaml --accetta-proposte
```

- `--metadata` è la cartella ISA-Tab, oppure la tabella dei campioni (`/dati/campioni.tsv`). Se la
  cartella ha più tabelle di assay (16S e ITS) `init` le nomina e si aggiunge `--assay NOME`, con
  il nome di quella del 16S: lo script del passo 2 lo scrive già nel comando che stampa. Il file
  del lotto, se c'è, si indica con `--batch-table /dati/lotto.tsv`.
- `init` legge letture e metadati (da decine di secondi ad alcuni minuti) e scrive
  `/lavoro/configurazione.yaml`: ciò che ha verificato e ciò che propone sono righe attive, le
  soglie di qualità sono al valore che non ferma (si fissano al passo 8), ciò che non può
  decidere è una riga in commento. Sopra ogni riga un commento dice da dove viene il valore; una
  riga attiva a `null` o `[]` dichiara un'assenza, e resta così se lo studio non ha quel dato.
- **Il riepilogo va letto**, sulla console o in testa al file: soprattutto le righe «Nota:». Una
  nota che nomina un parametro dice che `init` ha lasciato quella scelta a chi conosce i dati.
  **Come si risolve:** si apre il file alla riga del parametro; il commento sopra la riga riporta
  ciò che `init` ha visto e la sua raccomandazione; si segue quella, salvo sapere di più sullo
  studio. Una nota che non nomina parametri («resta da verificare») informa, e non chiede nulla.
- **Il primer.** Se è proposto dai metadati e il commento dice che «non è nelle letture», è già
  stato tolto: non c'è nulla da fare. Se la nota dice che non è stato trovato e lo si conosce, si
  rilancia con `--primer SEQUENZA` (con le coppie anche `--primer-reverse SEQUENZA`); se non lo
  si conosce si prosegue, e `qc.primer_sequence` resta `null`.
- `init` non scrive su un file che esiste: per rilanciarlo si cancella il file di prima.

## Passo 5. Le scelte da scrivere

```bash
p16s validate --config /lavoro/configurazione.yaml
```

`validate` elenca i parametri che mancano (`E-G15-10`) e, sotto «dove scriverli», **la riga del
file** di ciascuno. A quella riga c'è il parametro in commento, sotto il suo gruppo: per
`filt.exclude_taxa` la riga `  # exclude_taxa:` sotto `filt:`. Si toglie `# ` lasciando i due
spazi e si scrive il valore dopo i due punti: le stringhe fra virgolette, gli elenchi fra
parentesi quadre (`["a", "b"]`), i logici `true` o `false`.

**La chiave del campione.** Se `validate` indica per `meta.accession_column` e
`io.accession_regex` più righe «pronte», più colonne dei metadati corrispondono ai file. Nel file
ogni candidata è una coppia di righe già scritte, numerata («coppia 1», «coppia 2»): si toglie il
commento alle due righe della stessa coppia, una sotto `meta:` e una sotto `io:`, scegliendo
quella marcata «consigliata». Non c'è nessuna espressione da scrivere.

| Parametro | Che cosa decide | Regola semplice |
|---|---|---|
| `tax.ref_name`, `tax.ref_version` | nome e versione del riferimento, riportati con i risultati | `"SILVA"` e `"138"` per il file del passo 2 (di solito `init` li scrive da sé) |
| `filt.remove_na_phylum` | se togliere le varianti a cui il riferimento non dà un phylum | `true`: ciò che il riferimento non riconosce di solito non è un batterio |
| `filt.exclude_taxa` | i taxa da togliere | `["Chloroplast", "Mitochondria", "Eukaryota"]`: non sono batteri del campione |
| `prev.apply` | se togliere le varianti presenti in pochi campioni | `false` con poche decine di campioni, o se interessa ogni singolo campione; con `true` vanno scritti anche `prev.min_fraction` (con 20 campioni e `0.2` una variante deve comparire in 4) e `prev.min_count` (`1`) |
| `qc.min_reads_final` | le letture sotto le quali un campione esce dal risultato | `0` alla prima esecuzione: il valore giusto lo propone il passo 8 |
| `out.study_columns` | le colonne dei metadati da portare nel risultato: senza, le analisi ecologiche non hanno nulla da confrontare. `validate` non la chiede: va controllata comunque | le variabili dello studio, con il nome esatto dell'intestazione. In ISA-Tab sono le colonne «Factor Value[...]», e `init` le propone già; altrimenti la riga è attiva con `[]`, e si scrive l'elenco al posto di `[]`: `["Trattamento", "Tempo"]`. Il commento elenca le colonne disponibili |
| `ctrl.column` e le tre etichette `ctrl.*_values` | quali campioni sono controlli negativi (bianchi) o positivi | se `init` le ha proposte, controllare nel commento quali campioni porta ogni valore; `null` e tre elenchi vuoti se non ci sono controlli. Colture pure e standard senza diluizioni si lasciano fra i biologici |
| `decontam.threshold`, `decontam.min_blanks` | solo con controlli negativi: quanto è severa la decontaminazione, e quanti bianchi servono | `0.1`, il valore del metodo (più alto toglie di più, e si usa con bassa biomassa); il numero dei bianchi, non più di quelli che ci sono |
| `katharoseq.target_taxon`, `katharoseq.target_sensitivity`, `katharoseq.min_r2`, `ctrl.min_positive_pass_frac` | solo con controlli positivi in diluizione: il genere atteso e come si calibra la profondità minima | il genere del ceppo usato; gli altri tre sono scelte dello studio, fra 0 e 1: RIFERIMENTO.md 6.5 |
| `filter.truncLen` (con le coppie anche `pair.truncLen_reverse`, `pair.len_min`, `pair.len_max`) | dove si troncano le letture, e quanto sono lunghe le coppie fuse | di solito proposti. Se mancano, il commento riporta la raccomandazione; se consiglia le sole forward, vedi il passo 2 |

**Il classificatore.** Non c'è nulla da scrivere per quello predefinito, il bayesiano ingenuo
(`tax.classifier: naive_bayes`): è la scelta consigliata per un'analisi al livello di genere,
perché assegna più generi. L'alternativa è IdTaxa (`tax.classifier: idtaxa`), più prudente: assegna
meno generi ma sbaglia meno quando il batterio non è nel riferimento. È il classificatore della
pipeline di NASA GeneLab, e serve per confrontarsi con le sue tabelle: il file `2019` (SILVA r138)
corrisponde alle versioni precedenti di quella pipeline, il `138.2` a quella attuale. Per usarlo:

```bash
python3 dati/riferimento/scarica_addestrato.py 2019     # oppure 138.2; stampa percorso e MD5
```

Nel file, sotto `tax:`, si mettono in commento (`# ` davanti) le righe `ref_fasta`, `ref_md5` e
`try_rc`, e si aggiungono queste, con l'MD5 stampato dallo script (per il `138.2` il file è
`SILVA_SSU_r138_2_2024.RData`, e `ref_version` diventa `"138.2"`):

```yaml
  classifier: idtaxa
  idtaxa_trained: "/app/dati/riferimento/SILVA_SSU_r138_2019.RData"
  idtaxa_trained_md5: "<MD5 stampato dallo script>"
  idtaxa_threshold: 60
  idtaxa_strand: both
```

Con il file `138.2` i taxa da togliere si scrivono `["Chloroplast_2", "Mitochondria_2",
"Eukaryota"]`: quel classificatore li nomina con il suffisso.

## Passo 6. `validate`, finché passa

Si rilancia il comando del passo 5 dopo ogni correzione, finché scrive «Validazione superata». Gli
avvisi (`WARNING`) non fermano. Sono attesi, qui e durante `run`: i controlli di qualità spenti
(passo 8); `E-S0-17` ed `E-S12-03` senza controlli negativi o positivi; `E-S0-31`, le righe di un
altro amplicone, escluse. `E-S0-32` va corretto: quel taxon di `filt.exclude_taxa` non toglie nulla.

## Passo 7. `run` e `status`

```bash
rm -r "$LAVORO/esecuzione"      # validate vi ha scritto, e run parte solo da una cartella vuota
p16s run --config /lavoro/configurazione.yaml
p16s status --out /lavoro/esecuzione
```

`status` si può lanciare da un altro terminale, dopo le righe del passo 3. Se `run` si ferma dice
la fase, il codice e che cosa fare (passo 11); corretta la causa si riparte con `p16s resume
--config /lavoro/configurazione.yaml`. Il risultato è in `/lavoro/esecuzione/12_final/`.

## Passo 8. `soglie` e `resume`

```bash
p16s soglie --config /lavoro/configurazione.yaml --output /lavoro/con_soglie.yaml --accetta-proposte
```

`soglie` legge le misure dell'esecuzione e scrive un file nuovo, uguale al primo salvo le soglie
di qualità: proteggono dai peggioramenti, non dicono se la misura è buona. Due cose da guardare:

- **Quante letture restano** (`qc.min_frac_reads_retained`, «letture finali su letture senza
  chimere»). Con campioni di piante o di un ospite è normale che ne resti meno di un terzo: il
  resto sono cloroplasti e mitocondri, tolti dal filtro. Che cosa è stato tolto, con le letture, è
  in `/lavoro/esecuzione/report/tabelle/varianti_rimosse.tsv`, dopo il passo 9.
- **Le letture minime per campione.** Se alcuni campioni hanno molte meno letture degli altri
  `soglie` li elenca in fondo. Sotto 1.000 letture raccomanda di toglierli, perché le analisi
  ecologiche portano ogni campione alla profondità del più povero: si scrive il valore indicato in
  `con_soglie.yaml`, nella riga `min_reads_final:`. Se hanno più di 1.000 letture di solito si
  tengono. In ogni caso, se `soglie` scrive «ATTENZIONE» un gruppo dello studio resterebbe senza
  campioni: allora si tengono, salvo rinunciare a quel confronto.

```bash
p16s resume --config /lavoro/con_soglie.yaml
```

`resume` attiva le soglie. **Ripete tutti i calcoli: dura quanto `run`**; se `qc.min_reads_final`
non è cambiato i risultati sono gli stessi. Da qui in poi la configurazione è `con_soglie.yaml`;
un `resume` successivo, dopo aver cambiato solo `qc.min_reads_final`, dura meno di un minuto.

## Passo 9. `report`

```bash
p16s report --config /lavoro/con_soglie.yaml
```

Scrive `$LAVORO/esecuzione/report/report.html`, che si apre con un browser, e le stesse tabelle in
testo in `report/tabelle/`. Da guardare per primi:

1. **Segnalazioni in apertura**: i controlli di qualità spenti e le «degradazioni», cioè gli
   avvisi con cui l'esecuzione ha proseguito rinunciando a qualcosa;
2. **Tracciamento delle letture**: quante ne restano dopo ogni fase; una perdita forte dice dove
   guardare (con le coppie, la quota di coppie fuse);
3. **Risultato finale** e **Decisioni prese automaticamente**: campioni e varianti rimasti, quali
   campioni sono usciti e perché, le varianti tolte dai filtri.

## Passo 10. Analisi ecologiche

```bash
p16s eco init --object /lavoro/esecuzione/12_final/ps_final.rds --output /lavoro/eco.yaml --accetta-proposte
p16s eco validate --object /lavoro/esecuzione/12_final/ps_final.rds --config /lavoro/eco.yaml
```

`eco init` scrive `eco.yaml` e, sotto «Rilevato sull'oggetto», elenca le colonne dei campioni:
quelle portate con `out.study_columns` hanno il nome in minuscolo con i trattini bassi (`Factor
Value[Tempo]` diventa `factor_value_tempo`). Tre righe vanno lette:

- **«Gruppi delle colonne candidate»**: per ogni colonna, i gruppi e quanti campioni hanno. La
  variabile da confrontare è una di queste: quella della domanda dello studio, che si legge in
  `metadati/i_Investigation.txt` (righe «Study Title» e «Study Description») o sulla pagina
  dello studio.
- **«Variabili confuse»**: compare se due colonne dividono i campioni quasi allo stesso modo (per
  esempio tutti i campioni a terra sono radici). Un confronto sull'una misura anche l'altra: lo
  si dichiara nei risultati, oppure si restringe l'analisi con `design.subset` ai campioni di un
  solo gruppo dell'altra colonna (le sole radici), se ne restano almeno due per gruppo.
- **«Rarefazione»**: compare se il campione più povero ha meno di un decimo delle letture della
  mediana. La riga stessa dice che cosa fare: sopra 1.000 letture di solito nulla.

Dopo ogni `resume` si rilanciano `report` ed `eco init` (con un altro `--output`).
`eco validate` nomina le scelte che mancano (`E-ECO-02`), senza il numero della riga: ogni
parametro è nel file sotto il suo gruppo (`design.variable` è `# variable:` sotto `design:`).

| Parametro | Che cosa decide | Regola semplice |
|---|---|---|
| `design.variable` | la colonna con i gruppi da confrontare | quella della domanda dello studio, fra le colonne elencate; se manca dall'oggetto si aggiunge a `out.study_columns` in `con_soglie.yaml` e si rilanciano `resume` ed `eco init` |
| `run.seed` | il seme di rarefazione e permutazioni | un intero qualunque, per esempio `1` |
| `comp.rank`, `comp.top_n` | il rango della composizione e quanti taxa distinguere nei grafici | il rango è di solito proposto (`"Genus"`); `10` |
| `beta.distances` | le distanze fra campioni | `["bray", "jaccard"]` |
| `ord.methods`, `ord.distances` | le ordinazioni, e su quali distanze (è la riga `distances:` sotto `ord:`) | `["pcoa"]` e `["bray", "jaccard"]`; con `pcoa` si toglie il commento anche a `ord.pcoa_correction`, che ha già il valore |
| `stat.min_group_size`, `stat.significance_level` | i campioni che un gruppo deve avere per entrare nei test, e il livello di significatività | `2` (`3` se ogni gruppo da confrontare ha almeno tre campioni) e `0.05` |
| `design.subset` | facoltativo: quali campioni entrano nelle analisi | una colonna e i valori da tenere: `subset: {factor_value_tempo: ["Giorno 1", "Giorno 7"]}` |

```bash
p16s eco run --object /lavoro/esecuzione/12_final/ps_final.rds --config /lavoro/eco.yaml --out /lavoro/analisi
```

`eco run` vuole una cartella di uscita assente o vuota. In `/lavoro/analisi` scrive tabelle
(`.tsv`) e grafici (`.png`) in `alfa`, `composizione`, `beta`, `ordinazione` e `test`, e
`riepilogo.json` con gli avvisi (`E-ECO-`): non fermano, dicono che cosa non si è potuto calcolare
o va letto con cautela. Variabili continue, misure ripetute, UniFrac, NMDS: RIFERIMENTO.md, 4.

## Passo 11. Se qualcosa va storto

Ogni arresto porta un codice, il dettaglio e che cosa fare: va letto per intero. I più comuni:

| Codice o sintomo | Che cosa vuol dire | Che cosa fare |
|---|---|---|
| `E-G15-10` | mancano parametri obbligatori | scriverli alle righe che `validate` indica (passo 5) |
| `E-G15-99` | un valore non è ammesso, una chiave è scritta male o ripetuta | correggere la riga che il messaggio nomina; un parametro compare una volta sola nel suo gruppo |
| `E-S0-06`, `E-S0-04` | file e metadati non corrispondono: un file senza riga, una riga senza file | togliere dalla cartella il file estraneo; se la chiave è sbagliata, scegliere l'altra coppia (passo 5) |
| `E-S0-30`, o `init` che non propone la chiave per «righe collegate a un file di un altro amplicone» | i metadati del deposito sono sbagliati | copiare le tabelle in una cartella propria, correggere la copia (ogni campione al suo file) e rilanciare `init` sulla copia |
| `E-S0-11` | un campione non ricade in nessuna etichetta | aggiungere il suo valore a una delle tre etichette `ctrl.*_values` |
| `E-S0-24`, `E-S0-23` | due file per campione con la modalità single-end; oppure a un file manca il compagno | rilanciare `init`, o lasciare nella cartella i soli file forward (passo 2) |
| `E-S0-09`, `E-S1-02` | troppe letture sono più corte di `filter.truncLen` | abbassare `filter.truncLen`: il commento della riga, nel file, riporta le lunghezze delle letture |
| `E-S1-07`, `E-S1-09` | un file non contiene le letture che il nome dichiara (inverse, o scambiate) | il dettaglio dice il caso e il file: sostituirlo, o scambiare i nomi dei due file |
| `E-S7-02` | con le coppie: nessuna sequenza fusa ha una lunghezza fra `pair.len_min` e `pair.len_max` | guardare `05_pair_merge/lunghezze_fuse.tsv` e correggere l'intervallo |
| `E-S13-04` | nessun campione resta nel risultato | guardare `12_final/intermedi/esclusioni.tsv`; abbassare `qc.min_reads_final` |
| `run` dice che la cartella non è vuota | `validate` o un'esecuzione vi ha già scritto | cancellarla (passo 7), oppure usare `resume` |
| l'esecuzione finisce la memoria (il processo è terminato dal sistema) | troppi processi per la memoria disponibile | scrivere `threads: 4` sotto `run:` e rilanciare `resume` |
