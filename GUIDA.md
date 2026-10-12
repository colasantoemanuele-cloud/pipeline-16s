# Guida: dal proprio dataset al report, passo per passo

Per chi ha un dataset 16S, single-end o appaiato, e vuole scrivere la configurazione ed eseguire lo
strumento. Qui c'è solo il percorso; il significato di ogni parametro, i casi particolari e i
formati JSON sono in [RIFERIMENTO.md](RIFERIMENTO.md).

I messaggi dello strumento citano sezioni numerate (per esempio «GUIDA.md, 6.3») o un titolo
(«Letture appaiate»): sono in RIFERIMENTO.md, con lo stesso numero e lo stesso titolo, non nei passi
di questa guida; dove un argomento è trattato qui, RIFERIMENTO.md rimanda al passo giusto.

## Passo 1. Requisiti

- **Docker**. Dove non c'è (un cluster di calcolo) si usa Apptainer: RIFERIMENTO.md, 3.1.
- **L'immagine pubblicata**, scaricata per digest (1,86 GB da scaricare, 7,6 GB su disco):

  ```bash
  docker pull ghcr.io/colasantoemanuele-cloud/amplicon16s@sha256:930dda008581b24cc3a3fea718c0fed3690b011958bf1716c643f92056970154
  ```

- **Un clone di questo repository**: ogni comando si lancia dalla sua radice.

  ```bash
  git clone https://github.com/colasantoemanuele-cloud/pipeline-16s.git
  cd pipeline-16s
  ```

- **Memoria e spazio** in proporzione ai dati. Per orientarsi: 960 campioni e 21 milioni di letture
  chiedono 13 GB di memoria, 2,3 GB di disco e 50 minuti con 12 processori; qualche decina di
  campioni, una decina di minuti.

## Passo 2. Preparare i dati

Una cartella dei dati, fuori dal clone, con dentro:

- **Le letture**, in una cartella senza sottocartelle: file FASTQ (`.fastq.gz`, `.fq.gz`, `.fastq`,
  `.fq`), meglio se tutti con la stessa estensione. Il nome di ogni file contiene una parte unica
  per campione (un accession, un codice), che compare anche nei metadati.
  - Single-end: un file per campione.
  - Appaiate: due file per campione, con lo stesso nome salvo un marcatore in fondo, prima
    dell'estensione: `_R1` e `_R2`, `_1` e `_2`, `.1` e `.2` (anche `_R1_001`).
- **I metadati**: tabelle di testo UTF-8 separate da tabulazioni, con una riga di intestazione. O
  una sola tabella con una riga per campione (una colonna con il nome del file o la sua parte unica,
  una con il nome del campione, una con il tipo se ci sono controlli), o le tabelle ISA-Tab di assay
  e di studio.
- **Il file del lotto**, se c'è: una tabella con piastra di estrazione e corsa di sequenziamento per
  campione. È facoltativo; che cosa deve contenere è in RIFERIMENTO.md, 1.4.
- **Il riferimento tassonomico.** Dal clone, `python3 dati/riferimento/scarica_riferimento.py`
  scarica SILVA 138 (138 MB) in `dati/riferimento/`: nel container è
  `/app/dati/riferimento/silva_nr99_v138_train_set.fa.gz`. È il riferimento del classificatore
  predefinito (bayesiano ingenuo), la prima scelta per un'analisi al livello di genere: l'altro
  classificatore, IdTaxa, con SILVA 138.2 non ha superato la prova sulle comunità a composizione
  nota. Per un altro riferimento, o per IdTaxa, vedi RIFERIMENTO.md, 1.5.

Che cosa non mettere nella cartella delle letture:
- file di campioni che non sono nei metadati: la validazione si ferma (`E-S0-06`);
- file di un altro amplicone (ITS, 18S): se i metadati non dichiarano per riga l'amplicone o
  il primer, nulla li ferma;
- per usare le sole forward di una corsa appaiata, i file delle seconde letture: con i due file di
  ogni campione `init` sceglie la modalità appaiata;
- file con le due letture di ogni coppia insieme, o seconde letture senza il compagno: `init` si
  rifiuta e non scrive nulla. I primi vanno separati in due file.

## Passo 3. Un modo solo di lanciare i comandi

Si definisce una volta, nel terminale, dalla radice del clone, il comando `p16s`:

```bash
IMMAGINE=ghcr.io/colasantoemanuele-cloud/amplicon16s@sha256:930dda008581b24cc3a3fea718c0fed3690b011958bf1716c643f92056970154
DATI=/percorso/dei/dati        # la cartella dei dati del passo 2
LAVORO=/percorso/di/lavoro     # una cartella nuova, per configurazione e uscite
mkdir -p "$LAVORO"
p16s() { docker run --rm -u "$(id -u):$(id -g)" -v "$(pwd)":/app -v "$DATI":/dati:ro \
           -v "$LAVORO":/lavoro -w /app "$IMMAGINE" python3 scripts/esegui.py "$@"; }
```

Da qui ogni passo è una riga che comincia con `p16s`, lanciata dalla radice del clone. Nel container
la cartella dei dati si chiama `/dati` (in sola lettura), quella di lavoro `/lavoro`, il clone
`/app`: nei comandi e nella configurazione si scrivono questi nomi, non i percorsi del proprio
disco. La cartella di lavoro va creata prima (`mkdir -p`): se manca la crea Docker, di un altro
utente, e `init` non può scrivervi. Docker va lanciato sempre con `-u "$(id -u):$(id -g)"`, come fa
`p16s`: senza, i file scritti sono di root. In un terminale nuovo queste righe si ripetono. Dove un
messaggio dello strumento dice di lanciare `python3 scripts/esegui.py ...`, si scrive `p16s ...`.

## Passo 4. `init --accetta-proposte`

```bash
p16s init --reads /dati/letture --metadata /dati/campioni.tsv \
  --reference /app/dati/riferimento/silva_nr99_v138_train_set.fa.gz \
  --out-root /lavoro/esecuzione --output /lavoro/configurazione.yaml --accetta-proposte
```

- Con ISA-Tab `--metadata` è la cartella con le tabelle `a_*.txt` e `s_*.txt`; se le due tabelle
  hanno altri nomi si indicano entrambe, e il file del lotto con `--batch-table`: `--metadata
  /dati/assay.txt --study-table /dati/studio.txt --batch-table /dati/lotto.tsv`.
- `init` legge letture e metadati e scrive `/lavoro/configurazione.yaml`: ciò che ha rilevato e ciò
  che propone sono righe attive, le soglie di qualità senza un valore di partenza sono al valore che
  non ferma (si fissano al passo 8), le scelte che restano sono righe in commento. La modalità
  (single-end o appaiata) la riconosce dai nomi dei file.
- Dura qualche decina di secondi con decine di campioni, alcuni minuti con centinaia. Sulla console
  e in testa al file c'è il riepilogo: va letto. Ogni riga ha un commento che dice da dove viene il
  valore e come sceglierlo. Le righe lasciate da scegliere con una nota nel riepilogo vanno lette
  per prime: `init` le ha lasciate perché i dati erano ambigui, e il commento dice che cosa raccomanda.
- **Il primer.** Se `init` dice di non aver trovato il primer e il primer è noto, anche se è già
  stato tolto dalle letture, si rilancia con `--primer SEQUENZA` (in modalità appaiata anche
  `--primer-reverse SEQUENZA`) e un altro `--output`: serve a riconoscere la regione amplificata e i
  file di letture inverse. Se se ne conosce solo il nome, il commento di `qc.primer_sequence` nel
  file elenca le sequenze di quel nome; fra due forme si prende quella con più codici degenerati,
  che comprende l'altra (515F `GTGYCAGCMGCCGCGGTAA`, 806R `GGACTACNVGGGTWTCTAAT`; RIFERIMENTO.md,
  6.3 e 2, «In modalità appaiata»). Se non è noto si prosegue: `qc.primer_sequence` resta `null`.
- `init` non scrive su un file che esiste. Dopo un secondo `init` si cancella il primo file, e nei
  comandi che seguono si usa il nome del nuovo. Tutti gli argomenti: RIFERIMENTO.md, 2, «Il comando».

## Passo 5. Le scelte da scrivere

```bash
p16s validate --config /lavoro/configurazione.yaml
```

`validate` elenca i parametri che mancano (`E-G15-10`). Ciascuno è nel file come riga in commento
senza valore, sotto il suo gruppo: `filt.exclude_taxa` è la riga `  # exclude_taxa:` sotto `filt:`.
Si toglie `# `, lasciando i due spazi, e si scrive il valore: le stringhe fra virgolette, gli
elenchi in linea (`["a", "b"]`), i logici `true` o `false`. Le sezioni dell'ultima colonna sono di
RIFERIMENTO.md.

| Parametro | Che cosa decide | Regola semplice | Sezione |
|---|---|---|---|
| `tax.ref_name`, `tax.ref_version` | nome e versione del riferimento, riportati con i risultati | quelli della fonte: `"SILVA"` e `"138"` per il file del passo 2. `init` li scrive da sé se si leggono dal nome del file | 6.3 |
| `filt.remove_na_phylum` | se togliere le varianti a cui il riferimento non dà un phylum | `true` se ciò che il riferimento non riconosce è fuori bersaglio (il caso comune); `false` per tenerle | 6.4 |
| `filt.exclude_taxa` | i taxa da togliere, a qualunque rango | di solito `["Chloroplast", "Mitochondria", "Eukaryota"]`, scritti come li scrive il riferimento (così in SILVA 138); `[]` per non togliere nulla | 6.4 |
| `prev.apply` | se togliere le varianti presenti in pochi campioni | `false` se interessa la composizione di ogni singolo campione; `true` per togliere le varianti viste in un campione solo, spesso errori, sapendo che escono anche i taxa rari | 6.4 |
| `prev.min_fraction`, `prev.min_count` | solo con `prev.apply: true`: in quale frazione dei campioni, e con quante letture, una variante deve comparire | frazione più bassa con campioni numerosi o diversi fra loro (con 20 campioni e `0.2` servono 4 campioni); conteggio `1` conta ogni lettura, `2` ignora le letture singole | 6.4 |
| `qc.min_reads_final` | le letture sotto le quali un campione esce dall'oggetto finale | `0` non toglie nessuno: si parte da lì, si guardano nel report le letture per campione, si alza il valore e si lancia `resume` | 6.4 |
| `katharoseq.target_taxon` | solo con controlli positivi: il taxon atteso nei controlli | il genere del ceppo usato, come lo scrive il riferimento: lo sa chi ha allestito i controlli | 6.5 |
| `katharoseq.target_sensitivity` | solo con controlli positivi: la fedeltà (frazione di letture sul taxon atteso) a cui si legge la profondità minima dei campioni | strettamente fra 0 e 1; più alta, più campioni poveri esclusi. Non c'è un valore buono per tutti: si sceglie quanta contaminazione accettare, e nel report si guardano la soglia per piastra e i campioni esclusi | 6.5 |
| `katharoseq.min_r2` | solo con controlli positivi: la bontà minima (R²) perché una curva di calibrazione sia valida | fra 0 e 1; è una scelta: nel report si guardano gli R² delle curve, e una piastra sotto il minimo usa una soglia di ripiego | 6.5 |
| `ctrl.min_positive_pass_frac` | solo con controlli positivi: la frazione minima di controlli conformi al loro livello | fra 0 e 1; sotto, un avviso (`E-S11-04`), o un arresto se `ctrl.positive_gate` è `true`. `0` spegne il controllo: va bene finché la conformità non è nota | 6.5, 6.6 |
| `decontam.threshold` | solo con controlli negativi: la soglia sotto la quale una variante è un contaminante | strettamente fra 0 e 1. `0.1` è il predefinito di `decontam`; un valore più alto toglie di più e si usa con bassa biomassa | 6.5 |
| `decontam.min_blanks` | solo con controlli negativi: quanti controlli negativi con letture servono per decontaminare, e perché una piastra abbia un confronto proprio | non più dei controlli negativi del dataset. Una piastra che ne ha meno dà un avviso (`E-S0-15`) e si confronta con i negativi di tutte | 6.5 |

Le righe attive da controllare sempre, perché `init` può averle lasciate a un valore che dichiara
un'assenza (`null` o `[]`; il riepilogo le elenca):

| Parametro | Che cosa succede se resta così | Che cosa scrivere |
|---|---|---|
| `io.accession_regex`, `meta.accession_column` | se `init` non le ha rilevate restano in commento, e il riepilogo dice perché (nessuna colonna corrisponde ai file, o righe collegate al file sbagliato: passo 11) | la colonna che riporta il file di ogni campione, e un'espressione che estrae dal nome del file la parte unica per campione: `"^(CAMP[0-9]+)_"` (RIFERIMENTO.md, 1.1 e 6.1) |
| `meta.assay_keep_column`, `meta.assay_keep_values` | proposte solo se la tabella descrive anche un altro amplicone: tengono le sole righe 16S | controllare che i valori siano quelli delle righe da analizzare (RIFERIMENTO.md, 1.3) |
| `ctrl.column` e le tre etichette (`ctrl.blank_values`, `ctrl.positive_values`, `ctrl.biological_values`) | `null`: tutti i campioni sono biologici. La proposta guarda le parole dei valori e può sbagliare: il commento elenca i campioni di ogni valore | la colonna del tipo di campione e, per ogni etichetta, i valori di quella colonna; ogni campione deve ricadere in una sola. Colture pure e standard senza diluizioni vanno fra i biologici (RIFERIMENTO.md, 2.3, 6.1 e 6.5) |
| `ctrl.blank_override_column`, `ctrl.blank_override_values` | `null`: nessun campione è riclassificato | se alcuni campioni dichiarati biologici sono in realtà controlli negativi (materiale mai usato), la colonna e i valori che li riconoscono (RIFERIMENTO.md, 6.1) |
| `meta.module_column` | `null`: l'oggetto finale non ha la posizione di prelievo | la colonna della posizione, se lo studio ne ha una: nell'oggetto diventa `posizione`, e non va ripetuta in `out.study_columns` (RIFERIMENTO.md, 6.1) |
| `out.study_columns` | `[]`: l'oggetto finale non ha le altre variabili dello studio, e le analisi ecologiche non possono confrontarle | le colonne dei metadati con le variabili da confrontare, con il nome esatto dell'intestazione: `["Trattamento", "Tempo"]` (RIFERIMENTO.md, 6.1) |
| `qc.primer_sequence` | `null`: nessun primer è cercato nelle letture; se le letture lo portano entra nelle varianti | vedi il passo 4 (RIFERIMENTO.md, 6.3) |
| `filter.truncLen`, in modalità appaiata anche `pair.truncLen_reverse`, `pair.len_min`, `pair.len_max` | proposte dalla lunghezza più frequente e dalla qualità; il commento riporta lunghezze e qualità. Se mancano, vedi il passo 11 | (RIFERIMENTO.md, 6.3 e 6.7) |

Le altre righe a `null` o `[]` si lasciano così se lo studio non ha quell'informazione: il commento
di ciascuna dice che cosa dichiara (RIFERIMENTO.md, 6.1 e 6.2).

## Passo 6. `validate`, finché passa

Si rilancia il comando del passo 5 dopo ogni correzione, finché scrive «Validazione superata».
L'avviso sui controlli di qualità spenti è atteso: si chiude al passo 8. Gli altri avvisi
(`WARNING`) non fermano, ma vanno letti: senza controlli negativi o positivi è atteso `E-S0-17`, e
l'esecuzione dirà poi che non decontamina o non calibra; `E-S0-32` dice che un taxon di
`filt.exclude_taxa` non è nel riferimento con quel nome, e quindi non toglie nulla; `E-S0-31`
elenca le righe dei metadati escluse perché di un altro amplicone.

## Passo 7. `run` e `status`

```bash
rm -r "$LAVORO/esecuzione"      # validate vi ha scritto, e run parte solo da una cartella vuota
p16s run --config /lavoro/configurazione.yaml
p16s status --out /lavoro/esecuzione
```

Al posto di cancellare la cartella e lanciare `run` si può lanciare `resume`. `status` si può
lanciare durante `run` da un altro terminale, dopo avervi ripetuto le righe del passo 3. Se `run` si
ferma dice la fase, il codice e che cosa fare (passo 11). Corretta la causa, si riparte con `p16s
resume --config /lavoro/configurazione.yaml`, che rifà solo ciò che serve. I risultati sono in
`/lavoro/esecuzione/12_final/` (`ps_final.rds` e gli stessi dati in testo).

## Passo 8. `soglie --accetta-proposte` e `resume`

```bash
p16s soglie --config /lavoro/configurazione.yaml --output /lavoro/con_soglie.yaml --accetta-proposte
p16s resume --config /lavoro/con_soglie.yaml
```

`soglie` legge le misure dell'esecuzione conclusa e scrive un nuovo file, uguale al primo salvo le
soglie di qualità, ciascuna con la sua misura accanto. `resume` con il nuovo file le attiva: avvisa
che la configurazione è cambiata, ripete i calcoli (dura quanto `run`) e non cambia i risultati. Una
soglia proposta protegge dai peggioramenti, non giudica se la misura è buona: una misura che
sorprende va capita prima di accettarla (che cosa misura ogni soglia: RIFERIMENTO.md, 6.6; per le
coppie fuse, 6.7). Da qui in poi la configurazione è `con_soglie.yaml`.

## Passo 9. `report`

```bash
p16s report --config /lavoro/con_soglie.yaml
```

Scrive `$LAVORO/esecuzione/report/report.html`, che si apre con un browser. Come si legge è nel
README, «Come si legge il report»; le tre cose da guardare per prime:

1. **Segnalazioni in apertura**: i controlli di qualità spenti e il numero delle degradazioni, cioè
   degli avvisi con cui l'esecuzione ha proseguito rinunciando a qualcosa; il loro elenco è in
   «Decisioni prese automaticamente».
2. **Tracciamento delle letture**: quante letture restano dopo ogni fase; una perdita forte in una
   fase dice dove guardare (in modalità appaiata, la quota di coppie fuse).
3. **Risultato finale**: campioni e varianti rimasti. Quali campioni sono usciti, e perché, è in
   «Decisioni prese automaticamente», alla voce dei filtri finali.

## Passo 10. Analisi ecologiche

```bash
p16s eco init --object /lavoro/esecuzione/12_final/ps_final.rds --output /lavoro/eco.yaml --accetta-proposte
p16s eco validate --object /lavoro/esecuzione/12_final/ps_final.rds --config /lavoro/eco.yaml
p16s eco run --object /lavoro/esecuzione/12_final/ps_final.rds --config /lavoro/eco.yaml --out /lavoro/analisi
```

`eco init` scrive `eco.yaml` ed elenca le colonne dell'oggetto con i loro valori distinti: fra le
altre `posizione`, `piastra` e `corsa` (con 0 valori se non dichiarate) e quelle portate con
`out.study_columns`, con il nome in minuscolo e i trattini bassi (`Factor Value[Tempo]` diventa
`factor_value_tempo`). Le proposte attive (variabili tecniche, rango) si controllano. `eco validate`
nomina le scelte che mancano (`E-ECO-02`); si scrivono come al passo 5:

| Parametro | Che cosa decide | Regola semplice |
|---|---|---|
| `design.variable` | la colonna con i gruppi da confrontare | una delle colonne elencate da `eco init`, che la propone solo se una sola ha gruppi ripetuti. Va dichiarata anche se nessuna ne ha: i test non si eseguono, tabelle e grafici sì. Se la colonna giusta manca dall'oggetto, la si aggiunge a `out.study_columns` in `con_soglie.yaml` e si rilanciano `resume` (passo 8) ed `eco init` con un altro `--output` |
| `run.seed` | il seme di rarefazione e permutazioni | un intero qualunque |
| `comp.rank` | il rango a cui si riassume la composizione | di solito proposto: il più fine con almeno metà delle varianti assegnate; altrimenti uno dei ranghi elencati da `eco init` |
| `comp.top_n` | quanti taxa distinguere nei grafici di composizione | una decina: di più non si distinguono i colori |
| `beta.distances` | le distanze fra campioni | `["bray", "jaccard"]`. Con `"aitchison"` si scrive anche `beta.clr_pseudocount` (`1` è la scelta più diffusa); per le UniFrac serve un albero: `config/eco.example.yaml` |
| `ord.methods`, `ord.distances` | le ordinazioni, e su quali distanze | `["pcoa"]`, `["nmds"]` o entrambe; le distanze fra quelle di `beta.distances`. Con `pcoa` si toglie il commento a `ord.pcoa_correction`, con `nmds` a `ord.nmds_trymax`, che hanno già un valore |
| `stat.min_group_size` | i campioni che un gruppo deve avere per entrare nei test | non meno di `2` |
| `stat.significance_level` | il livello di significatività dello studio | di consueto `0.05` |
| `design.subset` | facoltativo: quali campioni entrano nelle analisi | toglie i campioni fuori dal confronto, per esempio una coltura pura tenuta come biologica: una colonna e i suoi valori da tenere (RIFERIMENTO.md, 4) |

`eco run` vuole una cartella di uscita assente o vuota. In `/lavoro/analisi` scrive tabelle (`.tsv`)
e grafici (`.png`) nelle sottocartelle `alfa`, `composizione`, `beta`, `ordinazione` e `test`, e
`riepilogo.json` con gli avvisi (`E-ECO-`): non fermano, ma dicono che cosa non si è potuto
calcolare (un gruppo troppo piccolo per i test, una NMDS con pochi campioni). Per ogni altra
scelta, le variabili continue e le misure ripetute: RIFERIMENTO.md, sezione 4.

## Passo 11. Se qualcosa va storto

Ogni arresto porta un codice, il dettaglio del caso e che cosa fare: va letto per intero, e basta
anche per i codici che qui non ci sono. I più comuni (le sezioni sono di RIFERIMENTO.md):

| Codice | Che cosa vuol dire | Che cosa fare | Sezione |
|---|---|---|---|
| `E-G15-10` | mancano parametri obbligatori | scrivere le scelte elencate (passo 5) | 6 |
| `E-G15-99` | un valore non è ammesso, o una chiave è scritta male | correggere la riga che il messaggio nomina | 6 |
| `E-S0-06`, `E-S0-04` | file e metadati non corrispondono: un file senza riga, una riga senza file, una chiave che non si estrae dal nome | togliere il file estraneo o completare i metadati; controllare `io.accession_regex` e `meta.accession_column` | 1.1, 6.1 |
| `E-S0-30`, o `init` che non propone la chiave | i metadati del deposito sono sbagliati: una riga è collegata al file di un altro amplicone, o a un file che non c'è | copiare le tabelle in una cartella propria e correggere la copia, mai gli originali (ogni campione al suo file, via le righe estranee), poi rilanciare `init` sulla copia. Va verificata ogni riga, non solo quelle segnalate: uno scambio fra campioni dello stesso amplicone non si vede dalle letture | 1.3 |
| `E-S0-11` | un campione non ricade in nessuna etichetta | aggiungere il suo valore a una delle tre etichette `ctrl.*_values` | 6.1 |
| `E-S0-12` | il riferimento manca o il suo MD5 non è quello dichiarato | controllare il percorso; `tax.ref_md5` è l'uscita di `md5sum FILE` | 1.5 |
| `E-S0-24`, `E-S0-23` | la cartella ha due file per campione ma la modalità è single-end; oppure, in modalità appaiata, a un file manca il compagno | rilanciare `init` (riconosce le coppie), o lasciare nella cartella i soli file forward; aggiungere il file che manca | 1.2 |
| `E-S0-09`, `E-S1-02` | troppe letture sono più corte di `filter.truncLen` | abbassare `filter.truncLen`, oppure alzare `qc.max_frac_short_reads` se la perdita è accettata | 6.3 |
| `E-S1-07` | single-end: un file non contiene le sole letture forward | il dettaglio dice il caso: sostituire il file con quello delle forward; togliere il campione se non ha quasi nulla del bersaglio; riportare a un ordine non legato alla sequenza un file riordinato | 1.2 |
| `E-S1-09` | appaiata: i file di una corsa non hanno l'orientamento dichiarato | il dettaglio dice il caso; se la corsa è invertita, scambiare i nomi dei suoi due file, poi `resume` | 1.2 |
| `E-S7-02` | appaiata: nessuna sequenza fusa ha una lunghezza fra `pair.len_min` e `pair.len_max` | guardare `05_pair_merge/lunghezze_fuse.tsv` e correggere l'intervallo, o i troncamenti | 6.7 |
| `E-S13-04` | nessun campione resta nell'oggetto finale | guardare `12_final/intermedi/esclusioni.tsv`; abbassare `qc.min_reads_final`, rivedere `prev.min_fraction` | 6.4 |
| `E-RUN-01` | un'altra esecuzione usa la stessa cartella di uscita | aspettare che finisca (`status`), o cambiare `io.out_root` | 3 |
| nessun codice: `run` dice che la cartella non è vuota | `validate` o un'esecuzione vi ha già scritto | cancellare la cartella, oppure usare `resume` (passo 7) | 3 |
| nessun codice: «Permission denied», o cartella non scrivibile, scrivendo nella cartella di lavoro o di uscita | la cartella è stata creata da Docker, o da un comando lanciato senza `-u`: appartiene a root | crearla prima con `mkdir -p` (passo 3); se esiste già, toglierla con `sudo rm -r` e ricrearla | 2, «Il comando» |
| nessun codice: in modalità appaiata `init` non propone `filter.truncLen`, `pair.truncLen_reverse`, `pair.len_min`, `pair.len_max` | non ha potuto misurare la sovrapposizione delle due letture né riconoscere i due primer (code con molti errori, o letture troppo corte per la regione) | se le letture coprono la regione, scegliere a mano i due troncamenti dalle lunghezze e dalla qualità che il commento delle due righe riporta, rispettando il vincolo della fusione, e l'intervallo dalla lunghezza della regione; altrimenti usare le sole forward | 2, «In modalità appaiata»; 6.7 |
| nessun codice: `init` dice che le coppie non si fondono senza basi di qualità bassa | le letture reverse sono corte o di qualità bassa: per fonderle i troncamenti dovrebbero entrare in basi di qualità mediana sotto 20 | usare le sole forward: togliere i file reverse dalla cartella delle letture e rilanciare `init` con un altro `--output` | 2, «In modalità appaiata»; 6.3 |
