# OSD-276: il secondo dataset

Questa cartella contiene cio' che serve a eseguire la pipeline su OSD-276, "Swabs
from the International Space Station" (NASA GeneLab/OSDR): 15 tamponi di superfici
della Stazione Spaziale Internazionale, regione V4 del gene 16S rRNA (primer
515F-806R), Illumina MiSeq. E' il dataset su cui si verifica che la pipeline non
presupponga la struttura di OSD-734, il dataset di riferimento del progetto.

## Letture appaiate in un solo file per corsa: prima vanno separate

OSD-276 e' depositato come PAIRED, e ENA distribuisce per ciascuna delle 15 corse **un
solo file** che contiene le due letture di ogni coppia, in due blocchi: i record la cui
intestazione termina con `/1` e quelli che terminano con `/2`; in 8 corse viene prima il
blocco `/1`, in 7 il blocco `/2`. La pipeline non accetta un file cosi' (vedi sotto): le
due letture vanno prima separate, e lo fa `scarica_letture.py`. La configurazione di base
di questa cartella usa le **sole letture forward**, in modalita' single-end; la modalita'
appaiata (`io.layout: paired-end`) si usa sui due file per corsa che lo script ricava con
`--coppie`. Le sole forward sono un uso legittimo (la lettura forward di un amplicone V4
e' cio' che un sequenziamento single-end produrrebbe) ma va dichiarato: i risultati non
sono confrontabili con un'analisi appaiata dello stesso dataset, che ricostruisce
l'intero amplicone.

`scarica_letture.py` scarica il file di ENA, lo verifica con l'MD5 dichiarato da ENA e
ne ricava le forward (i record `/1`) senza modificarle. Prima conta i record `/1`, i
record `/2` e gli altri, e rifiuta il file se ce ne sono di altri, se `/1` e `/2` non
sono in numero uguale o se la somma non e' il numero di letture dichiarato dal
deposito; un file forward con un'impronta diversa dall'attesa viene messo da parte
(`.md5_errato`) e non resta in `fastq/`; il guasto di una corsa finisce nel riepilogo
finale senza fermare le altre.

Con `--coppie` lo script ricava entrambe le letture, per la modalita' appaiata della
pipeline (`io.layout: paired-end`), in una cartella separata: per ogni corsa
`fastq_coppie/<corsa>_<campione>_R1.fastq.gz` (i record `/1`) e `..._R2.fastq.gz` (i
record `/2`), verificati contro `letture_forward.md5` e `letture_reverse.md5`. I due
file prendono il nome dall'etichetta del deposito e non vengono scambiati; `fastq/` non
viene toccata.

```bash
python3 dati/osd276/scarica_letture.py --coppie   # 30 file in dati/osd276/fastq_coppie/
```

### I file pubblicati da OSDR

OSDR pubblica per ogni campione due file,
`GLDS-276_GAmplicon_<campione>_R1_raw.fastq.gz` e `..._R2_raw.fastq.gz`. Ogni blocco dei
file di ENA coincide con uno dei due file di OSDR, ma **il nome R1 non individua la
lettura forward**: l'etichetta `/1` di ENA individua in ogni corsa la lettura che parte
dal lato del primer 515F, quella che la pipeline si attende (lo conferma il motivo
conservato a valle del primer, presente nel solo blocco `/1`); il file R1 di OSDR in
alcune corse e' la lettura forward e in altre la reverse. `corse_invertite_osdr.txt`
elenca i campioni in cui i due file di OSDR sono scambiati: usare per tutti il file R1
di OSDR mescolerebbe le due estremita' dell'amplicone, e per quei campioni la tabella
elaborata da GeneLab non e' confrontabile (`dati/confronto_genelab/`).

I file di ENA cosi' come sono **non vanno dati alla pipeline**, che li respinge
tutti e 15. La fase di validazione riconosce un file con le due letture di ogni
coppia dalle intestazioni delle letture che ispeziona (le prime `qc.head_reads`,
10.000 nella configurazione di questa cartella): in questi file ogni blocco e' piu' lungo (da 35.489 a 105.423
letture), e le letture ispezionate sono tutte dello stesso blocco. I 7 file che
cominciano con il blocco `/2` sono respinti dalla validazione (`E-S0-07`: sole
seconde letture); gli 8 che cominciano con il blocco `/1` sembrano, nelle letture
ispezionate, file di sole forward, e li respinge il profilo di qualita', che conta i
segni di coppia su tutte le letture di ogni file (`E-S1-04`, con i conteggi in
`02_qc_profiles/coppie.tsv`). Accettati, la catena mescolerebbe letture forward e
inverse.

## Perche' questo dataset

Criteri di scelta, fissati prima di guardare i dati: un dataset 16S di NASA GeneLab,
preferibilmente depositato come SINGLE, diverso da OSD-734 nel maggior numero di
caratteristiche critiche (lunghezza delle letture, primer nelle letture, assenza di
controlli positivi, assenza del file del lotto, formato dei nomi dei file e degli
accession) e, a parita', il piu' piccolo.

Candidati esaminati su OSDR e su ENA:

| Che cosa | OSD-276 | OSD-146 |
|---|---|---|
| Deposito | PAIRED; un file per corsa con le due letture | PAIRED; due file per corsa |
| Strumento | Illumina MiSeq | Illumina HiSeq 2500 |
| Lunghezza delle letture | 151 basi, tutte | 148 e 150 basi |
| Primer in testa alle letture | no (0 su 1.000 letture) | no (0 su 1.000 letture) |
| Corse e dimensione su ENA | 15 corse, 205 MB | 80 corse, 2,1 GB |
| Metadati ISA | si' | si' |
| Controlli positivi e negativi | nessuno | nessuno |
| File del lotto | nessuno | nessuno |
| Nomi dei file e accession | per corsa (`SRR`) | per corsa (`SRR`); su OSDR per campione, senza accession |

Nessuno dei due e' depositato come SINGLE, e nessuno ha il primer nelle letture: si
usano quindi le letture forward di un dataset paired-end. I due candidati differiscono
da OSD-734 nelle stesse caratteristiche (lunghezza delle letture, nessun controllo
positivo, nessun file del lotto, nomi e accession per corsa invece che per
esperimento): a parita' e' stato scelto OSD-276, il piu' piccolo.

In che cosa OSD-276 differisce da OSD-734:

| Caratteristica | OSD-734 | OSD-276 |
|---|---|---|
| Campioni | 960 (770 biologici, 80 positivi, 110 negativi) | 15, tutti biologici |
| Lunghezza delle letture | da 137 a 151 basi | 151 basi |
| Controlli | serie KatharoSeq e bianchi | nessuno |
| File del lotto, piastre, corse | 10 piastre, 2 corse | nessun file; una corsa |
| Chiave dei campioni | accession dell'esperimento (`ERX`) | accession della corsa (`SRR`) |
| Materiale dei campioni biologici | `Surface swab` | `Cells` |
| Posizione | sigla con il modulo (`LAB1P3`) | descrizione libera (`lab fwd air vent`) |

## Contenuto

| File | Nel repository | Che cosa e' |
|---|---|---|
| `metadati/` | si' | le tabelle ISA di OSDR: campioni di studio e assay 16S |
| `letture_ena.tsv` | si' | le 15 corse con indirizzo, MD5 e dimensione dichiarati da ENA, e il numero di letture forward attese |
| `letture_forward.md5` | si' | l'MD5 del contenuto non compresso di ogni file forward ricavato |
| `letture_reverse.md5` | si' | l'MD5 del contenuto non compresso di ogni file reverse ricavato con `--coppie` |
| `scarica_letture.py` | si' | scarica, verifica e ricava le forward; con `--coppie`, le due letture |
| `config_osd276.yaml` | si' | la configurazione della pipeline per questo dataset |
| `config_osd276_forward_idtaxa.yaml` | si' | la stessa con IdTaxa e il classificatore SILVA SSU r138 |
| `eco_osd276.yaml` | si' | la configurazione d'esempio delle analisi ecologiche di base |
| `corse_invertite_osdr.txt` | si' | i campioni in cui i file R1 e R2 pubblicati da OSDR sono scambiati rispetto alle letture forward e reverse |
| `abbinamento_genelab.tsv` | si' | per ogni campione, il nome con cui compare nella tabella elaborata da GeneLab (per `dati/confronto_genelab/`) |
| `FONTI.tsv` | si' | fonte, licenza e controllo di ogni file |
| `ena/`, `fastq/` | no | i file di ENA (196 MB) e le forward ricavate (91 MB) |
| `fastq_coppie/` | no | le due letture di ogni corsa, ricavate con `--coppie` |

Il riferimento tassonomico (SILVA 138) e' comune ai dataset: sta in
`dati/riferimento/` e si ottiene con `dati/riferimento/scarica_riferimento.py`.

## Come eseguire

Dalla radice del repository, con `<immagine>` l'immagine indicata nel README principale:

```bash
python3 dati/osd276/scarica_letture.py            # circa un minuto
python3 dati/riferimento/scarica_riferimento.py   # se il riferimento non c'e' gia'

docker run --rm \
  --memory=24g \
  --memory-swap=24g \
  -u "$(id -u):$(id -g)" \
  -v "$(pwd)":/app \
  -w /app \
  <immagine> \
  python3 scripts/esegui.py run --config dati/osd276/config_osd276.yaml
```

Le uscite vanno in `output/osd276/`.

La configurazione si puo' anche riscrivere da zero con `init`, come per un dataset
proprio (GUIDA.md): sulle letture e sulle tabelle di questa cartella rileva il modello
dei nomi, la colonna e l'espressione della chiave (`SRR[0-9]+`), la colonna del nome
del campione e l'MD5 del riferimento, e propone la colonna del tipo con la sua sola
etichetta; le scelte restano da scrivere, e quelle di questo dataset sono in
`config_osd276.yaml`. Una proposta differisce dalla configurazione di questa cartella:
`init` propone `filter.truncLen: 151`, la lunghezza di tutte le letture, mentre qui si
tronca a 150 per togliere l'ultima base (una scelta di questo esempio, non degli autori
dello studio). I due valori danno oggetti finali diversi: per riprodurre il risultato di
questa cartella il troncamento e' 150. Il primer non e' nelle letture, e i metadati ne
riportano solo il nome (515F-806R): la sequenza va scritta, o indicata con `--primer`.

Le analisi ecologiche di base sull'oggetto finale:

```bash
docker run --rm -u "$(id -u):$(id -g)" -v "$(pwd)":/app -w /app <immagine> \
  python3 scripts/esegui.py eco run --object output/osd276/12_final/ps_final.rds \
    --config dati/osd276/eco_osd276.yaml --out output/eco_osd276
```

## La configurazione

`config_osd276.yaml` dichiara tutti i parametri che descrivono il dataset, come la
pipeline richiede. Quelli che per OSD-276 non hanno contenuto sono dichiarati vuoti o
nulli: le etichette dei controlli (`ctrl.blank_values: []`, `ctrl.positive_values: []`),
il taxon atteso e la colonna delle cellule dei controlli positivi, il file del lotto e
le sue colonne, l'espressione che deriva il modulo dalla posizione. La fase di
validazione lo registra: dichiara subito che il dataset non ha controlli positivi ne'
negativi (`E-S0-17`), cioe' che la soglia di profondita' non puo' venire da una curva e
che i contaminanti non sono stimabili. Senza controlli positivi non si applica alcuna
soglia di profondita' (`E-S11-05`): sui campioni resta il solo `qc.min_reads_final`,
sulle letture dell'oggetto finale.

`config_osd276_forward_idtaxa.yaml` differisce da `config_osd276.yaml` solo nel
classificatore: assegna la tassonomia con IdTaxa e `SILVA_SSU_r138_2019.RData`, il
classificatore che GeneLab ha usato per la tabella elaborata di questo dataset
(`python3 dati/riferimento/scarica_addestrato.py 2019`). Dichiara la regola rigorosa
sulla provenienza (`run.strict_provenance: true`) e scrive in una cartella propria sotto
`output/`.

Gli avvisi attesi con questa configurazione sono quelli di un dataset senza controlli e
senza file del lotto: due `E-S0-17`, `E-S11-05` ed `E-S12-03`.
