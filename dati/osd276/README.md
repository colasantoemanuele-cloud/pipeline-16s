# OSD-276: il secondo dataset

Questa cartella contiene cio' che serve a eseguire la pipeline su OSD-276, "Swabs
from the International Space Station" (NASA GeneLab/OSDR): 15 tamponi di superfici
della Stazione Spaziale Internazionale, regione V4 del gene 16S rRNA (primer
515F-806R), Illumina MiSeq. E' il dataset su cui si verifica che la pipeline non
presupponga la struttura di OSD-734, il dataset di riferimento del progetto.

## Le letture sono paired-end: si usano le sole forward

OSD-276 non e' un dataset single-end. E' depositato come PAIRED, e ENA distribuisce
per ciascuna delle 15 corse **un solo file** che contiene le due letture di ogni
coppia: prima tutte le inverse (record il cui nome termina con `/2`), poi tutte le
forward (`/1`). La pipeline tratta dati single-end: **qui si usano le sole letture
forward come dato single-end**. E' un uso legittimo (la lettura forward di un
amplicone V4 e' cio' che un sequenziamento single-end produrrebbe) ma va dichiarato:
i risultati non sono confrontabili con un'analisi paired-end dello stesso dataset, che
ricostruisce l'intero amplicone.

`scarica_letture.py` scarica il file di ENA, lo verifica con l'MD5 dichiarato da ENA e
ne ricava le forward senza modificarle. Le forward cosi' ricavate coincidono con i file
`GLDS-276_GAmplicon_<campione>_R1_raw.fastq.gz` pubblicati da OSDR (verificato sul
campione SP9: stesse sequenze e stesse qualita').

I file di ENA cosi' come sono **non vanno dati alla pipeline**: la fase di validazione
li accetterebbe, perche' giudica il layout dai nomi dei file, e la catena mescolerebbe
letture forward e inverse.

## Perche' questo dataset

Criteri di scelta, fissati prima di guardare i dati: un dataset 16S di NASA GeneLab,
preferibilmente depositato come SINGLE, diverso da OSD-734 nel maggior numero di
caratteristiche critiche (lunghezza delle letture, primer nelle letture, assenza di
controlli positivi, assenza del file del lotto, formato dei nomi dei file e degli
accession) e, a parita', il piu' piccolo.

Candidati esaminati il 7 ottobre 2026 su OSDR e su ENA:

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
| `scarica_letture.py` | si' | scarica, verifica e ricava le forward |
| `config_osd276.yaml` | si' | la configurazione della pipeline per questo dataset |
| `FONTI.tsv` | si' | fonte, licenza e controllo di ogni file |
| `ena/`, `fastq/` | no | i file di ENA (196 MB) e le forward ricavate (91 MB) |

Il riferimento tassonomico e' lo stesso di OSD-734 (SILVA 138) e si ottiene con
`dati/osd734/scarica_riferimento.py`.

## Come eseguire

Dalla radice del repository, con `<immagine>` l'immagine indicata nel README principale:

```bash
python3 dati/osd276/scarica_letture.py        # circa un minuto
python3 dati/osd734/scarica_riferimento.py    # se il riferimento non c'e' gia'

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
  amplicon16s run --config dati/osd276/config_osd276.yaml
```

Le uscite vanno in `output/osd276/`.

## La configurazione

`config_osd276.yaml` dichiara tutti i parametri che descrivono il dataset, come la
pipeline richiede. Quelli che per OSD-276 non hanno contenuto sono dichiarati vuoti o
nulli: le etichette dei controlli (`ctrl.blank_values: []`, `ctrl.positive_values: []`),
il taxon atteso e la colonna delle cellule dei controlli positivi, il file del lotto e
le sue colonne, l'espressione che deriva il modulo dalla posizione. La fase di
validazione lo registra: dichiara subito che il dataset non ha controlli positivi ne'
negativi (`E-S0-17`), cioe' che la soglia di profondita' non puo' venire da una curva e
che i contaminanti non sono stimabili.

## Stato della prova

L'esito della prova della pipeline su questo dataset (dove si fermava con i valori
predefiniti di OSD-734, e fin dove arriva con questa configurazione) e' descritto
nella sezione "Stato dell'implementazione" del README principale.
