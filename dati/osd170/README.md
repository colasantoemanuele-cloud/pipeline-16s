# OSD-170: dati, configurazione ed esecuzione

OSD-170 (NASA GeneLab/OSDR, «NASA Aircraft Bioaerosol Collector») e' uno dei dataset su
cui la pipeline e' eseguita: 30 campioni sequenziati con Illumina MiSeq a letture
appaiate, di cui 24 di bioaerosol raccolti a terra e in quota in quattro voli, 4 campioni
dell'attrezzatura di prelievo e 2 controlli negativi (filtri bianchi). Non ha controlli
positivi ne' un file del lotto. Differisce dagli altri due dataset per lunghezza delle
letture (251 basi), ambiente e disegno dei controlli.

## Contenuto della cartella

| File | Che cos'e' |
|---|---|
| `metadati/` | le tabelle di assay e di studio dell'archivio ISA di OSDR |
| `letture_osdr.tsv` | per ogni campione gli indirizzi e le dimensioni dei due file grezzi, il lato tenuto e il suo MD5 |
| `scarica_letture.py` | scarica il file forward di ogni campione in `fastq/` |
| `config_osd170.yaml` | la configurazione della pipeline, commentata |
| `eco_osd170.yaml` | la configurazione delle analisi ecologiche di base |
| `abbinamento_genelab.tsv` | per ogni campione, il nome con cui compare nella tabella elaborata da GeneLab |
| `FONTI.tsv` | fonte, licenza e controllo di ogni file |

Le letture (`fastq/`, 966 MB) non sono nel repository.

## Le letture

OSDR pubblica per ogni campione due file, R1 e R2. La pipeline tratta dati single-end e
usa il solo file forward:

```bash
python3 dati/osd170/scarica_letture.py                  # 30 file in dati/osd170/fastq/
python3 dati/osd170/scarica_letture.py --solo-verifica  # controlla i file presenti
```

Lo script (`dati/letture_forward.py`) legge le prime letture dei due file di ogni corsa e
tiene quello che inizia con il primer forward; una corsa in cui lo porta il secondo file
e' dichiarata invertita. **Le letture di questo dataset sono depositate senza primer**: lo
script lo dichiara per tutte le corse e tiene l'R1, senza poter verificare l'inversione
con questa regola. Il file tenuto e' salvato come `<campione>.fastq.gz`, senza marcatori
di coppia nel nome. OSDR dichiara di ogni file la sola dimensione: lo scarico la
verifica, e verifica l'MD5 registrato in `letture_osdr.tsv` al primo scarico.

**I metadati non descrivono le letture depositate.** La tabella di assay dichiara il
primer forward `TCCTACGGGAGGCAGCAGT` e letture di 220 basi; le letture non portano alcun
primer, sono di 251 basi e iniziano con il motivo conservato a valle del primer 515F: la
regione letta e' la V4. Nella configurazione `qc.primer_sequence` e' il primer dichiarato
(la validazione ne verifica l'assenza) e `qc.conserved_motif` il motivo ricavato dalle
letture; `filter.trimLeft` e' 0.

## Esecuzione

Dalla radice del clone, con le letture e il riferimento scaricati
(`python3 dati/riferimento/scarica_riferimento.py`):

```bash
docker run --rm -u "$(id -u):$(id -g)" -v "$(pwd)":/app -w /app <immagine> \
  python3 scripts/esegui.py run --config dati/osd170/config_osd170.yaml
docker run --rm -u "$(id -u):$(id -g)" -v "$(pwd)":/app -w /app <immagine> \
  python3 scripts/esegui.py report --config dati/osd170/config_osd170.yaml
docker run --rm -u "$(id -u):$(id -g)" -v "$(pwd)":/app -w /app <immagine> \
  python3 scripts/esegui.py eco run --object output/osd170/12_final/ps_final.rds \
    --config dati/osd170/eco_osd170.yaml --out output/osd170_eco
```

`<immagine>` e' quella indicata nel README principale. La catena impiega una decina di
minuti. I due controlli negativi sono consegnati in `ps_controlli.rds`; i quattro
campioni dell'attrezzatura restano fra i campioni (sono tamponi di superfici, non
bianchi). Senza controlli positivi nessuna soglia di profondita' si applica, e la
pipeline lo dichiara.

`config_osd170.yaml` dichiara anche un parametro che su questo dataset non ha effetto e
che si potrebbe omettere (`katharoseq.target_sensitivity`: non ci sono controlli
positivi): e' dichiarato perche' la configurazione dia esattamente il suo digest, che
identifica le esecuzioni.

## Confronto con l'elaborazione di GeneLab

GeneLab pubblica per questo dataset la tabella `GLDS-170_GAmplicon_taxonomy-and-counts.tsv`
(pagina del dataset su OSDR). Scaricata, si confronta con l'oggetto finale a livello di
genere con `dati/confronto_genelab/confronta.py`, che dice che cosa misura e con quali
soglie:

```bash
python3 dati/confronto_genelab/confronta.py --finale output/osd170/12_final \
  --genelab GLDS-170_GAmplicon_taxonomy-and-counts.tsv \
  --abbinamento dati/osd170/abbinamento_genelab.tsv --uscita output/osd170_confronto
```

## Sottoinsieme di prova

`tests/fixtures/osd170/ridotto/` contiene le prime 3.000 letture di sei campioni e dei due
controlli negativi, con un riferimento tassonomico sintetico e la propria configurazione:
la catena vi si esegue senza scaricare nulla
(`python3 scripts/esegui.py run --config tests/fixtures/osd170/ridotto/config.yaml`).
