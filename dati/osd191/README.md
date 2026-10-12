# OSD-191: dati e scarico

OSD-191 (NASA GeneLab/OSDR, «Temporal dynamics of the gut microbiota in people sharing a
confined environment, a 520-day ground-based space simulation, MARS500») raccoglie 159
campioni di feci di sei soggetti, prelevati prima, durante e dopo 520 giorni di
isolamento e sequenziati con Illumina MiSeq a letture appaiate di 300 basi. La regione
amplificata e' la V3-V4 del gene 16S rRNA (primer 341F `CCTACGGGNGGCWGCAG` e 805R
`GACTACHVGGGTATCTAATCC`): circa 465 basi con i primer, piu' di una lettura. Non ha
controlli ne' un file del lotto.

Questa cartella contiene i metadati e cio' che serve a ottenere e verificare le letture.
**Non contiene una configurazione della pipeline**: su queste letture `init` non arriva
a proporre i troncamenti (vedi sotto), che vanno scelti da chi conduce l'analisi.

## Contenuto della cartella

| File | Che cos'e' |
|---|---|
| `metadati/` | le tabelle di assay e di studio dell'archivio ISA di OSDR |
| `letture_osdr.tsv` | per ogni campione gli indirizzi, le dimensioni e gli MD5 dei due file grezzi |
| `scarica_letture.py` | scarica i due file di ogni campione in `fastq_coppie/` |
| `abbinamento_genelab.tsv` | per ogni campione, il nome con cui compare nella tabella elaborata da GeneLab |
| `FONTI.tsv` | fonte, licenza e controllo di ogni file |

Le letture (`fastq_coppie/`, 6,6 GB) non sono nel repository.

## Le letture

OSDR pubblica per ogni campione due file, R1 e R2:

```bash
python3 dati/osd191/scarica_letture.py --coppie                  # 318 file in dati/osd191/fastq_coppie/
python3 dati/osd191/scarica_letture.py --coppie --solo-verifica  # controlla i file presenti
```

I due file di ogni campione si salvano come `<campione>_R1.fastq.gz` e
`<campione>_R2.fastq.gz`, dove `<campione>` e' il `Sample Name` della tabella di assay,
nell'ordine del deposito e senza scambiarli. OSDR dichiara di ogni file la sola
dimensione: lo scarico la verifica, e verifica l'MD5 registrato in `letture_osdr.tsv`.
Senza `--coppie` lo script scarica in `fastq/` il solo file forward di ogni campione.

OSDR pubblica gli stessi file anche con un altro nome
(`GLDS-191_GAmplicon_<campione>_R1_raw.fastq.gz`): hanno le stesse dimensioni, e qui si
usano quelli indicati nella colonna `Raw Data File` della tabella di assay.

## Che cosa propone `init` su queste letture

Le letture portano i primer: il 341F in testa alle forward e l'805R in testa alle
reverse, che `init` rileva (tagli di 17 e 21 basi). Le lunghezze vanno da 35 a 301
basi, con la moda a 301: una piccola parte delle letture e' stata accorciata prima del
deposito.

- **In modalita' appaiata `init` non stima la lunghezza dell'inserto**, quindi non
  propone `filter.truncLen`, `pair.truncLen_reverse`, `pair.len_min` e `pair.len_max`:
  misura la sovrapposizione sulle letture intere, e nelle code di queste letture di 300
  basi gli errori sono troppi perche' la riconosca (poche coppie su cento). Le due
  letture coprono comunque la regione: i due troncamenti si scelgono dal profilo di
  qualita' rispettando il vincolo della fusione, e l'intervallo delle lunghezze dalla
  regione amplificata (GUIDA.md, sezione 6.7).
- **Sulle sole forward `init` propone un troncamento di 35 basi**, la lunghezza della
  lettura piu' corta: la sua proposta e' la lunghezza maggiore che non scarta alcuna
  lettura, e qui va sostituita con una scelta fatta sul profilo di qualita', accettando
  di scartare le poche letture piu' corte (GUIDA.md, sezione 6.3).

## Confronto con l'elaborazione di GeneLab

GeneLab pubblica per questo dataset la tabella `GLDS-191_GAmplicon_taxonomy-and-counts.tsv`
(pagina del dataset su OSDR), ottenuta togliendo i primer con cutadapt, troncando le
letture a 250 e 200 basi, fondendo le coppie e assegnando la tassonomia con IdTaxa
(archivio `GLDS-191_GAmplicon_processing_info.tar`). L'oggetto finale di un'esecuzione
si confronta con quella tabella, a livello di genere, con
`dati/confronto_genelab/confronta.py`. Nella tabella di GeneLab i campioni portano il
`Sample Name` della tabella di assay: `abbinamento_genelab.tsv` ha lo stesso nome nelle
due colonne.
