# Riferimento tassonomico comune ai dataset

Questa cartella contiene cio' che serve a ottenere e verificare il riferimento
tassonomico usato dalle configurazioni del repository (`dati/osd734/config_osd734.yaml`,
`dati/osd276/config_osd276.yaml`): SILVA 138 nel formato per `assignTaxonomy` di dada2.
Il riferimento non appartiene a un dataset, e sta in una cartella propria perche' la
configurazione di un dataset non dipenda dalla cartella di un altro.

| File | Nel repository | Contenuto |
|---|---|---|
| `scarica_riferimento.py` | si' | scarica e verifica i tre file in questa cartella |
| `FONTI.tsv` | si' | fonte, licenza e controllo di ciascun file |
| `impronte.md5` | si' | gli MD5 dei tre file, nel formato di `md5sum -c` |
| `silva_nr99_v138_train_set.fa.gz` | no | il training set (`tax.ref_fasta`; l'MD5 e' `tax.ref_md5`), 137.973.851 byte |
| `SILVA_LICENSE.txt` | no | la licenza dichiarata dalla fonte |
| `silva_138_v2_bad-taxa.csv` | no | i taxa con un rango mancante nella versione 2 (`tax.ref_bad_taxa`) |

```bash
python3 dati/riferimento/scarica_riferimento.py                  # scarica, circa 138 MB da Zenodo
python3 dati/riferimento/scarica_riferimento.py --solo-verifica  # controlla e basta
(cd dati/riferimento && md5sum -c impronte.md5)                  # verifica indipendente dallo script
```

Lo script usa solo la libreria standard (`dati/scaricamento.py`), e' ripetibile e
riprende uno scarico interrotto; un file con il nome giusto e il contenuto sbagliato e'
messo da parte, non sovrascritto. I file scaricati sono esclusi da git (`.gitignore`).

## Fonti e licenze

- SILVA 138 per dada2: record Zenodo 3986799, versione 2 (DOI 10.5281/zenodo.3986799),
  licenza CC BY 4.0 (`SILVA_LICENSE.txt`).
- Elenco dei taxa difettosi: mikemc/dada2-reference-databases, licenza MIT, fissato al
  commit indicato in `FONTI.tsv`.

## Il classificatore addestrato per IdTaxa

Con `tax.classifier: idtaxa` la pipeline usa, al posto del riferimento FASTA, un
classificatore addestrato nel formato di DECIPHER (RIFERIMENTO.md, 1.5):

```bash
python3 dati/riferimento/scarica_addestrato.py 138.2                 # SILVA SSU r138.2, circa 299 MB
python3 dati/riferimento/scarica_addestrato.py 2019                  # SILVA SSU r138, circa 152 MB
python3 dati/riferimento/scarica_addestrato.py 2019 --solo-verifica  # controlla e basta
```

| File | Classificatore | Quando | Byte | MD5 (`tax.idtaxa_trained_md5`) |
|---|---|---|---|---|
| `SILVA_SSU_r138_2_2024.RData` | SILVA SSU r138.2 | per un'analisi nuova | 298.719.756 | `4e272e39c2d71f5d3e7a31b00dbb1df4` |
| `SILVA_SSU_r138_2019.RData` | SILVA SSU r138 | per confrontarsi con tabelle elaborate da NASA GeneLab con le versioni A e B della sua pipeline (GL-DPPD-7104) | 151.849.504 | `cb983b6a5e8cdb46f8c88b5afae21f66` |

I due classificatori sono addestrati dagli autori di DECIPHER. Il sito di DECIPHER non
serve il file del 2019 e non pubblica checksum: lo script scarica le copie archiviate su
figshare dal gruppo di elaborazione dati di GeneLab (DOI 10.6084/m9.figshare.28550993 e
10.6084/m9.figshare.25799980, licenza CC BY 4.0) e ne verifica la dimensione e l'MD5
pubblicati da figshare. Che la copia del 2019 sia il file allora distribuito da DECIPHER
lo dichiara chi l'ha archiviata. I file scaricati sono esclusi da git.

Chi ha gia' i file nella collocazione precedente (`dati/osd734/riferimento/`) puo'
spostarli qui: lo script li riconosce dall'MD5 e non li scarica di nuovo.
