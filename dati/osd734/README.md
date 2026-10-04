# Dati di OSD-734

Questa cartella rende recuperabili da un solo punto i dati con cui la pipeline si
esegue sul dataset di riferimento, OSD-734 (NASA OSDR; ENA PRJEB73327): 960 campioni di
sequenziamento 16S single-end dalla Stazione Spaziale Internazionale. Nel repository
stanno solo i file piccoli e necessari; i file grandi si scaricano con gli script, che
verificano ogni file contro il checksum della sua fonte. Fonte, licenza e checksum di
ogni file sono in `FONTI.tsv`.

## Che cosa c'e'

| Percorso | Nel repository | Che cos'e' |
|---|---|---|
| `metadati/` | si' | le tabelle ISA usate dalla pipeline: assay del 16S (`io.assay_table`) e campioni di studio (`io.study_table`), copiate senza modifiche dall'archivio ISA di OSDR |
| `letture_ena.tsv` | si' | le 960 corse: accession, nome del file locale, indirizzo, MD5 e dimensione dichiarati da ENA |
| `config_osd734.yaml` | si' | la configurazione per OSD-734, con i percorsi relativi alla radice del repository |
| `scarica_letture.py` | si' | scarica le letture in `fastq/` |
| `scarica_riferimento.py` | si' | scarica il riferimento tassonomico in `riferimento/` |
| `ricostruisci_lotto.py` | si' | ricostruisce il file del lotto in `lotto/` dalla fonte pubblica |
| `scaricamento.py` | si' | lo scarico verificato e riprendibile comune ai due script |
| `fastq/`, `riferimento/`, `lotto/` | no | i file ottenuti dagli script, esclusi da `.gitignore` |

## Come recuperare i dati

Dalla radice del repository, con Python 3.11 e la sola libreria standard:

```bash
python3 dati/osd734/scarica_letture.py        # 960 file, circa 2,4 GB, da ENA
python3 dati/osd734/scarica_riferimento.py    # SILVA 138 per dada2, circa 138 MB, da Zenodo
python3 dati/osd734/ricostruisci_lotto.py     # il file del lotto, dalla fonte degli autori
```

Gli script sono ripetibili: a ogni avvio controllano tutti i file e scaricano solo cio'
che manca o non corrisponde al checksum; uno scarico interrotto riprende dal punto a cui
era arrivato. Un file presente con il nome giusto ma il contenuto sbagliato non viene
sovrascritto: e' messo da parte con il suffisso `.md5_errato`. Con `--solo-verifica`
controllano senza scaricare ne' spostare nulla; con `--cartella` lavorano su un'altra
cartella, per esempio una che contiene gia' i file. L'esito e' 0 se tutti i file sono
presenti e integri.

Poi la pipeline, sempre dalla radice del repository, dopo aver indicato in
`config_osd734.yaml` il digest dell'immagine (`run.container`):

```bash
amplicon16s run --config dati/osd734/config_osd734.yaml
```

## Il file del lotto

`lotto/plate_well_map_960.tsv` (`io.batch_table`) associa a ogni campione la piastra
di estrazione, il pozzetto, la corsa di sequenziamento, il modulo della stazione e, per
i controlli positivi, le cellule seminate. Senza, la pipeline procede in modalita' a
corsa singola e senza lotto: un solo modello d'errore e nessuna soglia di profondita'
per piastra. I risultati di riferimento lo richiedono.

Deriva dai metadati Qiita (studio 14542) che gli autori dello studio hanno pubblicato
nel loro repository (`RodolfoSalido/3DMM`). Quel repository non dichiara una licenza,
quindi il file non e' ridistribuito qui: `ricostruisci_lotto.py` lo ricostruisce dalla
fonte, fissata a un commit e verificata con l'impronta git del contenuto, e ne documenta
la provenienza colonna per colonna. L'unica colonna che non viene dalla fonte e'
l'accession dell'esperimento, ricavato dalla tabella di assay: e' la chiave del join,
perche' il nome del campione si ripete fra le repliche.

## Licenze

- Metadati ISA: NASA OSDR, distribuiti secondo la Scientific Information Policy del NASA
  Science Mission Directorate (SPD-41a), che rende i dati pubblici senza restrizioni e
  raccomanda la licenza CC0. Citazione: OSD-734, NASA Open Science Data Repository.
- Letture ed elenco delle corse: ENA, secondo la politica INSDC, senza restrizioni d'uso.
- SILVA 138 per dada2: CC BY 4.0 (record Zenodo 3986799, `SILVA_LICENSE.txt`).
- Elenco dei taxa difettosi del riferimento: MIT (mikemc/dada2-reference-databases).
- Metadati Qiita degli autori: nessuna licenza dichiarata; usati solo come fonte della
  ricostruzione, non ridistribuiti.
