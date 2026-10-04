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
python3 dati/osd734/scarica_letture.py        # 960 file, 2,52 GB, da ENA
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

## Come eseguire la pipeline

I calcoli richiedono R 4.5.2 con Bioconductor e dada2 1.36.0.1 (la 1.36.0 con la
correzione dei pareggi di `assignTaxonomy`), presenti solo nell'immagine del container:
la pipeline non si esegue sulla macchina dell'utente ma nel container, con il
repository montato. Nei comandi `<immagine>` e' il nome che si sceglie per l'immagine
costruita, per esempio `amplicon16s:locale`.

**1. Costruire l'immagine**, dalla radice del repository (serve la rete: le immagini di
partenza, i pacchetti R dall'istantanea datata di CRAN e i sorgenti di dada2):

```bash
docker build -f container/Dockerfile -t <immagine> .
```

**2. Adattare la configurazione** `config_osd734.yaml`, prima dell'esecuzione:

- `run.threads`, per primo. La configurazione non lo imposta e vale 16: se la macchina
  ha meno processori utilizzabili, i controlli di avvio fermano l'esecuzione con
  `E-S0-14` prima di qualunque fase. Va indicato il numero di processori della macchina
  (`nproc`) o meno. Non cambia i risultati: e' fra i parametri senza effetto, fuori
  dall'impronta delle fasi.
- `run.container`, l'immagine con cui si esegue, ancorata per digest
  (`nome@sha256:<64 cifre esadecimali>`; un tag non e' accettato, perche' puo' essere
  riassegnato a un'immagine diversa). Un'immagine costruita in locale non ha un digest
  di registro: si indica il suo identificativo di contenuto, che si ottiene con

  ```bash
  echo "amplicon16s@$(docker image inspect <immagine> --format '{{.Id}}')"
  ```

  Il valore non sceglie l'immagine da eseguire, che e' quella passata a `docker run`:
  dichiara con quale immagine i risultati sono stati prodotti. E' registrato nella
  configurazione salvata in `00_config/` e nella provenienza del manifesto di ogni
  fase; una ripresa con un'immagine dichiarata diversa lo segnala con un avviso di
  provenienza. L'identificativo locale vale solo sulla macchina che ha costruito
  l'immagine: per un'immagine pubblicata in un registro si indica il digest con cui la
  si recupera (vedi "Identificare l'immagine per digest" nel README principale).

**3. Eseguire**, dalla radice del repository:

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
  amplicon16s run --config dati/osd734/config_osd734.yaml
```

- `-v "$(pwd)":/app -w /app`: il repository montato e' la cartella di lavoro, quindi i
  percorsi relativi della configurazione trovano dati e configurazione, e i risultati
  vanno in `output/osd734/` del repository;
- `PYTHONPATH` e `AMPLICON16S_R_DIR`: si eseguono il codice Python e gli script R del
  repository montato; senza la seconda variabile l'immagine eseguirebbe gli script R
  copiati quando e' stata costruita;
- `-u` e `HOME`: i file prodotti appartengono all'utente e non a root;
- `--memory=24g --memory-swap=24g`: un esaurimento della memoria ferma il container e
  non la macchina (vedi i requisiti sotto).

Un'esecuzione interrotta si riprende con lo stesso comando e `amplicon16s resume` al
posto di `amplicon16s run`; `amplicon16s validate` esegue solo i controlli di avvio e la
fase S0, in meno di un minuto.

Questo comando, con `amplicon16s validate`, e' stato provato sui file gia' presenti
sulla macchina di sviluppo: con la configurazione cosi' com'e' si ferma con `E-S0-14`
(12 processori utilizzabili), con `run.threads: 12` conclude S0 e i suoi artefatti
sono identici byte per byte a quelli dell'esecuzione di riferimento. L'esecuzione
completa da dati appena scaricati non e' ancora stata provata.

## Requisiti misurati

Misure sulla macchina di sviluppo (Intel Xeon W-10855M, 12 processori logici, 30 GB di
memoria), nel container della settimana 21, con `run.threads: 12`.

| Che cosa | Misura | Fonte |
|---|---|---|
| Letture FASTQ | 2,52 GB (2.515.101.383 byte) | somma delle dimensioni dichiarate da ENA in `letture_ena.tsv` |
| Riferimento SILVA 138 | 138 MB (137.973.851 byte) | `scarica_riferimento.py` |
| Immagine del container | 7,59 GB su disco, strati dell'immagine di partenza di Bioconductor compresi | `docker image ls` (colonna DISK USAGE) e `docker system df -v` (colonna SIZE), immagine della settimana 21 |
| Uscite dell'esecuzione completa | 2,27 GB, di cui 2,24 GB di letture filtrate (`03_filtered/`) e 32 MB di tutto il resto | esecuzione di riferimento S0-S14, settimana 21 |
| Durata dell'esecuzione completa | 55,5 e 64 minuti, in due esecuzioni dell'intera catena S0-S14; la fase S4 (denoising) ne occupa 26 e 32 | esecuzioni di riferimento, settimana 21 |
| Memoria di picco della pipeline | non misurata | |

Lo spazio dell'immagine e' misurato con Docker 29.1.3, che usa l'archivio di immagini
di containerd. Dei 7,59 GB, 5,67 sono strati condivisi con altre immagini della
macchina (quella di Bioconductor da cui l'immagine e' costruita, 6,55 GB da sola, e le
versioni precedenti della pipeline); la parte propria dell'immagine e' 1,92 GB.
`docker image inspect --format '{{.Size}}'`, il comando usato in precedenza, riporta
invece 1,86 GB: con l'archivio di containerd e' la dimensione del contenuto, cioe'
degli strati compressi che si scaricano da un registro (colonna CONTENT SIZE di
`docker image ls`), non lo spazio occupato su disco. La cache di
costruzione non e' compresa, e su una costruzione isolata non e' stata misurata; non e'
stato misurato nemmeno lo spazio complessivo su una macchina che parte da zero.

Le letture filtrate si conservano con `run.keep_filtered_fastq: true`, il predefinito;
con `false` si rimuovono quando tutte le fasi sono concluse, quindi lo spazio serve
comunque durante l'esecuzione. I controlli di avvio richiedono, sul volume di
`io.out_root`, almeno lo spazio occupato dalle letture. La durata con meno processori
non e' stata misurata.

La memoria di picco della sola pipeline non e' stata misurata. L'unica misura
disponibile e' quella della suite di test con i dati reali, che sul dataset completo
esegue la catena S0-S14 e altri test, uno dopo l'altro: 16,3 GiB di picco del container
(`memory.peak`), con il limite a 24 GB e senza esaurimenti.

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
