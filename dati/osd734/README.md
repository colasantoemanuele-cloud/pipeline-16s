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
| `confronta_risultati.py` | si' | confronta i risultati di un'esecuzione con i checksum attesi |
| `checksum_finali.sha256` | si' | i checksum attesi dei file consegnati di `12_final/`, nella forma di `sha256sum` |
| `checksum_artefatti.tsv` | si' | i checksum attesi di tutti i 1.031 artefatti dei manifesti di fase, con la fase che li produce |
| `fastq/`, `riferimento/`, `lotto/` | no | i file ottenuti dagli script, esclusi da `.gitignore` |

## Come recuperare i dati

Dalla radice del repository, con Python 3.11 o successivo e la sola libreria standard:

```bash
python3 dati/osd734/scarica_letture.py        # 960 file, 2,52 GB, da ENA
python3 dati/osd734/scarica_riferimento.py    # SILVA 138 per dada2, circa 138 MB, da Zenodo
python3 dati/osd734/ricostruisci_lotto.py     # il file del lotto, dalla fonte degli autori
```

Lo scarico delle letture e' il passo lungo: i file si scaricano uno alla volta, e la
durata dipende dal collegamento con ENA. Nella riproduzione del 5 ottobre 2026 ha
impiegato 2 ore e 48 minuti, senza interruzioni; il riferimento 44 secondi, il file del
lotto meno di un secondo. Lo si puo' interrompere e rilanciare: riprende dal punto a cui
era arrivato.

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

**2. Adattare la configurazione** `config_osd734.yaml`, prima dell'esecuzione. I suoi
percorsi sono relativi alla radice del repository, da cui il comando va lanciato: la
riga di comando li rende assoluti al caricamento, e la configurazione registrata in
`00_config/` riporta i percorsi dei file letti (nel container, sotto `/app`). Vanno
adattati due valori:

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

**La procedura e' stata eseguita da zero** il 5 ottobre 2026: repository clonato da
GitHub (commit `59ef39d`), immagine costruita dal suo `container/Dockerfile` senza
cache (4 minuti e 37 secondi), dati scaricati con i tre script in una cartella vuota,
catena S0-S14 con il comando qui sopra e `run.threads: 12`. Tutti i 1.031 artefatti
dei manifesti di fase sono risultati identici, byte per byte, a quelli dell'esecuzione
di riferimento precedente, prodotta su un'altra copia dei dati e con un'immagine
costruita in un altro momento. Quell'esecuzione e' ora il riferimento, e i suoi
checksum sono quelli pubblicati qui.

La prova ha fatto emergere un difetto, corretto: al commit `59ef39d` i percorsi
relativi di `config_osd734.yaml` superavano S0 ma fermavano S1 con `E-R-02`, perche' il
processo R di una fase parte nella cartella della fase e non trovava i propri
ingressi. La riproduzione e' stata completata indicando nella configurazione i
percorsi assoluti del container (`/app/dati/osd734/...`, `/app/output/osd734`); dalla
versione successiva la riga di comando rende assoluti i percorsi relativi, che danno
la stessa configurazione registrata, e il caso e' coperto da un test.

## Come verificare i risultati

Dalla radice del repository, a esecuzione conclusa:

```bash
python3 dati/osd734/confronta_risultati.py
```

Lo script ricalcola l'impronta SHA-256 di ogni artefatto elencato nei manifesti di
fase di `output/osd734/` (un'altra cartella si indica con `--uscita`) e la confronta
con `checksum_artefatti.tsv`. Riporta l'esito fase per fase, nell'ordine di esecuzione,
e termina con `RISULTATI IDENTICI` ed esito 0, oppure con `RISULTATI DIVERSI`, la prima
fase in cui compare una differenza ed esito 1. Impiega pochi secondi. Chi vuole
verificare solo i file consegnati, senza Python:

```bash
(cd output/osd734/12_final && sha256sum -c ../../../dati/osd734/checksum_finali.sha256)
```

Non si confrontano i manifesti, la configurazione registrata in `00_config/`, i log e
il report: portano date, identificativi dell'esecuzione e percorsi, e differiscono per
costruzione fra due esecuzioni.

**A quale ambiente si riferiscono i checksum.** All'immagine costruita da
`container/Dockerfile` (immagine di partenza ancorata per digest, pacchetti R
dell'istantanea datata e di `renv.lock`, dada2 1.36.0.1), con il codice del repository
e la configurazione `config_osd734.yaml`, filogenesi disattivata; `run.threads` non
incide sui risultati. Un'immagine pubblicata in un registro, con il suo digest, non
c'e' ancora: finche' non c'e', chi riproduce costruisce l'immagine dal Dockerfile, come
nella prova descritta sopra, in cui due immagini costruite in momenti diversi hanno
dato gli stessi byte. Fuori dal container, con un'altra versione di R o dei pacchetti,
i risultati possono differire senza che l'esecuzione sia sbagliata: il confronto prova
l'identita', non la correttezza.

### Se i risultati non coincidono

1. Guardare la prima fase indicata dallo script: le fasi successive dipendono da
   quella, quindi le loro differenze sono conseguenze, non cause.
2. Se la prima fase e' S0, S1 o S2, la differenza e' nei dati di ingresso: rilanciare i
   tre script di scarico con `--solo-verifica`, che devono terminare con esito 0, e
   confrontare la configurazione registrata in `output/osd734/00_config/resolved.yaml`
   con `config_osd734.yaml` (devono differire solo i percorsi, `run.threads` e
   `run.container`).
3. Se e' una fase successiva, la differenza e' nell'ambiente di calcolo: controllare
   di aver eseguito nell'immagine costruita dal Dockerfile del repository, con
   `AMPLICON16S_R_DIR=/app/R` (senza, l'immagine esegue gli script R copiati alla
   costruzione), e che `amplicon16s report --config dati/osd734/config_osd734.yaml`
   non riporti avvisi di provenienza ne' artefatti non integri. Il report indica anche
   la versione di dada2 dichiarata dal file di blocco: deve essere 1.36.0.1, con la
   correzione dei pareggi, senza la quale la tassonomia (S8) non e' riproducibile.
4. Una differenza che resta va segnalata con l'uscita dello script, il report e la
   configurazione registrata: e' un difetto di riproducibilita' della pipeline, non un
   errore di chi la esegue.

## Requisiti misurati

Misure sulla macchina di sviluppo (Intel Xeon W-10855M, 12 processori logici, 30 GB di
memoria), con `run.threads: 12`. Dove non e' detto altrimenti vengono dalla
riproduzione da zero del 5 ottobre 2026 (codice del commit `59ef39d`, immagine costruita
dal suo Dockerfile).

| Che cosa | Misura | Fonte |
|---|---|---|
| Letture FASTQ | 2,52 GB (2.515.101.383 byte), 960 file | somma delle dimensioni dichiarate da ENA in `letture_ena.tsv`, uguale ai byte scaricati |
| Riferimento SILVA 138 | 138 MB (137.973.851 byte) | `scarica_riferimento.py` |
| Durata dello scarico | 2 ore e 48 minuti le letture, 44 secondi il riferimento, meno di un secondo il lotto | i tre script, in una cartella vuota; dipende dal collegamento con ENA |
| Costruzione dell'immagine | 4 minuti e 37 secondi senza cache | `docker build --no-cache`; dipende dalla rete |
| Immagine del container | 7,59 GB su disco, strati dell'immagine di partenza di Bioconductor compresi | `docker image ls` (colonna DISK USAGE) |
| Uscite dell'esecuzione completa | 2,27 GB, di cui 2,24 GB di letture filtrate (`03_filtered/`) e 32 MB di tutto il resto | `output/osd734/` a catena conclusa |
| Durata dell'esecuzione completa | 53 minuti e 20 secondi per la catena S0-S14; la fase S4 (denoising) ne occupa 26, S2 9,5, S8 6, S1 5 | il comando di esecuzione qui sopra; nelle esecuzioni precedenti 55,5 e 64 minuti |
| Memoria di picco della pipeline | 12,3 GB (11,5 GiB; 12.304.121.856 byte) | `memory.peak` del cgroup del container, letto al termine della catena, con il limite a 24 GB |

Lo spazio dell'immagine e' misurato con Docker 29.1.3, che usa l'archivio di immagini
di containerd. Dei 7,59 GB, 5,67 sono strati condivisi con altre immagini della
macchina (quella di Bioconductor da cui l'immagine e' costruita, 6,55 GB da sola, e le
versioni precedenti della pipeline); la parte propria dell'immagine e' 1,92 GB (misura
della settimana 21). `docker image inspect --format '{{.Size}}'` riporta
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

La memoria di picco e' quella dell'intero container durante la sola catena, cache dei
file compresa: 12,3 GB su una macchina da 30 GB. Il limite di 24 GB del comando lascia
quindi margine; con meno di 16 GB di memoria l'esecuzione non e' stata provata. La
suite di test con i dati reali, che esegue la catena e altri test uno dopo l'altro,
arriva a 16,9 GB.

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
