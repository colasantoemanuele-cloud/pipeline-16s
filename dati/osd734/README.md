# Dati di OSD-734

Questa cartella rende recuperabili da un solo punto i dati con cui la pipeline si
esegue sul dataset di riferimento, OSD-734 (NASA OSDR; ENA PRJEB73327): 960 campioni di
sequenziamento 16S single-end dalla Stazione Spaziale Internazionale. Nel repository
stanno solo i file piccoli e necessari; i file grandi si scaricano con gli script, che
verificano ogni file contro il checksum della sua fonte. Fonte, licenza e checksum di
ogni file sono in `FONTI.tsv`. Il riferimento tassonomico, comune ai dataset, sta in
`dati/riferimento/`, con il suo script, le sue fonti e le sue impronte.

## Che cosa c'e'

| Percorso | Nel repository | Che cos'e' |
|---|---|---|
| `metadati/` | si' | le tabelle ISA usate dalla pipeline: assay del 16S (`io.assay_table`) e campioni di studio (`io.study_table`), copiate senza modifiche dall'archivio ISA di OSDR |
| `letture_ena.tsv` | si' | le 960 corse: accession, nome del file locale, indirizzo, MD5 e dimensione dichiarati da ENA |
| `config_osd734.yaml` | si' | la configurazione per OSD-734, con i percorsi relativi alla radice del repository |
| `eco_osd734.yaml` | si' | la configurazione d'esempio delle analisi ecologiche di base sull'oggetto finale |
| `FONTI.tsv` | si' | fonte, licenza e controllo di ogni file; per il file del lotto anche l'identificativo dell'archivio Software Heritage |
| `scarica_letture.py` | si' | scarica le letture in `fastq/` |
| `ricostruisci_lotto.py` | si' | ricostruisce il file del lotto in `lotto/` dalla fonte pubblica |
| `confronta_risultati.py` | si' | confronta i risultati di un'esecuzione con i checksum attesi |
| `checksum_finali.sha256` | si' | i checksum attesi dei file consegnati di `12_final/`, nella forma di `sha256sum` |
| `checksum_artefatti.tsv` | si' | i checksum attesi di tutti i 1.033 artefatti dei manifesti di fase, con la fase che li produce |
| `fastq/`, `lotto/` | no | i file ottenuti dagli script, esclusi da `.gitignore` |
| `../riferimento/` | in parte | lo script `scarica_riferimento.py`, `FONTI.tsv` e `impronte.md5` del riferimento tassonomico; i file scaricati sono esclusi da `.gitignore` |
| `../scaricamento.py` | si' | lo scarico verificato e riprendibile comune agli script |

## Come recuperare i dati

Dalla radice del repository, con Python 3.11 o successivo e la sola libreria standard:

```bash
python3 dati/osd734/scarica_letture.py            # 960 file, 2,52 GB, da ENA
python3 dati/riferimento/scarica_riferimento.py   # SILVA 138 per dada2, circa 138 MB, da Zenodo
python3 dati/osd734/ricostruisci_lotto.py         # il file del lotto, dalla fonte degli autori
```

Lo scarico delle letture e' il passo lungo: i file si scaricano quattro alla volta
(`--connessioni` cambia il numero), e la durata dipende dal collegamento con ENA: puo'
richiedere alcune ore. Lo si puo' interrompere e rilanciare: riprende dal punto a cui
era arrivato.

Gli script sono ripetibili: a ogni avvio controllano tutti i file e scaricano solo cio'
che manca o non corrisponde al checksum; uno scarico interrotto riprende dal punto a cui
era arrivato. Un file presente con il nome giusto ma il contenuto sbagliato non viene
sovrascritto: e' messo da parte con il suffisso `.md5_errato`. I due script di scarico
(`scarica_letture.py` e `dati/riferimento/scarica_riferimento.py`) con `--solo-verifica`
controllano senza scaricare ne' spostare nulla, e con `--cartella` lavorano su un'altra
cartella, per esempio una che contiene gia' i file; l'esito e' 0 se tutti i file sono
presenti e integri. `ricostruisci_lotto.py` non ha queste opzioni: ricostruisce il
file dalla fonte (`--fonte` per una copia locale, `--uscita` per la destinazione), e
rilanciarlo riscrive lo stesso file. Prima di consegnarlo confronta il file ricostruito
con l'MD5 dichiarato in `FONTI.tsv`: se non coincide si ferma e non scrive nulla.

## Come eseguire la pipeline

I calcoli richiedono R 4.5.2 con Bioconductor e dada2 1.36.0.1 (la 1.36.0 con la
correzione dei pareggi di `assignTaxonomy`), presenti solo nell'immagine del container:
la pipeline non si esegue sulla macchina dell'utente ma nel container, con il
repository montato.

**1. Procurarsi l'immagine.** La via principale e' scaricare l'immagine pubblicata,
identificata per digest:

```bash
docker pull ghcr.io/colasantoemanuele-cloud/amplicon16s@sha256:930dda008581b24cc3a3fea718c0fed3690b011958bf1716c643f92056970154
```

Il digest (`sha256:930dda008581b24cc3a3fea718c0fed3690b011958bf1716c643f92056970154`) identifica i byte dell'immagine: il nome con cui e' pubblicata potra'
cambiare, il digest no, e chi la scarica con quel digest ottiene la stessa immagine. Lo
scarico e' di 1,86 GB compressi; su disco l'immagine occupa 7,6 GB. Nei comandi che
seguono `<immagine>` e' quel riferimento completo.

In alternativa l'immagine si costruisce dal Dockerfile, dalla radice del repository
(serve la rete: le immagini di partenza, i pacchetti R dall'istantanea datata di CRAN e
i sorgenti di dada2):

```bash
docker build -f container/Dockerfile -t amplicon16s:locale .
```

e `<immagine>` e' allora `amplicon16s:locale`. **Due costruzioni non sono identiche
byte per byte**: date dei file, ordine di compressione e metadati cambiano, quindi
un'immagine costruita in locale ha un altro identificativo e non e' l'immagine
pubblicata. Contiene pero' gli stessi pacchetti alle stesse versioni, verificati contro
`renv.lock` durante la costruzione.

**2. Adattare la configurazione** `config_osd734.yaml`, prima dell'esecuzione. I suoi
percorsi sono relativi alla radice del repository, da cui il comando va lanciato: la
riga di comando li rende assoluti al caricamento, e la configurazione registrata in
`00_config/` riporta i percorsi dei file letti (nel container, sotto `/app`).

- `run.threads` non e' impostato e vale i processori utilizzabili dal processo: va
  indicato solo per usarne meno. Non cambia i risultati: e' fra i parametri senza
  effetto, fuori dall'impronta delle fasi.
- I parametri che dipendono dai dati o dallo studio (colonne e etichette dei metadati,
  colonne del file del lotto, primer, motivo conservato, troncamento, riferimento,
  scelte sui controlli, soglie di qualita') sono dichiarati tutti: la pipeline non ha
  per essi alcun valore predefinito, e una configurazione che ne omette uno e' respinta
  prima di qualunque calcolo (`E-G15-10`). La configurazione dichiara anche i parametri
  che un predefinito lo hanno, cosi' che i risultati non dipendano da esso.
- `run.strict_provenance` e' vero: e' la configurazione congelata, e l'esecuzione parte
  solo se il codice e l'ambiente sono quelli dichiarati. Prima di qualunque fase la
  pipeline verifica che il repository sia un clone git leggibile, che `src/` e `R/`
  non abbiano modifiche non committate ne' file non tracciati, e che la versione di R
  e quelle dei pacchetti installati siano quelle di `renv.lock`; altrimenti si ferma
  con `E-PROV-01`, `E-PROV-02` o `E-PROV-03`, e il messaggio dice che cosa manca. Va
  quindi eseguita da un clone (non da un archivio scaricato senza `.git`), con il
  comando qui sotto e nell'immagine pubblicata. Modificare questo file di
  configurazione, o scaricare i dati in `dati/`, non conta come modifica del codice.
  Chi vuole provare una modifica al codice porta il parametro a `false`: i risultati
  non saranno certificati rispetto al codice, e le fasi gia' calcolate con la regola
  vengono rifatte.
- Cinque valori decidono quali campioni e quali varianti restano, e cambiarli cambia i
  risultati. Il criterio di ciascuno, in una riga:
  - `katharoseq.target_sensitivity: 0.90`: nei campioni piu' poveri si accetta al piu'
    un decimo di letture estranee al campione;
  - `decontam.threshold: 0.5`: bassa biomassa, con controlli negativi di profondita'
    paragonabile ai biologici; e' contaminante cio' che e' piu' prevalente nei negativi;
  - `decontam.mode: aggregate`: le piastre hanno da 5 a 8 controlli negativi, e per
    piastra, con questa soglia, la decontaminazione toglierebbe una parte troppo grande
    delle letture dei biologici;
  - `decontam.batch_combine: minimum`: il predefinito di `decontam`, dichiarato;
  - `prev.min_fraction: 0.01`: i campioni vengono da oltre cento posizioni distinte, e
    una variante propria di poche posizioni e' segnale, non rumore.
- `run.container` dichiara gia' l'immagine pubblicata, con il suo digest di registro, e
  con quella non va toccato. Il valore non sceglie l'immagine da eseguire, che e'
  quella passata a `docker run`: dichiara con quale immagine i risultati sono stati
  prodotti. E' registrato nella configurazione salvata in `00_config/` e nella
  provenienza del manifesto di ogni fase; una ripresa con un'immagine dichiarata
  diversa lo segnala con un avviso di provenienza. Solo chi ha costruito l'immagine dal
  Dockerfile lo sostituisce con l'identificativo locale, che vale soltanto sulla sua
  macchina e si ottiene con

  ```bash
  echo "amplicon16s@$(docker image inspect amplicon16s:locale --format '{{.Id}}')"
  ```

  (un tag non e' accettato, perche' puo' essere riassegnato a un'immagine diversa).

**3. Eseguire**, dalla radice del repository:

```bash
docker run --rm \
  --memory=24g \
  --memory-swap=24g \
  -u "$(id -u):$(id -g)" \
  -v "$(pwd)":/app \
  -w /app \
  <immagine> \
  python3 scripts/esegui.py run --config dati/osd734/config_osd734.yaml
```

- `-v "$(pwd)":/app -w /app`: il repository montato e' la cartella di lavoro, quindi i
  percorsi relativi della configurazione trovano dati e configurazione, e i risultati
  vanno in `output/osd734/` del repository;
- `python3 scripts/esegui.py`: esegue il codice Python e gli script R del repository
  montato, senza variabili d'ambiente. Il comando `amplicon16s` dell'immagine
  eseguirebbe invece il codice copiato quando l'immagine e' stata costruita, e la
  regola rigorosa sulla provenienza, attiva in questa configurazione, lo rifiuterebbe
  dicendolo (`E-PROV-01`);
- `-u`: i file prodotti appartengono all'utente e non a root, con i permessi della sua
  umask;
- `--memory=24g --memory-swap=24g`: un esaurimento della memoria ferma il container e
  non la macchina (vedi i requisiti sotto).

Un'esecuzione interrotta si riprende con lo stesso comando e `resume` al posto di
`run` (un arresto stampa il comando di ripresa); `validate` esegue solo i controlli di
avvio e la fase S0, in meno di un minuto.

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

**A che cosa si riferiscono i checksum.** All'esecuzione di riferimento, con la
configurazione e l'immagine indicate: `config_osd734.yaml` cosi' com'e' (filogenesi
disattivata, regola rigorosa sulla provenienza attiva), il codice del repository e
l'immagine pubblicata,
`ghcr.io/colasantoemanuele-cloud/amplicon16s@sha256:930dda008581b24cc3a3fea718c0fed3690b011958bf1716c643f92056970154`.
`run.threads` non incide sui risultati. Un'immagine costruita da `container/Dockerfile`
(immagine di partenza ancorata per digest, pacchetti R dell'istantanea datata e di
`renv.lock`, dada2 1.36.0.1) non e' la stessa immagine: ha gli stessi pacchetti alle
stesse versioni, ma le impronte si riferiscono a quella pubblicata. Fuori dal container,
con un'altra versione di R o dei pacchetti, i risultati possono differire senza che
l'esecuzione sia sbagliata: il confronto prova l'identita', non la correttezza.

L'oggetto finale di riferimento ha 451 campioni e 1.737 varianti; gli artefatti dei
manifesti di fase sono 1.033, e i file consegnati di `12_final/` sei.

### Se i risultati non coincidono

1. Guardare la prima fase indicata dallo script: le fasi successive dipendono da
   quella, quindi le loro differenze sono conseguenze, non cause.
2. Se la prima fase e' S0, S1 o S2, la differenza e' nei dati di ingresso: rilanciare i
   due script di scarico con `--solo-verifica`, che devono terminare con esito 0,
   ricostruire di nuovo il file del lotto con `ricostruisci_lotto.py`, e
   confrontare la configurazione registrata in `output/osd734/00_config/resolved.yaml`
   con `config_osd734.yaml` (devono differire solo i percorsi).
3. Se e' una fase successiva, la differenza e' nell'ambiente di calcolo: controllare
   di aver eseguito nell'immagine pubblicata, scaricata con il suo digest, con
   `python3 scripts/esegui.py` (il comando `amplicon16s` dell'immagine esegue il codice
   copiato alla costruzione), e che il sottocomando `report` non riporti avvisi di
   provenienza ne' artefatti non integri. Il report indica la versione di dada2
   caricata da R in S8, accanto a quella dichiarata dal file di blocco: deve essere
   1.36.0.1, con la correzione dei pareggi presente, senza la quale la tassonomia (S8)
   non e' ripetibile e l'esecuzione lo dichiara con `E-S8-03`.
4. Una differenza che resta va segnalata con l'uscita dello script, il report e la
   configurazione registrata: e' un difetto di riproducibilita' della pipeline, non un
   errore di chi la esegue.

## Requisiti

Misure indicative, con dodici processori e l'immagine pubblicata.

| Che cosa | Misura |
|---|---|
| Letture FASTQ | 2,52 GB, 960 file (le dimensioni dichiarate da ENA in `letture_ena.tsv`) |
| Riferimento SILVA 138 | 138 MB |
| Immagine pubblicata | 1,86 GB da scaricare, 7,6 GB su disco |
| Uscite dell'esecuzione completa | 2,3 GB, quasi tutti di letture filtrate (`03_filtered/`); 1.033 artefatti nei manifesti di fase |
| Durata dell'esecuzione completa | circa 50 minuti per la catena S0-S14; l'inferenza delle varianti (S4) ne occupa circa meta', il filtro (S2) da 8 a 10 |
| Memoria di picco | 13 GB, per l'intero container, cache dei file compresa |

Le letture filtrate si conservano con `run.keep_filtered_fastq: true`, il predefinito;
con `false` si rimuovono quando tutte le fasi sono concluse, quindi lo spazio serve
comunque durante l'esecuzione. I controlli di avvio richiedono, sul volume di
`io.out_root`, almeno lo spazio occupato dalle letture. Il limite di 24 GB del comando di
esecuzione lascia margine sulla memoria di picco; con meno di 16 GB di memoria, e con
meno processori, l'esecuzione non e' stata misurata.

## Le analisi ecologiche d'esempio

Sull'oggetto finale di un'esecuzione, con la configurazione d'esempio (il modulo della
stazione come variabile biologica, corsa e piastra come variabili tecniche):

```bash
docker run --rm -u "$(id -u):$(id -g)" -v "$(pwd)":/app -w /app <immagine> \
  python3 scripts/esegui.py eco run --object output/osd734/12_final/ps_final.rds \
    --config dati/osd734/eco_osd734.yaml --out output/eco_osd734
```

La cartella di uscita deve essere assente o vuota. Un'esecuzione dura circa un minuto e
mezzo; due esecuzioni danno gli stessi byte in ogni file.

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

La fonte, al commit fissato, e' archiviata su Software Heritage: l'identificativo
(SWHID) del file e' in `FONTI.tsv`. Se il repository degli autori non fosse piu'
raggiungibile, il file si ottiene dall'archivio, con la stessa impronta, e si passa
allo script con `--fonte`.

## Licenze

- Metadati ISA: NASA OSDR, distribuiti secondo la Scientific Information Policy del NASA
  Science Mission Directorate (SPD-41a), che rende i dati pubblici senza restrizioni e
  raccomanda la licenza CC0. Citazione: OSD-734, NASA Open Science Data Repository.
- Letture ed elenco delle corse: ENA, secondo la politica INSDC, senza restrizioni d'uso.
- SILVA 138 per dada2 ed elenco dei taxa difettosi del riferimento: fonti e licenze
  in `dati/riferimento/` (CC BY 4.0 e MIT).
- Metadati Qiita degli autori: nessuna licenza dichiarata; usati solo come fonte della
  ricostruzione, non ridistribuiti.
