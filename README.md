# amplicon16s

Pipeline per l'analisi di dati 16S rRNA da sequenziamento Illumina single-end.

## Che cosa fa

Dalle letture grezze di un esperimento di sequenziamento dell'amplicone 16S e dai suoi
metadati produce un oggetto `phyloseq` pronto per le analisi ecologiche: tabella dei
conteggi delle varianti esatte di sequenza (ASV), tassonomia, sequenze e metadati dei
campioni, con accanto gli stessi dati in formato testo.

Le fasi, nell'ordine:

| Fase | Che cosa fa |
|---|---|
| S0 | valida configurazione, metadati e file delle letture, prima di ogni calcolo |
| S1 | profilo di qualità e di lunghezza di tutte le letture |
| S2 | filtro e troncamento delle letture (`dada2::filterAndTrim`) |
| S3 | modello d'errore, uno per corsa di sequenziamento |
| S4-S7 | inferenza delle varianti, tabella delle sequenze, rimozione delle chimere, selezione per lunghezza |
| S8 | assegnazione tassonomica (`dada2::assignTaxonomy`) |
| S10 | oggetto integrato con tutti i campioni |
| S11 | curva di calibrazione dai controlli positivi in diluizione (KatharoSeq) e soglia di profondità |
| S12 | riconoscimento dei contaminanti dai controlli negativi (`decontam`) |
| S13 | filtri finali: profondità, tassonomia, prevalenza, letture minime |
| S9 | albero filogenetico delle varianti finali (facoltativa) |
| S14 | oggetto finale, export e checksum |

Python orchestra l'esecuzione; R esegue i calcoli (dada2, phyloseq, DECIPHER, phangorn,
decontam) come processi separati. Ogni fase scrive i propri artefatti con un manifesto
che ne riporta dimensione e SHA-256: un'esecuzione interrotta si riprende dal punto in
cui si è fermata, e a parità di configurazione, codice e immagine due esecuzioni danno
gli stessi byte.

## Dati supportati

- Letture Illumina **single-end** dell'amplicone 16S, un file FASTQ per campione,
  compresso con gzip o no. Di una corsa appaiata si usano le sole letture forward; un
  file che contiene le due letture di ogni coppia viene respinto.
- Letture **senza primer in testa**, oppure con il primer a posizione fissa, che si
  toglie con `filter.trimLeft`.
- Metadati in tabelle di testo UTF-8 separate da tabulazioni: una tabella di assay (un
  campione per riga, con il file delle letture) e, facoltativa, una tabella di studio
  con il tipo di ogni campione. Il formato ISA-Tab degli archivi OSDR è letto così com'è.
- Facoltativi: un file del lotto (piastra di estrazione e corsa di sequenziamento per
  campione), controlli negativi, controlli positivi in diluizione seriale. Senza di
  essi la pipeline procede e dichiara che cosa non ha potuto fare.
- Un riferimento tassonomico nel formato dei training set di dada2.

Nel repository ci sono due dataset su cui la pipeline è eseguita, con dati,
configurazione e istruzioni: `dati/osd734/` (960 campioni della Stazione Spaziale
Internazionale, con controlli e file del lotto) e `dati/osd276/` (15 campioni, senza
controlli).

## Requisiti

- **Docker**, e l'immagine pubblicata (1,86 GB da scaricare, 7,6 GB su disco). È
  l'unica via supportata: l'immagine porta Python 3.11, R 4.5.2, Bioconductor 3.21 e i
  pacchetti alle versioni di `renv.lock`, fra cui dada2 1.36.0.1, cioè dada2 1.36.0 con
  una correzione che rende ripetibile la tassonomia nei pareggi fra generi
  (`container/dada2/`).
- Un clone di questo repository: il codice si esegue dal clone montato nel container.
- Memoria e spazio in proporzione ai dati. Sul dataset di 960 campioni: 13 GB di
  memoria di picco, 2,3 GB di uscite, 50 minuti con 12 processori.

## Esecuzione

L'immagine si scarica per digest: il digest identifica i byte dell'immagine, un nome no.

```bash
docker pull ghcr.io/colasantoemanuele-cloud/amplicon16s@sha256:930dda008581b24cc3a3fea718c0fed3690b011958bf1716c643f92056970154
```

I comandi si lanciano dalla radice del clone, montato nel container. Nei comandi che
seguono `<immagine>` è il riferimento qui sopra.

```bash
docker run --rm -u "$(id -u):$(id -g)" -v "$(pwd)":/app -w /app <immagine> \
  python3 scripts/esegui.py validate --config <configurazione.yaml>
```

| Sottocomando | Che cosa fa |
|---|---|
| `validate` | esegue la sola validazione (S0) |
| `run` | esegue dall'inizio, in una cartella di uscita nuova |
| `resume` | riprende un'esecuzione dagli artefatti esistenti, rifacendo solo ciò che serve |
| `report` | genera il report dell'esecuzione |

Non servono variabili d'ambiente. I dati devono essere raggiungibili dal container:
sotto la radice del clone, oppure montati con un altro `-v`. I percorsi relativi della
configurazione valgono rispetto alla cartella da cui il comando è lanciato.
`scripts/esegui.py` esegue il codice e gli script R del clone; il comando `amplicon16s`
dentro l'immagine esegue invece il codice copiato quando l'immagine è stata costruita.

`run` non sovrascrive mai un'esecuzione: se la cartella di uscita non è vuota si
rifiuta, e si sceglie fra `resume` e una cartella nuova.

## Configurazione

Un file YAML. `config/config.example.yaml` è il modello: elenca ogni parametro e dice,
per ciascuno, che cosa significa e come si sceglie.

- **I parametri che dipendono dai dati o dallo studio sono obbligatori** e non hanno un
  valore predefinito: colonne ed etichette dei metadati, primer e motivo conservato,
  troncamento, riferimento, scelte sui controlli, soglie di qualità. Se ne manca uno la
  pipeline si ferma prima di ogni calcolo ed elenca tutti quelli che mancano
  (`E-G15-10`). Un parametro non pertinente si dichiara nullo o vuoto, dove il modello
  lo indica.
- **Hanno un predefinito** solo i parametri il cui valore è quello standard del metodo
  (il predefinito di dada2 o di decontam, con la fonte indicata nel modello), quelli di
  cui la pipeline realizza una sola forma e quelli che governano l'esecuzione.
- Lo schema rifiuta le chiavi che non conosce: un refuso ferma l'esecuzione.

Per cominciare da un dataset simile a uno dei due del repository conviene partire dalla
sua configurazione (`dati/osd734/config_osd734.yaml`, `dati/osd276/config_osd276.yaml`)
e rivedere ogni valore: quelli scritti lì sono stati scelti per quei dati.

`run.strict_provenance: true` fa partire l'esecuzione solo da un clone senza modifiche
non committate in `src/` e `R/` e in un ambiente R identico a `renv.lock`: è il modo di
produrre risultati da consegnare.

## Ingressi

| Ingresso | Parametro | Note |
|---|---|---|
| Cartella delle letture | `io.fastq_dir`, `io.fastq_glob` | un file per campione, senza sottocartelle |
| Tabella di assay | `io.assay_table` | nome del campione e file delle letture |
| Tabella di studio | `io.study_table` | facoltativa: tipo del campione, posizione, variabili |
| File del lotto | `io.batch_table` | facoltativo: piastra, corsa, cellule dei controlli positivi |
| Riferimento tassonomico | `tax.ref_fasta`, `tax.ref_md5` | l'MD5 è verificato a ogni avvio |

Il campione si riconosce da una chiave estratta dal nome del file con
`io.accession_regex`, la stessa che si ricava dalla tabella di assay: deve essere unica
per campione.

## Uscite

Sotto `io.out_root`, una cartella per fase:

| Cartella | Contenuto |
|---|---|
| `00_config/` | la configurazione usata, con il suo digest e l'elenco dei parametri dichiarati |
| `01_input_validation/` | esito dei controlli di S0, inventario dei campioni |
| `02_qc_profiles/` | profili di qualità e di lunghezza, conteggi delle letture |
| `03_filtered/` | letture filtrate |
| `04_error_models/` ... `08_taxonomy/` | modelli d'errore, varianti, tabella delle sequenze, tassonomia |
| `10_phyloseq/` | oggetto integrato con tutti i campioni |
| `11_controls/` | curve di calibrazione, soglie di profondità, contaminanti |
| `12_final/` | **i risultati da usare** |
| `99_logs/` | log strutturato, registro degli avvii, punto di ripresa |
| `report/` | il report, dopo `report` |

In `12_final/`:

| File | Contenuto |
|---|---|
| `ps_final.rds` | oggetto `phyloseq` con i soli campioni biologici che hanno superato i filtri |
| `ps_controlli.rds` | i controlli, se il dataset ne ha |
| `conteggi.tsv`, `tassonomia.tsv`, `metadati.tsv`, `sequenze.fasta` | gli stessi dati in testo, sufficienti a ricostruire l'oggetto |
| `checksum.sha256` | SHA-256 dei file consegnati |

Ogni cartella ha un manifesto con dimensione e SHA-256 di ogni artefatto. Le varianti si
chiamano `ASV1`, `ASV2`, ... per abbondanza decrescente nei campioni biologici; i nomi
non cambiano dopo i filtri, quindi la numerazione finale ha dei vuoti.

## Come si legge il report

`report` scrive `report/report.html` e, in `report/tabelle/`, ogni sua tabella come file
TSV. Lo ricava da ciò che l'esecuzione ha lasciato su disco: non ricalcola nulla e si può
generare anche su un'esecuzione conclusa da tempo.

Le sezioni, nell'ordine in cui compaiono:

- **Segnalazioni in apertura**: ciò che va visto per primo. Arresti, avvisi dichiarati
  dalle fasi con il loro codice, tentativi ripetuti, differenze di provenienza.
- **Stato della catena**: quali fasi sono concluse.
- **Gate di validazione e controlli di avvio**: l'esito di ogni controllo di S0.
- **Configurazione**: ogni parametro con il suo valore e l'origine, cioè *dichiarato*
  nel file, *predefinito*, *derivato* da altri parametri, oppure *aggiustato* da un
  tentativo ripetuto.
- **Provenienza delle fasi**: versione del codice, impronta del sorgente e immagine.
- **Decisioni prese automaticamente**: curve di calibrazione e soglia di profondità per
  piastra, contaminanti rimossi, campioni riclassificati, ripieghi dichiarati.
- **Tracciamento delle letture**: quante letture ha ogni campione dopo ogni fase, e
  dove si perdono.
- **Risultato finale**: campioni e varianti dell'oggetto finale, campioni esclusi con
  il motivo di ciascuno.

Ogni avviso o arresto porta un codice (`E-S11-02`, `E-G15-10`, ...): il messaggio riporta
il dettaglio del caso e che cosa fare.

## Codici di uscita

| Codice | Significato |
|---|---|
| 0 | successo |
| 1 | errore imprevisto: un difetto del programma; la traccia è nel log in `99_logs/` |
| 2 | riga di comando non valida |
| 3 | configurazione non valida, oppure `run` su una cartella di uscita già usata; nessuna fase è partita |
| 4 | arresto con punto di ripresa dichiarato, stampato e scritto in `99_logs/punto_di_ripresa.json` |
| 5 | una fase prevista dal grafo non è realizzata; non si presenta con le fasi attuali |

Un arresto con codice 4 dice quale fase si è fermata, perché e che cosa fare; corretta
la causa, `resume` riparte da lì. Gli avvisi non fermano l'esecuzione: restano nel log,
nel manifesto della fase e nel report.

## Riprodurre il riferimento

Il dataset di riferimento è OSD-734 (NASA GeneLab/OSDR). `dati/osd734/` contiene i
metadati, l'elenco delle letture con i checksum della fonte, la configurazione e i
checksum attesi di tutti gli artefatti dell'esecuzione pubblicata.

```bash
python3 dati/osd734/scarica_letture.py            # 960 file FASTQ da ENA, 2,52 GB
python3 dati/riferimento/scarica_riferimento.py   # riferimento tassonomico da Zenodo
python3 dati/osd734/ricostruisci_lotto.py         # file del lotto dalla fonte degli autori

docker run --rm -u "$(id -u):$(id -g)" -v "$(pwd)":/app -w /app <immagine> \
  python3 scripts/esegui.py run --config dati/osd734/config_osd734.yaml

python3 dati/osd734/confronta_risultati.py        # 1.033 artefatti, fase per fase
```

Il confronto termina con `RISULTATI IDENTICI` oppure indica il primo artefatto diverso.
La configurazione ha la regola rigorosa attiva: va eseguita da un clone senza modifiche.
I dettagli, i valori da adattare alla propria macchina e i requisiti misurati sono in
`dati/osd734/README.md`; il secondo dataset è in `dati/osd276/README.md`.

## Limiti noti

Ingressi e configurazione:

- Un file con le due letture di ogni coppia intercalate è riconosciuto dai marcatori
  `/1` e `/2`, dal commento `1:N:` e `2:N:`, dal nome ripetuto, e dal marcatore `1` e `2`
  o `R1` e `R2` in fondo all'identificativo dopo `/`, `.`, `_` o `-`. Non è riconosciuto
  se contiene una sola coppia, se usa un marcatore diverso da questi, o se usa i
  marcatori in fondo all'identificativo ma riporta prima tutte le prime letture e poi
  tutte le seconde, o le sole seconde.
- Un file il cui nome finisce con `R2` prima dell'estensione, o con `_2` quando nella
  cartella c'è lo stesso nome con `_1`, è preso per il secondo di una coppia e respinto.
- Un dataset single-end le cui letture portano tutte il marcatore `/2` è respinto.
- Due file `_R1` e `_R2` con la stessa chiave sono respinti come chiave ripetuta
  (`E-S0-05`), non come dati appaiati.
- Un file di letture vuoto ferma la validazione (`E-S0-13`), anche se è di un controllo
  negativo.
- `ctrl.column` non può essere nullo: un dataset senza controlli indica comunque una
  colonna con il tipo dei campioni.
- Con `filter.trimLeft` maggiore di zero `qc.conserved_motif` si cerca dopo le basi
  tagliate: un'espressione che guarda all'indietro non trova nulla.
- `qc.conserved_motif` e `qc.primer_sequence` vanno scritti in maiuscolo.
- I parametri reali accettano `true` e `false`, letti come 1 e 0.
- `qc.primer_sequence` non può essere nullo.
- Il file di `tax.ref_bad_taxa` deve avere le colonne `name` e `rank`: con altre S8 si
  ferma con un errore senza codice.
- Per alcuni codici di S0 e per `E-S11-03` il testo «che cosa fare» descrive il caso più
  frequente; il dettaglio dell'errore è sempre quello del caso.
- Una riga finale di soli spazi in un FASTQ ferma S2 (`E-S2-03`); byte estranei dopo la
  fine dell'archivio gzip fermano S0 o S1.

Risultati:

- S8 cerca il rango `Phylum`: un riferimento con ranghi nominati diversamente dà un
  errore senza codice.
- La ricerca dell'albero filogenetico (S9) è locale: l'albero serve alle distanze fra
  campioni, non come risultato filogenetico.
- La tassonomia arriva al genere: l'assegnazione della specie non è realizzata.

## Test

La suite verifica la pipeline su dataset sintetici e su un sottoinsieme ridotto di
OSD-734 in `tests/fixtures/`. Nell'immagine, dalla radice del clone:

```bash
docker run --rm -e PYTHONPATH=/app/src -e AMPLICON16S_R_DIR=/app/R \
  -v "$(pwd)":/app -w /app <immagine> pytest -o cache_dir=/tmp/.pytest_cache
```

In un ambiente Python 3.11 senza R (`pip install -e ".[dev]"`, poi `pytest`) i test che
richiedono R o Bioconductor si saltano e dicono perché.
