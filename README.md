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
| S3-S7 | modello d'errore per corsa, inferenza delle varianti, tabella delle sequenze, rimozione delle chimere, selezione per lunghezza |
| S8 | assegnazione tassonomica (`dada2::assignTaxonomy`) |
| S10 | oggetto integrato con tutti i campioni |
| S11 | curva di calibrazione dai controlli positivi in diluizione (KatharoSeq) e soglia di profondità |
| S12 | riconoscimento dei contaminanti dai controlli negativi (`decontam`) |
| S13 | filtri finali: profondità, tassonomia, prevalenza, letture minime |
| S9, S14 | albero filogenetico delle varianti finali (facoltativo); oggetto finale, export e checksum |

Python orchestra l'esecuzione; R esegue i calcoli (dada2, phyloseq, DECIPHER, phangorn,
decontam) come processi separati. Ogni fase scrive i propri artefatti con un manifesto
che ne riporta dimensione e SHA-256: un'esecuzione interrotta si riprende dal punto in
cui si è fermata, e a parità di configurazione, codice e immagine due esecuzioni danno
gli stessi byte.

Sull'oggetto finale lavora un secondo pacchetto, le **analisi ecologiche di base**
(più sotto). **Per un dataset proprio la strada è in [GUIDA.md](GUIDA.md)**: come
preparare i dati, come scrivere la configurazione con `init`, che cosa significa ogni
parametro.
## Dati supportati

- Letture Illumina **single-end** dell'amplicone 16S, un file FASTQ per campione,
  compresso con gzip o no. Di una corsa appaiata si usano le sole letture forward; un
  file che contiene le due letture di ogni coppia viene respinto.
- Letture **senza primer in testa**, oppure con il primer a posizione fissa, che si
  toglie con `filter.trimLeft`.
- Metadati in tabelle di testo UTF-8 separate da tabulazioni: una tabella dei campioni,
  oppure le tabelle di assay e di studio del formato ISA-Tab degli archivi OSDR.
- Facoltativi: un file del lotto (piastra di estrazione e corsa di sequenziamento per
  campione), controlli negativi, controlli positivi in diluizione seriale. Senza di
  essi la pipeline procede e dichiara che cosa non ha potuto fare. Un controllo
  negativo senza alcuna lettura è ammesso e attraversa le fasi come tale.
- Un riferimento tassonomico nel formato dei training set di dada2, con il suo MD5,
  verificato a ogni avvio.

File e metadati si agganciano con una chiave estratta dal nome di ogni file, unica per
campione (`io.accession_regex`). Come preparare letture e tabelle è in GUIDA.md,
sezione 1.

Nel repository ci sono due dataset su cui la pipeline è eseguita, con dati,
configurazione e istruzioni: `dati/osd734/` (960 campioni della Stazione Spaziale
Internazionale, con controlli e file del lotto) e `dati/osd276/` (15, senza controlli).

## Requisiti

- **Docker**, e l'immagine pubblicata (1,86 GB da scaricare, 7,6 GB su disco). È
  l'unica via supportata: l'immagine porta Python 3.11, R 4.5.2, Bioconductor 3.21 e i
  pacchetti alle versioni di `renv.lock`, fra cui dada2 1.36.0.1, cioè dada2 1.36.0 con
  una correzione che rende ripetibile la tassonomia nei pareggi fra generi
  (`container/dada2/`).
- Un clone di questo repository: il codice si esegue dal clone montato nel container.
- Memoria e spazio in proporzione ai dati: sul dataset di 960 campioni 13 GB di memoria
  di picco, 2,3 GB di uscite, 50 minuti con 12 processori.

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
| `init` | scrive una configurazione di partenza a partire dai dati |
| `validate` | esegue la sola validazione (S0) |
| `run` | esegue dall'inizio, in una cartella di uscita nuova |
| `status` | dice a che punto è un'esecuzione, anche mentre è in corso |
| `resume` | riprende un'esecuzione dagli artefatti esistenti, rifacendo solo ciò che serve |
| `report` | genera il report dell'esecuzione |
| `schema` | esporta in JSON lo schema della configurazione |
| `eco ...` | le analisi ecologiche di base: `eco init`, `eco validate`, `eco run`, `eco status`, `eco schema` |

Non servono variabili d'ambiente. I dati devono essere raggiungibili dal container:
sotto la radice del clone, oppure montati con un altro `-v` (GUIDA.md, sezione 2). I
percorsi relativi della configurazione valgono rispetto alla cartella da cui il comando
è lanciato. `scripts/esegui.py` esegue il codice e gli script R del clone; il comando
`amplicon16s` dentro l'immagine esegue il codice copiato alla sua costruzione. `run`
non sovrascrive: se la cartella di uscita non è vuota si rifiuta, e si sceglie fra
`resume` e una cartella nuova.

### Avvio rapido sul sottoinsieme di prova

Il repository contiene un sottoinsieme di 28 campioni con poche migliaia di letture
ciascuno e un riferimento tassonomico sintetico, con la loro configurazione completa.
La catena si conclude in pochi minuti e non richiede di scaricare nulla oltre
all'immagine:

```bash
docker run --rm -u "$(id -u):$(id -g)" -v "$(pwd)":/app -w /app <immagine> \
  python3 scripts/esegui.py run --config tests/fixtures/osd734/ridotto/config_avvio_rapido.yaml
docker run --rm -u "$(id -u):$(id -g)" -v "$(pwd)":/app -w /app <immagine> \
  python3 scripts/esegui.py status --out output/avvio_rapido
docker run --rm -u "$(id -u):$(id -g)" -v "$(pwd)":/app -w /app <immagine> \
  python3 scripts/esegui.py report --config tests/fixtures/osd734/ridotto/config_avvio_rapido.yaml
```

Le uscite sono in `output/avvio_rapido/` (oggetto finale in `12_final/`, report in
`report/report.html`): servono a provare lo strumento, non a trarre conclusioni.

## Configurazione

Un file YAML. `config/config.example.yaml` è il modello: elenca ogni parametro e dice,
per ciascuno, che cosa significa e come si sceglie. Il modo consigliato di cominciare
è **`init`**, che scrive il file a partire dalle proprie letture e dai propri metadati:

```bash
docker run --rm -u "$(id -u):$(id -g)" -v "$(pwd)":/app -w /app <immagine> \
  python3 scripts/esegui.py init --reads <cartella delle letture> \
    --metadata <tabella dei campioni o cartella ISA-Tab> --output <configurazione.yaml>
```

Nel file ogni valore dichiara da dove viene: **rilevato** sui file (riga attiva),
**proposto** da un'euristica (riga in commento con il valore già scritto: per
accettarlo si toglie il commento), oppure **da scegliere** (riga in commento, con un
valore di partenza e il suo motivo, o senza valore). Nessun valore entra nella
configurazione senza un gesto di chi la scrive; `validate` elenca ciò che resta da
dichiarare. Argomenti, regole e significato di ogni parametro sono in GUIDA.md.

- I parametri che dipendono dai dati o dallo studio sono **obbligatori** e senza
  predefinito: se ne manca uno la pipeline si ferma prima di ogni calcolo e li elenca
  (`E-G15-10`). Uno non pertinente si dichiara nullo o vuoto. Hanno un predefinito solo
  i valori standard del metodo, le forme uniche e i parametri di esecuzione.
- Lo schema rifiuta le chiavi che non conosce: un refuso ferma l'esecuzione.
- `run.strict_provenance: true` fa partire l'esecuzione solo da un clone senza
  modifiche non committate in `src/` e `R/` e in un ambiente R identico a `renv.lock`:
  è il modo di produrre risultati da consegnare.

Le configurazioni dei due dataset del repository, in `dati/`, sono esempi commentati.

## Uscite

Sotto `io.out_root` c'è una cartella per fase (`00_config/`, `01_input_validation/`,
... `11_controls/`), ciascuna con un manifesto che riporta dimensione e SHA-256 di ogni
artefatto; `99_logs/` con il log, il registro degli avvii e il punto di ripresa;
`report/` dopo `report`. I risultati da usare sono in `12_final/`:

| File | Contenuto |
|---|---|
| `ps_final.rds` | oggetto `phyloseq` con i soli campioni biologici che hanno superato i filtri |
| `ps_controlli.rds` | i controlli, se il dataset ne ha |
| `conteggi.tsv`, `tassonomia.tsv`, `metadati.tsv`, `sequenze.fasta` | gli stessi dati in testo, sufficienti a ricostruire l'oggetto |
| `checksum.sha256` | SHA-256 dei file consegnati |

Le varianti si chiamano `ASV1`, `ASV2`, ... per abbondanza decrescente nei biologici;
i nomi non cambiano dopo i filtri, quindi la numerazione finale ha dei vuoti.

## Come si legge il report

`report` scrive `report/report.html` e, in `report/tabelle/`, ogni sua tabella in TSV.
Lo ricava da ciò che l'esecuzione ha lasciato su disco, senza ricalcolare nulla. Le
sezioni, nell'ordine:

- **Segnalazioni in apertura**: arresti, tentativi ripetuti e *degradazioni*, cioè gli
  avvisi con cui la pipeline prosegue rinunciando a qualcosa, ciascuno con il suo codice.
  Poi **Stato della catena** e **Gate di validazione e controlli di avvio**.
- **Configurazione** (ogni parametro con il valore e l'origine: *dichiarato*,
  *predefinito*, *derivato*, *aggiustato*) e **Provenienza delle fasi**.
- **Decisioni prese automaticamente**: curve di calibrazione e soglia di profondità per
  piastra, contaminanti rimossi, campioni riclassificati, troncamento suggerito.
- **Tracciamento delle letture**, **Risultato finale** (campioni e varianti, esclusi
  con il motivo), **Integrità di ciò che è stato letto**.

Ogni avviso o arresto porta un codice (`E-S11-02`, `E-G15-10`, ...): il messaggio
riporta il dettaglio del caso e che cosa fare.

## Codici di uscita

| Codice | Significato |
|---|---|
| 0 | successo |
| 1 | errore imprevisto: un difetto del programma; la traccia è nel log in `99_logs/` |
| 2 | riga di comando non valida |
| 3 | configurazione non valida, `run` su una cartella già usata, o `init` che si rifiuta; nessuna fase è partita |
| 4 | arresto con punto di ripresa, stampato e scritto in `99_logs/punto_di_ripresa.json` |
| 5 | una fase prevista dal grafo non è realizzata; non si presenta con le fasi attuali |

Corretta la causa di un arresto, `resume` riparte da lì. Gli avvisi non fermano.

## Analisi ecologiche di base

Il pacchetto `amplicon16s_eco` lavora sull'oggetto finale, che non modifica, con una
configurazione propria (`config/eco.example.yaml`): alfa diversità su conteggi
rarefatti, composizione a un rango, distanze (Bray-Curtis, Jaccard, Aitchison, UniFrac
con un albero), PCoA e NMDS, confronto fra gruppi (Wilcoxon o Kruskal-Wallis, PERMANOVA
con il test delle dispersioni).

```bash
docker run --rm -u "$(id -u):$(id -g)" -v "$(pwd)":/app -w /app <immagine> \
  python3 scripts/esegui.py eco init --object output/avvio_rapido/12_final/ps_final.rds \
    --output output/eco_avvio_rapido.yaml
docker run --rm -u "$(id -u):$(id -g)" -v "$(pwd)":/app -w /app <immagine> \
  python3 scripts/esegui.py eco run --object <ps_final.rds> --config <eco.yaml> --out <cartella>
```

`eco init` scrive una configurazione di partenza con le proposte e le scelte da fare
(completate quelle, `eco validate` la controlla contro l'oggetto); `eco run` scrive
tabelle, grafici, un riepilogo con gli avvisi (`E-ECO-`) e un manifesto con lo SHA-256
di ogni file, identici byte per byte a parità di oggetto, configurazione e seme. Esempi
in `dati/osd734/eco_osd734.yaml` e `dati/osd276/eco_osd276.yaml`; dettagli in GUIDA.md.

## Integrazione con altri programmi

Tre comandi scrivono sullo standard output un documento JSON con il campo
`versione_formato`, per la pipeline e, con `eco` davanti, per le analisi ecologiche:

| Comando | Contenuto |
|---|---|
| `schema` | ogni parametro della configurazione: tipo, obbligatorietà, valori ammessi, predefinito, descrizione |
| `validate --json` | per ogni problema: parametro, codice, messaggio, azione suggerita |
| `status --json` | fasi concluse, fase in corso, arresto con il codice e il punto di ripresa |

`status` legge soltanto: si può usare durante un'esecuzione. Il formato dei tre
documenti è in GUIDA.md, sezione 5.

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
La regola rigorosa è attiva: si esegue da un clone senza modifiche. Dettagli e requisiti
misurati sono in `dati/osd734/README.md`; il secondo dataset è in `dati/osd276/README.md`.

## Limiti noti

- Un file con le due letture di ogni coppia è riconosciuto dai marcatori `/1` e `/2`,
  dal commento `1:N:` e `2:N:`, dal nome ripetuto, e dal marcatore `1` e `2` o `R1` e
  `R2` in fondo all'identificativo. Non lo è se contiene una sola coppia, se usa un
  altro marcatore, o se con i marcatori in fondo all'identificativo riporta prima tutte
  le prime letture e poi le seconde, o le sole seconde.
- Un file il cui nome finisce con `R2`, o con `_2` quando c'è lo stesso nome con `_1`,
  è preso per il secondo di una coppia e respinto; così un dataset single-end le cui
  letture portano tutte il marcatore `/2`. Due file `_R1` e `_R2` con la stessa chiave
  sono respinti come chiave ripetuta (`E-S0-05`).
- `init` riconosce il primer solo se è in un elenco di primer 16S comuni e sta in testa
  alle letture; la chiave del campione solo se una colonna corrisponde uno a uno ai
  file. Negli altri casi lascia la riga da scegliere, e lo dice.
- Con `filter.trimLeft` maggiore di zero `qc.conserved_motif` si cerca dopo le basi
  tagliate: un'espressione che guarda all'indietro non trova nulla. Il motivo va
  scritto in maiuscolo; `qc.primer_sequence` non può essere nullo.
- I parametri reali accettano `true` e `false`, letti come 1 e 0. Il file di
  `tax.ref_bad_taxa` deve avere le colonne `name` e `rank`. Una riga finale di soli
  spazi in un FASTQ ferma S2 (`E-S2-03`); byte estranei dopo l'archivio gzip, S0 o S1.
- S8 cerca il rango `Phylum`: un riferimento con ranghi nominati diversamente dà un
  errore senza codice. La tassonomia arriva al genere; l'albero di S9 serve alle
  distanze fra campioni, non come risultato filogenetico.
- Le analisi ecologiche trattano ogni variabile come categorica: una variabile che
  sembra continua dà un avviso (`E-ECO-29`), non un'analisi di regressione. Un'analisi
  respinta non lascia nulla nella cartella di uscita: `eco status` non la riporta, il
  codice lo dà il comando.

## Test

La suite verifica pipeline e analisi su dataset sintetici e sul sottoinsieme di prova.
Nell'immagine, dalla radice del clone (in un ambiente Python 3.11 senza R i test che
richiedono R si saltano e dicono perché):

```bash
docker run --rm -e PYTHONPATH=/app/src -e AMPLICON16S_R_DIR=/app/R \
  -v "$(pwd)":/app -w /app <immagine> pytest -o cache_dir=/tmp/.pytest_cache
```
