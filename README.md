# amplicon16s

Pipeline per l'analisi di dati 16S rRNA da sequenziamento Illumina, single-end o appaiato.

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
| S2 | taglio del primer, se cercato lungo la lettura; filtro e troncamento (`dada2::filterAndTrim`) |
| S3-S7 | modello d'errore per corsa, inferenza delle varianti, in modalità appaiata fusione delle coppie (SM, `dada2::mergePairs`), tabella delle sequenze, rimozione delle chimere, selezione per lunghezza |
| S8 | assegnazione tassonomica (`dada2::assignTaxonomy`, oppure `DECIPHER::IdTaxa`) |
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

- Letture Illumina dell'amplicone 16S, in FASTQ compresso con gzip o no: **single-end**, un
  file per campione, o **appaiate** (`io.layout: paired-end`), due file per campione, fuse dopo
  l'inferenza. Un file con le due letture di ogni coppia insieme viene respinto (GUIDA.md, 1.2).
- Letture **senza primer**, con il primer a posizione fissa (`filter.trimLeft`), o con il primer dopo un
  numero variabile di basi, cercato lungo la lettura (`filter.primer_removal: search`; GUIDA.md, 6.3).
- Metadati in tabelle di testo UTF-8 separate da tabulazioni: una tabella dei campioni,
  oppure le tabelle di assay e di studio del formato ISA-Tab degli archivi OSDR.
- Facoltativi: un file del lotto (piastra di estrazione e corsa di sequenziamento per
  campione), controlli negativi, controlli positivi in diluizione seriale. Senza di
  essi la pipeline procede e dichiara che cosa non ha potuto fare. Un controllo
  negativo senza alcuna lettura è ammesso e attraversa le fasi come tale.
- Un riferimento tassonomico con il suo MD5, verificato a ogni avvio: un training set di dada2 (SILVA, RDP, GTDB,
  Greengenes) per `assignTaxonomy`, o un classificatore addestrato di DECIPHER per IdTaxa (`tax.classifier`; GUIDA.md, 1.5).

File e metadati si agganciano con una chiave estratta dal nome di ogni file, unica per
campione (`io.accession_regex`). Come preparare letture e tabelle è in GUIDA.md,
sezione 1.

Nel repository ci sono i dataset su cui la pipeline è eseguita, con dati, configurazione e
istruzioni: `dati/osd734/` (960 campioni, controlli e file del lotto), `dati/osd276/` (15, senza
controlli), `dati/osd170/` (30, soli controlli negativi), `dati/comunita_nota/` (uno standard noto).

## Requisiti

- **Docker** (o Apptainer, GUIDA.md 3.1) e l'immagine pubblicata (1,86 GB da scaricare, 7,6 GB su
  disco). È l'unica via supportata: l'immagine porta Python 3.11, R 4.5.2, Bioconductor 3.21 e i
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
| `soglie` | da un'esecuzione conclusa, propone le soglie di qualità dalle sue misure |
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

Nel file ogni valore dichiara da dove viene: **rilevato** sui file (riga attiva), **proposto** da
un'euristica o **da scegliere** (righe in commento, con il valore proposto o di partenza quando
c'è). Nessun valore entra senza un gesto: si toglie il commento riga per riga, oppure si accettano
tutti con `--accetta-proposte`. Il percorso breve: `init ... --accetta-proposte` (proposte e
valori di partenza attivi, soglie di qualità al valore che non ferma); si scrivono le scelte che
`validate` elenca; `run`; `soglie --config <configurazione.yaml> --output <nuova.yaml>
--accetta-proposte` propone le soglie dalle misure di quell'esecuzione, e `resume --config
<nuova.yaml>` le attiva. Una soglia proposta protegge dai peggioramenti, non giudica la misura.
Dettagli in GUIDA.md (sezioni 2 e 6.6); le configurazioni in `dati/` sono esempi commentati.

- I parametri che dipendono dai dati o dallo studio sono **obbligatori**, senza predefinito
  (lo hanno solo standard del metodo, forme uniche e parametri di esecuzione): se uno manca
  la pipeline li elenca e si ferma (`E-G15-10`). I **condizionali** si omettono se non
  pertinenti: omessi o dichiarati con un valore danno gli stessi conteggi, tassonomia, campioni
  e varianti ma un digest diverso; sette, se omessi, sono nulli nei riepiloghi (GUIDA.md, 2.5).
- Lo schema rifiuta le chiavi che non conosce: un refuso ferma l'esecuzione.
- `run.strict_provenance: true` fa partire l'esecuzione solo da un clone senza
  modifiche non committate in `src/` e `R/` e in un ambiente R identico a `renv.lock`:
  è il modo di produrre risultati da consegnare.

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

- **Segnalazioni in apertura**: arresti, tentativi ripetuti, *degradazioni* (gli avvisi con cui
  la pipeline prosegue rinunciando a qualcosa), integrità di ciò che è stato letto e controlli di
  qualità spenti. Poi **Stato della catena** e **Gate di validazione e controlli di avvio**.
- **Configurazione** (ogni parametro con il valore e l'origine: *dichiarato*,
  *predefinito*, *derivato*, *aggiustato*) e **Provenienza delle fasi**.
- **Decisioni prese automaticamente**: curve di calibrazione e soglia di profondità per
  piastra, contaminanti rimossi, campioni riclassificati, troncamento suggerito.
- **Tracciamento delle letture** e **Risultato finale** (campioni e varianti, esclusi
  con il motivo).

## Codici di uscita

| Codice | Significato |
|---|---|
| 0 | successo |
| 1 | errore imprevisto: un difetto del programma; la traccia è nel log in `99_logs/` |
| 2 | riga di comando non valida |
| 3 | configurazione non valida, `run` su una cartella già usata, o `init` che si rifiuta; nessuna fase è partita |
| 4 | arresto con punto di ripresa, stampato e scritto in `99_logs/punto_di_ripresa.json` |
| 5 | una fase prevista dal grafo non è realizzata; non si presenta con le fasi attuali |
| 6 | la cartella di uscita è tenuta da un'altra esecuzione (`E-RUN-01`): nulla è stato scritto |

Corretta la causa di un arresto, `resume` riparte da lì. Gli avvisi non fermano. Ogni messaggio porta un codice (`E-S11-02`, `E-G15-10`, ...), il dettaglio del caso e che cosa fare.

## Analisi ecologiche di base

Il pacchetto `amplicon16s_eco` lavora sull'oggetto finale, che non modifica, con una
configurazione propria (`config/eco.example.yaml`): alfa diversità su conteggi
rarefatti, composizione a un rango, distanze (Bray-Curtis, Jaccard, Aitchison, UniFrac
con un albero), PCoA e NMDS, confronto fra gruppi (Wilcoxon o Kruskal-Wallis, PERMANOVA
con il test delle dispersioni) e, se dichiarati, variabili continue (Spearman) e misure ripetute.

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
in `dati/osd734/`, `dati/osd276/` e `dati/osd170/` (`eco_*.yaml`); dettagli in GUIDA.md.

## Integrazione con altri programmi

Tre comandi scrivono sullo standard output un documento JSON con il campo
`versione_formato`, per la pipeline e, con `eco` davanti, per le analisi ecologiche:

| Comando | Contenuto |
|---|---|
| `schema` | ogni parametro della configurazione: tipo, obbligatorietà, valori ammessi, predefinito, descrizione |
| `validate --json` | per ogni problema: parametro, codice, messaggio, azione suggerita |
| `status --json` | fasi concluse, fase in corso, arresto con il codice e il punto di ripresa |

`status` legge soltanto, anche durante un'esecuzione. Formati, e `soglie --json`: GUIDA.md, 5.

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
misurati sono in `dati/osd734/README.md`; ogni altro dataset ha il suo README in `dati/`.

## Limiti noti

- Le letture inverse prese per forward si riconoscono dal contenuto (`E-S1-07`), con il motivo
  conservato o con il primer nelle basi tagliate: non con il motivo nullo e senza taglio, non in un
  file con meno di 500 letture o di un controllo negativo, non se il motivo è stato ricavato da
  file in maggioranza inversi, non se la libreria mescola i due versi (GUIDA.md, 1.2).
- In modalità single-end un file con il nome di una seconda lettura (`R2`, `_2`) è respinto per il
  nome solo dove il contenuto non si può giudicare (GUIDA.md, 1.2); lo è un dataset le cui letture
  portano tutte il marcatore `/2`. In modalità appaiata i taxa con la regione più lunga di quanto
  le letture troncate coprono non si fondono, senza alcun errore (GUIDA.md, 6.7).
- `init` riconosce il primer solo se è in un elenco di primer 16S comuni; la chiave del campione solo se una
  colonna corrisponde uno a uno ai file. Il primer cercato lungo la lettura ammette sole sostituzioni, non
  inserzioni né delezioni, e non toglie il primer dell'altra estremità in coda alle letture (GUIDA.md, 6.3).
- Con `filter.trimLeft` maggiore di zero `qc.conserved_motif` si cerca dopo le basi tagliate
  (un'espressione che guarda all'indietro non trova nulla) e il primer non si cerca: in
  `letture_ispezionate.tsv` la colonna `frazione_primer` vale 0.
- Byte estranei dopo la fine dell'archivio gzip di un FASTQ fermano S0 o S1, anche se
  altri strumenti leggono il file.
- Il testo generale «che cosa fare» di alcuni codici descrive solo il caso più frequente (quelli
  sui controlli negativi citano `decontam.min_blanks` anche se omesso): fa fede il dettaglio.
- Un riferimento senza cloroplasti, mitocondri ed eucarioti (GTDB) non permette di toglierli per nome. La tassonomia
  arriva al rango più fine del riferimento, con IdTaxa al genere. IdTaxa richiede un classificatore già addestrato, e può
  lasciare vuoto un rango sopra ranghi assegnati (GUIDA.md, 6.4). L'albero di S9 serve alle distanze fra campioni.
- Nelle analisi ecologiche una variabile è continua, e un disegno è a misure ripetute, solo se
  dichiarati (`design.continuous_variables`, `design.subject`: GUIDA.md, 4.1 e 4.2). Con una
  distanza non euclidea la PERMANOVA entro soggetto può non essere valutabile. Un'analisi respinta
  non lascia nulla nella cartella di uscita: `eco status` non la riporta, il codice lo dà il comando.

## Test

La suite verifica pipeline e analisi su dataset sintetici e sul sottoinsieme di prova.
Nell'immagine, dalla radice del clone (senza R i test che lo richiedono si saltano):

```bash
docker run --rm -u "$(id -u):$(id -g)" -e HOME=/tmp -e PYTHONPATH=/app/src \
  -e AMPLICON16S_R_DIR=/app/R -v "$(pwd)":/app -w /app <immagine> pytest -o cache_dir=/tmp/.pytest_cache
```
