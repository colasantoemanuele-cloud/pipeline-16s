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
(più sotto). **Per un dataset proprio la strada è in [GUIDA.md](GUIDA.md)**, passo per
passo, dai dati al report; che cosa significa ogni parametro, e i casi particolari, è
in [RIFERIMENTO.md](RIFERIMENTO.md).

## Dati supportati

- Letture Illumina dell'amplicone 16S, in FASTQ compresso con gzip o no: **single-end**, un
  file per campione, o **appaiate** (`io.layout: paired-end`), due file per campione, fuse dopo
  l'inferenza. Un file con le due letture di ogni coppia insieme viene respinto (RIFERIMENTO.md, 1.2).
- Letture **senza primer**, con il primer a posizione fissa (`filter.trimLeft`), o con il primer dopo un
  numero variabile di basi, cercato lungo la lettura (`filter.primer_removal: search`; RIFERIMENTO.md, 6.3).
- Metadati in tabelle di testo UTF-8 separate da tabulazioni: una tabella dei campioni,
  oppure le tabelle di assay e di studio del formato ISA-Tab degli archivi OSDR.
- Facoltativi: un file del lotto (piastra di estrazione e corsa di sequenziamento per
  campione), controlli negativi, controlli positivi in diluizione seriale. Senza di
  essi la pipeline procede e dichiara che cosa non ha potuto fare. Un controllo
  negativo senza alcuna lettura è ammesso e attraversa le fasi come tale.
- Un riferimento tassonomico con il suo MD5, verificato a ogni avvio: un training set di dada2 (SILVA, RDP, GTDB,
  Greengenes) per `assignTaxonomy`, o un classificatore addestrato di DECIPHER per IdTaxa (`tax.classifier`; RIFERIMENTO.md, 1.5).

File e metadati si agganciano con una chiave estratta dal nome di ogni file, unica per
campione (`io.accession_regex`). Come preparare letture e tabelle, e come scaricare uno studio di OSDR
(`scripts/scarica_osdr.py`), è in GUIDA.md, passo 2, e in RIFERIMENTO.md, sezione 1.

Nel repository ci sono i dataset su cui la pipeline è eseguita, tutti da Illumina MiSeq, ciascuno con il suo README:

| Cartella | Che cosa contiene |
|---|---|
| `dati/osd734/` | il dataset di riferimento: 960 campioni, controlli negativi e positivi, file del lotto, impronte attese |
| `dati/osd276/` | 15 campioni, senza controlli; configurazioni forward e con IdTaxa |
| `dati/osd170/` | 30 campioni, soli controlli negativi; configurazioni forward, appaiata e con IdTaxa |
| `dati/osd191/` | letture lunghe (2 x 300) della regione V3-V4: soli dati e scarico, nessuna configurazione |
| `dati/comunita_nota/`, `dati/comunita_nota_v3v4/` | uno standard di composizione nota, regioni V4 e V3-V4, con lo script che ne verifica la composizione |
| `dati/confronto_genelab/` | lo script che confronta, al livello di genere, un'esecuzione con le tabelle elaborate da NASA GeneLab |
| `dati/riferimento/` | scarico e impronte del riferimento tassonomico (SILVA) e dei classificatori addestrati per IdTaxa |
| `dati/scaricamento.py`, `dati/letture_forward.py` | funzioni comuni agli script di scarico: scarico verificato e riprendibile, scelta del file forward dal primer |

## Requisiti

- **Docker** (o Apptainer, RIFERIMENTO.md 3.1) e l'immagine pubblicata (1,86 GB da scaricare, 7,6 GB su
  disco), per l'architettura `linux/amd64`. È l'unica via supportata: l'immagine porta Python 3.11, R 4.5.2, Bioconductor 3.21 e i
  pacchetti alle versioni di `renv.lock`, fra cui dada2 1.36.0.1, cioè dada2 1.36.0 con
  una correzione che rende ripetibile la tassonomia nei pareggi fra generi
  (`container/dada2/`).
- Un clone di questo repository (`git clone https://github.com/colasantoemanuele-cloud/pipeline-16s.git`):
  il codice si esegue dal clone montato nel container.
- Memoria e tempo in proporzione ai dati e ai processori: con 12 processori, 960 campioni
  single-end circa 14 GB di memoria di picco e 45-50 minuti; 28 campioni appaiati circa 15 GB
  e 25 minuti.
- Spazio su disco: per un dataset proprio, nella cartella di uscita deve essere libero almeno
  quanto occupano i file delle letture; è ciò che la validazione verifica (`E-S0-14`), perché
  gli artefatti intermedi sono dello stesso ordine. Sul dataset di 960 campioni (2,5 GB di
  letture) le uscite occupano 2,3 GB.

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
sotto la radice del clone, oppure montati con un altro `-v` (RIFERIMENTO.md, 2). I
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

Un file YAML, che di solito scrive **`init`** a partire dalle proprie letture e dai propri
metadati: lo strumento propone, chi conduce l'analisi conferma, e nessun valore entra senza
un gesto. Il percorso, con i comandi, è in GUIDA.md, passi 4-8; ogni parametro, i marcatori
del file e i parametri con un predefinito sono in RIFERIMENTO.md (sezioni 2 e 6).
`config/config.example.yaml` è il modello con tutti i parametri, e le configurazioni in
`dati/` sono esempi commentati.

- I parametri che dipendono dai dati o dallo studio sono **obbligatori**, senza predefinito:
  se uno manca la pipeline li elenca e si ferma (`E-G15-10`). I **condizionali** si omettono
  se non pertinenti (RIFERIMENTO.md, 2.5).
- Lo schema rifiuta le chiavi che non conosce: un refuso ferma l'esecuzione.
- `run.strict_provenance: true` fa partire l'esecuzione solo da un clone senza modifiche non
  committate in `src/` e `R/` e in un ambiente R identico a `renv.lock`: è il modo di
  produrre risultati da consegnare.

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
Gli artefatti di ogni fase e a che cosa servono: RAZIONALE.md, sezione 3.

## Come si legge il report

`report` scrive `report/report.html` e, in `report/tabelle/`, ogni sua tabella in TSV.
Lo ricava da ciò che l'esecuzione ha lasciato su disco, senza ricalcolare nulla. Le
sezioni, nell'ordine:

- **Segnalazioni in apertura**: arresti, tentativi ripetuti, *degradazioni* (gli avvisi con cui
  la pipeline prosegue rinunciando a qualcosa), integrità di ciò che è stato letto e controlli di
  qualità spenti. Poi **Stato della catena** e **Gate di validazione e controlli di avvio**.
- **Configurazione** (ogni parametro con il valore e l'origine: *dichiarato*,
  *predefinito*, *derivato*, *aggiustato*) e **Provenienza delle fasi**.
- **Decisioni prese automaticamente**: troncamento suggerito, copertura del phylum e varianti
  assegnate per rango, curve e soglia di profondità per piastra, contaminanti e campioni esclusi.
- **Tracciamento delle letture** (segnala un controllo negativo con più letture dei
  biologici) e **Risultato finale** (campioni e varianti, esclusi con il motivo).

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
in `dati/osd734/`, `dati/osd276/` e `dati/osd170/` (`eco_*.yaml`); GUIDA.md (passo 10), RIFERIMENTO.md (4).

## Integrazione con altri programmi

Tre comandi scrivono sullo standard output un documento JSON con il campo
`versione_formato`, per la pipeline e, con `eco` davanti, per le analisi ecologiche:

| Comando | Contenuto |
|---|---|
| `schema` | ogni parametro della configurazione: tipo, obbligatorietà, valori ammessi, predefinito, descrizione |
| `validate --json` | per ogni problema: parametro, codice, messaggio, azione suggerita |
| `status --json` | fasi concluse, fase in corso, arresto con il codice e il punto di ripresa |

`status` legge soltanto, anche durante un'esecuzione. Formati, e `soglie --json`: RIFERIMENTO.md, 5.

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

- Le letture inverse prese per forward e i file scambiati si riconoscono dal contenuto solo con un motivo conservato o un primer da cercare.
- Le coppie si abbinano dal marcatore in fondo al nome del file; in modalità appaiata i taxa con la regione più lunga di quanto le letture coprono non si fondono, senza errore.
- `init` propone solo ciò che i dati sostengono senza ambiguità: riconosce i primer di un elenco, e il resto lo lascia da scegliere.
- I controlli positivi si valutano solo come serie di diluizione di un taxon; la tassonomia non assegna la specie, e i filtri tolgono solo i taxa che il riferimento nomina.
- Una soglia di qualità stretta a mano si controlla solo al calcolo successivo; l'unico ambiente supportato è l'immagine pubblicata.
- Il dettaglio di ogni limite è in RIFERIMENTO.md, sezione 7; i limiti scientifici dei metodi in RAZIONALE.md, sezione 5.

## Test

La suite verifica pipeline e analisi su dataset sintetici e sul sottoinsieme di prova.
Nell'immagine, dalla radice del clone (senza R i test che lo richiedono si saltano):

```bash
docker run --rm -u "$(id -u):$(id -g)" -e HOME=/tmp -e PYTHONPATH=/app/src \
  -e AMPLICON16S_R_DIR=/app/R -v "$(pwd)":/app -w /app <immagine> pytest -o cache_dir=/tmp/.pytest_cache
```

## Segnalare un problema

Si apre una issue su GitHub, nel repository da cui viene il clone, allegando l'uscita di
`status --out <cartella> --json`, il file `99_logs/punto_di_ripresa.txt` (se c'è) e
`00_config/resolved.yaml`.

## Versioni

Una versione dello strumento è un commit (o un tag) di questo repository insieme al digest
della sua immagine: i due insieme identificano i risultati.

## Documenti

- `README.md`: che cos'è lo strumento.
- [GUIDA.md](GUIDA.md): l'uso passo per passo.
- [RIFERIMENTO.md](RIFERIMENTO.md): ogni parametro, i formati, i casi particolari, il glossario.
- [RAZIONALE.md](RAZIONALE.md): il perché scientifico di ogni metodo e i suoi limiti.
- [SVILUPPO.md](SVILUPPO.md): il codice.
