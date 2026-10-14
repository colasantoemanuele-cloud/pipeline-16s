# Guida per chi sviluppa amplicon16s

Questa guida serve a chi deve **leggere o modificare il codice**. Per usare lo strumento c'è `GUIDA.md`; per i parametri e i formati, `RIFERIMENTO.md`. Qui c'è solo ciò che serve per orientarsi: dove sta ogni cosa, come scorre un'esecuzione, come si lavora senza rompere nulla. I termini (gate, manifesto, artefatto, degradazione, digest...) sono nel glossario di `RIFERIMENTO.md`, sezione 8; in fondo a questa guida c'è l'indice dei file.

## 1. In tre righe

- **Python orchestra, R calcola.** Python legge la configurazione, controlla, decide quale fase eseguire e la lancia; ogni calcolo scientifico è uno script R eseguito come processo separato. Fra i due passano solo file.
- **Lo stato è il disco.** Ogni fase scrive i suoi file in una cartella con un manifesto (impronte SHA-256). Non c'è un database: «a che punto siamo» si rilegge dai manifesti.
- **Nessun valore di un dataset nel codice.** Ogni scelta è un parametro della configurazione. Un predefinito esiste solo se è lo standard del metodo (con la fonte), l'unica forma realizzata o un parametro dello strumento (risorse, ripresa, seme); tutto il resto è obbligatorio.

## 2. Prima di cominciare: l'ambiente

- **L'ambiente è l'immagine** costruita da `container/Dockerfile` (Python, R, Bioconductor e i pacchetti alle versioni di `renv.lock`): è l'unico supportato (README, «Requisiti»). Il codice si esegue dal clone montato nel container.
- **Prova veloce.** La catena completa sul sottoinsieme di prova gira in pochi minuti: `python3 scripts/esegui.py run --config tests/fixtures/osd734/ridotto/config_avvio_rapido.yaml` (README, «Avvio rapido»). È il modo più rapido per vedere tutte le cartelle di uscita.
- **Test.** La suite completa si lancia nell'immagine con il comando del README, sezione «Test», che punta `PYTHONPATH` e `AMPLICON16S_R_DIR` al clone (senza, si proverebbe la copia installata nell'immagine). Per la sola parte Python basta `pip install -e ".[dev]"` e `pytest`: i test che chiedono R si saltano. Un solo modulo, senza i quattro processi: `pytest tests/<modulo>.py -n 0`.
- **Convenzioni.** Nomi, commenti e messaggi sono in italiano. Negli script R le vocali accentate si scrivono con l'apostrofo (`e'`, `piu'`); nel Python prevale l'apostrofo; nei documenti `.md` si usa l'accento. Ogni file comincia con un commento che dice che cosa fa e perché.

## 3. La mappa

| Percorso | Che cosa contiene | Quando lo tocchi |
|---|---|---|
| `src/amplicon16s/` | la pipeline in Python (orchestrazione, controlli, report) | quasi sempre |
| `R/` | uno script per fase (`01_profile.R` … `14_finale.R`), due script di servizio `00_*.R` che non sono fasi (`00_ambiente.R` lo lanciano l'esecutore e S8, `00_classificatore.R` lo lancia S0 con IdTaxa), e le funzioni condivise in `R/lib/` | quando cambia un calcolo |
| `src/amplicon16s_eco/` | le analisi ecologiche: pacchetto separato, con i suoi script R in `R/` | quando cambiano le analisi ecologiche |
| `scripts/` | `esegui.py` (il punto d'ingresso: esegue `src/` e `R/` di questo clone); per lo sviluppo `registro_sorgente.py`, `build_test_subset.py`, `costruisci_riferimento_sintetico.py`, `genera_renv_lock.sh`; `scarica_osdr.py` per i dati | raramente |
| `tests/` | la suite (`pytest`), i dati di prova in `tests/fixtures/` e gli script R finti per provare il ponte in `tests/r_doppioni/` | sempre, insieme al codice |
| `config/` | le configurazioni di esempio della pipeline (`config.example.yaml`) e delle analisi ecologiche (`eco.example.yaml`), con ogni parametro commentato | quando aggiungi un parametro |
| `dati/` | i dataset di riferimento (configurazioni, metadati, impronte attese) e il riferimento tassonomico | quando cambia un risultato atteso |
| `container/` | `Dockerfile`, installazione e verifica dei pacchetti R; `dada2/` contiene lo script e la correzione con cui dada2 viene compilato nell'immagine | quando cambia l'ambiente |
| `renv.lock` | le versioni esatte dei pacchetti R | solo con `scripts/genera_renv_lock.sh` |
| `.github/workflows/ci.yml` | i due job della CI (il secondo in due parti) | quando aggiungi un modulo di test lento |

Dentro `src/amplicon16s/`:

| Modulo | Ruolo |
|---|---|
| `cli.py` | i sottocomandi (`validate`, `run`, `resume`, `status`, `init`, `soglie`, `report`, `schema`) e, in testa, i codici di uscita |
| `config/` | `schema.py` (parametri, pydantic, chiavi sconosciute respinte), `defaults.py` (predefiniti con la fonte), `resolve.py` (parametri derivati), `vista.py` (ciò che una fase può leggere) |
| `runner/` | `graph.py` (fasi e dipendenze), `project.py` (che cosa è concluso, letto dal disco), `executor.py` (esegue, ritenta, si ferma), `retry.py`, `provenienza.py`, `tracciamento.py` |
| `steps/` | una classe per fase, `s00_validate.py` … `s14_finale.py`, sulla base comune `base.py` |
| `rbridge/` | il ponte verso R: `payload.py` (i due file JSON), `runner.py` (lancio ed esito) |
| `gates/` | i 15 controlli di S0 (`g01_g15.py`) e il loro ordine (`registry.py`, G15 per primo) |
| `metadata/` | lettura delle tabelle, inventario file-campioni (`crosswalk.py`), classi dei campioni (`controls_map.py`) |
| `io_layer/` | cartelle di uscita e manifesti (`artifacts.py`), impronte (`checksums.py`), lettura dei FASTQ e riconoscimento delle coppie (`reads.py`) |
| `errors/` | il catalogo dei codici (`catalog.py`) e le eccezioni |
| `logging/` | il log strutturato (`logger.py`) |
| `init/`, `soglie/`, `report/` | i tre comandi che preparano la configurazione, propongono le soglie e scrivono il report |
| `interfaccia.py` | le uscite JSON per programmi esterni (`schema`, `validate --json`, `status --json`) |

## 4. Il filo rosso: che cosa succede con `run`

Segui questi passi aprendo i file indicati: è il modo più rapido per capire tutto il resto.

1. **L'avvio.** `scripts/esegui.py run --config X.yaml` mette `src/` del clone nel percorso di Python, punta `AMPLICON16S_R_DIR` a `R/` del clone e chiama `amplicon16s.cli`.
2. **La configurazione.** `cli.py` carica il file con `config/schema.py`: un parametro sconosciuto o fuori dominio ferma tutto qui. `config/resolve.py` calcola i parametri derivati.
3. **I controlli di avvio.** A ogni avvio, prima della prima fase, `runner/executor.py` esegue G15 (coerenza della configurazione), G14 (processori e disco) e G12 (impronta del riferimento tassonomico), più la verifica di commit e ambiente R se `run.strict_provenance` è attivo. Poi registra la configurazione risolta in `00_config/resolved.yaml`.
4. **Il grafo.** `runner/graph.py` dichiara le fasi, le loro cartelle e le dipendenze **di dato** (una fase dipende solo da quelle di cui legge i file). `runner/project.py` rilegge i manifesti e dice quali fasi sono concluse e quale viene dopo.
5. **Una fase.** Ogni fase è una classe in `steps/` che eredita da `PipelineStep` (`steps/base.py`). La base fa sempre le stesse quattro cose: verifica che le fasi a monte siano concluse; chiama `calcola()`, l'unica parte scritta dalla fase; controlla che i file dichiarati esistano e siano integri; scrive il manifesto.
6. **Il passaggio a R.** Dentro `calcola()` la fase prepara una richiesta e la passa a `rbridge/runner.py`: Python scrive un JSON di richiesta, lancia `Rscript` sullo script della fase, lo script scrive i suoi file e per ultimo un JSON con l'esito (successo o codice d'errore del catalogo).
7. **Errori e tentativi.** Un errore R diventa un'eccezione Python con il codice del catalogo. `runner/retry.py` ritenta solo i codici in `retry.whitelist` (che può solo restringere l'elenco ammesso dal catalogo), entro `retry.max_attempts`, applicando l'azione correttiva che la fase dichiara. Altrimenti l'esecuzione si ferma e scrive il punto di ripresa.
8. **La fine.** S14 scrive in `12_final/` l'oggetto `ps_final.rds`, le tabelle di testo e `checksum.sha256`. Il comando `report` (`report/builder.py`) legge la cartella e scrive il report HTML, senza ricalcolare nulla.
9. **Le analisi ecologiche.** `scripts/esegui.py eco run …` passa a `amplicon16s_eco`: copia l'oggetto finale (l'originale non si tocca), lancia gli script di `src/amplicon16s_eco/R/` e scrive le uscite (elenco e impronte in `manifest.json`, sintesi in `riepilogo.json`) solo ad analisi conclusa.

**Gli altri comandi.** `run` e `resume` usano lo stesso esecutore: `run` vuole una cartella di uscita vuota, `resume` rivaluta lo stato e rifà solo ciò che serve. `validate` si ferma dopo S0. `status` legge soltanto `99_logs/` e i manifesti, quindi si può lanciare durante un'esecuzione.

**Dove guardare quando qualcosa va storto.** In `99_logs/`: `pipeline.jsonl` (il log), `avvii.jsonl` (un evento per avvio), `punto_di_ripresa.json` e `.txt` (dopo un arresto, con il codice e che cosa fare), `stato_esecuzione.json` (durante l'esecuzione). Ogni fase scrive inoltre `letture_<passo>.tsv` nella sua cartella: è il tracciamento delle letture, ricomposto da `runner/tracciamento.py`.

## 5. Le fasi, in ordine di esecuzione

| Fase | Che cosa fa | Python (`steps/`) | R (`R/`) | Cartella |
|---|---|---|---|---|
| S0 | controlli sugli ingressi (gate G01-G15), inventario | `s00_validate.py` | `00_classificatore.R` (solo con IdTaxa) | `01_input_validation` |
| S1 | profilo di qualità e lunghezza | `s01_profile.py` | `01_profile.R` | `02_qc_profiles` |
| S2 | taglio del primer, filtro e troncamento | `s02_filter.py` | `02_filter.R` | `03_filtered` |
| S3 | modello d'errore per corsa | `s03_learn_errors.py` | `03_learn_errors.R` | `04_error_models` |
| S4 | inferenza delle varianti (dada2) | `s04_dada.py` | `04_dada.R` | `05_asv_inference` |
| SM | fusione delle coppie (solo appaiate) | `s04m_merge.py` | `04m_merge.R` | `05_pair_merge` |
| S5 | tabella delle sequenze | `s05_seqtab.py` | `05_seqtab.R` | `06_seqtab` |
| S6, S7 | chimere, filtro di lunghezza | `s06_chimera.py`, `s07_asv_length.py` | `06_chimera.R`, `07_asv_length.R` | `07_chimera` |
| S8 | tassonomia | `s08_taxonomy.py` | `08_taxonomy.R` | `08_taxonomy` |
| S10 | oggetto integrato (phyloseq) | `s10_phyloseq.py` | `10_phyloseq.R` | `10_phyloseq` |
| S11 | controlli positivi, soglia di profondità | `s11_controls.py` | `11_controls.R` | `11_controls` |
| S12 | decontaminazione (decontam) | `s12_decontam.py` | `12_decontam.R` | `11_controls` |
| S13 | filtri finali | `s13_filtri.py` | `13_filtri.R` | `12_final/intermedi` |
| S9 | albero filogenetico (facoltativo), dopo S13 | `s09_phylogeny.py` | `09_phylogeny.R` | `09_phylogeny` |
| S14 | oggetto finale ed export | `s14_finale.py` | `14_finale.R` | `12_final` |

Ogni file Python e R comincia con un commento che dice che cosa fa la fase, che cosa scrive e perché: leggilo prima del codice. Le funzioni R condivise (scrittura atomica, esito, lettura delle letture, KatharoSeq, decontaminazione) stanno in `R/lib/`.

## 6. I concetti da avere chiari

- **Fase conclusa.** Una fase è conclusa se le fasi a monte lo sono, il suo manifesto esiste, i suoi file corrispondono alle impronte e `calcolata_su` (parametri dichiarati, fasi a monte, versione) coincide con lo stato attuale. Altrimenti `resume` la rifà, con tutto ciò che ne dipende.
- **Parametri dichiarati.** Ogni fase elenca in `parametri` quelli da cui dipende (più gli elenchi per modalità: `parametri_appaiata`, `parametri_taglio_cercato`, `parametri_bayesiano`, `parametri_idtaxa`) e li legge attraverso una `VistaConfig` (`config/vista.py`): leggerne uno non dichiarato è un errore. È ciò che garantisce che cambiare un parametro rifaccia esattamente le fasi giuste.
- **Versione e registro del sorgente.** Ogni fase ha una `versione`; se cambi ciò che una fase calcola, la incrementi. Il registro `steps/registro_sorgente.json` tiene versione e impronta del sorgente di ogni fase (il suo modulo e i moduli che importa, lo script R, i file di `R/lib` che carica), e un test fallisce se non corrisponde al codice: toccare un modulo condiviso sposta il registro di più fasi. Dopo una modifica si aggiorna con `python3 scripts/registro_sorgente.py --aggiorna`; se la modifica non cambia i risultati (commenti, riorganizzazione) si dichiara con `--aggiorna --senza-effetto S4`, senza toccare la versione. **Lo script deve leggere il codice del clone**: fuori dall'immagine `PYTHONPATH=src python3 scripts/registro_sorgente.py --aggiorna` (o con il pacchetto installato con `pip install -e .`); nell'immagine con `-e PYTHONPATH=/app/src -e AMPLICON16S_R_DIR=/app/R`, come nel comando dei test del README. Senza, nell'immagine aggiornerebbe in silenzio il registro della copia installata, non quello del clone.
- **Catalogo degli errori.** Ogni codice (`E-S2-03`, `E-G15-10`…) sta in `errors/catalog.py` con il messaggio operativo (che cosa fare) e la categoria di gestione: fermarsi e far decidere l'operatore, ritentare, oppure proseguire con un ripiego registrato. I codici delle analisi ecologiche (`E-ECO-*`) stanno in `amplicon16s_eco/catalogo.py`.
- **Riproducibilità.** A parità di configurazione, codice e immagine, gli artefatti delle fasi sono identici byte per byte (manifesti, log, report e configurazione registrata no: portano date e percorsi).

## 7. Le modifiche più comuni

**Aggiungere un parametro.** Sei passi, in quest'ordine:

1. `config/schema.py`: tipo, vincoli e `description` (è il testo che il comando `schema` espone ai programmi esterni);
2. `config/defaults.py`: il predefinito con la fonte, se è lo standard del metodo; altrimenti l'elenco `OBBLIGATORI` o `CONDIZIONALI`;
3. `init/comando.py`, funzione `prepara`: la riga che `init` scrive (`metti(...)`). Un `assert` verso la fine della funzione verifica che ci sia una riga per ogni obbligatorio e condizionale, e senza `init` si ferma; i parametri della modalità appaiata hanno un controllo proprio (`PARAMETRI_DELLE_COPPIE`) e, se condizionali, vanno negli elenchi `CONDIZIONALI_DELLE_COPPIE` o `CONDIZIONALI_DI_MODALITA` di `defaults.py`;
4. `config/config.example.yaml` e `RIFERIMENTO.md`: le sezioni 6.1-6.7 per un obbligatorio o un condizionale (un test lo verifica per gli obbligatori), la tabella della sezione 6.8 per un parametro con un predefinito (con la fonte o il motivo); per le analisi ecologiche la sezione 4 e, se ha un predefinito, la 4.4; se è una scelta dell'utente con un valore di partenza, anche la tabella del passo 5 di `GUIDA.md`; se è una soglia che ferma o avvisa, anche `soglie/tabella.py`;
5. la dichiarazione nella fase che lo usa (`parametri` o l'elenco della sua modalità); se non incide sui risultati (per esempio i processori), `PARAMETRI_SENZA_EFFETTO` in `config/resolve.py`;
6. il test.

**Cambiare il calcolo di una fase.** Modifichi lo script R o la classe Python; incrementi `versione`; aggiorni il registro del sorgente; esegui i test della fase e il confronto sul dataset di riferimento (sezione 8) per vedere quali file cambiano.

**Aggiungere un codice d'errore.** Va in `errors/catalog.py` e nell'elenco dei codici attesi di `tests/test_w05_errors_catalog.py` (un test respinge i codici non elencati); poi il controllo nella fase o nel gate, e il test che verifica che il codice esca nel caso previsto.

**Aggiungere una fase.** `Passo` e `Nodo` in `runner/graph.py`; la cartella in `Fase` (`io_layer/artifacts.py`); la classe in `steps/`, registrata in `passi_realizzati()` (`runner/project.py`); lo script R; i codici nel catalogo e nel test; il registro del sorgente; `ORDINE` e `CARTELLE` in `dati/osd734/confronta_risultati.py` (un test li tiene allineati al grafo).

**Cambiare le analisi ecologiche.** Tutto sta in `src/amplicon16s_eco/`: parametri in `config.py` (e in `config/eco.example.yaml`, `RIFERIMENTO.md` sezione 4; un file di uscita nuovo va nella tabella della 4.3), comandi in `cli.py`, codici in `catalogo.py`, calcoli in `R/` (`valida.R`, `analisi.R`, `descrivi.R` e `R/lib/`). Test in `tests/test_eco_*.py`, dati in `tests/fixtures/eco/`. Le uscite della pipeline non si toccano.

**Aggiornare l'ambiente** (un pacchetto R, la data dell'istantanea CRAN, l'immagine di base). Otto passi:

1. modifichi `container/install_r_packages.R`, oppure nel `Dockerfile` `CRAN_SNAPSHOT` (la data dell'istantanea dei pacchetti CRAN) o il digest dell'immagine di base Bioconductor;
2. sposti `renv.lock` fuori dal clone: la build confronta i pacchetti con il file e fallirebbe, mentre senza file la verifica si salta;
3. costruisci: `docker build -f container/Dockerfile -t amplicon16s:nuova .`;
4. `scripts/genera_renv_lock.sh amplicon16s:nuova` rigenera `renv.lock` dall'immagine appena costruita;
5. ricostruisci l'immagine, che ora si verifica contro il file nuovo;
6. esegui la suite nell'immagine e il confronto sul riferimento (sezione 8); se i risultati cambiano per scelta, `confronta_risultati.py --aggiorna` riscrive le impronte attese;
7. pubblichi l'immagine nel registro (`docker push`) e ne annoti il digest;
8. sostituisci il digest vecchio con il nuovo ovunque compaia: README, GUIDA, `dati/osd734/README.md` e la chiave `run.container` delle configurazioni in `dati/` (`grep -rl <digest vecchio> .` li elenca).

## 8. Test e verifiche

- **La suite:** `pytest` (sezione 2) gira su 4 processi. I test che richiedono il dataset completo sono marcati `dati_reali` e `riferimento` e si attivano con `AMPLICON16S_CONFIG_DATI_REALI=<configurazione>`; quelli lunghi sono marcati `lenta` e si attivano con `AMPLICON16S_PROVE_LENTE=1`. I moduli di test per area sono nella tabella qui sotto; per un modulo del codice si cercano i test che lo nominano (`grep -l s04_dada tests/`).
- **I dati di prova** (`tests/fixtures/`): `osd734/ridotto/` (versione piccola di OSD-734, ricostruibile con `scripts/build_test_subset.py`, su cui gira la catena in CI), `osd276/` (letture appaiate), `osd170/`, `riferimento_sintetico/`, `idtaxa/`, `eco/`.
- **La CI** (`.github/workflows/ci.yml`): un job rapido con la suite Python e R con il solo `jsonlite` (senza Bioconductor, quindi i test di calcolo si saltano) e un job dentro l'immagine, diviso in due parti; un nuovo modulo di test lungo va aggiunto a `MODULI_DI_CALCOLO`.
- **Il confronto con il riferimento:** `dati/osd734/confronta_risultati.py` confronta un'esecuzione completa di OSD-734 con le impronte pubblicate (`checksum_artefatti.tsv`, `checksum_finali.sha256`) e, se qualcosa differisce, indica la prima fase in cui compare la differenza. Servono i dati completi (README, «Riprodurre il riferimento») e un clone senza modifiche non salvate in un commit, perché la configurazione di riferimento usa `strict_provenance`. Se un cambiamento è voluto, `--aggiorna` riscrive le impronte attese.

**I moduli di test per area:**

| Area | Moduli |
|---|---|
| configurazione e schema | `test_w03_w04_config_schema`, `test_config_condizionali`, `test_dipendenze_parametri`, `test_w26_congelamento` |
| ingressi, metadati e S0 | `test_w07_gates`, `test_w06_crosswalk`, `test_w06_join_restricted`, `test_w27_ingressi`, `test_righe_e_file`, `test_modalita_appaiata` |
| esecutore, grafo, ripresa, ponte R, log | `test_w05_artifacts`, `test_w05_logging`, `test_w05_errors_catalog`, `test_w08_rbridge`, `test_w08_w09_graph_resume`, `test_w09_w10_retry_policy`, `test_w16_recupero` |
| S1-S7 | `test_w11_s01_profile`, `test_s01_orientamento`, `test_w12_s02_filter`, `test_taglio_primer`, `test_w13_s03_batch`, `test_w14_s04_s07_denoising`, `test_w15_consolidamento` |
| S8, S10 | `test_w17_tassonomia`, `test_idtaxa`, `test_w18_oggetto` |
| S11, S12 | `test_w19_controlli`, `test_w29_soglia`, `test_w20_decontam`, `test_controlli_assenti` |
| S13, S14, S9, report | `test_w21_finale`, `test_w23_filogenesi`, `test_w22_report`, `test_w28_calcolo` |
| catene su sottoinsiemi | `test_catena_appaiata`, `test_catena_osd170_ridotto` |
| `init`, `soglie`, JSON | `test_init`, `test_soglie`, `test_interfaccia_json` |
| analisi ecologiche | `test_eco_configurazione`, `test_eco_diversita`, `test_eco_ordinazione`, `test_eco_statistica` |
| dati e script | `test_dati_dei_dataset`, `test_w21_dati_osd734`, `test_scarica_osdr`, `test_letture_forward`, `test_comunita_nota`, `test_confronto_genelab` |
| ambiente, riproducibilità, dataset completo | `test_w01_w02_import`, `test_w25_immagine`, `test_w30_ambiente`, `test_w30_correzioni`, `test_w30_indipendenza`, `test_w24_riproduzione`, `test_w07_dati_reali` |

## 9. Le regole che non si rompono

1. Nessun valore di un dataset nel codice o nei test: va in configurazione o in `dati/`.
2. Una fase scrive solo nella sua cartella e non modifica mai i file di un'altra.
3. Chi cambia un risultato cambia la `versione` della fase e lo verifica sul riferimento.
4. Le dipendenze (Python in `pyproject.toml`, R in `renv.lock`) sono a versione esatta.
5. Ogni errore ha un codice del catalogo e un messaggio che dice che cosa fare.

## 10. Dove trovare il resto

- `README.md`: che cos'è lo strumento, requisiti, avvio rapido, uscite, test.
- `GUIDA.md`: come si usa, passo per passo.
- `RIFERIMENTO.md`: ogni parametro (sezione 6: obbligatori e condizionali nelle 6.1-6.7, quelli con un predefinito nella 6.8), le analisi ecologiche (sezione 4, con le uscite nella 4.3 e i parametri con un predefinito nella 4.4), i formati JSON (sezione 5), i casi particolari, i limiti noti (sezione 7) e il glossario dei termini usati qui e nel codice (sezione 8: gate, manifesto, artefatto, degradazione, digest, lotto...).
- `RAZIONALE.md`: il perché scientifico di ogni metodo e, per ciascun risultato, i parametri e i file da toccare per cambiarlo.
- I commenti in testa a ogni file: il perché di ogni scelta.

## Appendice. Indice dei file

Una riga per file; il commento in testa a ciascuno dice il resto. I file `__init__.py` dichiarano il pacchetto e ne descrivono il contenuto, senza codice proprio (quello di `amplicon16s` porta la versione).

**`src/amplicon16s/`**

| File | Che cosa fa |
|---|---|
| `cli.py` | la riga di comando: sottocomandi, codici di uscita, blocco della cartella di uscita |
| `interfaccia.py` | i documenti JSON di `schema`, `validate --json` e `status --json`, e lo stato dell'esecuzione in corso |

**`config/`**

| File | Che cosa fa |
|---|---|
| `schema.py` | i gruppi e i parametri (pydantic), i loro vincoli e le descrizioni; caricamento del file YAML |
| `defaults.py` | i predefiniti con la fonte, e gli elenchi di obbligatori, condizionali e parametri di modalità |
| `resolve.py` | i parametri derivati, il digest della configurazione e la sua registrazione in `00_config/` |
| `vista.py` | la configurazione come la vede una fase: solo i parametri che dichiara |

**`runner/`**

| File | Che cosa fa |
|---|---|
| `graph.py` | il grafo: fasi, dipendenze di dato, fasi attive secondo la configurazione |
| `project.py` | lo stato dell'esecuzione riletto dai manifesti: fasi concluse, da eseguire, da rifare |
| `executor.py` | percorre il grafo: controlli di avvio, esecuzione delle fasi, tentativi, punto di ripresa |
| `retry.py` | la politica dei tentativi ripetuti e le azioni correttive |
| `provenienza.py` | versione e impronta del sorgente di ogni fase; la regola rigorosa (git, ambiente R) |
| `tracciamento.py` | la tabella delle letture per campione, ricomposta dalle tabelle delle fasi |

**`steps/`**

| File | Che cosa fa |
|---|---|
| `base.py` | `PipelineStep`: lo scheletro comune (fasi a monte, `calcola()`, verifica degli artefatti, manifesto) |
| `s00_validate.py` … `s14_finale.py` | una fase per file (tabella della sezione 5): parametri dichiarati, richiesta allo script R, controlli sull'esito |
| `registro_sorgente.json` | versione e impronta del sorgente di ogni fase; lo aggiorna solo `scripts/registro_sorgente.py` |

**`rbridge/`, `gates/`, `errors/`, `logging/`**

| File | Che cosa fa |
|---|---|
| `rbridge/payload.py` | i due file JSON, richiesta ed esito, che passano fra Python e R |
| `rbridge/runner.py` | l'unico punto che lancia `Rscript`: processo isolato, tempo massimo, lettura dell'esito |
| `gates/g01_g15.py` | i quindici controlli della validazione, e i controlli di avvio G12, G14, G15 |
| `gates/registry.py` | l'ordine dei gate e la loro esecuzione in sequenza |
| `errors/catalog.py` | il catalogo dei codici: messaggio, che cosa fare, categoria di gestione |
| `errors/exceptions.py` | le eccezioni della pipeline, legate ai codici |
| `logging/logger.py` | il log: console per l'operatore, file strutturato in `99_logs/` |

**`metadata/`, `io_layer/`**

| File | Che cosa fa |
|---|---|
| `metadata/tabelle.py` | lettura delle tabelle separate da tabulazioni |
| `metadata/crosswalk.py` | l'inventario: file, righe dei metadati e file del lotto messi in corrispondenza |
| `metadata/controls_map.py` | la classe di ogni campione (biologico, controllo negativo, positivo) dalle etichette dichiarate |
| `metadata/ampliconi.py` | il confronto fra ciò che una riga dichiara (amplicone, primer) e il primer del file collegato |
| `metadata/models.py` | i modelli dell'inventario (`Campione`, `Inventario`) |
| `metadata/lettura_inventario.py` | la lettura dell'inventario scritto da S0, per le fasi successive |
| `io_layer/artifacts.py` | l'albero delle cartelle di uscita, i manifesti, la scrittura atomica |
| `io_layer/checksums.py` | calcolo e verifica degli SHA-256 |
| `io_layer/conteggi.py` | lettura delle tabelle del tracciamento |
| `io_layer/reads.py` | lettura dei FASTQ: ispezione delle prime letture, coppie e marcatori, verso delle letture, ricerca del primer |

**`init/`, `soglie/`, `report/`**

| File | Che cosa fa |
|---|---|
| `init/comando.py` | il sottocomando `init`: la funzione `prepara` costruisce una voce per parametro |
| `init/rilevamento.py` | ciò che `init` ricava dai file: fatti verificati ed euristiche |
| `init/vocabolario.py` | le parole e i primer con cui `init` riconosce colonne, controlli e regioni |
| `init/scrittura.py` | il file YAML che `init` scrive, con marcatori e commenti |
| `soglie/tabella.py` | le soglie che fermano o avvisano: direzione, limite, margine, controlli spenti |
| `soglie/misure.py` | la misura di ogni controllo, letta dagli artefatti di un'esecuzione conclusa |
| `soglie/comando.py` | il sottocomando `soglie`: proposte, file riscritto, documento JSON |
| `report/builder.py` | il report HTML e le sue tabelle, dagli artefatti verificati sui manifesti |
| `report/templates/` | il modello HTML e lo stile del report |

**`R/`**

| File | Che cosa fa |
|---|---|
| `00_ambiente.R` | non è una fase: le versioni di R e dei pacchetti, per la regola rigorosa e per S8 |
| `00_classificatore.R` | non è una fase: verifica del classificatore addestrato di IdTaxa, per S0 |
| `01_profile.R` … `14_finale.R` | uno script per fase (tabella della sezione 5); l'intestazione elenca i file che scrive |
| `lib/errors.R` | l'esito verso Python: codice del catalogo, uscita 3 per un errore dichiarato |
| `lib/io_json.R` | lettura della richiesta, scritture atomiche di JSON e RDS |
| `lib/letture.R` | la tabella delle letture per campione che ogni fase registra |
| `lib/risorse.R` | la memoria di picco del processo, per il log |
| `lib/primer.R` | ricerca e taglio del primer lungo la lettura (S1, S2) |
| `lib/errore_loess.R` | le funzioni di stima dei tassi d'errore, compresa la variante monotona (S3, S4) |
| `lib/idtaxa.R` | la classificazione con `DECIPHER::IdTaxa` (S8, `00_classificatore.R`) |
| `lib/oggetto.R` | le verifiche dell'oggetto integrato sugli identificativi (S10) |
| `lib/katharoseq.R` | la curva di calibrazione e la conformità dei controlli positivi (S11) |
| `lib/decontaminazione.R` | la decontaminazione per prevalenza, aggregata e per piastra (S12) |
| `lib/albero.R` | scrittura, lettura e verifica dell'albero (S9, S14) |
| `lib/export.R` | la ricostruzione dell'oggetto dagli export di testo (S14) |

**`src/amplicon16s_eco/`**

| File | Che cosa fa |
|---|---|
| `cli.py`, `__main__.py` | la riga di comando `eco`: `init`, `validate`, `run`, `status`, `schema` |
| `config.py` | lo schema dei parametri, il caricamento del file YAML e il digest |
| `catalogo.py` | i codici `E-ECO-`: rifiuti e avvisi |
| `esecuzione.py` | l'orchestrazione: copia dell'oggetto, validazione, calcolo, manifesto |
| `ponte.py` | il lancio degli script R del pacchetto attraverso il ponte della pipeline |
| `init.py` | `eco init`: proposte e avvisi dalla descrizione dell'oggetto |
| `interfaccia.py` | i documenti JSON di `eco schema`, `eco validate --json` ed `eco status --json` |
| `R/descrivi.R` | la descrizione dell'oggetto (colonne, ranghi, albero) per `eco init` |
| `R/valida.R` | la validazione dell'analisi sull'oggetto, prima di ogni calcolo |
| `R/analisi.R` | le analisi, nell'ordine: campioni, alfa, composizione, distanze, ordinazione, test |
| `R/lib/uscite.R` | la scrittura di tabelle e grafici, identici byte per byte |
| `R/lib/campioni.R` | lettura dell'oggetto, campioni di ogni analisi, esame delle variabili |
| `R/lib/alfa.R`, `composizione.R`, `beta.R`, `ordinazione.R` | rarefazione e indici; agglomerazione e proporzioni; distanze; PCoA e NMDS |
| `R/lib/test.R` | Wilcoxon e Kruskal-Wallis, PERMANOVA e test delle dispersioni |
| `R/lib/continua.R`, `ripetute.R` | le variabili continue; i test entro soggetto |

Il flusso di `eco run`:

1. `cli.py` legge gli argomenti, carica la configurazione con `config.py` e chiama `esecuzione.py`;
2. `esecuzione.py` copia l'oggetto in una cartella di lavoro, prepara la richiesta e lancia, attraverso `ponte.py`, prima `R/valida.R` e poi `R/analisi.R`;
3. `ponte.py` usa il ponte della pipeline (`amplicon16s.rbridge`): stesso contratto a file JSON, stessi processi R isolati;
4. `R/analisi.R` carica le librerie di `R/lib/` e scrive tabelle, grafici e `riepilogo.json` in una cartella provvisoria; solo ad analisi conclusa `esecuzione.py` vi aggiunge `configurazione.json` e `manifest.json` e sposta tutto nella cartella indicata.

**`scripts/`, `container/`**

| File | Che cosa fa |
|---|---|
| `scripts/esegui.py` | il punto d'ingresso: esegue `src/` e `R/` di questo clone, per la pipeline e per `eco` |
| `scripts/scarica_osdr.py` | scarica da OSDR metadati e letture 16S di uno studio |
| `scripts/registro_sorgente.py` | verifica e aggiorna il registro del sorgente delle fasi |
| `scripts/build_test_subset.py` | ricostruisce il sottoinsieme di prova dal dataset di riferimento |
| `scripts/costruisci_riferimento_sintetico.py` | costruisce il riferimento tassonomico sintetico dei test |
| `scripts/genera_renv_lock.sh` | rigenera `renv.lock` dall'immagine costruita |
| `container/Dockerfile` | l'immagine: basi ancorate per digest, pacchetti R, dada2, Python, verifica contro `renv.lock` |
| `container/install_r_packages.R` | l'installazione dei pacchetti R dall'istantanea datata |
| `container/snapshot_renv_lock.R` | la scrittura di `renv.lock` dalle versioni installate |
| `container/verify_renv_lock.R` | la verifica dei pacchetti installati contro `renv.lock` |
| `container/dada2/` | lo script e la correzione con cui dada2 è compilato nell'immagine |
