# Razionale scientifico di amplicon16s

Questo documento descrive **tutti i risultati** che lo strumento produce, nell'ordine in cui li produce. Per ciascuno dice:

- **Che cosa produce**: i file di uscita, nella cartella dell'esecuzione;
- **Dove nasce**: il codice che lo calcola;
- **Metodo**: lo strumento scientifico usato e la sua fonte;
- **Perché**: il razionale della scelta;
- **A che cosa serve**: che cosa permette di sapere o di fare a chi usa il risultato;
- **Per cambiarlo**: i parametri della configurazione, se basta una scelta diversa, e i file su cui intervenire, se cambia il metodo.

Serve a chi deve valutare la validità scientifica dello strumento e a chi deve cambiarne un metodo. I percorsi della sezione 3 sono relativi a `src/amplicon16s/` per il Python e alla radice per `R/`; quelli della sezione 4 a `src/amplicon16s_eco/`. Il significato di ogni parametro è in `RIFERIMENTO.md` (sezioni 4 e 6) e in `config/config.example.yaml`; i termini (gate, manifesto, artefatto, degradazione, fedeltà, lotto, piastra...) sono nel glossario di `RIFERIMENTO.md`, sezione 8; come si lavora sul codice è in `SVILUPPO.md`.

**Indice.** 1. Il problema scientifico e il flusso · 2. I principi e le scelte proprie dello strumento · 3. Dalle letture all'oggetto finale (3.1 S0 … 3.17 risultati trasversali) · 4. Le analisi ecologiche (4.1-4.7) · 5. Limiti dichiarati e procedure di verifica · 6. Fonti

## 1. Il problema scientifico e il flusso

**Il gene 16S.** Il gene dell'RNA ribosomiale 16S è presente in tutti i batteri e gli archei. Alterna regioni conservate, molto simili in quasi tutte le specie, e regioni variabili (V1-V9), che cambiano da un gruppo all'altro. Si amplifica per PCR una regione variabile (per esempio V4) con due primer che si legano alle regioni conservate ai suoi lati, e si sequenziano i prodotti: ogni lettura è un frammento del 16S di un organismo presente nel campione, o portato dal laboratorio (vedi sotto). Dalle letture si vuole sapere **quali organismi ci sono** (la tassonomia) e **in che proporzioni** (una tabella di conteggi per campione), per confrontare i campioni fra loro.

**Le fonti di errore.** Fra il campione e la tabella ci sono diverse fonti di errore; lo strumento ne tratta tre, una per una:

- il **sequenziamento** introduce sostituzioni di basi: una sequenza vera genera molte sequenze quasi uguali. L'inferenza delle varianti (S3-S4) le separa;
- la **PCR** produce chimere, sequenze ibride di due sequenze genitrici (S6);
- la **contaminazione** da reagenti e ambiente di laboratorio aggiunge DNA che non viene dal campione. Nei campioni a **bassa biomassa**, con poco DNA di partenza, può essere una parte grande delle letture (Salter et al., 2014). Per questo i controlli sono centrali: i **controlli negativi** (bianchi, senza campione) dicono che cosa porta il laboratorio, e servono alla decontaminazione (S12); i **controlli positivi** (un organismo noto, in diluizioni note) dicono da quale profondità di sequenziamento la quota di letture del DNA del campione raggiunge la sensibilità dichiarata (S11).

**Varianti esatte, non raggruppamenti.** I metodi più vecchi raggruppavano le sequenze simili al 97% in unità operative (OTU). Lo strumento usa invece **varianti di sequenza esatte** (ASV): sequenze risolte fino alla singola base, distinte dagli errori con un modello statistico. Sono riproducibili e confrontabili fra studi che usano la stessa regione, con gli stessi primer e lo stesso troncamento, e non dipendono da una soglia di similarità arbitraria (Callahan et al., 2017).

**Le proporzioni sono relative.** Il numero di letture di un campione dipende dal sequenziamento, non dalla quantità di DNA: i conteggi stimano la quota relativa delle copie del 16S amplificate di ogni taxon, non la sua quantità assoluta (Gloor et al., 2017). Quella quota è distorta anche dal numero di copie del gene per genoma e dall'efficienza diversa di estrazione e primer fra taxa (Kembel et al., 2012; McLaren et al., 2019).

**Il flusso, in cinque blocchi:**

1. **Controllo degli ingressi** (S0, S1): i dati sono quelli dichiarati e sono utilizzabili?
2. **Dalle letture alle varianti** (S2-S7): filtro, modello d'errore, inferenza, fusione delle coppie, chimere, lunghezza.
3. **Nome e assemblaggio** (S8, S10): la tassonomia di ogni variante e l'oggetto con conteggi, tassonomia e metadati.
4. **Controlli e filtri** (S11, S12, S13): la soglia di profondità dai controlli positivi, la decontaminazione dai negativi, i filtri finali.
5. **Albero e consegna** (S9, S14): l'albero filogenetico delle varianti finali e l'oggetto finale con i suoi export. Poi, separate, le analisi ecologiche.

**Perché quest'ordine.** La tassonomia viene dopo l'inferenza perché classifica le sequenze esatte, una sola volta per variante. I controlli vengono dopo l'oggetto integrato perché richiedono la tassonomia (la quota del taxon atteso nei positivi) e i metadati (classe e piastra di ogni campione). La decontaminazione precede i filtri perché usa i controlli negativi, che l'oggetto finale non contiene, e l'insieme completo delle varianti. L'albero viene dopo i filtri perché serve solo alle varianti consegnate.

## 2. I principi e le scelte proprie dello strumento

- **Metodi pubblicati.** I calcoli principali sono affidati a pacchetti R pubblicati (dada2, DECIPHER, decontam, phyloseq, phangorn, ape, vegan, stats). Le parti scritte dallo strumento (ricerca del primer, stima dell'errore per qualità raggruppate, curva KatharoSeq, pseudo-pooling a lotti) dichiarano nel commento la regola e la fonte, e così ogni scostamento dal predefinito di un pacchetto.
- **Nessuna scelta ereditata in silenzio.** Un parametro ha un predefinito solo se è lo standard del metodo (con la fonte in `config/defaults.py`), l'unica forma realizzata o un parametro dello strumento (risorse, seme, formato delle uscite). Tutto ciò che dipende dai dati o dallo studio (regione, primer, lunghezza delle letture, controlli, soglie) va dichiarato; `init` ne propone una parte, da accettare esplicitamente.
- **Un controllo che fallisce ferma o dichiara, non nasconde.** Ogni condizione anomala ha un codice nel catalogo (`errors/catalog.py`): l'esecuzione si ferma, oppure prosegue e registra la degradazione nel manifesto della fase e nel report.
- **Le classi dei campioni si dichiarano.** Biologici, controlli negativi e positivi si leggono da una colonna dei metadati con etichette dichiarate (`metadata/controls_map.py`), mai dal nome del campione. I controlli di qualità si applicano alle classi per cui hanno senso: un controllo negativo con poche letture è un bianco pulito, non un difetto.
- **Riproducibilità.** Ogni passo casuale usa il seme dichiarato, il numero di processi non cambia i risultati, e a parità di dati, configurazione e immagine gli artefatti delle fasi sono identici byte per byte.

**Le scelte proprie dello strumento.** Dove lo strumento non applica un metodo pubblicato così com'è, ma lo adatta o aggiunge una regola, la scelta è qui, con il perché e la sezione in cui è descritta. Sono i punti su cui un revisore deve esercitare il proprio giudizio.

| Scelta | Perché | Dove |
|---|---|---|
| Ricerca del primer con sole sostituzioni e discordanze dichiarate | posizione del taglio non ambigua | 3.3 |
| Troncamento proposto da `init`: qualità mediana 30, non più del 5% di letture scartate per classe; per le coppie, allungamento fino a qualità 20 se serve alla fusione | regola scritta e uguale per ogni dataset; il valore resta una scelta dell'utente | 3.3 |
| Stima dell'errore monotona per le qualità raggruppate | evitare tassi d'errore che crescono con la qualità | 3.4 |
| Pseudo-pooling in due passate esplicite, a lotti, con le informazioni a priori raccolte su tutte le corse | memoria indipendente dal numero di campioni; il modello d'errore resta per corsa | 3.5 |
| Identificativi delle varianti per abbondanza nei soli campioni biologici | il ceppo dei controlli positivi non va in testa alla numerazione | 3.11 |
| KatharoSeq con curva aggregata e curve per piastra sugli stessi punti, scelta per AIC, ripieghi dichiarati | una soglia per piastra quando i dati la sostengono, confrontabile fra piastre | 3.12 |
| Conformità dei controlli positivi con lo z modificato, per livello di diluizione | un controllo si giudica rispetto a quelli della stessa concentrazione | 3.12 |
| Soglia di profondità misurata e applicata sulle letture senza chimere | una sola grandezza per stima e filtro | 3.12, 3.14 |
| Decontam calcolato sempre nelle due modalità; per piastra, combinazione riprodotta, e una piastra con pochi negativi confrontata con i negativi di tutte le piastre | la modalità si dichiara prima; l'altra resta come diagnostica | 3.13 |
| Albero con ricerca deterministica e radice al punto medio | lo stesso albero a ogni esecuzione; radice indipendente da scelte casuali | 3.15 |
| Analisi ecologiche: rarefazione senza reinserimento, categoria dei non assegnati, Jaccard in presenza/assenza, p non sotto 1/D con poche disposizioni | evitare scostamenti silenziosi dai predefiniti dei pacchetti | 4.2-4.6 |

## 3. Dalle letture all'oggetto finale

### 3.1 Validazione degli ingressi e inventario (S0)

- **Che cosa produce.** In `01_input_validation/`: `gates.json` (l'esito dei 15 controlli G01-G15), `crosswalk.tsv` e `inventario.json` (ogni file abbinato al suo campione e alla sua classe), `letture_ispezionate.tsv` (e `_reverse` per le appaiate: primer e motivo conservato nelle prime letture di ogni file); con IdTaxa `classificatore_addestrato.json`.
- **Dove nasce.** `steps/s00_validate.py`, `gates/g01_g15.py`, `gates/registry.py`, `metadata/crosswalk.py`, `metadata/controls_map.py`, `io_layer/reads.py` (riconoscimento delle coppie dai nomi dei file, comprese le forme rinominate da GeneLab, `_R1_raw`/`_R2_raw`); `R/00_classificatore.R` solo con IdTaxa.
- **Metodo.** Controlli deterministici su file e metadati: leggibilità, corrispondenza uno a uno fra file e campioni, classe di ogni campione (G11), integrità degli archivi (G13), formato single-end o appaiato (G07), coerenza del lotto (G08), compatibilità del troncamento con le lunghezze (G09), primer e segnale della regione (G10: il primer non deve essere in testa alle letture quando il filtro non lo toglie, e il motivo conservato della regione dichiarata deve comparire nella mediana dei biologici e dei positivi), impronta del riferimento tassonomico (G12), risorse (G14). Le verifiche sulle sequenze si fanno sulle prime `qc.head_reads` letture.
- **Perché.** Un errore negli ingressi (un file corrotto, una colonna sbagliata, un troncamento incompatibile) produrrebbe risultati plausibili e sbagliati dopo ore di calcolo; qui si scopre in pochi minuti. Il valore che `init` propone per `qc.head_reads`, 9.604 letture, dà a una frazione un intervallo di confidenza al 95% di più o meno un punto percentuale se le letture fossero estratte a caso; le prime letture di un file non lo sono, e per questo S1 ripete sui file completi i controlli di lunghezza e di coppia. Un valore della colonna dei controlli che non ricade in alcuna etichetta ferma l'esecuzione: non viene attribuito a una classe di ripiego.
- **A che cosa serve.** `gates.json` dice, prima di qualunque calcolo, se i dati sono utilizzabili e, per ogni controllo fallito, il codice e il dettaglio; che cosa fare è nella voce del codice nel catalogo (`errors/catalog.py`), che il messaggio sulla console riporta. L'inventario è la corrispondenza file-campione-classe che tutte le fasi successive usano. `letture_ispezionate.tsv` permette di verificare che primer e regione siano quelli dichiarati.
- **Per cambiarlo.** `io.*`, `meta.*`, `ctrl.*` (colonna ed etichette delle classi, riclassificazione con `ctrl.blank_override_values`), `qc.head_reads`, `qc.primer_sequence`, `qc.conserved_motif`, `qc.max_primer_hit_frac`, `qc.min_motif_frac`, `tax.ref_md5` (o `tax.idtaxa_trained_md5`). Codice: il gate in `gates/g01_g15.py`; i nomi delle coppie in `io_layer/reads.py`.

### 3.2 Profilo di qualità e di lunghezza (S1)

- **Che cosa produce.** In `02_qc_profiles/`: `lunghezze.tsv`, `qualita.tsv` (qualità per posizione: media e percentili), `letture_grezze.tsv`, `valori_qualita.tsv`, `riepilogo.json`, `coppie.tsv`; per le appaiate gli stessi profili con il suffisso `_reverse`.
- **Dove nasce.** `steps/s01_profile.py`, `R/01_profile.R`, `io_layer/reads.py`.
- **Metodo.** Lettura completa di ogni file, a blocchi.
- **Perché.** Chiude i limiti dei controlli di S0, che guardano solo le prime letture: la quota di letture più corte del troncamento si misura per classe (biologici e positivi separatamente, perché in una media unica una classe piccola potrebbe perdere tutte le letture senza farsi notare); si riconoscono file con le due letture della coppia mescolate o letture inverse prese per forward; si segnalano le qualità raggruppate in pochi valori (tipiche di NovaSeq e NextSeq), che richiedono una stima diversa del modello d'errore (3.4).
- **A che cosa serve.** `qualita.tsv` e `lunghezze.tsv` sono la base su cui si sceglie il troncamento (`init` lo propone da qui, 3.3); `valori_qualita.tsv` mostra se le qualità sono raggruppate, da cui la scelta di `err.error_function`; `coppie.tsv` conta i segni di coppia di ogni file, per riconoscere file con le due letture mescolate.
- **Per cambiarlo.** `filter.truncLen`, `qc.max_frac_short_reads`, `filter.truncLen_shortfall_warn`. Codice: `R/01_profile.R`.

### 3.3 Taglio del primer, filtro e troncamento (S2)

- **Che cosa produce.** In `03_filtered/`: le letture filtrate `<accession>_filt.fastq.gz` (e `_filt_reverse` per le appaiate), `letture_prefiltro.tsv`, `letture_filtrate.tsv`; con il primer cercato lungo la lettura anche `taglio_primer.tsv`, `posizioni_primer.tsv`, `letture_taglio_primer.tsv`.
- **Dove nasce.** `steps/s02_filter.py`, `R/02_filter.R`, `R/lib/primer.R`.
- **Metodo.** `dada2::filterAndTrim` (Callahan et al., 2016): troncamento a lunghezza fissa e scarto delle letture con più di `maxEE` errori attesi, con basi indeterminate o di PhiX. Il primer si toglie a lunghezza fissa (`filter.trimLeft`) oppure si cerca lungo la lettura, quando il protocollo lo fa precedere da un numero variabile di basi: sole sostituzioni (nessuna inserzione, così la posizione del taglio non è ambigua), fino a `filter.primer_max_mismatch` discordanze (`init` propone un decimo della lunghezza del primer, la frazione d'errore predefinita di cutadapt, qui applicata alle sole sostituzioni).
- **Perché.** dada2 tratta lunghezze diverse come sequenze diverse: il troncamento fisso fa sì che la stessa variante dia la stessa sequenza e toglie le code di bassa qualità; le basi indeterminate non sono ammesse dall'inferenza; il primer va tolto perché la sua sequenza è sintetica, non del campione. `maxEE = 2` è il valore dei parametri standard del tutorial di dada2; `truncQ = 2`, `maxN = 0` e `rm.phix` sono i predefiniti del pacchetto. Il troncamento che `init` propone è il minore fra la lunghezza oltre la quale la qualità mediana scende sotto 30 (un errore atteso ogni mille basi) e la lunghezza che scarta non più del 5% delle letture di ciascuna classe (`RIFERIMENTO.md`, 6.3); per le letture appaiate, se la somma dei due troncamenti non basta alla fusione, `init` li allunga solo dentro posizioni con qualità mediana almeno 20, e altrimenti li lascia da scegliere. La perdita al filtro si controlla sui biologici e sui positivi, non sui negativi, che hanno poca biomassa e qualità spesso peggiore senza che i parametri siano sbagliati.
- **A che cosa serve.** Le letture filtrate sono l'ingresso dell'inferenza. `letture_filtrate.tsv` dice quanto ogni campione perde al filtro, cioè se i parametri sono adatti; `posizioni_primer.tsv` permette di verificare la ricerca del primer.
- **Per cambiarlo.** `filter.*` (`truncLen`, `trimLeft`, `maxEE`, `truncQ`, `primer_removal` e i parametri della ricerca), `qc.max_zeroed_samples`, `qc.max_frac_lost_filter`. Codice: `R/02_filter.R`, `R/lib/primer.R` (la stessa regola è in `io_layer/reads.py`).

### 3.4 Modello d'errore di sequenziamento (S3)

- **Che cosa produce.** In `04_error_models/`: `modello_<nome>.rds` e il grafico diagnostico `modello_<nome>.png` (un modello per corsa, o `tutti`), `convergenza.json`, `corrispondenza.tsv` (quale modello usa ogni campione), `modelli.json` (su quali campioni e quante basi è stata fatta la stima).
- **Dove nasce.** `steps/s03_learn_errors.py`, `R/03_learn_errors.R`, `R/lib/errore_loess.R`.
- **Metodo.** `dada2::learnErrors`: stima dei tassi di sostituzione in funzione della qualità, con l'algoritmo di auto-consistenza di dada2, su `err.nbases` basi (predefinito del pacchetto, 10^8). La curva si adatta con il loess standard di dada2 (`loessErrfun`) oppure, per le qualità raggruppate, con la variante `loess_monotono`: loess di primo grado con span 2 e pesi logaritmici, come nelle modifiche discusse dalla comunità di dada2 per i dati NovaSeq (issue 938 e 1307 del repository di dada2), più un vincolo di monotonia (tasso non crescente con la qualità) che è una scelta dello strumento, più forte di quello dell'issue.
- **Perché.** Un modello per corsa, perché corse diverse hanno profili d'errore diversi e un modello unico li medierebbe. I campioni su cui si stima li sceglie la fase, mescolati con il seme se `err.randomize` è vero, così la stima è riproducibile e documentata. Con le qualità raggruppate il loess standard può stimare un tasso d'errore che cresce con la qualità, che non ha senso fisico. Un modello che non converge ferma l'esecuzione: aumentare le basi cambierebbe i risultati, ed è una decisione dell'operatore.
- **A che cosa serve.** I modelli sono l'ingresso dell'inferenza. Il grafico permette di giudicare a occhio se la curva stimata segue i tassi osservati: è il controllo visivo raccomandato dagli autori di dada2. `convergenza.json` e `modelli.json` documentano la stima.
- **Per cambiarlo.** `err.batch_column` (la colonna della corsa), `err.error_function`, `err.nbases`, `err.max_consist`, `err.randomize`. Codice: `R/03_learn_errors.R`, `R/lib/errore_loess.R`.

### 3.5 Inferenza delle varianti di sequenza (S4)

- **Che cosa produce.** In `05_asv_inference/`: `varianti_per_campione.rds`, `varianti.tsv`, `letture_denoised.tsv`, `inferenza.json`; con il pseudo-pooling `priori.tsv` e i modelli d'errore della seconda passata; per le appaiate gli stessi per la lettura reverse e i file per la fusione.
- **Dove nasce.** `steps/s04_dada.py`, `R/04_dada.R`.
- **Metodo.** `dada2::dada` (Callahan et al., 2016): separa le sequenze vere dagli errori di sequenziamento con il modello d'errore della corsa del campione, e restituisce varianti di sequenza esatte. Con `dada.pool = pseudo` il pseudo-pooling è realizzato in due passate esplicite, a lotti, con la regola del sorgente di dada2 1.36.0 (informazioni a priori: le varianti presenti in almeno due campioni); la seconda passata usa, come fa dada2, il modello d'errore ristimato dalle transizioni della prima, separatamente per ciascun modello di S3. `OMEGA_A` è il predefinito del pacchetto.
- **Perché.** Il pseudo-pooling aumenta la sensibilità per le varianti rare presenti in più campioni senza il costo del pooling completo; le due passate esplicite permettono di lavorare su centinaia di campioni senza tenerli tutti in memoria, e il risultato non dipende dalla dimensione del lotto (verificato dai test). Le informazioni a priori si raccolgono su tutte le corse insieme, una scelta diversa da chiamare `dada(pool = "pseudo")` corsa per corsa: il modello d'errore descrive la macchina, la presenza di una variante in più campioni no. Le informazioni a priori abbassano solo la soglia per riconoscere una sequenza già presente nel campione, non ne aggiungono; una variante condivisa può comunque venire anche da contaminazione incrociata, che le fasi successive trattano (3.13).
- **A che cosa serve.** Le varianti sono le unità di tutte le analisi successive. `varianti.tsv` dice quante varianti ha ogni campione; `priori.tsv` elenca le sequenze usate come informazione a priori nella seconda passata, con prevalenza e abbondanza nella prima.
- **Per cambiarlo.** `dada.pool`, `dada.omega_a`. Codice: `R/04_dada.R`.

### 3.6 Fusione delle coppie (SM, solo letture appaiate)

- **Che cosa produce.** In `05_pair_merge/`: `varianti_per_campione.rds`, `fusione.tsv`, `fusione.json`, `lunghezze_fuse.tsv`, `letture_fuse.tsv`.
- **Dove nasce.** `steps/s04m_merge.py`, `R/04m_merge.R`.
- **Metodo.** `dada2::mergePairs`: la lettura forward e il complemento inverso della reverse si fondono se si sovrappongono per almeno `pair.min_overlap` basi (predefinito 12) con al più `pair.max_mismatch` discordanze (predefinito 0); le code oltre l'inizio dell'altra lettura si tolgono sempre (`trimOverhang`).
- **Perché.** Ricompone l'intera regione amplificata. Le code tolte sono, con letture più lunghe dell'amplicone, il complemento inverso del primer opposto e l'adattatore. Una quota bassa di coppie fuse ferma l'esecuzione: indica troncamenti troppo corti per la regione o file che non corrispondono, e proseguire consegnerebbe i soli campioni fortunati.
- **A che cosa serve.** `fusione.tsv` dà la quota di coppie fuse per campione, per giudicare i troncamenti; `lunghezze_fuse.tsv` permette di verificare che le sequenze fuse abbiano la lunghezza della regione attesa.
- **Per cambiarlo.** `pair.*` (sovrapposizione, discordanze, quote minime, lunghezza attesa). Codice: `R/04m_merge.R`.

### 3.7 Tabella delle sequenze (S5)

- **Che cosa produce.** In `06_seqtab/`: `tabella.rds` (campioni × varianti, conteggi interi), `letture_tabella.tsv`, `tabella.json`.
- **Dove nasce.** `steps/s05_seqtab.py`, `R/05_seqtab.R`.
- **Metodo.** `dada2::makeSequenceTable`.
- **Perché.** Il numero di varianti si controlla prima di allocare la tabella (`qc.max_asv_count`): la memoria cresce con campioni × varianti, e un numero anomalo di varianti segnala un problema a monte.
- **A che cosa serve.** È la tabella su cui lavorano tutte le fasi successive; `tabella.json` ne riporta dimensioni e distribuzione delle lunghezze.
- **Per cambiarlo.** `qc.max_asv_count`. Codice: `R/05_seqtab.R`.

### 3.8 Rimozione delle chimere (S6)

- **Che cosa produce.** In `07_chimera/`: `tabella_senza_chimere.rds`, `chimere.tsv` (ogni variante chimerica), `chimere_per_campione.tsv`, `chimere_per_gruppo.tsv`, `letture_senza_chimere.tsv`, `chimere.json`.
- **Dove nasce.** `steps/s06_chimera.py`, `R/06_chimera.R`.
- **Metodo.** `dada2::removeBimeraDenovo` con il metodo `chimera.method` (predefinito `consensus`): una variante è chimerica se si ricostruisce dall'unione di due varianti almeno `chimera.min_fold_parent_over_abundance` volte più abbondanti (parametro da dichiarare, perché il predefinito del pacchetto dipende da metodo e versione). Con `consensus` la decisione si prende campione per campione, e la variante si toglie se è chimerica in almeno il 90% dei campioni in cui compare; gli altri parametri sono i predefiniti di dada2.
- **Perché.** Le chimere sono artefatti della PCR, non organismi. La quota chimerica si misura sulle letture, non sulle varianti (le chimere sono tante varianti rare, e una chimera abbondante pesa più di cento chimere da una lettura), ed è giudicata sui biologici e sui positivi.
- **A che cosa serve.** `chimere.tsv` dice quali varianti sono state tolte e quante letture portavano; `chimere_per_campione.tsv` fa vedere i campioni con una PCR anomala.
- **Per cambiarlo.** `chimera.method`, `chimera.min_fold_parent_over_abundance`, gli altri `chimera.*`, `qc.warn_frac_chimeric`, `qc.stop_frac_chimeric`. Codice: `R/06_chimera.R`.

### 3.9 Filtro di lunghezza delle varianti (S7)

- **Che cosa produce.** In `07_chimera/`: `tabella_asv.rds`, `lunghezze.tsv`, `letture_lunghezza.tsv`.
- **Dove nasce.** `steps/s07_asv_length.py`, `R/07_asv_length.R`; l'intervallo è calcolato in `config/resolve.py`.
- **Metodo.** Si tengono le varianti con lunghezza fra `asv.len_min` e `asv.len_max`, ricavati da troncamento e taglio (single-end) o dichiarati (`pair.len_min`, `pair.len_max`, appaiate).
- **Perché.** Una variante di lunghezza molto diversa da quella della regione è di norma un'amplificazione aspecifica. In single-end, con letture troncate a lunghezza fissa, il filtro non toglie nulla per costruzione.
- **A che cosa serve.** `lunghezze.tsv` mostra la distribuzione delle lunghezze e fa vedere le varianti di lunghezza inattesa; le letture a questo stadio sono la grandezza su cui si misura la soglia di profondità (3.12).
- **Per cambiarlo.** `asv.len_tol`, `pair.len_min`, `pair.len_max`. Codice: `R/07_asv_length.R`.

### 3.10 Assegnazione tassonomica (S8)

- **Che cosa produce.** In `08_taxonomy/`: `tassonomia.tsv` (un rango per colonna), `bootstrap.tsv` (o `confidenza.tsv` con IdTaxa), `tassonomia.rds`, `copertura_per_gruppo.tsv`, `riepilogo.json` (con la versione di dada2 usata); `difetto_riferimento.tsv` se si dichiarano i taxa con un difetto noto del riferimento; `ranghi_scartati.json` con IdTaxa.
- **Dove nasce.** `steps/s08_taxonomy.py`, `R/08_taxonomy.R`, `R/lib/idtaxa.R`.
- **Metodo.** Due classificatori, a scelta con `tax.classifier`. Il classificatore bayesiano naive di dada2, `assignTaxonomy` (Wang et al., 2007), predefinito, su un riferimento FASTA: i ranghi sono quelli del riferimento (da due a sette, nominati per posizione da Kingdom a Species), e un rango si assegna solo se il bootstrap raggiunge `tax.min_boot` (predefinito del pacchetto, 50). Oppure `DECIPHER::IdTaxa` (Murali et al., 2018) su un classificatore addestrato, con la soglia di confidenza `tax.idtaxa_threshold`, da dichiarare (60 è il predefinito di `DECIPHER::IdTaxa`; più bassa assegna di più con più errori), da Kingdom a Genus. L'assegnazione esatta della specie (`assignSpecies`) non è realizzata.
- **Perché.** Il bayesiano naive è il metodo standard del flusso dada2; IdTaxa è un'alternativa pubblicata che, nel confronto dei suoi autori, assegna meno ranghi errati (meno sovraclassificazioni). Come si comportano i due su uno standard di composizione nota si può verificare con le comunità note del repository (sezione 5). Il **riferimento** (SILVA, GTDB, RDP, Greengenes nel formato del classificatore) decide i nomi possibili e la loro versione: per questo se ne dichiara l'impronta MD5, che S0 verifica, e i risultati riportano nome e versione. Il seme rende ripetibile il bootstrap; la versione di dada2 nell'immagine rende ripetibile anche la scelta fra generi a pari probabilità, e la fase dichiara (`E-S8-03`) quando gira con una versione che non lo è. Il riferimento non si modifica: i taxa con un difetto noto sono elencati, non corretti. La copertura (frazione di varianti con il phylum assegnato) si controlla per classe di campioni.
- **A che cosa serve.** La tassonomia dà il nome alle varianti, per la composizione e per i filtri. Il bootstrap (o la confidenza) dice quanto fidarsi di ciascun rango; la copertura dice se il riferimento è adatto ai campioni (una copertura bassa indica una regione o un ambiente che il riferimento conosce poco); `difetto_riferimento.tsv` elenca le assegnazioni toccate da un difetto noto.
- **Per cambiarlo.** `tax.classifier`, `tax.ref_fasta`, `tax.ref_md5`, `tax.ref_name`, `tax.ref_version`, `tax.min_boot`, `tax.try_rc`, `tax.idtaxa_*`, `tax.ref_bad_taxa`, `qc.min_frac_phylum`. Codice: `R/08_taxonomy.R`, `R/lib/idtaxa.R`.

### 3.11 Oggetto integrato (S10)

- **Che cosa produce.** In `10_phyloseq/`: `ps_integrato.rds` (tutti i campioni dell'inventario, controlli compresi), `varianti.tsv`, `varianti_accessorie.tsv`, `metadati_campioni.tsv`, `colonne_metadati.tsv`, `riepilogo.json`.
- **Dove nasce.** `steps/s10_phyloseq.py`, `R/10_phyloseq.R`.
- **Metodo.** Oggetto `phyloseq` (McMurdie e Holmes, 2013) con conteggi, tassonomia, metadati e sequenze, allineati e verificati sugli identificativi.
- **Perché.** È il formato standard per le analisi di microbioma in R. Gli identificativi `ASV1`, `ASV2`… seguono le letture decrescenti **nei soli campioni biologici** (i positivi metterebbero in testa il ceppo di calibrazione) e non si rinumerano dopo i filtri. L'orientamento della tabella si verifica sui nomi, non sulle dimensioni, e i nomi delle colonne dei metadati si scelgono esplicitamente: `phyloseq` scarterebbe in silenzio campioni non allineati e R rinominerebbe le colonne senza avvisare.
- **A che cosa serve.** È la base di controlli e filtri (S11-S13). `varianti.tsv` permette di seguire una variante lungo tutta la pipeline, per identificativo e per sequenza; `colonne_metadati.tsv` dà i nomi delle colonne da usare nelle analisi ecologiche.
- **Per cambiarlo.** `out.study_columns`, `out.batch_columns`. Codice: `R/10_phyloseq.R`.

### 3.12 Controlli positivi e soglia di profondità (S11)

- **Che cosa produce.** In `11_controls/`: `positivi.tsv` (fedeltà e conformità di ogni controllo positivo), `curve.tsv`, `soglia.json` (la soglia di ogni piastra, con origine e motivo), `profondita_campioni.tsv`, `riepilogo.json`.
- **Dove nasce.** `steps/s11_controls.py`, `R/11_controls.R`, `R/lib/katharoseq.R`.
- **Metodo.** KatharoSeq (Minich et al., 2018). I controlli positivi contengono un solo organismo in diluizione seriale; la **fedeltà** è la quota delle loro letture assegnate a quel taxon. La curva f = x^h / (k' + x^h), con x = log10 della profondità, adattata con `stats::nls`, lega la fedeltà alla profondità; la soglia è la profondità a cui la fedeltà raggiunge `katharoseq.target_sensitivity`. Si adattano la curva aggregata e le curve per piastra sugli stessi punti e si sceglie il modello con l'AIC minore; una curva vale solo se converge, ha abbastanza punti, R² sufficiente e una soglia sostenuta da osservazioni. Un controllo è non conforme se, rispetto ai controlli dello stesso livello nelle altre piastre, ha fedeltà o profondità anomale (z modificato oltre 3,5; Iglewicz e Hoaglin, 1993).
- **Perché.** In bassa biomassa un campione con poche letture può essere dominato dai contaminanti. La soglia è la profondità a cui, nei controlli a un solo organismo della stessa piastra (o di tutte, con la curva aggregata), la quota di letture del taxon atteso raggiunge la sensibilità dichiarata: alla soglia, i controlli hanno al più 1 meno la sensibilità di letture estranee; applicarla ai campioni assume che questi abbiano un carico di contaminanti simile. **La sensibilità obiettivo** dice quanta contaminazione si accetta nei campioni più poveri: più è alta, più la soglia sale e più campioni escono. La soglia è un risultato, non un parametro. Si misura sulle letture senza chimere, perché una soglia su un'altra grandezza non sarebbe confrontabile fra piastre. La conformità tiene conto della concentrazione: a forte diluizione una fedeltà bassa è attesa. Quando una piastra non ha una curva propria valida usa l'aggregata o la mediana delle altre, e lo dichiara (`E-S11-02`).
- **A che cosa serve.** `soglia.json` è la profondità minima che S13 applica ai campioni di ogni piastra; `positivi.tsv` permette di giudicare la corsa (quanti controlli sono conformi); `curve.tsv` la qualità delle curve; `profondita_campioni.tsv` dice quali campioni cadranno sotto soglia.
- **Per cambiarlo.** `katharoseq.*` (taxon atteso, rango, sensibilità, colonna dei livelli, R² minimo), `ctrl.min_positives`, `ctrl.min_positive_pass_frac`, `ctrl.positive_gate`, `qc.min_reads_mode`. Codice: `R/11_controls.R`, `R/lib/katharoseq.R`.

### 3.13 Decontaminazione (S12)

- **Che cosa produce.** In `11_controls/`: `decontam_varianti.tsv`, `decontam_per_piastra.tsv`, `decontam_contaminanti.tsv`, `decontam_riepilogo.json` (le due modalità a confronto), `ps_decontaminato.rds`, `letture_decontaminate.tsv`.
- **Dove nasce.** `steps/s12_decontam.py`, `R/12_decontam.R`, `R/lib/decontaminazione.R`.
- **Metodo.** `decontam::isContaminant` per prevalenza (Davis et al., 2018): una variante è un contaminante se il suo punteggio è sotto `decontam.threshold`. Si calcolano sempre due modalità: aggregata (tutti i biologici contro tutti i negativi) e per piastra (punteggi combinati con la regola di decontam, predefinito il minimo; una piastra con meno di `decontam.min_blanks` negativi confronta i propri biologici con i negativi di tutte le piastre); decide quella dichiarata in `decontam.mode`, l'altra resta come diagnostica.
- **Perché.** In bassa biomassa la contaminazione da reagente può essere una quota rilevante del segnale. Il metodo per prevalenza è l'unico applicabile senza una misura della concentrazione di DNA. **La soglia**: il predefinito di decontam è 0,1; la sua documentazione usa 0,5 come soglia più aggressiva per il metodo per prevalenza, adatta quando i contaminanti pesano molto, come in bassa biomassa. A 0,5 è contaminante, nella modalità aggregata, ogni variante più prevalente nei negativi che nei biologici; nella modalità per piastra con la regola del minimo basta che lo sia in una piastra. I positivi restano fuori dal confronto, perché contengono un organismo aggiunto e, ai livelli diluiti, contaminanti amplificati. La modalità si dichiara prima, non si sceglie dopo aver visto quanto rimuove. Una quota rimossa dai biologici oltre `qc.max_frac_contaminant` ferma l'esecuzione. La decontaminazione precede i filtri di S13 perché usa i controlli negativi, che l'oggetto finale non contiene, e l'insieme completo delle varianti; il filtro di prevalenza da solo non toglierebbe i contaminanti da reagente, che sono prevalenti.
- **A che cosa serve.** `decontam_contaminanti.tsv` dice quali varianti sono state tolte, con punteggio e prevalenze; `decontam_riepilogo.json` confronta le due modalità e dice quale quota di letture è stata tolta da ogni classe. È ciò su cui si giudica se la decontaminazione è proporzionata.
- **Per cambiarlo.** `decontam.threshold`, `decontam.mode`, `decontam.batch_column`, `decontam.batch_combine`, `decontam.min_blanks`, `qc.max_frac_contaminant`; i negativi sono definiti da `ctrl.blank_values` e `ctrl.blank_override_values`. Codice: `R/12_decontam.R`, `R/lib/decontaminazione.R`.

### 3.14 Filtri finali (S13)

- **Che cosa produce.** In `12_final/intermedi/`: `ps_filtrato.rds`, `esclusioni.tsv` (campioni esclusi, con filtro e motivo), `varianti_rimosse.tsv`, `filtri_riepilogo.json`, `letture_finali.tsv`.
- **Dove nasce.** `steps/s13_filtri.py` (filtro per profondità), `R/13_filtri.R` (gli altri).
- **Metodo.** Quattro filtri in ordine dichiarato: (1) profondità, con la soglia di S11 della piastra del campione, confrontata con le letture del campione senza chimere dopo il filtro di lunghezza, la grandezza su cui è stata stimata; (2) tassonomico, le varianti senza phylum e i taxa esclusi (per esempio mitocondri e cloroplasti); (3) prevalenza, una variante resta se ha almeno `prev.min_count` letture in almeno una frazione `prev.min_fraction` dei campioni tenuti; (4) letture minime per campione, `qc.min_reads_final`. Restano solo i campioni biologici.
- **Perché.** L'oggetto finale è destinato alle analisi ecologiche, e un controllo al suo interno sarebbe trattato come un campione ambientale (i controlli si consegnano a parte, 3.16). Mitocondri e cloroplasti portano un 16S di origine batterica, ma sono organelli di cellule eucariotiche, non batteri del campione. Il filtro di prevalenza toglie il rumore delle varianti presenti in pochissimi campioni, a costo di perdere i batteri davvero rari; si lascia spento quando conta la composizione di ogni singolo campione. Il denominatore della prevalenza sono i campioni tenuti dopo profondità e filtro tassonomico, perché la composizione di un campione sotto soglia non è attendibile. Le letture minime rendono confrontabili i campioni nelle analisi a valle.
- **A che cosa serve.** `esclusioni.tsv` dice quali campioni mancano dall'oggetto finale e perché; `varianti_rimosse.tsv` che cosa si è perso e per quale filtro; `filtri_riepilogo.json` quanto ha tolto ciascun filtro.
- **Per cambiarlo.** `filt.remove_na_phylum`, `filt.exclude_taxa`, `prev.apply`, `prev.min_fraction`, `prev.min_count`, `qc.min_reads_final`. Codice: `R/13_filtri.R`, `steps/s13_filtri.py`.

### 3.15 Albero filogenetico (S9, dopo S13)

- **Che cosa produce.** In `09_phylogeny/`: `albero.nwk`, `allineamento.fasta`, `filogenesi.json`; nessun file se la fase si conclude senza albero (quelli di un'esecuzione precedente sono tolti).
- **Dove nasce.** `steps/s09_phylogeny.py`, `R/09_phylogeny.R`, `R/lib/albero.R`.
- **Metodo.** Allineamento multiplo con `DECIPHER::AlignSeqs` (Wright, 2015 e 2020); albero di partenza neighbor-joining, poi massima verosimiglianza con `phangorn` (Schliep, 2011) sul modello GTR+G+I, con ricerca per scambi fra rami vicini (NNI); radicamento al punto medio.
- **Perché.** Si costruisce sulle sole varianti finali: su quelle di S7 conterrebbe contaminanti e organelli che i filtri tolgono. La ricerca è deterministica per scelta: le alternative stocastiche di phangorn dipendono dal seme e costano molte volte tanto; la ricerca NNI può però fermarsi a un massimo locale. Il punto medio è una regola che dipende solo dall'albero, mentre una variante scelta a caso come radice renderebbe le distanze dipendenti da quella scelta; assume tassi simili nei due rami principali, ed è dichiarato in `filogenesi.json`. **`init` propone di costruire sempre l'albero** (`phylo.enabled` vero): su letture corte di un solo tratto del gene le relazioni profonde sono poco sostenute dai dati, e questo è un limite da tenere presente leggendo il risultato, non un motivo per non calcolarlo. Se il dataset non lo permette (varianti finali oltre `phylo.max_seqs`, `E-S9-01`, o meno di quattro, `E-S9-02`) la fase si conclude senza albero, lo dichiara come degradazione nel manifesto e nel report, e l'oggetto finale è consegnato senza albero, con gli stessi conteggi, tassonomia, metadati e sequenze.
- **A che cosa serve.** Alle misure che tengono conto della parentela fra i taxa: le distanze UniFrac (4.4) e altre misure filogenetiche (per esempio la diversità filogenetica di Faith), con i limiti detti sopra.
- **Per cambiarlo.** `phylo.enabled`, `phylo.max_seqs`; allineatore e modello hanno oggi una sola forma realizzata. Codice: `R/09_phylogeny.R`, `steps/s09_phylogeny.py`.

### 3.16 Oggetto finale ed export (S14)

- **Che cosa produce.** In `12_final/`: `ps_final.rds`; `ps_controlli.rds` (i controlli, se l'inventario ne ha); gli export di testo `conteggi.tsv`, `tassonomia.tsv`, `metadati.tsv`, `sequenze.fasta` (e `albero.nwk`, se c'è l'albero); `checksum.sha256`.
- **Dove nasce.** `steps/s14_finale.py`, `R/14_finale.R`.
- **Metodo.** Serializzazione dell'oggetto e validazione strutturale: componenti allineati, solo campioni biologici, nessun campione o variante vuoti, oggetto riletto identico a quello scritto, oggetto ricostruito dai soli export identico a quello serializzato. Con la filogenesi attiva, l'albero di S9 entra nell'oggetto; se S9 ha dichiarato di non poterlo costruire, l'oggetto è consegnato senza. La quota delle letture trattenute dai campioni finali si confronta con `qc.min_frac_reads_retained`.
- **Perché.** Il risultato deve essere usabile anche fuori da R e verificabile da chiunque: gli export non dipendono dalle impostazioni locali e sono identici byte per byte fra due esecuzioni. Una quota trattenuta bassa indica che decontaminazione e filtri hanno tolto la maggior parte del segnale.
- **A che cosa serve.** `ps_final.rds` è l'ingresso delle analisi ecologiche e di qualunque analisi a valle in R; gli export servono a chi lavora con altri programmi; `ps_controlli.rds` al controllo di qualità dei controlli; `checksum.sha256` a verificare l'integrità dei file consegnati (non è una firma: chi altera un file può rifare anche il checksum).
- **Per cambiarlo.** `out.export_flat`, `qc.min_frac_reads_retained`. Codice: `R/14_finale.R`.

### 3.17 Risultati trasversali

- **Tracciamento delle letture.** Ogni fase scrive `letture_<passo>.tsv` nella sua cartella; `runner/tracciamento.py` ricompone la tabella completa, campione per campione, dalle letture grezze all'oggetto finale. Serve a vedere in quale fase un campione perde le letture.
- **Provenienza.** `00_config/resolved.yaml` (la configurazione usata, con ogni valore derivato; una ripresa con valori diversi ne aggiunge un'altra), il manifesto di ogni cartella (impronte dei file, parametri e versione della fase, sorgente del codice) e, con `run.strict_provenance`, `99_logs/ambiente_r.json` (le versioni dei pacchetti R). Servono a verificare che un risultato venga da quella configurazione e da quel codice.
- **Report.** `report/report.html`, con le tabelle anche in `report/tabelle/`, da `report/builder.py`: riassume l'esecuzione leggendo gli artefatti verificati sui manifesti, senza ricalcolare nulla. Le decisioni automatiche (soglie, contaminanti, esclusioni) e le degradazioni sono lette dai file che le registrano. Serve a leggere un'esecuzione in un solo documento.
- **Proposte di configurazione.** `init` (`init/`) propone i valori che i dati sostengono senza ambiguità (troncamento, primer, colonne di controlli e lotto, filogenesi attiva); `soglie` (`soglie/`) propone le soglie di qualità dalle misure di un'esecuzione conclusa, mai più strette della misura. Una soglia proposta così non giudica se la misura è buona: protegge solo dai peggioramenti (campioni aggiunti, una scelta cambiata, lo strumento aggiornato). Nessun valore entra nella configurazione senza un gesto dell'utente.

## 4. Le analisi ecologiche

Lavorano su una copia di `12_final/ps_final.rds`, con un file di configurazione proprio, e scrivono le uscite nella cartella indicata: `campioni.tsv`, le cartelle descritte qui sotto, `configurazione.json` (la configurazione usata), `riepilogo.json` (sintesi e avvisi) e `manifest.json` (impronta di ogni file). Codice: orchestrazione in `esecuzione.py`, parametri in `config.py`, calcoli in `R/analisi.R` e `R/lib/`. Nei nomi dei file, `NN` è la posizione della variabile nella configurazione. Rispondono alle domande di base di uno studio di microbioma: quanto è diverso ogni campione al suo interno (alfa), chi c'è (composizione), quanto i campioni differiscono fra loro (beta, ordinazione) e se i gruppi dello studio differiscono (test).

### 4.1 Campioni e disegno

- **Che cosa produce.** `campioni.tsv`: i campioni di ogni analisi e le esclusioni, con il motivo.
- **Dove nasce.** `R/lib/campioni.R`, `R/valida.R`.
- **Perché.** Indici, abbondanze e distanze si calcolano su tutti i campioni del sottoinsieme; un valore mancante in una variabile toglie il campione dalle sole analisi di quella variabile; un gruppo sotto la dimensione minima esce dai test, non dalle descrizioni. `eco init` segnala le variabili confuse fra loro (due colonne che dividono i campioni negli stessi gruppi, o gruppi interi di una che cadono in un solo gruppo dell'altra): un confronto sull'una misurerebbe anche l'altra.
- **A che cosa serve.** Dice, per ogni analisi, quali campioni sono entrati e perché gli altri no.
- **Per cambiarlo.** `design.variable`, `design.other_variables`, `design.subset`, `design.technical_variables`, `stat.min_group_size`.

### 4.2 Alfa diversità

- **Che cosa produce.** In `alfa/`: `alfa_diversita.tsv` (ricchezza osservata, Shannon, Gini-Simpson per campione), `conteggi_rarefatti.tsv`, e i grafici a scatole per gruppo `gruppi_NN_<variabile>.png` (con le variabili continue, `continua_NN_<variabile>.png`).
- **Dove nasce.** `R/lib/alfa.R` (grafici in `R/lib/test.R` e `R/lib/continua.R`).
- **Metodo.** Una rarefazione senza reinserimento alla profondità `alpha.rarefy_depth` (o, se non dichiarata, a quella del campione più povero), con il seme del file delle analisi ecologiche; poi i tre indici.
- **Perché.** Gli indici dipendono dalla profondità di campionamento, e la rarefazione li rende confrontabili fra campioni. È una scelta discussa: rarefare scarta letture, e c'è chi la sconsiglia (McMurdie e Holmes, 2014) e chi la ritiene oggi il modo più affidabile di controllare la profondità (Schloss, 2024). Lo strumento la applica solo all'alfa diversità; le distanze beta si calcolano sui conteggi non rarefatti (4.4). Senza reinserimento, perché la rarefazione è un sottocampione delle letture osservate: con il reinserimento (il predefinito di `phyloseq::rarefy_even_depth`) una variante può ricevere più letture di quante ne abbia, e si aggiunge una variabilità di campionamento che il dato non ha. **La profondità**: senza un valore dichiarato si rarefà al campione più povero; `eco init` avvisa se questo è sotto un decimo della mediana, nomina i campioni che la abbassano e, se hanno poche letture, suggerisce di alzare `qc.min_reads_final` nella pipeline o di toglierli con `design.subset`; altrimenti, di norma, di lasciarli. Con una profondità dichiarata, i campioni più poveri escono dall'alfa diversità.
- **A che cosa serve.** A confrontare la diversità interna dei campioni: quante varianti hanno (ricchezza) e, insieme a questa, quanto equamente si distribuiscono le letture (Shannon; Gini-Simpson, più sensibile alle varianti abbondanti).
- **Per cambiarlo.** `alpha.rarefy_depth`, `run.seed` (del file delle analisi ecologiche). Codice: `R/lib/alfa.R`.

### 4.3 Composizione tassonomica

- **Che cosa produce.** In `composizione/`: `abbondanze_relative.tsv`, e per ogni variabile `gruppi_NN_<variabile>.tsv` e il grafico `.png`.
- **Dove nasce.** `R/lib/composizione.R`.
- **Metodo.** Agglomerazione al rango `comp.rank` e proporzioni per campione; i `comp.top_n` taxa più abbondanti nei grafici.
- **Perché.** Con letture corte molte varianti non hanno un'assegnazione ai ranghi fini: con il predefinito di `phyloseq::tax_glom` sparirebbero e le proporzioni degli altri taxa risulterebbero gonfiate. Qui restano in una categoria dichiarata, e le proporzioni sommano a 1. `eco init` propone il rango più fine a cui almeno metà delle varianti è assegnata.
- **A che cosa serve.** A vedere quali taxa hanno le quote maggiori in ogni campione e in ogni gruppo. Le proporzioni sono relative e distorte dal numero di copie del 16S e dall'efficienza di amplificazione (sezione 1): non sono quantità assolute.
- **Per cambiarlo.** `comp.rank`, `comp.top_n`, `design.other_variables`. Codice: `R/lib/composizione.R`.

### 4.4 Beta diversità: le distanze

- **Che cosa produce.** In `beta/`: `distanza_<metodo>.tsv` per ogni distanza scelta.
- **Dove nasce.** `R/lib/beta.R`.
- **Metodo.** Sui conteggi non rarefatti: Bray-Curtis sulle proporzioni (Bray e Curtis, 1957); Jaccard in presenza/assenza; Aitchison, cioè la distanza euclidea sulla trasformata log-rapporto centrata, con pseudoconteggio dichiarato (Aitchison, 1982; Gloor et al., 2017); UniFrac non pesato (Lozupone e Knight, 2005) e pesato normalizzato (Lozupone et al., 2007) sull'albero, radicato come dichiarato.
- **Perché.** Le distanze rispondono a domande diverse: Bray-Curtis pesa i taxa abbondanti, Jaccard la sola presenza, Aitchison tratta i dati come composizionali (proporzioni, non quantità), UniFrac tiene conto della parentela fra i taxa e richiede l'albero (3.15). Si calcolano sui conteggi non rarefatti, così non dipendono dal seme; le misure di presenza/assenza (Jaccard, UniFrac non pesato) dipendono però dalla profondità, e differenze di profondità fra gruppi possono apparire come differenze di composizione (Weiss et al., 2017). Jaccard è dichiarata in presenza/assenza perché il predefinito di `vegan::vegdist` ne calcola la forma quantitativa, che è un'altra misura. Lo pseudoconteggio di Aitchison (di norma 1) serve perché lo zero non ha logaritmo; un valore più piccolo rende gli zeri più estremi nel log-rapporto, e dà più peso agli zeri e alle varianti con pochi conteggi.
- **A che cosa serve.** Quanto i campioni differiscono fra loro nella composizione: sono l'ingresso dell'ordinazione e dei test multivariati.
- **Per cambiarlo.** `beta.distances`, `beta.clr_pseudocount`, `beta.unifrac_root`. Codice: `R/lib/beta.R`.

### 4.5 Ordinazione

- **Che cosa produce.** In `ordinazione/`: `<metodo>_<distanza>_coordinate.tsv` (per la NMDS lo stress è nell'intestazione e in `riepilogo.json`), `pcoa_<distanza>_autovalori.tsv`, e i grafici `<metodo>_<distanza>_NN_<variabile>.png`.
- **Dove nasce.** `R/lib/ordinazione.R`.
- **Metodo.** PCoA (Gower, 1966; `ape::pcoa`, Paradis e Schliep, 2019, con correzione dichiarata per gli autovalori negativi) e NMDS (Kruskal, 1964; `vegan::metaMDS`, sulla matrice senza trasformazioni, con `ord.nmds_trymax` tentativi e il seme).
- **Perché.** Le distanze non euclidee danno alla PCoA autovalori negativi, e la correzione è una scelta dichiarata; il segno di ogni asse è fissato, così le coordinate non cambiano fra macchine. Uno stress della NMDS sopra 0,2 è segnalato (`E-ECO-17`): la rappresentazione deforma le distanze.
- **A che cosa serve.** A vedere in due dimensioni come si dispongono i campioni: raggruppamenti, gradienti, campioni anomali. È una rappresentazione, non un test: la differenza fra gruppi si giudica con i test (4.6).
- **Per cambiarlo.** `ord.methods`, `ord.distances`, `ord.pcoa_correction`, `ord.nmds_trymax`. Codice: `R/lib/ordinazione.R`.

### 4.6 Test fra gruppi

- **Che cosa produce.** In `test/`: `alfa_diversita.tsv` (test sugli indici) e `permanova.tsv` (PERMANOVA e test delle dispersioni).
- **Dove nasce.** `R/lib/test.R`.
- **Metodo.** Alfa: Wilcoxon con due gruppi; con più gruppi Kruskal-Wallis (1952) e Wilcoxon per ogni coppia; correzione di Benjamini e Hochberg (1995) per ogni variabile, sull'insieme di tutti gli indici e di tutti i confronti. Beta: PERMANOVA (`vegan::adonis2`; Anderson, 2001) con somme dei quadrati sequenziali, le variabili tecniche come primi termini, `stat.permanova_permutations` permutazioni; accanto, sempre, il test di omogeneità delle dispersioni (`vegan::betadisper`; Anderson, 2006) sugli stessi campioni e sulla stessa matrice.
- **Perché.** Test non parametrici, che non assumono la normalità degli indici. Benjamini-Hochberg controlla la quota attesa di falsi positivi fra i risultati dichiarati significativi (FDR): è meno severa di una correzione che controlla la probabilità di anche un solo falso positivo (FWER, come Bonferroni) ed è adatta a molti confronti esplorativi; la garanzia vale con test indipendenti o positivamente dipendenti. Con gruppi molto piccoli un test per ranghi non può dare una p piccola, comunque i gruppi differiscano; `stat.min_group_size` fissa sotto quale dimensione un gruppo esce dai test. Una PERMANOVA significativa può riflettere una differenza di composizione media o una differenza di variabilità fra gruppi (Anderson e Walsh, 2013): il test delle dispersioni dice se una differenza di variabilità può contribuire al risultato, e se è significativo la PERMANOVA va letta con cautela. Con le variabili tecniche (piastra, corsa) come primi termini, l'effetto della variabile studiata è stimato al netto del lotto; la confusione fra lotto e variabile è dichiarata. **Le p con pochi campioni**: con 999 permutazioni la p minima è 0,001, ma con pochi campioni le disposizioni distinte D possono essere meno delle permutazioni; allora, se le permutazioni sono estratte a caso, la p è riportata a non meno di 1/D (con gruppi di uguale dimensione il minimo effettivo è più alto), con la nota e l'avviso `E-ECO-21`; quando le disposizioni sono così poche da essere enumerate tutte, la p è già esatta. Per il test delle dispersioni D è quello delle sole etichette dei gruppi, per la PERMANOVA quello delle combinazioni del modello.
- **A che cosa serve.** A dire se i gruppi dello studio differiscono nella diversità interna (alfa) e nella composizione (PERMANOVA), e se la differenza di composizione può essere in parte una differenza di variabilità (dispersioni).
- **Per cambiarlo.** `design.variable`, `design.other_variables`, `design.technical_variables`, `stat.min_group_size`, `stat.significance_level`, `stat.permanova_permutations`, `stat.permanova_strata`. Codice: `R/lib/test.R`.

### 4.7 Variabili continue e misure ripetute (facoltative)

- **Che cosa produce.** In `test/`: `alfa_continua.tsv` e `permanova_continua.tsv` per le variabili continue; `alfa_appaiata.tsv`, `permanova_appaiata.tsv` e `soggetti.tsv` per le misure ripetute.
- **Dove nasce.** `R/lib/continua.R`, `R/lib/ripetute.R`.
- **Metodo.** Variabili continue: correlazione di Spearman con gli indici alfa (con Benjamini-Hochberg) e PERMANOVA con la variabile come termine numerico. Misure ripetute: Wilcoxon dei ranghi con segno o Friedman (1937), e PERMANOVA con il soggetto come primo termine e permutazioni entro il soggetto.
- **Perché.** Una misura divisa in gruppi perde l'ordine dei valori; la correlazione per ranghi non presume una relazione lineare. Quando gli stessi soggetti sono campionati in più condizioni, la differenza fra soggetti è di norma maggiore di quella fra condizioni: un test entro soggetto la toglie, confrontando ogni soggetto con sé stesso.
- **A che cosa serve.** A studiare gradienti (una distanza, un tempo, una concentrazione) e disegni con campionamenti ripetuti degli stessi soggetti, che i confronti fra gruppi indipendenti tratterebbero male.
- **Per cambiarlo.** `design.continuous_variables`, `design.subject`, `stat.significance_level` (che fissa anche il minimo di campioni o soggetti), `stat.permanova_permutations`, `stat.permanova_strata`. Codice: `R/lib/continua.R`, `R/lib/ripetute.R`.

## 5. Limiti dichiarati e procedure di verifica

- Le abbondanze del 16S sono **relative** e distorte dal numero di copie del gene per genoma e dall'efficienza di estrazione e amplificazione: lo strumento non stima quantità assolute né corregge queste distorsioni.
- Con IdTaxa la tassonomia si ferma al **genere**; con il bayesiano arriva al rango più fine del riferimento, ma su letture corte il bootstrap è più basso a parità di correttezza e un rango fine (la specie) è poco sostenuto.
- La soglia di profondità richiede **controlli positivi** in diluizione seriale, e assume un carico di contaminanti simile fra controlli e campioni; senza controlli positivi resta solo `qc.min_reads_final`, e lo strumento lo dichiara.
- La decontaminazione richiede **controlli negativi**; con meno di `decontam.min_blanks` negativi l'oggetto non è decontaminato, e lo strumento lo dichiara. Il metodo per prevalenza può marcare come contaminanti taxa veri passati nei bianchi per contaminazione incrociata, e non vede i contaminanti assenti dai bianchi.
- Il filtro di prevalenza toglie anche i **batteri davvero rari**.
- L'alfa diversità si calcola su **una sola rarefazione**; le distanze di presenza/assenza dipendono dalla profondità.
- L'albero su letture corte è **poco risolto** nelle relazioni profonde, la ricerca NNI può fermarsi a un massimo locale e la radice al punto medio è un'approssimazione.
- La riproducibilità byte per byte dimostra che il calcolo è **deterministico**, non che è corretto: la correttezza si giudica sui metodi, descritti qui, e sulle procedure di verifica qui sotto.

**Le procedure di verifica.** Il repository contiene tre procedure, ciascuna con i propri criteri scritti, che chiunque può ripetere:

- **composizione nota** (`dati/comunita_nota/` e `dati/comunita_nota_v3v4/`): standard commerciali di composizione nota; lo script di verifica controlla che ogni genere atteso sia trovato e che i generi non attesi restino sotto una quota dichiarata;
- **confronto con le tabelle di GeneLab** (`dati/confronto_genelab/`): sugli stessi dati, i generi più abbondanti in comune e la correlazione delle abbondanze con le tabelle elaborate dalla pipeline di NASA GeneLab;
- **identità sul riferimento** (`dati/osd734/`): un'esecuzione completa di OSD-734 confrontata, file per file, con le impronte pubblicate, ricavate da un'esecuzione di riferimento.

I criteri sono nei README delle rispettive cartelle.

## 6. Fonti

- Aitchison J (1982). The statistical analysis of compositional data. *Journal of the Royal Statistical Society B* 44: 139-177.
- Anderson MJ (2001). A new method for non-parametric multivariate analysis of variance. *Austral Ecology* 26: 32-46.
- Anderson MJ (2006). Distance-based tests for homogeneity of multivariate dispersions. *Biometrics* 62: 245-253.
- Anderson MJ, Walsh DCI (2013). PERMANOVA, ANOSIM, and the Mantel test in the face of heterogeneous dispersions. *Ecological Monographs* 83: 557-574.
- Benjamini Y, Hochberg Y (1995). Controlling the false discovery rate: a practical and powerful approach to multiple testing. *Journal of the Royal Statistical Society B* 57: 289-300.
- Bray JR, Curtis JT (1957). An ordination of the upland forest communities of southern Wisconsin. *Ecological Monographs* 27: 325-349.
- Callahan BJ et al. (2016). DADA2: high-resolution sample inference from Illumina amplicon data. *Nature Methods* 13: 581-583.
- Callahan BJ, McMurdie PJ, Holmes SP (2017). Exact sequence variants should replace operational taxonomic units in marker-gene data analysis. *The ISME Journal* 11: 2639-2643.
- Davis NM et al. (2018). Simple statistical identification and removal of contaminant sequences in marker-gene and metagenomics data. *Microbiome* 6: 226.
- Friedman M (1937). The use of ranks to avoid the assumption of normality implicit in the analysis of variance. *Journal of the American Statistical Association* 32: 675-701.
- Gloor GB et al. (2017). Microbiome datasets are compositional: and this is not optional. *Frontiers in Microbiology* 8: 2224.
- Gower JC (1966). Some distance properties of latent root and vector methods used in multivariate analysis. *Biometrika* 53: 325-338.
- Iglewicz B, Hoaglin DC (1993). *How to Detect and Handle Outliers*. ASQC Quality Press.
- Kembel SW, Wu M, Eisen JA, Green JL (2012). Incorporating 16S gene copy number information improves estimates of microbial diversity and abundance. *PLoS Computational Biology* 8: e1002743.
- Kruskal JB (1964). Multidimensional scaling by optimizing goodness of fit to a nonmetric hypothesis. *Psychometrika* 29: 1-27.
- Kruskal WH, Wallis WA (1952). Use of ranks in one-criterion variance analysis. *Journal of the American Statistical Association* 47: 583-621.
- Lozupone C, Knight R (2005). UniFrac: a new phylogenetic method for comparing microbial communities. *Applied and Environmental Microbiology* 71: 8228-8235.
- Lozupone CA, Hamady M, Kelley ST, Knight R (2007). Quantitative and qualitative β diversity measures lead to different insights into factors that structure microbial communities. *Applied and Environmental Microbiology* 73: 1576-1585.
- Martin M (2011). Cutadapt removes adapter sequences from high-throughput sequencing reads. *EMBnet.journal* 17: 10-12.
- McLaren MR, Willis AD, Callahan BJ (2019). Consistent and correctable bias in metagenomic sequencing experiments. *eLife* 8: e46923.
- McMurdie PJ, Holmes S (2013). phyloseq: an R package for reproducible interactive analysis and graphics of microbiome census data. *PLoS ONE* 8: e61217.
- McMurdie PJ, Holmes S (2014). Waste not, want not: why rarefying microbiome data is inadmissible. *PLoS Computational Biology* 10: e1003531.
- Minich JJ et al. (2018). KatharoSeq enables high-throughput microbiome analysis from low-biomass samples. *mSystems* 3: e00218-17.
- Murali A, Bhargava A, Wright ES (2018). IDTAXA: a novel approach for accurate taxonomic classification of microbiome sequences. *Microbiome* 6: 140.
- Oksanen J et al. vegan: Community Ecology Package. Pacchetto R, versione 2.7-6 (quella dell'immagine, in `renv.lock`).
- Paradis E, Schliep K (2019). ape 5.0: an environment for modern phylogenetics and evolutionary analyses in R. *Bioinformatics* 35: 526-528.
- Salter SJ et al. (2014). Reagent and laboratory contamination can critically impact sequence-based microbiome analyses. *BMC Biology* 12: 87.
- Schliep KP (2011). phangorn: phylogenetic analysis in R. *Bioinformatics* 27: 592-593.
- Schloss PD (2024). Rarefaction is currently the best approach to control for uneven sequencing effort in amplicon sequence analyses. *mSphere* 9: e00354-23.
- Wang Q, Garrity GM, Tiedje JM, Cole JR (2007). Naive Bayesian classifier for rapid assignment of rRNA sequences into the new bacterial taxonomy. *Applied and Environmental Microbiology* 73: 5261-5267.
- Weiss S et al. (2017). Normalization and microbial differential abundance strategies depend upon data characteristics. *Microbiome* 5: 27.
- Wright ES (2015). DECIPHER: harnessing local sequence context to improve protein multiple sequence alignment. *BMC Bioinformatics* 16: 322.
- Wright ES (2020). RNAconTest: comparing tools for noncoding RNA multiple sequence alignment based on structural consistency. *RNA* 26: 531-540.
- Tutorial di dada2 e issue 938 e 1307 del repository benjjneb/dada2 (parametri standard del filtro; stima dell'errore con qualità raggruppate).
