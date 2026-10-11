# Comunita' nota: uno standard di composizione dichiarata

Cinque corse dello standard ZymoBIOMICS Microbial Community DNA Standard (DNA genomico di
otto batteri in parti uguali, cataloghi D6305 e D6306), sequenziate come controlli
positivi nello studio PRJNA758928 (Illumina MiSeq, letture appaiate 2 x 250, regione V4
con i primer 515F e 806R). Servono a verificare che la pipeline ricostruisca una
composizione che si conosce in anticipo.

## Contenuto della cartella

| File | Che cos'e' |
|---|---|
| `campioni.tsv` | la tabella dei campioni: una riga per corsa (un TSV semplice, non ISA-Tab) |
| `letture_ena.tsv` | per ogni corsa indirizzi, dimensioni e MD5 dei due file dichiarati da ENA, e il lato tenuto |
| `scarica_letture.py` | scarica il file forward di ogni corsa in `fastq/` |
| `config_comunita_nota.yaml` | la configurazione della pipeline, commentata |
| `attesi.tsv` | i generi attesi con la percentuale teorica e la fonte |
| `verifica_composizione.py` | verifica la composizione ottenuta contro quella attesa |
| `FONTI.tsv` | fonte, licenza e controllo di ogni file |

## Le letture

```bash
python3 dati/comunita_nota/scarica_letture.py                  # 5 file in dati/comunita_nota/fastq/
python3 dati/comunita_nota/scarica_letture.py --solo-verifica  # controlla i file presenti
```

Con `--coppie` lo script scarica invece i due file di ogni corsa in `fastq_coppie/`
(`<nome>_R1` il primo file del deposito, `<nome>_R2` il secondo), nell'ordine del
deposito e senza scambiarli: sono l'ingresso della modalita' appaiata della pipeline
(`io.layout: paired-end`), che riconosce da se' una corsa invertita.

Di ogni corsa lo script (`dati/letture_forward.py`) tiene il file le cui letture iniziano
con il primer 515F, dopo aver verificato dimensione e MD5 dichiarati da ENA, e lo salva
come `<corsa>.fastq.gz`. In queste corse e' sempre il primo file.

## Esecuzione

Dalla radice del clone, con le letture e il riferimento scaricati:

```bash
docker run --rm -u "$(id -u):$(id -g)" -v "$(pwd)":/app -w /app <immagine> \
  python3 scripts/esegui.py run --config dati/comunita_nota/config_comunita_nota.yaml
```

Nella tabella dei campioni ogni corsa e' un «positive control» dello studio di origine:
qui lo standard e' l'oggetto dell'analisi, e la configurazione dichiara `ctrl.column:
null` (ogni corsa e' un campione). I parametri condizionali che servono solo con i
controlli sono omessi. La configurazione dichiara anche parametri che su questo dataset
non hanno effetto e che si potrebbero omettere (`katharoseq.target_sensitivity`,
`decontam.threshold`, `decontam.min_blanks`, `decontam.mode` e `qc.max_frac_contaminant`,
perche' non ci sono controlli; `prev.min_fraction` e `prev.min_count`, perche' il filtro
di prevalenza e' spento): sono dichiarati perche' la configurazione dia esattamente il
suo digest, che identifica le esecuzioni. I filtri finali sono quelli di una verifica di composizione:
nessun filtro di prevalenza, nessuna rimozione delle varianti senza phylum, nessuna
soglia sulle letture finali, perche' filtri piu' stretti toglierebbero cio' che la
verifica deve contare.

## Composizione attesa

`attesi.tsv` riporta la composizione teorica per il gene 16S dichiarata dal produttore
(colonna «16S Only» della tabella 1 del manuale dello standard), con il genere e la
famiglia di ogni specie come li scrive il riferimento tassonomico (SILVA 138:
*Escherichia coli* e' `Escherichia-Shigella`, *Lactobacillus fermentum* e'
`Lactobacillus`). I due lieviti dello standard non hanno il 16S.

## Verifica

```bash
python3 dati/comunita_nota/verifica_composizione.py \
  --finale output/comunita_nota/12_final --attesi dati/comunita_nota/attesi.tsv \
  --campioni dati/comunita_nota/campioni.tsv --uscita output/comunita_nota_verifica
```

Ogni corsa si valuta da sola, e lo standard supera un criterio se lo supera ogni corsa;
una corsa di `campioni.tsv` che non e' nell'oggetto finale non li supera. Le quote sono
sulle letture del campione nell'oggetto finale: le letture scartate prima (filtro,
chimere) non hanno un genere e non entrano nel conto.

| Criterio | Soglia |
|---|---|
| 1. Presenza: ogni genere atteso e' trovato | almeno lo 0,1% delle letture del campione |
| 2. Letture su generi non attesi | non oltre il 2% delle letture del campione |

Una variante senza genere conta secondo il rango piu' fine che ha: se e' una famiglia
che contiene un genere atteso resta «senza genere» e la sua quota e' riportata a parte;
altrimenti conta fra le non attese. Lo script riporta anche, senza criterio, il rapporto
fra abbondanza misurata e teorica di ogni genere (sulle sole letture dei generi attesi,
anche in log2): il criterio 1 dice che un genere c'e', il rapporto dice quanto la sua
abbondanza si discosta dall'attesa. Esce con 0 se ogni corsa supera i due criteri, con 1
altrimenti, e scrive `composizione.tsv`, `campioni.tsv` e `riepilogo.json`.
