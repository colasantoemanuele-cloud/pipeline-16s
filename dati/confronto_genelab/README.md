# Confronto con le tabelle elaborate da GeneLab

GeneLab pubblica, per molti dataset 16S di OSDR, una tabella con le varianti ottenute
dalla propria elaborazione, la loro tassonomia e i conteggi per campione
(`*taxonomy-and-counts.tsv`). `confronta.py` la confronta, a livello di genere, con
l'oggetto finale della pipeline per lo stesso dataset.

```bash
python3 dati/confronto_genelab/confronta.py --finale <cartella 12_final> \
  --genelab <tabella di GeneLab> --abbinamento <abbinamento.tsv> \
  [--escludi <campioni.txt>] --uscita <cartella>
```

| Ingresso | Che cos'e' |
|---|---|
| `--finale` | la cartella `12_final` di un'esecuzione (`conteggi.tsv`, `tassonomia.tsv`, `metadati.tsv`) |
| `--genelab` | la tabella di GeneLab, scaricata dalla pagina del dataset |
| `--abbinamento` | due colonne con intestazione: il nome del campione nei metadati ISA e l'identificativo con cui compare nella tabella di GeneLab |
| `--escludi` | facoltativo: un nome di campione per riga, da togliere dal confronto |

Gli abbinamenti dei dataset del repository sono accanto ai loro dati
(`dati/osd170/abbinamento_genelab.tsv`, `dati/osd276/abbinamento_genelab.tsv`,
`dati/osd191/abbinamento_genelab.tsv`). Per
OSD-276 vanno escluse le corse in cui i due file pubblicati da OSDR sono scambiati
(`dati/osd276/corse_invertite_osdr.txt`): l'elaborazione di GeneLab di quelle corse non
e' confrontabile.

## Che cosa misura

- I conteggi delle varianti si sommano per genere. Il nome si confronta come stringa
  esatta, tolti gli spazi ai bordi e un prefisso di rango (`g__`): nessun'altra
  equivalenza fra nomi (due scritture dello stesso genere restano due generi).
- Un genere vuoto, scritto come assente (`NA`, `N/A`, `NaN`, `None`, `Unassigned`,
  senza distinguere le maiuscole) o con un'etichetta non informativa (`uncultured`,
  `unidentified`, `unknown`, `incertae sedis`, `metagenome`, `unclassified`) e' «non
  assegnato»: entra nel totale del campione, non nella classifica ne' nella
  correlazione.
- L'abbondanza relativa e' il conteggio del genere sul totale del campione.
- Si confrontano i campioni presenti in entrambe le tabelle, abbinati per nome. I
  controlli e i campioni esclusi dai filtri della pipeline non sono nell'oggetto finale:
  lo script li elenca, con il motivo dell'esclusione.

| Criterio | Soglia |
|---|---|
| 1. Generi in comune fra i 10 piu' abbondanti di ciascuna tabella (abbondanza relativa media; a pari merito, ordine alfabetico) | almeno 8 |
| 2. Mediana, sui campioni, della correlazione di Spearman fra le abbondanze relative dei generi condivisi | almeno 0,7 |

Un genere e' condiviso se ha letture in almeno un campione comune in tutte e due le
tabelle. In ogni campione entrano nella correlazione i generi condivisi diversi da zero
in almeno una delle due tabelle; un campione con meno di 5 di questi generi e' escluso
dalla mediana, e se gli esclusi superano il 20% dei campioni comuni il criterio 2 non e'
superato. Una correlazione non definita vale 0. Accanto al criterio 2 lo script riporta
la quota delle letture che cade nei generi condivisi, per le due tabelle: dice quanta
parte dei dati la correlazione non guarda.

## Uscite

`generi.tsv` (abbondanza media di ogni genere nelle due tabelle), `solo_in_una.tsv`
(generi presenti in una sola), `campioni.tsv` (per campione: generi usati, correlazione,
quote dei non assegnati e dei condivisi) e `riepilogo.json`. Sulla console l'esito dei
due criteri. Esce con 0 se sono superati entrambi, con 1 se almeno uno non lo e', con 2
se il confronto non si puo' eseguire: ingressi non leggibili, una riga con un numero di
campi diverso dall'intestazione, un identificativo dell'abbinamento che la tabella non
ha, un campione da escludere che l'abbinamento non conosce, nessun campione in comune.
Le virgolette nelle tabelle sono testo, non delimitatori.

Le due elaborazioni differiscono per costruzione: GeneLab fonde le due letture di ogni
coppia, assegna la tassonomia con IdTaxa e non applica filtri finali; la pipeline toglie
contaminanti, organelli e varianti rare e, secondo la configurazione, usa le sole forward
o fonde le coppie, e assegna la tassonomia con il classificatore bayesiano ingenuo o con
IdTaxa. Il confronto misura quanto, a dispetto delle differenze che restano, le due
descrivono la stessa comunita'. I nomi si confrontano come stringhe esatte: due
classificatori che scrivono lo stesso genere in due modi danno due generi.

Per confrontare con la stessa tabella due configurazioni della pipeline (le sole forward
e le coppie fuse, o i due classificatori), ciascuna si confronta sugli stessi campioni:
quelli che una delle due ha escluso con i filtri si tolgono dal confronto dell'altra con
`--escludi`.
