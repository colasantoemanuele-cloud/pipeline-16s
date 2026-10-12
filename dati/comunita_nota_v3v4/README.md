# Comunita' nota sulla regione V3-V4

Sei corse dello standard ZymoBIOMICS Microbial Community DNA Standard (DNA genomico di
otto batteri in parti uguali, cataloghi D6305 e D6306), sequenziate come comunita' di
controllo nello studio PRJEB30924 (Illumina MiSeq, letture appaiate 2 x 300, regione
V3-V4 con i primer 341F `CCTACGGGNGGCWGCAG` e 806R `GGACTACNVGGGTWTCTAAT`). La regione,
circa 465 basi con i primer, e' piu' lunga di una lettura: servono a verificare che la
modalita' appaiata della pipeline ricostruisca una composizione che si conosce in
anticipo. La comunita' nota sulla regione V4 e' in `dati/comunita_nota/`.

## Contenuto della cartella

| File | Che cos'e' |
|---|---|
| `campioni.tsv` | la tabella dei campioni: una riga per corsa (un TSV semplice, non ISA-Tab) |
| `letture_ena.tsv` | per ogni corsa indirizzi, dimensioni e MD5 dei due file dichiarati da ENA |
| `scarica_letture.py` | scarica i due file di ogni corsa in `fastq_coppie/` |
| `config_comunita_nota_v3v4_appaiata_naive_bayes.yaml` | modalita' appaiata, classificatore bayesiano ingenuo su SILVA 138 |
| `config_comunita_nota_v3v4_appaiata_idtaxa.yaml` | modalita' appaiata, IdTaxa con il classificatore SILVA SSU r138.2 |
| `FONTI.tsv` | fonte, licenza e controllo di ogni file |

Le letture (`fastq_coppie/`, circa 80 MB) non sono nel repository.

## Le corse

Lo studio confronta due procedure di amplificazione su piu' regioni: le sei corse sono
quelle dello standard sulla regione V3-V4, tre per ciascuna delle due librerie (colonna
`libreria` di `campioni.tsv`). Due delle sei hanno poco piu' di 5.000 coppie di letture.
Lo studio ha depositato le stesse sei corse una seconda volta con altri accession: qui
si usano le prime.

```bash
python3 dati/comunita_nota_v3v4/scarica_letture.py --coppie                  # 12 file in fastq_coppie/
python3 dati/comunita_nota_v3v4/scarica_letture.py --coppie --solo-verifica  # controlla i file presenti
```

I due file di ogni corsa si salvano come `<corsa>_R1.fastq.gz` e `<corsa>_R2.fastq.gz`,
nell'ordine del deposito e senza scambiarli, dopo la verifica di dimensione e MD5
dichiarati da ENA. Le letture portano i primer; nelle reverse il primer e' preceduto da
un numero variabile di basi, e la configurazione lo cerca lungo la lettura
(`filter.primer_removal: search`).

## Le due configurazioni

Nella tabella dei campioni ogni corsa e' una «Mock community» dello studio di origine:
qui lo standard e' l'oggetto dell'analisi, e le configurazioni dichiarano `ctrl.column:
null` (ogni corsa e' un campione). I filtri finali sono quelli di una verifica di
composizione: nessun filtro di prevalenza, nessuna rimozione delle varianti senza
phylum, nessuna soglia sulle letture finali, perche' filtri piu' stretti toglierebbero
cio' che la verifica deve contare. Troncamenti, taglio dei primer, intervallo delle
lunghezze fuse e soglie di qualita' sono quelli scritti da `init --accetta-proposte`;
le soglie sono al valore che non ferma.

Le due configurazioni differiscono solo nel classificatore. Con IdTaxa il classificatore
e' `SILVA_SSU_r138_2_2024.RData` (`python3 dati/riferimento/scarica_addestrato.py 138.2`),
e i taxa da escludere sono scritti come li scrive quel classificatore (`Chloroplast_2`,
`Mitochondria_2`). Entrambe dichiarano la regola rigorosa sulla provenienza
(`run.strict_provenance: true`): si eseguono da un clone senza modifiche in `src/` e
`R/`, nell'immagine indicata in `run.container`.

```bash
docker run --rm -u "$(id -u):$(id -g)" -v "$(pwd)":/app -w /app <immagine> \
  python3 scripts/esegui.py run \
    --config dati/comunita_nota_v3v4/config_comunita_nota_v3v4_appaiata_naive_bayes.yaml
```

## Verifica

La composizione attesa e lo script sono quelli di `dati/comunita_nota/`, il cui README
dice che cosa misura la verifica e con quali soglie. La tabella degli attesi e' quella
del riferimento con cui la tassonomia e' stata assegnata:

```bash
python3 dati/comunita_nota/verifica_composizione.py \
  --finale output/comunita_nota_v3v4_appaiata_naive_bayes/12_final \
  --attesi dati/comunita_nota/attesi.tsv \
  --campioni dati/comunita_nota_v3v4/campioni.tsv --uscita output/comunita_nota_v3v4_verifica
python3 dati/comunita_nota/verifica_composizione.py \
  --finale output/comunita_nota_v3v4_appaiata_idtaxa/12_final \
  --attesi dati/comunita_nota/attesi_idtaxa_silva_138_2.tsv \
  --campioni dati/comunita_nota_v3v4/campioni.tsv --uscita output/comunita_nota_v3v4_idtaxa_verifica
```
