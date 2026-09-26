# Sottoinsieme di prova di OSD-734

Ventotto campioni del dataset di riferimento, scelti perché le fasi della pipeline
possano essere sviluppate e verificate su dati veri, fino alla validazione finale.

- `selezione.tsv` — la selezione: accession, run ENA, classe, piastra, corsa, letture e
  lunghezza minima del file originale, MD5, e il motivo di ciascuna scelta. È la fonte
  di verità del sottoinsieme.
- `ridotto/` — la versione ridotta: le letture sottocampionate (`fastq/`), le righe dei
  metadati che riguardano i campioni scelti (`metadati/`) e il `manifesto.tsv` con
  letture conservate, lunghezza minima e MD5 di ogni file. È la versione su cui girano i
  test e la catena di integrazione continua.

Entrambe si rigenerano con `scripts/build_test_subset.py`: `completo` ricostruisce il
sottoinsieme con le letture intere, verificando gli MD5; `ridotto` rigenera questa
versione, in modo deterministico, con il seme fissato nello script.

I dati di origine sono pubblici: OSD-734 in NASA GeneLab/OSDR, letture nell'archivio
ENA. La colonna `run_ena` indica il run da cui riscaricare ciascun file; gli MD5 della
selezione coincidono con quelli pubblicati da ENA.
