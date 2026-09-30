"""Costruisce il riferimento tassonomico sintetico per i test di S8.

Il riferimento reale (SILVA 138, 138 MB) non sta nel repository ne'
nell'immagine: i test in CI usano questo, e il reale serve solo ai test sui
dati reali. Le sequenze sono le varianti piu' abbondanti della versione ridotta
del sottoinsieme di prova (gia' nel repository come letture), dopo S7; le linee
tassonomiche sono inventate, con nomi che non esistono (``Phylum_A``,
``Genere_07``), coerenti come albero: ogni genere ha una famiglia, ogni famiglia
un ordine, e cosi' via.

Una famiglia, con i suoi due generi, riproduce il difetto noto della versione 2
di SILVA 138: il rango dell'ordine manca dal percorso, e il nome della famiglia
compare nella colonna dell'ordine. Il file dei taxa difettosi elenca quella
famiglia e quei generi nella forma del file degli autori
(``rank,name,path,ranks_present,ranks_expected``).

Uso, con la tabella di S7 di un'esecuzione di S0-S7 sulla versione ridotta:

    python scripts/costruisci_riferimento_sintetico.py <07_chimera/tabella_asv.rds>

L'archivio si scrive senza nome ne' istante nell'intestazione gzip, cosi' da
essere identico a ogni rigenerazione.
"""

from __future__ import annotations

import gzip
import io
import subprocess
import sys
from pathlib import Path

from amplicon16s.rbridge.runner import trova_rscript

DESTINAZIONE = Path(__file__).resolve().parents[1] / "tests" / "fixtures" / "riferimento_sintetico"
#: Quante varianti entrano nel riferimento.
VARIANTI = 24
#: Il genere i appartiene alla famiglia i // 2, all'ordine i // 4, alla classe
#: i // 8 e al phylum i // 12.
FAMIGLIA_DIFETTOSA = 10


def varianti_per_abbondanza(tabella: Path) -> list[tuple[str, int]]:
    """Le varianti della tabella con le loro letture, dalla piu' abbondante."""
    codice = (
        f"t <- readRDS('{tabella}'); l <- colSums(t); "
        "cat(paste(names(l), format(l, scientific = FALSE, trim = TRUE), sep = '\\t'), sep = '\\n')"
    )
    uscita = subprocess.run(
        [str(trova_rscript()), "--vanilla", "-e", codice],
        capture_output=True, text=True, check=True,
    ).stdout
    coppie = [(s, int(float(n))) for s, n in (r.split("\t") for r in uscita.splitlines() if r)]
    return sorted(coppie, key=lambda c: (-c[1], c[0]))


def linea(genere: int) -> list[str]:
    """La linea tassonomica sintetica del genere indicato, dal regno al genere."""
    famiglia = genere // 2
    ranghi = [
        "Bacteria",
        f"Phylum_{'AB'[genere // 12]}",
        f"Classe_{genere // 8}",
        f"Ordine_{genere // 4}",
        f"Famiglia_{famiglia}",
        f"Genere_{genere:02d}",
    ]
    if famiglia == FAMIGLIA_DIFETTOSA:
        del ranghi[3]  # manca l'ordine: la famiglia scala nella colonna dell'ordine
    return ranghi


def main(argv: list[str]) -> int:
    """Scrive il riferimento e il file dei taxa difettosi in ``tests/fixtures``."""
    if len(argv) != 1:
        print(__doc__, file=sys.stderr)
        return 2
    scelte = varianti_per_abbondanza(Path(argv[0]))[:VARIANTI]
    if len(scelte) < VARIANTI:
        raise SystemExit(f"servono {VARIANTI} varianti, la tabella ne ha {len(scelte)}")

    testo = "".join(f">{';'.join(linea(i))};\n{s}\n" for i, (s, _) in enumerate(scelte))
    DESTINAZIONE.mkdir(parents=True, exist_ok=True)
    buffer = io.BytesIO()
    with gzip.GzipFile(filename="", mode="wb", fileobj=buffer, mtime=0) as archivio:
        archivio.write(testo.encode("ascii"))
    (DESTINAZIONE / "riferimento.fa.gz").write_bytes(buffer.getvalue())

    righe = ["rank,name,path,ranks_present,ranks_expected"]
    difettosi = [i for i in range(VARIANTI) if i // 2 == FAMIGLIA_DIFETTOSA]
    percorso = linea(difettosi[0])[:4]
    righe.append(f"family,Famiglia_{FAMIGLIA_DIFETTOSA},{';'.join(percorso)};,4,5")
    for i in difettosi:
        righe.append(f"genus,Genere_{i:02d},{';'.join(linea(i))};,5,6")
    (DESTINAZIONE / "taxa_difettosi.csv").write_text("\n".join(righe) + "\n", encoding="ascii")
    print(f"{VARIANTI} sequenze in {DESTINAZIONE / 'riferimento.fa.gz'}, "
          f"{len(righe) - 1} taxa difettosi")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
