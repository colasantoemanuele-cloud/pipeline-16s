#!/usr/bin/env python3
"""Confronta i risultati di un'esecuzione su OSD-734 con i checksum attesi.

Chi riproduce i risultati dal repository e dalle fonti pubbliche verifica qui di
aver ottenuto gli stessi byte dell'esecuzione di riferimento. I checksum attesi
stanno accanto a questo script:

* ``checksum_finali.sha256``: i file consegnati di ``12_final/`` (l'oggetto
  finale e i suoi export), nella forma di ``sha256sum``. E' cio' che interessa
  a chi usa i risultati;
* ``checksum_artefatti.tsv``: ogni artefatto elencato nei manifesti di fase, con
  la fase che lo produce. Serve quando un risultato finale non coincide: dice
  qual e' la prima fase, nell'ordine di esecuzione, in cui la divergenza
  compare, ed e' da li' che va cercata la causa.

**Che cosa si confronta.** I soli artefatti dei manifesti di fase, letti dal
disco e ricalcolati: non i manifesti, la configurazione registrata, i log e il
report, che differiscono per costruzione fra due esecuzioni (portano date,
identificativi dell'esecuzione e percorsi).

**A quale ambiente si riferiscono.** All'immagine costruita da
``container/Dockerfile`` del commit indicato nel README della cartella, con la
configurazione ``config_osd734.yaml``. Con un'altra versione di R o dei
pacchetti i risultati possono differire negli ultimi decimali senza che
l'esecuzione sia sbagliata: il confronto prova l'identita', non la correttezza.

    python3 dati/osd734/confronta_risultati.py                     # output/osd734
    python3 dati/osd734/confronta_risultati.py --uscita CARTELLA   # un'altra cartella di output
    python3 dati/osd734/confronta_risultati.py --aggiorna          # riscrive i checksum attesi

L'esito e' 0 se tutti gli artefatti coincidono, 1 altrimenti. Richiede Python
3.11 e la sola libreria standard.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import sys
from pathlib import Path

QUI = Path(__file__).resolve().parent
ATTESI_FINALI = QUI / "checksum_finali.sha256"
ATTESI_ARTEFATTI = QUI / "checksum_artefatti.tsv"
#: La cartella di output indicata da config_osd734.yaml, dalla radice del repository.
USCITA_PREDEFINITA = QUI.parents[1] / "output" / "osd734"

#: Le fasi nell'ordine di esecuzione: S9, la filogenesi, viene dopo S13. Un test
#: tiene questo elenco allineato al grafo della pipeline.
ORDINE = (
    "S0", "S1", "S2", "S3", "S4", "S5", "S6", "S7", "S8",
    "S10", "S11", "S12", "S13", "S9", "S14",
)
#: Le cartelle in cui le fasi scrivono il proprio manifesto: quelle del grafo,
#: e la sottocartella degli intermedi dei filtri finali. I manifesti si cercano
#: solo qui: una ricerca in tutte le sottocartelle leggerebbe anche quelli di
#: un'uscita annidata o di una copia lasciata in una cartella di fase. Un test
#: tiene questo elenco allineato alle cartelle della pipeline.
CARTELLE = (
    "01_input_validation", "02_qc_profiles", "03_filtered", "04_error_models",
    "05_asv_inference", "06_seqtab", "07_chimera", "08_taxonomy", "09_phylogeny",
    "10_phyloseq", "11_controls", "12_final", "12_final/intermedi",
)
#: La cartella e il manifesto dei file consegnati.
CARTELLA_FINALE = "12_final"
MANIFESTO_FINALE = "checksum.sha256"
BLOCCO = 1024 * 1024


def sha256(percorso: Path) -> str:
    """L'impronta SHA-256 di un file, letto a blocchi."""
    impronta = hashlib.sha256()
    with open(percorso, "rb") as file:
        while blocco := file.read(BLOCCO):
            impronta.update(blocco)
    return impronta.hexdigest()


def artefatti_dei_manifesti(uscita: Path) -> list[dict[str, str]]:
    """Gli artefatti elencati nei manifesti di fase, nell'ordine di esecuzione.

    Per ciascuno la fase, la cartella, il nome e i byte dichiarati dal
    manifesto. L'impronta non si prende dal manifesto: si ricalcola dal file.
    I manifesti si leggono nelle sole cartelle di :data:`CARTELLE`, e una fase
    con due manifesti ferma il confronto: non si saprebbe quale dei due
    descrive l'esecuzione.
    """
    per_fase: dict[str, tuple[str, list[dict[str, object]]]] = {}
    origine: dict[str, Path] = {}
    for cartella in CARTELLE:
        for manifesto in sorted((uscita / cartella).glob("manifest_S*.json")):
            documento = json.loads(manifesto.read_text(encoding="utf-8"))
            fase = documento["passo"]
            if fase in per_fase:
                raise SystemExit(
                    f"la fase {fase} ha due manifesti, {origine[fase].relative_to(uscita)} e "
                    f"{manifesto.relative_to(uscita)}: uno dei due e' di un'altra esecuzione, "
                    "toglilo prima di confrontare")
            per_fase[fase] = (documento["cartella"], documento["artefatti"])
            origine[fase] = manifesto
    sconosciute = sorted(set(per_fase) - set(ORDINE))
    if sconosciute:
        raise SystemExit(f"fasi non previste nei manifesti: {', '.join(sconosciute)}")
    righe = []
    for fase in ORDINE:
        cartella, artefatti = per_fase.get(fase, ("", []))
        for voce in sorted(artefatti, key=lambda v: str(v["nome"])):
            righe.append({"fase": fase, "cartella": cartella, "nome": str(voce["nome"]),
                          "byte": str(voce["byte"])})
    return righe


def leggi_attesi() -> list[dict[str, str]]:
    """I checksum attesi di tutti gli artefatti di fase."""
    with open(ATTESI_ARTEFATTI, encoding="utf-8", newline="") as file:
        return list(csv.DictReader(file, delimiter="\t"))


def aggiorna(uscita: Path) -> int:
    """Riscrive i due file dei checksum attesi dall'esecuzione in ``uscita``."""
    righe = artefatti_dei_manifesti(uscita)
    if not righe:
        raise SystemExit(f"nessun manifesto di fase in {uscita}")
    testo = io.StringIO()
    scrittore = csv.writer(testo, delimiter="\t", lineterminator="\n")
    scrittore.writerow(["fase", "cartella", "nome", "sha256", "byte"])
    for riga in righe:
        percorso = uscita / riga["cartella"] / riga["nome"]
        scrittore.writerow([riga["fase"], riga["cartella"], riga["nome"], sha256(percorso),
                            percorso.stat().st_size])
    ATTESI_ARTEFATTI.write_text(testo.getvalue(), encoding="utf-8")
    finali = (uscita / CARTELLA_FINALE / MANIFESTO_FINALE).read_text(encoding="utf-8")
    for riga in finali.splitlines():
        impronta, nome = riga.split(maxsplit=1)
        if sha256(uscita / CARTELLA_FINALE / nome) != impronta:
            raise SystemExit(f"{nome} non corrisponde a {MANIFESTO_FINALE}: esecuzione non integra")
    ATTESI_FINALI.write_text(finali, encoding="utf-8")
    print(f"scritti {ATTESI_ARTEFATTI.name} ({len(righe)} artefatti) e {ATTESI_FINALI.name} "
          f"({len(finali.splitlines())} file)")
    return 0


def confronta(uscita: Path) -> int:
    """Confronta gli artefatti in ``uscita`` con quelli attesi, fase per fase."""
    if not uscita.is_dir():
        print(f"la cartella di output {uscita} non esiste: indicala con --uscita")
        return 1
    attesi = leggi_attesi()
    diversi_per_fase: dict[str, list[str]] = {}
    contati: dict[str, int] = {}
    for riga in attesi:
        fase = riga["fase"]
        contati[fase] = contati.get(fase, 0) + 1
        percorso = uscita / riga["cartella"] / riga["nome"]
        if not percorso.is_file():
            diversi_per_fase.setdefault(fase, []).append(f"{riga['nome']}: manca")
        elif sha256(percorso) != riga["sha256"]:
            diversi_per_fase.setdefault(fase, []).append(f"{riga['nome']}: contenuto diverso")
    # Un artefatto in piu' e' una differenza quanto uno in meno.
    noti = {(r["fase"], r["nome"]) for r in attesi}
    for riga in artefatti_dei_manifesti(uscita):
        if (riga["fase"], riga["nome"]) not in noti:
            diversi_per_fase.setdefault(riga["fase"], []).append(f"{riga['nome']}: non atteso")

    for fase in ORDINE:
        if fase not in contati and fase not in diversi_per_fase:
            continue
        diversi = diversi_per_fase.get(fase, [])
        esito = "identica" if not diversi else f"{len(diversi)} differenze"
        print(f"{fase:<4} {contati.get(fase, 0):>4} artefatti attesi: {esito}")
        for voce in diversi[:10]:
            print(f"       {voce}")
        if len(diversi) > 10:
            print(f"       ... e altre {len(diversi) - 10}")

    finali = [r.split(maxsplit=1)[1] for r in ATTESI_FINALI.read_text(encoding="utf-8").splitlines()]
    print()
    if not diversi_per_fase:
        print(f"RISULTATI IDENTICI: {len(attesi)} artefatti di {len(contati)} fasi, compresi i "
              f"{len(finali)} file consegnati di {CARTELLA_FINALE}/.")
        return 0
    prima = next(f for f in ORDINE if f in diversi_per_fase)
    print(f"RISULTATI DIVERSI. La prima fase in cui compare una differenza e' {prima}: le fasi "
          "successive ne dipendono, quindi la causa va cercata li'. Vedi il README della "
          "cartella, 'Se i risultati non coincidono'.")
    return 1


def main(argomenti: list[str] | None = None) -> int:
    """Punto di ingresso della riga di comando."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--uscita", type=Path, default=USCITA_PREDEFINITA,
                        help="la cartella di output dell'esecuzione (io.out_root)")
    parser.add_argument("--aggiorna", action="store_true",
                        help="riscrive i checksum attesi dall'esecuzione indicata")
    opzioni = parser.parse_args(argomenti)
    return aggiorna(opzioni.uscita) if opzioni.aggiorna else confronta(opzioni.uscita)


if __name__ == "__main__":
    sys.exit(main())
