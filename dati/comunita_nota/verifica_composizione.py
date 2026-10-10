#!/usr/bin/env python3
"""Verifica la composizione di una comunita' nota contro quella attesa.

Ingressi:

* ``--finale``: la cartella ``12_final`` di un'esecuzione (``conteggi.tsv``,
  ``tassonomia.tsv``, ``metadati.tsv``);
* ``--attesi``: la tabella dei generi attesi (``attesi.tsv``: genere e
  famiglia come li scrive il riferimento tassonomico, specie, percentuale
  teorica, fonte);
* ``--campioni``: la tabella dei campioni dati alla pipeline; la prima
  colonna e' il nome del campione. Un campione che vi compare e non e'
  nell'oggetto finale e' riportato come escluso e non supera i criteri: senza
  l'elenco dei campioni attesi una corsa persa passerebbe inosservata.

Uscite, nella cartella ``--uscita``: ``composizione.tsv``, ``campioni.tsv`` e
``riepilogo.json``; sulla console l'esito dei due criteri per ogni campione.

Definizioni. Ogni campione si valuta da solo; la comunita' supera un criterio
se lo supera ogni campione, e un campione escluso dalla pipeline non li supera.

* Criterio 1, presenza: ogni genere atteso ha almeno :data:`QUOTA_MINIMA`
  delle letture del campione.
* Criterio 2: le letture su generi non attesi non superano
  :data:`MASSIMO_NON_ATTESI` del campione. Le letture del campione sono quelle
  dell'oggetto finale: cio' che la pipeline ha scartato prima (letture
  filtrate, chimere) non ha un genere e non entra nel conto. Una variante senza genere (vuoto,
  ``NA`` o etichetta non informativa) conta secondo il rango piu' fine che ha:
  se e' una famiglia che contiene un genere atteso resta «senza genere», e la
  sua quota e' riportata a parte; se si ferma sopra la famiglia, o non e'
  assegnata, conta come non attesa.
* Rapporto fra misurato e teorico, senza criterio: il misurato e'
  l'abbondanza relativa sulle sole letture dei generi attesi, cosi' che
  misurato e teorico sommino entrambi a 100; e' riportato anche in log2.

I nomi si confrontano come nel confronto fra tabelle
(``dati/confronto_genelab/confronta.py``): stringhe esatte, tolti gli spazi ai
bordi e un prefisso di rango.

Esce con 0 se ogni campione supera i due criteri, con 1 altrimenti, con 2 se
gli ingressi non sono leggibili.
"""

from __future__ import annotations

import argparse
import csv
import importlib.util
import json
import math
import sys
from pathlib import Path

QUOTA_MINIMA = 0.001
MASSIMO_NON_ATTESI = 0.02


def _confronto():
    """Il modulo del confronto fra tabelle, da cui viene la normalizzazione dei
    nomi: una sola definizione per i due script."""
    percorso = Path(__file__).resolve().parents[1] / "confronto_genelab" / "confronta.py"
    specifica = importlib.util.spec_from_file_location("confronta_genelab", percorso)
    modulo = importlib.util.module_from_spec(specifica)
    specifica.loader.exec_module(modulo)
    return modulo


_CG = _confronto()
genere_normalizzato = _CG.genere_normalizzato
NON_ASSEGNATO = _CG.NON_ASSEGNATO
IngressoNonValido = _CG.IngressoNonValido

ATTESO, SENZA_GENERE, NON_ATTESO = "atteso", "senza genere, famiglia attesa", "non atteso"


def leggi_attesi(percorso: Path) -> tuple[dict[str, float], set[str]]:
    """Genere atteso -> percentuale teorica, e le famiglie dei generi attesi."""
    try:
        with open(percorso, encoding="utf-8", newline="") as file:
            righe = list(csv.DictReader(file, delimiter="\t"))
    except OSError as errore:
        raise IngressoNonValido(f"{percorso}: {errore}") from errore
    if not righe or not {"genere", "famiglia", "percentuale"} <= set(righe[0]):
        raise IngressoNonValido(f"{percorso}: servono le colonne genere, famiglia e percentuale")
    attesi: dict[str, float] = {}
    for riga in righe:
        genere = genere_normalizzato(riga["genere"])
        if genere == NON_ASSEGNATO or genere in attesi:
            raise IngressoNonValido(f"{percorso}: genere vuoto o ripetuto ({riga['genere']!r})")
        attesi[genere] = float(riga["percentuale"])
    famiglie = {genere_normalizzato(r["famiglia"]) for r in righe} - {NON_ASSEGNATO}
    return attesi, famiglie


def categoria(genere: str | None, famiglia: str | None, attesi: dict[str, float],
              famiglie: set[str]) -> tuple[str, str]:
    """La categoria di una variante e il nome sotto cui si conta."""
    g = genere_normalizzato(genere)
    if g in attesi:
        return ATTESO, g
    if g != NON_ASSEGNATO:
        return NON_ATTESO, g
    f = genere_normalizzato(famiglia)
    if f in famiglie:
        return SENZA_GENERE, f
    return NON_ATTESO, NON_ASSEGNATO


def valuta(conteggi: dict[str, dict[str, int]], attesi: dict[str, float]) -> dict:
    """I due criteri e i rapporti per un campione, dai conteggi per categoria
    (``{categoria: {nome: letture}}``).
    """
    per_genere = conteggi.get(ATTESO, {})
    totale = sum(sum(v.values()) for v in conteggi.values())
    su_attesi = sum(per_genere.values())
    non_attese = sum(conteggi.get(NON_ATTESO, {}).values())
    senza = sum(conteggi.get(SENZA_GENERE, {}).values())
    teorico_totale = sum(attesi.values())
    generi = []
    for genere, percentuale in attesi.items():
        letture = per_genere.get(genere, 0)
        teorico = 100 * percentuale / teorico_totale
        misurato = 100 * letture / su_attesi if su_attesi else None
        rapporto = misurato / teorico if misurato is not None else None
        generi.append({
            "genere": genere, "letture": letture,
            "quota_sul_campione": letture / totale if totale else None,
            "misurato": misurato, "teorico": teorico, "rapporto": rapporto,
            "log2": math.log2(rapporto) if rapporto else None,
            "trovato": bool(totale) and letture / totale >= QUOTA_MINIMA,
        })
    mancanti = [g["genere"] for g in generi if not g["trovato"]]
    quota_non_attese = non_attese / totale if totale else None
    return {
        "letture": totale,
        "quota_attesi": su_attesi / totale if totale else None,
        "quota_non_attesi": quota_non_attese,
        "quota_senza_genere_in_famiglia_attesa": senza / totale if totale else None,
        "criterio_1": {"superato": not mancanti, "mancanti": mancanti},
        "criterio_2": {"superato": quota_non_attese is not None
                       and quota_non_attese <= MASSIMO_NON_ATTESI, "valore": quota_non_attese},
        "generi": generi,
        "non_attesi": dict(sorted(conteggi.get(NON_ATTESO, {}).items(), key=lambda x: (-x[1], x[0]))),
    }


def leggi_finale(finale: Path, attesi: dict[str, float], famiglie: set[str]) -> dict[str, dict]:
    """Per nome di campione, i conteggi per categoria dell'oggetto finale."""
    righe = _CG._righe
    _, metadati = righe(finale / "metadati.tsv")
    intestazione, conteggi = righe(finale / "conteggi.tsv")
    _, tassonomia = righe(finale / "tassonomia.tsv")
    if not metadati or not {"accession", "sample_name"} <= set(metadati[0]):
        raise IngressoNonValido(f"{finale / 'metadati.tsv'}: servono accession e sample_name")
    if tassonomia and not {"Family", "Genus"} <= set(tassonomia[0]):
        raise IngressoNonValido(f"{finale / 'tassonomia.tsv'}: servono le colonne Family e Genus")
    nome_di = {r["accession"]: r["sample_name"] for r in metadati}
    if len(set(nome_di.values())) != len(nome_di):
        raise IngressoNonValido(f"{finale / 'metadati.tsv'}: nomi di campione ripetuti")
    di_variante = {r["asv_id"]: categoria(r["Genus"], r["Family"], attesi, famiglie) for r in tassonomia}
    campioni: dict[str, dict] = {nome_di[c]: {} for c in intestazione[1:]}
    for riga in conteggi:
        specie, nome = di_variante[riga[intestazione[0]]]
        for colonna in intestazione[1:]:
            letture = _CG._intero(riga[colonna], f"conteggi.tsv, colonna {colonna}")
            if letture:
                voce = campioni[nome_di[colonna]].setdefault(specie, {})
                voce[nome] = voce.get(nome, 0) + letture
    return campioni


def _numero(valore: float | None) -> str:
    return "" if valore is None else f"{valore:.6g}"


def main(argomenti: list[str] | None = None) -> int:
    analizzatore = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    analizzatore.add_argument("--finale", type=Path, required=True)
    analizzatore.add_argument("--attesi", type=Path, required=True)
    analizzatore.add_argument("--campioni", type=Path, required=True)
    analizzatore.add_argument("--uscita", type=Path, required=True)
    args = analizzatore.parse_args(argomenti)
    try:
        attesi, famiglie = leggi_attesi(args.attesi)
        conteggi = leggi_finale(args.finale, attesi, famiglie)
        intestazione, righe = _CG._righe(args.campioni)
        dichiarati = [r[intestazione[0]].strip() for r in righe]
        if not dichiarati:
            raise IngressoNonValido(f"{args.campioni}: nessun campione")
        motivi: dict[str, str] = {}
        esclusioni = args.finale / "intermedi" / "esclusioni.tsv"
        if esclusioni.is_file():
            motivi = {r["sample_name"]: f"{r['filtro']}: {r['motivo']}"
                      for r in _CG._righe(esclusioni)[1]}
    except (IngressoNonValido, OSError, KeyError, ValueError, TypeError, OverflowError) as errore:
        print(f"verifica non eseguita: {errore}")
        return 2

    esiti = {nome: valuta(c, attesi) for nome, c in sorted(conteggi.items())}
    esclusi = {
        nome: motivi.get(nome, "non e' nell'oggetto finale")
        for nome in sorted((set(dichiarati) | set(motivi)) - set(esiti))}
    args.uscita.mkdir(parents=True, exist_ok=True)
    with open(args.uscita / "composizione.tsv", "w", encoding="utf-8", newline="") as file:
        scrittore = csv.writer(file, delimiter="\t", lineterminator="\n")
        scrittore.writerow(["campione", "genere", "letture", "quota_sul_campione", "misurato",
                            "teorico", "rapporto", "log2", "trovato"])
        for nome, esito in esiti.items():
            for g in esito["generi"]:
                scrittore.writerow([nome, g["genere"], g["letture"], _numero(g["quota_sul_campione"]),
                                    _numero(g["misurato"]), _numero(g["teorico"]),
                                    _numero(g["rapporto"]), _numero(g["log2"]),
                                    "si" if g["trovato"] else "no"])
    with open(args.uscita / "campioni.tsv", "w", encoding="utf-8", newline="") as file:
        scrittore = csv.writer(file, delimiter="\t", lineterminator="\n")
        scrittore.writerow(["campione", "letture", "quota_attesi", "quota_non_attesi",
                            "quota_senza_genere_in_famiglia_attesa", "criterio_1", "generi_mancanti",
                            "criterio_2", "nota"])
        for nome, e in esiti.items():
            scrittore.writerow([nome, e["letture"], _numero(e["quota_attesi"]),
                                _numero(e["quota_non_attesi"]),
                                _numero(e["quota_senza_genere_in_famiglia_attesa"]),
                                "superato" if e["criterio_1"]["superato"] else "non superato",
                                ", ".join(e["criterio_1"]["mancanti"]),
                                "superato" if e["criterio_2"]["superato"] else "non superato", ""])
        for nome, motivo in esclusi.items():
            scrittore.writerow([nome, "", "", "", "", "non superato", "", "non superato",
                                f"escluso dalla pipeline: {motivo}"])
    tutti_1 = bool(esiti) and not esclusi and all(e["criterio_1"]["superato"] for e in esiti.values())
    tutti_2 = bool(esiti) and not esclusi and all(e["criterio_2"]["superato"] for e in esiti.values())
    riepilogo = {
        "criterio_1": {"descrizione": f"ogni genere atteso con almeno {QUOTA_MINIMA:.1%} delle letture",
                       "superato": tutti_1},
        "criterio_2": {"descrizione": f"non oltre {MASSIMO_NON_ATTESI:.0%} di letture su generi non attesi",
                       "superato": tutti_2},
        "campioni": esiti, "esclusi_dalla_pipeline": esclusi,
        "attesi": attesi, "famiglie_ammesse": sorted(famiglie),
    }
    with open(args.uscita / "riepilogo.json", "w", encoding="utf-8") as file:
        json.dump(riepilogo, file, ensure_ascii=False, indent=2, sort_keys=True)
        file.write("\n")
    for nome, e in esiti.items():
        print(f"{nome}: criterio 1 {'superato' if e['criterio_1']['superato'] else 'NON superato'}"
              + (f" (mancano: {', '.join(e['criterio_1']['mancanti'])})" if e["criterio_1"]["mancanti"] else "")
              + f"; criterio 2 {'superato' if e['criterio_2']['superato'] else 'NON superato'} "
              f"(non attesi {_numero(e['criterio_2']['valore'])}, senza genere in famiglia attesa "
              f"{_numero(e['quota_senza_genere_in_famiglia_attesa'])})")
    for nome, motivo in esclusi.items():
        print(f"{nome}: escluso dalla pipeline ({motivo}): criteri NON superati")
    print(f"comunita': criterio 1 {'superato' if tutti_1 else 'NON superato'}, "
          f"criterio 2 {'superato' if tutti_2 else 'NON superato'}; tabelle in {args.uscita}")
    return 0 if tutti_1 and tutti_2 else 1


if __name__ == "__main__":
    sys.exit(main())
