#!/usr/bin/env python3
"""Confronta, a livello di genere, i risultati della pipeline con la tabella
elaborata da GeneLab per lo stesso dataset.

Ingressi:

* ``--finale``: la cartella ``12_final`` di un'esecuzione (``conteggi.tsv``,
  ``tassonomia.tsv``, ``metadati.tsv``);
* ``--genelab``: la tabella di GeneLab con tassonomia e conteggi
  (``*taxonomy-and-counts.tsv``: una riga per variante, le colonne dei ranghi
  e una colonna per campione);
* ``--abbinamento``: una tabella a due colonne, con intestazione, che per ogni
  nome di campione dei metadati ISA da' l'identificativo con cui il campione
  compare nella tabella di GeneLab;
* ``--escludi`` (facoltativo): un file con un nome di campione per riga, da
  togliere dal confronto (per esempio le corse in cui i due file della coppia
  sono scambiati nell'elaborazione di GeneLab).

Uscite, nella cartella ``--uscita``: ``generi.tsv``, ``campioni.tsv``,
``solo_in_una.tsv`` e ``riepilogo.json``; sulla console l'esito dei due criteri.

Definizioni (le stesse per le due tabelle):

* i conteggi delle varianti si sommano per genere; il nome si confronta come
  stringa esatta dopo aver tolto gli spazi ai bordi e un prefisso di rango
  (``g__``); nessun'altra equivalenza fra nomi;
* un genere vuoto, scritto come assente (:data:`ASSENTI`: ``NA`` e simili) o
  con un'etichetta non informativa (:data:`NON_INFORMATIVE`) e' «non assegnato»: entra nel totale del campione,
  non nella classifica ne' nella correlazione;
* l'abbondanza relativa e' il conteggio del genere sul totale del campione,
  non assegnati compresi;
* si confrontano i campioni presenti in entrambe le tabelle, abbinati per nome;
* criterio 1: almeno :data:`MINIMO_COMUNI` generi in comune fra i
  :data:`QUANTI_PRIMI` con l'abbondanza relativa media piu' alta in ciascuna
  tabella (a pari merito, l'ordine alfabetico);
* criterio 2: mediana, sui campioni, della correlazione di Spearman (ranghi
  medi nei pareggi) fra le abbondanze relative dei generi condivisi, almeno
  :data:`MEDIANA_MINIMA`. Un genere e' condiviso se ha un conteggio maggiore
  di zero in almeno un campione comune in tutte e due le tabelle. In ogni
  campione entrano i generi condivisi diversi da zero in almeno una delle due
  tabelle. Un campione con meno di :data:`MINIMO_GENERI` di questi generi e'
  escluso dalla mediana; se gli esclusi superano :data:`MASSIMO_ESCLUSI` dei
  campioni comuni il criterio non e' superato. Una correlazione non definita
  (valori tutti uguali in una tabella) vale 0.

Esce con 0 se i due criteri sono superati, con 1 se almeno uno non lo e', con
2 se gli ingressi non sono leggibili o non hanno campioni in comune.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import statistics
import sys
from pathlib import Path

#: Parole che, dentro il nome di un genere, ne fanno un'etichetta non
#: informativa (confronto senza distinguere le maiuscole).
NON_INFORMATIVE = (
    "uncultured", "unidentified", "unknown", "incertae_sedis", "incertae sedis",
    "metagenome", "unclassified",
)
#: Scritture di un genere assente, confrontate per intero e senza distinguere
#: le maiuscole: nessuna e' il nome di un genere.
ASSENTI = ("na", "n/a", "nan", "none", "unassigned")
NON_ASSEGNATO = "non assegnato"
QUANTI_PRIMI = 10
MINIMO_COMUNI = 8
MEDIANA_MINIMA = 0.7
MINIMO_GENERI = 5
MASSIMO_ESCLUSI = 0.20

_PREFISSO = re.compile(r"^[a-z]__")


class IngressoNonValido(Exception):
    """Un ingresso manca, non ha le colonne attese o non da' campioni comuni."""


def genere_normalizzato(nome: str | None) -> str:
    """Il nome del genere come si confronta: senza spazi ai bordi e senza il
    prefisso di rango; :data:`NON_ASSEGNATO` se vuoto, ``NA`` o non informativo.
    """
    testo = _PREFISSO.sub("", (nome or "").strip()).strip()
    minuscolo = testo.lower()
    if testo == "" or minuscolo in ASSENTI:
        return NON_ASSEGNATO
    if any(parola in minuscolo for parola in NON_INFORMATIVE):
        return NON_ASSEGNATO
    return testo


def _righe(percorso: Path) -> tuple[list[str], list[dict[str, str]]]:
    """Intestazione e righe di una tabella TSV.

    Le virgolette non hanno un significato: un campo e' cio' che sta fra due
    tabulazioni. Una riga con piu' o meno campi dell'intestazione e' un errore,
    non una riga da completare o da troncare.
    """
    try:
        with open(percorso, encoding="utf-8", newline="") as file:
            lettore = csv.DictReader(file, delimiter="\t", quoting=csv.QUOTE_NONE)
            righe = list(lettore)
            intestazione = list(lettore.fieldnames or [])
    except (OSError, UnicodeDecodeError, csv.Error) as errore:
        raise IngressoNonValido(f"{percorso}: {errore}") from errore
    if len(set(intestazione)) != len(intestazione):
        raise IngressoNonValido(f"{percorso}: colonne con lo stesso nome nell'intestazione")
    for numero, riga in enumerate(righe, start=2):
        if None in riga or None in riga.values():
            raise IngressoNonValido(
                f"{percorso}, riga {numero}: numero di campi diverso dall'intestazione")
    return intestazione, righe


def _intero(valore: str, dove: str) -> int:
    """Un conteggio come intero (le tabelle possono scriverlo con il punto)."""
    try:
        numero = float(valore)
    except ValueError as errore:
        raise IngressoNonValido(f"{dove}: conteggio non numerico {valore!r}") from errore
    if numero != numero or numero in (float("inf"), float("-inf")):
        raise IngressoNonValido(f"{dove}: conteggio non finito {valore!r}")
    if numero < 0 or numero != int(numero):
        raise IngressoNonValido(f"{dove}: conteggio non intero o negativo {valore!r}")
    return int(numero)


def _somma_per_genere(
    righe: list[dict[str, str]], colonna_genere: str, colonne: dict[str, str], dove: str
) -> dict[str, dict[str, int]]:
    """Per campione (nome), i conteggi sommati per genere normalizzato."""
    somme: dict[str, dict[str, int]] = {nome: {} for nome in colonne.values()}
    for riga in righe:
        genere = genere_normalizzato(riga.get(colonna_genere))
        for colonna, nome in colonne.items():
            conteggio = _intero(riga[colonna], f"{dove}, colonna {colonna}")
            if conteggio:
                somme[nome][genere] = somme[nome].get(genere, 0) + conteggio
    return somme


def tabella_nostra(finale: Path) -> tuple[dict[str, dict[str, int]], dict[str, str]]:
    """I conteggi per genere dei campioni dell'oggetto finale, per nome del
    campione, e per ogni nome il motivo di un'eventuale esclusione dai filtri
    (campioni che non sono nell'oggetto finale).
    """
    _, metadati = _righe(finale / "metadati.tsv")
    intestazione, conteggi = _righe(finale / "conteggi.tsv")
    _, tassonomia = _righe(finale / "tassonomia.tsv")
    if not metadati or "sample_name" not in metadati[0] or "accession" not in metadati[0]:
        raise IngressoNonValido(f"{finale / 'metadati.tsv'}: servono le colonne accession e sample_name")
    if tassonomia and "Genus" not in tassonomia[0]:
        raise IngressoNonValido(f"{finale / 'tassonomia.tsv'}: manca la colonna Genus")
    nome_di = {r["accession"]: r["sample_name"] for r in metadati}
    if len(set(nome_di.values())) != len(nome_di):
        raise IngressoNonValido(f"{finale / 'metadati.tsv'}: nomi di campione ripetuti")
    colonne = {c: nome_di[c] for c in intestazione[1:] if c in nome_di}
    if len(colonne) != len(intestazione) - 1:
        raise IngressoNonValido(f"{finale / 'conteggi.tsv'}: campioni assenti dai metadati")
    genere_di = {r["asv_id"]: r["Genus"] for r in tassonomia}
    chiave = intestazione[0]
    for riga in conteggi:
        if riga[chiave] not in genere_di:
            raise IngressoNonValido(f"variante {riga[chiave]} senza tassonomia")
        riga["__genere"] = genere_di[riga[chiave]]
    somme = _somma_per_genere(conteggi, "__genere", colonne, "conteggi.tsv")

    motivi: dict[str, str] = {}
    esclusioni = finale / "intermedi" / "esclusioni.tsv"
    if esclusioni.is_file():
        for riga in _righe(esclusioni)[1]:
            motivi[riga["sample_name"]] = f"{riga['filtro']}: {riga['motivo']}"
    return somme, motivi


def tabella_genelab(percorso: Path, abbinamento: dict[str, str]) -> dict[str, dict[str, int]]:
    """I conteggi per genere dei campioni della tabella di GeneLab, per nome
    del campione: le sole colonne che l'abbinamento riconosce. Un
    identificativo dell'abbinamento che la tabella non ha e' un errore: il
    campione sparirebbe dal confronto senza che nessuno lo veda.
    """
    intestazione, righe = _righe(percorso)
    genere = next((c for c in intestazione if c.strip().lower() == "genus"), None)
    if genere is None:
        raise IngressoNonValido(f"{percorso}: manca la colonna del genere (genus)")
    assenti = sorted(i for i in abbinamento.values() if i not in intestazione)
    if assenti:
        raise IngressoNonValido(
            f"{percorso}: l'abbinamento indica colonne che la tabella non ha: {', '.join(assenti)}")
    colonne = {identificativo: nome for nome, identificativo in abbinamento.items()}
    return _somma_per_genere(righe, genere, colonne, percorso.name)


def leggi_abbinamento(percorso: Path) -> dict[str, str]:
    """Nome del campione -> identificativo nella tabella di GeneLab."""
    intestazione, righe = _righe(percorso)
    if len(intestazione) != 2:
        raise IngressoNonValido(f"{percorso}: servono due colonne (nome del campione, identificativo)")
    nome, identificativo = intestazione
    abbinamento = {r[nome].strip(): r[identificativo].strip() for r in righe}
    if len(abbinamento) != len(righe) or len(set(abbinamento.values())) != len(righe):
        raise IngressoNonValido(f"{percorso}: nomi o identificativi ripetuti")
    return abbinamento


def _ranghi(valori: list[float]) -> list[float]:
    """I ranghi dei valori, con il rango medio nei pareggi."""
    ordine = sorted(range(len(valori)), key=lambda i: valori[i])
    ranghi = [0.0] * len(valori)
    i = 0
    while i < len(ordine):
        j = i
        while j + 1 < len(ordine) and valori[ordine[j + 1]] == valori[ordine[i]]:
            j += 1
        medio = (i + j) / 2 + 1
        for k in range(i, j + 1):
            ranghi[ordine[k]] = medio
        i = j + 1
    return ranghi


def spearman(a: list[float], b: list[float]) -> tuple[float, bool]:
    """La correlazione di Spearman e se e' definita: con i valori tutti uguali
    in una delle due serie non lo e', e vale 0.
    """
    ra, rb = _ranghi(a), _ranghi(b)
    ma, mb = sum(ra) / len(ra), sum(rb) / len(rb)
    va = sum((x - ma) ** 2 for x in ra)
    vb = sum((y - mb) ** 2 for y in rb)
    if va == 0 or vb == 0:
        return 0.0, False
    covarianza = sum((x - ma) * (y - mb) for x, y in zip(ra, rb))
    return covarianza / (va * vb) ** 0.5, True


def _relative(somme: dict[str, int]) -> dict[str, float]:
    """Le abbondanze relative di un campione, sul totale con i non assegnati."""
    totale = sum(somme.values())
    return {g: c / totale for g, c in somme.items()} if totale else {}


def _primi(medie: dict[str, float]) -> list[str]:
    """I generi con l'abbondanza media piu' alta; a pari merito, alfabetico."""
    return [g for g, _ in sorted(medie.items(), key=lambda x: (-x[1], x[0]))[:QUANTI_PRIMI]]


def _quartili(valori: list[float]) -> dict[str, float | None]:
    """Minimo, quartili e massimo (quartili per interpolazione lineare fra i
    valori ordinati, estremi compresi)."""
    if not valori:
        return dict.fromkeys(("minimo", "q1", "mediana", "q3", "massimo"))
    if len(valori) == 1:
        return dict.fromkeys(("minimo", "q1", "mediana", "q3", "massimo"), valori[0])
    q1, mediana, q3 = statistics.quantiles(valori, n=4, method="inclusive")
    return {"minimo": min(valori), "q1": q1, "mediana": mediana, "q3": q3, "massimo": max(valori)}


def confronta(
    nostra: dict[str, dict[str, int]],
    genelab: dict[str, dict[str, int]],
    escludi: set[str] | None = None,
) -> dict:
    """Applica le definizioni alle due tabelle (per campione, conteggi per
    genere gia' normalizzato) e restituisce il riepilogo del confronto.
    """
    escludi = escludi or set()
    comuni = sorted((set(nostra) & set(genelab)) - escludi)
    if not comuni:
        raise IngressoNonValido("nessun campione in comune fra le due tabelle")
    tabelle = {"nostra": nostra, "genelab": genelab}
    relative = {t: {c: _relative(tab[c]) for c in comuni} for t, tab in tabelle.items()}
    generi = {
        t: {g for c in comuni for g, n in tab[c].items() if n > 0 and g != NON_ASSEGNATO}
        for t, tab in tabelle.items()
    }
    medie = {
        t: {g: sum(relative[t][c].get(g, 0.0) for c in comuni) / len(comuni) for g in generi[t]}
        for t in tabelle
    }
    primi = {t: _primi(medie[t]) for t in tabelle}
    in_comune = sorted(set(primi["nostra"]) & set(primi["genelab"]))
    condivisi = generi["nostra"] & generi["genelab"]

    campioni = []
    for c in comuni:
        usati = sorted(g for g in condivisi
                       if tabelle["nostra"][c].get(g, 0) > 0 or tabelle["genelab"][c].get(g, 0) > 0)
        voce = {"campione": c, "generi": len(usati)}
        for t in tabelle:
            totale = sum(tabelle[t][c].values())
            voce[f"letture_{t}"] = totale
            voce[f"quota_non_assegnato_{t}"] = (
                tabelle[t][c].get(NON_ASSEGNATO, 0) / totale if totale else None)
            voce[f"quota_condivisi_{t}"] = (
                sum(tabelle[t][c].get(g, 0) for g in condivisi) / totale if totale else None)
        if len(usati) < MINIMO_GENERI:
            voce.update(correlazione=None, stato="escluso: meno di %d generi" % MINIMO_GENERI)
        else:
            rho, definita = spearman([relative["nostra"][c].get(g, 0.0) for g in usati],
                                     [relative["genelab"][c].get(g, 0.0) for g in usati])
            voce.update(correlazione=rho,
                        stato="nella mediana" if definita else "nella mediana: non definita, vale 0")
        campioni.append(voce)

    correlazioni = [v["correlazione"] for v in campioni if v["correlazione"] is not None]
    esclusi = [v["campione"] for v in campioni if v["correlazione"] is None]
    troppi = len(esclusi) > MASSIMO_ESCLUSI * len(comuni)
    mediana = statistics.median(correlazioni) if correlazioni else None
    quote = {}
    for t, tab in tabelle.items():
        totale = sum(sum(tab[c].values()) for c in comuni)
        quote[t] = {
            "letture": totale,
            "quota_condivisi": sum(tab[c].get(g, 0) for c in comuni for g in condivisi) / totale
            if totale else None,
            "quota_non_assegnato": sum(tab[c].get(NON_ASSEGNATO, 0) for c in comuni) / totale
            if totale else None,
        }
    return {
        "campioni_comuni": comuni,
        "criterio_1": {
            "descrizione": f"almeno {MINIMO_COMUNI} generi in comune fra i {QUANTI_PRIMI} piu' abbondanti",
            "primi_nostra": primi["nostra"], "primi_genelab": primi["genelab"],
            "in_comune": in_comune, "valore": len(in_comune),
            "superato": len(in_comune) >= MINIMO_COMUNI,
        },
        "criterio_2": {
            "descrizione": f"mediana delle correlazioni di Spearman almeno {MEDIANA_MINIMA}",
            "valore": mediana, "campioni_nella_mediana": len(correlazioni),
            "esclusi": esclusi, "esclusi_oltre_il_massimo": troppi,
            "superato": mediana is not None and mediana >= MEDIANA_MINIMA and not troppi,
            "distribuzione": _quartili(correlazioni),
            "quota_letture_nei_generi_condivisi": {t: quote[t]["quota_condivisi"] for t in tabelle},
        },
        "generi": {
            "nostra": len(generi["nostra"]), "genelab": len(generi["genelab"]),
            "condivisi": len(condivisi),
        },
        "quota_non_assegnato": {t: quote[t]["quota_non_assegnato"] for t in tabelle},
        "letture": {t: quote[t]["letture"] for t in tabelle},
        "_medie": medie, "_campioni": campioni, "_condivisi": sorted(condivisi),
    }


def _numero(valore: float | None) -> str:
    return "" if valore is None else f"{valore:.6g}"


def scrivi(esito: dict, uscita: Path, in_piu: dict) -> None:
    """Scrive le tabelle e il riepilogo nella cartella di uscita."""
    uscita.mkdir(parents=True, exist_ok=True)
    medie, condivisi = esito.pop("_medie"), set(esito.pop("_condivisi"))
    campioni = esito.pop("_campioni")
    tutti = sorted(set(medie["nostra"]) | set(medie["genelab"]))
    with open(uscita / "generi.tsv", "w", encoding="utf-8", newline="") as file:
        scrittore = csv.writer(file, delimiter="\t", lineterminator="\n")
        scrittore.writerow(["genere", "media_nostra", "media_genelab", "condiviso"])
        for g in tutti:
            scrittore.writerow([g, _numero(medie["nostra"].get(g)), _numero(medie["genelab"].get(g)),
                                "si" if g in condivisi else "no"])
    with open(uscita / "solo_in_una.tsv", "w", encoding="utf-8", newline="") as file:
        scrittore = csv.writer(file, delimiter="\t", lineterminator="\n")
        scrittore.writerow(["genere", "tabella", "abbondanza_media"])
        for t in ("nostra", "genelab"):
            for g, m in sorted(medie[t].items(), key=lambda x: (-x[1], x[0])):
                if g not in condivisi:
                    scrittore.writerow([g, t, _numero(m)])
    colonne = ["campione", "generi", "correlazione", "stato", "letture_nostra", "letture_genelab",
               "quota_non_assegnato_nostra", "quota_non_assegnato_genelab",
               "quota_condivisi_nostra", "quota_condivisi_genelab"]
    with open(uscita / "campioni.tsv", "w", encoding="utf-8", newline="") as file:
        scrittore = csv.writer(file, delimiter="\t", lineterminator="\n")
        scrittore.writerow(colonne)
        for v in campioni:
            scrittore.writerow([_numero(v[c]) if isinstance(v[c], float) or v[c] is None else v[c]
                                for c in colonne])
    with open(uscita / "riepilogo.json", "w", encoding="utf-8") as file:
        json.dump({**esito, **in_piu}, file, ensure_ascii=False, indent=2, sort_keys=True)
        file.write("\n")


def main(argomenti: list[str] | None = None) -> int:
    analizzatore = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    analizzatore.add_argument("--finale", type=Path, required=True)
    analizzatore.add_argument("--genelab", type=Path, required=True)
    analizzatore.add_argument("--abbinamento", type=Path, required=True)
    analizzatore.add_argument("--escludi", type=Path)
    analizzatore.add_argument("--uscita", type=Path, required=True)
    args = analizzatore.parse_args(argomenti)
    try:
        abbinamento = leggi_abbinamento(args.abbinamento)
        nostra, motivi = tabella_nostra(args.finale)
        genelab = tabella_genelab(args.genelab, abbinamento)
        escludi: set[str] = set()
        if args.escludi:
            escludi = {r.strip() for r in args.escludi.read_text(encoding="utf-8").splitlines()
                       if r.strip() and not r.startswith("#")}
            ignoti = sorted(escludi - set(abbinamento))
            if ignoti:
                raise IngressoNonValido(
                    f"{args.escludi}: campioni che l'abbinamento non conosce: {', '.join(ignoti)}")
        esito = confronta(nostra, genelab, escludi)
    except (IngressoNonValido, OSError, KeyError, TypeError, ValueError, OverflowError) as errore:
        print(f"confronto non eseguito: {errore}")
        return 2
    comuni = set(esito["campioni_comuni"])
    in_piu = {
        "esclusi_dal_confronto": sorted(escludi & set(nostra) & set(genelab)),
        "solo_nella_tabella_di_genelab": {
            nome: motivi.get(nome, "non e' fra i campioni biologici dell'oggetto finale")
            for nome in sorted(set(genelab) - set(nostra))},
        "solo_nella_nostra": sorted(set(nostra) - set(genelab)),
    }
    scrivi(esito, args.uscita, in_piu)
    c1, c2 = esito["criterio_1"], esito["criterio_2"]
    print(f"campioni comuni: {len(comuni)}; esclusi dal confronto: {len(in_piu['esclusi_dal_confronto'])}; "
          f"solo nella tabella di GeneLab: {len(in_piu['solo_nella_tabella_di_genelab'])}; "
          f"solo nella nostra: {len(in_piu['solo_nella_nostra'])}")
    if in_piu["solo_nella_nostra"]:
        print("campioni dell'oggetto finale senza una colonna abbinata, fuori dal confronto: "
              + ", ".join(in_piu["solo_nella_nostra"]))
    print(f"criterio 1 ({c1['descrizione']}): {c1['valore']} in comune: "
          f"{'superato' if c1['superato'] else 'NON superato'}")
    quote = c2["quota_letture_nei_generi_condivisi"]
    print(f"criterio 2 ({c2['descrizione']}): mediana {_numero(c2['valore']) or 'non calcolabile'} "
          f"su {c2['campioni_nella_mediana']} campioni, {len(c2['esclusi'])} esclusi; quota delle "
          f"letture nei generi condivisi: nostra {_numero(quote['nostra'])}, GeneLab "
          f"{_numero(quote['genelab'])}: {'superato' if c2['superato'] else 'NON superato'}")
    print(f"tabelle e riepilogo in {args.uscita}")
    return 0 if c1["superato"] and c2["superato"] else 1


if __name__ == "__main__":
    sys.exit(main())
