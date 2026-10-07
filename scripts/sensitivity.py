#!/usr/bin/env python3
"""Analisi di sensibilita' dei parametri a valle dell'oggetto integrato.

Varia un parametro alla volta attorno alla configurazione di un'esecuzione
conclusa, misura l'effetto sull'oggetto finale e applica la regola di decisione
di ``docs/decision_log.md`` (sezione 1), fissata prima dei risultati.

**Come lavora.** I parametri variati incidono solo sulle fasi da S11 in poi:
non serve rieseguire la catena. Per ogni valore lo script copia la cartella di
output dell'esecuzione di partenza con collegamenti fisici, che non occupano
spazio, e la riprende con il parametro cambiato: la valutazione dello stato
rifa' le sole fasi che ne dipendono. L'esecuzione di partenza non viene
toccata, perche' ogni scrittura della pipeline sostituisce il file per
rinomina; i log, che invece crescono in aggiunta, non vengono copiati.

**Che cosa misura**, per ogni valore: se la catena si conclude o quale arresto
scatta; campioni biologici conservati e quanti cambiano stato rispetto alla
configurazione corrente; varianti finali e similarita' di Jaccard del loro
insieme con quello corrente; letture conservate; soglie di profondita' per
piastra e ripieghi; contaminanti rimossi.

**Che cosa scrive**, nella cartella indicata con ``--uscita``:

* ``misure.tsv``: una riga per valore di ogni parametro;
* ``passi.tsv``: i cambiamenti di passo fra valori adiacenti della griglia;
* ``decisioni.json``: l'applicazione della regola a ciascun parametro.

E' rieseguibile: una variante gia' calcolata, riconosciuta dal file
``misure.json`` nella sua copia, non viene rifatta (``--rifai`` la ricalcola).
Va eseguito nell'ambiente della pipeline, il container, con l'esecuzione di
partenza e la cartella delle copie sotto lo stesso punto di montaggio: un
collegamento fisico non attraversa due montaggi.

    python scripts/sensitivity.py --config CONFIG --lavoro CARTELLA --uscita docs/sensibilita
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import shutil
import statistics
import sys
from pathlib import Path
from typing import Any

#: Le griglie della regola di decisione: parametro, valore corrente, valori.
GRIGLIE: dict[str, tuple[float, tuple[float, ...]]] = {
    "katharoseq.target_sensitivity": (0.90, (0.70, 0.75, 0.80, 0.85, 0.90, 0.95)),
    "decontam.threshold": (0.5, (0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7)),
    "prev.min_fraction": (0.01, (0.005, 0.010, 0.015, 0.020, 0.025, 0.030)),
}
#: Valori riportati fuori dalla griglia di decisione.
FUORI_GRIGLIA: dict[str, tuple[float, ...]] = {"prev.min_fraction": (0.05,)}
#: Le modalita' di decontaminazione: nome, parametri del gruppo decontam.
MODALITA: dict[str, dict[str, str]] = {
    "aggregata": {"mode": "aggregate", "batch_combine": "minimum"},
    "per piastra, minimum": {"mode": "batch", "batch_combine": "minimum"},
    "per piastra, fisher": {"mode": "batch", "batch_combine": "fisher"},
}
MODALITA_CORRENTE = "aggregata"

#: Le soglie della regola (docs/decision_log.md, 1.3).
SOGLIE_ASSOLUTE = {"varianti": 0.10, "campioni": 0.05, "letture": 0.05}
FATTORE_SPROPORZIONE = 3.0
MISURE = tuple(SOGLIE_ASSOLUTE)


# --------------------------------------------------------------------------- #
# La regola: funzioni pure, senza la pipeline                                  #
# --------------------------------------------------------------------------- #


def jaccard(a: set[str], b: set[str]) -> float:
    """La similarita' di Jaccard fra due insiemi; 1 se sono entrambi vuoti."""
    return len(a & b) / len(a | b) if a | b else 1.0


def cambiamento(a: dict[str, Any], b: dict[str, Any], biologici: int, letture_correnti: int) -> dict[str, float] | None:
    """I tre cambiamenti di passo fra due configurazioni, o ``None`` se una
    delle due non e' ammissibile.
    """
    if not (a["ammissibile"] and b["ammissibile"]):
        return None
    return {
        "varianti": 1.0 - jaccard(set(a["varianti"]), set(b["varianti"])),
        "campioni": len(set(a["campioni"]) ^ set(b["campioni"])) / biologici,
        "letture": abs(a["letture"] - b["letture"]) / letture_correnti,
    }


def applica_regola(valori: list[float], corrente: float, passi: list[dict[str, float] | None]) -> dict[str, Any]:
    """Applica la regola dei parametri numerici.

    ``valori`` e' la griglia ordinata; ``passi[i]`` il cambiamento fra
    ``valori[i]`` e ``valori[i + 1]``, ``None`` se uno dei due non e'
    ammissibile. Restituisce la decisione con le grandezze che l'hanno
    determinata.
    """
    assert len(passi) == len(valori) - 1
    i = valori.index(corrente)
    assert 0 < i < len(valori) - 1, "il valore corrente deve avere due vicini"
    definiti = [p for p in passi if p is not None]
    mediane = {m: statistics.median(p[m] for p in definiti) if definiti else 0.0 for m in MISURE}
    sproporzione = {m: FATTORE_SPROPORZIONE * mediane[m] for m in MISURE}

    motivi = []
    for lato, passo in (("inferiore", passi[i - 1]), ("superiore", passi[i])):
        if passo is None:
            motivi.append(f"vicino {lato} non ammissibile")
            continue
        for m in MISURE:
            if passo[m] > SOGLIE_ASSOLUTE[m] and passo[m] > sproporzione[m]:
                motivi.append(
                    f"{m} verso il vicino {lato}: {passo[m]:.4f}, oltre la soglia assoluta "
                    f"{SOGLIE_ASSOLUTE[m]} e la soglia di sproporzione {sproporzione[m]:.4f}"
                )
    stabili = [
        valori[j] for j in range(1, len(valori) - 1)
        if j != i and passi[j - 1] is not None and passi[j] is not None
        and all(max(passi[j - 1][m], passi[j][m]) < SOGLIE_ASSOLUTE[m] for m in MISURE)
    ]
    decisione: dict[str, Any] = {
        "corrente": corrente,
        "mediane_dei_passi": mediane,
        "soglie_di_sproporzione": sproporzione,
        "instabile": bool(motivi),
        "motivi_di_instabilita": motivi,
        "alternative_in_zona_stabile": stabili,
    }
    if not motivi:
        decisione.update(scelto=corrente, esito="mantenuto: il valore corrente non e' instabile")
    elif not stabili:
        decisione.update(
            scelto=corrente,
            esito="mantenuto: instabile, ma nessuna alternativa e' in una zona stabile (limite noto)",
        )
    else:
        distanza = min(abs(v - corrente) for v in stabili)
        decisione.update(
            candidati=[v for v in stabili if abs(abs(v - corrente) - distanza) < 1e-12],
            esito="sostituito: instabile, con un'alternativa in zona stabile",
        )
    return decisione


def scegli_modalita(corrente: str, ammissibili: dict[str, bool], somiglianza: dict[str, float]) -> dict[str, Any]:
    """Applica la regola delle modalita' di decontaminazione: si mantiene la
    corrente se ammissibile, altrimenti la piu' simile fra le ammissibili.
    """
    if ammissibili[corrente]:
        return {"corrente": corrente, "scelto": corrente,
                "esito": "mantenuta: la modalita' corrente e' ammissibile"}
    candidate = [m for m, ok in ammissibili.items() if ok]
    if not candidate:
        return {"corrente": corrente, "scelto": corrente,
                "esito": "mantenuta: nessuna modalita' e' ammissibile (limite noto)"}
    return {"corrente": corrente, "scelto": max(candidate, key=lambda m: somiglianza[m]),
            "esito": "sostituita: la modalita' corrente non e' ammissibile"}


# --------------------------------------------------------------------------- #
# Le esecuzioni                                                                #
# --------------------------------------------------------------------------- #


def _copia(origine: Path, destinazione: Path) -> None:
    """Copia l'esecuzione con collegamenti fisici, senza i log e senza il report."""
    shutil.copytree(
        origine, destinazione, copy_function=os.link,
        ignore=shutil.ignore_patterns("99_logs", "report"),
    )


def _tsv(percorso: Path) -> list[dict[str, str]]:
    with open(percorso, encoding="utf-8", newline="") as file:
        return list(csv.DictReader(file, delimiter="\t"))


def misura(uscita: Path, esito: str) -> dict[str, Any]:
    """Le misure di un'esecuzione, lette dai suoi artefatti."""
    misure: dict[str, Any] = {"esito": esito, "ammissibile": esito == "completata"}
    controlli = uscita / "11_controls"
    soglia = json.loads((controlli / "soglia.json").read_text(encoding="utf-8"))
    piastre = soglia["per_piastra"]
    misure["soglie"] = {p: v["valore"] for p, v in piastre.items()}
    misure["ripieghi"] = sorted(p for p, v in piastre.items() if v["origine"].startswith("ripiego"))
    misure["modello_soglia"] = soglia["scelta"]
    # Dopo un arresto il riepilogo della decontaminazione descrive un oggetto a
    # cui non e' stato tolto nulla: le sue cifre non si riportano.
    if not misure["ammissibile"]:
        return misure
    decontam = json.loads((controlli / "decontam_riepilogo.json").read_text(encoding="utf-8"))
    misure["contaminanti"] = decontam["contaminanti_rimossi"]
    misure["letture_rimosse_biologici"] = decontam["letture_rimosse"]["biologico"]
    # I file consegnati stanno in 12_final, gli intermedi dei filtri sotto di essa.
    finale = uscita / "12_final"
    filtri = json.loads(
        (finale / "intermedi" / "filtri_riepilogo.json").read_text(encoding="utf-8"))
    misure["varianti"] = [r["asv_id"] for r in _tsv(finale / "tassonomia.tsv")]
    intestazione = (finale / "conteggi.tsv").read_text(encoding="utf-8").split("\n", 1)[0]
    misure["campioni"] = intestazione.split("\t")[1:]
    misure["letture"] = filtri["letture"]["finali"]
    misure["esclusi"] = filtri["campioni"]["esclusi"]
    misure["varianti_rimosse"] = filtri["varianti"]["rimosse"]
    misure["biologici"] = filtri["campioni"]["biologici"]
    assert len(misure["varianti"]) == filtri["varianti"]["finali"]
    assert len(misure["campioni"]) == filtri["campioni"]["finali"]
    return misure


def esegui_variante(nome: str, base: dict[str, Any], origine: Path, lavoro: Path,
                    cambi: dict[str, dict[str, Any]], rifai: bool) -> dict[str, Any]:
    """Copia l'esecuzione di partenza, la riprende con i parametri cambiati e la
    misura; riusa una variante gia' calcolata.
    """
    from amplicon16s.config.schema import valida
    from amplicon16s.runner.executor import Conclusione, Esecutore
    from amplicon16s.runner.project import ProjectRun

    destinazione = lavoro / nome
    fatto = destinazione / "misure.json"
    if fatto.is_file() and not rifai:
        return json.loads(fatto.read_text(encoding="utf-8"))
    if destinazione.exists():
        shutil.rmtree(destinazione)
    _copia(origine, destinazione)
    dati = json.loads(json.dumps(base))
    dati["io"]["out_root"] = str(destinazione)
    for gruppo, valori in cambi.items():
        dati.setdefault(gruppo, {}).update(valori)
    run = ProjectRun(valida(dati))
    esito = Esecutore(run).esegui()
    if esito.conclusione is Conclusione.COMPLETATA:
        descrizione = "completata"
    else:
        assert esito.punto is not None
        descrizione = f"arresto in {esito.punto.passo} con {esito.punto.codice}"
    misure = misura(destinazione, descrizione)
    misure["fasi_rieseguite"] = [str(r.passo) for r in esito.eseguite]
    fatto.write_text(json.dumps(misure), encoding="utf-8")
    return misure


def _riga(parametro: str, valore: Any, m: dict[str, Any], corrente: dict[str, Any], nota: str = "") -> list[Any]:
    """Una riga di ``misure.tsv``."""
    if not m["ammissibile"]:
        vuote = [""] * 8
    else:
        vuote = [
            len(m["campioni"]), len(set(m["campioni"]) ^ set(corrente["campioni"])),
            len(m["varianti"]), f"{jaccard(set(m['varianti']), set(corrente['varianti'])):.4f}",
            m["letture"], f"{m['letture'] / corrente['letture']:.4f}",
            m["esclusi"]["profondita"], m["varianti_rimosse"]["prevalenza"],
        ]
    soglie = "; ".join(f"{p}: {v}" for p, v in sorted(m["soglie"].items(), key=lambda x: int(x[0])))
    return [
        parametro, valore, m["esito"], *vuote, m["contaminanti"] if m["ammissibile"] else "",
        f"{m['letture_rimosse_biologici']:.4f}" if m["ammissibile"] else "",
        m["modello_soglia"], ", ".join(sorted(m["ripieghi"], key=int)) or "nessuno", soglie, nota,
    ]


def main(argomenti: list[str] | None = None) -> int:
    """Esegue le varianti, scrive le misure e applica la regola."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--config", type=Path, required=True,
                        help="la configurazione dell'esecuzione di partenza, conclusa")
    parser.add_argument("--lavoro", type=Path, required=True,
                        help="la cartella delle copie, sotto lo stesso montaggio dell'esecuzione")
    parser.add_argument("--uscita", type=Path, required=True, help="dove scrivere le tabelle")
    parser.add_argument("--rifai", action="store_true", help="ricalcola anche le varianti gia' fatte")
    opzioni = parser.parse_args(argomenti)

    from amplicon16s.config.schema import carica

    config = carica(opzioni.config)
    base = config.model_dump(mode="json")
    origine = Path(config.io.out_root)
    opzioni.lavoro.mkdir(parents=True, exist_ok=True)
    opzioni.uscita.mkdir(parents=True, exist_ok=True)

    corrente = misura(origine, "completata")
    biologici, letture = corrente["biologici"], corrente["letture"]
    righe: list[list[Any]] = []
    passi_tsv: list[list[Any]] = []
    decisioni: dict[str, Any] = {}

    for parametro, (valore_corrente, griglia) in GRIGLIE.items():
        gruppo, chiave = parametro.split(".")
        assert getattr(getattr(config, gruppo), chiave) == valore_corrente, parametro
        misure = {}
        for valore in (*griglia, *FUORI_GRIGLIA.get(parametro, ())):
            if valore == valore_corrente:
                misure[valore] = corrente
            else:
                nome = f"{gruppo}_{chiave}_{valore}"
                print(f"{parametro} = {valore}", flush=True)
                misure[valore] = esegui_variante(nome, base, origine, opzioni.lavoro,
                                                 {gruppo: {chiave: valore}}, opzioni.rifai)
            nota = ("corrente" if valore == valore_corrente
                    else "fuori dalla griglia di decisione" if valore not in griglia else "")
            righe.append(_riga(parametro, valore, misure[valore], corrente, nota))
        passi = [cambiamento(misure[a], misure[b], biologici, letture)
                 for a, b in zip(griglia, griglia[1:])]
        for (a, b), passo in zip(zip(griglia, griglia[1:]), passi):
            passi_tsv.append([parametro, a, b, *(
                f"{passo[m]:.4f}" if passo else "non definito" for m in MISURE)])
        decisione = applica_regola(list(griglia), valore_corrente, passi)
        if "candidati" in decisione:
            candidati = decisione.pop("candidati")
            decisione["scelto"] = max(candidati, key=lambda v: jaccard(
                set(misure[v]["varianti"]), set(corrente["varianti"])))
        decisioni[parametro] = decisione

    ammissibili, somiglianza = {}, {}
    for nome, parametri in MODALITA.items():
        if nome == MODALITA_CORRENTE:
            assert (config.decontam.mode, config.decontam.batch_combine) == tuple(parametri.values())
            m = corrente
        else:
            print(f"modalita' di decontaminazione: {nome}", flush=True)
            etichetta = "decontam_" + nome.replace(" ", "_").replace(",", "")
            m = esegui_variante(etichetta, base, origine, opzioni.lavoro, {"decontam": parametri},
                                opzioni.rifai)
        ammissibili[nome] = m["ammissibile"]
        somiglianza[nome] = (jaccard(set(m["varianti"]), set(corrente["varianti"]))
                             if m["ammissibile"] else 0.0)
        righe.append(_riga("modalita' di decontaminazione", nome, m, corrente,
                           "corrente" if nome == MODALITA_CORRENTE else ""))
        if not m["ammissibile"]:
            # Diagnostica: la stessa modalita' con il vincolo sospeso, per
            # misurarne l'effetto a valle. Non e' una configurazione candidata.
            sospesa = esegui_variante(
                etichetta + "_vincolo_sospeso", base, origine, opzioni.lavoro,
                {"decontam": parametri, "qc": {"max_frac_contaminant": 1.0}}, opzioni.rifai)
            righe.append(_riga("modalita' di decontaminazione", nome, sospesa, corrente,
                               "diagnostica: qc.max_frac_contaminant sospeso, non candidata"))
    decisioni["modalita' di decontaminazione"] = {
        **scegli_modalita(MODALITA_CORRENTE, ammissibili, somiglianza), "ammissibili": ammissibili,
    }

    intestazione = [
        "parametro", "valore", "esito", "campioni_conservati", "campioni_che_cambiano",
        "varianti_finali", "jaccard_con_la_corrente", "letture_conservate",
        "letture_rispetto_alla_corrente", "esclusi_per_profondita", "varianti_rimosse_prevalenza",
        "contaminanti_rimossi", "letture_biologici_rimosse_come_contaminanti", "modello_soglia",
        "piastre_con_ripiego", "soglie_per_piastra", "nota",
    ]
    for nome, testa, tabella in (
        ("misure.tsv", intestazione, righe),
        ("passi.tsv", ["parametro", "da", "a", *MISURE], passi_tsv),
    ):
        with open(opzioni.uscita / nome, "w", encoding="utf-8", newline="") as file:
            scrittore = csv.writer(file, delimiter="\t", lineterminator="\n")
            scrittore.writerow(testa)
            scrittore.writerows(tabella)
    (opzioni.uscita / "decisioni.json").write_text(
        json.dumps({"soglie_assolute": SOGLIE_ASSOLUTE, "fattore_di_sproporzione": FATTORE_SPROPORZIONE,
                    "campioni_biologici": biologici, "letture_correnti": letture,
                    "decisioni": decisioni}, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    for parametro, d in decisioni.items():
        print(f"{parametro}: {d['esito']} (scelto: {d['scelto']})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
