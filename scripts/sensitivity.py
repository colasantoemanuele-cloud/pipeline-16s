#!/usr/bin/env python3
"""Analisi di sensibilita' dei parametri a valle dell'oggetto integrato.

Varia un parametro alla volta attorno alla configurazione di un'esecuzione
conclusa, misura l'effetto sull'oggetto finale e applica la regola di decisione
di ``docs/decision_log.md`` (sezione 1), fissata prima dei risultati.

**Nulla di un dataset sta nel codice.** Le griglie vengono da un file indicato
con ``--griglie`` (per il dataset di riferimento ``docs/sensibilita/griglie.yaml``);
il valore corrente di ogni parametro, e la modalita' corrente, si leggono
dalla configurazione. Il file ha due voci:

* ``numerici``: per parametro (``gruppo.chiave``), la ``griglia`` ordinata e a
  passo uniforme, in cui il valore della configurazione deve avere due vicini,
  e facoltativamente ``fuori_griglia``, valori riportati ma non decisi;
* ``modalita``: un ``nome`` e, in ``valori``, per ogni modalita' i parametri
  che la definiscono; e' corrente quella i cui parametri coincidono con la
  configurazione.

**Come lavora.** I parametri variati incidono solo sulle fasi a valle
dell'oggetto integrato: non serve rieseguire la catena. Per ogni valore lo
script copia la cartella di output dell'esecuzione di partenza con
collegamenti fisici, che non occupano spazio, e la riprende con il parametro
cambiato: la valutazione dello stato rifa' le sole fasi che ne dipendono.
L'esecuzione di partenza non viene toccata, perche' ogni scrittura della
pipeline sostituisce il file per rinomina; i log, che invece crescono in
aggiunta, non vengono copiati. I percorsi relativi della configurazione sono
resi assoluti come fa la riga di comando.

**Ammissibilita'.** Un valore non e' ammissibile solo se produce una delle
condizioni che la regola prevede (``ARRESTI_PREVISTI``, ``SVUOTATI_DALLA_PREVALENZA``).
Ogni altro arresto e ogni errore fermano lo script: un guasto non e' una
misura, e contarlo come valore non ammissibile lo nasconderebbe nella tabella.

**Che cosa misura**, per ogni valore: se la catena si conclude o quale arresto
scatta; campioni biologici conservati e quanti cambiano stato rispetto alla
configurazione corrente; varianti finali e similarita' di Jaccard del loro
insieme con quello corrente; letture conservate; soglie di profondita' per
piastra con la loro origine; contaminanti rimossi.

**Che cosa scrive**, nella cartella indicata con ``--uscita``:

* ``misure.tsv``: una riga per valore di ogni parametro;
* ``passi.tsv``: i cambiamenti di passo fra valori adiacenti della griglia;
* ``decisioni.json``: l'applicazione della regola a ciascun parametro.

E' rieseguibile: una variante gia' calcolata, riconosciuta dal file
``misure.json`` nella sua copia, non viene rifatta (``--rifai`` la ricalcola).
Va eseguito nell'ambiente della pipeline, il container, con l'esecuzione di
partenza e la cartella delle copie sotto lo stesso punto di montaggio: un
collegamento fisico non attraversa due montaggi.

    python scripts/sensitivity.py --config CONFIG --griglie GRIGLIE \
        --lavoro CARTELLA --uscita CARTELLA_DELLE_TABELLE
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

#: Le soglie della regola (docs/decision_log.md, 1.3).
SOGLIE_ASSOLUTE = {"varianti": 0.10, "campioni": 0.05, "letture": 0.05}
FATTORE_SPROPORZIONE = 3.0
MISURE = tuple(SOGLIE_ASSOLUTE)

#: Gli arresti che la regola prevede (1.3, punto 1): le letture rimosse come
#: contaminanti oltre qc.max_frac_contaminant e la frazione di letture
#: trattenute sotto qc.min_frac_reads_retained. Un valore che li produce non e'
#: ammissibile; qualunque altro arresto ferma lo script.
ARRESTI_PREVISTI = {
    "E-S12-02": "letture rimosse come contaminanti oltre qc.max_frac_contaminant",
    "E-S14-01": "frazione di letture trattenute sotto qc.min_frac_reads_retained",
}
#: Un codice che la pipeline usa anche per altre cause e' un arresto previsto
#: solo se il suo dettaglio nomina il parametro della condizione: E-S14-01 e'
#: anche l'oggetto finale non valido o un ingresso alterato, che sono guasti.
CAUSA_RICHIESTA = {"E-S14-01": "qc.min_frac_reads_retained"}
#: La terza condizione della regola, i campioni svuotati dal filtro di
#: prevalenza: la pipeline la dichiara senza fermarsi, e la regola la tratta
#: comunque come non ammissibile.
SVUOTATI_DALLA_PREVALENZA = "E-S13-02"
#: Il vincolo che un arresto per i contaminanti lascia sospendere, per la sola
#: diagnostica: il parametro e il valore che non pone alcun limite.
VINCOLO_SOSPESO = {"E-S12-02": {"qc": {"max_frac_contaminant": 1.0}}}


class ErroreSensibilita(RuntimeError):
    """Un guasto dell'analisi: non una misura, e lo script si ferma."""


def leggi_griglie(percorso: Path) -> tuple[dict[str, dict[str, tuple[float, ...]]], str, dict[str, dict[str, Any]]]:
    """Le griglie numeriche, il nome delle modalita' e le modalita' dal file.

    Restituisce, per parametro, ``griglia`` e ``fuori_griglia``; il nome con
    cui le modalita' compaiono nelle tabelle; e per modalita' i parametri
    ``gruppo.chiave`` che la definiscono.
    """
    import yaml

    dati = yaml.safe_load(percorso.read_text(encoding="utf-8")) or {}
    numerici = {
        parametro: {"griglia": tuple(voce["griglia"]),
                    "fuori_griglia": tuple(voce.get("fuori_griglia", ()))}
        for parametro, voce in (dati.get("numerici") or {}).items()
    }
    modalita = dati.get("modalita") or {}
    return numerici, modalita.get("nome", "modalita"), dict(modalita.get("valori") or {})


def valore_corrente(config: Any, parametro: str) -> Any:
    """Il valore di ``gruppo.chiave`` nella configurazione."""
    gruppo, chiave = parametro.split(".")
    return getattr(getattr(config, gruppo), chiave)


def modalita_corrente(config: Any, modalita: dict[str, dict[str, Any]]) -> str:
    """La modalita' i cui parametri coincidono con la configurazione."""
    correnti = [nome for nome, parametri in modalita.items()
                if all(valore_corrente(config, p) == v for p, v in parametri.items())]
    if len(correnti) != 1:
        raise ErroreSensibilita(
            f"la configurazione corrisponde a {len(correnti)} delle modalita' del file "
            f"delle griglie ({', '.join(modalita)}): ne serve una sola")
    return correnti[0]


def _per_gruppo(parametri: dict[str, Any]) -> dict[str, dict[str, Any]]:
    """Da ``gruppo.chiave: valore`` a ``gruppo: {chiave: valore}``."""
    cambi: dict[str, dict[str, Any]] = {}
    for parametro, valore in parametri.items():
        gruppo, chiave = parametro.split(".")
        cambi.setdefault(gruppo, {})[chiave] = valore
    return cambi


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


def _degradazioni(uscita: Path) -> set[str]:
    """I codici delle degradazioni dichiarate nei manifesti di fase dell'esecuzione."""
    codici: set[str] = set()
    for manifesto in uscita.rglob("manifest_S*.json"):
        contenuto = json.loads(manifesto.read_text(encoding="utf-8"))
        codici |= {d["codice"] for d in contenuto.get("degradazioni", [])}
    return codici


def misura(uscita: Path, arresto: str | None = None) -> dict[str, Any]:
    """Le misure di un'esecuzione, lette dai suoi artefatti.

    ``arresto`` e' il codice dell'arresto previsto, se la catena si e' fermata.
    Un'esecuzione conclusa che dichiara campioni svuotati dal filtro di
    prevalenza non e' ammissibile: la regola la esclude.
    """
    if arresto is not None:
        esito = f"arresto: {ARRESTI_PREVISTI[arresto]} ({arresto})"
    elif SVUOTATI_DALLA_PREVALENZA in _degradazioni(uscita):
        esito = f"non ammissibile: campioni svuotati dal filtro di prevalenza ({SVUOTATI_DALLA_PREVALENZA})"
    else:
        esito = "completata"
    misure: dict[str, Any] = {"esito": esito, "ammissibile": esito == "completata",
                              "codice_arresto": arresto}
    controlli = uscita / "11_controls"
    soglia = json.loads((controlli / "soglia.json").read_text(encoding="utf-8"))
    piastre = soglia["per_piastra"] or {}
    misure["soglie"] = {p: v["valore"] for p, v in piastre.items()}
    misure["origini"] = {p: v["origine"] for p, v in piastre.items()}
    misure["modello_soglia"] = soglia["scelta"]
    # Dopo un arresto il riepilogo della decontaminazione descrive un oggetto a
    # cui non e' stato tolto nulla: le sue cifre non si riportano.
    if arresto is not None:
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
    arresto = None
    if esito.conclusione is not Conclusione.COMPLETATA:
        previsto = (
            esito.punto is not None and esito.punto.codice in ARRESTI_PREVISTI
            and CAUSA_RICHIESTA.get(esito.punto.codice, "") in esito.punto.dettaglio
        )
        if not previsto:
            dove = "senza punto di ripresa" if esito.punto is None else (
                f"in {esito.punto.passo} con {esito.punto.codice} ({esito.punto.dettaglio})")
            raise ErroreSensibilita(
                f"la variante {nome} si e' fermata {dove}: non e' un arresto previsto dalla "
                f"regola ({', '.join(ARRESTI_PREVISTI)}). Log in {destinazione / '99_logs'}")
        arresto = esito.punto.codice
    misure = misura(destinazione, arresto)
    misure["fasi_rieseguite"] = [str(r.passo) for r in esito.eseguite]
    fatto.write_text(json.dumps(misure), encoding="utf-8")
    return misure


def _ordine(piastra: str) -> tuple[bool, int, str]:
    """Le piastre numerate in ordine numerico, le altre dopo."""
    return (not piastra.isdigit(), int(piastra) if piastra.isdigit() else 0, piastra)


def _riga(parametro: str, valore: Any, m: dict[str, Any], corrente: dict[str, Any], nota: str = "") -> list[Any]:
    """Una riga di ``misure.tsv``."""
    if "campioni" not in m:
        vuote = [""] * 10
    else:
        vuote = [
            len(m["campioni"]), len(set(m["campioni"]) ^ set(corrente["campioni"])),
            len(m["varianti"]), f"{jaccard(set(m['varianti']), set(corrente['varianti'])):.4f}",
            m["letture"], f"{m['letture'] / corrente['letture']:.4f}",
            m["esclusi"]["profondita"], m["varianti_rimosse"]["prevalenza"],
            m["contaminanti"], f"{m['letture_rimosse_biologici']:.4f}",
        ]
    piastre = sorted(m["soglie"], key=_ordine)
    soglie = "; ".join(f"{p}: {m['soglie'][p]}" for p in piastre)
    non_proprie = "; ".join(f"{p}: {m['origini'][p]}" for p in piastre
                            if m["origini"][p] != "propria")
    return [parametro, valore, m["esito"], *vuote, m["modello_soglia"],
            non_proprie or "nessuna", soglie, nota]


def main(argomenti: list[str] | None = None) -> int:
    """Esegue le varianti, scrive le misure e applica la regola."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--config", type=Path, required=True,
                        help="la configurazione dell'esecuzione di partenza, conclusa")
    parser.add_argument("--griglie", type=Path, required=True,
                        help="il file YAML con le griglie numeriche e le modalita'")
    parser.add_argument("--lavoro", type=Path, required=True,
                        help="la cartella delle copie, sotto lo stesso montaggio dell'esecuzione")
    parser.add_argument("--uscita", type=Path, required=True, help="dove scrivere le tabelle")
    parser.add_argument("--rifai", action="store_true", help="ricalcola anche le varianti gia' fatte")
    opzioni = parser.parse_args(argomenti)

    from amplicon16s.cli import _percorsi_assoluti
    from amplicon16s.config.schema import carica

    config = _percorsi_assoluti(carica(opzioni.config))
    numerici, nome_modalita, modalita = leggi_griglie(opzioni.griglie)
    base = config.model_dump(mode="json")
    origine = Path(config.io.out_root)
    lavoro = opzioni.lavoro.resolve()
    lavoro.mkdir(parents=True, exist_ok=True)
    opzioni.uscita.mkdir(parents=True, exist_ok=True)

    corrente = misura(origine)
    if not corrente["ammissibile"]:
        raise ErroreSensibilita(f"l'esecuzione di partenza non e' ammissibile: {corrente['esito']}")
    biologici, letture = corrente["biologici"], corrente["letture"]
    righe: list[list[Any]] = []
    passi_tsv: list[list[Any]] = []
    decisioni: dict[str, Any] = {}

    for parametro, voce in numerici.items():
        griglia = voce["griglia"]
        valore_della_config = valore_corrente(config, parametro)
        if valore_della_config not in griglia:
            raise ErroreSensibilita(
                f"{parametro} vale {valore_della_config} nella configurazione, che non e' "
                f"nella griglia {griglia}")
        if griglia.index(valore_della_config) in (0, len(griglia) - 1):
            raise ErroreSensibilita(
                f"{parametro} vale {valore_della_config}, a un estremo della griglia {griglia}: "
                "la regola confronta il valore corrente con due vicini")
        misure = {}
        for valore in (*griglia, *voce["fuori_griglia"]):
            if valore == valore_della_config:
                misure[valore] = corrente
            else:
                nome = f"{parametro.replace('.', '_')}_{valore}"
                print(f"{parametro} = {valore}", flush=True)
                misure[valore] = esegui_variante(nome, base, origine, lavoro,
                                                 _per_gruppo({parametro: valore}), opzioni.rifai)
            nota = ("corrente" if valore == valore_della_config
                    else "fuori dalla griglia di decisione" if valore not in griglia else "")
            righe.append(_riga(parametro, valore, misure[valore], corrente, nota))
        passi = [cambiamento(misure[a], misure[b], biologici, letture)
                 for a, b in zip(griglia, griglia[1:])]
        for (a, b), passo in zip(zip(griglia, griglia[1:]), passi):
            passi_tsv.append([parametro, a, b, *(
                f"{passo[m]:.4f}" if passo else "non definito" for m in MISURE)])
        decisione = applica_regola(list(griglia), valore_della_config, passi)
        if "candidati" in decisione:
            candidati = decisione.pop("candidati")
            decisione["scelto"] = max(candidati, key=lambda v: jaccard(
                set(misure[v]["varianti"]), set(corrente["varianti"])))
        decisioni[parametro] = decisione

    if modalita:
        attuale = modalita_corrente(config, modalita)
        ammissibili, somiglianza = {}, {}
        for nome, parametri in modalita.items():
            etichetta = "modalita_" + "".join(c if c.isalnum() else "_" for c in nome)
            if nome == attuale:
                m = corrente
            else:
                print(f"{nome_modalita}: {nome}", flush=True)
                m = esegui_variante(etichetta, base, origine, lavoro, _per_gruppo(parametri),
                                    opzioni.rifai)
            ammissibili[nome] = m["ammissibile"]
            somiglianza[nome] = (jaccard(set(m["varianti"]), set(corrente["varianti"]))
                                 if m["ammissibile"] else 0.0)
            righe.append(_riga(nome_modalita, nome, m, corrente,
                               "corrente" if nome == attuale else ""))
            sospeso = VINCOLO_SOSPESO.get(m.get("codice_arresto"))
            if sospeso:
                # Diagnostica: la stessa modalita' con il vincolo sospeso, per
                # misurarne l'effetto a valle. Non e' una configurazione candidata.
                cambi = _per_gruppo(parametri)
                for gruppo, valori in sospeso.items():
                    cambi.setdefault(gruppo, {}).update(valori)
                diagnostica = esegui_variante(etichetta + "_vincolo_sospeso", base, origine,
                                              lavoro, cambi, opzioni.rifai)
                righe.append(_riga(
                    nome_modalita, nome, diagnostica, corrente,
                    "diagnostica: " + ", ".join(f"{g}.{k}" for g, v in sospeso.items() for k in v)
                    + " sospeso, non candidata"))
        decisioni[nome_modalita] = {
            **scegli_modalita(attuale, ammissibili, somiglianza), "ammissibili": ammissibili,
        }

    intestazione = [
        "parametro", "valore", "esito", "campioni_conservati", "campioni_che_cambiano",
        "varianti_finali", "jaccard_con_la_corrente", "letture_conservate",
        "letture_rispetto_alla_corrente", "esclusi_per_profondita", "varianti_rimosse_prevalenza",
        "contaminanti_rimossi", "letture_biologici_rimosse_come_contaminanti", "modello_soglia",
        "piastre_senza_curva_propria", "soglie_per_piastra", "nota",
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
    try:
        sys.exit(main())
    except ErroreSensibilita as guasto:
        sys.exit(f"analisi di sensibilita' interrotta: {guasto}")
