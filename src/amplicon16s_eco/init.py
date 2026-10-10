"""Configurazione assistita delle analisi ecologiche: ``eco init``.

Dall'oggetto finale della pipeline scrive una configurazione di partenza delle
analisi, con le tre specie di valore della pipeline (:mod:`amplicon16s.init`):

* **rilevato**: cio' che l'oggetto contiene (ranghi tassonomici, albero e sua
  radice). Non sono parametri: stanno in testa al file, come fatti;
* **proposto**: l'esito di un'euristica, in commento con il valore gia' scritto
  (la variabile biologica, le variabili tecniche, il rango);
* **da scegliere**: le scelte di analisi, in commento, con un valore di
  partenza solo dove c'e' il predefinito di un pacchetto.

Il file non e' eseguibile cosi' com'e': ``validate`` nomina i parametri che
mancano. L'oggetto si legge da una copia e non viene modificato.

Con ``--accetta-proposte`` le proposte e i valori di partenza sono righe
attive, con la loro origine nel commento. Restano in commento, con il loro
valore, i parametri condizionali: sono richiesti da un metodo e respinti senza,
e quali metodi usare e' una scelta che ``init`` non fa.
"""

from __future__ import annotations

import json
import tempfile
from dataclasses import replace
import textwrap
from pathlib import Path
from typing import Any, Final

from amplicon16s.init.scrittura import SENZA_VALORE, Specie, Voce, valore_yaml
from amplicon16s.io_layer.artifacts import scrivi_atomico
from amplicon16s.metadata.crosswalk import _VALORI_ASSENTI
from amplicon16s.metadata.tabelle import COLONNE_INVENTARIO

from amplicon16s_eco.config import CONDIZIONALI, FACOLTATIVI, STANDARD_DEL_METODO, ConfigEco
from amplicon16s_eco.esecuzione import _copia, _impronta_oggetto
from amplicon16s_eco.ponte import esegui_r

__all__ = ["QUOTA_MINIMA_ASSEGNATA", "accetta", "descrivi", "prepara", "scrivi", "testo_yaml"]

#: La quota di varianti assegnate sotto la quale un rango non si propone: sotto
#: la meta', la categoria «non assegnato» sarebbe la piu' grande del grafico.
QUOTA_MINIMA_ASSEGNATA: Final = 0.5

#: Le colonne dell'inventario della pipeline che sono variabili tecniche (la
#: piastra di estrazione e la corsa di sequenziamento, le ultime due), e quelle
#: che identificano il campione o ne dicono la classe (le prime tre), che non
#: sono variabili da confrontare. I nomi sono quelli del contratto
#: dell'oggetto, non scritti qui.
_TECNICHE: Final = COLONNE_INVENTARIO[-2:]
_NON_VARIABILI: Final = COLONNE_INVENTARIO[:3]

PROPOSTO, SCELTA = Specie.PROPOSTO, Specie.DA_SCEGLIERE


def descrivi(oggetto: Path | str) -> dict[str, Any]:
    """Che cosa contiene l'oggetto: campioni, ranghi, albero, colonne.

    Solleva :class:`~amplicon16s_eco.catalogo.ErroreEco` se l'oggetto non e'
    utilizzabile (``E-ECO-04``) o se R non si esegue (``E-ECO-90``).
    """
    oggetto = Path(oggetto)
    impronta = _impronta_oggetto(oggetto)
    with tempfile.TemporaryDirectory(prefix="amplicon16s-eco-") as temporanea:
        lavoro = Path(temporanea)
        copia = _copia(oggetto, lavoro, impronta)
        esegui_r("descrivi.R", {"oggetto": str(copia),
                                "valori_mancanti": sorted(_VALORI_ASSENTI)}, lavoro, "descrizione")
        return json.loads((lavoro / "descrizione.json").read_text(encoding="utf-8"))


def prepara(descrizione: dict[str, Any]) -> tuple[dict[str, Voce], list[str]]:
    """Le voci della configurazione e i fatti rilevati sull'oggetto."""
    voci: dict[str, Voce] = {}

    def metti(chiave: str, specie: Specie, valore: Any, origine: str, *note: str) -> None:
        voci[chiave] = Voce(chiave, specie, valore, origine, tuple(n for n in note if n))

    colonne = descrizione["colonne"]
    ranghi = descrizione["ranghi"]
    albero = descrizione["albero"]
    fatti = [
        f"Oggetto: {descrizione['campioni']} campioni, {descrizione['varianti']} varianti.",
        "Ranghi tassonomici (quota di varianti assegnate): "
        + ", ".join(f"{r['rango']} {r['quota_assegnata']:.2f}" for r in ranghi) + ".",
        "Albero filogenetico: " + (
            ("presente, " + ("radicato" if albero["radicato"] else "senza radice")
             + ("" if albero["con_lunghezze"] else ", senza lunghezze dei rami"))
            if albero["presente"] else "assente (le distanze UniFrac non sono disponibili)") + ".",
    ]

    # Variabile biologica: proposta solo se una sola colonna ha gruppi ripetuti.
    # Il numero di gruppi non dice quale colonna lo studio confronta: una
    # coordinata di laboratorio si ripete quanto un fattore sperimentale.
    candidate = sorted(
        (c for c in colonne
         if c["nome"] not in _TECNICHE + _NON_VARIABILI and c["gruppi_ripetuti"] >= 2),
        key=lambda c: c["nome"])
    if len(candidate) == 1:
        sola = candidate[0]
        metti("design.variable", PROPOSTO, sola["nome"],
              f"la sola colonna con almeno due gruppi di almeno due campioni "
              f"({sola['gruppi_ripetuti']} su {sola['distinti']} valori)",
              *(["attenzione: i valori sono tutti numerici e quasi tutti diversi; le analisi "
                 "la tratterebbero come categorica"] if _sembra_continua(sola) else []))
    elif candidate:
        metti("design.variable", SCELTA, SENZA_VALORE,
              "piu' colonne hanno almeno due gruppi di almeno due campioni, e init non sa "
              "quale lo studio confronta: il numero di gruppi non lo dice",
              "colonne con gruppi ripetuti (numero di gruppi): "
              + ", ".join(f"{c['nome']} ({c['gruppi_ripetuti']})" for c in candidate),
              *(["attenzione: " + ", ".join(c["nome"] for c in candidate if _sembra_continua(c))
                 + " ha valori tutti numerici e quasi tutti diversi; le analisi la "
                 "tratterebbero come categorica"]
                if any(_sembra_continua(c) for c in candidate) else []))
    else:
        metti("design.variable", SCELTA, SENZA_VALORE,
              "nessuna colonna ha almeno due gruppi di almeno due campioni: i test fra "
              "gruppi non si potranno eseguire, restano le uscite descrittive",
              "colonne dell'oggetto: " + ", ".join(c["nome"] for c in colonne))
    tecniche = sorted(
        (c for c in colonne if c["nome"] in _TECNICHE and c["distinti"] > 1),
        key=lambda c: c["distinti"])
    metti("design.technical_variables", PROPOSTO, [c["nome"] for c in tecniche],
          ("le colonne di corsa e di piastra dell'inventario che hanno piu' di un valore, "
           "dalla piu' grossolana alla piu' annidata (meno valori prima)") if tecniche else
          "ne' la corsa ne' la piastra hanno piu' di un valore fra i campioni")

    # Dal rango piu' generale in giu', finche' la quota regge: un rango fine
    # non si propone attraverso uno piu' generale che non la raggiunge.
    proponibili = []
    for rango in ranghi:
        if rango["quota_assegnata"] < QUOTA_MINIMA_ASSEGNATA:
            break
        proponibili.append(rango)
    if proponibili:
        rango = proponibili[-1]
        metti("comp.rank", PROPOSTO, rango["rango"],
              f"il rango piu' fine a cui almeno il {QUOTA_MINIMA_ASSEGNATA:.0%} delle "
              f"varianti ha un'assegnazione ({rango['quota_assegnata']:.2f} delle varianti; "
              "la quota e' sul numero delle varianti, non sulle letture)")
    else:
        metti("comp.rank", SCELTA, SENZA_VALORE,
              f"nessun rango ha almeno il {QUOTA_MINIMA_ASSEGNATA:.0%} delle varianti assegnate")
    if albero["presente"] and albero["radicato"] and albero["con_lunghezze"]:
        metti("beta.unifrac_root", PROPOSTO, "existing",
              "l'albero dell'oggetto e' gia' radicato; serve solo se beta.distances "
              "comprende una UniFrac, ed e' respinto altrimenti")

    metti("run.seed", SCELTA, SENZA_VALORE,
          "un intero qualunque, da dichiarare: fissa la rarefazione e le permutazioni")
    metti("comp.top_n", SCELTA, SENZA_VALORE,
          "quanti taxa mostrare nei grafici prima di riunire gli altri")
    metti("beta.distances", SCELTA, SENZA_VALORE,
          "fra bray, jaccard, aitchison" + (
              ", unifrac_weighted, unifrac_unweighted" if albero["presente"] else
              " (le due UniFrac richiedono un albero, che l'oggetto non ha)"))
    metti("beta.clr_pseudocount", SCELTA, SENZA_VALORE,
          "solo con aitchison: lo pseudoconteggio sommato ai conteggi prima del CLR")
    metti("ord.methods", SCELTA, SENZA_VALORE, "fra pcoa e nmds; [] per nessuna ordinazione")
    metti("ord.distances", SCELTA, SENZA_VALORE, "fra quelle di beta.distances")
    metti("ord.pcoa_correction", SCELTA, "none",
          "solo con pcoa. Partenza: il predefinito di ape::pcoa (correction = \"none\"); "
          "cailliez o lingoes correggono gli autovalori negativi delle distanze non euclidee")
    metti("ord.nmds_trymax", SCELTA, 20,
          "solo con nmds. Partenza: il predefinito di vegan::metaMDS (trymax = 20)")
    metti("stat.min_group_size", SCELTA, SENZA_VALORE,
          "quanti campioni deve avere un gruppo per entrare nei test (almeno 2)")
    metti("stat.significance_level", SCELTA, SENZA_VALORE,
          "il livello sotto il quale le dispersioni dei gruppi si dichiarano diverse")
    return voci, fatti


def accetta(voci: dict[str, Voce]) -> dict[str, Voce]:
    """Le voci con proposte e valori di partenza accettati in blocco.

    Un condizionale non si attiva: e' respinto se il metodo che lo richiede
    non e' fra quelli scelti, e gli elenchi dei metodi non hanno un valore di
    partenza. Resta in commento con il suo valore, e lo dice.
    """
    accettate = {}
    for chiave, voce in voci.items():
        if voce.valore is SENZA_VALORE:
            pass
        elif chiave in CONDIZIONALI:
            elenco, metodi = CONDIZIONALI[chiave]
            voce = replace(voce, note=(
                *voce.note,
                f"non attivato da --accetta-proposte: va dichiarato solo se {elenco} "
                f"comprende {' o '.join(metodi)}, ed e' respinto altrimenti"))
        elif voce.specie is PROPOSTO:
            voce = replace(voce, accettata=True,
                           origine=voce.origine + "; accettato con --accetta-proposte")
        else:
            voce = replace(voce, specie=Specie.PARTENZA,
                           origine=voce.origine + "; accettato con --accetta-proposte")
        accettate[chiave] = voce
    return accettate


def _sembra_continua(colonna: dict[str, Any]) -> bool:
    """La stessa regola dell'avviso sulle variabili che sembrano continue."""
    return bool(colonna["numerica"]) and colonna["distinti"] >= 10 and (
        colonna["distinti"] >= 0.8 * colonna["con_valore"])


def testo_yaml(
    voci: dict[str, Voce], fatti: list[str], oggetto: str, accettate: bool = False
) -> str:
    """Il testo del file di configurazione delle analisi."""
    da_scegliere = [v for v in voci.values() if v.specie is SCELTA]
    proposti = [v for v in voci.values() if v.specie is PROPOSTO]
    attive = [v for v in voci.values() if v.attiva]

    def commento(testo: str, rientro: str = "  ##") -> list[str]:
        # Dentro i gruppi le spiegazioni portano ##, le righe dei parametri in
        # commento un solo #: non si scambiano fra loro.
        return [f"{rientro} {r}" for r in textwrap.wrap(
            testo, 77 - len(rientro), break_long_words=False, break_on_hyphens=False)]

    righe = [
        "# " + "=" * 76,
        "# Configurazione di partenza delle analisi ecologiche, scritta da eco init",
        "# " + "=" * 76,
        "#",
        f"# Oggetto esaminato: {oggetto}",
        "# [RILEVATO] sull'oggetto:",
        *[riga for fatto in fatti for riga in commento(fatto, "#  ")],
        "#",
        *([
            "#   [PROPOSTO]      l'esito di un'euristica: riga in commento con il valore",
            "#                   gia' scritto. Per accettarlo, togli il commento.",
            "#   [DA SCEGLIERE]  una scelta di analisi: riga in commento, con un valore di",
            "#                   partenza solo dove c'e' il predefinito di un pacchetto.",
        ] if not accettate else [
            "#   [PROPOSTO]      l'esito di un'euristica: riga attiva, accettata con",
            "#                   --accetta-proposte. Controllala comunque.",
            "#   [PARTENZA]      il valore di partenza di una scelta, accettato con",
            "#                   --accetta-proposte: riga attiva.",
            "#   [DA SCEGLIERE]  una scelta di analisi: riga in commento, da scrivere. Un",
            "#                   condizionale resta in commento anche se ha un valore.",
        ]),
        "#   [CONDIZIONALE]  va dichiarato solo se il metodo indicato e' richiesto.",
        "#",
        f"# Proposti: {len(proposti)}; da scegliere: {len(da_scegliere)}."
        + (f" Righe attive, accettate con --accetta-proposte: {len(attive)}." if accettate
           else ""),
        "# Da scegliere:",
        *commento(", ".join(v.chiave for v in da_scegliere), "#"),
        "# Le spiegazioni cominciano con ##, le righe dei parametri in commento con #.",
        "",
    ]
    for gruppo, campo_gruppo in ConfigEco.model_fields.items():
        modello = campo_gruppo.annotation
        righe.append(f"# {(modello.__doc__ or '').strip().splitlines()[0]}")
        righe.append(f"{gruppo}:")
        for nome, campo in modello.model_fields.items():
            chiave = f"{gruppo}.{nome}"
            righe += commento(campo.description or "")
            voce = voci.get(chiave)
            if chiave in CONDIZIONALI:
                elenco, metodi = CONDIZIONALI[chiave]
                righe += commento(f"[CONDIZIONALE] solo se {elenco} comprende "
                                  + " o ".join(metodi))
            if voce is None:
                predefinito = campo.get_default(call_default_factory=True)
                if chiave in STANDARD_DEL_METODO:
                    righe += commento(f"[STANDARD] {STANDARD_DEL_METODO[chiave][1]}")
                elif chiave in FACOLTATIVI:
                    righe += commento("[FACOLTATIVO] il predefinito e' l'assenza")
                valore = {} if isinstance(predefinito, dict) else (
                    list(predefinito) if isinstance(predefinito, tuple) else predefinito)
                scritto = "{}" if valore == {} else valore_yaml(valore)
                righe.append(f"  # {nome}: {scritto}")
                continue
            righe += commento(f"[{voce.specie.value}] {voce.origine}")
            for nota in voce.note:
                righe += commento(nota)
            scritto = "" if voce.valore is SENZA_VALORE else f" {valore_yaml(voce.valore)}"
            righe.append(f"  {'' if voce.attiva else '# '}{nome}:{scritto}")
        righe.append("")
    return "\n".join(righe)


def scrivi(
    oggetto: Path | str, uscita: Path | str, accettate: bool = False
) -> tuple[dict[str, Voce], list[str]]:
    """Esamina l'oggetto e scrive la configurazione di partenza.

    Con ``accettate`` le proposte e i valori di partenza sono righe attive
    (:func:`accetta`). Solleva ``FileExistsError`` se il file esiste gia': non
    lo sovrascrive.
    """
    uscita = Path(uscita)
    if uscita.exists():
        raise FileExistsError(uscita)
    voci, fatti = prepara(descrivi(oggetto))
    if accettate:
        voci = accetta(voci)
    scrivi_atomico(uscita, testo_yaml(voci, fatti, Path(oggetto).name, accettate) + "\n")
    return voci, fatti
