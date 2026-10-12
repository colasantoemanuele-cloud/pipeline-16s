"""Il sottocomando ``soglie``: dalle misure di un'esecuzione alle soglie.

Uso, dalla radice di un clone e nel container::

    python3 scripts/esegui.py soglie --config configurazione.yaml \\
        --output nuova.yaml [--out CARTELLA] [--accetta-proposte] [--json]

Legge la configurazione usata e la cartella di un'esecuzione conclusa, senza
ricalcolare nulla e senza modificarla, e scrive un **nuovo** file di
configurazione: lo stesso testo, salvo le soglie che fermano o avvisano
(:mod:`amplicon16s.soglie.tabella`), ognuna con la misura osservata e il valore
proposto. Senza ``--accetta-proposte`` le righe proposte sono in commento, e il
parametro manca finche' chi legge non toglie il commento; con l'opzione sono
attive. Lo strumento propone, l'utente conferma: l'opzione e' il suo gesto.

**La regola.** Per una soglia massima il suggerimento e' la misura piu' il
margine, per una minima la misura meno il margine; non e' mai piu' stretto della
misura, cosi' che la ripresa sugli stessi dati non avvisi ne' si fermi per una
soglia appena scritta. Se la soglia della configurazione e' attiva, e' rispettata
dalla misura e non e' piu' larga del suggerimento, si tiene quella. Se il valore
finale coincide con il limite della soglia il controllo resta spento, e il
comando lo dichiara.

**Che cosa una soglia suggerita non fa.** Non giudica se la misura e' buona: la
protegge dai peggioramenti (campioni aggiunti, una scelta cambiata, lo
strumento aggiornato). Per questo ogni soglia e' riportata accanto alla sua
misura, sulla console, nel documento JSON e nel file.

**Codici di uscita**: ``0`` file scritto; ``3`` rifiuto (file gia' esistente,
configurazione non valida o diversa da quella dell'esecuzione, esecuzione non
conclusa, artefatto mancante, file non riscrivibile): nulla e' stato scritto.
"""

from __future__ import annotations

import argparse
import copy
import json
import math
import os
import re
import statistics
import textwrap
from collections.abc import Callable, Mapping
from dataclasses import dataclass, replace
from decimal import Decimal
from fractions import Fraction
from pathlib import Path
from typing import Any, Final

import yaml

from amplicon16s import interfaccia
from amplicon16s.config.resolve import risolvi
from amplicon16s.config.schema import Config, ErroreConfigurazione, carica, valida
from amplicon16s.init.scrittura import valore_yaml
from amplicon16s.io_layer.artifacts import nome_manifesto_passo, scrivi_atomico
from amplicon16s.runner.graph import GRAFO
from amplicon16s.soglie.misure import (
    ArtefattoNonLeggibile,
    LettureFinali,
    Misura,
    campioni_controllati,
    gruppi_dei_campioni_finali,
    leggi,
    letture_finali,
)
from amplicon16s.steps.base import impronta_parametri
from amplicon16s.config import defaults
from amplicon16s.soglie.tabella import (
    LETTURE_FINALI_POCHE,
    QUOTA_DELLA_MEDIANA_DELLE_LETTURE,
    SOGLIE_SENZA_LIMITE_AMMESSO,
    TRONCAMENTO_DELLO_SCARTO,
    Direzione,
    Soglia,
    Unita,
    controlli_spenti,
    limite,
    oltre_il_limite,
    soglie_di,
    valori_piatti,
)

__all__ = ["USCITA_RIFIUTO", "Esito", "LettureMinime", "Rifiuto", "aggiungi_argomenti",
           "esegui", "letture_minime", "proponi", "riscrivi", "suggerisci"]

USCITA_RIFIUTO: Final = 3

#: Le frazioni proposte si scrivono con quattro decimali, arrotondate verso
#: l'esterno: per eccesso le massime, per difetto le minime.
_DECIMALI: Final = 10_000

#: Gli esiti possibili per una soglia.
PROPOSTA: Final = "proposta"
TENUTA: Final = "tenuta"
SENZA_MISURA: Final = "senza_misura"
NON_LEGGIBILE: Final = "esclusa"
NON_PROPOSTA: Final = "non_proposta"
OMESSA: Final = "omessa"

_AVVERTENZA: Final = (
    "Una soglia ricavata dalla prima esecuzione non giudica se la misura e' buona: la "
    "protegge dai peggioramenti (campioni aggiunti, una scelta cambiata, lo strumento "
    "aggiornato). Il giudizio sulla misura resta a chi conduce l'analisi: guarda ogni misura "
    "accanto al valore proposto."
)


class Rifiuto(Exception):
    """``soglie`` non scrive alcun file, e dice perche'."""


@dataclass(frozen=True)
class Esito:
    """Che cosa ``soglie`` dice di una soglia."""

    voce: Soglia
    esito: str
    misura: Misura | None
    #: Il valore della soglia nella configurazione usata (``None`` se omessa).
    configurazione: Any
    #: Misura piu' o meno il margine, riportata nell'intervallo ammesso.
    suggerito: Any = None
    #: Il valore che il file nuovo dichiara: il suggerimento, o quello tenuto.
    valore: Any = None
    regola: str = ""
    spento: bool = False
    note: tuple[str, ...] = ()

    @property
    def chiave(self) -> str:
        return self.voce.chiave

    @property
    def cambia(self) -> bool:
        """Vero se il file nuovo porta per questa soglia una riga proposta."""
        return self.esito == PROPOSTA

    def documento(self) -> dict[str, Any]:
        """La soglia nella forma del documento JSON."""
        misura = self.misura
        return {
            "parametro": self.chiave,
            "esito": self.esito,
            "direzione": self.voce.direzione.value,
            "unita": self.voce.unita.value,
            "ferma": self.voce.ferma,
            "codici": list(self.voce.codici),
            "misura": None if misura is None or misura.valore is None else _numero(
                misura.valore, self.voce),
            "statistica": misura.statistica if misura else None,
            "origine_della_misura": misura.origine if misura else None,
            "valore_configurazione": self.configurazione,
            "valore_suggerito": self.suggerito,
            "valore": self.valore,
            "margine": self.voce.testo_margine,
            "regola": self.regola,
            "controllo_spento": self.spento,
            "note": list(self.note),
        }


def _numero(valore: Fraction, voce: Soglia) -> float | int:
    """Un valore esatto nel tipo del parametro; una frazione con sei decimali,
    quanti ne scrivono gli artefatti piu' precisi.
    """
    return round(float(valore), 6) if voce.unita is Unita.FRAZIONE else int(valore)


def _esatto(valore: Any) -> Fraction:
    return Fraction(Decimal(str(valore)))


def _verso_l_esterno(valore: Fraction, voce: Soglia) -> Fraction:
    """Il valore nella risoluzione con cui si scrive, senza renderlo piu' stretto."""
    if voce.unita is not Unita.FRAZIONE:
        return Fraction(math.ceil(valore))
    arrotonda = math.ceil if voce.direzione is Direzione.MASSIMA else math.floor
    return Fraction(arrotonda(valore * _DECIMALI), _DECIMALI)


def proponi(voce: Soglia, misura: Misura, configurazione: Any, al_limite: int | None) -> Esito:
    """Applica la regola del suggerimento a una soglia.

    ``configurazione`` e' il valore della soglia nella configurazione usata,
    ``al_limite`` il valore oltre il quale il controllo non puo' scattare
    (``None`` se non e' definito).
    """
    if configurazione is None:
        return Esito(voce, OMESSA, misura, None,
                     regola=("facoltativa e non dichiarata: resta omessa, e il suo controllo "
                             "non si esegue"
                             if voce.chiave in defaults.FACOLTATIVI_ASSENTI_SE_OMESSI else
                             "condizionale non pertinente: resta omessa"))
    if voce.non_proposta is not None:
        return Esito(voce, NON_PROPOSTA, misura, configurazione, valore=configurazione,
                     regola=f"non proposta: {voce.non_proposta}",
                     note=(misura.assente,) if misura.assente else ())
    if misura.valore is None:
        if misura.non_leggibile is not None:
            return Esito(voce, NON_LEGGIBILE, misura, configurazione, valore=configurazione,
                         regola="esclusa dai suggerimenti: la misura non si legge "
                                "dall'esecuzione; la soglia resta com'e'",
                         note=(misura.non_leggibile,))
        return Esito(voce, SENZA_MISURA, misura, configurazione, valore=configurazione,
                     regola="senza misura in questa esecuzione: la soglia resta com'e'",
                     note=(misura.assente or "",))

    massima = voce.direzione is Direzione.MASSIMA
    segno = 1 if massima else -1
    # La misura peggiorata della tolleranza di lettura: il valore vero non e'
    # peggiore di cosi'.
    peggiore = misura.valore + segno * misura.tolleranza
    usata = _esatto(configurazione)
    usata_spenta = oltre_il_limite(voce, configurazione, al_limite)
    note: list[str] = []
    # Un controllo che ferma, in un'esecuzione arrivata in fondo, non e'
    # scattato: il valore vero rispetta la soglia usata. Se questa cade nella
    # tolleranza di lettura, la misura peggiore possibile e' la soglia stessa.
    letta_rispettata = usata >= misura.valore if massima else usata <= misura.valore
    if voce.ferma and not usata_spenta and letta_rispettata:
        peggiore = min(peggiore, usata) if massima else max(peggiore, usata)
    margine = Fraction(voce.margine)
    suggerito = _verso_l_esterno(peggiore + segno * margine, voce)
    if al_limite is not None and oltre_il_limite(voce, suggerito, al_limite):
        if suggerito != al_limite:
            note.append(f"misura e margine escono dall'intervallo ammesso: riportata al "
                        f"limite, {al_limite}")
        suggerito = Fraction(al_limite)
    verso = "piu'" if massima else "meno"
    regola = f"misura {verso} il margine ({voce.testo_margine})"

    basso, alto = (peggiore, suggerito) if massima else (suggerito, peggiore)
    if not usata_spenta and basso <= usata <= alto:
        return Esito(
            voce, TENUTA, misura, configurazione, suggerito=_numero(suggerito, voce),
            valore=configurazione,
            regola=f"tenuta la soglia della configurazione: e' rispettata dalla misura e non "
                   f"e' piu' larga del suggerimento ({regola})",
            spento=False)
    if usata_spenta:
        note.append(f"nella configurazione usata il controllo era spento ({configurazione})")
    elif (usata < peggiore) if massima else (usata > peggiore):
        note.append(f"la soglia della configurazione usata ({configurazione}) era gia' piu' "
                    "stretta della misura"
                    + ("" if voce.ferma else ": il controllo ha avvisato"))
    valore = _numero(suggerito, voce)
    spento = oltre_il_limite(voce, valore, al_limite)
    if spento and voce.chiave in SOGLIE_SENZA_LIMITE_AMMESSO:
        ancora = (suggerito < peggiore) if massima else (suggerito > peggiore)
        note.append(
            "CONTROLLO SPENTO: il valore proposto e' il massimo che lo schema ammette "
            f"({TRONCAMENTO_DELLO_SCARTO[voce.chiave]} meno 1)" + (
                f"; la misura ({_numero(misura.valore, voce)}) lo supera: l'avviso "
                f"{', '.join(voce.codici)} continuera' a scattare, perche' lo schema non "
                "ammette un valore piu' largo" if ancora else ""))
    elif spento:
        note.append("CONTROLLO SPENTO: il valore proposto coincide con il limite oltre il "
                    "quale il controllo non puo' scattare")
    return Esito(voce, PROPOSTA, misura, configurazione, suggerito=valore, valore=valore,
                 regola=regola, spento=spento, note=tuple(note))


def suggerisci(
    valori: Mapping[str, Any], misure: Mapping[str, Misura], controllati: int | None
) -> list[Esito]:
    """L'esito di ogni soglia, nell'ordine della tabella.

    Avviso e arresto delle chimere si calcolano ciascuno con il suo margine;
    poi l'avviso non supera mai l'arresto: se lo supererebbe (l'arresto della
    configurazione e' stato tenuto, ed e' sotto la misura piu' il margine
    dell'avviso) viene riportato al valore dell'arresto.
    """
    esiti = {
        voce.chiave: proponi(voce, misure[voce.chiave], valori.get(voce.chiave),
                             limite(voce, valori, controllati))
        for voce in soglie_di(valori)
    }
    # La quota di coppie fuse, in modalita' appaiata: soglie minime, quindi
    # l'avviso non scende mai sotto l'arresto.
    if "pair.min_merged_frac_warn" in esiti:
        avviso, arresto = esiti["pair.min_merged_frac_warn"], esiti["pair.min_merged_frac_stop"]
        if (avviso.valore is not None and arresto.valore is not None
                and avviso.valore < arresto.valore):
            esiti[avviso.chiave] = replace(
                avviso, esito=PROPOSTA, valore=arresto.valore, spento=arresto.spento,
                note=(*avviso.note,
                      f"riportato al valore dell'arresto ({arresto.valore}): l'avviso non "
                      "puo' stare sotto pair.min_merged_frac_stop"))
    avviso, arresto = esiti["qc.warn_frac_chimeric"], esiti["qc.stop_frac_chimeric"]
    if avviso.valore is not None and arresto.valore is not None and avviso.valore > arresto.valore:
        esiti[avviso.chiave] = replace(
            avviso, esito=PROPOSTA, valore=arresto.valore, spento=arresto.spento,
            note=(*avviso.note,
                  f"riportato al valore dell'arresto ({arresto.valore}): l'avviso non puo' "
                  "superare qc.stop_frac_chimeric"))
    # Spenta e' ogni soglia che, con il valore che il file nuovo dichiara, non
    # puo' far scattare il suo controllo: anche una rimasta com'era.
    finali = {**valori, **{c: e.valore for c, e in esiti.items() if e.esito != OMESSA}}
    spente = {s.chiave for s in controlli_spenti(finali, controllati)}
    return [replace(esiti[voce.chiave], spento=voce.chiave in spente)
            for voce in soglie_di(valori)]


# --------------------------------------------------------------------------- #
# Il file                                                                      #
# --------------------------------------------------------------------------- #

_GRUPPO: Final = re.compile(r"^([A-Za-z_]\w*):\s*(#.*)?$")
_MARCATORE: Final = re.compile(r"^\s*## \[[A-Z ]+\]")
_NOSTRO: Final = "[SOGLIE]"
_LARGHEZZA: Final = 78


def _commento(testo: str, rientro: str) -> list[str]:
    return [f"{rientro}## {riga}" for riga in textwrap.wrap(
        testo, _LARGHEZZA - len(rientro) - 3, break_long_words=False, break_on_hyphens=False)]


def _testo_misura(esito: Esito) -> str:
    misura = esito.misura
    if misura is None or misura.valore is None:
        return "nessuna misura"
    return f"misura {_numero(misura.valore, esito.voce)} ({misura.statistica}; {misura.origine})"


def _righe_di_commento(esito: Esito, attive: bool, rientro: str) -> list[str]:
    """Il commento che precede la riga di una soglia nel file nuovo."""
    if esito.cambia:
        testo = (f"[PROPOSTO] da soglie{', accettato con --accetta-proposte' if attive else ''}: "
                 f"{_testo_misura(esito)}; regola: {esito.regola}; valore nella configurazione "
                 f"di partenza: {valore_yaml(esito.configurazione)}.")
    elif esito.esito == TENUTA:
        testo = (f"{_NOSTRO} {esito.regola}; {_testo_misura(esito)}; suggerimento "
                 f"{valore_yaml(esito.suggerito)}.")
    elif esito.esito == NON_PROPOSTA:
        testo = f"{_NOSTRO} {esito.regola}; {_testo_misura(esito)}."
    else:
        testo = f"{_NOSTRO} {esito.regola}."
    righe = _commento(testo, rientro)
    for nota in esito.note:
        if nota:
            righe += _commento(nota, rientro)
    return righe


def riscrivi(testo: str, esiti: list[Esito], attive: bool) -> str:
    """Il testo della configurazione con le soglie proposte al posto delle loro righe.

    Ogni soglia con una proposta prende il suo valore, in una riga attiva o in
    commento, preceduta dal commento con origine, misura e regola; le altre
    conservano la riga e ricevono un commento che ne dichiara l'esito. Il
    resto del file non cambia. Solleva :class:`Rifiuto` se una riga non si
    trova: il file non ha la forma a gruppi e righe che il comando sa riscrivere.
    """
    righe = testo.split("\n")
    per_gruppo: dict[str, dict[str, Esito]] = {}
    for esito in esiti:
        if esito.esito != OMESSA:
            gruppo, nome = esito.chiave.split(".")
            per_gruppo.setdefault(gruppo, {})[nome] = esito

    uscita: list[str] = []
    trovate: set[str] = set()
    gruppo: str | None = None
    for riga in righe:
        intestazione = _GRUPPO.match(riga)
        if intestazione:
            gruppo = intestazione.group(1)
        elif riga and not riga[0].isspace() and not riga.startswith("#"):
            gruppo = None
        esito = None
        if gruppo in per_gruppo and riga[:1].isspace():
            trovata = re.match(r"^(\s+)([A-Za-z_]\w*):(\s|$)", riga)
            if trovata and trovata.group(2) in per_gruppo[gruppo]:
                esito = per_gruppo[gruppo][trovata.group(2)]
        if esito is None:
            uscita.append(riga)
            continue
        rientro, nome = trovata.group(1), trovata.group(2)
        trovate.add(esito.chiave)
        # Le righe di commento che precedono: da un marcatore in poi dicono da
        # dove veniva il valore di prima. Una proposta le sostituisce tutte; gli
        # altri esiti tolgono solo il commento di un passaggio precedente.
        inizio = len(uscita)
        while inizio and uscita[inizio - 1].lstrip().startswith("##"):
            inizio -= 1
        marcatori = [i for i in range(inizio, len(uscita)) if _MARCATORE.match(uscita[i])]
        if esito.cambia and marcatori:
            del uscita[marcatori[0]:]
        elif marcatori and _NOSTRO in uscita[marcatori[-1]]:
            del uscita[marcatori[-1]:]
        uscita += _righe_di_commento(esito, attive, rientro)
        if esito.cambia:
            in_linea = re.search(r"\s#\s*(.*\S)\s*$", riga)
            if in_linea:
                # Il commento che l'autore aveva messo accanto al valore di prima.
                uscita += _commento(
                    f"commento sulla riga di partenza: {in_linea.group(1)}", rientro)
            uscita.append(f"{rientro}{'' if attive else '# '}{nome}: {valore_yaml(esito.valore)}")
        else:
            uscita.append(riga)
    mancanti = sorted(e.chiave for e in esiti if e.esito != OMESSA and e.chiave not in trovate)
    if mancanti:
        raise Rifiuto(
            "nel file di configurazione non trovo la riga di " + ", ".join(mancanti)
            + ": il comando riscrive solo file con un gruppo per blocco e un parametro per riga")
    return "\n".join(uscita)


def _intestazione(esiti: list[Esito], attive: bool, cartella: Path) -> str:
    """Le righe in testa al file nuovo: che cosa e' stato fatto, e l'avvertenza."""
    proposte = [e for e in esiti if e.cambia]
    spente = [e.chiave for e in esiti if e.spento]
    altre = [e for e in esiti if e.esito in (SENZA_MISURA, NON_LEGGIBILE)]
    barra = "# " + "=" * (_LARGHEZZA - 2)

    def blocco(testo: str) -> list[str]:
        return [f"# {r}" for r in textwrap.wrap(testo, _LARGHEZZA - 2, break_long_words=False,
                                                break_on_hyphens=False)]

    righe = [barra, "# Soglie proposte da soglie, dalle misure di un'esecuzione conclusa", barra]
    righe += blocco(f"Esecuzione letta: {cartella}.")
    righe += blocco(
        f"Soglie proposte: {len(proposte)} ("
        + ("righe attive, accettate con --accetta-proposte" if attive else
           "righe in commento: togli il commento per accettare, o correggi il valore; finche' "
           "restano in commento validate le elenca come mancanti")
        + "): " + (", ".join(e.chiave for e in proposte) or "nessuna") + ".")
    tenute = [e.chiave for e in esiti if e.esito == TENUTA]
    if tenute:
        righe += blocco("Tenute come sono: " + ", ".join(tenute) + ".")
    if altre:
        righe += blocco("Senza misura o escluse, lasciate come sono: "
                        + ", ".join(e.chiave for e in altre) + ".")
    if spente:
        righe += blocco("CONTROLLI SPENTI anche dopo la proposta: " + ", ".join(spente) + ".")
    righe += blocco(_AVVERTENZA)
    righe += blocco("Le intestazioni che seguono descrivono il file com'era prima di questo "
                    "passaggio: per le soglie vale cio' che e' scritto qui e accanto a "
                    "ciascuna riga.")
    righe += ["#", ""]
    return "\n".join(righe)


def _senza_intestazione(testo: str) -> str:
    """Il testo senza l'intestazione di un passaggio precedente del comando."""
    righe = testo.split("\n")
    if len(righe) > 2 and righe[1].startswith("# Soglie proposte da soglie"):
        # L'intestazione finisce con la sua sola riga di cancelletto vuoto.
        for i, riga in enumerate(righe):
            if riga == "#":
                return "\n".join(righe[i + 1:])
    return testo


def _dati(testo: str) -> dict[str, Any]:
    dati = yaml.safe_load(testo)
    if not isinstance(dati, dict):
        raise Rifiuto("il file di configurazione non contiene una mappa di gruppi")
    return dati


def _verifica(originale: str, nuovo: str, accettato: str, esiti: list[Esito]) -> None:
    """Il file nuovo, riletto, ha gli stessi valori di quello di partenza salvo
    le soglie proposte; con le proposte attive e' una configurazione valida.
    """
    attesi = _dati(originale)
    senza = copy.deepcopy(attesi)
    con = copy.deepcopy(attesi)
    for esito in esiti:
        if esito.cambia:
            gruppo, nome = esito.chiave.split(".")
            senza[gruppo].pop(nome)
            con[gruppo][nome] = esito.valore
    try:
        if _dati(accettato) != con or _dati(nuovo) not in (con, senza):
            raise Rifiuto(
                "il file riscritto, riletto, non ha gli stessi valori del file di partenza "
                "salvo le soglie: il comando non sa riscrivere questo file")
        valida(_dati(accettato), origine="configurazione con le soglie proposte")
    except yaml.YAMLError as guasto:
        raise Rifiuto(f"il file riscritto non e' YAML leggibile: {guasto}") from guasto
    except ErroreConfigurazione as guasto:
        raise Rifiuto("con le soglie proposte la configurazione non sarebbe valida: "
                      + "; ".join(guasto.problemi)) from guasto


# --------------------------------------------------------------------------- #
# Letture minime per campione                                                  #
# --------------------------------------------------------------------------- #

#: Il parametro che toglie dall'oggetto finale i campioni con poche letture.
_LETTURE_MINIME: Final = "qc.min_reads_final"
#: Il marcatore del commento che ``soglie`` scrive accanto a quella riga.
_DA_SCEGLIERE: Final = "[DA SCEGLIERE] da soglie, letture minime per campione:"


@dataclass(frozen=True)
class LettureMinime:
    """La proposta di ``qc.min_reads_final``: non una soglia di qualita', ma
    una scelta che cambia l'oggetto finale, e che per questo resta a chi
    conduce l'analisi.
    """

    #: Il valore della configurazione usata.
    configurazione: int
    #: Il valore con cui i campioni elencati escono.
    raccomandato: int
    #: La mediana delle letture finali dei campioni rimasti.
    mediana: float
    #: Quanti campioni ha l'oggetto finale.
    campioni: int
    #: I campioni con molte meno letture degli altri, dal piu' povero.
    poveri: tuple[LettureFinali, ...]
    #: Vero se toglierli e' raccomandato: hanno meno letture di quante ne
    #: servono a descrivere un campione. Falso se ne hanno molte meno degli
    #: altri ma abbastanza in assoluto: si nominano, e la scelta resta aperta.
    consigliato: bool = True
    #: I gruppi dello studio che, togliendoli, resterebbero con meno di due
    #: campioni: (colonna, valore, campioni che restano, campioni di prima).
    gruppi_svuotati: tuple[tuple[str, str, int, int], ...] = ()

    def testo(self) -> str:
        """La proposta a parole, per la console e per il file."""
        elenco = ", ".join(
            f"{c.nome or c.accession} ({c.letture})" for c in self.poveri[:10]) + (
            f" e altri {len(self.poveri) - 10}" if len(self.poveri) > 10 else "")
        inizio = (
            f"{len(self.poveri)} campioni su {self.campioni} dell'oggetto finale hanno molte "
            f"meno letture degli altri (meno di {self.raccomandato}; mediana "
            f"{self.mediana:.0f}): {elenco}. Le analisi che rarefanno alla profondita' del "
            f"campione piu' povero scenderebbero a {self.poveri[0].letture} letture per tutti. ")
        if self.consigliato:
            consiglio = (
                f"Raccomandato {_LETTURE_MINIME}: {self.raccomandato}, con cui i campioni "
                "elencati escono dall'oggetto finale")
        else:
            consiglio = (
                f"Hanno pero' almeno {LETTURE_FINALI_POCHE} letture, che bastano a "
                "descriverli: di solito si tengono, e non c'e' un valore raccomandato. Per "
                f"toglierli comunque: {_LETTURE_MINIME}: {self.raccomandato}")
        gruppi = ""
        if self.gruppi_svuotati:
            gruppi = ". ATTENZIONE: togliendoli resterebbero con meno di due campioni " + "; ".join(
                f"il gruppo {valore!r} di {colonna} ({restano} su {prima})"
                for colonna, valore, restano, prima in self.gruppi_svuotati[:6]
            ) + ": un confronto che lo riguarda non si potrebbe piu' fare"
        return (inizio + consiglio + gruppi + ". Non e' una soglia di qualita': cambia i "
                "risultati, e per questo non si attiva con --accetta-proposte")

    def documento(self) -> dict[str, Any]:
        """La proposta nella forma del documento JSON."""
        return {
            "parametro": _LETTURE_MINIME,
            "valore_configurazione": self.configurazione,
            "valore_raccomandato": self.raccomandato,
            "raccomandato": self.consigliato,
            "mediana": self.mediana,
            "campioni": self.campioni,
            "campioni_che_uscirebbero": [
                {"accession": c.accession, "nome": c.nome, "letture": c.letture}
                for c in self.poveri],
            "gruppi_con_meno_di_due_campioni": [
                {"colonna": colonna, "valore": valore, "restano": restano, "campioni": prima}
                for colonna, valore, restano, prima in self.gruppi_svuotati],
        }


def letture_minime(
    campioni: list[LettureFinali] | None, configurazione: Any,
    gruppi: Mapping[str, Mapping[str, str]] | None = None,
) -> LettureMinime | None:
    """La proposta di ``qc.min_reads_final``, o ``None`` se nessun campione
    dell'oggetto finale ha molte meno letture degli altri.

    Poveri sono i campioni sotto :data:`LETTURE_FINALI_POCHE` letture: toglierli
    e' raccomandato, perche' non bastano a descrivere un campione. Se nessuno e'
    sotto quel numero, ma alcuni sono sotto un decimo della mediana, si
    nominano senza raccomandare di toglierli: hanno abbastanza letture, e
    possono essere un intero gruppo dello studio. Se sotto il limite cadrebbe
    almeno la meta' dei campioni il dataset e' povero nel suo insieme: vale
    allora il solo decimo della mediana, e se anche cosi' sono almeno la meta'
    non c'e' proposta. ``gruppi`` da', per ogni colonna dello studio, il valore
    di ogni campione: serve a dire quali gruppi resterebbero con meno di due.
    """
    if not campioni or len(campioni) < 2 or not isinstance(configurazione, int):
        return None
    mediana = statistics.median(c.letture for c in campioni)
    decimo = math.ceil(Fraction(mediana) * Fraction(QUOTA_DELLA_MEDIANA_DELLE_LETTURE))
    scelto: tuple[int, list[LettureFinali], bool] | None = None
    for limite_inferiore, consigliato in (
            (LETTURE_FINALI_POCHE, True), (decimo, decimo < LETTURE_FINALI_POCHE)):
        poveri = sorted((c for c in campioni if c.letture < limite_inferiore),
                        key=lambda c: (c.letture, c.accession))
        if poveri and 2 * len(poveri) < len(campioni):
            scelto = (limite_inferiore, poveri, consigliato)
            break
    if scelto is None:
        return None
    limite_inferiore, poveri, consigliato = scelto
    fuori = {c.accession for c in poveri}
    presenti = {c.accession for c in campioni}
    svuotati = []
    for colonna, valori in (gruppi or {}).items():
        per_valore: dict[str, list[str]] = {}
        for accession, valore in valori.items():
            if accession in presenti and valore:
                per_valore.setdefault(valore, []).append(accession)
        for valore, membri in sorted(per_valore.items()):
            restano = sum(a not in fuori for a in membri)
            if len(membri) >= 2 and restano < 2:
                svuotati.append((colonna, valore, restano, len(membri)))
    return LettureMinime(configurazione, limite_inferiore, float(mediana), len(campioni),
                         tuple(poveri), consigliato, tuple(svuotati))


def _con_le_letture_minime(testo: str, proposta: LettureMinime | None) -> str:
    """Il testo con la proposta scritta sopra la riga di ``qc.min_reads_final``.

    Toglie il commento di un passaggio precedente del comando; senza una
    proposta, e senza quel commento, il testo resta com'e'. La riga del
    parametro non cambia: chi accetta scrive il valore raccomandato al posto
    del suo.
    """
    gruppo_cercato, nome = _LETTURE_MINIME.split(".")
    uscita: list[str] = []
    gruppo: str | None = None
    for riga in testo.split("\n"):
        intestazione = _GRUPPO.match(riga)
        if intestazione:
            gruppo = intestazione.group(1)
        trovata = re.match(rf"^(\s+){nome}:(\s|$)", riga) if gruppo == gruppo_cercato else None
        if trovata:
            inizio = len(uscita)
            while inizio and uscita[inizio - 1].lstrip().startswith("##"):
                inizio -= 1
            nostre = [i for i in range(inizio, len(uscita)) if _DA_SCEGLIERE in uscita[i]]
            if nostre:
                del uscita[nostre[0]:]
            if proposta is not None:
                uscita += _commento(
                    f"{_DA_SCEGLIERE} {proposta.testo()}. Per toglierli scrivi "
                    f"{proposta.raccomandato} al posto del valore nella riga qui sotto, poi "
                    "resume.", trovata.group(1))
        uscita.append(riga)
    return "\n".join(uscita)


# --------------------------------------------------------------------------- #
# L'esecuzione letta                                                           #
# --------------------------------------------------------------------------- #


def _esecuzione_conclusa(radice: Path, config: Config) -> None:
    """Rifiuta se la cartella non contiene un'esecuzione conclusa con questa
    configurazione, coerente da una fase all'altra.
    """
    stato = interfaccia.stato(radice)
    if stato["stato"] != "conclusa":
        mancanti = [f["passo"] for f in stato["fasi"] if f["stato"] in ("da_eseguire", "in_corso")]
        dettaglio = {
            "non_avviata": "non contiene un'esecuzione",
            "in_corso": "l'esecuzione e' in corso",
            "interrotta": "l'esecuzione e' stata interrotta: riprendila con resume",
            "ferma": "l'esecuzione si e' fermata con un codice ("
                     + str((stato.get("arresto") or {}).get("codice")) + "): correggi la causa "
                     "e riprendila con resume",
            "incompleta": "l'esecuzione non e' arrivata in fondo",
        }[stato["stato"]]
        raise Rifiuto(
            f"la cartella {radice} non contiene un'esecuzione conclusa: {dettaglio}"
            + (f" (fasi da eseguire: {', '.join(mancanti)})" if mancanti else "")
            + ". Le soglie si propongono dalle misure di un'esecuzione arrivata in fondo")
    registrata = (stato.get("configurazione") or {}).get("digest")
    digest = risolvi(config).digest
    if registrata != digest:
        raise Rifiuto(
            f"la configurazione indicata ({digest}) non e' quella registrata dall'esecuzione "
            f"di {radice} ({registrata}): le misure non sarebbero le sue. Indica il file "
            "con cui l'esecuzione e' stata fatta, dagli stessi percorsi")
    # Ogni fase deve essere calcolata sulle fasi a monte oggi su disco: una
    # fase rifatta da sola (per esempio da validate con un'altra
    # configurazione) lascia le successive calcolate su artefatti superati.
    manifesti: dict[str, dict[str, Any]] = {}
    for passo in GRAFO.ordine():
        percorso = radice / GRAFO.nodo(passo).cartella.value / nome_manifesto_passo(str(passo))
        if percorso.is_file():
            try:
                manifesto = json.loads(percorso.read_text(encoding="utf-8"))
            except (OSError, ValueError) as guasto:
                raise Rifiuto(f"il manifesto della fase {passo} non si legge: {guasto}") from guasto
            if not isinstance(manifesto, dict):
                raise Rifiuto(f"il manifesto della fase {passo} non ha la forma attesa")
            manifesti[str(passo)] = manifesto
    risolta = risolvi(config)
    for passo, manifesto in manifesti.items():
        calcolata = manifesto.get("calcolata_su")
        calcolata = calcolata if isinstance(calcolata, dict) else {}
        # La fase deve essere stata calcolata con i valori che la configurazione
        # indicata da' ai parametri da cui dipende: una configurazione
        # registrata dopo, da una validazione che non ha rifatto la fase, ha lo
        # stesso digest dell'ultima registrata ma non e' quella delle misure.
        dei_parametri = calcolata.get("configurazione")
        dei_parametri = dei_parametri if isinstance(dei_parametri, dict) else {}
        dichiarati = dei_parametri.get("parametri")
        try:
            attesa = impronta_parametri(
                risolta, None if dichiarati is None else tuple(dichiarati))
        except (KeyError, TypeError, ValueError):
            attesa = None
        if attesa is None or dei_parametri.get("impronta") != attesa:
            raise Rifiuto(
                f"la fase {passo} non e' stata calcolata con i valori che la configurazione "
                "indicata da' ai suoi parametri: le sue misure non sarebbero quelle di questa "
                "configurazione. Lancia resume con questa configurazione, poi soglie")
        a_monte = calcolata.get("a_monte") or {}
        for monte, impronta in a_monte.items():
            if monte in manifesti and manifesti[monte].get("impronta") != impronta:
                raise Rifiuto(
                    f"la fase {passo} e' calcolata su artefatti di {monte} che non sono piu' "
                    "quelli su disco: l'esecuzione non e' coerente. Lancia resume, poi soglie")


# --------------------------------------------------------------------------- #
# Il comando                                                                   #
# --------------------------------------------------------------------------- #


def aggiungi_argomenti(parser: argparse.ArgumentParser) -> None:
    """Gli argomenti del sottocomando."""
    parser.add_argument("--config", "-c", required=True, type=Path, metavar="FILE",
                        help="la configurazione con cui l'esecuzione e' stata fatta")
    parser.add_argument("--out", type=Path, metavar="CARTELLA",
                        help="la cartella di uscita dell'esecuzione conclusa (senza: "
                             "io.out_root della configurazione)")
    parser.add_argument("--output", required=True, type=Path, metavar="FILE",
                        help="il file di configurazione da scrivere: non deve esistere")
    parser.add_argument("--accetta-proposte", action="store_true",
                        help="scrive le soglie proposte come righe attive (senza: in commento)")
    parser.add_argument("--json", action="store_true",
                        help="scrive sullo standard output l'esito come documento JSON")


def _documento(args: argparse.Namespace, cartella: Path | None, esiti: list[Esito],
               uscita: int, rifiuto: str | None,
               minime: LettureMinime | None = None) -> dict[str, Any]:
    return {
        # Solo se c'e' una proposta: senza, il documento e' quello di sempre.
        **({"letture_minime": minime.documento()} if minime is not None else {}),
        "versione_formato": interfaccia.VERSIONE_FORMATO,
        "comando": "soglie",
        "strumento": "amplicon16s",
        "configurazione": str(args.config),
        "cartella": str(cartella) if cartella is not None else None,
        "file_scritto": str(args.output) if uscita == 0 else None,
        "proposte_accettate": bool(args.accetta_proposte),
        "codice_uscita": uscita,
        "rifiuto": rifiuto,
        "soglie": [e.documento() for e in esiti],
        "controlli_spenti": [e.chiave for e in esiti if e.spento],
        "avvertenza": _AVVERTENZA if esiti else None,
    }


def _riga_console(esito: Esito) -> list[str]:
    testa = f"  {esito.chiave}: "
    if esito.esito == OMESSA:
        return [testa + "omessa (condizionale non pertinente): resta omessa"]
    if esito.esito in (SENZA_MISURA, NON_LEGGIBILE):
        che = "senza misura" if esito.esito == SENZA_MISURA else "esclusa"
        return [testa + f"{che}, resta {valore_yaml(esito.configurazione)}",
                *[f"      {n}" for n in esito.note if n]]
    if esito.esito == NON_PROPOSTA:
        misurata = "" if esito.misura is None or esito.misura.valore is None else (
            f"misura {_numero(esito.misura.valore, esito.voce)}, ")
        return [testa + f"{misurata}resta {valore_yaml(esito.configurazione)}"
                + (" [CONTROLLO SPENTO]" if esito.spento else ""),
                f"      {esito.regola}", *[f"      {n}" for n in esito.note if n]]
    misura = esito.misura
    assert misura is not None and misura.valore is not None
    righe = [
        testa + f"misura {_numero(misura.valore, esito.voce)}, configurazione "
        f"{valore_yaml(esito.configurazione)}, "
        + (f"proposta {valore_yaml(esito.valore)}" if esito.cambia else
           f"tenuta (suggerimento {valore_yaml(esito.suggerito)})")
        + (" [CONTROLLO SPENTO]" if esito.spento else ""),
        f"      {misura.statistica}; letta da {misura.origine}",
        f"      regola: {esito.regola}",
    ]
    return righe + [f"      {n}" for n in esito.note if n]


def prepara(args: argparse.Namespace, in_cartella: Callable[[Config], Config]) -> tuple[
    Path, list[Esito], str, LettureMinime | None
]:
    """Legge esecuzione e configurazione e restituisce la cartella letta, gli
    esiti, il testo del file nuovo e la proposta delle letture minime per
    campione, se c'e'. Solleva :class:`Rifiuto`.
    """
    uscita: Path = args.output
    if uscita.exists():
        raise Rifiuto(f"{uscita} esiste gia': soglie non scrive su un file esistente. Indica "
                      "con --output un file nuovo")
    if not uscita.parent.is_dir():
        raise Rifiuto(f"la cartella {uscita.parent} non esiste")
    if not os.access(uscita.parent, os.W_OK | os.X_OK):
        raise Rifiuto(f"la cartella {uscita.parent} non e' scrivibile dall'utente che esegue")
    try:
        config = in_cartella(carica(args.config))
        grezzo = Path(args.config).read_bytes()
        testo = grezzo.decode("utf-8").replace("\r\n", "\n")
    except ErroreConfigurazione as guasto:
        raise Rifiuto("la configurazione non e' valida: " + "; ".join(guasto.problemi)) from guasto
    except (OSError, UnicodeDecodeError) as guasto:
        raise Rifiuto(f"la configurazione non si legge: {guasto}") from guasto
    radice = Path(args.out) if args.out is not None else Path(config.io.out_root)
    _esecuzione_conclusa(radice, config)

    valori = valori_piatti(config)
    try:
        misure = leggi(radice, valori)
    except ArtefattoNonLeggibile as guasto:
        raise Rifiuto(f"un artefatto dell'esecuzione manca o non si legge ({guasto}): la "
                      "cartella non e' quella di un'esecuzione conclusa integra") from guasto
    esiti = suggerisci(valori, misure, campioni_controllati(radice))

    minime = letture_minime(letture_finali(radice), valori.get(_LETTURE_MINIME),
                            gruppi_dei_campioni_finali(radice))

    base = _senza_intestazione(testo)
    nuovo = _con_le_letture_minime(riscrivi(base, esiti, bool(args.accetta_proposte)), minime)
    _verifica(testo, nuovo, riscrivi(base, esiti, True), esiti)
    scritto = _intestazione(esiti, bool(args.accetta_proposte), radice) + nuovo
    # Il file nuovo ha il fine riga di quello di partenza.
    return (radice, esiti, scritto.replace("\n", "\r\n") if b"\r\n" in grezzo else scritto,
            minime)


def esegui(
    args: argparse.Namespace,
    in_cartella: Callable[[Config], Config],
    stampa: Callable[[str], None] = print,
    stampa_json: Callable[[dict[str, Any]], None] | None = None,
) -> int:
    """Esegue ``soglie`` e restituisce il codice di uscita.

    ``in_cartella`` rende assoluti i percorsi della configurazione come fanno
    gli altri comandi, perche' il digest sia quello registrato dall'esecuzione.
    """
    try:
        radice, esiti, testo, minime = prepara(args, in_cartella)
    except Rifiuto as rifiuto:
        stampa(f"soglie: {rifiuto}. Nulla e' stato scritto")
        if args.json and stampa_json is not None:
            stampa_json(_documento(args, None, [], USCITA_RIFIUTO, str(rifiuto)))
        return USCITA_RIFIUTO
    scrivi_atomico(args.output, testo)

    proposte = [e for e in esiti if e.cambia]
    stampa(f"Soglie dalle misure dell'esecuzione in {radice}")
    for esito in esiti:
        for riga in _riga_console(esito):
            stampa(riga)
    stampa(f"Configurazione scritta in {args.output}: {len(proposte)} soglie proposte, "
           + ("righe attive (--accetta-proposte)." if args.accetta_proposte else
              "righe in commento: togli il commento per accettarle."))
    spente = [e.chiave for e in esiti if e.spento]
    if spente:
        stampa("Avviso: controlli spenti anche dopo la proposta: " + ", ".join(spente) + ".")
    for riga in textwrap.wrap(_AVVERTENZA, 96):
        stampa(riga)
    if minime is not None:
        for riga in textwrap.wrap(
                f"Letture minime per campione, da scegliere: {minime.testo()}. Per toglierli "
                f"scrivi {minime.raccomandato} nella riga di {_LETTURE_MINIME} in "
                f"{args.output}, prima di resume.", 96):
            stampa(riga)
    stampa(f"Poi: resume --config {args.output}" if args.accetta_proposte or not proposte else
           f"Poi: togli il commento alle soglie che accetti in {args.output} (finche' restano "
           "in commento validate le elenca come mancanti), e lancia resume --config "
           f"{args.output}")
    if args.json and stampa_json is not None:
        stampa_json(_documento(args, radice, esiti, 0, None, minime))
    return 0
