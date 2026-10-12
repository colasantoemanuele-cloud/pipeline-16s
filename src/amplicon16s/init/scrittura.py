"""Il file di configurazione scritto da ``init``: tre specie di valore.

Ogni riga dichiara da dove viene il suo valore, e la specie decide se la riga
e' attiva:

* **rilevato**: un fatto verificato sui file. Riga attiva;
* **indicato**: un percorso dato sulla riga di comando. Riga attiva;
* **proposto**: l'esito di un'euristica. Riga in commento con il valore gia'
  scritto: per accettarlo si toglie il commento;
* **da scegliere**: una scelta di chi conduce l'analisi. Riga in commento, con
  un valore di partenza e il suo motivo, oppure senza valore.

I parametri che hanno un predefinito seguono il modello di configurazione: il
predefinito in un commento, con la sua specie e la fonte, senza riga attiva.
Il file e' quindi una configurazione valida solo quando chi lo legge ha tolto
il commento a cio' che accetta e ha scritto cio' che manca: nessun valore vi
entra senza un suo gesto.

Il gesto puo' essere uno solo per tutto il file: con ``--accetta-proposte``
le proposte e i valori di partenza sono righe attive, ciascuna con la sua
origine nel commento. Compaiono allora due specie in piu':

* **partenza**: il valore di partenza di una scelta (il predefinito del
  metodo, o un valore ricavato dai dati), accettato. Riga attiva;
* **esplorativa**: una soglia di qualita' senza valore di partenza, al valore
  che non ferma: il suo controllo e' spento. Serve alla prima esecuzione; poi
  il sottocomando ``soglie`` propone il valore dalle misure. Riga attiva.

Razionale sistemistico: le descrizioni dei parametri sono quelle dello schema
(:mod:`amplicon16s.config.schema`), le stesse che ``schema`` esporta: una sola
fonte per il modulo di un programma esterno, per questo file e per la guida.
"""

from __future__ import annotations

import json
import textwrap
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from enum import StrEnum
from typing import Any, Final

from amplicon16s.config import defaults
from amplicon16s.config.schema import Config

__all__ = ["SENZA_VALORE", "Specie", "Voce", "testo_yaml", "valore_yaml"]


class Specie(StrEnum):
    """Da dove viene il valore di una riga."""

    RILEVATO = "RILEVATO"
    INDICATO = "INDICATO"
    PROPOSTO = "PROPOSTO"
    DA_SCEGLIERE = "DA SCEGLIERE"
    CONDIZIONALE = "CONDIZIONALE"
    PARTENZA = "PARTENZA"
    ESPLORATIVA = "ESPLORATIVA"

    @property
    def attiva(self) -> bool:
        """Vero se la riga entra nella configurazione senza altri gesti: un
        fatto, un percorso indicato, o un valore accettato con l'opzione.
        """
        return self in (Specie.RILEVATO, Specie.INDICATO, Specie.PARTENZA, Specie.ESPLORATIVA)


class _SenzaValore:
    """Il segnaposto di una riga senza valore (``None`` e' un valore: nullo)."""

    def __repr__(self) -> str:
        return "SENZA_VALORE"


SENZA_VALORE: Final = _SenzaValore()


@dataclass(frozen=True)
class Voce:
    """Un parametro nel file scritto da ``init``."""

    chiave: str
    specie: Specie
    valore: Any
    #: Da dove viene il valore, o perche' non ce n'e' uno: finisce nel commento.
    origine: str
    #: Altre righe di commento: alternative, misure, che cosa guardare.
    note: tuple[str, ...] = ()
    #: Vero per una proposta accettata con l'opzione: la riga e' attiva.
    accettata: bool = False

    @property
    def attiva(self) -> bool:
        """Vero se la riga del parametro e' scritta attiva."""
        return self.specie.attiva or self.accettata


#: I gruppi dello schema riservati alle analisi a valle, senza parametri: non
#: si scrivono.
_GRUPPI_RISERVATI: Final = ("norm", "glom", "beta", "ord", "stat")

_LARGHEZZA: Final = 78


def valore_yaml(valore: Any) -> str:
    """Il valore in una forma YAML su una riga, senza ambiguita'.

    Le stringhe sono sempre fra virgolette doppie, con le sequenze di escape di
    JSON, che YAML legge allo stesso modo: una colonna con le parentesi quadre,
    un'espressione regolare o un numero scritto come testo restano cio' che sono.
    """
    if valore is None:
        return "null"
    if isinstance(valore, bool):
        return "true" if valore else "false"
    if isinstance(valore, int):
        return repr(valore)
    if isinstance(valore, float):
        # YAML legge come numero solo una mantissa con il punto e un esponente
        # con il segno: "1e-40" sarebbe una stringa.
        testo = repr(valore)
        if "e" in testo:
            mantissa, esponente = testo.split("e")
            if "." not in mantissa:
                mantissa += ".0"
            if esponente[0] not in "+-":
                esponente = "+" + esponente
            testo = f"{mantissa}e{esponente}"
        return testo
    if isinstance(valore, (list, tuple)):
        return "[" + ", ".join(valore_yaml(x) for x in valore) + "]"
    return json.dumps(str(valore), ensure_ascii=False)


def _commento(testo: str, rientro: str = "  ", segno: str = "##") -> list[str]:
    """Il testo in righe di commento, mandate a capo.

    Le spiegazioni dentro i gruppi portano un doppio cancelletto, le righe dei
    parametri in commento uno solo: chi toglie i commenti ai parametri, a mano
    o con uno strumento, non scambia una spiegazione per un parametro.
    """
    return [
        f"{rientro}{segno} {riga}" for riga in textwrap.wrap(
            testo, _LARGHEZZA - len(rientro) - len(segno) - 1, break_long_words=False,
            break_on_hyphens=False)
    ]


#: Nel file le lettere accentate si scrivono con l'apostrofo, come nel modello
#: di configurazione.
_SENZA_ACCENTI: Final = str.maketrans({
    "à": "a'", "è": "e'", "é": "e'", "ì": "i'", "ò": "o'", "ù": "u'", "È": "E'",
})


def _primo_capoverso(testo: str | None) -> str:
    """Il primo capoverso di una docstring, su una riga."""
    capoverso = (testo or "").strip().split("\n\n")[0]
    return " ".join(capoverso.split()).translate(_SENZA_ACCENTI)


_MARCATORI: Final = {
    "unica_forma": "[UNICA FORMA] la pipeline realizza questo solo valore",
    "facoltativo": "[FACOLTATIVO] il predefinito e' l'assenza",
    "strumento": "[STRUMENTO] governa come si esegue, non che cosa si conclude dai dati",
    "modalita": (
        "[MODALITA'] omesso vale naive_bayes. idtaxa e' disponibile: richiede un "
        "classificatore addestrato nel formato di DECIPHER (tax.idtaxa_trained con il suo "
        "MD5, tax.idtaxa_threshold, tax.idtaxa_strand) al posto del riferimento FASTA. "
        "dati/riferimento/scarica_addestrato.py scarica quelli di SILVA: 138.2 per "
        "un'analisi nuova, 2019 per confrontarsi con tabelle elaborate con SILVA r138"),
}

#: I parametri che esistono solo con una modalita' che init non ha proposto:
#: non si scrivono, come il gruppo ``pair`` in un file single-end.
_DI_UNA_MODALITA: Final = frozenset((
    "filter.primer_removal", *defaults.PARAMETRI_DEL_TAGLIO_CERCATO,
    *defaults.PARAMETRI_DI_IDTAXA,
    # I facoltativi che omessi non esistono: si scrivono solo se init li tratta.
    *defaults.FACOLTATIVI_ASSENTI_SE_OMESSI))


def _marcatori_dei_predefiniti() -> dict[str, str]:
    """Per ogni parametro con un predefinito, il marcatore della sua specie:
    la stessa classificazione che ``schema`` esporta.
    """
    from amplicon16s.interfaccia import schema

    return {
        p["nome"]: (f"[STANDARD] {p['fonte']}" if p["specie"] == "standard"
                    else _MARCATORI[p["specie"]])
        for p in schema()["parametri"]
        if not p["obbligatorio"] and p["specie"] not in ("condizionale", "obbligatorio")
    }


def _guida_alle_specie(accettate: bool) -> list[str]:
    """Le righe d'intestazione che spiegano i marcatori."""
    if not accettate:
        return [
            "#   [PROPOSTO]      l'esito di un'euristica, che puo' sbagliare: riga in",
            "#                   commento con il valore gia' scritto. Per accettarlo,",
            "#                   togli il commento; altrimenti correggilo.",
            "#   [DA SCEGLIERE]  una scelta di chi conduce l'analisi: riga in commento,",
            "#                   con un valore di partenza e il suo motivo, o senza valore.",
        ]
    return [
        "#   [PROPOSTO]      l'esito di un'euristica, che puo' sbagliare: riga attiva,",
        "#                   accettata con --accetta-proposte. Controllala comunque.",
        "#   [PARTENZA]      il valore di partenza di una scelta (il predefinito del",
        "#                   metodo, o un valore ricavato dai dati), accettato con",
        "#                   --accetta-proposte: riga attiva.",
        "#   [ESPLORATIVA]   una soglia di qualita' al valore che non ferma: riga attiva,",
        "#                   ma il controllo e' spento. Serve alla prima esecuzione;",
        "#                   poi il sottocomando soglie propone il valore dalle misure.",
        "#   [DA SCEGLIERE]  una scelta di chi conduce l'analisi, senza un valore di",
        "#                   partenza: riga in commento. Va scritta.",
    ]


def testo_yaml(
    voci: Mapping[str, Voce], intestazione: Sequence[str], non_indovinati: Sequence[str],
    accettate: bool = False,
) -> str:
    """Il testo del file di configurazione.

    ``voci`` da' la riga di ogni parametro che ``init`` tratta; gli altri
    parametri dello schema, che hanno tutti un predefinito, sono scritti in
    commento con quello. ``intestazione`` sono le righe sul dataset esaminato,
    ``non_indovinati`` cio' che ``init`` dichiara di non poter stabilire.
    """
    # Nell'ordine dello schema, che e' quello in cui le righe compaiono nel file.
    ordine = [
        f"{gruppo}.{nome}" for gruppo, campo in Config.model_fields.items()
        for nome in campo.annotation.model_fields
    ]
    voci = {chiave: voci[chiave] for chiave in ordine if chiave in voci}
    per_specie = {s: [v for v in voci.values() if v.specie is s] for s in Specie}
    da_scegliere = per_specie[Specie.DA_SCEGLIERE]
    barra = "# " + "=" * (_LARGHEZZA - 2)
    righe = [
        barra,
        "# Configurazione di partenza scritta da init: da rivedere prima di eseguire"
        if not accettate else
        "# Configurazione di partenza scritta da init con --accetta-proposte",
        barra,
        "#",
        *[f"# {r}" for r in intestazione],
        "#",
        "# Ogni riga dice da dove viene il suo valore:",
        "#   [RILEVATO]      un fatto verificato sui file: riga attiva.",
        "#   [INDICATO]      un percorso dato sulla riga di comando: riga attiva.",
        *_guida_alle_specie(accettate),
        "#   [CONDIZIONALE]  un parametro che serve solo se una condizione e' vera, e",
        "#                   con i valori di questo file non lo e': riga in commento,",
        "#                   senza valore. Non va scelto; validate lo chiede se la",
        "#                   condizione diventa vera.",
        "# I parametri con un predefinito (standard del metodo, una sola forma,",
        "# dello strumento) sono in commento con quello: si possono lasciare cosi'.",
        "# Le spiegazioni cominciano con ##, le righe dei parametri in commento con #.",
        "#",
        f"# Rilevati: {len(per_specie[Specie.RILEVATO])}; indicati: "
        f"{len(per_specie[Specie.INDICATO])}; proposti: {len(per_specie[Specie.PROPOSTO])}; "
        f"da scegliere: {len(da_scegliere)}; non richiesti: "
        f"{len(per_specie[Specie.CONDIZIONALE])}.",
        *([
            "#",
            f"# Accettati con --accetta-proposte: {len(per_specie[Specie.PROPOSTO])} proposti, "
            f"{len(per_specie[Specie.PARTENZA])} valori di partenza,",
            f"# {len(per_specie[Specie.ESPLORATIVA])} soglie esplorative, al valore che non ferma:",
            *_commento(", ".join(v.chiave for v in per_specie[Specie.ESPLORATIVA]) or "nessuna",
                       "", "#"),
            "#",
            "# Valori di partenza che dichiarano un'assenza (nullo o elenco vuoto): sono",
            "# validi, ma se lo studio ha quell'informazione va scritta.",
            *_commento(", ".join(
                v.chiave for v in per_specie[Specie.PARTENZA] if v.valore in (None, [])
            ) or "nessuno", "", "#"),
        ] if accettate else []),
        "#",
        "# Finche' una riga [PROPOSTO] o [DA SCEGLIERE] resta in commento il parametro",
        "# non e' dichiarato: validate elenca quelli che mancano (E-G15-10).",
        "#",
        "# Da scegliere:",
        *[riga for testo in (
            ", ".join(v.chiave for v in da_scegliere) or "nulla",)
          for riga in _commento(testo, "", "#")],
    ]
    righe += [
        "#",
        "# Non richiesti con le proposte attuali (condizionali):",
        *_commento(", ".join(v.chiave for v in per_specie[Specie.CONDIZIONALE]) or "nessuno",
                   "", "#"),
    ]
    if non_indovinati:
        righe += ["#", "# Che cosa init non stabilisce:"]
        for voce in non_indovinati:
            capo, *resto = textwrap.wrap(voce, _LARGHEZZA - 6, break_long_words=False,
                                         break_on_hyphens=False)
            righe.append(f"#   - {capo}")
            righe += [f"#     {r}" for r in resto]
    righe.append("")

    marcatori = _marcatori_dei_predefiniti()
    for gruppo, campo_gruppo in Config.model_fields.items():
        if gruppo in _GRUPPI_RISERVATI:
            continue
        if gruppo == "pair" and not any(chiave.startswith("pair.") for chiave in voci):
            # Il gruppo della modalita' appaiata non esiste in un file single-end.
            continue
        modello = campo_gruppo.annotation
        righe += _commento(_primo_capoverso(modello.__doc__), "", "#")
        righe.append(f"{gruppo}:")
        for nome, campo in modello.model_fields.items():
            chiave = f"{gruppo}.{nome}"
            voce = voci.get(chiave)
            if voce is None and chiave in _DI_UNA_MODALITA:
                continue
            righe += _commento(campo.description or "")
            if voce is None:
                if campo.is_required() or chiave in defaults.PARAMETRI_DEL_BAYESIANO:
                    # Un percorso obbligatorio non indicato sulla riga di comando.
                    righe += _commento("[DA INDICARE] non indicato sulla riga di comando")
                    righe.append(f"  # {nome}:")
                else:
                    predefinito = campo.get_default(call_default_factory=True)
                    righe += _commento(marcatori[chiave])
                    righe.append(f"  # {nome}: {valore_yaml(predefinito)}")
                continue
            righe += _commento(f"[{voce.specie.value}] {voce.origine}")
            for nota in voce.note:
                righe += _commento(nota)
            scritto = "" if voce.valore is SENZA_VALORE else f" {valore_yaml(voce.valore)}"
            righe.append(f"  {'' if voce.attiva else '# '}{nome}:{scritto}")
        righe.append("")
    return "\n".join(righe)
