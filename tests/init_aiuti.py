"""Aiuti per i test della configurazione assistita.

Il file scritto da ``init`` va completato da una persona: toglie il commento
alle proposte che accetta e scrive le scelte che mancano. :func:`completa` fa
lo stesso gesto su un testo, cosi' i test verificano che il file, una volta
completato, sia una configurazione valida.
"""

from __future__ import annotations

import re
from typing import Any

from amplicon16s.config.schema import Config
from amplicon16s.init.scrittura import valore_yaml

_GRUPPO = re.compile(r"^([a-z]+):\s*$")
_RIGA_IN_COMMENTO = re.compile(r"^  # ([A-Za-z_0-9]+):(.*)$")
_RIGA_ATTIVA = re.compile(r"^  ([A-Za-z_0-9]+):(.*)$")


def _e_un_parametro(gruppo: str | None, nome: str, modello: Any = Config) -> bool:
    """Vero se ``gruppo.nome`` e' un parametro dello schema: una riga di
    commento in prosa puo' cominciare con una parola seguita da due punti.
    """
    campo = modello.model_fields.get(gruppo or "")
    return campo is not None and nome in campo.annotation.model_fields


def righe_dei_parametri(testo: str, modello: Any = Config) -> dict[str, tuple[bool, str]]:
    """Per ogni parametro scritto nel file, se la riga e' attiva e il valore
    che porta (testo dopo i due punti, vuoto se senza valore).
    """
    trovate: dict[str, tuple[bool, str]] = {}
    gruppo = None
    for riga in testo.splitlines():
        if (m := _GRUPPO.match(riga)):
            gruppo = m.group(1)
        elif (m := _RIGA_ATTIVA.match(riga)) and _e_un_parametro(gruppo, m.group(1), modello):
            trovate[f"{gruppo}.{m.group(1)}"] = (True, m.group(2).strip())
        elif (m := _RIGA_IN_COMMENTO.match(riga)) and _e_un_parametro(gruppo, m.group(1), modello):
            trovate.setdefault(f"{gruppo}.{m.group(1)}", (False, m.group(2).strip()))
    return trovate


def completa(
    testo: str, scelte: dict[str, Any], accetta: set[str] | None = None, modello: Any = Config
) -> str:
    """Il testo del file con le righe in commento completate.

    Per ogni parametro in ``scelte`` la riga diventa attiva con quel valore,
    anche al posto di un valore rilevato. Le altre righe in commento che
    portano un valore diventano attive se il parametro e' in ``accetta``
    (``None``: tutte quelle con un valore): e' accettare la proposta o il
    valore di partenza. Le righe senza valore e non in ``scelte`` restano come sono.
    """
    righe = []
    gruppo = None
    for riga in testo.splitlines():
        if (m := _GRUPPO.match(riga)):
            gruppo = m.group(1)
        attiva, commento = _RIGA_ATTIVA.match(riga), _RIGA_IN_COMMENTO.match(riga)
        trovata = attiva or commento
        if trovata and _e_un_parametro(gruppo, trovata.group(1), modello):
            chiave = f"{gruppo}.{trovata.group(1)}"
            if chiave in scelte:
                riga = f"  {trovata.group(1)}: {valore_yaml(scelte[chiave])}"
            elif commento and commento.group(2).strip() and (accetta is None or chiave in accetta):
                riga = f"  {commento.group(1)}:{commento.group(2)}"
        righe.append(riga)
    return "\n".join(righe) + "\n"
