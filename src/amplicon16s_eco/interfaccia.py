"""Uscite in JSON delle analisi ecologiche, per i programmi esterni.

Gli stessi tre documenti della pipeline (:mod:`amplicon16s.interfaccia`), con
lo stesso campo ``versione_formato``: lo schema della configurazione, l'esito
della validazione e lo stato di una cartella di uscita.

Un'analisi ecologica e' un solo calcolo, senza fasi e senza ripresa: lo stato
di una cartella e' ``non_avviata``, ``in_corso``, ``interrotta`` o
``conclusa``. Un rifiuto non lascia nulla nella cartella di uscita: lo riporta
il comando stesso, con il suo codice, e ``validate --json`` lo anticipa.
"""

from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Any

from amplicon16s.interfaccia import (
    VERSIONE_FORMATO,
    descrivi_parametri,
    documento_schema,
    esecuzione_attiva,
)

from amplicon16s_eco import __version__
from amplicon16s_eco.catalogo import CATALOGO
from amplicon16s_eco.config import CONDIZIONALI, OBBLIGATORI, STANDARD_DEL_METODO, ConfigEco

__all__ = ["esito_validazione", "problema", "schema", "stato"]

STRUMENTO = "amplicon16s_eco"


def schema() -> dict[str, Any]:
    """Lo schema della configurazione delle analisi, per un programma esterno.

    La specie di un parametro e' ``obbligatorio``, ``condizionale``
    (obbligatorio se uno dei metodi di ``richiesto_da`` e' dichiarato, respinto
    altrimenti), ``standard`` (il predefinito e' quello del pacchetto R, con la
    fonte) oppure ``facoltativo``.
    """
    parametri = descrivi_parametri(ConfigEco)
    for voce in parametri:
        nome = voce["nome"]
        if nome in OBBLIGATORI:
            voce["specie"] = "obbligatorio"
        elif nome in CONDIZIONALI:
            elenco, metodi = CONDIZIONALI[nome]
            voce["specie"] = "condizionale"
            voce["richiesto_da"] = {"parametro": elenco, "valori": list(metodi)}
        elif nome in STANDARD_DEL_METODO:
            voce["specie"] = "standard"
            voce["fonte"] = STANDARD_DEL_METODO[nome][1]
        else:
            voce["specie"] = "facoltativo"
    return documento_schema(STRUMENTO, __version__, ConfigEco, parametri)


_NOME = re.compile(r"\b([a-z]+\.[a-z_]+)\b")


def problema(codice: str, dettaglio: str) -> dict[str, Any]:
    """Un rifiuto o un avviso nella forma del documento."""
    noti = {
        f"{gruppo}.{nome}" for gruppo, campo in ConfigEco.model_fields.items()
        for nome in campo.annotation.model_fields
    }
    parametri = [n for n in dict.fromkeys(_NOME.findall(dettaglio)) if n in noti]
    voce = CATALOGO.get(codice)
    return {
        "parametro": parametri[0] if parametri else None,
        "parametri": parametri,
        "codice": codice,
        "categoria": voce.specie.value if voce else None,
        "messaggio": f"{voce.sintesi} {dettaglio}".strip() if voce else dettaglio,
        "azione": voce.azione if voce else None,
    }


def esito_validazione(
    *, configurazione: Path | str, oggetto: Path | str | None, codice_uscita: int,
    problemi: list[dict[str, Any]], avvisi: list[dict[str, Any]] | None = None,
) -> dict[str, Any]:
    """Il documento di ``validate --json``."""
    return {
        "versione_formato": VERSIONE_FORMATO,
        "comando": "validate",
        "strumento": STRUMENTO,
        "configurazione": str(configurazione),
        "oggetto": str(oggetto) if oggetto is not None else None,
        "valida": codice_uscita == 0 and not problemi,
        "codice_uscita": codice_uscita,
        "problemi": problemi,
        "avvisi": avvisi or [],
    }


def _json(percorso: Path) -> Any | None:
    try:
        return json.loads(percorso.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return None


def stato(cartella: Path | str, adesso: float | None = None) -> dict[str, Any]:
    """Lo stato dell'analisi nella cartella di uscita indicata, senza scrivervi."""
    radice = Path(cartella)
    attiva = esecuzione_attiva(radice, adesso)
    manifesto = _json(radice / "manifest.json")
    riepilogo = _json(radice / "riepilogo.json")
    if attiva is not None:
        complessivo = "in_corso" if attiva["battito_recente"] else "interrotta"
    elif isinstance(manifesto, dict):
        complessivo = "conclusa"
    else:
        complessivo = "non_avviata"
    documento: dict[str, Any] = {
        "versione_formato": VERSIONE_FORMATO,
        "comando": "status",
        "strumento": STRUMENTO,
        "cartella": str(radice),
        "stato": complessivo,
        "esecuzione": attiva,
        "configurazione": None,
        "oggetto_di_partenza": None,
        "file": None,
        "avvisi": [],
        # Un rifiuto non lascia nulla nella cartella: lo riporta il comando.
        "arresto": None,
    }
    if isinstance(manifesto, dict) and complessivo == "conclusa":
        documento["configurazione"] = manifesto.get("configurazione")
        documento["oggetto_di_partenza"] = manifesto.get("oggetto_di_partenza")
        # Come il comando run: i file elencati dal manifesto e il manifesto.
        documento["file"] = len(manifesto.get("file") or []) + 1
        if isinstance(riepilogo, dict):
            documento["avvisi"] = [a.get("codice") for a in riepilogo.get("avvisi") or []]
    return documento
