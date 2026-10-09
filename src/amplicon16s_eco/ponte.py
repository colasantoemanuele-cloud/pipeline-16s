"""Esecuzione degli script R del pacchetto, attraverso il ponte della pipeline.

Gli script stanno in ``R/`` accanto a questo modulo e girano come la pipeline
esegue i propri: processo separato, ``Rscript --vanilla``, sessione propria,
``LANGUAGE=en`` e ``LC_ALL=C.UTF-8``, richiesta e dichiarazione d'esito in due
file JSON con l'identificativo dell'invocazione. Di :mod:`amplicon16s.rbridge`
si riusano il contratto dei due file, il lancio del processo e la
classificazione dell'esito; non l'albero degli artefatti e il catalogo della
pipeline, che appartengono alle sue fasi.

Gli script caricano le funzioni condivise della pipeline (scrittura atomica,
dichiarazione dell'esito) dalla cartella indicata da ``AMPLICON16S_R_LIB`` e
quelle del pacchetto da ``AMPLICON16S_ECO_R_LIB``: le due variabili le imposta
questo modulo nel solo processo figlio.
"""

from __future__ import annotations

import os
import uuid
from collections.abc import Mapping
from pathlib import Path
from typing import Any, Final

from amplicon16s.rbridge.payload import (
    DichiarazioneNonValida,
    leggi_dichiarazione,
    scrivi_richiesta,
)
from amplicon16s.rbridge.runner import (
    VARIABILE_LIB_R,
    Condizione,
    _lancia,
    cartella_r,
    classifica,
    trova_rscript,
)

from amplicon16s_eco.catalogo import CATALOGO, ErroreEco, Specie

__all__ = ["CARTELLA_R", "VARIABILE_LIB_ECO", "esegui_r"]

#: Gli script R del pacchetto e, in ``lib``, le loro funzioni condivise.
CARTELLA_R: Final = Path(__file__).resolve().parent / "R"

#: Variabile che indica al processo figlio le funzioni R del pacchetto.
VARIABILE_LIB_ECO: Final = "AMPLICON16S_ECO_R_LIB"

#: Quanto dell'uscita di errore di R entra nel dettaglio di un rifiuto.
_CODA_STDERR: Final = 2000


def esegui_r(script: str, parametri: Mapping[str, Any], cartella: Path, passo: str) -> None:
    """Esegue uno script R del pacchetto nella cartella di lavoro indicata.

    Ritorna se lo script dichiara il successo. Solleva
    :class:`~amplicon16s_eco.catalogo.ErroreEco` con il codice dichiarato dallo
    script, se e' un rifiuto del catalogo, oppure con ``E-ECO-90`` in ogni
    altro caso (interprete assente, errore non previsto, processo morto senza
    dichiarare l'esito), con il messaggio di R nel dettaglio.
    """
    interprete = trova_rscript()
    if interprete is None:
        raise ErroreEco([("E-ECO-90", "Rscript non trovato nel PATH.")])
    lib_pipeline = cartella_r() / "lib"
    percorso = CARTELLA_R / script
    for richiesto in (percorso, lib_pipeline, CARTELLA_R / "lib"):
        if not richiesto.exists():
            raise ErroreEco([("E-ECO-90", f"Non trovato: {richiesto}.")])

    invocazione = uuid.uuid4().hex
    richiesta, percorso_esito = scrivi_richiesta(cartella, invocazione, parametri, passo)
    ambiente = {
        **os.environ,
        "LANGUAGE": "en",
        "LC_ALL": "C.UTF-8",
        VARIABILE_LIB_R: str(lib_pipeline),
        VARIABILE_LIB_ECO: str(CARTELLA_R / "lib"),
    }
    comando = [str(interprete), "--vanilla", str(percorso), str(richiesta)]
    try:
        processo = _lancia(comando, cartella, ambiente, None)
    except OSError as e:
        raise ErroreEco([("E-ECO-90", f"L'interprete {interprete} non si avvia: {e}.")]) from e

    try:
        dichiarazione = leggi_dichiarazione(percorso_esito, invocazione)
    except DichiarazioneNonValida:
        dichiarazione = None
    condizione = classifica(processo.codice_uscita, dichiarazione, processo.stderr)
    if condizione is Condizione.RIUSCITO:
        return

    if condizione is Condizione.ERRORE_DICHIARATO and dichiarazione is not None:
        v = CATALOGO.get(dichiarazione.codice or "")
        if v is not None and v.specie is Specie.RIFIUTO:
            raise ErroreEco([(v.codice, dichiarazione.messaggio)])
    messaggio = dichiarazione.messaggio if dichiarazione and dichiarazione.messaggio else ""
    coda = processo.stderr.strip()[-_CODA_STDERR:]
    dettaglio = (
        f"Script {script}, esito del ponte: {condizione.value}, codice di uscita "
        f"{processo.codice_uscita}. {messaggio or coda}"
    )
    raise ErroreEco([("E-ECO-90", dettaglio.strip())])
