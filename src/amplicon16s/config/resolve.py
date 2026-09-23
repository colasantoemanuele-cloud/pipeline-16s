"""Risoluzione dei parametri derivati e registrazione della configurazione usata.

Alcuni parametri non si scrivono a mano: discendono da altri e calcolarli è
l'unico modo di garantire che restino coerenti. Scriverli nella configurazione
significherebbe poterli mettere in contraddizione con i parametri da cui
dipendono, ed è per questo che tentarlo è un errore.

**La risoluzione avviene in due momenti**, perché i derivati non dipendono
tutti dalle stesse informazioni:

1. :func:`risolvi` — risoluzione *statica*. Calcola tutto ciò che discende
   dalla sola configurazione: ``filter.minLen``, ``asv.len_min``,
   ``asv.len_max``. Può essere eseguita prima di toccare qualunque dato, ed è
   quello che serve al gate G15.
2. :meth:`ConfigRisolta.con_campioni_biologici` — risoluzione *dipendente dai
   dati*. Calcola ``prev.min_samples``, che richiede il numero di campioni
   biologici: un dato che non esiste finché i metadati non sono stati letti.

La separazione non è una comodità: senza di essa il gate G15 non potrebbe
girare prima della lettura dei metadati, e l'intero scopo di un controllo che
precede qualunque calcolo verrebbe meno.

**Il digest** identifica la combinazione di parametri impiegata ed è calcolato
sulla configurazione dichiarata più i derivati statici. ``prev.min_samples``
ne è deliberatamente escluso: dipende da quanti campioni biologici contiene il
dataset, cioè descrive i dati e non la configurazione. Includerlo farebbe
cambiare il digest fra la prima e la seconda fase della risoluzione, e due
esecuzioni sulla stessa configurazione non risulterebbero più confrontabili.
"""

from __future__ import annotations

import hashlib
import json
import math
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Any, Final

import yaml

from amplicon16s.config.schema import PARAMETRI_DERIVATI, Config
from amplicon16s.io_layer.artifacts import AlberoOutput, Fase

__all__ = [
    "ConfigRisolta",
    "Derivati",
    "NOME_FILE_RISOLTO",
    "PARAMETRI_DERIVATI",
    "risolvi",
    "scrivi_risolta",
]

#: Nome del file con la configurazione effettivamente usata. La cartella in cui
#: finisce è ``Fase.CONFIG``: la creazione delle cartelle di output passa dal
#: servizio degli artefatti, non da qui.
NOME_FILE_RISOLTO: Final = "resolved.yaml"


@dataclass(frozen=True)
class Derivati:
    """Parametri calcolati, non letti dalla configurazione.

    ``prev_min_samples`` vale ``None`` finché il numero di campioni biologici
    non è noto: è l'unico derivato che dipende dai dati e non dalla sola
    configurazione.
    """

    filter_minLen: int
    asv_len_min: int
    asv_len_max: int
    prev_min_samples: int | None = None

    @property
    def completi(self) -> bool:
        """Vero quando anche i derivati dipendenti dai dati sono stati calcolati."""
        return self.prev_min_samples is not None

    def come_chiavi(self) -> dict[str, Any]:
        """I derivati nella forma ``gruppo.parametro``, per confronti e report."""
        return {
            "filter.minLen": self.filter_minLen,
            "asv.len_min": self.asv_len_min,
            "asv.len_max": self.asv_len_max,
            "prev.min_samples": self.prev_min_samples,
        }


def _deriva_statici(config: Config) -> Derivati:
    """Calcola i derivati che discendono dalla sola configurazione."""
    # La lunghezza minima accettata coincide con quella di troncamento: una
    # lettura piu' corta non ha superato il troncamento.
    minLen = config.filter.truncLen

    # Le ASV attese hanno la lunghezza delle letture troncate, meno quanto e'
    # stato tagliato in testa. La tolleranza allarga l'intervallo nei due versi.
    base = config.filter.truncLen - config.filter.trimLeft
    len_min = base - config.asv.len_tol
    len_max = base + config.asv.len_tol

    return Derivati(filter_minLen=minLen, asv_len_min=len_min, asv_len_max=len_max)


def _digest(config: Config, derivati: Derivati) -> str:
    """Digest della combinazione di parametri impiegata.

    Calcolato su una forma canonica: chiavi ordinate, separatori fissi,
    nessuno spazio. Due configurazioni identiche danno lo stesso digest anche
    se scritte con ordine o formattazione diversi.
    """
    statici = derivati.come_chiavi()
    statici.pop("prev.min_samples")  # dipende dai dati, non dalla configurazione

    canonico = json.dumps(
        {"config": config.model_dump(mode="json"), "derivati": statici},
        sort_keys=True,
        ensure_ascii=False,
        separators=(",", ":"),
    )
    return "sha256:" + hashlib.sha256(canonico.encode("utf-8")).hexdigest()


@dataclass(frozen=True)
class ConfigRisolta:
    """Configurazione validata insieme ai suoi parametri derivati."""

    config: Config
    derivati: Derivati

    @property
    def digest(self) -> str:
        """Identificatore della combinazione di parametri impiegata."""
        return _digest(self.config, self.derivati)

    @property
    def completa(self) -> bool:
        return self.derivati.completi

    def con_campioni_biologici(self, quanti: int) -> ConfigRisolta:
        """Seconda fase: calcola i derivati che dipendono dai dati.

        ``prev.min_samples`` è il numero minimo di campioni in cui una variante
        deve comparire, arrotondato per eccesso perché una soglia frazionaria
        non ha significato: i campioni si contano interi.
        """
        if quanti < 0:
            raise ValueError(
                f"il numero di campioni biologici non puo' essere negativo: {quanti}"
            )

        min_samples = math.ceil(self.config.prev.min_fraction * quanti)
        return replace(self, derivati=replace(self.derivati, prev_min_samples=min_samples))

    def come_mappa(self) -> dict[str, Any]:
        """Configurazione e derivati in una sola mappa, pronta per il file.

        I derivati compaiono dentro il gruppo a cui appartengono, così il file
        mostra i parametri effettivamente in uso e non solo quelli dichiarati.
        """
        mappa = self.config.model_dump(mode="json")
        for chiave, valore in self.derivati.come_chiavi().items():
            gruppo, parametro = chiave.split(".")
            mappa[gruppo][parametro] = valore
        return mappa


def risolvi(config: Config) -> ConfigRisolta:
    """Prima fase: risolve i derivati che non dipendono dai dati."""
    return ConfigRisolta(config=config, derivati=_deriva_statici(config))


_INTESTAZIONE = """\
# =========================================================================== #
# Configurazione effettivamente usata — generata, non da modificare           #
# =========================================================================== #
#
# Registra i parametri con cui l'esecuzione e' stata avviata, compresi quelli
# calcolati dalla pipeline e non dichiarati nel file di ingresso.
#
# Il digest identifica la combinazione di parametri: due esecuzioni con lo
# stesso digest hanno usato la stessa configurazione. E' calcolato sulla
# configurazione dichiarata piu' i derivati statici, in forma canonica
# (chiavi ordinate, separatori fissi).
#
# prev.min_samples e' escluso dal digest: dipende da quanti campioni biologici
# contiene il dataset, quindi descrive i dati e non la configurazione. Vale
# null finche' i metadati non sono stati letti.
#
# Questo file non e' riutilizzabile come configurazione di ingresso: i
# parametri derivati non sono ammessi in ingresso proprio perche' calcolati.
"""


def scrivi_risolta(risolta: ConfigRisolta, out_root: Path | str) -> Path:
    """Scrive la configurazione risolta in ``00_config``.

    Passa dal servizio degli artefatti, quindi il file entra nel manifesto
    della fase con il proprio checksum come qualunque altro artefatto: la
    configurazione usata è un artefatto dell'esecuzione, non un file a parte.
    """
    documento = {
        "digest": risolta.digest,
        "derivati": list(risolta.derivati.come_chiavi()),
        "parametri": risolta.come_mappa(),
    }

    corpo = yaml.safe_dump(
        documento,
        sort_keys=False,
        allow_unicode=True,
        default_flow_style=False,
    )

    albero = AlberoOutput(out_root)
    artefatto = albero.scrivi_testo(Fase.CONFIG, NOME_FILE_RISOLTO, _INTESTAZIONE + corpo)
    return artefatto.percorso
