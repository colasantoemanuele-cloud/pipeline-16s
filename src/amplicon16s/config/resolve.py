"""Risoluzione dei parametri derivati e registrazione della configurazione usata.

Alcuni parametri non si scrivono a mano: discendono da altri e calcolarli è
l'unico modo di garantire che restino coerenti. Scriverli nella configurazione
significherebbe poterli mettere in contraddizione con i parametri da cui
dipendono, ed è per questo che tentarlo è un errore.

:func:`risolvi` calcola ``asv.len_min`` e ``asv.len_max``, che discendono
dalla sola configurazione: la risoluzione può avvenire prima di toccare qualunque dato, ed è quello che serve al gate G15.
Il numero minimo di campioni del filtro di prevalenza non è un derivato: S13 lo
calcola sui campioni biologici che tiene, un dato che esiste solo dentro la
fase (``R/13_filtri.R``).

**Digest e impronta dei risultati sono due cose diverse.** Il digest
identifica la configurazione ed è calcolato su tutta, compresi i parametri che
non toccano i risultati. L'impronta dei risultati
(:attr:`ConfigRisolta.impronta_risultati`) serve a decidere se una fase
conclusa è ancora valida, e ne esclude i parametri elencati in
:data:`PARAMETRI_SENZA_EFFETTO`: altrimenti cambiare il numero di thread, o
spostare la cartella di output, farebbe rifare ore di calcolo che darebbero
gli stessi risultati.

**Il digest** identifica la combinazione di parametri impiegata ed è calcolato
sulla configurazione dichiarata più i derivati.

**L'origine dei valori resta agli atti.** La configurazione registrata elenca
in ``dichiarati`` i parametri impostati nel file di ingresso
(:func:`parametri_dichiarati`): gli altri valgono il predefinito. Senza, a
esecuzione conclusa un valore predefinito sarebbe indistinguibile da una scelta
deliberata. L'elenco non entra nel digest, che identifica i valori e non il
modo in cui sono stati dati.
"""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Final

import re

import yaml

from amplicon16s.config.schema import PARAMETRI_DERIVATI, Config
from amplicon16s.io_layer.artifacts import AlberoOutput, Fase

__all__ = [
    "ConfigRisolta",
    "Derivati",
    "NOME_FILE_RISOLTO",
    "Registrazione",
    "PARAMETRI_DERIVATI",
    "PARAMETRI_SENZA_EFFETTO",
    "parametri_dichiarati",
    "registra_risolta",
    "risolvi",
    "scrivi_risolta",
    "versioni_registrate",
]

#: Nome del file con la configurazione effettivamente usata. La cartella in cui
#: finisce è ``Fase.CONFIG``: la creazione delle cartelle di output passa dal
#: servizio degli artefatti, non da qui.
NOME_FILE_RISOLTO: Final = "resolved.yaml"

#: Parametri che non incidono sui risultati, esclusi dall'impronta con cui si
#: decide se una fase conclusa è ancora valida. E' l'unico punto in cui si
#: dichiarano, e l'elenco e' volutamente prudente: nel dubbio un parametro
#: resta nell'impronta, perche' escluderne uno che incide farebbe consegnare
#: risultati calcolati con un valore diverso da quello dichiarato, senza alcun
#: errore. Per questo ``run.lockfile``, che fissa le versioni dei pacchetti
#: di calcolo, resta dentro.
#:
#: * ``run.threads``: quanti processori usare, non che cosa calcolare;
#: * ``io.out_root``: dove scrivere; spostare la cartella di un'esecuzione
#:   conclusa non deve renderla incompleta;
#: * ``retry.enabled``, ``retry.max_attempts``: se e quante volte ritentare un
#:   errore ammesso al retry, la cui azione correttiva per definizione non
#:   cambia alcuna assunzione metodologica;
#: * ``run.keep_filtered_fastq``: se conservare le letture filtrate dopo che
#:   tutte le fasi le hanno usate, non come sono state calcolate;
#: * ``run.batch_size``: quanti campioni elaborare insieme, cioe' quanta
#:   memoria chiedere. S2 filtra ogni file da solo; S4 da' gli stessi byte con
#:   lotti diversi, e lo stesso risultato di ``dada(pool = "pseudo")`` in una
#:   sola chiamata (verificati nel container). E' anche la condizione perche'
#:   il retry di E-S4-02, che lo dimezza, sia legittimo.
PARAMETRI_SENZA_EFFETTO: Final[tuple[str, ...]] = (
    "run.threads",
    "io.out_root",
    "retry.enabled",
    "retry.max_attempts",
    "run.keep_filtered_fastq",
    "run.batch_size",
)


@dataclass(frozen=True)
class Derivati:
    """Parametri calcolati dalla configurazione, non letti da essa."""

    asv_len_min: int
    asv_len_max: int

    def come_chiavi(self) -> dict[str, Any]:
        """I derivati nella forma ``gruppo.parametro``, per confronti e report."""
        return {
            "asv.len_min": self.asv_len_min,
            "asv.len_max": self.asv_len_max,
        }


def _deriva_statici(config: Config) -> Derivati:
    """Calcola i derivati che discendono dalla sola configurazione."""
    # Le ASV attese hanno la lunghezza delle letture troncate, meno quanto e'
    # stato tagliato in testa. La tolleranza allarga l'intervallo nei due versi.
    base = config.filter.truncLen - config.filter.trimLeft
    len_min = base - config.asv.len_tol
    len_max = base + config.asv.len_tol

    return Derivati(asv_len_min=len_min, asv_len_max=len_max)


def _digest(
    config: Config, derivati: Derivati, escludi: tuple[str, ...] = ()
) -> str:
    """Digest della combinazione di parametri impiegata.

    Calcolato su una forma canonica: chiavi ordinate, separatori fissi,
    nessuno spazio. Due configurazioni identiche danno lo stesso digest anche
    se scritte con ordine o formattazione diversi. ``escludi`` toglie i
    parametri indicati nella forma ``gruppo.parametro``.
    """
    parametri = config.model_dump(mode="json")
    for chiave in escludi:
        gruppo, nome = chiave.split(".")
        del parametri[gruppo][nome]

    canonico = json.dumps(
        {"config": parametri, "derivati": derivati.come_chiavi()},
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
    def impronta_risultati(self) -> str:
        """Impronta dei soli parametri che incidono sui risultati.

        E' il digest senza :data:`PARAMETRI_SENZA_EFFETTO`: due configurazioni
        che differiscono solo in quelli danno gli stessi risultati.
        """
        return _digest(self.config, self.derivati, PARAMETRI_SENZA_EFFETTO)

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


def parametri_dichiarati(config: Config) -> tuple[str, ...]:
    """I parametri impostati nel file di configurazione, come ``gruppo.parametro``.

    Lo dice lo schema, che tiene traccia dei campi ricevuti in ingresso
    (``model_fields_set``): confrontare i valori con i predefiniti non
    distinguerebbe un predefinito ereditato dallo stesso valore scritto di
    proposito. Un gruppo assente dal file non ha parametri dichiarati.
    """
    dichiarati = []
    for gruppo in Config.model_fields:
        if gruppo not in config.model_fields_set:
            continue
        modello = getattr(config, gruppo)
        dichiarati += [
            f"{gruppo}.{chiave}"
            for chiave in type(modello).model_fields
            if chiave in modello.model_fields_set
        ]
    return tuple(dichiarati)


def risolvi(config: Config) -> ConfigRisolta:
    """Risolve i parametri derivati della configurazione."""
    return ConfigRisolta(config=config, derivati=_deriva_statici(config))


_INTESTAZIONE = """\
# =========================================================================== #
# Configurazione effettivamente usata: generata, non da modificare            #
# =========================================================================== #
#
# Registra i parametri con cui l'esecuzione e' stata avviata, compresi quelli
# calcolati dalla pipeline e non dichiarati nel file di ingresso. L'elenco
# "dichiarati" dice quali parametri il file di ingresso ha impostato: gli altri
# valgono il predefinito della pipeline.
#
# Il digest identifica la combinazione di parametri: due esecuzioni con lo
# stesso digest hanno usato la stessa configurazione. E' calcolato sulla
# configurazione dichiarata piu' i derivati, in forma canonica (chiavi
# ordinate, separatori fissi).
#
# Questo file non e' riutilizzabile come configurazione di ingresso: i
# parametri derivati non sono ammessi in ingresso proprio perche' calcolati.
"""


def _testo_risolta(risolta: ConfigRisolta, aggiunte: dict[str, Any] | None = None) -> str:
    """Il testo YAML della configurazione risolta, preceduto dall'intestazione che ne
    vieta il riuso.
    """
    documento: dict[str, Any] = {
        "digest": risolta.digest,
        **(aggiunte or {}),
        "derivati": list(risolta.derivati.come_chiavi()),
        # Quali valori vengono dal file di ingresso: gli altri sono predefiniti.
        "dichiarati": list(parametri_dichiarati(risolta.config)),
        "parametri": risolta.come_mappa(),
    }
    corpo = yaml.safe_dump(
        documento,
        sort_keys=False,
        allow_unicode=True,
        default_flow_style=False,
    )
    return _INTESTAZIONE + corpo


def scrivi_risolta(risolta: ConfigRisolta, out_root: Path | str) -> Path:
    """Scrive la configurazione risolta in ``00_config``.

    Passa dal servizio degli artefatti, quindi il file entra nel manifesto
    della fase con il proprio checksum come qualunque altro artefatto: la
    configurazione usata è un artefatto dell'esecuzione, non un file a parte.
    Scrive sempre ``resolved.yaml``: chi deve conservare le versioni
    precedenti usa :func:`registra_risolta`.
    """
    albero = AlberoOutput(out_root)
    artefatto = albero.scrivi_testo(Fase.CONFIG, NOME_FILE_RISOLTO, _testo_risolta(risolta))
    return artefatto.percorso


#: Le versioni successive alla prima: resolved_2.yaml, resolved_3.yaml, ...
_VERSIONE: Final = re.compile(r"^resolved_(\d+)\.yaml$")


def versioni_registrate(out_root: Path | str) -> list[Path]:
    """Le configurazioni registrate in ``00_config``, dalla prima all'ultima."""
    cartella = AlberoOutput(out_root).cartella(Fase.CONFIG)
    versioni: list[tuple[int, Path]] = []
    if (cartella / NOME_FILE_RISOLTO).is_file():
        versioni.append((1, cartella / NOME_FILE_RISOLTO))
    if cartella.is_dir():
        for percorso in cartella.iterdir():
            trovato = _VERSIONE.match(percorso.name)
            if trovato:
                versioni.append((int(trovato.group(1)), percorso))
    return [p for _, p in sorted(versioni)]


@dataclass(frozen=True)
class Registrazione:
    """Con quale configurazione si esegue, e se è stata registrata ora."""

    percorso: Path
    digest: str
    #: Vero se il file è stato scritto adesso.
    nuova: bool
    #: Per una nuova versione, i parametri che differiscono dalla precedente.
    differenze: tuple[str, ...] = ()


def _appiattisci(mappa: dict[str, Any]) -> dict[str, Any]:
    """Da ``{gruppo: {chiave: valore}}`` a ``{"gruppo.chiave": valore}``."""
    return {
        f"{gruppo}.{chiave}": valore
        for gruppo, parametri in mappa.items()
        for chiave, valore in parametri.items()
    }


def registra_risolta(risolta: ConfigRisolta, out_root: Path | str) -> Registrazione:
    """Registra la configurazione con cui un'esecuzione parte, senza sovrascrivere.

    * nessuna configurazione registrata: scrive ``resolved.yaml``;
    * l'ultima registrata ha lo stesso digest e gli stessi parametri
      dichiarati: non scrive nulla;
    * altrimenti la conserva e scrive la nuova accanto, ``resolved_2.yaml``,
      ``resolved_3.yaml`` e cosi' via, con cio' che differisce dalla
      precedente: i valori, e i parametri dichiarati o non piu' dichiarati.

    Il confronto è con l'ultima e non con una qualunque: tornare a una
    configurazione gia' usata e' a sua volta un cambiamento, e ricostruire
    con quale configurazione ha girato ciascuna ripresa richiede di vederlo.
    Una versione si distingue dalla precedente per i valori, cioe' per il
    digest, e per l'elenco dei parametri dichiarati: dichiarare un parametro
    al valore che aveva per difetto non cambia alcun risultato, ma lo rende
    una scelta e non piu' un valore ereditato, e il report deve poterlo leggere.
    """
    versioni = versioni_registrate(out_root)
    if not versioni:
        return Registrazione(scrivi_risolta(risolta, out_root), risolta.digest, True)

    ultima = versioni[-1]
    registrata = yaml.safe_load(ultima.read_text(encoding="utf-8"))
    dichiarati_prima = registrata.get("dichiarati")
    dichiarati_ora = list(parametri_dichiarati(risolta.config))
    stessi_dichiarati = dichiarati_prima is None or set(dichiarati_prima) == set(dichiarati_ora)
    if registrata.get("digest") == risolta.digest and stessi_dichiarati:
        return Registrazione(ultima, risolta.digest, False)

    prima = _appiattisci(registrata.get("parametri", {}))
    ora = _appiattisci(risolta.come_mappa())
    differenze = tuple(
        f"{chiave}: {prima.get(chiave)!r} -> {ora.get(chiave)!r}"
        for chiave in sorted(set(prima) | set(ora))
        if prima.get(chiave) != ora.get(chiave)
    )
    if not stessi_dichiarati:
        # Solo per i parametri il cui valore non e' cambiato: per gli altri la
        # differenza di valore dice gia' tutto.
        invariati = {chiave for chiave in ora if prima.get(chiave) == ora.get(chiave)}
        differenze += tuple(
            f"{chiave}: ora dichiarato"
            for chiave in sorted((set(dichiarati_ora) - set(dichiarati_prima)) & invariati)
        ) + tuple(
            f"{chiave}: non piu' dichiarato"
            for chiave in sorted((set(dichiarati_prima) - set(dichiarati_ora)) & invariati)
        )
    nome = f"resolved_{len(versioni) + 1}.yaml"
    testo = _testo_risolta(
        risolta,
        {"precedente": ultima.name, "differenze_dalla_precedente": list(differenze)},
    )
    artefatto = AlberoOutput(out_root).scrivi_testo(Fase.CONFIG, nome, testo)
    return Registrazione(artefatto.percorso, risolta.digest, True, differenze)
