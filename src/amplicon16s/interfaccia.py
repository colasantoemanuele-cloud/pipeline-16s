"""Uscite in JSON per i programmi che guidano la pipeline dall'esterno.

Un programma esterno (un'interfaccia grafica, uno scheduler) non deve leggere
il testo pensato per una persona ne' conoscere il codice della pipeline: gli
bastano tre documenti JSON, ciascuno con il campo ``versione_formato``.

* :func:`schema`: che cosa si puo' scrivere in un file di configurazione. Per
  ogni parametro il tipo, se e' obbligatorio, i valori ammessi, il predefinito
  con la sua specie e la descrizione; accanto, lo schema JSON generato dal
  modello.
* :func:`esito_validazione`: l'esito di ``validate``. Per ogni problema il
  parametro, il codice del catalogo, il messaggio e l'azione suggerita.
* :func:`stato`: a che punto e' un'esecuzione, letto dalla sua cartella di
  uscita senza scrivervi nulla: fasi concluse (dai manifesti), fase in corso
  (dal file di stato che l'esecutore aggiorna), arresto con il codice e il
  punto di ripresa.

``versione_formato`` cambia solo se un campo esistente cambia significato o
sparisce: un campo aggiunto non la cambia, e chi legge deve ignorare i campi
che non conosce.

**Lo stato non ricalcola i checksum.** Una fase risulta conclusa se il suo
manifesto esiste e si legge: e' cio' che l'esecuzione ha dichiarato, e basta a
seguire un'esecuzione in corso senza contendere il disco con essa. Se quegli
artefatti siano ancora integri e validi per una data configurazione lo
stabilisce ``resume``, che li rilegge.
"""

from __future__ import annotations

import json
import os
import re
import socket
import threading
import time
from collections.abc import Iterable, Mapping
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Final

import yaml

from amplicon16s import __version__
from amplicon16s.config import defaults
from amplicon16s.config.schema import Config
from amplicon16s.errors.catalog import CATALOGO
from amplicon16s.io_layer.artifacts import Fase, nome_manifesto_passo
from amplicon16s.runner.graph import GRAFO

__all__ = [
    "INTERVALLO_BATTITO_S",
    "NOME_STATO",
    "StatoInCorso",
    "descrivi_parametri",
    "esecuzione_attiva",
    "VERSIONE_FORMATO",
    "chiavi_dello_schema",
    "esito_validazione",
    "parametri_nel_testo",
    "problema",
    "schema",
    "stato",
]

#: La versione del formato dei tre documenti.
VERSIONE_FORMATO: Final = 1

#: Il file di stato dell'esecuzione in corso, in ``99_logs``: lo scrive
#: l'esecutore all'avvio, lo aggiorna a ogni fase e con un battito periodico, e
#: lo toglie quando l'esecuzione finisce.
NOME_STATO: Final = "stato_esecuzione.json"

#: Ogni quanti secondi l'esecutore aggiorna il battito. Un'esecuzione il cui
#: battito e' fermo da oltre tre intervalli e' riportata come interrotta: fra
#: container diversi l'identificativo del processo non si puo' verificare,
#: l'eta' di un file si'.
INTERVALLO_BATTITO_S: Final = 20.0

#: I parametri che lo schema ammette in una sola forma pur non essendo un
#: insieme chiuso di un solo valore: due valori logici con un solo valore
#: realizzato.
_UNICA_FORMA_LOGICA: Final = (
    "tax.assign_species", "out.taxa_are_rows", "pair.trim_overhang")


class StatoInCorso:
    """Il file di stato dell'esecuzione in corso, per chi la osserva da fuori.

    Sta nella cartella indicata (``99_logs`` per la pipeline) e dice quale fase
    e' in esecuzione e da quando; un
    battito periodico ne aggiorna l'istante, cosi' chi lo legge distingue
    un'esecuzione viva da una uccisa senza dichiarazione anche da un altro
    container, dove l'identificativo del processo non e' verificabile. Non e'
    un artefatto di fase: non entra in alcun manifesto, e sparisce quando
    l'esecuzione finisce, comunque finisca.

    Razionale sistemistico: il battito gira in un thread che scrive un file e
    lo rinomina, senza usare il log ne' altre risorse condivise con le fasi:
    una fase che avvia processi figli non puo' ereditarne un blocco.
    """

    def __init__(self, cartella: Path, comando: str, intervallo_s: float = INTERVALLO_BATTITO_S) -> None:
        self._percorso = cartella / NOME_STATO
        self._intervallo = intervallo_s
        self._blocco = threading.Lock()
        self._fine = threading.Event()
        self._thread: threading.Thread | None = None
        self._documento: dict[str, Any] = {
            "pid": os.getpid(),
            "host": socket.gethostname(),
            "comando": comando,
            "avviata": self._adesso(),
            "passo": None,
            "fase_avviata": None,
            "intervallo_battito_s": intervallo_s,
        }

    @staticmethod
    def _adesso() -> str:
        return datetime.now(timezone.utc).isoformat(timespec="seconds")

    def _scrivi(self) -> None:
        with self._blocco:
            self._documento["aggiornato"] = self._adesso()
            self._documento["aggiornato_epoch"] = round(time.time(), 3)
            temporaneo = self._percorso.with_name(f".{NOME_STATO}.scrittura")
            temporaneo.write_text(
                json.dumps(self._documento, ensure_ascii=False) + "\n", encoding="utf-8")
            os.replace(temporaneo, self._percorso)

    def _batti(self) -> None:
        while not self._fine.wait(self._intervallo):
            try:
                self._scrivi()
            except OSError:
                # Lo stato e' un'informazione per chi osserva: non poterlo
                # scrivere non deve fermare il calcolo.
                pass

    def avvia(self) -> None:
        """Scrive lo stato e avvia il battito."""
        self._percorso.parent.mkdir(parents=True, exist_ok=True)
        self._scrivi()
        self._thread = threading.Thread(target=self._batti, name="battito", daemon=True)
        self._thread.start()

    def fase(self, passo: Any) -> None:
        """Registra la fase che sta per essere eseguita."""
        self._documento["passo"] = str(passo)
        self._documento["fase_avviata"] = self._adesso()
        try:
            self._scrivi()
        except OSError:
            pass

    def da_rifare(self, passi: Iterable[Any]) -> None:
        """Registra le fasi che l'esecuzione deve ancora calcolare: chi osserva
        distingue cosi' una fase conclusa da una che ha solo gli artefatti di
        un'esecuzione precedente.
        """
        self._documento["da_rifare"] = [str(p) for p in passi]
        try:
            self._scrivi()
        except OSError:
            pass

    def chiudi(self) -> None:
        """Ferma il battito e toglie il file."""
        self._fine.set()
        if self._thread is not None:
            self._thread.join(timeout=5)
        with self._blocco:
            self._percorso.unlink(missing_ok=True)


# --------------------------------------------------------------------------- #
# Schema                                                                       #
# --------------------------------------------------------------------------- #


def chiavi_dello_schema() -> tuple[str, ...]:
    """Tutti i parametri dello schema, nella forma ``gruppo.parametro``."""
    return tuple(
        f"{gruppo}.{nome}"
        for gruppo, campo in Config.model_fields.items()
        for nome in campo.annotation.model_fields
    )


def _senza_nullo(proprieta: Mapping[str, Any]) -> tuple[dict[str, Any], bool]:
    """La descrizione del tipo senza l'alternativa nulla, e se il nullo e' ammesso."""
    alternative = proprieta.get("anyOf")
    if not alternative:
        return dict(proprieta), False
    non_nulle = [a for a in alternative if a.get("type") != "null"]
    ammette_nullo = len(non_nulle) < len(alternative)
    if len(non_nulle) == 1:
        return {**proprieta, **non_nulle[0]}, ammette_nullo
    return {**proprieta, "alternative": non_nulle}, ammette_nullo


_TIPI: Final = {
    "string": "stringa", "integer": "intero", "number": "reale", "boolean": "logico",
    "object": "mappa",
}


def _tipo(tipo: Mapping[str, Any]) -> str:
    """Il tipo di un parametro in una parola, per chi costruisce un modulo."""
    if "alternative" in tipo:
        return " oppure ".join(_tipo(a) for a in tipo["alternative"])
    if "enum" in tipo or "const" in tipo:
        return "scelta"
    if tipo.get("type") == "array":
        return "elenco di " + _tipo(tipo.get("items", {})).replace("stringa", "stringhe")
    if tipo.get("format") == "path":
        return "percorso"
    return _TIPI.get(str(tipo.get("type")), "valore")


def _valori_ammessi(tipo: Mapping[str, Any]) -> list[Any] | None:
    """I valori di un insieme chiuso, o ``None`` se il parametro non ne ha uno."""
    if "alternative" in tipo:
        valori: list[Any] = []
        for alternativa in tipo["alternative"]:
            if alternativa.get("type") == "boolean":
                valori += [True, False]
            else:
                valori += _valori_ammessi(alternativa) or []
        return valori or None
    if "enum" in tipo:
        return list(tipo["enum"])
    if "const" in tipo:
        return [tipo["const"]]
    if tipo.get("type") == "array" and isinstance(tipo.get("items"), dict):
        # Un elenco i cui elementi vengono da un insieme chiuso.
        return _valori_ammessi(tipo["items"])
    return None


def _specie(chiave: str, obbligatorio: bool, predefinito: Any, ammessi: list[Any] | None) -> str:
    """Di che specie e' un parametro, con le stesse parole del modello."""
    if chiave in defaults.OBBLIGATORI:
        return "obbligatorio"
    if chiave in defaults.CONDIZIONALI or chiave in defaults.CONDIZIONALI_RESPINTI:
        return "condizionale"
    if obbligatorio:
        return "percorso"
    if chiave in defaults.SCELTE_DI_MODALITA:
        return "modalita"
    if chiave in defaults.STANDARD_DEL_METODO:
        return "standard"
    if (ammessi is not None and len(ammessi) == 1) or chiave in _UNICA_FORMA_LOGICA:
        return "unica_forma"
    if predefinito is None:
        return "facoltativo"
    return "strumento"


def descrivi_parametri(radice: Any) -> list[dict[str, Any]]:
    """I parametri di un modello di configurazione a gruppi, uno per voce:
    nome, tipo, obbligatorieta', valori ammessi, predefinito, vincoli e
    descrizione. La specie la aggiunge chi conosce le regole dello strumento.
    """
    definizioni = radice.model_json_schema().get("$defs", {})
    parametri = []
    for gruppo, campo_gruppo in radice.model_fields.items():
        modello = campo_gruppo.annotation
        proprieta = definizioni[modello.__name__]["properties"]
        for nome, campo in modello.model_fields.items():
            tipo, ammette_nullo = _senza_nullo(proprieta[nome])
            obbligatorio = campo.is_required()
            voce: dict[str, Any] = {
                "nome": f"{gruppo}.{nome}",
                "gruppo": gruppo,
                "chiave": nome,
                "tipo": _tipo(tipo),
                "obbligatorio": obbligatorio,
                "ammette_nullo": ammette_nullo,
                "valori_ammessi": _valori_ammessi(tipo),
                "specie": None,
                "descrizione": campo.description,
            }
            if not obbligatorio:
                predefinito = campo.get_default(call_default_factory=True)
                voce["predefinito"] = json.loads(json.dumps(predefinito, default=str))
            for limite in ("minimum", "maximum", "exclusiveMinimum", "exclusiveMaximum",
                           "minLength", "minItems", "pattern"):
                if limite in tipo:
                    voce.setdefault("vincoli", {})[limite] = tipo[limite]
            parametri.append(voce)
    return parametri


def _primo_capoverso(testo: str | None) -> str:
    """Il primo capoverso di una docstring su una riga sola: la frase che
    descrive un gruppo puo' andare a capo, e la sola prima riga la troncherebbe."""
    capoverso = (testo or "").strip().split("\n\n")[0]
    return " ".join(capoverso.split())


def documento_schema(strumento: str, versione: str, radice: Any,
                     parametri: list[dict[str, Any]]) -> dict[str, Any]:
    """Il documento di ``schema`` per uno strumento."""
    return {
        "versione_formato": VERSIONE_FORMATO,
        "comando": "schema",
        "strumento": strumento,
        "versione_strumento": versione,
        "gruppi": {
            gruppo: _primo_capoverso(campo.annotation.__doc__)
            for gruppo, campo in radice.model_fields.items()
        },
        "parametri": parametri,
        "json_schema": radice.model_json_schema(),
    }


#: Il valore che rende vera la condizione di una modalita'.
_VALORE_DELLA_MODALITA: Final = {
    defaults.APPAIATA: defaults.LAYOUT_APPAIATO,
    defaults.CERCATO: defaults.PRIMER_CERCATO,
    defaults.IDTAXA: defaults.CLASSIFICATORE_IDTAXA,
    defaults.BAYESIANO: defaults.TAX_CLASSIFIER,
}


def schema() -> dict[str, Any]:
    """Lo schema della configurazione della pipeline, per un programma esterno."""
    parametri = descrivi_parametri(Config)
    for voce in parametri:
        voce["specie"] = _specie(
            voce["nome"], voce["obbligatorio"], voce.get("predefinito"), voce["valori_ammessi"])
        if voce["nome"] in defaults.STANDARD_DEL_METODO:
            voce["fonte"] = defaults.STANDARD_DEL_METODO[voce["nome"]]
        condizioni = {**defaults.CONDIZIONALI, **defaults.CONDIZIONALI_RESPINTI}
        if voce["nome"] in condizioni:
            # Obbligatorio se la condizione e' vera, omissibile altrimenti (e,
            # dichiarato comunque, validato come ogni valore ma senza effetto
            # sui risultati). La condizione puo' avere piu' alternative, di cui
            # ne basta una: la prima e' richiesto_da, le altre il suo "oppure".
            # I parametri di una modalita' (appaiata, primer cercato, uno dei
            # due classificatori) con la condizione falsa non si dichiarano:
            # sono respinti.
            prima, *altre = [
                {"parametro": parametro, "condizione": forma,
                 **({"valori": [True]} if forma == defaults.VERO else {}),
                 **({"valori": [_VALORE_DELLA_MODALITA[forma]]}
                    if forma in _VALORE_DELLA_MODALITA else {})}
                for parametro, forma in condizioni[voce["nome"]]
            ]
            voce["richiesto_da"] = {**prima, **({"oppure": altre} if altre else {})}
    return documento_schema("amplicon16s", __version__, Config, parametri)


# --------------------------------------------------------------------------- #
# Validazione                                                                  #
# --------------------------------------------------------------------------- #

_NOME_DI_PARAMETRO: Final = re.compile(r"\b([a-z]+\.[A-Za-z_0-9]+)\b")


def parametri_nel_testo(testo: str, chiavi: Iterable[str] | None = None) -> list[str]:
    """I parametri dello schema nominati in un messaggio, nell'ordine in cui
    compaiono e senza ripetizioni.
    """
    note = set(chiavi_dello_schema() if chiavi is None else chiavi)
    trovati: list[str] = []
    for nome in _NOME_DI_PARAMETRO.findall(testo):
        if nome in note and nome not in trovati:
            trovati.append(nome)
    return trovati


def problema(
    codice: str,
    messaggio: str,
    *,
    parametri: Iterable[str] | None = None,
    gate: str | None = None,
    azione: str | None = None,
) -> dict[str, Any]:
    """Un problema nella forma del documento: parametro, codice, messaggio, azione.

    ``parametro`` e' il primo dei parametri coinvolti, o nullo se il problema
    non riguarda un parametro (un file che non si apre, una risorsa che manca);
    ``parametri`` li elenca tutti. L'azione e' quella del catalogo dei codici.
    """
    coinvolti = list(parametri) if parametri is not None else parametri_nel_testo(messaggio)
    voce = CATALOGO.get(codice)
    return {
        "parametro": coinvolti[0] if coinvolti else None,
        "parametri": coinvolti,
        "codice": codice,
        "categoria": voce.categoria.value if voce else None,
        "messaggio": messaggio,
        "azione": azione if azione is not None else (voce.azione if voce else None),
        "gate": gate,
    }


def esito_validazione(
    *,
    configurazione: Path | str,
    codice_uscita: int,
    problemi: list[dict[str, Any]],
    avvisi: list[dict[str, Any]] | None = None,
    cartella: Path | str | None = None,
    controlli_spenti: list[dict[str, Any]] | None = None,
) -> dict[str, Any]:
    """Il documento di ``validate --json``.

    ``controlli_spenti`` elenca le soglie il cui controllo non puo' scattare
    (:mod:`amplicon16s.soglie.tabella`): la configurazione e' valida, ma chi la
    esegue deve saperlo. ``None`` se non sono stati valutati, perche' la
    configurazione e' stata respinta prima.
    """
    return {
        "versione_formato": VERSIONE_FORMATO,
        "comando": "validate",
        "strumento": "amplicon16s",
        "configurazione": str(configurazione),
        "cartella": str(cartella) if cartella is not None else None,
        "valida": codice_uscita == 0 and not problemi,
        "codice_uscita": codice_uscita,
        "problemi": problemi,
        "avvisi": avvisi or [],
        "controlli_spenti": controlli_spenti,
    }


# --------------------------------------------------------------------------- #
# Stato di un'esecuzione                                                       #
# --------------------------------------------------------------------------- #


def _json(percorso: Path) -> Any | None:
    """Il contenuto di un file JSON, o ``None`` se non c'e' o non si legge: un
    file che l'esecuzione sta sostituendo in questo momento non e' un errore.
    """
    try:
        return json.loads(percorso.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return None


def _configurazione_registrata(radice: Path) -> dict[str, Any] | None:
    """L'ultima configurazione registrata in ``00_config``, se c'e'."""
    cartella = radice / Fase.CONFIG.value
    if not cartella.is_dir():
        return None

    def ordine(percorso: Path) -> int:
        trovato = re.fullmatch(r"resolved_(\d+)\.yaml", percorso.name)
        return int(trovato.group(1)) if trovato else 1

    versioni = sorted(cartella.glob("resolved*.yaml"), key=ordine)
    if not versioni:
        return None
    try:
        registrata = yaml.safe_load(versioni[-1].read_text(encoding="utf-8"))
    except (OSError, yaml.YAMLError):
        return None
    if not isinstance(registrata, dict):
        return None
    return {"file": versioni[-1].name, "digest": registrata.get("digest"),
            "parametri": registrata.get("parametri") or {}}


def esecuzione_attiva(cartella: Path, adesso: float | None = None) -> dict[str, Any] | None:
    """Che cosa dice il file di stato della cartella, o ``None`` se non c'e'.

    ``battito_recente`` e' vero se l'ultimo aggiornamento non e' piu' vecchio di
    tre intervalli: l'esecuzione e' viva.
    """
    esecuzione = _json(cartella / NOME_STATO)
    if not isinstance(esecuzione, dict):
        return None
    adesso = time.time() if adesso is None else adesso
    eta = adesso - float(esecuzione.get("aggiornato_epoch", 0))
    intervallo = float(esecuzione.get("intervallo_battito_s", INTERVALLO_BATTITO_S))
    return {
        "passo": esecuzione.get("passo"),
        "comando": esecuzione.get("comando"),
        "avviata": esecuzione.get("avviata"),
        "fase_avviata": esecuzione.get("fase_avviata"),
        "ultimo_battito": esecuzione.get("aggiornato"),
        "secondi_dall_ultimo_battito": round(max(eta, 0.0), 1),
        "battito_recente": eta <= 3 * intervallo,
        "da_rifare": [str(p) for p in esecuzione.get("da_rifare") or []],
        "pid": esecuzione.get("pid"),
        "host": esecuzione.get("host"),
    }


def stato(cartella: Path | str, adesso: float | None = None) -> dict[str, Any]:
    """Lo stato dell'esecuzione nella cartella di uscita indicata.

    Non scrive nulla e non apre il log: legge i manifesti di fase, il file di
    stato dell'esecutore e il punto di ripresa. ``stato`` vale:

    * ``non_avviata``: la cartella non c'e' o non contiene un'esecuzione;
    * ``in_corso``: un'esecuzione e' attiva (battito recente);
    * ``interrotta``: il file di stato c'e' ma il battito e' fermo: il processo
      e' stato ucciso senza poter dichiarare nulla;
    * ``ferma``: l'esecuzione si e' fermata con un codice e un punto di ripresa;
    * ``conclusa``: ogni fase attiva ha il suo manifesto;
    * ``incompleta``: alcune fasi sono concluse e nessuna esecuzione e' attiva
      (per esempio dopo ``validate``, che esegue la sola prima fase).
    """
    radice = Path(cartella)
    adesso = time.time() if adesso is None else adesso
    registrata = _configurazione_registrata(radice)
    parametri = (registrata or {}).get("parametri", {})
    # Senza una configurazione registrata non si sa se la fase facoltativa e'
    # attiva: la si riporta come sconosciuta, non come mancante.
    filogenesi = (parametri.get("phylo") or {}).get("enabled")

    attiva = esecuzione_attiva(radice / Fase.LOGS.value, adesso)
    # Le fasi che l'esecuzione in corso (o interrotta) deve ancora calcolare.
    da_rifare = set(attiva["da_rifare"]) if attiva is not None else set()

    # La fase che esiste solo in modalita' appaiata si elenca se la
    # configurazione registrata dichiara quella modalita'.
    appaiata = (parametri.get("io") or {}).get("layout") == "paired-end"

    fasi = []
    for passo in GRAFO.ordine():
        nodo = GRAFO.nodo(passo)
        voce: dict[str, Any] = {
            "passo": str(passo), "descrizione": nodo.descrizione,
            "cartella": nodo.cartella.value,
        }
        if not nodo.elencata_se_inattiva and not appaiata:
            continue
        manifesto = _json(radice / nodo.cartella.value / nome_manifesto_passo(str(passo)))
        # La sola fase che una configurazione puo' disattivare e' la filogenesi:
        # quella della modalita' appaiata, se elencata, e' attiva.
        facoltativa = nodo.facoltativa and nodo.elencata_se_inattiva
        if facoltativa and filogenesi is False:
            voce["stato"] = "disattivata"
        elif attiva is not None and attiva["battito_recente"] and attiva["passo"] == str(passo):
            voce["stato"] = "in_corso"
        elif isinstance(manifesto, dict):
            voce["stato"] = "conclusa"
            voce["artefatti"] = len(manifesto.get("artefatti") or [])
            voce["conclusa_il"] = manifesto.get("conclusa")
            voce["versione"] = (manifesto.get("calcolata_su") or {}).get("versione")
            voce["degradazioni"] = [
                d.get("codice") for d in manifesto.get("degradazioni") or []
            ]
            # Il manifesto e' quello di un'esecuzione precedente: lo stato resta
            # quello che il manifesto dice, e il campo lo distingue.
            if str(passo) in da_rifare:
                voce["da_rifare"] = True
        elif facoltativa and filogenesi is None:
            voce["stato"] = "sconosciuta"
        else:
            voce["stato"] = "da_eseguire"
        fasi.append(voce)

    punto = _json(radice / Fase.LOGS.value / "punto_di_ripresa.json")
    concluse = [f["passo"] for f in fasi if f["stato"] == "conclusa"]
    mancanti = [f["passo"] for f in fasi if f["stato"] in ("da_eseguire", "in_corso")]
    if attiva is not None and attiva["battito_recente"]:
        complessivo = "in_corso"
    elif attiva is not None:
        complessivo = "interrotta"
    elif isinstance(punto, dict):
        complessivo = "ferma"
    elif concluse and not mancanti:
        complessivo = "conclusa"
    elif concluse or registrata is not None:
        complessivo = "incompleta"
    else:
        complessivo = "non_avviata"

    return {
        "versione_formato": VERSIONE_FORMATO,
        "comando": "status",
        "strumento": "amplicon16s",
        "cartella": str(radice),
        "stato": complessivo,
        "configurazione": (
            {"file": registrata["file"], "digest": registrata["digest"]} if registrata else None
        ),
        "fasi": fasi,
        "fasi_concluse": concluse,
        "fasi_da_rifare": [f["passo"] for f in fasi if f.get("da_rifare")],
        "fase_in_corso": attiva["passo"] if complessivo == "in_corso" and attiva else None,
        "esecuzione": attiva,
        # Il punto di ripresa di un arresto resta sul disco finche' un'esecuzione
        # non arriva in fondo: durante una ripresa in corso descrive l'arresto
        # precedente, e lo stato complessivo lo dice.
        "arresto": punto if isinstance(punto, dict) else None,
    }
