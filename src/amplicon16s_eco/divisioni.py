"""Le divisioni dei metadati: il catalogo dei modi noti di dividere i campioni.

Una **divisione** e' un modo di dividere i campioni in gruppi a partire dalla
tabella dei campioni dell'oggetto: una colonna cosi' com'e', il prefisso dei
suoi valori, una mappatura nota di valori in gruppi, l'ordine di tempo, fasce
di un valore numerico, oppure l'incrocio di due divisioni per una griglia di
pannelli. Il catalogo (``divisioni.yaml``, versionato con il pacchetto) elenca
le divisioni note con il nome, la spiegazione, le colonne a cui si applicano e
l'origine dei nomi.

Chi lo usa:

* ``eco init`` (:mod:`amplicon16s_eco.init`) applica il catalogo ai metadati
  dell'oggetto e propone solo le divisioni che funzionano su quel dataset
  (:func:`proponi`), piu' le divisioni generiche: ogni colonna con gruppi
  ripetuti e il prefisso dei valori di una colonna (``<colonna>_prefisso``);
* ``eco validate`` e ``eco run`` (:mod:`amplicon16s_eco.esecuzione`) ricavano
  la regola dal nome scritto nella configurazione e la espandono in cio' che
  le analisi conoscono gia': una colonna, un gruppo derivato di
  ``design.groupings``, una griglia ``{rows, cols}`` (:func:`espandi`).

**Come si legge un nome** nei campi che accettano una colonna
(:data:`CAMPI`), nell'ordine, e vale la prima risposta:

1. un gruppo derivato dichiarato in ``design.groupings``;
2. una colonna dell'oggetto con almeno un valore;
3. una divisione del catalogo, o una divisione generica ``<colonna>_prefisso``;
4. una colonna dell'oggetto senza valori (come nella versione 0.2.0: le
   analisi la respingono piu' avanti);
5. altrimenti il nome e' sconosciuto (``E-ECO-47``).

Una colonna con valori vince quindi sempre su una divisione con lo stesso nome,
e una configurazione della versione 0.2.0 si legge come prima. La colonna
vuota dell'inventario della pipeline (``modulo`` quando la pipeline non ha
derivato il modulo) non e' una variabile: il nome passa alla divisione.

**Quando una divisione funziona** su un dataset (:func:`giudica`): la regola
da' almeno :data:`MINIMO_GRUPPI` gruppi e non piu' di :data:`MASSIMO_GRUPPI`,
i campioni con un gruppo sono almeno :data:`QUOTA_COPERTA` di quelli
dell'oggetto, e almeno :data:`MINIMO_GRUPPI` gruppi hanno abbastanza campioni
per i test (``stat.min_group_size`` in validazione,
:data:`DIMENSIONE_DELLE_PROPOSTE` per le proposte, quando la dimensione non e'
ancora scelta). La colonna su cui si applica la divisione si sceglie senza la
dimensione dei gruppi: la prima, nell'ordine delle fonti del catalogo e poi
delle colonne dell'oggetto, su cui la regola da' un numero di gruppi ammesso
coprendo abbastanza campioni. Cosi' ``init`` e ``validate`` scelgono sempre la
stessa colonna. I gruppi si contano su tutti i campioni dell'oggetto, prima di
``design.subset``.

Tutto e' deterministico: stesse proposte, nello stesso ordine, a ogni
esecuzione; i gruppi in ordine di nome, come li ordinano le analisi.
"""

from __future__ import annotations

import hashlib
import re
from dataclasses import dataclass
from functools import cache
from importlib import resources
from typing import Any, Final

import yaml

from amplicon16s.metadata.crosswalk import _VALORI_ASSENTI


__all__ = [
    "CAMPI", "DIMENSIONE_DELLE_PROPOSTE", "MASSIMO_GRUPPI", "MINIMO_GRUPPI", "PREFISSI_ISA",
    "PREFISSO", "QUOTA_COPERTA", "REGOLE", "SUFFISSO_DEL_PREFISSO", "Applicazione",
    "Divisione", "Espansione", "Fonte", "Proposta", "Tabella", "applica", "carica_catalogo",
    "catalogo", "chiave", "espandi", "giudica", "impronta_del_catalogo", "proponi",
    "stessa_partizione", "traduci",
]

#: I prefissi delle colonne ISA-Tab come li scrive la pipeline nell'oggetto
#: (``Factor Value[Time]`` diventa ``factor_value_time``): la chiave di una
#: colonna non li comprende.
PREFISSI_ISA: Final = ("factor_value_", "characteristics_", "comment_", "parameter_value_")

#: Le regole che una divisione puo' avere.
REGOLE: Final = ("colonna", "prefisso", "mappatura", "ordine", "fasce", "incrocio")

#: Il criterio delle divisioni generiche per prefisso: lettere iniziali seguite
#: da cifre (il modulo di una posizione come «NOD1D4»), primo gruppo di cattura.
PREFISSO: Final = r"^([A-Za-z]+[0-9]+)"
#: Il nome di una divisione generica per prefisso: la colonna e questo suffisso.
SUFFISSO_DEL_PREFISSO: Final = "_prefisso"
#: Una divisione generica per prefisso si propone solo per una colonna con piu'
#: di tanti valori distinti, e se riunisce in media almeno due valori per gruppo:
#: altrimenti la colonna stessa basta.
_VALORI_PER_IL_PREFISSO: Final = 8

MINIMO_GRUPPI: Final = 2
MASSIMO_GRUPPI: Final = 12
QUOTA_COPERTA: Final = 0.5
#: La dimensione dei gruppi con cui si giudicano le proposte di ``eco init``,
#: quando ``stat.min_group_size`` non e' ancora scelto: il valore di mezzo di
#: quelli che init elenca (2, 5, 10).
DIMENSIONE_DELLE_PROPOSTE: Final = 5

#: I campi della configurazione che accettano il nome di una divisione, come
#: accettano una colonna; ``comp.facet_grid`` accetta un incrocio, o una
#: divisione in ``rows`` e ``cols``.
CAMPI: Final = ("design.variable", "design.other_variables", "comp.group_by", "ord.color_by",
                "comp.facet_wrap", "comp.facet_grid")
#: I campi che nominano variabili gia' nominate nei precedenti: seguono la stessa
#: sostituzione.
_CAMPI_CHE_SEGUONO: Final = ("constrained.variables", "constrained.conditions")

_NOME_DEL_FILE: Final = "divisioni.yaml"


def chiave(nome: str) -> str:
    """La chiave con cui il catalogo riconosce una colonna.

    Minuscole, ogni sequenza di segni diversi da lettere e cifre come un
    trattino basso (la regola con cui la pipeline nomina le colonne
    nell'oggetto), senza uno dei prefissi :data:`PREFISSI_ISA`: «Factor
    Value[Sample Location]», «factor_value_sample_location» e
    «Sample-Location» hanno la stessa chiave ``sample_location``.
    """
    testo = re.sub(r"[^a-z0-9]+", "_", nome.casefold()).strip("_")
    for prefisso in PREFISSI_ISA:
        if testo.startswith(prefisso) and len(testo) > len(prefisso):
            return testo[len(prefisso):]
    return testo


def _nome_canonico(nome: str) -> str:
    """Un nome di divisione per il confronto: spazi ridotti, senza maiuscole."""
    return " ".join(nome.split()).casefold()


def _presente(valore: str | None) -> bool:
    """Vero se il valore non e' nullo, vuoto, ne' un modo di dire «assente»."""
    return valore is not None and valore.strip().casefold() not in _VALORI_ASSENTI


# --------------------------------------------------------------------------- #
# Il catalogo                                                                  #
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class Fonte:
    """Una colonna a cui si applica una divisione, con la regola."""

    #: La chiave della colonna (:func:`chiave`).
    colonna: str
    regola: str
    #: prefisso: l'espressione regolare, il cui primo gruppo da' il gruppo.
    espressione: str | None = None
    #: mappatura: per gruppo, nell'ordine, le espressioni che lo riconoscono.
    gruppi: tuple[tuple[str, tuple[str, ...]], ...] = ()
    #: mappatura: il numero d'ordine davanti ai nomi dei gruppi.
    ordinata: bool = False
    #: fasce: i limiti crescenti e le etichette (una in piu' dei limiti).
    limiti: tuple[float, ...] = ()
    etichette: tuple[str, ...] = ()

    def descrizione(self, colonna: str) -> str:
        """La regola in parole, sulla colonna dell'oggetto data."""
        if self.regola == "colonna":
            return f"colonna {colonna} cosi' com'e'"
        if self.regola == "prefisso":
            return f"prefisso {self.espressione} dei valori di {colonna}"
        if self.regola == "mappatura":
            return f"mappatura nota dei valori di {colonna} in {', '.join(_nomi_mappati(self))}"
        if self.regola == "ordine":
            return f"valori di {colonna} in ordine di tempo"
        return (f"fasce di {colonna} (limiti {', '.join(_numero_testo(x) for x in self.limiti)})")


@dataclass(frozen=True)
class Divisione:
    """Una voce del catalogo."""

    nome: str
    spiegazione: str
    origine: str
    #: Le colonne, nell'ordine in cui si provano; vuoto per un incrocio.
    fonti: tuple[Fonte, ...] = ()
    #: Per un incrocio: le divisioni delle righe e delle colonne della griglia.
    righe: str | None = None
    colonne: str | None = None

    @property
    def incrocio(self) -> bool:
        """Vero per una griglia di due divisioni."""
        return self.righe is not None


def _numero_testo(x: float) -> str:
    """Un limite come testo, senza decimali se e' intero."""
    return str(int(x)) if float(x).is_integer() else str(x)


def _nomi_mappati(fonte: Fonte) -> list[str]:
    """I nomi dei gruppi di una mappatura, con il numero d'ordine se e' ordinata."""
    return [f"{i} {nome}" if fonte.ordinata else nome
            for i, (nome, _) in enumerate(fonte.gruppi, start=1)]


def _etichette_delle_fasce(fonte: Fonte) -> list[str]:
    """I nomi dei gruppi di una divisione in fasce: il numero d'ordine e l'etichetta."""
    return [f"{i} {e}" for i, e in enumerate(fonte.etichette, start=1)]


def _fonte(dati: Any, dove: str) -> Fonte:
    """Una fonte del catalogo letta e verificata."""
    if not isinstance(dati, dict):
        raise ValueError(f"{dove}: una fonte e' una mappa")
    ammesse = {"colonna", "regola", "espressione", "gruppi", "ordinata", "limiti", "etichette"}
    estranee = sorted(set(dati) - ammesse)
    if estranee:
        raise ValueError(f"{dove}: chiavi sconosciute {estranee}")
    colonna, regola = dati.get("colonna"), dati.get("regola")
    if not isinstance(colonna, str) or chiave(colonna) != colonna or not colonna:
        raise ValueError(f"{dove}: colonna {colonna!r} non e' una chiave (minuscole e _)")
    if regola not in REGOLE or regola == "incrocio":
        raise ValueError(f"{dove}: regola {regola!r} non ammessa per una fonte")
    espressione = dati.get("espressione")
    gruppi: tuple[tuple[str, tuple[str, ...]], ...] = ()
    limiti: tuple[float, ...] = ()
    etichette: tuple[str, ...] = ()
    if regola == "prefisso":
        if not isinstance(espressione, str) or re.compile(espressione).groups < 1:
            raise ValueError(f"{dove}: un prefisso chiede un'espressione con un gruppo di cattura")
    elif espressione is not None:
        raise ValueError(f"{dove}: espressione vale solo per la regola prefisso")
    if regola == "mappatura":
        grezzi = dati.get("gruppi")
        if not isinstance(grezzi, dict) or len(grezzi) < 2:
            raise ValueError(f"{dove}: una mappatura chiede almeno due gruppi")
        for nome, espressioni in grezzi.items():
            if not isinstance(nome, str) or not isinstance(espressioni, list) or not espressioni:
                raise ValueError(f"{dove}: il gruppo {nome!r} non ha espressioni")
            for e in espressioni:
                re.compile(e)
        gruppi = tuple((nome, tuple(e)) for nome, e in grezzi.items())
    elif "gruppi" in dati or "ordinata" in dati:
        raise ValueError(f"{dove}: gruppi e ordinata valgono solo per la regola mappatura")
    if regola == "fasce":
        limiti = tuple(float(x) for x in dati.get("limiti") or ())
        etichette = tuple(str(e) for e in dati.get("etichette") or ())
        if not limiti or list(limiti) != sorted(set(limiti)):
            raise ValueError(f"{dove}: i limiti delle fasce sono crescenti e distinti")
        if len(etichette) != len(limiti) + 1:
            raise ValueError(f"{dove}: le etichette delle fasce sono una piu' dei limiti")
    elif "limiti" in dati or "etichette" in dati:
        raise ValueError(f"{dove}: limiti ed etichette valgono solo per la regola fasce")
    return Fonte(colonna, regola, espressione, gruppi, bool(dati.get("ordinata", False)),
                 limiti, etichette)


def carica_catalogo(testo: str) -> tuple[Divisione, ...]:
    """Le divisioni di un catalogo, dal testo YAML, verificate.

    Solleva ``ValueError`` per un catalogo malformato: nomi ripetuti o vuoti,
    un nome che finisce come una divisione generica, una spiegazione o
    un'origine mancante, una regola sconosciuta, un incrocio che non nomina
    due divisioni del catalogo diverse e non incrociate.
    """
    dati = yaml.safe_load(testo)
    if not isinstance(dati, dict) or not isinstance(dati.get("divisioni"), list):
        raise ValueError("il catalogo e' una mappa con l'elenco divisioni")
    divisioni: list[Divisione] = []
    for i, voce in enumerate(dati["divisioni"]):
        dove = f"divisione {i + 1}"
        if not isinstance(voce, dict):
            raise ValueError(f"{dove}: una divisione e' una mappa")
        nome = voce.get("nome")
        if not isinstance(nome, str) or not nome.strip() or nome != " ".join(nome.split()):
            raise ValueError(f"{dove}: nome {nome!r} vuoto o con spazi doppi")
        dove = f"divisione {nome!r}"
        if nome.endswith(SUFFISSO_DEL_PREFISSO):
            raise ValueError(f"{dove}: il suffisso {SUFFISSO_DEL_PREFISSO} e' delle divisioni "
                             "generiche")
        if any(_nome_canonico(d.nome) == _nome_canonico(nome) for d in divisioni):
            raise ValueError(f"{dove}: nome ripetuto")
        for campo in ("spiegazione", "origine"):
            if not isinstance(voce.get(campo), str) or not voce[campo].strip():
                raise ValueError(f"{dove}: {campo} mancante")
        ammesse = {"nome", "spiegazione", "origine", "fonti", "regola", "righe", "colonne"}
        estranee = sorted(set(voce) - ammesse)
        if estranee:
            raise ValueError(f"{dove}: chiavi sconosciute {estranee}")
        if voce.get("regola") == "incrocio":
            righe, colonne = voce.get("righe"), voce.get("colonne")
            semplici = {d.nome: d for d in divisioni if not d.incrocio}
            if righe not in semplici or colonne not in semplici or righe == colonne:
                raise ValueError(f"{dove}: un incrocio nomina due divisioni diverse, non "
                                 "incroci, che lo precedono nel catalogo")
            divisioni.append(Divisione(nome, voce["spiegazione"], voce["origine"],
                                       righe=righe, colonne=colonne))
            continue
        if "regola" in voce or "righe" in voce or "colonne" in voce:
            raise ValueError(f"{dove}: regola, righe e colonne valgono solo per un incrocio "
                             "(le altre regole stanno nelle fonti)")
        fonti = voce.get("fonti")
        if not isinstance(fonti, list) or not fonti:
            raise ValueError(f"{dove}: nessuna fonte")
        divisioni.append(Divisione(
            nome, voce["spiegazione"], voce["origine"],
            tuple(_fonte(f, f"{dove}, fonte {j + 1}") for j, f in enumerate(fonti))))
    return tuple(divisioni)


def _testo_del_catalogo() -> str:
    return resources.files("amplicon16s_eco").joinpath(_NOME_DEL_FILE).read_text(encoding="utf-8")


@cache
def catalogo() -> tuple[Divisione, ...]:
    """Il catalogo del pacchetto."""
    return carica_catalogo(_testo_del_catalogo())


@cache
def impronta_del_catalogo() -> str:
    """Lo SHA-256 del file del catalogo: identifica le regole usate."""
    return "sha256:" + hashlib.sha256(_testo_del_catalogo().encode("utf-8")).hexdigest()


def _del_catalogo(nome: str) -> Divisione | None:
    """La divisione del catalogo con quel nome (spazi e maiuscole non contano)."""
    return next((d for d in catalogo() if _nome_canonico(d.nome) == _nome_canonico(nome)), None)


# --------------------------------------------------------------------------- #
# I valori di tempo                                                            #
# --------------------------------------------------------------------------- #

_GIORNI_PER_UNITA: Final = {"hour": 1 / 24, "day": 1, "week": 7, "month": 30, "year": 365}
_TERMINE: Final = re.compile(
    r"\s*(minus|plus|-|\+)?\s*(\d+(?:\.\d+)?)\s*(hour|day|week|month|year)s?\b", re.IGNORECASE)
_MESI: Final = {m: i for i, m in enumerate(
    ("jan", "feb", "mar", "apr", "may", "jun", "jul", "aug", "sep", "oct", "nov", "dec"), 1)}


def durata_in_giorni(valore: str) -> float | None:
    """Una durata scritta in parole, in giorni; ``None`` se non lo e'.

    Riconosce uno o piu' termini «numero unita'» (hour, day, week, month,
    year, al singolare o al plurale), ciascuno con un segno facoltativo
    (minus, plus, -, +): «7 days», «minus10 days», «520 days plus 3 month».
    Un mese vale 30 giorni, un anno 365.
    """
    totale, posizione, trovati = 0.0, 0, 0
    testo = valore.strip()
    while posizione < len(testo):
        termine = _TERMINE.match(testo, posizione)
        if termine is None:
            return None
        segno = -1 if (termine.group(1) or "").casefold() in ("minus", "-") else 1
        totale += segno * float(termine.group(2)) * _GIORNI_PER_UNITA[termine.group(3).casefold()]
        posizione, trovati = termine.end(), trovati + 1
    return totale if trovati else None


def _data(valore: str) -> tuple[int, int, int] | None:
    """Una data «01-Nov-2017» o «2017-11-01», come (anno, mese, giorno)."""
    if (m := re.fullmatch(r"(\d{1,2})-([A-Za-z]{3})-(\d{4})", valore.strip())):
        mese = _MESI.get(m.group(2).casefold())
        return (int(m.group(3)), mese, int(m.group(1))) if mese else None
    if (m := re.fullmatch(r"(\d{4})-(\d{2})-(\d{2})", valore.strip())):
        return int(m.group(1)), int(m.group(2)), int(m.group(3))
    return None


def _numero(valore: str) -> float | None:
    """Un numero, o una durata in giorni."""
    try:
        return float(valore)
    except ValueError:
        return durata_in_giorni(valore)


# --------------------------------------------------------------------------- #
# L'applicazione a una tabella                                                 #
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class Tabella:
    """La tabella dei campioni dell'oggetto: un valore per campione e colonna."""

    campioni: tuple[str, ...]
    #: Colonna -> valori, nell'ordine dei campioni; ``None`` se mancante.
    valori: dict[str, tuple[str | None, ...]]

    @classmethod
    def dalla_descrizione(cls, descrizione: dict[str, Any]) -> "Tabella":
        """La tabella che ``descrivi.R`` scrive in ``descrizione.json``."""
        tabella = descrizione.get("tabella") or {}
        valori = tabella.get("valori") or {}
        if not isinstance(valori, dict):
            valori = {}
        return cls(tuple(tabella.get("campioni") or ()),
                   {nome: tuple(v) for nome, v in valori.items()})

    def con_valore(self, colonna: str) -> bool:
        """Vero se la colonna esiste e almeno un campione vi ha un valore."""
        return any(_presente(v) for v in self.valori.get(colonna, ()))


@dataclass(frozen=True)
class Applicazione:
    """Che cosa da' una divisione su una tabella."""

    nome: str
    #: La regola in parole, con la colonna (o le due divisioni dell'incrocio).
    regola: str
    #: La colonna dell'oggetto su cui si applica; vuota per un incrocio.
    colonna: str
    #: Il gruppo di ogni campione, nell'ordine della tabella; ``None`` se nessuno.
    gruppi_dei_campioni: tuple[str | None, ...]
    #: Ogni gruppo con i suoi campioni, in ordine di nome.
    gruppi: tuple[tuple[str, int], ...]
    #: Che cosa diventa nella configurazione: ``{"colonna": nome}``,
    #: ``{"gruppo": regola di design.groupings}`` o ``{"griglia": (righe, colonne)}``.
    espansione: dict[str, Any]
    #: Per un incrocio, le applicazioni delle righe e delle colonne.
    componenti: tuple["Applicazione", "Applicazione"] | None = None
    #: Perche' la regola non si applica a questa colonna (una regola «ordine»
    #: su valori che non sono tempi), o ``None``.
    non_applicabile: str | None = None

    @property
    def senza_gruppo(self) -> int:
        """Quanti campioni la divisione lascia senza gruppo."""
        return sum(g is None for g in self.gruppi_dei_campioni)


def _conta(gruppi: tuple[str | None, ...]) -> tuple[tuple[str, int], ...]:
    """I gruppi con i loro campioni, in ordine di nome (codepoint: l'ordine radix di R)."""
    conti: dict[str, int] = {}
    for g in gruppi:
        if g is not None:
            conti[g] = conti.get(g, 0) + 1
    return tuple(sorted(conti.items()))


def _gruppi_della_mappa(mappa: dict[str, str]) -> dict[str, list[str]]:
    """Da valore -> gruppo, i gruppi con i loro valori: la forma di design.groupings."""
    gruppi: dict[str, list[str]] = {}
    for valore, gruppo in sorted(mappa.items()):
        gruppi.setdefault(gruppo, []).append(valore)
    return dict(sorted(gruppi.items()))


def _applica_fonte(nome: str, fonte: Fonte, colonna: str,
                   valori: tuple[str | None, ...]) -> Applicazione:
    """La regola di una fonte sui valori di una colonna dell'oggetto."""
    presenti = sorted({v for v in valori if _presente(v)})
    descritta = fonte.descrizione(colonna)
    mappa: dict[str, str] = {}
    espansione: dict[str, Any]
    if fonte.regola == "prefisso":
        for v in presenti:
            trovato = re.match(fonte.espressione or PREFISSO, v)
            if trovato and trovato.group(1):
                mappa[v] = trovato.group(1)
        espansione = {"gruppo": {"from": colonna, "pattern": fonte.espressione}}
    elif fonte.regola == "mappatura":
        nomi = _nomi_mappati(fonte)
        for v in presenti:
            # Il primo gruppo, nell'ordine del catalogo, che riconosce il valore.
            trovato = next((nomi[i] for i, (_, espressioni) in enumerate(fonte.gruppi)
                            if any(re.search(e, v, re.IGNORECASE) for e in espressioni)), None)
            if trovato is not None:
                mappa[v] = trovato
        espansione = {"gruppo": {"from": colonna, "groups": _gruppi_della_mappa(mappa)}}
    elif fonte.regola == "fasce":
        etichette = _etichette_delle_fasce(fonte)
        for v in presenti:
            x = _numero(v)
            if x is not None:
                mappa[v] = etichette[sum(x >= limite for limite in fonte.limiti)]
        espansione = {"gruppo": {"from": colonna, "groups": _gruppi_della_mappa(mappa)}}
    elif fonte.regola == "ordine":
        # Solo valori che sono tutti durate, o tutti date: altrimenti l'ordine
        # di tempo non c'e', e la fonte non vale (si prova la successiva).
        tempi = [durata_in_giorni(v) for v in presenti]
        date = [_data(v) for v in presenti]
        if all(t is not None for t in tempi):
            chiavi: list[Any] = tempi
        elif all(date):
            chiavi = date
        else:
            estranei = [v for v, t, g in zip(presenti, tempi, date) if t is None and g is None]
            return Applicazione(
                nome, descritta, colonna, (None,) * len(valori), (), {},
                non_applicabile=("i valori non sono tutti durate o date (per esempio "
                                 + ", ".join(repr(v) for v in estranei[:3]) + ")"))
        ordinati = [v for _, v in sorted(zip(chiavi, presenti))]
        cifre = len(str(len(ordinati)))
        mappa = {v: f"{i:0{cifre}d} {v}" for i, v in enumerate(ordinati, start=1)}
        espansione = {"gruppo": {"from": colonna, "groups": _gruppi_della_mappa(mappa)}}
    else:
        mappa = {v: v for v in presenti}
        espansione = {"colonna": colonna}
    gruppi = tuple(mappa.get(v) if _presente(v) else None for v in valori)
    return Applicazione(nome, descritta, colonna, gruppi, _conta(gruppi), espansione)


def _struttura(app: Applicazione, campioni: int) -> str | None:
    """Perche' i gruppi di un'applicazione non sono utilizzabili (numero di
    gruppi e campioni coperti), o ``None``."""
    if app.non_applicabile:
        return app.non_applicabile
    if len(app.gruppi) < MINIMO_GRUPPI:
        return f"{len(app.gruppi)} gruppi, ne servono almeno {MINIMO_GRUPPI}"
    if len(app.gruppi) > MASSIMO_GRUPPI:
        return f"{len(app.gruppi)} gruppi, piu' di {MASSIMO_GRUPPI}"
    coperti = campioni - app.senza_gruppo
    if coperti < QUOTA_COPERTA * campioni:
        return f"un gruppo per {coperti} campioni su {campioni}, meno della meta'"
    return None


def giudica(app: Applicazione, campioni: int, dimensione: int) -> str | None:
    """Perche' una divisione applicata non funziona, o ``None`` se funziona:
    numero di gruppi, campioni coperti e gruppi grandi abbastanza per i test.
    Un incrocio funziona se funzionano le sue due divisioni."""
    if app.componenti is not None:
        for parte in app.componenti:
            motivo = giudica(parte, campioni, dimensione)
            if motivo:
                return f"{parte.nome}: {motivo}"
        return None
    motivo = _struttura(app, campioni)
    if motivo:
        return motivo
    grandi = sum(n >= dimensione for _, n in app.gruppi)
    if grandi < MINIMO_GRUPPI:
        return (f"{grandi} gruppi con almeno {dimensione} campioni, ne servono almeno "
                f"{MINIMO_GRUPPI}")
    return None


def _partizione(gruppi: tuple[str | None, ...]) -> tuple[int, ...]:
    """Per ogni campione l'indice del suo gruppo in ordine di apparizione; -1
    senza gruppo: due divisioni con la stessa partizione dividono i campioni
    negli stessi gruppi, con qualunque nome."""
    indice: dict[str, int] = {}
    return tuple(-1 if g is None else indice.setdefault(g, len(indice)) for g in gruppi)


def stessa_partizione(a: Applicazione, b: Applicazione) -> bool:
    """Vero se le due divisioni dividono i campioni negli stessi gruppi."""
    return _partizione(a.gruppi_dei_campioni) == _partizione(b.gruppi_dei_campioni)


def applica(divisione: Divisione, tabella: Tabella) -> tuple[Applicazione | None, list[str]]:
    """La divisione sulla tabella: la prima fonte che da' gruppi utilizzabili
    (:func:`_struttura`), e perche' le colonne provate prima non li danno.

    Le fonti si provano nell'ordine del catalogo; per ciascuna, le colonne
    dell'oggetto con quella chiave e almeno un valore, nell'ordine
    dell'oggetto. Per un incrocio, le due divisioni.
    """
    campioni = len(tabella.campioni)
    if divisione.incrocio:
        parti = []
        motivi: list[str] = []
        for nome in (divisione.righe, divisione.colonne):
            parte, perche = applica(_del_catalogo(nome), tabella)
            if parte is None:
                return None, [f"{nome}: " + ("; ".join(perche) or "nessuna colonna")]
            parti.append(parte)
        righe, colonne = parti
        if stessa_partizione(righe, colonne):
            return None, [f"{righe.nome} e {colonne.nome} dividono i campioni negli stessi "
                          "gruppi"]
        gruppi = tuple(None if r is None or c is None else f"{r} × {c}"
                       for r, c in zip(righe.gruppi_dei_campioni, colonne.gruppi_dei_campioni))
        return Applicazione(
            divisione.nome, f"griglia: righe {righe.nome} ({righe.regola}), colonne "
            f"{colonne.nome} ({colonne.regola})", "", gruppi, _conta(gruppi),
            {"griglia": (righe.nome, colonne.nome)}, (righe, colonne)), motivi
    motivi = []
    for fonte in divisione.fonti:
        for colonna, valori in tabella.valori.items():
            if chiave(colonna) != fonte.colonna or not any(_presente(v) for v in valori):
                continue
            app = _applica_fonte(divisione.nome, fonte, colonna, valori)
            motivo = _struttura(app, campioni)
            if motivo is None:
                return app, motivi
            motivi.append(f"{app.regola}: {motivo}")
    return None, motivi


def _annidate(righe: Applicazione, colonne: Applicazione) -> bool:
    """Vero se una delle due divisioni e' annidata nell'altra: nessun suo gruppo
    ha campioni in due gruppi dell'altra. Una griglia di due divisioni
    annidate ha una sola cella piena per riga, o per colonna: non le incrocia,
    e ``eco init`` non la propone (``validate`` la accetta, se richiesta)."""
    def dentro(a: Applicazione, b: Applicazione) -> bool:
        gruppi: dict[str, set[str]] = {}
        for x, y in zip(a.gruppi_dei_campioni, b.gruppi_dei_campioni):
            if x is not None and y is not None:
                gruppi.setdefault(x, set()).add(y)
        return all(len(v) == 1 for v in gruppi.values())
    return dentro(righe, colonne) or dentro(colonne, righe)


def _generica_per_prefisso(colonna: str, tabella: Tabella) -> Applicazione:
    """La divisione generica ``<colonna>_prefisso``."""
    fonte = Fonte(chiave(colonna), "prefisso", PREFISSO)
    return _applica_fonte(colonna + SUFFISSO_DEL_PREFISSO, fonte, colonna,
                          tabella.valori[colonna])


def _generica_per_colonna(colonna: str, tabella: Tabella) -> Applicazione:
    """La divisione generica di una colonna cosi' com'e'."""
    return _applica_fonte(colonna, Fonte(chiave(colonna), "colonna"), colonna,
                          tabella.valori[colonna])


# --------------------------------------------------------------------------- #
# Le proposte di eco init                                                      #
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class Proposta:
    """Una divisione che funziona sul dataset, per il file di configurazione."""

    applicazione: Applicazione
    #: La riga di spiegazione, con la colonna a cui si applica.
    spiegazione: str
    #: «catalogo» o «generica».
    origine: str

    @property
    def nome(self) -> str:
        return self.applicazione.nome


@dataclass(frozen=True)
class Proposte:
    """Le proposte di init con cio' che e' rimasto fuori."""

    proposte: tuple[Proposta, ...]
    #: Divisioni tolte perche' danno la stessa partizione di una proposta:
    #: (nome, nome della proposta uguale).
    doppioni: tuple[tuple[str, str], ...]
    #: Divisioni del catalogo con una colonna nell'oggetto che non funzionano:
    #: (nome, motivo).
    scartate: tuple[tuple[str, str], ...]


def proponi(tabella: Tabella, selezionabili: list[str],
            dimensione: int = DIMENSIONE_DELLE_PROPOSTE) -> Proposte:
    """Le divisioni che funzionano sulla tabella, senza doppioni, nell'ordine:
    il catalogo, poi per ogni colonna ``selezionabili`` (nell'ordine dato) la
    colonna stessa e il suo prefisso, poi gli incroci del catalogo fra le
    divisioni proposte che non sono annidate (:func:`_annidate`). Fra due
    divisioni con la stessa partizione resta la prima.
    """
    campioni = len(tabella.campioni)
    proposte: list[Proposta] = []
    doppioni: list[tuple[str, str]] = []
    scartate: list[tuple[str, str]] = []

    def aggiungi(app: Applicazione, spiegazione: str, origine: str) -> None:
        uguale = next((p for p in proposte if p.applicazione.componenti is None
                       and stessa_partizione(p.applicazione, app)), None)
        if uguale is not None:
            doppioni.append((app.nome, uguale.nome))
        else:
            proposte.append(Proposta(app, spiegazione, origine))

    for divisione in catalogo():
        if divisione.incrocio:
            continue
        app, motivi = applica(divisione, tabella)
        if app is None:
            if motivi:
                scartate.append((divisione.nome, "; ".join(motivi)))
            continue
        motivo = giudica(app, campioni, dimensione)
        if motivo:
            scartate.append((divisione.nome, f"{app.regola}: {motivo}"))
            continue
        if tabella.con_valore(divisione.nome) and app.espansione != {"colonna": divisione.nome}:
            # Il nome si leggerebbe come la colonna, che vince sulla divisione.
            scartate.append((divisione.nome, f"il nome e' quello di una colonna dell'oggetto "
                             f"con valori, che vince sulla divisione ({app.regola})"))
            continue
        aggiungi(app, f"{divisione.spiegazione} ({app.colonna})", "catalogo")
    for colonna in selezionabili:
        app = _generica_per_colonna(colonna, tabella)
        if giudica(app, campioni, dimensione) is None:
            aggiungi(app, "la colonna dell'oggetto cosi' com'e'", "generica")
        distinti = len({v for v in tabella.valori[colonna] if _presente(v)})
        if distinti <= _VALORI_PER_IL_PREFISSO or tabella.con_valore(
                colonna + SUFFISSO_DEL_PREFISSO):
            continue
        app = _generica_per_prefisso(colonna, tabella)
        if giudica(app, campioni, dimensione) is None and 2 * len(app.gruppi) <= len(
                {v for v, g in zip(tabella.valori[colonna], app.gruppi_dei_campioni)
                 if g is not None}):
            aggiungi(app, f"il prefisso dei valori di {colonna} (lettere iniziali seguite da "
                     "cifre)", "generica")
    nomi = {p.nome for p in proposte}
    for divisione in catalogo():
        if divisione.incrocio and {divisione.righe, divisione.colonne} <= nomi:
            app, _ = applica(divisione, tabella)
            if app is None or giudica(app, campioni, dimensione) is not None:
                continue
            if _annidate(*app.componenti):
                scartate.append((divisione.nome, f"{divisione.righe} e {divisione.colonne} sono "
                                 "annidate: la griglia non le incrocia"))
                continue
            proposte.append(Proposta(app, divisione.spiegazione, "catalogo"))
    return Proposte(tuple(proposte), tuple(doppioni), tuple(scartate))


# --------------------------------------------------------------------------- #
# L'espansione della configurazione                                            #
# --------------------------------------------------------------------------- #


def _nomi(dati: dict[str, Any]) -> list[tuple[str, str]]:
    """I nomi scritti nei campi :data:`CAMPI`, con il campo, nell'ordine."""
    design, comp, ord_ = (dati.get(g) or {} for g in ("design", "comp", "ord"))
    nomi = [("design.variable", design.get("variable"))]
    nomi += [("design.other_variables", n) for n in design.get("other_variables") or ()]
    nomi += [("comp.group_by", n) for n in comp.get("group_by") or ()]
    nomi += [("ord.color_by", n) for n in ord_.get("color_by") or ()]
    nomi.append(("comp.facet_wrap", comp.get("facet_wrap")))
    griglia = comp.get("facet_grid")
    if isinstance(griglia, dict):
        nomi += [("comp.facet_grid.rows", griglia.get("rows")),
                 ("comp.facet_grid.cols", griglia.get("cols"))]
    elif isinstance(griglia, str):
        nomi.append(("comp.facet_grid", griglia))
    return [(campo, nome) for campo, nome in nomi if isinstance(nome, str)]


def _da_espandere(nome: str, tabella: Tabella, gruppi_dichiarati: set[str]) -> bool:
    """Vero se il nome e' da leggere come divisione (regola 3 di :mod:`divisioni`)."""
    if nome in gruppi_dichiarati or tabella.con_valore(nome):
        return False
    if _del_catalogo(nome) is not None:
        return True
    radice = nome.removesuffix(SUFFISSO_DEL_PREFISSO)
    return nome.endswith(SUFFISSO_DEL_PREFISSO) and tabella.con_valore(radice)


@dataclass(frozen=True)
class Espansione:
    """La configurazione con le divisioni espanse, e che cosa ne e' venuto."""

    #: La configurazione espansa, in forma di dati.
    dati: dict[str, Any]
    #: Vero se ``dati`` differisce dalla configurazione di partenza.
    cambiata: bool
    #: Che cosa e' diventata ogni divisione nominata.
    scelte: tuple[dict[str, Any], ...]
    #: I rifiuti E-ECO-47 ed E-ECO-48 trovati qui.
    rifiuti: tuple[tuple[str, str], ...]
    #: I nomi respinti qui: le analisi non li trovano fra le colonne, e il loro
    #: rifiuto di colonna mancante (E-ECO-05) e' gia' detto.
    respinti: frozenset[str]
    #: Ogni nome scritto nei campi delle divisioni, con i campi.
    campi_dei_nomi: dict[str, tuple[str, ...]]


def espandi(dati: dict[str, Any], tabella: Tabella) -> Espansione:
    """La configurazione con le divisioni nominate espanse, e che cosa e'
    diventata ognuna.

    ``dati`` e' la configurazione validata, in forma di dati
    (``model_dump(by_alias=True)``). Una divisione con la regola «colonna»
    diventa la colonna dell'oggetto, in ogni campo in cui e' nominata; ogni
    altra diventa un gruppo derivato con il suo nome, aggiunto in fondo a
    ``design.groupings``; un incrocio in ``comp.facet_grid`` diventa
    ``{rows, cols}``. Senza divisioni i dati restano quelli di partenza.

    I rifiuti non si sollevano qui, perche' si dicono insieme a quelli delle
    analisi: ``E-ECO-48`` per una divisione che sul dataset non da' gruppi
    validi, ``E-ECO-47`` per un incrocio fuori da ``comp.facet_grid`` o un
    nome in ``comp.facet_grid`` che non e' un incrocio. Un nome che non e'
    nulla di noto resta com'e': lo respingono le analisi come colonna mancante,
    e :func:`traduci` ne fa un ``E-ECO-47``.
    """
    design = dati.get("design") or {}
    gruppi_dichiarati = set(design.get("groupings") or {})
    dimensione = (dati.get("stat") or {}).get("min_group_size") or DIMENSIONE_DELLE_PROPOSTE
    campioni = len(tabella.campioni)
    rifiuti: list[tuple[str, str]] = []
    respinti: set[str] = set()
    sostituzioni: dict[str, str] = {}
    nuovi_gruppi: dict[str, dict[str, Any]] = {}
    griglia: dict[str, str] | None = None
    togli_griglia = False
    scelte: dict[str, dict[str, Any]] = {}
    campi_dei_nomi: dict[str, list[str]] = {}

    def registra(nome: str, app: Applicazione, origine: str) -> None:
        """L'espansione di un'applicazione valida."""
        if "colonna" in app.espansione:
            sostituzioni[nome] = app.espansione["colonna"]
            esito = f"la colonna {app.espansione['colonna']}"
        else:
            sostituzioni[nome] = nome
            nuovi_gruppi[nome] = app.espansione["gruppo"]
            esito = f"il gruppo derivato {nome}"
        scelte.setdefault(nome, {"nome": nome, "origine": origine, "campi": [],
                                 "regola": app.regola, "colonna": app.colonna,
                                 "diventa": esito, "gruppi": dict(app.gruppi),
                                 "campioni_senza_gruppo": app.senza_gruppo})

    def valida_app(divisione: Divisione) -> tuple[Applicazione | None, str | None]:
        app, motivi = applica(divisione, tabella)
        if app is None:
            return None, "; ".join(motivi) or "nessuna colonna dell'oggetto a cui applicarla"
        return app, giudica(app, campioni, dimensione)

    for campo, nome in _nomi(dati):
        campi_dei_nomi.setdefault(nome, []).append(campo)
        divisione = _del_catalogo(nome)
        if campo == "comp.facet_grid":
            togli_griglia = True
            if nome in gruppi_dichiarati or tabella.con_valore(nome) or divisione is None \
                    or not divisione.incrocio:
                rifiuti.append(("E-ECO-47", (
                    f"comp.facet_grid: {nome!r} non e' un incrocio del catalogo; una griglia si "
                    "scrive {rows: <nome>, cols: <nome>} oppure con il nome di un incrocio. "
                    "Incroci del catalogo: " + ", ".join(
                        d.nome for d in catalogo() if d.incrocio) + ".")))
                continue
            app, motivo = valida_app(divisione)
            if motivo:
                rifiuti.append(("E-ECO-48", _non_valida(campo, divisione, motivo, dimensione)))
                continue
            for parte in app.componenti:
                if _da_espandere(parte.nome, tabella, gruppi_dichiarati):
                    registra(parte.nome, parte, "catalogo")
                    scelte[parte.nome]["campi"].append(f"{campo} (in {nome})")
                else:
                    sostituzioni.setdefault(parte.nome, parte.nome)
            righe, colonne = app.componenti
            griglia = {"rows": sostituzioni[righe.nome], "cols": sostituzioni[colonne.nome]}
            togli_griglia = False
            scelte[nome] = {"nome": nome, "origine": "catalogo", "campi": [campo],
                            "regola": app.regola, "colonna": "",
                            "diventa": f"la griglia rows {griglia['rows']}, cols "
                                       f"{griglia['cols']}",
                            "gruppi": dict(app.gruppi),
                            "campioni_senza_gruppo": app.senza_gruppo}
            continue
        if nome in scelte:
            scelte[nome]["campi"].append(campo)
            continue
        if nome in respinti or not _da_espandere(nome, tabella, gruppi_dichiarati):
            continue
        if divisione is None:
            radice = nome.removesuffix(SUFFISSO_DEL_PREFISSO)
            app = _generica_per_prefisso(radice, tabella)
            motivo, origine = giudica(app, campioni, dimensione), "generica"
            if motivo:
                respinti.add(nome)
                rifiuti.append(("E-ECO-48", (
                    f"{campo}: la divisione generica {nome!r} ({app.regola}) sul dataset non da' "
                    f"gruppi validi: {motivo}. Rimedio: nomina la colonna {radice} o dichiara un "
                    "gruppo derivato in design.groupings.")))
                continue
        elif divisione.incrocio:
            respinti.add(nome)
            rifiuti.append(("E-ECO-47", (
                f"{campo}: {nome!r} e' un incrocio di due divisioni e si usa solo come valore "
                f"di comp.facet_grid; qui va una delle due: {divisione.righe} o "
                f"{divisione.colonne}.")))
            continue
        else:
            app, motivo = valida_app(divisione)
            origine = "catalogo"
            if motivo:
                respinti.add(nome)
                rifiuti.append(("E-ECO-48", _non_valida(campo, divisione, motivo, dimensione)))
                continue
        registra(nome, app, origine)
        scelte[nome]["campi"].append(campo)

    # I campi del modello vincolato nominano variabili del disegno: seguono la
    # sostituzione, e si registrano fra i campi della divisione.
    vincolata = dati.get("constrained") or {}
    for campo in _CAMPI_CHE_SEGUONO:
        for nome in vincolata.get(campo.split(".")[1]) or ():
            if nome in scelte and campo not in scelte[nome]["campi"]:
                scelte[nome]["campi"].append(campo)
    cambiata = bool(scelte) or togli_griglia
    espansi = dati
    if cambiata:
        espansi = _sostituiti(dati, sostituzioni, nuovi_gruppi, griglia, togli_griglia)
    return Espansione(espansi, cambiata, tuple(scelte.values()), tuple(rifiuti),
                      frozenset(respinti), {n: tuple(c) for n, c in campi_dei_nomi.items()})


def _sostituiti(dati: dict[str, Any], sostituzioni: dict[str, str],
                nuovi_gruppi: dict[str, dict[str, Any]], griglia: dict[str, str] | None,
                togli_griglia: bool) -> dict[str, Any]:
    """I dati della configurazione con i nomi sostituiti, i gruppi derivati
    aggiunti in fondo e la griglia espansa (o tolta, se respinta)."""
    def via(valore: Any) -> Any:
        return sostituzioni.get(valore, valore) if isinstance(valore, str) else valore

    espansi = {g: (dict(v) if isinstance(v, dict) else v) for g, v in dati.items()}
    design = espansi["design"] = dict(espansi["design"])
    design["variable"] = via(design["variable"])
    design["other_variables"] = [via(n) for n in design.get("other_variables") or ()]
    design["groupings"] = {**(design.get("groupings") or {}), **nuovi_gruppi}
    comp = espansi["comp"] = dict(espansi["comp"])
    comp["group_by"] = [via(n) for n in comp.get("group_by") or ()]
    comp["facet_wrap"] = via(comp.get("facet_wrap"))
    if griglia is not None:
        comp["facet_grid"] = griglia
    elif togli_griglia:
        comp["facet_grid"] = None
    elif isinstance(comp.get("facet_grid"), dict):
        comp["facet_grid"] = {k: via(v) for k, v in comp["facet_grid"].items()}
    ord_ = espansi["ord"] = dict(espansi["ord"])
    ord_["color_by"] = [via(n) for n in ord_.get("color_by") or ()]
    if isinstance(espansi.get("constrained"), dict):
        vincolata = espansi["constrained"] = dict(espansi["constrained"])
        for nome in ("variables", "conditions"):
            if vincolata.get(nome):
                vincolata[nome] = [via(n) for n in vincolata[nome]]
    return espansi


_COLONNA_MANCANTE: Final = re.compile(r"^Colonna '(.*)' non trovata\. Colonne dell'oggetto: (.*)\.$",
                                      re.DOTALL)


def traduci(rifiuti: list[tuple[str, str]], espansione: Espansione) -> list[tuple[str, str]]:
    """I rifiuti delle analisi letti alla luce dell'espansione: la colonna
    mancante (``E-ECO-05``) di un nome scritto in un campo delle divisioni e'
    un nome sconosciuto (``E-ECO-47``), con le colonne e le divisioni; quella
    di un nome gia' respinto qui non si ripete."""
    tradotti = []
    for codice, dettaglio in rifiuti:
        trovato = _COLONNA_MANCANTE.match(dettaglio) if codice == "E-ECO-05" else None
        if trovato and trovato.group(1) in espansione.respinti:
            continue
        if trovato and trovato.group(1) in espansione.campi_dei_nomi:
            nome = trovato.group(1)
            campi = ", ".join(dict.fromkeys(espansione.campi_dei_nomi[nome]))
            tradotti.append(("E-ECO-47", (
                f"{campi}: {nome!r} non e' una colonna dell'oggetto, ne' un gruppo derivato di "
                "design.groupings, ne' una divisione del catalogo. Colonne dell'oggetto: "
                f"{trovato.group(2)}. Divisioni del catalogo: "
                f"{', '.join(d.nome for d in catalogo())}; per ogni colonna con valori anche "
                f"<colonna>{SUFFISSO_DEL_PREFISSO}.")))
            continue
        tradotti.append((codice, dettaglio))
    return tradotti


def _non_valida(campo: str, divisione: Divisione, motivo: str, dimensione: int) -> str:
    """Il dettaglio di E-ECO-48 per una divisione del catalogo."""
    return (f"{campo}: la divisione {divisione.nome!r} ({divisione.spiegazione}) sul dataset "
            f"non da' gruppi validi: {motivo}. Servono da {MINIMO_GRUPPI} a {MASSIMO_GRUPPI} "
            f"gruppi, almeno {MINIMO_GRUPPI} con {dimensione} campioni (stat.min_group_size), "
            "e un gruppo per almeno meta' dei campioni. Rimedio: scegli un'altra divisione fra "
            "quelle che eco init propone, nomina una colonna dell'oggetto, o dichiara un "
            "gruppo derivato in design.groupings.")
