"""Il report di esecuzione: un documento HTML ricavato dalla cartella di output.

Il report non è una fase: non entra nel grafo e non ha un manifesto. Si ricava
da ciò che un'esecuzione ha lasciato su disco (manifesti di fase, artefatti,
configurazioni registrate in ``00_config``, registro degli avvii in ``99_logs``),
senza rieseguire alcun calcolo, senza R e senza leggere i dati grezzi.

**È funzione del contenuto della cartella.** Non riporta date di generazione né
valori dell'ambiente: generato due volte dalla stessa esecuzione dà gli stessi
byte, e le date che vi compaiono sono quelle dei manifesti e del registro. L'unico
confronto con qualcosa che sta fuori dalla cartella è quello di provenienza: il
sorgente registrato nei manifesti contro il sorgente della pipeline che genera
il documento, perché un avviso di provenienza è per definizione quel confronto.

**Non altera lo stato.** Scrive solo nella cartella ``report`` sotto
``io.out_root``, che nessuna fase legge e nessun manifesto elenca; non scrive
nel log. Lo stato delle fasi che riporta è quello dei manifesti: la validità
rispetto alla configurazione e ai dati di adesso la giudica ``resume``.

**Ogni artefatto letto è verificato** contro il checksum del manifesto della
sua fase: un artefatto alterato non viene riportato, e l'alterazione compare in
apertura. Le decisioni automatiche sono lette dagli artefatti che le
registrano (manifesti, ``soglia.json``, ``decontam_riepilogo.json``,
``esclusioni.tsv``), non ricostruite.

Ogni tabella del documento è scritta anche come file TSV in ``report/tabelle``.
"""

from __future__ import annotations

import base64
import csv
import html
import io
import json
import statistics
from collections import defaultdict
from collections.abc import Iterable, Mapping
from dataclasses import dataclass, field
from pathlib import Path
from string import Template
from typing import Any, Final

import yaml

from amplicon16s import __version__
from amplicon16s.config.defaults import FATTI_OSD734
from amplicon16s.config.resolve import versioni_registrate
from amplicon16s.gates.registry import REGISTRO
from amplicon16s.io_layer.artifacts import (
    AlberoOutput,
    Fase,
    ManifestoPasso,
    ManifestoPassoNonValido,
    scrivi_atomico,
)
from amplicon16s.io_layer.checksums import corrisponde
from amplicon16s.logging.logger import NOME_FILE_AVVII
from amplicon16s.metadata.lettura_inventario import NOME_CROSSWALK
from amplicon16s.metadata.models import CLASSI_CONTROLLATE, ClasseCampione
from amplicon16s.rbridge.runner import cartella_r
from amplicon16s.runner.executor import (
    EVENTO_CONCLUSIONE,
    EVENTO_CONFIGURAZIONE,
    EVENTO_CONTROLLI,
    EVENTO_FASE,
    NOME_PUNTO_DI_RIPRESA,
)
from amplicon16s.runner.graph import GRAFO, Passo
from amplicon16s.runner.project import avviso_di_provenienza, passi_realizzati
from amplicon16s.runner.provenienza import file_del_sorgente, impronta_sorgente
from amplicon16s.runner.tracciamento import PREFISSO, Tracciamento, componi
from amplicon16s.steps.base import PipelineStep
from amplicon16s.steps.s01_profile import NOME_LUNGHEZZE, NOME_QUALITA

__all__ = [
    "CARTELLA_REPORT",
    "CARTELLA_TABELLE",
    "NOME_REPORT",
    "Report",
    "Tabella",
    "costruisci",
    "genera",
    "scrivi",
    "troncamento_suggerito",
]

#: La cartella del report sotto ``io.out_root``. Non è una cartella di fase:
#: nessun manifesto la elenca e la valutazione dello stato non la guarda.
CARTELLA_REPORT: Final = "report"
NOME_REPORT: Final = "report.html"
#: Le tabelle del documento come file TSV, sotto la cartella del report.
CARTELLA_TABELLE: Final = "tabelle"

_MODELLI: Final = Path(__file__).resolve().parent / "templates"
_CLASSI: Final = tuple(c.value for c in ClasseCampione)
#: Le due regole del troncamento suggerito: la quota di letture che un
#: troncamento puo' scartare perche' piu' corte, e la qualita' mediana sotto
#: la quale una posizione non conviene tenerla.
QUOTA_CORTE_SUGGERITA: Final = 0.05
QUALITA_MINIMA_SUGGERITA: Final = 30


# --------------------------------------------------------------------------- #
# Forma del testo                                                              #
# --------------------------------------------------------------------------- #


def _e(valore: Any) -> str:
    """Il testo di un valore, pronto per l'HTML."""
    return html.escape(str(valore), quote=True)


def _n(valore: Any) -> str:
    """Un numero con il punto come separatore delle migliaia; un valore non
    intero, come una mediana fra due conteggi, con una cifra decimale.
    """
    if valore == int(valore):
        return f"{int(valore):,}".replace(",", ".")
    return f"{valore:,.1f}".replace(",", " ").replace(".", ",").replace(" ", ".")


def _righe(numero: int) -> str:
    """Il numero di righe di una tabella, concordato."""
    return "1 riga" if numero == 1 else f"{_n(numero)} righe"


def _ordine_piastre(piastra: str) -> tuple[bool, int, str]:
    """Le piastre numerate in ordine numerico, le altre dopo."""
    return (not piastra.isdigit(), int(piastra) if piastra.isdigit() else 0, piastra)


def _pct(frazione: Any) -> str:
    """Una frazione come percentuale con una cifra decimale e la virgola; un
    trattino se la frazione non esiste (una classe senza campioni o senza
    letture: un dataset senza controlli non ha la loro quota).
    """
    if frazione is None:
        return "-"
    return f"{100 * float(frazione):.1f}%".replace(".", ",")


def _valore(valore: Any) -> str:
    """Il valore di un parametro come lo si scriverebbe nella configurazione."""
    return valore if isinstance(valore, str) else json.dumps(valore, ensure_ascii=False)


def _corta(impronta: Any) -> str:
    """Le prime dodici cifre di un'impronta, senza il nome dell'algoritmo."""
    return str(impronta).removeprefix("sha256:")[:12] if impronta else ""


def _si_no(valore: Any) -> str:
    """Un valore di verità in parole; vuoto se non registrato."""
    return "" if valore is None else ("sì" if valore else "no")


def _esito(voce: Mapping[str, Any]) -> str:
    """L'esito di un controllo: non eseguito non equivale a superato."""
    if not voce.get("eseguito", True):
        return "non eseguito"
    return "superato" if voce["superato"] else "NON SUPERATO"


def _elenco(voci: Iterable[Mapping[str, Any]]) -> str:
    """Violazioni o avvisi di un controllo su una riga: codice e dettaglio."""
    return "; ".join(f"{v['codice']}: {v['dettaglio']}" for v in voci)


# --------------------------------------------------------------------------- #
# Tabelle                                                                      #
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class Tabella:
    """Una tabella del report: compare nel documento ed è scritta come TSV."""

    #: Nome del file TSV, senza estensione.
    nome: str
    titolo: str
    colonne: tuple[str, ...]
    #: I conteggi restano numeri: nel documento prendono il separatore delle
    #: migliaia, nel file TSV no, perché lì servono a chi li rielabora.
    righe: tuple[tuple[Any, ...], ...]
    #: Nel documento sta ripiegata: è lunga, e serve a chi cerca un caso.
    ripiegata: bool = False
    #: Righe mostrate nel documento; il file TSV le contiene tutte.
    limite: int | None = None

    def tsv(self) -> str:
        """La tabella come testo separato da tabulazioni, con fine riga POSIX."""
        testo = io.StringIO()
        scrittore = csv.writer(testo, delimiter="\t", lineterminator="\n")
        scrittore.writerow(self.colonne)
        scrittore.writerows(self.righe)  # i numeri senza separatore delle migliaia
        return testo.getvalue()

    def html(self) -> str:
        """La tabella nel documento, con il rimando al file TSV."""
        righe = self.righe if self.limite is None else self.righe[: self.limite]
        corpo = "\n".join(
            "<tr>" + "".join(f"<td>{_e(_cella(c))}</td>" for c in riga) + "</tr>" for riga in righe
        )
        testa = "".join(f"<th>{_e(c)}</th>" for c in self.colonne)
        mostrate = (
            f"le prime {_n(len(righe))} di {_righe(len(self.righe))}"
            if len(righe) < len(self.righe)
            else _righe(len(self.righe))
        )
        percorso = f"{CARTELLA_TABELLE}/{self.nome}.tsv"
        tabella = (
            f'<div class="scorre"><table>\n<caption>{_e(self.titolo)}</caption>\n'
            f"<thead><tr>{testa}</tr></thead>\n<tbody>\n{corpo}\n</tbody>\n</table></div>\n"
            f'<p class="sorgente">{mostrate}; tabella completa in '
            f'<a href="{_e(percorso)}"><code>{_e(percorso)}</code></a></p>'
        )
        if not self.ripiegata:
            return tabella
        return (
            f"<details><summary>{_e(self.titolo)} ({_righe(len(self.righe))})</summary>\n"
            f"{tabella}\n</details>"
        )


def _cella(valore: Any) -> str:
    """Il contenuto di una cella nel documento: i numeri con il separatore."""
    if isinstance(valore, (int, float)) and not isinstance(valore, bool):
        return _n(valore)
    return str(valore)


def _tabella(nome: str, titolo: str, colonne: Iterable[str], righe: Iterable[Iterable[Any]],
             **opzioni: Any) -> Tabella:
    """Una :class:`Tabella` dalle sue righe; un valore assente è una cella vuota."""
    return Tabella(
        nome, titolo, tuple(colonne),
        tuple(tuple("" if c is None else c for c in riga) for riga in righe),
        **opzioni,
    )


@dataclass(frozen=True)
class Report:
    """Il report costruito: il documento e le sue tabelle."""

    html: str
    tabelle: tuple[Tabella, ...]
    #: Le segnalazioni in apertura, una riga ciascuna, per chi lo genera.
    segnalazioni: tuple[str, ...]


@dataclass
class _Sezione:
    """Una sezione del documento in costruzione."""

    ident: str
    titolo: str
    blocchi: list[str] = field(default_factory=list)
    tabelle: list[Tabella] = field(default_factory=list)

    def testo(self, paragrafo: str, classe: str | None = None) -> None:
        """Aggiunge un paragrafo; il testo è già HTML."""
        attributo = f' class="{classe}"' if classe else ""
        self.blocchi.append(f"<p{attributo}>{paragrafo}</p>")

    def sottotitolo(self, titolo: str) -> None:
        """Aggiunge un titolo di sottosezione."""
        self.blocchi.append(f"<h3>{_e(titolo)}</h3>")

    def sintesi(self, voci: Iterable[tuple[str, Any]]) -> None:
        """Aggiunge un elenco di coppie nome e valore."""
        righe = "".join(f"<dt>{_e(n)}</dt><dd>{_e(v)}</dd>" for n, v in voci)
        self.blocchi.append(f'<dl class="sintesi">{righe}</dl>')

    def tabella(self, tabella: Tabella) -> None:
        """Aggiunge una tabella al documento e ai file TSV."""
        self.tabelle.append(tabella)
        self.blocchi.append(tabella.html())

    def html(self) -> str:
        """La sezione nel documento."""
        return (
            f'<section id="{self.ident}">\n<h2>{_e(self.titolo)}</h2>\n'
            + "\n".join(self.blocchi) + "\n</section>"
        )


# --------------------------------------------------------------------------- #
# Lettura della cartella di output                                             #
# --------------------------------------------------------------------------- #


class _Esecuzione:
    """Ciò che un'esecuzione ha lasciato su disco, letto senza modificarlo.

    Una fase conta come conclusa se il suo manifesto esiste e la sua impronta
    corrisponde al contenuto. Un artefatto si legge solo attraverso il
    manifesto della sua fase, che ne dà il checksum.
    """

    def __init__(self, radice: Path, passi: Mapping[Passo, PipelineStep]) -> None:
        self.albero = AlberoOutput(radice)
        self.passi = dict(passi)
        self.manifesti: dict[Passo, ManifestoPasso] = {}
        self.non_validi: dict[Passo, str] = {}
        #: Artefatti che non corrispondono al manifesto: chiave per non ripeterli.
        self.alterati: dict[str, None] = {}
        for passo, fase in self.passi.items():
            try:
                manifesto = self.albero.manifesto_passo(passo, fase.cartella)
            except ManifestoPassoNonValido as errore:
                self.non_validi[passo] = str(errore)
                continue
            if manifesto is not None:
                self.manifesti[passo] = manifesto
        self.versioni = [
            {"file": p.name, **yaml.safe_load(p.read_text(encoding="utf-8"))}
            for p in versioni_registrate(radice)
        ]
        #: I parametri della configurazione in uso, come ``gruppo.parametro``.
        self.parametri: dict[str, Any] = {
            f"{gruppo}.{chiave}": valore
            for gruppo, valori in (self.versioni[-1]["parametri"] if self.versioni else {}).items()
            for chiave, valore in valori.items()
        }
        eventi = _eventi(self.albero.cartella(Fase.LOGS))
        #: Falso se l'esecuzione non ha il registro degli avvii: è stata
        #: prodotta prima che esistesse, o senza la riga di comando.
        self.registro_avvii = eventi is not None
        self.avvii = _avvii(eventi or [])

    def byte(self, passo: Passo, nome: str) -> bytes | None:
        """Un artefatto di una fase conclusa, se corrisponde al suo checksum."""
        manifesto = self.manifesti.get(passo)
        voce = next((v for v in manifesto.artefatti if v["nome"] == nome), None) if manifesto else None
        if manifesto is None or voce is None:
            return None
        percorso = self.albero.cartella(manifesto.fase) / nome
        if not percorso.is_file():
            if nome not in self.albero.rimossi_del_passo(manifesto):
                self.alterati[f"{passo}: {nome} manca dalla cartella {manifesto.fase.value}"] = None
            return None
        if not corrisponde(percorso, str(voce["checksum"])):
            self.alterati[f"{passo}: {nome} non corrisponde al checksum del manifesto"] = None
            return None
        return percorso.read_bytes()

    def json(self, passo: Passo, nome: str) -> Any:
        """Un artefatto JSON, o ``None`` se la fase non lo ha o non è integro."""
        dati = self.byte(passo, nome)
        return None if dati is None else json.loads(dati)

    def tsv(self, passo: Passo, nome: str) -> list[dict[str, str]] | None:
        """Le righe di un artefatto TSV, o ``None`` se non c'è o non è integro."""
        dati = self.byte(passo, nome)
        if dati is None:
            return None
        return list(csv.DictReader(io.StringIO(dati.decode("utf-8"), newline=""), delimiter="\t"))

    def tracciamento(self) -> Tracciamento:
        """Il tracciamento ricomposto dalle fasi concluse con tabelle integre."""
        concluse = [
            fase for passo, fase in self.passi.items()
            if passo in self.manifesti
            and all(self.byte(passo, f"{PREFISSO}{p}.tsv") is not None for p in fase.passi_tracciamento)
        ]
        return componi(self.albero, concluse)

    def coerenza(self, passo: Passo) -> list[str]:
        """Le fasi a monte ricalcolate dopo che ``passo`` le aveva lette.

        Il manifesto registra l'impronta di ogni fase di cui la fase ha
        consumato gli artefatti: se oggi quella fase ha un'altra impronta, o
        non è conclusa, questo risultato è stato calcolato su artefatti che
        non sono più quelli su disco.
        """
        registrate = self.manifesti[passo].calcolata_su.get("a_monte", {})
        return [
            monte for monte, impronta in registrate.items()
            if Passo(monte) not in self.manifesti or self.manifesti[Passo(monte)].impronta != impronta
        ]


def _eventi(cartella: Path) -> list[dict[str, Any]] | None:
    """Gli eventi del registro degli avvii, dal più vecchio; ``None`` se manca.

    Il registro è solo in aggiunta e non ruota: contiene ogni avvio
    dell'esecuzione. Il log, che ruota, non viene letto al suo posto: darebbe
    una storia che può essere incompleta senza poterlo dire.
    """
    registro = cartella / NOME_FILE_AVVII
    if not registro.is_file():
        return None
    eventi = []
    for riga in registro.read_text(encoding="utf-8").splitlines():
        try:
            evento = json.loads(riga)
        except ValueError:
            continue  # una riga troncata da un'interruzione non ferma il report
        if isinstance(evento, dict):
            eventi.append(evento)
    return eventi


def _avvii(eventi: Iterable[Mapping[str, Any]]) -> list[dict[str, Any]]:
    """Gli avvii dell'esecutore nel registro degli avvii, ciascuno con i suoi controlli,
    la configurazione con cui ha eseguito, le fasi eseguite e la conclusione.
    """
    avvii: list[dict[str, Any]] = []

    def nuovo(evento: Mapping[str, Any]) -> dict[str, Any]:
        avvii.append({
            "istante": evento.get("istante", ""), "controlli": None, "fino_a": None,
            "blocco_r": None, "rigorosa": None, "configurazione": None, "fasi": [],
            "conclusione": None,
        })
        return avvii[-1]

    for evento in eventi:
        tipo = evento.get("evento")
        aperto = avvii[-1] if avvii and avvii[-1]["conclusione"] is None else None
        if tipo == EVENTO_CONTROLLI:
            avvio = nuovo(evento)
            avvio.update(
                controlli=evento["esiti"], fino_a=evento.get("fino_a"),
                blocco_r=evento.get("blocco_r"), rigorosa=evento.get("provenienza_rigorosa"),
            )
        elif tipo == EVENTO_CONFIGURAZIONE:
            avvio = aperto if aperto and aperto["configurazione"] is None else nuovo(evento)
            avvio["configurazione"] = evento
        elif tipo == EVENTO_FASE:
            (aperto or nuovo(evento))["fasi"].append((evento["passo"], evento["esito"]))
        elif tipo == EVENTO_CONCLUSIONE:
            conclusione = evento["conclusione"]
            if evento.get("codice"):
                conclusione += f" ({evento['codice']})"
            (aperto or nuovo(evento))["conclusione"] = conclusione
    return avvii


# --------------------------------------------------------------------------- #
# Sezioni                                                                      #
# --------------------------------------------------------------------------- #


def _rigorosa(esito: Mapping[str, Any] | None) -> str:
    """La regola rigorosa sulla provenienza in un avvio, in una riga: che cosa
    e' stato verificato e che cosa e' solo dichiarato.
    """
    if esito is None:
        return "non registrata"
    if not esito.get("attiva"):
        return "non attiva"
    git, ambiente = esito.get("git", {}), esito.get("ambiente_r", {})
    if not git.get("verificato"):
        return "rifiutata: codice non identificato da un commit"
    if not ambiente.get("verificato"):
        return "rifiutata: ambiente R diverso dal file di blocco"
    return (
        f"verificati il codice (commit {str(git.get('commit'))[:12]}, senza modifiche) e "
        f"l'ambiente R ({ambiente.get('r')}, {ambiente.get('pacchetti_confrontati')} pacchetti "
        "uguali al file di blocco); immagine dichiarata, non verificabile"
    )


def _origine(chiave: str, versione: Mapping[str, Any]) -> str:
    """Da dove viene il valore di un parametro nella configurazione registrata."""
    if chiave in versione.get("derivati", ()):
        return "derivato"
    dichiarati = versione.get("dichiarati")
    if dichiarati is None:
        return "non registrata"
    return "dichiarato" if chiave in dichiarati else "predefinito"


def _aggiustamenti(esecuzione: _Esecuzione) -> list[tuple[Passo, dict[str, Any]]]:
    """I tentativi ripetuti registrati nei manifesti, con la fase."""
    return [(p, a) for p, m in esecuzione.manifesti.items() for a in m.aggiustamenti]


def _valori_osd734(esecuzione: _Esecuzione) -> Tabella:
    """I parametri che in questa esecuzione valgono il predefinito tarato su
    OSD-734 senza essere stati dichiarati, con il fatto che lo giustificava.

    I parametri obbligatori non possono essere ereditati (la configurazione
    senza di essi non parte) e non vi compaiono. Un parametro tarato dichiarato
    nel file, anche con lo stesso valore, e' una scelta e non compare. Se la
    configurazione registrata non dice quali parametri sono stati dichiarati
    (esecuzioni di versioni precedenti), compaiono tutti quelli che hanno il
    valore di OSD-734.
    """
    versione = esecuzione.versioni[-1] if esecuzione.versioni else {}
    righe = []
    for chiave, riferimento in FATTI_OSD734.items():
        if chiave not in esecuzione.parametri:
            continue
        in_uso = esecuzione.parametri[chiave]
        # La configurazione registrata è YAML: una tupla vi diventa un elenco.
        atteso = list(riferimento.valore) if isinstance(riferimento.valore, tuple) else riferimento.valore
        origine = _origine(chiave, versione)
        if origine == "dichiarato" or in_uso != atteso:
            continue
        righe.append((chiave, _valore(in_uso), origine, riferimento.fatto))
    return _tabella(
        "valori_osd734",
        "Parametri non dichiarati che valgono il predefinito tarato su OSD-734",
        ("parametro", "valore in uso", "origine", "fatto accertato su OSD-734"),
        righe,
    )


def _sezione_stato(esecuzione: _Esecuzione) -> tuple[_Sezione, bool, list[str]]:
    """Lo stato di ogni fase secondo i manifesti; se la catena è completa; e le
    fasi calcolate su artefatti a monte che non sono più quelli su disco.
    """
    sezione = _Sezione("stato", "Stato della catena")
    ultima_versione = {
        passo: avvio["configurazione"]["file"]
        for avvio in esecuzione.avvii if avvio["configurazione"]
        for passo, _ in avvio["fasi"]
    }
    righe, completa, superate = [], True, []
    for passo in GRAFO.ordine():
        nodo = GRAFO.nodo(passo)
        manifesto = esecuzione.manifesti.get(passo)
        conclusa: Any = ""
        artefatti: Any = ""
        if nodo.parametro_attivazione and not esecuzione.parametri.get(nodo.parametro_attivazione):
            stato = f"disattivata ({nodo.parametro_attivazione} falso)"
        elif passo not in esecuzione.passi:
            stato, completa = "non realizzata", False
        elif passo in esecuzione.non_validi:
            stato, completa = f"manifesto non valido: {esecuzione.non_validi[passo]}", False
        elif manifesto is None:
            stato, completa = "non conclusa", False
        else:
            conclusa, artefatti = manifesto.conclusa, len(manifesto.artefatti)
            ricalcolate = esecuzione.coerenza(passo)
            rimossi = len(esecuzione.albero.rimossi_del_passo(manifesto))
            stato = "conclusa" + (f"; {_n(rimossi)} artefatti rimossi di proposito" if rimossi else "")
            if ricalcolate:
                stato = f"conclusa su artefatti superati: ricalcolate dopo {', '.join(ricalcolate)}"
                completa = False
                superate.append(f"{passo} (a monte: {', '.join(ricalcolate)})")
        righe.append((
            passo, nodo.descrizione, nodo.cartella.value, stato, conclusa, artefatti,
            ultima_versione.get(str(passo), ""),
        ))
    sezione.testo(
        "Lo stato è quello registrato nei manifesti di fase: una fase è conclusa se il suo "
        "manifesto esiste e la sua impronta corrisponde al contenuto. La validità rispetto "
        "alla configurazione e ai dati correnti è giudicata dal sottocomando "
        "<code>resume</code>, che questo documento non sostituisce."
    )
    sezione.tabella(_tabella(
        "fasi", "Fasi del grafo",
        ("fase", "descrizione", "cartella", "stato", "conclusa (manifesto)", "artefatti",
         "configurazione dell'ultima esecuzione (registro degli avvii)"),
        righe,
    ))
    punto = esecuzione.albero.cartella(Fase.LOGS) / f"{NOME_PUNTO_DI_RIPRESA}.json"
    if punto.is_file():
        documento = json.loads(punto.read_text(encoding="utf-8"))
        sezione.sottotitolo("Punto di ripresa dichiarato")
        sezione.sintesi((nome, documento.get(nome) or "") for nome in (
            "origine", "passo", "codice", "categoria", "sintesi", "dettaglio", "motivo",
            "tentativi", "azione", "comando",
        ))
    return sezione, completa, superate


def _sezione_gate(esecuzione: _Esecuzione) -> _Sezione:
    """I quindici gate come li ha visti S0, e i controlli di avvio di ogni avvio."""
    sezione = _Sezione("gate", "Gate di validazione e controlli di avvio")
    gate = esecuzione.json(Passo.S0, "gates.json")
    if gate is None:
        sezione.testo("S0 non è conclusa: l'esito dei quindici gate non è fra gli artefatti.")
    else:
        superati = sum(g["superato"] for g in gate["gate"])
        sezione.testo(
            f"Esito di S0: {superati} gate superati su {len(gate['gate'])}, nell'ordine di "
            "esecuzione (G15, la coerenza della configurazione, apre la sequenza)."
        )
        sezione.tabella(_tabella(
            "gate", "Esito dei quindici gate in S0 (01_input_validation/gates.json)",
            ("gate", "descrizione", "esito", "violazioni", "avvisi"),
            ((g["gate"], g["descrizione"], _esito(g), _elenco(g["violazioni"]), _elenco(g["avvisi"]))
             for g in gate["gate"]),
        ))

    sezione.sottotitolo("Controlli di avvio di ogni esecuzione")
    sezione.testo(
        "G15 (coerenza della configurazione), G14 (risorse della macchina) e G12 "
        "(riferimento tassonomico contro <code>tax.ref_md5</code>) sono ripetuti "
        "dall'esecutore a ogni avvio e a ogni ripresa, anche con S0 già conclusa: "
        "verificano la macchina e i file in uso in quel momento, non i dati."
    )
    descrizioni = {v.nome: v.descrizione for v in REGISTRO}
    registrati = [(i, a) for i, a in enumerate(esecuzione.avvii, 1) if a["controlli"] is not None]
    if not registrati:
        sezione.testo(
            "Questa esecuzione non ha il registro degli avvii "
            f"(<code>{Fase.LOGS.value}/{NOME_FILE_AVVII}</code>): è stata prodotta da una "
            "versione della pipeline che non lo scriveva, o senza la riga di comando. I "
            "controlli di avvio non sono ricostruiti dal log, che ruota e potrebbe darne una "
            "storia incompleta."
            if not esecuzione.registro_avvii else
            "Il registro degli avvii non riporta alcun controllo di avvio.",
            "nota",
        )
    else:
        sezione.tabella(_tabella(
            "controlli_di_avvio", "Controlli di avvio, per avvio (99_logs/avvii.jsonl)",
            ("avvio", "istante", "gate", "descrizione", "esito", "violazioni"),
            ((numero, avvio["istante"], c["gate"], descrizioni.get(c["gate"], ""), _esito(c),
              _elenco(c["violazioni"]))
             for numero, avvio in registrati for c in avvio["controlli"]),
        ))
    return sezione


def _sezione_configurazione(esecuzione: _Esecuzione) -> _Sezione:
    """Avvii, versioni della configurazione registrata e origine di ogni parametro."""
    sezione = _Sezione("configurazione", "Configurazione: esecuzioni, versioni e origine dei valori")
    # Gli avvii prima delle versioni: un avvio rifiutato dai controlli non
    # registra alcuna configurazione, e deve restare leggibile lo stesso.
    sezione.sottotitolo("Avvii dell'esecuzione")
    if esecuzione.avvii:
        sezione.tabella(_tabella(
            "esecuzioni", "Avvii dell'esecutore (99_logs/avvii.jsonl)",
            ("avvio", "istante", "fino a", "controlli di avvio",
             "regola rigorosa sulla provenienza", "configurazione", "digest",
             "fasi eseguite", "conclusione"),
            ((numero, a["istante"], a["fino_a"] or "",
              "non registrati" if a["controlli"] is None
              else " ".join(f"{c['gate']}: {_esito(c)};" for c in a["controlli"]).rstrip(";"),
              _rigorosa(a["rigorosa"]),
              a["configurazione"]["file"] if a["configurazione"] else "non registrata",
              _corta(a["configurazione"]["digest"]) if a["configurazione"] else "",
              ", ".join(p for p, _ in a["fasi"]) or "nessuna",
              a["conclusione"] or "non registrata")
             for numero, a in enumerate(esecuzione.avvii, 1)),
        ))
    else:
        sezione.testo(
            "Il registro degli avvii è assente o vuoto: gli avvii, e le fasi eseguite con "
            "ciascuna versione della configurazione, non sono determinabili, e non vengono "
            "ricostruiti dal log.", "nota",
        )

    if not esecuzione.versioni:
        sezione.testo("Nessuna configurazione registrata in <code>00_config</code>.")
        return sezione

    sezione.sottotitolo("Versioni della configurazione registrata")
    sezione.testo(
        "Ogni avvio con un digest diverso dall'ultimo registrato scrive una nuova versione "
        "accanto alle precedenti, senza sovrascriverle. La configurazione in uso è l'ultima."
    )
    eseguite: dict[str, list[str]] = {}
    for numero, avvio in enumerate(esecuzione.avvii, 1):
        if avvio["configurazione"] and avvio["fasi"]:
            eseguite.setdefault(avvio["configurazione"]["file"], []).append(
                f"avvio {numero}: {', '.join(p for p, _ in avvio['fasi'])}"
            )
    sezione.tabella(_tabella(
        "configurazioni", "Versioni registrate in 00_config",
        ("file", "digest", "precedente", "differenze dalla precedente", "fasi eseguite (registro degli avvii)"),
        ((v["file"], v["digest"], v.get("precedente", ""),
          "; ".join(v.get("differenze_dalla_precedente", [])),
          "; ".join(eseguite.get(v["file"], [])) or ("nessuna" if esecuzione.avvii else "non determinabili"))
         for v in esecuzione.versioni),
    ))

    versione = esecuzione.versioni[-1]
    aggiustati: dict[str, list[str]] = {}
    for passo, a in _aggiustamenti(esecuzione):
        if a["parametro"]:
            aggiustati.setdefault(a["parametro"], []).append(
                f"{passo}: {_valore(a['precedente'])} -> {_valore(a['usato'])} dopo {a['codice']}"
            )
    righe = [
        (chiave, _valore(valore), _origine(chiave, versione), "; ".join(aggiustati.get(chiave, [])))
        for chiave, valore in esecuzione.parametri.items()
    ]
    conteggi = {o: sum(r[2] == o for r in righe) for o in ("dichiarato", "predefinito", "derivato", "non registrata")}
    sezione.sottotitolo("Parametri effettivamente usati")
    sezione.sintesi((
        ("configurazione in uso", versione["file"]),
        ("digest", versione["digest"]),
        ("dichiarati nel file di configurazione", _n(conteggi["dichiarato"])),
        ("presi dal predefinito della pipeline", _n(conteggi["predefinito"])),
        ("derivati da altri parametri", _n(conteggi["derivato"])),
        ("aggiustati da un tentativo ripetuto", _n(len(aggiustati))),
    ))
    sezione.testo(
        "L'origine viene dallo schema, che registra quali campi il file di configurazione "
        "ha impostato: non da un confronto fra valori. Un parametro <em>dichiarato</em> è "
        "scritto nel file; uno <em>predefinito</em> vi è assente e vale il predefinito della "
        "pipeline; uno <em>derivato</em> è calcolato da altri parametri e non è ammesso nel "
        "file. Un valore <em>aggiustato</em> è quello usato da una fase dopo un tentativo "
        "ripetuto: il valore dichiarato resta nella colonna del valore. Dichiarato non "
        "significa scelto: una configurazione copiata da un'istanza completa dichiara anche "
        "i valori che nessuno ha rivisto, ed è il motivo della segnalazione in apertura."
    )
    if conteggi["non registrata"]:
        sezione.testo(
            "Questa configurazione è stata registrata da una versione della pipeline che non "
            "annotava i parametri dichiarati: per i parametri non derivati l'origine non è "
            "determinabile dal contenuto della cartella, e non viene ricostruita.", "nota",
        )
    sezione.tabella(_tabella(
        "parametri", f"Parametri della configurazione in uso ({versione['file']})",
        ("parametro", "valore", "origine", "aggiustamento"), righe, ripiegata=True,
    ))
    return sezione


def _sezione_provenienza(esecuzione: _Esecuzione) -> tuple[_Sezione, list[tuple[str, str]]]:
    """La provenienza registrata di ogni fase conclusa, e gli avvisi."""
    sezione = _Sezione("provenienza", "Provenienza delle fasi")
    immagine = esecuzione.parametri.get("run.container")
    righe, avvisi = [], []
    for passo, manifesto in esecuzione.manifesti.items():
        registrata = manifesto.provenienza
        try:
            sorgente, per_file = impronta_sorgente(
                file_del_sorgente(esecuzione.passi[passo], cartella_r())
            )
            avviso = avviso_di_provenienza(
                registrata, {"sorgente": sorgente, "file": per_file, "immagine": immagine}
            )
        except (OSError, ValueError, ImportError) as errore:
            avviso = f"confronto non eseguibile: {errore}"
        if avviso:
            avvisi.append((str(passo), avviso))
        righe.append((
            passo, registrata.get("versione", ""), _corta(registrata.get("sorgente")),
            registrata.get("commit") or "non registrato",
            _si_no(registrata.get("modifiche_non_committate")),
            registrata.get("immagine", ""), manifesto.conclusa, avviso or "",
        ))
    sezione.testo(
        "Per ogni fase conclusa il manifesto registra la versione del calcolo, l'impronta "
        "del sorgente effettivamente eseguito (moduli Python, script R e funzioni "
        "condivise), il commit del repository e l'immagine dichiarata in "
        "<code>run.container</code>. Il commit non è registrato quando il codice non sta in "
        "un repository git, come nell'immagine. L'immagine è una dichiarazione: dall'interno "
        "del container il suo digest non è verificabile. Con la regola rigorosa sulla "
        "provenienza (<code>run.strict_provenance</code>, riportata per ogni avvio nella "
        'sezione <a href="#configurazione">Configurazione</a>) il codice e l\'ambiente R '
        "sono invece verificati a ogni avvio. L'avviso confronta il sorgente e l'immagine "
        f"registrati con il sorgente della pipeline che ha generato questo documento "
        f"(amplicon16s {_e(__version__)}) e con l'immagine della configurazione in uso: a "
        "parità di versione una provenienza diversa non rende la fase da rifare, ma va "
        "conosciuta."
    )
    sezione.tabella(_tabella(
        "provenienza", "Provenienza registrata nei manifesti di fase",
        ("fase", "versione", "sorgente", "commit", "modifiche non committate",
         "immagine dichiarata", "conclusa", "avviso di provenienza"),
        righe,
    ))

    sezione.sottotitolo("dada2 e correzione dei pareggi di assignTaxonomy")
    blocchi = [a["blocco_r"] for a in esecuzione.avvii if a["blocco_r"]]
    if blocchi:
        blocco = blocchi[-1]
        correzione = blocco["dada2"].get("correzione") or {}
        sezione.sintesi((
            ("versione di dada2 dichiarata", blocco["dada2"].get("versione") or "non dichiarata"),
            ("versione di partenza", correzione.get("Base", "nessuna correzione dichiarata")),
            ("correzione", correzione.get("File", "")),
            ("SHA-256 della correzione", correzione.get("SHA256", "")),
            ("file di blocco", f"{blocco['file']} (SHA-256 {blocco['sha256']})"),
        ))
        sezione.testo(
            "È quanto dichiara il file di blocco dei pacchetti R letto all'ultimo avvio che "
            "lo ha registrato, come l'immagine è quella dichiarata in configurazione: la "
            "corrispondenza fra il file di blocco e le librerie installate è verificata alla "
            "costruzione dell'immagine, non misurata dalle fasi.", "nota",
        )
    else:
        sezione.testo(
            "Nessun avvio di questa esecuzione ha registrato il file di blocco dei pacchetti "
            f"R (<code>run.lockfile</code>: {_e(esecuzione.parametri.get('run.lockfile', ''))}): "
            "il registro degli avvii è assente, oppure il file non era "
            "raggiungibile all'avvio. La versione di dada2 non è determinabile dal contenuto "
            "della cartella.", "nota",
        )
    return sezione, avvisi


def _decisioni_controlli(esecuzione: _Esecuzione, sezione: _Sezione) -> None:
    """Soglia di profondità per piastra (S11) e modalità di decontaminazione (S12)."""
    soglia = esecuzione.json(Passo.S11, "soglia.json")
    if soglia is not None:
        sezione.sottotitolo("Soglia di profondità per piastra (S11)")
        sezione.sintesi((
            ("modo dichiarato", soglia.get("modo", "")),
            ("modello scelto", soglia.get("scelta", "")),
            ("motivo della scelta", soglia.get("motivo_scelta", "")),
        ))
        piastre = dict(soglia.get("per_piastra", {}))
        if "senza_piastra" in soglia:
            piastre["senza piastra"] = soglia["senza_piastra"]
        sezione.tabella(_tabella(
            "soglie_profondita", "Soglia di profondità e ripieghi (11_controls/soglia.json)",
            ("piastra", "soglia (letture)", "stadio delle letture", "origine", "motivo del ripiego"),
            ((p, s["valore"], s.get("stadio", ""), s.get("origine", ""), s.get("motivo", ""))
             for p, s in sorted(piastre.items(), key=lambda voce: _ordine_piastre(voce[0]))),
        ))
    positivi = esecuzione.json(Passo.S11, "riepilogo.json")
    if positivi is not None:
        sezione.testo(
            f"Controlli positivi conformi al proprio livello di diluizione: "
            f"{_e(positivi.get('conformi'))} su {_e(positivi.get('positivi'))}; non conformi: "
            f"{_e(', '.join(positivi.get('non_conformi_elenco', [])) or 'nessuno')}."
        )

    decontam = esecuzione.json(Passo.S12, "decontam_riepilogo.json")
    if decontam is not None:
        sezione.sottotitolo("Modalità di decontaminazione (S12)")
        sezione.sintesi((
            ("modalità dichiarata", decontam.get("modalita_dichiarata", "")),
            ("esito", decontam.get("esito", "")),
            ("metodo e soglia", f"{decontam.get('metodo', '')}, {decontam.get('soglia', '')}"),
            ("contaminanti rimossi", _n(decontam.get("contaminanti_rimossi", 0))),
        ))
        sezione.tabella(_tabella(
            "decontaminazione",
            "Le due modalità a confronto: una decide, l'altra è diagnostica "
            "(11_controls/decontam_riepilogo.json)",
            ("modalità", "contaminanti", *(f"letture rimosse: {c}" for c in _CLASSI)),
            # Una modalita' non calcolata (senza negativi a sufficienza) e' nulla;
            # una classe senza campioni non ha la sua quota.
            ((nome, m["contaminanti"], *(_pct(m["letture_rimosse"].get(c)) for c in _CLASSI))
             for nome, m in decontam.get("modalita", {}).items() if m is not None),
        ))


def _decisioni_campioni(esecuzione: _Esecuzione, sezione: _Sezione) -> None:
    """Campioni riclassificati dalla configurazione (S0) ed esclusioni dei filtri
    finali (S13).
    """
    inventario = esecuzione.tsv(Passo.S0, NOME_CROSSWALK)
    if inventario is not None:
        dichiarati = set(esecuzione.parametri.get("ctrl.blank_values", []))
        # Il crosswalk conserva il materiale dichiarato accanto alla classe
        # assegnata: un controllo negativo con un altro materiale è stato
        # riclassificato dalla regola ctrl.blank_override_*.
        riclassificati = [
            r for r in inventario
            if r["classe"] == ClasseCampione.CONTROLLO_NEGATIVO.value and r["materiale"] not in dichiarati
        ]
        sezione.sottotitolo("Campioni riclassificati dalla regola di configurazione (S0)")
        sezione.testo(
            f"{_n(len(riclassificati))} campioni sono trattati come controlli negativi pur "
            "dichiarando un altro materiale, per la regola "
            f"<code>ctrl.blank_override_column</code> = "
            f"{_e(_valore(esecuzione.parametri.get('ctrl.blank_override_column')))}, "
            f"<code>ctrl.blank_override_values</code> = "
            f"{_e(_valore(esecuzione.parametri.get('ctrl.blank_override_values')))}. Il "
            "materiale dichiarato resta nell'inventario."
        )
        if riclassificati:
            sezione.tabella(_tabella(
                "riclassificati", "Campioni riclassificati (01_input_validation/crosswalk.tsv)",
                ("accession", "campione", "materiale dichiarato", "posizione", "piastra", "classe assegnata"),
                ((r["accession"], r["sample_name"], r["materiale"], r["posizione"], r["piastra"],
                  r["classe"]) for r in riclassificati),
                ripiegata=True,
            ))

    filtri = esecuzione.json(Passo.S13, "filtri_riepilogo.json")
    esclusi = esecuzione.tsv(Passo.S13, "esclusioni.tsv")
    rimosse = esecuzione.tsv(Passo.S13, "varianti_rimosse.tsv")
    if filtri is None:
        return
    sezione.sottotitolo("Campioni e varianti esclusi dai filtri finali (S13)")
    ordine = filtri["ordine"]
    sezione.tabella(_tabella(
        "filtri_finali", "Esclusioni per filtro, nell'ordine di applicazione (12_final/intermedi/filtri_riepilogo.json)",
        ("filtro", "campioni biologici esclusi", "varianti rimosse"),
        ((f, filtri["campioni"]["esclusi"].get(f, 0),
          filtri["varianti"]["rimosse"].get(f, "non si applica"))
         for f in ordine),
    ))
    prevalenza = filtri.get("prevalenza", {})
    sezione.testo(
        f"Campioni biologici: {_n(filtri['campioni']['biologici'])} in ingresso, "
        f"{_n(filtri['campioni']['finali'])} nell'oggetto finale; i "
        f"{_n(filtri['campioni']['controlli_a_parte'])} controlli sono conservati a parte. "
        f"Varianti: {_n(filtri['varianti']['iniziali'])} in ingresso, "
        f"{_n(filtri['varianti']['finali'])} finali. Prevalenza "
        f"{'applicata' if prevalenza.get('applicata') else 'non applicata'}: denominatore "
        f"{_e(prevalenza.get('denominatore'))} campioni tenuti, almeno "
        f"{_e(prevalenza.get('min_campioni'))} campioni per variante."
    )
    if esclusi:
        sezione.tabella(_tabella(
            "campioni_esclusi", "Campioni esclusi, con il motivo (12_final/intermedi/esclusioni.tsv)",
            ("accession", "campione", "piastra", "filtro", "motivo"),
            ((r["accession"], r["sample_name"], r["piastra"], r["filtro"], r["motivo"]) for r in esclusi),
            ripiegata=True,
        ))
    if rimosse:
        sezione.tabella(_tabella(
            "varianti_rimosse", "Varianti rimosse, con il motivo (12_final/intermedi/varianti_rimosse.tsv)",
            ("variante", "filtro", "motivo", "letture nei biologici tenuti"),
            ((r["asv_id"], r["filtro"], r["motivo"], r["letture_biologici_tenuti"]) for r in rimosse),
            ripiegata=True, limite=50,
        ))


def troncamento_suggerito(
    lunghezze: Iterable[Mapping[str, str]],
    qualita: Iterable[Mapping[str, str]],
    classi: Mapping[str, str],
) -> dict[str, int | None]:
    """Il troncamento che le letture suggeriscono: un'indicazione, non una scelta.

    E' il minore fra due valori. ``per_lunghezza``: il troncamento piu' lungo
    che scarta, perche' piu' corte, non oltre ``QUOTA_CORTE_SUGGERITA`` delle
    letture dei campioni biologici e dei controlli positivi (i negativi non
    contano: le loro poche letture sono spesso dimeri corti). ``per_qualita``:
    l'ultima posizione prima che la qualita' mediana dei campioni biologici (la
    mediana, fra i campioni, della mediana di ciascuno) scenda sotto
    ``QUALITA_MINIMA_SUGGERITA``; la lunghezza massima se non scende mai.
    ``lunghezze`` e ``qualita`` sono le righe delle tabelle di S1, ``classi``
    la classe di ogni campione. Un valore e' ``None`` se mancano le letture
    per calcolarlo.
    """
    controllate = {c.value for c in CLASSI_CONTROLLATE}
    per_lunghezza: dict[int, int] = defaultdict(int)
    for riga in lunghezze:
        if classi.get(riga["campione"]) in controllate:
            per_lunghezza[int(riga["lunghezza"])] += int(riga["letture"])
    totale = sum(per_lunghezza.values())
    da_lunghezza = None
    piu_corte = 0
    for lunghezza in sorted(per_lunghezza):
        # Troncare a questa lunghezza scarta le letture piu' corte di essa.
        if piu_corte > QUOTA_CORTE_SUGGERITA * totale:
            break
        da_lunghezza = lunghezza
        piu_corte += per_lunghezza[lunghezza]

    mediane: dict[int, list[float]] = defaultdict(list)
    for riga in qualita:
        if classi.get(riga["campione"]) == ClasseCampione.BIOLOGICO.value:
            mediane[int(riga["posizione"])].append(float(riga["mediana"]))
    da_qualita = None
    for posizione in sorted(mediane):
        if statistics.median(mediane[posizione]) < QUALITA_MINIMA_SUGGERITA:
            break
        da_qualita = posizione

    presenti = [v for v in (da_lunghezza, da_qualita) if v is not None]
    return {
        "per_lunghezza": da_lunghezza,
        "per_qualita": da_qualita,
        "suggerito": min(presenti) if len(presenti) == 2 else None,
    }


def _decisioni_troncamento(esecuzione: _Esecuzione, sezione: _Sezione) -> None:
    """Il troncamento suggerito dalle letture (S1), accanto a quello dichiarato."""
    lunghezze = esecuzione.tsv(Passo.S1, NOME_LUNGHEZZE)
    qualita = esecuzione.tsv(Passo.S1, NOME_QUALITA)
    inventario = esecuzione.tsv(Passo.S0, NOME_CROSSWALK)
    if lunghezze is None or qualita is None or inventario is None:
        return
    esito = troncamento_suggerito(
        lunghezze, qualita, {r["accession"]: r["classe"] for r in inventario}
    )
    sezione.sottotitolo("Troncamento suggerito dalle letture (S1)")
    sezione.sintesi((
        ("troncamento dichiarato (filter.truncLen)",
         _valore(esecuzione.parametri.get("filter.truncLen"))),
        ("troncamento suggerito", _e(esito["suggerito"] if esito["suggerito"] is not None
                                     else "non calcolabile")),
        (f"per lunghezza: il più lungo che scarta non oltre "
         f"{_pct(QUOTA_CORTE_SUGGERITA)} delle letture di biologici e positivi",
         _e(esito["per_lunghezza"] if esito["per_lunghezza"] is not None else "non calcolabile")),
        (f"per qualità: l'ultima posizione con qualità mediana dei biologici almeno "
         f"{QUALITA_MINIMA_SUGGERITA}",
         _e(esito["per_qualita"] if esito["per_qualita"] is not None else "non calcolabile")),
    ))
    sezione.testo(
        "Il valore suggerito è un'indicazione ricavata da 02_qc_profiles/lunghezze.tsv e "
        "qualita.tsv: il troncamento applicato resta quello dichiarato in "
        "<code>filter.truncLen</code>, che è una scelta di chi conduce l'analisi."
    )


def _sezione_decisioni(esecuzione: _Esecuzione) -> tuple[_Sezione, int, int]:
    """Le decisioni prese senza intervento umano, lette dagli artefatti che le
    registrano; con il numero di tentativi ripetuti e di degradazioni.
    """
    sezione = _Sezione("decisioni", "Decisioni prese automaticamente")
    sezione.testo(
        "Ogni decisione è letta dall'artefatto che la registra, con il motivo scritto dalla "
        "fase che l'ha presa: nulla in questa sezione è ricalcolato."
    )
    tentativi = _aggiustamenti(esecuzione)
    sezione.sottotitolo("Tentativi ripetuti")
    if tentativi:
        sezione.tabella(_tabella(
            "tentativi_ripetuti", "Tentativi ripetuti, con l'aggiustamento applicato (manifesti di fase)",
            ("fase", "codice", "tentativo fallito", "parametro", "valore dichiarato",
             "valore precedente", "valore usato", "azione"),
            ((p, a["codice"], a["tentativo_fallito"], a["parametro"] or "nessuno: configurazione invariata",
              _valore(a["dichiarato"]) if a["parametro"] else "",
              _valore(a["precedente"]) if a["parametro"] else "",
              _valore(a["usato"]) if a["parametro"] else "", a["azione"])
             for p, a in tentativi),
        ))
    else:
        sezione.testo("Nessuna fase conclusa ha richiesto un tentativo ripetuto.")

    degradazioni = [(p, d) for p, m in esecuzione.manifesti.items() for d in m.degradazioni]
    sezione.sottotitolo("Degradazioni")
    if degradazioni:
        sezione.tabella(_tabella(
            "degradazioni", "Degradazioni: la fase ha proseguito con un ripiego (manifesti di fase)",
            ("fase", "codice", "sintesi", "dettaglio", "ripiego adottato"),
            ((p, d["codice"], d.get("sintesi", ""), d.get("dettaglio", ""), d.get("azione", ""))
             for p, d in degradazioni),
        ))
    else:
        sezione.testo("Nessuna fase conclusa ha registrato una degradazione.")

    _decisioni_troncamento(esecuzione, sezione)
    _decisioni_controlli(esecuzione, sezione)
    _decisioni_campioni(esecuzione, sezione)
    return sezione, len(tentativi), len(degradazioni)


def _sezione_tracciamento(esecuzione: _Esecuzione) -> tuple[_Sezione, Tracciamento, dict[str, str]]:
    """Le letture di ogni campione fase per fase, e il riassunto per classe."""
    sezione = _Sezione("tracciamento", "Tracciamento delle letture")
    tracciamento = esecuzione.tracciamento()
    inventario = esecuzione.tsv(Passo.S0, NOME_CROSSWALK) or []
    classe = {r["accession"]: r["classe"] for r in inventario}
    if not tracciamento.passi:
        sezione.testo("Nessuna fase conclusa ha registrato un passo del tracciamento.")
        return sezione, tracciamento, classe
    sezione.testo(
        "Ogni fase registra quante letture restano a ciascun campione dopo il proprio passo; "
        "la tabella è ricomposta dai file delle sole fasi concluse, nell'ordine del grafo. "
        "I controlli negativi sono riassunti a parte: con poca o nessuna biomassa perdono "
        "letture per ragioni diverse dai campioni biologici, e una media sull'insieme le "
        "confonderebbe."
    )
    righe = []
    for passo in tracciamento.passi:
        riga: list[Any] = [passo, tracciamento.origine[passo]]
        for c in _CLASSI:
            valori = [v[passo] for a, v in tracciamento.letture.items() if classe.get(a) == c and passo in v]
            mediana = statistics.median(valori) if valori else ""
            # Con un numero pari di campioni la mediana puo' cadere a meta' fra due conteggi.
            riga += [len(valori), sum(valori), int(mediana) if mediana == int(mediana or 0) else mediana]
        righe.append(riga)
    sezione.tabella(_tabella(
        "tracciamento_per_classe", "Letture per classe di campioni, passo per passo",
        ("passo", "fase", *(f"{c}: {m}" for c in _CLASSI for m in ("campioni", "letture", "mediana"))),
        righe,
    ))
    nome = {r["accession"]: r for r in inventario}
    sezione.tabella(_tabella(
        "tracciamento_per_campione", "Letture per campione, passo per passo",
        ("accession", "campione", "classe", "piastra", *tracciamento.passi),
        ((a, nome.get(a, {}).get("sample_name", ""), classe.get(a, ""), nome.get(a, {}).get("piastra", ""),
          *(v.get(p, "") for p in tracciamento.passi))
         for a, v in sorted(tracciamento.letture.items())),
        ripiegata=True,
    ))
    return sezione, tracciamento, classe


def _sezione_risultato(
    esecuzione: _Esecuzione, tracciamento: Tracciamento, classe: Mapping[str, str]
) -> tuple[_Sezione, str]:
    """L'oggetto finale: dimensioni, letture per classe, file consegnati."""
    sezione = _Sezione("risultato", "Risultato finale")
    finale = esecuzione.manifesti.get(Passo.S14)
    dimensioni = "non prodotto"
    if finale is None:
        sezione.testo("S14 non è conclusa: l'oggetto finale non è stato prodotto.")
    else:
        metriche = finale.metriche
        dimensioni = (
            f"{_n(metriche['campioni'])} campioni x {_n(metriche['varianti'])} varianti"
            if {"campioni", "varianti"} <= set(metriche)
            else "dimensioni non registrate nel manifesto di S14"
        )
        filtri = esecuzione.json(Passo.S13, "filtri_riepilogo.json") or {}
        sezione.sintesi((
            ("oggetto finale", f"ps_final.rds, {dimensioni}"),
            ("letture nell'oggetto finale", _n(filtri["letture"]["finali"]) if filtri else ""),
            ("frazione delle letture trattenute dai campioni finali, da senza chimere a finali",
             _pct(metriche["frazione_letture_trattenute"])
             if "frazione_letture_trattenute" in metriche else "non registrata"),
            ("soglia qc.min_frac_reads_retained", _valore(esecuzione.parametri.get("qc.min_frac_reads_retained"))),
            ("conclusa", finale.conclusa),
        ))
        sezione.testo(
            "L'oggetto finale contiene i soli campioni biologici che hanno superato i filtri; "
            + ("i controlli sono consegnati a parte in <code>ps_controlli.rds</code>."
               if metriche.get("controlli", True) else
               "l'inventario non ha controlli, e <code>ps_controlli.rds</code> non esiste.")
        )
        checksum = esecuzione.byte(Passo.S14, "checksum.sha256")
        if checksum is not None:
            sezione.tabella(_tabella(
                "file_finali", "File consegnati e loro impronta (12_final/checksum.sha256)",
                ("file", "SHA-256"),
                (reversed(riga.split(maxsplit=1)) for riga in checksum.decode("utf-8").splitlines() if riga),
            ))
    if tracciamento.passi:
        sezione.tabella(_tabella(
            "letture_per_classe", "Letture totali per classe lungo la catena",
            ("passo", "fase", *_CLASSI, "totale"),
            ((p, tracciamento.origine[p],
              *(sum(v.get(p, 0) for a, v in tracciamento.letture.items() if classe.get(a) == c)
                for c in _CLASSI),
              tracciamento.totali()[p])
             for p in tracciamento.passi),
        ))
    modello = esecuzione.manifesti.get(Passo.S3)
    for nome in (n for n in (modello.nomi if modello else ()) if n.endswith(".png")):
        immagine = esecuzione.byte(Passo.S3, nome)
        if immagine is not None:
            sezione.blocchi.append(
                f'<figure><img alt="{_e(nome)}" src="data:image/png;base64,'
                f'{base64.b64encode(immagine).decode("ascii")}">'
                f"<figcaption>Modello d'errore stimato da S3: {_e(nome)} "
                f"(04_error_models)</figcaption></figure>"
            )
    return sezione, dimensioni


def _sezione_evidenza(
    esecuzione: _Esecuzione, avvisi: list[tuple[str, str]], superate: list[str],
    tentativi: int, degradazioni: int,
) -> tuple[_Sezione, tuple[str, ...]]:
    """Ciò che chi apre il report deve vedere per primo."""
    sezione = _Sezione("evidenza", "Segnalazioni in apertura")
    segnalazioni: list[str] = []

    tabella = _valori_osd734(esecuzione)
    ereditati = len(tabella.righe)
    sezione.sottotitolo("Parametri tarati su OSD-734 e non dichiarati")
    if ereditati:
        segnalazioni.append(
            f"{ereditati} parametri non dichiarati valgono il predefinito tarato su OSD-734"
        )
        sezione.testo(
            f"<strong>{ereditati} parametri</strong> non sono stati dichiarati nella "
            "configurazione e valgono quindi il predefinito, che è stato tarato sul dataset di "
            "riferimento (NASA GeneLab OSD-734). La tabella li elenca con il fatto, accertato "
            "su OSD-734, che giustificava ciascun valore. <strong>Per ognuno va controllato se "
            "quel fatto vale anche per il dataset in uso</strong>: se vale, il parametro va "
            "dichiarato nella configurazione con lo stesso valore, e non comparirà più qui; se "
            "non vale, va dichiarato con il valore adatto al dataset. I parametri che descrivono "
            "il dataset (formato dei metadati, etichette, primer, troncamento) non compaiono "
            "perché sono obbligatori: la configurazione li dichiara tutti.", "evidenza",
        )
    elif not esecuzione.versioni:
        segnalazioni.append("nessuna configurazione registrata")
        sezione.testo(
            "Nessuna configurazione è registrata in <code>00_config</code>: i valori dei "
            "parametri non sono verificabili.", "evidenza",
        )
    else:
        sezione.testo(
            "Nessun parametro vale un predefinito tarato su OSD-734 senza essere stato "
            "dichiarato nella configurazione.", "regolare",
        )
    if tabella.righe:
        sezione.tabella(tabella)

    sezione.sottotitolo("Avvisi di provenienza")
    if avvisi:
        segnalazioni.append(f"avvisi di provenienza su {', '.join(p for p, _ in avvisi)}")
        sezione.testo(
            f"<strong>{len(avvisi)} fasi</strong> sono state calcolate da un sorgente o con "
            "un'immagine diversi da quelli della pipeline che ha generato questo documento, a "
            "parità di versione del calcolo. I risultati non sono invalidati: la differenza "
            'va conosciuta e, se non è attesa, chiarita (sezione <a href="#provenienza">'
            "Provenienza delle fasi</a>).", "evidenza",
        )
        sezione.tabella(_tabella(
            "avvisi_provenienza", "Avvisi di provenienza",
            ("fase", "avviso"), avvisi,
        ))
    else:
        sezione.testo(
            "Nessun avviso: il sorgente e l'immagine registrati dalle fasi concluse "
            "coincidono con quelli della pipeline che ha generato questo documento.", "regolare",
        )

    sezione.sottotitolo("Integrità di ciò che è stato letto")
    problemi = (
        [f"manifesto non valido: {p} ({m})" for p, m in esecuzione.non_validi.items()]
        + [f"artefatto non integro: {a}" for a in esecuzione.alterati]
        + [f"fase calcolata su artefatti superati: {s}" for s in superate]
    )
    if problemi:
        segnalazioni += problemi
        voci = "".join(f"<li>{_e(p)}</li>" for p in problemi)
        sezione.blocchi.append(
            f'<div class="evidenza"><p>Il contenuto della cartella non corrisponde in tutto a '
            f"quanto registrato nei manifesti; ciò che non è integro non è riportato nel "
            f"documento.</p><ul>{voci}</ul></div>"
        )
    else:
        sezione.testo(
            "I manifesti di fase sono validi, ogni fase conclusa è calcolata sugli artefatti a "
            "monte oggi su disco, e ogni artefatto letto per questo documento corrisponde al "
            "checksum del suo manifesto.", "regolare",
        )

    sezione.sottotitolo("Decisioni automatiche")
    if tentativi or degradazioni:
        segnalazioni.append(f"tentativi ripetuti: {tentativi}; degradazioni: {degradazioni}")
    sezione.testo(
        f"Tentativi ripetuti: {tentativi}. Degradazioni: {degradazioni}. Il dettaglio, con "
        'i motivi, è nella sezione <a href="#decisioni">Decisioni prese automaticamente</a>.',
        "evidenza" if tentativi or degradazioni else "regolare",
    )
    return sezione, tuple(segnalazioni)


# --------------------------------------------------------------------------- #
# Costruzione e scrittura                                                      #
# --------------------------------------------------------------------------- #


def costruisci(
    out_root: Path | str, passi: Mapping[Passo, PipelineStep] | None = None
) -> Report:
    """Costruisce il report dalla cartella di output, senza scrivere nulla.

    ``passi`` sono le fasi realizzate, per difetto quelle della pipeline: ne
    servono la cartella, i passi del tracciamento e il sorgente.
    """
    esecuzione = _Esecuzione(Path(out_root), passi_realizzati() if passi is None else passi)
    stato, completa, superate = _sezione_stato(esecuzione)
    gate = _sezione_gate(esecuzione)
    configurazione = _sezione_configurazione(esecuzione)
    provenienza, avvisi = _sezione_provenienza(esecuzione)
    decisioni, tentativi, degradazioni = _sezione_decisioni(esecuzione)
    tracciamento, letture, classe = _sezione_tracciamento(esecuzione)
    risultato, dimensioni = _sezione_risultato(esecuzione, letture, classe)
    # Per ultima, perché riassume le altre; nel documento apre.
    evidenza, segnalazioni = _sezione_evidenza(esecuzione, avvisi, superate, tentativi, degradazioni)

    sezioni = (evidenza, stato, gate, configurazione, provenienza, decisioni, tracciamento, risultato)
    inventario = esecuzione.json(Passo.S0, "inventario.json") or {}
    versione = esecuzione.versioni[-1] if esecuzione.versioni else {}
    concluse = [m.conclusa for m in esecuzione.manifesti.values()]
    intestazione = _Sezione("", "")
    intestazione.sintesi((
        ("catena", "completa" if completa else "non completa: vedi lo stato delle fasi"),
        ("fasi concluse", f"{len(esecuzione.manifesti)} su {len(GRAFO)}"),
        ("ultima fase conclusa", max(concluse) if concluse else "nessuna"),
        ("campioni", ", ".join(f"{_n(n)} {c}" for c, n in inventario.get("conteggi", {}).items())
         or "inventario non disponibile"),
        ("oggetto finale", dimensioni),
        ("configurazione in uso", f"{versione.get('file', 'nessuna')} {versione.get('digest', '')}".strip()),
        ("versioni della configurazione", len(esecuzione.versioni)),
        ("avvii registrati", len(esecuzione.avvii) if esecuzione.registro_avvii
         else "registro degli avvii assente"),
    ))
    documento = Template((_MODELLI / "report.html").read_text(encoding="utf-8")).substitute(
        titolo="Report di esecuzione della pipeline 16S rRNA single-end",
        stile=(_MODELLI / "stile.css").read_text(encoding="utf-8"),
        intestazione=intestazione.blocchi[0],
        indice="<ol>" + "".join(
            f'<li><a href="#{s.ident}">{_e(s.titolo)}</a></li>' for s in sezioni
        ) + "</ol>",
        corpo="\n".join(s.html() for s in sezioni),
        piede=(
            "<p>Documento ricavato dal solo contenuto della cartella di output: manifesti di "
            "fase, artefatti, configurazioni registrate e registro degli avvii. Non riporta la data "
            "in cui è stato generato: le date sono quelle dei manifesti e del registro. Generarlo "
            "di nuovo dalla stessa cartella dà lo stesso documento.</p>"
        ),
    )
    return Report(documento, tuple(t for s in sezioni for t in s.tabelle), segnalazioni)


def scrivi(report: Report, out_root: Path | str) -> Path:
    """Scrive il report nella cartella ``report`` e ne restituisce il percorso.

    Ogni file passa per una scrittura atomica. I file di una generazione
    precedente che questa non produce più vengono rimossi: la cartella contiene
    sempre e solo il report dello stato letto.
    """
    cartella = Path(out_root) / CARTELLA_REPORT
    tabelle = cartella / CARTELLA_TABELLE
    tabelle.mkdir(parents=True, exist_ok=True)
    attesi = {f"{t.nome}.tsv" for t in report.tabelle}
    if len(attesi) != len(report.tabelle):
        raise ValueError("due tabelle del report hanno lo stesso nome")
    for tabella in report.tabelle:
        scrivi_atomico(tabelle / f"{tabella.nome}.tsv", tabella.tsv())
    for superato in tabelle.iterdir():
        if superato.name not in attesi:
            superato.unlink()
    scrivi_atomico(cartella / NOME_REPORT, report.html)
    return cartella / NOME_REPORT


def genera(out_root: Path | str, passi: Mapping[Passo, PipelineStep] | None = None) -> Path:
    """Costruisce il report dalla cartella di output e lo scrive."""
    return scrivi(costruisci(out_root, passi), out_root)
