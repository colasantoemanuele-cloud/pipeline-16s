"""Suite di verifica del grafo DAG delle fasi, della propagazione dell'invalidità e della ripresa su filesystem.

Inquadramento nel Piano Operativo
---------------------------------
* **Settimane di riferimento**: **Settimana 8 e Settimana 9 (W8/W9 — Fase F3:
  Grafo DAG delle dipendenze, propagazione dell'invalidità e ripresa/resume su
  filesystem)**.
* **Moduli sorgente coperti**:
  - ``src/amplicon16s/steps/base.py``
  - ``src/amplicon16s/runner/graph.py``
  - ``src/amplicon16s/runner/project.py``
  - ``src/amplicon16s/steps/s00_validate.py``

Scopo sperimentale e razionale scientifico/sistemistico
-------------------------------------------------------
La Fase S0 utilizzata nei test è quella reale, eseguita su uno scenario
sintetico scritto su disco; le fasi da ``S1`` a ``S14`` sono rappresentate da
doppioni deterministici (``Doppione``) che scrivono su disco un artefatto
dipendente dalle impronte a monte e dai propri parametri, annotando ogni
effettiva esecuzione. In questo modo i test misurano **quale lavoro la ripresa
rifà davvero su disco**, certificando sette proprietà architetturali:

1. **Ordine topologico guidato dalle dipendenze di dato (DAG)**: il grafo esegue
   e invalida le fasi seguendo l'effettivo flusso degli artefatti e non una
   sequenza monolitica rigida (es. ricalcolare ``S1`` non invalida ``S2``;
   modificare ``tax.min_boot`` in ``S8`` non ricalcola l'albero filogenetico ``S9``).
2. **Precedenza metodologica obbligatoria ``S12 -> S13`` (``E-S13-01``)**:
   garantisce che la decontaminazione ``decontam`` (``S12``) preceda sempre il
   filtro di prevalenza (``S13``), impedendo che contaminanti da reagente
   superino la soglia di prevalenza prima di essere rimossi.
3. **Filogenesi ``S9`` come unica fase facoltativa (``phylo.enabled``)**: se
   disattivata, l'assemblaggio ``phyloseq`` in ``S10`` adatta dinamicamente le
   proprie dipendenze (``S0, S7, S8``) e la corsa risulta completa senza ``S9``.
4. **Isolamento dei manifesti nelle cartelle condivise (``07_chimera``,
   ``11_controls``, ``12_final``)**: ogni passo scrive il proprio manifesto
   atomico ``manifest_<passo>.json`` (es. ``manifest_s11.json`` e
   ``manifest_s12.json``), impedendo che il completamento di ``S11`` faccia
   apparire già conclusa ``S12``.
5. **Propagazione selettiva dell'invalidità**: la cancellazione o alterazione di
   un artefatto (es. in ``S4``) invalida ``S4`` e le sole fasi discendenti,
   preservando intatto il lavoro già svolto in ``S0..S3``.
6. **Invarianza rispetto ai parametri operativi (``run.threads``, ``io.out_root``,
   ``retry.*``)**: variare il numero di CPU o spostare la cartella di output non
   invalida ore di calcolo DADA2/DECIPHER già concluse, poiché non altera il
   risultato biologico.
7. **Caching dei checksum per valutazione**: durante una valutazione dello stato
   (``run.valuta()``), l'hash SHA-256 di ciascun artefatto viene calcolato una
   sola volta su disco per non degradare le prestazioni.
"""

from __future__ import annotations

import json
import os
import shutil
from typing import Any, ClassVar

import pytest
from conftest import NEGATIVO, POSITIVO, Campione, crea_scenario

from amplicon16s.config.resolve import risolvi
from amplicon16s.config.schema import Config, valida
from amplicon16s.errors.exceptions import ErroreRevisioneUmana
from amplicon16s.gates.g01_g15 import Contesto
from amplicon16s.io_layer.artifacts import (
    AlberoOutput,
    Fase,
    nome_manifesto_passo,
)
from amplicon16s.runner.graph import GRAFO, Grafo, Nodo, Passo
from amplicon16s.runner.project import (
    ProjectRun,
    StatoPasso,
    passi_realizzati,
)
from amplicon16s.steps.base import (
    Esito,
    PipelineStep,
    Produzione,
    StepContext,
    impronta_parametri,
)
from amplicon16s.steps.s00_validate import ValidazioneIngressi, esegui_s0

TUTTE = tuple(Passo)


# --------------------------------------------------------------------------- #
# Scenario e doppioni                                                          #
# --------------------------------------------------------------------------- #


def _campioni() -> list[Campione]:
    return [
        Campione("ERX3000001", "NOD1D4.L1", piastra="1"),
        Campione("ERX3000002", "NOD1D4.L2", piastra="1"),
        Campione("ERX3000003", "POS.P1.1", materiale=POSITIVO,
                 posizione="Not Applicable", piastra="1"),
        Campione("ERX3000004", "BLANK.P1.1", materiale=NEGATIVO,
                 posizione="Not Applicable", piastra="1"),
    ]


def _variante(config: Config, **valori: Any) -> Config:
    """La stessa configurazione con alcuni parametri cambiati: ``gruppo__chiave``."""
    dati = config.model_dump(mode="python")
    for nome, valore in valori.items():
        gruppo, chiave = nome.split("__")
        dati[gruppo][chiave] = valore
    return valida(dati)


class Doppione(PipelineStep):
    """Una fase finta che scrive un artefatto derivato dai propri ingressi."""

    def __init__(self, registro: list[Passo]) -> None:
        self.registro = registro

    @property
    def nome_artefatto(self) -> str:
        return f"{str(self.passo).lower()}.json"

    def calcola(self, contesto: StepContext) -> Produzione:
        self.registro.append(self.passo)
        contenuto = {
            "passo": str(self.passo),
            "a_monte": {str(p): i for p, i in contesto.a_monte.items()},
            "configurazione": impronta_parametri(contesto.risolta, self.parametri),
        }
        artefatto = contesto.albero.scrivi_testo(
            self.cartella, self.nome_artefatto, json.dumps(contenuto, sort_keys=True)
        )
        return Produzione((artefatto,), {"righe": 1})


class S0Contata(ValidazioneIngressi):
    """La S0 vera, che annota le proprie esecuzioni."""

    def __init__(self, registro: list[Passo]) -> None:
        super().__init__()
        self.registro = registro

    def calcola(self, contesto: StepContext) -> Produzione:
        self.registro.append(self.passo)
        return super().calcola(contesto)


def _passi(
    registro: list[Passo], parametri: dict[Passo, tuple[str, ...] | None] | None = None
) -> dict[Passo, PipelineStep]:
    """S0 vera piu' un doppione per ogni altra fase, con i parametri indicati."""
    parametri = parametri or {}
    classe_s0 = type(
        "S0Parametrica", (S0Contata,), {"parametri": parametri.get(Passo.S0)}
    )
    passi: dict[Passo, PipelineStep] = {Passo.S0: classe_s0(registro)}
    for passo in TUTTE[1:]:
        classe = type(
            f"Doppione{passo}",
            (Doppione,),
            {"passo": passo, "parametri": parametri.get(passo)},
        )
        passi[passo] = classe(registro)
    return passi


#: Parametri ristretti: ogni fase dipende solo dal proprio gruppo, cosi' si
#: vede quali fasi invalida un singolo cambiamento.
RISTRETTI: dict[Passo, tuple[str, ...]] = {
    Passo.S0: ("io", "meta", "ctrl", "qc"),
    Passo.S1: ("qc",),
    Passo.S2: ("filter",),
    Passo.S3: ("err",),
    Passo.S4: ("asv",),
    Passo.S5: ("asv",),
    Passo.S6: ("chimera",),
    Passo.S7: ("asv",),
    Passo.S8: ("tax",),
    Passo.S9: ("phylo",),
    Passo.S10: ("out",),
    Passo.S11: ("katharoseq",),
    Passo.S12: ("decontam",),
    Passo.S13: ("prev",),
    Passo.S14: ("out",),
}


def riprendi(run: ProjectRun, registro: list[Passo]) -> list[Passo]:
    """Esegue in ordine le fasi da eseguire e restituisce quelle eseguite."""
    registro.clear()
    for _ in range(2 * len(TUTTE)):
        valutazione = run.valuta()
        passo = valutazione.prossima()
        if passo is None:
            return list(registro)
        run.fase(passo).esegui(run.contesto(passo, valutazione))
    raise AssertionError("la ripresa non converge")


@pytest.fixture
def scenario(tmp_path):
    return crea_scenario(tmp_path, _campioni(), con_letture=True)


@pytest.fixture
def registro() -> list[Passo]:
    return []


@pytest.fixture
def eseguita(scenario, registro):
    """Un'esecuzione completa, con la filogenesi disattivata come per difetto."""
    run = ProjectRun(scenario.config, passi=_passi(registro))
    riprendi(run, registro)
    return run


def _attive_da(passo: Passo, config: Config) -> list[Passo]:
    """La fase data e tutte quelle attive che ne dipendono, in ordine."""
    coinvolte = {passo, *GRAFO.discendenti(passo)}
    return [p for p in GRAFO.attive(config) if p in coinvolte]


def _file(run: ProjectRun, passo: Passo):
    fase = run.fase(passo)
    return run.albero.cartella(fase.cartella) / fase.nome_artefatto


# --------------------------------------------------------------------------- #
# Il grafo                                                                     #
# --------------------------------------------------------------------------- #


def test_il_grafo_ordina_le_quindici_fasi():
    """
    **Obiettivo**: Verificare che ``GRAFO`` contenga tutte le 15 fasi (``S0``–``S14``)
    in ordine topologico coerente con ``TUTTE``.

    **Razionale Scientifico/Sistemistico**: Garantisce che la definizione statica
    del DAG di pipeline copra l'intero ciclo di vita bioinformatico dalla
    validazione degli ingressi (``S0``) alla generazione dei report finali (``S14``).
    """
    assert GRAFO.ordine() == TUTTE
    assert len(GRAFO) == 15


def test_la_decontaminazione_precede_il_filtro_di_prevalenza():
    """
    **Obiettivo**: Verificare che nel grafo ``Passo.S12`` preceda ``Passo.S13``,
    compaia nelle sue dipendenze dirette e sia vincolato dalla precedenza ``E-S13-01``.

    **Razionale Scientifico/Sistemistico**: La decontaminazione statistica con
    ``decontam`` (``S12``) deve tassativamente precedere il filtro di prevalenza
    ecologica (``S13``): se si filtrasse prima per prevalenza, i contaminanti da
    reagente presenti in quasi tutti i pozzetti supererebbero la soglia mentre il
    conteggio totale delle letture e la composizione su cui lavora ``decontam``
    risulterebbero distorti.
    """
    ordine = GRAFO.ordine()
    assert ordine.index(Passo.S12) < ordine.index(Passo.S13)
    assert Passo.S12 in GRAFO.nodo(Passo.S13).dipendenze
    assert (Passo.S12, Passo.S13, "E-S13-01") in GRAFO.precedenze


def test_ogni_dipendenza_precede_la_fase_che_la_usa():
    """
    **Obiettivo**: Verificare che per ogni nodo di ``GRAFO`` tutte le sue
    dipendenze abbiano un indice inferiore nell'ordinamento ``GRAFO.ordine()``.

    **Razionale Scientifico/Sistemistico**: Certifica che il grafo delle
    dipendenze è un DAG aciclico provvisto di ordinamento topologico valido,
    escludendo cicli o riferimenti in avanti impossibili da soddisfare.
    """
    ordine = GRAFO.ordine()
    for nodo in GRAFO:
        for dipendenza in nodo.dipendenze:
            assert ordine.index(dipendenza) < ordine.index(nodo.passo)


def test_la_sola_fase_facoltativa_e_la_filogenesi():
    """
    **Obiettivo**: Verificare che l'unico nodo con ``facoltativa is True`` in
    ``GRAFO`` sia ``Passo.S9``, governato da ``phylo.enabled``.

    **Razionale Scientifico/Sistemistico**: Tutte le fasi di controllo qualità,
    inferenza ASV, rimozione chimere, tassonomia e decontaminazione sono
    obbligatorie per la validità scientifica dello studio; solo l'albero
    filogenetico (``S9``, computazionalmente oneroso con DECIPHER + phangorn) è
    opzionale quando non sono richieste metriche UniFrac.
    """
    assert [n.passo for n in GRAFO if n.facoltativa] == [Passo.S9]
    assert GRAFO.nodo(Passo.S9).parametro_attivazione == "phylo.enabled"


def test_fasi_che_condividono_una_cartella():
    """
    **Obiettivo**: Verificare che le sole cartelle fisiche condivise da due fasi
    distinte siano ``07_chimera`` (``S6, S7``), ``11_controls`` (``S11, S12``) e
    ``12_final`` (``S13, S14``).

    **Razionale Scientifico/Sistemistico**: Poiché 15 fasi scrivono in 12
    cartelle di output, individuare esplicitamente le 3 coppie che condividono la
    cartella giustifica l'adozione di manifesti atomici per-passo
    (``manifest_<passo>.json``) anziché di un unico manifesto di directory.
    """
    per_cartella: dict[Fase, list[Passo]] = {}
    for nodo in GRAFO:
        per_cartella.setdefault(nodo.cartella, []).append(nodo.passo)
    condivise = {f: p for f, p in per_cartella.items() if len(p) > 1}
    assert condivise == {
        Fase.CHIMERA: [Passo.S6, Passo.S7],
        Fase.CONTROLS: [Passo.S11, Passo.S12],
        Fase.FINAL: [Passo.S13, Passo.S14],
    }


def _nodi_base() -> list[Nodo]:
    return [
        Nodo(Passo.S12, "a", Fase.CONTROLS),
        Nodo(Passo.S13, "b", Fase.FINAL, (Passo.S12,)),
    ]


def test_il_grafo_rifiuta_una_precedenza_obbligatoria_non_dichiarata():
    """
    **Obiettivo**: Verificare che istanziare ``Grafo`` senza che ``S13`` dichiari
    ``S12`` tra le sue dipendenze sollevi ``ValueError("S12 deve precedere S13")``.

    **Razionale Scientifico/Sistemistico**: Impedisce che una modifica futura
    alla definizione dei nodi del grafo possa accidentalmente rimuovere il
    vincolo di dipendenza tra decontaminazione e filtro di prevalenza.
    """
    with pytest.raises(ValueError, match="S12 deve precedere S13"):
        Grafo([Nodo(Passo.S12, "a", Fase.CONTROLS), Nodo(Passo.S13, "b", Fase.FINAL)])


def test_il_grafo_rifiuta_una_precedenza_obbligatoria_facoltativa():
    """
    **Obiettivo**: Verificare che ``Grafo`` sollevi ``ValueError`` se un passo
    che compare come prerequisito in una precedenza obbligatoria (``S12``) viene
    reso facoltativo (``attiva_se=...``).

    **Razionale Scientifico/Sistemistico**: Una fase obbligatoria per metodo
    scientifico (come la decontaminazione ``S12`` prima di ``S13``) non può
    essere disattivabile dall'utente, altrimenti ``S13`` girerebbe su dati non
    decontaminati.
    """
    nodi = _nodi_base()
    nodi[0] = Nodo(Passo.S12, "a", Fase.CONTROLS, attiva_se=lambda c: True)
    with pytest.raises(ValueError, match="facoltativa"):
        Grafo(nodi)


def test_il_grafo_rifiuta_una_dipendenza_che_non_precede():
    """
    **Obiettivo**: Verificare che ``Grafo`` sollevi ``ValueError`` se un nodo
    (``S13``) è posizionato prima della propria dipendenza (``S12``) nella lista dei nodi.

    **Razionale Scientifico/Sistemistico**: Preserva la proprietà fondamentale
    dell'ordinamento topologico per cui ogni produttore di artefatti deve essere
    valutato ed eseguito prima dei suoi consumatori.
    """
    with pytest.raises(ValueError, match="non la precede"):
        Grafo(
            [Nodo(Passo.S13, "b", Fase.FINAL, (Passo.S12,)),
             Nodo(Passo.S12, "a", Fase.CONTROLS)],
            precedenze=(),
        )


def test_con_la_filogenesi_disattivata_s10_non_ne_dipende(scenario):
    """
    **Obiettivo**: Verificare che con ``phylo.enabled = False`` le dipendenze
    attive di ``S10`` siano ``(S0, S7, S8)``, mentre con ``phylo.enabled = True``
    includano anche ``S9``.

    **Razionale Scientifico/Sistemistico**: L'oggetto ``phyloseq`` costruito in
    ``S10`` assembla tabella ASV (``S7``), tassonomia (``S8``), metadati (``S0``)
    ed eventualmente l'albero filogenetico (``S9``); quando ``S9`` è disattivata,
    ``S10`` non deve restare bloccata in attesa di un albero che non verrà prodotto.
    """
    config = scenario.config
    assert not config.phylo.enabled
    assert Passo.S9 not in GRAFO.attive(config)
    assert GRAFO.dipendenze_attive(Passo.S10, config) == (Passo.S0, Passo.S7, Passo.S8)

    attiva = _variante(config, phylo__enabled=True)
    assert Passo.S9 in GRAFO.dipendenze_attive(Passo.S10, attiva)


# --------------------------------------------------------------------------- #
# Prima esecuzione e ripresa di un'esecuzione completa                         #
# --------------------------------------------------------------------------- #


def test_la_prima_esecuzione_esegue_ogni_fase_attiva_una_volta(scenario, registro):
    """
    **Obiettivo**: Verificare che partendo da una cartella di output vuota
    ``riprendi`` esegua esattamente una volta ciascuna delle 14 fasi attive
    (tutte tranne ``S9``) in ordine topologico, portando ``run.completa`` a ``True``.

    **Razionale Scientifico/Sistemistico**: Certifica la convergenza del ciclo
    di esecuzione del DAG dall'inizio alla fine senza esecuzioni ridondanti.
    """
    run = ProjectRun(scenario.config, passi=_passi(registro))
    assert run.prossima() is Passo.S0

    eseguite = riprendi(run, registro)

    assert eseguite == [p for p in TUTTE if p is not Passo.S9]
    assert run.completa
    assert run.prossima() is None


def test_un_esecuzione_completa_non_ripete_alcun_lavoro(eseguita, scenario, registro):
    """
    **Obiettivo**: Verificare che invocare ``riprendi`` su una corsa già
    completa (anche istanziando un nuovo ``ProjectRun`` da zero) restituisca
    ``[]`` senza rieseguire alcuna fase.

    **Razionale Scientifico/Sistemistico**: Garantisce l'idempotenza assoluta
    della ripresa (*resume*) basata esclusivamente sui manifesti e sugli hash
    SHA-256 persistiti su filesystem.
    """
    assert riprendi(eseguita, registro) == []

    # Anche un oggetto nuovo, che sa solo cio' che trova su disco.
    nuova = ProjectRun(scenario.config, passi=_passi(registro))
    assert nuova.completa
    assert riprendi(nuova, registro) == []


def test_la_filogenesi_disattivata_non_e_lavoro_mancante(eseguita):
    """
    **Obiettivo**: Verificare che con ``phylo.enabled = False`` lo stato di
    ``Passo.S9`` sia ``StatoPasso.DISATTIVATA``, non figuri in ``da_eseguire`` e
    non impedisca a ``eseguita.completa`` di essere ``True``.

    **Razionale Scientifico/Sistemistico**: Distingue nettamente una fase
    *disattivata per scelta sperimentale* (``DISATTIVATA``) da una fase *non
    ancora eseguita* (``DA_ESEGUIRE``), evitando che l'assenza dell'albero
    filogenetico venga segnalata come esecuzione incompleta.
    """
    situazione = eseguita.situazione()
    assert situazione[Passo.S9].stato is StatoPasso.DISATTIVATA
    assert "phylo.enabled" in situazione[Passo.S9].motivo
    assert eseguita.disattivate == (Passo.S9,)
    assert Passo.S9 not in eseguita.da_eseguire
    assert eseguita.completa
    assert set(eseguita.completate) == set(TUTTE) - {Passo.S9}


def test_senza_filogenesi_s10_e_calcolata_senza_s9(eseguita):
    """
    **Obiettivo**: Verificare che il manifesto ``manifest_s10.json`` registri in
    ``calcolata_su["a_monte"]`` esattamente ``{"S0", "S7", "S8"}`` quando ``S9`` è disattivata.

    **Razionale Scientifico/Sistemistico**: Traccia con precisione nel manifesto
    di ``S10`` quali fasi hanno contribuito alla costruzione dell'oggetto
    ``phyloseq``, così che un'eventuale attivazione successiva di ``S9`` venga
    immediatamente rilevata come variazione dell'insieme delle dipendenze.
    """
    manifesto = eseguita.albero.manifesto_passo(Passo.S10, Fase.PHYLOSEQ)
    assert set(manifesto.calcolata_su["a_monte"]) == {"S0", "S7", "S8"}


# --------------------------------------------------------------------------- #
# Artefatti cancellati e alterati                                              #
# --------------------------------------------------------------------------- #


def test_cancellando_un_artefatto_si_riparte_da_quella_fase(eseguita, scenario, registro):
    """
    **Obiettivo**: Verificare che cancellando ``s4.json`` la ripresa riesegua
    ``S4`` e tutti i suoi discendenti attivi, lasciando intatte ``S0, S1, S2, S3``.

    **Razionale Scientifico/Sistemistico**: In una corsa su 960 campioni, le
    fasi ``S2`` (``filterAndTrim``) e ``S3`` (``learnErrors``) richiedono tempo
    computazionale significativo; se un file di ``S4`` viene rimosso o va
    ricalcolato, il DAG preserva tutto il lavoro a monte (``S0..S3``) e propaga
    il ricalcolo solo alle fasi a valle che dipendono dalle ASV di ``S4``.
    """
    _file(eseguita, Passo.S4).unlink()

    situazione = eseguita.situazione()
    assert situazione[Passo.S4].stato is StatoPasso.DA_ESEGUIRE
    assert "s4.json" in situazione[Passo.S4].motivo
    assert situazione[Passo.S5].motivo == "a monte da eseguire: S4"
    assert eseguita.prossima() is Passo.S4

    eseguite = riprendi(eseguita, registro)
    assert eseguite == _attive_da(Passo.S4, scenario.config)
    # A monte e di lato non si rifa' nulla.
    assert not {Passo.S0, Passo.S1, Passo.S2, Passo.S3} & set(eseguite)
    assert eseguita.completa


def test_cancellando_un_artefatto_di_s0_si_rifa_tutto(eseguita, scenario, registro):
    """
    **Obiettivo**: Verificare che la cancellazione di ``crosswalk.tsv`` in
    ``01_input_validation/`` provochi la riesecuzione di tutte le fasi attive
    della pipeline a partire da ``S0``.

    **Razionale Scientifico/Sistemistico**: Poiché ``S0`` è la radice del DAG da
    cui dipendono direttamente o transitivamente tutte le fasi successive, la
    perdita del crosswalk invalida l'intera catena analitica.
    """
    cartella = eseguita.albero.cartella(Fase.INPUT_VALIDATION)
    (cartella / "crosswalk.tsv").unlink()
    assert riprendi(eseguita, registro) == list(GRAFO.attive(scenario.config))


def test_ricalcolare_una_fase_senza_dipendenti_non_tocca_le_altre(eseguita, registro):
    """
    **Obiettivo**: Verificare che cancellando l'artefatto di ``S1`` la ripresa
    riesegua esclusivamente ``[Passo.S1]`` senza toccare ``S2..S14``.

    **Razionale Scientifico/Sistemistico**: Dimostra il vantaggio di un vero DAG
    rispetto a una catena sequenziale: ``S1`` produce i profili di qualità QC
    che servono al ricercatore ma non sono letti da ``S2`` (che legge da ``S0``).
    Ricalcolare ``S1`` non deve quindi invalidare il denoising DADA2 in ``S2..S14``.
    """
    _file(eseguita, Passo.S1).unlink()
    assert riprendi(eseguita, registro) == [Passo.S1]


def test_un_artefatto_alterato_fa_rieseguire_la_fase(eseguita, scenario, registro):
    """
    **Obiettivo**: Verificare che aggiungendo anche un solo spazio a ``s8.json``
    lo stato di ``S8`` passi a ``DA_ESEGUIRE`` con motivo ``alterati: s8.json`` e
    la ripresa ricalcoli ``S8`` e i suoi discendenti.

    **Razionale Scientifico/Sistemistico**: La verifica sistematica dei checksum
    SHA-256 impedisce che un file di output modificato manualmente o corrotto su
    disco venga silenziosamente riutilizzato da ``phyloseq`` in ``S10``.
    """
    percorso = _file(eseguita, Passo.S8)
    percorso.write_text(percorso.read_text() + " ", encoding="utf-8")

    situazione = eseguita.situazione()
    assert situazione[Passo.S8].stato is StatoPasso.DA_ESEGUIRE
    assert "alterati: s8.json" in situazione[Passo.S8].motivo

    assert riprendi(eseguita, registro) == _attive_da(Passo.S8, scenario.config)


def test_un_manifesto_alterato_non_vale_come_conclusione(eseguita, registro):
    """
    **Obiettivo**: Verificare che manomettere ``manifest_s5.json`` (alterando
    una metrica) invalidi l'impronta del manifesto e porti ``S5`` a ``DA_ESEGUIRE``.

    **Razionale Scientifico/Sistemistico**: Anche il manifesto stesso è
    autosigillato tramite la propria ``impronta`` crittografica; una modifica non
    autorizzata ai metadati di esecuzione di una fase ne impone il ricalcolo.
    """
    percorso = eseguita.albero.percorso_manifesto_passo(Passo.S5, Fase.SEQTAB)
    documento = json.loads(percorso.read_text())
    documento["metriche"]["righe"] = 999
    percorso.write_text(json.dumps(documento))

    situazione = eseguita.situazione()
    assert situazione[Passo.S5].stato is StatoPasso.DA_ESEGUIRE
    assert "alterato" in situazione[Passo.S5].motivo


def test_una_fase_interrotta_non_risulta_conclusa(eseguita, scenario, registro):
    """
    **Obiettivo**: Verificare che se ``S6`` viene interrotta a metà da
    ``KeyboardInterrupt`` dopo aver scritto un file parziale, il manifesto
    precedente di ``S6`` venga rimosso all'avvio e ``prossima()`` riparta da ``S6``.

    **Razionale Scientifico/Sistemistico**: ``PipelineStep.esegui`` rimuove il
    vecchio ``manifest_<passo>.json`` **prima** di iniziare il calcolo e scrive
    quello nuovo solo al completamento atomico: un'interruzione (Ctrl+C, SIGTERM,
    OOM) non lascia mai un vecchio manifesto valido accanto a file parziali.
    """

    class Interrotta(Doppione):
        passo: ClassVar[Passo] = Passo.S6

        def calcola(self, contesto):
            contesto.albero.scrivi_testo(self.cartella, "parziale.tsv", "meta'")
            raise KeyboardInterrupt

    with pytest.raises(KeyboardInterrupt):
        Interrotta(registro).esegui(eseguita.contesto(Passo.S6))

    assert eseguita.albero.manifesto_passo(Passo.S6, Fase.CHIMERA) is None
    assert eseguita.prossima() is Passo.S6
    assert riprendi(eseguita, registro) == _attive_da(Passo.S6, scenario.config)


# --------------------------------------------------------------------------- #
# Su che cosa e' stata calcolata una fase                                      #
# --------------------------------------------------------------------------- #


def test_cambiando_un_parametro_nessuna_fase_resta_conclusa(eseguita, scenario, registro):
    """
    **Obiettivo**: Verificare che quando una fase non restringe i propri
    parametri (``parametri = None``), la modifica di ``decontam.threshold``
    invalidi tutte le fasi con motivo ``configurazione cambiata``.

    **Razionale Scientifico/Sistemistico**: Il comportamento predefinito è
    massimamente cautelativo: se un passo non dichiara esplicitamente il
    sottoinsieme di sezioni YAML da cui dipende, qualsiasi modifica ai parametri
    scientifici ne forza il ricalcolo.
    """
    cambiata = _variante(scenario.config, decontam__threshold=0.4)
    run = ProjectRun(cambiata, passi=_passi(registro))

    situazione = run.situazione()
    assert run.completate == ()
    assert situazione[Passo.S0].motivo == "configurazione cambiata"
    assert situazione[Passo.S9].stato is StatoPasso.DISATTIVATA

    assert riprendi(run, registro) == list(GRAFO.attive(cambiata))

    # Tornare al valore precedente e' a sua volta un cambiamento.
    tornata = ProjectRun(scenario.config, passi=_passi(registro))
    assert tornata.completate == ()


def test_con_parametri_ristretti_si_rifa_solo_cio_che_ne_dipende(scenario, registro):
    """
    **Obiettivo**: Verificare che con ``RISTRETTI`` la modifica di ``tax.min_boot``
    invalidi ``S8`` e i suoi discendenti ``S10..S14``, lasciando intatte ``S0..S7``
    e la filogenesi ``S9``.

    **Razionale Scientifico/Sistemistico**: Quando ogni fase dichiara il proprio
    gruppo di parametri (es. ``S8 -> ("tax",)``, ``S9 -> ("phylo",)``), cambiare
    la soglia di bootstrap tassonomico ricalcola solo l'assegnazione SILVA in
    ``S8`` e le matrici a valle, risparmiando sia il denoising DADA2 (``S2..S7``)
    sia l'allineamento filogenetico (``S9``).
    """
    config = _variante(scenario.config, phylo__enabled=True)
    run = ProjectRun(config, passi=_passi(registro, RISTRETTI))
    assert riprendi(run, registro) == list(TUTTE)

    cambiata = _variante(config, tax__min_boot=80)
    run = ProjectRun(cambiata, passi=_passi(registro, RISTRETTI))
    situazione = run.situazione()
    assert situazione[Passo.S8].motivo == "configurazione cambiata"
    assert situazione[Passo.S10].motivo == "a monte da eseguire: S8"

    eseguite = riprendi(run, registro)
    assert eseguite == [Passo.S8, Passo.S10, Passo.S11, Passo.S12, Passo.S13, Passo.S14]
    assert Passo.S9 not in eseguite  # la filogenesi non legge la tassonomia


def test_attivare_la_filogenesi_rifa_s9_e_cio_che_ne_dipende(scenario, registro):
    """
    **Obiettivo**: Verificare che passando da ``phylo.enabled = False`` a
    ``True`` vengano eseguite solo ``S9..S14`` (con ``S8`` che resta ``COMPLETATA``),
    e che riportando ``phylo.enabled = False`` vengano rieseguite ``S10..S14``
    per ``dipendenze cambiate``.

    **Razionale Scientifico/Sistemistico**: Attivare l'albero filogenetico in un
    secondo momento deve calcolare ``S9`` e aggiornare l'oggetto ``phyloseq`` in
    ``S10`` (e i controlli/filtri successivi) senza ricalcolare DADA2 o SILVA;
    viceversa, disattivarlo deve rigenerare ``S10`` affinché non contenga più un
    albero obsoleto.
    """
    run = ProjectRun(scenario.config, passi=_passi(registro, RISTRETTI))
    riprendi(run, registro)

    attiva = _variante(scenario.config, phylo__enabled=True)
    run = ProjectRun(attiva, passi=_passi(registro, RISTRETTI))
    situazione = run.situazione()
    assert situazione[Passo.S9].motivo == "nessun manifesto: mai conclusa"
    assert situazione[Passo.S8].stato is StatoPasso.COMPLETATA

    assert riprendi(run, registro) == [
        Passo.S9, Passo.S10, Passo.S11, Passo.S12, Passo.S13, Passo.S14
    ]

    # E disattivarla di nuovo rifa' S10, che non deve piu' contenere l'albero.
    run = ProjectRun(scenario.config, passi=_passi(registro, RISTRETTI))
    assert "dipendenze cambiate" in run.situazione()[Passo.S10].motivo
    assert riprendi(run, registro) == [
        Passo.S10, Passo.S11, Passo.S12, Passo.S13, Passo.S14
    ]


def test_dati_di_ingresso_sostituiti_invalidano_s0_e_il_resto(eseguita, scenario, registro):
    """
    **Obiettivo**: Verificare che l'aggiornamento del timestamp ``mtime`` (o
    della dimensione) di un file FASTQ di ingresso faccia passare ``S0`` a
    ``DA_ESEGUIRE`` con motivo ``dati di ingresso cambiati`` e ricalcoli l'intera pipeline.

    **Razionale Scientifico/Sistemistico**: Se un file ``*.fastq.gz`` o una
    tabella ISA-Tab viene sostituita su disco mantenendo lo stesso nome e la
    stessa configurazione YAML, l'impronta degli ingressi di ``S0`` rileva la
    modifica senza dover ricalcolare l'hash di tutti i 960 FASTQ ad ogni avvio,
    prevenendo l'uso di risultati derivati da file di input non più attuali.
    """
    fastq = sorted(scenario.config.io.fastq_dir.glob(scenario.config.io.fastq_glob))[0]
    stato = fastq.stat()
    os.utime(fastq, ns=(stato.st_atime_ns, stato.st_mtime_ns + 1_000_000_000))

    situazione = eseguita.situazione()
    assert situazione[Passo.S0].motivo == "dati di ingresso cambiati"
    assert riprendi(eseguita, registro) == list(GRAFO.attive(scenario.config))


def test_un_parametro_dichiarato_male_e_un_errore(scenario):
    """
    **Obiettivo**: Verificare che ``impronta_parametri`` sollevi ``KeyError`` su
    nomi di parametri o sezioni inesistenti (``tax.min_bootstrap``, ``tassonomia``)
    e ``ValueError`` se si tenta di includere un parametro escluso (``run.threads``).

    **Razionale Scientifico/Sistemistico**: Evita che un errore di battitura
    nella dichiarazione ``parametri`` di una sottoclasse di ``PipelineStep``
    venga ignorato silenziosamente, impedendo alla fase di accorgersi quando il
    parametro cambia.
    """
    risolta = risolvi(scenario.config)
    with pytest.raises(KeyError, match="tax.min_bootstrap"):
        impronta_parametri(risolta, ("tax.min_bootstrap",))
    with pytest.raises(KeyError, match="tassonomia"):
        impronta_parametri(risolta, ("tassonomia",))
    assert impronta_parametri(risolta, None) == risolta.impronta_risultati
    assert risolta.impronta_risultati != risolta.digest
    with pytest.raises(ValueError, match="run.threads"):
        impronta_parametri(risolta, ("run.threads",))


# --------------------------------------------------------------------------- #
# Parametri che non incidono sui risultati                                     #
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize(
    "variazione",
    [
        {"run__threads": 2},
        {"retry__enabled": False},
        {"retry__max_attempts": 7},
        {"run__threads": 2, "retry__enabled": False, "retry__max_attempts": 7},
    ],
    ids=["threads", "retry.enabled", "retry.max_attempts", "tutti"],
)
def test_un_parametro_escluso_non_invalida_alcuna_fase(
    eseguita, scenario, registro, variazione
):
    """
    **Obiettivo**: Verificare che modificare ``run.threads``, ``retry.enabled``
    o ``retry.max_attempts`` mantenga ``run.completa is True`` e ``riprendi == []``,
    pur cambiando il ``digest`` globale in ``00_config``.

    **Razionale Scientifico/Sistemistico**: Quando si riprende una corsa su una
    macchina o coda HPC con un numero diverso di core (``run.threads``) o con
    una diversa politica di retry, i risultati matematici e biologici delle fasi
    già concluse sono identici: invalidare ore di calcolo DADA2 per un cambio di
    thread sarebbe uno spreco ingiustificato.
    """
    cambiata = _variante(scenario.config, **variazione)
    run = ProjectRun(cambiata, passi=_passi(registro))
    assert run.completa
    assert riprendi(run, registro) == []

    # Il digest scritto in 00_config identifica comunque la configurazione.
    assert risolvi(cambiata).digest != risolvi(scenario.config).digest


def test_spostare_la_cartella_di_output_non_rende_incompleta_l_esecuzione(
    eseguita, scenario, registro, tmp_path
):
    """
    **Obiettivo**: Verificare che copiare l'albero di output in una nuova
    directory e aggiornare ``io.out_root`` mantenga ``run.completa is True`` senza
    rieseguire alcuna fase.

    **Razionale Scientifico/Sistemistico**: L'esclusione di ``io.out_root``
    dall'impronta dei risultati rende l'albero degli artefatti **rilocabile**
    tra workstation, server HPC e archivi di revisione scientifica senza perdere
    lo stato di completamento delle fasi.
    """
    nuova_radice = tmp_path / "spostata"
    shutil.copytree(scenario.config.io.out_root, nuova_radice)
    spostata = _variante(scenario.config, io__out_root=nuova_radice)

    run = ProjectRun(spostata, passi=_passi(registro))
    assert run.completa
    assert riprendi(run, registro) == []


@pytest.mark.parametrize(
    "variazione",
    [{"run__batch_size": 7}, {"run__lockfile": "altro.lock"}, {"decontam__threshold": 0.4}],
    ids=["batch_size", "lockfile", "decontam.threshold"],
)
def test_un_parametro_incluso_invalida_le_fasi(eseguita, scenario, registro, variazione):
    """
    **Obiettivo**: Verificare che la modifica di ``run.batch_size``,
    ``run.lockfile`` o ``decontam.threshold`` invalidi le fasi portando
    ``run.completate`` a ``()``.

    **Razionale Scientifico/Sistemistico**: A differenza di ``run.threads``, il
    parametro ``run.batch_size`` altera il partizionamento dei campioni nei lotti
    di ``learnErrors``/``dada``, ``run.lockfile`` altera le versioni delle
    librerie R/Bioconductor e ``decontam.threshold`` altera la soglia di
    decontaminazione: tutti e tre incidono sui risultati e devono invalidare il
    calcolo precedente.
    """
    cambiata = _variante(scenario.config, **variazione)
    run = ProjectRun(cambiata, passi=_passi(registro))
    assert run.completate == ()
    assert run.situazione()[Passo.S0].motivo == "configurazione cambiata"


def test_da_un_gruppo_i_parametri_esclusi_restano_fuori(scenario):
    """
    **Obiettivo**: Verificare che richiedendo l'impronta dell'intero gruppo
    ``("run",)``, una variazione di ``run.threads`` produca la stessa impronta
    mentre una variazione di ``run.batch_size`` produca un'impronta diversa.

    **Razionale Scientifico/Sistemistico**: Garantisce che una fase che dichiara
    dipendenza dalla sezione ``run`` erediti automaticamente l'esclusione di
    ``run.threads`` senza dover elencare a mano ogni singola chiave di ``run``.
    """
    prima = impronta_parametri(risolvi(scenario.config), ("run",))
    dopo = impronta_parametri(risolvi(_variante(scenario.config, run__threads=2)), ("run",))
    assert prima == dopo
    lotti = impronta_parametri(
        risolvi(_variante(scenario.config, run__batch_size=7)), ("run",)
    )
    assert lotti != prima


# --------------------------------------------------------------------------- #
# Costo della valutazione                                                      #
# --------------------------------------------------------------------------- #


@pytest.fixture
def conteggio_checksum(monkeypatch):
    """Conta i checksum calcolati durante la valutazione dello stato."""
    import amplicon16s.runner.project as progetto

    calcolati: list = []
    originale = progetto.checksum_file

    def contato(percorso):
        calcolati.append(percorso)
        return originale(percorso)

    monkeypatch.setattr(progetto, "checksum_file", contato)
    return calcolati


def test_una_valutazione_calcola_ogni_checksum_una_volta(eseguita, conteggio_checksum):
    r"""
    **Obiettivo**: Verificare che una chiamata a ``eseguita.valuta()`` calcoli
    lo SHA-256 su disco esattamente una volta per ciascuno dei 17 artefatti (4 di
    ``S0`` + 13 dei doppioni) e che interrogazioni successive sulla stessa
    ``valutazione`` non ricalcolino alcun checksum.

    **Razionale Scientifico/Sistemistico**: In una pipeline con artefatti di
    grandi dimensioni (matrici ASV, oggetti ``phyloseq``), ricalcolare gli hash
    SHA-256 ad ogni interrogazione delle proprietà (``prossima()``, ``completate``,
    ``da_eseguire``, ``contesto()``) introdurrebbe un collo di bottiglia I/O
    quadratico; il caching dentro ``Valutazione`` garantisce complessità lineare
    $O(N_{\text{artefatti}})$.
    """
    valutazione = eseguita.valuta()
    artefatti = {
        eseguita.albero.cartella(eseguita.fase(p).cartella) / nome
        for p in valutazione.completate
        for nome in eseguita.albero.manifesto_passo(p, eseguita.fase(p).cartella).nomi
    }
    assert len(artefatti) == 17  # 4 di S0 e uno per ciascuna delle 13 fasi doppione
    assert len(conteggio_checksum) == len(artefatti)
    assert set(conteggio_checksum) == artefatti

    # Tutte le domande poste alla stessa valutazione non ne calcolano altri.
    for _ in range(3):
        valutazione.prossima()
        valutazione.completate
        valutazione.disattivate
        valutazione.da_eseguire
        valutazione.completa
        eseguita.contesto(Passo.S14, valutazione)
    assert len(conteggio_checksum) == len(artefatti)


def test_le_scorciatoie_compiono_una_valutazione_nuova(eseguita, conteggio_checksum):
    """
    **Obiettivo**: Verificare che i metodi scorciatoia su ``ProjectRun`` (come
    ``eseguita.prossima()``) effettuino una nuova lettura del disco intercettando
    la cancellazione immediata di ``s3.json``.

    **Razionale Scientifico/Sistemistico**: Mentre un oggetto ``Valutazione``
    fotografa uno stato coerente in un dato istante, interrogare direttamente
    ``ProjectRun`` dopo una modifica sul filesystem deve riflettere lo stato
    attuale del disco senza rischiare cache stantie (*stale cache*).
    """
    eseguita.prossima()
    primo = len(conteggio_checksum)
    _file(eseguita, Passo.S3).unlink()
    assert eseguita.prossima() is Passo.S3
    assert len(conteggio_checksum) > primo


# --------------------------------------------------------------------------- #
# Cartelle condivise                                                           #
# --------------------------------------------------------------------------- #


def test_la_cartella_completa_non_rende_completa_la_fase_accanto(scenario, registro):
    """
    **Obiettivo**: Verificare che dopo aver eseguito ``S11`` (che scrive in
    ``11_controls/``), ``S12`` (che condivide ``11_controls/``) risulti comunque
    ``StatoPasso.DA_ESEGUIRE`` con motivo ``nessun manifesto: mai conclusa``.

    **Razionale Scientifico/Sistemistico**: Risolve la collisione nelle cartelle
    condivise (``07_chimera``, ``11_controls``, ``12_final``): se il runner si
    basasse sul manifesto di directory ``manifest.json``, il completamento di
    ``S11`` (KatharoSeq) farebbe apparire già completata ``S12`` (``decontam``)
    prima ancora che ``S12`` abbia mai girato, saltando silenziosamente la
    decontaminazione! L'uso di ``manifest_s11.json`` e ``manifest_s12.json``
    elimina alla radice questo bug.
    """
    run = ProjectRun(scenario.config, passi=_passi(registro))
    for passo in GRAFO.attive(scenario.config):
        if passo is Passo.S12:
            break
        run.fase(passo).esegui(run.contesto(passo))

    # Il manifesto della cartella direbbe che 11_controls e' completa...
    assert run.albero.fase_completa(Fase.CONTROLS)
    # ...ma S12 non ha mai girato, e il suo criterio e' il suo manifesto.
    situazione = run.situazione()
    assert situazione[Passo.S11].stato is StatoPasso.COMPLETATA
    assert situazione[Passo.S12].stato is StatoPasso.DA_ESEGUIRE
    assert situazione[Passo.S12].motivo == "nessun manifesto: mai conclusa"


@pytest.mark.parametrize(
    ("prima", "dopo", "cartella"),
    [
        (Passo.S6, Passo.S7, Fase.CHIMERA),
        (Passo.S11, Passo.S12, Fase.CONTROLS),
        (Passo.S13, Passo.S14, Fase.FINAL),
    ],
)
def test_due_fasi_nella_stessa_cartella_si_concludono_separatamente(
    eseguita, registro, prima, dopo, cartella
):
    """
    **Obiettivo**: Verificare che per tutte e tre le cartelle condivise
    (``CHIMERA``, ``CONTROLS``, ``FINAL``) la cancellazione dell'artefatto della
    seconda fase (``dopo``) faccia rieseguire solo ``dopo`` lasciando ``prima``
    in stato ``COMPLETATA`` con impronta invariata.

    **Razionale Scientifico/Sistemistico**: Garantisce che i due passi
    coinquilini della stessa directory abbiano cicli di vita e invalidazioni
    completamente ortogonali.
    """
    albero = eseguita.albero
    assert albero.manifesto_passo(prima, cartella).nomi == (f"{prima.lower()}.json",)
    assert albero.manifesto_passo(dopo, cartella).nomi == (f"{dopo.lower()}.json",)

    # Rifare la seconda non tocca la conclusione della prima.
    impronta_prima = albero.manifesto_passo(prima, cartella).impronta
    _file(eseguita, dopo).unlink()
    situazione = eseguita.situazione()
    assert situazione[prima].stato is StatoPasso.COMPLETATA
    assert situazione[dopo].stato is StatoPasso.DA_ESEGUIRE

    eseguite = riprendi(eseguita, registro)
    assert prima not in eseguite and dopo in eseguite
    assert albero.manifesto_passo(prima, cartella).impronta == impronta_prima


def test_una_fase_che_riscrive_il_file_dell_altra_la_invalida(eseguita):
    """
    **Obiettivo**: Verificare che se un processo sovrascrive ``s11.json`` dentro
    ``11_controls/``, lo stato di ``Passo.S11`` diventi immediatamente ``DA_ESEGUIRE``.

    **Razionale Scientifico/Sistemistico**: Poiché ``manifest_s11.json``
    custodisce lo SHA-256 originale di ``s11.json``, qualsiasi sovrascrittura
    accidentale da parte di un'altra fase nella medesima cartella viene
    immediatamente intercettata come violazione di integrità.
    """
    albero = eseguita.albero
    albero.scrivi_testo(Fase.CONTROLS, "s11.json", "sovrascritto da S12")
    assert eseguita.situazione()[Passo.S11].stato is StatoPasso.DA_ESEGUIRE


def test_i_manifesti_di_fase_hanno_nomi_distinti():
    """
    **Obiettivo**: Verificare che ``nome_manifesto_passo("S11")`` sia distinto
    sia da ``nome_manifesto_passo("S12")`` sia dal manifesto cumulativo ``manifest.json``.

    **Razionale Scientifico/Sistemistico**: Impedisce collisioni di nome su
    filesystem tra i manifesti atomici dei singoli passi e il manifesto di directory.
    """
    assert nome_manifesto_passo("S11") != nome_manifesto_passo("S12")
    assert nome_manifesto_passo("S11") != "manifest.json"


# --------------------------------------------------------------------------- #
# Prerequisiti                                                                 #
# --------------------------------------------------------------------------- #


def test_s13_prima_di_s12_e_rifiutata(eseguita):
    """
    **Obiettivo**: Verificare che invocare ``S13.esegui`` dopo aver rimosso il
    manifesto di ``S12`` sollevi ``ErroreRevisioneUmana`` con codice specifico
    ``E-S13-01`` prima ancora di rimuovere il manifesto esistente di ``S13``.

    **Razionale Scientifico/Sistemistico**: Impone a livello di esecuzione del
    singolo passo la guardia metodologica ``E-S13-01`` (divieto assoluto di
    filtrare per prevalenza senza aver concluso la decontaminazione ``S12``),
    proteggendo al contempo gli artefatti preesistenti di ``S13``.
    """
    eseguita.albero.rimuovi_manifesto_passo(Passo.S12, Fase.CONTROLS)

    with pytest.raises(ErroreRevisioneUmana) as info:
        eseguita.fase(Passo.S13).esegui(eseguita.contesto(Passo.S13))
    assert info.value.codice == "E-S13-01"
    assert info.value.contesto["mancanti"] == ["S12"]
    # Il rifiuto avviene prima di toccare alcunche'.
    assert eseguita.albero.manifesto_passo(Passo.S13, Fase.FINAL) is not None


def test_nessuna_fase_gira_prima_delle_sue_dipendenze(scenario, registro):
    """
    **Obiettivo**: Verificare che tentare di eseguire ``S2`` quando ``S0`` non è
    ancora stata eseguita sollevi ``ErroreRevisioneUmana`` con codice ``E-GRAFO-01``.

    **Razionale Scientifico/Sistemistico**: Garantisce che nessun passo della
    pipeline possa mai essere invocato fuori ordine topologico quando le sue
    dipendenze a monte non sono state completate e certificate.
    """
    run = ProjectRun(scenario.config, passi=_passi(registro))
    with pytest.raises(ErroreRevisioneUmana, match="S0") as info:
        run.fase(Passo.S2).esegui(run.contesto(Passo.S2))
    assert info.value.codice == "E-GRAFO-01"
    assert info.value.fase == "GRAFO"
    assert registro == []


def test_s13_senza_s11_ne_s12_porta_il_codice_della_precedenza(eseguita):
    """
    **Obiettivo**: Verificare che quando mancano sia ``S11`` sia ``S12``, il
    tentativo di eseguire ``S13`` sollevi il codice metodologico specifico
    ``E-S13-01`` (e non quello generico ``E-GRAFO-01``).

    **Razionale Scientifico/Sistemistico**: Dà priorità diagnostica alla
    violazione metodologica specifica (``S12 -> S13``, ``E-S13-01``) rispetto
    alla violazione generica del grafo, fornendo all'operatore il messaggio
    scientifico più pertinente.
    """
    for passo in (Passo.S11, Passo.S12):
        eseguita.albero.rimuovi_manifesto_passo(passo, Fase.CONTROLS)
    contesto = eseguita.contesto(Passo.S13)
    with pytest.raises(ErroreRevisioneUmana) as info:
        eseguita.fase(Passo.S13).esegui(contesto)
    assert info.value.codice == "E-S13-01"


def test_s12_senza_s11_porta_il_codice_generico(eseguita):
    """
    **Obiettivo**: Verificare che tentare di eseguire ``S12`` in assenza di
    ``S11`` sollevi ``ErroreRevisioneUmana`` con il codice generico ``E-GRAFO-01``.

    **Razionale Scientifico/Sistemistico**: Conferma che ``E-S13-01`` è
    riservato esclusivamente alla coppia ``(S12, S13)``, mentre tutte le altre
    violazioni di dipendenza nel DAG ricadono sotto ``E-GRAFO-01``.
    """
    eseguita.albero.rimuovi_manifesto_passo(Passo.S11, Fase.CONTROLS)
    with pytest.raises(ErroreRevisioneUmana) as info:
        eseguita.fase(Passo.S12).esegui(eseguita.contesto(Passo.S12))
    assert info.value.codice == "E-GRAFO-01"


def test_una_precedenza_con_un_codice_inesistente_e_rifiutata():
    """
    **Obiettivo**: Verificare che dichiarare in ``Grafo`` una precedenza
    associata a un codice non censito in ``CATALOGO`` (``E-S99-01``) sollevi ``KeyError``.

    **Razionale Scientifico/Sistemistico**: Impedisce di registrare nel grafo
    regole di precedenza che al momento della violazione fallirebbero per
    mancanza del codice nel catalogo degli errori.
    """
    with pytest.raises(KeyError, match="E-S99-01"):
        Grafo(_nodi_base(), precedenze=((Passo.S12, Passo.S13, "E-S99-01"),))


# --------------------------------------------------------------------------- #
# S0 sulla classe base, e l'esecuzione reale di oggi                           #
# --------------------------------------------------------------------------- #


def test_s0_invocata_da_sola_risulta_conclusa_per_l_esecuzione(scenario):
    """
    **Obiettivo**: Verificare che l'invocazione diretta di ``esegui_s0(scenario.config)``
    scriva anche il manifesto del passo ``manifest_s0.json``, facendo risultare
    ``Passo.S0`` in stato ``COMPLETATA`` per ``ProjectRun``.

    **Razionale Scientifico/Sistemistico**: Garantisce l'interoperabilità tra il
    comando standalone ``amplicon16s validate`` (che invoca ``esegui_s0``) e una
    successiva esecuzione con ``amplicon16s resume``, che riconosce ``S0`` come
    già conclusa senza ripeterla.
    """
    esegui_s0(scenario.config)
    run = ProjectRun(scenario.config)
    assert run.situazione()[Passo.S0].stato is StatoPasso.COMPLETATA


def test_s0_non_superata_non_viene_registrata_come_conclusa(tmp_path):
    """
    **Obiettivo**: Verificare che quando ``ValidazioneIngressi(solleva=False)``
    termina con un gate fallito, restituisca ``Esito.NON_SUPERATA``, non scriva
    ``manifest_s0.json`` e mantenga ``run.prossima() is Passo.S0``.

    **Razionale Scientifico/Sistemistico**: Anche se ``gates.json`` viene scritto
    a scopo diagnostico, una validazione S0 fallita non deve mai sigillare il
    manifesto di completamento, altrimenti una ripresa successiva salterebbe S0
    procedendo su dati invalidi.
    """
    scenario = crea_scenario(
        tmp_path, _campioni(), con_letture=True, sovrascrivi={"run": {"threads": 100000}}
    )
    fase = ValidazioneIngressi(solleva=False)
    run = ProjectRun(scenario.config)
    risultato = fase.esegui(run.contesto(Passo.S0))

    assert risultato.esito is Esito.NON_SUPERATA
    assert risultato.impronta is None
    assert run.albero.manifesto_passo(Passo.S0, Fase.INPUT_VALIDATION) is None
    assert run.prossima() is Passo.S0


def test_il_risultato_porta_artefatti_con_checksum_e_metriche(scenario):
    """
    **Obiettivo**: Verificare che ``run.fase(Passo.S0).esegui(...)`` restituisca
    un ``RisultatoPasso`` completo dei 4 artefatti con checksum ``sha256:...``,
    metrica ``campioni == 4`` e ``impronta`` coincidente con ``manifest_s0.json``.

    **Razionale Scientifico/Sistemistico**: Certifica che il Template Method di
    ``PipelineStep`` calcoli e restituisca gli hash crittografici e le metriche
    quantitative esattamente come vengono salvati su disco.
    """
    run = ProjectRun(scenario.config)
    risultato = run.fase(Passo.S0).esegui(run.contesto(Passo.S0))

    assert risultato.completata
    assert {a.nome for a in risultato.artefatti} == {
        "gates.json", "crosswalk.tsv", "inventario.json", "letture_ispezionate.tsv"
    }
    assert all(a.checksum.startswith("sha256:") for a in risultato.artefatti)
    assert risultato.metriche["campioni"] == 4
    manifesto = run.albero.manifesto_passo(Passo.S0, Fase.INPUT_VALIDATION)
    assert manifesto.impronta == risultato.impronta


def test_l_inventario_e_riletto_dall_artefatto_di_s0(scenario):
    r"""
    **Obiettivo**: Verificare che dopo l'esecuzione di ``S0``, ``run.inventario``
    venga ricostruito automaticamente da ``01_input_validation/inventario.json``
    valorizzando ``run.risolta.derivati.prev_min_samples``.

    **Razionale Scientifico/Sistemistico**: Il parametro derivato dinamico
    ``prev_min_samples`` ($\lceil \text{prev.min\_fraction} \times N_{\text{bio}} \rceil$)
    dipende dal conteggio dei campioni biologici scoperto in ``S0``; rileggerlo
    da ``inventario.json`` permette alle fasi successive (e alle riprese su
    nuovi processi) di conoscere il denominatore di prevalenza senza rieseguire ``S0``.
    """
    run = ProjectRun(scenario.config)
    assert run.inventario is None
    esegui_s0(scenario.config)

    assert run.inventario == Contesto(scenario.config).inventario
    assert run.risolta.derivati.prev_min_samples is not None


def test_oggi_esiste_solo_s0(scenario):
    """
    **Obiettivo**: Verificare che allo stato della Settimana 9 ``passi_realizzati()``
    contenga ``{Passo.S0}``, mentre ``Passo.S1`` sia marcata ``StatoPasso.NON_REALIZZATA``
    e sollevi ``LookupError`` se richiesta a ``run.fase(Passo.S1)``.

    **Razionale Scientifico/Sistemistico**: Separa in modo trasparente le fasi
    già implementate nel codice di produzione (``S0`` al termine della Fase F3)
    dalle fasi successive (``S1..S14``, pianificate nelle Fasi F4–F7), evitando
    falsi stati di completamento.
    """
    run = ProjectRun(scenario.config)
    assert set(passi_realizzati()) == {Passo.S0}
    esegui_s0(scenario.config)

    situazione = run.situazione()
    assert situazione[Passo.S0].stato is StatoPasso.COMPLETATA
    assert situazione[Passo.S1].stato is StatoPasso.NON_REALIZZATA
    assert run.prossima() is Passo.S1
    assert not run.completa
    with pytest.raises(LookupError, match="S1"):
        run.fase(Passo.S1)


def test_l_albero_e_quello_della_configurazione(scenario):
    """
    **Obiettivo**: Verificare che ``run.albero.radice`` coincida con il percorso
    risolto di ``AlberoOutput(scenario.config.io.out_root).radice``.

    **Razionale Scientifico/Sistemistico**: Garantisce che ``ProjectRun`` operi
    sempre sulla directory di output dichiarata nella configurazione validata.
    """
    run = ProjectRun(scenario.config)
    assert run.albero.radice == AlberoOutput(scenario.config.io.out_root).radice
