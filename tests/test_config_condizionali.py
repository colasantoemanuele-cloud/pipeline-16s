"""Parametri condizionali della configurazione.

Un parametro condizionale (``defaults.CONDIZIONALI``) e' obbligatorio solo
quando la condizione che lo rende pertinente e' vera: senza controlli positivi
non serve dire quanti ne occorrono per una curva, senza controlli negativi non
serve la soglia della decontaminazione, senza filogenesi non serve il limite
delle varianti. Il modulo verifica, per ogni condizione:

* condizione vera e parametro mancante: G15 respinge con ``E-G15-10`` e il
  messaggio dice quale condizione lo richiede;
* condizione falsa e parametro omesso: la configurazione e' accettata, nessun
  calcolo usa il parametro e la catena arriva in fondo con gli stessi conteggi,
  tassonomia, campioni e varianti;
* condizione falsa e parametro dichiarato: la configurazione e' accettata
  senza effetto sui risultati, con il digest di sempre.

Alcuni condizionali una fase li riporta in un riepilogo anche con la
condizione falsa (``defaults.CONDIZIONALI_RIPORTATI``): dichiarati, il
riepilogo ne porta il valore; omessi, il campo resta e vale nullo. Gli altri
senza valore di assenza non si leggono affatto: il modulo verifica la guardia
che lo impedisce, l'origine ``non pertinente`` nel report e la specie
``condizionale`` nello schema esportato.

I test che eseguono la catena richiedono R con i pacchetti di calcolo:

    docker run --rm -e PYTHONPATH=/app/src -e AMPLICON16S_R_DIR=/app/R \\
      -v "$(pwd)":/app -w /app <immagine> pytest tests/test_config_condizionali.py
"""

from __future__ import annotations

import copy
import csv
import functools
import json
import os
from pathlib import Path
from typing import Any

import pytest
import sintetico
import yaml
from conftest import esegui_g15
from sottoinsieme import motivo_pacchetti_r_assenti

from amplicon16s import cli, interfaccia
from amplicon16s.config import defaults
from amplicon16s.config.resolve import parametri_dichiarati, risolvi
from amplicon16s.config.schema import (
    INTESTAZIONE_MANCANTI,
    ErroreConfigurazione,
    obbligatori_mancanti,
    valida,
)
from amplicon16s.config.vista import CondizionaleNonPertinente, VistaConfig
from amplicon16s.gates.g01_g15 import ErroreGate
from amplicon16s.io_layer.artifacts import Fase
from amplicon16s.logging.logger import chiudi
from amplicon16s.metadata.models import ClasseCampione
from amplicon16s.report.builder import CARTELLA_REPORT, CARTELLA_TABELLE, genera
from amplicon16s.runner.executor import Conclusione, Esecutore
from amplicon16s.runner.project import ProjectRun
from amplicon16s.runner.retry import applica
from amplicon16s.steps.s00_validate import esegui_s0

#: Il digest della configurazione sintetica senza colonna del tipo, con ogni
#: condizionale dichiarato, calcolato prima che i condizionali esistessero: e'
#: scritto per esteso perche' il test deve dimostrare che dichiararli non ha
#: cambiato l'identita' di una configurazione. Cambia, legittimamente, solo
#: quando lo schema acquista o perde un parametro.
DIGEST_SENZA_TIPO = "sha256:ef247324c1983224689b3e895ca8d4bdebd36b6aaa0554f8666e40e2b452ffa6"
#: Lo stesso, per la configurazione con la colonna del tipo e senza controlli
#: positivi.
DIGEST_CON_TIPO = "sha256:4caf7e25c720caf27023f7029330389a18eec60528f198afbe38d58184a1a6e6"


@pytest.fixture(autouse=True)
def uscite_pulite():
    """Chiude le uscite del log prima e dopo ogni test."""
    chiudi()
    yield
    chiudi()


@functools.cache
def _sonda_bioc() -> str | None:
    """Perche' l'ambiente R richiesto non c'e', o ``None``."""
    return motivo_pacchetti_r_assenti(
        "dada2", "ShortRead", "jsonlite", "phyloseq", "Biostrings", "decontam"
    )


@pytest.fixture
def bioc():
    """Richiede R con i pacchetti di tutte le fasi: salta senza, ma in CI fallisce."""
    if _sonda_bioc() is not None:
        if os.environ.get("AMPLICON16S_RICHIEDI_BIOC") == "1":
            pytest.fail(f"Bioconductor e' richiesto in questo ambiente: {_sonda_bioc()}")
        pytest.skip(_sonda_bioc())


#: Un percorso fisso del riferimento: quello vero dipende da dove sta il clone,
#: ed entrerebbe nel digest.
_RIFERIMENTO = {"ref_fasta": "/dati/riferimento.fa.gz"}


def _dataset_fittizio(con_tipo: bool) -> sintetico.DatasetSintetico:
    """La descrizione di un dataset che non esiste su disco: basta alla sola
    validazione della configurazione, e ha percorsi fissi (il digest non
    dipende dalla cartella temporanea del test).
    """
    return sintetico.DatasetSintetico(
        Path("/dati"), Path("/dati/letture"), Path("/dati/campioni.tsv"), None, {}, (),
        None, 159, con_tipo)


def _senza_tipo() -> dict[str, Any]:
    """La configurazione completa di un dataset senza colonna del tipo, con il
    filtro di prevalenza attivo: quella dei digest scritti per esteso.
    """
    return sintetico.configurazione(_dataset_fittizio(False), Path("/uscita"), tax=_RIFERIMENTO)


def _tutte_false() -> dict[str, Any]:
    """Una configurazione completa in cui ogni condizione e' falsa e ogni
    condizionale e' dichiarato: nessuna colonna del tipo, nessun controllo,
    nessuna riclassificazione, filogenesi, prevalenza e modulo derivato spenti.
    """
    return sintetico.configurazione(
        _dataset_fittizio(False), Path("/uscita"), tax=_RIFERIMENTO, prev={"apply": False})


def _tutte_vere() -> dict[str, Any]:
    """Una configurazione completa in cui ogni condizione e' vera."""
    return sintetico.configurazione(
        _dataset_fittizio(True), Path("/uscita"), tax=_RIFERIMENTO,
        meta={"derive_module": True, "module_regex": "^([A-Z]+)"},
        phylo={"enabled": True},
        ctrl={"positive_values": ["positivo"], "blank_values": ["bianco"],
              "blank_override_column": sintetico.COLONNA_SITO,
              "blank_override_values": ["mai aperto"]},
        katharoseq={"target_taxon": "Bacillus", "cell_count_column": "cellule"},
    )


def _senza(dati: dict[str, Any], *chiavi: str) -> dict[str, Any]:
    """Una copia della configurazione senza i parametri indicati."""
    copia = copy.deepcopy(dati)
    for chiave in chiavi:
        gruppo, nome = chiave.split(".")
        del copia[gruppo][nome]
    return copia


# --------------------------------------------------------------------------- #
# L'elenco                                                                     #
# --------------------------------------------------------------------------- #


def test_ogni_condizionale_ha_una_condizione_su_un_parametro_dello_schema():
    """
    **Obiettivo**: Verificare che ogni condizionale sia un parametro dello
    schema fuori dagli obbligatori, che i parametri che decidono la condizione
    esistano, che i condizionali senza un valore di assenza e quelli riportati
    nei riepiloghi siano due insiemi distinti di condizionali, e che omessi
    valgano tutti nullo (gli altri: nullo o elenco vuoto).

    **Razionale scientifico e sistemistico**: L'elenco e' la sola fonte da cui
    validazione, schema esportato, report e ``init`` leggono le condizioni: un
    nome sbagliato renderebbe un parametro omissibile sempre.
    """
    chiavi = interfaccia.chiavi_dello_schema()
    for chiave, alternative in defaults.CONDIZIONALI.items():
        assert chiave in chiavi and chiave not in defaults.OBBLIGATORI and alternative, chiave
        for parametro, forma in alternative:
            assert parametro in chiavi, chiave
            assert forma in (defaults.NON_NULLO, defaults.NON_VUOTO, defaults.VERO)
            assert parametro in defaults.testo_condizione(chiave)
        assert defaults.testo_condizione(chiave).count(" oppure ") == len(alternative) - 1
    assert defaults.CONDIZIONALI_SENZA_VALORE < set(defaults.CONDIZIONALI)
    assert defaults.CONDIZIONALI_RIPORTATI < set(defaults.CONDIZIONALI)
    assert not defaults.CONDIZIONALI_SENZA_VALORE & defaults.CONDIZIONALI_RIPORTATI
    assert defaults.CONDIZIONALI_NULLI_SE_OMESSI == (
        defaults.CONDIZIONALI_SENZA_VALORE | defaults.CONDIZIONALI_RIPORTATI)
    config = valida(_senza(_tutte_false(), *defaults.CONDIZIONALI))
    for chiave in defaults.CONDIZIONALI:
        gruppo, nome = chiave.split(".")
        valore = getattr(getattr(config, gruppo), nome)
        if chiave in defaults.CONDIZIONALI_NULLI_SE_OMESSI:
            assert valore is None, chiave
        else:
            assert valore in (None, []), chiave


def test_la_decontaminazione_e_richiesta_da_una_delle_due_vie_dei_negativi():
    """
    **Obiettivo**: Verificare che i quattro parametri della decontaminazione
    siano richiesti quando i controlli negativi sono dichiarati con le
    etichette, quando lo sono con la sola riclassificazione (con e senza la
    colonna del tipo), e non lo siano quando nessuna delle due vie ne dichiara;
    che il messaggio di ``E-G15-10`` riporti la disgiunzione; e che la
    validazione iniziale conti come controlli negativi i campioni
    riclassificati, cioe' che le fasi vedano dei negativi esattamente quando
    la condizione e' vera.

    **Razionale scientifico e sistemistico**: Un campione diventa controllo
    negativo per la sua etichetta o per la regola di riclassificazione. Se la
    condizione guardasse le sole etichette, un dataset con i negativi dati
    dalla riclassificazione potrebbe omettere la soglia, e la decontaminazione
    partirebbe senza.
    """
    quattro = ["decontam.threshold", "decontam.min_blanks", "decontam.mode",
               "qc.max_frac_contaminant"]
    riclassificati = {"blank_override_column": sintetico.COLONNA_SITO,
                      "blank_override_values": ["sito 1"]}
    casi = {
        "etichette": (True, {}),
        "riclassificazione con la colonna del tipo": (True, {"blank_values": [], **riclassificati}),
        "riclassificazione senza la colonna del tipo": (False, riclassificati),
    }
    for nome, (con_tipo, ctrl) in casi.items():
        completa = sintetico.configurazione(
            _dataset_fittizio(con_tipo), Path("/uscita"), tax=_RIFERIMENTO,
            ctrl={"blank_values": ["bianco"] if con_tipo else [], **ctrl})
        esegui_g15(completa)
        for chiave in quattro:
            assert defaults.condizione_vera(chiave, {
                f"ctrl.{n}": v for n, v in completa["ctrl"].items()}), (nome, chiave)
        assert obbligatori_mancanti(_senza(completa, *quattro)) == quattro, nome
        with pytest.raises(ErroreGate) as rifiuto:
            esegui_g15(_senza(completa, *quattro))
        assert rifiuto.value.codice == "E-G15-10", nome
        assert ("decontam.threshold perche' ctrl.blank_values non e' vuoto oppure "
                "ctrl.blank_override_values non e' vuoto") in rifiuto.value.violazioni[0].dettaglio
    # Nessuna delle due vie: i quattro non sono richiesti.
    for con_tipo in (True, False):
        nessuno = sintetico.configurazione(
            _dataset_fittizio(con_tipo), Path("/uscita"), tax=_RIFERIMENTO,
            ctrl={"blank_values": []})
        assert obbligatori_mancanti(_senza(nessuno, *quattro)) == []
        esegui_g15(_senza(nessuno, *quattro))


def test_i_campioni_riclassificati_sono_i_negativi_che_le_fasi_vedono(tmp_path):
    """
    **Obiettivo**: Verificare su un dataset sintetico senza colonna del tipo
    che, con la sola regola di riclassificazione, la validazione iniziale
    classifichi dei controlli negativi nell'inventario, e che senza la regola
    non ne classifichi alcuno: con i quattro parametri della decontaminazione
    omessi la fase si conclude e l'avviso sui controlli assenti non cita un
    minimo che non e' dichiarato.

    **Razionale scientifico e sistemistico**: La condizione letta dalla
    configurazione e i campioni che le fasi trattano da controlli negativi
    devono coincidere: e' cio' che rende lecito, con la condizione falsa,
    passare alla decontaminazione dei parametri nulli.
    """
    dataset = sintetico.genera(tmp_path / "dati", biologici=6, letture=40, con_tipo=False)
    quattro = ["decontam.threshold", "decontam.min_blanks", "decontam.mode",
               "qc.max_frac_contaminant"]
    con_regola = sintetico.configurazione(
        dataset, tmp_path / "con", ctrl={"blank_override_column": sintetico.COLONNA_SITO,
                                         "blank_override_values": ["sito 1"]})
    esito = esegui_s0(valida(con_regola))
    negativi = [c for c in esito.inventario if c.classe is ClasseCampione.CONTROLLO_NEGATIVO]
    assert esito.superata and len(negativi) == 2
    chiudi()
    senza_regola = _senza(sintetico.configurazione(dataset, tmp_path / "senza"), *quattro)
    esito = esegui_s0(valida(senza_regola))
    assert esito.superata
    assert {c.classe for c in esito.inventario} == {ClasseCampione.BIOLOGICO}
    (avviso,) = [a for a in esito.avvisi if "controlli negativi" in a.dettaglio]
    assert avviso.codice == "E-S0-17"
    assert avviso.dettaglio == (
        "il dataset ha 0 controlli negativi (ctrl.column e' nullo): i contaminanti non "
        "sono stimabili dai controlli")
    # Dichiarato il minimo, il messaggio e' quello di sempre e lo cita.
    chiudi()
    esito = esegui_s0(valida(sintetico.configurazione(dataset, tmp_path / "dichiarato")))
    (avviso,) = [a for a in esito.avvisi if "controlli negativi" in a.dettaglio]
    assert avviso.dettaglio == (
        "il dataset ha 0 controlli negativi, meno dei 2 di decontam.min_blanks "
        "(ctrl.column e' nullo): i contaminanti non sono stimabili dai controlli")



# --------------------------------------------------------------------------- #
# Condizione vera: il parametro e' obbligatorio                                #
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize("chiave", sorted(defaults.CONDIZIONALI))
def test_con_la_condizione_vera_il_condizionale_mancante_e_respinto(chiave):
    """
    **Obiettivo**: Verificare, per ogni condizionale, che ometterlo da una
    configurazione in cui la sua condizione e' vera faccia respingere la
    configurazione da G15 con ``E-G15-10``, in una sola violazione che lo
    elenca fra i mancanti e dice quale condizione lo richiede; e che il
    documento JSON di ``validate`` riporti la stessa condizione.

    **Razionale scientifico e sistemistico**: Un condizionale e' un
    obbligatorio a tutti gli effetti quando serve: senza il rifiuto la fase
    che lo usa riceverebbe un valore assente a calcolo iniziato.
    """
    completa = _tutte_vere()
    esegui_g15(completa)
    dati = _senza(completa, chiave)
    assert obbligatori_mancanti(dati) == [chiave]
    with pytest.raises(ErroreGate) as rifiuto:
        esegui_g15(dati)
    assert rifiuto.value.gate == "G15" and rifiuto.value.codice == "E-G15-10"
    (violazione,) = rifiuto.value.violazioni
    assert violazione.dettaglio.startswith(f"{INTESTAZIONE_MANCANTI} (1): {chiave}.")
    assert f"{chiave} perche' {defaults.testo_condizione(chiave)}" in violazione.dettaglio
    (problema,) = cli._problemi_di([violazione], gate="G15")
    assert problema["parametro"] == chiave and problema["codice"] == "E-G15-10"
    assert defaults.testo_condizione(chiave) in problema["messaggio"]


@pytest.mark.parametrize("chiave", sorted(defaults.CONDIZIONALI_NULLI_SE_OMESSI))
def test_il_nullo_di_un_condizionale_senza_valore_vale_come_l_omissione(chiave):
    """
    **Obiettivo**: Verificare che un condizionale senza valore di assenza, o
    riportato nei riepiloghi, scritto nullo, sia respinto come mancante quando
    la condizione e' vera e accettato come non dichiarato quando e' falsa.

    **Razionale scientifico e sistemistico**: Una configurazione riletta da
    una sua copia (un nuovo tentativo, un programma esterno) riporta nulli i
    parametri omessi: il nullo non puo' ne' diventare un valore ne' essere un
    errore dove il parametro non serve.
    """
    gruppo, nome = chiave.split(".")
    vera = _tutte_vere()
    vera[gruppo][nome] = None
    with pytest.raises(ErroreConfigurazione, match=INTESTAZIONE_MANCANTI) as rifiuto:
        valida(vera)
    assert chiave in str(rifiuto.value)
    falsa = _tutte_false()
    falsa[gruppo][nome] = None
    assert getattr(getattr(valida(falsa), gruppo), nome) is None


def test_i_mancanti_si_elencano_tutti_insieme_con_le_loro_condizioni():
    """
    **Obiettivo**: Verificare che obbligatori e condizionali mancanti escano
    in una sola voce, i primi senza condizione e i secondi con la propria, e
    che un condizionale la cui condizione dipende da un parametro a sua volta
    mancante non sia richiesto.

    **Razionale scientifico e sistemistico**: Chi compila una configurazione
    deve vedere in una volta tutto cio' che manca; e chiedere i parametri dei
    controlli positivi a chi non ha ancora detto se ne ha sarebbe chiedere un
    numero inventato.
    """
    dati = _senza(_tutte_vere(), "filter.truncLen", "ctrl.min_positives", "phylo.max_seqs")
    with pytest.raises(ErroreConfigurazione) as rifiuto:
        valida(dati)
    (voce,) = [p for p in rifiuto.value.problemi if p.startswith(INTESTAZIONE_MANCANTI)]
    assert voce.startswith(
        f"{INTESTAZIONE_MANCANTI} (3): filter.truncLen, ctrl.min_positives, phylo.max_seqs.")
    assert "ctrl.min_positives perche' ctrl.positive_values non e' vuoto" in voce
    assert "phylo.max_seqs perche' phylo.enabled e' vero" in voce
    assert "filter.truncLen perche'" not in voce
    assert len(rifiuto.value.problemi) == 1

    # Senza ctrl.column le etichette non sono richieste, e senza etichette dei
    # positivi non lo sono i parametri dei controlli positivi.
    soli = {"io": {"fastq_dir": "a", "assay_table": "b", "out_root": "c"},
            "tax": {"ref_fasta": "d", "ref_md5": "0" * 32}}
    assert obbligatori_mancanti(soli) == list(defaults.OBBLIGATORI)


# --------------------------------------------------------------------------- #
# Condizione falsa: omesso o dichiarato                                        #
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize("chiave", sorted(defaults.CONDIZIONALI))
def test_con_la_condizione_falsa_il_condizionale_si_omette(chiave):
    """
    **Obiettivo**: Verificare, per ogni condizionale, che ometterlo da una
    configurazione in cui la sua condizione e' falsa lasci la configurazione
    valida per G15, che il parametro non risulti dichiarato, e che la vista di
    una fase sollevi ``CondizionaleNonPertinente`` alla lettura di un
    condizionale senza valore di assenza (gli altri si leggono vuoti o nulli:
    fra questi i riportati nei riepiloghi, che la fase legge per scriverli).

    **Razionale scientifico e sistemistico**: Un parametro che nessuna fase
    usa non va chiesto. La guardia rende verificabile che nessuna fase lo usi:
    una lettura senza aver guardato la condizione e' un errore, non un nullo
    che scorre nel calcolo.
    """
    dati = _senza(_tutte_false(), chiave)
    assert obbligatori_mancanti(dati) == []
    config = esegui_g15(dati).config
    assert chiave not in parametri_dichiarati(config)
    gruppo, nome = chiave.split(".")
    for dichiarati in ((gruppo,), (chiave,)):
        vista = VistaConfig(config, dichiarati)
        if chiave in defaults.CONDIZIONALI_SENZA_VALORE:
            with pytest.raises(CondizionaleNonPertinente, match=chiave):
                getattr(getattr(vista, gruppo), nome)
        else:
            assert getattr(getattr(vista, gruppo), nome) in (None, [])


def test_un_condizionale_dichiarato_con_la_condizione_falsa_non_cambia_il_digest():
    """
    **Obiettivo**: Verificare che una configurazione che dichiara ogni
    condizionale con la condizione falsa sia accettata e abbia il digest che
    aveva prima che i condizionali esistessero (scritto per esteso), con e
    senza la colonna del tipo; che ometterli dia un digest diverso solo dove
    il valore dichiarato non era quello di assenza; e che la vista di una fase
    li lasci leggere.

    **Razionale scientifico e sistemistico**: Le configurazioni gia' scritte
    dichiarano tutto e identificano esecuzioni gia' consegnate: la regola
    nuova non puo' cambiarne l'identita' ne' respingerle.
    """
    completa = _senza_tipo()
    config = valida(completa)
    assert risolvi(config).digest == DIGEST_SENZA_TIPO
    con_tipo = sintetico.configurazione(_dataset_fittizio(True), Path("/uscita"), tax=_RIFERIMENTO)
    assert risolvi(valida(con_tipo)).digest == DIGEST_CON_TIPO
    assert set(defaults.CONDIZIONALI) <= set(parametri_dichiarati(config))
    assert VistaConfig(config, ("ctrl", "phylo")).ctrl.min_positives == 3
    assert VistaConfig(config, ("ctrl", "phylo")).phylo.max_seqs == 1000

    # I valori di assenza coincidono con cio' che si dichiarava: stesso digest.
    assenti = sorted(set(defaults.CONDIZIONALI) - defaults.CONDIZIONALI_NULLI_SE_OMESSI)
    assert risolvi(valida(_senza(completa, *assenti))).digest == DIGEST_SENZA_TIPO
    # Gli altri, omessi, sono nulli: il digest lo riporta, uno per uno. Quelli
    # della prevalenza qui hanno la condizione vera e non si possono omettere.
    falsi = [c for c in sorted(defaults.CONDIZIONALI_NULLI_SE_OMESSI)
             if not c.startswith("prev.")]
    for chiave in falsi:
        assert risolvi(valida(_senza(completa, chiave))).digest != DIGEST_SENZA_TIPO, chiave


def test_la_copia_per_un_nuovo_tentativo_conserva_i_condizionali_omessi():
    """
    **Obiettivo**: Verificare che la configurazione rivalidata per un nuovo
    tentativo, costruita da una copia in memoria, resti valida quando i
    condizionali sono omessi, con gli stessi valori.

    **Razionale scientifico e sistemistico**: Il nuovo tentativo rilegge la
    configurazione dalla sua forma serializzata, dove un parametro omesso e'
    nullo: se il nullo fosse un errore, ogni retry su un dataset senza
    controlli fallirebbe sulla configurazione invece che sul calcolo.
    """
    config = valida(_senza(_tutte_false(), *defaults.CONDIZIONALI))
    nuova = applica(config, "run.batch_size", 3)
    assert nuova.run.batch_size == 3
    assert nuova.ctrl.min_positives is None and nuova.phylo.max_seqs is None
    assert nuova.model_dump(exclude={"run"}) == config.model_dump(exclude={"run"})


# --------------------------------------------------------------------------- #
# Schema esportato e report                                                    #
# --------------------------------------------------------------------------- #


def test_lo_schema_esportato_riporta_la_specie_e_la_condizione():
    """
    **Obiettivo**: Verificare che ``schema`` dia a ogni condizionale la specie
    ``condizionale`` e ``richiesto_da`` con il parametro e la forma della
    condizione (e ``valori`` per le condizioni su un valore vero), e in piu'
    ``oppure`` con la seconda via per i parametri della decontaminazione; che
    non lo dia agli altri parametri e che la versione del formato sia
    invariata.

    **Razionale scientifico e sistemistico**: Un programma esterno che
    costruisce un modulo di configurazione deve sapere quali campi chiedere
    solo in presenza di un altro, anche quando le vie per renderlo pertinente
    sono due: la seconda e' un campo in piu', che non cambia il significato di
    quelli esistenti, quindi la versione del formato non cambia.
    """
    documento = interfaccia.schema()
    assert documento["versione_formato"] == 1
    per_nome = {p["nome"]: p for p in documento["parametri"]}
    for chiave, alternative in defaults.CONDIZIONALI.items():
        voce = per_nome[chiave]
        assert voce["specie"] == "condizionale" and voce["obbligatorio"] is False
        prima, *altre = [
            {"parametro": parametro, "condizione": forma,
             **({"valori": [True]} if forma == defaults.VERO else {})}
            for parametro, forma in alternative
        ]
        assert voce["richiesto_da"] == {**prima, **({"oppure": altre} if altre else {})}
    assert per_nome["decontam.threshold"]["richiesto_da"] == {
        "parametro": "ctrl.blank_values", "condizione": "non_vuoto",
        "oppure": [{"parametro": "ctrl.blank_override_values", "condizione": "non_vuoto"}]}
    assert "oppure" not in per_nome["phylo.max_seqs"]["richiesto_da"]
    # I parametri della modalita' appaiata sono condizionali della modalita'.
    for chiave in defaults.CONDIZIONALI_DELLE_COPPIE:
        assert per_nome[chiave]["specie"] == "condizionale"
        assert per_nome[chiave]["richiesto_da"] == {
            "parametro": "io.layout", "condizione": "appaiata", "valori": ["paired-end"]}
    assert {n for n, p in per_nome.items() if "richiesto_da" in p} == (
        set(defaults.CONDIZIONALI) | set(defaults.CONDIZIONALI_RESPINTI))
    # I parametri di una modalita' dicono quale valore la sceglie.
    assert per_nome["filter.primer_max_offset"]["richiesto_da"] == {
        "parametro": "filter.primer_removal", "condizione": "cercato", "valori": ["search"]}
    assert per_nome["tax.idtaxa_trained"]["richiesto_da"] == {
        "parametro": "tax.classifier", "condizione": "idtaxa", "valori": ["idtaxa"]}
    assert per_nome["tax.ref_fasta"]["richiesto_da"] == {
        "parametro": "tax.classifier", "condizione": "bayesiano", "valori": ["naive_bayes"]}
    assert per_nome["tax.classifier"]["specie"] == "modalita"
    assert per_nome["tax.classifier"]["valori_ammessi"] == ["naive_bayes", "idtaxa"]
    assert {p["specie"] for p in documento["parametri"]} == {
        "obbligatorio", "condizionale", "percorso", "standard", "unica_forma", "facoltativo",
        "strumento", "modalita"}
    json.dumps(documento)


def _origini(radice: Path) -> dict[str, str]:
    """L'origine di ogni parametro nella tabella del report."""
    percorso = Path(radice) / CARTELLA_REPORT / CARTELLA_TABELLE / "parametri.tsv"
    with open(percorso, encoding="utf-8", newline="") as file:
        return {r["parametro"]: r["origine"] for r in csv.DictReader(file, delimiter="\t")}


def test_il_report_distingue_non_pertinente_e_dichiarato_senza_effetto(tmp_path, capsys):
    """
    **Obiettivo**: Verificare, dopo la sola validazione di un dataset senza
    colonna del tipo, che il report dia l'origine ``non pertinente`` ai
    condizionali omessi e ``dichiarato, senza effetto`` a quelli dichiarati con
    la condizione falsa, che li conti nella sintesi della configurazione e che
    le due diciture non compaiano in altre sezioni.

    **Razionale scientifico e sistemistico**: Chi legge il report deve
    distinguere un valore scelto, uno che nessuno ha dovuto scegliere e uno
    scritto ma inutile: mostrarli tutti come dichiarati attribuirebbe
    all'analisi scelte che non ha fatto.
    """
    dataset = sintetico.genera(tmp_path / "dati", biologici=3, letture=40, con_tipo=False)
    omessi = sorted(defaults.CONDIZIONALI)
    for nome, tolti in (("omessi", omessi), ("dichiarati", [])):
        dati = _senza(
            sintetico.configurazione(dataset, tmp_path / nome, prev={"apply": False}), *tolti)
        file_config = tmp_path / f"{nome}.yaml"
        file_config.write_text(yaml.safe_dump(dati, sort_keys=False), encoding="utf-8")
        assert cli.main(["validate", "--config", str(file_config)]) == 0
        chiudi()
        pagina = genera(tmp_path / nome).read_text(encoding="utf-8")
        origini = _origini(tmp_path / nome)
        attesa = "non pertinente" if tolti else "dichiarato, senza effetto"
        assert {c: origini[c] for c in defaults.CONDIZIONALI} == dict.fromkeys(
            defaults.CONDIZIONALI, attesa)
        assert origini["ctrl.column"] == "dichiarato"
        assert origini["prev.apply"] == "dichiarato"
        quanti = len(defaults.CONDIZIONALI)
        sintesi = {"non pertinenti (condizionali omessi)": quanti if tolti else 0,
                   "dichiarati senza effetto (condizionali con la condizione falsa)":
                       0 if tolti else quanti}
        for voce, numero in sintesi.items():
            assert voce in pagina
            assert pagina.split(voce)[1].split("</")[1].endswith(f">{numero}"), voce
        fuori = pagina.split('id="configurazione"')[0]
        assert "senza effetto" not in fuori and "non pertinent" not in fuori
    capsys.readouterr()


# --------------------------------------------------------------------------- #
# La catena                                                                    #
# --------------------------------------------------------------------------- #


def _sha_degli_artefatti(run: ProjectRun) -> dict[str, str]:
    """Lo sha256 di ogni artefatto delle fasi, letto dai manifesti, senza la
    configurazione registrata (che dice quali parametri sono dichiarati)."""
    sha: dict[str, str] = {}
    for manifesto in sorted(Path(run.config.io.out_root).rglob("manifest_S*.json")):
        if manifesto.parent.name == Fase.CONFIG.value:
            continue
        for artefatto in json.loads(manifesto.read_text(encoding="utf-8"))["artefatti"]:
            cartella = manifesto.parent.relative_to(run.config.io.out_root).as_posix()
            sha[f"{cartella}/{artefatto['nome']}"] = artefatto["checksum"]
    return sha


#: I riepiloghi che riportano un condizionale anche con la condizione falsa,
#: con i campi che lo riportano: nulli se il parametro e' omesso.
RIPORTATI_NEI_RIEPILOGHI = {
    "11_controls/soglia.json": {"sensibilita": "katharoseq.target_sensitivity"},
    "11_controls/decontam_riepilogo.json": {
        "soglia": "decontam.threshold", "min_negativi": "decontam.min_blanks",
        "modalita_dichiarata": "decontam.mode", "max_frazione": "qc.max_frac_contaminant"},
    "12_final/intermedi/filtri_riepilogo.json": {
        "prevalenza.min_conteggio": "prev.min_count",
        # Derivato da prev.min_fraction: il minimo di campioni per variante.
        "prevalenza.min_campioni": "prev.min_fraction"},
}
#: I testi registrati che citano il valore di decontam.min_blanks: omesso, lo
#: tacciono.
TESTI_CHE_CITANO_IL_MINIMO = {
    "11_controls/decontam_riepilogo.json": {"esito"},
    "01_input_validation/gates.json": None,
}


def _foglie(valore: Any, percorso: str = "") -> dict[str, Any]:
    """Le foglie di un documento JSON, come ``percorso.puntato`` -> valore."""
    if isinstance(valore, dict):
        voci = valore.items()
    elif isinstance(valore, list):
        voci = enumerate(valore)
    else:
        return {percorso: valore}
    foglie: dict[str, Any] = {}
    for chiave, sotto in voci:
        foglie.update(_foglie(sotto, f"{percorso}.{chiave}" if percorso else str(chiave)))
    return foglie


def test_la_catena_arriva_in_fondo_senza_i_condizionali_non_pertinenti(bioc, tmp_path):
    """
    **Obiettivo**: Verificare su un dataset sintetico senza colonna del tipo,
    con filogenesi, prevalenza, riclassificazione e modulo derivato spenti, che
    la catena S0-S14 si concluda con ogni condizionale omesso. Rispetto alla
    catena con ogni condizionale dichiarato: conteggi, tassonomia, metadati e
    sequenze finali sono identici byte per byte; differiscono soltanto i
    riepiloghi che riportano un condizionale, e in essi soltanto i campi di
    quei parametri (nulli se omessi, con il valore dichiarato altrimenti) e i
    testi che citavano il minimo dei controlli negativi; il report mostra i
    campi nulli come ``non pertinente``.

    **Razionale scientifico e sistemistico**: E' la definizione di
    condizionale resa verificabile: con la condizione falsa il valore non entra
    in alcun calcolo. La vista delle fasi solleva un errore alla lettura di un
    condizionale senza valore di assenza, e gli script R si fermano se un
    valore nullo servisse a un calcolo: la catena conclusa dimostra che
    nessuna fase ne usa uno.
    """
    dataset = sintetico.genera(tmp_path / "dati", biologici=8, con_tipo=False)
    sha, pagine, valori = {}, {}, {}
    for nome, tolti in (("dichiarati", ()), ("omessi", tuple(defaults.CONDIZIONALI))):
        completa = sintetico.configurazione(
            dataset, tmp_path / nome, qc={"min_reads_mode": "katharoseq_if_available"},
            prev={"apply": False})
        valori[nome] = {f"{g}.{n}": v for g, c in completa.items() for n, v in c.items()}
        run = ProjectRun(valida(_senza(completa, *tolti)))
        esito = Esecutore(run).esegui()
        assert esito.conclusione is Conclusione.COMPLETATA, esito.punto
        chiudi()
        sha[nome] = _sha_degli_artefatti(run)
        pagine[nome] = genera(run.config.io.out_root).read_text(encoding="utf-8")
        assert "Risultato finale" in pagine[nome]
    assert set(sha["omessi"]) == set(sha["dichiarati"]) and len(sha["omessi"]) > 20
    diversi = {n for n in sha["omessi"] if sha["omessi"][n] != sha["dichiarati"][n]}
    # Differiscono i soli riepiloghi: anche l'oggetto finale e il suo elenco
    # di checksum sono identici.
    riepiloghi = set(RIPORTATI_NEI_RIEPILOGHI) | set(TESTI_CHE_CITANO_IL_MINIMO)
    assert diversi == riepiloghi, sorted(diversi ^ riepiloghi)
    for nome in ("conteggi.tsv", "tassonomia.tsv", "metadati.tsv", "sequenze.fasta"):
        assert (tmp_path / "omessi" / "12_final" / nome).read_bytes() == (
            tmp_path / "dichiarati" / "12_final" / nome).read_bytes(), nome

    # Campo per campo: nei riepiloghi cambiano i soli campi dei parametri omessi.
    for relativo in sorted(riepiloghi):
        foglie = {
            nome: _foglie(json.loads((tmp_path / nome / relativo).read_text(encoding="utf-8")))
            for nome in ("dichiarati", "omessi")}
        assert set(foglie["omessi"]) == set(foglie["dichiarati"]), relativo
        cambiati = {c for c in foglie["omessi"] if foglie["omessi"][c] != foglie["dichiarati"][c]}
        campi = RIPORTATI_NEI_RIEPILOGHI.get(relativo, {})
        testi = TESTI_CHE_CITANO_IL_MINIMO.get(relativo, set())
        if testi is None:
            # Gli avvisi dei gate: i soli testi diversi sono quelli che citavano il minimo.
            testi = cambiati
            assert cambiati, relativo
        assert cambiati == set(campi) | testi, (relativo, sorted(cambiati))
        for campo, parametro in campi.items():
            assert foglie["omessi"][campo] is None, (relativo, campo)
            if campo != "prevalenza.min_campioni":
                assert foglie["dichiarati"][campo] == valori["dichiarati"][parametro], campo
        for campo in testi:
            con, senza = foglie["dichiarati"][campo], foglie["omessi"][campo]
            assert "decontam.min_blanks" in con and "decontam.min_blanks" not in senza, campo
            assert "None" not in senza and "NULL" not in senza, campo
    minimo = json.loads((tmp_path / "dichiarati" / "12_final/intermedi/filtri_riepilogo.json")
                        .read_text(encoding="utf-8"))["prevalenza"]["min_campioni"]
    assert isinstance(minimo, int) and minimo >= 1

    # Il report, dopo la sezione della configurazione (che ha la sua colonna
    # dell'origine): i campi nulli dei riepiloghi sono mostrati come non
    # pertinenti, mai come un nullo.
    dopo = {nome: pagina.split('id="provenienza"')[1] for nome, pagina in pagine.items()}
    for testo in ("prevalence, non pertinente", "<dd>non pertinente</dd>",
                  "minimo di campioni per variante non pertinente"):
        assert testo in dopo["omessi"], testo
    assert "None" not in dopo["omessi"]
    assert "non pertinente" not in dopo["dichiarati"]
