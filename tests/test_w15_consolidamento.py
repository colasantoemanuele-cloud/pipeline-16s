r"""Suite di consolidamento delle fasi realizzate, da S0 a S7.

1. Inquadramento nel Piano Operativo
------------------------------------
Settimana 15 (W15), Fase F4 (consolidamento delle fasi di calcolo): nessuna
funzionalita' nuova, si chiude il lavoro aperto sulle fasi gia' realizzate.

2. Moduli sorgente coperti
--------------------------
* ``src/amplicon16s/steps/s00_validate.py`` (parametri dichiarati, G15 sulla
  configurazione completa), ``src/amplicon16s/gates/registry.py``
* ``src/amplicon16s/runner/executor.py``, ``src/amplicon16s/runner/project.py``
  (ripresa, valutazione dello stato)
* ``src/amplicon16s/errors/catalog.py`` (copertura dei codici), e per i codici
  che nessun test provocava ``src/amplicon16s/gates/g01_g15.py`` (E-G15-01,
  E-G15-99), ``src/amplicon16s/steps/s02_filter.py`` (E-S2-01) e
  ``src/amplicon16s/steps/s06_chimera.py`` (E-S6-01, E-S6-02)

3. Cosa valuta questo file
--------------------------
- che ogni codice del catalogo sollevabile dal codice realizzato abbia un
  test che lo provoca e ne verifica la categoria di gestione;
- che S0 non dipenda dai parametri che legge solo per G15: cambiare
  ``qc.warn_frac_chimeric`` rifa' soltanto S6 e S7, e una configurazione
  incoerente e' comunque respinta all'avvio;
- la ripresa da ogni fase, da S1 a S7: alterato o rimosso un suo artefatto,
  si rifanno esattamente quella fase e quelle che ne dipendono, con gli
  stessi byte.

4. Comandi Bash e scenari di esecuzione
---------------------------------------
    ``<immagine>`` e' l'immagine del container della pipeline; quella corrente
    e' indicata in ``test.txt``, sezione 1.3.

    1. Modalita locale standard (senza Bioconductor R):
       pytest tests/test_w15_consolidamento.py -v
       Risultato atteso: 16 test (5 passed in Python puro, 11 skipped per
       assenza di dada2 in R locale) in ~5s.

    2. Modalita container Docker standard (subset ridotto con Bioconductor):
       docker run --rm \
         -e PYTHONPATH=/app/src \
         -e AMPLICON16S_R_DIR=/app/R \
         -v "$(pwd)":/app \
         -w /app \
         <immagine> \
         pytest -o cache_dir=/tmp/.pytest_cache tests/test_w15_consolidamento.py -v
       Risultato atteso: 16 passed in ~4 minuti (le riprese rifanno sulla
       versione ridotta fino a sei fasi ciascuna).

    3. Modalita container Docker completa (con dati reali OSD-734; nessun test
       del modulo li legge, il risultato e' lo stesso: 16 passed):
       docker run --rm \
         --memory=24g \
         --memory-swap=24g \
         -e PYTHONPATH=/app/src \
         -e AMPLICON16S_R_DIR=/app/R \
         -e AMPLICON16S_CONFIG_DATI_REALI="$HOME/ASI/config_osd734.yaml" \
         -v "$(pwd)":/app \
         -v "$HOME/ASI":"$HOME/ASI" \
         -w /app \
         <immagine> \
         pytest -o cache_dir=/tmp/.pytest_cache tests/test_w15_consolidamento.py -v

5. Risultato atteso
-------------------
16 test totali (10 funzioni di test, di cui 1 parametrizzata su 7 fasi):
- 5 passed, 11 skipped in ambiente locale privo di Bioconductor;
- 16 passed nel container CI e con i dati reali OSD-734.

6. Razionale scientifico e sistemistico
---------------------------------------
- Un codice mai provocato da un test e' un ramo di gestione mai visto
  funzionare: la sua categoria (arresto, retry, degradazione) e' una promessa
  non verificata.
- La ripresa e' cio' che rende sostenibile un'esecuzione di un'ora: deve
  rifare tutto cio' che dipende da un artefatto perso, e nient'altro.
"""

from __future__ import annotations

import hashlib
import os
from pathlib import Path

import pytest
from conftest import copia_esecuzione
from sottoinsieme import config_ridotta, motivo_pacchetti_r_assenti

from amplicon16s.gates.g01_g15 import ErroreGate
from amplicon16s.io_layer.artifacts import Fase
from amplicon16s.logging.logger import chiudi
from amplicon16s.runner.executor import Conclusione, Esecutore
from amplicon16s.runner.graph import Passo
from amplicon16s.runner.project import StatoPasso, passi_realizzati
from amplicon16s.steps.s00_validate import esegui_s0


@pytest.fixture(autouse=True)
def uscite_pulite():
    """Chiude le uscite del log prima e dopo ogni test, perché nessun handler resti
    aperto sulla cartella temporanea.
    """
    chiudi()
    yield
    chiudi()


_MOTIVO_ASSENTI = motivo_pacchetti_r_assenti("dada2", "ggplot2", "ShortRead", "jsonlite")


@pytest.fixture
def dada2():
    """Richiede R con dada2: salta senza, ma in CI fallisce."""
    if _MOTIVO_ASSENTI is not None:
        if os.environ.get("AMPLICON16S_RICHIEDI_BIOC") == "1":
            pytest.fail(f"Bioconductor e' richiesto in questo ambiente: {_MOTIVO_ASSENTI}")
        pytest.skip(_MOTIVO_ASSENTI)


def _impronte(cartella: Path) -> dict[str, str]:
    """Calcola il dizionario MD5 degli artefatti nella cartella, escludendo i file temporanei rbridge_ e i manifesti."""
    return {
        p.name: hashlib.md5(p.read_bytes()).hexdigest()
        for p in sorted(cartella.iterdir())
        if not p.name.startswith(("rbridge_", "manifest"))
    }


_FASI = (Fase.QC_PROFILES, Fase.FILTERED, Fase.ERROR_MODELS, Fase.ASV_INFERENCE,
         Fase.SEQTAB, Fase.CHIMERA)


# --------------------------------------------------------------------------- #
# S0 e i parametri di G15                                                      #
# --------------------------------------------------------------------------- #


def test_una_configurazione_incoerente_e_respinta_da_s0_da_sola(tmp_path):
    """
    **Obiettivo**: Verificare che S0 eseguita da sola (``esegui_s0``, il
    comando ``validate``) respinga con G15 una ``retry.whitelist`` che
    ammette al retry un codice di revisione umana, anche se S0 non dichiara
    ``retry.whitelist`` fra i propri parametri.

    **Razionale scientifico e sistemistico**: G15 si esegue sulla configurazione
    completa, prima che la fase la veda ristretta: restringere l'impronta di
    S0 non deve togliere il controllo di coerenza.
    """
    config = config_ridotta(tmp_path, retry={"whitelist": ["E-S2-03", "E-S6-01"]})
    with pytest.raises(ErroreGate) as info:
        esegui_s0(config)
    assert info.value.gate == "G15"
    assert [v.codice for v in info.value.violazioni] == ["E-G15-09"]
    assert info.value.categoria.value == "revisione_umana"


def test_cambiare_qc_warn_frac_chimeric_rifa_solo_s6_e_s7(dada2, catena_calcolata, tmp_path):
    """
    **Obiettivo**: Verificare che, con la catena S0-S7 conclusa, cambiare
    ``qc.warn_frac_chimeric`` lasci concluse S0-S5 e rifaccia soltanto S6 e
    la fase che ne dipende, S7.

    **Razionale scientifico e sistemistico**: S0 leggeva l'intero gruppo qc per
    G15, e una soglia di S6 rifaceva la catena intera, S4 compresa (25 minuti
    sul dataset completo). G15 si ripete a ogni avvio: i parametri che legge
    solo lui non determinano i risultati di S0.
    """
    run = copia_esecuzione(catena_calcolata, tmp_path, qc={"warn_frac_chimeric": 0.3})
    situazione = run.valuta().situazioni
    for passo in (Passo.S0, Passo.S1, Passo.S2, Passo.S3, Passo.S4, Passo.S5):
        assert situazione[passo].stato is StatoPasso.COMPLETATA, passo
    assert situazione[Passo.S6].motivo == "configurazione cambiata"
    esito = Esecutore(run, fino_a=Passo.S7).esegui()
    assert esito.conclusione is Conclusione.COMPLETATA
    assert [r.passo for r in esito.eseguite] == [Passo.S6, Passo.S7]


def test_una_configurazione_incoerente_e_respinta_all_avvio(dada2, catena_calcolata, tmp_path):
    """
    **Obiettivo**: Verificare che, con la catena conclusa, una configurazione
    che G15 respinge fermi l'esecutore all'avvio, prima di qualunque fase, e
    che S0 resti conclusa: la whitelist non e' nella sua impronta.

    **Razionale scientifico e sistemistico**: Togliere i parametri di G15
    dall'impronta di S0 e' sicuro solo perche' G15 si ripete a ogni avvio.
    """
    run = copia_esecuzione(
        catena_calcolata, tmp_path, retry={"whitelist": ["E-S2-03", "E-S6-01"]}
    )
    assert run.valuta().situazioni[Passo.S0].stato is StatoPasso.COMPLETATA
    manifesto = run.albero.percorso_manifesto_passo(Passo.S6, Fase.CHIMERA).read_bytes()
    with pytest.raises(ErroreGate) as info:
        Esecutore(run, fino_a=Passo.S7).esegui()
    assert info.value.gate == "G15"
    assert [v.codice for v in info.value.violazioni] == ["E-G15-09"]
    assert run.albero.percorso_manifesto_passo(Passo.S6, Fase.CHIMERA).read_bytes() == manifesto


# --------------------------------------------------------------------------- #
# Ripresa da ogni fase                                                         #
# --------------------------------------------------------------------------- #


def _altera(percorso: Path) -> None:
    """Manomette deterministicamente un file aggiungendo un singolo byte newline per collaudare il ricalcolo selettivo."""
    percorso.write_bytes(percorso.read_bytes() + b"\n")


_RIPRESE = [
    # fase, cartella, artefatto, azione, fasi che si rifanno
    (Passo.S2, Fase.FILTERED, "letture_filtrate.tsv", "rimuovi",
     [Passo.S2, Passo.S3, Passo.S4, Passo.S5, Passo.S6, Passo.S7]),
    (Passo.S3, Fase.ERROR_MODELS, "corrispondenza.tsv", "altera",
     [Passo.S3, Passo.S4, Passo.S5, Passo.S6, Passo.S7]),
    (Passo.S4, Fase.ASV_INFERENCE, "varianti_per_campione.rds", "rimuovi",
     [Passo.S4, Passo.S5, Passo.S6, Passo.S7]),
    (Passo.S5, Fase.SEQTAB, "tabella.json", "altera", [Passo.S5, Passo.S6, Passo.S7]),
    (Passo.S6, Fase.CHIMERA, "chimere.tsv", "rimuovi", [Passo.S6, Passo.S7]),
    (Passo.S7, Fase.CHIMERA, "lunghezze.tsv", "altera", [Passo.S7]),
]


@pytest.mark.parametrize(
    ("passo", "cartella", "nome", "azione", "rifatte"),
    _RIPRESE,
    ids=[f"{r[0]}-{r[3]}-{r[2]}" for r in _RIPRESE],
)
def test_la_ripresa_rifa_la_fase_e_solo_quelle_che_ne_dipendono(
    dada2, catena_calcolata, tmp_path, passo, cartella, nome, azione, rifatte
):
    """
    **Obiettivo**: Verificare, per ogni fase da S2 a S7, che rimosso o
    alterato un suo artefatto la valutazione la dia da rifare, che la ripresa
    riesegua esattamente quella fase e le fasi che ne dipendono, nell'ordine
    del grafo, e che gli artefatti rifatti siano identici agli originali.

    **Razionale scientifico e sistemistico**: Una ripresa che rifacesse meno
    consegnerebbe risultati calcolati su un artefatto che non c'e' piu'; una
    che rifacesse di piu' costerebbe ore sul dataset completo. S2 alimenta S3,
    S4 e S7; S7 non alimenta nulla. S1 non e' fra gli antenati di S7, e
    l'esecuzione fino a S7 non la rifarebbe: la sua ripresa e' verificata con i
    test di S1 (``tests/test_w11_s01_profile.py``).
    """
    originale, _ = catena_calcolata
    run = copia_esecuzione(catena_calcolata, tmp_path)
    artefatto = run.albero.cartella(cartella) / nome
    if azione == "rimuovi":
        artefatto.unlink()
    else:
        _altera(artefatto)

    situazione = run.valuta().situazioni
    assert situazione[passo].stato is StatoPasso.DA_ESEGUIRE
    assert nome in situazione[passo].motivo
    for altro in (Passo.S0, Passo.S1, Passo.S2, Passo.S3, Passo.S4, Passo.S5, Passo.S6, Passo.S7):
        if altro not in rifatte:
            assert situazione[altro].stato is StatoPasso.COMPLETATA, altro

    esito = Esecutore(run, fino_a=Passo.S7).esegui()
    assert esito.conclusione is Conclusione.COMPLETATA
    assert [r.passo for r in esito.eseguite] == rifatte
    for fase in _FASI:
        assert _impronte(run.albero.cartella(fase)) == _impronte(originale.albero.cartella(fase)), fase


# --------------------------------------------------------------------------- #
# E-S6-01 ed E-S6-02 dalla fase vera                                           #
# --------------------------------------------------------------------------- #


def test_s6_oltre_l_avviso_registra_la_degradazione_e_prosegue(dada2, catena_calcolata, tmp_path):
    """
    **Obiettivo**: Verificare che S6, con ``qc.warn_frac_chimeric`` sotto la
    frazione chimerica della versione ridotta (0,0008), si concluda
    registrando E-S6-02 nel proprio manifesto, e che S7 prosegua.

    **Razionale scientifico e sistemistico**: E-S6-02 e' una degradazione
    automatica: non ferma nulla, ma deve restare agli atti della fase.
    """
    run = copia_esecuzione(
        catena_calcolata, tmp_path, qc={"warn_frac_chimeric": 0.0001, "stop_frac_chimeric": 0.5}
    )
    esito = Esecutore(run, fino_a=Passo.S7).esegui()
    assert esito.conclusione is Conclusione.COMPLETATA
    assert [r.passo for r in esito.eseguite] == [Passo.S6, Passo.S7]
    manifesto = run.albero.manifesto_passo(Passo.S6, Fase.CHIMERA)
    (degradazione,) = manifesto.degradazioni
    assert degradazione["codice"] == "E-S6-02"
    assert degradazione["categoria"] == "degradazione_automatica"


def test_s6_oltre_l_arresto_si_ferma_per_la_revisione(dada2, catena_calcolata, tmp_path):
    """
    **Obiettivo**: Verificare che S6, con ``qc.stop_frac_chimeric`` sotto la
    frazione chimerica della versione ridotta, si fermi con E-S6-01 senza
    tentativi, lasciando le misure in ``chimere.json`` e nessun manifesto.

    **Razionale scientifico e sistemistico**: E-S6-01 e' di revisione umana: un
    eccesso di chimere si capisce, non si ritenta.
    """
    run = copia_esecuzione(
        catena_calcolata, tmp_path, qc={"warn_frac_chimeric": 0.0001, "stop_frac_chimeric": 0.0002}
    )
    esito = Esecutore(run, fino_a=Passo.S7).esegui()
    assert esito.conclusione is Conclusione.ARRESTATA
    assert esito.punto.passo is Passo.S6
    assert esito.punto.codice == "E-S6-01"
    assert esito.punto.categoria == "revisione_umana"
    assert esito.punto.tentativi == 1
    assert run.albero.manifesto_passo(Passo.S6, Fase.CHIMERA) is None
    assert (run.albero.cartella(Fase.CHIMERA) / "chimere.json").exists()


# --------------------------------------------------------------------------- #
# Codici senza un test che li provocasse dal codice realizzato                 #
# --------------------------------------------------------------------------- #


def test_g15_sorveglia_la_regola_di_derivazione(tmp_path):
    """
    **Obiettivo**: Verificare che G15 respinga con E-G15-01, di revisione
    umana, una configurazione risolta in cui ``filter.minLen`` supera
    ``filter.truncLen``.

    **Razionale scientifico e sistemistico**: ``filter.minLen`` discende da
    ``filter.truncLen``, quindi nessuna configurazione puo' violare la
    disuguaglianza: il controllo sorveglia la regola di derivazione. Si
    provoca costruendo a mano i derivati incoerenti, come li produrrebbe una
    regola cambiata per errore.
    """
    from dataclasses import replace

    from amplicon16s.config.resolve import risolvi
    from amplicon16s.errors.catalog import Categoria
    from amplicon16s.gates.g01_g15 import _controlla_coerenza

    risolta = risolvi(config_ridotta(tmp_path))
    derivati = replace(risolta.derivati, filter_minLen=risolta.config.filter.truncLen + 1)
    (violazione,) = _controlla_coerenza(replace(risolta, derivati=derivati))
    assert violazione.codice == "E-G15-01"
    assert violazione.voce.categoria is Categoria.REVISIONE_UMANA


def test_un_problema_di_schema_non_attribuito_e_e_g15_99(tmp_path):
    """
    **Obiettivo**: Verificare che una chiave sconosciuta, che lo schema
    respinge senza che un controllo specifico la reclami, fermi G15 con
    E-G15-99, di revisione umana.

    **Razionale scientifico e sistemistico**: Le chiavi sconosciute sono
    respinte perche' un errore di battitura non diventi un parametro ignorato
    in silenzio; il codice generico garantisce che nessun rifiuto resti senza
    codice.
    """
    from sottoinsieme import dati_config

    from amplicon16s.errors.catalog import Categoria
    from amplicon16s.gates.g01_g15 import esegui_g15

    dati = dati_config(tmp_path)
    dati["filter"]["truncLenn"] = 137
    with pytest.raises(ErroreGate) as info:
        esegui_g15(dati)
    assert [v.codice for v in info.value.violazioni] == ["E-G15-99"]
    assert info.value.categoria is Categoria.REVISIONE_UMANA


def test_s2_ferma_su_un_biologico_azzerato_e_non_su_un_negativo(tmp_path):
    """
    **Obiettivo**: Verificare che il controllo di S2 sollevi E-S2-01, di
    revisione umana, per un campione biologico azzerato dal filtro, e non per
    un controllo negativo azzerato.

    **Razionale scientifico e sistemistico**: Il test sulla fase vera con un
    filtro troppo severo accetta E-S2-01 o E-S2-02, perche' scattano insieme:
    questo li separa. Un bianco azzerato e' un bianco pulito.
    """
    from amplicon16s.errors.exceptions import ErroreRevisioneUmana
    from amplicon16s.metadata.models import ClasseCampione
    from amplicon16s.steps.s02_filter import controlla_filtro

    qc = config_ridotta(tmp_path).qc
    classi = {"B1": ClasseCampione.BIOLOGICO, "B2": ClasseCampione.BIOLOGICO,
              "N1": ClasseCampione.CONTROLLO_NEGATIVO}
    ingresso = {"B1": 1000, "B2": 1000, "N1": 100}
    metriche = controlla_filtro(ingresso, {"B1": 990, "B2": 990, "N1": 0}, classi, qc)
    assert metriche["azzerati"]["controllo_negativo"] == ["N1"]
    with pytest.raises(ErroreRevisioneUmana) as info:
        controlla_filtro(ingresso, {"B1": 990, "B2": 0, "N1": 90}, classi, qc)
    assert info.value.codice == "E-S2-01"
    assert info.value.contesto["campioni"] == ["B2"]


# --------------------------------------------------------------------------- #
# Copertura dei codici del catalogo                                            #
# --------------------------------------------------------------------------- #

#: Le fasi del catalogo realizzate come codice, ricavate dalle fasi registrate
#: nell'esecuzione (non scritte a mano: una fase nuova entra da sola), con G15,
#: il grafo e il ponte. I codici di una fase non ancora realizzata restano fuori.
_REALIZZATE = frozenset({str(p) for p in passi_realizzati()} | {"G15", "GRAFO", "R"})

#: Per ogni codice sollevabile dal codice realizzato, il test che lo provoca
#: attraverso quel codice (non costruendo l'errore a mano) e ne verifica la
#: categoria di gestione. L'elenco e' stato ricavato registrando, durante
#: l'intera suite nel container, ogni codice creato con il test e il punto del
#: codice di produzione che lo ha sollevato.
COPERTURA: dict[str, str] = {
    "E-G15-01": "test_w15_consolidamento.py::test_g15_sorveglia_la_regola_di_derivazione",
    "E-G15-02": "test_w03_w04_config_schema.py::test_g15_respinge_le_configurazioni_incoerenti",
    "E-G15-03": "test_w03_w04_config_schema.py::test_g15_respinge_le_configurazioni_incoerenti",
    "E-G15-04": "test_w03_w04_config_schema.py::test_g15_respinge_le_configurazioni_incoerenti",
    "E-G15-05": "test_w03_w04_config_schema.py::test_g15_respinge_le_configurazioni_incoerenti",
    "E-G15-06": "test_w03_w04_config_schema.py::test_g15_respinge_le_configurazioni_incoerenti",
    "E-G15-07": "test_w03_w04_config_schema.py::test_g15_respinge_le_configurazioni_incoerenti",
    "E-G15-08": "test_w03_w04_config_schema.py::test_derivato_impostato_a_mano_e_un_errore",
    "E-G15-09": "test_w09_w10_retry_policy.py::test_una_whitelist_incoerente_e_respinta_prima_di_ogni_fase",
    "E-G15-99": "test_w15_consolidamento.py::test_un_problema_di_schema_non_attribuito_e_e_g15_99",
    "E-GRAFO-01": "test_w08_w09_graph_resume.py::test_nessuna_fase_gira_prima_delle_sue_dipendenze",
    "E-R-01": "test_w08_rbridge.py::test_senza_interprete_l_errore_e_del_catalogo",
    "E-R-02": "test_w08_rbridge.py::test_processo_bloccato_ucciso_allo_scadere_del_tempo",
    "E-R-03": "test_w08_rbridge.py::test_errore_r_non_catalogato",
    "E-R-04": "test_w08_rbridge.py::test_memoria_in_una_fase_senza_codice",
    "E-S0-01": "test_w07_gates.py::test_g01_fallisce_se_non_ci_sono_file_di_letture",
    "E-S0-02": "test_w07_gates.py::test_g02_fallisce_su_una_colonna_assente",
    "E-S0-03": "test_w06_join_restricted.py::test_g03_fallisce_se_un_campione_dell_assay_manca_nello_studio",
    "E-S0-04": "test_w06_crosswalk.py::test_g04_fallisce_se_l_accession_e_ambiguo",
    "E-S0-05": "test_w06_crosswalk.py::test_g05_fallisce_su_accession_duplicati_fra_i_file",
    "E-S0-06": "test_w06_crosswalk.py::test_g06_fallisce_su_un_file_senza_riga_nell_assay",
    "E-S0-07": "test_w07_gates.py::test_g07_fallisce_su_un_file_di_lettura_inversa",
    "E-S0-08": "test_w07_gates.py::test_g08_fallisce_se_manca_la_colonna_del_lotto",
    "E-S0-09": "test_w07_gates.py::test_g09_fallisce_se_il_troncamento_supera_le_letture",
    "E-S0-10": "test_w07_gates.py::test_g10_fallisce_se_il_primer_e_in_testa",
    "E-S0-11": "test_w06_crosswalk.py::test_g11_fallisce_su_un_materiale_non_mappato",
    "E-S0-12": "test_w07_gates.py::test_g12_fallisce_su_un_checksum_diverso",
    "E-S0-13": "test_w07_gates.py::test_g13_fallisce_se_la_qualita_ha_lunghezza_diversa",
    "E-S0-14": "test_w07_gates.py::test_g14_fallisce_se_i_thread_eccedono_le_cpu",
    "E-S0-15": "test_w07_gates.py::test_g08_avvisa_su_una_piastra_con_pochi_controlli_negativi",
    "E-S1-01": "test_w11_s01_profile.py::test_con_truncLen_120_lo_scarto_e_registrato_da_s1_non_da_s0",
    "E-S1-02": "test_w11_s01_profile.py::test_s1_ferma_se_il_troncamento_supera_il_minimo_vero",
    "E-S2-01": "test_w15_consolidamento.py::test_s2_ferma_su_un_biologico_azzerato_e_non_su_un_negativo",
    "E-S2-02": "test_w12_s02_filter.py::test_una_perdita_media_oltre_il_30_per_cento_ferma",
    "E-S2-03": "test_w12_s02_filter.py::test_un_archivio_corrotto_si_ferma_dopo_i_tentativi",
    "E-S3-01": "test_w13_s03_batch.py::test_la_mancata_convergenza_con_basi_in_piu_si_ritenta",
    "E-S3-02": "test_w13_s03_batch.py::test_un_campione_senza_corsa_con_la_colonna_attiva_ferma",
    "E-S4-02": "test_w14_s04_s07_denoising.py::test_con_memoria_ridotta_scatta_e_s4_02_e_il_retry_dimezza_il_lotto",
    "E-S5-01": "test_w14_s04_s07_denoising.py::test_oltre_qc_max_asv_count_s5_si_ferma_prima_della_tabella",
    "E-S6-01": "test_w15_consolidamento.py::test_s6_oltre_l_arresto_si_ferma_per_la_revisione",
    "E-S6-02": "test_w15_consolidamento.py::test_s6_oltre_l_avviso_registra_la_degradazione_e_prosegue",
    "E-S7-01": "test_w16_recupero.py::test_s7_senza_varianti_si_ferma_con_e_s7_01",
    "E-S8-02": "test_w17_tassonomia.py::test_e_s8_02_ferma_con_una_copertura_insufficiente",
    "E-S10-01": "test_w18_oggetto.py::test_s10_intercetta_una_tabella_di_s7_quadrata_trasposta",
    "E-S10-02": "test_w18_oggetto.py::test_una_colonna_ambigua_ferma_con_e_s10_02",
    "E-S11-02": "test_w19_controlli.py::test_il_ripiego_e_registrato_con_il_motivo",
    "E-S11-03": "test_w19_controlli.py::test_controlli_non_conformi_fanno_scattare_e_s11_03_o_l_avviso",
    "E-S11-04": "test_w19_controlli.py::test_controlli_non_conformi_fanno_scattare_e_s11_03_o_l_avviso",
    "E-S12-02": "test_w20_decontam.py::test_oltre_la_frazione_massima_la_modalita_dichiarata_si_ferma",
    "E-S13-01": "test_w20_decontam.py::test_il_filtro_prima_della_decontaminazione_fallisce_con_e_s13_01",
    "E-S13-02": "test_w21_finale.py::test_un_campione_svuotato_dal_filtro_di_prevalenza_ferma_con_e_s13_02",
    "E-S13-03": "test_w21_finale.py::test_un_campione_svuotato_dal_filtro_tassonomico_esce_con_e_s13_03",
    "E-S14-01": "test_w21_finale.py::test_una_frazione_trattenuta_sotto_la_soglia_ferma_con_e_s14_01",
}

_RADICE = Path(__file__).resolve().parents[1]


def _sollevabili() -> set[str]:
    """I codici che compaiono nel codice realizzato, delle fasi realizzate."""
    import re

    from amplicon16s.errors.catalog import voce

    sorgenti = [
        p for p in (_RADICE / "src" / "amplicon16s").rglob("*.py") if p.name != "catalog.py"
    ] + list((_RADICE / "R").rglob("*.R"))
    codici = set()
    for percorso in sorgenti:
        codici |= set(re.findall(r'"(E-[A-Z0-9]+-\d{2})"', percorso.read_text(encoding="utf-8")))
    return {c for c in codici if voce(c).fase in _REALIZZATE}


def test_ogni_codice_sollevabile_ha_un_test_che_lo_provoca():
    """
    **Obiettivo**: Verificare che l'insieme dei codici sollevabili dal codice
    realizzato, ricavato dal sorgente Python e R, coincida con quello della
    tabella :data:`COPERTURA`, e che ogni test indicato esista.

    **Razionale scientifico e sistemistico**: Un codice nuovo in una fase
    realizzata fa fallire questo test finche' non gli si associa un test che
    lo provochi; un codice che il codice non solleva piu' esce dalla tabella.
    La categoria di gestione e' verificata dai test indicati.
    """
    assert _sollevabili() == set(COPERTURA)
    for codice, riferimento in COPERTURA.items():
        modulo, nome = riferimento.split("::")
        testo = (_RADICE / "tests" / modulo).read_text(encoding="utf-8")
        assert f"def {nome}(" in testo, f"{codice}: {riferimento} non esiste"
