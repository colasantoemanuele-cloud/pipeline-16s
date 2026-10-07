r"""Suite di test per le fasi S4-S7 (il nucleo del denoising).

1. Inquadramento nel Piano Operativo
------------------------------------
Settimana 14 (W14), Fase F4 (Fasi di calcolo R/Bioconductor: S04 Denoising
DADA2 in ``05_asv_inference/``, S05 Tabella delle sequenze in ``06_seqtab/``,
S06 Rimozione delle chimere e S07 Filtro di lunghezza ASV in ``07_chimera/``).

2. Moduli sorgente coperti
--------------------------
* ``src/amplicon16s/steps/s04_dada.py``, ``R/04_dada.R``
* ``src/amplicon16s/steps/s05_seqtab.py``, ``R/05_seqtab.R``
* ``src/amplicon16s/steps/s06_chimera.py``, ``R/06_chimera.R``
* ``src/amplicon16s/steps/s07_asv_length.py``, ``R/07_asv_length.R``
* ``R/lib/risorse.R`` (``memoria_di_picco``), ``R/lib/errors.R`` (``richiedi_pacchetti``)
* ``src/amplicon16s/runner/retry.py``, ``src/amplicon16s/runner/tracciamento.py``
* ``src/amplicon16s/rbridge/runner.py``

3. Cosa valuta questo file
--------------------------
- inferenza delle varianti ASV con ``dada2::dada`` e il pseudo-pooling a lotti
  (``run.batch_size``): equivalenza esatta con ``dada(pool = "pseudo")`` in una
  sola chiamata su campioni della stessa corsa, stima di seconda passata con
  ``loessErrfun`` e invarianza byte per byte rispetto a ``run.batch_size``;
- aggregazione della matrice campioni per ASV (``dada2::makeSequenceTable``) e
  controllo preventivo ``qc.max_asv_count`` prima dell'allocazione (``E-S5-01``);
- rimozione de novo delle chimere bimeriche (``dada2::removeBimeraDenovo``),
  mancata duplicazione di ``tabella.rds`` di S5, calcolo delle frazioni
  chimeriche per classe (su letture e su varianti) e applicazione dei controlli
  ``E-S6-01`` ed ``E-S6-02`` alle sole classi controllate (esclusione dei
  controlli negativi);
- filtraggio posizionale della lunghezza dell'amplicone V4 (``asv.len_min`` e
  ``asv.len_max``) su varianti sintetiche e invarianza sui dati troncati a
  lunghezza fissa;
- gestione dei retry (dimezzamento di ``run.batch_size`` su ``E-S4-02``) e
  soppressioni motivate (``RITENTARE_INUTILE`` per ``E-S5-01`` e per ``E-S4-02``
  con ``dada.pool`` vero);
- ricomposizione completa del tracciamento delle letture da S1 a S7 (la
  riproducibilita' byte per byte della catena e' verificata da
  ``tests/test_w16_recupero.py`` su una seconda esecuzione completa S0-S7).

4. Comandi Bash e scenari di esecuzione
---------------------------------------
    ``<immagine>`` e' l'immagine del container della pipeline; quella corrente
    e' indicata in ``test.txt``, sezione 1.3.

    1. Modalita locale standard (senza Bioconductor R):
       pytest tests/test_w14_s04_s07_denoising.py -v
       Risultato atteso: 16 test (4 passed in Python puro, 12 skipped per
       assenza di dada2/ShortRead in R locale e dei dati reali in ~0.60s).

    2. Modalita container Docker standard (subset ridotto con Bioconductor):
       docker run --rm \
         -e PYTHONPATH=/app/src \
         -e AMPLICON16S_R_DIR=/app/R \
         -v "$(pwd)":/app \
         -w /app \
         <immagine> \
         pytest -o cache_dir=/tmp/.pytest_cache tests/test_w14_s04_s07_denoising.py -v
       Risultato atteso: 15 passed, 1 skipped in ~80s (resta saltato solo il
       test sui 960 campioni reali OSD-734).

    3. Modalita container Docker completa (con dati reali OSD-734):
       docker run --rm \
         --memory=24g \
         --memory-swap=24g \
         -e PYTHONPATH=/app/src \
         -e AMPLICON16S_R_DIR=/app/R \
         -e AMPLICON16S_CONFIG_DATI_REALI=$HOME/ASI/config_osd734.yaml \
         -v "$(pwd)":/app \
         -v $HOME/ASI:$HOME/ASI \
         -w /app \
         <immagine> \
         pytest -o cache_dir=/tmp/.pytest_cache tests/test_w14_s04_s07_denoising.py -v
       Risultato atteso: 16 passed (100% verde).

5. Risultato atteso
-------------------
16 test totali (15 funzioni di test, di cui 1 parametrizzata su 2 casi):
- 4 passed, 12 skipped in ambiente locale privo di Bioconductor (~0.60s);
- 15 passed, 1 skipped nel container CI sul sottoinsieme ridotto (~80s);
- 16 passed nel container Docker con il dataset completo OSD-734.

6. Razionale scientifico e sistemistico
---------------------------------------
- ``dada(pool = "pseudo")`` in una sola chiamata trattiene in memoria gli
  oggetti di tutti i campioni; l'esecuzione in due passate esplicite a lotti di
  ``run.batch_size`` riproduce il medesimo risultato biologico e gli stessi byte
  vincolando il picco di memoria alla dimensione del singolo lotto. Poiche' il
  lotto non altera il risultato, l'azione correttiva di ``E-S4-02`` che dimezza
  ``run.batch_size`` preserva integralmente le assunzioni metodologiche.
- In S6 la frazione chimerica si valuta sulle letture (e non sul conteggio
  grezzo delle varianti, dove le chimere rare peserebbero quanto le ASV
  dominanti) ed esclude i controlli negativi a bassa biomassa.
"""

from __future__ import annotations

import csv
import hashlib
import json
import os
import shutil
import subprocess
from pathlib import Path

import pytest
from sottoinsieme import config_ridotta, motivo_pacchetti_r_assenti

from amplicon16s.config.schema import valida
from amplicon16s.errors.exceptions import ErroreRevisioneUmana, ErroreRitentabileConRevisione
from amplicon16s.io_layer.artifacts import AlberoOutput, Fase
from amplicon16s.logging.logger import chiudi
from amplicon16s.rbridge.runner import cartella_r, esegui_script, trova_rscript
from amplicon16s.runner.executor import Conclusione, Esecutore
from amplicon16s.runner.graph import Passo
from amplicon16s.runner.project import ProjectRun, StatoPasso
from amplicon16s.runner.retry import Motivo, PoliticaRetry
from amplicon16s.runner.tracciamento import ricomponi
from amplicon16s.steps.s04_dada import InferenzaVarianti
from amplicon16s.steps.s05_seqtab import TabellaSequenze
from amplicon16s.steps.s06_chimera import controlla_chimere, misura_chimere


@pytest.fixture(autouse=True)
def uscite_pulite():
    """Chiude le uscite del log prima e dopo ogni test, perché nessun handler resti
    aperto sulla cartella temporanea.
    """
    chiudi()
    yield
    chiudi()


# --------------------------------------------------------------------------- #
# Retry: E-S4-02 dimezza il lotto, E-S5-01 e' dichiarato inutile               #
# --------------------------------------------------------------------------- #


def test_e_s4_02_si_ritenta_con_il_lotto_dimezzato(tmp_path):
    """
    **Obiettivo**: Verificare che la politica ritenti E-S4-02 con
    ``run.batch_size`` dimezzato, e smetta quando il lotto e' gia' 1.

    **Razionale scientifico e sistemistico**: Il lotto decide la memoria
    dell'inferenza e non il suo risultato: ridurlo non cambia alcuna
    assunzione metodologica.
    """
    config = config_ridotta(tmp_path, run={"batch_size": 24})
    politica = PoliticaRetry.da_config(config)
    azione = InferenzaVarianti.aggiustamenti["E-S4-02"]
    decisione = politica.decidi("E-S4-02", 1, azione, config)
    assert decisione.ritenta and decisione.nuovo_valore == 12
    uno = config_ridotta(tmp_path, run={"batch_size": 1})
    assert politica.decidi("E-S4-02", 1, azione, uno).motivo is Motivo.AZIONE_ESAURITA


def test_e_s5_01_dichiarato_inutile_ferma_con_il_motivo(tmp_path):
    """
    **Obiettivo**: Verificare che S5 non dichiari azioni correttive e che un
    E-S5-01 che dichiara il retry inutile si fermi con il motivo della fase,
    anche senza un'azione correttiva.

    **Razionale scientifico e sistemistico**: La tabella e' una matrice densa
    campioni x varianti: ridurre il lotto non ne cambia la memoria. Il motivo
    specifico dice all'operatore perche' non si ritenta.
    """
    assert TabellaSequenze.aggiustamenti == {}
    config = config_ridotta(tmp_path)
    politica = PoliticaRetry.da_config(config)
    decisione = politica.decidi("E-S5-01", 1, None, config, "matrice densa")
    assert decisione.motivo is Motivo.AZIONE_INUTILE
    assert politica.decidi("E-S5-01", 1, None, config).motivo is Motivo.NESSUNA_AZIONE


# --------------------------------------------------------------------------- #
# Controlli sulla frazione chimerica                                          #
# --------------------------------------------------------------------------- #


def _classe(letture: int, chimeriche: int, varianti: int = 10, var_chim: int = 1) -> dict[str, int]:
    """Le misure delle chimere di una classe con un solo campione."""
    return {
        "campioni": 1, "letture": letture, "letture_chimeriche": chimeriche,
        "varianti": varianti, "varianti_chimeriche": var_chim,
    }


def test_la_frazione_chimerica_si_misura_sulle_letture_delle_classi_controllate(tmp_path):
    """
    **Obiettivo**: Verificare che la frazione controllata sia quella delle
    letture di biologici e positivi insieme, che i negativi ne restino fuori,
    e che per ogni classe si riportino letture e varianti.

    **Razionale scientifico e sistemistico**: Una chimera abbondante pesa piu'
    di cento chimere da una lettura; i bianchi hanno poche letture e una
    frazione che non dice nulla sulla PCR dei campioni.
    """
    misure = misura_chimere({
        "biologico": _classe(900, 90, 100, 30),
        "controllo_positivo": _classe(100, 10),
        "controllo_negativo": _classe(10, 9),
    })
    assert misure["frazione_controllata"] == 0.1
    assert misure["frazione_chimerica"]["biologico"] == {"letture": 0.1, "varianti": 0.3}
    assert misure["frazione_chimerica"]["controllo_negativo"]["letture"] == 0.9
    assert controlla_chimere(misure, config_ridotta(tmp_path).qc) is None


def test_oltre_l_avviso_si_degrada_oltre_l_arresto_si_ferma(tmp_path):
    """
    **Obiettivo**: Verificare che una frazione fra ``qc.warn_frac_chimeric`` e
    ``qc.stop_frac_chimeric`` restituisca l'avviso di E-S6-02, e che oltre
    l'arresto si sollevi E-S6-01, di revisione umana.

    **Razionale scientifico e sistemistico**: Un avviso e un arresto non
    condividono mai lo stesso codice.
    """
    qc = config_ridotta(tmp_path).qc  # 0.25 e 0.50
    avviso = controlla_chimere(misura_chimere({"biologico": _classe(100, 30)}), qc)
    assert avviso is not None and "qc.warn_frac_chimeric" in avviso
    with pytest.raises(ErroreRevisioneUmana) as info:
        controlla_chimere(misura_chimere({"biologico": _classe(100, 51)}), qc)
    assert info.value.codice == "E-S6-01"
    # Solo i negativi oltre la soglia: nessun controllo scatta.
    solo_negativi = misura_chimere({"biologico": _classe(100, 1), "controllo_negativo": _classe(100, 99)})
    assert controlla_chimere(solo_negativi, qc) is None


# --------------------------------------------------------------------------- #
# Le fasi vere, sulla versione ridotta                                        #
# --------------------------------------------------------------------------- #

_MOTIVO_ASSENTI = motivo_pacchetti_r_assenti("dada2", "ggplot2", "ShortRead", "jsonlite")


@pytest.fixture
def dada2():
    """Richiede R con dada2: salta senza, ma in CI fallisce."""
    if _MOTIVO_ASSENTI is not None:
        if os.environ.get("AMPLICON16S_RICHIEDI_BIOC") == "1":
            pytest.fail(f"Bioconductor e' richiesto in questo ambiente: {_MOTIVO_ASSENTI}")
        pytest.skip(_MOTIVO_ASSENTI)


def _copia(base, cartella: Path, **sovrascrivi):
    """L'albero di base, copiato, con i parametri indicati cambiati."""
    from conftest import copia_esecuzione

    return copia_esecuzione(base, cartella, **sovrascrivi)


@pytest.fixture(scope="module")
def catena(catena_calcolata):
    """S4-S7 sulla versione ridotta, condivisa con gli altri moduli, in sola lettura."""
    return catena_calcolata


def _impronte(cartella: Path) -> dict[str, str]:
    """L'MD5 di ogni file della cartella, esclusi i file del ponte e i manifesti."""
    return {
        p.name: hashlib.md5(p.read_bytes()).hexdigest()
        for p in sorted(cartella.iterdir())
        if not p.name.startswith(("rbridge_", "manifest"))
    }


def _r(codice: str) -> str:
    """Esegue un'espressione con Rscript e ne restituisce l'output; fallisce se R esce
    con errore.
    """
    esito = subprocess.run(
        [str(trova_rscript()), "--vanilla", "-e", codice],
        capture_output=True, text=True, check=False,
    )
    assert esito.returncode == 0, esito.stderr
    return esito.stdout


def _tsv(percorso: Path) -> list[dict[str, str]]:
    """Le righe di una tabella separata da tabulazioni."""
    with open(percorso, encoding="utf-8", newline="") as file:
        return list(csv.DictReader(file, delimiter="\t"))


def test_la_catena_s4_s7_si_conclude(dada2, catena):
    """
    **Obiettivo**: Verificare che S4, S5, S6 e S7 si concludano sulla versione
    ridotta, ciascuna con il proprio manifesto, e le metriche attese.

    **Razionale scientifico e sistemistico**: E' il nucleo del denoising: dalle
    letture filtrate alla tabella delle varianti senza chimere.
    """
    run, esito = catena
    assert esito.conclusione is Conclusione.COMPLETATA
    assert [r.passo for r in esito.eseguite] == [Passo.S4, Passo.S5, Passo.S6, Passo.S7]
    for passo in (Passo.S4, Passo.S5, Passo.S6, Passo.S7):
        assert run.valuta().situazioni[passo].stato is StatoPasso.COMPLETATA, passo
    s4 = next(r for r in esito.eseguite if r.passo is Passo.S4)
    assert s4.metriche["campioni"] == 28
    assert s4.metriche["priori"] > 0
    s5 = next(r for r in esito.eseguite if r.passo is Passo.S5)
    assert s5.metriche["letture"] == s4.metriche["letture_denoised"]
    assert s5.metriche["varianti"] == s4.metriche["varianti_distinte"]


def test_il_pseudo_pooling_a_lotti_e_quello_di_dada2(dada2, ridotta_calcolata, tmp_path):
    """
    **Obiettivo**: Verificare che lo script di S4, a lotti di 3 campioni, dia
    per 8 campioni della stessa corsa esattamente le varianti di
    ``dada(pool = "pseudo")`` in una sola chiamata con lo stesso modello.

    **Razionale scientifico e sistemistico**: Le due passate esplicite devono
    riprodurre la regola di dada2 sulle informazioni a priori (prevalenza
    almeno PSEUDO_PREVALENCE) e il modello d'errore che dada ricalcola per la
    seconda passata; se una sola delle due differisse, i risultati no.
    """
    run, _ = ridotta_calcolata
    modelli = run.albero.cartella(Fase.ERROR_MODELS)
    righe = [r for r in _tsv(modelli / "corrispondenza.tsv") if "rerun" in r["modello"]]
    gruppo = sorted(r["campione"] for r in righe)[:8]
    modello = modelli / righe[0]["modello"]
    filtrati = run.albero.cartella(Fase.FILTERED)
    campioni = {a: str(filtrati / f"{a}_filt.fastq.gz") for a in gruppo}
    albero = AlberoOutput(tmp_path / "out")
    esegui_script(
        cartella_r() / "04_dada.R",
        {
            "campioni": campioni, "modelli": {a: str(modello) for a in gruppo},
            "senza_letture": [], "pool": "pseudo", "omega_a": 1e-40, "lotto": 3, "processi": 2,
            "funzione_errore": "loess",
        },
        albero, Fase.ASV_INFERENCE,
        passo="S4",
    )
    file = ", ".join(f'"{a}" = "{p}"' for a, p in campioni.items())
    uscita = _r(
        "suppressPackageStartupMessages(library(dada2)); "
        f"a <- dada(c({file}), err = readRDS('{modello}'), pool = 'pseudo', "
        "OMEGA_A = 1e-40, multithread = 2, verbose = 0); "
        "a <- lapply(a, getUniques); "
        f"b <- readRDS('{albero.cartella(Fase.ASV_INFERENCE) / 'varianti_per_campione.rds'}'); "
        "cat(identical(a[names(b)], b), length(b), sum(lengths(b)))"
    )
    identici, quanti, varianti = uscita.split()
    assert identici == "TRUE"
    assert int(quanti) == 8 and int(varianti) > 0


def test_il_lotto_non_cambia_i_byte_dell_inferenza(dada2, ridotta_calcolata, catena, tmp_path):
    """
    **Obiettivo**: Verificare che S4 con ``run.batch_size`` 5 (sei lotti) dia
    gli stessi byte di S4 con il valore d'esempio (24, due lotti), per ogni
    artefatto.

    **Razionale scientifico e sistemistico**: E' la condizione per escludere
    ``run.batch_size`` dall'impronta e per rendere legittimo il retry di
    E-S4-02 che lo dimezza.
    """
    run, _ = catena
    assert run.config.run.batch_size == 24
    altro = _copia(ridotta_calcolata, tmp_path, run={"batch_size": 5})
    esito = Esecutore(altro, fino_a=Passo.S4).esegui()
    assert esito.conclusione is Conclusione.COMPLETATA
    assert [r.passo for r in esito.eseguite] == [Passo.S4]
    prima = _impronte(run.albero.cartella(Fase.ASV_INFERENCE))
    assert "priori.tsv" in prima and "varianti_per_campione.rds" in prima
    assert _impronte(altro.albero.cartella(Fase.ASV_INFERENCE)) == prima


def test_con_memoria_ridotta_scatta_e_s4_02_e_il_retry_dimezza_il_lotto(
    dada2, ridotta_calcolata, tmp_path, monkeypatch
):
    """
    **Obiettivo**: Verificare che con la memoria del processo R limitata a
    900 MB S4 fallisca con E-S4-02, riconosciuto dal ponte, e che il retry
    riesegua l'inferenza con ``run.batch_size`` dimezzato, fino ai tentativi
    ammessi.

    **Razionale scientifico e sistemistico**: Il limite e' sotto quanto chiede
    anche un lotto di un campione (misurato: circa 1,3 GB di memoria
    virtuale), quindi fallisce ogni tentativo: il test verifica il
    riconoscimento e l'azione correttiva, non il loro successo.
    """
    import amplicon16s.steps.s04_dada as s04

    lotti: list[int] = []
    originale = s04.esegui_script

    def con_limite(script, parametri, *argomenti, **opzioni):
        lotti.append(parametri["lotto"])
        opzioni.update(limite_memoria_byte=900 << 20, tempo_massimo_s=600)
        return originale(script, parametri, *argomenti, **opzioni)

    monkeypatch.setattr(s04, "esegui_script", con_limite)
    run = _copia(ridotta_calcolata, tmp_path)
    esito = Esecutore(run, fino_a=Passo.S4).esegui()
    assert esito.conclusione is Conclusione.ARRESTATA
    assert esito.punto.passo is Passo.S4
    assert esito.punto.codice == "E-S4-02"
    assert esito.punto.tentativi == 2
    assert lotti == [24, 12]
    assert run.albero.manifesto_passo(Passo.S4, Fase.ASV_INFERENCE) is None


def test_oltre_qc_max_asv_count_s5_si_ferma_prima_della_tabella(dada2, catena, tmp_path):
    """
    **Obiettivo**: Verificare che con ``qc.max_asv_count`` sotto le varianti
    della versione ridotta S5 sollevi E-S5-01 prima di costruire la tabella,
    dichiarando inutile il retry, e che la politica non ritenti.

    **Razionale scientifico e sistemistico**: Il numero di varianti e' una
    proprieta' dei dati: un nuovo tentativo darebbe lo stesso numero. La fase
    si chiama direttamente: attraverso l'esecutore, cambiare un parametro qc
    rifarebbe anche S0-S4, perche' S0 dichiara l'intero gruppo qc.
    """
    import logging

    from amplicon16s.config.resolve import risolvi
    from amplicon16s.runner.retry import RITENTARE_INUTILE
    from amplicon16s.steps.base import StepContext

    run, _ = catena
    dati = run.config.model_dump(mode="python")
    dati["io"]["out_root"] = str(tmp_path / "out")
    dati["qc"]["max_asv_count"] = 10
    shutil.copytree(run.config.io.out_root, tmp_path / "out")
    shutil.rmtree(tmp_path / "out" / Fase.SEQTAB.value)
    config = valida(dati)
    fase = TabellaSequenze()
    contesto = StepContext(
        risolta=risolvi(config), albero=AlberoOutput(config.io.out_root),
        logger=logging.getLogger("test"), inventario=run.inventario,
    ).ristretto(fase.parametri)
    with pytest.raises(ErroreRitentabileConRevisione) as info:
        fase.calcola(contesto)
    assert info.value.codice == "E-S5-01"
    assert "oltre le 10 di qc.max_asv_count" in info.value.dettaglio
    motivo = info.value.contesto[RITENTARE_INUTILE]
    assert "matrice densa" in motivo
    assert not (tmp_path / "out" / Fase.SEQTAB.value / "tabella.rds").exists()
    decisione = PoliticaRetry.da_config(config).decidi("E-S5-01", 1, None, config, motivo)
    assert decisione.motivo is Motivo.AZIONE_INUTILE


def test_la_tabella_prima_delle_chimere_resta_in_s5_e_si_ricostruisce(dada2, catena):
    """
    **Obiettivo**: Verificare che S6 non ricopi la tabella di S5, ne registri
    il riferimento con il checksum del manifesto di S5, elenchi le varianti
    chimeriche, e che la tabella senza chimere sia quella di S5 meno le
    colonne elencate, con gli stessi conteggi.

    **Razionale scientifico e sistemistico**: La tabella prima delle chimere
    resta recuperabile e distinta senza raddoppiare su disco l'artefatto
    piu' grande della pipeline.
    """
    run, _ = catena
    chimera = run.albero.cartella(Fase.CHIMERA)
    s6 = run.albero.manifesto_passo(Passo.S6, Fase.CHIMERA)
    s5 = run.albero.manifesto_passo(Passo.S5, Fase.SEQTAB)
    assert "tabella.rds" not in s6.nomi and not (chimera / "tabella.rds").exists()
    riepilogo = json.loads((chimera / "chimere.json").read_text())
    riferimento = riepilogo["tabella_prima_delle_chimere"]
    voce = next(v for v in s5.artefatti if v["nome"] == "tabella.rds")
    assert riferimento == {
        "fase": "S5", "cartella": "06_seqtab", "nome": "tabella.rds", "checksum": voce["checksum"],
    }
    chimere = _tsv(chimera / "chimere.tsv")
    assert all(c["rimossa"] == "si" for c in chimere)  # metodo consensus
    prima = run.albero.cartella(Fase.SEQTAB) / "tabella.rds"
    uscita = _r(
        f"p <- readRDS('{prima}'); s <- readRDS('{chimera / 'tabella_senza_chimere.rds'}'); "
        f"c <- read.delim('{chimera / 'chimere.tsv'}')$sequenza; "
        "cat(identical(p[, !colnames(p) %in% c, drop = FALSE], s), length(c), ncol(p), ncol(s))"
    )
    identiche, chimeriche, prima_n, dopo_n = uscita.split()
    assert identiche == "TRUE"
    assert int(chimeriche) == len(chimere) == s6.metriche["varianti_chimeriche"]
    assert int(prima_n) == int(chimeriche) + int(dopo_n)


def test_la_frazione_chimerica_e_riportata_per_classe(dada2, catena):
    """
    **Obiettivo**: Verificare che S6 riporti per ogni classe la frazione
    chimerica su letture e varianti, e la frazione controllata su biologici e
    positivi, coerente con i conteggi per campione.

    **Razionale scientifico e sistemistico**: Il controllo guarda le letture
    delle classi controllate; le altre misure restano per chi legge.
    """
    run, _ = catena
    chimera = run.albero.cartella(Fase.CHIMERA)
    riepilogo = json.loads((chimera / "chimere.json").read_text())
    assert set(riepilogo["frazione_chimerica"]) == {
        "biologico", "controllo_positivo", "controllo_negativo",
    }
    for misure in riepilogo["frazione_chimerica"].values():
        assert set(misure) == {"letture", "varianti"}
    classi = {r["accession"]: r["classe"] for r in _tsv(
        run.albero.cartella(Fase.INPUT_VALIDATION) / "crosswalk.tsv")}
    per_campione = _tsv(chimera / "chimere_per_campione.tsv")
    controllati = [r for r in per_campione if classi[r["campione"]] != "controllo_negativo"]
    attesa = sum(float(r["letture_chimeriche"]) for r in controllati) / sum(
        float(r["letture"]) for r in controllati)
    assert riepilogo["frazione_controllata"] == round(attesa, 4)


@pytest.mark.parametrize(
    ("tolleranza", "attese"),
    [(0, [137]), (3, [135, 137, 140])],
    ids=["tolleranza-0", "tolleranza-3"],
)
def test_il_filtro_di_lunghezza_toglie_le_varianti_fuori_intervallo(dada2, tmp_path, tolleranza, attese):
    """
    **Obiettivo**: Verificare su varianti sintetiche di lunghezze 120, 135,
    137, 140 e 150 che S7 tenga solo quelle fra ``137 - tolleranza`` e
    ``137 + tolleranza``, con le letture tolte registrate.

    **Razionale scientifico e sistemistico**: Sul dataset di riferimento le
    varianti hanno tutte 137 basi e il filtro non toglie nulla: l'efficacia
    si verifica solo su varianti costruite apposta.
    """
    tabella = tmp_path / "sintetica.rds"
    _r(
        "set.seed(1); l <- c(120, 135, 137, 140, 150); "
        "s <- vapply(l, function(n) paste(sample(c('A','C','G','T'), n, TRUE), collapse = ''), ''); "
        "m <- matrix(c(5L, 1L, 10L, 2L, 3L, 0L, 4L, 7L, 1L, 9L), nrow = 2, "
        "dimnames = list(c('X1', 'X2'), s)); "
        f"saveRDS(m, '{tabella}')"
    )
    albero = AlberoOutput(tmp_path / "out")
    esegui_script(
        cartella_r() / "07_asv_length.R",
        {"tabella": str(tabella), "senza_letture": ["X3"],
         "len_min": 137 - tolleranza, "len_max": 137 + tolleranza},
        albero, Fase.CHIMERA,
        passo="S7",
    )
    cartella = albero.cartella(Fase.CHIMERA)
    tenute = _r(f"cat(sort(nchar(colnames(readRDS('{cartella / 'tabella_asv.rds'}')))))")
    assert [int(x) for x in tenute.split()] == attese
    lunghezze = {int(r["lunghezza"]): r for r in _tsv(cartella / "lunghezze.tsv")}
    assert set(lunghezze) == {120, 135, 137, 140, 150}
    assert {l for l, r in lunghezze.items() if r["ammessa"] == "si"} == set(attese)
    tracciate = {r["campione"]: int(r["letture"]) for r in _tsv(cartella / "letture_lunghezza.tsv")}
    colonne = {120: (5, 1), 135: (10, 2), 137: (3, 0), 140: (4, 7), 150: (1, 9)}
    assert tracciate == {
        "X1": sum(colonne[l][0] for l in attese),
        "X2": sum(colonne[l][1] for l in attese),
        "X3": 0,
    }


def test_sui_dati_reali_il_filtro_di_lunghezza_non_toglie_nulla(dada2, catena):
    """
    **Obiettivo**: Verificare che sulla versione ridotta S7 non tolga alcuna
    variante, perche' hanno tutte la lunghezza del troncamento.

    **Razionale scientifico e sistemistico**: Letture troncate a 137 basi e
    senza fusione di coppie danno varianti di 137 basi: il filtro non ha
    effetto, e non gliene va attribuito.
    """
    run, esito = catena
    s7 = next(r for r in esito.eseguite if r.passo is Passo.S7)
    assert s7.metriche["len_min"] == s7.metriche["len_max"] == 137
    assert s7.metriche["varianti_escluse"] == 0
    assert s7.metriche["letture_escluse"] == 0
    righe = _tsv(run.albero.cartella(Fase.CHIMERA) / "lunghezze.tsv")
    assert [r["lunghezza"] for r in righe] == ["137"]


def test_il_tracciamento_e_completo_un_passo_per_fase(dada2, catena):
    """
    **Obiettivo**: Verificare che il tracciamento ricomposto abbia un passo
    per ciascuna fase da S1 a S7, nell'ordine, per tutti i 28 campioni, con
    letture che non crescono mai da un passo al successivo.

    **Razionale scientifico e sistemistico**: Un buco nel tracciamento
    nasconderebbe dove si perdono le letture di un campione.
    """
    run, _ = catena
    tracciamento = ricomponi(run)
    assert tracciamento.passi == (
        "grezze", "prefiltro", "filtrate", "denoised", "tabella", "senza_chimere", "lunghezza",
    )
    assert tracciamento.origine["denoised"] == "S4"
    assert tracciamento.origine["lunghezza"] == "S7"
    assert len(tracciamento.letture) == 28
    for campione, letture in tracciamento.letture.items():
        valori = [letture[p] for p in tracciamento.passi]
        assert all(a >= b for a, b in zip(valori, valori[1:])), campione
        assert letture["tabella"] == letture["denoised"], campione


# --------------------------------------------------------------------------- #
# Dataset completo, nel container                                              #
# --------------------------------------------------------------------------- #


@pytest.mark.dati_reali
def test_s4_s7_sul_dataset_completo(dada2, catena_reale):
    """
    **Obiettivo**: Verificare che sul dataset completo la catena si concluda
    fino a S7, e riportarne durate, varianti e frazioni chimeriche per classe.

    **Razionale scientifico e sistemistico**: Le misure sul dataset completo
    decidono se i controlli sono tarati sulle classi giuste. La catena e'
    quella condivisa ``catena_reale`` (S0-S14), calcolata una volta per sessione.
    """
    run, esito = catena_reale
    assert esito.conclusione is Conclusione.COMPLETATA
    for risultato in esito.eseguite:
        print(f"\n{risultato.passo}: {risultato.secondi} s, {dict(risultato.metriche)}")
    s7 = next(r for r in esito.eseguite if r.passo is Passo.S7)
    assert s7.metriche["varianti_escluse"] == 0
    tracciamento = ricomponi(run)
    assert len(tracciamento.letture) == 960
