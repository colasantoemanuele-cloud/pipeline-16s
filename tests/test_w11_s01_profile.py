r"""Suite di test per il sottoinsieme di prova OSD-734 e la Fase S1 (Profilo di qualità Phred e lunghezze).

1. Inquadramento nel Piano Operativo
------------------------------------
Settimana 11 (W11), Fase F4 (fase S1: profilo di qualita' Phred e delle
lunghezze; sottoinsieme di prova OSD-734).

2. Moduli sorgente coperti
--------------------------
* ``src/amplicon16s/steps/s01_profile.py``
* ``R/01_profile.R``
* ``scripts/build_test_subset.py``
* ``tests/sottoinsieme.py``

3. Cosa valuta questo file
--------------------------
1. **Integrità e rappresentatività del micro-dataset di test (28 campioni)**:
   verifica che la selezione stratificata e la versione ridotta
   (``tests/fixtures/osd734/ridotto/``, < 3 MB) conservino i 5 controlli
   negativi per piastra richiesti da ``decontam.min_blanks``, gli 8 livelli
   di diluizione dei controlli positivi, la lunghezza minima globale di
   137 bp, i casi estremi di profondità (44 vs 228.431 letture) e il
   tracciamento Git.
2. **Esecuzione di ``R/01_profile.R`` tramite ``rbridge`` e Fase S1**:
   verifica la lettura integrale dei file FASTQ tramite Bioconductor
   ``ShortRead``, la generazione in ``02_qc_profiles/`` dei profili per
   campione e aggregati (``lunghezze.tsv``, ``qualita.tsv``,
   ``letture_grezze.tsv``, ``riepilogo.json``) e la scrittura atomica di
   ``manifest_S1.json``.
3. **Chiusura del limite di G09 e gestione dei troncamenti (137 bp vs 120 bp)**:
   verifica l'assenza di scarto con ``truncLen = 137``, la registrazione
   della degradazione ``E-S1-01`` in S1 (e non in S0) con ``truncLen = 120``
   e l'arresto con ``E-S1-02`` quando ``truncLen`` supera il minimo vero
   dell'intero file non intercettato dal campionamento ``qc.head_reads``.
4. **Riproducibilità dello script e validazione su scala reale (960 campioni)**:
   verifica la ricostruzione deterministica bit-a-bit del subset ridotto e
   l'esecuzione end-to-end di S1 sull'intero dataset NASA GeneLab OSD-734.
5. **Ripresa di S1**: alterato un suo artefatto, si rifa' soltanto S1, con gli
   stessi byte.

4. Comandi Bash e scenari di esecuzione
---------------------------------------
    ``<immagine>`` e' l'immagine del container della pipeline; quella corrente
    e' indicata in ``test.txt``, sezione 1.3.

    1. Modalita' locale standard (senza dati reali ne' Bioconductor R):
       pytest tests/test_w11_s01_profile.py -v

    2. Modalita' locale con dati reali (senza Bioconductor R):
       AMPLICON16S_CONFIG_DATI_REALI=$HOME/ASI/config_osd734.yaml pytest tests/test_w11_s01_profile.py -v

    3. Modalita' container Docker standard (sottoinsieme ridotto con Bioconductor):
       docker run --rm \
         -e PYTHONPATH=/app/src \
         -e AMPLICON16S_R_DIR=/app/R \
         -v "$(pwd)":/app \
         -w /app \
         <immagine> \
         pytest -o cache_dir=/tmp/.pytest_cache tests/test_w11_s01_profile.py -v

    4. Modalita' container Docker completa (con i 960 FASTQ reali):
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
         bash -c "git config --global --add safe.directory /app && pytest -o cache_dir=/tmp/.pytest_cache tests/test_w11_s01_profile.py -v"

    Accorgimenti operativi per il container Docker:
    - Impostare '-e PYTHONPATH=/app/src' per caricare la versione corrente di amplicon16s.
    - Impostare '-e AMPLICON16S_R_DIR=/app/R': l'immagine punta agli script R copiati
      al momento della costruzione, e senza la variabile non userebbe quelli montati.
    - Montare sia il repository ('-v $(pwd):/app') sia i dati reali ('-v $HOME/ASI:$HOME/ASI').
    - Usare '-o cache_dir=/tmp/.pytest_cache' per proteggere i permessi della cartella locale.
    - Aggiungere 'safe.directory /app' in Git per abilitare il collaudo di
      test_i_fastq_di_prova_non_sono_ignorati_da_git, che altrimenti si salta.

5. Risultato atteso
-------------------
16 test totali (16 funzioni di test):
- modalita' 1: 7 passed, 9 skipped in ~1.2s (8 test attendono
  ShortRead/Biostrings, 1 test attende i dati completi);
- modalita' 2: 8 passed, 8 skipped in ~20s (test_lo_script... passa a verde);
- modalita' 3: 13 passed, 3 skipped in ~2 minuti (restano saltati i 2 test sui
  dati reali e il controllo git);
- modalita' 4: 16 passed in ~10 minuti (zero test saltati).

6. Razionale scientifico e sistemistico
---------------------------------------
1. **Profilatura Phred e controllo esatto delle lunghezze a monte di DADA2**:
   il Gate G09 in S0 ispeziona solo le prime ``qc.head_reads`` letture per
   offrire una barriera rapida; S1 legge invece l'intero volume dei file
   FASTQ tramite ``ShortRead``, garantendo che nessuna lettura più corta di
   ``filter.truncLen`` sfugga al controllo prima del filtraggio e
   dell'apprendimento parametrico degli errori di DADA2.
2. **Dataset miniaturizzato e portabile in CI**: il sottoinsieme di 28
   campioni reali sottocampionati con seme fisso permette di eseguire i
   calcoli R/Bioconductor in Continuous Integration in pochi secondi,
   mantenendo intatte tutte le condizioni al contorno biologiche e
   sperimentali di OSD-734.
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import itertools
import json
import logging
import os
import subprocess
import sys
from collections import Counter
from pathlib import Path

import pytest
from conftest import copia_esecuzione
from sottoinsieme import RIDOTTO, config_ridotta, manifesto, selezione

from amplicon16s.config.schema import carica, valida
from amplicon16s.io_layer.artifacts import Fase
from amplicon16s.logging.logger import NOME_FILE_LOG, chiudi, configura
from amplicon16s.rbridge.runner import trova_rscript
from amplicon16s.runner.executor import Conclusione, Esecutore
from amplicon16s.runner.graph import Passo
from amplicon16s.runner.project import ProjectRun

RADICE = Path(__file__).resolve().parents[1]


@pytest.fixture(autouse=True)
def uscite_pulite():
    """Chiude le uscite del log prima e dopo ogni test, perché nessun handler resti
    aperto sulla cartella temporanea.
    """
    chiudi()
    yield
    chiudi()


# --------------------------------------------------------------------------- #
# La selezione e la versione ridotta                                           #
# --------------------------------------------------------------------------- #


def test_la_selezione_e_motivata_campione_per_campione():
    """Verifica la completezza formale di ``selezione.tsv`` sui 28 campioni scelti.

    * **Obiettivo**: controllare che la tabella di selezione contenga esattamente
      28 accession distinti e che ciascuna riga dichiari un motivo esplicito,
      un'impronta MD5 a 32 caratteri esadecimali e un identificativo ``run_ena``
      con prefisso ``ERR``.
    * **Razionale scientifico e sistemistico**: garantisce la tracciabilità verso
      l'archivio pubblico ENA e impedisce l'inclusione nel sottoinsieme di prova
      di campioni privi di giustificazione metodologica o checksum verificabile.
    """
    righe = selezione()
    assert len(righe) == 28
    assert len({r["accession"] for r in righe}) == 28
    for riga in righe:
        assert riga["motivo"].strip(), riga["accession"]
        assert len(riga["md5"]) == 32
        assert riga["run_ena"].startswith("ERR")


def test_la_selezione_esercita_decontaminazione_e_calibrazione():
    """Verifica il bilanciamento sperimentale di classi, corse e piastre nella selezione.

    * **Obiettivo**: accertare che i 28 campioni coprano 9 biologici, 9 controlli
      positivi e 10 controlli negativi distribuiti su 2 corse e 2 piastre (piastre
      4 e 10, ciascuna con 5 blank e campioni biologici associati), oltre alla
      serie completa di 8 diluizioni di controllo positivo nella piastra 10.
    * **Razionale scientifico e sistemistico**: assicura che il sottoinsieme
      soddisfi il requisito statistico ``decontam.min_blanks = 5`` per batch e
      fornisca l'intera curva di calibrazione della mock community per le fasi
      successive (S11 e S12).
    """
    righe = selezione()
    # Le tre classi e le due corse.
    assert Counter(r["classe"] for r in righe) == {
        "biologico": 9, "controllo_positivo": 9, "controllo_negativo": 10,
    }
    assert len({r["corsa"] for r in righe}) == 2
    # Due piastre con i cinque negativi che decontam.min_blanks richiede, e
    # biologici della stessa piastra con cui confrontarli.
    for piastra in ("4", "10"):
        della = [r for r in righe if r["piastra"] == piastra]
        assert sum(r["classe"] == "controllo_negativo" for r in della) == 5
        assert any(r["classe"] == "biologico" for r in della)
    # Una serie di calibrazione completa, otto livelli nella stessa piastra.
    serie = [r for r in righe if r["classe"] == "controllo_positivo" and r["piastra"] == "10"]
    assert len(serie) == 8


def test_la_selezione_contiene_i_casi_noti():
    """Verifica la presenza dei campioni con casi limite noti di OSD-734.

    * **Obiettivo**: verificare che nella selezione figurino l'accession
      ``ERX12083125`` (con la lunghezza minima globale di 137 bp), il campione
      ``JLP1A1.L4`` (``ERX12083170``) e gli estremi di profondità di lettura
      (da 44 a 228.431 letture).
    * **Razionale scientifico e sistemistico**: preserva nel subset di test i
      casi limite reali che determinano le soglie di troncamento e la tenuta
      della pipeline su librerie a bassissima o altissima copertura.
    """
    per_accession = {r["accession"]: r for r in selezione()}
    assert per_accession["ERX12083125"]["lunghezza_minima"] == "137"  # BLANK.3DMM.P4.8C
    assert per_accession["ERX12083170"]["campione"] == "JLP1A1.L4"
    letture = sorted(int(r["letture"]) for r in per_accession.values())
    assert letture[0] == 44 and letture[-1] == 228431  # il piu' povero e il piu' profondo


def test_la_versione_ridotta_e_integra_e_piccola():
    """Verifica l'integrità MD5 e il peso complessivo (< 3 MB) dei FASTQ ridotti.

    * **Obiettivo**: controllare la corrispondenza esatta tra ``manifesto.tsv``
      e i file in ``ridotto/fastq/``, verificare l'impronta MD5 di ogni archivio
      ``.fastq.gz`` e accertare che l'ingombro totale resti sotto i 3 MB.
    * **Razionale scientifico e sistemistico**: garantisce che gli archivi
      versionati nel repository non siano corrotti e mantengano un'impronta su
      disco compatibile con la clonazione rapida nei runner di Continuous
      Integration.
    """
    voci = manifesto()
    cartella = RIDOTTO / "fastq"
    assert {v["file"] for v in voci} == {p.name for p in cartella.iterdir()}
    for voce in voci:
        percorso = cartella / voce["file"]
        assert hashlib.md5(percorso.read_bytes()).hexdigest() == voce["md5"], voce["file"]
    totale = sum(p.stat().st_size for p in cartella.iterdir())
    assert totale < 3 * 1024 * 1024


def _lunghezze(percorso: Path) -> list[int]:
    """Le lunghezze delle sequenze di un file FASTQ compresso, nell'ordine del file."""
    with gzip.open(percorso, "rt") as file:
        return [len(r.rstrip("\n")) for i, r in enumerate(file) if i % 4 == 1]


def test_la_versione_ridotta_conserva_le_proprieta_della_selezione():
    """Verifica che il sottocampionamento conservi lunghezze minime e rapporto di profondità.

    * **Obiettivo**: verificare che ciascun file FASTQ ridotto contenga il numero
      di letture dichiarato nel manifesto, conservi invariata la propria
      lunghezza minima originale (incluso il minimo globale di 137 bp) e
      mantenga un rapporto > 50x tra il campione più profondo e quello più
      povero (che conserva tutte le 44 letture).
    * **Razionale scientifico e sistemistico**: evita che il sottocampionamento
      casuale elimini le rare letture corte (137 bp) che governano il controllo
      di troncamento in Fase S1 o appiattisca le differenze di copertura tra
      campioni.
    """
    per_accession = {r["accession"]: r for r in selezione()}
    for voce in manifesto():
        accession = voce["file"].split("_")[0]
        lunghezze = _lunghezze(RIDOTTO / "fastq" / voce["file"])
        assert len(lunghezze) == int(voce["letture"])
        # Ogni file conserva la propria lunghezza minima, 137 compreso.
        assert min(lunghezze) == int(per_accession[accession]["lunghezza_minima"])
    letture = {v["file"].split("_")[0]: int(v["letture"]) for v in manifesto()}
    # Il contrasto fra il campione piu' profondo e il piu' povero resta.
    assert letture["ERX12084006"] > 50 * letture["ERX12083285"]
    assert letture["ERX12083285"] == 44  # sotto il minimo: tutte le letture


def test_i_metadati_ridotti_coprono_esattamente_la_selezione():
    """Verifica l'allineamento tra ``metadati/lotti.tsv`` e i 28 accession selezionati.

    * **Obiettivo**: accertare che l'insieme degli ``experiment_accession``
      presenti nella tabella ``lotti.tsv`` del subset ridotto coincida
      esattamente con i 28 accession di ``selezione.tsv``.
    * **Razionale scientifico e sistemistico**: garantisce la coerenza
      relazionale del crosswalk per il sottoinsieme di prova senza righe orfane
      o campioni mancanti.
    """
    accession = {r["accession"] for r in selezione()}
    with open(RIDOTTO / "metadati" / "lotti.tsv", encoding="utf-8") as file:
        assert {r["experiment_accession"] for r in csv.DictReader(file, delimiter="\t")} == accession


def test_i_fastq_di_prova_non_sono_ignorati_da_git():
    """Verifica che le regole di ``.gitignore`` non escludano i FASTQ del subset ridotto.

    * **Obiettivo**: eseguire ``git check-ignore`` su un file di
      ``tests/fixtures/osd734/ridotto/fastq/`` e verificare il codice di uscita 1
      (file non ignorato).
    * **Razionale scientifico e sistemistico**: mentre i FASTQ grezzi di
      produzione sono ignorati per evitare di versionare gigabyte di dati,
      l'eccezione in ``.gitignore`` per ``tests/fixtures/osd734/ridotto/`` è
      indispensabile affinché i runner CI ricevano i file di test.
    """
    dentro = subprocess.run(
        ["git", "rev-parse", "--is-inside-work-tree"],
        cwd=RADICE, capture_output=True, check=False,
    )
    if dentro.returncode != 0:
        pytest.skip("non e' un repository git: nell'immagine .git non viene copiato")
    esempio = next((RIDOTTO / "fastq").iterdir())
    esito = subprocess.run(
        ["git", "check-ignore", "-q", str(esempio)], cwd=RADICE, check=False
    )
    assert esito.returncode == 1  # 1: il file non e' ignorato


# --------------------------------------------------------------------------- #
# S1 sulla versione ridotta                                                    #
# --------------------------------------------------------------------------- #


def _motivo_bioconductor_assente() -> str | None:
    """Il motivo per saltare i test di S1 (Rscript, ShortRead, Biostrings o jsonlite
    assenti), o ``None``.
    """
    rscript = trova_rscript()
    if rscript is None:
        return "Rscript non disponibile"
    prova = subprocess.run(
        [str(rscript), "--vanilla", "-e",
         "quit(status = !all(vapply(c('ShortRead', 'Biostrings', 'jsonlite'), "
         "requireNamespace, logical(1), quietly = TRUE)))"],
        capture_output=True, check=False,
    )
    if prova.returncode != 0:
        return "ShortRead, Biostrings o jsonlite non installati"
    return None


_MOTIVO_BIOC_ASSENTE = _motivo_bioconductor_assente()


@pytest.fixture
def bioconductor():
    """Richiede R con Bioconductor: salta senza, ma in CI fallisce."""
    if _MOTIVO_BIOC_ASSENTE is not None:
        if os.environ.get("AMPLICON16S_RICHIEDI_BIOC") == "1":
            pytest.fail(f"Bioconductor e' richiesto in questo ambiente: {_MOTIVO_BIOC_ASSENTE}")
        pytest.skip(_MOTIVO_BIOC_ASSENTE)


def _fino_a_s1(config):
    """Esegue il grafo fino a S1 e restituisce l'esecuzione e il suo esito."""
    run = ProjectRun(config)
    return run, Esecutore(run, fino_a=Passo.S1).esegui()


def _profili(run) -> Path:
    """La cartella dei profili di qualità dell'esecuzione."""
    return run.albero.cartella(Fase.QC_PROFILES)


@pytest.fixture(scope="module")
def eseguita(ridotta_calcolata):
    """S0 e S1 sulla versione ridotta, una volta sola per i test che leggono.

    E' l'esecuzione condivisa, che esegue S1 a parte (``ridotta_calcolata``).
    Caricare ShortRead costa da solo alcuni secondi: ripeterlo per ogni test
    non aggiungerebbe nulla alla verifica.
    """
    return ridotta_calcolata


def test_s1_scrive_i_profili_con_il_proprio_manifesto(bioconductor, eseguita):
    """Verifica che la Fase S1 produca i 4 artefatti attesi e registri ``manifest_S1.json``.

    * **Obiettivo**: controllare che l'esecuzione di S0 e S1 termini con stato
      ``COMPLETATA``, che ``manifest_S1.json`` in ``02_qc_profiles/`` elenchi
      ``lunghezze.tsv``, ``qualita.tsv``, ``letture_grezze.tsv`` e
      ``riepilogo.json`` e che ``ProjectRun.valuta()`` riconosca S1 come
      completata.
    * **Razionale scientifico e sistemistico**: garantisce l'integrità
      crittografica SHA-256 degli output di profilatura prodotti da R e la
      possibilità di ripresa deterministica (``resume``) senza ricalcolare S1.
    """
    run, esito = eseguita
    assert esito.conclusione is Conclusione.COMPLETATA

    manifesto_s1 = run.albero.manifesto_passo(Passo.S1, Fase.QC_PROFILES)
    assert set(manifesto_s1.nomi) == {
        "lunghezze.tsv", "qualita.tsv", "letture_grezze.tsv", "riepilogo.json",
        "valori_qualita.tsv",
    }
    assert run.valuta().situazioni[Passo.S1].stato.value == "completata"


def test_la_ripresa_rifa_solo_s1_se_un_suo_artefatto_e_alterato(bioconductor, eseguita, tmp_path):
    """Verifica che, alterato ``lunghezze.tsv`` di S1, la ripresa rifaccia soltanto S1.

    * **Obiettivo**: in una copia dell'esecuzione, alterare un artefatto di S1
      e controllare che la valutazione dia S1 da rifare, per quell'artefatto,
      lasciando concluse S0, S2 e S3; che l'esecuzione fino a S1 rifaccia
      soltanto S1; e che gli artefatti rifatti siano identici agli originali.
    * **Razionale scientifico e sistemistico**: S1 non ha fasi realizzate a
      valle: alterarne un profilo non deve rifare il filtro e i modelli
      d'errore, e la profilatura rifatta deve dare gli stessi byte.
    """
    originale, _ = eseguita
    run = copia_esecuzione(eseguita, tmp_path)
    artefatto = run.albero.cartella(Fase.QC_PROFILES) / "lunghezze.tsv"
    artefatto.write_bytes(artefatto.read_bytes() + b"\n")

    situazione = run.valuta().situazioni
    assert situazione[Passo.S1].stato.value == "da_eseguire"
    assert "lunghezze.tsv" in situazione[Passo.S1].motivo
    for altro in (Passo.S0, Passo.S2, Passo.S3):
        assert situazione[altro].stato.value == "completata", altro

    esito = Esecutore(run, fino_a=Passo.S1).esegui()
    assert esito.conclusione is Conclusione.COMPLETATA
    assert [r.passo for r in esito.eseguite] == [Passo.S1]

    def impronte(cartella):
        return {
            f.name: hashlib.md5(f.read_bytes()).hexdigest()
            for f in sorted(cartella.iterdir()) if not f.name.startswith(("rbridge_", "manifest"))
        }

    assert impronte(run.albero.cartella(Fase.QC_PROFILES)) == impronte(
        originale.albero.cartella(Fase.QC_PROFILES)
    )


def test_s1_legge_tutte_le_letture(bioconductor, eseguita):
    """Verifica che ``R/01_profile.R`` ispezioni l'intero contenuto di ogni file FASTQ.

    * **Obiettivo**: confrontare i conteggi di ``letture_grezze.tsv`` e di
      ``riepilogo.json`` con i valori esatti del manifesto dei 28 campioni,
      verificando anche ``lunghezza_minima == 137`` e ``moda == 151``.
    * **Razionale scientifico e sistemistico**: dimostra che S1 supera il
      campionamento iniziale ``qc.head_reads`` di S0 e censisce il 100% delle
      letture grezze prima del filtraggio.
    """
    run, _ = eseguita
    attese = {v["file"].split("_")[0]: int(v["letture"]) for v in manifesto()}
    with open(_profili(run) / "letture_grezze.tsv", encoding="utf-8") as file:
        contate = {r["campione"]: int(r["letture"]) for r in csv.DictReader(file, delimiter="\t")}
    assert contate == attese

    riepilogo = json.loads((_profili(run) / "riepilogo.json").read_text())
    assert riepilogo["campioni"] == 28
    assert riepilogo["letture"] == sum(attese.values())
    assert riepilogo["lunghezza_minima"] == 137
    assert riepilogo["moda"] == 151


def test_s1_scrive_distribuzione_e_profilo_per_campione(bioconductor, eseguita):
    """Verifica l'esattezza delle distribuzioni di lunghezza e dei quantili Phred per posizione.

    * **Obiettivo**: controllare per il campione ``ERX12083125`` che la
      distribuzione in ``lunghezze.tsv`` coincida con il conteggio diretto sul
      FASTQ e che ``qualita.tsv`` riporti le posizioni 1..151 con l'ordinamento
      monotono dei quantili ``p10 <= mediana <= p90`` e la copertura esatta alle
      posizioni 137 e 151.
    * **Razionale scientifico e sistemistico**: valida la correttezza del calcolo
      posizionale Phred in ``R/01_profile.R`` anche in presenza di letture di
      lunghezza variabile (tra 137 e 151 bp).
    """
    run, _ = eseguita
    with open(_profili(run) / "lunghezze.tsv", encoding="utf-8") as file:
        lunghezze = list(csv.DictReader(file, delimiter="\t"))
    blank = {int(r["lunghezza"]): int(r["letture"]) for r in lunghezze if r["campione"] == "ERX12083125"}
    attese = Counter(_lunghezze(next((RIDOTTO / "fastq").glob("ERX12083125_*"))))
    assert blank == dict(attese)

    with open(_profili(run) / "qualita.tsv", encoding="utf-8") as file:
        qualita = [r for r in csv.DictReader(file, delimiter="\t") if r["campione"] == "ERX12083125"]
    assert [int(r["posizione"]) for r in qualita] == list(range(1, 152))
    # Alla posizione 137 ci sono tutte le letture, alla 151 solo le complete.
    assert int(qualita[136]["letture"]) == sum(attese.values())
    assert int(qualita[150]["letture"]) == attese[151]
    for riga in qualita:
        assert int(riga["p10"]) <= int(riga["mediana"]) <= int(riga["p90"])


def test_con_truncLen_137_non_c_e_scarto(bioconductor, eseguita):
    """Verifica l'assenza di degradazioni quando ``filter.truncLen`` coincide con il minimo (137 bp).

    * **Obiettivo**: accertare che con ``truncLen = 137`` la tupla
      ``degradazioni`` di ``manifest_S1.json`` sia vuota.
    * **Razionale scientifico e sistemistico**: conferma che la parametrizzazione
      di riferimento per OSD-734 (137 bp) conserva tutte le letture senza
      attivare avvisi di perdita di basi (``E-S1-01``).
    """
    run, _ = eseguita
    assert run.albero.manifesto_passo(Passo.S1, Fase.QC_PROFILES).degradazioni == ()


def test_con_truncLen_120_lo_scarto_e_registrato_da_s1_non_da_s0(bioconductor, tmp_path):
    """Verifica che la degradazione ``E-S1-01`` per troncamento conservativo sia emessa da S1.

    * **Obiettivo**: eseguire S0 e S1 con ``filter.truncLen = 120`` (scarto di
      17 bp rispetto a 137 bp, superiore a ``truncLen_shortfall_warn = 10``) e
      verificare che ``E-S1-01`` compaia esclusivamente nel manifesto e nei log
      della Fase S1 e non in S0.
    * **Razionale scientifico e sistemistico**: attribuisce la segnalazione di
      perdita di risoluzione tassonomica alla fase che dispone del minimo reale
      su tutte le letture (S1), anziché alla stima preliminare di S0.
    """
    configura(tmp_path / "out")
    run, esito = _fino_a_s1(config_ridotta(tmp_path, filter={"truncLen": 120}))
    assert esito.conclusione is Conclusione.COMPLETATA

    (degradazione,) = run.albero.manifesto_passo(Passo.S1, Fase.QC_PROFILES).degradazioni
    assert degradazione["codice"] == "E-S1-01"
    assert degradazione["scarto"] == 17
    s0 = run.albero.manifesto_passo(Passo.S0, Fase.INPUT_VALIDATION).degradazioni
    assert "E-S1-01" not in {d["codice"] for d in s0}

    righe = [json.loads(r) for r in (tmp_path / "out" / "99_logs" / NOME_FILE_LOG).read_text().splitlines()]
    eventi = [r for r in righe if r.get("codice") == "E-S1-01"]
    assert len(eventi) == 1 and eventi[0]["fase"] == "S1"


def test_s1_ferma_se_il_troncamento_supera_il_minimo_vero(bioconductor, tmp_path):
    """Verifica la chiusura del limite di G09 mediante arresto ``E-S1-02`` in Fase S1.

    * **Obiettivo**: impostare ``qc.head_reads = 20`` e ``filter.truncLen`` pari
      al minimo delle prime 20 letture (> 137 bp), verificando che il Gate G09 in
      S0 passi ma che S1 individui la lettura da 137 bp oltre la testa del file,
      arresti l'esecuzione con codice ``E-S1-02`` e lasci ``riepilogo.json`` su
      disco senza scrivere ``manifest_S1.json``.
    * **Razionale scientifico e sistemistico**: impedisce che un troncamento
      eccessivo scarti silenziosamente letture durante ``filterAndTrim`` in S2,
      conservando al contempo su disco i profili di lunghezza per consentire
      all'operatore di calibrare ``filter.truncLen``.
    """
    def prime(percorso, k):
        with gzip.open(percorso, "rt") as file:
            return [len(r.rstrip()) for i, r in enumerate(itertools.islice(file, 4 * k)) if i % 4 == 1]

    stima = min(min(prime(p, 20)) for p in (RIDOTTO / "fastq").iterdir())
    assert stima > 137
    # Dalla settimana 28 S1 giudica la frazione di letture piu' corte: con la
    # frazione ammessa a zero, una sola lettura corta oltre la testa ferma.
    config = config_ridotta(tmp_path, qc={"head_reads": 20, "max_frac_short_reads": 0.0},
                            filter={"truncLen": stima})
    run, esito = _fino_a_s1(config)

    assert esito.conclusione is Conclusione.ARRESTATA
    assert [r.passo for r in esito.eseguite] == [Passo.S0]  # G09 e' passato
    assert esito.punto.passo is Passo.S1
    assert esito.punto.codice == "E-S1-02"
    assert esito.punto.categoria == "revisione_umana"
    assert esito.punto.tentativi == 1
    assert "137 bp" in esito.punto.dettaglio
    # I profili restano per decidere il troncamento, ma S1 non e' conclusa.
    assert (_profili(run) / "riepilogo.json").exists()
    assert run.albero.manifesto_passo(Passo.S1, Fase.QC_PROFILES) is None


# --------------------------------------------------------------------------- #
# Dataset completo, nel container                                              #
# --------------------------------------------------------------------------- #


@pytest.mark.dati_reali
def test_s1_sul_dataset_completo(bioconductor, tmp_path):
    """Verifica l'esecuzione end-to-end di S0 e S1 sui 960 file FASTQ reali di OSD-734.

    * **Obiettivo**: eseguire la pipeline fino a S1 sull'intero dataset reale
      indicato da ``AMPLICON16S_CONFIG_DATI_REALI`` e verificare in
      ``riepilogo.json`` i valori globali attesi: 960 campioni, lunghezza minima
      137 bp, moda 151 bp e qualità mediana minima >= 25.
    * **Razionale scientifico e sistemistico**: valida empiricamente su scala
      reale (2.4 GB compressi) le prestazioni, la tenuta di memoria dello
      streaming Bioconductor e le proprietà biofisiche del dataset NASA GeneLab
      OSD-734.
    """
    percorso = os.environ.get("AMPLICON16S_CONFIG_DATI_REALI")
    if not percorso or not Path(percorso).expanduser().is_file():
        pytest.skip("AMPLICON16S_CONFIG_DATI_REALI non impostata o file assente")
    dati = carica(Path(percorso).expanduser()).model_dump(mode="python")
    dati["io"]["out_root"] = str(tmp_path / "out")
    logging.getLogger("amplicon16s").setLevel(logging.WARNING)

    run, esito = _fino_a_s1(valida(dati))
    assert esito.conclusione is Conclusione.COMPLETATA
    riepilogo = json.loads((_profili(run) / "riepilogo.json").read_text())
    assert riepilogo["campioni"] == 960
    assert riepilogo["lunghezza_minima"] == 137
    assert riepilogo["moda"] == 151
    assert riepilogo["qualita_mediana_minima"] >= 25


@pytest.mark.dati_reali
def test_lo_script_ricostruisce_la_versione_ridotta_identica(tmp_path, monkeypatch, capsys):
    """Verifica la riproducibilità bit-a-bit di ``scripts/build_test_subset.py ridotto``.

    * **Obiettivo**: invocare ``build_test_subset.py ridotto`` a partire dai dati
      reali locali verso una directory temporanea e verificare che
      ``manifesto.tsv`` (con gli MD5 dei file ``.fastq.gz``) e le tabelle in
      ``metadati/`` siano identici byte per byte a quelli versionati nel
      repository.
    * **Razionale scientifico e sistemistico**: garantisce che il seme fisso e la
      compressione gzip priva di timestamp producano artefatti di test
      integralmente riproducibili su qualunque macchina.
    """
    percorso = os.environ.get("AMPLICON16S_CONFIG_DATI_REALI")
    if not percorso or not Path(percorso).expanduser().is_file():
        pytest.skip("AMPLICON16S_CONFIG_DATI_REALI non impostata o file assente")
    import importlib.util

    spec = importlib.util.spec_from_file_location(
        "build_test_subset", RADICE / "scripts" / "build_test_subset.py"
    )
    modulo = importlib.util.module_from_spec(spec)
    # Le dataclass del modulo cercano il modulo in sys.modules.
    monkeypatch.setitem(sys.modules, "build_test_subset", modulo)
    spec.loader.exec_module(modulo)

    modulo.main(["ridotto", "--config", str(Path(percorso).expanduser()),
                 "--destinazione", str(tmp_path)])
    rigenerato = (tmp_path / "manifesto.tsv").read_text(encoding="utf-8")
    assert rigenerato == (RIDOTTO / "manifesto.tsv").read_text(encoding="utf-8")
    for nome in ("assay.txt", "studio.txt", "lotti.tsv"):
        assert (tmp_path / "metadati" / nome).read_bytes() == (RIDOTTO / "metadati" / nome).read_bytes()
