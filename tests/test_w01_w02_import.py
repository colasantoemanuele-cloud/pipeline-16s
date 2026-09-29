r"""Smoke test di installazione dei pacchetti, assenza di import circolari e interfaccia CLI.

1. Inquadramento nel Piano Operativo
------------------------------------
Settimane 1 e 2 (W1/W2), Fase F1 (inizializzazione del repository, namespace
dei pacchetti, smoke test di integrazione continua e interfaccia CLI di base).

2. Moduli sorgente coperti
--------------------------
* ``src/amplicon16s/__init__.py``
* ``src/amplicon16s_eco/__init__.py``
* ``src/amplicon16s/cli.py``

3. Cosa valuta questo file
--------------------------
- importabilità dei due pacchetti del progetto (``amplicon16s``, pipeline di
  produzione, e ``amplicon16s_eco``, analisi ecologiche) senza ``ImportError``
  né cicli di importazione, ciascuno con un attributo ``__version__`` non vuoto;
- registrazione nel parser ``argparse`` dei quattro sottocomandi ``run``,
  ``resume``, ``validate`` e ``report``, ciascuno dei quali risponde a
  ``--help`` con codice di uscita 0;
- uscita con codice 2 (riga di comando non valida) di ciascun sottocomando
  invocato senza l'argomento obbligatorio ``--config``.

4. Comandi Bash e scenari di esecuzione
---------------------------------------
    1. Modalità locale standard (R di base con jsonlite, senza Bioconductor né
       dati reali):
       pytest tests/test_w01_w02_import.py -v

    2. Modalità container Docker standard (sottoinsieme ridotto con
       R/Bioconductor):
       docker run --rm \
         -e PYTHONPATH=/app/src \
         -v "$(pwd)":/app \
         -w /app \
         amplicon16s:dev \
         pytest -o cache_dir=/tmp/.pytest_cache tests/test_w01_w02_import.py -v

    3. Modalità container Docker completa (con i 2.4 GB di dati reali OSD-734;
       la configurazione e i percorsi che contiene devono stare nella cartella
       montata):
       docker run --rm \
         -e PYTHONPATH=/app/src \
         -e AMPLICON16S_CONFIG_DATI_REALI="$HOME/ASI/config_osd734.yaml" \
         -v "$(pwd)":/app \
         -v "$HOME/ASI":"$HOME/ASI" \
         -w /app \
         amplicon16s:dev \
         pytest -o cache_dir=/tmp/.pytest_cache tests/test_w01_w02_import.py -v

5. Risultato atteso
-------------------
10 test totali:
- 10 passed in ambiente locale standard (~0.1s);
- 10 passed nel container Docker standard sul sottoinsieme ridotto (~0.1s);
- 10 passed nel container Docker con i dati reali OSD-734 (~0.1s).

6. Razionale scientifico e sistemistico
---------------------------------------
- È il primo presidio di integrazione continua: un pacchetto non importabile
  o un ciclo di importazione renderebbero inutilizzabile l'intera pipeline, e
  vanno intercettati prima di qualunque test funzionale.
- Il contratto della riga di comando (sottocomandi e codici di uscita) è
  l'interfaccia su cui si appoggiano gli script dell'operatore: un sottocomando
  scomparso o un codice di uscita cambiato romperebbero l'automazione a valle
  senza alcun errore visibile.
"""

import importlib

import pytest


def test_amplicon16s_importabile():
    """
    **Obiettivo**: Verificare che il pacchetto principale ``amplicon16s`` sia
    importabile tramite ``importlib.import_module`` ed esponga una stringa
    ``__version__`` non vuota.

    **Razionale scientifico e sistemistico**: Garantisce che l'installazione
    ``pyproject.toml`` e la catena dei moduli interni di ``amplicon16s`` siano
    prive di errori di sintassi o di *circular imports*, e che la versione del
    software sia disponibile per essere registrata nei manifesti di riproducibilità.
    """
    modulo = importlib.import_module("amplicon16s")
    assert modulo.__version__


def test_amplicon16s_eco_importabile():
    """
    **Obiettivo**: Verificare che il pacchetto complementare ``amplicon16s_eco``
    sia importabile ed esponga l'attributo ``__version__``.

    **Razionale scientifico e sistemistico**: Certifica fin dalla Settimana 1 la
    corretta separazione architetturale tra il motore di processamento primario
    16S (``amplicon16s``) e il namespace dedicato alle analisi ecologiche e
    multivariate a valle (``amplicon16s_eco``).
    """
    modulo = importlib.import_module("amplicon16s_eco")
    assert modulo.__version__


@pytest.mark.parametrize("comando", ["run", "resume", "validate", "report"])
def test_sottocomandi_rispondono(comando, capsys):
    """
    **Obiettivo**: Verificare che ciascuno dei 4 sottocomandi della CLI
    (``run``, ``resume``, ``validate``, ``report``) invocato con ``--help``
    termini con ``SystemExit(0)`` stampando il nome del comando su ``stdout``.

    **Razionale scientifico e sistemistico**: Accerta che l'interfaccia a riga di
    comando esponga tutti e quattro i punti di ingresso operativi previsti dal
    protocollo (esecuzione completa, ripresa da checkpoint, sola validazione S0
    e generazione del report) e che la documentazione ``--help`` sia accessibile
    dentro il container Docker/Apptainer.
    """
    from amplicon16s.cli import main

    with pytest.raises(SystemExit) as uscita:
        main([comando, "--help"])
    assert uscita.value.code == 0
    assert comando in capsys.readouterr().out


@pytest.mark.parametrize("comando", ["run", "resume", "validate", "report"])
def test_sottocomandi_chiedono_la_configurazione(comando):
    """
    **Obiettivo**: Verificare che invocare ciascuno dei 4 sottocomandi senza
    argomenti provochi l'uscita immediata di ``argparse`` con codice ``SystemExit(2)``.

    **Razionale scientifico e sistemistico**: Impedisce che la pipeline possa
    essere avviata con percorsi o parametri impliciti non dichiarati: ogni
    esecuzione deve ricevere esplicitamente il file di configurazione YAML dello
    studio (es. ``--config config_osd734.yaml``).
    """
    from amplicon16s.cli import main

    with pytest.raises(SystemExit) as uscita:
        main([comando])
    assert uscita.value.code == 2
