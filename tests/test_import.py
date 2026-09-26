"""Smoke test di installazione dei pacchetti, assenza di import circolari e interfaccia CLI.

Inquadramento nel Piano Operativo
---------------------------------
* **Settimane di riferimento**: **Settimana 1 e Settimana 2 (W1/W2 — Fase F1:
  Inizializzazione repository, package namespace, CI smoke tests e interfaccia
  CLI base)**.
* **Moduli sorgente coperti**:
  - ``src/amplicon16s/__init__.py``
  - ``src/amplicon16s_eco/__init__.py``
  - ``src/amplicon16s/cli.py``

Scopo sperimentale e razionale scientifico/sistemistico
-------------------------------------------------------
Questo modulo rappresenta il primo presidio di integrazione continua (CI)
introdotto fin dalla Settimana 1 per garantire due proprietà fondamentali:

1. **Integrità del packaging e assenza di import circolari**: verifica che i due
   namespace Python del progetto — ``amplicon16s`` (core di validazione,
   orchestrazione, gate e ponte R) e ``amplicon16s_eco`` (estensioni di analisi
   ecologica e modellistica) — siano correttamente installati nell'ambiente
   virtuale/container, importabili dinamicamente senza ``ImportError`` o cicli
   di dipendenza e provvisti di attributo ``__version__``.
2. **Contratto sintattico della CLI (``amplicon16s``)**: verifica che tutti e 4
   i sottocomandi previsti dal Piano Operativo (``run``, ``resume``,
   ``validate``, ``report``) siano registrati nel parser ``argparse``, rispondano
   con codice di uscita ``0`` all'opzione ``--help`` ed escano con codice POSIX
   ``2`` quando invocati senza gli argomenti obbligatori (es. ``--config``).
"""

import importlib

import pytest


def test_amplicon16s_importabile():
    """
    **Obiettivo**: Verificare che il pacchetto principale ``amplicon16s`` sia
    importabile tramite ``importlib.import_module`` ed esponga una stringa
    ``__version__`` non vuota.

    **Razionale Scientifico/Sistemistico**: Garantisce che l'installazione
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

    **Razionale Scientifico/Sistemistico**: Certifica fin dalla Settimana 1 la
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

    **Razionale Scientifico/Sistemistico**: Accerta che l'interfaccia a riga di
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

    **Razionale Scientifico/Sistemistico**: Impedisce che la pipeline possa
    essere avviata con percorsi o parametri impliciti non dichiarati: ogni
    esecuzione deve ricevere esplicitamente il file di configurazione YAML dello
    studio (es. ``--config config_osd734.yaml``).
    """
    from amplicon16s.cli import main

    with pytest.raises(SystemExit) as uscita:
        main([comando])
    assert uscita.value.code == 2
