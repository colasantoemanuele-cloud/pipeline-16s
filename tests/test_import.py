"""Test segnaposto.

Verifica che i pacchetti siano importabili e che la riga di comando risponda. Serve
perché la catena di integrazione continua abbia qualcosa da eseguire davvero: una
cartella di test vuota fa uscire pytest con un codice di errore e produce un
fallimento fuorviante alla prima esecuzione.
"""

import importlib

import pytest


def test_amplicon16s_importabile():
    modulo = importlib.import_module("amplicon16s")
    assert modulo.__version__


def test_amplicon16s_eco_importabile():
    modulo = importlib.import_module("amplicon16s_eco")
    assert modulo.__version__


@pytest.mark.parametrize("comando", ["run", "resume", "validate", "report"])
def test_sottocomandi_rispondono(comando, capsys):
    from amplicon16s.cli import main

    with pytest.raises(SystemExit) as uscita:
        main([comando, "--help"])
    assert uscita.value.code == 0
    assert comando in capsys.readouterr().out


@pytest.mark.parametrize("comando", ["run", "resume", "validate", "report"])
def test_sottocomandi_chiedono_la_configurazione(comando):
    from amplicon16s.cli import main

    with pytest.raises(SystemExit) as uscita:
        main([comando])
    assert uscita.value.code == 2
