"""Scelta e scarico del file forward di una corsa a letture appaiate.

Verifica ``dati/letture_forward.py`` su depositi finti, serviti da file locali:
il file tenuto e' quello le cui letture iniziano con il primer forward, anche
quando nel deposito e' il secondo (corsa invertita); senza primer nelle letture
si tiene il primo e lo si dichiara; il file tenuto e' verificato contro la
dimensione e l'MD5 della fonte e salvato con il nome dell'elenco, senza
marcatori di coppia. Un file inverso scambiato per forward darebbe alla
pipeline l'altra estremita' dell'amplicone senza alcun errore.
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import importlib.util
import sys
from pathlib import Path

import pytest

RADICE = Path(__file__).resolve().parents[1]
SCRIPT = RADICE / "dati" / "letture_forward.py"

pytestmark = pytest.mark.skipif(not SCRIPT.is_file(), reason="dati/letture_forward.py non presente")

PRIMER = "GTGYCAGCMGCCGCGGTAA"
CON_PRIMER = "GTGCCAGCAGCCGCGGTAATACGTAGGGTGCAAGCGTTAATCGG"
SENZA_PRIMER = "TACGTAGGGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCG"
INVERSA = "GGACTACTAGGGTTTCTAATCCTGTTTGCTCCCCACGCTTTCG"


@pytest.fixture(scope="module")
def lf():
    """Il modulo, caricato dal suo file."""
    specifica = importlib.util.spec_from_file_location("letture_forward", SCRIPT)
    modulo = importlib.util.module_from_spec(specifica)
    sys.modules["letture_forward"] = modulo
    specifica.loader.exec_module(modulo)
    return modulo


def _fastq(percorso: Path, sequenza: str, quante: int = 50, altre: str | None = None) -> dict[str, str]:
    """Scrive un FASTQ compresso e ne restituisce URL, dimensione e MD5. Con
    ``altre``, una lettura su dieci porta quella sequenza."""
    with gzip.open(percorso, "wt", encoding="ascii") as file:
        for i in range(quante):
            s = altre if altre and i % 10 == 0 else sequenza
            file.write(f"@r{i}\n{s}\n+\n{'I' * len(s)}\n")
    dati = percorso.read_bytes()
    return {"url": percorso.as_uri(), "byte": str(len(dati)), "md5": hashlib.md5(dati).hexdigest()}


def _elenco(percorso: Path, corse: dict[str, tuple[dict, dict]], con_md5: bool = True) -> None:
    with open(percorso, "w", encoding="utf-8", newline="") as file:
        scrittore = csv.writer(file, delimiter="\t", lineterminator="\n")
        scrittore.writerow(["file", "url_1", "byte_1", "md5_1", "url_2", "byte_2", "md5_2", "lato"])
        for nome, (uno, due) in corse.items():
            scrittore.writerow([nome, uno["url"], uno["byte"], uno["md5"] if con_md5 else "",
                                due["url"], due["byte"], due["md5"] if con_md5 else "", ""])


def test_il_lato_forward_si_decide_dal_primer_in_testa(lf):
    """
    **Obiettivo**: Verificare che la frazione di letture con il primer in
    testa conti il primer in codici IUPAC solo all'inizio della lettura, e che
    la scelta del lato sia: il primo file se lo porta, il secondo (corsa
    invertita) se lo porta solo quello, il primo dichiarando l'assenza se non
    lo porta nessuno, un errore se lo portano entrambi.

    **Razionale scientifico e sistemistico**: Il nome del file nel deposito e'
    una dichiarazione, il primer in testa alle letture e' un fatto: dove i due
    discordano vale il fatto, e dove il fatto manca lo si dice.
    """
    assert lf.frazione_con_primer([CON_PRIMER, SENZA_PRIMER, "A" + CON_PRIMER, CON_PRIMER], PRIMER) == 0.5
    assert lf.frazione_con_primer([], PRIMER) == 0.0
    assert lf.lato_forward(0.94, 0.0) == ("1", "forward")
    assert lf.lato_forward(0.0, 0.9) == ("2", "invertita")
    assert lf.lato_forward(0.2, 0.1) == ("1", "senza primer")
    assert lf.lato_forward(0.5, 0.49) == ("1", "forward")
    with pytest.raises(lf.CoppiaAmbigua):
        lf.lato_forward(0.8, 0.7)
    # L'ultima riga, interrotta dal taglio dei byte, non e' una lettura.
    dati = gzip.compress(f"@a\n{CON_PRIMER}\n+\nIII\n@b\n{SENZA_PRIMER[:9]}".encode())
    assert lf.letture_di_testa(dati) == [CON_PRIMER]


def test_lo_scarico_tiene_il_file_forward_anche_nelle_corse_invertite(lf, tmp_path, capsys):
    """
    **Obiettivo**: Verificare su tre corse (regolare, con i due file scambiati,
    senza primer nelle letture) che in ``fastq/`` arrivi per ciascuna il file
    forward con il nome dell'elenco, identico byte per byte all'originale; che
    la corsa invertita e quella senza primer siano dichiarate; che un secondo
    lancio non scarichi nulla; che un file alterato sia riconosciuto dalla sola
    verifica; e che un MD5 della fonte diverso dal file faccia fallire la corsa
    senza lasciare il file.

    **Razionale scientifico e sistemistico**: Lo scarico e' il solo punto in
    cui un file inverso puo' entrare al posto del forward, e un file corrotto
    al posto di quello pubblicato: le due cose si verificano li', non a valle.
    """
    deposito = tmp_path / "deposito"
    deposito.mkdir()
    corse = {
        "regolare.fastq.gz": (_fastq(deposito / "a_1.fq.gz", CON_PRIMER, altre=SENZA_PRIMER),
                              _fastq(deposito / "a_2.fq.gz", INVERSA)),
        "scambiata.fastq.gz": (_fastq(deposito / "b_1.fq.gz", INVERSA),
                               _fastq(deposito / "b_2.fq.gz", CON_PRIMER)),
        "senza.fastq.gz": (_fastq(deposito / "c_1.fq.gz", SENZA_PRIMER),
                           _fastq(deposito / "c_2.fq.gz", INVERSA)),
    }
    elenco = tmp_path / "letture.tsv"
    _elenco(elenco, corse)
    lavoro = tmp_path / "lavoro"
    argomenti = ["--cartella", str(lavoro), "--pausa", "0"]
    assert lf.main(elenco, PRIMER, lavoro, "prova", [*argomenti, "--registra"]) == 0
    stampato = capsys.readouterr().out
    assert "corse invertite (tenuto il secondo file): scambiata.fastq.gz" in stampato
    assert "letture senza il primer forward in testa in 1 corse" in stampato
    attesi = {"regolare.fastq.gz": "a_1.fq.gz", "scambiata.fastq.gz": "b_2.fq.gz",
              "senza.fastq.gz": "c_1.fq.gz"}
    for nome, origine in attesi.items():
        assert (lavoro / "fastq" / nome).read_bytes() == (deposito / origine).read_bytes(), nome
    assert sorted(p.name for p in (lavoro / "fastq").iterdir()) == sorted(attesi)
    assert {r["file"]: r["lato"] for r in lf.leggi_elenco(elenco)} == {
        "regolare.fastq.gz": "1", "scambiata.fastq.gz": "2", "senza.fastq.gz": "1"}

    # Secondo lancio: tutto conforme, nessuno scarico.
    assert lf.main(elenco, PRIMER, lavoro, "prova", argomenti) == 0
    assert "file forward conformi 3, mancanti o non conformi 0" in capsys.readouterr().out
    # Un file alterato: la sola verifica lo riconosce.
    (lavoro / "fastq" / "senza.fastq.gz").write_bytes(b"alterato")
    assert lf.main(elenco, PRIMER, lavoro, "prova", [*argomenti, "--solo-verifica"]) == 1
    assert "senza.fastq.gz: non conforme" in capsys.readouterr().out

    # Un MD5 della fonte che non corrisponde al file: la corsa fallisce.
    rotto = dict(corse["regolare.fastq.gz"][0], md5="0" * 32)
    _elenco(tmp_path / "rotto.tsv", {"rotta.fastq.gz": (rotto, corse["regolare.fastq.gz"][1])})
    assert lf.main(tmp_path / "rotto.tsv", PRIMER, tmp_path / "altro", "prova",
                   ["--cartella", str(tmp_path / "altro"), "--pausa", "0", "--tentativi", "1"]) == 1
    assert not (tmp_path / "altro" / "fastq" / "rotta.fastq.gz").exists()


def test_un_md5_non_dichiarato_dalla_fonte_si_registra_al_primo_scarico(lf, tmp_path, capsys):
    """
    **Obiettivo**: Verificare che con un elenco senza MD5 lo scarico si
    rifiuti senza ``--registra``, che con ``--registra`` verifichi la
    dimensione dichiarata e scriva nell'elenco l'MD5 del file tenuto e il lato,
    e che da quel momento un file diverso sia riconosciuto.

    **Razionale scientifico e sistemistico**: Una fonte che dichiara la sola
    dimensione non permette di verificare il contenuto al primo scarico: lo si
    registra allora, e ogni scarico successivo deve dare gli stessi byte.
    """
    deposito = tmp_path / "deposito"
    deposito.mkdir()
    uno, due = _fastq(deposito / "x_1.fq.gz", CON_PRIMER), _fastq(deposito / "x_2.fq.gz", INVERSA)
    elenco = tmp_path / "letture.tsv"
    _elenco(elenco, {"campione.fastq.gz": (uno, due)}, con_md5=False)
    lavoro = tmp_path / "lavoro"
    argomenti = ["--cartella", str(lavoro), "--pausa", "0", "--tentativi", "1"]
    assert lf.main(elenco, PRIMER, lavoro, "prova", argomenti) == 1
    assert "registralo con --registra" in capsys.readouterr().out
    assert lf.main(elenco, PRIMER, lavoro, "prova", [*argomenti, "--registra"]) == 0
    (riga,) = lf.leggi_elenco(elenco)
    assert riga["md5_1"] == uno["md5"] and riga["md5_2"] == "" and riga["lato"] == "1"
    assert lf.main(elenco, PRIMER, lavoro, "prova", [*argomenti, "--solo-verifica"]) == 0

    # Una dimensione dichiarata diversa da quella del file: rifiutato.
    _elenco(tmp_path / "corto.tsv", {"corto.fastq.gz": (dict(uno, byte="10"), due)}, con_md5=False)
    assert lf.main(tmp_path / "corto.tsv", PRIMER, tmp_path / "altro", "prova",
                   ["--cartella", str(tmp_path / "altro"), "--pausa", "0", "--tentativi", "1",
                    "--registra"]) == 1
    assert not (tmp_path / "altro" / "fastq" / "corto.fastq.gz").exists()
