"""Letture inverse prese per forward: il controllo di S1 sul contenuto dei file.

Verifica che un file con le sole seconde letture di una corsa appaiata, o con
le due letture a blocchi, sia respinto anche quando le intestazioni numerano le
letture come in un file single-end e nessun segno di coppia lo rivela; e che
non siano respinti un file buono, un file con poche letture, un file con la
quota del segnale naturalmente bassa e un controllo negativo senza il
bersaglio. Il segnale di una lettura forward e' il motivo conservato dopo le
basi tagliate oppure, con il taglio, il primer in testa
(``amplicon16s/io_layer/reads.py``); la fase che giudica e' S1
(``amplicon16s/steps/s01_profile.py``), che si ferma con ``E-S1-07``.

I file sono costruiti dal dataset sintetico (``tests/sintetico.py``): le
letture inverse sono il complemento inverso delle letture generate, cioe' cio'
che il sequenziatore legge dall'altra estremita' dell'amplicone. Il motivo
conservato si ricava qui dal riferimento sintetico, non e' scritto a mano.

I due test che arrivano al profilo in R richiedono i pacchetti di calcolo:

    docker run --rm -e PYTHONPATH=/app/src -e AMPLICON16S_R_DIR=/app/R \\
      -v "$(pwd)":/app -w /app <immagine> pytest tests/test_s01_orientamento.py
"""

from __future__ import annotations

import functools
import gzip
import json
import logging
import os
from collections import Counter
from pathlib import Path

import pytest
import sintetico
from sottoinsieme import motivo_pacchetti_r_assenti

from amplicon16s.config.schema import valida
from amplicon16s.errors.catalog import CATALOGO
from amplicon16s.io_layer import reads
from amplicon16s.io_layer.reads import (
    LETTURE_MINIME_PER_LA_QUOTA,
    SEGNALE_MINIMO_DI_UNA_PARTE,
    StatisticheFile,
    conta_coppie,
    segnale_di_orientamento,
)
from amplicon16s.logging.logger import chiudi
from amplicon16s.runner.executor import Conclusione, Esecutore
from amplicon16s.runner.graph import Passo
from amplicon16s.runner.project import ProjectRun
from amplicon16s.steps.s01_profile import NOME_COPPIE, NOME_RIEPILOGO


@pytest.fixture(autouse=True)
def uscite_pulite():
    """Chiude le uscite del log prima e dopo ogni test."""
    chiudi()
    yield
    chiudi()


@functools.cache
def _sonda_r() -> str | None:
    """Perche' il profilo in R non si puo' eseguire, o ``None``."""
    return motivo_pacchetti_r_assenti("ShortRead", "Biostrings", "jsonlite")


@pytest.fixture
def con_r():
    """Richiede R con i pacchetti del profilo: salta senza, ma in CI fallisce."""
    if _sonda_r() is not None:
        if os.environ.get("AMPLICON16S_RICHIEDI_BIOC") == "1":
            pytest.fail(f"Bioconductor e' richiesto in questo ambiente: {_sonda_r()}")
        pytest.skip(_sonda_r())


_COMPLEMENTO = str.maketrans("ACGT", "TGCA")


def _inversa(sequenza: str) -> str:
    """La lettura dall'altra estremita': il complemento inverso."""
    return sequenza.translate(_COMPLEMENTO)[::-1]


def _motivo(lunghezza: int = 19) -> str:
    """Le prime basi delle sequenze del riferimento sintetico come espressione:
    la base dove quasi tutte concordano, un punto dove variano.
    """
    sequenze = [s for _, s in sintetico.sequenze_di_riferimento()]
    simboli = []
    for posizione in range(lunghezza):
        base, quante = Counter(s[posizione] for s in sequenze).most_common(1)[0]
        simboli.append(base if quante >= 0.9 * len(sequenze) else ".")
    assert sum(1 for s in simboli if s != ".") >= 8, simboli
    return "".join(simboli)


def _record(percorso: Path) -> list[tuple[str, str, str]]:
    with gzip.open(percorso, "rt", encoding="ascii") as file:
        righe = file.read().splitlines()
    return [(righe[i][1:], righe[i + 1], righe[i + 3]) for i in range(0, len(righe), 4)]


def _file_di(dataset: sintetico.DatasetSintetico, chiave: str) -> Path:
    (percorso,) = dataset.letture.glob(f"{chiave}_*.fq.gz")
    return percorso


def _riscrivi(percorso: Path, forward: int, inverse: int) -> None:
    """Riscrive il file con le prime ``forward`` letture seguite da ``inverse``
    letture inverse, numerate di seguito come in un file single-end.
    """
    record = _record(percorso)
    chiave = record[0][0].split(".")[0]
    nuovi = [(s, q) for _, s, q in record[:forward]]
    nuovi += [(_inversa(s), q[::-1]) for _, s, q in record[:inverse]]
    sintetico._scrivi_fastq(
        percorso, [(f"{chiave}.{n + 1}", s, q) for n, (s, q) in enumerate(nuovi)])


def _statistiche(con_1: int, in_1: int, con_2: int, in_2: int) -> StatisticheFile:
    return StatisticheFile("f.fq.gz", in_1 + in_2, 100, 150, 0, 0, esaurito=True,
                           segnale_per_tratto=((con_1, in_1), (con_2, in_2)))


def _config(dataset, uscita: Path, **gruppi):
    qc = {"conserved_motif": _motivo(), "min_motif_frac": 0.0, **gruppi.pop("qc", {})}
    return valida(sintetico.configurazione(dataset, uscita, qc=qc, **gruppi))


def _arresto(config) -> tuple[Conclusione, str | None, str]:
    """Esegue S0 e S1 e restituisce la conclusione, il codice e il dettaglio."""
    esito = Esecutore(ProjectRun(config), fino_a=Passo.S1).esegui()
    chiudi()
    punto = esito.punto
    return esito.conclusione, punto.codice if punto else None, (punto.dettaglio or "") if punto else ""


# --------------------------------------------------------------------------- #
# Il giudizio                                                                  #
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize(("meta", "con_quota", "respinto"), [
    # Sole letture inverse: nessuna lettura con il segnale.
    ((0, 400, 0, 400), True, "quota"),
    ((0, 250, 0, 250), True, "quota"),          # il minimo di letture, compreso
    ((0, 249, 0, 250), True, None),             # sotto il minimo non si giudica
    ((0, 400, 0, 400), False, None),            # un controllo negativo senza bersaglio
    ((1, 5000, 2, 5000), True, "quota"),        # tre su diecimila: sotto uno per mille
    ((6, 5000, 6, 5000), True, None),           # sopra uno per mille: un file povero ma forward
    # A blocchi: una meta' con il segnale, l'altra quasi nulla.
    ((180, 200, 0, 200), True, "blocchi"),
    ((0, 200, 180, 200), True, "blocchi"),      # prima le inverse, poi le forward
    ((180, 200, 0, 200), False, "blocchi"),     # anche in un controllo negativo
    ((30, 40, 2, 40), True, "blocchi"),         # il minimo di letture con il segnale
    ((29, 40, 0, 40), True, None),              # sotto il minimo non si giudica
    ((30, 40, 3, 40), True, None),              # un decimo: non e' sotto il rapporto
    # File buoni.
    ((180, 200, 175, 200), True, None),
    ((15, 22, 11, 22), True, None),             # poche letture, meta' disuguali
    ((2400, 41526, 2562, 41526), True, None),   # quota bassa ma uguale nelle due meta'
    ((14679, 52711, 12683, 52712), True, None), # una differenza reale ma piccola
])
def test_il_giudizio_distingue_inverse_blocchi_e_file_buoni(meta, con_quota, respinto):
    """
    **Obiettivo**: Verificare il giudizio sulle letture con il segnale nelle
    due meta' di un file, ai suoi estremi: respinto per la quota quando su
    almeno il minimo di letture il segnale e' sotto uno per mille (e solo dove
    la quota si giudica); respinto per la struttura a blocchi quando una meta'
    ha almeno il minimo di letture con il segnale e l'altra meno di un decimo
    della sua quota, in qualunque ordine e in qualunque classe; non respinto un
    file buono, con poche letture, con la quota bassa o con una differenza
    piccola fra le meta'.

    **Razionale scientifico e sistemistico**: Un file di letture inverse non
    porta il segnale forward se non per caso; un file a blocchi lo porta in
    una sola meta'. Le soglie del giudizio sono costanti di metodo: i casi al
    confine fissano che cosa vuol dire «almeno» e «sotto».
    """
    assert LETTURE_MINIME_PER_LA_QUOTA == 500 and SEGNALE_MINIMO_DI_UNA_PARTE == 30
    motivo = _statistiche(*meta).orientamento_inatteso(con_quota)
    if respinto is None:
        assert motivo is None
    elif respinto == "quota":
        assert motivo is not None and "quota" in motivo and "blocchi" not in motivo
    else:
        assert motivo is not None and "a blocchi" in motivo
    # Senza un segnale cercato non c'e' nulla da giudicare.
    assert StatisticheFile("f", 900, 100, 150, 0, 0).orientamento_inatteso() is None
    # Le letture alterne: stessa quota nei tratti, tutta nelle letture di un posto.
    alterno = StatisticheFile("f", 800, 100, 150, 0, 0, segnale_per_tratto=((180, 400),) * 2,
                              segnale_alterno=((360, 400), (0, 400)))
    assert "alterna letture dei due orientamenti" in alterno.orientamento_inatteso(False)


def test_il_segnale_e_il_motivo_dopo_il_taglio_oppure_il_primer_nelle_basi_tagliate():
    """
    **Obiettivo**: Verificare che con il taglio a zero il segnale sia il solo
    motivo (il primer non si cerca, perche' deve essere assente); che con il
    taglio comprenda il primer in testa, anche senza motivo; che senza motivo
    e senza un primer da cercare nelle basi tagliate non ci sia alcun segnale;
    e che la funzione valga per un primer qualunque, quindi anche per quello
    inverso di un file di letture inverse.

    **Razionale scientifico e sistemistico**: Con le letture depositate senza
    primer una regola che lo cercasse in testa respingerebbe ogni file; con il
    primer ancora nelle letture e' il segno piu' diretto dell'orientamento. Il
    giudizio riceve primer, motivo e taglio: non e' legato al verso forward.
    """
    senza_taglio = segnale_di_orientamento("GTGYCAGC", "TAC.GA", 0)
    assert (senza_taglio.primer, senza_taglio.motivo, senza_taglio.taglio) == (None, "TAC.GA", 0)
    con_taglio = segnale_di_orientamento("GTGYCAGC", None, 8)
    assert con_taglio.primer == reads.espandi_iupac("GTGYCAGC") and con_taglio.motivo is None
    assert segnale_di_orientamento("GTGYCAGC", None, 0) is None
    assert segnale_di_orientamento(None, None, 8) is None
    assert segnale_di_orientamento(None, "TAC", 8).primer is None


def test_la_scansione_conta_il_segnale_nelle_due_meta_e_unisce_primer_e_motivo(tmp_path):
    """
    **Obiettivo**: Verificare su file scritti a mano che la scansione conti
    per meta' di file le letture con il segnale; che con il taglio una lettura
    conti se ha il primer in testa oppure il motivo dopo il taglio (unione per
    lettura), e non conti se non ha ne' l'uno ne' l'altro; e che senza
    ``orientamento`` le statistiche non portino alcun conteggio.

    **Razionale scientifico e sistemistico**: Una lettura con un errore nel
    primer e il motivo integro e' forward quanto una che li ha entrambi:
    contarla assente abbasserebbe la quota dei file buoni senza alzare quella
    dei file di letture inverse, che non hanno nessuno dei due.
    """
    primer, motivo = "ACGTAC", "TTGG"
    letture = [
        primer + motivo + "AAAA",      # entrambi
        "ACCTAC" + motivo + "AAAA",    # primer con un errore, motivo integro
        primer + "TAGG" + "AAAA",      # primer integro, motivo con un errore
        "CCCCCC" + "CCCC" + "AAAA",    # nessuno dei due
    ]
    percorso = tmp_path / "f.fastq"
    percorso.write_text("".join(f"@l.{i}\n{s}\n+\n{'I' * len(s)}\n" for i, s in enumerate(letture, 1)))
    con_taglio = conta_coppie(percorso, segnale_di_orientamento(primer, motivo, len(primer)))
    assert con_taglio.segnale_alterno == ((2, 2), (1, 2))
    assert [k for k, n in con_taglio.segnale_per_tratto if n] == [1, 1, 1, 0]
    assert sum(n for _, n in con_taglio.segnale_per_tratto) == 4
    solo_motivo = conta_coppie(percorso, segnale_di_orientamento(primer, "^" + primer + motivo, 0))
    assert solo_motivo.segnale_alterno == ((1, 2), (0, 2))
    assert conta_coppie(percorso).segnale_per_tratto is None
    # Un numero dispari di letture: la seconda meta' ne ha una in piu'.
    percorso.write_text("".join(f"@l.{i}\n{s}\n+\n{'I' * len(s)}\n"
                                for i, s in enumerate(letture[:3], 1)))
    dispari = conta_coppie(percorso, segnale_di_orientamento(primer, motivo, len(primer)))
    assert dispari.segnale_alterno == ((2, 2), (1, 1))
    assert sum(k for k, _ in dispari.segnale_per_tratto) == 3
    assert len(dispari.segnale_per_tratto) == reads.TRATTI


# --------------------------------------------------------------------------- #
# La fase                                                                      #
# --------------------------------------------------------------------------- #


def test_un_file_di_sole_letture_inverse_ferma_s1_con_il_suo_nome_e_la_quota(tmp_path):
    """
    **Obiettivo**: Verificare su un dataset sintetico che, sostituite le
    letture di un campione biologico con le loro inverse (numerate come in un
    file single-end, senza alcun segno di coppia), la validazione passi, S1 si
    fermi con ``E-S1-07`` e il messaggio nomini il file, le letture, la quota
    misurata e il segnale cercato; e che il codice sia a revisione umana e dica
    che cosa fare.

    **Razionale scientifico e sistemistico**: Le seconde letture di una corsa
    appaiata, prese per forward, danno varianti dell'altra estremita'
    dell'amplicone senza alcun errore: e' il solo difetto che falsa un
    risultato in silenzio. Tutti gli altri file sono buoni: la regola non
    dipende da loro.
    """
    dataset = sintetico.genera(tmp_path / "dati", biologici=4, negativi=1, letture=600)
    inverso = _file_di(dataset, "LIB002")
    _riscrivi(inverso, forward=0, inverse=600)
    conclusione, codice, dettaglio = _arresto(_config(dataset, tmp_path / "out"))
    assert (conclusione, codice) == (Conclusione.ARRESTATA, "E-S1-07")
    assert inverso.name in dettaglio and "LIB002" in dettaglio
    assert "0 letture su 600" in dettaglio and "quota 0.0000" in dettaglio
    assert "1 file su 5" in dettaglio and _motivo() in dettaglio
    voce = CATALOGO["E-S1-07"]
    assert voce.categoria.value == "revisione_umana" and not voce.ammette_retry
    assert "letture forward" in voce.azione and "qc.conserved_motif" in voce.azione
    # Il conteggio dei segni di coppia non lo rivela: e' scritto, e dice "no".
    coppie = (tmp_path / "out" / "02_qc_profiles" / NOME_COPPIE).read_text("utf-8")
    assert "\tsi\n" not in coppie


def test_un_file_a_blocchi_ferma_s1_anche_se_s0_ne_vede_solo_la_testa(tmp_path):
    """
    **Obiettivo**: Verificare che un file con tutte le letture forward seguite
    dalle loro inverse superi S0, che ne ispeziona le sole prime letture (tutte
    forward), e fermi S1 con ``E-S1-07`` dichiarando in quale meta' sta il
    segnale; lo stesso con i blocchi in ordine opposto e in un controllo
    negativo.

    **Razionale scientifico e sistemistico**: Un file a blocchi ha circa meta'
    della quota di un file buono, come un file buono con molte sequenze fuori
    bersaglio: lo distingue il contrasto fra le sue due meta', che non dipende
    dagli altri file ne' dalla classe del campione.
    """
    dataset = sintetico.genera(tmp_path / "dati", biologici=3, negativi=1, letture=400,
                               letture_negativi=200)
    blocchi = _file_di(dataset, "LIB001")
    _riscrivi(blocchi, forward=400, inverse=400)
    config = _config(dataset, tmp_path / "out", qc={"head_reads": 100})
    conclusione, codice, dettaglio = _arresto(config)
    assert (conclusione, codice) == (Conclusione.ARRESTATA, "E-S1-07")
    assert blocchi.name in dettaglio and "nel tratto " in dettaglio
    assert "a blocchi" in dettaglio and "0 su 80" in dettaglio
    # S0 e' conclusa: nelle prime letture il file e' forward.
    assert (tmp_path / "out" / "01_input_validation" / "manifest_S0.json").is_file()

    # Prima le inverse, in un controllo negativo: stessa struttura, stesso esito.
    altro = sintetico.genera(tmp_path / "dati2", biologici=3, negativi=1, letture=400,
                             letture_negativi=200)
    negativo = _file_di(altro, "LIB004")
    record = _record(negativo)
    sintetico._scrivi_fastq(negativo, [
        (f"LIB004.{n + 1}", s, q) for n, (s, q) in enumerate(
            [(_inversa(s), q[::-1]) for _, s, q in record] + [(s, q) for _, s, q in record])])
    conclusione, codice, dettaglio = _arresto(_config(altro, tmp_path / "out2"))
    assert (conclusione, codice) == (Conclusione.ARRESTATA, "E-S1-07")
    assert negativo.name in dettaglio and "quasi assente nel tratto 1 di 10" in dettaglio


def test_un_dataset_di_sole_letture_inverse_e_respinto_anche_con_le_soglie_esplorative(tmp_path):
    """
    **Obiettivo**: Verificare che, quando ogni file contiene le sole letture
    inverse, S0 passi con ``qc.min_motif_frac`` al valore che non ferma e S1 si
    fermi comunque con ``E-S1-07``, elencando tutti i file dei campioni
    biologici e non il controllo negativo.

    **Razionale scientifico e sistemistico**: Una regola relativa alla mediana
    del dataset non vede il caso in cui tutti i file sono sbagliati allo stesso
    modo, e il controllo sulla mediana e' spento nella prima esecuzione
    esplorativa: la soglia per file e' assoluta e non dipende dalla
    configurazione delle soglie.
    """
    dataset = sintetico.genera(tmp_path / "dati", biologici=3, negativi=1, letture=500,
                               letture_negativi=500)
    for chiave in dataset.campioni:
        _riscrivi(_file_di(dataset, chiave), forward=0, inverse=500)
    conclusione, codice, dettaglio = _arresto(_config(dataset, tmp_path / "out"))
    assert (conclusione, codice) == (Conclusione.ARRESTATA, "E-S1-07")
    assert "3 file su 4" in dettaglio
    assert all(c in dettaglio for c in ("LIB001", "LIB002", "LIB003"))
    assert "LIB004" not in dettaglio


def test_senza_un_segnale_il_controllo_dichiara_di_non_applicarsi(tmp_path, caplog, con_r):
    """
    **Obiettivo**: Verificare che con il motivo nullo e il taglio a zero S1
    dichiari nel log che il controllo sulle letture inverse non si applica, non
    si fermi nemmeno con un file di sole letture inverse, e non scriva nulla al
    riguardo negli artefatti della fase.

    **Razionale scientifico e sistemistico**: Senza motivo e senza un primer
    nelle basi tagliate nulla distingue un orientamento dall'altro: il
    controllo non puo' giudicare, e chi esegue deve saperlo invece di credere
    che i file siano stati verificati.
    """
    dataset = sintetico.genera(tmp_path / "dati", biologici=3, negativi=1, letture=600)
    _riscrivi(_file_di(dataset, "LIB002"), forward=0, inverse=600)
    config = valida(sintetico.configurazione(dataset, tmp_path / "out"))
    assert config.qc.conserved_motif is None and config.filter.trimLeft == 0
    with caplog.at_level(logging.WARNING, logger="amplicon16s"):
        logging.getLogger("amplicon16s").propagate = True
        esito = Esecutore(ProjectRun(config), fino_a=Passo.S1).esegui()
    chiudi()
    assert esito.conclusione is Conclusione.COMPLETATA
    assert any("letture inverse prese per forward non si applica" in r.getMessage()
               for r in caplog.records)
    cartella = tmp_path / "out" / "02_qc_profiles"
    for artefatto in (NOME_COPPIE, NOME_RIEPILOGO, "manifest_S1.json"):
        assert "orientament" not in (cartella / artefatto).read_text("utf-8")


def test_un_dataset_buono_passa_e_il_controllo_non_lascia_traccia_negli_artefatti(tmp_path, con_r):
    """
    **Obiettivo**: Verificare che su un dataset sintetico buono, con il motivo
    dichiarato, S1 si concluda; che ``coppie.tsv`` abbia le colonne di sempre;
    che ne' il riepilogo ne' le metriche del manifesto portino un campo del
    nuovo controllo; e che il manifesto dichiari fra i parametri della fase
    primer, motivo e taglio, da cui il giudizio dipende.

    **Razionale scientifico e sistemistico**: Un controllo che passa non deve
    cambiare alcun artefatto: gli artefatti di un'esecuzione gia' pubblicata
    restano confrontabili byte per byte. I parametri dichiarati garantiscono
    che cambiare il motivo rifaccia la fase che lo usa.
    """
    dataset = sintetico.genera(tmp_path / "dati", biologici=3, negativi=1, letture=300)
    config = _config(dataset, tmp_path / "out")
    run = ProjectRun(config)
    assert Esecutore(run, fino_a=Passo.S1).esegui().conclusione is Conclusione.COMPLETATA
    chiudi()
    cartella = tmp_path / "out" / "02_qc_profiles"
    intestazione = (cartella / NOME_COPPIE).read_text("utf-8").splitlines()[0]
    assert intestazione == ("campione\tletture\tprime_di_coppia\tseconde_di_coppia\t"
                            "nomi_ripetuti\tnomi_oltre_due\tcoppie_nello_stesso_file")
    manifesto = json.loads((cartella / "manifest_S1.json").read_text("utf-8"))
    assert set(manifesto["metriche"]) == {
        "campioni", "letture", "lunghezza_minima", "lunghezza_massima", "moda",
        "qualita_mediana_minima", "letture_piu_corte_del_troncamento",
        "valori_di_qualita_distinti"}
    assert "segnale" not in (cartella / NOME_RIEPILOGO).read_text("utf-8")
    assert {"qc.primer_sequence", "qc.conserved_motif", "filter.trimLeft"} <= set(
        run.fase(Passo.S1).parametri)


@pytest.mark.parametrize("blocchi", [
    "FRFR", "FRF", "RFFR", "FFFR", "RFFFFFFFFF", "FFFFFFRRRR",
])
def test_i_blocchi_si_riconoscono_anche_disuguali_o_ripetuti_e_le_letture_alterne(
    tmp_path, blocchi
):
    """
    **Obiettivo**: Verificare sulla scansione di file costruiti dal dataset
    sintetico che siano riconosciuti i file con piu' di due blocchi, con
    blocchi di lunghezza diversa (fino a un decimo del file) e in qualunque
    ordine; e, a parte, un file in cui le letture dei due versi si alternano
    senza alcun segno di coppia nelle intestazioni.

    **Razionale scientifico e sistemistico**: Un file che concatena piu' corse
    appaiate ha un blocco per corsa e per verso, non due meta': il contrasto si
    cerca fra tratti del file, e fra letture di posto dispari e pari, non in un
    punto fissato.
    """
    dataset = sintetico.genera(tmp_path / "dati", biologici=1, negativi=0, letture=400)
    percorso = _file_di(dataset, "LIB001")
    record = _record(percorso)
    parte = len(record) // len(blocchi)
    nuovi = []
    for i, verso in enumerate(blocchi):
        for _, sequenza, qualita in record[i * parte:(i + 1) * parte]:
            nuovi.append((sequenza, qualita) if verso == "F" else (_inversa(sequenza), qualita[::-1]))
    sintetico._scrivi_fastq(percorso, [(f"LIB001.{n + 1}", s, q) for n, (s, q) in enumerate(nuovi)])
    segnale = segnale_di_orientamento(None, _motivo(), 0)
    esito = conta_coppie(percorso, segnale)
    assert esito.coppie_nello_stesso_file is None
    assert "a blocchi" in (esito.orientamento_inatteso() or "")
    if blocchi == "FRFR":
        alterni = [(s, q) if n % 2 == 0 else (_inversa(s), q[::-1])
                   for n, (_, s, q) in enumerate(record)]
        sintetico._scrivi_fastq(
            percorso, [(f"LIB001.{n + 1}", s, q) for n, (s, q) in enumerate(alterni)])
        esito = conta_coppie(percorso, segnale)
        assert esito.coppie_nello_stesso_file is None
        assert "alterna letture dei due orientamenti" in (esito.orientamento_inatteso() or "")
