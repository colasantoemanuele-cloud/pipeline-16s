"""Configurazione assistita: il sottocomando ``init``.

Verifica che ``init`` distingua cio' che rileva sui file, cio' che propone e
cio' che lascia da scegliere; che il file scritto non sia eseguibile finche'
proposte e scelte non sono confermate, e lo diventi dopo; che i casi in cui
non c'e' una risposta sola (nessuna colonna, due colonne, letture appaiate) non
vengano risolti in silenzio; e che nel codice di ``init`` non ci siano valori
dei dataset distribuiti con il repository.

I dataset sono sintetici (``tests/sintetico.py``): una sola tabella dei
campioni, nomi di colonne e di file propri, il primer in testa alle letture.
I test che calcolano il profilo di qualita' richiedono R con ShortRead:

    docker run --rm -e PYTHONPATH=/app/src -e AMPLICON16S_R_DIR=/app/R \\
      -v "$(pwd)":/app -w /app <immagine> pytest tests/test_init.py
"""

from __future__ import annotations

import functools
import json
import os
import re
import shutil
from pathlib import Path

import pytest
import sintetico
import yaml
from init_aiuti import completa, righe_dei_parametri
from sottoinsieme import motivo_pacchetti_r_assenti

from amplicon16s import interfaccia
from amplicon16s_eco import interfaccia as interfaccia_eco
from amplicon16s.cli import main
from amplicon16s.config import defaults
from amplicon16s.config.schema import carica, obbligatori_mancanti
from amplicon16s.init import rilevamento, vocabolario
from amplicon16s.init.scrittura import SENZA_VALORE, Specie
from amplicon16s.logging.logger import chiudi
from amplicon16s.metadata.crosswalk import chiave_dalla_tabella, estrai_accession
from amplicon16s.metadata.tabelle import COLONNE_INVENTARIO
from amplicon16s.steps.s00_validate import esegui_s0

RADICE = Path(__file__).parents[1]


@pytest.fixture(autouse=True)
def uscite_pulite():
    """Chiude le uscite del log prima e dopo ogni test."""
    chiudi()
    yield
    chiudi()


@functools.cache
def _sonda_profilo() -> str | None:
    """Perche' il profilo di S1 non si puo' calcolare qui, o ``None``."""
    return motivo_pacchetti_r_assenti("ShortRead", "Biostrings", "jsonlite")


@pytest.fixture
def profilo():
    """Richiede R con i pacchetti del profilo: salta senza, ma in CI fallisce."""
    if _sonda_profilo() is not None:
        if os.environ.get("AMPLICON16S_RICHIEDI_BIOC") == "1":
            pytest.fail(f"Bioconductor e' richiesto in questo ambiente: {_sonda_profilo()}")
        pytest.skip(_sonda_profilo())


def _init(dataset, uscita: Path, *altri: str, profilo: bool = False) -> int:
    """Lancia ``init`` sul dataset sintetico dalla riga di comando."""
    argomenti = ["init", "--reads", str(dataset.letture), "--metadata", str(dataset.tabella),
                 "--output", str(uscita), *altri]
    if dataset.lotto is not None:
        argomenti += ["--batch-table", str(dataset.lotto)]
    if not profilo:
        argomenti.append("--skip-profile")
    return main(argomenti)


#: Le scelte che completano la configurazione del dataset sintetico: cio' che
#: ``init`` lascia senza valore, deciso da chi conosce i dati.
def _scelte(dataset) -> dict:
    return {
        "filter.truncLen": dataset.troncamento, "filter.truncLen_shortfall_warn": 30,
        "qc.max_frac_short_reads": 0.02, "tax.ref_name": "sintetico", "tax.ref_version": "1",
        "qc.head_reads": 200, "filt.remove_na_phylum": True, "filt.exclude_taxa": [],
        "phylo.max_seqs": 1000, "prev.apply": True, "prev.min_fraction": 0.2,
        "prev.min_count": 1, "ctrl.min_positive_pass_frac": 0.5,
        "katharoseq.target_sensitivity": 0.8, "katharoseq.min_r2": 0.5,
        "decontam.threshold": 0.1, "decontam.min_blanks": 2, "qc.max_primer_hit_frac": 0.1,
        "qc.min_motif_frac": 0.2, "qc.max_frac_lost_filter": 0.5, "qc.warn_frac_chimeric": 0.2,
        "qc.stop_frac_chimeric": 0.6, "qc.min_frac_phylum": 0.5,
        "qc.max_frac_contaminant": 0.5, "qc.min_reads_final": 100,
        "qc.min_frac_reads_retained": 0.2,
    }


def test_init_distingue_rilevato_proposto_e_da_scegliere(tmp_path, capsys):
    """
    **Obiettivo**: Verificare su un dataset sintetico con una sola tabella dei
    campioni, un file del lotto e il primer in testa alle letture che ``init``
    scriva come righe attive i soli fatti verificati (modello dei nomi, colonna
    ed espressione della chiave, chiave del file del lotto, primer e taglio
    pari alla sua lunghezza, MD5 del riferimento) e i percorsi indicati; che
    proposte e scelte siano righe in commento; che l'intestazione conti le tre
    specie come le righe; e che l'espressione rilevata estragga da ogni file
    una chiave unica, la stessa che si ricava dalla tabella.

    **Razionale scientifico e sistemistico**: Un valore attivo entra nella
    configurazione senza che nessuno lo guardi: puo' esserlo solo cio' che e'
    stato verificato su tutti i file. Tutto il resto deve chiedere un gesto.
    """
    dataset = sintetico.genera(tmp_path / "dati", biologici=5, negativi=2, letture=120,
                               primer=sintetico.PRIMER_341F, con_lotto=True)
    uscita = tmp_path / "config.yaml"
    assert _init(dataset, uscita, "--reference", str(sintetico.RIFERIMENTO),
                 "--out-root", str(tmp_path / "out")) == 0
    testo = uscita.read_text(encoding="utf-8")
    righe = righe_dei_parametri(testo)
    attive = {chiave: valore for chiave, (attiva, valore) in righe.items() if attiva}

    assert yaml.safe_load(attive["io.fastq_glob"]) == "*.fq.gz"
    assert yaml.safe_load(attive["meta.accession_column"]) == sintetico.COLONNA_FILE
    assert yaml.safe_load(attive["meta.batch_key_column"]) == sintetico.COLONNA_CHIAVE_LOTTO
    assert yaml.safe_load(attive["qc.primer_sequence"]) == sintetico.PRIMER_341F
    assert yaml.safe_load(attive["filter.trimLeft"]) == len(sintetico.PRIMER_341F)
    assert len(yaml.safe_load(attive["tax.ref_md5"])) == 32
    assert set(attive) == {
        "io.fastq_dir", "io.assay_table", "io.study_table", "io.out_root", "io.batch_table",
        "io.fastq_glob", "io.accession_regex", "meta.accession_column",
        "meta.batch_key_column", "qc.primer_sequence", "filter.trimLeft", "tax.ref_fasta",
        "tax.ref_md5",
    }
    # La chiave: unica per file, e la stessa dalla tabella.
    espressione = re.compile(yaml.safe_load(attive["io.accession_regex"]))
    file = sorted(p.name for p in dataset.letture.iterdir())
    dai_file = [estrai_accession(n, espressione) for n in file]
    assert len(set(dai_file)) == len(file)
    assert sorted(chiave_dalla_tabella(n, espressione) for n in file) == sorted(dai_file)
    assert sorted(dai_file) == sorted(dataset.campioni)  # la stessa del file del lotto

    # Proposte e scelte: in commento, con l'origine nel commento che le precede.
    for chiave in ("ctrl.column", "ctrl.blank_values", "ctrl.biological_values",
                   "meta.sample_id_column", "decontam.batch_column", "err.batch_column",
                   "qc.conserved_motif"):
        assert righe[chiave][0] is False and righe[chiave][1], chiave
    assert yaml.safe_load(righe["ctrl.column"][1]) == sintetico.COLONNA_TIPO
    assert yaml.safe_load(righe["ctrl.blank_values"][1]) == [sintetico.ETICHETTA_NEGATIVO]
    assert yaml.safe_load(righe["decontam.batch_column"][1]) == sintetico.COLONNA_PIASTRA
    assert yaml.safe_load(righe["err.batch_column"][1]) == sintetico.COLONNA_CORSA
    assert righe["decontam.threshold"] == (False, "")
    assert righe["filter.truncLen"] == (False, "")
    for marcatore in ("[RILEVATO]", "[PROPOSTO]", "[DA SCEGLIERE]", "[INDICATO]", "[STANDARD]"):
        assert marcatore in testo
    conteggio = re.search(r"Rilevati: (\d+); indicati: (\d+); proposti: (\d+); da scegliere: (\d+)", testo)
    assert conteggio is not None
    rilevati, indicati, proposti, scelte = map(int, conteggio.groups())
    assert rilevati + indicati == len(attive)
    assert rilevati + proposti + scelte >= len(defaults.OBBLIGATORI)
    stampato = capsys.readouterr().out
    assert "Che cosa init non stabilisce:" in stampato and "decontam.threshold" in stampato
    # I parametri con un predefinito non hanno una riga attiva.
    con_predefinito = {p["nome"] for p in interfaccia.schema()["parametri"] if not p["obbligatorio"]}
    assert not (set(attive) & con_predefinito) - {"io.study_table", "io.batch_table"}


def test_il_file_di_init_non_e_eseguibile_finche_non_e_completato(tmp_path, capsys):
    """
    **Obiettivo**: Verificare che sul file appena scritto ``validate`` si fermi
    con ``E-G15-10`` elencando proprio i parametri obbligatori che ``init`` ha
    lasciato in commento (proposte e scelte), e che accettate le proposte e
    scritte le scelte mancanti la validazione superi tutti i gate.

    **Razionale scientifico e sistemistico**: E' la regola «lo strumento
    propone, l'utente conferma» resa verificabile: l'elenco di ``validate`` e'
    la lista delle cose da fare, e si svuota solo con i gesti di chi configura.
    """
    dataset = sintetico.genera(tmp_path / "dati", biologici=5, negativi=2, letture=120,
                               primer=sintetico.PRIMER_341F, con_lotto=True)
    uscita = tmp_path / "config.yaml"
    assert _init(dataset, uscita, "--reference", str(sintetico.RIFERIMENTO),
                 "--out-root", str(tmp_path / "out")) == 0
    testo = uscita.read_text(encoding="utf-8")
    righe = righe_dei_parametri(testo)
    in_commento = {c for c, (attiva, _) in righe.items() if not attiva}
    mancanti = obbligatori_mancanti(yaml.safe_load(testo))
    assert set(mancanti) == in_commento & set(defaults.OBBLIGATORI) and mancanti
    capsys.readouterr()
    assert main(["validate", "--config", str(uscita)]) == 3
    stampato = capsys.readouterr().out
    assert "E-G15-10" in stampato and all(m in stampato for m in mancanti)
    assert not (tmp_path / "out").exists()

    completo = tmp_path / "completo.yaml"
    completo.write_text(completa(testo, _scelte(dataset)), encoding="utf-8")
    config = carica(completo)
    assert obbligatori_mancanti(yaml.safe_load(completo.read_text(encoding="utf-8"))) == []
    risultato = esegui_s0(config)
    assert risultato.superata, [str(v) for e in risultato.falliti for v in e.violazioni]
    assert len(risultato.inventario) == len(dataset.campioni)


def test_init_non_chiede_i_condizionali_non_pertinenti(tmp_path, capsys):
    """
    **Obiettivo**: Verificare che ``init`` scriva come ``[CONDIZIONALE]``, in
    commento, senza valore e con la condizione, i parametri la cui condizione
    e' falsa con i valori rilevati, proposti o di partenza; che non li conti
    fra i «da scegliere» e li elenchi in testa fra i non richiesti; che
    ``validate`` non li chieda; e che li chieda, con la condizione, appena chi
    configura rende vera la condizione.

    **Razionale scientifico e sistemistico**: Chiedere di scegliere un valore
    che nessuna fase usera' costringe a inventarlo, e il report lo mostrerebbe
    come una scelta. La condizione si valuta sulle proposte, che l'utente puo'
    correggere: per questo il controllo vero resta quello di ``validate``.
    """
    # Con la colonna del tipo e i soli controlli negativi: le etichette e i
    # parametri della decontaminazione sono richiesti, quelli dei controlli
    # positivi no; quelli della prevalenza dipendono da prev.apply, che init
    # lascia da scegliere.
    dataset = sintetico.genera(tmp_path / "dati", biologici=4, negativi=2, letture=80)
    uscita = tmp_path / "config.yaml"
    assert _init(dataset, uscita, "--reference", str(sintetico.RIFERIMENTO),
                 "--out-root", str(tmp_path / "out")) == 0
    stampato = capsys.readouterr().out
    testo = uscita.read_text(encoding="utf-8")
    righe = righe_dei_parametri(testo)
    non_richiesti = {
        "ctrl.blank_override_values", "ctrl.min_positives", "ctrl.min_positive_pass_frac",
        "ctrl.positive_gate", "katharoseq.target_taxon", "katharoseq.collapse_rank",
        "katharoseq.min_r2", "katharoseq.target_sensitivity", "phylo.max_seqs",
        "meta.module_regex", "prev.min_fraction", "prev.min_count",
    }
    for chiave in non_richiesti:
        assert righe[chiave] == (False, ""), chiave
        nome = chiave.split(".")[1]
        prima = testo.split(f"\n  # {nome}:")[0].rsplit("\n\n", 1)[-1]
        condizione = defaults.testo_condizione(chiave).split(" ")[0]
        assert "[CONDIZIONALE] serve solo se" in prima and condizione in prima, chiave
    for chiave in ("ctrl.blank_values", "ctrl.positive_values", "ctrl.biological_values"):
        assert righe[chiave][1], chiave
    intestazione = testo.split("\nio:")[0]
    elenco = " ".join(
        intestazione.split("# Non richiesti con le proposte attuali (condizionali):")[1]
        .split("#\n")[0].replace("#", " ").split())
    assert set(elenco.split(", ")) == non_richiesti
    da_scegliere = " ".join(
        intestazione.split("# Da scegliere:")[1].split("#\n")[0].replace("#", " ").split())
    assert not set(da_scegliere.split(", ")) & non_richiesti
    assert {"decontam.threshold", "decontam.min_blanks", "decontam.mode",
            "qc.max_frac_contaminant", "prev.apply"} <= set(da_scegliere.split(", "))
    assert f"non richiesti: {len(non_richiesti)}." in intestazione
    assert "Non richiesti con le proposte attuali (condizionali):" in stampato
    assert f"non richiesti: {len(non_richiesti)}" in stampato

    # validate non li chiede; completate le altre scelte, la configurazione e' valida.
    mancanti = obbligatori_mancanti(yaml.safe_load(testo))
    assert mancanti and not set(mancanti) & non_richiesti
    scelte = {c: v for c, v in _scelte(dataset).items() if c not in non_richiesti}
    # Il primer non e' nelle letture di questo dataset: lo dichiara chi lo conosce.
    scelte["qc.primer_sequence"] = sintetico.PRIMER_341F
    # La prevalenza resta spenta: i suoi due parametri non servono.
    scelte["prev.apply"] = False
    completo = tmp_path / "completo.yaml"
    completo.write_text(completa(testo, scelte), encoding="utf-8")
    dichiarati = yaml.safe_load(completo.read_text(encoding="utf-8"))
    assert obbligatori_mancanti(dichiarati) == []
    assert not non_richiesti & {f"{g}.{n}" for g, v in dichiarati.items() for n in (v or {})}
    assert esegui_s0(carica(completo)).superata

    # Resa vera la condizione, validate chiede il condizionale e dice perche'.
    con_filogenesi = tmp_path / "filogenesi.yaml"
    con_filogenesi.write_text(
        completa(testo, {**scelte, "phylo.enabled": True, "prev.apply": True}), encoding="utf-8")
    capsys.readouterr()
    assert main(["validate", "--config", str(con_filogenesi)]) == 3
    stampato = capsys.readouterr().out
    assert "E-G15-10" in stampato and "phylo.max_seqs perche' phylo.enabled e' vero" in stampato
    assert "prev.min_fraction perche' prev.apply e' vero" in stampato
    assert "prev.min_count perche' prev.apply e' vero" in stampato

    # Senza colonna del tipo anche le tre etichette sono non richieste.
    senza = sintetico.genera(tmp_path / "senza", biologici=3, letture=40, con_tipo=False)
    assert _init(senza, tmp_path / "senza.yaml") == 0
    righe = righe_dei_parametri((tmp_path / "senza.yaml").read_text(encoding="utf-8"))
    assert all(righe[c] == (False, "") for c in defaults.CONDIZIONALI)


def test_init_non_sceglie_in_silenzio_quando_la_risposta_non_e_una(tmp_path, capsys):
    """
    **Obiettivo**: Verificare i casi senza una risposta sola. Nessuna colonna
    che corrisponda ai file: colonna ed espressione della chiave restano da
    scegliere, senza valore, e il riepilogo lo dice. Due colonne con i nomi dei
    file: nessuna e' scelta, entrambe sono elencate. Nessun controllo e nessuna
    colonna del tipo: ``ctrl.column`` proposto nullo e le etichette non richieste.

    **Razionale scientifico e sistemistico**: Un'euristica che sceglie fra due
    risposte plausibili sbaglia meta' delle volte senza dirlo: qui la scelta
    resta a chi conosce i dati, e ``init`` porta gli elementi per farla.
    """
    # Nessuna colonna corrisponde ai file.
    senza = sintetico.genera(tmp_path / "senza", biologici=3, letture=40, con_tipo=False)
    tabella = senza.tabella.read_text(encoding="utf-8").replace(".fq.gz", "").replace("LIB", "X")
    senza.tabella.write_text(re.sub(r"X(\d+)_", r"Y\1-", tabella), encoding="utf-8")
    assert _init(senza, tmp_path / "senza.yaml") == 0
    righe = righe_dei_parametri((tmp_path / "senza.yaml").read_text(encoding="utf-8"))
    assert righe["meta.accession_column"] == (False, "") and righe["io.accession_regex"] == (False, "")
    assert "Chiave del campione non rilevata" in capsys.readouterr().out
    # Senza colonna del tipo: nessun controllo, proposto e non attivo.
    assert righe["ctrl.column"] == (False, "null")
    # Le etichette sono condizionali: senza colonna non sono richieste.
    assert righe["ctrl.biological_values"] == (False, "")

    # Due colonne con i nomi dei file.
    due = sintetico.genera(tmp_path / "due", biologici=3, letture=40)
    intestazione, *resto = due.tabella.read_text(encoding="utf-8").splitlines()
    colonne = intestazione.split("\t")
    posizione = colonne.index(sintetico.COLONNA_FILE)
    doppia = ["\t".join(colonne + ["archivio"])] + [
        r + "\t" + r.split("\t")[posizione] for r in resto]
    due.tabella.write_text("\n".join(doppia) + "\n", encoding="utf-8")
    assert _init(due, tmp_path / "due.yaml") == 0
    testo = (tmp_path / "due.yaml").read_text(encoding="utf-8")
    righe = righe_dei_parametri(testo)
    assert righe["meta.accession_column"] == (False, "")
    assert "2 colonne corrispondono uno a uno ai file" in testo
    assert f"'{sintetico.COLONNA_FILE}'" in testo and "'archivio'" in testo


def test_init_rifiuta_letture_appaiate_e_un_file_gia_esistente(tmp_path, capsys):
    """
    **Obiettivo**: Verificare che ``init`` non scriva alcuna configurazione, ed
    esca con 3, quando le letture risultano appaiate: per il nome dei file
    (``_R1`` e ``_R2``) e per il contenuto (le due letture di ogni coppia nello
    stesso file); che il messaggio dica di usare le sole forward; e che un file
    di configurazione gia' esistente non venga toccato.

    **Razionale scientifico e sistemistico**: La pipeline tratta dati
    single-end: proporre una configurazione per letture appaiate porterebbe a
    un arresto in S0, o a varianti di due regioni mescolate. E un comando che
    riscrive una configurazione gia' rivista cancellerebbe delle scelte.
    """
    per_nome = sintetico.genera(tmp_path / "nome", biologici=2, letture=30)
    for file in sorted(per_nome.letture.iterdir()):
        chiave = file.name.split("_")[0]
        shutil.copy(file, per_nome.letture / f"{chiave}_R2.fq.gz")
        file.rename(per_nome.letture / f"{chiave}_R1.fq.gz")
    assert _init(per_nome, tmp_path / "a.yaml") == 3
    stampato = capsys.readouterr().out
    assert "appaiate" in stampato and "sole letture forward" in stampato
    assert not (tmp_path / "a.yaml").exists()

    intercalato = sintetico.genera(tmp_path / "contenuto", biologici=2, letture=30)
    for file in sorted(intercalato.letture.iterdir()):
        import gzip

        record = gzip.decompress(file.read_bytes()).decode().splitlines()
        coppie = []
        for i in range(0, len(record), 4):
            nome = record[i].split(".")[0] + f":{i}"
            coppie += [f"{nome}/1", *record[i + 1:i + 4], f"{nome}/2", *record[i + 1:i + 4]]
        file.write_bytes(gzip.compress(("\n".join(coppie) + "\n").encode(), mtime=0))
    assert _init(intercalato, tmp_path / "b.yaml") == 3
    assert "non contengono le sole prime letture" in capsys.readouterr().out
    assert not (tmp_path / "b.yaml").exists()

    buono = sintetico.genera(tmp_path / "buono", biologici=2, letture=30)
    esistente = tmp_path / "c.yaml"
    esistente.write_text("io: {}\n", encoding="utf-8")
    assert _init(buono, esistente) == 3
    assert esistente.read_text(encoding="utf-8") == "io: {}\n"
    assert "esiste gia'" in capsys.readouterr().out


def test_con_due_tabelle_il_nome_del_campione_e_la_colonna_che_le_aggancia(tmp_path):
    """
    **Obiettivo**: Verificare che con una tabella di assay e una di studio
    ``init`` rilevi ``meta.sample_id_column`` come la colonna con lo stesso
    nome che aggancia ogni riga dell'assay a una sola riga dello studio; che
    con nomi diversi rilevi anche ``meta.study_sample_id_column``; che le
    etichette si propongano dalle sole righe dei campioni dell'assay; e che il
    primer dichiarato dai metadati, assente dalle letture, sia una proposta
    con il taglio proposto a zero.

    **Razionale scientifico e sistemistico**: La tabella di studio e'
    condivisa fra piu' assay: una colonna costante compare in tutte le righe e
    non aggancia nulla, e i valori di altri assay non sono etichette di questo.
    """
    dataset = sintetico.genera(tmp_path / "dati", biologici=4, negativi=1, letture=40)
    intestazione, *righe = dataset.tabella.read_text(encoding="utf-8").splitlines()
    colonne = intestazione.split("\t")
    campione, file, tipo = (colonne.index(c) for c in (
        sintetico.COLONNA_CAMPIONE, sintetico.COLONNA_FILE, sintetico.COLONNA_TIPO))
    primer = "F, 5'-GTGCCAGCAGCCGCGGTAA-3', R, 5'-GGACTACHVGGGTWTCTAAT-3'"
    assay = tmp_path / "dati" / "assay.tsv"
    assay.write_text("campione\tfile_letture\tprotocollo\tprimer\n" + "".join(
        f"{r.split(chr(9))[campione]}\t{r.split(chr(9))[file]}\tuno\t{primer}\n" for r in righe),
        encoding="utf-8")
    studio = tmp_path / "dati" / "studio.tsv"
    for nome_colonna in ("campione", "nome_nello_studio"):
        studio.write_text(f"{nome_colonna}\ttipo\tprotocollo\n" + "".join(
            f"{r.split(chr(9))[campione]}\t{r.split(chr(9))[tipo]}\tuno\n" for r in righe)
            + "ALTRO1\tcontrollo positivo (mock) di un altro assay\tuno\n", encoding="utf-8")
        uscita = tmp_path / f"{nome_colonna}.yaml"
        assert main(["init", "--reads", str(dataset.letture), "--metadata", str(assay),
                     "--study-table", str(studio), "--output", str(uscita),
                     "--skip-profile"]) == 0
        righe_file = righe_dei_parametri(uscita.read_text(encoding="utf-8"))
        assert righe_file["meta.sample_id_column"] == (True, '"campione"')
        if nome_colonna == "campione":
            assert righe_file["meta.study_sample_id_column"][0] is False
        else:
            assert righe_file["meta.study_sample_id_column"] == (True, '"nome_nello_studio"')
        assert yaml.safe_load(righe_file["ctrl.positive_values"][1]) == []
        assert yaml.safe_load(righe_file["ctrl.blank_values"][1]) == [sintetico.ETICHETTA_NEGATIVO]
        assert righe_file["qc.primer_sequence"] == (False, '"GTGCCAGCAGCCGCGGTAA"')
        assert righe_file["filter.trimLeft"] == (False, "0")
        assert righe_file["qc.conserved_motif"][0] is False and righe_file["qc.conserved_motif"][1]


def test_le_regole_di_rilevamento_sui_casi_semplici():
    """
    **Obiettivo**: Verificare le funzioni di rilevamento su casi minimi: la
    chiave ricavata quando la colonna riporta un altro nome di file che
    condivide con quello vero un solo elemento; il valore intero delimitato
    («S1» non corrisponde a «S10»); nome e versione del riferimento letti dal
    nome del file; il primer forward fra piu' sequenze dei metadati; le
    colonne del lotto riconosciute dal nome; il motivo dal consenso.

    **Razionale scientifico e sistemistico**: Ogni regola ha un caso che la
    distingue da una scorciatoia: una chiave che compare come prefisso di
    un'altra, un nome di file con due elementi numerici, due sequenze di cui
    una sola e' il primer forward.
    """
    file = [f"RUN{n}_camp{n % 3}_L001.fastq.gz" for n in (101, 102, 103)]
    righe = [{"nome": f"c{n}", "archivio": f"PRJ_X_RUN{n}_raw.fastq.gz", "gruppo": "a"}
             for n in (101, 102, 103)]
    (chiave,) = rilevamento.chiave_del_campione(file, righe, ["nome", "archivio", "gruppo"])
    assert (chiave.colonna, chiave.espressione, chiave.di_file) == ("archivio", "RUN[0-9]+", True)
    assert chiave.file_per_chiave == {f"RUN{n}": f for n, f in zip((101, 102, 103), file)}

    file = ["S1_a.fq", "S10_a.fq", "S2_a.fq"]
    righe = [{"id": "S1"}, {"id": "S10"}, {"id": "S2"}]
    (chiave,) = rilevamento.chiave_del_campione(file, righe, ["id"])
    assert chiave.file_per_chiave == {"S1": "S1_a.fq", "S10": "S10_a.fq", "S2": "S2_a.fq"}
    assert rilevamento.chiave_del_lotto(["S1", "S10", "S2"], chiave.espressione,
                                        [{"k": "S1"}, {"k": "S10"}, {"k": "S2"}, {"k": "S3"}],
                                        ["k"]) == ["k"]
    assert rilevamento.chiave_del_lotto(["S1", "S10", "S2"], chiave.espressione,
                                        [{"k": "S1"}, {"k": "S1"}, {"k": "S10"}, {"k": "S2"}],
                                        ["k"]) == []

    assert rilevamento.nome_e_versione_del_riferimento("gtdb_bac120_r220.1_train.fa.gz") == ("GTDB", None)
    assert rilevamento.nome_e_versione_del_riferimento("rdp_train_set_18.fa.gz") == ("RDP", "18")
    assert rilevamento.nome_e_versione_del_riferimento("pr2_v5.0.1_dada2.fasta") == ("PR2", "5.0.1")
    assert rilevamento.nome_e_versione_del_riferimento("riferimento.fa.gz") == (None, None)

    due = [{"p": "fwd 5'-CCTACGGGNGGCWGCAG-3'; rev 5'-GACTACHVGGGTATCTAATCC-3'"}]
    assert rilevamento.primer_dai_metadati([due]) == (
        "CCTACGGGNGGCWGCAG", ["CCTACGGGNGGCWGCAG", "GACTACHVGGGTATCTAATCC"], True)
    senza_etichetta = [{"p": "CCTACGGGNGGCWGCAG / GACTACHVGGGTATCTAATCC"}]
    assert rilevamento.primer_dai_metadati([senza_etichetta]) == (
        None, ["CCTACGGGNGGCWGCAG", "GACTACHVGGGTATCTAATCC"], False)
    assert rilevamento.primer_dai_metadati([[{"p": "515F-806R"}]]) == (None, [], False)
    # L'etichetta puo' stare nel nome della colonna; un reverse non si propone;
    # una sequenza sola senza etichetta si propone dicendo che l'indicazione manca.
    per_colonna = [{"forward_primer": "CCTACGGGNGGCWGCAG", "reverse_primer": "GACTACHVGGGTATCTAATCC"}]
    assert rilevamento.primer_dai_metadati([per_colonna])[::2] == ("CCTACGGGNGGCWGCAG", True)
    assert rilevamento.primer_dai_metadati([[{"primer": "806R GGACTACHVGGGTWTCTAAT"}]])[0] is None
    assert rilevamento.primer_dai_metadati([[{"adattatore": "AATGATACGGCGACCACCGAGATCTACAC"}]]) == (
        "AATGATACGGCGACCACCGAGATCTACAC", ["AATGATACGGCGACCACCGAGATCTACAC"], False)

    lotto = rilevamento.colonne_del_lotto(
        ["id", "primer_plate", "extraction_plate", "seq_run", "cells_per_well", "note"])
    assert lotto == {"piastra": ["extraction_plate", "primer_plate"], "corsa": ["seq_run"],
                     "cellule": ["cells_per_well"]}
    # Con i valori: la colonna della chiave e quelle con un valore per riga non
    # raggruppano, e "flowcell" non e' la colonna delle cellule.
    righe = [{"run_accession": f"R{i}", "piastra": "p1", "flowcell": "FC1",
              "cellule": "" if i else "200"} for i in range(4)]
    assert rilevamento.colonne_del_lotto(
        ["run_accession", "piastra", "flowcell", "cellule"], righe, "run_accession") == {
            "piastra": ["piastra"], "corsa": ["flowcell"], "cellule": ["cellule"]}
    # Una colonna che non lascia biologici viene dopo quella del tipo.
    tipi = rilevamento.colonne_del_tipo(
        [{"tipo": t, "esito": e} for t, e in (("tampone", "positive"), ("blank", "negative"),
                                              ("mock community", "positive"))],
        ["esito", "tipo"])
    assert [t.colonna for t in tipi] == ["tipo", "esito"] and tipi[0].biologici == ("tampone",)
    # Una forma degenerata corrisponde anche alle letture della forma originale:
    # a parita' vale la meno degenerata.
    originale = [p.nome for p in vocabolario.PRIMER_16S].index("515F (originale)")
    conteggi = tuple(100 if p.nome.startswith("515F") else 0 for p in vocabolario.PRIMER_16S)

    class _Statistiche:
        valido, letture_esaminate = True, 100

    osservato, _ = rilevamento.primer_nelle_letture(
        {"a": rilevamento.Ispezione(_Statistiche(), conteggi, (), 0)})  # type: ignore[arg-type]
    assert osservato.primer is vocabolario.PRIMER_16S[originale] and osservato.in_ogni_file

    basi = tuple({b: 100} for b in "ACGT") + ({"A": 60, "G": 38, "C": 2}, {"A": 30, "C": 30, "G": 40})
    ispezione = rilevamento.Ispezione(None, (), basi, 100)  # type: ignore[arg-type]
    assert rilevamento.motivo_dal_consenso([ispezione], 0, 6) == "ACGT[AG]"
    assert rilevamento.motivo_dal_consenso([ispezione], 4, 2) is None


def test_un_fatto_e_attivo_solo_se_vale_per_ogni_file(tmp_path, capsys):
    """
    **Obiettivo**: Verificare che il primer e il taglio siano righe attive solo
    se ogni file con letture porta il primer in testa e sono state lette tutte
    le letture: con una parte dei file senza primer, o con ``--head-reads``,
    sono proposte, con i file che non lo portano nominati; che una chiave di
    sole cifre che corrisponde ai file sia proposta e non rilevata, con
    l'abbinamento mostrato; che una tabella senza righe e un numero di letture
    non positivo siano rifiutati.

    **Razionale scientifico e sistemistico**: Con il taglio dichiarato la
    validazione non cerca piu' il primer: un taglio attivo su un dataset in cui
    alcuni file non lo portano toglierebbe loro basi vere senza alcun avviso. E
    una colonna di numeri che coincide con i nomi dei file scambierebbe i
    campioni fra loro superando ogni controllo, perche' la corrispondenza e'
    uno a uno.
    """
    import gzip

    misto = sintetico.genera(tmp_path / "misto", biologici=7, letture=60,
                             primer=sintetico.PRIMER_341F, con_tipo=False)
    for file in sorted(misto.letture.iterdir())[:3]:
        record = gzip.decompress(file.read_bytes()).decode().splitlines()
        for i in range(1, len(record), 4):
            taglio = len(sintetico.PRIMER_341F)
            record[i], record[i + 2] = record[i][taglio:], record[i + 2][taglio:]
        file.write_bytes(gzip.compress(("\n".join(record) + "\n").encode(), mtime=0))
    assert _init(misto, tmp_path / "misto.yaml") == 0
    testo = (tmp_path / "misto.yaml").read_text(encoding="utf-8")
    righe = righe_dei_parametri(testo)
    assert righe["filter.trimLeft"] == (False, "17")
    assert righe["qc.primer_sequence"] == (False, f'"{sintetico.PRIMER_341F}"')
    assert "3 file su 7 non lo portano in testa" in " ".join(testo.replace("#", " ").split())
    assert "proposti, non rilevati" in capsys.readouterr().out

    intero = sintetico.genera(tmp_path / "intero", biologici=4, letture=60,
                              primer=sintetico.PRIMER_341F, con_tipo=False)
    assert _init(intero, tmp_path / "intero.yaml") == 0
    assert righe_dei_parametri((tmp_path / "intero.yaml").read_text())["filter.trimLeft"] == (True, "17")
    assert _init(intero, tmp_path / "testa.yaml", "--head-reads", "10") == 0
    testo = (tmp_path / "testa.yaml").read_text(encoding="utf-8")
    assert righe_dei_parametri(testo)["filter.trimLeft"] == (False, "17")
    assert "sole prime letture (--head-reads)" in " ".join(testo.replace("#", " ").split())
    for sbagliato in ("0", "-5"):
        with pytest.raises(SystemExit):
            _init(intero, tmp_path / "no.yaml", "--head-reads", sbagliato)
    with pytest.raises(SystemExit):
        _init(intero, tmp_path / "no.yaml", "--threads", "-1")
    capsys.readouterr()

    # Una colonna di numeri che coincide con i nomi dei file.
    numeri = tmp_path / "numeri"
    (numeri / "letture").mkdir(parents=True)
    for n, file in enumerate(sorted(intero.letture.iterdir()), start=1):
        shutil.copy(file, numeri / "letture" / f"{n}.fq.gz")
    (numeri / "campioni.tsv").write_text(
        "campione\tprofondita_cm\n" + "".join(f"sito{n}\t{5 - n}\n" for n in range(1, 5)),
        encoding="utf-8")
    assert main(["init", "--reads", str(numeri / "letture"), "--metadata",
                 str(numeri / "campioni.tsv"), "--output", str(tmp_path / "numeri.yaml"),
                 "--skip-profile"]) == 0
    testo = (tmp_path / "numeri.yaml").read_text(encoding="utf-8")
    righe = righe_dei_parametri(testo)
    assert righe["meta.accession_column"] == (False, '"profondita_cm"')
    assert righe["io.accession_regex"][0] is False
    piatto = " ".join(testo.replace("#", " ").split())
    assert "puo' essere una coincidenza" in piatto
    assert re.search(r"abbinamento, per esempio: (\d)\.fq\.gz -> \1;", piatto)
    assert "proposta e non rilevata" in capsys.readouterr().out

    # Una tabella senza righe non e' una tabella dei campioni.
    (numeri / "vuota.tsv").write_text("campione\tfile\n", encoding="utf-8")
    assert main(["init", "--reads", str(numeri / "letture"), "--metadata",
                 str(numeri / "vuota.tsv"), "--output", str(tmp_path / "vuota.yaml"),
                 "--skip-profile"]) == 3
    assert "alcuna riga di dati" in capsys.readouterr().out
    assert not (tmp_path / "vuota.yaml").exists()

    # Le spiegazioni non hanno la forma di un parametro in commento, e un
    # numero in notazione esponenziale si rilegge come numero.
    from amplicon16s.init.scrittura import valore_yaml

    assert not [r for r in testo.splitlines() if re.match(r"^  # (?![a-zA-Z_0-9]+:)", r)]
    for numero in (1e-40, 1.0e8, 0.25, 3):
        assert yaml.safe_load(valore_yaml(numero)) == numero
        assert type(yaml.safe_load(valore_yaml(numero))) is type(numero)


def test_nessun_valore_dei_dataset_distribuiti_nel_codice_di_init(tmp_path):
    """
    **Obiettivo**: Verificare che nessun valore testuale delle configurazioni
    dei dataset distribuiti con il repository (colonne, etichette, espressioni,
    motivo, taxa, valori delle analisi ecologiche) compaia nel codice di
    ``init``, salvo il vocabolario pubblico dichiarato (il primer nell'elenco
    dei primer 16S, ciascuno con la pubblicazione, e il nome del riferimento
    nell'elenco dei riferimenti noti) e le parole dello schema; e che ogni
    valore proposto o di partenza su un dataset sintetico che coincide con
    quello di un dataset distribuito sia una coincidenza dichiarata, con la
    sua fonte.

    **Razionale scientifico e sistemistico**: Un valore di un dataset usato
    come proposta o come punto di partenza tornerebbe a essere un predefinito
    preso da quei dati, sotto un altro nome. La verifica e' sul sorgente e sul
    comportamento: un valore puo' arrivare da una costante o da una regola.
    """
    if not (RADICE / "dati" / "osd734" / "config_osd734.yaml").is_file():
        pytest.skip("la cartella dati/ non e' disponibile in questo ambiente")
    import amplicon16s.init
    import amplicon16s_eco.init

    sorgente = "\n".join(
        p.read_text(encoding="utf-8")
        for p in [*Path(amplicon16s.init.__file__).parent.glob("*.py"),
                  Path(amplicon16s_eco.init.__file__)])
    dello_schema = {
        str(v) for p in interfaccia.schema()["parametri"] for v in (p["valori_ammessi"] or [])
    } | set(COLONNE_INVENTARIO) | {
        str(v) for p in interfaccia_eco.schema()["parametri"] for v in (p["valori_ammessi"] or [])
    }
    pubblico = {p.sequenza for p in vocabolario.PRIMER_16S} | set(vocabolario.RIFERIMENTI.values())
    assert all(p.fonte and re.search(r"(19|20)\d\d", p.fonte) for p in vocabolario.PRIMER_16S)

    def stringhe(valore):
        if isinstance(valore, str):
            yield valore
        elif isinstance(valore, list):
            for elemento in valore:
                yield from stringhe(elemento)
        elif isinstance(valore, dict):
            for elemento in valore.values():
                yield from stringhe(elemento)

    pubblicate: dict[str, dict] = {}
    trovate = []
    for nome in ("osd734/config_osd734.yaml", "osd276/config_osd276.yaml",
                 "osd734/eco_osd734.yaml", "osd276/eco_osd276.yaml"):
        dati = yaml.safe_load((RADICE / "dati" / nome).read_text(encoding="utf-8"))
        pubblicate[nome] = dati
        for gruppo, valori in dati.items():
            if gruppo in ("io", "tax", "run") and "eco" not in nome:
                # I percorsi e l'immagine non sono valori da proporre; del
                # riferimento contano nome e versione, verificati sotto.
                valori = {k: v for k, v in valori.items()
                          if k in ("accession_regex", "ref_name", "ref_version")}
            for testo in stringhe(valori):
                if len(testo) < 4 or testo in dello_schema | pubblico or testo.startswith("E-S"):
                    continue
                # Come parola a se', non come parte di un nome piu' lungo.
                if re.search(rf"(?<![A-Za-z0-9_]){re.escape(testo)}(?![A-Za-z0-9_])", sorgente):
                    trovate.append((nome, gruppo, testo))
    assert trovate == []

    # Sul comportamento: le coincidenze fra cio' che init propone su un
    # dataset sintetico e i valori pubblicati devono essere dichiarate.
    coincidenze_dichiarate = {
        "meta.batch_module_column": "assenza: regola generale",
        "meta.module_regex": "assenza: regola generale",
        "ctrl.blank_override_column": "assenza: regola generale",
        "ctrl.blank_override_values": "assenza: regola generale",
        "meta.non_surface_positions": "assenza: regola generale",
        "out.batch_columns": "assenza: regola generale",
        "ctrl.positive_values": "ricavato dai dati: il dataset sintetico non ha positivi",
        "katharoseq.target_taxon": "nullo, coerente con l'assenza di controlli positivi",
        "katharoseq.cell_count_column": "assenza: regola generale",
        "phylo.enabled": "fase facoltativa: disattivata",
        "ctrl.positive_gate": "avviso e non arresto finche' la conformita' non e' nota",
        "ctrl.min_positives": "un punto piu' dei due parametri della curva",
        "katharoseq.collapse_rank": "il rango piu' fine che la pipeline assegna",
        "decontam.mode": "aggregate: l'altra modalita' si calcola comunque",
        "qc.max_zeroed_samples": "nessun campione biologico perso senza vederlo",
        "filter.trimLeft": "ricavato dai dati",
        "err.batch_column": "ricavato dai dati",
        "decontam.batch_column": "ricavato dai dati",
    }
    dataset = sintetico.genera(tmp_path / "dati", biologici=4, negativi=2, letture=60)
    import argparse

    from amplicon16s.init.comando import aggiungi_argomenti, prepara

    parser = argparse.ArgumentParser()
    aggiungi_argomenti(parser)
    voci, *_ = prepara(parser.parse_args([
        "--reads", str(dataset.letture), "--metadata", str(dataset.tabella),
        "--output", str(tmp_path / "x.yaml"), "--skip-profile"]))
    non_dichiarate = []
    for nome in ("osd734/config_osd734.yaml", "osd276/config_osd276.yaml"):
        for chiave, voce in voci.items():
            if voce.specie not in (Specie.PROPOSTO, Specie.DA_SCEGLIERE):
                continue
            if voce.valore is SENZA_VALORE:
                continue
            gruppo, parametro = chiave.split(".")
            if pubblicate[nome][gruppo].get(parametro, SENZA_VALORE) == voce.valore:
                if chiave not in coincidenze_dichiarate:
                    non_dichiarate.append((nome, chiave, voce.valore))
    assert non_dichiarate == []


def test_il_profilo_propone_troncamento_e_funzione_d_errore(profilo, tmp_path):
    """
    **Obiettivo**: Verificare, con il profilo calcolato dallo script di S1, che
    ``init`` proponga come troncamento la lunghezza che non scarta alcuna
    lettura (qui tutte le letture hanno la stessa lunghezza e la qualita' non
    scende sotto la soglia del report), con lo scarto di partenza a zero e la
    frazione di letture corte di partenza a zero; che con tre soli valori di
    qualita' proponga ``err.error_function: loess_monotono``; che le righe
    restino in commento; e che senza il profilo (``--skip-profile``) non ci sia
    alcuna proposta di troncamento.

    **Razionale scientifico e sistemistico**: Il troncamento proposto e' quello
    che il report suggerirebbe dopo S1, calcolato con la stessa funzione sulle
    stesse tabelle: chi lo accetta non trova poi nel report un valore diverso.
    """
    dataset = sintetico.genera(tmp_path / "dati", biologici=4, negativi=1, letture=200,
                               primer=sintetico.PRIMER_341F)
    assert _init(dataset, tmp_path / "a.yaml", profilo=True) == 0
    righe = righe_dei_parametri((tmp_path / "a.yaml").read_text(encoding="utf-8"))
    assert righe["filter.truncLen"] == (False, str(dataset.lunghezza))
    assert righe["filter.truncLen_shortfall_warn"] == (False, "0")
    assert righe["qc.max_frac_short_reads"] == (False, "0.0")
    assert righe["err.error_function"] == (False, '"loess"')  # il predefinito, in commento

    raggruppate = sintetico.genera(tmp_path / "raggruppate", biologici=4, letture=200,
                                   qualita_raggruppate=True)
    assert _init(raggruppate, tmp_path / "b.yaml", profilo=True) == 0
    testo = (tmp_path / "b.yaml").read_text(encoding="utf-8")
    assert righe_dei_parametri(testo)["err.error_function"] == (False, '"loess_monotono"')
    assert "3 valori di qualita' distinti (11, 25, 37)" in testo

    assert _init(dataset, tmp_path / "c.yaml") == 0
    assert righe_dei_parametri((tmp_path / "c.yaml").read_text(encoding="utf-8"))[
        "filter.truncLen"] == (False, "")


def test_la_guida_copre_ogni_parametro_obbligatorio_e_ogni_argomento_di_init():
    """
    **Obiettivo**: Verificare che la guida per un nuovo dataset nomini, nella
    sezione dei parametri, ogni parametro obbligatorio dello schema; che elenchi
    ogni argomento del sottocomando ``init`` e nessuno che non esiste; che ogni
    sottocomando che il README e la guida mostrano esista nella riga di
    comando; e che il README resti sotto le trecento righe.

    **Razionale scientifico e sistemistico**: Un parametro obbligatorio che la
    guida non spiega, o un argomento che la guida nomina e il programma non
    accetta, ferma chi configura un dataset nuovo al primo passo: la
    documentazione d'uso si tiene allineata al codice con un controllo.
    """
    import argparse

    from amplicon16s.cli import build_parser
    from amplicon16s.init.comando import aggiungi_argomenti
    from amplicon16s_eco.cli import build_parser as build_parser_eco

    if not (RADICE / "GUIDA.md").is_file():
        pytest.skip("la guida non e' disponibile in questo ambiente (non e' copiata nell'immagine)")
    guida = (RADICE / "GUIDA.md").read_text(encoding="utf-8")
    readme = (RADICE / "README.md").read_text(encoding="utf-8")
    parametri = guida[guida.index("## 6. "):]
    assert [c for c in defaults.OBBLIGATORI if f"`{c}`" not in parametri] == []

    parser = argparse.ArgumentParser()
    aggiungi_argomenti(parser)
    argomenti = {o for azione in parser._actions for o in azione.option_strings if o.startswith("--")}
    argomenti.discard("--help")
    tabella = guida[guida.index("| Argomento | Significato |"):guida.index("`init` legge tutte le letture")]
    assert set(re.findall(r"`(--[a-z-]+)", tabella)) == argomenti

    def sottocomandi(parser_radice) -> set[str]:
        (azione,) = [a for a in parser_radice._actions if isinstance(a, argparse._SubParsersAction)]
        return set(azione.choices)

    pipeline, eco = sottocomandi(build_parser()), sottocomandi(build_parser_eco())
    for testo in (guida, readme):
        for comando in re.findall(r"scripts/esegui\.py (eco )?([a-z]+)", testo):
            assert comando[1] in (eco if comando[0] else pipeline), comando
    assert {"init", "schema", "status", "validate", "run", "resume", "report"} <= pipeline
    assert {"init", "schema", "status", "validate", "run"} <= eco
    assert len(readme.splitlines()) < 300


def test_il_primer_indicato_a_mano_e_una_cartella_non_scrivibile(tmp_path, capsys):
    """
    **Obiettivo**: Verificare che ``--primer`` dichiari il primer come valore
    indicato (riga attiva), che lo si cerchi nelle letture come quelli
    dell'elenco, e che con il primer cosi' noto il motivo conservato sia
    proposto; che un primer nominato dai metadati senza la sequenza faccia
    comparire nel commento le forme che ``init`` conosce con quel nome; che una
    sequenza non valida sia respinta dalla riga di comando; e che una cartella
    di destinazione non scrivibile sia un rifiuto dichiarato, non un errore
    imprevisto a lavoro finito.

    **Razionale scientifico e sistemistico**: Il primer lo conosce chi ha
    amplificato: poterlo indicare evita di lasciare senza proposta il motivo
    della regione. E un comando che legge tutte le letture per poi fallire
    sulla scrittura fa perdere il lavoro e restituisce un codice che la
    documentazione riserva ai difetti del programma.
    """
    dataset = sintetico.genera(tmp_path / "dati", biologici=3, letture=60, con_tipo=False)
    tabella = dataset.tabella.read_text(encoding="utf-8").splitlines()
    dataset.tabella.write_text(
        "\n".join([tabella[0] + "\tprimer"] + [r + "\t341F-785R" for r in tabella[1:]]) + "\n",
        encoding="utf-8")
    assert _init(dataset, tmp_path / "a.yaml") == 0
    testo = (tmp_path / "a.yaml").read_text(encoding="utf-8")
    righe = righe_dei_parametri(testo)
    assert righe["qc.primer_sequence"] == (False, "")
    piatto = " ".join(testo.replace("#", " ").split())
    assert "i metadati nominano un primer" in piatto and "341F CCTACGGGNGGCWGCAG" in piatto
    assert righe["qc.conserved_motif"] == (False, "null")

    assert _init(dataset, tmp_path / "b.yaml", "--primer", "cctacgggnggcwgcag") == 0
    righe = righe_dei_parametri((tmp_path / "b.yaml").read_text(encoding="utf-8"))
    assert righe["qc.primer_sequence"] == (True, '"CCTACGGGNGGCWGCAG"')
    assert righe["filter.trimLeft"] == (False, "0")
    assert righe["qc.conserved_motif"][0] is False and righe["qc.conserved_motif"][1] not in ("", "null")
    with pytest.raises(SystemExit):
        _init(dataset, tmp_path / "c.yaml", "--primer", "ACGTXX")
    capsys.readouterr()

    chiusa = tmp_path / "chiusa"
    chiusa.mkdir()
    chiusa.chmod(0o555)
    try:
        if os.access(chiusa, os.W_OK):
            pytest.skip("l'utente che esegue i test scrive ovunque (amministratore)")
        assert _init(dataset, chiusa / "d.yaml") == 3
        assert "non e' scrivibile" in capsys.readouterr().out
    finally:
        chiusa.chmod(0o755)


#: Le scelte che restano dopo ``--accetta-proposte`` sul dataset sintetico:
#: cio' che non ha ne' una proposta ne' un valore di partenza.
def _scelte_dopo_l_accettazione(dataset, con_profilo: bool) -> dict:
    scelte = {
        "tax.ref_name": "sintetico", "tax.ref_version": "1", "qc.head_reads": 200,
        "filt.remove_na_phylum": True, "filt.exclude_taxa": [], "prev.apply": False,
        "decontam.threshold": 0.1, "decontam.min_blanks": 2, "qc.min_reads_final": 100,
    }
    if not con_profilo:
        scelte |= {"filter.truncLen": dataset.troncamento, "filter.truncLen_shortfall_warn": 30}
    return scelte


def test_init_senza_l_opzione_scrive_gli_stessi_byte(tmp_path, monkeypatch):
    """
    **Obiettivo**: Verificare che ``init`` senza ``--accetta-proposte`` scriva,
    su un dataset sintetico fissato e con la memoria della macchina fissata,
    un file identico byte per byte a quello registrato in
    ``tests/fixtures/init/senza_opzione.yaml`` (i soli percorsi della cartella
    di prova e del repository sono sostituiti da segnaposto).

    **Razionale scientifico e sistemistico**: L'opzione e' un gesto in piu',
    non un altro comportamento: chi non la usa deve ottenere lo stesso file,
    con ogni proposta in commento. Un confronto sui byte vede anche una riga
    d'intestazione cambiata.
    """
    from amplicon16s.init import comando

    monkeypatch.setattr(comando, "_memoria_della_macchina", lambda: 8 * 2**30)
    dataset = sintetico.genera(tmp_path / "dati", biologici=5, negativi=2, letture=120,
                               primer=sintetico.PRIMER_341F, con_lotto=True)
    uscita = tmp_path / "c.yaml"
    assert _init(dataset, uscita, "--reference", str(sintetico.RIFERIMENTO),
                 "--out-root", str(tmp_path / "out"), "--threads", "1") == 0
    testo = uscita.read_text(encoding="utf-8").replace(str(tmp_path), "<DATI>").replace(
        str(RADICE), "<RADICE>")
    atteso = (RADICE / "tests" / "fixtures" / "init" / "senza_opzione.yaml").read_text("utf-8")
    assert testo == atteso


def test_init_accetta_proposte_lascia_da_scrivere_le_sole_scelte(tmp_path, capsys):
    """
    **Obiettivo**: Verificare che con ``--accetta-proposte`` le proposte siano
    righe attive con il marcatore ``[PROPOSTO]``, i valori di partenza righe
    attive con ``[PARTENZA]``, e le soglie che fermano o avvisano senza valore
    di partenza righe attive con ``[ESPLORATIVA]`` e il valore che non ferma
    (non ``qc.min_reads_final``, che cambia il risultato); che ogni riga attiva
    sia preceduta dalla sua origine; che restino in commento le scelte senza
    valore e i condizionali non pertinenti; che l'intestazione e la console
    dicano che cosa e' stato accettato e che cosa resta; che ``validate``
    elenchi come mancanti proprio le scelte rimaste; e che, scritte quelle, la
    validazione superi tutti i gate dichiarando i controlli spenti.

    **Razionale scientifico e sistemistico**: L'opzione e' il gesto con cui
    l'utente accetta in blocco cio' che lo strumento propone. Non deve
    inventare valori: una scelta senza proposta ne' partenza resta sua.
    """
    from amplicon16s.soglie import tabella

    dataset = sintetico.genera(tmp_path / "dati", biologici=5, negativi=2, letture=120,
                               primer=sintetico.PRIMER_341F, con_lotto=True)
    uscita, senza = tmp_path / "config.yaml", tmp_path / "senza.yaml"
    comuni = ("--reference", str(sintetico.RIFERIMENTO), "--out-root", str(tmp_path / "out"))
    assert _init(dataset, senza, *comuni) == 0
    capsys.readouterr()
    assert _init(dataset, uscita, *comuni, "--accetta-proposte") == 0
    stampato = capsys.readouterr().out
    testo = uscita.read_text(encoding="utf-8")
    righe, prima = righe_dei_parametri(testo), righe_dei_parametri(senza.read_text("utf-8"))

    # Ogni riga che senza l'opzione portava un valore in commento ora e' attiva
    # con lo stesso valore, salvo i condizionali; nessun valore e' cambiato.
    da_dichiarare = {*defaults.OBBLIGATORI, *defaults.CONDIZIONALI}
    for chiave, (attiva, valore) in prima.items():
        if chiave not in da_dichiarare:
            # I parametri con un predefinito restano in commento con quello.
            assert righe[chiave] == prima[chiave], chiave
        if attiva or valore == "" or chiave not in da_dichiarare:
            continue
        assert righe[chiave] == (True, valore), chiave
    esplorative = {c for c, (attiva, _) in righe.items() if attiva and prima[c] == (False, "")}
    attese = {s.chiave: tabella.valore_esplorativo(s) for s in tabella.SOGLIE
              if s.unita is tabella.Unita.FRAZIONE}
    assert esplorative == set(attese)
    for chiave, valore in attese.items():
        assert yaml.safe_load(righe[chiave][1]) == valore
    assert righe["qc.min_reads_final"] == (False, "")
    assert yaml.safe_load(righe["qc.max_zeroed_samples"][1]) == 0

    # L'origine di ogni riga attiva e' nel commento che la precede.
    blocchi = re.findall(r"## \[([A-Z ]+)\][^\n]*\n(?:  ##[^\n]*\n)*  (\w+):", testo)
    per_nome = {nome: marcatore for marcatore, nome in blocchi}
    assert per_nome["column"] == "PROPOSTO" and per_nome["pool"] == "PARTENZA"
    assert per_nome["max_frac_lost_filter"] == "ESPLORATIVA"
    assert per_nome["max_zeroed_samples"] == "PARTENZA"
    assert "accettato con --accetta-proposte" in testo
    assert "il controllo (E-S2-02) e' spento" in testo
    # Con il primer tolto dal filtro la soglia sul primer non ha effetto, e lo dice.
    assert "la validazione non cerca il primer" in testo
    assert "scritta da init con --accetta-proposte" in testo
    assert re.search(r"Accettati con --accetta-proposte: \d+ proposti, \d+ valori di partenza,"
                     r"\n# 9 soglie esplorative", testo)
    # I valori di partenza che dichiarano un'assenza sono segnalati, da rivedere.
    assenze = re.search(r"dichiarano un'assenza \(nullo o elenco vuoto\): sono\n"
                        r"# validi[^\n]*\n((?:# [^\n]*\n)+?)#\n", testo)
    assert assenze and "meta.module_column" in assenze.group(1)
    assert "out.study_columns" in assenze.group(1) and "dada.pool" not in assenze.group(1)
    assert "dichiarano un'assenza" in stampato and "out.study_columns" in stampato
    assert "Soglie esplorative, al valore che non ferma" in stampato
    assert "Restano da scegliere" in stampato

    mancanti = obbligatori_mancanti(yaml.safe_load(testo))
    scelte = _scelte_dopo_l_accettazione(dataset, con_profilo=False)
    assert set(mancanti) == set(scelte)
    assert main(["validate", "--config", str(uscita)]) == 3
    stampato = capsys.readouterr().out
    assert all(m in stampato for m in mancanti)

    completo = tmp_path / "completo.yaml"
    completo.write_text(completa(testo, scelte, accetta=set()), encoding="utf-8")
    assert main(["validate", "--config", str(completo), "--json"]) == 0
    catturato = capsys.readouterr()
    documento = json.loads(catturato.out)
    spenti = {s["parametro"] for s in documento["controlli_spenti"]}
    # Il primer con il taglio non si giudica: non e' fra gli spenti.
    assert spenti == esplorative - {"qc.max_primer_hit_frac"}
    assert "controlli di qualita' sono spenti" in catturato.err


def test_init_accetta_proposte_con_il_profilo_non_chiede_il_troncamento(profilo, tmp_path, capsys):
    """
    **Obiettivo**: Verificare con il profilo di qualita' calcolato che
    ``--accetta-proposte`` scriva attivi il troncamento proposto e, come
    valori di partenza, lo scarto ammesso e la frazione di letture corte pari
    a zero (non il valore che non ferma: per quelle soglie ``init`` ha un
    valore ricavato dai dati); e che le scelte rimaste non comprendano alcun
    parametro del filtro.

    **Razionale scientifico e sistemistico**: Dove lo strumento ha un valore
    difendibile, la soglia parte attiva: il valore esplorativo serve solo dove
    non c'e' nulla da cui partire.
    """
    dataset = sintetico.genera(tmp_path / "dati", biologici=5, negativi=2, letture=120,
                               primer=sintetico.PRIMER_341F, con_lotto=True)
    uscita = tmp_path / "config.yaml"
    assert _init(dataset, uscita, "--reference", str(sintetico.RIFERIMENTO),
                 "--out-root", str(tmp_path / "out"), "--accetta-proposte", profilo=True) == 0
    capsys.readouterr()
    testo = uscita.read_text(encoding="utf-8")
    righe = righe_dei_parametri(testo)
    assert righe["filter.truncLen"][0] and righe["filter.truncLen_shortfall_warn"][0]
    assert righe["qc.max_frac_short_reads"] == (True, "0.0")
    assert re.search(r"## \[PARTENZA\][^\n]*\n(?:  ##[^\n]*\n)*  max_frac_short_reads:", testo)
    scelte = _scelte_dopo_l_accettazione(dataset, con_profilo=True)
    assert set(obbligatori_mancanti(yaml.safe_load(testo))) == set(scelte)
    completo = tmp_path / "completo.yaml"
    completo.write_text(completa(testo, scelte, accetta=set()), encoding="utf-8")
    assert esegui_s0(carica(completo)).superata
