"""Uscite in JSON per i programmi esterni: ``schema``, ``validate --json``,
``status --json``, per la pipeline e per le analisi ecologiche.

Verifica che i tre documenti abbiano la forma dichiarata (``versione_formato``,
un parametro, un codice, un messaggio e un'azione per ogni problema), che lo
standard output porti il solo documento, che lo stato si legga da una cartella
di uscita senza scrivervi, anche mentre un'esecuzione e' in corso, e che le
descrizioni dei parametri non entrino nel digest della configurazione.

Le esecuzioni sono costruite: una validazione vera su uno scenario sintetico,
e i file che un'esecuzione conclusa, in corso o ferma lascia nella cartella.
"""

from __future__ import annotations

import json
import os
import time
from pathlib import Path

import pytest
import yaml
from conftest import BIOLOGICO, NEGATIVO, Campione, crea_scenario, dichiarazione_minima

from amplicon16s import interfaccia
from amplicon16s.cli import main
from amplicon16s.config import defaults
from amplicon16s.config.resolve import risolvi
from amplicon16s.config.schema import Config
from amplicon16s.errors.catalog import CATALOGO
from amplicon16s.io_layer.artifacts import nome_manifesto_passo
from amplicon16s.logging.logger import chiudi
from amplicon16s.runner.executor import Conclusione, Esecutore
from amplicon16s.runner.graph import GRAFO, Passo
from amplicon16s.runner.project import ProjectRun
from amplicon16s_eco import interfaccia as interfaccia_eco
from amplicon16s_eco.cli import main as main_eco
from amplicon16s_eco.config import CONDIZIONALI, OBBLIGATORI, ConfigEco


@pytest.fixture(autouse=True)
def uscite_pulite():
    """Chiude le uscite del log prima e dopo ogni test."""
    chiudi()
    yield
    chiudi()


def _campioni() -> list[Campione]:
    """Tre biologici e un controllo negativo."""
    return [Campione(f"ERX000000{i}", f"MOD1.C{i}", BIOLOGICO) for i in (1, 2, 3)] + [
        Campione("ERX0000009", "BIANCO.1", NEGATIVO)]


def _scrivi_config(scenario, cartella: Path, **modifiche) -> Path:
    """La configurazione dello scenario in un file YAML, con le modifiche date
    (``gruppo.chiave`` -> valore; ``...`` toglie la chiave).
    """
    dati = dichiarazione_minima(scenario.config)
    for chiave, valore in modifiche.items():
        gruppo, nome = chiave.split("__")
        if valore is ...:
            dati[gruppo].pop(nome, None)
        else:
            dati.setdefault(gruppo, {})[nome] = valore
    percorso = cartella / "config.yaml"
    percorso.write_text(yaml.safe_dump(dati, allow_unicode=True), encoding="utf-8")
    return percorso


def _documento(capsys) -> dict:
    """L'unico documento JSON sullo standard output."""
    return json.loads(capsys.readouterr().out)


def _istantanea(cartella: Path) -> dict[str, tuple[int, int]]:
    """Ogni file della cartella con dimensione e istante di modifica."""
    return {
        str(p.relative_to(cartella)): (p.stat().st_size, p.stat().st_mtime_ns)
        for p in sorted(cartella.rglob("*")) if p.is_file()
    }


def test_lo_schema_descrive_ogni_parametro(capsys):
    """
    **Obiettivo**: Verificare che ``schema`` scriva un documento JSON con
    ``versione_formato``, un elenco con ogni parametro della configurazione
    (tipo, obbligatorieta', valori ammessi, specie, descrizione non vuota), che
    gli obbligatori siano quelli senza predefinito, che gli standard portino la
    fonte, e che lo schema JSON del modello accompagni l'elenco; lo stesso per
    le analisi ecologiche, con i parametri condizionali e il metodo che li
    richiede.

    **Razionale scientifico e sistemistico**: Un programma che costruisce un
    modulo di configurazione deve poterlo fare dal solo documento, senza
    leggere il codice: un parametro senza descrizione o senza i valori ammessi
    lo obbligherebbe a indovinare.
    """
    assert main(["schema"]) == 0
    documento = _documento(capsys)
    assert documento["versione_formato"] == interfaccia.VERSIONE_FORMATO == 1
    assert (documento["comando"], documento["strumento"]) == ("schema", "amplicon16s")
    parametri = {p["nome"]: p for p in documento["parametri"]}
    assert set(parametri) == set(interfaccia.chiavi_dello_schema())
    for nome, p in parametri.items():
        assert p["descrizione"] and len(p["descrizione"]) > 15, nome
        assert p["tipo"] != "valore", nome
        assert ("predefinito" in p) != p["obbligatorio"], nome
    assert {n for n, p in parametri.items() if p["specie"] == "obbligatorio"} == set(defaults.OBBLIGATORI)
    assert {n for n, p in parametri.items() if p["specie"] == "standard"} == set(
        defaults.STANDARD_DEL_METODO)
    assert all(parametri[n]["fonte"] for n in defaults.STANDARD_DEL_METODO)
    assert parametri["decontam.mode"]["valori_ammessi"] == ["aggregate", "batch"]
    assert parametri["dada.pool"]["valori_ammessi"] == [True, False, "pseudo"]
    assert parametri["ctrl.column"]["ammette_nullo"] and parametri["ctrl.column"]["obbligatorio"]
    assert parametri["io.fastq_dir"]["tipo"] == "percorso"
    assert parametri["ctrl.blank_values"]["tipo"] == "elenco di stringhe"
    assert parametri["decontam.threshold"]["vincoli"] == {"minimum": 0.0, "maximum": 1.0}
    assert parametri["run.threads"]["predefinito"] is None
    assert "$defs" in documento["json_schema"] and set(documento["gruppi"]) == set(Config.model_fields)
    assert all(d.endswith(".") and "\n" not in d for d in documento["gruppi"].values())

    assert main_eco(["schema"]) == 0
    eco = _documento(capsys)
    assert (eco["versione_formato"], eco["strumento"]) == (1, "amplicon16s_eco")
    parametri = {p["nome"]: p for p in eco["parametri"]}
    assert set(parametri) == {
        f"{g}.{n}" for g, c in ConfigEco.model_fields.items() for n in c.annotation.model_fields}
    assert all(p["descrizione"] for p in parametri.values())
    assert {n for n, p in parametri.items() if p["specie"] == "obbligatorio"} == set(OBBLIGATORI)
    assert {n for n, p in parametri.items() if p["specie"] == "condizionale"} == set(CONDIZIONALI)
    assert parametri["beta.unifrac_root"]["richiesto_da"] == {
        "parametro": "beta.distances", "valori": ["unifrac_weighted", "unifrac_unweighted"]}
    assert "aitchison" in parametri["beta.distances"]["valori_ammessi"]


def test_le_descrizioni_non_entrano_nel_digest(tmp_path, monkeypatch):
    """
    **Obiettivo**: Verificare che il digest della configurazione e l'impronta
    dei risultati non cambino quando cambia la descrizione di un parametro, e
    che nessuna descrizione compaia nella forma su cui il digest e' calcolato.

    **Razionale scientifico e sistemistico**: Il digest identifica i valori
    con cui un'esecuzione e' stata calcolata: riscrivere il testo che spiega un
    parametro non deve far risultare diversa una configurazione, ne' da rifare
    una fase conclusa.
    """
    scenario = crea_scenario(tmp_path, _campioni())
    risolta = risolvi(scenario.config)
    prima = (risolta.digest, risolta.impronta_risultati)
    canonico = json.dumps(risolta.come_mappa(), ensure_ascii=False)
    for gruppo, campo_gruppo in Config.model_fields.items():
        for nome, campo in campo_gruppo.annotation.model_fields.items():
            assert campo.description not in canonico, f"{gruppo}.{nome}"
            monkeypatch.setattr(campo, "description", "un altro testo")
    dopo = risolvi(scenario.config)
    assert (dopo.digest, dopo.impronta_risultati) == prima


def test_validate_json_riporta_parametro_codice_messaggio_e_azione(tmp_path, capsys):
    """
    **Obiettivo**: Verificare ``validate --json`` su configurazioni sbagliate e
    su una giusta: un parametro obbligatorio mancante e' un problema a se' con
    ``E-G15-10`` e il suo nome; un valore fuori dominio porta il parametro e il
    codice del controllo; una colonna inesistente ferma S0 e il problema porta
    il codice del gate, il gate e il parametro; una configurazione valida da'
    ``valida`` vero con gli avvisi; lo standard output contiene il solo
    documento, e il codice di uscita e' quello del comando senza ``--json``.

    **Razionale scientifico e sistemistico**: Il controllo e' lo stesso del
    comando per l'operatore: cambia la forma dell'esito. Un programma esterno
    deve poter mostrare ogni problema accanto al campo che lo causa, con
    l'azione del catalogo, senza interpretare un testo.
    """
    scenario = crea_scenario(tmp_path / "a", _campioni(), con_letture=True)

    mancante = _scrivi_config(scenario, tmp_path / "a", filter__truncLen=...,
                              decontam__threshold=...)
    assert main(["validate", "--config", str(mancante), "--json"]) == 3
    documento = _documento(capsys)
    assert (documento["versione_formato"], documento["comando"]) == (1, "validate")
    assert documento["valida"] is False and documento["codice_uscita"] == 3
    assert [(p["parametro"], p["codice"]) for p in documento["problemi"]] == [
        ("filter.truncLen", "E-G15-10"), ("decontam.threshold", "E-G15-10")]
    for p in documento["problemi"]:
        assert p["azione"] == CATALOGO["E-G15-10"].azione and p["parametro"] in p["messaggio"]

    fuori = _scrivi_config(scenario, tmp_path / "a", prev__min_fraction=1.5,
                           qc__warn_frac_chimeric=0.9, qc__stop_frac_chimeric=0.5)
    assert main(["validate", "--config", str(fuori), "--json"]) == 3
    problemi = {p["codice"]: p for p in _documento(capsys)["problemi"]}
    assert problemi["E-G15-04"]["parametro"] == "prev.min_fraction"
    assert set(problemi["E-G15-05"]["parametri"]) == {
        "qc.warn_frac_chimeric", "qc.stop_frac_chimeric"}
    assert all(p["azione"] and p["messaggio"] for p in problemi.values())

    colonna = _scrivi_config(scenario, tmp_path / "a", meta__sample_id_column="Colonna Assente")
    assert main(["validate", "--config", str(colonna), "--json"]) == 4
    documento = _documento(capsys)
    assert documento["valida"] is False and documento["codice_uscita"] == 4
    # Una violazione per tabella in cui la colonna manca.
    assert len(documento["problemi"]) == 2
    for problema in documento["problemi"]:
        assert (problema["codice"], problema["gate"], problema["parametro"]) == (
            "E-S0-02", "G02", "meta.sample_id_column")
        assert "Colonna Assente" in problema["messaggio"] and problema["azione"]

    chiudi()
    altro = crea_scenario(tmp_path / "b", _campioni(), con_letture=True)
    valida = _scrivi_config(altro, tmp_path / "b")
    assert main(["validate", "--config", str(valida), "--json"]) == 0
    documento = _documento(capsys)
    assert documento["valida"] is True and documento["problemi"] == []
    assert {a["codice"] for a in documento["avvisi"]} >= {"E-S0-17"}
    assert all(a["messaggio"] and a["azione"] for a in documento["avvisi"])


def test_status_legge_un_esecuzione_conclusa_in_corso_o_ferma(tmp_path, capsys):
    """
    **Obiettivo**: Verificare ``status --json`` su cartelle costruite: assente
    (``non_avviata``); con la sola prima fase conclusa (``incompleta``, con la
    fase facoltativa disattivata letta dalla configurazione registrata); con il
    file di stato dell'esecutore e un battito recente (``in_corso``, con la
    fase in esecuzione); con il battito fermo (``interrotta``); con ogni fase
    attiva conclusa (``conclusa``); ferma su un gate (``ferma``, con codice,
    fase e comando di ripresa). E che il comando non scriva nulla nella
    cartella.

    **Razionale scientifico e sistemistico**: Chi osserva un'esecuzione da un
    altro processo, o da un altro container, non deve poterla disturbare ne'
    dover ricalcolare checksum di gigabyte: legge cio' che l'esecuzione
    dichiara, e distingue un processo vivo da uno ucciso dall'eta' del battito.
    """
    assert main(["status", "--out", str(tmp_path / "assente"), "--json"]) == 0
    documento = _documento(capsys)
    assert (documento["versione_formato"], documento["comando"]) == (1, "status")
    assert documento["stato"] == "non_avviata" and documento["fasi_concluse"] == []
    assert not (tmp_path / "assente").exists()

    scenario = crea_scenario(tmp_path / "a", _campioni(), con_letture=True)
    radice = Path(scenario.config.io.out_root)
    esito = Esecutore(ProjectRun(scenario.config), fino_a=Passo.S0).esegui()
    assert esito.conclusione is Conclusione.COMPLETATA
    chiudi()
    prima = _istantanea(radice)
    assert main(["status", "--out", str(radice), "--json"]) == 0
    documento = _documento(capsys)
    assert documento["stato"] == "incompleta" and documento["fasi_concluse"] == ["S0"]
    assert documento["configurazione"]["digest"] == risolvi(scenario.config).digest
    fasi = {f["passo"]: f for f in documento["fasi"]}
    assert fasi["S0"]["stato"] == "conclusa" and fasi["S0"]["artefatti"] >= 3
    assert fasi["S9"]["stato"] == "disattivata" and fasi["S1"]["stato"] == "da_eseguire"
    assert [f["passo"] for f in documento["fasi"]] == [str(p) for p in GRAFO.ordine()]
    assert documento["esecuzione"] is None and documento["arresto"] is None
    assert _istantanea(radice) == prima

    # In corso: il file di stato dell'esecutore, con un battito recente.
    stato = radice / "99_logs" / interfaccia.NOME_STATO
    adesso = time.time()
    stato.write_text(json.dumps({
        "pid": 1, "host": "altro-container", "comando": "resume", "passo": "S4",
        "avviata": "2000-01-01T00:00:00+00:00", "fase_avviata": "2000-01-01T00:10:00+00:00",
        "intervallo_battito_s": 20.0, "aggiornato": "x", "aggiornato_epoch": adesso - 5,
    }), encoding="utf-8")
    documento = interfaccia.stato(radice, adesso)
    assert (documento["stato"], documento["fase_in_corso"]) == ("in_corso", "S4")
    assert {f["passo"]: f["stato"] for f in documento["fasi"]}["S4"] == "in_corso"
    assert documento["esecuzione"]["battito_recente"] is True
    # Lo stesso file, un'ora dopo: il processo non ha piu' dato segno.
    documento = interfaccia.stato(radice, adesso + 3600)
    assert documento["stato"] == "interrotta" and documento["fase_in_corso"] is None
    assert documento["esecuzione"]["secondi_dall_ultimo_battito"] > 3000
    stato.unlink()

    # Conclusa: ogni fase attiva ha il suo manifesto.
    manifesto = (radice / "01_input_validation" / nome_manifesto_passo("S0")).read_text()
    for passo in GRAFO.ordine():
        nodo = GRAFO.nodo(passo)
        if passo is Passo.S9:
            continue
        (radice / nodo.cartella.value).mkdir(parents=True, exist_ok=True)
        (radice / nodo.cartella.value / nome_manifesto_passo(str(passo))).write_text(manifesto)
    documento = interfaccia.stato(radice)
    assert documento["stato"] == "conclusa"
    assert documento["fasi_concluse"] == [str(p) for p in GRAFO.ordine() if p is not Passo.S9]

    # Ferma: un gate respinge, e l'arresto porta codice e punto di ripresa.
    chiudi()
    altro = crea_scenario(tmp_path / "b", _campioni(), con_letture=True)
    sbagliata = _scrivi_config(altro, tmp_path / "b", meta__sample_id_column="Colonna Assente")
    assert main(["validate", "--config", str(sbagliata)]) == 4
    capsys.readouterr()
    assert main(["status", "--config", str(sbagliata), "--json"]) == 0
    documento = _documento(capsys)
    assert documento["stato"] == "ferma"
    assert (documento["arresto"]["codice"], documento["arresto"]["passo"]) == ("E-S0-02", "S0")
    assert "resume --config" in documento["arresto"]["comando"]
    assert {f["passo"]: f["stato"] for f in documento["fasi"]}["S0"] == "da_eseguire"
    # La forma per l'operatore dice le stesse cose.
    assert main(["status", "--out", str(altro.config.io.out_root)]) == 0
    testo = capsys.readouterr().out
    assert "ferma con un punto di ripresa" in testo and "E-S0-02" in testo


def test_l_esecutore_dichiara_la_fase_in_corso_e_toglie_lo_stato_alla_fine(tmp_path, monkeypatch):
    """
    **Obiettivo**: Verificare che durante l'esecuzione di una fase il file di
    stato esista e ``stato`` riporti l'esecuzione in corso con quella fase; che
    il battito aggiorni l'istante senza altre scritture; e che alla fine,
    riuscita o ferma, il file non resti nella cartella.

    **Razionale scientifico e sistemistico**: Lo stato in corso e'
    un'informazione sull'esecuzione, non un suo artefatto: non entra in alcun
    manifesto, e un file lasciato dopo la fine farebbe credere viva
    un'esecuzione conclusa finche' il battito non invecchia.
    """
    scenario = crea_scenario(tmp_path, _campioni(), con_letture=True)
    radice = Path(scenario.config.io.out_root)
    run = ProjectRun(scenario.config)
    fase = run.fase(Passo.S0)
    visti: list[dict] = []
    originale = type(fase).calcola

    def osserva(self, contesto):
        visti.append(interfaccia.stato(radice))
        return originale(self, contesto)

    monkeypatch.setattr(type(fase), "calcola", osserva)
    assert Esecutore(run, fino_a=Passo.S0).esegui().conclusione is Conclusione.COMPLETATA
    (durante,) = visti
    assert (durante["stato"], durante["fase_in_corso"]) == ("in_corso", "S0")
    assert durante["esecuzione"]["pid"] == os.getpid()
    assert not (radice / "99_logs" / interfaccia.NOME_STATO).exists()
    assert interfaccia.stato(radice)["stato"] == "incompleta"

    # Il battito: l'istante avanza da solo, e chiudere toglie il file.
    cartella = tmp_path / "battito"
    in_corso = interfaccia.StatoInCorso(cartella, "prova", intervallo_s=0.05)
    in_corso.avvia()
    primo = json.loads((cartella / interfaccia.NOME_STATO).read_text())["aggiornato_epoch"]
    time.sleep(0.4)
    secondo = json.loads((cartella / interfaccia.NOME_STATO).read_text())["aggiornato_epoch"]
    assert secondo > primo
    in_corso.fase("S3")
    assert interfaccia.esecuzione_attiva(cartella)["passo"] == "S3"
    in_corso.chiudi()
    in_corso.chiudi()
    assert list(cartella.iterdir()) == []


def test_validate_e_status_delle_analisi_ecologiche_in_json(tmp_path, capsys):
    """
    **Obiettivo**: Verificare per le analisi ecologiche che ``validate --json``
    riporti ogni rifiuto della configurazione con il parametro, il codice
    ``E-ECO-``, il messaggio e l'azione, sul solo standard output; e che
    ``status --json`` distingua una cartella senza analisi, una con l'analisi
    in corso (file di stato con battito recente), una interrotta e una
    conclusa (manifesto presente), senza scrivervi.

    **Razionale scientifico e sistemistico**: Gli stessi tre comandi valgono
    per i due strumenti: un programma esterno li usa allo stesso modo. Un
    rifiuto delle analisi non lascia nulla nella cartella di uscita, quindi lo
    riporta il comando stesso.
    """
    config = tmp_path / "eco.yaml"
    config.write_text(yaml.safe_dump({
        "design": {"variable": "gruppo"}, "run": {"seed": 1},
        "comp": {"rank": "Genus"},
        "beta": {"distances": ["bray", "manhattan"]},
        "ord": {"methods": ["pcoa"], "distances": ["bray"]},
        "stat": {"min_group_size": 3, "significance_level": 0.05},
    }), encoding="utf-8")
    assert main_eco(["validate", "--object", str(tmp_path / "x.rds"), "--config", str(config),
                     "--json"]) == 3
    documento = _documento(capsys)
    assert (documento["versione_formato"], documento["strumento"]) == (1, "amplicon16s_eco")
    assert documento["valida"] is False and documento["codice_uscita"] == 3
    trovati = {(p["parametro"], p["codice"]) for p in documento["problemi"]}
    assert ("comp.top_n", "E-ECO-02") in trovati and ("beta.distances", "E-ECO-03") in trovati
    assert all(p["azione"] and p["messaggio"] for p in documento["problemi"])

    uscita = tmp_path / "analisi"
    assert main_eco(["status", "--out", str(uscita), "--json"]) == 0
    assert _documento(capsys)["stato"] == "non_avviata" and not uscita.exists()
    uscita.mkdir()
    adesso = time.time()
    (uscita / interfaccia.NOME_STATO).write_text(json.dumps({
        "pid": 1, "host": "h", "comando": "analisi", "passo": "analisi",
        "intervallo_battito_s": 20.0, "aggiornato_epoch": adesso - 2}), encoding="utf-8")
    assert interfaccia_eco.stato(uscita, adesso)["stato"] == "in_corso"
    assert interfaccia_eco.stato(uscita, adesso + 600)["stato"] == "interrotta"
    (uscita / interfaccia.NOME_STATO).unlink()
    (uscita / "manifest.json").write_text(json.dumps({
        "configurazione": {"digest": "abc", "file": "configurazione.json"},
        "oggetto_di_partenza": {"sha256": "sha256:0", "byte": 1},
        "file": [{"nome": "a.tsv"}, {"nome": "b.png"}]}), encoding="utf-8")
    (uscita / "riepilogo.json").write_text(json.dumps(
        {"avvisi": [{"codice": "E-ECO-11", "dettaglio": "x"}]}), encoding="utf-8")
    prima = _istantanea(uscita)
    assert main_eco(["status", "--out", str(uscita), "--json"]) == 0
    documento = _documento(capsys)
    assert (documento["stato"], documento["file"], documento["avvisi"]) == (
        "conclusa", 3, ["E-ECO-11"])
    assert documento["configurazione"]["digest"] == "abc" and documento["arresto"] is None
    assert _istantanea(uscita) == prima
