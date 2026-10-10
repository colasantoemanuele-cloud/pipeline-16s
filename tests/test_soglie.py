"""Soglie che fermano o avvisano: controlli spenti e sottocomando ``soglie``.

Verifica la tabella delle soglie (:mod:`amplicon16s.soglie.tabella`): ogni
valore esplorativo e' ammesso dallo schema, e una soglia al limite o oltre e'
riconosciuta come controllo spento; la regola con cui ``soglie`` propone un
valore dalla misura di un'esecuzione conclusa; che una soglia proposta lasci
passare la misura da cui viene e faccia scattare il controllo, con il suo
codice, appena la misura peggiora oltre; i rifiuti del comando; il file e il
documento JSON che scrive.

I test senza R lavorano su una cartella di uscita costruita a mano, con i soli
artefatti che il comando legge. Quelli che eseguono la catena su un dataset
sintetico (``tests/sintetico.py``) richiedono R con i pacchetti di calcolo:

    docker run --rm -e PYTHONPATH=/app/src -e AMPLICON16S_R_DIR=/app/R \\
      -v "$(pwd)":/app -w /app <immagine> pytest tests/test_soglie.py
"""

from __future__ import annotations

import copy
import functools
import gzip
import hashlib
import json
import os
import re
import shutil
from fractions import Fraction
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import pytest
import sintetico
import yaml
from sottoinsieme import motivo_pacchetti_r_assenti

from amplicon16s.cli import main
from amplicon16s.config.resolve import registra_risolta, risolvi
from amplicon16s.config.schema import ErroreConfigurazione, valida
from amplicon16s.errors.exceptions import ErrorePipeline
from amplicon16s.gates.g01_g15 import _g09_troncamento_compatibile, _g10_primer_assente
from amplicon16s.io_layer.artifacts import nome_manifesto_passo
from amplicon16s.io_layer.reads import StatisticheFile
from amplicon16s.logging.logger import chiudi
from amplicon16s.metadata.models import Campione, ClasseCampione
from amplicon16s.report.builder import costruisci
from amplicon16s.runner.graph import GRAFO
from amplicon16s.runner.project import passi_realizzati
from amplicon16s.steps.base import impronta_parametri
from amplicon16s.soglie import comando, misure, tabella
from amplicon16s.soglie.tabella import SOGLIE, Direzione, Unita
from amplicon16s.steps.s01_profile import classi_oltre_soglia
from amplicon16s.steps.s02_filter import controlla_filtro
from amplicon16s.steps.s06_chimera import controlla_chimere
from amplicon16s.steps.s08_taxonomy import controlla_copertura

RADICE = Path(__file__).parents[1]

pytestmark = pytest.mark.xdist_group("soglie")


@pytest.fixture(autouse=True)
def uscite_pulite():
    """Chiude le uscite del log prima e dopo ogni test."""
    chiudi()
    yield
    chiudi()


@functools.cache
def _sonda_bioc() -> str | None:
    """Perche' l'ambiente R richiesto non c'e', o ``None``."""
    return motivo_pacchetti_r_assenti(
        "dada2", "ShortRead", "jsonlite", "phyloseq", "Biostrings", "decontam"
    )


@pytest.fixture(scope="module")
def bioc():
    """Richiede R con i pacchetti di tutte le fasi: salta senza, ma in CI fallisce."""
    if _sonda_bioc() is not None:
        if os.environ.get("AMPLICON16S_RICHIEDI_BIOC") == "1":
            pytest.fail(f"Bioconductor e' richiesto in questo ambiente: {_sonda_bioc()}")
        pytest.skip(_sonda_bioc())


# --------------------------------------------------------------------------- #
# Una cartella di uscita costruita a mano                                      #
# --------------------------------------------------------------------------- #

#: Le soglie al valore che non ferma, o a un valore di partenza: la
#: configurazione di una prima esecuzione esplorativa.
ESPLORATIVE: dict[str, Any] = {
    "qc.max_primer_hit_frac": 1.0, "qc.min_motif_frac": 0.0, "qc.max_frac_short_reads": 1.0,
    "filter.truncLen_shortfall_warn": 0, "qc.max_zeroed_samples": 0,
    "qc.max_frac_lost_filter": 1.0, "qc.max_asv_count": 100000,
    "qc.warn_frac_chimeric": 1.0, "qc.stop_frac_chimeric": 1.0, "qc.min_frac_phylum": 0.0,
    "qc.max_frac_contaminant": 1.0, "qc.min_frac_reads_retained": 0.0,
}

#: Le misure della cartella costruita a mano, una per soglia.
MISURE: dict[str, Any] = {
    "primer": ["0.002000", "0.010000", "0.000000"],   # per file dei biologici
    "motivo": ["0.900000", "0.800000", "0.700000"],
    "corte_s1": (30, 3000),                            # letture corte e letture, biologici
    "lunghezza_minima": 160,
    "azzerati": [],
    "perdita": 0.2,
    "varianti": 400,
    "chimerica": 0.12,
    "phylum": 0.9,
    "rimossa": 0.03,
    "negativi": 3,
    "finali": (700, 1000),
}

_BIOLOGICI = ("LIB001", "LIB002", "LIB003")
_NEGATIVI = ("LIB004", "LIB005", "LIB006")


def _dataset_fittizio() -> sintetico.DatasetSintetico:
    """La descrizione di un dataset senza i suoi file: basta alla configurazione."""
    campioni = {c: sintetico.ETICHETTA_BIOLOGICO for c in _BIOLOGICI}
    campioni |= {c: sintetico.ETICHETTA_NEGATIVO for c in _NEGATIVI}
    cartella = Path("/dati/fittizi")
    return sintetico.DatasetSintetico(
        cartella, cartella / "letture", cartella / "campioni.tsv", None, campioni, (),
        None, 137 + len(sintetico.CODA), True, {})


def _dati(uscita: Path, soglie: dict[str, Any] | None = None, **gruppi: dict[str, Any]) -> dict:
    """La configurazione del dataset fittizio, con le soglie indicate."""
    dati = sintetico.configurazione(_dataset_fittizio(), uscita, **gruppi)
    dati["qc"]["conserved_motif"] = "ACGT"
    for chiave, valore in {**ESPLORATIVE, **(soglie or {})}.items():
        gruppo, nome = chiave.split(".")
        if valore is None:
            dati[gruppo].pop(nome, None)
        else:
            dati[gruppo][nome] = valore
    return dati


def _scrivi_json(percorso: Path, documento: Any) -> None:
    percorso.parent.mkdir(parents=True, exist_ok=True)
    percorso.write_text(json.dumps(documento, indent=2), encoding="utf-8")


def _esecuzione(uscita: Path, dati: dict, **misure_diverse: Any) -> None:
    """Scrive in ``uscita`` cio' che un'esecuzione conclusa lascia e che
    ``soglie`` legge: la configurazione registrata, un manifesto per fase e gli
    artefatti con le misure.
    """
    m = {**MISURE, **misure_diverse}
    config = valida(copy.deepcopy(dati))
    risolta = risolvi(config)
    registra_risolta(risolta, uscita)
    fasi = passi_realizzati()
    troncamento = config.filter.truncLen
    metriche: dict[str, dict[str, Any]] = {
        "S1": {"letture_piu_corte_del_troncamento": {
            "biologico": {"corte": m["corte_s1"][0], "letture": m["corte_s1"][1],
                          "frazione": round(m["corte_s1"][0] / m["corte_s1"][1], 6)}}},
        "S2": {"azzerati": {"biologico": list(m["azzerati"]), "controllo_negativo": [],
                            "controllo_positivo": []},
               **({} if m["perdita"] is None else {"perdita_media_controllata": m["perdita"]})},
    }
    for passo in GRAFO.ordine():
        nodo = GRAFO.nodo(passo)
        if nodo.facoltativa:
            continue
        _scrivi_json(uscita / nodo.cartella.value / nome_manifesto_passo(str(passo)), {
            "passo": str(passo), "impronta": f"sha256:{passo}", "artefatti": [],
            "calcolata_su": {
                "configurazione": {
                    "parametri": list(fasi[passo].parametri),
                    "impronta": impronta_parametri(risolta, fasi[passo].parametri)},
                "a_monte": {str(d): f"sha256:{d}" for d in nodo.dipendenze
                            if not GRAFO.nodo(d).facoltativa}},
            "metriche": metriche.get(str(passo), {}), "degradazioni": [],
        })
    s0 = uscita / "01_input_validation"
    classi = {**{c: "biologico" for c in _BIOLOGICI}, **{c: "controllo_negativo" for c in _NEGATIVI}}
    (s0 / "crosswalk.tsv").write_text(
        "accession\tsample_name\tclasse\n"
        + "".join(f"{c}\t{c}\t{classe}\n" for c, classe in classi.items()), encoding="utf-8")
    _scrivi_json(s0 / "inventario.json", {"campioni": 6, "conteggi": {
        "biologico": 3, "controllo_positivo": 0, "controllo_negativo": 3}})
    righe = ["accession\tfile\tletture_esaminate\tfile_esaurito\tlunghezza_minima\t"
             "lunghezza_massima\tfrazione_primer\tfrazione_motivo\terrore"]
    minime = m.get("minime", [m["lunghezza_minima"]] * 3)
    esauriti = m.get("esauriti", ["si"] * 3)
    for i, campione in enumerate(_BIOLOGICI):
        righe.append(f"{campione}\t{campione}.fq.gz\t500\t{esauriti[i]}\t{minime[i]}\t170\t"
                     f"{m['primer'][i]}\t{m['motivo'][i]}\t")
    for campione in _NEGATIVI:
        # I controlli negativi non contano: valori che, contati, cambierebbero ogni misura.
        righe.append(f"{campione}\t{campione}.fq.gz\t50\tsi\t20\t170\t0.900000\t0.000000\t")
    (s0 / "letture_ispezionate.tsv").write_text("\n".join(righe) + "\n", encoding="utf-8")
    s1 = uscita / "02_qc_profiles"
    _scrivi_json(s1 / "riepilogo.json", {"lunghezza_minima": m["lunghezza_minima"]})
    (s1 / "lunghezze.tsv").write_text(
        "campione\tlunghezza\tletture\n"
        + "".join(f"{c}\t{troncamento - 1}\t{m.get('corte_per_campione', 10)}\n"
                  f"{c}\t170\t490\n" for c in _BIOLOGICI), encoding="utf-8")
    _scrivi_json(uscita / "06_seqtab" / "tabella.json", {"varianti": m["varianti"]})
    _scrivi_json(uscita / "07_chimera" / "chimere.json", {"frazione_controllata": m["chimerica"]})
    _scrivi_json(uscita / "08_taxonomy" / "riepilogo.json", {"frazione_con_phylum": {
        "biologico": {"letture": 0.99, "varianti": m["phylum"]},
        **({"controllo_positivo": {"letture": 0.99, "varianti": m["phylum_positivi"]}}
           if "phylum_positivi" in m else {}),
        "controllo_negativo": {"letture": 0.1, "varianti": 0.1}}})
    _scrivi_json(uscita / "11_controls" / "decontam_riepilogo.json", {
        "min_negativi": config.decontam.min_blanks, "confronto": {"negativi": m["negativi"]},
        "modalita_dichiarata": m.get("modalita", "aggregate"),
        "modalita": {"aggregata": {"letture_rimosse": {"biologico": m["rimossa"]}},
                     "per_piastra": {"letture_rimosse": {"biologico": m["rimossa_piastra"]}}
                     if "rimossa_piastra" in m else None}})
    _scrivi_json(uscita / "12_final" / "intermedi" / "filtri_riepilogo.json", {"letture": {
        "nonchimeric_finali": m["finali"][1], "finali": m["finali"][0]}})


def _soglie(tmp_path: Path, dati: dict, *argomenti: str, capsys=None, **misure_diverse: Any):
    """Scrive configurazione ed esecuzione, lancia ``soglie --json`` e
    restituisce codice di uscita, documento e file scritto.
    """
    uscita = Path(dati["io"]["out_root"])
    if not uscita.exists():
        _esecuzione(uscita, dati, **misure_diverse)
    file = tmp_path / "configurazione.yaml"
    file.write_text(yaml.safe_dump(dati, sort_keys=False, allow_unicode=True), encoding="utf-8")
    nuovo = tmp_path / f"nuova_{len(list(tmp_path.glob('nuova_*')))}.yaml"
    capsys.readouterr()
    codice = main(["soglie", "--config", str(file), "--output", str(nuovo), "--json", *argomenti])
    catturato = capsys.readouterr()
    return codice, json.loads(catturato.out), nuovo, catturato.err


def _per_parametro(documento: dict) -> dict[str, dict]:
    return {s["parametro"]: s for s in documento["soglie"]}


# --------------------------------------------------------------------------- #
# La tabella e il riconoscimento dei controlli spenti                          #
# --------------------------------------------------------------------------- #


def test_la_tabella_ha_le_dodici_soglie_e_i_valori_esplorativi_sono_ammessi(tmp_path):
    """
    **Obiettivo**: Verificare che la tabella elenchi dodici parametri dello
    schema, ciascuno con un codice del catalogo; che il valore esplorativo di
    ogni soglia in frazione sia l'estremo che non ferma (1 per le massime, 0
    per le minime) e che una configurazione con tutti i valori esplorativi sia
    accettata dallo schema; che le soglie in conteggio e in basi non abbiano un
    valore esplorativo.

    **Razionale scientifico e sistemistico**: Il valore scritto da ``init``
    per una prima esecuzione deve essere accettato da ``validate``, altrimenti
    il percorso si ferma al secondo passo; e deve essere proprio il limite, non
    un valore qualunque abbastanza largo.
    """
    from amplicon16s import interfaccia
    from amplicon16s.errors.catalog import CATALOGO

    chiavi = set(interfaccia.chiavi_dello_schema())
    assert len(SOGLIE) == 12 and {s.chiave for s in SOGLIE} <= chiavi
    assert all(c in CATALOGO for s in SOGLIE for c in s.codici)
    dati = _dati(tmp_path / "out")
    for voce in SOGLIE:
        valore = tabella.valore_esplorativo(voce)
        if voce.unita is Unita.FRAZIONE:
            assert valore == (1.0 if voce.direzione is Direzione.MASSIMA else 0.0)
            gruppo, nome = voce.chiave.split(".")
            dati[gruppo][nome] = valore
        else:
            assert valore is None
    config = valida(dati)
    spenti = {s.chiave for s in tabella.controlli_spenti(tabella.valori_piatti(config), 3)}
    assert spenti == {s.chiave for s in SOGLIE if s.unita is Unita.FRAZIONE}


@pytest.mark.parametrize(("chiave", "valore", "spento"), [
    ("qc.max_frac_lost_filter", 1.0, True),
    ("qc.max_frac_lost_filter", 0.9999, False),
    ("qc.min_frac_phylum", 0.0, True),
    ("qc.min_frac_phylum", 0.0001, False),
    ("qc.max_zeroed_samples", 3, True),      # quanti i campioni controllati
    ("qc.max_zeroed_samples", 999, True),    # oltre il limite, non solo uguale
    ("qc.max_zeroed_samples", 2, False),
    ("filter.truncLen_shortfall_warn", 136, True),   # filter.truncLen - 1
    ("filter.truncLen_shortfall_warn", 135, False),
    ("qc.max_asv_count", 10**12, False),     # senza limite: mai spenta
    ("qc.warn_frac_chimeric", 1.0, True),
])
def test_una_soglia_al_limite_o_oltre_e_un_controllo_spento(tmp_path, chiave, valore, spento):
    """
    **Obiettivo**: Verificare che una soglia sia riconosciuta spenta quando e'
    uguale al limite della tabella o lo supera nella sua direzione (1 come
    frazione massima, 0 come minima, un numero di campioni azzerati pari o
    superiore ai campioni controllati, lo scarto massimo ammesso), e non lo sia
    un passo prima del limite; e che la soglia sulle varianti, che non ha un
    limite, non sia mai riportata come spenta.

    **Razionale scientifico e sistemistico**: Un controllo e' spento perche'
    non puo' scattare, non perche' porta il valore che ``init`` scrive: chi
    scrive 999 campioni azzerati ha lo stesso controllo spento di chi ne
    scrive quanti sono i campioni.
    """
    soglie = {chiave: valore}
    if chiave == "qc.warn_frac_chimeric":
        soglie["qc.stop_frac_chimeric"] = 1.0
    attive = {"qc.max_primer_hit_frac": 0.1, "qc.min_motif_frac": 0.5,
              "qc.max_frac_short_reads": 0.1, "qc.max_frac_lost_filter": 0.5,
              "qc.warn_frac_chimeric": 0.2, "qc.stop_frac_chimeric": 0.6,
              "qc.min_frac_phylum": 0.5, "qc.max_frac_contaminant": 0.5,
              "qc.min_frac_reads_retained": 0.2}
    config = valida(_dati(tmp_path / "out", {**attive, **soglie}))
    spenti = {s.chiave for s in tabella.controlli_spenti(tabella.valori_piatti(config), 3)}
    assert (chiave in spenti) is spento
    assert spenti <= {chiave, "qc.stop_frac_chimeric"}


def test_un_controllo_che_non_si_esegue_non_e_elencato_fra_gli_spenti(tmp_path):
    """
    **Obiettivo**: Verificare che non siano elencate fra i controlli spenti la
    soglia sul primer con ``filter.trimLeft`` maggiore di zero, quella sul
    motivo con ``qc.conserved_motif`` nullo e la frazione di contaminanti
    omessa, anche al valore che non ferma; e che la soglia sui campioni
    azzerati non sia giudicata finche' il numero dei campioni non e' noto.

    **Razionale scientifico e sistemistico**: L'avviso deve nominare cio' che
    l'utente puo' accendere scrivendo una soglia. Un controllo che non si
    esegue qualunque sia la soglia, elencato fra gli spenti, chiederebbe un
    valore che non cambierebbe nulla.
    """
    dati = _dati(tmp_path / "out", {"qc.max_frac_contaminant": None},
                 ctrl={"blank_values": [], "positive_values": []})
    dati["filter"]["trimLeft"] = 17
    dati["qc"]["conserved_motif"] = None
    valori = tabella.valori_piatti(valida(dati))
    spenti = {s.chiave for s in tabella.controlli_spenti(valori, None)}
    assert not spenti & {"qc.max_primer_hit_frac", "qc.min_motif_frac",
                         "qc.max_frac_contaminant", "qc.max_zeroed_samples"}
    assert "qc.max_frac_lost_filter" in spenti
    # Dichiarata comunque, senza controlli negativi dichiarati: non si esegue.
    dichiarata = {**valori, "qc.max_frac_contaminant": 1.0}
    assert "qc.max_frac_contaminant" not in {
        s.chiave for s in tabella.controlli_spenti(dichiarata, 3)}
    # Con zero campioni controllati nessun campione puo' essere azzerato.
    assert "qc.max_zeroed_samples" in {
        s.chiave for s in tabella.controlli_spenti(valori, 0)}


def test_validate_e_il_report_dichiarano_i_controlli_spenti(tmp_path, capsys):
    """
    **Obiettivo**: Verificare su un dataset sintetico che ``validate`` con le
    soglie al valore che non ferma riesca (uscita 0), elenchi sulla console le
    soglie spente e le riporti nel campo ``controlli_spenti`` del documento
    JSON; che con soglie attive non dica nulla; che l'avviso non entri in
    ``gates.json``; e che il report le riporti fra le segnalazioni in apertura.

    **Razionale scientifico e sistemistico**: Una configurazione esplorativa e'
    valida ma non protegge: chi la esegue deve saperlo da ogni comando, senza
    che l'avviso cambi un artefatto di fase.
    """
    dataset = sintetico.genera(tmp_path / "dati", biologici=4, negativi=2, letture=60)
    esplorativa = {c: v for c, v in ESPLORATIVE.items() if not c.startswith("filter.")}
    dati = sintetico.configurazione(dataset, tmp_path / "out")
    attivi = sintetico.configurazione(dataset, tmp_path / "out_attivi")
    for chiave, valore in esplorativa.items():
        gruppo, nome = chiave.split(".")
        dati[gruppo][nome] = valore
    for nome, contenuto in (("esplorativa", dati), ("attiva", attivi)):
        (tmp_path / f"{nome}.yaml").write_text(yaml.safe_dump(contenuto), encoding="utf-8")

    assert main(["validate", "--config", str(tmp_path / "esplorativa.yaml"), "--json"]) == 0
    catturato = capsys.readouterr()
    documento = json.loads(catturato.out)
    spenti = {s["parametro"] for s in documento["controlli_spenti"]}
    # conserved_motif e' nullo nel dataset sintetico: il motivo non si giudica.
    assert spenti == {s.chiave for s in SOGLIE if s.unita is Unita.FRAZIONE} - {"qc.min_motif_frac"}
    assert documento["valida"] and "controlli di qualita' sono spenti" in catturato.err
    assert all(p in catturato.err for p in spenti)
    assert "spent" not in (tmp_path / "out/01_input_validation/gates.json").read_text("utf-8")
    chiudi()

    report = costruisci(tmp_path / "out")
    assert any("controlli di qualità spenti" in s for s in report.segnalazioni)
    assert {t.nome for t in report.tabelle} >= {"controlli_spenti"}

    assert main(["validate", "--config", str(tmp_path / "attiva.yaml"), "--json"]) == 0
    catturato = capsys.readouterr()
    assert json.loads(catturato.out)["controlli_spenti"] == [] and "spent" not in catturato.err
    chiudi()
    assert not any("spenti" in s for s in costruisci(tmp_path / "out_attivi").segnalazioni)


# --------------------------------------------------------------------------- #
# La regola del suggerimento                                                   #
# --------------------------------------------------------------------------- #


def _misura(valore: Any, tolleranza: Fraction = Fraction(0)) -> misure.Misura:
    return misure.Misura(None if valore is None else Fraction(str(valore)), "statistica",
                         "artefatto", tolleranza=tolleranza)


@pytest.mark.parametrize(("chiave", "misura", "usata", "esito", "valore", "spento"), [
    # Massima spenta: misura piu' margine.
    ("qc.max_frac_lost_filter", "0.2", 1.0, "proposta", 0.25, False),
    # Minima spenta: misura meno margine.
    ("qc.min_frac_phylum", "0.9", 0.0, "proposta", 0.85, False),
    # Attiva, fra misura e suggerimento: si tiene.
    ("qc.max_frac_lost_filter", "0.2", 0.22, "tenuta", 0.22, False),
    ("qc.max_frac_lost_filter", "0.2", 0.2, "tenuta", 0.2, False),
    ("qc.min_frac_phylum", "0.9", 0.88, "tenuta", 0.88, False),
    # Attiva e uguale al suggerimento: estremo compreso, si tiene.
    ("qc.max_frac_lost_filter", "0.2", 0.25, "tenuta", 0.25, False),
    ("qc.min_frac_phylum", "0.9", 0.85, "tenuta", 0.85, False),
    # Attiva ma piu' larga del suggerimento: si propone il suggerimento.
    ("qc.max_frac_lost_filter", "0.2", 0.6, "proposta", 0.25, False),
    ("qc.min_frac_phylum", "0.9", 0.5, "proposta", 0.85, False),
    # Misura e margine escono dall'intervallo: riportata al limite, spenta.
    ("qc.max_frac_lost_filter", "0.97", 1.0, "proposta", 1.0, True),
    ("qc.min_frac_phylum", "0.03", 0.0, "proposta", 0.0, True),
    # La misura e' gia' al limite.
    ("qc.max_frac_lost_filter", "1", 1.0, "proposta", 1.0, True),
    # Arrotondamento verso l'esterno a quattro decimali.
    ("qc.max_frac_lost_filter", "0.123456", 1.0, "proposta", 0.1735, False),
    ("qc.min_frac_reads_retained", "0.803667", 0.0, "proposta", 0.7536, False),
    # Avviso gia' piu' stretto della misura.
    ("qc.warn_frac_chimeric", "0.3", 0.1, "proposta", 0.35, False),
    # Arresto delle chimere: margine doppio.
    ("qc.stop_frac_chimeric", "0.3", 1.0, "proposta", 0.4, False),
    # Conteggi: nessun margine sugli azzerati, la meta' in piu' sulle varianti.
    ("qc.max_zeroed_samples", "1", 0, "proposta", 1, False),
    ("qc.max_zeroed_samples", "0", 0, "tenuta", 0, False),
    ("qc.max_zeroed_samples", "3", 3, "proposta", 3, True),
    ("qc.max_asv_count", "401", 100000, "proposta", 602, False),
    ("qc.max_asv_count", "401", 500, "tenuta", 500, False),
    # Basi: lo scarto osservato, non oltre il massimo ammesso.
    ("filter.truncLen_shortfall_warn", "5", 0, "proposta", 5, False),
    ("filter.truncLen_shortfall_warn", "500", 0, "proposta", 136, True),
])
def test_la_regola_propone_tiene_e_riporta_al_limite(chiave, misura, usata, esito, valore, spento):
    """
    **Obiettivo**: Verificare la regola del suggerimento caso per caso: misura
    piu' il margine per una soglia massima e meno per una minima; la soglia
    della configurazione tenuta quando e' attiva, rispettata dalla misura e non
    piu' larga del suggerimento; il suggerimento al posto di una soglia piu'
    larga o gia' superata; il riporto al limite, dichiarato come controllo
    spento, quando misura e margine escono dall'intervallo; l'arrotondamento
    verso l'esterno; i margini dei conteggi e delle basi.

    **Razionale scientifico e sistemistico**: Una sola regola per tipo di
    soglia, senza eccezioni nascoste: chi legge il valore proposto deve poterlo
    ricalcolare dalla misura e dal margine dichiarato.
    """
    voce = tabella.soglia(chiave)
    al_limite = {"qc.max_zeroed_samples": 3, "filter.truncLen_shortfall_warn": 136}.get(
        chiave, tabella.limite(voce, {}))
    risultato = comando.proponi(voce, _misura(misura), usata, al_limite)
    assert (risultato.esito, risultato.valore) == (esito, valore)
    assert ("CONTROLLO SPENTO" in " ".join(risultato.note)) is spento
    # Mai piu' stretto della misura.
    if voce.direzione is Direzione.MASSIMA:
        assert Fraction(str(risultato.valore)) >= Fraction(misura) or al_limite is not None
    else:
        assert Fraction(str(risultato.valore)) <= Fraction(misura)


def test_la_tolleranza_di_lettura_allarga_il_suggerimento():
    """
    **Obiettivo**: Verificare che una misura letta arrotondata (quattro
    decimali) sia peggiorata di meta' dell'ultima cifra prima del margine, per
    una soglia massima e per una minima; che la soglia di un arresto che cade
    dentro quella tolleranza sia tenuta (l'esecuzione e' arrivata in fondo: il
    valore vero la rispetta) e quella di un avviso no.

    **Razionale scientifico e sistemistico**: Il controllo confronta il valore
    non arrotondato: una soglia calcolata sul valore scritto potrebbe stare
    sotto quello vero, e la ripresa sugli stessi dati si fermerebbe per una
    soglia appena proposta.
    """
    mezza = Fraction(1, 20_000)
    massima = comando.proponi(tabella.soglia("qc.max_frac_lost_filter"),
                              _misura("0.2", mezza), 1.0, 1)
    assert massima.valore == 0.2501
    # Un controllo che ferma non e' scattato: una soglia usata che cade nella
    # tolleranza e' rispettata dal valore vero, e si tiene.
    tenuta = comando.proponi(tabella.soglia("qc.max_frac_lost_filter"),
                             _misura("0.2", mezza), 0.20002, 1)
    assert tenuta.esito == "tenuta" and tenuta.valore == 0.20002
    assert not tenuta.note
    # Un avviso puo' essere scattato: la soglia dentro la tolleranza non si tiene.
    avviso = comando.proponi(tabella.soglia("qc.warn_frac_chimeric"),
                             _misura("0.2", mezza), 0.20002, 1)
    assert avviso.esito == "proposta" and avviso.valore == 0.2501
    minima = comando.proponi(tabella.soglia("qc.min_motif_frac"), _misura("0.9", mezza), 0.0, 0)
    assert minima.valore == 0.8499


def test_l_avviso_delle_chimere_non_supera_mai_l_arresto(tmp_path, capsys):
    """
    **Obiettivo**: Verificare che, quando l'arresto della configurazione viene
    tenuto ed e' sotto la misura piu' il margine dell'avviso, l'avviso proposto
    sia riportato al valore dell'arresto e la nota lo dichiari; che nel caso
    ordinario avviso e arresto abbiano i due margini; e che il file scritto
    sia una configurazione valida in entrambi i casi.

    **Razionale scientifico e sistemistico**: Lo schema respinge un avviso
    superiore all'arresto: una coppia proposta in modo indipendente darebbe un
    file che ``validate`` rifiuta.
    """
    codice, documento, _, _ = _soglie(
        tmp_path / "a", _dati(tmp_path / "a/out"), "--accetta-proposte", capsys=capsys)
    soglie = _per_parametro(documento)
    assert codice == 0
    assert soglie["qc.warn_frac_chimeric"]["valore"] == 0.17
    assert soglie["qc.stop_frac_chimeric"]["valore"] == 0.22

    stretta = _dati(tmp_path / "b/out", {"qc.stop_frac_chimeric": 0.14,
                                         "qc.warn_frac_chimeric": 0.05})
    codice, documento, nuovo, _ = _soglie(tmp_path / "b", stretta, "--accetta-proposte",
                                          capsys=capsys)
    soglie = _per_parametro(documento)
    assert codice == 0
    assert soglie["qc.stop_frac_chimeric"]["esito"] == "tenuta"
    assert soglie["qc.warn_frac_chimeric"]["valore"] == 0.14
    assert any("riportato al valore dell'arresto" in n for n in soglie["qc.warn_frac_chimeric"]["note"])
    qc = valida(yaml.safe_load(nuovo.read_text(encoding="utf-8"))).qc
    assert qc.warn_frac_chimeric == qc.stop_frac_chimeric == 0.14


# --------------------------------------------------------------------------- #
# Il comando: misure, file, documento, rifiuti                                 #
# --------------------------------------------------------------------------- #


def test_soglie_legge_la_statistica_che_ogni_controllo_confronta(tmp_path, capsys):
    """
    **Obiettivo**: Verificare su una cartella costruita a mano che ``soglie``
    legga per ogni soglia la statistica del suo controllo e non un'altra: il
    massimo fra i file dei biologici per il primer (non i controlli negativi,
    non la media), la mediana per il motivo, la frazione di classe per le
    letture corte, lo scarto in basi, i campioni azzerati delle classi
    controllate, la perdita media, le varianti, la frazione chimerica, il
    minimo fra le classi controllate per il phylum, la frazione rimossa dalla
    modalita' dichiarata, le letture finali su quelle senza chimere; e che il
    documento riporti per ciascuna misura, artefatto, valore della
    configurazione, valore proposto e regola.

    **Razionale scientifico e sistemistico**: Una soglia proposta su una
    statistica diversa da quella del controllo non protegge: la mediana non
    dice nulla di un controllo che giudica il file peggiore.
    """
    codice, documento, _, errore = _soglie(tmp_path, _dati(tmp_path / "out"), capsys=capsys)
    assert codice == 0 and documento["comando"] == "soglie" and documento["versione_formato"] == 1
    soglie = _per_parametro(documento)
    attese = {
        "qc.max_primer_hit_frac": (0.01, 0.06),
        "qc.min_motif_frac": (0.8, 0.75),
        "qc.max_frac_short_reads": (0.01, 0.06),
        "filter.truncLen_shortfall_warn": (160 - 137, 23),
        "qc.max_zeroed_samples": (0, 0),
        "qc.max_frac_lost_filter": (0.2, 0.2501),
        "qc.max_asv_count": (400, 600),
        "qc.warn_frac_chimeric": (0.12, 0.17),
        "qc.stop_frac_chimeric": (0.12, 0.22),
        "qc.min_frac_phylum": (0.9, 0.85),
        "qc.max_frac_contaminant": (0.03, 0.08),
        "qc.min_frac_reads_retained": (0.7, 0.65),
    }
    assert {c: (s["misura"], s["valore"]) for c, s in soglie.items()} == attese
    for voce in soglie.values():
        assert voce["origine_della_misura"] and voce["statistica"] and voce["regola"]
        assert voce["valore_configurazione"] == ESPLORATIVE[voce["parametro"]]
    assert soglie["qc.max_zeroed_samples"]["esito"] == "tenuta"
    assert documento["controlli_spenti"] == []
    # La console riporta la misura accanto al valore, e l'avvertenza.
    assert "misura 0.12" in errore and "proposta 0.17" in errore
    assert "non giudica se la misura e' buona" in errore
    assert "non giudica" in documento["avvertenza"]

    # Con due classi controllate conta la peggiore; con la modalita' per
    # piastra la frazione rimossa e' quella della modalita' dichiarata.
    dati = _dati(tmp_path / "due/out", decontam={"mode": "batch", "batch_column": "piastra"},
                 io={"batch_table": "/dati/fittizi/lotto.tsv"},
                 meta={"batch_key_column": "libreria"}, err={"batch_column": None})
    codice, documento, _, _ = _soglie(
        tmp_path / "due", dati, capsys=capsys, phylum_positivi=0.7, modalita="batch",
        rimossa_piastra=0.2)
    assert codice == 0, documento
    soglie = _per_parametro(documento)
    assert (soglie["qc.min_frac_phylum"]["misura"], soglie["qc.min_frac_phylum"]["valore"]) == (
        0.7, 0.65)
    assert soglie["qc.max_frac_contaminant"]["misura"] == 0.2


def test_il_file_scritto_cambia_solo_le_soglie_e_porta_origine_e_misura(tmp_path, capsys):
    """
    **Obiettivo**: Verificare che senza ``--accetta-proposte`` le soglie
    proposte siano righe in commento (il file, riletto, ha gli stessi valori
    di quello di partenza meno quelle soglie, e la validazione le elenca come
    mancanti) e con l'opzione siano attive (stessi valori salvo le soglie, e
    configurazione valida); che ogni riga proposta sia preceduta dal marcatore
    ``[PROPOSTO]`` con misura, artefatto, regola e valore di partenza; che una
    soglia tenuta conservi la sua riga; che i commenti del file di partenza
    restino; e che rilanciare il comando sul file scritto non accumuli
    intestazioni ne' commenti.

    **Razionale scientifico e sistemistico**: Lo strumento propone e l'utente
    conferma: senza il gesto nessuna soglia entra nella configurazione, e con
    il gesto ogni valore dichiara da dove viene. Il resto del file e' del suo
    autore e non va toccato.
    """
    dati = _dati(tmp_path / "out")
    codice, documento, nuovo, _ = _soglie(tmp_path, dati, capsys=capsys)
    assert codice == 0 and documento["proposte_accettate"] is False
    testo = nuovo.read_text(encoding="utf-8")
    proposte = {s["parametro"]: s["valore"] for s in documento["soglie"] if s["esito"] == "proposta"}
    riletto = yaml.safe_load(testo)
    atteso = copy.deepcopy(dati)
    for chiave in proposte:
        gruppo, nome = chiave.split(".")
        atteso[gruppo].pop(nome)
    assert riletto == atteso
    with pytest.raises(ErroreConfigurazione) as rifiuto:
        valida(riletto)
    assert all(c in " ".join(rifiuto.value.problemi) for c in proposte)
    assert "  # max_frac_lost_filter: 0.2501" in testo
    assert re.search(r"## \[PROPOSTO\] da soglie: misura 0\.2 \(perdita media.*?"
                     r"03_filtered/manifest_S2\.json\).*?misura piu' il margine \(0\.05\).*?"
                     r"valore nella configurazione di partenza: 1\.0\.", testo, re.S)
    assert "  max_zeroed_samples: 0" in testo and "## [SOGLIE] tenuta" in testo

    codice, documento, accettato, _ = _soglie(tmp_path, dati, "--accetta-proposte", capsys=capsys)
    assert codice == 0 and documento["proposte_accettate"] is True
    riletto = yaml.safe_load(accettato.read_text(encoding="utf-8"))
    atteso = copy.deepcopy(dati)
    for chiave, valore in proposte.items():
        gruppo, nome = chiave.split(".")
        atteso[gruppo][nome] = valore
    assert riletto == atteso and valida(riletto)
    assert "accettato con --accetta-proposte" in accettato.read_text(encoding="utf-8")

    # Un file con i propri commenti, e un secondo passaggio sul file scritto.
    commentato = tmp_path / "commentato.yaml"
    commentato.write_text(
        "# la mia configurazione\n" + yaml.safe_dump(dati, sort_keys=False).replace(
            "  max_frac_lost_filter: 1.0",
            "  # deciso a mano\n  max_frac_lost_filter: 1.0   # limite largo").replace("\n", "\r\n"),
        encoding="utf-8", newline="")
    primo, secondo = tmp_path / "primo.yaml", tmp_path / "secondo.yaml"
    assert main(["soglie", "--config", str(commentato), "--output", str(primo),
                 "--accetta-proposte"]) == 0
    grezzo = primo.read_bytes()
    # Il fine riga del file di partenza e' conservato, e cosi' i suoi commenti:
    # quello accanto al valore sostituito passa nel commento che precede la riga.
    assert grezzo.count(b"\r\n") == grezzo.count(b"\n") > 100
    testo = grezzo.decode("utf-8").replace("\r\n", "\n")
    assert "# la mia configurazione" in testo and "  # deciso a mano" in testo
    assert "## commento sulla riga di partenza: limite largo" in testo
    shutil.rmtree(tmp_path / "out")
    _esecuzione(tmp_path / "out", yaml.safe_load(testo))
    assert main(["soglie", "--config", str(primo), "--output", str(secondo),
                 "--accetta-proposte"]) == 0
    testo = secondo.read_bytes().decode("utf-8").replace("\r\n", "\n")
    assert testo.count("# Soglie proposte da soglie") == 1
    assert len(re.findall(r"^\s*#? ?max_frac_lost_filter:", testo, re.M)) == 1
    # Le soglie ora attive stanno fra misura e suggerimento: tenute. Resta il
    # commento con l'origine del loro valore, e uno solo con l'esito di ora.
    assert testo.count("[PROPOSTO] da soglie") == primo.read_bytes().decode("utf-8").count(
        "[PROPOSTO] da soglie") == 11
    assert len(re.findall(r"## \[SOGLIE\]", testo)) == 12
    capsys.readouterr()


def test_soglie_dichiara_controllo_spento_senza_misura_esclusa_e_omessa(tmp_path, capsys):
    """
    **Obiettivo**: Verificare che il comando dichiari, sulla console, nel
    documento JSON e nel file: «controllo spento» quando il valore proposto
    coincide con il limite (frazione chimerica misurata oltre il margine dal
    limite; tutti i campioni azzerati); «senza misura» per il motivo nullo,
    per il primer con il taglio e per la decontaminazione senza abbastanza
    controlli negativi; «esclusa» per le letture corte quando un file non
    letto per intero ne ha fra le prime; e che una soglia condizionale omessa
    resti omessa. Le soglie senza misura o escluse conservano la loro riga.

    **Razionale scientifico e sistemistico**: Una soglia che resta al valore
    che non ferma senza che nessuno lo dica e' un controllo creduto attivo: il
    comando deve dire di ogni soglia che cosa ne ha fatto, e perche'.
    """
    # Controllo spento: misura al limite.
    codice, documento, nuovo, errore = _soglie(
        tmp_path / "a", _dati(tmp_path / "a/out"), "--accetta-proposte", capsys=capsys,
        chimerica=0.97, azzerati=list(_BIOLOGICI))
    soglie = _per_parametro(documento)
    assert codice == 0
    for chiave in ("qc.warn_frac_chimeric", "qc.stop_frac_chimeric", "qc.max_zeroed_samples"):
        assert soglie[chiave]["controllo_spento"] is True, chiave
        assert chiave in documento["controlli_spenti"]
    assert soglie["qc.stop_frac_chimeric"]["valore"] == 1.0
    assert soglie["qc.max_zeroed_samples"]["valore"] == 3
    testo = nuovo.read_text(encoding="utf-8")
    assert "CONTROLLO SPENTO" in testo and "CONTROLLI SPENTI anche dopo la proposta" in testo
    assert "[CONTROLLO SPENTO]" in errore and "controlli spenti anche dopo la proposta" in errore

    # Senza misura: motivo nullo, primer con il taglio, pochi controlli negativi.
    dati = _dati(tmp_path / "b/out")
    dati["qc"]["conserved_motif"] = None
    dati["filter"]["trimLeft"] = 17
    codice, documento, nuovo, errore = _soglie(tmp_path / "b", dati, "--accetta-proposte",
                                               capsys=capsys, negativi=1)
    soglie = _per_parametro(documento)
    testo = nuovo.read_text(encoding="utf-8")
    for chiave in ("qc.min_motif_frac", "qc.max_primer_hit_frac", "qc.max_frac_contaminant"):
        assert soglie[chiave]["esito"] == "senza_misura" and soglie[chiave]["misura"] is None
        assert soglie[chiave]["valore"] == ESPLORATIVE[chiave] and soglie[chiave]["note"][0]
        assert f"{chiave}: senza misura" in errore
        assert f"  {chiave.split('.')[1]}: {ESPLORATIVE[chiave]}" in testo
    assert testo.count("## [SOGLIE] senza misura in questa esecuzione") == 3
    # Il primer e il motivo non si eseguono: non sono controlli spenti dal valore.
    # La decontaminazione non ha misura qui, ma la sua soglia resta al limite.
    assert documento["controlli_spenti"] == ["qc.max_frac_contaminant"]

    # Esclusa: letture corte fra le prime di un file non letto per intero.
    troncamento = valida(_dati(tmp_path / "x")).filter.truncLen
    codice, documento, nuovo, errore = _soglie(
        tmp_path / "c", _dati(tmp_path / "c/out"), "--accetta-proposte", capsys=capsys,
        minime=[troncamento - 1, 160, 160], esauriti=["no", "si", "si"])
    corte = _per_parametro(documento)["qc.max_frac_short_reads"]
    assert corte["esito"] == "esclusa" and corte["valore"] == 1.0
    assert "non e' stato letto per intero" in corte["note"][0]
    assert "qc.max_frac_short_reads: esclusa" in errore
    assert "esclusa dai suggerimenti" in nuovo.read_text(encoding="utf-8")
    assert "qc.max_frac_short_reads" in documento["controlli_spenti"]

    # Letto per intero: le letture corte delle prime sono quelle del profilo.
    codice, documento, _, _ = _soglie(
        tmp_path / "d", _dati(tmp_path / "d/out"), capsys=capsys,
        minime=[troncamento - 1, 160, 160], esauriti=["si", "si", "si"], corte_per_campione=60,
        corte_s1=(60, 3000))
    corte = _per_parametro(documento)["qc.max_frac_short_reads"]
    # 60 corte su 1500 prime letture dei tre biologici: piu' della frazione su tutte.
    assert (corte["misura"], corte["valore"]) == (0.04, 0.09)

    # Omessa: senza controlli negativi dichiarati.
    senza = _dati(tmp_path / "e/out", {"qc.max_frac_contaminant": None},
                  ctrl={"blank_values": []})
    for nome in ("threshold", "min_blanks", "mode"):
        senza["decontam"].pop(nome)
    codice, documento, nuovo, errore = _soglie(tmp_path / "e", senza, "--accetta-proposte",
                                               capsys=capsys)
    assert codice == 0
    assert _per_parametro(documento)["qc.max_frac_contaminant"]["esito"] == "omessa"
    assert "max_frac_contaminant" not in nuovo.read_text(encoding="utf-8").split("qc:")[1]
    assert "qc.max_frac_contaminant: omessa" in errore


def test_soglie_si_rifiuta_senza_scrivere(tmp_path, capsys):
    """
    **Obiettivo**: Verificare che ``soglie`` esca con 3 senza scrivere nulla, e
    con il motivo sulla console e nel documento JSON, quando: il file di uscita
    esiste gia'; la cartella non contiene un'esecuzione; l'esecuzione e'
    incompleta (manca il manifesto di una fase: il motivo nomina le fasi da
    eseguire); la configurazione non e' quella registrata dall'esecuzione; una
    fase e' calcolata su artefatti a monte che non sono piu' quelli su disco;
    un artefatto con una misura manca. E che la cartella dell'esecuzione non
    sia mai modificata.

    **Razionale scientifico e sistemistico**: Una soglia proposta da misure che
    non sono quelle della configurazione indicata, o da un'esecuzione a meta',
    sembrerebbe fondata senza esserlo. Il comando legge soltanto: non deve
    lasciare tracce nella cartella che legge.
    """
    dati = _dati(tmp_path / "out")
    uscita = tmp_path / "out"

    def impronta() -> dict[str, str]:
        return {str(p.relative_to(uscita)): hashlib.sha256(p.read_bytes()).hexdigest()
                for p in sorted(uscita.rglob("*")) if p.is_file()}

    # Nessuna esecuzione.
    file = tmp_path / "configurazione.yaml"
    file.write_text(yaml.safe_dump(dati), encoding="utf-8")
    capsys.readouterr()
    assert main(["soglie", "--config", str(file), "--output", str(tmp_path / "n.yaml"),
                 "--json"]) == 3
    documento = json.loads(capsys.readouterr().out)
    assert documento["codice_uscita"] == 3 and "non contiene un'esecuzione" in documento["rifiuto"]
    assert documento["file_scritto"] is None and not (tmp_path / "n.yaml").exists()
    assert not uscita.exists()

    _esecuzione(uscita, dati)
    prima = impronta()
    codice, documento, nuovo, _ = _soglie(tmp_path, dati, capsys=capsys)
    assert codice == 0 and impronta() == prima

    # File di uscita esistente.
    capsys.readouterr()
    assert main(["soglie", "--config", str(file), "--output", str(nuovo), "--json"]) == 3
    catturato = capsys.readouterr()
    assert "esiste gia'" in json.loads(catturato.out)["rifiuto"]
    assert "Nulla e' stato scritto" in catturato.err

    def rifiuto(nome: str, configurazione: Path = file) -> str:
        capsys.readouterr()
        destinazione = tmp_path / f"rifiuto_{nome}.yaml"
        assert main(["soglie", "--config", str(configurazione), "--output", str(destinazione),
                     "--json"]) == 3
        assert not destinazione.exists()
        return json.loads(capsys.readouterr().out)["rifiuto"]

    # Configurazione diversa da quella dell'esecuzione.
    altra = tmp_path / "altra.yaml"
    altra.write_text(yaml.safe_dump(_dati(uscita, {"qc.max_frac_lost_filter": 0.9})),
                     encoding="utf-8")
    assert "non e' quella registrata dall'esecuzione" in rifiuto("altra", altra)

    # Stesso digest dell'ultima registrata, ma le fasi sono calcolate con
    # altri valori: una validazione ha registrato la configurazione senza
    # rifare le fasi che dipendono dal parametro cambiato.
    cambiata = _dati(uscita)
    cambiata["filter"]["maxEE"] = 9.0
    registra_risolta(risolvi(valida(copy.deepcopy(cambiata))), uscita)
    diversa = tmp_path / "diversa.yaml"
    diversa.write_text(yaml.safe_dump(cambiata), encoding="utf-8")
    motivo = rifiuto("fasi", diversa)
    assert "non e' stata calcolata con i valori" in motivo and "S2" in motivo
    for registrata in (uscita / "00_config").glob("resolved_*.yaml"):
        registrata.unlink()

    # Un artefatto che si apre ma non ha la forma attesa.
    guasto = uscita / "06_seqtab" / "tabella.json"
    contenuto = guasto.read_bytes()
    guasto.write_text('{"varianti": null}', encoding="utf-8")
    assert "non ha la forma attesa" in rifiuto("forma")
    guasto.write_bytes(contenuto)

    # Artefatto mancante.
    tolto = uscita / "07_chimera" / "chimere.json"
    contenuto = tolto.read_bytes()
    tolto.unlink()
    assert "07_chimera/chimere.json" in rifiuto("artefatto")
    tolto.write_bytes(contenuto)

    # Fase calcolata su artefatti superati.
    manifesto = uscita / "01_input_validation" / nome_manifesto_passo("S0")
    originale = manifesto.read_text(encoding="utf-8")
    manifesto.write_text(originale.replace("sha256:S0", "sha256:rifatta"), encoding="utf-8")
    assert "non sono piu' quelli su disco" in rifiuto("superata")
    manifesto.write_text(originale, encoding="utf-8")

    # Esecuzione incompleta.
    (uscita / "12_final" / nome_manifesto_passo("S14")).unlink()
    motivo = rifiuto("incompleta")
    assert "non contiene un'esecuzione conclusa" in motivo and "S14" in motivo


def test_nel_codice_delle_soglie_non_ci_sono_valori_dei_dataset():
    """
    **Obiettivo**: Verificare sul sorgente del pacchetto ``soglie`` che non vi
    compaia alcuna stringa delle configurazioni dei dataset distribuiti con il
    repository (colonne, etichette, espressioni, taxa), ne' il loro nome; e che
    ogni margine sia una costante con un nome, uguale per ogni dataset.

    **Razionale scientifico e sistemistico**: Un margine o un limite preso dai
    dati di un dataset tornerebbe a essere un predefinito di quel dataset. I
    limiti vengono dallo schema e dai confronti, i margini sono dichiarati.
    """
    if not (RADICE / "dati" / "osd734" / "config_osd734.yaml").is_file():
        pytest.skip("la cartella dati/ non e' disponibile in questo ambiente")
    sorgente = "\n".join(
        p.read_text(encoding="utf-8") for p in Path(tabella.__file__).parent.glob("*.py"))
    assert not re.search(r"osd[-_ ]?\d+|comunita_nota", sorgente, re.I)
    from amplicon16s import interfaccia

    dello_schema = {str(v) for p in interfaccia.schema()["parametri"]
                    for v in (p["valori_ammessi"] or [])}

    def stringhe(valore):
        if isinstance(valore, str):
            yield valore
        elif isinstance(valore, list):
            for elemento in valore:
                yield from stringhe(elemento)

    trovate = []
    for nome in ("osd734", "osd276", "osd170", "comunita_nota"):
        dati = yaml.safe_load(
            (RADICE / "dati" / nome / f"config_{nome}.yaml").read_text(encoding="utf-8"))
        for gruppo, valori in dati.items():
            if gruppo in ("io", "run", "tax"):
                continue
            for testo in (t for v in valori.values() for t in stringhe(v)):
                # "campione" e' anche una parola dell'italiano, usata nei messaggi.
                if len(testo) >= 4 and testo not in dello_schema | {"campione"} and re.search(
                        rf"(?<![A-Za-z0-9_]){re.escape(testo)}(?![A-Za-z0-9_])", sorgente):
                    trovate.append((nome, gruppo, testo))
    assert trovate == []
    assert {s.margine for s in SOGLIE} == {
        tabella.MARGINE_FRAZIONE, tabella.MARGINE_ARRESTO_CHIMERE, tabella.MARGINE_VARIANTI,
        tabella.MARGINE_AZZERATI, tabella.MARGINE_BASI}


# --------------------------------------------------------------------------- #
# Una soglia proposta ferma un peggioramento: i controlli in Python            #
# --------------------------------------------------------------------------- #


def _qc(tmp_path: Path, documento: dict) -> Any:
    """La configurazione con le soglie che ``soglie`` ha proposto o tenuto."""
    valori = {s["parametro"]: s["valore"] for s in documento["soglie"]}
    return valida(_dati(tmp_path / "qc", valori))


def _contesto(config: Any, scansione: dict[str, StatisticheFile]) -> SimpleNamespace:
    """Cio' che i gate sulle letture leggono dal contesto."""
    inventario = [
        Campione(c, c, ClasseCampione.BIOLOGICO if c in _BIOLOGICI
                 else ClasseCampione.CONTROLLO_NEGATIVO, "x") for c in scansione]
    return SimpleNamespace(config=config, inventario=inventario, scansione=scansione)


def _file(nome: str, letture: int, primer: int = 0, motivo: int = 0, corte: int = 0):
    return StatisticheFile(nome, letture, 100, 170, primer, motivo, esaurito=True,
                           piu_corte=corte)


def test_ogni_soglia_proposta_lascia_passare_la_misura_e_ferma_appena_oltre(tmp_path, capsys):
    """
    **Obiettivo**: Verificare, per le soglie il cui controllo e' una funzione
    Python (primer e motivo in G10, letture corte in G09 e S1, campioni
    azzerati e perdita in S2, avviso e arresto delle chimere in S6, phylum in
    S8), che con la soglia proposta da ``soglie`` la misura da cui viene non
    faccia scattare il controllo, e che una misura appena oltre la soglia lo
    faccia scattare con il codice del catalogo dichiarato nella tabella.

    **Razionale scientifico e sistemistico**: Una soglia proposta ha senso
    solo se e' davvero il confine del controllo: larga abbastanza da non
    fermare la ripresa sugli stessi dati, stretta abbastanza da fermare il
    primo peggioramento oltre il margine dichiarato.
    """
    codice, documento, _, _ = _soglie(tmp_path, _dati(tmp_path / "out"), capsys=capsys)
    assert codice == 0
    soglie = _per_parametro(documento)
    config = _qc(tmp_path, documento)
    qc = config.qc

    # G10, primer: per file. La misura era 0,01; la soglia 0,06.
    assert qc.max_primer_hit_frac == 0.06
    dentro = {c: _file(c, 1000, primer=p, motivo=900)
              for c, p in zip(_BIOLOGICI, (2, 10, 60), strict=True)}
    assert _g10_primer_assente(_contesto(config, dentro))[0] == []
    oltre = {**dentro, "LIB003": _file("LIB003", 1000, primer=61, motivo=900)}
    assert [v.codice for v in _g10_primer_assente(_contesto(config, oltre))[0]] == ["E-S0-10"]
    # Un controllo negativo con il primer non conta.
    con_negativo = {**dentro, "LIB004": _file("LIB004", 100, primer=90)}
    assert _g10_primer_assente(_contesto(config, con_negativo))[0] == []

    # G10, motivo: mediana. La misura era 0,8; la soglia 0,75.
    assert qc.min_motif_frac == 0.75
    for motivi, codici in (((900, 750, 700), []), ((900, 749, 700), ["E-S0-16"])):
        scansione = {c: _file(c, 1000, motivo=m)
                     for c, m in zip(_BIOLOGICI, motivi, strict=True)}
        assert [v.codice for v in _g10_primer_assente(_contesto(config, scansione))[0]] == codici

    # G09 e S1, letture corte: per classe. La misura era 0,01; la soglia 0,06.
    assert qc.max_frac_short_reads == 0.06
    for corte, codici in ((180, []), (181, ["E-S0-09"])):
        scansione = {c: _file(c, 1000, corte=k)
                     for c, k in zip(_BIOLOGICI, (corte, 0, 0), strict=True)}
        esito = _g09_troncamento_compatibile(_contesto(config, scansione))[0]
        assert [v.codice for v in esito] == codici
    assert classi_oltre_soglia({"biologico": {"frazione": 0.06}}, qc.max_frac_short_reads) == []
    assert classi_oltre_soglia(
        {"biologico": {"frazione": 0.060001}}, qc.max_frac_short_reads) == ["biologico"]

    # S2, campioni azzerati (misura 0, soglia 0) e perdita media (0,2; 0,2501).
    classi = {c: ClasseCampione.BIOLOGICO for c in _BIOLOGICI}
    classi |= {c: ClasseCampione.CONTROLLO_NEGATIVO for c in _NEGATIVI}
    ingresso = {c: 10000 for c in classi}
    assert qc.max_zeroed_samples == 0 and qc.max_frac_lost_filter == 0.2501
    uscita = {**{c: 8000 for c in _BIOLOGICI}, **{c: 0 for c in _NEGATIVI}}
    assert controlla_filtro(ingresso, uscita, classi, qc)["perdita_media_controllata"] == 0.2
    uscita["LIB001"] = 7496      # perdita media 0,2168: ancora dentro
    controlla_filtro(ingresso, uscita, classi, qc)
    with pytest.raises(ErrorePipeline) as errore:
        controlla_filtro(ingresso, {**uscita, "LIB001": 6496}, classi, qc)   # 0,2501(3)
    assert errore.value.codice == "E-S2-02"
    with pytest.raises(ErrorePipeline) as errore:
        controlla_filtro(ingresso, {**uscita, "LIB002": 0}, classi, qc)
    assert errore.value.codice == "E-S2-01"

    # S6: avviso a 0,17 e arresto a 0,22; la misura era 0,12.
    assert (qc.warn_frac_chimeric, qc.stop_frac_chimeric) == (0.17, 0.22)

    def chimere(frazione: float) -> dict[str, Any]:
        return {"frazione_controllata": frazione, "letture_chimeriche_controllate": 1,
                "letture_controllate": 2}

    assert controlla_chimere(chimere(0.12), qc) is None
    assert controlla_chimere(chimere(0.17), qc) is None
    assert "qc.warn_frac_chimeric" in controlla_chimere(chimere(0.1701), qc)
    assert "qc.warn_frac_chimeric" in controlla_chimere(chimere(0.22), qc)
    with pytest.raises(ErrorePipeline) as errore:
        controlla_chimere(chimere(0.2201), qc)
    assert errore.value.codice == "E-S6-01"

    # S8: phylum, per classe controllata. La misura era 0,9; la soglia 0,85.
    assert qc.min_frac_phylum == 0.85
    controlla_copertura({"biologico": {"varianti": 0.85},
                         "controllo_negativo": {"varianti": 0.1}}, qc)
    with pytest.raises(ErrorePipeline) as errore:
        controlla_copertura({"biologico": {"varianti": 0.8499}}, qc)
    assert errore.value.codice == "E-S8-02"
    # I codici sono quelli che la tabella dichiara per ogni soglia.
    assert soglie["qc.min_frac_phylum"]["codici"] == ["E-S8-02"]
    assert soglie["qc.stop_frac_chimeric"]["codici"] == ["E-S6-01"]


# --------------------------------------------------------------------------- #
# Catene complete sul dataset sintetico                                        #
# --------------------------------------------------------------------------- #

_FINALI = ("ps_final.rds", "conteggi.tsv", "tassonomia.tsv", "metadati.tsv", "sequenze.fasta")


def _sha(percorso: Path) -> str:
    return hashlib.sha256(percorso.read_bytes()).hexdigest()


def _scrivi_config(percorso: Path, dati: dict) -> Path:
    percorso.write_text(yaml.safe_dump(dati, sort_keys=False), encoding="utf-8")
    return percorso


def _codice_di_arresto(uscita: Path) -> str:
    return json.loads((uscita / "99_logs" / "punto_di_ripresa.json").read_text("utf-8"))["codice"]


@pytest.fixture(scope="module")
def esplorativa(bioc, tmp_path_factory):
    """Un'esecuzione esplorativa conclusa sul dataset sintetico, con le soglie
    al valore che non ferma, e il file che ``soglie`` ne ricava.
    """
    base = tmp_path_factory.mktemp("soglie")
    dataset = sintetico.genera(base / "dati", biologici=8, negativi=3)
    # Un taxon del riferimento sintetico escluso: i filtri finali tolgono
    # letture, e la frazione trattenuta non e' al suo massimo.
    dati = sintetico.configurazione(dataset, base / "out", filt={"exclude_taxa": ["Classe_0"]})
    for chiave, valore in ESPLORATIVE.items():
        if not chiave.startswith("filter."):
            gruppo, nome = chiave.split(".")
            dati[gruppo][nome] = valore
    partenza = _scrivi_config(base / "esplorativa.yaml", dati)
    chiudi()
    assert main(["run", "--config", str(partenza)]) == 0
    chiudi()
    finali = {f: _sha(base / "out" / "12_final" / f) for f in _FINALI}
    proposta = base / "con_soglie.yaml"
    assert main(["soglie", "--config", str(partenza), "--output", str(proposta),
                 "--accetta-proposte"]) == 0
    return SimpleNamespace(base=base, dataset=dataset, dati=dati, partenza=partenza,
                           proposta=proposta, finali=finali)


def test_la_ripresa_con_le_soglie_proposte_non_cambia_i_risultati(esplorativa, capsys):
    """
    **Obiettivo**: Verificare sulla catena completa che l'esecuzione
    esplorativa dichiari i controlli spenti; che ``resume`` con il file
    scritto da ``soglie --accetta-proposte`` si concluda senza arresti ne'
    avvisi dovuti alle soglie appena scritte; che oggetto finale, conteggi,
    tassonomia, metadati e sequenze siano identici byte per byte a quelli
    dell'esecuzione esplorativa; e che dopo la ripresa nessun controllo con una
    misura resti spento.

    **Razionale scientifico e sistemistico**: Le soglie fermano o avvisano, non
    calcolano: attivarle non deve cambiare alcun risultato, e una soglia
    appena proposta non deve fermare i dati da cui e' stata ricavata.
    """
    base = esplorativa.base
    capsys.readouterr()
    assert main(["validate", "--config", str(esplorativa.partenza), "--json"]) == 0
    spenti = {s["parametro"] for s in json.loads(capsys.readouterr().out)["controlli_spenti"]}
    assert "qc.max_frac_lost_filter" in spenti and "qc.min_frac_reads_retained" in spenti
    chiudi()

    assert main(["resume", "--config", str(esplorativa.proposta)]) == 0
    testo = capsys.readouterr().out
    chiudi()
    assert "Esecuzione completata" in testo
    assert {f: _sha(base / "out" / "12_final" / f) for f in _FINALI} == esplorativa.finali
    for passo in ("S1", "S6"):
        nodo = GRAFO.nodo(next(p for p in GRAFO.ordine() if str(p) == passo))
        manifesto = json.loads(
            (base / "out" / nodo.cartella.value / nome_manifesto_passo(passo)).read_text("utf-8"))
        assert not [d for d in manifesto["degradazioni"]
                    if d["codice"] in ("E-S1-01", "E-S6-02")], manifesto["degradazioni"]
    assert main(["validate", "--config", str(esplorativa.proposta), "--json"]) == 0
    documento = json.loads(capsys.readouterr().out)
    # conserved_motif e' nullo e trimLeft e' zero nel dataset sintetico.
    assert documento["controlli_spenti"] == []
    chiudi()


@pytest.mark.parametrize(("chiave", "passo", "codice"), [
    ("qc.min_frac_reads_retained", "S14", "E-S14-01"),
    ("qc.max_frac_contaminant", "S12", "E-S12-02"),
    ("qc.max_asv_count", "S5", "E-S5-01"),
])
def test_una_soglia_appena_piu_stretta_della_misura_ferma_la_fase(
    esplorativa, capsys, chiave, passo, codice
):
    """
    **Obiettivo**: Verificare, per le soglie il cui controllo sta dentro il
    calcolo di una fase (letture trattenute in S14, contaminanti in S12,
    varianti in S5), che con la soglia proposta la catena si concluda e che
    con la stessa soglia spostata di un passo oltre la misura la fase si fermi
    con il codice dichiarato nella tabella.

    **Razionale scientifico e sistemistico**: Spostare la soglia di un passo
    oltre la misura equivale a peggiorare la misura di quel passo: se la fase
    non si ferma, il parametro non governa il controllo che dichiara.
    """
    base = esplorativa.base
    proposte = yaml.safe_load(esplorativa.proposta.read_text(encoding="utf-8"))
    documento_misure = misure.leggi(base / "out", tabella.valori_piatti(valida(
        copy.deepcopy(proposte))))
    misura = documento_misure[chiave]
    if misura.valore is None:
        pytest.skip(f"{chiave}: senza misura in questa esecuzione ({misura.assente})")
    voce = tabella.soglia(chiave)
    passo_stretto = Fraction(1, 10000) if voce.unita is Unita.FRAZIONE else 1
    stretta = misura.valore + (-passo_stretto if voce.direzione is Direzione.MASSIMA
                               else passo_stretto)
    if stretta < 0 or (voce.unita is Unita.FRAZIONE and stretta > 1):
        pytest.skip(f"{chiave}: la misura e' all'estremo dell'intervallo, non c'e' una "
                    "soglia piu' stretta")
    dati = copy.deepcopy(proposte)
    dati["io"]["out_root"] = str(base / f"out_{passo}")
    gruppo, nome = chiave.split(".")
    dati[gruppo][nome] = float(stretta) if voce.unita is Unita.FRAZIONE else int(stretta)
    # Stessi artefatti a monte: la fase da fermare si esegue su una copia.
    shutil.copytree(base / "out", base / f"out_{passo}")
    chiudi()
    assert main(["resume", "--config", str(_scrivi_config(base / f"{passo}.yaml", dati))]) == 4
    chiudi()
    capsys.readouterr()
    assert _codice_di_arresto(base / f"out_{passo}") == codice == voce.codici[0]


def _peggiora_la_qualita(letture: Path, quota: float) -> None:
    """Porta al minimo la qualita' della coda di una quota delle letture di
    ogni file: il filtro le tronca sotto ``filter.truncLen`` e le scarta.
    """
    for percorso in sorted(letture.glob("*.fq.gz")):
        with gzip.open(percorso, "rt", encoding="ascii") as file:
            righe = file.read().splitlines()
        for i in range(3, len(righe), 4):
            if (i // 4) % 100 < quota * 100:
                righe[i] = righe[i][:40] + "#" * (len(righe[i]) - 40)
        with open(percorso, "wb") as file:
            with gzip.GzipFile(fileobj=file, mode="wb", mtime=0, filename="") as archivio:
                archivio.write(("\n".join(righe) + "\n").encode("ascii"))


def test_letture_di_bassa_qualita_fermano_l_esecuzione_con_le_soglie_proposte(
    esplorativa, capsys
):
    """
    **Obiettivo**: Verificare con una catena completa che, peggiorando i dati
    in modo controllato (la coda di meta' delle letture di ogni campione
    portata alla qualita' minima), l'esecuzione con le soglie proposte
    dall'esecuzione buona si fermi nel filtro con ``E-S2-02``, mentre con le
    soglie esplorative gli stessi dati arrivano oltre il filtro.

    **Razionale scientifico e sistemistico**: E' il caso per cui la soglia
    esiste: una corsa di sequenziamento peggiore della prima. Con il controllo
    spento passerebbe inosservata, con la soglia proposta no.
    """
    base = esplorativa.base
    peggiori = base / "dati_peggiori"
    shutil.copytree(esplorativa.dataset.cartella, peggiori)
    _peggiora_la_qualita(peggiori / "letture", 0.5)

    def configurazione(file: Path, nome: str) -> Path:
        dati = yaml.safe_load(file.read_text(encoding="utf-8"))
        dati["io"]["fastq_dir"] = str(peggiori / "letture")
        dati["io"]["assay_table"] = str(peggiori / "campioni.tsv")
        dati["io"]["out_root"] = str(base / nome)
        return _scrivi_config(base / f"{nome}.yaml", dati)

    chiudi()
    assert main(["run", "--config", str(configurazione(esplorativa.proposta, "peggiori"))]) == 4
    chiudi()
    assert _codice_di_arresto(base / "peggiori") == "E-S2-02"
    # Con il controllo spento il filtro non ferma: l'esecuzione va oltre S2.
    from amplicon16s.runner.executor import Esecutore
    from amplicon16s.runner.graph import Passo
    from amplicon16s.runner.project import ProjectRun
    from amplicon16s.config.schema import carica

    run = ProjectRun(carica(configurazione(esplorativa.partenza, "peggiori_spenta")))
    esito = Esecutore(run, fino_a=Passo.S2).esegui()
    chiudi()
    assert [str(r.passo) for r in esito.eseguite][-1] == "S2" and esito.punto is None
    capsys.readouterr()


def test_un_riferimento_che_lascia_varianti_senza_phylum_ferma_l_esecuzione(
    esplorativa, capsys
):
    """
    **Obiettivo**: Verificare con una catena completa che, sostituendo il
    riferimento tassonomico con uno in cui meta' delle sequenze porta il solo
    regno, l'esecuzione con le soglie proposte dall'esecuzione buona si fermi
    nell'assegnazione tassonomica con ``E-S8-02``.

    **Razionale scientifico e sistemistico**: Un riferimento sbagliato, o di
    un'altra versione, lascia le varianti senza phylum senza alcun errore: la
    frazione minima proposta dall'esecuzione buona e' cio' che lo rivela.
    """
    base = esplorativa.base
    povero = base / "riferimento_povero.fa.gz"
    with gzip.open(sintetico.RIFERIMENTO, "rt", encoding="ascii") as file:
        righe = file.read().splitlines()
    intestazioni = 0
    for i, riga in enumerate(righe):
        if riga.startswith(">"):
            intestazioni += 1
            if intestazioni % 2:
                righe[i] = riga.split(";")[0] + ";"
    with open(povero, "wb") as file:
        with gzip.GzipFile(fileobj=file, mode="wb", mtime=0, filename="") as archivio:
            archivio.write(("\n".join(righe) + "\n").encode("ascii"))
    dati = yaml.safe_load(esplorativa.proposta.read_text(encoding="utf-8"))
    dati["tax"]["ref_fasta"] = str(povero)
    dati["tax"]["ref_md5"] = hashlib.md5(povero.read_bytes()).hexdigest()
    dati["io"]["out_root"] = str(base / "out_S8")
    shutil.copytree(base / "out", base / "out_S8")
    chiudi()
    assert main(["resume", "--config", str(_scrivi_config(base / "S8.yaml", dati))]) == 4
    chiudi()
    capsys.readouterr()
    assert _codice_di_arresto(base / "out_S8") == "E-S8-02"
