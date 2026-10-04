r"""Suite di test della settimana 21: filtri finali e serializzazione, fasi S13 e S14.

1. Inquadramento nel Piano Operativo
------------------------------------
Settimana 21 (W21), Fase F5 (filtri per profondita', tassonomici e di
prevalenza S13; serializzazione, export e validazione S14; ``12_final/``).

2. Moduli sorgente coperti
--------------------------
* ``src/amplicon16s/steps/s13_filtri.py``, ``R/13_filtri.R``
* ``src/amplicon16s/steps/s14_finale.py``, ``R/14_finale.R``, ``R/lib/export.R``
* ``src/amplicon16s/config/schema.py`` (``qc.min_frac_reads_retained``),
  ``src/amplicon16s/errors/catalog.py`` (``E-S13-02``, ``E-S13-03``, ``E-S14-01``)

3. Cosa valuta questo file
--------------------------
- il filtro per profondita' confronta ogni campione con la soglia della sua
  piastra allo stadio dichiarato da S11, sulle letture del tracciamento;
- S13 e S14 sul sottoinsieme di prova: l'oggetto finale con i soli biologici,
  i controlli a parte, le esclusioni e le varianti rimosse con il motivo,
  gli identificativi non rinumerati, il denominatore della prevalenza;
- gli export ricostruiscono l'oggetto, letti da codice indipendente da quello
  di S14; il manifesto dei checksum e' completo; gli artefatti sono identici
  fra due esecuzioni e con impostazioni locali diverse;
- un campione svuotato dal filtro tassonomico esce con ``E-S13-03``; uno
  svuotato dal filtro di prevalenza ferma con ``E-S13-02``; i campioni con
  poche letture finali escono senza fermare; una frazione trattenuta sotto
  ``qc.min_frac_reads_retained`` ferma con ``E-S14-01``;
- i test sui dati reali girano tutti sullo stesso processo di pytest-xdist;
- sul dataset completo: la catena arriva alla fine.

4. Comandi Bash e scenari di esecuzione
---------------------------------------
    ``<immagine>`` e' l'immagine del container della pipeline; quella corrente
    e' indicata in ``test.txt``, sezione 1.3.

    1. Modalità locale standard (R di base con jsonlite, senza Bioconductor né
       dati reali):
       pytest tests/test_w21_finale.py -v

    2. Modalità container Docker standard (sottoinsieme ridotto con
       R/Bioconductor):
       docker run --rm \
         -e PYTHONPATH=/app/src \
         -e AMPLICON16S_R_DIR=/app/R \
         -v "$(pwd)":/app \
         -w /app \
         <immagine> \
         pytest -o cache_dir=/tmp/.pytest_cache tests/test_w21_finale.py -v

    3. Modalità container Docker completa (con i dati reali OSD-734):
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
         pytest -o cache_dir=/tmp/.pytest_cache tests/test_w21_finale.py -v

5. Risultato atteso
-------------------
Vedi ``test.txt``, scheda W21.

6. Razionale scientifico e sistemistico
---------------------------------------
- La soglia di profondita' vale sulla grandezza con cui e' stata stimata: le
  letture dell'oggetto dopo decontaminazione e filtri sono un'altra cosa.
- L'oggetto finale e' destinato alle analisi ecologiche: un controllo al suo
  interno sarebbe trattato come un campione ambientale.
- Gli export servono a chi non usa lo stesso ambiente di calcolo: devono
  bastare a ricostruire l'oggetto, e non dipendere dalla macchina.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
import os
import re
import subprocess
import sys
from pathlib import Path

import pytest
from conftest import copia_esecuzione
from sottoinsieme import config_ridotta, motivo_pacchetti_r_assenti

from amplicon16s.config.resolve import risolvi
from amplicon16s.errors.exceptions import ErrorePipeline
from amplicon16s.io_layer.artifacts import Fase
from amplicon16s.logging.logger import chiudi
from amplicon16s.metadata.models import ClasseCampione
from amplicon16s.rbridge.payload import PREFISSO
from amplicon16s.rbridge.runner import trova_rscript
from amplicon16s.runner.executor import Conclusione, Esecutore
from amplicon16s.runner.graph import Passo
from amplicon16s.runner.project import StatoPasso, passi_realizzati
from amplicon16s.io_layer.conteggi import leggi_conteggi
from amplicon16s.steps.s02_filter import NOME_PREFILTRO
from amplicon16s.steps.s10_phyloseq import NOME_OGGETTO
from amplicon16s.steps.s13_filtri import NOME_LUNGHEZZA, esclusi_per_profondita


@pytest.fixture(autouse=True)
def uscite_pulite():
    """Chiude le uscite del log prima e dopo ogni test, perché nessun handler resti
    aperto sulla cartella temporanea.
    """
    chiudi()
    yield
    chiudi()


_MOTIVO_BIOC = motivo_pacchetti_r_assenti(
    "dada2", "ggplot2", "ShortRead", "jsonlite", "phyloseq", "Biostrings", "decontam"
)


@pytest.fixture
def bioc():
    """Richiede R con phyloseq, decontam e i pacchetti delle fasi a monte: salta
    senza, ma in CI fallisce.
    """
    if _MOTIVO_BIOC is not None:
        if os.environ.get("AMPLICON16S_RICHIEDI_BIOC") == "1":
            pytest.fail(f"Bioconductor e' richiesto in questo ambiente: {_MOTIVO_BIOC}")
        pytest.skip(_MOTIVO_BIOC)


def _tsv(percorso: Path) -> list[dict[str, str]]:
    """Le righe di una tabella separata da tabulazioni."""
    with open(percorso, encoding="utf-8", newline="") as file:
        return list(csv.DictReader(file, delimiter="\t"))


def _r(codice: str, cartella: Path) -> object:
    """Esegue codice R che scrive un oggetto JSON in ``uscita``."""
    uscita = cartella / "uscita.json"
    script = cartella / "codice.R"
    script.write_text(f"uscita <- {json.dumps(str(uscita))}\n{codice}\n", encoding="utf-8")
    subprocess.run([str(trova_rscript()), "--vanilla", str(script)], check=True,
                   capture_output=True)
    return json.loads(uscita.read_text(encoding="utf-8"))


def _impronte(cartella: Path) -> dict[str, str]:
    """L'MD5 di ogni file della cartella, esclusi i file del ponte e i manifesti."""
    return {
        p.name: hashlib.md5(p.read_bytes()).hexdigest()
        for p in sorted(cartella.iterdir())
        if not p.name.startswith((PREFISSO, "manifest"))
    }


def _riepilogo(run) -> dict:
    """``filtri_riepilogo.json`` di S13."""
    return json.loads((run.albero.cartella(Fase.FINAL) / "filtri_riepilogo.json").read_text())


# --------------------------------------------------------------------------- #
# 1. Il filtro per profondita' e i parametri                                   #
# --------------------------------------------------------------------------- #


SOGLIA = {
    "per_piastra": {
        "1": {"valore": 1000, "stadio": "raw", "origine": "ripiego: qc.min_reads_raw"},
        "3": {"valore": 9000, "stadio": "nonchimeric", "origine": "curva della piastra"},
    },
    "senza_piastra": {"valore": 1000, "stadio": "raw", "origine": "ripiego: qc.min_reads_raw"},
}


def test_il_filtro_per_profondita_usa_lo_stadio_della_piastra():
    """
    **Obiettivo**: Verificare che ogni campione si confronti con la soglia della
    sua piastra sulle letture dello stadio dichiarato: nella piastra 3 le
    letture senza chimere, nella piastra 1 e senza piastra le grezze; e che il
    motivo dell'esclusione dica letture, stadio, soglia e piastra.

    **Razionale scientifico e sistemistico**: Le letture grezze e quelle senza
    chimere non sono confrontabili: un campione con 9.500 letture grezze e
    8.000 senza chimere e' sotto la soglia di 9.000 che vale sulle seconde.
    """
    letture = {
        "raw": {"A": 1500, "B": 900, "C": 9500, "D": 9500, "E": 800},
        "nonchimeric": {"A": 700, "B": 600, "C": 8000, "D": 9200, "E": 700},
    }
    campioni = [("A", "1"), ("B", "1"), ("C", "3"), ("D", "3"), ("E", None)]
    tenuti, esclusi = esclusi_per_profondita(campioni, SOGLIA, letture)
    assert tenuti == ["A", "D"]
    assert [e["accession"] for e in esclusi] == ["B", "C", "E"]
    assert esclusi[1]["motivo"] == (
        "8000 letture allo stadio nonchimeric, meno della soglia 9000 della piastra 3 "
        "(curva della piastra)"
    )
    assert "piastra non nota" in esclusi[2]["motivo"]
    with pytest.raises(RuntimeError, match="stadio"):
        esclusi_per_profondita([("A", "1")], {**SOGLIA, "per_piastra": {
            "1": {"valore": 1, "stadio": "decontaminate", "origine": "x"}}}, letture)


def test_i_parametri_di_s13_e_s14(tmp_path):
    """
    **Obiettivo**: Verificare i parametri dei filtri e della validazione
    (``qc.min_frac_reads_retained`` 0,40, nuovo), che i filtri cambino
    l'impronta di S13 e non quella delle fasi a monte, e che gli export
    cambino l'impronta di S14 e non quella di S13.

    **Razionale scientifico e sistemistico**: Cambiare la soglia di prevalenza
    deve rifare S13 e S14, non la catena.
    """
    config = config_ridotta(tmp_path)
    assert (config.prev.min_fraction, config.prev.min_count, config.qc.min_reads_final,
            config.qc.min_frac_reads_retained) == (0.01, 2, 1000, 0.40)
    prima = risolvi(config)
    passi = passi_realizzati()
    for variazione, cambiata in (
        ({"prev": {"min_fraction": 0.05}}, Passo.S13),
        ({"filt": {"exclude_taxa": ["Chloroplast"]}}, Passo.S13),
        ({"out": {"export_flat": False}}, Passo.S14),
        ({"qc": {"min_frac_reads_retained": 0.5}}, Passo.S14),
    ):
        dopo = risolvi(config_ridotta(tmp_path, **variazione))
        for passo, fase in passi.items():
            diversa = fase.calcolata_su(prima, {}) != fase.calcolata_su(dopo, {})
            assert diversa is (passo is cambiata), (variazione, passo)


# --------------------------------------------------------------------------- #
# 2. S13 e S14 sul sottoinsieme di prova                                       #
# --------------------------------------------------------------------------- #


def test_l_oggetto_finale_ha_solo_i_biologici_e_i_controlli_stanno_a_parte(
    bioc, finale_calcolata, tmp_path
):
    """
    **Obiettivo**: Verificare che S13 e S14 si concludano dopo S12, che
    ``ps_final.rds`` contenga solo campioni biologici, quelli dell'inventario
    meno gli esclusi di ``esclusioni.tsv``, che ``ps_controlli.rds`` contenga
    tutti i controlli positivi e negativi, e che gli identificativi delle
    varianti finali siano un sottoinsieme ordinato di quelli di S10.

    **Razionale scientifico e sistemistico**: Un controllo nell'oggetto finale
    sarebbe trattato come un campione ambientale; un campione escluso senza
    traccia sparirebbe in silenzio.
    """
    run, esito = finale_calcolata
    assert esito.conclusione is Conclusione.COMPLETATA
    assert [r.passo for r in esito.eseguite] == [Passo.S13, Passo.S14]
    inventario = run.valuta().inventario
    cartella = run.albero.cartella(Fase.FINAL)
    letti = _r(f"""
suppressPackageStartupMessages(library(phyloseq))
f <- readRDS({json.dumps(str(cartella / 'ps_final.rds'))})
c <- readRDS({json.dumps(str(cartella / 'ps_controlli.rds'))})
i <- readRDS({json.dumps(str(run.albero.cartella(Fase.PHYLOSEQ) / NOME_OGGETTO))})
jsonlite::write_json(list(
  finale = sample_names(f), classi = unique(as(sample_data(f), "data.frame")$classe),
  controlli = sample_names(c), varianti = taxa_names(f), iniziali = taxa_names(i)
), uscita)
""", tmp_path)
    esclusi = [e["accession"] for e in _tsv(cartella / "esclusioni.tsv")]
    biologici = [c.accession for c in inventario.di_classe(ClasseCampione.BIOLOGICO)]
    assert letti["classi"] == ["biologico"]
    assert letti["finale"] == [a for a in biologici if a not in set(esclusi)]
    assert set(esclusi) <= set(biologici)
    assert set(letti["controlli"]) == {c.accession for c in inventario if c.classe.e_controllo}
    posizioni = [letti["iniziali"].index(v) for v in letti["varianti"]]
    assert posizioni == sorted(posizioni)
    for e in _tsv(cartella / "esclusioni.tsv"):
        assert e["filtro"] and e["motivo"]
    print(f"\nS13: {dict(esito.eseguite[0].metriche)}\nS14: {dict(esito.eseguite[1].metriche)}")


def test_il_filtro_per_profondita_usa_le_letture_del_tracciamento(bioc, finale_calcolata):
    """
    **Obiettivo**: Verificare che i campioni esclusi per profondita' siano
    esattamente quelli che il tracciamento, letto allo stadio della soglia di
    ciascuna piastra in ``soglia.json``, mette sotto soglia; e che sulla
    versione ridotta ci siano piastre con entrambi gli stadi.

    **Razionale scientifico e sistemistico**: La soglia si applica alla
    grandezza con cui S11 l'ha stimata, non alla profondita' dell'oggetto dopo
    decontaminazione e filtri.
    """
    run, _ = finale_calcolata
    soglia = json.loads((run.albero.cartella(Fase.CONTROLS) / "soglia.json").read_text())
    assert {v["stadio"] for v in soglia["per_piastra"].values()} == {"raw", "nonchimeric"}
    letture = {
        "raw": leggi_conteggi(run.albero.cartella(Fase.FILTERED) / NOME_PREFILTRO),
        "nonchimeric": leggi_conteggi(run.albero.cartella(Fase.CHIMERA) / NOME_LUNGHEZZA),
    }
    attesi = set()
    for c in run.valuta().inventario.di_classe(ClasseCampione.BIOLOGICO):
        voce = soglia["per_piastra"][c.piastra]
        if letture[voce["stadio"]][c.accession] < voce["valore"]:
            attesi.add(c.accession)
    esclusi = {e["accession"] for e in _tsv(run.albero.cartella(Fase.FINAL) / "esclusioni.tsv")
               if e["filtro"] == "profondita"}
    assert esclusi == attesi


def test_la_prevalenza_si_calcola_sui_campioni_tenuti(bioc, finale_calcolata, tmp_path):
    """
    **Obiettivo**: Verificare che il denominatore della prevalenza siano i
    biologici tenuti dopo i filtri per profondita' e tassonomico, che
    ``min_campioni`` sia ``ceil(prev.min_fraction * denominatore)``, che ogni
    variante finale sia presente con almeno ``prev.min_count`` letture in
    almeno ``min_campioni`` campioni finali, e che le rimosse lo dichiarino.

    **Razionale scientifico e sistemistico**: Numeratore e denominatore della
    prevalenza si riferiscono allo stesso insieme di campioni, quelli
    dell'oggetto finale.
    """
    run, _ = finale_calcolata
    riepilogo = _riepilogo(run)
    prevalenza = riepilogo["prevalenza"]
    esclusi = riepilogo["campioni"]["esclusi"]
    assert prevalenza["denominatore"] == (
        riepilogo["campioni"]["biologici"] - esclusi["profondita"] - esclusi["tassonomico"])
    assert prevalenza["min_campioni"] == math.ceil(0.01 * prevalenza["denominatore"])
    cartella = run.albero.cartella(Fase.FINAL)
    presenze = _r(f"""
suppressPackageStartupMessages(library(phyloseq))
f <- readRDS({json.dumps(str(cartella / 'ps_final.rds'))})
jsonlite::write_json(as.list(rowSums(as(otu_table(f), "matrix") >= 2)), uscita, auto_unbox = TRUE)
""", tmp_path)
    assert min(presenze.values()) >= prevalenza["min_campioni"]
    rimosse = [r for r in _tsv(cartella / "varianti_rimosse.tsv") if r["filtro"] == "prevalenza"]
    assert rimosse and all("campioni su" in r["motivo"] or "nessuna lettura" in r["motivo"]
                           for r in rimosse)


#: Ricostruisce l'oggetto dagli export con codice scritto qui, indipendente da
#: R/lib/export.R, e lo confronta componente per componente con ps_final.rds.
RICOSTRUISCI = r"""
suppressPackageStartupMessages(library(phyloseq))
leggi <- function(nome) {
  righe <- readLines(file.path(cartella, nome), encoding = "UTF-8")
  campi <- strsplit(righe, "\t", fixed = TRUE)
  intestazione <- campi[[1]]
  corpo <- do.call(rbind, lapply(campi[-1], function(c) c(c, rep("", length(intestazione) - length(c)))))
  colnames(corpo) <- intestazione
  corpo
}
cont <- leggi("conteggi.tsv")
m <- matrix(as.integer(cont[, -1]), nrow = nrow(cont), dimnames = list(cont[, 1], colnames(cont)[-1]))
tx <- leggi("tassonomia.tsv")
t <- tx[, -1, drop = FALSE]; t[t == ""] <- NA; rownames(t) <- tx[, 1]
md <- leggi("metadati.tsv")
d <- as.data.frame(md, stringsAsFactors = FALSE); d[d == ""] <- NA; rownames(d) <- d$accession
fa <- readLines(file.path(cartella, "sequenze.fasta"))
s <- setNames(fa[seq(2, length(fa), 2)], sub("^>", "", fa[seq(1, length(fa), 2)]))
f <- readRDS(file.path(cartella, "ps_final.rds"))
fd <- as(sample_data(f), "data.frame")
jsonlite::write_json(list(
  conteggi = identical(m, as(otu_table(f), "matrix")),
  tassonomia = identical(t, as(tax_table(f), "matrix")),
  metadati = identical(as.matrix(d), as.matrix(fd)) && identical(colnames(d), colnames(fd)),
  sequenze = identical(s, as.character(refseq(f))),
  dimensioni = dim(m)
), uscita, auto_unbox = TRUE)
"""


def test_gli_export_ricostruiscono_l_oggetto_finale(bioc, finale_calcolata, tmp_path):
    """
    **Obiettivo**: Verificare che i quattro export piatti, letti come testo da
    codice indipendente da quello di S14, ricostruiscano conteggi, tassonomia,
    metadati e sequenze di ``ps_final.rds`` identici.

    **Razionale scientifico e sistemistico**: Gli export servono a chi non usa
    lo stesso ambiente di calcolo: devono bastare, da soli, a riavere l'oggetto.
    """
    run, _ = finale_calcolata
    esito = _r(f"cartella <- {json.dumps(str(run.albero.cartella(Fase.FINAL)))}\n" + RICOSTRUISCI,
               tmp_path)
    assert {k: esito[k] for k in ("conteggi", "tassonomia", "metadati", "sequenze")} == {
        "conteggi": True, "tassonomia": True, "metadati": True, "sequenze": True}
    assert esito["dimensioni"][1] == _riepilogo(run)["campioni"]["finali"]


def test_il_manifesto_dei_checksum_e_completo(bioc, finale_calcolata):
    """
    **Obiettivo**: Verificare che ogni file di ``12_final``, salvo i file del
    ponte e i manifesti, sia registrato nel manifesto di S13 o di S14 con il
    suo checksum, e che ``checksum.sha256`` dia l'impronta di ogni file
    consegnato da S14.

    **Razionale scientifico e sistemistico**: Un file consegnato senza
    checksum non si puo' verificare dopo il trasferimento.
    """
    run, _ = finale_calcolata
    cartella = run.albero.cartella(Fase.FINAL)
    registrati = {}
    for passo in (Passo.S13, Passo.S14):
        for voce in run.albero.manifesto_passo(passo, Fase.FINAL).artefatti:
            registrati[voce["nome"]] = voce["checksum"]
    presenti = {p.name for p in cartella.iterdir() if not p.name.startswith((PREFISSO, "manifest"))}
    assert presenti == set(registrati)
    for nome, checksum in registrati.items():
        assert checksum == "sha256:" + hashlib.sha256((cartella / nome).read_bytes()).hexdigest()
    somme = dict(reversed(r.split("  ")) for r in (cartella / "checksum.sha256").read_text().splitlines())
    assert set(somme) == {"ps_final.rds", "conteggi.tsv", "tassonomia.tsv", "metadati.tsv",
                          "sequenze.fasta"}
    for nome, impronta in somme.items():
        assert impronta == hashlib.sha256((cartella / nome).read_bytes()).hexdigest()


@pytest.mark.parametrize("lingua", ["C", "en_US.UTF-8"])
def test_stessi_byte_in_due_esecuzioni_e_con_impostazioni_locali_diverse(
    bioc, finale_calcolata, tmp_path, monkeypatch, lingua
):
    """
    **Obiettivo**: Verificare che S13 e S14, rieseguite con ``LC_ALL`` diverso
    da quello della prima esecuzione, diano in ``12_final`` gli stessi byte,
    oggetto ed export compresi.

    **Razionale scientifico e sistemistico**: Separatore decimale, codifica e
    ordinamento dipendono dalle impostazioni locali: gli export non devono.
    """
    run, _ = finale_calcolata
    attese = _impronte(run.albero.cartella(Fase.FINAL))
    copia = copia_esecuzione(finale_calcolata, tmp_path)
    for passo in (Passo.S14, Passo.S13):
        copia.albero.rimuovi_manifesto_passo(passo, Fase.FINAL)
    monkeypatch.setenv("LC_ALL", lingua)
    esito = Esecutore(copia, fino_a=Passo.S14).esegui()
    assert [r.passo for r in esito.eseguite] == [Passo.S13, Passo.S14]
    assert _impronte(copia.albero.cartella(Fase.FINAL)) == attese


# --------------------------------------------------------------------------- #
# 3. Campioni svuotati, campioni poveri, validazione                           #
# --------------------------------------------------------------------------- #


def _s13_su_decontaminato_modificato(base, cartella: Path, modifica_r: str, **sovrascrivi):
    """Una copia della catena fino a S12 con l'oggetto decontaminato modificato da
    ``modifica_r`` (che lavora su ``ps``), ed S13 calcolata direttamente.
    """
    run = copia_esecuzione(base, cartella / "copia", **sovrascrivi)
    oggetto = run.albero.cartella(Fase.CONTROLS) / "ps_decontaminato.rds"
    _r(f"""
suppressPackageStartupMessages(library(phyloseq))
ps <- readRDS({json.dumps(str(oggetto))})
{modifica_r}
saveRDS(ps, {json.dumps(str(oggetto))})
jsonlite::write_json(list(), uscita)
""", cartella)
    fase = run.fase(Passo.S13)
    contesto = run.contesto(Passo.S13, run.valuta()).ristretto(fase.parametri)
    return run, fase, contesto


def _campione_e_varianti(run, tmp_path) -> dict:
    """Un campione biologico dell'oggetto finale, la sua variante piu' abbondante
    con il genere, e una variante senza letture nei campioni finali.
    """
    finale = run.albero.cartella(Fase.FINAL) / "ps_final.rds"
    decontaminato = run.albero.cartella(Fase.CONTROLS) / "ps_decontaminato.rds"
    return _r(f"""
suppressPackageStartupMessages(library(phyloseq))
f <- readRDS({json.dumps(str(finale))}); d <- readRDS({json.dumps(str(decontaminato))})
a <- sample_names(f)[1]
m <- as(otu_table(d), "matrix")
v <- rownames(m)[which.max(m[, a])]
assente <- rownames(m)[rowSums(m[, sample_names(f), drop = FALSE]) == 0][1]
jsonlite::write_json(list(campione = a, variante = v, genere = as(tax_table(d), "matrix")[v, "Genus"],
                          assente = assente), uscita, auto_unbox = TRUE)
""", tmp_path)


def test_un_campione_svuotato_dal_filtro_tassonomico_esce_con_e_s13_03(
    bioc, finale_calcolata, decontam_calcolata, tmp_path
):
    """
    **Obiettivo**: Verificare che un campione che, dopo il filtro tassonomico,
    resta senza letture (costruito lasciandogli una sola variante e mettendone
    il genere in ``filt.exclude_taxa``) esca dall'oggetto finale con il filtro
    ``tassonomico`` in ``esclusioni.tsv``, e che la fase registri ``E-S13-03``
    e prosegua.

    **Razionale scientifico e sistemistico**: Un campione fatto solo di
    organelli o di varianti senza phylum non ha segnale batterico: escluderlo e'
    un fatto, non una scelta di soglia.
    """
    scelta = _campione_e_varianti(finale_calcolata[0], tmp_path)
    run, fase, contesto = _s13_su_decontaminato_modificato(
        decontam_calcolata, tmp_path,
        f"m <- otu_table(ps); m[rownames(m) != {json.dumps(scelta['variante'])}, "
        f"{json.dumps(scelta['campione'])}] <- 0L; otu_table(ps) <- m",
        filt={"exclude_taxa": ["Chloroplast", "Mitochondria", "Eukaryota", scelta["genere"]]},
    )
    fase.calcola(contesto)
    assert [d.codice for d in contesto.degradazioni] == ["E-S13-03"]
    esclusione = [e for e in _tsv(run.albero.cartella(Fase.FINAL) / "esclusioni.tsv")
                  if e["accession"] == scelta["campione"]]
    assert [e["filtro"] for e in esclusione] == ["tassonomico"]


def test_un_campione_svuotato_dal_filtro_di_prevalenza_ferma_con_e_s13_02(
    bioc, finale_calcolata, decontam_calcolata, tmp_path
):
    """
    **Obiettivo**: Verificare che un campione la cui sola variante non e'
    presente con almeno ``prev.min_count`` letture in nessun campione tenuto
    resti vuoto dopo il filtro di prevalenza, e che S13 si fermi con
    ``E-S13-02``, di revisione umana, lasciando ``esclusioni.tsv``.

    **Razionale scientifico e sistemistico**: Il campione aveva letture
    batteriche, solo rare: toglierlo dipende dalle soglie, e la decisione non
    e' automatica.
    """
    scelta = _campione_e_varianti(finale_calcolata[0], tmp_path)
    run, fase, contesto = _s13_su_decontaminato_modificato(
        decontam_calcolata, tmp_path,
        f"m <- otu_table(ps); m[, {json.dumps(scelta['campione'])}] <- 0L; "
        f"m[{json.dumps(scelta['assente'])}, {json.dumps(scelta['campione'])}] <- 1L; otu_table(ps) <- m",
    )
    with pytest.raises(ErrorePipeline) as info:
        fase.calcola(contesto)
    assert info.value.codice == "E-S13-02"
    assert info.value.categoria.value == "revisione_umana"
    assert scelta["campione"] in info.value.dettaglio
    righe = _tsv(run.albero.cartella(Fase.FINAL) / "esclusioni.tsv")
    assert [e["filtro"] for e in righe if e["accession"] == scelta["campione"]] == ["prevalenza"]


def test_i_campioni_con_poche_letture_finali_escono_senza_fermare(bioc, finale_calcolata, tmp_path):
    """
    **Obiettivo**: Verificare che, con ``qc.min_reads_final`` alla mediana delle
    letture finali, la catena si concluda e i campioni sotto la soglia escano
    dall'oggetto finale con il filtro ``letture_finali`` e il motivo.

    **Razionale scientifico e sistemistico**: Un singolo campione povero non
    deve fermare l'intera esecuzione, come non la ferma un bianco azzerato: si
    gestisce come filtro, per campione.
    """
    run, _ = finale_calcolata
    finali = sorted(int(r["letture"]) for r in _tsv(run.albero.cartella(Fase.FINAL) / "letture_finali.tsv")
                    if int(r["letture"]) > 0)
    mediana = finali[len(finali) // 2]
    copia = copia_esecuzione(finale_calcolata, tmp_path, qc={"min_reads_final": mediana})
    assert copia.valuta().situazioni[Passo.S13].stato is StatoPasso.DA_ESEGUIRE
    esito = Esecutore(copia, fino_a=Passo.S14).esegui()
    assert esito.conclusione is Conclusione.COMPLETATA
    poveri = [e for e in _tsv(copia.albero.cartella(Fase.FINAL) / "esclusioni.tsv")
              if e["filtro"] == "letture_finali"]
    assert len(poveri) == sum(1 for n in finali if n < mediana)
    assert all("qc.min_reads_final" in e["motivo"] for e in poveri)


def test_una_frazione_trattenuta_sotto_la_soglia_ferma_con_e_s14_01(bioc, finale_calcolata, tmp_path):
    """
    **Obiettivo**: Verificare che, con ``qc.min_frac_reads_retained`` sopra la
    frazione delle letture trattenute dall'insieme dei campioni finali, S14 si
    fermi con ``E-S14-01`` senza rifare S13.

    **Razionale scientifico e sistemistico**: La frazione misura l'insieme: se
    i filtri hanno tolto la maggior parte del segnale dei campioni tenuti, il
    difetto e' della catena, non di un campione.
    """
    copia = copia_esecuzione(finale_calcolata, tmp_path, qc={"min_frac_reads_retained": 0.9999})
    assert copia.valuta().situazioni[Passo.S13].stato is StatoPasso.COMPLETATA
    esito = Esecutore(copia, fino_a=Passo.S14).esegui()
    assert esito.conclusione is Conclusione.ARRESTATA
    assert (esito.punto.passo, esito.punto.codice) == (Passo.S14, "E-S14-01")
    assert "qc.min_frac_reads_retained" in esito.punto.dettaglio


def test_un_oggetto_filtrato_alterato_ferma_s14_con_e_s14_01(bioc, finale_calcolata, tmp_path):
    """
    **Obiettivo**: Verificare che, con ``ps_filtrato.rds`` diverso da quello
    registrato nel manifesto di S13, il calcolo di S14 si fermi con
    ``E-S14-01`` prima di serializzare.

    **Razionale scientifico e sistemistico**: S14 consegna l'oggetto che S13 ha
    prodotto e registrato, non un file qualunque con lo stesso nome.
    """
    copia = copia_esecuzione(finale_calcolata, tmp_path)
    filtrato = copia.albero.cartella(Fase.FINAL) / "ps_filtrato.rds"
    filtrato.write_bytes(filtrato.read_bytes() + b"\0")
    fase = copia.fase(Passo.S14)
    contesto = copia.contesto(Passo.S14, copia.valuta()).ristretto(fase.parametri)
    with pytest.raises(ErrorePipeline) as info:
        fase.calcola(contesto)
    assert info.value.codice == "E-S14-01"
    assert "ps_filtrato.rds" in info.value.dettaglio


# --------------------------------------------------------------------------- #
# 4. I test sui dati reali su un solo processo                                 #
# --------------------------------------------------------------------------- #


def test_i_test_sui_dati_reali_girano_su_un_solo_processo(tmp_path):
    """
    **Obiettivo**: Verificare che, con le opzioni predefinite della suite
    (``addopts``: quattro processi, ``--dist loadgroup``), i 22 test sui dati
    reali di ``test_w07_dati_reali.py`` girino tutti sullo stesso processo, nel
    gruppo ``dati_reali``, insieme a test di un altro modulo che non vi
    appartengono.

    **Razionale scientifico e sistemistico**: Sul dataset completo i test reali
    lanciano catene con decine di processi R: distribuiti, giravano in
    parallelo fra loro e la suite saturava la memoria della macchina.
    """
    ambiente = {k: v for k, v in os.environ.items() if k != "AMPLICON16S_CONFIG_DATI_REALI"}
    radice = Path(__file__).resolve().parents[1]
    esito = subprocess.run(
        [sys.executable, "-m", "pytest", "-p", "no:cacheprovider", "-v", "--color=no",
         "--basetemp", str(tmp_path / "base"),
         "tests/test_w07_dati_reali.py", "tests/test_w01_w02_import.py"],
        cwd=radice, env=ambiente, capture_output=True, text=True, check=False,
    )
    assert esito.returncode == 0, esito.stdout[-2000:]
    righe = re.findall(r"^\[(gw\d+)\] \[ *\d+%\] (\w+) (\S+)", esito.stdout, re.MULTILINE)
    reali = {nodo: processo for processo, _, nodo in righe if "test_w07_dati_reali.py" in nodo}
    altri = [nodo for _, _, nodo in righe if "test_w01_w02_import.py" in nodo]
    assert len(reali) == 22 and altri
    assert all(nodo.endswith("@dati_reali") for nodo in reali)
    assert len(set(reali.values())) == 1, reali
    assert not any(nodo.endswith("@dati_reali") for nodo in altri)


# --------------------------------------------------------------------------- #
# Dataset completo, nel container                                              #
# --------------------------------------------------------------------------- #


@pytest.mark.dati_reali
def test_la_catena_arriva_alla_fine_sul_dataset_completo(bioc, catena_reale):
    """
    **Obiettivo**: Verificare che sul dataset completo la catena si concluda con
    S14, e riportare le misure dei filtri e le dimensioni dell'oggetto finale.

    **Razionale scientifico e sistemistico**: E' il risultato dell'intera
    pipeline.
    """
    run, esito = catena_reale
    assert esito.conclusione is Conclusione.COMPLETATA
    for r in esito.eseguite:
        if r.passo in (Passo.S13, Passo.S14):
            print(f"\n{r.passo}: {r.secondi} s, {dict(r.metriche)}")
    riepilogo = _riepilogo(run)
    assert riepilogo["campioni"]["classi_finali"]["controllo_negativo"] == 0
    assert riepilogo["campioni"]["finali"] > 0
