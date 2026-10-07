r"""Suite di test della settimana 19: validazione della corsa dai controlli positivi, fase S11.

1. Inquadramento nel Piano Operativo
------------------------------------
Settimana 19 (W19), Fase F5 (validazione della corsa dai controlli positivi
S11 in ``11_controls/``: curva KatharoSeq, soglia di profondita' e conformita'
dei controlli).

2. Moduli sorgente coperti
--------------------------
* ``src/amplicon16s/steps/s11_controls.py``, ``R/11_controls.R``, ``R/lib/katharoseq.R``
* ``src/amplicon16s/config/schema.py`` (gruppo ``katharoseq``, ``ctrl.min_positives``,
  ``ctrl.min_positive_pass_frac``, ``ctrl.positive_gate``, ``qc.min_reads_mode``)
* ``src/amplicon16s/errors/catalog.py`` (``E-S11-02``, ``E-S11-03``, ``E-S11-04``)
* ``src/amplicon16s/steps/s00_validate.py`` (le sole chiavi di ``ctrl`` che S0 legge)

3. Cosa valuta questo file
--------------------------
- la sigmoide allosterica ritrova i parametri da cui sono generati i punti, la
  bonta' e' ``1 - SS_res / SS_tot`` sulla fedelta', la soglia e' quella
  dell'equazione dichiarata; una soglia senza osservazioni fra il punto medio
  della curva e la soglia non e' determinata;
- la conformita' dipende dal livello di concentrazione: a bassa concentrazione
  un controllo dominato dai contaminanti come gli altri del suo livello e'
  conforme; un controllo con fedelta' bassa o profondita' anomala per il suo
  livello no; un livello con meno di ``ctrl.min_positives`` controlli non e'
  valutabile;
- S10 non e' rifatta dai parametri di S11, e S0 non dipende da essi;
- S11 sul sottoinsieme di prova: curva aggregata e per piastra con la bonta'
  di entrambe, scelta motivata, soglia per piastra con lo stadio a cui si
  applica, ripiego con ``E-S11-02`` e motivo per le piastre senza curva,
  stessi byte in due esecuzioni;
- forzando ``katharoseq.min_r2`` sopra ogni bonta' tutte le piastre ripiegano
  su ``qc.min_reads_raw``, sulle letture grezze, con ``E-S11-02`` e il motivo;
  senza la colonna dei livelli anche; con ``qc.min_reads_mode: fixed`` si usa
  il ripiego senza degradazione;
- su un oggetto costruito: controlli conformi alla loro concentrazione non
  fanno scattare nulla; una piastra di controlli non conformi fa scattare
  ``E-S11-03`` con ``ctrl.positive_gate`` vero e l'avviso ``E-S11-04`` con
  falso;
- sul dataset completo: soglia, curve, conformita' e campioni biologici sotto
  la soglia.

4. Comandi Bash e scenari di esecuzione
---------------------------------------
    ``<immagine>`` e' l'immagine del container della pipeline; quella corrente
    e' indicata in ``test.txt``, sezione 1.3.

    1. Modalità locale standard (R di base con jsonlite, senza Bioconductor né
       dati reali):
       pytest tests/test_w19_controlli.py -v

    2. Modalità container Docker standard (sottoinsieme ridotto con
       R/Bioconductor):
       docker run --rm \
         -e PYTHONPATH=/app/src \
         -e AMPLICON16S_R_DIR=/app/R \
         -v "$(pwd)":/app \
         -w /app \
         <immagine> \
         pytest -o cache_dir=/tmp/.pytest_cache tests/test_w19_controlli.py -v

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
         pytest -o cache_dir=/tmp/.pytest_cache tests/test_w19_controlli.py -v

5. Risultato atteso
-------------------
Vedi ``test.txt``, scheda W19.

6. Razionale scientifico e sistemistico
---------------------------------------
- La soglia di profondita' e' un risultato dei controlli, non una scelta: va
  scritta con lo stadio delle letture a cui si applica, perche' il ripiego
  sulle letture grezze e la soglia sulle letture senza chimere non sono
  confrontabili.
- Ai livelli piu' diluiti della serie i contaminanti prevalgono per
  costruzione: giudicare un controllo senza la sua concentrazione darebbe un
  falso allarme proprio dove la curva prevede la fedelta' bassa.
- Un ripiego non deve passare inosservato: il motivo resta negli artefatti e
  nel manifesto.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
import os
import subprocess
from pathlib import Path

import pytest
from conftest import copia_esecuzione
from sottoinsieme import TAXON_SINTETICO, config_ridotta, motivo_pacchetti_r_assenti

from amplicon16s.config import defaults
from amplicon16s.config.resolve import risolvi
from amplicon16s.errors.exceptions import ErrorePipeline
from amplicon16s.gates.g01_g15 import Contesto
from amplicon16s.gates.registry import esegui_gate
from amplicon16s.io_layer.artifacts import Fase
from amplicon16s.logging.logger import chiudi
from amplicon16s.rbridge.payload import PREFISSO
from amplicon16s.rbridge.runner import cartella_r, trova_rscript
from amplicon16s.runner.executor import Conclusione, Esecutore
from amplicon16s.runner.graph import Passo
from amplicon16s.runner.project import StatoPasso, passi_realizzati
from amplicon16s.steps.s10_phyloseq import NOME_OGGETTO
from amplicon16s.steps.s11_controls import colonna_dei_livelli


@pytest.fixture(autouse=True)
def uscite_pulite():
    """Chiude le uscite del log prima e dopo ogni test, perché nessun handler resti
    aperto sulla cartella temporanea.
    """
    chiudi()
    yield
    chiudi()


_MOTIVO_R = motivo_pacchetti_r_assenti("jsonlite")
_MOTIVO_BIOC = motivo_pacchetti_r_assenti(
    "dada2", "ggplot2", "ShortRead", "jsonlite", "phyloseq", "Biostrings"
)


@pytest.fixture
def r():
    """Richiede R con jsonlite: salta senza, ma in CI fallisce."""
    if _MOTIVO_R is not None:
        if os.environ.get("AMPLICON16S_RICHIEDI_R") == "1":
            pytest.fail(f"R e' richiesto in questo ambiente: {_MOTIVO_R}")
        pytest.skip(_MOTIVO_R)


@pytest.fixture
def bioc():
    """Richiede R con phyloseq e i pacchetti delle fasi a monte: salta senza, ma in
    CI fallisce.
    """
    if _MOTIVO_BIOC is not None:
        if os.environ.get("AMPLICON16S_RICHIEDI_BIOC") == "1":
            pytest.fail(f"Bioconductor e' richiesto in questo ambiente: {_MOTIVO_BIOC}")
        pytest.skip(_MOTIVO_BIOC)


def _tsv(percorso: Path) -> list[dict[str, str]]:
    """Le righe di una tabella separata da tabulazioni."""
    with open(percorso, encoding="utf-8", newline="") as file:
        return list(csv.DictReader(file, delimiter="\t"))


def _impronte(cartella: Path) -> dict[str, str]:
    """L'MD5 dei file della cartella, esclusi i file del ponte e i manifesti."""
    return {
        p.name: hashlib.md5(p.read_bytes()).hexdigest()
        for p in sorted(cartella.iterdir())
        if not p.name.startswith((PREFISSO, "manifest"))
    }


def _r(codice: str, cartella: Path) -> object:
    """Esegue codice R con le funzioni di ``R/lib/katharoseq.R`` caricate; il codice
    scrive un oggetto JSON in ``uscita``.
    """
    uscita = cartella / "uscita.json"
    script = cartella / "prova.R"
    libreria = cartella_r() / "lib" / "katharoseq.R"
    script.write_text(
        f"source({json.dumps(str(libreria))})\nuscita <- {json.dumps(str(uscita))}\n{codice}\n",
        encoding="utf-8",
    )
    subprocess.run([str(trova_rscript()), "--vanilla", str(script)], check=True,
                   capture_output=True)
    return json.loads(uscita.read_text(encoding="utf-8"))


def _soglia(run) -> dict:
    """``soglia.json`` di S11."""
    return json.loads((run.albero.cartella(Fase.CONTROLS) / "soglia.json").read_text())


def _riepilogo(run) -> dict:
    """``riepilogo.json`` di S11."""
    return json.loads((run.albero.cartella(Fase.CONTROLS) / "riepilogo.json").read_text())


# --------------------------------------------------------------------------- #
# 1. Parametri e dipendenze                                                    #
# --------------------------------------------------------------------------- #


def test_i_parametri_di_s11_hanno_i_valori_del_piano(tmp_path):
    """
    **Obiettivo**: Verificare che lo schema porti i parametri di S11 con i
    valori del piano, e che la colonna dei livelli sia un parametro.

    **Razionale scientifico e sistemistico**: Il livello di diluizione viene da
    una colonna dei dati indicata in configurazione, non dai nomi dei campioni.
    """
    config = config_ridotta(tmp_path)
    k = config.katharoseq
    assert (k.cell_count_column, k.collapse_rank, k.curve_model) == (
        "katharoseq_cell_count", "Genus", "allosteric_sigmoid"
    )
    # Nel riferimento sintetico Variovorax ha un nome inventato.
    assert k.target_taxon == TAXON_SINTETICO
    assert defaults.ESEMPIO_OSD734["katharoseq.target_taxon"].valore == "Variovorax"
    assert (k.target_sensitivity, k.min_r2, k.read_stage) == (0.90, 0.80, "nonchimeric")
    assert config.qc.min_reads_mode == "katharoseq_if_available"
    assert (config.ctrl.min_positives, config.ctrl.min_positive_pass_frac,
            config.ctrl.positive_gate) == (3, 0.75, False)


@pytest.mark.parametrize(
    "variazione",
    [{"ctrl": {"min_positives": 5}}, {"ctrl": {"positive_gate": True}},
     {"katharoseq": {"min_r2": 0.5}}, {"qc": {"min_reads_mode": "fixed"}}],
    ids=["min_positives", "positive_gate", "min_r2", "min_reads_mode"],
)
def test_i_parametri_di_s11_non_toccano_s0_ne_s10(tmp_path, variazione):
    """
    **Obiettivo**: Verificare che cambiare un parametro di S11 lasci invariata
    l'impronta di configurazione di S0 e di S10, e cambi quella di S11.

    **Razionale scientifico e sistemistico**: S0 dichiarava l'intero gruppo
    ``ctrl``: le chiavi nuove di S11 avrebbero rifatto la catena intera. Ora S0
    dichiara le sole quattro chiavi di ``ctrl`` che legge.
    """
    prima = risolvi(config_ridotta(tmp_path))
    dopo = risolvi(config_ridotta(tmp_path, **variazione))
    passi = passi_realizzati()
    for passo in (Passo.S0, Passo.S10):
        assert passi[passo].calcolata_su(prima, {}) == passi[passo].calcolata_su(dopo, {})
    s11 = passi[Passo.S11]
    assert s11.calcolata_su(prima, {}) != s11.calcolata_su(dopo, {})


def test_la_colonna_dei_livelli_si_cerca_fra_quelle_del_lotto():
    """
    **Obiettivo**: Verificare che il nome nell'oggetto della colonna dei livelli
    si trovi per nome originale fra le colonne venute dal file di arricchimento,
    e che una colonna assente dia ``None``.

    **Razionale scientifico e sistemistico**: Il nome nell'oggetto e' scelto da
    S10; la corrispondenza scritta e' l'unico legame con la colonna originale.
    """
    corrispondenza = [
        {"colonna": "materiale", "origine": "tabella campioni di studio (io.study_table)",
         "colonna_originale": "katharoseq_cell_count"},
        {"colonna": "katharoseq_cell_count", "origine": "file di arricchimento (io.batch_table)",
         "colonna_originale": "katharoseq_cell_count"},
    ]
    assert colonna_dei_livelli(corrispondenza, "katharoseq_cell_count") == "katharoseq_cell_count"
    assert colonna_dei_livelli(corrispondenza[:1], "katharoseq_cell_count") is None


# --------------------------------------------------------------------------- #
# 2. La curva e la conformita', in R                                           #
# --------------------------------------------------------------------------- #


def test_la_curva_ritrova_i_parametri_e_la_soglia_e_quella_dichiarata(r, tmp_path):
    """
    **Obiettivo**: Verificare che ``adatta_curva`` ritrovi h e x50 da punti
    generati con la sigmoide allosterica piu' un rumore fisso, che la bonta'
    sia ``1 - SS_res / SS_tot`` calcolata sulla fedelta', che ``k' = x50^h``
    e che la soglia sia ``10^(x50 (s / (1 - s))^(1 / h))``.

    **Razionale scientifico e sistemistico**: L'equazione e la bonta' sono
    dichiarate nel codice: il test le ricalcola da fuori.
    """
    esito = _r("""
n <- 10^seq(2.5, 5.3, length.out = 16)
rumore <- rep(c(0.02, -0.03, 0.01, -0.01), 4)
f <- pmin(pmax(1 / (1 + (3.9 / log10(n))^18) + rumore, 0), 1)
cv <- adatta_curva(n, f, 0.9)
jsonlite::write_json(c(cv, list(n = n, f = f)), uscita, auto_unbox = TRUE, digits = NA)
""", tmp_path)
    assert esito["converge"] is True
    assert abs(esito["h"] - 18) < 3 and abs(esito["x50"] - 3.9) < 0.05
    assert esito["k"] == pytest.approx(esito["x50"] ** esito["h"], rel=1e-9)
    previste = [1 / (1 + (esito["x50"] / math.log10(n)) ** esito["h"]) for n in esito["n"]]
    media = sum(esito["f"]) / len(esito["f"])
    r2 = 1 - sum((f - p) ** 2 for f, p in zip(esito["f"], previste)) / sum(
        (f - media) ** 2 for f in esito["f"])
    assert esito["r2"] == pytest.approx(r2, rel=1e-9)
    assert esito["soglia"] == pytest.approx(10 ** (esito["x50"] * 9 ** (1 / esito["h"])), rel=1e-9)


def test_con_meno_di_tre_punti_o_fedelta_costante_nessuna_curva(r, tmp_path):
    """
    **Obiettivo**: Verificare che con due punti, o con la fedelta' identica in
    tutti, ``adatta_curva`` non converga e dica perche'.

    **Razionale scientifico e sistemistico**: Due parametri non si stimano da due
    punti, ne' una sigmoide da una retta orizzontale: sono motivi di ripiego.
    """
    esito = _r("""
jsonlite::write_json(list(
  pochi = adatta_curva(c(100, 1000), c(0.1, 0.9), 0.9)$motivo,
  piatta = adatta_curva(c(100, 1000, 10000), c(0.5, 0.5, 0.5), 0.9)$motivo
), uscita, auto_unbox = TRUE)
""", tmp_path)
    assert "2 punti" in esito["pochi"]
    assert "identica" in esito["piatta"]


def test_una_soglia_in_un_tratto_senza_osservazioni_non_e_determinata(r, tmp_path):
    """
    **Obiettivo**: Verificare che, con i controlli della piastra 1 del dataset di
    riferimento (un controllo a fedelta' 0,96 e 99.521 letture, gli altri sotto
    0,03 e 49.417 letture), la curva converga con bonta' alta ma la soglia
    risulti non determinata, perche' nessuna osservazione cade fra il punto
    medio della curva e la soglia; e che una serie con osservazioni dentro la
    transizione dia una soglia determinata.

    **Razionale scientifico e sistemistico**: Nel tratto vuoto qualunque
    posizione del gradino darebbe un adattamento equivalente: la soglia
    verrebbe dalla forma del modello, non dai dati.
    """
    esito = _r("""
gradino_n <- c(99521, 49417, 26423, 24186, 30195)
gradino_f <- c(0.961, 0.023, 0.008, 0.002, 0.004)
cv <- adatta_curva(gradino_n, gradino_f, 0.9)
serie_n <- c(100106, 64074, 29682, 8981, 4644, 2946, 2027, 2803)
serie_f <- c(0.978, 0.967, 0.910, 0.817, 0.431, 0.085, 0.002, 0.000)
cs <- adatta_curva(serie_n, serie_f, 0.9)
jsonlite::write_json(list(r2 = cv$r2, gradino = soglia_determinata(gradino_n, cv),
                          serie = soglia_determinata(serie_n, cs)),
                     uscita, auto_unbox = TRUE, digits = NA)
""", tmp_path)
    assert esito["r2"] > 0.99
    assert esito["gradino"]["determinata"] is False
    assert esito["gradino"]["osservazioni"] == 0
    assert "non determinata dai dati" in esito["gradino"]["motivo"]
    assert "49417" in esito["gradino"]["motivo"]
    assert esito["serie"]["determinata"] is True
    assert esito["serie"]["osservazioni"] >= 1


#: Dieci piastre con la serie di otto livelli, come sul dataset di riferimento:
#: la fedelta' scende con la concentrazione e ai livelli piu' diluiti prevalgono
#: i contaminanti; la profondita' scende con la concentrazione.
SERIE = """
cellule <- 2e6 / 5^(0:7)
fed_livello <- c(0.98, 0.96, 0.91, 0.78, 0.40, 0.08, 0.02, 0.003)
prof_livello <- c(120000, 64000, 35000, 9000, 4500, 3000, 2200, 2500)
piastre <- 1:10
livello <- rep(cellule, times = length(piastre))
scarto <- rep(seq(-0.02, 0.02, length.out = length(piastre)), each = 8)
fedelta <- pmax(rep(fed_livello, times = length(piastre)) + scarto * rep(fed_livello, times = 10), 0)
profondita <- round(rep(prof_livello, times = length(piastre)) * (1 + 5 * scarto))
"""


def test_a_bassa_concentrazione_i_contaminanti_non_rendono_non_conforme(r, tmp_path):
    """
    **Obiettivo**: Verificare che, in dieci piastre che si comportano come
    previsto (fedelta' sotto 0,1 ai tre livelli piu' diluiti), ogni controllo
    sia conforme: 80 su 80, mentre per dominanza del taxon atteso lo sarebbero
    meno di 5 su 8 livelli.

    **Razionale scientifico e sistemistico**: E' il falso allarme da evitare: ai
    livelli diluiti la curva stessa prevede che il taxon atteso non domini.
    """
    esito = _r(SERIE + """
conf <- conformita(livello, fedelta, profondita, 3L)
jsonlite::write_json(list(esiti = as.list(table(conf$esito)), dominati = sum(fedelta > 0.5)),
                     uscita, auto_unbox = TRUE)
""", tmp_path)
    assert esito["esiti"] == {"conforme": 80}
    assert esito["dominati"] == 40


@pytest.mark.parametrize(
    ("modifica", "motivo"),
    [
        ("fedelta[17] <- 0.2", "fedelta' bassa per il livello"),       # piastra 3, 2.000.000 cellule
        ("profondita[2] <- 259", "profondita' bassa per il livello"),  # come KATHARO.P4.2
        ("profondita[8] <- 30000", "profondita' alta per il livello"),  # 25,6 cellule, 30.000 letture
    ],
    ids=["fedelta_bassa", "profondita_bassa", "profondita_alta"],
)
def test_un_controllo_costruito_per_non_essere_conforme_e_intercettato(r, tmp_path, modifica, motivo):
    """
    **Obiettivo**: Verificare che un solo controllo che si discosta dagli altri
    del suo livello (fedelta' bassa a 2.000.000 di cellule, 259 letture a
    400.000 cellule, 30.000 letture a 25,6 cellule) sia l'unico non conforme,
    con il motivo.

    **Razionale scientifico e sistemistico**: Non conforme vuol dire diverso da
    cio' che ci si attende per la sua concentrazione, in fedelta' o in
    profondita'.
    """
    esito = _r(SERIE + modifica + """
conf <- conformita(livello, fedelta, profondita, 3L)
jsonlite::write_json(list(non = which(conf$esito == "non conforme"),
                          motivo = conf$motivo[conf$esito == "non conforme"]),
                     uscita, auto_unbox = FALSE)
""", tmp_path)
    indice = int(modifica.split("[")[1].split("]")[0])
    assert esito["non"] == [indice]
    assert esito["motivo"] == [motivo]


def test_un_livello_con_pochi_controlli_non_e_valutabile(r, tmp_path):
    """
    **Obiettivo**: Verificare che con due soli controlli per livello, sotto
    ``ctrl.min_positives``, nessun controllo sia giudicato, e che uno senza
    livello sia non valutabile.

    **Razionale scientifico e sistemistico**: Un confronto fra due valori non
    dice quale dei due e' anomalo.
    """
    esito = _r("""
conf <- conformita(c(1e6, 1e6, NA), c(0.9, 0.1, 0.5), c(1e4, 1e4, 1e4), 3L)
jsonlite::write_json(conf, uscita)
""", tmp_path)
    assert [c["esito"] for c in esito] == ["non valutabile"] * 3
    assert "meno di 3" in esito[0]["motivo"]
    assert esito[2]["motivo"] == "senza livello di concentrazione"


# --------------------------------------------------------------------------- #
# 3. S11 sul sottoinsieme di prova, nel container                               #
# --------------------------------------------------------------------------- #


def test_s11_adatta_le_curve_e_scrive_la_soglia_con_il_suo_stadio(bioc, controlli_calcolati):
    """
    **Obiettivo**: Verificare che S11 sulla versione ridotta si concluda dopo
    S10, che ``curve.tsv`` registri la curva aggregata e quelle per piastra con
    la bonta', che la scelta sia motivata in ``soglia.json``, e che ogni
    soglia porti lo stadio: la piastra 10, con la serie completa, ha la soglia
    della propria curva sulle letture ``nonchimeric``; le altre, senza una
    serie, ripiegano su ``qc.min_reads_raw`` sulle letture ``raw``, con
    ``E-S11-02`` e il motivo nel manifesto.

    **Razionale scientifico e sistemistico**: La versione ridotta ha la serie
    di otto livelli della sola piastra 10 e un positivo della piastra 4: e' il
    caso misto, curva dove c'e' e ripiego dichiarato dove manca.
    """
    run, esito = controlli_calcolati
    assert esito.conclusione is Conclusione.COMPLETATA
    assert [r.passo for r in esito.eseguite] == [Passo.S11]
    curve = {c["modello"]: c for c in _tsv(run.albero.cartella(Fase.CONTROLS) / "curve.tsv")}
    assert {"aggregato", "piastra 10"} <= set(curve)
    assert curve["piastra 10"]["valida"] == "si"
    assert float(curve["piastra 10"]["r2"]) >= 0.8
    soglia = _soglia(run)
    assert soglia["scelta"] in ("per_piastra", "aggregato")
    assert soglia["motivo_scelta"]
    assert set(soglia["stadi"]) == {"nonchimeric", "raw"}
    per_piastra = soglia["per_piastra"]
    assert set(per_piastra) == {"2", "4", "8", "10"}
    for piastra, voce in per_piastra.items():
        assert voce["stadio"] in ("nonchimeric", "raw")
        if voce["stadio"] == "raw":
            assert voce["valore"] == 1000 and voce["motivo"]
        else:
            assert voce["origine"].startswith("curva") and voce["valore"] > 0
    if soglia["scelta"] == "per_piastra":
        assert per_piastra["10"]["valore"] == int(curve["piastra 10"]["soglia"])
        assert per_piastra["10"]["stadio"] == "nonchimeric"
    manifesto = run.albero.manifesto_passo(Passo.S11, Fase.CONTROLS)
    assert [d["codice"] for d in manifesto.degradazioni] == (
        ["E-S11-02"] if soglia["ripiego"]["degradazione"] else []
    )
    print(f"\nS11 sulla versione ridotta: {dict(esito.eseguite[0].metriche)}")
    print(f"scelta: {soglia['scelta']} ({soglia['motivo_scelta']})")


def test_la_misura_sui_campioni_usa_la_soglia_del_suo_stadio(bioc, controlli_calcolati):
    """
    **Obiettivo**: Verificare che ``profondita_campioni.tsv`` abbia una riga per
    campione dell'inventario, con la soglia e lo stadio della sua piastra, e
    che ``sotto_soglia`` confronti la soglia con le letture di quello stadio:
    grezze per il ripiego, senza chimere per la soglia derivata.

    **Razionale scientifico e sistemistico**: Le due grandezze non sono
    confrontabili: una soglia applicata alle letture sbagliate darebbe un
    filtro plausibile e sbagliato.
    """
    run, _ = controlli_calcolati
    righe = _tsv(run.albero.cartella(Fase.CONTROLS) / "profondita_campioni.tsv")
    assert [r["accession"] for r in righe] == list(run.valuta().inventario.accessioni)
    per_piastra = _soglia(run)["per_piastra"]
    for r in righe:
        voce = per_piastra[r["piastra"]]
        assert (int(r["soglia"]), r["stadio"]) == (voce["valore"], voce["stadio"])
        letture = int(r["letture_grezze"] if r["stadio"] == "raw" else r["letture_nonchimeric"])
        assert r["sotto_soglia"] == ("si" if letture < voce["valore"] else "no")
    riepilogo = _riepilogo(run)
    assert riepilogo["biologici_sotto_soglia"] == sum(
        1 for r in righe if r["classe"] == "biologico" and r["sotto_soglia"] == "si")


def test_s11_da_gli_stessi_byte_in_due_esecuzioni(bioc, controlli_calcolati, tmp_path):
    """
    **Obiettivo**: Verificare che S11, rieseguita sugli stessi artefatti, dia
    in ``11_controls`` gli stessi byte.

    **Razionale scientifico e sistemistico**: L'adattamento non usa numeri
    casuali e parte da valori iniziali ricavati dai dati: due esecuzioni danno
    la stessa curva e la stessa soglia.
    """
    run, _ = controlli_calcolati
    attese = _impronte(run.albero.cartella(Fase.CONTROLS))
    assert "soglia.json" in attese
    copia = copia_esecuzione(controlli_calcolati, tmp_path)
    copia.albero.rimuovi_manifesto_passo(Passo.S11, Fase.CONTROLS)
    esito = Esecutore(copia, fino_a=Passo.S11).esegui()
    assert [r.passo for r in esito.eseguite] == [Passo.S11]
    assert _impronte(copia.albero.cartella(Fase.CONTROLS)) == attese


def test_la_colonna_dei_livelli_assente_e_respinta_in_s0(tmp_path):
    """
    **Obiettivo**: Verificare che una ``katharoseq.cell_count_column`` che non
    e' fra le colonne portate nell'oggetto sia respinta da G15 con
    ``E-G15-12``; che, elencata in ``out.batch_columns`` ma assente dal file
    del lotto, sia respinta da G08 con ``E-S0-08``; e che la colonna si trovi
    nell'oggetto da qualunque delle due tabelle venga.

    **Razionale scientifico e sistemistico**: Senza la colonna dei livelli S11
    ripiegava sulla soglia fissa dopo l'intera catena di calcolo: una colonna
    dichiarata e inesistente e' un errore della configurazione, e si scopre
    leggendo le intestazioni.
    """
    config = config_ridotta(tmp_path, katharoseq={"cell_count_column": "cellule_inesistenti"})
    esito = esegui_gate("G15", Contesto(config))
    assert not esito.superato and {v.codice for v in esito.violazioni} == {"E-G15-12"}
    assert "katharoseq.cell_count_column" in esito.violazioni[0].dettaglio
    assert "cellule_inesistenti" in esito.violazioni[0].dettaglio

    colonne = [*config.out.batch_columns, "cellule_inesistenti"]
    elencata = config_ridotta(tmp_path, katharoseq={"cell_count_column": "cellule_inesistenti"},
                              out={"batch_columns": colonne})
    assert esegui_gate("G15", Contesto(elencata)).superato
    esito = esegui_gate("G08", Contesto(elencata))
    assert not esito.superato and {v.codice for v in esito.violazioni} == {"E-S0-08"}
    assert "cellule_inesistenti" in esito.violazioni[0].dettaglio

    dal_lotto = {"colonna": "cellule", "origine": "file di arricchimento (io.batch_table)",
                 "colonna_originale": "Cellule seminate", "valore": "valore originale"}
    dallo_studio = {**dal_lotto, "origine": "tabella campioni di studio (io.study_table)"}
    inventario = {**dal_lotto, "colonna": "piastra", "colonna_originale": "Cellule seminate"}
    assert colonna_dei_livelli([dal_lotto], "Cellule seminate") == "cellule"
    assert colonna_dei_livelli([dallo_studio], "Cellule seminate") == "cellule"
    assert colonna_dei_livelli([inventario], "Cellule seminate") is None
    assert colonna_dei_livelli([dal_lotto], None) is None


@pytest.mark.parametrize(
    ("sovrascrivi", "nel_motivo", "degrada"),
    [
        ({"katharoseq": {"min_r2": 1.0}}, "katharoseq.min_r2", True),
        ({"ctrl": {"min_positives": 50}}, "ctrl.min_positives", True),
        ({"qc": {"min_reads_mode": "fixed"}}, "fixed", False),
    ],
    ids=["bonta_insufficiente", "controlli_insufficienti", "modo_fixed"],
)
def test_il_ripiego_e_registrato_con_il_motivo(bioc, oggetto_calcolato, tmp_path,
                                               sovrascrivi, nel_motivo, degrada):
    """
    **Obiettivo**: Verificare che con ``katharoseq.min_r2`` sopra ogni bonta'
    o con meno controlli di ``ctrl.min_positives`` ogni piastra ripieghi su ``qc.min_reads_raw`` sulle
    letture grezze, con ``E-S11-02`` nel manifesto e il motivo in
    ``soglia.json``; e che con ``qc.min_reads_mode: fixed`` si usi il ripiego
    senza degradazione.

    **Razionale scientifico e sistemistico**: Il ripiego e' una degradazione
    automatica: non ferma, ma resta scritto, con il suo motivo e il suo stadio.
    """
    run = copia_esecuzione(oggetto_calcolato, tmp_path, **sovrascrivi)
    assert run.valuta().situazioni[Passo.S10].stato is StatoPasso.COMPLETATA
    esito = Esecutore(run, fino_a=Passo.S11).esegui()
    assert esito.conclusione is Conclusione.COMPLETATA
    soglia = _soglia(run)
    assert soglia["scelta"] == "nessuno"
    assert nel_motivo in soglia["motivo_scelta"]
    assert {v["stadio"] for v in soglia["per_piastra"].values()} == {"raw"}
    assert {v["valore"] for v in soglia["per_piastra"].values()} == {1000}
    assert soglia["ripiego"]["degradazione"] is degrada
    manifesto = run.albero.manifesto_passo(Passo.S11, Fase.CONTROLS)
    assert [d["codice"] for d in manifesto.degradazioni] == (["E-S11-02"] if degrada else [])
    if degrada:
        assert nel_motivo in manifesto.degradazioni[0]["dettaglio"]


#: Un oggetto integrato costruito: tre piastre con la serie di otto livelli che
#: segue la curva, sei campioni biologici; ``contaminate`` sono le posizioni
#: (1-24) dei controlli resi non conformi, con il taxon atteso quasi assente.
COSTRUISCI_OGGETTO = r"""
suppressPackageStartupMessages(library(phyloseq))
cellule <- 2e6 / 5^(0:7)
fed <- c(0.98, 0.96, 0.91, 0.78, 0.40, 0.08, 0.02, 0.003)
prof <- c(120000, 64000, 35000, 9000, 4500, 3000, 2200, 2500)
pos <- sprintf("ERXP%02d", 1:24); bio <- sprintf("ERXB%02d", 1:6)
campioni <- c(pos, bio)
piastra <- c(rep(c("1", "2", "3"), each = 8), rep(c("1", "2", "3"), 2))
n <- c(round(rep(prof, 3) * rep(c(0.97, 1, 1.03), each = 8)), c(500, 4000, 9000, 20000, 40000, 80000))
f <- c(rep(fed, 3) * rep(c(0.99, 1, 1.01), each = 8), rep(0, 6))
f[contaminate] <- 0.01
otu <- rbind(Var = round(n * f), Altro = n - round(n * f))
colnames(otu) <- campioni
tax <- matrix(c("Bacteria", "Bacteria", "Proteobacteria", "Cyanobacteria", "C", "C", "O", "O",
                "Comamonadaceae", "F", taxon, "G"), nrow = 2,
              dimnames = list(c("Var", "Altro"), c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")))
dati <- data.frame(accession = campioni, sample_name = campioni,
                   classe = c(rep("controllo_positivo", 24), rep("biologico", 6)),
                   piastra = piastra,
                   katharoseq_cell_count = c(format(rep(cellule, 3)), rep(NA, 6)),
                   row.names = campioni, stringsAsFactors = FALSE)
ps <- phyloseq(otu_table(otu, taxa_are_rows = TRUE), tax_table(tax), sample_data(dati))
saveRDS(ps, oggetto)
"""


def _s11_su_oggetto_costruito(run, tmp_path, contaminate: list[int]):
    """Sostituisce l'oggetto di S10 con quello costruito ed esegue il calcolo di S11.

    Restituisce il contesto, con le degradazioni registrate.
    """
    oggetto = run.albero.cartella(Fase.PHYLOSEQ) / NOME_OGGETTO
    script = tmp_path / "costruisci.R"
    script.write_text(
        f"oggetto <- {json.dumps(str(oggetto))}\n"
        f"contaminate <- c({', '.join(str(i) for i in contaminate)})\n"
        f"taxon <- {json.dumps(run.config.katharoseq.target_taxon)}\n" + COSTRUISCI_OGGETTO,
        encoding="utf-8",
    )
    subprocess.run([str(trova_rscript()), "--vanilla", str(script)], check=True,
                   capture_output=True)
    fase = run.fase(Passo.S11)
    contesto = run.contesto(Passo.S11, run.valuta())
    ristretto = contesto.ristretto(fase.parametri)
    fase.calcola(ristretto)
    return ristretto


def test_controlli_conformi_alla_concentrazione_non_fanno_scattare_nulla(
    bioc, oggetto_calcolato, tmp_path
):
    """
    **Obiettivo**: Verificare che su un oggetto costruito con tre piastre la cui
    serie segue la curva, con il taxon atteso minoritario ai livelli diluiti,
    S11 dia 24 controlli conformi su 24, nessun ``E-S11-03`` ne' ``E-S11-04``,
    anche con ``ctrl.positive_gate`` vero, e una soglia derivata.

    **Razionale scientifico e sistemistico**: Per dominanza del taxon atteso
    sarebbero conformi 15 su 24, sotto 0,75: la conformita' per concentrazione
    non scatta su un comportamento che la curva prevede.
    """
    run = copia_esecuzione(oggetto_calcolato, tmp_path / "copia", ctrl={"positive_gate": True})
    contesto = _s11_su_oggetto_costruito(run, tmp_path, [])
    assert [d.codice for d in contesto.degradazioni] == []
    riepilogo = _riepilogo(run)
    assert (riepilogo["conformi"], riepilogo["non_conformi"]) == (24, 0)
    assert _soglia(run)["scelta"] in ("per_piastra", "aggregato")


@pytest.mark.parametrize("cancello", [True, False])
def test_controlli_non_conformi_fanno_scattare_e_s11_03_o_l_avviso(
    bioc, oggetto_calcolato, tmp_path, cancello
):
    """
    **Obiettivo**: Verificare che su un oggetto costruito con i sette livelli piu'
    concentrati della piastra 3 quasi privi del taxon atteso (17 conformi su
    24, 0,708), S11 si fermi con ``E-S11-03`` se ``ctrl.positive_gate`` e'
    vero, e prosegua registrando l'avviso ``E-S11-04`` se e' falso; e che i
    non conformi siano proprio quei sette, scritti in ``positivi.tsv``.

    **Razionale scientifico e sistemistico**: Un avviso e un arresto non
    condividono il codice: ``ctrl.positive_gate`` sceglie fra i due.
    """
    run = copia_esecuzione(oggetto_calcolato, tmp_path / "copia", ctrl={"positive_gate": cancello})
    contaminate = list(range(17, 24))
    if cancello:
        with pytest.raises(ErrorePipeline) as info:
            _s11_su_oggetto_costruito(run, tmp_path, contaminate)
        assert info.value.codice == "E-S11-03"
        assert "17 controlli positivi conformi su 24" in info.value.dettaglio
    else:
        contesto = _s11_su_oggetto_costruito(run, tmp_path, contaminate)
        assert [d.codice for d in contesto.degradazioni if d.codice != "E-S11-02"] == ["E-S11-04"]
    positivi = _tsv(run.albero.cartella(Fase.CONTROLS) / "positivi.tsv")
    assert [p["accession"] for p in positivi if p["conformita"] == "non conforme"] == [
        f"ERXP{i:02d}" for i in contaminate
    ]


# --------------------------------------------------------------------------- #
# Dataset completo, nel container                                              #
# --------------------------------------------------------------------------- #


@pytest.mark.dati_reali
def test_s11_sul_dataset_completo(bioc, catena_reale):
    """
    **Obiettivo**: Verificare che sul dataset completo S11 si concluda dopo
    S10 con una soglia derivata dalle curve, senza ``E-S11-03`` ne'
    ``E-S11-04``, e riportarne curve, scelta, conformita' e campioni biologici
    sotto la soglia.

    **Razionale scientifico e sistemistico**: Sono le misure che decidono la
    soglia dei filtri finali di S13.
    """
    run, esito = catena_reale
    assert esito.conclusione is Conclusione.COMPLETATA
    s11 = next(r for r in esito.eseguite if r.passo is Passo.S11)
    print(f"\nS11: {s11.secondi} s, {dict(s11.metriche)}")
    soglia = _soglia(run)
    assert soglia["scelta"] in ("per_piastra", "aggregato")
    riepilogo = _riepilogo(run)
    assert riepilogo["frazione_conformi"] >= run.config.ctrl.min_positive_pass_frac
    manifesto = run.albero.manifesto_passo(Passo.S11, Fase.CONTROLS)
    assert "E-S11-04" not in [d["codice"] for d in manifesto.degradazioni]
