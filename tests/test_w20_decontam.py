r"""Suite di test della settimana 20: decontaminazione dai controlli negativi, fase S12.

1. Inquadramento nel Piano Operativo
------------------------------------
Settimana 20 (W20), Fase F4 (decontaminazione per prevalenza S12 in
``11_controls/``, accanto a S11, con il confronto fra modalita' per piastra e
aggregata).

2. Moduli sorgente coperti
--------------------------
* ``src/amplicon16s/steps/s12_decontam.py``, ``R/12_decontam.R``,
  ``R/lib/decontaminazione.R``
* ``src/amplicon16s/config/schema.py`` (``decontam.batch_combine``,
  ``qc.max_frac_contaminant``), ``src/amplicon16s/errors/catalog.py`` (``E-S12-02``)
* ``src/amplicon16s/runner/graph.py`` (precedenza S12 prima di S13, ``E-S13-01``)

3. Cosa valuta questo file
--------------------------
- la regola di riclassificazione di ``ctrl``: un campione riconosciuto dalla
  colonna e dal valore dichiarati e' un controllo negativo, con il materiale
  originale; la regola cambia l'impronta di S0 e il denominatore della
  prevalenza;
- i parametri di S12, con ``batch_combine`` dichiarato al predefinito di
  decontam e ``mode`` che dichiara la modalita', e la loro indipendenza dalle
  fasi a monte;
- la combinazione delle probabilita' delle piastre coincide con quella di
  decontam (``batch.combine``) per le tre regole, quando ogni piastra ha il
  suo confronto;
- su un oggetto costruito: i controlli positivi non entrano nel confronto (le
  probabilita' non cambiano qualunque cosa contengano), un reagente presente
  nei negativi e' un contaminante e una variante dei soli biologici no, una
  piastra con meno di ``decontam.min_blanks`` negativi si confronta con i
  negativi di tutte le piastre, la modalita' che decide e' quella dichiarata
  in ``decontam.mode`` e oltre ``qc.max_frac_contaminant`` la fase si ferma
  con ``E-S12-02``;
- S12 sul sottoinsieme di prova: le due modalita', quella dichiarata,
  l'oggetto ripulito con gli identificativi delle altre varianti invariati, il
  tracciamento, gli stessi byte in due esecuzioni;
- l'ordine invertito, il filtro prima della decontaminazione vera, fallisce
  con ``E-S13-01``;
- sul dataset completo: modalita' dichiarata, contaminanti, riscontri e ASV1.

4. Comandi Bash e scenari di esecuzione
---------------------------------------
    1. Modalità locale standard (R di base con jsonlite, senza Bioconductor né
       dati reali):
       pytest tests/test_w20_decontam.py -v

    2. Modalità container Docker standard (sottoinsieme ridotto con
       R/Bioconductor):
       docker run --rm \
         -e PYTHONPATH=/app/src \
         -e AMPLICON16S_R_DIR=/app/R \
         -v "$(pwd)":/app \
         -w /app \
         amplicon16s:dev \
         pytest -o cache_dir=/tmp/.pytest_cache tests/test_w20_decontam.py -v

    3. Modalità container Docker completa (con i dati reali OSD-734):
       docker run --rm \
         -e PYTHONPATH=/app/src \
         -e AMPLICON16S_R_DIR=/app/R \
         -e AMPLICON16S_CONFIG_DATI_REALI="$HOME/ASI/config_osd734.yaml" \
         -v "$(pwd)":/app \
         -v "$HOME/ASI":"$HOME/ASI" \
         -w /app \
         amplicon16s:dev \
         pytest -o cache_dir=/tmp/.pytest_cache tests/test_w20_decontam.py -v

5. Risultato atteso
-------------------
Vedi ``test.txt``, scheda W20.

6. Razionale scientifico e sistemistico
---------------------------------------
- I controlli positivi contengono un organismo aggiunto e contaminanti di
  reagente amplificati: nel confronto falserebbero la prevalenza.
- Il minimo delle probabilita' di molte piastre con pochi negativi e' una
  regola permissiva per costruzione: la modalita' si dichiara in
  configurazione, non si sceglie dopo aver visto quanto rimuove.
- La decontaminazione precede il filtro di prevalenza: i contaminanti da
  reagente sono prevalenti per costruzione.
"""

from __future__ import annotations

import csv
import hashlib
import json
import os
import subprocess
from pathlib import Path
from typing import ClassVar

import pytest
from conftest import Campione, copia_esecuzione, crea_scenario
from sottoinsieme import config_ridotta, motivo_pacchetti_r_assenti

from amplicon16s.config.resolve import risolvi
from amplicon16s.config.schema import ErroreConfigurazione, valida
from amplicon16s.metadata.crosswalk import analizza
from amplicon16s.metadata.models import ClasseCampione
from amplicon16s.errors.exceptions import ErrorePipeline
from amplicon16s.io_layer.artifacts import Fase
from amplicon16s.logging.logger import chiudi
from amplicon16s.rbridge.payload import PREFISSO
from amplicon16s.rbridge.runner import cartella_r, trova_rscript
from amplicon16s.runner.executor import Conclusione, Esecutore
from amplicon16s.runner.graph import Passo
from amplicon16s.runner.project import ProjectRun, StatoPasso, passi_realizzati
from amplicon16s.runner.tracciamento import ricomponi
from amplicon16s.steps.base import PipelineStep, Produzione, StepContext
from amplicon16s.steps.s10_phyloseq import NOME_OGGETTO


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
    "dada2", "ggplot2", "ShortRead", "jsonlite", "phyloseq", "Biostrings", "decontam"
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
    """Esegue codice R con ``R/lib/decontaminazione.R`` caricato; il codice scrive un
    oggetto JSON in ``uscita``.
    """
    uscita = cartella / "uscita.json"
    script = cartella / "prova.R"
    libreria = cartella_r() / "lib" / "decontaminazione.R"
    script.write_text(
        f"source({json.dumps(str(libreria))})\nuscita <- {json.dumps(str(uscita))}\n{codice}\n",
        encoding="utf-8",
    )
    subprocess.run([str(trova_rscript()), "--vanilla", str(script)], check=True,
                   capture_output=True)
    return json.loads(uscita.read_text(encoding="utf-8"))


def _riepilogo(run) -> dict:
    """``decontam_riepilogo.json`` di S12."""
    return json.loads(
        (run.albero.cartella(Fase.CONTROLS) / "decontam_riepilogo.json").read_text()
    )


# --------------------------------------------------------------------------- #
# 0. La riclassificazione dei controlli di campo                               #
# --------------------------------------------------------------------------- #


def _tre_campioni() -> list[Campione]:
    """Due superfici e un tampone mai aperto, tutti dichiarati "Surface swab"."""
    return [
        Campione("ERX4000001", "NOD1D4.L1"),
        Campione("ERX4000002", "NOD1D4.L2"),
        Campione("ERX4000003", "NOD1.N5", posizione="Unopened 3DMM Swab Tube"),
    ]


def test_un_tampone_mai_aperto_diventa_controllo_negativo(tmp_path):
    """
    **Obiettivo**: Verificare che, con la regola predefinita di ``ctrl``
    (colonna ``Factor Value[Sample Location]``, valore ``Unopened 3DMM Swab
    Tube``), il tampone mai aperto sia un controllo negativo con il materiale
    originale ``Surface swab``, che le superfici restino biologiche, che il
    denominatore della prevalenza conti solo queste, e che il confronto non
    distingua maiuscole e spazi ai bordi.

    **Razionale scientifico e sistemistico**: Un tampone mai aperto non ha
    campionato alcuna superficie: tenerlo fra i biologici lo porterebbe nelle
    analisi ecologiche come campione ambientale. Il materiale originale resta
    perche' la riclassificazione sia tracciabile.
    """
    inventario = analizza(crea_scenario(tmp_path, _tre_campioni()).config).inventario()
    classi = {c.nome: (c.classe, c.materiale) for c in inventario}
    assert classi["NOD1.N5"] == (ClasseCampione.CONTROLLO_NEGATIVO, "Surface swab")
    assert classi["NOD1D4.L1"] == (ClasseCampione.BIOLOGICO, "Surface swab")
    assert inventario.denominatore_prevalenza() == 2
    campioni = _tre_campioni()
    campioni[2].posizione = "  unopened 3dmm swab tube "
    inventario = analizza(crea_scenario(tmp_path / "b", campioni).config).inventario()
    assert inventario["ERX4000003"].classe is ClasseCampione.CONTROLLO_NEGATIVO


def test_senza_regola_il_tampone_resta_biologico_e_s0_cambia_impronta(tmp_path):
    """
    **Obiettivo**: Verificare che con l'elenco dei valori vuoto il tampone resti
    biologico, che la regola entri nell'impronta di S0 (cambiarla rifa' la
    catena), e che dichiarare i valori senza la colonna sia respinto.

    **Razionale scientifico e sistemistico**: La classe determina il confronto
    della decontaminazione, il denominatore della prevalenza e la numerazione
    delle varianti: la regola e' un ingresso di S0.
    """
    scenario = crea_scenario(tmp_path, _tre_campioni(), sovrascrivi={"ctrl": {"blank_override_values": []}})
    inventario = analizza(scenario.config).inventario()
    assert inventario["ERX4000003"].classe is ClasseCampione.BIOLOGICO
    s0 = passi_realizzati()[Passo.S0]
    con_regola = risolvi(crea_scenario(tmp_path / "b", _tre_campioni()).config)
    assert s0.calcolata_su(risolvi(scenario.config), {}) != s0.calcolata_su(con_regola, {})
    dati = scenario.config.model_dump(mode="python")
    dati["ctrl"]["blank_override_column"] = None
    dati["ctrl"]["blank_override_values"] = ["Unopened 3DMM Swab Tube"]
    with pytest.raises(ErroreConfigurazione, match="blank_override_column"):
        valida(dati)


# --------------------------------------------------------------------------- #
# 1. Parametri e combinazione delle piastre                                    #
# --------------------------------------------------------------------------- #


def test_i_parametri_di_s12(tmp_path):
    """
    **Obiettivo**: Verificare i parametri di S12: prevalenza, soglia 0,5, cinque
    negativi minimi per piastra, la combinazione delle piastre dichiarata al
    predefinito di decontam (``minimum``), la modalita' aggregata e la frazione
    massima 0,40; e che cambiarli non tocchi l'impronta di S0, S10 e S11.

    **Razionale scientifico e sistemistico**: Un predefinito lasciato implicito
    cambierebbe in silenzio con una versione nuova di decontam.
    """
    config = config_ridotta(tmp_path)
    d = config.decontam
    assert (d.method, d.threshold, d.min_blanks, d.batch_combine, d.mode) == (
        "prevalence", 0.5, 5, "minimum", "aggregate")
    assert config.qc.max_frac_contaminant == 0.40
    prima = risolvi(config)
    dopo = risolvi(config_ridotta(tmp_path, decontam={"batch_combine": "fisher", "mode": "batch"},
                                  qc={"max_frac_contaminant": 0.2}))
    passi = passi_realizzati()
    for passo in (Passo.S0, Passo.S10, Passo.S11):
        assert passi[passo].calcolata_su(prima, {}) == passi[passo].calcolata_su(dopo, {})
    assert passi[Passo.S12].calcolata_su(prima, {}) != passi[Passo.S12].calcolata_su(dopo, {})


def test_la_combinazione_segue_le_regole_dichiarate(r, tmp_path):
    """
    **Obiettivo**: Verificare che ``combina_probabilita`` dia il minimo delle
    probabilita' non mancanti (NA se mancano tutte), il prodotto, e il metodo
    di Fisher con le mancanti a 0,5.

    **Razionale scientifico e sistemistico**: E' la regola di ``batch.combine``
    di decontam 1.28.0, riprodotta per poter dare a una piastra con pochi
    negativi un confronto diverso.
    """
    esito = _r("""
p <- matrix(c(0.2, 0.6, NA, 0.9, NA, NA), nrow = 3, byrow = TRUE,
            dimnames = list(c("A", "B", "C"), c("1", "2")))
jsonlite::write_json(list(
  minimum = as.list(combina_probabilita(p, "minimum")),
  product = as.list(combina_probabilita(p, "product")),
  fisher = as.list(combina_probabilita(p, "fisher"))
), uscita, auto_unbox = TRUE, digits = NA, na = "null")
""", tmp_path)
    import math

    assert esito["minimum"] == {"A": 0.2, "B": 0.9, "C": None}
    assert esito["product"]["A"] == pytest.approx(0.12)
    assert esito["product"]["B"] == pytest.approx(0.9)
    # Fisher con k = 2 probabilita': P(chi2_4 > -2 ln p) = p (1 - ln p).
    for v, attesa in (("A", 0.2 * 0.6), ("B", 0.9 * 0.5), ("C", 0.25)):
        assert esito["fisher"][v] == pytest.approx(attesa * (1 - math.log(attesa)))


def test_la_combinazione_coincide_con_quella_di_decontam(bioc, tmp_path):
    """
    **Obiettivo**: Verificare che, con ogni piastra nel suo confronto, le
    probabilita' per piastra combinate con ``combina_probabilita`` coincidano
    con quelle di ``decontam::isContaminant(batch = , batch.combine = )`` per
    ``minimum``, ``product`` e ``fisher``.

    **Razionale scientifico e sistemistico**: La riproduzione serve solo a
    trattare le piastre con pochi negativi: per le altre deve dare esattamente
    cio' che darebbe decontam.
    """
    esito = _r("""
set.seed(1)
n <- 60; v <- 40
m <- matrix(rpois(n * v, 3) * rbinom(n * v, 1, 0.4), nrow = n,
            dimnames = list(NULL, paste0("V", seq_len(v))))
neg <- rep(c(rep(TRUE, 6), rep(FALSE, 14)), 3)
lotto <- rep(c("1", "2", "3"), each = 20)
uguali <- sapply(c("minimum", "product", "fisher"), function(regola) {
  attesa <- decontam::isContaminant(m, neg = neg, method = "prevalence", batch = lotto,
                                    batch.combine = regola, threshold = 0.5)$p
  p <- sapply(c("1", "2", "3"), function(l) probabilita_prevalenza(m[lotto == l, ], neg[lotto == l]))
  isTRUE(all.equal(unname(combina_probabilita(p, regola)), attesa))
})
jsonlite::write_json(as.list(uguali), uscita, auto_unbox = TRUE)
""", tmp_path)
    assert esito == {"minimum": True, "product": True, "fisher": True}


# --------------------------------------------------------------------------- #
# 2. S12 su un oggetto costruito                                               #
# --------------------------------------------------------------------------- #


#: Tre piastre: le prime due con sei negativi, la terza con due; dieci
#: biologici per piastra, quattro positivi. Varianti: un reagente presente in
#: tutti i negativi e in pochi biologici, una variante dei soli biologici, una
#: variante dei biologici e di qualche negativo, e una dei soli positivi.
#: ``riempi_positivi`` decide che cosa contengono i positivi.
COSTRUISCI = r"""
suppressPackageStartupMessages(library(phyloseq))
piastra <- c(rep("1", 16), rep("2", 16), rep("3", 12), rep("1", 4))
classe <- c(rep(c(rep("controllo_negativo", 6), rep("biologico", 10)), 2),
            rep("controllo_negativo", 2), rep("biologico", 10), rep("controllo_positivo", 4))
campioni <- sprintf("ERX%03d", seq_along(classe))
n <- length(campioni)
neg <- classe == "controllo_negativo"; bio <- classe == "biologico"; pos <- classe == "controllo_positivo"
otu <- rbind(
  Reagente = ifelse(neg, 50, ifelse(bio & seq_len(n) %% 5 == 0, 3, 0)),
  Ambiente = ifelse(bio, 400, 0),
  Comune = ifelse(bio, 200, ifelse(neg & seq_len(n) %% 3 == 0, 5, 0)),
  Positivo = ifelse(pos, 900, 0)
)
if (riempi_positivi) otu[c("Reagente", "Ambiente", "Comune"), pos] <- c(500, 500, 500)
colnames(otu) <- campioni
storage.mode(otu) <- "integer"
tax <- matrix(c("Bacteria", "Bacteria", "Bacteria", "Bacteria", "Reag", "Amb", "Com", "Pos"),
              ncol = 2, dimnames = list(rownames(otu), c("Kingdom", "Genus")))
dati <- data.frame(accession = campioni, classe = classe, piastra = piastra,
                   row.names = campioni, stringsAsFactors = FALSE)
saveRDS(phyloseq(otu_table(otu, taxa_are_rows = TRUE), tax_table(tax), sample_data(dati)), oggetto)
"""


def _s12_su_oggetto_costruito(run, cartella: Path, *, riempi_positivi: bool = False):
    """Sostituisce l'oggetto di S10 con quello costruito ed esegue il calcolo di S12."""
    oggetto = run.albero.cartella(Fase.PHYLOSEQ) / NOME_OGGETTO
    script = cartella / "costruisci.R"
    script.write_text(
        f"oggetto <- {json.dumps(str(oggetto))}\n"
        f"riempi_positivi <- {'TRUE' if riempi_positivi else 'FALSE'}\n" + COSTRUISCI,
        encoding="utf-8",
    )
    subprocess.run([str(trova_rscript()), "--vanilla", str(script)], check=True,
                   capture_output=True)
    fase = run.fase(Passo.S12)
    contesto = run.contesto(Passo.S12, run.valuta())
    return fase.calcola(contesto.ristretto(fase.parametri))


def test_s12_su_un_oggetto_costruito(bioc, oggetto_calcolato, tmp_path):
    """
    **Obiettivo**: Verificare su un oggetto costruito che il reagente presente in
    tutti i negativi sia un contaminante e la variante dei soli biologici no,
    che i quattro positivi siano esclusi dal confronto, e che la piastra 3, con
    due negativi, si confronti con i negativi di tutte le piastre (14).

    **Razionale scientifico e sistemistico**: Sono le regole dichiarate del
    confronto, verificate dove l'esito e' noto in anticipo.
    """
    run = copia_esecuzione(oggetto_calcolato, tmp_path / "copia")
    _s12_su_oggetto_costruito(run, tmp_path)
    varianti = {v["asv_id"]: v for v in _tsv(run.albero.cartella(Fase.CONTROLS) / "decontam_varianti.tsv")}
    assert varianti["Reagente"]["contaminante_aggregata"] == "si"
    assert varianti["Ambiente"]["contaminante_aggregata"] == "no"
    assert varianti["Ambiente"]["contaminante_per_piastra"] == "no"
    riepilogo = _riepilogo(run)
    assert riepilogo["confronto"] == {
        "negativi": 14, "biologici": 30, "positivi_esclusi": 4, "senza_letture_esclusi": 0}
    assert riepilogo["piastre"]["1"]["confronto"] == "negativi della piastra"
    assert riepilogo["piastre"]["3"]["negativi"] == 14
    assert "negativi di tutte le piastre" in riepilogo["piastre"]["3"]["confronto"]
    assert riepilogo["modalita_dichiarata"] == "aggregate"
    assert varianti["Reagente"]["rimossa"] == "si"


def test_la_modalita_dichiarata_decide_i_contaminanti_rimossi(bioc, oggetto_calcolato, tmp_path):
    """
    **Obiettivo**: Verificare che con ``decontam.mode: batch`` le varianti
    rimosse siano i contaminanti per piastra e con ``aggregate`` quelli
    aggregati, e che entrambe le modalita' siano comunque registrate.

    **Razionale scientifico e sistemistico**: La modalita' e' una scelta di
    metodo dichiarata, non un esito: l'altra resta come diagnostica.
    """
    for modo, colonna in (("aggregate", "contaminante_aggregata"), ("batch", "contaminante_per_piastra")):
        run = copia_esecuzione(oggetto_calcolato, tmp_path / modo, decontam={"mode": modo})
        _s12_su_oggetto_costruito(run, tmp_path)
        righe = _tsv(run.albero.cartella(Fase.CONTROLS) / "decontam_varianti.tsv")
        assert [r["rimossa"] for r in righe] == [r[colonna] for r in righe], modo
        assert set(_riepilogo(run)["modalita"]) == {"aggregata", "per_piastra"}


def test_i_controlli_positivi_non_entrano_nel_confronto(bioc, oggetto_calcolato, tmp_path):
    """
    **Obiettivo**: Verificare che le probabilita' di decontam, aggregate e per
    piastra, siano identiche quando i positivi sono vuoti di reagente e quando
    lo contengono insieme alle altre varianti.

    **Razionale scientifico e sistemistico**: Se i positivi entrassero fra i
    campioni, la loro composizione cambierebbe la prevalenza di ogni variante
    che contengono.
    """
    colonne = ("asv_id", "p_aggregata", "p_per_piastra")
    esiti = []
    for riempi in (False, True):
        run = copia_esecuzione(oggetto_calcolato, tmp_path / f"copia_{riempi}")
        _s12_su_oggetto_costruito(run, tmp_path, riempi_positivi=riempi)
        righe = _tsv(run.albero.cartella(Fase.CONTROLS) / "decontam_varianti.tsv")
        esiti.append([tuple(r[c] for c in colonne) for r in righe])
    assert esiti[0] == esiti[1]


@pytest.mark.parametrize("modo", ["aggregate", "batch"])
def test_oltre_la_frazione_massima_la_modalita_dichiarata_si_ferma(
    bioc, oggetto_calcolato, tmp_path, modo
):
    """
    **Obiettivo**: Verificare che, con ``qc.max_frac_contaminant`` a zero e un
    reagente da rimuovere che compare anche in alcuni biologici, S12 si fermi
    con ``E-S12-02``, di revisione umana, qualunque sia la modalita' dichiarata,
    lasciando le misure di entrambe le modalita'.

    **Razionale scientifico e sistemistico**: Una decontaminazione che toglie
    quasi tutto il segnale dei campioni indica piu' spesso un confronto
    sbagliato che una contaminazione reale: la decisione non e' automatica.
    """
    run = copia_esecuzione(oggetto_calcolato, tmp_path / "copia",
                           qc={"max_frac_contaminant": 0.0}, decontam={"mode": modo})
    with pytest.raises(ErrorePipeline) as info:
        _s12_su_oggetto_costruito(run, tmp_path)
    assert info.value.codice == "E-S12-02"
    assert info.value.categoria.value == "revisione_umana"
    riepilogo = _riepilogo(run)
    assert riepilogo["modalita_dichiarata"] == modo
    assert riepilogo["entro_max_frazione"] is False
    assert riepilogo["modalita"]["aggregata"]["letture_rimosse"]["biologico"] > 0
    assert riepilogo["modalita"]["per_piastra"]["letture_rimosse"]["biologico"] > 0


# --------------------------------------------------------------------------- #
# 3. S12 sul sottoinsieme di prova                                             #
# --------------------------------------------------------------------------- #


def test_s12_confronta_le_modalita_e_rimuove_senza_rinumerare(bioc, decontam_calcolata, tmp_path):
    """
    **Obiettivo**: Verificare che S12 sulla versione ridotta si concluda dopo
    S11, che il riepilogo registri le due modalita' con la frazione rimossa per
    classe e la modalita' dichiarata, che l'oggetto ripulito abbia tutti i campioni e
    le varianti non rimosse con gli stessi identificativi, nello stesso ordine,
    e che l'elenco dei contaminanti sia quello delle varianti rimosse.

    **Razionale scientifico e sistemistico**: Gli identificativi si assegnano in
    S10 e non si rinumerano: dopo la rimozione la numerazione ha dei vuoti.
    """
    run, esito = decontam_calcolata
    assert esito.conclusione is Conclusione.COMPLETATA
    assert [r.passo for r in esito.eseguite] == [Passo.S12]
    riepilogo = _riepilogo(run)
    assert set(riepilogo["modalita"]) == {"aggregata", "per_piastra"}
    for m in riepilogo["modalita"].values():
        assert set(m["letture_rimosse"]) == {"biologico", "controllo_negativo", "controllo_positivo"}
    assert riepilogo["modalita_dichiarata"] == "aggregate"
    assert riepilogo["entro_max_frazione"] is True and riepilogo["esito"].startswith("aggregata")
    assert riepilogo["contaminanti_rimossi"] == riepilogo["modalita"]["aggregata"]["contaminanti"]
    assert riepilogo["confronto"]["positivi_esclusi"] == 9
    cartella = run.albero.cartella(Fase.CONTROLS)
    varianti = _tsv(cartella / "decontam_varianti.tsv")
    rimosse = [v["asv_id"] for v in varianti if v["rimossa"] == "si"]
    assert [c["asv_id"] for c in _tsv(cartella / "decontam_contaminanti.tsv")] == rimosse
    assert riepilogo["contaminanti_rimossi"] == len(rimosse)
    uscita = tmp_path / "uscita.json"
    script = tmp_path / "leggi.R"
    script.write_text(
        "suppressPackageStartupMessages(library(phyloseq))\n"
        f"a <- readRDS({json.dumps(str(run.albero.cartella(Fase.PHYLOSEQ) / NOME_OGGETTO))})\n"
        f"b <- readRDS({json.dumps(str(cartella / 'ps_decontaminato.rds'))})\n"
        "jsonlite::write_json(list(prima = taxa_names(a), dopo = taxa_names(b), "
        "campioni_prima = sample_names(a), campioni_dopo = sample_names(b)), "
        f"{json.dumps(str(uscita))})\n",
        encoding="utf-8",
    )
    subprocess.run([str(trova_rscript()), "--vanilla", str(script)], check=True, capture_output=True)
    nomi = json.loads(uscita.read_text())
    assert nomi["campioni_dopo"] == nomi["campioni_prima"]
    assert nomi["dopo"] == [v for v in nomi["prima"] if v not in set(rimosse)]
    tracciamento = ricomponi(run)
    assert tracciamento.passi[-1] == "decontaminate"
    assert tracciamento.origine["decontaminate"] == "S12"
    for campione, letture in tracciamento.letture.items():
        assert letture["decontaminate"] <= letture["lunghezza"], campione
    print(f"\nS12 sulla versione ridotta: {dict(esito.eseguite[0].metriche)}")


def test_s12_da_gli_stessi_byte_in_due_esecuzioni(bioc, decontam_calcolata, tmp_path):
    """
    **Obiettivo**: Verificare che S12, rieseguita sugli stessi artefatti, dia
    gli stessi byte nei propri artefatti di ``11_controls``, oggetto ripulito
    compreso, senza toccare quelli di S11.

    **Razionale scientifico e sistemistico**: decontam per prevalenza non usa
    numeri casuali; la cartella e' condivisa, e ogni fase riscrive solo i suoi
    file.
    """
    run, _ = decontam_calcolata

    def impronte(cartella: Path) -> dict[str, str]:
        return {
            p.name: hashlib.md5(p.read_bytes()).hexdigest()
            for p in sorted(cartella.iterdir())
            if not p.name.startswith((PREFISSO, "manifest"))
        }

    attese = impronte(run.albero.cartella(Fase.CONTROLS))
    assert "ps_decontaminato.rds" in attese and "soglia.json" in attese
    copia = copia_esecuzione(decontam_calcolata, tmp_path)
    copia.albero.rimuovi_manifesto_passo(Passo.S12, Fase.CONTROLS)
    esito = Esecutore(copia, fino_a=Passo.S12).esegui()
    assert [r.passo for r in esito.eseguite] == [Passo.S12]
    assert impronte(copia.albero.cartella(Fase.CONTROLS)) == attese
    assert copia.valuta().situazioni[Passo.S11].stato is StatoPasso.COMPLETATA


class _FiltroDoppione(PipelineStep):
    """S13 finta: il filtro di prevalenza, che non deve girare prima di S12."""

    passo: ClassVar[Passo] = Passo.S13
    parametri: ClassVar[tuple[str, ...]] = ("prev",)

    def calcola(self, contesto: StepContext) -> Produzione:
        """Scrive un artefatto vuoto."""
        return Produzione((contesto.albero.scrivi_testo(self.cartella, "filtro.txt", "x"),))


def test_il_filtro_prima_della_decontaminazione_fallisce_con_e_s13_01(
    bioc, controlli_calcolati, tmp_path
):
    """
    **Obiettivo**: Verificare che, con S0-S11 concluse e la S12 vera non ancora
    eseguita, avviare S13 sollevi ``E-S13-01``; e che dopo la S12 vera S13 si
    concluda.

    **Razionale scientifico e sistemistico**: I contaminanti da reagente sono
    prevalenti per costruzione: filtrare prima li renderebbe indistinguibili
    dal segnale. Il vincolo e' nel grafo, e vale con la fase vera.
    """
    copia = copia_esecuzione(controlli_calcolati, tmp_path)
    passi = {**passi_realizzati(), Passo.S13: _FiltroDoppione()}
    run = ProjectRun(copia.config, passi=passi)
    assert run.valuta().situazioni[Passo.S12].stato is StatoPasso.DA_ESEGUIRE
    with pytest.raises(ErrorePipeline) as info:
        run.fase(Passo.S13).esegui(run.contesto(Passo.S13))
    assert info.value.codice == "E-S13-01"
    esito = Esecutore(run, fino_a=Passo.S13).esegui()
    assert [r.passo for r in esito.eseguite] == [Passo.S12, Passo.S13]


# --------------------------------------------------------------------------- #
# Dataset completo, nel container                                              #
# --------------------------------------------------------------------------- #


@pytest.mark.dati_reali
def test_s12_sul_dataset_completo(bioc, catena_reale):
    """
    **Obiettivo**: Verificare che sul dataset completo S12 si concluda dopo
    S11, e riportarne la modalita' dichiarata, le due modalita' a confronto, i
    contaminanti principali e lo stato di ASV1, che non deve essere rimossa.

    **Razionale scientifico e sistemistico**: ASV1, un Pseudomonas, raccoglie
    oltre 4 milioni di letture dei biologici: la sua rimozione toglierebbe una
    frazione rilevante del dataset.
    """
    run, esito = catena_reale
    assert esito.conclusione is Conclusione.COMPLETATA
    s12 = next(r for r in esito.eseguite if r.passo is Passo.S12)
    print(f"\nS12: {s12.secondi} s, {dict(s12.metriche)}")
    riepilogo = _riepilogo(run)
    print(json.dumps(riepilogo["modalita"], indent=1))
    varianti = {v["asv_id"]: v for v in _tsv(run.albero.cartella(Fase.CONTROLS) / "decontam_varianti.tsv")}
    print(f"ASV1: {varianti['ASV1']}")
    assert varianti["ASV1"]["rimossa"] == "no"
