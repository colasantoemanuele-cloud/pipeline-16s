"""Suite di verifica dei 15 Gate di validazione pre-analitica (G01–G15) e della Fase S0.

Inquadramento nel Piano Operativo
---------------------------------
* **Settimana di riferimento**: **Settimana 7 (W7 — Fase F2: Validazione degli
  ingressi, i 15 Gate G01–G15 e Fase S0)**.
* **Moduli sorgente coperti**:
  - ``src/amplicon16s/gates/g01_g15.py``
  - ``src/amplicon16s/gates/registry.py``
  - ``src/amplicon16s/steps/s00_validate.py``

Scopo sperimentale e razionale scientifico/sistemistico
-------------------------------------------------------
Ogni gate possiede un caso nominale che lo supera e almeno un caso avverso che
lo fa fallire con il codice di catalogo atteso (``E-G15-*`` o ``E-S0-*``). Gli
scenari scrivono su disco veri archivi FASTQ compressi ``gzip``, tabelle ISA-Tab
di metadati e file FASTA di riferimento tassonomico, collaudando l'intero
percorso I/O reale della Fase S0:

1. **G15 valutato per primo a costo zero**: la coerenza matematica della
   configurazione YAML viene verificata prima di aprire qualsiasi file su disco.
2. **Esecuzione Fail-Fast**: il registro interrompe la valutazione al primo gate
   fallito, impedendo che un errore strutturale a monte (es. file mancante o
   tabella corrotta) generi una cascata di falsi errori derivati nei gate successivi.
3. **G07 (Layout Single-End)**: rifiuta file di lettura inversa (``_R2``) poiché
   l'intera pipeline per OSD-734 è calibrata rigorosamente su corse single-end.
4. **G08 (Lotti e controlli negativi per piastra)**: verifica la colonna di
   lotto e, qualora una piastra contenga meno di ``decontam.min_blanks`` controlli
   negativi, emette l'avviso di degradazione ``E-S0-15`` senza bloccare la corsa.
5. **G09 (Compatibilità di ``truncLen``)**: garantisce che ``filter.truncLen``
   non superi la lunghezza minima osservata nei FASTQ, poiché ``filterAndTrim``
   di DADA2 **scarta** (e non accorcia) le letture più corte di ``truncLen``.
6. **G10 (Assenza del primer 515F e controllo positivo del motivo V4)**:
   accerta che il primer non sia ancora attaccato in 5' e che la regione
   conservata V4 sia presente; esclude dal requisito i controlli negativi (da
   cui non si pretende amplificazione) e calcola la **mediana** sui campioni
   biologici per tollerare singoli campioni a bassa resa.
7. **G12 (Integrità MD5 del riferimento SILVA)**: verifica il checksum MD5 del
   database FASTA prima dell'assegnazione tassonomica in S8.
8. **G13 (Grammatica FASTQ e integrità Gzip)**: convalida la decompressione
   ``gzip`` e la struttura a 4 righe per record con pari lunghezza sequenza/Phred.
9. **G14 (Risorse di calcolo e spazio disco)**: confronta ``run.threads`` con le
   CPU effettivamente assegnate al processo (``os.sched_getaffinity``) e verifica
   lo spazio libero su disco rispetto alla mole dei FASTQ in ingresso.
10. **Fase S0 completa**: certifica la scrittura atomica dei 4 artefatti in
    ``01_input_validation/`` (``gates.json``, ``crosswalk.tsv``, ``inventario.json``,
    ``letture_ispezionate.tsv``) e il calcolo del denominatore di prevalenza
    sui soli campioni biologici (803 nel dataset reale OSD-734).
"""

from __future__ import annotations

import gzip
import json

import pytest
from conftest import (
    INIZIO_CON_PRIMER,
    INIZIO_MUTO,
    NEGATIVO,
    POSITIVO,
    Campione,
    crea_scenario,
    lettura,
)

from amplicon16s.gates.g01_g15 import Contesto, ErroreGate
from amplicon16s.gates.registry import REGISTRO, esegui_tutti, nomi_dei_gate
from amplicon16s.metadata.models import ClasseCampione
from amplicon16s.steps.s00_validate import esegui_s0


def _campioni() -> list[Campione]:
    """Quattro campioni: due biologici, un positivo, un negativo."""
    return [
        Campione("ERX3000001", "NOD1D4.L1", piastra="1"),
        Campione("ERX3000002", "NOD1D4.L2", piastra="1"),
        Campione("ERX3000003", "POS.P1.1", materiale=POSITIVO,
                 posizione="Not Applicable", piastra="1"),
        Campione("ERX3000004", "BLANK.P1.1", materiale=NEGATIVO,
                 posizione="Not Applicable", piastra="1"),
    ]


def _scenario(tmp_path, campioni=None, **extra):
    opzioni = {"con_letture": True}
    opzioni.update(extra)
    return crea_scenario(tmp_path, campioni or _campioni(), **opzioni)


def _esegui(scenario) -> dict:
    """Esegue tutti i gate e restituisce gli esiti indicizzati per nome."""
    return {e.gate: e for e in esegui_tutti(Contesto(scenario.config))}


def _fallito(esiti: dict) -> str | None:
    for nome, esito in esiti.items():
        if esito.eseguito and not esito.superato:
            return nome
    return None


# --------------------------------------------------------------------------- #
# Il registro                                                                  #
# --------------------------------------------------------------------------- #


def test_i_gate_sono_quindici():
    """
    **Obiettivo**: Verificare che il registro ``REGISTRO`` contenga esattamente
    15 gate con nomi univoci (``G01``–``G15``).

    **Razionale Scientifico/Sistemistico**: Assicura che nessuno dei 15 gate
    previsti dal protocollo di validazione pre-analitica della Fase S0 sia stato
    omesso o duplicato nel registro di esecuzione.
    """
    assert len(REGISTRO) == 15
    assert len(set(nomi_dei_gate())) == 15


def test_la_configurazione_si_controlla_per_prima():
    """
    **Obiettivo**: Verificare che il primo gate nella sequenza di esecuzione
    ``nomi_dei_gate()`` sia ``G15``.

    **Razionale Scientifico/Sistemistico**: Il Gate G15 valida la coerenza
    matematica dei parametri in RAM a costo zero (senza toccare il disco).
    Valutarlo per primo evita di ispezionare inutilmente centinaia di file FASTQ
    quando la configurazione contiene già una contraddizione logica.
    """
    assert nomi_dei_gate()[0] == "G15"


def test_tutti_i_gate_passano_su_uno_scenario_corretto(tmp_path):
    """
    **Obiettivo**: Verificare che su uno scenario completo e conforme tutti i 15
    gate risultino eseguiti (``eseguito is True``) e superati (``superato is True``).

    **Razionale Scientifico/Sistemistico**: Certifica l'assenza di falsi allarmi
    quando i file di ingresso, le tabelle ISA-Tab, le letture FASTQ e il
    riferimento SILVA soddisfano i requisiti dello studio.
    """
    esiti = _esegui(_scenario(tmp_path))
    assert all(e.eseguito and e.superato for e in esiti.values()), _fallito(esiti)


def test_l_esecuzione_si_ferma_al_primo_gate_fallito(tmp_path):
    """
    **Obiettivo**: Verificare la politica *fail-fast* di ``esegui_tutti``: se
    ``G04`` fallisce, i gate successivi (es. ``G11``) non vengono eseguiti e
    risultano ``superato is False``.

    **Razionale Scientifico/Sistemistico**: Proseguire la validazione su tabelle
    o elenchi di file già incoerenti produrrebbe una cascata di errori derivati
    fuorvianti, rendendo difficile per il biologo computazionale individuare la
    causa radice del problema.
    """
    campioni = _campioni()
    campioni[0].materiale = "solvent control"   # fa fallire G11
    campioni[1].file = "senza_codice.fastq.gz"  # farebbe fallire G04, prima
    esiti = _esegui(_scenario(tmp_path, campioni))

    assert esiti["G04"].eseguito and not esiti["G04"].superato
    assert not esiti["G11"].eseguito
    assert not esiti["G11"].superato, "un gate non eseguito non e' superato"


# --------------------------------------------------------------------------- #
# G01 — gli ingressi esistono e sono leggibili                                 #
# --------------------------------------------------------------------------- #


def test_g01_passa(tmp_path):
    """
    **Obiettivo**: Verificare che il Gate G01 sia superato quando cartella FASTQ,
    Assay Table, Study Table e riferimento esistono e sono leggibili.

    **Razionale Scientifico/Sistemistico**: Conferma il controllo preliminare di
    esistenza fisica dei percorsi di input prima della lettura delle tabelle.
    """
    assert _esegui(_scenario(tmp_path))["G01"].superato


@pytest.mark.parametrize("quale", ["assay.txt", "studio.txt"])
def test_g01_fallisce_su_un_ingresso_mancante(quale, tmp_path):
    """
    **Obiettivo**: Verificare che la rimozione di ``assay.txt`` o ``studio.txt``
    faccia fallire G01 con codice ``E-S0-01`` indicando il nome del file assente.

    **Razionale Scientifico/Sistemistico**: Intercetta immediatamente percorsi
    errati nel file YAML o file ISA-Tab non scaricati, fornendo una diagnosi
    puntuale sul file mancante.
    """
    scenario = _scenario(tmp_path)
    (scenario.radice / quale).unlink()

    esito = _esegui(scenario)["G01"]
    assert not esito.superato
    assert {v.codice for v in esito.violazioni} == {"E-S0-01"}
    assert quale in str(esito.violazioni[0])


def test_g01_fallisce_se_non_ci_sono_file_di_letture(tmp_path):
    """
    **Obiettivo**: Verificare che G01 fallisca con ``E-S0-01`` se la directory
    ``io.fastq_dir`` esiste ma non contiene alcun file ``*.fastq.gz``.

    **Razionale Scientifico/Sistemistico**: Impedisce che una directory montata
    vuota o un pattern ``io.fastq_glob`` errato lascino procedere la pipeline
    con zero campioni.
    """
    scenario = _scenario(tmp_path)
    for percorso in (scenario.radice / "fastq").glob("*.fastq.gz"):
        percorso.unlink()

    esito = _esegui(scenario)["G01"]
    assert not esito.superato
    assert "non contiene alcun file" in str(esito.violazioni[0])


# --------------------------------------------------------------------------- #
# G02 — le tabelle si aprono e hanno le colonne attese                         #
# --------------------------------------------------------------------------- #


def test_g02_passa(tmp_path):
    """
    **Obiettivo**: Verificare che il Gate G02 sia superato quando Assay Table e
    Study Table espongono tutte le colonne dichiarate nella sezione ``meta`` e ``ctrl``.

    **Razionale Scientifico/Sistemistico**: Accerta la conformità dello schema
    tabellare ISA-Tab prima di avviare il join del crosswalk.
    """
    assert _esegui(_scenario(tmp_path))["G02"].superato


def test_g02_fallisce_su_una_colonna_assente(tmp_path):
    """
    **Obiettivo**: Verificare che la rinomina della colonna ``Characteristics[Material Type]``
    in ``studio.txt`` faccia fallire G02 con codice ``E-S0-02`` citando ``ctrl.column``.

    **Razionale Scientifico/Sistemistico**: Previene errori silenziosi di
    classificazione dei controlli causati da discrepanze nei nomi delle colonne
    tra versioni diverse dei metadati GeneLab.
    """
    scenario = _scenario(tmp_path)
    percorso = scenario.radice / "studio.txt"
    righe = percorso.read_text(encoding="utf-8").splitlines()
    righe[0] = righe[0].replace("Characteristics[Material Type]", "Altro Nome")
    percorso.write_text("\n".join(righe) + "\n", encoding="utf-8")

    esito = _esegui(scenario)["G02"]
    assert not esito.superato
    assert {v.codice for v in esito.violazioni} == {"E-S0-02"}
    assert "ctrl.column" in str(esito.violazioni[0])


# --------------------------------------------------------------------------- #
# G07 — layout single-end                                                      #
# --------------------------------------------------------------------------- #


def test_g07_passa(tmp_path):
    """
    **Obiettivo**: Verificare che il Gate G07 sia superato quando tutti i file
    FASTQ seguono la nomenclatura single-end.

    **Razionale Scientifico/Sistemistico**: Conferma l'assenza di file mate ``_R2``
    nello scenario standard single-end.
    """
    assert _esegui(_scenario(tmp_path))["G07"].superato


def test_g07_fallisce_su_un_file_di_lettura_inversa(tmp_path):
    """
    **Obiettivo**: Verificare che la presenza di un file con suffisso ``_R2.fastq.gz``
    faccia fallire il Gate G07 con codice ``E-S0-07``.

    **Razionale Scientifico/Sistemistico**: La pipeline è strettamente calibrata
    su letture single-end (senza step ``mergePairs`` di DADA2). Se venissero
    introdotti file reverse ``_R2``, questi verrebbero erroneamente trattati come
    campioni indipendenti raddoppiando artificialmente il dataset o perdendo
    l'accoppiamento paired-end.
    """
    campioni = _campioni()
    campioni[0].file = "ERX3000001_NOD1D4.L1_R2.fastq.gz"
    esito = _esegui(_scenario(tmp_path, campioni))["G07"]

    assert not esito.superato
    assert {v.codice for v in esito.violazioni} == {"E-S0-07"}
    assert "letture inverse" in str(esito.violazioni[0])


# --------------------------------------------------------------------------- #
# G08 — lotto coerente e piastre plausibili                                    #
# --------------------------------------------------------------------------- #


def test_g08_passa_senza_arricchimento(tmp_path):
    """
    **Obiettivo**: Verificare che il Gate G08 passi senza errori quando la
    tabella di arricchimento opzionale dei lotti non è attiva.

    **Razionale Scientifico/Sistemistico**: Garantisce che la validazione dei
    lotti funzioni sia con metadati nativi sia con tabella di lotto esterna.
    """
    assert _esegui(_scenario(tmp_path))["G08"].superato


def test_g08_fallisce_se_manca_la_colonna_del_lotto(tmp_path):
    """
    **Obiettivo**: Verificare che l'assenza della colonna ``extraction_plate_num``
    (``decontam.batch_column``) nel file ``lotti.tsv`` faccia fallire G08 con ``E-S0-08``.

    **Razionale Scientifico/Sistemistico**: ``decontam`` in S12 richiede la
    stratificazione per piastra di estrazione; se la colonna indicata da
    ``decontam.batch_column`` non esiste, la decontaminazione per lotto è impossibile.
    """
    scenario = _scenario(tmp_path, con_arricchimento=True)
    percorso = scenario.radice / "lotti.tsv"
    righe = percorso.read_text(encoding="utf-8").splitlines()
    righe[0] = righe[0].replace("extraction_plate_num", "altro")
    percorso.write_text("\n".join(righe) + "\n", encoding="utf-8")

    esito = _esegui(scenario)["G08"]
    assert not esito.superato
    assert {v.codice for v in esito.violazioni} == {"E-S0-08"}
    assert "decontam.batch_column" in str(esito.violazioni[0])


def test_g08_avvisa_su_una_piastra_con_pochi_controlli_negativi(tmp_path):
    """
    **Obiettivo**: Verificare che quando una piastra ha meno controlli negativi
    di ``decontam.min_blanks`` (1 invece di 3), G08 risulti comunque superato
    (``superato is True``) ma registri l'avviso di degradazione ``E-S0-15``.

    **Razionale Scientifico/Sistemistico**: Una piastra con pochi blank (come
    ``Plate 4`` in OSD-734 quando si alza ``min_blanks``) non è un errore fatale
    che deve abortire l'intero studio: G08 emette l'avviso ``E-S0-15`` per
    attivare in S12 la degradazione controllata verso la stima globale di
    ``decontam`` su quella piastra.
    """
    scenario = _scenario(
        tmp_path, con_arricchimento=True,
        sovrascrivi={"decontam": {"min_blanks": 3}},
    )
    esito = _esegui(scenario)["G08"]

    assert esito.superato, "un avviso non blocca"
    assert len(esito.avvisi) == 1
    assert esito.avvisi[0].codice == "E-S0-15"
    assert "1 controlli negativi" in str(esito.avvisi[0])
    assert "decontam.min_blanks" in str(esito.avvisi[0])


def test_g08_non_avvisa_quando_i_negativi_bastano(tmp_path):
    """
    **Obiettivo**: Verificare che quando il numero di controlli negativi sulla
    piastra eguaglia ``decontam.min_blanks`` (4 su 4), G08 passi con zero avvisi.

    **Razionale Scientifico/Sistemistico**: Conferma che la soglia ``min_blanks``
    è valutata in modo assoluto per-piastra (criterio metodologico di potenza
    statistica per ``decontam``) e non per confronto relativo tra piastre.
    """
    campioni = _campioni()
    for campione in campioni:
        campione.materiale = NEGATIVO
        campione.posizione = "Not Applicable"
    scenario = _scenario(
        tmp_path, campioni, con_arricchimento=True,
        sovrascrivi={"decontam": {"min_blanks": 4}},
    )
    esito = _esegui(scenario)["G08"]

    assert esito.superato
    assert esito.avvisi == ()


# --------------------------------------------------------------------------- #
# G09 — troncamento compatibile con le lunghezze                               #
# --------------------------------------------------------------------------- #


def test_g09_passa_con_il_troncamento_predefinito(tmp_path):
    """
    **Obiettivo**: Verificare che il Gate G09 sia superato con ``truncLen = 137``
    su letture lunghe 151 bp.

    **Razionale Scientifico/Sistemistico**: Con letture di 151 bp e troncamento
    a 137 bp, tutte le read possiedono lunghezza sufficiente per essere troncate
    senza essere scartate da DADA2.
    """
    assert _esegui(_scenario(tmp_path))["G09"].superato


def test_g09_fallisce_se_il_troncamento_supera_le_letture(tmp_path):
    """
    **Obiettivo**: Verificare che impostare ``filter.truncLen = 152`` su letture
    di 151 bp faccia fallire G09 con codice ``E-S0-09``.

    **Razionale Scientifico/Sistemistico**: In DADA2 ``filterAndTrim``, le
    letture più corte di ``truncLen`` vengono **scartate integralmente, non
    accorciate**. Un ``truncLen`` superiore alla lunghezza osservata causerebbe
    la perdita del 100% delle letture nella Fase S2.
    """
    scenario = _scenario(tmp_path, sovrascrivi={"filter": {"truncLen": 152}})
    esito = _esegui(scenario)["G09"]

    assert not esito.superato
    assert {v.codice for v in esito.violazioni} == {"E-S0-09"}
    assert "scartate, non accorciate" in str(esito.violazioni[0])


def test_g09_fallisce_su_un_singolo_campione_piu_corto(tmp_path):
    """
    **Obiettivo**: Verificare che G09 fallisca se anche un solo campione ha
    letture di 120 bp inferiori a ``truncLen = 137`` bp.

    **Razionale Scientifico/Sistemistico**: Impedisce che un singolo campione
    con ciclo di sequenziamento più corto venga silenziosamente azzerato da
    ``filterAndTrim`` in S2.
    """
    campioni = _campioni()
    campioni[2].lunghezza_letture = 120
    esito = _esegui(_scenario(tmp_path, campioni))["G09"]

    assert not esito.superato
    assert "120 bp" in str(esito.violazioni[0])


def test_g09_avvisa_quando_si_tronca_molto_sotto_il_minimo(tmp_path):
    """
    **Obiettivo**: Verificare che impostare ``filter.truncLen = 100`` su letture
    di 151 bp (scarto di 51 bp) superi G09 ma emetta l'avviso ``E-S1-01``.

    **Razionale Scientifico/Sistemistico**: Troncare eccessivamente non fa
    perdere campioni ma sacrifica 51 basi utili della regione ipervariabile V4,
    riducendo la risoluzione tassonomica delle ASV fino al livello di genere/specie.
    """
    scenario = _scenario(tmp_path, sovrascrivi={"filter": {"truncLen": 100}})
    esito = _esegui(scenario)["G09"]

    assert esito.superato
    assert len(esito.avvisi) == 1
    assert esito.avvisi[0].codice == "E-S1-01"
    assert "51 bp" in str(esito.avvisi[0])


# --------------------------------------------------------------------------- #
# G10 — primer assente, con controllo positivo                                 #
# --------------------------------------------------------------------------- #


def test_g10_passa(tmp_path):
    """
    **Obiettivo**: Verificare che G10 sia superato quando le letture non hanno
    il primer 515F in testa e contengono il motivo conservato V4 atteso.

    **Razionale Scientifico/Sistemistico**: Valida simultaneamente che i FASTQ
    siano già stati demultiplexati/trimmati dal primer e che contengano vero
    amplicone 16S V4.
    """
    assert _esegui(_scenario(tmp_path))["G10"].superato


def test_g10_fallisce_se_il_primer_e_in_testa(tmp_path):
    """
    **Obiettivo**: Verificare che la presenza della sequenza del primer forward
    515F in 5' (``INIZIO_CON_PRIMER``) faccia fallire G10 con codice ``E-S0-10``
    suggerendo ``filter.trimLeft``.

    **Razionale Scientifico/Sistemistico**: Se il primer PCR degenerato (con basi
    IUPAC ambigue) rimane attaccato in 5' con ``trimLeft = 0``, DADA2 interpreta
    le degenerazioni del primer come vere variazioni biologiche, generando
    decine di ASV chimeriche o spurie per ogni specie reale.
    """
    campioni = _campioni()
    for campione in campioni:
        campione.inizio_letture = INIZIO_CON_PRIMER
    esito = _esegui(_scenario(tmp_path, campioni))["G10"]

    assert not esito.superato
    codici = {v.codice for v in esito.violazioni}
    assert codici == {"E-S0-10"}
    assert "filter.trimLeft" in str(esito.violazioni[0])


def test_g10_fallisce_se_manca_il_segnale_atteso(tmp_path):
    """
    **Obiettivo**: Verificare che G10 fallisca con ``E-S0-10`` se le letture non
    hanno il primer ma sono prive anche del motivo conservato 16S V4 (``INIZIO_MUTO``).

    **Razionale Scientifico/Sistemistico**: Costituisce il *controllo positivo*
    del Gate G10: la sola assenza del primer 515F non basta (potrebbe trattarsi
    di sequenze casuali, adattatori o di un'altra regione genica come ITS/18S).
    La verifica del motivo conservato V4 garantisce che le letture appartengano
    davvero al gene 16S rRNA.
    """
    campioni = _campioni()
    for campione in campioni:
        campione.inizio_letture = INIZIO_MUTO
    esito = _esegui(_scenario(tmp_path, campioni))["G10"]

    assert not esito.superato
    assert "motivo conservato" in str(esito.violazioni[0])
    assert "qc.min_motif_frac" in str(esito.violazioni[0])


def test_g10_non_pretende_il_segnale_dai_controlli_negativi(tmp_path):
    """
    **Obiettivo**: Verificare che G10 sia superato anche quando i controlli
    negativi (``NEGATIVO``) hanno letture prive del motivo conservato V4 (``INIZIO_MUTO``).

    **Razionale Scientifico/Sistemistico**: Da un bianco di estrazione o di
    reagente perfettamente pulito non ci si attende DNA batterico 16S ma solo
    rumore di fondo o dimeri di adattatori; pretendere il motivo V4 nei controlli
    negativi penalizzerebbe paradossalmente i lotti con i blank più puliti.
    """
    campioni = _campioni()
    for campione in campioni:
        if campione.materiale == NEGATIVO:
            campione.inizio_letture = INIZIO_MUTO
    esito = _esegui(_scenario(tmp_path, campioni))["G10"]

    assert esito.superato


def test_g10_tollera_un_singolo_campione_biologico_muto(tmp_path):
    """
    **Obiettivo**: Verificare che G10 sia superato se un singolo campione
    biologico è privo del motivo V4 ma la mediana sui campioni biologici/positivi
    supera ``qc.min_motif_frac``.

    **Razionale Scientifico/Sistemistico**: In studi ad alto campionamento su
    superfici a bassissima biomassa come la ISS (OSD-734), è fisiologico che
    singoli tamponi abbiano resa PCR quasi nulla. L'uso della **mediana**
    anziché del minimo evita che un singolo campione fallito blocchi la
    validazione dell'intera coorte di 960 campioni.
    """
    campioni = _campioni()
    campioni[0].inizio_letture = INIZIO_MUTO
    esito = _esegui(_scenario(tmp_path, campioni))["G10"]

    assert esito.superato


# --------------------------------------------------------------------------- #
# G12 — il riferimento tassonomico è quello dichiarato                         #
# --------------------------------------------------------------------------- #


def test_g12_passa(tmp_path):
    """
    **Obiettivo**: Verificare che il Gate G12 sia superato quando l'impronta MD5
    del file ``riferimento.fa.gz`` coincide con ``tax.ref_md5``.

    **Razionale Scientifico/Sistemistico**: Garantisce l'identità crittografica
    del database tassonomico SILVA prima di iniziare il calcolo.
    """
    assert _esegui(_scenario(tmp_path))["G12"].superato


def test_g12_fallisce_su_un_checksum_diverso(tmp_path):
    """
    **Obiettivo**: Verificare che alterare il contenuto di ``riferimento.fa.gz``
    faccia fallire G12 con codice ``E-S0-12`` citando ``tax.ref_md5``.

    **Razionale Scientifico/Sistemistico**: Impedisce di eseguire l'assegnazione
    tassonomica Bayesiana (Fase S8) su un archivio FASTA SILVA corrotto,
    incompleto o appartenente a una release diversa da quella dichiarata.
    """
    scenario = _scenario(tmp_path)
    (scenario.radice / "riferimento.fa.gz").write_bytes(b">seq1\nTTTT\n")

    esito = _esegui(scenario)["G12"]
    assert not esito.superato
    assert {v.codice for v in esito.violazioni} == {"E-S0-12"}
    assert "tax.ref_md5" in str(esito.violazioni[0])


def test_g12_fallisce_se_il_riferimento_non_esiste(tmp_path):
    """
    **Obiettivo**: Verificare che ``_g12_riferimento_verificato`` invocato su un
    file ``riferimento.fa.gz`` rimosso restituisca la violazione ``E-S0-12``.

    **Razionale Scientifico/Sistemistico**: Verifica l'autonomia difensiva del
    Gate G12 anche quando testato isolatamente rispetto a G01.
    """
    scenario = _scenario(tmp_path)
    (scenario.radice / "riferimento.fa.gz").unlink()
    # G01 lo intercetta per primo: qui si prova il gate da solo.
    from amplicon16s.gates.g01_g15 import _g12_riferimento_verificato

    violazioni, _ = _g12_riferimento_verificato(Contesto(scenario.config))
    assert violazioni and violazioni[0].codice == "E-S0-12"


# --------------------------------------------------------------------------- #
# G13 — gli archivi sono validi                                                #
# --------------------------------------------------------------------------- #


def test_g13_passa(tmp_path):
    """
    **Obiettivo**: Verificare che il Gate G13 sia superato su archivi
    ``*.fastq.gz`` ben formati con record FASTQ a 4 righe.

    **Razionale Scientifico/Sistemistico**: Certifica l'integrità dello stream
    ``gzip`` e della sintassi FASTQ su tutti i campioni.
    """
    assert _esegui(_scenario(tmp_path))["G13"].superato


def test_g13_fallisce_su_un_archivio_non_decomprimibile(tmp_path):
    """
    **Obiettivo**: Verificare che un file ``*.fastq.gz`` contenente byte non
    compressi in formato ``gzip`` faccia fallire G13 con codice ``E-S0-13``.

    **Razionale Scientifico/Sistemistico**: Intercetta in S0 file FASTQ corrotti
    da trasferimenti di rete interrotti prima che facciano andare in crash
    ``filterAndTrim`` di DADA2 a metà della Fase S2.
    """
    campioni = _campioni()
    campioni[0].contenuto_grezzo = b"questo non e' un archivio gzip"
    esito = _esegui(_scenario(tmp_path, campioni))["G13"]

    assert not esito.superato
    assert {v.codice for v in esito.violazioni} == {"E-S0-13"}
    assert "archivio non leggibile" in str(esito.violazioni[0])


def test_g13_fallisce_su_un_record_troncato(tmp_path):
    """
    **Obiettivo**: Verificare che un archivio ``gzip`` valido ma contenente un
    record FASTQ incompleto (2 righe invece di 4) faccia fallire G13 con ``E-S0-13``.

    **Razionale Scientifico/Sistemistico**: Garantisce il rispetto della
    grammatica FASTQ a 4 righe per record (``@header``, sequenza, ``+``, qualità
    Phred), intercettando troncamenti avvenuti prima della compressione.
    """
    campioni = _campioni()
    campioni[0].contenuto_grezzo = gzip.compress(
        f"@r1\n{lettura()}\n+\n{'I' * 151}\n@r2\n{lettura()}\n".encode()
    )
    esito = _esegui(_scenario(tmp_path, campioni))["G13"]

    assert not esito.superato
    assert "record troncato" in str(esito.violazioni[0])


def test_g13_fallisce_se_la_qualita_ha_lunghezza_diversa(tmp_path):
    """
    **Obiettivo**: Verificare che un record FASTQ in cui la stringa di qualità
    Phred (10 caratteri) ha lunghezza diversa dalla sequenza nucleotidica (151 bp)
    faccia fallire G13 con ``E-S0-13``.

    **Razionale Scientifico/Sistemistico**: In DADA2 il modello parametrico di
    errore (Fase S3) associa la probabilità di sostituzione base per base alla
    relativa qualità Phred; un disallineamento di lunghezza corromperebbe
    l'inferenza delle ASV.
    """
    campioni = _campioni()
    campioni[0].contenuto_grezzo = gzip.compress(
        f"@r1\n{lettura()}\n+\n{'I' * 10}\n".encode()
    )
    esito = _esegui(_scenario(tmp_path, campioni))["G13"]

    assert not esito.superato
    assert "lunghezze diverse" in str(esito.violazioni[0])


# --------------------------------------------------------------------------- #
# G14 — le risorse richieste sono disponibili                                  #
# --------------------------------------------------------------------------- #


def test_g14_passa(tmp_path):
    """
    **Obiettivo**: Verificare che il Gate G14 sia superato quando ``run.threads``
    e lo spazio disco disponibile in ``io.out_root`` soddisfano i requisiti.

    **Razionale Scientifico/Sistemistico**: Accerta la disponibilità delle
    risorse hardware prima di avviare fasi multithread e scritture pesanti.
    """
    assert _esegui(_scenario(tmp_path))["G14"].superato


def test_g14_fallisce_se_i_thread_eccedono_le_cpu(tmp_path):
    """
    **Obiettivo**: Verificare che richiedere ``run.threads = 100000`` (superiore
    alle CPU affini del processo restituite da ``sched_getaffinity``) faccia
    fallire G14 con codice ``E-S0-14``.

    **Razionale Scientifico/Sistemistico**: Previene l'*oversubscription*
    massiva dei thread OpenMP/Rcpp in DADA2 e DECIPHER dentro container o code
    HPC con cpuset limitati, che degraderebbe drasticamente le prestazioni per
    context-switching e consumo di RAM per-thread.
    """
    scenario = _scenario(tmp_path, sovrascrivi={"run": {"threads": 100000}})
    esito = _esegui(scenario)["G14"]

    assert not esito.superato
    assert {v.codice for v in esito.violazioni} == {"E-S0-14"}
    assert "run.threads" in str(esito.violazioni[0])


# --------------------------------------------------------------------------- #
# G15 — coerenza della configurazione                                          #
# --------------------------------------------------------------------------- #


def test_g15_passa(tmp_path):
    """
    **Obiettivo**: Verificare che il Gate G15 sia superato sulla configurazione
    predefinita dello scenario.

    **Razionale Scientifico/Sistemistico**: Conferma che tutte le relazioni
    algebriche inter-parametro del Gate G15 sono rispettate.
    """
    assert _esegui(_scenario(tmp_path))["G15"].superato


def test_g15_fallisce_su_una_soglia_degenere(tmp_path):
    """
    **Obiettivo**: Verificare che ``decontam.threshold = 1.0`` faccia fallire
    G15 con ``E-G15-02`` senza nemmeno eseguire ``G01`` (``esiti["G01"].eseguito is False``).

    **Razionale Scientifico/Sistemistico**: Una soglia ``decontam`` pari a 1.0
    classificherebbe come contaminante il 100% delle ASV presenti nei campioni;
    l'arresto immediato in G15 dimostra che la pipeline non apre alcun file su
    disco quando la configurazione è matematicamente degenere.
    """
    scenario = _scenario(tmp_path, sovrascrivi={"decontam": {"threshold": 1.0}})
    esiti = _esegui(scenario)

    assert not esiti["G15"].superato
    assert {v.codice for v in esiti["G15"].violazioni} == {"E-G15-02"}
    assert not esiti["G01"].eseguito, "la configurazione si controlla per prima"


# --------------------------------------------------------------------------- #
# La fase S0                                                                   #
# --------------------------------------------------------------------------- #


def test_s0_scrive_gli_artefatti(tmp_path):
    """
    **Obiettivo**: Verificare che ``esegui_s0`` produca esattamente i 4 file
    ``gates.json``, ``crosswalk.tsv``, ``inventario.json`` e
    ``letture_ispezionate.tsv`` nella cartella ``01_input_validation/``.

    **Razionale Scientifico/Sistemistico**: Questi 4 artefatti costituiscono il
    contratto di dati tra la validazione Python in S0 e tutte le fasi successive
    della pipeline (S1–S14).
    """
    scenario = _scenario(tmp_path, con_arricchimento=True)
    risultato = esegui_s0(scenario.config)

    assert risultato.superata
    prodotti = {p.name for p in risultato.artefatti}
    assert prodotti == {
        "gates.json", "crosswalk.tsv", "inventario.json", "letture_ispezionate.tsv"
    }
    for percorso in risultato.artefatti:
        assert percorso.is_file()
        assert percorso.parent.name == "01_input_validation"


def test_s0_registra_l_esito_di_ogni_gate(tmp_path):
    """
    **Obiettivo**: Verificare che ``gates.json`` contenga ``superata: true`` e
    l'elenco completo dei 15 gate con ``superato: true``.

    **Razionale Scientifico/Sistemistico**: Assicura la tracciabilità completa
    dell'audit pre-analitico per il report finale della Fase S14.
    """
    scenario = _scenario(tmp_path)
    risultato = esegui_s0(scenario.config)

    documento = json.loads(
        next(p for p in risultato.artefatti if p.name == "gates.json").read_text()
    )
    assert documento["superata"] is True
    assert len(documento["gate"]) == 15
    assert {g["gate"] for g in documento["gate"]} == set(nomi_dei_gate())
    assert all(g["superato"] for g in documento["gate"])


def test_s0_registra_il_denominatore_di_prevalenza(tmp_path):
    """
    **Obiettivo**: Verificare che ``inventario.json`` distingua il numero totale
    di campioni (``campioni == 4``) dal ``denominatore_prevalenza == 2``
    calcolato sui soli campioni biologici.

    **Razionale Scientifico/Sistemistico**: Sul dataset reale OSD-734 (960
    campioni totali), il filtro di prevalenza dell'1% in S13 deve essere
    calcolato sui soli **803 campioni biologici** (escludendo i controlli
    positivi e negativi); questo test verifica l'esattezza della formula di
    conteggio registrata in ``inventario.json``.
    """
    scenario = _scenario(tmp_path)
    risultato = esegui_s0(scenario.config)

    documento = json.loads(
        next(p for p in risultato.artefatti if p.name == "inventario.json").read_text()
    )
    assert documento["campioni"] == 4
    assert documento["denominatore_prevalenza"] == 2  # i soli biologici


def test_s0_registra_il_limite_della_stima_sulle_lunghezze(tmp_path):
    """
    **Obiettivo**: Verificare che ``letture_ispezionate.tsv`` includa la colonna
    booleana ``file_esaurito`` e una riga per ciascuno dei 4 campioni (più l'intestazione).

    **Razionale Scientifico/Sistemistico**: Documenta per ogni singolo file
    FASTQ se le statistiche sulle lunghezze e sui primer (G09/G10) coprono
    l'intero archivio (``file_esaurito == True``) oppure un campione in
    streaming delle prime ``qc.head_reads`` letture.
    """
    scenario = _scenario(tmp_path)
    risultato = esegui_s0(scenario.config)

    righe = next(
        p for p in risultato.artefatti if p.name == "letture_ispezionate.tsv"
    ).read_text(encoding="utf-8").splitlines()
    assert righe[0].split("\t")[3] == "file_esaurito"
    assert len(righe) == 5


def test_s0_gli_artefatti_entrano_nel_manifesto(tmp_path):
    """
    **Obiettivo**: Verificare che al termine di ``esegui_s0`` la cartella
    ``01_input_validation`` risulti ``fase_completa == True`` e ``non_integri == ()``.

    **Razionale Scientifico/Sistemistico**: Garantisce che i checksum SHA-256
    dei 4 file prodotti in S0 siano registrati nel manifesto per abilitare il
    meccanismo di ripresa (*resume*).
    """
    from amplicon16s.io_layer.artifacts import AlberoOutput, Fase

    scenario = _scenario(tmp_path)
    esegui_s0(scenario.config)

    albero = AlberoOutput(scenario.config.io.out_root)
    assert albero.fase_completa(Fase.INPUT_VALIDATION)
    assert albero.non_integri(Fase.INPUT_VALIDATION) == ()


def test_s0_solleva_quando_un_gate_fallisce(tmp_path):
    """
    **Obiettivo**: Verificare che per difetto (``solleva=True``) ``esegui_s0``
    sollevi ``ErroreGate`` indicando il gate fallito (``G14``).

    **Razionale Scientifico/Sistemistico**: Impedisce che la pipeline prosegua
    verso S1 quando la validazione pre-analitica non è superata.
    """
    scenario = _scenario(tmp_path, sovrascrivi={"run": {"threads": 100000}})
    with pytest.raises(ErroreGate) as errore:
        esegui_s0(scenario.config)
    assert errore.value.gate == "G14"


def test_s0_puo_riportare_senza_sollevare(tmp_path):
    """
    **Obiettivo**: Verificare che invocando ``esegui_s0(..., solleva=False)`` su
    uno scenario invalido venga comunque scritto ``gates.json`` e restituito
    ``risultato.superata is False`` senza sollevare eccezioni.

    **Razionale Scientifico/Sistemistico**: Consente al comando CLI
    ``amplicon16s validate`` e ai report diagnostici di salvare su disco il
    referto JSON del fallimento anche quando un gate blocca la corsa.
    """
    scenario = _scenario(tmp_path, sovrascrivi={"run": {"threads": 100000}})
    risultato = esegui_s0(scenario.config, solleva=False)

    assert not risultato.superata
    assert [e.gate for e in risultato.falliti] == ["G14"]
    assert any(p.name == "gates.json" for p in risultato.artefatti)


def test_s0_registra_gli_avvisi_senza_fallire(tmp_path):
    """
    **Obiettivo**: Verificare che gli avvisi non bloccanti (es. ``E-S0-15`` dal
    Gate G08) siano conservati in ``risultato.avvisi`` mantenendo ``risultato.superata is True``.

    **Razionale Scientifico/Sistemistico**: Assicura che le condizioni di
    degradazione controllata vengano trasmesse alle fasi a valle e al registro
    degli eventi senza interrompere l'esecuzione di S0.
    """
    scenario = _scenario(
        tmp_path, con_arricchimento=True,
        sovrascrivi={"decontam": {"min_blanks": 3}},
    )
    risultato = esegui_s0(scenario.config)

    assert risultato.superata
    da_g08 = [a for a in risultato.avvisi if a.codice == "E-S0-15"]
    assert len(da_g08) == 1
    assert "controlli negativi" in str(da_g08[0])
