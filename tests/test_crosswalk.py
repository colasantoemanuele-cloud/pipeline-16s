"""Suite di verifica per il crosswalk dei metadati, l'arricchimento lotti e i Gate G04, G05, G06 e G11.

Inquadramento nel Piano Operativo:
    - **Settimana di riferimento**: **Settimana 6 (W6 — Fase F2: Crosswalk,
      inventario dei campioni e classificazione dei controlli)**.
    - **Scopo del modulo**: Verifica su scenari sintetici controllati la
      costruzione dell'inventario campioni dall'incrocio fra gli archivi
      ``.fastq.gz`` su disco, l'Assay Table 16S, la Study Sample Table ISA-Tab e
      la tabella opzionale dei lotti (``batch_table``). Collauda specificamente:
        * L'estrazione non ambigua dell'accession ENA/SRA tramite ``re.finditer``
          (Gate **G04** → ``E-S0-04``);
        * L'univocità degli accession su disco e nella tabella di assay
          (Gate **G05** → ``E-S0-05``);
        * La corrispondenza biunivoca 1:1 (insiemi simmetrici senza orfani) fra
          file FASTQ e righe dell'Assay Table (Gate **G06** → ``E-S0-06``);
        * La classificazione esaustiva in ``BIOLOGICO``, ``CONTROLLO_POSITIVO``
          e ``CONTROLLO_NEGATIVO`` (Gate **G11** → ``E-S0-11``);
        * Il vincolo architetturale per cui la chiave primaria è sempre
          l'``accession`` e mai il ``Sample Name`` (che collasserebbe le repliche
          tecniche), e la coerenza nell'attribuzione del modulo spaziale ISS.
    - **Moduli sorgente coperti**:
        * ``src/amplicon16s/metadata/models.py``
        * ``src/amplicon16s/metadata/crosswalk.py``
        * ``src/amplicon16s/metadata/controls_map.py``
        * ``src/amplicon16s/gates/g01_g15.py`` (limitatamente a ``G04, G05, G06, G11``)
"""

from __future__ import annotations

import re

import pytest
from conftest import BIOLOGICO, NEGATIVO, POSITIVO, Campione, crea_scenario

from amplicon16s.gates.g01_g15 import ErroreGate, esegui_gate_metadati
from amplicon16s.metadata.controls_map import MappaControlli
from amplicon16s.metadata.crosswalk import (
    AccessionNonEstraibile,
    analizza,
    estrai_accession,
)
from amplicon16s.metadata.models import ClasseCampione

ACCESSION = re.compile(r"(E|S|D)RX[0-9]{4,}")


def _tre_campioni() -> list[Campione]:
    """Costruisce una terna minima bilanciata: 1 biologico, 1 controllo positivo e 1 negativo."""
    return [
        Campione("ERX1000001", "NOD1D4.L1"),
        Campione("ERX1000002", "POS.P1.1", materiale=POSITIVO, posizione="Not Applicable"),
        Campione("ERX1000003", "BLANK.P1.1", materiale=NEGATIVO, posizione="Not Applicable"),
    ]


# --------------------------------------------------------------------------- #
# Estrazione dell'accession                                                    #
# --------------------------------------------------------------------------- #


def test_estrae_l_accession_anche_con_gruppi_di_cattura():
    """
    **Obiettivo**: Verificare che ``estrai_accession()`` restituisca l'intero
    identificativo ``"ERX12083297"`` anche quando l'espressione regolare contiene
    gruppi di cattura parentetici come ``(E|S|D)RX[0-9]{4,}``.

    **Razionale Scientifico/Sistemistico**: In Python, ``re.findall()`` su una
    regex dotata di gruppi di cattura restituisce solo la sottostringa catturata
    dalle parentesi (``"E"`` anziché ``"ERX12083297"``), facendo collassare tutti
    i 960 file di OSD-734 sull'unica chiave ``"E"``. L'uso di ``re.finditer()``
    e ``match.group(0)`` garantisce l'estrazione dell'accession completo.
    """
    assert estrai_accession("GLDS-653_Amplicon_ERX12083297_raw.fastq.gz", ACCESSION) == (
        "ERX12083297"
    )


def test_nessun_accession_e_un_errore():
    """
    **Obiettivo**: Verificare che ``estrai_accession()`` sollevi
    ``AccessionNonEstraibile`` quando la stringa non contiene alcun codice ENA/SRA.

    **Razionale Scientifico/Sistemistico**: Impedisce che file estranei presenti
    nella cartella ``fastq_dir`` (es. file di riepilogo o archivi rinominati
    manualmente) vengano associati a una chiave vuota o ``None``.
    """
    with pytest.raises(AccessionNonEstraibile, match="nessun accession"):
        estrai_accession("campione_senza_codice.fastq.gz", ACCESSION)


def test_due_accession_sono_un_errore_quanto_nessuno():
    """
    **Obiettivo**: Verificare che la presenza di due codici accession distinti
    nello stesso nome di file (``ERX1000001_ERX1000002.fastq.gz``) sollevi
    ``AccessionNonEstraibile`` per ambiguità.

    **Razionale Scientifico/Sistemistico**: Se l'estrattore utilizzasse
    ``re.search()`` limitandosi a prendere la prima corrispondenza, un file nato
    dalla concatenazione accidentale di due campioni verrebbe silenziosamente
    attribuito al primo accession, fondendo due profili microbici distinti.
    """
    with pytest.raises(AccessionNonEstraibile, match="ambiguo"):
        estrai_accession("ERX1000001_ERX1000002.fastq.gz", ACCESSION)


# --------------------------------------------------------------------------- #
# Lo scenario corretto passa tutti i gate                                      #
# --------------------------------------------------------------------------- #


def test_scenario_corretto_produce_l_inventario(tmp_path):
    """
    **Obiettivo**: Verificare che uno scenario ben formato superi tutti i gate
    dei metadati (`G03, G04, G05, G06, G11`) producendo l'oggetto ``Inventario``
    con i conteggi esatti per ciascuna delle tre classi.

    **Razionale Scientifico/Sistemistico**: Certifica il percorso nominale di
    costruzione dell'inventario che alimenta tutte le fasi successive della
    pipeline (da `S1` a `S14`).
    """
    scenario = crea_scenario(tmp_path, _tre_campioni())
    inventario = esegui_gate_metadati(scenario.config)

    assert len(inventario) == 3
    assert inventario.conteggi() == {
        ClasseCampione.BIOLOGICO: 1,
        ClasseCampione.CONTROLLO_POSITIVO: 1,
        ClasseCampione.CONTROLLO_NEGATIVO: 1,
    }
    assert sorted(inventario.accessioni) == ["ERX1000001", "ERX1000002", "ERX1000003"]


def test_l_inventario_conserva_il_valore_grezzo_del_materiale(tmp_path):
    """
    **Obiettivo**: Verificare che ogni voce dell'``Inventario`` conservi, oltre
    alla classe normalizzata ``ClasseCampione``, anche la stringa originale
    ``materiale`` letta dalla Study Table.

    **Razionale Scientifico/Sistemistico**: Mantiene la tracciabilità di audit
    in ``crosswalk.tsv`` e nei metadati dell'oggetto ``phyloseq`` finale (`S10`),
    permettendo di distinguere sottotipi di campioni biologici (es. ``Surface swab``
    contro ``Air Sample``) o di controlli negativi.
    """
    scenario = crea_scenario(tmp_path, _tre_campioni())
    inventario = esegui_gate_metadati(scenario.config)
    assert inventario["ERX1000002"].materiale == POSITIVO


def test_ogni_campione_e_legato_al_proprio_file(tmp_path):
    """
    **Obiettivo**: Verificare che ogni ``CampioneInventario`` punti al percorso
    fisico del proprio file ``.fastq.gz`` contenente il relativo accession.

    **Razionale Scientifico/Sistemistico**: Garantisce che non avvengano
    trasposizioni di indice o riordinamenti spuri durante l'accoppiamento tra
    tabella di assay e filesystem.
    """
    scenario = crea_scenario(tmp_path, _tre_campioni())
    for campione in esegui_gate_metadati(scenario.config):
        assert campione.file is not None
        assert campione.accession in campione.file.name


# --------------------------------------------------------------------------- #
# G04 — l'accession è estraibile e non ambiguo                                 #
# --------------------------------------------------------------------------- #


def test_g04_passa_sullo_scenario_corretto(tmp_path):
    """
    **Obiettivo**: Verificare il superamento del Gate G04 quando tutti i file e
    le righe di assay presentano un accession univoco ed estraibile.

    **Razionale Scientifico/Sistemistico**: Conferma l'assenza di falsi positivi
    del Gate G04 su nomenclature standard NASA GeneLab.
    """
    scenario = crea_scenario(tmp_path, _tre_campioni())
    assert len(esegui_gate_metadati(scenario.config)) == 3


def test_g04_fallisce_se_un_nome_di_file_non_contiene_accession(tmp_path):
    """
    **Obiettivo**: Verificare che la presenza di un file FASTQ privo di codice
    accession faccia fallire il Gate G04 con codice ``E-S0-04`` nominando il file.

    **Razionale Scientifico/Sistemistico**: Blocca immediatamente l'esecuzione
    se nella directory di input sono presenti file con naming non conforme che
    non potrebbero essere collegati ai metadati clinici o ambientali.
    """
    campioni = _tre_campioni()
    campioni[0].file = "campione_senza_codice.fastq.gz"
    scenario = crea_scenario(tmp_path, campioni)

    with pytest.raises(ErroreGate) as errore:
        esegui_gate_metadati(scenario.config)

    assert errore.value.gate == "G04"
    assert errore.value.codice == "E-S0-04"
    assert "campione_senza_codice.fastq.gz" in str(errore.value)


def test_g04_fallisce_se_l_accession_e_ambiguo(tmp_path):
    """
    **Obiettivo**: Verificare che un file il cui nome contiene due accession
    distinti faccia fallire il Gate G04 con codice ``E-S0-04`` e menzione di ambiguità.

    **Razionale Scientifico/Sistemistico**: Impedisce che file con nomi ambigui
    superino la validazione iniziale S0.
    """
    campioni = _tre_campioni()
    campioni[0].file = "ERX1000001_ERX1000009_doppio.fastq.gz"
    scenario = crea_scenario(tmp_path, campioni)

    with pytest.raises(ErroreGate) as errore:
        esegui_gate_metadati(scenario.config)

    assert errore.value.codice == "E-S0-04"
    assert "ambiguo" in str(errore.value)


# --------------------------------------------------------------------------- #
# G05 — gli accession sono univoci                                             #
# --------------------------------------------------------------------------- #


def test_g05_fallisce_su_accession_duplicati_fra_i_file(tmp_path):
    """
    **Obiettivo**: Verificare che la presenza su disco di due file FASTQ aventi
    lo stesso accession (``ERX1000001``) faccia fallire il Gate G05 con ``E-S0-05``.

    **Razionale Scientifico/Sistemistico**: Previene il rischio che copie di
    backup o file rielaborati nella stessa directory facciano processare due
    volte la stessa corsa o sovrascrivano silenziosamente le letture in S2.
    """
    scenario = crea_scenario(
        tmp_path, _tre_campioni(), file_in_piu=["ERX1000001_copia.fastq.gz"]
    )

    with pytest.raises(ErroreGate) as errore:
        esegui_gate_metadati(scenario.config)

    assert errore.value.gate == "G05"
    assert errore.value.codice == "E-S0-05"
    assert "ERX1000001" in str(errore.value)


def test_g05_fallisce_su_accession_duplicati_nell_assay(tmp_path):
    """
    **Obiettivo**: Verificare che due righe dell'Assay Table che puntano allo
    stesso accession ``ERX1000001`` facciano fallire il Gate G05 con ``E-S0-05``.

    **Razionale Scientifico/Sistemistico**: Impedisce che uno stesso file di
    sequenziamento venga associato per errore a due campioni biologici diversi
    nella tabella dei metadati, generando un profilo duplicato artificiale.
    """
    campioni = _tre_campioni()
    campioni.append(Campione("ERX1000001", "ALTRO.NOME", file=""))
    scenario = crea_scenario(tmp_path, campioni)

    with pytest.raises(ErroreGate) as errore:
        esegui_gate_metadati(scenario.config)

    assert errore.value.codice == "E-S0-05"
    assert "tabella di assay" in str(errore.value)


# --------------------------------------------------------------------------- #
# G06 — gli insiemi coincidono                                                 #
# --------------------------------------------------------------------------- #


def test_g06_fallisce_su_un_file_senza_riga_nell_assay(tmp_path):
    """
    **Obiettivo**: Verificare che un file FASTQ presente su disco ma non
    dichiarato nell'Assay Table (``ERX1000099_orfano.fastq.gz``) faccia fallire
    il Gate G06 con codice ``E-S0-06``.

    **Razionale Scientifico/Sistemistico**: Garantisce la simmetria 1:1 tra
    storage e metadati: un file FASTQ orfano privo di metadati non saprebbe se
    essere trattato come campione biologico o controllo negativo durante `decontam`.
    """
    scenario = crea_scenario(
        tmp_path, _tre_campioni(), file_in_piu=["ERX1000099_orfano.fastq.gz"]
    )

    with pytest.raises(ErroreGate) as errore:
        esegui_gate_metadati(scenario.config)

    assert errore.value.gate == "G06"
    assert errore.value.codice == "E-S0-06"
    assert "ERX1000099" in str(errore.value)
    assert "ma non nella tabella di assay" in str(errore.value)


def test_g06_fallisce_su_una_riga_dell_assay_senza_file(tmp_path):
    """
    **Obiettivo**: Verificare che una riga dichiarata nell'Assay Table il cui
    file FASTQ manca dalla cartella ``fastq_dir`` (``ERX1000004``) faccia
    fallire il Gate G06 con codice ``E-S0-06``.

    **Razionale Scientifico/Sistemistico**: Intercetta download incompleti dal
    repository NASA GeneLab / ENA prima di iniziare la corsa: se mancasse ad
    esempio il file FASTQ di un controllo negativo, una piastra potrebbe scendere
    sotto la soglia ``decontam.min_blanks`` senza che il ricercatore se ne accorga.
    """
    campioni = _tre_campioni()
    campioni.append(Campione("ERX1000004", "SENZA.FILE", file=""))
    scenario = crea_scenario(tmp_path, campioni)

    with pytest.raises(ErroreGate) as errore:
        esegui_gate_metadati(scenario.config)

    assert errore.value.codice == "E-S0-06"
    assert "ERX1000004" in str(errore.value)
    assert "ma non fra i file" in str(errore.value)


# --------------------------------------------------------------------------- #
# G11 — ogni campione ricade in una classe dichiarata                          #
# --------------------------------------------------------------------------- #


def test_g11_fallisce_su_un_materiale_non_mappato(tmp_path):
    """
    **Obiettivo**: Verificare che un campione dell'assay 16S avente un
    ``Material Type`` non elencato in ``ctrl.*_values`` (es. ``"solvent control"``)
    faccia fallire il Gate G11 con codice ``E-S0-11``.

    **Razionale Scientifico/Sistemistico**: Ogni campione processato dalla
    pipeline deve avere un ruolo noto (biologico, positivo o negativo) per
    calcolare correttamente il denominatore di prevalenza e i modelli di
    decontaminazione `decontam` e `KatharoSeq`.
    """
    campioni = _tre_campioni()
    campioni[0].materiale = "solvent control"
    scenario = crea_scenario(tmp_path, campioni)

    with pytest.raises(ErroreGate) as errore:
        esegui_gate_metadati(scenario.config)

    assert errore.value.gate == "G11"
    assert errore.value.codice == "E-S0-11"
    assert "solvent control" in str(errore.value)
    assert "ctrl.biological_values" in str(errore.value)


def test_g11_fallisce_su_un_materiale_vuoto(tmp_path):
    """
    **Obiettivo**: Verificare che una cella vuota nella colonna ``ctrl.column``
    della Study Table faccia fallire il Gate G11 con ``E-S0-11``.

    **Razionale Scientifico/Sistemistico**: Impedisce che campioni privi di
    annotazione nella Study Table vengano silenziosamente ignorati o assegnati
    a una categoria predefinita arbitraria.
    """
    campioni = _tre_campioni()
    campioni[0].materiale = ""
    scenario = crea_scenario(tmp_path, campioni)

    with pytest.raises(ErroreGate, match="E-S0-11"):
        esegui_gate_metadati(scenario.config)


def test_la_classificazione_tollera_maiuscole_e_spazi(tmp_path):
    """
    **Obiettivo**: Verificare che ``MappaControlli`` classifichi correttamente
    ``"  positive control  "`` come ``CONTROLLO_POSITIVO`` ignorando differenze
    di maiuscole/minuscole e spazi bianchi periferici.

    **Razionale Scientifico/Sistemistico**: Le tabelle ISA-Tab sono compilate
    manualmente da curatori diversi e spesso presentano variazioni tipografiche
    (``"Positive Control"`` vs ``"positive control "``) che non devono causare
    falsi blocchi del Gate G11.
    """
    campioni = _tre_campioni()
    campioni[1].materiale = "  positive control  "
    scenario = crea_scenario(tmp_path, campioni)
    inventario = esegui_gate_metadati(scenario.config)
    assert inventario["ERX1000002"].classe is ClasseCampione.CONTROLLO_POSITIVO


def test_la_mappa_non_inventa_una_categoria_di_ripiego(tmp_path):
    """
    **Obiettivo**: Verificare che ``mappa.classifica()`` restituisca ``None``
    su etichette sconosciute o valori ``None``, senza mai ripiegare su ``BIOLOGICO``.

    **Razionale Scientifico/Sistemistico**: Un fallback automatico a ``BIOLOGICO``
    per le stringhe non riconosciute farebbe entrare eventuali nuovi tipi di
    controllo negativo tra i campioni biologici, falsando il denominatore di
    prevalenza e nascondendo i contaminanti a `decontam`.
    """
    scenario = crea_scenario(tmp_path, _tre_campioni())
    mappa = MappaControlli.da_configurazione(scenario.config.ctrl)
    assert mappa.classifica("materiale sconosciuto") is None
    assert mappa.classifica(None) is None


# --------------------------------------------------------------------------- #
# La chiave è l'accession, non il nome del campione                            #
# --------------------------------------------------------------------------- #


def test_il_nome_del_campione_si_ripete_fra_repliche_l_accession_no(tmp_path):
    """
    **Obiettivo**: Verificare che due repliche tecniche con lo stesso
    ``Sample Name`` (``LAB1P3.L1``) ma accession distinti (``ERX1000001`` ed
    ``ERX1000002``) rimangano due entità separate nel crosswalk mentre
    collasserebbero in una sola se indicizzate per nome.

    **Razionale Scientifico/Sistemistico**: In studi come OSD-734 alcune
    superfici spaziali sono state sequenziate in replicato tecnico (stesso
    campione biologico nella Study Table, due corse distinte nell'Assay Table).
    Usare il ``Sample Name`` come chiave primaria farebbe collassare le due
    repliche su un'unica voce del dizionario, scartando metà dei dati di
    sequenziamento o accoppiando il file FASTQ sbagliato alla piastra sbagliata.
    """
    campioni = [
        Campione("ERX1000001", "LAB1P3.L1", posizione="LAB1P3"),
        Campione("ERX1000002", "LAB1P3.L1", posizione="LAB1P3"),
    ]
    scenario = crea_scenario(tmp_path, campioni)

    analisi = analizza(scenario.config)
    assert analisi.accession_ripetuti_nell_assay == {}

    per_nome = {c.nome for c in campioni}
    per_accession = {c.accession for c in campioni}
    assert len(per_nome) == 1, "il nome non distingue le repliche"
    assert len(per_accession) == 2, "l'accession le distingue"


def test_appaiare_per_nome_produce_una_corrispondenza_errata(tmp_path):
    """
    **Obiettivo**: Dimostrare sperimentalmente che un dizionario indicizzato per
    ``Sample Name`` sovrascrive la classe del primo campione (``BIOLOGICO``) con
    quella del secondo (``POSITIVO``), mentre l'indicizzazione per ``accession``
    preserva l'identità e la classe di entrambi.

    **Razionale Scientifico/Sistemistico**: Prova formalmente perché l'intera
    architettura di ``crosswalk.py`` utilizza l'``accession`` come chiave
    primaria: un appaiamento per nome attribuirebbe la classe di un controllo
    al file FASTQ di un campione biologico (o viceversa), corrompendo
    silenziosamente la decontaminazione S12 e la calibrazione KatharoSeq S11.
    """
    campioni = [
        Campione("ERX1000001", "RIPETUTO", materiale=BIOLOGICO),
        Campione("ERX1000002", "RIPETUTO", materiale=POSITIVO, file="ERX1000002_b.fastq.gz"),
    ]
    scenario = crea_scenario(tmp_path, campioni)

    per_nome = {}
    for campione in campioni:
        per_nome[campione.nome] = campione.materiale
    assert len(per_nome) == 1
    assert per_nome["RIPETUTO"] == POSITIVO

    per_accession = {c.accession: c.materiale for c in campioni}
    assert per_accession["ERX1000001"] == BIOLOGICO
    assert per_accession["ERX1000002"] == POSITIVO


# --------------------------------------------------------------------------- #
# Derivazione del modulo                                                       #
# --------------------------------------------------------------------------- #


def test_il_modulo_si_deriva_dalla_posizione(tmp_path):
    """
    **Obiettivo**: Verificare che dalla stringa di posizione ``"NOD3O2 (ARED)"``
    venga estratto tramite ``meta.module_regex`` il prefisso del modulo ``"NOD3"``.

    **Razionale Scientifico/Sistemistico**: Permette di ricostruire il modulo
    abitativo della Stazione Spaziale Internazionale (ISS) direttamente dalla
    Study Table standard quando non viene fornito il file opzionale ``batch_table``.
    """
    campioni = [Campione("ERX1000001", "QUALUNQUE", posizione="NOD3O2 (ARED)")]
    scenario = crea_scenario(tmp_path, campioni)
    assert esegui_gate_metadati(scenario.config)["ERX1000001"].modulo == "NOD3"


@pytest.mark.parametrize("posizione", ["Air Sample", "Not Applicable", ""])
def test_le_posizioni_che_non_sono_superfici_non_hanno_modulo(posizione, tmp_path):
    """
    **Obiettivo**: Verificare che i campioni con posizione ``"Air Sample"``,
    ``"Not Applicable"`` o vuota ricevano ``modulo = None`` anche se il loro
    ``Sample Name`` inizia con il prefisso di un modulo (``"NOD1.F3"``).

    **Razionale Scientifico/Sistemistico**: I campioni d'aria e i controlli non
    rappresentano superfici di un modulo ISS; attribuire loro un modulo per
    via del prefisso nel nome mescolerebbe aerobioma e microbioma di superficie
    nelle analisi ecologiche stratificate per modulo (`amplicon16s_eco`).
    """
    campioni = [Campione("ERX1000001", "NOD1.F3", posizione=posizione)]
    scenario = crea_scenario(tmp_path, campioni)
    assert esegui_gate_metadati(scenario.config)["ERX1000001"].modulo is None


def test_il_modulo_non_si_deriva_dal_nome(tmp_path):
    """
    **Obiettivo**: Verificare che l'estrazione del modulo legga esclusivamente
    la colonna ``meta.module_column`` (``Factor Value[Sample Location]``) e mai
    il ``Sample Name``.

    **Razionale Scientifico/Sistemistico**: In OSD-734 alcuni filtri d'aria o
    controlli sul campo portano nel nome il codice del nodo dove è stato
    posizionato il campionatore, ma la loro natura non-superficiale è dichiarata
    nella colonna Location (``Air Sample``): ignorare il nome evita false
    attribuzioni spaziali.
    """
    campioni = [Campione("ERX1000001", "NOD1.F3", posizione="Air Sample")]
    scenario = crea_scenario(tmp_path, campioni)
    campione = esegui_gate_metadati(scenario.config)["ERX1000001"]
    assert campione.nome.startswith("NOD1")
    assert campione.modulo is None


def test_la_derivazione_si_puo_disattivare(tmp_path):
    """
    **Obiettivo**: Verificare che impostando ``meta.derive_module = False``
    tutti i campioni abbiano ``modulo = None`` in assenza di ``batch_table``.

    **Razionale Scientifico/Sistemistico**: Consente di riutilizzare la pipeline
    su dataset 16S non spaziali in cui i codici di posizione non rappresentano
    moduli architettonici.
    """
    scenario = crea_scenario(
        tmp_path, _tre_campioni(), sovrascrivi={"meta": {"derive_module": False}}
    )
    assert all(c.modulo is None for c in esegui_gate_metadati(scenario.config))


# --------------------------------------------------------------------------- #
# Arricchimento del lotto: facoltativo                                         #
# --------------------------------------------------------------------------- #


def test_senza_arricchimento_piastra_e_corsa_restano_nulle(tmp_path):
    """
    **Obiettivo**: Verificare che omettendo ``io.batch_table`` l'inventario sia
    valido e riporti ``piastre = ()``, ``corse = ()`` e tutti i campioni in ``senza_lotto``.

    **Razionale Scientifico/Sistemistico**: Dimostra che ``batch_table`` è un
    arricchimento facoltativo: in sua assenza la pipeline degrada in modo
    controllato alla stima globale del modello d'errore DADA2 (S3) e alla
    decontaminazione globale (S12).
    """
    scenario = crea_scenario(tmp_path, _tre_campioni(), con_arricchimento=False)
    inventario = esegui_gate_metadati(scenario.config)

    assert len(inventario) == 3
    assert inventario.piastre == ()
    assert inventario.corse == ()
    assert len(inventario.senza_lotto) == 3


def test_con_arricchimento_piastra_e_corsa_sono_valorizzate(tmp_path):
    """
    **Obiettivo**: Verificare che fornendo ``io.batch_table`` i campi ``piastra``
    (per `decontam`) e ``corsa`` (per `learnErrors`) vengano popolati per ogni campione.

    **Razionale Scientifico/Sistemistico**: Permette a S3 di stimare i profili
    di errore Phred separatamente per ciascuna corsa di sequenziamento MiSeq e
    a S12 di applicare `decontam` piastra per piastra (batch-specific decontamination).
    """
    campioni = _tre_campioni()
    campioni[2].piastra = "2"
    campioni[2].corsa = "corsa_B"
    scenario = crea_scenario(tmp_path, campioni, con_arricchimento=True)
    inventario = esegui_gate_metadati(scenario.config)

    assert inventario.piastre == ("1", "2")
    assert inventario.corse == ("corsa_A", "corsa_B")
    assert inventario.senza_lotto == ()


def test_un_campione_non_coperto_dall_arricchimento_resta_senza_lotto(tmp_path):
    """
    **Obiettivo**: Verificare che un campione assente da ``batch_table`` resti
    con ``piastra = None`` e venga elencato in ``senza_riga_di_arricchimento``.

    **Razionale Scientifico/Sistemistico**: Traccia esplicitamente la copertura
    parziale del file di arricchimento senza inventare assegnazioni di piastra.
    """
    campioni = _tre_campioni()
    campioni[0].chiave_arricchimento = "ERX9999999"
    scenario = crea_scenario(tmp_path, campioni, con_arricchimento=True)

    inventario = esegui_gate_metadati(scenario.config)
    assert len(inventario) == 3
    assert len(inventario.senza_lotto) == 1
    assert analizza(scenario.config).senza_riga_di_arricchimento == ["ERX1000001"]


def test_una_chiave_ambigua_nell_arricchimento_non_viene_indovinata(tmp_path):
    """
    **Obiettivo**: Verificare che se un accession compare due volte nella
    ``batch_table``, ``inventario[acc].piastra`` resti ``None`` e l'accession
    compaia in ``arricchimento_ambiguo``.

    **Razionale Scientifico/Sistemistico**: Scegliere arbitrariamente la prima
    o l'ultima riga duplicata nella ``batch_table`` assegnerebbe il campione
    ai reagenti e ai blank della piastra sbagliata in S12.
    """
    campioni = _tre_campioni()
    campioni[1].chiave_arricchimento = campioni[0].accession
    scenario = crea_scenario(tmp_path, campioni, con_arricchimento=True)

    inventario = esegui_gate_metadati(scenario.config)
    assert inventario["ERX1000001"].piastra is None
    assert analizza(scenario.config).arricchimento_ambiguo == {"ERX1000001": 2}


# --------------------------------------------------------------------------- #
# La chiave dell'arricchimento è l'accession, non il nome                      #
# --------------------------------------------------------------------------- #


def test_l_arricchimento_si_aggancia_per_accession(tmp_path):
    """
    **Obiettivo**: Verificare che due repliche (`_rep1` su piastra 3/corsa A e
    `_rep2` su piastra 9/corsa B) ricevano ciascuna il proprio lotto distinto
    tramite join su ``accession``.

    **Razionale Scientifico/Sistemistico**: Nel dataset reale OSD-734 le repliche
    tecniche sono state estratte e sequenziate su piastre e corse diverse per
    controllare l'effetto batch: agganciare la ``batch_table`` per ``accession``
    è l'unico modo per associare a ciascuna replica il modello d'errore della
    sua vera corsa e i controlli negativi della sua vera piastra di estrazione.
    """
    campioni = [
        Campione("ERX1000001", "RIPETUTO_rep1", piastra="3", corsa="corsa_A"),
        Campione("ERX1000002", "RIPETUTO_rep2", piastra="9", corsa="corsa_B"),
    ]
    scenario = crea_scenario(tmp_path, campioni, con_arricchimento=True)
    inventario = esegui_gate_metadati(scenario.config)

    assert inventario["ERX1000001"].piastra == "3"
    assert inventario["ERX1000001"].corsa == "corsa_A"
    assert inventario["ERX1000002"].piastra == "9"
    assert inventario["ERX1000002"].corsa == "corsa_B"
    assert inventario.senza_lotto == ()


def test_un_file_di_arricchimento_senza_colonna_dell_accession_e_respinto(tmp_path):
    """
    **Obiettivo**: Verificare che una ``batch_table`` priva di colonna con
    codici accession faccia fallire G04 (`E-S0-04`) con istruzioni operative
    dettagliate su come ricavarla da ``io.assay_table``.

    **Razionale Scientifico/Sistemistico**: Impedisce l'uso di tabelle batch
    indicizzate solo per ``Sample Name`` e guida l'operatore nella correzione.
    """
    scenario = crea_scenario(
        tmp_path,
        _tre_campioni(),
        con_arricchimento=True,
        colonna_chiave_arricchimento="sample_name_ena",
    )

    with pytest.raises(ErroreGate) as errore:
        esegui_gate_metadati(scenario.config)

    testo = str(errore.value)
    assert errore.value.codice == "E-S0-04"
    assert "io.batch_table" in testo
    assert "experiment_accession" in testo
    assert "ricavala dalla tabella di assay" in testo
    assert "io.assay_table" in testo
    assert "togli io.batch_table" in testo


def test_una_chiave_senza_accession_e_respinta(tmp_path):
    """
    **Obiettivo**: Verificare che una riga della ``batch_table`` la cui chiave
    non contiene un accession ENA/SRA valido sollevi ``ErroreGate`` (`E-S0-04`).

    **Razionale Scientifico/Sistemistico**: Evita che righe spurie o malformate
    nella tabella dei lotti vengano ignorate silenziosamente.
    """
    campioni = _tre_campioni()
    campioni[0].chiave_arricchimento = "nome_senza_codice"
    scenario = crea_scenario(tmp_path, campioni, con_arricchimento=True)

    with pytest.raises(ErroreGate) as errore:
        esegui_gate_metadati(scenario.config)

    assert errore.value.codice == "E-S0-04"
    assert "nome_senza_codice" in str(errore.value)


# --------------------------------------------------------------------------- #
# Precedenza del modulo dichiarato                                             #
# --------------------------------------------------------------------------- #


def test_il_modulo_dichiarato_ha_la_precedenza_sulla_derivazione(tmp_path):
    """
    **Obiettivo**: Verificare che quando ``batch_table`` dichiara esplicitamente
    il modulo (``"Node 3"``), questo prevalga sulla derivazione via regex (``"NOD3"``).

    **Razionale Scientifico/Sistemistico**: Il dato curato nella ``batch_table``
    ha priorità sulla ricostruzione euristica per prefisso (e consente di
    recuperare moduli come ``"Airlock"`` il cui codice ``A/L1`` non rispetta il
    pattern a 3 lettere).
    """
    campioni = [
        Campione("ERX1000001", "QUALUNQUE", posizione="NOD3O2",
                 modulo_arricchimento="Node 3")
    ]
    scenario = crea_scenario(tmp_path, campioni, con_arricchimento=True)
    assert esegui_gate_metadati(scenario.config)["ERX1000001"].modulo == "Node 3"


def test_senza_arricchimento_il_modulo_si_deriva(tmp_path):
    """
    **Obiettivo**: Verificare che senza ``batch_table`` il modulo venga
    derivato dal prefisso della posizione (``"NOD3O2" -> "NOD3"``).

    **Razionale Scientifico/Sistemistico**: Garantisce l'autonomia della
    pipeline sui soli file ISA-Tab ufficiali di GeneLab.
    """
    campioni = [Campione("ERX1000001", "QUALUNQUE", posizione="NOD3O2")]
    scenario = crea_scenario(tmp_path, campioni, con_arricchimento=False)
    assert esegui_gate_metadati(scenario.config)["ERX1000001"].modulo == "NOD3"


def test_un_modulo_dichiarato_non_applicabile_resta_assente(tmp_path):
    """
    **Obiettivo**: Verificare che se la ``batch_table`` dichiara ``"not applicable"``
    per il modulo, il valore finale resti ``None`` senza cadere nel fallback regex.

    **Razionale Scientifico/Sistemistico**: Un ``"not applicable"`` esplicito nel
    file di arricchimento è un'asserzione negativa deliberata che non deve
    essere sovrascritta dall'espressione regolare.
    """
    campioni = [
        Campione("ERX1000001", "QUALUNQUE", posizione="NOD3O2",
                 modulo_arricchimento="not applicable")
    ]
    scenario = crea_scenario(tmp_path, campioni, con_arricchimento=True)
    assert esegui_gate_metadati(scenario.config)["ERX1000001"].modulo is None


def test_senza_la_colonna_del_modulo_si_ricade_sulla_derivazione(tmp_path):
    """
    **Obiettivo**: Verificare che se la ``batch_table`` contiene le colonne di
    piastra e corsa ma non la colonna ``module``, il modulo venga comunque
    derivato dalla posizione nella Study Table.

    **Razionale Scientifico/Sistemistico**: Rende indipendenti l'arricchimento
    tecnico dei lotti di sequenziamento e la stratificazione spaziale dei campioni.
    """
    campioni = [Campione("ERX1000001", "QUALUNQUE", posizione="NOD3O2")]
    scenario = crea_scenario(
        tmp_path, campioni, con_arricchimento=True, colonna_modulo_arricchimento=None
    )
    assert esegui_gate_metadati(scenario.config)["ERX1000001"].modulo == "NOD3"


@pytest.mark.parametrize("posizione", ["Air Sample", "Unopened 3DMM Swab Tube", "Not Applicable", ""])
def test_una_posizione_non_di_superficie_non_ha_modulo_in_nessuna_modalita(
    posizione, tmp_path
):
    """
    **Obiettivo**: Verificare che le posizioni non di superficie (``"Air Sample"``,
    ``"Unopened 3DMM Swab Tube"``, ``"Not Applicable"``, ``""``) ricevano sempre
    ``modulo = None`` sia con ``batch_table`` attiva sia senza ``batch_table``,
    anche quando la ``batch_table`` riporta ``"Node 1"``.

    **Razionale Scientifico/Sistemistico**: Nella ``plate_well_map_960.tsv`` di
    OSD-734 i filtri d'aria e i tamponi non aperti portano talvolta nella colonna
    ``module`` il nodo in cui sono stati raccolti. Se questa regola si applicasse
    solo alla derivazione regex e non alla ``batch_table``, lo stesso identico
    campione d'aria verrebbe incluso nelle analisi di superficie quando si usa
    ``batch_table`` ed escluso quando non la si usa, producendo due risultati
    biologici discordanti sullo stesso dataset.
    """
    campioni = [
        Campione("ERX1000001", "NOD1.F3", posizione=posizione,
                 modulo_arricchimento="Node 1")
    ]

    con = crea_scenario(tmp_path / "con", campioni, con_arricchimento=True)
    senza = crea_scenario(tmp_path / "senza", campioni, con_arricchimento=False)

    assert esegui_gate_metadati(con.config)["ERX1000001"].modulo is None
    assert esegui_gate_metadati(senza.config)["ERX1000001"].modulo is None


def test_una_posizione_di_superficie_ha_il_modulo_in_entrambe(tmp_path):
    """
    **Obiettivo**: Verificare che un vero campione di superficie (``"NOD1D4"``)
    riceva un modulo valido in entrambe le modalità (``"Node 1"`` da ``batch_table``,
    ``"NOD1"`` per derivazione regex).

    **Razionale Scientifico/Sistemistico**: Conferma la simmetria del filtro
    sulle posizioni non di superficie: non sottrae alcun campione di superficie
    legittimo dalla stratificazione spaziale.
    """
    campioni = [
        Campione("ERX1000001", "QUALUNQUE", posizione="NOD1D4",
                 modulo_arricchimento="Node 1")
    ]

    con = crea_scenario(tmp_path / "con", campioni, con_arricchimento=True)
    senza = crea_scenario(tmp_path / "senza", campioni, con_arricchimento=False)

    assert esegui_gate_metadati(con.config)["ERX1000001"].modulo == "Node 1"
    assert esegui_gate_metadati(senza.config)["ERX1000001"].modulo == "NOD1"


def test_le_posizioni_non_di_superficie_sono_configurabili(tmp_path):
    """
    **Obiettivo**: Verificare che l'elenco delle posizioni non di superficie sia
    governato dal parametro ``meta.non_surface_positions`` e non cablato nel codice.

    **Razionale Scientifico/Sistemistico**: Garantisce la portabilità del
    crosswalk verso altri studi ambientali senza dover modificare il codice Python.
    """
    campioni = [
        Campione("ERX1000001", "QUALUNQUE", posizione="Vuoto Spaziale",
                 modulo_arricchimento="Node 1")
    ]
    scenario = crea_scenario(
        tmp_path, campioni, con_arricchimento=True,
        sovrascrivi={"meta": {"non_surface_positions": ["Vuoto Spaziale"]}},
    )
    assert esegui_gate_metadati(scenario.config)["ERX1000001"].modulo is None


def test_un_campione_senza_riga_di_lotto_usa_la_derivazione(tmp_path):
    """
    **Obiettivo**: Verificare che se un singolo campione non è presente nella
    ``batch_table``, il suo modulo venga comunque ricavato tramite ``meta.module_regex``.

    **Razionale Scientifico/Sistemistico**: Garantisce una degradazione morbida
    per-campione anche in presenza di una tabella di arricchimento incompleta.
    """
    campioni = [Campione("ERX1000001", "QUALUNQUE", posizione="NOD3O2",
                         chiave_arricchimento="ERX9999999")]
    scenario = crea_scenario(tmp_path, campioni, con_arricchimento=True)
    assert esegui_gate_metadati(scenario.config)["ERX1000001"].modulo == "NOD3"
