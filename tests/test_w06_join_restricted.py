r"""Suite di verifica del Gate G03 e della restrizione del join all'Assay Table.

1. Inquadramento nel Piano Operativo
------------------------------------
Settimana 6 (W6), Fase F2 (crosswalk e metadati, Gate G03: join ristretto).

2. Moduli sorgente coperti
--------------------------
* ``src/amplicon16s/metadata/crosswalk.py``
* ``src/amplicon16s/gates/g01_g15.py`` (Gate G03)

3. Cosa valuta questo file
--------------------------
- restrizione per costruzione alla tabella di assay: le righe di altri saggi
  presenti nella tabella di studio non entrano nell'inventario, e il loro
  materiale non fa fallire G11; l'inventario ha la dimensione dell'assay e non
  dello studio;
- Gate G03 (``E-S0-03``): fallimento quando un nome della tabella di assay
  compare più volte nella tabella di studio o vi manca, con una diagnosi che
  distingue ambiguità e assenza;
- un nome ripetuto nella tabella di studio ma estraneo all'assay non disturba.

4. Comandi Bash e scenari di esecuzione
---------------------------------------
    ``<immagine>`` e' l'immagine del container della pipeline; quella corrente
    e' indicata in ``test.txt``, sezione 1.3.

    1. Modalità locale standard (R di base con jsonlite, senza Bioconductor né
       dati reali):
       pytest tests/test_w06_join_restricted.py -v

    2. Modalità container Docker standard (sottoinsieme ridotto con
       R/Bioconductor):
       docker run --rm \
         -e PYTHONPATH=/app/src \
         -e AMPLICON16S_R_DIR=/app/R \
         -v "$(pwd)":/app \
         -w /app \
         <immagine> \
         pytest -o cache_dir=/tmp/.pytest_cache tests/test_w06_join_restricted.py -v

    3. Modalità container Docker completa (con i 2.4 GB di dati reali OSD-734;
       la configurazione e i percorsi che contiene devono stare nella cartella
       montata):
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
         pytest -o cache_dir=/tmp/.pytest_cache tests/test_w06_join_restricted.py -v

5. Risultato atteso
-------------------
7 test totali:
- 7 passed in ambiente locale standard (~0.1s);
- 7 passed nel container Docker standard sul sottoinsieme ridotto (~0.1s);
- 7 passed nel container Docker con i dati reali OSD-734 (~0.1s).

6. Razionale scientifico e sistemistico
---------------------------------------
Negli studi multi-omici depositati su NASA GeneLab, come OSD-734, la tabella
dei campioni di studio (``s_OSD-734.txt``) è condivisa fra più saggi
(ampliconi 16S, metagenomica shotgun, metabolomica con campioni di tipo
``"solvent control"``) e contiene righe che non appartengono all'assay 16S. Il
crosswalk itera solo sulla tabella di assay e usa quella di studio come
dizionario: così i campioni di altri saggi non entrano mai nell'inventario. G03
sorveglia i due soli casi in cui il join ristretto può rompersi: un nome
ambiguo, che attribuirebbe metadati sbagliati, e un nome assente, che
lascerebbe il campione senza classe né posizione.
"""

from __future__ import annotations

import pytest
from conftest import BIOLOGICO, NEGATIVO, POSITIVO, Campione, crea_scenario

from amplicon16s.errors.catalog import Categoria
from amplicon16s.gates.g01_g15 import ErroreGate, esegui_gate_metadati
from amplicon16s.metadata.crosswalk import analizza

#: Righe che appartengono ad altri assay dello stesso studio: nomi che
#: nell'assay dell'amplicone non compaiono, e un materiale che le tre categorie
#: dichiarate non prevedono.
ALTRI_ASSAY = [
    ("MS.SOLV.01", "solvent control", "Not Applicable"),
    ("MS.SOLV.02", "solvent control", "Not Applicable"),
    ("WGS.SWAB.01", "Surface swab", "NOD2D1"),
]


def _campioni_dell_assay() -> list[Campione]:
    """Quattro campioni della tabella di assay: due biologici, un positivo, un negativo.
    """
    return [
        Campione("ERX2000001", "NOD1D4.L1"),
        Campione("ERX2000002", "NOD1D4.L2"),
        Campione("ERX2000003", "POS.P1.1", materiale=POSITIVO, posizione="Not Applicable"),
        Campione("ERX2000004", "BLANK.P1.1", materiale=NEGATIVO, posizione="Not Applicable"),
    ]


# --------------------------------------------------------------------------- #
# Il caso che passa: la tabella di studio è più grande, l'inventario no        #
# --------------------------------------------------------------------------- #


def test_le_righe_di_altri_assay_non_entrano_nell_inventario(tmp_path):
    """
    **Obiettivo**: Verificare che aggiungendo 3 righe di altri saggi
    (``ALTRI_ASSAY``) alla Study Table, l'inventario restituito da
    ``esegui_gate_metadati`` contenga esclusivamente i 4 campioni dell'Assay Table.

    **Razionale scientifico e sistemistico**: Impedisce che campioni di
    spettrometria di massa o metagenomica shotgun presenti nella Study Table
    condivisa di NASA GeneLab (OSD-734) vengano inclusi nell'inventario 16S,
    evitando falsi allarmi per file FASTQ mancanti o alterazioni del
    denominatore di prevalenza.
    """
    campioni = _campioni_dell_assay()
    scenario = crea_scenario(tmp_path, campioni, righe_studio_extra=ALTRI_ASSAY)

    inventario = esegui_gate_metadati(scenario.config)

    assert len(inventario) == len(campioni) == 4
    nomi = {c.nome for c in inventario}
    assert nomi == {c.nome for c in campioni}
    assert not nomi & {nome for nome, _, _ in ALTRI_ASSAY}


def test_il_materiale_di_un_altro_assay_non_fa_fallire_g11(tmp_path):
    """
    **Obiettivo**: Verificare che la presenza del materiale ``"solvent control"``
    nelle righe extra della Study Table non produca ``materiali_non_mappati`` in
    ``analizza`` e venga contabilizzata in ``righe_studio_estranee == 3``.

    **Razionale scientifico e sistemistico**: Se la classificazione dei ruoli (G11)
    venisse applicata all'intera Study Table prima di restringere il dominio
    all'Assay Table 16S, la pipeline fallirebbe segnalando ``"solvent control"``
    come materiale non riconosciuto, generando un falso positivo bloccante su
    campioni che non appartengono al sequenziamento 16S.
    """
    scenario = crea_scenario(
        tmp_path, _campioni_dell_assay(), righe_studio_extra=ALTRI_ASSAY
    )
    analisi = analizza(scenario.config)

    assert analisi.materiali_non_mappati == {}
    assert analisi.righe_studio_estranee == len(ALTRI_ASSAY)


def test_l_inventario_ha_la_dimensione_dell_assay_non_dello_studio(tmp_path):
    """
    **Obiettivo**: Verificare che con 7 righe nella Study Table e 4 nell'Assay
    Table la cardinalità dell'inventario sia esattamente 4.

    **Razionale scientifico e sistemistico**: Certifica formalmente che l'insieme
    universo dello studio 16S è definito dall'Assay Table (960 campioni in
    OSD-734) e mai dalla Study Sample Table.
    """
    campioni = _campioni_dell_assay()
    scenario = crea_scenario(tmp_path, campioni, righe_studio_extra=ALTRI_ASSAY)

    inventario = esegui_gate_metadati(scenario.config)

    righe_di_studio = len(campioni) + len(ALTRI_ASSAY)
    assert righe_di_studio == 7
    assert len(inventario) == 4, "l'insieme di riferimento è la tabella di assay"


# --------------------------------------------------------------------------- #
# Il caso che fallisce: un nome riusato da un altro assay                      #
# --------------------------------------------------------------------------- #


def test_g03_fallisce_se_un_nome_compare_due_volte_nello_studio(tmp_path):
    """
    **Obiettivo**: Verificare che se il campione ``NOD1D4.L1`` dell'Assay Table
    compare 2 volte nella Study Table, ``esegui_gate_metadati`` sollevi
    ``ErroreGate`` sul Gate ``G03`` con codice ``E-S0-03``.

    **Razionale scientifico e sistemistico**: Una chiave duplicata nella Study
    Table per un campione dell'assay 16S causerebbe un prodotto cartesiano (1:N)
    durante il join o l'attribuzione arbitraria dei metadati di una delle due
    righe omonime.
    """
    campioni = _campioni_dell_assay()
    scenario = crea_scenario(
        tmp_path, campioni, righe_studio_ripetute=["NOD1D4.L1"]
    )

    with pytest.raises(ErroreGate) as errore:
        esegui_gate_metadati(scenario.config)

    assert errore.value.gate == "G03"
    assert errore.value.codice == "E-S0-03"
    assert "NOD1D4.L1" in str(errore.value)
    assert "2 righe" in str(errore.value)


def test_g03_fallisce_se_un_campione_dell_assay_manca_nello_studio(tmp_path):
    """
    **Obiettivo**: Verificare che rimuovendo ``NOD1D4.L2`` da ``studio.txt`` il
    Gate G03 fallisca con codice ``E-S0-03`` segnalando il campione ``senza classe``.

    **Razionale scientifico e sistemistico**: Senza la riga corrispondente nella
    Study Table, un campione presente nell'Assay Table resterebbe privo di
    ``Characteristics[Material Type]`` e ``Factor Value``, rendendo impossibile
    stabilirne il ruolo (biologico, blank o mock) e il modulo spaziale ISS.
    """
    campioni = _campioni_dell_assay()
    scenario = crea_scenario(tmp_path, campioni)

    # Si riscrive la tabella di studio senza uno dei campioni dell'assay.
    righe = (scenario.radice / "studio.txt").read_text(encoding="utf-8").splitlines()
    tenute = [r for r in righe if not r.startswith("NOD1D4.L2\t")]
    (scenario.radice / "studio.txt").write_text("\n".join(tenute) + "\n", encoding="utf-8")

    with pytest.raises(ErroreGate) as errore:
        esegui_gate_metadati(scenario.config)

    assert errore.value.codice == "E-S0-03"
    assert errore.value.categoria is Categoria.REVISIONE_UMANA
    assert "NOD1D4.L2" in str(errore.value)
    assert "senza classe" in str(errore.value)


def test_la_diagnosi_distingue_ambiguita_e_assenza(tmp_path):
    """
    **Obiettivo**: Verificare che la struttura diagnostica di ``analizza``
    separi in campi distinti ``nomi_ambigui_nello_studio``,
    ``nomi_assenti_nello_studio`` e ``righe_studio_estranee``.

    **Razionale scientifico e sistemistico**: Fornisce all'operatore un referto
    diagnostico non ambiguo che distingue immediatamente se un fallimento di G03
    dipende da una riga duplicata o da una riga mancante nella Study Table.
    """
    campioni = _campioni_dell_assay()
    scenario = crea_scenario(
        tmp_path, campioni, righe_studio_ripetute=["POS.P1.1"], righe_studio_extra=ALTRI_ASSAY
    )
    analisi = analizza(scenario.config)

    assert analisi.nomi_ambigui_nello_studio == {"POS.P1.1": 2}
    assert analisi.nomi_assenti_nello_studio == []
    assert analisi.righe_studio_estranee == len(ALTRI_ASSAY)


def test_un_nome_riusato_ma_estraneo_all_assay_non_disturba(tmp_path):
    """
    **Obiettivo**: Verificare che se il nome duplicato nella Study Table
    appartiene a un altro assay (``MS.SOLV.01``) e non all'Assay Table 16S, il
    Gate G03 passi regolarmente restituendo i 4 campioni dell'assay.

    **Razionale scientifico e sistemistico**: Conferma che G03 non segnala
    anomalie interne ad altri saggi dello studio che non toccano alcuno dei
    campioni dell'amplicone 16S.
    """
    campioni = _campioni_dell_assay()
    extra = ALTRI_ASSAY + [("MS.SOLV.01", "solvent control", "Not Applicable")]
    scenario = crea_scenario(tmp_path, campioni, righe_studio_extra=extra)

    inventario = esegui_gate_metadati(scenario.config)
    assert len(inventario) == 4
