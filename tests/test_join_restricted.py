"""Suite di verifica del Gate G03 e della restrizione del join all'Assay Table.

Inquadramento nel Piano Operativo
---------------------------------
* **Settimana di riferimento**: **Settimana 6 (W6 — Fase F2: Crosswalk e
  metadati, Gate G03: join ristretto)**.
* **Moduli sorgente coperti**:
  - ``src/amplicon16s/metadata/crosswalk.py``
  - ``src/amplicon16s/gates/g01_g15.py`` (Gate G03)

Scopo sperimentale e razionale scientifico/sistemistico
-------------------------------------------------------
Negli studi multi-omici depositati su **NASA GeneLab** (come **OSD-734**), la
*Study Sample Table* ISA-Tab (``s_OSD-734.txt``) è condivisa fra più saggi dello
stesso studio (es. sequenziamento amplicone 16S rRNA, metagenomica shotgun WGS e
spettrometria di massa con campioni di tipo ``"solvent control"``) e contiene
quindi righe che non appartengono all'assay 16S:

1. **Restrizione per costruzione all'Assay Table**: il crosswalk itera
   esclusivamente sulle righe dell'*Assay Table 16S* e consulta la *Study Table*
   come dizionario di lookup. In questo modo i campioni estranei di altri saggi
   non entrano mai nell'inventario 16S e non fanno fallire il Gate G11 per
   materiali sconosciuti (come ``"solvent control"``).
2. **Sorveglianza del Gate G03 (``E-S0-03``)**: G03 presidia gli unici due casi
   in cui l'integrità relazionale del join ristretto può rompersi:
   - un nome di campione dell'Assay Table che compare **più volte** nella Study
     Table (ambiguità 1:N che moltiplicherebbe le righe dell'assay o
     attribuirebbe metadati sbagliati);
   - un campione dell'Assay Table **assente** nella Study Table (che resterebbe
     privo di classe biologica/controllo e di posizione).
"""

from __future__ import annotations

import pytest
from conftest import BIOLOGICO, NEGATIVO, POSITIVO, Campione, crea_scenario

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

    **Razionale Scientifico/Sistemistico**: Impedisce che campioni di
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

    **Razionale Scientifico/Sistemistico**: Se la classificazione dei ruoli (G11)
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

    **Razionale Scientifico/Sistemistico**: Certifica formalmente che l'insieme
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

    **Razionale Scientifico/Sistemistico**: Una chiave duplicata nella Study
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

    **Razionale Scientifico/Sistemistico**: Senza la riga corrispondente nella
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
    assert "NOD1D4.L2" in str(errore.value)
    assert "senza classe" in str(errore.value)


def test_la_diagnosi_distingue_ambiguita_e_assenza(tmp_path):
    """
    **Obiettivo**: Verificare che la struttura diagnostica di ``analizza``
    separi in campi distinti ``nomi_ambigui_nello_studio``,
    ``nomi_assenti_nello_studio`` e ``righe_studio_estranee``.

    **Razionale Scientifico/Sistemistico**: Fornisce all'operatore un referto
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

    **Razionale Scientifico/Sistemistico**: Conferma che G03 non segnala
    anomalie interne ad altri saggi dello studio che non toccano alcuno dei
    campioni dell'amplicone 16S.
    """
    campioni = _campioni_dell_assay()
    extra = ALTRI_ASSAY + [("MS.SOLV.01", "solvent control", "Not Applicable")]
    scenario = crea_scenario(tmp_path, campioni, righe_studio_extra=extra)

    inventario = esegui_gate_metadati(scenario.config)
    assert len(inventario) == 4
