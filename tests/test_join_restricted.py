"""Test del gate G03: il join resta ristretto alla tabella di assay.

La tabella campioni di studio è condivisa fra più assay dello stesso studio e
contiene righe che all'assay dell'amplicone non appartengono. Se l'insieme dei
campioni venisse preso da lì, la pipeline analizzerebbe campioni che non ha —
e, nel caso peggiore, li classificherebbe con valori che in questo assay non
esistono nemmeno.

La restrizione è garantita per costruzione: si itera sulle righe dell'assay e
la tabella di studio viene solo consultata. Resta un modo in cui può rompersi,
ed è quello che G03 sorveglia: un nome di campione che nella tabella di studio
compare più volte, perché riusato da un altro assay.
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
    campioni = _campioni_dell_assay()
    scenario = crea_scenario(tmp_path, campioni, righe_studio_extra=ALTRI_ASSAY)

    inventario = esegui_gate_metadati(scenario.config)

    assert len(inventario) == len(campioni) == 4
    nomi = {c.nome for c in inventario}
    assert nomi == {c.nome for c in campioni}
    assert not nomi & {nome for nome, _, _ in ALTRI_ASSAY}


def test_il_materiale_di_un_altro_assay_non_fa_fallire_g11(tmp_path):
    """"solvent control" non è fra le categorie dichiarate.

    Se il join non fosse ristretto, quelle righe entrerebbero e G11
    fallirebbe segnalando un materiale non mappato che a questo assay non
    appartiene: un errore vero, ma con la causa sbagliata.
    """
    scenario = crea_scenario(
        tmp_path, _campioni_dell_assay(), righe_studio_extra=ALTRI_ASSAY
    )
    analisi = analizza(scenario.config)

    assert analisi.materiali_non_mappati == {}
    assert analisi.righe_studio_estranee == len(ALTRI_ASSAY)


def test_l_inventario_ha_la_dimensione_dell_assay_non_dello_studio(tmp_path):
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
    """Senza riga di studio il campione resterebbe senza classe."""
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
    campioni = _campioni_dell_assay()
    scenario = crea_scenario(
        tmp_path, campioni, righe_studio_ripetute=["POS.P1.1"], righe_studio_extra=ALTRI_ASSAY
    )
    analisi = analizza(scenario.config)

    assert analisi.nomi_ambigui_nello_studio == {"POS.P1.1": 2}
    assert analisi.nomi_assenti_nello_studio == []
    assert analisi.righe_studio_estranee == len(ALTRI_ASSAY)


def test_un_nome_riusato_ma_estraneo_all_assay_non_disturba(tmp_path):
    """La ripetizione conta solo se riguarda un campione di questo assay."""
    campioni = _campioni_dell_assay()
    extra = ALTRI_ASSAY + [("MS.SOLV.01", "solvent control", "Not Applicable")]
    scenario = crea_scenario(tmp_path, campioni, righe_studio_extra=extra)

    inventario = esegui_gate_metadati(scenario.config)
    assert len(inventario) == 4
