"""Test del crosswalk e dei gate G04, G05, G06 e G11, su fixture sintetiche.

Ogni gate ha un caso che lo supera e uno che lo fa fallire con il codice
atteso. Il dato reale è corretto per costruzione, quindi i casi di fallimento
possono esistere solo qui.
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
    return [
        Campione("ERX1000001", "NOD1D4.L1"),
        Campione("ERX1000002", "POS.P1.1", materiale=POSITIVO, posizione="Not Applicable"),
        Campione("ERX1000003", "BLANK.P1.1", materiale=NEGATIVO, posizione="Not Applicable"),
    ]


# --------------------------------------------------------------------------- #
# Estrazione dell'accession                                                    #
# --------------------------------------------------------------------------- #


def test_estrae_l_accession_anche_con_gruppi_di_cattura():
    """findall restituirebbe il gruppo, non la corrispondenza intera."""
    assert estrai_accession("GLDS-653_Amplicon_ERX12083297_raw.fastq.gz", ACCESSION) == (
        "ERX12083297"
    )


def test_nessun_accession_e_un_errore():
    with pytest.raises(AccessionNonEstraibile, match="nessun accession"):
        estrai_accession("campione_senza_codice.fastq.gz", ACCESSION)


def test_due_accession_sono_un_errore_quanto_nessuno():
    """Prendere il primo sarebbe la scelta comoda e silenziosamente sbagliata."""
    with pytest.raises(AccessionNonEstraibile, match="ambiguo"):
        estrai_accession("ERX1000001_ERX1000002.fastq.gz", ACCESSION)


# --------------------------------------------------------------------------- #
# Lo scenario corretto passa tutti i gate                                      #
# --------------------------------------------------------------------------- #


def test_scenario_corretto_produce_l_inventario(tmp_path):
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
    """Un risultato deve restare riconducibile al dato che l'ha prodotto."""
    scenario = crea_scenario(tmp_path, _tre_campioni())
    inventario = esegui_gate_metadati(scenario.config)
    assert inventario["ERX1000002"].materiale == POSITIVO


def test_ogni_campione_e_legato_al_proprio_file(tmp_path):
    scenario = crea_scenario(tmp_path, _tre_campioni())
    for campione in esegui_gate_metadati(scenario.config):
        assert campione.file is not None
        assert campione.accession in campione.file.name


# --------------------------------------------------------------------------- #
# G04 — l'accession è estraibile e non ambiguo                                 #
# --------------------------------------------------------------------------- #


def test_g04_passa_sullo_scenario_corretto(tmp_path):
    scenario = crea_scenario(tmp_path, _tre_campioni())
    assert len(esegui_gate_metadati(scenario.config)) == 3


def test_g04_fallisce_se_un_nome_di_file_non_contiene_accession(tmp_path):
    campioni = _tre_campioni()
    campioni[0].file = "campione_senza_codice.fastq.gz"
    scenario = crea_scenario(tmp_path, campioni)

    with pytest.raises(ErroreGate) as errore:
        esegui_gate_metadati(scenario.config)

    assert errore.value.gate == "G04"
    assert errore.value.codice == "E-S0-04"
    assert "campione_senza_codice.fastq.gz" in str(errore.value)


def test_g04_fallisce_se_l_accession_e_ambiguo(tmp_path):
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
    scenario = crea_scenario(
        tmp_path, _tre_campioni(), file_in_piu=["ERX1000001_copia.fastq.gz"]
    )

    with pytest.raises(ErroreGate) as errore:
        esegui_gate_metadati(scenario.config)

    assert errore.value.gate == "G05"
    assert errore.value.codice == "E-S0-05"
    assert "ERX1000001" in str(errore.value)


def test_g05_fallisce_su_accession_duplicati_nell_assay(tmp_path):
    campioni = _tre_campioni()
    # Due campioni diversi che dichiarano lo stesso file grezzo.
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
    campioni = _tre_campioni()
    campioni[0].materiale = ""
    scenario = crea_scenario(tmp_path, campioni)

    with pytest.raises(ErroreGate, match="E-S0-11"):
        esegui_gate_metadati(scenario.config)


def test_la_classificazione_tollera_maiuscole_e_spazi(tmp_path):
    """Le tabelle di metadati sono compilate a mano."""
    campioni = _tre_campioni()
    campioni[1].materiale = "  positive control  "
    scenario = crea_scenario(tmp_path, campioni)
    inventario = esegui_gate_metadati(scenario.config)
    assert inventario["ERX1000002"].classe is ClasseCampione.CONTROLLO_POSITIVO


def test_la_mappa_non_inventa_una_categoria_di_ripiego(tmp_path):
    scenario = crea_scenario(tmp_path, _tre_campioni())
    mappa = MappaControlli.da_configurazione(scenario.config.ctrl)
    assert mappa.classifica("materiale sconosciuto") is None
    assert mappa.classifica(None) is None


# --------------------------------------------------------------------------- #
# La chiave è l'accession, non il nome del campione                            #
# --------------------------------------------------------------------------- #


def test_il_nome_del_campione_si_ripete_fra_repliche_l_accession_no(tmp_path):
    """Usare il nome come chiave appaierebbe le due repliche fra loro.

    Due repliche dello stesso campione biologico condividono il nome e hanno
    accession diversi. Chi usasse il nome come chiave troverebbe due file per
    una chiave sola e dovrebbe sceglierne uno: metà delle volte sbaglierebbe, e
    il risultato resterebbe plausibile.
    """
    campioni = [
        Campione("ERX1000001", "LAB1P3.L1", posizione="LAB1P3"),
        Campione("ERX1000002", "LAB1P3.L1", posizione="LAB1P3"),
    ]
    scenario = crea_scenario(tmp_path, campioni)

    # Con l'accession come chiave i due restano distinti e il gate passa:
    # l'assay dichiara due righe con lo stesso nome ma accession diversi.
    analisi = analizza(scenario.config)
    assert analisi.accession_ripetuti_nell_assay == {}

    # Con il nome come chiave i due collassano in uno.
    per_nome = {c.nome for c in campioni}
    per_accession = {c.accession for c in campioni}
    assert len(per_nome) == 1, "il nome non distingue le repliche"
    assert len(per_accession) == 2, "l'accession le distingue"


def test_appaiare_per_nome_produce_una_corrispondenza_errata(tmp_path):
    """Dimostrazione concreta dell'appaiamento sbagliato.

    Lo scenario ha due repliche con lo stesso nome e materiali diversi: la
    prima biologica, la seconda un controllo. Appaiando per nome, al file
    della seconda verrebbe attribuita la classe della prima.
    """
    campioni = [
        Campione("ERX1000001", "RIPETUTO", materiale=BIOLOGICO),
        Campione("ERX1000002", "RIPETUTO", materiale=POSITIVO, file="ERX1000002_b.fastq.gz"),
    ]
    scenario = crea_scenario(tmp_path, campioni)

    # Appaiamento per nome: una sola chiave, l'ultima scrittura vince.
    per_nome = {}
    for campione in campioni:
        per_nome[campione.nome] = campione.materiale
    assert len(per_nome) == 1
    assert per_nome["RIPETUTO"] == POSITIVO  # il primo campione ha perso la sua classe

    # Appaiamento per accession: ciascuno conserva la propria.
    per_accession = {c.accession: c.materiale for c in campioni}
    assert per_accession["ERX1000001"] == BIOLOGICO
    assert per_accession["ERX1000002"] == POSITIVO


# --------------------------------------------------------------------------- #
# Derivazione del modulo                                                       #
# --------------------------------------------------------------------------- #


def test_il_modulo_si_deriva_dalla_posizione(tmp_path):
    campioni = [Campione("ERX1000001", "QUALUNQUE", posizione="NOD3O2 (ARED)")]
    scenario = crea_scenario(tmp_path, campioni)
    assert esegui_gate_metadati(scenario.config)["ERX1000001"].modulo == "NOD3"


@pytest.mark.parametrize("posizione", ["Air Sample", "Not Applicable", ""])
def test_le_posizioni_che_non_sono_superfici_non_hanno_modulo(posizione, tmp_path):
    """Non vanno forzate in un modulo: restano una categoria a parte."""
    campioni = [Campione("ERX1000001", "NOD1.F3", posizione=posizione)]
    scenario = crea_scenario(tmp_path, campioni)
    assert esegui_gate_metadati(scenario.config)["ERX1000001"].modulo is None


def test_il_modulo_non_si_deriva_dal_nome(tmp_path):
    """Un campione d'aria puo' chiamarsi come una superficie."""
    campioni = [Campione("ERX1000001", "NOD1.F3", posizione="Air Sample")]
    scenario = crea_scenario(tmp_path, campioni)
    campione = esegui_gate_metadati(scenario.config)["ERX1000001"]
    assert campione.nome.startswith("NOD1")
    assert campione.modulo is None


def test_la_derivazione_si_puo_disattivare(tmp_path):
    scenario = crea_scenario(
        tmp_path, _tre_campioni(), sovrascrivi={"meta": {"derive_module": False}}
    )
    assert all(c.modulo is None for c in esegui_gate_metadati(scenario.config))


# --------------------------------------------------------------------------- #
# Arricchimento del lotto: facoltativo                                         #
# --------------------------------------------------------------------------- #


def test_senza_arricchimento_piastra_e_corsa_restano_nulle(tmp_path):
    scenario = crea_scenario(tmp_path, _tre_campioni(), con_arricchimento=False)
    inventario = esegui_gate_metadati(scenario.config)

    assert len(inventario) == 3
    assert inventario.piastre == ()
    assert inventario.corse == ()
    assert len(inventario.senza_lotto) == 3


def test_con_arricchimento_piastra_e_corsa_sono_valorizzate(tmp_path):
    campioni = _tre_campioni()
    campioni[2].piastra = "2"
    campioni[2].corsa = "corsa_B"
    scenario = crea_scenario(tmp_path, campioni, con_arricchimento=True)
    inventario = esegui_gate_metadati(scenario.config)

    assert inventario.piastre == ("1", "2")
    assert inventario.corse == ("corsa_A", "corsa_B")
    assert inventario.senza_lotto == ()


def test_un_campione_non_coperto_dall_arricchimento_resta_senza_lotto(tmp_path):
    """Copertura parziale: non è un errore, il file è un ingresso in più."""
    campioni = _tre_campioni()
    campioni[0].chiave_arricchimento = "ERX9999999"
    scenario = crea_scenario(tmp_path, campioni, con_arricchimento=True)

    inventario = esegui_gate_metadati(scenario.config)
    assert len(inventario) == 3
    assert len(inventario.senza_lotto) == 1
    assert analizza(scenario.config).senza_riga_di_arricchimento == ["ERX1000001"]


def test_una_chiave_ambigua_nell_arricchimento_non_viene_indovinata(tmp_path):
    """Attribuirne una a caso metterebbe il campione nel lotto sbagliato."""
    campioni = _tre_campioni()
    # Due righe del file dichiarano lo stesso accession.
    campioni[1].chiave_arricchimento = campioni[0].accession
    scenario = crea_scenario(tmp_path, campioni, con_arricchimento=True)

    inventario = esegui_gate_metadati(scenario.config)
    assert inventario["ERX1000001"].piastra is None
    assert analizza(scenario.config).arricchimento_ambiguo == {"ERX1000001": 2}


# --------------------------------------------------------------------------- #
# La chiave dell'arricchimento è l'accession, non il nome                      #
# --------------------------------------------------------------------------- #


def test_l_arricchimento_si_aggancia_per_accession(tmp_path):
    """Due repliche dello stesso campione finiscono nei lotti giusti.

    È la forma che il dato reale ha: la tabella di assay distingue le repliche
    con un suffisso, e i due file stanno su piastre e corse diverse. Un file di
    arricchimento che le identificasse con il nome biologico condiviso
    aggancerebbe una delle due a caso; identificandole per accession ciascuna
    riceve il proprio lotto.
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
    """Il messaggio deve dire cosa fare, non solo che manca."""
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
    """Il dato originale batte la ricostruzione per espressione regolare."""
    campioni = [
        Campione("ERX1000001", "QUALUNQUE", posizione="NOD3O2",
                 modulo_arricchimento="Node 3")
    ]
    scenario = crea_scenario(tmp_path, campioni, con_arricchimento=True)
    assert esegui_gate_metadati(scenario.config)["ERX1000001"].modulo == "Node 3"


def test_senza_arricchimento_il_modulo_si_deriva(tmp_path):
    campioni = [Campione("ERX1000001", "QUALUNQUE", posizione="NOD3O2")]
    scenario = crea_scenario(tmp_path, campioni, con_arricchimento=False)
    assert esegui_gate_metadati(scenario.config)["ERX1000001"].modulo == "NOD3"


def test_un_modulo_dichiarato_non_applicabile_resta_assente(tmp_path):
    """Un "non applicabile" esplicito è un'informazione, non un vuoto da colmare.

    Ricadere sulla regex contraddirebbe il dato dichiarato.
    """
    campioni = [
        Campione("ERX1000001", "QUALUNQUE", posizione="NOD3O2",
                 modulo_arricchimento="not applicable")
    ]
    scenario = crea_scenario(tmp_path, campioni, con_arricchimento=True)
    assert esegui_gate_metadati(scenario.config)["ERX1000001"].modulo is None


def test_senza_la_colonna_del_modulo_si_ricade_sulla_derivazione(tmp_path):
    campioni = [Campione("ERX1000001", "QUALUNQUE", posizione="NOD3O2")]
    scenario = crea_scenario(
        tmp_path, campioni, con_arricchimento=True, colonna_modulo_arricchimento=None
    )
    assert esegui_gate_metadati(scenario.config)["ERX1000001"].modulo == "NOD3"


@pytest.mark.parametrize("posizione", ["Air Sample", "Unopened 3DMM Swab Tube", "Not Applicable", ""])
def test_una_posizione_non_di_superficie_non_ha_modulo_in_nessuna_modalita(
    posizione, tmp_path
):
    """La regola vale in entrambe le vie di attribuzione.

    Se valesse solo per la derivazione, lo stesso campione avrebbe un modulo
    con il file di arricchimento e non l'avrebbe senza: un'analisi raggruppata
    per modulo lo includerebbe o lo escluderebbe a seconda della presenza di un
    ingresso facoltativo, e lo stesso dataset darebbe due risultati.
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
    """Il controcampo: su una superficie entrambe attribuiscono un modulo.

    I nomi differiscono — il file dichiara il nome originale, la derivazione
    estrae il prefisso — ma nessuna delle due lascia il campione senza.
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
    """Il criterio non è scritto nel codice con i valori di questo dataset."""
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
    """Nessun dato dichiarato per quella riga: il ripiego resta disponibile."""
    campioni = [Campione("ERX1000001", "QUALUNQUE", posizione="NOD3O2",
                         chiave_arricchimento="ERX9999999")]
    scenario = crea_scenario(tmp_path, campioni, con_arricchimento=True)
    assert esegui_gate_metadati(scenario.config)["ERX1000001"].modulo == "NOD3"
