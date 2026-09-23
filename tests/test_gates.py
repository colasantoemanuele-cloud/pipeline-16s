"""Test dei quindici gate di validazione e della fase S0.

Ogni gate ha un caso che lo supera e uno che lo fa fallire con il codice
atteso. I casi di fallimento vivono solo qui: il dato reale è corretto e non
li contiene.

Gli scenari scrivono su disco veri file FASTQ compressi, tabelle di metadati e
riferimento tassonomico, così il percorso provato dai test è quello che userà
la pipeline e non una simulazione.
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
    assert len(REGISTRO) == 15
    assert len(set(nomi_dei_gate())) == 15


def test_la_configurazione_si_controlla_per_prima():
    """Un errore di configurazione va scoperto prima di aprire un file."""
    assert nomi_dei_gate()[0] == "G15"


def test_tutti_i_gate_passano_su_uno_scenario_corretto(tmp_path):
    esiti = _esegui(_scenario(tmp_path))
    assert all(e.eseguito and e.superato for e in esiti.values()), _fallito(esiti)


def test_l_esecuzione_si_ferma_al_primo_gate_fallito(tmp_path):
    """Proseguire produrrebbe errori derivati che confondono la diagnosi."""
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
    assert _esegui(_scenario(tmp_path))["G01"].superato


@pytest.mark.parametrize("quale", ["assay.txt", "studio.txt"])
def test_g01_fallisce_su_un_ingresso_mancante(quale, tmp_path):
    scenario = _scenario(tmp_path)
    (scenario.radice / quale).unlink()

    esito = _esegui(scenario)["G01"]
    assert not esito.superato
    assert {v.codice for v in esito.violazioni} == {"E-S0-01"}
    assert quale in str(esito.violazioni[0])


def test_g01_fallisce_se_non_ci_sono_file_di_letture(tmp_path):
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
    assert _esegui(_scenario(tmp_path))["G02"].superato


def test_g02_fallisce_su_una_colonna_assente(tmp_path):
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
    assert _esegui(_scenario(tmp_path))["G07"].superato


def test_g07_fallisce_su_un_file_di_lettura_inversa(tmp_path):
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
    assert _esegui(_scenario(tmp_path))["G08"].superato


def test_g08_fallisce_se_manca_la_colonna_del_lotto(tmp_path):
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
    """Segnala senza bloccare: una composizione inconsueta non è un errore."""
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
    """Il criterio è di metodo, non un confronto con le altre piastre."""
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
    assert _esegui(_scenario(tmp_path))["G09"].superato


def test_g09_fallisce_se_il_troncamento_supera_le_letture(tmp_path):
    scenario = _scenario(tmp_path, sovrascrivi={"filter": {"truncLen": 152}})
    esito = _esegui(scenario)["G09"]

    assert not esito.superato
    assert {v.codice for v in esito.violazioni} == {"E-S0-09"}
    assert "scartate, non accorciate" in str(esito.violazioni[0])


def test_g09_fallisce_su_un_singolo_campione_piu_corto(tmp_path):
    campioni = _campioni()
    campioni[2].lunghezza_letture = 120
    esito = _esegui(_scenario(tmp_path, campioni))["G09"]

    assert not esito.superato
    assert "120 bp" in str(esito.violazioni[0])


def test_g09_avvisa_quando_si_tronca_molto_sotto_il_minimo(tmp_path):
    """Uno scarto ampio non perde campioni, perde basi utili."""
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
    assert _esegui(_scenario(tmp_path))["G10"].superato


def test_g10_fallisce_se_il_primer_e_in_testa(tmp_path):
    campioni = _campioni()
    for campione in campioni:
        campione.inizio_letture = INIZIO_CON_PRIMER
    esito = _esegui(_scenario(tmp_path, campioni))["G10"]

    assert not esito.superato
    codici = {v.codice for v in esito.violazioni}
    assert codici == {"E-S0-10"}
    assert "filter.trimLeft" in str(esito.violazioni[0])


def test_g10_fallisce_se_manca_il_segnale_atteso(tmp_path):
    """Il controllo positivo: l'assenza del primer da sola non basta."""
    campioni = _campioni()
    for campione in campioni:
        campione.inizio_letture = INIZIO_MUTO
    esito = _esegui(_scenario(tmp_path, campioni))["G10"]

    assert not esito.superato
    assert "motivo conservato" in str(esito.violazioni[0])
    assert "qc.min_motif_frac" in str(esito.violazioni[0])


def test_g10_non_pretende_il_segnale_dai_controlli_negativi(tmp_path):
    """Da un bianco non ci si attende il segnale del bersaglio.

    È un principio, non una constatazione: sul dataset di riferimento i
    bianchi il motivo ce l'hanno, per via dei contaminanti che amplificano.
    Qui lo scenario costruisce il caso in cui non ce l'hanno, e il gate deve
    passare lo stesso.
    """
    campioni = _campioni()
    for campione in campioni:
        if campione.materiale == NEGATIVO:
            campione.inizio_letture = INIZIO_MUTO
    esito = _esegui(_scenario(tmp_path, campioni))["G10"]

    assert esito.superato


def test_g10_tollera_un_singolo_campione_biologico_muto(tmp_path):
    """La mediana: un campione mal riuscito non ferma l'intera validazione."""
    campioni = _campioni()
    campioni[0].inizio_letture = INIZIO_MUTO
    esito = _esegui(_scenario(tmp_path, campioni))["G10"]

    assert esito.superato


# --------------------------------------------------------------------------- #
# G12 — il riferimento tassonomico è quello dichiarato                         #
# --------------------------------------------------------------------------- #


def test_g12_passa(tmp_path):
    assert _esegui(_scenario(tmp_path))["G12"].superato


def test_g12_fallisce_su_un_checksum_diverso(tmp_path):
    scenario = _scenario(tmp_path)
    (scenario.radice / "riferimento.fa.gz").write_bytes(b">seq1\nTTTT\n")

    esito = _esegui(scenario)["G12"]
    assert not esito.superato
    assert {v.codice for v in esito.violazioni} == {"E-S0-12"}
    assert "tax.ref_md5" in str(esito.violazioni[0])


def test_g12_fallisce_se_il_riferimento_non_esiste(tmp_path):
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
    assert _esegui(_scenario(tmp_path))["G13"].superato


def test_g13_fallisce_su_un_archivio_non_decomprimibile(tmp_path):
    campioni = _campioni()
    campioni[0].contenuto_grezzo = b"questo non e' un archivio gzip"
    esito = _esegui(_scenario(tmp_path, campioni))["G13"]

    assert not esito.superato
    assert {v.codice for v in esito.violazioni} == {"E-S0-13"}
    assert "archivio non leggibile" in str(esito.violazioni[0])


def test_g13_fallisce_su_un_record_troncato(tmp_path):
    campioni = _campioni()
    campioni[0].contenuto_grezzo = gzip.compress(
        f"@r1\n{lettura()}\n+\n{'I' * 151}\n@r2\n{lettura()}\n".encode()
    )
    esito = _esegui(_scenario(tmp_path, campioni))["G13"]

    assert not esito.superato
    assert "record troncato" in str(esito.violazioni[0])


def test_g13_fallisce_se_la_qualita_ha_lunghezza_diversa(tmp_path):
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
    assert _esegui(_scenario(tmp_path))["G14"].superato


def test_g14_fallisce_se_i_thread_eccedono_le_cpu(tmp_path):
    scenario = _scenario(tmp_path, sovrascrivi={"run": {"threads": 100000}})
    esito = _esegui(scenario)["G14"]

    assert not esito.superato
    assert {v.codice for v in esito.violazioni} == {"E-S0-14"}
    assert "run.threads" in str(esito.violazioni[0])


# --------------------------------------------------------------------------- #
# G15 — coerenza della configurazione                                          #
# --------------------------------------------------------------------------- #


def test_g15_passa(tmp_path):
    assert _esegui(_scenario(tmp_path))["G15"].superato


def test_g15_fallisce_su_una_soglia_degenere(tmp_path):
    scenario = _scenario(tmp_path, sovrascrivi={"decontam": {"threshold": 1.0}})
    esiti = _esegui(scenario)

    assert not esiti["G15"].superato
    assert {v.codice for v in esiti["G15"].violazioni} == {"E-G15-02"}
    assert not esiti["G01"].eseguito, "la configurazione si controlla per prima"


# --------------------------------------------------------------------------- #
# La fase S0                                                                   #
# --------------------------------------------------------------------------- #


def test_s0_scrive_gli_artefatti(tmp_path):
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
    scenario = _scenario(tmp_path)
    risultato = esegui_s0(scenario.config)

    documento = json.loads(
        next(p for p in risultato.artefatti if p.name == "inventario.json").read_text()
    )
    assert documento["campioni"] == 4
    assert documento["denominatore_prevalenza"] == 2  # i soli biologici


def test_s0_registra_il_limite_della_stima_sulle_lunghezze(tmp_path):
    """Il file dice per quali file la statistica copre tutto e per quali no."""
    scenario = _scenario(tmp_path)
    risultato = esegui_s0(scenario.config)

    righe = next(
        p for p in risultato.artefatti if p.name == "letture_ispezionate.tsv"
    ).read_text(encoding="utf-8").splitlines()
    assert righe[0].split("\t")[3] == "file_esaurito"
    assert len(righe) == 5


def test_s0_gli_artefatti_entrano_nel_manifesto(tmp_path):
    from amplicon16s.io_layer.artifacts import AlberoOutput, Fase

    scenario = _scenario(tmp_path)
    esegui_s0(scenario.config)

    albero = AlberoOutput(scenario.config.io.out_root)
    assert albero.fase_completa(Fase.INPUT_VALIDATION)
    assert albero.non_integri(Fase.INPUT_VALIDATION) == ()


def test_s0_solleva_quando_un_gate_fallisce(tmp_path):
    scenario = _scenario(tmp_path, sovrascrivi={"run": {"threads": 100000}})
    with pytest.raises(ErroreGate) as errore:
        esegui_s0(scenario.config)
    assert errore.value.gate == "G14"


def test_s0_puo_riportare_senza_sollevare(tmp_path):
    scenario = _scenario(tmp_path, sovrascrivi={"run": {"threads": 100000}})
    risultato = esegui_s0(scenario.config, solleva=False)

    assert not risultato.superata
    assert [e.gate for e in risultato.falliti] == ["G14"]
    assert any(p.name == "gates.json" for p in risultato.artefatti)


def test_s0_registra_gli_avvisi_senza_fallire(tmp_path):
    scenario = _scenario(
        tmp_path, con_arricchimento=True,
        sovrascrivi={"decontam": {"min_blanks": 3}},
    )
    risultato = esegui_s0(scenario.config)

    assert risultato.superata
    da_g08 = [a for a in risultato.avvisi if a.codice == "E-S0-15"]
    assert len(da_g08) == 1
    assert "controlli negativi" in str(da_g08[0])
