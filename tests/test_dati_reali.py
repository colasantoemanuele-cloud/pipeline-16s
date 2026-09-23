"""Test sul dataset di riferimento OSD-734, se disponibile sul filesystem.

**Sono separati dai test su fixture e non girano da soli.** I dati pesano
2,4 GB e non risiedono nel repository: nella catena di integrazione continua
non ci sono, e la suite deve restare verde lo stesso. Questi test si saltano
quindi in modo pulito quando i dati mancano.

Per eseguirli, imposta ``AMPLICON16S_CONFIG_DATI_REALI`` sul percorso di una
configurazione che punti al dataset::

    AMPLICON16S_CONFIG_DATI_REALI=<percorso/della/configurazione.yaml> \
        pytest -m dati_reali

Quella configurazione è un file dell'utente, non del repository, e vive dove
vivono i dati. Deve essere un'istanza valida dello schema, con ``io.fastq_dir``
sulla cartella delle 960 letture, ``io.assay_table`` sulla tabella di assay
dell'amplicone, ``io.study_table`` sulla tabella campioni di studio e, se si
vuole coprire anche i test sul lotto, ``io.batch_table`` sul file di
arricchimento. ``config/config.example.yaml`` ne è il modello.

I numeri attesi sono quelli accertati del dataset: 960 campioni, 803
biologici, 80 controlli positivi, 77 controlli negativi.
"""

from __future__ import annotations

import csv
import json
import os
from collections import Counter
from pathlib import Path

import pytest
import yaml

from amplicon16s.config.schema import carica, valida
from amplicon16s.gates.g01_g15 import esegui_gate_metadati
from amplicon16s.metadata.crosswalk import analizza
from amplicon16s.metadata.models import ClasseCampione

VARIABILE = "AMPLICON16S_CONFIG_DATI_REALI"

CAMPIONI_ATTESI = 960
BIOLOGICI_ATTESI = 803
POSITIVI_ATTESI = 80
NEGATIVI_ATTESI = 77

pytestmark = pytest.mark.dati_reali


def _percorso_configurazione() -> Path | None:
    grezzo = os.environ.get(VARIABILE)
    if not grezzo:
        return None
    percorso = Path(grezzo).expanduser()
    return percorso if percorso.is_file() else None


def _sorgenti_presenti(dati: dict) -> bool:
    obbligatorie = [dati["io"]["fastq_dir"], dati["io"]["assay_table"], dati["io"]["study_table"]]
    return all(Path(p).expanduser().exists() for p in obbligatorie)


@pytest.fixture(scope="module")
def configurazione():
    percorso = _percorso_configurazione()
    if percorso is None:
        pytest.skip(f"{VARIABILE} non impostata o file assente: dati reali non disponibili")
    dati = yaml.safe_load(percorso.read_text(encoding="utf-8"))
    if not _sorgenti_presenti(dati):
        pytest.skip("le sorgenti indicate dalla configurazione non esistono")
    return carica(percorso)


@pytest.fixture(scope="module")
def inventario(configurazione):
    return esegui_gate_metadati(configurazione)


# --------------------------------------------------------------------------- #
# Crosswalk                                                                    #
# --------------------------------------------------------------------------- #


def test_960_accession_univoci(inventario):
    assert len(inventario) == CAMPIONI_ATTESI
    assert len(set(inventario.accessioni)) == CAMPIONI_ATTESI


def test_nessun_orfano_da_nessuna_delle_due_parti(configurazione):
    analisi = analizza(configurazione)
    assert analisi.solo_nei_file == ()
    assert analisi.solo_nell_assay == ()
    assert analisi.file_senza_accession == []
    assert analisi.assay_senza_accession == []


def test_ogni_campione_ha_il_proprio_file(inventario):
    assert all(c.file is not None and c.file.is_file() for c in inventario)


# --------------------------------------------------------------------------- #
# Inventario e classi                                                          #
# --------------------------------------------------------------------------- #


def test_conteggi_delle_classi(inventario):
    assert inventario.conteggi() == {
        ClasseCampione.BIOLOGICO: BIOLOGICI_ATTESI,
        ClasseCampione.CONTROLLO_POSITIVO: POSITIVI_ATTESI,
        ClasseCampione.CONTROLLO_NEGATIVO: NEGATIVI_ATTESI,
    }


def test_il_join_ristretto_esclude_il_materiale_di_altri_assay(configurazione):
    """La tabella di studio contiene "solvent control", che a questo assay
    non appartiene: se il join non fosse ristretto, G11 fallirebbe."""
    righe_studio = list(
        csv.DictReader(
            open(configurazione.io.study_table, encoding="utf-8"), delimiter="\t"
        )
    )
    materiali = Counter(
        (r.get(configurazione.ctrl.column) or "").strip().strip('"') for r in righe_studio
    )
    assert len(righe_studio) > CAMPIONI_ATTESI, "la tabella di studio è più grande"
    assert materiali.get("solvent control", 0) > 0, "contiene righe di un altro assay"
    assert analizza(configurazione).materiali_non_mappati == {}


# --------------------------------------------------------------------------- #
# Arricchimento del lotto                                                      #
# --------------------------------------------------------------------------- #


def test_dieci_piastre_da_96_e_due_corse(inventario, configurazione):
    if configurazione.io.batch_table is None:
        pytest.skip("la configurazione non indica il file di arricchimento")
    assert len(inventario.piastre) == 10
    assert len(inventario.corse) == 2
    assert set(inventario.campioni_per_piastra().values()) == {96}
    assert inventario.senza_lotto == ()


def test_il_file_di_arricchimento_copre_dieci_piastre_da_96(configurazione):
    if configurazione.io.batch_table is None:
        pytest.skip("la configurazione non indica il file di arricchimento")
    righe = list(
        csv.DictReader(
            open(configurazione.io.batch_table, encoding="utf-8"), delimiter="\t"
        )
    )
    per_piastra = Counter(r[configurazione.decontam.batch_column] for r in righe)
    assert len(per_piastra) == 10
    assert set(per_piastra.values()) == {96}


def test_i_campioni_replicati_ricevono_lotti_distinti(inventario, configurazione):
    """Due campioni sono stati sequenziati due volte.

    Nella tabella di assay compaiono come ``_rep1`` e ``_rep2``, con accession
    distinti. Agganciandoli per accession ciascuna replica riceve il proprio
    lotto; con una chiave per nome le due righe sarebbero indistinguibili e
    una delle due verrebbe attribuita a caso.
    """
    if configurazione.io.batch_table is None:
        pytest.skip("la configurazione non indica il file di arricchimento")

    repliche = {c.nome: c for c in inventario if "_rep" in c.nome}
    assert set(repliche) == {
        "LAB1P3.L1_rep1", "LAB1P3.L1_rep2",
        "NOD2S4.R6_rep1", "NOD2S4.R6_rep2",
    }
    assert all(c.piastra is not None and c.corsa is not None for c in repliche.values())

    # Le due repliche di uno stesso campione stanno su piastre e corse diverse.
    for radice in ("LAB1P3.L1", "NOD2S4.R6"):
        prima, seconda = repliche[f"{radice}_rep1"], repliche[f"{radice}_rep2"]
        assert prima.piastra != seconda.piastra
        assert prima.corsa != seconda.corsa


def test_senza_arricchimento_l_esecuzione_completa_lo_stesso(configurazione):
    dati = configurazione.model_dump(mode="python")
    dati["io"]["batch_table"] = None
    inventario = esegui_gate_metadati(valida(dati))

    assert len(inventario) == CAMPIONI_ATTESI
    assert inventario.piastre == ()
    assert inventario.corse == ()
    assert len(inventario.senza_lotto) == CAMPIONI_ATTESI


# --------------------------------------------------------------------------- #
# Particolarità note del dataset                                               #
# --------------------------------------------------------------------------- #


def test_la_piastra_a_composizione_diversa_non_produce_avvisi(inventario, configurazione):
    """Una piastra ha più campioni biologici e meno controlli negativi.

    È un dato reale dei metadati originali, non un errore: l'inventario si
    costruisce senza segnalazioni, e nessun gate la intercetta.
    """
    if configurazione.io.batch_table is None:
        pytest.skip("la configurazione non indica il file di arricchimento")

    composizione = {}
    for campione in inventario:
        if campione.piastra:
            composizione.setdefault(campione.piastra, Counter())[campione.classe] += 1

    negativi = {p: c[ClasseCampione.CONTROLLO_NEGATIVO] for p, c in composizione.items()}
    anomale = [p for p, n in negativi.items() if n != 8]
    assert len(anomale) == 1, "una sola piastra ha una composizione diversa"

    piastra = anomale[0]
    assert negativi[piastra] < 8
    assert composizione[piastra][ClasseCampione.BIOLOGICO] > 80
    # L'inventario è stato costruito: nessuna eccezione, nessun gate fallito.
    assert len(inventario) == CAMPIONI_ATTESI


def test_con_arricchimento_i_moduli_sono_nove_con_i_nomi_originali(
    inventario, configurazione
):
    """Il dato dichiarato dal file, non una ricostruzione per prefisso."""
    if configurazione.io.batch_table is None:
        pytest.skip("la configurazione non indica il file di arricchimento")
    assert inventario.moduli == (
        "Airlock", "Columbus", "JLP", "JPM",
        "Node 1", "Node 2", "Node 3", "PMM", "US Lab",
    )


def test_senza_arricchimento_i_moduli_sono_otto_e_manca_l_airlock(configurazione):
    """Limite noto della derivazione per prefisso.

    Il prefisso dell'Airlock è ``A/L1``: una lettera, una barra, una lettera,
    una cifra. L'espressione predefinita pretende tre lettere consecutive e
    non lo cattura, quindi i 16 tamponi di superficie dell'Airlock restano
    senza modulo. È un limite dichiarato, non un difetto nascosto: quando il
    file di arricchimento c'è, il modulo viene da lì.
    """
    dati = configurazione.model_dump(mode="python")
    dati["io"]["batch_table"] = None
    inventario = esegui_gate_metadati(valida(dati))

    assert inventario.moduli == (
        "COL1", "JLP1", "JPM1", "LAB1", "NOD1", "NOD2", "NOD3", "PMM1",
    )
    assert "A/L1" not in inventario.moduli

    senza_modulo = [c for c in inventario if c.modulo is None]
    per_classe = Counter(c.classe for c in senza_modulo)
    assert per_classe[ClasseCampione.CONTROLLO_POSITIVO] == POSITIVI_ATTESI
    assert per_classe[ClasseCampione.CONTROLLO_NEGATIVO] == NEGATIVI_ATTESI
    # I 16 tamponi dell'Airlock sono fra i biologici rimasti senza modulo.
    airlock = [c for c in senza_modulo if (c.posizione or "").startswith("A/L")]
    assert len(airlock) == 16


def test_le_due_modalita_differiscono_solo_per_l_airlock(inventario, configurazione):
    """La regola sulle posizioni non di superficie vale in entrambe.

    L'unica differenza deve restare la copertura: il file riconosce l'Airlock,
    la derivazione per prefisso no. Per ogni altro campione le due modalità
    devono concordare sul fatto che un modulo ci sia o non ci sia — i nomi
    differiscono, l'attribuzione no.
    """
    if configurazione.io.batch_table is None:
        pytest.skip("la configurazione non indica il file di arricchimento")

    dati = configurazione.model_dump(mode="python")
    dati["io"]["batch_table"] = None
    derivato = esegui_gate_metadati(valida(dati))

    con_modulo_solo_col_file = [
        c for c in inventario
        if c.modulo is not None and derivato[c.accession].modulo is None
    ]
    con_modulo_solo_senza_file = [
        c for c in derivato
        if c.modulo is not None and inventario[c.accession].modulo is None
    ]

    assert con_modulo_solo_senza_file == []
    assert len(con_modulo_solo_col_file) == 16
    assert all(
        (c.posizione or "").startswith("A/L") for c in con_modulo_solo_col_file
    )

    # Stesso numero di campioni senza modulo, a meno dei 16 dell'Airlock.
    senza_col_file = {c.accession for c in inventario if c.modulo is None}
    senza_derivato = {c.accession for c in derivato if c.modulo is None}
    assert len(senza_col_file) == 223
    assert len(senza_derivato) == 239
    assert senza_col_file < senza_derivato
    assert len(senza_derivato - senza_col_file) == 16


# --------------------------------------------------------------------------- #
# La fase S0 sul dataset completo                                              #
# --------------------------------------------------------------------------- #


@pytest.fixture(scope="module")
def risultato_s0(configurazione, tmp_path_factory):
    """Esegue S0 una volta sola, scrivendo in una cartella temporanea."""
    from amplicon16s.steps.s00_validate import esegui_s0

    dati = configurazione.model_dump(mode="python")
    dati["io"]["out_root"] = str(tmp_path_factory.mktemp("s0"))
    return esegui_s0(valida(dati), solleva=False)


def test_s0_supera_tutti_i_gate_sul_dataset_reale(risultato_s0):
    assert risultato_s0.superata, [e.gate for e in risultato_s0.falliti]
    assert len(risultato_s0.esiti) == 15
    assert all(e.eseguito for e in risultato_s0.esiti)


def test_s0_costa_minuti_non_ore(risultato_s0):
    """Il vincolo di costo: nessun gate legge un file per intero."""
    assert risultato_s0.secondi < 300, f"S0 ha impiegato {risultato_s0.secondi:.0f} s"


def test_s0_produce_gli_artefatti_con_checksum(risultato_s0, configurazione):
    from amplicon16s.io_layer.artifacts import AlberoOutput, Fase

    prodotti = {p.name for p in risultato_s0.artefatti}
    assert prodotti == {
        "gates.json", "crosswalk.tsv", "inventario.json", "letture_ispezionate.tsv"
    }
    albero = AlberoOutput(risultato_s0.artefatti[0].parent.parent)
    assert albero.fase_completa(Fase.INPUT_VALIDATION)
    assert albero.non_integri(Fase.INPUT_VALIDATION) == ()


def test_s0_registra_il_denominatore_di_prevalenza(risultato_s0):
    documento = json.loads(
        next(p for p in risultato_s0.artefatti if p.name == "inventario.json")
        .read_text(encoding="utf-8")
    )
    assert documento["campioni"] == CAMPIONI_ATTESI
    assert documento["denominatore_prevalenza"] == BIOLOGICI_ATTESI


def test_g08_non_avvisa_sulle_dieci_piastre_reali(risultato_s0, configurazione):
    """La piastra a composizione diversa non produce alcun avviso.

    Ha cinque controlli negativi contro gli otto delle altre, e
    decontam.min_blanks ne chiede cinque: passa per una ragione di metodo, non
    per un'eccezione.
    """
    if configurazione.io.batch_table is None:
        pytest.skip("la configurazione non indica il file di arricchimento")
    g08 = next(e for e in risultato_s0.esiti if e.gate == "G08")
    assert g08.superato
    assert g08.avvisi == ()


def test_g10_passa_anche_sui_controlli_negativi(risultato_s0, configurazione):
    """Il controllo positivo non pretende il segnale del bersaglio dai bianchi."""
    g10 = next(e for e in risultato_s0.esiti if e.gate == "G10")
    assert g10.superato
    assert g10.violazioni == ()


def test_il_motivo_non_distingue_segnale_e_contaminazione(risultato_s0, configurazione):
    """Misura che fissa cosa il controllo positivo misura davvero.

    Il motivo conservato compare nei controlli negativi quanto nei campioni
    biologici: i bianchi amplificano contaminanti, e i contaminanti sono
    batteri con lo stesso 16S. Il controllo positivo verifica quindi la
    presenza della regione amplificata, non la qualità del campione, e questo
    test esiste perché nessuno lo usi per ciò che non fa.
    """
    import statistics

    righe = {
        r["accession"]: float(r["frazione_motivo"])
        for r in csv.DictReader(
            (next(p for p in risultato_s0.artefatti if p.name == "letture_ispezionate.tsv"))
            .open(encoding="utf-8"),
            delimiter="\t",
        )
    }
    per_classe = {}
    for campione in risultato_s0.inventario:
        per_classe.setdefault(campione.classe, []).append(righe[campione.accession])

    negativi = per_classe[ClasseCampione.CONTROLLO_NEGATIVO]
    biologici = per_classe[ClasseCampione.BIOLOGICO]

    assert min(negativi) > 0.70, "i bianchi portano il motivo quanto i biologici"
    assert statistics.median(negativi) > statistics.median(biologici) * 0.9

    # Ed e' per questo che la mediana serve: qualche biologico legittimo sta
    # sotto la soglia, e una verifica file per file lo farebbe fallire.
    sotto_soglia = [f for f in biologici if f < configurazione.qc.min_motif_frac]
    assert sotto_soglia, "ci sono biologici sotto la soglia"
    assert min(sotto_soglia) < 0.10


def test_g09_fallisce_se_il_troncamento_supera_il_minimo(configurazione, tmp_path):
    """Con 152 il troncamento supera le letture piu' corte, che sono di 137 bp."""
    from amplicon16s.steps.s00_validate import esegui_s0

    dati = configurazione.model_dump(mode="python")
    dati["io"]["out_root"] = str(tmp_path / "out")
    dati["filter"]["truncLen"] = 152
    risultato = esegui_s0(valida(dati), solleva=False)

    g09 = next(e for e in risultato.esiti if e.gate == "G09")
    assert not g09.superato
    assert {v.codice for v in g09.violazioni} == {"E-S0-09"}
    assert "137 bp" in str(g09.violazioni[0])


def test_le_letture_non_contengono_il_primer(risultato_s0, configurazione):
    """Dato misurato: le letture iniziano a valle del sito 515F."""
    righe = list(
        csv.DictReader(
            (next(p for p in risultato_s0.artefatti if p.name == "letture_ispezionate.tsv"))
            .open(encoding="utf-8"),
            delimiter="\t",
        )
    )
    assert len(righe) == CAMPIONI_ATTESI
    assert max(float(r["frazione_primer"]) for r in righe) == 0.0
