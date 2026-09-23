"""Test dello schema di configurazione.

Il file di configurazione è il punto di controllo dell'intera esecuzione: un
errore al suo interno deve emergere prima che venga allocato qualunque calcolo.
Questi test verificano che emerga davvero, e che il messaggio dica quale chiave
è il problema.
"""

from __future__ import annotations

import copy
from pathlib import Path

import pytest
import yaml

from amplicon16s.config import defaults
from amplicon16s.config.resolve import (
    NOME_FILE_RISOLTO,
    ConfigRisolta,
    risolvi,
    scrivi_risolta,
)
from amplicon16s.config.schema import (
    PARAMETRI_DERIVATI,
    Config,
    ErroreConfigurazione,
    carica,
    chiavi_schema,
    gruppi_schema,
    valida,
)
from amplicon16s.gates.g01_g15 import CONTROLLI, ErroreGate, esegui_g15
from amplicon16s.io_layer.artifacts import Fase

RADICE = Path(__file__).resolve().parents[1]
ESEMPIO = RADICE / "config" / "config.example.yaml"

#: Parametri obbligatori dichiarati senza valore predefinito.
OBBLIGATORI = (
    "io.fastq_dir",
    "io.assay_table",
    "io.study_table",
    "io.out_root",
    "tax.ref_fasta",
    "tax.ref_name",
    "tax.ref_version",
    "tax.ref_md5",
    "run.container",
)


@pytest.fixture
def dati_esempio() -> dict:
    """Contenuto del file di esempio, ricaricato per ogni test."""
    return yaml.safe_load(ESEMPIO.read_text(encoding="utf-8"))


def _senza(dati: dict, chiave: str) -> dict:
    gruppo, campo = chiave.split(".")
    copia = copy.deepcopy(dati)
    del copia[gruppo][campo]
    return copia


# --------------------------------------------------------------------------- #
# Il file di esempio è un'istanza valida                                       #
# --------------------------------------------------------------------------- #


def test_esempio_esiste():
    assert ESEMPIO.is_file(), f"file di esempio non trovato in {ESEMPIO}"


def test_esempio_si_carica_e_valida():
    config = carica(ESEMPIO)
    assert isinstance(config, Config)


def test_esempio_conserva_i_tipi_numerici(dati_esempio):
    """In YAML l'esponente richiede il segno: "1e8" sarebbe letto come stringa."""
    assert isinstance(dati_esempio["err"]["nbases"], float)
    assert isinstance(dati_esempio["dada"]["omega_a"], float)


def test_esempio_dichiara_versione_come_stringa(dati_esempio):
    assert isinstance(dati_esempio["tax"]["ref_version"], str)


# --------------------------------------------------------------------------- #
# Corrispondenza biunivoca fra schema e file di esempio                        #
# --------------------------------------------------------------------------- #


def _chiavi_del_file(dati: dict) -> set[str]:
    return {
        f"{gruppo}.{campo}"
        for gruppo, contenuto in dati.items()
        if isinstance(contenuto, dict)
        for campo in contenuto
    }


def test_nessuna_chiave_dello_schema_manca_nel_file(dati_esempio):
    mancanti = chiavi_schema() - _chiavi_del_file(dati_esempio)
    assert not mancanti, f"chiavi dello schema assenti dal file: {sorted(mancanti)}"


def test_nessuna_chiave_del_file_e_sconosciuta_allo_schema(dati_esempio):
    estranee = _chiavi_del_file(dati_esempio) - chiavi_schema()
    assert not estranee, f"chiavi del file assenti dallo schema: {sorted(estranee)}"


def test_tutti_i_gruppi_compaiono_nel_file(dati_esempio):
    """I gruppi vuoti non hanno foglie: vanno confrontati a parte."""
    mancanti = set(gruppi_schema()) - set(dati_esempio)
    assert not mancanti, f"gruppi assenti dal file: {sorted(mancanti)}"
    estranei = set(dati_esempio) - set(gruppi_schema())
    assert not estranei, f"gruppi del file assenti dallo schema: {sorted(estranei)}"


# --------------------------------------------------------------------------- #
# Parametri obbligatori                                                        #
# --------------------------------------------------------------------------- #


def test_elenco_obbligatori_coincide_con_lo_schema():
    dallo_schema = {
        f"{gruppo}.{campo}"
        for gruppo, descrittore in Config.model_fields.items()
        for campo, sotto in descrittore.annotation.model_fields.items()
        if sotto.is_required()
    }
    assert dallo_schema == set(OBBLIGATORI)


@pytest.mark.parametrize("chiave", OBBLIGATORI)
def test_obbligatorio_mancante_viene_segnalato_col_suo_nome(chiave, dati_esempio):
    with pytest.raises(ErroreConfigurazione) as errore:
        valida(_senza(dati_esempio, chiave))

    messaggio = str(errore.value)
    assert chiave in messaggio, f"il messaggio non nomina {chiave}: {messaggio}"
    assert "obbligatorio" in messaggio


def test_gruppo_obbligatorio_mancante_viene_segnalato(dati_esempio):
    del dati_esempio["io"]
    with pytest.raises(ErroreConfigurazione, match="io"):
        valida(dati_esempio)


# --------------------------------------------------------------------------- #
# Tipi e vincoli di dominio                                                    #
# --------------------------------------------------------------------------- #


def test_tipo_errato_viene_respinto_con_messaggio_leggibile(dati_esempio):
    dati_esempio["filter"]["maxEE"] = "due"

    with pytest.raises(ErroreConfigurazione) as errore:
        valida(dati_esempio)

    messaggio = str(errore.value)
    assert "filter.maxEE" in messaggio
    assert "'due'" in messaggio


@pytest.mark.parametrize(
    ("gruppo", "campo", "valore"),
    [
        ("run", "threads", 0),               # deve essere positivo
        ("run", "threads", -4),
        ("filter", "truncLen", 0),
        ("tax", "min_boot", 101),            # percentuale fuori scala
        ("decontam", "threshold", 1.5),      # frazione fuori da [0, 1]
        ("qc", "warn_frac_chimeric", -0.1),
        ("chimera", "min_sample_fraction", 2.0),
    ],
)
def test_valore_fuori_dominio_viene_respinto(gruppo, campo, valore, dati_esempio):
    dati_esempio[gruppo][campo] = valore
    with pytest.raises(ErroreConfigurazione, match=f"{gruppo}.{campo}"):
        valida(dati_esempio)


@pytest.mark.parametrize(
    ("gruppo", "campo", "valore"),
    [
        ("chimera", "method", "consenso"),
        ("decontam", "method", "prevalenza"),
        ("out", "serialization", "qs"),
        ("dada", "pool", "pseudoo"),
    ],
)
def test_valore_fuori_vocabolario_viene_respinto(gruppo, campo, valore, dati_esempio):
    dati_esempio[gruppo][campo] = valore
    with pytest.raises(ErroreConfigurazione, match=f"{gruppo}.{campo}"):
        valida(dati_esempio)


def test_dada_pool_accetta_il_booleano(dati_esempio):
    """dada2 ammette TRUE, FALSE o "pseudo"."""
    dati_esempio["dada"]["pool"] = True
    assert valida(dati_esempio).dada.pool is True


def test_booleano_scritto_come_stringa_viene_respinto(dati_esempio):
    """Una stringa al posto di un booleano è un errore, non un valore da convertire."""
    dati_esempio["phylo"]["enabled"] = "true"
    with pytest.raises(ErroreConfigurazione, match="phylo.enabled"):
        valida(dati_esempio)


def test_regex_malformata_viene_respinta(dati_esempio):
    dati_esempio["meta"]["module_regex"] = "^([A-Z]"
    with pytest.raises(ErroreConfigurazione, match="meta.module_regex"):
        valida(dati_esempio)


def test_immagine_senza_digest_viene_respinta(dati_esempio):
    """Un tag può essere riassegnato a un'immagine diversa, un digest no."""
    dati_esempio["run"]["container"] = "amplicon16s:dev"
    with pytest.raises(ErroreConfigurazione, match="run.container"):
        valida(dati_esempio)


def test_md5_malformato_viene_respinto(dati_esempio):
    dati_esempio["tax"]["ref_md5"] = "non-un-md5"
    with pytest.raises(ErroreConfigurazione, match="tax.ref_md5"):
        valida(dati_esempio)


# --------------------------------------------------------------------------- #
# Chiavi sconosciute                                                           #
# --------------------------------------------------------------------------- #


def test_refuso_in_un_parametro_viene_respinto(dati_esempio):
    """Senza questo vincolo il refuso passerebbe e varrebbe il predefinito."""
    dati_esempio["filter"]["maxEEE"] = 2
    with pytest.raises(ErroreConfigurazione, match="filter.maxEEE"):
        valida(dati_esempio)


def test_gruppo_sconosciuto_viene_respinto(dati_esempio):
    dati_esempio["inventato"] = {"x": 1}
    with pytest.raises(ErroreConfigurazione, match="inventato"):
        valida(dati_esempio)


@pytest.mark.parametrize("gruppo", ["norm", "glom", "beta", "ord", "stat"])
def test_gruppi_ecologici_non_accettano_parametri(gruppo, dati_esempio):
    """Finché nessun codice sa interpretarli, vanno respinti."""
    dati_esempio[gruppo] = {"metodo": "qualcosa"}
    with pytest.raises(ErroreConfigurazione, match=f"{gruppo}.metodo"):
        valida(dati_esempio)


# --------------------------------------------------------------------------- #
# Coerenza fra campi                                                           #
# --------------------------------------------------------------------------- #


def test_avviso_non_puo_superare_arresto(dati_esempio):
    dati_esempio["qc"]["warn_frac_chimeric"] = 0.8
    dati_esempio["qc"]["stop_frac_chimeric"] = 0.5
    with pytest.raises(ErroreConfigurazione, match="warn_frac_chimeric"):
        valida(dati_esempio)


def test_tolleranza_non_puo_raggiungere_il_troncamento(dati_esempio):
    dati_esempio["filter"]["truncLen_shortfall_warn"] = 137
    with pytest.raises(ErroreConfigurazione, match="truncLen_shortfall_warn"):
        valida(dati_esempio)


def test_categorie_di_controllo_devono_essere_disgiunte(dati_esempio):
    dati_esempio["ctrl"]["blank_values"] = ["Surface swab"]
    with pytest.raises(ErroreConfigurazione, match="ctrl"):
        valida(dati_esempio)


# --------------------------------------------------------------------------- #
# Segnalazione degli errori                                                    #
# --------------------------------------------------------------------------- #


def test_vengono_riportati_tutti_i_problemi_non_solo_il_primo(dati_esempio):
    """Chi corregge un file vuole l'elenco completo, non un errore alla volta."""
    dati_esempio["filter"]["maxEE"] = "due"
    dati_esempio["run"]["threads"] = -1
    del dati_esempio["io"]["fastq_dir"]

    with pytest.raises(ErroreConfigurazione) as errore:
        valida(dati_esempio)

    assert len(errore.value.problemi) == 3
    messaggio = str(errore.value)
    for atteso in ("filter.maxEE", "run.threads", "io.fastq_dir"):
        assert atteso in messaggio


def test_errore_di_caricamento_indica_il_file(tmp_path):
    percorso = tmp_path / "rotto.yaml"
    percorso.write_text("io: [questo non e' una mappa]\n", encoding="utf-8")
    with pytest.raises(ErroreConfigurazione, match="rotto.yaml"):
        carica(percorso)


def test_file_vuoto_viene_segnalato(tmp_path):
    percorso = tmp_path / "vuoto.yaml"
    percorso.write_text("", encoding="utf-8")
    with pytest.raises(ErroreConfigurazione, match="vuoto"):
        carica(percorso)


def test_yaml_malformato_viene_segnalato(tmp_path):
    percorso = tmp_path / "malformato.yaml"
    percorso.write_text("io: {non chiuso\n", encoding="utf-8")
    with pytest.raises(ErroreConfigurazione, match="YAML"):
        carica(percorso)


# --------------------------------------------------------------------------- #
# Predefiniti derivati dal dataset di riferimento                              #
# --------------------------------------------------------------------------- #


def test_chiavi_derivate_dal_dataset_esistono_nello_schema():
    ignote = set(defaults.DERIVATI_DAL_DATASET) - chiavi_schema()
    assert not ignote, f"chiavi inesistenti nello schema: {sorted(ignote)}"


def test_esempio_usa_i_valori_del_dataset_di_riferimento(dati_esempio):
    """Il file di esempio descrive OSD-734: deve restare allineato a defaults."""
    assert dati_esempio["err"]["batch_column"] == defaults.ERR_BATCH_COLUMN
    assert dati_esempio["decontam"]["batch_column"] == defaults.DECONTAM_BATCH_COLUMN
    assert dati_esempio["ctrl"]["column"] == defaults.CTRL_COLUMN
    assert dati_esempio["tax"]["ref_name"] == defaults.TAX_REF_NAME
    assert dati_esempio["tax"]["ref_version"] == defaults.TAX_REF_VERSION
    assert dati_esempio["katharoseq"]["target_taxon"] == defaults.KATHAROSEQ_TARGET_TAXON


def test_whitelist_dei_ritentativi_non_dipende_dal_dataset():
    """L'elenco raccoglie codici di errore, non campioni.

    E' chiuso per scelta metodologica — si ritenta solo dove l'azione
    correttiva non modifica alcuna assunzione dell'analisi — quindi non va
    annoverato fra i valori da rivedere cambiando dataset.
    """
    assert "retry.whitelist" not in defaults.DERIVATI_DAL_DATASET


def test_esempio_resta_allineato_alla_whitelist_predefinita(dati_esempio):
    assert list(dati_esempio["retry"]["whitelist"]) == list(defaults.RETRY_WHITELIST)


# --------------------------------------------------------------------------- #
# Immutabilità                                                                 #
# --------------------------------------------------------------------------- #


def test_la_configurazione_non_si_modifica_dopo_il_caricamento():
    config = carica(ESEMPIO)
    with pytest.raises(Exception):
        config.run.threads = 1


# =========================================================================== #
# Gate G15 — coerenza interna della configurazione                            #
# =========================================================================== #

#: Numero di campioni biologici del dataset di riferimento.
CAMPIONI_BIOLOGICI_OSD734 = 803


def _viola(dati_esempio: dict, modifica) -> ErroreGate:
    """Applica una modifica incoerente e restituisce l'errore sollevato."""
    modifica(dati_esempio)
    with pytest.raises(ErroreGate) as errore:
        esegui_g15(dati_esempio)
    return errore.value


# --------------------------------------------------------------------------- #
# Il gate accetta una configurazione coerente                                  #
# --------------------------------------------------------------------------- #


def test_g15_accetta_il_file_di_esempio(dati_esempio):
    risolta = esegui_g15(dati_esempio)
    assert isinstance(risolta, ConfigRisolta)


def test_registro_dei_controlli_e_coerente():
    codici = [c.codice for c in CONTROLLI]
    assert len(codici) == len(set(codici)), "codici duplicati"
    for controllo in CONTROLLI:
        assert controllo.implementato_da in ("schema", "gate")
        assert controllo.descrizione


# --------------------------------------------------------------------------- #
# Configurazioni incoerenti                                                    #
# --------------------------------------------------------------------------- #
# Una per tipo di controllo. Ogni caso indica il codice atteso e i parametri
# che il messaggio deve nominare: non basta sapere che qualcosa non va, serve
# sapere quali valori sono in conflitto fra loro.

INCOERENTI = [
    pytest.param(
        lambda d: d["qc"].update(warn_frac_chimeric=0.8, stop_frac_chimeric=0.5),
        "E-G15-05",
        ("qc.warn_frac_chimeric", "qc.stop_frac_chimeric"),
        id="avviso-oltre-arresto",
    ),
    pytest.param(
        lambda d: d["decontam"].__setitem__("threshold", 1.0),
        "E-G15-02",
        ("decontam.threshold",),
        id="soglia-decontam-degenere",
    ),
    pytest.param(
        lambda d: d["asv"].__setitem__("len_tol", 200),
        "E-G15-03",
        ("asv.len_min", "asv.len_tol", "filter.truncLen"),
        id="tolleranza-oltre-lunghezza",
    ),
    pytest.param(
        lambda d: d["prev"].__setitem__("min_fraction", 1.5),
        "E-G15-04",
        ("prev.min_fraction",),
        id="frazione-fuori-scala",
    ),
    pytest.param(
        lambda d: d["retry"].__setitem__("max_attempts", -1),
        "E-G15-06",
        ("retry.max_attempts",),
        id="tentativi-negativi",
    ),
    pytest.param(
        lambda d: d["ctrl"].__setitem__("blank_values", []),
        "E-G15-07",
        ("decontam.method", "ctrl.blank_values"),
        id="nessun-controllo-negativo",
    ),
]


@pytest.mark.parametrize(("modifica", "codice", "parametri"), INCOERENTI)
def test_g15_respinge_le_configurazioni_incoerenti(
    modifica, codice, parametri, dati_esempio
):
    errore = _viola(dati_esempio, modifica)

    codici = {v.codice for v in errore.violazioni}
    assert codice in codici, f"atteso {codice}, ottenuti {sorted(codici)}"

    testo = str(errore)
    for parametro in parametri:
        assert parametro in testo, f"il messaggio non nomina {parametro}: {testo}"


def test_g15_termina_con_un_codice_di_errore(dati_esempio):
    errore = _viola(dati_esempio, lambda d: d["decontam"].__setitem__("threshold", 0.0))
    assert errore.codice == "E-G15-02"
    assert errore.gate == "G15"


def test_g15_raccoglie_tutte_le_violazioni(dati_esempio):
    errore = _viola(
        dati_esempio,
        lambda d: (
            d["qc"].update(warn_frac_chimeric=0.9, stop_frac_chimeric=0.1),
            d["retry"].__setitem__("max_attempts", -1),
        ),
    )
    assert len(errore.violazioni) >= 2


# --------------------------------------------------------------------------- #
# Parametri derivati                                                           #
# --------------------------------------------------------------------------- #


def test_derivati_statici_assumono_i_valori_attesi(dati_esempio):
    derivati = esegui_g15(dati_esempio).derivati
    assert derivati.filter_minLen == 137
    assert derivati.asv_len_min == 137
    assert derivati.asv_len_max == 137


def test_derivati_seguono_i_parametri_da_cui_discendono(dati_esempio):
    """Non sono costanti: cambiando l'origine cambiano anche loro."""
    dati_esempio["filter"]["truncLen"] = 150
    dati_esempio["filter"]["trimLeft"] = 10
    dati_esempio["asv"]["len_tol"] = 2

    derivati = esegui_g15(dati_esempio).derivati
    assert derivati.filter_minLen == 150
    assert derivati.asv_len_min == 138  # 150 - 10 - 2
    assert derivati.asv_len_max == 142  # 150 - 10 + 2


def test_la_prima_fase_e_utilizzabile_senza_il_gate(dati_esempio):
    """La risoluzione statica è indipendente dal gate che la usa."""
    risolta = risolvi(valida(dati_esempio))
    assert risolta.derivati.filter_minLen == 137
    assert risolta.derivati.prev_min_samples is None


def test_min_samples_non_e_noto_prima_dei_metadati(dati_esempio):
    risolta = esegui_g15(dati_esempio)
    assert risolta.derivati.prev_min_samples is None
    assert not risolta.completa


def test_min_samples_calcolato_quando_il_dato_diventa_disponibile(dati_esempio):
    completa = esegui_g15(dati_esempio).con_campioni_biologici(
        CAMPIONI_BIOLOGICI_OSD734
    )
    assert completa.derivati.prev_min_samples == 9  # ceil(0.01 * 803)
    assert completa.completa


def test_min_samples_arrotonda_per_eccesso(dati_esempio):
    """Una soglia frazionaria non ha significato: i campioni si contano interi."""
    risolta = esegui_g15(dati_esempio)
    assert risolta.con_campioni_biologici(100).derivati.prev_min_samples == 1
    assert risolta.con_campioni_biologici(101).derivati.prev_min_samples == 2
    assert risolta.con_campioni_biologici(0).derivati.prev_min_samples == 0


def test_numero_di_campioni_negativo_viene_respinto(dati_esempio):
    with pytest.raises(ValueError, match="negativo"):
        esegui_g15(dati_esempio).con_campioni_biologici(-1)


@pytest.mark.parametrize("chiave", PARAMETRI_DERIVATI)
def test_derivato_impostato_a_mano_e_un_errore(chiave, dati_esempio):
    gruppo, campo = chiave.split(".")
    dati_esempio[gruppo][campo] = 1

    with pytest.raises(ErroreGate) as errore:
        esegui_g15(dati_esempio)

    assert "E-G15-08" in {v.codice for v in errore.value.violazioni}
    assert chiave in str(errore.value)


def test_elenco_dei_derivati_coincide_con_quelli_calcolati(dati_esempio):
    calcolati = set(esegui_g15(dati_esempio).derivati.come_chiavi())
    assert calcolati == set(PARAMETRI_DERIVATI)


# --------------------------------------------------------------------------- #
# Digest                                                                       #
# --------------------------------------------------------------------------- #


def test_configurazioni_identiche_hanno_lo_stesso_digest(dati_esempio):
    primo = esegui_g15(copy.deepcopy(dati_esempio)).digest
    secondo = esegui_g15(copy.deepcopy(dati_esempio)).digest
    assert primo == secondo


def test_un_solo_parametro_diverso_cambia_il_digest(dati_esempio):
    originale = esegui_g15(copy.deepcopy(dati_esempio)).digest
    dati_esempio["run"]["seed"] = 101
    assert esegui_g15(dati_esempio).digest != originale


def test_il_digest_non_dipende_dall_ordine_delle_chiavi(dati_esempio):
    """La forma canonica ordina le chiavi: la formattazione non conta."""
    originale = esegui_g15(copy.deepcopy(dati_esempio)).digest
    rovesciato = {k: dati_esempio[k] for k in reversed(list(dati_esempio))}
    assert esegui_g15(rovesciato).digest == originale


def test_il_digest_cambia_se_cambia_un_derivato(dati_esempio):
    originale = esegui_g15(copy.deepcopy(dati_esempio)).digest
    dati_esempio["asv"]["len_tol"] = 3  # sposta len_min e len_max
    assert esegui_g15(dati_esempio).digest != originale


def test_il_digest_non_cambia_fra_le_due_fasi(dati_esempio):
    """Identifica la configurazione, non i dati su cui viene applicata."""
    risolta = esegui_g15(dati_esempio)
    completa = risolta.con_campioni_biologici(CAMPIONI_BIOLOGICI_OSD734)
    assert completa.digest == risolta.digest


def test_il_digest_ha_la_forma_attesa(dati_esempio):
    digest = esegui_g15(dati_esempio).digest
    algoritmo, _, valore = digest.partition(":")
    assert algoritmo == "sha256"
    assert len(valore) == 64
    assert set(valore) <= set("0123456789abcdef")


# --------------------------------------------------------------------------- #
# Registrazione della configurazione usata                                     #
# --------------------------------------------------------------------------- #


def test_resolved_viene_scritto_con_il_digest(dati_esempio, tmp_path):
    risolta = esegui_g15(dati_esempio).con_campioni_biologici(
        CAMPIONI_BIOLOGICI_OSD734
    )
    percorso = scrivi_risolta(risolta, tmp_path)

    assert percorso == tmp_path / Fase.CONFIG.value / NOME_FILE_RISOLTO
    assert percorso.is_file()

    documento = yaml.safe_load(percorso.read_text(encoding="utf-8"))
    assert documento["digest"] == risolta.digest


def test_resolved_contiene_i_parametri_derivati(dati_esempio, tmp_path):
    risolta = esegui_g15(dati_esempio).con_campioni_biologici(
        CAMPIONI_BIOLOGICI_OSD734
    )
    documento = yaml.safe_load(
        scrivi_risolta(risolta, tmp_path).read_text(encoding="utf-8")
    )

    parametri = documento["parametri"]
    assert parametri["filter"]["minLen"] == 137
    assert parametri["asv"]["len_min"] == 137
    assert parametri["asv"]["len_max"] == 137
    assert parametri["prev"]["min_samples"] == 9
    assert set(documento["derivati"]) == set(PARAMETRI_DERIVATI)


def test_resolved_registra_il_derivato_non_ancora_noto(dati_esempio, tmp_path):
    """Prima della lettura dei metadati il valore manca, e il file lo dice."""
    documento = yaml.safe_load(
        scrivi_risolta(esegui_g15(dati_esempio), tmp_path).read_text(encoding="utf-8")
    )
    assert documento["parametri"]["prev"]["min_samples"] is None


def test_resolved_registra_tutti_i_parametri_dichiarati(dati_esempio, tmp_path):
    documento = yaml.safe_load(
        scrivi_risolta(esegui_g15(dati_esempio), tmp_path).read_text(encoding="utf-8")
    )
    scritte = {
        f"{gruppo}.{campo}"
        for gruppo, contenuto in documento["parametri"].items()
        for campo in contenuto
    }
    assert chiavi_schema() <= scritte
    assert set(PARAMETRI_DERIVATI) <= scritte


def test_resolved_e_deterministico(dati_esempio, tmp_path):
    """Due scritture della stessa configurazione danno lo stesso file."""
    risolta = esegui_g15(dati_esempio)
    primo = scrivi_risolta(risolta, tmp_path / "a").read_text(encoding="utf-8")
    secondo = scrivi_risolta(risolta, tmp_path / "b").read_text(encoding="utf-8")
    assert primo == secondo


def test_la_cartella_di_destinazione_viene_creata(dati_esempio, tmp_path):
    out_root = tmp_path / "non" / "ancora" / "esistente"
    percorso = scrivi_risolta(esegui_g15(dati_esempio), out_root)
    assert percorso.is_file()
