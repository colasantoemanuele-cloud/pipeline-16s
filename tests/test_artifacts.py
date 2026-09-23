"""Test dell'albero degli artefatti e dei manifesti dei checksum."""

from __future__ import annotations

import json

import pytest

from amplicon16s.io_layer.artifacts import NOME_MANIFESTO, AlberoOutput, Fase
from amplicon16s.io_layer.checksums import (
    ALGORITMO,
    checksum_bytes,
    checksum_file,
    corrisponde,
)

CARTELLE_ATTESE = (
    "00_config", "01_input_validation", "02_qc_profiles", "03_filtered",
    "04_error_models", "05_asv_inference", "06_seqtab", "07_chimera",
    "08_taxonomy", "09_phylogeny", "10_phyloseq", "11_controls", "12_final",
    "99_logs",
)


@pytest.fixture
def albero(tmp_path) -> AlberoOutput:
    return AlberoOutput(tmp_path)


# --------------------------------------------------------------------------- #
# Checksum                                                                     #
# --------------------------------------------------------------------------- #


def test_il_checksum_dichiara_il_proprio_algoritmo(tmp_path):
    file = tmp_path / "x.txt"
    file.write_text("contenuto")
    assert checksum_file(file).startswith(f"{ALGORITMO}:")


def test_checksum_di_file_e_di_byte_coincidono(tmp_path):
    dati = b"contenuto qualunque"
    file = tmp_path / "x.bin"
    file.write_bytes(dati)
    assert checksum_file(file) == checksum_bytes(dati)


def test_contenuti_diversi_hanno_checksum_diversi():
    assert checksum_bytes(b"a") != checksum_bytes(b"b")


def test_un_file_grande_non_viene_caricato_in_memoria(tmp_path):
    """La lettura è a blocchi: il contenuto supera la dimensione del blocco."""
    dati = b"x" * (3 * 1024 * 1024 + 7)
    file = tmp_path / "grande.bin"
    file.write_bytes(dati)
    assert checksum_file(file) == checksum_bytes(dati)


def test_un_file_assente_non_corrisponde(tmp_path):
    """Non è un errore da sollevare: è la risposta «la fase non è completa»."""
    assert corrisponde(tmp_path / "mai_scritto", checksum_bytes(b"")) is False


# --------------------------------------------------------------------------- #
# Albero delle cartelle                                                        #
# --------------------------------------------------------------------------- #


def test_le_fasi_sono_quattordici():
    assert len(list(Fase)) == 14


def test_i_nomi_delle_cartelle_sono_quelli_previsti():
    assert tuple(f.value for f in Fase) == CARTELLE_ATTESE


def test_creare_l_albero_produce_tutte_le_cartelle(albero, tmp_path):
    albero.crea()
    presenti = sorted(p.name for p in tmp_path.iterdir() if p.is_dir())
    assert presenti == sorted(CARTELLE_ATTESE)


def test_creare_l_albero_e_ripetibile(albero, tmp_path):
    albero.crea()
    albero.scrivi_testo(Fase.FINAL, "nota.txt", "contenuto")
    albero.crea()  # non deve cancellare nulla
    assert (tmp_path / Fase.FINAL.value / "nota.txt").is_file()


def test_costruire_l_albero_non_tocca_il_filesystem(tmp_path):
    radice = tmp_path / "non_creata"
    AlberoOutput(radice)
    assert not radice.exists()


def test_la_cartella_nasce_alla_prima_scrittura(albero, tmp_path):
    assert not (tmp_path / Fase.SEQTAB.value).exists()
    albero.scrivi_testo(Fase.SEQTAB, "x.txt", "y")
    assert (tmp_path / Fase.SEQTAB.value).is_dir()


# --------------------------------------------------------------------------- #
# Manifesto                                                                    #
# --------------------------------------------------------------------------- #


def test_l_artefatto_scritto_finisce_nel_manifesto(albero):
    artefatto = albero.scrivi_testo(Fase.TAXONOMY, "taxa.tsv", "a\tb\n")

    manifesto = albero.manifesto(Fase.TAXONOMY)
    assert manifesto["taxa.tsv"]["checksum"] == artefatto.checksum
    assert manifesto["taxa.tsv"]["byte"] == len("a\tb\n")


def test_il_manifesto_e_interrogabile(albero):
    albero.scrivi_testo(Fase.QC_PROFILES, "profili.tsv", "x")
    documento = json.loads(
        (albero.cartella(Fase.QC_PROFILES) / NOME_MANIFESTO).read_text(encoding="utf-8")
    )
    assert documento["fase"] == "02_qc_profiles"
    assert [v["nome"] for v in documento["artefatti"]] == ["profili.tsv"]


def test_il_manifesto_non_dipende_dall_ordine_di_scrittura(tmp_path):
    uno, due = AlberoOutput(tmp_path / "a"), AlberoOutput(tmp_path / "b")
    for nome in ("z.txt", "a.txt", "m.txt"):
        uno.scrivi_testo(Fase.FILTERED, nome, nome)
    for nome in ("m.txt", "z.txt", "a.txt"):
        due.scrivi_testo(Fase.FILTERED, nome, nome)

    assert (uno.cartella(Fase.FILTERED) / NOME_MANIFESTO).read_text() == (
        due.cartella(Fase.FILTERED) / NOME_MANIFESTO
    ).read_text()


def test_riscrivere_un_artefatto_aggiorna_la_sua_voce(albero):
    albero.scrivi_testo(Fase.SEQTAB, "tab.tsv", "primo")
    secondo = albero.scrivi_testo(Fase.SEQTAB, "tab.tsv", "secondo")

    manifesto = albero.manifesto(Fase.SEQTAB)
    assert len(manifesto) == 1
    assert manifesto["tab.tsv"]["checksum"] == secondo.checksum


def test_il_manifesto_di_una_fase_mai_eseguita_e_vuoto(albero):
    assert albero.manifesto(Fase.PHYLOGENY) == {}


# --------------------------------------------------------------------------- #
# Integrità                                                                    #
# --------------------------------------------------------------------------- #


def test_un_artefatto_alterato_non_e_integro(albero):
    artefatto = albero.scrivi_testo(Fase.FINAL, "ps_final.txt", "contenuto originale")
    assert artefatto.integro

    artefatto.percorso.write_text("contenuto manomesso", encoding="utf-8")

    assert not artefatto.integro
    assert albero.non_integri(Fase.FINAL) == ("ps_final.txt",)
    assert not albero.fase_completa(Fase.FINAL)


def test_un_artefatto_rimosso_non_e_integro(albero):
    albero.scrivi_testo(Fase.CHIMERA, "chimere.tsv", "x")
    (albero.cartella(Fase.CHIMERA) / "chimere.tsv").unlink()

    assert albero.non_integri(Fase.CHIMERA) == ("chimere.tsv",)
    assert not albero.fase_completa(Fase.CHIMERA)


def test_un_troncamento_viene_rilevato(albero):
    """Un file interrotto a meta' scrittura esiste ma non e' quello prodotto."""
    albero.scrivi_testo(Fase.ERROR_MODELS, "modello.txt", "a" * 1000)
    (albero.cartella(Fase.ERROR_MODELS) / "modello.txt").write_text("a" * 400)
    assert albero.non_integri(Fase.ERROR_MODELS) == ("modello.txt",)


def test_una_fase_con_artefatti_integri_risulta_completa(albero):
    albero.scrivi_testo(Fase.CONTROLS, "katharoseq.tsv", "x")
    albero.scrivi_testo(Fase.CONTROLS, "positivi.tsv", "y")
    assert albero.fase_completa(Fase.CONTROLS)
    assert albero.non_integri(Fase.CONTROLS) == ()


def test_una_fase_mai_eseguita_non_risulta_completa(albero):
    assert not albero.fase_completa(Fase.ASV_INFERENCE)


def test_un_manifesto_parziale_non_significa_fase_conclusa(albero):
    """Alla ripresa serve sapere se ci sono *quegli* artefatti, non «qualcosa»."""
    albero.scrivi_testo(Fase.SEQTAB, "seqtab.tsv", "x")
    assert albero.fase_completa(Fase.SEQTAB)
    assert not albero.fase_completa(Fase.SEQTAB, attesi=["seqtab.tsv", "conteggi.tsv"])


# --------------------------------------------------------------------------- #
# Artefatti prodotti fuori da Python                                           #
# --------------------------------------------------------------------------- #


def test_un_file_scritto_da_altri_puo_essere_registrato(albero):
    """Le fasi di calcolo sono processi R: scrivono i propri file da se'."""
    cartella = albero.prepara(Fase.PHYLOSEQ)
    (cartella / "ps.rds").write_bytes(b"finto contenuto rds")

    artefatto = albero.registra(Fase.PHYLOSEQ, "ps.rds")

    assert artefatto.integro
    assert albero.manifesto(Fase.PHYLOSEQ)["ps.rds"]["checksum"] == artefatto.checksum
    assert albero.fase_completa(Fase.PHYLOSEQ)


def test_registrare_un_file_inesistente_e_un_errore(albero):
    albero.prepara(Fase.PHYLOSEQ)
    with pytest.raises(FileNotFoundError, match="non trovato"):
        albero.registra(Fase.PHYLOSEQ, "mai_scritto.rds")


def test_artefatti_binari(albero):
    dati = b"\x00\x01\x02 contenuto binario"
    artefatto = albero.scrivi_bytes(Fase.SEQTAB, "seqtab.bin", dati)
    assert artefatto.percorso.read_bytes() == dati
    assert artefatto.checksum == checksum_bytes(dati)
