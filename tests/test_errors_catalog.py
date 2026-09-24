"""Test del catalogo degli errori e delle eccezioni che vi si appoggiano.

Il catalogo è la sola fonte di verità dei codici: se una voce è incompleta,
l'utente si ritrova davanti a un codice che non gli dice cosa fare. Questi
test impediscono che una voce possa restare a metà.
"""

from __future__ import annotations

import pytest

from amplicon16s.config import defaults
from amplicon16s.errors.catalog import (
    CATALOGO,
    Categoria,
    codici_con_retry,
    codici_di_fase,
    voce,
)
from amplicon16s.errors.exceptions import (
    DegradazioneRichiesta,
    ErrorePipeline,
    ErroreRevisioneUmana,
    ErroreRitentabile,
    ErroreRitentabileConRevisione,
    errore,
)
from amplicon16s.gates.g01_g15 import CONTROLLI

#: I codici delle fasi della pipeline, come previsti dalla specifica. L'elenco
#: è ripetuto qui apposta: se una voce sparisse dal catalogo per errore, il
#: catalogo da solo resterebbe coerente e nessuno se ne accorgerebbe.
CODICI_DI_FASE = (
    "E-S0-01", "E-S0-02", "E-S0-03", "E-S0-04", "E-S0-05", "E-S0-06", "E-S0-07",
    "E-S0-08", "E-S0-09", "E-S0-10", "E-S0-11", "E-S0-12",
    "E-S0-13", "E-S0-14", "E-S0-15",
    "E-S1-01",
    "E-S2-01", "E-S2-02", "E-S2-03",
    "E-S3-01",
    "E-S4-02",
    "E-S5-01",
    "E-S6-01", "E-S6-02",
    "E-S8-02",
    "E-S9-01",
    "E-S10-01",
    "E-S11-02", "E-S11-03",
    "E-S12-02",
    "E-S13-01", "E-S13-02",
    "E-S14-01",
)

#: I codici del ponte verso R: non appartengono a una fase, descrivono cio' che
#: una fase non puo' dichiarare da se'. Anch'essi elencati apposta.
CODICI_DEL_PONTE = ("E-R-01", "E-R-02", "E-R-03", "E-R-04")

#: Il codice del grafo delle fasi: una fase avviata prima delle sue dipendenze.
CODICI_DEL_GRAFO = ("E-GRAFO-01",)

#: Elenco chiuso dei codici ammessi al retry automatico.
AMMESSI_AL_RETRY = ("E-S2-03", "E-S3-01", "E-S4-02", "E-S5-01")


# --------------------------------------------------------------------------- #
# Completezza del catalogo                                                     #
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize("codice", sorted(CATALOGO))
def test_ogni_voce_ha_messaggio_e_categoria(codice):
    """Nessuna voce può restare a metà: sarebbe un codice muto."""
    v = CATALOGO[codice]

    assert v.codice == codice
    assert v.fase, f"{codice}: fase non compilata"
    assert v.sintesi.strip(), f"{codice}: sintesi vuota"
    assert v.azione.strip(), f"{codice}: azione vuota"
    assert isinstance(v.categoria, Categoria)


@pytest.mark.parametrize("codice", sorted(CATALOGO))
def test_il_messaggio_dice_anche_cosa_fare(codice):
    """La sintesi constata, l'azione prescrive: devono essere cose diverse."""
    v = CATALOGO[codice]
    assert v.azione != v.sintesi
    assert len(v.azione) > 30, f"{codice}: azione troppo breve per essere operativa"
    assert v.sintesi in v.messaggio and v.azione in v.messaggio


@pytest.mark.parametrize("codice", CODICI_DI_FASE + CODICI_DEL_PONTE + CODICI_DEL_GRAFO)
def test_tutti_i_codici_previsti_sono_catalogati(codice):
    assert codice in CATALOGO


def test_la_fase_corrisponde_al_codice():
    for codice, v in CATALOGO.items():
        atteso = codice.removeprefix("E-").rsplit("-", 1)[0]
        assert v.fase == atteso, f"{codice}: fase {v.fase}, attesa {atteso}"


def test_nessun_codice_inatteso_nel_catalogo():
    dai_gate = {c.codice for c in CONTROLLI}
    assert set(CATALOGO) == (
        set(CODICI_DI_FASE) | set(CODICI_DEL_PONTE) | set(CODICI_DEL_GRAFO) | dai_gate
    )


def test_codice_sconosciuto_solleva_un_errore_esplicito():
    with pytest.raises(KeyError, match="non presente nel catalogo"):
        voce("E-S99-99")


def test_codici_di_fase_raggruppa():
    assert codici_di_fase("S6") == ("E-S6-01", "E-S6-02")
    assert codici_di_fase("S7") == ()


# --------------------------------------------------------------------------- #
# Il retry è un elenco chiuso                                                  #
# --------------------------------------------------------------------------- #


def test_i_codici_ammessi_al_retry_sono_esattamente_quattro():
    ammessi = codici_con_retry()
    assert len(ammessi) == 4
    assert ammessi == frozenset(AMMESSI_AL_RETRY)


@pytest.mark.parametrize("codice", AMMESSI_AL_RETRY)
def test_ogni_codice_ammesso_al_retry_esiste_nel_catalogo(codice):
    assert codice in CATALOGO
    assert CATALOGO[codice].ammette_retry


@pytest.mark.parametrize(
    "codice", [c for c in sorted(CATALOGO) if c not in AMMESSI_AL_RETRY]
)
def test_nessun_altro_codice_ammette_il_retry(codice):
    """Il verso opposto: fuori dall'elenco non si ritenta."""
    assert not CATALOGO[codice].ammette_retry, (
        f"{codice} risulta ritentabile ma non e' nell'elenco chiuso"
    )


def test_l_elenco_del_catalogo_coincide_con_retry_whitelist():
    """La configurazione e il catalogo devono dire la stessa cosa.

    retry.whitelist è ciò che la pipeline leggerà a runtime; il catalogo è ciò
    che dichiara quali codici siano ritentabili. Se divergessero, il
    comportamento seguirebbe la configurazione e la documentazione mentirebbe.
    """
    assert set(defaults.RETRY_WHITELIST) == codici_con_retry()


# --------------------------------------------------------------------------- #
# Le eccezioni seguono la categoria del codice                                 #
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize(
    ("codice", "classe"),
    [
        ("E-S0-01", ErroreRevisioneUmana),
        ("E-S2-03", ErroreRitentabile),
        ("E-S3-01", ErroreRitentabileConRevisione),
        ("E-S6-02", DegradazioneRichiesta),
    ],
)
def test_la_classe_discende_dalla_categoria(codice, classe):
    assert isinstance(errore(codice), classe)


@pytest.mark.parametrize("codice", sorted(CATALOGO))
def test_ogni_codice_produce_un_errore_coerente(codice):
    e = errore(codice)
    assert isinstance(e, ErrorePipeline)
    assert e.codice == codice
    assert e.categoria is CATALOGO[codice].categoria
    assert e.ammette_retry == CATALOGO[codice].ammette_retry


def test_sollevare_un_codice_con_la_classe_sbagliata_e_respinto():
    """Chi solleva dichiara cosa e' successo, non come va gestito."""
    with pytest.raises(TypeError, match="revisione_umana"):
        ErroreRitentabile("E-S0-01")


def test_il_messaggio_contiene_codice_dettaglio_e_azione():
    e = errore("E-S2-01", "campione X1 azzerato")
    testo = str(e)
    assert "E-S2-01" in testo
    assert "campione X1 azzerato" in testo
    assert CATALOGO["E-S2-01"].azione in testo


def test_l_errore_si_traduce_in_evento_strutturato():
    e = errore("E-S4-02", "lotto 7", lotto=7, memoria_mb=32000)
    evento = e.come_evento()

    assert evento["codice"] == "E-S4-02"
    assert evento["fase"] == "S4"
    assert evento["categoria"] == "retry_automatico"
    assert evento["lotto"] == 7
    assert evento["memoria_mb"] == 32000


# --------------------------------------------------------------------------- #
# Convivenza con i codici dei gate                                             #
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize("controllo", CONTROLLI, ids=lambda c: c.codice)
def test_ogni_controllo_del_gate_e_nel_catalogo(controllo):
    """Un solo insieme di codici, non due meccanismi paralleli."""
    assert controllo.codice in CATALOGO
    assert controllo.descrizione == CATALOGO[controllo.codice].sintesi
    assert controllo.categoria is CATALOGO[controllo.codice].categoria


@pytest.mark.parametrize("codice", CODICI_DEL_PONTE + CODICI_DEL_GRAFO)
def test_i_codici_del_ponte_e_del_grafo_chiedono_revisione_umana(codice):
    """Un guasto del processo R, o un ordine violato, va capito prima di ritentare."""
    assert CATALOGO[codice].categoria is Categoria.REVISIONE_UMANA


def test_nessun_controllo_di_configurazione_e_ritentabile():
    """Ritentare con la stessa configurazione darebbe lo stesso esito."""
    for controllo in CONTROLLI:
        assert controllo.categoria is Categoria.REVISIONE_UMANA


# --------------------------------------------------------------------------- #
# Proprietà delle categorie                                                    #
# --------------------------------------------------------------------------- #


#: Codici che proseguono con un ripiego. Il ripiego non inventa mai una
#: decisione scientifica: usa un valore gia' scelto dall'operatore, oppure
#: salta una fase che era dichiarata opzionale.
DEGRADAZIONI = ("E-S0-15", "E-S1-01", "E-S6-02", "E-S11-02")


def test_le_degradazioni_sono_quelle_previste():
    degradano = {
        c for c, v in CATALOGO.items()
        if v.categoria is Categoria.DEGRADAZIONE_AUTOMATICA
    }
    assert degradano == set(DEGRADAZIONI)


def test_la_guardia_sulla_filogenesi_ferma_l_esecuzione():
    """E-S9-01 non e' una degradazione.

    La fase viene saltata in silenzio quando phylo.enabled e' false, e in quel
    caso non si genera alcun codice. Il codice si presenta solo se l'albero e'
    stato richiesto esplicitamente: saltarlo allora consegnerebbe un oggetto
    privo proprio della componente domandata, quindi la guardia ferma
    l'esecuzione e lascia la scelta all'operatore.
    """
    v = CATALOGO["E-S9-01"]
    assert v.categoria is Categoria.REVISIONE_UMANA
    assert "phylo.max_seqs" in v.azione
    assert "phylo.enabled" in v.azione


def test_le_categorie_sono_quattro():
    assert len(Categoria) == 4


def test_ogni_categoria_e_usata_almeno_una_volta():
    """Una categoria senza codici sarebbe una distinzione senza differenza."""
    usate = {v.categoria for v in CATALOGO.values()}
    assert usate == set(Categoria)


@pytest.mark.parametrize(
    ("categoria", "retry", "ferma"),
    [
        (Categoria.REVISIONE_UMANA, False, True),
        (Categoria.RETRY_AUTOMATICO, True, False),
        (Categoria.RETRY_POI_REVISIONE, True, True),
        (Categoria.DEGRADAZIONE_AUTOMATICA, False, False),
    ],
)
def test_proprieta_delle_categorie(categoria, retry, ferma):
    assert categoria.ammette_retry is retry
    assert categoria.ferma_esecuzione is ferma
