"""Suite di verifica del catalogo degli errori, della whitelist di retry e delle eccezioni tipizzate.

Inquadramento nel Piano Operativo:
    - **Settimana di riferimento**: **Settimana 5 (W5 — Fase F1: Catalogo degli
      errori, gestione degli artefatti e logging)**.
    - **Scopo del modulo**: Verifica l'integrità formale e metodologica del
      catalogo centralizzato dei codici di errore (``CATALOGO``) e della
      gerarchia di eccezioni (``ErrorePipeline``). In particolare certifica:
        * Che ogni codice possieda sia una sintesi diagnostica sia un'azione
          operativa prescrittiva in italiano (> 30 caratteri) che indichi
          all'operatore *cosa fare* e non solo *cosa è fallito*;
        * La chiusura ermetica della whitelist dei tentativi ripetuti a
          **esattamente 4 codici** (``E-S2-03``, ``E-S3-01``, ``E-S4-02``,
          ``E-S5-01``), poiché il retry automatico è ammesso solo dove l'azione
          correttiva (es. riduzione di ``run.batch_size`` per OOM) non altera
          alcuna assunzione scientifica dell'analisi;
        * Che tutti i guasti del ponte verso R (``E-R-01`` .. ``E-R-04``),
          del grafo delle precedenze (``E-GRAFO-01``) e dei controlli di
          configurazione (``CONTROLLI`` di G15) appartengano tassativamente a
          ``Categoria.REVISIONE_UMANA``.
    - **Moduli sorgente coperti**:
        * ``src/amplicon16s/errors/catalog.py``
        * ``src/amplicon16s/errors/exceptions.py``
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

#: Elenco esplicito di controllo dei codici di fase (S0–S14): mantenuto nel test
#: per intercettare qualsiasi rimozione accidentale dal dizionario ``CATALOGO``.
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

#: Codici infrastrutturali del ponte verso R (`rbridge`).
CODICI_DEL_PONTE = ("E-R-01", "E-R-02", "E-R-03", "E-R-04")

#: Codice di violazione dell'ordine topologico nel grafo delle fasi.
CODICI_DEL_GRAFO = ("E-GRAFO-01",)

#: Elenco chiuso dei 4 soli codici ammessi al retry automatico nella pipeline.
AMMESSI_AL_RETRY = ("E-S2-03", "E-S3-01", "E-S4-02", "E-S5-01")


# --------------------------------------------------------------------------- #
# Completezza del catalogo                                                     #
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize("codice", sorted(CATALOGO))
def test_ogni_voce_ha_messaggio_e_categoria(codice):
    """
    **Obiettivo**: Verificare che ogni voce del ``CATALOGO`` abbia ``codice``,
    ``fase``, ``sintesi``, ``azione`` non vuoti e una ``categoria`` tipizzata
    ``Categoria``.

    **Razionale Scientifico/Sistemistico**: Impedisce che venga registrato nel
    codice un identificativo d'errore "muto" o incompleto che lascerebbe
    l'operatore e il file ``punto_di_ripresa.json`` privi di istruzioni.
    """
    v = CATALOGO[codice]

    assert v.codice == codice
    assert v.fase, f"{codice}: fase non compilata"
    assert v.sintesi.strip(), f"{codice}: sintesi vuota"
    assert v.azione.strip(), f"{codice}: azione vuota"
    assert isinstance(v.categoria, Categoria)


@pytest.mark.parametrize("codice", sorted(CATALOGO))
def test_il_messaggio_dice_anche_cosa_fare(codice):
    """
    **Obiettivo**: Verificare che per ciascuno dei 47 codici del catalogo il
    campo ``azione`` sia distinto dalla ``sintesi``, abbia lunghezza ``> 30``
    caratteri e compaia nel messaggio completo dell'eccezione.

    **Razionale Scientifico/Sistemistico**: Un messaggio d'errore che si limita
    a constatare il sintomo (es. *"troncamento superiore alle letture"*) non è
    operativo. Imporre un campo ``azione`` distinto e articolato (> 30 caratteri)
    obbliga ogni voce del catalogo a prescrivere il parametro YAML o l'intervento
    esatto richiesto al biologo computazionale per risolvere l'arresto.
    """
    v = CATALOGO[codice]
    assert v.azione != v.sintesi
    assert len(v.azione) > 30, f"{codice}: azione troppo breve per essere operativa"
    assert v.sintesi in v.messaggio and v.azione in v.messaggio


@pytest.mark.parametrize("codice", CODICI_DI_FASE + CODICI_DEL_PONTE + CODICI_DEL_GRAFO)
def test_tutti_i_codici_previsti_sono_catalogati(codice):
    """
    **Obiettivo**: Verificare che ciascuno dei codici di fase, del ponte R e
    del grafo previsti dalla specifica sia effettivamente presente in ``CATALOGO``.

    **Razionale Scientifico/Sistemistico**: Evita che durante il refactoring
    un codice referenziato da una fase R o da un gate venga cancellato dal
    catalogo causando un ``KeyError`` proprio nel momento in cui si verifica il guasto.
    """
    assert codice in CATALOGO


def test_la_fase_corrisponde_al_codice():
    """
    **Obiettivo**: Verificare che l'attributo ``v.fase`` di ogni voce coincida
    con il segmento centrale del codice ``E-<FASE>-<NN>`` (es. ``E-S12-02 -> S12``).

    **Razionale Scientifico/Sistemistico**: Garantisce la coerenza dei filtri di
    ricerca nei log strutturati ``99_logs/pipeline.jsonl`` e nei report di fase.
    """
    for codice, v in CATALOGO.items():
        atteso = codice.removeprefix("E-").rsplit("-", 1)[0]
        assert v.fase == atteso, f"{codice}: fase {v.fase}, attesa {atteso}"


def test_nessun_codice_inatteso_nel_catalogo():
    """
    **Obiettivo**: Verificare che l'insieme delle chiavi di ``CATALOGO`` coincida
    esattamente con l'unione di ``CODICI_DI_FASE``, ``CODICI_DEL_PONTE``,
    ``CODICI_DEL_GRAFO`` e dei ``CONTROLLI`` del Gate G15.

    **Razionale Scientifico/Sistemistico**: Impedisce l'accumulo di codici
    orfani o non documentati nella specifica dell'architettura.
    """
    dai_gate = {c.codice for c in CONTROLLI}
    assert set(CATALOGO) == (
        set(CODICI_DI_FASE) | set(CODICI_DEL_PONTE) | set(CODICI_DEL_GRAFO) | dai_gate
    )


def test_codice_sconosciuto_solleva_un_errore_esplicito():
    """
    **Obiettivo**: Verificare che ``voce("E-S99-99")`` sollevi ``KeyError`` con
    il messaggio esplicito ``"non presente nel catalogo"``.

    **Razionale Scientifico/Sistemistico**: Intercetta immediatamente eventuali
    errori di battitura nei codici d'errore passati a ``errore()`` o ``codice_memoria``.
    """
    with pytest.raises(KeyError, match="non presente nel catalogo"):
        voce("E-S99-99")


def test_codici_di_fase_raggruppa():
    """
    **Obiettivo**: Verificare che ``codici_di_fase("S6")`` restituisca tutti e
    soli i codici appartenenti alla fase indicata (``("E-S6-01", "E-S6-02")``).

    **Razionale Scientifico/Sistemistico**: Supporta l'introspezione per-fase
    del catalogo da parte del generatore di report e della documentazione.
    """
    assert codici_di_fase("S6") == ("E-S6-01", "E-S6-02")
    assert codici_di_fase("S7") == ()


# --------------------------------------------------------------------------- #
# Il retry è un elenco chiuso                                                  #
# --------------------------------------------------------------------------- #


def test_i_codici_ammessi_al_retry_sono_esattamente_quattro():
    """
    **Obiettivo**: Verificare che ``codici_con_retry()`` restituisca esattamente
    i 4 codici ``{"E-S2-03", "E-S3-01", "E-S4-02", "E-S5-01"}``.

    **Razionale Scientifico/Sistemistico**: La chiusura ermetica della whitelist
    di retry a questi 4 soli codici è un pilastro metodologico della pipeline:
    ritentare automaticamente è lecito solo per esaurimenti di memoria in
    `filterAndTrim` (`E-S2-03`), inferenza DADA2 (`E-S4-02`), costruzione
    `seqtab` (`E-S5-01`) o mancata convergenza di `learnErrors` (`E-S3-01`), dove
    l'azione correttiva agisce sul partizionamento dei lotti o su ``err.nbases``
    senza mai rilassare soglie biologiche (come ``maxEE`` o ``decontam.threshold``).
    """
    ammessi = codici_con_retry()
    assert len(ammessi) == 4
    assert ammessi == frozenset(AMMESSI_AL_RETRY)


@pytest.mark.parametrize("codice", AMMESSI_AL_RETRY)
def test_ogni_codice_ammesso_al_retry_esiste_nel_catalogo(codice):
    """
    **Obiettivo**: Verificare che ciascuno dei 4 codici in ``AMMESSI_AL_RETRY``
    esista in ``CATALOGO`` ed esponga ``ammette_retry is True``.

    **Razionale Scientifico/Sistemistico**: Garantisce la consistenza interna
    del flag ``ammette_retry`` sui 4 codici della whitelist.
    """
    assert codice in CATALOGO
    assert CATALOGO[codice].ammette_retry


@pytest.mark.parametrize(
    "codice", [c for c in sorted(CATALOGO) if c not in AMMESSI_AL_RETRY]
)
def test_nessun_altro_codice_ammette_il_retry(codice):
    """
    **Obiettivo**: Verificare che tutti i restanti 43 codici del catalogo abbiano
    tassativamente ``ammette_retry is False``.

    **Razionale Scientifico/Sistemistico**: Impedisce che un nuovo codice
    d'errore biologico o di validazione venga accidentalmente classificato in
    una categoria che abilita il retry automatico.
    """
    assert not CATALOGO[codice].ammette_retry, (
        f"{codice} risulta ritentabile ma non e' nell'elenco chiuso"
    )


def test_l_elenco_del_catalogo_coincide_con_retry_whitelist():
    """
    **Obiettivo**: Verificare l'identità insiemistica tra ``defaults.RETRY_WHITELIST``
    e ``codici_con_retry()``.

    **Razionale Scientifico/Sistemistico**: Allinea la costante predefinita
    dello schema di configurazione (`defaults.py`) con la definizione formale
    del catalogo degli errori (`catalog.py`), garantendo che il Gate G15
    (`E-G15-09`) validi contro lo stesso identico insieme chiuso.
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
    """
    **Obiettivo**: Verificare che la factory ``errore(codice)`` istanzi la
    sottoclasse di ``ErrorePipeline`` corrispondente alla ``Categoria`` del codice.

    **Razionale Scientifico/Sistemistico**: Consente all'``Esecutore`` e a
    ``PipelineStep.esegui()`` di intercettare il comportamento di gestione
    tramite ``except DegradazioneRichiesta:`` o ``isinstance(..., ErroreRitentabile)``
    senza catene di ``if`` sui codici stringa.
    """
    assert isinstance(errore(codice), classe)


@pytest.mark.parametrize("codice", sorted(CATALOGO))
def test_ogni_codice_produce_un_errore_coerente(codice):
    """
    **Obiettivo**: Verificare che ``errore(codice)`` produca un'istanza valida
    di ``ErrorePipeline`` per tutti i 47 codici del catalogo.

    **Razionale Scientifico/Sistemistico**: Assicura che nessuna voce del
    catalogo abbia una categoria non mappata nella tabella di dispatch di
    ``exceptions.py``.
    """
    e = errore(codice)
    assert isinstance(e, ErrorePipeline)
    assert e.codice == codice
    assert e.categoria is CATALOGO[codice].categoria
    assert e.ammette_retry == CATALOGO[codice].ammette_retry


def test_sollevare_un_codice_con_la_classe_sbagliata_e_respinto():
    """
    **Obiettivo**: Verificare che istanziare manualmente ``ErroreRitentabile("E-S0-01")``
    su un codice di ``REVISIONE_UMANA`` sollevi immediatamente ``TypeError``.

    **Razionale Scientifico/Sistemistico**: Chi solleva un errore nel codice
    dichiara *quale condizione si è verificata* (il codice), ma è solo il
    catalogo a decidere *come va gestita* (la categoria): questo controllo
    impedisce a qualsiasi modulo di forzare un retry su un codice non autorizzato.
    """
    with pytest.raises(TypeError, match="revisione_umana"):
        ErroreRitentabile("E-S0-01")


def test_il_messaggio_contiene_codice_dettaglio_e_azione():
    """
    **Obiettivo**: Verificare che la rappresentazione testuale ``str(e)`` di
    ``ErrorePipeline`` includa il codice, il dettaglio contestuale e l'azione
    prescrittiva del catalogo.

    **Razionale Scientifico/Sistemistico**: Garantisce che sia sulla console
    ``stderr`` sia nei traceback l'operatore legga immediatamente quale campione
    ha causato il problema e quale azione correttiva adottare.
    """
    e = errore("E-S2-01", "campione X1 azzerato")
    testo = str(e)
    assert "E-S2-01" in testo
    assert "campione X1 azzerato" in testo
    assert CATALOGO["E-S2-01"].azione in testo


def test_l_errore_si_traduce_in_evento_strutturato():
    """
    **Obiettivo**: Verificare che ``e.come_evento()`` serializzi l'eccezione in
    un dizionario JSON-compatibile comprensivo dei metadati contestuali
    (``lotto=7``, ``memoria_mb=32000``).

    **Razionale Scientifico/Sistemistico**: Fornisce il payload strutturato per
    ``registra_errore()`` in ``99_logs/pipeline.jsonl`` e per i manifesti di fase.
    """
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
    """
    **Obiettivo**: Verificare che ogni controllo del Gate G15 (`CONTROLLI`) sia
    registrato in ``CATALOGO`` con la medesima descrizione e categoria.

    **Razionale Scientifico/Sistemistico**: Unifica i codici dei gate `G01–G15`
    e i codici delle fasi `S1–S14` in un'unica fonte di verità.
    """
    assert controllo.codice in CATALOGO
    assert controllo.descrizione == CATALOGO[controllo.codice].sintesi
    assert controllo.categoria is CATALOGO[controllo.codice].categoria


@pytest.mark.parametrize("codice", CODICI_DEL_PONTE + CODICI_DEL_GRAFO)
def test_i_codici_del_ponte_e_del_grafo_chiedono_revisione_umana(codice):
    """
    **Obiettivo**: Verificare che tutti i codici del ponte R (``E-R-01`` ..
    ``E-R-04``) e del grafo (``E-GRAFO-01``) abbiano ``Categoria.REVISIONE_UMANA``.

    **Razionale Scientifico/Sistemistico**: Se l'interprete ``Rscript`` non è
    installato (`E-R-01`), se il processo R muore per segfault o timeout senza
    dichiarare l'esito (`E-R-02`), se solleva un bug R imprevisto (`E-R-03`),
    se esaurisce la memoria in una fase che non prevede riduzione del lotto
    (`E-R-04`), o se una fase viene invocata violando l'ordine delle dipendenze
    (`E-GRAFO-01`), ritentare automaticamente sarebbe inutile o dannoso: l'esecutore
    deve fermarsi subito ed esigere l'intervento umano.
    """
    assert CATALOGO[codice].categoria is Categoria.REVISIONE_UMANA


def test_nessun_controllo_di_configurazione_e_ritentabile():
    """
    **Obiettivo**: Verificare che tutti i controlli di configurazione in
    ``CONTROLLI`` (Gate G15) appartengano a ``Categoria.REVISIONE_UMANA``.

    **Razionale Scientifico/Sistemistico**: Un errore formale o logico nel file
    YAML dell'utente è deterministico: rilanciare la pipeline con la stessa
    configurazione invalida fallirebbe all'infinito.
    """
    for controllo in CONTROLLI:
        assert controllo.categoria is Categoria.REVISIONE_UMANA


# --------------------------------------------------------------------------- #
# Proprietà delle categorie                                                    #
# --------------------------------------------------------------------------- #


#: Elenco chiuso dei codici autorizzati alla degradazione automatica controllata.
DEGRADAZIONI = ("E-S0-15", "E-S1-01", "E-S6-02", "E-S11-02")


def test_le_degradazioni_sono_quelle_previste():
    """
    **Obiettivo**: Verificare che i codici con ``Categoria.DEGRADAZIONE_AUTOMATICA``
    siano esattamente ``{"E-S0-15", "E-S1-01", "E-S6-02", "E-S11-02"}``.

    **Razionale Scientifico/Sistemistico**: Limita i comportamenti di ripiego
    non bloccanti ai soli 4 casi scientificamente previsti dal protocollo (piastra
    con pochi blank che passa a `decontam` globale in G08/S12, avviso di scarto
    `truncLen` in G09, warning chimere intermedio in S6, e controllo positivo
    con composizione deviante annotata in S11).
    """
    degradano = {
        c for c, v in CATALOGO.items()
        if v.categoria is Categoria.DEGRADAZIONE_AUTOMATICA
    }
    assert degradano == set(DEGRADAZIONI)


def test_la_guardia_sulla_filogenesi_ferma_l_esecuzione():
    """
    **Obiettivo**: Verificare che ``E-S9-01`` (superamento di ``phylo.max_seqs``
    con filogenesi abilitata) appartenga a ``Categoria.REVISIONE_UMANA`` e non a
    ``DEGRADAZIONE_AUTOMATICA``.

    **Razionale Scientifico/Sistemistico**: Quando l'utente abilita esplicitamente
    ``phylo.enabled = True`` per calcolare metriche UniFrac, saltare silenziosamente
    la costruzione dell'albero in S9 per l'elevato numero di ASV consegnerebbe
    in ``10_phyloseq`` un oggetto privo dell'albero filogenetico richiesto: la
    pipeline deve invece fermarsi e chiedere all'operatore se alzare ``phylo.max_seqs``
    o disattivare ``phylo.enabled``.
    """
    v = CATALOGO["E-S9-01"]
    assert v.categoria is Categoria.REVISIONE_UMANA
    assert "phylo.max_seqs" in v.azione
    assert "phylo.enabled" in v.azione


def test_le_categorie_sono_quattro():
    """
    **Obiettivo**: Verificare che l'enumerazione ``Categoria`` contenga
    esattamente 4 stati di gestione.

    **Razionale Scientifico/Sistemistico**: Mantiene chiusa e ortogonale la
    tassonomia decisionale dell'esecutore della pipeline.
    """
    assert len(Categoria) == 4


def test_ogni_categoria_e_usata_almeno_una_volta():
    """
    **Obiettivo**: Verificare che ciascuna delle 4 categorie di ``Categoria``
    sia assegnata ad almeno un codice in ``CATALOGO``.

    **Razionale Scientifico/Sistemistico**: Assicura che non esistano rami morti
    nella logica di gestione delle eccezioni.
    """
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
    """
    **Obiettivo**: Verificare la tabella di verità delle proprietà booleane
    ``ammette_retry`` e ``ferma_esecuzione`` per ciascuna ``Categoria``.

    **Razionale Scientifico/Sistemistico**: ``PoliticaRetry`` ed ``Esecutore``
    decidono se ritentare e se arrestare la pipeline interrogando queste due
    proprietà; questo test certifica il comportamento esatto delle 4 categorie.
    """
    assert categoria.ammette_retry is retry
    assert categoria.ferma_esecuzione is ferma
