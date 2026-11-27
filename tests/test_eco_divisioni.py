"""Test del catalogo delle divisioni dei metadati (``amplicon16s_eco.divisioni``).

Una divisione e' un modo noto di dividere i campioni in gruppi: ``eco init``
propone quelle che funzionano sul dataset, e ``eco validate`` e ``eco run``
ricavano la regola dal nome scritto nella configurazione. I test coprono il
catalogo e la sua verifica, il riconoscimento delle colonne, le regole, i
criteri, i doppioni, la lettura dei nomi nella configurazione (con i rifiuti
E-ECO-47 ed E-ECO-48) e la compatibilita' con le configurazioni della
versione 0.2.0; quelli con R confrontano un'esecuzione con la divisione scelta
per nome e una con il gruppo derivato scritto a mano.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest
import yaml

from amplicon16s_eco import config as schema
from amplicon16s_eco import divisioni as d
from amplicon16s_eco.catalogo import CATALOGO, ErroreEco, Specie

from eco_aiuti import RADICE, configurazione, costruisci, eco_r  # noqa: F401


def _tabella(**colonne: list[str | None]) -> d.Tabella:
    n = len(next(iter(colonne.values())))
    return d.Tabella(tuple(f"c{i:02d}" for i in range(n)),
                     {nome: tuple(v) for nome, v in colonne.items()})


# --------------------------------------------------------------------------- #
# 1. Il catalogo                                                               #
# --------------------------------------------------------------------------- #


def test_il_catalogo_del_pacchetto_e_valido_e_ogni_voce_ha_la_sua_origine():
    """
    **Obiettivo**: il catalogo del pacchetto si carica; i nomi sono unici; ogni
    voce ha una spiegazione e un'origine che nomina un dataset scaricato (un
    file ``s_OSD-``), una fonte pubblica (la guida o l'API di OSDR),
    l'inventario della pipeline o, per un incrocio, le due divisioni; ogni
    fonte nomina la colonna con la sua chiave; l'impronta identifica il file.

    **Razionale scientifico e sistemistico**: nessun nome di colonna e'
    inventato: chi legge il catalogo sa da dove viene ogni voce.
    """
    catalogo = d.catalogo()
    nomi = [v.nome for v in catalogo]
    assert len(nomi) == len(set(nomi)) >= 20
    for voce in catalogo:
        assert voce.spiegazione and voce.origine
        assert any(f in voce.origine for f in ("s_OSD-", "guida OSDR", "API OSDR",
                                                  "inventario della pipeline",
                                                  "incrocio di due divisioni")), voce.nome
        for fonte in voce.fonti:
            assert d.chiave(fonte.colonna) == fonte.colonna
    assert {"posizione", "modulo", "fase", "volo o terra", "soggetto", "tipo di materiale",
            "tempo di raccolta", "soggetto per fase"} <= set(nomi)
    assert d.impronta_del_catalogo().startswith("sha256:")
    # Il pacchetto installato porta il file del catalogo.
    pyproject = (RADICE / "pyproject.toml").read_text(encoding="utf-8")
    assert '"divisioni.yaml"' in pyproject


@pytest.mark.parametrize("voce, errore", [
    ({"nome": "a", "spiegazione": "x", "origine": "y",
      "fonti": [{"colonna": "a", "regola": "sconosciuta"}]}, "regola"),
    ({"nome": "a", "spiegazione": "x", "origine": "y",
      "fonti": [{"colonna": "Sample Location", "regola": "colonna"}]}, "chiave"),
    ({"nome": "a", "spiegazione": "x", "fonti": [{"colonna": "a", "regola": "colonna"}]},
     "origine"),
    ({"nome": "a_prefisso", "spiegazione": "x", "origine": "y",
      "fonti": [{"colonna": "a", "regola": "colonna"}]}, "suffisso"),
    ({"nome": "a", "spiegazione": "x", "origine": "y",
      "fonti": [{"colonna": "a", "regola": "prefisso", "espressione": "^[A-Z]+"}]},
     "gruppo di cattura"),
    ({"nome": "a", "spiegazione": "x", "origine": "y",
      "fonti": [{"colonna": "a", "regola": "mappatura", "gruppi": {"g": ["x"]}}]},
     "due gruppi"),
    ({"nome": "a", "spiegazione": "x", "origine": "y",
      "fonti": [{"colonna": "a", "regola": "fasce", "limiti": [5, 1], "etichette": [1, 2, 3]}]},
     "crescenti"),
    ({"nome": "a", "spiegazione": "x", "origine": "y",
      "fonti": [{"colonna": "a", "regola": "fasce", "limiti": [1], "etichette": ["x"]}]},
     "etichette"),
    ({"nome": "a", "spiegazione": "x", "origine": "y", "regola": "incrocio",
      "righe": "b", "colonne": "c"}, "incrocio"),
    ({"nome": "a", "spiegazione": "x", "origine": "y", "sconosciuta": 1,
      "fonti": [{"colonna": "a", "regola": "colonna"}]}, "sconosciute"),
])
def test_un_catalogo_malformato_e_respinto(voce, errore):
    """
    **Obiettivo**: il caricamento respinge una regola sconosciuta, una colonna
    che non e' una chiave, un'origine mancante, un nome con il suffisso delle
    divisioni generiche, un prefisso senza gruppo di cattura, una mappatura con
    un solo gruppo, fasce con limiti non crescenti o con etichette sbagliate,
    un incrocio di divisioni che non esistono, una chiave sconosciuta.

    **Razionale scientifico e sistemistico**: il catalogo e' un dato del
    pacchetto; un errore di scrittura deve fermarsi al caricamento, non
    diventare una divisione che divide male.
    """
    with pytest.raises(ValueError, match=errore):
        d.carica_catalogo(yaml.safe_dump({"versione": 1, "divisioni": [voce]}))


def test_un_nome_ripetuto_nel_catalogo_e_respinto():
    """
    **Obiettivo**: due voci con lo stesso nome (maiuscole e spazi non contano)
    sono respinte.

    **Razionale scientifico e sistemistico**: un nome della configurazione
    deve indicare una sola regola.
    """
    voce = {"spiegazione": "x", "origine": "y", "fonti": [{"colonna": "a", "regola": "colonna"}]}
    with pytest.raises(ValueError, match="ripetuto"):
        d.carica_catalogo(yaml.safe_dump({"divisioni": [{"nome": "Fase", **voce},
                                                         {"nome": "fase", **voce}]}))


@pytest.mark.parametrize("nome, chiave", [
    ("Factor Value[Sample Location]", "sample_location"),
    ("factor_value_sample_location", "sample_location"),
    ("FACTOR VALUE [ Sample-Location ]", "sample_location"),
    ("Characteristics[Host Organism]", "host_organism"),
    ("characteristics_subject", "subject"),
    ("Comment[Sample Description]", "sample_description"),
    ("comment_sample_description", "sample_description"),
    ("Parameter Value[Sample Preservation Method]", "sample_preservation_method"),
    ("modulo", "modulo"),
    ("factor_value_", "factor_value"),
])
def test_le_colonne_si_riconoscono_senza_maiuscole_spazi_parentesi_e_prefissi(nome, chiave):
    """
    **Obiettivo**: la chiave di una colonna e' la stessa nella forma ISA-Tab
    (``Factor Value[...]``), in quella dell'oggetto (``factor_value_...``) e
    con maiuscole, spazi, parentesi o trattini diversi; un nome fatto del solo
    prefisso resta com'e'.

    **Razionale scientifico e sistemistico**: lo stesso metadato arriva
    scritto in modi diversi dai file ISA e dalla pipeline.
    """
    assert d.chiave(nome) == chiave


# --------------------------------------------------------------------------- #
# 2. Le regole e i criteri                                                     #
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize("testo, giorni", [
    ("7 days", 7), ("1 day", 1), ("minus10 days", -10), ("-3 days", -3),
    ("520 days plus 3 month", 610), ("520 days plus 6 months", 700), ("2 weeks", 14),
    ("1 year", 365), ("Pre-flight", None), ("", None), ("7 days dopo", None),
])
def test_le_durate_in_parole_si_leggono_in_giorni(testo, giorni):
    """
    **Obiettivo**: una durata in parole («7 days», «minus10 days», «520 days
    plus 3 month») si legge in giorni; un testo che non e' una durata no.

    **Razionale scientifico e sistemistico**: e' la forma dei tempi di
    raccolta nei metadati OSDR (OSD-191); l'ordine alfabetico metterebbe
    «120 days» prima di «14 days».
    """
    assert d.durata_in_giorni(testo) == giorni


def test_le_regole_del_catalogo_dividono_i_valori_come_dichiarato():
    """
    **Obiettivo**: il prefisso riunisce le posizioni del modulo; la mappatura
    ordinata riconosce le fasi e mette il numero d'ordine davanti ai nomi;
    l'ordine di tempo numera i valori dal primo; le fasce dividono le durate
    ai limiti; un valore mancante o non riconosciuto resta senza gruppo; il
    gruppo derivato prodotto elenca solo i valori presenti.

    **Razionale scientifico e sistemistico**: la regola scritta nel catalogo
    e' quella che le analisi applicano.
    """
    tabella = _tabella(
        factor_value_sample_location=["NOD1A1", "NOD1B2", "LAB1C3", "LAB1D1", "Air Sample",
                                      None],
        comment_sample_description=["samples collected before entering the module"] * 2
        + ["samples collected while inside the module"] * 2
        + ["samples collected after exiting the module", "NA"],
        factor_value_time=["minus10 days", "1 day", "14 days", "120 days",
                           "520 days plus 3 month", "not applicable"],
    )
    per_nome = {v.nome: v for v in d.catalogo()}
    modulo, _ = d.applica(per_nome["modulo"], tabella)
    assert modulo.gruppi == (("LAB1", 2), ("NOD1", 2)) and modulo.senza_gruppo == 2
    assert modulo.espansione == {"gruppo": {"from": "factor_value_sample_location",
                                            "pattern": "^([A-Za-z]+[0-9]+)"}}
    fase, motivi = d.applica(per_nome["fase"], tabella)
    assert motivi and "factor_value_time" in motivi[0]
    assert fase.colonna == "comment_sample_description"
    assert fase.gruppi == (("1 prima", 2), ("2 durante", 2), ("3 dopo", 1))
    assert fase.espansione["gruppo"]["groups"]["1 prima"] == [
        "samples collected before entering the module"]
    tempo, _ = d.applica(per_nome["tempo di raccolta"], tabella)
    assert [g for g, _ in tempo.gruppi] == [
        "1 minus10 days", "2 1 day", "3 14 days", "4 120 days", "5 520 days plus 3 month"]
    fasce, _ = d.applica(per_nome["fascia di tempo"], tabella)
    assert fasce.gruppi_dei_campioni == (
        "1 prima", "2 0-29 giorni", "2 0-29 giorni", "4 90-179 giorni", "6 da 365 giorni", None)
    date = _tabella(characteristics_collection_date=["02-Nov-2017", "30-Oct-2017",
                                                     "01-Nov-2017", "30-Oct-2017"])
    raccolta, _ = d.applica(per_nome["tempo di raccolta"], date)
    assert [g for g, _ in raccolta.gruppi] == ["1 30-Oct-2017", "2 01-Nov-2017", "3 02-Nov-2017"]
    # Valori che non sono tutti tempi: la fonte non vale e si passa alla
    # successiva, le date di raccolta; senza un'altra fonte la divisione non c'e'.
    testo = ["Pre-flight", "Post-flight", "Pre-flight", "Post-flight"]
    poi_date = _tabella(factor_value_time=testo, characteristics_collection_date=[
        "31-Oct-2017", "02-Nov-2017", "30-Oct-2017", "02-Nov-2017"])
    raccolta, motivi = d.applica(per_nome["tempo di raccolta"], poi_date)
    assert raccolta.colonna == "characteristics_collection_date"
    assert [g for g, _ in raccolta.gruppi] == ["1 30-Oct-2017", "2 31-Oct-2017", "3 02-Nov-2017"]
    assert len(motivi) == 1 and "non sono tutti durate o date" in motivi[0]
    assert "'Post-flight'" in motivi[0]
    solo_testo, motivi = d.applica(per_nome["tempo di raccolta"], _tabella(factor_value_time=testo))
    assert solo_testo is None and "non sono tutti durate o date" in motivi[0]
    misti = _tabella(factor_value_time=["7 days", "14 days", "Pre-flight", "7 days"])
    assert d.applica(per_nome["tempo di raccolta"], misti)[0] is None


def test_volo_o_terra_riconosce_i_valori_documentati_e_ignora_gli_altri():
    """
    **Obiettivo**: «Space Flight», «Flight» e «in-flight» sono il volo,
    «Ground Control» e «basal control» la terra, senza distinguere le
    maiuscole; «Not Applicable» resta senza gruppo.

    **Razionale scientifico e sistemistico**: sono i valori dei dataset
    scaricati e della documentazione OSDR citati nel catalogo.
    """
    tabella = _tabella(factor_value_spaceflight=["Space Flight", "space flight", "Flight",
                                                 "in-flight", "Ground Control",
                                                 "basal control", "Not Applicable"])
    app, _ = d.applica(next(v for v in d.catalogo() if v.nome == "volo o terra"), tabella)
    assert app.gruppi_dei_campioni == ("volo",) * 4 + ("terra",) * 2 + (None,)
    # Un valore che due gruppi riconoscono va al primo, nell'ordine del catalogo.
    fasi = _tabella(factor_value_time=["Pre-flight, then in-flight", "Post-flight"])
    app, _ = d.applica(next(v for v in d.catalogo() if v.nome == "fase"), fasi)
    assert app.gruppi_dei_campioni == ("1 prima", "3 dopo")


@pytest.mark.parametrize("valori, dimensione, motivo", [
    (["a"] * 10, 2, "1 gruppi, ne servono almeno 2"),
    ([f"v{i}" for i in range(13)] * 2, 2, "13 gruppi, piu' di 12"),
    (["a", "a", "b", "b"] + [None] * 6, 2, "meno della meta'"),
    (["a"] * 6 + ["b"] * 4, 5, "1 gruppi con almeno 5 campioni"),
    (["a"] * 5 + ["b"] * 5, 5, None),
])
def test_i_criteri_dicono_quando_una_divisione_funziona(valori, dimensione, motivo):
    """
    **Obiettivo**: una divisione funziona con da 2 a 12 gruppi, un gruppo per
    almeno meta' dei campioni e almeno due gruppi con la dimensione richiesta;
    altrimenti il motivo lo dice.

    **Razionale scientifico e sistemistico**: una divisione con un solo
    gruppo, con un gruppo per campione o che lascia fuori quasi tutti i
    campioni non sostiene un confronto.
    """
    tabella = _tabella(colonna=valori)
    app = d._generica_per_colonna("colonna", tabella)
    trovato = d.giudica(app, len(valori), dimensione)
    assert (trovato is None) if motivo is None else (motivo in trovato)


def test_le_proposte_tolgono_i_doppioni_e_sono_le_stesse_a_ogni_esecuzione():
    """
    **Obiettivo**: fra due divisioni che dividono i campioni negli stessi
    gruppi resta la prima (il catalogo prima delle generiche); le proposte,
    nello stesso ordine, sono le stesse a ogni chiamata; un incrocio si
    propone se lo sono le sue due divisioni; una divisione del catalogo il cui
    nome e' quello di una colonna con valori, ma che userebbe un'altra regola,
    non si propone.

    **Razionale scientifico e sistemistico**: due nomi per la stessa
    partizione confondono chi sceglie; e risultati riproducibili chiedono
    proposte riproducibili.
    """
    soggetti = ["S1", "S2", "S3"] * 6
    fasi = (["samples collected before entering the module"] * 6
            + ["samples collected while inside the module"] * 6
            + ["samples collected after exiting the module"] * 6)
    tabella = _tabella(characteristics_subject=soggetti, comment_sample_description=fasi,
                       copia_del_soggetto=[s.lower() for s in soggetti])
    selezionabili = ["characteristics_subject", "comment_sample_description",
                     "copia_del_soggetto"]
    prima = d.proponi(tabella, selezionabili)
    assert prima == d.proponi(tabella, selezionabili)
    assert [p.nome for p in prima.proposte] == ["fase", "soggetto", "soggetto per fase"]
    assert ("descrizione del campione", "fase") in prima.doppioni
    assert ("copia_del_soggetto", "soggetto") in prima.doppioni
    griglia = prima.proposte[-1].applicazione
    assert griglia.espansione == {"griglia": ("soggetto", "fase")}
    assert ("S1 × 1 prima", 2) in griglia.gruppi
    # Un incrocio di due divisioni annidate (ogni fase in una sola quota) non
    # si propone: la griglia non le incrocia.
    annidate = _tabella(factor_value_altitude=["Ground"] * 6 + ["Atmosphere"] * 6
                        + ["Ground"] * 6,
                        factor_value_time=["Pre-flight"] * 6 + ["Cruise, in-flight"] * 6
                        + ["Post-flight"] * 6)
    proposte = d.proponi(annidate, [])
    assert [p.nome for p in proposte.proposte] == ["fase", "quota"]
    assert any(n == "quota per fase" and "annidate" in m for n, m in proposte.scartate)
    # Una colonna di tempi che non sono durate ne' date non e' «tempo di
    # raccolta»: si propone la colonna stessa, come divisione generica.
    tempi = ["Pre-flight", "Post-flight", "Cruise, in-flight", "Ascent/Descent, in-flight"] * 5
    proposte = d.proponi(_tabella(factor_value_time=tempi), ["factor_value_time"])
    assert [(p.nome, p.origine) for p in proposte.proposte] == [
        ("fase", "catalogo"), ("factor_value_time", "generica")]
    assert any(n == "tempo di raccolta" and "durate o date" in m for n, m in proposte.scartate)
    # Una colonna «soggetto» con valori vince sul nome: la divisione non si propone.
    con_colonna = _tabella(characteristics_subject=soggetti, soggetto=["x", "y"] * 9)
    proposte = d.proponi(con_colonna, ["soggetto"])
    assert "soggetto" in [p.nome for p in proposte.proposte]
    assert next(p for p in proposte.proposte if p.nome == "soggetto").origine == "generica"
    assert any(n == "soggetto" and "vince" in m for n, m in proposte.scartate)


# --------------------------------------------------------------------------- #
# 3. I nomi nella configurazione                                               #
# --------------------------------------------------------------------------- #


def _dati(**modifica) -> dict:
    return schema.valida(configurazione(**modifica)).model_dump(mode="json", by_alias=True)


_POSIZIONI = ["NOD1A1", "NOD1A1", "NOD1B2", "NOD1B2", "LAB1C3", "LAB1C3", "LAB1D1", "LAB1D1"]


def test_una_divisione_per_nome_diventa_un_gruppo_derivato_o_una_colonna():
    """
    **Obiettivo**: una divisione con una regola diventa un gruppo derivato con
    il suo nome, aggiunto dopo quelli dichiarati; una divisione «colonna»
    diventa la colonna in ogni campo in cui e' nominata, anche nel modello
    vincolato; la configurazione espansa e' valida; che cosa e' diventata ogni
    divisione si registra con i campi.

    **Razionale scientifico e sistemistico**: le analisi conoscono gia'
    colonne e gruppi derivati: la divisione e' un modo di scriverli, non un
    calcolo nuovo.
    """
    tabella = _tabella(gruppo=["a", "b"] * 4, factor_value_sample_location=_POSIZIONI,
                       modulo=[None] * 8)
    dati = _dati(design={"variable": "modulo", "other_variables": ["posizione"],
                         "groupings": {"altro": {"from": "gruppo", "pattern": "^(.)"}}},
                 comp={"group_by": ["modulo", "gruppo"], "facet_wrap": "posizione"},
                 ord={"color_by": ["gruppo"]},
                 constrained={"method": "dbrda", "variables": ["modulo", "posizione"],
                              "distances": ["bray"], "max_df_fraction": 0.5})
    espansione = d.espandi(dati, tabella)
    assert espansione.cambiata and not espansione.rifiuti
    design = espansione.dati["design"]
    assert design["variable"] == "modulo"
    assert design["other_variables"] == ["factor_value_sample_location"]
    assert list(design["groupings"]) == ["altro", "modulo"]
    assert design["groupings"]["modulo"] == {"from": "factor_value_sample_location",
                                             "pattern": "^([A-Za-z]+[0-9]+)"}
    assert espansione.dati["comp"]["facet_wrap"] == "factor_value_sample_location"
    assert espansione.dati["constrained"]["variables"] == ["modulo",
                                                           "factor_value_sample_location"]
    schema.valida(espansione.dati)
    scelte = {s["nome"]: s for s in espansione.scelte}
    assert scelte["modulo"]["campi"] == ["design.variable", "comp.group_by",
                                         "constrained.variables"]
    assert scelte["posizione"]["campi"] == ["design.other_variables", "comp.facet_wrap",
                                            "constrained.variables"]
    assert scelte["modulo"]["diventa"] == "il gruppo derivato modulo"
    assert scelte["posizione"]["diventa"] == "la colonna factor_value_sample_location"


def test_un_incrocio_in_facet_grid_diventa_righe_e_colonne():
    """
    **Obiettivo**: il nome di un incrocio come valore di ``comp.facet_grid``
    diventa ``{rows, cols}`` con le due divisioni espanse; uno fuori da
    ``comp.facet_grid``, o un nome che non e' un incrocio come valore di
    ``comp.facet_grid``, e' ``E-ECO-47``; lo schema accetta il nome.

    **Razionale scientifico e sistemistico**: una griglia ha due dimensioni, e
    un incrocio ne nomina due in una volta.
    """
    tabella = _tabella(
        gruppo=["a", "b"] * 9, characteristics_subject=["S1", "S2", "S3"] * 6,
        comment_sample_description=["before entering"] * 6 + ["while inside"] * 6
        + ["after exiting"] * 6)
    dati = _dati(comp={"facet_grid": "soggetto per fase"})
    assert schema.valida(configurazione(comp={"facet_grid": "soggetto per fase"}))
    espansione = d.espandi(dati, tabella)
    assert espansione.dati["comp"]["facet_grid"] == {"rows": "characteristics_subject",
                                                     "cols": "fase"}
    assert list(espansione.dati["design"]["groupings"]) == ["fase"]
    schema.valida(espansione.dati)
    fuori = d.espandi(_dati(design={"variable": "soggetto per fase"}), tabella)
    assert [c for c, _ in fuori.rifiuti] == ["E-ECO-47"] and "soggetto per fase" in fuori.respinti
    sbagliata = d.espandi(_dati(comp={"facet_grid": "fase"}), tabella)
    assert [c for c, _ in sbagliata.rifiuti] == ["E-ECO-47"]
    assert sbagliata.dati["comp"]["facet_grid"] is None


def test_la_lettura_dei_nomi_segue_una_regola_univoca():
    """
    **Obiettivo**: un gruppo derivato dichiarato vince su una divisione con lo
    stesso nome; una colonna con valori vince su una divisione; una colonna
    vuota cede il nome alla divisione; ``<colonna>_prefisso`` e' la divisione
    generica per prefisso; un nome che non e' nulla di noto resta per le
    analisi, che lo respingono come colonna mancante, e la traduzione ne fa un
    ``E-ECO-47`` con le colonne e le divisioni.

    **Razionale scientifico e sistemistico**: un nome deve indicare una cosa
    sola, e una configurazione della versione 0.2.0, che nomina colonne,
    deve leggersi come prima.
    """
    tabella = _tabella(gruppo=["a", "b"] * 4, factor_value_sample_location=_POSIZIONI,
                       modulo=[None] * 8, posizione=["p", "q"] * 4)
    dichiarato = d.espandi(_dati(design={"variable": "modulo", "groupings": {
        "modulo": {"from": "gruppo", "groups": {"x": ["a"], "y": ["b"]}}}}), tabella)
    assert not dichiarato.cambiata
    colonna = d.espandi(_dati(design={"variable": "posizione"}), tabella)
    assert not colonna.cambiata and colonna.dati["design"]["variable"] == "posizione"
    vuota = d.espandi(_dati(design={"variable": "modulo"}), tabella)
    assert vuota.dati["design"]["groupings"]["modulo"]["pattern"] == d.PREFISSO
    generica = d.espandi(_dati(design={"variable": "factor_value_sample_location_prefisso"}),
                         tabella)
    assert generica.dati["design"]["groupings"] == {"factor_value_sample_location_prefisso": {
        "from": "factor_value_sample_location", "pattern": d.PREFISSO}}
    ignoto = d.espandi(_dati(design={"variable": "grupo"}), tabella)
    assert not ignoto.cambiata and not ignoto.rifiuti
    tradotti = d.traduci([("E-ECO-05", "Colonna 'grupo' non trovata. Colonne dell'oggetto: "
                                       "gruppo, modulo."),
                          ("E-ECO-05", "Colonna 'lotto' non trovata. Colonne dell'oggetto: "
                                       "gruppo, modulo.")], ignoto)
    assert [c for c, _ in tradotti] == ["E-ECO-47", "E-ECO-05"]
    assert "design.variable: 'grupo'" in tradotti[0][1] and "volo o terra" in tradotti[0][1]


def test_una_divisione_senza_gruppi_validi_e_respinta_con_il_motivo():
    """
    **Obiettivo**: una divisione che sul dataset non da' gruppi validi (nessuna
    colonna a cui applicarla, troppo pochi gruppi, gruppi sotto
    ``stat.min_group_size``) e' ``E-ECO-48``, con il motivo e il rimedio; il
    nome respinto non si ripete come colonna mancante.

    **Razionale scientifico e sistemistico**: la regola di una divisione vale
    su un dataset solo se i gruppi che produce sostengono i test; dirlo prima
    del calcolo evita un'analisi senza confronti.
    """
    tabella = _tabella(gruppo=["a", "b"] * 4, factor_value_sample_location=_POSIZIONI)
    nessuna = d.espandi(_dati(design={"variable": "soggetto"}), tabella)
    assert [c for c, _ in nessuna.rifiuti] == ["E-ECO-48"]
    assert "nessuna colonna" in nessuna.rifiuti[0][1] and "Rimedio" in nessuna.rifiuti[0][1]
    piccoli = d.espandi(_dati(design={"variable": "modulo"}, stat={"min_group_size": 5}),
                        tabella)
    assert [c for c, _ in piccoli.rifiuti] == ["E-ECO-48"]
    assert "con almeno 5 campioni" in piccoli.rifiuti[0][1]
    assert d.traduci([("E-ECO-05", "Colonna 'modulo' non trovata. Colonne dell'oggetto: "
                                   "gruppo.")], piccoli) == []
    # Una divisione in ordine di tempo su valori che non sono tempi: E-ECO-48
    # con il motivo, e la colonna resta nominabile cosi' com'e'.
    fasi = _tabella(gruppo=["a", "b"] * 4,
                    factor_value_time=["Pre-flight", "Post-flight"] * 4)
    tempo = d.espandi(_dati(design={"variable": "tempo di raccolta"}), fasi)
    assert [c for c, _ in tempo.rifiuti] == ["E-ECO-48"]
    assert "non sono tutti durate o date" in tempo.rifiuti[0][1]
    assert "tempo di raccolta" in tempo.respinti
    for codice in ("E-ECO-47", "E-ECO-48"):
        assert CATALOGO[codice].specie is Specie.RIFIUTO


def test_le_configurazioni_dei_dataset_non_cambiano_con_l_espansione():
    """
    **Obiettivo**: ogni ``dati/*/eco_*.yaml`` (configurazioni della versione
    0.2.0) letto su una tabella che ha le sue colonne resta identico: nessun
    nome e' una divisione, e le analisi ricevono gli stessi parametri e lo
    stesso digest.

    **Razionale scientifico e sistemistico**: il catalogo aggiunge un modo di
    scrivere i gruppi senza cambiare il significato dei nomi gia' in uso.
    """
    esempi = sorted((RADICE / "dati").glob("*/eco_*.yaml"))
    if not esempi:
        pytest.skip("la cartella dati/ non e' disponibile in questo ambiente")
    for esempio in esempi:
        config = schema.carica(esempio)
        dati = config.model_dump(mode="json", by_alias=True)
        nomi = {n for _, n in d._nomi(dati)} | {
            g["from"] for g in (dati["design"].get("groupings") or {}).values()}
        tabella = _tabella(**{n: ["x", "y"] * 5 for n in nomi})
        espansione = d.espandi(dati, tabella)
        assert not espansione.cambiata and espansione.dati is dati, esempio.name
        assert schema.digest(schema.valida(espansione.dati)) == schema.digest(config)


# --------------------------------------------------------------------------- #
# 4. Con R: la divisione per nome e il gruppo scritto a mano                   #
# --------------------------------------------------------------------------- #

_CONTEGGI = [[10 + (7 * i + 3 * j) % 23 for j in range(12)] for i in range(6)]
_METADATI = {
    "factor_value_sample_location": ["NOD1A1", "NOD1A2", "NOD1B1", "NOD1B2", "LAB1A1",
                                     "LAB1A2", "LAB1B1", "LAB1B2", "COL1A1", "COL1A2",
                                     "Air Sample", "Air Sample"],
    "modulo": [None] * 12,
    "gruppo": ["a", "b"] * 6,
    "characteristics_subject": ["S1", "S2"] * 6,
    "comment_sample_description": ["before entering the module"] * 4
    + ["while inside the module"] * 4 + ["after exiting the module"] * 4,
}


def test_la_divisione_per_nome_da_gli_stessi_file_del_gruppo_scritto_a_mano(eco_r, tmp_path):
    """
    **Obiettivo**: l'analisi con ``design.variable: modulo`` (la divisione del
    catalogo, mentre la colonna ``modulo`` dell'oggetto e' vuota) e quella con
    il gruppo derivato ``modulo`` scritto a mano con la stessa regola danno
    file identici byte per byte, salvo ``configurazione.json``, che nella prima
    ha in piu' il campo descrittivo ``divisioni`` (stessi parametri e stesso
    digest), e il manifesto, che ne riporta l'impronta. Lo stesso vale per la
    griglia «soggetto per fase» in validazione.

    **Razionale scientifico e sistemistico**: scegliere una divisione per nome
    e' solo un modo di scrivere la regola; i risultati devono essere gli
    stessi.
    """
    from eco_aiuti import impronte

    from amplicon16s_eco.esecuzione import esegui, valida_sull_oggetto

    oggetto = costruisci(tmp_path, _CONTEGGI, _METADATI)
    per_nome = configurazione(design={"variable": "modulo"}, comp={"group_by": ["modulo"]})
    a_mano = configurazione(design={"variable": "modulo", "groupings": {"modulo": {
        "from": "factor_value_sample_location", "pattern": d.PREFISSO}}},
        comp={"group_by": ["modulo"]})
    uscite = {}
    for nome, dati in (("per_nome", per_nome), ("a_mano", a_mano)):
        esito = esegui(oggetto, schema.valida(dati), tmp_path / nome)
        uscite[nome] = (esito, impronte(tmp_path / nome))
    (esito_n, file_n), (esito_m, file_m) = uscite["per_nome"], uscite["a_mano"]
    assert esito_n.digest_configurazione == esito_m.digest_configurazione
    assert set(file_n) == set(file_m)
    diversi = sorted(f for f in file_n if file_n[f] != file_m[f])
    assert diversi == ["configurazione.json", "manifest.json"]
    conf_n = json.loads((tmp_path / "per_nome" / "configurazione.json").read_text("utf-8"))
    conf_m = json.loads((tmp_path / "a_mano" / "configurazione.json").read_text("utf-8"))
    assert conf_n["parametri"] == conf_m["parametri"] and conf_n["digest"] == conf_m["digest"]
    assert set(conf_n) - set(conf_m) == {"divisioni"}
    (scelta,) = conf_n["divisioni"]["scelte"]
    assert scelta["nome"] == "modulo" and scelta["gruppi"] == {"COL1": 2, "LAB1": 4, "NOD1": 4}
    assert esito_n.riepilogo == esito_m.riepilogo

    griglia_n = valida_sull_oggetto(oggetto, schema.valida(configurazione(
        comp={"facet_grid": "soggetto per fase"})))
    griglia_m = valida_sull_oggetto(oggetto, schema.valida(configurazione(
        design={"variable": "gruppo", "groupings": {"fase": {
            "from": "comment_sample_description", "groups": {
                "1 prima": ["before entering the module"],
                "2 durante": ["while inside the module"],
                "3 dopo": ["after exiting the module"]}}}},
        comp={"facet_grid": {"rows": "characteristics_subject", "cols": "fase"}})))
    assert griglia_n.digest_configurazione == griglia_m.digest_configurazione
    assert griglia_n.riepilogo == griglia_m.riepilogo


def test_validate_respinge_i_nomi_sconosciuti_e_le_divisioni_senza_gruppi(eco_r, tmp_path):
    """
    **Obiettivo**: sull'oggetto, un nome sconosciuto in un campo delle
    divisioni e' ``E-ECO-47`` insieme agli altri rifiuti delle analisi (una
    variabile tecnica mancante resta ``E-ECO-05``); una divisione sotto
    ``stat.min_group_size`` e' ``E-ECO-48``; un gruppo derivato non puo'
    prendere il nome di una colonna con valori (``E-ECO-40``), ma sostituisce
    una colonna vuota.

    **Razionale scientifico e sistemistico**: ogni problema si dice prima del
    calcolo, tutti nello stesso giro, con un codice che dice che cosa
    correggere.
    """
    from amplicon16s_eco.esecuzione import valida_sull_oggetto

    oggetto = costruisci(tmp_path, _CONTEGGI, _METADATI)

    def codici(**modifica):
        with pytest.raises(ErroreEco) as rifiuto:
            valida_sull_oggetto(oggetto, schema.valida(configurazione(**modifica)))
        return rifiuto.value.codici, str(rifiuto.value)

    trovati, testo = codici(design={"variable": "modullo", "technical_variables": ["lotto"]})
    assert trovati == ["E-ECO-47", "E-ECO-05"] and "Divisioni del catalogo" in testo
    trovati, testo = codici(design={"variable": "modulo"}, stat={"min_group_size": 5})
    assert trovati == ["E-ECO-48"] and "con almeno 5 campioni" in testo
    trovati, _ = codici(design={"variable": "gruppo", "groupings": {
        "characteristics_subject": {"from": "gruppo", "pattern": "^(.)"}}})
    assert trovati == ["E-ECO-40"]
    esito = valida_sull_oggetto(oggetto, schema.valida(configurazione(design={
        "variable": "modulo", "groupings": {"modulo": {
            "from": "factor_value_sample_location", "pattern": d.PREFISSO}}})))
    (variabile,) = esito.riepilogo["variabili"]
    assert variabile["nome"] == "modulo" and variabile["gruppi"] == {"COL1": 2, "LAB1": 4,
                                                                     "NOD1": 4}


def test_eco_init_scrive_il_file_dei_gruppi_accanto_alla_configurazione(eco_r, tmp_path):
    """
    **Obiettivo**: ``eco init`` scrive accanto alla configurazione il file dei
    gruppi delle divisioni proposte; nel file di configurazione ci sono il
    numero e l'elenco delle divisioni, non i gruppi; se uno dei due file esiste
    gia', init non scrive nulla.

    **Razionale scientifico e sistemistico**: il file di configurazione resta
    leggibile, e due esecuzioni di init non si sovrascrivono.
    """
    from amplicon16s_eco.cli import USCITA_RIFIUTO, USCITA_SUCCESSO, main

    oggetto = costruisci(tmp_path, _CONTEGGI, _METADATI)
    uscita = tmp_path / "eco.yaml"
    assert main(["init", "--object", str(oggetto), "--output", str(uscita)]) == USCITA_SUCCESSO
    testo = uscita.read_text(encoding="utf-8")
    tsv = (tmp_path / "eco_divisioni_proposte.tsv").read_text(encoding="utf-8")
    assert "# Divisioni suggerite dai metadati: " in testo
    assert "eco_divisioni_proposte.tsv" in " ".join(testo.replace("#", " ").split())
    assert tsv.startswith("divisione\tregola\tgruppo\tcampioni\tcampioni_senza_gruppo\n")
    # Il modulo ha gruppi da 2 e 4 campioni, sotto la soglia delle proposte; il
    # soggetto due da 6.
    assert "\nmodulo\t" not in tsv
    assert "soggetto\tcolonna characteristics_subject cosi' com'e'\tS1\t6\t0\n" in tsv
    assert "#   soggetto" in testo and "S1" not in testo
    (tmp_path / "altro.yaml").unlink(missing_ok=True)
    (tmp_path / "altro_divisioni_proposte.tsv").write_text("", encoding="utf-8")
    assert main(["init", "--object", str(oggetto), "--output",
                 str(tmp_path / "altro.yaml")]) == USCITA_RIFIUTO
    assert not (tmp_path / "altro.yaml").exists()
