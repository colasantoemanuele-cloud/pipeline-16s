"""Gruppi scelti dall'utente: gruppi derivati, colori descrittivi, palette,
proposte di raggruppamento di ``eco init``.

Un gruppo derivato (``design.groupings``) e' una nuova variabile categorica
ricavata da una colonna, con un elenco di valori per gruppo o con
un'espressione regolare: si verifica, su tabelle scritte qui, che ogni
campione vada nel gruppo atteso, che chi resta senza gruppo esca dalle analisi
di quella variabile con il motivo, e che ogni regola sbagliata sia respinta
prima del calcolo. ``ord.color_by`` colora le ordinazioni senza alcun test;
``plot.palette`` sceglie i colori di tutti i grafici.
"""

from __future__ import annotations

import pytest

from amplicon16s_eco import config as schema
from amplicon16s_eco.catalogo import CATALOGO, ErroreEco, Specie
from amplicon16s_eco.esecuzione import esegui, valida_sull_oggetto

from eco_aiuti import (  # noqa: F401
    RADICE,
    configurazione,
    costruisci,
    eco_r,
    esegui_rscript,
    impronte,
    leggi_tsv,
    righe_permanova,
)


def _codici(dati) -> list[str]:
    with pytest.raises(ErroreEco) as rifiuto:
        schema.valida(dati)
    return rifiuto.value.codici


# --------------------------------------------------------------------------- #
# 1. Lo schema                                                                 #
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize("regola", [
    {"from": "luogo", "groups": {"a": ["x"]}, "pattern": "^(.)"},
    {"from": "luogo"},
    {"from": "luogo", "groups": {}},
    {"from": "luogo", "groups": {"a": []}},
    {"from": "luogo", "groups": {"a": ["x", "y"], "b": ["y"]}},
    {"from": "luogo", "pattern": "^([A-Z"},
    {"from": "luogo", "pattern": "^[A-Z]+"},
    {"groups": {"a": ["x"]}},
    {"from": "luogo", "groups": {"a": ["x"]}, "colonna": "x"},
])
def test_una_regola_di_raggruppamento_sbagliata_e_respinta_con_il_suo_codice(regola):
    """
    **Obiettivo**: un gruppo derivato con entrambe le regole o con nessuna, con
    un gruppo vuoto, con un valore in due gruppi, con un'espressione non
    valida o senza gruppo di cattura, senza la colonna di partenza o con una
    chiave sconosciuta e' respinto con ``E-ECO-40``, prima di leggere
    l'oggetto.

    **Razionale scientifico e sistemistico**: una regola ambigua metterebbe un
    campione in un gruppo a caso, o in nessuno, senza che chi legge i
    risultati lo sappia.
    """
    codici = _codici(configurazione(design={"groupings": {"zona": regola}}))
    assert set(codici) == {"E-ECO-40"}
    assert CATALOGO["E-ECO-40"].specie is Specie.RIFIUTO


def test_i_gruppi_derivati_e_i_colori_entrano_nella_configurazione_solo_se_dichiarati():
    """
    **Obiettivo**: senza gruppi derivati ne' colori la configurazione risolta
    non ha quelle chiavi; dichiarati, vi compaiono (``from`` con il suo nome)
    e cambiano il digest. Un gruppo derivato non si ricava da un altro;
    ``ord.color_by`` non ripete una colonna ne' una variabile del disegno; la
    palette e' un nome o un elenco di colori esadecimali (``E-ECO-41``).

    **Razionale scientifico e sistemistico**: il digest identifica l'analisi;
    una figura descrittiva in piu' non deve cambiare quello delle
    configurazioni che non la chiedono.
    """
    base = schema.risolta(schema.valida(configurazione()))
    assert "groupings" not in base["design"] and "color_by" not in base["ord"]
    assert base["plot"] == {"palette": "Dark 3"}
    regola = {"from": "luogo", "pattern": "^([A-Z]+)"}
    dichiarata = schema.valida(configurazione(design={"groupings": {"zona": regola}},
                                              ord={"color_by": ["zona"]}))
    risolta = schema.risolta(dichiarata)
    assert risolta["design"]["groupings"] == {"zona": {**regola, "groups": None}}
    assert risolta["ord"]["color_by"] == ["zona"]
    assert schema.digest(dichiarata) != schema.digest(schema.valida(configurazione()))
    assert _codici(configurazione(design={"groupings": {
        "zona": regola, "macro": {"from": "zona", "pattern": "^(.)"}}})) == ["E-ECO-40"]
    assert _codici(configurazione(ord={"color_by": ["luogo", "luogo"]})) == ["E-ECO-01"]
    assert _codici(configurazione(ord={"color_by": ["gruppo"]})) == ["E-ECO-01"]
    for palette in (["#12345"], ["rosso"], [], 3):
        assert _codici(configurazione(plot={"palette": palette})) == ["E-ECO-41"], palette
    for palette in ("Okabe-Ito", ["#1B9E77", "#D95F02AA"]):
        schema.valida(configurazione(plot={"palette": palette}))


# --------------------------------------------------------------------------- #
# 2. I gruppi derivati sull'oggetto                                            #
# --------------------------------------------------------------------------- #

#: Dieci campioni con un luogo ciascuno; tre varianti, il terzo luogo diverso.
LUOGHI = ["NOD1A1", "NOD1A2", "NOD2B1", "NOD2B2", "LAB1C1", "LAB1C2", "LAB1C3",
          "Aria", "NOD1A3", "Tubo"]
CONTEGGI = [
    [30, 32, 5, 6, 10, 11, 12, 20, 31, 15],
    [5, 6, 30, 31, 10, 12, 11, 20, 4, 15],
    [10, 9, 10, 8, 30, 29, 31, 20, 11, 15],
]


def _oggetto(tmp_path, nome="dati"):
    cartella = tmp_path / nome
    cartella.mkdir(exist_ok=True)
    return costruisci(cartella, CONTEGGI, {"gruppo": ["a", "b"] * 5, "luogo": LUOGHI})


def test_un_gruppo_derivato_per_espressione_assegna_i_campioni_attesi(eco_r, tmp_path):
    """
    **Obiettivo**: con ``pattern: "^([A-Z]+[0-9])"`` i campioni vanno nei
    gruppi NOD1 (3), NOD2 (2), LAB1 (3), calcolati a mano; i valori che
    l'espressione non riconosce restano senza gruppo, escono dalle analisi
    della variabile con il motivo in ``campioni.tsv``, e ``validate`` riporta i
    gruppi ottenuti con il numero di campioni. Il gruppo derivato si usa come
    variabile: PERMANOVA, composizione e grafici lo trattano come una colonna.

    **Razionale scientifico e sistemistico**: la regola si decide nelle analisi
    ecologiche, senza rifare la pipeline; il suo esito va controllato prima
    del calcolo.
    """
    oggetto = _oggetto(tmp_path)
    modulo = {"from": "luogo", "pattern": "^([A-Z]+[0-9])"}
    config = schema.valida(configurazione(
        design={"variable": "modulo", "groupings": {"modulo": modulo}},
        ord={"methods": ["pcoa"], "distances": ["bray"], "pcoa_correction": "none"}))
    esito = valida_sull_oggetto(oggetto, config)
    gruppi = esito.riepilogo["gruppi_derivati"]
    assert gruppi == [{"nome": "modulo", "da": "luogo",
                       "regola": "pattern: ^([A-Z]+[0-9]) (primo gruppo di cattura)",
                       "gruppi": {"LAB1": 3, "NOD1": 3, "NOD2": 2},
                       "campioni_senza_gruppo": 2}]
    uscita = tmp_path / "uscita"
    esito = esegui(oggetto, config, uscita)
    _, campioni = leggi_tsv(uscita / "campioni.tsv")
    per = {r["campione"]: r for r in campioni}
    assert per["c01"]["valore:modulo"] == "NOD1" and per["c05"]["valore:modulo"] == "LAB1"
    assert per["c08"]["valore:modulo"] == "NA"
    assert per["c08"]["stato:modulo"] == (
        "escluso: il valore di 'luogo' non rientra in alcun gruppo di 'modulo'")
    assert per["c10"]["stato:modulo"] == per["c08"]["stato:modulo"]
    assert per["c08"]["stato"] == "analizzato"
    righe = righe_permanova(uscita, variabile="modulo")
    assert righe[("permanova", "modulo")]["campioni"] == "8"
    assert righe[("permanova", "modulo")]["gruppi"] == "3"
    assert (uscita / "composizione" / "gruppi_01_modulo.png").is_file()
    assert (uscita / "ordinazione" / "pcoa_bray_01_modulo.png").is_file()


def test_un_gruppo_derivato_per_elenco_e_le_sue_regole_sull_oggetto(eco_r, tmp_path):
    """
    **Obiettivo**: con ``groups`` ogni valore elencato va al suo gruppo e gli
    altri restano senza; un valore elencato che la colonna non ha, un nome
    che e' gia' una colonna o una colonna di partenza che non esiste sono
    respinti dalla validazione con ``E-ECO-40``, che nomina il problema.

    **Razionale scientifico e sistemistico**: un refuso in un valore
    toglierebbe in silenzio dei campioni dal loro gruppo.
    """
    oggetto = _oggetto(tmp_path)
    elenco = {"from": "luogo", "groups": {"nodi": ["NOD1A1", "NOD1A2", "NOD1A3", "NOD2B1"],
                                          "laboratorio": ["LAB1C1", "LAB1C2"]}}
    config = schema.valida(configurazione(design={"variable": "zona",
                                                  "groupings": {"zona": elenco}}))
    esito = valida_sull_oggetto(oggetto, config)
    assert esito.riepilogo["gruppi_derivati"][0]["gruppi"] == {"laboratorio": 2, "nodi": 4}
    assert esito.riepilogo["gruppi_derivati"][0]["campioni_senza_gruppo"] == 4

    for regola, parola in (
        ({"from": "luogo", "groups": {"nodi": ["NOD9"]}}, "NOD9"),
        ({"from": "colonna_che_non_c_e", "pattern": "^(.)"}, "colonna_che_non_c_e"),
    ):
        with pytest.raises(ErroreEco) as rifiuto:
            valida_sull_oggetto(oggetto, schema.valida(configurazione(
                design={"groupings": {"zona": regola}})))
        assert rifiuto.value.codici == ["E-ECO-40"] and parola in str(rifiuto.value)
    with pytest.raises(ErroreEco) as rifiuto:
        valida_sull_oggetto(oggetto, schema.valida(configurazione(
            design={"groupings": {"luogo": {"from": "gruppo", "pattern": "^(.)"}}})))
    assert rifiuto.value.codici == ["E-ECO-40"] and "colonna dell'oggetto" in str(rifiuto.value)


def test_i_colori_descrittivi_non_hanno_test_e_non_cambiano_le_altre_uscite(eco_r, tmp_path):
    """
    **Obiettivo**: ``ord.color_by`` aggiunge per ogni ordinazione un grafico
    ``<metodo>_<distanza>_colore_NN_<colonna>.png`` e nessuna riga nelle
    tabelle dei test; tutti gli altri file restano identici byte per byte a
    quelli di un'analisi senza.

    **Razionale scientifico e sistemistico**: un raggruppamento scelto dopo
    aver visto i dati non va testato sugli stessi dati: e' una figura, e non
    deve toccare i test.
    """
    oggetto = _oggetto(tmp_path)
    ordinazione = {"methods": ["pcoa"], "distances": ["bray"], "pcoa_correction": "none"}
    senza = tmp_path / "senza"
    esegui(oggetto, schema.valida(configurazione(ord=ordinazione)), senza)
    con = tmp_path / "con"
    esegui(oggetto, schema.valida(configurazione(
        design={"groupings": {"modulo": {"from": "luogo", "pattern": "^([A-Z]+[0-9])"}}},
        ord={**ordinazione, "color_by": ["modulo", "luogo"]})), con)
    prima, dopo = impronte(senza), impronte(con)
    nuovi = set(dopo) - set(prima)
    assert nuovi == {"ordinazione/pcoa_bray_colore_01_modulo.png",
                     "ordinazione/pcoa_bray_colore_02_luogo.png"}
    diversi = {n for n in prima if prima[n] != dopo[n]}
    assert diversi == {"configurazione.json", "manifest.json", "riepilogo.json", "campioni.tsv"}
    _, test = leggi_tsv(con / "test" / "permanova.tsv")
    assert {r["variabile"] for r in test} == {"gruppo"}


# --------------------------------------------------------------------------- #
# 3. La palette                                                                #
# --------------------------------------------------------------------------- #


def test_la_palette_si_estende_per_interpolazione_e_un_nome_sconosciuto_e_respinto(
    eco_r, tmp_path
):
    """
    **Obiettivo**: in R, con ``"Dark 3"`` i colori sono quelli di sempre
    (``hcl.colors``); con ``"Okabe-Ito"`` i primi otto sono quelli di
    ``palette.colors`` senza il grigio (riservato ad «altri»), e con dodici gruppi la palette si estende per
    interpolazione; un elenco esadecimale vale cosi' com'e'; un nome di
    RColorBrewer da' i colori di ColorBrewer. Un nome che nessun pacchetto
    conosce e' respinto dalla validazione con ``E-ECO-41``; con un'altra
    palette i grafici cambiano e le tabelle no.

    **Razionale scientifico e sistemistico**: i colori servono a distinguere i
    gruppi, anche a chi non distingue rosso e verde; non devono cambiare alcun
    numero.
    """
    script = tmp_path / "palette.R"
    script.write_text("\n".join([
        f'source("{RADICE / "src" / "amplicon16s_eco" / "R" / "lib" / "uscite.R"}")',
        'stopifnot(identical(colori_gruppi(4), grDevices::hcl.colors(4, "Dark 3")))',
        'imposta_palette("Okabe-Ito")',
        'oi <- unname(grDevices::palette.colors(NULL, "Okabe-Ito"))',
        'oi <- oi[oi != "#999999"]',
        'stopifnot(length(oi) == 8, !"#999999" %in% colori_gruppi(20))',
        'stopifnot(identical(colori_gruppi(5), oi[1:5]))',
        'dodici <- colori_gruppi(12)',
        'stopifnot(length(dodici) == 12, length(unique(dodici)) == 12,',
        '          identical(dodici, grDevices::colorRampPalette(oi)(12)))',
        'imposta_palette(c("#112233", "#445566"))',
        'stopifnot(identical(colori_gruppi(2), c("#112233", "#445566")),',
        '          length(colori_gruppi(3)) == 3)',
        'imposta_palette("Set2")',
        'stopifnot(identical(colori_gruppi(4), RColorBrewer::brewer.pal(4, "Set2")))',
        'stopifnot(is.na(origine_palette("Inesistente")))',
    ]), encoding="utf-8")
    esegui_rscript(script)

    oggetto = _oggetto(tmp_path)
    with pytest.raises(ErroreEco) as rifiuto:
        valida_sull_oggetto(oggetto, schema.valida(configurazione(
            plot={"palette": "Inesistente"})))
    assert rifiuto.value.codici == ["E-ECO-41"]
    prima = tmp_path / "prima"
    esegui(oggetto, schema.valida(configurazione()), prima)
    dopo = tmp_path / "dopo"
    esegui(oggetto, schema.valida(configurazione(plot={"palette": "Okabe-Ito"})), dopo)
    a, b = impronte(prima), impronte(dopo)
    cambiati = {n for n in a if a[n] != b[n]}
    assert cambiati == {"configurazione.json", "manifest.json", "composizione/gruppi_01_gruppo.png",
                        "alfa/gruppi_01_gruppo.png"}


def test_una_palette_estesa_per_interpolazione_e_dichiarata_con_e_eco_46(eco_r, tmp_path):
    """
    **Obiettivo**: con Okabe-Ito (otto colori, senza il grigio) e
    ``ord.color_by`` su una colonna di dieci valori, la palette si estende per
    interpolazione: ``validate`` e ``run`` danno l'avviso ``E-ECO-46``, il
    riepilogo dice quale palette, di quanti colori e per quali grafici, e
    l'intestazione della tabella delle coordinate lo scrive per quel grafico.
    Con una palette di due colori esadecimali e tre taxa lo stesso vale per le
    tabelle della composizione (barre per gruppo e pannelli) e, con tre gruppi
    nel confronto, per la tabella dei test dell'alfa diversita'. Con la palette
    predefinita (HCL, tanti colori quanti ne servono) nessun avviso e nessuna
    nota.

    **Razionale scientifico e sistemistico**: i colori di partenza di una
    palette adatta ai daltonici sono scelti per essere distinti; quelli
    interpolati no, e chi legge il grafico deve saperlo.
    """
    from amplicon16s_eco.catalogo import CATALOGO

    assert CATALOGO["E-ECO-46"].specie is Specie.AVVISO
    nota = ("palette estesa per interpolazione: i colori aggiunti non sono garantiti per chi "
            "non distingue i colori")
    oggetto = _oggetto(tmp_path)
    pcoa = {"methods": ["pcoa"], "distances": ["bray"], "pcoa_correction": "none",
            "color_by": ["luogo"]}
    config = schema.valida(configurazione(ord=pcoa, plot={"palette": "Okabe-Ito"}))
    validato = valida_sull_oggetto(oggetto, config)
    avvisi = {a.codice: a.dettaglio for a in validato.avvisi}
    assert "E-ECO-46" in avvisi and "ord.color_by 'luogo' (10 colori)" in avvisi["E-ECO-46"]
    uscita = tmp_path / "okabe"
    esito = esegui(oggetto, config, uscita)
    assert "E-ECO-46" in {a.codice for a in esito.avvisi}
    assert esito.riepilogo["palette"] == {
        "nome": "'Okabe-Ito'", "colori": 8, "nota": nota,
        "estesa_per_interpolazione": [{"grafici": "ord.color_by 'luogo'", "colori": 10}]}
    commenti, _ = leggi_tsv(uscita / "ordinazione" / "pcoa_bray_coordinate.tsv")
    assert any(c.startswith("grafico pcoa_bray_colore_01_luogo.png: " + nota) for c in commenti)
    assert not any("gruppi_01_gruppo" in c for c in commenti)

    due = schema.valida(configurazione(
        design={"variable": "zona",
                "groupings": {"zona": {"from": "luogo", "pattern": "^([A-Z]+[0-9])"}}},
        comp={"facet_wrap": "zona"}, plot={"palette": ["#1B9E77", "#D95F02"]}))
    uscita = tmp_path / "due"
    esito = esegui(oggetto, due, uscita)
    assert [g["grafici"] for g in esito.riepilogo["palette"]["estesa_per_interpolazione"]] == [
        "variabile 'zona'", "taxa della composizione (comp.top_n)"]
    for tabella in ("gruppi_01_zona.tsv", "pannelli_zona.tsv"):
        commenti, _ = leggi_tsv(uscita / "composizione" / tabella)
        assert any(c.startswith(nota) and "usa 3 colori" in c for c in commenti), tabella
    # I tre gruppi NOD1, NOD2 e LAB1 entrano tutti nel confronto.
    commenti, _ = leggi_tsv(uscita / "test" / "alfa_diversita.tsv")
    assert any(c.startswith("grafico alfa/gruppi_01_zona.png: " + nota) for c in commenti)

    predefinita = tmp_path / "predefinita"
    esito = esegui(oggetto, schema.valida(configurazione(ord=pcoa)), predefinita)
    assert "E-ECO-46" not in {a.codice for a in esito.avvisi}
    assert "palette" not in esito.riepilogo
    for tabella in predefinita.rglob("*.tsv"):
        assert "E-ECO-46" not in tabella.read_text(encoding="utf-8"), tabella


# --------------------------------------------------------------------------- #
# 4. Le proposte di eco init                                                   #
# --------------------------------------------------------------------------- #


def _colonna(nome, valori: dict[str, int], numerica=False):
    return {"nome": nome, "distinti": len(valori), "con_valore": sum(valori.values()),
            "gruppi_ripetuti": sum(n >= 2 for n in valori.values()), "numerica": numerica,
            "valori": [{"valore": v, "campioni": n} for v, n in sorted(valori.items())]}


def test_eco_init_elenca_le_colonne_e_propone_in_commento_il_gruppo_per_prefisso():
    """
    **Obiettivo**: sopra ``design`` il file elenca ogni colonna
    selezionabile con il numero di gruppi, quanti hanno almeno 2, 5 e 10
    campioni e, se sono pochi, i valori; le altre con il motivo (tecnica, un
    solo valore, tutti diversi, identificativo). Per una colonna con molti
    valori il cui prefisso (lettere seguite da cifre) li divide in pochi gruppi
    ripetuti, ``design.groupings`` porta la proposta in commento con i gruppi e
    i campioni che produrrebbe; ``--accetta-proposte`` non la attiva; tolto il
    commento, e' una configurazione valida. Una colonna senza un prefisso
    comune non riceve proposte.

    **Razionale scientifico e sistemistico**: i gruppi li sceglie chi conduce
    lo studio; init mostra che cosa i dati permettono, e la regola la decide
    l'utente.
    """
    from init_aiuti import completa, righe_dei_parametri

    from amplicon16s_eco import init

    posizioni = {f"{m}{s}{i}": 3 for m in ("NOD1", "NOD2", "LAB1") for s in "AB"
                 for i in range(1, 4)}
    posizioni.update({"Aria": 4, "Tubo": 2})
    descrizione = {
        "campioni": sum(posizioni.values()), "varianti": 10,
        "ranghi": [{"rango": "Regno", "quota_assegnata": 1.0}],
        "albero": {"presente": False, "radicato": None, "con_lunghezze": None},
        "colonne": [
            {"nome": "accession", "distinti": 60, "con_valore": 60, "gruppi_ripetuti": 0,
             "numerica": False, "valori": []},
            _colonna("posizione_prelievo", posizioni),
            _colonna("stagione", {"estate": 30, "inverno": 30}),
            _colonna("piastra", {"p1": 30, "p2": 30}),
            _colonna("costante", {"x": 60}),
            _colonna("codice", {f"Z{i}": 2 for i in range(30)}),
        ],
    }
    voci, fatti = init.prepara(descrizione)
    colonne = init.colonne_per_gruppi(descrizione, None)
    testo = init.testo_yaml(voci, fatti, "ps.rds", False, colonne)
    intestazione = " ".join(testo.split("design:")[0].replace("#", " ").split())
    assert ("- posizione_prelievo: 20 gruppi; con almeno 2, 5, 10 campioni: 20, 0, 0"
            in intestazione)
    assert "- stagione: 2 gruppi" in intestazione and "'estate' 30, 'inverno' 30" in intestazione
    assert "- piastra: variabile tecnica" in intestazione
    assert "- costante: un solo valore" in intestazione
    assert "- accession: identifica il campione" in intestazione
    assert "Colonne derivate dalla pipeline: non note" in intestazione
    righe = righe_dei_parametri(testo, schema.ConfigEco)
    attiva, valore = righe["design.groupings"]
    assert not attiva and "posizione_prelievo_prefisso" in valore
    assert "codice_prefisso" not in valore
    spiegazioni = " ".join(testo.replace("#", " ").split())
    assert "LAB1 18, NOD1 18, NOD2 18; senza gruppo: 6 campioni" in spiegazioni
    assert init.PREFISSO in spiegazioni

    accettato = init.testo_yaml(init.accetta(voci), fatti, "ps.rds", True, colonne)
    assert righe_dei_parametri(accettato, schema.ConfigEco)["design.groupings"][0] is False
    completo = completa(testo, {"run.seed": 1, "comp.top_n": 3, "beta.distances": ["bray"],
                                "ord.methods": [], "ord.distances": [],
                                "stat.min_group_size": 2, "stat.significance_level": 0.05,
                                "design.variable": "posizione_prelievo_prefisso",
                                "comp.rank": "Regno"},
                        accetta={"design.groupings", "alpha.method"}, modello=schema.ConfigEco)
    import yaml
    config = schema.valida(yaml.safe_load(completo))
    regola = config.design.groupings["posizione_prelievo_prefisso"]
    assert (regola.da, regola.pattern) == ("posizione_prelievo", init.PREFISSO)


# --------------------------------------------------------------------------- #
# 5. Composizione per gruppi scelti                                            #
# --------------------------------------------------------------------------- #


def test_la_composizione_si_disegna_per_le_variabili_scelte_senza_cambiare_i_test(
    eco_r, tmp_path
):
    """
    **Obiettivo**: con ``comp.group_by`` le barre della composizione sono per le
    variabili elencate (anche un gruppo derivato), con il numero nell'ordine
    dell'elenco, e quelle delle variabili del disegno non si disegnano; i test
    e le altre uscite restano identici byte per byte.

    **Razionale scientifico e sistemistico**: la composizione per modulo con
    una variabile dei test diversa e' una figura, non un'analisi: non deve
    cambiare i test.
    """
    oggetto = _oggetto(tmp_path)
    modulo = {"modulo": {"from": "luogo", "pattern": "^([A-Z]+[0-9])"}}
    senza = tmp_path / "senza"
    esegui(oggetto, schema.valida(configurazione(design={"groupings": modulo})), senza)
    con = tmp_path / "con"
    esegui(oggetto, schema.valida(configurazione(design={"groupings": modulo},
                                                 comp={"group_by": ["modulo", "luogo"]})), con)
    prima, dopo = impronte(senza), impronte(con)
    assert set(prima) - set(dopo) == {"composizione/gruppi_01_gruppo.tsv",
                                      "composizione/gruppi_01_gruppo.png"}
    assert set(dopo) - set(prima) == {f"composizione/gruppi_{n}.{e}" for n in (
        "01_modulo", "02_luogo") for e in ("tsv", "png")}
    for nome in ("test/alfa_diversita.tsv", "test/permanova.tsv", "campioni.tsv",
                 "composizione/abbondanze_relative.tsv"):
        assert prima[nome] == dopo[nome], nome
    commenti, righe = leggi_tsv(con / "composizione" / "gruppi_01_modulo.tsv")
    assert list(righe[0])[1:] == ["LAB1", "NOD1", "NOD2"]
    assert any("LAB1 = 3; NOD1 = 3; NOD2 = 2" in c for c in commenti)
    assert _codici(configurazione(comp={"group_by": ["luogo", "luogo"]})) == ["E-ECO-01"]


def test_la_media_per_gruppo_coincide_con_tax_glom_tranne_i_non_assegnati(eco_r, tmp_path):
    """
    **Obiettivo**: con ``comp.rank`` Famiglia le proporzioni medie per gruppo
    dei taxa assegnati coincidono con quelle di ``transform_sample_counts``
    (proporzioni per campione), ``tax_glom(ps, taxrank = "Famiglia")`` e media
    per gruppo, calcolate in R in modo indipendente; le varianti senza
    famiglia, che ``tax_glom`` toglie (``NArm = TRUE``), restano nella voce
    «non assegnato» e le colonne sommano a 1.

    **Razionale scientifico e sistemistico**: e' il calcolo che i postdoc
    farebbero a mano: deve dare gli stessi numeri, e la sola differenza deve
    essere dichiarata.
    """
    import json

    conteggi = [[20, 5, 30, 10], [10, 15, 0, 25], [5, 5, 10, 5], [15, 25, 10, 10], [0, 10, 5, 0]]
    tassonomia = {"Regno": ["B"] * 5, "Famiglia": ["F1", "F1", "F2", None, "F3"],
                  "Genere": ["G1", "G2", "G3", None, None]}
    cartella = tmp_path / "dati"
    cartella.mkdir()
    oggetto = costruisci(cartella, conteggi, {"gruppo": ["a", "a", "b", "b"]}, tassonomia)
    uscita = tmp_path / "uscita"
    esegui(oggetto, schema.valida(configurazione(comp={"rank": "Famiglia", "top_n": 10})), uscita)
    _, righe = leggi_tsv(uscita / "composizione" / "gruppi_01_gruppo.tsv")
    nostre = {r["taxon"]: (float(r["a"]), float(r["b"])) for r in righe}

    script = tmp_path / "glom.R"
    atteso = tmp_path / "atteso.json"
    script.write_text(f"""
ps <- readRDS("{oggetto}")
rel <- phyloseq::transform_sample_counts(ps, function(x) x / sum(x))
fam <- phyloseq::tax_glom(rel, taxrank = "Famiglia")
m <- methods::as(phyloseq::otu_table(fam), "matrix")
nomi <- as.character(phyloseq::tax_table(fam)[rownames(m), "Famiglia"])
g <- phyloseq::sample_data(ps)$gruppo
media <- sapply(c("a", "b"), function(x) rowMeans(m[, g == x, drop = FALSE]))
rownames(media) <- nomi
jsonlite::write_json(as.list(as.data.frame(t(media))), "{atteso}", digits = NA)
""", encoding="utf-8")
    esegui_rscript(script)
    di_tax_glom = json.loads(atteso.read_text(encoding="utf-8"))
    assert set(di_tax_glom) == {"F1", "F2", "F3"}
    for famiglia, valori in di_tax_glom.items():
        assert nostre[famiglia] == pytest.approx(tuple(valori), abs=1e-9), famiglia
    # La variante senza famiglia: 10/50, 25/60 nel gruppo a; 5/55, 10/50 nel b.
    assert nostre["non assegnato"] == pytest.approx(
        ((15 / 50 + 25 / 60) / 2, (10 / 55 + 10 / 50) / 2), abs=1e-9)
    assert nostre["altri"] == pytest.approx((0.0, 0.0), abs=1e-9)


# --------------------------------------------------------------------------- #
# 6. Grafici a pannelli                                                        #
# --------------------------------------------------------------------------- #


def test_i_grafici_a_pannelli_disegnano_i_valori_della_loro_tabella(eco_r, tmp_path):
    """
    **Obiettivo**: con ``comp.facet_wrap`` il grafico ha un pannello per gruppo
    e la tabella accanto riporta, per ogni campione con un valore, le
    abbondanze dei primi taxa e di «altri» di ``abbondanze_relative.tsv``; con
    ``comp.facet_grid`` la tabella riporta per ogni cella (riga, colonna,
    gruppo della variabile) la media calcolata qui a mano e il numero di
    campioni. Due esecuzioni danno gli stessi byte; una colonna che non esiste
    e' respinta con ``E-ECO-05``, righe e colonne uguali con ``E-ECO-01``.

    **Razionale scientifico e sistemistico**: un grafico va letto con i suoi
    numeri accanto; i pannelli sono figure della stessa composizione, non un
    altro calcolo.
    """
    oggetto = _oggetto(tmp_path)
    pannelli = {"top_n": 2, "facet_wrap": "modulo", "facet_grid": {"rows": "modulo",
                                                                   "cols": "gruppo"}}
    modulo = {"modulo": {"from": "luogo", "pattern": "^([A-Z]+[0-9])"}}
    config = schema.valida(configurazione(design={"groupings": modulo}, comp=pannelli))
    esegui(oggetto, config, tmp_path / "a")
    esegui(oggetto, config, tmp_path / "b")
    assert impronte(tmp_path / "a") == impronte(tmp_path / "b")
    cartella = tmp_path / "a" / "composizione"
    for nome in ("pannelli_modulo", "griglia_modulo_gruppo"):
        assert (cartella / f"{nome}.png").stat().st_size > 1000, nome
    _, relative = leggi_tsv(cartella / "abbondanze_relative.tsv")
    primi = [r["taxon"] for r in relative[:2]]
    per_campione = {(r["taxon"], c): float(v) for r in relative for c, v in r.items()
                    if c.startswith("c")}
    commenti, wrap = leggi_tsv(cartella / "pannelli_modulo.tsv")
    assert {r["campione"] for r in wrap} == {f"c{i:02d}" for i in (1, 2, 3, 4, 5, 6, 7, 9)}
    for riga in wrap:
        if riga["taxon"] == "altri":
            atteso = 1 - sum(per_campione[(t, riga["campione"])] for t in primi)
        else:
            atteso = per_campione[(riga["taxon"], riga["campione"])]
        assert float(riga["abbondanza_relativa"]) == pytest.approx(atteso, abs=1e-9)
    assert any("campioni senza valore in 'modulo', non disegnati: 2" in c for c in commenti)
    # Nel pannello, i campioni in ordine decrescente del taxon piu' abbondante
    # nel complesso (il primo di abbondanze_relative.tsv), a parita' per nome.
    assert any(f"per abbondanza decrescente di '{primi[0]}'" in c for c in commenti)
    ordine = list(dict.fromkeys((r["gruppo"], r["campione"]) for r in wrap))
    gruppi_in_ordine = list(dict.fromkeys(g for g, _ in ordine))
    assert gruppi_in_ordine == sorted(gruppi_in_ordine)
    for gruppo in gruppi_in_ordine:
        nel_gruppo = [c for g, c in ordine if g == gruppo]
        assert nel_gruppo == sorted(nel_gruppo, key=lambda c: (-per_campione[(primi[0], c)], c))
    _, griglia = leggi_tsv(cartella / "griglia_modulo_gruppo.tsv")
    # Cella NOD1 / a / a: i campioni NOD1 del gruppo a sono c01 e c09.
    cella = [r for r in griglia if (r["riga"], r["colonna"], r["gruppo"]) == ("NOD1", "a", "a")]
    assert {r["campioni"] for r in cella} == {"2"}
    for riga in cella:
        if riga["taxon"] != "altri":
            atteso = (per_campione[(riga["taxon"], "c01")] + per_campione[(riga["taxon"], "c09")]) / 2
            assert float(riga["abbondanza_relativa_media"]) == pytest.approx(atteso, abs=1e-9)
    with pytest.raises(ErroreEco) as rifiuto:
        valida_sull_oggetto(oggetto, schema.valida(configurazione(comp={"facet_wrap": "stanza"})))
    assert rifiuto.value.codici == ["E-ECO-05"]
    assert _codici(configurazione(comp={"facet_grid": {"rows": "luogo", "cols": "luogo"}})) == [
        "E-ECO-01"]


def _script_dei_pannelli(cartella, corpo: str) -> str:
    """Uno script R che carica le librerie dei grafici della composizione e
    definisce ``pannelli(m, gruppi, sovrapposizione)``: il grafico per campione
    di una matrice taxa x campioni, con un gruppo per campione, e le colonne di
    pixel bianchi fra due pixel colorati lungo piu' di 600 righe (una riga
    bianca verticale dentro un pannello)."""
    lib = RADICE / "src" / "amplicon16s_eco" / "R" / "lib"
    return "\n".join([
        f'source("{RADICE / "R" / "lib" / "io_json.R"}")',
        f'source("{lib / "uscite.R"}")',
        f'source("{lib / "composizione.R"}")',
        "pannelli <- function(m, gruppi, sovrapposizione = SOVRAPPOSIZIONE_BARRE) {",
        "  vecchia <- SOVRAPPOSIZIONE_BARRE",
        "  SOVRAPPOSIZIONE_BARRE <<- sovrapposizione",
        "  on.exit(SOVRAPPOSIZIONE_BARRE <<- vecchia)",
        "  dimensioni <- table(gruppi)",
        "  v <- list(nome = 'g', valori = stats::setNames(gruppi, colnames(m)),",
        "            dimensioni = stats::setNames(as.integer(dimensioni), names(dimensioni)),",
        "            mancanti = character())",
        f'  d <- tempfile(tmpdir = "{cartella}"); dir.create(d)',
        "  esito <- .pannelli_per_campione(m, v, 'R', d, 'x')",
        "  img <- png::readPNG(file.path(d, 'pannelli_g.png'))[, , 1:3]",
        "  b <- img[, , 1] == 1 & img[, , 2] == 1 & img[, , 3] == 1",
        "  W <- ncol(b)",
        "  f <- b[, 2:(W - 1)] & !b[, 1:(W - 2)] & !b[, 3:W]",
        "  esito$fessure <- sum(colSums(f) > 600)",
        "  esito$cartella <- d",
        "  esito",
        "}",
        "casuale <- function(n) {",
        "  set.seed(1)",
        "  m <- matrix(stats::runif(3 * n, 0, 0.3), nrow = 3)",
        "  m <- rbind(m, 1 - colSums(m))",
        "  dimnames(m) <- list(c(sprintf('taxon_%02d', 1:3), 'altri'), sprintf('c%03d', seq_len(n)))",
        "  m",
        "}",
        corpo,
    ])


def test_i_pannelli_non_lasciano_righe_bianche_fra_barre_vicine(eco_r, tmp_path):
    """
    **Obiettivo**: nel grafico per campione (``comp.facet_wrap``) fra due barre
    vicine non resta una riga bianca verticale. Con 78 campioni in un pannello
    e barre larghe esattamente 1 la riga compare (cosi' si sa che la prova la
    vede); con la sovrapposizione di un millesimo dello strumento no, ne'
    con 78 ne' con 118, 132 o 156 campioni, e il disegno non da' avvisi.

    **Razionale scientifico e sistemistico**: una riga bianca sembra un
    campione senza letture; nella dimostrazione su OSD-734 ne compariva una
    nel pannello NOD2, dove nessun campione aveva valori mancanti: era un
    difetto del disegno, i bordi di due barre calcolati separatamente ai due
    lati del centro di un pixel.
    """
    script = tmp_path / "fessure.R"
    script.write_text(_script_dei_pannelli(tmp_path, "\n".join([
        "stopifnot(pannelli(casuale(78), rep('G1', 78), 0)$fessure > 0)",
        "for (n in c(78, 118, 132, 156)) {",
        "  avvisi <- NULL",
        "  esito <- withCallingHandlers(pannelli(casuale(n), rep('G1', n)),",
        "    warning = function(w) { avvisi <<- c(avvisi, conditionMessage(w)); invokeRestart('muffleWarning') })",
        "  if (esito$fessure != 0 || length(avvisi)) stop(n, ': ', esito$fessure, ' ', paste(avvisi, collapse = '; '))",
        "}",
    ])), encoding="utf-8")
    esegui_rscript(script)


def test_un_campione_senza_abbondanze_e_dichiarato_nei_pannelli(eco_r, tmp_path):
    """
    **Obiettivo**: un campione le cui abbondanze non sono calcolabili (nessuna
    lettura al rango: 0/0) non e' disegnato e non lascia uno spazio vuoto; la
    tabella del grafico lo nomina, il numero del suo gruppo non lo conta e la
    funzione lo restituisce per il riepilogo.

    **Razionale scientifico e sistemistico**: uno spazio vuoto in un grafico
    non dice perche' c'e'; una nota si'.
    """
    uscita = tmp_path / "esito.json"
    script = tmp_path / "senza.R"
    script.write_text(_script_dei_pannelli(tmp_path, "\n".join([
        "m <- casuale(10)",
        "m[, 'c004'] <- NaN",
        "esito <- pannelli(m, rep(c('A', 'B'), each = 5))",
        "stopifnot(esito$campioni_disegnati == 9, esito$campioni_senza_abbondanze == 1)",
        "righe <- readLines(file.path(esito$cartella, 'pannelli_g.tsv'))",
        "stopifnot(any(grepl('non disegnati: 1 (c004)', righe, fixed = TRUE)),",
        "          any(grepl('A = 4; B = 5', righe, fixed = TRUE)),",
        "          !any(grepl('c004\\t', righe, fixed = TRUE)))",
        f'jsonlite::write_json(esito[c("campioni_disegnati", "campioni_senza_abbondanze")], "{uscita}", auto_unbox = TRUE)',
    ])), encoding="utf-8")
    esegui_rscript(script)
    import json

    assert json.loads(uscita.read_text(encoding="utf-8")) == {
        "campioni_disegnati": 9, "campioni_senza_abbondanze": 1}
