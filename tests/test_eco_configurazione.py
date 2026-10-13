"""Validazione delle analisi ecologiche: configurazione, oggetto, riga di comando.

Verifica che un'analisi mal dichiarata sia respinta prima di ogni calcolo, con
un codice del catalogo del pacchetto che dice la causa e che cosa correggere:
il ruolo che la validazione della configurazione ha per la pipeline. Verifica
inoltre che il pacchetto non porti con se' valori dei dataset distribuiti e che
il comando di esecuzione dal clone non cambi comportamento per la pipeline.
"""

from __future__ import annotations

import re
import subprocess
import sys
from pathlib import Path

import pytest
import yaml

from amplicon16s_eco import config as schema
from amplicon16s_eco.catalogo import CATALOGO, ErroreEco, Specie, descrivi
from amplicon16s_eco.cli import USCITA_RIFIUTO, USCITA_SUCCESSO, main
from amplicon16s_eco.esecuzione import esegui, valida_sull_oggetto

from eco_aiuti import RADICE, configurazione, costruisci, eco_r, impronte  # noqa: F401

PACCHETTO = RADICE / "src" / "amplicon16s_eco"
MODELLO = RADICE / "config" / "eco.example.yaml"

#: Sei campioni in due gruppi da tre, quattro varianti.
CONTEGGI = [
    [10, 12, 9, 2, 1, 3],
    [5, 4, 6, 8, 9, 7],
    [3, 2, 4, 6, 5, 6],
    [0, 1, 0, 9, 8, 7],
]
METADATI = {
    "gruppo": ["a", "a", "a", "b", "b", "b"],
    "lotto": ["x", "y", "x", "y", "x", "y"],
}


def _codici(dati) -> list[str]:
    with pytest.raises(ErroreEco) as rifiuto:
        schema.valida(dati)
    return rifiuto.value.codici


# --------------------------------------------------------------------------- #
# 1. Il file di configurazione                                                 #
# --------------------------------------------------------------------------- #


def test_una_configurazione_completa_e_accettata_e_il_digest_dipende_dai_valori():
    """
    **Obiettivo**: una configurazione con tutti gli obbligatori e' accettata; il
    digest non cambia con l'ordine delle chiavi e cambia con un valore.

    **Razionale scientifico e sistemistico**: il digest registrato con i
    risultati deve identificare i valori con cui sono stati calcolati, non i
    byte del file.
    """
    dati = configurazione()
    rovesciata = {k: dict(reversed(list(v.items()))) for k, v in reversed(list(dati.items()))}
    assert schema.digest(schema.valida(dati)) == schema.digest(schema.valida(rovesciata))
    altra = configurazione(run={"seed": 8})
    assert schema.digest(schema.valida(dati)) != schema.digest(schema.valida(altra))


def test_senza_la_variabile_biologica_l_analisi_non_parte():
    """
    **Obiettivo**: senza ``design.variable`` la configurazione e' respinta con
    ``E-ECO-02`` e il messaggio nomina il parametro.

    **Razionale scientifico e sistemistico**: la variabile biologica dipende
    dallo studio e non puo' avere un valore predefinito.
    """
    with pytest.raises(ErroreEco) as rifiuto:
        schema.valida(configurazione(design={"variable": ...}))
    assert rifiuto.value.codici == ["E-ECO-02"]
    assert "design.variable" in str(rifiuto.value)


@pytest.mark.parametrize("chiave", schema.OBBLIGATORI)
def test_ogni_parametro_obbligatorio_mancante_e_respinto_e_nominato(chiave):
    """
    **Obiettivo**: tolto un parametro obbligatorio, il rifiuto e' ``E-ECO-02``
    e nomina proprio quel parametro.

    **Razionale scientifico e sistemistico**: nessun parametro che dipende
    dallo studio o che e' una scelta di analisi puo' valere per difetto.
    """
    gruppo, nome = chiave.split(".")
    with pytest.raises(ErroreEco) as rifiuto:
        schema.valida(configurazione(**{gruppo: {nome: ...}}))
    assert rifiuto.value.codici == ["E-ECO-02"]
    assert chiave in str(rifiuto.value)


def test_un_metodo_non_supportato_e_respinto_con_l_elenco_degli_ammessi():
    """
    **Obiettivo**: una distanza che il pacchetto non realizza e' respinta con
    ``E-ECO-03`` e il messaggio elenca quelle ammesse.

    **Razionale scientifico e sistemistico**: un metodo sconosciuto non deve
    essere ignorato ne' sostituito in silenzio con un altro.
    """
    with pytest.raises(ErroreEco) as rifiuto:
        schema.valida(configurazione(beta={"distances": ["bray", "manhattan"]}))
    assert rifiuto.value.codici == ["E-ECO-03"]
    testo = str(rifiuto.value)
    assert "manhattan" in testo and all(d in testo for d in schema.DISTANZE)


def test_le_chiavi_sconosciute_sono_respinte():
    """
    **Obiettivo**: una chiave non prevista, in un gruppo o alla radice, e'
    respinta con ``E-ECO-01`` che la nomina.

    **Razionale scientifico e sistemistico**: un refuso nel nome di un
    parametro lascerebbe valere il predefinito senza che nessuno se ne accorga.
    """
    assert _codici(configurazione(alpha={"rarefy_dept": 100})) == ["E-ECO-01"]
    assert _codici({**configurazione(), "extra": {}}) == ["E-ECO-01"]
    with pytest.raises(ErroreEco) as rifiuto:
        schema.valida(configurazione(alpha={"rarefy_dept": 100}))
    assert "alpha.rarefy_dept" in str(rifiuto.value)


def test_lo_pseudoconteggio_e_obbligatorio_con_aitchison_e_respinto_senza():
    """
    **Obiettivo**: ``beta.clr_pseudocount`` manca con Aitchison richiesta
    (``E-ECO-02``); e' dichiarato senza Aitchison (``E-ECO-01``).

    **Razionale scientifico e sistemistico**: lo pseudoconteggio cambia la
    distanza e va scelto; dichiarato senza il metodo, entrerebbe nel digest
    senza avere effetto.
    """
    assert _codici(configurazione(beta={"clr_pseudocount": ...})) == ["E-ECO-02"]
    assert _codici(configurazione(beta={"distances": ["bray"]})) == ["E-ECO-01"]
    schema.valida(configurazione(beta={"distances": ["bray"], "clr_pseudocount": ...}))


def test_i_parametri_dell_ordinazione_seguono_il_metodo_richiesto():
    """
    **Obiettivo**: ``ord.pcoa_correction`` e ``ord.nmds_trymax`` mancano con il
    metodo richiesto (``E-ECO-02``) e sono respinti senza (``E-ECO-01``); un
    metodo di ordinazione o una correzione non realizzati sono respinti con
    ``E-ECO-03`` e l'elenco degli ammessi.

    **Razionale scientifico e sistemistico**: la correzione degli autovalori
    negativi e il numero di avvii cambiano il risultato e non hanno un valore
    che valga per ogni studio.
    """
    pcoa = {"methods": ["pcoa"], "distances": ["bray"]}
    nmds = {"methods": ["nmds"], "distances": ["bray"]}
    assert _codici(configurazione(ord=pcoa)) == ["E-ECO-02"]
    assert _codici(configurazione(ord=nmds)) == ["E-ECO-02"]
    assert _codici(configurazione(ord={**pcoa, "pcoa_correction": "none", "nmds_trymax": 5})) == ["E-ECO-01"]
    assert _codici(configurazione(ord={"pcoa_correction": "none"})) == ["E-ECO-01"]
    schema.valida(configurazione(ord={**pcoa, "pcoa_correction": "lingoes"}))
    schema.valida(configurazione(ord={**nmds, "nmds_trymax": 5}))
    for modifica, parola in (
        ({**pcoa, "pcoa_correction": "radice"}, "cailliez"),
        ({"methods": ["tsne"], "distances": ["bray"]}, "nmds"),
    ):
        with pytest.raises(ErroreEco) as rifiuto:
            schema.valida(configurazione(ord=modifica))
        assert rifiuto.value.codici == ["E-ECO-03"] and parola in str(rifiuto.value)
    assert schema.valida(configurazione()).stat.permanova_permutations == 999


def test_la_radice_dell_albero_va_dichiarata_solo_con_unifrac():
    """
    **Obiettivo**: ``beta.unifrac_root`` manca (``E-ECO-02``) se e' richiesta
    una delle due UniFrac, pesata o non pesata; e' respinto (``E-ECO-01``) se
    non ne e' richiesta nessuna; un valore diverso da ``midpoint`` ed
    ``existing`` e' un metodo non realizzato (``E-ECO-03``).

    **Razionale scientifico e sistemistico**: UniFrac dipende dalla radice, e
    il modo di ottenerla cambia la distanza: non ha un valore che valga per
    ogni albero, e dichiarato senza UniFrac entrerebbe nel digest senza effetto.
    """
    for distanza in ("unifrac_weighted", "unifrac_unweighted"):
        beta = {"distances": ["bray", distanza], "clr_pseudocount": ...}
        assert _codici(configurazione(beta=beta)) == ["E-ECO-02"]
        for radice in schema.RADICI_UNIFRAC:
            config = schema.valida(configurazione(beta={**beta, "unifrac_root": radice}))
            assert config.beta.unifrac_root == radice
        with pytest.raises(ErroreEco) as rifiuto:
            schema.valida(configurazione(beta={**beta, "unifrac_root": "outgroup"}))
        assert rifiuto.value.codici == ["E-ECO-03"] and "midpoint" in str(rifiuto.value)
    assert _codici(configurazione(beta={"unifrac_root": "midpoint"})) == ["E-ECO-01"]
    assert schema.valida(configurazione()).beta.unifrac_root is None


@pytest.mark.parametrize("modifica", [
    {"run": {"seed": "sette"}},
    {"run": {"seed": 1.5}},
    {"comp": {"top_n": 0}},
    {"stat": {"min_group_size": 1}},
    {"alpha": {"rarefy_depth": 0}},
    {"beta": {"distances": []}},
    {"beta": {"distances": ["bray", "bray"], "clr_pseudocount": ...}},
    {"design": {"other_variables": ["gruppo"]}},
    {"design": {"technical_variables": ["gruppo"]}},
    {"design": {"subset": {"lotto": []}}},
    {"design": {"subset": {"lotto": [1]}}},
    {"stat": {"significance_level": 0}},
    {"stat": {"significance_level": 1}},
    {"stat": {"permanova_permutations": 0}},
    {"stat": {"permanova_strata": ""}},
    {"ord": {"methods": ["pcoa"], "distances": [], "pcoa_correction": "none"}},
    {"ord": {"methods": ["pcoa"], "distances": ["unifrac_weighted"], "pcoa_correction": "none"}},
    {"ord": {"methods": ["pcoa", "pcoa"], "distances": ["bray"], "pcoa_correction": "none"}},
    {"ord": {"methods": ["nmds"], "distances": ["bray"], "nmds_trymax": 0}},
])
def test_i_valori_fuori_dominio_sono_respinti(modifica):
    """
    **Obiettivo**: tipi sbagliati, valori fuori dai limiti, ripetizioni e
    sottoinsiemi vuoti sono respinti con ``E-ECO-01``.

    **Razionale scientifico e sistemistico**: un valore privo di senso deve
    fermare l'analisi prima del calcolo, non produrre un errore di R a meta'.
    """
    assert _codici(configurazione(**modifica)) == ["E-ECO-01"]


def test_tutti_i_problemi_sono_riportati_insieme():
    """
    **Obiettivo**: due problemi indipendenti danno due rifiuti nello stesso
    messaggio.

    **Razionale scientifico e sistemistico**: chi corregge una configurazione
    non deve scoprire gli errori uno per esecuzione.
    """
    codici = _codici(configurazione(design={"variable": ...}, beta={"distances": ["x"]}))
    assert sorted(codici) == ["E-ECO-02", "E-ECO-03"]


def test_il_modello_lascia_da_dichiarare_ogni_obbligatorio_e_i_marcatori_coincidono():
    """
    **Obiettivo**: il modello ``config/eco.example.yaml`` non assegna alcun
    parametro obbligatorio, e ogni parametro vi compare con il marcatore della
    sua specie secondo il codice.

    **Razionale scientifico e sistemistico**: il modello non deve suggerire
    valori, e la sua classificazione dei parametri non deve divergere dallo
    schema.
    """
    testo = MODELLO.read_text(encoding="utf-8")
    dati = yaml.safe_load(testo)
    with pytest.raises(ErroreEco) as rifiuto:
        schema.valida(dati)
    assert set(rifiuto.value.codici) == {"E-ECO-02"}

    marcatore: dict[str, str] = {}
    gruppo, corrente = None, None
    for riga in testo.splitlines():
        if m := re.match(r"^(\w+):\s*$", riga):
            gruppo = m.group(1)
        elif m := re.match(r"^  # \[(\w+)", riga):
            corrente = m.group(1)
        elif (m := re.match(r"^  (?:# )?(\w+):", riga)) and gruppo and corrente:
            marcatore[f"{gruppo}.{m.group(1)}"] = corrente
            corrente = None
    attesi = {
        **{c: "OBBLIGATORIO" for c in schema.OBBLIGATORI},
        **{c: "CONDIZIONALE" for c in schema.CONDIZIONALI},
        **{c: "FACOLTATIVO" for c in schema.FACOLTATIVI},
        **{c: "STANDARD" for c in schema.STANDARD_DEL_METODO},
    }
    assert marcatore == attesi
    campi = {
        f"{g}.{n}" for g, modello in schema.ConfigEco.model_fields.items()
        for n in modello.annotation.model_fields
    }
    assert campi == set(attesi)


def test_ogni_predefinito_e_l_assenza_o_un_valore_del_metodo_con_la_fonte():
    """
    **Obiettivo**: i soli campi con un predefinito sono i facoltativi (vuoti) e
    quelli elencati fra gli standard del metodo, con il valore dichiarato e una
    fonte.

    **Razionale scientifico e sistemistico**: un predefinito e' ammesso solo se
    e' quello del pacchetto che realizza il metodo.
    """
    for gruppo, modello in schema.ConfigEco.model_fields.items():
        for nome, campo in modello.annotation.model_fields.items():
            chiave = f"{gruppo}.{nome}"
            if campo.is_required():
                assert chiave in schema.OBBLIGATORI
            elif chiave in schema.STANDARD_DEL_METODO:
                valore, fonte = schema.STANDARD_DEL_METODO[chiave]
                assert campo.get_default(call_default_factory=True) == valore
                assert "::" in fonte
            elif chiave in schema.CONDIZIONALI:
                assert campo.get_default(call_default_factory=True) is None
            else:
                assert chiave in schema.FACOLTATIVI
                assert not campo.get_default(call_default_factory=True)


# --------------------------------------------------------------------------- #
# 2. Il catalogo                                                               #
# --------------------------------------------------------------------------- #


def test_i_codici_usati_e_quelli_del_catalogo_coincidono():
    """
    **Obiettivo**: ogni codice ``E-ECO-`` scritto nel sorgente Python o R del
    pacchetto e' nel catalogo, e ogni codice del catalogo e' usato.

    **Razionale scientifico e sistemistico**: un codice sollevato senza voce
    non avrebbe una causa ne' un'azione; una voce senza uso sarebbe un codice
    morto.
    """
    usati: set[str] = set()
    for percorso in [*PACCHETTO.rglob("*.py"), *PACCHETTO.rglob("*.R")]:
        if percorso.name != "catalogo.py":
            usati |= set(re.findall(r"E-ECO-\d{2}", percorso.read_text(encoding="utf-8")))
    assert usati == set(CATALOGO)


def test_ogni_voce_dice_la_causa_e_che_cosa_fare():
    """
    **Obiettivo**: ogni voce ha sintesi e azione non vuote, il codice nella
    forma attesa, e un avviso non puo' fermare l'analisi.

    **Razionale scientifico e sistemistico**: ogni rifiuto deve dire la causa e
    che cosa correggere; un avviso e un arresto non condividono un codice.
    """
    for codice, v in CATALOGO.items():
        assert re.fullmatch(r"E-ECO-\d{2}", codice)
        assert len(v.sintesi) > 10 and len(v.azione) > 10
        assert v.sintesi in descrivi(codice) and v.azione in descrivi(codice)
    avviso = next(c for c, v in CATALOGO.items() if v.specie is Specie.AVVISO)
    with pytest.raises(TypeError):
        ErroreEco([(avviso, "")])


# --------------------------------------------------------------------------- #
# 3. La validazione sull'oggetto                                               #
# --------------------------------------------------------------------------- #


def _rifiuto(oggetto, **modifica) -> ErroreEco:
    with pytest.raises(ErroreEco) as rifiuto:
        valida_sull_oggetto(oggetto, schema.valida(configurazione(**modifica)))
    return rifiuto.value


def test_una_colonna_inesistente_e_respinta_con_le_colonne_disponibili(eco_r, tmp_path):
    """
    **Obiettivo**: una colonna assente dalla tabella dei campioni (variabile
    biologica, tecnica o del sottoinsieme) e' respinta con ``E-ECO-05``, una
    volta per colonna, con l'elenco di quelle esistenti.

    **Razionale scientifico e sistemistico**: il nome di una colonna e' il
    punto in cui la configurazione incontra i metadati dello studio, e un
    refuso va detto prima del calcolo.
    """
    oggetto = costruisci(tmp_path, CONTEGGI, METADATI)
    rifiuto = _rifiuto(oggetto, design={"variable": "grupo", "technical_variables": ["loto"],
                                        "subset": {"assente": ["x"]}})
    assert rifiuto.codici == ["E-ECO-05"] * 3
    assert "gruppo, lotto" in str(rifiuto)


def test_un_gruppo_troppo_piccolo_e_dichiarato_prima_dell_esecuzione(eco_r, tmp_path):
    """
    **Obiettivo**: la sola validazione dichiara con ``E-ECO-11`` il gruppo
    sotto ``stat.min_group_size`` e con ``E-ECO-12`` che restano meno di due
    gruppi, senza scrivere nulla.

    **Razionale scientifico e sistemistico**: un gruppo con pochi campioni non
    sostiene un test; esce dai test con una dichiarazione, mentre le uscite
    descrittive restano valide e l'analisi non si ferma.
    """
    metadati = {**METADATI, "gruppo": ["a", "a", "a", "a", "a", "b"]}
    oggetto = costruisci(tmp_path, CONTEGGI, metadati)
    prima = impronte(tmp_path)
    esito = valida_sull_oggetto(oggetto, schema.valida(configurazione()))
    codici = [a.codice for a in esito.avvisi]
    assert codici[:2] == ["E-ECO-11", "E-ECO-12"]
    assert "b: 1" in esito.avvisi[0].dettaglio
    (variabile,) = esito.riepilogo["variabili"]
    assert variabile["gruppi_esclusi_dai_test"] == ["b"]
    assert variabile["gruppi_nei_test"] == ["a"]
    assert impronte(tmp_path) == prima


def test_unifrac_senza_albero_e_respinta(eco_r, tmp_path):
    """
    **Obiettivo**: una distanza UniFrac su un oggetto senza albero e' respinta
    con ``E-ECO-06`` qualunque radice sia dichiarata; su un albero senza radice
    e' respinta con ``existing`` e accettata con ``midpoint``; su un albero
    radicato e' accettata con entrambe; un albero senza lunghezze dei rami e'
    respinto.

    **Razionale scientifico e sistemistico**: senza albero la distanza non
    esiste; senza radice phyloseq ne sceglierebbe una a caso e il risultato non
    sarebbe riproducibile, quindi la radice o c'e' gia' o si ottiene in un modo
    dichiarato; UniFrac e il punto medio sono definiti sulle lunghezze dei rami.
    """
    def beta(radice):
        return {"distances": ["unifrac_unweighted"], "clr_pseudocount": ...,
                "unifrac_root": radice}

    senza = costruisci(tmp_path, CONTEGGI, METADATI)
    for radice in schema.RADICI_UNIFRAC:
        assert _rifiuto(senza, beta=beta(radice)).codici == ["E-ECO-06"]
    non_radicato = costruisci(tmp_path, CONTEGGI, METADATI, nome="non_radicato.rds",
                              albero="(v01:1,v02:1,(v03:1,v04:1):1);")
    rifiuto = _rifiuto(non_radicato, beta=beta("existing"))
    assert rifiuto.codici == ["E-ECO-06"] and "midpoint" in str(rifiuto)
    valida_sull_oggetto(non_radicato, schema.valida(configurazione(beta=beta("midpoint"))))
    radicato = costruisci(tmp_path, CONTEGGI, METADATI, nome="radicato.rds",
                          albero="((v01:1,v02:1):1,(v03:1,v04:1):1);")
    for radice in schema.RADICI_UNIFRAC:
        valida_sull_oggetto(radicato, schema.valida(configurazione(beta=beta(radice))))
    senza_rami = costruisci(tmp_path, CONTEGGI, METADATI, nome="senza_rami.rds",
                            albero="(v01,v02,(v03,v04));")
    for radice in schema.RADICI_UNIFRAC:
        rifiuto = _rifiuto(senza_rami, beta=beta(radice))
        assert rifiuto.codici == ["E-ECO-06"] and "lunghezze dei rami" in str(rifiuto)


def test_un_rango_assente_e_respinto_con_i_ranghi_dell_oggetto(eco_r, tmp_path):
    """
    **Obiettivo**: un ``comp.rank`` che l'oggetto non ha e' respinto con
    ``E-ECO-07`` e il messaggio elenca i ranghi disponibili.

    **Razionale scientifico e sistemistico**: i ranghi dipendono dal
    riferimento tassonomico usato dalla pipeline.
    """
    oggetto = costruisci(tmp_path, CONTEGGI, METADATI)
    rifiuto = _rifiuto(oggetto, comp={"rank": "Specie"})
    assert rifiuto.codici == ["E-ECO-07"]
    assert "Regno, Famiglia, Genere" in str(rifiuto)


def test_il_sottoinsieme_si_applica_prima_e_un_valore_assente_e_respinto(eco_r, tmp_path):
    """
    **Obiettivo**: il sottoinsieme restringe i campioni analizzati; un valore
    ammesso che non compare, o un sottoinsieme con meno di due campioni, e'
    respinto con ``E-ECO-08``.

    **Razionale scientifico e sistemistico**: un valore che non compare e' piu'
    probabilmente un refuso che un'intenzione, e analizzerebbe in silenzio un
    insieme diverso da quello voluto.
    """
    oggetto = costruisci(tmp_path, CONTEGGI, METADATI)
    esito = valida_sull_oggetto(
        oggetto, schema.valida(configurazione(design={"subset": {"lotto": ["x"]}})))
    assert esito.riepilogo["campioni"]["analizzati"] == 3
    assert esito.riepilogo["campioni"]["esclusi"]["fuori_dal_sottoinsieme"] == 3
    assert _rifiuto(oggetto, design={"subset": {"lotto": ["x", "z"]}}).codici == ["E-ECO-08"]
    metadati = {**METADATI, "lotto": ["x", "y", "y", "y", "y", "y"]}
    uno = costruisci(tmp_path, CONTEGGI, metadati, nome="uno.rds")
    assert _rifiuto(uno, design={"subset": {"lotto": ["x"]}}).codici == ["E-ECO-08"]


def test_i_valori_mancanti_escono_dalla_variabile_e_sono_registrati(eco_r, tmp_path):
    """
    **Obiettivo**: vuoto, ``NA`` e un marcatore di «non applicabile» nella
    variabile analizzata tolgono il campione dalle analisi di quella variabile,
    con ``E-ECO-13`` e il conteggio nel riepilogo.

    **Razionale scientifico e sistemistico**: un valore mancante non e' un
    gruppo: trattarlo come tale creerebbe un gruppo fittizio nei test.
    """
    metadati = {**METADATI, "gruppo": ["a", "a", None, "b", "b", "Not Applicable"]}
    oggetto = costruisci(tmp_path, CONTEGGI, metadati)
    esito = valida_sull_oggetto(oggetto, schema.valida(configurazione()))
    mancanti = [a for a in esito.avvisi if a.codice == "E-ECO-13"]
    assert len(mancanti) == 1 and "c03, c06" in mancanti[0].dettaglio
    (variabile,) = esito.riepilogo["variabili"]
    assert variabile["campioni_con_valore_mancante"] == 2
    assert variabile["gruppi"] == {"a": 2, "b": 2}


def test_un_oggetto_che_non_e_phyloseq_e_respinto(eco_r, tmp_path):
    """
    **Obiettivo**: un file assente, un file che non e' RDS e un RDS che non
    contiene un oggetto phyloseq sono respinti con ``E-ECO-04``.

    **Razionale scientifico e sistemistico**: il pacchetto lavora solo
    sull'oggetto consegnato dalla pipeline; qualunque altra cosa va respinta
    con la causa, non con un errore di R.
    """
    config = schema.valida(configurazione())
    testo = tmp_path / "testo.rds"
    testo.write_text("non sono un oggetto", encoding="utf-8")
    for percorso in (tmp_path / "assente.rds", testo):
        with pytest.raises(ErroreEco) as rifiuto:
            valida_sull_oggetto(percorso, config)
        assert rifiuto.value.codici == ["E-ECO-04"]


def test_la_cartella_di_uscita_deve_essere_vuota_e_fuori_da_quella_dell_oggetto(eco_r, tmp_path):
    """
    **Obiettivo**: una cartella di uscita dentro quella dell'oggetto, o non
    vuota, e' respinta con ``E-ECO-09`` senza scrivere nulla.

    **Razionale scientifico e sistemistico**: le uscite dell'analisi non devono
    mescolarsi con gli artefatti della pipeline ne' sovrascrivere un'analisi
    precedente.
    """
    oggetto = costruisci(tmp_path, CONTEGGI, METADATI)
    config = schema.valida(configurazione())
    piena = tmp_path / "piena"
    piena.mkdir()
    (piena / "resto.txt").write_text("x", encoding="utf-8")
    prima = impronte(tmp_path)
    for uscita in (oggetto.parent / "eco", oggetto.parent, piena):
        with pytest.raises(ErroreEco) as rifiuto:
            esegui(oggetto, config, uscita)
        assert rifiuto.value.codici == ["E-ECO-09"]
    assert impronte(tmp_path) == prima


# --------------------------------------------------------------------------- #
# 4. La riga di comando                                                        #
# --------------------------------------------------------------------------- #


def _scrivi_config(cartella: Path, dati) -> Path:
    percorso = cartella / "eco.yaml"
    percorso.write_text(yaml.safe_dump(dati), encoding="utf-8")
    return percorso


def test_la_riga_di_comando_respinge_con_il_codice_e_uscita_3(tmp_path, capsys):
    """
    **Obiettivo**: una configurazione senza la variabile biologica fa uscire
    ``validate`` e ``run`` con 3, il codice ``E-ECO-02`` sull'uscita di errore
    e nessuna cartella di uscita creata.

    **Razionale scientifico e sistemistico**: chi lancia l'analisi da uno
    script distingue un rifiuto da un difetto dal codice di uscita.
    """
    config = _scrivi_config(tmp_path, configurazione(design={"variable": ...}))
    uscita = tmp_path / "uscita"
    comune = ["--object", str(tmp_path / "x.rds"), "--config", str(config)]
    assert main(["validate", *comune]) == USCITA_RIFIUTO
    assert main(["run", *comune, "--out", str(uscita)]) == USCITA_RIFIUTO
    assert "[E-ECO-02]" in capsys.readouterr().err
    assert not uscita.exists()


def test_validate_indica_la_riga_del_file_di_ogni_parametro_mancante(tmp_path, capsys):
    """
    **Obiettivo**: sul modello di configurazione, che lascia in commento ogni
    parametro obbligatorio, ``validate`` stampa dopo il rifiuto, per ciascun
    parametro mancante, il numero della riga del file che lo porta in commento
    sotto il suo gruppo; per un parametro senza riga indica la riga del gruppo
    sotto cui aggiungerla; il messaggio del rifiuto e il documento JSON restano
    quelli di prima.

    **Razionale scientifico e sistemistico**: l'elenco dei parametri mancanti
    e' la lista delle cose da scrivere, e il numero della riga evita di
    scriverne uno sotto il gruppo sbagliato. La regola e' quella della
    pipeline, applicata ai parametri dello schema delle analisi: un parametro
    della pipeline con lo stesso nome non e' una riga di questo file.
    """
    import json

    from amplicon16s.cli import righe_nel_file
    from amplicon16s_eco.interfaccia import chiavi_dello_schema

    comune = ["--object", str(tmp_path / "x.rds"), "--config", str(MODELLO)]
    assert main(["validate", *comune]) == USCITA_RIFIUTO
    rifiuto = capsys.readouterr().err
    assert all(f"Parametro mancante: {chiave}" in rifiuto for chiave in schema.OBBLIGATORI)
    assert f"dove scriverli, in {MODELLO} (parametro: riga del file):" in rifiuto
    # Le indicazioni seguono il messaggio del rifiuto.
    assert rifiuto.index("dove scriverli") > rifiuto.rindex("[E-ECO-02]")
    righe = MODELLO.read_text(encoding="utf-8").splitlines()
    indicate = dict(re.findall(r"^    ([a-z_.0-9]+): (.+)$", rifiuto, flags=re.M))
    assert set(schema.OBBLIGATORI) <= set(indicate)
    for chiave in schema.OBBLIGATORI:
        gruppo, parametro = chiave.split(".")
        numeri = [int(n) for n in re.findall(r"\d+", indicate[chiave].split("(")[0])]
        assert numeri, chiave
        for numero in numeri:
            assert righe[numero - 1].startswith(f"  # {parametro}:"), (chiave, righe[numero - 1])
            sopra = [r for r in righe[:numero] if re.match(r"^[a-z_]+:\s*$", r)]
            assert sopra[-1] == f"{gruppo}:", chiave

    # Senza la riga lo dice, e una riga di spiegazione non e' un parametro.
    breve = tmp_path / "breve.yaml"
    breve.write_text("design:\n  ## variable: come si sceglie\n  technical_variables: []\n",
                     encoding="utf-8")
    assert main(["validate", "--object", str(tmp_path / "x.rds"), "--config", str(breve)]) == (
        USCITA_RIFIUTO)
    rifiuto = capsys.readouterr().err
    assert ("design.variable: nessuna riga nel file: aggiungila sotto «design:» (riga 1)"
            in rifiuto)

    # I parametri noti sono quelli dello schema indicato: per difetto la pipeline.
    testo = breve.read_text(encoding="utf-8")
    assert righe_nel_file(testo, chiavi_dello_schema())[0] == {
        "design.technical_variables": [(3, True)]}
    assert righe_nel_file(testo) == ({}, {"design": 1})

    # Con --json il documento resta sull'uscita standard, senza le righe.
    assert main(["validate", *comune, "--json"]) == USCITA_RIFIUTO
    catturato = capsys.readouterr()
    documento = json.loads(catturato.out)
    assert all("dove scriverli" not in p["messaggio"] for p in documento["problemi"])
    assert any("Parametro mancante: design.variable" in p["messaggio"]
               for p in documento["problemi"])
    assert "dove scriverli" in catturato.err and "dove scriverli" not in catturato.out


def test_la_riga_di_comando_esegue_le_analisi(eco_r, tmp_path, capsys):
    """
    **Obiettivo**: ``run`` con oggetto, configurazione e cartella di uscita
    esce con 0, scrive le uscite e stampa gli avvisi con il loro codice.

    **Razionale scientifico e sistemistico**: le analisi di base si eseguono
    indicando tre sole cose, come la pipeline si esegue indicando la sua
    configurazione.
    """
    oggetto = costruisci(tmp_path, CONTEGGI, METADATI)
    config = _scrivi_config(tmp_path, configurazione(stat={"min_group_size": 4}))
    uscita = tmp_path / "uscita"
    assert main(["run", "--object", str(oggetto), "--config", str(config),
                 "--out", str(uscita)]) == USCITA_SUCCESSO
    stampato = capsys.readouterr().out
    assert "avviso: [E-ECO-11]" in stampato and "avviso: [E-ECO-12]" in stampato
    assert (uscita / "manifest.json").is_file()
    assert (uscita / "alfa" / "alfa_diversita.tsv").is_file()


def _esegui_py(*argomenti: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [sys.executable, str(RADICE / "scripts" / "esegui.py"), *argomenti],
        capture_output=True, text=True, check=False, cwd=RADICE,
    )


def test_esegui_py_porta_alle_analisi_senza_cambiare_la_pipeline(tmp_path):
    """
    **Obiettivo**: ``scripts/esegui.py eco`` esegue la riga di comando delle
    analisi ecologiche del clone; senza ``eco`` gli argomenti vanno alla
    pipeline come prima, con gli stessi sottocomandi e lo stesso rifiuto di un
    sottocomando sconosciuto.

    **Razionale scientifico e sistemistico**: il comando con cui si esegue la
    pipeline dal clone e' quello documentato per riprodurre i risultati, e non
    deve cambiare significato.
    """
    aiuto_eco = _esegui_py("eco", "--help")
    assert aiuto_eco.returncode == 0
    assert "scripts/esegui.py eco" in aiuto_eco.stdout and "--object" not in aiuto_eco.stdout
    assert "--object" in _esegui_py("eco", "run", "--help").stdout

    aiuto = _esegui_py("--help")
    assert aiuto.returncode == 0
    for sottocomando in ("validate", "run", "resume", "report"):
        assert sottocomando in aiuto.stdout
    assert "amplicon16s_eco" not in aiuto.stdout and "analisi ecologiche" not in aiuto.stdout
    assert "--config" in _esegui_py("run", "--help").stdout
    assert "--object" not in _esegui_py("run", "--help").stdout
    assert _esegui_py("sconosciuto").returncode == 2

    config = _scrivi_config(tmp_path, configurazione(design={"variable": ...}))
    rifiuto = _esegui_py("eco", "validate", "--object", str(tmp_path / "x.rds"),
                         "--config", str(config))
    assert rifiuto.returncode == USCITA_RIFIUTO and "[E-ECO-02]" in rifiuto.stderr


# --------------------------------------------------------------------------- #
# 5. Nessun valore dei dataset nel pacchetto                                   #
# --------------------------------------------------------------------------- #

#: Le configurazioni di esempio dei dataset distribuiti con il repository.
ESEMPI = sorted((RADICE / "dati").glob("*/eco_*.yaml"))

#: Parametri il cui valore e' una voce di un vocabolario chiuso dello schema.
_VOCABOLARI = ("beta.distances", "beta.unifrac_root", "ord.methods", "ord.distances",
               "ord.pcoa_correction")


def _stringhe_dei_dataset() -> set[str]:
    """I valori testuali delle configurazioni di esempio dei dataset: nomi di
    colonna, ranghi, valori di sottoinsieme.
    """
    stringhe: set[str] = set()

    def raccogli(valore) -> None:
        if isinstance(valore, str):
            stringhe.add(valore)
        elif isinstance(valore, dict):
            for chiave, v in valore.items():
                stringhe.add(str(chiave))
                raccogli(v)
        elif isinstance(valore, list):
            for v in valore:
                raccogli(v)

    for esempio in ESEMPI:
        dati = yaml.safe_load(esempio.read_text(encoding="utf-8"))
        for gruppo, valori in dati.items():
            for nome, valore in valori.items():
                if f"{gruppo}.{nome}" not in _VOCABOLARI:
                    raccogli(valore)
    return stringhe


def test_gli_esempi_dei_dataset_sono_configurazioni_valide():
    """
    **Obiettivo**: ogni ``dati/*/eco_*.yaml`` supera la validazione del file e
    dichiara tutti i parametri, anche quelli con un predefinito; fanno
    eccezione i facoltativi che, non dichiarati, non entrano nella
    configurazione risolta, e che un esempio dichiara solo se li usa.

    **Razionale scientifico e sistemistico**: la configurazione di un dataset
    e' il luogo in cui stanno i suoi valori, e deve bastare a se stessa.
    """
    if not ESEMPI:
        pytest.skip("la cartella dati/ non e' disponibile in questo ambiente")
    for esempio in ESEMPI:
        dati = yaml.safe_load(esempio.read_text(encoding="utf-8"))
        config = schema.carica(esempio)
        dichiarati = {f"{g}.{n}" for g, valori in dati.items() for n in valori}
        campi = {
            f"{g}.{n}" for g, modello in schema.ConfigEco.model_fields.items()
            for n in modello.annotation.model_fields
        }
        usati = {c for c in campi
                 if c not in (*schema.CONDIZIONALI, *schema.ASSENTI_SE_NON_DICHIARATI)
                 or getattr(getattr(config, c.split(".")[0]), c.split(".")[1])}
        assert dichiarati == usati, esempio.name


def test_il_pacchetto_non_contiene_valori_dei_dataset():
    """
    **Obiettivo**: nessun valore testuale delle configurazioni di esempio dei
    dataset (nomi di colonna, ranghi, valori di sottoinsieme) compare come
    letterale nel sorgente Python o R del pacchetto, ne' come parola nel
    modello di configurazione; nessun campo dello schema lo ha per predefinito.

    **Razionale scientifico e sistemistico**: il pacchetto deve valere per
    qualunque oggetto consegnato dalla pipeline; i valori di un dataset stanno
    solo nella sua cartella sotto ``dati/``.
    """
    if not ESEMPI:
        pytest.skip("la cartella dati/ non e' disponibile in questo ambiente")
    stringhe = _stringhe_dei_dataset()
    assert len(stringhe) >= 4
    modello = MODELLO.read_text(encoding="utf-8")
    for stringa in stringhe:
        assert not re.search(rf"\b{re.escape(stringa)}\b", modello), stringa
    for percorso in [*PACCHETTO.rglob("*.py"), *PACCHETTO.rglob("*.R")]:
        testo = percorso.read_text(encoding="utf-8")
        for stringa in stringhe:
            for virgolette in "\"'":
                assert f"{virgolette}{stringa}{virgolette}" not in testo, (percorso.name, stringa)
    # I soli predefiniti non vuoti sono valori del metodo, con la fonte, e nessuno
    # e' una stringa.
    for gruppo, modello_gruppo in schema.ConfigEco.model_fields.items():
        for nome, campo in modello_gruppo.annotation.model_fields.items():
            if not campo.is_required():
                predefinito = campo.get_default(call_default_factory=True)
                if predefinito not in (None, (), {}):
                    assert schema.STANDARD_DEL_METODO[f"{gruppo}.{nome}"][0] == predefinito
                    assert not isinstance(predefinito, str)


# --------------------------------------------------------------------------- #
# 6. Ingressi legittimi ma scomodi                                             #
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize("modifica", [
    {"beta": {"clr_pseudocount": float("inf")}},
    {"beta": {"clr_pseudocount": True}},
    {"beta": {"clr_pseudocount": "0.5"}},
    {"stat": {"significance_level": "0.05"}},
    {"comp": {"top_n": 3_000_000_000}},
    {"stat": {"min_group_size": 3_000_000_000}},
    {"stat": {"permanova_permutations": 3_000_000_000}},
    {"alpha": {"rarefy_depth": 3_000_000_000}},
])
def test_i_numeri_non_rappresentabili_sono_respinti_dallo_schema(modifica):
    """
    **Obiettivo**: un reale infinito, un booleano o un testo al posto di un
    numero, e un intero oltre il massimo che R rappresenta, sono respinti con
    ``E-ECO-01``.

    **Razionale scientifico e sistemistico**: oltre lo schema quei valori
    arriverebbero a R come ``NA`` o non attraverserebbero il confine, e
    l'analisi si fermerebbe con un errore che non dice la causa.
    """
    assert _codici(configurazione(**modifica)) == ["E-ECO-01"]


def test_una_chiave_ripetuta_nel_file_e_respinta(tmp_path):
    """
    **Obiettivo**: un file in cui una chiave compare due volte nello stesso
    gruppo e' respinto con ``E-ECO-01`` che la nomina.

    **Razionale scientifico e sistemistico**: il lettore YAML terrebbe in
    silenzio l'ultima dichiarazione, e l'analisi userebbe un valore diverso da
    quello che chi legge il file vede per primo.
    """
    percorso = tmp_path / "eco.yaml"
    testo = yaml.safe_dump(configurazione())
    percorso.write_text(testo.replace("variable: gruppo", "variable: gruppo\n  variable: lotto"),
                        encoding="utf-8")
    with pytest.raises(ErroreEco) as rifiuto:
        schema.carica(percorso)
    assert rifiuto.value.codici == ["E-ECO-01"] and "variable" in str(rifiuto.value)


def test_le_colonne_con_nomi_non_sintattici_si_dichiarano_con_il_loro_nome(eco_r, tmp_path):
    """
    **Obiettivo**: colonne i cui nomi hanno spazi, cifre iniziali o segni si
    dichiarano con il nome che hanno nell'oggetto, e un nome sbagliato e'
    respinto con l'elenco dei nomi veri; gli spazi attorno ai valori non
    creano gruppi distinti.

    **Razionale scientifico e sistemistico**: i metadati di uno studio hanno
    nomi come «Factor Value[...]»; R li riscriverebbe, e la configurazione
    dovrebbe usare nomi che nell'oggetto non esistono.
    """
    metadati = {
        "Factor Value[Sito n°1]": ["a", "a ", " a", "b", "b", "b"],
        "1": ["x", "y", "x", "y", "x", "y"],
    }
    oggetto = costruisci(tmp_path, CONTEGGI, metadati)
    esito = valida_sull_oggetto(oggetto, schema.valida(configurazione(
        design={"variable": "Factor Value[Sito n°1]", "technical_variables": ["1"]})))
    (variabile,) = esito.riepilogo["variabili"]
    assert variabile["nome"] == "Factor Value[Sito n°1]"
    assert variabile["gruppi"] == {"a": 3, "b": 3}
    rifiuto = _rifiuto(oggetto, design={"variable": "Factor.Value.Sito.n.1."})
    assert rifiuto.codici == ["E-ECO-05"] and "Factor Value[Sito n°1], 1" in str(rifiuto)


def test_una_variabile_senza_valori_non_ferma_l_analisi(eco_r, tmp_path):
    """
    **Obiettivo**: una variabile con tutti i valori mancanti e' dichiarata
    (``E-ECO-13``, ``E-ECO-12``) e l'analisi si conclude con le uscite
    descrittive, senza test per quella variabile.

    **Razionale scientifico e sistemistico**: una colonna vuota fra le altre
    variabili non deve impedire l'analisi delle altre.
    """
    metadati = {**METADATI, "vuota": [None, "", "NA", "n/a", None, ""]}
    oggetto = costruisci(tmp_path, CONTEGGI, metadati)
    config = schema.valida(configurazione(design={"other_variables": ["vuota"]}))
    esito = esegui(oggetto, config, tmp_path / "uscita")
    codici = [a.codice for a in esito.avvisi if "'vuota'" in a.dettaglio]
    assert codici == ["E-ECO-13", "E-ECO-12"]
    vuota = esito.riepilogo["variabili"][1]
    assert vuota["campioni_con_valore"] == 0 and vuota["gruppi_nei_test"] == []
    assert (tmp_path / "uscita" / "test" / "permanova.tsv").is_file()


def test_conteggi_che_non_sono_conteggi_e_oggetti_minimi_sono_respinti(eco_r, tmp_path):
    """
    **Obiettivo**: un oggetto con una sola variante e' respinto con
    ``E-ECO-04`` e il motivo; lo stesso vale, per costruzione del controllo,
    per conteggi negativi o non interi.

    **Razionale scientifico e sistemistico**: rarefazione e indici sono
    definiti su conteggi di letture; proporzioni o valori trasformati darebbero
    numeri senza significato o un errore di R senza causa.
    """
    una = costruisci(tmp_path, [[5, 6, 7, 8, 9, 10]], METADATI,
                     tassonomia={"Regno": ["Batteri"], "Famiglia": ["F"], "Genere": ["G"]})
    rifiuto = _rifiuto(una)
    assert rifiuto.codici == ["E-ECO-04"] and "1 varianti" in str(rifiuto)
    sorgente = (PACCHETTO / "R" / "lib" / "campioni.R").read_text(encoding="utf-8")
    assert "valori negativi o non interi" in sorgente and "valori mancanti o non numerici" in sorgente


def test_una_cartella_di_uscita_non_utilizzabile_e_respinta_prima_del_calcolo(eco_r, tmp_path):
    """
    **Obiettivo**: un collegamento che non porta a nulla e un percorso sotto
    un file sono respinti con ``E-ECO-09`` prima del calcolo; un rifiuto in
    validazione non lascia la cartella di uscita appena creata.

    **Razionale scientifico e sistemistico**: un percorso non scrivibile deve
    emergere subito e con un codice, non come errore del programma a calcolo
    concluso.
    """
    oggetto = costruisci(tmp_path, CONTEGGI, METADATI)
    config = schema.valida(configurazione())
    pendente = tmp_path / "pendente"
    pendente.symlink_to(tmp_path / "non_esiste")
    file = tmp_path / "file.txt"
    file.write_text("x", encoding="utf-8")
    for uscita in (pendente, file / "sotto"):
        with pytest.raises(ErroreEco) as rifiuto:
            esegui(oggetto, config, uscita)
        assert rifiuto.value.codici == ["E-ECO-09"]
    nuova = tmp_path / "nuova" / "annidata"
    with pytest.raises(ErroreEco) as rifiuto:
        esegui(oggetto, schema.valida(configurazione(comp={"rank": "Specie"})), nuova)
    assert rifiuto.value.codici == ["E-ECO-07"] and not nuova.exists()


def test_eco_init_rileva_propone_e_lascia_da_scegliere(eco_r, tmp_path, capsys):
    """
    **Obiettivo**: ``init`` sull'oggetto scrive una configurazione di partenza
    in cui nessuna riga e' attiva: riporta come fatti i ranghi con la quota di
    varianti assegnate e l'albero con la sua radice; propone la variabile
    biologica solo se una sola colonna ha gruppi ripetuti (con piu' candidate
    le elenca e lascia la scelta senza valore), come variabili tecniche le
    colonne di piastra e di corsa che hanno piu' di un valore (la piu'
    grossolana per prima), come rango il piu' fine con almeno meta' delle
    varianti assegnate, e la radice esistente per UniFrac; lascia da scegliere
    distanze, metodi, seme, dimensione dei gruppi e significativita'. Il file
    cosi' com'e' e' respinto nominando i parametri mancanti; completato, la
    validazione sull'oggetto passa. Un file esistente non viene riscritto.

    **Razionale scientifico e sistemistico**: quale variabile lo studio
    confronta, e con quali metodi, lo decide chi conduce l'analisi: ``init``
    porta cio' che l'oggetto contiene e un punto di partenza, e nessun valore
    entra nella configurazione senza un gesto.
    """
    from init_aiuti import completa, righe_dei_parametri

    from amplicon16s_eco.init import QUOTA_MINIMA_ASSEGNATA

    conteggi = [[10 + (3 * i + j) % 9 for j in range(12)] for i in range(4)]
    metadati = {
        "accession": [f"c{i}" for i in range(12)],
        "trattamento": ["a", "b", "c"] * 4,
        "sito": ["x"] * 6 + ["y"] * 6,
        "piastra": ["p1", "p2", "p3", "p4"] * 3,
        "corsa": ["r1"] * 6 + ["r2"] * 6,
        "misura": [str(i * 1.5) for i in range(12)],
    }
    tassonomia = {
        "Regno": ["Batteri"] * 4,
        "Famiglia": ["F1", "F1", "F2", None],
        "Genere": ["G1", None, None, None],
    }
    oggetto = costruisci(tmp_path, conteggi, metadati, tassonomia,
                         albero="((v01:1,v02:1):1,(v03:1,v04:1):1);")
    uscita = tmp_path / "eco.yaml"
    assert main(["init", "--object", str(oggetto), "--output", str(uscita)]) == USCITA_SUCCESSO
    stampato = capsys.readouterr().out
    testo = uscita.read_text(encoding="utf-8")
    # I fatti stanno nell'intestazione, in righe di commento mandate a capo.
    intestazione = " ".join(testo.replace("#", " ").split())
    assert "Regno 1.00, Famiglia 0.75, Genere 0.25" in intestazione
    assert "presente, radicato" in intestazione and "Regno 1.00" in stampato
    assert "Da scegliere:" in stampato and "design.variable" in stampato

    righe = righe_dei_parametri(testo, schema.ConfigEco)
    assert not [chiave for chiave, (attiva, _) in righe.items() if attiva]
    assert righe["design.variable"] == (False, "")
    assert "sito (2), trattamento (3)" in intestazione
    assert righe["design.technical_variables"] == (False, '["corsa", "piastra"]')
    assert QUOTA_MINIMA_ASSEGNATA == 0.5 and righe["comp.rank"] == (False, '"Famiglia"')
    assert righe["beta.unifrac_root"] == (False, '"existing"')
    for chiave in ("run.seed", "beta.distances", "ord.methods", "ord.distances",
                   "stat.min_group_size", "stat.significance_level", "comp.top_n"):
        assert righe[chiave] == (False, ""), chiave
    assert righe["ord.nmds_trymax"] == (False, "20") and "vegan::metaMDS" in testo
    assert righe["stat.permanova_permutations"] == (False, "999")

    # Cosi' com'e' non e' una configurazione: mancano tutti gli obbligatori.
    assert main(["validate", "--object", str(oggetto), "--config", str(uscita)]) == USCITA_RIFIUTO
    rifiuto = capsys.readouterr().err
    assert all(f"Parametro mancante: {chiave}" in rifiuto for chiave in schema.OBBLIGATORI)

    completo = tmp_path / "completo.yaml"
    completo.write_text(completa(
        testo,
        {"run.seed": 3, "comp.top_n": 3, "beta.distances": ["bray", "unifrac_weighted"],
         "ord.methods": ["pcoa"], "ord.distances": ["bray"], "stat.min_group_size": 2,
         "stat.significance_level": 0.05, "design.variable": "trattamento"},
        accetta={"design.technical_variables", "comp.rank",
                 "beta.unifrac_root", "ord.pcoa_correction"},
        modello=schema.ConfigEco), encoding="utf-8")
    config = schema.carica(completo)
    assert (config.design.variable, config.comp.rank) == ("trattamento", "Famiglia")
    assert valida_sull_oggetto(oggetto, config).riepilogo["campioni"]["analizzati"] == 12

    prima = uscita.read_bytes()
    assert main(["init", "--object", str(oggetto), "--output", str(uscita)]) == USCITA_RIFIUTO
    assert uscita.read_bytes() == prima

    # Con una sola colonna a gruppi ripetuti la variabile e' una proposta.
    (tmp_path / "sola").mkdir()
    oggetto = costruisci(tmp_path / "sola", conteggi,
                         {k: v for k, v in metadati.items() if k != "sito"}, tassonomia)
    sola = tmp_path / "sola.yaml"
    assert main(["init", "--object", str(oggetto), "--output", str(sola)]) == USCITA_SUCCESSO
    righe = righe_dei_parametri(sola.read_text(encoding="utf-8"), schema.ConfigEco)
    assert righe["design.variable"] == (False, '"trattamento"')


def test_eco_init_accetta_proposte_attiva_proposte_e_partenze_non_i_condizionali():
    """
    **Obiettivo**: Verificare, sulla descrizione di un oggetto con una sola
    variabile a gruppi ripetuti e un albero radicato, che con
    ``--accetta-proposte`` siano righe attive la variabile biologica, le
    variabili tecniche e il rango proposti, con la loro origine; che i
    parametri condizionali che hanno un valore (la radice per UniFrac, la
    correzione della PCoA, gli avvii della NMDS) restino in commento con quel
    valore e con la nota che lo spiega; che le scelte senza valore restino in
    commento; e che senza l'opzione nessuna riga sia attiva.

    **Razionale scientifico e sistemistico**: Un condizionale e' respinto se
    il metodo che lo richiede non e' fra quelli scelti: attivarlo prima che
    distanze e ordinazioni siano scelte darebbe una configurazione che la
    validazione rifiuta. Quali metodi usare e' una scelta di chi analizza.
    """
    from init_aiuti import righe_dei_parametri

    from amplicon16s_eco import init

    descrizione = {
        "campioni": 12, "varianti": 4,
        "ranghi": [{"rango": "Regno", "quota_assegnata": 1.0},
                   {"rango": "Famiglia", "quota_assegnata": 0.75},
                   {"rango": "Genere", "quota_assegnata": 0.25}],
        "albero": {"presente": True, "radicato": True, "con_lunghezze": True},
        "colonne": [
            {"nome": "trattamento", "distinti": 3, "gruppi_ripetuti": 3, "numerica": False,
             "con_valore": 12},
            {"nome": "piastra", "distinti": 4, "gruppi_ripetuti": 4, "numerica": False,
             "con_valore": 12},
            {"nome": "corsa", "distinti": 2, "gruppi_ripetuti": 2, "numerica": False,
             "con_valore": 12},
        ],
    }
    voci, fatti = init.prepara(descrizione)
    senza = righe_dei_parametri(init.testo_yaml(voci, fatti, "ps.rds"), schema.ConfigEco)
    assert not [c for c, (attiva, _) in senza.items() if attiva]

    testo = init.testo_yaml(init.accetta(voci), fatti, "ps.rds", True)
    righe = righe_dei_parametri(testo, schema.ConfigEco)
    attive = {c for c, (attiva, _) in righe.items() if attiva}
    assert attive == {"design.variable", "design.technical_variables", "comp.rank"}
    assert righe["design.variable"] == (True, '"trattamento"')
    assert righe["comp.rank"] == (True, '"Famiglia"')
    for chiave, valore in (("beta.unifrac_root", '"existing"'), ("ord.pcoa_correction", '"none"'),
                           ("ord.nmds_trymax", "20")):
        assert righe[chiave] == (False, valore), chiave
    assert testo.count("non attivato da --accetta-proposte") == 3
    assert "accettato con --accetta-proposte" in testo
    # Nessuna riga di parametro porta [PARTENZA]: i soli valori di partenza
    # delle analisi sono condizionali, e restano in commento.
    assert not re.search(r"^  ## \[PARTENZA\]", testo, re.M)
    # Una scelta non condizionale con un valore di partenza diventerebbe attiva.
    from amplicon16s.init.scrittura import Specie as SpecieInit
    from amplicon16s.init.scrittura import Voce

    finta = init.accetta({"run.seed": Voce("run.seed", SpecieInit.DA_SCEGLIERE, 7, "prova")})
    assert finta["run.seed"].specie is SpecieInit.PARTENZA and finta["run.seed"].attiva
    for chiave in schema.OBBLIGATORI:
        if chiave not in attive:
            assert righe[chiave][0] is False, chiave
    dati = yaml.safe_load(testo)
    assert dati["design"]["variable"] == "trattamento" and dati["comp"] == {"rank": "Famiglia"}


def test_i_condizionali_dei_metodi_gia_dichiarati_si_chiedono_nello_stesso_giro():
    """
    **Obiettivo**: Verificare che, quando la configurazione ha gia' l'elenco
    delle distanze e dei metodi di ordinazione ma le manca un altro parametro
    obbligatorio, il rifiuto nomini insieme a quello anche i condizionali che i
    metodi dichiarati richiedono (lo pseudoconteggio con Aitchison, la
    correzione con la PCoA); che un condizionale gia' respinto dallo schema per
    il suo valore non risulti anche mancante; che con l'elenco dei metodi non
    ancora dichiarato nessun condizionale sia chiesto; e che il testo «che
    cosa fare» di un parametro mancante valga per ogni parametro.

    **Razionale scientifico e sistemistico**: Chi completa la configurazione
    deve vedere in un solo giro tutto cio' che manca e che si puo' gia'
    conoscere. Un condizionale dipende dai metodi scelti: senza l'elenco non
    lo si puo' chiedere senza inventare una scelta.
    """
    from amplicon16s_eco.catalogo import CATALOGO

    completa = configurazione()
    completa["beta"]["distances"] = ["bray", "aitchison"]
    completa["ord"] = {"methods": ["pcoa"], "distances": ["bray"]}
    completa["beta"].pop("clr_pseudocount", None)
    del completa["stat"]["significance_level"]
    with pytest.raises(ErroreEco) as rifiuto:
        schema.valida(completa)
    testo = str(rifiuto.value)
    assert rifiuto.value.codici == ["E-ECO-02"] * 3
    for parametro in ("stat.significance_level", "beta.clr_pseudocount", "ord.pcoa_correction"):
        assert parametro in testo

    # Un valore non ammesso e' un solo problema, non anche un parametro mancante.
    completa["beta"]["clr_pseudocount"] = -1
    with pytest.raises(ErroreEco) as rifiuto:
        schema.valida(completa)
    assert sorted(rifiuto.value.codici) == ["E-ECO-01", "E-ECO-02", "E-ECO-02"]

    # Senza l'elenco dei metodi i suoi condizionali non si conoscono.
    senza = configurazione()
    del senza["beta"]["distances"]
    senza["ord"] = {"methods": [], "distances": []}
    with pytest.raises(ErroreEco) as rifiuto:
        schema.valida(senza)
    assert "clr_pseudocount" not in str(rifiuto.value)

    assert "design.variable" not in CATALOGO["E-ECO-02"].azione


def test_eco_init_non_scrive_un_valore_per_un_condizionale_senza_proposta():
    """
    **Obiettivo**: Verificare che, per un oggetto senza albero, la riga di
    ``beta.unifrac_root`` scritta da ``eco init`` sia in commento e senza
    valore, con la condizione che la richiede, e non riporti il nullo dello
    schema.

    **Razionale scientifico e sistemistico**: Il nullo e' l'assenza del
    parametro, non una scelta: scritto accanto al nome fa pensare a un valore
    da confermare, mentre senza albero il parametro non va dichiarato.
    """
    from init_aiuti import righe_dei_parametri

    from amplicon16s_eco import init

    descrizione = {
        "campioni": 6, "varianti": 4,
        "ranghi": [{"rango": "Regno", "quota_assegnata": 1.0}],
        "albero": {"presente": False, "radicato": False, "con_lunghezze": False},
        "colonne": [{"nome": "trattamento", "distinti": 2, "gruppi_ripetuti": 2,
                     "numerica": False, "con_valore": 6}],
    }
    voci, fatti = init.prepara(descrizione)
    testo = init.testo_yaml(voci, fatti, "ps.rds")
    assert righe_dei_parametri(testo, schema.ConfigEco)["beta.unifrac_root"] == (False, "")
    assert "unifrac_root: null" not in testo
    assert re.search(r"\[CONDIZIONALE\] solo se beta.distances comprende[^\n]*\n(?:  ##[^\n]*\n)*"
                     r"  # unifrac_root:\n", testo)


def test_eco_init_propone_come_continua_la_colonna_che_sembra_una_misura():
    """
    **Obiettivo**: Verificare, sulla descrizione di un oggetto con una colonna
    numerica dai valori quasi tutti diversi, che ``eco init`` la proponga in
    ``design.continuous_variables`` e non fra le candidate a variabile
    categorica; che una colonna numerica con pochi valori ripetuti (una dose a
    due livelli) resti candidata categorica; che senza colonne di quel tipo il
    parametro resti la riga facoltativa del modello; e che il soggetto non sia
    mai proposto.

    **Razionale scientifico e sistemistico**: una misura continua trattata
    come categorica diventa un gruppo per campione. La regola della proposta e'
    quella dell'avviso sulle variabili che sembrano continue: numerica, almeno
    dieci valori distinti, distinti in almeno quattro quinti dei campioni. Chi
    sia il soggetto di un disegno a misure ripetute lo sa solo chi lo ha fatto.
    """
    from init_aiuti import righe_dei_parametri

    from amplicon16s_eco import init

    def descrizione(*colonne):
        return {
            "campioni": 12, "varianti": 4,
            "ranghi": [{"rango": "Regno", "quota_assegnata": 1.0}],
            "albero": {"presente": False, "radicato": False, "con_lunghezze": False},
            "colonne": list(colonne),
        }

    gruppo = {"nome": "trattamento", "distinti": 2, "gruppi_ripetuti": 2, "numerica": False,
              "con_valore": 12}
    dose = {"nome": "dose", "distinti": 2, "gruppi_ripetuti": 2, "numerica": True,
            "con_valore": 12}
    misura = {"nome": "distanza_m", "distinti": 11, "gruppi_ripetuti": 1, "numerica": True,
              "con_valore": 12}
    voci, fatti = init.prepara(descrizione(gruppo, misura))
    assert voci["design.continuous_variables"].valore == ["distanza_m"]
    # Tutte le colonne dell'oggetto sono elencate, con cio' che serve a scegliere.
    assert any("trattamento (2)" in f and "distanza_m (11, numerica)" in f for f in fatti)
    assert voci["design.variable"].valore == "trattamento"
    testo = init.testo_yaml(voci, fatti, "ps.rds")
    assert righe_dei_parametri(testo, schema.ConfigEco)["design.continuous_variables"] == (
        False, '["distanza_m"]')
    accettato = init.testo_yaml(init.accetta(voci), fatti, "ps.rds", True)
    assert righe_dei_parametri(accettato, schema.ConfigEco)["design.continuous_variables"][0]
    assert righe_dei_parametri(accettato, schema.ConfigEco)["design.subject"] == (False, "null")

    # Una dose a due livelli e' numerica ma non una misura: resta categorica.
    voci, _ = init.prepara(descrizione(gruppo, dose))
    assert "design.continuous_variables" not in voci
    assert voci["design.variable"].valore is not None
    assert "design.subject" not in voci


def _consegna(cartella: Path, letture: dict[str, int], colonne: dict[str, list[str]]) -> Path:
    """Scrive in ``cartella`` un oggetto finto con gli export in testo e il file
    dei checksum, come la pipeline li consegna; restituisce il percorso
    dell'oggetto. ``eco init`` non apre l'oggetto per leggere gli export: ne
    verifica l'impronta.
    """
    from amplicon16s.io_layer.checksums import checksum_file

    cartella.mkdir(parents=True, exist_ok=True)
    campioni = list(letture)
    (cartella / "ps_final.rds").write_bytes(b"oggetto di prova")
    (cartella / "conteggi.tsv").write_text(
        "asv_id\t" + "\t".join(campioni) + "\n"
        + "ASV1\t" + "\t".join(str(letture[c] - letture[c] // 2) for c in campioni) + "\n"
        + "ASV2\t" + "\t".join(str(letture[c] // 2) for c in campioni) + "\n", encoding="utf-8")
    (cartella / "metadati.tsv").write_text(
        "accession\tsample_name\t" + "\t".join(colonne) + "\n" + "".join(
            f"{c}\tcampione {c}\t" + "\t".join(v[i] for v in colonne.values()) + "\n"
            for i, c in enumerate(campioni)), encoding="utf-8")
    (cartella / "checksum.sha256").write_text("".join(
        f"{checksum_file(cartella / n).removeprefix('sha256:')}  {n}\n"
        for n in ("ps_final.rds", "conteggi.tsv", "metadati.tsv")), encoding="utf-8")
    return cartella / "ps_final.rds"


def _descrizione(colonne: dict[str, list[str]]) -> dict:
    """La descrizione che ``descrivi`` darebbe di una tabella dei campioni."""
    def una(nome: str, valori: list[str]) -> dict:
        presenti = [v for v in valori if v]
        return {"nome": nome, "con_valore": len(presenti), "distinti": len(set(presenti)),
                "gruppi_ripetuti": sum(presenti.count(v) >= 2 for v in set(presenti)),
                "numerica": False}
    return {"campioni": len(next(iter(colonne.values()))), "varianti": 2,
            "ranghi": [{"rango": "Regno", "quota_assegnata": 1.0}],
            "albero": {"presente": False, "radicato": None, "con_lunghezze": None},
            "colonne": [una(n, v) for n, v in colonne.items()]}


def test_eco_init_avvisa_se_la_rarefazione_scende_molto_sotto_la_mediana(tmp_path):
    """
    **Obiettivo**: Verificare che, quando il campione piu' povero ha meno di un
    decimo della mediana delle letture, ``eco init`` lo dica fra i fatti
    rilevati (in testa al file e nel riepilogo), con la profondita' a cui le
    analisi rarefarebbero, i campioni che la abbassano nominati dal piu'
    povero con le loro letture, e il rimedio (alzare ``qc.min_reads_final`` e
    ``resume``, o ``design.subset``); che con campioni di profondita'
    paragonabile l'avviso non compaia; e che gli export si usino solo se il
    file dei checksum li dichiara per questo oggetto: con l'oggetto o un
    export cambiati dopo la consegna, o senza export, l'avviso non si valuta e
    il riepilogo lo dichiara.

    **Razionale scientifico e sistemistico**: Senza una profondita' dichiarata
    le analisi rarefanno al campione piu' povero: un campione fallito porta
    tutti gli altri a poche centinaia di letture, e chi legge i grafici non lo
    vede. I conteggi in testo valgono per l'oggetto solo se la loro impronta e'
    quella consegnata con lui: altrimenti si dichiara di non sapere.
    """
    from amplicon16s_eco import init

    colonne = {"gruppo": ["a", "a", "a", "b", "b", "b", "b"]}
    letture = {"c1": 60000, "c2": 250000, "c3": 357, "c4": 1584, "c5": 102, "c6": 44000,
               "c7": 8700}
    oggetto = _consegna(tmp_path / "poveri", letture, colonne)
    trovati = init.esportati(oggetto)
    assert trovati is not None and trovati.letture == letture
    voci, fatti = init.prepara(_descrizione(colonne), trovati)
    avviso = next(f for f in fatti if f.startswith("Rarefazione:"))
    assert "la profondita' sarebbe 102 letture" in avviso and "mediana (8700)" in avviso
    assert "La abbassano: campione c5 (102), campione c3 (357)." in avviso
    assert "c4" not in avviso
    assert "alzare qc.min_reads_final a 1000" in avviso
    assert "design.subset" in avviso and "resume" in avviso
    assert "non siano un intero gruppo dello studio" in avviso
    assert any(f.startswith("Gruppi delle colonne candidate") and "gruppo: 'a' 3, 'b' 4" in f
               for f in fatti)
    # Con almeno mille letture nel campione piu' povero il consiglio e' tenerlo.
    profondi = init.rarefazione_bassa({"c1": 9000, "c2": 150000, "c3": 160000, "c4": 170000})
    assert "di solito si lascia cosi'" in profondi and "qc.min_reads_final" not in profondi
    testo = " ".join(init.testo_yaml(voci, fatti, "ps_final.rds").replace("#", " ").split())
    assert "la profondita' sarebbe 102 letture" in testo

    pari = _consegna(tmp_path / "pari", {"c1": 9000, "c2": 12000, "c3": 30000, "c4": 4000},
                     {"gruppo": ["a", "a", "b", "b"]})
    _, fatti = init.prepara(_descrizione({"gruppo": ["a", "a", "b", "b"]}), init.esportati(pari))
    assert not [f for f in fatti if "Rarefazione" in f or "non valutate" in f]

    # Un export cambiato dopo la consegna, l'oggetto sostituito, nessun export.
    (tmp_path / "pari" / "conteggi.tsv").write_text("asv_id\tc1\nASV1\t5\n", encoding="utf-8")
    assert init.esportati(pari) is None
    oggetto.write_bytes(b"un altro oggetto")
    assert init.esportati(oggetto) is None
    (tmp_path / "solo").mkdir()
    (tmp_path / "solo" / "ps.rds").write_bytes(b"oggetto")
    assert init.esportati(tmp_path / "solo" / "ps.rds") is None
    _, fatti = init.prepara(_descrizione(colonne), None)
    assert any(f.startswith("Profondita' di rarefazione e variabili confuse: non valutate")
               for f in fatti)


def test_eco_init_avvisa_sulle_variabili_confuse(tmp_path):
    """
    **Obiettivo**: Verificare che ``eco init`` riconosca fra le colonne
    candidate a variabile biologica: due colonne che dividono i campioni negli
    stessi gruppi (la stessa variabile, detta una volta); una colonna i cui
    gruppi cadono ciascuno in un solo gruppo dell'altra (confuse); e il caso
    parziale, in cui un gruppo intero di una cade in un solo gruppo dell'altra
    senza che valga per tutti (i campioni a terra sono tutti di un tessuto, e
    un tessuto e' tutto in volo); che l'avviso stia nel riepilogo e nel
    commento di ``design.variable``, con il rimedio; e che due colonne
    incrociate, in cui ogni gruppo dell'una ha campioni di piu' gruppi
    dell'altra, non diano alcun avviso.

    **Razionale scientifico e sistemistico**: Se tutti i campioni di una
    condizione sono anche di un solo tessuto, il confronto fra le condizioni
    misura insieme la differenza fra i tessuti, e nessun test le separa: chi
    sceglie la variabile deve saperlo prima di leggere un valore di p, e i
    metadati da soli non lo mostrano.
    """
    from amplicon16s_eco import init

    colonne = {
        "condizione": ["volo", "volo", "volo", "volo", "volo", "volo", "terra", "terra"],
        "tessuto": ["radice", "radice", "foglia", "foglia", "foglia", "foglia", "radice",
                    "radice"],
        "materiale": ["r", "r", "f", "f", "f", "f", "r", "r"],
        "replica": ["1", "2", "1", "2", "1", "2", "1", "2"],
    }
    letture = {f"c{i}": 5000 + i for i in range(8)}
    trovati = init.esportati(_consegna(tmp_path / "parziale", letture, colonne))
    voci, fatti = init.prepara(_descrizione(colonne), trovati)
    avvisi = [f.removeprefix("Variabili confuse: ") for f in fatti
              if f.startswith("Variabili confuse: ")]
    assert len(avvisi) == 2
    assert avvisi[0].startswith("materiale e tessuto dividono i campioni negli stessi gruppi")
    assert avvisi[1].startswith(
        "condizione e materiale sono confuse in parte: i campioni 'terra' di condizione sono "
        "tutti 'r' di materiale; i campioni 'f' di materiale sono tutti 'volo' di condizione.")
    assert "design.subset" in avvisi[1]
    # La replica incrocia ogni altra colonna: nessun avviso la nomina.
    assert not [a for a in avvisi if "replica" in a]
    note = " ".join(voci["design.variable"].note)
    assert "attenzione: condizione e materiale sono confuse in parte" in note
    testo = " ".join(init.testo_yaml(voci, fatti, "ps_final.rds").replace("#", " ").split())
    assert testo.count("condizione e materiale sono confuse in parte") == 2

    # Ogni gruppo di una colonna in un solo gruppo dell'altra: confuse del tutto.
    annidate = {"soggetto": ["s1", "s1", "s2", "s2", "s3", "s3", "s4", "s4"],
                "dieta": ["x", "x", "x", "x", "y", "y", "y", "y"]}
    assert init.confuse(
        [dict(zip(annidate, valori)) for valori in zip(*annidate.values())], list(annidate)
    )[0].startswith("soggetto e dieta sono confuse: ogni gruppo di soggetto cade in un solo "
                    "gruppo di dieta.")
    # Un valore assente non conta come gruppo.
    con_assenti = [{"a": "p", "b": "u"}, {"a": "p", "b": "v"}, {"a": "q", "b": "u"},
                   {"a": "q", "b": "v"}, {"a": "NA", "b": "u"}, {"a": "NA", "b": "u"}]
    assert init.confuse(con_assenti, ["a", "b"]) == []
