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
    dichiara tutti i parametri, anche quelli con un predefinito.

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
                 if c not in schema.CONDIZIONALI
                 or getattr(getattr(config, c.split(".")[0]), c.split(".")[1]) is not None}
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
