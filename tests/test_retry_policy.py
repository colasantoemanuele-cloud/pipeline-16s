"""Suite di verifica della politica di retry automatico, delle degradazioni, del punto di ripresa e della CLI.

Inquadramento nel Piano Operativo
---------------------------------
* **Settimane di riferimento**: **Settimana 9 e Settimana 10 (W9/W10 — Fase F3:
  Politica di retry automatico, aggiustamenti correttivi, punto di ripresa e
  CLI ``run`` / ``resume`` / ``validate`` / ``report``)**.
* **Moduli sorgente coperti**:
  - ``src/amplicon16s/runner/retry.py``
  - ``src/amplicon16s/runner/executor.py``
  - ``src/amplicon16s/runner/project.py``
  - ``src/amplicon16s/cli.py``
  - ``src/amplicon16s/config/resolve.py``

Scopo sperimentale e razionale scientifico/sistemistico
-------------------------------------------------------
Questo modulo verifica il comportamento dell'``Esecutore`` e dell'interfaccia a
riga di comando ``amplicon16s`` davanti a guasti transitori, arresti critici,
degradazioni controllate e riprese successive:

1. **Chiusura ermetica del retry ai 4 codici tecnici della whitelist
   (``E-S2-03``, ``E-S3-01``, ``E-S4-02``, ``E-S5-01``)**: il retry automatico è
   consentito esclusivamente dove l'azione correttiva (dimezzare ``run.batch_size``
   per ridurre l'impronta RAM o raddoppiare ``err.nbases`` per far convergere
   ``learnErrors``) agisce su parametri operativi/computazionali senza mai
   alterare soglie biologiche (come ``filter.maxEE`` o ``decontam.threshold``).
2. **Tracciabilità degli aggiustamenti**: ogni modifica applicata da un retry
   viene annotata sia in ``manifest_<passo>.json`` (con valore ``dichiarato`` e
   valore ``usato``) sia nel log strutturato ``99_logs/pipeline.jsonl``, senza
   alterare il file ``config.yaml`` dell'utente né la configurazione delle fasi
   successive.
3. **Punto di ripresa persistente (``punto_di_ripresa.json`` e ``.txt``)**:
   quando i tentativi si esauriscono o interviene un errore a revisione umana,
   l'esecutore salva in ``99_logs/`` la diagnosi completa e il comando esatto
   ``amplicon16s resume --config ...`` da eseguire dopo aver rimosso la causa.
4. **Degradazioni automatiche controllate (``E-S6-02``, ``E-S11-02``,
   ``E-S0-15``, ``E-S1-01``)**: verifica che le condizioni di ripiego previste
   dal protocollo registrino l'avviso nel manifesto e proseguano l'esecuzione.
5. **Regole di sicurezza della CLI e versionamento di ``resolved.yaml``**:
   - ``run`` rifiuta di sovrascrivere una cartella di output già esistente o non vuota;
   - ``resume`` riparte dalla prima fase non conclusa senza ripetere il lavoro a
     monte e, se cambiano parametri non computazionali (es. ``run.threads`` o
     ``retry.max_attempts``), affianca una nuova versione progressiva
     (``resolved_2.yaml``, ``resolved_3.yaml``, ...) con il diff esplicito
     ``differenze_dalla_precedente`` senza sovrascrivere ``resolved.yaml``;
   - ``validate`` esegue solo la Fase S0; ``report`` sintetizza stato,
     aggiustamenti e degradazioni.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, ClassVar

import pytest
import yaml
from conftest import NEGATIVO, POSITIVO, Campione, crea_scenario

import amplicon16s.cli as cli
from amplicon16s.config.resolve import risolvi, versioni_registrate
from amplicon16s.config.schema import Config, valida
from amplicon16s.errors.exceptions import errore
from amplicon16s.gates.g01_g15 import ErroreGate
from amplicon16s.io_layer.artifacts import Fase
from amplicon16s.logging.logger import NOME_FILE_LOG, chiudi, configura
from amplicon16s.runner.executor import Conclusione, Esecutore
from amplicon16s.runner.graph import GRAFO, Passo
from amplicon16s.runner.project import ProjectRun, StatoPasso
from amplicon16s.runner.retry import (
    Aggiustamento,
    PoliticaRetry,
    dimezza,
    raddoppia,
)
from amplicon16s.steps.base import PipelineStep, Produzione, StepContext
from amplicon16s.steps.s00_validate import ValidazioneIngressi

TUTTE = tuple(Passo)


# --------------------------------------------------------------------------- #
# Scenario e doppioni                                                          #
# --------------------------------------------------------------------------- #


def _campioni() -> list[Campione]:
    return [
        Campione("ERX3000001", "NOD1D4.L1", piastra="1"),
        Campione("ERX3000002", "NOD1D4.L2", piastra="1"),
        Campione("ERX3000003", "POS.P1.1", materiale=POSITIVO,
                 posizione="Not Applicable", piastra="1"),
        Campione("ERX3000004", "BLANK.P1.1", materiale=NEGATIVO,
                 posizione="Not Applicable", piastra="1"),
    ]


def _variante(config: Config, **valori: Any) -> Config:
    dati = config.model_dump(mode="python")
    for nome, v in valori.items():
        gruppo, chiave = nome.split("__")
        dati[gruppo][chiave] = v
    return valida(dati)


class Doppione(PipelineStep):
    """Scrive un artefatto e annota i parametri che ha visto."""

    #: Per i fragili: codice sollevato e quanti tentativi falliscono.
    codice: ClassVar[str | None] = None
    fallimenti: ClassVar[int] = 0
    #: "degrada": registra una degradazione e prosegue; "solleva": solleva
    #: una DegradazioneRichiesta; con ``ripiego`` la fase sa proseguire.
    degradazione: ClassVar[str | None] = None
    codice_degradazione: ClassVar[str] = "E-S6-02"
    ripiego: ClassVar[bool] = False

    def __init__(self, registro: list[tuple[Passo, int, float]]) -> None:
        self.registro = registro
        self.tentativi = 0

    def _scrivi(self, contesto: StepContext, testo: str = "") -> Produzione:
        artefatto = contesto.albero.scrivi_testo(
            self.cartella, f"{str(self.passo).lower()}.txt", f"{self.passo} {testo}"
        )
        return Produzione((artefatto,))

    def calcola(self, contesto: StepContext) -> Produzione:
        config = contesto.config
        self.registro.append((self.passo, config.run.batch_size, config.err.nbases))
        self.tentativi += 1
        if self.codice and self.tentativi <= self.fallimenti:
            raise errore(self.codice, f"doppione, tentativo {self.tentativi}")
        if self.degradazione == "degrada":
            contesto.degrada(self.codice_degradazione, "frazione chimerica 0,31")
        elif self.degradazione == "solleva":
            raise errore(self.codice_degradazione, "soglia non attendibile")
        return self._scrivi(contesto)

    def ripiega(self, contesto, degradazione):
        if not self.ripiego:
            return super().ripiega(contesto, degradazione)
        return self._scrivi(contesto, "ripiego")


def _passi(registro, **speciali: dict[str, Any]) -> dict[Passo, PipelineStep]:
    """S0 vera e un doppione per ogni altra fase; ``S4={...}`` ne cambia uno."""
    passi: dict[Passo, PipelineStep] = {Passo.S0: ValidazioneIngressi()}
    for passo in TUTTE[1:]:
        attributi = {"passo": passo, **speciali.get(str(passo), {})}
        passi[passo] = type(f"Doppione{passo}", (Doppione,), attributi)(registro)
    return passi


FRAGILE_S4 = {"codice": "E-S4-02", "fallimenti": 1,
              "aggiustamenti": {"E-S4-02": dimezza("run.batch_size")}}


@pytest.fixture(autouse=True)
def uscite_pulite():
    chiudi()
    yield
    chiudi()


@pytest.fixture
def scenario(tmp_path):
    return crea_scenario(tmp_path, _campioni(), con_letture=True)


@pytest.fixture
def registro() -> list:
    return []


def _esegui(config, registro, **speciali):
    run = ProjectRun(config, passi=_passi(registro, **speciali))
    return run, Esecutore(run, comando_ripresa="amplicon16s resume --config prova.yaml").esegui()


def _visti(registro, passo) -> list[int]:
    return [batch for p, batch, _ in registro if p is passo]


# --------------------------------------------------------------------------- #
# La politica                                                                  #
# --------------------------------------------------------------------------- #


def test_max_attempts_conta_i_tentativi_totali_compreso_il_primo(scenario):
    """
    **Obiettivo**: Verificare che con ``retry.max_attempts = 2`` la politica
    consenta il retry dopo il 1° tentativo fallito (``ritenta is True``) e lo
    neghi dopo il 2° tentativo fallito (``ritenta is False``).

    **Razionale Scientifico/Sistemistico**: Chiarisce senza ambiguità la
    semantica di ``max_attempts``: indica il numero totale di esecuzioni
    ammesse per la fase (1° tentativo iniziale + 1 ritentativo con parametro
    aggiustato), evitando cicli di retry superiori a quanto preventivato.
    """
    politica = PoliticaRetry.da_config(scenario.config)
    assert politica.tentativi_massimi == 2
    config = scenario.config
    adatta = dimezza("run.batch_size")
    assert politica.decidi("E-S4-02", 1, adatta, config).ritenta
    assert not politica.decidi("E-S4-02", 2, adatta, config).ritenta


def test_un_aggiustamento_su_un_parametro_non_ammesso_e_rifiutato():
    """
    **Obiettivo**: Verificare che tentare di costruire un ``Aggiustamento`` su
    un parametro biologico come ``filter.maxEE`` sollevi immediatamente ``ValueError``.

    **Razionale Scientifico/Sistemistico**: È il vincolo metodologico cardine:
    un retry automatico può modificare solo parametri di partizionamento/campionamento
    tecnico (``run.batch_size``, ``err.nbases``). Alzare automaticamente
    ``filter.maxEE`` per far passare un campione a bassa qualità falserebbe
    silenziosamente il tasso di errore delle ASV.
    """
    with pytest.raises(ValueError, match="filter.maxEE"):
        Aggiustamento("filter.maxEE", lambda v: v * 2, "soglia alzata")


def test_una_fase_non_puo_aggiustare_un_codice_non_ripetibile(scenario, registro):
    """
    **Obiettivo**: Verificare che dichiarare un aggiustamento sul codice non
    ritentabile ``E-S2-01`` in una fase faccia sollevare ``ValueError``
    all'istanziazione di ``Esecutore``.

    **Razionale Scientifico/Sistemistico**: Impedisce a chi implementa una
    classe ``PipelineStep`` di agganciare azioni di retry a codici che il
    catalogo classifica come ``REVISIONE_UMANA``.
    """
    passi = _passi(registro, S2={"aggiustamenti": {"E-S2-01": dimezza("run.batch_size")}})
    with pytest.raises(ValueError, match="E-S2-01"):
        Esecutore(ProjectRun(scenario.config, passi=passi))


# --------------------------------------------------------------------------- #
# Ritentare con l'aggiustamento                                                #
# --------------------------------------------------------------------------- #


def test_un_codice_ammesso_viene_ritentato_con_l_aggiustamento(scenario, registro):
    """
    **Obiettivo**: Verificare che quando ``S4`` fallisce al 1° tentativo con
    ``E-S4-02``, l'``Esecutore`` la rilanci dimezzando ``run.batch_size`` da
    ``24`` a ``12``, completi la pipeline e ripristini ``batch_size = 24`` per ``S5``.

    **Razionale Scientifico/Sistemistico**: Dimostra che l'aggiustamento di
    emergenza per superare un picco di memoria in ``S4`` resta strettamente
    confinato a ``S4`` senza inquinare i parametri delle fasi successive.
    """
    run, esito = _esegui(scenario.config, registro, S4=FRAGILE_S4)

    assert esito.conclusione is Conclusione.COMPLETATA
    assert _visti(registro, Passo.S4) == [24, 12]
    # L'esecuzione prosegue dopo il tentativo riuscito, con il valore dichiarato.
    assert _visti(registro, Passo.S5) == [24]
    assert run.valuta().completa


def test_l_aggiustamento_e_nel_manifesto_della_fase(scenario, registro):
    """
    **Obiettivo**: Verificare che ``manifest_s4.json`` registri l'aggiustamento
    applicato (``codice="E-S4-02"``, ``parametro="run.batch_size"``,
    ``dichiarato=24``, ``usato=12``, ``tentativo_fallito=1``), mentre il
    manifesto di ``S5`` abbia ``aggiustamenti == ()``.

    **Razionale Scientifico/Sistemistico**: Garantisce la trasparenza
    scientifica: anche se ``batch_size`` non cambia il risultato biologico, il
    manifesto della fase documenta formalmente per la tracciabilità dello studio che ``S4`` è stata completata al
    2° tentativo usando lotti da 12 campioni anziché 24.
    """
    run, _ = _esegui(scenario.config, registro, S4=FRAGILE_S4)

    manifesto = run.albero.manifesto_passo(Passo.S4, Fase.ASV_INFERENCE)
    (aggiustamento,) = manifesto.aggiustamenti
    assert aggiustamento["codice"] == "E-S4-02"
    assert aggiustamento["parametro"] == "run.batch_size"
    assert aggiustamento["dichiarato"] == 24
    assert aggiustamento["usato"] == 12
    assert aggiustamento["tentativo_fallito"] == 1
    # Le altre fasi non ne hanno.
    assert run.albero.manifesto_passo(Passo.S5, Fase.SEQTAB).aggiustamenti == ()


def test_l_aggiustamento_e_nel_log(scenario, registro):
    """
    **Obiettivo**: Verificare che l'evento di retry con ``parametro="run.batch_size"``
    e ``usato=12`` venga scritto in ``99_logs/pipeline.jsonl`` insieme al
    ``WARNING`` del codice ``E-S4-02``.

    **Razionale Scientifico/Sistemistico**: Assicura che ogni intervento
    automatico della politica di retry lasci traccia cronologica nel log JSONL.
    """
    configura(scenario.config.io.out_root)
    _esegui(scenario.config, registro, S4=FRAGILE_S4)
    righe = [
        json.loads(r)
        for r in (Path(scenario.config.io.out_root) / Fase.LOGS.value / NOME_FILE_LOG)
        .read_text().splitlines()
    ]
    ritento = [r for r in righe if r.get("parametro") == "run.batch_size"]
    assert len(ritento) == 1
    assert ritento[0]["codice"] == "E-S4-02"
    assert ritento[0]["usato"] == 12
    assert any(r.get("codice") == "E-S4-02" and r["livello"] == "WARNING" for r in righe)


def test_la_ripresa_non_rifa_la_fase_riuscita_con_l_aggiustamento(scenario, registro):
    """
    **Obiettivo**: Verificare che dopo il completamento di ``S4`` al 2°
    tentativo (con ``batch_size=12``), una successiva ripresa con la
    configurazione originale (``batch_size=24``) riconosca ``S4`` come
    ``COMPLETATA`` senza rieseguirla.

    **Razionale Scientifico/Sistemistico**: Se ``S4`` venisse invalidata perché
    il suo aggiustamento temporaneo (``12``) differisce dal valore nel file YAML
    (``24``), ogni ``resume`` successivo ricalcolerebbe ``S4`` andando
    nuovamente in OOM a ``24``! Il manifesto conserva l'impronta della
    configurazione dichiarata così che la fase resti validamente conclusa.
    """
    _esegui(scenario.config, registro, S4=FRAGILE_S4)
    registro.clear()

    run, esito = _esegui(scenario.config, registro)
    situazione = run.valuta().situazioni[Passo.S4]
    assert situazione.stato is StatoPasso.COMPLETATA
    assert esito.conclusione is Conclusione.COMPLETATA
    assert esito.eseguite == ()
    assert registro == []


def test_tre_tentativi_totali_con_due_fallimenti(scenario, registro):
    """
    **Obiettivo**: Verificare che con ``retry.max_attempts = 3`` e 2 fallimenti
    consecutivi ``run.batch_size`` venga dimezzato progressivamente ``24 -> 12 -> 6``
    registrando entrambi gli aggiustamenti ``[12, 6]`` nel manifesto.

    **Razionale Scientifico/Sistemistico**: Conferma che la funzione di
    aggiustamento viene applicata cumulativamente al valore del tentativo
    precedente quando sono permessi più ritentativi.
    """
    config = _variante(scenario.config, retry__max_attempts=3)
    fragile = {**FRAGILE_S4, "fallimenti": 2}
    run, esito = _esegui(config, registro, S4=fragile)

    assert esito.conclusione is Conclusione.COMPLETATA
    assert _visti(registro, Passo.S4) == [24, 12, 6]
    usati = [a["usato"] for a in run.albero.manifesto_passo(Passo.S4, Fase.ASV_INFERENCE).aggiustamenti]
    assert usati == [12, 6]


def test_il_modello_d_errore_si_ritenta_con_piu_basi(scenario, registro):
    r"""
    **Obiettivo**: Verificare che un fallimento di convergenza ``E-S3-01`` in
    ``S3`` attivi l'aggiustamento ``raddoppia("err.nbases")`` passando da
    ``1e8`` a ``2e8`` basi campionate.

    **Razionale Scientifico/Sistemistico**: In DADA2 ``learnErrors`` (Fase S3),
    la mancata convergenza del modello parametrico Phred si risolve fornendo un
    campione più ampio di nucleotidi (``nbases`` da $10^8$ a $2 \times 10^8$).
    """
    fragile = {"codice": "E-S3-01", "fallimenti": 1,
               "aggiustamenti": {"E-S3-01": raddoppia("err.nbases")}}
    _, esito = _esegui(scenario.config, registro, S3=fragile)
    assert esito.conclusione is Conclusione.COMPLETATA
    assert [n for p, _, n in registro if p is Passo.S3] == [1e8, 2e8]


# --------------------------------------------------------------------------- #
# Fermarsi                                                                     #
# --------------------------------------------------------------------------- #


def test_fallendo_due_volte_si_ferma_con_il_punto_di_ripresa(scenario, registro):
    """
    **Obiettivo**: Verificare che quando ``S4`` fallisce tutti i ``max_attempts=2``
    tentativi, l'``Esecutore`` si fermi con ``Conclusione.ARRESTATA``, non
    esegua ``S5`` e scriva su disco ``99_logs/punto_di_ripresa.json`` e
    ``punto_di_ripresa.txt`` con il comando ``amplicon16s resume``.

    **Razionale Scientifico/Sistemistico**: Garantisce che in caso di arresto
    notturno o su cluster HPC l'operatore trovi sul disco un referto leggibile
    (``.txt``) e uno strutturato (``.json``) con la fase esatta di arresto, i
    tentativi effettuati, l'azione consigliata dal catalogo e il comando pronto
    da copiare per riprendere l'esecuzione.
    """
    fragile = {**FRAGILE_S4, "fallimenti": 99}
    run, esito = _esegui(scenario.config, registro, S4=fragile)

    assert esito.conclusione is Conclusione.ARRESTATA
    punto = esito.punto
    assert punto.passo is Passo.S4
    assert punto.codice == "E-S4-02"
    assert punto.tentativi == 2 and punto.tentativi_massimi == 2
    assert "run.batch_size ridotto" in punto.azione  # il messaggio del catalogo
    assert punto.comando == "amplicon16s resume --config prova.yaml"
    assert _visti(registro, Passo.S4) == [24, 12]
    assert _visti(registro, Passo.S5) == []

    # Anche su disco, per dopo la chiusura del terminale.
    cartella = run.albero.cartella(Fase.LOGS)
    documento = json.loads((cartella / "punto_di_ripresa.json").read_text())
    assert documento["passo"] == "S4"
    assert documento["codice"] == "E-S4-02"
    assert documento["tentativi"] == 2
    assert documento["comando"] == punto.comando
    testo = (cartella / "punto_di_ripresa.txt").read_text()
    assert "ESECUZIONE FERMATA nella fase S4" in testo and punto.comando in testo


def test_un_esecuzione_completata_rimuove_il_vecchio_punto_di_ripresa(scenario, registro):
    """
    **Obiettivo**: Verificare che quando una corsa precedentemente arrestata
    viene ripresa e portata a termine con ``Conclusione.COMPLETATA``, il file
    ``99_logs/punto_di_ripresa.json`` venga rimosso.

    **Razionale Scientifico/Sistemistico**: Evita che nella directory ``99_logs/``
    di uno studio completato con successo rimanga un vecchio file
    ``punto_di_ripresa.json`` che farebbe sembrare la corsa ancora in stato di arresto.
    """
    _esegui(scenario.config, registro, S4={**FRAGILE_S4, "fallimenti": 99})
    cartella = Path(scenario.config.io.out_root) / Fase.LOGS.value
    assert (cartella / "punto_di_ripresa.json").exists()

    _, esito = _esegui(scenario.config, registro)
    assert esito.conclusione is Conclusione.COMPLETATA
    assert not (cartella / "punto_di_ripresa.json").exists()


def test_un_codice_a_revisione_umana_si_ferma_al_primo_tentativo(scenario, registro):
    """
    **Obiettivo**: Verificare che un errore di categoria ``REVISIONE_UMANA``
    (``E-S2-01``) arresti immediatamente la corsa al 1° tentativo
    (``tentativi == 1``) senza alcun retry.

    **Razionale Scientifico/Sistemistico**: Ribadisce che i codici di revisione
    umana (es. azzeramento delle letture di un campione dopo il filtro) non
    subiscono mai tentativi automatici.
    """
    fragile = {"codice": "E-S2-01", "fallimenti": 99}
    _, esito = _esegui(scenario.config, registro, S2=fragile)
    assert esito.punto.tentativi == 1
    assert esito.punto.categoria == "revisione_umana"
    assert "revisione umana" in esito.punto.motivo
    assert _visti(registro, Passo.S2) == [24]


def test_un_codice_ripetibile_fuori_whitelist_si_ferma_al_primo(scenario, registro):
    """
    **Obiettivo**: Verificare che se l'utente restringe ``retry.whitelist`` a
    ``["E-S2-03"]``, il codice ``E-S4-02`` (pur essendo ritentabile per catalogo)
    si fermi al 1° tentativo citando ``retry.whitelist``.

    **Razionale Scientifico/Sistemistico**: Consente al ricercatore di
    disabilitare selettivamente il retry automatico per singole fasi tramite il
    file YAML.
    """
    config = _variante(scenario.config, retry__whitelist=["E-S2-03"])
    _, esito = _esegui(config, registro, S4={**FRAGILE_S4, "fallimenti": 99})
    assert esito.punto.tentativi == 1
    assert "retry.whitelist" in esito.punto.motivo
    assert _visti(registro, Passo.S4) == [24]


def test_con_retry_disattivato_nessun_codice_viene_ritentato(scenario, registro):
    """
    **Obiettivo**: Verificare che con ``retry.enabled = False`` qualsiasi errore
    ritentabile (``E-S3-01``) fermi la pipeline al 1° tentativo citando ``retry.enabled``.

    **Razionale Scientifico/Sistemistico**: Rispetta l'interruttore globale
    ``retry.enabled`` quando si desidera un'esecuzione strettamente a singolo tentativo.
    """
    config = _variante(scenario.config, retry__enabled=False)
    fragile_s3 = {"codice": "E-S3-01", "fallimenti": 1,
                  "aggiustamenti": {"E-S3-01": raddoppia("err.nbases")}}
    _, esito = _esegui(config, registro, S3=fragile_s3, S4=FRAGILE_S4)
    assert esito.conclusione is Conclusione.ARRESTATA
    assert esito.punto.passo is Passo.S3 and esito.punto.tentativi == 1
    assert "retry.enabled" in esito.punto.motivo


def test_senza_azione_correttiva_non_si_ritenta_identico(scenario, registro):
    """
    **Obiettivo**: Verificare che se una fase solleva ``E-S4-02`` ma non ha
    dichiarato alcun ``Aggiustamento`` per quel codice, l'esecutore si fermi al
    1° tentativo citando ``azione correttiva``.

    **Razionale Scientifico/Sistemistico**: Impedisce i "retry ciechi": rilanciare
    lo stesso identico calcolo deterministico con gli stessi identici parametri
    produrrebbe lo stesso fallimento sprecando risorse.
    """
    fragile = {"codice": "E-S4-02", "fallimenti": 1}
    _, esito = _esegui(scenario.config, registro, S4=fragile)
    assert esito.punto.tentativi == 1
    assert "azione correttiva" in esito.punto.motivo


def test_un_azione_correttiva_al_limite_non_si_ripete(scenario, registro):
    """
    **Obiettivo**: Verificare che se ``run.batch_size`` è già pari a ``1`` (il
    minimo fisico indivisibile), ``dimezza("run.batch_size")`` non possa
    ridurlo ulteriormente e l'esecutore si fermi al 1° tentativo citando ``limite``.

    **Razionale Scientifico/Sistemistico**: Evita che il dimezzamento porti
    ``batch_size`` a ``0`` o ripeta inutilmente il tentativo con ``batch_size = 1``.
    """
    config = _variante(scenario.config, run__batch_size=1)
    _, esito = _esegui(config, registro, S4=FRAGILE_S4)
    assert esito.punto.tentativi == 1
    assert "limite" in esito.punto.motivo


def test_il_messaggio_distingue_retry_e_retry_poi_revisione(scenario, registro):
    """
    **Obiettivo**: Verificare che all'esaurimento dei tentativi il motivo in
    ``punto_di_ripresa`` indichi ``serve la revisione umana`` per ``E-S3-01``
    (``retry_poi_revisione``) e ``ammette il retry`` per ``E-S4-02`` (``retry_automatico``).

    **Razionale Scientifico/Sistemistico**: Comunica fedelmente all'operatore la
    differenza tra un mancato raggiungimento della convergenza statistica in S3
    (che dopo il retry richiede l'ispezione dei profili Phred da parte del
    biologo) e un limite puramente hardware di RAM in S4.
    """
    s3 = {"codice": "E-S3-01", "fallimenti": 99,
          "aggiustamenti": {"E-S3-01": raddoppia("err.nbases")}}
    _, dopo_s3 = _esegui(scenario.config, registro, S3=s3)
    assert dopo_s3.punto.categoria == "retry_poi_revisione"
    assert "serve la revisione umana" in dopo_s3.punto.motivo

    _, dopo_s4 = _esegui(
        _variante(scenario.config, io__out_root=Path(scenario.radice) / "altra"),
        registro, S4={**FRAGILE_S4, "fallimenti": 99},
    )
    assert dopo_s4.punto.categoria == "retry_automatico"
    assert "serve la revisione umana" not in dopo_s4.punto.motivo
    assert "ammette il retry" in dopo_s4.punto.motivo


# --------------------------------------------------------------------------- #
# Whitelist incoerente con il catalogo                                         #
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize(
    "whitelist",
    [["E-S4-02", "E-S2-01"], ["E-S99-01"]],
    ids=["codice-non-ripetibile", "codice-inesistente"],
)
def test_una_whitelist_incoerente_e_respinta_prima_di_ogni_fase(
    scenario, registro, whitelist
):
    """
    **Obiettivo**: Verificare che inserire in ``retry.whitelist`` un codice non
    ritentabile (``E-S2-01``) o inesistente (``E-S99-01``) venga bloccato dal
    Gate G15 con codice ``E-G15-09`` prima di creare qualsiasi cartella di fase.

    **Razionale Scientifico/Sistemistico**: Impedisce che una configurazione
    YAML manomessa allarghi la whitelist di retry oltre i 4 codici autorizzati
    dal catalogo.
    """
    config = _variante(scenario.config, retry__whitelist=whitelist)
    with pytest.raises(ErroreGate) as info:
        _esegui(config, registro)
    assert info.value.gate == "G15"
    assert "E-G15-09" in {v.codice for v in info.value.violazioni}
    assert registro == []
    assert not (Path(config.io.out_root) / Fase.INPUT_VALIDATION.value).exists()


def test_una_whitelist_ristretta_e_ammessa(scenario, registro):
    """
    **Obiettivo**: Verificare che un sottoinsieme proprio dei 4 codici ammessi
    (es. ``["E-S4-02"]``) superi G15 e completi regolarmente l'esecuzione.

    **Razionale Scientifico/Sistemistico**: L'utente non può aggiungere codici
    alla whitelist, ma ha sempre facoltà di restringerla.
    """
    config = _variante(scenario.config, retry__whitelist=["E-S4-02"])
    _, esito = _esegui(config, registro)
    assert esito.conclusione is Conclusione.COMPLETATA


# --------------------------------------------------------------------------- #
# Degradazioni                                                                 #
# --------------------------------------------------------------------------- #


def test_una_degradazione_non_interrompe_ed_e_registrata(scenario, registro):
    """
    **Obiettivo**: Verificare che chiamare ``contesto.degrada("E-S6-02", ...)``
    durante ``S6`` non interrompa la corsa e salvi l'evento nella lista
    ``degradazioni`` di ``manifest_s6.json``.

    **Razionale Scientifico/Sistemistico**: Quando la frazione di sequenze
    chimeriche rimosse in ``S6`` supera la soglia di attenzione ma resta sotto
    la soglia bloccante (``E-S6-02``), la pipeline prosegue annotando l'anomalia
    nel manifesto per il report finale.
    """
    run, esito = _esegui(scenario.config, registro, S6={"degradazione": "degrada"})

    assert esito.conclusione is Conclusione.COMPLETATA
    (degradazione,) = run.albero.manifesto_passo(Passo.S6, Fase.CHIMERA).degradazioni
    assert degradazione["codice"] == "E-S6-02"
    assert degradazione["categoria"] == "degradazione_automatica"
    assert degradazione["dettaglio"] == "frazione chimerica 0,31"


def test_una_degradazione_sollevata_con_ripiego_non_interrompe(scenario, registro):
    """
    **Obiettivo**: Verificare che se ``S11`` solleva ``DegradazioneRichiesta("E-S11-02")``
    e implementa il metodo ``ripiega``, il passo completi usando il ripiego e
    registri ``E-S11-02`` nel manifesto.

    **Razionale Scientifico/Sistemistico**: Se in ``S11`` la curva KatharoSeq sui
    controlli positivi a diluizione seriale non converge o presenta una
    composizione deviante (``E-S11-02``), ``PipelineStep`` invoca ``ripiega()``
    per applicare la soglia di profondità di riserva senza abortire la pipeline.
    """
    s11 = {"degradazione": "solleva", "codice_degradazione": "E-S11-02", "ripiego": True}
    run, esito = _esegui(scenario.config, registro, S11=s11)

    assert esito.conclusione is Conclusione.COMPLETATA
    manifesto = run.albero.manifesto_passo(Passo.S11, Fase.CONTROLS)
    assert [d["codice"] for d in manifesto.degradazioni] == ["E-S11-02"]
    assert "ripiego" in (run.albero.cartella(Fase.CONTROLS) / "s11.txt").read_text()


def test_una_degradazione_sollevata_senza_ripiego_e_un_difetto(scenario, registro):
    """
    **Obiettivo**: Verificare che sollevare ``DegradazioneRichiesta`` in una
    fase che non ha ridefinito ``ripiega()`` sollevi ``RuntimeError("senza dichiarare un ripiego")``.

    **Razionale Scientifico/Sistemistico**: Impedisce che una fase sollevi
    un'eccezione di degradazione senza aver effettivamente implementato la
    logica alternativa di calcolo.
    """
    with pytest.raises(RuntimeError, match="senza dichiarare un ripiego"):
        _esegui(scenario.config, registro, S11={"degradazione": "solleva",
                                                "codice_degradazione": "E-S11-02"})


def test_degradare_con_un_codice_che_ferma_e_rifiutato(scenario):
    """
    **Obiettivo**: Verificare che passare a ``contesto.degrada(...)`` un codice
    di ``REVISIONE_UMANA`` (``E-S2-01``) sollevi ``ValueError``.

    **Razionale Scientifico/Sistemistico**: Impedisce che un errore bloccante
    venga declassato a semplice avviso di degradazione tramite ``contesto.degrada``.
    """
    run = ProjectRun(scenario.config)
    contesto = run.contesto(Passo.S0)
    with pytest.raises(ValueError, match="revisione_umana"):
        contesto.degrada("E-S2-01", "non e' una degradazione")


def test_le_degradazioni_di_s0_finiscono_nel_suo_manifesto(tmp_path, registro):
    """
    **Obiettivo**: Verificare che gli avvisi di degradazione emessi dai gate di
    ``S0`` (``E-S0-15`` da G08 per piastra con pochi blank ed ``E-S1-01`` da G09
    per scarto di troncamento) vengano registrati in ``manifest_s0.json``.

    **Razionale Scientifico/Sistemistico**: Collega gli avvisi pre-analitici dei
    gate in S0 al sistema dei manifesti di fase, così che ``S12`` e il report
    finale ``S14`` possano consultarli in modo uniforme.
    """
    scenario = crea_scenario(
        tmp_path, _campioni(), con_letture=True, con_arricchimento=True,
        sovrascrivi={"decontam": {"min_blanks": 3}},
    )
    run, esito = _esegui(scenario.config, registro)
    assert esito.conclusione is Conclusione.COMPLETATA
    codici = [d["codice"] for d in run.albero.manifesto_passo(Passo.S0, Fase.INPUT_VALIDATION).degradazioni]
    assert sorted(codici) == ["E-S0-15", "E-S1-01"]


# --------------------------------------------------------------------------- #
# Controlli di risorse a ogni avvio                                            #
# --------------------------------------------------------------------------- #


def test_i_controlli_di_risorse_girano_anche_con_s0_conclusa(scenario, registro):
    """
    **Obiettivo**: Verificare che se ``S0`` è già ``COMPLETATA`` su disco ma si
    riprende l'esecuzione con ``run.threads = 100000``, l'``Esecutore`` arresti
    subito la corsa ai ``controlli di avvio`` con ``E-S0-14`` senza invalidare
    né rieseguire ``S0``.

    **Razionale Scientifico/Sistemistico**: Poiché ``run.threads`` è escluso
    dall'impronta dei risultati, ``S0`` resta giustamente ``COMPLETATA`` quando
    si cambia macchina o numero di core; tuttavia l'``Esecutore`` deve sempre
    rieseguire il Gate G14 all'avvio del processo per verificare che i nuovi
    ``threads`` richiesti esistano davvero sulla macchina corrente.
    """
    _esegui(scenario.config, registro)
    registro.clear()

    # run.threads non incide sui risultati: S0 resta conclusa.
    troppi = _variante(scenario.config, run__threads=100000)
    run = ProjectRun(troppi, passi=_passi(registro))
    assert run.valuta().situazioni[Passo.S0].stato is StatoPasso.COMPLETATA
    impronta_s0 = run.valuta().situazioni[Passo.S0].impronta

    esito = Esecutore(run).esegui()
    assert esito.conclusione is Conclusione.ARRESTATA
    assert esito.punto.codice == "E-S0-14"
    assert esito.punto.origine == "controlli di avvio"
    assert "100000" in esito.punto.dettaglio
    # Nessuna fase eseguita, S0 compresa.
    assert registro == [] and esito.eseguite == ()
    assert run.valuta().situazioni[Passo.S0].impronta == impronta_s0


# --------------------------------------------------------------------------- #
# Riga di comando                                                              #
# --------------------------------------------------------------------------- #


@pytest.fixture
def file_config(scenario, tmp_path) -> Path:
    percorso = tmp_path / "config.yaml"
    percorso.write_text(
        yaml.safe_dump(scenario.config.model_dump(mode="json"), sort_keys=False),
        encoding="utf-8",
    )
    return percorso


@pytest.fixture
def con_doppioni(monkeypatch, registro):
    """Sostituisce le fasi realizzate con S0 vera piu' i doppioni."""
    speciali: dict[str, Any] = {}
    monkeypatch.setattr(cli, "_passi", lambda: _passi(registro, **speciali))
    return speciali


def _cli(*argomenti, capsys) -> tuple[int, str]:
    codice = cli.main(list(argomenti))
    return codice, capsys.readouterr().out


def test_cli_run_completa_con_uscita_zero(file_config, con_doppioni, capsys):
    """
    **Obiettivo**: Verificare che ``amplicon16s run --config config.yaml`` su
    una directory di output nuova completi tutte le fasi restituendo exit code ``0``.

    **Razionale Scientifico/Sistemistico**: Certifica il funzionamento nominale
    del comando primario di esecuzione da terminale.
    """
    codice, uscita = _cli("run", "--config", str(file_config), capsys=capsys)
    assert codice == cli.USCITA_SUCCESSO == 0
    assert "Esecuzione completata" in uscita


def test_cli_run_non_sovrascrive_un_esecuzione_esistente(
    file_config, scenario, con_doppioni, capsys
):
    """
    **Obiettivo**: Verificare che lanciare una seconda volta ``amplicon16s run``
    sulla stessa ``out_root`` fallisca con codice ``3`` senza modificare alcun
    file su disco e suggerendo ``amplicon16s resume``.

    **Razionale Scientifico/Sistemistico**: Previene la distruzione accidentale
    di un'elaborazione già conclusa o parziale se l'operatore invoca per errore
    ``run`` anziché ``resume``.
    """
    assert _cli("run", "--config", str(file_config), capsys=capsys)[0] == 0
    radice = Path(scenario.config.io.out_root)
    prima = {p: p.stat().st_mtime_ns for p in radice.rglob("*") if p.is_file()}

    codice, uscita = _cli("run", "--config", str(file_config), capsys=capsys)
    assert codice == cli.USCITA_CONFIGURAZIONE == 3
    assert "non sovrascrive" in uscita
    assert f"amplicon16s resume --config {file_config.resolve()}" in uscita
    dopo = {p: p.stat().st_mtime_ns for p in radice.rglob("*") if p.is_file()}
    assert dopo == prima  # nessun file toccato, nemmeno i log


def test_cli_run_rifiuta_anche_una_cartella_con_altri_file(file_config, scenario, capsys):
    """
    **Obiettivo**: Verificare che ``amplicon16s run`` rifiuti con codice ``3``
    anche una cartella ``out_root`` contenente un singolo file estraneo (``appunti.txt``).

    **Razionale Scientifico/Sistemistico**: Garantisce che ``run`` parta sempre
    da una directory immacolata, preservando eventuali file dell'utente.
    """
    radice = Path(scenario.config.io.out_root)
    radice.mkdir(parents=True)
    (radice / "appunti.txt").write_text("miei")
    assert _cli("run", "--config", str(file_config), capsys=capsys)[0] == 3
    assert (radice / "appunti.txt").read_text() == "miei"


def test_cli_arresto_e_ripresa(file_config, con_doppioni, registro, capsys):
    """
    **Obiettivo**: Verificare il flusso operativo completo: ``run`` si arresta
    in ``S4`` con exit code ``4`` stampando il comando di ripresa; rimossa la
    causa del guasto, ``resume`` riparte esattamente da ``S4`` (senza rifare
    ``S0..S3``) ed esce con codice ``0``.

    **Razionale Scientifico/Sistemistico**: Collauda lo scenario reale di
    recupero da incidente durante l'elaborazione dei 960 campioni di OSD-734.
    """
    con_doppioni["S4"] = {**FRAGILE_S4, "fallimenti": 99}
    codice, uscita = _cli("run", "--config", str(file_config), capsys=capsys)
    assert codice == cli.USCITA_ARRESTO == 4
    assert "ESECUZIONE FERMATA nella fase S4" in uscita
    assert f"amplicon16s resume --config {file_config.resolve()}" in uscita

    # Rimossa la causa, resume riparte da S4 senza rifare S0-S3.
    con_doppioni.clear()
    registro.clear()
    codice, _ = _cli("resume", "--config", str(file_config), capsys=capsys)
    assert codice == 0
    assert [p for p, _, _ in registro][0] is Passo.S4


def test_cli_validate_esegue_solo_s0(file_config, scenario, con_doppioni, registro, capsys):
    """
    **Obiettivo**: Verificare che ``amplicon16s validate`` esegua unicamente la
    Fase ``S0`` scrivendo ``manifest_S0.json`` (senza avviare ``S1..S14``) e che
    una seconda chiamata rilevi ``S0 gia' conclusa``.

    **Razionale Scientifico/Sistemistico**: Permette al biologo di verificare
    rapidamente (in ~24 secondi su OSD-734) tutti i 15 gate di ingresso prima di
    sottomettere il job completo di calcolo R/DADA2.
    """
    codice, uscita = _cli("validate", "--config", str(file_config), capsys=capsys)
    assert codice == 0 and "S0 conclusa" in uscita
    assert registro == []  # nessun doppione eseguito
    radice = Path(scenario.config.io.out_root)
    assert (radice / Fase.INPUT_VALIDATION.value / "manifest_S0.json").exists()

    codice, uscita = _cli("validate", "--config", str(file_config), capsys=capsys)
    assert codice == 0 and "gia' conclusa" in uscita


def test_cli_report_e_un_resoconto_provvisorio(file_config, con_doppioni, capsys):
    """
    **Obiettivo**: Verificare che ``amplicon16s report`` stampi lo stato di
    tutte le fasi (inclusa ``S9 disattivata``), gli aggiustamenti applicati
    (``[E-S4-02]: run.batch_size 24 -> 12``) e le degradazioni (``[E-S6-02]``).

    **Razionale Scientifico/Sistemistico**: Fornisce un cruscotto immediato da
    riga di comando per ispezionare lo stato di avanzamento e le decisioni
    automatiche prese dalla pipeline.
    """
    con_doppioni["S4"] = FRAGILE_S4
    con_doppioni["S6"] = {"degradazione": "degrada"}
    assert _cli("run", "--config", str(file_config), capsys=capsys)[0] == 0

    codice, uscita = _cli("report", "--config", str(file_config), capsys=capsys)
    assert codice == 0
    assert "RESOCONTO PROVVISORIO" in uscita
    assert "S9   disattivata" in uscita
    assert "aggiustamento [E-S4-02]: run.batch_size 24 -> 12" in uscita
    assert "degradazione [E-S6-02]" in uscita
    assert "Esecuzione completa." in uscita


def test_cli_report_senza_esecuzione(file_config, capsys):
    """
    **Obiettivo**: Verificare che invocare ``amplicon16s report`` su una
    cartella ``out_root`` non ancora esistente restituisca exit code ``3``.

    **Razionale Scientifico/Sistemistico**: Segnala chiaramente che non esiste
    ancora alcuna esecuzione da rendicontare.
    """
    assert _cli("report", "--config", str(file_config), capsys=capsys)[0] == 3


def test_cli_oggi_si_ferma_alla_prima_fase_non_realizzata(file_config, capsys):
    """
    **Obiettivo**: Verificare che senza i doppioni di test ``amplicon16s run``
    completi ``S0`` e si fermi ordinatamente su ``S1`` con codice
    ``USCITA_FASE_NON_REALIZZATA == 5``.

    **Razionale Scientifico/Sistemistico**: Verifica il comportamento onesto del
    codice di produzione al termine della Settimana 9/10 prima dell'innesto
    degli step R reali delle fasi successive.
    """
    codice, uscita = _cli("run", "--config", str(file_config), capsys=capsys)
    assert codice == cli.USCITA_FASE_NON_REALIZZATA == 5
    assert "S1, non e' ancora realizzata" in uscita


def test_cli_errore_di_configurazione(tmp_path, capsys):
    """
    **Obiettivo**: Verificare che un file YAML sintatticamente corrotto o
    inesistente faccia uscire ``run`` e ``resume`` con codice ``3``.

    **Razionale Scientifico/Sistemistico**: Mappa tutti gli errori di lettura e
    validazione dello schema YAML sul codice di uscita dedicato alla
    configurazione (``3``).
    """
    rotta = tmp_path / "rotta.yaml"
    rotta.write_text("filter: [\n")
    assert _cli("run", "--config", str(rotta), capsys=capsys)[0] == 3
    assert _cli("resume", "--config", str(tmp_path / "assente.yaml"), capsys=capsys)[0] == 3


def test_cli_whitelist_incoerente_e_errore_di_configurazione(
    scenario, tmp_path, con_doppioni, registro, capsys
):
    """
    **Obiettivo**: Verificare che una violazione del Gate G15 (es. ``E-G15-09``
    sulla whitelist) invocata da CLI restituisca codice ``3`` stampando ``E-G15-09``.

    **Razionale Scientifico/Sistemistico**: Garantisce che gli errori statici di
    configurazione (G15) siano distinti tramite exit code ``3`` dagli arresti
    operativi di esecuzione (exit code ``4``).
    """
    dati = scenario.config.model_dump(mode="json")
    dati["retry"]["whitelist"] = ["E-S4-02", "E-S14-01"]
    percorso = tmp_path / "incoerente.yaml"
    percorso.write_text(yaml.safe_dump(dati))
    codice, uscita = _cli("resume", "--config", str(percorso), capsys=capsys)
    assert codice == 3 and "E-G15-09" in uscita
    assert registro == []


def test_cli_controlli_di_risorse_nella_ripresa(scenario, tmp_path, con_doppioni, registro, capsys):
    """
    **Obiettivo**: Verificare che un fallimento di G14 (``run.threads = 100000``)
    durante ``amplicon16s resume`` esca con codice di arresto ``4`` citando
    ``E-S0-14`` e ``controlli di avvio``.

    **Razionale Scientifico/Sistemistico**: Accerta che i controlli sulle
    risorse hardware della macchina vengano applicati anche quando si riprende
    da CLI una corsa già iniziata.
    """
    percorso = tmp_path / "c.yaml"
    dati = scenario.config.model_dump(mode="json")
    percorso.write_text(yaml.safe_dump(dati))
    assert _cli("run", "--config", str(percorso), capsys=capsys)[0] == 0

    dati["run"]["threads"] = 100000
    percorso.write_text(yaml.safe_dump(dati))
    registro.clear()
    codice, uscita = _cli("resume", "--config", str(percorso), capsys=capsys)
    assert codice == 4 and "E-S0-14" in uscita and "controlli di avvio" in uscita
    assert registro == []


def test_il_file_di_configurazione_non_viene_mai_modificato(
    file_config, con_doppioni, capsys
):
    """
    **Obiettivo**: Verificare che l'esecuzione di ``validate``, ``run``,
    ``resume`` e ``report`` (anche in presenza di aggiustamenti automatici di
    retry) lasci inalterati byte per byte e per ``mtime`` il file ``config.yaml`` dell'utente.

    **Razionale Scientifico/Sistemistico**: Il file YAML sorgente fornito dal
    ricercatore è strettamente *read-only*; tutte le risoluzioni e gli
    aggiustamenti vengono scritti esclusivamente dentro ``out_root/00_config/``
    e nei manifesti di fase.
    """
    contenuto = file_config.read_bytes()
    istante = file_config.stat().st_mtime_ns
    con_doppioni["S4"] = FRAGILE_S4
    for comando in ("validate", "run", "resume", "report"):
        _cli(comando, "--config", str(file_config), capsys=capsys)
    assert file_config.read_bytes() == contenuto
    assert file_config.stat().st_mtime_ns == istante


def test_l_aggiustamento_non_modifica_la_configurazione_in_memoria(scenario, registro):
    """
    **Obiettivo**: Verificare che al termine di un retry su ``S4`` gli oggetti
    ``run.config.run.batch_size`` e ``scenario.config.run.batch_size`` in
    memoria conservino il valore originale ``24``.

    **Razionale Scientifico/Sistemistico**: Garantisce l'immutabilità della
    configurazione di progetto in memoria durante l'esecuzione.
    """
    run, _ = _esegui(scenario.config, registro, S4=FRAGILE_S4)
    assert run.config.run.batch_size == 24
    assert scenario.config.run.batch_size == 24


def test_un_codice_ripetibile_senza_aggiustamento_e_un_arresto_non_un_errore(
    scenario, registro
):
    """
    **Obiettivo**: Verificare che se ``S5`` solleva ``E-S5-01`` senza avere un
    aggiustamento configurato, ``Esecutore.esegui()`` restituisca ordinatamente
    ``Conclusione.ARRESTATA`` su ``Passo.S5`` senza far propagare un'eccezione non gestita.

    **Razionale Scientifico/Sistemistico**: Assicura che qualsiasi fallimento di
    fase si traduca sempre in un arresto governato con scrittura del punto di
    ripresa.
    """
    _, esito = _esegui(scenario.config, registro,
                       S5={"codice": "E-S5-01", "fallimenti": 99})
    assert esito.conclusione is Conclusione.ARRESTATA
    assert esito.punto.passo is Passo.S5


def test_le_fasi_eseguite_sono_quelle_attive(scenario, registro):
    """
    **Obiettivo**: Verificare che ``esito.eseguite`` contenga esattamente la
    sequenza restituita da ``GRAFO.attive(scenario.config)``.

    **Razionale Scientifico/Sistemistico**: Conferma l'allineamento tra la
    pianificazione statica del grafo e la lista dei risultati prodotta dall'esecutore.
    """
    _, esito = _esegui(scenario.config, registro)
    assert [r.passo for r in esito.eseguite] == list(GRAFO.attive(scenario.config))


# --------------------------------------------------------------------------- #
# La configurazione registrata in 00_config                                   #
# --------------------------------------------------------------------------- #


def _scrivi_config(percorso: Path, config: Config, **valori: Any) -> Path:
    dati = config.model_dump(mode="json")
    for nome, v in valori.items():
        gruppo, chiave = nome.split("__")
        dati[gruppo][chiave] = v
    percorso.write_text(yaml.safe_dump(dati, sort_keys=False), encoding="utf-8")
    return percorso


def _registrate(config: Config) -> list[dict]:
    return [
        yaml.safe_load(p.read_text()) | {"_nome": p.name}
        for p in versioni_registrate(config.io.out_root)
    ]


def test_run_registra_la_configurazione_con_il_digest(
    file_config, scenario, con_doppioni, capsys
):
    """
    **Obiettivo**: Verificare che ``amplicon16s run`` salvi in ``00_config/resolved.yaml``
    la configurazione risolta con ``digest`` coincidente con ``risolvi(config).digest``.

    **Razionale Scientifico/Sistemistico**: Sigilla all'avvio della corsa la
    fotografia esatta di tutti i parametri espliciti, predefiniti e derivati.
    """
    assert _cli("run", "--config", str(file_config), capsys=capsys)[0] == 0
    (registrata,) = _registrate(scenario.config)
    assert registrata["_nome"] == "resolved.yaml"
    assert registrata["digest"] == risolvi(scenario.config).digest


def test_resume_con_lo_stesso_digest_non_scrive_nulla(
    file_config, scenario, con_doppioni, capsys
):
    """
    **Obiettivo**: Verificare che eseguire ``amplicon16s resume`` con la
    medesima configurazione non crei nuovi file in ``00_config/`` né modifichi i
    timestamp dei file esistenti.

    **Razionale Scientifico/Sistemistico**: Evita la proliferazione di file
    ``resolved_<N>.yaml`` identici quando si riprende una corsa senza aver
    toccato il file YAML.
    """
    _cli("run", "--config", str(file_config), capsys=capsys)
    cartella = Path(scenario.config.io.out_root) / Fase.CONFIG.value
    prima = {p.name: p.stat().st_mtime_ns for p in cartella.iterdir()}

    assert _cli("resume", "--config", str(file_config), capsys=capsys)[0] == 0
    assert {p.name: p.stat().st_mtime_ns for p in cartella.iterdir()} == prima


def test_resume_con_un_digest_diverso_registra_accanto_senza_sovrascrivere(
    file_config, scenario, con_doppioni, registro, tmp_path, capsys
):
    """
    **Obiettivo**: Verificare che riprendendo con ``run.threads=2`` (e poi con
    ``retry.max_attempts=3`` e infine tornando alla configurazione iniziale) non
    venga ricalcolata alcuna fase, ``resolved.yaml`` resti intatto e vengano
    affiancati ``resolved_2.yaml``, ``resolved_3.yaml`` e ``resolved_4.yaml``
    con i campi ``precedente`` e ``differenze_dalla_precedente``.

    **Razionale Scientifico/Sistemistico**: Coniuga il risparmio computazionale
    (non ricalcolare le fasi se cambiano solo ``threads`` o ``retry``) con la
    tracciabilità forense: la catena ``resolved.yaml -> resolved_2.yaml -> ...``
    documenta ogni singola variazione operativa avvenuta tra una ripresa e l'altra.
    """
    _cli("run", "--config", str(file_config), capsys=capsys)
    originale = versioni_registrate(scenario.config.io.out_root)[0].read_bytes()

    # run.threads e' fuori dall'impronta dei risultati ma non dal digest.
    registro.clear()
    diversa = _scrivi_config(tmp_path / "c2.yaml", scenario.config, run__threads=2)
    codice, uscita = _cli("resume", "--config", str(diversa), capsys=capsys)
    assert codice == 0 and registro == []  # nessuna fase rifatta
    assert "registrata in resolved_2.yaml" in uscita

    prima, seconda = _registrate(scenario.config)
    assert versioni_registrate(scenario.config.io.out_root)[0].read_bytes() == originale
    assert seconda["_nome"] == "resolved_2.yaml"
    assert seconda["precedente"] == "resolved.yaml"
    assert seconda["differenze_dalla_precedente"] == ["run.threads: 1 -> 2"]
    assert seconda["digest"] != prima["digest"]

    # Un parametro di retry e' un'altra configurazione ancora.
    terza = _scrivi_config(tmp_path / "c3.yaml", scenario.config, run__threads=2,
                           retry__max_attempts=3)
    _cli("resume", "--config", str(terza), capsys=capsys)
    assert _registrate(scenario.config)[-1]["differenze_dalla_precedente"] == [
        "retry.max_attempts: 2 -> 3"
    ]

    # Tornare alla prima e' a sua volta un cambiamento, e resta ricostruibile.
    _cli("resume", "--config", str(file_config), capsys=capsys)
    versioni = _registrate(scenario.config)
    assert [v["_nome"] for v in versioni] == [
        "resolved.yaml", "resolved_2.yaml", "resolved_3.yaml", "resolved_4.yaml"
    ]
    assert versioni[-1]["digest"] == versioni[0]["digest"]


def test_validate_su_cartella_nuova_segue_la_regola_di_run(
    file_config, scenario, con_doppioni, capsys
):
    """
    **Obiettivo**: Verificare che ``amplicon16s validate`` su una directory
    nuova scriva ``00_config/resolved.yaml`` esattamente come ``run``.

    **Razionale Scientifico/Sistemistico**: Garantisce che anche una corsa
    iniziata con ``validate`` inizializzi correttamente ``00_config/resolved.yaml``.
    """
    assert _cli("validate", "--config", str(file_config), capsys=capsys)[0] == 0
    assert [v["_nome"] for v in _registrate(scenario.config)] == ["resolved.yaml"]


def test_validate_con_s0_conclusa_segue_la_regola_di_resume(
    file_config, scenario, con_doppioni, tmp_path, capsys
):
    """
    **Obiettivo**: Verificare che rieseguire ``validate`` con ``S0`` già
    conclusa non duplichi ``resolved.yaml`` a parità di digest, ma crei
    ``resolved_2.yaml`` se viene modificato ``retry.enabled``.

    **Razionale Scientifico/Sistemistico**: Mantiene coerente il versionamento
    di ``00_config/`` tra i comandi ``validate`` e ``resume``.
    """
    _cli("validate", "--config", str(file_config), capsys=capsys)
    _cli("validate", "--config", str(file_config), capsys=capsys)
    assert len(_registrate(scenario.config)) == 1

    diversa = _scrivi_config(tmp_path / "c2.yaml", scenario.config, retry__enabled=False)
    codice, uscita = _cli("validate", "--config", str(diversa), capsys=capsys)
    assert codice == 0 and "gia' conclusa" in uscita
    assert [v["_nome"] for v in _registrate(scenario.config)] == [
        "resolved.yaml", "resolved_2.yaml"
    ]


def test_un_arresto_ai_controlli_di_avvio_non_registra_nulla(
    file_config, scenario, con_doppioni, tmp_path, capsys
):
    """
    **Obiettivo**: Verificare che se ``amplicon16s resume`` viene respinto ai
    controlli di avvio (``run.threads = 100000``, codice ``4``), non venga
    registrato alcun ``resolved_2.yaml`` in ``00_config/``.

    **Razionale Scientifico/Sistemistico**: Impedisce che una configurazione
    hardware impossibile (rifiutata da G14 prima di eseguire alcunché) entri
    nella storia delle configurazioni effettivamente utilizzate.
    """
    _cli("run", "--config", str(file_config), capsys=capsys)
    troppi = _scrivi_config(tmp_path / "c2.yaml", scenario.config, run__threads=100000)
    assert _cli("resume", "--config", str(troppi), capsys=capsys)[0] == 4
    assert len(_registrate(scenario.config)) == 1


def test_il_log_dice_con_quale_configurazione_si_esegue(
    file_config, scenario, con_doppioni, capsys
):
    """
    **Obiettivo**: Verificare che ad ogni invocazione di ``run`` e ``resume`` il
    log JSONL registri l'evento ``"configurazione in uso"`` specificando il nome
    del file ``resolved*.yaml`` e il flag booleano ``registrata_ora``.

    **Razionale Scientifico/Sistemistico**: Lega ogni sessione di log in
    ``99_logs/pipeline.jsonl`` all'esatta versione di ``00_config/resolved*.yaml``
    in vigore in quel momento.
    """
    _cli("run", "--config", str(file_config), capsys=capsys)
    _cli("resume", "--config", str(file_config), capsys=capsys)
    righe = [
        json.loads(r)
        for r in (Path(scenario.config.io.out_root) / Fase.LOGS.value / NOME_FILE_LOG)
        .read_text().splitlines()
    ]
    in_uso = [r for r in righe if r["messaggio"] == "configurazione in uso"]
    assert [(r["file"], r["registrata_ora"]) for r in in_uso] == [
        ("resolved.yaml", True), ("resolved.yaml", False)
    ]


def test_le_versioni_si_ordinano_per_numero(tmp_path):
    """
    **Obiettivo**: Verificare che ``versioni_registrate`` restituisca i file in
    ordine numerico intero (``resolved.yaml``, ``resolved_2.yaml``,
    ``resolved_9.yaml``, ``resolved_10.yaml``) ignorando ``manifest.json`` e nomi non numerici.

    **Razionale Scientifico/Sistemistico**: Previene il classico bug di
    ordinamento lessicografico ASCII (in cui ``resolved_10.yaml`` precederebbe
    ``resolved_2.yaml``), garantendo che ``versioni_registrate[-1]`` punti
    sempre all'ultima configurazione registrata anche oltre la 9ª ripresa.
    """
    cartella = tmp_path / Fase.CONFIG.value
    cartella.mkdir()
    for nome in ("resolved_10.yaml", "resolved.yaml", "resolved_2.yaml", "resolved_9.yaml",
                 "manifest.json", "resolved_x.yaml"):
        (cartella / nome).write_text("digest: x\n")
    assert [p.name for p in versioni_registrate(tmp_path)] == [
        "resolved.yaml", "resolved_2.yaml", "resolved_9.yaml", "resolved_10.yaml"
    ]
