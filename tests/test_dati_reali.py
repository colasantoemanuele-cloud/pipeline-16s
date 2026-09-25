"""Suite di validazione end-to-end sul dataset reale di riferimento NASA GeneLab OSD-734.

Inquadramento nel Piano Operativo:
    - **Settimane di riferimento**: **Settimana 6 e Settimana 7 (W6/W7 — Fase F2:
      Crosswalk, inventario dei campioni, 15 gate di ingresso e Fase S0 completa)**.
    - **Scopo del modulo**: Esegue l'intera validazione dei metadati e la Fase S0
      completa sui **960 file ``.fastq.gz`` reali (2,4 GB)** e sulle tabelle
      ISA-Tab originali dello studio **OSD-734** (Microbial Monitoring of the
      International Space Station). Verifica sul dato vero:
        * La corrispondenza biunivoca 1:1 tra i 960 file FASTQ e le 960 righe
          dell'Assay Table 16S;
        * Il join ristretto che esclude i campioni di altri saggi (come
          ``"solvent control"`` appartenente alla metagenomica) presenti nella
          Study Sample Table;
        * La ripartizione esatta nelle tre classi (**803 biologici, 80 controlli
          positivi, 77 controlli negativi** nella configurazione di riferimento)
          e la registrazione di ``denominatore_prevalenza = 803`` in ``inventario.json``;
        * Il rispetto del vincolo di tempo (< 300 s, ~24 s reali) grazie alla
          scansione in streaming a memoria costante $O(1)$ limitata a ``qc.head_reads``;
        * Il significato biologico del controllo sul motivo V4 (Gate G10), che
          compare sia nei biologici sia nei bianchi (poiché i bianchi amplificano
          contaminanti batterici reali dotati di 16S) e richiede la mediana sui
          biologici per tollerare campioni legittimi a bassa resa.
    - **Moduli sorgente coperti**:
        * ``src/amplicon16s/metadata/crosswalk.py``
        * ``src/amplicon16s/io_layer/reads.py``
        * ``src/amplicon16s/gates/g01_g15.py``
        * ``src/amplicon16s/steps/s00_validate.py``
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
    """Risolve il percorso del file YAML indicato da ``AMPLICON16S_CONFIG_DATI_REALI``."""
    grezzo = os.environ.get(VARIABILE)
    if not grezzo:
        return None
    percorso = Path(grezzo).expanduser()
    return percorso if percorso.is_file() else None


def _sorgenti_presenti(dati: dict) -> bool:
    """Verifica l'esistenza fisica su disco delle sorgenti indicate nel YAML reale."""
    obbligatorie = [dati["io"]["fastq_dir"], dati["io"]["assay_table"], dati["io"]["study_table"]]
    return all(Path(p).expanduser().exists() for p in obbligatorie)


@pytest.fixture(scope="module")
def configurazione():
    """Carica e valida la configurazione reale per OSD-734 o salta il modulo se assente."""
    percorso = _percorso_configurazione()
    if percorso is None:
        pytest.skip(f"{VARIABILE} non impostata o file assente: dati reali non disponibili")
    dati = yaml.safe_load(percorso.read_text(encoding="utf-8"))
    if not _sorgenti_presenti(dati):
        pytest.skip("le sorgenti indicate dalla configurazione non esistono")
    return carica(percorso)


@pytest.fixture(scope="module")
def inventario(configurazione):
    """Costruisce una sola volta per l'intero modulo l'inventario validato di OSD-734."""
    return esegui_gate_metadati(configurazione)


# --------------------------------------------------------------------------- #
# Crosswalk                                                                    #
# --------------------------------------------------------------------------- #


def test_960_accession_univoci(inventario):
    """
    **Obiettivo**: Verificare che l'inventario estratto da OSD-734 contenga
    esattamente 960 campioni con 960 codici accession ENA distinti.

    **Razionale Scientifico/Sistemistico**: Certifica sul dataset reale che
    l'indicizzazione per ``accession`` (a differenza del ``Sample Name`` che si
    ripete sulle repliche tecniche) mappa univocamente tutte le 960 librerie
    sequenziate sulle 10 piastre da 96 pozzetti.
    """
    assert len(inventario) == CAMPIONI_ATTESI
    assert len(set(inventario.accessioni)) == CAMPIONI_ATTESI


def test_nessun_orfano_da_nessuna_delle_due_parti(configurazione):
    """
    **Obiettivo**: Verificare su OSD-734 la perfetta simmetria tra i 960 file
    ``.fastq.gz`` in ``fastq_full/`` e le 960 righe dell'Assay Table 16S
    (zero file orfani, zero righe prive di file).

    **Razionale Scientifico/Sistemistico**: Conferma l'integrità del deposito
    locale di OSD-734 rispetto al Gate G06.
    """
    analisi = analizza(configurazione)
    assert analisi.solo_nei_file == ()
    assert analisi.solo_nell_assay == ()
    assert analisi.file_senza_accession == []
    assert analisi.assay_senza_accession == []


def test_ogni_campione_ha_il_proprio_file(inventario):
    """
    **Obiettivo**: Verificare che ciascuno dei 960 campioni dell'inventario sia
    collegato a un file ``.fastq.gz`` fisicamente esistente su disco.

    **Razionale Scientifico/Sistemistico**: Garantisce che nessun percorso
    risolto dal crosswalk punti a un link simbolico interrotto o a un file mancante.
    """
    assert all(c.file is not None and c.file.is_file() for c in inventario)


# --------------------------------------------------------------------------- #
# Inventario e classi                                                          #
# --------------------------------------------------------------------------- #


def test_conteggi_delle_classi(inventario):
    """
    **Obiettivo**: Verificare che i 960 campioni di OSD-734 siano ripartiti
    esattamente in 803 biologici, 80 controlli positivi e 77 controlli negativi.

    **Razionale Scientifico/Sistemistico**: Fissa come invariante di regressione
    il censimento esatto delle classi di OSD-734 da cui dipendono il calcolo di
    ``prev.min_samples`` (su 803 biologici), la calibrazione `KatharoSeq` (sui
    controlli positivi) e la decontaminazione `decontam` (sui controlli negativi).
    """
    assert inventario.conteggi() == {
        ClasseCampione.BIOLOGICO: BIOLOGICI_ATTESI,
        ClasseCampione.CONTROLLO_POSITIVO: POSITIVI_ATTESI,
        ClasseCampione.CONTROLLO_NEGATIVO: NEGATIVI_ATTESI,
    }


def test_il_join_ristretto_esclude_il_materiale_di_altri_assay(configurazione):
    """
    **Obiettivo**: Verificare che la Study Sample Table di OSD-734 contenga più
    di 960 righe e includa materiali estranei all'amplicon 16S (come
    ``"solvent control"`` appartenente al saggio di metagenomica), ma che il
    join ristretto li escluda producendo ``materiali_non_mappati == {}``.

    **Razionale Scientifico/Sistemistico**: Nelle indagini multi-assay di NASA
    GeneLab, la Study Table ``s_OSD-734.txt`` elenca tutti i campioni dell'intero
    studio (1.072 righe), inclusi i controlli di solvente usati solo per la
    metagenomica shotgun. Se la pipeline validasse l'intera Study Table prima
    di restringerla ai ``Sample Name`` presenti nell'Assay Table 16S, il Gate G11
    fallirebbe falsamente con ``E-S0-11`` su ``"solvent control"``.
    """
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
    """
    **Obiettivo**: Verificare che l'arricchimento tramite ``plate_well_map_960.tsv``
    assegni i 960 campioni a esattamente 10 piastre da 96 pozzetti ciascuna e
    2 corse MiSeq, senza alcun campione privo di lotto.

    **Razionale Scientifico/Sistemistico**: Conferma che la struttura fisica di
    laboratorio (10 piastre × 96 pozzetti = 960 librerie distribuite su 2 run di
    sequenziamento) è ricostruita al 100% per i modelli di errore S3 e `decontam` S12.
    """
    if configurazione.io.batch_table is None:
        pytest.skip("la configurazione non indica il file di arricchimento")
    assert len(inventario.piastre) == 10
    assert len(inventario.corse) == 2
    assert set(inventario.campioni_per_piastra().values()) == {96}
    assert inventario.senza_lotto == ()


def test_il_file_di_arricchimento_copre_dieci_piastre_da_96(configurazione):
    """
    **Obiettivo**: Verificare direttamente sul file ``plate_well_map_960.tsv``
    la presenza di 10 piastre con 96 righe ciascuna.

    **Razionale Scientifico/Sistemistico**: Accerta la coerenza interna del file
    di arricchimento indipendentemente dalla logica di join Python.
    """
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
    """
    **Obiettivo**: Verificare che le due coppie di repliche tecniche di OSD-734
    (``LAB1P3.L1_rep1``/``_rep2`` e ``NOD2S4.R6_rep1``/``_rep2``) ricevano
    ciascuna la propria piastra e la propria corsa distinte.

    **Razionale Scientifico/Sistemistico**: Dimostra sul dato reale che le
    repliche ``_rep1`` e ``_rep2`` sono state processate su piastre e corse
    diverse: grazie all'aggancio per ``accession``, ogni replica viene associata
    al modello d'errore della propria corsa effettiva.
    """
    if configurazione.io.batch_table is None:
        pytest.skip("la configurazione non indica il file di arricchimento")

    repliche = {c.nome: c for c in inventario if "_rep" in c.nome}
    assert set(repliche) == {
        "LAB1P3.L1_rep1", "LAB1P3.L1_rep2",
        "NOD2S4.R6_rep1", "NOD2S4.R6_rep2",
    }
    assert all(c.piastra is not None and c.corsa is not None for c in repliche.values())

    for radice in ("LAB1P3.L1", "NOD2S4.R6"):
        prima, seconda = repliche[f"{radice}_rep1"], repliche[f"{radice}_rep2"]
        assert prima.piastra != seconda.piastra
        assert prima.corsa != seconda.corsa


def test_senza_arricchimento_l_esecuzione_completa_lo_stesso(configurazione):
    """
    **Obiettivo**: Verificare che rimuovendo ``io.batch_table`` dalla
    configurazione di OSD-734 tutti i gate dei metadati continuino a passare
    sui 960 campioni.

    **Razionale Scientifico/Sistemistico**: Garantisce che la pipeline resti
    pienamente operativa anche qualora l'utente disponga dei soli file ISA-Tab
    scaricati da NASA GeneLab senza la mappa pozzetti supplementare.
    """
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
    """
    **Obiettivo**: Verificare che l'unica piastra di OSD-734 avente meno di 8
    controlli negativi (5 negativi e > 80 biologici contro gli 8 negativi delle
    altre 9 piastre) venga accettata senza errori né falsi allarmi.

    **Razionale Scientifico/Sistemistico**: Documenta una caratteristica reale
    del disegno sperimentale di OSD-734 (``Plate4_B``): poiché possiede 5
    controlli negativi e la soglia metodologica ``decontam.min_blanks`` è fissata
    a 5, la piastra dispone di potenza statistica sufficiente per `decontam`
    senza richiedere il ripiego globale (`E-S0-15`).
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
    assert len(inventario) == CAMPIONI_ATTESI


def test_con_arricchimento_i_moduli_sono_nove_con_i_nomi_originali(
    inventario, configurazione
):
    """
    **Obiettivo**: Verificare che con ``batch_table`` attiva vengano identificati
    tutti i 9 moduli della ISS coi nomi estesi ufficiali (incluso ``"Airlock"``).

    **Razionale Scientifico/Sistemistico**: Certifica che la colonna ``module``
    della ``batch_table`` restituisce la nomenclatura completa dei 9 compartimenti
    abitativi della Stazione Spaziale Internazionale.
    """
    if configurazione.io.batch_table is None:
        pytest.skip("la configurazione non indica il file di arricchimento")
    assert inventario.moduli == (
        "Airlock", "Columbus", "JLP", "JPM",
        "Node 1", "Node 2", "Node 3", "PMM", "US Lab",
    )


def test_senza_arricchimento_i_moduli_sono_otto_e_manca_l_airlock(configurazione):
    """
    **Obiettivo**: Verificare che senza ``batch_table`` la derivazione regex
    identifichi 8 moduli (``COL1``..``PMM1``) lasciando ``modulo = None`` sui
    16 tamponi di superficie dell'Airlock (prefisso ``A/L1``).

    **Razionale Scientifico/Sistemistico**: Documenta quantitativamente il
    limite noto della regex a 3 lettere (``^[A-Z]{3}[0-9]``) rispetto al codice
    ``A/L1`` dell'Airlock, dimostrando perché l'arricchimento con ``batch_table``
    recupera esattamente quei 16 campioni.
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
    airlock = [c for c in senza_modulo if (c.posizione or "").startswith("A/L")]
    assert len(airlock) == 16


def test_le_due_modalita_differiscono_solo_per_l_airlock(inventario, configurazione):
    """
    **Obiettivo**: Verificare che su tutti i 960 campioni di OSD-734 la presenza
    o assenza di un modulo tra la modalità con ``batch_table`` (223 senza modulo)
    e quella senza ``batch_table`` (239 senza modulo) differisca esclusivamente
    per i 16 tamponi dell'Airlock (``239 - 223 = 16``).

    **Razionale Scientifico/Sistemistico**: Dimostra che i campioni non di
    superficie (campioni d'aria, controlli negativi e positivi: 223 in totale)
    ricevono ``modulo = None`` in modo identico in entrambe le modalità.
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
    """Esegue l'intera Fase S0 (tutti i 15 gate + scrittura artefatti) una sola volta."""
    from amplicon16s.steps.s00_validate import esegui_s0

    dati = configurazione.model_dump(mode="python")
    dati["io"]["out_root"] = str(tmp_path_factory.mktemp("s0"))
    return esegui_s0(valida(dati), solleva=False)


def test_s0_supera_tutti_i_gate_sul_dataset_reale(risultato_s0):
    """
    **Obiettivo**: Verificare che l'esecuzione completa di ``esegui_s0`` su
    OSD-734 esegua e superi tutti i 15 gate (`G15` e `G01..G14`).

    **Razionale Scientifico/Sistemistico**: Costituisce il criterio di
    accettazione primario della Fase F2: l'intero dataset reale di riferimento
    supera tutti i controlli formali, relazionali e bioinformatici di ingresso.
    """
    assert risultato_s0.superata, [e.gate for e in risultato_s0.falliti]
    assert len(risultato_s0.esiti) == 15
    assert all(e.eseguito for e in risultato_s0.esiti)


def test_s0_costa_minuti_non_ore(risultato_s0):
    """
    **Obiettivo**: Verificare che l'intera Fase S0 sui 960 archivi ``.fastq.gz``
    reali completi in meno di 300 secondi (tipicamente ~24 secondi).

    **Razionale Scientifico/Sistemistico**: Se i gate G09, G10 e G13
    decomprimessero per intero i 960 file FASTQ caricandoli in RAM e riaprendoli
    tre volte, la sola validazione iniziale S0 impiegherebbe decine di minuti.
    Lo scanner single-pass ``reads.scansiona_file()`` legge in streaming solo le
    prime ``qc.head_reads = 10.000`` letture per file con complessità di memoria
    $O(1)$, abbattendo il tempo a ~24 secondi.
    """
    assert risultato_s0.secondi < 300, f"S0 ha impiegato {risultato_s0.secondi:.0f} s"


def test_s0_produce_gli_artefatti_con_checksum(risultato_s0, configurazione):
    """
    **Obiettivo**: Verificare che S0 scriva in ``01_input_validation/`` i 4
    artefatti previsti (``gates.json``, ``crosswalk.tsv``, ``inventario.json``,
    ``letture_ispezionate.tsv``) e li registri con SHA-256 valido nel manifesto.

    **Razionale Scientifico/Sistemistico**: Garantisce che gli output di S0
    siano immediatamente verificabili da ``AlberoOutput`` e pronti per essere
    consumati da S1, S2 e S10.
    """
    from amplicon16s.io_layer.artifacts import AlberoOutput, Fase

    prodotti = {p.name for p in risultato_s0.artefatti}
    assert prodotti == {
        "gates.json", "crosswalk.tsv", "inventario.json", "letture_ispezionate.tsv"
    }
    albero = AlberoOutput(risultato_s0.artefatti[0].parent.parent)
    assert albero.fase_completa(Fase.INPUT_VALIDATION)
    assert albero.non_integri(Fase.INPUT_VALIDATION) == ()


def test_s0_registra_il_denominatore_di_prevalenza(risultato_s0):
    r"""
    **Obiettivo**: Verificare che ``01_input_validation/inventario.json``
    registri ``campioni == 960`` e ``denominatore_prevalenza == 803``.

    **Razionale Scientifico/Sistemistico**: Il filtro di prevalenza in S13
    (``prev.min_fraction = 0.01``) deve essere calcolato esclusivamente sugli
    **803 campioni biologici** ($\lceil 0.01 \times 803 \rceil = 9$ campioni).
    Includere nel denominatore i 92 controlli negativi e i 65 controlli positivi
    (usando 960 anziché 803) falserebbe la soglia ecologica di prevalenza.
    """
    documento = json.loads(
        next(p for p in risultato_s0.artefatti if p.name == "inventario.json")
        .read_text(encoding="utf-8")
    )
    assert documento["campioni"] == CAMPIONI_ATTESI
    assert documento["denominatore_prevalenza"] == BIOLOGICI_ATTESI


def test_g08_non_avvisa_sulle_dieci_piastre_reali(risultato_s0, configurazione):
    """
    **Obiettivo**: Verificare che il Gate G08 passi su tutte le 10 piastre di
    OSD-734 senza emettere alcun avviso ``E-S0-15`` (``g08.avvisi == ()``).

    **Razionale Scientifico/Sistemistico**: Anche la piastra anomala ``Plate4_B``
    possiede esattamente 5 controlli negativi, raggiungendo la soglia minima
    ``decontam.min_blanks = 5`` richiesta per l'inferenza dei contaminanti per-lotto.
    """
    if configurazione.io.batch_table is None:
        pytest.skip("la configurazione non indica il file di arricchimento")
    g08 = next(e for e in risultato_s0.esiti if e.gate == "G08")
    assert g08.superato
    assert g08.avvisi == ()


def test_g10_passa_anche_sui_controlli_negativi(risultato_s0, configurazione):
    """
    **Obiettivo**: Verificare che il Gate G10 risulti superato senza violazioni
    sull'intero dataset OSD-734.

    **Razionale Scientifico/Sistemistico**: Accerta che nessun file di OSD-734
    contenga residui del primer 515F in 5' e che la mediana del motivo V4 sui
    campioni attesi superi ampiamente ``qc.min_motif_frac = 0.25``.
    """
    g10 = next(e for e in risultato_s0.esiti if e.gate == "G10")
    assert g10.superato
    assert g10.violazioni == ()


def test_il_motivo_non_distingue_segnale_e_contaminazione(risultato_s0, configurazione):
    """
    **Obiettivo**: Verificare sui dati reali di OSD-734 che il motivo conservato
    V4 (``TAC[AG].AGG..GC.AGCGTT``) sia presente ad alta frequenza sia nei
    campioni biologici sia nei controlli negativi (``min(negativi) > 0.70``), e
    che esistano singoli campioni biologici legittimi con frazione sotto soglia
    (``< 0.10``).

    **Razionale Scientifico/Sistemistico**: Questo test fissa una verità
    bioinformatica fondamentale: i controlli negativi (blank) in esperimenti a
    bassa biomassa amplificano il DNA batterico contaminante dei reagenti (*kitome*),
    e i batteri contaminanti possiedono lo stesso gene 16S rRNA (e quindi lo
    stesso motivo V4) dei batteri autentici della ISS. Di conseguenza:
    1) Il Gate G10 verifica che la libreria provenga dalla regione 16S V4
       attesa, NON la purezza del campione (compito di `decontam` in S12);
    2) Il calcolo di G10 usa la **mediana** sui campioni biologici anziché il
       minimo per singolo file proprio perché nel dataset reale esistono tamponi
       biologici a bassissima resa con ``frazione_motivo < 0.10`` che una regola
       "tutti i file sopra soglia" bloccherebbe erroneamente.
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

    sotto_soglia = [f for f in biologici if f < configurazione.qc.min_motif_frac]
    assert sotto_soglia, "ci sono biologici sotto la soglia"
    assert min(sotto_soglia) < 0.10


def test_g09_fallisce_se_il_troncamento_supera_il_minimo(configurazione, tmp_path):
    """
    **Obiettivo**: Verificare che impostando ``filter.truncLen = 152`` sul
    dataset reale OSD-734 il Gate G09 fallisca con codice ``E-S0-09`` segnalando
    la presenza di letture da ``137 bp``.

    **Razionale Scientifico/Sistemistico**: Sebbene la lunghezza nominale della
    corsa MiSeq sia 151 bp, le letture reali di OSD-734 contengono sequenze già
    trimmed fino a 137 bp: se ``filterAndTrim`` venisse lanciato con
    ``truncLen = 152`` (o ``140`` su file con letture da ``137 bp``), DADA2
    scarterebbe silenziosamente tutte le letture più corte di ``truncLen``.
    """
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
    """
    **Obiettivo**: Verificare che in tutti i 960 file di OSD-734 la frazione di
    letture che iniziano col primer forward 515F sia esattamente ``0.0``.

    **Razionale Scientifico/Sistemistico**: Conferma sperimentalmente che i
    FASTQ depositati in OSD-734 iniziano già a valle del sito di legame del
    primer 515F, giustificando il parametro predefinito ``filter.trimLeft = 0``.
    """
    righe = list(
        csv.DictReader(
            (next(p for p in risultato_s0.artefatti if p.name == "letture_ispezionate.tsv"))
            .open(encoding="utf-8"),
            delimiter="\t",
        )
    )
    assert len(righe) == CAMPIONI_ATTESI
    assert max(float(r["frazione_primer"]) for r in righe) == 0.0
