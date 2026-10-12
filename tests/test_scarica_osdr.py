"""Scarico di uno studio da OSDR: ``scripts/scarica_osdr.py``.

Verifica, senza rete e con risposte simulate dell'archivio, che lo script
prepari la cartella che ``init`` chiede: le tabelle ISA-Tab estratte, i soli
file di letture elencati nella tabella di assay del 16S (anche quando una
cella nomina i due file di una coppia), le sole prime letture con
``--solo-forward``; che dichiari quanti file e quanti GB prima di scaricare;
che un file gia' presente con la dimensione dichiarata non si riscarichi e uno
rimasto a meta' riprenda; e che con piu' tabelle candidate chieda ``--assay``
invece di scegliere.
"""

from __future__ import annotations

import importlib.util
import io
import json
import zipfile
from pathlib import Path

import pytest

RADICE = Path(__file__).parents[1]

_spec = importlib.util.spec_from_file_location(
    "scarica_osdr", RADICE / "scripts" / "scarica_osdr.py")
scarica_osdr = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(scarica_osdr)

STUDIO = "4242"
ARCHIVIO = f"OSD-{STUDIO}_metadata_OSD-{STUDIO}-ISA.zip"


def _assay(righe: list[tuple[str, str]]) -> str:
    """Una tabella di assay con il nome del campione e la cella dei file grezzi."""
    return "Sample Name\tRaw Data File\tAltro\n" + "".join(
        f"{nome}\t{file}\tx\n" for nome, file in righe)


class Archivio:
    """Un archivio simulato: l'elenco dei file di uno studio e il loro contenuto.

    Conta le richieste per indirizzo e registra da quale byte ciascuna parte,
    cosi' i test verificano che cosa e' stato chiesto davvero.
    """

    def __init__(self, tabelle: dict[str, str], letture: dict[str, bytes]) -> None:
        zip_in_memoria = io.BytesIO()
        with zipfile.ZipFile(zip_in_memoria, "w") as z:
            for nome, testo in tabelle.items():
                z.writestr(nome, testo)
        self.contenuti = {ARCHIVIO: zip_in_memoria.getvalue(), **letture}
        self.richieste: list[tuple[str, int]] = []
        #: Per nome di file, dopo quanti byte la prima risposta si interrompe.
        self.interrompi: dict[str, int] = {}

    def elenco(self) -> bytes:
        return json.dumps({"studies": {f"OSD-{STUDIO}": {"study_files": [
            {"file_name": nome, "file_size": len(dati),
             "remote_url": f"/scarico?file={nome}"}
            for nome, dati in self.contenuti.items()]}}}).encode()

    def apri(self, url: str, da: int = 0):
        if url.endswith(f"/files/{STUDIO}"):
            return _Risposta(self.elenco(), 200)
        if "/files/" in url:
            return _Risposta(json.dumps({"studies": {}}).encode(), 200)
        nome = url.split("file=")[1]
        self.richieste.append((nome, da))
        dati = self.contenuti[nome][da:]
        taglio = self.interrompi.pop(nome, None)
        return _Risposta(dati, 206 if da else 200, taglio)


class _Risposta:
    """La risposta a una richiesta: un flusso, che puo' cadere a meta'."""

    def __init__(self, dati: bytes, stato: int, cade_dopo: int | None = None) -> None:
        self._flusso = io.BytesIO(dati if cade_dopo is None else dati[:cade_dopo])
        self._cade = cade_dopo is not None
        self.status = stato

    def read(self, quanti: int = -1) -> bytes:
        blocco = self._flusso.read(quanti)
        if not blocco and self._cade:
            raise ConnectionResetError("connessione caduta")
        return blocco

    def __enter__(self):
        return self

    def __exit__(self, *_):
        return False


@pytest.fixture
def senza_attese(monkeypatch):
    """Toglie le pause fra un tentativo e l'altro."""
    monkeypatch.setattr(scarica_osdr.time, "sleep", lambda _: None)


def _lancia(archivio: Archivio, cartella: Path, *argomenti: str) -> tuple[int, list[str]]:
    righe: list[str] = []
    args = scarica_osdr.build_parser().parse_args([f"OSD-{STUDIO}", str(cartella), *argomenti])
    return scarica_osdr.esegui(args, righe.append, archivio.apri), righe


def _studio_appaiato() -> Archivio:
    letture = {f"S{n}_R{l}_001.fastq.gz": f"letture {n} {l} ".encode() * (n + l)
               for n in (1, 2, 3) for l in (1, 2)}
    return Archivio(
        {"i_Investigation.txt": "x\n", f"s_OSD-{STUDIO}.txt": "Sample Name\nA\nB\nC\n",
         "a_studio_amplicon-sequencing_16s_illumina.txt": _assay(
             [(f"C{n}", f"S{n}_R1_001.fastq.gz,S{n}_R2_001.fastq.gz") for n in (1, 2, 3)]),
         "a_studio_amplicon-sequencing_its_illumina.txt": _assay(
             [("F1", "ITS1_R1_001.fastq.gz,ITS1_R2_001.fastq.gz")])},
        {**letture, "ITS1_R1_001.fastq.gz": b"its", "ITS1_R2_001.fastq.gz": b"its",
         "tabella_elaborata.tsv": b"altro"})


def test_scarica_metadati_e_le_sole_letture_della_tabella_del_16s(tmp_path, senza_attese):
    """
    **Obiettivo**: Verificare che da uno studio con due tabelle di assay (16S
    e ITS) lo script estragga le tabelle in ``metadati/``, riconosca dal nome
    quella del 16S e scarichi in ``letture/`` i soli file che essa elenca, i
    due di ogni coppia nominati nella stessa cella, con il contenuto
    dell'archivio; che dichiari prima quanti file e quanti GB; che stampi il
    comando di ``init`` con la tabella scelta; e che un secondo lancio non
    chieda di nuovo alcun file di letture.

    **Razionale scientifico e sistemistico**: Uno studio pubblica insieme piu'
    ampliconi e i file elaborati: nella cartella delle letture devono entrare
    i soli file grezzi del 16S, altrimenti la validazione si ferma su file
    senza riga o, peggio, analizza un altro amplicone. Lo scarico rilanciato
    non deve ripetere cio' che e' gia' integro sul disco.
    """
    archivio = _studio_appaiato()
    esito, righe = _lancia(archivio, tmp_path / "dati")
    assert esito == 0
    assert sorted(p.name for p in (tmp_path / "dati" / "metadati").iterdir()) == sorted([
        ARCHIVIO, "i_Investigation.txt", f"s_OSD-{STUDIO}.txt",
        "a_studio_amplicon-sequencing_16s_illumina.txt",
        "a_studio_amplicon-sequencing_its_illumina.txt"])
    scaricati = sorted((tmp_path / "dati" / "letture").iterdir())
    assert [p.name for p in scaricati] == sorted(
        f"S{n}_R{l}_001.fastq.gz" for n in (1, 2, 3) for l in (1, 2))
    assert all(p.read_bytes() == archivio.contenuti[p.name] for p in scaricati)
    testo = "\n".join(righe)
    assert "Letture: 6 file, 0,00 GB; gia' presenti 0, da scaricare 6" in testo
    assert testo.index("Letture: 6 file") < testo.index("S1_R1_001.fastq.gz")
    assert "--assay a_studio_amplicon-sequencing_16s_illumina.txt" in testo
    assert "p16s init --reads /dati/letture --metadata /dati/metadati" in testo

    archivio.richieste.clear()
    esito, righe = _lancia(archivio, tmp_path / "dati")
    assert esito == 0
    assert [n for n, _ in archivio.richieste if n != ARCHIVIO] == []
    assert "gia' presenti 6, da scaricare 0" in "\n".join(righe)


def test_solo_forward_scarica_la_prima_lettura_di_ogni_coppia(tmp_path, senza_attese):
    """
    **Obiettivo**: Verificare che con ``--solo-forward`` di ogni cella con i
    due file di una coppia si scarichi il solo file della prima lettura, che
    il file della seconda non sia mai chiesto all'archivio, e che la regola
    riconosca le coppie dal marcatore della lettura nel nome (1 e 2 dopo un
    separatore o una R), lasciando com'e' un file senza compagno e non
    scambiando per marcatore il numero del campione o della corsia.

    **Razionale scientifico e sistemistico**: L'analisi delle sole forward e'
    la scelta quando le reverse sono corte o di qualita' bassa: la cartella
    deve contenere un file per campione, perche' con i due file ``init``
    sceglie la modalita' appaiata.
    """
    archivio = _studio_appaiato()
    esito, righe = _lancia(archivio, tmp_path / "dati", "--solo-forward")
    assert esito == 0
    assert sorted(p.name for p in (tmp_path / "dati" / "letture").iterdir()) == [
        f"S{n}_R1_001.fastq.gz" for n in (1, 2, 3)]
    assert not any("_R2_" in nome for nome, _ in archivio.richieste)
    assert "3 file" in "\n".join(righe) and "le sole prime letture" in "\n".join(righe)

    assert scarica_osdr.forward_di(["a_1.fq.gz", "a_2.fq.gz"]) == ["a_1.fq.gz"]
    assert scarica_osdr.forward_di(["x.R2.fastq", "x.R1.fastq"]) == ["x.R1.fastq"]
    assert scarica_osdr.forward_di(["solo_R2.fastq.gz"]) == ["solo_R2.fastq.gz"]
    # Il numero del campione o della corsia non e' il marcatore della lettura.
    assert scarica_osdr.forward_di(["c1_R1.fq", "c2_R1.fq"]) == ["c1_R1.fq", "c2_R1.fq"]
    corsie = ["S1_L001_R1_001.fastq.gz", "S1_L002_R1_001.fastq.gz",
              "S1_L001_R2_001.fastq.gz", "S1_L002_R2_001.fastq.gz"]
    assert scarica_osdr.forward_di(corsie) == corsie[:2]


def test_avverte_se_le_coppie_non_hanno_il_marcatore_in_fondo_al_nome(tmp_path, senza_attese):
    """
    **Obiettivo**: Verificare che, quando i due file di ogni campione portano
    il marcatore della lettura seguito da un'altra parte del nome (``_R1_raw``),
    lo script li scarichi e avverta che lo strumento non li riconosce come
    coppie, indicando ``--solo-forward`` al posto del comando di ``init``; che
    con ``--solo-forward`` sugli stessi nomi stampi il comando di ``init``
    senza avvertenze; e che il marcatore in fondo, anche seguito dal numero di
    blocco, sia riconosciuto come tale.

    **Razionale scientifico e sistemistico**: Le due letture di un campione si
    abbinano per nome, con il marcatore in fondo. Con un'altra convenzione
    ``init`` rifiuta la cartella dopo aver letto ogni file: chi scarica deve
    saperlo subito, con la strada che funziona.
    """
    letture = {f"S{n}_R{l}_raw.fastq.gz": b"x" * (n + l) for n in (1, 2) for l in (1, 2)}
    archivio = Archivio(
        {"a_studio_16s.txt": _assay(
            [(f"C{n}", f"S{n}_R1_raw.fastq.gz,S{n}_R2_raw.fastq.gz") for n in (1, 2)])},
        letture)
    esito, righe = _lancia(archivio, tmp_path / "coppie")
    testo = "\n".join(righe)
    assert esito == 0 and len(list((tmp_path / "coppie" / "letture").iterdir())) == 4
    assert "ATTENZIONE" in testo and "--solo-forward" in testo and "p16s init" not in testo
    esito, righe = _lancia(archivio, tmp_path / "forward", "--solo-forward")
    testo = "\n".join(righe)
    assert sorted(p.name for p in (tmp_path / "forward" / "letture").iterdir()) == [
        "S1_R1_raw.fastq.gz", "S2_R1_raw.fastq.gz"]
    assert "ATTENZIONE" not in testo and "p16s init" in testo
    for nome, atteso in (("a_R1.fastq.gz", True), ("a_R2_001.fastq.gz", True),
                         ("a_1.fq", True), ("a.2.fq.gz", True), ("a_R1_raw.fastq.gz", False),
                         ("campione1.fastq.gz", False)):
        assert scarica_osdr.marcatore_in_fondo(nome) is atteso, nome


def test_uno_scarico_interrotto_riprende_dal_punto_in_cui_era(tmp_path, senza_attese):
    """
    **Obiettivo**: Verificare che, se la risposta dell'archivio cade a meta'
    di un file, lo script ritenti chiedendo il solo resto (dal byte gia'
    scritto), che il file finale sia identico a quello dell'archivio, e che un
    file parziale non porti mai il nome definitivo; e che un file presente ma
    di dimensione diversa da quella dichiarata sia scaricato di nuovo.

    **Razionale scientifico e sistemistico**: Un FASTQ compresso troncato si
    legge senza errori fino al punto del taglio: con il nome giusto entrerebbe
    nell'analisi con meno letture. Il nome definitivo arriva solo con la
    dimensione dichiarata dall'archivio.
    """
    grande = bytes(range(256)) * 40
    archivio = Archivio(
        {"a_studio_16s.txt": _assay([("C1", "grande.fastq.gz"), ("C2", "altro.fastq.gz")])},
        {"grande.fastq.gz": grande, "altro.fastq.gz": b"contenuto giusto"})
    archivio.interrompi["grande.fastq.gz"] = 1000
    letture = tmp_path / "dati" / "letture"
    letture.mkdir(parents=True)
    (letture / "altro.fastq.gz").write_bytes(b"troncato")
    esito, _ = _lancia(archivio, tmp_path / "dati", "--connessioni", "1")
    assert esito == 0
    assert ("grande.fastq.gz", 0) in archivio.richieste
    assert ("grande.fastq.gz", 1000) in archivio.richieste
    assert (letture / "grande.fastq.gz").read_bytes() == grande
    assert (letture / "altro.fastq.gz").read_bytes() == b"contenuto giusto"
    assert not list(letture.glob("*.parziale"))


def test_con_piu_tabelle_candidate_chiede_assay_e_non_sceglie(tmp_path, senza_attese):
    """
    **Obiettivo**: Verificare che con due tabelle di assay che hanno «16S» nel
    nome, o con piu' tabelle di cui nessuna lo ha, lo script si rifiuti
    nominandole e chiedendo ``--assay``, senza scaricare letture, e che con
    ``--assay`` usi quella indicata; che uno studio che l'archivio non elenca,
    un identificativo malformato e una tabella senza la colonna dei file
    grezzi siano rifiutati con un messaggio che dice che cosa fare.

    **Razionale scientifico e sistemistico**: Scegliere in silenzio fra due
    tabelle scaricherebbe le letture di un altro esperimento dello stesso
    studio: la scelta resta a chi conosce lo studio, con i nomi per farla.
    """
    tabelle = {
        "a_studio_16s_v4.txt": _assay([("C1", "v4.fastq.gz")]),
        "a_studio_16s_v1v3.txt": _assay([("C1", "v1v3.fastq.gz")])}
    archivio = Archivio(tabelle, {"v4.fastq.gz": b"quattro", "v1v3.fastq.gz": b"uno"})
    with pytest.raises(scarica_osdr.Rifiuto) as rifiuto:
        _lancia(archivio, tmp_path / "due")
    assert "--assay" in str(rifiuto.value)
    assert all(nome in str(rifiuto.value) for nome in tabelle)
    assert not (tmp_path / "due" / "letture").exists()
    esito, _ = _lancia(archivio, tmp_path / "due", "--assay", "a_studio_16s_v4.txt")
    assert esito == 0
    assert [p.name for p in (tmp_path / "due" / "letture").iterdir()] == ["v4.fastq.gz"]

    senza_nome = Archivio({"a_uno.txt": _assay([("C1", "v4.fastq.gz")]),
                           "a_due.txt": _assay([("C1", "v4.fastq.gz")])},
                          {"v4.fastq.gz": b"quattro"})
    with pytest.raises(scarica_osdr.Rifiuto, match="--assay"):
        _lancia(senza_nome, tmp_path / "senza")

    with pytest.raises(scarica_osdr.Rifiuto, match="non elenca alcun file"):
        scarica_osdr.file_dello_studio("1", archivio.apri)
    with pytest.raises(scarica_osdr.Rifiuto, match="OSD-NNN"):
        scarica_osdr.numero_dello_studio("studio 12b")
    assert scarica_osdr.numero_dello_studio("osd-0042") == "42"
    senza_colonna = Archivio({"a_studio_16s.txt": "Sample Name\tAltro\nC1\tx\n"}, {})
    with pytest.raises(scarica_osdr.Rifiuto, match="Raw Data File"):
        _lancia(senza_colonna, tmp_path / "colonna")
