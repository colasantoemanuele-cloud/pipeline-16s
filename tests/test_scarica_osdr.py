"""Scarico di uno studio da OSDR: ``scripts/scarica_osdr.py``.

Verifica, senza rete e con risposte simulate dell'archivio, che lo script
prepari la cartella che ``init`` chiede: le tabelle ISA-Tab estratte, i soli
file di letture elencati nella tabella di assay del 16S (anche quando una
cella nomina i due file di una coppia), le sole prime letture con
``--solo-forward``; che dichiari quanti file e quanti GB prima di scaricare;
che un file gia' presente con la dimensione dichiarata non si riscarichi e uno
rimasto a meta' riprenda; che con piu' tabelle candidate chieda ``--assay``
invece di scegliere; che verifichi ogni file con l'MD5 pubblicato dall'archivio,
quando c'e', riscaricando una volta un file diverso; e che scriva l'elenco
degli SHA-256 dei file scaricati.
"""

from __future__ import annotations

import hashlib
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
        #: Per nome di file, quante risposte portano un byte alterato, con la
        #: dimensione giusta: un trasferimento guasto che la dimensione non vede.
        self.guasta: dict[str, int] = {}

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
        if self.guasta.get(nome):
            self.guasta[nome] -= 1
            dati = bytes([dati[0] ^ 0xFF]) + dati[1:]
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
    il marcatore della lettura seguito da un'altra parte del nome (``_R1_L001``),
    lo script li scarichi e avverta che lo strumento non li riconosce come
    coppie, indicando ``--solo-forward`` al posto del comando di ``init``; che
    con ``--solo-forward`` sugli stessi nomi stampi il comando di ``init``
    senza avvertenze; e che il marcatore in fondo, anche seguito dal numero di
    blocco o da ``_raw``, sia riconosciuto come tale.

    **Razionale scientifico e sistemistico**: Le due letture di un campione si
    abbinano per nome, con il marcatore in fondo. Con un'altra convenzione
    ``init`` rifiuta la cartella dopo aver letto ogni file: chi scarica deve
    saperlo subito, con la strada che funziona.
    """
    letture = {f"S{n}_R{l}_L001.fastq.gz": b"x" * (n + l) for n in (1, 2) for l in (1, 2)}
    archivio = Archivio(
        {"a_studio_16s.txt": _assay(
            [(f"C{n}", f"S{n}_R1_L001.fastq.gz,S{n}_R2_L001.fastq.gz") for n in (1, 2)])},
        letture)
    esito, righe = _lancia(archivio, tmp_path / "coppie")
    testo = "\n".join(righe)
    assert esito == 0 and len(list((tmp_path / "coppie" / "letture").iterdir())) == 4
    assert "ATTENZIONE" in testo and "--solo-forward" in testo and "p16s init" not in testo
    esito, righe = _lancia(archivio, tmp_path / "forward", "--solo-forward")
    testo = "\n".join(righe)
    assert sorted(p.name for p in (tmp_path / "forward" / "letture").iterdir()) == [
        "S1_R1_L001.fastq.gz", "S2_R1_L001.fastq.gz"]
    assert "ATTENZIONE" not in testo and "p16s init" in testo
    for nome, atteso in (("a_R1.fastq.gz", True), ("a_R2_001.fastq.gz", True),
                         ("a_1.fq", True), ("a.2.fq.gz", True), ("a_R1_raw.fastq.gz", True),
                         ("a_R2_raw.fastq.gz", True), ("a_R1_RAW.fastq.gz", True),
                         ("a_R2_RAW.fastq.gz", True), ("a_R1_L001.fastq.gz", False),
                         ("a_R1_trimmed.fastq.gz", False), ("campione1.fastq.gz", False)):
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


# --------------------------------------------------------------------------- #
# Le impronte dei file: l'MD5 pubblicato e l'elenco degli SHA-256              #
# --------------------------------------------------------------------------- #

PREFISSO = f"GLDS-{STUDIO}_Amplicon_"
ELENCO_MD5 = f"GLDS-{STUDIO}_Amplicon_raw_md5sum_GLAmpSeq.txt"


def _md5(dati: bytes) -> str:
    return hashlib.md5(dati).hexdigest()


def _studio_con_md5(letture: dict[str, bytes], md5: dict[str, str]) -> Archivio:
    """Uno studio i cui file portano il prefisso dell'archivio, con un elenco
    MD5 che li nomina senza prefisso, come quelli pubblicati da OSDR."""
    elenco = "".join(f"{valore}  {nome}\n" for nome, valore in md5.items())
    righe = [(f"C{i}", PREFISSO + nome) for i, nome in enumerate(letture)]
    return Archivio({"a_studio_16s.txt": _assay(righe)},
                    {**{PREFISSO + n: d for n, d in letture.items()},
                     ELENCO_MD5: elenco.encode()})


def _impronte(cartella: Path) -> list[dict[str, str]]:
    righe = (cartella / scarica_osdr.ELENCO_IMPRONTE).read_text(encoding="utf-8").splitlines()
    intestazione = righe[0].split("\t")
    return [dict(zip(intestazione, r.split("\t"))) for r in righe[1:]]


def test_con_l_md5_pubblicato_ogni_file_e_verificato_e_l_elenco_ha_lo_sha256(
        tmp_path, senza_attese):
    """
    **Obiettivo**: Verificare che, con un elenco MD5 fra i file dello studio
    che nomina i file senza il prefisso dell'archivio, lo script lo legga,
    abbini ogni file al suo MD5, dica quanti file ne hanno uno, scarichi i
    file una volta sola e scriva ``impronte.tsv`` con percorso, byte, SHA-256
    calcolato, MD5 pubblicato ed elenco di provenienza per l'archivio ISA-Tab
    e per ogni lettura; che un file che l'elenco non nomina abbia le colonne
    dell'MD5 vuote; che un secondo lancio non chieda alcun file di letture e
    riscriva lo stesso elenco, byte per byte.

    **Razionale scientifico e sistemistico**: La dimensione non vede un
    trasferimento guasto che ne lascia invariata la lunghezza: l'MD5 che
    l'archivio pubblica e' la sola verifica del contenuto. L'elenco degli
    SHA-256 permette di verificare i file in seguito senza l'archivio.
    """
    letture = {"A_raw.fastq.gz": b"letture a " * 50, "B_raw.fastq.gz": b"letture b " * 30,
               "C_raw.fastq.gz": b"letture c " * 10}
    archivio = _studio_con_md5(
        letture, {n: _md5(d) for n, d in letture.items() if n != "C_raw.fastq.gz"})
    cartella = tmp_path / "dati"
    esito, righe = _lancia(archivio, cartella)
    assert esito == 0
    testo = "\n".join(righe)
    assert f"Impronte: 2 file su 3 con l'MD5 pubblicato da OSDR ({ELENCO_MD5})" in testo
    assert "per 1 l'archivio non pubblica un MD5" in testo
    assert sorted(n for n, _ in archivio.richieste if n.endswith(".fastq.gz")) == sorted(
        PREFISSO + n for n in letture)
    attese = {f"letture/{PREFISSO}{n}": d for n, d in letture.items()}
    attese[f"metadati/{ARCHIVIO}"] = archivio.contenuti[ARCHIVIO]
    voci = _impronte(cartella)
    assert [v["file"] for v in voci] == sorted(attese)
    for voce in voci:
        dati = attese[voce["file"]]
        assert voce["byte"] == str(len(dati))
        assert voce["sha256"] == hashlib.sha256(dati).hexdigest()
        con_md5 = voce["file"].startswith("letture/") and not voce["file"].endswith("C_raw.fastq.gz")
        assert voce["md5_pubblicato"] == (_md5(dati) if con_md5 else "")
        assert voce["elenco_md5"] == (ELENCO_MD5 if con_md5 else "")
    assert not list(cartella.rglob("*.parziale"))

    primo = (cartella / scarica_osdr.ELENCO_IMPRONTE).read_bytes()
    archivio.richieste.clear()
    esito, righe = _lancia(archivio, cartella)
    assert esito == 0
    assert not [n for n, _ in archivio.richieste if n.endswith(".fastq.gz")]
    assert "gia' presenti 3, da scaricare 0" in "\n".join(righe)
    assert (cartella / scarica_osdr.ELENCO_IMPRONTE).read_bytes() == primo


def test_un_md5_diverso_si_riscarica_una_volta_poi_e_un_errore(tmp_path, senza_attese):
    """
    **Obiettivo**: Verificare che un file arrivato con la dimensione giusta ma
    un MD5 diverso da quello pubblicato sia cancellato e scaricato di nuovo da
    capo, una volta, e tenuto se il secondo trasferimento coincide; che se
    anche il secondo e' diverso il file non resti sul disco ne' nell'elenco
    delle impronte, lo script lo dica nominando l'elenco MD5 e chiuda con il
    codice del guasto, senza fermare gli altri file; e che un file gia'
    presente con la dimensione giusta ma un MD5 diverso non conti come
    completo e venga riscaricato.

    **Razionale scientifico e sistemistico**: Un FASTQ alterato si legge
    spesso senza errori e cambia i conteggi delle letture: un file che non
    corrisponde all'impronta dichiarata non deve entrare nell'analisi con il
    nome giusto. Un solo nuovo tentativo distingue un guasto del trasferimento
    da un file pubblicato diverso dal suo MD5, che rilanciare non corregge.
    """
    letture = {"A_raw.fastq.gz": b"letture a " * 50, "B_raw.fastq.gz": b"letture b " * 30}
    archivio = _studio_con_md5(letture, {n: _md5(d) for n, d in letture.items()})
    archivio.guasta[PREFISSO + "A_raw.fastq.gz"] = 1
    cartella = tmp_path / "una_volta"
    esito, _ = _lancia(archivio, cartella)
    assert esito == 0
    assert [d for n, d in archivio.richieste if n == PREFISSO + "A_raw.fastq.gz"] == [0, 0]
    assert (cartella / "letture" / (PREFISSO + "A_raw.fastq.gz")).read_bytes() == (
        letture["A_raw.fastq.gz"])

    sbagliato = _studio_con_md5(letture, {"A_raw.fastq.gz": "0" * 32,
                                          "B_raw.fastq.gz": _md5(letture["B_raw.fastq.gz"])})
    cartella = tmp_path / "sempre"
    esito, righe = _lancia(sbagliato, cartella)
    testo = "\n".join(righe)
    assert esito == scarica_osdr.USCITA_GUASTO
    assert [d for n, d in sbagliato.richieste if n == PREFISSO + "A_raw.fastq.gz"] == [0, 0]
    assert not (cartella / "letture" / (PREFISSO + "A_raw.fastq.gz")).exists()
    assert (cartella / "letture" / (PREFISSO + "B_raw.fastq.gz")).is_file()
    assert ELENCO_MD5 in testo and "1 file con un MD5 diverso" in testo
    assert [v["file"] for v in _impronte(cartella)] == [
        f"letture/{PREFISSO}B_raw.fastq.gz", f"metadati/{ARCHIVIO}"]
    assert not list(cartella.rglob("*.parziale"))

    cartella = tmp_path / "presente"
    (cartella / "letture").mkdir(parents=True)
    alterato = bytearray(letture["B_raw.fastq.gz"])
    alterato[0] ^= 0xFF
    (cartella / "letture" / (PREFISSO + "B_raw.fastq.gz")).write_bytes(bytes(alterato))
    esito, righe = _lancia(_studio_con_md5(letture, {n: _md5(d) for n, d in letture.items()}),
                           cartella)
    assert esito == 0 and "gia' presenti 0, da scaricare 2" in "\n".join(righe)
    assert (cartella / "letture" / (PREFISSO + "B_raw.fastq.gz")).read_bytes() == (
        letture["B_raw.fastq.gz"])


def test_senza_md5_pubblicato_si_verifica_la_dimensione_e_l_elenco_ha_lo_sha256(
        tmp_path, senza_attese, monkeypatch):
    """
    **Obiettivo**: Verificare che, con uno studio che non pubblica elenchi MD5,
    lo scarico proceda verificando la sola dimensione, lo dica, e scriva
    comunque ``impronte.tsv`` con lo SHA-256 di ogni file e le colonne
    dell'MD5 vuote, senza un file rimasto a meta' per un guasto; e che le regole di abbinamento e di lettura degli elenchi
    reggano i casi limite: nome identico, suffisso piu' lungo dopo un trattino
    basso, nome senza separatore (non abbinato), nome con due MD5 diversi
    (ambiguo, quindi nessun MD5), forma binaria di ``md5sum``, righe non valide
    ignorate, elenco riconosciuto dalla sottocategoria dell'archivio.

    **Razionale scientifico e sistemistico**: Non tutti gli studi pubblicano
    l'MD5: lo scarico non deve fermarsi, ma dichiarare che cosa ha potuto
    verificare. Un abbinamento sbagliato fra nome e MD5 rifiuterebbe file
    integri; uno ambiguo non deve decidere a caso.
    """
    archivio = _studio_appaiato()
    cartella = tmp_path / "dati"
    esito, righe = _lancia(archivio, cartella)
    assert esito == 0
    assert "Impronte: 0 file su 6 con l'MD5 pubblicato da OSDR; per 6 l'archivio non " in (
        "\n".join(righe))
    voci = _impronte(cartella)
    assert len(voci) == 7
    for voce in voci:
        dati = (cartella / voce["file"]).read_bytes()
        assert voce["sha256"] == hashlib.sha256(dati).hexdigest()
        assert voce["md5_pubblicato"] == voce["elenco_md5"] == ""

    # Un file rimasto a meta' per un guasto di rete non entra nell'elenco.
    monkeypatch.setattr(scarica_osdr, "TENTATIVI", 1)
    guasto = _studio_appaiato()
    guasto.interrompi["S1_R1_001.fastq.gz"] = 3
    esito, _ = _lancia(guasto, tmp_path / "guasto")
    assert esito == scarica_osdr.USCITA_GUASTO
    assert (tmp_path / "guasto" / "letture" / "S1_R1_001.fastq.gz.parziale").is_file()
    assert "letture/S1_R1_001.fastq.gz" not in [v["file"] for v in _impronte(tmp_path / "guasto")]
    assert len(_impronte(tmp_path / "guasto")) == 6

    pubblicati = {"ERX1_raw.fastq.gz": ("1" * 32, "e1.txt"),
                  "X_ERX1_raw.fastq.gz": ("2" * 32, "e1.txt"),
                  "doppio.fastq.gz": (None, "e2.txt")}
    trova = scarica_osdr.md5_del_file
    assert trova("ERX1_raw.fastq.gz", pubblicati) == ("1" * 32, "e1.txt")
    assert trova("GLDS-1_Amplicon_ERX1_raw.fastq.gz", pubblicati) == ("1" * 32, "e1.txt")
    assert trova("GLDS-1_Amplicon_X_ERX1_raw.fastq.gz", pubblicati) == ("2" * 32, "e1.txt")
    assert trova("GLDS-1_AmpliconERX1_raw.fastq.gz", pubblicati) is None
    assert trova("GLDS-1_doppio.fastq.gz", pubblicati) is None
    assert scarica_osdr.leggi_elenco_md5(
        f"{'A' * 32} *a.fq.gz\nnon e' una riga\n{'b' * 32}  b c.fq.gz\n") == {
        "a.fq.gz": ["a" * 32], "b c.fq.gz": ["b" * 32]}
    due_elenchi = Archivio({}, {"a_md5sum_1.txt": f"{'1' * 32}  f.fq.gz\n".encode(),
                                "b_md5sum_2.txt": f"{'2' * 32}  f.fq.gz\n{'3' * 32}  g.fq.gz\n"
                                .encode()})
    letti = scarica_osdr.md5_pubblicati(
        STUDIO, scarica_osdr.file_dello_studio(STUDIO, due_elenchi.apri), due_elenchi.apri)
    assert letti == {"f.fq.gz": (None, "b_md5sum_2.txt"), "g.fq.gz": ("3" * 32, "b_md5sum_2.txt")}
    assert scarica_osdr.md5_del_file("f.fq.gz", letti) is None
    assert scarica_osdr.elenchi_md5({
        "x.txt": {"subcategory": "MD5 Checksums"}, "y_md5sum_z.txt": {},
        "letture.fastq.gz": {"subcategory": "Raw sequence data"}}) == ["x.txt", "y_md5sum_z.txt"]
