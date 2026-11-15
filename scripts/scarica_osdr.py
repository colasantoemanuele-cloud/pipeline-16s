#!/usr/bin/env python3
"""Scarica da OSDR i metadati e le letture 16S di uno studio.

Uso, dalla radice di un clone::

    python3 scripts/scarica_osdr.py OSD-NNN CARTELLA [--assay NOME] [--solo-forward]

Nella cartella indicata prepara cio' che ``init`` chiede:

* ``metadati/``: le tabelle ISA-Tab dello studio (investigazione, studio,
  assay), estratte dall'archivio pubblicato;
* ``letture/``: i file di letture elencati nella tabella di assay del 16S.

Razionale sistemistico:

* **solo libreria standard**: lo script gira prima che l'ambiente della
  pipeline esista, sul sistema o dentro l'immagine;
* **l'elenco dei file viene dall'archivio** (dimensione compresa): prima di
  scaricare si dichiara quanti file e quanti GB, e un file gia' presente con
  la dimensione dichiarata non si riscarica. Lo script si puo' rilanciare dopo
  un'interruzione: riprende i file rimasti a meta' (``<nome>.parziale``,
  intestazione HTTP ``Range``) e salta quelli completi;
* **il nome definitivo si assegna solo a trasferimento concluso**, con la
  dimensione verificata: un'interruzione non lascia un file incompleto con il
  nome giusto;
* **l'impronta si verifica con l'MD5 che OSDR pubblica**: l'elenco JSON dei
  file non porta checksum, ma molti studi pubblicano fra i loro file uno o piu'
  elenchi MD5 dei file grezzi (sottocategoria ``MD5 Checksums``, nomi come
  ``GLDS-NNN_Amplicon_raw_md5sum_GLAmpSeq.txt``), in cui ogni file e' nominato
  senza il prefisso ``GLDS-NNN_..._`` che OSDR aggiunge. Un file con l'MD5
  pubblicato e' completo solo se dimensione e MD5 coincidono; uno scaricato
  con un MD5 diverso si cancella e si riscarica una volta, poi e' un errore.
  Senza MD5 pubblicato resta la verifica della dimensione. Non si usa l'MD5 di
  ENA: GeneLab ricomprime i file, che non coincidono byte per byte con quelli
  di ENA;
* **un elenco delle impronte** (``impronte.tsv`` nella cartella dei dati):
  l'SHA-256 calcolato di ogni file scaricato, con la dimensione e l'MD5
  pubblicato se c'e', per verificare i file in seguito senza l'archivio;
* **poche connessioni alla volta**: l'archivio e' pubblico, e aprirne molte in
  parallelo e' un modo di abusarne.

Razionale sperimentale: uno studio puo' avere piu' tabelle di assay (16S, ITS,
metagenomica). Le letture da analizzare sono quelle che la tabella del 16S
elenca, riga per riga, nella colonna dei file grezzi; una cella puo' nominare
i due file di una coppia. Con ``--solo-forward`` si scarica di ogni coppia il
solo file della prima lettura.
"""

from __future__ import annotations

import argparse
import csv
import functools
import hashlib
import http.client
import json
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
import zipfile
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Callable, Dict, List, Optional, Sequence, Tuple

#: L'archivio: l'elenco dei file di uno studio e lo scarico di un file.
SITO = "https://osdr.nasa.gov"
ELENCO = SITO + "/osdr/data/osd/files/{numero}"
SCARICO = SITO + "/geode-py/ws/studies/OSD-{numero}/download?source=datamanager&file={nome}"

#: La colonna ISA-Tab che elenca i file grezzi di ogni riga di una tabella di assay.
COLONNA_DEI_FILE = "Raw Data File"
#: Le estensioni dei file di letture.
ESTENSIONI = (".fastq.gz", ".fq.gz", ".fastq", ".fq")

#: La sottocategoria con cui l'elenco dell'archivio marca gli elenchi MD5 dei
#: file grezzi; li si riconosce anche dal nome, che contiene ``md5sum``.
SOTTOCATEGORIA_MD5 = "MD5 Checksums"
#: L'elenco delle impronte dei file scaricati, nella cartella dei dati.
ELENCO_IMPRONTE = "impronte.tsv"
#: Quante volte si scarica un file il cui MD5 non coincide con quello pubblicato.
SCARICHI_PER_MD5 = 2

BLOCCO = 1024 * 1024
TEMPO_MASSIMO = 60
CONNESSIONI = 4
TENTATIVI = 5
GUASTI_DI_RETE = (urllib.error.URLError, OSError, http.client.HTTPException)

USCITA_RIFIUTO = 2
USCITA_GUASTO = 1


class Rifiuto(Exception):
    """La richiesta non si puo' soddisfare: il messaggio dice perche' e che cosa fare."""


class ImprontaDiversa(OSError):
    """Un file scaricato due volte ha un MD5 diverso da quello pubblicato."""


def numero_dello_studio(testo: str) -> str:
    """Il numero di uno studio scritto come ``OSD-NNN``, ``osdNNN`` o il solo numero."""
    trovato = re.fullmatch(r"(?i)(?:OSD[-_ ]?)?([0-9]+)", testo.strip())
    if trovato is None:
        raise Rifiuto(f"{testo!r} non e' l'identificativo di uno studio: si scrive OSD-NNN")
    return str(int(trovato.group(1)))


def apri(url: str, da: int = 0):
    """Apre l'indirizzo, dal byte indicato se non e' zero."""
    richiesta = urllib.request.Request(url, headers={"User-Agent": "scarica_osdr"})
    if da:
        richiesta.add_header("Range", f"bytes={da}-")
    return urllib.request.urlopen(richiesta, timeout=TEMPO_MASSIMO)


def file_dello_studio(numero: str, apri_url: Callable = apri) -> Dict[str, dict]:
    """I file pubblicati dello studio, per nome, dall'elenco dell'archivio."""
    with apri_url(ELENCO.format(numero=numero)) as risposta:
        documento = json.load(risposta)
    studi = documento.get("studies") or {}
    studio = studi.get(f"OSD-{numero}")
    if not studio:
        raise Rifiuto(f"l'archivio non elenca alcun file per OSD-{numero}: controlla il numero")
    return {f["file_name"]: f for f in studio.get("study_files", [])}


def archivio_dei_metadati(file: Dict[str, dict]) -> str:
    """Il nome dell'archivio ISA-Tab fra i file dello studio."""
    trovati = sorted(n for n in file if n.lower().endswith(".zip") and "isa" in n.lower()
                     and "metadata" in n.lower())
    if len(trovati) != 1:
        raise Rifiuto(
            "fra i file dello studio non c'e' un solo archivio ISA-Tab (trovati: "
            + (", ".join(trovati) or "nessuno") + ")")
    return trovati[0]


def elenchi_md5(file: Dict[str, dict]) -> List[str]:
    """I nomi degli elenchi MD5 fra i file dello studio, in ordine."""
    return sorted(n for n, voce in file.items()
                  if voce.get("subcategory") == SOTTOCATEGORIA_MD5
                  or ("md5sum" in n.lower() and n.lower().endswith(".txt")))


def leggi_elenco_md5(testo: str) -> Dict[str, List[str]]:
    """Gli MD5 delle righe ``<md5>  <nome>`` (o ``<md5> *<nome>``, la forma
    binaria di ``md5sum``) di un elenco, per nome: piu' d'uno se il nome
    compare piu' volte. Le altre righe si ignorano."""
    voci: Dict[str, List[str]] = {}
    for riga in testo.splitlines():
        trovata = re.fullmatch(r"\s*([0-9a-fA-F]{32})\s+\*?(.+?)\s*", riga)
        if trovata:
            voci.setdefault(trovata.group(2), []).append(trovata.group(1).lower())
    return voci


def md5_pubblicati(numero: str, file: Dict[str, dict],
                   apri_url: Callable = apri) -> Dict[str, Tuple[Optional[str], str]]:
    """Gli MD5 che lo studio pubblica, per nome nell'elenco: valore ed elenco
    di provenienza. Un nome con due MD5 diversi vale ``None``: e' ambiguo, e
    per quel file resta la sola verifica della dimensione."""
    pubblicati: Dict[str, Tuple[Optional[str], str]] = {}
    for elenco in elenchi_md5(file):
        with apri_url(indirizzo(numero, file[elenco])) as risposta:
            testo = risposta.read().decode("utf-8", errors="replace")
        for nome, valori in leggi_elenco_md5(testo).items():
            precedenti = [pubblicati[nome][0]] if nome in pubblicati else []
            distinti = set(valori + precedenti)
            pubblicati[nome] = (distinti.pop() if len(distinti) == 1 else None, elenco)
    return pubblicati


def md5_del_file(nome: str, pubblicati: Dict[str, Tuple[Optional[str], str]]
                 ) -> Optional[Tuple[str, str]]:
    """L'MD5 pubblicato di un file dell'archivio e l'elenco che lo dichiara, o
    ``None``.

    Gli elenchi nominano i file senza il prefisso che OSDR aggiunge
    (``GLDS-653_Amplicon_ERX1_raw.fastq.gz`` e' ``ERX1_raw.fastq.gz``): vale il
    nome identico, altrimenti il nome piu' lungo con cui il file finisce dopo
    un trattino basso.
    """
    candidati = [n for n in pubblicati if n == nome or nome.endswith("_" + n)]
    if not candidati:
        return None
    valore, elenco = pubblicati[max(candidati, key=len)]
    return None if valore is None else (valore, elenco)


@functools.lru_cache(maxsize=None)
def _impronte(percorso: str, byte: int, modificato: int, inode: int) -> Tuple[str, str]:
    """MD5 e SHA-256 di un file, letto una volta a blocchi; la chiave comprende
    dimensione, istante di modifica e inode, quindi un file riscritto si
    rilegge."""
    md5, sha256 = hashlib.md5(), hashlib.sha256()
    with open(percorso, "rb") as ingresso:
        for blocco in iter(lambda: ingresso.read(BLOCCO), b""):
            md5.update(blocco)
            sha256.update(blocco)
    return md5.hexdigest(), sha256.hexdigest()


def impronte(percorso: Path) -> Tuple[str, str]:
    """MD5 e SHA-256 del file, calcolati una volta per versione del file."""
    stato = percorso.stat()
    return _impronte(str(percorso), stato.st_size, stato.st_mtime_ns, stato.st_ino)


def completo(percorso: Path, byte: int, md5: Optional[str] = None) -> bool:
    """Vero se il file c'e' con la dimensione dichiarata e, se l'MD5 e'
    pubblicato, con quell'MD5."""
    if not percorso.is_file() or percorso.stat().st_size != byte:
        return False
    return md5 is None or impronte(percorso)[0] == md5


def indirizzo(numero: str, voce: dict) -> str:
    """L'indirizzo di scarico di un file dell'elenco."""
    relativo = voce.get("remote_url")
    if relativo:
        return SITO + relativo
    return SCARICO.format(numero=numero, nome=urllib.parse.quote(voce["file_name"]))


def scarica(url: str, destinazione: Path, byte: int, apri_url: Callable = apri,
            md5: Optional[Tuple[str, str]] = None) -> bool:
    """Porta il file a destinazione; restituisce falso se c'era gia', integro.

    ``md5`` e' l'MD5 pubblicato con l'elenco che lo dichiara, o ``None``. Il
    file gia' presente e' integro se ha la dimensione dichiarata e, con
    ``md5``, quell'MD5. Scaricato un file con un MD5 diverso, lo si cancella
    e lo si scarica di nuovo, da capo: un trasferimento puo' essersi guastato
    senza cambiare la dimensione. Se anche il secondo e' diverso il file non
    si tiene, e l'errore e' :class:`ImprontaDiversa`.
    """
    atteso = None if md5 is None else md5[0]
    if completo(destinazione, byte, atteso):
        return False
    for _ in range(SCARICHI_PER_MD5):
        trasferisci(url, destinazione, byte, apri_url)
        calcolato = None if atteso is None else impronte(destinazione)[0]
        if calcolato == atteso:
            return True
        destinazione.unlink()
    raise ImprontaDiversa(
        f"{destinazione.name}: l'MD5 del file scaricato ({calcolato}) non coincide con "
        f"quello pubblicato da OSDR in {md5[1]} ({atteso}), anche scaricandolo "
        f"{SCARICHI_PER_MD5} volte: il file non e' stato tenuto. Se rilanciando l'esito "
        "non cambia, l'archivio serve un file diverso da quello di cui pubblica l'MD5: "
        "va segnalato a OSDR")


def trasferisci(url: str, destinazione: Path, byte: int, apri_url: Callable = apri) -> None:
    """Trasferisce il file su ``<nome>.parziale``, riprendendo da dove si era
    interrotto, e gli da' il nome definitivo solo con la dimensione dichiarata.
    """
    parziale = destinazione.with_name(destinazione.name + ".parziale")
    ultimo: Optional[BaseException] = None
    for tentativo in range(TENTATIVI):
        if tentativo:
            time.sleep(2 ** tentativo)
        presenti = parziale.stat().st_size if parziale.is_file() else 0
        if presenti > byte:
            parziale.unlink()
            presenti = 0
        try:
            if presenti < byte:
                with apri_url(url, presenti) as risposta:
                    # Senza il 206 il server manda il file intero: si riparte da capo.
                    ripresa = presenti > 0 and getattr(risposta, "status", 200) == 206
                    with open(parziale, "ab" if ripresa else "wb") as uscita:
                        while True:
                            blocco = risposta.read(BLOCCO)
                            if not blocco:
                                break
                            uscita.write(blocco)
        except GUASTI_DI_RETE as guasto:
            ultimo = guasto
            continue
        if parziale.stat().st_size == byte:
            parziale.replace(destinazione)
            return
        ultimo = OSError(f"ricevuti {parziale.stat().st_size} byte su {byte}")
    raise OSError(f"{destinazione.name}: non scaricato dopo {TENTATIVI} tentativi ({ultimo}). "
                  "Rilancia lo stesso comando: riprende da dove si e' fermato")


def estrai(archivio: Path, cartella: Path) -> List[str]:
    """Estrae le tabelle dell'archivio nella cartella, senza sottocartelle."""
    estratti = []
    try:
        with zipfile.ZipFile(archivio) as contenuto:
            for voce in contenuto.infolist():
                # Il solo nome del file: un percorso scritto nell'archivio non
                # porta mai fuori dalla cartella dei metadati.
                nome = Path(voce.filename).name
                if voce.is_dir() or not nome:
                    continue
                (cartella / nome).write_bytes(contenuto.read(voce))
                estratti.append(nome)
    except zipfile.BadZipFile as guasto:
        raise Rifiuto(f"{archivio.name} non si legge come archivio zip ({guasto}): "
                      "cancellalo e rilancia") from guasto
    return sorted(estratti)


def leggi_tabella(percorso: Path) -> Tuple[List[str], List[List[str]]]:
    """Intestazione e righe di una tabella separata da tabulazioni."""
    with open(percorso, encoding="utf-8-sig", newline="") as ingresso:
        righe = [r for r in csv.reader(ingresso, delimiter="\t") if any(c.strip() for c in r)]
    if not righe:
        raise Rifiuto(f"{percorso.name} e' vuota")
    return righe[0], righe[1:]


def tabella_del_16s(cartella: Path, scelta: Optional[str]) -> Path:
    """La tabella di assay del 16S fra quelle estratte.

    Con piu' tabelle candidate non sceglie: le nomina e chiede ``--assay``.
    """
    tabelle = sorted(cartella.glob("a_*.txt"))
    if not tabelle:
        raise Rifiuto("l'archivio dei metadati non contiene tabelle di assay (a_*.txt)")
    if scelta is not None:
        indicate = [t for t in tabelle if t.name == scelta]
        if not indicate:
            raise Rifiuto(f"--assay: {scelta} non e' fra le tabelle di assay: "
                          + ", ".join(t.name for t in tabelle))
        return indicate[0]
    del_16s = [t for t in tabelle if re.search(r"(?i)16s", t.name)]
    if len(del_16s) == 1:
        return del_16s[0]
    candidate = del_16s or tabelle
    if len(candidate) == 1:
        return candidate[0]
    raise Rifiuto(
        ("piu' tabelle di assay hanno «16S» nel nome" if del_16s else
         "nessuna tabella di assay ha «16S» nel nome, e ce n'e' piu' d'una")
        + ": rilancia con --assay NOME, scegliendo quella del sequenziamento 16S fra "
        + ", ".join(t.name for t in candidate))


def nomi_di_una_cella(cella: str) -> List[str]:
    """I nomi di file elencati in una cella, separati da virgole."""
    return [n.strip() for n in cella.split(",") if n.strip()]


#: Il marcatore della lettura nel nome di un file: 1 o 2 dopo un separatore o
#: dopo una R, e prima di un separatore o dell'estensione (``_R1_001``, ``_1``,
#: ``.R2``). Il numero della corsia (``L001``) o del campione non lo e'.
_MARCATORE = re.compile(r"(?:(?<=[_.\-])|(?<=[_.\-]R)|(?<=[_.\-]r))[12](?=[_.\-]|$)")


def forward_di(nomi: Sequence[str]) -> List[str]:
    """Dei file di una riga, quelli della prima lettura.

    Due file che differiscono per il solo marcatore della lettura, 1 nell'uno
    e 2 nell'altro, sono le due letture di una coppia: resta il primo. Un file
    senza compagno resta com'e'.
    """
    def radice(nome: str) -> str:
        """Il nome del file senza l'estensione riconosciuta delle letture."""
        for estensione in ESTENSIONI:
            if nome.lower().endswith(estensione):
                return nome[: -len(estensione)]
        return nome

    scartati = set()
    for a in nomi:
        for trovato in _MARCATORE.finditer(radice(a)):
            if trovato.group() != "1":
                continue
            compagno = a[: trovato.start()] + "2" + a[trovato.end():]
            if compagno in nomi:
                scartati.add(compagno)
    return [n for n in nomi if n not in scartati]


def marcatore_in_fondo(nome: str) -> bool:
    """Vero se il marcatore della lettura chiude il nome del file, prima
    dell'estensione (anche seguito dal numero di blocco, ``_001``, o da
    ``_raw``): e' la sola forma in cui lo strumento abbina i due file di una coppia.
    """
    radice = nome
    for estensione in ESTENSIONI:
        if nome.lower().endswith(estensione):
            radice = nome[: -len(estensione)]
            break
    return re.search(r"[_.][Rr]?[12](?:_[0-9]+|_[Rr][Aa][Ww])?$", radice) is not None


def letture_della_tabella(tabella: Path, solo_forward: bool) -> List[str]:
    """I file di letture che la tabella di assay elenca, nell'ordine delle righe."""
    intestazione, righe = leggi_tabella(tabella)
    colonne = [i for i, nome in enumerate(intestazione) if nome.strip() == COLONNA_DEI_FILE]
    if not colonne:
        raise Rifiuto(f"{tabella.name} non ha la colonna «{COLONNA_DEI_FILE}»: i file delle "
                      "letture non si possono elencare")
    nomi: List[str] = []
    for riga in righe:
        della_riga = [n for i in colonne if i < len(riga) for n in nomi_di_una_cella(riga[i])
                      if n.lower().endswith(ESTENSIONI)]
        nomi += forward_di(della_riga) if solo_forward else della_riga
    nomi = list(dict.fromkeys(nomi))
    if not nomi:
        raise Rifiuto(f"la colonna «{COLONNA_DEI_FILE}» di {tabella.name} non nomina file di "
                      "letture (" + ", ".join(ESTENSIONI) + ")")
    return nomi


def _gb(byte: int) -> str:
    """La dimensione in gigabyte, con due decimali e la virgola."""
    return f"{byte / 1e9:.2f}".replace(".", ",") + " GB"


def scrivi_impronte(cartella: Path, dimensioni: Dict[Path, int],
                    md5: Dict[Path, Optional[Tuple[str, str]]]) -> Path:
    """Scrive ``impronte.tsv`` nella cartella dei dati: per ogni file completo
    (dimensione dichiarata e, se pubblicato, MD5) il percorso relativo alla
    cartella, i byte, l'SHA-256 calcolato e l'MD5 pubblicato con l'elenco che
    lo dichiara (vuoti se l'archivio non lo pubblica). Un file rimasto
    incompleto per un guasto non vi compare. Le righe sono in ordine di
    percorso, e la scrittura passa da un file provvisorio: rilanciare lo
    script riscrive lo stesso elenco.
    """
    righe = ["file\tbyte\tsha256\tmd5_pubblicato\telenco_md5"]
    for percorso in sorted(dimensioni):
        pubblicato = md5.get(percorso)
        if not completo(percorso, dimensioni[percorso], pubblicato[0] if pubblicato else None):
            continue
        righe.append("\t".join((
            percorso.relative_to(cartella).as_posix(), str(percorso.stat().st_size),
            impronte(percorso)[1], pubblicato[0] if pubblicato else "",
            pubblicato[1] if pubblicato else "")))
    elenco = cartella / ELENCO_IMPRONTE
    provvisorio = elenco.with_name(elenco.name + ".parziale")
    provvisorio.write_text("\n".join(righe) + "\n", encoding="utf-8")
    provvisorio.replace(elenco)
    return elenco


def esegui(args: argparse.Namespace, stampa: Callable[[str], None] = print,
           apri_url: Callable = apri) -> int:
    """Scarica metadati e letture, verifica e registra le impronte e stampa il
    comando di ``init``."""
    numero = numero_dello_studio(args.studio)
    cartella: Path = args.cartella
    metadati, letture = cartella / "metadati", cartella / "letture"
    metadati.mkdir(parents=True, exist_ok=True)

    stampa(f"OSD-{numero}: elenco dei file dall'archivio...")
    file = file_dello_studio(numero, apri_url)
    pubblicati = md5_pubblicati(numero, file, apri_url)
    nome_archivio = archivio_dei_metadati(file)
    archivio = metadati / nome_archivio
    md5 = {archivio: md5_del_file(nome_archivio, pubblicati)}
    scarica(indirizzo(numero, file[nome_archivio]), archivio,
            int(file[nome_archivio]["file_size"]), apri_url, md5[archivio])
    estratti = estrai(archivio, metadati)
    stampa(f"Metadati in {metadati}: " + ", ".join(estratti))

    tabella = tabella_del_16s(metadati, args.assay)
    nomi = letture_della_tabella(tabella, args.solo_forward)
    assenti = [n for n in nomi if n not in file]
    if assenti:
        raise Rifiuto(
            f"{len(assenti)} file elencati in {tabella.name} non sono fra quelli pubblicati "
            "dello studio: " + ", ".join(assenti[:5]) + ("..." if len(assenti) > 5 else ""))
    letture.mkdir(exist_ok=True)
    dimensioni = {n: int(file[n]["file_size"]) for n in nomi}
    for n in nomi:
        md5[letture / n] = md5_del_file(n, pubblicati)
    con_md5 = [n for n in nomi if md5[letture / n] is not None]
    presenti = [n for n in nomi if completo(
        letture / n, dimensioni[n], None if md5[letture / n] is None else md5[letture / n][0])]
    da_scaricare = [n for n in nomi if n not in presenti]
    stampa(f"Tabella di assay: {tabella.name}")
    stampa(f"Letture: {len(nomi)} file, {_gb(sum(dimensioni.values()))}"
           + (" (le sole prime letture di ogni coppia)" if args.solo_forward else "")
           + f"; gia' presenti {len(presenti)}, da scaricare {len(da_scaricare)} "
           f"({_gb(sum(dimensioni[n] for n in da_scaricare))})")
    elenchi = sorted({md5[letture / n][1] for n in con_md5})
    stampa(f"Impronte: {len(con_md5)} file su {len(nomi)} con l'MD5 pubblicato da OSDR"
           + (f" ({', '.join(elenchi)})" if elenchi else "")
           + (", verificati con dimensione e MD5" if con_md5 else "")
           + (f"; per {len(nomi) - len(con_md5)} l'archivio non pubblica un MD5, e si "
              "verifica la sola dimensione" if len(con_md5) < len(nomi) else ""))

    falliti: List[str] = []
    diversi: List[str] = []
    with ThreadPoolExecutor(max_workers=args.connessioni) as gruppo:
        in_corso = {
            gruppo.submit(scarica, indirizzo(numero, file[n]), letture / n, dimensioni[n],
                          apri_url, md5[letture / n]): n for n in da_scaricare}
        try:
            for fatti, futuro in enumerate(as_completed(in_corso), 1):
                nome = in_corso[futuro]
                try:
                    futuro.result()
                    stampa(f"  [{fatti}/{len(da_scaricare)}] {nome}")
                except OSError as guasto:
                    (diversi if isinstance(guasto, ImprontaDiversa) else falliti).append(nome)
                    stampa(f"  [{fatti}/{len(da_scaricare)}] NON SCARICATO: {guasto}")
        except KeyboardInterrupt:
            for futuro in in_corso:
                futuro.cancel()
            stampa("Interrotto: rilancia lo stesso comando per riprendere.")
            raise
    elenco = scrivi_impronte(
        cartella, {archivio: int(file[nome_archivio]["file_size"]),
                   **{letture / n: dimensioni[n] for n in nomi}}, md5)
    stampa(f"SHA-256 di ogni file scaricato e verificato in {elenco}.")
    if falliti:
        stampa(f"{len(falliti)} file non scaricati: rilancia lo stesso comando, che riprende "
               "da dove si e' fermato.")
    if diversi:
        stampa(f"{len(diversi)} file con un MD5 diverso da quello pubblicato da OSDR, anche "
               "scaricati due volte, e non tenuti: " + ", ".join(diversi[:5])
               + ("..." if len(diversi) > 5 else ""))
    if falliti or diversi:
        return USCITA_GUASTO

    appaiate = not args.solo_forward and len(forward_di(nomi)) < len(nomi)
    stampa(f"Fatto. In {letture}: {len(nomi)} file di letture"
           + (", due per campione (letture appaiate)" if appaiate else "") + ".")
    if appaiate and not all(marcatore_in_fondo(n) for n in nomi):
        # Lo strumento abbina i due file di una coppia solo con il marcatore in
        # fondo al nome: va detto prima che init rifiuti la cartella.
        stampa("ATTENZIONE: nei nomi di questi file il marcatore della lettura (R1, R2) non "
               "e' in fondo al nome, e lo strumento non li riconosce come coppie: init "
               "rifiuterebbe questa cartella. Usa le sole prime letture: rilancia con "
               "--solo-forward in un'altra cartella, e lancia init su quella.")
        return 0
    if appaiate:
        stampa("Per usare le sole prime letture: rilancia con --solo-forward in un'altra "
               "cartella.")
    stampa("Il comando di init, con la cartella montata nel container come /dati:")
    stampa(f"  p16s init --reads /dati/letture --metadata /dati/metadati --assay {tabella.name} \\")
    stampa("    --reference /app/dati/riferimento/silva_nr99_v138_train_set.fa.gz \\")
    stampa("    --out-root /lavoro/esecuzione --output /lavoro/configurazione.yaml "
           "--accetta-proposte")
    return 0


def build_parser() -> argparse.ArgumentParser:
    """Gli argomenti della riga di comando: lo studio, la cartella dei dati e
    le opzioni sulla tabella di assay, sulle sole forward e sulle connessioni."""
    parser = argparse.ArgumentParser(
        description="Scarica da OSDR i metadati ISA-Tab e le letture 16S di uno studio.")
    parser.add_argument("studio", metavar="OSD-NNN", help="l'identificativo dello studio")
    parser.add_argument("cartella", type=Path, metavar="CARTELLA",
                        help="la cartella dei dati: vi si creano metadati/ e letture/")
    parser.add_argument("--assay", metavar="NOME",
                        help="il nome della tabella di assay del 16S, se lo studio ne ha "
                             "piu' d'una candidata")
    parser.add_argument("--solo-forward", action="store_true",
                        help="di ogni coppia scarica il solo file della prima lettura")
    parser.add_argument("--connessioni", type=int, default=CONNESSIONI, metavar="N",
                        help=f"connessioni aperte insieme (predefinito {CONNESSIONI})")
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    """Punto d'ingresso: esegue lo scarico e traduce gli esiti in codici di
    uscita.

    Un argomento non valido o un rifiuto danno il codice del rifiuto; un
    guasto di rete quello del guasto, con l'invito a rilanciare (lo scarico
    riprende), e cosi' un MD5 diverso da quello pubblicato, con il suo
    messaggio; un'interruzione da tastiera 130, senza messaggio.
    """
    args = build_parser().parse_args(argv)
    if args.connessioni < 1:
        print("scarica_osdr: --connessioni deve essere almeno 1", file=sys.stderr)
        return USCITA_RIFIUTO
    try:
        return esegui(args)
    except Rifiuto as rifiuto:
        print(f"scarica_osdr: {rifiuto}", file=sys.stderr)
        return USCITA_RIFIUTO
    except ImprontaDiversa as diversa:
        print(f"scarica_osdr: {diversa}", file=sys.stderr)
        return USCITA_GUASTO
    except GUASTI_DI_RETE as guasto:
        print(f"scarica_osdr: l'archivio non risponde ({guasto}). Controlla la rete e "
              "rilancia lo stesso comando: riprende da dove si e' fermato", file=sys.stderr)
        return USCITA_GUASTO
    except KeyboardInterrupt:
        return 130


if __name__ == "__main__":
    raise SystemExit(main())
