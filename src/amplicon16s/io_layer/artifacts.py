"""Albero degli artefatti di output e manifesti dei checksum.

Sotto ``io.out_root`` la pipeline scrive una cartella per fase, numerata in
modo che l'ordine sul filesystem coincida con l'ordine di esecuzione.

Ogni artefatto scritto viene registrato nel **manifesto della sua cartella**,
con il checksum e la dimensione. Quel manifesto risponde alla domanda «i file
di questa cartella sono ancora quelli scritti?», e ogni scrittura lo aggiorna.

**Il completamento di una fase ha un manifesto proprio.** Le cartelle sono
quattordici e le fasi quindici, e alcune fasi condividono una cartella: S11 e
S12 scrivono in ``11_controls``, S13 e S14 in ``12_final``. Con il solo
manifesto di cartella, la conclusione di S11 renderebbe «completa» la cartella
e con essa S12, che non ha mai girato. Ogni fase conclusa scrive quindi il
proprio :class:`ManifestoPasso` (``manifest_S11.json``, ``manifest_S12.json``)
con i propri artefatti e ciò su cui è stata calcolata. Lo scrive per ultimo e
in modo atomico: una fase interrotta a metà non ne ha uno, e una fase che
riscrive un file di un'altra ne rompe il checksum invece di sembrare conclusa.

Gli artefatti non nascono tutti qui: le fasi di calcolo sono processi R che
scrivono i propri file da sé. Per questo :meth:`AlberoOutput.registra` esiste
accanto ai metodi di scrittura — un file prodotto altrove entra nel manifesto
allo stesso modo di uno scritto da Python, e la verifica di completezza non
distingue i due casi.
"""

from __future__ import annotations

import hashlib
import json
import os
import tempfile
import uuid
from collections.abc import Callable, Iterable, Iterator, Mapping
from dataclasses import dataclass
from datetime import datetime, timezone
from enum import StrEnum
from pathlib import Path
from typing import Any, Final

from amplicon16s.io_layer.checksums import checksum_file, corrisponde

__all__ = [
    "AlberoOutput",
    "Artefatto",
    "Fase",
    "ManifestoPasso",
    "ManifestoPassoNonValido",
    "NOME_MANIFESTO",
    "nome_manifesto_passo",
]

#: Nome del manifesto dentro ogni cartella di fase.
NOME_MANIFESTO: Final = "manifest.json"


def nome_manifesto_passo(passo: str) -> str:
    """Nome del manifesto di completamento di una fase, nella sua cartella."""
    return f"manifest_{passo}.json"


class Fase(StrEnum):
    """Cartelle di output, una per fase.

    Il valore è il nome della cartella: la numerazione fa sì che l'ordine
    alfabetico sul filesystem sia l'ordine di esecuzione.
    """

    CONFIG = "00_config"
    INPUT_VALIDATION = "01_input_validation"
    QC_PROFILES = "02_qc_profiles"
    FILTERED = "03_filtered"
    ERROR_MODELS = "04_error_models"
    ASV_INFERENCE = "05_asv_inference"
    SEQTAB = "06_seqtab"
    CHIMERA = "07_chimera"
    TAXONOMY = "08_taxonomy"
    PHYLOGENY = "09_phylogeny"
    PHYLOSEQ = "10_phyloseq"
    CONTROLS = "11_controls"
    FINAL = "12_final"
    LOGS = "99_logs"


@dataclass(frozen=True)
class Artefatto:
    """Un file prodotto dalla pipeline, con la sua impronta."""

    fase: Fase
    nome: str
    percorso: Path
    checksum: str
    byte: int

    def come_voce(self) -> dict[str, object]:
        """Forma registrata nel manifesto."""
        return {"nome": self.nome, "checksum": self.checksum, "byte": self.byte}

    @property
    def integro(self) -> bool:
        """Se il file su disco corrisponde ancora al checksum registrato."""
        from amplicon16s.io_layer.checksums import corrisponde

        return corrisponde(self.percorso, self.checksum)


def _canonico(documento: Mapping[str, Any]) -> str:
    return json.dumps(
        documento, sort_keys=True, ensure_ascii=False, separators=(",", ":"),
        allow_nan=False,
    )


def _impronta(documento: Mapping[str, Any]) -> str:
    return "sha256:" + hashlib.sha256(_canonico(documento).encode("utf-8")).hexdigest()


class ManifestoPassoNonValido(ValueError):
    """Il manifesto di una fase esiste ma è illeggibile o è stato alterato."""


@dataclass(frozen=True)
class ManifestoPasso:
    """La traccia che una fase lascia quando si conclude.

    ``calcolata_su`` dice su che cosa la fase è stata calcolata: la
    configurazione, le fasi a monte, i dati esterni. ``impronta`` identifica
    questo calcolo ed è ciò che le fasi a valle registrano come proprio
    ingresso: ricalcolare la fase ne produce una nuova, e chi aveva registrato
    la vecchia non risulta piu' calcolato sugli ingressi di adesso.
    """

    passo: str
    fase: Fase
    impronta: str
    calcolata_su: dict[str, Any]
    artefatti: tuple[dict[str, Any], ...]
    metriche: dict[str, Any]
    conclusa: str
    esecuzione: str
    #: Parametri cambiati da un'azione correttiva: codice che l'ha causata,
    #: parametro, valore dichiarato e valore usato.
    aggiustamenti: tuple[dict[str, Any], ...] = ()
    #: Degradazioni registrate: la fase si e' conclusa con un ripiego.
    degradazioni: tuple[dict[str, Any], ...] = ()

    def _contenuto(self) -> dict[str, Any]:
        return {
            "passo": self.passo,
            "cartella": self.fase.value,
            "calcolata_su": self.calcolata_su,
            "artefatti": list(self.artefatti),
            "metriche": self.metriche,
            "aggiustamenti": list(self.aggiustamenti),
            "degradazioni": list(self.degradazioni),
            "conclusa": self.conclusa,
            # Identifica l'esecuzione: due calcoli della stessa fase non hanno
            # mai la stessa impronta, anche se producono gli stessi byte.
            "esecuzione": self.esecuzione,
        }

    def come_documento(self) -> dict[str, Any]:
        return {"impronta": self.impronta, **self._contenuto()}

    @property
    def nomi(self) -> tuple[str, ...]:
        return tuple(str(v["nome"]) for v in self.artefatti)


class AlberoOutput:
    """L'albero delle cartelle di output sotto una radice.

    Creare l'oggetto non tocca il filesystem: le cartelle nascono con
    :meth:`crea`, oppure da sé alla prima scrittura in una fase.
    """

    def __init__(self, radice: Path | str) -> None:
        self.radice = Path(radice)

    # ----------------------------------------------------------------- #
    # Cartelle                                                           #
    # ----------------------------------------------------------------- #

    def cartella(self, fase: Fase) -> Path:
        """Percorso della cartella di una fase, senza crearla."""
        return self.radice / fase.value

    def prepara(self, fase: Fase) -> Path:
        """Percorso della cartella di una fase, creandola se manca."""
        percorso = self.cartella(fase)
        percorso.mkdir(parents=True, exist_ok=True)
        return percorso

    def crea(self) -> tuple[Path, ...]:
        """Crea l'albero completo. Ripetibile senza effetti."""
        return tuple(self.prepara(fase) for fase in Fase)

    # ----------------------------------------------------------------- #
    # Scrittura                                                          #
    # ----------------------------------------------------------------- #

    def scrivi_testo(self, fase: Fase, nome: str, contenuto: str) -> Artefatto:
        """Scrive un artefatto di testo e lo registra nel manifesto."""
        percorso = self.prepara(fase) / nome
        percorso.write_text(contenuto, encoding="utf-8")
        return self.registra(fase, nome)

    def scrivi_bytes(self, fase: Fase, nome: str, contenuto: bytes) -> Artefatto:
        """Scrive un artefatto binario e lo registra nel manifesto."""
        percorso = self.prepara(fase) / nome
        percorso.write_bytes(contenuto)
        return self.registra(fase, nome)

    def registra(self, fase: Fase, nome: str) -> Artefatto:
        """Registra nel manifesto un file già presente nella cartella di fase.

        È la via per gli artefatti prodotti dai processi R: il file esiste
        già, e qui se ne calcola l'impronta.
        """
        percorso = self.cartella(fase) / nome
        if not percorso.is_file():
            raise FileNotFoundError(
                f"artefatto non trovato, non si puo' registrare: {percorso}"
            )

        artefatto = Artefatto(
            fase=fase,
            nome=nome,
            percorso=percorso,
            checksum=checksum_file(percorso),
            byte=percorso.stat().st_size,
        )
        self._aggiorna_manifesto(artefatto)
        return artefatto

    # ----------------------------------------------------------------- #
    # Manifesto                                                          #
    # ----------------------------------------------------------------- #

    def percorso_manifesto(self, fase: Fase) -> Path:
        return self.cartella(fase) / NOME_MANIFESTO

    def manifesto(self, fase: Fase) -> dict[str, dict[str, object]]:
        """Contenuto del manifesto di una fase, vuoto se non esiste."""
        percorso = self.percorso_manifesto(fase)
        if not percorso.is_file():
            return {}
        documento = json.loads(percorso.read_text(encoding="utf-8"))
        return {v["nome"]: v for v in documento.get("artefatti", [])}

    def _aggiorna_manifesto(self, artefatto: Artefatto) -> None:
        voci = self.manifesto(artefatto.fase)
        voci[artefatto.nome] = artefatto.come_voce()

        documento = {
            "fase": artefatto.fase.value,
            # Ordinate per nome: il manifesto non deve cambiare solo perche' e'
            # cambiato l'ordine di scrittura.
            "artefatti": [voci[n] for n in sorted(voci)],
        }
        self.percorso_manifesto(artefatto.fase).write_text(
            json.dumps(documento, indent=2, ensure_ascii=False, sort_keys=False) + "\n",
            encoding="utf-8",
        )

    # ----------------------------------------------------------------- #
    # Verifica                                                           #
    # ----------------------------------------------------------------- #

    def artefatti(self, fase: Fase) -> Iterator[Artefatto]:
        """Gli artefatti registrati per una fase, come oggetti."""
        for nome, voce in self.manifesto(fase).items():
            yield Artefatto(
                fase=fase,
                nome=nome,
                percorso=self.cartella(fase) / nome,
                checksum=str(voce["checksum"]),
                byte=int(voce["byte"]),
            )

    def non_integri(self, fase: Fase) -> tuple[str, ...]:
        """Nomi degli artefatti registrati che mancano o sono stati alterati."""
        return tuple(a.nome for a in self.artefatti(fase) if not a.integro)

    def fase_completa(self, fase: Fase, attesi: Iterable[str] | None = None) -> bool:
        """Se la cartella contiene artefatti registrati e sono ancora integri.

        Senza ``attesi`` la domanda è «tutto ciò che risulta prodotto è ancora
        integro?»; con ``attesi`` è «ci sono anche questi?».

        Riguarda la cartella, non una fase del grafo: con due fasi nella
        stessa cartella non distingue l'una dall'altra, e non sa su che cosa i
        file sono stati calcolati. Il completamento di una fase si stabilisce
        con :meth:`manifesto_passo` e :meth:`non_integri_del_passo`.
        """
        registrati = self.manifesto(fase)
        if not registrati:
            return False
        if attesi is not None and not set(attesi) <= set(registrati):
            return False
        return not self.non_integri(fase)

    # ----------------------------------------------------------------- #
    # Manifesto di una fase                                              #
    # ----------------------------------------------------------------- #

    def percorso_manifesto_passo(self, passo: str, fase: Fase) -> Path:
        return self.cartella(fase) / nome_manifesto_passo(passo)

    def rimuovi_manifesto_passo(self, passo: str, fase: Fase) -> None:
        """Una fase che viene ricalcolata smette subito di risultare conclusa."""
        self.percorso_manifesto_passo(passo, fase).unlink(missing_ok=True)

    def concludi_passo(
        self,
        passo: str,
        fase: Fase,
        artefatti: Iterable[Artefatto],
        calcolata_su: Mapping[str, Any],
        metriche: Mapping[str, Any] | None = None,
        aggiustamenti: Iterable[Mapping[str, Any]] = (),
        degradazioni: Iterable[Mapping[str, Any]] = (),
    ) -> ManifestoPasso:
        """Scrive il manifesto di una fase conclusa, per ultimo e atomicamente."""
        voci = []
        for artefatto in artefatti:
            if artefatto.fase is not fase:
                raise ValueError(
                    f"{passo} scrive in {fase.value}, non in {artefatto.fase.value}: "
                    f"{artefatto.nome}"
                )
            voci.append(artefatto.come_voce())
        voci.sort(key=lambda v: str(v["nome"]))

        bozza = ManifestoPasso(
            passo=passo,
            fase=fase,
            impronta="",
            calcolata_su=dict(calcolata_su),
            artefatti=tuple(voci),
            metriche=dict(metriche or {}),
            conclusa=datetime.now(timezone.utc).isoformat(timespec="seconds"),
            esecuzione=uuid.uuid4().hex,
            aggiustamenti=tuple(dict(a) for a in aggiustamenti),
            degradazioni=tuple(dict(d) for d in degradazioni),
        )
        # Il passaggio per JSON e ritorno fissa la forma che verra' riletta:
        # l'impronta si calcola su quella, non sugli oggetti in memoria.
        contenuto = json.loads(_canonico(bozza._contenuto()))
        manifesto = ManifestoPasso(
            passo=passo,
            fase=fase,
            impronta=_impronta(contenuto),
            calcolata_su=contenuto["calcolata_su"],
            artefatti=tuple(contenuto["artefatti"]),
            metriche=contenuto["metriche"],
            conclusa=contenuto["conclusa"],
            esecuzione=contenuto["esecuzione"],
            aggiustamenti=tuple(contenuto["aggiustamenti"]),
            degradazioni=tuple(contenuto["degradazioni"]),
        )

        self.prepara(fase)
        percorso = self.percorso_manifesto_passo(passo, fase)
        descrittore, temporaneo = tempfile.mkstemp(
            prefix=".scrittura-", suffix=".json", dir=percorso.parent
        )
        try:
            with os.fdopen(descrittore, "w", encoding="utf-8") as file:
                json.dump(
                    manifesto.come_documento(), file, indent=2, ensure_ascii=False
                )
                file.write("\n")
            os.replace(temporaneo, percorso)
        except BaseException:
            Path(temporaneo).unlink(missing_ok=True)
            raise
        return manifesto

    def manifesto_passo(self, passo: str, fase: Fase) -> ManifestoPasso | None:
        """Il manifesto di una fase, o ``None`` se la fase non si è conclusa.

        Solleva :class:`ManifestoPassoNonValido` se il file esiste ma non è
        leggibile, appartiene a un'altra fase o è stato modificato dopo la
        scrittura: l'impronta non corrisponderebbe al contenuto.
        """
        percorso = self.percorso_manifesto_passo(passo, fase)
        if not percorso.is_file():
            return None
        try:
            documento = json.loads(percorso.read_text(encoding="utf-8"))
            impronta = documento.pop("impronta")
            manifesto = ManifestoPasso(
                passo=documento["passo"],
                fase=Fase(documento["cartella"]),
                impronta=impronta,
                calcolata_su=documento["calcolata_su"],
                artefatti=tuple(documento["artefatti"]),
                metriche=documento["metriche"],
                conclusa=documento["conclusa"],
                esecuzione=documento["esecuzione"],
                aggiustamenti=tuple(documento.get("aggiustamenti", ())),
                degradazioni=tuple(documento.get("degradazioni", ())),
            )
        except (OSError, ValueError, KeyError, TypeError, AttributeError) as e:
            raise ManifestoPassoNonValido(f"manifesto di {passo} illeggibile: {e}") from e

        if manifesto.passo != passo or manifesto.fase is not fase:
            raise ManifestoPassoNonValido(
                f"il manifesto in {percorso} appartiene a {manifesto.passo}"
            )
        if _impronta(documento) != impronta:
            raise ManifestoPassoNonValido(f"il manifesto di {passo} e' stato alterato")
        return manifesto

    def non_integri_del_passo(
        self,
        manifesto: ManifestoPasso,
        checksum: Callable[[Path], str | None] | None = None,
    ) -> tuple[str, ...]:
        """Artefatti della fase che mancano o non corrispondono al checksum.

        ``checksum`` restituisce il checksum di un file, o ``None`` se il file
        manca; permette a chi valuta piu' fasi di calcolare ciascun checksum
        una volta sola. Senza, il checksum si calcola qui.
        """
        cartella = self.cartella(manifesto.fase)
        if checksum is None:
            return tuple(
                str(v["nome"])
                for v in manifesto.artefatti
                if not corrisponde(cartella / str(v["nome"]), str(v["checksum"]))
            )
        return tuple(
            str(v["nome"])
            for v in manifesto.artefatti
            if checksum(cartella / str(v["nome"])) != str(v["checksum"])
        )
