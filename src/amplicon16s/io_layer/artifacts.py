"""Albero degli artefatti di output e manifesti dei checksum.

Sotto ``io.out_root`` la pipeline scrive una cartella per fase, numerata in
modo che l'ordine sul filesystem coincida con l'ordine di esecuzione.

Ogni artefatto scritto viene registrato nel **manifesto della sua fase**, con
il checksum e la dimensione. Il manifesto è per fase e non unico perché il suo
uso principale è la domanda «questa fase è già stata completata?»: con un
manifesto per cartella la risposta si ottiene leggendo un file solo, e una
fase che non è mai partita non ha manifesto da distinguere da una che è andata
a metà.

Gli artefatti non nascono tutti qui: le fasi di calcolo sono processi R che
scrivono i propri file da sé. Per questo :meth:`AlberoOutput.registra` esiste
accanto ai metodi di scrittura — un file prodotto altrove entra nel manifesto
allo stesso modo di uno scritto da Python, e la verifica di completezza non
distingue i due casi.
"""

from __future__ import annotations

import json
from collections.abc import Iterable, Iterator
from dataclasses import dataclass
from enum import StrEnum
from pathlib import Path
from typing import Final

from amplicon16s.io_layer.checksums import checksum_file

__all__ = [
    "AlberoOutput",
    "Artefatto",
    "Fase",
    "NOME_MANIFESTO",
]

#: Nome del manifesto dentro ogni cartella di fase.
NOME_MANIFESTO: Final = "manifest.json"


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
        """Se la fase ha prodotto i suoi artefatti e sono ancora integri.

        Senza ``attesi`` la domanda è «tutto ciò che risulta prodotto è ancora
        integro?»; con ``attesi`` è «ci sono anche questi?». La distinzione
        serve alla ripresa di un'esecuzione interrotta, dove un manifesto
        parziale non significa fase conclusa.
        """
        registrati = self.manifesto(fase)
        if not registrati:
            return False
        if attesi is not None and not set(attesi) <= set(registrati):
            return False
        return not self.non_integri(fase)
