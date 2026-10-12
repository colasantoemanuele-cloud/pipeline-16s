"""Fase S0: validazione iniziale.

È la barriera che gira prima di allocare qualunque calcolo costoso. Il suo
scopo è far emergere in pochi minuti un problema che altrimenti si scoprirebbe
dopo ore di elaborazione: un file corrotto, una colonna sbagliata, un
troncamento incompatibile con le lunghezze reali.

**Nessun gate legge un file intero.** Le verifiche sulle sequenze si fermano
alle prime ``qc.head_reads`` letture; l'integrità completa dei file è già
garantita dai checksum. È il vincolo che tiene la validazione nell'ordine dei
minuti.

La fase scrive in ``01_input_validation/`` l'esito di ogni gate, il crosswalk
e l'inventario, passando dal servizio di scrittura degli artefatti: ogni file
prodotto entra nel manifesto con il proprio checksum, come qualunque altro
artefatto della pipeline.

È realizzata come :class:`ValidazioneIngressi`, sulla classe base comune a
tutte le fasi; :func:`esegui_s0` resta il modo di invocarla da sola, con la
stessa interfaccia di prima. S0 è l'unica fase che legge dati fuori
dall'albero di output, e ne registra un'impronta: letture, tabelle e
riferimento sostituiti sotto lo stesso percorso la rendono non piu' conclusa,
anche a configurazione invariata.
"""

from __future__ import annotations

import csv
import hashlib
import json
import time
from dataclasses import dataclass, replace
from io import StringIO
from pathlib import Path
from typing import Any, ClassVar, Final

from amplicon16s.config.resolve import risolvi
from amplicon16s.config.schema import Config
from amplicon16s.errors.catalog import Categoria, voce
from amplicon16s.errors.exceptions import ErrorePipeline
from amplicon16s.gates.g01_g15 import Contesto, ErroreGate, Violazione
from amplicon16s.gates.registry import EsitoGate, esegui_gate, esegui_tutti
from amplicon16s.io_layer.artifacts import AlberoOutput, Fase
from amplicon16s.logging.logger import ottieni
from amplicon16s.metadata.lettura_inventario import NOME_CROSSWALK
from amplicon16s.metadata.models import Inventario
from amplicon16s.rbridge.runner import cartella_r, esegui_script
from amplicon16s.runner.graph import Passo
from amplicon16s.steps.base import PipelineStep, Produzione, StepContext, StepResult

__all__ = [
    "RisultatoS0",
    "ValidazioneIngressi",
    "esegui_s0",
    "impronta_dati_grezzi",
]

#: Lo script che verifica il classificatore addestrato di IdTaxa.
NOME_SCRIPT_CLASSIFICATORE: Final = "00_classificatore.R"

NOME_ESITI: Final = "gates.json"
NOME_INVENTARIO: Final = "inventario.json"
NOME_LETTURE: Final = "letture_ispezionate.tsv"
#: In modalita' appaiata, le stesse statistiche per i file reverse.
NOME_LETTURE_REVERSE: Final = "letture_ispezionate_reverse.tsv"


@dataclass(frozen=True)
class RisultatoS0:
    """Esito della fase, con gli artefatti prodotti."""

    esiti: tuple[EsitoGate, ...]
    inventario: Inventario | None
    artefatti: tuple[Path, ...]
    secondi: float

    @property
    def superata(self) -> bool:
        """Vero se tutti i gate sono superati."""
        return all(e.superato for e in self.esiti)

    @property
    def falliti(self) -> tuple[EsitoGate, ...]:
        """I gate eseguiti e non superati."""
        return tuple(e for e in self.esiti if e.eseguito and not e.superato)

    @property
    def avvisi(self) -> tuple[Any, ...]:
        """Gli avvisi di tutti i gate, nell'ordine di esecuzione."""
        return tuple(a for e in self.esiti for a in e.avvisi)


def _crosswalk_tsv(inventario: Inventario, appaiata: bool = False) -> str:
    """Il crosswalk in forma tabellare: una riga per campione. In modalita'
    appaiata ha in fondo la colonna del file reverse.
    """
    buffer = StringIO()
    scrittore = csv.writer(buffer, delimiter="\t", lineterminator="\n")
    scrittore.writerow(
        ["accession", "sample_name", "classe", "materiale", "file",
         "posizione", "modulo", "piastra", "corsa"]
        + (["file_reverse"] if appaiata else [])
    )
    for campione in inventario:
        scrittore.writerow(
            [campione.accession, campione.nome, campione.classe.value,
             campione.materiale, campione.file.name if campione.file else "",
             campione.posizione or "", campione.modulo or "",
             campione.piastra or "", campione.corsa or ""]
            + ([campione.file_reverse.name if campione.file_reverse else ""]
               if appaiata else [])
        )
    return buffer.getvalue()


def _inventario_json(inventario: Inventario) -> str:
    """L'inventario come documento JSON: conteggi per classe, denominatore di
    prevalenza, moduli, piastre e corse.
    """
    conteggi = inventario.conteggi()
    documento = {
        "campioni": len(inventario),
        "conteggi": {classe.value: numero for classe, numero in conteggi.items()},
        # Il denominatore e' esplicito perche' a valle non lo si confonda con
        # il totale dei campioni: sono entrambi plausibili.
        "denominatore_prevalenza": inventario.denominatore_prevalenza(),
        "moduli": list(inventario.moduli),
        "piastre": list(inventario.piastre),
        "corse": list(inventario.corse),
        "campioni_per_piastra": inventario.campioni_per_piastra(),
        "senza_lotto": len(inventario.senza_lotto),
    }
    return json.dumps(documento, indent=2, ensure_ascii=False) + "\n"


def _letture_tsv(
    scansione: dict[str, Any], head_reads: int, con_motivo: bool = True, con_primer: bool = True
) -> str:
    """Le statistiche per file, con il limite della stima messo agli atti.

    Le lunghezze sono quelle delle prime ``head_reads`` letture, non di tutto
    il file: il minimo puo' quindi essere piu' alto del vero. La colonna
    ``file_esaurito`` dice per quali file la statistica copre tutto, e per
    quali e' una stima. Con ``con_motivo`` falso il motivo conservato non e'
    stato cercato (``qc.conserved_motif`` nullo) e la sua colonna resta vuota:
    uno zero direbbe che nessuna lettura lo porta. Lo stesso per il primer con
    ``con_primer`` falso (``qc.primer_sequence`` nullo).
    """
    buffer = StringIO()
    scrittore = csv.writer(buffer, delimiter="\t", lineterminator="\n")
    scrittore.writerow(
        ["accession", "file", "letture_esaminate", "file_esaurito",
         "lunghezza_minima", "lunghezza_massima", "frazione_primer",
         "frazione_motivo", "errore"]
    )
    for accession in sorted(scansione):
        s = scansione[accession]
        scrittore.writerow(
            [accession, s.nome, s.letture_esaminate, "si" if s.esaurito else "no",
             s.lunghezza_minima if s.lunghezza_minima is not None else "",
             s.lunghezza_massima if s.lunghezza_massima is not None else "",
             f"{s.frazione_primer:.6f}" if con_primer else "",
             f"{s.frazione_motivo:.6f}" if con_motivo else "",
             s.errore or ""]
        )
    return buffer.getvalue()


def _esiti_json(esiti: tuple[EsitoGate, ...]) -> str:
    """L'esito di ogni gate come documento JSON.

    Le durate non ci sono: cambiano a ogni esecuzione, e l'artefatto deve
    avere lo stesso checksum fra due esecuzioni sugli stessi ingressi. Vanno
    nel log (:meth:`ValidazioneIngressi.calcola`).
    """
    documento = {
        "superata": all(e.superato for e in esiti),
        "gate": [e.come_voce() for e in esiti],
    }
    return json.dumps(documento, indent=2, ensure_ascii=False) + "\n"


def impronta_dati_grezzi(config: Config) -> str:
    """Impronta dei dati letti da S0: letture e tabelle.

    Si basa su nome, dimensione e istante di modifica di ogni file, non sul
    contenuto: rileggere gigabyte di letture a ogni ripresa costerebbe piu' della
    validazione stessa. Un file sostituito sotto lo stesso nome cambia quasi
    sempre dimensione o istante di modifica.

    Il riferimento tassonomico non ne fa parte: non e' un ingresso dei
    risultati di S0, e cambiarlo deve rifare S8 e le fasi successive, non la
    catena intera. Il suo contenuto entra nell'impronta di S8 con tax.ref_md5,
    e G12 lo verifica contro quel checksum a ogni avvio.
    """

    def descrivi(percorso: Path) -> list[Any]:
        try:
            stato = percorso.stat()
        except OSError:
            return [str(percorso), "assente"]
        return [str(percorso), stato.st_size, stato.st_mtime_ns]

    io = config.io
    file: list[Path] = sorted(Path(io.fastq_dir).glob(io.fastq_glob))
    file.append(Path(io.assay_table))
    if io.study_table is not None:
        file.append(Path(io.study_table))
    if io.batch_table is not None:
        file.append(Path(io.batch_table))

    canonico = json.dumps([descrivi(f) for f in file], separators=(",", ":"))
    return "sha256:" + hashlib.sha256(canonico.encode("utf-8")).hexdigest()


class ValidazioneIngressi(PipelineStep):
    """S0: i quindici gate, e l'esito in ``01_input_validation/``.

    Con ``solleva`` falso restituisce il risultato anche quando un gate
    fallisce, invece di sollevare: la fase risulta allora non superata e non
    viene registrata come conclusa.
    """

    passo: ClassVar[Passo] = Passo.S0
    #: 6: G10 cerca il primer nelle sole classi controllate; G08 dichiara le
    #: righe del lotto senza campione (E-S0-19); il modulo viene dal lotto
    #: anche senza meta.module_column. 7: G09 giudica ogni classe controllata
    #: da sola; G07 riconosce i file con le due letture di ogni coppia; G02 e
    #: G08 respingono i valori con tabulazioni o a capo. 8: G15, i cui esiti
    #: sono in gates.json, respinge anche meta.non_surface_positions senza la
    #: colonna della posizione. 9: G07 giudica i segni di coppia sulla
    #: frazione delle letture esaminate, sull'intera intestazione, e respinge un
    #: file di sole seconde letture; G02 guarda nella tabella di studio le sole
    #: righe dei campioni dell'assay. 10: un nome ripetuto piu' di due volte
    #: non e' un segno di coppia, e servono almeno due letture con il segno.
    #: 11: i nomi oltre le due copie si tollerano entro la stessa frazione degli
    #: altri segni, e le sole seconde letture tollerano almeno un record.
    #: 12: le righe vuote in fondo a un file di letture sono tollerate; un file
    #: del lotto non UTF-8 ferma G08 (E-S0-08) invece di un errore imprevisto;
    #: G02 riconosce la riga di studio con una tabulazione nell'identificativo.
    #: 13: G07 riconosce le due letture di una coppia anche dal marcatore in
    #: fondo all'identificativo, in letture consecutive; G02 respinge le righe
    #: con piu' valori dell'intestazione e ogni identificativo con una
    #: tabulazione o un a capo. 14: ctrl.column puo' essere nullo (nessun
    #: controllo: ogni campione e' biologico); G13 ammette il file senza
    #: letture di un controllo negativo e lo dichiara (E-S0-20). 15: con il
    #: motivo conservato nullo la colonna frazione_motivo e' vuota. 16: il
    #: primer puo' essere nullo (non si cerca, E-S0-21); una riga finale di
    #: soli spazi in un file di letture e' un difetto di formato (E-S0-13).
    #: 17: due file dello stesso campione che sono le due letture di una coppia
    #: hanno un codice proprio (E-S0-24); un file con il nome di una seconda
    #: lettura non e' respinto per il nome se S1 ne puo' giudicare il contenuto.
    versione: ClassVar[int] = 17
    #: I parametri con cui i gate producono i risultati di S0: l'inventario, il
    #: crosswalk, la scansione delle letture, gli esiti e le degradazioni.
    #: Ingressi e metadati per intero (io, meta); di ctrl la colonna, le tre
    #: etichette della classificazione e la regola di riclassificazione in
    #: controllo negativo, non le chiavi di S11; di qc le cinque chiavi
    #: di G10 e della scansione; di decontam la colonna e il minimo di bianchi
    #: per piastra (G08, E-S0-15); la colonna della corsa (crosswalk, G08); il
    #: troncamento (G09) e il taglio iniziale, da cui dipende dove G10 cerca
    #: primer e motivo; le colonne da portare nell'oggetto e quella delle
    #: cellule dei controlli positivi, che G02 e G08 verificano sulle
    #: intestazioni. L'elenco e' stato ricavato registrando i parametri letti
    #: da un'esecuzione vera, e confrontato con il codice dei gate.
    #:
    #: G15 e G12 sono precondizioni, come G14: non contribuiscono ai risultati
    #: di S0 (se falliscono S0 non produce nulla) e l'esecutore li ripete a ogni
    #: avvio, prima di qualunque fase. G15 legge l'intera configurazione per
    #: verificarne la coerenza; G12 verifica il riferimento tassonomico contro
    #: tax.ref_md5. I parametri che servono solo a loro (asv, prev, le altre
    #: chiavi di decontam e di qc, retry.whitelist, tax)
    #: restano quindi fuori dall'impronta di S0: cambiare una soglia di S6 o il
    #: riferimento tassonomico non la rende da rifare, e con lei la catena
    #: intera; il riferimento entra nell'impronta di S8. G15 e G12 si eseguono
    #: in :meth:`esegui`, sulla configurazione completa, prima che la fase la
    #: veda ristretta, e figurano in gates.json fra gli altri gate.
    parametri: ClassVar[tuple[str, ...]] = (
        "io", "meta",
        "ctrl.column", "ctrl.blank_values", "ctrl.positive_values", "ctrl.biological_values",
        "ctrl.blank_override_column", "ctrl.blank_override_values",
        "qc.primer_sequence", "qc.conserved_motif", "qc.head_reads",
        "qc.max_primer_hit_frac", "qc.min_motif_frac", "qc.max_frac_short_reads",
        "decontam.batch_column", "decontam.min_blanks", "err.batch_column",
        "filter.truncLen", "filter.trimLeft",
        "out.study_columns", "out.batch_columns", "katharoseq.cell_count_column",
    )
    #: In modalita' appaiata: i marcatori con cui G07 abbina i file; il
    #: troncamento, il taglio, il primer e il motivo della lettura reverse, da
    #: cui dipendono la scansione dei file reverse, G09 e G10; la soglia del
    #: primer reverse in testa.
    parametri_appaiata: ClassVar[tuple[str, ...]] = (
        "pair.file_markers", "pair.truncLen_reverse", "pair.trimLeft_reverse",
        "pair.primer_sequence_reverse", "pair.conserved_motif_reverse",
        "pair.max_primer_hit_frac_reverse",
    )

    #: Con il primer cercato lungo la lettura: la regola di ricerca, con cui la
    #: scansione conta il primer, cerca il motivo e misura le lunghezze.
    parametri_taglio_cercato: ClassVar[tuple[str, ...]] = (
        "filter.primer_max_mismatch", "filter.primer_max_offset", "filter.primer_untrimmed",
    )

    def __init__(self, *, solleva: bool = True) -> None:
        self.solleva = solleva
        #: Gli esiti di G15 e G12, calcolati in :meth:`esegui` sulla
        #: configurazione completa e registrati da :meth:`calcola` fra gli altri gate.
        self._precondizioni: dict[str, EsitoGate] | None = None
        #: Con il classificatore IdTaxa, il file addestrato da verificare:
        #: letto in :meth:`esegui` dalla configurazione completa, come G12. Il
        #: riferimento non e' fra i parametri della fase (entra in S8): cambiare
        #: classificatore non rifa' la validazione ne' l'inferenza.
        self._addestrato: Path | None = None

    def impronta_dati_esterni(self, config: Config) -> str | None:
        """L'impronta dei dati grezzi da nome, dimensione e istante di modifica dei
        file.
        """
        return impronta_dati_grezzi(config)

    def esegui(self, contesto: StepContext) -> StepResult:
        """G15 e G12 sulla configurazione completa, poi la fase sulla sua vista."""
        completa = Contesto(contesto.config)
        self._precondizioni = {g: esegui_gate(g, completa) for g in ("G15", "G12")}
        if contesto.config.idtaxa:
            self._addestrato = Path(contesto.config.tax.idtaxa_trained)
        try:
            return super().esegui(contesto)
        finally:
            self._precondizioni = None
            self._addestrato = None

    def _classificatore_addestrato(
        self, contesto: StepContext, esiti: tuple[EsitoGate, ...]
    ) -> tuple[tuple[EsitoGate, ...], list[Any]]:
        """Con il classificatore IdTaxa, e dopo che G12 ne ha verificato il
        checksum: carica il file in R e verifica che contenga un classificatore
        addestrato con la tabella dei ranghi. Se non lo contiene, l'esito di
        G12 porta la violazione ``E-S0-28``. Restituisce gli esiti e gli
        artefatti scritti.

        Scoprirlo qui costa qualche secondo; all'assegnazione tassonomica
        costerebbe l'intera inferenza delle varianti. Se il classificatore
        cambia a validazione gia' conclusa, la stessa verifica la fa S8 quando
        lo carica.
        """
        g12 = next(e for e in esiti if e.gate == "G12")
        if not g12.superato:
            return esiti, []
        try:
            esito = esegui_script(
                cartella_r() / NOME_SCRIPT_CLASSIFICATORE,
                {"addestrato": str(self._addestrato)},
                contesto.albero, Fase.INPUT_VALIDATION, passo="S0_CLASSIFICATORE",
                tempo_massimo_s=contesto.config.run.r_timeout_s, logger=contesto.logger,
            )
        except ErrorePipeline as guasto:
            if guasto.codice != "E-S0-28":
                raise
            respinto = replace(g12, superato=False, violazioni=(
                Violazione("E-S0-28", f"tax.idtaxa_trained: {guasto.dettaglio}"),))
            return tuple(respinto if e.gate == "G12" else e for e in esiti), []
        return esiti, list(esito.artefatti)

    def calcola(self, contesto: StepContext) -> Produzione:
        """Esegue i gate G15 e G01-G14 e scrive esiti, crosswalk, inventario e
        statistiche delle letture.
        """
        config = contesto.config
        log = contesto.logger
        inizio = time.perf_counter()

        gate = Contesto(config)
        if self._precondizioni is None:
            raise RuntimeError(
                "S0 si esegue con esegui(): G15 e G12 leggono la configurazione completa"
            )
        esiti = tuple(esegui_tutti(gate, self._precondizioni))
        del_classificatore: list[Any] = []
        if config.idtaxa:
            esiti, del_classificatore = self._classificatore_addestrato(contesto, esiti)
        secondi = time.perf_counter() - inizio

        # L'inventario esiste solo se i gate che lo costruiscono sono passati.
        inventario: Inventario | None = None
        if all(e.superato for e in esiti if e.gate in ("G03", "G04", "G05", "G06", "G11")):
            try:
                inventario = gate.inventario
            except Exception:  # noqa: BLE001 - l'inventario e' accessorio al referto
                inventario = None

        albero = contesto.albero
        artefatti = [
            albero.scrivi_testo(
                Fase.INPUT_VALIDATION, NOME_ESITI, _esiti_json(esiti)
            ),
            *del_classificatore,
        ]
        if inventario is not None:
            artefatti.append(
                albero.scrivi_testo(
                    Fase.INPUT_VALIDATION, NOME_CROSSWALK,
                    _crosswalk_tsv(inventario, config.appaiata),
                )
            )
            artefatti.append(
                albero.scrivi_testo(
                    Fase.INPUT_VALIDATION, NOME_INVENTARIO, _inventario_json(inventario)
                )
            )
        if gate.scansione_gia_fatta:
            artefatti.append(
                albero.scrivi_testo(
                    Fase.INPUT_VALIDATION,
                    NOME_LETTURE,
                    _letture_tsv(gate.scansione, config.qc.head_reads,
                                 con_motivo=config.qc.conserved_motif is not None,
                                 con_primer=config.qc.primer_sequence is not None),
                )
            )
            if config.appaiata:
                artefatti.append(
                    albero.scrivi_testo(
                        Fase.INPUT_VALIDATION,
                        NOME_LETTURE_REVERSE,
                        _letture_tsv(
                            gate.scansione_reverse, config.qc.head_reads,
                            con_motivo=config.pair.conserved_motif_reverse is not None,
                            con_primer=config.pair.primer_sequence_reverse is not None),
                    )
                )

        # Gli avvisi di degradazione (E-S0-15 ed E-S0-19 da G08, E-S0-17 da G11,
        # E-S0-18, E-S0-21 ed E-S0-29 da G10, E-S0-20 da G13) non fermano la
        # fase e finiscono nel suo manifesto; gli altri restano segnalazioni
        # nel log.
        for esito in esiti:
            for avviso in esito.avvisi:
                if voce(avviso.codice).categoria is Categoria.DEGRADAZIONE_AUTOMATICA:
                    contesto.degrada(avviso.codice, avviso.dettaglio, gate=esito.gate)
                else:
                    log.warning(
                        avviso.dettaglio,
                        extra={"gate": esito.gate, "codice": avviso.codice, "fase": "S0"},
                    )

        risultato = RisultatoS0(
            esiti=esiti,
            inventario=inventario,
            artefatti=tuple(a.percorso for a in artefatti),
            secondi=secondi,
        )
        for esito in esiti:
            if esito.eseguito:
                log.info(
                    f"{esito.gate} {'superato' if esito.superato else 'non superato'}",
                    extra={"fase": "S0", "gate": esito.gate, "secondi": round(esito.secondi, 3)},
                )
        # La durata complessiva la registra nel log PipelineStep.esegui, come
        # per ogni fase; nelle metriche, che finiscono nel manifesto, non entra.
        metriche = {
            "campioni": len(inventario) if inventario else 0,
            "avvisi": len(risultato.avvisi),
        }

        if risultato.superata:
            log.info("S0 superata", extra={"fase": "S0", **metriche})
            return Produzione(tuple(artefatti), metriche, True, risultato)

        falliti = risultato.falliti
        for esito in falliti:
            for violazione in esito.violazioni:
                log.error(
                    violazione.dettaglio,
                    extra={"gate": esito.gate, "codice": violazione.codice, "fase": "S0"},
                )

        if self.solleva:
            primo = falliti[0]
            raise ErroreGate(primo.gate, primo.violazioni)
        return Produzione(tuple(artefatti), metriche, False, risultato)


def esegui_s0(config: Config, *, solleva: bool = True) -> RisultatoS0:
    """Esegue i quindici gate e registra l'esito in ``01_input_validation/``.

    Con ``solleva`` falso restituisce il risultato anche quando un gate
    fallisce, invece di sollevare: serve a chi vuole ispezionare l'esito
    completo, per esempio il comando ``validate``.
    """
    contesto = StepContext(
        risolta=risolvi(config),
        albero=AlberoOutput(config.io.out_root),
        logger=ottieni("s0"),
        inventario=None,
    )
    return ValidazioneIngressi(solleva=solleva).esegui(contesto).dettaglio
