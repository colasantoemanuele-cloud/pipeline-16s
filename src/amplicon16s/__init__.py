"""Pipeline 16S rRNA — single-end Illumina.

Pacchetto di produzione: orchestrazione in Python delle 15 fasi (S0-S14) che
terminano producendo l'oggetto phyloseq serializzato ``ps_final.rds``. I calcoli
scientifici sono delegati a R come processi separati; il confine fra i due
linguaggi passa per il filesystem (parametri JSON in ingresso, artefatti su disco
in uscita).

Le analisi ecologiche a valle vivono nel pacchetto separato ``amplicon16s_eco``.
"""

__version__ = "0.1.0"

__all__ = ["__version__"]
