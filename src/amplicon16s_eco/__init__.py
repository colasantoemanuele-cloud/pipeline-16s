"""Analisi ecologiche a valle della pipeline 16S.

Pacchetto separato dalla pipeline di produzione: opera sempre su una copia
dell'oggetto phyloseq prodotto da ``amplicon16s`` (``ps_final.rds``), mai
sull'originale.
"""

__version__ = "0.1.0"

__all__ = ["__version__"]
