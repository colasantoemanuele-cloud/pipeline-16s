"""Configurazione assistita: una configurazione di partenza scritta dai dati.

Il sottocomando ``init`` guarda le letture e i metadati di un dataset e scrive
un file di configurazione in cui ogni valore dichiara da dove viene: rilevato
sui file, proposto da un'euristica, oppure lasciato da scegliere a chi conduce
l'analisi. Il file non e' eseguibile cosi' com'e': ``validate`` elenca cio' che
resta da confermare e da decidere.
"""
