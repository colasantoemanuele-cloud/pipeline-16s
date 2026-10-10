"""Le soglie che fermano o avvisano: quando un controllo e' spento, e come si
propone una soglia dalle misure di un'esecuzione conclusa.

* :mod:`amplicon16s.soglie.tabella`: per ogni soglia la direzione, il limite
  oltre il quale il suo controllo non puo' scattare, il margine del
  suggerimento; il riconoscimento dei controlli spenti;
* :mod:`amplicon16s.soglie.misure`: la statistica che ogni controllo confronta
  con la soglia, letta dagli artefatti di un'esecuzione;
* :mod:`amplicon16s.soglie.comando`: il sottocomando ``soglie``.
"""
