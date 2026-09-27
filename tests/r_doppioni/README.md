# Script R doppioni

Script che esistono solo per i test del ponte verso R (`tests/test_w08_rbridge.py`).
Non sono script di fase e non vanno mai invocati dalla pipeline: stanno qui, fuori da
`R/`, perché nessuno possa scambiarli per quelli veri.

- `successo.R`: riesce, scrive un artefatto e un file di tracciamento delle letture.
- `fallimento_dichiarato.R`: si interrompe dichiarando un codice del catalogo; con il
  parametro `modo` può invece interrompersi con un errore non catalogato, uscire senza
  dichiarare, morire per un errore di segmentazione o restare bloccato.
- `memoria.R`: esaurisce la memoria: con un'allocazione che fallisce sotto il limite
  imposto dal ponte al processo figlio, oppure facendosi uccidere con `SIGKILL`, come
  farebbe il sistema operativo.
- `parametro_non_ricevuto.R` — legge un parametro che la fase non gli ha passato: la
  lettura deve fallire, invece di restituire `NULL` e proseguire con un valore predefinito.
