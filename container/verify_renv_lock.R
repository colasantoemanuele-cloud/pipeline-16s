#!/usr/bin/env Rscript
# Verifica che i pacchetti R installati nell'immagine corrispondano a renv.lock.
#
# renv.lock viene generato dall'immagine già costruita, quindi alla primissima
# build non esiste ancora e la verifica viene saltata. Una volta versionato, il
# file diventa vincolante: ogni build successiva fallisce se una versione non
# coincide.

lock <- "renv.lock"

if (!file.exists(lock)) {
  cat("renv.lock assente: verifica saltata (atteso solo alla prima build)\n")
  quit(status = 0)
}

registrati <- jsonlite::fromJSON(lock, simplifyVector = FALSE)$Packages

discordanze <- character()
for (nome in names(registrati)) {
  attesa <- registrati[[nome]]$Version

  if (!requireNamespace(nome, quietly = TRUE)) {
    discordanze <- c(discordanze, sprintf("%s: atteso %s, non installato", nome, attesa))
    next
  }

  # Il confronto avviene fra oggetti package_version e non fra stringhe: R
  # tratta "-" e "." come separatori equivalenti, quindi as.character() di
  # packageVersion("1.90.0-1") restituisce "1.90.0.1" e un confronto testuale
  # segnalerebbe come discordante ogni versione che contiene un trattino.
  trovata <- packageVersion(nome)
  if (trovata != package_version(attesa)) {
    discordanze <- c(
      discordanze,
      sprintf("%s: atteso %s, trovato %s", nome, attesa, as.character(trovata))
    )
  }
}

if (length(discordanze)) {
  # Le discordanze vanno su stderr una per riga: il messaggio di stop() viene
  # troncato da R oltre una certa lunghezza e nasconderebbe parte dell'elenco.
  message("l'immagine non corrisponde a renv.lock:")
  for (d in discordanze) message("  ", d)
  stop(sprintf("%d pacchetti discordanti rispetto a renv.lock", length(discordanze)),
       call. = FALSE)
}

cat(sprintf("renv.lock verificato: %d pacchetti alle versioni registrate\n",
            length(registrati)))
