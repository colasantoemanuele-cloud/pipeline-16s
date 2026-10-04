#!/usr/bin/env bash
# Rigenera renv.lock dalle versioni installate nell'immagine.
#
# Va eseguito dopo ogni modifica a container/Dockerfile che tocchi i pacchetti
# R: il file di lock deve sempre riflettere l'immagine, non precederla.
#
#   scripts/genera_renv_lock.sh <immagine>
#
# <immagine> e' l'immagine appena costruita da cui leggere le versioni: e'
# obbligatoria, perche' un valore predefinito potrebbe indicare un'immagine
# vecchia e produrre un lock che non corrisponde a quella nuova.

set -euo pipefail

if [ "$#" -ne 1 ] || [ -z "$1" ]; then
  echo "uso: scripts/genera_renv_lock.sh <immagine>" >&2
  echo "  <immagine>: il nome dell'immagine costruita da container/Dockerfile, per" >&2
  echo "  esempio quella indicata in test.txt (sezione 1.3)" >&2
  exit 2
fi
IMMAGINE="$1"
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

echo "Immagine: ${IMMAGINE}"
echo "Destinazione: ${REPO}/renv.lock"

docker run --rm \
  -v "${REPO}:/out" \
  -e RENV_LOCK_OUT=/out/renv.lock \
  --entrypoint Rscript \
  "${IMMAGINE}" /opt/amplicon16s/container/snapshot_renv_lock.R

# Il container gira come root: il file appena scritto appartiene a root.
if [ ! -w "${REPO}/renv.lock" ]; then
  echo "renv.lock non scrivibile dall'utente corrente, correggo i permessi" >&2
  docker run --rm -v "${REPO}:/out" --entrypoint chown \
    "${IMMAGINE}" "$(id -u):$(id -g)" /out/renv.lock
fi

echo "Fatto."
