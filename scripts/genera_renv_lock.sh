#!/usr/bin/env bash
# Rigenera renv.lock dalle versioni installate nell'immagine.
#
# Va eseguito dopo ogni modifica a container/Dockerfile che tocchi i pacchetti
# R: il file di lock deve sempre riflettere l'immagine, non precederla.
#
#   scripts/genera_renv_lock.sh [tag]

set -euo pipefail

IMMAGINE="${1:-amplicon16s:dev}"
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
