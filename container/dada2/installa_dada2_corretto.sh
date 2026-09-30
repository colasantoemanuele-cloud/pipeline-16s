#!/usr/bin/env bash
# Installa dada2 1.36.0 con la correzione dei pareggi, al posto della versione
# ufficiale installata dai binari di Bioconductor.
#
# Perché. In dada2 1.36.0 (src/taxonomy.cpp, get_best_genus) assignTaxonomy
# sceglie fra generi a pari probabilità con std::random_device, un seme preso dal
# sistema che set.seed() non controlla: due esecuzioni con lo stesso seme danno
# risultati diversi. Sul dataset di riferimento cambiano 3-4 varianti su 12.045
# e il bootstrap di un centinaio. La correzione (dada2-1.36.0-pareggi.patch)
# inizializza il generatore dei pareggi con un seme che dipende dal seme di R,
# dall'indice della sequenza e dal caso (classificazione diretta, complementare
# inversa, replica di bootstrap): ogni pareggio ha un seme fisso qualunque thread
# lo esegua. La scelta resta uniforme fra i generi a pari probabilità; fuori dai
# pareggi nulla cambia.
#
# Come. I sorgenti ufficiali si scaricano a versione fissa e si verificano per
# checksum prima di applicare la correzione. La versione compilata è 1.36.0.1 e
# porta nel DESCRIPTION la versione di partenza, il nome della correzione e il
# suo SHA-256: renv.lock li registra, e container/verify_renv_lock.R li confronta.

set -euo pipefail

BASE="1.36.0"
CORRETTA="1.36.0.1"
SHA256_SORGENTI="314658ecbad252326b6a2b0eafedf5ea5389f4156f5225f3bcac0fb82950bb53"
CORREZIONE="dada2-${BASE}-pareggi.patch"
QUI="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# Il file resta al primo indirizzo finché è la versione corrente della release,
# poi passa all'archivio: il checksum garantisce che sia comunque lo stesso.
INDIRIZZI=(
  "https://bioconductor.org/packages/${BIOCONDUCTOR_VERSION}/bioc/src/contrib/dada2_${BASE}.tar.gz"
  "https://bioconductor.org/packages/${BIOCONDUCTOR_VERSION}/bioc/src/contrib/Archive/dada2/dada2_${BASE}.tar.gz"
)

lavoro="$(mktemp -d)"
trap 'rm -rf "${lavoro}"' EXIT
cd "${lavoro}"

for indirizzo in "${INDIRIZZI[@]}"; do
  if curl -fsSL --retry 3 -o sorgenti.tar.gz "${indirizzo}"; then
    break
  fi
done
echo "${SHA256_SORGENTI}  sorgenti.tar.gz" | sha256sum -c -

tar -xzf sorgenti.tar.gz
patch -p1 -d dada2 --forward < "${QUI}/${CORREZIONE}"

# La versione compilata si distingue da quella ufficiale nel pacchetto: numero
# di versione, versione di partenza, correzione e suo checksum. Il campo
# Repository si toglie perché questa build non viene dal repository Bioconductor.
sha_correzione="$(sha256sum "${QUI}/${CORREZIONE}" | cut -d' ' -f1)"
sed -i -e "s/^Version: ${BASE}\$/Version: ${CORRETTA}/" -e '/^Repository:/d' dada2/DESCRIPTION
cat >> dada2/DESCRIPTION <<EOF
Amplicon16sBase: ${BASE}
Amplicon16sPatch: ${CORREZIONE}
Amplicon16sPatchSHA256: ${sha_correzione}
EOF

# Nella stessa libreria della versione ufficiale, che viene sostituita: non
# devono restarne due copie, una delle quali verrebbe caricata per prima.
libreria="$(Rscript -e 'cat(dirname(find.package("dada2")))')"
R CMD INSTALL -l "${libreria}" dada2

Rscript -e "
  v <- as.character(packageVersion('dada2'))
  d <- packageDescription('dada2')
  if (v != '${CORRETTA}' || !identical(d[['Amplicon16sPatch']], '${CORREZIONE}')) {
    stop('dada2 installato ', v, ', correzione ', d[['Amplicon16sPatch']], ': atteso ${CORRETTA}')
  }
  cat('dada2', v, 'con', d[['Amplicon16sPatch']], '\n')
"
