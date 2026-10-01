# syntax=docker/dockerfile:1.7

# Ambiente di esecuzione della pipeline: Python e R nella stessa immagine, con le
# versioni bloccate su entrambi i fronti.
#
# Le immagini di partenza sono ancorate per digest e non per tag: un tag può
# essere riassegnato a un'immagine diversa, un digest no. Aggiornare l'ambiente
# significa quindi cambiare esplicitamente un digest qui sotto, mai subire un
# cambiamento deciso a monte.

# bioconductor/bioconductor_docker:RELEASE_3_21 : Ubuntu 24.04, R 4.5, Bioconductor 3.21
ARG BIOCONDUCTOR_IMAGE=bioconductor/bioconductor_docker@sha256:359702d482e70e343d2a436731b7c6deffd1ec184a1f1700de5d64b0c97cacba
# ghcr.io/astral-sh/uv:0.12.5 : serve solo a procurare l'interprete Python
ARG UV_IMAGE=ghcr.io/astral-sh/uv@sha256:e85be844203885286c60ffad8a858d48afb6c5a5c237ca0e67f12e74b8f174b1

FROM ${UV_IMAGE} AS uv

FROM ${BIOCONDUCTOR_IMAGE}

LABEL org.opencontainers.image.title="amplicon16s" \
      org.opencontainers.image.description="Ambiente di esecuzione della pipeline 16S rRNA: Python e R" \
      org.opencontainers.image.licenses="Proprietary"

# Versioni bloccate: cambiarle qui è l'unico modo di cambiare l'ambiente.
ARG PYTHON_VERSION=3.11.16
ARG BIOCONDUCTOR_VERSION=3.21
# Istantanea datata del repository CRAN da cui si installano le dipendenze che
# l'immagine di partenza non contiene. L'immagine Bioconductor punta a
# p3m.dev/.../latest, che cambia nel tempo: una build da zero prenderebbe le
# versioni del giorno, diverse da renv.lock (igraph 2.3.4 contro 2.3.3 il 1
# ottobre 2026). La data e' quella in cui e' stato generato renv.lock.
ARG CRAN_SNAPSHOT=https://p3m.dev/cran/__linux__/noble/2026-09-25

ENV DEBIAN_FRONTEND=noninteractive \
    LANG=C.UTF-8 \
    PYTHON_VERSION=${PYTHON_VERSION} \
    BIOCONDUCTOR_VERSION=${BIOCONDUCTOR_VERSION} \
    UV_PYTHON_INSTALL_DIR=/opt/python \
    UV_LINK_MODE=copy \
    VIRTUAL_ENV=/opt/venv \
    PATH=/opt/venv/bin:$PATH

# --------------------------------------------------------------------------- #
# Python 3.11                                                                  #
# --------------------------------------------------------------------------- #
# L'immagine Bioconductor porta il Python di sistema di Ubuntu 24.04, che è il
# 3.12. La pipeline richiede il 3.11, quindi l'interprete viene procurato a
# parte e bloccato sulla versione esatta: la stessa in sviluppo, in CI e qui.
COPY --from=uv /uv /usr/local/bin/uv

RUN uv python install "${PYTHON_VERSION}" \
 && uv venv --python "${PYTHON_VERSION}" "${VIRTUAL_ENV}" \
 && python --version

# --------------------------------------------------------------------------- #
# R: pacchetti dei calcoli scientifici                                         #
# --------------------------------------------------------------------------- #
COPY container/install_r_packages.R /opt/amplicon16s/container/
RUN CRAN_SNAPSHOT="${CRAN_SNAPSHOT}" Rscript /opt/amplicon16s/container/install_r_packages.R

# dada2 1.36.0 con la correzione dei pareggi di assignTaxonomy, compilato dai
# sorgenti ufficiali verificati per checksum, al posto della versione dei binari.
# In get_best_genus la scelta fra generi a pari probabilità usa
# std::random_device, che set.seed() non controlla: senza la correzione due
# esecuzioni con lo stesso seme non danno la stessa tassonomia. La versione
# compilata è 1.36.0.1 e renv.lock la registra con la correzione e il suo
# checksum (container/dada2/installa_dada2_corretto.sh).
COPY container/dada2/ /opt/amplicon16s/container/dada2/
RUN bash /opt/amplicon16s/container/dada2/installa_dada2_corretto.sh

# --------------------------------------------------------------------------- #
# Pacchetto Python del repository                                              #
# --------------------------------------------------------------------------- #
WORKDIR /opt/amplicon16s

COPY pyproject.toml README.md ./
COPY src/ src/
COPY tests/ tests/
# La configurazione di esempio e' un'istanza dello schema e i test la
# verificano: senza, la suite non potrebbe girare dentro l'immagine.
COPY config/ config/
# Gli strumenti del repository che i test esercitano (per esempio
# scripts/registro_sorgente.py).
COPY scripts/ scripts/
# Gli script R e le loro funzioni condivise. Il pacchetto Python e' installato
# in modo non modificabile e non sta piu' accanto a R/: la variabile dice al
# ponte dove trovarli.
COPY R/ R/
ENV AMPLICON16S_R_DIR=/opt/amplicon16s/R

RUN uv pip install --python "${VIRTUAL_ENV}/bin/python" ".[dev]"

# --------------------------------------------------------------------------- #
# Verifica contro il file di lock R                                            #
# --------------------------------------------------------------------------- #
# La parentesi in "renv.loc[k]" rende facoltativa la COPY: alla primissima build
# il file non esiste ancora, perché viene generato proprio da questa immagine.
# Il glob è accoppiato a pyproject.toml, che esiste sempre: un glob che non
# corrisponde a nulla fa fallire la COPY sul builder legacy, e l'immagine deve
# potersi costruire anche senza BuildKit.
COPY container/snapshot_renv_lock.R container/verify_renv_lock.R container/
COPY pyproject.toml renv.loc[k] ./

RUN Rscript container/verify_renv_lock.R

CMD ["amplicon16s", "--help"]
