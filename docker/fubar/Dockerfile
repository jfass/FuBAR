FROM mambaorg/micromamba
SHELL ["/bin/bash", "-o", "pipefail", "-c"]

RUN micromamba install -y -n base -c conda-forge -c bioconda \
      r-base=4.1.2 \
      bioconductor-rsamtools \
      bioconductor-shortread \
      bioconductor-rtracklayer \
      bioconductor-rsubread \
      samtools \
      tabix \
      picard \
      seqtk \
      igv-reports \
    && micromamba clean --all --yes
COPY --chmod=644 fubar.R /opt/

COPY Dockerfile /opt/

