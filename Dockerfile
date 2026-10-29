FROM mambaorg/micromamba:1.5.8

LABEL org.opencontainers.image.title="Velvet Worm Slime Proteomics & MD"
LABEL org.opencontainers.image.description="Proteomics and SMD analysis toolkit for Onychophora slime proteins"
LABEL org.opencontainers.image.authors="Gretta Y. Rozenberg"

USER root
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app
COPY pyproject.toml README.md LICENSE ./
COPY src/ ./src/
COPY data/ ./data/

RUN micromamba create -y -n slime-md -c conda-forge \
    python=3.11 \
    biopython \
    numpy \
    scipy \
    pandas \
    matplotlib \
    openmm \
    pdbfixer \
    pytest \
    && micromamba clean --all --yes

ARG MAMBA_DOCKERFILE_ACTIVATE=1
RUN python -m pip install --no-deps .

WORKDIR /app
ENTRYPOINT ["slime-proteo-md"]
CMD ["--help"]
