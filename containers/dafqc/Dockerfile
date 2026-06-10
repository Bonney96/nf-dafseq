# nf-dafseq-dafqc : self-contained image for step 2 (DAFQC_SMK).
# Bundles pixi + a pinned DAF-QC-SMK checkout + Snakemake 8.21, with the workflow's two
# conda envs (cmd.yaml, python.yaml) PRE-BUILT into a fixed prefix so the step needs no
# network or env-solving at runtime. This is what removes the `maxForks 1` race.
#
# Build:  docker build -t ghcr.io/bonney96/nf-dafseq-dafqc:0.1.0 containers/dafqc
FROM debian:bookworm-slim

ARG DAFQC_COMMIT=43184be6a44f436d4a32be310fbb43df990a82b2
ARG PIXI_VERSION=0.56.0

LABEL org.opencontainers.image.source="https://github.com/Bonney96/nf-dafseq" \
      org.opencontainers.image.description="nf-dafseq step 2: pixi + DAF-QC-SMK (pinned) + prebuilt Snakemake conda envs"

ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y --no-install-recommends \
        git ca-certificates curl bash procps && \
    rm -rf /var/lib/apt/lists/*

# --- pixi (pinned) ---
RUN curl -fsSL https://pixi.sh/install.sh | PIXI_VERSION=v${PIXI_VERSION} PIXI_HOME=/opt/pixi bash
ENV PATH="/opt/pixi/bin:${PATH}"

# --- DAF-QC-SMK pinned at a known-good commit ---
RUN git clone https://github.com/StergachisLab/DAF-QC-SMK.git /opt/DAF-QC-SMK && \
    git -C /opt/DAF-QC-SMK checkout ${DAFQC_COMMIT}

# --- pixi env: snakemake 8.21 + conda (from the repo's pixi.lock) ---
RUN pixi install --manifest-path /opt/DAF-QC-SMK/pixi.toml

# --- pre-build the snakemake conda envs into a fixed, image-local prefix ---
# Uses a throwaway config with placeholder inputs so Snakemake can build the DAG and create
# every rule's conda env, then exits (--conda-create-envs-only) without running any job.
# Snakemake 8 keys envs by file CONTENT, so the same env.yaml resolves to the same prefix at
# runtime -> no rebuild. If the DAG proves too strict about placeholder inputs on a future
# DAF-QC-SMK version, swap in the repo's real test data (`pixi run ... test-data`); see README.
ENV SNAKEMAKE_CONDA_PREFIX=/opt/snakemake-conda-envs
RUN mkdir -p /opt/snakemake-conda-envs /opt/prebake && cd /opt/prebake && \
    printf '>chr1\n%s\n' "$(head -c 2000 /dev/zero | tr '\0' 'A')" > ref.fa && \
    : > reads.bam && \
    printf 'sample\tfile\tregs\n'                              >  config.tbl && \
    printf 'prebake\t/opt/prebake/reads.bam\tchr1:1-100\n'     >> config.tbl && \
    printf 'ref: /opt/prebake/ref.fa\nmanifest: /opt/prebake/config.tbl\nplatform: ont\nis_fastq: False\n' > config.yaml && \
    pixi run --manifest-path /opt/DAF-QC-SMK/pixi.toml \
        snakemake --configfile config.yaml \
        --software-deployment-method conda \
        --conda-prefix /opt/snakemake-conda-envs \
        --conda-create-envs-only --cores 1 && \
    rm -rf /opt/prebake

# Defaults the pipeline points at (overridable by the washu profile back to cluster paths).
ENV DAFQC_REPO=/opt/DAF-QC-SMK
