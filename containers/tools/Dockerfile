# nf-dafseq-tools : minimap2 + samtools + the python stack for dedup & phasing.
# Used by steps 1 (MINIMAP_ALIGN), 3 (MARK_DUPLICATES) and 4 (PHASE_READS).
#
# Build:   docker build -t ghcr.io/bonney96/nf-dafseq-tools:0.1.0 containers/tools
# Versions are pinned to match the DAF-QC-SMK reference (samtools 1.22.1, minimap2 2.30,
# pysam 0.23.3, matplotlib 3.10.3); see containers/tools/env.yml.
FROM mambaorg/micromamba:1.5.10

LABEL org.opencontainers.image.source="https://github.com/Bonney96/nf-dafseq" \
      org.opencontainers.image.description="nf-dafseq tools image: minimap2, samtools, pysam + dedup/phasing python deps"

# Solve the env into base so every binary is on PATH without `micromamba activate`.
COPY --chown=$MAMBA_USER:$MAMBA_USER env.yml /tmp/env.yml
RUN micromamba install -y -n base -f /tmp/env.yml && \
    micromamba clean --all --yes

# Nextflow invokes .command.sh directly (no entrypoint activation), so put the env on PATH.
ENV PATH="/opt/conda/bin:$PATH"

# Sanity: fail the build if any required tool is missing.
RUN minimap2 --version && samtools --version | head -1 && \
    bedGraphToBigWig 2>&1 | head -1; \
    python -c "import pysam, pandas, numpy, scipy, sklearn, matplotlib, seaborn; print('py deps OK')"
