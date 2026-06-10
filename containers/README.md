# Container images for nf-dafseq

Tier 1 portability ships the pipeline's tools as containers instead of `module load`. Two
images cover all four steps:

| Image | Steps | Build dir |
|-------|-------|-----------|
| `ghcr.io/bonney96/nf-dafseq-tools` | 1 (align), 3 (dedup), 4 (phasing) | `containers/tools` |
| `ghcr.io/bonney96/nf-dafseq-dafqc` | 2 (DAF-QC-SMK wrap) | `containers/dafqc` |

The pipeline references these via `container` directives in `modules/local/*.nf`. A container
engine is only used when you select a profile that enables one (`-profile docker|singularity|
apptainer`); under `-profile washu` the directives are inert and the cluster `module load`s are
used instead.

## Build & push

```bash
TAG=0.1.0
docker build -t ghcr.io/bonney96/nf-dafseq-tools:$TAG containers/tools
docker build -t ghcr.io/bonney96/nf-dafseq-dafqc:$TAG containers/dafqc
docker push ghcr.io/bonney96/nf-dafseq-tools:$TAG
docker push ghcr.io/bonney96/nf-dafseq-dafqc:$TAG
```

On the cluster (no docker daemon) the same OCI images are pulled by Apptainer automatically
when you run `-profile apptainer`; no separate build is required. To pre-pull/convert:

```bash
apptainer pull docker://ghcr.io/bonney96/nf-dafseq-tools:0.1.0
```

Keep the `:<tag>` in the module `container` directives in sync with what you push.

## Notes on the `dafqc` image (step 2)

- DAF-QC-SMK is pinned by `ARG DAFQC_COMMIT` (currently `43184be`, the validated commit) from
  `https://github.com/StergachisLab/DAF-QC-SMK`. Bump that arg to update.
- Its two Snakemake conda envs (`workflow/envs/cmd.yaml`, `python.yaml`) are **pre-built** at
  image-build into `/opt/snakemake-conda-envs` via `snakemake --conda-create-envs-only`, so the
  step does no env-solving or network access at runtime. Snakemake 8 keys envs by file content,
  so the runtime prefix matches the prebuilt one.
- The prebake uses placeholder inputs just to let Snakemake construct the DAG. If a future
  DAF-QC-SMK revision rejects the placeholders during DAG building, swap the prebake `RUN` to
  use the repo's real test data instead:
  ```dockerfile
  RUN pixi run --manifest-path /opt/DAF-QC-SMK/pixi.toml test-data && \
      cd /opt/DAF-QC-SMK/dafqc-test-data && \
      pixi run --manifest-path /opt/DAF-QC-SMK/pixi.toml snakemake \
        --configfile test.yaml --software-deployment-method conda \
        --conda-prefix /opt/snakemake-conda-envs --conda-create-envs-only --cores 1
  ```

## Version pins

Mirrors the DAF-QC-SMK repo (`workflow/envs/*.yaml`): `samtools==1.22.1`, `minimap2==2.30`,
`pysam==0.23.3`, `matplotlib==3.10.3`. The `tools` image uses pandas 2.x (our dedup script is
validated there); the `dafqc` image keeps DAF-QC-SMK's internal `pandas==1.4` pin.
