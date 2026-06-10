# nf-dafseq

A [Nextflow](https://www.nextflow.io/) pipeline for processing **DAF-seq** (Direct
Accessibility Footprinting sequencing) data, end to end:

| Step | Process | What it does | Output |
|------|---------|--------------|--------|
| 1 | `MINIMAP_ALIGN`   | Map ONT reads with `minimap2` + add read group | `{sample}.bam` |
| 2 | `DAFQC_SMK`       | QC + strand decoration by wrapping the [DAF-QC-SMK](https://github.com/StephanieBohaczuk/DAF-QC-SMK) Snakemake | `{sample}.decorated.reads.bam`, QC dashboard |
| 3 | `MARK_DUPLICATES` | Clustering-based combinatorial dedup (deamination-pattern uniqueness) | `{sample}.mkdup.bam` |
| 4 | `PHASE_READS`     | Phase reads by a known SNP or deletion (optional, per sample) | phased BAM + hap1/hap2/non_phased bigWigs |

This pipeline replaces the older four-script manual workflow (with hardcoded paths). All
per-run inputs now live in a single samplesheet.

## Requirements

You need Nextflow (`>=24.04`, tested with 25.10.4) and **one** way to supply tools:

- **Anywhere** — a container runtime: Docker (laptops/cloud) or Apptainer/Singularity (HPC).
  Tools ship as two images: `dhspence/docker-dafseq` (steps 1/3/4; the lab image from
  [dhslab-docker-images](https://github.com/dhslab/dhslab-docker-images)) and
  `ghcr.io/bonney96/nf-dafseq-dafqc` (step 2 — the wrapped DAF-QC-SMK Snakemake, with its conda
  envs prebuilt; built from `containers/dafqc`). Select with
  `-profile docker|singularity|apptainer`.
- **WashU RIS in-house** — `-profile washu` uses Lmod (`module load labtools`/`pixi`) + SLURM
  instead of containers, and restores the on-cluster `--ref`/`--dafqc_repo` defaults. On first
  run Snakemake builds its conda envs into `params.conda_prefix` (slow once, then cached).

## Samplesheet (`assets/samplesheet.tsv`)

Tab-separated, one row per sample:

| column | required | description |
|--------|----------|-------------|
| `sample`      | yes | unique sample name |
| `fastq`       | yes | path to the ONT FASTQ |
| `regions`     | yes | comma-separated target regions, `chr:start-end[,chr:start-end...]` |
| `phase_type`  | no  | `snp`, `del`, or empty (no phasing) |
| `phase_value` | no  | for `snp`: `chr:pos,ref,alt` (e.g. `chr3:128491529,C,T`); for `del`: `chr:start-end` |

TSV is used (not CSV) so SNP values containing commas don't clash with the delimiter.

## Running

```bash
# Dry wiring check — no tools or data needed:
nextflow run . -stub -profile test

# Anywhere with a container runtime (--ref is required):
nextflow run . -profile docker    --input sheet.tsv --ref genome.fa --outdir results
nextflow run . -profile apptainer --input sheet.tsv --ref genome.fa --outdir results

# WashU RIS in-house (Lmod + SLURM; --ref/--dafqc_repo default to cluster paths):
nextflow run . -profile washu --input assets/samplesheet.tsv --outdir results
```

`--ref` is required for every real (non-`washu`) run. Other common overrides: `--dafqc_repo`,
`--conda_prefix`, `--outdir`, `--decorated_samplesize`, `--del_overlap_frac`. Profiles compose,
e.g. `-profile docker,local`.

## Outputs

```
results/<sample>/align/    <sample>.bam, <sample>.decorated.reads.bam, QC dashboard
results/<sample>/dedup/     <sample>.mkdup.bam, cluster table + plots
results/<sample>/phasing/   <sample>.filtered_phased.bam, hap1/hap2/non_phased .bw + .bedgraph
results/pipeline_info/      timeline, report, trace, DAG
```

## Notes / TODO

- Step 2 wraps the DAF-QC-SMK Snakemake as one process (pinned + baked into the `nf-dafseq-dafqc`
  image). The layout keeps it isolated so its rules could be ported to native Nextflow processes
  later (Tier 3 in `docs/PORTABILITY_HANDOFF.md`).
- nf-core conventions (schema validation, lint/CI, institutional configs) are the next tier and
  not yet applied.
