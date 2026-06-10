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

Run on the lab cluster. The pipeline uses environment modules + pixi (no per-process
containers in this version):

- `module load nextflow`  (tested with 25.10.4)
- `module load labtools`  — provides minimap2, samtools, bedGraphToBigWig, and a Python with
  pysam/numpy/pandas/scipy/scikit-learn/seaborn (used by steps 1, 3, 4). Loaded automatically
  per process via `beforeScript`.
- `module load pixi`      — used by step 2 to run the wrapped DAF-QC-SMK Snakemake. Loaded
  automatically per process.

Step 2 calls the **DAF-QC-SMK** repo set by `params.dafqc_repo`. On first real run, Snakemake
builds its conda environments into `params.conda_prefix` (slow once, then cached).

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
module load nextflow

# Dry wiring check — no tools or data needed:
nextflow run . -stub -profile test

# Real run on the cluster:
nextflow run . -profile slurm --input assets/samplesheet.tsv --outdir results
```

Common overrides: `--ref`, `--dafqc_repo`, `--conda_prefix`, `--outdir`,
`--decorated_samplesize`, `--del_overlap_frac`.

## Outputs

```
results/<sample>/align/    <sample>.bam, <sample>.decorated.reads.bam, QC dashboard
results/<sample>/dedup/     <sample>.mkdup.bam, cluster table + plots
results/<sample>/phasing/   <sample>.filtered_phased.bam, hap1/hap2/non_phased .bw + .bedgraph
results/pipeline_info/      timeline, report, trace, DAG
```

## Notes / TODO

- `params.dafqc_repo` currently defaults to a personal path; clone DAF-QC-SMK to a shared
  location and update the default.
- This version wraps the DAF-QC-SMK Snakemake as one process. The layout keeps step 2 isolated
  so its rules could be ported to native Nextflow processes later.
- A container (`apptainer`) profile can be added later for cross-cluster reproducibility.
