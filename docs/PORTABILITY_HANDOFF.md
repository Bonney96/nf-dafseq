# nf-dafseq — Portability Handoff (Road to nf-core)

**Audience:** an agent/developer tasked with planning the work to make `nf-dafseq`
runnable outside WashU RIS infrastructure, and (optionally) bringing it toward nf-core
standards. This is a **planning brief**, not an implementation spec — produce a tiered
plan, surface decisions, then implement once approved.

**Branch:** `road-to-nf-core` (branched from `main` @ commit `3db95c8`).

---

## 1. What this pipeline is (current state)

A lightweight Nextflow DSL2 pipeline that ports a previously-manual 4-step DAF-seq workflow.
Repo: `github.com/Bonney96/nf-dafseq`. Structure:

```
main.nf                         # parses TSV samplesheet, chains 4 processes
nextflow.config                 # params + profiles (standard/slurm/test)
conf/{base,slurm}.config        # resource labels; SLURM executor
modules/local/
  minimap_align.nf              # 1. minimap2 + read group
  dafqc_smk.nf                  # 2. WRAPS the DAF-QC-SMK Snakemake via pixi
  mark_duplicates.nf            # 3. clustering combinatorial dedup
  phase_reads.nf                # 4. phase by SNP/deletion (optional per sample)
bin/                            # step 3 & 4 python scripts (on PATH automatically)
assets/samplesheet*.tsv         # sample,fastq,regions,phase_type,phase_value
```

**Pipeline flow:** `MINIMAP_ALIGN → DAFQC_SMK → MARK_DUPLICATES → (filter variant rows) → PHASE_READS`.

**Validation status:** VALIDATED on a real SLURM run (3 MOLM13 samples). Phasing hap1/hap2
counts reproduce the original manual workflow **exactly** (GATA2 del 213/331; 287XL snp
301/151; identical totals/written/strand). pandas 2.x worked for the dedup script. So any
portability refactor must preserve these numbers — use them as a regression check.

**Design decisions already made by the user (do not relitigate without asking):**
- Step 2 **wraps** the published DAF-QC-SMK Snakemake rather than re-porting it natively.
- Lightweight DSL2, *not* the full nf-core template (yet).
- Current dependency model: environment modules + pixi.

---

## 2. Why it is not portable today (concrete blockers)

| # | Blocker | Where | Impact |
|---|---------|-------|--------|
| 1 | `beforeScript 'module load labtools'` (steps 1/3/4) and `'module load pixi'` (step 2) | `modules/local/*.nf` | Lmod modules are WashU-RIS-only. **No `conda`/`container` directives exist on any process**, so there is no alternative way to resolve tools. This is the #1 blocker. |
| 2 | Hardcoded filesystem paths | `nextflow.config`: `params.ref` (`/storage2/.../hg38_mgi_patch.fa`), `params.dafqc_repo` (`/storage2/.../mohamed/.../DAF-QC-SMK`) | External users have neither the reference nor the Snakemake repo at those paths. |
| 3 | SLURM + partition hardwired | `conf/slurm.config` (`executor='slurm'`, `queue='general-cpu'`) | Assumes SLURM and a WashU partition name. No abstraction for local/LSF/PBS/SGE/k8s/cloud. |
| 4 | Nested workflow engine (step 2) | `modules/local/dafqc_smk.nf` | Nextflow → pixi → Snakemake → conda. Portable only insofar as pixi can be installed and the DAF-QC-SMK repo is present. The hardest part to make portable. |
| 5 | Tool envs are implicit | steps 3/4 rely on `labtools` python having pysam/scipy/scikit-learn/seaborn; step 4 needs UCSC `bedGraphToBigWig` | No declared environment ⇒ not reproducible off-cluster. |

---

## 3. Key facts to design against (verified on this cluster, 2026-06)

**Available here as modules:** `nextflow/25.10.4`, `apptainer/1.4.5`, `pixi/0.56.0`.
Apptainer means we **can build container images on-cluster**. No docker/singularity binary
directly (apptainer is the singularity successor).

**Cluster specifics (for the WashU profile that outside users won't select):**
SLURM partition `general-cpu`; default account `compute2-dspencer` (no `--account` flag needed).

**Exact tool versions are already pinned in the DAF-QC-SMK repo** (reuse these for container/conda specs):
- `workflow/envs/cmd.yaml`: `samtools==1.22.1`, `minimap2==2.30`, `pbmarkdup==1.2.0`, `rustybam==0.1.34`
- `workflow/envs/python.yaml`: `pysam==0.23.3`, `numpy==1.24`, `pandas==1.4`, `matplotlib==3.10.3`, `pyabpoa==1.5.4`, `panel==1.7.4`
- `pixi.toml`: `snakemake==8.21`, plus conda/awscli
- `workflow/profiles/default/config.yaml`: `software-deployment-method: [apptainer, conda]` (Snakemake already supports conda/containers internally — a portability asset)

**What steps 1/3/4 actually need (for writing their conda/container specs):**
- Step 1 (`minimap_align.nf`): `minimap2`, `samtools`.
- Step 3 (`bin/mark_duplicates.decorated.py`): python with `pysam, pandas, numpy, scipy, scikit-learn, seaborn, matplotlib`; plus `samtools` for indexing.
- Step 4 (`bin/phase_reads_del.py`): python with `pysam`; plus UCSC `bedGraphToBigWig` (bioconda: `ucsc-bedgraphtobigwig`).

**Gotchas already discovered (must be preserved through any refactor):**
- `bin/mark_duplicates.decorated.py` needs its `#!/usr/bin/env python3` shebang (added) — without it, calling by name runs under /bin/sh.
- `DAFQC_SMK` has `maxForks 1` to avoid concurrent Snakemake conda-env builds racing into the shared `SNAKEMAKE_CONDA_PREFIX`. A containerized step 2 (envs baked in) would remove this constraint.
- The pixi `snakemake` task already injects `-s <repo>/workflow/Snakefile`; don't pass `-s` again.

---

## 4. The three tiers to plan

For each tier, produce: concrete task list, files to touch, key decisions to put to the user,
and a verification step. **Regression check for every tier:** re-run the 3 MOLM13 samples and
confirm hap1/hap2 still = 213/331 and 301/151.

### Tier 1 — "Runs anywhere" (portability; high payoff, modest effort)
Goal: `nextflow run Bonney96/nf-dafseq -profile docker --input sheet.tsv --ref genome.fa`
works on an arbitrary machine with a container runtime, no Lmod/WashU assumptions.

Plan should cover:
- **Software packaging:** add `conda` + `container` directives to all four processes; remove
  the `module load` `beforeScript`s (or keep them only inside a `washu` profile). Reuse the
  pinned versions in §3. Decide: per-process biocontainers vs a small number of custom images.
- **Step 2 containerization (the crux):** build one image bundling pixi + DAF-QC-SMK + its
  Snakemake conda envs, OR vendor DAF-QC-SMK as a git submodule/pinned release and document a
  pixi bootstrap. Recommend the container approach (apptainer can build it here). Decide how
  DAF-QC-SMK is obtained/pinned (submodule vs release tag vs baked into image).
- **Profiles:** add `docker`, `singularity`, `apptainer`, `conda`, `local`, `test` profiles
  for *software/execution*; move SLURM+`general-cpu`+account into a separate `washu` profile.
- **References/inputs:** make `--ref` required (drop the `/storage2` default); keep samplesheet
  as the only other required input. Consider an optional URL/igenomes fetch.
- **Decision to surface:** which container runtimes to officially support first (Docker for
  laptops/cloud, Apptainer/Singularity for HPC).

### Tier 2 — nf-core conventions (standardization; optional, à la carte)
Goal: recognizable, validated, sharable — without necessarily the full template rewrite.
Plan should cover, in increasing commitment:
- **nf-schema:** JSON-schema validation for the samplesheet and params (`nextflow_schema.json`,
  `assets/schema_input.json`). Catches bad input pre-launch.
- **Config layout & linting:** align to nf-core directory/config conventions; add `nf-core lint`,
  GitHub Actions CI (stub test on PRs), versioned releases, `CHANGELOG`, `CITATIONS`.
- **nf-core/configs:** contribute a `washu`/`ris` institutional profile so other sites get a
  tested profile for free; outside users select theirs.
- **MultiQC / docs** conventions.
- **Decision to surface:** adopt nf-core conventions on the existing lightweight pipeline vs
  regenerate from the `nf-core create` template (bigger, more opinionated). Note the step-2
  Snakemake wrap does **not** fit the nf-core-modules model cleanly — strict compliance pushes
  toward Tier 3.

### Tier 3 — Native port of DAF-QC-SMK (only if needed)
Goal: remove the nested-engine dependency; reimplement DAF-QC-SMK's ~12 rules as native nf
processes with their own container directives. Largest lift; enables true nf-core-modules
compliance. Plan should map each Snakemake rule (`workflow/rules/rules.smk`: deduplicate/align/
targeting_qc/sequence_qc/decorate_strands/build_consensus/make_dashboard/...) to an nf process,
reusing the existing `workflow/scripts/*.py`. Only pursue if the wrap becomes a maintenance
burden or full nf-core compliance is required. Flag the consensus/PacBio path as extra scope
(this lab's DAF-seq is ONT; consensus is PacBio-only).

---

## 5. Recommended sequencing & open items

1. **Tier 1 first** — it is what actually answers "someone outside WashU wants to run it," and
   it makes the pending HOXA run reproducible off-cluster too.
2. Tier 2 only if community sharing/citation matters.
3. Tier 3 only on demand.

**Do not** change the four-step science or the wrap decision without asking the user.
**Always** end with the MOLM13 regression check (§4).

**Unrelated pending task (not portability):** a real `HOXA_POOL` run is still owed — the user
has provided target regions `chr7:27092000-27200000`; that work can proceed independently of
this branch on `main`.

**Repo hygiene note:** never add AI/Claude co-author trailers to commits (standing user rule);
`.gitignore` covers `results*/`, `work/`, `*.log`, `.snakemake-conda-envs/`.
