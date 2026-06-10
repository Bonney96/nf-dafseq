// Step 2: QC + strand decoration, by wrapping the published DAF-QC-SMK Snakemake.
// Invoked per-sample: writes a 1-row config.tbl + a config.yaml, then runs Snakemake
// LOCALLY inside this job (--cores, not its own slurm executor) via pixi.
// Conda envs are forced (--sdm conda) and cached in a persistent prefix so they build once.

process DAFQC_SMK {
    tag { sample }
    label 'process_high'
    publishDir "${params.outdir}/${sample}", mode: 'copy',
        saveAs: { fn -> fn.startsWith("results/${sample}/") ? fn.substring("results/${sample}/".length()) : fn }

    beforeScript 'module load pixi'

    input:
    tuple val(sample), path(bam), path(bai), val(regions), val(phase_type), val(phase_value)

    output:
    tuple val(sample),
          path("results/${sample}/align/${sample}.decorated.reads.bam"),
          path("results/${sample}/align/${sample}.decorated.reads.bam.bai"),
          val(regions), val(phase_type), val(phase_value),                  emit: decorated
    path "results/${sample}/qc",                                            emit: qc, optional: true

    script:
    """
    # --- per-sample manifest + config for DAF-QC-SMK ---
    printf 'sample\\tfile\\tregs\\n'        >  config.tbl
    printf '%s\\t%s\\t%s\\n' "${sample}" "\$(realpath ${bam})" "${regions}" >> config.tbl

    cat > config.yaml <<EOF
ref: ${params.ref}
manifest: config.tbl
platform: ${params.platform}
is_fastq: ${params.is_fastq ? 'True' : 'False'}
chimera_cutoff: ${params.chimera_cutoff}
min_deamination_count: ${params.min_deamination_count}
end_tolerance: ${params.end_tolerance}
decorated_samplesize: ${params.decorated_samplesize}
benchmark: ${params.benchmark ? 'True' : 'False'}
EOF

    # --- persistent conda env cache so snakemake builds tool envs only once ---
    mkdir -p "${params.conda_prefix}"
    export SNAKEMAKE_CONDA_PREFIX="${params.conda_prefix}"

    pixi run --manifest-path ${params.dafqc_repo}/pixi.toml \\
        snakemake \\
        -s ${params.dafqc_repo}/workflow/Snakefile \\
        --configfile config.yaml \\
        --software-deployment-method conda \\
        --conda-prefix "${params.conda_prefix}" \\
        --cores ${task.cpus} \\
        -k
    """

    stub:
    """
    mkdir -p results/${sample}/align results/${sample}/qc
    touch results/${sample}/align/${sample}.decorated.reads.bam
    touch results/${sample}/align/${sample}.decorated.reads.bam.bai
    touch results/${sample}/qc/${sample}.dashboard.html
    """
}
