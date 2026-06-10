// Step 3: clustering-based combinatorial dedup on the decorated BAM.
// Calls bin/mark_duplicates.decorated.py (auto-on-PATH via Nextflow bin/).

process MARK_DUPLICATES {
    tag { sample }
    label 'process_medium'
    publishDir "${params.outdir}/${sample}/dedup", mode: 'copy'

    container 'dhspence/docker-dafseq:latest'

    input:
    tuple val(sample), path(decorated_bam), path(decorated_bai),
          val(regions), val(phase_type), val(phase_value)

    output:
    tuple val(sample), path("${sample}.mkdup.bam"), path("${sample}.mkdup.bam.bai"),
          val(phase_type), val(phase_value),                       emit: mkdup
    path "${sample}.clusters.txt",                                 emit: clusters, optional: true
    path "*.clusters.png",                                         emit: plots,    optional: true

    script:
    """
    # regions string (comma-separated chr:start-end) -> BED
    echo "${regions}" | tr ',' '\\n' | awk -F'[:-]' 'NF>=3 {print \$1"\\t"\$2"\\t"\$3}' > regions.bed

    mark_duplicates.decorated.py \\
        -i ${decorated_bam} \\
        -o ${sample} \\
        -b regions.bed

    samtools index ${sample}.mkdup.bam
    """

    stub:
    """
    touch ${sample}.mkdup.bam ${sample}.mkdup.bam.bai ${sample}.clusters.txt
    touch ${sample}.stub_0_1.clusters.png
    """
}
