// Step 1: map ONT reads with minimap2, add a read group, sort + index.
// Ports scripts/minimap.sh, parameterized and without the config.tbl side effect.

process MINIMAP_ALIGN {
    tag { sample }
    label 'process_medium'
    publishDir "${params.outdir}/${sample}/align", mode: 'copy'

    container 'dhspence/docker-dafseq:latest'

    input:
    tuple val(sample), path(fastq), val(regions), val(phase_type), val(phase_value)
    path ref

    output:
    tuple val(sample), path("${sample}.bam"), path("${sample}.bam.bai"),
          val(regions), val(phase_type), val(phase_value), emit: bam

    script:
    """
    minimap2 \\
        -ax map-ont \\
        -t ${task.cpus} \\
        -R "@RG\\tID:${sample}\\tSM:${sample}\\tLB:${sample}" \\
        ${ref} \\
        ${fastq} \\
    | samtools view -bS - \\
    | samtools sort -@ ${task.cpus} -o ${sample}.bam
    samtools index ${sample}.bam
    """

    stub:
    """
    touch ${sample}.bam ${sample}.bam.bai
    """
}
