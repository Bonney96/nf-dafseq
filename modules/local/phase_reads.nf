// Step 4: phase reads by a known SNP or deletion, emit phased BAM + hap1/hap2/non_phased bigWigs.
// Calls bin/phase_reads_del.py. Only runs for samples whose phase_type is 'snp' or 'del'.
// bedGraphToBigWig is provided by `module load labtools`.

process PHASE_READS {
    tag { "${sample}:${phase_type}" }
    label 'process_low'
    publishDir "${params.outdir}/${sample}/phasing", mode: 'copy'

    beforeScript 'module load labtools'

    input:
    tuple val(sample), path(mkdup_bam), path(mkdup_bai), val(phase_type), val(phase_value)

    output:
    tuple val(sample), path("${sample}.filtered_phased.bam"), path("${sample}.filtered_phased.bam.bai"), emit: phased
    path "${sample}.*.bw",                                                                                emit: bigwigs
    path "${sample}.*.bedgraph",                                                                          emit: bedgraphs

    script:
    def phase_flag = phase_type == 'snp' ? '--snp' : '--del'
    def extra      = phase_type == 'del' ? "--del-overlap-frac ${params.del_overlap_frac}" : ''
    """
    phase_reads_del.py \\
        --bam ${mkdup_bam} \\
        ${phase_flag} '${phase_value}' \\
        ${extra} \\
        --results-path . \\
        --sample-prefix ${sample}
    """

    stub:
    """
    touch ${sample}.filtered_phased.bam ${sample}.filtered_phased.bam.bai
    touch ${sample}.non_phased.bw ${sample}.hap1.bw ${sample}.hap2.bw
    touch ${sample}.non_phased.bedgraph ${sample}.hap1.bedgraph ${sample}.hap2.bedgraph
    """
}
