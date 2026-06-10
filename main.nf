#!/usr/bin/env nextflow

/*
 * nf-dafseq : DAF-seq processing pipeline
 *   1) MINIMAP_ALIGN   minimap2 + read group        -> {sample}.bam
 *   2) DAFQC_SMK       wrap DAF-QC-SMK (Snakemake)  -> {sample}.decorated.reads.bam + QC dashboard
 *   3) MARK_DUPLICATES clustering combinatorial dedup-> {sample}.mkdup.bam
 *   4) PHASE_READS     phase by known SNP/deletion  -> hap1/hap2/non_phased bigWigs (optional)
 */

nextflow.enable.dsl = 2

include { MINIMAP_ALIGN   } from './modules/local/minimap_align'
include { DAFQC_SMK       } from './modules/local/dafqc_smk'
include { MARK_DUPLICATES } from './modules/local/mark_duplicates'
include { PHASE_READS     } from './modules/local/phase_reads'

workflow {

    // Reference (skip existence check during -stub so it can run with no data)
    ch_ref = file(params.ref, checkIfExists: !workflow.stubRun)

    // Parse the TSV samplesheet: sample, fastq, regions, phase_type, phase_value
    ch_samples = Channel
        .fromPath(params.input, checkIfExists: true)
        .splitCsv(sep: '\t', header: true, strip: true)
        .map { row ->
            def sample  = row.sample?.trim()
            def fastq   = row.fastq?.trim()
            def regions = row.regions?.trim()
            def ptype   = (row.phase_type  ?: '').trim().toLowerCase()
            def pval    = (row.phase_value ?: '').trim()
            // treat common "no value" placeholders as no-phasing
            if( ptype in ['none','na','n/a','.','-'] ) { ptype = ''; pval = '' }
            if( pval  in ['none','na','n/a','.','-'] ) { pval = '' }
            if( !sample || !fastq || !regions )
                error "Samplesheet row is missing required sample/fastq/regions: ${row}"
            if( ptype && !(ptype in ['snp','del']) )
                error "phase_type for '${sample}' must be 'snp', 'del', or empty (got '${ptype}')"
            if( ptype && !pval )
                error "phase_type '${ptype}' set for '${sample}' but phase_value is empty"
            tuple(sample, file(fastq, checkIfExists: !workflow.stubRun), regions, ptype, pval)
        }

    // 1 -> 2 -> 3
    MINIMAP_ALIGN( ch_samples, ch_ref )
    DAFQC_SMK( MINIMAP_ALIGN.out.bam )
    MARK_DUPLICATES( DAFQC_SMK.out.decorated )

    // 4 (only samples with a phasing variant defined)
    ch_to_phase = MARK_DUPLICATES.out.mkdup.filter { it[3] in ['snp','del'] }
    PHASE_READS( ch_to_phase )
}

workflow.onComplete {
    log.info ( workflow.success
        ? "\nnf-dafseq finished. Results in: ${params.outdir}\n"
        : "\nnf-dafseq failed. See ${params.outdir}/pipeline_info/ and .nextflow.log\n" )
}
