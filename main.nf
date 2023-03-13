#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// include definitions
include  { helpMessage; Version } from './modules/messages.nf'

// include processes
include { MINIMAP2_INDEX; MINIMAP2_ALIGN; EXTRACT_MICROBIAL_READS; BAM_STATISTICS;  COMBINE_BAM_STATISTICS } from './modules/processes.nf'

log.info """\
    ==================================================================
    HUMAN DNA REMOVAL AND MAPPING STATISTICS  - TAPIR  P I P E L I N E
    ==================================================================
    output_dir      : ${params.output_dir}
    """
    .stripIndent()


workflow  {
         fasta_ch = channel
                          .fromPath( params.reference_fasta, checkIfExists: true )
                          .map { file -> tuple(file.baseName, file) }
         align_ch = channel
                          .fromPath( params.reads, checkIfExists: true )
                          .map { file -> tuple(file.simpleName, file) }


         MINIMAP2_INDEX(fasta_ch)

         MINIMAP2_ALIGN(align_ch, MINIMAP2_INDEX.out.index_ch)

         EXTRACT_MICROBIAL_READS(MINIMAP2_ALIGN.out.bam_ch)

         BAM_STATISTICS(MINIMAP2_ALIGN.out.bam_ch)

         collected_bam_statistics_ch = BAM_STATISTICS.out.bam_stats_ch.collect( sort: {a, b -> a[0].getBaseName() <=> b[0].getBaseName()} )

         COMBINE_BAM_STATISTICS(collected_bam_statistics_ch)
}
