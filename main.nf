#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// include non-process modules
include { help_message; version_message; complete_message; error_message; pipeline_start_message } from './modules/messages.nf'
include { default_params; check_params } from './modules/params_parser.nf'
include { help_or_version } from './modules/params_utilities.nf'

version = '1.0dev'

// setup default params
default_params = default_params()
// merge defaults with user params
merged_params = default_params + params

// help and version messages
help_or_version(merged_params, version)
final_params = check_params(merged_params)
// starting pipeline
pipeline_start_message(version, final_params)

// include processes
include { MINIMAP2_INDEX; MINIMAP2_ALIGN; EXTRACT_MICROBIAL_READS; BAM_STATISTICS;  COMBINE_BAM_STATISTICS } from './modules/processes.nf' addParams(final_params)

workflow  {
         fasta_ch = channel
                          .fromPath( final_params.reference_fasta, checkIfExists: true )
                          .map { file -> tuple(file.baseName, file) }
         align_ch = channel
                          .fromPath( final_params.reads, checkIfExists: true )
                          .map { file -> tuple(file.simpleName, file) }


         MINIMAP2_INDEX(fasta_ch)

         MINIMAP2_ALIGN(align_ch, MINIMAP2_INDEX.out.index_ch)

         EXTRACT_MICROBIAL_READS(MINIMAP2_ALIGN.out.bam_ch)

         BAM_STATISTICS(MINIMAP2_ALIGN.out.bam_ch)

         collected_bam_statistics_ch = BAM_STATISTICS.out.bam_stats_ch.collect( sort: {a, b -> a[0].getBaseName() <=> b[0].getBaseName()} )

         COMBINE_BAM_STATISTICS(collected_bam_statistics_ch)
}

workflow.onComplete {
    complete_message(final_params, workflow, version)
}

workflow.onError {
    error_message(workflow)
}
