/*
 * pipeline input parameters
 */
params.reference_fasta = "$projectDir/data/ggal/gut.fasta"

params.reads = "$projectDir/data/ggal/gut_{1,2}.fastq"
params.output_dir = "$PWD/results"
nextflow.enable.dsl=2


log.info """\
    MAPPING  - TAPIR   P I P E L I N E
    ============================================
    output_dir       : ${params.output_dir}
    """
    .stripIndent()

/*
 * define the `index` process that creates a binary index
 * given the transcriptome file
 */
 
process MINIMAP2_INDEX {
//    publishDir "${params.output_dir}/minimap", mode:'copy'
    
    // Note: the versions here need to match the versions used in minimap2/align
    conda "bioconda::minimap2=2.24"
    
    input:
    tuple val(meta), path(reference_fasta)

    output:
    path("*.mmi"),                  emit: index_ch
    path "versions.yml"           , emit: versions_ch

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    minimap2 \\
        -t $task.cpus \\
        -d ${reference_fasta.baseName}.mmi \\
        $args \\
        $reference_fasta
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$(minimap2 --version 2>&1)
    END_VERSIONS
    """
}


process MINIMAP2_ALIGN {
    publishDir "${params.output_dir}", mode:'copy'
    tag "minimap alignment on $sample_id"
    
    // Note: the versions here need to match the versions used in the mulled container below and minimap2/index
    conda "bioconda::minimap2=2.24 bioconda::samtools=1.14"
    
    input:
    tuple val(sample_id), path(reads)
    each reference
    
    

    output:
    tuple val(sample_id), path("*.bam"), emit: bam_ch
    path "versions.yml"                , emit: versions_ch

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${sample_id}"
    
    """	 
    minimap2 $args -ax map-ont -t $task.cpus ${reference} ${reads} | samtools sort > ${prefix}_alignment.sorted.bam
    
    samtools index -b ${prefix}_alignment.sorted.bam	 
	 
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$(minimap2 --version 2>&1)
    END_VERSIONS
    """
}

workflow  {
         fasta_ch = channel
                          .fromPath( params.reference_fasta, checkIfExists: true )
			  .map { file -> tuple(file.baseName, file) }
	 align_ch = channel
                          .fromPath( params.reads, checkIfExists: true )
			  .map { file -> tuple(file.simpleName, file) }		
	 
       
	     MINIMAP2_INDEX(fasta_ch)
	     align_ch.view()
//	     MINIMAP2_INDEX.out.index_ch.first().view()
	  mmi_value_ch = MINIMAP2_INDEX.out.index_ch.view()     
	  
	     MINIMAP2_ALIGN(align_ch, MINIMAP2_INDEX.out.index_ch)
//             MINIMAP2_ALIGN(align_ch, mmi_value_ch)
	     MINIMAP2_ALIGN.out.bam_ch.view()
	
	   
}

