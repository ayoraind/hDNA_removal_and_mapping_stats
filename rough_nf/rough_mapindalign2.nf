/*
 * pipeline input parameters
 */
params.reference_fasta = "$projectDir/data/ggal/gut.fasta"
params.bam_format = true

params.reads = "$projectDir/data/ggal/gut_{1,2}.fastq"
params.output_dir = "$PWD/results"


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
    publishDir "${params.output_dir}/minimap", mode:'copy'
    
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
    path reference
    val bam_format
    
    

    output:
    tuple val(sample_id), path("*.paf"), optional: true, emit: paf_ch
    tuple val(sample_id), path("*.bam"), optional: true, emit: bam_ch
    path "versions.yml"           , emit: versions_ch

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${sample_id}"
    def bam_output = bam_format ? "-a | samtools sort | samtools view -@ ${task.cpus} -b -h -o ${prefix}.bam" : "-o ${prefix}.paf"
    """
    minimap2 \\
        $args \\
        -t $task.cpus \\
        "${reference ?: reads}" \\
        "$reads" \\
         $bam_output
	 
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
	 bam_define_ch = channel.of( params.bam_format )	  
        	  
	 MINIMAP2_INDEX(fasta_ch)
	 //align_ch.view()
	 //MINIMAP2_INDEX.out.index_ch.view()     
	  
	 MINIMAP2_ALIGN(align_ch, MINIMAP2_INDEX.out.index_ch, bam_define_ch)
	   
}
