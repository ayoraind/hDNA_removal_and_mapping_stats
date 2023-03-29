params.reference_fasta = "$projectDir/data/ggal/gut.fasta"

params.reads = "$projectDir/data/ggal/gut_{1,2}.fastq"
params.output_dir = "$PWD/results"
nextflow.enable.dsl=2


process MINIMAP2_INDEX {
    
    input:
    tuple val(meta), path(reference_fasta)

    output:
    path("*.mmi"),       emit: index_ch
    path "versions.yml", emit: versions_ch

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    minimap2 \\
        -t $task.cpus \\
        -d ${reference_fasta.baseName}.mmi \\
        $reference_fasta
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$(minimap2 --version 2>&1)
    END_VERSIONS
    """
}

process MINIMAP2_ALIGN {
    tag "$meta"
    
    publishDir "${params.output_dir}", mode:'copy'
    
    input:
    tuple val(meta), path(reads), path(index)

    output:
    tuple val(meta), path("*.bam"), emit: bam_ch
    path "versions.yml" , emit: versions_ch

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    minimap2 \\
        -ax map-ont \\
        -t $task.cpus \\
        $index \\
        $reads | samtools sort > ${meta}.alignment.sorted.bam
	
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$(minimap2 --version 2>&1)
    END_VERSIONS
    """
}

process EXTRACT_MICROBIAL_READS {
    publishDir "${params.output_dir}", mode:'copy'
    tag "extract microbial DNA from $sample_id"

    // Note: the versions here need to match the versions used in the mulled container below and minimap2/index
    conda "/MIGE/01_DATA/07_TOOLS_AND_SOFTWARE/nextflow_pipelines/hDNA_removal_and_map_stats/mapping_env.yml"

    input:
    tuple val(sample_id), path(bam)


    output:
    tuple val(sample_id), path("*.fastq.gz"), emit: microbial_reads_ch
    path "versions.yml"                     , emit: versions_ch

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${sample_id}"

    """
    samtools view --threads ${task.cpus} -b -f 4 ${prefix}.alignment.sorted.bam > ${prefix}.bam

    samtools sort $args -@ $task.cpus ${prefix}.bam -o ${prefix}.sorted.bam

    samtools index -@ ${task.cpus} $args -b ${prefix}.sorted.bam

    samtools fastq $args --threads ${task.cpus} ${prefix}.sorted.bam > ${prefix}.fastq

    gzip ${prefix}.fastq

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
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
	// MINIMAP2_INDEX.out.index_ch.view()
	
	ch_index  = MINIMAP2_INDEX.out.index_ch
	
	combine_ch = align_ch.combine(ch_index)
	
	
	MINIMAP2_ALIGN ( combine_ch )
	
	collect_align = MINIMAP2_ALIGN.out.bam_ch.collect(){it -> [ it[1] ] }
	// collect_align.view()		
	
	collect_align
	            .flatten()
		    .map { file -> tuple(file.simpleName, file) }
		    .set { flatten_collect_align_ch }
		    
       // flatten_collect_align_ch.view()
		    
    //    EXTRACT_MICROBIAL_READS(flatten_collect_align_ch)		    
}
