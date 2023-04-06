process MINIMAP2_INDEX {
    
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


process MINIMAP2_SAM {
    tag "$meta"
    memory { 4.GB * task.attempt }
    errorStrategy { task.attempt <= 5 ? "retry" : "finish" }
    maxRetries 5

    input:
    tuple val(meta), path(reads), path(index)

    output:
    tuple val(meta), path("*.sam"), emit: sam_ch
    path "versions.yml" , emit: versions_ch

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    minimap2 -ax map-ont -t $task.cpus $index $reads  > ${meta}.sam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$(minimap2 --version 2>&1)
    END_VERSIONS
    """
}

process SAM_SORT_AND_INDEX {
    tag "$meta"

    input:
    tuple val(meta), path(sam)

    output:
    tuple val(meta), path("*.bam"), emit: bam_ch
    path "versions.yml" , emit: versions_ch

    when:
    task.ext.when == null || task.ext.when

    script:
    """

    samtools sort $sam > ${meta}.alignment.sorted.bam

    samtools index -b ${meta}.alignment.sorted.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}

process EXTRACT_MICROBIAL_READS {
    publishDir "${params.output_dir}", mode:'copy'
    tag "extract microbial DNA from $sample_id"
    
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


process BAM_STATISTICS {
    publishDir "${params.output_dir}", mode:'copy'
    tag "$sample_id"
    
    
    input:
    tuple val(sample_id), path(bam)
    

    output:
    path("*.mappedstats.txt"), emit: bam_stats_ch

    
    script:
    """ 
    bam_stats.sh ${sample_id}

    """
}


process COMBINE_BAM_STATISTICS {
    publishDir "${params.output_dir}", mode:'copy'
    tag { 'combine bam statistics files'} 
    
    
    input:
    path(bam_statistics_files)
    

    output:
    path("combined_bam_stats.txt"), emit: bam_comb_stats_ch

    
    script:
    """ 
    BAM_STATISTICS_FILES=(${bam_statistics_files})
    
    for index in \${!BAM_STATISTICS_FILES[@]}; do
    BAM_STATISTICS_FILE=\${BAM_STATISTICS_FILES[\$index]}
    
    # add header line if first file
    if [[ \$index -eq 0 ]]; then
      echo "\$(head -1 \${BAM_STATISTICS_FILE})" >> combined_bam_stats.txt
    fi
    echo "\$(awk 'FNR==2 {print}' \${BAM_STATISTICS_FILE})" >> combined_bam_stats.txt
    done

    """
}
