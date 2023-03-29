params.reads = ""
params.reference = ""
params.output = ""
params.reads_rna = true

process MINIMAP2_NANO {
  
  
  publishDir "${params.output}/${name}/minimap2", mode: 'copy', pattern: "*.contamination.sorted.bam*"

  input: 
    tuple val(name), path(reads)
    path db

  output:
    tuple val(name), path ('idxstats.tsv'), emit: idxstats
    tuple val(name), val('clean'), path('*clean.fastq'), emit: cleaned_reads
    tuple val(name), val('contamination'), path('*contamination.fastq'), emit: contaminated_reads
    path '*.contamination.sorted.bam*'

  script:
  """
  PARAMS="-ax map-ont"
  if [[ ${params.reads_rna} != 'false' ]]; then
    PARAMS="-ax splice -k14"
  fi
  minimap2 \$PARAMS -N 5 --secondary=no -t ${task.cpus} -o ${name}.sam ${db} ${reads}
  samtools fastq -f 4 -0 ${reads.baseName}.clean.fastq ${name}.sam
  samtools fastq -F 4 -0 ${reads.baseName}.contamination.fastq ${name}.sam
  samtools view -b -F 2052 ${name}.sam | samtools sort -o ${name}.contamination.sorted.bam --threads ${task.cpus}
  samtools index ${name}.contamination.sorted.bam
  samtools idxstats  ${name}.contamination.sorted.bam > idxstats.tsv
  """
}

workflow {
         reads_ch = channel
                          .fromPath( params.reads, checkIfExists: true )
                          .map { file -> tuple(file.simpleName, file) }
         ref_ch = channel.fromPath( params.reference )

    MINIMAP2_NANO(reads_ch, ref_ch)
    

}
