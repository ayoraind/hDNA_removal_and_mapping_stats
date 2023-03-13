## Workflow to remove Human DNA contaminants from microbial reads.
### Usage

```

=======================================================================
 HUMAN DNA REMOVAL AND MAPPING STATISTICS TAPIR Pipeline version 1.0dev
=======================================================================
 The typical command for running the pipeline is as follows:
        nextflow run main.nf --reads "PathToReadFile(s)" --output_dir "PathToOutputDir" --reference_fasta "PathToRefFasta" 

        Mandatory arguments:
         --reads                        Query fastqz file of sequences you wish to supply as input (e.g., "/MIGE/01_DATA/01_FASTQ/T055-8-*.fastq.gz")
         --reference_fasta              fasta file to be used as reference (e.g., /path/to/GRCh38.primary_assembly.genome.fa)
         --output_dir                   Output directory to place output (e.g., "/MIGE/01_DATA/03_ASSEMBLY")
         
        Optional arguments:
         --help                         This usage statement.
         --version                      Version statement

```


## Introduction
This pipeline maps reads against the human genome, removes human DNA contaminants from reads, estimates the proportion of reads that align with (and do not align with) the human genome, and calculates a few more descriptive statistics. A certain percentage of this pipeline was adapted from the NF Core's [Minimap2 module](https://github.com/nf-core/modules/tree/master/modules/nf-core/minimap2).  


## Sample command
An example of a command to run this pipeline is:

```
nextflow run main.nf --reads "Sample_files/*.fastq.gz" --output_dir "test2" --reference_fasta "PathToRefFasta"
```

## Word of Note
This is an ongoing project at the Microbial Genome Analysis Group, Institute for Infection Prevention and Hospital Epidemiology, Üniversitätsklinikum, Freiburg. The project is funded by BMBF, Germany, and is led by [Dr. Sandra Reuter](https://www.uniklinik-freiburg.de/iuk-en/infection-prevention-and-hospital-epidemiology/research-group-reuter.html).


## Authors and acknowledgment
The TAPIR (Track Acquisition of Pathogens In Real-time) team.
