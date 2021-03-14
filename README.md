# illumina-preprocess-assembly:

preprocessing and assembly workflows for Illumina MiSeq data

## IlluminaPreprocessAssembly.wdl
The IlluminaPreprocessAssembly.wdl workflow was developed for the preprocessing and assembly of Illumina MiSeq sequencing of SARS-CoV-2 to be run on the GCP Terra platform. It takes paired-end fastq files as input.

The columns for the input data table in Terra should be arranged as:

1. entity:sample_id (column of sample names/ids). If there is more than one data table in the Terra Workspace, add a number after the word sample (e.g. entity:sample2_id).
2. fastq_1 - This is the R1 fastq file.
3. fastq_2 - This is the R2 fastq file.

Also needed as input (these should be saved as Workspace data and included in Terra input field for the workflow): 
1. covid_genome - the path to the google bucket directory containing the SARS-CoV-2 reference genome fasta
2. covid_annotation - the path to the google bucket directory containing the SARS-CoV-2 reference genome gff annotation file
3. adapter_and_contaminants - the path to the google bucket directory containing a fasta file of adapter sequences as potential contaminants (e.g. PhiX)
4. V3arctic - the path to the google bucket directory containing a bed file with the primers used for amplicon sequencing

The IlluminaPreprocessAssembly.wdl workflow will:

1. Use Seqyclean to quality filter and trim raw fastq files
2. Run FastQC on both the raw and cleaned reads
3. Align reads to the reference genome using bwa and then sort the bam by coordinates using Samtools
4. Use iVar trim to trim primer regions and then sort the trimmed bam by coordinates using Samtools
5. Use iVar variants to call variants from the trimmed and sorted bam
6. Use iVar consensus to call the consensus genome sequence from the trimmed and sorted bam
7. Use Samtools flagstat, stats, and coverage to output statistics from the bam
8. Renames consensus sequences to CO-CDPHE-{sample_id}
