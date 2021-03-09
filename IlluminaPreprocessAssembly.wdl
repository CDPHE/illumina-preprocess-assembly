version 1.0

workflow IlluminaPreprocessAssembly {

    input {
        String    sample_id
        File    fastq_1
        File    fastq_2
        File    V3arctic
        File    adapters_and_contaminants
        File    covid_genome
        File    covid_gff
        String    out_dir
    }

    call seqyclean {
        input:
            contam = adapters_and_contaminants,
            sample_id = sample_id,
            fastq_1 = fastq_1,
            fastq_2 = fastq_2
    }
    
    call fastqc as fastqc_raw {
        input:
           fastq_1 = fastq_1,
           fastq_2 = fastq_2
    }

    call fastqc as fastqc_cleaned {
        input:
            fastq_1 = seqyclean.cleaned_1,
            fastq_2 = seqyclean.cleaned_2
    }

    call align_reads {
        input:
            sample_id = sample_id,
            ref = covid_genome,
            fastq_1 = seqyclean.cleaned_1,
            fastq_2 = seqyclean.cleaned_2
    }

    call ivar_trim {
        input:
            sample_id = sample_id,
            primers = V3arctic,
            bam = align_reads.out_bam
    }

    call ivar_var {
        input:
            sample_id = sample_id,
            ref = covid_genome,
            gff = covid_gff,
            bam = ivar_trim.trimsort_bam
    }

    call ivar_consensus {
        input:
            sample_id = sample_id,
            ref = covid_genome,
            bam = ivar_trim.trimsort_bam
    }

    call bam_stats {
        input:
            sample_id = sample_id,
            bam = ivar_trim.trimsort_bam
    }

    call pangolin {
        input:
            sample_id = sample_id,
            fasta = ivar_consensus.consensus_out
    }
    
    call nextclade_one_sample {
        input:
             genome_fasta = ivar_consensus.consensus_out
    }
    
    call transfer_outputs {
        input:
            out_dir = out_dir,
            data = [seqyclean.cleaned_1, seqyclean.cleaned_2, seqyclean.seqyclean_summary, fastqc_raw.fastqc1_html, fastqc_raw.fastqc1_zip, fastqc_raw.fastqc2_html, fastqc_raw.fastqc2_zip, fastqc_cleaned.fastqc1_html, fastqc_cleaned.fastqc1_zip, fastqc_cleaned.fastqc2_html, fastqc_cleaned.fastqc2_zip, align_reads.out_bam, align_reads.out_bamindex, ivar_trim.trim_bam, ivar_trim.trimsort_bam, ivar_trim.trimsort_bamindex, ivar_var.var_out, ivar_consensus.consensus_out, bam_stats.flagstat_out, bam_stats.stats_out, bam_stats.covhist_out, bam_stats.cov_out, pangolin.lineage, nextclade_one_sample.nextclade_json, nextclade_one_sample.auspice_json, nextclade_one_sample.nextclade_csv]
    }

    output {
        File filtered_reads_1 = seqyclean.cleaned_1
        File filtered_reads_2 = seqyclean.cleaned_2
        File seqyclean_summary = seqyclean.seqyclean_summary
        File fastqc_raw1_html = fastqc_raw.fastqc1_html
        File fastqc_raw1_zip = fastqc_raw.fastqc1_zip
        File fastqc_raw2_html = fastqc_raw.fastqc2_html
        File fastqc_raw2_zip = fastqc_raw.fastqc2_zip
        File fastqc_clean1_html = fastqc_cleaned.fastqc1_html
        File fastqc_clean1_zip = fastqc_cleaned.fastqc1_zip
        File fastqc_clean2_html = fastqc_cleaned.fastqc2_html
        File fastqc_clean2_zip = fastqc_cleaned.fastqc2_zip
        File out_bam = align_reads.out_bam
        File out_bamindex = align_reads.out_bamindex
        File trim_bam = ivar_trim.trim_bam
        File trimsort_bam = ivar_trim.trimsort_bam
        File trimsort_bamindex = ivar_trim.trimsort_bamindex
        File variants = ivar_var.var_out
        File consensus = ivar_consensus.consensus_out
        File flagstat_out = bam_stats.flagstat_out
        File stats_out = bam_stats.stats_out
        File covhist_out = bam_stats.covhist_out
        File cov_out = bam_stats.cov_out
        String pangolin_version = pangolin.pangolin_version
        File pangolin_lineage = pangolin.lineage
        String nextclade_version = nextclade_one_sample.nextclade_version
        File nextclade_json = nextclade_one_sample.nextclade_json
        File auspice_json = nextclade_one_sample.auspice_json
        File nextclade_csv = nextclade_one_sample.nextclade_csv
        String transfer_date = transfer_outputs.transfer_date
    }
}

task seqyclean {
    input {
        File contam
        String sample_id
        File fastq_1
        File fastq_2
    }

    command {

        seqyclean -minlen 70 -qual 30 30 -gz -1 ${fastq_1} -2 ${fastq_2} -c ${contam} -o ${sample_id}_clean

    }

    output {

        File cleaned_1 = "${sample_id}_clean_PE1.fastq.gz"
        File cleaned_2 = "${sample_id}_clean_PE2.fastq.gz"
        File seqyclean_summary = "${sample_id}_clean_SummaryStatistics.tsv"

    }

    runtime {
        cpu:    2
        memory:    "4 GiB"
        disks:    "local-disk 1 HDD"
        bootDiskSizeGb:    10
        preemptible:    0
        maxRetries:    0
        docker:    "staphb/seqyclean:1.10.09"
    }
}

task fastqc {
    input {

        File fastq_1
        File fastq_2
    }
    
    String fastq1_name = basename(basename(basename(fastq_1, ".gz"), ".fastq"), ".fq")
    String fastq2_name = basename(basename(basename(fastq_2, ".gz"), ".fastq"), ".fq")

    command {

        fastqc --outdir $PWD ${fastq_1} ${fastq_2}

    }

    output {

        File fastqc1_html = "${fastq1_name}_fastqc.html"
        File fastqc1_zip = "${fastq1_name}_fastqc.zip"
        File fastqc2_html = "${fastq2_name}_fastqc.html"
        File fastqc2_zip = "${fastq2_name}_fastqc.zip"

    }

    runtime {
        cpu:    1
        memory:    "2 GiB"
        disks:    "local-disk 1 HDD"
        bootDiskSizeGb:    10
        preemptible:    0
        maxRetries:    0
        docker:    "staphb/fastqc:0.11.9"
    }
}

task align_reads {

    input {

        File fastq_1
        File fastq_2
        File ref
        String sample_id
    }

    command {

        bwa index -p reference.fasta -a is ${ref}
        bwa mem -t 2 reference.fasta ${fastq_1} ${fastq_2} | \
        samtools sort | \
        samtools view -u -h -F 4 -o ./${sample_id}_aln.sorted.bam
        samtools index ./${sample_id}_aln.sorted.bam

    }

    output {

        File out_bam = "${sample_id}_aln.sorted.bam"
        File out_bamindex = "${sample_id}_aln.sorted.bam.bai"

    }

    runtime {
        cpu:    2
        memory:    "8 GiB"
        disks:    "local-disk 1 HDD"
        bootDiskSizeGb:    10
        preemptible:    0
        maxRetries:    0
        docker:    "broadinstitute/viral-core:latest"
    }
}

task ivar_trim {

    input {

        File primers
        File bam
        String sample_id
    }

    command {

        ivar trim -e -i ${bam} -b ${primers} -p ${sample_id}_trim.bam
        samtools sort ${sample_id}_trim.bam -o ${sample_id}_trim.sort.bam
        samtools index ${sample_id}_trim.sort.bam

    }

    output {

        File trim_bam = "${sample_id}_trim.bam"
        File trimsort_bam = "${sample_id}_trim.sort.bam"
        File trimsort_bamindex = "${sample_id}_trim.sort.bam.bai"

    }

    runtime {
        cpu:    2
        memory:    "8 GiB"
        disks:    "local-disk 1 HDD"
        bootDiskSizeGb:    10
        preemptible:    0
        maxRetries:    0
        docker:    "andersenlabapps/ivar:1.3.1"
    }
}

task ivar_var {

    input {

        String sample_id
        File ref
        File gff
        File bam
    }

    command {

        samtools faidx ${ref}
        samtools mpileup -A -aa -d 600000 -B -Q 20 -q 20 -f ${ref} ${bam} | \
        ivar variants -p ${sample_id}_variants -q 20 -t 0.6 -m 10 -r ${ref} -g ${gff}

    }

    output {

        File var_out = "${sample_id}_variants.tsv"

    }

    runtime {
        cpu:    2
        memory:    "8 GiB"
        disks:    "local-disk 1 HDD"
        bootDiskSizeGb:    10
        preemptible:    0
        maxRetries:    0
        docker:    "andersenlabapps/ivar:1.3.1"
    }
}

task ivar_consensus {

    input {

        String sample_id
        File ref
        File bam
    }

    command {

        samtools faidx ${ref}
        samtools mpileup -A -aa -d 600000 -B -Q 20 -q 20 -f ${ref} ${bam} | \
        ivar consensus -p ${sample_id}_consensus -q 20 -t 0.6 -m 10

    }

    output {

        File consensus_out = "${sample_id}_consensus.fa"

    }

    runtime {
        cpu:    2
        memory:    "8 GiB"
        disks:    "local-disk 1 HDD"
        bootDiskSizeGb:    10
        preemptible:    0
        maxRetries:    0
        docker:    "andersenlabapps/ivar:1.3.1"
    }
}

task bam_stats {

    input {

        String sample_id
        File bam
    }

    command {

        samtools flagstat ${bam} > ${sample_id}_flagstat.txt
        samtools stats ${bam} > ${sample_id}_stats.txt
        samtools coverage -m -o ${sample_id}_coverage_hist.txt ${bam}
        samtools coverage -o ${sample_id}_coverage.txt ${bam}

    }

    output {

        File flagstat_out  = "${sample_id}_flagstat.txt"
        File stats_out  = "${sample_id}_stats.txt"
        File covhist_out  = "${sample_id}_coverage_hist.txt"
        File cov_out  = "${sample_id}_coverage.txt"

    }

    runtime {
        cpu:    2
        memory:    "8 GiB"
        disks:    "local-disk 1 HDD"
        bootDiskSizeGb:    10
        preemptible:    0
        maxRetries:    0
        docker:    "staphb/samtools:1.10"
    }
}

task pangolin {

    input {

        String sample_id
        File fasta
    }

    command {

        pangolin --version > VERSION
        pangolin --outfile ${sample_id}_lineage_report.csv ${fasta}

    }

    output {

        String pangolin_version = read_string("VERSION")
        File lineage = "${sample_id}_lineage_report.csv"

    }

    runtime {
        cpu:    2
        memory:    "8 GiB"
        disks:    "local-disk 1 HDD"
        bootDiskSizeGb:    10
        preemptible:    0
        maxRetries:    0
        docker:    "staphb/pangolin"
    }
}

task nextclade_one_sample {

    input {
        File genome_fasta
    }
    
    String basename = basename(genome_fasta, ".fa")
    
    command {
        nextclade --version > VERSION
        nextclade --input-fasta "${genome_fasta}" --output-json "${basename}".nextclade.json --output-csv "${basename}".nextclade.csv --output-tree "${basename}".nextclade.auspice.json
    }
    
    output {
        String nextclade_version = read_string("VERSION")
        File nextclade_json = "${basename}.nextclade.json"
        File auspice_json = "${basename}.nextclade.auspice.json"
        File nextclade_csv = "${basename}.nextclade.csv"
    }
    
    runtime {
        docker: "nextstrain/nextclade:0.13.0"
        memory: "3 GB"
        cpu: 2
        disks: "local-disk 50 HDD"
        dx_instance_type: "mem1_ssd1_v2_x2"
    }
}

task transfer_outputs {
    input {
        String out_dir
        Array[File] data
    }
    
    String outdir = sub(out_dir, "/$", "")
    
    parameter_meta {
        data: {
            localization_optional: true
        }
    }

    command <<<
        gsutil -m cp ~{sep=' ' data} ~{outdir}
        transferdate=`date`
        echo $transferdate | tee TRANSFERDATE
    >>>

    output {
        String transfer_date = read_string("TRANSFERDATE")
    }

    runtime {
        docker: "theiagen/utility:1.0"
        memory: "1 GB"
        cpu: 1
        disks: "local-disk 10 SSD"
    }
}
