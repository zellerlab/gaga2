process fastqc {
    publishDir params.output_dir, mode: params.publish_mode, pattern: "raw_counts/*.txt"

    input:
    tuple val(sample), path(reads)

    output:
    tuple val(sample), path("fastqc/*/*fastqc_data.txt"), emit: reports
    tuple val(sample), path("raw_counts/${sample.id}.txt"), emit: counts

    script:
    def process_r2 = (sample.is_paired) ? "fastqc -t $task.cpus --extract --outdir=fastqc ${sample.id}_R2.fastq.gz && mv fastqc/${sample.id}_R2_fastqc/fastqc_data.txt fastqc/${sample.id}_R2_fastqc/${sample.id}_R2_fastqc_data.txt" : "";

    """
    mkdir -p fastqc
    mkdir -p raw_counts
    fastqc -t $task.cpus --extract --outdir=fastqc ${sample.id}_R1.fastq.gz && mv fastqc/${sample.id}_R1_fastqc/fastqc_data.txt fastqc/${sample.id}_R1_fastqc/${sample.id}_R1_fastqc_data.txt
    ${process_r2}
    grep "Total Sequences" fastqc/*/*data.txt > seqcount.txt
    echo \$(wc -l seqcount.txt)\$'\t'\$(head -n1 seqcount.txt | cut -f 2) > raw_counts/${sample.id}.txt
    """
}
