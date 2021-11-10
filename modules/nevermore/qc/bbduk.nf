
process qc_bbduk {
	publishDir path: params.output_dir, mode: params.publish_mode, pattern: "${sample.id}/${sample.id}.bbduk_stats.txt"

    input:
    tuple val(sample), path(reads)
	path(adapters)
	path(run_sentinel)

    output:
    tuple val(sample), path("${sample.id}/${sample.id}_R*.fastq.gz"), emit: reads
    tuple val(sample), path("${sample.id}/${sample.id}_O.fastq.gz"), emit: orphans, optional: true
    path("${sample.id}/${sample.id}.bbduk_stats.txt")

    script:
    def maxmem = task.memory.toString().replace(/ GB/, "g")
    //def read1 = "in1=${sample.id}_R1.fastq.gz out1=${sample.id}/${sample.id}_R1.fastq.gz"
    def read1 = "in1=${sample.id}_R1.fastq.gz out1=${sample.id}/${sample.id}_R1.fastq.gz"
    //read2 = sample.is_paired ? "in2=${sample.id}_R2.fastq.gz out2=${sample.id}/${sample.id}_R2.fastq.gz outs=${sample.id}/${sample.id}_O.fastq.gz" : ""
    read2 = sample.is_paired ? "in2=${sample.id}_R2.fastq.gz out2=${sample.id}/${sample.id}_R2.fastq.gz outs=${sample.id}/${sample.id}_O.fastq.gz" : ""

	if (params.primers) {
		qc_params = params.qc_params_primers
		trim_params = "literal=${params.primers} minlen=${params.qc_minlen}"
	} else {
		qc_params = params.qc_params_adapters
		trim_params = "ref=${adapters} minlen=${params.qc_minlen}"
	}

    """
    mkdir -p ${sample.id}
    bbduk.sh -Xmx${maxmem} t=${task.cpus} ordered=t ${qc_params} ${trim_params} stats=${sample.id}/${sample.id}.bbduk_stats.txt ${read1} ${read2}
    """
}
