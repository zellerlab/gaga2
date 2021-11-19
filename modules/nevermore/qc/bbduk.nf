process qc_bbduk_amplicon {
	label 'bbduk'
	publishDir path: params.output_dir, mode: params.publish_mode, pattern: "${sample.id}/${sample.id}.*bbduk_stats.txt"

    input:
    tuple val(sample), path(reads)
	path(adapters)

    output:
    tuple val(sample), path("${sample.id}/${sample.id}_R*.fastq.gz"), emit: reads
    tuple val(sample), path("${sample.id}/${sample.id}_O.fastq.gz"), emit: orphans, optional: true
    path("${sample.id}/${sample.id}.*bbduk_stats.txt")

    script:
    def maxmem = task.memory.toString().replace(/ GB/, "g")

	def read1_p5 = "in1=${sample.id}_R1.fastq.gz out1=${sample.id}_R1.p5.fastq.gz"
	def read1_p3 = "in1=${sample.id}_R1.p5.fastq.gz out1=${sample.id}/${sample.id}_R1.fastq.gz"

	def read2_p5 = ""
	def read2_p3 = ""

	if (sample.is_paired) {
		read2_p5 = "in2=${sample.id}_R2.fastq.gz out2=${sample.id}_R2.p5.fastq.gz"
		read2_p3 = "in2=${sample.id}_R2.p5.fastq.gz out2=${sample.id}/${sample.id}_R2.fastq.gz"
	}

	if (params.primers) {
		trim_params = "literal=${params.primers} minlength=${params.qc_minlen}"
	} else {
		trim_params = "ref=${adapters} minlength=${params.qc_minlen}"
	}

    """
    mkdir -p ${sample.id}

    bbduk.sh -Xmx${maxmem} t=${task.cpus} ordered=t ${params.p5_primer_params} ${trim_params} stats=${sample.id}/${sample.id}.fwd_bbduk_stats.txt ${read1_p5} ${read2_p5}
    bbduk.sh -Xmx${maxmem} t=${task.cpus} ordered=t ${params.p3_primer_params} ${trim_params} stats=${sample.id}/${sample.id}.rev_bbduk_stats.txt ${read1_p3} ${read2_p3}
    """
}
