process qc_bbduk_amplicon {
	label 'bbduk'
	publishDir path: params.output_dir, mode: params.publish_mode, pattern: "${sample.id}/${sample.id}.*{txt,lhist}" //bbduk_stats.txt"

    input:
    tuple val(sample), path(reads)
	path(adapters)

    output:
    tuple val(sample), path("${sample.id}/${sample.id}_R*.fastq.gz"), emit: reads
    tuple val(sample), path("${sample.id}/${sample.id}_O.fastq.gz"), emit: orphans, optional: true
    path("${sample.id}/${sample.id}.*bbduk_stats.txt")
    path("${sample.id}/${sample.id}*lhist")

    script:
    def maxmem = task.memory.toString().replace(/ GB/, "g")

	def read1_p5 = "in1=${sample.id}_R1.fastq.gz out1=${sample.id}_R1.p5.fastq.gz"
	def read1_p3 = "in1=${sample.id}_R1.p5.fastq.gz out1=${sample.id}/${sample.id}_R1.fastq.gz"

	def read2_p5 = ""
	def read2_p3 = ""
	def read2_hist = ""

	if (sample.is_paired) {
		read2_p5 = "in2=${sample.id}_R2.fastq.gz out2=${sample.id}_R2.p5.fastq.gz"
		read2_p3 = "in2=${sample.id}_R2.p5.fastq.gz out2=${sample.id}/${sample.id}_R2.fastq.gz"
		read2_hist = "bbduk.sh -Xmx${maxmem} in1=${sample.id}/${sample.id}_R2.fastq.gz lhist=${sample.id}/${sample.id}_R2.post_lhist"
	}

	if (params.primers) {
		trim_params = "literal=${params.primers} minlength=${params.qc_minlen}"
	} else {
		trim_params = "ref=${adapters} minlength=${params.qc_minlen}"
	}

    """
    mkdir -p ${sample.id}
    bbduk.sh -Xmx${maxmem} t=${task.cpus} ordered=t ${params.p5_primer_params} ${trim_params} stats=${sample.id}/${sample.id}.fwd_bbduk_stats.txt ${read1_p5} ${read2_p5} lhist=${sample.id}/${sample.id}.p5_lhist
    bbduk.sh -Xmx${maxmem} t=${task.cpus} ordered=t ${params.p3_primer_params} ${trim_params} stats=${sample.id}/${sample.id}.rev_bbduk_stats.txt ${read1_p3} ${read2_p3} lhist=${sample.id}/${sample.id}.p3_lhist
	bbduk.sh -Xmx${maxmem} t=${task.cpus} in1=${sample.id}/${sample.id}_R1.fastq.gz \$(echo ${read2_p3} | cut -f 1,3 -d =) lhist=${sample.id}/${sample.id}.post_lhist

	bbduk.sh -Xmx${maxmem} in1=${sample.id}/${sample.id}_R1.fastq.gz lhist=${sample.id}/${sample.id}_R1.post_lhist
	${read2_hist}

    """
	// bbduk.sh -Xmx${maxmem} t=${task.cpus} in1=${sample.id}_R1.fastq.gz \$(echo ${read2_p5} | cut -f 1 -d " ") lhist=${sample.id}/${sample.id}.pre_lhist
}

process qc_bbduk_stepwise_amplicon {
	label 'bbduk'
	publishDir path: params.output_dir, mode: params.publish_mode, pattern: "${sample.id}/${sample.id}*txt"

    input:
    tuple val(sample), path(reads)
	path(adapters)

    output:
    tuple val(sample), path("${sample.id}/${sample.id}_R*.fastq.gz"), emit: reads
    tuple val(sample), path("${sample.id}/${sample.id}_O.fastq.gz"), emit: orphans, optional: true
    path("${sample.id}/${sample.id}.*bbduk_stats.txt"), optional: true
    path("${sample.id}/${sample.id}*lhist.txt"), emit: read_lengths, optional: true

    script:
    def maxmem = task.memory.toString().replace(/ GB/, "g")

	if (params.primers) {
		trim_params = "literal=${params.primers} minlength=${params.qc_minlen}"
	} else {
		trim_params = "ref=${adapters} minlength=${params.qc_minlen}"
	}

	def bbduk_call = "bbduk.sh -Xmx${maxmem} t=${task.cpus} ordered=t trd=t"

	ref_p5_r1 = (params.primers) ? "literal=" + params.primers.split(",")[0] : "ref=${adapters}"
	ref_p5_r2 = (params.primers && !params.single_end) ? "literal=" + params.primers.split(",")[1] : "ref=${adapters}"
	ref_p3_r1 = ref_p5_r2
	ref_p3_r2 = ref_p5_r1

	if (params.single_end) {
		"""
	    mkdir -p ${sample.id}
		${bbduk_call} ${trim_params} in1=${sample.id}_R1.fastq.gz out1=${sample.id}/${sample.id}_R1.fastq.gz stats=${sample.id}/${sample.id}.fwd_bbduk_stats.txt lhist=${sample.id}/${sample.id}.p5_lhist.txt
		${bbduk_call} in1=${sample.id}/${sample.id}_R1.fastq.gz lhist=${sample.id}/${sample.id}_R1.post_lhist.txt
		"""
	} else if (params.long_reads) {
		"""
	    mkdir -p ${sample.id}
		${bbduk_call} ${ref_p5_r1} minlength=${params.qc_minlen} ${params.p5_primer_params} in1=${sample.id}_R1.fastq.gz out1=fwd_p5.fastq.gz
		${bbduk_call} ${ref_p5_r2} minlength=${params.qc_minlen} ${params.p5_primer_params} in1=${sample.id}_R2.fastq.gz out1=rev_p5.fastq.gz
		${bbduk_call} ${ref_p3_r1} minlength=${params.qc_minlen} ${params.p3_primer_params} in1=fwd_p5.fastq.gz out1=fwd.fastq.gz
		${bbduk_call} ${ref_p3_r2} minlength=${params.qc_minlen} ${params.p3_primer_params} in1=rev_p5.fastq.gz out1=rev.fastq.gz
		gzip -dc fwd.fastq.gz | awk 'NR%4==1' | sed 's/^@//' | sed 's/\\/1//' | sort > fwd.txt
        gzip -dc rev.fastq.gz | awk 'NR%4==1' | sed 's/^@//' | sed 's/\\/2//' | sort > rev.txt
		join -1 1 -2 1 fwd.txt rev.txt > both.txt
		seqtk subseq fwd.fastq.gz both.txt | gzip -c - > ${sample.id}/${sample.id}_R1.fastq.gz
		seqtk subseq rev.fastq.gz both.txt | gzip -c - > ${sample.id}/${sample.id}_R2.fastq.gz
		${bbduk_call} in1=${sample.id}/${sample.id}_R1.fastq.gz lhist=${sample.id}/${sample.id}_R1.post_lhist.txt
		${bbduk_call} in1=${sample.id}/${sample.id}_R2.fastq.gz lhist=${sample.id}/${sample.id}_R2.post_lhist.txt
		"""
	} else {
		"""
	    mkdir -p ${sample.id}
		${bbduk_call} ${ref_p5_r1} minlength=${params.qc_minlen} ${params.p5_primer_params} in1=${sample.id}_R1.fastq.gz out1=fwd.fastq.gz
		${bbduk_call} ${ref_p5_r2} minlength=${params.qc_minlen} ${params.p5_primer_params} in1=${sample.id}_R2.fastq.gz out1=rev.fastq.gz
		gzip -dc fwd.fastq.gz | awk 'NR%4==1' | sed 's/^@//' | sed 's/\\/1//' | sort > fwd.txt
        gzip -dc rev.fastq.gz | awk 'NR%4==1' | sed 's/^@//' | sed 's/\\/2//' | sort > rev.txt
		join -1 1 -2 1 fwd.txt rev.txt > both.txt
		seqtk subseq fwd.fastq.gz both.txt | gzip -c - > ${sample.id}/${sample.id}_R1.fastq.gz
		seqtk subseq rev.fastq.gz both.txt | gzip -c - > ${sample.id}/${sample.id}_R2.fastq.gz
		${bbduk_call} in1=${sample.id}/${sample.id}_R1.fastq.gz lhist=${sample.id}/${sample.id}_R1.post_lhist.txt
		${bbduk_call} in1=${sample.id}/${sample.id}_R2.fastq.gz lhist=${sample.id}/${sample.id}_R2.post_lhist.txt
		"""
	}
}
