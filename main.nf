#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { qc_bbduk } from "./modules/nevermore/qc/bbduk"
include { fastqc } from "./modules/nevermore/qc/fastqc"
include { classify_sample } from "./modules/nevermore/functions"


//def helpMessage() {
//	log.info """
//
//	This is the gaga2 (figaro/dada2) 16S amplicon sequencing processing pipeline of the Zeller lab @ EMBL.
//
//	Usage:
//
//	The typical command for running the pipeline is as follows:
//
//        nextflow run gaga2.nf -c <config_file> --input_dir </path/to/input_dir> --output_dir </path/to/output_dir> ...
//
//    Mandatory arguments:
//
//            --input_dir               path to input read files (needs to be absolute, 1 subdirectory per sample)
//                                      reads need to be Illumina paired end amplicon seq (e.g. MiSeq)
//            --output_dir              path to output directory (absolute)
//            --amplicon_length         expected amplicon length
//            --left_primer             length of left primer
//            --right_primer            length of right primer
//
//    Optional arguments:
//
//            --min_overlap             Minimum read pair overlap [bp] (default=20)
//            --nthreads                Number of threads used for dada2 (default=8)
//            -work-dir, -w             Path to working directory
//            --help                    Show this help.
//
//	""".stripIndent()
//}
//
//if (params.help) {
//	helpMessage()
//	exit 0
//}


process figaro {
	publishDir "${params.output_dir}", mode: params.publish_mode

	input:
	path input_reads
	val is_paired_end

	output:
	path("figaro/trimParameters.json"), emit: trim_params
	path("figaro/*.png")

	script:
	def paired_params = (is_paired_end == true) ? "-r ${params.right_primer} -m ${params.min_overlap}" : ""

	"""
	figaro -i . -o figaro/ -a ${params.amplicon_length} -f ${params.left_primer} ${paired_params}
	"""
}


process extract_trim_parameters {
	input:
	path(trim_params)

	output:
	path("trim_params.txt"), emit: trim_params

	script:
	"""
	python $projectDir/scripts/trim_params.py $trim_params > trim_params.txt
	"""
}


process dada2_preprocess {
	label 'dada2'
	publishDir "${params.output_dir}", mode: params.publish_mode

	input:
	path input_reads
	path trim_params
	path dada2_script
	val is_paired_end

	output:
	path("read_quality.pdf")
	path("read_quality_postqc.pdf")
	path("filter_trim_table.tsv")
	path("filter_trim_table.final.tsv"), emit: trim_table
	path("dada2_preprocess.log")
	path("dada2_preprocess/*.{fastq,fq,fastq.gz,fq.gz}"), emit: filtered_reads

	script:
	"""
	mkdir -p dada2_in/
	for f in \$(find . -maxdepth 1 -type l); do ln -s ../\$f dada2_in/; done
	rm dada2_in/trim_params.txt dada2_in/*.R
	Rscript --vanilla ${dada2_script} dada2_in/ dada2_preprocess/ \$(cat trim_params.txt) $task.cpus > dada2_preprocess.log
	"""
}


process dada2_analysis {
	label 'dada2'
	publishDir "${params.output_dir}", mode: params.publish_mode

	input:
	path input_reads
	path filter_trim_table
	path dada2_script
	val is_paired_end

	output:
	path("dada2_analysis.log")
	path("error_model.pdf")
	path("summary_table.tsv")
	path("result.RData")
	path("dada2_figures.pdf")
	path("ASVs.tsv")
	path("asv_table.tsv")

	script:
	"""
	mkdir -p dada2_in/
	for f in \$(find . -maxdepth 1 -type l); do ln -s ../\$f dada2_in/; done
	rm dada2_in/*.R dada2_in/filter_trim_table.final.tsv
	Rscript --vanilla ${dada2_script} dada2_in/ dada2_analysis/ filter_trim_table.final.tsv $task.cpus > dada2_analysis.log
	"""
}


process assess_read_length_distribution {
	input:
	path(fastq_reports)

	output:
	path("read_length_thresholds.txt"), emit: read_length
	path("READSET_HOMOGENEOUS"), emit: hom_reads_marker, optional: true

	script:
	"""
	python ${projectDir}/scripts/assess_readlengths.py . > read_length_thresholds.txt
	"""
}


process homogenise_readlengths {
	label 'bbduk'

	input:
	tuple val(sample), path(reads)
	path(read_lengths)

	output:
	path("${sample.id}/*.{fastq,fq,fastq.gz,fq.gz}"), emit: reads

	script:
	def maxmem = task.memory.toString().replace(/ GB/, "g")

	if (sample.is_paired) {
		"""
		mkdir -p ${sample.id}
		r1len=\$(head -n 1 ${read_lengths} | cut -f 1)
		r2len=\$(head -n 1 ${read_lengths} | cut -f 4)
		bbduk.sh -Xmx${maxmem} t=${task.cpus} ordered=t minlength=\$((r1len-1)) ftr=\$((r1len-1)) stats=${sample.id}/${sample.id}.homr_stats_1.txt in=${sample.id}_R1.fastq.gz out=${sample.id}/${sample.id}_R1.fastq.gz
		bbduk.sh -Xmx${maxmem} t=${task.cpus} ordered=t minlength=\$((r2len-1)) ftr=\$((r2len-1)) stats=${sample.id}/${sample.id}.homr_stats_2.txt in=${sample.id}_R2.fastq.gz out=${sample.id}/${sample.id}_R2.fastq.gz
		"""
	} else {
		"""
		mkdir -p ${sample.id}
		r1len=\$(head -n 1 ${read_lengths} | cut -f 1)
		bbduk.sh -Xmx${maxmem} t=${task.cpus} ordered=t minlength=\$((r1len-1)) ftr=\$((r1len-1)) stats=${sample.id}/${sample.id}.homr_stats_1.txt in=${sample.id}_R1.fastq.gz out=${sample.id}/${sample.id}_R1.fastq.gz
		"""
	}
}


workflow raw_reads_figaro {

	take:
		reads
		run_figaro

	main:
		qc_bbduk(reads, "${projectDir}/assets/adapters.fa", run_figaro)

		fastqc(qc_bbduk.out.reads)
		fastqc_ch = fastqc.out.reports
			.map { sample, report -> return report }
			.collect()

		assess_read_length_distribution(fastqc_ch)
		homogenise_readlengths(qc_bbduk.out.reads, assess_read_length_distribution.out.read_length)

		hom_reads = homogenise_readlengths.out.reads.collect()
		figaro(hom_reads, !params.single_end)
		extract_trim_parameters(figaro.out.trim_params)

	emit:
		reads = hom_reads
		trim_params = extract_trim_parameters.out.trim_params

}


workflow check_for_preprocessing {
	take:
		reads

	main:
		reads.view()
		fastqc(reads)
		assess_read_length_distribution(
			fastqc.out.reports
				.map { sample, report -> return report }
				.collect()
		)

	emit:
		readlen_dist = assess_read_length_distribution.out.read_length
		hom_reads_marker = assess_read_length_distribution.out.hom_reads_marker
}


process prepare_fastqs {
	input:
	tuple val(sample), path(fq)

	output:
	tuple val(sample), path("fastq/${sample.id}/${sample.id}_R*.fastq.gz"), emit: reads

	script:
	if (sample.is_paired) {
		"""
		mkdir -p fastq/${sample.id}
		ln -sf ../../${fq[0]} fastq/${sample.id}/${sample.id}_R1.fastq.gz
		ln -sf ../../${fq[1]} fastq/${sample.id}/${sample.id}_R2.fastq.gz
		"""
	} else {
		"""
		mkdir -p fastq/${sample.id}
		ln -sf ../../${fq[0]} fastq/${sample.id}/${sample.id}_R1.fastq.gz
		"""
	}
}


workflow {
	fastq_ch = Channel
		.fromPath(params.input_dir + "/**_*[12].{fastq,fq,fastq.gz,fq.gz}")
		.map { file ->
			def sample = file.name.replaceAll(/.(fastq|fq)(.gz)?$/, "")
			sample = sample.replaceAll(/_R?[12]$/, "")
			return tuple(sample, file)
		}
		.groupTuple(sort: true)
		.map { classify_sample(it[0], it[1]) }

	prepare_fastqs(fastq_ch)

	if (params.single_end) {
		library_layout = "SINGLE";
		dada2_preprocess_script = "$projectDir/R_scripts/dada2_preprocess_single.R"
		dada2_analysis_script = "$projectDir/R_scripts/dada2_analysis_single.R"
		is_paired_end = false
	} else {
		library_layout = "PAIRED";
		dada2_preprocess_script = "$projectDir/R_scripts/dada2_preprocess_paired.R"
		dada2_analysis_script = "$projectDir/R_scripts/dada2_analysis_paired.R"
		is_paired_end = true
	}

	print library_layout

	trim_params_ch = Channel.empty()
	dada_reads_ch = Channel.empty()

	if (!params.preprocessed) {

		/* check if dataset was preprocessed */

		check_for_preprocessing(prepare_fastqs.out.reads)

		raw_reads_figaro(prepare_fastqs.out.reads, check_for_preprocessing.out.hom_reads_marker)
		trim_params_ch = trim_params_ch
			.concat(raw_reads_figaro.out.trim_params)

		dada_reads_ch = dada_reads_ch.concat(raw_reads_figaro.out.reads)

	}

	trim_params = file("${workDir}/trim_params.txt")
	trim_params.text = "-1 -1\n"

	trim_params_ch = trim_params_ch
		.concat(Channel.fromPath("${workDir}/trim_params.txt"))

	dada_reads_ch = dada_reads_ch.concat(
		prepare_fastqs.out.reads
			.map { sample, reads -> reads }
			.collect()
	)

	dada2_preprocess(dada_reads_ch.first(), trim_params_ch.first(), dada2_preprocess_script, is_paired_end)
	dada2_analysis(dada2_preprocess.out.filtered_reads, dada2_preprocess.out.trim_table, dada2_analysis_script, is_paired_end)
}
