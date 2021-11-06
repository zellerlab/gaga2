#!/usr/bin/env nextflow

nextflow.enable.dsl = 2


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
//
if ( !params.min_overlap ) {
	params.min_overlap = 20
}

if ( !params.left_primer ) {
	params.left_primer = 0
}

if ( !params.right_primer ) {
	params.right_primer = 0
}


process ltrim {
	publishDir "${params.output_dir}", mode: params.publish_mode

	input:
	path input_reads

	output:
	path("ltrim/*.{fastq,fq,fastq.gz,fq.gz}"), emit: ltrimmed_reads

	script:
	"""
	mkdir -p ltrim/
	python ${projectDir}/scripts/ltrim.py . ltrim/
	"""

}

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
	//  -r ${params.right_primer} -m ${params.min_overlap}
}

process extract_trim_parameters {
	input:
	path(trim_params)

	output:
	//tuple val(left), val(right), emit: trim_params
	path("trim_params.txt"), emit: trim_params

	script:
	"""
	python $projectDir/scripts/trim_params.py $trim_params > trim_params.txt
	"""

}

process dada2_preprocess {
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


workflow {
	fastq_ch = Channel
		.fromPath(params.input_dir + "/**_*[12].{fastq,fq,fastq.gz,fq.gz}")
		.map { file ->
			def sample = file.name.replaceAll(/.(fastq|fq)(.gz)?$/, "")
			sample = sample.replaceAll(/_R?[12]$/, "")
			return tuple(sample, file)
		}
		.groupTuple(sort: true)

	if (!params.single_end) {
		library_layout = "PAIRED";
		dada2_preprocess_script = "$projectDir/R_scripts/dada2_preprocess_paired.R"
		dada2_analysis_script = "$projectDir/R_scripts/dada2_analysis_paired.R"
		is_paired_end = true
	} else {
		library_layout = "SINGLE";
		dada2_preprocess_script = "$projectDir/R_scripts/dada2_preprocess_single.R"
		dada2_analysis_script = "$projectDir/R_scripts/dada2_analysis_single.R"
		is_paired_end = false
	}

	print library_layout

	fastq_ch.collect().view()

	files_only_ch = fastq_ch
		.map { sample, files -> return files }
		.collect()

	files_only_ch.view()

	ltrim(files_only_ch)

	figaro(ltrim.out.ltrimmed_reads, is_paired_end)
	extract_trim_parameters(figaro.out.trim_params)

	dada2_preprocess(ltrim.out.ltrimmed_reads, extract_trim_parameters.out.trim_params, dada2_preprocess_script, is_paired_end)
	dada2_analysis(dada2_preprocess.out.filtered_reads, dada2_preprocess.out.trim_table, dada2_analysis_script, is_paired_end)
}
