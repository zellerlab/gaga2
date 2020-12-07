#!/usr/bin/env nextflow

params.figaro = "/g/scb/bork/schudoma/16S/figaro/figaro.py"
params.scripts = "/g/scb/bork/schudoma/16S/scripts"
params.amplicon_length = 465
params.forward_primer = 21
params.reverse_primer = 17
params.min_overlap = 20

params.output_dir = "output"

Channel
	.fromFilePairs(params.input_dir + "/*/*_*[12].{fastq,fq,fastq.gz,fq.gz}")
	.set { samples_ch }

samples_ch.into { run_figaro_input_ch; run_trimming_input_ch; }


process run_figaro {
	conda "anaconda::numpy>=1.13.1 anaconda::scipy>=1.2.1 anaconda::matplotlib>=3.0.2"
	publishDir "${params.output_dir}/figaro/${sampleId}"

	input:
	set sampleId, file(reads) from run_figaro_input_ch

	output:
	stdout run_figaro_stdout
	file("figaro_out/forwardExpectedError.png")
	file("figaro_out/reverseExpectedError.png")
	tuple val(sampleId), file("figaro_out/trimParameters.json") into run_figaro_output_ch

	script:
	"""
	mkdir figaro_in
	cd figaro_in
	ln -s ../${reads[0]}
	ln -s ../${reads[1]}
	cd ..
	python ${params.figaro} -i figaro_in -o figaro_out -a ${params.amplicon_length} -f ${params.forward_primer} -r ${params.reverse_primer} -m ${params.min_overlap}
	"""

}

//run_figaro_stdout.view { it }

run_trimming_input_ch = run_trimming_input_ch.join(run_figaro_output_ch)

process run_trimming {
	// dada2's dependencies are not right
	//conda "python>=3.7 conda-forge::gcc r-essentials r-base bioconda::bioconductor-dada2 conda-forge::r-rcpp"
	conda "r-essentials r-base bioconda::bioconductor-dada2 conda-forge::r-rcpp"
	input:
	set sampleId, file(reads), file(trim_params) from run_trimming_input_ch

	output:
	stdout run_trimming_stdout

	script:
	"""
	# python ${params.scripts}/trim_params.py $trim_params ${reads[0]}
	Rscript ${params.scripts}/test.R
	"""

}

run_trimming_stdout.view { it }


