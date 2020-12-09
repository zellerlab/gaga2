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

samples_ch.into { run_figaro_input_ch; run_dada2_input_ch; }


process run_figaro_all {
	//conda "anaconda::numpy anaconda::scipy anaconda::matplotlib"
	conda "anaconda::numpy>=1.13.1 anaconda::scipy>=1.2.1 anaconda::matplotlib>=3.0.2"
	//conda "etc/figaro.yml"
	publishDir "${params.output_dir}/figaro"

	input:
	run_figaro_input_ch.collect()

	output:
	stdout run_figaro_all_stdout
	file("figaro_out/forwardExpectedError.png")
	file("figaro_out/reverseExpectedError.png")
	file("figaro_out/trimParameters.json") into run_figaro_all_output_ch

	shell:
	"""
	python ${params.scripts}/gather_fastq_files.py ${params.input_dir} figaro_in
	python ${params.figaro} -i figaro_in -o figaro_out -a ${params.amplicon_length} -f ${params.forward_primer} -r ${params.reverse_primer} -m ${params.min_overlap}
	"""
}


process run_dada2 {
	// dada2's dependencies are not right
	//conda "python>=3.7 conda-forge::gcc r-essentials r-base bioconda::bioconductor-dada2 conda-forge::r-rcpp"
	//conda "r-essentials r-base bioconda::bioconductor-dada2 conda-forge::r-rcpp"
	publishDir "${params.output_dir}/qc"
	input:
	file(trim_params) from run_figaro_all_output_ch
	run_dada2_input_ch.collect()

	output:
	stdout run_dada2_stdout
	file("read_quality.pdf")

	script:
	"""
	ls -l
	tparams=\$(python ${params.scripts}/trim_params.py $trim_params)
	echo \$tparams
	module load R/3.5.0-foss-2017b-X11-20171023
	Rscript --vanilla ${params.scripts}/dada2_qc.R ${params.input_dir} ${params.output_dir} \$tparams
	"""
}

run_dada2_stdout.view { it }


