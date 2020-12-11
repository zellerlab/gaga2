#!/usr/bin/env nextflow

/*
params.figaro = "/home/schudoma/gaga2_test/gaga2/figaro/figaro.py"
params.scripts = "/home/schudoma/gaga2_test/gaga2/scripts"
params.envs = "/home/schudoma/gaga2_test/gaga2/etc"
params.output_dir = "/home/schudoma/gaga2_test/output"
*/
/*
Test params for Bullman 2017
params.amplicon_length = 465
params.forward_primer = 21
params.reverse_primer = 17
params.min_overlap = 20
*/


def helpMessage() {
	log.info """
	
	This is the gaga2 (figaro/dada2) 16S amplicon sequencing processing pipeline of the Zeller lab @ EMBL.

	Usage:

	The typical command for running the pipeline is as follows:

        nextflow run gaga2.nf -c <config_file> --input_dir </path/to/input_dir> --output_dir </path/to/output_dir> ...

    Mandatory arguments:

            --input_dir               path to input read files (needs to be absolute, 1 subdirectory per sample)
                                      reads need to be Illumina paired end amplicon seq (e.g. MiSeq)
            --output_dir              path to output directory (absolute)
            --amplicon_length         expected amplicon length
            --left_primer             length of left primer
            --right_primer            length of right primer

    Optional arguments:

            --min_overlap             Minimum read pair overlap [bp] (default=20)
            --help                    Show this help.

	""".stripIndent()
}

if (params.help) {
	helpMessage()
	exit 0
}

if ( !params.min_overlap ) {
	params.min_overlap = 20
}


Channel
	.fromFilePairs(params.input_dir + "/*/*_*[12].{fastq,fq,fastq.gz,fq.gz}")
	.set { samples_ch }

samples_ch.into { run_figaro_input_ch; run_dada2_input_ch; }


process run_figaro_all {
	//conda "anaconda::numpy anaconda::scipy anaconda::matplotlib"
	//	conda "anaconda::numpy>=1.13.1 anaconda::scipy>=1.2.1 anaconda::matplotlib>=3.0.2"
	conda "${params.envs}/figaro.yml"
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
	which python
	python ${params.scripts}/check_readsets.py ${params.input_dir} ${params.output_dir}
	python ${params.scripts}/gather_fastq_files.py ${params.input_dir} figaro_in
	if [[ ! -f ${params.output_dir}/SKIP_FIGARO ]]; then
	python ${params.figaro} -i figaro_in -o figaro_out -a ${params.amplicon_length} -f ${params.forward_primer} -r ${params.reverse_primer} -m ${params.min_overlap}
	fi

	mkdir -p figaro_out
	touch figaro_out/forwardExpectedError.png
	touch figaro_out/reverseExpectedError.png
	touch figaro_out/trimParameters.json
	"""
}


process run_dada2 {
	// dada2's dependencies are not right
	//conda "python>=3.7 conda-forge::gcc r-essentials r-base bioconda::bioconductor-dada2 conda-forge::r-rcpp"
	//conda "r-essentials r-base bioconda::bioconductor-dada2 conda-forge::r-rcpp"
	publishDir "${params.output_dir}/dada2"
	input:
	file(trim_params) from run_figaro_all_output_ch
	run_dada2_input_ch.collect()

	output:
	stdout run_dada2_stdout
	file("read_quality.pdf")
	file("filter_trim_table.tsv")
	file("error_model.pdf")
	file("summary_table.tsv")
	file("result.RData")
	script:
	"""
	ls -l
	tparams=\$(python ${params.scripts}/trim_params.py $trim_params)
	echo \$tparams
	module load R/3.5.0-foss-2017b-X11-20171023
	Rscript --vanilla ${params.scripts}/dada2.R ${params.input_dir} ${params.output_dir} \$tparams
	"""
}

run_dada2_stdout.view { it }


