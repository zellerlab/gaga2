#!/usr/bin/env nextflow

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
            --nthreads                Number of threads used for dada2 (default=8)
            -work-dir, -w             Path to working directory
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

//TODO: get the default from config? -> but needs to be provided to dada2 call (is there a way to get "threads" from params?)
if ( !params.nthreads ) {
	params.nthreads = 8
}


Channel
	.fromFilePairs(params.input_dir + "/*/*_*[12].{fastq,fq,fastq.gz,fq.gz}")
	.set { samples_ch }

samples_ch.into { run_figaro_input_ch; run_dada2_input_ch; }


process run_figaro_all {
	//conda "${params.envs}/figaro.yml"
	conda "${params.envs}/figaro_env"
	publishDir "${params.output_dir}/figaro", mode: "link"

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
	figaro -i figaro_in -o figaro_out -a ${params.amplicon_length} -f ${params.forward_primer} -r ${params.reverse_primer} -m ${params.min_overlap}
	fi

	mkdir -p figaro_out
	touch figaro_out/forwardExpectedError.png
	touch figaro_out/reverseExpectedError.png
	touch figaro_out/trimParameters.json
	"""
}

process dada2_preprocess {
	publishDir "${params.output_dir}/dada2", mode: "link"

	input:
	file(trim_params) from run_figaro_all_output_ch
	run_dada2_input_ch.collect()

	output:
	stdout dada2_preprocess_stdout
	file("read_quality.pdf")
    file("read_quality_postqc.pdf")
    file("filter_trim_table.tsv")
    file("filter_trim_table.final.tsv") into dada2_preprocess_output_ch
	file("dada2_preprocess.log")

	script:
	"""
	tparams=\$(python ${params.scripts}/trim_params.py $trim_params)
	echo \$tparams
	dada2_preprocess.R ${params.input_dir} ${params.output_dir} \$tparams ${params.nthreads} > dada2_preprocess.log
	"""
	//Rscript --vanilla ${params.scripts}/dada2_preprocess.R ${params.input_dir} ${params.output_dir} \$tparams ${params.nthreads} > dada2_preprocess.log
}

process dada2_analysis {
	publishDir "${params.output_dir}/dada2", mode: "link"

	input:
	file(filter_trim_table) from dada2_preprocess_output_ch

	output:
	file("dada2_analysis.log")
	file("error_model.pdf")
	file("summary_table.tsv")
	file("result.RData")
	file("dada2_figures.pdf")
	file("ASVs.tsv")
	file("asv_table.tsv")

	script:
	"""
	dada2_analysis.R ${params.output_dir}/filtered ${params.output_dir} ${filter_trim_table} ${params.nthreads} > dada2_analysis.log
	"""
	//Rscript --vanilla ${params.scripts}/dada2_analysis.R ${params.output_dir}/filtered ${params.output_dir} ${filter_trim_table} ${params.nthreads} > dada2_analysis.log

}
