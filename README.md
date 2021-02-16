# gaga2 - automated 16S amplicon analysis with Figaro/DADA2

## Installation instructions
`TBD`

## Usage instructions

### Input
`gaga2` takes as input Illumina paired-end 16S amplicon sequences (e.g. sequenced on a MiSeq).

#### Prerequisites
* Read files need to be named according the typical pattern `<prefix=sample_id>_*[12].{fastq,fq,fastq.gz,fq.gz}`
and need to be arranged in a sample-based directory structure:

```
<project_dir> (aka "input_dir")
     |___ <sample_1>
     |        |____ <sample_1_forward_reads>
     |        |____ <sample_2_reverse_reads>
     |
     |___ <sample_2>
     |        |____ <empty samples will be ignored>
     |        
     |___ <sample_n>
              |____ <sample_n_forward_reads>
              |____ <sample_n_reverse_reads>
```

* If input reads have been "invasively" preprocessed (as evidenced by length differences between reads), 
`gaga2` will generate an empty sentinel file `<output_dir>/SKIP_FIGARO`.
This will prevent `Figaro` from executing and should automatically advance to `dada2` processing.

### Running gaga2
The typical command for running `gaga2` is

`gaga2.nf --input_dir <input_dir> --output_dir <output_dir> --amplicon_length <int> --left_primer <int:length_left_primer> --right_primer <int:length_right_primer>`

#### Mandatory arguments
* `input_dir` is the project directory mentioned above. **Recommended: absolute path**
* `output_dir` will be created automatically. **Recommended: absolute path**
* `amplicon_length`, `left_primer`, and `right_primer` are derived from the experiment parameters 

#### Optional arguments
* `--min_overlap` of read pairs is 20bp by default
* `-w, -work-dir` should be set to `<output_dir>/work`, otherwise it will be automatically placed in the current directory.


### internal beta-testing instructions
* `source /g/scb2/zeller/schudoma/software/wrappers/gaga2_wrapper` **before** submitting job to cluster
* please report issues/requests/feedback in the github issue tracker 

