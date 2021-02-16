#!/usr/bin/env Rscript
library("dada2");packageVersion("dada2")
library("tidyverse")
library("cowplot")

print("XXX")

# handle command line args
args = commandArgs(trailingOnly=TRUE)
input_dir = args[1]
output_dir = args[2]

length_cut = as.numeric(c(args[3], args[4]))
print(length_cut)
if (length_cut[1] == -1 || length_cut[2] == -1) {
	length_cut = 0
}
exp_err = as.numeric(c(2,2)) # get from cargs? apparently not.

if (length(args) >= 5) {
	nthreads = as.numeric(c(args[5]))
} else {
	nthreads = TRUE
}

# get the read files and sample names
list.files(input_dir)
sample_ids = basename(list.files(input_dir))
print(sample_ids)

r1_raw = sort(list.files(file.path(input_dir, sample_ids), pattern = "_R?1.(fastq|fq)(.gz)?", full.names = TRUE))
r2_raw = sort(list.files(file.path(input_dir, sample_ids), pattern = "_R?2.(fastq|fq)(.gz)?", full.names = TRUE))
print(r1_raw)
print(r2_raw)

# get the quality profiles
x.f = plotQualityProfile(r1_raw, aggregate = TRUE)
x.r = plotQualityProfile(r2_raw, aggregate = TRUE)
g = plot_grid(x.f, x.r, labels = c('forward', 'reverse'))
ggsave(g, filename = "read_quality.pdf",
       width = 12, height = 7, useDingbats=FALSE)

# perform filtering and trimming
r1_filtered = file.path(output_dir, "filtered", basename(r1_raw))
r2_filtered = file.path(output_dir, "filtered", basename(r2_raw))

print(r1_filtered)
print(r2_filtered)

out = filterAndTrim(r1_raw, r1_filtered, r2_raw, r2_filtered, truncLen=length_cut,
                    maxN=0, maxEE=exp_err, truncQ=2, rm.phix=TRUE,
                    compress=TRUE, multithread=nthreads)
# colnames(out) = c("sample", "reads.in", "reads.out")
rownames(out) = str_replace(rownames(out), "_R?[12].(fastq|fq)(.gz)?", "")
write.table(out, file="filter_trim_table.tsv", sep="\t")

# remove 
MIN_READ_THRESHOLD = 100
keep = file.exists(r1_filtered) & file.exists(r2_filtered) & out[,"reads.out"] >= MIN_READ_THRESHOLD
r1_remove = r1_filtered[!keep] #file.exists(r1_filtered) & out["reads.out"] < MIN_READ_THRESHOLD]
r2_remove = r2_filtered[!keep] #file.exists(r2_filtered) & out["reads.out"] < MIN_READ_THRESHOLD]
sapply(r1_remove, file.remove)
sapply(r2_remove, file.remove)

out = out[keep,]
write.table(out, file="filter_trim_table.final.tsv", sep="\t")

# get the quality profiles post-qc
x.f = plotQualityProfile(r1_filtered, aggregate = TRUE)
x.r = plotQualityProfile(r2_filtered, aggregate = TRUE)
g = plot_grid(x.f, x.r, labels = c('forward', 'reverse'))
ggsave(g, filename = "read_quality_postqc.pdf",
       width = 12, height = 7, useDingbats=FALSE)

# update read files
#keep = file.exists(r1_filtered) & file.exists(r2_filtered) & out[,"reads.out"] > 100
#r1_filtered = r1_filtered[keep]
#r2_filtered = r2_filtered[keep]
#sample_ids = sample_ids[keep]
