#!/usr/bin/env Rscript
.libPaths(c('/g/scb2/zeller/SHARED/software/R/3.5', .libPaths()))
library("dada2");packageVersion("dada2")

args = commandArgs(trailingOnly=TRUE)
input_dir = args[1]
output_dir = args[2]

length_cut = as.numeric(c(args[3], args[4]))
print(length_cut)
if (length_cut[1] == -1 || length_cut[2] == -1) {
	length_cut = 0
}
exp_err = as.numeric(c(2,2)) # get from cargs? apparently not.

list.files(input_dir)
sample_ids = basename(list.files(input_dir))
print(sample_ids)

r1_raw = sort(list.files(file.path(input_dir, sample_ids), pattern="_R?1.(fastq|fq)(.gz)?", full.names=TRUE))
r2_raw = sort(list.files(file.path(input_dir, sample_ids), pattern="_R?2.(fastq|fq)(.gz)?", full.names=TRUE))

print(r1_raw)
print(r2_raw)

#TODO: this only generates one plot 
pdf("read_quality.pdf", useDingbats=FALSE, width=12, height=6)
for (i in 1:length(r1_raw)) {
	print(plotQualityProfile(c(r1_raw[i], r2_raw[i])))
}
dev.off()

r1_filtered = file.path(output_dir, "filtered", basename(r1_raw))
r2_filtered = file.path(output_dir, "filtered", basename(r2_raw))

print(r1_filtered)
print(r2_filtered)


out = filterAndTrim(r1_raw, r1_filtered, r2_raw, r2_filtered, truncLen=length_cut,
                    maxN=0, maxEE=exp_err, truncQ=2, rm.phix=TRUE,
                    compress=TRUE, multithread=TRUE)
write.table(out, file="filter_trim_table.tsv", sep="\t")

# update read files
keep = file.exists(r1_filtered) & file.exists(r2_filtered) & out[,"reads.out"] > 100
r1_filtered = r1_filtered[keep]
r2_filtered = r2_filtered[keep]
sample_ids = sample_ids[keep]

r1_error = learnErrors(r1_filtered, multithread=TRUE)
r2_error = learnErrors(r2_filtered, multithread=TRUE)

pdf("error_model.pdf", width=12, height=12, paper='special')
print(plotErrors(r1_error, nominalQ=TRUE))
print(plotErrors(r2_error, nominalQ=TRUE))
dev.off()


r1_uniq = derepFastq(r1_filtered, verbose=TRUE)
r2_uniq = derepFastq(r2_filtered, verbose=TRUE)
names(r1_uniq) = sample_ids
names(r2_uniq) = sample_ids

r1_dada = dada(r1_uniq, err=r1_error, multithread=TRUE)
r2_dada = dada(r2_uniq, err=r2_error, multithread=TRUE)

#Inspecting the returned dada-class object:
r1_dada[[1]]


merged = mergePairs(r1_dada, r1_uniq, r2_dada, r2_uniq, verbose=TRUE)

# Inspect the merger data.frame from the first sample
head(merged[[1]])

seqtab = makeSequenceTable(merged)
n_OTU = dim(seqtab)[2]
# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

seqtab.nochim = removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
n_OTU_after_removing_chimeras = dim(seqtab.nochim)[2]
# relative abundance of the chimeric sequences
rel_ab_chimeras = 1 - (sum(seqtab.nochim)/sum(seqtab))
table(nchar(getSequences(seqtab.nochim)))

getN = function(x) sum(getUniques(x))
track = cbind(out[keep,], sapply(r1_dada, getN), sapply(r2_dada, getN), sapply(merged, getN), rowSums(seqtab.nochim))
colnames(track) = c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) = sample_ids
head(track)
write.table(track, file = "summary_table.tsv", sep="\t")

save(seqtab, seqtab.nochim, r1_error, r2_error, track, file = "result.RData")
