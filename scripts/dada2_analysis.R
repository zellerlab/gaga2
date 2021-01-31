#!/usr/bin/env Rscript
library("dada2");packageVersion("dada2")
library("tidyverse")
library("cowplot")

print("YYY")

# handle command line args
args = commandArgs(trailingOnly=TRUE)
input_dir = args[1]
output_dir = args[2]

filter_table = read.table(args[3])

if (length(args) >= 4) {
	nthreads = as.numeric(c(args[4]))
} else {
	nthreads = TRUE
}

# get the read files and sample names
list.files(input_dir)
sample_ids = basename(list.files(input_dir))
print(sample_ids)

r1_filtered = sort(list.files(file.path(input_dir, sample_ids), pattern = "_R?1.(fastq|fq)(.gz)?", full.names = TRUE))
r2_filtered = sort(list.files(file.path(input_dir, sample_ids), pattern = "_R?2.(fastq|fq)(.gz)?", full.names = TRUE))
r1_filtered = sort(list.files(file.path(input_dir), pattern = "_R?1.(fastq|fq)(.gz)?", full.names = TRUE))
r2_filtered = sort(list.files(file.path(input_dir), pattern = "_R?2.(fastq|fq)(.gz)?", full.names = TRUE))

print(r1_filtered)
print(r2_filtered)

# learn error rate
r1_error = learnErrors(r1_filtered, multithread=nthreads)
r2_error = learnErrors(r2_filtered, multithread=nthreads)

pdf("error_model.pdf", width=12, height=12, paper='special')
print(plotErrors(r1_error, nominalQ=TRUE))
print(plotErrors(r2_error, nominalQ=TRUE))
dev.off()

# remove replicates
r1_uniq = derepFastq(r1_filtered, verbose=TRUE)
r2_uniq = derepFastq(r2_filtered, verbose=TRUE)
print(r1_uniq[[1]])
print(c("LENGTHCHECK:", length(r1_uniq), length(r2_uniq), length(sample_ids)))
#names(r1_uniq) = sample_ids
#names(r2_uniq) = sample_ids
r1_dada = dada(r1_uniq, err=r1_error, multithread=nthreads)
r2_dada = dada(r2_uniq, err=r2_error, multithread=nthreads)

#Inspecting the returned dada-class object:
r1_dada[[1]]

# merge read pairs
print("merging")
merged = mergePairs(r1_dada, r1_uniq, r2_dada, r2_uniq, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(merged[[1]])

# construct sequence table
print("making sequence table")
seqtab = makeSequenceTable(merged)
n_OTU = dim(seqtab)[2]
# Inspect distribution of sequence lengths
print(dim(seqtab))
print("ASV length")
print(table(nchar(getSequences(seqtab))))

# remove chimeras
seqtab.nochim = removeBimeraDenovo(seqtab, method="consensus", multithread=nthreads, verbose=TRUE)
n_OTU_after_removing_chimeras = dim(seqtab.nochim)[2]
print(dim(seqtab.nochim))
asv.table = t(seqtab.nochim)

# export asv table
asvs = tibble(id=paste0('ASV_', seq_len(nrow(asv.table))),
               ASV=rownames(asv.table))
write_tsv(asvs, 'ASVs.tsv')
rownames(asv.table) = asvs$id
write.table(asv.table, file = 'asv_table.tsv',
            sep='\t', quote = FALSE, row.names = TRUE, col.names = TRUE)

# relative abundance of the chimeric sequences
rel_ab_chimeras = 1 - (sum(seqtab.nochim)/sum(seqtab))
table(nchar(getSequences(seqtab.nochim)))

# track reads
getN = function(x) sum(getUniques(x))
track = cbind(filter_table, sapply(r1_dada, getN), sapply(r2_dada, getN), sapply(merged, getN), rowSums(seqtab.nochim))
print(c("TRACK", length(track), length(sample_ids)))
colnames(track) = c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
# rownames(track) = sample_ids
print(head(track))
print("Filtering")
print(summary(track[,6] / track[,1]))
write.table(track, file = "summary_table.tsv", sep="\t")

save(seqtab, seqtab.nochim, r1_error, r2_error, track, file = "result.RData")

# prevalance sanity check
pdf('dada2_figures.pdf')
temp = prop.table(asv.table, 2)
print(temp)
hist(log10(temp), 100, main='abundance')
prev = rowMeans(asv.table!=0)
hist(prev, 50, main='prevalence')
mean.ab = rowMeans(log10(temp + 1e-05))
hist(mean.ab, 50, main='log.abundance')
print("Prevalence")
print(summary(prev))
print("Mean abundance")
print(summary(mean.ab))
plot(prev, mean.ab, main='prevalence vs abundance')
plot(nchar(rownames(asv.table)), prev, main='ASV length vs prevalence')
dev.off()
