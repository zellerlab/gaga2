#!/usr/bin/env Rscript
.libPaths(c('/g/scb2/zeller/SHARED/software/R/3.5', .libPaths()))
library("dada2");packageVersion("dada2")

args = commandArgs(trailingOnly=TRUE)
input_dir = args[1]
#output_dir = args[2]

list.files(input_dir)
sample_ids = basename(list.files(input_dir))
print(sample_ids)

r1_files = sort(list.files(file.path(input_dir, sample_ids), pattern="_R?1.(fastq|fq)(.gz)?", full.names=TRUE))
r2_files = sort(list.files(file.path(input_dir, sample_ids), pattern="_R?2.(fastq|fq)(.gz)?", full.names=TRUE))

print(r1_files)
print(r2_files)

pdf("read_quality.pdf", useDingbats=FALSE, width=12, height=6)
for (i in 1:length(r1_files)) {
	print(plotQualityProfile(c(r1_files[i], r2_files[i])))
}
dev.off()

w_dir = "/g/scb2/zeller/milanese/other_projects/helena_2019_3/"
path_data = "/g/scb2/zeller/milanese/DATA/private_DATASET/helena_2019_3/"
# forward and reverse, where to truncate
length_cut = c(180,180)

# YOU NEED TO CREATE A DIRECTORY IN path_data WITH NAME filtered

# check how the samples names are separated

# -----------------------------------------------------------------
# Load DADA 2
# -----------------------------------------------------------------
cat("1. load dada2\n-----------------\n")
.libPaths(c('/g/scb2/zeller/SHARED/software/R/3.5', .libPaths()))
library("dada2");packageVersion("dada2")

# -----------------------------------------------------------------
# Find fastq files
# -----------------------------------------------------------------
cat("2. find files\n-----------------\n")
list.files(path_data)
# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path_data, pattern="1_sequence.txt", full.names = TRUE))
fnRs <- sort(list.files(path_data, pattern="2_sequence.txt", full.names = TRUE))

# Extract sample names, assuming filenames have format: SAMPLENAME.XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_1_sequence"), `[`, 1)

cat("\nnumber of samples:\n")
cat(length(sample.names))
cat("\n")



# -----------------------------------------------------------------
# Filter and trim reads
# -----------------------------------------------------------------
cat("3. filter\n-------------\n")
# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path_data, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path_data, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=length_cut,
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
					                      compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE
										  head(out)

										  # remove the one that did not pass the filter &
										  # remove also the one with few reads
										  not.lost <- file.exists(filtFs) & out[,"reads.out"] > 100
										  filtFs <- filtFs[not.lost]
										  filtRs <- filtRs[not.lost]
										  sample.names <- sample.names[not.lost]

										  # remove also the one with few reads
										  #keep <- out[,"reads.out"] > 100
										  #filtFs <- file.path(filtpathF, fastqFs)[keep]
										  #filtRs <- file.path(filtpathR, fastqRs)[keep]

										  # -----------------------------------------------------------------
										  # Learn error rates
										  # -----------------------------------------------------------------
										  cat("4. learn error\n------------------\n")
										  errF <- learnErrors(filtFs, multithread=TRUE)
										  errR <- learnErrors(filtRs, multithread=TRUE)

										  pdf(paste0(w_dir,"error_model.pdf"),width=12,height=12,paper='special')
										  print(plotErrors(errF, nominalQ=TRUE))
										  print(plotErrors(errR, nominalQ=TRUE))
										  dev.off()

										  # -----------------------------------------------------------------
										  # Dereplication
										  # -----------------------------------------------------------------
										  cat("5. dereplication\n---------------------\n")
										  derepFs <- derepFastq(filtFs, verbose=TRUE)
										  derepRs <- derepFastq(filtRs, verbose=TRUE)
										  # Name the derep-class objects by the sample names
										  names(derepFs) <- sample.names
										  names(derepRs) <- sample.names

										  # -----------------------------------------------------------------
										  # Core of dada2 
										  # -----------------------------------------------------------------
										  cat("6. dada core\n----------------\n")
										  dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
										  dadaRs <- dada(derepRs, err=errR, multithread=TRUE)

										  #Inspecting the returned dada-class object:
										  dadaFs[[1]]

										  # -----------------------------------------------------------------
										  # Merge paired reads
										  # -----------------------------------------------------------------
										  cat("7. merge pairs\n----------------\n")
										  mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
										  # Inspect the merger data.frame from the first sample
										  head(mergers[[1]])

										  # -----------------------------------------------------------------
										  # Construct sequence table
										  # -----------------------------------------------------------------
										  cat("8. construct the table\n-----------------------\n")
										  seqtab <- makeSequenceTable(mergers)

										  n_OTU = dim(seqtab)[2]

										  # Inspect distribution of sequence lengths
										  table(nchar(getSequences(seqtab)))

										  # -----------------------------------------------------------------
										  # Remove chimeras
										  # -----------------------------------------------------------------
										  cat("9. remove chimeras\n---------------------\n")
										  seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
										  n_OTU_after_removing_chimeras = dim(seqtab.nochim)[2]

										  # relative abundance of the chimeric sequences
										  rel_ab_chimeras = 1 - (sum(seqtab.nochim)/sum(seqtab))

										  table(nchar(getSequences(seqtab.nochim)))

										  # -----------------------------------------------------------------
										  # Track reads through the pipeline
										  # -----------------------------------------------------------------
										  cat("10. track reads\n-------------------------------\n")
										  getN <- function(x) sum(getUniques(x))
										  track <- cbind(out[not.lost,], sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
										  colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
										  rownames(track) <- sample.names
										  head(track)

										  # -----------------------------------------------------------------
										  # Save file
										  # -----------------------------------------------------------------
										  cat("11. save file\n------------------------------------------------------\n")
										  save(seqtab,seqtab.nochim, errF, errR,track , file = paste0(w_dir,"result.RData"))

										  #load(file = paste0(w_dir,"result.RData"))




