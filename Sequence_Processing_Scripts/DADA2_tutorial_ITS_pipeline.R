# This is pulled from Dr. Benjamin Callehan's DADA2 tutorial
## here is the link: https://benjjneb.github.io/dada2/tutorial.html
## here is a link to the publication describing the DADA2 algorithm: https://www.nature.com/articles/nmeth.3869#Sec2


# **** If running R on the cluster, do the following before you run DADA2:
## srun --partition=aronsonlab --mem=400gb --cpus-per-task 4 --ntasks 1 --time 3-00:00:00 --pty bash -l

## module load tmux
## tmux
## R or nvim DADA2_tutorial_pipeline.R
## then load libraries and check DADA2 package

getwd()
#install.packages("devtools")
#library("devtools")
#devtools::install_github("benjjneb/dada2", ref="taxonomy-segfault-fix")
packageVersion("dada2") # Should be 1.26
suppressPackageStartupMessages({
  library(dada2)
  library(tidyr)
  library(ggpubr)
  library(decontam)
})

### Set the Path ####
path <- "/bigdata/aronsonlab/shared/dstev/seqsLT/Trimmed_Seqs" # CHANGE ME to the directory containing the fastq files after unzipping.
# path <- "/bigdata/aronsonlab/shared/XCZO/ITS1/XCZO_ITS1_HFreund/sequences" # cluster path
# path <- "/bigdata/aronsonlab/shared/HannahFreund/MtStHelens_ITS_MM/Trimmed_Seqs/its2_R1_reads"
list.files(path)

#### Read in sample/file names ####

# Now we read in the names of the fastq files, and perform some string manipulation to get matched lists of the forward and reverse fastq files.
# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq

fnFs <- sort(list.files(path, pattern="_R1_clean.fastq", full.names = TRUE))
#fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_clean.fastq", full.names = TRUE)); save.image(file = "mydada_its2.Rdata")

fnFs

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_R1"), `[`, 1); save.image(file = "mydada_its2.Rdata") #pattern where you want to split the name
sample.names

#### Check Read Quality Profiles ####

plot1<-plotQualityProfile(fnFs[1:2])
plot2<-plotQualityProfile(fnRs[1:2])
ggsave(plot1,filename = "its2_pretrim_DADA2_F_quality.pdf", width=15, height=12, dpi=600)
ggsave(plot2,filename = "its2_pretrim_DADA2_R_quality.pdf", width=15, height=12, dpi=600)

#### Filter + Trim ####

## assign filenames for filtered fastq.gz files
# Place filtered files in filtered/ subdirectory
path
filtFs <- file.path(path, "Filtered", paste0(sample.names, "_F_filtered.fastq.gz"))
filtRs <- file.path(path, "Filtered", paste0(sample.names, "_R_filtered.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names; save.image(file = "mydada_its2.Rdata")

## Standard filtering parameters: maxN=0 (DADA2 requires no Ns), truncQ=2, rm.phix=TRUE and maxEE=2.
## The maxEE parameter sets the maximum number of “expected errors” allowed in a read, which is a better filter than simply averaging quality scores.
# * if you are only doing F reads, remove the "truncLen" command - truncLen=c(240,160) [for PE reads]
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, trimLeft=0,
                     maxN=0, maxEE = c(5, 10), truncQ=2, minLen=50, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE); save.image(file = "mydada_its2.Rdata") # On Windows set multithread=FALSE

## If too few reads are passing the filter, consider relaxing maxEE, perhaps especially on the reverse reads (eg. maxEE=c(2,5))
# and reducing the truncLen to remove low quality tails. Remember though, when choosing truncLen for paired-end reads you must maintain overlap after truncation in order to merge them later
head(out)

# removes files that were not included in output because 0 reads passed filter step
filtFs <- filtFs[file.exists(filtFs)]
filtRs <- filtRs[file.exists(filtRs)] ; save.image(file = "mydada_its2.Rdata")

#### Learn the error rates ####

errF <- learnErrors(filtFs, multithread=TRUE); save.image(file = "mydada_its2.Rdata")
errR <- learnErrors(filtRs, multithread=TRUE); save.image(file = "mydada_its2.Rdata")
# The learnErrors method learns this error model from the data by alternating estimation of the error rates and inference of sample composition until they converge on a jointly consistent solution.
# As in many machine-learning (ML) problems, the algorithm must begin with an initial guess, for which the maximum possible error rates in this data are used (the error rates if only the most abundant sequence is correct and all the rest are errors)

plot_error<-plotErrors(errF, nominalQ=TRUE)## sanity check by visualizing estimated error rates -- should see error rates drop w/ increased quality
ggsave(plot_error,filename = "its2_errormodel_DADA2.pdf", width=15, height=15, dpi=600)

## Error rates for each possible transition: A2A (A -> A), A2C (A -> C), etc
## Points are observed error rates for each consensus quality score
## Black line shows estimated error rates after convergence of the ML algorithm
## Red line shows error rates expected under the nominal definition of the Q-score


#### Sample Inference! ####
# DADA2's selfConsist mode alternates sample inference (conditional on the parameters of the error model) with parameter estimation (conditional on the inferred sample composition) until convergence, at which point jointly consistent estimates of the error parameters and sample composition are reported.
# https://www.nature.com/articles/nmeth.3869#methods
# this step (w/ dada command) won't work if you've lost samples due to not surpassing the filtering step
filtFs <- filtFs[file.exists(filtFs)] # removes files that were not included in output because 0 reads passed filter step
filtRs <- filtRs[file.exists(filtRs)]
dadaFs <- dada(filtFs, err=errF, multithread=TRUE, pool=TRUE); save.image(file = "mydada_its2.Rdata") # pseudo pooling is computationally more efficient but similar in results to pooling; pool = True will pool samples together before sample inference
dadaRs <- dada(filtRs, err=errR, multithread=TRUE, pool=TRUE); save.image(file = "mydada_its2.Rdata")

dadaFs[1] # Returns first section of dada-class object {one sample}
dadaRs[1]
#dadaRs[[1]]

#### Merge Paired End Reads (obtain full denoised seqs) ####

## Merging is performed by aligning the denoised forward reads with the reverse-complement of the corresponding denoised reverse reads, and then constructing the merged “contig” sequences.
## By default, merged sequences are only output if the forward and reverse reads overlap by at least 12 bases, and are identical to each other in the overlap region (but these conditions can be changed via function arguments)
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE); save.image(file = "mydada_its2.Rdata")

# Inspect the merger data.frame from the first sample
head(mergers[[1]])

# NOTE: If sequences do not merge, or if reverse reads are too low quality, use only Forward reads for rest of workflow.
# replace mergers object with dadaFs object to create seqtab below

#### Construct Sequence Table (ASVs)! ####

## sequence table is a matrix of rows (samples) and columns (ASVs)
## NOTE: reads were not able to merge, going to use only F reads for ITS2#
seqtab <- makeSequenceTable(dadaFs); save.image(file = "mydada_its2.Rdata")
dim(seqtab)

#### Inspect distribution of sequence lengths ####
table(nchar(getSequences(seqtab)))

# Multiple length reads - ITS1 are of variable length (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5309391/)

## Look at merged sequences in a plot -- see their distribution and frequency of sequences of certain length
compare_reads_plot = ggdensity(rowSums(seqtab), fill = "blue4", alpha = 0.7); save.image(file = "mydada_its2.Rdata")
pdf(file = "its2_compare_plots.pdf", height = 10, width = 12) # Make empty pdf
compare_reads_plot # Add plot to empty pdf
dev.off()


#### Remove Chimeras ####
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE); save.image(file = "mydada_its2.Rdata")
# Chimeric sequences are identified if they can be exactly reconstructed by combining a left-segment and a right-segment from two more abundant “parent” sequences
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab) # comparing reads after chimera removal over total reads (after filtering)
# Most of your reads should remain after chimera removal (it is not uncommon for a majority of sequence variants to be removed though)

#### Track reads through the pipeline: see how many reads made it through each step of pipeline ####
getN <- function(x) sum(getUniques(x))
#track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim)); save.image(file = "mydada_its2.Rdata")
track <- cbind(out, sapply(dadaFs, getN),sapply(dadaRs, getN),rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs) ;  sapply(dadaRs, getN), sapply(mergers, getN),
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "nonchim") # remove whichever labels you didn't include
rownames(track) <- sample.names
head(track)
write.table(track,"ITS2_tracking_reads_dada2.txt",sep="\t",row.names=TRUE,col.names=TRUE); save.image(file = "mydada_its2.Rdata")

# If a majority of reads failed to merge, you may need to revisit the truncLen parameter used in the filtering step and make sure that the truncated reads span your amplicon.
# If a majority of reads were removed as chimeric, you may need to revisit the removal of primers, as the ambiguous nucleotides in unremoved primers interfere with chimera identification

#### Assign Taxonomy! May have to untar/zip file to run this step ####
## download trainig zip file: https://zenodo.org/record/3986799

# DADA2 package provides a native implementation of the naive Bayesian classifier method for this purpose.
# The assignTaxonomy function takes as input a set of sequences to be classified and a training set of reference sequences with known taxonomy, and outputs taxonomic assignments with at least minBoot bootstrap confidence.
# For fungal taxonomy, the General Fasta release files from the UNITE ITS database can be downloaded and used as the reference.
# NOTE: include tryRC=TRUE for ITS2 regions because primers can be in reverse, and thus forward ITS2 reads do not align to UNITE db
taxa <- assignTaxonomy(seqtab.nochim, "/bigdata/aronsonlab/shared/RefDBs/UNITE/sh_general_release_dynamic_all_29.11.2022.fasta", multithread=TRUE,tryRC=TRUE)

#### Looking at the taxonomic assingments ####

taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
taxa.print<-as.data.frame(apply(taxa.print,2, function(x) gsub("[^.]__", "", x))) # remove leading letters and __ with gsub
head(taxa.print)

#### Extracting the DADA2 data from R ####

# giving our seq headers more manageable names (ASV_1, ASV_2...)
asv_seqs <- colnames(seqtab.nochim)
asv_headers <- vector(dim(seqtab.nochim)[2], mode="character")
head(seqtab.nochim)
for (i in 1:dim(seqtab.nochim)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
}
asv_fasta <- c(rbind(asv_headers, asv_seqs))

# making and writing out a fasta of our final ASV seqs:
write(asv_fasta, "ITS2_ASVs_dada2.fa") # write fasta file
write.table(asv_fasta,"ITS2_ASVs_dada2.txt",sep="\t",row.names=TRUE,col.names=TRUE)

# count table:
asv_tab <- t(seqtab.nochim)
row.names(asv_tab) <- sub(">", "", asv_headers)
write.table(asv_tab, "ITS2_ASVs_Counts_dada2.tsv", sep="\t", quote=F, col.names=NA)
write.table(asv_tab,"ITS2_ASVs_Counts_dada2.txt",sep="\t",row.names=TRUE,col.names=TRUE)

# tax table:
asv_tax <- taxa
row.names(asv_tax) <- sub(">", "", asv_headers)
write.table(asv_tax, "ITS2_ASVs_Taxonomy_dada2.tsv", sep="\t", quote=F, col.names=NA)
write.table(asv_tax,"ITS2_ASVs_Taxonomy_dada2.txt",sep="\t",row.names=TRUE,col.names=TRUE)

#### Save as R objects ####
saveRDS(asv_tax, file = "ITS2_ASVs_Taxonomy_dada2_Robject.rds", ascii = FALSE, version = NULL,
        compress = TRUE, refhook = NULL)
saveRDS(asv_tab, file = "ITS2_ASVs_Counts_dada2_Robject.rds", ascii = FALSE, version = NULL,
        compress = TRUE, refhook = NULL)
saveRDS(asv_fasta, file = "ITS2_ASV_Sequences_dada2_Robject.rds", ascii = FALSE, version = NULL,
        compress = TRUE, refhook = NULL)

## how to read in R objects for future
readRDS(file, refhook = NULL)
infoRDS(file)
