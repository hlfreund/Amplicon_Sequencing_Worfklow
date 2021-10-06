# This is pulled from Dr. Benjamin Callehan's DADA2 tutorial
## here is the link: https://benjjneb.github.io/dada2/tutorial.html
## here is a link to the publication describing the DADA2 algorithm: https://www.nature.com/articles/nmeth.3869#Sec2


# **** If running R on the cluster, do the following before you run DADA2:
## srun --partition=aronsonlab --mem=400gb --cpus-per-task 4 --ntasks 1 --time 1-00:00:00 --pty bash -l

## module load tmux
## tmux
## R or nvim DADA2_tutorial_pipeline.R
## then load libraries and check DADA2 package

getwd()
#install.packages("devtools")
#library("devtools")
#devtools::install_github("benjjneb/dada2", ref="taxonomy-segfault-fix")
packageVersion("dada2") # Should be 1.16

suppressPackageStartupMessages({
  library(dada2)
  library(tidyr)
  library(ggpubr)
  library(decontam)
})
### Set the Path ####
path <- "/PATH/TO/SEQUENCE/DATA/" # CHANGE ME to the directory containing the fastq files after unzipping.
list.files(path)

## If you are picking up where you left off, load your mydada_16S.Rdata file now
load("/PATH/TO/SEQUENCE/DATA/mydada_16S.V3V4.Rdata")

#### Read in sample/file names ####

# Now we read in the names of the fastq files, and perform some string manipulation to get matched lists of the forward and reverse fastq files.
# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq

fnFs <- sort(list.files(path, pattern="_R1_clean.fastq", full.names = TRUE))
#fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_clean.fastq", full.names = TRUE)); save.image(file = "mydada_16S.V3V4.Rdata")

fnFs

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_R1"), `[`, 1); save.image(file = "mydada_16S.V3V4.Rdata") #pattern where you want to split the name
sample.names

#### Check Read Quality Profiles ####

plot1<-plotQualityProfile(fnFs[1:2])
plot2<-plotQualityProfile(fnRs[1:2])
ggsave(plot1,filename = "16S.V3V4_pretrim_DADA2_F_quality.pdf", width=15, height=12, dpi=600)
ggsave(plot2,filename = "16S.V3V4_pretrim_DADA2_R_quality.pdf", width=15, height=12, dpi=600)

#### Filter + Trim ####

## assign filenames for filtered fastq.gz files
# Place filtered files in filtered/ subdirectory
path
filtFs <- file.path(path, "Filtered", paste0(sample.names, "_F_filtered.fastq.gz"))
filtRs <- file.path(path, "Filtered", paste0(sample.names, "_R_filtered.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names; save.image(file = "mydada_16S.V3V4.Rdata")

## Standard filtering parameters: maxN=0 (DADA2 requires no Ns), truncQ=2, rm.phix=TRUE and maxEE=2.
## The maxEE parameter sets the maximum number of “expected errors” allowed in a read, which is a better filter than simply averaging quality scores.
## * Refer to eestats2 results to see ideal read length for given expected error frequencies
# * if you are only doing F reads, remove the "truncLen" command - truncLen=c(240,160) [for PE reads]
# sometimes there is a trimLeft=15 argument here, but I removed this because I already trimmed my sequences with bbduk
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(280,230),
                     maxN=0, maxEE=c(2,3), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE); save.image(file = "mydada_16S.V3V4.Rdata") # On Windows set multithread=FALSE
# truncLen=c(240,230) -- trim F reads to 240 bp, trim R reads to 230 bp
# ^^^ check eestats results for truncLen limits
## Notes for trunc length of 2x300 PE reads: https://github.com/benjjneb/dada2/issues/236
## If too few reads are passing the filter, consider relaxing maxEE, perhaps especially on the reverse reads (eg. maxEE=c(2,5))
# and reducing the truncLen to remove low quality tails. Remember though, when choosing truncLen for paired-end reads you must maintain overlap after truncation in order to merge them later
head(out)

#### Learn the error rates ####

errF <- learnErrors(filtFs, multithread=TRUE); save.image(file = "mydada_16S.V3V4.Rdata")
errR <- learnErrors(filtRs, multithread=TRUE); save.image(file = "mydada_16S.V3V4.Rdata")
# The learnErrors method learns this error model from the data by alternating estimation of the error rates and inference of sample composition until they converge on a jointly consistent solution.
# As in many machine-learning (ML) problems, the algorithm must begin with an initial guess, for which the maximum possible error rates in this data are used (the error rates if only the most abundant sequence is correct and all the rest are errors)

plot_error<-plotErrors(errF, nominalQ=TRUE)## sanity check by visualizing estimated error rates -- should see error rates drop w/ increased quality
ggsave(plot_error,filename = "16S.V3V4_errormodel_DADA2.pdf", width=15, height=15, dpi=600)

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
dadaFs <- dada(filtFs, err=errF, multithread=TRUE, pool=TRUE); save.image(file = "mydada_16S.V3V4.Rdata") # pseudo pooling is computationally more efficient but similar in results to pooling; pool = True will pool samples together before sample inference
dadaRs <- dada(filtRs, err=errR, multithread=TRUE, pool=TRUE); save.image(file = "mydada_16S.V3V4.Rdata")

dadaFs[1] # Returns first section of dada-class object {one sample}
dadaRs[[1]]

#### Merge Paired End Reads (obtain full denoised seqs) ####

## Merging is performed by aligning the denoised forward reads with the reverse-complement of the corresponding denoised reverse reads, and then constructing the merged “contig” sequences.
## By default, merged sequences are only output if the forward and reverse reads overlap by at least 12 bases, and are identical to each other in the overlap region (but these conditions can be changed via function arguments)
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE); save.image(file = "mydada_16S.V3V4.Rdata")
# pick up where you left off vvv
load("/PATH/TO/SEQUENCE/DATA/mydada_16S.V3V4.Rdata")
# Inspect the merger data.frame from the first sample
head(mergers[[1]])

#### Construct Sequence Table (ASVs)! ####

## sequence table is a matrix of rows (samples) and columns (ASVs)
## PAIRED-END READS ##
seqtab <- makeSequenceTable(mergers); save.image(file = "mydada_16S.V3V4.Rdata")
dim(seqtab)
seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% 390:470]
# note on trimming V3-V4 regions: https://github.com/benjjneb/dada2/issues/1033#issuecomment-641281945
# *** note on how long V3-V4 region is: https://www.nature.com/articles/sdata20197
# another article mentioning V3V4 region length: https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-019-0743-1

write.table(seqtab,"16S.V3V4_Read_Length_Counts_dada2.txt",sep="\t",row.names=TRUE,col.names=TRUE); save.image(file = "mydada_16S.V3V4.Rdata")

## Forward Reads Only! ##
#seqtab <- makeSequenceTable(dadaFs)
#dim(seqtab)
# *** can replace dadaFs w "mergers" object if you are using PE reads that have been merged


#### Inspect distribution of sequence lengths ####
table(nchar(getSequences(seqtab)))
table(nchar(getSequences(seqtab2)))
dim(seqtab2)
#table(nchar(getSequences(seqtab2))) # Multiple length reads after cutadapt - ITS1 are of variable length (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5309391/)

# Sequences that are much longer or shorter than expected may be the result of non-specific priming...#
# You can remove non-target-length sequences from your sequence table vvvv
#seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% 250:256] # specifying length between 250 - 256

## Look at merged sequences in a plot -- see their distribution and frequency of sequences of certain length
compare_reads_plot = ggdensity(rowSums(seqtab2), fill = "blue4", alpha = 0.7); save.image(file = "mydada_16S.V3V4.Rdata")
pdf(file = "16S.V3V4_compare_plots.pdf", height = 10, width = 12) # Make empty pdf
compare_reads_plot # Add plot to empty pdf
dev.off()

# ## Look at F reads only
# compare_reads_F_its2=ggdensity(rowSums(seqtab2), fill = "red4", alpha = 0.7)
# pdf(file = "compare_plots2_F_its2_7.9.20.pdf", height = 10, width = 12) # Make empty pdf
# compare_reads_F_its2 # Add plot pdf
# dev.off()
#
# ## Compare Merged + F Reads in one plot
# pdf(file = "compare_plots2.pdf", height = 10, width = 12) # Make empty pdf
# ggdensity(gather(data.frame(Merged = rowSums(seqtab), FwdOnly = rowSums(seqtab2)), key = "Method", value = "ReadCount"), x = "ReadCount", fill = "Method", alpha = 0.7, palette = "tron", ggtheme = theme_dark())
# dev.off()
table(nchar(getSequences(seqtab2)))

#### Remove Chimeras ####
seqtab.nochim <- removeBimeraDenovo(seqtab2, method="consensus", multithread=TRUE, verbose=TRUE); save.image(file = "mydada_16S.V3V4.Rdata")
# Chimeric sequences are identified if they can be exactly reconstructed by combining a left-segment and a right-segment from two more abundant “parent” sequences
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab) # comparing reads after chimera removal over total reads (after filtering)
# Most of your reads should remain after chimera removal (it is not uncommon for a majority of sequence variants to be removed though)

#filtFs <- filtFs[file.exists(filtFs)] # removes files that were not included in output because 0 reads passed filter step

#### Track reads through the pipeline: see how many reads made it through each step of pipeline ####
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim)); save.image(file = "mydada_16S.V3V4.Rdata")
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs) ;  sapply(dadaRs, getN), sapply(mergers, getN),
head(track)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim") # remove whichever labels you didn't include
rownames(track) <- sample.names
head(track)
write.table(track,"16S.V3V4_tracking_reads_dada2.txt",sep="\t",row.names=TRUE,col.names=TRUE); save.image(file = "mydada_16S.V3V4.Rdata")

# If a majority of reads failed to merge, you may need to revisit the truncLen parameter used in the filtering step and make sure that the truncated reads span your amplicon.
# If a majority of reads were removed as chimeric, you may need to revisit the removal of primers, as the ambiguous nucleotides in unremoved primers interfere with chimera identification

#### Assign Taxonomy! May have to untar/zip file to run this step ####
## download trainig zip file: https://zenodo.org/record/3986799

# DADA2 package provides a native implementation of the naive Bayesian classifier method for this purpose.
# The assignTaxonomy function takes as input a set of sequences to be classified and a training set of reference sequences with known taxonomy, and outputs taxonomic assignments with at least minBoot bootstrap confidence.
taxa <- assignTaxonomy(seqtab.nochim, "/bigdata/aronsonlab/shared/DADA2_Silva_Files/silva_nr99_v138.1_wSpecies_train_set.fa.gz", multithread=TRUE)

# dada2 package also implements a method to make species level assignments based on exact matching between ASVs and sequenced reference strains. Recent analysis suggests that exact matching (or 100% identity) is the only appropriate way to assign species to 16S.V3V4 gene fragments.
# Currently, species-assignment training fastas are available for the Silva and RDP 16S.V3V4 databases. To follow the optional species addition step, download the silva_species_assignment_v132.fa.gz file, and place it in the directory with the fastq files.
#taxa <- addSpecies(taxa, "/bigdata/aronsonlab/shared/DADA2_Silva_Files/silva_species_assignment_v138.fa.gz"); save.image(file = "mydada_16S.V3V4.Rdata")

#### Looking at the taxonomic assingments ####

taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)
head
#### Evaluate Accuracy w/ mock communities if mock dataset available ####

unqs.mock <- seqtab.nochim["Mock",]
unqs.mock <- sort(unqs.mock[unqs.mock>0], decreasing=TRUE) # Drop ASVs absent in the Mock
cat("DADA2 inferred", length(unqs.mock), "sample sequences present in the Mock community.\n")

mock.ref <- getSequences(file.path(path, "HMP_MOCK.v35.fasta"))
match.ref <- sum(sapply(names(unqs.mock), function(x) any(grepl(x, mock.ref))))
cat("Of those,", sum(match.ref), "were exact matches to the expected reference sequences.\n")
save.image(file = "mydada_16S.V3V4.Rdata")

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
write(asv_fasta, "16S.V3V4_ASVs_dada2.fa") # write fasta file
write.table(asv_fasta,"16S.V3V4_ASVs_dada2.txt",sep="\t",row.names=TRUE,col.names=TRUE)

# count table:
asv_tab <- t(seqtab.nochim)
row.names(asv_tab) <- sub(">", "", asv_headers)
write.table(asv_tab, "16S.V3V4_ASVs_Counts_dada2.tsv", sep="\t", quote=F, col.names=NA)
write.table(asv_tab,"16S.V3V4_ASVs_Counts_dada2.txt",sep="\t",row.names=TRUE,col.names=TRUE)

# tax table:
asv_tax <- taxa
row.names(asv_tax) <- sub(">", "", asv_headers)
write.table(asv_tax, "16S.V3V4_ASVs_Taxonomy_dada2.tsv", sep="\t", quote=F, col.names=NA)
write.table(asv_tax,"16S.V3V4_ASVs_Taxonomy_dada2.txt",sep="\t",row.names=TRUE,col.names=TRUE)

#### Save as R objects ####
saveRDS(asv_tax, file = "16S.V3V4_ASVs_Taxonomy_dada2_Robject.rds", ascii = FALSE, version = NULL,
        compress = TRUE, refhook = NULL)
saveRDS(asv_tab, file = "16S.V3V4_ASVs_Counts_dada2_Robject.rds", ascii = FALSE, version = NULL,
        compress = TRUE, refhook = NULL)
saveRDS(asv_fasta, file = "16S.V3V4_ASV_Sequences_dada2_Robject.rds", ascii = FALSE, version = NULL,
        compress = TRUE, refhook = NULL)

## how to read in R objects for future
readRDS(file, refhook = NULL)
infoRDS(file)