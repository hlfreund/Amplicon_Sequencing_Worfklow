# This is pulled from Dr. Benjamin Callehan's DADA2 tutorial
## here is the link: https://benjjneb.github.io/dada2/tutorial.html
## here is a link to the publication describing the DADA2 algorithm: https://www.nature.com/articles/nmeth.3869#Sec2

setwd("/Volumes/HLF_SSD/Aronson_Lab_Data/XCZO_ITS1/xczo_fastq")
getwd()
#install.packages("devtools")
library("devtools")
#devtools::install_github("benjjneb/dada2", ref="taxonomy-segfault-fix")
packageVersion("dada2") # Should be 1.15.1
library(dada2)

path <- "/Volumes/HLF_SSD/Aronson_Lab_Data/XCZO_ITS1/xczo_fastq/cut_sequences" # CHANGE ME to the directory containing the fastq files after unzipping.
path <- "/bigdata/aronsonlab/shared/HannahFreund/MtStHelens_ITS_MM/Trimmed_Seqs/its2_R1_reads"
#list.files(path)

# Now we read in the names of the fastq files, and perform some string manipulation to get matched lists of the forward and reverse fastq files.
# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
#fnFs <- sort(list.files(path, pattern="_R1_001_cut.fastq", full.names = TRUE))


fnFs_cut2 <- sort(list.files(path2, pattern="_R1_001.fastq", full.names = TRUE))
#fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))
#fnFs
fnFs_cut2
saveRDS(fnFs_cut2, file = "ITS1_ASVs_unchanged_c2_Robject.rds", ascii = FALSE, version = NULL,
        compress = TRUE, refhook = NULL)
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
#sample.names <- sapply(strsplit(basename(fnFs), "_L"), `[`, 1) #pattern where you want to split the name
sample.names1 <- sapply(strsplit(basename(fnFs_cut2), "_L001"), `[`, 1) #pattern where you want to split the name
sample.names1
# Check Read Quality Profiles
plotQualityProfile(fnFs_cut2[1:2])

# Filter + Trim
## assign filenames for filtered fastq.gz files
# Place filtered files in filtered/ subdirectory
path2
#filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtFs_cut2 <- file.path(path2, "filtered", paste0(sample.names1, "_F_filt2.fastq.gz"))
#filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs_cut2) <- sample.names1
#names(filtRs) <- sample.names

## Standard filtering parameters: maxN=0 (DADA2 requires no Ns), truncQ=2, rm.phix=TRUE and maxEE=2. 
## The maxEE parameter sets the maximum number of “expected errors” allowed in a read, which is a better filter than simply averaging quality scores.
# * if you are only doing F reads, remove the "truncLen" command - truncLen=c(240,160) [for PE reads]
out <- filterAndTrim(fnFs_cut2, filtFs_cut2, maxN=0, maxEE=c(2), truncQ=2, rm.phix=TRUE,compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE
## If too few reads are passing the filter, consider relaxing maxEE, perhaps especially on the reverse reads (eg. maxEE=c(2,5))
# and reducing the truncLen to remove low quality tails. Remember though, when choosing truncLen for paired-end reads you must maintain overlap after truncation in order to merge them later
head(out)
saveRDS(out, file = "ITS1_ASVs_filtertrim_c2_Robject.rds", ascii = FALSE, version = NULL,
        compress = TRUE, refhook = NULL)
# Learn the error rates
errF_cut2 <- learnErrors(filtFs_cut2, multithread=TRUE)
# The learnErrors method learns this error model from the data by alternating estimation of the error rates and inference of sample composition until they converge on a jointly consistent solution. 
# As in many machine-learning (ML) problems, the algorithm must begin with an initial guess, for which the maximum possible error rates in this data are used (the error rates if only the most abundant sequence is correct and all the rest are errors)
#errR <- learnErrors(filtRs, multithread=TRUE)
saveRDS(errF_cut2, file = "ITS1_ASVs_errorrates_c2_Robject.rds", ascii = FALSE, version = NULL,
        compress = TRUE, refhook = NULL)
plotErrors(errF_cut2, nominalQ=TRUE) ## sanity check by visualizing estimated error rates -- should see error rates drop w/ increased quality
## Error rates for each possible transition: A2A (A -> A), A2C (A -> C), etc
## Points are observed error rates for each consensus quality score
## Black line shows estimated error rates after convergence of the ML algorithm
## Red line shows error rates expected under the nominal definition of the Q-score

# Sample Inference!

# this step (w/ dada command) won't work if you've lost samples due to not surpassing the filtering step
filtFs_cut2 <- filtFs_cut2[file.exists(filtFs_cut2)] # removes files that were not included in output because 0 reads passed filter step
saveRDS(filtFs_cut2, file = "ITS1_ASVs_filtered_c2_Robject.rds", ascii = FALSE, version = NULL,
        compress = TRUE, refhook = NULL)
dadaFs_cut2 <- dada(filtFs_cut2, err=errF_cut2, multithread=TRUE) # pool = True will pool samples together before sample inference
#dadaRs <- dada(filtRs, err=errR, multithread=TRUE)
dadaFs_cut2[[1]]
dadaFs_cut2[1] # Returns first section of dada-class object {one sample}
saveRDS(dadaFs_cut2, file = "ITS1_ASVs_sampleinference.dada_c2_Robject.rds", ascii = FALSE, version = NULL,
        compress = TRUE, refhook = NULL)
# Merge Paired End Reads (obtain full denoised seqs)
## Merging is performed by aligning the denoised forward reads with the reverse-complement of the corresponding denoised reverse reads, and then constructing the merged “contig” sequences. 
## By default, merged sequences are only output if the forward and reverse reads overlap by at least 12 bases, and are identical to each other in the overlap region (but these conditions can be changed via function arguments)
#mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)

# Inspect the merger data.frame from the first sample
#head(mergers[[1]])

# Construct Sequence Table (ASVs)!
## sequence table is a matrix of rows (samples) and columns (ASVs)
seqtab_c2 <- makeSequenceTable(dadaFs_cut2)
dim(seqtab_c2)
#can replace dadaFs w "mergers" object if you are using PE reads that have been merged

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab_c2))) # Multiple length reads after cutadapt - ITS1 are of variable length (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5309391/)
dim(seqtab_c2)
# Sequences that are much longer or shorter than expected may be the result of non-specific priming...#
# You can remove non-target-length sequences from your sequence table vvvv
seqtab2_c2 <- seqtab_c2[,nchar(colnames(seqtab_c2)) %in% 210:260] # specifying length between 210 - 2 ## see paper above, average ITS1 length on UNITE is 227
saveRDS(seqtab2_c2, file = "ITS1_ASVs_withchims_c2_Robject.rds", ascii = FALSE, version = NULL,
        compress = TRUE, refhook = NULL)

# Remove chimeras
seqtab.nochim_c2 <- removeBimeraDenovo(seqtab2_c2, method="consensus", multithread=TRUE, verbose=TRUE)
# Chimeric sequences are identified if they can be exactly reconstructed by combining a left-segment and a right-segment from two more abundant “parent” sequences
dim(seqtab.nochim_c2)
sum(seqtab.nochim_c2)/sum(seqtab2_c2) # comparing reads after chimera removal over total reads (after filtering)
# Most of your reads should remain after chimera removal (it is not uncommon for a majority of sequence variants to be removed though)

#filtFs <- filtFs[file.exists(filtFs)] # removes files that were not included in output because 0 reads passed filter step

# Track reads through the pipeline: see how many reads made it through each step of pipeline
getN <- function(x) sum(getUniques(x))
track_c2 <- cbind(out, sapply(dadaFs_cut2, getN), rowSums(seqtab.nochim_c2))
head(track_c2)
#track <- cbind(out, sapply(dadaFs, getN),rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs) ;  sapply(dadaRs, getN), sapply(mergers, getN), 
colnames(track_c2) <- c("input", "filtered", "denoisedF", "nonchimera") # remove whichever labels you didn't include
rownames(track_c2) <- sample.names1
head(track_c2)
beep()
saveRDS(seqtab.nochim_c2, file = "ITS1_ASVs_nochims_c2_Robject.rds", ascii = FALSE, version = NULL,
        compress = TRUE, refhook = NULL)
# If a majority of reads failed to merge, you may need to revisit the truncLen parameter used in the filtering step and make sure that the truncated reads span your amplicon. 
# If a majority of reads were removed as chimeric, you may need to revisit the removal of primers, as the ambiguous nucleotides in unremoved primers interfere with chimera identification

# Assign Taxonomy! May have to untar/zip file to run this step

taxa <- assignTaxonomy(seqtab.nochim, "/Volumes/HLF_SSD/Aronson_Lab_Data/sh_general_release_dynamic_04.02.2020.fasta", multithread=TRUE)
#taxa <- assignTaxonomy(seqtab.nochim, "/bigdata/aronsonlab/shared/RefDBs/UNITE/sh_general_release_dynamic_04.02.2020.fasta", multithread=TRUE)

# The dada2 package also implements a method to make species level assignments based on exact matching between ASVs and sequenced reference strains..
# Currently, species-assignment training fastas are available for the Silva and RDP 16S databases. To follow the optional species addition step, download the silva_species_assignment_v132.fa.gz file, and place it in the directory with the fastq files.
#taxa <- addSpecies(taxa, "~/tax/silva_species_assignment_v132.fa.gz")

# Looking at the taxonomic assingments

taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)
head
# Evaluate Accuracy w/ mock communities if mock dataset available

unqs.mock <- seqtab.nochim["Mock",]
unqs.mock <- sort(unqs.mock[unqs.mock>0], decreasing=TRUE) # Drop ASVs absent in the Mock
cat("DADA2 inferred", length(unqs.mock), "sample sequences present in the Mock community.\n")

mock.ref <- getSequences(file.path(path, "HMP_MOCK.v35.fasta"))
match.ref <- sum(sapply(names(unqs.mock), function(x) any(grepl(x, mock.ref))))
cat("Of those,", sum(match.ref), "were exact matches to the expected reference sequences.\n")


# Extracting the DADA2 data from R

# giving our seq headers more manageable names (ASV_1, ASV_2...)
asv_seqs <- colnames(seqtab.nochim)
asv_headers <- vector(dim(seqtab.nochim)[2], mode="character")

for (i in 1:dim(seqtab.nochim)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
}

# making and writing out a fasta of our final ASV seqs:
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, "ITS1_ASVs_dada2.fa")
write.table(asv_fasta,"ITS1_ASVs_dada2.txt",sep="\t",row.names=TRUE,col.names=TRUE)

# count table:
asv_tab <- t(seqtab.nochim)
row.names(asv_tab) <- sub(">", "", asv_headers)
write.table(asv_tab, "ASVs_counts.tsv", sep="\t", quote=F, col.names=NA)
write.table(asv_tab,"ITS1_ASV_Counts_dada2.txt",sep="\t",row.names=TRUE,col.names=TRUE)

# tax table:
asv_tax <- taxa
row.names(asv_tax) <- sub(">", "", asv_headers)
write.table(asv_tax, "ASVs_taxonomy.tsv", sep="\t", quote=F, col.names=NA)
write.table(asv_tax,"ITS1_ASVs_taxID_dada2.txt",sep="\t",row.names=TRUE,col.names=TRUE)

# Save as R objects
saveRDS(asv_tax, file = "ITS1_ASVs_taxID_Robject.rds", ascii = FALSE, version = NULL,
        compress = TRUE, refhook = NULL)
saveRDS(asv_tab, file = "ITS1_ASVs_counts_Robject.rds", ascii = FALSE, version = NULL,
        compress = TRUE, refhook = NULL)
saveRDS(asv_fasta, file = "ITS1_ASVs_seqs_Robject.rds", ascii = FALSE, version = NULL,
        compress = TRUE, refhook = NULL)
readRDS(file, refhook = NULL)
infoRDS(file)