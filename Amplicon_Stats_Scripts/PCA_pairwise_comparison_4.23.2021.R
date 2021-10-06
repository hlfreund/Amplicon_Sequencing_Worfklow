#if (!requireNamespace("BiocManager", quietly = TRUE))
#install.packages("BiocManager")
#BiocManager::install(version = "3.10")
#BiocManager::install(c("phyloseq"))
library(phyloseq)
library(devtools)
#install.packages("ggplot2")
library (ggplot2)
library (vegan)
library(reshape2)
#install.packages("ggpubr")
library(ggpubr)
#install.packages("scales")
library(scales)
library(grid)
library(ape)
#install.packages("dplyr")
library(dplyr)
#install.packages("viridis")
library(viridis)
#install.packages("readxl")
library(readxl)
#BiocManager::install("metagenomeSeq")
library(metagenomeSeq)
#if (!requireNamespace("BiocManager", quietly = TRUE))
#install.packages("BiocManager")
#BiocManager::install("DESeq2")
library(DESeq2)
library(dplyr)
library(magrittr)
library(MASS)
library(ade4)
library("dendextend")
library("tidyr")
library("viridis")
library("reshape")
# Install
#install.packages("wesanderson")
# Load
library(wesanderson)
#devtools::install_github("katiejolly/nationalparkcolors")
library(nationalparkcolors)
library(shades)
#install.packages("tidyverse")
#install.packages("tidygraph")
#install.packages("ggraph")
#install.packages("RAM")
#install.packages("FactoMineR")
library(FactoMineR)
#library(tidyverse)
#library(tidygraph)
#library(ggraph)
#library(RAM)

getwd()

#### Import and reformat the data ####

### ITS1 now
its1_otus_counts_taxa<- as.data.frame(read.csv("data/filtered_table_FUNGI.csv"))

dim(its1_otus_counts_taxa)
tail(its1_otus_counts_taxa) ## make sure that OTUs are rows and samples are columns!!

empty_to_unknown<-function(df){
  for(i in 1:nrow(df)){
    for(j in 1:ncol(df)){
      if(df[i,j]==""){
        df[i,j]<-sub("^$", "Unknown", df[i,j])
      }
    }
  }
  return(df)
}

its1_otus_counts_taxa<-empty_to_unknown(its1_otus_counts_taxa)
tail(its1_otus_counts_taxa) ## make sure that the function worked
head(its1_otus_counts_taxa)

its1_otus_counts_taxa<-subset(its1_otus_counts_taxa, Kingdom!="Unknown") ## keep only bacteria and archaean -- drop Unknowns 

'Unknown' %in% its1_otus_counts_taxa # check if Chloroplast counts are still in df, should be false because they've been removed

dim(its1_otus_counts_taxa)

head(its1_otus_counts_taxa)

rownames(its1_otus_counts_taxa)<-its1_otus_counts_taxa$OTU_ID
head(its1_otus_counts_taxa)

its1_otu_tax_ID<-subset(its1_otus_counts_taxa, select=c(Kingdom, Phylum, Class, Order, Family, Genus, Species)) # subset only part of the metadata we need
head(its1_otu_tax_ID)
its1_otu_counts<-subset(its1_otus_counts_taxa, select=-c(OTU_ID, Kingdom, Phylum, Class, Order, Family, Genus, Species)) # subset only part of the metadata we need
head(its1_otu_counts) ## these are RAW otu counts!

its1_otu_table<-data.frame(t(its1_otu_counts)) ## a work around to get a table of samples as rows and OTU as columns for when you melt the data
head(its1_otu_table)
min(colSums(its1_otu_table)) ## if min is 1 or 0, drop singletons before moving forward
#its1_otu_table<-its1_otu_table[, which(colSums(its1_otu_table) > 1)] # drop 0s and singletons 

dim(its1_otu_table)
## For Vegan analyses, use sample/site x species/OTUs/ASVs matrix (its1_otu_table)
# in vegan: rows are samples or site IDs, and columns are OTUs/ASVs/Species

### 16S
bac_otus_counts_taxa<- as.data.frame(read.csv("data/filtered_table_16S_191117_HF.2.csv"))
dim(bac_otus_counts_taxa)
tail(bac_otus_counts_taxa) ## making sure that OTUs are rows and samples are columns

bac_otus_counts_taxa<-empty_to_unknown(bac_otus_counts_taxa)

tail(bac_otus_counts_taxa)

#bac_otus_counts_taxa<-subset(bac_otus_counts_taxa, Class!="Chloroplast") ## keep only bacteria -- exclude Chloroplast sequences
'Unknown' %in% bac_otus_counts_taxa # check if UKs counts are still in df, should be false because they've been removed

dim(bac_otus_counts_taxa)
### creating new OTU IDs for these OTUs that are less confusing
prefix<-"OTU" # prefix of names will be OTU
suffix<-seq(1:37785) # creating sequence of all OTUs to assign them each a unique #
length(suffix)
bac_otu_ids_clean<-paste(prefix, suffix, sep="") # create vector where you attach prefix and suffix together with no separator
rownames(bac_otus_counts_taxa)<-bac_otu_ids_clean # set new OTU IDs as rownames
head(bac_otus_counts_taxa)

bac_otus_counts_taxa<-cbind(bac_otu_ids_clean, bac_otus_counts_taxa) # attach new OTU IDs to dataframe
names(bac_otus_counts_taxa)[names(bac_otus_counts_taxa) == "bac_otu_ids_clean"] <- "OTU_ID" # rename new column

bac_otu_tax_ID<-subset(bac_otus_counts_taxa, select=c(OTU_ID, Kingdom, Phylum, Class, Order, Family, Genus, Species)) # subset only part of the metadata we need
head(bac_otu_tax_ID)
bac_otu_counts<-subset(bac_otus_counts_taxa, select=-c(OTUID, OTU_ID, Kingdom, Phylum, Class, Order, Family, Genus, Species)) # subset only part of the metadata we need
head(bac_otu_counts) ## these are RAW otu counts!

bac_otu_table<-data.frame(t(bac_otu_counts)) ## a work around to get a table of samples as rows and OTU as columns for when you melt the data
head(bac_otu_table) # should see long list of OTU IDs

dim(bac_otu_table)

min(colSums(bac_otu_table)) ## if min is 1 or 0, drop singletons before moving forward
#bac_otu_table<-bac_otu_table[, which(colSums(bac_otu_table) > 1)] # drop 0s and singletons 

## For Vegan analysis, use sample x species matrix (bac_otu_table)


#### Import metadata ####
metadata<-as.data.frame(read_excel("data/SierraMapEle2SH.xls"))

head(metadata)
tail(metadata)
metadata<-subset(metadata, select=-c(DateCode, BarcodeSequence, LinkerPrimerSequence, SiteCode, RepNum, SiteRep, TubeID, DateCode, Description)) # subset only part of the metadata we need
head(metadata)

rownames(metadata)<-metadata$SampleID

# checking that we have data for all 31 samples
dim(metadata)
dim(its1_otu_table)
dim(bac_otu_table)
metadata$SampleID %in% row.names(its1_otu_table)
metadata$SampleID %in% row.names(bac_otu_table)

row.names(metadata)
row.names(its1_otu_table)
row.names(bac_otu_table)

dust_comp<-read.csv("data/DCSierra.csv", header=TRUE, sep=",")
head(dust_comp)
dust_comp<-subset(dust_comp, select=c(SampleID, DustComp1,DustComp2)) # subset only part of the metadata we need
rownames(dust_comp)<-dust_comp$SampleID

metadata<-merge(metadata, dust_comp, by.x="SampleID", by.y="SampleID")
rownames(metadata)<-metadata$SampleID

## Reorder rows of metadata and all dfs
metadata=metadata[rownames(its1_otu_table),] ## reorder metadata to have same rows as original OTU table
# ^ this indexing method will only work if the two dfs have the same # of rows AND the same row names!
head(metadata)
bac_otu_table=bac_otu_table[rownames(metadata),] ## reorder bacterial OTU table to have same rows as ITS1 OTU table + metadata
row.names(bac_otu_table)

#### Relativize the OTU/ASV counts ####

## fungi first
class(its1_otu_table) # raw counts
dim(its1_otu_table)
row.names(its1_otu_table)

its1_table_RA<-data.frame(decostand(its1_otu_table, method="total", MARGIN=1, na.rm=TRUE))  
# relative abundance of taxa data where everything is divided by margin total (default MARGIN = 1 = rows) -- rows = samples
rowSums(its1_table_RA)


## bacterial/archaeal next
class(bac_otu_table) # raw counts
dim(bac_otu_table)
row.names(bac_otu_table)

bac_table_RA<-data.frame(decostand(bac_otu_table, method="total", MARGIN=1, na.rm=TRUE))  
# relative abundance of taxa data where everything is divided by margin total (default MARGIN = 1 = rows) -- rows = samples
rowSums(bac_table_RA)

#### Subsetting Metadata ####
# subset metadata based on categorical variables and continuous variables
# metadata<-lapply(metadata, function(x) if(is.numeric(x)){
#   scale(x, center=TRUE, scale=TRUE)
# } else x)
head(metadata)
meta_category<-subset(metadata, select=c(SampleID, Year, Month, Site, Elevation, DustComp1, DustComp2)) # subset only part of the metadata we need
head(meta_category)
meta_quant<-subset(metadata, select=-c(SampleID, Year, NdEp, A_87Sr86Sr, Month, Site, SampType, Elevation, DustComp1, DustComp2)) # subset only part of the metadata we need
head(meta_quant)
meta_sites<-subset(metadata, select=c(SampleID, Site))
head(meta_sites)

meta_quant<-sapply(meta_quant, as.numeric) # turns all columns in df into numeric data (for pcoa)
tail(meta_quant)
meta_quant[,1] %in% metadata$Cu # sanity check
# match(meta_quant[,1],metadata$Cu, nomatch = NA_integer_)
# match(x, y, nomatch=NA): returns a vector showing the POSITION of x's match in table y
# in this case since we are comparing just two columns, any repeats in x[i] found in y[j] will just show the same index
# e.g. Cu val is 45.0 throughout the column in both dfs, so anytime there is a 45 in x found it y, it will return the first index where that match is found
# ^^ can use this vector to then order a df or array, like below
meta_quant<-meta_quant[order(match(meta_quant[,1],metadata$Cu)),]
# ordering meta_quant based on position order of matches from meta_quant[,1] found in metadata$Cu
# cannot call column using $ for array, only df
class(meta_quant)
rownames(meta_quant)<-rownames(metadata) # can only do this if you're sure the arrays/dfs are in the same exact order
tail(meta_quant)

meta_quant.scale<-as.data.frame(scale(meta_quant))
head(meta_quant.scale)
meta_quant.scale$SampleID<-rownames(meta_quant.scale)
meta_qdat_scale<-merge(meta_quant.scale,meta_sites, by.x="SampleID", by.y="SampleID")
head(meta_qdat_scale)

#### PCoA from comp data ####
## PCoA. It should be reserved for cases where no Euclidean measure 
## is appropriate (e.g., lots of ## shared zeros) and you prefer to use a similarity measure.
## If the metric is non-euclidean (as in our case), then the PCoA may produce several negative eigenvalues in
## addition to the positive ones. In most applications, this does not affect the representation of the first
## several axes. You will still receive a warning message, though, when this occurs. You also will receive a
## warning about species scores not being available; there is a way to project weighted averages of species
## abundances on a PCoA plot using the function wascores

## each ordination axis has an eigenvalue, it is a measure of the strength of
## an axis, the amount of variation explained by the axis. We typically consider the first several axes, as they
## have greater eigenvalues and thus explain much of the variation in a dataset

### ITS1 PCOA first...
## compute a matrix of Bray-Curtis similarities among sites
its1.bray<-vegdist(its1_otu_table, method="bray")
# pcoa of Bray-Curtis dissimilarity matrix made w/ relativized abundance (site x species) data
its1.sq.bray<-sqrt(its1.bray)
its1.pcoa = pcoa(its1.sq.bray) 
its1.pcoa ## Check to see if negative eigenvalues affect the interpretation of the first several axes
# The proportion of variances explained is in its element values$Relative_eig

str(its1.pcoa)
biplot(its1.pcoa)
## project (scaled) chem data on the PCoA ordination
biplot(its1.pcoa, meta_qdat_scale[,2:9], main="ITS1 PCoA + Scaled Chemical Data")

### 16S PCOA next...
## compute a matrix of Bray-Curtis similarities among sites
bac.bray<-vegdist(bac_otu_table, method="bray")
# pcoa of Bray-Curtis dissimilarity matrix made w/ relativized abundance (site x species) data
bac.sq.bray<-sqrt(bac.bray)
is.euclid(bac.sq.bray)
bac.pcoa = pcoa(bac.sq.bray) 
bac.pcoa ## Check to see if negative eigenvalues affect the interpretation of the first several axes
# The proportion of variances explained is in its element values$Relative_eig
#cmdscale(bac.bray, k=(nrow(bac_otu_table)-1), eig=TRUE)

biplot(bac.pcoa)
## project (scaled) chem data on the PCoA ordination
biplot(bac.pcoa, meta_qdat_scale[,2:9], main="16S PCoA + Scaled Chemical Data")

## Save for comparisons of ITS1 vs 16S
#png(file="figures/its1_pcoa.biplot_5.2.21.png",width=800, height=600)
biplot(its1.pcoa, meta_qdat_scale[,2:9], main="ITS1 PCoA + Scaled Chemical Data")
dev.off()

#png(file="figures/bac_pcoa.biplot_5.2.21.png",width=800, height=600)
biplot(bac.pcoa, meta_qdat_scale[,2:9], main="16S PCoA + Scaled Chemical Data")
dev.off()

#### Extract axes from PCoAs for comparison ####
## unlike NMDS, we can use PCoA sores for other analyses

# ITS1 first
str(its1.pcoa)
sum(abs(its1.pcoa$vectors[,1]))
its1.pcoa$vectors[,1:2]
its1_pcoa.scores=as.data.frame(its1.pcoa$vectors[,1:2])
its1_pcoa.scores$SampleID<-rownames(its1_pcoa.scores)

## plot to check
plot(its1_pcoa.scores$Axis.2 ~ its1_pcoa.scores$Axis.1 )

# 16S next
str(bac.pcoa)
sum(abs(bac.pcoa$vectors[,1]))
bac.pcoa$vectors[,1:2]
bac_pcoa.scores<-as.data.frame(bac.pcoa$vectors[,1:2])
bac_pcoa.scores$SampleID<-rownames(bac_pcoa.scores)

## plot to check
plot(bac_pcoa.scores$Axis.2 ~ bac_pcoa.scores$Axis.1)
## these can be used as a predictor variable in other analyses. NMDS scores cannot!!!

#### Pairwise comparison of PCoA axes vs env. data ####
## Note to self: to analyze the PC axes that capture the top 50% of variation

## ITS1 first

its1_sc.meta<-merge(meta_qdat_scale, its1_pcoa.scores, by.x="SampleID", by.y="SampleID")
pairs(Axis.1~Cu+Fe+Mg+Mn+Ni+P+S+Zn,data=its1_sc.meta, lower.panel=NULL, main="ITS1 Axis 1 (PCoA) vs. Scaled Chem Data")
## save ITS1 pairwise comparisons
#png(file="figures/its1_pcoa_Axis1_pairwise_5.2.21.png",width=1000, height=800)
pairs(Axis.1~Cu+Fe+Mg+Mn+Ni+P+S+Zn,data=its1_sc.meta, lower.panel=NULL, main="ITS1 Axis 1 (PCoA) vs. Scaled Chem Data")
dev.off()

adonis(its1.bray~Cu*Fe*Mg*Mn, data = meta_qdat_scale, strata=meta_qdat_scale$Site, permutations = 999)
adonis(its1.bray~Ni*P*S*Zn, data = meta_qdat_scale, strata=meta_qdat_scale$Site, permutations = 999)

adonis2(its1.bray~Cu*Fe*Mg*Mn, data = meta_qdat_scale, permutations = 999, by=NULL)
adonis2(its1.bray~Ni*P*S*Zn, data = meta_qdat_scale, permutations = 999, by=NULL)
adonis2(its1.bray~Cu*Fe*Mg*Mn*Ni*P*S*Zn, data = meta_qdat_scale, permutations = 999, by=NULL)

## Models to consider...
adonis2(its1.bray~Cu*Fe*Mn*P, data = meta_qdat_scale, permutations = 999, by=NULL)
adonis2(its1.bray~Cu*Mn*P, data = meta_qdat_scale, permutations = 999, by=NULL)
adonis2(its1.bray~Cu*Fe*P, data = meta_qdat_scale, permutations = 999, by=NULL)
adonis2(its1.bray~Cu*P, data = meta_qdat_scale, permutations = 999, by=NULL)


## 16S next...
bac_sc.meta<-merge(meta_qdat_scale, bac_pcoa.scores, by.x="SampleID", by.y="SampleID")
pairs(Axis.1~Cu*Fe*Mg*Mn*Ni*P*S*Zn,data=bac_sc.meta, lower.panel=NULL, main="16S Axis 1 (PCoA) vs. Scaled Chem Data")
# save 16S pairwise comparisons
png(file="figures/16s_pcoa_Axis1_pairwise_5.2.21.png",width=1000, height=800)
pairs(Axis.1~Cu*Fe*Mg*Mn*Ni*P*S*Zn,data=bac_sc.meta, lower.panel=NULL, main="16S Axis 1 (PCoA) vs. Scaled Chem Data")
dev.off()
png(file="figures/16s_pcoa_Axis2_pairwise_5.2.21.png",width=1000, height=800)
pairs(Axis.2~Cu*Fe*Mg*Mn*Ni*P*S*Zn,data=bac_sc.meta, lower.panel=NULL, main="16S Axis 1 (PCoA) vs. Scaled Chem Data")
dev.off()

adonis(bac.bray~Cu*Fe*Mg*Mn, data = meta_qdat_scale, strata=meta_qdat_scale$Site, permutations = 999)
adonis(bac.bray~Ni*P*S*Zn, data = meta_qdat_scale, strata=meta_qdat_scale$Site, permutations = 999)

adonis2(bac.bray~Cu*Fe*Mg*Mn, data = meta_qdat_scale, permutations = 999, by=NULL)
adonis2(bac.bray~Ni*P*S*Zn, data = meta_qdat_scale, permutations = 999,by=NULL)

# Models to consider...
adonis2(bac.bray~Cu, data = meta_qdat_scale, permutations = 999, by=NULL)


#### PCA ####
# Notes...
### PCoA w/ Euclidean distances = PCA!
### For categorical variables, since you can't do PCA...
### PCoA w/ Gower's dissimilarity -- can handle categorical and continuous data in same DF!! 
prcomp(its1_table_RA)
its1.pca <- prcomp(its1_table_RA)
summary(its1.pca)
## cumulative proportion is calculating the axis' or axes' proportion of variance
# PC1 + PC2 = 65.18% of variation
its1.PC1 = data.frame(scores(its1.pca, choices=c(1), display=c("sites")))
head(its1.PC1)
##changing the choices gives us different axes
its1.PC2 = data.frame(scores(its1.pca, choices=c(2), display=c("sites")))
head(its1.PC2)

# ## changing the display gives us the "rotation" or loadings
# its1.PC1s = scores(its1.pca, choices=c(1), display=c("species"))
# its1.PC2s = scores(its1.pca, choices=c(2), display=c("species"))

bac.pca <- prcomp(bac_table_RA)
summary(bac.pca)
## cumulative proportion is calculating the axis' or axes' proportion of variance
# PC1 + PC2 + PC3 + PC4 = 57.99% variation total

bac.PC1 = data.frame(scores(bac.pca, choices=c(1), display=c("sites")))
head(bac.PC1)
##changing the choices gives us different axes
bac.PC2 = data.frame(scores(bac.pca, choices=c(2), display=c("sites")))
head(bac.PC2)
bac.PC3 = data.frame(scores(bac.pca, choices=c(3), display=c("sites")))
head(bac.PC3)
bac.PC4 = data.frame(scores(bac.pca, choices=c(4), display=c("sites")))
head(bac.PC4)

#### Pairwise comparison of ITS1 PC axes vs env. data ####
## Note to self: to analyze the PC axes that capture the top 50% of variation

## ITS1 first
its1.PCs<-cbind(its1.PC1, its1.PC2)
its1.PCs$SampleID<-rownames(its1.PCs)
its1_PC.meta<-merge(meta_qdat_scale, its1.PCs, by.x="SampleID", by.y="SampleID")
pairs(PC1~Cu+Fe+Mg+Mn+Ni+P+S+Zn,data=its1_PC.meta, lower.panel=NULL, main="ITS1 PC1 (PCA) vs. Scaled Chem Data")
pairs(PC2~Cu+Fe+Mg+Mn+Ni+P+S+Zn,data=its1_PC.meta, lower.panel=NULL, main="ITS1 PC2 (PCA) vs. Scaled Chem Data")

## save ITS1 pairwise comparisons
#png(file="figures/its1_PC1_pairwise_5.2.21.png",width=1000, height=800)
pairs(PC1~Cu+Fe+Mg+Mn+Ni+P+S+Zn,data=its1_PC.meta, lower.panel=NULL, main="ITS1 PC1 (PCA) vs. Scaled Chem Data")
dev.off()
#png(file="figures/its1_PC1_pairwise_5.2.21.png",width=1000, height=800)
pairs(PC2~Cu+Fe+Mg+Mn+Ni+P+S+Zn,data=its1_PC.meta, lower.panel=NULL, main="ITS1 PC2 (PCA) vs. Scaled Chem Data")
dev.off()

# PERMANOVA check
adonis(PC1~Cu+Fe+Mg+Mn, data = its1_PC.meta, strata=its1_PC.meta$Site, permutations = 999)
adonis(its1.PC2~Ni+P+S+Zn, data = its1_PC.meta, strata=its1_PC.meta$Site, permutations = 999)

adonis2(its1.PC1~Cu+Fe+Mg+Mn, data = its1_PC.meta, permutations = 999, by="terms")
adonis2(its1.PC~Ni+P+S+Zn, data = its1_PC.meta, permutations = 999, by="terms")

install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(pairwiseAdonis)

its1.pair.mod<-pairwise.adonis(its1.bray,meta_qdat_scale$Site, p.adjust.m='bonferroni') # shows us variation for each sample to see which ones are different
its1.pair.mod

### SELF REMINDER FOR R^2
### Coefficient of Determination, denoted R2 or r2 and pronounced "R squared"
### is the proportion of the variance in the dependent variable that is predictable from the independent variable(s)

### Pseudo F stat for PERMANOVA
### pseudo F-ratio: It compares the total sum of squared dissimilarities (or ranked dissimilarities) among objects belonging to different groups to that of objects belonging to the same group. 
### Larger F-ratios indicate more pronounced group separation, however, the significance of this ratio is usually of more interest than its magnitude.

#### Pairwise comparison of 16S PC axes vs env. data ####
## Note to self: to analyze the PC axes that capture the top 50% of variation

bac.PCs<-cbind(bac.PC1, bac.PC2, bac.PC3, bac.PC4)
bac.PCs$SampleID<-rownames(bac.PCs)
bac_PC.meta<-merge(meta_qdat_scale, bac.PCs, by.x="SampleID", by.y="SampleID")
pairs(PC1~Cu+Fe+Mg+Mn+Ni+P+S+Zn,data=bac_PC.meta, lower.panel=NULL, main="16S PC1 (PCA) vs. Scaled Chem Data")
pairs(PC2~Cu+Fe+Mg+Mn+Ni+P+S+Zn,data=bac_PC.meta, lower.panel=NULL, main="16S PC2 (PCA) vs. Scaled Chem Data")
pairs(PC3~Cu+Fe+Mg+Mn+Ni+P+S+Zn,data=bac_PC.meta, lower.panel=NULL, main="16S PC3 (PCA) vs. Scaled Chem Data")
pairs(PC4~Cu+Fe+Mg+Mn+Ni+P+S+Zn,data=bac_PC.meta, lower.panel=NULL, main="16S PC4 (PCA) vs. Scaled Chem Data")

## save 16S pairwise comparisons
#png(file="figures/16s_PC1_pairwise_5.2.21.png",width=1000, height=800)
pairs(PC1~Cu+Fe+Mg+Mn+Ni+P+S+Zn,data=bac_PC.meta, lower.panel=NULL, main="16S PC1 (PCA) vs. Scaled Chem Data")
dev.off()
#png(file="figures/16s_PC2_pairwise_5.2.21.png",width=1000, height=800)
pairs(PC2~Cu+Fe+Mg+Mn+Ni+P+S+Zn,data=bac_PC.meta, lower.panel=NULL, main="16S PC2 (PCA) vs. Scaled Chem Data")
dev.off()

# PERMANOVA check
adonis(bac.bray~Cu+Fe+Mg+Mn, data = meta_qdat_scale, strata=meta_qdat_scale$Site, permutations = 999)
adonis(bac.bray~Ni+P+S+Zn, data = meta_qdat_scale, strata=meta_qdat_scale$Site, permutations = 999)

adonis2(bac.bray~Cu+Fe+Mg+Mn, data = meta_qdat_scale, permutations = 999, by="terms")
adonis2(bac.bray~Ni+P+S+Zn, data = meta_qdat_scale, permutations = 999, by="terms")

install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(pairwiseAdonis)

bac.pair.mod<-pairwise.adonis(bac.bray,meta_qdat_scale$Site, p.adjust.m='bonferroni') # shows us variation for each sample to see which ones are different
bac.pair.mod

### SELF REMINDER FOR R^2
### Coefficient of Determination, denoted R2 or r2 and pronounced "R squared"
### is the proportion of the variance in the dependent variable that is predictable from the independent variable(s)

### Pseudo F stat for PERMANOVA
### pseudo F-ratio: It compares the total sum of squared dissimilarities (or ranked dissimilarities) among objects belonging to different groups to that of objects belonging to the same group. 
### Larger F-ratios indicate more pronounced group separation, however, the significance of this ratio is usually of more interest than its magnitude.
