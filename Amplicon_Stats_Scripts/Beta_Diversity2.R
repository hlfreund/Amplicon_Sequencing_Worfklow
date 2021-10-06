#if (!requireNamespace("BiocManager", quietly = TRUE))
#install.packages("BiocManager")
#BiocManager::install(version = "3.10")
#BiocManager::install(c("phyloseq"))
library(phyloseq)
#install.packages("ggplot2")
library (ggplot2)
library (vegan)
library(reshape2)
#install.packages("ggpubr")
library(ggpubr)
library(scales)
library(grid)
library(ape)
library(dplyr)
library(viridis)
library(readxl)
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
library(devtools)
library(wesanderson)
#devtools::install_github("katiejolly/nationalparkcolors")
library(nationalparkcolors)
library(shades)

library(FactoMineR)
#library(tidyverse)
#library(tidygraph)
#library(ggraph)
#library(RAM)

#setwd("/Volumes/HLF_SSD/Aronson_Lab_Data/SierraDust/New_Results_September2020")
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

tail(its1_otus_counts_taxa)

its1_otus_counts_taxa<-subset(its1_otus_counts_taxa, Kingdom!="Unknown") ## keep only bacteria and archaean -- drop Unknowns

'Unknown' %in% its1_otus_counts_taxa # check if Chloroplast counts are still in df, should be false because they've been removed

dim(its1_otus_counts_taxa)

head(its1_otus_counts_taxa)

rownames(its1_otus_counts_taxa)<-its1_otus_counts_taxa$OTU_ID
head(its1_otus_counts_taxa)

its1_otu_tax_ID<-subset(its1_otus_counts_taxa, select=c(OTU_ID, Kingdom, Phylum, Class, Order, Family, Genus, Species)) # subset only part of the metadata we need
head(its1_otu_tax_ID)
its1_otu_counts<-subset(its1_otus_counts_taxa, select=-c(OTU_ID, Kingdom, Phylum, Class, Order, Family, Genus, Species)) # subset only part of the metadata we need
head(its1_otu_counts) ## these are RAW otu counts!

its1_otu_table<-data.frame(t(its1_otu_counts)) ## a work around to get a table of samples as rows and OTU as columns for when you melt the data
head(its1_otu_table)

dim(its1_otu_table)

its1_otu_table<-its1_otu_table[,which(colSums(its1_otu_table) > 1)] # drop 0s and singletons
dim(its1_otu_table)

## For Vegan analyses, use sample/site x species/OTUs/ASVs matrix (its1_otu_table)
# in vegan: rows are samples or site IDs, and columns are OTUs/ASVs/Species

### 16S
bac_otus_counts_taxa<- as.data.frame(read.csv("data/filtered_table_16S_191117_HF.2.csv"))

dim(bac_otus_counts_taxa)
tail(bac_otus_counts_taxa) ## making sure that OTUs are rows and samples are columns

bac_otus_counts_taxa<-empty_to_unknown(bac_otus_counts_taxa)

tail(bac_otus_counts_taxa)

bac_otus_counts_taxa<-subset(bac_otus_counts_taxa, Kingdom!="Unknown") ## keep only bacteria and archaean -- drop Unknowns
#bac_otus_counts_taxa<-subset(bac_otus_counts_taxa, Phylum!="Unknown") ## keep only bacteria and archaean -- drop Unknowns

#bac_otus_counts_taxa<-subset(bac_otus_counts_taxa, Class!="Chloroplast") ## keep only bacteria -- exclude Chloroplast sequences
'Unknown' %in% bac_otus_counts_taxa # check if UKs counts are still in df, should be false because they've been removed

dim(bac_otus_counts_taxa)
### creating new OTU IDs for these OTUs that are less confusing
prefix<-"OTU" # prefix of names will be OTU
suffix<-seq(1:37785) # creating sequence of all OTUs to assign them each a unique #
length(suffix)
bac_otu_ids_clean<-paste(prefix, suffix, sep="") # create vector where you attach prefix and suffix together with no separator
rownames(bac_otus_counts_taxa)<-bac_otu_ids_clean # set new OTU IDs as rownames
head(bac_otu_ids_clean)
head(bac_otus_counts_taxa)

bac_otus_counts_taxa<-cbind(bac_otu_ids_clean, bac_otus_counts_taxa) # attach new OTU IDs to dataframe
names(bac_otus_counts_taxa)[names(bac_otus_counts_taxa) == "bac_otu_ids_clean"] <- "OTU_ID" # rename new column added to dataframe

bac_otu_tax_ID<-subset(bac_otus_counts_taxa, select=c(OTU_ID, Kingdom, Phylum, Class, Order, Family, Genus, Species)) # subset only part of the metadata we need
head(bac_otu_tax_ID)
bac_otu_counts<-subset(bac_otus_counts_taxa, select=-c(OTUID, OTU_ID, Kingdom, Phylum, Class, Order, Family, Genus, Species)) # subset only part of the metadata we need
head(bac_otu_counts) ## these are RAW otu counts!

bac_otu_table<-data.frame(t(bac_otu_counts)) ## a work around to get a table of samples as rows and OTU as columns for when you melt the data
head(bac_otu_table) # should see long list of OTU IDs

dim(bac_otu_table)

bac_otu_table<-bac_otu_table[,which(colSums(bac_otu_table) > 1)] # drop 0s and singletons
dim(bac_otu_table)

## For Vegan analysis, use sample x species matrix (bac_otu_table)

#### Import metadata ####
metadata<-as.data.frame(read_excel("data/SierraMapEle2SH.xls"))

head(metadata)
tail(metadata)
metadata<-subset(metadata, select=-c(DateCode, BarcodeSequence, LinkerPrimerSequence, SiteCode, RepNum, SiteRep, TubeID, DateCode, Description)) # subset only part of the metadata we need
head(metadata)

rownames(metadata)<-metadata$SampleID

dim(metadata)
dim(its1_otu_table)
dim(bac_otu_table)


# create colors for later figures

fair_cols <- melt(c('400'="#BF1B0B", '1100'="#FFC465", '2000'="#66ADE5", '2700'="#252A52"))
fair_cols$Elevation<-rownames(fair_cols)
colnames(fair_cols)[which(names(fair_cols) == "value")] <- "color"
fair_cols

metadata<-merge(metadata, fair_cols, by="Elevation")
head(metadata)
metadata$color <- as.character(metadata$color)
rownames(metadata)<-metadata$SampleID

## Reorder rows of metadata and all dfs
metadata=metadata[rownames(its1_otu_table),] ## reorder metadata to have same rows as original OTU table
# ^ this indexing method will only work if the two dfs have the same # of rows AND the same row names!
bac_otu_table=bac_otu_table[rownames(metadata),] ## reorder bacterial OTU table to have same rows as ITS1 OTU table + metadata

meta_category<-subset(metadata, select=c(SampleID, Year, Month, Site, Elevation)) # subset only part of the metadata we need
meta_quant<-subset(metadata, select=-c(Year, Month, Site, Elevation)) # subset only part of the metadata we need

#### Relative Abundance ####
## ITS1 first
class(its1_otu_table) # raw counts
dim(its1_otu_table)
row.names(its1_otu_table)

its1_RA_table<-data.frame(decostand(its1_otu_table, method="total", MARGIN=1, na.rm=TRUE))
# relative abundance of taxa data where everything is divided by margin total (default MARGIN = 1 = rows) -- rows = samples
rowSums(its1_RA_table)

## 16S next
class(bac_otu_table) # raw counts
dim(bac_otu_table)
row.names(bac_otu_table)

bac_RA_table<-data.frame(decostand(bac_otu_table, method="total", MARGIN=1, na.rm=TRUE))
# relative abundance of taxa data where everything is divided by margin total (default MARGIN = 1 = rows) -- rows = samples
rowSums(bac_RA_table)

#### Merge OTU tables, metadata, + taxonomy data ####

### ITS1 first
head(its1_otus_counts_taxa)
its1.all_melt<-melt(its1_otus_counts_taxa, id.vars = c("OTU_ID", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"))
head(its1.all_melt)
names(its1.all_melt)[which(names(its1.all_melt) == "variable")] <- "SampleID"
names(its1.all_melt)[which(names(its1.all_melt) == "value")] <- "Counts"
head(its1.all_melt)

head(metadata)

all_its1<-merge(its1.all_melt, metadata, by = "SampleID")

head(all_its1) # this dataframe contains EVERYTHING: OTU_IDs, Sample IDs, counts, metadata, and taxonomy info
# ^ use this for calculating relative abundance by x (i.e., elevation, month, year)

### 16S next
head(bac_otus_counts_taxa)
bac_otus_counts_taxa<-subset(bac_otus_counts_taxa, select=-c(OTUID)) # getting rid of original OTU column because it's irrelevant, renamed OTUs
bac.all_melt<-melt(bac_otus_counts_taxa, id.vars = c("OTU_ID", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"))
head(bac.all_melt)
names(bac.all_melt)[which(names(bac.all_melt) == "variable")] <- "SampleID"
names(bac.all_melt)[which(names(bac.all_melt) == "value")] <- "Counts"
head(bac.all_melt)

head(metadata)

all_bac<-merge(bac.all_melt, metadata, by = "SampleID")

head(all_bac) # this dataframe contains EVERYTHING: OTU_IDs, Sample IDs, counts, metadata, and taxonomy info
# ^ use this for calculating relative abundance by x (i.e., elevation, month, year)


#### ITS1: Relative Abundance by Taxa Level ####
head(all_its1)
### * below we use the dcast() function to "cast" the data into a wide format based on given elements (column names), taking sum of "Counts"
### * decostand(df, method="total) is the function (with argument total) used to get relative abundance of OTU table

## phylum ....
f.phyla_counts <- as.data.frame(dcast(all_its1, SampleID~Phylum, value.var="Counts", fun.aggregate=sum)) ###
head(f.phyla_counts) # counts by phyla per sample
rownames(f.phyla_counts)<-f.phyla_counts$SampleID
f.phyla_counts$SampleID<-NULL
f.phyla_RelAb<-data.frame(decostand(f.phyla_counts, method="total", MARGIN=1, na.rm=TRUE))  # relative abundance of taxa data where everything is divided by col total (b/c Margin=1 meaning rows == SAMPLES in this case)
rowSums(f.phyla_RelAb) # sanity check to make sure the transformation worked!
#f.phyla_RelAb$SampleID<-rownames(f.phyla_RelAb)
head(f.phyla_RelAb)
head(metadata)
f.phyla_RelAb=f.phyla_RelAb[rownames(metadata),] # reorder both dfs by row names

#write.csv(f.phyla_RelAb,"data/Sierra_ITS1_Phyla_Relative_Abundance.csv",row.names=TRUE)

## class...

f.class_counts <- as.data.frame(dcast(all_its1, SampleID~Class, value.var="Counts", fun.aggregate=sum)) ###
head(f.class_counts) # counts by class per sample
rownames(f.class_counts)<-f.class_counts$SampleID
f.class_counts$SampleID<-NULL
f.class_RelAb<-data.frame(decostand(f.class_counts, method="total", MARGIN=1, na.rm=TRUE))  # relative abundance of taxa data where everything is divided by col total (b/c Margin=1 meaning rows == SAMPLES in this case)
rowSums(f.class_RelAb) # sanity check to make sure the transformation worked!
#f.class_RelAb$SampleID<-rownames(f.class_RelAb)
head(f.class_RelAb)
#write.csv(f.class_RelAb,"data/Sierra_ITS1_Class_Relative_Abundance.csv",row.names=TRUE)
f.class_RelAb=f.class_RelAb[rownames(metadata),] # reorder both dfs by row names

## order ...

f.order_counts <- as.data.frame(dcast(all_its1, SampleID~Order, value.var="Counts", fun.aggregate=sum)) ###
head(f.order_counts) # counts by order per sample
rownames(f.order_counts)<-f.order_counts$SampleID
f.order_counts$SampleID<-NULL
f.order_RelAb<-data.frame(decostand(f.order_counts, method="total", MARGIN=1, na.rm=TRUE))  # relative abundance of taxa data where everything is divided by col total (b/c Margin=1 meaning rows == SAMPLES in this case)
rowSums(f.order_RelAb) # sanity check to make sure the transformation worked!
#f.order_RelAb$SampleID<-rownames(f.order_RelAb)
head(f.order_RelAb)
f.order_RelAb=f.order_RelAb[rownames(metadata),] # reorder both dfs by row names

#f.order_m<-melt(f.order_RelAb)

#### 16S: Relative Abundance by Taxa Level ####
head(all_bac)
### * below we use the dcast() function to "cast" the data into a wide format based on given elements (column names), taking sum of "Counts"
### * decostand(df, method="total) is the function (with argument total) used to get relative abundance of OTU table
## phylum ....
b.phyla_counts <- as.data.frame(dcast(all_bac, SampleID~Phylum, value.var="Counts", fun.aggregate=sum)) ###
head(b.phyla_counts) # counts by phyla per sample
rownames(b.phyla_counts)<-b.phyla_counts$SampleID
b.phyla_counts$SampleID<-NULL
b.phyla_RelAb<-data.frame(decostand(b.phyla_counts, method="total", MARGIN=1, na.rm=TRUE))  # relative abundance of taxa data where everything is divided by col total (b/c Margin=1 meaning rows == SAMPLES in this case)
rowSums(b.phyla_RelAb) # sanity check to make sure the transformation worked!
#b.phyla_RelAb$SampleID<-rownames(b.phyla_RelAb)
head(b.phyla_RelAb)
b.phyla_RelAb=b.phyla_RelAb[rownames(metadata),] # reorder both dfs by row names

#write.csv(b.phyla_RelAb,"data/Sierra_16S_Phyla_Relative_Abundance.csv",row.names=TRUE)

## class...

b.class_counts <- as.data.frame(dcast(all_bac, SampleID~Class, value.var="Counts", fun.aggregate=sum)) ###
head(b.class_counts) # counts by class per sample
rownames(b.class_counts)<-b.class_counts$SampleID
b.class_counts$SampleID<-NULL
b.class_RelAb<-data.frame(decostand(b.class_counts, method="total", MARGIN=1, na.rm=TRUE))  # relative abundance of taxa data where everything is divided by col total (b/c Margin=1 meaning rows == SAMPLES in this case)
rowSums(b.class_RelAb) # sanity check to make sure the transformation worked!
#b.class_RelAb$SampleID<-rownames(b.class_RelAb)
head(b.class_RelAb)
b.class_RelAb=b.class_RelAb[rownames(metadata),] # reorder both dfs by row names

#write.csv(b.class_RelAb,"data/Sierra_16S_Class_Relative_Abundance.csv",row.names=TRUE)

## order ...

b.order_counts <- as.data.frame(dcast(all_bac, SampleID~Order, value.var="Counts", fun.aggregate=sum)) ###
head(b.order_counts) # counts by order per sample
rownames(b.order_counts)<-b.order_counts$SampleID
b.order_counts$SampleID<-NULL
b.order_RelAb<-data.frame(decostand(b.order_counts, method="total", MARGIN=1, na.rm=TRUE))  # relative abundance of taxa data where everything is divided by col total (b/c Margin=1 meaning rows == SAMPLES in this case)
rowSums(b.order_RelAb) # sanity check to make sure the transformation worked!
#b.order_RelAb$SampleID<-rownames(b.order_RelAb)
head(b.order_RelAb)
b.order_RelAb=b.order_RelAb[rownames(metadata),] # reorder both dfs by row names

#### Bray-Curtis Dissimilarity Matrices w/ Relativized Count Data ####
head(metadata)
metadata=metadata[rownames(its1_RA_table),] # reorder both dfs by row names
bac_RA_table=bac_RA_table[rownames(metadata),] # reorder both dfs by row names

# ITS1 first
its1_RA_bray<-vegdist(its1_RA_table,method="bray")     ####### distance matrix with RELATIVIZED data!!!
# **** in vegan: ROWS need to be SITES/samples; COLUMNS are SPECIES (OTUs, ASVs) -- vegdist!
class(its1_RA_bray) # needs to be class dist for pcoa()

# 16S next
bac_RA_bray<-vegdist(bac_RA_table,method="bray")     ####### distance matrix with RELATIVIZED data!!!
# **** in vegan: ROWS need to be SITES/samples; COLUMNS are SPECIES (OTUs, ASVs) -- vegdist!
class(bac_RA_bray) # needs to be class dist for pcoa()

#### Merge OTU tables, metadata, + taxonomy data ####

### ITS1 first
head(its1_otus_counts_taxa)
its1.all_melt<-melt(its1_otus_counts_taxa, id.vars = c("OTU_ID", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"))
head(its1.all_melt)
names(its1.all_melt)[which(names(its1.all_melt) == "variable")] <- "SampleID"
names(its1.all_melt)[which(names(its1.all_melt) == "value")] <- "Counts"
head(its1.all_melt)

head(metadata)

all_its1<-merge(its1.all_melt, metadata, by = "SampleID")

head(all_its1) # this dataframe contains EVERYTHING: OTU_IDs, Sample IDs, counts, metadata, and taxonomy info
# ^ use this for calculating relative abundance by x (i.e., elevation, month, year)

### 16S next
head(bac_otus_counts_taxa)
bac_otus_counts_taxa<-subset(bac_otus_counts_taxa, select=-c(OTUID)) # getting rid of original OTU column because it's irrelevant, renamed OTUs
bac.all_melt<-melt(bac_otus_counts_taxa, id.vars = c("OTU_ID", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"))
head(bac.all_melt)
names(bac.all_melt)[which(names(bac.all_melt) == "variable")] <- "SampleID"
names(bac.all_melt)[which(names(bac.all_melt) == "value")] <- "Counts"
head(bac.all_melt)

head(metadata)

all_bac<-merge(bac.all_melt, metadata, by = "SampleID")

head(all_bac) # this dataframe contains EVERYTHING: OTU_IDs, Sample IDs, counts, metadata, and taxonomy info
# ^ use this for calculating relative abundance by x (i.e., elevation, month, year)


#### ITS1 BetaDisper + PERMANOVA ####

# First by binned month
head(metadata)
metadata=metadata[rownames(its1_RA_table),] # reorder both dfs by row names
bac_RA_table=bac_RA_table[rownames(metadata),] # reorder both dfs by row names

# Evaluate homogeneity of group dispersion
#(analogous to Levene's test for equal variances). The vegan function for this test is “betadisper”:
anova(betadisper(its1_RA_bray, metadata$MonthBin)) # p = 0.06877
# This tells us that there is NOT a difference in between-group dispersions, which means that we CAN trust the results of an adonis (permutational anova) test on this
## aka not significant betadisper result -- groups are homogeneous in dispersion rather than being heterogeneous
## the more heterogeneous within-group results exhibit, the more likely that the significant PERMANOVA results are due to the within-group differences (aka heterogenous dispersion) rather than between-group differences
# ^ more info on this here: https://www.researchgate.net/post/Betadisper_and_adonis_in_R_Am_I_interpreting_my_output_correctly

## Significant differences in homogeneities can be tested using either parametric or permutational tests,
##and parametric post hoc contrasts can also be investigated:
permutest(betadisper(its1_RA_bray, metadata$MonthBin))
its1.hsd<-TukeyHSD(betadisper(its1_RA_bray, metadata$MonthBin))

plot(its1.hsd)
plot(betadisper(its1_RA_bray, metadata$MonthBin))
#plot(betadisper(its1_RA_bray, metadata$Site))
boxplot(betadisper(its1_RA_bray, metadata$MonthBin))

# PERMANOVA assumption of homogenous within-group disperions is met - meaning we can run PERMANOVA!
adonis2(its1_RA_table ~ MonthBin, data = metadata, permutations = 999, method="bray", by='terms') # looks for interactions between predictor variables
#adonis2(its1_RA_table ~ MonthBin, data = metadata, permutations = 999, method="bray", by=NULL) # looks at entire model
pf1<-ggtexttable(adonis2(its1_RA_table ~ MonthBin, data = metadata, permutations = 999, method="bray", by='terms'))

pf2<-ggtexttable(adonis2(its1_RA_table ~ MonthBin*Elevation, data = metadata, permutations = 999, method="bray", by='terms')) # looks for interactions between predictor variables
adonis2(its1_RA_table ~ MonthBin*Elevation, data = metadata, permutations = 999, method="bray", by=NULL) # looks at entire model

pf3<-ggtexttable(adonis2(its1_RA_table ~ MonthBin*Year, data = metadata, permutations = 999, method="bray", by='terms')) # looks for interactions between predictor variables
adonis2(its1_RA_table ~ MonthBin*Year, data = metadata, permutations = 999, method="bray", by=NULL) # looks at entire model

p.its1.table<-ggarrange(pf1, pf2, pf3, legend = "none", common.legend = FALSE, labels = c("A. MonthBin", "B. MonthBin * Elevation", "C. MonthBin * Year"),align="h", ncol = 1, nrow = 3)
ggsave(p.its1.table,filename = "results/ITS1_PERMANOVAs_MonthBin_Comparisons_8.11.21.pdf", width=10, height=10, dpi=600)

#These models provide the significance for statistical interactions and the main effects.
## An advantage of adonis2 is that we can also test for overall model fit, setting by=NULL

##random person on the internet to the rescue!
install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(pairwiseAdonis)

pair.mod<-pairwise.adonis(its1_RA_table,metadata$MonthBin, p.adjust.m='bonferroni') # shows us variation for each sample to see which ones are different
pair.mod # in this case, all sites are different (F, R2)

## PERMANOVA w/ Rel.Ab by phyla, class
pf.phy1<-ggtexttable(adonis2(f.phyla_RelAb ~ MonthBin, data = metadata, permutations = 999, method="bray", by="terms")) # looks at entire model
pf.phy2<-ggtexttable(adonis2(f.phyla_RelAb ~ MonthBin*Elevation, data = metadata, permutations = 999, method="bray", by='terms')) # looks for interactions between predictor variables
pf.phy3<-ggtexttable(adonis2(f.phyla_RelAb ~ MonthBin*Year, data = metadata, permutations = 999, method="bray", by='terms')) # looks for interactions between predictor variables

p.its1.phy.table<-ggarrange(pf.phy1, pf.phy2, pf.phy3, legend = "none", common.legend = FALSE, labels = c("A. (Phyla Rel.Ab) MonthBin", "B. (Phyla Rel.Ab) MonthBin * Elevation", "C. (Phyla Rel.Ab) MonthBin * Year"),align="h", ncol = 1, nrow = 3)
ggsave(p.its1.phy.table,filename = "results/ITS1_Phyla_PERMANOVAs_MonthBin_Comparisons_8.11.21.pdf", width=10, height=15, dpi=600)

pf.cls1<-ggtexttable(adonis2(f.class_RelAb ~ MonthBin, data = metadata, permutations = 999, method="bray", by="terms")) # looks at entire model
pf.cls2<-ggtexttable(adonis2(f.class_RelAb ~ MonthBin*Elevation, data = metadata, permutations = 999, method="bray", by='terms')) # looks for interactions between predictor variables
pf.cls3<-ggtexttable(adonis2(f.class_RelAb ~ MonthBin*Year, data = metadata, permutations = 999, method="bray", by='terms')) # looks for interactions between predictor variables

p.its1.cls.table<-ggarrange(pf.cls1, pf.cls2, pf.cls3, legend = "none", common.legend = FALSE, labels = c("A. (Class Rel.Ab) MonthBin", "B. (Class Rel.Ab) MonthBin * Elevation", "C. (Class Rel.Ab) MonthBin * Year"),align="h", ncol = 1, nrow = 3)
ggsave(p.its1.cls.table,filename = "results/ITS1_Class_PERMANOVAs_MonthBin_Comparisons_8.11.21.pdf", width=10, height=10, dpi=600)

## Elevation...

# Evaluate homogeneity of group dispersion
#(analogous to Levene's test for equal variances). The vegan function for this test is “betadisper”:
anova(betadisper(its1_RA_bray, metadata$Elevation)) # p = 0.5069
# This tells us that there is NOT a difference in between-group dispersions, which means that we CAN trust the results of an adonis (permutational anova) test on this
## aka not significant betadisper result -- groups are homogeneous in dispersion rather than being heterogeneous
## the more heterogeneous within-group results exhibit, the more likely that the significant PERMANOVA results are due to the within-group differences (aka heterogenous dispersion) rather than between-group differences
# ^ more info on this here: https://www.researchgate.net/post/Betadisper_and_adonis_in_R_Am_I_interpreting_my_output_correctly

## Significant differences in homogeneities can be tested using either parametric or permutational tests,
##and parametric post hoc contrasts can also be investigated:
permutest(betadisper(its1_RA_bray, metadata$Elevation)) # p = 0.519
its1.hsd1<-TukeyHSD(betadisper(its1_RA_bray, metadata$Elevation))

plot(its1.hsd1)
plot(betadisper(its1_RA_bray, metadata$Elevation))
#plot(betadisper(its1_RA_bray, metadata$Site))
boxplot(betadisper(its1_RA_bray, metadata$Elevation))

# PERMANOVA assumption of homogenous within-group disperions is not met - meaning we can run PERMANOVA!
adonis2(its1_RA_table ~ Elevation, data = metadata, permutations = 999, method="bray", by='terms') # looks for interactions between predictor variables
#These models provide the significance for statistical interactions and the main effects.
## An advantage of adonis2 is that we can also test for overall model fit, setting by=NULL



#### 16S BetaDisper + PERMANOVA ####

# First by binned month
head(metadata)
metadata=metadata[rownames(its1_RA_table),] # reorder both dfs by row names
bac_RA_table=bac_RA_table[rownames(metadata),] # reorder both dfs by row names

# Evaluate homogeneity of group dispersion
#(analogous to Levene's test for equal variances). The vegan function for this test is “betadisper”:
anova(betadisper(bac_RA_bray, metadata$MonthBin)) # p = 0.06877
# This tells us that there is NOT a difference in between-group dispersions, which means that we CAN trust the results of an adonis (permutational anova) test on this
## aka not significant betadisper result -- groups are homogeneous in dispersion rather than being heterogeneous
## the more heterogeneous within-group results exhibit, the more likely that the significant PERMANOVA results are due to the within-group differences (aka heterogenous dispersion) rather than between-group differences
# ^ more info on this here: https://www.researchgate.net/post/Betadisper_and_adonis_in_R_Am_I_interpreting_my_output_correctly

## Significant differences in homogeneities can be tested using either parametric or permutational tests,
##and parametric post hoc contrasts can also be investigated:
permutest(betadisper(bac_RA_bray, metadata$MonthBin))
bac.hsd<-TukeyHSD(betadisper(bac_RA_bray, metadata$MonthBin))

plot(bac.hsd)
plot(betadisper(bac_RA_bray, metadata$MonthBin))
#plot(betadisper(bac_RA_bray, metadata$Site))
boxplot(betadisper(bac_RA_bray, metadata$MonthBin))

# PERMANOVA assumption of homogenous within-group disperions is not met - meaning we can run PERMANOVA!
adonis2(bac_RA_table ~ MonthBin, data = metadata, permutations = 999, method="bray", by='terms') # looks for interactions between predictor variables
#adonis2(bac_RA_table ~ MonthBin, data = metadata, permutations = 999, method="bray", by=NULL) # looks at entire model
pb1<-ggtexttable(adonis2(bac_RA_table ~ MonthBin, data = metadata, permutations = 999, method="bray", by='terms'))

pb2<-ggtexttable(adonis2(bac_RA_table ~ MonthBin*Elevation, data = metadata, permutations = 999, method="bray", by='terms')) # looks for interactions between predictor variables
adonis2(bac_RA_table ~ MonthBin*Elevation, data = metadata, permutations = 999, method="bray", by=NULL) # looks at entire model

pb3<-ggtexttable(adonis2(bac_RA_table ~ MonthBin*Year, data = metadata, permutations = 999, method="bray", by='terms')) # looks for interactions between predictor variables
adonis2(bac_RA_table ~ MonthBin*Year, data = metadata, permutations = 999, method="bray", by=NULL) # looks at entire model

p.bac.table<-ggarrange(pb1, pb2, pb3, legend = "none", common.legend = FALSE, labels = c("A. MonthBin", "B. MonthBin * Elevation", "C. MonthBin * Year"),align="h", ncol = 1, nrow = 3)
ggsave(p.bac.table,filename = "results/16S_PERMANOVAs_MonthBin_Comparisons_8.11.21.pdf", width=10, height=10, dpi=600)

#These models provide the significance for statistical interactions and the main effects.
## An advantage of adonis2 is that we can also test for overall model fit, setting by=NULL

##random person on the internet to the rescue!
install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(pairwiseAdonis)

pair.mod<-pairwise.adonis(bac_RA_table,metadata$MonthBin, p.adjust.m='bonferroni') # shows us variation for each sample to see which ones are different
pair.mod # in this case, all sites are different (F, R2)

## PERMANOVA w/ Rel.Ab by phyla, class
pb.phy1<-ggtexttable(adonis2(b.phyla_RelAb ~ MonthBin, data = metadata, permutations = 999, method="bray", by="terms")) # looks at entire model
pb.phy2<-ggtexttable(adonis2(b.phyla_RelAb ~ MonthBin*Elevation, data = metadata, permutations = 999, method="bray", by='terms')) # looks for interactions between predictor variables
pb.phy3<-ggtexttable(adonis2(b.phyla_RelAb ~ MonthBin*Year, data = metadata, permutations = 999, method="bray", by='terms')) # looks for interactions between predictor variables

p.bac.phy.table<-ggarrange(pb.phy1, pb.phy2, pb.phy3, legend = "none", common.legend = FALSE, labels = c("A. (Phyla Rel.Ab) MonthBin", "B. (Phyla Rel.Ab) MonthBin * Elevation", "C. (Phyla Rel.Ab) MonthBin * Year"),align="h", ncol = 1, nrow = 3)
ggsave(p.bac.phy.table,filename = "results/16S_Phyla_PERMANOVAs_MonthBin_Comparisons_8.11.21.pdf", width=10, height=15, dpi=600)

pb.cls1<-ggtexttable(adonis2(b.class_RelAb ~ MonthBin, data = metadata, permutations = 999, method="bray", by="terms")) # looks at entire model
pb.cls2<-ggtexttable(adonis2(b.class_RelAb ~ MonthBin*Elevation, data = metadata, permutations = 999, method="bray", by='terms')) # looks for interactions between predictor variables
pb.cls3<-ggtexttable(adonis2(b.class_RelAb ~ MonthBin*Year, data = metadata, permutations = 999, method="bray", by='terms')) # looks for interactions between predictor variables

p.bac.cls.table<-ggarrange(pb.cls1, pb.cls2, pb.cls3, legend = "none", common.legend = FALSE, labels = c("A. (Class Rel.Ab) MonthBin", "B. (Class Rel.Ab) MonthBin * Elevation", "C. (Class Rel.Ab) MonthBin * Year"),align="h", ncol = 1, nrow = 3)
ggsave(p.bac.cls.table,filename = "results/16S_Class_PERMANOVAs_MonthBin_Comparisons_8.11.21.pdf", width=10, height=10, dpi=600)

## Elevation...

# Evaluate homogeneity of group dispersion
#(analogous to Levene's test for equal variances). The vegan function for this test is “betadisper”:
anova(betadisper(bac_RA_bray, metadata$Elevation)) # p = 0.5069
# This tells us that there is NOT a difference in between-group dispersions, which means that we CAN trust the results of an adonis (permutational anova) test on this
## aka not significant betadisper result -- groups are homogeneous in dispersion rather than being heterogeneous
## the more heterogeneous within-group results exhibit, the more likely that the significant PERMANOVA results are due to the within-group differences (aka heterogenous dispersion) rather than between-group differences
# ^ more info on this here: https://www.researchgate.net/post/Betadisper_and_adonis_in_R_Am_I_interpreting_my_output_correctly

## Significant differences in homogeneities can be tested using either parametric or permutational tests,
##and parametric post hoc contrasts can also be investigated:
permutest(betadisper(bac_RA_bray, metadata$Elevation)) # p = 0.519
bac.hsd1<-TukeyHSD(betadisper(bac_RA_bray, metadata$Elevation))

plot(bac.hsd1)
plot(betadisper(bac_RA_bray, metadata$Elevation))
#plot(betadisper(bac_RA_bray, metadata$Site))
boxplot(betadisper(bac_RA_bray, metadata$Elevation))

# PERMANOVA assumption of homogenous within-group disperions is not met - meaning we can run PERMANOVA!
adonis2(bac_RA_table ~ Elevation, data = metadata, permutations = 999, method="bray", by='terms') # looks for interactions between predictor variables
#These models provide the significance for statistical interactions and the main effects.
## An advantage of adonis2 is that we can also test for overall model fit, setting by=NULL


#### mantel correlation test ####

man1<-mantel(f.euc_dist,b.euc_dist)
man1$statistic
man1$signif
