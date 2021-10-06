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
'Unknown' %in% bac_otus_counts_taxa # check if UKs counts are still in df, should be false because they've been removed

'Chloroplast' %in% bac_otus_counts_taxa$Class
#bac_otus_counts_taxa<-subset(bac_otus_counts_taxa, Class!="Chloroplast") ## keep only bacteria -- exclude Chloroplast sequences from Cyano counts

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
write.csv(its1_RA_table,"data/Sierra_ITS1_Relative_Abundance.csv",row.names=TRUE)

## 16S next
class(bac_otu_table) # raw counts
dim(bac_otu_table)
row.names(bac_otu_table)

bac_RA_table<-data.frame(decostand(bac_otu_table, method="total", MARGIN=1, na.rm=TRUE))  
# relative abundance of taxa data where everything is divided by margin total (default MARGIN = 1 = rows) -- rows = samples
rowSums(bac_RA_table)

write.csv(bac_RA_table,"data/Sierra_16S_Relative_Abundance.csv",row.names=TRUE)

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
f.phyla_RelAb$SampleID<-rownames(f.phyla_RelAb)
head(f.phyla_RelAb)

write.csv(f.phyla_RelAb,"data/Sierra_ITS1_Phyla_Relative_Abundance.csv",row.names=TRUE)

f.phyla_m<-melt(f.phyla_RelAb)

head(f.phyla_m)
colnames(f.phyla_m)[which(names(f.phyla_m) == "variable")] <- "Phylum"
colnames(f.phyla_m)[which(names(f.phyla_m) == "value")] <- "Counts"
head(f.phyla_m) ## relative abundance based on sum of counts by phyla!

f.phyla_RA_meta<-merge(f.phyla_m,metadata, by="SampleID")
head(f.phyla_RA_meta) ## relative abundance based on sum of counts by phyla!

## class...

f.class_counts <- as.data.frame(dcast(all_its1, SampleID~Class, value.var="Counts", fun.aggregate=sum)) ###
head(f.class_counts) # counts by class per sample
rownames(f.class_counts)<-f.class_counts$SampleID
f.class_counts$SampleID<-NULL
f.class_RelAb<-data.frame(decostand(f.class_counts, method="total", MARGIN=1, na.rm=TRUE))  # relative abundance of taxa data where everything is divided by col total (b/c Margin=1 meaning rows == SAMPLES in this case)
rowSums(f.class_RelAb) # sanity check to make sure the transformation worked!
f.class_RelAb$SampleID<-rownames(f.class_RelAb)
head(f.class_RelAb)
write.csv(f.class_RelAb,"data/Sierra_ITS1_Class_Relative_Abundance.csv",row.names=TRUE)

f.class_m<-melt(f.class_RelAb)

head(f.class_m)
colnames(f.class_m)[which(names(f.class_m) == "variable")] <- "Class"
colnames(f.class_m)[which(names(f.class_m) == "value")] <- "Counts"
head(f.class_m) ## relative abundance based on sum of counts by class!

f.class_RA_meta<-merge(f.class_m,metadata, by="SampleID")
head(f.class_RA_meta) ## relative abundance based on sum of counts by class!

## order ...

f.order_counts <- as.data.frame(dcast(all_its1, SampleID~Order, value.var="Counts", fun.aggregate=sum)) ###
head(f.order_counts) # counts by order per sample
rownames(f.order_counts)<-f.order_counts$SampleID
f.order_counts$SampleID<-NULL
f.order_RelAb<-data.frame(decostand(f.order_counts, method="total", MARGIN=1, na.rm=TRUE))  # relative abundance of taxa data where everything is divided by col total (b/c Margin=1 meaning rows == SAMPLES in this case)
rowSums(f.order_RelAb) # sanity check to make sure the transformation worked!
f.order_RelAb$SampleID<-rownames(f.order_RelAb)
head(f.order_RelAb)
f.order_m<-melt(f.order_RelAb)

head(f.order_m)
colnames(f.order_m)[which(names(f.order_m) == "variable")] <- "Order"
colnames(f.order_m)[which(names(f.order_m) == "value")] <- "Counts"
head(f.order_m) ## relative abundance based on sum of counts by order!

f.order_RA_meta<-merge(f.order_m,metadata, by="SampleID")
head(f.order_RA_meta) ## relative abundance based on sum of counts by order!

## family ...

f.family_counts <- as.data.frame(dcast(all_its1, SampleID~Family, value.var="Counts", fun.aggregate=sum)) ###
head(f.family_counts) # counts by family per sample
rownames(f.family_counts)<-f.family_counts$SampleID
f.family_counts$SampleID<-NULL
f.family_RelAb<-data.frame(decostand(f.family_counts, method="total", MARGIN=1, na.rm=TRUE))  # relative abundance of taxa data where everything is divided by col total (b/c Margin=1 meaning rows == SAMPLES in this case)
rowSums(f.family_RelAb) # sanity check to make sure the transformation worked!
f.family_RelAb$SampleID<-rownames(f.family_RelAb)
head(f.family_RelAb)
f.family_m<-melt(f.family_RelAb)

head(f.family_m)
colnames(f.family_m)[which(names(f.family_m) == "variable")] <- "Family"
colnames(f.family_m)[which(names(f.family_m) == "value")] <- "Counts"
head(f.family_m) ## relative abundance based on sum of counts by family!

f.family_RA_meta<-merge(f.family_m,metadata, by="SampleID")
head(f.family_RA_meta) ## relative abundance based on sum of counts by family!

## genus ...

f.genus_counts <- as.data.frame(dcast(all_its1, SampleID~Genus, value.var="Counts", fun.aggregate=sum)) ###
head(f.genus_counts) # counts by genus per sample
rownames(f.genus_counts)<-f.genus_counts$SampleID
f.genus_counts$SampleID<-NULL
f.genus_RelAb<-data.frame(decostand(f.genus_counts, method="total", MARGIN=1, na.rm=TRUE))  # relative abundance of taxa data where everything is divided by col total (b/c Margin=1 meaning rows == SAMPLES in this case)
rowSums(f.genus_RelAb) # sanity check to make sure the transformation worked!
f.genus_RelAb$SampleID<-rownames(f.genus_RelAb)
head(f.genus_RelAb)
f.genus_m<-melt(f.genus_RelAb)

head(f.genus_m)
colnames(f.genus_m)[which(names(f.genus_m) == "variable")] <- "Genus"
colnames(f.genus_m)[which(names(f.genus_m) == "value")] <- "Counts"
head(f.genus_m) ## relative abundance based on sum of counts by genus!

f.genus_RA_meta<-merge(f.genus_m,metadata, by="SampleID")
head(f.genus_RA_meta) ## relative abundance based on sum of counts by genus!

## species ...

f.species_counts <- as.data.frame(dcast(all_its1, SampleID~Genus+Species, value.var="Counts", fun.aggregate=sum)) ###
head(f.species_counts) # counts by species per sample
rownames(f.species_counts)<-f.species_counts$SampleID
f.species_counts$SampleID<-NULL
f.species_RelAb<-data.frame(decostand(f.species_counts, method="total", MARGIN=1, na.rm=TRUE))  # relative abundance of taxa data where everything is divided by col total (b/c Margin=1 meaning rows == SAMPLES in this case)
rowSums(f.species_RelAb) # sanity check to make sure the transformation worked!
f.species_RelAb$SampleID<-rownames(f.species_RelAb)
head(f.species_RelAb)
f.species_m<-melt(f.species_RelAb)

head(f.species_m)
colnames(f.species_m)[which(names(f.species_m) == "variable")] <- "Genus_Species"
colnames(f.species_m)[which(names(f.species_m) == "value")] <- "Counts"
f.species_m$Genus_Species<-gsub("_", " ", f.species_m$Genus_Species) ## gsub is global sub (does not just remove first instance of pattern, but multiple)

head(f.species_m) ## relative abundance based on sum of counts by species!

f.species_RA_meta<-merge(f.species_m,metadata, by="SampleID")
head(f.species_RA_meta) ## relative abundance based on sum of counts by species!


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
b.phyla_RelAb$SampleID<-rownames(b.phyla_RelAb)
head(b.phyla_RelAb)
write.csv(b.phyla_RelAb,"data/Sierra_16S_Phyla_Relative_Abundance.csv",row.names=TRUE)

b.phyla_m<-melt(b.phyla_RelAb)

head(b.phyla_m)
colnames(b.phyla_m)[which(names(b.phyla_m) == "variable")] <- "Phylum"
colnames(b.phyla_m)[which(names(b.phyla_m) == "value")] <- "Counts"
head(b.phyla_m) ## relative abundance based on sum of counts by phyla!
b.phyla_m<-subset(b.phyla_m, SampleID!="Undetermined_S0") ## keep only bacteria -- exclude Chloroplast sequences

b.phyla_RA_meta<-merge(b.phyla_m,metadata, by="SampleID")
head(b.phyla_RA_meta) ## relative abundance based on sum of counts by phyla!

## class...

b.class_counts <- as.data.frame(dcast(all_bac, SampleID~Class, value.var="Counts", fun.aggregate=sum)) ###
head(b.class_counts) # counts by class per sample
rownames(b.class_counts)<-b.class_counts$SampleID
b.class_counts$SampleID<-NULL
b.class_RelAb<-data.frame(decostand(b.class_counts, method="total", MARGIN=1, na.rm=TRUE))  # relative abundance of taxa data where everything is divided by col total (b/c Margin=1 meaning rows == SAMPLES in this case)
rowSums(b.class_RelAb) # sanity check to make sure the transformation worked!
b.class_RelAb$SampleID<-rownames(b.class_RelAb)
head(b.class_RelAb)
write.csv(b.class_RelAb,"data/Sierra_16S_Class_Relative_Abundance.csv",row.names=TRUE)

b.class_m<-melt(b.class_RelAb)

head(b.class_m)
colnames(b.class_m)[which(names(b.class_m) == "variable")] <- "Class"
colnames(b.class_m)[which(names(b.class_m) == "value")] <- "Counts"
head(b.class_m) ## relative abundance based on sum of counts by class!
b.class_m<-subset(b.class_m, SampleID!="Undetermined_S0") ## keep only bacteria -- exclude Chloroplast sequences

b.class_RA_meta<-merge(b.class_m,metadata, by="SampleID")
head(b.class_RA_meta) ## relative abundance based on sum of counts by class!

## order ...

b.order_counts <- as.data.frame(dcast(all_bac, SampleID~Order, value.var="Counts", fun.aggregate=sum)) ###
head(b.order_counts) # counts by order per sample
rownames(b.order_counts)<-b.order_counts$SampleID
b.order_counts$SampleID<-NULL
b.order_RelAb<-data.frame(decostand(b.order_counts, method="total", MARGIN=1, na.rm=TRUE))  # relative abundance of taxa data where everything is divided by col total (b/c Margin=1 meaning rows == SAMPLES in this case)
rowSums(b.order_RelAb) # sanity check to make sure the transformation worked!
b.order_RelAb$SampleID<-rownames(b.order_RelAb)
head(b.order_RelAb)
b.order_m<-melt(b.order_RelAb)

head(b.order_m)
colnames(b.order_m)[which(names(b.order_m) == "variable")] <- "Order"
colnames(b.order_m)[which(names(b.order_m) == "value")] <- "Counts"
head(b.order_m) ## relative abundance based on sum of counts by order!
b.order_m<-subset(b.order_m, SampleID!="Undetermined_S0") ## keep only bacteria -- exclude Chloroplast sequences

b.order_RA_meta<-merge(b.order_m,metadata, by="SampleID")
head(b.order_RA_meta) ## relative abundance based on sum of counts by order!

## family ...

b.family_counts <- as.data.frame(dcast(all_bac, SampleID~Family, value.var="Counts", fun.aggregate=sum)) ###
head(b.family_counts) # counts by family per sample
rownames(b.family_counts)<-b.family_counts$SampleID
b.family_counts$SampleID<-NULL
b.family_RelAb<-data.frame(decostand(b.family_counts, method="total", MARGIN=1, na.rm=TRUE))  # relative abundance of taxa data where everything is divided by col total (b/c Margin=1 meaning rows == SAMPLES in this case)
rowSums(b.family_RelAb) # sanity check to make sure the transformation worked!
b.family_RelAb$SampleID<-rownames(b.family_RelAb)
head(b.family_RelAb)
b.family_m<-melt(b.family_RelAb)

head(b.family_m)
colnames(b.family_m)[which(names(b.family_m) == "variable")] <- "Family"
colnames(b.family_m)[which(names(b.family_m) == "value")] <- "Counts"
head(b.family_m) ## relative abundance based on sum of counts by family!
b.family_m<-subset(b.family_m, SampleID!="Undetermined_S0") ## keep only bacteria -- exclude Chloroplast sequences

b.family_RA_meta<-merge(b.family_m,metadata, by="SampleID")
head(b.family_RA_meta) ## relative abundance based on sum of counts by family!

## genus ...

b.genus_counts <- as.data.frame(dcast(all_bac, SampleID~Genus, value.var="Counts", fun.aggregate=sum)) ###
head(b.genus_counts) # counts by genus per sample
rownames(b.genus_counts)<-b.genus_counts$SampleID
b.genus_counts$SampleID<-NULL
b.genus_RelAb<-data.frame(decostand(b.genus_counts, method="total", MARGIN=1, na.rm=TRUE))  # relative abundance of taxa data where everything is divided by col total (b/c Margin=1 meaning rows == SAMPLES in this case)
rowSums(b.genus_RelAb) # sanity check to make sure the transformation worked!
b.genus_RelAb$SampleID<-rownames(b.genus_RelAb)
head(b.genus_RelAb)
b.genus_m<-melt(b.genus_RelAb)

head(b.genus_m)
colnames(b.genus_m)[which(names(b.genus_m) == "variable")] <- "Genus"
colnames(b.genus_m)[which(names(b.genus_m) == "value")] <- "Counts"
head(b.genus_m) ## relative abundance based on sum of counts by genus!
b.genus_m<-subset(b.genus_m, SampleID!="Undetermined_S0") ## keep only bacteria -- exclude Chloroplast sequences

b.genus_RA_meta<-merge(b.genus_m,metadata, by="SampleID")
head(b.genus_RA_meta) ## relative abundance based on sum of counts by genus!

## species ...

b.species_counts <- as.data.frame(dcast(all_bac, SampleID~Genus+Species, value.var="Counts", fun.aggregate=sum)) ###
head(b.species_counts) # counts by species per sample
rownames(b.species_counts)<-b.species_counts$SampleID
b.species_counts$SampleID<-NULL
b.species_RelAb<-data.frame(decostand(b.species_counts, method="total", MARGIN=1, na.rm=TRUE))  # relative abundance of taxa data where everything is divided by col total (b/c Margin=1 meaning rows == SAMPLES in this case)
rowSums(b.species_RelAb) # sanity check to make sure the transformation worked!
b.species_RelAb$SampleID<-rownames(b.species_RelAb)
head(b.species_RelAb)
b.species_m<-melt(b.species_RelAb)

head(b.species_m)
colnames(b.species_m)[which(names(b.species_m) == "variable")] <- "Genus_Species"
colnames(b.species_m)[which(names(b.species_m) == "value")] <- "Counts"
b.species_m$Genus_Species<-gsub("_", " ", b.species_m$Genus_Species) ## gsub is global sub (does not just remove first instance of pattern, but multiple)
head(b.species_m) ## relative abundance based on sum of counts by species!
b.species_m<-subset(b.species_m, SampleID!="Undetermined_S0") ## keep only bacteria -- exclude Chloroplast sequences

b.species_RA_meta<-merge(b.species_m,metadata, by="SampleID")
head(b.species_RA_meta) ## relative abundance based on sum of counts by species!

#### ITS1 Relative Abundance: by Elevation ####

head(all_its1)

# by phylum + elevation
its1_phy_elv<- as.data.frame(dcast(all_its1,Elevation~Phylum, value.var="Counts", fun.aggregate=sum)) ### 
head(its1_phy_elv) # counts by phyla + elevation
rownames(its1_phy_elv)<-as.factor(its1_phy_elv$Elevation)

its1_RA_phy.elv<-data.frame(decostand(its1_phy_elv[,-1], method="total", MARGIN=1, na.rm=TRUE))  
# relative abundance of taxa data where everything is divided by margin total (default MARGIN = 1 = rows) -- rows = samples
rowSums(its1_RA_phy.elv)
head(its1_RA_phy.elv)

# by class + elevation
its1_cls_elv <- as.data.frame(dcast(all_its1,Elevation~Class, value.var="Counts", fun.aggregate=sum)) ### 
head(its1_cls_elv) # counts by class + elevation
rownames(its1_cls_elv)<-its1_cls_elv$Elevation

its1_RA_cls.elv<-data.frame(decostand(its1_cls_elv[,-1], method="total", MARGIN=1, na.rm=TRUE))  
# relative abundance of taxa data where everything is divided by margin total (default MARGIN = 1 = rows) -- rows = samples
rowSums(its1_RA_cls.elv)
head(its1_RA_cls.elv)

#### 16s Relative Abundance: by Elevation ####
head(all_bac)

# by phylum + elevation
bac_phy_elv<- as.data.frame(dcast(all_bac,Elevation~Phylum, value.var="Counts", fun.aggregate=sum)) ### 
head(bac_phy_elv) # counts by phyla + elevation
rownames(bac_phy_elv)<-as.factor(bac_phy_elv$Elevation)

bac_RA_phy.elv<-data.frame(decostand(bac_phy_elv[,-1], method="total", MARGIN=1, na.rm=TRUE))  
# relative abundance of taxa data where everything is divided by margin total (default MARGIN = 1 = rows) -- rows = samples
rowSums(bac_RA_phy.elv)
head(bac_RA_phy.elv)

# by class + elevation
bac_cls_elv <- as.data.frame(dcast(all_bac,Elevation~Class, value.var="Counts", fun.aggregate=sum)) ### 
head(bac_cls_elv) # counts by class + elevation
rownames(bac_cls_elv)<-bac_cls_elv$Elevation

bac_RA_cls.elv<-data.frame(decostand(bac_cls_elv[,-1], method="total", MARGIN=1, na.rm=TRUE))  
# relative abundance of taxa data where everything is divided by margin total (default MARGIN = 1 = rows) -- rows = samples
rowSums(bac_RA_cls.elv)
head(bac_RA_cls.elv)


#### ITS1 Relative Abundance: by Elevation + Month ####

head(all_its1)

# by phylum + elevation
its1_phy_e.m<- as.data.frame(dcast(all_its1,Elevation+Month~Phylum, value.var="Counts", fun.aggregate=sum)) ### 
head(its1_phy_e.m)
rownames(its1_phy_e.m)<-interaction(its1_phy_e.m$Month, its1_phy_e.m$Elevation, sep=".")
head(its1_phy_e.m) 

its1_RA_phy_e.m<-data.frame(decostand(its1_phy_e.m[,-c(1:2)], method="total", MARGIN=1, na.rm=TRUE))  
# relative abundance of taxa data where everything is divided by margin total (default MARGIN = 1 = rows) -- rows = samples
rowSums(its1_RA_phy_e.m)
head(its1_RA_phy_e.m)

# by class + elevation
its1_cls_e.m <- as.data.frame(dcast(all_its1,Elevation+Month~Class, value.var="Counts", fun.aggregate=sum)) ### 
head(its1_cls_e.m) # counts by class + elevation
rownames(its1_cls_e.m)<-interaction(its1_cls_e.m$Month, its1_cls_e.m$Elevation, sep=".")

its1_RA_cls_e.m<-data.frame(decostand(its1_cls_e.m[,-c(1:2)], method="total", MARGIN=1, na.rm=TRUE))  
# relative abundance of taxa data where everything is divided by margin total (default MARGIN = 1 = rows) -- rows = samples
rowSums(its1_RA_cls_e.m)
head(its1_RA_cls_e.m)

#### 16s Relative Abundance: by Elevation + Month ####
head(all_bac)

# by phylum + elevation
bac_phy_e.m<- as.data.frame(dcast(all_bac,Elevation+Month~Phylum, value.var="Counts", fun.aggregate=sum)) ### 
head(bac_phy_e.m) # counts by phyla + elevation
rownames(bac_phy_e.m)<-interaction(bac_phy_e.m$Month, bac_phy_e.m$Elevation, sep=".")

bac_RA_phy_e.m<-data.frame(decostand(bac_phy_e.m[,-c(1:2)], method="total", MARGIN=1, na.rm=TRUE))  
# relative abundance of taxa data where everything is divided by margin total (default MARGIN = 1 = rows) -- rows = samples
rowSums(bac_RA_phy_e.m)
head(bac_RA_phy_e.m)

#### 16s Relative Abundance: by Month + Year####
head(all_bac)

# by phylum + elevation
bac_phy_m.y<- as.data.frame(dcast(all_bac,Month+Year~Phylum, value.var="Counts", fun.aggregate=sum)) ### 
head(bac_phy_m.y) # counts by phyla + elevation
rownames(bac_phy_m.y)<-interaction(bac_phy_m.y$Month, bac_phy_m.y$Year, sep=".")

bac_phy_RA_m.y<-data.frame(decostand(bac_phy_m.y[,-c(1:2)], method="total", MARGIN=1, na.rm=TRUE))  
# relative abundance of taxa data where everything is divided by margin total (default MARGIN = 1 = rows) -- rows = samples
rowSums(bac_phy_RA_m.y)
head(bac_phy_RA_m.y)


#### 16s Relative Abundance: by Elevation + Month + Year####
head(all_bac)

# by phylum + elevation
bac_phy_e.m.y<- as.data.frame(dcast(all_bac,Elevation+Month+Year~Phylum, value.var="Counts", fun.aggregate=sum)) ### 
head(bac_phy_e.m.y) # counts by phyla + elevation
head(bac_phy_e.m.y)

# 2014 first

bac_phy_em2014<-subset(bac_phy_e.m.y, Year=="2014")
head(bac_phy_em2014)
rownames(bac_phy_em2014)<-interaction(bac_phy_em2014$Month, bac_phy_em2014$Elevation, sep=".")
dim(bac_phy_em2014)

bac_phy_RA_em2014<-data.frame(decostand(bac_phy_em2014[,-c(1:3)], method="total", MARGIN=1, na.rm=TRUE))  
# relative abundance of taxa data where everything is divided by margin total (default MARGIN = 1 = rows) -- rows = samples
rowSums(bac_phy_RA_em2014)
head(bac_phy_RA_em2014)
bac_phy_RA_em2014$EvMo<-rownames(bac_phy_RA_em2014)

# 2015 next

bac_phy_em2015<-subset(bac_phy_e.m.y, Year=="2015")
head(bac_phy_em2015)
rownames(bac_phy_em2015)<-interaction(bac_phy_em2015$Month, bac_phy_em2015$Elevation, sep=".")
dim(bac_phy_em2015)

bac_phy_RA_em2015<-data.frame(decostand(bac_phy_em2015[,-c(1:3)], method="total", MARGIN=1, na.rm=TRUE))  
# relative abundance of taxa data where everything is divided by margin total (default MARGIN = 1 = rows) -- rows = samples
rowSums(bac_phy_RA_em2015)
head(bac_phy_RA_em2015)
bac_phy_RA_em2015$EvMo<-rownames(bac_phy_RA_em2015)

#### ITS1 Relative Abundance: by Binned Month ####

head(all_its1)

# by phylum + MonthBin
its1_phy_mbin<- as.data.frame(dcast(all_its1,MonthBin~Phylum, value.var="Counts", fun.aggregate=sum)) ### 
head(its1_phy_mbin) # counts by phyla + MonthBin
rownames(its1_phy_mbin)<-as.factor(its1_phy_mbin$MonthBin)

its1_RA_phy.mbin<-data.frame(decostand(its1_phy_mbin[,-1], method="total", MARGIN=1, na.rm=TRUE))  
# relative abundance of taxa data where everything is divided by margin total (default MARGIN = 1 = rows) -- rows = samples
rowSums(its1_RA_phy.mbin)
head(its1_RA_phy.mbin)

# by class + MonthBin
its1_cls_mbin <- as.data.frame(dcast(all_its1,MonthBin~Class, value.var="Counts", fun.aggregate=sum)) ### 
head(its1_cls_mbin) # counts by class + MonthBin
rownames(its1_cls_mbin)<-its1_cls_mbin$MonthBin

its1_RA_cls.mbin<-data.frame(decostand(its1_cls_mbin[,-1], method="total", MARGIN=1, na.rm=TRUE))  
# relative abundance of taxa data where everything is divided by margin total (default MARGIN = 1 = rows) -- rows = samples
rowSums(its1_RA_cls.mbin)
head(its1_RA_cls.mbin)

#### 16s Relative Abundance: by Binned Month ####
head(all_bac)

# by phylum + MonthBin
bac_phy_mbin<- as.data.frame(dcast(all_bac,MonthBin~Phylum, value.var="Counts", fun.aggregate=sum)) ### 
head(bac_phy_mbin) # counts by phyla + MonthBin
rownames(bac_phy_mbin)<-as.factor(bac_phy_mbin$MonthBin)

bac_RA_phy.mbin<-data.frame(decostand(bac_phy_mbin[,-1], method="total", MARGIN=1, na.rm=TRUE))  
# relative abundance of taxa data where everything is divided by margin total (default MARGIN = 1 = rows) -- rows = samples
rowSums(bac_RA_phy.mbin)
head(bac_RA_phy.mbin)

# by class + MonthBin
bac_cls_mbin <- as.data.frame(dcast(all_bac,MonthBin~Class, value.var="Counts", fun.aggregate=sum)) ### 
head(bac_cls_mbin) # counts by class + MonthBin
rownames(bac_cls_mbin)<-bac_cls_mbin$MonthBin

bac_RA_cls.mbin<-data.frame(decostand(bac_cls_mbin[,-1], method="total", MARGIN=1, na.rm=TRUE))  
# relative abundance of taxa data where everything is divided by margin total (default MARGIN = 1 = rows) -- rows = samples
rowSums(bac_RA_cls.mbin)
head(bac_RA_cls.mbin)


#### ITS1 Relative Abundance: by Binned Month + Year ####

head(all_its1)

# by phylum + MonthBin
its1_phy_mbin.yr<- as.data.frame(dcast(all_its1,MonthBin+Year~Phylum, value.var="Counts", fun.aggregate=sum)) ### 
head(its1_phy_mbin.yr)
rownames(its1_phy_mbin.yr)<-interaction(its1_phy_mbin.yr$Year, its1_phy_mbin.yr$MonthBin, sep=".")
head(its1_phy_mbin.yr) 

its1_RA_phy_mbin.yr<-data.frame(decostand(its1_phy_mbin.yr[,-c(1:2)], method="total", MARGIN=1, na.rm=TRUE))  
# relative abundance of taxa data where everything is divided by margin total (default MARGIN = 1 = rows) -- rows = samples
rowSums(its1_RA_phy_mbin.yr)
head(its1_RA_phy_mbin.yr)

# by class + MonthBin
its1_cls_mbin.yr <- as.data.frame(dcast(all_its1,MonthBin+Year~Class, value.var="Counts", fun.aggregate=sum)) ### 
head(its1_cls_mbin.yr) # counts by class + MonthBin
rownames(its1_cls_mbin.yr)<-interaction(its1_cls_mbin.yr$Year, its1_cls_mbin.yr$MonthBin, sep=".")

its1_RA_cls_mbin.yr<-data.frame(decostand(its1_cls_mbin.yr[,-c(1:2)], method="total", MARGIN=1, na.rm=TRUE))  
# relative abundance of taxa data where everything is divided by margin total (default MARGIN = 1 = rows) -- rows = samples
rowSums(its1_RA_cls_mbin.yr)
head(its1_RA_cls_mbin.yr)

#### 16S Relative Abundance: by Binned Month + Year ####
head(all_bac)

# by phylum + MonthBin
bac_phy_mbin.yr<- as.data.frame(dcast(all_bac,MonthBin+Year~Phylum, value.var="Counts", fun.aggregate=sum)) ### 
head(bac_phy_mbin.yr) # counts by phyla + MonthBin
rownames(bac_phy_mbin.yr)<-interaction(bac_phy_mbin.yr$Year, bac_phy_mbin.yr$MonthBin, sep=".")

bac_RA_phy_mbin.yr<-data.frame(decostand(bac_phy_mbin.yr[,-c(1:2)], method="total", MARGIN=1, na.rm=TRUE))  
# relative abundance of taxa data where everything is divided by margin total (default MARGIN = 1 = rows) -- rows = samples
rowSums(bac_RA_phy_mbin.yr)
head(bac_RA_phy_mbin.yr)

#### Checking out some color pallettes before we visualize ####
wes1<-wes_palette("Chevalier1")
wes2<-wes_palette("Moonrise3")
wes3<-wes_palette("IsleofDogs1")
wes4<-wes_palette("GrandBudapest1")
wes5<-wes_palette("GrandBudapest2")
#scale_fill_manual(values = wes_palette("IsleofDogs1"))

SM_pal <- park_palette("SmokyMountains") # create a palette and specify # of colors youw ant
Arc_pal <- park_palette("Arches") # create a palette and specify # of colors youw ant
CL_pal <- park_palette("CraterLake") # create a palette and specify # of colors youw ant
Sag_pal <- park_palette("Saguaro") # create a palette and specify # of colors youw ant
Aca_pal <- park_palette("Acadia") # create a palette and specify # of colors youw ant
DV_pal <- park_palette("DeathValley") # create a palette and specify # of colors youw ant
CI_pal <- park_palette("ChannelIslands") # create a palette and specify # of colors youw ant
Bad_pal <- park_palette("Badlands") # create a palette and specify # of colors youw ant
MR_pal <- park_palette("MtRainier") # create a palette and specify # of colors youw ant
HI_pal <- park_palette("Hawaii") # create a palette and specify # of colors youw ant

fair_cols <- c("#BF1B0B", "#FFC465", "#66ADE5", "#252A52")
names(fair_cols) <- letters[1:4]
fair_ramp <- scales::colour_ramp(fair_cols)
fair_sat <- saturation(fair_ramp, 1)

#### ITS1 Specific Taxa by Elevation ####

## Phylum First
head(its1_RA_phy.elv)
its1_RA_phy.elv$Elevation<-rownames(its1_RA_phy.elv)
its1.RA.phy.ev.m<-melt(its1_RA_phy.elv)
head(its1.RA.phy.ev.m)
names(its1.RA.phy.ev.m)[which(names(its1.RA.phy.ev.m) == "variable")] <- "Phylum"
names(its1.RA.phy.ev.m)[which(names(its1.RA.phy.ev.m) == "value")] <- "Rel.Ab"
its1.RA.phy.ev.m$Elevation<-factor(its1.RA.phy.ev.m$Elevation, levels = c("400","1100","2000","2700"))
its1.RA.elv.all<-merge(its1.RA.phy.ev.m, metadata, by = "Elevation")

its1_phy_RA.elev<-ggplot(its1.RA.phy.ev.m, aes(x=factor(Elevation), y=Rel.Ab, fill=Phylum))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  labs(title = "Fungal Phylum Relative Abundance by Elevation", x="Elevation", y="Relative Abundance", fill="Phylum")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15))
ggsave(its1_phy_RA.elev,filename = "figures/fungal_phyla_RA_elev_6.9.21.pdf", width=8, height=6, dpi=600)

## Class next
head(its1_RA_cls.elv)
its1_RA_cls.elv$Elevation<-rownames(its1_RA_cls.elv)
its1.RA.cls.ev.m<-melt(its1_RA_cls.elv)
head(its1.RA.cls.ev.m)
names(its1.RA.cls.ev.m)[which(names(its1.RA.cls.ev.m) == "variable")] <- "Class"
names(its1.RA.cls.ev.m)[which(names(its1.RA.cls.ev.m) == "value")] <- "Rel.Ab"
its1.RA.cls.ev.m$Elevation<-factor(its1.RA.cls.ev.m$Elevation, levels = c("400","1100","2000","2700"))

its1_cls_RA.elev<-ggplot(its1.RA.cls.ev.m, aes(x=factor(Elevation), y=Rel.Ab, fill=Class))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  labs(title = "Fungal Class Relative Abundance by Elevation", x="Elevation", y="Relative Abundance", fill="Class")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15))
ggsave(its1_cls_RA.elev,filename = "figures/fungal_class_RA_elev_6.9.21.pdf", width=8, height=6, dpi=600)


#### 16S Specific Taxa by Elevation ####

## Phylum First
head(bac_RA_phy.elv)
bac_RA_phy.elv$Elevation<-rownames(bac_RA_phy.elv)
bac.RA.phy.ev.m<-melt(bac_RA_phy.elv)
head(bac.RA.phy.ev.m)
names(bac.RA.phy.ev.m)[which(names(bac.RA.phy.ev.m) == "variable")] <- "Phylum"
names(bac.RA.phy.ev.m)[which(names(bac.RA.phy.ev.m) == "value")] <- "Rel.Ab"
bac.RA.phy.ev.m$Elevation<-factor(bac.RA.phy.ev.m$Elevation, levels = c("400","1100","2000","2700"))

bac_phy_RA.elev<-ggplot(bac.RA.phy.ev.m, aes(x=factor(Elevation), y=Rel.Ab, fill=Phylum))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  labs(title = "Microbial Phylum Relative Abundance by Elevation", x="Elevation", y="Relative Abundance", fill="Phylum")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15))
ggsave(bac_phy_RA.elev,filename = "figures/microb_phyla_RA_elev_6.9.21.pdf", width=10, height=8, dpi=600)

Proteo_RA.elev<-ggplot(bac.RA.phy.ev.m[which(bac.RA.phy.ev.m$Phylum=='Proteobacteria'),], aes(x=factor(Elevation), y=Rel.Ab, fill=Phylum))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+
  theme_classic()+labs(title = "Proteobacteria Relative Abundance by Elevation", x="Elevation", y="Relative Abundance", fill="Phylum")+scale_fill_manual(values = "cornflowerblue")+
  theme(legend.position="none",axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15)) +
  
ggsave(Proteo_RA.elev,filename = "figures/Proteobacteria_RA_Elevation_6.9.21.pdf", width=8, height=6, dpi=600)

Firm_RA.elev<-ggplot(bac.RA.phy.ev.m[which(bac.RA.phy.ev.m$Phylum=='Firmicutes'),], aes(x=factor(Elevation), y=Rel.Ab, fill=Phylum))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+
  theme_classic()+labs(title = "Firmicutes Relative Abundance by Elevation", x="Elevation", y="Relative Abundance", fill="Phylum")+scale_fill_manual(values = brightness("seagreen3", 0.8))+
  theme(legend.position="none",axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15))
ggsave(Firm_RA.elev,filename = "figures/Firmicutes_RA_Elevation_6.9.21.pdf", width=8, height=6, dpi=600)

Actino_RA.elev<-ggplot(bac.RA.phy.ev.m[which(bac.RA.phy.ev.m$Phylum=='Actinobacteria'),], aes(x=factor(Elevation), y=Rel.Ab, fill=Phylum))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+
  theme_classic()+labs(title = "Actinobacteria Relative Abundance by Elevation", x="Elevation", y="Relative Abundance", fill="Phylum")+scale_fill_manual(values = brightness("coral2", 1))+
  theme(legend.position="none",axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15))
ggsave(Actino_RA.elev,filename = "figures/Actinobacteria_RA_Elevation_6.9.21.pdf", width=8, height=6, dpi=600)

Cyano_RA.elev<-ggplot(bac.RA.phy.ev.m[which(bac.RA.phy.ev.m$Phylum=='Cyanobacteria'),], aes(x=factor(Elevation), y=Rel.Ab, fill=Phylum))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+
  theme_classic()+labs(title = "Cyanobacteria Relative Abundance by Elevation", x="Elevation", y="Relative Abundance", fill="Phylum")+scale_fill_manual(values = brightness("olivedrab4", 0.7))+
  theme(legend.position="none",axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15))
ggsave(Cyano_RA.elev,filename = "figures/Cyanobacteria_RA_Elevation_6.9.21.pdf", width=8, height=6, dpi=600)

Acido_RA.elev<-ggplot(bac.RA.phy.ev.m[which(bac.RA.phy.ev.m$Phylum=='Acidobacteria'),], aes(x=factor(Elevation), y=Rel.Ab, fill=Phylum))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+
  theme_classic()+labs(title = "Acidobacteria Relative Abundance by Elevation", x="Elevation", y="Relative Abundance", fill="Phylum")+scale_fill_manual(values = saturation("tomato2", 0.6))+
  theme(legend.position="none",axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15))
ggsave(Acido_RA.elev,filename = "figures/Acidobacteria_RA_Elevation_6.9.21.pdf", width=8, height=6, dpi=600)

## Class next
head(bac_RA_cls.elv)
bac_RA_cls.elv$Elevation<-rownames(bac_RA_cls.elv)
bac.RA.cls.ev.m<-melt(bac_RA_cls.elv)
head(bac.RA.cls.ev.m)
names(bac.RA.cls.ev.m)[which(names(bac.RA.cls.ev.m) == "variable")] <- "Class"
names(bac.RA.cls.ev.m)[which(names(bac.RA.cls.ev.m) == "value")] <- "Rel.Ab"
bac.RA.cls.ev.m$Elevation<-factor(bac.RA.cls.ev.m$Elevation, levels = c("400","1100","2000","2700"))

bac_cls_RA.elev<-ggplot(bac.RA.cls.ev.m, aes(x=factor(Elevation), y=Rel.Ab, fill=Class))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  labs(title = "Fungal Class Relative Abundance by Elevation", x="Elevation", y="Relative Abundance", fill="Class")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15))
ggsave(bac_cls_RA.elev,filename = "figures/fungal_RA_elev_6.9.21.pdf", width=8, height=6, dpi=600)

#### 16S Specific Taxa by Elevation + Month ####

## Phylum 
head(bac_RA_phy_e.m)
bac_RA_phy_e.m$Elev.Month<-rownames(bac_RA_phy_e.m)
bac_RA_phy_e.m.m<-melt(bac_RA_phy_e.m)
head(bac_RA_phy_e.m.m)
names(bac_RA_phy_e.m.m)[which(names(bac_RA_phy_e.m.m) == "variable")] <- "Phylum"
names(bac_RA_phy_e.m.m)[which(names(bac_RA_phy_e.m.m) == "value")] <- "Rel.Ab"
bac_RA_phy_e.m.m$Elev.Month<-factor(bac_RA_phy_e.m.m$Elev.Month, levels = c("July.400","August.400", "October.400","July.1100","August.1100", "October.1100", 
                                                                            "July.2000","August.2000","October.2000","July.2700", "August.2700", "October.2700"))

bac_phy_RA.em1<-ggplot(bac_RA_phy_e.m.m, aes(x=Elev.Month, y=Rel.Ab, fill=Phylum))+geom_bar(stat="identity",colour="black")+scale_x_discrete(labels=c("July-400ft","Aug-400ft", "Oct-400ft","July-1100ft","Aug-1100ft", "Oct-1100ft", "July-2000ft", "Aug-2000ft", "Oct-2000ft", "July-2700ft", "Aug-2700ft", "Oct-2700ft"))+
  theme_classic()+labs(title = "Microbial Phylum Relative Abundance by Month + Elevation", x="Month & Elevation (ft)", y="Relative Abundance", fill="Phylum")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15))
ggsave(bac_phy_RA.em1,filename = "figures/microb_phyla_RA_Elev.month_6.10.21.pdf", width=20, height=10, dpi=600)

Proteo_RA.e.m.<-ggplot(bac_RA_phy_e.m.m[which(bac_RA_phy_e.m.m$Phylum=='Proteobacteria'),], aes(x=Elev.Month, y=Rel.Ab, fill=Phylum))+geom_bar(stat="identity",colour="black")+
  theme_classic()+labs(title = "Proteobacteria Relative Abundance by Month & Elevation", x="Month & Elevation (ft)", y="Relative Abundance", fill="Phylum")+scale_fill_manual(values = "cornflowerblue")+
  theme(legend.position="none",axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15)) +
  scale_x_discrete(labels=c("July-400ft","Aug-400ft", "Oct-400ft","July-1100ft","Aug-1100ft", "Oct-1100ft", "July-2000ft", "Aug-2000ft", "Oct-2000ft", "July-2700ft", "Aug-2700ft", "Oct-2700ft"))
ggsave(Proteo_RA.e.m.,filename = "figures/Proteobacteria_RA_Elev.Month_6.10.21.pdf", width=12, height=8, dpi=600)

Firm_RA.e.m<-ggplot(bac_RA_phy_e.m.m[which(bac_RA_phy_e.m.m$Phylum=='Firmicutes'),], aes(x=Elev.Month, y=Rel.Ab, fill=Phylum))+geom_bar(stat="identity",colour="black")+
  theme_classic()+labs(title = "Firmicutes Relative Abundance by Month & Elevation", x="Month & Elevation (ft)", y="Relative Abundance", fill="Phylum")+scale_fill_manual(values = brightness("seagreen3", 0.8))+
  theme(legend.position="none",axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15)) +
  scale_x_discrete(labels=c("July-400ft","Aug-400ft", "Oct-400ft","July-1100ft","Aug-1100ft", "Oct-1100ft", "July-2000ft", "Aug-2000ft", "Oct-2000ft", "July-2700ft", "Aug-2700ft", "Oct-2700ft"))
ggsave(Firm_RA.e.m,filename = "figures/Firmicutes_RA_Elev.Month_6.10.21.pdf", width=12, height=8, dpi=600)

Actino_RA.e.m<-ggplot(bac_RA_phy_e.m.m[which(bac_RA_phy_e.m.m$Phylum=='Actinobacteria'),], aes(x=Elev.Month, y=Rel.Ab, fill=Phylum))+geom_bar(stat="identity",colour="black")+
  theme_classic()+labs(title = "Actinobacteria Relative Abundance by Month & Elevation", x="Month & Elevation (ft)", y="Relative Abundance", fill="Phylum")+scale_fill_manual(values = brightness("coral2", 1))+
  theme(legend.position="none",axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15)) +
  scale_x_discrete(labels=c("July-400ft","Aug-400ft", "Oct-400ft","July-1100ft","Aug-1100ft", "Oct-1100ft", "July-2000ft", "Aug-2000ft", "Oct-2000ft", "July-2700ft", "Aug-2700ft", "Oct-2700ft"))
ggsave(Actino_RA.e.m,filename = "figures/Actinobacteria_RA_Elev.Month_6.10.21.pdf", width=12, height=8, dpi=600)

Cyano_RA.e.m<-ggplot(bac_RA_phy_e.m.m[which(bac_RA_phy_e.m.m$Phylum=='Cyanobacteria'),], aes(x=Elev.Month, y=Rel.Ab, fill=Phylum))+geom_bar(stat="identity",colour="black")+
  theme_classic()+labs(title = "Cyanobacteria Relative Abundance by Month & Elevation", x="Month & Elevation (ft)", y="Relative Abundance", fill="Phylum")+scale_fill_manual(values = brightness("olivedrab4", 0.7))+
  theme(legend.position="none",axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15)) +
  scale_x_discrete(labels=c("July-400ft","Aug-400ft", "Oct-400ft","July-1100ft","Aug-1100ft", "Oct-1100ft", "July-2000ft", "Aug-2000ft", "Oct-2000ft", "July-2700ft", "Aug-2700ft", "Oct-2700ft"))
ggsave(Cyano_RA.e.m,filename = "figures/Cyanobacteria_RA_Elev.Month_6.10.21.pdf", width=12, height=8, dpi=600)

Acido_RA.e.m<-ggplot(bac_RA_phy_e.m.m[which(bac_RA_phy_e.m.m$Phylum=='Acidobacteria'),], aes(x=Elev.Month, y=Rel.Ab, fill=Phylum))+geom_bar(stat="identity",colour="black")+
  theme_classic()+labs(title = "Acidobacteria Relative Abundance by Month & Elevation", x="Month & Elevation (ft)", y="Relative Abundance", fill="Phylum")+scale_fill_manual(values = saturation("tomato2", 0.6))+
  theme(legend.position="none",axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15)) +
  scale_x_discrete(labels=c("July-400ft","Aug-400ft", "Oct-400ft","July-1100ft","Aug-1100ft", "Oct-1100ft", "July-2000ft", "Aug-2000ft", "Oct-2000ft", "July-2700ft", "Aug-2700ft", "Oct-2700ft"))
ggsave(Acido_RA.e.m,filename = "figures/Acidobacteria_RA_Elev.Month_6.10.21.pdf", width=12, height=8, dpi=600)


#### 16S Specific Taxa by Month + Year ####

## Phylum 
head(bac_phy_RA_m.y)
bac_phy_RA_m.y$Month.Yr<-rownames(bac_phy_RA_m.y)
bac_phy_RA_m.y.m<-melt(bac_phy_RA_m.y)
head(bac_phy_RA_m.y.m)
names(bac_phy_RA_m.y.m)[which(names(bac_phy_RA_m.y.m) == "variable")] <- "Phylum"
names(bac_phy_RA_m.y.m)[which(names(bac_phy_RA_m.y.m) == "value")] <- "Rel.Ab"
bac_phy_RA_m.y.m$Month.Yr<-factor(bac_phy_RA_m.y.m$Month.Yr, levels = c("July.2014","August.2014", "October.2014","July.2015","October.2015"))

bac_phy_RA.my1<-ggplot(bac_phy_RA_m.y.m, aes(x=Month.Yr, y=Rel.Ab, fill=Phylum))+geom_bar(stat="identity",colour="black")+scale_x_discrete(labels=c("July 2014","August 2014", "October 2014","July 2015","October 2015"))+
  theme_classic()+labs(title = "Microbial Phylum Relative Abundance by Month & Year", x="Month & Year", y="Relative Abundance", fill="Phylum")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15))
ggsave(bac_phy_RA.my1,filename = "figures/microb_phyla_RA_month.year_6.11.21.pdf", width=15, height=8, dpi=600)

Proteo_RA.my1<-ggplot(bac_phy_RA_m.y.m[which(bac_phy_RA_m.y.m$Phylum=='Proteobacteria'),], aes(x=Month.Yr, y=Rel.Ab, fill=Phylum))+geom_bar(stat="identity",colour="black")+
  theme_classic()+labs(title = "Proteobacteria Relative Abundance by Month & Year", x="Month & Year", y="Relative Abundance", fill="Phylum")+scale_fill_manual(values = "cornflowerblue")+
  theme(legend.position="none",axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15)) +
  scale_x_discrete(labels=c("July 2014","August 2014", "October 2014","July 2015","October 2015"))
ggsave(Proteo_RA.my1,filename = "figures/Proteobacteria_RA_Month.Yr_6.11.21.pdf", width=12, height=8, dpi=600)

Firm_RA.my1<-ggplot(bac_phy_RA_m.y.m[which(bac_phy_RA_m.y.m$Phylum=='Firmicutes'),], aes(x=Month.Yr, y=Rel.Ab, fill=Phylum))+geom_bar(stat="identity",colour="black")+
  theme_classic()+labs(title = "Firmicutes Relative Abundance by Month & Year", x="Month & Year", y="Relative Abundance", fill="Phylum")+scale_fill_manual(values = brightness("seagreen3", 0.8))+
  theme(legend.position="none",axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15)) +
  scale_x_discrete(labels=c("July 2014","August 2014", "October 2014","July 2015","October 2015"))
ggsave(Firm_RA.my1,filename = "figures/Firmicutes_RA_Month.Yr_6.11.21.pdf", width=12, height=8, dpi=600)

Actino_RA.my1<-ggplot(bac_phy_RA_m.y.m[which(bac_phy_RA_m.y.m$Phylum=='Actinobacteria'),], aes(x=Month.Yr, y=Rel.Ab, fill=Phylum))+geom_bar(stat="identity",colour="black")+
  theme_classic()+labs(title = "Actinobacteria Relative Abundance by Month & Year", x="Month & Year", y="Relative Abundance", fill="Phylum")+scale_fill_manual(values = brightness("coral2", 1))+
  theme(legend.position="none",axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15)) +
  scale_x_discrete(labels=c("July 2014","August 2014", "October 2014","July 2015", "October 2015"))
ggsave(Actino_RA.my1,filename = "figures/Actinobacteria_RA_Month.Yr_6.11.21.pdf", width=12, height=8, dpi=600)

Cyano_RA.my1<-ggplot(bac_phy_RA_m.y.m[which(bac_phy_RA_m.y.m$Phylum=='Cyanobacteria'),], aes(x=Month.Yr, y=Rel.Ab, fill=Phylum))+geom_bar(stat="identity",colour="black")+
  theme_classic()+labs(title = "Cyanobacteria Relative Abundance by Month & Year", x="Month & Year", y="Relative Abundance", fill="Phylum")+scale_fill_manual(values = brightness("olivedrab4", 0.7))+
  theme(legend.position="none",axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15)) +
  scale_x_discrete(labels=c("July 2014","August 2014", "October 2014","July 2015","October 2015"))
ggsave(Cyano_RA.my1,filename = "figures/Cyanobacteria_RA_Month.Yr_6.11.21.pdf", width=12, height=8, dpi=600)

Acido_RA.my1<-ggplot(bac_phy_RA_m.y.m[which(bac_phy_RA_m.y.m$Phylum=='Acidobacteria'),], aes(x=Month.Yr, y=Rel.Ab, fill=Phylum))+geom_bar(stat="identity",colour="black")+
  theme_classic()+labs(title = "Acidobacteria Relative Abundance by Month & Year", x="Month & Year", y="Relative Abundance", fill="Phylum")+scale_fill_manual(values = saturation("tomato2", 0.6))+
  theme(legend.position="none",axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15)) +
  scale_x_discrete(labels=c("July 2014","August 2014", "October 2014","July 2015","October 2015"))
ggsave(Acido_RA.my1,filename = "figures/Acidobacteria_RA_Month.Yr_6.11.21.pdf", width=12, height=8, dpi=600)


#### 16S Specific Taxa by Elevation + Month + Year ####

## 2014 
head(bac_phy_RA_em2014)
bac_phy_RA_em2014.m<-melt(bac_phy_RA_em2014, id.vars="EvMo")
head(bac_phy_RA_em2014.m)
names(bac_phy_RA_em2014.m)[which(names(bac_phy_RA_em2014.m) == "variable")] <- "Phylum"
names(bac_phy_RA_em2014.m)[which(names(bac_phy_RA_em2014.m) == "value")] <- "Rel.Ab"
bac_phy_RA_em2014.m$EvMo<-factor(bac_phy_RA_em2014.m$EvMo, levels = c("July.400","August.400", "October.400","July.1100","August.1100", "October.1100", 
                                                                      "July.2000","August.2000","October.2000","July.2700", "August.2700", "October.2700"))

bac_phy_RA.em.2014.1<-ggplot(bac_phy_RA_em2014.m, aes(x=EvMo, y=Rel.Ab, fill=Phylum))+geom_bar(stat="identity",colour="black")+scale_x_discrete(labels=c("July-400ft","Aug-400ft", "Oct-400ft","July-1100ft","Aug-1100ft", "Oct-1100ft", "July-2000ft", "Aug-2000ft", "Oct-2000ft", "July-2700ft", "Aug-2700ft", "Oct-2700ft"))+
  theme_classic()+labs(title = "Microbial Phylum Relative Abundance by Month & Elevation (2014)", x="Month & Elevation (ft)", y="Relative Abundance", fill="Phylum")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15))
ggsave(bac_phy_RA.em.2014.1,filename = "figures/microb_phyla_RA_elev.mon.2014_6.11.21.pdf", width=20, height=10, dpi=600)

Proteo_RA.em.2014.1<-ggplot(bac_phy_RA_em2014.m[which(bac_phy_RA_em2014.m$Phylum=='Proteobacteria'),], aes(x=EvMo, y=Rel.Ab, fill=Phylum))+geom_bar(stat="identity",colour="black")+
  theme_classic()+labs(title = "Proteobacteria Relative Abundance by Month & Elevation (2014)", x="Month & Elevation (ft)", y="Relative Abundance", fill="Phylum")+scale_fill_manual(values = "cornflowerblue")+
  theme(legend.position="none",axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15)) +
  scale_x_discrete(labels=c("July-400ft","Aug-400ft", "Oct-400ft","July-1100ft","Aug-1100ft", "Oct-1100ft", "July-2000ft", "Aug-2000ft", "Oct-2000ft", "July-2700ft", "Aug-2700ft", "Oct-2700ft"))
ggsave(Proteo_RA.em.2014.1,filename = "figures/Proteobacteria_RA_Elev.Mon.2014_6.11.21.pdf", width=12, height=8, dpi=600)

Firm_RA.em.2014.1<-ggplot(bac_phy_RA_em2014.m[which(bac_phy_RA_em2014.m$Phylum=='Firmicutes'),], aes(x=EvMo, y=Rel.Ab, fill=Phylum))+geom_bar(stat="identity",colour="black")+
  theme_classic()+labs(title = "Firmicutes Relative Abundance by Month & Elevation (2014)", x="Month & Elevation (ft)", y="Relative Abundance", fill="Phylum")+scale_fill_manual(values = brightness("seagreen3", 0.8))+
  theme(legend.position="none",axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15)) +
  scale_x_discrete(labels=c("July-400ft","Aug-400ft", "Oct-400ft","July-1100ft","Aug-1100ft", "Oct-1100ft", "July-2000ft", "Aug-2000ft", "Oct-2000ft", "July-2700ft", "Aug-2700ft", "Oct-2700ft"))
ggsave(Firm_RA.em.2014.1,filename = "figures/Firmicutes_RA_Elev.Mon.2014_6.11.21.pdf", width=12, height=8, dpi=600)

Actino_RA.em.2014.1<-ggplot(bac_phy_RA_em2014.m[which(bac_phy_RA_em2014.m$Phylum=='Actinobacteria'),], aes(x=EvMo, y=Rel.Ab, fill=Phylum))+geom_bar(stat="identity",colour="black")+
  theme_classic()+labs(title = "Actinobacteria Relative Abundance by Month & Elevation (2014)", x="Month & Elevation (ft)", y="Relative Abundance", fill="Phylum")+scale_fill_manual(values = brightness("coral2", 1))+
  theme(legend.position="none",axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15)) +
  scale_x_discrete(labels=c("July-400ft","Aug-400ft", "Oct-400ft","July-1100ft","Aug-1100ft", "Oct-1100ft", "July-2000ft", "Aug-2000ft", "Oct-2000ft", "July-2700ft", "Aug-2700ft", "Oct-2700ft"))
ggsave(Actino_RA.em.2014.1,filename = "figures/Actinobacteria_RA_Elev.Mon.2014_6.11.21.pdf", width=12, height=8, dpi=600)

Cyano_RA.em.2014.1<-ggplot(bac_phy_RA_em2014.m[which(bac_phy_RA_em2014.m$Phylum=='Cyanobacteria'),], aes(x=EvMo, y=Rel.Ab, fill=Phylum))+geom_bar(stat="identity",colour="black")+
  theme_classic()+labs(title = "Cyanobacteria Relative Abundance by Month & Elevation (2014)", x="Month & Elevation (ft)", y="Relative Abundance", fill="Phylum")+scale_fill_manual(values = brightness("olivedrab4", 0.7))+
  theme(legend.position="none",axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15)) +
  scale_x_discrete(labels=c("July-400ft","Aug-400ft", "Oct-400ft","July-1100ft","Aug-1100ft", "Oct-1100ft", "July-2000ft", "Aug-2000ft", "Oct-2000ft", "July-2700ft", "Aug-2700ft", "Oct-2700ft"))
ggsave(Cyano_RA.em.2014.1,filename = "figures/Cyanobacteria_RA_Elev.Mon.2014_6.11.21.pdf", width=12, height=8, dpi=600)

Acido_RA.em.2014.1<-ggplot(bac_phy_RA_em2014.m[which(bac_phy_RA_em2014.m$Phylum=='Acidobacteria'),], aes(x=EvMo, y=Rel.Ab, fill=Phylum))+geom_bar(stat="identity",colour="black")+
  theme_classic()+labs(title = "Acidobacteria Relative Abundance by Month & Elevation (2014)", x="Month & Elevation (ft)", y="Relative Abundance", fill="Phylum")+scale_fill_manual(values = saturation("tomato2", 0.6))+
  theme(legend.position="none",axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15)) +
  scale_x_discrete(labels=c("July-400ft","Aug-400ft", "Oct-400ft","July-1100ft","Aug-1100ft", "Oct-1100ft", "July-2000ft", "Aug-2000ft", "Oct-2000ft", "July-2700ft", "Aug-2700ft", "Oct-2700ft"))
ggsave(Acido_RA.em.2014.1,filename = "figures/Acidobacteria_RA_Elev.Mon.2014_6.11.21.pdf", width=12, height=8, dpi=600)

## 2015 
head(bac_phy_RA_em2015)
bac_phy_RA_em2015.m<-melt(bac_phy_RA_em2015, id.vars="EvMo")
head(bac_phy_RA_em2015.m)
names(bac_phy_RA_em2015.m)[which(names(bac_phy_RA_em2015.m) == "variable")] <- "Phylum"
names(bac_phy_RA_em2015.m)[which(names(bac_phy_RA_em2015.m) == "value")] <- "Rel.Ab"
bac_phy_RA_em2015.m$EvMo<-factor(bac_phy_RA_em2015.m$EvMo, levels = c("July.400","October.400","July.1100", "October.1100", 
                                                                      "July.2000","October.2000","July.2700", "October.2700"))

bac_phy_RA.em.2015.1<-ggplot(bac_phy_RA_em2015.m, aes(x=EvMo, y=Rel.Ab, fill=Phylum))+geom_bar(stat="identity",colour="black")+scale_x_discrete(labels=c("July-400ft","Oct-400ft","July-1100ft","Oct-1100ft", "July-2000ft", "Oct-2000ft", "July-2700ft", "Oct-2700ft"))+
  theme_classic()+labs(title = "Microbial Phylum Relative Abundance by Month & Elevation (2015)", x="Month & Elevation (ft)", y="Relative Abundance", fill="Phylum")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15))
ggsave(bac_phy_RA.em.2015.1,filename = "figures/microb_phyla_RA_elev.mon.2015_6.11.21.pdf", width=20, height=10, dpi=600)

Proteo_RA.em.2015.1<-ggplot(bac_phy_RA_em2015.m[which(bac_phy_RA_em2015.m$Phylum=='Proteobacteria'),], aes(x=EvMo, y=Rel.Ab, fill=Phylum))+geom_bar(stat="identity",colour="black")+
  theme_classic()+labs(title = "Proteobacteria Relative Abundance by Month & Elevation (2015)", x="Month & Elevation (ft)", y="Relative Abundance", fill="Phylum")+scale_fill_manual(values = "cornflowerblue")+
  theme(legend.position="none",axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15)) +
  scale_x_discrete(labels=c("July-400ft","Oct-400ft","July-1100ft","Oct-1100ft", "July-2000ft", "Oct-2000ft", "July-2700ft", "Oct-2700ft"))
ggsave(Proteo_RA.em.2015.1,filename = "figures/Proteobacteria_RA_Elev.Mon.2015_6.11.21.pdf", width=12, height=8, dpi=600)

Firm_RA.em.2015.1<-ggplot(bac_phy_RA_em2015.m[which(bac_phy_RA_em2015.m$Phylum=='Firmicutes'),], aes(x=EvMo, y=Rel.Ab, fill=Phylum))+geom_bar(stat="identity",colour="black")+
  theme_classic()+labs(title = "Firmicutes Relative Abundance by Month & Elevation (2015)", x="Month & Elevation (ft)", y="Relative Abundance", fill="Phylum")+scale_fill_manual(values = brightness("seagreen3", 0.8))+
  theme(legend.position="none",axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15)) +
  scale_x_discrete(labels=c("July-400ft","Oct-400ft","July-1100ft","Oct-1100ft", "July-2000ft", "Oct-2000ft", "July-2700ft", "Oct-2700ft"))
ggsave(Firm_RA.em.2015.1,filename = "figures/Firmicutes_RA_Elev.Mon.2015_6.11.21.pdf", width=12, height=8, dpi=600)

Actino_RA.em.2015.1<-ggplot(bac_phy_RA_em2015.m[which(bac_phy_RA_em2015.m$Phylum=='Actinobacteria'),], aes(x=EvMo, y=Rel.Ab, fill=Phylum))+geom_bar(stat="identity",colour="black")+
  theme_classic()+labs(title = "Actinobacteria Relative Abundance by Month & Elevation (2015)", x="Month & Elevation (ft)", y="Relative Abundance", fill="Phylum")+scale_fill_manual(values = brightness("coral2", 1))+
  theme(legend.position="none",axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15)) +
  scale_x_discrete(labels=c("July-400ft","Oct-400ft","July-1100ft","Oct-1100ft", "July-2000ft", "Oct-2000ft", "July-2700ft", "Oct-2700ft"))
ggsave(Actino_RA.em.2015.1,filename = "figures/Actinobacteria_RA_Elev.Mon.2015_6.11.21.pdf", width=12, height=8, dpi=600)

Cyano_RA.em.2015.1<-ggplot(bac_phy_RA_em2015.m[which(bac_phy_RA_em2015.m$Phylum=='Cyanobacteria'),], aes(x=EvMo, y=Rel.Ab, fill=Phylum))+geom_bar(stat="identity",colour="black")+
  theme_classic()+labs(title = "Cyanobacteria Relative Abundance by Month & Elevation (2015)", x="Month & Elevation (ft)", y="Relative Abundance", fill="Phylum")+scale_fill_manual(values = brightness("olivedrab4", 0.7))+
  theme(legend.position="none",axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15)) +
  scale_x_discrete(labels=c("July-400ft", "Oct-400ft","July-1100ft", "Oct-1100ft", "July-2000ft", "Oct-2000ft", "July-2700ft",  "Oct-2700ft"))
ggsave(Cyano_RA.em.2015.1,filename = "figures/Cyanobacteria_RA_Elev.Mon.2015_6.11.21.pdf", width=12, height=8, dpi=600)

Acido_RA.em.2015.1<-ggplot(bac_phy_RA_em2015.m[which(bac_phy_RA_em2015.m$Phylum=='Acidobacteria'),], aes(x=EvMo, y=Rel.Ab, fill=Phylum))+geom_bar(stat="identity",colour="black")+
  theme_classic()+labs(title = "Acidobacteria Relative Abundance by Month & Elevation (2015)", x="Month & Elevation (ft)", y="Relative Abundance", fill="Phylum")+scale_fill_manual(values = saturation("tomato2", 0.6))+
  theme(legend.position="none",axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15)) +
  scale_x_discrete(labels=c("July-400ft", "Oct-400ft","July-1100ft", "Oct-1100ft", "July-2000ft", "Oct-2000ft", "July-2700ft",  "Oct-2700ft"))
ggsave(Acido_RA.em.2015.1,filename = "figures/Acidobacteria_RA_Elev.Mon.2015_6.11.21.pdf", width=12, height=8, dpi=600)

#### Combine 2014 + 2015 Figs together

Prot_comp<-ggarrange(Proteo_RA.em.2014.1, Proteo_RA.em.2015.1, legend="none", ncol=2, nrow=1)
ggsave(Prot_comp,filename = "figures/Proteobacteria_ElevMonth_both_yrs_6.11.21.pdf", width=25, height=10, dpi=600)

Firm_comp<-ggarrange(Firm_RA.em.2014.1, Firm_RA.em.2015.1, legend="none", ncol=2, nrow=1)
ggsave(Prot_comp,filename = "figures/Firmicutes_ElevMonth_both_yrs_6.11.21.pdf", width=25, height=10, dpi=600)

Actino_comp<-ggarrange(Actino_RA.em.2014.1, Actino_RA.em.2015.1, legend="none", ncol=2, nrow=1)
ggsave(Actino_comp,filename = "figures/Actinobacteria_ElevMonth_both_yrs_6.11.21.pdf", width=25, height=10, dpi=600)

Cyano_comp<-ggarrange(Cyano_RA.em.2014.1, Cyano_RA.em.2015.1, legend="none", ncol=2, nrow=1)
ggsave(Cyano_comp,filename = "figures/Cyanobacteria_ElevMonth_both_yrs_6.11.21.pdf", width=25, height=10, dpi=600)

Acido_comp<-ggarrange(Acido_RA.em.2014.1, Acido_RA.em.2015.1, legend="none", ncol=2, nrow=1)
ggsave(Acido_comp,filename = "figures/Acidobacteria_ElevMonth_both_yrs_6.11.21.pdf", width=25, height=10, dpi=600)



#### ITS1 Specific Taxa by Binned Month ####

## Phylum First
head(its1_RA_phy.mbin)
its1_RA_phy.mbin$MonthBin<-rownames(its1_RA_phy.mbin)
its1.RA.phy.mbin.m<-melt(its1_RA_phy.mbin)
head(its1.RA.phy.mbin.m)
names(its1.RA.phy.mbin.m)[which(names(its1.RA.phy.mbin.m) == "variable")] <- "Phylum"
names(its1.RA.phy.mbin.m)[which(names(its1.RA.phy.mbin.m) == "value")] <- "Rel.Ab"
its1.RA.phy.mbin.m$MonthBin<-factor(its1.RA.phy.mbin.m$MonthBin, levels = c("Early","Late"))
its1.RA.mbin.all<-merge(its1.RA.phy.mbin.m, metadata, by = "MonthBin")

its1_phy_RA.MB<-ggplot(its1.RA.phy.mbin.m, aes(x=factor(MonthBin), y=Rel.Ab, fill=Phylum))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  labs(title = "Fungal Phylum Relative Abundance by Time in Dry Season", x="Time in Dry Season", y="Relative Abundance", fill="Phylum")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15))
ggsave(its1_phy_RA.MB,filename = "figures/fungal_phyla_RA_monthbin_8.9.21.pdf", width=10, height=6, dpi=600)

Basidio_RA.mbin<-ggplot(its1.RA.phy.mbin.m[which(its1.RA.phy.mbin.m$Phylum=='Basidiomycota'),], aes(x=MonthBin, y=Rel.Ab, fill=Phylum))+geom_bar(stat="identity",colour="black")+
  theme_classic()+labs(title = "Proteobacteria Relative Abundance by Time in Dry Season", x="Time in Dry Season", y="Relative Abundance", fill="Phylum")+scale_fill_manual(values = brightness("orange2", 0.83))+
  theme(legend.position="none",axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15))
ggsave(Basidio_RA.mbin,filename = "figures/Basidiomycota_RA_monthbin_8.9.21.pdf", width=10, height=8, dpi=600)

## Class next
head(its1_RA_cls.mbin)
its1_RA_cls.mbin$MonthBin<-rownames(its1_RA_cls.mbin)
its1.RA.cls.mbin.m<-melt(its1_RA_cls.mbin)
head(its1.RA.cls.mbin.m)
names(its1.RA.cls.mbin.m)[which(names(its1.RA.cls.mbin.m) == "variable")] <- "Class"
names(its1.RA.cls.mbin.m)[which(names(its1.RA.cls.mbin.m) == "value")] <- "Rel.Ab"
its1.RA.cls.mbin.m$MonthBin<-factor(its1.RA.cls.mbin.m$MonthBin, levels = c("Early","Late"))

its1_cls_RA.mbin<-ggplot(its1.RA.cls.mbin.m, aes(x=factor(MonthBin), y=Rel.Ab, fill=Class))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  labs(title = "Fungal Class Relative Abundance by Time in Dry Season", x="Time in Dry Season", y="Relative Abundance", fill="Class")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15))
ggsave(its1_cls_RA.mbin,filename = "figures/fungal_class_RA_monthbin_7.6.21.pdf", width=12, height=6, dpi=600)

Trem_RA.mbiny<-ggplot(its1.RA.cls.mbin.m[which(its1.RA.cls.mbin.m$Class=='Tremellomycetes'),], aes(x=factor(MonthBin), y=Rel.Ab, fill=Class))+geom_bar(stat="identity",colour="black")+
  theme_classic()+labs(title = "Tremellomycetes Relative Abundance by Time in Dry Season", x="Time in Dry Season", y="Relative Abundance", fill="Class")+scale_fill_manual(values = brightness("violet", 0.91))+
  theme(legend.position="none",axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15)) +
  stat_compare_means(method="kruskal.test", hide.ns = TRUE,label = "p.signif")
ggsave(Trem_RA.mbiny,filename = "figures/Tremellomycetes_RA_monthbin_8.9.21.pdf", width=12, height=8, dpi=600)

Dothideo_RA.mbiny<-ggplot(its1.RA.cls.mbin.m[which(its1.RA.cls.mbin.m$Class=='Dothideomycetes'),], aes(x=MonthBin, y=Rel.Ab, fill=Class))+geom_bar(stat="identity",colour="black")+
  theme_classic()+labs(title = "Dothideomycetes Relative Abundance by Time in Dry Season", x="Time in Dry Season", y="Relative Abundance", fill="Class")+scale_fill_manual(values = brightness("olivedrab", 0.65))+
  theme(legend.position="none",axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15)) + 
  stat_compare_means(method="kruskal.test", hide.ns = TRUE,label = "p.signif")
ggsave(Dothideo_RA.mbiny,filename = "figures/Dothideomycetes_RA_monthbin_8.9.21.pdf", width=12, height=8, dpi=600)

Eurotio_RA.mbiny<-ggplot(its1.RA.cls.mbin.m[which(its1.RA.cls.mbin.m$Class=='Eurotiomycetes'),], aes(x=MonthBin, y=Rel.Ab, fill=Class))+geom_bar(stat="identity",colour="black")+
  theme_classic()+labs(title = "Eurotiomycetes Relative Abundance by Time in Dry Season", x="Time in Dry Season", y="Relative Abundance", fill="Class")+scale_fill_manual(values = brightness("green2", 0.7))+
  theme(legend.position="none",axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15)) + 
  stat_compare_means(method="kruskal.test", hide.ns = TRUE,label = "p.signif")
ggsave(Eurotio_RA.mbiny,filename = "figures/Eurotiomycetes_RA_monthbin_8.11.21.pdf", width=12, height=8, dpi=600)

Arthionio_RA.mbiny<-ggplot(its1.RA.cls.mbin.m[which(its1.RA.cls.mbin.m$Class=='Arthoniomycetes'),], aes(x=MonthBin, y=Rel.Ab, fill=Class))+geom_bar(stat="identity",colour="black")+
  theme_classic()+labs(title = "Arthoniomycetes Relative Abundance by Time in Dry Season", x="Time in Dry Season", y="Relative Abundance", fill="Class")+scale_fill_manual(values = brightness("orange3", 0.8))+
  theme(legend.position="none",axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15)) + 
  stat_compare_means(method="kruskal.test", hide.ns = TRUE,label = "p.signif")
ggsave(Arthionio_RA.mbiny,filename = "figures/Arthoniomycetes_RA_monthbin_8.11.21.pdf", width=12, height=8, dpi=600)

Cystobasidio_RA.mbiny<-ggplot(its1.RA.cls.mbin.m[which(its1.RA.cls.mbin.m$Class=='Cystobasidiomycetes'),], aes(x=MonthBin, y=Rel.Ab, fill=Class))+geom_bar(stat="identity",colour="black")+
  theme_classic()+labs(title = "Cystobasidiomycetes Relative Abundance by Time in Dry Season", x="Time in Dry Season", y="Relative Abundance", fill="Class")+scale_fill_manual(values = brightness("peru", 0.75))+
  theme(legend.position="none",axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15)) + 
  stat_compare_means(method="kruskal.test", hide.ns = TRUE,label = "p.signif")
ggsave(Cystobasidio_RA.mbiny,filename = "figures/Cystobasidiomycetes_RA_monthbin_8.11.21.pdf", width=12, height=8, dpi=600)


#### 16S Specific Taxa by Binned Month ####

## Phylum First
head(bac_RA_phy.mbin)
bac_RA_phy.mbin$MonthBin<-rownames(bac_RA_phy.mbin)
bac.RA.phy.mbin.m<-melt(bac_RA_phy.mbin)
head(bac.RA.phy.mbin.m)
names(bac.RA.phy.mbin.m)[which(names(bac.RA.phy.mbin.m) == "variable")] <- "Phylum"
names(bac.RA.phy.mbin.m)[which(names(bac.RA.phy.mbin.m) == "value")] <- "Rel.Ab"
bac.RA.phy.mbin.m$MonthBin<-factor(bac.RA.phy.mbin.m$MonthBin, levels = c("Early","Late"))
bac.RA.mbin.all<-merge(bac.RA.phy.mbin.m, metadata, by = "MonthBin")

bac_phy_RA.MB<-ggplot(bac.RA.phy.mbin.m, aes(x=factor(MonthBin), y=Rel.Ab, fill=Phylum))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  labs(title = "Microbial Phylum Relative Abundance by Time in Dry Season", x="Time in Dry Season", y="Relative Abundance", fill="Phylum")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15))
# stat_compare_means(method="kruskal.test")
ggsave(bac_phy_RA.MB,filename = "figures/microb_phyla_RA_mbin_8.9.21.pdf", width=8, height=6, dpi=600)

Proteo_RA.MB<-ggplot(bac.RA.phy.mbin.m[which(bac.RA.phy.mbin.m$Phylum=='Proteobacteria'),], aes(x=factor(MonthBin), y=Rel.Ab, fill=Phylum))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+
  theme_classic()+labs(title = "Proteobacteria Relative Abundance by Time in Dry Season", x="Time in Dry Season", y="Relative Abundance", fill="Phylum")+scale_fill_manual(values = "cornflowerblue")+
  theme(legend.position="none",axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15))
ggsave(Proteo_RA.MB,filename = "figures/Proteobacteria_RA_BinnedMonth_8.9.21.pdf", width=8, height=6, dpi=600)

Firm_RA.MB<-ggplot(bac.RA.phy.mbin.m[which(bac.RA.phy.mbin.m$Phylum=='Firmicutes'),], aes(x=factor(MonthBin), y=Rel.Ab, fill=Phylum))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+
  theme_classic()+labs(title = "Firmicutes Relative Abundance by Time in Dry Season", x="Time in Dry Season", y="Relative Abundance", fill="Phylum")+scale_fill_manual(values = brightness("seagreen3", 0.8))+
  theme(legend.position="none",axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15))
ggsave(Firm_RA.MB,filename = "figures/Firmicutes_RA_BinnedMonth_8.9.21.pdf", width=8, height=6, dpi=600)

Actino_RA.MB<-ggplot(bac.RA.phy.mbin.m[which(bac.RA.phy.mbin.m$Phylum=='Actinobacteria'),], aes(x=factor(MonthBin), y=Rel.Ab, fill=Phylum))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+
  theme_classic()+labs(title = "Actinobacteria Relative Abundance by Time in Dry Season", x="Time in Dry Season", y="Relative Abundance", fill="Phylum")+scale_fill_manual(values = brightness("coral2", 1))+
  theme(legend.position="none",axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15))
ggsave(Actino_RA.MB,filename = "figures/Actinobacteria_RA_BinnedMonth_8.9.21.pdf", width=8, height=6, dpi=600)

Cyano_RA.MB<-ggplot(bac.RA.phy.mbin.m[which(bac.RA.phy.mbin.m$Phylum=='Cyanobacteria'),], aes(x=factor(MonthBin), y=Rel.Ab, fill=Phylum))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+
  theme_classic()+labs(title = "Cyanobacteria Relative Abundance by Time in Dry Season", x="Time in Dry Season", y="Relative Abundance", fill="Phylum")+scale_fill_manual(values = brightness("olivedrab4", 0.7))+
  theme(legend.position="none",axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15))
ggsave(Cyano_RA.MB,filename = "figures/Cyanobacteria_RA_BinnedMonth_8.9.21.pdf", width=8, height=6, dpi=600)

Acido_RA.MB<-ggplot(bac.RA.phy.mbin.m[which(bac.RA.phy.mbin.m$Phylum=='Acidobacteria'),], aes(x=factor(MonthBin), y=Rel.Ab, fill=Phylum))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+
  theme_classic()+labs(title = "Acidobacteria Relative Abundance by Time in Dry Season", x="Time in Dry Season", y="Relative Abundance", fill="Phylum")+scale_fill_manual(values = saturation("tomato2", 0.6))+
  theme(legend.position="none",axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15))
ggsave(Acido_RA.MB,filename = "figures/Acidobacteria_RA_BinnedMonth_8.9.21.pdf", width=8, height=6, dpi=600)

## Class next
head(bac_RA_cls.mbin)
bac_RA_cls.mbin$MonthBin<-rownames(bac_RA_cls.mbin)
bac.RA.cls.mbin.m<-melt(bac_RA_cls.mbin)
head(bac.RA.cls.mbin.m)
names(bac.RA.cls.mbin.m)[which(names(bac.RA.cls.mbin.m) == "variable")] <- "Class"
names(bac.RA.cls.mbin.m)[which(names(bac.RA.cls.mbin.m) == "value")] <- "Rel.Ab"
bac.RA.cls.mbin.m$MonthBin<-factor(bac.RA.cls.mbin.m$MonthBin, levels = c("Early","Late"))

bac_cls_RA.MB<-ggplot(bac.RA.cls.mbin.m, aes(x=factor(MonthBin), y=Rel.Ab, fill=Class))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  labs(title = "Microbial Class Relative Abundance by Time in Dry Season", x="Time in Dry Season", y="Relative Abundance", fill="Class")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15))
#ggsave(bac_cls_RA.MB,filename = "figures/microb_class_RA_monthbin_7.6.21.pdf", width=8, height=6, dpi=600)


#### 16S (Specific) Taxa by Binned Month + Year ####

## Phylum 
head(bac_RA_phy_mbin.yr)
bac_RA_phy_mbin.yr$BinYr<-rownames(bac_RA_phy_mbin.yr)
bac_RA_phy_mbin.yr.m<-melt(bac_RA_phy_mbin.yr)
head(bac_RA_phy_mbin.yr.m)
names(bac_RA_phy_mbin.yr.m)[which(names(bac_RA_phy_mbin.yr.m) == "variable")] <- "Phylum"
names(bac_RA_phy_mbin.yr.m)[which(names(bac_RA_phy_mbin.yr.m) == "value")] <- "Rel.Ab"
bac_RA_phy_mbin.yr.m$BinYr<-factor(bac_RA_phy_mbin.yr.m$BinYr, levels = c("2014.Early","2014.Late", "2015.Early","2015.Late"))
head(bac_RA_phy_mbin.yr.m)

bac_phy_RA.mbiny<-ggplot(bac_RA_phy_mbin.yr.m, aes(x=BinYr, y=Rel.Ab, fill=Phylum))+geom_bar(stat="identity",colour="black")+scale_x_discrete(labels=c("Early 2014","Late 2014", "Early 2015","Late 2015"))+
  theme_classic()+labs(title = "Microbial Phylum Relative Abundance by Time in Dry Season & Year", x="Time in Season & Year", y="Relative Abundance", fill="Phylum")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15))
ggsave(bac_phy_RA.mbiny,filename = "figures/microb_phyla_RA_bin.month.year_8.9.21.pdf", width=15, height=8, dpi=600)

Proteo_RA.mbiny<-ggplot(bac_RA_phy_mbin.yr.m[which(bac_RA_phy_mbin.yr.m$Phylum=='Proteobacteria'),], aes(x=BinYr, y=Rel.Ab, fill=Phylum))+geom_bar(stat="identity",colour="black")+
  theme_classic()+labs(title = "Proteobacteria Relative Abundance by Time in Dry Season & Year", x="Time in Dry Season & Year", y="Relative Abundance", fill="Phylum")+scale_fill_manual(values = "cornflowerblue")+
  theme(legend.position="none",axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15)) +
  scale_x_discrete(labels=c("Early 2014","Late 2014", "Early 2015","Late 2015"))
ggsave(Proteo_RA.mbiny,filename = "figures/Proteobacteria_RA_bin.month.yr_8.9.21.pdf", width=12, height=8, dpi=600)

Firm_RA.mbiny<-ggplot(bac_RA_phy_mbin.yr.m[which(bac_RA_phy_mbin.yr.m$Phylum=='Firmicutes'),], aes(x=BinYr, y=Rel.Ab, fill=Phylum))+geom_bar(stat="identity",colour="black")+
  theme_classic()+labs(title = "Firmicutes Relative Abundance by Time in Dry Season & Year", x="Time in Dry Season & Year", y="Relative Abundance", fill="Phylum")+scale_fill_manual(values = brightness("seagreen3", 0.8))+
  theme(legend.position="none",axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15)) +
  scale_x_discrete(labels=c("Early 2014","Late 2014", "Early 2015","Late 2015"))
ggsave(Firm_RA.mbiny,filename = "figures/Firmicutes_RA_bin.month.yr_8.9.21.pdf", width=12, height=8, dpi=600)

Actino_RA.mbiny<-ggplot(bac_RA_phy_mbin.yr.m[which(bac_RA_phy_mbin.yr.m$Phylum=='Actinobacteria'),], aes(x=BinYr, y=Rel.Ab, fill=Phylum))+geom_bar(stat="identity",colour="black")+
  theme_classic()+labs(title = "Actinobacteria Relative Abundance by Time in Dry Season & Year", x="Time in Dry Season & Year", y="Relative Abundance", fill="Phylum")+scale_fill_manual(values = brightness("coral2", 1))+
  theme(legend.position="none",axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15)) +
  scale_x_discrete(labels=c("Early 2014","Late 2014", "Early 2015","Late 2015"))
ggsave(Actino_RA.mbiny,filename = "figures/Actinobacteria_RA_bin.month.yr_8.9.21.pdf", width=12, height=8, dpi=600)

Cyano_RA.mbiny<-ggplot(bac_RA_phy_mbin.yr.m[which(bac_RA_phy_mbin.yr.m$Phylum=='Cyanobacteria'),], aes(x=BinYr, y=Rel.Ab, fill=Phylum))+geom_bar(stat="identity",colour="black")+
  theme_classic()+labs(title = "Cyanobacteria Relative Abundance by Time in Dry Season & Year", x="Time in Dry Season & Year", y="Relative Abundance", fill="Phylum")+scale_fill_manual(values = brightness("olivedrab4", 0.7))+
  theme(legend.position="none",axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15)) +
  scale_x_discrete(labels=c("Early 2014","Late 2014", "Early 2015","Late 2015"))
ggsave(Cyano_RA.mbiny,filename = "figures/Cyanobacteria_RA_bin.month.yr_8.9.21.pdf", width=12, height=8, dpi=600)

Acido_RA.mbiny<-ggplot(bac_RA_phy_mbin.yr.m[which(bac_RA_phy_mbin.yr.m$Phylum=='Acidobacteria'),], aes(x=BinYr, y=Rel.Ab, fill=Phylum))+geom_bar(stat="identity",colour="black")+
  theme_classic()+labs(title = "Acidobacteria Relative Abundance by Time in Dry Season & Year", x="Time in Dry Season & Year", y="Relative Abundance", fill="Phylum")+scale_fill_manual(values = saturation("tomato2", 0.6))+
  theme(legend.position="none",axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15)) +
  scale_x_discrete(labels=c("Early 2014","Late 2014", "Early 2015","Late 2015"))
ggsave(Acido_RA.mbiny,filename = "figures/Acidobacteria_RA_bin.month.yr_8.9.21.pdf", width=12, height=8, dpi=600)



#### ITS1 (Specific) Taxa by Binned Month + Year ####

## Phylum 
head(its1_RA_phy_mbin.yr)
its1_RA_phy_mbin.yr$BinYr<-rownames(its1_RA_phy_mbin.yr)
its1_RA_phy_mbin.yr.m<-melt(its1_RA_phy_mbin.yr)
head(its1_RA_phy_mbin.yr.m)
names(its1_RA_phy_mbin.yr.m)[which(names(its1_RA_phy_mbin.yr.m) == "variable")] <- "Phylum"
names(its1_RA_phy_mbin.yr.m)[which(names(its1_RA_phy_mbin.yr.m) == "value")] <- "Rel.Ab"
its1_RA_phy_mbin.yr.m$BinYr<-factor(its1_RA_phy_mbin.yr.m$BinYr, levels = c("2014.Early","2014.Late", "2015.Early","2015.Late"))

its1_phy_RA.mbiny<-ggplot(its1_RA_phy_mbin.yr.m, aes(x=BinYr, y=Rel.Ab, fill=Phylum))+geom_bar(stat="identity",colour="black")+scale_x_discrete(labels=c("Early 2014","Late 2014", "Early 2015","Late 2015"))+
  theme_classic()+labs(title = "Fungal Phylum Relative Abundance by Time in Dry Season & Year", x="Time in Season & Year", y="Relative Abundance", fill="Phylum")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15))
ggsave(its1_phy_RA.mbiny,filename = "figures/fungal_phyla_RA_bin.month.year_8.9.21.pdf", width=15, height=8, dpi=600)

Basidio_RA.mbiny<-ggplot(its1_RA_phy_mbin.yr.m[which(its1_RA_phy_mbin.yr.m$Phylum=='Basidiomycota'),], aes(x=BinYr, y=Rel.Ab, fill=Phylum))+geom_bar(stat="identity",colour="black")+
  theme_classic()+labs(title = "Proteobacteria Relative Abundance by Time in Dry Season & Year", x="Time in Dry Season & Year", y="Relative Abundance", fill="Phylum")+scale_fill_manual(values = brightness("orange2", 0.83))+
  theme(legend.position="none",axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15)) +
  scale_x_discrete(labels=c("Early 2014","Late 2014", "Early 2015","Late 2015"))
ggsave(Basidio_RA.mbiny,filename = "figures/Basidiomycota_RA_bin.month.yr_8.9.21.pdf", width=12, height=8, dpi=600)

# class

head(its1_RA_cls_mbin.yr)
its1_RA_cls_mbin.yr$BinYr<-rownames(its1_RA_cls_mbin.yr)
its1_RA_cls_mbin.yr.m<-melt(its1_RA_cls_mbin.yr)
head(its1_RA_cls_mbin.yr.m)
names(its1_RA_cls_mbin.yr.m)[which(names(its1_RA_cls_mbin.yr.m) == "variable")] <- "Class"
names(its1_RA_cls_mbin.yr.m)[which(names(its1_RA_cls_mbin.yr.m) == "value")] <- "Rel.Ab"
its1_RA_cls_mbin.yr.m$BinYr<-factor(its1_RA_cls_mbin.yr.m$BinYr, levels = c("2014.Early","2014.Late", "2015.Early","2015.Late"))
head(its1_RA_cls_mbin.yr.m)

its1_cls_RA.mbiny<-ggplot(its1_RA_cls_mbin.yr.m, aes(x=BinYr, y=Rel.Ab, fill=Class))+geom_bar(stat="identity",colour="black")+scale_x_discrete(labels=c("Early 2014","Late 2014", "Early 2015","Late 2015"))+
  theme_classic()+labs(title = "Fungal Class Relative Abundance by Time in Dry Season & Year", x="Time in Season & Year", y="Relative Abundance", fill="Class")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15))
ggsave(its1_cls_RA.mbiny,filename = "figures/fungal_cls_RA_bin.month.year_8.9.21.pdf", width=15, height=8, dpi=600)

# Check classes w/ greater than 10% relative abundance
its1_cls_mbin.yr.10<-subset(its1_RA_cls_mbin.yr.m, its1_RA_cls_mbin.yr.m$Rel.Ab>0.1)

its1_cls_RA.mbiny2<-ggplot(its1_cls_mbin.yr.10, aes(x=BinYr, y=Rel.Ab, fill=Class))+geom_bar(stat="identity",colour="black")+scale_x_discrete(labels=c("Early 2014","Late 2014", "Early 2015","Late 2015"))+
  theme_classic()+labs(title = "Fungal Class Relative Abundance by Time in Dry Season & Year", x="Time in Season & Year", y="Relative Abundance", fill="Class")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15))
ggsave(its1_cls_RA.mbiny,filename = "figures/fungal_cls_RA_bin.month.year_8.9.21.pdf", width=15, height=8, dpi=600)

Trem_RA.mbiny<-ggplot(its1_RA_cls_mbin.yr.m[which(its1_RA_cls_mbin.yr.m$Class=='Tremellomycetes'),], aes(x=BinYr, y=Rel.Ab, fill=Class))+geom_bar(stat="identity",colour="black")+
  theme_classic()+labs(title = "Tremellomycetes Relative Abundance by Time in Dry Season & Year", x="Time in Dry Season & Year", y="Relative Abundance", fill="Class")+scale_fill_manual(values = brightness("violet", 0.91))+
  theme(legend.position="none",axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15)) +
  scale_x_discrete(labels=c("Early 2014","Late 2014", "Early 2015","Late 2015"))+stat_compare_means(comparisons = list(c(1,2), c(2,3),  c(1,3),  c(3,4),  c(2,4),  c(1,4)), method="kruskal.test", hide.ns = TRUE,label = "p.signif")
ggsave(Trem_RA.mbiny,filename = "figures/Tremellomycetes_RA_bin.month.yr_8.9.21.pdf", width=12, height=8, dpi=600)

Dothideo_RA.mbiny<-ggplot(its1_RA_cls_mbin.yr.m[which(its1_RA_cls_mbin.yr.m$Class=='Dothideomycetes'),], aes(x=BinYr, y=Rel.Ab, fill=Class))+geom_bar(stat="identity",colour="black")+
  theme_classic()+labs(title = "Dothideomycetes Relative Abundance by Time in Dry Season & Year", x="Time in Dry Season & Year", y="Relative Abundance", fill="Class")+scale_fill_manual(values = brightness("olivedrab", 0.65))+
  theme(legend.position="none",axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15)) +
  scale_x_discrete(labels=c("Early 2014","Late 2014", "Early 2015","Late 2015"))
ggsave(Dothideo_RA.mbiny,filename = "figures/Dothideomycetes_RA_bin.month.yr_8.9.21.pdf", width=12, height=8, dpi=600)
