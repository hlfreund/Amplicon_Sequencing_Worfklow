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

season_col <- melt(c('Early'="#e0607e", 'Late'="#721121"))
season_col$MonthBin<-rownames(season_col)
names(season_col)[which(names(season_col) == "value")] <- "color2"
season_col

metadata<-merge(metadata, season_col, by="MonthBin")
head(metadata)
metadata$color2 <- as.character(metadata$color2)
rownames(metadata)<-metadata$SampleID
head(metadata)

## Reorder rows of metadata and all dfs
metadata=metadata[rownames(its1_otu_table),] ## reorder metadata to have same rows as original OTU table
# ^ this indexing method will only work if the two dfs have the same # of rows AND the same row names!
bac_otu_table=bac_otu_table[rownames(metadata),] ## reorder bacterial OTU table to have same rows as ITS1 OTU table + metadata

#meta_category<-subset(metadata, select=c(SampleID, Year, Month, Site, Elevation)) # subset only part of the metadata we need
#meta_quant<-subset(metadata, select=-c(Year, Month, Site, Elevation)) # subset only part of the metadata we need

#### Prep count data for taxa summaries ####
# bacteria first
bac_otu_counts$OTU_ID<-rownames(bac_otu_counts)
bac.tax<-merge(bac_otu_counts, bac_otu_tax_ID, by="OTU_ID")
bac.tax.m<-melt(bac.tax)
head(bac.tax.m)
colnames(bac.tax.m)[which(names(bac.tax.m) == "value")] <- "Count"
colnames(bac.tax.m)[which(names(bac.tax.m) == "variable")] <- "SampleID"

# fungi next
its1_otu_counts$OTU_ID<-rownames(its1_otu_counts)
its1.tax<-merge(its1_otu_counts, its1_otu_tax_ID, by="OTU_ID")
its1.tax.m<-melt(its1.tax)
head(its1.tax.m)
colnames(its1.tax.m)[which(names(its1.tax.m) == "value")] <- "Count"
colnames(its1.tax.m)[which(names(its1.tax.m) == "variable")] <- "SampleID"

#### 16S Taxonomic Summaries by Elevation ####
head(bac.tax.m)
#major_taxa_proportions_tab <- apply(major_taxa_counts_tab, 2, function(x) x/sum(x)*100)

# Phylum
bac_phy<- as.data.frame(dcast(bac.tax.m,SampleID~Phylum, value.var="Count", fun.aggregate=sum)) ### 
head(bac_phy) # counts by phyla + elevation
rownames(bac_phy)<-bac_phy$SampleID

bac.phy_RA<-data.frame(decostand(bac_phy[,-c(1)], method="total", MARGIN=1, na.rm=TRUE))  
# relative abundance of taxa data where everything is divided by margin total (default MARGIN = 1 = rows) -- rows = samples
rowSums(bac.phy_RA)
head(bac.phy_RA)
bac.phy_RA$SampleID<-rownames(bac.phy_RA)
head(bac.phy_RA)
bac.phy_m<-melt(bac.phy_RA)
head(bac.phy_m)

names(bac.phy_m)[which(names(bac.phy_m) == "variable")] <- "Phylum"
names(bac.phy_m)[which(names(bac.phy_m) == "value")] <- "Count"
head(bac.phy_m) ## relative abundance based on sum of counts by phyla!
b.phy.dat<-merge(bac.phy_m, metadata, by="SampleID")
head(b.phy.dat)
b.phy.dat$Elevation<-factor(b.phy.dat$Elevation, levels=c("2700", "2000", "1100", "400"))

b.phy_1<-subset(b.phy.dat, (Count)>(1/100)) ## DROP BACTERIAL PHYLA that are less than 1% abundant!!!!!!1\
head(b.phy_1)

# 16S phyla

b1<-ggplot(b.phy.dat, aes(Phylum, Count)) +
  geom_jitter(aes(color=factor(Elevation)), size=2, width=0.15, height=0) +
  scale_color_manual(name ="Elevation", labels=c("400"="400 ft", "1100"="1100 ft", "2000"="2000 ft", "2700"="2700 ft"),
                     values=unique(b.phy.dat$color[order(b.phy.dat$Elevation)])) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  labs(x="Microbial Phyla", y="Relative Abundance", title="Bacteria/Archaea & Elevation")
ggsave(b1,filename = "figures/16S_taxa.summary_phyla_elev_6.22.21.pdf", width=10, height=6, dpi=600)

b2<-ggplot(b.phy_1, aes(Phylum, Count)) +
  geom_jitter(aes(color=factor(Elevation)), size=2, width=0.15, height=0) +
  scale_color_manual(name ="Elevation", labels=c("400"="400 ft", "1100"="1100 ft", "2000"="2000 ft", "2700"="2700 ft"),
                     values=unique(b.phy_1$color[order(b.phy_1$Elevation)])) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  labs(x="Microbial Phyla", y="Relative Abundance", title="Bacteria/Archaea & Elevation", subtitle="Includes Phyla > 1% Relative Abundance")
ggsave(b2,filename = "figures/16S_taxa.summary_phyla.1perc_elev_6.22.21.pdf", width=10, height=6, dpi=600)

# Class
bac_cls<- as.data.frame(dcast(bac.tax.m,SampleID~Class, value.var="Count", fun.aggregate=sum)) ### 
head(bac_cls) # counts by clsla + elevation
rownames(bac_cls)<-bac_cls$SampleID

bac.cls_RA<-data.frame(decostand(bac_cls[,-c(1)], method="total", MARGIN=1, na.rm=TRUE))  
# relative abundance of taxa data where everything is divided by margin total (default MARGIN = 1 = rows) -- rows = samples
rowSums(bac.cls_RA)
head(bac.cls_RA)
bac.cls_RA$SampleID<-rownames(bac.cls_RA)
head(bac.cls_RA)
bac.cls_m<-melt(bac.cls_RA)
head(bac.cls_m)

names(bac.cls_m)[which(names(bac.cls_m) == "variable")] <- "Class"
names(bac.cls_m)[which(names(bac.cls_m) == "value")] <- "Count"
head(bac.cls_m) ## relative abundance based on sum of counts by clsla!
b.cls.dat<-merge(bac.cls_m, metadata, by="SampleID")
head(b.cls.dat)
b.cls.dat$Elevation<-factor(b.cls.dat$Elevation, levels=c("2700", "2000", "1100", "400"))

b.cls_1<-subset(b.cls.dat, (Count)>(1/100)) ## DROP BACTERIAL clsLA that are less than 1% abundant!!!!!!1\
head(b.cls_1)

# 16S Class

b3<-ggplot(b.cls.dat, aes(Class, Count)) +
  geom_jitter(aes(color=factor(Elevation)), size=2, width=0.15, height=0) +
  scale_color_manual(name ="Elevation", labels=c("400"="400 ft", "1100"="1100 ft", "2000"="2000 ft", "2700"="2700 ft"),
                     values=unique(b.cls.dat$color[order(b.cls.dat$Elevation)])) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  labs(x="Microbial Classes", y="Relative Abundance", title="Bacteria/Archaea & Elevation")
ggsave(b3,filename = "figures/16S_taxa.summary_class_elev_6.22.21.pdf", width=13, height=6, dpi=600)

b4<-ggplot(b.cls_1, aes(Class, Count)) +
  geom_jitter(aes(color=factor(Elevation)), size=2, width=0.15, height=0) +
  scale_color_manual(name ="Elevation", labels=c("400"="400 ft", "1100"="1100 ft", "2000"="2000 ft", "2700"="2700 ft"),
                     values=unique(b.cls_1$color[order(b.cls_1$Elevation)])) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  labs(x="Microbial Classes", y="Relative Abundance", title="Bacteria/Archaea & Elevation", subtitle="Includes Classes > 1% Relative Abundance")
ggsave(b4,filename = "figures/16S_taxa.summary_class.1perc_elev_6.22.21.pdf", width=10, height=6, dpi=600)

# Genera
bac_gen<- as.data.frame(dcast(bac.tax.m,SampleID~Genus, value.var="Count", fun.aggregate=sum)) ### 
head(bac_gen) # counts by genla + elevation
rownames(bac_gen)<-bac_gen$SampleID

bac.gen_RA<-data.frame(decostand(bac_gen[,-c(1)], method="total", MARGIN=1, na.rm=TRUE))  
# relative abundance of taxa data where everything is divided by margin total (default MARGIN = 1 = rows) -- rows = samples
rowSums(bac.gen_RA)
head(bac.gen_RA)
bac.gen_RA$SampleID<-rownames(bac.gen_RA)
head(bac.gen_RA)
bac.gen_m<-melt(bac.gen_RA)
head(bac.gen_m)

names(bac.gen_m)[which(names(bac.gen_m) == "variable")] <- "Genus"
names(bac.gen_m)[which(names(bac.gen_m) == "value")] <- "Count"
head(bac.gen_m) ## relative abundance based on sum of counts by genla!
b.gen.dat<-merge(bac.gen_m, metadata, by="SampleID")
head(b.gen.dat)
b.gen.dat$Elevation<-factor(b.gen.dat$Elevation, levels=c("2700", "2000", "1100", "400"))

b.gen_5<-subset(b.gen.dat, (Count)>(5/100)) ## DROP BACTERIAL genLA that are less than 1% abundant!!!!!!1\
head(b.gen_5)

b.gen_10<-subset(b.gen.dat, (Count)>(10/100)) ## DROP BACTERIAL genLA that are less than 1% abundant!!!!!!1\
head(b.gen_10)

# 16S Genera

b5<-ggplot(b.gen.dat, aes(Genus, Count)) +
  geom_jitter(aes(color=factor(Elevation)), size=2, width=0.15, height=0) +
  scale_color_manual(name ="Elevation", labels=c("400"="400 ft", "1100"="1100 ft", "2000"="2000 ft", "2700"="2700 ft"),
                     values=unique(b.gen.dat$color[order(b.gen.dat$Elevation)])) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  labs(x="Microbial Genera", y="Relative Abundance", title="Bacteria/Archaea & Elevation")
ggsave(b5,filename = "figures/16S_taxa.summary_Genus_elev_6.22.21.pdf", width=35, height=6, dpi=600)

b6<-ggplot(b.gen_5, aes(Genus, Count)) +
  geom_jitter(aes(color=factor(Elevation)), size=2, width=0.15, height=0) +
  scale_color_manual(name ="Elevation", labels=c("400"="400 ft", "1100"="1100 ft", "2000"="2000 ft", "2700"="2700 ft"),
                     values=unique(b.gen_5$color[order(b.gen_5$Elevation)])) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  labs(x="Microbial Genera", y="Relative Abundance", title="Bacteria/Archaea & Elevation", subtitle="Includes Genera > 5% Relative Abundance")
ggsave(b6,filename = "figures/16S_taxa.summary_Genus.5perc_elev_6.22.21.pdf", width=10, height=6, dpi=600)

# Grouped boxplots, no jitter points: do the following commented out block
b6.1<-ggplot(b.gen_5, aes(Genus, Count), fill=Elevation) +
  scale_fill_manual(name ="Elevation",
                    values=unique(b.gen_5$color[order(b.gen_5$Elevation)])) +
  geom_boxplot(aes(fill=factor(Elevation))) + theme_classic() +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  labs(x="Microbial Genera", y="Relative Abundance", title="Bacteria/Archaea & Elevation", subtitle="Includes Genera > 5% Relative Abundance")
ggsave(b6.1,filename = "figures/16S_taxa.summary_Genus.5perc.2_Elev_6.22.21.pdf", width=12, height=6, dpi=600)

b7<-ggplot(b.gen_10, aes(Genus, Count)) +
  geom_jitter(aes(color=factor(Elevation)), size=2, width=0.15, height=0) +
  scale_color_manual(name ="Elevation", labels=c("400"="400 ft", "1100"="1100 ft", "2000"="2000 ft", "2700"="2700 ft"),
                     values=unique(b.gen_10$color[order(b.gen_10$Elevation)])) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  labs(x="Microbial Genera", y="Relative Abundance", title="Bacteria/Archaea & Elevation", subtitle="Includes Genera > 10% Relative Abundance")
ggsave(b7,filename = "figures/16S_taxa.summary_Genus.10perc_elev_6.22.21.pdf", width=10, height=6, dpi=600)

# Grouped boxplots, no jitter points: do the following commented out block
b7.1<-ggplot(b.gen_10, aes(Genus, Count), fill=Elevation) +
  scale_fill_manual(name ="Elevation",
                    values=unique(b.gen_10$color[order(b.gen_10$Elevation)])) +
  geom_boxplot(aes(fill=factor(Elevation))) + theme_classic() +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  labs(x="Microbial Genera", y="Relative Abundance", title="Bacteria/Archaea & Elevation", subtitle="Includes Genera > 10% Relative Abundance")
ggsave(b7.1,filename = "figures/16S_taxa.summary_Genus.10perc.2_Elev_6.22.21.pdf", width=12, height=6, dpi=600)

# Species
bac.spec<- as.data.frame(dcast(bac.tax.m,SampleID~Genus+Species, value.var="Count", fun.aggregate=sum)) ### 
head(bac.spec) # counts by genla + elevation
rownames(bac.spec)<-bac.spec$SampleID

bac.spec_RA<-data.frame(decostand(bac.spec[,-c(1)], method="total", MARGIN=1, na.rm=TRUE))  
# relative abundance of taxa data where everything is divided by margin total (default MARGIN = 1 = rows) -- rows = samples
rowSums(bac.spec_RA)
head(bac.spec_RA)
bac.spec_RA$SampleID<-rownames(bac.spec_RA)
head(bac.spec_RA)
bac.spec_m<-melt(bac.spec_RA)
head(bac.spec_m)

names(bac.spec_m)[which(names(bac.spec_m) == "variable")] <- "Genus_species"
names(bac.spec_m)[which(names(bac.spec_m) == "value")] <- "Count"
head(bac.spec_m) ## relative abundance based on sum of counts by specla!
bac.spec_m$Genus_species<-gsub('_', ' ', bac.spec_m$Genus_species)
bac.spec_m$Genus_species<-gsub('Unknown', 'unknown', bac.spec_m$Genus_species)
bac.spec_m$Genus_species<-gsub('unknown unknown', 'Unknown', bac.spec_m$Genus_species)

head(bac.spec_m)

b.spec.dat<-merge(bac.spec_m, metadata, by="SampleID")
head(b.spec.dat)
b.spec.dat$Elevation<-factor(b.spec.dat$Elevation, levels=c("2700", "2000", "1100", "400"))

b.spec_5<-subset(b.spec.dat, (Count)>(5/100)) ## DROP BACTERIAL specLA that are less than 1% abundant!!!!!!1\
head(b.spec_5)

b.spec_10<-subset(b.spec.dat, (Count)>(10/100)) ## DROP BACTERIAL specLA that are less than 1% abundant!!!!!!1\
head(b.spec_10)

# 16S Species

b8<-ggplot(b.spec.dat, aes(Genus_species, Count)) +
  geom_jitter(aes(color=factor(Elevation)), size=2, width=0.15, height=0) +
  scale_color_manual(name ="Elevation", labels=c("400"="400 ft", "1100"="1100 ft", "2000"="2000 ft", "2700"="2700 ft"),
                     values=unique(b.spec.dat$color[order(b.spec.dat$Elevation)])) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  labs(x="Microbial Species", y="Relative Abundance", title="Bacteria/Archaea & Elevation")
ggsave(b8,filename = "figures/16S_taxa.summary_Species_elev_6.22.21.pdf", width=45, height=6, dpi=600)

b9<-ggplot(b.spec_5, aes(Genus_species, Count)) +
  geom_jitter(aes(color=factor(Elevation)), size=2, width=0.15, height=0) +
  scale_color_manual(name ="Elevation", labels=c("400"="400 ft", "1100"="1100 ft", "2000"="2000 ft", "2700"="2700 ft"),
                     values=unique(b.spec_5$color[order(b.spec_5$Elevation)])) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  labs(x="Microbial Species", y="Relative Abundance", title="Bacteria/Archaea & Elevation", subtitle="Includes Species > 5% Relative Abundance")
ggsave(b9,filename = "figures/16S_taxa.summary_Species.5perc_elev_6.22.21.pdf", width=10, height=6, dpi=600)

# Grouped boxplots, no jitter points: do the following commented out block
b9.1<-ggplot(b.spec_5, aes(Genus_species, Count), fill=Elevation) +
  scale_fill_manual(name ="Elevation",
                    values=unique(b.spec_5$color[order(b.spec_5$Elevation)])) +
  geom_boxplot(aes(fill=factor(Elevation))) + theme_classic() +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  labs(x="Microbial Species", y="Relative Abundance", title="Bacteria/Archaea & Elevation", subtitle="Includes Species > 5% Relative Abundance")
ggsave(b9.1,filename = "figures/16S_taxa.summary_Species.5perc.2_Elev_6.22.21.pdf", width=12, height=6, dpi=600)

b10<-ggplot(b.spec_10, aes(Genus_species, Count)) +
  geom_jitter(aes(color=factor(Elevation)), size=2, width=0.15, height=0) +
  scale_color_manual(name ="Elevation", labels=c("400"="400 ft", "1100"="1100 ft", "2000"="2000 ft", "2700"="2700 ft"),
                     values=unique(b.spec_10$color[order(b.spec_10$Elevation)])) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  labs(x="Microbial Species", y="Relative Abundance", title="Bacteria/Archaea & Elevation", subtitle="Includes Species > 10% Relative Abundance")
ggsave(b10,filename = "figures/16S_taxa.summary_Species.10perc_elev_6.22.21.pdf", width=10, height=6, dpi=600)

b10.1<-ggplot(b.spec_10, aes(Genus_species, Count), fill=Elevation) +
  scale_fill_manual(name ="Elevation", labels=c("400"="400 ft", "1100"="1100 ft", "2000"="2000 ft", "2700"="2700 ft"),
                    values=unique(b.spec_10$color[order(b.spec_10$Elevation)])) +
  geom_boxplot(aes(fill=factor(Elevation))) + theme_classic() +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  labs(x="Microbial Species", y="Relative Abundance", title="Bacteria/Archaea & Elevation", subtitle="Includes Species > 10% Relative Abundance")
ggsave(b10.1,filename = "figures/16S_taxa.summary_Species.10perc.2_Elev_6.22.21.pdf", width=10, height=6, dpi=600)


#### 16S Taxonomic Summaries - specific taxa ####
head(bac.tax.m)

## Acidobacteria
acido.dat<-subset(bac.tax.m, Phylum=="Acidobacteria")

acido_c<- as.data.frame(dcast(acido.dat,SampleID~Class, value.var="Count", fun.aggregate=sum)) ### 
head(acido_c) # counts by phyla + elevation
rownames(acido_c)<-acido_c$SampleID

acido_RA<-data.frame(decostand(acido_c[,-c(1)], method="total", MARGIN=1, na.rm=TRUE))  
# relative abundance of taxa data where everything is divided by margin total (default MARGIN = 1 = rows) -- rows = samples
rowSums(acido_RA)
head(acido_RA)
acido_RA$SampleID<-rownames(acido_RA)
head(acido_RA)
acido_RA.m<-melt(acido_RA)
head(acido_RA.m)

names(acido_RA.m)[which(names(acido_RA.m) == "variable")] <- "Class"
names(acido_RA.m)[which(names(acido_RA.m) == "value")] <- "Count"
head(acido_RA.m) ## relative abundance based on sum of counts by phyla!
acido.meta<-merge(acido_RA.m, metadata, by="SampleID")
head(acido.meta)
acido.meta$Elevation<-factor(acido.meta$Elevation, levels=c("2700", "2000", "1100", "400"))
acido.meta$MonthBin<-factor(acido.meta$MonthBin, levels=c("Early", "Late"))

# 16S phyla

ac1<-ggplot(acido.meta, aes(Class, Count)) +
  geom_jitter(aes(color=factor(Elevation)), size=2, width=0.15, height=0) +
  scale_color_manual(name ="Elevation", labels=c("400"="400 ft", "1100"="1100 ft", "2000"="2000 ft", "2700"="2700 ft"),
                     values=unique(acido.meta$color[order(acido.meta$Elevation)])) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  labs(x="Acidobacteria Classes", y="Relative Abundance", title="Acidobacteria & Elevation")
ggsave(ac1,filename = "figures/Acidobacteria_taxa.summary_elev_8.11.21.pdf", width=10, height=6, dpi=600)

ac2<-ggplot(acido.meta, aes(Class, Count)) +
  geom_jitter(aes(color=factor(MonthBin)), size=2, width=0.15, height=0) +
  scale_color_manual(name ="Time in Dry Season",
                     values=unique(acido.meta$color2[order(acido.meta$MonthBin)])) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  labs(x="Acidobacteria Classes", y="Relative Abundance", title="Acidobacteria & Time in Dry Season")
ggsave(ac2,filename = "figures/Acidobacteria_taxa.summary_monthbin_8.11.21.pdf", width=10, height=6, dpi=600)

## Proteobacteria

proteo.dat<-subset(bac.tax.m, Phylum=="Proteobacteria")

proteo_c<- as.data.frame(dcast(proteo.dat,SampleID~Class, value.var="Count", fun.aggregate=sum)) ### 
head(proteo_c) # counts by phyla + elevation
rownames(proteo_c)<-proteo_c$SampleID

proteo_RA<-data.frame(decostand(proteo_c[,-c(1)], method="total", MARGIN=1, na.rm=TRUE))  
# relative abundance of taxa data where everything is divided by margin total (default MARGIN = 1 = rows) -- rows = samples
rowSums(proteo_RA)
head(proteo_RA)
proteo_RA$SampleID<-rownames(proteo_RA)
head(proteo_RA)
proteo_RA.m<-melt(proteo_RA)
head(proteo_RA.m)

names(proteo_RA.m)[which(names(proteo_RA.m) == "variable")] <- "Class"
names(proteo_RA.m)[which(names(proteo_RA.m) == "value")] <- "Count"
head(proteo_RA.m) ## relative abundance based on sum of counts by phyla!
proteo.meta<-merge(proteo_RA.m, metadata, by="SampleID")
head(proteo.meta)
proteo.meta$Elevation<-factor(proteo.meta$Elevation, levels=c("2700", "2000", "1100", "400"))
proteo.meta$MonthBin<-factor(proteo.meta$MonthBin, levels=c("Early", "Late"))

# 16S phyla

pr1<-ggplot(proteo.meta, aes(Class, Count)) +
  geom_jitter(aes(color=factor(Elevation)), size=2, width=0.15, height=0) +
  scale_color_manual(name ="Elevation", labels=c("400"="400 ft", "1100"="1100 ft", "2000"="2000 ft", "2700"="2700 ft"),
                     values=unique(proteo.meta$color[order(proteo.meta$Elevation)])) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  labs(x="Proteobacteria Classes", y="Relative Abundance", title="Proteobacteria & Elevation")
ggsave(pr1,filename = "figures/Proteobacteria_taxa.summary_elev_8.11.21.pdf", width=10, height=6, dpi=600)

pr2<-ggplot(proteo.meta, aes(Class, Count)) +
  geom_jitter(aes(color=factor(MonthBin)), size=2, width=0.15, height=0) +
  scale_color_manual(name ="Time in Dry Season",
                     values=unique(proteo.meta$color2[order(proteo.meta$MonthBin)])) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  labs(x="Proteobacteria Classes", y="Relative Abundance", title="Proteobacteria & Time in Dry Season")
ggsave(pr2,filename = "figures/Proteobacteria_taxa.summary_monthbin_8.11.21.pdf", width=10, height=6, dpi=600)

#### ITS1 Taxonomic Summaries by Elevation ####

#major_taxa_proportions_tab <- apply(major_taxa_counts_tab, 2, function(x) x/sum(x)*100)

# Phylum
its1_phy<- as.data.frame(dcast(its1.tax.m,SampleID~Phylum, value.var="Count", fun.aggregate=sum)) ### 
head(its1_phy) # counts by phyla + elevation
rownames(its1_phy)<-its1_phy$SampleID

its1.phy_RA<-data.frame(decostand(its1_phy[,-c(1)], method="total", MARGIN=1, na.rm=TRUE))  
# relative abundance of taxa data where everything is divided by margin total (default MARGIN = 1 = rows) -- rows = samples
rowSums(its1.phy_RA)
head(its1.phy_RA)
its1.phy_RA$SampleID<-rownames(its1.phy_RA)
head(its1.phy_RA)
its1.phy_m<-melt(its1.phy_RA)
head(its1.phy_m)

names(its1.phy_m)[which(names(its1.phy_m) == "variable")] <- "Phylum"
names(its1.phy_m)[which(names(its1.phy_m) == "value")] <- "Count"
head(its1.phy_m) ## relative abundance based on sum of counts by phyla!
f.phy.dat<-merge(its1.phy_m, metadata, by="SampleID")
head(f.phy.dat)
f.phy.dat$Elevation<-factor(f.phy.dat$Elevation, levels=c("2700", "2000", "1100", "400"))

f.phy_1<-subset(f.phy.dat, (Count)>(1/100)) ## DROP FUNGAL PHYLA that are less than 1% abundant!!!!!!1\
head(f.phy_1)

# ITS1 phyla

f1<-ggplot(f.phy.dat, aes(Phylum, Count)) +
  geom_jitter(aes(color=factor(Elevation)), size=2, width=0.15, height=0) +
  scale_color_manual(name ="Elevation", labels=c("400"="400 ft", "1100"="1100 ft", "2000"="2000 ft", "2700"="2700 ft"),
                     values=unique(f.phy.dat$color[order(f.phy.dat$Elevation)])) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  labs(x="Fungal Phyla", y="Relative Abundance", title="Fungi & Elevation")
ggsave(f1,filename = "figures/ITS1_taxa.summary_phyla_elev_6.22.21.pdf", width=10, height=6, dpi=600)

f2<-ggplot(f.phy_1, aes(Phylum, Count)) +
  geom_jitter(aes(color=factor(Elevation)), size=2, width=0.15, height=0) +
  scale_color_manual(name ="Elevation", labels=c("400"="400 ft", "1100"="1100 ft", "2000"="2000 ft", "2700"="2700 ft"),
                     values=unique(f.phy_1$color[order(f.phy_1$Elevation)])) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  labs(x="Fungal Phyla", y="Relative Abundance", title="Fungi & Elevation", subtitle="Includes Phyla > 1% Relative Abundance")
ggsave(f2,filename = "figures/ITS1_taxa.summary_phyla.1perc_elev_6.22.21.pdf", width=10, height=6, dpi=600)

# Class
its1_cls<- as.data.frame(dcast(its1.tax.m,SampleID~Class, value.var="Count", fun.aggregate=sum)) ### 
head(its1_cls) # counts by clsla + elevation
rownames(its1_cls)<-its1_cls$SampleID

its1.cls_RA<-data.frame(decostand(its1_cls[,-c(1)], method="total", MARGIN=1, na.rm=TRUE))  
# relative abundance of taxa data where everything is divided by margin total (default MARGIN = 1 = rows) -- rows = samples
rowSums(its1.cls_RA)
head(its1.cls_RA)
its1.cls_RA$SampleID<-rownames(its1.cls_RA)
head(its1.cls_RA)
its1.cls_m<-melt(its1.cls_RA)
head(its1.cls_m)

names(its1.cls_m)[which(names(its1.cls_m) == "variable")] <- "Class"
names(its1.cls_m)[which(names(its1.cls_m) == "value")] <- "Count"
head(its1.cls_m) ## relative abundance based on sum of counts by clsla!
f.cls.dat<-merge(its1.cls_m, metadata, by="SampleID")
head(f.cls.dat)
f.cls.dat$Elevation<-factor(f.cls.dat$Elevation, levels=c("2700", "2000", "1100", "400"))

f.cls_1<-subset(f.cls.dat, (Count)>(1/100)) ## DROP FUNGAL clsLA that are less than 1% abundant!!!!!!1\
head(f.cls_1)

# ITS1 Class

f3<-ggplot(f.cls.dat, aes(Class, Count)) +
  geom_jitter(aes(color=factor(Elevation)), size=2, width=0.15, height=0) +
  scale_color_manual(name ="Elevation", labels=c("400"="400 ft", "1100"="1100 ft", "2000"="2000 ft", "2700"="2700 ft"),
                     values=unique(f.cls.dat$color[order(f.cls.dat$Elevation)])) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  labs(x="Fungal Classes", y="Relative Abundance", title="Fungi & Elevation")
ggsave(f3,filename = "figures/ITS1_taxa.summary_class_elev_6.22.21.pdf", width=13, height=6, dpi=600)

f4<-ggplot(f.cls_1, aes(Class, Count)) +
  geom_jitter(aes(color=factor(Elevation)), size=2, width=0.15, height=0) +
  scale_color_manual(name ="Elevation", labels=c("400"="400 ft", "1100"="1100 ft", "2000"="2000 ft", "2700"="2700 ft"),
                     values=unique(f.cls_1$color[order(f.cls_1$Elevation)])) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  labs(x="Fungal Classes", y="Relative Abundance", title="Fungi & Elevation", subtitle="Includes Classes > 1% Relative Abundance")
ggsave(f4,filename = "figures/ITS1_taxa.summary_class.1perc_elev_6.22.21.pdf", width=10, height=6, dpi=600)

# Genera
its1_gen<- as.data.frame(dcast(its1.tax.m,SampleID~Genus, value.var="Count", fun.aggregate=sum)) ### 
head(its1_gen) # counts by genla + elevation
rownames(its1_gen)<-its1_gen$SampleID

its1.gen_RA<-data.frame(decostand(its1_gen[,-c(1)], method="total", MARGIN=1, na.rm=TRUE))  
# relative abundance of taxa data where everything is divided by margin total (default MARGIN = 1 = rows) -- rows = samples
rowSums(its1.gen_RA)
head(its1.gen_RA)
its1.gen_RA$SampleID<-rownames(its1.gen_RA)
head(its1.gen_RA)
its1.gen_m<-melt(its1.gen_RA)
head(its1.gen_m)

names(its1.gen_m)[which(names(its1.gen_m) == "variable")] <- "Genus"
names(its1.gen_m)[which(names(its1.gen_m) == "value")] <- "Count"
head(its1.gen_m) ## relative abundance based on sum of counts by genla!
f.gen.dat<-merge(its1.gen_m, metadata, by="SampleID")
head(f.gen.dat)
f.gen.dat$Elevation<-factor(f.gen.dat$Elevation, levels=c("2700", "2000", "1100", "400"))

f.gen_5<-subset(f.gen.dat, (Count)>(5/100)) ## DROP FUNGAL genLA that are less than 1% abundant!!!!!!1\
head(f.gen_5)

f.gen_10<-subset(f.gen.dat, (Count)>(10/100)) ## DROP FUNGAL genLA that are less than 1% abundant!!!!!!1\
head(f.gen_10)

# ITS1 Genera

f5<-ggplot(f.gen.dat, aes(Genus, Count)) +
  geom_jitter(aes(color=factor(Elevation)), size=2, width=0.15, height=0) +
  scale_color_manual(name ="Elevation", labels=c("400"="400 ft", "1100"="1100 ft", "2000"="2000 ft", "2700"="2700 ft"),
                     values=unique(f.gen.dat$color[order(f.gen.dat$Elevation)])) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  labs(x="Fungal Genera", y="Relative Abundance", title="Fungi & Elevation")
ggsave(f5,filename = "figures/ITS1_taxa.summary_Genus_elev_6.22.21.pdf", width=30, height=6, dpi=600)

f6<-ggplot(f.gen_5, aes(Genus, Count)) +
  geom_jitter(aes(color=factor(Elevation)), size=2, width=0.15, height=0) +
  scale_color_manual(name ="Elevation", labels=c("400"="400 ft", "1100"="1100 ft", "2000"="2000 ft", "2700"="2700 ft"),
                     values=unique(f.gen_5$color[order(f.gen_5$Elevation)])) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  labs(x="Fungal Genera", y="Relative Abundance", title="Fungi & Elevation", subtitle="Includes Genera > 5% Relative Abundance")
ggsave(f6,filename = "figures/ITS1_taxa.summary_Genus.5perc_elev_6.22.21.pdf", width=10, height=6, dpi=600)

f6.1<-ggplot(f.gen_5, aes(Genus, Count), fill=Elevation) +
  scale_fill_manual(name ="Elevation", labels=c("400"="400 ft", "1100"="1100 ft", "2000"="2000 ft", "2700"="2700 ft"),
                    values=unique(f.gen_5$color[order(f.gen_5$Elevation)])) +
  geom_boxplot(aes(fill=factor(Elevation))) + theme_classic() +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  labs(x="Fungal Genera", y="Relative Abundance", title="Fungi & Elevation", subtitle="Includes Genera > 5% Relative Abundance")
ggsave(f6.1,filename = "figures/ITS1_taxa.summary_Genus.5perc.2_Elev_6.22.21.pdf", width=10, height=6, dpi=600)

f7<-ggplot(f.gen_10, aes(Genus, Count)) +
  geom_jitter(aes(color=factor(Elevation)), size=2, width=0.15, height=0) +
  scale_color_manual(name ="Elevation", labels=c("400"="400 ft", "1100"="1100 ft", "2000"="2000 ft", "2700"="2700 ft"),
                     values=unique(f.gen_10$color[order(f.gen_10$Elevation)])) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  labs(x="Fungal Genera", y="Relative Abundance", title="Fungi & Elevation", subtitle="Includes Genera > 10% Relative Abundance")
ggsave(f7,filename = "figures/ITS1_taxa.summary_Genus.10perc_elev_6.22.21.pdf", width=10, height=6, dpi=600)

f7.1<-ggplot(f.gen_10, aes(Genus, Count), fill=Elevation) +
  scale_fill_manual(name ="Elevation", labels=c("400"="400 ft", "1100"="1100 ft", "2000"="2000 ft", "2700"="2700 ft"),
                    values=unique(f.gen_10$color[order(f.gen_10$Elevation)])) +
  geom_boxplot(aes(fill=factor(Elevation))) + theme_classic() +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  labs(x="Fungal Genera", y="Relative Abundance", title="Fungi & Elevation", subtitle="Includes Genera > 10% Relative Abundance")
ggsave(f7.1,filename = "figures/ITS1_taxa.summary_Genus.10perc.2_Elev_6.22.21.pdf", width=10, height=6, dpi=600)

# Species
its1.spec<- as.data.frame(dcast(its1.tax.m,SampleID~Genus+Species, value.var="Count", fun.aggregate=sum)) ### 
head(its1.spec) # counts by genla + elevation
rownames(its1.spec)<-its1.spec$SampleID

its1.spec_RA<-data.frame(decostand(its1.spec[,-c(1)], method="total", MARGIN=1, na.rm=TRUE))  
# relative abundance of taxa data where everything is divided by margin total (default MARGIN = 1 = rows) -- rows = samples
rowSums(its1.spec_RA)
head(its1.spec_RA)
its1.spec_RA$SampleID<-rownames(its1.spec_RA)
head(its1.spec_RA)
its1.spec_m<-melt(its1.spec_RA)
head(its1.spec_m)

names(its1.spec_m)[which(names(its1.spec_m) == "variable")] <- "Genus_species"
names(its1.spec_m)[which(names(its1.spec_m) == "value")] <- "Count"
head(its1.spec_m) ## relative abundance based on sum of counts by specla!
its1.spec_m$Genus_species<-gsub('_', ' ', its1.spec_m$Genus_species)
its1.spec_m$Genus_species<-gsub('Unknown', 'unknown', its1.spec_m$Genus_species)
its1.spec_m$Genus_species<-gsub('unknown unknown', 'Unknown', its1.spec_m$Genus_species)

head(its1.spec_m)

f.spec.dat<-merge(its1.spec_m, metadata, by="SampleID")
head(f.spec.dat)
f.spec.dat$Elevation<-factor(f.spec.dat$Elevation, levels=c("2700", "2000", "1100", "400"))

f.spec_5<-subset(f.spec.dat, (Count)>(5/100)) ## DROP FUNGAL specLA that are less than 1% abundant!!!!!!1\
head(f.spec_5)

f.spec_10<-subset(f.spec.dat, (Count)>(10/100)) ## DROP FUNGAL specLA that are less than 1% abundant!!!!!!1\
head(f.spec_10)

# ITS1 Species

f8<-ggplot(f.spec.dat, aes(Genus_species, Count)) +
  geom_jitter(aes(color=factor(Elevation)), size=2, width=0.15, height=0) +
  scale_color_manual(name ="Elevation", labels=c("400"="400 ft", "1100"="1100 ft", "2000"="2000 ft", "2700"="2700 ft"),
                     values=unique(f.spec.dat$color[order(f.spec.dat$Elevation)])) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  labs(x="Fungal Species", y="Relative Abundance", title="Fungi & Elevation")
ggsave(f8,filename = "figures/ITS1_taxa.summary_Species_elev_6.22.21.pdf", width=20, height=6, dpi=600)

f9<-ggplot(f.spec_5, aes(Genus_species, Count)) +
  geom_jitter(aes(color=factor(Elevation)), size=2, width=0.15, height=0) +
  scale_color_manual(name ="Elevation", labels=c("400"="400 ft", "1100"="1100 ft", "2000"="2000 ft", "2700"="2700 ft"),
                     values=unique(f.spec_5$color[order(f.spec_5$Elevation)])) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  labs(x="Fungal Species", y="Relative Abundance", title="Fungi & Elevation", subtitle="Includes Species > 5% Relative Abundance")
ggsave(f9,filename = "figures/ITS1_taxa.summary_Species.5perc_elev_6.22.21.pdf", width=10, height=6, dpi=600)

f9.1<-ggplot(f.spec_5, aes(Genus_species, Count), fill=Elevation) +
  scale_fill_manual(name ="Elevation", labels=c("400"="400 ft", "1100"="1100 ft", "2000"="2000 ft", "2700"="2700 ft"),
                    values=unique(f.spec_5$color[order(f.spec_5$Elevation)])) +
  geom_boxplot(aes(fill=factor(Elevation))) + theme_classic() +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  labs(x="Fungal Species", y="Relative Abundance", title="Fungi & Elevation", subtitle="Includes Species > 5% Relative Abundance")
ggsave(f9.1,filename = "figures/ITS1_taxa.summary_Species.5perc.2_Elev_6.22.21.pdf", width=10, height=6, dpi=600)

f10<-ggplot(f.spec_10, aes(Genus_species, Count)) +
  geom_jitter(aes(color=factor(Elevation)), size=2, width=0.15, height=0) +
  scale_color_manual(name ="Elevation", labels=c("400"="400 ft", "1100"="1100 ft", "2000"="2000 ft", "2700"="2700 ft"),
                     values=unique(f.spec_10$color[order(f.spec_10$Elevation)])) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  labs(x="Fungal Species", y="Relative Abundance", title="Fungi & Elevation", subtitle="Includes Species > 10% Relative Abundance")
ggsave(f10,filename = "figures/ITS1_taxa.summary_Species.10perc_elev_6.22.21.pdf", width=10, height=6, dpi=600)

f10.1<-ggplot(f.spec_10, aes(Genus_species, Count), fill=Elevation) +
  scale_fill_manual(name ="Elevation", labels=c("400"="400 ft", "1100"="1100 ft", "2000"="2000 ft", "2700"="2700 ft"),
                    values=unique(f.spec_10$color[order(f.spec_10$Elevation)])) +
  geom_boxplot(aes(fill=factor(Elevation))) + theme_classic() +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  labs(x="Fungal Species", y="Relative Abundance", title="Fungi & Elevation", subtitle="Includes Species > 10% Relative Abundance")
ggsave(f10.1,filename = "figures/ITS1_taxa.summary_Species.10perc.2_Elev_6.22.21.pdf", width=10, height=6, dpi=600)

#### 16S Taxonomic Summaries by Season ####

head(bac_otu_counts)
#major_taxa_proportions_tab <- apply(major_taxa_counts_tab, 2, function(x) x/sum(x)*100)
#bac_otu_counts$OTU_ID<-rownames(bac_otu_counts)
#bac.tax<-merge(bac_otu_counts, bac_otu_tax_ID, by="OTU_ID")
#bac.tax.m<-melt(bac.tax)
head(bac.tax.m)
#colnames(bac.tax.m)[which(names(bac.tax.m) == "value")] <- "Count"
#colnames(bac.tax.m)[which(names(bac.tax.m) == "variable")] <- "SampleID"

# Phylum
# bac_phy<- as.data.frame(dcast(bac.tax.m,SampleID~Phylum, value.var="Count", fun.aggregate=sum)) ### 
# head(bac_phy) # counts by phyla + MonthBin
# rownames(bac_phy)<-bac_phy$SampleID
# 
# bac.phy_RA<-data.frame(decostand(bac_phy[,-c(1)], method="total", MARGIN=1, na.rm=TRUE))  
# # relative abundance of taxa data where everything is divided by margin total (default MARGIN = 1 = rows) -- rows = samples
# rowSums(bac.phy_RA)
# head(bac.phy_RA)
# bac.phy_RA$SampleID<-rownames(bac.phy_RA)
# head(bac.phy_RA)
# bac.phy_m<-melt(bac.phy_RA)
# head(bac.phy_m)
# 
# names(bac.phy_m)[which(names(bac.phy_m) == "variable")] <- "Phylum"
# names(bac.phy_m)[which(names(bac.phy_m) == "value")] <- "Count"
# head(bac.phy_m) ## relative abundance based on sum of counts by phyla!
# b.phy.dat<-merge(bac.phy_m, metadata, by="SampleID")
head(b.phy.dat)
b.phy.dat$MonthBin<-factor(b.phy.dat$MonthBin, levels=c("Summer", "Fall"))

b.phy_1<-subset(b.phy.dat, (Count)>(1/100)) ## DROP BACTERIAL PHYLA that are less than 1% abundant!!!!!!1\
head(b.phy_1)

# 16S phyla

b1a<-ggplot(b.phy.dat, aes(Phylum, Count)) +
  geom_jitter(aes(color=factor(MonthBin)), size=2, width=0.15, height=0) +
  scale_color_manual(name ="Time in Dry Season",
                     values=unique(b.phy.dat$color2[order(b.phy.dat$MonthBin)])) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  labs(x="Microbial Phyla", y="Relative Abundance", title="Bacteria/Archaea & Time in Dry Season")
ggsave(b1a,filename = "figures/16S_taxa.summary_phy.MonthBin_6.22.21.pdf", width=10, height=6, dpi=600)

# Grouped boxplots, no jitter points: do the following commented out block
# b1b<-ggplot(b.phy.dat, aes(Phylum, Count), fill=MonthBin) +
#   scale_fill_manual(name ="Time in Dry Season",
#                      values=unique(b.phy.dat$color2[order(b.phy.dat$MonthBin)])) +
#   geom_boxplot(aes(fill=factor(MonthBin))) + theme_classic() +
#   theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
#   theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
#   labs(x="Microbial Phyla", y="Relative Abundance", title="Bacteria/Archaea & Time in Dry Season")
# ggsave(b1b,filename = "figures/16S_taxa.summary.2_phy.MonthBin_6.22.21.pdf", width=10, height=6, dpi=600)

b2a<-ggplot(b.phy_1, aes(Phylum, Count)) +
  geom_jitter(aes(color=factor(MonthBin)), size=2, width=0.15, height=0) +
  scale_color_manual(name ="Time in Dry Season",
                     values=unique(b.phy_1$color2[order(b.phy_1$MonthBin)])) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  labs(x="Microbial Phyla", y="Relative Abundance", title="Bacteria/Archaea & Time in Dry Season", subtitle="Includes Phyla > 1% Relative Abundance")
ggsave(b2a,filename = "figures/16S_taxa.summary_phyla.1perc_MonthBin_6.22.21.pdf", width=10, height=6, dpi=600)

# Class
# bac_cls<- as.data.frame(dcast(bac.tax.m,SampleID~Class, value.var="Count", fun.aggregate=sum)) ### 
# head(bac_cls) # counts by clsla + MonthBin
# rownames(bac_cls)<-bac_cls$SampleID
# 
# bac.cls_RA<-data.frame(decostand(bac_cls[,-c(1)], method="total", MARGIN=1, na.rm=TRUE))  
# # relative abundance of taxa data where everything is divided by margin total (default MARGIN = 1 = rows) -- rows = samples
# rowSums(bac.cls_RA)
# head(bac.cls_RA)
# bac.cls_RA$SampleID<-rownames(bac.cls_RA)
# head(bac.cls_RA)
# bac.cls_m<-melt(bac.cls_RA)
# head(bac.cls_m)
# 
# names(bac.cls_m)[which(names(bac.cls_m) == "variable")] <- "Class"
# names(bac.cls_m)[which(names(bac.cls_m) == "value")] <- "Count"
# head(bac.cls_m) ## relative abundance based on sum of counts by clsla!
# b.cls.dat<-merge(bac.cls_m, metadata, by="SampleID")
head(b.cls.dat)
b.cls.dat$MonthBin<-factor(b.cls.dat$MonthBin, levels=c("Summer", "Fall"))

b.cls_1<-subset(b.cls.dat, (Count)>(1/100)) ## DROP BACTERIAL clsLA that are less than 1% abundant!!!!!!1\
head(b.cls_1)

# 16S Class

b3a<-ggplot(b.cls.dat, aes(Class, Count)) +
  geom_jitter(aes(color=factor(MonthBin)), size=2, width=0.15, height=0) +
  scale_color_manual(name ="Time in Dry Season",
                     values=unique(b.cls.dat$color2[order(b.cls.dat$MonthBin)])) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  labs(x="Microbial Classes", y="Relative Abundance", title="Bacteria/Archaea & Time in Dry Season")
ggsave(b3a,filename = "figures/16S_taxa.summary_class_MonthBin_6.22.21.pdf", width=13, height=6, dpi=600)

b4a<-ggplot(b.cls_1, aes(Class, Count)) +
  geom_jitter(aes(color=factor(MonthBin)), size=2, width=0.15, height=0) +
  scale_color_manual(name ="Time in Dry Season",
                     values=unique(b.cls_1$color2[order(b.cls_1$MonthBin)])) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  labs(x="Microbial Classes", y="Relative Abundance", title="Bacteria/Archaea & Time in Dry Season", subtitle="Includes Classes > 1% Relative Abundance")
ggsave(b4a,filename = "figures/16S_taxa.summary_class.1perc_MonthBin_6.22.21.pdf", width=10, height=6, dpi=600)

# Genera
# bac_gen<- as.data.frame(dcast(bac.tax.m,SampleID~Genus, value.var="Count", fun.aggregate=sum)) ### 
# head(bac_gen) # counts by genla + MonthBin
# rownames(bac_gen)<-bac_gen$SampleID
# 
# bac.gen_RA<-data.frame(decostand(bac_gen[,-c(1)], method="total", MARGIN=1, na.rm=TRUE))  
# # relative abundance of taxa data where everything is divided by margin total (default MARGIN = 1 = rows) -- rows = samples
# rowSums(bac.gen_RA)
# head(bac.gen_RA)
# bac.gen_RA$SampleID<-rownames(bac.gen_RA)
# head(bac.gen_RA)
# bac.gen_m<-melt(bac.gen_RA)
# head(bac.gen_m)
# 
# names(bac.gen_m)[which(names(bac.gen_m) == "variable")] <- "Genus"
# names(bac.gen_m)[which(names(bac.gen_m) == "value")] <- "Count"
# head(bac.gen_m) ## relative abundance based on sum of counts by genla!
# b.gen.dat<-merge(bac.gen_m, metadata, by="SampleID")
head(b.gen.dat)
b.gen.dat$MonthBin<-factor(b.gen.dat$MonthBin, levels=c("Summer", "Fall"))

b.gen_5<-subset(b.gen.dat, (Count)>(5/100)) ## DROP BACTERIAL genLA that are less than 1% abundant!!!!!!1\
head(b.gen_5)

b.gen_10<-subset(b.gen.dat, (Count)>(10/100)) ## DROP BACTERIAL genLA that are less than 1% abundant!!!!!!1\
head(b.gen_10)

# 16S Genera

b5a<-ggplot(b.gen.dat, aes(Genus, Count)) +
  geom_jitter(aes(color=factor(MonthBin)), size=2, width=0.15, height=0) +
  scale_color_manual(name ="Time in Dry Season",
                     values=unique(b.gen.dat$color2[order(b.gen.dat$MonthBin)])) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  labs(x="Microbial Genera", y="Relative Abundance", title="Bacteria/Archaea & Time in Dry Season")
ggsave(b5a,filename = "figures/16S_taxa.summary_Genus_MonthBin_6.22.21.pdf", width=35, height=6, dpi=600)

b6a<-ggplot(b.gen_5, aes(Genus, Count)) +
  geom_jitter(aes(color=factor(MonthBin)), size=2, width=0.15, height=0) +
  scale_color_manual(name ="Time in Dry Season",
                     values=unique(b.gen_5$color2[order(b.gen_5$MonthBin)])) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  labs(x="Microbial Genera", y="Relative Abundance", title="Bacteria/Archaea & Time in Dry Season", subtitle="Includes Genera > 5% Relative Abundance")
ggsave(b6a,filename = "figures/16S_taxa.summary_Genus.5perc_MonthBin_6.22.21.pdf", width=10, height=6, dpi=600)

# Grouped boxplots, no jitter points: do the following commented out block
b6b<-ggplot(b.gen_5, aes(Genus, Count), fill=MonthBin) +
  scale_fill_manual(name ="Time in Dry Season",
                     values=unique(b.gen_5$color2[order(b.gen_5$MonthBin)])) +
  geom_boxplot(aes(fill=factor(MonthBin))) + theme_classic() +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  labs(x="Microbial Genera", y="Relative Abundance", title="Bacteria/Archaea & Time in Dry Season", subtitle="Includes Genera > 5% Relative Abundance")
ggsave(b6b,filename = "figures/16S_taxa.summary_Genus.5perc.2_MonthBin_6.22.21.pdf", width=10, height=6, dpi=600)

b7a<-ggplot(b.gen_10, aes(Genus, Count)) +
  geom_jitter(aes(color=factor(MonthBin)), size=2, width=0.15, height=0) +
  scale_color_manual(name ="Time in Dry Season",
                     values=unique(b.gen_10$color2[order(b.gen_10$MonthBin)])) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  labs(x="Microbial Genera", y="Relative Abundance", title="Bacteria/Archaea & Time in Dry Season", subtitle="Includes Genera > 10% Relative Abundance")
ggsave(b7a,filename = "figures/16S_taxa.summary_Genus.10perc_MonthBin_6.22.21.pdf", width=10, height=6, dpi=600)

b7b<-ggplot(b.gen_10, aes(Genus, Count), fill=MonthBin) +
  scale_fill_manual(name ="Time in Dry Season",
                    values=unique(b.gen_10$color2[order(b.gen_10$MonthBin)])) +
  geom_boxplot(aes(fill=factor(MonthBin))) + theme_classic() +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  labs(x="Microbial Genera", y="Relative Abundance", title="Bacteria/Archaea & Time in Dry Season", subtitle="Includes Genera > 10% Relative Abundance")
ggsave(b7b,filename = "figures/16S_taxa.summary_Genus.10perc.2_MonthBin_6.22.21.pdf", width=10, height=6, dpi=600)

# Species
# bac.spec<- as.data.frame(dcast(bac.tax.m,SampleID~Genus+Species, value.var="Count", fun.aggregate=sum)) ### 
# head(bac.spec) # counts by genla + MonthBin
# rownames(bac.spec)<-bac.spec$SampleID
# 
# bac.spec_RA<-data.frame(decostand(bac.spec[,-c(1)], method="total", MARGIN=1, na.rm=TRUE))  
# # relative abundance of taxa data where everything is divided by margin total (default MARGIN = 1 = rows) -- rows = samples
# rowSums(bac.spec_RA)
# head(bac.spec_RA)
# bac.spec_RA$SampleID<-rownames(bac.spec_RA)
# head(bac.spec_RA)
# bac.spec_m<-melt(bac.spec_RA)
# head(bac.spec_m)
# 
# names(bac.spec_m)[which(names(bac.spec_m) == "variable")] <- "Genus_species"
# names(bac.spec_m)[which(names(bac.spec_m) == "value")] <- "Count"
# head(bac.spec_m) ## relative abundance based on sum of counts by specla!
# bac.spec_m$Genus_species<-gsub('_', ' ', bac.spec_m$Genus_species)
# bac.spec_m$Genus_species<-gsub('Unknown', 'unknown', bac.spec_m$Genus_species)
# bac.spec_m$Genus_species<-gsub('unknown unknown', 'Unknown', bac.spec_m$Genus_species)
# 
# head(bac.spec_m)
# 
# b.spec.dat<-merge(bac.spec_m, metadata, by="SampleID")
head(b.spec.dat)
b.spec.dat$MonthBin<-factor(b.spec.dat$MonthBin, levels=c("Summer", "Fall"))

b.spec_5<-subset(b.spec.dat, (Count)>(5/100)) ## DROP BACTERIAL specLA that are less than 1% abundant!!!!!!1\
head(b.spec_5)

b.spec_10<-subset(b.spec.dat, (Count)>(10/100)) ## DROP BACTERIAL specLA that are less than 1% abundant!!!!!!1\
head(b.spec_10)

# 16S Species

b8<-ggplot(b.spec.dat, aes(Genus_species, Count)) +
  geom_jitter(aes(color=factor(MonthBin)), size=2, width=0.15, height=0) +
  scale_color_manual(name ="Time in Dry Season",
                     values=unique(b.spec.dat$color2[order(b.spec.dat$MonthBin)])) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  labs(x="Microbial Species", y="Relative Abundance", title="Bacteria/Archaea & Time in Dry Season")
ggsave(b8,filename = "figures/16S_taxa.summary_Species_MonthBin_6.22.21.pdf", width=45, height=6, dpi=600)

b9<-ggplot(b.spec_5, aes(Genus_species, Count)) +
  geom_jitter(aes(color=factor(MonthBin)), size=2, width=0.15, height=0) +
  scale_color_manual(name ="Time in Dry Season",
                     values=unique(b.spec_5$color2[order(b.spec_5$MonthBin)])) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  labs(x="Microbial Species", y="Relative Abundance", title="Bacteria/Archaea & Time in Dry Season", subtitle="Includes Species > 5% Relative Abundance")
ggsave(b9,filename = "figures/16S_taxa.summary_Species.5perc_MonthBin_6.22.21.pdf", width=10, height=6, dpi=600)

b9b<-ggplot(b.spec_5, aes(Genus_species, Count), fill=MonthBin) +
  scale_fill_manual(name ="Time in Dry Season",
                    values=unique(b.spec_5$color2[order(b.spec_5$MonthBin)])) +
  geom_boxplot(aes(fill=factor(MonthBin))) + theme_classic() +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  labs(x="Microbial Species", y="Relative Abundance", title="Bacteria/Archaea & Time in Dry Season", subtitle="Includes Species > 5% Relative Abundance")
ggsave(b9b,filename = "figures/16S_taxa.summary_Species.5perc.2_MonthBin_6.22.21.pdf", width=10, height=6, dpi=600)

b10<-ggplot(b.spec_10, aes(Genus_species, Count)) +
  geom_jitter(aes(color=factor(MonthBin)), size=2, width=0.15, height=0) +
  scale_color_manual(name ="Time in Dry Season",
                     values=unique(b.spec_10$color2[order(b.spec_10$MonthBin)])) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  labs(x="Microbial Species", y="Relative Abundance", title="Bacteria/Archaea & Time in Dry Season", subtitle="Includes Species > 10% Relative Abundance")
ggsave(b10,filename = "figures/16S_taxa.summary_Species.10perc_MonthBin_6.22.21.pdf", width=10, height=6, dpi=600)

b10b<-ggplot(b.spec_10, aes(Genus_species, Count), fill=MonthBin) +
  scale_fill_manual(name ="Time in Dry Season",
                    values=unique(b.spec_10$color2[order(b.spec_10$MonthBin)])) +
  geom_boxplot(aes(fill=factor(MonthBin))) + theme_classic() +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  labs(x="Microbial Species", y="Relative Abundance", title="Bacteria/Archaea & Time in Dry Season", subtitle="Includes Species > 10% Relative Abundance")
ggsave(b10b,filename = "figures/16S_taxa.summary_Species.10perc.2_MonthBin_6.22.21.pdf", width=10, height=6, dpi=600)

#### ITS1 Taxonomic Summaries by Season ####

head(its1_otu_counts)
#major_taxa_proportions_tab <- apply(major_taxa_counts_tab, 2, function(x) x/sum(x)*100)
#bac_otu_counts$OTU_ID<-rownames(bac_otu_counts)
#bac.tax<-merge(bac_otu_counts, bac_otu_tax_ID, by="OTU_ID")
#bac.tax.m<-melt(bac.tax)
head(its1.tax.m)
#colnames(bac.tax.m)[which(names(bac.tax.m) == "value")] <- "Count"
#colnames(bac.tax.m)[which(names(bac.tax.m) == "variable")] <- "SampleID"

# Phylum
# bac_phy<- as.data.frame(dcast(bac.tax.m,SampleID~Phylum, value.var="Count", fun.aggregate=sum)) ### 
# head(bac_phy) # counts by phyla + MonthBin
# rownames(bac_phy)<-bac_phy$SampleID
# 
# bac.phy_RA<-data.frame(decostand(bac_phy[,-c(1)], method="total", MARGIN=1, na.rm=TRUE))  
# # relative abundance of taxa data where everything is divided by margin total (default MARGIN = 1 = rows) -- rows = samples
# rowSums(bac.phy_RA)
# head(bac.phy_RA)
# bac.phy_RA$SampleID<-rownames(bac.phy_RA)
# head(bac.phy_RA)
# bac.phy_m<-melt(bac.phy_RA)
# head(bac.phy_m)
# 
# names(bac.phy_m)[which(names(bac.phy_m) == "variable")] <- "Phylum"
# names(bac.phy_m)[which(names(bac.phy_m) == "value")] <- "Count"
# head(bac.phy_m) ## relative abundance based on sum of counts by phyla!
# b.phy.dat<-merge(bac.phy_m, metadata, by="SampleID")
head(f.phy.dat)
f.phy.dat$MonthBin<-factor(f.phy.dat$MonthBin, levels=c("Summer", "Fall"))

f.phy_1<-subset(f.phy.dat, (Count)>(1/100)) ## DROP Fungal PHYLA that are less than 1% abundant!!!!!!1\
head(f.phy_1)

# ITS1 phyla

f1a<-ggplot(f.phy.dat, aes(Phylum, Count)) +
  geom_jitter(aes(color=factor(MonthBin)), size=2, width=0.15, height=0) +
  scale_color_manual(name ="Time in Dry Season",
                     values=unique(f.phy.dat$color2[order(f.phy.dat$MonthBin)])) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  labs(x="Fungal Phyla", y="Relative Abundance", title="Fungi & Time in Dry Season")
ggsave(f1a,filename = "figures/ITS1_taxa.summary_phy.MonthBin_6.22.21.pdf", width=10, height=6, dpi=600)

# Grouped boxplots, no jitter points: do the following commented out block
# f1b<-ggplot(f.phy.dat, aes(Phylum, Count), fill=MonthBin) +
#   scale_fill_manual(name ="Time in Dry Season",
#                      values=unique(f.phy.dat$color2[order(f.phy.dat$MonthBin)])) +
#   geom_boxplot(aes(fill=factor(MonthBin))) + theme_classic() +
#   theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
#   theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
#   labs(x="Fungal Phyla", y="Relative Abundance", title="Fungi & Time in Dry Season")
# ggsave(f1b,filename = "figures/ITS1_taxa.summary.2_phy.MonthBin_6.22.21.pdf", width=10, height=6, dpi=600)

f2a<-ggplot(f.phy_1, aes(Phylum, Count)) +
  geom_jitter(aes(color=factor(MonthBin)), size=2, width=0.15, height=0) +
  scale_color_manual(name ="Time in Dry Season",
                     values=unique(f.phy_1$color2[order(f.phy_1$MonthBin)])) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  labs(x="Fungal Phyla", y="Relative Abundance", title="Fungi & Time in Dry Season", subtitle="Includes Phyla > 1% Relative Abundance")
ggsave(f2a,filename = "figures/ITS1_taxa.summary_phyla.1perc_MonthBin_6.22.21.pdf", width=10, height=6, dpi=600)

# Class
# bac_cls<- as.data.frame(dcast(bac.tax.m,SampleID~Class, value.var="Count", fun.aggregate=sum)) ### 
# head(bac_cls) # counts by clsla + MonthBin
# rownames(bac_cls)<-bac_cls$SampleID
# 
# bac.cls_RA<-data.frame(decostand(bac_cls[,-c(1)], method="total", MARGIN=1, na.rm=TRUE))  
# # relative abundance of taxa data where everything is divided by margin total (default MARGIN = 1 = rows) -- rows = samples
# rowSums(bac.cls_RA)
# head(bac.cls_RA)
# bac.cls_RA$SampleID<-rownames(bac.cls_RA)
# head(bac.cls_RA)
# bac.cls_m<-melt(bac.cls_RA)
# head(bac.cls_m)
# 
# names(bac.cls_m)[which(names(bac.cls_m) == "variable")] <- "Class"
# names(bac.cls_m)[which(names(bac.cls_m) == "value")] <- "Count"
# head(bac.cls_m) ## relative abundance based on sum of counts by clsla!
# f.cls.dat<-merge(bac.cls_m, metadata, by="SampleID")
head(f.cls.dat)
f.cls.dat$MonthBin<-factor(f.cls.dat$MonthBin, levels=c("Summer", "Fall"))

f.cls_1<-subset(f.cls.dat, (Count)>(1/100)) ## DROP Fungal clsLA that are less than 1% abundant!!!!!!1\
head(f.cls_1)

# ITS1 Class

f3a<-ggplot(f.cls.dat, aes(Class, Count)) +
  geom_jitter(aes(color=factor(MonthBin)), size=2, width=0.15, height=0) +
  scale_color_manual(name ="Time in Dry Season",
                     values=unique(f.cls.dat$color2[order(f.cls.dat$MonthBin)])) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  labs(x="Fungal Classes", y="Relative Abundance", title="Fungi & Time in Dry Season")
ggsave(f3a,filename = "figures/ITS1_taxa.summary_class_MonthBin_6.22.21.pdf", width=13, height=6, dpi=600)

f4a<-ggplot(f.cls_1, aes(Class, Count)) +
  geom_jitter(aes(color=factor(MonthBin)), size=2, width=0.15, height=0) +
  scale_color_manual(name ="Time in Dry Season",
                     values=unique(f.cls_1$color2[order(f.cls_1$MonthBin)])) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  labs(x="Fungal Classes", y="Relative Abundance", title="Fungi & Time in Dry Season", subtitle="Includes Classes > 1% Relative Abundance")
ggsave(f4a,filename = "figures/ITS1_taxa.summary_class.1perc_MonthBin_6.22.21.pdf", width=10, height=6, dpi=600)

# Genera
# bac_gen<- as.data.frame(dcast(bac.tax.m,SampleID~Genus, value.var="Count", fun.aggregate=sum)) ### 
# head(bac_gen) # counts by genla + MonthBin
# rownames(bac_gen)<-bac_gen$SampleID
# 
# bac.gen_RA<-data.frame(decostand(bac_gen[,-c(1)], method="total", MARGIN=1, na.rm=TRUE))  
# # relative abundance of taxa data where everything is divided by margin total (default MARGIN = 1 = rows) -- rows = samples
# rowSums(bac.gen_RA)
# head(bac.gen_RA)
# bac.gen_RA$SampleID<-rownames(bac.gen_RA)
# head(bac.gen_RA)
# bac.gen_m<-melt(bac.gen_RA)
# head(bac.gen_m)
# 
# names(bac.gen_m)[which(names(bac.gen_m) == "variable")] <- "Genus"
# names(bac.gen_m)[which(names(bac.gen_m) == "value")] <- "Count"
# head(bac.gen_m) ## relative abundance based on sum of counts by genla!
# f.gen.dat<-merge(bac.gen_m, metadata, by="SampleID")
head(f.gen.dat)
f.gen.dat$MonthBin<-factor(f.gen.dat$MonthBin, levels=c("Summer", "Fall"))

f.gen_5<-subset(f.gen.dat, (Count)>(5/100)) ## DROP Fungal genLA that are less than 1% abundant!!!!!!1\
head(f.gen_5)

f.gen_10<-subset(f.gen.dat, (Count)>(10/100)) ## DROP Fungal genLA that are less than 1% abundant!!!!!!1\
head(f.gen_10)

# ITS1 Genera

f5a<-ggplot(f.gen.dat, aes(Genus, Count)) +
  geom_jitter(aes(color=factor(MonthBin)), size=2, width=0.15, height=0) +
  scale_color_manual(name ="Time in Dry Season",
                     values=unique(f.gen.dat$color2[order(f.gen.dat$MonthBin)])) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  labs(x="Fungal Genera", y="Relative Abundance", title="Fungi & Time in Dry Season")
ggsave(f5a,filename = "figures/ITS1_taxa.summary_Genus_MonthBin_6.22.21.pdf", width=35, height=6, dpi=600)

f6a<-ggplot(f.gen_5, aes(Genus, Count)) +
  geom_jitter(aes(color=factor(MonthBin)), size=2, width=0.15, height=0) +
  scale_color_manual(name ="Time in Dry Season",
                     values=unique(f.gen_5$color2[order(f.gen_5$MonthBin)])) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  labs(x="Fungal Genera", y="Relative Abundance", title="Fungi & Time in Dry Season", subtitle="Includes Genera > 5% Relative Abundance")
ggsave(f6a,filename = "figures/ITS1_taxa.summary_Genus.5perc_MonthBin_6.22.21.pdf", width=10, height=6, dpi=600)

# Grouped boxplots, no jitter points: do the following commented out block
f6b<-ggplot(f.gen_5, aes(Genus, Count), fill=MonthBin) +
  scale_fill_manual(name ="Time in Dry Season",
                    values=unique(f.gen_5$color2[order(f.gen_5$MonthBin)])) +
  geom_boxplot(aes(fill=factor(MonthBin))) + theme_classic() +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  labs(x="Fungal Genera", y="Relative Abundance", title="Fungi & Time in Dry Season", subtitle="Includes Genera > 5% Relative Abundance")
ggsave(f6b,filename = "figures/ITS1_taxa.summary_Genus.5perc.2_MonthBin_6.22.21.pdf", width=10, height=6, dpi=600)

f7a<-ggplot(f.gen_10, aes(Genus, Count)) +
  geom_jitter(aes(color=factor(MonthBin)), size=2, width=0.15, height=0) +
  scale_color_manual(name ="Time in Dry Season",
                     values=unique(f.gen_10$color2[order(f.gen_10$MonthBin)])) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  labs(x="Fungal Genera", y="Relative Abundance", title="Fungi & Time in Dry Season", subtitle="Includes Genera > 10% Relative Abundance")
ggsave(f7a,filename = "figures/ITS1_taxa.summary_Genus.10perc_MonthBin_6.22.21.pdf", width=10, height=6, dpi=600)

f7b<-ggplot(f.gen_10, aes(Genus, Count), fill=MonthBin) +
  scale_fill_manual(name ="Time in Dry Season",
                    values=unique(f.gen_10$color2[order(f.gen_10$MonthBin)])) +
  geom_boxplot(aes(fill=factor(MonthBin))) + theme_classic() +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  labs(x="Fungal Genera", y="Relative Abundance", title="Fungi & Time in Dry Season", subtitle="Includes Genera > 10% Relative Abundance")
ggsave(f7b,filename = "figures/ITS1_taxa.summary_Genus.10perc.2_MonthBin_6.22.21.pdf", width=10, height=6, dpi=600)

# Species
# bac.spec<- as.data.frame(dcast(bac.tax.m,SampleID~Genus+Species, value.var="Count", fun.aggregate=sum)) ### 
# head(bac.spec) # counts by genla + MonthBin
# rownames(bac.spec)<-bac.spec$SampleID
# 
# bac.spec_RA<-data.frame(decostand(bac.spec[,-c(1)], method="total", MARGIN=1, na.rm=TRUE))  
# # relative abundance of taxa data where everything is divided by margin total (default MARGIN = 1 = rows) -- rows = samples
# rowSums(bac.spec_RA)
# head(bac.spec_RA)
# bac.spec_RA$SampleID<-rownames(bac.spec_RA)
# head(bac.spec_RA)
# bac.spec_m<-melt(bac.spec_RA)
# head(bac.spec_m)
# 
# names(bac.spec_m)[which(names(bac.spec_m) == "variable")] <- "Genus_species"
# names(bac.spec_m)[which(names(bac.spec_m) == "value")] <- "Count"
# head(bac.spec_m) ## relative abundance based on sum of counts by specla!
# bac.spec_m$Genus_species<-gsub('_', ' ', bac.spec_m$Genus_species)
# bac.spec_m$Genus_species<-gsub('Unknown', 'unknown', bac.spec_m$Genus_species)
# bac.spec_m$Genus_species<-gsub('unknown unknown', 'Unknown', bac.spec_m$Genus_species)
# 
# head(bac.spec_m)
# 
# f.spec.dat<-merge(bac.spec_m, metadata, by="SampleID")
head(f.spec.dat)
f.spec.dat$MonthBin<-factor(f.spec.dat$MonthBin, levels=c("Summer", "Fall"))

f.spec_5<-subset(f.spec.dat, (Count)>(5/100)) ## DROP Fungal specLA that are less than 1% abundant!!!!!!1\
head(f.spec_5)

f.spec_10<-subset(f.spec.dat, (Count)>(10/100)) ## DROP Fungal specLA that are less than 1% abundant!!!!!!1\
head(f.spec_10)

# ITS1 Species

f8<-ggplot(f.spec.dat, aes(Genus_species, Count)) +
  geom_jitter(aes(color=factor(MonthBin)), size=2, width=0.15, height=0) +
  scale_color_manual(name ="Time in Dry Season",
                     values=unique(f.spec.dat$color2[order(f.spec.dat$MonthBin)])) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  labs(x="Fungal Species", y="Relative Abundance", title="Fungi & Time in Dry Season")
ggsave(f8,filename = "figures/ITS1_taxa.summary_Species_MonthBin_6.22.21.pdf", width=45, height=6, dpi=600)

f9<-ggplot(f.spec_5, aes(Genus_species, Count)) +
  geom_jitter(aes(color=factor(MonthBin)), size=2, width=0.15, height=0) +
  scale_color_manual(name ="Time in Dry Season",
                     values=unique(f.spec_5$color2[order(f.spec_5$MonthBin)])) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  labs(x="Fungal Species", y="Relative Abundance", title="Fungi & Time in Dry Season", subtitle="Includes Species > 5% Relative Abundance")
ggsave(f9,filename = "figures/ITS1_taxa.summary_Species.5perc_MonthBin_6.22.21.pdf", width=10, height=6, dpi=600)

f9b<-ggplot(f.spec_5, aes(Genus_species, Count), fill=MonthBin) +
  scale_fill_manual(name ="Time in Dry Season",
                    values=unique(f.spec_5$color2[order(f.spec_5$MonthBin)])) +
  geom_boxplot(aes(fill=factor(MonthBin))) + theme_classic() +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  labs(x="Fungal Species", y="Relative Abundance", title="Fungi & Time in Dry Season", subtitle="Includes Species > 5% Relative Abundance")
ggsave(f9b,filename = "figures/ITS1_taxa.summary_Species.5perc.2_MonthBin_6.22.21.pdf", width=10, height=6, dpi=600)

f10<-ggplot(f.spec_10, aes(Genus_species, Count)) +
  geom_jitter(aes(color=factor(MonthBin)), size=2, width=0.15, height=0) +
  scale_color_manual(name ="Time in Dry Season",
                     values=unique(f.spec_10$color2[order(f.spec_10$MonthBin)])) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  labs(x="Fungal Species", y="Relative Abundance", title="Fungi & Time in Dry Season", subtitle="Includes Species > 10% Relative Abundance")
ggsave(f10,filename = "figures/ITS1_taxa.summary_Species.10perc_MonthBin_6.22.21.pdf", width=10, height=6, dpi=600)

f10b<-ggplot(f.spec_10, aes(Genus_species, Count), fill=MonthBin) +
  scale_fill_manual(name ="Time in Dry Season",
                    values=unique(f.spec_10$color2[order(f.spec_10$MonthBin)])) +
  geom_boxplot(aes(fill=factor(MonthBin))) + theme_classic() +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1)) +
  labs(x="Fungal Species", y="Relative Abundance", title="Fungi & Time in Dry Season", subtitle="Includes Species > 10% Relative Abundance")
ggsave(f10b,filename = "figures/ITS1_taxa.summary_Species.10perc.2_MonthBin_6.22.21.pdf", width=10, height=6, dpi=600)

