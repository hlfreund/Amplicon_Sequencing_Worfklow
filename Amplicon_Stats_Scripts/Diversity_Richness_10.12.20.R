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

its1_otu_tax_ID<-subset(its1_otus_counts_taxa, select=c(Kingdom, Phylum, Class, Order, Family, Genus, Species)) # subset only part of the metadata we need
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

#### Import binary OTU tables ####

# 16S first
bac_binary1<- as.data.frame(read.csv("data/Bac_OTU_binarytable.csv"))
head(bac_binary1)
bac_binary1<-subset(bac_binary1, select=-c(taxonomy)) # subset only part of the metadata we need
head(bac_binary1)
dim(bac_binary1)
rownames(bac_binary1)<-bac_otu_ids_clean # set new OTU IDs as rownames
head(bac_binary1)
bac_binary1$OTUID<-NULL

bac_binary_table<-as.data.frame(t(bac_binary1))
head(bac_binary_table)
dim(bac_binary_table)

## binary ITS table

its1_binary1<- as.data.frame(read.csv("data/SierraDustFungiBinary.csv"))
head(its1_binary1)
dim(its1_binary1)
rownames(its1_binary1)<-its1_binary1$OTUID # set new OTU IDs as rownames
head(its1_binary1)
its1_binary1$OTUID<-NULL

its1_binary_table<-as.data.frame(t(its1_binary1))
head(its1_binary_table)
dim(its1_binary_table)
dim(its1_otu_table)

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


## Reorder rows of metadata and all dfs
metadata=metadata[rownames(its1_otu_table),] ## reorder metadata to have same rows as original OTU table
# ^ this indexing method will only work if the two dfs have the same # of rows AND the same row names!
bac_otu_table=bac_otu_table[rownames(metadata),] ## reorder bacterial OTU table to have same rows as ITS1 OTU table + metadata

#bac_binary_table=bac_binary_table[rownames(metadata),] ## reorder binary tables to have same rows as original OTU tables/metadata
#its1_binary_table=its1_binary_table[rownames(metadata),] ## reorder binary tables to have same rows as original OTU tables/metadata

meta_category<-subset(metadata, select=c(SampleID, Year, Month, MonthBin, Site, Elevation)) # subset only part of the metadata we need
meta_quant<-subset(metadata, select=-c(Year, Month, Site, Elevation)) # subset only part of the metadata we need

#### Species Accumulation Curves ####

# fungi
sc1<-specaccum(its1_otu_table,"random")
plot(sc1, ci.type="poly", col="blue", lwd=2, ci.lty=0, ci.col="lightblue")
boxplot(sc1, col="yellow", add=TRUE, pch=20)

# bacteria/archaea
sc2<-specaccum(bac_otu_table,"random")
plot(sc2, ci.type="poly", col="darkgreen", lwd=2, ci.lty=0, ci.col="lightgreen")
boxplot(sc2, col="yellow", add=TRUE, pch=20)

#### Rarefaction Curves ####
# determine read abundance and distribution across samples

sort(colSums(its1_otu_table))
rarecurve(as.matrix(its1_otu_table), step=1000) # fungi

sort(colSums(bac_otu_table))
rarecurve(as.matrix(bac_otu_table), step=1000) # bacteria/archaea

#### Rarefaction ####
# RAREFACTION in vegan: ROWS need to be SITES/samples; COLUMNS are SPECIES (OTUs, ASVs)
# see "Numerical Ecology with R", page 13-14
min_its1<-min(rowSums(its1_otu_table))
min_its1

its_otu_rar<-rrarefy(its1_otu_table,min_its1) 
head(its_otu_rar)
max(rowSums(its_otu_rar)) # should be equal to the min used for rarefaction

min_16s<-min(rowSums(bac_otu_table))
min_16s

bac_otu_rar<-rrarefy(bac_otu_table,min_16s) 
head(bac_otu_rar)
max(rowSums(bac_otu_rar))

#### ITS1 - Alpha Diversity + Species Richness ####
# ** using rarefied data

# ## H <- diversity() -- Shannon entropy function
## exp(H) -- Shannon diversity (exp of entropy)

# Shannon diversity - used rarefied data
# Species richness - used raw data as to not lose rare species
Shan_ent.its1<-vegan::diversity(its_otu_rar, index="shannon") # Shannon entropy
Shan_div.its1<- exp(Shan_ent.its1) # Shannon Diversity aka Hill number 1
div_its1<-data.frame(ITS1_Shannon_Entropy=Shan_ent.its1,ITS1_Shannon_Diversity=Shan_div.its1)
class(div_its1)
div_its1$SampleID<-rownames(div_its1)

S_its1<-specnumber(its1_otu_table) # finds # of species per sample using RAW count data; if MARGIN = 2 it finds frequencies of species
#S_its1 # tells you number of OTUs in every sample
S_its1<-data.frame(ITS1_Species_Richness=specnumber(its1_otu_table), SampleID=rownames(its1_otu_table)) # finds # of species per sample using RAW count data; if MARGIN = 2 it finds frequencies of species
d_r_its1<-merge(div_its1, S_its1, by.x="SampleID", by.y="SampleID")
head(d_r_its1)

S.freq_its1<-specnumber(its1_otu_table, MARGIN = 2) # # finds how many times each ASV appeared in samples (frequency)
S.freq_its1 # # of OTUs that appear across ALL sites/samples

#S_its1a<-data.frame(its1_Species_Richness=specnumber(its1_binary_table)) # finds # of species per sample using RAW count data; if MARGIN = 2 it finds frequencies of species
S_its1 # tells you number of OTUs in every sample


#### 16S - Alpha Diversity + Species Richness ####
# ** using rarefied data

# ## H <- diversity() -- Shannon entropy function
## exp(H) -- Shannon diversity (exp of entropy)

# Shannon diversity - used rarefied data
# Species richness - used raw data as to not lose rare species
Shan_ent.16s<-vegan::diversity(bac_otu_rar, index="shannon") # Shannon entropy
Shan_div.16s<- exp(Shan_ent.16s) # Shannon Diversity aka Hill number 1
div_16s<-data.frame(Bac_Shannon_Entropy=Shan_ent.16s,Bac_Shannon_Diversity=Shan_div.16s)
class(div_16s)
div_16s$SampleID<-rownames(div_16s)

S_16s<-data.frame(Bac_Species_Richness=specnumber(bac_otu_table), SampleID=rownames(bac_otu_table)) # finds # of species per sample using RAW count data; if MARGIN = 2 it finds frequencies of species
d_r_16s<-merge(div_16s, S_16s, by.x="SampleID", by.y="SampleID")
head(d_r_16s)

S.freq_16s<-specnumber(bac_otu_table, MARGIN = 2) # # finds how many times each ASV appeared in samples (frequency)
S.freq_16s # # of OTUs that appear across ALL sites/samples

#S_16s<-data.frame(bac_Species_Richness=specnumber(bac_binary_table)) # finds # of species per sample using RAW count data; if MARGIN = 2 it finds frequencies of species
S_16s # tells you number of OTUs in every sample

#### Merge Alpha div w/ count/taxa data + metadata ####
head(d_r_its1)
head(d_r_16s)

its1.div.metadat <- merge(d_r_its1,metadata, by.x="SampleID", by.y="SampleID")
head(its1.div.metadat)
class(its1.div.metadat)

bac.div.metadat <- merge(d_r_16s,metadata, by.x="SampleID", by.y="SampleID")
head(bac.div.metadat)
class(bac.div.metadat)

## Merging all diversity data together for export table

all_div<-merge(d_r_its1, d_r_16s, by.x="SampleID", by.y="SampleID")
head(all_div)

#write.csv(all_div,"results/Sierra_ALL_Diversity_Richness_4.18.21.csv",row.names=FALSE)

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

season_col <- c("#ffba08", "#9b2226")
names(season_col) <- c("Early", "Late")

#### ITS1 Alpha Diversity Visualization ####
head(its1.div.metadat)
its1.div.metadat$Month2 <- factor(its1.div.metadat$Month, levels = c("July","August","October")) ## reeorder month column so that it will be listed in chronological order in x axis
its1.div.metadat$MonthBin <- factor(its1.div.metadat$MonthBin, levels = c("Early","Late")) 

# shannon entropy by year
its1_ent.yr<-ggplot(its1.div.metadat, aes(x=factor(Year), y=ITS1_Shannon_Entropy, fill=factor(Year))) +geom_boxplot(color="black")+scale_x_discrete()+theme_bw()+scale_fill_manual(values = saturation(wes1, 1))+theme_classic()+
  labs(title = "Fungal Shannon Entropy by Sampling Year", x="Year", y="Shannon Entropy", fill="Year")+theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15))
ggsave(its1_ent.yr,filename = "figures/Fungal_shannon_entropy_rarefied_by_year_4.18.21.pdf", width=8, height=6, dpi=600)

# shannon diversity by year
its1_div.yr<-ggplot(its1.div.metadat, aes(x=factor(Year), y=ITS1_Shannon_Diversity, fill=factor(Year))) +geom_boxplot(color="black")+scale_x_discrete()+theme_bw()+scale_fill_manual(values = saturation(wes1, 1))+theme_classic()+
  labs(title = "Fungal Shannon Diversity by Sampling Year", x="Year", y="Shannon Diversity", fill="Year")+theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15))
ggsave(its1_div.yr,filename = "figures/Fungal_shannon_diversity_rarefied_by_year_4.18.21.pdf", width=8, height=6, dpi=600)

# shannon entropy by month
its1_ent.month<-ggplot(its1.div.metadat, aes(x=Month2, y=ITS1_Shannon_Entropy, fill=Month2)) +geom_boxplot(color="black")+scale_x_discrete()+theme_bw()+scale_fill_manual(values = saturation(wes1, 1))+theme_classic()+
  labs(title = "Fungal Shannon Entropy by Sampling Month", x="Month", y="Shannon Entropy", fill="Month")+theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15))
ggsave(its1_ent.month,filename = "figures/fungal_shannon_entropy_rarefied_by_month_4.18.21.pdf", width=8, height=6, dpi=600)

# shannon diversity by month
its1_div.month<-ggplot(its1.div.metadat, aes(x=Month2, y=ITS1_Shannon_Diversity, fill=Month2)) +geom_boxplot(color="black")+scale_x_discrete()+theme_bw()+scale_fill_manual(values = saturation(wes1, 1))+theme_classic()+
  labs(title = "Fungal Shannon Diversity by Sampling Month", x="Month", y="Shannon Diversity", fill="Month")+theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15))
ggsave(its1_div.month,filename = "figures/fungal_shannon_diversity_rarefied_by_month_4.18.21.pdf", width=8, height=6, dpi=600)

# shannon entropy across elevation
its1_ent.elev<-ggplot(its1.div.metadat, aes(x=factor(Elevation), y=ITS1_Shannon_Entropy, fill=factor(Elevation))) +geom_boxplot(color="black")+scale_x_discrete()+theme_bw()+scale_fill_manual(values = saturation(wes1, 1))+theme_classic()+
  labs(title = "Fungal Shannon Entropy by Elevation", x="Elevation", y="Shannon Entropy", fill="Elevation")+theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15))
ggsave(its1_ent.elev,filename = "figures/fungal_shannon_entropy_rarefied_by_elevation_4.18.21.pdf", width=8, height=6, dpi=600)

# shannon diversity across elevation
its1_div.elev<-ggplot(its1.div.metadat, aes(x=factor(Elevation), y=ITS1_Shannon_Diversity, fill=factor(Elevation))) +geom_boxplot(color="black")+scale_x_discrete()+theme_bw()+scale_fill_manual(values = saturation(wes1, 1))+theme_classic()+
  labs(title = "Fungal Shannon Diversity by Elevation", x="Elevation", y="Shannon Diversity", fill="Elevation")+theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15))
ggsave(its1_div.elev,filename = "figures/fungal_shannon_diversity_rarefied_by_elevation_4.18.21.pdf", width=8, height=6, dpi=600)

# shannon entropy across elevation + year
its1_ent.elev.yr<-ggplot(its1.div.metadat, aes(x=factor(Elevation), y=ITS1_Shannon_Entropy, fill=factor(Year)))+geom_boxplot(color="black")+scale_x_discrete()+theme_bw()+scale_fill_manual(values = saturation(wes1, 1))+theme_classic()+
  labs(title = "Fungal Shannon Entropy by Elevation and Year", x="Elevation", y="Shannon Entropy", fill="Year")+theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15))
ggsave(its1_ent.elev.yr,filename = "figures/fungal_shan_entropy_rarefied_by_elevation.year_4.18.21.pdf", width=8, height=6, dpi=600)

# shannon diversity across elevation + year
its1_div.elev.yr<-ggplot(its1.div.metadat, aes(x=factor(Elevation), y=ITS1_Shannon_Diversity, fill=factor(Year)))+geom_boxplot(color="black")+scale_x_discrete()+theme_bw()+scale_fill_manual(values = saturation(wes1, 1))+theme_classic()+
  labs(title = "Fungal Shannon Diversity by Elevation and Year", x="Elevation", y="Shannon Diversity", fill="Year")+theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15))
ggsave(its1_div.elev.yr,filename = "figures/fungal_shan_diversity_rarefied_by_elevation.year_4.18.21.pdf", width=8, height=6, dpi=600)

# shannon entropy across elevation + year (no box fill, color outlines only)
its_ent.elev.yr1<-ggplot(its1.div.metadat, aes(x=factor(Elevation), y=ITS1_Shannon_Entropy, col=factor(Year)))+geom_boxplot()+scale_x_discrete()+theme_bw()+scale_colour_manual(values = saturation(wes1, 1))+theme_classic()+
  labs(title = "Fungal Shannon Entropy by Elevation and Year", x="Elevation", y="Shannon Entropy", col="Year")+theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15))
ggsave(its_ent.elev.yr1,filename = "figures/fungal_shannon_entropy_rarefied_by_elevation.year.2.4.18.21.pdf", width=8, height=6, dpi=600)

# shannon diversity across elevation + year (no box fill, color outlines only)
its_div.elev.yr1<-ggplot(its1.div.metadat, aes(x=factor(Elevation), y=ITS1_Shannon_Diversity, col=factor(Year)))+geom_boxplot()+scale_x_discrete()+theme_bw()+scale_colour_manual(values = saturation(wes1, 1))+theme_classic()+
  labs(title = "Fungal Shannon Diversity by Elevation and Year", x="Elevation", y="Shannon Diversity", col="Year")+theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15))
ggsave(its_div.elev.yr1,filename = "figures/fungal_shannon_diversity_rarefied_by_elevation.year.2_4.18.21.pdf", width=8, height=6, dpi=600)

dev.off()

### ITS1 Alpha Diversity within Elevations
its1.div_400<-subset(its1.div.metadat, Elevation=='400')
head(its1.div_400)
its1.div_1100<-subset(its1.div.metadat, Elevation=='1100')
its1.div_2000<-subset(its1.div.metadat, Elevation=='2000')
its1.div_2700<-subset(its1.div.metadat, Elevation=='2700')

## its1 diversity w/in 400 ft.
f400.its1.1<-ggplot(its1.div_400, aes(x=factor(Month), y=ITS1_Shannon_Diversity, fill=factor(Month)))+geom_boxplot(color="black")+scale_x_discrete()+theme_bw()+scale_fill_manual(values = saturation(Sag_pal, 0.8))+theme_classic()+
  labs(title = "Fungal Shannon Diversity by Month (400 ft)", x="Month", y="Shannon Diversity", fill="Month")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15))+
  stat_compare_means(method = "anova")
ggsave(f400.its1.1,filename = "figures/400ft_fungal_shan_div_month_6.9.21.pdf", width=8, height=6, dpi=600)

f400.its1.1a<-ggplot(its1.div_400, aes(x=MonthBin, y=ITS1_Shannon_Diversity, fill=MonthBin))+geom_boxplot(color="black")+scale_x_discrete()+theme_bw()+scale_fill_manual(values = saturation(season_col, 1))+theme_classic()+
  labs(title = "Fungal Shannon Diversity (400 ft)", x="Time in Dry Season", y="Shannon Diversity", fill="Time in Dry Season")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15))+
  stat_compare_means(method = "anova",label.y=45)+stat_compare_means(comparisons = list(c(1,2)), method="t.test", hide.ns = FALSE,label = "p.signif")
ggsave(f400.its1.1a,filename = "figures/400ft_fungal_shan_div_binned_month_8.3.21.pdf", width=8, height=6, dpi=600)

f400.its1.2<-ggplot(its1.div_400, aes(x=factor(Year), y=ITS1_Shannon_Diversity, fill=factor(Month)))+geom_boxplot(color="black")+scale_x_discrete()+theme_bw()+scale_fill_manual(values = saturation(Sag_pal, 1))+theme_classic()+
  labs(title = "Fungal Shannon Diversity (400 ft) by Year & Month", x="Year", y="Shannon Diversity", fill="Month")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15))+
  stat_compare_means(method = "anova") +stat_compare_means(comparisons = list(c(1,2), c(2,3),  c(1,3),  c(3,4),  c(2,4),  c(1,4), c(1,5), c(2,5), c(3,5), c(4,5)), method="t.test", hide.ns = TRUE,label = "p.signif")
ggsave(f400.its1.2,filename = "figures/400ft_fungal_shan_div_yr.mon_6.9.21.pdf", width=8, height=6, dpi=600)


## its1 diversity w/in 1100 ft.
f1100.its1.1<-ggplot(its1.div_1100, aes(x=factor(Month), y=ITS1_Shannon_Diversity, fill=factor(Month)))+geom_boxplot(color="black")+scale_x_discrete()+theme_bw()+scale_fill_manual(values = saturation(Sag_pal, 0.8))+theme_classic()+
  labs(title = "Fungal Shannon Diversity by Month (1100 ft)", x="Month", y="Shannon Diversity", fill="Month")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15)) +
  stat_compare_means(method = "anova") +stat_compare_means(comparisons = list(c(1,2), c(2,3),  c(1,3),  c(3,4),  c(2,4),  c(1,4), c(1,5), c(2,5), c(3,5), c(4,5)), method="t.test", hide.ns = TRUE,label = "p.signif")
ggsave(f1100.its1.1,filename = "figures/1100ft_fungal_shan_div_month_6.9.21.pdf", width=8, height=6, dpi=600)

f1100.its1.1a<-ggplot(its1.div_1100, aes(x=MonthBin, y=ITS1_Shannon_Diversity, fill=MonthBin))+geom_boxplot(color="black")+scale_x_discrete()+theme_bw()+scale_fill_manual(values = saturation(season_col, 1))+theme_classic()+
  labs(title = "Fungal Shannon Diversity (1100 ft)", x="Time in Dry Season", y="Shannon Diversity", fill="Time in Dry Season")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15)) +
  stat_compare_means(method = "anova", label.y=110) +stat_compare_means(comparisons = list(c(1,2)), method="t.test", hide.ns = FALSE,label = "p.signif")
ggsave(f1100.its1.1a,filename = "figures/1100ft_fungal_shan_div_binned_month_8.3.21.pdf", width=8, height=6, dpi=600)

f1100.its1.2<-ggplot(its1.div_1100, aes(x=factor(Year), y=ITS1_Shannon_Diversity, fill=factor(Month)))+geom_boxplot(color="black")+scale_x_discrete()+theme_bw()+scale_fill_manual(values = saturation(Sag_pal, 1))+theme_classic()+
  labs(title = "Fungal Shannon Diversity (1100 ft) by Year & Month", x="Year", y="Shannon Diversity", fill="Month")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15)) +
  stat_compare_means(method = "anova") +stat_compare_means(comparisons = list(c(1,2), c(2,3),  c(1,3),  c(3,4),  c(2,4),  c(1,4), c(1,5), c(2,5), c(3,5), c(4,5)), method="t.test", hide.ns = TRUE,label = "p.signif")
ggsave(f1100.its1.2,filename = "figures/1100ft_fungal_shan_div_yr.mon_6.9.21.pdf", width=8, height=6, dpi=600)

## its1 diversity w/in 2000 ft.
f2000.its1.1<-ggplot(its1.div_2000, aes(x=factor(Month), y=ITS1_Shannon_Diversity, fill=factor(Month)))+geom_boxplot(color="black")+scale_x_discrete()+theme_bw()+scale_fill_manual(values = saturation(Sag_pal, 0.8))+theme_classic()+
  labs(title = "Fungal Shannon Diversity by Month (2000 ft)", x="Month", y="Shannon Diversity", fill="Month")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15)) +
  stat_compare_means(method = "anova") +stat_compare_means(comparisons = list(c(1,2), c(2,3),  c(1,3),  c(3,4),  c(2,4),  c(1,4), c(1,5), c(2,5), c(3,5), c(4,5)), method="t.test", hide.ns = TRUE,label = "p.signif")
ggsave(f2000.its1.1,filename = "figures/2000ft_fungal_shan_div_month_6.9.21.pdf", width=8, height=6, dpi=600)

f2000.its1.1a<-ggplot(its1.div_2000, aes(x=MonthBin, y=ITS1_Shannon_Diversity, fill=MonthBin))+geom_boxplot(color="black")+scale_x_discrete()+theme_bw()+scale_fill_manual(values = saturation(season_col, 1))+theme_classic()+
  labs(title = "Fungal Shannon Diversity (2000 ft)", x="Time in Dry Season", y="Shannon Diversity", fill="Time in Dry Season")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15)) +
  stat_compare_means(method = "anova",label.y=35) +stat_compare_means(comparisons = list(c(1,2)), method="t.test", hide.ns = FALSE,label = "p.signif")
ggsave(f2000.its1.1a,filename = "figures/2000ft_fungal_shan_div_binned_month_8.3.21.pdf", width=8, height=6, dpi=600)

f2000.its1.2<-ggplot(its1.div_2000, aes(x=factor(Year), y=ITS1_Shannon_Diversity, fill=factor(Month)))+geom_boxplot(color="black")+scale_x_discrete()+theme_bw()+scale_fill_manual(values = saturation(Sag_pal, 1))+theme_classic()+
  labs(title = "Fungal Shannon Diversity (2000 ft) by Year & Month", x="Year", y="Shannon Diversity", fill="Month")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15)) +
  stat_compare_means(method = "anova") +stat_compare_means(comparisons = list(c(1,2), c(2,3),  c(1,3),  c(3,4),  c(2,4),  c(1,4), c(1,5), c(2,5), c(3,5), c(4,5)), method="t.test", hide.ns = TRUE,label = "p.signif")
ggsave(f2000.its1.2,filename = "figures/2000ft_fungal_shan_div_yr.mon_6.9.21.pdf", width=8, height=6, dpi=600)

## its1 diversity w/in 2700 ft.
f2700.its1.1<-ggplot(its1.div_2700, aes(x=factor(Month), y=ITS1_Shannon_Diversity, fill=factor(Month)))+geom_boxplot(color="black")+scale_x_discrete()+theme_bw()+scale_fill_manual(values = saturation(Sag_pal, 0.8))+theme_classic()+
  labs(title = "Fungal Shannon Diversity by Month (2700 ft)", x="Month", y="Shannon Diversity", fill="Month")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15)) +
  stat_compare_means(method = "anova") +stat_compare_means(comparisons = list(c(1,2), c(2,3),  c(1,3),  c(3,4),  c(2,4),  c(1,4), c(1,5), c(2,5), c(3,5), c(4,5)), method="t.test", hide.ns = TRUE,label = "p.signif")
ggsave(f2700.its1.1,filename = "figures/2700ft_fungal_shan_div_month_6.9.21.pdf", width=8, height=6, dpi=600)

f2700.its1.1a<-ggplot(its1.div_2700, aes(x=MonthBin, y=ITS1_Shannon_Diversity, fill=MonthBin))+geom_boxplot(color="black")+scale_x_discrete()+theme_bw()+scale_fill_manual(values = saturation(season_col, 1))+theme_classic()+
  labs(title = "Fungal Shannon Diversity (2700 ft)", x="Time in Dry Season", y="Shannon Diversity", fill="Time in Dry Season")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15)) +
  stat_compare_means(method = "anova",label.y=19) +stat_compare_means(comparisons = list(c(1,2)), method="t.test", hide.ns = FALSE,label = "p.signif")
ggsave(f2700.its1.1a,filename = "figures/2700ft_fungal_shan_div_binned_month_8.3.21.pdf", width=8, height=6, dpi=600)

f2700.its1.2<-ggplot(its1.div_2700, aes(x=factor(Year), y=ITS1_Shannon_Diversity, fill=factor(Month)))+geom_boxplot(color="black")+scale_x_discrete()+theme_bw()+scale_fill_manual(values = saturation(Sag_pal, 1))+theme_classic()+
  labs(title = "Fungal Shannon Diversity (2700 ft) by Year & Month", x="Year", y="Shannon Diversity", fill="Month")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15)) +
  stat_compare_means(method = "anova") +stat_compare_means(comparisons = list(c(1,2), c(2,3),  c(1,3),  c(3,4),  c(2,4),  c(1,4), c(1,5), c(2,5), c(3,5), c(4,5)), method="t.test", hide.ns = TRUE,label = "p.signif")
ggsave(f2700.its1.2,filename = "figures/2700ft_fungal_shan_div_month.yr_6.9.21.pdf", width=8, height=6, dpi=600)


## ITS1 Alpha Diversity within Elevations + Years by MonthBin
its1.div_2014<-subset(its1.div.metadat, Year=='2014')
its1.div_2015<-subset(its1.div.metadat, Year=='2015')

# Pull out 2014 metadata
its1.div.14_400<-subset(its1.div_2014, Elevation=='400')
head(its1.div.14_400)
its1.div.14_1100<-subset(its1.div_2014, Elevation=='1100')
its1.div.14_2000<-subset(its1.div_2014, Elevation=='2000')
its1.div.14_2700<-subset(its1.div_2014, Elevation=='2700')

# Pull out 2015 metadata
its1.div.15_400<-subset(its1.div_2015, Elevation=='400')
head(its1.div.15_400)
its1.div.15_1100<-subset(its1.div_2015, Elevation=='1100')
its1.div.15_2000<-subset(its1.div_2015, Elevation=='2000')
its1.div.15_2700<-subset(its1.div_2015, Elevation=='2700')

## ITS1 Diversity by Binned Month 2014

its1.mbin.1a<-ggplot(its1.div_2014, aes(x=MonthBin, y=ITS1_Shannon_Diversity, fill=MonthBin))+geom_boxplot(color="black")+scale_x_discrete()+theme_bw()+scale_fill_manual(values = saturation(season_col, 1))+theme_classic()+
  labs(title = "Fungal Shannon Diversity 2014", x="Time in Dry Season", y="Shannon Diversity", fill="Time in Dry Season")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15))+
  stat_compare_means(method = "anova",label.y=86)+stat_compare_means(comparisons = list(c(1,2)), method="t.test", hide.ns = FALSE,label = "p.signif")
ggsave(its1.mbin.1a,filename = "figures/fungal_shan_div_2014_binned_month_8.3.21.pdf", width=8, height=6, dpi=600)

## its1 diversity w/in 400 ft.

f400.its1.1a<-ggplot(its1.div.14_400, aes(x=MonthBin, y=ITS1_Shannon_Diversity, fill=MonthBin))+geom_boxplot(color="black")+scale_x_discrete()+theme_bw()+scale_fill_manual(values = saturation(season_col, 1))+theme_classic()+
  labs(title = "Fungal Shannon Diversity (400 ft) 2014", x="Time in Dry Season", y="Shannon Diversity", fill="Time in Dry Season")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15))+
  stat_compare_means(method = "anova",label.y=45)+stat_compare_means(comparisons = list(c(1,2)), method="t.test", hide.ns = FALSE,label = "p.signif")
ggsave(f400.its1.1a,filename = "figures/400ft_fungal_shan_div_2014_binned_month_8.3.21.pdf", width=8, height=6, dpi=600)

## its1 diversity w/in 1100 ft.

f1100.its1.1a<-ggplot(its1.div.14_1100, aes(x=MonthBin, y=ITS1_Shannon_Diversity, fill=MonthBin))+geom_boxplot(color="black")+scale_x_discrete()+theme_bw()+scale_fill_manual(values = saturation(season_col, 1))+theme_classic()+
  labs(title = "Fungal Shannon Diversity (1100 ft) 2014", x="Time in Dry Season", y="Shannon Diversity", fill="Time in Dry Season")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15)) +
  stat_compare_means(method = "anova", label.y=110) +stat_compare_means(comparisons = list(c(1,2)), method="t.test", hide.ns = FALSE,label = "p.signif")
ggsave(f1100.its1.1a,filename = "figures/1100ft_fungal_shan_div_2014_binned_month_8.3.21.pdf", width=8, height=6, dpi=600)

## its1 diversity w/in 2000 ft.

f2000.its1.1a<-ggplot(its1.div.14_2000, aes(x=MonthBin, y=ITS1_Shannon_Diversity, fill=MonthBin))+geom_boxplot(color="black")+scale_x_discrete()+theme_bw()+scale_fill_manual(values = saturation(season_col, 1))+theme_classic()+
  labs(title = "Fungal Shannon Diversity (2000 ft) 2014", x="Time in Dry Season", y="Shannon Diversity", fill="Time in Dry Season")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15)) +
  stat_compare_means(method = "anova",label.y=35) +stat_compare_means(comparisons = list(c(1,2)), method="t.test", hide.ns = FALSE,label = "p.signif")
ggsave(f2000.its1.1a,filename = "figures/2000ft_fungal_shan_div_2014_binned_month_8.3.21.pdf", width=8, height=6, dpi=600)

## its1 diversity w/in 2700 ft.

f2700.its1.1a<-ggplot(its1.div.14_2700, aes(x=MonthBin, y=ITS1_Shannon_Diversity, fill=MonthBin))+geom_boxplot(color="black")+scale_x_discrete()+theme_bw()+scale_fill_manual(values = saturation(season_col, 1))+theme_classic()+
  labs(title = "Fungal Shannon Diversity (2700 ft) 2014", x="Time in Dry Season", y="Shannon Diversity", fill="Time in Dry Season")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15)) +
  stat_compare_means(method = "anova",label.y=19) +stat_compare_means(comparisons = list(c(1,2)), method="t.test", hide.ns = FALSE,label = "p.signif")
ggsave(f2700.its1.1a,filename = "figures/2700ft_fungal_shan_div_2014_binned_month_8.3.21.pdf", width=8, height=6, dpi=600)


## ITS1 Diversity by Binned Month 2015

its1.mbin.2a<-ggplot(its1.div_2015, aes(x=MonthBin, y=ITS1_Shannon_Diversity, fill=MonthBin))+geom_boxplot(color="black")+scale_x_discrete()+theme_bw()+scale_fill_manual(values = saturation(season_col, 1))+theme_classic()+
  labs(title = "Fungal Shannon Diversity 2015", x="Time in Dry Season", y="Shannon Diversity", fill="Time in Dry Season")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15))+
  stat_compare_means(method = "anova",label.y=105)+stat_compare_means(comparisons = list(c(1,2)), method="t.test", hide.ns = FALSE,label = "p.signif")
ggsave(its1.mbin.2a,filename = "figures/fungal_shan_div_2015_binned_month_8.3.21.pdf", width=8, height=6, dpi=600)

## its1 diversity w/in 400 ft.

f400.its1.15.1a<-ggplot(its1.div.15_400, aes(x=MonthBin, y=ITS1_Shannon_Diversity, fill=MonthBin))+geom_boxplot(color="black")+scale_x_discrete()+theme_bw()+scale_fill_manual(values = saturation(season_col, 1))+theme_classic()+
  labs(title = "Fungal Shannon Diversity (400 ft) 2015", x="Time in Dry Season", y="Shannon Diversity", fill="Time in Dry Season")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15))+
  stat_compare_means(method = "anova",label.y=45)+stat_compare_means(comparisons = list(c(1,2)), method="t.test", hide.ns = FALSE,label = "p.signif")
ggsave(f400.its1.15.1a,filename = "figures/400ft_fungal_shan_div_2015_binned_month_8.3.21.pdf", width=8, height=6, dpi=600)

## its1 diversity w/in 1100 ft.

f1100.its1.15.1a<-ggplot(its1.div.15_1100, aes(x=MonthBin, y=ITS1_Shannon_Diversity, fill=MonthBin))+geom_boxplot(color="black")+scale_x_discrete()+theme_bw()+scale_fill_manual(values = saturation(season_col, 1))+theme_classic()+
  labs(title = "Fungal Shannon Diversity (1100 ft) 2015", x="Time in Dry Season", y="Shannon Diversity", fill="Time in Dry Season")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15)) +
  stat_compare_means(method = "anova", label.y=110) +stat_compare_means(comparisons = list(c(1,2)), method="t.test", hide.ns = FALSE,label = "p.signif")
ggsave(f1100.its1.15.1a,filename = "figures/1100ft_fungal_shan_div_2015_binned_month_8.3.21.pdf", width=8, height=6, dpi=600)

## its1 diversity w/in 2000 ft.

f2000.its1.15.1a<-ggplot(its1.div.15_2000, aes(x=MonthBin, y=ITS1_Shannon_Diversity, fill=MonthBin))+geom_boxplot(color="black")+scale_x_discrete()+theme_bw()+scale_fill_manual(values = saturation(season_col, 1))+theme_classic()+
  labs(title = "Fungal Shannon Diversity (2000 ft) 2015", x="Time in Dry Season", y="Shannon Diversity", fill="Time in Dry Season")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15)) +
  stat_compare_means(method = "anova",label.y=11.2) +stat_compare_means(comparisons = list(c(1,2)), method="t.test", hide.ns = FALSE,label = "p.signif")
ggsave(f2000.its1.15.1a,filename = "figures/2000ft_fungal_shan_div_2015_binned_month_8.3.21.pdf", width=8, height=6, dpi=600)

## its1 diversity w/in 2700 ft.

f2700.its1.15.1a<-ggplot(its1.div.15_2700, aes(x=MonthBin, y=ITS1_Shannon_Diversity, fill=MonthBin))+geom_boxplot(color="black")+scale_x_discrete()+theme_bw()+scale_fill_manual(values = saturation(season_col, 1))+theme_classic()+
  labs(title = "Fungal Shannon Diversity (2700 ft) 2015", x="Time in Dry Season", y="Shannon Diversity", fill="Time in Dry Season")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15)) +
  stat_compare_means(method = "anova",label.y=15) +stat_compare_means(comparisons = list(c(1,2)), method="t.test", hide.ns = FALSE,label = "p.signif")
ggsave(f2700.its1.15.1a,filename = "figures/2700ft_fungal_shan_div_2015_binned_month_8.3.21.pdf", width=8, height=6, dpi=600)

#### 16S Alpha Diversity Visualization ####
head(bac.div.metadat)
bac.div.metadat$Month2 <- factor(bac.div.metadat$Month, levels = c("July","August","October")) ## reeorder month column so that it will be listed in chronological order in x axis
bac.div.metadat$MonthBin<-factor(bac.div.metadat$MonthBin, levels=c("Early", "Late"))

# shannon entropy by year
bac_ent.yr<-ggplot(bac.div.metadat, aes(x=factor(Year), y=Bac_Shannon_Entropy, fill=factor(Year))) +geom_boxplot(color="black")+scale_x_discrete()+theme_bw()+scale_fill_manual(values = saturation(wes1, 1))+theme_classic()+
  labs(title = "Microbial Shannon Entropy by Sampling Year", x="Year", y="Shannon Entropy", fill="Year")+theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15))
ggsave(bac_ent.yr,filename = "figures/Microbial_shannon_entropy_rarefied_by_year_4.18.21.pdf", width=8, height=6, dpi=600)

# shannon diversity by year
bac_div.yr<-ggplot(bac.div.metadat, aes(x=factor(Year), y=Bac_Shannon_Diversity, fill=factor(Year))) +geom_boxplot(color="black")+scale_x_discrete()+theme_bw()+scale_fill_manual(values = saturation(wes1, 1))+theme_classic()+
  labs(title = "Microbial Shannon Diversity by Sampling Year", x="Year", y="Shannon Diversity", fill="Year")+theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15))
ggsave(bac_div.yr,filename = "figures/Microbial_shannon_diversity_rarefied_by_year_4.18.21.pdf", width=8, height=6, dpi=600)

# shannon entropy by month
bac_ent.month<-ggplot(bac.div.metadat, aes(x=Month2, y=Bac_Shannon_Entropy, fill=Month2)) +geom_boxplot(color="black")+scale_x_discrete()+theme_bw()+scale_fill_manual(values = saturation(wes1, 1))+theme_classic()+
  labs(title = "Bacterial Shannon Entropy by Sampling Month", x="Month", y="Shannon Entropy", fill="Month")+theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15))
ggsave(bac_ent.month,filename = "figures/Microbial_shannon_entropy_rarefied_by_month_4.18.21.pdf", width=8, height=6, dpi=600)

# shannon diversity by month
bac_div.month<-ggplot(bac.div.metadat, aes(x=Month2, y=Bac_Shannon_Diversity, fill=Month2)) +geom_boxplot(color="black")+scale_x_discrete()+theme_bw()+scale_fill_manual(values = saturation(wes1, 1))+theme_classic()+
  labs(title = "Microbial Shannon Diversity by Sampling Month", x="Month", y="Shannon Diversity", fill="Month")+theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15))
ggsave(bac_div.month,filename = "figures/Microbial_shannon_diversity_rarefied_by_month_4.18.21.pdf", width=8, height=6, dpi=600)

# shannon entropy across elevation
bac_ent.elev<-ggplot(bac.div.metadat, aes(x=factor(Elevation), y=Bac_Shannon_Entropy, fill=factor(Elevation))) +geom_boxplot(color="black")+scale_x_discrete()+theme_bw()+scale_fill_manual(values = saturation(wes1, 1))+theme_classic()+
  labs(title = "Microbial Shannon Entropy by Elevation", x="Elevation", y="Shannon Entropy", fill="Elevation")+theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15))
ggsave(bac_ent.elev,filename = "figures/Microbial_shannon_entropy_rarefied_by_elevation_4.18.21.pdf", width=8, height=6, dpi=600)

# shannon diversity across elevation
bac_div.elev<-ggplot(bac.div.metadat, aes(x=factor(Elevation), y=Bac_Shannon_Diversity, fill=factor(Elevation))) +geom_boxplot(color="black")+scale_x_discrete()+theme_bw()+scale_fill_manual(values = saturation(wes1, 1))+theme_classic()+
  labs(title = "Microbial Shannon Diversity by Elevation", x="Elevation", y="Shannon Diversity", fill="Elevation")+theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15))
ggsave(bac_div.elev,filename = "figures/Microbial_shannon_diversity_rarefied_by_elevation_4.18.21.pdf", width=8, height=6, dpi=600)

# shannon entropy across elevation + year
bac_ent.elev.yr<-ggplot(bac.div.metadat, aes(x=factor(Elevation), y=Bac_Shannon_Entropy, fill=factor(Year)))+geom_boxplot(color="black")+scale_x_discrete()+theme_bw()+scale_fill_manual(values = saturation(wes1, 1))+theme_classic()+
  labs(title = "Microbial Shannon Entropy by Elevation and Year", x="Elevation", y="Shannon Entropy", fill="Year")+theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15))
ggsave(bac_ent.elev.yr,filename = "figures/Microbial_shan_entropy_rarefied_by_elevation.year_4.18.21.pdf", width=8, height=6, dpi=600)

# shannon diversity across elevation + year
bac_div.elev.yr<-ggplot(bac.div.metadat, aes(x=factor(Elevation), y=Bac_Shannon_Diversity, fill=factor(Year)))+geom_boxplot(color="black")+scale_x_discrete()+theme_bw()+scale_fill_manual(values = saturation(wes1, 1))+theme_classic()+
  labs(title = "Microbial Shannon Diversity by Elevation and Year", x="Elevation", y="Shannon Diversity", fill="Year")+theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15))
ggsave(bac_div.elev.yr,filename = "figures/Microbial_shan_diversity_rarefied_by_elevation.year_4.18.21.pdf", width=8, height=6, dpi=600)

# shannon entropy across elevation + year (no box fill, color outlines only)
bac_ent.elev.yr1<-ggplot(bac.div.metadat, aes(x=factor(Elevation), y=Bac_Shannon_Entropy, col=factor(Year)))+geom_boxplot()+scale_x_discrete()+theme_bw()+scale_colour_manual(values = saturation(wes1, 1))+theme_classic()+
  labs(title = "Microbial Shannon Entropy by Elevation and Year", x="Elevation", y="Shannon Entropy", col="Year")+theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15))
ggsave(bac_ent.elev.yr1,filename = "figures/Microbial_shannon_entropy_rarefied_by_elevation.year.2.4.18.21.pdf", width=8, height=6, dpi=600)

# shannon diversity across elevation + year (no box fill, color outlines only)
bac_div.elev.yr1<-ggplot(bac.div.metadat, aes(x=factor(Elevation), y=Bac_Shannon_Diversity, col=factor(Year)))+geom_boxplot()+scale_x_discrete()+theme_bw()+scale_colour_manual(values = saturation(wes1, 1))+theme_classic()+
  labs(title = "Microbial Shannon Diversity by Elevation and Year", x="Elevation", y="Shannon Diversity", col="Year")+theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15))
ggsave(bac_div.elev.yr1,filename = "figures/Microbial_shannon_diversity_rarefied_by_elevation.year.2_4.18.21.pdf", width=8, height=6, dpi=600)

dev.off()

### Bac Alpha Diversity within Elevations
bac.div_400<-subset(bac.div.metadat, Elevation=='400')
head(bac.div_400)
bac.div_1100<-subset(bac.div.metadat, Elevation=='1100')
bac.div_2000<-subset(bac.div.metadat, Elevation=='2000')
bac.div_2700<-subset(bac.div.metadat, Elevation=='2700')

## bac diversity w/in 400 ft.
f400.bac.1<-ggplot(bac.div_400, aes(x=factor(Month), y=Bac_Shannon_Diversity, fill=factor(Month)))+geom_boxplot(color="black")+scale_x_discrete()+theme_bw()+scale_fill_manual(values = saturation(Sag_pal, 0.8))+theme_classic()+
  labs(title = "Microbial Shannon Diversity by Month (400 ft)", x="Month", y="Shannon Diversity", fill="Month")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15))+
  stat_compare_means(method = "anova")
ggsave(f400.bac.1,filename = "figures/400ft_bacterial_shan_div_month_6.9.21.pdf", width=8, height=6, dpi=600)

f400.bac.1a<-ggplot(bac.div_400, aes(x=MonthBin, y=Bac_Shannon_Diversity, fill=MonthBin))+geom_boxplot(color="black")+scale_x_discrete()+theme_bw()+scale_fill_manual(values = saturation(season_col, 1))+theme_classic()+
  labs(title = "Microbial Shannon Diversity (400 ft)", x="Time in Dry Season", y="Shannon Diversity", fill="Time in Dry Season")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15))+
  stat_compare_means(method = "anova",label.y=505)+stat_compare_means(comparisons = list(c(1,2)), method="t.test", hide.ns = FALSE,label = "p.signif")
ggsave(f400.bac.1a,filename = "figures/400ft_bacterial_shan_div_binned_month_8.3.21.pdf", width=8, height=6, dpi=600)

f400.bac.2<-ggplot(bac.div_400, aes(x=factor(Year), y=Bac_Shannon_Diversity, fill=factor(Month)))+geom_boxplot(color="black")+scale_x_discrete()+theme_bw()+scale_fill_manual(values = saturation(Sag_pal, 1))+theme_classic()+
  labs(title = "Microbial Shannon Diversity (400 ft) by Year & Month", x="Year", y="Shannon Diversity", fill="Month")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15))+
  stat_compare_means(method = "anova") +stat_compare_means(comparisons = list(c(1,2), c(2,3),  c(1,3),  c(3,4),  c(2,4),  c(1,4), c(1,5), c(2,5), c(3,5), c(4,5)), method="t.test", hide.ns = TRUE,label = "p.signif")
ggsave(f400.bac.2,filename = "figures/400ft_bacterial_shan_div_yr.mon_6.9.21.pdf", width=8, height=6, dpi=600)

## bac diversity w/in 1100 ft.
f1100.bac.1<-ggplot(bac.div_1100, aes(x=factor(Month), y=Bac_Shannon_Diversity, fill=factor(Month)))+geom_boxplot(color="black")+scale_x_discrete()+theme_bw()+scale_fill_manual(values = saturation(Sag_pal, 0.8))+theme_classic()+
  labs(title = "Microbial Shannon Diversity by Month (1100 ft)", x="Month", y="Shannon Diversity", fill="Month")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15)) +
  stat_compare_means(method = "anova") +stat_compare_means(comparisons = list(c(1,2), c(2,3),  c(1,3),  c(3,4),  c(2,4),  c(1,4), c(1,5), c(2,5), c(3,5), c(4,5)), method="t.test", hide.ns = TRUE,label = "p.signif")
ggsave(f1100.bac.1,filename = "figures/1100ft_bacterial_shan_div_month_6.9.21.pdf", width=8, height=6, dpi=600)

f1100.bac.1a<-ggplot(bac.div_1100, aes(x=MonthBin, y=Bac_Shannon_Diversity, fill=MonthBin))+geom_boxplot(color="black")+scale_x_discrete()+theme_bw()+scale_fill_manual(values = saturation(season_col, 1))+theme_classic()+
  labs(title = "Microbial Shannon Diversity (1100 ft)", x="Time in Dry Season", y="Shannon Diversity", fill="Time in Dry Season")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15)) +
  stat_compare_means(method = "anova",label.y=280) +stat_compare_means(comparisons = list(c(1,2)), method="t.test", hide.ns = FALSE,label = "p.signif")
ggsave(f1100.bac.1a,filename = "figures/1100ft_bacterial_shan_div_binned_month_8.3.21.pdf", width=8, height=6, dpi=600)

f1100.bac.2<-ggplot(bac.div_1100, aes(x=factor(Year), y=Bac_Shannon_Diversity, fill=factor(Month)))+geom_boxplot(color="black")+scale_x_discrete()+theme_bw()+scale_fill_manual(values = saturation(Sag_pal, 1))+theme_classic()+
  labs(title = "Microbial Shannon Diversity (1100 ft) by Year & Month", x="Year", y="Shannon Diversity", fill="Month")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15)) +
  stat_compare_means(method = "anova") +stat_compare_means(comparisons = list(c(1,2), c(2,3),  c(1,3),  c(3,4),  c(2,4),  c(1,4), c(1,5), c(2,5), c(3,5), c(4,5)), method="t.test", hide.ns = TRUE,label = "p.signif")
ggsave(f1100.bac.2,filename = "figures/1100ft_bacterial_shan_div_yr.mon_6.9.21.pdf", width=8, height=6, dpi=600)

## bac diversity w/in 2000 ft.
f2000.bac.1<-ggplot(bac.div_2000, aes(x=factor(Month), y=Bac_Shannon_Diversity, fill=factor(Month)))+geom_boxplot(color="black")+scale_x_discrete()+theme_bw()+scale_fill_manual(values = saturation(Sag_pal, 0.8))+theme_classic()+
  labs(title = "Microbial Shannon Diversity by Month (2000 ft)", x="Month", y="Shannon Diversity", fill="Month")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15)) +
  stat_compare_means(method = "anova", label.y=310) +stat_compare_means(comparisons = list(c(1,2), c(2,3),  c(1,3),  c(3,4),  c(2,4),  c(1,4), c(1,5), c(2,5), c(3,5), c(4,5)), method="t.test", hide.ns = TRUE,label = "p.signif")
ggsave(f2000.bac.1,filename = "figures/2000ft_bacterial_shan_div_month_6.9.21.pdf", width=8, height=6, dpi=600)

f2000.bac.1a<-ggplot(bac.div_2000, aes(x=MonthBin, y=Bac_Shannon_Diversity, fill=MonthBin))+geom_boxplot(color="black")+scale_x_discrete()+theme_bw()+scale_fill_manual(values = saturation(season_col, 1))+theme_classic()+
  labs(title = "Microbial Shannon Diversity by Month (2000 ft)", x="Time in Dry Season", y="Shannon Diversity", fill="Time in Dry Season")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15)) +
  stat_compare_means(method = "anova", label.y=340) +stat_compare_means(comparisons = list(c(1,2)), method="t.test", hide.ns = FALSE,label = "p.signif")
ggsave(f2000.bac.1a,filename = "figures/2000ft_bacterial_shan_div_binned_month_8.3.21.pdf", width=8, height=6, dpi=600)

f2000.bac.2<-ggplot(bac.div_2000, aes(x=factor(Year), y=Bac_Shannon_Diversity, fill=factor(Month)))+geom_boxplot(color="black")+scale_x_discrete()+theme_bw()+scale_fill_manual(values = saturation(Sag_pal, 1))+theme_classic()+
  labs(title = "Microbial Shannon Diversity (2000 ft) by Year & Month", x="Year", y="Shannon Diversity", fill="Month")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15)) +
  stat_compare_means(method = "anova", label.y=310) +stat_compare_means(comparisons = list(c(1,2), c(2,3),  c(1,3),  c(3,4),  c(2,4),  c(1,4), c(1,5), c(2,5), c(3,5), c(4,5)), method="t.test", hide.ns = TRUE,label = "p.signif")
ggsave(f2000.bac.2,filename = "figures/2000ft_bacterial_shan_div_yr.mon_6.9.21.pdf", width=8, height=6, dpi=600)

## bac richness w/in 2700 ft.
f2700.bac.1<-ggplot(bac.div_2700, aes(x=factor(Month), y=Bac_Shannon_Diversity, fill=factor(Month)))+geom_boxplot(color="black")+scale_x_discrete()+theme_bw()+scale_fill_manual(values = saturation(Sag_pal, 0.8))+theme_classic()+
  labs(title = "Microbial Shannon Diversity by Month (2700 ft)", x="Month", y="Shannon Diversity", fill="Month")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15)) +
  stat_compare_means(method = "anova") +stat_compare_means(comparisons = list(c(1,2), c(2,3),  c(1,3),  c(3,4),  c(2,4),  c(1,4), c(1,5), c(2,5), c(3,5), c(4,5)), method="t.test", hide.ns = TRUE,label = "p.signif")
ggsave(f2700.bac.1,filename = "figures/2700ft_bacterial_shan_div_month_6.9.21.pdf", width=8, height=6, dpi=600)

f2700.bac.1a<-ggplot(bac.div_2700, aes(x=MonthBin, y=Bac_Shannon_Diversity, fill=MonthBin))+geom_boxplot(color="black")+scale_x_discrete()+theme_bw()+scale_fill_manual(values = saturation(season_col, 1))+theme_classic()+
  labs(title = "Microbial Shannon Diversity by Month (2700 ft)", x="Time in Dry Season", y="Shannon Diversity", fill="Time in Dry Season")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15)) +
  stat_compare_means(method = "anova",label.y=62) +stat_compare_means(comparisons = list(c(1,2)), method="t.test", hide.ns = FALSE,label = "p.signif")
ggsave(f2700.bac.1a,filename = "figures/2700ft_bacterial_shan_div_binned_month_8.3.21.pdf", width=8, height=6, dpi=600)

f2700.bac.2<-ggplot(bac.div_2700, aes(x=factor(Year), y=Bac_Shannon_Diversity, fill=factor(Month)))+geom_boxplot(color="black")+scale_x_discrete()+theme_bw()+scale_fill_manual(values = saturation(Sag_pal, 1))+theme_classic()+
  labs(title = "Microbial Shannon Diversity (2700 ft) by Year & Month", x="Year", y="Shannon Diversity", fill="Month")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15)) +
  stat_compare_means(method = "anova") +stat_compare_means(comparisons = list(c(1,2), c(2,3),  c(1,3),  c(3,4),  c(2,4),  c(1,4), c(1,5), c(2,5), c(3,5), c(4,5)), method="t.test", hide.ns = TRUE,label = "p.signif")
ggsave(f2700.bac.2,filename = "figures/2700ft_bacterial_shan_div_month.yr_6.9.21.pdf", width=8, height=6, dpi=600)

## Bac Alpha Diversity within Elevations + Years by MonthBin
bac.div_2014<-subset(bac.div.metadat, Year=='2014')
bac.div_2015<-subset(bac.div.metadat, Year=='2015')

# Pull out 2014 metadata
bac.div.14_400<-subset(bac.div_2014, Elevation=='400')
head(bac.div.14_400)
bac.div.14_1100<-subset(bac.div_2014, Elevation=='1100')
bac.div.14_2000<-subset(bac.div_2014, Elevation=='2000')
bac.div.14_2700<-subset(bac.div_2014, Elevation=='2700')

# Pull out 2015 metadata
bac.div.15_400<-subset(bac.div_2015, Elevation=='400')
head(bac.div.15_400)
bac.div.15_1100<-subset(bac.div_2015, Elevation=='1100')
bac.div.15_2000<-subset(bac.div_2015, Elevation=='2000')
bac.div.15_2700<-subset(bac.div_2015, Elevation=='2700')

## BAC Diversity by Binned Month 2014

bac.mbin.1a<-ggplot(bac.div_2014, aes(x=MonthBin, y=Bac_Shannon_Diversity, fill=MonthBin))+geom_boxplot(color="black")+scale_x_discrete()+theme_bw()+scale_fill_manual(values = saturation(season_col, 1))+theme_classic()+
  labs(title = "Microbial Shannon Diversity 2014", x="Time in Dry Season", y="Shannon Diversity", fill="Time in Dry Season")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15))+
  stat_compare_means(method = "anova",label.y=320)+stat_compare_means(comparisons = list(c(1,2)), method="t.test", hide.ns = FALSE,label = "p.signif")
ggsave(bac.mbin.1a,filename = "figures/microb_shan_div_2014_binned_month_8.3.21.pdf", width=8, height=6, dpi=600)

## bac diversity w/in 400 ft.

f400.bac.1a<-ggplot(bac.div.14_400, aes(x=MonthBin, y=Bac_Shannon_Diversity, fill=MonthBin))+geom_boxplot(color="black")+scale_x_discrete()+theme_bw()+scale_fill_manual(values = saturation(season_col, 1))+theme_classic()+
  labs(title = "Microbial Shannon Diversity (400 ft) 2014", x="Time in Dry Season", y="Shannon Diversity", fill="Time in Dry Season")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15))+
  stat_compare_means(method = "anova",label.y=290)+stat_compare_means(comparisons = list(c(1,2)), method="t.test", hide.ns = FALSE,label = "p.signif")
ggsave(f400.bac.1a,filename = "figures/400ft_microb_shan_div_2014_binned_month_8.3.21.pdf", width=8, height=6, dpi=600)

## bac diversity w/in 1100 ft.

f1100.bac.1a<-ggplot(bac.div.14_1100, aes(x=MonthBin, y=Bac_Shannon_Diversity, fill=MonthBin))+geom_boxplot(color="black")+scale_x_discrete()+theme_bw()+scale_fill_manual(values = saturation(season_col, 1))+theme_classic()+
  labs(title = "Microbial Shannon Diversity (1100 ft) 2014", x="Time in Dry Season", y="Shannon Diversity", fill="Time in Dry Season")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15)) +
  stat_compare_means(method = "anova", label.y=260) +stat_compare_means(comparisons = list(c(1,2)), method="t.test", hide.ns = FALSE,label = "p.signif")
ggsave(f1100.bac.1a,filename = "figures/1100ft_microb_shan_div_2014_binned_month_8.3.21.pdf", width=8, height=6, dpi=600)

## bac diversity w/in 2000 ft.

f2000.bac.1a<-ggplot(bac.div.14_2000, aes(x=MonthBin, y=Bac_Shannon_Diversity, fill=MonthBin))+geom_boxplot(color="black")+scale_x_discrete()+theme_bw()+scale_fill_manual(values = saturation(season_col, 1))+theme_classic()+
  labs(title = "Microbial Shannon Diversity (2000 ft) 2014", x="Time in Dry Season", y="Shannon Diversity", fill="Time in Dry Season")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15)) +
  stat_compare_means(method = "anova",label.y=300) +stat_compare_means(comparisons = list(c(1,2)), method="t.test", hide.ns = FALSE,label = "p.signif")
ggsave(f2000.bac.1a,filename = "figures/2000ft_microb_shan_div_2014_binned_month_8.3.21.pdf", width=8, height=6, dpi=600)

## bac diversity w/in 2700 ft.

f2700.bac.1a<-ggplot(bac.div.14_2700, aes(x=MonthBin, y=Bac_Shannon_Diversity, fill=MonthBin))+geom_boxplot(color="black")+scale_x_discrete()+theme_bw()+scale_fill_manual(values = saturation(season_col, 1))+theme_classic()+
  labs(title = "Microbial Shannon Diversity (2700 ft) 2014", x="Time in Dry Season", y="Shannon Diversity", fill="Time in Dry Season")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15)) +
  stat_compare_means(method = "anova",label.y=19) +stat_compare_means(comparisons = list(c(1,2)), method="t.test", hide.ns = FALSE,label = "p.signif")
ggsave(f2700.bac.1a,filename = "figures/2700ft_microb_shan_div_2014_binned_month_8.3.21.pdf", width=8, height=6, dpi=600)


## Bac Diversity by Binned Month 2015

bac.mbin.2a<-ggplot(bac.div_2015, aes(x=MonthBin, y=Bac_Shannon_Diversity, fill=MonthBin))+geom_boxplot(color="black")+scale_x_discrete()+theme_bw()+scale_fill_manual(values = saturation(season_col, 1))+theme_classic()+
  labs(title = "Microbial Shannon Diversity 2015", x="Time in Dry Season", y="Shannon Diversity", fill="Time in Dry Season")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15))+
  stat_compare_means(method = "anova",label.y=500)+stat_compare_means(comparisons = list(c(1,2)), method="t.test", hide.ns = FALSE,label = "p.signif")
ggsave(bac.mbin.2a,filename = "figures/microb_shan_div_2015_binned_month_8.3.21.pdf", width=8, height=6, dpi=600)

## bac diversity w/in 400 ft.

f400.bac.15.1a<-ggplot(bac.div.15_400, aes(x=MonthBin, y=Bac_Shannon_Diversity, fill=MonthBin))+geom_boxplot(color="black")+scale_x_discrete()+theme_bw()+scale_fill_manual(values = saturation(season_col, 1))+theme_classic()+
  labs(title = "Microbial Shannon Diversity (400 ft) 2015", x="Time in Dry Season", y="Shannon Diversity", fill="Time in Dry Season")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15))+
  stat_compare_means(method = "anova",label.y=500)+stat_compare_means(comparisons = list(c(1,2)), method="t.test", hide.ns = FALSE,label = "p.signif")
ggsave(f400.bac.15.1a,filename = "figures/400ft_microb_shan_div_2015_binned_month_8.3.21.pdf", width=8, height=6, dpi=600)

## bac diversity w/in 1100 ft.

f1100.bac.15.1a<-ggplot(bac.div.15_1100, aes(x=MonthBin, y=Bac_Shannon_Diversity, fill=MonthBin))+geom_boxplot(color="black")+scale_x_discrete()+theme_bw()+scale_fill_manual(values = saturation(season_col, 1))+theme_classic()+
  labs(title = "Microbial Shannon Diversity (1100 ft) 2015", x="Time in Dry Season", y="Shannon Diversity", fill="Time in Dry Season")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15)) +
  stat_compare_means(method = "anova", label.y=200) +stat_compare_means(comparisons = list(c(1,2)), method="t.test", hide.ns = FALSE,label = "p.signif")
ggsave(f1100.bac.15.1a,filename = "figures/1100ft_microb_shan_div_2015_binned_month_8.3.21.pdf", width=8, height=6, dpi=600)

## bac diversity w/in 2000 ft.

f2000.bac.15.1a<-ggplot(bac.div.15_2000, aes(x=MonthBin, y=Bac_Shannon_Diversity, fill=MonthBin))+geom_boxplot(color="black")+scale_x_discrete()+theme_bw()+scale_fill_manual(values = saturation(season_col, 1))+theme_classic()+
  labs(title = "Microbial Shannon Diversity (2000 ft) 2015", x="Time in Dry Season", y="Shannon Diversity", fill="Time in Dry Season")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15)) +
  stat_compare_means(method = "anova",label.y=125) +stat_compare_means(comparisons = list(c(1,2)), method="t.test", hide.ns = FALSE,label = "p.signif")
ggsave(f2000.bac.15.1a,filename = "figures/2000ft_microb_shan_div_2015_binned_month_8.3.21.pdf", width=8, height=6, dpi=600)

## bac diversity w/in 2700 ft.

f2700.bac.15.1a<-ggplot(bac.div.15_2700, aes(x=MonthBin, y=Bac_Shannon_Diversity, fill=MonthBin))+geom_boxplot(color="black")+scale_x_discrete()+theme_bw()+scale_fill_manual(values = saturation(season_col, 1))+theme_classic()+
  labs(title = "Microbial Shannon Diversity (2700 ft) 2015", x="Time in Dry Season", y="Shannon Diversity", fill="Time in Dry Season")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15)) +
  stat_compare_means(method = "anova",label.y=60) +stat_compare_means(comparisons = list(c(1,2)), method="t.test", hide.ns = FALSE,label = "p.signif")
ggsave(f2700.bac.15.1a,filename = "figures/2700ft_microb_shan_div_2015_binned_month_8.3.21.pdf", width=8, height=6, dpi=600)

#### ITS1 Species Richness Visualization ####
head(its1.div.metadat)

its1.div.meta_14<-subset(its1.div.metadat, its1.div.metadat$Year=='2014')
its1.div.meta_15<-subset(its1.div.metadat, its1.div.metadat$Year=='2015')

#its1.div.metadat$Month2 <- factor(its1.div.metadat$Month, levels = c("July","August","October")) ## reeorder month column so that it will be listed in chronological order in x axis
summary(aov(ITS1_Species_Richness~Year, data=its1.div.metadat))
t.test(ITS1_Species_Richness~Year, data=its1.div.metadat)

# species richness by year
its1_SR.yr<-ggplot(its1.div.metadat, aes(x=factor(Year), y=ITS1_Species_Richness, fill=factor(Year))) +geom_boxplot(color="black")+scale_x_discrete()+theme_bw()+scale_fill_manual(values = saturation(wes1, 1))+theme_classic()+
  labs(title = "Fungal Species Richness by Sampling Year", x="Year", y="Species Richness", fill="Year")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15)) + 
  stat_compare_means(comparisons = list(c(1,2)), method="t.test", hide.ns = FALSE,label = "p.signif")+stat_compare_means(method = "anova",label.y=1.5) 
  # ANOVA is significant but t test isn't...
ggsave(its1_SR.yr,filename = "figures/Fungal_SR_by_Year_5.27.21.pdf", width=8, height=6, dpi=600)

# species richness by year + elevation
its1_SR.yr.elev<-ggplot(its1.div.metadat, aes(x = factor(Year), y = ITS1_Species_Richness, fill=factor(Elevation))) + 
  geom_boxplot() + theme_classic() + guides(fill = guide_legend(reverse=TRUE)) + 
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11)) +
  labs(fill="Elevation (ft)")+ylab("Species Richness")+xlab("Year")+scale_fill_manual(values=saturation(fair_cols, 0.9))+stat_compare_means(method = "anova",label.y=1.5) +stat_compare_means(comparisons = list(c(1,2), c(2,3),  c(1,3),  c(3,4),  c(2,4),  c(1,4)), method="t.test", hide.ns = TRUE,label = "p.signif")

ggsave(its1_SR.yr.elev,filename = "figures/ITS1_Spec_Richness_Elev_Year_5.27.21.pdf", width=10, height=8, dpi=600)

# species richness + elevation 2014

its1_SR.14.elev<-ggplot(its1.div.meta_14, aes(x = factor(Elevation), y = ITS1_Species_Richness, fill=factor(Elevation))) + 
  geom_boxplot() + theme_classic() + guides(fill = guide_legend(reverse=TRUE)) + 
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11)) +
  labs(title="Fungal Species Richness x Elevation (2014)",fill="Elevation (ft)")+ylab("Species Richness")+xlab("Elevation (ft)")+scale_fill_manual(values=saturation(fair_cols, 0.9))+stat_compare_means(method = "anova",label.y=1.5) +stat_compare_means(comparisons = list(c(1,2), c(2,3),  c(1,3),  c(3,4),  c(2,4),  c(1,4)), method="t.test", hide.ns = FALSE,label = "p.signif")

ggsave(its1_SR.14.elev,filename = "figures/ITS1_Spec_Richness_Elev_2014_sigbars_5.27.21.pdf", width=10, height=8, dpi=600)

its1_SR.14.elev2<-ggplot(its1.div.meta_14, aes(x = factor(Elevation), y = ITS1_Species_Richness, fill=factor(Elevation))) + 
  geom_boxplot() + theme_classic() + guides(fill = guide_legend(reverse=TRUE)) + 
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11)) +
  labs(title="Fungal Species Richness x Elevation (2014)",fill="Elevation (ft)")+ylab("Species Richness")+xlab("Elevation (ft)")+scale_fill_manual(values=saturation(fair_cols, 0.9))+stat_compare_means(method = "anova",label.y=1.5) +stat_compare_means(comparisons = list(c(1,2), c(2,3),  c(1,3),  c(3,4),  c(2,4),  c(1,4)), method="t.test", label = "p.format")

ggsave(its1_SR.14.elev2,filename = "figures/ITS1_Spec_Richness_Elev_2014_pvals_5.27.21.pdf", width=10, height=8, dpi=600)

# species richness + elevation 2015

its1_SR.15.elev<-ggplot(its1.div.meta_15, aes(x = factor(Elevation), y = ITS1_Species_Richness, fill=factor(Elevation))) + 
  geom_boxplot() + theme_classic() + guides(fill = guide_legend(reverse=TRUE)) + 
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11)) +
  labs(title="Fungal Species Richness x Elevation (2015)",fill="Elevation (ft)")+ylab("Species Richness")+xlab("Elevation (ft)")+scale_fill_manual(values=saturation(fair_cols, 0.9))+stat_compare_means(method = "anova",label.y=1.5) +stat_compare_means(comparisons = list(c(1,2), c(2,3),  c(1,3),  c(3,4),  c(2,4),  c(1,4)), method="t.test", hide.ns = FALSE,label = "p.signif")

ggsave(its1_SR.15.elev,filename = "figures/ITS1_Spec_Richness_Elev_2015_sigbars_5.27.21.pdf", width=10, height=8, dpi=600)

its1_SR.15.elev2<-ggplot(its1.div.meta_15, aes(x = factor(Elevation), y = ITS1_Species_Richness, fill=factor(Elevation))) + 
  geom_boxplot() + theme_classic() + guides(fill = guide_legend(reverse=TRUE)) + 
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11)) +
  labs(title="Fungal Species Richness x Elevation (2015)",fill="Elevation (ft)")+ylab("Species Richness")+xlab("Elevation (ft)")+scale_fill_manual(values=saturation(fair_cols, 0.9))+stat_compare_means(method = "anova",label.y=1.5) +stat_compare_means(comparisons = list(c(1,2), c(2,3),  c(1,3),  c(3,4),  c(2,4),  c(1,4)), method="t.test", label = "p.format")

ggsave(its1_SR.14.elev2,filename = "figures/ITS1_Spec_Richness_Elev_2015_pvals_5.27.21.pdf", width=10, height=8, dpi=600)

itssr_yr_1<-ggarrange(its1_SR.14.elev, its1_SR.15.elev, legend="right", common.legend=TRUE,ncol=2, nrow=1)
ggsave(itssr_yr_1,filename = "figures/ITS1_Spec_Richness_Elev_both_yrs_5.27.21.pdf", width=20, height=8, dpi=600)

itssr_yr_2<-ggarrange(its1_SR.14.elev2, its1_SR.15.elev2, legend="right", common.legend=TRUE,ncol=2, nrow=1)
ggsave(itssr_yr_2,filename = "figures/ITS1_Spec_Richness_Elev_both_yrs_pvals_5.27.21.pdf", width=20, height=8, dpi=600)

# species richness by month
its1_SR.month<-ggplot(its1.div.metadat, aes(x=Month2, y=ITS1_Species_Richness, fill=Month2)) +geom_boxplot(color="black")+scale_x_discrete()+theme_bw()+scale_fill_manual(values = saturation(wes1, 1))+theme_classic()+
  labs(title = "Fungal Species Richness by Sampling Month", x="Month", y="Species Richness", fill="Month")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15)) + 
  stat_compare_means(method = "anova",label.y=1.5) +stat_compare_means(comparisons = list(c(1,2), c(2,3),  c(1,3),  c(3,4),  c(2,4),  c(1,4)), method="t.test", hide.ns = TRUE,label = "p.signif")
ggsave(its1_SR.month,filename = "figures/Fungal_Species_Richness_by_Month_Sierra_3.29.21.pdf", width=8, height=6, dpi=600)

# species richness across elevation
its1_SR.elev<-ggplot(its1.SR, aes(x=factor(Elevation), y=its1_Species_Richness, fill=factor(Elevation))) +geom_boxplot(color="black")+scale_x_discrete()+theme_bw()+scale_fill_manual(values = saturation(wes1, 1))+theme_classic()+
  labs(title = "Fungal Species Richness by Elevation", x="Elevation", y="Species Richness", fill="Elevation")+theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15))
ggsave(its1_SR.elev,filename = "figures/Fungal_Species_Richness_by_Elevation_Sierra_3.29.21.pdf", width=8, height=6, dpi=600)

# species richness across elevation + year
its1_SR.elev.yr<-ggplot(its1.SR, aes(x=factor(Elevation), y=its1_Species_Richness, fill=factor(Year)))+geom_boxplot(color="black")+scale_x_discrete()+theme_bw()+scale_fill_manual(values = saturation(wes1, 1))+theme_classic()+
  labs(title = "Fungal Species Richness by Elevation and Year", x="Elevation", y="Species Richness", fill="Year")+theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15))
ggsave(its1_SR.elev.yr,filename = "figures/Fungal_Species_Richness_by_Elevation_x_Year_Sierra_3.29.21.pdf", width=8, height=6, dpi=600)

### ITS1 Species Richness within Elevations
## its1 richness w/in 400 ft.
f400.its1.3<-ggplot(its1.div_400, aes(x=factor(Month), y=ITS1_Species_Richness, fill=factor(Month)))+geom_boxplot(color="black")+scale_x_discrete()+theme_bw()+scale_fill_manual(values = saturation(Sag_pal, 0.8))+theme_classic()+
  labs(title = "Fungal Species Richness by Month (400 ft)", x="Month", y="Species Richness", fill="Month")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15))+
  stat_compare_means(method = "anova")
ggsave(f400.its1.3,filename = "figures/400ft_fungal_spec_ric_month_6.9.21.pdf", width=8, height=6, dpi=600)

f400.its1.3a<-ggplot(its1.div_400, aes(x=MonthBin, y=ITS1_Species_Richness, fill=MonthBin))+geom_boxplot(color="black")+scale_x_discrete()+theme_bw()+scale_fill_manual(values = saturation(season_col, 1))+theme_classic()+
  labs(title = "Fungal Species Richness (400 ft)", x="Time in Dry Season", y="Species Richness", fill="Time in Dry Season")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15))+
  stat_compare_means(method = "anova",label.y=730)+stat_compare_means(comparisons = list(c(1,2)), method="t.test", hide.ns = FALSE,label = "p.signif")
ggsave(f400.its1.3a,filename = "figures/400ft_fungal_spec_ric_binned_month_8.3.21.pdf", width=8, height=6, dpi=600)

f400.its1.4<-ggplot(its1.div_400, aes(x=factor(Year), y=ITS1_Species_Richness, fill=factor(Month)))+geom_boxplot(color="black")+scale_x_discrete()+theme_bw()+scale_fill_manual(values = saturation(Sag_pal, 1))+theme_classic()+
  labs(title = "Fungal Species Richness (400 ft) by Year & Month", x="Year", y="Species Richness", fill="Month")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15))+
  stat_compare_means(method = "anova") +stat_compare_means(comparisons = list(c(1,2), c(2,3),  c(1,3),  c(3,4),  c(2,4),  c(1,4), c(1,5), c(2,5), c(3,5), c(4,5)), method="t.test", hide.ns = TRUE,label = "p.signif")
ggsave(f400.its1.4,filename = "figures/400ft_fungal_spec_ric_yr.mon_6.9.21.pdf", width=8, height=6, dpi=600)

## its1 Richness w/in 1100 ft.
f1100.its1.3<-ggplot(its1.div_1100, aes(x=factor(Month), y=ITS1_Species_Richness, fill=factor(Month)))+geom_boxplot(color="black")+scale_x_discrete()+theme_bw()+scale_fill_manual(values = saturation(Sag_pal, 0.8))+theme_classic()+
  labs(title = "Fungal Species Richness by Month (1100 ft)", x="Month", y="Species Richness", fill="Month")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15)) +
  stat_compare_means(method = "anova") +stat_compare_means(comparisons = list(c(1,2), c(2,3),  c(1,3),  c(3,4),  c(2,4),  c(1,4), c(1,5), c(2,5), c(3,5), c(4,5)), method="t.test", hide.ns = TRUE,label = "p.signif")
ggsave(f1100.its1.3,filename = "figures/1100ft_fungal_spec_ric_month_6.9.21.pdf", width=8, height=6, dpi=600)

f1100.its1.3a<-ggplot(its1.div_1100, aes(x=MonthBin, y=ITS1_Species_Richness, fill=MonthBin))+geom_boxplot(color="black")+scale_x_discrete()+theme_bw()+scale_fill_manual(values = saturation(season_col, 1))+theme_classic()+
  labs(title = "Fungal Species Richness (1100 ft)", x="Time in Dry Season", y="Species Richness", fill="Time in Dry Season")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15)) +
  stat_compare_means(method = "anova",label.y=950) +stat_compare_means(comparisons = list(c(1,2)), method="t.test", hide.ns = FALSE,label = "p.signif")
ggsave(f1100.its1.3a,filename = "figures/1100ft_fungal_spec_ric_binned_month_8.3.21.pdf", width=8, height=6, dpi=600)

f1100.its1.4<-ggplot(its1.div_1100, aes(x=factor(Year), y=ITS1_Species_Richness, fill=factor(Month)))+geom_boxplot(color="black")+scale_x_discrete()+theme_bw()+scale_fill_manual(values = saturation(Sag_pal, 1))+theme_classic()+
  labs(title = "Fungal Species Richness (1100 ft) by Year & Month", x="Year", y="Species Richness", fill="Month")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15)) +
  stat_compare_means(method = "anova") +stat_compare_means(comparisons = list(c(1,2), c(2,3),  c(1,3),  c(3,4),  c(2,4),  c(1,4), c(1,5), c(2,5), c(3,5), c(4,5)), method="t.test", hide.ns = TRUE,label = "p.signif")
ggsave(f1100.its1.4,filename = "figures/1100ft_fungal_spec_ric_yr.mon_6.9.21.pdf", width=8, height=6, dpi=600)

## its1 Richness w/in 2000 ft.
f2000.its1.3<-ggplot(its1.div_2000, aes(x=factor(Month), y=ITS1_Species_Richness, fill=factor(Month)))+geom_boxplot(color="black")+scale_x_discrete()+theme_bw()+scale_fill_manual(values = saturation(Sag_pal, 0.8))+theme_classic()+
  labs(title = "Fungal Species Richness by Month (2000 ft)", x="Month", y="Species Richness", fill="Month")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15)) +
  stat_compare_means(method = "anova", label.y=700) +stat_compare_means(comparisons = list(c(1,2), c(2,3),  c(1,3),  c(3,4),  c(2,4),  c(1,4), c(1,5), c(2,5), c(3,5), c(4,5)), method="t.test", hide.ns = TRUE,label = "p.signif")
ggsave(f2000.its1.3,filename = "figures/2000ft_fungal_spec_ric_month_6.9.21.pdf", width=8, height=6, dpi=600)

f2000.its1.3a<-ggplot(its1.div_2000, aes(x=MonthBin, y=ITS1_Species_Richness, fill=MonthBin))+geom_boxplot(color="black")+scale_x_discrete()+theme_bw()+scale_fill_manual(values = saturation(season_col, 1))+theme_classic()+
  labs(title = "Fungal Species Richness (2000 ft)", x="Time in Dry Season", y="Species Richness", fill="Time in Dry Season")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15)) +
  stat_compare_means(method = "anova", label.y=700) +stat_compare_means(comparisons = list(c(1,2)), method="t.test", hide.ns = FALSE,label = "p.signif")
ggsave(f2000.its1.3a,filename = "figures/2000ft_fungal_spec_ric_binned_month_8.3.21.pdf", width=8, height=6, dpi=600)

f2000.its1.4<-ggplot(its1.div_2000, aes(x=factor(Year), y=ITS1_Species_Richness, fill=factor(Month)))+geom_boxplot(color="black")+scale_x_discrete()+theme_bw()+scale_fill_manual(values = saturation(Sag_pal, 1))+theme_classic()+
  labs(title = "Fungal Species Richness (2000 ft) by Year & Month", x="Year", y="Species Richness", fill="Month")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15)) +
  stat_compare_means(method = "anova", label.y=700) +stat_compare_means(comparisons = list(c(1,2), c(2,3),  c(1,3),  c(3,4),  c(2,4),  c(1,4), c(1,5), c(2,5), c(3,5), c(4,5)), method="t.test", hide.ns = TRUE,label = "p.signif")
ggsave(f2000.its1.4,filename = "figures/2000ft_fungal_spec_ric_yr.mon_6.9.21.pdf", width=8, height=6, dpi=600)

## its1 Richness w/in 2700 ft.
f2700.its1.3<-ggplot(its1.div_2700, aes(x=factor(Month), y=ITS1_Species_Richness, fill=factor(Month)))+geom_boxplot(color="black")+scale_x_discrete()+theme_bw()+scale_fill_manual(values = saturation(Sag_pal, 0.8))+theme_classic()+
  labs(title = "Fungal Species Richness by Month (2700 ft)", x="Month", y="Species Richness", fill="Month")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15)) +
  stat_compare_means(method = "anova") +stat_compare_means(comparisons = list(c(1,2), c(2,3),  c(1,3),  c(3,4),  c(2,4),  c(1,4), c(1,5), c(2,5), c(3,5), c(4,5)), method="t.test", hide.ns = TRUE,label = "p.signif")
ggsave(f2700.its1.3,filename = "figures/2700ft_fungal_spec_ric_month_6.9.21.pdf", width=8, height=6, dpi=600)

f2700.its1.3a<-ggplot(its1.div_2700, aes(x=MonthBin, y=ITS1_Species_Richness, fill=MonthBin))+geom_boxplot(color="black")+scale_x_discrete()+theme_bw()+scale_fill_manual(values = saturation(season_col, 1))+theme_classic()+
  labs(title = "Fungal Species Richness (2700 ft)", x="Time in Dry Season", y="Species Richness", fill="Time in Dry Season")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15)) +
  stat_compare_means(method = "anova",label.y=400) +stat_compare_means(comparisons = list(c(1,2)), method="t.test", hide.ns = FALSE,label = "p.signif")
ggsave(f2700.its1.3a,filename = "figures/2700ft_fungal_spec_ric_binned_month_8.3.21.pdf", width=8, height=6, dpi=600)

f2700.its1.4<-ggplot(its1.div_2700, aes(x=factor(Year), y=ITS1_Species_Richness, fill=factor(Month)))+geom_boxplot(color="black")+scale_x_discrete()+theme_bw()+scale_fill_manual(values = saturation(Sag_pal, 1))+theme_classic()+
  labs(title = "Fungal Species Richness (2700 ft) by Year & Month", x="Year", y="Species Richness", fill="Month")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15)) +
  stat_compare_means(method = "anova") +stat_compare_means(comparisons = list(c(1,2), c(2,3),  c(1,3),  c(3,4),  c(2,4),  c(1,4), c(1,5), c(2,5), c(3,5), c(4,5)), method="t.test", hide.ns = TRUE,label = "p.signif")
ggsave(f2700.its1.4,filename = "figures/2700ft_fungal_spec_ric_month.yr_6.9.21.pdf", width=8, height=6, dpi=600)


## ITS1 Species Richness within Elevations + Years by MonthBin
# sanity check to see if these dataframes are in the global env.
head(its1.div_2014)
head(its1.div.14_400)

## ITS1 Richness by Binned Month 2014

its1.mbin.2a<-ggplot(its1.div_2014, aes(x=MonthBin, y=ITS1_Species_Richness, fill=MonthBin))+geom_boxplot(color="black")+scale_x_discrete()+theme_bw()+scale_fill_manual(values = saturation(season_col, 1))+theme_classic()+
  labs(title = "Fungal Species Richness 2014", x="Time in Dry Season", y="Species Richness", fill="Time in Dry Season")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15))+
  stat_compare_means(method = "anova",label.y=1000)+stat_compare_means(comparisons = list(c(1,2)), method="t.test", hide.ns = FALSE,label = "p.signif")
ggsave(its1.mbin.2a,filename = "figures/fungal_SR_2014_binned_month_8.3.21.pdf", width=8, height=6, dpi=600)

## its1 species richness w/in 400 ft.

f400.its1.2a<-ggplot(its1.div.14_400, aes(x=MonthBin, y=ITS1_Species_Richness, fill=MonthBin))+geom_boxplot(color="black")+scale_x_discrete()+theme_bw()+scale_fill_manual(values = saturation(season_col, 1))+theme_classic()+
  labs(title = "Fungal Species Richness (400 ft) 2014", x="Time in Dry Season", y="Species Richness", fill="Time in Dry Season")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15))+
  stat_compare_means(method = "anova",label.y=650)+stat_compare_means(comparisons = list(c(1,2)), method="t.test", hide.ns = FALSE,label = "p.signif")
ggsave(f400.its1.2a,filename = "figures/400ft_fungal_SR_2014_binned_month_8.3.21.pdf", width=8, height=6, dpi=600)

## its1 species richness w/in 1100 ft.

f1100.its1.2a<-ggplot(its1.div.14_1100, aes(x=MonthBin, y=ITS1_Species_Richness, fill=MonthBin))+geom_boxplot(color="black")+scale_x_discrete()+theme_bw()+scale_fill_manual(values = saturation(season_col, 1))+theme_classic()+
  labs(title = "Fungal Species Richness (1100 ft) 2014", x="Time in Dry Season", y="Species Richness", fill="Time in Dry Season")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15)) +
  stat_compare_means(method = "anova", label.y=900) +stat_compare_means(comparisons = list(c(1,2)), method="t.test", hide.ns = FALSE,label = "p.signif")
ggsave(f1100.its1.2a,filename = "figures/1100ft_fungal_SR_2014_binned_month_8.3.21.pdf", width=8, height=6, dpi=600)

## its1 species richness w/in 2000 ft.

f2000.its1.2a<-ggplot(its1.div.14_2000, aes(x=MonthBin, y=ITS1_Species_Richness, fill=MonthBin))+geom_boxplot(color="black")+scale_x_discrete()+theme_bw()+scale_fill_manual(values = saturation(season_col, 1))+theme_classic()+
  labs(title = "Fungal Species Richness (2000 ft) 2014", x="Time in Dry Season", y="Species Richness", fill="Time in Dry Season")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15)) +
  stat_compare_means(method = "anova",label.y=650) +stat_compare_means(comparisons = list(c(1,2)), method="t.test", hide.ns = FALSE,label = "p.signif")
ggsave(f2000.its1.2a,filename = "figures/2000ft_fungal_SR_2014_binned_month_8.3.21.pdf", width=8, height=6, dpi=600)

## its1 species richness w/in 2700 ft.

f2700.its1.2a<-ggplot(its1.div.14_2700, aes(x=MonthBin, y=ITS1_Species_Richness, fill=MonthBin))+geom_boxplot(color="black")+scale_x_discrete()+theme_bw()+scale_fill_manual(values = saturation(season_col, 1))+theme_classic()+
  labs(title = "Fungal Species Richness (2700 ft) 2014", x="Time in Dry Season", y="Species Richness", fill="Time in Dry Season")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15)) +
  stat_compare_means(method = "anova",label.y=370) +stat_compare_means(comparisons = list(c(1,2)), method="t.test", hide.ns = FALSE,label = "p.signif")
ggsave(f2700.its1.2a,filename = "figures/2700ft_fungal_SR_2014_binned_month_8.3.21.pdf", width=8, height=6, dpi=600)


## ITS1 Species Richness by Binned Month 2015

its1.mbin.3a<-ggplot(its1.div_2015, aes(x=MonthBin, y=ITS1_Species_Richness, fill=MonthBin))+geom_boxplot(color="black")+scale_x_discrete()+theme_bw()+scale_fill_manual(values = saturation(season_col, 1))+theme_classic()+
  labs(title = "Fungal Species Richness 2015", x="Time in Dry Season", y="Species Richness", fill="Time in Dry Season")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15))+
  stat_compare_means(method = "anova",label.y=860)+stat_compare_means(comparisons = list(c(1,2)), method="t.test", hide.ns = FALSE,label = "p.signif")
ggsave(its1.mbin.3a,filename = "figures/fungal_SR_2015_binned_month_8.3.21.pdf", width=8, height=6, dpi=600)

## its1 species richness w/in 400 ft.

f400.its1.15.2a<-ggplot(its1.div.15_400, aes(x=MonthBin, y=ITS1_Species_Richness, fill=MonthBin))+geom_boxplot(color="black")+scale_x_discrete()+theme_bw()+scale_fill_manual(values = saturation(season_col, 1))+theme_classic()+
  labs(title = "Fungal Species Richness (400 ft) 2015", x="Time in Dry Season", y="Species Richness", fill="Time in Dry Season")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15))+
  stat_compare_means(method = "anova",label.y=720)+stat_compare_means(comparisons = list(c(1,2)), method="t.test", hide.ns = FALSE,label = "p.signif")
ggsave(f400.its1.15.2a,filename = "figures/400ft_fungal_SR_2015_binned_month_8.3.21.pdf", width=8, height=6, dpi=600)

## its1 species richness w/in 1100 ft.

f1100.its1.15.2a<-ggplot(its1.div.15_1100, aes(x=MonthBin, y=ITS1_Species_Richness, fill=MonthBin))+geom_boxplot(color="black")+scale_x_discrete()+theme_bw()+scale_fill_manual(values = saturation(season_col, 1))+theme_classic()+
  labs(title = "Fungal Species Richness (1100 ft) 2015", x="Time in Dry Season", y="Species Richness", fill="Time in Dry Season")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15)) +
  stat_compare_means(method = "anova", label.y=825) +stat_compare_means(comparisons = list(c(1,2)), method="t.test", hide.ns = FALSE,label = "p.signif")
ggsave(f1100.its1.15.2a,filename = "figures/1100ft_fungal_SR_2015_binned_month_8.3.21.pdf", width=8, height=6, dpi=600)

## its1 species richness w/in 2000 ft.

f2000.its1.15.2a<-ggplot(its1.div.15_2000, aes(x=MonthBin, y=ITS1_Species_Richness, fill=MonthBin))+geom_boxplot(color="black")+scale_x_discrete()+theme_bw()+scale_fill_manual(values = saturation(season_col, 1))+theme_classic()+
  labs(title = "Fungal Species Richness (2000 ft) 2015", x="Time in Dry Season", y="Species Richness", fill="Time in Dry Season")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15)) +
  stat_compare_means(method = "anova",label.y=315) +stat_compare_means(comparisons = list(c(1,2)), method="t.test", hide.ns = FALSE,label = "p.signif")
ggsave(f2000.its1.15.2a,filename = "figures/2000ft_fungal_SR_2015_binned_month_8.3.21.pdf", width=8, height=6, dpi=600)

## its1 species richness w/in 2700 ft.

f2700.its1.15.2a<-ggplot(its1.div.15_2700, aes(x=MonthBin, y=ITS1_Species_Richness, fill=MonthBin))+geom_boxplot(color="black")+scale_x_discrete()+theme_bw()+scale_fill_manual(values = saturation(season_col, 1))+theme_classic()+
  labs(title = "Fungal Species Richness (2700 ft) 2015", x="Time in Dry Season", y="Species Richness", fill="Time in Dry Season")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15)) +
  stat_compare_means(method = "anova",label.y=310) +stat_compare_means(comparisons = list(c(1,2)), method="t.test", hide.ns = FALSE,label = "p.signif")
ggsave(f2700.its1.15.2a,filename = "figures/2700ft_fungal_SR_2015_binned_month_8.3.21.pdf", width=8, height=6, dpi=600)

#### 16S Species Richness VIsualization ####

# ## using binary table for species richness visualization
# bac.SR$Month2 <- factor(bac.SR$Month, levels = c("July","August","October")) ## reeorder month column so that it will be listed in chronological order in x axis
# summary(aov(Bac_Species_Richness~Year, data=bac.div.metadat))
# t.test(Bac_Species_Richness~Year, data=bac.div.metadat)
# 
# bac.div.meta_14<-subset(bac.div.metadat, bac.div.metadat$Year=='2014')
# bac.div.meta_15<-subset(bac.div.metadat, bac.div.metadat$Year=='2015')

# species richness by year
bac_SR.yr<-ggplot(bac.div.metadat, aes(x=factor(Year), y=Bac_Species_Richness, fill=factor(Year))) +geom_boxplot(color="black")+scale_x_discrete()+theme_bw()+scale_fill_manual(values = saturation(wes1, 1))+theme_classic()+
  labs(title = "Microbial Species Richness by Sampling Year", x="Year", y="Species Richness", fill="Year")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15)) + 
  stat_compare_means(comparisons = list(c(1,2)), method="t.test", hide.ns = FALSE,label = "p.signif")+stat_compare_means(method = "anova",label.y=1.5) 

ggsave(bac_SR.yr,filename = "figures/Bacterial_SR_by_Year_5.27.21.pdf", width=8, height=6, dpi=600)

# species richness + elevation 2014

bac_SR.14.elev<-ggplot(bac.div.meta_14, aes(x = factor(Elevation), y = Bac_Species_Richness, fill=factor(Elevation))) + 
  geom_boxplot() + theme_classic() + guides(fill = guide_legend(reverse=TRUE)) + 
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11)) +
  labs(title="Microbial Species Richness x Elevation (2014)",fill="Elevation (ft)")+ylab("Species Richness")+xlab("Elevation (ft)")+scale_fill_manual(values=saturation(fair_cols, 0.9))+stat_compare_means(method = "anova",label.y=1.5) +stat_compare_means(comparisons = list(c(1,2), c(2,3),  c(1,3),  c(3,4),  c(2,4),  c(1,4)), method="t.test", hide.ns = FALSE,label = "p.signif")

ggsave(bac_SR.14.elev,filename = "figures/Bac_Spec_Richness_Elev_2014_sigbars_5.27.21.pdf", width=10, height=8, dpi=600)

bac_SR.14.elev2<-ggplot(bac.div.meta_14, aes(x = factor(Elevation), y = Bac_Species_Richness, fill=factor(Elevation))) + 
  geom_boxplot() + theme_classic() + guides(fill = guide_legend(reverse=TRUE)) + 
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11)) +
  labs(title="Microbial Species Richness x Elevation (2014)",fill="Elevation (ft)")+ylab("Species Richness")+xlab("Elevation (ft)")+scale_fill_manual(values=saturation(fair_cols, 0.9))+stat_compare_means(method = "anova",label.y=1.5) +stat_compare_means(comparisons = list(c(1,2), c(2,3),  c(1,3),  c(3,4),  c(2,4),  c(1,4)), method="t.test", label = "p.format")

ggsave(bac_SR.14.elev2,filename = "figures/Bac_Spec_Richness_Elev_2014_pvals_5.27.21.pdf", width=10, height=8, dpi=600)

# species richness + elevation 2015

bac_SR.15.elev<-ggplot(bac.div.meta_15, aes(x = factor(Elevation), y = Bac_Species_Richness, fill=factor(Elevation))) + 
  geom_boxplot() + theme_classic() + guides(fill = guide_legend(reverse=TRUE)) + 
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11)) +
  labs(title="Microbial Species Richness x Elevation (2015)",fill="Elevation (ft)")+ylab("Species Richness")+xlab("Elevation (ft)")+scale_fill_manual(values=saturation(fair_cols, 0.9))+stat_compare_means(method = "anova",label.y=1.5) +stat_compare_means(comparisons = list(c(1,2), c(2,3),  c(1,3),  c(3,4),  c(2,4),  c(1,4)), method="t.test", hide.ns = FALSE,label = "p.signif")

ggsave(bac_SR.15.elev,filename = "figures/Bac_Spec_Richness_Elev_2015_sigbars_5.27.21.pdf", width=10, height=8, dpi=600)

bac_SR.15.elev2<-ggplot(bac.div.meta_15, aes(x = factor(Elevation), y = Bac_Species_Richness, fill=factor(Elevation))) + 
  geom_boxplot() + theme_classic() + guides(fill = guide_legend(reverse=TRUE)) + 
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11)) +
  labs(title="Microbial Species Richness x Elevation (2015)",fill="Elevation (ft)")+ylab("Species Richness")+xlab("Elevation (ft)")+scale_fill_manual(values=saturation(fair_cols, 0.9))+stat_compare_means(method = "anova",label.y=1.5) +stat_compare_means(comparisons = list(c(1,2), c(2,3),  c(1,3),  c(3,4),  c(2,4),  c(1,4)), method="t.test", label = "p.format")

ggsave(bac_SR.15.elev2,filename = "figures/Bac_Spec_Richness_Elev_2015_pvals_5.27.21.pdf", width=10, height=8, dpi=600)

bacsr_yr_1<-ggarrange(bac_SR.14.elev, bac_SR.15.elev, legend="right", common.legend=TRUE,ncol=2, nrow=1)
ggsave(bacsr_yr_1,filename = "figures/Bac_Spec_Richness_Elev_both_yrs_5.27.21.pdf", width=20, height=8, dpi=600)

bacsr_yr_2<-ggarrange(bac_SR.14.elev2, bac_SR.15.elev2, legend="right", common.legend=TRUE, ncol=2, nrow=1)
ggsave(bacsr_yr_2,filename = "figures/Bac_Spec_Richness_Elev_both_yrs_pvals_5.27.21.pdf", width=20, height=8, dpi=600)


bac_SR.yr<-ggplot(bac.SR, aes(x=factor(Year), y=bac_Species_Richness, fill=factor(Year))) +geom_boxplot(color="black")+scale_x_discrete()+theme_bw()+scale_fill_manual(values = saturation(wes1, 1))+theme_classic()+
  labs(title = "Microbial Species Richness by Sampling Year", x="Year", y="Species Richness", fill="Year")+theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15))
ggsave(bac_SR.yr,filename = "figures/Microbial_Species_Richness_by_Year_Sierra_3.29.21.pdf", width=8, height=6, dpi=600)

# species richness by month
bac_SR.month<-ggplot(bac.SR, aes(x=Month2, y=bac_Species_Richness, fill=Month2)) +geom_boxplot(color="black")+scale_x_discrete()+theme_bw()+scale_fill_manual(values = saturation(wes1, 1))+theme_classic()+
  labs(title = "Microbial Species Richness by Sampling Month", x="Month", y="Species Richness", fill="Month")+theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15))
ggsave(bac_SR.month,filename = "figures/Microbial_Species_Richness_by_Month_Sierra_3.29.21_TEST.pdf", width=8, height=6, dpi=600)

# species richness across elevation
bac_SR.elev<-ggplot(bac.SR, aes(x=factor(Elevation), y=bac_Species_Richness, fill=factor(Elevation))) +geom_boxplot(color="black")+scale_x_discrete()+theme_bw()+scale_fill_manual(values = saturation(wes1, 1))+theme_classic()+
  labs(title = "Microbial Species Richness by Elevation", x="Elevation", y="Species Richness", fill="Elevation")+theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15))
ggsave(bac_SR.elev,filename = "figures/Microbial_Species_Richness_by_Elevation_Sierra_3.29.21.pdf", width=8, height=6, dpi=600)

# species richness across elevation + year
bac_SR.elev.yr<-ggplot(bac.SR, aes(x=factor(Elevation), y=bac_Species_Richness, fill=factor(Year)))+geom_boxplot(color="black")+scale_x_discrete()+theme_bw()+scale_fill_manual(values = saturation(wes1, 1))+theme_classic()+
  labs(title = "Microbial Species Richness by Elevation and Year", x="Elevation", y="Species Richness", fill="Year")+theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15))
ggsave(bac_SR.elev.yr,filename = "figures/Microbial_Species_Richness_by_Elevation_x_Year_Sierra_3.29.21.pdf", width=8, height=6, dpi=600)

### bac Species Richness within Elevations
## bac richness w/in 400 ft.
f400.bac.3<-ggplot(bac.div_400, aes(x=factor(Month), y=Bac_Species_Richness, fill=factor(Month)))+geom_boxplot(color="black")+scale_x_discrete()+theme_bw()+scale_fill_manual(values = saturation(Sag_pal, 0.8))+theme_classic()+
  labs(title = "Microbial Species Richness by Month (400 ft)", x="Month", y="Species Richness", fill="Month")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15))+
  stat_compare_means(method = "anova")
ggsave(f400.bac.3,filename = "figures/400ft_bacterial_spec_ric_month_6.9.21.pdf", width=8, height=6, dpi=600)

f400.bac.3a<-ggplot(bac.div_400, aes(x=MonthBin, y=Bac_Species_Richness, fill=MonthBin))+geom_boxplot(color="black")+scale_x_discrete()+theme_bw()+scale_fill_manual(values = saturation(season_col, 1))+theme_classic()+
  labs(title = "Microbial Species Richness (400 ft)", x="Time in Dry Season", y="Species Richness", fill="Time in Dry Season")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15))+
  stat_compare_means(method = "anova",label.y=9250)+stat_compare_means(comparisons = list(c(1,2)), method="t.test", hide.ns = FALSE,label = "p.signif")
ggsave(f400.bac.3a,filename = "figures/400ft_bacterial_spec_ric_binned_month_8.3.21.pdf", width=8, height=6, dpi=600)

f400.bac.4<-ggplot(bac.div_400, aes(x=factor(Year), y=Bac_Species_Richness, fill=factor(Month)))+geom_boxplot(color="black")+scale_x_discrete()+theme_bw()+scale_fill_manual(values = saturation(Sag_pal, 1))+theme_classic()+
  labs(title = "Microbial Species Richness (400 ft) by Year & Month", x="Year", y="Species Richness", fill="Month")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15))+
  stat_compare_means(method = "anova") +stat_compare_means(comparisons = list(c(1,2), c(2,3),  c(1,3),  c(3,4),  c(2,4),  c(1,4), c(1,5), c(2,5), c(3,5), c(4,5)), method="t.test", hide.ns = TRUE,label = "p.signif")
ggsave(f400.bac.4,filename = "figures/400ft_bacterial_spec_ric_yr.mon_6.9.21.pdf", width=8, height=6, dpi=600)

## bac Richness w/in 1100 ft.
f1100.bac.3<-ggplot(bac.div_1100, aes(x=factor(Month), y=Bac_Species_Richness, fill=factor(Month)))+geom_boxplot(color="black")+scale_x_discrete()+theme_bw()+scale_fill_manual(values = saturation(Sag_pal, 0.8))+theme_classic()+
  labs(title = "Microbial Species Richness by Month (1100 ft)", x="Month", y="Species Richness", fill="Month")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15)) +
  stat_compare_means(method = "anova") +stat_compare_means(comparisons = list(c(1,2), c(2,3),  c(1,3),  c(3,4),  c(2,4),  c(1,4), c(1,5), c(2,5), c(3,5), c(4,5)), method="t.test", hide.ns = TRUE,label = "p.signif")
ggsave(f1100.bac.3,filename = "figures/1100ft_bacterial_spec_ric_month_6.9.21.pdf", width=8, height=6, dpi=600)

f1100.bac.3a<-ggplot(bac.div_1100, aes(x=MonthBin, y=Bac_Species_Richness, fill=MonthBin))+geom_boxplot(color="black")+scale_x_discrete()+theme_bw()+scale_fill_manual(values = saturation(season_col, 1))+theme_classic()+
  labs(title = "Microbial Species Richness (1100 ft)", x="Time in Dry Season", y="Species Richness", fill="Time in Dry Season")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15)) +
  stat_compare_means(method = "anova",label.y=7000) +stat_compare_means(comparisons = list(c(1,2)), method="t.test", hide.ns = FALSE,label = "p.signif")
ggsave(f1100.bac.3a,filename = "figures/1100ft_bacterial_spec_ric_binned_month_8.3.21.pdf", width=8, height=6, dpi=600)

f1100.bac.4<-ggplot(bac.div_1100, aes(x=factor(Year), y=Bac_Species_Richness, fill=factor(Month)))+geom_boxplot(color="black")+scale_x_discrete()+theme_bw()+scale_fill_manual(values = saturation(Sag_pal, 1))+theme_classic()+
  labs(title = "Microbial Species Richness (1100 ft) by Year & Month", x="Year", y="Species Richness", fill="Month")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15)) +
  stat_compare_means(method = "anova") +stat_compare_means(comparisons = list(c(1,2), c(2,3),  c(1,3),  c(3,4),  c(2,4),  c(1,4), c(1,5), c(2,5), c(3,5), c(4,5)), method="t.test", hide.ns = TRUE,label = "p.signif")
ggsave(f1100.bac.4,filename = "figures/1100ft_bacterial_spec_ric_yr.mon_6.9.21.pdf", width=8, height=6, dpi=600)

## bac Richness w/in 2000 ft.
f2000.bac.3<-ggplot(bac.div_2000, aes(x=factor(Month), y=Bac_Species_Richness, fill=factor(Month)))+geom_boxplot(color="black")+scale_x_discrete()+theme_bw()+scale_fill_manual(values = saturation(Sag_pal, 0.8))+theme_classic()+
  labs(title = "Microbial Species Richness by Month (2000 ft)", x="Month", y="Species Richness", fill="Month")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15)) +
  stat_compare_means(method = "anova", label.y=7300) +stat_compare_means(comparisons = list(c(1,2), c(2,3),  c(1,3),  c(3,4),  c(2,4),  c(1,4), c(1,5), c(2,5), c(3,5), c(4,5)), method="t.test", hide.ns = TRUE,label = "p.signif")
ggsave(f2000.bac.3,filename = "figures/2000ft_bacterial_spec_ric_month_6.9.21.pdf", width=8, height=6, dpi=600)

f2000.bac.3a<-ggplot(bac.div_2000, aes(x=MonthBin, y=Bac_Species_Richness, fill=MonthBin))+geom_boxplot(color="black")+scale_x_discrete()+theme_bw()+scale_fill_manual(values = saturation(season_col, 1))+theme_classic()+
  labs(title = "Microbial Species Richness (2000 ft)", x="Time in Dry Season", y="Species Richness", fill="Time in Dry Season")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15)) +
  stat_compare_means(method = "anova", label.y=7800) +stat_compare_means(comparisons = list(c(1,2)), method="t.test", hide.ns = FALSE,label = "p.signif")
ggsave(f2000.bac.3a,filename = "figures/2000ft_bacterial_spec_ric_binned_month_8.3.21.pdf", width=8, height=6, dpi=600)

f2000.bac.4<-ggplot(bac.div_2000, aes(x=factor(Year), y=Bac_Species_Richness, fill=factor(Month)))+geom_boxplot(color="black")+scale_x_discrete()+theme_bw()+scale_fill_manual(values = saturation(Sag_pal, 1))+theme_classic()+
  labs(title = "Microbial Species Richness (2000 ft) by Year & Month", x="Year", y="Species Richness", fill="Month")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15)) +
  stat_compare_means(method = "anova", label.y=7300) +stat_compare_means(comparisons = list(c(1,2), c(2,3),  c(1,3),  c(3,4),  c(2,4),  c(1,4), c(1,5), c(2,5), c(3,5), c(4,5)), method="t.test", hide.ns = TRUE,label = "p.signif")
ggsave(f2000.bac.4,filename = "figures/2000ft_bacterial_spec_ric_yr.mon_6.9.21.pdf", width=8, height=6, dpi=600)

## bac Richness w/in 2700 ft.
f2700.bac.3<-ggplot(bac.div_2700, aes(x=factor(Month), y=Bac_Species_Richness, fill=factor(Month)))+geom_boxplot(color="black")+scale_x_discrete()+theme_bw()+scale_fill_manual(values = saturation(Sag_pal, 0.8))+theme_classic()+
  labs(title = "Microbial Species Richness by Month (2700 ft)", x="Month", y="Species Richness", fill="Month")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15)) +
  stat_compare_means(method = "anova") +stat_compare_means(comparisons = list(c(1,2), c(2,3),  c(1,3),  c(3,4),  c(2,4),  c(1,4), c(1,5), c(2,5), c(3,5), c(4,5)), method="t.test", hide.ns = TRUE,label = "p.signif")
ggsave(f2700.bac.3,filename = "figures/2700ft_bacterial_spec_ric_month_6.9.21.pdf", width=8, height=6, dpi=600)

f2700.bac.3a<-ggplot(bac.div_2700, aes(x=MonthBin, y=Bac_Species_Richness, fill=MonthBin))+geom_boxplot(color="black")+scale_x_discrete()+theme_bw()+scale_fill_manual(values = saturation(season_col, 1))+theme_classic()+
  labs(title = "Microbial Species Richness (2700 ft)", x="Time in Dry Season", y="Species Richness", fill="Time in Dry Season")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15)) +
  stat_compare_means(method = "anova",label.y=4100) +stat_compare_means(comparisons = list(c(1,2)), method="t.test", hide.ns = FALSE,label = "p.signif")
ggsave(f2700.bac.3a,filename = "figures/2700ft_bacterial_spec_ric_binned_month_8.3.21.pdf", width=8, height=6, dpi=600)

f2700.bac.4<-ggplot(bac.div_2700, aes(x=factor(Year), y=Bac_Species_Richness, fill=factor(Month)))+geom_boxplot(color="black")+scale_x_discrete()+theme_bw()+scale_fill_manual(values = saturation(Sag_pal, 1))+theme_classic()+
  labs(title = "Microbial Species Richness (2700 ft) by Year & Month", x="Year", y="Species Richness", fill="Month")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15)) +
  stat_compare_means(method = "anova") +stat_compare_means(comparisons = list(c(1,2), c(2,3),  c(1,3),  c(3,4),  c(2,4),  c(1,4), c(1,5), c(2,5), c(3,5), c(4,5)), method="t.test", hide.ns = TRUE,label = "p.signif")
ggsave(f2700.bac.4,filename = "figures/2700ft_bacterial_spec_ric_month.yr_6.9.21.pdf", width=8, height=6, dpi=600)


## 16S Species Richness within Elevations + Years by MonthBin
# sanity check to see if these dataframes are in the global env.
head(bac.div_2014)
head(bac.div.14_400)

## 16S Richness by Binned Month 2014

bac.mbin.2a<-ggplot(bac.div_2014, aes(x=MonthBin, y=Bac_Species_Richness, fill=MonthBin))+geom_boxplot(color="black")+scale_x_discrete()+theme_bw()+scale_fill_manual(values = saturation(season_col, 1))+theme_classic()+
  labs(title = "Microbial Species Richness 2014", x="Time in Dry Season", y="Species Richness", fill="Time in Dry Season")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15))+
  stat_compare_means(method = "anova",label.y=7500)+stat_compare_means(comparisons = list(c(1,2)), method="t.test", hide.ns = FALSE,label = "p.signif")
ggsave(bac.mbin.2a,filename = "figures/microb_SR_2014_binned_month_8.3.21.pdf", width=8, height=6, dpi=600)

## bac species richness w/in 400 ft.

f400.bac.2a<-ggplot(bac.div.14_400, aes(x=MonthBin, y=Bac_Species_Richness, fill=MonthBin))+geom_boxplot(color="black")+scale_x_discrete()+theme_bw()+scale_fill_manual(values = saturation(season_col, 1))+theme_classic()+
  labs(title = "Microbial Species Richness (400 ft) 2014", x="Time in Dry Season", y="Species Richness", fill="Time in Dry Season")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15))+
  stat_compare_means(method = "anova",label.y=5700)+stat_compare_means(comparisons = list(c(1,2)), method="t.test", hide.ns = FALSE,label = "p.signif")
ggsave(f400.bac.2a,filename = "figures/400ft_microb_SR_2014_binned_month_8.3.21.pdf", width=8, height=6, dpi=600)

## bac species richness w/in 1100 ft.

f1100.bac.2a<-ggplot(bac.div.14_1100, aes(x=MonthBin, y=Bac_Species_Richness, fill=MonthBin))+geom_boxplot(color="black")+scale_x_discrete()+theme_bw()+scale_fill_manual(values = saturation(season_col, 1))+theme_classic()+
  labs(title = "Microbial Species Richness (1100 ft) 2014", x="Time in Dry Season", y="Species Richness", fill="Time in Dry Season")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15)) +
  stat_compare_means(method = "anova", label.y=4250) +stat_compare_means(comparisons = list(c(1,2)), method="t.test", hide.ns = FALSE,label = "p.signif")
ggsave(f1100.bac.2a,filename = "figures/1100ft_microb_SR_2014_binned_month_8.3.21.pdf", width=8, height=6, dpi=600)

## bac species richness w/in 2000 ft.

f2000.bac.2a<-ggplot(bac.div.14_2000, aes(x=MonthBin, y=Bac_Species_Richness, fill=MonthBin))+geom_boxplot(color="black")+scale_x_discrete()+theme_bw()+scale_fill_manual(values = saturation(season_col, 1))+theme_classic()+
  labs(title = "Microbial Species Richness (2000 ft) 2014", x="Time in Dry Season", y="Species Richness", fill="Time in Dry Season")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15)) +
  stat_compare_means(method = "anova",label.y=7100) +stat_compare_means(comparisons = list(c(1,2)), method="t.test", hide.ns = FALSE,label = "p.signif")
ggsave(f2000.bac.2a,filename = "figures/2000ft_microb_SR_2014_binned_month_8.3.21.pdf", width=8, height=6, dpi=600)

## bac species richness w/in 2700 ft.

f2700.bac.2a<-ggplot(bac.div.14_2700, aes(x=MonthBin, y=Bac_Species_Richness, fill=MonthBin))+geom_boxplot(color="black")+scale_x_discrete()+theme_bw()+scale_fill_manual(values = saturation(season_col, 1))+theme_classic()+
  labs(title = "Microbial Species Richness (2700 ft) 2014", x="Time in Dry Season", y="Species Richness", fill="Time in Dry Season")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15)) +
  stat_compare_means(method = "anova",label.y=2250) +stat_compare_means(comparisons = list(c(1,2)), method="t.test", hide.ns = FALSE,label = "p.signif")
ggsave(f2700.bac.2a,filename = "figures/2700ft_microb_SR_2014_binned_month_8.3.21.pdf", width=8, height=6, dpi=600)


## 16S Species Richness by Binned Month 2015

bac.mbin.3a<-ggplot(bac.div_2015, aes(x=MonthBin, y=Bac_Species_Richness, fill=MonthBin))+geom_boxplot(color="black")+scale_x_discrete()+theme_bw()+scale_fill_manual(values = saturation(season_col, 1))+theme_classic()+
  labs(title = "Microbial Species Richness 2015", x="Time in Dry Season", y="Species Richness", fill="Time in Dry Season")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15))+
  stat_compare_means(method = "anova",label.y=9100)+stat_compare_means(comparisons = list(c(1,2)), method="t.test", hide.ns = FALSE,label = "p.signif")
ggsave(bac.mbin.3a,filename = "figures/microb_SR_2015_binned_month_8.3.21.pdf", width=8, height=6, dpi=600)

## bac species richness w/in 400 ft.

f400.bac.15.2a<-ggplot(bac.div.15_400, aes(x=MonthBin, y=Bac_Species_Richness, fill=MonthBin))+geom_boxplot(color="black")+scale_x_discrete()+theme_bw()+scale_fill_manual(values = saturation(season_col, 1))+theme_classic()+
  labs(title = "Microbial Species Richness (400 ft) 2015", x="Time in Dry Season", y="Species Richness", fill="Time in Dry Season")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15))+
  stat_compare_means(method = "anova",label.y=9000)+stat_compare_means(comparisons = list(c(1,2)), method="t.test", hide.ns = FALSE,label = "p.signif")
ggsave(f400.bac.15.2a,filename = "figures/400ft_microb_SR_2015_binned_month_8.3.21.pdf", width=8, height=6, dpi=600)

## bac species richness w/in 1100 ft.

f1100.bac.15.2a<-ggplot(bac.div.15_1100, aes(x=MonthBin, y=Bac_Species_Richness, fill=MonthBin))+geom_boxplot(color="black")+scale_x_discrete()+theme_bw()+scale_fill_manual(values = saturation(season_col, 1))+theme_classic()+
  labs(title = "Microbial Species Richness (1100 ft) 2015", x="Time in Dry Season", y="Species Richness", fill="Time in Dry Season")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15)) +
  stat_compare_means(method = "anova", label.y=6700) +stat_compare_means(comparisons = list(c(1,2)), method="t.test", hide.ns = FALSE,label = "p.signif")
ggsave(f1100.bac.15.2a,filename = "figures/1100ft_microb_SR_2015_binned_month_8.3.21.pdf", width=8, height=6, dpi=600)

## bac species richness w/in 2000 ft.

f2000.bac.15.2a<-ggplot(bac.div.15_2000, aes(x=MonthBin, y=Bac_Species_Richness, fill=MonthBin))+geom_boxplot(color="black")+scale_x_discrete()+theme_bw()+scale_fill_manual(values = saturation(season_col, 1))+theme_classic()+
  labs(title = "Microbial Species Richness (2000 ft) 2015", x="Time in Dry Season", y="Species Richness", fill="Time in Dry Season")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15)) +
  stat_compare_means(method = "anova",label.y=5060) +stat_compare_means(comparisons = list(c(1,2)), method="t.test", hide.ns = FALSE,label = "p.signif")
ggsave(f2000.bac.15.2a,filename = "figures/2000ft_microb_SR_2015_binned_month_8.3.21.pdf", width=8, height=6, dpi=600)

## bac species richness w/in 2700 ft.

f2700.bac.15.2a<-ggplot(bac.div.15_2700, aes(x=MonthBin, y=Bac_Species_Richness, fill=MonthBin))+geom_boxplot(color="black")+scale_x_discrete()+theme_bw()+scale_fill_manual(values = saturation(season_col, 1))+theme_classic()+
  labs(title = "Microbial Species Richness (2700 ft) 2015", x="Time in Dry Season", y="Species Richness", fill="Time in Dry Season")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15)) +
  stat_compare_means(method = "anova",label.y=4050) +stat_compare_means(comparisons = list(c(1,2)), method="t.test", hide.ns = FALSE,label = "p.signif")
ggsave(f2700.bac.15.2a,filename = "figures/2700ft_microb_SR_2015_binned_month_8.3.21.pdf", width=8, height=6, dpi=600)
