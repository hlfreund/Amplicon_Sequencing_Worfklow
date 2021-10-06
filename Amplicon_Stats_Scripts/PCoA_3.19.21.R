#if (!requireNamespace("BiocManager", quietly = TRUE))
#install.packages("BiocManager")
#BiocManager::install(version = "3.10")
#BiocManager::install(c("phyloseq"))
library(phyloseq)
library(devtools)
install_github("vqv/ggbiplot")
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

## binary 16S table

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

#### Import metadata ####
metadata<-as.data.frame(read_excel("data/SierraMapEle2SH.xls"))

# Import dust metadata
dust_comp<-read.csv("data/BinnedHiMidLow.csv", header=TRUE, sep=",")
head(dust_comp)
dust_comp<-subset(dust_comp, select=c(SampleID, DustComplexity, LoMidHi)) # subset only part of the metadata we need
rownames(dust_comp)<-dust_comp$SampleID
metadata<-merge(metadata, dust_comp, by.x="SampleID", by.y="SampleID")
head(metadata)
rownames(metadata)<-metadata$SampleID

tail(metadata)
metadata<-subset(metadata, select=-c(DateCode, BarcodeSequence, LinkerPrimerSequence, SiteCode, RepNum, SiteRep, TubeID, DateCode, Description)) # subset only part of the metadata we need
head(metadata)

# some sanity checks
dim(metadata)
dim(its1_otu_table)
dim(bac_otu_table)

row.names(metadata)
row.names(its1_otu_table)
row.names(bac_otu_table)

## Reorder rows of metadata and all dfs
metadata=metadata[rownames(its1_otu_table),] ## reorder metadata to have same rows as original OTU table
# ^ this indexing method will only work if the two dfs have the same # of rows AND the same row names!
bac_otu_table=bac_otu_table[rownames(metadata),] ## reorder bacterial OTU table to have same rows as ITS1 OTU table + metadata

bac_binary_table=bac_binary_table[rownames(metadata),] ## reorder binary tables to have same rows as original OTU tables/metadata
its1_binary_table=its1_binary_table[rownames(metadata),] ## reorder binary tables to have same rows as original OTU tables/metadata


#### Subset metadata ####
meta_category<-subset(metadata, select=c(SampleID, Year, Month, Site, Elevation, LoMidHi)) # subset only part of the metadata we need
meta_quant<-subset(metadata, select=-c(SampleID, Year, NdEp, A_87Sr86Sr, Month, Site, SampType, Elevation, DustComplexity, LoMidHi)) # subset only part of the metadata we need

meta_quant<-sapply(meta_quant, as.numeric) # turns all columns in df into numeric data (for pcoa)
tail(meta_quant)
class(meta_quant)
meta_quant<-as.data.frame(meta_quant); rownames(meta_quant)<-rownames(metadata) # can only do this if you're sure the dfs are in the same exact order
head(meta_quant)

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

# scale continuous (numeric) data
meta_quant.scale<-as.data.frame(scale(meta_quant))
head(meta_quant.scale)
meta_quant.scale$SampleID<-rownames(meta_quant.scale)

site_dust<-data.frame(SampleID=metadata$SampleID, DustComplex=metadata$DustComplexity, Site=metadata$Site)
meta_qdat<-merge(meta_quant.scale, site_dust, by.x="SampleID", by.y="SampleID")
head(meta_qdat)

meta_all_scaled<-merge(meta_qdat[,-11], meta_category, by.x="SampleID", by.y="SampleID")
head(meta_all_scaled)

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

fair_cols <- c("#38170B","#BF1B0B", "#FFC465", "#66ADE5", "#252A52")
names(fair_cols) <- letters[1:5]
fair_ramp <- scales::colour_ramp(fair_cols)
fair_sat <- saturation(fair_ramp, 1)

# custom color palette
head(metadata)

warm2cold1<-get_palette(paste0("#", c("720026", "D14D60", "0077B6","03045E")),k=4)
names(warm2cold1) <- levels(metadata$Elevation)

warm2cold2<-get_palette(paste0("#", c("720026", "D14D60", "0077B6","03045E")),k=4)
names(warm2cold2) <- levels(metadata$Site)

warm2cold2<-get_palette(paste0("#", c("720026", "D14D60", "0077B6","03045E")),k=4)
names(warm2cold2) <- levels(metadata$Site)

#### Bray-Curtis Distance Matrices w/ relativized data ####
# *** using relativized data for Bray-Curtis distance matrices
## ITS1 first
class(its1_otu_table) # raw counts
dim(its1_otu_table)
row.names(its1_otu_table)

its1_RA_table<-data.frame(decostand(its1_otu_table, method="total", MARGIN=1, na.rm=TRUE))  
# relative abundance of taxa data where everything is divided by margin total (default MARGIN = 1 = rows) -- rows = samples
rowSums(its1_RA_table)

its1_RA_bray<-vegdist(its1_RA_table,method="bray")     ####### distance matrix with RELATIVIZED data!!! 
# **** in vegan: ROWS need to be SITES/samples; COLUMNS are SPECIES (OTUs, ASVs) -- vegdist!
class(its1_RA_bray) # needs to be class dist for pcoa()

## 16S next
class(bac_otu_table) # raw counts
dim(bac_otu_table)
row.names(bac_otu_table)

bac_RA_table<-data.frame(decostand(bac_otu_table, method="total", MARGIN=1, na.rm=TRUE))  
# relative abundance of taxa data where everything is divided by margin total (default MARGIN = 1 = rows) -- rows = samples
rowSums(bac_RA_table)
bac_RA_bray<-vegdist(bac_RA_table,method="bray")     ####### distance matrix with RELATIVIZED data!!! 
# **** in vegan: ROWS need to be SITES/samples; COLUMNS are SPECIES (OTUs, ASVs) -- vegdist!
class(bac_RA_bray) # needs to be class dist for pcoa()

#### PCoA w/ relativized, BrayCurtis matrix ####
## fungi first
its1_pcoa1 = pcoa(its1_RA_bray) # pcoa of Bray-Curtis dissimilarity matrix made w/ relativized abundance (site x species) data
# The proportion of variances explained is in its element values$Relative_eig
#cmdscale(bac.bray, k=(nrow(bac_otu_table)-1), eig=TRUE)

str(its1_pcoa1)
biplot(its1_pcoa1)
## project (scaled) chem data on the PCoA ordination
biplot(its1_pcoa1, meta_qdat[,2:9], main="ITS1 PCoA + Scaled Chemical Data", xlab = "Axis 1", ylab = "Axis 2")
biplot(its1_pcoa1, meta_qdat[,2:9],display = c("sites", "species"),type = c("text","points"), main="ITS1 PCoA + Scaled Chemical Data", xlab = "Axis 1", ylab = "Axis 2")

## bacteria + archaea next

bac_pcoa1 = pcoa(bac_RA_bray) # pcoa of Bray-Curtis dissimilarity matrix made w/ relativized abundance (site x species) data
# The proportion of variances explained is in its element values$Relative_eig
#cmdscale(bac.bray, k=(nrow(bac_otu_table)-1), eig=TRUE)

biplot(bac_pcoa1)
## project (scaled) chem data on the PCoA ordination
biplot(bac_pcoa1, meta_qdat[,2:9], main="16S PCoA + Scaled Chemical Data", xlab = "Axis 1", ylab = "Axis 2")
biplot(bac_pcoa1, meta_qdat[,2:9],display = c("sites", "species"),type = c("text","points"))

#### Function for custom biplot ####
PCbiplot <- function(PC, x="PC1", y="PC2") {
  # PC being a prcomp object
  data <- data.frame(obsnames=row.names(PC$x), PC$x)
  plot <- ggplot(data, aes_string(x=x, y=y)) + geom_text(alpha=.4, size=3, aes(label=obsnames))
  plot <- plot + geom_hline(aes(0), size=.2) + geom_vline(aes(0), size=.2)
  datapc <- data.frame(varnames=rownames(PC$rotation), PC$rotation)
  mult <- min(
    (max(data[,y]) - min(data[,y])/(max(datapc[,y])-min(datapc[,y]))),
    (max(data[,x]) - min(data[,x])/(max(datapc[,x])-min(datapc[,x])))
  )
  datapc <- transform(datapc,
                      v1 = .7 * mult * (get(x)),
                      v2 = .7 * mult * (get(y))
  )
  plot <- plot + coord_equal() + geom_text(data=datapc, aes(x=v1, y=v2, label=varnames), size = 5, vjust=1, color="red")
  plot <- plot + geom_segment(data=datapc, aes(x=0, y=0, xend=v1, yend=v2), arrow=arrow(length=unit(0.2,"cm")), alpha=0.75, color="red")
  plot
}

fit <- prcomp(USArrests, scale=T)
PCbiplot(fit)
#### Fungal PCoA Visualization (relativized Bray-Curtis counts) ####

### fungi first

#tiff(file="figures/its1_pcoa_biplot.tiff",width=10, height=10, units="in", res=600)
biplot(its1_pcoa1, scale(meta_quant))
dev.off()

# setting up the dataframe for ggplot2 visualization

its1_pcoa1.vectors<-data.frame(its1_pcoa1$vectors)
its1_pcoa1.vectors$SampleID<-rownames(its1_pcoa1.vectors)

its1_pcoa1_meta<-merge(its1_pcoa1.vectors, meta_all_scaled, by.x="SampleID", by.y="SampleID")
its1_pcoa1_df=data.frame(SampleID = its1_pcoa1_meta$SampleID, Axis1 = its1_pcoa1_meta$Axis.1, Axis2 = its1_pcoa1_meta$Axis.2, Site = its1_pcoa1_meta$Site, Elevation = its1_pcoa1_meta$Elevation, Month=its1_pcoa1_meta$Month, Year= its1_pcoa1_meta$Year, DustComp=its1_pcoa1_meta$DustComplex, DustBin=its1_pcoa1_meta$LoMidHi) #create dataframe with PCoA axes and some metadata

its1_pcoa1_df$Month2 <- factor(its1_pcoa1_df$Month, levels = c("July","August","October")) 
its1_pcoa1_df$Elevation2 <- factor(its1_pcoa1_df$Elevation, levels = c("400","1100","2000","2700")) 
its1_pcoa1_df$Site2 <- factor(its1_pcoa1_df$Site, levels = c("SJER","SOAPROOT","PROVIDENCE","SHORTHAIR")) 

head(its1_pcoa1_df)

# plot time

its1.pcoa.fig.1<-ggplot(its1_pcoa1_df, aes(x=Axis1, y=Axis2, col=Month2, shape=factor(Year))) +geom_point(alpha=0.5,size=5)+theme_bw()+scale_colour_manual(values =  saturation(SM_pal, 1))+
  coord_cartesian(xlim=c(-0.5,0.5),ylim=c(-0.5,0.5))+labs(xlab="Axis 1", ylab="Axis 2",color="Month", shape="Year")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))

ggsave(its1.pcoa.fig.1,filename = "figures/fungal_pcoa1_BC_Sierra_4.18.21.pdf", width=8, height=6, dpi=600)

its1.pcoa.fig.2<-ggplot(its1_pcoa1_df, aes(x=Axis1, y=Axis2, shape=factor(Year), col=factor(Elevation))) +geom_point(alpha=0.5, size=5)+theme_bw()+scale_colour_manual(values = saturation(warm2cold1, 2))+
  coord_cartesian(xlim=c(-0.5,0.5),ylim=c(-0.5,0.5))+labs(xlab="Axis 1", ylab="Axis 2",color="Elevation", shape="Year",size="Elevation")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(color = guide_legend(override.aes = list(size = 5)),shape = guide_legend(override.aes = list(size = 5)))

ggsave(its1.pcoa.fig.2,filename = "figures/fungal_pcoa2_BC_Sierra_4.18.21.pdf", width=8, height=6, dpi=600)

its1.pcoa.fig.3<-ggplot(its1_pcoa1_df, aes(x=Axis1, y=Axis2, col=Site2, shape=Month2)) +geom_point(alpha=0.5,size=5)+theme_bw()+scale_colour_manual(values =  saturation(warm2cold2, 1))+
  coord_cartesian(xlim=c(-0.5,0.5),ylim=c(-0.5,0.5))+labs(xlab="Axis 1", ylab="Axis 2",color="Site", shape="Month")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))

ggsave(its1.pcoa.fig.3,filename = "figures/fungal_pcoa3_BC_Sierra_4.18.21.pdf", width=8, height=6, dpi=600)

its1.pcoa.fig.4<-ggplot(its1_pcoa1_df, aes(x=Axis1, y=Axis2, col=Site2, shape=Elevation2)) +geom_point(alpha=0.5,size=5)+theme_bw()+scale_colour_manual(values =  saturation(warm2cold2, 1))+
  coord_cartesian(xlim=c(-0.5,0.5),ylim=c(-0.5,0.5))+labs(xlab="Axis 1", ylab="Axis 2",color="Site", shape="Elevation")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))

ggsave(its1.pcoa.fig.4,filename = "figures/fungal_pcoa4_BC_Sierra_4.18.21.pdf", width=8, height=6, dpi=600)

its1.pcoa.fig.5a<-ggplot(its1_pcoa1_df, aes(x=Axis1, y=Axis2, col=Elevation2, shape=Month2)) +geom_point(alpha=0.5,size=5)+theme_bw()+scale_colour_manual(values =  saturation(warm2cold1, 1))+
  coord_cartesian(xlim=c(-0.5,0.5),ylim=c(-0.5,0.5))+labs(xlab="Axis 1", ylab="Axis 2",color="Elevation", shape="Month")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))

ggsave(its1.pcoa.fig.5a,filename = "figures/fungal_pcoa5a_BC_Sierra_4.18.21.pdf", width=8, height=6, dpi=600)

# Dust Copmlexity PCoAs

f.dust.pcoa1<-ggplot(its1_pcoa1_df, aes(x=Axis1, y=Axis2, col=DustComp, shape=factor(Year))) +geom_point(alpha=0.5,size=5)+theme_bw()+scale_colour_gradientn(colours=fair_cols)+
  coord_cartesian(xlim=c(-0.5,0.5),ylim=c(-0.5,0.5))+labs(xlab="Axis 1", ylab="Axis 2",color="Dust Complexity", shape="Year")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))
ggsave(f.dust.pcoa1,filename = "figures/Sierra_ITS1_Dust_PCOA_5.3.21.pdf", width=8, height=6, dpi=600)

f.dust.pcoa2<-ggplot(its1_pcoa1_df, aes(x=Axis1, y=Axis2, col=DustComp, shape=factor(Elevation))) +geom_point(alpha=0.5,size=5)+theme_bw()+scale_shape_manual(values = c(18, 16, 17,15))+scale_colour_gradientn(colours=fair_cols)+
  coord_cartesian(xlim=c(-0.5,0.5),ylim=c(-0.5,0.5))+labs(xlab="Axis 1", ylab="Axis 2",color="Dust Complexity", shape="Elevation (ft)")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))
ggsave(f.dust.pcoa2,filename = "figures/Sierra_ITS1_Dust_PCOA2_5.3.21.pdf", width=8, height=6, dpi=600)

f.dust.pcoa3<-ggplot(its1_pcoa1_df, aes(x=Axis1, y=Axis2, col=DustComp, shape=factor(interaction(Month,Year)))) +geom_point(alpha=0.5,size=5)+theme_bw()+scale_shape_manual(name="Month & Year", labels=c("August 2014","July 2014","October 2014","July 2015","October 2015"),values = c(16, 17, 15, 2, 0))+saturation(scale_colour_gradientn(colours=fair_cols), 0.9)+
  coord_cartesian(xlim=c(-0.5,0.5),ylim=c(-0.5,0.5))+labs(xlab="Axis 1", ylab="Axis 2",color="Dust Complexity", shape="Month & Year")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))
ggsave(f.dust.pcoa3,filename = "figures/Sierra_ITS1_Dust_PCOA3_5.3.21.pdf", width=8, height=6, dpi=600)

f.dust.pcoa4<-ggplot(its1_pcoa1_df, aes(x=Axis1, y=Axis2, col=DustComp, shape=Site2)) +geom_point(alpha=0.5,size=5)+theme_bw()+scale_shape_manual(values = c(15, 16, 17, 18))+saturation(scale_colour_gradientn(colours=fair_cols), 0.9)+
  coord_cartesian(xlim=c(-0.5,0.5),ylim=c(-0.5,0.5))+labs(xlab="Axis 1", ylab="Axis 2",color="Dust Complexity", shape="Site")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))
ggsave(f.dust.pcoa4,filename = "figures/Sierra_ITS1_Dust_PCOA4_5.3.21.pdf", width=8, height=6, dpi=600)

f.dust.pcoa5<-ggplot(its1_pcoa1_df, aes(x=Axis1, y=Axis2, col=DustComp, shape=factor(Year), size=factor(Elevation))) +geom_point(alpha=0.5)+theme_bw()+scale_shape_manual(values = c(15, 16))+saturation(scale_colour_gradientn(colours=fair_cols), 0.9)+
  coord_cartesian(xlim=c(-0.5,0.5),ylim=c(-0.5,0.5))+scale_size_manual(values=c(3,5,7,9))+labs(xlab="Axis 1", ylab="Axis 2",color="Dust Complexity", size = "Elevation (ft)",shape="Year")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
    guides(shape = guide_legend(override.aes = list(size = 5)))
ggsave(f.dust.pcoa5,filename = "figures/Sierra_ITS1_Dust_PCOA5_5.3.21.pdf", width=8, height=6, dpi=600)

f.dust.pcoa5a<-ggplot(its1_pcoa1_df, aes(x=Axis1, y=Axis2, col=DustComp, shape=factor(Year), size=factor(Elevation))) +geom_point(alpha=0.5)+theme_bw()+scale_shape_manual(values = c(15, 1))+saturation(scale_colour_gradientn(colours=fair_cols), 0.9)+
  coord_cartesian(xlim=c(-0.5,0.5),ylim=c(-0.5,0.5))+scale_size_manual(values=c(3,5,7,9))+labs(xlab="Axis 1", ylab="Axis 2",color="Dust Complexity", size = "Elevation (ft)",shape="Year")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))
ggsave(f.dust.pcoa5a,filename = "figures/Sierra_ITS1_Dust_PCOA5a_5.3.21.pdf", width=8, height=6, dpi=600)

f.dust.pcoa6<-ggplot(its1_pcoa1_df, aes(x=Axis1, y=Axis2, col=DustComp, shape=interaction(Site,Year))) +geom_point(alpha=0.5,size=5)+theme_bw()+scale_shape_manual(name="Site & Year", labels=c("Providence 2014","Shorthair 2014","SJER 2014","Soaproot 2014","Providence 2015","Shorthair 2015","SJER 2015","Soaproot 2015"),values = c(15, 16, 17, 18,0,1,2,5))+saturation(scale_colour_gradientn(colours=fair_cols), 0.9)+
  coord_cartesian(xlim=c(-0.5,0.5),ylim=c(-0.5,0.5))+labs(xlab="Axis 1", ylab="Axis 2",color="Dust Complexity", shape="Site & Year")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))
ggsave(f.dust.pcoa6,filename = "figures/Sierra_ITS1_Dust_PCOA6_5.3.21.pdf", width=8, height=6, dpi=600)

f.dust.pcoa7<-ggplot(its1_pcoa1_df, aes(x=Axis1, y=Axis2, col=Elevation, shape=DustBin)) +geom_point(alpha=0.5,size=5)+theme_bw()+scale_shape_manual(name="Dust Complexity", labels=c("Low (< 0.5)","Mid (0.5-0.65)","High (> 0.65)"),values = c(15, 16, 17))+saturation(scale_colour_gradientn(colours=fair_cols), 0.9)+
  coord_cartesian(xlim=c(-0.5,0.5),ylim=c(-0.5,0.5))+labs(xlab="Axis 1", ylab="Axis 2",color="Elevation (ft)", shape="Dust Complexity")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))
ggsave(f.dust.pcoa7,filename = "figures/Sierra_ITS1_Dust_PCOA7_5.4.21.pdf", width=8, height=6, dpi=600)

f.dust.pcoa8<-ggplot(its1_pcoa1_df, aes(x=Axis1, y=Axis2, col=Elevation, shape=interaction(DustBin,Year))) +geom_point(alpha=0.5,size=5)+theme_bw()+scale_shape_manual(name="Dust Complexity & Year", labels=c("Low 2014","Mid 2014","High 2014","Low 2015","Mid 2015","High 2015"),values = c(15, 16, 17, 0, 1, 2))+saturation(scale_colour_gradientn(colours=fair_cols), 0.9)+
  coord_cartesian(xlim=c(-0.5,0.5),ylim=c(-0.5,0.5))+labs(xlab="Axis 1", ylab="Axis 2",color="Elevation (ft)", shape="Dust Complexity & Year")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))
ggsave(f.dust.pcoa8,filename = "figures/Sierra_ITS1_Dust_PCOA8_5.4.21.pdf", width=8, height=6, dpi=600)

# its1.pcoa.fig.5b<-ggplot(its1_pcoa1_df, aes(x=Axis1, y=Axis2, col=Month2, shape=Elevation2)) +geom_point(alpha=0.5,size=5)+theme_bw()+scale_colour_manual(values =  saturation(Sag_pal, 1))+
#   coord_cartesian(xlim=c(-0.5,0.5),ylim=c(-0.5,0.5))+labs(xlab="Axis 1", ylab="Axis 2",color="Month", shape="Elevation")+theme_classic()+
#   theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
#   guides(shape = guide_legend(override.aes = list(size = 5)))
# 
# ggsave(its1.pcoa.fig.5b,filename = "figures/fungal_pcoa5b_BC_Sierra_4.18.21.pdf", width=8, height=6, dpi=600)


#### Microbial PCoA Visualization (relativized Bray-Curtis counts) ####


#tiff(file="figures/bac_pcoa_biplot.tiff",width=10, height=10, units="in", res=600)
biplot(bac_pcoa1, scale(meta_quant))
dev.off()

# setting up the dataframe for ggplot2 visualization

bac_pcoa1.vectors<-data.frame(bac_pcoa1$vectors)
bac_pcoa1.vectors$SampleID<-rownames(bac_pcoa1.vectors)

bac_pcoa1_meta<-merge(bac_pcoa1.vectors, meta_all_scaled, by="SampleID")
bac_pcoa1_df=data.frame(SampleID = bac_pcoa1_meta$SampleID, Axis1 = bac_pcoa1_meta$Axis.1, Axis2 = bac_pcoa1_meta$Axis.2, Site = bac_pcoa1_meta$Site,Elevation = bac_pcoa1_meta$Elevation, Month=bac_pcoa1_meta$Month, Year= bac_pcoa1_meta$Year, DustComp=bac_pcoa1_meta$DustComplex, DustBin=its1_pcoa1_meta$LoMidHi) #create dataframe with PCoA axes and some metadata

bac_pcoa1_df$Month2 <- factor(bac_pcoa1_df$Month, levels = c("July","August","October")) 
bac_pcoa1_df$Elevation2 <- factor(bac_pcoa1_df$Elevation, levels = c("400","1100","2000","2700")) 
bac_pcoa1_df$Site2 <- factor(bac_pcoa1_df$Site, levels = c("SJER","SOAPROOT","PROVIDENCE","SHORTHAIR")) 

head(bac_pcoa1_df)

bac.pcoa.fig.1<-ggplot(bac_pcoa1_df, aes(x=Axis1, y=Axis2, col=Month2, shape=factor(Year))) +geom_point(alpha=0.5,size=5)+theme_bw()+scale_colour_manual(values =  saturation(SM_pal, 1))+
  coord_cartesian(xlim=c(-0.5,0.5),ylim=c(-0.5,0.5))+labs(xlab="Axis 1", ylab="Axis 2",color="Month", shape="Year")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))

ggsave(bac.pcoa.fig.1,filename = "figures/microbial_pcoa1_BC_Sierra_4.18.21.pdf", width=8, height=6, dpi=600)

bac.pcoa.fig.2<-ggplot(bac_pcoa1_df, aes(x=Axis1, y=Axis2, col=Elevation2, shape=factor(Year))) +geom_point(alpha=0.5,size=5)+theme_bw()+scale_colour_manual(values =  saturation(warm2cold1, 1))+
  coord_cartesian(xlim=c(-0.5,0.5),ylim=c(-0.5,0.5))+labs(xlab="Axis 1", ylab="Axis 2",color="Month", shape="Year",size="Elevation")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(color = guide_legend(override.aes = list(size = 5)),shape = guide_legend(override.aes = list(size = 5)))

ggsave(bac.pcoa.fig.2,filename = "figures/microbial_pcoa2_BC_Sierra_4.18.21.pdf", width=8, height=6, dpi=600)

bac.pcoa.fig.3<-ggplot(bac_pcoa1_df, aes(x=Axis1, y=Axis2, col=Site2, shape=Month2)) +geom_point(alpha=0.5,size=5)+theme_bw()+scale_colour_manual(values =  saturation(warm2cold2, 1))+
  coord_cartesian(xlim=c(-0.5,0.5),ylim=c(-0.5,0.5))+labs(xlab="Axis 1", ylab="Axis 2",color="Site", shape="Month")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))

ggsave(bac.pcoa.fig.3,filename = "figures/microbial_pcoa3_BC_Sierra_4.18.21.pdf", width=8, height=6, dpi=600)

bac.pcoa.fig.4<-ggplot(bac_pcoa1_df, aes(x=Axis1, y=Axis2, col=Site2, shape=Elevation2)) +geom_point(alpha=0.5,size=5)+theme_bw()+scale_colour_manual(values =  saturation(warm2cold2, 1))+
  coord_cartesian(xlim=c(-0.5,0.5),ylim=c(-0.5,0.5))+labs(xlab="Axis 1", ylab="Axis 2",color="Site", shape="Elevation")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))

ggsave(bac.pcoa.fig.4,filename = "figures/microbial_pcoa4_BC_Sierra_4.18.21.pdf", width=8, height=6, dpi=600)

bac.pcoa.fig.5a<-ggplot(bac_pcoa1_df, aes(x=Axis1, y=Axis2, col=Elevation2, shape=Month2)) +geom_point(alpha=0.5,size=5)+theme_bw()+scale_colour_manual(values =  saturation(warm2cold1, 1))+
  coord_cartesian(xlim=c(-0.5,0.5),ylim=c(-0.5,0.5))+labs(xlab="Axis 1", ylab="Axis 2",color="Elevation", shape="Month")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))

ggsave(bac.pcoa.fig.5a,filename = "figures/microbial_pcoa5a_BC_Sierra_4.18.21.pdf", width=8, height=6, dpi=600)

# Dust Copmlexity PCoAs

b.dust.pcoa1<-ggplot(bac_pcoa1_df, aes(x=Axis1, y=Axis2, col=DustComp, shape=factor(Year))) +geom_point(alpha=0.5,size=5)+theme_bw()+scale_colour_gradientn(colours=fair_cols)+
  coord_cartesian(xlim=c(-0.5,0.5),ylim=c(-0.5,0.5))+labs(xlab="Axis 1", ylab="Axis 2",color="Dust Complexity", shape="Year")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))
ggsave(b.dust.pcoa1,filename = "figures/Sierra_16S_Dust_PCOA_5.3.21.pdf", width=8, height=6, dpi=600)

b.dust.pcoa2<-ggplot(bac_pcoa1_df, aes(x=Axis1, y=Axis2, col=DustComp, shape=factor(Elevation))) +geom_point(alpha=0.5,size=5)+theme_bw()+scale_shape_manual(values = c(18, 16, 17,15))+scale_colour_gradientn(colours=fair_cols)+
  coord_cartesian(xlim=c(-0.5,0.5),ylim=c(-0.5,0.5))+labs(xlab="Axis 1", ylab="Axis 2",color="Dust Complexity", shape="Elevation (ft)")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))
ggsave(b.dust.pcoa2,filename = "figures/Sierra_16S_Dust_PCOA2_5.3.21.pdf", width=8, height=6, dpi=600)

b.dust.pcoa3<-ggplot(bac_pcoa1_df, aes(x=Axis1, y=Axis2, col=DustComp, shape=factor(interaction(Month,Year)))) +geom_point(alpha=0.5,size=5)+theme_bw()+scale_shape_manual(name="Month & Year", labels=c("August 2014","July 2014","October 2014","July 2015","October 2015"),values = c(16, 17, 15, 2, 0))+saturation(scale_colour_gradientn(colours=fair_cols), 0.9)+
  coord_cartesian(xlim=c(-0.5,0.5),ylim=c(-0.5,0.5))+labs(xlab="Axis 1", ylab="Axis 2",color="Dust Complexity", shape="Month & Year")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))
ggsave(b.dust.pcoa3,filename = "figures/Sierra_16S_Dust_PCOA3_5.3.21.pdf", width=8, height=6, dpi=600)

b.dust.pcoa4<-ggplot(bac_pcoa1_df, aes(x=Axis1, y=Axis2, col=DustComp, shape=Site2)) +geom_point(alpha=0.5,size=5)+theme_bw()+scale_shape_manual(values = c(15, 16, 17, 18))+saturation(scale_colour_gradientn(colours=fair_cols), 0.9)+
  coord_cartesian(xlim=c(-0.5,0.5),ylim=c(-0.5,0.5))+labs(xlab="Axis 1", ylab="Axis 2",color="Dust Complexity", shape="Site")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))
ggsave(b.dust.pcoa4,filename = "figures/Sierra_16S_Dust_PCOA4_5.3.21.pdf", width=8, height=6, dpi=600)

b.dust.pcoa5<-ggplot(bac_pcoa1_df, aes(x=Axis1, y=Axis2, col=DustComp, shape=factor(Year), size=factor(Elevation))) +geom_point(alpha=0.5)+theme_bw()+scale_shape_manual(values = c(15, 16))+saturation(scale_colour_gradientn(colours=fair_cols), 0.9)+
  coord_cartesian(xlim=c(-0.5,0.5),ylim=c(-0.5,0.5))+scale_size_manual(values=c(3,5,7,9))+labs(xlab="Axis 1", ylab="Axis 2",color="Dust Complexity", size = "Elevation (ft)",shape="Year")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))
ggsave(b.dust.pcoa5,filename = "figures/Sierra_16S_Dust_PCOA5_5.3.21.pdf", width=8, height=6, dpi=600)

b.dust.pcoa5a<-ggplot(bac_pcoa1_df, aes(x=Axis1, y=Axis2, col=DustComp, shape=factor(Year), size=factor(Elevation))) +geom_point(alpha=0.5)+theme_bw()+scale_shape_manual(values = c(15, 1))+saturation(scale_colour_gradientn(colours=fair_cols), 0.9)+
  coord_cartesian(xlim=c(-0.5,0.5),ylim=c(-0.5,0.5))+scale_size_manual(values=c(3,5,7,9))+labs(xlab="Axis 1", ylab="Axis 2",color="Dust Complexity", size = "Elevation (ft)",shape="Year")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))
ggsave(b.dust.pcoa5a,filename = "figures/Sierra_16S_Dust_PCOA5a_5.3.21.pdf", width=8, height=6, dpi=600)

b.dust.pcoa6<-ggplot(bac_pcoa1_df, aes(x=Axis1, y=Axis2, col=DustComp, shape=interaction(Site,Year))) +geom_point(alpha=0.5,size=5)+theme_bw()+scale_shape_manual(name="Site & Year", labels=c("Providence 2014","Shorthair 2014","SJER 2014","Soaproot 2014","Providence 2015","Shorthair 2015","SJER 2015","Soaproot 2015"),values = c(15, 16, 17, 18,0,1,2,5))+saturation(scale_colour_gradientn(colours=fair_cols), 0.9)+
  coord_cartesian(xlim=c(-0.5,0.5),ylim=c(-0.5,0.5))+labs(xlab="Axis 1", ylab="Axis 2",color="Dust Complexity", shape="Site & Year")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))
ggsave(b.dust.pcoa6,filename = "figures/Sierra_16S_Dust_PCOA6_5.3.21.pdf", width=8, height=6, dpi=600)

b.dust.pcoa7<-ggplot(bac_pcoa1_df, aes(x=Axis1, y=Axis2, col=Elevation, shape=DustBin)) +geom_point(alpha=0.5,size=5)+theme_bw()+scale_shape_manual(name="Dust Complexity", labels=c("Low (< 0.5)","Mid (0.5-0.65)","High (> 0.65)"),values = c(15, 16, 17))+saturation(scale_colour_gradientn(colours=fair_cols), 0.9)+
  coord_cartesian(xlim=c(-0.5,0.5),ylim=c(-0.5,0.5))+labs(xlab="Axis 1", ylab="Axis 2",color="Elevation (ft)", shape="Dust Complexity")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))
ggsave(b.dust.pcoa7,filename = "figures/Sierra_16S_Dust_PCOA7_5.4.21.pdf", width=8, height=6, dpi=600)

b.dust.pcoa8<-ggplot(bac_pcoa1_df, aes(x=Axis1, y=Axis2, col=Elevation, shape=interaction(DustBin,Year))) +geom_point(alpha=0.5,size=5)+theme_bw()+scale_shape_manual(name="Dust Complexity & Year", labels=c("Low 2014","Mid 2014","High 2014","Low 2015","Mid 2015","High 2015"),values = c(15, 16, 17, 0, 1, 2))+saturation(scale_colour_gradientn(colours=fair_cols), 0.9)+
  coord_cartesian(xlim=c(-0.5,0.5),ylim=c(-0.5,0.5))+labs(xlab="Axis 1", ylab="Axis 2",color="Elevation (ft)", shape="Dust Complexity & Year")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))
ggsave(b.dust.pcoa8,filename = "figures/Sierra_16S_Dust_PCOA8_5.4.21.pdf", width=8, height=6, dpi=600)

# bac.pcoa.fig.5b<-ggplot(bac_pcoa1_df, aes(x=Axis1, y=Axis2, col=Month2, shape=Elevation2)) +geom_point(alpha=0.5,size=5)+theme_bw()+scale_colour_manual(values =  saturation(Sag_pal, 1))+
#   coord_cartesian(xlim=c(-0.5,0.5),ylim=c(-0.5,0.5))+labs(xlab="Axis 1", ylab="Axis 2",color="Month", shape="Elevation")+theme_classic()+
#   theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
#   guides(shape = guide_legend(override.aes = list(size = 5)))
# 
# ggsave(bac.pcoa.fig.5b,filename = "figures/microbial_pcoa5b_BC_Sierra_4.18.21.pdf", width=8, height=6, dpi=600)



#### Jaccard Distance Matrices ####
## using binary tables for these Jaccard

## fungi first
class(its1_binary_table) # raw counts
dim(its1_binary_table)

its1_jac<-vegdist(its1_binary_table,method="jaccard")     ####### distance matrix with RELATIVIZED data!!! 
# **** in vegan: ROWS need to be SITES/samples; COLUMNS are SPECIES (OTUs, ASVs) -- vegdist!
class(its1_jac) # needs to be class dist for pcoa()

## bacterial/archaeal next

class(bac_binary_table) # raw counts
dim(bac_binary_table)

bac_jac<-vegdist(bac_binary_table,method="jaccard")     ####### distance matrix with RELATIVIZED data!!! 
# **** in vegan: ROWS need to be SITES/samples; COLUMNS are SPECIES (OTUs, ASVs) -- vegdist!
class(bac_jac) # needs to be class dist for pcoa()



#### PCoA -- w/ Jaccard matrix (from binary counts) ####

## fungi first

its1_pcoa2 = pcoa(its1_jac) # pcoa of Bray-Curtis dissimilarity matrix made w/ relativized abundance (site x species) data

biplot(its1_pcoa2, scale(meta_quant))

its1_pcoa2.vectors<-data.frame(its1_pcoa2$vectors)
its1_pcoa2.vectors$SampleID<-rownames(its1_pcoa2.vectors)

its1_pcoa2_meta<-merge(its1_pcoa2.vectors, metadata, by="SampleID")
head(its1_pcoa2_meta)

its1_pcoa2_df=data.frame(SampleID = its1_pcoa2_meta$SampleID, Axis1 = its1_pcoa2_meta$Axis.1, Axis2 = its1_pcoa2_meta$Axis.2, Site = its1_pcoa2_meta$Site, Elevation = its1_pcoa2_meta$Elevation, Month=its1_pcoa2_meta$Month, Year= its1_pcoa2_meta$Year) #create dataframe with PCoA axes and some metadata
head(its1_pcoa2_df)

its1_pcoa2_df$Month2 <- factor(its1_pcoa2_df$Month, levels = c("July","August","October")) 
its1_pcoa2_df$Elevation2 <- factor(its1_pcoa2_df$Elevation, levels = c("400","1100","2000","2700")) 

head(its1_pcoa2_df)

its1.pcoa2.fig.1<-ggplot(its1_pcoa2_df, aes(x=Axis1, y=Axis2, col=Month2, shape=factor(Year))) +geom_point(alpha=0.5,size=5)+theme_bw()+scale_colour_manual(values =  saturation(SM_pal, 1))+
  coord_cartesian(xlim=c(-0.5,0.5),ylim=c(-0.5,0.5))+labs(xlab="Axis 1", ylab="Axis 2",color="Month", shape="Year")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))

ggsave(its1.pcoa2.fig.1,filename = "figures/fungal_pcoa1_Jaccard_Sierra_3.29.21.pdf", width=8, height=6, dpi=600)

its1.pcoa2.fig.2<-ggplot(its1_pcoa2_df, aes(x=Axis1, y=Axis2, col=Month2, shape=factor(Year), size=Elevation2)) +geom_point(alpha=0.5)+theme_bw()+scale_colour_manual(values =  saturation(SM_pal, 1))+
  coord_cartesian(xlim=c(-0.5,0.5),ylim=c(-0.5,0.5))+labs(xlab="Axis 1", ylab="Axis 2",color="Month", shape="Year",size="Elevation")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(color = guide_legend(override.aes = list(size = 5)),shape = guide_legend(override.aes = list(size = 5)))

ggsave(its1.pcoa2.fig.2,filename = "figures/fungal_pcoa2_Jaccard_Sierra_3.29.21.pdf", width=8, height=6, dpi=600)

## bacteria + archaea next

bac_pcoa2 = pcoa(bac_jac) # pcoa of Bray-Curtis dissimilarity matrix made w/ relativized abundance (site x species) data

biplot(bac_pcoa2, scale(meta_quant))

bac_pcoa2.vectors<-data.frame(bac_pcoa2$vectors)
bac_pcoa2.vectors$SampleID<-rownames(bac_pcoa2.vectors)

bac_pcoa2_meta<-merge(bac_pcoa2.vectors, metadata, by="SampleID")
head(bac_pcoa2_meta)

bac_pcoa2_df=data.frame(SampleID = bac_pcoa2_meta$SampleID, Axis1 = bac_pcoa2_meta$Axis.1, Axis2 = bac_pcoa2_meta$Axis.2, Site = bac_pcoa2_meta$Site, Elevation = bac_pcoa2_meta$Elevation, Month=bac_pcoa2_meta$Month, Year= bac_pcoa2_meta$Year) #create dataframe with PCoA axes and some metadata
head(bac_pcoa2_df)

bac_pcoa2_df$Month2 <- factor(bac_pcoa2_df$Month, levels = c("July","August","October")) 
bac_pcoa2_df$Elevation2 <- factor(bac_pcoa2_df$Elevation, levels = c("400","1100","2000","2700")) 

head(bac_pcoa2_df)

bac.pcoa2.fig.1<-ggplot(bac_pcoa2_df, aes(x=Axis1, y=Axis2, col=Month2, shape=factor(Year))) +geom_point(alpha=0.5,size=5)+theme_bw()+scale_colour_manual(values =  saturation(SM_pal, 1))+
  coord_cartesian(xlim=c(-0.5,0.5),ylim=c(-0.5,0.5))+labs(xlab="Axis 1", ylab="Axis 2",color="Month", shape="Year")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))

ggsave(bac.pcoa2.fig.1,filename = "figures/microbial_pcoa1_Jaccard_Sierra_3.29.21.pdf", width=8, height=6, dpi=600)

bac.pcoa2.fig.2<-ggplot(bac_pcoa2_df, aes(x=Axis1, y=Axis2, col=Month2, shape=factor(Year), size=Elevation2)) +geom_point(alpha=0.5)+theme_bw()+scale_colour_manual(values =  saturation(SM_pal, 1))+
  coord_cartesian(xlim=c(-0.5,0.5),ylim=c(-0.5,0.5))+labs(xlab="Axis 1", ylab="Axis 2",color="Month", shape="Year",size="Elevation")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(color = guide_legend(override.aes = list(size = 5)),shape = guide_legend(override.aes = list(size = 5)))

ggsave(bac.pcoa2.fig.2,filename = "figures/microbial_pcoa2_Jaccard_Sierra_3.29.21.pdf", width=8, height=6, dpi=600)




##### Bray-Curtis dissimilarity matrix from RAW OTU counts ####

# ITS1 first
its1.bray<-vegdist(its1_otu_table, method="bray")
# pcoa of Bray-Curtis dissimilarity matrix made w/ relativized abundance (site x species) data
its1.sq.bray<-sqrt(its1.bray)
is.euclid(its1.sq.bray)

# 16S next
## compute a matrix of from raw OTU counts
bac.bray<-vegdist(bac_otu_table, method="bray")
# pcoa of Bray-Curtis dissimilarity matrix made w/ relativized abundance (site x species) data
bac.sq.bray<-sqrt(bac.bray)
is.euclid(bac.sq.bray)


#### PCoA w/ Bray-Curtis matrix from raw OTU counts####
## ** PCoA: build Bray-Curtis dissimilarity matrix from raw OTU counts, then take square root of matrix before PCoA...
### OR: build Bray-Curtis dissimilarity matrix from raw OTU counts, then use correction when running PCoA

## fungi first
its1.pcoa = pcoa(its1.sq.bray) 
its1.pcoa ## Check to see if negative eigenvalues affect the interpretation of the first several axes
# The proportion of variances explained is in its element values$Relative_eig

str(its1.pcoa)
biplot(its1.pcoa)
## project (scaled) chem data on the PCoA ordination
biplot(its1.pcoa, meta_qdat[,2:9], main="ITS1 PCoA + Scaled Chemical Data")
biplot(its1.pcoa, meta_qdat[,2:9],display = c("sites", "species"),type = c("text","points"))


## 16S next...

bac.pcoa = pcoa(bac.sq.bray) 
bac.pcoa ## Check to see if negative eigenvalues affect the interpretation of the first several axes
# The proportion of variances explained is in its element values$Relative_eig
#cmdscale(bac.bray, k=(nrow(bac_otu_table)-1), eig=TRUE)

biplot(bac.pcoa)
## project (scaled) chem data on the PCoA ordination
biplot(bac.pcoa, meta_qdat[,2:9], main="16S PCoA + Scaled Chemical Data")
biplot(its1.pcoa, meta_qdat[,2:9],display = c("sites", "species"),type = c("text","points"))



#### Fungal PCoA Visualization (raw counts) ####

### fungi first

#tiff(file="figures/its1_pcoa_biplot.tiff",width=10, height=10, units="in", res=600)
biplot(its1.pcoa, scale(meta_quant))
dev.off()

# setting up the dataframe for ggplot2 visualization

its1.pcoa.vectors<-data.frame(its1.pcoa$vectors)
its1.pcoa.vectors$SampleID<-rownames(its1.pcoa.vectors)

its1.pcoa_meta<-merge(its1.pcoa.vectors, meta_all_scaled, by.x="SampleID", by.y="SampleID")
its1.pcoa_df=data.frame(SampleID = its1.pcoa_meta$SampleID, Axis1 = its1.pcoa_meta$Axis.1, Axis2 = its1.pcoa_meta$Axis.2, Site = its1.pcoa_meta$Site, Elevation = its1.pcoa_meta$Elevation, Month=its1.pcoa_meta$Month, Year= its1.pcoa_meta$Year, DustComp=its1.pcoa_meta$DustComplex) #create dataframe with PCoA axes and some metadata

its1.pcoa_df$Month2 <- factor(its1.pcoa_df$Month, levels = c("July","August","October")) 
its1.pcoa_df$Elevation2 <- factor(its1.pcoa_df$Elevation, levels = c("400","1100","2000","2700")) 
its1.pcoa_df$Site2 <- factor(its1.pcoa_df$Site, levels = c("SJER","SOAPROOT","PROVIDENCE","SHORTHAIR")) 

head(its1.pcoa_df)

# plot time

its1.pcoa.fig.1<-ggplot(its1.pcoa_df, aes(x=Axis1, y=Axis2, col=Month2, shape=factor(Year))) +geom_point(alpha=0.5,size=5)+theme_bw()+scale_colour_manual(values =  saturation(SM_pal, 1))+
  coord_cartesian(xlim=c(-0.5,0.5),ylim=c(-0.5,0.5))+labs(xlab="Axis 1", ylab="Axis 2",color="Month", shape="Year")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))

ggsave(its1.pcoa.fig.1,filename = "figures/fungal_pcoa1_BC_Sierra_4.18.21.pdf", width=8, height=6, dpi=600)

its1.pcoa.fig.2<-ggplot(its1.pcoa_df, aes(x=Axis1, y=Axis2, shape=factor(Year), col=factor(Elevation))) +geom_point(alpha=0.5, size=5)+theme_bw()+scale_colour_manual(values = saturation(warm2cold1, 2))+
  coord_cartesian(xlim=c(-0.5,0.5),ylim=c(-0.5,0.5))+labs(xlab="Axis 1", ylab="Axis 2",color="Elevation", shape="Year",size="Elevation")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(color = guide_legend(override.aes = list(size = 5)),shape = guide_legend(override.aes = list(size = 5)))

ggsave(its1.pcoa.fig.2,filename = "figures/fungal_pcoa2_BC_Sierra_4.18.21.pdf", width=8, height=6, dpi=600)

its1.pcoa.fig.3<-ggplot(its1.pcoa_df, aes(x=Axis1, y=Axis2, col=Site2, shape=Month2)) +geom_point(alpha=0.5,size=5)+theme_bw()+scale_colour_manual(values =  saturation(warm2cold2, 1))+
  coord_cartesian(xlim=c(-0.5,0.5),ylim=c(-0.5,0.5))+labs(xlab="Axis 1", ylab="Axis 2",color="Site", shape="Month")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))

ggsave(its1.pcoa.fig.3,filename = "figures/fungal_pcoa3_BC_Sierra_4.18.21.pdf", width=8, height=6, dpi=600)

its1.pcoa.fig.4<-ggplot(its1.pcoa_df, aes(x=Axis1, y=Axis2, col=Site2, shape=Elevation2)) +geom_point(alpha=0.5,size=5)+theme_bw()+scale_colour_manual(values =  saturation(warm2cold2, 1))+
  coord_cartesian(xlim=c(-0.5,0.5),ylim=c(-0.5,0.5))+labs(xlab="Axis 1", ylab="Axis 2",color="Site", shape="Elevation")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))

ggsave(its1.pcoa.fig.4,filename = "figures/fungal_pcoa4_BC_Sierra_4.18.21.pdf", width=8, height=6, dpi=600)

its1.pcoa.fig.5a<-ggplot(its1.pcoa_df, aes(x=Axis1, y=Axis2, col=Elevation2, shape=Month2)) +geom_point(alpha=0.5,size=5)+theme_bw()+scale_colour_manual(values =  saturation(warm2cold1, 1))+
  coord_cartesian(xlim=c(-0.5,0.5),ylim=c(-0.5,0.5))+labs(xlab="Axis 1", ylab="Axis 2",color="Elevation", shape="Month")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))

ggsave(its1.pcoa.fig.5a,filename = "figures/fungal_pcoa5a_BC_Sierra_4.18.21.pdf", width=8, height=6, dpi=600)

# Dust Copmlexity PCoAs
f1.dust.pcoa1<-ggplot(its1.pcoa_df, aes(x=Axis1, y=Axis2, col=DustComp, shape=factor(Year))) +geom_point(alpha=0.5,size=5)+theme_bw()+scale_colour_gradientn(colours=fair_cols)+
  coord_cartesian(xlim=c(-0.5,0.5),ylim=c(-0.5,0.5))+labs(xlab="Axis 1", ylab="Axis 2",color="Dust Complexity", shape="Year")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))
ggsave(f1.dust.pcoa1,filename = "figures/Sierra_ITS1_Dust_PCOAraw_5.3.21.pdf", width=8, height=6, dpi=600)

f1.dust.pcoa2<-ggplot(its1.pcoa_df, aes(x=Axis1, y=Axis2, col=DustComp, shape=factor(Elevation))) +geom_point(alpha=0.5,size=5)+theme_bw()+scale_shape_manual(values = c(18, 16, 17,15))+scale_colour_gradientn(colours=fair_cols)+
  coord_cartesian(xlim=c(-0.5,0.5),ylim=c(-0.5,0.5))+labs(xlab="Axis 1", ylab="Axis 2",color="Dust Complexity", shape="Elevation (ft)")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))
ggsave(f1.dust.pcoa2,filename = "figures/Sierra_ITS1_Dust_PCOAraw2_5.3.21.pdf", width=8, height=6, dpi=600)

f1.dust.pcoa3<-ggplot(its1.pcoa_df, aes(x=Axis1, y=Axis2, col=DustComp, shape=factor(interaction(Month,Year)))) +geom_point(alpha=0.5,size=5)+theme_bw()+scale_shape_manual(name="Month & Year", labels=c("August 2014","July 2014","October 2014","July 2015","October 2015"),values = c(20, 18, 19, 17,15))+scale_colour_gradientn(colours=fair_cols)+
  coord_cartesian(xlim=c(-0.5,0.5),ylim=c(-0.5,0.5))+labs(xlab="Axis 1", ylab="Axis 2",color="Dust Complexity", shape="Elevation (ft)")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))
ggsave(f1.dust.pcoa3,filename = "figures/Sierra_ITS1_Dust_PCOAraw3_5.3.21.pdf", width=8, height=6, dpi=600)

# its1.pcoa.fig.5b<-ggplot(its1_pcoa1_df, aes(x=Axis1, y=Axis2, col=Month2, shape=Elevation2)) +geom_point(alpha=0.5,size=5)+theme_bw()+scale_colour_manual(values =  saturation(Sag_pal, 1))+
#   coord_cartesian(xlim=c(-0.5,0.5),ylim=c(-0.5,0.5))+labs(xlab="Axis 1", ylab="Axis 2",color="Month", shape="Elevation")+theme_classic()+
#   theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
#   guides(shape = guide_legend(override.aes = list(size = 5)))
# 
# ggsave(its1.pcoa.fig.5b,filename = "figures/fungal_pcoa5b_BC_Sierra_4.18.21.pdf", width=8, height=6, dpi=600)


#### Microbial PCoA Visualization (raw counts) ####


#tiff(file="figures/bac_pcoa_biplot.tiff",width=10, height=10, units="in", res=600)
biplot(bac.pcoa, scale(meta_quant))
dev.off()

# setting up the dataframe for ggplot2 visualization

bac.pcoa.vectors<-data.frame(bac.pcoa$vectors)
bac.pcoa.vectors$SampleID<-rownames(bac.pcoa.vectors)

bac.pcoa_meta<-merge(bac.pcoa.vectors, meta_all_scaled, by="SampleID")
bac.pcoa_df=data.frame(SampleID = bac.pcoa_meta$SampleID, Axis1 = bac.pcoa_meta$Axis.1, Axis2 = bac.pcoa_meta$Axis.2, Site = bac.pcoa_meta$Site,Elevation = bac.pcoa_meta$Elevation, Month=bac.pcoa_meta$Month, Year= bac.pcoa_meta$Year, DustComp=bac.pcoa_meta$DustComplex) #create dataframe with PCoA axes and some metadata

bac.pcoa_df$Month2 <- factor(bac.pcoa_df$Month, levels = c("July","August","October")) 
bac.pcoa_df$Elevation2 <- factor(bac.pcoa_df$Elevation, levels = c("400","1100","2000","2700")) 
bac.pcoa_df$Site2 <- factor(bac.pcoa_df$Site, levels = c("SJER","SOAPROOT","PROVIDENCE","SHORTHAIR")) 

head(bac.pcoa_df)

bac.pcoa.fig.1<-ggplot(bac.pcoa_df, aes(x=Axis1, y=Axis2, col=Month2, shape=factor(Year))) +geom_point(alpha=0.5,size=5)+theme_bw()+scale_colour_manual(values =  saturation(SM_pal, 1))+
  coord_cartesian(xlim=c(-0.5,0.5),ylim=c(-0.5,0.5))+labs(xlab="Axis 1", ylab="Axis 2",color="Month", shape="Year")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))

ggsave(bac.pcoa.fig.1,filename = "figures/microbial_pcoa1_BC_Sierra_4.18.21.pdf", width=8, height=6, dpi=600)

bac.pcoa.fig.2<-ggplot(bac.pcoa_df, aes(x=Axis1, y=Axis2, col=Elevation2, shape=factor(Year))) +geom_point(alpha=0.5,size=5)+theme_bw()+scale_colour_manual(values =  saturation(warm2cold1, 1))+
  coord_cartesian(xlim=c(-0.5,0.5),ylim=c(-0.5,0.5))+labs(xlab="Axis 1", ylab="Axis 2",color="Month", shape="Year",size="Elevation")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(color = guide_legend(override.aes = list(size = 5)),shape = guide_legend(override.aes = list(size = 5)))

ggsave(bac.pcoa.fig.2,filename = "figures/microbial_pcoa2_BC_Sierra_4.18.21.pdf", width=8, height=6, dpi=600)

bac.pcoa.fig.3<-ggplot(bac.pcoa_df, aes(x=Axis1, y=Axis2, col=Site2, shape=Month2)) +geom_point(alpha=0.5,size=5)+theme_bw()+scale_colour_manual(values =  saturation(warm2cold2, 1))+
  coord_cartesian(xlim=c(-0.5,0.5),ylim=c(-0.5,0.5))+labs(xlab="Axis 1", ylab="Axis 2",color="Site", shape="Month")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))

ggsave(bac.pcoa.fig.3,filename = "figures/microbial_pcoa3_BC_Sierra_4.18.21.pdf", width=8, height=6, dpi=600)

bac.pcoa.fig.4<-ggplot(bac.pcoa_df, aes(x=Axis1, y=Axis2, col=Site2, shape=Elevation2)) +geom_point(alpha=0.5,size=5)+theme_bw()+scale_colour_manual(values =  saturation(warm2cold2, 1))+
  coord_cartesian(xlim=c(-0.5,0.5),ylim=c(-0.5,0.5))+labs(xlab="Axis 1", ylab="Axis 2",color="Site", shape="Elevation")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))

ggsave(bac.pcoa.fig.4,filename = "figures/microbial_pcoa4_BC_Sierra_4.18.21.pdf", width=8, height=6, dpi=600)

bac.pcoa.fig.5a<-ggplot(bac.pcoa_df, aes(x=Axis1, y=Axis2, col=Elevation2, shape=Month2)) +geom_point(alpha=0.5,size=5)+theme_bw()+scale_colour_manual(values =  saturation(warm2cold1, 1))+
  coord_cartesian(xlim=c(-0.5,0.5),ylim=c(-0.5,0.5))+labs(xlab="Axis 1", ylab="Axis 2",color="Elevation", shape="Month")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))

ggsave(bac.pcoa.fig.5a,filename = "figures/microbial_pcoa5a_BC_Sierra_4.18.21.pdf", width=8, height=6, dpi=600)

# Dust Complexity PCoAs
b1.dust.pcoa1<-ggplot(bac.pcoa_df, aes(x=Axis1, y=Axis2, col=DustComp, shape=factor(Year))) +geom_point(alpha=0.5,size=5)+theme_bw()+scale_colour_gradientn(colours=fair_cols)+
  coord_cartesian(xlim=c(-0.5,0.5),ylim=c(-0.5,0.5))+labs(xlab="Axis 1", ylab="Axis 2",color="Dust Complexity", shape="Year")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))
ggsave(b1.dust.pcoa1,filename = "figures/Sierra_16S_Dust_PCOAraw_5.3.21.pdf", width=8, height=6, dpi=600)

b1.dust.pcoa2<-ggplot(bac.pcoa_df, aes(x=Axis1, y=Axis2, col=DustComp, shape=factor(Elevation))) +geom_point(alpha=0.5,size=5)+theme_bw()+scale_shape_manual(values = c(18, 16, 17,15))+scale_colour_gradientn(colours=fair_cols)+
  coord_cartesian(xlim=c(-0.5,0.5),ylim=c(-0.5,0.5))+labs(xlab="Axis 1", ylab="Axis 2",color="Dust Complexity", shape="Elevation (ft)")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))
ggsave(b1.dust.pcoa2,filename = "figures/Sierra_16S_Dust_PCOAraw2_5.3.21.pdf", width=8, height=6, dpi=600)

b1.dust.pcoa3<-ggplot(bac.pcoa_df, aes(x=Axis1, y=Axis2, col=DustComp, shape=factor(interaction(Month,Year)))) +geom_point(alpha=0.5,size=5)+theme_bw()+scale_shape_manual(name="Month & Year", labels=c("August 2014","July 2014","October 2014","July 2015","October 2015"),values = c(20, 18, 19, 17,15))+scale_colour_gradientn(colours=fair_cols)+
  coord_cartesian(xlim=c(-0.5,0.5),ylim=c(-0.5,0.5))+labs(xlab="Axis 1", ylab="Axis 2",color="Dust Complexity", shape="Elevation (ft)")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))
ggsave(b1.dust.pcoa3,filename = "figures/Sierra_16S_Dust_PCOAraw3_5.3.21.pdf", width=8, height=6, dpi=600)

# bac.pcoa.fig.5b<-ggplot(bac_pcoa1_df, aes(x=Axis1, y=Axis2, col=Month2, shape=Elevation2)) +geom_point(alpha=0.5,size=5)+theme_bw()+scale_colour_manual(values =  saturation(Sag_pal, 1))+
#   coord_cartesian(xlim=c(-0.5,0.5),ylim=c(-0.5,0.5))+labs(xlab="Axis 1", ylab="Axis 2",color="Month", shape="Elevation")+theme_classic()+
#   theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
#   guides(shape = guide_legend(override.aes = list(size = 5)))
# 
# ggsave(bac.pcoa.fig.5b,filename = "figures/microbial_pcoa5b_BC_Sierra_4.18.21.pdf", width=8, height=6, dpi=600)
