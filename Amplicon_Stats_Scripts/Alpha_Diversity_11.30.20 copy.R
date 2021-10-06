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
library(ALDEx2)
library(breakaway)
library(DivNet)
library(compositions)

setwd("/Volumes/HLF_SSD/Aronson_Lab_Data/Denise_Data")
getwd()


#### Import and reformat the data ####
otus_counts_taxa<- as.data.frame(read.csv("FilteredTable_clean.csv"))

dim(otus_counts_taxa)
tail(otus_counts_taxa) ## make sure that OTUs are rows and samples are columns!!

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

otus_counts_taxa<-empty_to_unknown(otus_counts_taxa)
tail(otus_counts_taxa) ## make sure that the function worked
head(otus_counts_taxa)

otu_counts_all<-as.data.frame(otus_counts_taxa) # making copy of OTU table that includes Unknown taxa (at Kingdom level)
head(otu_counts_all)

otus_counts_taxa<-subset(otus_counts_taxa, Kingdom!="Unknown") ## keep only bacteria and archaean -- drop Unknown at Kingdom level
'Unknown' %in% otus_counts_taxa # check if Chloroplast counts are still in df, should be false because they've been removed

otu_counts<-subset(otus_counts_taxa, select=-c(Kingdom, Phylum, Class, Order, Family, Genus, Species)) # subset only part of the metadata we need

head(otu_counts) # Need Site x OTU table; currently OTUs as rows, Samples as columns - must transpose thee dataframe!

rownames(otu_counts)<-otu_counts$OTUID
head(otu_counts)
otu_counts$OTUID<-NULL

otu_table<-data.frame(t(otu_counts)) # transposes df so that Sample IDs are row names, OTU IDs are column names
#head(otu_table)

otu_table<-otu_table[, which(colSums(otu_table) > 0)] # drop 0s

### Import metadata ####
metadata<-as.data.frame(read_xlsx("SPMapDJ.xlsx"))

head(metadata)
tail(metadata)
metadata<-subset(metadata, select=-c(Barcodesequences, LinkerPrimersequences)) # subset only part of the metadata we need

rownames(metadata)<-metadata$"#SampleID"
colnames(metadata)[which(names(metadata) == "#SampleID")] <- "SampleID"
tail(metadata)

metadata=metadata[rownames(otu_table),] ## reorder metadata to have same rows as original OTU table
# ** ^ this indexing method will only work if the two dfs have the same # of rows AND the same row names!

# are our rows the same in both dataframes?
checkrow<-(is.element(row.names(otu_table), row.names(metadata)))

ori<-as.factor(metadata$Origin)

ster<-as.factor(metadata$Sterility)

metadata$Ori_Ster<-interaction(metadata$Origin,metadata$Sterility)

#### CLR transformation ####

#clr.otu.table<-aldex.clr(otu_table, mc.samples = 128, verbose = FALSE, useMC=FALSE)

otu.table_clr<-microbiome::transform(otu_table, 'clr')


otu.table_log2<-as.data.frame(lapply(otu_table, function(x) log2(x+1)))
otu.table_log2<-otu.table_log2[order(match(otu.table_log2[,1],otu_table[,1])),]
# ordering meta_quant based on position order of matches from meta_quant[,1] found in metadata$Cu
# cannot call column using $ for array, only df
rownames(otu.table_log2)<-rownames(otu_table) # can only do this if you're sure the arrays/dfs are in the same exact order
rownames(otu.table_log2)

#### Alpha Diversity and Richness (raw data) ####

## We can also describe the data with diversity measures, such as:
SR <- rowSums(otu_table > 0) ## Species richness
specnumber(otu_table) ## you can also use this for Richness
H <- vegan::diversity(otu_table) ## Shannon entropy
Shan_div <- exp(H) ## Shannon's diversity (number of abundant species)
Simp_div <- vegan::diversity(otu_table, "inv") ## Simpson diversity (number of dominant species)
Eve <- H/log(SR) ## Pielou evenness
Shan_eve <- Shan_div/SR ## Shannon evenness (Hill's ratio)
Simp_eve <- Simp_div/SR ## Simpson evenness (Hill's ratio)
div <- data.frame(SR, H, Shan_div, Simp_div, Shan_eve, Simp_eve, Eve) ## create a dataframe for the above measures
head(div)
div$SampleID<-rownames(div)
dim(div)
div=div[rownames(metadata),] ## reorder metadata to have same rows as original OTU table

cor(div[,1:7])
pairs(div[,-8], pch = 19,  cex = 0.5,lower.panel=panel.smooth, main="Pearson Correlation Matrix")
library(psych)
pairs.panels(div[,-8],
             method = "pearson", # correlation method
             hist.col = "#00AFBB",
             density = TRUE,  # show density plots
             ellipses = TRUE # show correlation ellipses
)

otu_meta<-merge(metadata, div, by="SampleID")
rownames(otu_meta)<-otu_meta$SampleID
otu_meta=otu_meta[rownames(metadata),] ## reorder metadata to have same rows as original OTU table
head(otu_meta)

otu_meta$Ori_Ster<-interaction(otu_meta$Origin,otu_meta$Sterility)
head(otu_meta)
otu_meta$Ori_Ster <- factor(otu_meta$Ori_Ster,levels = c("wild.N", "wild.S", "lab.N", "lab.S"))
#metadata$Ori_Ster <- factor(metadata$Ori_Ster,levels = c("wild.N", "wild.S", "lab.N", "lab.S"))

labN<-subset(otu_meta, Ori_Ster=="lab.N")
wildN<-subset(otu_meta, Ori_Ster=="wild.N")
labS<-subset(otu_meta, Ori_Ster=="lab.S")
wildS<-subset(otu_meta, Ori_Ster=="wild.S")


#### Alpha Diversity and Richness (log2 transformed data) ####

## We can also describe the data with diversity measures, such as:
#SR <- rowSums(otu.table_log2 > 0) ## Species richness
#specnumber(otu.table_log2) ## you can also use this for Richness
H_l <- vegan::diversity(otu.table_log2) ## Shannon entropy
Shan_div_l <- exp(H_l) ## Shannon's diversity (number of abundant species)
Simp_div <- vegan::diversity(otu_table, "inv") ## Simpson diversity (number of dominant species)
Eve <- H/log(SR) ## Pielou evenness
Shan_eve <- Shan_div/SR ## Shannon evenness (Hill's ratio)
Simp_eve <- Simp_div/SR ## Simpson evenness (Hill's ratio)
div <- data.frame(SR, H, Shan_div, Simp_div, Shan_eve, Simp_eve, Eve) ## create a dataframe for the above measures
head(div)
div$SampleID<-rownames(div)
dim(div)
div=div[rownames(metadata),] ## reorder metadata to have same rows as original OTU table

otu_meta<-merge(metadata, div, by="SampleID")
rownames(otu_meta)<-otu_meta$SampleID
otu_meta=otu_meta[rownames(metadata),] ## reorder metadata to have same rows as original OTU table
head(otu_meta)

otu_meta$Ori_Ster<-interaction(otu_meta$Origin,otu_meta$Sterility)
head(otu_meta)
otu_meta$Ori_Ster <- factor(otu_meta$Ori_Ster,levels = c("wild.N", "wild.S", "lab.N", "lab.S"))
#metadata$Ori_Ster <- factor(metadata$Ori_Ster,levels = c("wild.N", "wild.S", "lab.N", "lab.S"))

labN<-subset(otu_meta, Ori_Ster=="lab.N")
wildN<-subset(otu_meta, Ori_Ster=="wild.N")
labS<-subset(otu_meta, Ori_Ster=="lab.S")
wildS<-subset(otu_meta, Ori_Ster=="wild.S")



#### Alpha Diversity and Richness (clr transformed data) ####

## We can also describe the data with diversity measures, such as:
#SR <- rowSums(otu.table_log2 > 0) ## Species richness
#specnumber(otu.table_log2) ## you can also use this for Richness
H_clr <- vegan::diversity(otu.table_clr) ## Shannon entropy
Shan_div_clr <- exp(H_clr) ## Shannon's diversity (number of abundant species)
Simp_div <- vegan::diversity(otu_table, "inv") ## Simpson diversity (number of dominant species)
Eve <- H/log(SR) ## Pielou evenness
Shan_eve <- Shan_div/SR ## Shannon evenness (Hill's ratio)
Simp_eve <- Simp_div/SR ## Simpson evenness (Hill's ratio)
div <- data.frame(SR, H, Shan_div, Simp_div, Shan_eve, Simp_eve, Eve) ## create a dataframe for the above measures
head(div)
div$SampleID<-rownames(div)
dim(div)
div=div[rownames(metadata),] ## reorder metadata to have same rows as original OTU table

otu_meta<-merge(metadata, div, by="SampleID")
rownames(otu_meta)<-otu_meta$SampleID
otu_meta=otu_meta[rownames(metadata),] ## reorder metadata to have same rows as original OTU table
head(otu_meta)

otu_meta$Ori_Ster<-interaction(otu_meta$Origin,otu_meta$Sterility)
head(otu_meta)
otu_meta$Ori_Ster <- factor(otu_meta$Ori_Ster,levels = c("wild.N", "wild.S", "lab.N", "lab.S"))
#metadata$Ori_Ster <- factor(metadata$Ori_Ster,levels = c("wild.N", "wild.S", "lab.N", "lab.S"))

labN<-subset(otu_meta, Ori_Ster=="lab.N")
wildN<-subset(otu_meta, Ori_Ster=="wild.N")
labS<-subset(otu_meta, Ori_Ster=="lab.S")
wildS<-subset(otu_meta, Ori_Ster=="wild.S")




#### Histogram to determine distribution of data ####
#png(file = "myplot.png")
hist(otu_meta$Shan_div, freq=FALSE)
hist(otu_meta$SR, freq=FALSE)

#### T - test or Wilcoxon to Compare Diversity by Group ####
## Using Shapiro-Wilk test for normality
# The null-hypothesis of this test is that the population is normally distributed.
# Thus, if the p value is < the chosen alpha level, then the null hypothesis is rejected and there is evidence that the data tested are not normally distributed.
# If the p value is > than the chosen alpha level, then the null hypothesis (that the data came from a normally distributed population) can not be rejected
# (e.g., for an alpha level of .05, a data set with a p value of less than .05 rejects the null hypothesis that the data are from a normally distributed population
# p-value > 0.05 implying that the distribution of the data are not significantly different from normal distribution. In other words, we can assume the normality

## Visualizing with Normal Q-Q Plots: orders data from smallest to largest and plots it vs. expected normal order statistics

shapiro.test(otu_meta$SR) # p-value = 0.09316
shapiro.test(otu_meta$Shan_div) # p-value = 0.0006323

qqnorm(otu_meta$SR, pch = 1, frame = FALSE)
qqline(otu_meta$SR, col = "steelblue", lwd = 2)

qqnorm(otu_meta$Shan_div, pch = 1, frame = FALSE)
qqline(otu_meta$Shan_div, col = "steelblue", lwd = 2)

wild_div<-subset(otu_meta, Origin=="wild") ## keep only bacteria and archaea!!!!!!!!!
lab_div<-subset(otu_meta, Origin=="lab") ## keep only bacteria and archaea!!!!!!!!!


## Shannon Diversity by Treatment Origin
t.test(otu_meta$Shan_div~otu_meta$Origin) # where y is numeric and x is a binary factor
wilcox.test(otu_meta$Shan_div~otu_meta$Origin)

## Shannon Diversity by Treatment Source
t.test(otu_meta$Shan_div~otu_meta$Sterility) # where y is numeric and x is a binary factor
wilcox.test(otu_meta$Shan_div~otu_meta$Sterility)

## Species Richness by Treatment Origin
t.test(otu_meta$SR~otu_meta$Origin) # where y is numeric and x is a binary factor
wilcox.test(otu_meta$SR~otu_meta$Origin)

## Species Richness by Treatment Source
t.test(otu_meta$SR~otu_meta$Sterility) # where y is numeric and x is a binary factor
wilcox.test(otu_meta$SR~otu_meta$Sterility)

t.test(otu_meta$Shan_div~otu_meta$Ori_Ster) # where y is numeric and x is a binary factor

t1<-ggtexttable(compare_means(Shan_div~Ori_Ster,data=otu_meta)) ### t test or Wilcoxon tests with compare_means!
t2<-ggtexttable(compare_means(SR~Ori_Ster,data=otu_meta)) ### t test or Wilcoxon tests with compare_means!
wil.table<-ggarrange(t1, t2, legend = "none", common.legend = FALSE, labels = c("Wilcoxon Test (Shannon Diversity)", "Wilcoxon Test (Species Richness)"),align="hv", ncol = 1, nrow = 2)
ggsave(wil.table,filename = "16S_Wilcoxon_Results_Comparisons_5.25.21.pdf", width=10, height=10, dpi=600)

#### ANOVAs based on Diversity, SR for Origin, Sterility ####

# Compute the analysis of variance
anova1<-aov(Shan_div~Origin, data = otu_meta)
# Summary of the analysis
summary(anova1)
TukeyHSD(anova1)

anova2<-aov(SR~Origin, data = otu_meta)
# Summary of the analysis
summary(anova2)
TukeyHSD(anova2)

anova3<-aov(Shan_div~Sterility, data = otu_meta)
# Summary of the analysis
summary(anova3)
TukeyHSD(anova3)

anova4<-aov(SR~Sterility, data = otu_meta)
# Summary of the analysis
summary(anova4)
TukeyHSD(anova4)


# res.aov2 <- aov(Shan_div ~ Origin + Sterility, data = otu_meta)
# summary(res.aov2)
#
# res.aov3 <- aov(SR ~ Origin + Sterility, data = otu_meta)
# summary(res.aov3)

res.aov4 <- aov(Shan_div ~ Origin*Sterility, data = otu_meta)
summary(res.aov4)

res.aov5 <- aov(SR ~ Origin*Sterility, data = otu_meta)
summary(res.aov5)


#### PERMANOVAs with relative abundance ####
min(rowSums(otu_table))
### * below we use the dcast() function to "cast" the data into a wide format based on given elements (column names), taking sum of "Count"
### * decostand(df, method="total) is the function (with argument total) used to get relative abundance of OTU table

otu_RA<-data.frame(decostand(otu_table, method="total", MARGIN=1, na.rm=TRUE))  # relative abundance of taxa data where everything is divided by col total (b/c Margin=1 meaning rows == SAMPLES in this case)
rowSums(otu_RA) # sanity check to make sure the transformation worked!
#otu_RA$SampleID<-rownames(otu_RA)
head(otu_RA)

adonis2(otu_RA ~ Origin*Sterility, data = metadata, permutations = 999, method = "bray")
adonis2(otu_RA ~ Ori_Ster, data = metadata, permutations = 999, method = "bray")
adonis2(Shan_div ~ Ori_Ster, data = otu_meta, permutations = 999, method = "bray")
adonis2(SR ~ Ori_Ster, data = otu_meta, permutations = 999, method = "bray")


install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(pairwiseAdonis)

div$SampleID<-NULL
pair.mod<-pairwise.adonis(otu_RA,metadata$Ori_Ster, p.adjust.m='bonferroni') # shows us variation for each sample to see which ones are different
pair.mod # in this case, all sites are different (F, R2)


#### Checking out some color pallettes before we visualize ####
wes1<-wes_palette("Chevalier1")
wes2<-wes_palette("Moonrise3")
wes3<-wes_palette("IsleofDogs1")
wes4<-wes_palette("GrandBudapest1")
wes5<-wes_palette("GrandBudapest2")

SM_pal <- park_palette("SmokyMountains") # create a palette and specify # of colors youw ant
Arc_pal <- park_palette("Arches") # create a palette and specify # of colors youw ant
CL_pal <- park_palette("CraterLake") # create a palette and specify # of colors youw ant
Sag_pal <- park_palette("Saguaro") # create a palette and specify # of colors youw ant
Aca_pal <- park_palette("Acadia") # create a palette and specify # of colors youw ant
DV_pal <- park_palette("DeathValley") # create a palette and specify # of colors youw ant
DV_pal2 <- park_palette("DeathValley", 1) # create a palette and specify # of colors youw ant
CI_pal <- park_palette("ChannelIslands") # create a palette and specify # of colors youw ant
Bad_pal <- park_palette("Badlands") # create a palette and specify # of colors youw ant
MR_pal2 <- park_palette("MtRainier", 1) # create a palette and specify # of colors youw ant
MR_pal <- park_palette("MtRainier") # create a palette and specify # of colors youw ant

HI_pal <- park_palette("Hawaii") # create a palette and specify # of colors youw ant

bright_pal<-get_palette(palette = paste0("#", c("FFBE0B","D03325","8338EC","3A86FF")), k = 4)

N_S_pal<-get_palette(palette = paste0("#", c("4C7DAE","F24C4C")), k = 2)

L_G_pal<-get_palette(palette = paste0("#", c("579E71","462255")), k = 2)

#fungal_pal<-get_palette(palette = paste0("#", c("222529","c72421","da7049","eedc2b","c5debf","b6f2fa","513a56","fafccf")), k = 8)
#fungal_pal2<-get_palette(palette = paste0("#", c("07060f","dd1f19","eda964","bba20d","44be45","74dabc","5e2479","5f3523")), k = 8)

#scale_fill_manual(values = wes_palette("IsleofDogs1"))

#### Alpha Diversity Visualization -- based on Raw Counts ####

# f1<-ggplot(otu_meta, aes(x=SampleID, y=SR)) + geom_bar(stat="identity",colour="black", fill=DV_pal2)+theme_bw()+theme_classic()+labs(title="Microbial Species Richness Across Samples", x="Sample ID", y="Species Richness")+theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5))
# ggsave(f1,filename = "16S_SpeciesRichness_Denise_11.30.20.pdf", width=20, height=10, dpi=600)
#
# f2<-ggplot(otu_meta, aes(x=SampleID, y=Shan_div)) + geom_bar(stat="identity",colour="black", fill=DV_pal2)+theme_bw()+theme_classic()+labs(title="Microbial Alpha Diversity Across Samples", x="Sample ID", y="Shannon Diversity")+theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5))
# ggsave(f2,filename = "16S_ShannonDiversity_Denise_11.30.20.pdf", width=20, height=10, dpi=600)
w.t1<-compare_means(Shan_div~Ori_Ster,data=otu_meta)
names(w.t1)[which(names(w.t1) == ".y.")] <- "Result"
write.table(w.t1,"Wilcoxon_Diversity_Results.csv",  sep = ",", row.names=FALSE)

w.t2<-compare_means(SR~Ori_Ster,data=otu_meta)
names(w.t2)[which(names(w.t2) == ".y.")] <- "Result"
write.table(w.t2,"Wilcoxon_Richness_Results.csv",  sep = ",", row.names=FALSE)

t1<-ggtexttable(compare_means(Shan_div~Ori_Ster,data=otu_meta))
t1.a<-ggtexttable(compare_means(Shan_div~Ori_Ster,data=otu_meta,method="t.test"))
t1.b<-ggtexttable(compare_means(Shan_div~Ori_Ster,data=otu_meta,method="anova"))
table1<-ggarrange(t1, t1.a, t1.b, legend = "none", common.legend = FALSE, labels = c("Wilcoxon Test", "T-test", "ANOVA"),align="hv", ncol = 1, nrow = 3)
ggsave(table1,filename = "16S_ShannonDiversity_Comparisons_OriginxSterility_6.21.21.pdf", width=15, height=15, dpi=600)

# f3<-ggplot(otu_meta, aes(x=Origin, y=Shan_div, fill=Origin)) + geom_boxplot(color="black") +labs(title = "Shannon Diversity by Treatment Source", x="Treatment Source", y="Shannon Diversity") + theme_classic() + scale_x_discrete(name ="Treatment Source", labels=c("lab"="Lab", "wild"="Garden")) +
#   theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15)) +
#   scale_fill_manual(values = brightness(L_G_pal, 0.6), name ="Treatment Source", labels=c("lab"="Lab", "wild"="Garden"))
# ggsave(f3,filename = "16S_ShannonDiversity_byTreatmentSource_Denise_5.8.21.pdf", width=15, height=10, dpi=600)
#
# f4<-ggplot(otu_meta, aes(x=Sterility, y=Shan_div, fill=Sterility)) + geom_boxplot(color="black") +labs(title = "Shannon Diversity by Treatment Condition", x="Treatment Condition", y="Shannon Diversity") + theme_classic() + scale_x_discrete(name ="Treatment Source", labels=c("N"="Natural", "S"="Sterile")) +
#   theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15)) +
#   scale_fill_manual(values = saturation(N_S_pal, 0.7), name ="Treatment Condition", labels=c("N"="Natural", "S"="Sterile"))
# ggsave(f4,filename = "16S_ShannonDiversity_byTreatmentCondition_Denise_5.8.21.pdf", width=15, height=10, dpi=600)

t2<-ggtexttable(compare_means(SR~Ori_Ster,data=otu_meta))
t2.a<-ggtexttable(compare_means(SR~Ori_Ster,data=otu_meta,method="t.test"))
t2.b<-ggtexttable(compare_means(SR~Ori_Ster,data=otu_meta,method="anova"))
table2<-ggarrange(t2, t2.a, t2.b, legend = "none", common.legend = FALSE, labels = c("Wilcoxon Test", "T-test", "ANOVA"),align="hv", ncol = 1, nrow = 3)
ggsave(table2,filename = "16S_SpeciesRichness_Comparisons_OriginxSterility_6.21.21.pdf", width=15, height=15, dpi=600)

# f5<-ggplot(otu_meta, aes(x=Origin, y=SR, fill=Origin)) + geom_boxplot(color="black") +labs(title = "Species Richness by Treatment Source", x="Treatment Source", y="Species Richness") + theme_classic() + scale_x_discrete(name ="Treatment Source", labels=c("lab"="Lab", "wild"="Garden")) +
#   theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15)) +
#   scale_fill_manual(values = brightness(L_G_pal, 0.6), name ="Treatment Source", labels=c("lab"="Lab", "wild"="Garden"))
# ggsave(f5,filename = "16S_SpeciesRichness_byTreatmentSource_Denise_5.8.21.pdf", width=15, height=10, dpi=600)
#
# f6<-ggplot(otu_meta, aes(x=Sterility, y=SR, fill=Sterility)) + geom_boxplot(color="black") +labs(title = "Species Richness by Treatment Condition", x="Treatment Condition", y="Species Richness") + theme_classic() + scale_x_discrete(name ="Treatment Source", labels=c("N"="Natural", "S"="Sterile")) +
#   theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15)) +
#   scale_fill_manual(values = saturation(N_S_pal, 0.7), name ="Treatment Condition", labels=c("N"="Natural", "S"="Sterile"))
# ggsave(f6,filename = "16S_SpeciesRichness_byTreatmentCondition_Denise_5.8.21.pdf", width=15, height=10, dpi=600)


## comment out multiple lines with ctrl+shift+C
#
# a<-ggplot(wild_div, aes(x=SampleID, y=SR, fill=Sterility)) + geom_bar(stat="identity",colour="black")+theme_bw()+theme_classic()+scale_x_discrete(name ="Treatment Group", labels=c("lab"="Lab", "wild"="Garden"))+
#   labs(title="Bacterial Species Richness in Garden Slugs", x="Sample ID", y="Species Richness")+theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5))+
#   scale_fill_manual(values = saturation(CL_pal, 1), name ="Treatment Source", labels=c("N"="Natural", "S"="Sterile"))
#
# b<-ggplot(lab_div, aes(x=SampleID, y=SR, fill=Sterility)) + geom_bar(stat="identity",colour="black")+theme_bw()+theme_classic()+scale_x_discrete(name ="Treatment Group", labels=c("lab"="Lab", "wild"="Garden"))+
#   labs(title="Bacterial Species Richness in Lab-Reared Slugs", x="Sample ID", y="Species Richness")+theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5))+
#   scale_fill_manual(values = saturation(CL_pal, 1), name ="Treatment Source", labels=c("N"="Natural", "S"="Sterile"))
#
# SR_otu<-ggarrange(a, b, legend = "right", common.legend = TRUE, labels = c("A", "B"),ncol = 2, nrow = 1)
#
# ggsave(bac_phy_ori.st,filename = "16S_phyla_byTreatment_Denise_10.13.20.pdf", width=20, height=10, dpi=600)

ggplot(otu_meta, aes(x=SR, y=Shan_div)) + geom_point(size=2)
otu_meta$Ori_Ster <- factor(otu_meta$Ori_Ster ,levels = c("lab.N", "wild.N", "lab.S", "wild.S"))

# Species Richness

f7<-ggplot(otu_meta, aes(x=Ori_Ster, y=SR, fill=Ori_Ster)) + geom_boxplot(color="black") +labs(x="Treatment Condition by Source", y="Species Richness") + theme_classic() + scale_x_discrete(name ="Treatment Condition by Source", labels=c("lab.N"="Lab x Natural", "wild.N"="Garden x Natural", "lab.S"="Lab x Sterile", "wild.S"="Garden x Sterile")) +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15))+
  stat_compare_means(comparisons = list(c(1,2), c(3,4), c(2,3), c(2,4), c(1,3), c(1,4)), hide.ns = FALSE,label = "p.signif")+scale_fill_grey(name ="Treatment Condition by Source", labels=c("lab.N"="Lab x Natural", "wild.N"="Garden x Natural", "lab.S"="Lab x Sterile", "wild.S"="Garden x Sterile"))
ggsave(f7,filename = "16S_SpeciesRichness_ALL_ConditionXSource_Denise_5.8.21.pdf", width=15, height=10, dpi=600)

f8<-ggplot(otu_meta, aes(x=Ori_Ster, y=SR, fill=Ori_Ster)) + geom_boxplot(color="black") +labs(x="Treatment Condition by Source", y="Species Richness") + theme_classic() + scale_x_discrete(name ="Treatment Condition by Source", labels=c("lab.N"="Lab x Natural", "wild.N"="Garden x Natural", "lab.S"="Lab x Sterile", "wild.S"="Garden x Sterile")) +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15)) +
  scale_fill_grey(name ="Treatment Condition by Source", labels=c("lab.N"="Lab x Natural", "wild.N"="Garden x Natural", "lab.S"="Lab x Sterile", "wild.S"="Garden x Sterile"))+stat_compare_means(comparisons = list(c(1,3), c(2,3)), hide.ns = FALSE,label = "p.signif")
ggsave(f8,filename = "16S_SpeciesRichness_SIGNIFICANT_ConditionXSource_Denise_5.8.21.pdf", width=15, height=10, dpi=600)

sr_nobars<-ggplot(otu_meta, aes(x=Ori_Ster, y=SR, fill=Ori_Ster)) + geom_boxplot(color="black") +labs(x="Treatment Condition by Source", y="Species Richness") + theme_classic() + scale_x_discrete(name ="Treatment Condition by Source", labels=c("lab.N"="Lab x Natural", "wild.N"="Garden x Natural", "lab.S"="Lab x Sterile", "wild.S"="Garden x Sterile")) +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15))+
  scale_fill_grey(name ="Treatment Condition by Source", labels=c("lab.N"="Lab x Natural", "wild.N"="Garden x Natural", "lab.S"="Lab x Sterile", "wild.S"="Garden x Sterile"))
ggsave(sr_nobars,filename = "16S_NoSigBars_SpeciesRichness_ConditionXSource_Denise_5.25.21.pdf", width=15, height=10, dpi=600)

## Shannon diversity
f9<-ggplot(otu_meta, aes(x=Ori_Ster, y=Shan_div, fill=Ori_Ster)) + geom_boxplot(color="black") +labs(x="Treatment Condition by Source", y="Shannon Diversity") + theme_classic() + scale_x_discrete(name ="Treatment Condition by Source", labels=c("lab.N"="Lab x Natural", "wild.N"="Garden x Natural", "lab.S"="Lab x Sterile", "wild.S"="Garden x Sterile")) +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15)) +
  scale_fill_grey(name ="Treatment Condition by Source", labels=c("lab.N"="Lab x Natural", "lab.S"="Lab x Sterile", "wild.N"="Garden x Natural", "wild.S"="Garden x Sterile"))+stat_compare_means(comparisons = list(c(1,2), c(3,4), c(2,3), c(2,4), c(1,3), c(1,4)), hide.ns = FALSE,label = "p.signif")
ggsave(f9,filename = "16S_ShannonDiversity_ALL_ConditionXSource_Denise_5.8.21.pdf", width=15, height=10, dpi=600)

f10<-ggplot(otu_meta, aes(x=Ori_Ster, y=Shan_div, fill=Ori_Ster)) + geom_boxplot(color="black") +labs(x="Treatment Condition by Source", y="Shannon Diversity") + theme_classic() + scale_x_discrete(name ="Treatment Condition by Source", labels=c("lab.N"="Lab x Natural", "wild.N"="Garden x Natural", "lab.S"="Lab x Sterile", "wild.S"="Garden x Sterile")) +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15)) +
  scale_fill_grey(name ="Treatment Condition by Source", labels=c("lab.N"="Lab x Natural", "lab.S"="Lab x Sterile", "wild.N"="Garden x Natural", "wild.S"="Garden x Sterile"))+stat_compare_means(comparisons = list(c(2,3), c(3,4)), hide.ns = FALSE,label = "p.signif")
ggsave(f10,filename = "16S_ShannonDiversity_SIGNIFICANT_ConditionXSource_Denise_5.8.21.pdf", width=15, height=10, dpi=600)

shan_nobars<-ggplot(otu_meta, aes(x=Ori_Ster, y=Shan_div, fill=Ori_Ster)) + geom_boxplot(color="black") +labs(x="Treatment Condition by Source", y="Shannon Diversity") + theme_classic() + scale_x_discrete(name ="Treatment Condition by Source", labels=c("lab.N"="Lab x Natural", "wild.N"="Garden x Natural", "lab.S"="Lab x Sterile", "wild.S"="Garden x Sterile")) +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15))+
  scale_fill_grey(name ="Treatment Condition by Source", labels=c("lab.N"="Lab x Natural", "wild.N"="Garden x Natural", "lab.S"="Lab x Sterile", "wild.S"="Garden x Sterile"))
ggsave(shan_nobars,filename = "16S_NoSigBars_ShanDiv_ConditionXSource_Denise_5.25.21.pdf", width=15, height=10, dpi=600)

## Combo figures

combo_f1<-ggarrange(f8, f10, legend = "none", common.legend = FALSE, labels = c("A", "B"), align="hv", ncol = 1, nrow = 2)
ggsave(combo_f1,filename = "16S_Significant_Comparisons_ConditionXSource_6.21.21.pdf", width=20, height=10, dpi=800)
ggsave(combo_f1,filename = "16S_Significant_Comparisons_ConditionXSource_6.21.21.tiff", width=20, height=10, dpi=400)

combo_f2<-ggarrange(sr_nobars, shan_nobars, legend = "none", common.legend = FALSE, labels = c("A", "B"), align="hv", ncol = 1, nrow = 2)
ggsave(combo_f2,filename = "16S_NoBars_Comparisons_ConditionXSource_6.21.21.pdf", width=10, height=25, dpi=800)
ggsave(combo_f2,filename = "16S_NoBars_Comparisons_ConditionXSource_6.21.21.tiff", width=20, height=10, dpi=400)

#### Standardize Environmental Data for later ####

## standardize the data ...
## accomplished with a z-transformation (mean = 0, SD = 1) using the function 'scale'.

head(chem_meta)
scale(chem_meta) ## look at new scaled values as a safety check
chem_dist<-vegdist(scale(chem_meta), "euclidean") ## may need to transform these data prior to this step b/c comparing elements at varying concentraitons

coldiss(chem_dist, nc=31, diag=TRUE) # nc = # of colors

dev.off()
## quickly lets look a the heat maps for bray-Curtis and ED of the environment at the same time.
coldiss(spe.dbrel, byrank=FALSE, diag=TRUE)
coldiss(chem_dist, nc=16, diag=TRUE)

dev.off()














#### Rarefaction ####
# RAREFACTION in vegan: ROWS need to be SITES/samples; COLUMNS are SPECIES (OTUs, ASVs)
# see "Numerical Ecology with R", page 13-14
head(otu_table)

min(rowSums(otu_table)) ## seeing min sum of OTUs in a site — use this min for rarefaction
max(rowSums(otu_table)) # checking to see the max # of OTUs in a site

min<-min(rowSums(otu_table)) ## get min # of OTUs in a site for rarefaction

# Rarefy! vvvv
otu_table.r<-rrarefy(otu_table,min) ## be cognizant of min for rarefaction

max(rowSums(otu_table.r)) # SANITY CHECK: this should now be equivalent to the min from pre-rarefied data!
head(otu_table.r)


#### Alpha Diversity and Richness (rarefied data) ####

## We can also describe the data with diversity measures, such as:
H.r <- diversity(otu_table.r) ## Shannon entropy
Shan_div.r <- exp(H.r) ## Shannon's diversity (number of abundant species)
Simp_div.r <- diversity(otu_table.r, "inv") ## Simpson diversity (number of dominant species)
Eve.r <- H.r/log(SR) ## Pielou evenness
Shan_eve.r <- Shan_div/SR ## Shannon evenness (Hill's ratio)
Simp_eve.r <- Simp_div/SR ## Simpson evenness (Hill's ratio)
div.rar <- data.frame(H.r, Shan_div.r, Simp_div.r, Shan_eve.r, Simp_eve.r, Eve.r, SR) ## create a dataframe for the above measures
head(div.rar)
div.rar$SampleID<-rownames(div.rar)
dim(div.rar)
head(div)
otu.rar_meta<-merge(metadata, div.rar, by="SampleID")
head(otu.rar_meta)
## so now let's see if diversity differs between soils types
## but how do we know if the data are going to match since they are in different data sets
## one workaround is to have everything in a single dataset and things will always be in the same order

## let's make a few plots

par(mfrow=c(2,2))
boxplot(SR~metadata$Origin, xlab="Slug Origin", ylab="Species Richness")
boxplot(Shan_div~metadata$Origin, xlab="Slug Origin", ylab="Shannon's diversity")
boxplot(Simp_div~metadata$Origin, xlab="Slug Origin", ylab="Simpson diversity")
boxplot(Eve~metadata$Origin, xlab="Slug Origin", ylab="Pielou evenness")

dev.off()

#### Alpha Diversity Visualization -- based on Rarefied Counts ####

t1r<-ggtexttable(compare_means(Shan_div.r~Ori_Ster,data=otu.rar_meta))
#t1.ar<-ggtexttable(compare_means(Shan_div.r~Ori_Ster,data=otu.rar_meta,method="t.test"))
#t1.br<-ggtexttable(compare_means(Shan_div.r~Ori_Ster,data=otu.rar_meta,method="anova"))
table1r<-ggarrange(t1r, t1.ar, t1.br, legend = "none", common.legend = FALSE, labels = c("Wilcoxon Test", "T-test", "ANOVA"),align="hv", ncol = 1, nrow = 3)
ggsave(table1r,filename = "16S_RAREFIED_ShannonDiversity_Comparisons_OriginxSterility_5.8.21.pdf", width=15, height=15, dpi=600)

t2r<-ggtexttable(compare_means(SR~Ori_Ster,data=otu.rar_meta))
t2.ar<-ggtexttable(compare_means(SR~Ori_Ster,data=otu.rar_meta,method="t.test"))
t2.br<-ggtexttable(compare_means(SR~Ori_Ster,data=otu.rar_meta,method="anova"))
table2r<-ggarrange(t2r, t2.ar, t2.br, legend = "none", common.legend = FALSE, labels = c("Wilcoxon Test", "T-test", "ANOVA"),align="hv", ncol = 1, nrow = 3)
ggsave(table2r,filename = "16S_RAREFIED_SpeciesRichness_Comparisons_OriginxSterility_5.8.21.pdf", width=15, height=15, dpi=600)

f_nobars<-ggplot(otu_meta, aes(x=Ori_Ster, y=SR, fill=Ori_Ster)) + geom_boxplot(color="black") +labs(title = "Species Richness - Treatment Condition by Source", x="Treatment Condition by Source", y="Species Richness") + theme_classic() + scale_x_discrete(name ="Treatment Condition by Source", labels=c("lab.N"="Lab x Natural", "wild.N"="Garden x Natural", "lab.S"="Lab x Sterile", "wild.S"="Garden x Sterile")) +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15))+
  scale_fill_manual(values = saturation(bright_pal, 0.8), name ="Treatment Condition by Source", labels=c("lab.N"="Lab x Natural", "wild.N"="Garden x Natural", "lab.S"="Lab x Sterile", "wild.S"="Garden x Sterile"))
ggsave(f_nobars,filename = "16S_RAREFIED_NoSigBars_SpeciesRichness_ConditionXSource_Denise_1.22.21.pdf", width=15, height=10, dpi=600)

f7r<-ggplot(otu.rar_meta, aes(x=Ori_Ster, y=SR, fill=Ori_Ster)) + geom_boxplot(color="black") +labs(title = "Species Richness - Treatment Condition by Source", x="Treatment Condition by Source", y="Species Richness") + theme_classic() + scale_x_discrete(name ="Treatment Condition by Source", labels=c("lab.N"="Lab x Natural", "wild.N"="Garden x Natural", "lab.S"="Lab x Sterile", "wild.S"="Garden x Sterile")) +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15)) +
  scale_fill_manual(values = saturation(bright_pal, 0.8), name ="Treatment Condition by Source", labels=c("lab.N"="Lab x Natural", "wild.N"="Garden x Natural", "lab.S"="Lab x Sterile", "wild.S"="Garden x Sterile"))+stat_compare_means(comparisons = list(c(1,2), c(3,4), c(2,3), c(2,4), c(1,3), c(1,4)), hide.ns = FALSE,label = "p.signif")
ggsave(f7r,filename = "16S_RAREFIED_SpeciesRichness_ALL_ConditionXSource_Denise_5.8.21.pdf", width=15, height=10, dpi=600)

f8r<-ggplot(otu.rar_meta, aes(x=Ori_Ster, y=SR, fill=Ori_Ster)) + geom_boxplot(color="black") +labs(title = "Species Richness - Treatment Condition by Source", x="Treatment Condition by Source", y="Species Richness") + theme_classic() + scale_x_discrete(name ="Treatment Condition by Source", labels=c("lab.N"="Lab x Natural", "wild.N"="Garden x Natural", "lab.S"="Lab x Sterile", "wild.S"="Garden x Sterile")) +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15)) +
  scale_fill_manual(values = saturation(bright_pal, 0.8), name ="Treatment Condition by Source", labels=c("lab.N"="Lab x Natural", "wild.N"="Garden x Natural", "lab.S"="Lab x Sterile", "wild.S"="Garden x Sterile"))+stat_compare_means(comparisons = list(c(1,2), c(1,4)), hide.ns = FALSE,label = "p.signif")
ggsave(f8r,filename = "16S_RAREFIED_SpeciesRichness_SIGNIFICANT_ConditionXSource_Denise_5.8.21.pdf", width=15, height=10, dpi=600)


f_nobars2<-ggplot(otu.rar_meta, aes(x=Ori_Ster, y=Shan_div.r, fill=Ori_Ster)) + geom_boxplot(color="black") +labs(title = "Shannon Diversity - Treatment Condition by Source", x="Treatment Condition by Source", y="Shannon Diversity") + theme_classic() + scale_x_discrete(name ="Treatment Condition by Source", labels=c("lab.N"="Lab x Natural", "wild.N"="Garden x Natural", "lab.S"="Lab x Sterile", "wild.S"="Garden x Sterile")) +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15)) +
  scale_fill_manual(values = saturation(bright_pal, 0.8), name ="Treatment Condition by Source", labels=c("lab.N"="Lab x Natural", "lab.S"="Lab x Sterile", "wild.N"="Garden x Natural", "wild.S"="Garden x Sterile"))
ggsave(f_nobars2,filename = "16S_RAREFIED_NoSigBars_ShannonDiversity_ConditionXSource_Denise_1.22.21.pdf", width=15, height=10, dpi=600)

f9r<-ggplot(otu.rar_meta, aes(x=Ori_Ster, y=Shan_div.r, fill=Ori_Ster)) + geom_boxplot(color="black") +labs(title = "Shannon Diversity - Treatment Condition by Source", x="Treatment Condition by Source", y="Shannon Diversity") + theme_classic() + scale_x_discrete(name ="Treatment Condition by Source", labels=c("lab.N"="Lab x Natural", "wild.N"="Garden x Natural", "lab.S"="Lab x Sterile", "wild.S"="Garden x Sterile")) +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15)) +
  scale_fill_manual(values = saturation(bright_pal, 0.8), name ="Treatment Condition by Source", labels=c("lab.N"="Lab x Natural", "lab.S"="Lab x Sterile", "wild.N"="Garden x Natural", "wild.S"="Garden x Sterile"))+stat_compare_means(comparisons = list(c(1,2), c(3,4), c(2,3), c(2,4), c(1,3), c(1,4)), hide.ns = FALSE,label = "p.signif")
ggsave(f9r,filename = "16S_RAREFIED_ShannonDiversity_ALL_ConditionXSource_Denise_5.8.21.pdf", width=15, height=10, dpi=600)

f10r<-ggplot(otu.rar_meta, aes(x=Ori_Ster, y=Shan_div.r, fill=Ori_Ster)) + geom_boxplot(color="black") +labs(title = "Shannon Diversity - Treatment Condition by Source", x="Treatment Condition by Source", y="Shannon Diversity") + theme_classic() + scale_x_discrete(name ="Treatment Condition by Source", labels=c("lab.N"="Lab x Natural", "wild.N"="Garden x Natural", "lab.S"="Lab x Sterile", "wild.S"="Garden x Sterile")) +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15)) +
  scale_fill_manual(values = saturation(bright_pal, 0.8), name ="Treatment Condition by Source", labels=c("lab.N"="Lab x Natural", "lab.S"="Lab x Sterile", "wild.N"="Garden x Natural", "wild.S"="Garden x Sterile"))+stat_compare_means(comparisons = list(c(1,3), c(1,4)), hide.ns = FALSE,label = "p.signif")
ggsave(f10r,filename = "16S_RAREFIED_ShannonDiversity_SIGNIFICANT_ConditionXSource_Denise_5.8.21.pdf", width=15, height=10, dpi=600)



