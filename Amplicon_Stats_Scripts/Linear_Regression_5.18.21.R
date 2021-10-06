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
library(AER)

getwd()

#### Import & reformat the data ####

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
dim(its1_otu_table)

#### Import & format metadata ####
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
metadata<-subset(metadata, select=-c(BarcodeSequence, LinkerPrimerSequence, NdEp, A_87Sr86Sr, SiteCode, RepNum, SiteRep, SampType, TubeID, DateCode, Description)) # subset only part of the metadata we need
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

#### Rarefaction ####
# RAREFACTION in vegan: ROWS need to be SITES/samples; COLUMNS are SPECIES (OTUs, ASVs)
# see "Numerical Ecology with R", page 13-14
min_its1<-min(rowSums(its1_otu_table))
min_its1

its1_rar<-rrarefy(its1_otu_table,min_its1) 
head(its1_rar)
max(rowSums(its1_rar)) # should be equal to the min used for rarefaction
all(min_its1 == max(rowSums(its1_rar)))

min_16s<-min(rowSums(bac_otu_table))
min_16s

bac_rar<-rrarefy(bac_otu_table,min_16s) 
head(bac_rar)
max(rowSums(bac_rar))
all(min_16s == max(rowSums(bac_rar)))


#### ITS1 - Alpha Diversity + Species Richness (Rarefied data) ####
# ## H <- diversity() -- Shannon entropy function
## exp(H) -- Shannon diversity (exp of entropy)

# ** using rarefied data
Shan_ent.its1<-vegan::diversity(its1_rar, index="shannon") # Shannon entropy
Shan_div.its1<- exp(Shan_ent.its1) # Shannon Diversity aka Hill number 1
div_its1<-data.frame(ITS1_Shannon_Entropy=Shan_ent.its1,ITS1_Shannon_Diversity=Shan_div.its1)
class(div_its1)
div_its1$SampleID<-rownames(div_its1)

specnumber(its1_otu_table) # finds # of species per sample using RAW count data; if MARGIN = 2 it finds frequencies of species
# vvvv S_its1 # tells you number of OTUs in every sample
S_its1<-data.frame(ITS1_Species_Richness=specnumber(its1_otu_table), SampleID=rownames(its1_otu_table)) # finds # of species per sample using RAW count data; if MARGIN = 2 it finds frequencies of species
d_r_its1<-merge(div_its1, S_its1, by.x="SampleID", by.y="SampleID")

S.freq_its1<-specnumber(its1_otu_table, MARGIN = 2) # # finds how many times each ASV appeared in samples (frequency)
S.freq_its1 # # of OTUs that appear across ALL sites/samples

#S_its1a<-data.frame(its1_Species_Richness=specnumber(its1_binary_table)) # finds # of species per sample using RAW count data; if MARGIN = 2 it finds frequencies of species
S_its1 # tells you number of OTUs in every sample


#### 16S - Alpha Diversity + Species Richness (Rarefied data) ####
# ## H <- diversity() -- Shannon entropy function
## exp(H) -- Shannon diversity (exp of entropy)

Shan_ent.16s<-vegan::diversity(bac_rar, index="shannon") # Shannon entropy
Shan_div.16s<- exp(Shan_ent.16s) # Shannon Diversity aka Hill number 1
div_16s<-data.frame(Bac_Shannon_Entropy=Shan_ent.16s,Bac_Shannon_Diversity=Shan_div.16s)
class(div_16s)
div_16s$SampleID<-rownames(div_16s)

S_16s<-data.frame(Bac_Species_Richness=specnumber(bac_otu_table), SampleID=rownames(bac_otu_table)) # finds # of species per sample using RAW count data; if MARGIN = 2 it finds frequencies of species
# ^^^ tells you number of OTUs in every sample
d_r_16s<-merge(div_16s, S_16s, by.x="SampleID", by.y="SampleID")

S.freq_16s<-specnumber(bac_otu_table, MARGIN = 2) # # finds how many times each ASV appeared in samples (frequency)
S.freq_16s # # of OTUs that appear across ALL sites/samples

#S_16s<-data.frame(bac_Species_Richness=specnumber(bac_binary_table)) # finds # of species per sample using RAW count data; if MARGIN = 2 it finds frequencies of species

#### Create dataframe to hold all data to be compaired pairwise ####

div_dat_all<-merge(d_r_its1,d_r_16s, by.x="SampleID", by.y="SampleID")
head(div_dat_all)
pw_comp<-merge(div_dat_all, metadata, by.x="SampleID", by.y="SampleID")
head(pw_comp)

pw_comp$Elevation2<-factor(pw_comp$Elevation, levels = c("400", "1100", "2000", "2700"))
rownames(pw_comp)<-pw_comp$SampleID

#### Histogram to determine distribution of data ####
#png(file = "myplot.png")
hist(pw_comp$Elevation, freq=FALSE)
hist(pw_comp$DustComplexity, freq=FALSE)

hist(pw_comp$ITS1_Shannon_Diversity, freq=FALSE)
hist(pw_comp$ITS1_Species_Richness, freq=FALSE)

hist(pw_comp$Bac_Shannon_Diversity, freq=FALSE)
hist(pw_comp$Bac_Species_Richness, freq=FALSE)

#### Subset comp dataframe by year ####

pw_2014<-subset(pw_comp, pw_comp$Year=="2014")
pw_2015<-subset(pw_comp, pw_comp$Year=="2015")

#### Histogram to determine distribution - by year ####
#png(file = "myplot.png")
### 2014
hist(pw_2014$Elevation, freq=FALSE)
hist(pw_2014$DustComplexity, freq=FALSE)

hist(pw_2014$ITS1_Shannon_Diversity, freq=FALSE)
hist(pw_2014$ITS1_Species_Richness, freq=FALSE)

hist(pw_2014$Bac_Shannon_Diversity, freq=FALSE)
hist(pw_2014$Bac_Species_Richness, freq=FALSE)

### 2015
hist(pw_2015$Elevation, freq=FALSE)
hist(pw_2015$DustComplexity, freq=FALSE)

hist(pw_2015$ITS1_Shannon_Diversity, freq=FALSE)
hist(pw_2015$ITS1_Species_Richness, freq=FALSE)

hist(pw_2015$Bac_Shannon_Diversity, freq=FALSE)
hist(pw_2015$Bac_Species_Richness, freq=FALSE)


#### Color Palettes ####

fair_cols <- c("#BF1B0B", "#FFC465", "#66ADE5", "#252A52")
names(fair_cols) <- letters[1:4]
fair_ramp <- scales::colour_ramp(fair_cols)
fair_sat <- saturation(fair_ramp, 1)

#### Are our variables normally distributed? ####
## Using Shapiro-Wilk test for normality
# The null-hypothesis of this test is that the population is normally distributed. 
# Thus, if the p value is < the chosen alpha level, then the null hypothesis is rejected and there is evidence that the data tested are not normally distributed. 
# If the p value is > than the chosen alpha level, then the null hypothesis (that the data came from a normally distributed population) can not be rejected 
# (e.g., for an alpha level of .05, a data set with a p value of less than .05 rejects the null hypothesis that the data are from a normally distributed population
# p-value > 0.05 implying that the distribution of the data are not significantly different from normal distribution. In other words, we can assume the normality

## Visualizing with Normal Q-Q Plots: orders data from smallest to largest and plots it vs. expected normal order statistics

shapiro.test(pw_comp$DustComplexity) # p-value = 0.1778
qqnorm(pw_comp$DustComplexity, pch = 1, frame = FALSE)
qqline(pw_comp$DustComplexity, col = "steelblue", lwd = 2)

#shapiro.test(pw_comp$Elevation) -- not normally distributed, categorical variable
shapiro.test(pw_comp$ITS1_Shannon_Diversity) # p-value = 9.155e-06 -- significantly different from normal
qqnorm(pw_comp$ITS1_Shannon_Diversity, pch = 1, frame = FALSE)
qqline(pw_comp$ITS1_Shannon_Diversity, col = "steelblue", lwd = 2)

shapiro.test(pw_comp$ITS1_Species_Richness) # p-value = 0.2861
qqnorm(pw_comp$ITS1_Species_Richness, pch = 1, frame = FALSE)
qqline(pw_comp$ITS1_Species_Richness, col = "steelblue", lwd = 2)

shapiro.test(pw_comp$Bac_Shannon_Diversity) # p-value = 7.466e-06 -- significantly different from normal
qqnorm(pw_comp$Bac_Shannon_Diversity, pch = 1, frame = FALSE)
qqline(pw_comp$Bac_Shannon_Diversity, col = "steelblue", lwd = 2)

shapiro.test(pw_comp$Bac_Species_Richness) # p-value = 0.4302
qqnorm(pw_comp$Bac_Species_Richness, pch = 1, frame = FALSE)
qqline(pw_comp$Bac_Species_Richness, col = "steelblue", lwd = 2)

#### Linear Regression Comparisons ####
## here the focus is comparing dust complexity to alpha diversity, species richness, & elevation
head(pw_comp)

# Dust Complexity vs Elevation

fit1<-aov(DustComplexity ~ Elevation2, data=pw_comp)
summary(fit1)
TukeyHSD(fit1, which="Elevation2")

plot(DustComplexity ~ Elevation, data=pw_comp)
abline(fit1)

# fit.0<-aov(DustComplexity ~ as.factor(Elevation), data=pw_comp)
# summary(fit.0)
# TukeyHSD(fit.0)
# Levene's test with one independent variable
fligner.test(DustComplexity ~ Elevation, data = pw_comp)
# Levenes Test for Homogeneity of Variance
#        Df  F value  Pr(>F)
# group  3   2.3415   0.0818
# Which shows that the data do not deviate significantly from homogeneity.
elev<-pw_comp$Elevation2
compare_means(DustComplexity ~ Elevation2, data=pw_comp, group.by=elev, method="anova")

fit.test<-ggplot(pw_comp, aes(x = Elevation2, y = DustComplexity, fill=Elevation2)) + 
  geom_boxplot() + theme_classic() + guides(fill = guide_legend(reverse=TRUE)) + 
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11)) +
  labs(title="Dust Complexity x Elevation",fill="Elevation (ft)")+ylab("Dust Complexity")+xlab("Elevation (ft)")+scale_fill_manual(values=saturation(fair_cols, 0.9))+stat_compare_means(method = "anova",label.y=1.5) +stat_compare_means(comparisons = list(c(1,2), c(2,3),  c(1,3),  c(3,4),  c(2,4),  c(1,4)), method="t.test", hide.ns = TRUE,label = "p.format")

ggsave(fit.test,filename = "figures/DustComp_by_Elevation_ALL_sigbars_5.24.21.pdf", width=10, height=8, dpi=600)

fit.test0<-ggplot(pw_comp, aes(x = Elevation2, y = DustComplexity, fill=Elevation2)) + 
  geom_boxplot() + theme_classic() + guides(fill = guide_legend(reverse=TRUE)) + 
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11)) +
  labs(title="Dust Complexity x Elevation",fill="Elevation (ft)")+ylab("Dust Complexity")+xlab("Elevation (ft)")+scale_fill_manual(values=saturation(fair_cols, 0.9))+stat_compare_means(method = "anova",label.y=1.5) +stat_compare_means(comparisons = list(c(1,2), c(2,3),  c(1,3),  c(3,4),  c(2,4),  c(1,4)), method="t.test", hide.ns = TRUE,label = "p.signif")

ggsave(fit.test,filename = "figures/DustComp_by_Elevation_ALL_sigbars_5.24.21.pdf", width=10, height=8, dpi=600)

fit.testa<-ggplot(pw_comp, aes(x = as.factor(Elevation), y = DustComplexity, fill=as.factor(Elevation))) + 
  geom_boxplot() + theme_classic() + guides(fill = guide_legend(reverse=TRUE)) + 
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11)) +
  labs(title="Dust Complexity x Elevation",fill="Elevation (ft)")+ylab("Dust Complexity")+xlab("Elevation (ft)")+scale_fill_manual(values=saturation(fair_cols, 0.9))

ggsave(fit.testa,filename = "figures/DustComp_by_Elevation_ALL_no.sigbars_5.24.21.pdf", width=10, height=8, dpi=600)

fit.testb<-ggplot(pw_comp, aes(x = as.factor(Elevation), y = DustComplexity, fill=as.factor(Elevation))) + 
  geom_boxplot() + theme_classic() + guides(fill = guide_legend(reverse=TRUE)) + 
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11)) +
  labs(title="Dust Complexity x Elevation",fill="Elevation (ft)")+ylab("Dust Complexity")+xlab("Elevation (ft)")+scale_fill_grey(start=0.8, end=0.3)

ggsave(fit.testb,filename = "figures/DustComp_by_Elevation_ALL_gray_5.24.21.pdf", width=10, height=8, dpi=600)

fit.testb.0<-ggplot(pw_comp, aes(x = Elevation2, y = DustComplexity, fill=Elevation2)) + 
  geom_boxplot() + theme_classic() + guides(fill = guide_legend(reverse=TRUE)) + 
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11)) +
  labs(title="Dust Complexity x Elevation",fill="Elevation (ft)")+ylab("Dust Complexity")+xlab("Elevation (ft)")+scale_fill_grey(start=0.8, end=0.3)+stat_compare_means(method = "anova",label.y=1.5) +stat_compare_means(comparisons = list(c(1,2), c(2,3),  c(1,3),  c(3,4),  c(2,4),  c(1,4)), method="t.test", hide.ns = TRUE,label = "p.signif")

ggsave(fit.testb.0,filename = "figures/DustComp_by_Elevation_ALL_gray_sigbars_5.24.21.pdf", width=10, height=8, dpi=600)

### Fungi comparisons first
# Dust Comp x ITS1 Shannon diversity
its1.fit1<-glm(DustComplexity ~ ITS1_Shannon_Diversity, data=pw_comp, family=poisson)
summary(its1.fit1)
dispersiontest(its1.fit1) # null hypothesis is that equidispersion exists; alternative hypothesis is overdispersion
# p val is 1, dispersion = 0.06...going to go with no overdispersion

plot(DustComplexity ~ ITS1_Shannon_Diversity, data=pw_comp)
abline(its1.fit1)

leveneTest(pw_comp$DustComplexity,
            pw_comp$ITS1_Shannon_Diversity,
            location = c("median"),
            trim.alpha = 0.25)
# Levenes Test for Homogeneity of Variance
#        Df  F value  Pr(>F)
# group  3   2.3415   0.0818
# Which shows that the data do not deviate significantly from homogeneity.

fig.its1.fit1<-ggplot(pw_comp, aes(x = ITS1_Shannon_Diversity, y = DustComplexity)) + 
  geom_point(aes(color=Elevation), size=3) + theme_classic() + saturation(scale_colour_gradientn(colours=fair_cols), 0.9) +
  stat_smooth(method = "lm", col = "black", se=FALSE, size=1)+ labs(title="Dust Complexity x ITS1 Shannon Diversity", color="Elevation (ft)")+ylab("Dust Complexity")+xlab("ITS1 Shannon Diversity")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))

ggsave(fig.its1.fit1,filename = "figures/DustComp_by_ITS1_Shan_Div_ALL_test_5.19.21.pdf", width=10, height=8, dpi=600)


# DustComp x ITS1 Species Richness
its1.fit2<-glm(DustComplexity ~ ITS1_Species_Richness, data=pw_comp, family=poisson)
summary(its1.fit2)

plot(DustComplexity ~ ITS1_Species_Richness, data=pw_comp)
abline(its1.fit2)

fig.its1.fit2<-ggplot(pw_comp, aes(x = ITS1_Species_Richness, y = DustComplexity)) + 
  geom_point(aes(color=Elevation), size=3) + theme_classic() + saturation(scale_colour_gradientn(colours=fair_cols), 0.9) +
  stat_smooth(method = "lm", col = "black", se=FALSE, size=1)+ labs(title="Dust Complexity x ITS1 Species Richness", color="Elevation (ft)")+ylab("Dust Complexity")+xlab("ITS1 Species Richness")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))

ggsave(fig.its1.fit2,filename = "figures/DustComp_by_ITS1_Spec_Richness_ALL_test_5.24.21.pdf", width=10, height=8, dpi=600)

### Bacteria next

# Dust Comp x Bac Shannon diversity
Bac.fit1<-lm(DustComplexity ~ Bac_Shannon_Diversity, data=pw_comp)
summary(Bac.fit1)

plot(DustComplexity ~ Bac_Shannon_Diversity, data=pw_comp)
abline(Bac.fit1)

fig.Bac.fit1<-ggplot(pw_comp, aes(x = Bac_Shannon_Diversity, y = DustComplexity, color=Elevation)) + 
  geom_point(size=3) + theme_classic() + saturation(scale_colour_gradientn(colours=fair_cols), 0.9) +
  stat_smooth(method = "lm", col = "black", se=FALSE, size=1)+ labs(title="Dust Complexity x 16S Shannon Diversity", color="Elevation (ft)")+ylab("Dust Complexity")+xlab("16S Shannon Diversity")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))

ggsave(fig.Bac.fit1,filename = "figures/DustComp_by_16S_Shan_Div_ALL_test_5.19.21.pdf", width=10, height=8, dpi=600)

### wanted to see what regression line looked at without outliers
mean(pw_comp[-7,]$Bac_Shannon_Diversity)
pw_comp.bacclean<-pw_comp[(pw_comp$Bac_Shannon_Diversity) < 150,]
max(pw_comp.bacclean$Bac_Shannon_Diversity)

ggplot(pw_comp.bacclean, aes(x = Bac_Shannon_Diversity, y = DustComplexity, color=Elevation)) + 
  geom_point(size=3) + theme_classic() + saturation(scale_colour_gradientn(colours=fair_cols), 0.9) +
  stat_smooth(method = "lm", col = "black", se=FALSE, size=1)+ labs(title="Dust Complexity x 16S Shannon Diversity", color="Elevation (ft)")+ylab("Dust Complexity")+xlab("16S Shannon Diversity")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))

# DustComp x Bac Species Richness
Bac.fit2<-lm(DustComplexity ~ Bac_Species_Richness, data=pw_comp)
summary(Bac.fit2)

plot(DustComplexity ~ Bac_Species_Richness, data=pw_comp)
abline(Bac.fit2)

fig.Bac.fit2<-ggplot(pw_comp, aes(x = Bac_Species_Richness, y = DustComplexity, color=Elevation)) + 
  geom_point(size=3) + theme_classic() + saturation(scale_colour_gradientn(colours=fair_cols), 0.9) +
  stat_smooth(method = "lm", col = "black", se=FALSE, size=1)+ labs(title="Dust Complexity x 16S Species Richness", color="Elevation (ft)")+ylab("Dust Complexity")+xlab("16S Species Richness")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))

ggsave(fig.Bac.fit2,filename = "figures/DustComp_by_16S_Spec_Richness_ALL_test_5.24.21.pdf", width=10, height=8, dpi=600)



#### Model Validation ####

# For the fitted model, validation requires verification of:
# 1. Homogeneity of variance.
# 2. Model misfit.
# 3. Normality of residuals.
# 4. Absence of influential observations.

# 1. Homogeneity of variance
## plotting model residual variance (the variance in the response variable that is not explained by the model) against model fitted values
##Ideally, the distribution of residuals around zero should be consistent along the horizontal axis.
Fitted <- fitted(fit1)
Resid  <- resid(fit1, type = "pearson")
par(mfrow = c(1,1), mar = c(5,5,2,2))
plot(x = Fitted, y = Resid,
     xlab = "Fitted values",
     ylab = "Pearson Residuals")
abline(h = 0, lty = 2)

# 2. Model Misfit
# Model misfit occurs if key covariates (including interactions) are missing from the model, or the model departs from linearity. 
# Model misfit can be recognised visually by plotting Pearson residuals against each covariate in the model, as well as those not included in the model
# distribution of residuals around zero should be consistent along the horizontal axis
plot(x = pw_comp$DustComplexity,
     y = Resid,
     xlab = "Dust Complexity",
     ylab = "Pearson residuals",
     pch = 16, cex = 1.5)
abline(h = 0, lty = 2)

# 3. Normality of Residuals
# are the residuals normally distributed?
histogram(Resid)

# 4. Absence of influential observations
# plotting Cook's distamnce
# Cook's distance is estimated by systemically dropping each observation and comparing hte fitted values with those when all observations are included in model
# Cook's distance > 1 indicates an influential data point

par(mfrow = c(1, 1))
plot(cooks.distance(fit1),
       xlab = "Observation",
       ylab = "Cook's distance",
       type = "h",
       ylim = c(0, 1.1),
       cex.lab =  1.5)
abline(h = 1, lty = 2)

#### Linear Regression Comparisons 2014 ####
## here the focus is comparing dust complexity to alpha diversity, species richness, & elevation
head(pw_2014)

# Dust Complexity vs Elevation

fit1.14<-aov(DustComplexity ~ Elevation, data=pw_2014)
summary(fit1.14)

plot(DustComplexity ~ Elevation, data=pw_2014)
abline(fit1.14)

fit.14<-ggplot(pw_2014, aes(x = Elevation2, y = DustComplexity, fill=Elevation2)) + 
  geom_boxplot() + theme_classic() + guides(fill = guide_legend(reverse=TRUE)) + 
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11)) +
  labs(title="Dust Complexity x Elevation (2014)",fill="Elevation (ft)")+ylab("Dust Complexity")+xlab("Elevation (ft)")+scale_fill_manual(values=saturation(fair_cols, 0.9))+stat_compare_means(method = "anova",label.y=1.5) +stat_compare_means(comparisons = list(c(1,2), c(2,3),  c(1,3),  c(3,4),  c(2,4),  c(1,4)), method="t.test", hide.ns = TRUE,label = "p.format")

ggsave(fit.14,filename = "figures/DustComp_by_Elevation_2014_pvals_sigbars_5.24.21.pdf", width=10, height=8, dpi=600)

fit.14.0<-ggplot(pw_2014, aes(x = Elevation2, y = DustComplexity, fill=Elevation2)) + 
  geom_boxplot() + theme_classic() + guides(fill = guide_legend(reverse=TRUE)) + 
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11)) +
  labs(title="Dust Complexity x Elevation (2014)",fill="Elevation (ft)")+ylab("Dust Complexity")+xlab("Elevation (ft)")+scale_fill_manual(values=saturation(fair_cols, 0.9))+stat_compare_means(method = "anova",label.y=1.5) +stat_compare_means(comparisons = list(c(1,2), c(2,3),  c(1,3),  c(3,4),  c(2,4),  c(1,4)), method="t.test", hide.ns = TRUE,label = "p.signif")

ggsave(fit.14.0,filename = "figures/DustComp_by_Elevation_2014_sigbars_5.24.21.pdf", width=10, height=8, dpi=600)

fit1.14a<-ggplot(pw_2014, aes(x = as.factor(Elevation), y = DustComplexity, fill=as.factor(Elevation))) + 
  geom_boxplot() + theme_classic() + guides(fill = guide_legend(reverse=TRUE)) + 
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11)) +
  labs(title="Dust Complexity x Elevation (2014)",fill="Elevation (ft)")+ylab("Dust Complexity")+xlab("Elevation (ft)")+scale_fill_manual(values=saturation(fair_cols, 0.9))

ggsave(fit1.14a,filename = "figures/DustComp_by_Elevation_2014_no.sigbars_5.24.21.pdf", width=10, height=8, dpi=600)

fit1.14b<-ggplot(pw_2014, aes(x = as.factor(Elevation), y = DustComplexity, fill=as.factor(Elevation))) + 
  geom_boxplot() + theme_classic() + guides(fill = guide_legend(reverse=TRUE)) + 
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11)) +
  labs(title="Dust Complexity x Elevation (2014)",fill="Elevation (ft)")+ylab("Dust Complexity")+xlab("Elevation (ft)")+scale_fill_grey(start=0.8, end=0.3)

ggsave(fit1.14b,filename = "figures/DustComp_by_Elevation_2014_gray_5.24.21.pdf", width=10, height=8, dpi=600)

fit.14b.0<-ggplot(pw_2014, aes(x = Elevation2, y = DustComplexity, fill=Elevation2)) + 
  geom_boxplot() + theme_classic() + guides(fill = guide_legend(reverse=TRUE)) + 
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11)) +
  labs(title="Dust Complexity x Elevation (2014)",fill="Elevation (ft)")+ylab("Dust Complexity")+xlab("Elevation (ft)")+scale_fill_grey(start=0.8, end=0.3)+stat_compare_means(method = "anova",label.y=1.5) +stat_compare_means(comparisons = list(c(1,2), c(2,3),  c(1,3),  c(3,4),  c(2,4),  c(1,4)), method="t.test", hide.ns = TRUE,label = "p.signif")

ggsave(fit.14b.0,filename = "figures/DustComp_by_Elevation_2014_gray_sigbars_5.24.21.pdf", width=10, height=8, dpi=600)

### Fungi comparisons first
# Dust Comp x ITS1 Shannon diversity
its1.fit1.14<-lm(DustComplexity ~ ITS1_Shannon_Diversity, data=pw_2014)
summary(its1.fit1.14)

plot(DustComplexity ~ ITS1_Shannon_Diversity, data=pw_2014)
abline(its1.fit1.14)

fig.its1.fit1.14<-ggplot(pw_2014, aes(x = ITS1_Shannon_Diversity, y = DustComplexity)) + 
  geom_point(aes(color=Elevation), size=3) + theme_classic() + saturation(scale_colour_gradientn(colours=fair_cols), 0.9) +
  stat_smooth(method = "lm", col = "black", se=FALSE, size=1)+ labs(title="Dust Complexity x ITS1 Shannon Diversity (2014)", color="Elevation (ft)")+ylab("Dust Complexity")+xlab("ITS1 Shannon Diversity")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))

ggsave(fig.its1.fit1.14,filename = "figures/DustComp_by_ITS1_Shan_Div_2014_5.24.21.pdf", width=10, height=8, dpi=600)


# DustComp x ITS1 Species Richness
its1.fit2.14<-lm(DustComplexity ~ ITS1_Species_Richness, data=pw_2014)
summary(its1.fit2.14)

plot(DustComplexity ~ ITS1_Species_Richness, data=pw_2014)
abline(its1.fit2.14)

fig.its1.fit2.14<-ggplot(pw_2014, aes(x = ITS1_Species_Richness, y = DustComplexity)) + 
  geom_point(aes(color=Elevation), size=3) + theme_classic() + saturation(scale_colour_gradientn(colours=fair_cols), 0.9) +
  stat_smooth(method = "lm", col = "black", se=FALSE, size=1)+ labs(title="Dust Complexity x ITS1 Species Richness (2014)", color="Elevation (ft)")+ylab("Dust Complexity")+xlab("ITS1 Species Richness")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))

ggsave(fig.its1.fit2.14,filename = "figures/DustComp_by_ITS1_Spec_Richness_2014_5.24.21.pdf", width=10, height=8, dpi=600)

### Bacteria next

# Dust Comp x Bac Shannon diversity
Bac.fit1.14<-lm(DustComplexity ~ Bac_Shannon_Diversity, data=pw_2014)
summary(Bac.fit1.14)

plot(DustComplexity ~ Bac_Shannon_Diversity, data=pw_2014)
abline(Bac.fit1.14)

fig.Bac.fit1.14<-ggplot(pw_2014, aes(x = Bac_Shannon_Diversity, y = DustComplexity, color=Elevation)) + 
  geom_point(size=3) + theme_classic() + saturation(scale_colour_gradientn(colours=fair_cols), 0.9) +
  stat_smooth(method = "lm", col = "black", se=FALSE, size=1)+ labs(title="Dust Complexity x 16S Shannon Diversity (2014)", color="Elevation (ft)")+ylab("Dust Complexity")+xlab("16S Shannon Diversity")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))

ggsave(fig.Bac.fit1.14,filename = "figures/DustComp_by_16S_Shan_Div_2014_5.24.21.pdf", width=10, height=8, dpi=600)


# DustComp x Bac Species Richness
Bac.fit2.14<-lm(DustComplexity ~ Bac_Species_Richness, data=pw_2014)
summary(Bac.fit2.14)

plot(DustComplexity ~ Bac_Species_Richness, data=pw_2014)
abline(Bac.fit2.14)

fig.Bac.fit2.14<-ggplot(pw_2014, aes(x = Bac_Species_Richness, y = DustComplexity, color=Elevation)) + 
  geom_point(size=3) + theme_classic() + saturation(scale_colour_gradientn(colours=fair_cols), 0.9) +
  stat_smooth(method = "lm", col = "black", se=FALSE, size=1)+ labs(title="Dust Complexity x 16S Species Richness (2014)", color="Elevation (ft)")+ylab("Dust Complexity")+xlab("16S Species Richness")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))

ggsave(fig.Bac.fit2.14,filename = "figures/DustComp_by_16S_Spec_Richness_2014_5.24.21.pdf", width=10, height=8, dpi=600)

#### Linear Regression Comparisons 2015 ####
## here the focus is comparing dust complexity to alpha diversity, species richness, & elevation
head(pw_2015)

# Dust Complexity vs Elevation

fit1.15<-aov(DustComplexity ~ Elevation, data=pw_2015)
summary(fit1.15)

plot(DustComplexity ~ Elevation, data=pw_2015)
abline(fit1.15)

fit.15<-ggplot(pw_2015, aes(x = Elevation2, y = DustComplexity, fill=Elevation2)) + 
  geom_boxplot() + theme_classic() + guides(fill = guide_legend(reverse=TRUE)) + 
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11)) +
  labs(title="Dust Complexity x Elevation (2015)",fill="Elevation (ft)")+ylab("Dust Complexity")+xlab("Elevation (ft)")+scale_fill_manual(values=saturation(fair_cols, 0.9))+stat_compare_means(method = "anova",label.y=1.5) +stat_compare_means(comparisons = list(c(1,2), c(2,3),  c(1,3),  c(3,4),  c(2,4),  c(1,4)), method="t.test", hide.ns = TRUE,label = "p.format")

ggsave(fit.15,filename = "figures/DustComp_by_Elevation_2015_pvals_sigbars_5.24.21.pdf", width=10, height=8, dpi=600)

fit.15.0<-ggplot(pw_2015, aes(x = Elevation2, y = DustComplexity, fill=Elevation2)) + 
  geom_boxplot() + theme_classic() + guides(fill = guide_legend(reverse=TRUE)) + 
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11)) +
  labs(title="Dust Complexity x Elevation (2015)",fill="Elevation (ft)")+ylab("Dust Complexity")+xlab("Elevation (ft)")+scale_fill_manual(values=saturation(fair_cols, 0.9))+stat_compare_means(method = "anova",label.y=1.5) +stat_compare_means(comparisons = list(c(1,2), c(2,3),  c(1,3),  c(3,4),  c(2,4),  c(1,4)), method="t.test", hide.ns = TRUE,label = "p.signif")

ggsave(fit.15.0,filename = "figures/DustComp_by_Elevation_2015_sigbars_5.24.21.pdf", width=10, height=8, dpi=600)

fit1.15a<-ggplot(pw_2015, aes(x = as.factor(Elevation), y = DustComplexity, fill=as.factor(Elevation))) + 
  geom_boxplot() + theme_classic() + guides(fill = guide_legend(reverse=TRUE)) + 
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11)) +
  labs(title="Dust Complexity x Elevation (2015)",fill="Elevation (ft)")+ylab("Dust Complexity")+xlab("Elevation (ft)")+scale_fill_manual(values=saturation(fair_cols, 0.9))

ggsave(fit1.15a,filename = "figures/DustComp_by_Elevation_2015_no.sigbars_5.24.21.pdf", width=10, height=8, dpi=600)

fit1.15b<-ggplot(pw_2015, aes(x = as.factor(Elevation), y = DustComplexity, fill=as.factor(Elevation))) + 
  geom_boxplot() + theme_classic() + guides(fill = guide_legend(reverse=TRUE)) + 
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11)) +
  labs(title="Dust Complexity x Elevation (2015)",fill="Elevation (ft)")+ylab("Dust Complexity")+xlab("Elevation (ft)")+scale_fill_grey(start=0.8, end=0.3)

ggsave(fit1.15b,filename = "figures/DustComp_by_Elevation_2015_gray_5.24.21.pdf", width=10, height=8, dpi=600)

fit.15b.0<-ggplot(pw_2015, aes(x = Elevation2, y = DustComplexity, fill=Elevation2)) + 
  geom_boxplot() + theme_classic() + guides(fill = guide_legend(reverse=TRUE)) + 
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11)) +
  labs(title="Dust Complexity x Elevation (2015)",fill="Elevation (ft)")+ylab("Dust Complexity")+xlab("Elevation (ft)")+scale_fill_grey(start=0.8, end=0.3)+stat_compare_means(method = "anova",label.y=1.5) +stat_compare_means(comparisons = list(c(1,2), c(2,3),  c(1,3),  c(3,4),  c(2,4),  c(1,4)), method="t.test", hide.ns = TRUE,label = "p.signif")

ggsave(fit.15b.0,filename = "figures/DustComp_by_Elevation_2015_gray_sigbars_5.24.21.pdf", width=10, height=8, dpi=600)

### Fungi comparisons first
# Dust Comp x ITS1 Shannon diversity
its1.fit1.15<-lm(DustComplexity ~ ITS1_Shannon_Diversity, data=pw_2015)
summary(its1.fit1.15)

plot(DustComplexity ~ ITS1_Shannon_Diversity, data=pw_2015)
abline(its1.fit1.15)

fig.its1.fit1.15<-ggplot(pw_2015, aes(x = ITS1_Shannon_Diversity, y = DustComplexity, color=Elevation)) + 
  geom_point(size=3) + theme_classic() + saturation(scale_colour_gradientn(colours=fair_cols), 0.9) +
  stat_smooth(method = "lm", col = "black", se=FALSE, size=1)+ labs(title="Dust Complexity x ITS1 Shannon Diversity (2015)", color="Elevation (ft)")+ylab("Dust Complexity")+xlab("ITS1 Shannon Diversity")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))

ggsave(fig.its1.fit1.15,filename = "figures/DustComp_by_ITS1_Shan_Div_2015_5.24.21.pdf", width=10, height=8, dpi=600)


# DustComp x ITS1 Species Richness
its1.fit2.15<-lm(DustComplexity ~ ITS1_Species_Richness, data=pw_2015)
summary(its1.fit2.15)

plot(DustComplexity ~ ITS1_Species_Richness, data=pw_2015)
abline(its1.fit2.15)

fig.its1.fit2.15<-ggplot(pw_2015, aes(x = ITS1_Species_Richness, y = DustComplexity, color=Elevation)) + 
  geom_point(size=3) + theme_classic() + saturation(scale_colour_gradientn(colours=fair_cols), 0.9) +
  stat_smooth(method = "lm", col = "black", se=FALSE, size=1)+ labs(title="Dust Complexity x ITS1 Species Richness (2015)", color="Elevation (ft)")+ylab("Dust Complexity")+xlab("ITS1 Species Richness")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))

ggsave(fig.its1.fit2.15,filename = "figures/DustComp_by_ITS1_Spec_Richness_2015_5.24.21.pdf", width=10, height=8, dpi=600)

### Bacteria next

# Dust Comp x Bac Shannon diversity
Bac.fit1.15<-lm(DustComplexity ~ Bac_Shannon_Diversity, data=pw_2015)
summary(Bac.fit1.15)

plot(DustComplexity ~ Bac_Shannon_Diversity, data=pw_2015)
abline(Bac.fit1.15)

fig.Bac.fit1.15<-ggplot(pw_2015, aes(x = Bac_Shannon_Diversity, y = DustComplexity, color=Elevation)) + 
  geom_point(size=3) + theme_classic() + saturation(scale_colour_gradientn(colours=fair_cols), 0.9) +
  stat_smooth(method = "lm", col = "black", se=FALSE, size=1)+ labs(title="Dust Complexity x 16S Shannon Diversity (2015)", color="Elevation (ft)")+ylab("Dust Complexity")+xlab("16S Shannon Diversity")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))

ggsave(fig.Bac.fit1.15,filename = "figures/DustComp_by_16S_Shan_Div_2015_5.24.21.pdf", width=10, height=8, dpi=600)


# DustComp x Bac Species Richness
Bac.fit2.15<-lm(DustComplexity ~ Bac_Species_Richness, data=pw_2015)
summary(Bac.fit2.15)

plot(DustComplexity ~ Bac_Species_Richness, data=pw_2015)
abline(Bac.fit2.15)

fig.Bac.fit2.15<-ggplot(pw_2015, aes(x = Bac_Species_Richness, y = DustComplexity, color=Elevation)) + 
  geom_point(size=3) + theme_classic() + saturation(scale_colour_gradientn(colours=fair_cols), 0.9) +
  stat_smooth(method = "lm", col = "black", se=FALSE, size=1)+ labs(title="Dust Complexity x 16S Species Richness (2015)", color="Elevation (ft)")+ylab("Dust Complexity")+xlab("16S Species Richness")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))

ggsave(fig.Bac.fit2.15,filename = "figures/DustComp_by_16S_Spec_Richness_2015_5.24.21.pdf", width=10, height=8, dpi=600)




#### Visual Comparison of Regressions ####

comp1<-ggarrange(fit.14.0, fit.15.0, ncol=2, nrow=1)
ggsave(comp1,filename = "figures/DustComp_by_Elevation_Comparison_5.24.21.pdf", width=20, height=8, dpi=600)

comp2<-ggarrange(fig.its1.fit1.14, fig.its1.fit1.15, ncol=2, nrow=1)
ggsave(comp2,filename = "figures/DustComp_by_ITS1_ShanDiv_Comparison_5.24.21.pdf", width=20, height=8, dpi=600)

comp3<-ggarrange(fig.its1.fit2.14, fig.its1.fit2.15, ncol=2, nrow=1)
ggsave(comp3,filename = "figures/DustComp_by_ITS1_SpecRich_Comparison_5.24.21.pdf", width=20, height=8, dpi=600)

comp4<-ggarrange(fig.Bac.fit1.14, fig.Bac.fit1.15, ncol=2, nrow=1)
ggsave(comp4,filename = "figures/DustComp_by_16S_ShanDiv_Comparison_5.24.21.pdf", width=20, height=8, dpi=600)

comp5<-ggarrange(fig.Bac.fit2.14, fig.Bac.fit2.15, ncol=2, nrow=1)
ggsave(comp5,filename = "figures/DustComp_by_16S_SpecRich_Comparison_5.24.21.pdf", width=20, height=8, dpi=600)

## comparison of each year and combo of both years

comp1a<-ggarrange(fit.14.0, fit.15.0, fit.test0, ncol=3, nrow=1)
ggsave(comp1a,filename = "figures/DustComp_by_Elevation_Comp2_5.24.21.pdf", width=25, height=8, dpi=600)

comp2a<-ggarrange(fig.its1.fit1.14, fig.its1.fit1.15, fig.its1.fit1, ncol=3, nrow=1)
ggsave(comp2a,filename = "figures/DustComp_by_ITS1_ShanDiv_Comp2_5.24.21.pdf", width=25, height=8, dpi=600)

comp3a<-ggarrange(fig.its1.fit2.14, fig.its1.fit2.15, fig.its1.fit2, ncol=3, nrow=1)
ggsave(comp3a,filename = "figures/DustComp_by_ITS1_SpecRich_Comp2_5.24.21.pdf", width=25, height=8, dpi=600)

comp4a<-ggarrange(fig.Bac.fit1.14, fig.Bac.fit1.15, fig.Bac.fit1, ncol=3, nrow=1)
ggsave(comp4a,filename = "figures/DustComp_by_16S_ShanDiv_Comp2_5.24.21.pdf", width=25, height=8, dpi=600)

comp5a<-ggarrange(fig.Bac.fit2.14, fig.Bac.fit2.15, fig.Bac.fit2, ncol=3, nrow=1)
ggsave(comp5a,filename = "figures/DustComp_by_16S_SpecRich_Comp2_5.24.21.pdf", width=20, height=8, dpi=600)

#### Multiple Lines per plot ####
head(pw_comp)

testplot1<-ggplot(pw_comp) + 
  geom_jitter(aes(ITS1_Shannon_Diversity,DustComplexity), colour="blue") + geom_smooth(aes(ITS1_Shannon_Diversity,DustComplexity), method=lm, se=FALSE) +
  geom_jitter(aes(Bac_Shannon_Diversity,DustComplexity), colour="green") + geom_smooth(aes(Bac_Shannon_Diversity,DustComplexity), method=lm, se=FALSE) +
  labs(x = "Dust Complexity", y = "Shannon Diversity")+theme_classic()
ggsave(testplot1,filename = "figures/ShanDiv_Compare_ALL_5.24.21.pdf", width=20, height=8, dpi=600)

### Side by side comparisons ####

pw_div<-subset(pw_comp, select=c(DustComplexity, ITS1_Shannon_Diversity, Bac_Shannon_Diversity))
pw_div.m = melt(pw_div, id.vars='DustComplexity')

testplot2<-ggplot(pw_div.m) +
  geom_jitter(aes(value,DustComplexity, colour=variable),) + geom_smooth(aes(value,DustComplexity, colour=variable), method=lm, se=FALSE) +
  facet_wrap(~variable, scales="free_x") +
  labs(x = "Dust Complexity", y = "Shannon Diversity", color="Shannon Diversity")
ggsave(testplot2,filename = "figures/ShanDiv_Compare_ALL_2_5.24.21.pdf", width=20, height=8, dpi=600)
