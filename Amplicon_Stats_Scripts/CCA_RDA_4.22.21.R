#if (!requireNamespace("BiocManager", quietly = TRUE))
#install.packages("BiocManager")
#BiocManager::install(version = "3.10")
#BiocManager::install(c("phyloseq"))
library(phyloseq)
library(ggvegan)
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

# Import dust metadata
dust_comp<-read.csv("data/BinnedHiMidLow.csv", header=TRUE, sep=",")
head(dust_comp)
dust_comp<-subset(dust_comp, select=c(SampleID, DustComplexity, LoMidHi)) # subset only part of the metadata we need
rownames(dust_comp)<-dust_comp$SampleID
metadata<-merge(metadata, dust_comp, by.x="SampleID", by.y="SampleID")
head(metadata)
rownames(metadata)<-metadata$SampleID

head(metadata)
tail(metadata)
metadata<-subset(metadata, select=-c(DateCode, BarcodeSequence, LinkerPrimerSequence, SiteCode, RepNum, SiteRep, TubeID, DateCode, Description)) # subset only part of the metadata we need
head(metadata)

rownames(metadata)<-metadata$SampleID

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

meta_quant.scale<-as.data.frame(scale(meta_quant))
head(meta_quant.scale)
meta_quant.scale$SampleID<-rownames(meta_quant.scale)

site_dust<-data.frame(SampleID=metadata$SampleID, DustComplex=metadata$DustComplexity, Site=metadata$Site, DustBin=metadata$LoMidHi)
meta_qdat<-merge(meta_quant.scale, site_dust, by.x="SampleID", by.y="SampleID")
head(meta_qdat)

meta_all_scaled<-merge(meta_qdat[,-11], meta_category, by.x="SampleID", by.y="SampleID")

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

#### Subset metadata by site ####
unique(metadata$Site)

site_list<-unique(metadata$Site) #define an array of string values

#create a list of dataframes containing the subsets of the data by list
site_subsets<-lapply(site_list, function(x) {subset(metadata, Site==x)})
# here the function(x) is using site_list aka x to subset metadata, when $Site column == site_list
site_subsets # sanity check1
site_subsets[[1]] # sanity check2
#rename the list elements
names(site_subsets)<-site_list # can do this because used site_list to populate each element of list, so order is maintained
site_subsets$SJER # sanity check3 - should be able to pull dataframes by names rather than index now

site_subsets$PROVIDENCE # same as site_subsets[[2]]
site_subsets[[2]][8:10] # example of subsetting
# ^ subsetting to [[second dataframe]], [columns 8-10]
site_subsets[[2]][[6,6]] # [[second dataframe]], [[6 row, 6 column]]


df_specific.subset<-function(var_vec,var_subsets){
  # var_vec = vector of variable elements from specific categorical variable;
  ## e.g. vector of names from Site categorical variable (metadata sites)
  # var_subsets = list of dataframes subsetted by column$element from original dataframe;
  ## e.g. list of dataframes (each df = element of list) subsetted from metadata using vector of metadata$Site names
  for(i in seq_along(var_vec)){
    # print(var_vec[i]) -- var_vec[i] = each element in var_vec
    # print(var_subsets[[i]]) -- var_subsets[[i]] = each sub
    df<-paste(var_vec[i])
    #print(df)
    assign(df, var_subsets[[i]], envir = .GlobalEnv)
    print(paste("Dataframe", var_vec[i] ,"done"))

  }

}
df_specific.subset(site_list, site_subsets) # used scaled metadata quantitative values

head(SJER) # sanity check

#### Subset metadata by year ####
unique(metadata$Year)

yr_list<-unique(metadata$Year) #define an array of string values

#create a list of dataframes containing the subsets of the data by list
yr_subsets<-lapply(yr_list, function(x) {subset(metadata, Year==x)})
# here the function(x) is using site_list aka x to subset metadata, when $Site column == site_list
yr_subsets # sanity check1
yr_subsets[[1]] # sanity check2

#rename the list elements
# change names of dataframes so that you can call them easily later
## harder to call objects that have numeric names
yr_list<-gsub('(^.)', 'yr_\\1', yr_list)
names(yr_subsets)<-yr_list # can do this because used site_list to populate each element of list, so order is maintained
yr_subsets$yr_2015# sanity check3 - should be able to pull dataframes by names rather than index now
# ^ same as site_subsets[[2]]
yr_subsets[[2]][8:10] # example of subsetting
# ^ subsetting to [[second dataframe]], [columns 8-10]
yr_subsets[[2]][[6,6]] # [[second dataframe]], [[6 row, 6 column]]

df_specific.subset<-function(var_vec,var_subsets){
  # var_vec = vector of variable elements from specific categorical variable;
  ## e.g. vector of names from Site categorical variable (metadata sites)
  # var_subsets = list of dataframes subsetted by column$element from original dataframe;
  ## e.g. list of dataframes (each df = element of list) subsetted from metadata using vector of metadata$Site names
  for(i in seq_along(var_vec)){
    # print(var_vec[i]) -- var_vec[i] = each element in var_vec
    # print(var_subsets[[i]]) -- var_subsets[[i]] = each sub
    df<-paste(var_vec[i])
    #print(df)
    assign(df, var_subsets[[i]], envir = .GlobalEnv)
    print(paste("Dataframe", var_vec[i] ,"done"))

  }

}
df_specific.subset(yr_list, yr_subsets) # used scaled metadata quantitative values

head(yr_2014) # sanity check

#### Match relativized comp data to subsetted metadata (by year) ####

# some sanity checks...
head(yr_2014) # sanity check
its1_table_RA[1:3,1:3] # sanity check to see if relativized OTU table has rownames for next function
bac_table_RA[1:3,1:3]

# matching data with user defined function
#
match_dat<-function(compdata, subset_metadata){
  subset_comp_data = pullrow<-(is.element(row.names(compdata), row.names(subset_metadata)))
  ### * comp data and metadata need to have row names - rownames should be Sample IDs
  subset_comp_data=compdata[pullrow,]
  return(subset_comp_data)
}

# 2014
its1_2014<-match_dat(its1_table_RA, yr_2014)
its1_2014[1:3,1:3]
bac_2014<-match_dat(bac_table_RA, yr_2014)
bac_2014[1:3,1:3]

# 2015
its1_2015<-match_dat(its1_table_RA, yr_2015)
its1_2015[1:3,1:3]
bac_2015<-match_dat(bac_table_RA, yr_2015)
bac_2015[1:3,1:3]

#### Run DCA to determine if you should use: CCA or RDA ####
## *** ONLY USE DCA To determine if you should use a CCA or an RDA

## The length of first DCA axis:
## > 4 indicates heterogeneous dataset on which unimodal methods should be used (CCA),
##  < 3 indicates homogeneous dataset for which linear methods are suitable (RDA)
## between 3 and 4 both linear and unimodal methods are OK.

### all comp data dca
its1.dca = decorana(its1_table_RA)
plot(its1.dca)
summary (its1.dca) #DCA1 axis length = 3.0669; can use RDA or CCA

bac.dca = decorana(bac_table_RA)
plot(bac.dca)
summary (bac.dca) #DCA1 axis length = 5.0610; use CCA


### 2014 DCA
its1.2014.dca = decorana(its1_2014)
plot(its1.2014.dca)
summary (its1.2014.dca) #DCA1 axis length = 3.2023; can use RDA or CCA

bac.2014.dca = decorana(bac_2014)
plot(bac.2014.dca)
summary (bac.2014.dca) #DCA1 axis length = 4.3789; CCA

### 2015 DCA
its1.2015.dca = decorana(its1_2015)
plot(its1.2015.dca)
summary (its1.2015.dca) #DCA1 axis length = 3.1265; can use RDA or CCA

bac.2015.dca = decorana(bac_2015)
plot(bac.2015.dca)
summary (bac.2015.dca) #DCA1 axis length = 5.3095; CCA


# ### SJER dca
# sjr.its1.dca = decorana(sjer_its1)
# plot(sjr.its1.dca)
# summary (sjr.its1.dca) #DCA1 axis length = 3.7816; can use RDA or CCA
#
# sjr.bac.dca = decorana(sjer_bac)
# plot(sjr.bac.dca)
# summary (sjr.bac.dca) #DCA1 axis length = 4.7962; use CCA
#
# ### PROVIDENCE dca
# prov_its1.dca = decorana(prov_its1)
# plot(prov_its1.dca)
# summary (prov_its1.dca) #DCA1 axis length = 2.0514; use RDA
#
# prov.bac.dca = decorana(prov_bac)
# plot(prov.bac.dca)
# summary (prov.bac.dca) #DCA1 axis length = 3.3769; RDA or CCA
#
# ### SHORTHAIR dca
# sh_its1.dca = decorana(sh_its1)
# plot(sh_its1.dca)
# summary (sh_its1.dca) #DCA1 axis length = 2.3452; use RDA
#
# sh.bac.dca = decorana(sh_bac)
# plot(sh.bac.dca)
# summary (sh.bac.dca) #DCA1 axis length = 4.9442; use CCA
#
# ### SOAPROOT dca
# sr_its1.dca = decorana(sr_its1)
# plot(sr_its1.dca)
# summary (sr_its1.dca) #DCA1 axis length = 2.0693; use RDA
#
# sr.bac.dca = decorana(sr_bac)
# plot(sr.bac.dca)
# summary (sr.bac.dca) #DCA1 axis length = 3.5193; use RDA or CCA


#### Fungi CCA - across ALL sites #####

## CCA does not try to display all variation in the data,
## but only the part that can be explained by the used constraints. *******
## Consequently, the results are strongly dependent on the set of constraints
## The shotgun method is to use all environmental variables as constraints,
## however, rely on non-metric multidimensional scaling (metaMDS) and environmental
## interpretation after analysis (envfit, ordisurf).
## CCA is a good choice if the user has clear and strong a priori hypotheses on constraints

###### Fungi first!
head(metadata)
## let's start with a CCA of several environmental variables
its1.cca = cca(its1_table_RA ~  Cu + Fe + Mg + Mn + Ni + P + S + Zn + Elevation + DustComplexity + Year + Month, data = metadata )

## what do we find?
its1.cca
plot(its1.cca)

## remember, Scaling
plot(its1.cca, scaling = 1) ## this emphasizes relationships among sites
plot(its1.cca, scaling = 2) ## this emphasizes relationships among species

## what do we find?
summary(its1.cca)
its1.cca
## 1) partition of scaled chi-squared Constrained = explained by variable,
##    Unconstrained = explained by residuals
## 2) Importance of componenents - which axes contribute
## 3) Scores: LC or WA? LC give best fit if you have no noise in constraining
##    variables, WA scores good with noisy data

## need to bootstrap to get the real r-square
## how much variation does it explain?
RsquareAdj(its1.cca)
## ^^ use this b/c chance correlations can inflat R^2
## reminder: R^2 = % of variation in dependent variable explained by model

vif.cca(its1.cca)

## we can then test for significance of the model by permutation
## if it is not signfiicant,
## it doesn't matter how much of the variation is explained
## Is it significant?
anova(its1.cca, permutations = how(nperm=999))

## we can also do a permutation test by axis ********
anova(its1.cca, by = "axis", permutations = how(nperm=999)) ### by CCA axis
## or by terms
anova(its1.cca, by = "terms", permutations = how(nperm=999)) ### by variables
## this will help us interpret our CCA and we can see some variable are not significant

## we can use model selection instead of picking variables we think are important

its.cca2 = ordistep(cca(its1_table_RA ~ 1, data = metadata ),
                         scope=formula(its1.cca),
                         direction = "forward",
                         permutations = how(nperm=999))

## this tells us that we really only need a few variables
## its1_table_RA ~ Elevation + P

## so lets rerun our model with these and see how it looks

its1.cca.2 = cca(its1_table_RA ~ Elevation + P + Month, data = metadata )
summary(its1.cca.2)
anova(its1.cca.2, permutations = how(nperm=999))
# Model: cca(formula = its1_table_RA ~ Elevation + P + Month, data = metadata)
# Df ChiSquare      F Pr(>F)
# Model     4   0.95862 2.2488  0.001 ***
#   Residual 26   2.77089
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
anova(its1.cca.2, by = "terms", permutations = how(nperm=999))
# Model: cca(formula = its1_table_RA ~ Elevation + P + Month, data = metadata)
# Df ChiSquare      F Pr(>F)
# Elevation  1   0.41786 3.9209  0.001 ***
#   P          1   0.24723 2.3199  0.001 ***
#   Month      2   0.29353 1.3771  0.038 *
#   Residual  26   2.77089
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

##What about r-squared

## original model
RsquareAdj(its1.cca)
# 0.1303113

## new model
RsquareAdj(its1.cca.2)
# 0.141622

## now let rerun our model selection
its1.cca3 = ordistep(cca(its1_table_RA ~ 1, data = metadata ),
                          scope=formula(its1.cca.2),
                          direction = "forward",
                          permutations = how(nperm=999))

## so our model that was selected is
# its1 ~ Elevation + P + Month

## so how do we interpret this?
## the env explains 14.16% of the variation
## ^^ use Rsquareadj() b/c chance correlations can inflat R^2
## reminder: R^2 = % of variation in dependent variable explained by model

# Visualization w/ ggvegan
# autoplot(its1.cca)
# autoplot(its1.cca.2, axes=c(1,2),geom = c("point","text"), layers = c("sites", "biplot"),
#          legend.position = "right", title = "ITS1 CCA", ylab="CCA 2", xlab = "CCA 1")

## Note: object names are recycled here for everytime these figures are generated.
## Rerun all commands including & after fortify() command to ensure correct values are used in figures!
fits1.cca2 <- fortify(its1.cca.2, axes = 1:2)  # fortify the ordination -- convert to dataframe

take <- c('RDA1', 'RDA2')  # which columns contain the scores we want
arrows <- subset(rda.ccs, Score == 'biplot')  # take only biplot arrow scores from cca df
samps_its1 <- subset(fits1.cca2, Score == 'sites')  # take only biplot arrow scores from cca df
samp_sites<-data.frame(Label=metadata$SampleID, Site=metadata$Site, Year=metadata$Year)
samps_its1<-merge(samps_its1, samp_sites, by="Label")
mul <- ggvegan:::arrowMul(arrows[, take], subset(fits1.cca2, select = take, Score == 'sites'))
## ^ multiplier for arrows to scale them to the plot range
arrows[, take] <- arrows[, take] * mul  # scale biplot arrows w/ multiplier
arrows$Label<-gsub("Month(.*)", "\\1", arrows$Label) # change labels for months to be legible

summary(its1.cca.2)
# Importance of components:
#                        CCA1    CCA2
# Eigenvalue            0.4241 0.26116
# Proportion Explained  0.1137 0.07003
# Cumulative Proportion 0.1137 0.18375
screeplot(its1.cca.2)

its1.cca.all1<-ggplot() +
  geom_point(data = subset(fits1.cca2, Score == 'sites'),aes(x = CCA1, y = CCA2, colour=samps_its1$Site, shape=as.factor(samps_its1$Year)), size=3) +
  geom_segment(data = arrows,aes(x = 0, y = 0, xend = CCA1, yend = CCA2),
               arrow = arrow(length = unit(0.01, "npc"))) + theme_bw() + xlab("CCA1 (11.37%)") + ylab("CCA2 (7.00%)") +
  geom_label(data = arrows,aes(label = Label, x = CCA1, y = CCA2, fontface="bold",hjust = 0.5*(1 - sign(CCA1)), vjust = 0.5*(1-sign(CCA2)))) +
  scale_x_continuous(expand = c(.1, .1)) +scale_y_continuous(expand = c(.1, .1)) + scale_colour_brewer("Site", palette = "Dark2")+scale_shape_manual(name="Year",values = c(16, 17))

ggsave(its1.cca.all1,filename = "figures/ITS1_CCA_AllSites_5.17.21.pdf", width=10, height=8, dpi=600)

its1.cca.all2<-ggplot() +
  geom_point(data = subset(fits1.cca2, Score == 'sites'),aes(x = CCA1, y = CCA2, colour=samps_its1$Site, shape=as.factor(samps_its1$Year)), size=3) +
  geom_segment(data = arrows,aes(x = 0, y = 0, xend = CCA1, yend = CCA2),
               arrow = arrow(length = unit(0.01, "npc"))) + theme_bw() + xlab("CCA1 (11.37%)") + ylab("CCA2 (7.00%)") +
  geom_label(data = arrows,aes(label = Label, x = CCA1, y = CCA2, fontface="bold",hjust = 0.5*(1 - sign(CCA1)), vjust = 0.5*(1-sign(CCA2)))) +
  scale_x_continuous(expand = c(.1, .1)) +scale_y_continuous(expand = c(.1, .1)) + scale_colour_brewer("Site", palette = "Dark2") +
  geom_text(data=samps_its1, aes(x = CCA1, y = CCA2, label=Label), fontface = "plain", nudge_x = 0.1, nudge_y = 0.15)+scale_shape_manual(name="Year",values = c(16, 17))

ggsave(its1.cca.all2,filename = "figures/ITS1_CCA_AllSites_Labeled_5.17.21.pdf", width=10, height=8, dpi=600)


#### Bacteria/Archaea CCA - across ALL sites ########

head(metadata)
## let's start with a CCA of several environmental variables
bac.cca = cca(bac_table_RA ~  Cu + Fe + Mg + Mn + Ni + P + S + Zn + Elevation + DustComplexity + Year + Month, data = metadata)

## what do we find?
bac.cca
plot(bac.cca)

## remember, Scaling
plot(bac.cca, scaling = 1)
## this emphasizes relationships among sites
plot(bac.cca, scaling = 2)
## this emphasizes relationships among species

## what do we find?
summary(bac.cca)
bac.cca
## 1) partition of scaled chi-squared Constrained = explained by variable,
##    Unconstrained = explained by residuals
## 2) Importance of componenents - which axes contribute
## 3) Scores: LC or WA? LC give best fit if you have no noise in constraining
##    variables, WA scores good with noisy data

## need to bootstrap to get the real r-square
## how much variation does it explain?
RsquareAdj(bac.cca)
## ^^ use this b/c chance correlations can inflat R^2
## reminder: R^2 = % of variation in dependent variable explained by model

vif.cca(bac.cca)

## we can then test for significance of the model by permutation
## if it is not signfiicant,
## it doesn't matter how much of the variation is explained
## Is it significant?
anova(bac.cca, permutations = how(nperm=999))

## we can also do a permutation test by axis ********
anova(bac.cca, by = "axis", permutations = how(nperm=999)) ### by CCA axis
## or by terms
anova(bac.cca, by = "terms", permutations = how(nperm=999)) ### by variables
## this will help us interpret our CCA and we can see some variable are not significant

## we can use model selection instead of picking variables we think are important

bac.cca2 = ordistep(cca(bac_table_RA ~ 1, data = metadata ),
                         scope=formula(bac.cca),
                         direction = "forward",
                         permutations = how(nperm=999))

## this tells us that we really only need a few variables
## bac_table_RA ~ Elevation + Month

## so lets rerun our model with these and see how it looks

bac.cca.2 = cca(bac_table_RA ~ Elevation + Month, data = metadata )
summary(bac.cca.2)
anova(bac.cca.2, permutations = how(nperm=999))
anova(bac.cca.2, by = "terms", permutations = how(nperm=999))

##What about r-squared

## original model
RsquareAdj(bac.cca)
# 0.07258813

## new model
RsquareAdj(bac.cca.2)
# 0.05428455

## so our model that was selected is
# bac ~ Elevation + Month; maybe Zn?
## so how do we interpret this?
## the env explains 5.5% of the variation
## ^^ use Rsquareadj() b/c chance correlations can inflat R^2
## reminder: R^2 = % of variation in dependent variable explained by model

# Visualization w/ ggvegan
# autoplot(bac.cca) # original CCA
# autoplot(bac.cca.2)

fbac.cca2 <- fortify(bac.cca.2, axes = 1:2)  # fortify the ordination
## Note: object names are recycled here for everytime these figures are generated.
## Rerun all commands including & after fortify() command to ensure correct values are used in figures!
take <- c('CCA1', 'CCA2')  # which columns contain the scores we want
arrows <- subset(fbac.cca2, Score == 'biplot')  # take only biplot arrow scores from cca df
samps_bac <- subset(fbac.cca2, Score == 'sites')  # take only biplot arrow scores from cca df
samp_sites<-data.frame(Label=metadata$SampleID, Site=metadata$Site, Year=metadata$Year)
samps_bac<-merge(samps_bac, samp_sites, by="Label")
mul <- ggvegan:::arrowMul(arrows[, take], subset(fbac.cca2, select = take, Score == 'sites'))
## ^ multiplier for arrows to scale them to the plot range
arrows[, take] <- arrows[, take] * mul  # scale biplot arrows w/ multiplier
arrows$Label<-gsub("Month(.*)", "\\1", arrows$Label) # change labels for months to be legible

summary(bac.cca.2)
# Importance of components:
#                       CCA1    CCA2
# Eigenvalue            0.51652 0.34480
# Proportion Explained  0.06807 0.04544
# Cumulative Proportion 0.06807 0.11350
screeplot(bac.cca.2)

bac.cca.all1<-ggplot() +
  geom_point(data = subset(fbac.cca2, Score == 'sites'),aes(x = CCA1, y = CCA2, colour=samp_sites$Site, shape=as.factor(samp_sites$Year)), size=3) +
  geom_segment(data = arrows,aes(x = 0, y = 0, xend = CCA1, yend = CCA2),
               arrow = arrow(length = unit(0.01, "npc"))) + theme_bw() + xlab("CCA1 (6.80%)") + ylab("CCA2 (4.54%)") +
  geom_label(data = arrows,aes(label = Label, x = CCA1, y = CCA2, fontface="bold",hjust = 0.5*(1 - sign(CCA1)), vjust = 0.5*(1-sign(CCA2)))) +
  scale_x_continuous(expand = c(.1, .1)) +scale_y_continuous(expand = c(.1, .1)) + scale_colour_brewer("Site", palette = "Dark2")+scale_shape_manual(name="Year",values = c(16, 17))

ggsave(bac.cca.all1,filename = "figures/16S_CCA_AllSites_5.17.21.pdf", width=10, height=8, dpi=600)

bac.cca.all2<-ggplot() +
  geom_point(data = subset(fits1.cca2, Score == 'sites'),aes(x = CCA1, y = CCA2, colour=samp_sites$Site, shape=as.factor(samp_sites$Year)), size=3) +
  geom_segment(data = arrows,aes(x = 0, y = 0, xend = CCA1, yend = CCA2),
               arrow = arrow(length = unit(0.01, "npc"))) + theme_bw() + xlab("CCA1 (6.80%)") + ylab("CCA2 (4.54%)") +
  geom_label(data = arrows,aes(label = Label, x = CCA1, y = CCA2, fontface="bold",hjust = 0.5*(1 - sign(CCA1)), vjust = 0.5*(1-sign(CCA2)))) +
  scale_x_continuous(expand = c(.1, .1)) +scale_y_continuous(expand = c(.1, .1)) + scale_colour_brewer("Site", palette = "Dark2") +
  geom_text(data=samps_its1, aes(x = CCA1, y = CCA2, label=Label), fontface = "plain", nudge_x = 0.1, nudge_y = 0.15)+scale_shape_manual(name="Year",values = c(16, 17))

ggsave(its1.cca.all2,filename = "figures/16S_CCA_AllSites_Labeled_5.17.21.pdf", width=10, height=8, dpi=600)


#### Fungi CCA - by year #####

## CCA does not try to display all variation in the data,
## but only the part that can be explained by the used constraints. *******
## Consequently, the results are strongly dependent on the set of constraints
## The shotgun method is to use all environmental variables as constraints,
## however, rely on non-metric multidimensional scaling (metaMDS) and environmental
## interpretation after analysis (envfit, ordisurf).
## CCA is a good choice if the user has clear and strong a priori hypotheses on constraints

###### 2014 first!

head(yr_2014)
## let's start with a CCA of several environmental variables
its1.2014.cca = cca(its1_2014 ~  Cu + Fe + Mg + Mn + Ni + P + S + Zn + Elevation + DustComplexity + Month, data = yr_2014 )

## what do we find?
its1.2014.cca
plot(its1.2014.cca)

## remember, Scaling
plot(its1.2014.cca, scaling = 1) ## this emphasizes relationships among sites
plot(its1.2014.cca, scaling = 2) ## this emphasizes relationships among species

## what do we find?
summary(its1.2014.cca)
its1.2014.cca
## 1) partition of scaled chi-squared Constrained = explained by variable,
##    Unconstrained = explained by residuals
## 2) Importance of componenents - which axes contribute
## 3) Scores: LC or WA? LC give best fit if you have no noise in constraining
##    variables, WA scores good with noisy data

## need to bootstrap to get the real r-square
## how much variation does it explain?
RsquareAdj(its1.2014.cca)
## ^^ use this b/c chance correlations can inflat R^2
## reminder: R^2 = % of variation in dependent variable explained by model

vif.cca(its1.2014.cca)

## we can then test for significance of the model by permutation
## if it is not signfiicant,
## it doesn't matter how much of the variation is explained
## Is it significant?
anova(its1.2014.cca, permutations = how(nperm=999))

## we can also do a permutation test by axis ********
anova(its1.2014.cca, permutations = how(nperm=999)) ### by CCA axis
## or by terms
anova(its1.2014.cca, by = "terms", permutations = how(nperm=999)) ### by variables
## this will help us interpret our CCA and we can see some variable are not significant
# shows that Cu and Elevation are significant; S is a litle over 0.05

## we can use model selection instead of picking variables we think are important

its.2014.cca2 = ordistep(cca(its1_2014 ~ 1, data = yr_2014 ),
                    scope=formula(its1.2014.cca),
                    direction = "forward",
                    permutations = how(nperm=999))

## this tells us that we really only need a few variables
## its1_2014 + Elevation

## so lets rerun our model with these and see how it looks

its1.2014.cca.2 = cca(its1_2014 ~ Elevation, data = yr_2014 )
summary(its1.2014.cca.2)
anova(its1.2014.cca.2, permutations = how(nperm=999))
anova(its1.2014.cca.2, by = "terms", permutations = how(nperm=999))

##What about r-squared

## original model
RsquareAdj(its1.2014.cca)
# 0.9659091

## new model
RsquareAdj(its1.2014.cca.2)
# 0.05933682

## now let rerun our model selection
its.2014.cca3 = ordistep(cca(its1_2014 ~ 1, data = yr_2014 ),
                         scope=formula(its1.2014.cca.2),
                         direction = "forward",
                         permutations = how(nperm=999))
## so our model that was selected is
# its1_2014 ~ Elevation

## so how do we interpret this?
## the env explains 5.93% of the variation
## ^^ use Rsquareadj() b/c chance correlations can inflat R^2
## reminder: R^2 = % of variation in dependent variable explained by model

# Visualization w/ ggvegan
# autoplot(its1.2014.cca.2, axes=c(1,2),geom = c("point","text"), layers = c("sites", "biplot"),
#          legend.position = "right", title = "ITS1 CCA", ylab="CCA 2", xlab = "CCA 1")

fits1.2014.cca2 <- fortify(its1.2014.cca.2, axes = 1:2)  # fortify the ordination
## Note: object names are recycled here for everytime these figures are generated.
## Rerun all commands including & after fortify() command to ensure correct values are used in figures!
head(fits1.2014.cca2)
take <- c('CCA1', 'CA1')  # which columns contain the scores we want
arrows <- subset(fits1.2014.cca2, Score == 'biplot')  # take only biplot arrow scores from cca df
samps_its1.2014 <- subset(fits1.2014.cca2, Score == 'sites')  # take only biplot arrow scores from cca df
samp_sites<-data.frame(Label=metadata$SampleID, Site=metadata$Site)
samps_its1.2014<-merge(samps_its1.2014, samp_sites, by="Label")
mul <- ggvegan:::arrowMul(arrows[, take], subset(fits1.2014.cca2, select = take, Score == 'sites'))
## ^ multiplier for arrows to scale them to the plot range
arrows[, take] <- arrows[, take] * mul  # scale biplot arrows w/ multiplier
arrows$Label<-gsub("Month(.*)", "\\1", arrows$Label) # change labels for months to be legible

summary(its1.2014.cca.2)
# Importance of components:
#                         CCA1    CA1
# Eigenvalue            0.4517 0.5226
# Proportion Explained  0.1528 0.1768
# Cumulative Proportion 0.1528 0.3296
screeplot(its1.2014.cca.2)

its1.2014.cca.all1<-ggplot() +
  geom_point(data = subset(fits1.2014.cca2, Score == 'sites'),aes(x = CCA1, y = CA1, colour=samps_its1.2014$Site), size=3) +
  geom_segment(data = arrows,aes(x = 0, y = 0, xend = CCA1, yend = CA1),
               arrow = arrow(length = unit(0.01, "npc"))) + theme_bw() + xlab("CCA1 (15.28%)") + ylab("CA1 (17.68%)") +
  geom_label(data = arrows,aes(label = Label, x = CCA1, y = CA1, fontface="bold",hjust = 0.5*(1 - sign(CCA1)), vjust = 0.5*(1-sign(CA1)))) +
  scale_x_continuous(expand = c(.1, .1)) +scale_y_continuous(expand = c(.1, .1)) + scale_colour_brewer("Site", palette = "Dark2")

ggsave(its1.2014.cca.all1,filename = "figures/ITS1_CCA_2014_5.17.21.pdf", width=10, height=8, dpi=600)

its1.2014.cca.all2<-ggplot() +
  geom_point(data = subset(fits1.2014.cca2, Score == 'sites'),aes(x = CCA1, y = CA1, colour=samps_its1.2014$Site), size=3) +
  geom_segment(data = arrows,aes(x = 0, y = 0, xend = CCA1, yend = CA1),
               arrow = arrow(length = unit(0.01, "npc"))) + theme_bw() + xlab("CCA1 (15.28%)") + ylab("CA1 (17.68%)") +
  geom_label(data = arrows,aes(label = Label, x = CCA1, y = CA1, fontface="bold",hjust = 0.5*(1 - sign(CCA1)), vjust = 0.5*(1-sign(CA1)))) +
  scale_x_continuous(expand = c(.1, .1)) +scale_y_continuous(expand = c(.1, .1)) + scale_colour_brewer("Site", palette = "Dark2")+
  geom_text(data=samps_its1.2014, aes(x = CCA1, y = CA1, label=Label), fontface = "plain", nudge_x = 0.1, nudge_y = 0.15)

ggsave(its1.2014.cca.all2,filename = "figures/ITS1_CCA_2014_Labeled_5.17.21.pdf", width=10, height=6, dpi=600)

### 2015 next

head(yr_2015)
## let's start with a CCA of several environmental variables
its1.2015.cca = cca(its1_2015 ~  Cu + Fe + Mg + Mn + Ni + P + S + Zn + Elevation + DustComplexity + Month, data = yr_2015 )

## what do we find?
its1.2015.cca
plot(its1.2015.cca)

## remember, Scaling
plot(its1.2015.cca, scaling = 1) ## this emphasizes relationships among sites
plot(its1.2015.cca, scaling = 2) ## this emphasizes relationships among species

## what do we find?
summary(its1.2015.cca)
its1.2015.cca
## 1) partition of scaled chi-squared Constrained = explained by variable,
##    Unconstrained = explained by residuals
## 2) Importance of componenents - which axes contribute
## 3) Scores: LC or WA? LC give best fit if you have no noise in constraining
##    variables, WA scores good with noisy data

## need to bootstrap to get the real r-square
## how much variation does it explain?
RsquareAdj(its1.2015.cca)
## ^^ use this b/c chance correlations can inflat R^2
## reminder: R^2 = % of variation in dependent variable explained by model

vif.cca(its1.2015.cca)

## we can then test for significance of the model by permutation
## if it is not signfiicant,
## it doesn't matter how much of the variation is explained
## Is it significant?
anova(its1.2015.cca, permutations = how(nperm=999))

## we can also do a permutation test by axis ********
anova(its1.2015.cca, by = "axis", permutations = how(nperm=999)) ### by CCA axis
anova(its1.2015.cca, by = "terms", permutations = how(nperm=999)) ### by variables
## this will help us interpret our CCA and we can see some variable are not significant
# shows that Cu and Elevation are significant; S is a litle over 0.05

## we can use model selection instead of picking variables we think are important

its.2015.cca2 = ordistep(cca(its1_2015 ~ 1, data = yr_2015 ),
                         scope=formula(its1.2015.cca),
                         direction = "forward",
                         permutations = how(nperm=999))

## this tells us that we really only need a few variables
## its1_2015 + Elevation + Month

## so lets rerun our model with these and see how it looks

its1.2015.cca.2 = cca(its1_2015 ~ Elevation+Month, data = yr_2015 )
summary(its1.2015.cca.2)
anova(its1.2015.cca.2, permutations = how(nperm=999))
anova(its1.2015.cca.2, by = "terms", permutations = how(nperm=999))

##What about r-squared

## original model
RsquareAdj(its1.2015.cca)
# 0.06365602

## new model
RsquareAdj(its1.2015.cca.2)
# 0.1906947

## now let rerun our model selection
its.2015.cca3 = ordistep(cca(its1_2015 ~ 1, data = yr_2015 ),
                         scope=formula(its1.2015.cca.2),
                         direction = "forward",
                         permutations = how(nperm=999))
## so our model that was selected is
# its1_2015 ~ Elevation + Month

## so how do we interpret this?
## the env explains 19.07% of the variation
## ^^ use Rsquareadj() b/c chance correlations can inflat R^2
## reminder: R^2 = % of variation in dependent variable explained by model

# Visualization w/ ggvegan
# autoplot(its1.cca)
# autoplot(its1.cca.2, axes=c(1,2),geom = c("point","text"), layers = c("sites", "biplot"),
#          legend.position = "right", title = "ITS1 CCA", ylab="CCA 2", xlab = "CCA 1")

fits1.2015.cca2 <- fortify(its1.2015.cca.2, axes = 1:2)  # fortify the ordination
## Note: object names are recycled here for everytime these figures are generated.
## Rerun all commands including & after fortify() command to ensure correct values are used in figures!
head(fits1.2015.cca2)
take <- c('CCA1', 'CCA2')  # which columns contain the scores we want
arrows <- subset(fits1.2015.cca2, Score == 'biplot')  # take only biplot arrow scores from cca df
samps_its1.2015 <- subset(fits1.2015.cca2, Score == 'sites')  # take only biplot arrow scores from cca df
samp_sites<-data.frame(Label=metadata$SampleID, Site=metadata$Site)
samps_its1.2015<-merge(samps_its1.2015, samp_sites, by="Label")
mul <- ggvegan:::arrowMul(arrows[, take], subset(fits1.2015.cca2, select = take, Score == 'sites'))
## ^ multiplier for arrows to scale them to the plot range
arrows[, take] <- arrows[, take] * mul  # scale biplot arrows w/ multiplier
arrows$Label<-gsub("Month(.*)", "\\1", arrows$Label) # change labels for months to be legible

summary(its1.2015.cca.2)
# Importance of components:
#                       CCA1   CCA2
# Eigenvalue            0.4842 0.1834
# Proportion Explained  0.1998 0.0757
# Cumulative Proportion 0.1998 0.2755
screeplot(its1.2015.cca.2)

its1.2015.cca.all1<-ggplot() +
  geom_point(data = subset(fits1.2015.cca2, Score == 'sites'),aes(x = CCA1, y = CCA2, colour=as.factor(samps_its1.2015$Site)), size=3) +
  geom_segment(data = arrows,aes(x = 0, y = 0, xend = CCA1, yend = CCA2),
               arrow = arrow(length = unit(0.01, "npc"))) + theme_bw() + xlab("CCA1 (19.98%)") + ylab("CCA2 (7.57%)") +
  geom_label(data = arrows,aes(label = Label, x = CCA1, y = CCA2, fontface="bold",hjust = 0.5*(1 - sign(CCA1)), vjust = 0.5*(1-sign(CCA2)))) +
  scale_x_continuous(expand = c(.1, .1)) +scale_y_continuous(expand = c(.1, .1)) + scale_colour_brewer("Site", palette = "Dark2")

ggsave(its1.2015.cca.all1,filename = "figures/ITS1_CCA_2015_5.17.21.pdf", width=10, height=8, dpi=600)

its1.2015.cca.all2<-ggplot() +
  geom_point(data = subset(fits1.2015.cca2, Score == 'sites'),aes(x = CCA1, y = CCA2, colour=as.factor(samps_its1.2015$Site)), size=3) +
  geom_segment(data = arrows,aes(x = 0, y = 0, xend = CCA1, yend = CCA2),
               arrow = arrow(length = unit(0.01, "npc"))) + theme_bw() + xlab("CCA1 (19.98%)") + ylab("CCA2 (7.57%)") +
  geom_label(data = arrows,aes(label = Label, x = CCA1, y = CCA2, fontface="bold",hjust = 0.5*(1 - sign(CCA1)), vjust = 0.5*(1-sign(CCA2)))) +
  scale_x_continuous(expand = c(.1, .1)) +scale_y_continuous(expand = c(.1, .1)) + scale_colour_brewer("Site", palette = "Dark2")+
  geom_text(data=samps_its1.2015, aes(x = CCA1, y = CCA2, label=Label), fontface = "plain", nudge_x = 0.1, nudge_y = 0.15)

ggsave(its1.2015.cca.all2,filename = "figures/ITS1_CCA_2015_Labeled_5.17.21.pdf", width=10, height=6, dpi=600)


#### Bacteria/Archaea CCA - by year ########

head(yr_2014)

###### 2014 first!

## let's start with a CCA of several environmental variables
bac.2014.cca = cca(bac_2014 ~  Cu + Fe + Mg + Mn + Ni + P + S + Zn + Elevation + DustComplexity + Month, data = yr_2014 )

## what do we find?
bac.2014.cca
plot(bac.2014.cca)

## remember, Scaling
plot(bac.2014.cca, scaling = 1) ## this emphasizes relationships among sites
plot(bac.2014.cca, scaling = 2) ## this emphasizes relationships among species

## what do we find?
summary(bac.2014.cca)
bac.2014.cca
## 1) partition of scaled chi-squared Constrained = explained by variable,
##    Unconstrained = explained by residuals
## 2) Importance of componenents - which axes contribute
## 3) Scores: LC or WA? LC give best fit if you have no noise in constraining
##    variables, WA scores good with noisy data

## need to bootstrap to get the real r-square
## how much variation does it explain?
RsquareAdj(bac.2014.cca)
## ^^ use this b/c chance correlations can inflat R^2
## reminder: R^2 = % of variation in dependent variable explained by model

vif.cca(bac.2014.cca)

## we can then test for significance of the model by permutation
## if it is not signfiicant,
## it doesn't matter how much of the variation is explained
## Is it significant?
anova(bac.2014.cca, permutations = how(nperm=999))

## we can also do a permutation test by axis ********
anova(bac.2014.cca, by = "axis", permutations = how(nperm=999)) ### by CCA axis
## or by terms
anova(bac.2014.cca, by = "terms", permutations = how(nperm=999)) ### by variables
## this will help us interpret our CCA and we can see some variable are not significant
# shows that Cu and Elevation are significant; S is a litle over 0.05

## we can use model selection instead of picking variables we think are important

bac.2014.cca2 = ordistep(cca(bac_2014 ~ 1, data = yr_2014 ),
                         scope=formula(bac.2014.cca),
                         direction = "forward",
                         permutations = how(nperm=999))

## this tells us that we really only need a few variables
## bac_2014 + S

## so lets rerun our model with these and see how it looks

bac.2014.cca.2 = cca(bac_2014 ~ S, data = yr_2014 )
summary(bac.2014.cca.2)
anova(bac.2014.cca.2, permutations = how(nperm=999))
anova(bac.2014.cca.2, by = "terms", permutations = how(nperm=999))

##What about r-squared

## original model
RsquareAdj(bac.2014.cca)

## new model
RsquareAdj(bac.2014.cca.2)

## so how do we interpret this?
## the env explains 3.32%% of the variation
## ^^ use Rsquareadj() b/c chance correlations can inflat R^2
## reminder: R^2 = % of variation in dependent variable explained by model

# Visualization w/ ggvegan
# autoplot(bac.cca)
# autoplot(bac.cca.2, axes=c(1,2),geom = c("point","text"), layers = c("sites", "biplot"),
#          legend.position = "right", title = "bac CCA", ylab="CCA 2", xlab = "CCA 1")

fbac.2014.cca2 <- fortify(bac.2014.cca.2, axes = 1:2)  # fortify the ordination
## Note: object names are recycled here for everytime these figures are generated.
## Rerun all commands including & after fortify() command to ensure correct values are used in figures!
head(fbac.2014.cca2)
take <- c('CCA1', 'CA1')  # which columns contain the scores we want
arrows <- subset(fbac.2014.cca2, Score == 'biplot')  # take only biplot arrow scores from cca df
samps_bac.2014 <- subset(fbac.2014.cca2, Score == 'sites')  # take only biplot arrow scores from cca df
samp_sites<-data.frame(Label=metadata$SampleID, Site=metadata$Site)
samps_bac.2014<-merge(samps_bac.2014, samp_sites, by="Label")
mul <- ggvegan:::arrowMul(arrows[, take], subset(fbac.2014.cca2, select = take, Score == 'sites'))
## ^ multiplier for arrows to scale them to the plot range
arrows[, take] <- arrows[, take] * mul  # scale biplot arrows w/ multiplier
arrows$Label<-gsub("Month(.*)", "\\1", arrows$Label) # change labels for months to be legible

summary(bac.2014.cca.2)
# Importance of components:
#                       CCA1    CA1
# Eigenvalue            0.5397 0.6454
# Proportion Explained  0.1299 0.1553
# Cumulative Proportion 0.1299 0.2852
screeplot(bac.2014.cca.2)

bac.2014.cca.all1<-ggplot() +
  geom_point(data = subset(fbac.2014.cca2, Score == 'sites'),aes(x = CCA1, y = CA1, colour=as.factor(samps_bac.2014$Site)), size=3) +
  geom_segment(data = arrows,aes(x = 0, y = 0, xend = CCA1, yend = CA1),
               arrow = arrow(length = unit(0.01, "npc"))) + theme_bw() + xlab("CCA1 (12.99%)") + ylab("CA1 (15.53%)") +
  geom_label(data = arrows,aes(label = Label, x = CCA1, y = CA1, fontface="bold",hjust = 0.5*(1 - sign(CCA1)), vjust = 0.5*(1-sign(CA1)))) +
  scale_x_continuous(expand = c(.1, .1)) +scale_y_continuous(expand = c(.1, .1)) + scale_colour_brewer("Site", palette = "Dark2")

ggsave(bac.2014.cca.all1,filename = "figures/16S_CCA_2014_5.17.21.pdf", width=10, height=8, dpi=600)

bac.2014.cca.all2<-ggplot() +
  geom_point(data = subset(fbac.2014.cca2, Score == 'sites'),aes(x = CCA1, y = CA1, colour=as.factor(samps_bac.2014$Site)), size=3) +
  geom_segment(data = arrows,aes(x = 0, y = 0, xend = CCA1, yend = CA1),
               arrow = arrow(length = unit(0.01, "npc"))) + theme_bw() + xlab("CCA1 (12.99%)") + ylab("CA1 (15.53%)") +
  geom_label(data = arrows,aes(label = Label, x = CCA1, y = CA1, fontface="bold",hjust = 0.5*(1 - sign(CCA1)), vjust = 0.5*(1-sign(CA1)))) +
  scale_x_continuous(expand = c(.1, .1)) +scale_y_continuous(expand = c(.1, .1)) + scale_colour_brewer("Site", palette = "Dark2")+
  geom_text(data=samps_bac.2014, aes(x = CCA1, y = CA1, label=Label), fontface = "plain", nudge_x = 0.1, nudge_y = 0.15)

ggsave(bac.2014.cca.all2,filename = "figures/16S_CCA_2014_Labeled_5.17.21.pdf", width=10, height=6, dpi=600)

### 2015 next

head(metadata)
## let's start with a CCA of several environmental variables
bac.2015.cca = cca(bac_2015 ~  Cu + Fe + Mg + Mn + Ni + P + S + Zn + Elevation + DustComplexity + Month, data = yr_2015 )

## what do we find?
bac.2015.cca
plot(bac.2015.cca)

## remember, Scaling
plot(bac.2015.cca, scaling = 1) ## this emphasizes relationships among sites
plot(bac.2015.cca, scaling = 2) ## this emphasizes relationships among species

## what do we find?
summary(bac.2015.cca)
bac.2015.cca
## 1) partition of scaled chi-squared Constrained = explained by variable,
##    Unconstrained = explained by residuals
## 2) Importance of componenents - which axes contribute
## 3) Scores: LC or WA? LC give best fit if you have no noise in constraining
##    variables, WA scores good with noisy data

## need to bootstrap to get the real r-square
## how much variation does it explain?
RsquareAdj(bac.2015.cca)
## ^^ use this b/c chance correlations can inflat R^2
## reminder: R^2 = % of variation in dependent variable explained by model

vif.cca(bac.2015.cca)

## we can then test for significance of the model by permutation
## if it is not signfiicant,
## it doesn't matter how much of the variation is explained
## Is it significant?
anova(bac.2015.cca, permutations = how(nperm=999))

## we can also do a permutation test by axis ********
anova(bac.2015.cca, by = "axis", permutations = how(nperm=999)) ### by CCA axis
anova(bac.2015.cca, by = "terms", permutations = how(nperm=999)) ### by variables
## this will help us interpret our CCA and we can see some variable are not significant
# shows that Cu and Elevation are significant; S is a litle over 0.05

## we can use model selection instead of picking variables we think are important

bac.2015.cca2 = ordistep(cca(bac_2015 ~ 1, data = yr_2015 ),
                         scope=formula(bac.2015.cca),
                         direction = "forward",
                         permutations = how(nperm=999))

## this tells us that we really only need a few variables
## bac_2015 + Elevation

## so lets rerun our model with these and see how it looks

bac.2015.cca.2 = cca(bac_2015 ~ Elevation+Month, data = yr_2015 )
summary(bac.2015.cca.2)
anova(bac.2015.cca.2, permutations = how(nperm=999))
anova(bac.2015.cca.2, by = "terms", permutations = how(nperm=999))

##What about r-squared

## original model
RsquareAdj(bac.2015.cca)
# 0.1180616

## new model
RsquareAdj(bac.2015.cca.2)
# 0.0886777

## now let rerun our model selection
bac.2015.cca3 = ordistep(cca(bac_2015 ~ 1, data = yr_2015 ),
                         scope=formula(bac.2015.cca.2),
                         direction = "forward",
                         permutations = how(nperm=999))
## so our model that was selected is
# bac_2015 ~ Elevation + Month

## so how do we interpret this?
## the env explains 8.87% of the variation
## ^^ use Rsquareadj() b/c chance correlations can inflat R^2
## reminder: R^2 = % of variation in dependent variable explained by model

# Visualization w/ ggvegan
# autoplot(bac.cca)
# autoplot(bac.cca.2, axes=c(1,2),geom = c("point","text"), layers = c("sites", "biplot"),
#          legend.position = "right", title = "bac CCA", ylab="CCA 2", xlab = "CCA 1")

fbac.2015.cca2 <- fortify(bac.2015.cca.2, axes = 1:2)  # fortify the ordination
## Note: object names are recycled here for everytime these figures are generated.
## Rerun all commands including & after fortify() command to ensure correct values are used in figures!
head(fbac.2015.cca2)
take <- c('CCA1', 'CCA2')  # which columns contain the scores we want
arrows <- subset(fbac.2015.cca2, Score == 'biplot')  # take only biplot arrow scores from cca df
samps_bac.2015 <- subset(fbac.2015.cca2, Score == 'sites')  # take only biplot arrow scores from cca df
samp_sites<-data.frame(Label=metadata$SampleID, Site=metadata$Site)
samps_bac.2015<-merge(samps_bac.2015, samp_sites, by="Label")
mul <- ggvegan:::arrowMul(arrows[, take], subset(fbac.2015.cca2, select = take, Score == 'sites'))
## ^ multiplier for arrows to scale them to the plot range
arrows[, take] <- arrows[, take] * mul  # scale biplot arrows w/ multiplier
arrows$Label<-gsub("Month(.*)", "\\1", arrows$Label) # change labels for months to be legible

summary(bac.2015.cca.2)
# Importance of components:
#                       CCA1    CCA2
# Eigenvalue            0.6322 0.40788
# Proportion Explained  0.1120 0.07223
# Cumulative Proportion 0.1120 0.18418
screeplot(bac.2015.cca.2)

bac.2015.cca.all1<-ggplot() +
  geom_point(data = subset(fbac.2015.cca2, Score == 'sites'),aes(x = CCA1, y = CCA2, colour=as.factor(samps_bac.2015$Site)), size=3) +
  geom_segment(data = arrows,aes(x = 0, y = 0, xend = CCA1, yend = CCA2),
               arrow = arrow(length = unit(0.01, "npc"))) + theme_bw() + xlab("CCA1 (11.20%)") + ylab("CCA2 (7.22%)") +
  geom_label(data = arrows,aes(label = Label, x = CCA1, y = CCA2, fontface="bold",hjust = 0.5*(1 - sign(CCA1)), vjust = 0.5*(1-sign(CCA2)))) +
  scale_x_continuous(expand = c(.1, .1)) +scale_y_continuous(expand = c(.1, .1)) + scale_colour_brewer("Site", palette = "Dark2")

ggsave(bac.2015.cca.all1,filename = "figures/16S_CCA_2015_5.17.21.pdf", width=10, height=8, dpi=600)

bac.2015.cca.all2<-ggplot() +
  geom_point(data = subset(fbac.2015.cca2, Score == 'sites'),aes(x = CCA1, y = CCA2, colour=as.factor(samps_bac.2015$Site)), size=3) +
  geom_segment(data = arrows,aes(x = 0, y = 0, xend = CCA1, yend = CCA2),
               arrow = arrow(length = unit(0.01, "npc"))) + theme_bw() + xlab("CCA1 (11.20%)") + ylab("CCA2 (7.22%)") +
  geom_label(data = arrows,aes(label = Label, x = CCA1, y = CCA2, fontface="bold",hjust = 0.5*(1 - sign(CCA1)), vjust = 0.5*(1-sign(CCA2)))) +
  scale_x_continuous(expand = c(.1, .1)) +scale_y_continuous(expand = c(.1, .1)) + scale_colour_brewer("Site", palette = "Dark2")+
  geom_text(data=samps_bac.2015, aes(x = CCA1, y = CCA2, label=Label), fontface = "plain", nudge_x = 0.1, nudge_y = 0.15)

ggsave(bac.2015.cca.all2,filename = "figures/16S_CCA_2015_Labeled_5.17.21.pdf", width=10, height=6, dpi=600)


#### Match relativized comp data to subsetted metadata (by site)####
# some sanity checks...
head(SJER) # sanity check
its1_table_RA[1:3,1:3] # sanity check to see if relativized OTU table has rownames for next function
bac_table_RA[1:3,1:3]

# matching data with user defined function
#
match_dat<-function(compdata, subset_metadata){
  subset_comp_data = pullrow<-(is.element(row.names(compdata), row.names(subset_metadata)))
  ### * comp data and metadata need to have row names - rownames should be Sample IDs
  subset_comp_data=compdata[pullrow,]
  return(subset_comp_data)
}

# SJER
sjer_its1<-match_dat(its1_table_RA, SJER)
sjer_its1[1:3,1:3]
sjer_bac<-match_dat(bac_table_RA, SJER)
sjer_bac[1:3,1:3]

# PROVIDENCE
prov_its1<-match_dat(its1_table_RA, PROVIDENCE)
prov_its1[1:3,1:3]
prov_bac<-match_dat(bac_table_RA, PROVIDENCE)
prov_bac[1:3,1:3]

# SHORTHAIR
sh_its1<-match_dat(its1_table_RA, SHORTHAIR)
sh_its1[1:3,1:3]
sh_bac<-match_dat(bac_table_RA, SHORTHAIR)
sh_bac[1:3,1:3]

# SOAPROOT
sr_its1<-match_dat(its1_table_RA, SOAPROOT)
sr_its1[1:3,1:3]
sr_bac<-match_dat(bac_table_RA, SOAPROOT)
sr_bac[1:3,1:3]


########  CCA/RDAs by SITE  ###############
site_list

## CCA does not try to display all variation in the data,
## but only the part that can be explained by the used constraints. *******
## Consequently, the results are strongly dependent on the set of constraints
## The shotgun method is to use all environmental variables as constraints,
## however, rely on non-metric multidimensional scaling (metaMDS) and environmental
## interpretation after analysis (envfit, ordisurf).
## CCA is a good choice if the user has clear and strong a priori hypotheses on constraints

### SJER CCA for ITS1
head(sjer_its1)
head(metadata)
## let's start with a CCA of several environmental variables
sjer_its.cca = cca(sjer_its1 ~  Cu + Fe + Mg + Mn + Ni + P + S + Zn + Elevation, data = SJER )

## what do we find?
sjer_its.cca
plot(sjer_its.cca)

## remember, Scaling
plot(sjer_its.cca, scaling = 1)
## this emphasizes relationships among sites
plot(sjer_its.cca, scaling = 2)
## this emphasizes relationships among species

## what do we find?
summary(sjer_its.cca)
sjer_its.cca
## 1) partition of scaled chi-squared Constrained = explained by variable,
##    Unconstrained = explained by residuals
## 2) Importance of componenents - which axes contribute
## 3) Scores: LC or WA? LC give best fit if you have no noise in constraining
##    variables, WA scores good with noisy data

## need to bootstrap to get the real r-square
## how much variation does it explain?
RsquareAdj(sjer_its.cca)
## ^^ use this b/c chance correlations can inflat R^2
## reminder: R^2 = % of variation in dependent variable explained by model

vif.cca(sjer_its.cca)

## we can then test for significance of the model by permutation
## if it is not signfiicant,
## it doesn't matter how much of the variation is explained
## Is it significant?
anova(sjer_its.cca, permutations = how(nperm=999))

## we can also do a permutation test by axis ********
anova(sjer_its.cca, by = "axis", permutations = how(nperm=999)) ### by CCA axis
## or by terms
anova(sjer_its.cca, by = "terms", permutations = how(nperm=999)) ### by variables
## this will help us interpret our CCA and we can see some variable are not significant

## we can use model selection instead of picking variables we think are important

sjer_its.cca2 = ordistep(cca(sjer_its1 ~ 1, data = SJER ),
                  scope=formula(sjer_its.cca),
                  direction = "forward",
                  permutations = how(nperm=999))

## this tells us that we really only need a few variables
## sjer_its1 ~ Fe

## so lets rerun our model with these and see how it looks

sjer_its.cca.2 = cca(sjer_its1 ~  Fe, data = SJER )
summary(sjer_its.cca.2)
anova(sjer_its.cca.2, permutations = how(nperm=999))
anova(sjer_its.cca.2, by = "terms", permutations = how(nperm=999))

##What about r-squared

## original model
RsquareAdj(sjer_its.cca)

## new model
RsquareAdj(sjer_its.cca.2)

## but we have limited our model selection using "scope" to just what was in our initial model

## the "." here mean everything in the dataset
KLcca3 = cca(sjer_its1 ~ ., data = SJER )

## now let rerun our model selection
sjer_its.cca.3 = ordistep(cca(sjer_its1 ~ 1, data = SJER ),
                  scope=formula(sjer_its.cca.2),
                  direction = "forward",
                  permutations = how(nperm=999))

## so our model that was selected is
# sjer_its1 ~ Fe
## so how do we interpret this?
## the env explains 9.8% of the variation
## ^^ use Rsquareadj() b/c chance correlations can inflat R^2
## reminder: R^2 = % of variation in dependent variable explained by model


### SJER CCA for 16S
head(sjer_bac)
head(metadata)
## let's start with a CCA of several environmental variables
sjer_bac.cca = cca(sjer_bac ~  Cu + Fe + Mg + Mn + Ni + P + S + Zn + Elevation, data = SJER )

## what do we find?
sjer_bac.cca
plot(sjer_bac.cca)

## remember, Scaling
plot(sjer_bac.cca, scaling = 1)
## this emphasizes relationships among sites
plot(sjer_bac.cca, scaling = 2)
## this emphasizes relationships among species

## what do we find?
summary(sjer_bac.cca)
sjer_bac.cca
## 1) partition of scaled chi-squared Constrained = explained by variable,
##    Unconstrained = explained by residuals
## 2) Importance of componenents - which axes contribute
## 3) Scores: LC or WA? LC give best fit if you have no noise in constraining
##    variables, WA scores good with noisy data

## need to bootstrap to get the real r-square
## how much variation does it explain?
RsquareAdj(sjer_bac.cca)
## ^^ use this b/c chance correlations can inflat R^2
## reminder: R^2 = % of variation in dependent variable explained by model

vif.cca(sjer_bac.cca)

## we can then test for significance of the model by permutation
## if it is not signfiicant,
## it doesn't matter how much of the variation is explained
## Is it significant?
anova(sjer_bac.cca, permutations = how(nperm=999))

## we can also do a permutation test by axis ********
anova(sjer_bac.cca, by = "axis", permutations = how(nperm=999)) ### by CCA axis
## or by terms
anova(sjer_bac.cca, by = "terms", permutations = how(nperm=999)) ### by variables
## this will help us interpret our CCA and we can see some variable are not significant

## we can use model selection instead of picking variables we think are important

sjer_bac.cca2 = ordistep(cca(sjer_bac ~ 1, data = SJER ),
                         scope=formula(sjer_bac.cca),
                         direction = "forward",
                         permutations = how(nperm=999))

## this tells us that we really only need a few variables
## sjer_bac ~ Fe

## so lets rerun our model with these and see how it looks

sjer_bac.cca.2 = cca(sjer_bac ~  Fe, data = SJER )
summary(sjer_bac.cca.2)
anova(sjer_bac.cca.2, permutations = how(nperm=999))
anova(sjer_bac.cca.2, by = "terms", permutations = how(nperm=999))

##What about r-squared

## original model
RsquareAdj(sjer_bac.cca)

## new model
RsquareAdj(sjer_bac.cca.2)

## so our model that was selected is
# sjer_bac ~ Fe
## so how do we interpret this?
## the env explains 5.7% of the variation
## ^^ use Rsquareadj() b/c chance correlations can inflat R^2
## reminder: R^2 = % of variation in dependent variable explained by model

##### Providence

### SJER CCA for ITS1
head(sjer_its1)
head(metadata)
## let's start with a CCA of several environmental variables
sjer_its.cca = cca(sjer_its1 ~  Cu + Fe + Mg + Mn + Ni + P + S + Zn + Elevation, data = SJER )

## what do we find?
sjer_its.cca
plot(sjer_its.cca)

## remember, Scaling
plot(sjer_its.cca, scaling = 1)
## this emphasizes relationships among sites
plot(sjer_its.cca, scaling = 2)
## this emphasizes relationships among species

## what do we find?
summary(sjer_its.cca)
sjer_its.cca
## 1) partition of scaled chi-squared Constrained = explained by variable,
##    Unconstrained = explained by residuals
## 2) Importance of componenents - which axes contribute
## 3) Scores: LC or WA? LC give best fit if you have no noise in constraining
##    variables, WA scores good with noisy data

## need to bootstrap to get the real r-square
## how much variation does it explain?
RsquareAdj(sjer_its.cca)
## ^^ use this b/c chance correlations can inflat R^2
## reminder: R^2 = % of variation in dependent variable explained by model

vif.cca(sjer_its.cca)

## we can then test for significance of the model by permutation
## if it is not signfiicant,
## it doesn't matter how much of the variation is explained
## Is it significant?
anova(sjer_its.cca, permutations = how(nperm=999))

## we can also do a permutation test by axis ********
anova(sjer_its.cca, by = "axis", permutations = how(nperm=999)) ### by CCA axis
## or by terms
anova(sjer_its.cca, by = "terms", permutations = how(nperm=999)) ### by variables
## this will help us interpret our CCA and we can see some variable are not significant

## we can use model selection instead of picking variables we think are important

sjer_its.cca2 = ordistep(cca(sjer_its1 ~ 1, data = SJER ),
                         scope=formula(sjer_its.cca),
                         direction = "forward",
                         permutations = how(nperm=999))

## this tells us that we really only need a few variables
## sjer_its1 ~ Fe

## so lets rerun our model with these and see how it looks

sjer_its.cca.2 = cca(sjer_its1 ~  Fe, data = SJER )
summary(sjer_its.cca.2)
anova(sjer_its.cca.2, permutations = how(nperm=999))
anova(sjer_its.cca.2, by = "terms", permutations = how(nperm=999))

##What about r-squared

## original model
RsquareAdj(sjer_its.cca)

## new model
RsquareAdj(sjer_its.cca.2)

## but we have limited our model selection using "scope" to just what was in our initial model

## the "." here mean everything in the dataset
KLcca3 = cca(sjer_its1 ~ ., data = SJER )

## now let rerun our model selection
sjer_its.cca.3 = ordistep(cca(sjer_its1 ~ 1, data = SJER ),
                          scope=formula(sjer_its.cca.2),
                          direction = "forward",
                          permutations = how(nperm=999))

## so our model that was selected is
# sjer_its1 ~ Fe
## so how do we interpret this?
## the env explains 9.8% of the variation
## ^^ use Rsquareadj() b/c chance correlations can inflat R^2
## reminder: R^2 = % of variation in dependent variable explained by model


### SJER CCA for 16S
head(sjer_bac)
head(metadata)
## let's start with a CCA of several environmental variables
sjer_bac.cca = cca(sjer_bac ~  Cu + Fe + Mg + Mn + Ni + P + S + Zn + Elevation, data = SJER )

## what do we find?
sjer_bac.cca
plot(sjer_bac.cca)

## remember, Scaling
plot(sjer_bac.cca, scaling = 1)
## this emphasizes relationships among sites
plot(sjer_bac.cca, scaling = 2)
## this emphasizes relationships among species

## what do we find?
summary(sjer_bac.cca)
sjer_bac.cca
## 1) partition of scaled chi-squared Constrained = explained by variable,
##    Unconstrained = explained by residuals
## 2) Importance of componenents - which axes contribute
## 3) Scores: LC or WA? LC give best fit if you have no noise in constraining
##    variables, WA scores good with noisy data

## need to bootstrap to get the real r-square
## how much variation does it explain?
RsquareAdj(sjer_bac.cca)
## ^^ use this b/c chance correlations can inflat R^2
## reminder: R^2 = % of variation in dependent variable explained by model

vif.cca(sjer_bac.cca)

## we can then test for significance of the model by permutation
## if it is not signfiicant,
## it doesn't matter how much of the variation is explained
## Is it significant?
anova(sjer_bac.cca, permutations = how(nperm=999))

## we can also do a permutation test by axis ********
anova(sjer_bac.cca, by = "axis", permutations = how(nperm=999)) ### by CCA axis
## or by terms
anova(sjer_bac.cca, by = "terms", permutations = how(nperm=999)) ### by variables
## this will help us interpret our CCA and we can see some variable are not significant

## we can use model selection instead of picking variables we think are important

sjer_bac.cca2 = ordistep(cca(sjer_bac ~ 1, data = SJER ),
                         scope=formula(sjer_bac.cca),
                         direction = "forward",
                         permutations = how(nperm=999))

## this tells us that we really only need a few variables
## sjer_bac ~ Fe

## so lets rerun our model with these and see how it looks

sjer_bac.cca.2 = cca(sjer_bac ~  Fe, data = SJER )
summary(sjer_bac.cca.2)
anova(sjer_bac.cca.2, permutations = how(nperm=999))
anova(sjer_bac.cca.2, by = "terms", permutations = how(nperm=999))

##What about r-squared

## original model
RsquareAdj(sjer_bac.cca)

## new model
RsquareAdj(sjer_bac.cca.2)

## so our model that was selected is
# sjer_bac ~ Fe
## so how do we interpret this?
## the env explains 5.7% of the variation
## ^^ use Rsquareadj() b/c chance correlations can inflat R^2
## reminder: R^2 = % of variation in dependent variable explained by model











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

# custom color palette
head(metadata)

warm2cold1<-get_palette(paste0("#", c("720026", "D14D60", "0077B6","03045E")),k=4)
names(warm2cold1) <- levels(metadata$Elevation)

warm2cold2<-get_palette(paste0("#", c("720026", "D14D60", "0077B6","03045E")),k=4)
names(warm2cold2) <- levels(metadata$Site)

