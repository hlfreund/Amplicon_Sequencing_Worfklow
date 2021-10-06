library("phyloseq")
library("vegan")
library("DESeq2")
library("ggplot2")
library("dendextend")
library("tidyr")
library("viridis")
library("reshape")
library(reshape2)
library(ggpubr)
library(dplyr)
library(ape)
library(shades)

#### Import metadata ####
metadata<-as.data.frame(read.csv("data/Mouse_Seqs/mapzym.csv"), header=TRUE)
head(metadata)
metadata<-na.omit(metadata)

metadata$Description<-gsub("\xca", "", metadata$Description) ## gsub is global sub (does not just remove first instance of pattern, but multiple)
rownames(metadata)<-metadata$SampleID
head(metadata)

metadata$lungTissue<-sub("yes", "lung", metadata$lungTissue)
metadata$lungTissue<-sub("no", "BALF", metadata$lungTissue)
head(metadata)

# create colors for later figures
colorset = melt(c(Alt="#540b0e",Con="#4c956c",Silica="#98c1d9"))
colorset$Exposed<-rownames(colorset)
colnames(colorset)[which(names(colorset) == "value")] <- "color"
colorset

metadata<-merge(metadata, colorset, by="Exposed")
head(metadata)
metadata$color <- as.character(metadata$color)
rownames(metadata)<-metadata$SampleID

#### Import and reformat the data ####
# bacteria first
bac.ASV_counts<-data.frame(readRDS("data/Mouse_Seqs/Results_16S/16S_ASVs_Counts_dada2_6.20.2021_Robject.rds", refhook = NULL))
dim(bac.ASV_counts)
head(bac.ASV_counts)
rownames(bac.ASV_counts)<-bac.ASV_counts$ASV_ID

# fungi next
its2.ASV_counts<-data.frame(readRDS("data/Mouse_Seqs/Results_ITS2/ITS2_ASVs_Counts_dada2_6.20.2021_Robject.rds", refhook = NULL))
dim(its2.ASV_counts)
head(its2.ASV_counts)
rownames(its2.ASV_counts)<-its2.ASV_counts$ASV_ID

### ^^^ *** Note to self - before you run stats, make sure your table is in a SAMPLE x ASV/SPECIES orientation!!

#### Import ASV to taxa sheet for ASV identification ####

# bacteria first
bac.ASV_tax<-data.frame(readRDS("data/Mouse_Seqs/Results_16S/ASVs_Taxonomy_dada2_6.16.2021_Robject.rds", refhook = NULL))
head(bac.ASV_tax)
bac.ASV_tax$Genus<-gsub("\\[(.*)\\]", "\\1", bac.ASV_tax$Genus) ## drop brackets around Eubacterium while not losing the string inside brackets

bac.ASV_tax[is.na(bac.ASV_tax)]<- "Unknown"
bac.ASV_tax$Species<-gsub("Unknown", "unknown", bac.ASV_tax$Species) ## drop brackets around Eubacterium while not losing the string inside brackets

head(bac.ASV_tax)
class(bac.ASV_tax)
bac.ASV_tax$ASV_ID<-rownames(bac.ASV_tax)
head(bac.ASV_tax)

# fungi next
its2.ASV_tax<-data.frame(readRDS("data/Mouse_Seqs/Results_ITS2/ITS2_ASVs_Taxonomy_dada2_6.16.2021_Robject.rds", refhook = NULL))
head(its2.ASV_tax)
its2.ASV_tax$ASV_ID<-rownames(its2.ASV_tax)
its2.ASV_tax<-as.data.frame(lapply(its2.ASV_tax, function(x) gsub("[a-z]__", "", x)))

its2.ASV_tax[is.na(its2.ASV_tax)]<- "Unknown"
its2.ASV_tax$Species<-gsub("Unknown", "unknown", its2.ASV_tax$Species) ## drop brackets around Eubacterium while not losing the string inside brackets

head(its2.ASV_tax)
rownames(its2.ASV_tax)<-its2.ASV_tax$ASV_ID
head(its2.ASV_tax)

#### Data Formatting and Transformation ####

# bacteria first
bac.ASV_dat<-merge(bac.ASV_counts,bac.ASV_tax, by="ASV_ID")
head(bac.ASV_dat)

bac.ASV_dat<-subset(bac.ASV_dat, Kingdom!="Unknown") ## keep only bacteria and archaean -- drop Unknowns
bac.ASV_dat<-subset(bac.ASV_dat, Phylum!="Unknown") ## keep only bacteria and archaean -- drop Unknowns=
head(bac.ASV_dat)
bac.ASV_dat<-subset(bac.ASV_dat, Class!="Chloroplast") ## keep only bacteria -- exclude Chloroplast sequences
bac.ASV_dat<-subset(bac.ASV_dat, Order!="Chloroplast") ## keep only bacteria -- exclude Chloroplast sequences
bac.ASV_dat<-subset(bac.ASV_dat, Family!="Mitochondria") ## keep only bacteria -- exclude Chloroplast sequences

'Chloroplast' %in% bac.ASV_dat # check if Chloroplast counts are still in df, should be false because they've been removed
'Mitochondria' %in% bac.ASV_dat # check if Chloroplast counts are still in df, should be false because they've been removed

head(bac.ASV_dat)
rownames(bac.ASV_dat)<-bac.ASV_dat$ASV_ID
head(bac.ASV_dat)

"Undetermined" %in% bac.ASV_dat
b.dat.m<-melt(bac.ASV_dat)
head(b.dat.m)
colnames(b.dat.m)[which(names(b.dat.m) == "variable")] <- "SampleID"
colnames(b.dat.m)[which(names(b.dat.m) == "value")] <- "Count"

# ASVs (rows) x SampleID (cols)
bac.counts<-as.data.frame(dcast(b.dat.m, ASV_ID~SampleID, value.var="Count", fun.aggregate=sum)) ###
rownames(bac.counts)<-bac.counts$ASV_ID
bac.counts<-subset(bac.counts, select=-c(ASV_ID))
head(bac.counts)

# Sample ID x ASVs
bac.table<-as.data.frame(dcast(b.dat.m, SampleID~ASV_ID, value.var="Count", fun.aggregate=sum)) ###
rownames(bac.table)<-bac.table$SampleID
bac.table<-subset(bac.table, select=-c(SampleID))
head(bac.table)

# fungi next
its2.ASV_dat<-merge(its2.ASV_counts,its2.ASV_tax, by="ASV_ID")
head(its2.ASV_dat)

its2.ASV_dat<-subset(its2.ASV_dat, Kingdom!="Unknown") ## keep only its2teria and archaean -- drop Unknowns
#its2.ASV_dat<-subset(its2.ASV_dat, Phylum!="Unknown") ## keep only its2teria and archaean -- drop Unknowns=
head(its2.ASV_dat)
its2.ASV_dat<-subset(its2.ASV_dat, Class!="Chloroplast") ## keep only its2teria -- exclude Chloroplast sequences
its2.ASV_dat<-subset(its2.ASV_dat, Order!="Chloroplast") ## keep only its2teria -- exclude Chloroplast sequences
its2.ASV_dat<-subset(its2.ASV_dat, Family!="Mitochondria") ## keep only its2teria -- exclude Chloroplast sequences

'Chloroplast' %in% its2.ASV_dat # check if Chloroplast counts are still in df, should be false because they've been removed
'Mitochondria' %in% its2.ASV_dat # check if Chloroplast counts are still in df, should be false because they've been removed

head(its2.ASV_dat)
rownames(its2.ASV_dat)<-its2.ASV_dat$ASV_ID
head(its2.ASV_dat)

"Undetermined" %in% its2.ASV_dat
f.dat.m<-melt(its2.ASV_dat)
head(f.dat.m)
colnames(f.dat.m)[which(names(f.dat.m) == "variable")] <- "SampleID"
colnames(f.dat.m)[which(names(f.dat.m) == "value")] <- "Count"

# ASVs (rows) x SampleID (cols)
its2.counts<-as.data.frame(dcast(f.dat.m, ASV_ID~SampleID, value.var="Count", fun.aggregate=sum)) ###
rownames(its2.counts)<-its2.counts$ASV_ID
its2.counts<-subset(its2.counts, select=-c(ASV_ID))
head(its2.counts)

# Sample ID x ASVs
its2.table<-as.data.frame(dcast(f.dat.m, SampleID~ASV_ID, value.var="Count", fun.aggregate=sum)) ###
rownames(its2.table)<-its2.table$SampleID
its2.table<-subset(its2.table, select=-c(SampleID))
head(its2.table)

## Reorder metadata to have same rows as ASV tables
bac.counts=bac.counts[,rownames(metadata)]
its2.counts=its2.counts[,rownames(metadata)]
# ** ^ this indexing method will only work if the two dfs have the same # of rows AND the same row names!

#### 16S: Transform Data (Normalize for Sampling Depth) ####
# The point of these two transformations, the VST and the rlog, is to remove the dependence of the variance on the mean, particularly the high variance of the logarithm of count data when the mean is low
# More info here: http://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#count-data-transformations

head(bac.counts)

# first we need to make a DESeq2 object
bac.deseq_counts <- DESeqDataSetFromMatrix(bac.counts, colData = metadata, design = ~Exposed)
# Rows of colData (metadata) correspond to columns of countData (where SampleIDs are col names; not ASV table in sample x species format)
# we have to include the "colData" and "design" arguments because they are
# required, as they are needed for further downstream processing by DESeq2,
# but for our purposes of simply transforming the data right now, they don't
# matter
bac.deseq_counts <- estimateSizeFactors(bac.deseq_counts, type = "poscounts")
plotSparsity(bac.deseq_counts)

# now followed by the transformation function (Variance Stabilizing Transformation):
b.deseq_counts_vst <- varianceStabilizingTransformation(bac.deseq_counts)
# b.deseq_counts_vst <- vst(bac.deseq_counts) # better to use vst() when you have over 1000 ASVs
# This function calculates a variance stabilizing transformation (VST) from the fitted dispersionmean relation(s)
# then transforms the count data (normalized by division by the size factors or normalization factors), yielding a matrix of values which are now approximately homoskedastic
# homoskedastic = having constant variance along the range of mean values.
# The transformation also normalizes with respect to library size

# NOTE: If you get this error here with your dataset: "Error in
# estimateSizeFactorsForMatrix(counts(object), locfunc =locfunc, : every
# gene contains at least one zero, cannot compute log geometric means", that
# can be because the count table is sparse with many zeroes, which is common
# with marker-gene surveys. In that case you'd need to use a specific
# function first that is equipped to deal with that.
# You could run:
# deseq_counts <- estimateSizeFactors(deseq_counts, type = "poscounts")
# now followed by the transformation function:
# deseq_counts_vst <- varianceStabilizingTransformation(deseq_counts)

# and here is pulling out our transformed table
b.vst <- assay(b.deseq_counts_vst)
head(b.vst)

# and calculating our Euclidean distance matrix
b.euc_dist <- dist(t(b.vst), method = "euclidean")

## Trying rlog transformation
b.rl_count_tab<-rlog(bac.deseq_counts, blind=FALSE)
# In the above function calls, we specified blind = FALSE, which means that differences between cell lines and treatment (the variables in the design) will not contribute to the expected variance-mean trend of the experiment.
# The experimental design is not used directly in the transformation, only in estimating the global amount of variability in the counts. For a fully unsupervised transformation, one can set blind = TRUE (which is the default).
#bac.deseq_counts <- estimateSizeFactors(bac.deseq_counts, type = "poscounts")

## CLR Transformation
bac.clr<-microbiome::transform(bac.ASV_table, 'clr')
head(bac.clr)

b.clr.dist<-dist(bac.clr, method = "euclidean")
b.clr.clust<-hclust(b.clr.dist)

b.clr_dend <- as.dendrogram(b.clr.clust, hang=0.1)
b.clr.dend_cols <- as.character(metadata$exp_color[order.dendrogram(b.clr_dend)])
labels_colors(b.clr_dend) <- b.clr.dend_cols

plot(b.clr_dend, ylab="CLR Euclidean Distance")
## log2 Transformation
bac.log2<-as.data.frame(lapply(bac.table, function(x) log2(x+1)))
bac.log2<-bac.log2[order(match(bac.log2[,1],bac.table[,1])),]
# ordering meta_quant based on position order of matches from meta_quant[,1] found in metadata$Cu
# cannot call column using $ for array, only df
rownames(bac.log2)<-rownames(bac.table) # can only do this if you're sure the arrays/dfs are in the same exact order
rownames(bac.log2)
head(bac.log2)

bac.l2.dist<-dist(bac.log2, method = "euclidean")
plot(hclust(bac.l2.dist))

#### ITS2: Transform Data (Normalize for Sampling Depth) ####
# The point of these two transformations, the VST and the rlog, is to remove the dependence of the variance on the mean, particularly the high variance of the logarithm of count data when the mean is low
# More info here: http://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#count-data-transformations

#its2.counts<-na.omit(its2.counts)

# first we need to make a DESeq2 object
its2.deseq_counts <- DESeqDataSetFromMatrix(its2.counts, colData = metadata, design = ~Exposed)
# Rows of colData (metadata) correspond to columns of countData (where SampleIDs are col names; not ASV table in sample x species format)
# we have to include the "colData" and "design" arguments because they are
# required, as they are needed for further downstream processing by DESeq2,
# but for our purposes of simply transforming the data right now, they don't
# matter
its2.deseq_counts <- estimateSizeFactors(its2.deseq_counts, type = "poscounts")
# now followed by the transformation function (Variance Stabilizing Transformation):
f.deseq_counts_vst <- varianceStabilizingTransformation(its2.deseq_counts)
# b.deseq_counts_vst <- vst(bac.deseq_counts) # better to use vst() when you have over 1000 ASVs
# This function calculates a variance stabilizing transformation (VST) from the fitted dispersionmean relation(s)
# then transforms the count data (normalized by division by the size factors or normalization factors), yielding a matrix of values which are now approximately homoskedastic
# homoskedastic = having constant variance along the range of mean values.
# The transformation also normalizes with respect to library size

# NOTE: If you get this error here with your dataset: "Error in
# estimateSizeFactorsForMatrix(counts(object), locfunc =locfunc, : every
# gene contains at least one zero, cannot compute log geometric means", that
# can be because the count table is sparse with many zeroes, which is common
# with marker-gene surveys. In that case you'd need to use a specific
# function first that is equipped to deal with that.
# You could run:
# deseq_counts <- estimateSizeFactors(deseq_counts, type = "poscounts")
# now followed by the transformation function:
# deseq_counts_vst <- varianceStabilizingTransformation(deseq_counts)

# and here is pulling out our transformed table
f.vst <- assay(f.deseq_counts_vst)
head(f.vst)

# and calculating our Euclidean distance matrix
f.euc_dist <- dist(t(f.vst), method = "euclidean")

## Trying rlog transformation
f.rl_count_tab<-rlog(its2.deseq_counts, blind=FALSE)
# and calculating our Euclidean distance matrix
# f.rl_dist <- dist(t(assay(f.rl_count_tab)))
# plot(hclust(f.rl_dist))
# In the above function calls, we specified blind = FALSE, which means that differences between cell lines and treatment (the variables in the design) will not contribute to the expected variance-mean trend of the experiment.
# The experimental design is not used directly in the transformation, only in estimating the global amount of variability in the counts. For a fully unsupervised transformation, one can set blind = TRUE (which is the default).
plotSparsity(its2.deseq_counts)

## CLR Transformation
its2.clr<-microbiome::transform(its2.table, 'clr')
head(its2.clr)

## log2 Transformation
its2.log2<-as.data.frame(lapply(its2.table, function(x) log2(x+1)))
its2.log2<-its2.log2[order(match(its2.log2[,1],its2.table[,1])),]
# ordering meta_quant based on position order of matches from meta_quant[,1] found in metadata$Cu
# cannot call column using $ for array, only df
rownames(its2.log2)<-rownames(its2.table) # can only do this if you're sure the arrays/dfs are in the same exact order
rownames(its2.log2)
head(its2.log2)
its2.l2.dist<-dist(its2.log2, method = "euclidean")
plot(hclust(its2.l2.dist))

#### 16S - Hierarchical Clustering ####

# bacteria
b.euc_dist
b.euc_clust <- hclust(b.euc_dist, method="ward.D2")

# hclust objects like this can be plotted with the generic plot() function
plot(b.euc_clust)
# but i like to change them to dendrograms for two reasons:
# 1) it's easier to color the dendrogram plot by groups
# 2) if wanted you can rotate clusters with the rotate()
#    function of the dendextend package

b.euc_dend <- as.dendrogram(b.euc_clust, hang=0.1)
b.dend_cols <- as.character(metadata$color[order.dendrogram(b.euc_dend)])
labels_colors(b.euc_dend) <- b.dend_cols

pdf(file="figures/16S_VST_cluster_by.ExposMat_6.21.21.pdf", width=8, height=8)
#tiff(file="figures/16s_cluster_by.ExposMat_6.20.21.tiff",width=10, height=10, units="in", res=600)
plot(b.euc_dend, ylab="VST Euc. dist.") # colorized by exposure material
title(main = "Hierarchical Clustering: Bacteria/Archaea", sub = "Colorized by Exposure Material",
      cex.main = 2,   font.main= 2,
      cex.sub = 0.75, font.sub = 3
)
dev.off()

b.dend_cols2 <- as.character(metadata$color2[order.dendrogram(b.euc_dend)])
labels_colors(b.euc_dend) <- b.dend_cols2
# color reminder --> lung="#a53860" aka purple-red, BALF="#788bff" aka periwinkle

pdf(file="figures/16S_VST_cluster_by.LungTissue_6.21.21.pdf", width=8, height=8)
#tiff(file="figures/16s_cluster_by.LungTissue_6.20.21.tiff",width=10, height=10, units="in", res=600)
plot(b.euc_dend, ylab="VST Euc. dist.") # colorized by exposure material
title(main = "Hierarchical Clustering: Bacteria/Archaea", sub = "Colorized by Lung Tissue Type",
      cex.main = 2,   font.main= 2,
      cex.sub = 0.75, font.sub = 3
)
dev.off()

#### ITS2 - Hierarchical Clustering ####

f.euc_clust <- hclust(f.euc_dist, method="ward.D2")

# hclust objects like this can be plotted with the generic plot() function
plot(f.euc_clust)
# but i like to change them to dendrograms for two reasons:
# 1) it's easier to color the dendrogram plot by groups
# 2) if wanted you can rotate clusters with the rotate()
#    function of the dendextend package

f.euc_dend <- as.dendrogram(f.euc_clust, hang=0.1)
f.dend_cols <- as.character(metadata$color[order.dendrogram(f.euc_dend)])
labels_colors(f.euc_dend) <- f.dend_cols

pdf(file="figures/its2_VST_cluster_by.ExposMat_6.21.21.pdf", width=8, height=8)
#tiff(file="figures/its2_cluster_by.ExposMat_6.20.21.tiff",width=10, height=6, units="in", res=600)
plot(f.euc_dend, ylab="VST Euc. dist.") # colorized by exposure material
title(main = "Hierarchical Clustering: Fungi", sub = "Colorized by Exposure Material",
      cex.main = 2,   font.main= 2,
      cex.sub = 0.75, font.sub = 3
)
dev.off()

f.dend_cols2 <- as.character(metadata$color2[order.dendrogram(f.euc_dend)])
labels_colors(f.euc_dend) <- f.dend_cols2
# color reminder --> lung="#a53860" aka purple-red, BALF="#788bff" aka periwinkle

pdf(file="figures/ITS2_VST_cluster_by.LungTissue_6.21.21.pdf", width=8, height=8)
#tiff(file="figures/ITS2_cluster_by.LungTissue_6.20.21.tiff",width=10, height=6, units="in", res=400)
plot(f.euc_dend, ylab="VST Euc. dist.") # colorized bt exposure material
title(main = "Hierarchical Clustering: Fungi", sub = "Colorized fy Lung Tissue Type",
      cex.main = 2,   font.main= 2,
      cex.sub = 0.75, font.sub = 3
)
dev.off()


#### 16S PCA ####
head(bac.clr)
b.pca<-rda(b.clr.dist, scale=FALSE) # PCA using transposed VST transformed data
barplot(as.vector(b.pca$CA$eig)/sum(b.pca$CA$eig))

# Calculate the percent of variance explained by first two axes
sum((as.vector(b.pca$CA$eig)/sum(b.pca$CA$eig))[1:2]) # 50.79%
## cumulative proportion is calculating the axis' or axes' proportion of variance
# PC1 + PC2 = 50.79% of variation

# Sanity check... ^ axes % variance is same if you use rda() or prcomp()
# b.pca2<-prcomp(t(b.vst))
# summary(b.pca2)

# Now, we`ll plot our results with the plot function
#plot(b.pca)
#plot(b.pca, display = "sites", type = "points")
plot(b.pca, display = "sites", type = "text")

# Pull out PC1 and PC2
b.PC1 = data.frame(scores(b.pca, choices=c(1), display=c("sites")))
head(b.PC1)
b.PC1$SampleID<-rownames(b.PC1)

##changing the choices gives us different axes
b.PC2 = data.frame(scores(b.pca, choices=c(2), display=c("sites")))
head(b.PC2)
b.PC2$SampleID<-rownames(b.PC2)

## Create dataframe of PC axes for visualization
head(metadata)

b.pcs<-merge(b.PC1,b.PC2, by="SampleID")
b.pca.df<-merge(b.pcs, metadata, by="SampleID")
head(b.pca.df)

#### Visualize 16S PCA ####
head(b.pca.df)
sum((as.vector(b.pca$CA$eig)/sum(b.pca$CA$eig))[1]) # PC1 - 37.07%
sum((as.vector(b.pca$CA$eig)/sum(b.pca$CA$eig))[2]) # PC2 - 13.72%

b.f1<-ggplot(b.pca.df, aes(x=PC1, y=PC2)) +geom_point(aes(color=factor(Exposed)), size=3)+theme_bw()+
  labs(title="PCA: Bacteria/Archaea",subtitle="Using Variance Stabilizing Transformed (VST) Data",xlab="Axis 1", ylab="Axis 2",color="Exposure Material")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Exposure Material", labels=c("Alt"="Alternaria", "Con"="Control", "Silica"="Silica"),
                     values=unique(b.pca.df$color[order(b.pca.df$Exposed)])) +
  xlab("PC 1 [37.07 %]") + ylab("PC 2 [13.72 %]")

ggsave(b.f1,filename = "figures/16S_PCA_VST.transformed_6.21.21.pdf", width=10, height=6, dpi=600)

b.f2<-ggplot(b.pca.df, aes(x=PC1, y=PC2, label=SampleID)) +geom_text(aes(color=factor(Exposed)), size=4)+theme_bw()+
  labs(title="PCA: Bacteria/Archaea",subtitle="Using Variance Stabilizing Transformed (VST) Data",xlab="Axis 1", ylab="Axis 2",color="Exposure Material")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Exposure Material", labels=c("Alt"="Alternaria", "Con"="Control", "Silica"="Silica"),
                     values=unique(b.pca.df$color[order(b.pca.df$Exposed)])) +
  xlab("PC 1 [37.07 %]") + ylab("PC 2 [13.72 %]")

ggsave(b.f2,filename = "figures/16S_PCA.labeled_VST.transformed_6.21.21.pdf", width=15, height=6, dpi=600)


#### ITS2 PCA ####
head(f.vst)
f.pca<-rda(t(f.vst), scale=FALSE) # PCA using transposed VST transformed data
barplot(as.vector(f.pca$CA$eig)/sum(f.pca$CA$eig))

# Calculate the percent of variance explained by first two axes
sum((as.vector(f.pca$CA$eig)/sum(f.pca$CA$eig))[1:2]) # 32.79%
## cumulative proportion is calculating the axis' or axes' proportion of variance
# PC1 + PC2 = 32.79% of variation

# Sanity check... ^ axes % variance is same if you use rda() or prcomp()
#f.pca2 <- prcomp(t(f.vst))
#summary(f.pca2)

# Now, we`ll plot our results with the plot function
#plot(f.pca)
#plot(f.pca, display = "sites", type = "points")
plot(f.pca, display = "sites", type = "text")

# Pull out PC1 and PC2
f.PC1 = data.frame(scores(f.pca, choices=c(1), display=c("sites")))
head(f.PC1)
f.PC1$SampleID<-rownames(f.PC1)

##changing the choices gives us different axes
f.PC2 = data.frame(scores(f.pca, choices=c(2), display=c("sites")))
head(f.PC2)
f.PC2$SampleID<-rownames(f.PC2)

## Create dataframe of PC axes for visualization
head(metadata)

f.pcs<-merge(f.PC1,f.PC2, by="SampleID")
f.pca.df<-merge(f.pcs, metadata, by="SampleID")
head(f.pca.df)

#### Visualize ITS2 PCA ####
head(f.pca.df)
sum((as.vector(f.pca$CA$eig)/sum(f.pca$CA$eig))[1]) # PC1 - 17.86%
sum((as.vector(f.pca$CA$eig)/sum(f.pca$CA$eig))[2]) # PC2 - 14.93%

f.f1<-ggplot(f.pca.df, aes(x=PC1, y=PC2)) +geom_point(aes(color=factor(Exposed)), size=3)+theme_bw()+
  labs(title="PCA: Fungi",subtitle="Using Variance Stabilizing Transformed (VST) Data",xlab="Axis 1", ylab="Axis 2",color="Exposure Material")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Exposure Material", labels=c("Alt"="Alternaria", "Con"="Control", "Silica"="Silica"),
                     values=unique(f.pca.df$color[order(f.pca.df$Exposed)])) +
  xlab("PC 1 [17.86 %]") + ylab("PC 2 [14.93 %]")

ggsave(f.f1,filename = "figures/ITS2_PCA_VST.transformed_6.21.21.pdf", width=10, height=6, dpi=600)

f.f2<-ggplot(f.pca.df, aes(x=PC1, y=PC2, label=SampleID)) +geom_text(aes(color=factor(Exposed)), size=4)+theme_bw()+
  labs(title="PCA: Fungi",subtitle="Using Variance Stabilizing Transformed (VST) Data",xlab="Axis 1", ylab="Axis 2",color="Exposure Material")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Exposure Material", labels=c("Alt"="Alternaria", "Con"="Control", "Silica"="Silica"),
                     values=unique(f.pca.df$color[order(f.pca.df$Exposed)])) +
  xlab("PC 1 [17.86 %]") + ylab("PC 2 [14.93 %]")

ggsave(f.f2,filename = "figures/ITS2_PCA.labeled_VST.transformed_6.21.21.pdf", width=15, height=6, dpi=600)


#### 16S PCoA ####
### PCoA w/ Euclidean distances = PCA!
### For categorical variables, since you can't do PCA...
### PCoA w/ Gower's dissimilarity -- can handle categorical and continuous data in same DF!!
b.euc_dist
pcoa(b.euc_dist)
cmdscale(b.euc_dist)

b.pcoa <- pcoa(b.clr.dist)
b.pcoa$values$Relative_eig

# The proportion of variances explained is in its element values$Relative_eig
#cmdscale(bac.bray, k=(nrow(bac.ASV_table)-1), eig=TRUE)

biplot(b.pcoa)

b.pcoa.vectors<-data.frame(b.pcoa$vectors)
b.pcoa.vectors$SampleID<-rownames(b.pcoa$vectors)

b.pcoa.meta<-merge(b.pcoa.vectors, metadata, by.x="SampleID", by.y="SampleID")
head(b.pcoa.meta)

#### Visualize 16S PCoA ####
head(b.pcoa.meta)
b.pcoa$values$Relative_eig

b.f1<-ggplot(b.pcoa.meta, aes(x=Axis.1, y=Axis.2)) +geom_point(aes(color=factor(ExposedTo)), size=3)+theme_bw()+
  labs(title="PCoA: Bacteria/Archaea",subtitle="Using Centered-Log Ratio Data",xlab="Axis 1", ylab="Axis 2",color="Exposure Material")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Exposure Material", labels=c("Alt"="Alternaria", "Con"="Control"),
                     values=unique(b.pcoa.meta$exp_color[order(b.pcoa.meta$ExposedTo)])) +
  xlab("Axis 1") + ylab("Axis 2")

ggsave(b.f1,filename = "figures/16S_pcoa_VST.transformed_6.20.21.pdf", width=10, height=6, dpi=600)

b.f2<-ggplot(b.pcoa.meta, aes(x=Axis.1, y=Axis.2, label=SampleID)) +geom_text(aes(color=factor(Exposed)), size=4)+theme_bw()+
  labs(title="PCoA: Bacteria/Archaea",subtitle="Using Variance Stabilizing Transformed (VST) Data",xlab="Axis 1", ylab="Axis 2",color="Exposure Material")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Exposure Material", labels=c("Alt"="Alternaria", "Con"="Control", "Silica"="Silica"),
                     values=unique(b.pcoa.meta$color[order(b.pcoa.meta$Exposed)])) +
  xlab("Axis 1") + ylab("Axis 2")

ggsave(b.f2,filename = "figures/16S_pcoa.labeled_VST.transformed_6.20.21.pdf", width=15, height=6, dpi=600)

#### ITS2 PCoA ####
### PCoA w/ Euclidean distances = PCA!
### For categorical variables, since you can't do PCA...
### PCoA w/ Gower's dissimilarity -- can handle categorical and continuous data in same DF!!
f.euc_dist
pcoa(f.euc_dist)

f.pcoa <- pcoa(f.euc_dist)
f.pcoa$values$Relative_eig
# The proportion of variances explained is in its element values$Relative_eig
#cmdscale(bac.bray, k=(nrow(bac.ASV_table)-1), eig=TRUE)

biplot(f.pcoa)

f.pcoa.vectors<-data.frame(f.pcoa$vectors)
f.pcoa.vectors$SampleID<-rownames(f.pcoa$vectors)

f.pcoa.meta<-merge(f.pcoa.vectors, metadata, by.x="SampleID", by.y="SampleID")
head(f.pcoa.meta)


#### Visualize ITS2 PCoA ####
head(f.pcoa.meta)
f.pcoa$values$Relative_eig

f.f1<-ggplot(f.pcoa.meta, aes(x=Axis.1, y=Axis.2)) +geom_point(aes(color=factor(Exposed)), size=3)+theme_bw()+
  labs(title="PCoA: Fungi",subtitle="Using Variance Stabilizing Transformed (VST) Data",xlab="Axis 1", ylab="Axis 2",color="Exposure Material")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Exposure Material", labels=c("Alt"="Alternaria", "Con"="Control", "Silica"="Silica"),
                     values=unique(f.pcoa.meta$color[order(f.pcoa.meta$Exposed)])) +
  xlab("Axis 1") + ylab("Axis 2")

ggsave(f.f1,filename = "figures/ITS2_pcoa_VST.transformed_6.20.21.pdf", width=8, height=6, dpi=600)

f.f1<-ggplot(f.pcoa.meta, aes(x=Axis.1, y=Axis.2, label=SampleID)) +geom_text(aes(color=factor(Exposed)), size=4)+theme_bw()+
  labs(title="PCoA: Fungi",subtitle="Using Variance Stabilizing Transformed (VST) Data",xlab="Axis 1", ylab="Axis 2",color="Exposure Material")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Exposure Material", labels=c("Alt"="Alternaria", "Con"="Control", "Silica"="Silica"),
                     values=unique(f.pcoa.meta$color[order(f.pcoa.meta$Exposed)])) +
  xlab("Axis 1") + ylab("Axis 2")

ggsave(f.f1,filename = "figures/ITS2_pcoa.labeled_VST.transformed_6.20.21.pdf", width=10, height=6, dpi=600)


#### 16S BetaDisper + PERMANOVA ####

b.vst.t<-t(b.vst)

b.vst.t=b.vst.t[rownames(metadata),] # reorder both dfs by column names

head(b.vst.t)

b.euc_dist <- dist(b.vst.t, method = "euclidean")

anova(betadisper(b.euc_dist, metadata$Exposed)) # p = 0.4647
# This tells us that there is NOT a difference between group dispersions, which means that we CAN trust the results of an adonis (permutational anova) test on this
# PERMANOVA assumption of homogenous within-group disperions is not met
## aka not significant betadisper result -- groups are homogeneous in dispersion rather than being heterogeneous
## the more heterogeneous within-group results exhibit, the more likely that the significant PERMANOVA results are due to the within-group differences (aka heterogenous dispersion) rather than between-group differences
# ^ more info on this here: https://www.researchgate.net/post/Betadisper_and_adonis_in_R_Am_I_interpreting_my_output_correctly

adonis2(b.euc_dist ~ Exposed*runDuration, data = metadata, permutations = 999, by='terms') # looks for interactions between predictor variables
#These models provide the significance for statistical interactions and the main effects.
## An advantage of adonis2 is that we can also test for overall model fit, setting by=NULL

#While PERMANOVA tests differences in group means (analogous to MANOVA),
## a related test called PERMDISP can be used to evaluate homogeneity of group dispersion
#(analogous to Levene's test for equal variances). The vegan function for this test is “betadisper”:
bac.disper <- betadisper(b.euc_dist, metadata$Exposed)
bac.disper

## Significant differences in homogeneities can be tested using either parametric or permutational tests,
##and parametric post hoc contrasts can also be investigated:
anova(bac.disper)
permutest(bac.disper)
TukeyHSD(bac.disper)

#### ITS2 BetaDisper + PERMANOVA ####

f.vst.t<-t(f.vst)

f.vst.t=f.vst.t[rownames(metadata),] # reorder both dfs by column names

head(f.vst.t)

f.euc_dist <- dist(f.vst.t, method = "euclidean")

anova(betadisper(f.euc_dist, metadata$Exposed)) # p = 0.4405
# This tells us that there is NOT a difference between group dispersions, which means that we CAN trust the results of an adonis (permutational anova) test on this
# PERMANOVA assumption of homogenous within-group disperions is not met
adonis2(f.euc_dist ~ Exposed*exposed_y.n, data = metadata, permutations = 999, by='terms') # looks for interactions between predictor variables
#These models provide the significance for statistical interactions and the main effects.
## An advantage of adonis2 is that we can also test for overall model fit, setting by=NULL

#While PERMANOVA tests differences in group means (analogous to MANOVA),
## a related test called PERMDISP can be used to evaluate homogeneity of group dispersion
#(analogous to Levene's test for equal variances). The vegan function for this test is “betadisper”:
its2.disper <- betadisper(f.euc_dist, metadata$Exposed)
its2.disper

## Significant differences in homogeneities can be tested using either parametric or permutational tests,
##and parametric post hoc contrasts can also be investigated:
anova(its2.disper)
permutest(its2.disper)
TukeyHSD(its2.disper)



#### Mantel Tests ####
head(b.vst)
head(f.vst)

f.vst=f.vst[,colnames(b.vst)] # reorder both dfs by column names

head(b.vst)
head(f.vst)

b.euc_dist <- dist(t(b.vst), method = "euclidean")
f.euc_dist <- dist(t(f.vst), method = "euclidean")

#mantel correlation test

man1<-mantel(f.euc_dist,b.euc_dist)
man1$statistic
man1$signif
