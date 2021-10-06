#### Network Analysis ####
# first with just OTU table ...
# following tutorial: https://kelseyandersen.github.io/NetworksPlantPathology/Microbiome_network_ICPP2018_v2.html
library(igraph)  
#install.packages("Hmisc") # uncomment this line in order to install this package
library(Hmisc)  
#install.packages("Matrix") # uncomment this line in order to install this package
library(Matrix) 
head(otu_counts_rar)
head(t.otu_counts)
head(otu_tax_ID)
head(m.otu_counts_tax)

phy_counts<- as.data.frame(dcast(m.otus_counts_taxa, variable~Phylum, value.var="value", fun.aggregate=sum)) ### 
rownames(phy_counts)<-phy_counts$variable
phy_counts$variable<-NULL
cla_counts<- as.data.frame(dcast(m.otus_counts_taxa, variable~Class, value.var="value", fun.aggregate=sum)) ### 
rownames(cla_counts)<-cla_counts$variable
cla_counts$variable<-NULL
ord_counts<- as.data.frame(dcast(m.otus_counts_taxa, variable~Order, value.var="value", fun.aggregate=sum)) ### 
rownames(ord_counts)<-ord_counts$variable
ord_counts$variable<-NULL
fam_counts<- as.data.frame(dcast(m.otus_counts_taxa, variable~Family, value.var="value", fun.aggregate=sum)) ### 
rownames(fam_counts)<-fam_counts$variable
fam_counts$variable<-NULL
gen_counts<- as.data.frame(dcast(m.otus_counts_taxa, variable~Genus, value.var="value", fun.aggregate=sum)) ### 
rownames(gen_counts)<-gen_counts$variable
gen_counts$variable<-NULL
spec_counts<- as.data.frame(dcast(m.otus_counts_taxa, variable~Species, value.var="value", fun.aggregate=sum)) ### 
rownames(spec_counts)<-spec_counts$variable
spec_counts$variable<-NULL
head(spec_counts)
head(cla_counts)

otu.table.filter <- t.otu_counts[ ,colSums(t.otu_counts) >= 10] # drop low-abundant OTUs
#print(c(ncol(t.otu_counts),"versus",ncol(otu.table.filter))) #before and after drop OTUs count results

## maybe we could use the phyla_counts (and other taxa count dfs) to do this? 8-7-2020
otu.cor <- rcorr(as.matrix(otu.table.filter), type="spearman") # Calculating the “Spearman” correlation coefficient between OTUs using the function rcorr()
## Spearman is a ranking correlation coefficient, so I think this will serve as "weights" in our network
otu.pval <- forceSymmetric(otu.cor$P) # forceSymmetric will assign Self-correlation as NA
head(otu.pval)
rownames(otu.pval) # see what rownames are for matching step in two lines
rownames(otu_tax_ID)<-otu_tax_ID$OTU_ID # change rownames to be OTU IDs in taxa table 
sel.tax <- otu_tax_ID[rownames(otu.pval),,drop=FALSE] # Select only the taxa for the filtered OTUs by using rownames of otu.pval
all.equal(rownames(sel.tax), rownames(otu.pval)) # sanity check

p.yes <- otu.pval<0.05 # Filter the association based on p-values and level of correlations
# vvvv Select the r values for the filter probality of < 0.5.
r.val = otu.cor$r # select all the correlation values 
p.yes.r <- r.val*p.yes # only select correlation values based on p-value criterion 

# Select OTus by level of correlation
p.yes.r <- abs(p.yes.r)>0.75 # output is logical vector
p.yes.rr <- p.yes.r*r.val # use logical vector for subscripting.

#Create an adjacency matrix
adjm <- as.matrix(p.yes.rr)

#Add taxonomic information from the metadata associated with adjacency matrix
colnames(adjm) <- as.vector(sel.tax$Phylum)
rownames(adjm) <- as.vector(sel.tax$Phylum)

# Creating graph object from adjacency matrix
net.grph=graph.adjacency(adjm,mode="undirected",weighted=TRUE,diag=FALSE)

# Obtaining edge weight based on the Spearman correlation
edgew<-E(net.grph)$weight

# Creating a vector to remove the isolated nodes (nodes with no interactions)
bad.vs<-V(net.grph)[degree(net.grph) == 0] 

# Removing the isolated nodes from the graph object using the function delete.vertices()
net.grph <-delete.vertices(net.grph, bad.vs)
class(net.grph)
plot(net.grph,
     vertex.size=8,
     vertex.frame.color="black",
     edge.curved=F,
     edge.width=1.5,
     layout=layout.fruchterman.reingold,
     edge.color=ifelse(edgew < 0,"red","blue"),
     #vertex.label=NA,
     vertex.label.color="black",
     vertex.label.family="Arial",
     vertex.label.font=0.1)

### Another way.... https://chrischizinski.github.io/rstats/igraph-ggplotll/
head(otus_counts_taxa)
t.otu_counts_tax<-as.data.frame(t(otus_counts_taxa))
m.otus_counts_taxa<-melt(otus_counts_taxa)
head(m.otus_counts_taxa)



nw_phy<-graph_from_incidence_matrix(phy_counts, mode='all', multiple=TRUE)
plot(nw_phy)

caught.inc <- graph.incidence(t.otu_counts_tax, weighted = TRUE)  #make data into a bipartite graph object
obs.parties.all <- bipartite.projection(caught.inc)[[1]]
obs.spp.all <- bipartite.projection(caught.inc)[[2]]

op <- par(mfrow = c(1, 2))
fr.all <- layout.fruchterman.reingold(obs.spp.all)
plot(obs.spp.all, layout = fr.all, edge.color = "black", edge.width = E(obs.spp.all)$weight * 
       0.1, vertex.label = V(obs.spp.all)$name)
obs.sg.all <- fastgreedy.community(obs.spp.all, weights = E(obs.spp.all)$weight)
plot(obs.sg.all, obs.spp.all, layout = fr.all, edge.width = E(obs.spp.all)$weight * 
       0.25, vertex.label = V(obs.spp.all)$name, vertex.label.color = "blue")
par(op)


fr.all.df <- as.data.frame(fr.all)  ## convert the layout to a data.frame
fr.all.df$Phylum <- colnames(t.otu_counts_tax)  ## add in the species codes

fr.all.df

ggplot() +
  geom_point(data=fr.all.df,aes(x=V1,y=V2),size=21,colour="black") +  # adds a black border around the nodes
  geom_point(data=fr.all.df,aes(x=V1,y=V2),size=20,colour="lightgrey") +
  geom_text(data=fr.all.df,aes(x=V1,y=V2,label=Phylum)) + # add the node labels
  scale_x_continuous(expand=c(0,1))+  # expand the x limits 
  scale_y_continuous(expand=c(0,1))+ # expand the y limits
  theme_bw()  # use the ggplot black and white theme
