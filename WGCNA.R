### ----------------------------- STEP I ----------------------------------------------------
##----------------------- Data Input and Preprocessing --------------------------------------
setwd("D:/UNI/4TH YEAR/RESEARCH/Codes/Arabidopsis")
library(WGCNA)
options(stringsAsFactors = FALSE)

#Reading data
exp = read.csv('data_processed.csv',header = T) #using the averages is too noisy

#Preprocessing
header = as.character(exp[,1])
data = exp[,-1]
data = t(data)
colnames(data) = header
data = as.data.frame(data)
#head(data)
rm(header,exp)

#QC - checking for missing values
gsg = goodSamplesGenes(data, verbose = 2);
gsg$allOK

#Cluster samples
sampleTree = hclust(dist(data), method = "average")
# Plot the sample tree as a dendrogram
library(grDevices)
pdf(file = "Plots/1-sampleClustering.pdf", width = 12, height = 9)
par(cex = 0.6)
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers in Zygotic Embryogensis", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
dev.off()
#No outliers, therefore no need of removing samples

collectGarbage()

#saving data
save(data, file = "dataInput.RData")


### ---------------------------------------- STEP II --------------------------------------------------
##-------------------------- Network Construction & Module Detection ----------------------------------

rm(list = ls(all.names = TRUE))

setwd("D:/UNI/4TH YEAR/RESEARCH/Codes/Arabidopsis")
library(WGCNA)
options(stringsAsFactors = FALSE)
enableWGCNAThreads()
library(grDevices)

# Load the data saved in the first part
load(file = "dataInput.RData")

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=100, by=5))
# Call the network topology analysis function
sft = pickSoftThreshold(data, powerVector = powers, verbose = 5)

# Plot the results:
pdf(file = "Plots/2-thresholding.pdf", width = 12, height = 9)
par(mfrow = c(1,2))
cex1 = 0.9
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()

##We choose the power 16, which is the lowest power for which the scale-free topology fit
#index curve flattens out upon reaching a high value (in this case, roughly 0.90)

## One-step network construction and module detection
net = blockwiseModules(data, power = 16, corType = "pearson", networkType = "signed",
                       TOMType = "signed", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "TOM",
                       verbose = 3)

#SIGNED NETWORK is selected as it is more accurate. 
#Rationale: https://peterlangfelder.com/2018/11/25/signed-or-unsigned-which-network-type-is-preferable/

#number of modules and module sizes
table(net$colors)

#plotting modules
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
pdf(file = "Plots/3-Modules.pdf", width = 12, height = 9)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

#saving the environment
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];
save(MEs, moduleLabels, moduleColors, geneTree,
     file = "networkConstruction.RData")

### ---------------------------------------- STEP III --------------------------------------------------
##----------------- Network analysis with functional annotation and gene ontology ----------------------

#ERROR NOT RESOLVED
rm(list = ls(all.names = TRUE))

setwd("D:/UNI/4TH YEAR/RESEARCH/Codes/Arabidopsis")
library(WGCNA)
options(stringsAsFactors = FALSE)
enableWGCNAThreads()

#load data
load(file = "dataInput.RData");

# Load network data 
load(file = "networkConstruction.RData");

#Microarray details 
#BiocManager::install('ath1121501.db')
library(ath1121501.db)

annot = ath1121501CHRLOC
mapped_probes = mappedkeys(annot)
map = as.data.frame(annot[mapped_probes])
# Match probes in the data set to the probe IDs in the annotation file
probes = names(data)
probes2annot = match(probes,map$probe_id)

# Get the corresponding Locuis Link IDs
allLLIDs = map$start_location[probes2annot]

# As background in the enrichment analysis, we will use all probes in the analysis.
fileName = paste("LocusLinkIDs-all.txt", sep="");
write.table(as.data.frame(cbind(probes,allLLIDs)), file = fileName,
            row.names = FALSE, col.names = c('Probe_ID','Locus'))

#BiocManager::install("GO.db")
#BiocManager::install("AnnotationDbi")
library('GO.db')
library('AnnotationDbi')

org = "Arabidopsis thaliana"
GOenr = GOenrichmentAnalysis(moduleColors, allLLIDs, organism = org, nBestP = 10)
##ERROR - limited to particular organisms (animals)
#GEOenrichment canbe applied for user-inserted GO data


### -------------------------------------- STEP IV -----------------------------------------------------
##--------------------------------- Network visualization-----------------------------------------------

rm(list = ls(all.names = TRUE))

setwd("D:/UNI/4TH YEAR/RESEARCH/Codes/Arabidopsis")
library(WGCNA)
options(stringsAsFactors = FALSE)
enableWGCNAThreads()
library(grDevices)

#load data
load(file = "dataInput.RData")

# Load network data saved in the second part.
load(file = "networkConstruction.RData")

#Number of genes and samples
nGenes = ncol(data)
nSamples = nrow(data)

#Visualization using Heatmaps

# Calculate topological overlap anew: this could be done more efficiently by saving the TOM
# calculated during module detection, but let us do it again here.
TOM = TOMsimilarityFromExpr(data, power = 16)
save(TOM,file = 'TOM_signed.RData')
dissTOM = 1-TOM

# For reproducibility, we set the random seed
set.seed(10)
select = sample(nGenes, size = nGenes)
selectTOM = dissTOM[select, select]
# There's no simple way of restricting a clustering tree to a subset of genes, so we must re-cluster.
selectTree = hclust(as.dist(selectTOM), method = "average")
selectColors = moduleColors[select]

# Taking the dissimilarity to a power makes the plot more informative by effectively changing
# the color palette; setting the diagonal to NA also improves the clarity of the plot
plotDiss = selectTOM^7
diag(plotDiss) = NA
pdf(file = "Plots/4-Heatmap.pdf", width = 12, height = 9)
TOMplot(plotDiss, selectTree, selectColors, main = "Network heatmap plot, selected genes")
dev.off()

#Visualizing the network of eigengenes
MEs = moduleEigengenes(data, moduleColors)$eigengenes
# Plot the relationships among the eigengenes and the trait
pdf(file = "Plots/5-Network_Eigengenes.pdf", width = 12, height = 9)
par(cex = 0.9)
plotEigengeneNetworks(MEs, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), cex.lab = 0.8, 
                      xLabelsAngle = 90)
dev.off()

### --------------------------------------- STEP V ------------------------------------------------------
##---------------------------------- Exporting the Network ----------------------------------------------

rm(list = ls(all.names = TRUE))

setwd("D:/UNI/4TH YEAR/RESEARCH/Codes/Arabidopsis")
library(WGCNA)
options(stringsAsFactors = FALSE)
enableWGCNAThreads()

#load data
load(file = "dataInput.RData")

# Load network data saved in the second part.
load(file = "networkConstruction.RData")

# Load topological overlap from the previous step
load('TOM_signed.RData')

#CHECKING MODULES OF KEY GENES
genes = read.delim('gene_list_meinke.txt')
key = toupper(genes$From)
probes = names(data)

gene_module = match(key,probes)
module = moduleColors[gene_module]
sum = table(module)

#modules with more than 30 key genes
#sum[sum>30]

write.table(sum,file = 'key gene modules.txt')

#CHECKING MODULES FOR SE MARKERS
SE_markers = read.csv('SE_markers.csv',header = T)
mrk = SE_markers$Gene_ID
marker_module = match(mrk,probes)
md = as.data.frame(cbind(probes,moduleColors))
markers = na.omit(md[marker_module,])
write.csv(markers,'SE-marker-modules.csv',row.names = F)

#=====================================================================
#EXPORTING TO CYTOSCAPE
setwd("D:/UNI/4TH YEAR/RESEARCH/Codes/Arabidopsis/Export")
clr = unique(moduleColors)
probes = names(data)

# Select modules
for (modules in clr) {
  # Select module probes
  inModule = is.finite(match(moduleColors, modules));
  modProbes = probes[inModule];
  #modGenes = annot$Probe_ID[match(modProbes, annot$Probe_ID)];
  
  # Select the corresponding Topological Overlap
  modTOM = TOM[inModule, inModule]
  dimnames(modTOM) = list(modProbes, modProbes)
  modTOM = as.table(modTOM)
  write.csv(modTOM,file = paste("TOM",modules,".csv"))
  
  # Export the network into edge and node list files Cytoscape can read
  cyt = exportNetworkToCytoscape(modTOM,
                                 edgeFile = paste("CytoscapeInput-edges-", paste(modules, collapse="-"), ".txt", sep=""),
                                 #nodeFile = paste("CytoscapeInput-nodes-", paste(modules, collapse="-"), ".txt", sep=""),
                                 weighted = TRUE,
                                 threshold = 0.01,
                                 nodeNames = modProbes,
                                 nodeAttr = moduleColors[inModule]);
}


#genes in the network after thresholding
count = numeric()
mod = character()

for (i in 1: length(clr)) {
  modules = clr[i]
  path = paste("CytoscapeInput-edges-",modules,'.txt',sep = "")
  x = read.delim(path)
  dt = unique(c(as.character(x$fromNode,x$toNode)))
  count[i] = length(dt)
  mod[i] = modules
}

summary = cbind(mod,count)
write.table(summary,file = 'Gene count.txt',row.names = F)
