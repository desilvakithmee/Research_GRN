### ----------------------------- STEP I --------------------------------------------
##----------------------- Data Input and Preprocessing ------------------------------

setwd("D:/UNI/4TH YEAR/RESEARCH/Codes/Arabidopsis")
library(WGCNA)
options(stringsAsFactors = FALSE)

#Reading data
exp = read.csv('data_processed.csv',header = T) 
#exp = read.csv('data_processed_avg.csv', header = T) #using the averages is too noisy

#Preprocessing
header = as.character(exp[,1])
data = exp[,-1]
data = t(data)
colnames(data) = header
data = as.data.frame(data)
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
plot(sampleTree, main = "Sample clustering to detect outliers in Zygotic Embryogensis",
     sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
dev.off()
#No outliers, therefore no need of removing samples

collectGarbage()

#saving data
save(data, file = "dataInput.RData")


### --------------------------------- STEP II ---------------------------------------
##-------------------- Network Construction & Module Detection -----------------------

rm(list = ls(all.names = TRUE))

setwd("D:/UNI/4TH YEAR/RESEARCH/Codes/Arabidopsis")
library(WGCNA)
options(stringsAsFactors = FALSE)
enableWGCNAThreads()
library(grDevices)

# Load the data saved in the first part
load(file = "dataInput.RData")

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=30, by=5))
# Call the network topology analysis function
sft = pickSoftThreshold(data, powerVector = powers, verbose = 5)

# Plot the results:
pdf(file = "Plots/2-thresholding.pdf", width = 12, height = 9)
par(mfrow = c(1,2))
cex1 = 0.9
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",
     type="n", main = paste("Scale independence"))
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

##We choose the power 17, which is the lowest power for which the scale-free topology 
#fit index curve flattens out upon reaching a high value (in this case, roughly 0.90)

## One-step network construction and module detection
net = blockwiseModules(data, power = 17, corType = "pearson", networkType = "signed",
                       TOMType = "signed", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "TOM",
                       verbose = 3)
 
#SIGNED NETWORK is selected as it is more accurate. 
#https://peterlangfelder.com/2018/11/25/signed-or-unsigned-which-network-type-is-preferable/

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
MEs0 = moduleEigengenes(data, moduleColors)$eigengenes
MEs = removeGreyME(MEs0,  greyMEName = paste(moduleColor.getMEprefix(), "grey", sep=""))
geneTree = net$dendrograms[[1]]
save(MEs, moduleLabels, moduleColors, geneTree,
     file = "networkConstruction.RData")

### --------------------------------- STEP III -----------------=======---------------
#------------------------------ Module-Trait Relationship ----------------------------

rm(list = ls(all.names = TRUE))

setwd("D:/UNI/4TH YEAR/RESEARCH/Codes/Arabidopsis")
library(WGCNA)
options(stringsAsFactors = FALSE)
library(grDevices)

#load data
load(file = "dataInput.RData")
load(file = 'dataInput_avg.RData')

# Load network data 
load(file = "networkConstruction.RData")

# Define numbers of genes and samples
nGenes = ncol(data)
nSamples = nrow(data)

# Recalculate MEs with color
MEs0 = moduleEigengenes(data, moduleColors)$eigengenes
MEs0 = removeGreyME(MEs0,greyMEName = "MEgrey")
MEs = orderMEs(MEs0)
ME_names = colnames(MEs)
Zygotic = colMeans(MEs[1:2,])
Day7 = colMeans(MEs[3:4,])
Day14 = colMeans(MEs[5:6,])
Somatic = colMeans(MEs[7:8,])
MEs = data.frame(rbind(Zygotic,Day7,Day14,Somatic))
colnames(MEs) = ME_names
rm(Zygotic,Day14,Day7,Somatic)
head(MEs)

#trait data
trait = matrix(c(1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1),byrow = T,nrow = 4,ncol = 4)
colnames(trait) = row.names(MEs)
rownames(trait) = row.names(MEs)
trait = as.data.frame(trait)
trait

moduleTraitCor = cor(MEs,trait, use = "complete.obs")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

pdf(file = "Plots/6-Module-trait relationship.pdf", width = 12, height = 9, 
    pagecentre = T, paper = "a4r")

# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(4,8,2,2))
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = c("Stage I","Stage II","Stage III","Stage IV"),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               xLabelsAngle = 45,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.8,
               zlim = c(-1,1),
               plotLegend = TRUE,
               main = paste("Module-trait relationships"))
dev.off()


# Define variable weight containing the weight column of datTrait
weight = as.data.frame(trait$Somatic)
names(weight) = "weight"
# names (colors) of the modules
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(data, MEs0, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")
geneTraitSignificance = as.data.frame(cor(dt, weight, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) = paste("GS.", names(weight), sep="")
names(GSPvalue) = paste("p.GS.", names(weight), sep="")
head(GSPvalue)

cor = abs(moduleTraitCor) >= 0.5
p = moduleTraitPvalue < 0.05
xx = cor & p
write.csv(xx, 'Module-Trait.csv', quote = F)

#annot = read.csv('annotation.csv')
annot = annot[,1:6]
clr = sort(unique(moduleColors))
clr = clr[clr != "grey"]
hubs = character()

for (i in 1:length(clr)) {
  dt = geneModuleMembership[,i]
  kk = order(dt, decreasing = T)
  nm = rownames(geneModuleMembership)[kk]
  top10 = nm[1:10]  
  hubs = cbind(hubs,top10)
  #mt = match(top10,annot$Probe.Set.ID)
  #kg = annot[mt,]
  #file = paste(clr[i],' hub genes annot.csv', sep = "")
  #write.csv(kg,file,row.names = F,quote = F)
  
}
hubs = as.data.frame(hubs)
colnames(hubs) = clr
write.csv(hubs,'Top 10 hub genes.csv', row.names = F, quote = F)


setwd("D:/UNI/4TH YEAR/RESEARCH/Codes/Arabidopsis/Export")
options(stringsAsFactors = FALSE)
hubs = read.csv('Top 10 hub genes.csv')
gene = read.csv("gene_description.csv", header = T)
clr = colnames(hubs)
sum = data.frame()

for (i in 1:length(clr)) {
  gg = hubs[,i]
  key = match(gg, gene$ï..Column1)
  des = gene[key,]
  des[,5] = clr[i]
  sum = rbind.data.frame(sum,des)
}
write.csv(sum,'10 hub genes.csv', row.names = F, quote = F)
write.table(sum,'10 hub genes.txt', row.names = F, quote = F, sep = "\t")
#================== Plots =======================================================

#Module-trait siginificant
rm(list = ls(all.names = TRUE))
setwd("D:/UNI/4TH YEAR/RESEARCH/Codes/Arabidopsis")
library(gplots)

mt = read.csv('Module-trait.csv',header = T)
rownames(mt) = mt$X
mt[,1] = NULL
mt = as.matrix(mt)
dt = as.integer(mt)
dt = matrix(dt, nrow = 14,ncol = 4)
colnames(dt) = colnames(mt)
rownames(dt) = rownames(mt)
#dt = t(dt)
dt[c(10,14),1] = -1
dt[c(6,7),4] = -1

par(mar=c(7,4,4,2))
heatmap.2(dt, col = cm.colors(3),  trace = "none", Colv = F, colsep=0:ncol(dt),
          rowsep=0:nrow(dt), sepcolor = c("lightgrey"),margins=c(8,8),
          dendrogram = "none", Rowv = F)
dev.off()

#================== Plots =======================================================

abs(moduleTraitCor) > 0.6
moduleTraitPvalue < 0.05

setwd("D:/UNI/4TH YEAR/RESEARCH/Codes/Arabidopsis/Plots/figures")

clr1 = c('green','magenta','yellow','brown','cyan','purple','red','salmon',
         'greenyellow','pink','turquoise')
clr2 = c('brown','cyan','tan','midnightblue')
clr3 = c('yellow','black','tan','greenyellow')
clr4 = c('blue','green','purple','red')

for (module in clr4){
file = paste("Stage IV - ",module,".jpg", sep = "")
column = match(module, modNames);
moduleGenes = moduleColors==module;
jpeg(file, width = 900, height = 600 )
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for Stage 4",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
dev.off()
}

library(ggplot2)
library(gridExtra)
module = rownames(moduleTraitCor)
mod = data.frame(module,moduleTraitCor)
rownames(mod) = NULL

clr = c('grey23','royalblue2','tan3','cyan2','limegreen','greenyellow','deeppink3',
        'hotpink1','mediumpurple2','red3','salmon','tan','turquoise3','gold')


 q1 = ggplot(data = mod,aes(x = module, y = Zygotic, fill = module))+
  geom_col()+
  scale_fill_manual(values = clr, guide = F) +
  scale_x_discrete(labels = NULL) +
  #theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5) ) +
  labs(title = "Module Expression in Stage I") +
  ylab("Module-trait relationship") +
  xlab("Modules")
 
q2 = ggplot(data = mod,aes(x = module, y = Day7, fill = module))+
  geom_col()+
  scale_fill_manual(values = clr, guide = F) +
  scale_x_discrete(labels = NULL) +
  labs(title = "Module Expression in Stage II") +
  ylab("Module-trait relationship") +
  xlab("Modules")
q3 = ggplot(data = mod,aes(x = module, y = Day14, fill = module))+
  geom_col()+
  scale_fill_manual(values = clr, guide = F) +
  scale_x_discrete(labels = NULL) +
  labs(title = "Module Expression in Stage III") +
  ylab("Module-trait relationship") +
  xlab("Modules")
q4 = ggplot(data = mod,aes(x = module, y = Somatic, fill = module))+
  geom_col()+
  scale_fill_manual(values = clr, guide = F) +
  scale_x_discrete(labels = NULL) +
  labs(title = "Module Expression in Stage IV") +
  ylab("Module-trait relationship") +
  xlab("Modules")

grid.arrange(q1, q2,q3,q4, nrow = 2)


### ------------------------------- STEP IV -----------------------------------------
##----------------------- Network visualization -------------------------------------

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

# Calculate topological overlap 
#TOM = TOMsimilarityFromExpr(data, power = 17)
#save(TOM,file = 'TOM_signed.RData')
load('TOM_signed.RData')
dissTOM = 1-TOM

# For reproducibility, we set the random seed
set.seed(10)
nSize = 1000
select = sample(nGenes, size = nSize)
selectTOM = dissTOM[select, select]
selectTree = hclust(as.dist(selectTOM), method = "average")
selectColors = moduleColors[select]

plotDiss = selectTOM^7           #makes the plot more informative
diag(plotDiss) = NA             #improves the clarity of the plot
pdf(file = "Plots/4-Heatmap.pdf", width = 12, height = 9)
TOMplot(plotDiss, selectTree, selectColors, 
        main = "Network heatmap plot, selected genes")
dev.off()

#Visualizing the network of eigengenes
MEs = moduleEigengenes(data, moduleColors)$eigengenes
MEs = removeGreyME(MEs,greyMEName = "MEgrey")
# Plot the relationships among the eigengenes and the trait
pdf(file = "Plots/5-Network_Eigengenes.pdf", width = 12, height = 9)
par(cex = 0.9)
plotEigengeneNetworks(MEs, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), 
                      cex.lab = 0.8, 
                      xLabelsAngle = 90)
dev.off()

### ------------------------------- STEP V -------------------------------------------
##------------------------ Exporting the Network -------------------------------------

rm(list = ls(all.names = TRUE))

setwd("D:/UNI/4TH YEAR/RESEARCH/Codes/Arabidopsis")
library(WGCNA)
options(stringsAsFactors = FALSE)
enableWGCNAThreads()

#load data
load(file = "dataInput.RData")
load(file = "networkConstruction.RData")
load('TOM_signed.RData')

#Hub genes
hub = chooseTopHubInEachModule(data, moduleColors, omitColors = "grey", 
                               power = 17, type = "signed")
write.table(hub,'Hub genes.txt',quote = F, col.names = F)



#EXPORTING TO CYTOSCAPE
setwd("D:/UNI/4TH YEAR/RESEARCH/Codes/Arabidopsis/Export")
clr = unique(moduleColors)
clr = sort(clr)
probes = names(data)

# Select modules
for (modules in clr) {
# Select module probes
inModule = is.finite(match(moduleColors, modules));
modProbes = probes[inModule];

# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule]
dimnames(modTOM) = list(modProbes, modProbes)
modTOM = as.table(modTOM)
write.csv(modTOM,file = paste("TOM",modules,".csv"))

# Export the network into edge and node list files that Cytoscape can read
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("CytoscapeInput-edges-", 
                                                paste(modules, collapse="-"), 
                                                ".txt", sep=""),
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
write.table(summary,file = 'Gene count.txt',row.names = F,quote = F)
head(summary)


### ------------------------------- STEP VI ------------------------------------------
##-------------------------- Validating the Network ----------------------------------

rm(list = ls(all.names = TRUE))

setwd("D:/UNI/4TH YEAR/RESEARCH/Codes/Arabidopsis")
library(WGCNA)
options(stringsAsFactors = FALSE)
enableWGCNAThreads()

#load data
load(file = "dataInput.RData")
load(file = "networkConstruction.RData")

softPower = 17
adjacency = adjacency(data, power = softPower)

setLabels = c("Network", "Test")
multiExpr = list(Network = list(data = adjacency), Test = list(data = data))
multiColor = list(Network = moduleColors, Test = moduleColors);
nSets = 2

system.time( {
  mp = modulePreservation(multiExpr, multiColor,
                          referenceNetworks = c(1:2),
                          nPermutations = 200,
                          randomSeed = 1,
                          verbose = 3)
} );
# Save the results
save(mp, file = "ModulePreservation.RData")

ref = 1
test = 2
statsObs = cbind(mp$quality$observed[[ref]][[test]][, -1], mp$preservation$observed[[ref]][[test]][, -1])
statsZ = cbind(mp$quality$Z[[ref]][[test]][, -1], mp$preservation$Z[[ref]][[test]][, -1]);
print( cbind(statsObs[, c("medianRank.pres", "medianRank.qual")],
             signif(statsZ[, c("Zsummary.pres", "Zsummary.qual")], 2)) )
mp_sum = cbind(statsObs[, c("medianRank.pres", "medianRank.qual")],
               signif(statsZ[, c("Zsummary.pres", "Zsummary.qual")], 2))
write.csv(mp_sum,'Module Preservation.csv',quote = F)

# Module labels and module sizes are also contained in the results
modColors = rownames(mp$preservation$observed[[ref]][[test]])
moduleSizes = mp$preservation$Z[[ref]][[test]][, 1];
# leave grey and gold modules out
plotMods = !(modColors %in% c("grey", "gold"));
# Text labels for points
text = modColors[plotMods];
# Auxiliary convenience variable
plotData = cbind(mp$preservation$observed[[ref]][[test]][, 2], mp$preservation$Z[[ref]][[test]][, 2])
# Main titles for the plot
mains = c("Preservation Median rank", "Preservation Zsummary");
# Start the plot
sizeGrWindow(10, 5);
library(grDevices)
pdf(file = "Plots/7-modulePreservation.pdf", width = 12, height = 9)
par(mfrow = c(1,2))
par(mar = c(4.5,4.5,2.5,1))
for (p in 1:2)
{
  min = min(plotData[, p], na.rm = TRUE);
  max = max(plotData[, p], na.rm = TRUE);
  # Adjust ploting ranges appropriately
  if (p==2)
  {
    if (min > -max/10) min = -max/10
    ylim = c(min - 0.1 * (max-min), max + 0.1 * (max-min))
  } else
    ylim = c(max + 0.1 * (max-min), min - 0.1 * (max-min))
  plot(moduleSizes[plotMods], plotData[plotMods, p], col = 1, bg = modColors[plotMods], pch = 21,
       main = mains[p],
       cex = 2.4,
       ylab = mains[p], xlab = "Module size", log = "x",
       ylim = ylim,
       xlim = c(10, 2000), cex.lab = 1.2, cex.axis = 1.2, cex.main =1.4)
  labelPoints(moduleSizes[plotMods], plotData[plotMods, p], text, cex = 1, offs = 0.08);
  # For Zsummary, add threshold lines
  if (p==2)
  {
    abline(h=0)
    abline(h=2, col = "blue", lty = 2)
    abline(h=10, col = "darkgreen", lty = 2)
  }
}
dev.off()
