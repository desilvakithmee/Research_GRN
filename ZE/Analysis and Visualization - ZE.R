
#ANALYSIS AND VISUALIZATION OF THE NETWORK MODULES

#=================================================================================

## Checking the presence of key genes and markers in the modules

rm(list = ls(all.names = TRUE))

setwd("D:/UNI/4TH YEAR/RESEARCH/Codes/Zygotic_Embryogenesis")

#load data
load(file = "dataInput.RData")
load(file = "networkConstruction.RData")

#CHECKING MODULES OF KEY GENES
genes = read.delim('Key genes.txt', header = F)
key = toupper(genes$V1)
probes = names(data)

gene_module = match(key,probes)
module = moduleColors[gene_module]
sum = table(module)

write.table(sum,file = 'key gene modules.txt', row.names = F, quote = F)

md = cbind(probes,moduleColors)
key_genes = na.omit(md[gene_module,])
write.csv(key_genes,'Keygenes-modules.csv',row.names = F)

#=================================================================================

#Correlational Analysis of Network Modules

setwd("D:/UNI/4TH YEAR/RESEARCH/Codes/Zygotic_Embryogenesis")
library(WGCNA)
# Load network data 
load(file = "dataInput.RData")
load(file = "networkConstruction.RData")

MEs = moduleEigengenes(data, moduleColors)$eigengenes
MEs = removeGreyME(MEs,greyMEName = "MEgrey")


library(corrr)
library(tidyverse)
res.cor <- correlate(MEs)
name = names(res.cor)
name = gsub('ME',"",name)
name = c(name[1],str_to_title(name[-1]))
colnames(res.cor) = name
res.cor[,1] = name[-1]

#rplot(res.cor, colors = c("red","green"),shape = 19,legend = T,print_cor = T)

res.cor = res.cor[,-1]
row.names(res.cor) = name[-1]
res.cor = as.matrix(res.cor)

flattenCorrMatrix = function(cormat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = colnames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut]
  )
}

mod_corr = flattenCorrMatrix(res.cor)
cor_flt = mod_corr[abs(mod_corr$cor) > 0.5,]
write.table(cor_flt,'Module network.txt', quote = F, row.names = F, sep = ",")

library(ggplot2)
library(tidyverse)
setwd("D:/UNI/4TH YEAR/RESEARCH/Codes/Zygotic_Embryogenesis/Export")
data = read.delim('Gene count.txt', header = T, sep = " ")
data$mod = str_to_title(data$mod)
colnames(data) = c('Module', 'Gene frequency')

ggplot(data)  +
  geom_col(aes(x = data$Module, y = data$`Gene frequency`, fill = Module)) +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5) ) +
  labs(title = "Gene Count in Modules") +  
  scale_fill_manual(values = data$Module) +
  xlab("Module") +
  ylab("Gene Frequency")
ggsave('Gene count.jpg')

#===========================================================================
setwd("D:/UNI/4TH YEAR/RESEARCH/Codes/Zygotic_Embryogenesis")
options(stringsAsFactors = F)
library(WGCNA)

dt = read.csv('data_filtered.csv')
load(file = "dataInput.RData")
load(file = "networkConstruction.RData")
clr = 

data2 = sample(data,6000, replace = F)
mp = modulePreservation(dt, moduleColors, networkType = "signed",
                        referenceNetworks = 1, verbose = 4, maxModuleSize = 500,
                        maxGoldModuleSize = 500, nPermutations = 200)


#===========================================================================

#To get the genes under each module and their functions

setwd("D:/UNI/4TH YEAR/RESEARCH/Codes/Zygotic_Embryogenesis/Export/Modular genes")
options(stringsAsFactors = F)

annot = read.csv('annotation.csv', header = T)
Gene_ID = as.character(toupper(annot$Probe.Set.ID))
Gene_name = as.character(annot$Gene.Symbol)
Gene_function = as.character(annot$Gene.Title)

gene_file = data.frame(Gene_ID,Gene_name,Gene_function,stringsAsFactors = F)
write.csv(gene_file,'Gene details.csv', quote = F, row.names = F)

#---------------------------------------------------------------------------

rm(list = ls(all.names = TRUE))

setwd("D:/UNI/4TH YEAR/RESEARCH/Codes/Zygotic_Embryogenesis/Export/Modular genes")

gene_file = read.csv('Gene descriptors.csv',header = T, stringsAsFactors = F,
                     sep = ",")
load(file = "networkConstruction.RData")

clr = unique(moduleColors)
clr = stringi::stri_trans_totitle(clr)
clr = sort(clr)
clr = clr[-6]

for (i in 1: length(clr)) {
  modules = clr[i]
  path = paste(modules,' genes.txt',sep = "")
  x = read.delim(path, header = F)
  key = match(x$V1,gene_file$Locus.Identifier)
  genes = na.omit(gene_file[key,])
  path2 = paste(modules,' gene table.csv',sep = "")
  write.csv(genes, path2, quote = F, row.names = F)
}

#=============================================================================

#Modular Analysis
setwd("D:/UNI/4TH YEAR/RESEARCH/Codes/Zygotic_Embryogenesis/Export/Modular genes")
library(tidyverse)

#Go_summ.csv is a summary of GO terms common to multiple modules
data = read.csv('GO_summ.csv', header = T)
data[,1] =str_to_title(as.character(data$module))
data$Term = str_to_sentence(as.character(data$go_term))


data$prop = data$num/data$total
dt = data[order(data$Term),]
#plot_data = dt[c(4:5,8:13,16:19,22:23,27:31),]
#clr = sort(unique(plot_data$X))
clr = c('lightskyblue','burlywood','palegreen','pink','aquamarine1','lightgoldenrod1')


library(ggplot2)
ggplot(data = dt, aes(y = dt$prop, x = dt$Term, 
                                  fill = dt$module, width = 0.7))+
  geom_bar(position = position_dodge2(width = 0.9, preserve = "single"), 
           stat="identity", colour = "gray50") +
  coord_flip() + 
  scale_fill_manual(values = clr) + 
  scale_linetype_manual(values = 'white') +
  theme(panel.background = element_rect(fill = "gray95", colour = "gray95", 
                                        size = 0.5, linetype = "solid"))+
  labs(title = "GO terms related to zygotic embryogenesis in modules") +
  ylab("Gene proportion") +
  xlab("Gene Ontology term")+
  labs(fill = "Module")

#============================================================================
#Hormone-related genes in the network modules

rm(list = ls(all.names = TRUE))

setwd("D:/UNI/4TH YEAR/RESEARCH/Codes/Zygotic_Embryogenesis")
options(stringsAsFactors = FALSE)

#load data
load(file = "dataInput.RData")
load(file = "networkConstruction.RData")

#CHECKING MODULES OF GENES
genes = read.csv('hormone_related_genes.csv',header = T)
key = toupper(genes$Gene)
probes = names(data)

dt = data.frame(probes, moduleColors)
gene_module = match(key,probes)
func = dt[gene_module,]
dt_fun = na.omit(cbind(func,genes))

write.csv(dt_fun,'Hormone related genes in modules.csv',row.names = F,quote = F)

#===========================================================================
rm(list = ls(all.names = TRUE))

setwd("D:/UNI/4TH YEAR/RESEARCH/Codes/Zygotic_Embryogenesis")
options(stringsAsFactors = FALSE)

#Heatmap of hormone related genes

dt = read.csv('Hormone-module.csv', header = T)
dt[is.na(dt)] = 0
dt = unique(dt)
rownames(dt) = dt$gene
colnames(dt) = stringr::str_to_title(colnames(dt))
dt = dt[,-1]
dt = as.matrix(dt)
dt1 = dt[1:31,]
dt2 = dt[c(32:51,59:65),]
dt3 = dt[c(52:58,65:77),]

library(pheatmap)
pheatmap(t(dt1), treeheight_row = 0, treeheight_col = 0,cellwidth = 15, 
         cellheight = 18, color = rev(brewer.pal(8,"PuBu")), scale = "none",
         cluster_rows = F, cluster_cols = F, angle_col = 45)
pheatmap(t(dt2), treeheight_row = 0, treeheight_col = 0,cellwidth = 15, 
         cellheight = 18, color = rev(brewer.pal(8,"PuBu")), scale = "none",
         cluster_rows = F, cluster_cols = F, angle_col = 45)
pheatmap(t(dt3), treeheight_row = 0, treeheight_col = 0,cellwidth = 18, 
         cellheight = 18, color = rev(brewer.pal(8,"PuBu")), scale = "none",
         cluster_rows = F, cluster_cols = F, angle_col = 45)

#Heatmap of enzyme related genes

dt = read.csv('Enzyme-module.csv', header = T)
dt[is.na(dt)] = 0
dt = unique(dt)
rownames(dt) = dt$gene
colnames(dt) = stringr::str_to_title(colnames(dt))
dt = dt[,-1]
dt = as.matrix(dt)

library(pheatmap)
pheatmap(t(dt), treeheight_row = 0, treeheight_col = 0,cellwidth = 20, 
         cellheight = 18, color = rev(brewer.pal(8,"PuBu")), scale = "none",
         cluster_rows = F, cluster_cols = F, angle_col = 45, legend = F)

#===========================================================================
rm(list = ls(all.names = TRUE))

setwd("D:/UNI/4TH YEAR/RESEARCH/Codes/Zygotic_Embryogenesis/Export/Modular genes")
options(stringsAsFactors = FALSE)

load(file = "dataInput.RData")
load(file = "networkConstruction.RData")
key = read.delim('Key genes.txt', header = F)

key = key$V1
probes = names(data)

gene_module = match(key,probes)
md = as.data.frame(cbind(probes,moduleColors))
key_genes = na.omit(md[gene_module,])
write.csv(key_genes,'Keygenes-modules.csv',row.names = F)

#------------------------------------------------------------------

rm(list = ls(all.names = TRUE))

dt = read.csv('Keygenes-modules.csv')
clr = sort(unique(dt$moduleColors))
sum = data.frame()

thresh = c(0.6,rep(0.4,7),0.62,0.4)

for (c in 1:length(clr)) {
  mod = clr[c]
  exp = dt[dt$moduleColors==mod,]
  edges = read.delim(paste('CytoscapeInput-edges-',mod,'.txt',sep = ""))
  th = thresh[c]
  
  for (i in 1:dim(exp)[1]) {
    gene = exp[i,1]
    rcd = edges[(edges$fromNode == gene | edges$toNode == gene)& edges$weight>th,1:3]
    if (dim(rcd)[1] != 0)
      rcd$module = mod
      sum = data.frame(rbind(sum,rcd))
    }
  
}
table(sum$module)
write.csv(sum,'Key genes summary.csv',quote = F,row.names = F)

#----------------------------------------------------------------------

rm(list = ls(all.names = TRUE))
setwd("D:/UNI/4TH YEAR/RESEARCH/Codes/Zygotic_Embryogenesis/Export/Modular genes")
options(stringsAsFactors = FALSE)

annot = read.csv('annotation.csv')
key = read.csv('Keygenes-modules.csv')

mt = match(key$probes,annot$Probe.Set.ID)
kg = annot[mt,]
write.csv(kg,'key gene - annot.csv',row.names = F,quote = F)

