
#ANALYSIS AND VISUALIZATION OF THE NETWORK MODULES

#=================================================================================

## Checking the presence of key genes and markers in the modules

rm(list = ls(all.names = TRUE))

setwd("D:/UNI/4TH YEAR/RESEARCH/Codes/Arabidopsis")

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

write.table(sum,file = 'key gene modules.txt')

#CHECKING MODULES FOR SE MARKERS
SE_markers = read.csv('SE_markers.csv',header = T)
mrk = SE_markers$Gene_ID
marker_module = match(mrk,probes)
md = as.data.frame(cbind(probes,moduleColors))
markers = na.omit(md[marker_module,])
write.csv(markers,'SE-marker-modules.csv',row.names = F)

key_genes = na.omit(md[gene_module,])
write.csv(key_genes,'Keygenes-modules.csv',row.names = F)

#=================================================================================

#Correlational Analysis of Network Modules

setwd("D:/UNI/4TH YEAR/RESEARCH/Codes/Arabidopsis")
library(WGCNA)
# Load network data 
load(file = "dataInput.RData")
load(file = "networkConstruction.RData")

MEs = moduleEigengenes(data, moduleColors)$eigengenes
MEs = removeGreyME(MEs,greyMEName = "MEgrey")

top = chooseTopHubInEachModule(data, moduleColors, power = 2, type = 'signed')
hub = chooseOneHubInEachModule(data, moduleColors, numGenes = 100, power = 2, 
                               type = "signed")

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
cor_flt = rbind(mod_corr[mod_corr$cor > 0.5,],mod_corr[mod_corr$cor < -0.5,])
write.table(cor_flt,'Module network.txt', quote = F, row.names = F, sep = ",")

library(ggplot2)
library(tidyverse)
setwd("D:/UNI/4TH YEAR/RESEARCH/Codes/Arabidopsis/Export")
data = read.delim('Gene count.csv', header = T, sep = ",")
data$mod = str_to_title(data$mod)
colnames(data) = c('Module', 'Gene frequency')

qq = ggplot(data)  +
  geom_col(aes(x = data$Module, y = data$`Gene frequency`, fill = Module)) +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5) ) +
  labs(title = "Gene Count in Modules") +  
  scale_fill_manual(values = data$Module) +
  xlab("Module") +
  ylab("Gene Frequency")
ggsave('Gene count.jpg')
#===========================================================================

#To get the genes under each module and their functions

setwd("D:/UNI/4TH YEAR/RESEARCH/Codes/Arabidopsis/Export/Modular genes")

annot = read.csv('annotation.csv',header = T,sep = ",")
key = read.delim('ID.txt', header = F)
mt = match(annot$Probe.Set.ID,key$V1)
gene_ID = toupper(key[mt,3])
annot$Probe.Set.ID = gene_ID
Gene_ID = as.character(toupper(annot$Representative.Public.ID))
Gene_name = as.character(annot$Gene.Symbol)
Gene_function = as.character(annot$Gene.Title)

gene_file = data.frame(Gene_ID,Gene_name,Gene_function,stringsAsFactors = F)
write.csv(gene_file,'Gene details.csv', quote = F, row.names = F)

#---------------------------------------------------------------------------

rm(list = ls(all.names = TRUE))

setwd("D:/UNI/4TH YEAR/RESEARCH/Codes/Arabidopsis/Export/Modular genes")

gene_file = read.csv('Gene descriptors.csv',header = T, stringsAsFactors = F,
                     sep = ",")
load(file = "networkConstruction.RData")

clr = unique(moduleColors)
clr = clr[-14]
clr = stringi::stri_trans_totitle(clr)

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
setwd("D:/UNI/4TH YEAR/RESEARCH/Codes/Arabidopsis/Export/Modular genes")
library(tidyverse)

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
  labs(title = "GO terms related to somatic embryogenesis in modules") +
  ylab("Gene proportion") +
  xlab("Gene Ontology term")+
  labs(fill = "Module")

#============================================================================
#Hormone-related genes in the network modules

rm(list = ls(all.names = TRUE))

setwd("D:/UNI/4TH YEAR/RESEARCH/Codes/Arabidopsis")
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

#Plots to visualize 
rm(list = ls(all.names = TRUE))

library(ggplot2)
library(tidyverse)

#Hormone related genes
dt = read.csv('Hormone related genes in modules.csv', header = T)
dt$moduleColors = str_to_title(dt$moduleColors)
dt$Related.Hormone.Type = str_to_title(dt$Related.Hormone.Type)
clr = c('grey23','royalblue2','tan3','limegreen','olivedrab1','deeppink3',
        'mediumpurple2','red3','turquoise3','gold')

ggplot(data = dt)+
  geom_bar(aes(x=dt$Related.Hormone.Type,group = dt$moduleColors, 
               fill = dt$moduleColors), 
           position = position_dodge2(width = 0.9, preserve = "single")) +
  scale_fill_manual(values = clr, name = "Module") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10.5), 
        legend.position = "bottom") +
  labs(title = "Plant growth regulator related genes in modules") +
  ylab("Gene count") +
  xlab("Plant growth regulator")

title = unique(dt$Related.Hormone.Type)
path = "D:/UNI/4TH YEAR/RESEARCH/Codes/Arabidopsis/Plots/figures"

for (i in 1:length(title)) {
  dth = dt[dt$Related.Hormone.Type == title[i],]
  
  ggplot(data = dth) +
    geom_dotplot(aes(x = dth$moduleColors, y = dth$Related.Hormone.Type, 
                 fill = dth$Function.category)) +
                   xlab("") +
                   ylab("")
  file = paste("Plot ",title[i],".jpg",sep = "")
  ggsave(filename = file,path = path)
  
}


#Enzyme related genes
dt_enz = read.csv('Enzyme related genes in modules.csv', header = T)
dt_enz$moduleColors = str_to_title(dt_enz$moduleColors)
clr2 = c('tan3','olivedrab1','turquoise3')

ggplot(data = dt_enz)+
  geom_bar(aes(x=dt_enz$Related.Hormone.Type,group = dt_enz$moduleColors, 
               fill = dt_enz$moduleColors), 
           position = position_dodge2(width = 0.5, preserve = "single")) +
  scale_fill_manual(values = clr2, name = "Module") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10.5), 
        legend.position = "right") +
  labs(title = "Enzyme related genes in modules") +
  ylab("Gene count") +
  xlab("Enzyme")

#============================================================
#All-in-one plot
rm(list = ls(all.names = TRUE))

library(ggplot2)
library(tidyverse)

setwd("D:/UNI/4TH YEAR/RESEARCH/Codes/Arabidopsis")
options(stringsAsFactors = FALSE)
load(file = "networkConstruction.RData")
modules = sort(str_to_title(unique(moduleColors)))
clr = modules[modules != "Grey"]
count = clr
path = "D:/UNI/4TH YEAR/RESEARCH/Codes/Arabidopsis/Plots/figures"

dt = read.csv('Hormone related genes in modules.csv', header = T)
dt$moduleColors = str_to_title(dt$moduleColors)

#Sorting out genes based on hormone/enzyme
title = unique(dt$Related.Hormone.Type)
mod = numeric()

for (i in 1:length(title)) {
  dth = dt[dt$Related.Hormone.Type == title[i],]
  
  
  for (j in 1:length(clr)) {
    mod[j] = sum(dth$moduleColors == clr[j])
  }
  
  count = data.frame(cbind(count ,mod))
}
colnames(count) = c("Module",title)
write.csv(count,'Hormones.csv',quote = F, row.names = F)


#=========================================================================
#Heatmap of hormone related genes

setwd("D:/UNI/4TH YEAR/RESEARCH/Codes/Arabidopsis")
dt = read.csv('Hormone-module.csv', header = T)
dt[is.na(dt)] = 0
dt = unique(dt)
rownames(dt) = dt$gene
colnames(dt) = stringr::str_to_title(colnames(dt))
dt = dt[,-1]
dt = as.matrix(dt)
dt1 = dt[1:42,]
dt2 = dt[43:79,]

library(pheatmap)
library( RColorBrewer)
pheatmap(t(dt1), treeheight_row = 0, treeheight_col = 0,cellwidth = 15, 
         cellheight = 18, color = rev(brewer.pal(8,"PuBu")), scale = "none",
         cluster_rows = F, cluster_cols = F, angle_col = 45)
pheatmap(t(dt2), treeheight_row = 0, treeheight_col = 0,cellwidth = 15, 
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
pheatmap(dt, treeheight_row = 0, treeheight_col = 0,cellwidth = 15, 
         cellheight = 18, color = rev(brewer.pal(8,"PuBu")), scale = "none",
         cluster_rows = F, cluster_cols = F, angle_col = 45, legend = F)

#===========================================================================

rm(list = ls(all.names = TRUE))

setwd("D:/UNI/4TH YEAR/RESEARCH/Codes/Arabidopsis/Export/Modular genes")
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

setwd("D:/UNI/4TH YEAR/RESEARCH/Codes/Arabidopsis/Export/Modular genes")
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

setwd("D:/UNI/4TH YEAR/RESEARCH/Codes/Arabidopsis/Export/Modular genes")
dt = read.csv('Key genes summary.csv')
annot = read.csv('annotation.csv')
key = read.csv('Keygenes-modules.csv')

mt = match(key$probes,annot$Probe.Set.ID)
kg = annot[mt,]
write.csv(kg,'key gene - annot.csv',row.names = F,quote = F)

