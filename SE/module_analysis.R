
#========================= Module Analysis ===================================
rm(list = ls(all.names = TRUE))

setwd("D:/UNI/4TH YEAR/RESEARCH/Codes/Arabidopsis/Export/Black")
library(WGCNA)

dt = read.delim('CytoscapeInput-edges-black.txt')
genes = unique(c(as.character(dt$fromNode,dt$toNode)))

exp = read.csv('data_processed.csv')
key = match(genes,exp$header)
data = exp[key,]
data = data[,-1]
data = t(data)
colnames(data) = genes
data = as.data.frame(data)

cor = corAndPvalue(data, alternative = 'two.sided')
sum(cor$p<0.05)

flattenCorrMatrix = function(cormat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = colnames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut]
  )
}

p_val = flattenCorrMatrix(cor$p)  
bgenes = p_val[p_val$cor<0.05,]  #significant correlation
Black = c(as.character(bgenes$row,bgenes$column))
Black_genes = unique(Black)
write.table(Black_genes,'Black genes.txt',row.names = F,col.names = F, sep = "\t",
            quote = F)
#Since it is only 164 genes, not necessary to further breakdown to modules

#=================================================================================

#Checking for key genes

markers = read.delim('Key genes.txt',header = F)
key = match(Black_genes,markers$V1)
mrk = markers[key,]
cat('Key genes in Black module:',mrk[!is.na(mrk)])

#=================================================================================


#PlantGSEA -> Job ID: 595490762
#agriGO -> Job ID: 107664841