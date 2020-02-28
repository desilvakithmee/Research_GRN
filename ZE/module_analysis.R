
#======================================= Module Analysis ==============================================
rm(list = ls(all.names = TRUE))

setwd("D:/UNI/4TH YEAR/RESEARCH/Codes/Zygotic_Embryogenesis/Export/Blue")
library(WGCNA)

dt = read.delim('CytoscapeInput-edges-blue.txt')
genes = unique(c(as.character(dt$fromNode,dt$toNode)))

exp = read.csv('data_filtered.csv')
key = match(genes,exp$ID)
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
Blue = c(as.character(bgenes$row,bgenes$column))
Blue_genes = unique(Blue)

#Saving for further clustering
key = match(Blue_genes,exp$ID)
data = exp[key,]
data = data[,-1]
data = t(data)
colnames(data) = Blue_genes
data = as.data.frame(data)
save(data,file = 'Blue_data.RData')

#If there is no further clustering
diff = setdiff(genes,Blue_genes)
dt_new = dt
for (i in 1:length(diff)) {
  df = diff[i]
  dt_new = dt_new[!(dt$fromNode==df | dt$toNode==df),]
}
write.table(dt_new,'Blue_edges.txt',row.names = F,sep = "\t", quote = F)


#==================================================================================================

#Checking for key genes

markers = read.delim('Key genes.txt',header = F)
key = match(Blue_genes,markers)
mrk = markers[key]
cat('Key genes in Blue module:',mrk[!is.na(mrk)])
write.table(Blue_genes,'Blue genes.txt',row.names = F,col.names = F, sep = "\t", quote = F)

#==================================================================================================
#Transcription factors

TF = read.csv('TF list.csv',header = T)
TF = TF$Locus
nodes = intersect(Blue_genes,TF)

write.table(nodes,'Blue_TFs.txt',row.names = F,col.names = F, sep = "\t", quote = F)

#==================================================================================================
##ANALYSIS OF THE Blue MODULE USING WGCNA is not needed

#LEA Genes
setwd("D:/UNI/4TH YEAR/RESEARCH/Codes/Zygotic_Embryogenesis/Export/Blue")
lea = read.csv('LEA.csv')
lea = lea[,1:6]
lea_genes = lea$ï..AGI.code
lea_genes = toupper(lea_genes)
blue = read.delim('Blue genes.txt', header = F)

key = match(blue$V1,lea_genes)
key = na.omit(key)
dt = lea[key,]
write.csv(dt,'LEA proteins in Mature embryo.csv',row.names = F, quote = F)