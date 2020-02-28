#========================= PRE-PROCESSING ==================================

setwd("D:/UNI/4TH YEAR/RESEARCH/Codes/Zygotic_Embryogenesis")

#ADDING TAIR IDS TO THE SE MARKER LIST

ID = read.delim('ID.txt',header = F)
markers = read.csv('SE_markers.csv')
mrk = as.character(markers$Gene_ID)
mrk = trimws(mrk)
genes = casefold(ID$V3,upper = T)
key = match(mrk,genes)
tair_ID = ID[key,1]
markers = data.frame(markers,tair_ID)
write.csv(markers,'SE_markers.csv',row.names = F)
head(markers)

##================= DIFFERENTIAL GENE EXPRESSION ===========================
rm(list = ls(all.names = TRUE))

##DEG using fold change

setwd("D:/UNI/4TH YEAR/RESEARCH/Codes/Zygotic_Embryogenesis")
options(stringsAsFactors = FALSE)
exp = read.delim('GSE123010_Normalized_counts.txt', header = T)
exp = na.omit(exp)

zygote = rowMeans(exp[,2:3])
octant = rowMeans(exp[,4:5])
globular = rowMeans(exp[,6:7])
heart = rowMeans(exp[,8:9])
torpedo = rowMeans(exp[,10:11])
bent  = rowMeans(exp[,12:13])
mature  = rowMeans(exp[,14:15])
data = data.frame(exp$ID,zygote,octant,globular,heart,torpedo,bent,mature)

fc_oct = octant/zygote
fc_glb = globular/zygote
fc_hrt = heart/zygote
fc_trp = torpedo/zygote
fc_bnt = bent/zygote
fc_mtr = mature/zygote

data = data.frame(data,fc_oct,fc_glb,fc_hrt,fc_trp,fc_bnt,fc_mtr)
log = log2(data[,9:14])
names(log) = c('oct_log2','glb_log2','hrt_log2','trp_log2',
               'bnt_log2','mtr_log2')
data = data.frame(data,log)

data[mapply(is.infinite, data)] = NA
data[is.na(data)] = 0

thresh = 6
deg_oct = data[abs(data$fc_oct)>thresh,1]
deg_glb = data[abs(data$fc_glb)>thresh,1]
deg_hrt = data[abs(data$fc_hrt)>thresh,1]
deg_trp = data[abs(data$fc_trp)>thresh,1]
deg_bnt = data[abs(data$fc_bnt)>thresh,1]
deg_mtr = data[abs(data$fc_mtr)>thresh,1]

id = c(as.character(deg_bnt),as.character(deg_glb),as.character(deg_hrt),
       as.character(deg_mtr),as.character(deg_oct),as.character(deg_trp))
id = unique(id)

write.table(id,'DEGs.txt',row.names = F, col.names = F, quote = F)

mrk = as.character(data$exp.ID)
key = match(id,mrk)
dt = data[key,]
write.csv(dt,'data_processed.csv',row.names = F, quote = F)

#===========================================================================
rm(list = ls(all.names = TRUE))
## p-value

setwd("D:/UNI/4TH YEAR/RESEARCH/Codes/Zygotic_Embryogenesis")
options(stringsAsFactors = FALSE)
data = read.csv('data_filtered.csv', header = T)
result = numeric()

#function to calculate the p-value
ttest <- function(df, grp1, grp2) {
  for (i in 1:dim(df)[1]) {
    x = df[i,grp1]
    y = df[i,grp2]
    rslt = t.test(x,y)
    result[i] = rslt$p.value
  }
  return(result)
}

pvalue_oct = ttest(df = data, grp1 = c(2:3), grp2 = c(4:5))
pvalue_glb = ttest(df = data, grp1 = c(2:3), grp2 = c(6:7))
pvalue_hrt = ttest(df = data, grp1 = c(2:3), grp2 = c(8:9))
pvalue_trp = ttest(df = data, grp1 = c(2:3), grp2 = c(10:11))
pvalue_bnt = ttest(df = data, grp1 = c(2:3), grp2 = c(12:13))
pvalue_mtr = ttest(df = data, grp1 = c(2:3), grp2 = c(14:15))

exp = read.csv('data_processed.csv')
key = match(data$ID, exp$exp.ID)
exp = exp[key,]
exp = data.frame(exp,pvalue_oct,pvalue_glb,pvalue_hrt,pvalue_trp,pvalue_bnt,
                  pvalue_mtr)
write.csv(exp,'data_processed.csv', row.names = F)

#===========================================================================
#filtering data based on p-value

rm(list = ls(all.names = TRUE))

data = read.csv('data_processed.csv', header = T)

thresh = 0.05
deg_oct = data[abs(data$pvalue_oct)<thresh,1]
deg_glb = data[abs(data$pvalue_glb)>thresh,1]
deg_hrt = data[abs(data$pvalue_hrt)>thresh,1]
deg_trp = data[abs(data$pvalue_trp)>thresh,1]
deg_bnt = data[abs(data$pvalue_bnt)>thresh,1]
deg_mtr = data[abs(data$pvalue_mtr)>thresh,1]

id = c(as.character(deg_bnt),as.character(deg_glb),as.character(deg_hrt),
       as.character(deg_mtr),as.character(deg_oct),as.character(deg_trp))
id = unique(id)

write.table(id,'DEGs.txt',row.names = F, col.names = F, quote = F)

#FILTERED INPUT DATASET

exp = read.delim('GSE123010_Normalized_counts.txt', header = T)
exp = na.omit(exp)
mrk = as.character(exp$ID)
key = match(id,mrk)
dt = exp[key,]
write.csv(dt,'data_filtered.csv',row.names = F, quote = F)

#===========================================================================
rm(list = ls(all.names = TRUE))

#Checking if key genes and SE markers are present in the filtered dataset

genes = read.delim('gene_list_meinke.txt')
key_genes = as.character(genes$From)
key_genes = toupper(key_genes)

markers = read.delim('markers.txt',sep = "")
mrk_genes = as.character(markers$mrk)

probes = read.delim('DEGs.txt', header = F)
probes = as.character(probes$V1)

tf = read.csv('TF list.csv')
tfs = tf$Locus

mt1 = match(key_genes,probes)
mt2 = match(tfs,probes)
mt3 = match(mrk_genes,probes)
count1 = sum(is.na(mt1) == F)
count2 = sum(is.na(mt2) == F)
count3 = sum(is.na(mt3) == F)
cat('Only',count1,'out of 483 key genes ',count2,' transcription factors and',
    count3,'out of 12 markers are present.')


#========================== VISUALIZATION ==================================

#Fold change
library('ggplot2')
setwd("D:/UNI/4TH YEAR/RESEARCH/Codes/Zygotic_Embryogenesis")
data = read.csv('fc_data.csv')
x = integer()
y = integer()
for (i in 15:20) {
  x = rbind(x,sum(data[,i]>(6)))
  y = c(y,sum(data[,i]<(-6)))
}

#Fold changes
q1 = ggplot(data, aes(x=day7_log2)) +
  geom_histogram(binwidth=.5, colour="black", fill="white") +
  ggtitle(label =  "Fold Change - Timepoint 2") +
  xlab(label =  "log2 Fold change") +
  ylab(label = "Frequency")
ggsave("FC_day7.jpg",plot = q1, path = "D:/UNI/4TH YEAR/RESEARCH/Codes/
       Zygotic_Embryogenesis/Plots")

q2 = ggplot(data, aes(x=day14_log2)) +
  geom_histogram(binwidth=.5, colour="black", fill="white") +
  ggtitle(label =  "Fold Change - Timepoint 3") +
  xlab(label =  "log2 Fold change") +
  ylab(label = "Frequency")
ggsave("FC_day14.jpg",plot = q2, path = "D:/UNI/4TH YEAR/RESEARCH/Codes/
       Zygotic_Embryogenesis/Plots")

q3 = ggplot(data, aes(x=som_log2)) +
  geom_histogram(binwidth=.5, colour="black", fill="white") +
  ggtitle(label =  "Fold Change - Timepoint 3") +
  xlab(label =  "log2 Fold change") +
  ylab(label = "Frequency")
ggsave("FC_som.jpg", plot = q3, path = "D:/UNI/4TH YEAR/RESEARCH/Codes/
       Zygotic_Embryogenesis/Plots")


#===========================================================================

#Volcano Plot

library(ggplot2)
library(grDevices)
setwd("D:/UNI/4TH YEAR/RESEARCH/Codes/Zygotic_Embryogenesis")
data = read.csv('data_processed.csv')

#day7
results = cbind(log = data$trp_log2, pval = data$pvalue_trp)
results = as.data.frame(results)

group = rep("No change",dim(results)[1])
for (i in 1:dim(results)[1]) {
  if (results$log[i] > 6){
    group[i] = "Upregulated"
  }
  if (results$log[i] < -6){
    group[i] = "Downregulated"
  }
  }
results = cbind(results,group)

p1 = ggplot(data = results, aes(x = results$log, y = -1*log2(results$pval)))  + 
  geom_point(aes(colour = group, shape = group)) + 
  ggtitle(label =  "Volcano Plot - Timepoint 2") +
  xlab(label =  "log2 Fold change") +
  ylab(label = "log10 p-value") +
  coord_cartesian(ylim = c(-2,17))
ggsave("Volcano_StageII.jpg", plot = p1, path = "D:/UNI/4TH YEAR/RESEARCH/Codes/
       Zygotic_Embryogenesis/Plots")

results = cbind(ID = data$tair_ID,results)
results = results[order(results$group),]
write.csv(results,'Plots/Day7 chart.csv',row.names = F)

#day14
results = cbind(log = data$day14_log2, pval = data$pvalue_day14)
results = as.data.frame(results)

group = rep("No change",dim(results)[1])
for (i in 1:dim(results)[1]) {
  if (results$log[i] > 2){
    group[i] = "Upregulated"
  }
  if (results$log[i] < -2){
    group[i] = "Downregulated"
  }
}
results = cbind(results,group)

p2 = ggplot(data = results, aes(x = results$log, y = -1*log10(results$pval)))  + 
  geom_point(aes(colour = group, shape = group)) +
  ggtitle(label =  "Volcano Plot - Timepoint 3") +
  xlab(label =  "log2 Fold change") +
  ylab(label = "log10 p-value") +
  coord_cartesian(ylim = c(-2,5))
ggsave("Volcano_day14.jpg", plot = p2, path = "D:/UNI/4TH YEAR/RESEARCH/Codes/
       Zygotic_Embryogenesis/Plots")

results = cbind(ID = data$tair_ID,results)
results = results[order(results$group),]
write.csv(results,'Plots/Day14 chart.csv',row.names = F)

#som
results = cbind(log = data$som_log2, pval = data$pvalue_som)
results = as.data.frame(results)

group = rep("No change",dim(results)[1])
for (i in 1:dim(results)[1]) {
  if (results$log[i] > 2){
    group[i] = "Upregulated"
  }
  if (results$log[i] < -2){
    group[i] = "Downregulated"
  }
}
results = cbind(results,group)

p3 = ggplot(data = results, aes(x = results$log, y = -1*log10(results$pval)))  + 
  geom_point(aes(colour = group, shape = group)) +
  ggtitle(label =  "Volcano Plot - Timepoint 4") +
  xlab(label =  "log2 Fold change") +
  ylab(label = "log10 p-value") +
  coord_cartesian(ylim = c(-2,5))
ggsave("Volcano_som.jpg", plot = p3, path = "D:/UNI/4TH YEAR/RESEARCH/Codes/
       Zygotic_Embryogenesis/Plots")

results = cbind(ID = data$tair_ID,results)
results = results[order(results$group),]
write.csv(results,'Plots/Som chart.csv',row.names = F)

#===============================================================================
setwd( "D:/UNI/4TH YEAR/RESEARCH/Codes/Zygotic_Embryogenesis")
dt = read.csv('fc_data.csv')

up_stII = dt[dt$log_octant>6,1]
up_stIII = dt[dt$log_globular>6,1]
up_stIV = dt[dt$log_heart>6,1]
up_stV = dt[dt$log_torpedo>6,1]
up_stVI = dt[dt$log_bent>6,1]
up_stVII = dt[dt$log_mature>6,1]

down_stII = dt[dt$log_octant<(-6),1]
down_stIII = dt[dt$log_globular<(-6),1]
down_stIV = dt[dt$log_heart<(-6),1]
down_stV = dt[dt$log_torpedo<(-6),1]
down_stVI = dt[dt$log_bent<(-6),1]
down_stVII = dt[dt$log_mature<(-6),1]

setwd( "D:/UNI/4TH YEAR/RESEARCH/Codes/Zygotic_Embryogenesis/DEG")
write.table(up_stII,file = "up_stII.txt",quote = F,row.names = F, col.names = F)
write.table(up_stIII,file = "up_stIII.txt",quote = F,row.names = F, col.names = F)
write.table(up_stIV,file = "up_stIV.txt",quote = F,row.names = F, col.names = F)
write.table(up_stV,file = "up_stV.txt",quote = F,row.names = F, col.names = F)
write.table(up_stVI,file = "up_stVI.txt",quote = F,row.names = F, col.names = F)
write.table(up_stVII,file = "up_stVII.txt",quote = F,row.names = F, col.names = F)

write.table(down_stII,file = "down_stII.txt",quote = F,row.names = F, col.names = F)
write.table(down_stIII,file = "down_stIII.txt",quote = F,row.names = F, col.names = F)
write.table(down_stIV,file = "down_stIV.txt",quote = F,row.names = F, col.names = F)
write.table(down_stV,file = "down_stV.txt",quote = F,row.names = F, col.names = F)
write.table(down_stVI,file = "down_stVI.txt",quote = F,row.names = F, col.names = F)
write.table(down_stVII,file = "down_stVII.txt",quote = F,row.names = F, col.names = F)

#==============================================
#TFs
setwd("D:/UNI/4TH YEAR/RESEARCH/Codes/Zygotic_Embryogenesis")
options(stringsAsFactors = F)

tf_list = read.csv('TF list.csv')
load('networkConstruction.RData')
load('dataInput.RData')

dt = data.frame(gene = colnames(data),moduleColors)
key = match(tf_list$Locus,dt$gene)
tf_gene = na.omit(dt[key,])
tf_gene = tf_gene[order(tf_gene$moduleColors),]
clr = unique(tf_gene$moduleColors)

sum = numeric()
for (i in 1:length(clr)) {
  cl = clr[i]
  kk = tf_gene[tf_gene$moduleColors == cl,]
  x = dim(kk)[1]
  sum = c(sum,x)
}

ex = data.frame(module = clr, count = sum)
write.csv(ex,'TF modules.csv', row.names = F, quote = F)


