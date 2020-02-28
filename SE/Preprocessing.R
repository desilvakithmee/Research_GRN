#========================= PRE-PROCESSING ==================================

setwd("D:/UNI/4TH YEAR/RESEARCH/Codes/Arabidopsis")

#ADDING TAIR IDS TO THE SE MARKER LIST

ID = read.delim('ID.txt',header = F)
head(ID)
markers = read.csv('SE_markers.csv')
head(markers)
mrk = as.character(markers$Gene_ID)
mrk = trimws(mrk)
genes = casefold(ID$V3,upper = T)
key = match(mrk,genes)
tair_ID = ID[key,1]
markers = data.frame(markers,tair_ID)
write.csv(markers,'SE_markers.csv',row.names = F)
head(markers)

##================= DIFFERENTIAL GENE EXPRESSION ===========================

##DEG using fold change

setwd("D:/UNI/4TH YEAR/RESEARCH/Codes/Arabidopsis")
exp = read.csv('data_SE.csv',header = T)
exp = exp[1:22746,] #remove AFFX chips

zyg = rowMeans(exp[,2:3])
day7 = rowMeans(exp[,4:5])
day14 = rowMeans(exp[,6:7])
som = rowMeans(exp[,8:9])
data = data.frame(exp$ID,zyg,day7,day14,som)

fc7 = day7/zyg
fc14 = day14/zyg
fcsom = som/zyg
data = data.frame(data,fc7,fc14,fcsom)
log = log2(data[,6:8])
names(log) = c('day7_log2','day14_log2','som_log2')
data = data.frame(data,log)
deg7 = data[data$day7_log2>2 |data$day7_log2 < -2,1]
deg14 = data[data$day14_log2>2 |data$day14_log2 < -2,1]
deg_som = data[data$som_log2>2 |data$som_log2 < -2,1]
id = c(as.character(deg7),as.character(deg14),as.character(deg_som))
id = unique(id)
write.table(id,'fold_change.txt',row.names = F, col.names = F)

up_day7 = data[data$day7_log2>2,1]
up_day14 = data[data$day14_log2>2,1]
up_zyg = data[data$som_log2>2,1]
down_day7 = data[data$day7_log2<(-2),1]
down_day14 = data[data$day14_log2<(-2),1]
down_zyg = data[data$som_log2<(-2),1]

setwd( "D:/UNI/4TH YEAR/RESEARCH/Codes/Arabidopsis/DEG")
write.table(up_day7,file = "up_day7.txt",quote = F,row.names = F, col.names = F)
write.table(up_day14,file = "up_day14.txt",quote = F,row.names = F, col.names = F)
write.table(up_zyg,file = "up_zyg.txt",quote = F,row.names = F, col.names = F)
write.table(down_day7,file = "down_day7.txt",quote = F,row.names = F, col.names = F)
write.table(down_day14,file = "down_day14.txt",quote = F,row.names = F, col.names = F)
write.table(down_zyg,file = "down_zyg.txt",quote = F,row.names = F, col.names = F)


fi_se = day7/zyg
fi_th = day14/zyg
fi_fo = som/zyg
se_th = day7/day14
se_fo = day7/zyg
th_fo = day14/zyg

data = data.frame(data,fi_se,fi_th,fi_fo,se_th,se_fo,th_fo)
log = log2(data[,6:11])
log = data.frame(data$exp.ID,log)

deg = character()
for (i in 2:dim(log)[2]) {
  dg = as.character(log[abs(log[,i])>2,1])
  deg = c(deg,dg)
}
id = unique(deg)

write.table(id,'fold_change2.txt',row.names = F, col.names = F,quote = F)

ID = read.delim('ID.txt',header = F)
mrk = as.character(data$exp.ID)
genes = ID$V1
key = match(mrk,genes)
tair_ID = toupper(ID[key,3])
data = data[,-1]
data = data.frame(tair_ID,data)
write.csv(data,'fc_data.csv',row.names = F)
head(data)

#===========================================================================

## p-value

setwd("D:/UNI/4TH YEAR/RESEARCH/Codes/Arabidopsis")
exp = read.csv('data.csv')

t.test(data[1,2:3], data[1,4:5])
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

pvalue_day7 = ttest(df = data, grp1 = c(2:3), grp2 = c(4:5))
pvalue_day14 = ttest(df = data, grp1 = c(2:3), grp2 = c(6:7))
pvalue_som = ttest(df = data, grp1 = c(2:3), grp2 = c(8:9))

data = read.csv('fc_data.csv')
data = data.frame(data,pvalue_day7,pvalue_day14,pvalue_som)
write.csv(data, 'fc_data.csv', row.names = F)
head(data)

#===========================================================================

exp = read.csv('fc_data.csv',header = T)

#filtering data based on fold change
flt = read.delim('fold_change.txt',header = F)
genes = as.character(flt$V1)
key = match(genes,exp$exp.ID)    #sum((key == 'NA') == TRUE)
data = exp[key,]
write.csv(data, 'fc_data_filter.csv', row.names = F)

#===========================================================================

#Checking if key genes and SE markers are present in the filtered dataset

genes = read.delim('gene_list_meinke.txt')
key_genes = as.character(genes$To)

markers = read.csv('SE_markers.csv', header = T)
mrk = as.character(markers$tair_ID)

markers = read.delim('markers.txt',sep = "")
mrk_genes = as.character(markers$affy_ID)

probes = read.delim('fold_change.txt', header = F)
probes = as.character(probes$V1)
mt1 = match(key_genes,probes)
mt2 = match(mrk,probes)
mt3 = match(mrk_genes,probes)
count1 = sum(is.na(mt1) == F)
count2 = sum(is.na(mt2) == F)
count3 = sum(is.na(mt3) == F)
cat('Only',count1,'out of 483 key genes',count2,'out of 42 SE markers and',
    count3,'out of 12 markers are present.')


#======================== FILTERED INPUT DATASET ===========================

rm(list = ls(all.names = TRUE))

#PREPROCESSING THE INPUT DATASET

#Reading data
exp = read.csv('data.csv',header = T)
exp = exp[1:22746,] #remove AFFX chips - used for QC

#exp = read.csv('fc_data.csv',header = T)
#exp = exp[,1:5] #for average of replicates

#exploring data
dim(exp)
head(exp)

#filtering data based on fold change
flt = read.delim('fold_change.txt',header = F)
genes = as.character(flt$V1)
key = match(genes,exp$ID)    #sum((key == 'NA') == TRUE)
data = exp[key,]

#Conversion of probe IDs to gene IDs
ID = read.delim('ID.txt',header = F)
key = match(data$ID,ID$V1)
header = casefold(ID$V3[key],upper = T)
data = data[,-1]  #removing gene IDs
data = as.data.frame(cbind(header,data))
write.csv(data,'data_processed_avg.csv',row.names = F)
head(data)

nm = as.character(data$header)
write.table(nm,'DEGs.txt',quote = F,row.names = F,col.names = F)

#===========================================================================
#Combining all marker and key genes for convenience

rm(list = ls(all.names = TRUE))

x = read.delim('gene_list_meinke.txt')
y = read.csv('SE_markers.csv')
z = read.delim('markers.txt',sep = "")

x = toupper(x$From)
y = y$Gene_ID
z = z$mrk
dt = c(x,y,z)
dt = unique(dt)
write.table(dt,'Key genes.txt',row.names = F,col.names = F)

#========================== VISUALIZATION ==================================

#Fold change
library('ggplot2')
setwd("D:/UNI/4TH YEAR/RESEARCH/Codes/Arabidopsis")
data = read.csv('fc_data.csv')

#Fold changes
q1 = ggplot(data, aes(x=day7_log2)) +
  geom_histogram(binwidth=.5, colour="black", fill="white") +
  ggtitle(label =  "Fold Change - Timepoint 2") +
  xlab(label =  "log2 Fold change") +
  ylab(label = "Frequency")
ggsave("FC_day7.jpg",plot = q1, path = "D:/UNI/4TH YEAR/RESEARCH/Codes/
       Arabidopsis/Plots")

q2 = ggplot(data, aes(x=day14_log2)) +
  geom_histogram(binwidth=.5, colour="black", fill="white") +
  ggtitle(label =  "Fold Change - Timepoint 3") +
  xlab(label =  "log2 Fold change") +
  ylab(label = "Frequency")
ggsave("FC_day14.jpg",plot = q2, path = "D:/UNI/4TH YEAR/RESEARCH/Codes/
       Arabidopsis/Plots")

q3 = ggplot(data, aes(x=som_log2)) +
  geom_histogram(binwidth=.5, colour="black", fill="white") +
  ggtitle(label =  "Fold Change - Timepoint 3") +
  xlab(label =  "log2 Fold change") +
  ylab(label = "Frequency")
ggsave("FC_som.jpg", plot = q3, path = "D:/UNI/4TH YEAR/RESEARCH/Codes/
       Arabidopsis/Plots")


#===========================================================================

#Volcano Plot

library(ggplot2)
setwd("D:/UNI/4TH YEAR/RESEARCH/Codes/Arabidopsis")
data = read.csv('fc_data.csv')

#day7
results = cbind(log = data$day7_log2, pval = data$pvalue_day7)
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

p1 = ggplot(data = results, aes(x = results$log, y = -1*log2(results$pval)))  + 
  geom_point(aes(colour = group, shape = group)) + 
  ggtitle(label =  "Volcano Plot - Timepoint 2") + 
  xlab(label =  "log2 Fold change") +
  ylab(label = "log10 p-value") +
  coord_cartesian(ylim = c(-2,10))
ggsave("Volcano_day7.jpg", plot = p1, path = "D:/UNI/4TH YEAR/RESEARCH/Codes/
       Arabidopsis/Plots")

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
       Arabidopsis/Plots")

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
       Arabidopsis/Plots")

results = cbind(ID = data$tair_ID,results)
results = results[order(results$group),]
write.csv(results,'Plots/Som chart.csv',row.names = F)

#==============================================
#TFs
setwd("D:/UNI/4TH YEAR/RESEARCH/Codes/Arabidopsis")
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
#========================================================================
setwd("D:/UNI/4TH YEAR/RESEARCH/Codes/Arabidopsis")

markers = read.delim('Key genes.txt', header = F)
tf_list = read.csv('TF list.csv')
mrk = as.character(markers$V1)
tf = as.character(tf_list$Locus)

#Original
dt = read.csv('Plots/Day7 chart.csv')
probes = as.character(dt$ID)
mt1 = match(mrk,probes)
mt2 = match(tf,probes)
count1 = sum(is.na(mt1) == F)
count2 = sum(is.na(mt2) == F)

cat(count1,'key genes and',count2,'transcription factors are present.')

#Filtered dataset
dt = read.csv('data_processed.csv')
probes = as.character(dt$header)
mt1 = match(mrk,probes)
mt2 = match(tf,probes)
count1 = sum(is.na(mt1) == F)
count2 = sum(is.na(mt2) == F)

cat(count1,'key genes and',count2,'transcription factors are present.')


dt = read.csv('fc_data.csv')
probes = as.character(dt$tair_ID)
mt1 = match(mrk,probes)
xx = dt[mt1,]
xx = na.omit(xx)
