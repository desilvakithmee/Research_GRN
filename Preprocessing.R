#=================================== PRE-PROCESSING ====================================================

setwd("D:/UNI/4TH YEAR/RESEARCH/Codes/Arabidopsis")

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

##=========================== DIFFERENTIAL GENE EXPRESSION =============================================

##DEG using fold change

setwd("D:/UNI/4TH YEAR/RESEARCH/Codes/Arabidopsis")
exp = read.csv('data.csv',header = T)
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
write.csv(data,'fc_data.csv',row.names = F)

#===================================================================================================

## p-value

setwd("D:/UNI/4TH YEAR/RESEARCH/Codes/Arabidopsis")
exp = read.csv('data.csv')

t.test(data[1,2:3], data[1,4:5])
result = numeric()

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

#=================================================================================================

exp = read.csv('fc_data.csv',header = T)

#filtering data based on fold change
flt = read.delim('fold_change.txt',header = F)
genes = as.character(flt$V1)
key = match(genes,exp$exp.ID)    #sum((key == 'NA') == TRUE)
data = exp[key,]
write.csv(data, 'fc_data_filter.csv', row.names = F)

#=================================================================================================

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
cat('Only',count1,'out of 483 key genes',count2,'out of 42 SE markers and',count3,'out of 12 markers are present.')


#=================================== FILTERED INPUT DATASET =========================================

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
key = match(genes,exp$exp.ID)    #sum((key == 'NA') == TRUE)
data = exp[key,]

#Conversion of probe IDs to gene IDs
ID = read.delim('ID.txt',header = F)
key = match(data$exp.ID,ID$V1)
header = casefold(ID$V3[key],upper = T)
data = data[,-1]  #removing gene IDs
data = as.data.frame(cbind(header,data))
write.csv(data,'data_processed_avg.csv',row.names = F)


#====================================== VISUALIZATION ==============================================

#===================================================================================================
#Fold change

setwd("D:/UNI/4TH YEAR/RESEARCH/Codes/Arabidopsis")
data = read.csv('fc_data.csv')

#Fold changes
ggplot(data, aes(x=day7_log2)) +
  geom_histogram(binwidth=.5, colour="black", fill="white") +
  ggtitle(label =  "Fold Change - Timepoint 2") +
  xlab(label =  "log2 Fold change") +
  ylab(label = "Frequency")
ggsave("FC_day7.jpg", path = "D:/UNI/4TH YEAR/RESEARCH/Codes/Arabidopsis/Plots")

ggplot(data, aes(x=day14_log2)) +
  geom_histogram(binwidth=.5, colour="black", fill="white") +
  ggtitle(label =  "Fold Change - Timepoint 3") +
  xlab(label =  "log2 Fold change") +
  ylab(label = "Frequency")
ggsave("FC_day14.jpg", path = "D:/UNI/4TH YEAR/RESEARCH/Codes/Arabidopsis/Plots")

ggplot(data, aes(x=som_log2)) +
  geom_histogram(binwidth=.5, colour="black", fill="white") +
  ggtitle(label =  "Fold Change - Timepoint 3") +
  xlab(label =  "log2 Fold change") +
  ylab(label = "Frequency")
ggsave("FC_som.jpg", path = "D:/UNI/4TH YEAR/RESEARCH/Codes/Arabidopsis/Plots")


#=====================================================================================================

#Volcano Plot

library(ggplot2)
setwd("D:/UNI/4TH YEAR/RESEARCH/Codes/Arabidopsis")
data = read.csv('fc_data_filter.csv')

#day7
results = cbind(log = data$day7_log2, pval = data$pvalue_day7)
results = as.data.frame(results)
results$probename = data$exp.ID

ggplot(data = results, aes(x = log, y = -1*log2(pval)))  + 
  geom_point(colour = "skyblue4") +
  ggtitle(label =  "Volcano Plot - Timepoint 2") +
  xlab(label =  "log2 Fold change") +
  ylab(label = "log10 p-value") +
  coord_cartesian(ylim = c(-2,10))
ggsave("Volcano_day7.jpg", path = "D:/UNI/4TH YEAR/RESEARCH/Codes/Arabidopsis/Plots")

#day14
results = cbind(log = data$day14_log2, pval = data$pvalue_day14)
results = as.data.frame(results)
results$probename = data$exp.ID

ggplot(data = results, aes(x = log, y = -1*log10(pval)))  + 
  geom_point(colour = "skyblue4") +
  ggtitle(label =  "Volcano Plot - Timepoint 3") +
  xlab(label =  "log2 Fold change") +
  ylab(label = "log10 p-value") +
  coord_cartesian(ylim = c(-2,5))
ggsave("Volcano_day14.jpg", path = "D:/UNI/4TH YEAR/RESEARCH/Codes/Arabidopsis/Plots")

#som
results = cbind(log = data$som_log2, pval = data$pvalue_som)
results = as.data.frame(results)
results$probename = data$exp.ID

ggplot(data = results, aes(x = log, y = -1*log10(pval)))  + 
  geom_point(colour = "skyblue4") +
  ggtitle(label =  "Volcano Plot - Timepoint 4") +
  xlab(label =  "log2 Fold change") +
  ylab(label = "log10 p-value") +
  coord_cartesian(ylim = c(-2,5))
ggsave("Volcano_som.jpg", path = "D:/UNI/4TH YEAR/RESEARCH/Codes/Arabidopsis/Plots")

