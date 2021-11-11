### PCA comparing response - non response for each drug
setwd('/Users/asurace/Documents/GC-PEAC/Desktop/R4RA')

## unsupervised script

library(M3C)
library(edgeR)
library(limma)
library(ComplexHeatmap)
library(dplyr)
library(car)
load("./R4RA_V2_data_trimmed.RData")

### parameters

# variability and expression strength filtering
filtering <- TRUE # set to TRUE for unsupervised clustering
cvx <- 0.075

### main script

## get the data out
data <- txiready$vst
colnames(data) <- gsub('\\.','-',colnames(data))
colnames(data) <- substring(colnames(data), 2)

## remove outliers
  outliers <- read.csv('./outliers.csv',sep=',')
  outliers$ID <- as.character(outliers$ID)
  data <- data[ , -which(names(data) %in% outliers$ID)] 

# get the metadata
meta <- read.csv('./meta_data.csv')
meta <- subset(meta, meta$Visit == 3)
completeFun <- function(data, desiredCols) {
  completeVec <- complete.cases(data[, desiredCols])
  return(data[completeVec, ])
}

meta$Seq_ID.V2 <- as.character(meta$Seq_ID.V2)
meta$Seq_ID.V2 <- gsub('\\.','-',meta$Seq_ID.V2)
meta$Seq_ID.V2 <- paste(meta$Seq_ID.V2,'-Baseline',sep='')


setdiff(colnames(data),meta$Seq_ID.V2)

meta <- subset(meta, meta$Seq_ID.V2 %in% colnames(data))
data <- data[,meta$Seq_ID.V2]

## coefficient of variation
## second order coefficient of variation
mfilter <- 'cv'
if (filtering){
  # filter for just the strongest expressed genes
  tokeep <- row.names(txiraw$txistd$counts[ rowSums((cpm(txiraw$txistd$counts))>1)>=10, ])
  data2 <- subset(data, row.names(data) %in% tokeep)
  # variable gene filter
  if (mfilter == 'cv'){
    # co efficient of variation
    CV <- apply(data2,1,sd)/(rowMeans(data2))
    names <- names(CV)[CV>cvx] 
  }else if (mfilter == 'a2'){
    CV <- apply(data2,1,sd)/(rowMeans(data2))
    CV2 <- CV^2
    A2 <- sqrt(CV2/(CV2+1))
    names <- names(A2)[A2>0.8]
  }else if (mfilter == 'var'){
    xxx <- featurefilter(data2,method='var',percentile=40)
    names <- row.names(xxx$filtered_data)
  }
  data3 <- subset(data2, row.names(data2) %in% names)
  nrow(data3)
}else{
  data3 <- data
}

colnames(meta)[colnames(meta)=='Seq_ID.V2'] <- 'ID'

### 

datax2 <- t(scale(t(data3)))

## 3D PCA
pca <- prcomp(t(datax2))
pca <- as.data.frame(pca$x[,c(1,2,3)])
saveRDS(pca, file='PCscores_R4RA_forDESeq2.rds')
# plot
meta$CDAI.response.status.V7 <- as.factor(meta$CDAI.response.status.V7)
meta$CDAI.response.status.V7 <- ifelse(meta$CDAI.response.status.V7 == 'Non.Responder', 'Non Responder', paste(meta$CDAI.response.status.V7))
meta$Randomised.medication <- as.factor(meta$Randomised.medication)
meta$col <- ifelse(meta$Randomised.medication == 'Rituximab', 'Rituximab Responder', 'Tocilizumab Responder')
meta$col <- ifelse(meta$CDAI.response.status.V7 == 'Non.Responder', 'Non Responder', meta$col )
meta_r <- meta %>% filter(Randomised.medication =='Rituximab')
meta_t <- meta %>% filter(Randomised.medication =='Tocilizumab')

pca_r <- subset(pca, row.names(pca) %in% meta_r$ID)
pca_t <- subset(pca, row.names(pca) %in% meta_t$ID)

library(ggplot2)

p1 <- ggplot(pca, aes(x=PC1, y= PC3, colour = factor(meta$Pathotype.V2))) +
  geom_point(size=4)+
  theme_classic()+
  scale_color_manual(values=c('skyblue1','firebrick1','plum1','grey'))+
  stat_ellipse(level = 0.8, aes(linetype = meta$Pathotype.V2))+
  scale_linetype_manual(values=c(1,1,0,0))+
  labs(colour = 'Pathotype')+
  guides(linetype = 'none')+
  theme(axis.title = element_text(size = 14), axis.text = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14))
p1


## RTX response PCA
p2 <-  ggplot(pca_r, aes(x=PC1, y= PC2, colour = factor(meta_r$CDAI.response.status.V7))) +
  geom_point(size=4)+
  scale_color_manual(values=c('red', 'blue'))+
  labs(colour = 'Response to RTX')+
  stat_ellipse(level=0.8)+
  theme_classic()+
  theme(axis.title = element_text(size = 15), axis.text = element_text(size = 15),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 15))

p2


p3 <-  ggplot(pca_t, aes(x=PC1, y= PC2, colour = factor(meta_t$CDAI.response.status.V7))) +
  geom_point(size=4)+
  scale_color_manual(values=c('red', 'gold3'))+
  labs(colour = 'Response to TOC')+
  stat_ellipse(level=0.8)+
  theme_classic()+
  theme(axis.title = element_text(size = 15), axis.text = element_text(size = 15),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 15))


p3

