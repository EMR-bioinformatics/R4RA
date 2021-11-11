### DESeq2 of clusters with PC1 as covariate 
library(M3C)
library(edgeR)
library(limma)
load("./R4RA_V2_data_trimmed.RData")

### parameters

# variability and expression strength filtering
filtering <- TRUE # set to TRUE for unsupervised clustering
cvx <- 0.075

data <- txiready$vst

## remove outliers
outliers <- read.csv('./outliers.csv',sep=',')
outliers$ID <- as.character(outliers$ID)
data <- data[ , -which(names(data) %in% outliers$ID)] # 12 initial outliers

#
data <- fix_colnames(data)

# get the metadata
meta <- read.csv('./meta_data_V14_imp.csv')
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

### 

data3 <- t(scale(t(data3)))

### clustering


set.seed(123)
res <- M3C(data3,method=2,clusteralg='pam',repsref=1000,repsreal=1000) 

# reorder data according to any clustering algorithm
ind <- sort(as.vector(res$assignments),index.return=TRUE)
datax <- data3[,ind$ix]
annonx <- meta[ind$ix,]
annonx$consensuscluster <- as.factor(as.character(ind$x))

### optional reordering of consensus clusters manually
reorder <- TRUE
if (reorder){
  ## reorder consensus clusters manually
  newv <- as.factor(c(2,1)) 
  annonx <- annonx[order(match(annonx$consensuscluster, newv)),] # reorder clinical
  annonx$consensuscluster <- as.character(1:length(newv))[ match(annonx$consensuscluster, as.character(unique(annonx$consensuscluster)) ) ]
  #colnames(datax) <- gsub('\\.','-',colnames(datax)) # reformatting
  datax <- datax[,annonx$Seq_ID.V2] # get in right order -- change for PEAC vs R4RA
}

datax2 <- datax # we manipulate this later so save a copy

library(ggplot2)
library(DESeq2)
library(qusage)
library(stringr)
library(qvalue)
library(dplyr)

txi <- readRDS("txi.RDS")

##filter low expressed genes
tokeep <- row.names(txi$counts[ rowSums((cpm(txi$counts))>1)>=10, ])

annonx$consensuscluster <- as.factor(annonx$consensuscluster)
rownames(annonx) <- annonx$Seq_ID.V2
coldata <- annonx # need coldata

## add PC1 as covariate
coldata <-right_join(coldata, PC1, by='Seq_ID.V2')

t_meta <- coldata %>% filter (Randomised.medication == 'Tocilizumab')
r_meta <- coldata %>% filter (Randomised.medication == 'Rituximab')

r_txi <- txi 
t_txi <- txi

##filter off low expressed genes

r_txi$abundance <- r_txi$abundance[which(rownames(r_txi$abundance) %in% tokeep), r_meta$Seq_ID.V2]
r_txi$counts    <- r_txi$counts[which(rownames(r_txi$counts) %in% tokeep), r_meta$Seq_ID.V2]
r_txi$length    <- r_txi$length[which(rownames(r_txi$length) %in% tokeep), r_meta$Seq_ID.V2]

t_txi$abundance <- t_txi$abundance[which(rownames(t_txi$abundance) %in% tokeep), t_meta$Seq_ID.V2]
t_txi$counts    <- t_txi$counts[which(rownames(t_txi$counts) %in% tokeep), t_meta$Seq_ID.V2]
t_txi$length    <- t_txi$length[which(rownames(t_txi$length) %in% tokeep), t_meta$Seq_ID.V2]


# DeSeq2
#nocovariate
dds_t <- DESeqDataSetFromTximport(t_txi, colData = t_meta, design= ~PC1 + consensuscluster)
dds_r <- DESeqDataSetFromTximport(r_txi, colData = r_meta, design= ~PC1 + consensuscluster)


dds_t <- DESeq(dds_t)
dds_r <- DESeq(dds_r)

res_t <- results(dds_t, contrast = c("consensuscluster", "2", "1"))  
res_t$qvalue <- NA
res_t$qvalue[!is.na(res_t[, "padj"])] <- qvalue(res_t[!is.na(res_t[, "padj"]), "pvalue"])$qvalues
res_t <- as.data.frame(res_t)


res_r <- results(dds_r, contrast = c("consensuscluster", "2", "1"))  # only two levels, no need to specify a contrast. The second level is the reference (consensuscluster1)
res_r$qvalue <- NA
res_r$qvalue[!is.na(res_r[, "padj"])] <- qvalue(res_r[!is.na(res_r[, "padj"]), "pvalue"])$qvalues
res_r <- as.data.frame(res_r)

nrow(subset(res_t, res_t$qvalue < 0.05))

nrow(subset(res_r, res_r$qvalue < 0.05))

