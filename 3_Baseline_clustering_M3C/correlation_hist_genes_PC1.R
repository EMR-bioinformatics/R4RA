### correlate gene expression, histology scores and PC1
library(edgeR)
library(dplyr)
library(corrplot)

meta <- read.csv('./meta_data.csv')
meta <- subset(meta, meta$Visit == 3)

PC1 <- readRDS('/PCscores_R4RA_forDESeq2.rds')

load(".//R4RA_V2_data_trimmed.RData")

### parameters

data <- txiready$vst

## remove outliers
outliers <- read.csv('./outliers.csv',sep=',')
outliers$ID <- as.character(outliers$ID)
data <- data[ , -which(names(data) %in% outliers$ID)] # 12 initial outliers

setdiff(colnames(data),meta$Seq_ID.V2)

meta <- subset(meta, meta$Seq_ID.V2 %in% colnames(data))
data <- data[,meta$Seq_ID.V2]

data_backup <- data # save for boxplots

meta <- meta %>% select(c("Seq_ID.V2", "CD20.V2", "CD138.V2", "CD68L.V2", "CD68SL.V2", "CD3.V2" )) 

## calculate PC and then add as covariate
PC1$Seq_ID.V2 <- rownames(PC1)
PC1 <- PC1 %>% select(c("Seq_ID.V2", "PC1"))
coldata <-right_join(meta, PC1, by='Seq_ID.V2')

data2 <- data %>% filter(rownames(data) %in% c('MS4A1', 'CD79A', 'CD79B', 'PIK3CA', 'BTK', 'SYK', 'IL6R', 'IL6', 'IL6ST', 'JAK1', 'JAK2', 'STAT3' ))
data2 <- as.data.frame(t(data2))
data2$Seq_ID.V2 <- rownames(data2)
coldata2 <-right_join(coldata, data2, by='Seq_ID.V2')
coldata2 <- coldata2 %>% rename(CD20 = CD20.V2, CD138 = CD138.V2, CD68L = CD68L.V2, CD68SL = CD68SL.V2, CD3 = CD3.V2)
rownames(coldata2) <- coldata2$Seq_ID.V2
coldata2 <- coldata2[,-1]
M <- cor(coldata2, use="pairwise.complete.obs") 
testRes <- cor.mtest(coldata2, conf.level = 0.95)
p1 <- corrplot(M, p.mat=testRes$p, type= 'lower', insig = 'p-value', method = 'circle', diag = FALSE, tl.col = 'black')


