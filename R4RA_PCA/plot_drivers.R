## explore drivers of variation in data

library(bioplotr)
library(dplyr)
library(edgeR)
library(ggplot2)
library(purrr)

load("./R4RA_V2_data_trimmed.RData")
source('./R4RA_PCA/plotdrivers_function_AS.R')

## load all synovium data, baseline (visit 3 clinical, visit 2 biopsy)

data <- txiready$vst

outliers <- read.csv('./outliers.csv',sep=',')
outliers$ID <- as.character(outliers$ID)
data <- data[ , -which(names(data) %in% outliers$ID)] # 12 initial outliers



## get meta data and select visit 3 (baseline)
meta <- read.csv('./meta_data.csv')

## select visit 3, reformat
meta <- subset(meta, meta$Visit == 3)

meta <- subset(meta, meta$Seq_ID.V2 %in% colnames(data))



meta2 <- meta[,c('TJC','SJC','ESR','CRP', 'Arthritis.Activity',
                 'CDAI.response.status.V7', 'Gender','Age','Ethnicity',
                 'CCP_status.V1','RF_status.V1', 
                'Cell.Type.V2','Pathotype.V2')]

colnames(meta2)[colnames(meta2)=='Pathotype.V2'] <- 'Pathotype'
colnames(meta2)[colnames(meta2)=='Cell.Type.V2'] <- 'Cell Type'
colnames(meta2)[colnames(meta2)=='CCP_status.V1'] <- 'CCP_status'
colnames(meta2)[colnames(meta2)=='RF_status.V1'] <- 'RF_status'
colnames(meta2)[colnames(meta2)=='CDAI.response.status.V7'] <- 'CDAI 50% improvement'


data <- data[,meta2$Seq_ID.V2]

p <- plotdrivers(data, meta2, label=TRUE, n.pc=10, p.adj = NULL, alpha = 0.05)
p ## non adjusted p-value


p2 <- plotdrivers(data, meta2, label=TRUE, n.pc=10, p.adj = 'fdr', alpha = 0.05)
p2 ## adjusted p-value


