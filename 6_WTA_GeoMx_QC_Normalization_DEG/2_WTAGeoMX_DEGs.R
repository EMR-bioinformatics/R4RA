rm(list=ls())

library(qvalue)
library(DESeq2)
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(reshape2)
library(ggrepel)
library(grid)
library(IHW)
library(glmmSeq)
library(easylabel)

myannot <- read.delim("./TSCOL_0257 Annotation file.csv", sep = "\t", stringsAsFactors = F)
mypathotypes <- read.delim("./GeoMx DSP Samples_Batch 2_Decoded_update_06082021_AS_patients_and_pathotypes.csv", sep = "\t", stringsAsFactors = F)

# Added separately 
# Otherwise, if I do >>> myannot <- cbind(myannot,mypathotypes[match(myannot$Scan.name, mypathotypes$Anon.ID),])
# DESeq2 converts colnames of counts to some numbers
idx <- match(myannot$Scan.name, mypathotypes$Anon.ID)
myannot$Anon.ID <- mypathotypes$Anon.ID[idx]
myannot$Pathotype <- mypathotypes$Pathotype[idx]
myannot$Sample.ID <- mypathotypes$Sample.ID[idx]
myannot$Status <- mypathotypes$Status[idx]
myannot$Diagnosis <- mypathotypes$Diagnosis[idx]
myannot$Treatment <- mypathotypes$Treatment[idx]

myData <- read.delim("./TSCOL_0257 Initial Dataset.csv", sep = "\t", stringsAsFactors = F, as.is = T, check.names = F)
sort(table(myData$TargetName))[sort(table(myData$TargetName))>1]
myData <- myData[myData$TargetName!="NegProbe-WTX",]
rownames(myData) <- myData$TargetName
myData <- myData[,-1]
myData <- myData[,match(myannot$SegmentDisplayName, colnames(myData))]

myannot$ROIs <- paste0("ROI_",myannot$ROI..label.)

####################################################
############ Normalization & QC Chekcs #############
####################################################

dds <- DESeqDataSetFromMatrix(countData = myData, colData = myannot, design = ~ Location)
dds <- estimateSizeFactors(dds)
dds <- DESeq(dds)

# ROI
boxplot(log2(myData[,order(myannot$ROI_ID)]+1), col = rainbow(length(unique(myannot$ROIs)))[as.factor(myannot$ROIs[order(myannot$ROI_ID)])])
samps4myMeans <- split(myannot$SegmentDisplayName, myannot$ROIs)
myDataMean <- cbind(sapply(samps4myMeans[1:6], function(x) rowMeans(myData[, colnames(myData) %in%x])), myData[,unlist(samps4myMeans[7:8])])
boxplot(log2(myDataMean+1), las=2, pars=list(par(mar=c(12,2,2,2))))

# Treatment
boxplot(log2(myData[,order(myannot$Treatment)]+1), col = rainbow(length(unique(myannot$Treatment)))[as.factor(myannot$Treatment[order(myannot$Treatment)])])
samps4myMeans <- split(myannot$SegmentDisplayName, myannot$Treatment)
myDataMean <- sapply(samps4myMeans, function(x) rowMeans(myData[, colnames(myData) %in%x]))
boxplot(log2(myDataMean+1), las=2, pars=list(par(mar=c(12,2,2,2))))

# Patient
boxplot(log2(myData[,order(myannot$Sample.ID)]+1), col = rainbow(length(unique(myannot$Sample.ID)))[as.factor(myannot$Sample.ID[order(myannot$Sample.ID)])])
samps4myMeans <- split(myannot$SegmentDisplayName, myannot$Sample.ID)
myDataMean <- sapply(samps4myMeans, function(x) rowMeans(myData[, colnames(myData) %in%x]))
boxplot(log2(myDataMean+1), las=2, pars=list(par(mar=c(12,2,2,2))))

# Pathotype
boxplot(log2(myData[,order(myannot$Pathotype)]+1), col = rainbow(length(unique(myannot$Pathotype)))[as.factor(myannot$Pathotype[order(myannot$Pathotype)])])
samps4myMeans<- split(myannot$SegmentDisplayName, myannot$Pathotype)
myDataMean <- sapply(samps4myMeans, function(x) rowMeans(myData[, colnames(myData) %in%x]))
boxplot(log2(myDataMean+1), las=2, pars=list(par(mar=c(8,2,2,2))))

split(myannot$Sample.ID, myannot$ROIs)

myDataNorm <- read.delim("./TSCOL_0257 Q3 Normalized.csv", sep = "\t", stringsAsFactors = F, as.is = T, check.names = F, row.names = 1)
normQ3 <- as.matrix(myDataNorm[,match(myannot$SegmentDisplayName, colnames(myDataNorm))])
# ROI
boxplot(log2(normQ3[,order(myannot$ROI_ID)]+1), col = rainbow(length(unique(myannot$ROIs)))[as.factor(myannot$ROIs[order(myannot$ROI_ID)])])
# Treatment
boxplot(log2(normQ3[,order(myannot$Treatment)]+1), col = rainbow(length(unique(myannot$Treatment)))[as.factor(myannot$Treatment[order(myannot$Treatment)])])
# Patient
boxplot(log2(normQ3[,order(myannot$Sample.ID)]+1), col = rainbow(length(unique(myannot$Sample.ID)))[as.factor(myannot$Sample.ID[order(myannot$Sample.ID)])])
# Pathotype
boxplot(log2(normQ3[,order(myannot$Pathotype)]+1), col = rainbow(length(unique(myannot$Pathotype)))[as.factor(myannot$Pathotype[order(myannot$Pathotype)])])

samps4myMeans <- split(myannot$SegmentDisplayName, myannot$Sample.ID)
myDataMean <- sapply(samps4myMeans, function(x) rowMeans(normQ3[, colnames(normQ3) %in%x]))
boxplot(log2(myDataMean+1), las=2, pars=list(par(mar=c(12,2,2,2))))

normDDS <- log2(counts(dds, normalized=TRUE) + 1)
boxplot(normDDS[,order(myannot$ROI_ID)], col = rainbow(length(unique(myannot$ROIs)))[as.factor(myannot$ROIs[order(myannot$ROI_ID)])])
boxplot(normDDS[,order(myannot$Treatment)], col = rainbow(length(unique(myannot$Treatment)))[as.factor(myannot$Treatment[order(myannot$Treatment)])])
boxplot(normDDS[,order(myannot$Sample.ID)], col = rainbow(length(unique(myannot$Sample.ID)))[as.factor(myannot$Sample.ID[order(myannot$Sample.ID)])])
boxplot(normDDS[,order(myannot$Pathotype)], col = rainbow(length(unique(myannot$Pathotype)))[as.factor(myannot$Pathotype[order(myannot$Pathotype)])])
samps4myMeans <- split(myannot$SegmentDisplayName, myannot$Sample.ID)
myDataMean <- sapply(samps4myMeans, function(x) rowMeans(normDDS[, colnames(normDDS) %in%x]))
boxplot(myDataMean, las=2, pars=list(par(mar=c(12,2,2,2))))

# blind=T sets  design(object) <- ~1, higher nsub makes the normalization worse
vst <- vst(dds, blind=T,  nsub = 1000)
vst <- assay(vst)
boxplot(vst[,order(myannot$ROI_ID)], col = rainbow(length(unique(myannot$ROIs)))[as.factor(myannot$ROIs[order(myannot$ROI_ID)])])
boxplot(vst[,order(myannot$Treatment)], col = rainbow(length(unique(myannot$Treatment)))[as.factor(myannot$Treatment[order(myannot$Treatment)])])
boxplot(vst[,order(myannot$Sample.ID)], col = rainbow(length(unique(myannot$Sample.ID)))[as.factor(myannot$Sample.ID[order(myannot$Sample.ID)])])
boxplot(vst[,order(myannot$Pathotype)], col = rainbow(length(unique(myannot$Pathotype)))[as.factor(myannot$Pathotype[order(myannot$Pathotype)])])
samps4myMeans <- split(myannot$SegmentDisplayName, myannot$Sample.ID)
myDataMean <- sapply(samps4myMeans, function(x) vst[, colnames(vst) %in% x])
boxplot(myDataMean, las=2, pars=list(par(mar=c(12,2,2,2))))

# dispersions for glmmSeq
dds <- DESeqDataSetFromMatrix(countData = myData, colData = myannot, design = ~ 1)
dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)
dispersionsWTA <- dispersions(dds)
names(dispersionsWTA) <- rownames(myData)
rm(dds)

####################################################
################### DEG Analysis ###################
####################################################

######## DESeq2 -  DEGs for Response ######## 

res_responders <- list()

# 1)
dds <- DESeqDataSetFromMatrix(countData = myData, colData = myannot, design = ~ Location + Response)
dds <- estimateSizeFactors(dds)
dds <- DESeq(dds)
dssDF <- as.data.frame(results(dds))
dssDF$qvalue[!is.na(dssDF$padj)] <- qvalue::qvalue(dssDF$pvalue[!is.na(dssDF$padj)])$qvalues
# The standard IHW method presented below controls the FDR by using a weighted Benjamini-Hochberg procedure with data-driven weights. 
ihwRes <- ihw(pvalue ~ baseMean,  data = dssDF, alpha = 0.1)
dssDF$padj_ihwBH <- adj_pvalues(ihwRes)
res_responders[["All_samples_Responder_vs_Refractory_with_location_as_covariate"]] <- dssDF

# 2) 
for(ml in unique(myannot$Location)){
  dds <- DESeqDataSetFromMatrix(countData = myData[,myannot$Location %in% ml], colData = myannot[myannot$Location %in% ml,], design = ~ Response)
  dds <- estimateSizeFactors(dds)
  dds <- DESeq(dds)
  dssDF <- as.data.frame(results(dds))
  dssDF$qvalue[!is.na(dssDF$padj)] <- qvalue::qvalue(dssDF$pvalue[!is.na(dssDF$padj)])$qvalues
  ihwRes <- ihw(pvalue ~ baseMean,  data = dssDF, alpha = 0.1)
  dssDF$padj_ihwBH <- adj_pvalues(ihwRes)
  res_responders[[paste0(ml,"_Responder_vs_Refractory", collapse = "")]] <- dssDF
  
}

####################################################
######################### PCAs  ####################
####################################################
pdf("./WTA_GeoMx_PCA_of_subset_genes.pdf",  height = 8, width = 8)

scale <- F
all(rownames(myData)==rownames(vst))
all(colnames(vst)==myannot$SegmentDisplayName)
vst.pca <- vst[which(!apply(myData,1, function(x) all(x<=2))),]
pc <- prcomp(t(vst.pca), scale. = scale)
varExp <- round(pc$sdev^2/sum(pc$sdev^2),digits = 2) * 100
# plot(pc$x[,c(1,2)], main = "PCA", frame = FALSE, col=outcols)

x <- 1
y <- 2
PCs <- paste0("PC",c(x,y))
pc_df <- as.data.frame(pc$x[,PCs])
colnames(pc_df) <- c("pc_x","pc_y")
pc_df$Response <- as.factor(myannot$Response)
pc_df$Location <- as.factor(myannot$Location)
pc_df$Scan.name <- as.integer(gsub("Refrac ","",gsub("Resp ","",gsub(" Blk 2","",myannot$Scan.name))))
p<-ggplot(pc_df,aes(x=pc_x,y=pc_y,color=Response, shape=Location, label = Scan.name))
p <- p+geom_point(size=3, alpha = 0.7) 

p <- p + theme_classic() + guides(colour = guide_legend(override.aes = list(size=6)))
p <- p + theme(strip.text.x = element_text(face = "italic"), 
               axis.title.x = element_text(size=14, face="bold"),
               axis.title.y = element_text(size=14, face="bold"),
               legend.title = element_text(size=14, face="bold"),
               legend.text = element_text(size=14, face="bold"))
p <- p + xlab(paste0(PCs[1], " variance %", varExp[1]))
p <- p + ylab(paste0(PCs[2], " variance %", varExp[2]))
p <- p + labs(col="", shape ="", subtitle="") + theme(plot.margin = unit(c(1,1,1,1), "cm")) +
  scale_color_manual(values = c("#e90000", "#4c4cca")) 
p <- p + geom_text_repel(box.padding = 0.2, max.overlaps = Inf, show.legend = F) +
  theme(legend.position = c(0.15, 0.75),  legend.key.size = unit(0.2, "cm")) +
  geom_label(x=-85, y=125, label="Patient IDs are shown in numbers", color="black", size=4 , angle=45, fontface="bold" )

grid.force()
p

dev.off()

####################################################
####################### Volcano  ###################
####################################################

selectedGenes <- c("ASPN", "SPP1", "DKK3", "S100A10", "MAMDC2",  "TMSB10",  "ABHD2", "MT2A", "PRELP","GAPDH")

# response
myDEGs <- res_responders$All_samples_Responder_vs_Refractory_with_location_as_covariate
easyVolcano(myDEGs, useQ = TRUE, fullname = T, x = "log2FoldChange", y = NULL, showOutliers = T,
            fccut = 0.5, fdrcutoff = 0.05, showCounts = T,
            markerSize = 8, outline_lwd = 0.5,  outline_col = "white", outlier_pch = 3, ylim = c(0,50),
            startLabels = selectedGenes,
            scheme=c('darkgrey', 'blue', 'lightblue', 'orange', 'red')) 

myDEGs <- res_responders$SubLining_Responder_vs_Refractory
easyVolcano(myDEGs, useQ = TRUE, fullname = T, x = "log2FoldChange", y = NULL, showOutliers = T,
            fccut = 0.5, fdrcutoff = 0.05, showCounts = T,
            markerSize = 8, outline_lwd = 0.5,  outline_col = "white", outlier_pch = 3, ylim = c(0,50),
            startLabels = selectedGenes,
            scheme=c('darkgrey', 'blue', 'lightblue', 'orange', 'red')) 


myDEGs <- res_responders$Lining_Responder_vs_Refractory
easyVolcano(myDEGs, useQ = TRUE, fullname = T, x = "log2FoldChange", y = NULL, showOutliers = T,
            fccut = 0.5, fdrcutoff = 0.05, showCounts = T,
            markerSize = 8, outline_lwd = 0.5,  outline_col = "white", outlier_pch = 3, ylim = c(0,50),
            startLabels = selectedGenes,
            scheme=c('darkgrey', 'blue', 'lightblue', 'orange', 'red')) 

myDEGs <- res_responders$Aggregate_Responder_vs_Refractory
easyVolcano(myDEGs, useQ = TRUE, fullname = T, x = "log2FoldChange", y = NULL, showOutliers = T,
            fccut = 0.5, fdrcutoff = 0.05, showCounts = T,
            markerSize = 8, outline_lwd = 0.5,  outline_col = "white", outlier_pch = 3, ylim = c(0,50),
            startLabels = selectedGenes,
            scheme=c('darkgrey', 'blue', 'lightblue', 'orange', 'red')) 


