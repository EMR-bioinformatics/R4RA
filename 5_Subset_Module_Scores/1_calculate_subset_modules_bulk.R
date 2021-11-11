rm(list=ls())

library(Seurat)
library(patchwork)
library(ggplot2)
library(ggpubr)
library(reshape2)
library(MASS) 
library(sfsmisc)
source("./plot_functions.R")

outdir <- "./output/"
dir.create(outdir, showWarnings = F, recursive = T)
median_based <- F

# gexData <- tmm_edgeR_log2
# ....
# ....
load("./subsetModuleData_BulkRNAseq.RData")

######################
##  module scores  ##
#####################

r4ra_seurat <- CreateSeuratObject(counts = gexData, min.cells = 0, min.features = 0, project = "bulkRNAseq")
for(ct in c(unique(mygenes$cell.type))){
  mysubset <- mygenes[mygenes$cell.type==ct,]
  if(ct=="Macrophage_Ale"){
    mysubset <- lapply(split(mysubset$gene, mysubset$clusterClean), function(x) x[1:length(x)])
  }else{
    mysubset <- lapply(split(mysubset$gene, mysubset$clusterClean), function(x) x[1:5])  
  }
  warning(ifelse(median_based,"Median based module score","AddModuleScore/Seurat based module score"))
  if(median_based){
    for(m in 1:length(mysubset)){
      mx <- paste0(names(mysubset)[m],m)
      mygex <- r4ra_seurat@assays$RNA@data
      mygex <- mygex[rownames(mygex) %in% mysubset[[m]],]
      r4ra_seurat[[mx]] <-  colMeans(mygex)
    }
  }else{
    r4ra_seurat <- AddModuleScore(r4ra_seurat, features = mysubset, k=F, ctrl.size = 5, name = names(mysubset))  
  }
}

names(r4ra_seurat@meta.data)[grep("SC_",names(r4ra_seurat@meta.data))] <- 
  sapply(grep("SC_",names(r4ra_seurat@meta.data), value = T), function(x) substr(x,1,nchar(x)-1))

names(r4ra_seurat@meta.data)[grep("Cluster",names(r4ra_seurat@meta.data))] <- 
  sapply(grep("Cluster",names(r4ra_seurat@meta.data), value = T), function(x) substr(x,1,nchar(x)-1))

# normalize the scores between 0-1
for(sset in unique(mygenes$clusterClean)){
  myscore <- as.vector(r4ra_seurat[[sset]])
  mymin <- min(myscore)
  mymin <- ifelse(mymin > 0, -1 * mymin, abs(mymin))
  myscore <- myscore + mymin
  mymax <- max(myscore)
  myscore <- myscore/mymax
  r4ra_seurat[[sset]] <- myscore
}

## baseline samples
r4ra.meta.bl <- r4ra.meta[which(r4ra.meta$Visit==3),]
gexData_bl <- gexData[,colnames(gexData) %in% r4ra.meta.bl$noPrefix]
r4ra.switch.bl <- r4ra.switch[match(r4ra.meta.bl$Patient.I.D., r4ra.switch$Patient.I.D.),]
all(r4ra.switch.bl$Patient.I.D.==r4ra.meta.bl$Patient.I.D.)
all(colnames(gexData)==r4ra.meta$noPrefix)

res_list_bl_anova <- list()
res_list_bl <- list()

for(ct in unique(mygenes$cell.type)){
  
  mysubset <- mygenes[mygenes$cell.type==ct,]
  mysubset <- split(mysubset$gene, mysubset$clusterClean)
  r4ra_seurat <- r4ra_seurat[,which(colnames(r4ra_seurat) %in% colnames(gexData_bl))]
  
  all(colnames(r4ra_seurat)==r4ra.meta.bl$noPrefix[match(r4ra.switch.bl$Patient.I.D.,r4ra.meta.bl$Patient.I.D.)])
  
  for(sset in names(mysubset)){
    
    g1 <- r4ra_seurat[[sset]][r4ra.switch.bl$Responses %in% c("ResponderRIT"),]
    g2 <- r4ra_seurat[[sset]][r4ra.switch.bl$Responses %in% c("ResponderTOC"),]
    g3 <- r4ra_seurat[[sset]][r4ra.switch.bl$Responses %in% "Refractory",]
    mydataKW <- data.frame("score"=c(g1,g2,g3),"group"=c(rep("proRIT",length(g1)), rep("proTOC",length(g2)), rep("proRefactory",length(g3))), stringsAsFactors = T)
    mydataKW <- mydataKW[!is.na(mydataKW$score),]
    res_kw<- kruskal.test(score ~ group, data = mydataKW)
    res_list_bl_anova[[paste0(ct,"-",sset,"-RTX_&_TOC_@@@_ProTOC/RIT_vs_Refractory")]]$pvalue <-res_kw$p.value
    res_list_bl_anova[[paste0(ct,"-",sset,"-RTX_&_TOC_@@@_ProTOC/RIT_vs_Refractory")]]$g1 <- g1
    res_list_bl_anova[[paste0(ct,"-",sset,"-RTX_&_TOC_@@@_ProTOC/RIT_vs_Refractory")]]$g2 <- g2
    res_list_bl_anova[[paste0(ct,"-",sset,"-RTX_&_TOC_@@@_ProTOC/RIT_vs_Refractory")]]$g3 <- g3
    res_list_bl_anova[[paste0(ct,"-",sset,"-RTX_&_TOC_@@@_ProTOC/RIT_vs_Refractory")]]$fc <- NA
    
    g1 <- r4ra_seurat[[sset]][r4ra.switch.bl$Responses %in% c("ResponderRIT","ResponderTOC"),]
    g2 <- r4ra_seurat[[sset]][r4ra.switch.bl$Responses %in% "Refractory",]
    res1 <- wilcox.test(g1, g2, alternative = "two.sided", paired = F)
    res_list_bl[[paste0(ct,"-",sset,"-RTX_&_TOC_@@@_Responders_vs_Refractory")]]$pvalue <-res1$p.value
    res_list_bl[[paste0(ct,"-",sset,"-RTX_&_TOC_@@@_Responders_vs_Refractory")]]$g1 <- g1
    res_list_bl[[paste0(ct,"-",sset,"-RTX_&_TOC_@@@_Responders_vs_Refractory")]]$g2 <- g2
    res_list_bl[[paste0(ct,"-",sset,"-RTX_&_TOC_@@@_Responders_vs_Refractory")]]$fc <- median(g1,na.rm = T)/median(g2, na.rm = T)  
    
    g1 <- r4ra_seurat[[sset]][r4ra.switch.bl$responderTOC %in% c("Responder"),]
    g2 <- r4ra_seurat[[sset]][r4ra.switch.bl$responderTOC %in% "Non.Responder",]
    res2 <- wilcox.test(g1, g2, alternative = "two.sided", paired = F)
    res_list_bl[[paste0(ct,"-",sset,"-TOC_@@@_Responders_vs_Non_responders")]]$pvalue <-res2$p.value
    res_list_bl[[paste0(ct,"-",sset,"-TOC_@@@_Responders_vs_Non_responders")]]$g1 <- g1
    res_list_bl[[paste0(ct,"-",sset,"-TOC_@@@_Responders_vs_Non_responders")]]$g2 <- g2
    res_list_bl[[paste0(ct,"-",sset,"-TOC_@@@_Responders_vs_Non_responders")]]$fc <- median(g1,na.rm = T)/median(g2, na.rm = T)
    
    g1 <- r4ra_seurat[[sset]][r4ra.switch.bl$responderRIT %in% c("Responder"),]
    g2 <- r4ra_seurat[[sset]][r4ra.switch.bl$responderRIT %in% "Non.Responder",]
    res3 <- wilcox.test(g1, g2, alternative = "two.sided", paired = F)
    res_list_bl[[paste0(ct,"-",sset,"-RTX_@@@_Responders_vs_Non_responders")]]$pvalue <- res3$p.value
    res_list_bl[[paste0(ct,"-",sset,"-RTX_@@@_Responders_vs_Non_responders")]]$g1 <- g1
    res_list_bl[[paste0(ct,"-",sset,"-RTX_@@@_Responders_vs_Non_responders")]]$g2 <- g2
    res_list_bl[[paste0(ct,"-",sset,"-RTX_@@@_Responders_vs_Non_responders")]]$fc <- median(g1,na.rm = T)/median(g2, na.rm = T)
  
  }
}

sapply(res_list_bl, function(x) x$pvalue)
sapply(res_list_bl_anova, function(x) x$pvalue)
r4ra_seurat_responders <- r4ra_seurat

## refractory

aaa <- intersect(grep("RTX_&_TOC_@@@_Responders_vs_Refractory",names(res_list_bl)), grep("Fibroblast",names(res_list_bl)))

p_response_fib1 <- do.ridgePlot(res_list = res_list_bl[aaa], 
                                mycolors = c("#4c4cca", "#e90000","#cbae01","#07ca07"), mylevels = c("Responders","Refractory"), position = NULL, violin = F, densityPlot = T)

pdf(file = paste0(outdir,"R4RA_module_score_fibro_RTX&TOC_refractory.pdf"),  height = 6, width = 8)
print(p_response_fib1)
dev.off()

aaa <- intersect(grep("RTX_&_TOC_@@@_Responders_vs_Refractory",names(res_list_bl)), grep("Macrophage-",names(res_list_bl)))

p_response_mono <- do.ridgePlot(res_list = res_list_bl[aaa], 
                                mycolors = c("#4c4cca", "#e90000","#cbae01","#07ca07"), mylevels = c("Responders","Refractory"), position = NULL, violin = F, densityPlot = T)

pdf(file = paste0(outdir,"R4RA_module_score_macrophage_RTX&TOC_refractory.pdf"),  height = 6, width = 9)
print(p_response_mono)
dev.off()

########################
## clinical variables ##
########################

subsetLabels <- list(
  "SC-F1" = "CD34+ sublining",
  "SC-F2" = "HLA+ sublining",
  "SC-F3" = "DKK3+ sublining",
  "SC-F4" = "CD55+ lining",

  "SC-M1" = "IL1B+ pro-inflammatory",
  "SC-M2" = "NUPR1+",
  "SC-M3" = "C1QA+",
  "SC-M4" = "IFN-activated",

  "SC-T1" = "CCR7+ CD4+",
  "SC-T2" = "FOXP3+ Tregs",
  "SC-T3" = "PD-1+ Tph/Tfh",
  "SC-T4" = "GZMK+ CD8+",
  "SC-T5" = "GNLY+ GZMB+",
  "SC-T6" = "GZMK+/GZMB+",

  "SC-B1" = "IGHD+ CD270 naive",
  "SC-B2" = "IGHG3+ CD27- memory",
  "SC-B3" = "Autoimmune associated",
  "SC-B4" = "Plasmablasts",

  "Cluster-0" = "TREM2low (0)",
  "Cluster-1" = "TREM2high (1)",
  "Cluster-2" = "FOLR2+   ID2+ (2)",
  "Cluster-3" = "FOLR2high LYVE1+ (3)",
  "Cluster-8" = "FOLR2+ ICAM1+ (8)",
  "Cluster-4" = "HLAhigh CLEC10A+ (4)",
  "Cluster-7" = "HLAhigh ISG15+ (7)",
  "Cluster-5" =  "CD52+  S100A12+ (5)",
  "Cluster-6" = "CD52+  SSP1+ (6)")

load("./clinicalVarCorPlots.RData")

mycors <- list()
myplot <- list()

for(s in mysubsets){
  myscore <- r4ra_seurat_responders[[s]]
  idx <-  match(meta$noPrefix,rownames(myscore))
  
  for(m in c(histology1, histology2, clinical1, clinical2, response1, response2)){
    corcoeff <- cor.test(as.numeric(meta[[m]]), as.numeric(myscore[,1][idx]), use = "na.or.complete", method = "s")
    mycors[[paste0(s,"@@@",m)]]$pvalue <- corcoeff$p.value
    mycors[[paste0(s,"@@@",m)]]$coef <- corcoeff$estimate
  }
}

mycorDF <- data.frame("Subset"=sapply(strsplit(names(mycors), split = "@@@"), function(x) x[1]),
                      "ClinicalVariable"=sapply(strsplit(names(mycors), split = "@@@"), function(x) x[2]),
                      "Correlationcoefficient"=sapply(mycors, function(x) x$coef), 
                      "Correlationpvalue"=sapply(mycors, function(x) x$pvalue), 
                      stringsAsFactors = F) 

# boxplot histology1, clinical1, response1, 
summary(meta[,histology1])
summary(meta[,clinical1])
summary(meta[,response1])

meta$Cell.Type.V2[meta$Cell.Type.V2=="GC"] <- NA 
selectedvars <- c("Pathotype.V2", histology2, response2, clinical2)

# correlplot
for(s in mysubsets[grep("cluster|T|B",mysubsets, invert = T, ignore.case = T)]){
  
  ss_label <- gsub("_","-", s)
  ss_label <- paste0(ss_label,": ",subsetLabels[[ss_label]])
  
  par(mfrow=c(2,2))
  myscore <- r4ra_seurat_responders[[s]]
  idx <-  match(meta$noPrefix,rownames(myscore))
  for(m in selectedvars){
    # corcoeff <- cor.test(as.numeric(meta[[m]]), as.numeric(myscore[,1][idx]), use = "na.or.complete")
    if(m=="Pathotype.V2"){    
      rsl <- rlm(as.numeric(meta[[m]]) ~  as.numeric(myscore[,1][idx]))
      correlplot(x = as.numeric(meta[[m]]), extraLabel = "", xaxt = "n", 
                 y=as.numeric(myscore[,1][idx]), 
                 prlm = signif(f.robftest(rsl)$p.value, digits = 1),
                 method="spearman", q=NA, xlab=gsub("[.]"," ",m), ylab=ss_label, cex.axis=1, cex.lab=1.2, cex.pch=0.9, cex.p = 1)
      axis(1, at = c(1,2,3), labels=levels(meta[[m]]))
    }else{
      rsl <- rlm(as.numeric(meta[[m]]) ~  as.numeric(myscore[,1][idx]))
      correlplot(x = as.numeric(meta[[m]]), extraLabel = "", 
                 y=as.numeric(myscore[,1][idx]), 
                 prlm = signif(f.robftest(rsl)$p.value, digits = 1),
                 method="spearman", q=NA, xlab=gsub("[.]"," ",m), ylab=ss_label, cex.axis=1, cex.lab=1.2, cex.pch=0.9, cex.p = 1)
    }
  }
}


mycorDF <- mycorDF[mycorDF$ClinicalVariable %in% selectedvars,]

mycolnames <- list(
  "SC-F1" = "CD34+ sublining",
  "SC-F2" = "HLA+ sublining",
  "SC-F3" = "DKK3+ sublining",
  "SC-F4" = "CD55+ lining",
  
  "SC-M1" = "IL1B+ pro-inflammatory",
  "SC-M2" = "NUPR1+",
  "SC-M3" = "C1QA+",
  "SC-M4" = "IFN-activated")

mycorDF$ClinicalVariable <- gsub("[.]"," ",mycorDF$ClinicalVariable)
mycorDF <- mycorDF[order(paste0(mycorDF$Subset,"-",mycorDF$ClinicalVariable)),]
mycorDF <- mycorDF[grep("cluster|T|B",mycorDF$Subset, invert = T, ignore.case = T),]
mycorDF$Subset <- gsub("_","-",mycorDF$Subset)
mycorDF$Subset <- paste0(mycorDF$Subset,": ",unlist(mycolnames[mycorDF$Subset]))
mycorDF$Correlationcoefficient <- round(mycorDF$Correlationcoefficient, 2)
mycorDF$Correlationpvalue <- signif(mycorDF$Correlationpvalue, digits = 2)
write.table(mycorDF, file = paste0(outdir,"moduleScores_vs_clinicalVariables_long.tsv"), sep = "\t", row.names = F)
write.table(dcast(mycorDF, formula = Subset~ClinicalVariable, value.var = "Correlationcoefficient"), file = paste0(outdir,"moduleScores_vs_clinicalVariables_corCoef_wide.tsv"), sep = "\t", row.names = T)
write.table(dcast(mycorDF, formula = Subset~ClinicalVariable, value.var = "Correlationpvalue"), file = paste0(outdir,"moduleScores_vs_clinicalVariables_pvalue_wide.tsv"), sep = "\t", row.names = T)

mycorDF$Correlationcoefficient2 <- mycorDF$Correlationcoefficient
mycorDF$Correlationcoefficient2 <- paste0(mycorDF$Correlationcoefficient2, ifelse(mycorDF$Correlationpvalue < 0.001,"***",ifelse(mycorDF$Correlationpvalue < 0.01,"**",ifelse(mycorDF$Correlationpvalue < 0.05,"*",""))))

myorder <- c("Pathotype V2", "Krenn Score",
             "CD3 V2","CD20 V2","CD138 V2","CD68L V2","CD68SL V2",
             "SJC","TJC", "VAS", "VAS pain","ESR","CRP",
             "CDAI","DAS28 ESR","DAS28 CRP", "DAS2c",
             "HAQ Score","Change in CDAI from baseline V7","Change in DAS28 ESR from baseline V7",
             "Change in DAS28 CRP from baseline V7","Change in DAS2c from baseline V7",
             "delta SySc","delta CD3","delta CD20","delta CD138","delta CD68L","delta CD68SL")
any(is.na(match(unique(mycorDF$ClinicalVariable), myorder)))
mycorDF$ClinicalVariable <- factor(mycorDF$ClinicalVariable, levels=rev(myorder))

pdf(file = paste0(outdir,"moduleScores_vs_clinicalVariables.pdf"),  height = 10, width = 12)
do.correlation.heatmap(mycorDF)
dev.off()

