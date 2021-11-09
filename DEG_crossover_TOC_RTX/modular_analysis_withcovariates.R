library(qusage)
library(DESeq2)
library(edgeR)
library(dplyr)
source('./DEG_crossover_TOC_RTX/load_mods_qmod.R')
source('./DEG_crossover_TOC_RTX/qmod_AS.R')

mod.list <- c('reduced_wgcna', 'wgcna',  'LI_reduced','LI','singlecell') #, LI = tmod
data <- dds_t 
qusageres_list_CDAI50_TOC <- lapply(mod.list, function(module){
  c2.indices <- loadmodsforqusage(load = module)
  qusageres <- qmod(fit=data, filter= c(1,1), coef = 2,  geneSets =  c2.indices)
  print(module)
  print(nrow(subset(qusageres, qusageres$qval < 0.05)))
  
  write.csv(qusageres, file = paste0("Qmod_",module, "_TOC_CDAI50.csv"))
  return(qusageres)
})

names(qusageres_list_CDAI50_TOC) <- unlist(mod.list)
saveRDS(qusageres_list_CDAI50_TOC, file = 'qusageres_list_CDAI50_TOC.rds')
TOC <- do.call(rbind, qusageres_list_CDAI50_TOC)




mod.list <- c('reduced_wgcna', 'wgcna',  'LI_reduced','LI') 
data <- dds_r

qusageres_list_CDAI50_RTX <- lapply(mod.list, function(module){
  c2.indices <- loadmodsforqusage(load = module)
  qusageres <- qmod(fit=data, filter= c(1,1), coef = 2,  geneSets =  c2.indices)
  print(module)
  print(nrow(subset(qusageres, qusageres$qval < 0.05)))
  
  write.csv(qusageres, file = paste0("Qmod_",module, "_RTX_CDAI50_PCA1cov.csv"))
  return(qusageres)
})

names(qusageres_list_CDAI50_RTX) <- unlist(mod.list)
saveRDS(qusageres_list_CDAI50_RTX, file = 'qusageres_list_CDAI50_RTX.rds')
RTX <- do.call(rbind, qusageres_list_CDAI50_RTX)

### plot modules on q values 
library(ggplot2)

RTX_select <- c('S80 CD8_T_cells',  'S43 M1_macrophage, chemokines and cytokines', 'S158 M2_macrophage', 'S185 Hox genes', 'S214 B_cells', 
                'S156 Dendritic_cells, MHC class I antigen presentation', 'S104 Mast_cells, IL-4 and IL-13 signalling', 'S108 CD4 and CD8_T_cells', 
                'S139 CD8 and Tph T_cells', 'S2 Interferon (IFN) signalling', 'S210 Fibroblast_1_CD55+', 'S224 PPAR signaling pathway', 'M71 enriched in antigen presentation (I)',
                'M29 proinflammatory cytokines and chemokines', 'M88.0 leukocyte migration','M89.1 putative targets of PAX3', 'M146 MHC-TLR7-TLR8 cluster', 
                'M111.1 viral sensing & immunity; IRF2 targets network (II)', 'M36 T cell surface, activation',  'S180 Fibroblast_2a_THY1+, ECM organization'
                )
TOC_select <- c('M140 extracellular matrix, complement', 'M74 transcriptional targets of glucocorticoid receptor', 'M28 antigen presentation (lipids and proteins)',
                'M78 myeloid cell cytokines, metallopeptidases and laminins', 'M46 cell division stimulated CD4+ T cells', 'S185 Hox genes', 'S99 Fibroblast_1_CD55+',
                'S200 Fibroblast_2a_THY1+', 'S51 Th17_cell differentiation', 'S214 B_cells', 'S202 Dendritic_cells', 'S104 Mast_cells, IL-4 and IL-13 signalling',
                'S224 PPAR signaling pathway')
#qusageres_list <- readRDS('./Tocllizumab/ACR20/ACR_TOCI_qusageres.rds')
masterresults <- RTX
masterresults <- masterresults[order(masterresults$p.Value),]
masterresults1 <- subset(masterresults, masterresults$pathway.name %in% RTX_select)
masterresults1 <- masterresults1[order(masterresults1$log.fold.change),]
masterresults1$pathway.name <- factor(masterresults1$pathway.name, levels=masterresults1$pathway.name)

masterresults2 <- TOC
masterresults2 <- subset(masterresults2, masterresults2$pathway.name %in% TOC_select)
masterresults2 <- masterresults2[order(masterresults2$log.fold.change),]
masterresults2$pathway.name <- factor(masterresults2$pathway.name, levels=masterresults2$pathway.name)





p1 <- ggplot(masterresults1, aes(pathway.name, log.fold.change)) +
  geom_point(aes(colour=qval),size=3) +
  scale_colour_gradient(low="red", high="black") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, colour = 'black', size = 12.5),
        axis.text.y = element_text(color = 'black', size = 12.5),
        axis.title.x = element_text(color = 'black', size = 12.5),
        axis.title.y = element_text(color = 'black', size = 12.5),
        legend.text=element_text(size=12.5),
        legend.title=element_text(size=12.5), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.margin = unit(c(1,1,1,3), "cm")) +
  ylab(expression('log'[2]*'(Responder/Non Responder)')) + 
  geom_hline(yintercept=0,linetype="dashed",color = "grey") +
  labs(colour='Q value') 



png('qusage_overall_RTX_CDAI50_combined_nocov.png', height = 16, width = 25, units = 'cm', 
    res = 900)
p1

dev.off()

p2 <- ggplot(masterresults2, aes(pathway.name, log.fold.change)) +
  geom_point(aes(colour=qval),size=3) +
  scale_colour_gradient(limits=c(0.001, 0.33), n.breaks = 6, low="red", high="black") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, colour = 'black', size = 12.5),
        axis.text.y = element_text(color = 'black', size = 12.5),
        axis.title.x = element_text(color = 'black', size = 12.5),
        axis.title.y = element_text(color = 'black', size = 12.5),
        legend.text=element_text(size=12.5),
        legend.title=element_text(size=12.5), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.margin = unit(c(1,1,1,3), "cm")) +
  ylab(expression('log'[2]*'(Responder/Non Responder)')) + 
  geom_hline(yintercept=0,linetype="dashed",color = "grey") +
  labs(colour='Q value') 


png('qusage_overall_TOC_CDAI50_combined_nocov.png', height = 16, width = 22, units = 'cm', 
    res = 900)
p2

dev.off()
