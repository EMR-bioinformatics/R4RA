rm(list = ls())
library(DESeq2)
library(qvalue)
library(volcano3D)
library(dplyr)
library(plotly)
library(writexl)

# load data
metadata <- readRDS("Data/R4RA.crossover.metadata.RDS")
metadata <- subset(metadata, !is.na(metadata$best.trt.V2))
txi      <- readRDS("txi.R4RA.gene.symbols.RDS")
vst      <- readRDS("txi.R4RA.vst.gene.symbols.RDS")

# NOTE: column "best.trt.V2" in the metadata is the column to be used to 
# distinguish patients in Pro-RTX, Pro-TOC and Refractory as explained in 
# Figure 2 g 

# Restrincting the analysis to protein coding genes and aligning data 
coding.genes <- read.csv("gene_type.txt", stringsAsFactors = FALSE)[, c(1,2)]
coding.genes <- subset(coding.genes, coding.genes$Gene.type == "protein_coding")
txi$abundance <- txi$abundance[which(rownames(txi$abundance) %in%
                                     coding.genes$Gene.name), metadata$RNAseq.ID]
txi$counts <- txi$counts[which(rownames(txi$counts) %in% 
                                 coding.genes$Gene.name), metadata$RNAseq.ID]
txi$length <- txi$length[which(rownames(txi$length) %in%
                                    coding.genes$Gene.name), metadata$RNAseq.ID]

# Compute likelihood ratio test on the three groups
dds = DESeqDataSetFromTximport(txi, metadata, ~best.trt.V2)
dds_DE <- DESeq(dds)
dds_LRT <- DESeq(dds, test = "LRT", reduced = ~1, parallel = TRUE)

LRT <- results(dds_LRT, parallel=TRUE)
LRT$qvalue <- qvalue(p = LRT$pvalue)$qvalues
LRT <- LRT[match(rownames(txi$counts), rownames(LRT)), 
           c("pvalue", "padj", "qvalue", "log2FoldChange")]

# Compute pairwise DEGs 
groups <- levels(metadata$best.trt.V2)
contrasts <- c(paste(groups[1], groups[2], sep="-"),
               paste(groups[2], groups[3], sep="-"),
               paste(groups[3], groups[1], sep="-"))

Pvals_DESeq_DE <- lapply(contrasts, function(x){
  vars <- unlist(strsplit(x, split="-"))
  out <- results(dds_DE, contrast=c("best.trt.V2", vars))
  out$qvalue <- qvalue(p = out$pvalue)$qvalues
  out <- out[,c("pvalue", "padj", "qvalue", "log2FoldChange")]
  out <- out[match(rownames(txi$counts), rownames(out)), ]
})


# Combine pairwise + LRT results 
my_pvalues <- cbind(Pvals_DESeq_DE[[1]],
                     Pvals_DESeq_DE[[2]],
                     Pvals_DESeq_DE[[3]],
                     LRT)
my_pvalues <- as.data.frame(my_pvalues[complete.cases(my_pvalues), ])


# Calculating the polar coordinates needed for the following plots 
colnames(metadata)[which(colnames(metadata) == "RNAseq.ID")] <- "ID"
vst <- vst[rownames(my_pvalues), metadata$ID]

polar <- polar_coords(sampledata = metadata,
                      contrast = "best.trt.V2",
                      pvalues = my_pvalues,
                      expression = vst,
                      p_col_suffix = "pvalue",
                      padj_col_suffix = "qvalue",
                      fc_col_suffix = "log2FoldChange",
                      multi_group_prefix = "LRT",
                      non_sig_name = "Not Significant",
                      significance_cutoff = 0.05,
                      fc_cutoff = 0)

# Mark mixed groups as "NS" (Not Specific)
polar@polar$sig <- as.character(polar@polar$sig)
polar@polar$sig[which(polar@polar$sig != "Refractory+" &
                      polar@polar$sig != "Pro-RTX+" &
                      polar@polar$sig != "Pro-TOC+" &
                      polar@polar$sig != "Pro-RTX+Pro-TOC+"  )] <- "NS"
polar@polar$sig <- as.factor(polar@polar$sig)
polar@non_sig_name <- "NS"
# NOTE: NS will stand for either "Not Specific" or "Not Significant" at all


# Producuing radial and volcano plots as in Figure 2 i and j
labels_polar <- read.csv("labels.polar.plot.genes.csv", stringsAsFactors = F, header = F)[,1]
radial_plotly(polar,
                 plot_height = 720, 
                 plot_width = 720,
                 axis_title_size = 20,
                 label_size = 12,
                 label_rows = labels_polar,
                 colours = c("red", "gray60", "blue", "green3", "gold3", "gray60"),
                 colour_code_labels = FALSE, 
                 label_colour = 'black',
                 axis_title_offset = 100,
                 axis_thickness = 3,
                 axes_on_origin = TRUE,
                 marker_outline_width = 0) %>%
  config(editable=T, doubleClick = "reset",
         edits = list(annotationPosition = F, annotationText = FALSE),
         toImageButtonOptions=list(format="svg"))


labels_volcano <- read.csv("labels.volcano3D.genes.csv", stringsAsFactors = F, header = F)[,1]
volcano3D(polar,
             plot_height = 4000,
             plot_width = 4000,
             colours = c("red", "cyan", "blue", "green3", "gold3", "gold2"),
             label_size = 23,
             axis_size = 50,
             tick_size = 30,
             legend_size = 0,
             axis_width = 7,
             arrow_width = 1.5, 
             grid_thickness = 4,
             marker_size = 4.5,
             label_rows = labels_volcano,
             fc_or_zscore = 'zscore',
             z_aspectratio = 0.5,
             colour_code_labels = FALSE,
             label_colour = 'black',
             axis_title_offset = 1.2,
             marker_alpha = 0.7,
             marker_outline_width = 0) %>% 
  config(editable=T, toImageButtonOptions=list(format="svg"))


# Supplementary table
table(polar@polar$sig)
refractory <- my_pvalues[polar@polar[which(polar@polar$sig == "Refractory+"), "Name"],
                         -grep("padj", colnames(my_pvalues))]
Pro.RTX <- my_pvalues[polar@polar[which(polar@polar$sig == "Pro-RTX+"), "Name"], 
                      -grep("padj", colnames(my_pvalues))]
Pro.TOC <- my_pvalues[polar@polar[which(polar@polar$sig == "Pro-TOC+"), "Name"], 
                      -grep("padj", colnames(my_pvalues))]
Pro.RTX.Pro.TOC <- my_pvalues[polar@polar[which(polar@polar$sig == "Pro-RTX+Pro-TOC+"), "Name"],
                              -grep("padj", colnames(my_pvalues))]

refractory$Gene      <- rownames(refractory     )    
Pro.RTX$Gene         <- rownames(Pro.RTX        )
Pro.TOC$Gene         <- rownames(Pro.TOC        )
Pro.RTX.Pro.TOC$Gene <- rownames(Pro.RTX.Pro.TOC)

# move the last created column to the first position
refractory      <- refractory     [,c(ncol(refractory     ),1:(ncol(refractory     )-1))]
Pro.RTX         <- Pro.RTX        [,c(ncol(Pro.RTX        ),1:(ncol(Pro.RTX        )-1))]
Pro.TOC         <- Pro.TOC        [,c(ncol(Pro.TOC        ),1:(ncol(Pro.TOC        )-1))]
Pro.RTX.Pro.TOC <- Pro.RTX.Pro.TOC[,c(ncol(Pro.RTX.Pro.TOC),1:(ncol(Pro.RTX.Pro.TOC)-1))]

# export as excel spredsheet with multiple sheets 
volcano.sig.genes <- list(refractory, Pro.RTX, Pro.TOC, Pro.RTX.Pro.TOC)
names(volcano.sig.genes) <- c("a.Refractory", "b.Pro-RTX", "c.Pro-TOC",
                              "d.Pro-RTX and Pro-TOC")
names(volcano.sig.genes) <- substr(names(volcano.sig.genes), start = 1, stop = 29)
write_xlsx(volcano.sig.genes, "Supplementary_Volcano3D_Genes.xlsx")
