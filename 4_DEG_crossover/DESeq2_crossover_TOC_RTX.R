## Differential gene expression script for switchover patients
# load packages
library(ggplot2)
library(limma)
library(edgeR)
library(DESeq2)
library(qusage)
library(stringr)
library(qvalue)
library(dplyr)

txi <- readRDS("txi.RDS")

#metadata


meta <- read.csv('./meta_data.csv') ## with overall response

coldata <- meta %>% 
  filter(!is.na(Seq_ID.V2))

rownames(coldata) <- coldata$Seq_ID.V2

## calculate PC to add as covariate
coldata <-right_join(coldata, PC1, by='Seq_ID.V2')

##filter low expressed genes
tokeep <- row.names(txi$counts[ rowSums((cpm(txi$counts))>1)>=10, ])

##DeSeq2

r_meta <- coldata %>% filter(
  overall_response_CDAI50 %in% c("Neither","RNR_NoInfo","RNR_TR","RR","TNR_RR","Both") ## RR = Rituximab responder, RNR = Rituximab non responder, both =  rituximab and toci responder
)


t_meta <-  coldata %>% filter(
  overall_response_CDAI50 %in% c("Neither","TNR_NoInfo","RNR_TR","TR","TNR_RR","Both") ## TR = Toci responder, TNR = Toci non responder
)

#defining responders/non-responders

r_meta$overall_response_CDAI50 <- as.character(r_meta$overall_response_CDAI50)
r_meta$overall_response_CDAI50[r_meta$overall_response_CDAI50 == 'Both'] <- 'Responder'
r_meta$overall_response_CDAI50[r_meta$overall_response_CDAI50 == 'RR'] <- 'Responder'
r_meta$overall_response_CDAI50[r_meta$overall_response_CDAI50 == 'TNR_RR'] <- 'Responder'
r_meta$overall_response_CDAI50[r_meta$overall_response_CDAI50 == 'RNR_TR'] <- 'Non.Responder'
r_meta$overall_response_CDAI50[r_meta$overall_response_CDAI50 == 'RNR_NoInfo'] <- 'Non.Responder'
r_meta$overall_response_CDAI50[r_meta$overall_response_CDAI50 == 'Neither'] <- 'Non.Responder'
r_meta$overall_response_CDAI50 <- as.factor(r_meta$overall_response_CDAI50)

t_meta$overall_response_CDAI50 <- as.character(t_meta$overall_response_CDAI50)
t_meta$overall_response_CDAI50[t_meta$overall_response_CDAI50 == 'Both'] <- 'Responder'
t_meta$overall_response_CDAI50[t_meta$overall_response_CDAI50 == 'RNR_TR'] <- 'Responder'
t_meta$overall_response_CDAI50[t_meta$overall_response_CDAI50 == 'TR'] <- 'Responder'
t_meta$overall_response_CDAI50[t_meta$overall_response_CDAI50 == 'TNR_NoInfo'] <- 'Non.Responder'
t_meta$overall_response_CDAI50[t_meta$overall_response_CDAI50 == 'TNR_RR'] <- 'Non.Responder'
t_meta$overall_response_CDAI50[t_meta$overall_response_CDAI50 == 'Neither'] <- 'Non.Responder'
t_meta$overall_response_CDAI50 <- as.factor(t_meta$overall_response_CDAI50)

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
#covariate PC1
dds_t <- DESeqDataSetFromTximport(t_txi, colData = t_meta, design= ~PC1 + overall_response_CDAI50) 
dds_r <- DESeqDataSetFromTximport(r_txi, colData = r_meta, design= ~PC1 + overall_response_CDAI50)

#no covariate
dds_t_nocov <- DESeqDataSetFromTximport(t_txi, colData = t_meta, design= ~ overall_response_CDAI50)
dds_r_nocov <- DESeqDataSetFromTximport(r_txi, colData = r_meta, design= ~ overall_response_CDAI50)

dds_t <- DESeq(dds_t) 
dds_r <- DESeq(dds_r)
dds_t_nocov <- DESeq(dds_t_nocov) ## save dds for modular analysis
dds_r_nocov <- DESeq(dds_r_nocov) ## save dds for modular analysis


### insert dds_x_nocov for analysis without covariates
res_t <- results(dds_t, contrast = c("overall_response_CDAI50", "Responder", "Non.Responder"))  # only two levels, no need to specify a contrast. The second level is the reference (Non.Responder)
res_t$qvalue <- NA
res_t$qvalue[!is.na(res_t[, "padj"])] <- qvalue(res_t[!is.na(res_t[, "padj"]), "pvalue"])$qvalues
res_t <- as.data.frame(res_t)

res_r <- results(dds_r, contrast = c("overall_response_CDAI50", "Responder", "Non.Responder"))  # only two levels, no need to specify a contrast. The second level is the reference (Non.Responder)
res_r$qvalue <- NA
res_r$qvalue[!is.na(res_r[, "padj"])] <- qvalue(res_r[!is.na(res_r[, "padj"]), "pvalue"])$qvalues
res_r <- as.data.frame(res_r)

mydf <- res_t # or res_r

## plotting code

# make volcanoplot

mydf$logq <- -log10(mydf$qvalue)
mydf$col <- cut(mydf$logq,breaks = c(-Inf, 1.30103, Inf),labels = c("Insignificant", "Significant"))
res4 <- mydf
#colnames(res4)[1] <- 'BaseMean'

labelsdf <- data.frame('gene'= genes_of_interest) 
res4$col <- as.character(res4$col)
res4$col[res4$gene%in%labelsdf$gene] <- 'red'

p2<-ggplot(res4) +
  geom_point(aes(x = log2FoldChange, y = logq, color = col)) +
  geom_point(aes(x = log2FoldChange, y = logq, color = col), 
             data = subset(res4, col == 'red')) +
  theme_bw() +
  theme(axis.text.x = element_text(color = 'black', size = 12),
        axis.text.y = element_text(color = 'black', size = 12),
        panel.background = element_rect(colour = "white", fill = "white"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.margin = unit(c(1,1,1,1), "cm"),
        legend.position="none") +
  scale_colour_manual(values = c("Insignificant"= "grey", "Significant"="sky blue", "red"="black")) +
  ylab(expression('-log'[10]*'q')) +
  xlab(expression('-log'[2]*'(Non.Responder/Responder)')) +
  scale_x_continuous(limits=c(-8,5)) + 
  scale_y_continuous(limits=c(-0.1,14)) +
  geom_text(data=subset(res4, row.names(res4) %in% labelsdf$gene),color='black',
            size=3.5,fontface='italic',
            aes(log2FoldChange, logq, label=row.names(subset(res4, row.names(res4)%in%labelsdf$gene))),
            hjust=0, vjust=0)

library(ggrepel)

p3 <- p2 + geom_text_repel(data=subset(res4, row.names(res4) %in% labelsdf$gene),color='black',
                           size=3.1,fontface='italic',
                           aes(log2FoldChange, logq, 
                               label=row.names(subset(res4, row.names(res4)%in%labelsdf$gene))))

p3