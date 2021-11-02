rm(list=ls())
library(RNASeqPower)
library(tidyr)
library(dplyr)
library(edgeR)
library(ggplot2)
library(ggsci)

drug.abbreviation <- "toci"  # ritux or toci
drug <- "Tocilizumab"  # Rituximab or Tocilizumab

# load and align count data
# using the metadata Anna used for her DEGs 
metadata <- readRDS("Data/R4RA_metadata.RDS")
metadata <- subset(metadata, Randomised.medication == drug)

counts <- readRDS("Data/txi.R4RA.gene.symbols.RDS")$counts 
tokeep <- row.names(counts[ rowSums((cpm(counts))>1)>=10, ])
counts <- as.data.frame(counts[tokeep, rownames(metadata)])

# # DEGs showed in Fig.2 a and b (no covariate)
DEGs <- read.csv(paste0(drug.abbreviation, ".DEGs.CDAI50.csv"), 
                 stringsAsFactors = FALSE) %>%
  filter(qvalue < 0.05) %>%
  pull(Gene)

# # DEGs showed in Fig.2 c and d (PC1 covariate added)
# # Results are similar to the previous case
# DEGs <- read.csv(paste0("CDAI50.", drug.abbreviation, 
#                         ".switchr.DeSeq2_nonrespref_PC1cov.csv"),
#                              stringsAsFactors = FALSE) %>% 
#   filter(qvalue < 0.05) %>%
#   pull(X)

DEGs <- DEGs[DEGs %in% rownames(counts)]
counts.DEGs <- counts[DEGs,]
AveCounts <- mean(rowMeans(counts.DEGs))

meta.resp <- subset(metadata, metadata$CDAI.response.status.V7 == "Responder" )
counts.resp <- counts[, rownames(meta.resp)]

meta.non.resp <- subset(metadata, metadata$CDAI.response.status.V7 == "Non.Responder" )
counts.non.resp <- counts[, rownames(meta.non.resp)]

dispersion.resp   <- estimateDisp(y=counts.resp, 
                                  design = model.matrix(~as.factor(meta.resp$Pathotype.V2))
)$common.dispersion

dispersion.non.resp <- estimateDisp(y=counts.non.resp, 
                                    design = model.matrix(~as.factor(meta.non.resp$Pathotype.V2))
)$common.dispersion

my.depth <- AveCounts

power.values <- as.data.frame(rnapower(depth = my.depth,
                                       cv = dispersion.resp, 
                                       cv2 = dispersion.non.resp,
                                       effect=c(1.25, 1.5, 1.75, 2),
                                       alpha= .05, power=c(.7, .8, .9)))
power.values$logFC <- rownames(power.values)
power.values.long <- reshape2::melt(power.values, "logFC", 
                                    measures = c("0.7", "0.8", "0.9"),
                                    variable.name = "Power", 
                                    value.name = "No.of.Samples")

sampleSizes <- seq(5, 250, 5)
sample.sizes <- as.data.frame(rnapower(depth=my.depth, cv=dispersion.resp,
                                       cv2 = dispersion.non.resp, 
                                       effect=c(1.25, 1.5, 1.75, 2),
                                       alpha= .05, n = sampleSizes,
                                       n2 = sampleSizes))
sample.sizes$SampleSizes <- as.numeric(rownames(sample.sizes))
sample.sizes.long <- reshape2::melt(sample.sizes, "SampleSizes", 
                                    measures = c("1.25", "1.5", "1.75", "2"),
                                    variable.name = "logFC", 
                                    value.name = "Power")

sample.sizes.long.1.5 <- subset(sample.sizes.long, logFC == 1.5)

powerPlot <- ggplot(power.values.long, aes(logFC, No.of.Samples, group = Power, color = Power)) + 
  geom_point() +
  geom_line() +
  scale_color_npg() +
  theme_bw() +
  theme(axis.text = element_text(family = "TT Arial", colour = "black"),
        axis.title= element_text(family = "TT Arial", colour = "black")) +
  labs(title = paste(" Responder group dispersion:", round(dispersion.resp, digits = 2), "\n",
                     "Non responder group dispersion:", round(dispersion.non.resp, digits = 2), "\n",
                     "Average depth of coverage of DEGs:", round(my.depth)),
       y = "No of Samples", 
       x = "LogFC")

sampleSizePlot <- ggplot(sample.sizes.long, aes(SampleSizes, Power,
                                                group = logFC, color = logFC)) + 
  geom_point() +
  geom_line() +
  scale_color_npg() +
  theme_bw() +
  theme(axis.text  = element_text(family = "TT Arial", colour = "black"),
        axis.title = element_text(family = "TT Arial", colour = "black"),
        plot.caption = element_text(family = "TT Arial", colour = "black", size = 11)) +
  labs(title = paste0("Dispersion of the responder group: ", 
                      round(dispersion.resp, digits = 2), "\n",
                      "Dispersion of the non responder group: ", 
                      round(dispersion.non.resp, digits = 2), "\n",
                      "Average depth of coverage of DEGs: ", round(my.depth)),
       caption = paste0("# of samples needed to reach 80% power (logFC=1.5): ",
                        sample.sizes.long.1.5[which.min(abs(sample.sizes.long.1.5$Power - 0.80)), "SampleSizes"], "\n",
                        "# of samples needed to reach 90% power (logFC=1.5): ",
                        sample.sizes.long.1.5[which.min(abs(sample.sizes.long.1.5$Power - 0.90)), "SampleSizes"]),
       x = "Samples")

sampleSizePlot
powerPlot

svg(paste0(drug.abbreviation, ".sampleSizePlot.firstDEGs.svg"), height = 5.5)
# pdf(paste0(drug.abbreviation, ".sampleSizePlot.svg"), height = 5.5)
sampleSizePlot
dev.off()

svg(paste0(drug.abbreviation, ".powerPlot.firstDEGs.svg"), height = 5.5)
# pdf(paste0(drug.abbreviation, ".powerPlot.svg"), height = 5.5)
powerPlot
dev.off()