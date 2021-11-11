### clustering and heatmap

# R4RA clustering for all and individual treatment groups

library(M3C)
library(edgeR)
library(limma)
library(ComplexHeatmap)
load("./R4RA_V2_data_trimmed.RData")


# variability and expression strength filtering
filtering <- TRUE # set to TRUE for unsupervised clustering
cvx <- 0.075

### main script

## get the data out

data <- txiready$vst

# get the metadata
meta <- read.csv('./meta_data.csv')
meta <- subset(meta, meta$Visit == 3)


setdiff(colnames(data),meta$Seq_ID.V2)

meta <- subset(meta, meta$Seq_ID.V2 %in% colnames(data))
data <- data[,meta$Seq_ID.V2]


data_backup <- data # save for boxplots

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
  }
  data3 <- subset(data2, row.names(data2) %in% names)
  nrow(data3)
}else{
  data3 <- data
}

colnames(meta)[colnames(meta)=='Seq_ID.V2'] <- 'ID'

### 

data3 <- t(scale(t(data3)))

### clustering

set.seed(123)
res <- M3C(data3,method=2,clusteralg='pam',printres=FALSE,repsref=1000,repsreal=1000) 

# reorder data according to any clustering algorithm
ind <- sort(as.vector(res$assignments),index.return=TRUE)
datax <- data3[,ind$ix]
annonx <- meta[ind$ix,]
annonx$consensuscluster <- as.factor(as.character(ind$x))

  ## reorder consensus clusters 
  newv <- as.factor(c(2,1)) 
  annonx <- annonx[order(match(annonx$consensuscluster, newv)),] # reorder clinical
  annonx$consensuscluster <- as.character(1:length(newv))[ match(annonx$consensuscluster, as.character(unique(annonx$consensuscluster)) ) ]
  datax <- datax[,annonx$ID]


datax2 <- datax # we manipulate this later so save a copy


### all heatmap

annonx$DAS28.CRP.EULARresp.V7 <- as.character(annonx$DAS28.CRP.EULARresp.V7)
annonx$DAS28.CRP.EULARresp.V7[annonx$DAS28.CRP.EULARresp.V7 == 'Good.Responder'] <- 'Good'
annonx$DAS28.CRP.EULARresp.V7[annonx$DAS28.CRP.EULARresp.V7 == 'Non.Responder'] <- 'Moderate_or_none'
annonx$DAS28.CRP.EULARresp.V7[annonx$DAS28.CRP.EULARresp.V7 == 'Moderate.Responder'] <- 'Moderate_or_none'

annonx$Cell.Type.V2 <- as.character(annonx$Cell.Type.V2)

cutoffval <- 2.5

datax <- apply(datax, 2, function(x) ifelse(x > cutoffval, cutoffval, x)) # compress data within range
datax <- apply(datax, 2, function(x) ifelse(x < -cutoffval, -cutoffval, x)) # compress data within range


## individual
colnames(annonx)[colnames(annonx)=='CDAI.response.status.V7'] <- 'CDAI.50.improvement'
annonx$Current.Medication <- as.character(annonx$Current.Medication)
r_annonx <- subset(annonx, annonx$Current.Medication == 'Rituximab')
t_annonx <- subset(annonx, annonx$Current.Medication == 'Tocilizumab')

# check sizes of clusters for both treatment groups
# need as similar sizes as possible really
table(r_annonx$consensuscluster)
table(t_annonx$consensuscluster)

## CDAI 50 improvement for consensuscluster

table(r_annonx[c('consensuscluster','CDAI.50.improvement')])
fisher.test(table(r_annonx[c('consensuscluster','CDAI.50.improvement')]))

table(t_annonx[c('consensuscluster','CDAI.50.improvement')])
fisher.test(table(t_annonx[c('consensuscluster','CDAI.50.improvement')]))


### heatmaps individual

library(RColorBrewer)

scale1 <- brewer.pal(6, "YlOrRd")
scale2 <- brewer.pal(6, "YlGnBu")
scale3 <- brewer.pal(6, "RdPu")
scale4 <- brewer.pal(6, "YlGn")
scale5 <- brewer.pal(6, "BuPu")

## RT

r_datax <- datax[,as.character(r_annonx$ID)]

ann_colors2 <- ggsci::pal_futurama()(4) # futurama palette
chrispal1 <- c('cornflowerblue', 'gold1', 'darkorchid', 'skyblue1', 'plum1', 'violetred', 'forestgreen', 'midnightblue')

ha2 = HeatmapAnnotation(df = r_annonx[,c('consensuscluster','Cell.Type.V2' ,'Pathotype.V2','CDAI.50.improvement',
                                         'DAS28.CRP.EULARresp.V7','CD20.V2','CD138.V2',
                                         'CD68L.V2','CD68SL.V2','CD3.V2'),drop=FALSE],
                        simple_anno_size = unit(4, "mm"),
                        annotation_name_gp = gpar(fontsize = 0),
                        col = list(consensuscluster = c("1" = chrispal1[1], "2" = chrispal1[2], "3" = chrispal1[3], 
                                                        "4" = chrispal1[4], "5" = chrispal1[5], "6" = chrispal1[6]),
                                   Cell.Type.V2 = c('Brich' = '#FF6F00FF', 'Bpoor' = '#C71000FF', 'Unknown' = 'grey', 'GC' = 'blue'),
                                   Pathotype.V2 = c("Lymphoid" = 'firebrick1', "Fibroid" = chrispal1[4], "Myeloid" = chrispal1[5], "Ungraded" = 'grey'),
                                   CDAI.50.improvement = c('Responder' = 'mediumblue', 'Non.Responder' = 'black'),
                                   DAS28.CRP.EULARresp.V7 = c('Good' = 'slateblue1', 'Moderate_or_none' = 'black'),
                                   CD3.V2 = c('0'=scale1[2],'1'=scale1[3],'2'=scale1[4],'3'=scale1[5],'4'=scale1[6]),
                                   CD20.V2 = c('0'=scale2[2],'1'=scale2[3],'2'=scale2[4],'3'=scale2[5],'4'=scale2[6]),
                                   CD68L.V2 = c('0'=scale3[2],'1'=scale3[3],'2'=scale3[4],'3'=scale3[5],'4'=scale3[6]),
                                   CD68SL.V2 = c('0'=scale4[2],'1'=scale4[3],'2'=scale4[4],'3'=scale4[5],'4'=scale4[6]),
                                   CD138.V2 = c('0'=scale5[2],'1'=scale5[3],'2'=scale5[4],'3'=scale5[5],'4'=scale5[6])),
                        annotation_legend_param = list(consensuscluster = list(title = "Consensus cluster", 
                                                                               title_gp = gpar(fontface = 'plain', fontsize = '10')),
                                                       Pathotype.V2 = list(title = "Pathotype", 
                                                                           title_gp = gpar(fontface = 'plain', fontsize = '10')),
                                                       Cell.Type.V2 = list(title = "Cell.type", 
                                                       title_gp = gpar(fontface = 'plain', fontsize = '10')),
                                                       CDAI.50.improvement = list(title = "CDAI.50%.improvement", 
                                                                                  title_gp = gpar(fontface = 'plain', fontsize = '10')),
                                                       DAS28.CRP.EULARresp.V7 = list(title = "EULAR3.response", 
                                                                                     title_gp = gpar(fontface = 'plain', fontsize = '10')),
                                                       CD3.V2 = list(title = "CD3.max",
                                                                     title_gp = gpar(fontface = 'plain', fontsize = '10')),
                                                       CD20.V2 = list(title = "CD20.max",
                                                                      title_gp = gpar(fontface = 'plain', fontsize = '10')),
                                                       CD68L.V2 = list(title = "CD68L.max",
                                                                       title_gp = gpar(fontface = 'plain', fontsize = '10')),
                                                       CD68SL.V2 = list(title = "CD68SL.max",
                                                                        title_gp = gpar(fontface = 'plain', fontsize = '10')),
                                                       CD138.V2 = list(title = "CD138.max",
                                                                       title_gp = gpar(fontface = 'plain', fontsize = '10'))))


  h <- Heatmap(r_datax, name = "Norm intensity",
               cluster_columns = FALSE, col = gplots::bluered(256),
               heatmap_legend_param = list(color_bar = "continuous",
                                           grid_border = NA, title_gp = gpar(fontface = 'plain', fontsize = '10'), heatmap_annotation_side = "bottom"),
               row_names_gp = gpar(fontsize = 10), top_annotation = ha2,
               column_names_gp = gpar(fontsize = 0),
               clustering_distance_rows = "pearson",
               clustering_method_rows = "complete",
               show_column_dend = FALSE,
               annotation_legend_side = "bottom",
               #row_dend_width = unit(20, "mm"),
               show_row_names = FALSE,
               show_heatmap_legend = TRUE
  )

h


## TOC

t_datax <- datax[,as.character(t_annonx$ID)]
ann_colors2 <- ggsci::pal_futurama()(4) # futurama palette
chrispal1 <- c('cornflowerblue', 'gold1', 'darkorchid', 'skyblue1', 'plum1', 'violetred', 'forestgreen', 'midnightblue')

colnames(t_annonx)[colnames(t_annonx)=='CDAI.response.status.V7'] <- 'CDAI.50.improvement'

ha2 = HeatmapAnnotation(df = t_annonx[,c('consensuscluster', 'Cell.Type.V2','Pathotype.V2','CDAI.50.improvement',
                                         'DAS28.CRP.EULARresp.V7','CD20.V2','CD138.V2',
                                         'CD68L.V2','CD68SL.V2','CD3.V2'),drop=FALSE],
                        simple_anno_size = unit(4, "mm"),
                        annotation_name_gp = gpar(fontsize = 0),
                        col = list(consensuscluster = c("1" = chrispal1[1], "2" = chrispal1[2], "3" = chrispal1[3], 
                                                        "4" = chrispal1[4], "5" = chrispal1[5], "6" = chrispal1[6]),
                                   Cell.Type.V2 = c('Brich' = '#FF6F00FF', 'Bpoor' = '#C71000FF', 'Unknown' = 'grey', 'GC' = 'blue'),
                                   Pathotype.V2 = c("Lymphoid" = 'firebrick1', "Fibroid" = chrispal1[4], "Myeloid" = chrispal1[5], "Ungraded" = 'grey'),
                                   CDAI.50.improvement = c('Responder' = 'mediumblue', 'Non.Responder' = 'black'),
                                   DAS28.CRP.EULARresp.V7 = c('Good' = 'slateblue1', 'Moderate_or_none' = 'black'),
                                   CD3.V2 = c('0'=scale1[2],'1'=scale1[3],'2'=scale1[4],'3'=scale1[5],'4'=scale1[6]),
                                   CD20.V2 = c('0'=scale2[2],'1'=scale2[3],'2'=scale2[4],'3'=scale2[5],'4'=scale2[6]),
                                   CD68L.V2 = c('0'=scale3[2],'1'=scale3[3],'2'=scale3[4],'3'=scale3[5],'4'=scale3[6]),
                                   CD68SL.V2 = c('0'=scale4[2],'1'=scale4[3],'2'=scale4[4],'3'=scale4[5],'4'=scale4[6]),
                                   CD138.V2 = c('0'=scale5[2],'1'=scale5[3],'2'=scale5[4],'3'=scale5[5],'4'=scale5[6])),
                        annotation_legend_param = list(consensuscluster = list(title = "Consensus cluster", 
                                                                               title_gp = gpar(fontface = 'plain', fontsize = '10')),
                                                       Pathotype.V2 = list(title = "Pathotype", 
                                                                           title_gp = gpar(fontface = 'plain', fontsize = '10')),
                                                       Cell.Type.V2 = list(title = "Cell.type", 
                                                       title_gp = gpar(fontface = 'plain', fontsize = '10')),
                                                       CDAI.50.improvement = list(title = "CDAI.50%.improvement", 
                                                                                  title_gp = gpar(fontface = 'plain', fontsize = '10')),
                                                       DAS28.CRP.EULARresp.V7 = list(title = "EULAR3.response", 
                                                                                     title_gp = gpar(fontface = 'plain', fontsize = '10')),
                                                       CD3.V2 = list(title = "CD3.max",
                                                                     title_gp = gpar(fontface = 'plain', fontsize = '10')),
                                                       CD20.V2 = list(title = "CD20.max",
                                                                      title_gp = gpar(fontface = 'plain', fontsize = '10')),
                                                       CD68L.V2 = list(title = "CD68L.max",
                                                                       title_gp = gpar(fontface = 'plain', fontsize = '10')),
                                                       CD68SL.V2 = list(title = "CD68SL.max",
                                                                        title_gp = gpar(fontface = 'plain', fontsize = '10')),
                                                       CD138.V2 = list(title = "CD138.max",
                                                                       title_gp = gpar(fontface = 'plain', fontsize = '10'))))


  h <- Heatmap(t_datax, name = "Norm intensity",
               cluster_columns = FALSE, col = gplots::bluered(256),
               heatmap_legend_param = list(color_bar = "continuous",
                                           grid_border = NA, title_gp = gpar(fontface = 'plain', fontsize = '10')),
               row_names_gp = gpar(fontsize = 10), top_annotation = ha2,
               column_names_gp = gpar(fontsize = 0),
               clustering_distance_rows = "pearson",
               clustering_method_rows = "complete",
               row_dend_width = unit(20, "mm"),
               show_row_names = FALSE,
               show_heatmap_legend = TRUE)

h

### ritux boxplots

source('./Baseline_clustering_M3C/flexible_boxplots.R')
mdata <- data_backup[,r_annonx$ID]
### genes

# IL6 markers

# B cell markers
library(ggpubr)

gene_boxplots <- list()
for (g in c('MS4A1','CD79A','CD79B','PIK3CA','BTK','SYK')){
  toadd <- g
  mdata2 <- mdata[row.names(mdata)==toadd,]
  r_annonx[[toadd]] <- as.numeric(mdata2[toadd,])
  gene_boxplots[[g]] <- flexibleboxplot(r_annonx,'consensuscluster',toadd,xlab='consensuscluster',
                                        theme='bw',fill=15,title=toadd)+
    stat_compare_means(label = "p.format", label.x = 1.35, label.y.npc = 0.92, method ='kruskal.test', size = 8)
  ## stat test
  toaddy <- paste(toadd,'~',sep='')
  t <- kruskal.test(as.formula(paste(toaddy, paste('consensuscluster', collapse="+"))), data = r_annonx)
  print(t$p.value)
}

### modules

source('./Baseline_clustering_M3C/mergedatausingmodules.R')
mdata <- mergemodules(data=data_backup,type='wgcna',method='zmean')
r_annonx$ID <- gsub('\\-','.',r_annonx$ID)
mdata <- mdata[,r_annonx$ID]
row.names(mdata) <- gsub(' ','_',row.names(mdata)) # collapse all spaces
row.names(mdata) <- gsub(',','',row.names(mdata)) # collapse all spaces
row.names(mdata) <- gsub('\\+','',row.names(mdata))
module_boxplots <- list()

for (g in c('S136_B_cells','S39_M1_macrophage_cytokine_signalling','S115_Fibroblast_2a_THY1+')){
  g <- gsub('\\+','',g) # remove the special character
  toadd <- g
  mdata2 <- mdata[row.names(mdata)==toadd,]
  r_annonx[[toadd]] <- as.numeric(mdata2[toadd,])
  module_boxplots[[g]] <- flexibleboxplot(r_annonx,'consensuscluster',toadd,xlab='consensuscluster',
                                          theme='bw',fill=10)+
    stat_compare_means(label = "p.format", label.x = 1.35, label.y.npc = 0.92, method ='kruskal.test', size = 8)
  ## stat test
  toaddy <- paste(toadd,'~',sep='')
  t <- kruskal.test(as.formula(paste(toaddy, paste('consensuscluster', collapse="+"))), data = r_annonx)
  print(t$p.value)
}

# print boxplots: 1

library(grid)
library(gridExtra)
library(cowplot)

p1 <- gene_boxplots[[1]]+theme(plot.title = element_text(size = 24))
p2 <- gene_boxplots[[2]]+theme(plot.title = element_text(size = 24))
p3 <- gene_boxplots[[3]]+theme(plot.title = element_text(size = 24))

p4 <- gene_boxplots[[4]]+theme(plot.title = element_text(size = 24))
p5 <- gene_boxplots[[5]]+theme(plot.title = element_text(size = 24))
p6 <- gene_boxplots[[6]]+theme(plot.title = element_text(size = 24))

p7 <- module_boxplots[[1]]+
  labs(title = 'S136 B cells' ) +
  theme( plot.title = element_text(size = 24)) 
p8 <- module_boxplots[[2]] +
  labs(title = 'S39 M1 macrophage \ncytokine signalling' )+
  theme( plot.title = element_text(size = 24))
p9 <- module_boxplots[[3]]+
  labs(title = 'S115 \nFibroblast 2a THY1+' )+
  theme(plot.title = element_text(size = 24))

#setwd("~/Desktop/")
png('all_ritux_genes_modules_unsupervised_clusters_trimmed_K2_V2.png', height = 26, width = 26, units = 'cm', 
    res = 900, typ = 'cairo')
plot_grid(p1 + theme(axis.text.x = element_blank(),
                     axis.ticks.x = element_blank(),
                     axis.title.x = element_blank(),
                     axis.title.y = element_blank(),
                     legend.position = 'none'),
          p2+ theme(axis.text.x = element_blank(),
                    axis.ticks.x = element_blank(),
                    axis.title.x = element_blank(),
                    axis.title.y = element_blank(),
                    legend.position = 'none'),
          p3+ theme(axis.text.x = element_blank(),
                    axis.ticks.x = element_blank(),
                    axis.title.x = element_blank(),
                    axis.title.y = element_blank(), legend.position = 'none'),#legend.text = element_text(size=15), legend.key.width = unit(1.2, 'cm'), legend.key.height = unit(1.2, 'cm') ),
          NULL,NULL, NULL,
          p4+ theme(axis.text.x = element_blank(),
                    axis.ticks.x = element_blank(),
                    axis.title.x = element_blank(),
                    axis.title.y = element_blank(),
                    legend.position = 'none'),
          p5+theme(axis.text.x = element_blank(),
                   axis.ticks.x = element_blank(),
                   axis.title.x = element_blank(),
                   axis.title.y = element_blank(),
                   legend.position = 'none'),
          p6+ theme(axis.text.x = element_blank(),
                    axis.ticks.x = element_blank(),
                    axis.title.x = element_blank(),
                    axis.title.y = element_blank(), legend.position = 'none'), 
          NULL,NULL, NULL,
          p7 + theme(legend.position = 'none',
                     axis.title.y = element_blank()),
          p8+ theme(axis.title.y = element_blank(),
                    legend.position = 'none'),
          p9+ theme(axis.title.y = element_blank(), legend.position = 'none'), 
          align='h', rel_widths = c(0.85,0.85,0.85), rel_heights = c(1,-0.13,1,-0.11,1),
          ncol = 3, label_size = 24)

dev.off()

### tocil boxplots
t_annonx$ID <- gsub('\\.','-',t_annonx$ID)
mdata <- data_backup[,t_annonx$ID]

### genes

# IL6 markers

gene_boxplots <- list()
for (g in c('IL6R','IL6','IL6ST','JAK1','JAK2','STAT3')){
  toadd <- g
  mdata2 <- mdata[row.names(mdata)==toadd,]
  t_annonx[[toadd]] <- as.numeric(mdata2[toadd,])
  gene_boxplots[[g]] <- flexibleboxplot(t_annonx,'consensuscluster',toadd,xlab='consensuscluster',
                                        theme='bw',fill=15,title=toadd)+
    stat_compare_means(label = "p.format", label.x = 1.35, label.y.npc = 0.92, method ='kruskal.test', size = 8)
  ## stat test
  toaddy <- paste(toadd,'~',sep='')
  t <- kruskal.test(as.formula(paste(toaddy, paste('consensuscluster', collapse="+"))), data = t_annonx)
  print(t$p.value)
}

### modules

mdata <- mergemodules(data=data_backup,type='wgcna',method='zmean')
t_annonx$ID <- gsub('\\-','.',t_annonx$ID)
mdata <- mdata[,t_annonx$ID]
row.names(mdata) <- gsub(' ','_',row.names(mdata)) # collapse all spaces
row.names(mdata) <- gsub(',','',row.names(mdata)) # collapse all spaces
row.names(mdata) <- gsub('\\+','',row.names(mdata))
module_boxplots <- list()

for (g in c('S136_B_cells','S39_M1_macrophage_cytokine_signalling','S115_Fibroblast_2a_THY1+')){
  g <- gsub('\\+','',g) # remove the special character
  toadd <- g
  mdata2 <- mdata[row.names(mdata)==toadd,]
  t_annonx[[toadd]] <- as.numeric(mdata2[toadd,])
  module_boxplots[[g]] <- flexibleboxplot(t_annonx,'consensuscluster',toadd,xlab='consensuscluster',
                                          theme='bw',fill=10)+
    stat_compare_means(label = "p.format", label.x = 1.35, label.y.npc = 0.92, method ='kruskal.test', size = 8)
  ## stat test
  toaddy <- paste(toadd,'~',sep='')
  t <- kruskal.test(as.formula(paste(toaddy, paste('consensuscluster', collapse="+"))), data = t_annonx)
  print(t$p.value)
}

# print boxplots: 1

p1 <- gene_boxplots[[1]]+theme(plot.title = element_text(size = 24))
p2 <- gene_boxplots[[2]]+theme(plot.title = element_text(size = 24))
p3 <- gene_boxplots[[3]]+theme(plot.title = element_text(size = 24))

p4 <- gene_boxplots[[4]]+theme(plot.title = element_text(size = 24))
p5 <- gene_boxplots[[5]]+theme(plot.title = element_text(size = 24))
p6 <- gene_boxplots[[6]]+theme(plot.title = element_text(size = 24))

p7 <- module_boxplots[[1]]+
  labs(title = 'S136 B cells' ) +
  theme( plot.title = element_text(size = 24)) 
p8 <- module_boxplots[[2]] +
  labs(title = 'S39 M1 macrophage \ncytokine signalling' )+
  theme( plot.title = element_text(size = 24))
p9 <- module_boxplots[[3]]+
  labs(title = 'S115 \nFibroblast 2a THY1+' )+
  theme(plot.title = element_text(size = 24))

plot_grid(p1 + theme(axis.text.x = element_blank(),
                     axis.ticks.x = element_blank(),
                     axis.title.x = element_blank(),
                     axis.title.y = element_blank(),
                     legend.position = 'none'),
          p2+ theme(axis.text.x = element_blank(),
                    axis.ticks.x = element_blank(),
                    axis.title.x = element_blank(),
                    axis.title.y = element_blank(),
                    legend.position = 'none'),
          p3+ theme(axis.text.x = element_blank(),
                    axis.ticks.x = element_blank(),
                    axis.title.x = element_blank(),
                    axis.title.y = element_blank(), legend.position = 'none'),#legend.text = element_text(size=15), legend.key.width = unit(1.2, 'cm'), legend.key.height = unit(1.2, 'cm') ),
          NULL,NULL, NULL,
          p4+ theme(axis.text.x = element_blank(),
                    axis.ticks.x = element_blank(),
                    axis.title.x = element_blank(),
                    axis.title.y = element_blank(),
                    legend.position = 'none'),
          p5+theme(axis.text.x = element_blank(),
                   axis.ticks.x = element_blank(),
                   axis.title.x = element_blank(),
                   axis.title.y = element_blank(),
                   legend.position = 'none'),
          p6+ theme(axis.text.x = element_blank(),
                    axis.ticks.x = element_blank(),
                    axis.title.x = element_blank(),
                    axis.title.y = element_blank(), legend.position = 'none'), 
          NULL,NULL, NULL,
          p7 + theme(legend.position = 'none',
                     axis.title.y = element_blank()),
          p8+ theme(axis.title.y = element_blank(),
                    legend.position = 'none'),
          p9+ theme(axis.title.y = element_blank(), legend.position = 'none'), 
          align='h', rel_widths = c(0.85,0.85,0.85), rel_heights = c(1,-0.13,1,-0.11,1),
          ncol = 3, label_size = 24)


## histology for all samples - no subsets
## boxplots
source('./Baseline_clustering_M3C/histology_boxplots.R')
histology_boxplots <- list()
for (g in c('CD3.V2','CD20.V2','CD68L.V2','CD68SL.V2','CD138.V2', 'CD79a')){
  toadd <- g
  
  histology_boxplots[[g]] <- histoboxplot(annonx,'consensuscluster',toadd,xlab='consensuscluster',
                                          theme='bw',fill=15)+
    stat_compare_means(label = "p.format", label.x = 1.35, label.y = 3.5, method ='kruskal.test', size = 8)
  ## stat test
  toaddy <- paste(toadd,'~',sep='')
  t <- kruskal.test(as.formula(paste(toaddy, paste('consensuscluster', collapse="+"))), data = annonx)
  print(t$p.value)
}

p1 <- histology_boxplots[[1]]+labs(title = 'CD3')+ theme(plot.title = element_text(size = 26))
p2 <- histology_boxplots[[2]]+labs(title='CD20') +theme(plot.title = element_text(size = 26))
p3 <- histology_boxplots[[3]]+labs(title='CD68L') +theme(plot.title = element_text(size = 26))
p4 <- histology_boxplots[[4]]+labs(title='CD68SL') +theme(plot.title = element_text(size = 26))
p5 <- histology_boxplots[[5]]+labs(title='CD138') +theme(plot.title = element_text(size = 26))
p6 <- histology_boxplots[[6]]+labs(title='CD79a') +theme(plot.title = element_text(size = 26))

a <- plot_grid(p1 + theme(axis.text.x = element_blank(),
                          axis.ticks.x = element_blank(),
                          axis.title.x = element_blank(),
                          axis.title.y = element_blank(),
                          axis.text.y = element_text(size= 20),
                          legend.position = 'none'),
               p2+ theme(axis.text.x = element_blank(),
                         axis.ticks.x = element_blank(),
                         axis.title.x = element_blank(),
                         axis.title.y = element_blank(),
                         axis.text.y = element_blank(),
                         legend.position = 'none'),
               p3+ theme(axis.text.x = element_blank(),
                         axis.ticks.x = element_blank(),
                         axis.title.x = element_blank(),
                         axis.title.y = element_blank(),
                         axis.text.y = element_blank(),
                         legend.position = 'none'),
               p4 + theme(axis.text.x = element_blank(),
                          axis.ticks.x = element_blank(),
                          axis.title.x = element_blank(),
                          axis.title.y = element_blank(),
                          axis.text.y = element_blank(),
                          legend.position = 'none'),
               p5 + theme(axis.text.x = element_blank(),
                          axis.ticks.x = element_blank(),
                          axis.title.x = element_blank(),
                          axis.title.y = element_blank(),
                          axis.text.y = element_blank(),
                          legend.position = 'none'),
               p6+ theme(axis.text.x = element_blank(),
                         axis.ticks.x = element_blank(),
                         axis.title.x = element_blank(),
                         axis.title.y = element_blank(),
                         axis.text.y = element_blank(),
                         legend.position = 'none'),
               labels = NA,
               ncol = 6, label_size = 20)