setwd("/media/gcpeac/Katriona/R4RA/MachineLearning")
library(caret)
library(edgeR)
library(ggplot2)
library(pROC)
library(plotROC)
library(dplyr)
library(M3C)
library(data.table)
library(ggrepel)
library(magrittr)
library(ggpubr)
library(gbm)
library(parallel)
library(sjPlot)
library(jpeg)
library(grid)
set.seed(5217)

####################################
# Data Curation
###################################

# Load in the outlier data
outliers = read.csv2("/media/gcpeac/Chris/R4RA_data_for_sharing/outliers_trimmed_1.csv", sep=",")
outliers = outliers[grepl("Baseline", outliers$ID), ]
outliers = unlist(lapply(strsplit(gsub("-Baseline", "", outliers$ID), split="-"), 
                         function(x){paste(x[(length(x)-2):length(x)], collapse="-")}))


# Reformat IDs
source('/media/gcpeac/Chris_Desktop/R4RA/SCRIPTS/reformat_IDs_functions.R')
load("/media/gcpeac/Chris/R4RA_data_for_sharing/R4RA_V2_data_trimmed.RData")
## load vst data and reformat
vst <- txiready$vst
vst <- fix_colnames(vst) 
vst = vst[, grepl("Baseline", colnames(vst))]
colnames(vst) = gsub("\\-Baseline", "", colnames(vst))
vst = vst[, ! colnames(vst) %in%  outliers]

# Tidy metadata
metadata = readRDS("/media/gcpeac/Katriona/R4RA/full_meta_data.rds")
metadata = metadata[as.character(metadata$seqID) %in% colnames(vst), ]
metadata = metadata[metadata$Visit == 3, ]
cat("Expression samples missing metadata info:", 
    paste0(paste0(colnames(vst)[!colnames(vst) %in% metadata$seqID], sep=""), collaspe=", "))


vst = vst[, unique(as.character(metadata$seqID))]
metadata = metadata[! metadata$seqID %in% outliers, ]


# Load in Giovanni's switched data
load("/media/gcpeac/Katriona/R4RA/MachineLearning/Data/UpdatedR4RASwitchData_250820.RData")
clinical = r4ra.switch.rna
clinical = clinical[as.character(clinical$Patient.I.D.) %in% gsub("\\.", "", metadata$Patient_ID), ]

resp_metric = gsub("\\.Before", "", colnames(clinical)[grepl("Before", colnames(clinical))])
resp_metric = resp_metric[! grepl("rem", resp_metric)]

clinical$DAS28.CRP.EULAR.bin2.After[is.na(clinical$DAS28.CRP.EULAR.bin2.After) & 
                                      clinical$DAS28.CRP.EULAR.bin.After == "Non.Responder" & 
                                      ! is.na(clinical$DAS28.CRP.EULAR.bin.After)] = "Non.Responder" 
clinical$DAS28.ESR.EULAR.bin2.After[is.na(clinical$DAS28.ESR.EULAR.bin2.After) & 
                                      clinical$DAS28.ESR.EULAR.bin.After == "Non.Responder" & 
                                      ! is.na(clinical$DAS28.ESR.EULAR.bin.After)] = "Non.Responder"


# Load in gene information for gencode
load("/media/gcpeac/Myles/R4RA mixed model/gencode.v29.genelength.RData")
gencode_filter = read.csv("/media/gcpeac/Katriona/R4RA/MachineLearning/Data/gencodefilter.csv") 

genetype <- genelength$gene_type[match(rownames(vst), genelength$geneID)]
names(genetype) <- rownames(vst)
genegroup <- gencode_filter$Analysis[match(genetype, gencode_filter$Gene_type)]
names(genegroup) <- rownames(vst)

# keep only protein coding
genelist <- names(genegroup)[genegroup =='Protein' & !is.na(genegroup)]
table(genegroup, useNA='always')
vst = vst[names(genegroup)[genegroup == "Protein" & ! is.na(genegroup)], ]



# Break down overall response to each treatment 
# (RR=ritux response, RNR=ritux non-resp, TR=toci response, TNR=toci non-resp, )
########################################
endpoint_counts = data.frame("RR"=c(), "RNR"=c(), "TR"=c(), "TNR"=c(),"NoInfo"=c(), "N/A"=c())
initial_counts = data.frame("RR"=c(), "RNR"=c(), "TR"=c(), "TNR"=c(),"NoInfo"=c(), "N/A"=c())
overall_counts = data.frame()
toc_resp = data.frame("0"=c(), "1"=c(), "N/A"=c())
rit_resp = data.frame("0"=c(), "1"=c(), "N/A"=c())
any_resp = data.frame("0"=c(), "1"=c(), "N/A"=c())

# for each response metric
for(i in resp_metric){
  
  # resp_ionse pre-switching
  resp_i <- rep(NA, nrow(clinical))
  resp_i[clinical$Randomized.Medication == "Tocilizumab" & clinical[, paste0(i, ".Before")] == "Responder"] <- "TR"
  resp_i[clinical$Randomized.Medication == "Rituximab" & clinical[, paste0(i, ".Before")] == "Responder"] <- "RR"
  resp_i[clinical$Randomized.Medication == "Tocilizumab" & clinical[, paste0(i, ".Before")] == "Non.Responder"] <- "TNR"
  resp_i[clinical$Randomized.Medication == "Rituximab" & clinical[, paste0(i, ".Before")] == "Non.Responder"] <- "RNR"
  resp_i[is.na(clinical[, paste0(i, ".Before")])] <- "NoInfo"
  resp_i <- as.factor(resp_i)
  
  out = t(data.frame(as.numeric(table(resp_i, useNA="always"))))
  cn <- names(table(resp_i, useNA="always"))
  cn[is.na(cn)] <- "N/A"
  dimnames(out) = list(i, cn)
  
  initial_counts = rbind(initial_counts, out[, match(c("RR", "RNR", "TR", "TNR","NoInfo", "N/A"), colnames(out)) ])
  
  clinical[, paste0("initial_response_", i)] <- resp_i
  
  # response post-switch
  resp_f <- rep(NA, nrow(clinical))
  resp_f[clinical$swtrt == "TOCtoRIT" & clinical[, paste0(i, ".After")] == "Responder"] <- "RR"
  resp_f[clinical$swtrt == "RITtoTOC" & clinical[, paste0(i, ".After")] == "Responder"] <- "TR"
  resp_f[clinical$swtrt == "TOCtoRIT" & clinical[, paste0(i, ".After")] == "Non.Responder"] <- "RNR"
  resp_f[clinical$swtrt == "RITtoTOC" & clinical[, paste0(i, ".After")] == "Non.Responder"] <- "TNR"
  resp_f[is.na(clinical[, paste0(i, ".After")])] <- "NoInfo"
  resp_f[which(clinical$swtrt == "RIT" & clinical[, paste0(i, ".Before")] == "Responder")] <- "RR" 
  resp_f[which(clinical$swtrt == "TOC" & clinical[, paste0(i, ".Before")] == "Responder")] <- "TR" 
  resp_f[clinical$resp.path.any == "Early Drop-out"] = "Dropped out" 
  resp_f <- as.factor(resp_f)
  
  out = t(data.frame(as.numeric(table(resp_f, useNA="always"))))
  cn <- names(table(resp_f, useNA="always"))
  cn[is.na(cn)] <- "N/A"
  dimnames(out) = list(i, cn)
  
  endpoint_counts = rbind(endpoint_counts, 
                          out[, match(c("RR", "RNR", "TR", "TNR","NoInfo", "N/A", "Dropped out"), 
                                      colnames(out)) ])
  
  clinical[, paste0("endpoint_response_", i)] <- resp_f
  
  overall_resp <- paste0(clinical[, paste0("initial_response_", i)], "_", 
                         clinical[, paste0("endpoint_response_", i)])
  overall_resp[overall_resp == "RR_RR"] <- "RR"
  overall_resp[overall_resp == "TR_TR"] <- "TR"
  overall_resp[overall_resp == "TNR_RNR"] <- "Neither"
  overall_resp[overall_resp == "RNR_TNR"] <- "Neither"
  overall_resp[overall_resp == "TR_RR"] <- "Both"
  overall_resp[overall_resp == "RR_TR"] <- "Both"
  out = t(data.frame(as.numeric(table(overall_resp, useNA="always"))))
  cn <- names(table(overall_resp, useNA="always"))
  cn[is.na(cn)] <- "N/A"
  dimnames(out) = list(i, cn)
  
  overall_counts = rbind(overall_counts, out[, match(c("Both", "RR", "TNR_RR", "TR", "RNR_TR",
                                                       "TNR_NoInfo", "RNR_NoInfo", 
                                                       "RNR_Dropped out", "TNR_Dropped out", 
                                                       "Neither", "N/A"), colnames(out)) ])
  
  clinical[, paste0("overall_response_", i)] <- overall_resp
  
  tr = rep(NA, nrow(clinical))
  tr[clinical[, paste0("overall_response_", i)] %in% c("RNR_TR", "TR", "Both")] = 1
  tr[clinical[, paste0("overall_response_", i)] %in% c("TNR_RR", "TNR_Dropped out", "TNR_NoInfo", "Neither")] = 0
  n = t(data.frame(as.numeric(table(tr, useNA="always"))))
  cn <- names(table(tr, useNA="always"))
  cn[is.na(cn)] <- "N/A"
  dimnames(n) = list(i, cn)
  toc_resp = rbind(toc_resp, n[, match(c("0", "1", "N/A"), colnames(n)) ])
  
  rr = rep(NA, nrow(clinical))
  rr[clinical[, paste0("overall_response_", i)] %in% c("TNR_RR", "RR", "Both")] = 1
  rr[clinical[, paste0("overall_response_", i)] %in% c("RNR_TR", "RNR_Dropped out", "RNR_NoInfo", "Neither")] = 0
  n = t(data.frame(as.numeric(table(rr, useNA="always"))))
  cn <- names(table(rr, useNA="always"))
  cn[is.na(cn)] <- "N/A"
  dimnames(n) = list(i, cn)
  rit_resp = rbind(rit_resp, n[, match(c("0", "1", "N/A"), colnames(n)) ])
  
  ar = rep(NA, nrow(clinical))
  ar[clinical[, paste0("overall_response_", i)] %in% c("TNR_RR", "RR", "TR", "Both", "RNR_TR", 
                                                       "RR_TNR", "TR_RNR", "TR_NoInfo", "RR_NoInfo")] = 1
  ar[clinical[, paste0("overall_response_", i)] %in% c("Neither")] = 0
  n = t(data.frame(as.numeric(table(ar, useNA="always"))))
  cn <- names(table(ar, useNA="always"))
  cn[is.na(cn)] <- "N/A"
  dimnames(n) = list(i, cn)
  any_resp = rbind(any_resp, n[, match(c("0", "1", "N/A"), colnames(n)) ])
  
  clinical[, paste0(i, "_toc_resp")]= tr 
  clinical[, paste0(i, "_rtx_resp")]= rr 
  clinical[, paste0(i, "_any_resp")]= ar 
}


# Organise overall response
rownames(endpoint_counts) = rownames(initial_counts) = rownames(overall_counts) = resp_metric
colnames(endpoint_counts) = colnames(initial_counts) = c("RR", "RNR", "TR", "TNR","NoInfo", "N/A")
colnames(overall_counts) = c("Both", "RR", "TNR_RR", "TR", "RNR_TR","TNR_NoInfo", "RNR_NoInfo", 
                             "RNR_Dropped out", "TNR_Dropped out", "Neither", "N/A")

colnames(rit_resp) = colnames(toc_resp) = colnames(any_resp) = c("0", "1", "N/A")
rownames(rit_resp) = rownames(toc_resp) = rownames(any_resp) = resp_metric


metadata$Patient_ID <- gsub("\\.", "", metadata$Patient_ID)

response_df <- merge(metadata, clinical, by.x="Patient_ID", by.y="Patient.I.D.", all=TRUE)
dropped.out = as.character(clinical$Patient.I.D.[clinical$best.trt.any == "Early Drop-out"])

response_df = response_df[, grepl("Patient|seq|super|swtrt|Randomised|toc_resp|overall|rtx_resp|any_resp", 
                                  colnames(response_df))]
response_df = response_df[! duplicated(response_df[, colnames(response_df) != "super_ID"]), ]
response_df = response_df[! is.na(response_df$Randomised.medication), ]

response_df$Initial.medication = response_df$Randomised.medication
response_df$Second.medication = response_df$swtrt
response_df$Second.medication[response_df$Second.medication %in% c("RITtoTOC", "TOC")] = "Tocilizumab"
response_df$Second.medication[response_df$Second.medication %in% c("TOCtoRIT", "RIT")] = "Rituximab"
response_df$Second.medication[response_df$Patient_ID %in% dropped.out] = "Dropped out"


med = data.frame(table(response_df$Initial.medication, response_df$Second.medication, useNA="ifany"))
colnames(med)[1:2] = c("Initial Medication", "Second Medication")


cat("Patients who dropped out after visit 2:", 
    paste0(response_df$Patient_ID[response_df$Second.medication == "Dropped out"], collapse=", "))


# Baseline response
baseline_df <- merge(metadata, clinical, by.x="Patient_ID", by.y="Patient.I.D.", all=TRUE)
dropped.out = as.character(clinical$Patient.I.D.[clinical$best.trt.any == "Early Drop-out"])

baseline_df = baseline_df[, grepl("Patient|seq|super|swtrt|Randomised|Before", 
                                  colnames(baseline_df))]
baseline_df = baseline_df[! duplicated(baseline_df[, colnames(baseline_df) != "super_ID"]), ]
baseline_df = baseline_df[! is.na(baseline_df$Randomised.medication), ]
baseline_df[baseline_df == "Responder"] = 1
baseline_df[baseline_df == "Non.Responder"] = 0

toc_resp = baseline_df[baseline_df$Randomised.medication == "Tocilizumab", ]
rtx_resp = baseline_df[baseline_df$Randomised.medication == "Rituximab", ]

for(i in colnames(baseline_df)[grepl("Before", colnames(baseline_df))]){
  baseline_df[, gsub("\\.Before", "_toc_resp", i)] = toc_resp[match(baseline_df$Patient_ID, toc_resp$Patient_ID), i]
  baseline_df[, gsub("\\.Before", "_rtx_resp", i)] = rtx_resp[match(baseline_df$Patient_ID, rtx_resp$Patient_ID), i]
  baseline_df[, gsub("\\.Before", "_any_resp", i)] =  baseline_df[, gsub("\\.Before", "_toc_resp", i)]
  baseline_df[is.na(baseline_df[, gsub("\\.Before", "_any_resp", i)]), gsub("\\.Before", "_any_resp", i)] = 
    baseline_df[is.na(baseline_df[, gsub("\\.Before", "_any_resp", i)]), gsub("\\.Before", "_rtx_resp", i)]
}

baseline_df = baseline_df[, !grepl("Before", colnames(baseline_df))]



# Feature filtering by variance
###########################
# Chris' function
# vstx <- featurefilter(vst, method='A', percentile=15) # v1 was A and percentile=20
# exp = vstx$filtered_data

# my method
vars = apply(vst, 1, var)
vars = sort(vars, decreasing = T)
keep = vars[1:(ceiling(length(vars)*0.10))] # keep the top 25 for ~2400 features or 10% for ~1400 features)
exp = vst[rownames(vst) %in% names(keep), ]


exp_sd_uncor = t(exp[, response_df$seqID])

exp_sd = cor(exp_sd_uncor)
exp_sd[is.na(exp_sd)] = 0
hc = findCorrelation(exp_sd, cutoff=0.9) # remove genes where r > 0.95 (was 0.9 for ~ 1400 features)
hc = sort(hc)


output_clustering = FALSE

if(output_clustering){
  # hc are those to be removed
  library(ComplexHeatmap)
  library(circlize)
  library(dendextend)
  
  # hc are the genes to be removed
  hm_df = exp_sd[hc, -c(hc)] # rows are those to be removed, columns are those to be kept
  
  keep = apply(hm_df, 2, function(x) any(abs(x) > 0.95) ) # only look at the genes which actually cluster with those to be removed
  table(keep)
  hm_df = hm_df[, keep]
  
  pa = cluster::pam(hm_df, k = 5)
  pa2 = cluster::pam(t(hm_df), k = 5)
  
  spl = list("A"=c(1, 5), "B"=c(2, 4), "C"=c(3, 3), "D"=c(4, 2), "E"=c(5, 1))
  
  ht_master = Heatmap(hm_df, name = "correlation", 
                      show_row_names = F, show_column_names = F,
                      layer_fun = function(j, i, x, y, width, height, fill, slice_r, slice_c) {
                        v = pindex(hm_df, i, j)
                        if(slice_c == spl[[slice_r]][2]) {
                          grid.rect(gp = gpar(lwd = 2, fill = "transparent"))
                        }
                      },
                      row_split = factor(pa$clustering),
                      column_split = factor(pa2$clustering)
  )
  draw(ht_master)
  
  spls = list(c(1, 2), c(2, 1), c(5, 3), c(3, 5), c(4, 4))
  col_fun = colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
  
  # zoom in on each cluster
  pdf("/media/gcpeac/R4RA_machine_learning/Outputs/Plots/corr_clusters.pdf", height=18, width=30)
  pushViewport(viewport(layout = grid.layout(nrow=3, ncol=2, 
                                             widths=unit(c(0.9,0.1), "null"), 
                                             heights=unit(c(0.1,0.9, 0.1), "null"))))
  pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
  grid.text("Correlation matrix for feature selection", gp=gpar(fontsize=20), check=TRUE)
  popViewport()
  pushViewport(viewport(layout.pos.col=1, layout.pos.row=2))
  grid.draw(grid.grabExpr(draw(ht_master, heatmap_legend_side = "left"))) 
  popViewport()
  pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 2))
  grid.text("Dropped Features", rot=270, gp=gpar(fontsize=20), check=TRUE)
  popViewport()
  pushViewport(viewport(layout.pos.row = 3, layout.pos.col = 1))
  grid.text("Remaining Features which Correlate", gp=gpar(fontsize=20), check=TRUE)
  grid.newpage()
  for(i in spls){
    pushViewport(viewport(layout = grid.layout(nrow=2, ncol=2, 
                                               height=unit(c(length(which(pa2$clustering==i[2]))+1000, 480-length(which(pa2$clustering==i[2]))), "null"), 
                                               widths=unit(c(length(which(pa$clustering==i[1]))+4800, 102-length(which(pa$clustering==i[1]))), "null"))))
    clust_df = hm_df[pa$clustering==i[1], pa2$clustering ==i[2]]
    ht = Heatmap(t(clust_df), name = "correlation", col=col_fun, row_names_gp = gpar(fontsize = 6),
                 column_names_gp = gpar(fontsize = 6),
                 cell_fun = function(j, i, x, y, width, height, fill) {
                   if(abs(clust_df[j, i]) > 0.95) {
                     grid.rect(x = x, y = y, width = width, height = height, 
                               gp = gpar(col = "black", fill = "red"))
                     
                   }},
                 column_title = paste0("Cluster", names(which(unlist(lapply(spl, "[[", 1)) == i[1])), ": row ", i[1], ", column ", i[2]),
                 row_title = "KEPT Features", row_title_side ="right",
                 show_row_names = T, show_column_names = T)
    
    pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
    grid.draw(grid.grabExpr(draw(ht))) # or draw(ht_list, newpage = FALSE)
    grid.newpage()
  }
  dev.off()
  
  
  # dropped genes
  paste(rownames(hm_df)[! grepl("LINC|AC1|AL1|AP0|AL3|AL2|AL4|AC0|AL0|AL2|AL3|AL5|AL6|AL8|AP4|AC2", rownames(hm_df))], collapse=", ")
  length(rownames(hm_df)[grepl("LINC|AC1|AL1|AP0|AL3|AL2|AL4|AC0|AL0|AL2|AL3|AL5|AL6|AL8|AP4|AC2", rownames(hm_df))])
}


########################


exp_sd = exp_sd_uncor[, -c(hc)]
cat('going from', ncol(exp_sd_uncor), 'genes to', ncol(exp_sd))

loop_resp = paste0(resp_metric[resp_metric != "summed"],
                   rep(c("_toc_resp", "_rtx_resp", "_any_resp"),
                       each=length(resp_metric[resp_metric != "summed"])))


# focus on CDAI50
loop_resp = loop_resp[grepl("CDAI50", loop_resp)]

clin_exp = data.frame(cbind("ID"=rownames(exp_sd), 
                            "Initial"=as.character(response_df$Initial.medication), 
                            "Second"=as.character(response_df$Second.medication)), 
                      stringsAsFactors = F)


# Load in more metadata
m= readRDS("/media/gcpeac/Katriona/R4RA/full_meta_data.rds")

mbl = m[m$Visit == 3 & m$seqID %in% clin_exp$ID & m$seqID %in% rownames(exp_sd), ]
mbl = mbl[, c("Patient_ID" , 'seqID', 'Visit', 'TJC', 'SJC', 'CDAI', 'DAS28.CRP', 'DAS28.ESR',  
              'Age',  'Gender', 'ESR', 'CRP')]
m2 =  m[m$Visit == 2 & m$seqID %in% clin_exp$ID & m$seqID %in% rownames(exp_sd), ]
m2 = m2[, c("Patient_ID" ,  "CD20", "CD138" ,  "CD68L",   "CD68SL", "CD3")]
mbl$Gender = as.numeric(mbl$Gender)

clin_exp = cbind(mbl, m2[match(mbl$Patient_ID, m2$Patient_ID), ])
clin_exp= clin_exp[match(rownames(exp_sd), clin_exp$seqID), ]
clin_exp = clin_exp[, ! colnames(clin_exp) %in% c("Patient_ID" , 'seqID', 'Visit')]
clin_vars = colnames(clin_exp)


save(loop_resp, exp_sd, exp, exp_sd_uncor, clinical, vst, metadata, 
     baseline_df, response_df, resp_metric, clin_exp, clin_vars, 
     file="/media/gcpeac/R4RA_machine_learning/Data/curated_input.rdata")
