# This script produces modules to predict CDAI50 response in the R4RA cohort. 
#
# There are 2 possible response variables: 
#     1. predicting the final response (post treatment cross-over)
#     2. predicting the initial response (at the point of cross-over)
# There are also 2 possible feature sets: 
#     1. clinical and histology variables only
#     2. clinical and gene expression variables
# There are three treatment options: 
#     1. Rituximab (rtx)
#     2. Tocilizumab (toc)
#     3. Refractory (any), response to either


# ``````````````````````````````````````
# to run this script as a background process with nohup: 
# cd /media/gcpeac/R4RA_machine_learning/2_Machine_Learning_binary/nested_cv 
# nohup Rscript R4RA_CDAI50_prediction.R 
# ``````````````````````````````````````

library(caret)
library(ggplot2)
library(pROC)
library(plotROC)
library(dplyr)
library(data.table)
library(ggrepel)
library(magrittr)
library(ggpubr)
library(gbm)
library(parallel)
library(pbmcapply)
library(xgboost)


#############################
# Set up the working environment
#############################
set.seed(123)
setwd("/media/gcpeac/R4RA_machine_learning/2_Machine_Learning_binary/nested_cv/")
output_dir <- "/home/kgoldmann/Documents/R4RA_ML/Nested/"

dir.create(paste0(output_dir, "/Initial/"), showWarnings = FALSE)
dir.create(paste0(output_dir, "/Final/"), showWarnings = FALSE)
#############################

#############################
# Load in the data and source 
# additional scripts
#############################
source("./nested_machine_learning_functions.R")
load("/media/gcpeac/R4RA_machine_learning/Data/curated_input.rdata")

# Set up the response metrics 
resp_names <- setNames(c("CDAI50", "CDAI MTR", "EULAR ESR (bin1)", 
                         "EULAR ESR (bin2)", "DAS28 ESR", "EULAR CRP (bin1)", 
                         "EULAR CRP (bin2)", "DAS28 CRP"), resp_metric)
drug_shorthands <- setNames(c("Tocilizumab", "Rituximab", "Any"), 
                            c("toc", "rtx", "any"))

# Set up parameters for machine learning
models <- c("glmnet", "rf",  "svmPoly",  "svmRadial", "gbm",  "mda", "pda") 
clin_vars <- c('TJC', 'SJC', 'CDAI', 'DAS28.CRP', 'DAS28.ESR', 'Age', 
               'Gender', 'ESR', 'CRP') 
both_exp <- cbind(exp_sd, clin_exp[, clin_vars])
k <- 10

#############################

#############################
# Final CDAI50 Prediction using
# gene expression and clinical features
#############################

# Loop through different feature filtering sizes
for(size in c(25, 30, 50, 100)){
  print(paste("SIZE =", size))
  
  # Loop through different models
  for(m in models){
    if(m %in% c("mda", "pda")) feat <- exp_sd else feat <- both_exp
    print(paste( m, Sys.time(), "( model", which(models == m), "of", 
                 length(models), ")"))
    if(m == "glmnet") funcs <- lrFuncs else funcs <- caretFuncs
    
    # univariate filtering
    print(paste("---", m, "univariate", Sys.time()))
    o <- lapply(loop_resp, function(x) {
      print(paste("----", x, "univariate"))
      invisible(try(r4ra_nested_ml(x, size=size, rfe_funcs = funcs, ml_method= m,
                                   features = feat, k=k, inner_folds=10,
                                   metadata = clinical, 
                                   filtering_method = "univariate")))
    })
    assign(paste0("t_exp_and_clin_uni_", m), o)
    
    # filtering via recursive feature elimination
    print(paste("---", m, "rfe", Sys.time()))
    o <- lapply(loop_resp, function(x) {
      print(paste("--", x, "rfe"))
      invisible(try(r4ra_nested_ml(x, size=size, rfe_funcs = funcs, ml_method= m, 
                                   metadata = clinical, verbose=T,
                                   features = feat, k=k, inner_folds=10, 
                                   filtering_method = "rfe")))
    })
    assign(paste0("t_exp_and_clin_rfe_", m), o)
    
    # No feature filtering
    print(paste("---", m, "none", Sys.time()))
    o <- lapply(loop_resp, function(x) {
      print(paste("--", x, "no filtering"))
      invisible(r4ra_nested_ml(x, size=50, rfe_funcs = funcs, ml_method= m,
                               metadata = clinical, 
                               features = both_exp, k=k, inner_folds=10,
                               filtering_method = "none"))
    })
    assign(paste0("t_exp_and_clin_none_", m), o)
    
    dir.create(paste0(output_dir, "Final/input_", ncol(exp_sd),
                      "/"), showWarnings = FALSE)
    dir.create(paste0(output_dir, "Final/input_", ncol(exp_sd),
                      "/filter_size_", size), showWarnings = FALSE)
    save(list=ls()[grepl("t_exp_and", ls()) & grepl(m, ls())],
         file=paste0(output_dir, "Final/input_", ncol(exp_sd),
                     "/filter_size_", size, "/all_visit9_ml_nested_", 
                     m, "_", k, ".rdata"))
  }
}

#############################

#############################
# Final CDAI50 Prediction using
# clinincal & histology only
#############################
print("Clinical and Histological Models")
dir.create(paste0(output_dir, "Final/input_clinical"), showWarnings = FALSE)

for(m in models){
  if(m %in% c("mda", "pda")) {
    feat = clin_exp[complete.cases(clin_exp), ] 
  }else {
    feat =  clin_exp[, c('TJC', 'SJC', 'Age', 'Gender', 'ESR', 'CRP', 'CD20', 
                         'CD138', 'CD68L', 'CD68SL', 'CD3')]
  }
  print(paste( m, Sys.time(), "( model", which(models == m), "of", 
               length(models), ")"))
  if(m == "glmnet") funcs = lrFuncs else funcs = caretFuncs
  
  print(paste("---", m, "none", Sys.time()))
  o = lapply(loop_resp, function(x) {
    print(paste("--", x, "no filtering"))
    invisible(r4ra_nested_ml(x, size=50, rfe_funcs = funcs, ml_method= m,
                             metadata = clinical,
                             features = feat, k=k, inner_folds=10, remove_na=T,
                             filtering_method = "none"))
  })
  assign(paste0("t_clin_and_hist_", m), o)
  
  
  save(list=ls()[grepl("t_clin_and_hist", ls()) & grepl(m, ls())],
       file=paste0(output_dir, "/Final/input_clinical/all_visit9_clin_hist_ml_nested_", 
                   m, "_", k, ".rdata"))
}


#############################

#############################
# Initial CDAI50 Prediction using
# gene expression and clinical features
#############################
print("Initial Response Prediction")
models <- c("glmnet", "rf",  "svmPoly",  "svmRadial", "gbm",  "mda", "pda") 

baseline_df$Patient.I.D. <- baseline_df$Patient_ID


# Loop through different feature filtering sizes
for(size in c(30, 50)){  
  print(paste("SIZE =", size))
  
  # Loop through different models
  for(m in models){
    print(paste( m, Sys.time(), "( model", which(models == m), "of", 
                 length(models), ")"))
    if(m %in% c("mda", "pda")) feat = exp_sd else feat =  both_exp
    
    
    # univariate filtering
    print(paste("---", m, "univariate", Sys.time()))
    if(m == "glmnet") funcs = lrFuncs else funcs = caretFuncs
    o = lapply(loop_resp[1:2], function(x) {
      print(paste("----", x, "univariate"))
      invisible(try(r4ra_nested_ml(x, size=size, rfe_funcs = funcs, ml_method= m,
                                   metadata=baseline_df,
                                   features = feat, k=k, inner_folds=10,
                                   filtering_method = "univariate")))
    })
    assign(paste0("t_exp_and_clin_uni_", m), o)
    
    # filtering via recursive feature elimination
    print(paste("---", m, "rfe", Sys.time()))
    o = lapply(loop_resp[1:2], function(x) {
      print(paste("--", x, "rfe"))
      invisible(try(r4ra_nested_ml(x, size=50, rfe_funcs = funcs, ml_method= m,
                                   features = feat, k=k, inner_folds=10,
                                   filtering_method = "rfe")))
    })
    assign(paste0("t_exp_and_clin_rfe_", m), o)
    
    # no filtering
    print(paste("---", m, "none", Sys.time()))
    o = lapply(loop_resp[1:2], function(x) {
      print(paste("--", x, "no filtering"))
      invisible(try(r4ra_nested_ml(x, size=size, rfe_funcs = funcs, ml_method= m,
                                   metadata=baseline_df,
                                   features = feat, k=k, inner_folds=10,
                                   filtering_method = "none")))
    })
    assign(paste0("t_exp_and_clin_none2_", m), o)
    
    dir.create(paste0(output_dir, "/Initial/input_", ncol(exp_sd)),
               showWarnings = FALSE)
    dir.create(paste0(output_dir, "/Initial/input_", ncol(exp_sd),
                      "/filter_size_", size), showWarnings = FALSE)
    save(list=ls()[grepl("t_exp_and", ls()) & grepl(m, ls())],
         file=paste0(output_dir, "/Initial/input_", ncol(exp_sd),
                     "/filter_size_", size, "/all_visit7_ml_nested3_", m, 
                     "_", k, ".rdata"))
  }
}

#############################

#############################
# Initial CDAI50 prediction using
# clinincal & histology only
#############################
print("Initial clinical and Histological Models")
dir.create(paste0(output_dir, "Initial/input_clinical"), showWarnings = FALSE)

for(m in models){
  if(m %in% c("mda", "pda")) {
    feat = clin_exp[complete.cases(clin_exp), ] 
  }else {
    feat =  clin_exp[, c('TJC', 'SJC', 'Age', 'Gender', 'ESR', 'CRP', 'CD20', 
                         'CD138', 'CD68L', 'CD68SL', 'CD3')]
  }
  print(paste( m, Sys.time(), "( model", which(models == m), "of", 
               length(models), ")"))
  if(m == "glmnet") funcs = lrFuncs else funcs = caretFuncs
  
  print(paste("---", m, "none", Sys.time()))
  o = lapply(loop_resp[1:2], function(x) {
    print(paste("--", x, "no filtering"))
    invisible(r4ra_nested_ml(x, size=50, rfe_funcs = funcs, ml_method= m,
                             metadata=baseline_df,
                             features = feat, k=k, inner_folds=10, remove_na=T,
                             filtering_method = "none"))
  })
  assign(paste0("t_clin_and_hist_", m), o)
  
  
  save(list=ls()[grepl("t_clin_and_hist", ls()) & grepl(m, ls())],
       file=paste0(output_dir, "/Initial/input_clinical/all_visit7_clin_hist_ml_nested_", 
                   m, "_", k, ".rdata"))
}


#############################

