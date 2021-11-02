# Functions to perform nested cross-validation with optional feature filtering

#' Function to perform recursive feature elimination
#' 
#' @param sdata Feature space data with response in column 'resp'
#' @param rfe_funcs Rfe function. See 
#' \code{\link[caret:rfeControl]{caret::rfeControl()}}
#' @param size The rfe sizes. See \code{\link[caret:rfe]{caret::rfe()}}
#' @param ml_method The machine learning method/model
#' @param eval_metric Variable to optimise. 
#' See \code{\link[caret:rfe]{caret::rfe()}}
#' @param train_control Caret trControl. 
#' See \code{\link[caret:rfe]{caret::rfe()}}
#' @return A subset feature space 
rfe_filtering = function(sdata, rfe_funcs, size, ml_method, 
                         eval_metric, train_control){
  rfe.ctrl = rfeControl(functions = rfe_funcs,
                        allowParallel =TRUE, 
                        verbose = FALSE)
  
  predictors = sdata[, colnames(sdata) != "resp"]
  
  rfe_out <- rfe(x = predictors, 
                 y = sdata$resp, 
                 sizes = size,
                 method=ml_method, 
                 metric = eval_metric, 
                 rfeControl = rfe.ctrl, 
                 trControl=train_control)
  
  sdata <- sdata[, c(predictors(rfe_out), "resp")]
  return(sdata)
}

#' Function to perform univariate filtering
#' 
#' @param sdata Feature space data with response in column 'resp'
#' @param size The target feature size. 
#' @param train_control Caret trControl. See 
#' \code{\link[caret:rfe]{caret::rfe()}}
#' @return A subset feature space 
univariate_filtering = function(sdata, train_control, size){
  filterCtrl <- sbfControl(functions = rfSBF)
  
  predictors = sdata[, colnames(sdata) != "resp"]
  
  rfe_out <- sbf(x = predictors, 
                 y = sdata$resp, 
                 sbfControl = filterCtrl, 
                 trControl=train_control)
  
  # subset by the number of times each feature is selected (min in size)
  nsize = min(size, predictors(rfe_out), na.rm=T)
  min_n = min(sort(table(unlist(rfe_out$variables)), decreasing = T)[1:nsize])
  features = names(table(unlist(rfe_out$variables)))[table(unlist(rfe_out$variables)) >= min_n]
  
  subdata <- sdata[, c(features, "resp")]
  return(subdata)
}

#' Function to predict response to treatment
#' 
#' @param variable The response variable column name in metadata
#' @param metadata Data frame containing ID and response columns 
#' (response column defined by variable)
#' @param features A data frame of the feature space
#' @param size The target size for feature filtering (default=c(20, 30, 50))
#' @param inner_folds The number of inner folds (default=10)
#' @param repeats The number of repeats (default=1)
#' @param eval_metric Metric to optimise for feature filtering if rfe 
#' (default="Accuracy"). See \code{\link[caret:rfe]{caret::rfe()}}
#' @param ml_method machine learning model/method (default="glmnet")
#' @param filtering_method Feature filtering method 
#' (options = c("univariate", "rfe", "none")
#' @param rfe_funcs Rfe function. See 
#' \code{\link[caret:rfeControl]{caret::rfeControl()}}
#' @param k Number of outer folds (default=5)
#' @param cores The number of cores (default=6)
#' @param remove_na Whether to remove samples with NA in feature space
#' @param verbose Logical for verbose
r4ra_nested_ml = function(variable, 
                          features,
                          size=c(20, 30, 50), 
                          inner_folds=10, 
                          repeats=1, 
                          eval_metric="Accuracy", 
                          ml_method="glmnet", 
                          filtering_method = "univariate",
                          rfe_funcs = caretFuncs, 
                          metadata = clinical,
                          k=5, 
                          cores=6, 
                          remove_na=F,
                          verbose=FALSE){
  if(verbose) print(variable)
  if(! filtering_method %in% c("univariate", "rfe", "none")) {
    stop('filtering_method must be in c("univariate", "rfe", "none")')
  }
  
  tvst <- data.frame(features)
  y = metadata[, variable]
  
  ## add response variable
  tvst$resp = factor(y, labels=c("NR", "R"))
  tvst = tvst[! is.na(tvst$resp), ]
  
  # Control the computational nuances of the train function
  train_control <- trainControl(method="repeatedcv", 
                                number=inner_folds,
                                repeats=repeats,
                                classProbs = T,
                                savePredictions = TRUE)
  
  # Create outer folds
  outer_folds = createFolds(tvst$resp, k)
  
  ptm <- proc.time()  # measure the time
  
  # Acquire a model for each outer fold (using multi-cores with a progress bar)
  models <- pbmclapply(1:k, function(k_fold) {
    
    fold = outer_folds[[k_fold]] 
    if(verbose) print(paste("fold", k_fold, ":", paste(fold, collapse=", ")))
    
    # Split into train and test set 
    testset <- tvst[fold, ]
    subdata <- tvst[-fold, ]
    
    # Create a tuning grid for glmnet
    tgrid=NULL
    lambda.seq <- exp(seq(log(1e-5), log(1e0), length.out = 20))
    alpha.seq <- 7:10/10
    if(ml_method == "glmnet"){
      tgrid = expand.grid(alpha = alpha.seq, lambda = lambda.seq)
    } 
    prep = NULL
    
    if(verbose) print("Performing filtering")
    if(filtering_method == "rfe") {
      subdata = rfe_filtering(subdata, rfe_funcs, size, ml_method, 
                              eval_metric, train_control)
    }
    if(filtering_method == "univariate") {
      subdata = univariate_filtering(subdata, train_control, size=size)
    }
    if(verbose) print("Filtered successfully")
    testset=testset[, colnames(subdata)]
    
    if(verbose) print(paste("Feature dim: ", dim(subdata)))
    
    # Remove samples containing NA in features if required
    if(remove_na){
      subdata = subdata[complete.cases(subdata), ]
      testset = testset[complete.cases(testset), ]
    }
    
    # Train the model
    if(verbose) print("Training the model")
    fit <- caret::train(resp ~., 
                        data=subdata,
                        method = ml_method,
                        metric = eval_metric,
                        trControl = train_control,
                        tuneGrid = tgrid, 
                        preProcess=prep)
    
    if(verbose) print("Fit successfully")
    
    # Calculate evalutation metrics
    if(verbose) print("Performing evalutation calculations")
    accuracy_lo = fit$results
    lo_df = fit$pred
    for(colu in colnames(fit$bestTune)){
      accuracy_lo = accuracy_lo[accuracy_lo[, colu] == fit$bestTune[, colu], ]
      lo_df = lo_df[lo_df[, colu] == fit$bestTune[, colu], ]
    }
    
    accuracy_lo = accuracy_lo$Accuracy
    auc_lo = auc(lo_df$obs, lo_df$R, levels=c("NR", "R"), direction="<")
    
    test_fit <- predict(fit, newdata=testset, type="prob")
    test_fit$predicted = predict(fit, newdata=testset)
    test_fit$observed = testset$"resp"
    
    genes = as.character(fit$coefnames) 
    n_out = length(genes)
    genes = paste(genes, collapse="; ")
    ptest = predict(fit, newdata=testset)
    names(ptest) = rownames(testset)
    
    
    conf <- confusionMatrix(test_fit$predicted, test_fit$observed)
    if(all(is.na(test_fit$R)) | all(is.nan(test_fit$R))){
      auc = NA
    } else {
      auc <- auc(test_fit$observed, test_fit$R, levels=c("NR", "R"), 
                 direction="<")
    }
    
    cm = data.frame(conf$table)
    cm = data.frame("TP"=cm$Freq[cm$Prediction == "R" & cm$Reference == "R"], 
                    "FP"=cm$Freq[cm$Prediction == "R" & cm$Reference == "NR"], 
                    "TN"=cm$Freq[cm$Prediction == "NR" & cm$Reference == "NR"], 
                    "FN"=cm$Freq[cm$Prediction == "NR" & cm$Reference == "R"]
    )
    
    os = ncol(subdata) - 1
    data_output <- data.frame(gsub(".*_", "", gsub("_resp", "", variable)), 
                              gsub("_.*", "", variable), 
                              fit$method,  eval_metric, auc, auc_lo, 
                              accuracy_lo, os,
                              t(data.frame(conf$overall)),
                              t(data.frame(conf$byClass)), cm, 
                              genes,  n_out)
    
    colnames(data_output)[1:8] <- c("Drug", "Response_metric", "model", 
                                    "Optimising", "AUC",  "AUC (LO)", 
                                    "Accuracy (LO)", "n_genes_in")
    colnames(data_output)[(ncol(data_output)-1):ncol(data_output)] = 
      c("genes_output",  "n_output_genes")
    
    data_output = data_output[, c('Drug', 'Response_metric', 'model',  
                                  "Optimising", 'Balanced.Accuracy',
                                  'AUC', 'n_genes_in', 'Accuracy', 
                                  "AUC (LO)", "Accuracy (LO)",
                                  'Kappa', 'AccuracyLower', 'AccuracyUpper',
                                  'AccuracyNull', 'AccuracyPValue', 
                                  'McnemarPValue',
                                  'Sensitivity', 'Specificity', 'Pos.Pred.Value',
                                  'Neg.Pred.Value', 'Precision', 'Recall', 'F1',
                                  'Prevalence', 'Detection.Rate', 
                                  'Detection.Prevalence',
                                  'TP', 'FP', 'TN', 'FN',  'genes_output',
                                  'n_output_genes')]
    
    return(
      list("df"=data_output, 
           "fit"=fit, 
           "filtering"=filtering_method,
           "test"=testset,
           "prob"=test_fit,
           "tc"=train_control, 
           "predict"=as.numeric(ptest)-1,
           "obs"=as.numeric(factor(testset$resp))-1, 
           "cm"=conf, 
           "selected_features"=colnames(subdata)[colnames(subdata) !="resp"], 
           "metric"=variable))
    
  }, mc.cores = cores-1)
  proc.time() - ptm
  
  dfs = lapply(models, function(x) x$df)
  dfs = do.call(rbind, dfs)
  
  # Calculate the *mean* AUC
  auc = mean(dfs$AUC)
  
  return(list("models"=models, "df"=dfs, "auc"=auc, 
              "filtering"=filtering_method))
}
