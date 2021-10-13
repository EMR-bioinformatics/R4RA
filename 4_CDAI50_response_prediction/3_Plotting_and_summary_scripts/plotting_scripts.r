# Functions for machine learning visualisations and summaries

#' ROC curves for nested cross-validation
#' 
#' @param modelList A list of models
#' @param metric Response metric (in modelList names to filter by)
#' @param drug Drug (in modelList names to filter by)
#' @param filtering Filtering method (in modelList names to filter by)
#' @param fs Font size for ROC curves
#' @param linesize Line size for ROC curves
#' @param legpos Legend position
#' @param k If numeric (the target number of features) filters on names of 
#' modelList
#' @param cs Colour scheme for models (named vector)
nested_rocs = function(modelList, 
                       metric, 
                       drug, 
                       filtering="deg", 
                       fs=10, 
                       linesize=1, 
                       legpos=c(0.7, 0.2), 
                       k = FALSE, 
                       cs = setNames(c("#999999", "#E69F00", "#56B4E9", "#009E73",
                                       "#F0E442", "#0072B2", "#D55E00", "#CC79A7"),  
                                     c("rf", "glmnet", "svmRadial", "nb", "pda", 
                                       "gbm", "svmPoly", "mda"))){
  
  # filter the model list
  models = modelList[grepl(metric, names(modelList))]
  models = models[grepl(drug, names(models))]
  models = models[grepl(filtering, names(models))]
  if(class(k) == "numeric") 
    models = models[grepl(paste0("\\_", k, "\\_"), names(models))]
  
  
  # create data frame of stats and all test sets
  combined = data.frame()
  stats_df = data.frame()
  for(j in names(models)){
    temp = models[[j]]
    mod.combined = data.frame() # combine model info for each outer fold 
    
    # For each outer fold
    for(f in 1:length(temp$models)){
      fold = temp$models[[f]]
      rfFit = fold$fit
      test_set <- fold$prob
      
      # The predicted and observed results on the test set (outer folds)
      test = data.frame("D"=test_set$observed,  "M1"=test_set$R, 
                        "Predicted"=test_set$predicted,
                        eval="test", model=fold$df$model, fold=f, 
                        final_k = length(fold$selected_features),
                        stringsAsFactors = FALSE)
      
      
      # The predicted and observed results on the left-out set (inner folds)
      accuracy_lo = rfFit$results
      lo_df = rfFit$pred
      for(colu in colnames(rfFit$bestTune)){
        accuracy_lo = accuracy_lo[accuracy_lo[, colu] == rfFit$bestTune[, colu], ]
        lo_df = lo_df[lo_df[, colu] == rfFit$bestTune[, colu], ]
      }
      D.ex = lo_df$obs
      M1 = lo_df$R
      
      s_lo = paste0(fold$fit$method, " (AUC=", 
                    format(as.numeric(as.character(fold$df$`AUC (LO)`)), digits=3), ")")
      
      lo <- data.frame(D = D.ex, M1 = M1, 
                       "Predicted"=lo_df$pred,
                       eval="lo", model=fold$df$model,fold=f,
                       final_k = length(fold$selected_features),
                       stringsAsFactors = FALSE)
      
      # The predicted and observed results on entire set
      training = fold$fit$trainingData
      colnames(training)[1] = "resp"
      alldat = rbind(training[, colnames(fold$test)], 
                     fold$test[, colnames(fold$test)])
      
      test_set <- predict(rfFit, newdata=alldat, type="prob")
      test_set$predicted = predict(rfFit, newdata=alldat)
      test_set$observed = alldat[, "resp"]
      
      all_test = data.frame("D"=test_set$observed,  "M1"=test_set$R,  
                            "Predicted"=test_set$predicted,
                            eval="all",model=fold$df$model, fold=f, 
                            final_k = length(fold$selected_features),
                            stringsAsFactors = FALSE)
      
      mod.combined = rbind( mod.combined, rbind(test, lo, all_test))
    }
    
    # summarize results for each set
    t = summary_stats(df = mod.combined[mod.combined$eval == "test", ], 
                      coln = "Test Stats")
    lo = summary_stats(df = mod.combined[mod.combined$eval == "lo", ], 
                       coln = "Left-out Stats")
    all = summary_stats(df = mod.combined[mod.combined$eval == "all", ], 
                        coln = "Overall Stats")
    
    # combine into one data frame
    df = cbind(t, lo, all, "model"=j)
    df$file = temp$file
    stats_df = rbind(stats_df, df)
    
    mod.combined$auc = NA
    mod.combined$auc[mod.combined$eval == "test"] = df["AUC", "Test Stats"]
    mod.combined$auc[mod.combined$eval == "lo"] = df["AUC", "Left-out Stats"]
    mod.combined$auc[mod.combined$eval == "all"] = df["AUC", "Overall Stats"]
    
    mod.combined$balanced_accuracy = NA
    mod.combined$balanced_accuracy[mod.combined$eval == "test"] = 
      df["Balanced Accuracy", "Test Stats"]
    mod.combined$balanced_accuracy[mod.combined$eval == "lo"] = 
      df["Balanced Accuracy", "Left-out Stats"]
    mod.combined$balanced_accuracy[mod.combined$eval == "all"] = 
      df["Balanced Accuracy", "Overall Stats"]
    
    mod.combined$set = paste0(mod.combined$model, " (AUC=", 
                              format(as.numeric(as.character(mod.combined$auc)), 
                                     digits=2), 
                              ", BA=", format(mod.combined$balanced_accuracy, digits=2), ")")
    
    combined = rbind(combined, mod.combined)
  }
  
  combined = combined[order(combined$auc, decreasing = T), ]
  combined$set = factor(combined$set, levels=unique(combined$set))
  combined$model = unlist(lapply(strsplit(as.character(combined$set), split=" "), 
                                 "[[", 1))
  
  cs_test = as.character(cs[match(unique(combined[combined$eval == "test", "model"]), names(cs))])
  cs_lo = as.character(cs[match(unique(combined[combined$eval == "lo", "model"]), names(cs))])
  cs_all = as.character(cs[match(unique(combined[combined$eval == "all", "model"]), names(cs))])
  
  title = paste(metric, 
                setNames(c("Tocilizumab", "Rituximab", "Refractory"), 
                         c("toc", "rtx", "any"))[drug], 
                "Response")
  
  # Create a ROC curve for the test set (outer folds)
  roc_test <- ggplot(combined[combined$eval == "test", ], 
                     aes(d = D, m = M1, group=set, color=set)) +
    geom_roc(labels=F, pointsize=0, linemitre=1, size=linesize) +
    geom_abline(slope=1, intercept=0, color="grey60", linetype="dashed") + 
    labs(x="1 - Specificity", y="Sensitivity", 
         title=title, subtitle="Test set", color="Model") + 
    theme_classic() + 
    theme(text = element_text(size=fs), legend.position=legpos) + 
    scale_color_manual(values = cs_test) 
  
  # Create a ROC curve for the left out set (inner folds)
  roc_lo <- ggplot(combined[combined$eval == "lo", ], 
                   aes(d = D, m = M1, group=set, color=set)) +
    geom_roc(labels=FALSE, pointsize=0, linemitre=1, size=linesize) +
    geom_abline(slope=1, intercept=0, color="grey60", linetype="dashed") + 
    labs(x="1 - Specificity", y="Sensitivity", title="",
         subtitle="Left out set", color="Model") + 
    theme_classic() + 
    theme(text = element_text(size=fs), legend.position=legpos) + 
    scale_color_manual(values = cs_lo) 
  
  # Create a ROC curve for the entire set
  roc_all <- ggplot(combined[combined$eval == "all", ], 
                    aes(d = D, m = M1, group=set, color=set)) +
    geom_roc(labels=FALSE, pointsize=0, linemitre=1, size=linesize) +
    geom_abline(slope=1, intercept=0, color="grey60", linetype="dashed") + 
    labs(x="1 - Specificity", y="Sensitivity", title="",
         subtitle="All samples", color="Model") + 
    theme_classic() + 
    theme(text = element_text(size=fs), legend.position=legpos) + 
    scale_color_manual(values = cs_all) 
  
  return(list("test"=roc_test, "lo"=roc_lo, "all"=roc_all, 
              "combined"=ggarrange(roc_test, roc_lo, ncol=2), 
              "roc_df" = combined, "stats_df"=stats_df))
}




#' Function to create a table of summary statistics
#' 
#' @param df data frame of statistics
#' @param coln Column name to use in output
summary_stats <- function(df, coln="Stats"){
  # calculate AUC and confusion matrix from set
  auc <-  auc(df$D, df$M1, levels=c("NR", "R"), direction="<")
  cm <- confusionMatrix(df$Predicted, df$D, positive="R")
  
  # calculate Youden's sensitivity and specificity
  rc <- roc(df$D, df$M1, levels=c("NR", "R"), direction="<")
  vals <- rbind("youden"=coords(rc, best.method="youden"))
  YI =  which(vals$specificity + vals$sensitivity - 1 == 
                max(vals$specificity + vals$sensitivity - 1))[1]
  
  # generate an output data frame
  output <- rbind(
    data.frame("Stat"=cm$byClass),
    data.frame("Stat"=c(cm$overall, c("AUC"=auc))), 
    data.frame("Stat"=c("Youden_specificity"=vals[YI, "specificity"],
                        "Youden_sensitivity"=vals[YI, "sensitivity"])),
    data.frame("Stat"=c("final_k"=paste(paste(names(table(df$final_k)), 
                                              table(df$final_k), sep =": "), 
                                        collapse=", "))), 
    data.frame("Stat"=c("mean_inner_k"=mean(df$final_k))))
  
  
  colnames(output)=coln
  return(output)
  
}

