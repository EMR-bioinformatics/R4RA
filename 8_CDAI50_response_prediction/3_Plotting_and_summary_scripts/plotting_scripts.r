# Plots used by KG for Machine Learning Results: 
#      ROCs, variable importance, calibration curves, coefficient plots

drug_shorthands = setNames(c("Tocilizumab", "Rituximab", "Any"), c("toc", "rtx", "any"))
resp_names = setNames(c("CDAI50", "CDAI MTR", "EULAR ESR (bin1)", "EULAR ESR (bin2)",
                        "DAS28 ESR", "EULAR CRP (bin1)", "EULAR CRP (bin2)", "DAS28 CRP"), resp_metric)

# test and left-out ROCs on one plot
test_and_lo = function(modelOutput, color, i){
  
  rfFit = modelOutput$fit
  
  accuracy_lo = rfFit$results
  lo_df = rfFit$pred
  for(colu in colnames(rfFit$bestTune)){
    accuracy_lo = accuracy_lo[accuracy_lo[, colu] == rfFit$bestTune[, colu], ]
    lo_df = lo_df[lo_df[, colu] == rfFit$bestTune[, colu], ]
  }
  
  D.ex = lo_df$obs
  M1 = lo_df$R
  
  lo <- data.frame(D = D.ex, M1 = M1, set="left-out", stringsAsFactors = FALSE)
  
  
  test_set <- predict(rfFit, newdata=modelOutput$probs, type="prob")
  test_set$predicted = predict(rfFit, newdata=modelOutput$probs)
  test_set$observed = modelOutput$probs[, "resp"]
  
  
  
  combined = rbind(lo, data.frame("D"=test_set$observed,  "M1"=test_set$R, "set"="test", stringsAsFactors = FALSE))
  
  roc_test <-  ggplot(combined[combined$set == "test", ], aes(d = D, m = M1)) +
    geom_roc(color=color) +
    geom_abline(slope=1, intercept=0, color="grey60", linetype="dashed")
  roc_lo <-  ggplot(combined[combined$set == "left-out", ], aes(d = D, m = M1)) +
    geom_roc(color=color) +
    geom_abline(slope=1, intercept=0, color="grey60", linetype="dashed")
  
  roc <- ggplot(combined, aes(d = D, m = M1, group=set, linetype=set, color=set)) +
    geom_roc(labels=FALSE) +
    scale_color_manual(values=c("black", color)) +
    scale_linetype_manual(values=c("dashed", "solid")) +
    geom_abline(slope=1, intercept=0, color="grey60", linetype="dashed")
  
  title = unlist(strsplit(i, "_"))
  title = paste(resp_names[title[1]], drug_shorthands[title[2]], "Response")
  
  roc <- roc +
    annotate("text", x=.75, y=.1, color=color, vjust=1,size=4,
             label = paste("AUC =", round(calc_auc(roc_test)$AUC, 2))) +
    annotate("text", x=.75, y=.1,color="black", vjust=-1,size=4,
             label = paste("AUC =", round(calc_auc(roc_lo)$AUC, 2)))+
    labs(x="1 - Specificity", y="Sensitivity", subtitle=modelOutput$fit$method,
         title=title, color="Evaluation Set", linetype="Evaluation Set") + 
    theme_classic() + 
    theme(text = element_text(size=15))
  
  
  return(roc)
}

# ROC curves for left-out data
lo_rocs = function(modelList, metric, drug, ngenes=20457, subtitle=""){
  
  models = modelList[grepl(metric, names(modelList))]
  models = models[grepl(drug, names(models))]
  models = models[grepl(paste(ngenes, collapse="|"), names(models))]
  
  combined = data.frame()
  for(i in names(models)){
    
    rfFit = models[[i]]$fit
    
    accuracy_lo = rfFit$results
    lo_df = rfFit$pred
    for(colu in colnames(rfFit$bestTune)){
      accuracy_lo = accuracy_lo[accuracy_lo[, colu] == rfFit$bestTune[, colu], ]
      lo_df = lo_df[lo_df[, colu] == rfFit$bestTune[, colu], ]
    }
    
    D.ex = lo_df$obs
    M1 = lo_df$R
    
    s = paste0(models[[i]]$fit$method, " (AUC=", 
               format(as.numeric(as.character(models[[i]]$df$`AUC (LO)`)), digits=3), ")")
    
    lo <- data.frame(D = D.ex, M1 = M1, set=s, 
                     auc = rep(models[[i]]$df$`AUC (LO)`, length(D.ex)), 
                     stringsAsFactors = FALSE)
    
    combined = rbind(combined, lo)
  }
  combined = combined[order(combined$auc, decreasing = T), ]
  combined$set = factor(combined$set, levels=unique(combined$set))
  
  
  title = paste(metric, 
                setNames(c("Tocilizumab", "Rituximab"), c("toc", "rtx"))[drug], 
                "Response")
  
  roc <- ggplot(combined, aes(d = D, m = M1, group=set, color=set)) +
    geom_roc(labels=FALSE, pointsize=0, linemitre=1) +
    geom_abline(slope=1, intercept=0, color="grey60", linetype="dashed") + 
    labs(x="1 - Specificity", y="Sensitivity", 
         title=title, subtitle=subtitle, color="Model") + 
    theme_classic() + 
    theme(text = element_text(size=15))
  
  
  return(roc)
}

# ROC curves for test data
test_rocs = function(modelList, metric, drug, ngenes=20457, subtitle=""){
  
  models = modelList[grepl(metric, names(modelList))]
  models = models[grepl(drug, names(models))]
  if( ! is.null(ngenes)) 
    models = models[grepl(paste(ngenes, collapse="|"), names(models))]
  
  combined = data.frame()
  for(i in names(models)){
    
    rfFit = models[[i]]$fit
    
    test_set <- predict(rfFit, newdata=models[[i]]$probs, type="prob")
    test_set$predicted = predict(rfFit, newdata=models[[i]]$probs)
    test_set$observed = models[[i]]$probs[, "resp"]
    
    s = paste0(models[[i]]$fit$method, " (AUC=", format(models[[i]]$df$AUC, digits=3), ")")
    
    combined = rbind(combined, 
                     data.frame("D"=test_set$observed,  "M1"=test_set$R, 
                                "set"=rep(s, nrow(test_set)), "auc" = rep(models[[i]]$df$AUC, nrow(test_set)), 
                                stringsAsFactors = FALSE))
  }
  combined = combined[order(combined$auc, decreasing = T), ]
  combined$set = factor(combined$set, levels=unique(combined$set))
  
  
  title = paste(metric, 
                setNames(c("Tocilizumab", "Rituximab"), c("toc", "rtx"))[drug], 
                "Response")
  
  roc <- ggplot(combined, aes(d = D, m = M1, group=set, color=set)) +
    geom_roc(labels=FALSE, pointsize=0, linemitre=1) +
    geom_abline(slope=1, intercept=0, color="grey60", linetype="dashed") + 
    labs(x="1 - Specificity", y="Sensitivity", 
         title=title, subtitle=subtitle, color="Model") + 
    theme_classic() + 
    theme(text = element_text(size=15))
  
  
  return(roc)
}


final_rocs = function(modelList, metric, drug, subtitle=""){
  
  models = modelList
 
  combined = data.frame()
  all_resp = data.frame()
  for(i in 1:length(models$models)){
    m=models$models[[i]]
    rfFit = models$models[[i]]$fit
    
    test_set <- models$models[[i]]$prob
    all_resp = rbind(all_resp, test_set)
    
    s = paste0(models$models[[i]]$fit$method, " (AUC=", 
               format(models$models[[i]]$df$AUC, digits=3), ")")
    
    combined = rbind(combined, 
                     data.frame("D"=test_set$observed,  "M1"=test_set$R, 
                                "set"=rep(s, nrow(test_set)), "auc" = rep(models$models[[i]]$df$AUC, nrow(test_set)), 
                                stringsAsFactors = FALSE))
  }
  combined = combined[order(combined$auc, decreasing = T), ]
  combined$set = factor(combined$set, levels=unique(combined$set))
  
  auc_o = auc(all_resp$obs, all_resp$R, levels=c("NR", "R"), direction="<")
  
  title = paste(metric, 
                setNames(c("Tocilizumab", "Rituximab"), c("toc", "rtx"))[drug], 
                "Response")
  
  roc <- ggplot(combined, aes(d = D, m = M1)) +
    geom_roc(labels=F, pointsize=0, linemitre=1) +
    geom_abline(slope=1, intercept=0, color="grey60", linetype="dashed") + 
    labs(x="1 - Specificity", y="Sensitivity", 
         title=title, subtitle=paste("Outer (Test) AUC =", format(auc_o, digits=3)), color="Model") + 
    theme_classic() + 
    theme(text = element_text(size=15))
  
  
  return(roc)
}

# ROC curves for both test and left-out
all_rocs = function(modelList, metric, drug, ngenes=20457, subtitle=""){
  
  models = modelList[grepl(metric, names(modelList))]
  models = models[grepl(drug, names(models))]
  models = models[grepl(paste(ngenes, collapse="|"), names(models))]
  
  combined = data.frame()
  for(i in names(models)){
    
    rfFit = models[[i]]$fit
    training = models[[i]]$fit$trainingData
    colnames(training)[1] = "resp"
    alldat = rbind(training[, colnames(models[[i]]$probs)], models[[i]]$probs[, colnames(models[[i]]$probs)])
    
    test_set <- predict(rfFit, newdata=alldat, type="prob")
    test_set$predicted = predict(rfFit, newdata=alldat)
    test_set$observed = alldat[, "resp"]
    
    auc = auc(test_set$observed, test_set$NR, levels=c("NR", "R"), direction=">")
    
    s = paste0(models[[i]]$fit$method, " (AUC=", format(auc, digits=3), ")")
    
    combined = rbind(combined, 
                     data.frame("D"=test_set$observed,  "M1"=test_set$R, 
                                "set"=rep(s, nrow(test_set)), "auc" = rep(models[[i]]$df$AUC, nrow(test_set)), 
                                stringsAsFactors = FALSE))
  }
  combined = combined[order(combined$auc, decreasing = T), ]
  combined$set = factor(combined$set, levels=unique(combined$set))
  
  
  title = paste(metric, 
                setNames(c("Tocilizumab", "Rituximab"), c("toc", "rtx"))[drug], 
                "Response")
  
  roc <- ggplot(combined, aes(d = D, m = M1, group=set, color=set)) +
    geom_roc(labels=FALSE, pointsize=0, linemitre=1) +
    geom_abline(slope=1, intercept=0, color="grey60", linetype="dashed") + 
    labs(x="1 - Specificity", y="Sensitivity", 
         title=title, subtitle=subtitle, color="Model") + 
    theme_classic() + 
    theme(text = element_text(size=15))
  
  
  return(roc)
}

# ROC curves for test and left-out side by side for all models
side_by_side_rocs = function(modelList, metric, drug, ngenes=20457, fs=10, linesize=1, dig=2,
                             cs = setNames(c("#999999", "#E69F00", "#56B4E9", "#009E73",
                                             "#F0E442", "#0072B2", "#D55E00", "#CC79A7"),  
                                           c("rf", "glmnet", "svmRadial", "nb", "pda", "gbm", "svmPoly", "mda"))){
  
  
  models = modelList[grepl(metric, names(modelList))]
  models = models[grepl(drug, names(models))]
  models = models[unlist(lapply(names(models), function(x) strsplit(x, " ")[[1]][4])) %in% ngenes]
  
  combined = data.frame()
  for(i in names(models)){
    
    rfFit = models[[i]]$fit
    
    test_set <- predict(rfFit, newdata=models[[i]]$probs, type="prob")
    test_set$predicted = predict(rfFit, newdata=models[[i]]$probs)
    test_set$observed = models[[i]]$probs[, "resp"]
    
    # auc(test_set$observed, test_set$R, levels=c("NR", "R"), direction=">")
    
    s = paste0(models[[i]]$fit$method, " (AUC=", format(models[[i]]$df$AUC, digits=dig), ")")
    test = data.frame("D"=test_set$observed,  "M1"=test_set$R, 
                      "set"=rep(s, nrow(test_set)), "auc" = rep(models[[i]]$df$AUC, nrow(test_set)), 
                      eval="test",
                      stringsAsFactors = FALSE)
    
    accuracy_lo = rfFit$results
    lo_df = rfFit$pred
    for(colu in colnames(rfFit$bestTune)){
      accuracy_lo = accuracy_lo[accuracy_lo[, colu] == rfFit$bestTune[, colu], ]
      lo_df = lo_df[lo_df[, colu] == rfFit$bestTune[, colu], ]
    }
    
    
    D.ex = lo_df$obs
    M1 = lo_df$R
    
    s_lo = paste0(models[[i]]$fit$method, " (AUC=", 
                  format(models[[i]]$df$`AUC (LO)`, digits=2), ")")
    
    lo <- data.frame(D = D.ex, M1 = M1, set=s_lo, 
                     auc = rep(models[[i]]$df$`AUC (LO)`, length(D.ex)), 
                     eval="lo",
                     stringsAsFactors = FALSE)
    
    
    training = models[[i]]$fit$trainingData
    colnames(training)[1] = "resp"
    alldat = rbind(training[, colnames(models[[i]]$probs)], models[[i]]$probs[, colnames(models[[i]]$probs)])
    
    test_set <- predict(rfFit, newdata=alldat, type="prob")
    test_set$predicted = predict(rfFit, newdata=alldat)
    test_set$observed = alldat[, "resp"]
    
    auc = auc(test_set$observed, test_set$NR, levels=c("NR", "R"), direction=">")
    s = paste0(models[[i]]$fit$method, " (AUC=", format(auc, digits=dig), ")")
    all_test = data.frame("D"=test_set$observed,  "M1"=test_set$R, 
                          "set"=rep(s, nrow(test_set)), "auc" = rep(models[[i]]$df$AUC, nrow(test_set)), 
                          eval="all",
                          stringsAsFactors = FALSE)
    
    combined = rbind(combined, rbind(test, lo, all_test))
  }
  combined = combined[order(combined$auc, decreasing = T), ]
  combined$set = factor(combined$set, levels=unique(combined$set))
  combined$model = unlist(lapply(strsplit(as.character(combined$set), split=" "), "[[", 1))
  cs_test = as.character(cs[match(unique(combined[combined$eval == "test", "model"]), names(cs))])
  cs_lo = as.character(cs[match(unique(combined[combined$eval == "lo", "model"]), names(cs))])
  cs_all = as.character(cs[match(unique(combined[combined$eval == "all", "model"]), names(cs))])
  
  title = paste(metric, 
                setNames(c("Tocilizumab", "Rituximab"), c("toc", "rtx"))[drug], 
                "Response")
  
  roc_test <- ggplot(combined[combined$eval == "test", ], aes(d = D, m = M1, group=set, color=set)) +
    geom_roc(labels=F, pointsize=0, linemitre=1, size=linesize) +
    geom_abline(slope=1, intercept=0, color="grey60", linetype="dashed") + 
    labs(x="1 - Specificity", y="Sensitivity", 
         title=title, subtitle="Test set", color="Model") + 
    theme_classic() + 
    theme(text = element_text(size=fs)) + 
    scale_color_manual(values = cs_test) 
  
  roc_lo <- ggplot(combined[combined$eval == "lo", ], aes(d = D, m = M1, group=set, color=set)) +
    geom_roc(labels=FALSE, pointsize=0, linemitre=1, size=linesize) +
    geom_abline(slope=1, intercept=0, color="grey60", linetype="dashed") + 
    labs(x="1 - Specificity", y="Sensitivity", title="",
         subtitle="Left out set", color="Model") + 
    theme_classic() + 
    theme(text = element_text(size=fs)) + 
    scale_color_manual(values = cs_lo) 
  
  roc_all <- ggplot(combined[combined$eval == "all", ], aes(d = D, m = M1, group=set, color=set)) +
    geom_roc(labels=FALSE, pointsize=0, linemitre=1, size=linesize) +
    geom_abline(slope=1, intercept=0, color="grey60", linetype="dashed") + 
    labs(x="1 - Specificity", y="Sensitivity", title="",
         subtitle="All samples", color="Model") + 
    theme_classic() + 
    theme(text = element_text(size=fs)) + 
    scale_color_manual(values = cs_all) 
  
  return(list("test"=roc_test, "lo"=roc_lo, "all"=roc_all, "combined"=ggarrange(roc_test, roc_lo, ncol=2)))
}

# Calibration curves for all models
calibration_curve = function(modelList, metric, drug, ngenes=20457, fs=10, cuts=7,linesize=1, n_sig=3, 
                             cs = setNames(c("#999999", "#E69F00", "#56B4E9", "#009E73",
                                             "#F0E442", "#0072B2", "#D55E00", "#CC79A7"),  
                                           c("rf", "glmnet", "svmRadial", "nb", "pda", 
                                             "gbm", "svmPoly", "mda"))){
  
  
  models = modelList[grepl(metric, names(modelList))]
  models = models[grepl(drug, names(models))]
  
  
  models = models[unlist(lapply(strsplit(names(models), " "), "[[", 4)) %in% ngenes]
  
  
  aucs = setNames(as.numeric(lapply(models, function(x) x$df$AUC)), names(models))
  names(aucs[order(aucs, decreasing=T)])
  
  models = models[unique(c(paste(drug, metric, "glmnet", ngenes[1]), 
                           names(aucs[order(aucs, decreasing=T)])[! grepl("glmnet", names(aucs[order(aucs, decreasing=T)]))]  ))[1:n_sig]]
  
  
  combined = data.frame()
  for(i in names(models)){
    print(i)
    
    
    training = models[[i]]$fit$trainingData
    colnames(training)[1] = "resp"
    alldat = rbind(training[, colnames(models[[i]]$probs)], models[[i]]$probs[, colnames(models[[i]]$probs)])
    
    predy = cbind("Predicted"=predict(models[[i]]$fit, alldat[, colnames(alldat) != "resp"]),
                  "Observed"=alldat$resp, 
                  predict(models[[i]]$fit, alldat, type="prob"))
    
    
    #cal = calibration(obs ~ NR, data = models[[i]]$fit$pred, class=1)$data
    cal = calibration(Observed ~ R, data=predy, class="R", cuts=cuts)$data
    
    cal$model = models[[i]]$df$model
    
    
    combined = rbind(combined, cal)
  }
  
  title = paste(metric, 
                setNames(c("Tocilizumab", "Rituximab"), c("toc", "rtx"))[drug], 
                "Response")
  
  ggplot(combined, aes(x=midpoint, y=Percent, ymin=Lower, max=Upper, color=model)) +
    # geom_point() + 
    geom_line(size=linesize) +
    #geom_errorbar(width=1)+
    geom_abline(intercept=0, slope=1, colour="grey60", linetype="dashed") +
    lims(x=c(0, 100), y=c(0, 100)) +
    theme_classic() +
    theme(text = element_text(size=fs)) + 
    labs(x="Bin Midpoint", y="Observed Event Percentage", 
         title=title) +
    scale_color_manual(values = cs) 
}

# variable importance plots
varImp_plots = function(x, deg.list, clin_vars=colnames(clin_exp), titles, drop=T, int.list=NULL) {
  model.test = x$fit
  model = x$df$model
  drug = as.character(x$df$Drug)
  vi = varImp(model.test)
  
  coefs =  data.frame(as.matrix(vi$importance))
  colnames(coefs)[1] = "Overall"
  coefs$gene = rownames(coefs)
  coefs = coefs[order(abs(coefs$Overall), decreasing = T), ]
  if(nrow(coefs) > 20) coefs= coefs[1:20, ]
  coefs = coefs[order(abs(coefs$Overall)), ]
  
  coefs = coefs[coefs$Overall != 0, ]
  
  coefs$source = ifelse(coefs$gene %in% clin_vars, "Clinical", "EEG")
  coefs$source[coefs$gene %in% c('CD20', 'CD138', 'CD68L', 'CD68SL', 'CD3')] = "Histological"
  
  
  deg = deg.list[[drug]]
  coefs$source[coefs$gene %in% deg] = "DEG"
  
  interaction = rownames(int.list[[drug]])
  coefs$source[coefs$gene %in% interaction] = "Gene interaction"
  
  coefs$source = factor(coefs$source, levels= c("Histological", "Clinical", "DEG", "EEG", "Gene interaction"))
  
  
  coefs$gene[coefs$gene %in% interaction] = gsub("and", "/", coefs$gene[coefs$gene %in% interaction])
  coefs$gene = factor(coefs$gene, levels=unique(coefs$gene))
  
  st = titles[names(titles) == as.character(x$df$n_input_genes)]
  
  ggplot(coefs, aes(y=Overall, x=gene, fill=source)) +
    geom_bar(stat="identity") +
    coord_flip() +
    labs(y="Importance", x="") +
    theme_classic() + 
    theme(text=element_text(size=10), panel.grid = element_line(colour = "grey92"),
          panel.grid.minor = element_line(size = rel(0.5)), 
          panel.grid.major = element_line(size = rel(0.7))) +
    
    scale_fill_manual(values=c("goldenrod1", "#5858da", "midnightblue", "mediumseagreen", "coral"),
                      breaks=c("Clinical", "EEG", "DEG", "Histological", "Gene interaction"), drop=drop, name="Feature Type") +
    labs(title=paste(drug_shorthands[as.character(x$df$Drug)],  model),
         subtitle=paste0(st, " (n=", x$df$n_input_genes, ")"))
}

# Coefficient plots
coef_plots = function(x, deg.list, clin_vars = colnames(clin_exp), titles, drop=T){
  model.test = x$fit
  model = x$df$model
  drug = as.character(x$df$Drug)
  if(model =="glmnet"){
    tmp_coeffs <- coef(model.test$finalModel, s = model.test$bestTune$lambda)
  }else(tmp_coeffs <- coef(model.test$finalModel))
  coefs = data.frame(as.matrix(tmp_coeffs))
  colnames(coefs)="coef"
  coefs$gene = rownames(coefs)
  coefs = coefs[order(abs(coefs$coef)), ]
  coefs$gene = factor(coefs$gene, levels=unique(coefs$gene))
  coefs$sign = factor(sign(coefs$coef))
  coefs = coefs[coefs$sign != 0, ]
  coefs = coefs[! coefs$gene %in% c("(Intercept)", "Intercept"), ]
  coefs$source = ifelse(coefs$gene %in% clin_vars, "Clinical", "EEG")
  coefs$source[coefs$gene %in% c('CD20', 'CD138', 'CD68L', 'CD68SL', 'CD3')] = "Histological"
  
  deg = deg.list[[drug]]
  coefs$source[coefs$gene %in% deg] = "DEG"
  
  coefs$source = factor(coefs$source, levels= c("Histological", "Clinical", "DEG", "EEG"))
  
  
  
  if(nrow(coefs) > 20) coefs = coefs[1:20, ]
  
  st = titles[names(titles) == as.character(x$df$n_input_genes)]
  
  
  ggplot(coefs, aes(y=coef, x=gene, fill=source)) +
    geom_bar(stat="identity") +
    coord_flip() +
    labs(y="Coefficient", x="") +
    scale_fill_manual(values=c("goldenrod1", "mediumseagreen", "midnightblue", "#5858da"), 
                      breaks=c("Clinical", "Histological", "DEG", "EEG"), drop=drop, name="Feature Type") +
    theme_classic() + theme(text=element_text(size=10), panel.grid = element_line(colour = "grey92"),
                            panel.grid.minor = element_line(size = rel(0.5)), 
                            panel.grid.major = element_line(size = rel(0.7))) +
    labs(title=paste0(drug_shorthands[as.character(x$df$Drug)], " ", model),
         subtitle=paste0(st, " (n=", x$df$n_input_genes, ")"))
}

# Confusion matrix
confusion_matrix <- function(model){
  cm = data.frame(model$cm$table)
  ggplot(cm, aes(fill=Prediction, y=Freq, x=Reference)) +
    geom_bar(position="fill", stat="identity") +
    scale_fill_manual(values=c("goldenrod1", "midnightblue")) +
    theme_classic() +
    labs(title = paste(model$df$Response_metric,
                       drug_shorthands[as.character(model$df$Drug)], 
                       model$df$model),
         subtitle=paste0("n=", model$df$n_input_genes)) +
    theme(text = element_text(size=10))
}

# ROC curves for nested cross-validation
nested_rocs = function(modelList, metric, drug, filtering="deg", fs=10, linesize=1, 
                       legpos=c(0.7, 0.2), k = FALSE, 
                       cs = setNames(c("#999999", "#E69F00", "#56B4E9", "#009E73",
                                       "#F0E442", "#0072B2", "#D55E00", "#CC79A7"),  
                                     c("rf", "glmnet", "svmRadial", "nb", "pda", 
                                       "gbm", "svmPoly", "mda"))){
  
  models = modelList[grepl(metric, names(modelList))]
  models = models[grepl(drug, names(models))]
  models = models[grepl(filtering, names(models))]
  if(class(k) == "numeric") models = models[grepl(paste0("\\_", k), names(models))]
  
  combined = data.frame()
  stats_df = data.frame()
  for(j in names(models)){
    temp = models[[j]]
    mod.combined = data.frame() # combine the info for each fold
    
    for(f in 1:length(temp$models)){
      fold = temp$models[[f]]
      
      rfFit = fold$fit
      
      test_set <- fold$prob
      
   
      
      test = data.frame("D"=test_set$observed,  "M1"=test_set$R, 
                        "Predicted"=test_set$predicted,
                        eval="test", model=fold$df$model, fold=f, 
                        final_k = length(fold$selected_features),
                        stringsAsFactors = FALSE)
      
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

    t = summary_stats(df = mod.combined[mod.combined$eval == "test", ], coln = "Test Stats")
    lo = summary_stats(df = mod.combined[mod.combined$eval == "lo", ], coln = "Left-out Stats")
    all = summary_stats(df = mod.combined[mod.combined$eval == "all", ], coln = "Overall Stats")
    
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
                              format(as.numeric(as.character(mod.combined$auc)), digits=2), 
                              ", BA=", format(mod.combined$balanced_accuracy, digits=2), ")")

    combined = rbind(combined, mod.combined)
  }
  
  combined = combined[order(combined$auc, decreasing = T), ]
  combined$set = factor(combined$set, levels=unique(combined$set))
  combined$model = unlist(lapply(strsplit(as.character(combined$set), split=" "), "[[", 1))
  
  cs_test = setNames(as.character(cs[match(gsub(" .*", "", levels(droplevels(combined$set[combined$eval=="test"]))), names(cs))]), 
                     levels(droplevels(combined$set[combined$eval=="test"])))
  cs_lo = setNames(as.character(cs[match(gsub(" .*", "", levels(droplevels(combined$set[combined$eval=="lo"]))), names(cs))]), 
                   levels(droplevels(combined$set[combined$eval=="lo"])))
  cs_all = setNames(as.character(cs[match(gsub(" .*", "", levels(droplevels(combined$set[combined$eval=="all"]))), names(cs))]), 
                    levels(droplevels(combined$set[combined$eval=="all"])))
  
  title = paste(metric, 
                setNames(c("Tocilizumab", "Rituximab", "Refractory"), 
                         c("toc", "rtx", "any"))[drug], 
                "Response")
  
  roc_test <- ggplot(combined[combined$eval == "test", ], 
                     aes(d = D, m = M1, group=set, color=set)) +
    geom_roc(labels=F, pointsize=0, linemitre=1, size=linesize) +
    geom_abline(slope=1, intercept=0, color="grey60", linetype="dashed") + 
    labs(x="1 - Specificity", y="Sensitivity", 
         title=title, subtitle="Test set", color="Model") + 
    theme_classic() + 
    theme(text = element_text(size=fs), legend.position=legpos) + 
    scale_color_manual(values = cs_test) 
  
  roc_lo <- ggplot(combined[combined$eval == "lo", ], 
                   aes(d = D, m = M1, group=set, color=set)) +
    geom_roc(labels=FALSE, pointsize=0, linemitre=1, size=linesize) +
    geom_abline(slope=1, intercept=0, color="grey60", linetype="dashed") + 
    labs(x="1 - Specificity", y="Sensitivity", title="",
         subtitle="Left out set", color="Model") + 
    theme_classic() + 
    theme(text = element_text(size=fs), legend.position=legpos) + 
    scale_color_manual(values = cs_lo) 
  
  roc_all <- ggplot(combined[combined$eval == "all", ], aes(d = D, m = M1, group=set, color=set)) +
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

summary_stats <- function(df, coln="Stats"){
  auc <-  auc(df$D, df$M1, levels=c("NR", "R"), direction="<")
 
  
  cm <- confusionMatrix(df$Predicted, df$D, positive="R")
  
  rc <- roc(df$D, df$M1, levels=c("NR", "R"), direction="<")
  vals <- rbind("youden"=coords(rc, best.method="youden"))
  YI =  which(vals$specificity + vals$sensitivity - 1 == 
                max(vals$specificity + vals$sensitivity - 1))[1]
  
  # print(vals)
  
  
  output <- rbind(data.frame("Stat"=cm$byClass),
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

