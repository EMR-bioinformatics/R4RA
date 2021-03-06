mergemodules <- function(data=NA,type,list=NA,method='zmean'){
  # for debugging
  #data <- vst
  #list <- dtdf
  #type <- 'singlecell'
  #method <- 'zmean' -> but script was set to mean -> I set it method='zmean'
  #
  if (type == 'wgcna'){
    df <- read.csv('../PEAC/8_modularframework/synovium_baseline_modules_annotatedv11.csv')
    df <- df[,-1]
    df <- df[,c('Symbol','Module')]
    df$Module <- as.character(df$Module)
    df$Symbol <- as.character(df$Symbol)
  }else if (type == 'KEGG'){
    df <- clusterProfiler::read.gmt('../Genesets/c2.cp.kegg.v6.1.symbols.gmt')
    colnames(df) <- c('Module','Symbol')
  }else if (type == 'REACTOME'){
    df <- clusterProfiler::read.gmt("../Genesets/c2.cp.reactome.v6.1.symbols.gmt")
    colnames(df) <- c('Module','Symbol')
  }else if (type == 'blood'){
    # modules (chaussabel)
    df <- read.csv('../PEAC/1_BloodDiseaseSubtypes/chassmodulesV3forPEAC.csv')
    df <- df[,-1]
    df <- df[,-1]
  }else if (type == 'cos_beloved_modules'){
    cbm <- readRDS("./R4RA/5_modules/Module_list_KG.rds")
    colnames(cbm) <- c('Symbol','Module')
    cbm$Symbol <- as.character(cbm$Symbol)
    cbm$Module <- as.character(cbm$Module)
    df <- cbm
  }else if (type == 'tmod'){
    library(tmod)
    data(tmod)
    c2.indices <- tmod$MODULES2GENES
    c2.indices <- c2.indices[1:334]
    names(c2.indices) <- paste(names(c2.indices),tmod$MODULES$Title[1:334])
    names(c2.indices) <- gsub('LI.','',names(c2.indices))
    m <- matrix(nrow=0,ncol=2) #6617
    for (module in names(c2.indices)){
      genes <- c2.indices[[module]] # get the names for the genes in the module
      len <- length(genes) # get the module length
      name <- rep(module,len) # make the names
      matrix <- matrix(nrow=length(genes),ncol=2)
      matrix[,1] <- name 
      matrix[,2] <- genes
      m <- rbind(m,matrix)
    }
    df <- data.frame(m)
    colnames(df) <- c('Module','Symbol')
    c2.indices <- list
    m <- matrix(nrow=0,ncol=2) #6617
    for (module in names(c2.indices)){
      genes <- c2.indices[[module]] # get the names for the genes in the module
      len <- length(genes) # get the module length
      name <- rep(module,len) # make the names
      matrix <- matrix(nrow=length(genes),ncol=2)
      matrix[,1] <- name 
      matrix[,2] <- genes
      m <- rbind(m,matrix)
    }
    df <- data.frame(m)
    colnames(df) <- c('Module','Symbol')
  }else if (type == 'r4ra_drugtargetmodules'){
    df <- read.csv('./R4RA/5_modules/drug_target_modulesV2_0.6.csv')
    df <- df[,-1]
    colnames(df) <- c('Symbol','Module')
    df <- df[,c('Symbol','Module')]
    df$Module <- as.character(df$Module)
    df$Symbol <- as.character(df$Symbol)
  }else if (type == 'dataframe'){
    df <- list
    colnames(df) <- c('Symbol','Module')
    df <- df[,c('Symbol','Module')]
    df$Module <- as.character(df$Module)
    df$Symbol <- as.character(df$Symbol)
  }
  if (!is.na(data)){
    if (method == 'mean'){ # mean of VST data
      rows.to.keep2 <- subset(data, row.names(data) %in% as.character(df$Symbol)) # enter name of data
      rows.to.keep2$Symbol <- row.names(rows.to.keep2)
      merged <- merge(df, rows.to.keep2, by = 'Symbol')
      merged$Symbol <- NULL
      merged2 <- aggregate(. ~ Module, merged, mean)
      row.names(merged2) <- merged2$Module
      merged2$Module <- NULL
    }else if (method == 'meanthenzscore'){
      rows.to.keep2 <- subset(data, row.names(data) %in% as.character(df$Symbol)) # enter name of data
      rows.to.keep2$Symbol <- row.names(rows.to.keep2)
      merged <- merge(df, rows.to.keep2, by = 'Symbol')
      merged$Symbol <- NULL
      merged2 <- aggregate(. ~ Module, merged, mean)
      row.names(merged2) <- merged2$Module
      merged2$Module <- NULL
      merged2 <- data.frame(t(scale(t(merged2))))
    }else if (method == 'zmean'){ # calculate z score means
      data <- data.frame(t(scale(t(data)))) # take z score of genes before aggregating by mean
      rows.to.keep2 <- subset(data, row.names(data) %in% as.character(df$Symbol)) # enter name of data
      rows.to.keep2$Symbol <- row.names(rows.to.keep2)
      merged <- merge(df, rows.to.keep2, by = 'Symbol')
      merged$Symbol <- NULL
      merged$Module <- as.character(merged$Module) # recent addition
      merged2 <- aggregate(. ~ Module, merged, sum)
      row.names(merged2) <- merged2$Module
      merged2$Module <- NULL
      numberineach <- data.frame(table(merged$Module))
      numberineach$Freq <- sqrt(numberineach$Freq) # from GSVA vignette (Lee 2008 inferring pathway activity plos comp bio)
      merged2 <- merged2/numberineach$Freq # sum of z scores/ sqrt(gene set no.)
    }else if (method == 'zmean2'){ # calculate z score means w/o sqrt
      data <- data.frame(t(scale(t(data)))) # take z score of genes before aggregating by mean
      rows.to.keep2 <- subset(data, row.names(data) %in% as.character(df$Symbol)) # enter name of data
      rows.to.keep2$Symbol <- row.names(rows.to.keep2)
      merged <- merge(df, rows.to.keep2, by = 'Symbol')
      merged$Symbol <- NULL
      merged$Module <- as.character(merged$Module) # recent addition
      merged2 <- aggregate(. ~ Module, merged, sum)
      row.names(merged2) <- merged2$Module
      merged2$Module <- NULL
      numberineach <- data.frame(table(merged$Module))
      #numberineach$Freq <- (numberineach$Freq) # from GSVA vignette
      merged2 <- merged2/numberineach$Freq # sum of z scores/ sqrt(gene set no.)
    }else if (method == 'zpca'){ # calcuate first PC of z score
      data <- data.frame(t(scale(t(data)))) # take z score of genes before aggregating by mean
      rows.to.keep2 <- subset(data, row.names(data) %in% as.character(df$Symbol)) # enter name of data
      rows.to.keep2$Symbol <- row.names(rows.to.keep2)
      merged <- merge(df, rows.to.keep2, by = 'Symbol')
      merged$Symbol <- NULL
      # for module
      results <- matrix(ncol=ncol(merged)-1,nrow=0)
      i = 1
      for (module in unique(merged$Module)){
        mdata <- subset(merged, merged$Module == module)
        mdata <- mdata[,-1] # remove module name
        pca <- prcomp(t(mdata))
        scores <- pca$x
        results <- rbind(results,t(as.matrix(as.numeric(scores[,1]))))
        row.names(results)[i] <- module
        i = i + 1
      }
      colnames(results) <- row.names(scores)
      merged2 <- data.frame(results)
    }
  }else{
    merged2 <- df
  }
  return(merged2)
}