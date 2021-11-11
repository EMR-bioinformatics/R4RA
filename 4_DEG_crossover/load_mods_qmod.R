loadmodsforqusage <- function(load=module){ # tmod, wgcna
  if (load == 'LI'){
    library(tmod)
    data(tmod)
    c2.indices <- tmod$MODULES2GENES
    c2.indices <- c2.indices[1:334]
    names(c2.indices) <- paste(names(c2.indices),tmod$MODULES$Title[1:334])
    names(c2.indices) <- gsub('LI.','',names(c2.indices))
    
  }else if (load == 'LI_reduced'){
    library(tmod)
    data(tmod)
    c2.indices <- tmod$MODULES2GENES
    c2.indices <- c2.indices[1:334]
    names(c2.indices) <- paste(names(c2.indices),tmod$MODULES$Title[1:334])
    names(c2.indices) <- gsub('LI.','',names(c2.indices))
    TBA <- grep('TBA', names(c2.indices))
    c2.indices <- c2.indices[-TBA]
    
  }else if (load == 'wgcna'){
    # wgcna
    df <- read.csv('~/PEAC/8_modularframework/synovium_baseline_modules_annotatedv11.csv')
    df <- df[,-1]
    final_list <- list()
    modules <- unique(df$Module)
    for (module in modules){
      tempdf <- subset(df, Module == module)
      new_list <- list(as.character(tempdf$Symbol))
      names(new_list)[1]<-module
      final_list <- c(final_list, new_list)
    }
    c2.indices <- final_list
  } else if (load == 'reduced_wgcna'){
    df <- read.csv('~/PEAC/reduced_wgcna_modules_V1.csv')
    final_list <- list()
    modules <- unique(df$Module)
    for (module in modules){
      tempdf <- subset(df, Module == module)
      new_list <- list(as.character(tempdf$Symbol))
      names(new_list)[1]<-module
      final_list <- c(final_list, new_list)
    }
    c2.indices <- final_list

  }
  return(c2.indices)
}