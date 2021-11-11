flexibleboxplot <- function(meta,variable,toadd,xlab='no',fill=10,theme='bw',colours=FALSE,
                            title='givetitle',rotatex=FALSE,fsize=12.5,simple=FALSE){ 
  
  # changed fsize from 13 to 12.5
  
  # meta = input data with gene/module to add
  # variable = variable of interest (e.g. EULAR3)
  # toadd = name of gene/module
  
  # debugging
  #meta <- annonx
  #variable <- 'consensuscluster'
  #toadd <- 'MS4A1'
  
  library(ggplot2)
  
  if (colours == FALSE){
    chrispal1 <- c('cornflowerblue', 'gold1', 'darkorchid', 'skyblue1', 'plum1', 'violetred', 'forestgreen', 'midnightblue')
  }else{
    chrispal1 <- colours
  }
  
  # initial reformatting
  # remove NAs from input data
  meta[[variable]] <- as.character(meta[[variable]]) # new
  
  meta2 <- meta[!is.na(meta[[variable]]),]
  meta2 <- meta2[!(meta2[[variable]]=='Unknown'),]
  # get range of variable
  rnge <- range(meta[[toadd]])
  toaddd <- ((rnge[2]-rnge[1])/100)*fill
  
  # plotting code
  if (theme == 'bw' && rotatex == FALSE){ # use black and white theme
    p <- ggplot(meta2, aes_string(x = variable, y = toadd)) +
      geom_boxplot(aes_string(fill = variable), lwd=0.5,colour = "black", outlier.shape = NA) +
      geom_point(aes_string(fill = variable), size = 2, shape = 21, position = position_jitterdodge()) +
      theme_bw() +
      labs(title = toadd) +
      scale_fill_manual(values=c(chrispal1[1], chrispal1[2], chrispal1[3], chrispal1[4], chrispal1[5], chrispal1[6])) +
      theme(plot.title = element_text(hjust = 0.5),
            #panel.border = element_rect(size=1, color = 'black'), ## commented out
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_rect(colour = "white", fill = "white"),
            axis.text.y = element_text(size = fsize, colour = 'black'),
            axis.text.x =  element_text(size = fsize, colour = 'black'),
            axis.title.x = element_text(size = fsize, colour = 'black'),
            axis.title.y = element_text(size = fsize, colour = 'black'),
            legend.text = element_text(size = fsize, colour = 'black'),
            legend.title=element_blank()
      ) + ggtitle(title) + scale_y_continuous(limits=c(rnge[1]-toaddd,rnge[2]+toaddd)) +
      ylab('Normalised expression') +
      if (xlab == 'no'){
        xlab(variable)
      }else if (xlab != 'no'){
        xlab(xlab)
      }
  }else if (theme == 'bw' && rotatex == TRUE){
    p <- ggplot(meta2, aes_string(x = variable, y = toadd)) +
      geom_boxplot(aes_string(fill = variable), lwd=0.5,colour = "black", outlier.shape = NA) +
      geom_point(aes_string(fill = variable), size = 2, shape = 21, position = position_jitterdodge()) +
      theme_bw() +
      labs(title = toadd) +
      scale_fill_manual(values=c(chrispal1[1], chrispal1[2], chrispal1[3], chrispal1[4], chrispal1[5], chrispal1[6])) +
      theme(plot.title = element_text(hjust = 0.5),
            #panel.border = element_rect(size=1, color = 'black'), ## commented out
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_rect(colour = "white", fill = "white"),
            axis.text.y = element_text(size = fsize, colour = 'black'),
            axis.text.x =  element_text(size = fsize, colour = 'black', angle = 45, vjust = 1, hjust=1),
            axis.title.x = element_text(size = fsize, colour = 'black'),
            axis.title.y = element_text(size = fsize, colour = 'black'),
            legend.text = element_text(size = fsize, colour = 'black'),
            legend.title=element_blank()
      ) + ggtitle(title) + scale_y_continuous(limits=c(rnge[1]-toaddd,rnge[2]+toaddd)) +
      ylab('Normalised expression') +
      if (xlab == 'no'){
        xlab(variable)
      }else if (xlab != 'no'){
        xlab(xlab)
      }
  }
  if (simple){
    p <- ggplot(meta2, aes_string(x = variable, y = toadd)) +
      geom_boxplot(aes_string(fill = variable), lwd=0.5,colour = "black", outlier.shape = NA) +
      geom_point(aes_string(fill = variable), size = 2, shape = 21, position = position_jitterdodge()) +
      theme_bw() +
      labs(title = toadd) +
      scale_fill_manual(values=c(chrispal1[1], chrispal1[2], chrispal1[3], chrispal1[4], chrispal1[5], chrispal1[6])) +
      theme(plot.title = element_text(hjust = 0.5),
            #panel.border = element_rect(size=1, color = 'black'), ## commented out
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_rect(colour = "white", fill = "white"),
            axis.text.y = element_text(size = fsize, colour = 'black'),
            axis.text.x =  element_blank(),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            legend.text = element_blank(),
            legend.title = element_blank(), legend.position="none"
      ) + ggtitle(title) + scale_y_continuous(limits=c(rnge[1]-toaddd,rnge[2]+toaddd)) +
      ylab('Normalised expression') +
      if (xlab == 'no'){
        xlab(variable)
      }else if (xlab != 'no'){
        xlab(xlab)
      }
  }
  #
  print(p)
  return(p)
}
  

