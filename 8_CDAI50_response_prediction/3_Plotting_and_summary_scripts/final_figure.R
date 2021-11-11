library(ggplot2)
library(ggpubr)

load("/media/gcpeac/R4RA_machine_learning/Data/best_models7.rdata")

names = c("gbm"="Gradient Boosting \nMachine", "rf"=" \nRandom forest", 
          "svmRadial"="\nSupport Vector Machine (Radial)",
          "svmPoly"="\nSupport Vector Machine (Polynomial)",
          "glmnet"="Lasso and Elastic-Net \nGeneralized Linear Model")

bio_facet = function(text, colour="grey80", outline_colour="white", outline_width=5,  
                     font_colour="black", angle=90, align="center", size=5){
  p = ggplot(df=NULL)  + 
    theme_void() + 
    annotate("text", x=0, y=1, label=text, angle=angle, colour=font_colour, size=size) + 
    theme(panel.background=element_rect(fill=colour, color=outline_colour, size=outline_width))
  if(align == "left") p = p + lims(x=c(0, 1))
  if(align == "right") p = p + lims(x=c(-1, 0))
  return(p)
}

facet_grid <- function(rtx_mod, toc_mod, any_mod, w = 1.7*8.25, h = 0.8*1.7*11.75){
  plot_list = list(
    # titles
    bio_facet(text = "", colour="white", angle=0), 
    bio_facet(text = "Rituximab", angle=0, size=7), 
    bio_facet(text = "Tocilizumab", angle=0, size=7), 
    bio_facet(text = "Refractory", angle=0, size=7),
    
    # Plots  
    bio_facet(text = "Test ROC", size=7), 
    rtx_mod[[1]] + labs(title=names[rtx_mod$fit$method]) + theme_tr, 
    toc_mod[[1]] + labs(title=names[toc_mod$fit$method]) + theme_tr, 
    any_mod[[1]] + labs(title=names[any_mod$fit$method]) + theme_tr, 
    bio_facet(text = "Left-Out Inner Fold ROC", size=7), 
    rtx_mod[[2]] + labs(title=NULL), 
    toc_mod[[2]] + labs(title=NULL), 
    any_mod[[2]] + labs(title=NULL), 
    bio_facet(text = "Variable Importance", size=7), 
    rtx_mod[[3]] + labs(title=NULL) + theme(axis.text = element_text(size=10)), 
    toc_mod[[3]] + labs(title=NULL) + theme(axis.text = element_text(size=10)), 
    any_mod[[3]] + labs(title=NULL) + theme(axis.text = element_text(size=10))
  )
  
  hs = c(0.15, 1, 1, 1.8)
  hs = hs/sum(hs)
  y = hs[1]*h/w
  ws = c(y, 1/(y+3), 1/(y+3), 1/(y+3))
  
  ggarrange(plotlist = plot_list, ncol=4, nrow=4, widths=ws, heights=hs)
}

theme_tr =  theme(plot.title = element_text(hjust = 0.5, size=14), 
                  axis.title = element_text(size=12))



# Secondary timepoint
facet_grid(rtx_other_filt_plots[[5]], toc_other_filt_plots[[1]], any_other_filt_plots[[1]])
#facet_grid(rtx_other_filt_plots[[5]], toc_filt_plots, any_filt_plots)
cairo_pdf(paste0("/media/gcpeac/R4RA_machine_learning/Data/final_fig_secondary_resp_",  Sys.Date(), ".pdf"), 
          height=0.8*1.7*11.75, width=1.7*8.25)
facet_grid(rtx_other_filt_plots[[5]], toc_other_filt_plots[[1]], any_other_filt_plots[[1]])
dev.off()


# Primary timepoint
#load("/media/gcpeac/R4RA_machine_learning/Data/best_models4.rdata")
facet_grid(rtx_bl_other_filt_plots[[1]], toc_bl_other_filt_plots[[1]], any_other_filt_plots[[3]])
#facet_grid(rtx_filt_plots_bl, toc_filt_plots_bl2, any_filt_plots)
cairo_pdf(paste0("/media/gcpeac/R4RA_machine_learning/Data/final_fig_primary_resp_",  Sys.Date(), "_fix.pdf"), 
          height=0.8*1.7*11.75, width=1.7*8.25)
facet_grid(rtx_bl_other_filt_plots[[1]], toc_bl_other_filt_plots[[1]], any_other_filt_plots[[1]])
dev.off()


load("/media/gcpeac/R4RA_machine_learning/Data/best_models4.rdata")
length(intersect(rtx_bl_other_filt_plots[[1]]$vi$gene[rtx_bl_other_filt_plots[[1]]$vi$Overall != 0], 
          rtx_filt_plots_bl$vi$gene[rtx_filt_plots_bl$vi$Overall != 0]))/
  length(rtx_bl_other_filt_plots[[1]]$vi$gene[rtx_bl_other_filt_plots[[1]]$vi$Overall != 0])*100

rtx_bl_other_filt_plots[[1]]$vi$gene[rtx_bl_other_filt_plots[[1]]$vi$Overall != 0][! rtx_bl_other_filt_plots[[1]]$vi$gene[rtx_bl_other_filt_plots[[1]]$vi$Overall != 0] %in% rtx_filt_plots_bl$vi$gene[rtx_filt_plots_bl$vi$Overall != 0]]


# overlap
intersect(rtx_bl_other_filt_plots[[1]]$vi$gene[rtx_bl_other_filt_plots[[1]]$vi$Overall != 0], 
          rtx_filt_plots_bl$vi$gene[rtx_filt_plots_bl$vi$Overall != 0])

# final fit genes not in initial fig
setdiff(rtx_bl_other_filt_plots[[1]]$vi$gene[rtx_bl_other_filt_plots[[1]]$vi$Overall != 0], 
          rtx_filt_plots_bl$vi$gene[rtx_filt_plots_bl$vi$Overall != 0])

# initial fig not in final fit genes
setdiff(rtx_filt_plots_bl$vi$gene[rtx_filt_plots_bl$vi$Overall != 0], 
        rtx_bl_other_filt_plots[[1]]$vi$gene[rtx_bl_other_filt_plots[[1]]$vi$Overall != 0])



# overlap
intersect(toc_bl_other_filt_plots[[1]]$vi$gene[toc_bl_other_filt_plots[[1]]$vi$Overall != 0], 
          toc_filt_plots_bl$vi$gene[toc_filt_plots_bl$vi$Overall != 0])

# initial fig not in final fit genes 
setdiff(toc_bl_other_filt_plots[[1]]$vi$gene[toc_bl_other_filt_plots[[1]]$vi$Overall != 0], 
        toc_filt_plots_bl$vi$gene[toc_filt_plots_bl$vi$Overall != 0])

# final fit genes not in initial fig
setdiff(toc_filt_plots_bl$vi$gene[toc_filt_plots_bl$vi$Overall != 0], 
        toc_bl_other_filt_plots[[1]]$vi$gene[toc_bl_other_filt_plots[[1]]$vi$Overall != 0])




length(intersect(toc_bl_other_filt_plots[[1]]$vi$gene[toc_bl_other_filt_plots[[1]]$vi$Overall != 0], 
                 toc_filt_plots_bl2$vi$gene[toc_filt_plots_bl2$vi$Overall != 0]))/
  length(toc_bl_other_filt_plots[[1]]$vi$gene[toc_bl_other_filt_plots[[1]]$vi$Overall != 0])*100

length(intersect(any_other_filt_plots[[1]]$vi$gene[any_other_filt_plots[[1]]$vi$Overall != 0], 
                 any_filt_plots$vi$gene[any_filt_plots$vi$Overall != 0]))/
  length(any_other_filt_plots[[1]]$vi$gene[any_other_filt_plots[[1]]$vi$Overall != 0])*100

toc_model = toc_bl01
rtx_model = rtx_bl01
ref_model = any01

toc_plots = toc_bl_other_filt_plots[[1]]
rtx_plots = rtx_bl_other_filt_plots[[1]]
ref_plots = any_other_filt_plots

save(toc_model, rtx_model, ref_model, toc_plots, rtx_plots, ref_plots, 
     file=paste0("/media/gcpeac/R4RA_machine_learning/Data/final_models_",  Sys.Date(), ".rdata"))



intersect(toc_filt_plots_bl2$vi$gene, rtx_filt_plots_bl$vi$gene)

temp = toc_filt_plots_bl2$vi[toc_filt_plots_bl2$vi$gene %in% 
                               intersect(toc_filt_plots_bl2$vi$gene, rtx_filt_plots_bl$vi$gene), ]


load("/media/gcpeac/R4RA_machine_learning/Data/curated_input.rdata")


corv = cor(t(vst))
rownames(corv) = gsub("\\-", ".", rownames(corv))
corv = corv[as.character(temp2$Gene[temp2$description == "" & temp2$Type == ""]), ]
corv = corv[, apply(corv, 2, function(x) any(x >= 0.9))]
corv = corv[, ! is.na(colnames(corv))]
corv = corv[, ! gsub("\\-", ".", colnames(corv)) %in% rownames(corv)]


library(ComplexHeatmap)

Heatmap(t(corv))


temp3 = gene_summary(colnames(corv), disease_api_token = api_key)
