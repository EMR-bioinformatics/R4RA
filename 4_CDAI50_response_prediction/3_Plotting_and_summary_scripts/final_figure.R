library(ggplot2)
library(ggpubr)

load("/media/gcpeac/R4RA_machine_learning/Data/best_models5.rdata")

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

theme_tr =  theme(plot.title = element_text(hjust = 0.5, size=14), 
                  axis.title = element_text(size=12))


############################
# Secondary timepoint
############################

toc_mod = toc_other_filt_plots[[1]]
rtx_mod = rtx_other_filt_plots[[5]]
any_mod = any_filt_plots

toc_mod$data =   toc_mod$data[! is.na(  toc_mod$data$Overall), ]

pl1 = list(
  # titles
  bio_facet(text = "", colour="white", angle=0), 
  bio_facet(text = "Rituximab", angle=0, size=7), 
  bio_facet(text = "Tocilizumab", angle=0, size=7), 
  bio_facet(text = "Resistant", angle=0, size=7),
  
  # Plots  
  bio_facet(text = "Test ROC", size=7), 
  rtx_mod[[1]] + labs(title=names[rtx_mod$fit$method]) + theme_tr, 
  toc_mod[[1]] + labs(title=names[toc_mod$fit$method]) + theme_tr, 
  any_mod[[1]] + labs(title=names[any_mod$fit$method]) + theme_tr, 
  bio_facet(text = "Left-Out ROC", size=7), 
  rtx_mod[[2]] + labs(title=NULL), 
  toc_mod[[2]] + labs(title=NULL), 
  any_mod[[2]] + labs(title=NULL), 
  bio_facet(text = "Variable Importance", size=7), 
  rtx_mod[[3]] + labs(title=NULL) + theme(axis.text = element_text(size=10)), 
  toc_mod[[3]] + labs(title=NULL) + theme(axis.text = element_text(size=10)), 
  any_mod[[3]] + labs(title=NULL) + theme(axis.text = element_text(size=10))
)


cairo_pdf("/media/gcpeac/R4RA_machine_learning/Data/final_fig.pdf", 
          height=0.8*1.7*11.75, width=1.7*8.25)
ggarrange(plotlist = pl1, ncol=4, nrow=4, #align="hv",
          widths=c(0.15, 1, 1, 1), heights=c(0.15,  1.2, 1, 1.8))
dev.off()

# pl_rtx = lapply(rtx_other_filt_plots, function(x) {
#   vi = x$variable_importance + theme(axis.text.y = element_text(size=6))
#   ggarrange( bio_facet(text = x$fit$method, colour="grey60", angle=90), 
#              x$test_roc, x$lo_roc, vi + 
#                labs(title=paste("n features=", nrow(x$variable_importance$data))), 
#              ncol=4, widths=c(0.1, 1, 1, 1))
# })
# 
# cairo_pdf("/media/gcpeac/R4RA_machine_learning/Data/final_fig_rtx_models.pdf", 
#           height=2*0.8*1.7*11.75, width=2*8.25)
# ggarrange(plotlist=pl_rtx, nrow=length(rtx_other_filt_plots), ncol=1)
# dev.off()

############################





############################
# Primary timepoint
############################
rtx_filt_plots_bl[[3]]$data = 
  rtx_filt_plots_bl[[3]]$data[! is.na(  rtx_filt_plots_bl[[3]]$data$Overall), ]


pl2 = list(
  # titles
  bio_facet(text = "", colour="white", angle=0), 
  bio_facet(text = "Rituximab", angle=0, size=7), 
  bio_facet(text = "Tocilizumab", angle=0, size=7),
  
  # Plots  
  bio_facet(text = "Test ROC", size=7), 
  rtx_filt_plots_bl[[1]] + labs(title=names[rtx_filt_plots_bl$fit$method]) + 
    theme_tr, 
  toc_filt_plots_bl[[1]] + labs(title=names[toc_filt_plots_bl$fit$method]) + 
    theme_tr,  
  bio_facet(text = "Left-Out ROC", size=7), 
  rtx_filt_plots_bl[[2]] + labs(title=NULL), 
  toc_filt_plots_bl[[2]] + labs(title=NULL), 
  
  bio_facet(text = "Variable Importance", size=7), 
  rtx_filt_plots_bl[[3]] + labs(title=NULL) + theme(axis.text = element_text(size=10)), 
  toc_filt_plots_bl[[3]] + labs(title=NULL) + theme(axis.text = element_text(size=10))#, 
)


cairo_pdf("/media/gcpeac/R4RA_machine_learning/Data/final_fig_bl.pdf", height=0.8*1.7*11.75, width=3*1.7*8.25/4)
ggarrange(plotlist = pl2, ncol=3, nrow=4, #align="hv",
          widths=c(0.15*3/4, 1*3/4, 1*3/4), heights=c(0.15,  1.2, 1, 1.8))
dev.off()
############################



############################
# primary response for rtx and toc, secondary for refractory
############################

toc_mod = toc_bl_other_filt_plots[[1]]
rtx_mod = rtx_filt_plots_bl
any_mod = any_filt_plots

pl2 = list(
  # titles
  bio_facet(text = "", colour="white", angle=0), 
  bio_facet(text = "Rituximab", angle=0, size=7), 
  bio_facet(text = "Tocilizumab", angle=0, size=7), 
  bio_facet(text = "Resistant", angle=0, size=7),
  
  # Plots  
  bio_facet(text = "Test ROC", size=7), 
  rtx_mod[[1]] + labs(title=names[rtx_mod$fit$method]) + theme_tr, 
  toc_mod[[1]] + labs(title=names[toc_mod$fit$method]) + theme_tr, 
  any_mod[[1]] + labs(title=names[any_mod$fit$method]) + theme_tr, 
  bio_facet(text = "Left-Out ROC", size=7), 
  rtx_mod[[2]] + labs(title=NULL), 
  toc_mod[[2]] + labs(title=NULL), 
  any_mod[[2]] + labs(title=NULL), 
  bio_facet(text = "Variable Importance", size=7), 
  rtx_mod[[3]] + labs(title=NULL) + theme(axis.text = element_text(size=10)), 
  toc_mod[[3]] + labs(title=NULL) + theme(axis.text = element_text(size=10)), 
  any_mod[[3]] + labs(title=NULL) + theme(axis.text = element_text(size=10))
)

w = 1.7*8.25
h = 0.8*1.7*11.75
hs = c(0.15, 1, 1, 1.8)
hs = hs/sum(hs)
y = hs[1]*h/w
ws = c(y, 1/(y+3), 1/(y+3), 1/(y+3))

cairo_pdf("/media/gcpeac/R4RA_machine_learning/Data/final_fig_all.pdf", 
          height=h, width=w)
ggarrange(plotlist = pl2, ncol=4, nrow=4, #align="hv",
          widths=ws, heights=hs)
dev.off()
