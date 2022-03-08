plotdrivers <- function(pcs,
         clin,
         block = NULL,
         unblock = NULL,
         kernel = NULL,
         kpar = NULL,
         top = NULL,
         n.pc = 5L,
         label = FALSE,
         alpha = 0.05,
         p.adj = NULL,
         title = 'Variation By Feature',
         legend = 'right',
         hover = FALSE) {

  pca = prcomp(t(data))
  pve <- seq_len(n.pc) %>% map_chr(function(pc) {
    p <- round(pca$sdev[pc]^2L / sum(pca$sdev^2L) * 100L, 1L)
    paste0('PC', pc, '\n(', round(p, 1L), '%)')
  }) 
  pcs = pca$x
  
  
  sig <- function(j, pc) {  # p-val fn
    mod <- lm(pcs[, pc] ~ clin[[j]])
    if_else(clin[[j]] %>% is.numeric, summary(mod)$coef[2L, 4L], 
            anova(mod)[1L, 5L])
  }
  
  # remove columns all identical or all different ###need this otherwise Age is removed as well
  #keep.cols = which(apply(clin, 2, function(x) {
  #  length(unique(x[! is.na(x)])) < nrow(clin) &
  #    length(unique(x[! is.na(x)])) > 1
  #}))
  #clin = clin[, keep.cols]
  
  
  df <- expand.grid(Feature = colnames(clin), PC = paste0('PC', seq_len(n.pc))) %>%
    rowwise(.) %>%
    dplyr::mutate(Association = sig(Feature, PC)) %>%   # Populate
    ungroup(.)
  if (!p.adj %>% is.null) {
    df <- df %>% mutate(Association = p.adjust(Association, method = p.adj))
  }
  df <- df %>% 
    mutate(Significant = if_else(Association <= alpha, TRUE, FALSE), 
           Association = -log(Association))
  
  
  # Build plot
  if (!p.adj %>% is.null && p.adj %in% c('holm', 'hochberg', 'hommel', 'bonferroni', 'BH', 'BY', 'fdr')) {
    leg_lab <- expression(~-log(italic(q)))
  } else {
    leg_lab <- expression(~-log(italic(p)))
  }
  p <- ggplot(df, aes(PC, Feature, fill = Association, text = Association,
                      color = Significant)) +
    geom_tile(size = 1L, width = 0.9, height = 0.9) +
    coord_equal() +
    scale_fill_gradientn(
      colors = c('white', 'pink', 'orange', 'red', 'darkred'), 
      name = leg_lab) +
    scale_color_manual(values = c('grey90', 'black')) +
    guides(color = FALSE) +
    #scale_x_discrete(labels = pve) +
    labs(title = title, x = 'Principal Component', y='Feature') +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5), 
          axis.text.x = element_text(hjust = 0.5))
  if (label) {
    p <- p + geom_text(aes(label = round(Association, 2L)), size =3)
  }
  return(p)
}

