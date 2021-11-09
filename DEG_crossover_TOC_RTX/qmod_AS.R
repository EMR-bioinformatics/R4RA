#' Run QuSAGE on limma or DESeq2 model objects
#'
#' This function is a wrapper for the QuSAGE algorithm, which tests for pathway
#' enrichment, designed for easy integration with limma and DESeq2.
#'
#' @param fit An object of class \code{limma::\link[limma]{MArrayLM}}, as
#'   created by a call to \code{\link[limma]{eBayes}}, or a \code{
#'   \link[DESeq2]{DESeqDataSet}} that has been fit with a negative binomial
#'   GLM. See Details.
#' @param dat An expression matrix or matrix-like object, with rows
#'   corresponding to probes and columns to samples. Only necessary if \code{
#'   fit} is an \code{MArrayLM} object.
#' @param filter Numeric vector of length 2 specifying the filter criterion.
#'   Each probe must have at least \code{filter[1]} log2-counts per million in
#'   at least \code{filter[2]} libraries to pass the expression threshold. Only
#'   relevant if \code{fit} is a \code{DESeqDataSet}, in which case the
#'   normality of transformed residuals at various values of \code{filter}
#'   should be checked prior to running \code{qmod}, for example using \code{
#'   \link{check_resid}}. See Details.
#' @param coef Column name or number specifying which coefficient of the model
#'   is of interest. Alternatively, a vector of three or more such strings or
#'   numbers, in which case pathways are ranked by the \emph{F}-statistic for
#'   that set of coefficients.
#' @param contrast Character or numeric vector of length two, specifying the
#'   column names or numbers to be contrasted. The first and second elements
#'   will be the numerator and denominator, respectively, of the fold change
#'   calculation. Exactly one of \code{coef} or \code{contrast} must be \code{
#'   NULL}.
#' @param geneSets A named list of one or several gene sets.
#' @param n.points The number of points at which to sample the convoluted
#'   \emph{t}-distribution. See Details.
#'
#' @details
#' QuSAGE
#'
#' ### INTRO TO QMOD? ###
#' \code{qmod} combines the
#'
#'
#' the statistical flexibility and empirical Bayes methods of
#' \code{limma} with the specificity and sensitivity of the QuSAGE algorithm for
#' detecting pathway enrichment. Simulations have shown that this pipeline
#' outperforms each package's independent methods for gene set analysis in
#' complex experimental designs, namely \code{\link[limma]{camera}} and
#' \code{link[qusage]{qgen}}. See Watson & John, forthcoming.
#'
#' ### SOMETHING ON MArrayLM vs. DESeqDataSet OBJECTS ###
#'
#' If fit is a voom object, make sure to run lcpm(dat).
#'
#' By default \code{n.points} is set to \code{2^14}, or 16,384 points, which will give very
#' accurate \emph{p}-values in most cases. Sampling at more points will increase the
#' accuracy of the resulting \emph{p}-values, but will also linearly increase the
#' amount of time needed to calculate the result. With larger sample sizes, as few
#' as 1/4 this number of points can be used without seriously affecting the accuracy
#' of the resulting \emph{p}-values. However, when there are a small number of samples
#' (i.e., fewer than 8 samples total), the \emph{t}-distribution must be sampled over
#' a much wider range, and the number of points needed for sampling should be
#' increased accordingly. It is recommended that when running \code{qmod} with fewer
#' than 8 samples, the number of points be increased to at least 2^15 or 2^16. It may
#' also be useful to test higher values of this parameter, as it can often result in
#' a much more significant \emph{p}-value with small sample sizes.
#'
#' @return A data frame with the following columns:
#' \itemize{
#'   \item{'Pathway'} Name of the pathway.
#'   \item{'logFC'} Average log2 fold change of genes in the pathway.
#'   \item{'p.value'} The gene set \emph{p}-value, as calculated using
#'     \code{pdf.pVal}.
#'   \item{'FDR'} The false discovery rate, as estimated using the Benjamini-Hochberg
#'     algorithm.
#' }
#'
#' @references
#' Yaari, G. Bolen, C.R., Thakar, J. & Kleinstein, S.H. (2013).
#' \href{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3794608/}{Quantitative set
#' analysis for gene expression: a method to quantify gene set differential
#' expression including gene-gene correlations}. \emph{Nucleic Acids Res.},
#' \emph{41}(18): e170.
#'
#' Turner, J.A., Bolen, C.R. & Blankenship, D.M. (2015).
#' \href{http://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-015-0707-9}{
#' Quantitative gene set analysis generalized for repeated measures, confounder
#' adjustment, and continuous covariates}. \emph{BMC Bioinformatics}, \emph{
#' 16}(1): 272.
#'
#' Smyth, G.K. (2004).
#' \href{http://www.statsci.org/smyth/pubs/ebayes.pdf}{Linear models and
#' empirical Bayes methods for assessing differential expression in microarray
#' experiments}. \emph{Stat. Appl. Genet. Molec. Biol.}, \emph{3}(1).
#'
#' Love, M., Huber, W., & Anders, S. (2014).
#' \href{https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0550-8}{
#' Moderated estimation of fold change and dispersion for RNA-seq data with
#' DESeq2}. \emph{Genome Biology}, \strong{15}:550.
#'
#' @examples
#' # Simulate data, fit limma model
#' library(limma)
#' eset <- matrix(rnorm(5000 * 10), nrow = 5000, ncol = 10,
#'                dimnames = list(seq_len(5000), paste0('S', seq_len(10))))
#' clin <- data.frame(treat = rep(c("ctl", "trt"), each = 5),
#'                    batch = rep(c("b1", "b2"), times = 5),
#'                       x1 = rnorm(10),
#'                       x2 = runif(10))
#' des <- model.matrix(~ treat + batch + x1 + x2, data = clin)
#' fit <- eBayes(lmFit(eset, des))
#'
#' # Create list of differentially expressed pathways
#' geneSets <- list()
#' for (i in 0:10) {
#'   genes <- ((30 * i) + 1):(30 * (i + 1))
#'   eset[genes, clin$treat == "trt"] <- eset[genes, clin$treat == "trt"] + rnorm(1)
#'   geneSets[[paste("Set", i)]] <- genes
#' }
#'
#' # Run qmod
#' top <- qmod(fit, dat = eset, coef = 2, geneSets = geneSets)
#'
#' @export
#' @importFrom limma eBayes getEAWP
#' @importFrom DESeq2 sizeFactors normalizationFactors counts results
#' @importFrom SummarizedExperiment assays
#' @importFrom edgeR DGEList calcNormFactors cpm
#' @import qusage
#' @import dplyr
#'

qmod <- function(fit,
                 dat = NULL,
                 filter = NULL,
                 coef,
                 contrast = NULL,
                 geneSets,
                 n.points = 2^14) {
  
  # Preliminaries
  if (nrow(fit) < 3L) {
    stop('fit must have at least three probes.')
  }
  if (is(fit, 'MArrayLM')) {
    coefs <- colnames(fit$coefficients)
    p <- length(coefs)
    if (is.null(dat)) {
      stop('dat must be provided when fit is an MArrayLM object.')
    }
    if (nrow(fit) != nrow(dat) | nrow(fit$design) != ncol(dat)) {
      stop('dat is not conformal with fit.')
    }
    if (!identical(rownames(fit), rownames(dat))) {
      stop('dat and fit must have identical rownames.')
    }
    if (!is.null(filter)) {
      warning('filter is ignored when fit is an MArrayLM object.')
    }
  } else if (is(fit, 'DESeqDataSet')) {
    coefs <- resultsNames(fit)
    p <- ncol(model.matrix(design(fit), colData(fit)))
    if (is.null(rownames(fit))) {
      rownames(fit) <- seq_len(nrow(fit))
    }
    if (!is.null(dat)) {
      warning('dat is ignored when fit is a DESeqDataSet.')
    }
    if (is.null(filter)) {
      stop('filter must be supplied when fit is a DESeqDataSet. See ',
           '?check_resid.')
    } else if (length(filter) != 2L) {
      stop('filter must be a vector of length 2.')
    }
  } else {
    stop('fit must be an object of class MArrayLM or DESeqDataSet.')
  }
  if (is.null(contrast)) {
    if (is.character(coef) & !coef %in% coefs) {
      stop(paste0("'", coef, "' not found in fit's design matrix."))
    }
    if (is.numeric(coef) & !coef %in% seq_len(p)) {
      stop(paste("No coef number", coef, "found in fit's design matrix."))
    }
  } else if (is.null(coef)) {
    if (length(contrast) != 2L) {
      stop('contrast must be a vector of length 2.')
    }
    if ((is.character(contrast) & any(!contrast %in% coefs)) ||
        (is.numeric(contrast) & any(!contrast %in% seq_len(p)))) {
      stop("Both coefficients passed to contrast must be in fit's design ',
           'matrix.")
    }
  } else {
    stop('Exactly one of coef or contrast must be NULL.')
  }
  if (is(fit, 'MArrayLM')) {
    if (is.null(fit$t) & is.null(fit$F) & is.null(contrast)) {
      fit <- eBayes(fit)
      warning('Standard errors for fit have not been moderated. Running ',
              'eBayes before testing for enrichment. See ?eBayes for more ',
              'info.')
    }
    if (!is.null(fit$t) & !is.null(fit$F) & is.null(coef)) {
      stop('Standard errors for fit must not be moderated when passing a ',
           'contrast to qmod. Use an lmFit output instead. The function will ',
           'internally create the appropriate contrast matrix and run eBayes ',
           'on that. See ?contrasts.fit for more info.')
    }
  }
  if (!is.list(geneSets)) {
    stop('geneSets must be a list.')
  }
  if (is.null(names(geneSets))) {
    stop('geneSets must be a named list.')
  }
  overlap <- sapply(geneSets, function(g) sum(g %in% rownames(fit)))
  geneSets <- geneSets[overlap > 1L]
  if (length(geneSets) == 0L) {
    stop('No overlap detected between the genes in fit and those in geneSets.')
  }
  
  # Prep data
  if (is(fit, 'MArrayLM')) {
    dat <- getEAWP(dat)$exprs
    n <- ncol(dat)
    resid_mat <- residuals(fit, dat)
    if (is.null(coef)) {
      coef <- 'Contrast'
      cm <- makeContrasts(contrasts = paste(contrast[1], '-', contrast[2]), 
                          levels = coefs)
      colnames(cm) <- coef
      fit <- contrasts.fit(fit, cm)
      fit <- eBayes(fit)
    }
    mean <- fit$coefficients[, coef]
    SD <- sqrt(fit$s2.post) * fit$stdev.unscaled[, coef]
    sd.alpha <- SD / (fit$sigma * fit$stdev.unscaled[, coef])
    sd.alpha[is.infinite(sd.alpha)] <- 1L
    dof <- fit$df.total
  } else {
    cnts <- counts(fit)
    n <- ncol(cnts)
    keep <- rowSums(cpm(cnts) > filter[1L]) >= filter[2L]
    fit <- fit[keep, , drop = FALSE]
    cnts <- fit %>%
      counts(normalized = TRUE) %>%
      cpm(log = TRUE, prior.count = 1L)
    if (sizeFactors(fit) %>% is.null) {
      signal_mat <- assays(fit)[['mu']] / normalizationFactors(fit)
    } else {
      signal_mat <- t(t(assays(fit)[['mu']]) / sizeFactors(fit))
    }
    signal_mat <- cpm(signal_mat, log = TRUE, prior.count = 1L)
    resid_mat <- cnts - signal_mat
    if (contrast %>% is.null) {
      dds_res <- results(fit, name = resultsNames(fit)[coef], independentFiltering = FALSE)
    } else {
      dds_res <- results(fit, contrast = list(contrast),
                         independentFiltering = FALSE)
    }
    mean <- dds_res$log2FoldChange
    SD <- dds_res$lfcSE
    sd.alpha <- rep(1L, times = nrow(fit))
    dof <- rep((ncol(fit) - p), times = nrow(fit))
  }
  names(mean) <- names(SD) <- names(sd.alpha) <- names(dof) <- rownames(fit)
  
  # Run QuSAGE functions
  res <- newQSarray(mean = mean, SD = SD, sd.alpha = sd.alpha, dof = dof,
                    labels = rep('resid', n))          # Create QSarray obj
  res <- aggregateGeneSet(res, geneSets, n.points)     # PDF per gene set
  res <- calcVIF(resid_mat, res, useCAMERA = FALSE)    # VIF on resid_mat
  
  # Export
  output <- as.data.frame(qsTable(res, number = Inf, sort.by = 'p')) 
  myqs <- qvalue::qvalue(output$p.Value)
  output$qval <- myqs$qvalues
   return(output)
  
}


# Extend to ANOVA F-tests/likelihood ratio tests?