#' High-dimensional Mediation Analysis
#' 
#' \code{hima_rec} is used to estimate and test high-dimensional mediation effects
#' by recursively using [hima]. 
#' 
#' The most significant markers are added as covariables in order to test the 
#' other markers. This has the benefit to remove markers that are associated to 
#' a true causal, without being causal for the mediation. It also makes room for
#' other marker with smaller effects. This procedure stops when there are no 
#' more significant markers that are detected.
#' 
#' @inheritParams hima
#' @param maxiter Number of maximum iterations.
#' @param threshold Threshold on p-values to be considered as significant.
#'
#' @return The list of markers that are considered significant.
#' @export
#'
hima_rec <- function(X, Y, M, COV = NULL, 
                     family = c("gaussian", "binomial"), 
                     penalty = c("MCP", "SCAD", "lasso"), 
                     topN = NULL, 
                     parallel = TRUE, 
                     ncore = NULL, 
                     verbose = TRUE, 
                     maxiter = 5,
                     threshold = 0.05,
                     keep.all.first = TRUE,
                     ...) {
  
  # retrieve the arguments
  args <- as.list(environment())
  args$maxiter <- NULL
  args$threshold <- NULL
  
  markers <- colnames(M)
  outliers.all <- character()
  
  it <- 1
  conv <- FALSE
  while (!conv && it <= maxiter) {
    obj.hima <- do.call("hima", args)
    if (it == 1 && keep.all.first) {
      outliers <- row.names(obj.hima)
    } else {
      outliers <- row.names(obj.hima)[obj.hima$`p-value` < threshold]
    }
    print(outliers)
    outliers.all <- c(outliers.all, outliers)
    is.outlier <- markers %in% outliers
    if (is.null(args$COV)) {
      args$COV <- args$M[, is.outlier, drop = FALSE]
    } else {
      args$COV <- cbind(args$COV, args$M[, is.outlier, drop = FALSE])
    }
    args$M <- args$M[, !is.outlier, drop = FALSE]
    conv <- (length(outliers) == 0)
    it <- it + 1
    markers <- setdiff(markers, outliers)
  }
  
  outliers.all
}