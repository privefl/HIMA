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
#' @param keep.all.first Do not use threshold on the first round?
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
  
  markers <- colnames(M)
  outliers.all <- character()
  
  it <- 1
  conv <- FALSE
  while (!conv && it <= maxiter) {
    obj.hima <- hima(
      X = X, Y = Y, M = M, COV = COV, 
      family = family, 
      penalty = penalty, 
      topN = topN, 
      parallel = parallel, 
      ncore = ncore, 
      verbose = verbose, 
      ...
    )
    if (it == 1 && keep.all.first) {
      outliers <- row.names(obj.hima)
    } else {
      outliers <- row.names(obj.hima)[obj.hima$`p-value` < threshold]
    }
    print(outliers)
    outliers.all <- c(outliers.all, outliers)
    is.outlier <- markers %in% outliers
    if (is.null(COV)) {
      COV <- M[, is.outlier, drop = FALSE]
    } else {
      COV <- cbind(COV, M[, is.outlier, drop = FALSE])
    }
    M <- M[, !is.outlier, drop = FALSE]
    conv <- (length(outliers) == 0)
    it <- it + 1
    markers <- setdiff(markers, outliers)
  }
  
  outliers.all
}
