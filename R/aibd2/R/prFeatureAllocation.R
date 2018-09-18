#' Evaluation a Probabilty Mass Function of a Feature Allocation Distribution
#'
#' This is great!!
#'
#' @param featureAllocation test
#' @param distribution A feature allocation distribution
#' @param log Should results be given on the log scale?
#'
#' @return The probability (or log of the probability) of the given feature allocation.
#' @export
#'
#' @examples
#' fa <- 1
#' d <- ibp(2,nrow(fa))
#' prFeatureAllocation(fa,d)
prFeatureAllocation <- function(featureAllocation, distribution, log=FALSE) {
  if ( inherits(distribution,"ibpFADistribution") ) stop("Unsupported distribution.")
  # Code goes here!
  if ( log ) 0.0 else 1.0
}
