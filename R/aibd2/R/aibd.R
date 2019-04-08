#' Define an Attraction Indian Buffet Distribution (AIBD) for Feature
#' Allocations
#'
#' @param mass The mass (a.k.a., concentration) parameter.
#' @param permutation A permutation, i.e., a vector of integers \code{1, 2, ...,
#'   n} whose length is \code{n} and whose elements are unique.
#' @param similarity A similarity matrix, i.e., a symmetric matrix whose
#'   \code{(i,j)} is large if items \code{i} and \code{j} as similar (i.e., have
#'   a small distance).  An object of class \code{"dist"} is also permissible,
#'   but is interpreted as a similarity (not as a distance).
#' @param permutationMethod A string indicating how the permutation is handled.
#'   \code{'fixed'} indicates that the permutation is not updated from the
#'   original value given by \code{permutation}.  \code{'mcmc'} indicates a
#'   uniform prior on all permutations and updating via Markov chain Monte
#'   Carlo. \code{'enumeration'} indicates a uniform prior on all permutations
#'   which is marginalized away through enumeration.  This last option is only
#'   feasible for very small datasets.
#'
#' @return An object representing an attraction Indian buffet distribution
#'   (AIBD) for feature allocations.
#' @export
#'
#' @examples
#' states <- c("California","Wisconsin","Nebraska","New York")
#' data <- USArrests[states,]
#' dist <- dist(scale(data))
#' similarity <- exp(-1.0*dist)
#' a1 <- aibd(1,seq_along(states),similarity)
#'
#' \dontshow{
#' rscala::scalaDisconnect(aibd2:::s)
#' }
#'
aibd <- function(mass, permutation, similarity, permutationMethod=c("fixed","mcmc","enumeration")[1]) {
  if ( missing(mass) || is.null(mass) || any(is.na(mass)) || any(is.nan(mass)) || !is.numeric(mass) || ( length(mass) != 1 ) ) stop("'mass' is misspecified.")
  mass <- as.double(mass)
  if ( missing(permutation) ) permutation <- NULL
  if ( ! is.null(permutation) ) {
    if ( any(is.na(permutation)) || any(is.nan(permutation)) || !is.numeric(permutation) ) stop("'permutation' is misspecified.")
    permutation <- as.integer(permutation)
    if ( ( min(permutation) < 1 ) || ( max(permutation) > length(permutation) ) || ( length(unique(permutation)) != length(permutation) ) ) stop("'permutation' is misspecified.")
  } else if ( permutationMethod == "fixed" ) stop("'permutation' must be specified if 'permutationMethod' is 'fixed'.")
  if ( ! permutationMethod %in% c("mcmc","fixed","enumeration") ) stop("'permutationMethod' is not recognized.")
  if ( missing(similarity) || is.null(similarity) || any(is.na(similarity)) || any(is.nan(similarity)) || !is.numeric(similarity) ) stop("'similarity' is misspecified.")
  if ( inherits(similarity,"dist") ) similarity <- as.matrix(similarity)
  else if ( is.matrix(similarity) ) {
    if ( !isSymmetric(similarity) ) stop("'similarity' must be symmetric.")
  } else stop("'similarity' should must be a matrix or of class 'dist'.")
  upper <- similarity[upper.tri(similarity)]
  if ( any( upper <= 0.0 ) ) stop("Elements of 'similarity' must be positive.")
  if ( !all(is.finite(upper)) ) stop("Elements of 'similarity' must be finite.")
  labels <- rownames(similarity)
  if ( is.null(labels) ) labels <- 1:nrow(similarity)
  structure(list(mass=mass, nItems=nrow(similarity), permutation=permutation, permutationMethod=permutationMethod, similarity=similarity, labels=labels), class="aibdFADistribution")
}
