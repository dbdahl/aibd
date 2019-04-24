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
aibd <- function(mass, permutation, similarity) {
  if ( missing(mass) || is.null(mass) || any(is.na(mass)) || any(is.nan(mass)) || !is.numeric(mass) || ( length(mass) != 1 ) ) stop("'mass' is misspecified.")
  mass <- as.double(mass)
  if ( missing(permutation) ) permutation <- NULL
  if ( ! is.null(permutation) ) {
    if ( any(is.na(permutation)) || any(is.nan(permutation)) || !is.numeric(permutation) ) stop("'permutation' is misspecified.")
    permutation <- as.integer(permutation)
    if ( ( min(permutation) < 1 ) || ( max(permutation) > length(permutation) ) || ( length(unique(permutation)) != length(permutation) ) ) stop("'permutation' is misspecified.")
  }
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
  structure(list(mass=mass, nItems=nrow(similarity), permutation=permutation, similarity=similarity, labels=labels), class="aibdFADistribution")
}
