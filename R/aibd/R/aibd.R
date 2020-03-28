#' Define an Attraction Indian Buffet Distribution (AIBD) for Feature Allocations
#'
#' This function specifies an Attraction Indian Buffet Distribution (AIBD), which is a
#' distribution over feature allocations.
#'
#' @param mass The mass (a.k.a., concentration) parameter of the AIBD.
#' @param permutation A permutation, i.e., a vector of integers \code{1, 2, ...,
#'   n} whose length is \code{n} and whose elements are unique.  Using the Indian
#'   buffet analogy, the permutation represents the order the customers enter the buffet.
#' @param temperature A nonnegative scalar which determines how influential the
#'   distance matrix is in the feature allocation distribution.  The AIBD reduces to
#'   the IBP when the temperature is zero and diverges from the IBP as the temperature increases.
#' @param distance A distance matrix, i.e., a symmetric matrix whose
#'   \code{(i,j)} entry is small if items \code{i} and \code{j} are similar.  An object
#'   of class \code{"dist"} is also permissible.
#' @param decayFunction One of the following strings: \code{"exponential"}
#'   (making \code{similarity = exp(-temperature*distance)}),
#'   \code{"reciprocal"} (making \code{similarity = 1/distance^temperature}), or
#'   \code{"identity"} (in which case \code{distance} is interpreted as a
#'   similarity instead of a distance).
#'
#' @return An object representing an Attraction Indian Buffet Distribution
#'   (AIBD) for feature allocations.
#' @export
#'
#' @examples
#' states <- c("California","Wisconsin","Nebraska","New York")
#' data <- USArrests[states,]
#' dist <- dist(scale(data))
#' aibd(1, seq_along(states), 1.0, dist)
#'
aibd <- function(mass, permutation, temperature, distance, decayFunction=c("exponential","reciprocal","identity")[1]) {
  if ( missing(mass) || is.null(mass) || any(is.na(mass)) || any(is.nan(mass)) || !is.numeric(mass) || ( length(mass) != 1 ) ) stop("'mass' is misspecified.")
  mass <- as.double(mass)
  if ( missing(permutation) ) permutation <- NULL
  if ( ! is.null(permutation) ) {
    if ( any(is.na(permutation)) || any(is.nan(permutation)) || !is.numeric(permutation) ) stop("'permutation' is misspecified.")
    permutation <- as.integer(permutation)
    if ( ( min(permutation) < 1 ) || ( max(permutation) > length(permutation) ) || ( length(unique(permutation)) != length(permutation) ) ) stop("'permutation' is misspecified.")
  }
  if ( missing(temperature) || is.null(temperature) || any(is.na(temperature)) || any(is.nan(temperature)) || !is.numeric(temperature) || ( length(temperature) != 1 ) ) stop("'temperature' is misspecified.")
  if ( missing(distance) || is.null(distance) || any(is.na(distance)) || any(is.nan(distance)) || !is.numeric(distance) ) stop("'distance' is misspecified.")
  if ( inherits(distance,"dist") ) distance <- as.matrix(distance) else if ( is.matrix(distance) ) {
    if ( !isSymmetric(distance) ) stop("'distance' must be symmetric.")
  } else stop("'distance' should must be a matrix or of class 'dist'.")
  upper <- distance[upper.tri(distance)]
  if ( any( upper <= 0.0 ) ) stop("Elements of 'distance' must be positive.")
  if ( !all(is.finite(upper)) ) stop("Elements of 'distance' must be finite.")
  labels <- rownames(distance)
  if ( is.null(labels) ) labels <- 1:nrow(distance)
  if ( ! ( is.character(decayFunction) && ( length(decayFunction) == 1 ) && ( decayFunction %in% c("exponential","reciprocal","identity") ) ) ) stop("'decayFunction' is misspecified.")
  structure(list(mass=mass, nItems=nrow(distance), permutation=permutation, temperature=temperature, distance=distance, decayFunction=decayFunction, labels=labels), class="aibdFADistribution")
}
