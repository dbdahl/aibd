#' Define an Indian Buffet Process (IBP) Distribution for Feature Allocations
#'
#' This function specifies an Indian Buffet Process (IBP), which is a
#' distribution over feature allocations.
#'
#' @param mass The mass (a.k.a., concentration) parameter.
#' @param x A character vector giving the labels of the items to place in a feature allocation, or an integer giving the number of items.
#'
#' @return An object representing an Indian Buffet Process (IBP) feature allocation distribution.
#' @export
#'
#' @examples
#' ibp(1,5)
#' ibp(1,c("CA","WI","NE","NY","UT"))
#'
ibp <- function(mass, x) {
  if ( missing(mass) || is.null(mass) || any(is.na(mass)) || any(is.nan(mass)) || !is.numeric(mass) || ( length(mass) != 1 ) ) stop("'mass' is misspecified.")
  mass <- as.double(mass)
  if ( missing(x) || is.null(x) || any(is.na(x)) || any(is.nan(x)) ) stop("'nItems' is misspecified.")
  if ( length(x) == 1L ) {
    nItems <- as.integer(x)
    labels <- as.character(1:nItems)
  } else {
    nItems <- length(x)
    if ( ! is.character(x) ) stop("A character vector was expected.")
    labels <- x
  }
  structure(list(mass=mass, nItems=nItems, labels=labels), class="ibpFADistribution")
}
