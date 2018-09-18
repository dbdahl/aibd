#' @import rscala
#'
.onLoad <- function(libname, pkgname) {
  s <- sdols:::s
  scalaLazy(function(s) {
    s + 'import org.ddahl.sdols.featureallocation.FeatureAllocation'
  },s)
  assign("s",s,envir=parent.env(environment()))
}
