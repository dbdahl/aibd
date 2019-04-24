#' @import rscala
#'
.onLoad <- function(libname, pkgname) {
  s <- sdols:::s
  scalaJARs(pkgname,s)
  scalaLazy(function(s) {
    s + '
      import org.ddahl.matrix._
      import org.ddahl.aibd.{FeatureAllocation => FA}
      import org.ddahl.aibd.distribution._
      import org.ddahl.aibd.model.lineargaussian.{LinearGaussianLatentFeatureModel => LGLFM, _}
      import org.ddahl.sdols.featureallocation.FeatureAllocation
    '
  },s)
  assign("s",s,envir=parent.env(environment()))
}
