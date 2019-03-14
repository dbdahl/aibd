#' @import rscala
#'
.onLoad <- function(libname, pkgname) {
  s <- sdols:::s
  scalaJARs(c("breezeJARs",pkgname),s)
  scalaLazy(function(s) {
    s + '
      import org.ddahl.aibd._
      import org.ddahl.aibd.model.lineargaussian.LinearGaussianSamplingModel
      import org.ddahl.aibd.model.lineargaussian.{LinearGaussianLatentFeatureModel => LGLFM}
      import org.ddahl.aibd.model.lineargaussian.{FeatureAllocationUtilities => FAU}
      import org.ddahl.sdols.featureallocation.FeatureAllocation
    '
  },s)
  assign("s",s,envir=parent.env(environment()))
}
