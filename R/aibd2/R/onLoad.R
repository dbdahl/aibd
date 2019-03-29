#' @import rscala
#'
.onLoad <- function(libname, pkgname) {
  s <- sdols:::s
  scalaJARs(pkgname,s)
  scalaLazy(function(s) {
    s + '
      import org.ddahl.aibd._
      import org.ddahl.aibd.model.lineargaussian.LinearGaussianSamplingModel
      import org.ddahl.aibd.model.lineargaussian.{LinearGaussianLatentFeatureModel => LGLFM}
      import org.ddahl.aibd.model.lineargaussian.{FeatureAllocation => FA}
      import org.ddahl.aibd.model.lineargaussian.PosteriorSimulation
      import org.ddahl.matrix._
      import org.ddahl.sdols.featureallocation.FeatureAllocation
    '
  },s)
  assign("s",s,envir=parent.env(environment()))
}
