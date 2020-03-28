#' @import utils
#'
.onLoad <- function(libname, pkgname) {
  assign("s", NULL, envir=parent.env(environment()))
  globalVariables("s")
}

#' @import rscala commonsMath
#'
scalaEnsure <- function() {
  if ( ! is.null(s) ) return()
  commonsMath:::.packageName   # So CRAN checks recognize that its being used.
  s <- scala("commonsMath")
  scalaLazy(function(s) s + '
    import org.apache.commons.math3.random.{ RandomDataGenerator => RDG }
    import org.ddahl.matrix._
    import org.ddahl.aibd.{FeatureAllocation => FA}
    import org.ddahl.aibd.distribution._
    import org.ddahl.aibd.model.lineargaussian.{LinearGaussianLatentFeatureModel => LGLFM, _}

    def rdg() = {
      val ints = R.evalI1("runif(2,-.Machine$integer.max,.Machine$integer.max)")
      val seed = ((ints(0).asInstanceOf[Long]) << 32) | (ints(1) & 0xffffffffL)
      val r = new RDG()
      r.reSeed(seed)
      r
    }
  ')
  env <- parent.env(environment())
  # unlockBinding("s", env)
  eval(parse(text=paste0('unlockBinding("s",env)')))
  assign("s", s, envir=env)
  lockBinding("s", env)
}

.onUnload <- function(libpath) {
  if ( ! is.null(s) ) close(s)
}
