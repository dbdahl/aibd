#' @import rscala
#'
.onLoad <- function(libname, pkgname) {
  pkgEnv <- parent.env(environment())
  flagEnv <- new.env(parent=pkgEnv)
  assign("flagEnv",flagEnv,envir=pkgEnv)
  assign("isLoaded",FALSE,envir=flagEnv)
  delayedAssign("s",{
    assign("isLoaded",TRUE,envir=flagEnv)
    s <- scala("commonsMath")
    s + '
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
    '
    s
  },assign.env=pkgEnv)
}

.onUnload <- function(libpath) {
  if ( get("isLoaded",envir=flagEnv) ) close(s)
}

