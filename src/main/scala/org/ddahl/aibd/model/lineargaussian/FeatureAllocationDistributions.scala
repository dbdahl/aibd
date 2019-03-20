package org.ddahl.aibd.model.lineargaussian

import org.ddahl.matrix._
import org.ddahl.aibd.Utils._
import org.apache.commons.math3.util.FastMath.log

object FeatureAllocationDistributions {

  def logProbabilityIBP(Z: Matrix, N: Int, mass: Double): Double = {
    val const1 = -mass * harmonicNumber(N)
    if ( Z == null ) return const1
    val fp = FeatureAllocationUtilities.toFingerprint(Z)
    val const2 = log(mass) - logFactorial(N)
    var sum = const1 + fp.length * const2
    val a = fp.groupBy(identity).map(_._2.length)
    sum -= a.foldLeft(0.0)((s, x) => s + logFactorial(x))
    sum += fp.aggregate(0.0)((s, f) => {
      val mk = f.size
      s + logFactorial(N - mk) + logFactorial(mk - 1)
    }, _ + _)
    sum
  }

}
