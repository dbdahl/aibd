package org.ddahl.aibd.model.lineargaussian

import org.ddahl.aibd.Utils._
import org.apache.commons.math3.util.FastMath.log

object FeatureAllocationDistributions {

  def logProbabilityIBP(fa: FeatureAllocation, mass: Double): Double = {
    val const1 = -mass * harmonicNumber(fa.nItems)
    if ( fa.nFeatures == 0 ) return const1
    val const2 = log(mass) - logFactorial(fa.nItems)
    var sum = const1 + fa.nFeatures * const2
    sum = -fa.array.groupBy(identity).map(_._2.length).foldLeft(0.0)((s, x) => s + logFactorial(x))
    sum += fa.sizes.foldLeft(0.0) { (s, mk) =>
      s + logFactorial(fa.nItems - mk) + logFactorial(mk - 1)
    }
    sum
  }

}
