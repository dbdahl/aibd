package org.ddahl.aibd.model.lineargaussian

import org.ddahl.aibd.Utils._
import org.apache.commons.math3.util.FastMath.log

object FeatureAllocationDistributions {

  def logProbabilityAIBD(i: Int, fa: FeatureAllocation, mass: Double, logMass: Double, permutation: Array[Int], invertedPermutation: Array[Int], similarity: Array[Array[Double]]): Double = {
    var index = invertedPermutation(i)
    val state = fa.remove(permutation.drop(index), true)
    var sum = 0.0
    while ( index < fa.nItems ) {
      val ii = permutation(index)
      val divisor = (0 until index).foldLeft(0.0) { (s,index2) => s + similarity(ii)(permutation(index2)) }
      var newFeatureCount = 0
      var j = 0
      while ( j < fa.nFeatures ) {
        if ( state.sizes(j) == 0 ) {
          if ( fa.features(j)(ii) ) {
            state.mutateAdd(ii,j)
            newFeatureCount += 1
          }
        } else {
          val p = index.toDouble / (index + 1) * state.featuresAsList(j).foldLeft(0.0) { (s, iPrime) => s + similarity(ii)(iPrime) } / divisor
          if ( fa.features(j)(ii) ) {
            state.mutateAdd(ii,j)
            sum += log(p)
          } else sum += log(1-p)
        }
        j += 1
      }
      index += 1
      sum += newFeatureCount * ( logMass - logOnInt(index) ) - mass/index
    }
    sum -= fa.computeRegardingTies
    sum
  }

  def logProbabilityIBP(fa: FeatureAllocation, mass: Double, logMass: Double): Double = {
    val const1 = -mass * harmonicNumber(fa.nItems)
    if ( fa.nFeatures == 0 ) return const1
    val const2 = logMass - logFactorial(fa.nItems)
    var sum = const1 + fa.nFeatures * const2
    sum -= fa.computeRegardingTies
    sum += fa.sizes.foldLeft(0.0) { (s, mk) =>
      s + logFactorial(fa.nItems - mk) + logFactorial(mk - 1)
    }
    sum
  }

}

