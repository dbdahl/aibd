package org.ddahl.aibd.model.lineargaussian

import org.ddahl.aibd.Utils._
import org.apache.commons.math3.util.FastMath.log

import scala.collection.mutable.BitSet

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
          if ( fa.array(j)(ii) ) {
            state.mutateAdd(ii,j)
            newFeatureCount += 1
          }
        } else {
          val p = index.toDouble / (index + 1) * state.featuresAsList(j).foldLeft(0.0) { (s, iPrime) => s + similarity(ii)(iPrime) } / divisor
          if ( fa.array(j)(ii) ) {
            state.mutateAdd(ii,j)
            sum += log(p)
          } else sum += log(1-p)
        }
        j += 1
      }
      index += 1
      sum += newFeatureCount * ( logMass - logOnInt(index) ) - mass/index
    }
    sum -= computeRegardingTies(fa)
    sum
  }

  def logProbabilityIBP(fa: FeatureAllocation, mass: Double, logMass: Double): Double = {
    val const1 = -mass * harmonicNumber(fa.nItems)
    if ( fa.nFeatures == 0 ) return const1
    val const2 = logMass - logFactorial(fa.nItems)
    var sum = const1 + fa.nFeatures * const2
    sum -= computeRegardingTies(fa)
    sum += fa.sizes.foldLeft(0.0) { (s, mk) =>
      s + logFactorial(fa.nItems - mk) + logFactorial(mk - 1)
    }
    sum
  }

  def computeRegardingTiesSlow(fa: FeatureAllocation): Double = {
    fa.array.groupBy(identity).map(_._2.length).foldLeft(0.0)((s, x) => s + logFactorial(x))
  }

  def computeRegardingTies(fa: FeatureAllocation): Double = {
    val aa = fa.array.zip(fa.sizes).sortWith(lessThan)
    var sum = 0.0
    var run = 1
    var j = 1
    while ( j < aa.length ) {
      if ( compare(aa(j-1),aa(j)) == 0 ) run += 1
      else {
        sum += logFactorial(run)
        run = 1
      }
      j += 1
    }
    sum += logFactorial(run)
    sum
  }

  def compare(x: (BitSet,Int), y: (BitSet, Int)): Int = {
    if ( x._2 < y._2 ) return -1
    else if ( x._2 > y._2 ) return 1
    else {
      val xi = x._1.iterator
      val yi = y._1.iterator
      while (xi.hasNext && yi.hasNext) {
        val xn = xi.next()
        val yn = yi.next()
        if (xn < yn) return -1
        if (xn > yn) return 1
      }
      if (xi.hasNext) return 1
      if (yi.hasNext) return -1
      0
    }
  }

  def lessThan(x: (BitSet,Int), y: (BitSet, Int)): Boolean = compare(x,y) <= 0

}

