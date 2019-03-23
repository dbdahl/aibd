package org.ddahl.aibd.model.lineargaussian

import org.ddahl.aibd.Utils._
import org.apache.commons.math3.util.FastMath.log

import scala.collection.mutable.BitSet

object FeatureAllocationDistributions {

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

  def logProbabilityIBP(fa: FeatureAllocation, mass: Double): Double = {
    val const1 = -mass * harmonicNumber(fa.nItems)
    if ( fa.nFeatures == 0 ) return const1
    val const2 = log(mass) - logFactorial(fa.nItems)
    var sum = const1 + fa.nFeatures * const2
    //sum -= fa.array.groupBy(identity).map(_._2.length).foldLeft(0.0)((s, x) => s + logFactorial(x))
    val aa = fa.array.zip(fa.sizes).sortWith(lessThan)
    var sum2 = 0.0
    var run = 1
    var j = 1
    while ( j < aa.length ) {
      if ( compare(aa(j-1),aa(j)) == 0 ) run += 1
      else {
        sum2 += logFactorial(run)
        run = 1
      }
      j += 1
    }
    sum2 += logFactorial(run)
    sum -= sum2
    sum += fa.sizes.foldLeft(0.0) { (s, mk) =>
      s + logFactorial(fa.nItems - mk) + logFactorial(mk - 1)
    }
    sum
  }

}
