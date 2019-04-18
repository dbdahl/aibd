package org.ddahl.aibd.model.lineargaussian

import org.ddahl.aibd.Utils.{harmonicNumber, logFactorial}
import org.apache.commons.math3.random.RandomDataGenerator
import org.apache.commons.math3.util.FastMath.log

class IndianBuffetProcess private (val mass: Double, val nItems: Int) extends FeatureAllocationDistribution with HasMass[IndianBuffetProcess] {

  val logMass = log(mass)

  def updateMass(mass: Double): IndianBuffetProcess = {
    if ( mass <= 0.0 ) throw new IllegalArgumentException("'mass' must be positive.")
    new IndianBuffetProcess(mass, nItems)
  }

  def logProbability(i: Int, fa: FeatureAllocation): Double = logProbability(fa)   // This could be more efficient.

  def logProbability(fa: FeatureAllocation): Double = {
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

  def sample(rdg: RandomDataGenerator): FeatureAllocation = {
    val nNewFeaturesPerItems = Array.tabulate(nItems) { i => rdg.nextPoisson(mass / (i+1)).toInt }
    val nNewFeaturesCumulant = nNewFeaturesPerItems.scan(0)(_+_)
    var fa = FeatureAllocation.empty(nItems)
    var i = 0
    while (i < nItems) {
      var j = 0
      while (j < nNewFeaturesCumulant(i)) {
        if ( rdg.nextUniform(0.0,1.0) <= fa.sizes(j).toDouble / (i+1) ) fa = fa.add(i,j)
        j += 1
      }
      while (j < nNewFeaturesCumulant(i+1) ) {
        fa = fa.add(i)
        j += 1
      }
      i += 1
    }
    fa
  }

}

object IndianBuffetProcess {

  def apply(mass: Double, nItems: Int): IndianBuffetProcess = {
    if ( mass <= 0.0 ) throw new IllegalArgumentException("'mass' must be positive.")
    new IndianBuffetProcess(mass, nItems)
  }

}

