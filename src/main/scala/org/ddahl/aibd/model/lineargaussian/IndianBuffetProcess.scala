package org.ddahl.aibd.model.lineargaussian

import org.ddahl.aibd.{IndianBuffetProcess => IndianBuffetProcessAlternative}
import org.ddahl.aibd.Utils.{harmonicNumber, logFactorial}
import org.apache.commons.math3.random.RandomDataGenerator
import org.apache.commons.math3.util.FastMath.log

class IndianBuffetProcess private (val mass: Double, val nItems: Int) extends FeatureAllocationDistribution {

  val logMass = log(mass)

  def logProbability(i: Int, fa: FeatureAllocation): Double = {
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
    val ibp = IndianBuffetProcessAlternative(mass, nItems)
    FeatureAllocation.convertFromAlternativeImplementation(ibp.sample(rdg))
  }

}

object IndianBuffetProcess {

  def apply(mass: Double, nItems: Int): IndianBuffetProcess = {
    if ( mass <= 0.0 ) throw new IllegalArgumentException("'mass' must be positive.")
    new IndianBuffetProcess(mass, nItems)
  }

}

