package org.ddahl.aibd.parameter

import org.apache.commons.math3.random.RandomDataGenerator
import org.apache.commons.math3.distribution.NormalDistribution
import org.apache.commons.math3.util.FastMath.sqrt

case class GaussianParameterDistribution(mean: Double, precision: Double) extends ParameterDistribution[Double] {

  require(!mean.isInfinity && !mean.isNaN)
  require(!precision.isInfinity && !precision.isNaN && precision > 0.0)

  val standardDeviation = 1 / sqrt(precision)

  private val nd = new NormalDistribution(null, mean, standardDeviation)

  def logDensity(x: Double) = nd.logDensity(x)

  def sample(rdg: RandomDataGenerator) = rdg.nextGaussian(mean, standardDeviation)

}
