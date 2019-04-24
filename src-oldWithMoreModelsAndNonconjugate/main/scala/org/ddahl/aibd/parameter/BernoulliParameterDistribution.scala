package org.ddahl.aibd.parameter

import org.apache.commons.math3.random.RandomDataGenerator

case class BernoulliParameterDistribution(probability: Double) extends ParameterDistribution[Int] {

  require(0.0 < probability && probability < 1.0)

  import org.apache.commons.math3.util.FastMath.log

  private val logProbability = log(probability)
  private val log1MinusProbability = log(1.0 - probability)

  def logDensity(x: Int) = if (x == 1) logProbability else if (x == 0) log1MinusProbability else Double.NegativeInfinity

  def sample(rdg: RandomDataGenerator) = if (rdg.nextUniform(0.0, 1.0) <= probability) 1 else 0

  override val discreteSupport = Some(List(0, 1))

}
