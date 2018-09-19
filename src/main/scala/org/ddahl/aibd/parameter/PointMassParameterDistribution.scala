package org.ddahl.aibd.parameter

import org.apache.commons.math3.random.RandomDataGenerator

case class PointMassParameterDistribution[A](at: A) extends ParameterDistribution[A] {

  def logDensity(x: A) = if (x == at) 0.0 else Double.NegativeInfinity

  def sample(rdg: RandomDataGenerator) = at

  override val discreteSupport = Some(List(at))

}
