package org.ddahl.aibd.parameter

import org.apache.commons.math3.random.RandomDataGenerator

trait ParameterDistribution[A] {

  type tipe = A

  def logDensity(x: A): Double

  def sample(rdg: RandomDataGenerator): A

  val discreteSupport: Option[List[A]] = None

}
