package org.ddahl.aibd
package parameter

import org.apache.commons.math3.random.RandomDataGenerator
import org.apache.commons.math3.util.FastMath.{sqrt,log,PI}

class IndependentNormalsParameterDistribution private (val mean: Vector[Double], varianceOption: Option[Vector[Double]], precisionOption: Option[Vector[Double]]) extends ParameterDistribution[Vector[Double]] {

  private val k = mean.length
  require(k > 0)

  require(varianceOption.isDefined || precisionOption.isDefined)

  lazy val variance: Vector[Double] = if (varianceOption.isDefined) varianceOption.get
  else precisionOption.get.map(1/_)

  lazy val precision: Vector[Double] = if (precisionOption.isDefined) precisionOption.get
  else varianceOption.get.map(1/_)

  private val sd = if ( varianceOption.isDefined ) varianceOption.get.map(sqrt)
  else precisionOption.get.map(y => sqrt(1/y))

  require(sd.length == k)

  private lazy val logNormalizingConstant = -0.5 * k * log(2 * PI) - sd.foldLeft(0.0)((sum,s) => sum + log(s))

  def logDensity(x: Vector[Double]): Double = {
    require(x.size == k)
    logNormalizingConstant - 0.5 * Range(0,k).foldLeft(0.0) { (sum,i) => val y = (x(i)-mean(i))/sd(i); sum + y*y }
  }

  def sample(rdg: RandomDataGenerator): Vector[Double] = {
    Vector.tabulate(k) { i =>
      rdg.nextGaussian(mean(i),sd(i))
    }
  }

}

object IndependentNormalsParameterDistribution {

  def usingVariance(mean: Vector[Double], variance: Vector[Double]) = {
    new IndependentNormalsParameterDistribution(mean, Some(variance), None)
  }

  def usingVariance(mean: Array[Double], variance: Array[Double]) = {
    new IndependentNormalsParameterDistribution(mean.toVector, Some(variance.toVector), None)
  }

  def usingPrecision(mean: Vector[Double], precision: Vector[Double]) = {
    new IndependentNormalsParameterDistribution(mean, None, Some(precision))
  }

  def usingPrecision(mean: Array[Double], precision: Array[Double]) = {
    new IndependentNormalsParameterDistribution(mean.toVector, None, Some(precision.toVector))
  }

}

