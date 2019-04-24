package org.ddahl.aibd
package model
package regression

import org.ddahl.aibd.parameter.MultivariateNormalParameterDistribution
import org.ddahl.aibd.TimeMonitor
import org.ddahl.commonsmath._
import org.apache.commons.math3.linear.CholeskyDecomposition
import org.apache.commons.math3.random.RandomDataGenerator
import org.apache.commons.math3.util.FastMath._

class RegressionSamplingModel private(val nItems: Int, val nCovariates: Int, val response: Vector[Double], val covariates: Vector[Vector[Double]], private val covariateArrays: Vector[Array[Double]]) {

  def updatedResponse(index: Int, elem: Double) = {
    new RegressionSamplingModel(nItems, nCovariates, response.updated(index, elem), covariates, covariateArrays)
  }

  def mean(i: Int, fa: FeatureAllocation[Vector[Double]]): Double = {
    RegressionSamplingModel.mean(i, fa, covariates(i))
  }

  def mean(fa: FeatureAllocation[Vector[Double]]): Vector[Double] = {
    Vector.tabulate(nItems)(i => RegressionSamplingModel.mean(i, fa, covariates(i)))
  }

  def coefficients(i: Int, fa: FeatureAllocation[Vector[Double]]): Vector[Double] = {
    RegressionSamplingModel.coefficients(i, fa, nCovariates)
  }

  def coefficients(fa: FeatureAllocation[Vector[Double]]): Vector[Vector[Double]] = {
    Vector.tabulate(nItems)(i => RegressionSamplingModel.coefficients(i, fa, nCovariates))
  }

  private final val normalizingConstant = -0.5 * log(2 * PI)

  def logLikelihood(i: Int, precision: Double, mean: Double): Double = {
    val r = response(i) - mean
    normalizingConstant - 0.5 * r * r * precision
  }

  val timer = TimeMonitor()

  def logLikelihood(i: Int, precision: Double, fa: FeatureAllocation[Vector[Double]]): Double = {
    timer {
      logLikelihood(i, precision, RegressionSamplingModel.mean(i, fa, covariates(i)))
    }
  }

  def logLikelihood(precision: Double, fa: FeatureAllocation[Vector[Double]]): Double = {
    var sum = 0.0
    var i = 0
    while (i < nItems) {
      sum += logLikelihood(i, precision, fa)
      i += 1
    }
    sum
  }

  def logLikelihood(precision: Double, mean: Vector[Double]): Double = {
    var sum = 0.0
    var i = 0
    while (i < nItems) {
      sum += logLikelihood(i, precision, mean(i))
      i += 1
    }
    sum
  }

  def mkLogLikelihood(precision: Double): (Int, FeatureAllocation[Vector[Double]]) => Double = (i: Int, fa: FeatureAllocation[Vector[Double]]) => {
    logLikelihood(i, precision, fa)
  }

  def updatePrecision(fa: FeatureAllocation[Vector[Double]], rdg: RandomDataGenerator, shape: Double, rate: Double): Double = {
    val m = mean(fa)
    var sum = 0.0
    var i = 0
    while (i < nItems) {
      val r = response(i) - m(i)
      sum += r * r
      i += 1
    }
    val newShape = shape + nItems / 2.0
    val newRate = rate + 0.5 * sum
    rdg.nextGamma(newShape, 1.0 / newRate)
  }

  def updateCoefficients(fa: FeatureAllocation[Vector[Double]], rdg: RandomDataGenerator, precision: Double, prior: MultivariateNormalParameterDistribution): FeatureAllocation[Vector[Double]] = {
    var faCurrent = fa
    val priorPrecisionTimesMean = prior.precisionMatrix.multiply(prior.meanMatrix)
    for (fCurrent <- fa) {
      faCurrent = faCurrent.remove(fCurrent)
      val which = fCurrent.set.toArray
      val y = which.map(response)
      val X = MatrixFactory(which.map(covariateArrays), false)
      val Xt = X.transpose
      val precisionMatrix = (Xt.multiply(X)).scalarMultiply(precision).add(prior.precisionMatrix)
      val r = Array.ofDim[Double](y.size)
      var i = 0
      while (i < y.size) {
        r(i) = y(i) - RegressionSamplingModel.mean(which(i), faCurrent, covariateArrays(which(i)))
        i += 1
      }
      val scaledXtX = Xt.operate(r).map(_ * precision)
      val numerator = priorPrecisionTimesMean.add(MatrixFactory(Array.tabulate(nCovariates, 1)((i, j) => scaledXtX(i)), false))
      val meanMatrix = new CholeskyDecomposition(precisionMatrix).getSolver.solve(numerator)
      val result = MultivariateNormalParameterDistribution.usingPrecision(meanMatrix, precisionMatrix).sample(rdg)
      faCurrent = faCurrent.add(fCurrent.replace(result))
    }
    faCurrent
  }

}

object RegressionSamplingModel {

  def apply(covariates: Array[Array[Double]], precision: Double, fa: FeatureAllocation[Vector[Double]], rdg: RandomDataGenerator): RegressionSamplingModel = {
    val stdDevOfError = 1.0 / sqrt(precision)
    val response = Vector.tabulate(covariates.length)(i => rdg.nextGaussian(mean(i, fa, covariates(i)), stdDevOfError))
    apply(response, covariates)
  }

  def apply(covariates: Vector[Vector[Double]], precision: Double, fa: FeatureAllocation[Vector[Double]], rdg: RandomDataGenerator): RegressionSamplingModel = {
    apply(covariates.map(_.toArray).toArray, precision, fa, rdg)
  }

  def apply(response: IndexedSeq[Double], covariates: Array[Array[Double]]): RegressionSamplingModel = {
    apply(response, covariates.map(_.toVector).toVector)
  }

  def apply(response: IndexedSeq[Double], covariates: Vector[Vector[Double]]): RegressionSamplingModel = {
    val nItems: Int = covariates.length
    require(nItems > 0)
    val nCovariates: Int = covariates(0).length
    require(covariates.map(_.length).toSet.size == 1)
    require(response.size == nItems)
    val covariateArrays = covariates.map(_.toArray)
    new RegressionSamplingModel(nItems, nCovariates, response.toVector, covariates, covariateArrays)
  }

  private[model] def mean(i: Int, fa: FeatureAllocation[Vector[Double]], covars: IndexedSeq[Double]): Double = {
    val nCovariates = covars.length
    var sum = 0.0
    var p = 0
    while (p < nCovariates) {
      var sum2 = 0.0
      for ( f <- fa ) {
        if (f.contains(i)) sum2 += f.parameter(p)
      }
      sum += covars(p) * sum2
      p += 1
    }
    sum
  }

  private[model] def coefficients(i: Int, fa: FeatureAllocation[Vector[Double]], nCovariates: Int): Vector[Double] = {
    val coef = Array.ofDim[Double](nCovariates)
    var p = 0
    while (p < nCovariates) {
      var sum2 = 0.0
      for ( f <- fa ) {
        if (f.contains(i)) sum2 += f.parameter(p)
      }
      coef(p) = sum2
      p += 1
    }
    coef.toVector
  }

}
