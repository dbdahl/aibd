package org.ddahl.aibd
package parameter

import org.apache.commons.math3.linear.{ArrayRealVector, CholeskyDecomposition, LUDecomposition, RealMatrix}
import org.apache.commons.math3.random.RandomDataGenerator
import org.apache.commons.math3.util.FastMath.{PI, log}

// Be careful:  Objects of this class are not immutable since RealMatrix is mutable.
class MultivariateNormalParameterDistribution private (val meanMatrix: RealMatrix, covarianceOption: Option[RealMatrix], precisionOption: Option[RealMatrix]) extends ParameterDistribution[Vector[Double]] {

  private val k = meanMatrix.getRowDimension
  require(k > 0)
  require(meanMatrix.getColumnDimension == 1)
  // require((0 until k).forall(i => !mean.getEntry(i, 0).isInfinity && !mean.getEntry(i, 0).isNaN))
  require(covarianceOption.isDefined || precisionOption.isDefined)
  if (covarianceOption.isDefined) require(covarianceMatrix.getRowDimension == k)
  if (precisionOption.isDefined) require(precisionMatrix.getRowDimension == k)

  lazy val covarianceCholesky = new CholeskyDecomposition(covarianceMatrix)
  lazy val precisionCholesky = new CholeskyDecomposition(precisionMatrix)

  lazy val covarianceSolver = covarianceCholesky.getSolver
  lazy val precisionSolver = new LUDecomposition(precisionCholesky.getLT).getSolver

  lazy val covarianceMatrix: RealMatrix = if (covarianceOption.isDefined) covarianceOption.get
  else precisionCholesky.getSolver.getInverse

  lazy val precisionMatrix: RealMatrix = if (precisionOption.isDefined) precisionOption.get
  else covarianceCholesky.getSolver.getInverse

  private lazy val logNormalizingConstant = if (precisionOption.isDefined) {
    -0.5 * k * log(2 * PI) + 0.5 * log(precisionCholesky.getDeterminant)
  } else {
    -0.5 * k * log(2 * PI) - 0.5 * log(covarianceCholesky.getDeterminant)
  }

  def logDensity(x: Vector[Double]): Double = {
    require(x.size == k)
    val centered = x.toArray
    var i = 0
    while (i < k) {
      centered(i) -= meanMatrix.getEntry(i, 0)
      i += 1
    }
    val dotProduct = if (precisionOption.isDefined) {
      val z = precisionCholesky.getLT.operate(centered)
      z.foldLeft(0.0)((sum, x) => sum + x * x)
    } else {
      val y = covarianceSolver.solve(new ArrayRealVector(centered, false))
      y.dotProduct(new ArrayRealVector(centered, false))
    }
    logNormalizingConstant - 0.5 * dotProduct
  }

  def sample(rdg: RandomDataGenerator, mean: Option[RealMatrix]): Vector[Double] = {
    val x = Array.fill(k)(rdg.nextGaussian(0, 1))
    if (covarianceOption.isDefined) {
      val result = covarianceCholesky.getLT.preMultiply(x)
      if ( mean.isDefined ) {
        val m = mean.get
        require(m.getRowDimension == k)
        require(m.getColumnDimension == 1)
        Vector.tabulate(k)(i => result(i) + m.getEntry(i, 0))
      } else {
        result.toVector
      }
    } else {
      val result = precisionSolver.solve(new ArrayRealVector(x, false))
      if ( mean.isDefined ) {
        val m = mean.get
        require(m.getRowDimension == k)
        require(m.getColumnDimension == 1)
        Vector.tabulate(k)(i => result.getEntry(i) + m.getEntry(i, 0))
      } else {
        result.toArray.toVector
      }
    }
  }

  def sample(rdg: RandomDataGenerator): Vector[Double] = {
    sample(rdg,Some(meanMatrix))
  }

}

object MultivariateNormalParameterDistribution {

  def usingCovariance(mean: RealMatrix, covariance: RealMatrix) = {
    new MultivariateNormalParameterDistribution(mean, Some(covariance), None)
  }

  def usingPrecision(mean: RealMatrix, precision: RealMatrix) = {
    new MultivariateNormalParameterDistribution(mean, None, Some(precision))
  }

}
