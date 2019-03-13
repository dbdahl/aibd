package org.ddahl.lglfm

import breeze.linalg._
import org.apache.commons.math3.util.FastMath.log

class LinearGaussianLatentFeatureModel private (val X: DenseMatrix[Double], val precisionX: Double, val precisionW: Double) {

  val N = X.rows
  val D = X.cols

  private val Xt = X.t
  private val I = DenseMatrix.eye[Double](N)
  private val Dhalf = D / 2.0
  private val const = -N * Dhalf * log(2 * math.Pi)
  private val XtXtrace = trace(Xt * X)
  private val DhalfTimeslogPrecisionX = Dhalf * log(precisionX)
  private val halfPrecisionX = precisionX / 2
  private val logPrecisionW = log(precisionW)
  private val ratioOfPrecisionsTimesI = I *:* precisionW / precisionX

  def logLikelihood(lc: LikelihoodComponents): Double = {
    if (lc.K == 0) {
      const + (N) * DhalfTimeslogPrecisionX - halfPrecisionX * XtXtrace
    } else {
      const + (N - lc.K) * DhalfTimeslogPrecisionX + lc.K * Dhalf * logPrecisionW + Dhalf * log(det(lc.M)) - halfPrecisionX * trace(Xt * (I - lc.Z * lc.M * lc.Zt) * X)
    }
  }

  def logLikelihood(Z: DenseMatrix[Double]): Double = logLikelihood(computeLikelihoodComponents(Z))

  class LikelihoodComponents private[LinearGaussianLatentFeatureModel] (val Z: DenseMatrix[Double], val Zt: DenseMatrix[Double], val M: DenseMatrix[Double]) {
    val K = Z.cols
  }

  def computeLikelihoodComponents(Z: DenseMatrix[Double]) = {
    if (Z.rows != N) throw new IllegalArgumentException("Feature allocation has " + Z.rows + " items, but " + N + " were expected.")
    val Zt = Z.t
    val M = inv(Zt * Z + ratioOfPrecisionsTimesI)
    new LikelihoodComponents(Z, Zt, M)
  }

  private def update(M: DenseMatrix[Double], z: Transpose[DenseVector[Double]], add: Boolean) = {
    val zt = z.t
    M - M * zt * z * M /:/ (z * M * zt + ( if (add) 1 else -1))
  }

  def dropFeaturesFor(i: Int, lc: LikelihoodComponents): LikelihoodComponents = {
    val Z = lc.Z.copy
    val zOld = Z(i,::)
    Z(i,::) := DenseVector.zeros[Double](lc.K).t
    new LikelihoodComponents(Z, Z.t, update(lc.M, zOld, false))
  }

  def addFeaturesFor(i: Int, lc: LikelihoodComponents, z: Transpose[DenseVector[Double]]): LikelihoodComponents = {
    val Z = lc.Z.copy
    // assert(!any(Z(i,::)))
    Z(i,::) := z
    new LikelihoodComponents(Z, Z.t, update(lc.M, z, true))
  }

  def dropAndAddFeaturesFor(i: Int, lc: LikelihoodComponents, z: Transpose[DenseVector[Double]]): LikelihoodComponents = {
    val newLC = dropFeaturesFor(i, lc)
    newLC.Z(i,::) := z
    new LikelihoodComponents(newLC.Z, newLC.Z.t, update(newLC.M, z, true))
  }

}


object LinearGaussianLatentFeatureModel {

  def usingPrecisions(X: DenseMatrix[Double], precisionX: Double, precisionW: Double): LinearGaussianLatentFeatureModel = {
    new LinearGaussianLatentFeatureModel(X.copy, precisionX, precisionW)
  }

  def usingStandardDeviations(X: DenseMatrix[Double], standardDeviationX: Double, standardDeviationW: Double): LinearGaussianLatentFeatureModel = {
    new LinearGaussianLatentFeatureModel(X.copy, 1/(standardDeviationX*standardDeviationX), 1/(standardDeviationW*standardDeviationW))
  }

  def usingVariances(X: DenseMatrix[Double], varianceX: Double, varianceW: Double): LinearGaussianLatentFeatureModel = {
    new LinearGaussianLatentFeatureModel(X.copy, 1/varianceX, 1/varianceW)
  }

  def main(args: Array[String]): Unit = {
    val Z = DenseMatrix.zeros[Int](6,2)
    println(Z)
    println(3+4)
  }

}

