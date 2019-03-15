package org.ddahl.aibd.model.lineargaussian

import org.ddahl.matrix._
import org.apache.commons.math3.util.FastMath.log

class LikelihoodComponents private[lineargaussian] (val Z: Matrix, val Zt: Matrix, val M: Matrix) {
  val K = Z.cols
}

class LinearGaussianLatentFeatureModel private (val X: Matrix, val precisionX: Double, val precisionW: Double) {

  val N = X.rows
  val D = X.cols

  private val Xt = X.t
  private val XtX = X.t * X
  private val traceXtX = trace(XtX)
  private val I = eye(N)
  private val Dhalf = D / 2.0
  private val const = -N * Dhalf * log(2 * math.Pi) - precisionX / 2 * traceXtX
  private val DhalfTimesLogPrecisionX = Dhalf * log(precisionX)
  private val DhalfTimesLogPrecisionW = Dhalf * log(precisionW)
  private val halfPrecisionX = precisionX / 2
  private val ratioOfPrecisions = precisionW / precisionX

  def logLikelihood(lc: LikelihoodComponents): Double = {
    if (lc.K == 0) {
      const + (N       ) * DhalfTimesLogPrecisionX
    } else {
      const + (N - lc.K) * DhalfTimesLogPrecisionX + lc.K * DhalfTimesLogPrecisionW + Dhalf * log(det(lc.M)) + halfPrecisionX * trace(Xt * lc.Z * lc.M * lc.Zt * X)
    }
  }

  def logLikelihood(Z: Matrix): Double = logLikelihood(computeLikelihoodComponents(Z))

  def computeLikelihoodComponents(Z: Matrix): LikelihoodComponents = {
    if (Z.rows != N) throw new IllegalArgumentException("Feature allocation has " + Z.rows + " items, but " + N + " were expected.")
    val Zt = Z.t
    val M = inv(Zt * Z + diag(Array.fill(Z.cols)(ratioOfPrecisions)))
    new LikelihoodComponents(Z, Zt, M)
  }

  private def update(M: Matrix, z: Array[Double], add: Boolean): Matrix = {
    M - M * z ** z * M / (z * M * z + ( if (add) 1 else -1))
  }

  def dropFeaturesFor(i: Int, lc: LikelihoodComponents): LikelihoodComponents = {
    val Z = lc.Z.copy
    val zOld = Z(i,::)
    Z(i,::) = Array.ofDim[Double](lc.K)
    new LikelihoodComponents(Z, Z.t, update(lc.M, zOld, false))
  }

  def addFeaturesFor(i: Int, lc: LikelihoodComponents, z: Array[Double]): LikelihoodComponents = {
    val Z = lc.Z.copy
    // assert(Z(i,::).forall(_ == 0.0))
    Z(i,::) = z
    new LikelihoodComponents(Z, Z.t, update(lc.M, z, true))
  }

  def dropAndAddFeaturesFor(i: Int, lc: LikelihoodComponents, z: Array[Double]): LikelihoodComponents = {
    val newLC = dropFeaturesFor(i, lc)
    newLC.Z(i,::) = z
    new LikelihoodComponents(newLC.Z, newLC.Z.t, update(newLC.M, z, true))
  }

}


object LinearGaussianLatentFeatureModel {

  def usingPrecisions(X: Matrix, precisionX: Double, precisionW: Double): LinearGaussianLatentFeatureModel = {
    new LinearGaussianLatentFeatureModel(X, precisionX, precisionW)
  }

  def usingStandardDeviations(X: Matrix, standardDeviationX: Double, standardDeviationW: Double): LinearGaussianLatentFeatureModel = {
    new LinearGaussianLatentFeatureModel(X, 1/(standardDeviationX*standardDeviationX), 1/(standardDeviationW*standardDeviationW))
  }

  def usingVariances(X: Matrix, varianceX: Double, varianceW: Double): LinearGaussianLatentFeatureModel = {
    new LinearGaussianLatentFeatureModel(X, 1/varianceX, 1/varianceW)
  }

  def usingPrecisions(X: Array[Array[Double]], precisionX: Double, precisionW: Double): LinearGaussianLatentFeatureModel = {
    usingPrecisions(wrap(X), precisionX, precisionW)
  }

  def usingStandardDeviations(X: Array[Array[Double]], standardDeviationX: Double, standardDeviationW: Double): LinearGaussianLatentFeatureModel = {
    usingStandardDeviations(wrap(X), standardDeviationX, standardDeviationW)
  }

  def usingVariances(X: Array[Array[Double]], varianceX: Double, varianceW: Double): LinearGaussianLatentFeatureModel = {
    usingVariances(wrap(X), varianceX, varianceW)
  }

}

