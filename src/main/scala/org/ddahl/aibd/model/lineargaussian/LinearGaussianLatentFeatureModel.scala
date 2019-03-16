package org.ddahl.aibd.model.lineargaussian

import org.ddahl.matrix._
import org.apache.commons.math3.util.FastMath.log

class LikelihoodComponents private[lineargaussian] (val Z: Matrix, val Zt: Matrix, val M: Matrix, val d: Double) {
  val K = Z.cols
}

class LinearGaussianLatentFeatureModel private (val X: Matrix, val precisionX: Double, val precisionW: Double) {

  val N = X.rows
  val D = X.cols

  private val Xt = X.t
  private val XtX = X.t * X
  private val traceXtX = trace(XtX)
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
      const + (N - lc.K) * DhalfTimesLogPrecisionX + lc.K * DhalfTimesLogPrecisionW + Dhalf * log(lc.d) + halfPrecisionX * trace(Xt * lc.Z * lc.M * lc.Zt * X)
    }
  }

  def logLikelihood(Z: Matrix): Double = logLikelihood(computeLikelihoodComponents(Z))

  def computeLikelihoodComponents(Z: Matrix): LikelihoodComponents = {
    if (Z.rows != N) throw new IllegalArgumentException("Feature allocation has " + Z.rows + " items, but " + N + " were expected.")
    val Zt = Z.t
    val W = Zt * Z + diag(Array.fill(Z.cols)(ratioOfPrecisions))
    val m = inv(W)
    val d = 1/det(W)
    new LikelihoodComponents(Z, Zt, m, d)
  }

  private def update(M: Matrix, d: Double, z: Array[Double], add: Boolean): (Matrix, Double) = {
    val sign= if (add) 1 else -1
    val zMz = z * M * z
    (M - M * z ** z * M / (zMz + sign), d / ( 1 + sign * zMz ))
  }

  def dropFeaturesFor(i: Int, lc: LikelihoodComponents): LikelihoodComponents = {
    val Z = lc.Z.copy
    val zOld = Z(i,::)
    Z(i,::) = Array.ofDim[Double](lc.K)
    val (m,d) = update(lc.M, lc.d, zOld, false)
    new LikelihoodComponents(Z, Z.t, m, d)
  }

  def addFeaturesFor(i: Int, lc: LikelihoodComponents, z: Array[Double]): LikelihoodComponents = {
    val Z = lc.Z.copy
    // assert(Z(i,::).forall(_ == 0.0))
    Z(i,::) = z
    val (m,d) = update(lc.M, lc.d, z, true)
    new LikelihoodComponents(Z, Z.t, m, d)
  }

  def dropAndAddFeaturesFor(i: Int, lc: LikelihoodComponents, z: Array[Double]): LikelihoodComponents = {
    val newLC = dropFeaturesFor(i, lc)
    newLC.Z(i,::) = z
    val (m,d) = update(newLC.M, newLC.d, z, true)
    new LikelihoodComponents(newLC.Z, newLC.Z.t, m, d)
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

