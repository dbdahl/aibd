package org.ddahl.aibd.model.lineargaussian

import org.ddahl.matrix._
import org.apache.commons.math3.linear.CholeskyDecomposition
import org.apache.commons.math3.util.FastMath.log

class LikelihoodComponents private[lineargaussian] (val Z: Matrix, val Zt: Matrix, val M: Matrix, val d: Double) {
  val K = if ( Z == null ) 0 else Z.cols
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
      const + (N - lc.K) * DhalfTimesLogPrecisionX + lc.K * DhalfTimesLogPrecisionW - Dhalf * log(lc.d) + halfPrecisionX * trace(Xt * lc.Z * lc.M * lc.Zt * X)
    }
  }

  def logLikelihood(lcs: Array[LikelihoodComponents]): Array[Double] = lcs.map(logLikelihood)

  def logLikelihood(Z: Matrix): Double = logLikelihood(computeLikelihoodComponents(Z))

  def logLikelihood(Zs: Array[Matrix]): Array[Double] = logLikelihood(computeLikelihoodComponents(Zs))

  def logLikelihood(featureAllocation: FeatureAllocation): Double = logLikelihood(computeLikelihoodComponents(featureAllocation))

  def logLikelihood(featureAllocations: Array[FeatureAllocation]): Array[Double] = logLikelihood(computeLikelihoodComponents(featureAllocations))

  def computeLikelihoodComponents(Z: Matrix): LikelihoodComponents = {
    if ( Z == null ) new LikelihoodComponents(null, null, null, 0.0) else {
      if (Z.rows != N) throw new IllegalArgumentException("Feature allocation has " + Z.rows + " items, but " + N + " were expected.")
      val Zt = Z.t
      val W = Zt * Z + diag(Array.fill(Z.cols)(ratioOfPrecisions))
      val chol = new CholeskyDecomposition(W)
      val m = chol.getSolver.getInverse
      val d = chol.getDeterminant
      new LikelihoodComponents(Z, Zt, m, d)
    }
  }

  def computeLikelihoodComponents(featureAllocation: FeatureAllocation): LikelihoodComponents = {
    computeLikelihoodComponents(wrap(featureAllocation.matrix))
  }

  def computeLikelihoodComponents(Zs: Array[Matrix]): Array[LikelihoodComponents] = Zs.map(computeLikelihoodComponents)

  def computeLikelihoodComponents(featureAllocations: Array[FeatureAllocation]): Array[LikelihoodComponents] = featureAllocations.map(computeLikelihoodComponents)

  private def update(M: Matrix, d: Double, z: Array[Double], add: Boolean): (Matrix, Double) = {
    val sign= if (add) 1 else -1
    val zMz = z * M * z
    (M - M * z ** z * M / (zMz + sign), d * ( 1 + sign * zMz ))
  }

  def dropFeaturesFor(i: Int, lc: LikelihoodComponents): LikelihoodComponents = {
    if ( lc.Z == null ) lc
    else {
      val Z = lc.Z.copy
      val zOld = Z(i, ::)
      Z(i, ::) = Array.ofDim[Double](lc.K)
      val (m, d) = update(lc.M, lc.d, zOld, false)
      new LikelihoodComponents(Z, Z.t, m, d)
    }
  }

  def addFeaturesFor(i: Int, lc: LikelihoodComponents, z: Array[Double]): LikelihoodComponents = {
    if ( lc.Z == null ) computeLikelihoodComponents(wrap(z))
    else {
      val Z = lc.Z.copy
      // assert(Z(i,::).forall(_ == 0.0))
      Z(i, ::) = z
      val (m, d) = update(lc.M, lc.d, z, true)
      new LikelihoodComponents(Z, Z.t, m, d)
    }
  }

  def dropAndAddFeaturesFor(i: Int, lc: LikelihoodComponents, z: Array[Double]): LikelihoodComponents = {
    if ( lc.Z == null ) computeLikelihoodComponents(wrap(z))
    else {
      val newLC = dropFeaturesFor(i, lc)
      newLC.Z(i, ::) = z
      val (m, d) = update(newLC.M, newLC.d, z, true)
      new LikelihoodComponents(newLC.Z, newLC.Z.t, m, d)
    }
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

