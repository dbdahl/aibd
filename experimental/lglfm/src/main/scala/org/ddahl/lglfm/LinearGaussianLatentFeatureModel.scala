package org.ddahl.lglfm

import breeze.linalg._
import org.apache.commons.math3.util.FastMath.log

class LinearGaussianLatentFeatureModel private (X: DenseMatrix[Double], precisionX: Double, precisionW: Double) {

  val N = X.rows
  val D = X.cols

  private val Xt = X.t
  private val I = DenseMatrix.eye[Double](N)
  private val Dhalf = D / 2.0
  private val const = -N * Dhalf * log(2 * math.Pi)
  private val XtXtrace = trace(Xt * X)

  def logLikelihood(Z: DenseMatrix[Double]): Double = {
    if (Z.rows != N) throw new IllegalArgumentException("Feature allocation has " + Z.rows + " items, but " + N + " were expected.")
    val K = Z.cols
    if (K == 0) {
      const + (N) * Dhalf * log(precisionX) + -precisionX / 2 * XtXtrace
    } else {
      val Zt = Z.t
      val ZtZplusRatioI = Zt * Z + I *:* (precisionW / precisionX)
      const + (N - K) * Dhalf * log(precisionX) + K * Dhalf * log(precisionW) - Dhalf * log(det(ZtZplusRatioI)) - precisionX / 2 * trace(Xt * (I - Z * inv(ZtZplusRatioI) * Zt) * X)
    }
  }

}

object LinearGaussianLatentFeatureModel {

  def main(args: Array[String]): Unit = {
    val Z = DenseMatrix.zeros[Int](6,2)
    println(Z)
    println(3+4)
  }

}

