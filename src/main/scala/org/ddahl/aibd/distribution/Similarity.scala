package org.ddahl.aibd.distribution

import org.apache.commons.math3.util.FastMath.{exp, pow}

trait Similarity {

  val nItems: Int

  val isUniform: Boolean

  def apply(i: Int, j: Int): Double

  def getData(): Array[Array[Double]] = {
    Array.tabulate(nItems) { i =>
      Array.tabulate(nItems) { j =>
        apply(i,j)
      }
    }
  }

  def toText(fmt: String) = {
    val r = new StringBuilder()
    for (i <- 0 until nItems) {
      for (j <- 0 until nItems) {
        r.append(apply(i, j).formatted(fmt) + " ")
      }
      r.append("\n")
    }
    r.dropRight(1).toString
  }

  override def toString: String = toText("%4.3f")

}

class FixedSimilarity(similarity: Array[Array[Double]]) extends Similarity {

  val nItems = similarity.length

  val isUniform = false

  def apply(i: Int, j: Int): Double = if (j < i) similarity(i)(j) else similarity(j)(i)

}

object Similarity {

  def apply(x: Array[Array[Double]]): Similarity = new FixedSimilarity(checkAndCopyMatrix(x))

  def checkAndCopyMatrix(x: Array[Array[Double]]): Array[Array[Double]] = {
    if (x == null) throw new IllegalArgumentException("Array cannot be null.")
    val nItems = x.length
    if (nItems == 0) throw new IllegalArgumentException("Dimension must be at least 1x1.")
    val y = new Array[Array[Double]](x.length)
    for (i <- x.indices) {
      if (x(i) == null || x(i).length != nItems) throw new IllegalArgumentException("Inconsistent dimensions.")
      y(i) = new Array[Double](i + 1)
      for (j <- 0 until i) {
        if (x(i)(j) == Double.PositiveInfinity || x(i)(j) == Double.NaN || x(i)(j) <= 0.0) throw new IllegalArgumentException("Distances must be positive.")
        if (x(i)(j) != x(j)(i)) throw new IllegalArgumentException("Distances must be symmetric.")
        y(i)(j) = x(i)(j)
      }
    }
    y
  }

  def checkTemperature(x: Double): Unit = {
    if ( x == Double.PositiveInfinity || x == Double.NaN || x < 0.0 ) throw new IllegalArgumentException("Temperature must be nonnegative.")
  }

}

trait HasTemperature[T] {

  val isUniform = false

  val temperature: Double

  def updateTemperature(temperature: Double): Similarity with HasTemperature[T]

}

class ExponentialSimilarity private (distances: Array[Array[Double]], val temperature: Double) extends Similarity with HasTemperature[ExponentialSimilarity] {

  val nItems = distances.length

  def apply(i: Int, j: Int): Double = if (j < i) y(i)(j) else y(j)(i)

  private val y = distances.map(_.map { z => exp(-temperature*z) })

  def updateTemperature(temperature: Double): ExponentialSimilarity = {
    new ExponentialSimilarity(distances, temperature)
  }

}

object ExponentialSimilarity {

  def apply(distances: Array[Array[Double]], temperature: Double): ExponentialSimilarity = {
    val d = Similarity.checkAndCopyMatrix(distances)
    Similarity.checkTemperature(temperature)
    new ExponentialSimilarity(d, temperature)
  }

}

class ReciprocalSimilarity private (distances: Array[Array[Double]], val temperature: Double) extends Similarity with HasTemperature[ReciprocalSimilarity] {

  val nItems = distances.length

  def apply(i: Int, j: Int): Double = if (j < i) y(i)(j) else y(j)(i)

  private val y = distances.map(_.map { z => pow(z,-temperature) })

  def updateTemperature(temperature: Double): ReciprocalSimilarity = {
    new ReciprocalSimilarity(distances, temperature)
  }

}

object ReciprocalSimilarity {

  def apply(distances: Array[Array[Double]], temperature: Double): ReciprocalSimilarity = {
    val d = Similarity.checkAndCopyMatrix(distances)
    Similarity.checkTemperature(temperature)
    new ReciprocalSimilarity(d, temperature)
  }

}

