package org.ddahl.aibd

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

}

object Similarity {

  def uniform(size: Int) = new Similarity() {
    val nItems = size
    val isUniform = true

    def apply(i: Int, j: Int) = if ( i != j ) 1.0 else 0.0
  }

  def apply(x: Array[Array[Double]]) = fromDistance(x, identity)

  def fromDistance(x: Array[Array[Double]], similarityFunction: Double => Double) = new Similarity() {
    if (x == null) throw new IllegalArgumentException("Array cannot be null.")
    val nItems = x.length
    if (x.length == 0) throw new IllegalArgumentException("Dimension must be at least 1x1.")
    if (x(0) == null || x(0).length != nItems) throw new IllegalArgumentException("Inconsistent dimensions.")
    private val y = new Array[Array[Double]](x.length)
    for (i <- x.indices) {
      y(i) = new Array[Double](i + 1)
      for (j <- 0 until i) {
        if (x(i)(j) == Double.PositiveInfinity || x(i)(j) == Double.NaN || x(i)(j) <= 0.0) throw new IllegalArgumentException("Distances must be positive.")
        if (x(i)(j) != x(j)(i)) throw new IllegalArgumentException("Distances must be symmetric.")
        y(i)(j) = similarityFunction(x(i)(j))
      }
    }
    val isUniform = false

    def apply(i: Int, j: Int) = if (j < i) y(i)(j) else y(j)(i)
  }

}

