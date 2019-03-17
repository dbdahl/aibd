package org.ddahl

object matrix {

  import org.apache.commons.math3.linear.LUDecomposition

  type Matrix = org.apache.commons.math3.linear.RealMatrix

  def wrap(X: Array[Array[Double]]): Matrix = {
    if ( ( X.length == 0 ) || ( X(0).length == 0 ) ) null
    else new org.apache.commons.math3.linear.Array2DRowRealMatrix(X,false)
  }
  def eye(n: Int): Matrix = org.apache.commons.math3.linear.MatrixUtils.createRealIdentityMatrix(n)
  def diag(x: Array[Double]): Matrix = org.apache.commons.math3.linear.MatrixUtils.createRealDiagonalMatrix(x)
  def inv(X: Matrix): Matrix = org.apache.commons.math3.linear.MatrixUtils.inverse(X: Matrix)
  def trace(X: Matrix): Double = X.getTrace
  def det(X: Matrix): Double = new LUDecomposition(X).getDeterminant

  import scala.language.implicitConversions

  implicit def realMatrix2RichMatrix(X: Matrix): RichMatrix = new RichMatrix(X)
  implicit def array2RichArray(x: Array[Double]): RichArray = new RichArray(x)
  implicit def scalar2RichScalar(x: Double): RichDouble = new RichDouble(x)

}

import matrix._

class RichDouble(x: Double) {

  def *(Y: Matrix): Matrix = Y scalarMultiply x

}

class RichArray(x: Array[Double]) {

  def *(Y: Matrix): Array[Double] = Y preMultiply x

  def *(y: Double): Array[Double] = Array.tabulate(x.length)(y*x(_))

  def *(y: Array[Double]): Double = {
    if ( x.length != y.length ) throw new RuntimeException("Incompatible lengths.")
    var sum = 0.0
    var i = 0
    while ( i < x.length ) {
      sum += x(i)*y(i)
      i += 1
    }
    sum
  }

  def **(y: Array[Double]): Matrix = {
    val result = Array.ofDim[Double](x.length, y.length)
    var i = 0
    while ( i < x.length ) {
      val row = result(i)
      val xi = x(i)
      var j = 0
      while ( j < y.length) {
        row(j) = xi*y(j)
        j += 1
      }
      i += 1
    }
    wrap(result)
  }

}

class RichMatrix(X: Matrix) {

  def t: Matrix = X.transpose
  def trace: Double = X.getTrace
  def cols: Int = X.getColumnDimension
  def rows: Int = X.getRowDimension

  def *(Y: Matrix): Matrix = X multiply Y
  def *(y: Double): Matrix = X scalarMultiply y
  def *(y: Array[Double]): Array[Double] = X operate y

  def /(y: Double): Matrix = X scalarMultiply (1/y)

  def +(Y: Matrix): Matrix = X add Y
  def +(y: Double): Matrix = X scalarAdd y

  def -(Y: Matrix): Matrix = X subtract Y
  def -(y: Double): Matrix = X scalarAdd -y

  def apply(i: Int, j: Int): Double = X.getEntry(i,j)
  def apply(i: Int, j: scala.collection.immutable.::.type): Array[Double] = X.getRow(i)
  def apply(i: scala.collection.immutable.::.type, j: Int): Array[Double] = X.getColumn(j)

  def update(i: Int, j: Int, x: Double): Unit = X.setEntry(i,j,x)
  def update(i: Int, j: scala.collection.immutable.::.type, y: Array[Double]): Unit = X.setRow(i, y)
  def update(i: scala.collection.immutable.::.type, j: Int, y: Array[Double]): Unit = X.setColumn(j, y)

}

