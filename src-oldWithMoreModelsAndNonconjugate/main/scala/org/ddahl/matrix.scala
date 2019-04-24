package org.ddahl

object matrix {

  import org.apache.commons.math3.linear.LUDecomposition

  type Matrix = org.apache.commons.math3.linear.RealMatrix

  def matrixOfDim(rows: Int, cols: Int): Matrix = {
    if ( ( rows == 0 ) || ( cols == 0 ) ) null
    else new org.apache.commons.math3.linear.Array2DRowRealMatrix(rows,cols)
  }

  def wrap(X: Array[Array[Double]]): Matrix = {
    if ( ( X.length == 0 ) || ( X(0).length == 0 ) ) null
    else new org.apache.commons.math3.linear.Array2DRowRealMatrix(X,false)
  }

  def wrap(X: Array[Array[Int]]): Matrix = {
    if ( ( X.length == 0 ) || ( X(0).length == 0 ) ) null
    else new org.apache.commons.math3.linear.Array2DRowRealMatrix(X.map(_.map(_.toDouble)),false)
  }

  def wrap(x: Array[Double]): Matrix = {
    if ( ( x.length == 0 ) ) null
    else new org.apache.commons.math3.linear.Array2DRowRealMatrix(x)
  }

  def eye(n: Int): Matrix = org.apache.commons.math3.linear.MatrixUtils.createRealIdentityMatrix(n)
  def diag(x: Array[Double]): Matrix = org.apache.commons.math3.linear.MatrixUtils.createRealDiagonalMatrix(x)
  def inv(X: Matrix): Matrix = org.apache.commons.math3.linear.MatrixUtils.inverse(X: Matrix)
  def trace(X: Matrix): Double = X.getTrace
  def det(X: Matrix): Double = new LUDecomposition(X).getDeterminant

  def getData(X: Matrix): Array[Array[Double]] = {
    X match {
      case x : org.apache.commons.math3.linear.Array2DRowRealMatrix => x.getDataRef  // Avoid copying
      case x => x.getData                                                            // Copy if needed
    }
  }

  def pretty(X: Matrix): String = {
    if ( X == null ) "" else getData(X).map(_.map { x => "%1.0f".format(x) }.mkString(" ")).mkString("\n")
  }

  import scala.language.implicitConversions

  implicit def realMatrix2RichMatrix(X: Matrix): RichMatrix = new RichMatrix(X)
  implicit def array2RichDoubleArray(x: Array[Double]): RichDoubleArray = new RichDoubleArray(x)
  implicit def array2RichBooleanArray(x: Array[Boolean]): RichBooleanArray = new RichBooleanArray(x)
  implicit def scalar2RichScalar(x: Double): RichDouble = new RichDouble(x)

}

import matrix._

class RichDouble(x: Double) {

  def *(Y: Matrix): Matrix = Y scalarMultiply x

}

class RichBooleanArray(x: Array[Boolean]) {

  def negate: Array[Boolean] = x map(!_)

}

class RichDoubleArray(x: Array[Double]) {

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

  def |(y: Array[Double]): Matrix = {
    if ( y == null ) X
    else if ( X == null ) new org.apache.commons.math3.linear.Array2DRowRealMatrix(y)
    else {
      var i = 0
      wrap(getData(X).map { row =>
        val w = Array.ofDim[Double](row.length+1)
        Array.copy(row,0,w,0,row.length)
        w(row.length) = y(i)
        i += 1
        w
      })
    }
  }

  def |(Y: Matrix): Matrix = {
    if ( Y == null ) X
    else if ( X == null ) Y
    else {
      wrap(getData(X).zip(getData(Y)).map { tuple =>
        val w = Array.ofDim[Double](tuple._1.length + tuple._2.length)
        Array.copy(tuple._1, 0, w, 0, tuple._1.length)
        Array.copy(tuple._2, 0, w, tuple._1.length, tuple._2.length)
        w
      })
    }
  }

  def apply(i: Int, j: Int): Double = X.getEntry(i,j)
  def apply(i: Int, j: scala.collection.immutable.::.type): Array[Double] = X.getRow(i)
  def apply(i: IndexedSeq[Int], j: scala.collection.immutable.::.type)(implicit d: DummyImplicit): Matrix = wrap(i.map(X.getRow).toArray)
  def apply(i: IndexedSeq[Boolean], j: scala.collection.immutable.::.type)(implicit d1: DummyImplicit, d2: DummyImplicit): Matrix = {
    wrap(i.zipWithIndex.withFilter(_._1).map(k => X.getRow(k._2)).toArray)
  }
  def apply(i: scala.collection.immutable.::.type, j: Int): Array[Double] = X.getColumn(j)
  def apply(i: scala.collection.immutable.::.type, j: IndexedSeq[Int])(implicit d: DummyImplicit): Matrix = {
    val Xt = wrap(j.map(X.getColumn).toArray)
    if ( Xt == null ) null else Xt.transpose
  }
  def apply(i: scala.collection.immutable.::.type, j: IndexedSeq[Boolean])(implicit d1: DummyImplicit, d2: DummyImplicit): Matrix = {
    val Xt = wrap(j.zipWithIndex.withFilter(_._1).map(k => X.getColumn(k._2)).toArray)
    if ( Xt == null ) null else Xt.transpose
  }

  def update(i: Int, j: Int, x: Double): Unit = X.setEntry(i,j,x)
  def update(i: Int, j: scala.collection.immutable.::.type, y: Array[Double]): Unit = X.setRow(i, y)
  def update(i: scala.collection.immutable.::.type, j: Int, y: Array[Double]): Unit = X.setColumn(j, y)

}

