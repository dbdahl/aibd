package org.ddahl.aibd.model.lineargaussian

class Bob

/*
import breeze.linalg._

object FeatureAllocationUtilities {

  def arrays2Matrix(Z: Array[Array[Double]]): DenseMatrix[Double] = {
    DenseMatrix.create(Z.size, Z(0).size, Z.flatten, 0, Z(0).size, true)
  }

  def arrays2Matrix(Z: Array[Array[Int]]): DenseMatrix[Double] = {
    arrays2Matrix(Z.map(_.map(_.toDouble)))
  }

  def arrays2Matrix(Z: Array[Array[Boolean]]): DenseMatrix[Double] = {
    arrays2Matrix(Z.map(_.map(x => if ( x ) 1.0 else 0.0)))
  }

  def matrix2Arrays(Z: DenseMatrix[Double]): Array[Array[Double]] = {
    Array.tabulate(Z.rows) { i =>
      Z(i,::).t.toArray
    }
  }

  def isValid(Z: DenseMatrix[Double]): Boolean = {
    Z.forall(x => ( x == 0.0 ) || ( x == 1.0 ) )
  }

  def toFingerprint(Z: DenseMatrix[Double]): Array[BigInt] = {
    val coefficients = Array.fill[BigInt](Z.rows-1)(BigInt(2)).scan(BigInt(1))(_*_)
    val result = Array.fill(Z.cols+1)(BigInt(0))
    result(0) = Z.rows
    var j = 0
    while ( j < Z.cols ) {
      var sum = result(j+1)
      val a = Z(::,j).map(_ != 0.0)
      var i = 0
      while ( i < Z.rows ) {
        if ( a(i) ) sum += coefficients(Z.rows-i-1)
        i += 1
      }
      j += 1
      result(j) = sum
    }
    result
  }

  def toFingerprintString(Z: DenseMatrix[Double]): String = {
    toFingerprint(Z).mkString(",")
  }

  def fromFingerprint(x: Array[BigInt]): DenseMatrix[Double] = {
    val Z = DenseMatrix.zeros[Double](x(0).toInt, x.length-1)
    var j = 0
    while ( j < Z.cols ) {
      val w = x(j+1)
      var i = 0
      while ( i < Z.rows ) {
        Z(i,j) = if ( (w & (1<<i)) != 0 ) 1.0 else 0.0
        i += 1
      }
      j += 1
    }
    Z
  }

  def fromFingerprintString(x: String): DenseMatrix[Double] = {
    fromFingerprint(x.split(",").map(BigInt.apply))
  }

 def leftOrderedForm(Z: DenseMatrix[Double]): Matrix[Double] = {
    val fp = toFingerprint(Z)
    val order = (0 until Z.cols).zip(fp.view(1,Z.cols+1)).sortWith( _._2 > _._2 ).map(_._1)
    Z(::,order).toDenseMatrix
  }

  def partitionBySingletonsOf(i: Int, Z: DenseMatrix[Double]): (Matrix[Double],Matrix[Double]) = {
    val a = Z(::,*).map { f =>
      ( f(i) != 0.0 ) && ( sum(f) == 1.0 )
    }
    (Z(::,a.t), Z(::,!a.t))
  }

  def enumerateCombinationsFor(i: Int, Z: DenseMatrix[Double]): List[DenseMatrix[Double]] = {
    if ( Z.cols > 31 ) throw new IllegalArgumentException("No more than 31 columns are supported.")
    val max = (math.pow(2, Z.cols) - 1).toInt
    var result = List[DenseMatrix[Double]]()
    var w = 0
    while ( w <= max ) {
      val newZ = Z.copy
      newZ(i,::) := DenseVector(Array.tabulate(Z.cols)(j => if ( (w & (1<<j)) != 0 ) 1.0 else 0.0)).t
      result = newZ :: result
      w += 1
    }
    result
  }

  def main(args: Array[String]): Unit = {
    val m = Array(Array(0,1,0,1),Array(1,0,1,0),Array(1,0,1,0),Array(0,0,0,1),Array(1,0,0,0))
    val Z = arrays2Matrix(m)
    println(isValid(Z))
    Z(0,2) = 7.0
    println(isValid(Z))
    Z(0,2) = 1.0
    println(Z)
    println("---")
    println(leftOrderedForm(Z))
    println(partitionBySingletonsOf(0,Z))
    println(enumerateCombinationsFor(0,Z).mkString("\n---------\n"))

    // val part = partitionBySingletonsOf(0,Z)
    // val Zempty = part._2
    // Zempty(0,::) = DenseVector.zeros[Double](Z.cols).t
    // val lc = new LinearGaussianLatentFeatureModel(X,1,1)
    println(enumerateCombinationsFor(0,Z).map(ZZ => toFingerprintString(ZZ)))
  }

}
*/
