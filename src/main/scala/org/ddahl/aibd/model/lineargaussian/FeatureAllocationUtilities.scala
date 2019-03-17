package org.ddahl.aibd.model.lineargaussian

import org.ddahl.matrix._

object FeatureAllocationUtilities {

  def getData(Z: Matrix): Array[Array[Double]] = {
    Z match {
      case x : org.apache.commons.math3.linear.Array2DRowRealMatrix => x.getDataRef  // Avoid copying
      case x => x.getData                                                            // Copy if needed
    }
  }

  def pretty(Z: Matrix): String = {
    if ( Z == null ) "" else getData(Z).map(_.map { x => "%1.0f".format(x) }.mkString(" ")).mkString("\n")
  }

  def isValid(Z: Matrix): Boolean = {
    val data = getData(Z)
    data.forall(_.forall(x => ( x == 0.0 ) || ( x == 1.0 )))
  }

  def toFingerprint(Z: Matrix): Array[BigInt] = {
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

  def toFingerprintString(Z: Matrix): String = {
    toFingerprint(Z).mkString(",")
  }

  def fromFingerprint(x: Array[BigInt]): Matrix = {
    val Z = matrixOfDim(x(0).toInt, x.length-1)
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

  def fromFingerprintString(x: String): Matrix = {
    fromFingerprint(x.split(",").map(BigInt.apply))
  }

 def leftOrderedForm(Z: Matrix): Matrix = {
    val fp = toFingerprint(Z)
    val order = (0 until Z.cols).zip(fp.view(1,Z.cols+1)).sortWith( _._2 > _._2 ).map(_._1)
    Z(::,order)
  }

  def partitionBySingletonsOf(i: Int, Z: Matrix): (Matrix,Matrix) = {
    val data = getData(Z.t)
    val a = data.map { f =>
      ( f(i) != 0.0 ) && ( f.sum == 1.0 )
    }
    (Z(::,a), Z(::,a.negate))
  }

  def enumerateCombinationsFor(i: Int, Z: Matrix): Array[Matrix] = {
    if ( Z.cols > 31 ) throw new IllegalArgumentException("No more than 31 columns are supported.")
    val max = (math.pow(2, Z.cols) - 1).toInt
    var result = List[Matrix]()
    var w = 0
    while ( w <= max ) {
      val newZ = Z.copy
      newZ(i,::) = Array.tabulate(Z.cols)(j => if ( (w & (1<<j)) != 0 ) 1.0 else 0.0)
      result = newZ :: result
      w += 1
    }
    result.toArray
  }

  def main(args: Array[String]): Unit = {
    val m = Array(Array[Double](0,1,0,1),Array[Double](1,0,1,0),Array[Double](1,0,1,0),Array[Double](0,0,0,1),Array[Double](1,0,0,0))
    val Z = wrap(m)
    println(isValid(Z))
    Z(0,2) = 7.0
    println(isValid(Z))
    Z(0,2) = 1.0
    println(pretty(Z))
    println("---")
    println(pretty(leftOrderedForm(Z)))
    println("<<")
    val (a,b) = partitionBySingletonsOf(1,Z)
    println(pretty(a))
    println(pretty(b))
    println(">>")
    println(enumerateCombinationsFor(0,Z).map(pretty).mkString("\n---------\n"))
  }

}
*/
