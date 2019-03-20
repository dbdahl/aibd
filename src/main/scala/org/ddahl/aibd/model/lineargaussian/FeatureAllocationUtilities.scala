package org.ddahl.aibd.model.lineargaussian

import org.ddahl.matrix._
import scala.collection.mutable.BitSet

object FeatureAllocationUtilities {

  def toFingerprintString(Z: Matrix): String = {
    toFingerprint(Z).map(_.mkString(",")).mkString(";")
  }

  def fromFingerprintString(x: String, rows: Int): Matrix = {
    fromFingerprint(x.split(";").map( y => {
      new BitSet(y.split(",").map(_.toLong))
    }), rows)
  }

  def isValid(Z: Matrix): Boolean = {
    val data = getData(Z)
    data.forall(_.forall(x => ( x == 0.0 ) || ( x == 1.0 )))
  }

  def toFingerprint(Z: Matrix): Array[BitSet] = {
    if ( Z == null ) Array[BitSet]()
    else {
      val result = Array.fill(Z.cols)(BitSet())
      var j = 0
      while (j < Z.cols) {
        val bs = result(j)
        var i = 0
        while (i < Z.rows) {
          if (Z(i,j) != 0.0) bs.add(i)
          i += 1
        }
        result(j) = bs
        j += 1
      }
      result
    }
  }


  def fromFingerprint(x: Array[BitSet], rows: Int): Matrix = {
    val Z = matrixOfDim(rows, x.length-1)
    var j = 0
    while ( j < x.length ) {
      val w = x(j)
      var i = 0
      while ( i < rows ) {
        Z(i,j) = if ( w(i) ) 1.0 else 0.0
        i += 1
      }
      j += 1
    }
    Z
  }

  implicit val ordering = new Ordering[BitSet] {
    def compare(x: BitSet, y: BitSet): Int = {
      val xi = x.iterator
      val yi = x.iterator
      while (xi.hasNext && yi.hasNext) {
        val xn = xi.next()
        val yn = yi.next()
        if (xn < yn) return (-1)
        if (xn > yn) return (1)
      }
      if (xi.hasNext) return (1)
      if (yi.hasNext) return (-1)
      0
    }
  }

  def leftOrderedForm(Z: Matrix): Matrix = leftOrderedForm(toFingerprint(Z),Z.rows)

  def leftOrderedForm(fingerprint: Array[BitSet], rows: Int): Matrix = {
    fromFingerprint(fingerprint.sorted, rows)
  }

  def partitionBySingletonsOf(i: Int, Z: Matrix): (Matrix,Matrix) = {
    if ( Z == null) (null,null) else {
      val data = getData(Z.t)
      val a = data.map { f =>
        (f(i) == 1.0) && (f.sum == 1.0)
      }
      (Z(::, a), Z(::, a.negate))
    }
  }

  def enumerateCombinationsFor(i: Int, Z: Matrix): Array[Matrix] = {
    if ( Z == null ) new Array[Matrix](0) else {
      if (Z.cols > 31) throw new IllegalArgumentException("No more than 31 columns are supported.")
      val max = (math.pow(2, Z.cols) - 1).toInt
      var result = List[Matrix]()
      var w = 0
      while (w <= max) {
        val newZ = Z.copy
        newZ(i, ::) = Array.tabulate(Z.cols) { j => if ((w & (1 << j)) != 0) 1.0 else 0.0 }
        result = newZ :: result
        w += 1
      }
      result.toArray
    }
  }

  /*
  def enumerateCombinationsFor2(i: Int, Z: Matrix): Array[Matrix] = {
    val data = getData(Z.copy)
    data(i) = Array.ofDim[Double](Z.cols)
    val fp = toFingerprint(wrap(data))
    val a = fp.groupBy(identity).mapValues(_.size).toArray.map { x =>
      val off = x._1
      val on  = x._1 + i
      val count = x._2
      List.tabulate(count+1) { n =>
        List.fill(n){on} ++ List.fill(count-n){off}
      }
    }
    val n = a.map(_.size).product
    println(a.mkString(" ||| "))
    var all = List[List[BigInt]]()
    def engine(toProcess: Array[List[List[BigInt]]], result: List[BigInt]): Unit = {
      if ( toProcess.isEmpty ) all = result :: all
      else toProcess.head.foreach { h =>
        engine(toProcess.tail, h ++ result)
      }
    }
    engine(a, Nil)
    println(all.size)
    println(all)
    null
  }
  */

  /*
  def tabulateFeatures(Z: Matrix): List[(Array[Double],Int)] = {
    if ( Z == null ) null else {
      val fp = toFingerprint(Z)


      val features = getData(Z.t).map { feature =>
        var bs = BitSet()
        feature.index
        feature.map(x => if ( x == 1.0 ) bs(1) else 0L))
      }



      if (Z.cols > 31) throw new IllegalArgumentException("No more than 31 columns are supported.")
      val max = (math.pow(2, Z.cols) - 1).toInt
      var result = List[Matrix]()
      var w = 0
      while (w <= max) {
        val newZ = Z.copy
        newZ(i, ::) = Array.tabulate(Z.cols) { j => if ((w & (1 << j)) != 0) 1.0 else 0.0 }
        result = newZ :: result
        w += 1
      }
      result.toArray
    }
  }
  */

  def main(args: Array[String]): Unit = {
    val m = Array(Array[Double](0,1,0,1),Array[Double](1,1,1,0),Array[Double](1,1,1,0),Array[Double](0,0,0,1),Array[Double](0,0,0,1))
    val Z = wrap(m)
    println(pretty(Z))
    //println(enumerateCombinationsFor2(0,Z))
  }

}
