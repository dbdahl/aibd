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
    val Z = matrixOfDim(rows, x.length)
    var j = 0
    while ( j < x.length ) {
      x(j).foreach { Z(_,j) = 1.0 }
      j += 1
    }
    Z
  }

  implicit val ordering = new Ordering[BitSet] {
    def compare(x: BitSet, y: BitSet): Int = {
      val xi = x.iterator
      val yi = y.iterator
      while (xi.hasNext && yi.hasNext) {
        val xn = xi.next()
        val yn = yi.next()
        if (xn < yn) return -1
        if (xn > yn) return 1
      }
      if (xi.hasNext) return 1
      if (yi.hasNext) return -1
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

  def enumerateCombinationsFor(i: Int, Zsingletons: Matrix, Zexisting: Matrix): Array[Matrix] = {
    val fpSingletons = toFingerprint(Zsingletons).toList
    val data = getData(Zexisting.copy)
    data(i) = Array.ofDim[Double](Zexisting.cols)
    val fp = toFingerprint(wrap(data))
    val a = fp.groupBy(identity).mapValues(_.size)
    val b = a.map { x =>
      val off = x._1
      val on  = x._1 + i
      val count = x._2
      List.tabulate(count+1) { n => List.tabulate(count) { k => if ( k < n ) on else off } }
    }
    val n = b.map(_.length).product
    var counter = 0
    val collector = Array.ofDim[List[BitSet]](n)
    def engine(toProcess: Iterable[List[List[BitSet]]], result: List[BitSet]): Unit = {
      if ( toProcess.isEmpty ) {
        collector(counter) = result
        counter += 1
      } else toProcess.head.foreach { h => engine(toProcess.tail, h ++ result) }
    }
    engine(b, fpSingletons)
    collector.map(x => fromFingerprint(x.toArray.sorted, Zexisting.rows))
  }

  def main(args: Array[String]): Unit = {
    val m = Array(Array[Double](0,1,0,1,1,1,0,1,0,0,1,1),Array[Double](1,1,1,0,0,1,1,1,1,1,0,0),Array[Double](1,1,1,0,1,0,0,0,0,0,0,0),Array[Double](0,0,0,1,0,1,1,1,0,1,0,1),Array[Double](0,0,0,1,0,0,1,0,1,0,1,1))
    val Z = wrap(m)
    println()
    println("-- original ---")
    println(pretty(Z))
    println("--")
    val Zs = enumerateCombinationsFor(0,null,Z)
    Zs.foreach { ZZ =>
      println(pretty(ZZ))
      println("--")
    }
    println(Zs.length)
    val fas = Zs.map { Z =>
      PosteriorSimulation.Z2fa(Z,Z.rows)
    }
    println(fas.length)
    println(fas.toSet.size)
  }

}
