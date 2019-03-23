package org.ddahl.aibd.model.lineargaussian

import org.ddahl.sdols.featureallocation.{FeatureAllocation => FeatureAllocationAlternative}
import org.ddahl.sdols.featureallocation.{Feature => FeatureAlternative}
import scala.collection.mutable.BitSet

sealed trait FeatureAllocation {

  val nItems: Int
  val nFeatures: Int
  val sizes: Array[Int]
  val array: Array[BitSet]
  val matrix: Array[Array[Double]]

  override def toString(): String = {
    if ( nItems == 0 ) "" else matrix.map(_.map { x => "%1.0f".format(x) }.mkString(" ")).mkString("\n")
  }

  def check(): Unit = {
    assert(nItems == matrix.length)
    assert(nFeatures == array.length)
    assert(nFeatures == sizes.length)
    assert( ( nFeatures == 0 ) || ( nFeatures == matrix(0).length ) )
    val that1 = FeatureAllocation(matrix)
    val that2 = new FeatureAllocationWithArray(nItems, array.map(_.clone))
    val that3 = new FeatureAllocationWithArrayAndSizes(nItems, array.map(_.clone), sizes.clone)
    val a = that1.matrix.flatten
    val b = that2.matrix.flatten
    val c = that3.matrix.flatten
    a.zip(b).forall { x => x._1 == x._2 }
    a.zip(c).forall { x => x._1 == x._2 }
    val aa = that1.array
    val bb = that2.array
    val cc = that3.array
    aa.zip(bb).forall { x => x._1 == x._2 }
    aa.zip(cc).forall { x => x._1 == x._2 }
    val aaa = that1.sizes
    val bbb = that2.sizes
    val ccc = that3.sizes
    aaa.zip(bbb).forall { x => x._1 == x._2 }
    aaa.zip(ccc).forall { x => x._1 == x._2 }
  }

  protected def computeSizes: Array[Int] = array.map(_.size)

  protected def computeArray: Array[BitSet] = {
    val result = Array.fill(nFeatures)(BitSet())
    var j = 0
    while (j < nFeatures) {
      val bs = result(j)
      var i = 0
      while (i < nItems) {
        if (matrix(i)(j) != 0.0) bs.add(i)
        i += 1
      }
      result(j) = bs
      bs.size
      j += 1
    }
    result
  }

  protected def computeMatrix: Array[Array[Double]] = {
    val result = Array.ofDim[Double](nItems, nFeatures)
    var j = 0
    while ( j < nFeatures ) {
      array(j).foreach { result(_)(j) = 1.0 }
      j += 1
    }
    result
  }

  def row(i: Int): Array[Double] = {
    if ( ( i < 0 ) || ( i >= nItems ) ) throw new IllegalArgumentException("Row index "+i+" is out of bounds [0"+(nItems-1)+"].")
    Array.tabulate(nFeatures) { j => if ( array(j)(i) ) 1.0 else 0.0 }
  }

  def add(i: Int): FeatureAllocation = {
    val newArray = array :+ BitSet(i)
    val newSizes = sizes :+ 1
    new FeatureAllocationWithArrayAndSizes(nItems, newArray, newSizes)
  }

  def remove(i: Int): FeatureAllocation = {
    if ( ( i < 0 ) || ( i >= nItems ) ) throw new IllegalArgumentException("Row index "+i+" is out of bounds [0"+(nItems-1)+"].")
    val newArray = array.map(_.clone)
    val newSizes = sizes.clone
    var nToDelete = 0
    var j = 0
    while ( j < nFeatures ) {
      if ( array(j)(i) ) {
        newSizes(j) -= 1
        newArray(j).remove(i)
        if ( newSizes(j) == 0 ) nToDelete += 1
      }
      j += 1
    }
    if ( nToDelete > 0 ) {
      val newArray2 = new Array[BitSet](nFeatures - nToDelete)
      val newSizes2 = new Array[Int](nFeatures - nToDelete)
      var j = 0
      var jj = 0
      while (jj < newSizes2.length) {
        while (newSizes(j) > 0) {
          newArray2(jj) = newArray(j)
          newSizes2(jj) = newSizes(j)
          j += 1
          jj += 1
        }
        j += 1
      }
      new FeatureAllocationWithArrayAndSizes(nItems, newArray2, newSizes2)
    } else {
      new FeatureAllocationWithArrayAndSizes(nItems, newArray, newSizes)
    }
  }

  def matrixWithout(i: Int): Array[Array[Double]] = {
    if ( ( i < 0 ) || ( i >= nItems ) ) throw new IllegalArgumentException("Row index "+i+" is out of bounds [0"+(nItems-1)+"].")
    val newArray = array.map(_.clone)
    var j = 0
    while ( j < nFeatures ) {
      if ( array(j)(i) ) newArray(j).remove(i)
      j += 1
    }
    new FeatureAllocationWithArray(nItems, newArray).matrix
  }

  def partitionBySingletonsOf(i: Int): (FeatureAllocation, FeatureAllocation) = {
    var sum = 0
    val sel = Array.tabulate(nFeatures) { j =>
      val result = ( sizes(j) == 1 ) && array(j)(i)
      if ( result ) sum += 1
      result
    }
    val leftArray = Array.ofDim[BitSet](sum)
    val leftSizes = Array.ofDim[Int](sum)
    val rightArray = Array.ofDim[BitSet](nFeatures-sum)
    val rightSizes = Array.ofDim[Int](nFeatures-sum)
    var j = 0
    var jl = 0
    var jr = 0
    while ( j < nFeatures ) {
      if ( sel(j) ) {
        leftArray(jl) = array(j).clone
        leftSizes(jl) = sizes(j)
        jl += 1
      } else {
        rightArray(jr) = array(j).clone
        rightSizes(jr) = sizes(j)
        jr += 1
      }
      j += 1
    }
    (new FeatureAllocationWithArrayAndSizes(nItems,leftArray,leftSizes), new FeatureAllocationWithArrayAndSizes(nItems,rightArray,rightSizes))
  }

  def enumerateCombinationsFor(i: Int): Array[FeatureAllocation] = {   // Careful, this results in tons of shared mutable instances.
    val (singletons, existing) = partitionBySingletonsOf(i)
    val a = existing.remove(i).array.groupBy(identity).mapValues(_.size)
    val b = a.map { x =>
      val off = x._1
      val on  = x._1 + i  // Leads to a clone, unlike the "add" method.
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
    engine(b, singletons.array.toList)
    if ( collector.length == 0 ) Array[FeatureAllocation]()
    else {
      val zeroedOutMatrix = new FeatureAllocationWithArray(nItems, collector.head.toArray).matrixWithout(i)
      collector.map { x =>
        val arr = x.toArray
        val mat = zeroedOutMatrix.clone // Shallow copy
        mat(i) = arr.map { f => if (f(i)) 1.0 else 0.0 }
        new FeatureAllocationWithMatrixAndArray(mat, arr)
      }
    }
  }

  def convertToAlternativeImplementation: FeatureAllocationAlternative[Null] = {
    FeatureAllocationAlternative(nItems, array.map { f =>
      FeatureAlternative(f.toArray:_*)
    }:_*)
  }

}

sealed class FeatureAllocationEmpty private[lineargaussian] (val nItems: Int) extends FeatureAllocation {

  val nFeatures = 0
  val sizes = Array[Int]()
  val array = Array[BitSet]()
  val matrix = Array.ofDim[Double](nItems,0)

}

sealed class FeatureAllocationWithMatrix private[lineargaussian] (val matrix: Array[Array[Double]]) extends FeatureAllocation {

  val nItems = matrix.length
  val nFeatures = matrix(0).length
  lazy val sizes = computeSizes
  lazy val array = computeArray

}

sealed class FeatureAllocationWithArray private[lineargaussian] (val nItems: Int, val array: Array[BitSet]) extends FeatureAllocation {

  val nFeatures = array.length
  lazy val sizes = computeSizes
  lazy val matrix = computeMatrix

}

sealed class FeatureAllocationWithArrayAndSizes private[lineargaussian] (val nItems: Int, val array: Array[BitSet], val sizes: Array[Int]) extends FeatureAllocation {

  val nFeatures = array.length
  lazy val matrix = computeMatrix

}

sealed class FeatureAllocationWithMatrixAndArray private[lineargaussian] (val matrix: Array[Array[Double]], val array: Array[BitSet]) extends FeatureAllocation {

  val nItems = matrix.length
  val nFeatures = array.length
  lazy val sizes = computeSizes

}

sealed class FeatureAllocationWithAll private[lineargaussian] (val matrix: Array[Array[Double]], val array: Array[BitSet], val sizes: Array[Int]) extends FeatureAllocation {

  val nItems = matrix.length
  val nFeatures = array.length

}

object FeatureAllocation {

  def apply(nItems: Int): FeatureAllocation = {
    if ( nItems < 0 ) throw new IllegalArgumentException("Number of items must be at least 0.")
    new FeatureAllocationEmpty(nItems)
  }

  def apply(matrix: Array[Array[Double]]): FeatureAllocation = {
    val rows = matrix.length
    if ( rows == 0 ) return new FeatureAllocationEmpty(rows)
    val cols = matrix(0).length
    if ( ! matrix.forall(_.length == cols) ) throw new IllegalArgumentException("Number of features must be at least 1.")
    if (cols == 0) return new FeatureAllocationEmpty(rows)
    if ( ! matrix.forall(_.forall(x => ( x == 0.0 ) || ( x == 1.0 ))) ) throw new IllegalArgumentException("Elements should be either 0 or 1.")
    new FeatureAllocationWithMatrix(matrix)
  }

  def apply(fa: FeatureAllocation): FeatureAllocation = {
    fa match {
      case e: FeatureAllocationEmpty => new FeatureAllocationEmpty(e.nItems)
      case e: FeatureAllocationWithMatrix => new FeatureAllocationWithMatrix(e.matrix.map(_.clone))
      case e: FeatureAllocationWithArray => new FeatureAllocationWithArray(e.nItems, e.array.map(_.clone))
      case e: FeatureAllocationWithArrayAndSizes => new FeatureAllocationWithArrayAndSizes(e.nItems, e.array.map(_.clone), e.sizes.clone)
      case e: FeatureAllocationWithMatrixAndArray => new FeatureAllocationWithMatrixAndArray(e.matrix.map(_.clone), e.array.map(_.clone))
      case e: FeatureAllocationWithAll => new FeatureAllocationWithAll(e.matrix.map(_.clone), e.array.map(_.clone), e.sizes.clone)
    }
  }

}

