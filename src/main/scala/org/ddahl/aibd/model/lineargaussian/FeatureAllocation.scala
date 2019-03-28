package org.ddahl.aibd.model.lineargaussian

import org.ddahl.sdols.featureallocation.{FeatureAllocation => FeatureAllocationAlternative}
import org.ddahl.sdols.featureallocation.{Feature => FeatureAlternative}
import scala.collection.mutable.BitSet

sealed trait FeatureAllocation {

  val nItems: Int
  val nFeatures: Int
  val sizes: Array[Int]
  val features: Array[BitSet]
  val featuresAsList: Array[List[Int]]
  val matrix: Array[Array[Double]]

  protected var matrixIsCached: Boolean = false

  override def toString(): String = {
    if (nItems == 0) "" else matrix.map(_.map { x => "%1.0f".format(x) }.mkString(" ")).mkString("\n")
  }

  override def equals(that: Any): Boolean = {
    that match {
      case that: FeatureAllocation => ( that.nItems == this.nItems ) &&
                                      ( that.nFeatures == this.nFeatures ) &&
                                      that.sizes.zip(this.sizes).forall(x => x._1 == x._2) &&
                                      that.features.zip(this.features).forall(x => x._1 == x._2)
      case _ => false
    }
  }

  def check(): Unit = {
    assert(nItems == matrix.length)
    assert(nFeatures == features.length)
    assert(nFeatures == sizes.length)
    assert( ( nFeatures == 0 ) || ( nFeatures == matrix(0).length ) )
    val that1 = FeatureAllocation(matrix)
    val that2 = new FeatureAllocationWithArray(nItems, features.map(_.clone))
    val that3 = new FeatureAllocationWithArrayAndSizes(nItems, features.map(_.clone), sizes.clone)
    val a = that1.matrix.flatten
    val b = that2.matrix.flatten
    val c = that3.matrix.flatten
    a.zip(b).forall { x => x._1 == x._2 }
    a.zip(c).forall { x => x._1 == x._2 }
    val aa = that1.features
    val bb = that2.features
    val cc = that3.features
    aa.zip(bb).forall { x => x._1 == x._2 }
    aa.zip(cc).forall { x => x._1 == x._2 }
    val aaa = that1.sizes
    val bbb = that2.sizes
    val ccc = that3.sizes
    aaa.zip(bbb).forall { x => x._1 == x._2 }
    aaa.zip(ccc).forall { x => x._1 == x._2 }
  }

  protected def computeSizes: Array[Int] = features.map(_.size)

  protected def computeFeatures: Array[BitSet] = {
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

  protected def computeFeaturesAsList: Array[List[Int]] = features.map(_.toList)

  protected def computeMatrix: Array[Array[Double]] = {
    matrixIsCached = true
    val result = Array.ofDim[Double](nItems, nFeatures)
    var j = 0
    while ( j < nFeatures ) {
      features(j).foreach { result(_)(j) = 1.0 }
      j += 1
    }
    result
  }

  def matrixRow(i: Int): Array[Double] = {
    if ( ( i < 0 ) || ( i >= nItems ) ) throw new IllegalArgumentException("Item index "+i+" is out of bounds [0"+(nItems-1)+"].")
    if ( matrixIsCached ) matrix(i) else Array.tabulate(nFeatures) { j => if ( features(j)(i) ) 1.0 else 0.0 }
  }

  def matrixColumn(j: Int): Array[Double] = {
    if ( ( j < 0 ) || ( j >= nFeatures ) ) throw new IllegalArgumentException("Feature index "+j+" is out of bounds [0"+(nFeatures-1)+"].")
    val f = features(j)
    Array.tabulate(nItems) { i => if ( f(i) ) 1.0 else 0.0 }
  }

  def itemsOf(j: Int): Array[Int] = {
    if ( ( j < 0 ) || ( j >= nFeatures ) ) throw new IllegalArgumentException("Feature index "+j+" is out of bounds [0"+(nFeatures-1)+"].")
    val f = features(j)
    val result = new Array[Int](sizes(j))
    var ii = 0
    var i = 0
    while ( ii < result.size ) {
      if ( f(i) ) {
        result(ii) = i
        ii += 1
      }
      i += 1
    }
    result
  }

  def featuresOf(i: Int): Array[Int] = {
    if ( ( i < 0 ) || ( i >= nItems ) ) throw new IllegalArgumentException("Item index "+i+" is out of bounds [0"+(nItems-1)+"].")
    var list = List[Int]()
    var size = 0
    var j = 0
    while ( j < nFeatures ) {
      if ( features(j)(i) ) {
        size += 1
        list = j :: list
      }
      j += 1
    }
    val result = new Array[Int](size)
    while ( size > 0 ) {
      size -= 1
      result(size) = list.head
      list = list.tail
    }
    result
  }

  def featuresOf(i: Array[Int]): Array[Array[Int]] = i.map(featuresOf)

  def add(i: Int): FeatureAllocation = {
    val newArray = features :+ BitSet(i)  // Shallow copy
    val newSizes = sizes :+ 1
    new FeatureAllocationWithArrayAndSizes(nItems, newArray, newSizes)
  }

  def add(i: Int, j: Int): FeatureAllocation = {
    if ( ( i < 0 ) || ( i >= nItems ) ) throw new IllegalArgumentException("Item index "+i+" is out of bounds [0"+(nItems-1)+"].")
    if ( ( j < 0 ) || ( j >= nFeatures ) ) throw new IllegalArgumentException("Feature index "+j+" is out of bounds [0"+(nFeatures-1)+"].")
    if ( features(j)(i) ) this
    else {
      val newArray = features.clone  // Shallow copy
      val newSizes = sizes.clone
      newArray(j) = newArray(j) + i  // Clones
      newSizes(j) += 1
      if ( matrixIsCached ) {
        val newMatrix = matrix.clone  // Shallow copy
        newMatrix(i) = new Array[Double](nFeatures)
        newMatrix(i)(j) = 1.0
        new FeatureAllocationWithAll(newMatrix, newArray, newSizes)
      } else {
        new FeatureAllocationWithArrayAndSizes(nItems, newArray, newSizes)
      }
    }
  }

  def mutateAdd(i: Int, j: Int): Unit = {
    if ( ( i < 0 ) || ( i >= nItems ) ) throw new IllegalArgumentException("Item index "+i+" is out of bounds [0"+(nItems-1)+"].")
    if ( ( j < 0 ) || ( j >= nFeatures ) ) throw new IllegalArgumentException("Feature index "+j+" is out of bounds [0"+(nFeatures-1)+"].")
    if ( ! features(j)(i) ) {
      sizes(j) += 1
      featuresAsList(j) = i :: featuresAsList(j)  // Lazy, so must be before next line!
      features(j).add(i)                             // Mutates
      if ( matrixIsCached ) matrix(i)(j) = 1.0
    }
  }

  def remove(i: Int): FeatureAllocation = remove(Array(i), false)

  def remove(i: Int, cloneAndKeepEmptyFeatures: Boolean): FeatureAllocation = remove(Array(i), cloneAndKeepEmptyFeatures)

  def remove(i: Array[Int]): FeatureAllocation = remove(i, false)

  def remove(i: Array[Int], cloneAndKeepEmptyFeatures: Boolean): FeatureAllocation = {
    val newArray = if ( cloneAndKeepEmptyFeatures ) features.map(_.clone) else features.clone
    val newSizes = sizes.clone
    var nToDelete = 0
    var iii = 0
    while ( iii < i.length ) {
      val ii = i(iii)
      if ( ( ii < 0 ) || ( ii >= nItems ) ) throw new IllegalArgumentException("Item index "+ii+" is out of bounds [0"+(nItems-1)+"].")
      var j = 0
      while (j < nFeatures) {
        if ( ( newSizes(j) > 0 ) && features(j)(ii) ) {
          if ( cloneAndKeepEmptyFeatures ) newArray(j).remove(ii) else newArray(j) = newArray(j) - ii
          newSizes(j) -= 1
          if (newSizes(j) == 0) nToDelete += 1
        }
        j += 1
      }
      iii += 1
    }
    if ( ( !cloneAndKeepEmptyFeatures ) && ( nToDelete > 0 ) ) {
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
      if ( matrixIsCached ) {
        val newMatrix = if ( cloneAndKeepEmptyFeatures ) matrix.map(_.clone) else matrix.clone
        i.foreach { newMatrix(_) = new Array[Double](nFeatures) }
        new FeatureAllocationWithAll(newMatrix, newArray, newSizes)
      } else {
        new FeatureAllocationWithArrayAndSizes(nItems, newArray, newSizes)
      }
    }
  }

  def matrixWithout(i: Int): Array[Array[Double]] = {
    if ( ( i < 0 ) || ( i >= nItems ) ) throw new IllegalArgumentException("Item index "+i+" is out of bounds [0"+(nItems-1)+"].")
    val newMatrix = matrix.clone  // Shallow copy
    newMatrix(i) = new Array[Double](nFeatures)
    newMatrix
  }

  def partitionBySingletonsOf(i: Int): (FeatureAllocation, FeatureAllocation) = {
    var sum = 0
    val sel = Array.tabulate(nFeatures) { j =>
      val result = ( sizes(j) == 1 ) && features(j)(i)
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
        leftArray(jl) = features(j).clone
        leftSizes(jl) = sizes(j)
        jl += 1
      } else {
        rightArray(jr) = features(j).clone
        rightSizes(jr) = sizes(j)
        jr += 1
      }
      j += 1
    }
    (new FeatureAllocationWithArrayAndSizes(nItems,leftArray,leftSizes), new FeatureAllocationWithArrayAndSizes(nItems,rightArray,rightSizes))
  }

  def enumerateCombinationsFor(i: Int): Array[FeatureAllocation] = {   // Careful, this results in tons of shared mutable instances.
    val (singletons, existing) = partitionBySingletonsOf(i)
    val a = existing.remove(i).features.groupBy(identity).mapValues(_.size)
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
    engine(b, singletons.features.toList)
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
    FeatureAllocationAlternative(nItems, features.map { f =>
      FeatureAlternative(f.toArray:_*)
    }:_*)
  }

}

sealed class FeatureAllocationEmpty private[lineargaussian] (override val nItems: Int) extends FeatureAllocation {

  override val nFeatures = 0
  override val sizes = Array[Int]()
  override val features = Array[BitSet]()
  override lazy val featuresAsList = computeFeaturesAsList
  override val matrix = Array.ofDim[Double](nItems,0)
  matrixIsCached = true

}

sealed class FeatureAllocationWithMatrix private[lineargaussian] (override val matrix: Array[Array[Double]]) extends FeatureAllocation {

  override val nItems = matrix.length
  override val nFeatures = matrix(0).length
  override lazy val sizes = computeSizes
  override lazy val features = computeFeatures
  override lazy val featuresAsList = computeFeaturesAsList
  matrixIsCached = true

}

sealed class FeatureAllocationWithArray private[lineargaussian] (override val nItems: Int, override val features: Array[BitSet]) extends FeatureAllocation {

  override val nFeatures = features.length
  override lazy val sizes = computeSizes
  override lazy val featuresAsList = computeFeaturesAsList
  override lazy val matrix = computeMatrix

}

sealed class FeatureAllocationWithArrayAndSizes private[lineargaussian] (override val nItems: Int, override val features: Array[BitSet], override val sizes: Array[Int]) extends FeatureAllocation {

  override val nFeatures = features.length
  override lazy val featuresAsList = computeFeaturesAsList
  override lazy val matrix = computeMatrix

}

sealed class FeatureAllocationWithMatrixAndArray private[lineargaussian] (override val matrix: Array[Array[Double]], override val features: Array[BitSet]) extends FeatureAllocation {

  override val nItems = matrix.length
  override val nFeatures = features.length
  override lazy val sizes = computeSizes
  override lazy val featuresAsList = computeFeaturesAsList
  matrixIsCached = true

}

sealed class FeatureAllocationWithAll private[lineargaussian] (override val matrix: Array[Array[Double]], override val features: Array[BitSet], override val sizes: Array[Int]) extends FeatureAllocation {

  override val nItems = matrix.length
  override val nFeatures = features.length
  override lazy val featuresAsList = computeFeaturesAsList
  matrixIsCached = true

}

object FeatureAllocation {

  def apply(nItems: Int): FeatureAllocation = {
    if ( nItems < 0 ) throw new IllegalArgumentException("Number of items must be at least 0.")
    new FeatureAllocationEmpty(nItems)
  }

  def apply(matrix: Array[Array[Int]]): FeatureAllocation = apply(matrix.map(_.map(_.toDouble)))

  def apply(matrix: Array[Array[Double]]): FeatureAllocation = {
    val rows = matrix.length
    if ( rows == 0 ) return new FeatureAllocationEmpty(rows)
    val cols = matrix(0).length
    if ( ! matrix.forall(_.length == cols) ) throw new IllegalArgumentException("Number of features must be at least 1.")
    if (cols == 0) return new FeatureAllocationEmpty(rows)
    if ( ! matrix.forall(_.forall(x => ( x == 0.0 ) || ( x == 1.0 ))) ) throw new IllegalArgumentException("Elements should be either 0 or 1.")
    new FeatureAllocationWithMatrix(matrix)
  }

  def apply(fa: FeatureAllocation): FeatureAllocation = {  // Shallow clones
    fa match {
      case e: FeatureAllocationEmpty => new FeatureAllocationEmpty(e.nItems)
      case e: FeatureAllocationWithMatrix => new FeatureAllocationWithMatrix(e.matrix.clone)
      case e: FeatureAllocationWithArray => new FeatureAllocationWithArray(e.nItems, e.features.clone)
      case e: FeatureAllocationWithArrayAndSizes => new FeatureAllocationWithArrayAndSizes(e.nItems, e.features.clone, e.sizes.clone)
      case e: FeatureAllocationWithMatrixAndArray => new FeatureAllocationWithMatrixAndArray(e.matrix.clone, e.features.clone)
      case e: FeatureAllocationWithAll => new FeatureAllocationWithAll(e.matrix.clone, e.features.clone, e.sizes.clone)
    }
  }

}

