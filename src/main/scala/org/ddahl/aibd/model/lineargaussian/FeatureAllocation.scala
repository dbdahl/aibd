package org.ddahl.aibd.model.lineargaussian

import org.ddahl.sdols.featureallocation.{FeatureAllocation => FeatureAllocationAlternative}
import org.ddahl.sdols.featureallocation.{Feature => FeatureAlternative}
import org.ddahl.aibd.Utils.logFactorial
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
      case that: FeatureAllocation => (that.nItems == this.nItems) &&
        (that.nFeatures == this.nFeatures) &&
        that.sizes.zip(this.sizes).forall(x => x._1 == x._2) &&
        that.features.zip(this.features).forall(x => x._1 == x._2)
      case _ => false
    }
  }

  def check(): Unit = {
    assert(nItems == matrix.length)
    assert(nFeatures == features.length)
    assert(nFeatures == sizes.length)
    assert((nFeatures == 0) || (nFeatures == matrix(0).length))
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

  def isEmpty(j: Int): Boolean = sizes(j) == 0

  def matrixRow(i: Int): Array[Double] = {
    if ( ( i < 0 ) || ( i >= nItems ) ) throw new IllegalArgumentException("Item index "+i+" is out of bounds [0,"+(nItems-1)+"].")
    if ( matrixIsCached ) matrix(i) else Array.tabulate(nFeatures) { j => if ( features(j)(i) ) 1.0 else 0.0 }
  }

  def matrixColumn(j: Int): Array[Double] = {
    if ( ( j < 0 ) || ( j >= nFeatures ) ) throw new IllegalArgumentException("Feature index "+j+" is out of bounds [0,"+(nFeatures-1)+"].")
    val f = features(j)
    Array.tabulate(nItems) { i => if ( f(i) ) 1.0 else 0.0 }
  }

  def itemsOf(j: Int): Array[Int] = {
    if ( ( j < 0 ) || ( j >= nFeatures ) ) throw new IllegalArgumentException("Feature index "+j+" is out of bounds [0,"+(nFeatures-1)+"].")
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
    if ( ( i < 0 ) || ( i >= nItems ) ) throw new IllegalArgumentException("Item index "+i+" is out of bounds [0,"+(nItems-1)+"].")
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

  def add(fa: FeatureAllocation): FeatureAllocation = {
    if ( fa.nFeatures == 0 ) this
    else {
      val newFeatures = new Array[BitSet](nFeatures + fa.nFeatures)
      Array.copy(   features, 0, newFeatures, 0,            nFeatures)
      Array.copy(fa.features, 0, newFeatures, nFeatures, fa.nFeatures)
      new FeatureAllocationWithArray(nItems, newFeatures)
    }
  }

  def add(i: Int): FeatureAllocation = {
    val newFeatures = features :+ BitSet(i)  // Shallow copy
    val newSizes = sizes :+ 1
    new FeatureAllocationWithArrayAndSizes(nItems, newFeatures, newSizes)
  }

  def addOld(i: Int, j: Int): FeatureAllocation = {
    if ( ( i < 0 ) || ( i >= nItems ) ) throw new IllegalArgumentException("Item index "+i+" is out of bounds [0,"+(nItems-1)+"].")
    if ( ( j < 0 ) || ( j >= nFeatures ) ) throw new IllegalArgumentException("Feature index "+j+" is out of bounds [0,"+(nFeatures-1)+"].")
    if ( features(j)(i) ) this
    else {
      val newFeatures = features.clone  // Shallow copy
      val newSizes = sizes.clone
      newFeatures(j) = newFeatures(j) + i  // Clones
      newSizes(j) += 1
      if ( matrixIsCached ) {
        val newMatrix = matrix.clone  // Shallow copy
        newMatrix(i) = matrix(i).clone
        newMatrix(i)(j) = 1.0
        new FeatureAllocationWithAll(newMatrix, newFeatures, newSizes)
      } else {
        new FeatureAllocationWithArrayAndSizes(nItems, newFeatures, newSizes)
      }
    }
  }

  def add(i: Int, j: Int): FeatureAllocation = {
    if ( ( i < 0 ) || ( i >= nItems ) ) throw new IllegalArgumentException("Item index "+i+" is out of bounds [0,"+(nItems-1)+"].")
    if ( ( j < 0 ) || ( j >= nFeatures ) ) throw new IllegalArgumentException("Feature index "+j+" is out of bounds [0,"+(nFeatures-1)+"].")
    val newMatrix = matrix.map(_.clone)
    newMatrix(i)(j) = 1.0
    new FeatureAllocationWithMatrix(newMatrix)
  }

  def mutateAdd(i: Int, j: Int): Unit = {
    if ( ( i < 0 ) || ( i >= nItems ) ) throw new IllegalArgumentException("Item index "+i+" is out of bounds [0,"+(nItems-1)+"].")
    if ( ( j < 0 ) || ( j >= nFeatures ) ) throw new IllegalArgumentException("Feature index "+j+" is out of bounds [0,"+(nFeatures-1)+"].")
    if ( ! features(j)(i) ) {
      sizes(j) += 1
      featuresAsList(j) = i :: featuresAsList(j)  // Could be lazy, so must be before next line!
      features(j).add(i)                          // Mutates
      if ( matrixIsCached ) matrix(i)(j) = 1.0
    }
  }

  def mutateRemove(i: Int, j: Int): Unit = {
    if ( ( i < 0 ) || ( i >= nItems ) ) throw new IllegalArgumentException("Item index "+i+" is out of bounds [0,"+(nItems-1)+"].")
    if ( ( j < 0 ) || ( j >= nFeatures ) ) throw new IllegalArgumentException("Feature index "+j+" is out of bounds [0,"+(nFeatures-1)+"].")
    if ( features(j)(i) ) {
      sizes(j) -= 1
      featuresAsList(j) = featuresAsList(j).diff(List(i))  // Could be lazy, so must be before next line!
      features(j).remove(i)                                // Mutates
      if ( matrixIsCached ) matrix(i)(j) = 0.0
    }
  }

  def removeOld(i: Int, j: Int): FeatureAllocation = {
    if ( ( i < 0 ) || ( i >= nItems ) ) throw new IllegalArgumentException("Item index "+i+" is out of bounds [0,"+(nItems-1)+"].")
    if ( ( j < 0 ) || ( j >= nFeatures ) ) throw new IllegalArgumentException("Feature index "+j+" is out of bounds [0,"+(nFeatures-1)+"].")
    if ( ! features(j)(i) ) this
    else {
      val newFeatures = features.clone  // Shallow copy
      val newSizes = sizes.clone
      newFeatures(j) = newFeatures(j) - i  // Clones
      newSizes(j) -= 1
      if ( matrixIsCached ) {
        val newMatrix = matrix.clone  // Shallow copy
        newMatrix(i) = matrix(i).clone
        newMatrix(i)(j) = 0.0
        new FeatureAllocationWithAll(newMatrix, newFeatures, newSizes)
      } else {
        new FeatureAllocationWithArrayAndSizes(nItems, newFeatures, newSizes)
      }
    }
  }

  def remove(i: Int, j: Int): FeatureAllocation = {
    if ( ( i < 0 ) || ( i >= nItems ) ) throw new IllegalArgumentException("Item index "+i+" is out of bounds [0,"+(nItems-1)+"].")
    if ( ( j < 0 ) || ( j >= nFeatures ) ) throw new IllegalArgumentException("Feature index "+j+" is out of bounds [0,"+(nFeatures-1)+"].")
    val newMatrix = matrix.map(_.clone)
    newMatrix(i)(j) = 0.0
    new FeatureAllocationWithMatrix(newMatrix)
  }

  def removeRow(i: Int): FeatureAllocation = {
    val newMatrix = matrix.map(_.clone)
    newMatrix(i) = new Array[Double](nFeatures)
    val counts = (0 until nFeatures).map { j =>
      (0 until nItems).count{ i => newMatrix(i)(j) == 1.0 }
    }
    val thoseToKeep = counts.zipWithIndex.partition { case (count,index) => count != 0 }._1.map(_._2)
    val newNewMatrix = Array.tabulate(nItems, thoseToKeep.sum) { case (i,j) =>
      newMatrix(i)(thoseToKeep(j))
    }
    new FeatureAllocationWithMatrix(newNewMatrix)
  }

  def remove(i: Int): FeatureAllocation = remove(Array(i), false)

  def remove(i: Int, cloneAndKeepEmptyFeatures: Boolean): FeatureAllocation = remove(Array(i), cloneAndKeepEmptyFeatures)

  def remove(i: Iterable[Int]): FeatureAllocation = remove(i, false)

  def remove(i: Iterable[Int], cloneAndKeepEmptyFeatures: Boolean): FeatureAllocation = {
    val newFeatures = if ( cloneAndKeepEmptyFeatures ) features.map(_.clone) else features.clone
    val newSizes = sizes.clone
    var nToDelete = 0
    for ( ii <- i ) {
      if ( ( ii < 0 ) || ( ii >= nItems ) ) throw new IllegalArgumentException("Item index "+ii+" is out of bounds [0"+(nItems-1)+"].")
      var j = 0
      while (j < nFeatures) {
        if ( ( newSizes(j) > 0 ) && features(j)(ii) ) {
          if ( cloneAndKeepEmptyFeatures ) newFeatures(j).remove(ii) else newFeatures(j) = newFeatures(j) - ii
          newSizes(j) -= 1
          if (newSizes(j) == 0) nToDelete += 1
        }
        j += 1
      }
    }
    if ( ( !cloneAndKeepEmptyFeatures ) && ( nToDelete > 0 ) ) {
      val newFeatures2 = new Array[BitSet](nFeatures - nToDelete)
      val newSizes2 = new Array[Int](nFeatures - nToDelete)
      var j = 0
      var jj = 0
      while (jj < newSizes2.length) {
        while (newSizes(j) > 0) {
          newFeatures2(jj) = newFeatures(j)
          newSizes2(jj) = newSizes(j)
          j += 1
          jj += 1
        }
        j += 1
      }
      new FeatureAllocationWithArrayAndSizes(nItems, newFeatures2, newSizes2)
    } else {
      if ( matrixIsCached ) {
        val newMatrix = if ( cloneAndKeepEmptyFeatures ) matrix.map(_.clone) else matrix.clone
        i.foreach { newMatrix(_) = new Array[Double](nFeatures) }
        new FeatureAllocationWithAll(newMatrix, newFeatures, newSizes)
      } else {
        new FeatureAllocationWithArrayAndSizes(nItems, newFeatures, newSizes)
      }
    }
  }

  def matrixWithout(i: Int): Array[Array[Double]] = {
    if ( ( i < 0 ) || ( i >= nItems ) ) throw new IllegalArgumentException("Item index "+i+" is out of bounds [0,"+(nItems-1)+"].")
    val newMatrix = matrix.clone  // Shallow copy
    newMatrix(i) = new Array[Double](nFeatures)
    newMatrix
  }

  def pairwiseAllocationMatrix: Array[Array[Int]] = {
    val result = Array.ofDim[Int](nItems,nItems)
    for ( k <- 0 until nFeatures ) {
      var i = 0
      var hits = 0
      while ( hits < sizes(k) ) {
        if ( features(k).contains(i) ) {
          result(i)(i) += 1
          hits += 1
          var j = i + 1
          while ( j < nItems ) {
            if ( features(k).contains(j) ) {
              result(i)(j) += 1
              result(j)(i) += 1
            }
            j += 1
          }
        }
        i += 1
      }
    }
    result
  }

  def isSingleton(i: Int, j: Int): Boolean = {
    if ( ( i < 0 ) || ( i >= nItems ) ) throw new IllegalArgumentException("Item index "+i+" is out of bounds [0,"+(nItems-1)+"].")
    if ( ( j < 0 ) || ( j >= nFeatures ) ) throw new IllegalArgumentException("Feature index "+j+" is out of bounds [0,"+(nFeatures-1)+"].")
    ( sizes(j) == 1 ) && features(j)(i)
  }

  def partitionBySingletonsOf(i: Int): (FeatureAllocation, FeatureAllocation) = {
    if ( ( i < 0 ) || ( i >= nItems ) ) throw new IllegalArgumentException("Item index "+i+" is out of bounds [0,"+(nItems-1)+"].")
    var sum = 0
    val sel = Array.tabulate(nFeatures) { j =>
      val result = isSingleton(i,j)
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

  def enumerateFor(i: Int): Array[FeatureAllocation] = {   // Careful, this results in tons of shared mutable instances.
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

  def computeRegardingTiesSlow(fa: FeatureAllocation): Double = {
    fa.features.groupBy(identity).map(_._2.length).foldLeft(0.0)((s, x) => s + logFactorial(x))
  }

  def tiesMapper[B](initial: => B, aggregator: (B,(BitSet,Int),Int) => B): B = {
    val aa = features.zip(sizes).sortWith(lessThan)
    var j = 1
    var sum = initial
    if ( nFeatures == 0 ) sum
    else {
      var run = 1
      while ( j < nFeatures ) {
        if ( compare(aa(j-1),aa(j)) == 0 ) run += 1
        else {
          sum = aggregator(sum,aa(j-1),run)
          run = 1
        }
        j += 1
      }
      aggregator(sum,aa(j-1),run)
    }
  }

  def asMap: Map[(BitSet,Int),Int] = tiesMapper(Map[(BitSet,Int),Int](), (sum: Map[(BitSet,Int),Int], pair, count) => {
    sum + { pair -> count }
  })

  def computeRegardingTies: Double = tiesMapper(0.0, (sum: Double, pair, count) => {
    sum + logFactorial(count)
  })

  private def compare(x: (Iterable[Int],Int), y: (Iterable[Int], Int)): Int = {
    if ( x._2 < y._2 ) return -1
    else if ( x._2 > y._2 ) return 1
    else {
      val xi = x._1.iterator
      val yi = y._1.iterator
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

  private def lessThan(x: (Iterable[Int],Int), y: (Iterable[Int], Int)): Boolean = compare(x,y) < 0

  def convertToAlternativeImplementation: FeatureAllocationAlternative[Null] = {
    FeatureAllocationAlternative(nItems, features.map { f =>
      FeatureAlternative(f.toArray:_*)
    }:_*)
  }

}

sealed class FeatureAllocationNone private[lineargaussian](override val nItems: Int) extends FeatureAllocation {

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

sealed class FeatureAllocationEmpty private[lineargaussian] (override val nItems: Int, override val nFeatures: Int) extends FeatureAllocation {

  override val sizes = Array.ofDim[Int](nFeatures)
  override val features = Array.fill(nFeatures)(BitSet())
  override val featuresAsList = Array.fill(nFeatures)(List[Int]())
  override val matrix = Array.ofDim[Double](nItems,nFeatures)
  matrixIsCached = true

}

object FeatureAllocation {

  def apply(nItems: Int): FeatureAllocation = {
    if ( nItems < 0 ) throw new IllegalArgumentException("Number of items must be at least 0.")
    new FeatureAllocationNone(nItems)
  }

  def apply(nItems: Int, nFeatures: Int): FeatureAllocation = {
    if ( nItems < 0 ) throw new IllegalArgumentException("Number of items must be at least 0.")
    if ( nFeatures < 0 ) throw new IllegalArgumentException("Number of features must be at least 0.")
    new FeatureAllocationEmpty(nItems,nFeatures)
  }

  def apply(matrix: Array[Array[Int]]): FeatureAllocation = apply(matrix.map(_.map(_.toDouble)))

  def apply(matrix: Array[Array[Double]]): FeatureAllocation = {
    val rows = matrix.length
    if ( rows == 0 ) return new FeatureAllocationNone(rows)
    val cols = matrix(0).length
    if ( ! matrix.forall(_.length == cols) ) throw new IllegalArgumentException("Number of features is not consistent.")
    if (cols == 0) return new FeatureAllocationNone(rows)
    if ( ! matrix.forall(_.forall(x => ( x == 0.0 ) || ( x == 1.0 ))) ) throw new IllegalArgumentException("Elements should be either 0 or 1.")
    val allColumnsNonempty = (0 until cols).forall { j =>
      (0 until rows).exists{ i => matrix(i)(j) == 1.0 }
    }
    if ( ! allColumnsNonempty) throw new IllegalArgumentException("Empty features are not permitted.")
    new FeatureAllocationWithMatrix(matrix)
  }

  def apply(fa: FeatureAllocation): FeatureAllocation = {  // Shallow clones
    fa match {
      case e: FeatureAllocationNone => new FeatureAllocationNone(e.nItems)
      case e: FeatureAllocationWithMatrix => new FeatureAllocationWithMatrix(e.matrix.clone)
      case e: FeatureAllocationWithArray => new FeatureAllocationWithArray(e.nItems, e.features.clone)
      case e: FeatureAllocationWithArrayAndSizes => new FeatureAllocationWithArrayAndSizes(e.nItems, e.features.clone, e.sizes.clone)
      case e: FeatureAllocationWithMatrixAndArray => new FeatureAllocationWithMatrixAndArray(e.matrix.clone, e.features.clone)
      case e: FeatureAllocationWithAll => new FeatureAllocationWithAll(e.matrix.clone, e.features.clone, e.sizes.clone)
      case e: FeatureAllocationEmpty => throw new IllegalArgumentException("This subtype cannot be shallow copied.")
    }
  }

  def enumerate(nItems: Int, maxNFeatures: Int): Array[FeatureAllocation] = {
    import scala.collection.parallel.mutable.ParArray
    def engine(i: Int, fa: FeatureAllocation): ParArray[FeatureAllocation] = {
      var state = fa
      var bag = state.enumerateFor(i).par
      for ( k <- state.nFeatures until maxNFeatures ) {
        state = state.add(i)
        bag ++= state.enumerateFor(i)
      }
      if ( i < nItems-1 ) bag.flatMap(fa => engine(i+1, fa)) else bag
    }
    engine(0,apply(nItems)).toArray
  }

  def convertFromAlternativeImplementation(faa: FeatureAllocationAlternative[Null]): FeatureAllocation = {
    apply(faa.toMatrix)
  }

}

