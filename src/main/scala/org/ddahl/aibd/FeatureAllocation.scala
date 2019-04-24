package org.ddahl.aibd

import util.Functions.logFactorial

import org.ddahl.sdols.featureallocation.{FeatureAllocation => FeatureAllocationAlternative}
import org.ddahl.sdols.featureallocation.{Feature => FeatureAlternative}
import scala.collection.immutable.BitSet

sealed trait FeatureAllocation {

  val nItems: Int
  val nFeatures: Int
  val sizes: Vector[Int]
  val features: Vector[BitSet]
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

  def id = {
    asCountList.flatMap { case (bs, size, count) =>
      val str = bs.toBitMask.mkString(",")
      Array.fill(count)(str)
    }.mkString(":")
  }

  def check(): Unit = {
    assert(nItems == matrix.length)
    assert(nFeatures == features.length)
    assert(nFeatures == sizes.length)
    assert((nFeatures == 0) || (nFeatures == matrix(0).length))
    val that1 = new FeatureAllocationWithMatrix(matrix)
    val that2 = new FeatureAllocationWithFeatures(nItems, features)
    val that3 = new FeatureAllocationWithFeaturesAndSizes(nItems, features, sizes)
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

  protected def computeSizes: Vector[Int] = features.map(_.size)

  protected def computeFeatures: Vector[BitSet] = {
    val result = Array.fill(nFeatures)(BitSet())
    var j = 0
    while (j < nFeatures) {
      var bs = result(j)
      var i = 0
      while (i < nItems) {
        if (matrix(i)(j) != 0.0) bs += i
        i += 1
      }
      result(j) = bs
      j += 1
    }
    result.toVector
  }

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

  def add(fa: FeatureAllocation): FeatureAllocation = {
    if ( fa.nFeatures == 0 ) this
    else if ( this.nFeatures == 0 ) fa
    else new FeatureAllocationWithFeaturesAndSizes(nItems, features ++ fa.features, sizes ++ fa.sizes)
  }

  def add(i: Int): FeatureAllocation = {
    val newFeatures = features :+ BitSet(i)
    val newSizes = sizes :+ 1
    new FeatureAllocationWithFeaturesAndSizes(nItems, newFeatures, newSizes)
  }

  private def updated(i: Int, j: Int, add: Boolean): FeatureAllocation = {
    if ( ( i < 0 ) || ( i >= nItems ) ) throw new IllegalArgumentException("Item index "+i+" is out of bounds [0,"+(nItems-1)+"].")
    if ( ( j < 0 ) || ( j >= nFeatures ) ) throw new IllegalArgumentException("Feature index "+j+" is out of bounds [0,"+(nFeatures-1)+"].")
    if ( add == features(j)(i) ) this
    else {
      val newFeatures = features.updated(j, if ( add ) features(j) + i else features(j) - i)
      val newSizes = sizes.updated(j, sizes(j) + { if ( add ) 1 else -1 })
      if ( matrixIsCached ) {
        val newMatrix = matrix.clone    // Shallow copy
        newMatrix(i) = matrix(i).clone  // Deep copy
        newMatrix(i)(j) = if ( add ) 1.0 else 0.0
        new FeatureAllocationWithAll(newMatrix, newFeatures, newSizes)
      } else {
        new FeatureAllocationWithFeaturesAndSizes(nItems, newFeatures, newSizes)
      }
    }
  }

  def add(i: Int, j: Int): FeatureAllocation = updated(i,j,true)

  def remove(i: Int, j: Int): FeatureAllocation = updated(i,j,false)

  def remove(i: Iterable[Int]): FeatureAllocation = {
    val newFeatures = features.toArray
    val newSizes = sizes.toArray
    var nToDelete = 0
    for ( ii <- i ) {
      if ( ( ii < 0 ) || ( ii >= nItems ) ) throw new IllegalArgumentException("Item index "+ii+" is out of bounds [0"+(nItems-1)+"].")
      var j = 0
      while (j < nFeatures) {
        if ( ( newSizes(j) > 0 ) && features(j)(ii) ) {
          newFeatures(j) -= ii
          newSizes(j) -= 1
          if (newSizes(j) == 0) nToDelete += 1
        }
        j += 1
      }
    }
    if ( nToDelete > 0 ) {
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
      new FeatureAllocationWithFeaturesAndSizes(nItems, newFeatures2.toVector, newSizes2.toVector)
    } else {
      if ( matrixIsCached ) {
        val newMatrix = matrix.clone
        i.foreach { newMatrix(_) = new Array[Double](nFeatures) }
        new FeatureAllocationWithAll(newMatrix, newFeatures.toVector, newSizes.toVector)
      } else {
        new FeatureAllocationWithFeaturesAndSizes(nItems, newFeatures.toVector, newSizes.toVector)
      }
    }
  }

  def remove(i: Int): FeatureAllocation = remove(Array(i))

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
        leftArray(jl) = features(j)
        leftSizes(jl) = sizes(j)
        jl += 1
      } else {
        rightArray(jr) = features(j)
        rightSizes(jr) = sizes(j)
        jr += 1
      }
      j += 1
    }
    (new FeatureAllocationWithFeaturesAndSizes(nItems,leftArray.toVector, leftSizes.toVector),
     new FeatureAllocationWithFeaturesAndSizes(nItems,rightArray.toVector,rightSizes.toVector))
  }

  def enumerateFor(i: Int): Vector[FeatureAllocation] = {
    val (singletons, existing) = partitionBySingletonsOf(i)
    val a = existing.remove(i).features.groupBy(identity).mapValues(_.size)
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
    engine(b, singletons.features.toList)
    if ( collector.length == 0 ) Vector[FeatureAllocation]()
    else {
      val zeroedOutMatrix = new FeatureAllocationWithFeatures(nItems, collector.head.toVector).matrixWithout(i)
      collector.toVector.map { x =>
        val arr = x.toArray
        val mat = zeroedOutMatrix.clone // Shallow copy
        mat(i) = arr.map { f => if (f(i)) 1.0 else 0.0 }
        new FeatureAllocationWithMatrixAndFeatures(mat, arr.toVector)
      }
    }
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

  def asCountMap: Map[(BitSet,Int),Int] = tiesMapper(Map[(BitSet,Int),Int](), (sum: Map[(BitSet,Int),Int], pair, count) => {
    sum + { pair -> count }
  })

  def asCountList: List[(BitSet,Int,Int)] = tiesMapper(List[(BitSet,Int,Int)](), (sum: List[(BitSet,Int,Int)], pair, count) => {
    (pair._1, pair._2, count) :: sum
  })

  def asLists: Vector[List[Int]] = features.map(_.toList)

  def asListsWithout(except: Iterable[Int]): Vector[List[Int]] = {
    features.map { f => ( f -- except ).toList }
  }

  def computeRegardingTies: Double = tiesMapper(0.0, (sum: Double, pair, count) => {
    sum + logFactorial(count)
  })

  def computeRegardingTiesSlow(fa: FeatureAllocation): Double = {
    fa.features.groupBy(identity).map(_._2.length).foldLeft(0.0)((s, x) => s + logFactorial(x))
  }

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

sealed class FeatureAllocationNone private[aibd] (override val nItems: Int) extends FeatureAllocation {

  override val nFeatures = 0
  override val sizes = Vector[Int]()
  override val features = Vector[BitSet]()
  override val matrix = Array.ofDim[Double](nItems,0)
  matrixIsCached = true

}

sealed class FeatureAllocationWithMatrix private[aibd] (override val matrix: Array[Array[Double]]) extends FeatureAllocation {

  override val nItems = matrix.length
  override val nFeatures = matrix(0).length
  override lazy val sizes = computeSizes
  override lazy val features = computeFeatures
  matrixIsCached = true

}

sealed class FeatureAllocationWithFeatures private[aibd] (override val nItems: Int, override val features: Vector[BitSet]) extends FeatureAllocation {

  override val nFeatures = features.length
  override lazy val sizes = computeSizes
  override lazy val matrix = computeMatrix

}

sealed class FeatureAllocationWithFeaturesAndSizes private[aibd] (override val nItems: Int, override val features: Vector[BitSet], override val sizes: Vector[Int]) extends FeatureAllocation {

  override val nFeatures = features.length
  override lazy val matrix = computeMatrix

}

sealed class FeatureAllocationWithMatrixAndFeatures private[aibd] (override val matrix: Array[Array[Double]], override val features: Vector[BitSet]) extends FeatureAllocation {

  override val nItems = matrix.length
  override val nFeatures = features.length
  override lazy val sizes = computeSizes
  matrixIsCached = true

}

sealed class FeatureAllocationWithAll private[aibd] (override val matrix: Array[Array[Double]], override val features: Vector[BitSet], override val sizes: Vector[Int]) extends FeatureAllocation {

  override val nItems = matrix.length
  override val nFeatures = features.length
  matrixIsCached = true

}

sealed class FeatureAllocationEmpty private[aibd] (override val nItems: Int, override val nFeatures: Int) extends FeatureAllocation {

  override val sizes = Vector.fill(nFeatures) { 0 }
  override val features = Vector.fill(nFeatures)(BitSet())
  override val matrix = Array.ofDim[Double](nItems,nFeatures)
  matrixIsCached = true

}

object FeatureAllocation {

  def empty(nItems: Int): FeatureAllocation = {
    if ( nItems < 0 ) throw new IllegalArgumentException("Number of items must be at least 0.")
    new FeatureAllocationNone(nItems)
  }

  def fromBitSets(nItems: Int, features: Iterable[BitSet]): FeatureAllocation = {
    if ( nItems < 0 ) throw new IllegalArgumentException("Number of items must be at least 0.")
    new FeatureAllocationWithFeatures(nItems, features.toVector)
  }

  def fromLists(nItems: Int, lists: Iterable[List[Int]]): FeatureAllocation = {
    if ( nItems < 0 ) throw new IllegalArgumentException("Number of items must be at least 0.")
    new FeatureAllocationWithFeatures(nItems, lists.map { _.foldLeft(BitSet()) { _ + _ } }.toVector)
  }

  def fromMatrix(matrix: Array[Array[Double]]): FeatureAllocation = {
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

  def fromMatrix(matrix: Array[Array[Int]]): FeatureAllocation = fromMatrix(matrix.map(_.map(_.toDouble)))

  def fromCountMap(nItems: Int, map: Map[(BitSet,Int),Int]): FeatureAllocation = {
    val pair = map.map { case (pair, count) =>
      (Seq.fill(count)(pair._1), Seq.fill(count)(pair._2))
    }.foldLeft((Vector[BitSet](),Vector[Int]())) { case (sum, x) =>
      (sum._1 ++ x._1, sum._2 ++ x._2)
    }
    new FeatureAllocationWithFeaturesAndSizes(nItems, pair._1, pair._2)
  }

  def enumerate(nItems: Int, maxNFeatures: Int): Vector[FeatureAllocation] = {
    def engine(i: Int, fa: FeatureAllocation): Vector[FeatureAllocation] = {
      var state = fa
      var bag = state.enumerateFor(i)
      for ( k <- state.nFeatures until maxNFeatures ) {
        state = state.add(i)
        bag ++= state.enumerateFor(i)
      }
      if ( i < nItems-1 ) bag.flatMap(fa => engine(i+1, fa)) else bag
    }
    engine(0,empty(nItems))
  }

  def convertFromAlternativeImplementation(faa: FeatureAllocationAlternative[Null]): FeatureAllocation = {
    fromMatrix(faa.toMatrix)
  }

}

