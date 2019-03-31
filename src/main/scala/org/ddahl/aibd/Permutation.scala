package org.ddahl.aibd

import org.apache.commons.math3.random.RandomDataGenerator

class Permutation private (private val x: Array[Int], val nPerShuffle: Int) extends Iterable[Int] {

  def iterator = x.iterator

  override def isEmpty = x.isEmpty

  override def size = x.size

  val nItems: Int = x.size

  def apply(index: Int): Int = x(index)

  override def drop(n: Int) = x.drop(n)

  lazy val inverse: Permutation = {
    val y = new Array[Int](nItems)
    for (i <- x.indices) y(x(i)) = i
    new Permutation(y, nPerShuffle)
  }

  def shuffle(rdg: RandomDataGenerator): Permutation = {
    val xx = x.clone
    val indices = rdg.nextPermutation(nItems, nPerShuffle)
    val sorted = indices.sorted
    for (i <- 0 until nPerShuffle) {
      xx(sorted(i)) = x(indices(i))
    }
    new Permutation(xx, nPerShuffle)
  }

  def nPerShuffle(n: Int): Permutation = new Permutation(x, n)

  def toArray = x.clone

}

object Permutation {

  def apply(x: Int*): Permutation = {
    // Check that it's a valid permutation.
    val visited = new Array[Boolean](x.length)
    for (i <- x.indices) {
      val j = x(i)
      if (j < 0 || j >= x.length || visited(j)) throw new IllegalArgumentException("Not a permutation.")
      visited(j) = true
    }
    new Permutation(Array(x: _*), x.length)
  }

  def apply(x: Array[Int]): Permutation = apply(x: _*)

  def natural(nItems: Int): Permutation = new Permutation(Array.tabulate(nItems)(identity), nItems)

  def random(nItems: Int, rdg: RandomDataGenerator): Permutation = new Permutation(rdg.nextPermutation(nItems, nItems), nItems)

  def enumerate(nItems: Int): Vector[Permutation] = {
    Array.tabulate(nItems)(identity).permutations.map(new Permutation(_,nItems)).toVector
  }

}
