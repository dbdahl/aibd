package org.ddahl.aibd.util

import Misc.repeat
import org.apache.commons.math3.random.RandomDataGenerator

class RandomSweep private(rdg: RandomDataGenerator) {

  private var blocks = Vector[(String, Double, TimeMonitor, () => Unit)]()
  private var cumweights = Vector[Double]()

  def add[B](label: String, weight: Double)(block: => B) = {
    val timer = TimeMonitor()
    val blockFunction = () => {
      timer {
        block
      }
      ()
    }
    blocks = blocks :+ ((label, weight, timer, blockFunction))
    cumweights = blocks.scanLeft(0.0)(_ + _._2)
  }

  def times = blocks.map(x => x._1 -> x._3.total).toMap

  def apply(nUpdates: Int = 1): Unit = {
    repeat(nUpdates) {
      val x = rdg.nextUniform(0.0, cumweights(cumweights.size - 1))
      var i = 1
      while (x > cumweights(i)) {
        i += 1
      }
      blocks(i - 1)._4()
    }
  }

  override def toString = {
    blocks.map(x => {
      x._1 + ": " + x._3
    }).mkString(", ")
  }

}

object RandomSweep {

  def apply(rdg: RandomDataGenerator) = new RandomSweep(rdg)

}
