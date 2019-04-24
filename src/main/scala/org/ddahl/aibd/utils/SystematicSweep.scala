package org.ddahl.aibd.util

import Misc.repeat

class SystematicSweep private () {

  private var blocks = Vector[(String, TimeMonitor, () => Unit)]()

  def add[B](label: String, reps: Int)(block: => B) = {
    val timer = TimeMonitor()
    val blockFunction = () => {
      timer {
        repeat(reps) {
          block
        }
      }
      ()
    }
    blocks = blocks :+ ((label, timer, blockFunction))
  }

  def times = blocks.map(x => x._1 -> x._2.total).toMap

  def apply(nUpdates: Int = 1): Unit = {
    repeat(nUpdates) {
      var i = 0
      while (i < blocks.size) {
        blocks(i)._3()
        i += 1
      }
    }
  }

  override def toString = {
    blocks.map(x => {
      x._1 + ": " + x._2
    }).mkString(", ")
  }

}

object SystematicSweep {

  def apply() = new SystematicSweep()

}
