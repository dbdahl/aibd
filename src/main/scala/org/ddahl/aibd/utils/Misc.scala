package org.ddahl.aibd.util

object Misc {

  @scala.annotation.tailrec
  def repeat(n: Long)(f: => Unit): Unit = {
    if (n > 0) {
      f
      repeat(n - 1)(f)
    }
  }

  def time[R](block: => R): R = {
    val t0 = System.nanoTime()
    val result = block // call-by-name
    val t1 = System.nanoTime()
    println("Elapsed time: " + (t1 - t0) + "ns")
    result
  }

}
