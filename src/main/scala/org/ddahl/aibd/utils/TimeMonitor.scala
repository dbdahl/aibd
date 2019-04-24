package org.ddahl.aibd.util

class TimeMonitor private(private var sum: Long, private var invocations: Long) {

  def total = sum

  def nInvocations = invocations

  def apply[B](block: => B): B = {
    invocations += 1
    val t0 = System.nanoTime()
    val result = block // call-by-name
    sum += System.nanoTime() - t0
    result
  }

  def /(d: TimeMonitor) = sum.toDouble / d.total

  def +(d: TimeMonitor) = new TimeMonitor(sum + d.sum, invocations + d.invocations)

  def rate = sum.toDouble / invocations

  def reset(): Unit = {
    sum = 0L
    invocations = 0L
  }

  override def toString = "%d invocations in %4.4f seconds".format(invocations, sum / 1e9)

}

object TimeMonitor {

  def apply() = new TimeMonitor(0L, 0L)

}
