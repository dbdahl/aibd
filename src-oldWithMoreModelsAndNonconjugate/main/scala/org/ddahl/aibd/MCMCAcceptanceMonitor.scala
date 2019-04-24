package org.ddahl.aibd

class MCMCAcceptanceMonitor1 {

  private var num1 = 0L
  private var den1 = 0L

  def nAcceptances: Long = num1

  def nAcceptances1: Long = num1

  def nAttempts: Long = den1

  def nAttempts1: Long = den1

  def rate: Double = num1.toDouble / den1

  def rate1: Double = num1.toDouble / den1

  def apply[B](func: => (B, Int, Int)): B = {
    val x = func
    num1 += x._2
    den1 += x._3
    x._1
  }

}

object MCMCAcceptanceMonitor1 {

  def apply() = new MCMCAcceptanceMonitor1()

}

class MCMCAcceptanceMonitor2 {

  private var num1 = 0L
  private var den1 = 0L
  private var num2 = 0L
  private var den2 = 0L

  def nAcceptances: Long = num1

  def nAcceptances1: Long = num1

  def nAcceptances2: Long = num2

  def nAttempts: Long = den1

  def nAttempts1: Long = den1

  def nAttempts2: Long = den1

  def rate: Double = num1.toDouble / den1

  def rate1: Double = num1.toDouble / den1

  def rate2: Double = num2.toDouble / den2

  def apply[B](func: => (B, Int, Int, Int, Int)): B = {
    val x = func
    num1 += x._2
    den1 += x._3
    num2 += x._4
    den2 += x._5
    x._1
  }

}

object MCMCAcceptanceMonitor2 {

  def apply() = new MCMCAcceptanceMonitor2()

}
