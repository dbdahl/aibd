package org.ddahl.aibd

import org.apache.commons.math3.util.FastMath.log

object Utils {

  val harmonicNumber = new CachedSumFunctionNonnegativeInt(1.0/_, 0.0, 1.2)
  val logFactorial = new CachedSumFunctionNonnegativeInt(i => log(i), 0.0, 1.2)
  val logOnInt = new CachedFunctionNonnegativeInt(i => log(i), 1.2)

  sealed class CachedSumFunctionNonnegativeInt private[Utils] (f: Int => Double, initialValue: Double, growFactor: Double) {

    final def apply(i: Int): Double = {
      if ( i >= seq.size ) synchronized {
        val max = (growFactor*(i max seq.size)).ceil.toInt
        seq = seq.init ++ (seq.size until max).scanLeft(seq.last)( (sum,i) => sum + f(i) )
      }
      seq(i)
    }

    private final var seq: Array[Double] = Array(initialValue)

  }

  sealed class CachedFunctionNonnegativeInt private[Utils] (f: Int => Double, growFactor: Double) {

    final def apply(i: Int): Double = {
      if ( i >= seq.size ) synchronized {
        val max = (growFactor*(i max seq.size)).ceil.toInt
        seq = seq ++ (seq.size until max).map(f)
      }
      seq(i)
    }

    private final var seq: Array[Double] = Array(f(0))

  }

}

