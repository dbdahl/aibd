package org.ddahl.aibd

import org.ddahl.aibd.parameter.ParameterDistribution
import org.ddahl.commonsmath._
import org.apache.commons.math3.linear.RealMatrix
import org.apache.commons.math3.random.RandomDataGenerator
import org.apache.commons.math3.util.FastMath.exp

abstract class FeatureAllocationDistribution[A] {

  val mass: Double

  require(mass>0,"Mass must be greater than 0.")

  val nItems: Int

  val permutation = Permutation.natural(nItems)

  val parameterDistribution: ParameterDistribution[A]

  def dropParameter: FeatureAllocationDistribution[Null]

  def sample(rdg: RandomDataGenerator): FeatureAllocation[A]

  def sample[B](rdg: RandomDataGenerator, nSamples: Int, nCores: Int = Runtime.getRuntime.availableProcessors, transop: FeatureAllocation[A] => B = (x: FeatureAllocation[A]) => x): Vector[B] = {
    require(nSamples > 0, "nSamples must be at least 1.")
    require(nCores >= 0, "nCores cannot be negative.")
    val nCores2 = if ( nCores == 0 ) Runtime.getRuntime.availableProcessors else nCores
    if (nCores2 > 1 ) {
      val nSamplesPerCore = (nSamples - 1) / nCores2 + 1
      rdg.nextRandomDataGenerators(nCores).flatMap(rdg => {
        Vector.fill(nSamplesPerCore)(transop(sample(rdg)))
      }).toVector
    } else {
      Vector.fill(nSamples)(transop(sample(rdg)))
    }
  }

  def logDensity(fa: FeatureAllocation[A]): Double = logDensity(fa, false)

  def logDensity(fa: FeatureAllocation[A], parallel: Boolean): Double

  def logDensityWithParameters(fa: FeatureAllocation[A], parallel: Boolean): Double = {
    logDensity(fa,parallel) + fa.foldLeft(0.0)( (sum,f) => sum + parameterDistribution.logDensity(f.parameter) )
  }

  def logDensityStartingFromIndex(index: Int, fa: FeatureAllocation[A], parallel: Boolean): Double = logDensity(fa, parallel)

  def update(nScans: Int, fa: FeatureAllocation[A], rdg: RandomDataGenerator, parallel: Boolean): (FeatureAllocationDistribution[A], Int, Int) = (this,0,0)

  def integrateEnumeration[B](f: (FeatureAllocation[A], Double) => B, combop: (B, B) => B, maxNFeatures: Int, parallel: Boolean = false): B = {
    var sum: Option[B] = None
    FeatureAllocation.foreach(nItems, parameterDistribution.discreteSupport.get, maxNFeatures) {fa =>
      val v = f(fa, exp(logDensity(fa, parallel)))
      if ( sum.isEmpty ) sum = Some(v)
      else sum = Some(combop(sum.get, v))
    }
    sum.get
  }

  def integrateEnumerationMatrix(f: FeatureAllocation[A] => RealMatrix, maxNFeatures: Int, parallel: Boolean = false): RealMatrix = {
    integrateEnumeration((x: FeatureAllocation[A], p: Double) => f(x) :* p, (x: RealMatrix, y: RealMatrix) => x+y, maxNFeatures, parallel)
  }

  def integrateMonteCarlo[B,C](f: (FeatureAllocation[A]) => B, combop: (B, B) => B, scaleop: (B, Int) => C, rdg: RandomDataGenerator, nSamples: Int, nCores: Int = Runtime.getRuntime.availableProcessors): C = {
    require(nSamples > 0, "nSamples must be at least 1.")
    require(nCores >= 0, "nCores cannot be negative.")
    val nCores2 = if ( nCores == 0 ) Runtime.getRuntime.availableProcessors else nCores
    if (nCores2 > 1) {
      val rdgs = rdg.nextRandomDataGenerators(nCores2)
      val nSamplesPerCore = (nSamples - 1) / nCores2 + 1
      val sum = rdgs.map(r => {
        var counter = 0
        var sum = f(sample(rdg))
        while ( counter < nSamplesPerCore-1 ) {
          sum = combop(sum, f(sample(rdg)))
          counter += 1
        }
        sum 
      }).reduce(combop)
      scaleop(sum, nCores2*nSamplesPerCore)
    } else {
      var counter = 0
      var sum = f(sample(rdg))
      while ( counter < nSamples-1 ) {
        sum = combop(sum, f(sample(rdg)))
        counter += 1
      }
      scaleop(sum, nSamples)
    }
  }

  def integrateMonteCarloMatrix(f: FeatureAllocation[A] => RealMatrix, rdg: RandomDataGenerator, nSamples: Int, nCores: Int = Runtime.getRuntime.availableProcessors): RealMatrix = {
    integrateMonteCarlo(f, (x: RealMatrix, y: RealMatrix) => x+y, (t: RealMatrix, n: Int) => t :/ n, rdg, nSamples, nCores)
  }

  def integrateMonteCarloArray(f: FeatureAllocation[A] => Array[Double], rdg: RandomDataGenerator, nSamples: Int, nCores: Int = Runtime.getRuntime.availableProcessors): Array[Double] = {
    integrateMonteCarlo(f, (x: Array[Double], y: Array[Double]) => Array.tabulate(x.length) { i => x(i)+y(i) }, (t: Array[Double], n: Int) => t.map(_/n), rdg, nSamples, nCores)
  }

  def integrateMonteCarloDouble(f: FeatureAllocation[A] => Double, rdg: RandomDataGenerator, nSamples: Int, nCores: Int = Runtime.getRuntime.availableProcessors): Double = {
    integrateMonteCarlo(f, (x: Double, y: Double) => x+y, (t: Double, n: Int) => t/n, rdg, nSamples, nCores)
  }

  def toEmpericalDistributionViaDirect(rdg: RandomDataGenerator, nSamples: Int, nCores: Int = Runtime.getRuntime.availableProcessors): EmpericalFeatureAllocationDistribution = {
    val dist = this.dropParameter
    val nSamplesPerCore = (nSamples - 1) / nCores + 1
    val efad = rdg.nextRandomDataGenerators(nCores).map(rdg => {
      var efad = EmpericalFeatureAllocationDistribution()
      repeat(nSamplesPerCore) {
        efad = efad.tally(dist.sample(rdg))
      }
      efad
    }).fold(EmpericalFeatureAllocationDistribution())((x, y) => {
      x ++ y
    })
    efad
  }

}
