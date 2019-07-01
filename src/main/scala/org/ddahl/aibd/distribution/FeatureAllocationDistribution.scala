package org.ddahl.aibd.distribution

import org.ddahl.aibd._
import org.ddahl.commonsmath.Implicits.RandomDataGeneratorImprovements
import org.apache.commons.math3.random.RandomDataGenerator

import scala.collection.parallel.immutable.ParVector
import scala.reflect.ClassTag

trait FeatureAllocationDistribution {

  val nItems: Int

  def logProbability(i: Int, fa: FeatureAllocation): Double

  def logProbability(i: Int, fa: Array[Array[Double]]): Double = logProbability(i, FeatureAllocation.fromMatrix(fa))

  def logProbability(i: Int, fa: Array[FeatureAllocation]): Array[Double] = fa.map(logProbability(i,_))

  def logProbability(i: Int, fa: Array[Array[Array[Double]]]): Array[Double] = fa.map(logProbability(i,_))

  def logProbability(fa: FeatureAllocation): Double

  def logProbability(fa: Array[Array[Double]]): Double = logProbability(FeatureAllocation.fromMatrix(fa))

  def logProbability(fa: Array[FeatureAllocation]): Array[Double] = fa.map(logProbability)

  def logProbability(fa: Array[Array[Array[Double]]]): Array[Double] = fa.map(logProbability)

  def sample(rdg: RandomDataGenerator): FeatureAllocation

  def sample[B: ClassTag](rdg: RandomDataGenerator, nSamples: Int, nCores: Int = 0, transop: FeatureAllocation => B = (x: FeatureAllocation) => x): Array[B] = {
    require(nSamples > 0, "nSamples must be at least 1.")
    require(nCores >= 0, "nCores cannot be negative.")
    val nCores2 = if ( nCores <= 0 ) Runtime.getRuntime.availableProcessors else nCores
    if (nCores2 > 1 ) {
      val nSamplesPerCore = (nSamples - 1) / nCores2 + 1
      rdg.nextRandomDataGenerators(nCores).flatMap(rdg => {
        Array.fill(nSamplesPerCore)(transop(sample(rdg)))
      }).toArray
    } else {
      Array.fill(nSamples)(transop(sample(rdg)))
    }
  }

  def expectedPairwiseAllocationMatrix(maxNFeatures: Int): Array[Array[Double]] = {
    import org.ddahl.matrix._
    import org.apache.commons.math3.util.FastMath.exp
    val allFeatures = ParVector(FeatureAllocation.enumerate(nItems, maxNFeatures):_*)
    val expectation = allFeatures.aggregate(matrixOfDim(nItems, nItems))( (sum,fa) =>
      sum add ( exp(logProbability(fa)) * wrap(fa.pairwiseAllocationMatrix) ), _ add _)
    expectation.getData()
  }

  def sum(maxNFeatures: Int): Double = {
    import org.apache.commons.math3.util.FastMath.exp
    FeatureAllocation.enumerate(nItems, maxNFeatures).foldLeft(0.0) { (sum,fa) =>
      sum + exp(logProbability(fa))
    }
  }

}


