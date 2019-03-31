package org.ddahl.aibd.model.lineargaussian

import org.ddahl.aibd.{AttractionIndianBuffetDistribution => AttractionIndianBuffetDistributionAlternative, IndianBuffetProcess => IndianBuffetProcessAlternative}
import org.apache.commons.math3.random.RandomDataGenerator
import org.ddahl.aibd.Utils._
import org.apache.commons.math3.util.FastMath.log
import org.ddahl.aibd.{Permutation, Similarity}
import scala.reflect.ClassTag
import org.ddahl.commonsmath._

trait FeatureAllocationDistribution {

  val nItems: Int

  def logProbability(i: Int, fa: FeatureAllocation): Double

  def logProbability(i: Int, fa: Array[Array[Double]]): Double = logProbability(i, FeatureAllocation(fa))

  def logProbability(i: Int, fa: Array[FeatureAllocation]): Array[Double] = fa.map(logProbability(i,_))

  def logProbability(i: Int, fa: Array[Array[Array[Double]]]): Array[Double] = fa.map(logProbability(i,_))

  def logProbability(fa: FeatureAllocation): Double = logProbability(0,fa)

  def logProbability(fa: Array[Array[Double]]): Double = logProbability(0,fa)

  def logProbability(fa: Array[FeatureAllocation]): Array[Double] = fa.map(logProbability)

  def logProbability(fa: Array[Array[Array[Double]]]): Array[Double] = fa.map(logProbability)

  def sample(rdg: RandomDataGenerator): FeatureAllocation

  def sample[B: ClassTag](rdg: RandomDataGenerator, nSamples: Int, nCores: Int = Runtime.getRuntime.availableProcessors, transop: FeatureAllocation => B = (x: FeatureAllocation) => x): Array[B] = {
    require(nSamples > 0, "nSamples must be at least 1.")
    require(nCores >= 0, "nCores cannot be negative.")
    val nCores2 = if ( nCores == 0 ) Runtime.getRuntime.availableProcessors else nCores
    if (nCores2 > 1 ) {
      val nSamplesPerCore = (nSamples - 1) / nCores2 + 1
      rdg.nextRandomDataGenerators(nCores).flatMap(rdg => {
        Array.fill(nSamplesPerCore)(transop(sample(rdg)))
      }).toArray
    } else {
      Array.fill(nSamples)(transop(sample(rdg)))
    }
  }

}

class AttractionIndianBuffetDistribution private (val mass: Double, val permutation: Permutation, val similarity: Similarity) extends FeatureAllocationDistribution {

  val nItems = permutation.nItems
  val logMass = log(mass)

  def logProbability(i: Int, fa: FeatureAllocation): Double = {
    var index = permutation.inverse(i)
    val state = fa.remove(permutation.drop(index), true)
    var sum = 0.0
    while ( index < fa.nItems ) {
      val ii = permutation(index)
      val divisor = (0 until index).foldLeft(0.0) { (s,iPrime) => s + similarity(ii,permutation(iPrime)) }
      var newFeatureCount = 0
      var j = 0
      while ( j < fa.nFeatures ) {
        if ( state.isEmpty(j) ) {
          if ( fa.features(j)(ii) ) {
            state.mutateAdd(ii,j)
            newFeatureCount += 1
          }
        } else {
          val p = index.toDouble / (index + 1) * state.featuresAsList(j).foldLeft(0.0) { (s, iPrime) => s + similarity(ii,iPrime) } / divisor
          if ( fa.features(j)(ii) ) {
            state.mutateAdd(ii,j)
            sum += log(p)
          } else sum += log(1-p)
        }
        j += 1
      }
      index += 1
      sum += newFeatureCount * ( logMass - logOnInt(index) ) - mass/index
    }
    sum -= fa.computeRegardingTies
    sum
  }

  def sample(rdg: RandomDataGenerator): FeatureAllocation = {
    val aibd = AttractionIndianBuffetDistributionAlternative(mass, permutation, similarity)
    FeatureAllocation.convertFromAlternativeImplementation(aibd.sample(rdg))
  }

}

object AttractionIndianBuffetDistribution {

  def apply(mass: Double, permutation: Permutation, similarity: Similarity): AttractionIndianBuffetDistribution = {
    if ( mass <= 0.0 ) throw new IllegalArgumentException("'mass' must be positive.")
    new AttractionIndianBuffetDistribution(mass, permutation, similarity)
  }

}

class IndianBuffetProcess private (val mass: Double, val nItems: Int) extends FeatureAllocationDistribution {

  val logMass = log(mass)

  def logProbability(i: Int, fa: FeatureAllocation): Double = {
    val const1 = -mass * harmonicNumber(fa.nItems)
    if ( fa.nFeatures == 0 ) return const1
    val const2 = logMass - logFactorial(fa.nItems)
    var sum = const1 + fa.nFeatures * const2
    sum -= fa.computeRegardingTies
    sum += fa.sizes.foldLeft(0.0) { (s, mk) =>
      s + logFactorial(fa.nItems - mk) + logFactorial(mk - 1)
    }
    sum
  }

  def sample(rdg: RandomDataGenerator): FeatureAllocation = {
    val ibp = IndianBuffetProcessAlternative(mass, nItems)
    FeatureAllocation.convertFromAlternativeImplementation(ibp.sample(rdg))
  }

}

object IndianBuffetProcess {

  def apply(mass: Double, nItems: Int): IndianBuffetProcess = {
    if ( mass <= 0.0 ) throw new IllegalArgumentException("'mass' must be positive.")
    new IndianBuffetProcess(mass, nItems)
  }

}

