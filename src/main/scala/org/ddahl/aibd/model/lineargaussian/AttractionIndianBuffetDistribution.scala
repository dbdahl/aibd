package org.ddahl.aibd.model.lineargaussian

import org.ddahl.aibd.{Permutation, Similarity}
import org.ddahl.aibd.Utils.logOnInt
import org.apache.commons.math3.random.RandomDataGenerator
import org.apache.commons.math3.util.FastMath.log

class AttractionIndianBuffetDistribution private (val mass: Double, val permutation: Permutation, val similarity: Similarity) extends FeatureAllocationDistribution with HasMass[AttractionIndianBuffetDistribution] {

  val nItems = similarity.nItems
  val logMass = log(mass)

  def updateMass(mass: Double): AttractionIndianBuffetDistribution = {
    if ( mass <= 0.0 ) throw new IllegalArgumentException("'mass' must be positive.")
    new AttractionIndianBuffetDistribution(mass, permutation, similarity)
  }

  def updatePermutation(permutation: Permutation): AttractionIndianBuffetDistribution = {
    if ( permutation.nItems != similarity.nItems ) throw new IllegalArgumentException("Inconsistent number of items.")
    new AttractionIndianBuffetDistribution(mass, permutation, similarity)
  }

  def updateSimilarity(similarity: Similarity): AttractionIndianBuffetDistribution = {
    if ( permutation.nItems != similarity.nItems ) throw new IllegalArgumentException("Inconsistent number of items.")
    new AttractionIndianBuffetDistribution(mass, permutation, similarity)
  }

  def logProbability(i: Int, fa: FeatureAllocation): Double = {
    var index = permutation.inverse(i)
    val state = fa.featuresAsListWithout(permutation.drop(index))
    var sum = 0.0
    while ( index < fa.nItems ) {
      val ii = permutation(index)
      val divisor = (0 until index).foldLeft(0.0) { (s,indexPrime) => s + similarity(ii,permutation(indexPrime)) }
      var newFeatureCount = 0
      var j = 0
      while ( j < fa.nFeatures ) {
        if ( state(j).isEmpty ) {
          if ( fa.features(j)(ii) ) {
            state(j) = ii :: state(j)
            newFeatureCount += 1
          }
        } else {
          val p = index / (index + 1.0) * state(j).foldLeft(0.0) { (s, iPrime) => s + similarity(ii,iPrime) } / divisor
          if ( fa.features(j)(ii) ) {
            state(j) = ii :: state(j)
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

  def logProbability(fa: FeatureAllocation): Double = logProbability(permutation(0), fa)

  def sample(rdg: RandomDataGenerator): FeatureAllocation = {
    val nNewFeaturesPerItems = Array.tabulate(nItems) { i => rdg.nextPoisson(mass / (i+1)).toInt }
    val nNewFeaturesCumulant = nNewFeaturesPerItems.scan(0)(_+_)
    var fa = FeatureAllocation(nItems)
    var index = 0
    while (index < nItems) {
      val ii = permutation(index)
      val divisor = (0 until index).foldLeft(0.0) { (s,indexPrime) => s + similarity(ii,permutation(indexPrime)) }
      var j = 0
      while (j < nNewFeaturesCumulant(index)) {
        val p = index / (index + 1.0) * fa.featuresAsList(j).foldLeft(0.0) { (s, iPrime) => s + similarity(ii,iPrime) } / divisor
        if ( rdg.nextUniform(0.0,1.0) <= p ) fa = fa.add(ii,j)
        j += 1
      }
      while (j < nNewFeaturesCumulant(index+1) ) {
        fa = fa.add(ii)
        j += 1
      }
      index += 1
    }
    fa
  }

}

object AttractionIndianBuffetDistribution {

  def apply(mass: Double, permutation: Permutation, similarity: Similarity): AttractionIndianBuffetDistribution = {
    if ( mass <= 0.0 ) throw new IllegalArgumentException("'mass' must be positive.")
    if ( permutation.nItems != similarity.nItems ) throw new IllegalArgumentException("Inconsistent number of items.")
    new AttractionIndianBuffetDistribution(mass, permutation, similarity)
  }

}

