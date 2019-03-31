package org.ddahl.aibd.model.lineargaussian

import org.ddahl.aibd.{AttractionIndianBuffetDistribution => AttractionIndianBuffetDistributionAlternative}
import org.ddahl.aibd.{Permutation, Similarity}
import org.ddahl.aibd.Utils.logOnInt
import org.apache.commons.math3.random.RandomDataGenerator
import org.apache.commons.math3.util.FastMath.log

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

