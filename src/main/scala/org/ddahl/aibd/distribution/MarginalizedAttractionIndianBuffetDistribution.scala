package org.ddahl.aibd.distribution

import org.ddahl.aibd._
import org.apache.commons.math3.random.RandomDataGenerator
import org.apache.commons.math3.util.FastMath.{exp, log}
import scala.collection.parallel.immutable.ParVector

class MarginalizedAttractionIndianBuffetDistribution private (val mass: Double, val similarity: Similarity) extends FeatureAllocationDistribution {

  implicit val ordering = CrossCompatibility.doubleOrdering

  val nItems = similarity.nItems

  def logProbability(i: Int, fa: FeatureAllocation): Double = logProbability(fa)

  def logProbability(fa: FeatureAllocation): Double = {
    val components = ParVector(Permutation.enumerate(nItems):_*).map(p => AttractionIndianBuffetDistribution(mass, p, similarity))
    val nPermutations = components.length
    val logProbs = components.map(_.logProbability(fa))
    val max = logProbs.max
    -log(nPermutations) + max + log(logProbs.map(lp => exp(lp-max)).sum)   // More numerically stable.
  }

  def sample(rdg: RandomDataGenerator): FeatureAllocation = {
    val component = AttractionIndianBuffetDistribution(mass, Permutation.random(nItems, rdg), similarity)
    component.sample(rdg)
  }

}

object MarginalizedAttractionIndianBuffetDistribution {

  def apply(mass: Double, similarity: Similarity): MarginalizedAttractionIndianBuffetDistribution = {
    if ( mass <= 0.0 ) throw new IllegalArgumentException("'mass' must be positive.")
    new MarginalizedAttractionIndianBuffetDistribution(mass, similarity)
  }

}

