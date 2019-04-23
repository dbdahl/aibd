package org.ddahl.aibd

import org.ddahl.aibd.parameter.{NullParameterDistribution, ParameterDistribution}
import org.apache.commons.math3.random.RandomDataGenerator
import org.apache.commons.math3.util.CombinatoricsUtils.{factorial, factorialLog}
import org.apache.commons.math3.util.FastMath.{exp, log}

class MarginalizedAttractionIndianBuffetDistribution[A] protected (val mass: Double, val similarity: Similarity, val parameterDistribution: ParameterDistribution[A], val enumeration: Boolean) extends FeatureAllocationDistribution[A] {

  val nItems = similarity.nItems

  private val nFactorial = if (enumeration) factorial(nItems).toInt else 0

  private val logNFactorial = if (enumeration) factorialLog(nItems) else Double.NegativeInfinity

  private val allAIBD = if (enumeration) Permutation.enumerate(nItems).map(AttractionIndianBuffetDistribution(mass, _, similarity, parameterDistribution)) else Vector[AttractionIndianBuffetDistribution[A]]()

  def dropParameter = new MarginalizedAttractionIndianBuffetDistribution[Null](mass, similarity, NullParameterDistribution, enumeration)

  def sample(rdg: RandomDataGenerator): FeatureAllocation[A] = {
    if (enumeration) {
      allAIBD(rdg.nextInt(0, nFactorial - 1)).sample(rdg)
    } else {
      AttractionIndianBuffetDistribution(mass, Permutation.random(nItems, rdg), similarity, parameterDistribution).sample(rdg)
    }
  }

  def logDensity(fa: FeatureAllocation[A], parallel: Boolean): Double = {
    if (enumeration) {
      if (parallel) {
        log(allAIBD.par.map(d => exp(d.logDensity(fa, parallel))).sum) - logNFactorial
      } else {
        log(allAIBD.map(d => exp(d.logDensity(fa, parallel))).sum) - logNFactorial
      }
    } else throw new IllegalStateException("This function is only implemented when enumeration is true.")
  }

}

object MarginalizedAttractionIndianBuffetDistribution {

  def apply(mass: Double, similarity: Similarity, enumeration: Boolean) = {
    new MarginalizedAttractionIndianBuffetDistribution(mass, similarity, NullParameterDistribution, enumeration)
  }

  def apply[A](mass: Double, similarity: Similarity, enumeration: Boolean, parameterDistribution: ParameterDistribution[A]) = {
    new MarginalizedAttractionIndianBuffetDistribution(mass, similarity, parameterDistribution, enumeration)
  }

}
