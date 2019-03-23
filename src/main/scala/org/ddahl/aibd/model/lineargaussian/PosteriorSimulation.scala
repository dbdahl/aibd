package org.ddahl.aibd.model.lineargaussian

import org.apache.commons.math3.random.RandomDataGenerator
import org.apache.commons.math3.util.FastMath.log
import org.ddahl.commonsmath.RandomDataGeneratorImprovements
import org.ddahl.aibd.TimeMonitor

object PosteriorSimulation {

  def updateFeatureAllocationViaNeighborhoods(featureAllocation: FeatureAllocation, mass: Double, lglfm: LinearGaussianLatentFeatureModel, nSamples: Int, thin: Int, width: Int, rdg: RandomDataGenerator, newFeaturesTruncationDivisor: Double = 1000): Array[FeatureAllocation] = {
    val nItems = lglfm.N
    var state = featureAllocation
    val tmAll = TimeMonitor()
    val tmPrior = TimeMonitor()
    val tmLikelihood = TimeMonitor()
    val tmEnumeration = TimeMonitor()
    val nIterations = thin*nSamples
    val rate = if ( width <= 0 ) 1
    else {
      print("[" + (" " * width) + "]" + ("\b" * (width + 1)))
      nSamples / width
    }
    val logNewFeaturesTruncationDivisor = log(newFeaturesTruncationDivisor)
    val results = Array.ofDim[FeatureAllocation](nSamples)
    var b = 1
    tmAll {
    while (b <= nIterations) {
      for (i <- 0 until nItems) {
        val proposals = tmEnumeration { state.enumerateCombinationsFor(i) }
        val logWeights = proposals.map(fa => (fa, tmPrior { FeatureAllocationDistributions.logProbabilityIBP(fa,mass) } + tmLikelihood { lglfm.logLikelihood(fa) } ) )
        state = rdg.nextItem(logWeights, onLogScale = true)._1
        state = state.partitionBySingletonsOf(i)._2
        @scala.annotation.tailrec
        def engine(weights: List[(FeatureAllocation, Double)], max: Double): List[(FeatureAllocation, Double)] = {
          val newCumState = FeatureAllocation(weights.head._1).add(i)
          val newLogWeight = tmPrior { FeatureAllocationDistributions.logProbabilityIBP(newCumState,mass) } + tmLikelihood { lglfm.logLikelihood(newCumState) }
          val expanded = (newCumState, newLogWeight) :: weights
          if (newLogWeight < max - logNewFeaturesTruncationDivisor) expanded
          else engine(expanded, if (newLogWeight > max) newLogWeight else max)
        }
        val setup = (state, tmPrior { FeatureAllocationDistributions.logProbabilityIBP(state,mass) } + tmLikelihood { lglfm.logLikelihood(state) } ) :: Nil
        val weights = engine(setup, Double.NegativeInfinity).toIndexedSeq
        state = rdg.nextItem(weights, onLogScale = true)._1
      }
      if (b % thin == 0) {
        val index = (b-1)/thin
        results(index) = state
        if ( ( width > 0 ) && ( index % rate == 0 ) ) print("*")
      }
      b += 1
    }
    }
    if ( width > 0 ) println("]")
    println(tmAll)
    println(tmLikelihood)
    println(tmPrior)
    println(tmEnumeration)
    results
  }

}
