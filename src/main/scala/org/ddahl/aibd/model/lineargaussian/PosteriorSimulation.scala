package org.ddahl.aibd.model.lineargaussian

import org.apache.commons.math3.random.RandomDataGenerator
import org.apache.commons.math3.util.FastMath.log
import org.ddahl.commonsmath.RandomDataGeneratorImprovements
import org.ddahl.aibd.TimeMonitor

object PosteriorSimulation {

  def updateFeatureAllocationViaNeighborhoods(featureAllocation: FeatureAllocation, mass: Double, lglfm: LinearGaussianLatentFeatureModel, nSamples: Int, thin: Int, width: Int, rdg: RandomDataGenerator, parallel: Boolean, newFeaturesTruncationDivisor: Double = 1000): Array[FeatureAllocation] = {
    val nItems = lglfm.N
    val logNewFeaturesTruncationDivisor = log(newFeaturesTruncationDivisor)
    var state = featureAllocation
    val tmAll = TimeMonitor()
    val tmPosterior1 = TimeMonitor()
    val tmPosterior2 = TimeMonitor()
    val tmEnumeration = TimeMonitor()
    val nIterations = thin*nSamples
    val rate = if ( width <= 0 ) 1
    else {
      print("[" + (" " * width) + "]" + ("\b" * (width + 1)))
      nSamples / width
    }
    def logPosterior(fa: FeatureAllocation): (FeatureAllocation, Double) = {
      (fa, FeatureAllocationDistributions.logProbabilityIBP(fa,mass) + lglfm.logLikelihood(fa))
    }
    val results = Array.ofDim[FeatureAllocation](nSamples)
    var b = 1
    tmAll {
    while (b <= nIterations) {
      for (i <- 0 until nItems) {
        val proposals = tmEnumeration { state.enumerateCombinationsFor(i) }
        val logWeights = tmPosterior1 { if ( parallel ) proposals.par.map(logPosterior).toArray else proposals.map(logPosterior) }
        state = rdg.nextItem(logWeights, onLogScale = true)._1
        state = state.partitionBySingletonsOf(i)._2
        @scala.annotation.tailrec
        def engine(weights: List[(FeatureAllocation, Double)], max: Double): List[(FeatureAllocation, Double)] = {
          val newCumState = FeatureAllocation(weights.head._1).add(i)
          val newLogWeight = logPosterior(newCumState)._2
          val expanded = (newCumState, newLogWeight) :: weights
          if (newLogWeight < max - logNewFeaturesTruncationDivisor) expanded
          else engine(expanded, if (newLogWeight > max) newLogWeight else max)
        }
        val weights = tmPosterior2 { engine(logPosterior(state) :: Nil, Double.NegativeInfinity).toIndexedSeq }
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
    println("Parallel: "+parallel)
    println(tmAll)
    println(tmPosterior1)
    println(tmPosterior2)
    println(tmEnumeration)
    results
  }

}
