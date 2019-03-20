package org.ddahl.aibd.model.lineargaussian

import org.ddahl.matrix._
import org.apache.commons.math3.random.RandomDataGenerator
import org.apache.commons.math3.util.FastMath.log
import org.ddahl.commonsmath.RandomDataGeneratorImprovements
import org.ddahl.aibd.TimeMonitor

object PosteriorSimulation {

  def updateFeatureAllocationViaNeighborhoods(Z: Matrix, mass: Double, lglfm: LinearGaussianLatentFeatureModel, nSamples: Int, thin: Int, rdg: RandomDataGenerator, newFeaturesTruncationDivisor: Double = 1000): Array[Matrix] = {
    val tm = TimeMonitor()
    val tm3 = TimeMonitor()
    val nItems = lglfm.N
    var state = Z
    val nIterations = thin*nSamples
    val logNewFeaturesTruncationDivisor = log(newFeaturesTruncationDivisor)
    val results = Array.ofDim[Matrix](nSamples)
    var b = 1
    tm {
    while (b <= nIterations) {
      for (i <- 0 until nItems) {
        val (singletons, existing) = FeatureAllocationUtilities.partitionBySingletonsOf(i, state)
        val proposals = if (existing == null) Array(singletons)
        else tm3 { FeatureAllocationUtilities.enumerateCombinationsFor(i, singletons, existing) }
        val logWeights = proposals.map(Z => (Z, FeatureAllocationDistributions.logProbabilityIBP(Z,nItems,mass) + lglfm.logLikelihood(Z)))
        state = rdg.nextItem(logWeights, onLogScale = true)._1
        state = FeatureAllocationUtilities.partitionBySingletonsOf(i, state)._2
        val featureWithOnlyI = Array.ofDim[Double](nItems)
        featureWithOnlyI(i) = 1.0
        @scala.annotation.tailrec
        def engine(weights: List[(Matrix, Double)], max: Double): List[(Matrix, Double)] = {
          val newCumState = weights.head._1 | featureWithOnlyI
          val newLogWeight = FeatureAllocationDistributions.logProbabilityIBP(newCumState,nItems,mass) + lglfm.logLikelihood(newCumState)
          val expanded = (newCumState, newLogWeight) :: weights
          if (newLogWeight < max - logNewFeaturesTruncationDivisor) expanded
          else engine(expanded, if (newLogWeight > max) newLogWeight else max)
        }

        val setup = (state, FeatureAllocationDistributions.logProbabilityIBP(state,nItems,mass) + lglfm.logLikelihood(state)) :: Nil
        val weights = engine(setup, Double.NegativeInfinity).toIndexedSeq
        state = rdg.nextItem(weights, onLogScale = true)._1
      }
      if (b % thin == 0) results((b - 1) / thin) = state
      b += 1
    }
    }
    println(lglfm.tm)
    println(tm3)
    println(tm)
    results
  }

}
