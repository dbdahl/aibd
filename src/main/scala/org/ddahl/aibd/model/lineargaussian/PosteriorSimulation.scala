package org.ddahl.aibd.model.lineargaussian

import org.apache.commons.math3.random.RandomDataGenerator
import org.apache.commons.math3.util.FastMath.log
import org.ddahl.commonsmath.RandomDataGeneratorImprovements
import org.ddahl.aibd.{AttractionIndianBuffetDistribution, Permutation, Similarity, TimeMonitor}

object PosteriorSimulation {

  def mkLogPriorProbabilityIBP(mass: Double): (Int, FeatureAllocation) => Double = (i: Int, fa: FeatureAllocation) => {
    FeatureAllocationDistributions.logProbabilityIBP(fa, mass)
  }

  def mkLogPriorProbabilityAIBD(mass: Double, permutation: Permutation, similarity: Similarity): (Int, FeatureAllocation) => Double = {
    val aibd = AttractionIndianBuffetDistribution(mass, permutation, similarity)
    (i: Int, fa: FeatureAllocation) => {
      val faa = fa.convertToAlternativeImplementation
      aibd.logDensityStartingFromIndex(i, faa, false)
    }
  }

  def updateFeatureAllocationViaNeighborhoods(featureAllocation: FeatureAllocation, logPriorProbability: (Int, FeatureAllocation) => Double, lglfm: LinearGaussianLatentFeatureModel, nSamples: Int, thin: Int, progressWidth: Int, rdg: RandomDataGenerator, parallel: Boolean, rankOneUpdates: Boolean, newFeaturesTruncationDivisor: Double = 1000): Array[FeatureAllocation] = {
    val nItems = lglfm.N
    val logNewFeaturesTruncationDivisor = log(newFeaturesTruncationDivisor)
    var state = featureAllocation
    val tmAll = TimeMonitor()
    val tmPosterior1 = TimeMonitor()
    val tmPosterior2 = TimeMonitor()
    val tmEnumeration = TimeMonitor()
    val nIterations = thin*nSamples
    val (width,rate) = if ( progressWidth <= 0 ) (0,1)
    else {
      val r = nSamples / progressWidth
      if ( r == 0 ) (nSamples, 1) else (progressWidth, r)
    }
    if ( width > 0 ) print("[" + (" " * width) + "]" + ("\b" * (width + 1)))
    def logPosterior1(i: Int, fa: FeatureAllocation): (FeatureAllocation, Double) = {
      (fa, logPriorProbability(i, fa) + lglfm.logLikelihood(fa))
    }
    def logPosterior2(i: Int, fa: FeatureAllocation, lc: LikelihoodComponents): (FeatureAllocation, Double) = {
      (fa, logPriorProbability(i, fa) + lglfm.logLikelihood(lglfm.addFeaturesFor(i,lc,fa.matrix(i))))
    }
    val results = Array.ofDim[FeatureAllocation](nSamples)
    var b = 1
    tmAll {
    while (b <= nIterations) {
      for (i <- 0 until nItems) {
        val proposals = tmEnumeration { state.enumerateCombinationsFor(i) }
        val logWeights = if ( rankOneUpdates ) {
          if ( proposals.isEmpty ) Array[(FeatureAllocation,Double)]()
          else {
            val lc1 = lglfm.computeLikelihoodComponents(proposals.head)
            val lc2 = lglfm.dropFeaturesFor(i,lc1)
            if ( parallel ) proposals.par.map(logPosterior2(i,_,lc2)).toArray else proposals.map(logPosterior2(i,_,lc2))
          }
        } else {
          tmPosterior1 {
            if (parallel) proposals.par.map(logPosterior1(i,_)).toArray else proposals.map(logPosterior1(i,_))
          }
        }
        state = rdg.nextItem(logWeights, onLogScale = true)._1
        state = state.partitionBySingletonsOf(i)._2
        @scala.annotation.tailrec
        def engine(weights: List[(FeatureAllocation, Double)], max: Double): List[(FeatureAllocation, Double)] = {
          val newCumState = FeatureAllocation(weights.head._1).add(i)
          val newLogWeight = logPosterior1(i,newCumState)._2
          val expanded = (newCumState, newLogWeight) :: weights
          if (newLogWeight < max - logNewFeaturesTruncationDivisor) expanded
          else engine(expanded, if (newLogWeight > max) newLogWeight else max)
        }
        val weights = tmPosterior2 { engine(logPosterior1(i,state) :: Nil, Double.NegativeInfinity).toIndexedSeq }
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
    println("Parallel:         "+parallel)
    println("Rank-one updates: "+rankOneUpdates)
    println(tmAll)
    println(tmPosterior1)
    println(tmPosterior2)
    println(tmEnumeration)
    results
  }

}
