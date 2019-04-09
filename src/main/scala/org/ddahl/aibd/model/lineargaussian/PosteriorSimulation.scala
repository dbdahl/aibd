package org.ddahl.aibd.model.lineargaussian

import org.ddahl.aibd.{MCMCAcceptanceMonitor1, TimeMonitor}
import org.ddahl.aibd.Utils.harmonicNumber
import org.apache.commons.math3.random.RandomDataGenerator
import org.apache.commons.math3.util.FastMath.log
import org.ddahl.commonsmath.RandomDataGeneratorImprovements

object PosteriorSimulation {

  def update4AIBD(featureAllocation: FeatureAllocation, featureAllocationPrior: FeatureAllocationDistribution, lglfm: LinearGaussianLatentFeatureModel, massPriorShape: Double, massPriorRate: Double, nSamples: Int, thin: Int, progressWidth: Int, rdg: RandomDataGenerator, nPerShuffle: Int, parallel: Boolean, rankOneUpdates: Boolean, newFeaturesTruncationDivisor: Double = 1000): (Array[FeatureAllocation], Array[Double]) = {
    var stateFA = featureAllocation
    var stateFAPrior = featureAllocationPrior
    val monitorFA = MCMCAcceptanceMonitor1()
    val monitorFAMass = MCMCAcceptanceMonitor1()
    val monitorFAPermutation = MCMCAcceptanceMonitor1()
    val tmAll = TimeMonitor()
    val tmAllocation = TimeMonitor()
    val tmMass = TimeMonitor()
    val tmPermutation = TimeMonitor()
    val nIterations = thin*nSamples
    val (width,rate) = if ( progressWidth <= 0 ) (0,1)
    else {
      val r = nSamples / progressWidth
      if ( r == 0 ) (nSamples, 1) else (progressWidth, r)
    }
    if ( width > 0 ) print("[" + (" " * width) + "]" + ("\b" * (width + 1)))
    val resultFA = Array.ofDim[FeatureAllocation](nSamples)
    val resultMass = Array.ofDim[Double](nSamples)
    var b = 1
    tmAll {
      while (b <= nIterations) {
        stateFA = monitorFA(tmAllocation(updateFeatureAllocationViaNeighborhoods(stateFA, stateFAPrior, lglfm, rdg, parallel, rankOneUpdates, newFeaturesTruncationDivisor)))
        stateFAPrior = stateFAPrior match {
          case faPrior: AttractionIndianBuffetDistribution =>
            monitorFAPermutation(tmPermutation(updatePermutation(stateFA, faPrior, rdg, nPerShuffle)))
          case _ =>
            stateFAPrior
        }
        stateFAPrior = stateFAPrior match {
          case faPrior: FeatureAllocationDistribution with HasMass[_] =>
            if ( ( massPriorShape <= 0 ) || ( massPriorRate <= 0 ) ) faPrior
            else monitorFAMass(tmMass(updateMass(stateFA, faPrior, rdg, massPriorShape, massPriorRate)))
          case _ =>
            stateFAPrior
        }
        if (b % thin == 0) {
          val index = (b-1)/thin
          resultFA(index) = stateFA
          resultMass(index) = stateFAPrior match {
            case faPrior: FeatureAllocationDistribution with HasMass[_] => faPrior.mass
            case _ => 0.0
          }
          if ( ( width > 0 ) && ( index % rate == 0 ) ) print("*")
        }
        b += 1
      }
    }
    if ( width > 0 ) println("]")
    println("Parallel: "+parallel)
    println("Rank-one updates: "+rankOneUpdates)
    println("Main lapse time: "+tmAll)
    println("Allocation lapse time: "+tmAllocation)
    println("Permutation lapse time: "+tmPermutation)
    println("Mass lapse time: "+tmMass)
    println("Permutation acceptance rate: "+monitorFAPermutation.rate)
    (resultFA, resultMass)
  }

  def updatePermutation(featureAllocation: FeatureAllocation, aibdPrior: AttractionIndianBuffetDistribution, rdg: RandomDataGenerator, nPerShuffle: Int): (AttractionIndianBuffetDistribution, Int, Int) = {
    val proposedPermutation = aibdPrior.permutation.nPerShuffle(nPerShuffle).shuffle(rdg)
    val proposedAIBDPrior = aibdPrior.updatePermutation(proposedPermutation)
    val diff = proposedAIBDPrior.logProbability(featureAllocation) - aibdPrior.logProbability(featureAllocation)
    if ( ( diff > 0 ) || ( log(rdg.nextUniform(0,1)) < diff ) ) (proposedAIBDPrior,1,1) else (aibdPrior,0,1)
  }

  def updateMass[T](featureAllocation: FeatureAllocation, hasMassPrior: FeatureAllocationDistribution with HasMass[T], rdg: RandomDataGenerator, priorShape: Double, priorRate: Double): (FeatureAllocationDistribution with HasMass[T], Int, Int) = {
    val posteriorShape = priorShape + featureAllocation.nFeatures
    val posteriorRate = priorRate + harmonicNumber(featureAllocation.nItems)
    (hasMassPrior.updateMass(rdg.nextGamma(posteriorShape, 1/posteriorRate)), 1, 1)
  }

  def updateFeatureAllocationViaNeighborhoods(featureAllocation: FeatureAllocation, featureAllocationPrior: FeatureAllocationDistribution, lglfm: LinearGaussianLatentFeatureModel, rdg: RandomDataGenerator, parallel: Boolean, rankOneUpdates: Boolean, newFeaturesTruncationDivisor: Double = 1000): (FeatureAllocation, Int, Int) = {
    val nItems = lglfm.N
    val logNewFeaturesTruncationDivisor = log(newFeaturesTruncationDivisor)
    var state = featureAllocation
    def logPosterior1(i: Int, fa: FeatureAllocation): (FeatureAllocation, Double) = {
      (fa, featureAllocationPrior.logProbability(i, fa) + lglfm.logLikelihood(fa) )
    }
    def logPosterior2(i: Int, fa: FeatureAllocation, lc: LikelihoodComponents): (FeatureAllocation, Double) = {
      (fa, featureAllocationPrior.logProbability(i, fa) + lglfm.logLikelihood(lglfm.allocateFeaturesFor(i,lc,fa.matrix(i))) )
    }
    for (i <- 0 until nItems) {
      val proposals = state.enumerateFor(i)
      val logWeights = {
        if (rankOneUpdates) {
          if (proposals.isEmpty) Array[(FeatureAllocation, Double)]()
          else {
            val lc1 = lglfm.computeLikelihoodComponents(proposals.head)
            val lc2 = lglfm.deallocateFeaturesFor(i, lc1)
            if (parallel) proposals.par.map(logPosterior2(i, _, lc2)).toArray else proposals.map(logPosterior2(i, _, lc2))
          }
        } else {
          if (parallel) proposals.par.map(logPosterior1(i, _)).toArray else proposals.map(logPosterior1(i, _))
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
      val weights = engine(logPosterior1(i,state) :: Nil, Double.NegativeInfinity).toIndexedSeq
      state = rdg.nextItem(weights, onLogScale = true)._1
    }
    (state, 1, 1)
  }

}

