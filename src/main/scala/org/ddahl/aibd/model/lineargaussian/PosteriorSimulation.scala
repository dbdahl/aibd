package org.ddahl.aibd.model.lineargaussian

import org.ddahl.aibd.repeat
import org.ddahl.aibd.{MCMCAcceptanceMonitor1, TimeMonitor}
import org.ddahl.aibd.Utils.harmonicNumber
import org.apache.commons.math3.random.RandomDataGenerator
import org.apache.commons.math3.util.FastMath.{sqrt, log}
import org.ddahl.commonsmath.RandomDataGeneratorImprovements

object PosteriorSimulation {

  def update4AIBD(featureAllocation: FeatureAllocation, featureAllocationPrior: FeatureAllocationDistribution, lglfm: LinearGaussianLatentFeatureModel, massPriorShape: Double, massPriorRate: Double, nPerShuffle: Int, maxStandardDeviationX: Double, maxStandardDeviationW: Double, sdProposedStandardDeviationX: Double, sdProposedStandardDeviationW: Double, corProposedSdXSdW: Double, nSamples: Int, thin: Int, progressWidth: Int, rdg: RandomDataGenerator, parallel: Boolean, rankOneUpdates: Boolean, newFeaturesTruncationDivisor: Double = 1000): (Array[FeatureAllocation], Array[Array[Double]]) = {
    var stateFA = featureAllocation
    var stateFAPrior = featureAllocationPrior
    var stateLGLFM = lglfm
    val monitorFA = MCMCAcceptanceMonitor1()
    val monitorFAMass = MCMCAcceptanceMonitor1()
    val monitorFAPermutation = MCMCAcceptanceMonitor1()
    val monitorLGLFM = MCMCAcceptanceMonitor1()
    val tmAll = TimeMonitor()
    val tmAllocation = TimeMonitor()
    val tmMass = TimeMonitor()
    val tmPermutation = TimeMonitor()
    val tmLGLFM = TimeMonitor()
    val nIterations = thin*nSamples
    val (width,rate) = if ( progressWidth <= 0 ) (0,1)
    else {
      val r = nSamples / progressWidth
      if ( r == 0 ) (nSamples, 1) else (nSamples / r, r)
    }
    if ( width > 0 ) print("[" + (" " * width) + "]" + ("\b" * (width + 1)) + "   ")
    val resultFA = Array.ofDim[FeatureAllocation](nSamples)
    val resultOthers = Array.ofDim[Double](nSamples,3)
    var b = 1
    while (b <= nIterations) {
      tmAll {
        stateFA = monitorFA(tmAllocation(updateFeatureAllocation(stateFA, stateFAPrior, stateLGLFM, rdg, parallel, rankOneUpdates, newFeaturesTruncationDivisor)))
        stateFAPrior = stateFAPrior match {
          case faPrior: AttractionIndianBuffetDistribution =>
            if ( nPerShuffle < 2 ) faPrior
            else monitorFAPermutation(tmPermutation(updatePermutation(stateFA, faPrior, rdg, nPerShuffle)))
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
        stateLGLFM = if ( ( sdProposedStandardDeviationX <= 0.0 ) || ( sdProposedStandardDeviationW <= 0.0 ) ) stateLGLFM
        else monitorLGLFM(tmLGLFM(updateSamplingModel(stateFA, stateFAPrior, stateLGLFM, rdg, maxStandardDeviationX, maxStandardDeviationW, sdProposedStandardDeviationX, sdProposedStandardDeviationW, corProposedSdXSdW)))
        if (b % thin == 0) {
          val index = (b-1)/thin
          resultFA(index) = stateFA
          resultOthers(index)(0) = stateFAPrior match {
            case faPrior: FeatureAllocationDistribution with HasMass[_] => faPrior.mass
            case _ => 0.0
          }
          resultOthers(index)(1) = stateLGLFM.standardDeviationX
          resultOthers(index)(2) = stateLGLFM.standardDeviationW
          if ( ( width > 0 ) && ( index % rate == 0 ) ) {
            print("\b"*3)
            print("*")
            print(" ")
            print("%2d".format(stateFA.nFeatures))
          }
        }
      }
      b += 1
    }
    if ( width > 0 ) {
      print("\b"*3)
      println("]  ")
      println("Parallel: " + parallel)
      println("Rank-one updates: " + rankOneUpdates)
      println("Lapse times:")
      println("  Total: " + tmAll)
      println("  Allocation: " + tmAllocation)
      println("  Permutation: " + tmPermutation)
      println("  Mass: " + tmMass)
      println("  Sampling model: " + tmLGLFM)
      println("Acceptance rates:")
      println("  Permutation: " + monitorFAPermutation.rate)
      println("  Sampling model: " + monitorLGLFM.rate)
    }
    (resultFA, resultOthers)
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

  def updateSamplingModel(featureAllocation: FeatureAllocation, featureAllocationPrior: FeatureAllocationDistribution, lglfm: LinearGaussianLatentFeatureModel, rdg: RandomDataGenerator, maxStandardDeviationX: Double, maxStandardDeviationW: Double, sdProposedStandardDeviationX: Double, sdProposedStandardDeviationW: Double, corProposedSdXSdW: Double): (LinearGaussianLatentFeatureModel, Int, Int) = {
    val m1 = lglfm.standardDeviationX
    val m2 = lglfm.standardDeviationW
    val proposedStandardDeviationX = rdg.nextGaussian(m1, sdProposedStandardDeviationX)
    if ( ( proposedStandardDeviationX < 0 ) || ( proposedStandardDeviationX > maxStandardDeviationX ) ) return (lglfm,0,1)
    val proposedStandardDeviationW = rdg.nextGaussian(m2 + (sdProposedStandardDeviationW/sdProposedStandardDeviationX) * corProposedSdXSdW * (proposedStandardDeviationX - m1), sqrt((1 - corProposedSdXSdW*corProposedSdXSdW)*sdProposedStandardDeviationW*sdProposedStandardDeviationW))
    if ( ( proposedStandardDeviationW < 0 ) || ( proposedStandardDeviationW > maxStandardDeviationW ) ) return (lglfm,0,1)
    val proposedLGLFM = lglfm.updateStandardDeviations(proposedStandardDeviationX, proposedStandardDeviationW)
    val diff = proposedLGLFM.logLikelihood(featureAllocation) - lglfm.logLikelihood(featureAllocation)
    if ( ( diff > 0 ) || ( log(rdg.nextUniform(0,1)) < diff ) ) (proposedLGLFM,1,1) else (lglfm,0,1)
  }

  def updateFeatureAllocation(featureAllocation: FeatureAllocation, featureAllocationPrior: FeatureAllocationDistribution, lglfm: LinearGaussianLatentFeatureModel, rdg: RandomDataGenerator, parallel: Boolean, rankOneUpdates: Boolean, newFeaturesTruncationDivisor: Double = 1000): (FeatureAllocation, Int, Int) = {
    val nItems = lglfm.N
    val logNewFeaturesTruncationDivisor = log(newFeaturesTruncationDivisor)
    var state = featureAllocation
    var accepts = 0
    var attempts = 0
    for (i <- 0 until nItems) {
      val (stateNew, n, d) = updateFeatureAllocationOfExistingByEnumeration(i, state, featureAllocationPrior, lglfm, rdg, parallel, rankOneUpdates)
      // val (stateNew, n, d) = updateFeatureAllocationOfExistingWholeRow(i, state, featureAllocationPrior, lglfm, rdg)
      state = stateNew
      accepts += n
      attempts += d
      state = updateFeatureAllocationOfSingletons(i, state, featureAllocationPrior, lglfm, rdg, logNewFeaturesTruncationDivisor)
    }
    (state, accepts, attempts)
  }

  private def logPosterior0(i: Int, fa: FeatureAllocation, featureAllocationPrior: FeatureAllocationDistribution, lglfm: LinearGaussianLatentFeatureModel): Double = {
    featureAllocationPrior.logProbability(i, fa) + lglfm.logLikelihood(fa)
  }

  private def logPosterior1(i: Int, fa: FeatureAllocation, featureAllocationPrior: FeatureAllocationDistribution, lglfm: LinearGaussianLatentFeatureModel): (FeatureAllocation, Double) = {
    (fa, featureAllocationPrior.logProbability(i, fa) + lglfm.logLikelihood(fa) )
  }

  private def logPosterior2(i: Int, fa: FeatureAllocation, featureAllocationPrior: FeatureAllocationDistribution, lglfm: LinearGaussianLatentFeatureModel, lc: LikelihoodComponents): (FeatureAllocation, Double) = {
    (fa, featureAllocationPrior.logProbability(i, fa) + lglfm.logLikelihood(lglfm.allocateFeaturesFor(i,lc,fa.matrix(i))) )
  }

  def updateFeatureAllocationOfSingletons(i: Int, featureAllocation: FeatureAllocation, featureAllocationPrior: FeatureAllocationDistribution, lglfm: LinearGaussianLatentFeatureModel, rdg: RandomDataGenerator, logNewFeaturesTruncationDivisor: Double) = {
    @scala.annotation.tailrec
    def engine(weights: List[(FeatureAllocation, Double)], max: Double): List[(FeatureAllocation, Double)] = {
      val newCumState = FeatureAllocation(weights.head._1).add(i)
      val newLogWeight = logPosterior1(i, newCumState, featureAllocationPrior, lglfm)._2
      val expanded = (newCumState, newLogWeight) :: weights
      if (newLogWeight < max - logNewFeaturesTruncationDivisor) expanded
      else engine(expanded, if (newLogWeight > max) newLogWeight else max)
    }
    val existing = featureAllocation.partitionBySingletonsOf(i)._2
    val weights = engine(logPosterior1(i, existing, featureAllocationPrior, lglfm) :: Nil, Double.NegativeInfinity).toIndexedSeq
    rdg.nextItem(weights, onLogScale = true)._1
  }

  def updateFeatureAllocationOfExistingByEnumeration(i: Int, featureAllocation: FeatureAllocation, featureAllocationPrior: FeatureAllocationDistribution, lglfm: LinearGaussianLatentFeatureModel, rdg: RandomDataGenerator, parallel: Boolean, rankOneUpdates: Boolean): (FeatureAllocation, Int, Int) = {
    val proposals = featureAllocation.enumerateFor(i)
    val logWeights = {
      if (rankOneUpdates) {
        if (proposals.isEmpty) Array[(FeatureAllocation, Double)]()
        else {
          val lc1 = lglfm.computeLikelihoodComponents(proposals.head)
          val lc2 = lglfm.deallocateFeaturesFor(i, lc1)
          if (parallel) proposals.par.map(logPosterior2(i, _, featureAllocationPrior, lglfm, lc2)).toArray else proposals.map(logPosterior2(i, _, featureAllocationPrior, lglfm, lc2))
        }
      } else {
        if (parallel) proposals.par.map(logPosterior1(i, _, featureAllocationPrior, lglfm)).toArray else proposals.map(logPosterior1(i, _, featureAllocationPrior, lglfm))
      }
    }
    (rdg.nextItem(logWeights, onLogScale = true)._1, 1, 1)
  }

  def updateFeatureAllocationOfExistingSimplyOff(i: Int, featureAllocation: FeatureAllocation, featureAllocationPrior: FeatureAllocationDistribution, lglfm: LinearGaussianLatentFeatureModel, rdg: RandomDataGenerator): (FeatureAllocation, Int, Int) = {
    println("---Start---")
    println("i: "+i)
    println("State:")
    println(featureAllocation)
    var (singletons, state) = featureAllocation.partitionBySingletonsOf(i)
    println("Singletons:")
    println(singletons)
    println("Existing:")
    println(state)
    var accepts = 0
    var attempts = 0
    repeat(state.nFeatures) {
      val j = rdg.nextInt(0,state.nFeatures-1)
      val proposal = if ( rdg.nextUniform(0,1) < 0.5 ) {
        println("add: "+i+" "+j)
        state.add(i,j)
      } else {
        println("remove: "+i+" "+j)
        state.remove(i,j)
      }
      val fullProposal = proposal add singletons
      val fullState = state add singletons
      println("Full State:\n"+fullState)
      println("Full Proposal:\n"+fullProposal)
      fullProposal.check()
      val diff = logPosterior0(i, fullProposal, featureAllocationPrior, lglfm) - logPosterior0(i, fullState, featureAllocationPrior, lglfm)
      attempts += 1
      if ( ( diff >= 0 ) || ( log(rdg.nextUniform(0.0,1.0)) < diff ) ) {
        println("Accept!")
        accepts += 1
        state = proposal
      }
    }
    (state add singletons, accepts, attempts)
  }

  def updateFeatureAllocationOfExistingSimply1(i: Int, featureAllocation: FeatureAllocation, featureAllocationPrior: FeatureAllocationDistribution, lglfm: LinearGaussianLatentFeatureModel, rdg: RandomDataGenerator): (FeatureAllocation, Int, Int) = {
    var (singletons, state) = featureAllocation.partitionBySingletonsOf(i)
    var accepts = 0
    var attempts = 0
    for ( j <- 0 until state.nFeatures ) {
      val proposal = if ( state.features(j).contains(i) ) state.remove(i,j) else state.add(i,j)
      val diff = logPosterior0(i, proposal add singletons, featureAllocationPrior, lglfm) - logPosterior0(i, state add singletons, featureAllocationPrior, lglfm)
      attempts += 1
      if ( ( diff >= 0 ) || ( log(rdg.nextUniform(0.0,1.0)) < diff ) ) {
        accepts += 1
        state = proposal
      }
    }
    (state add singletons, accepts, attempts)
  }

  def updateFeatureAllocationOfExistingSimplyBroken(i: Int, featureAllocation: FeatureAllocation, featureAllocationPrior: FeatureAllocationDistribution, lglfm: LinearGaussianLatentFeatureModel, rdg: RandomDataGenerator): (FeatureAllocation, Int, Int) = {
    var state = featureAllocation
    var accepts = 0
    var attempts = 0
    for ( j <- 0 until state.nFeatures; if ! state.isSingleton(i,j) ) {
      val proposal = if ( state.features(j).contains(i) ) state.remove(i,j) else state.add(i,j)
      val diff = logPosterior0(i, proposal, featureAllocationPrior, lglfm) - logPosterior0(i, state, featureAllocationPrior, lglfm)
      attempts += 1
      if ( ( diff >= 0 ) || ( log(rdg.nextUniform(0.0,1.0)) < diff ) ) {
        accepts += 1
        state = proposal
      }
    }
    (state, accepts, attempts)
  }

  def updateFeatureAllocationOfExistingSimplySwapBroken(i: Int, featureAllocation: FeatureAllocation, featureAllocationPrior: FeatureAllocationDistribution, lglfm: LinearGaussianLatentFeatureModel, rdg: RandomDataGenerator): (FeatureAllocation, Int, Int) = {
    if ( featureAllocation.partitionBySingletonsOf(i)._2.nFeatures < 2 ) return (featureAllocation,0,0)
    var state = featureAllocation
    var accepts = 0
    var attempts = 0
    repeat(10) {
      val (j1, j2) = {
        @scala.annotation.tailrec
        def get(taken: Int): Int = {
          val j = rdg.nextInt(0,state.nFeatures-1)
          if ( ( j == taken ) || state.isSingleton(i,j) ) get(taken) else j
        }
        val j1 = get(-1)
        (j1, get(j1))
      }
      val proposal = if ( state.features(j1).contains(i) == state.features(j2).contains(i) ) state
      else {
        if ( state.features(j1).contains(i) ) state.remove(i,j1).add(i,j2)
        else state.remove(i,j2).add(i,j1)
      }
      val diff = logPosterior0(i, proposal, featureAllocationPrior, lglfm) - logPosterior0(i, state, featureAllocationPrior, lglfm)
      attempts += 1
      if ( ( diff >= 0 ) || ( log(rdg.nextUniform(0.0,1.0)) < diff ) ) {
        accepts += 1
        state = proposal
      }
    }
    (state, accepts, attempts)
  }

  def updateFeatureAllocationOfExistingWholeRowBroken(i: Int, featureAllocation: FeatureAllocation, featureAllocationPrior: FeatureAllocationDistribution, lglfm: LinearGaussianLatentFeatureModel, rdg: RandomDataGenerator): (FeatureAllocation, Int, Int) = {
    var state = featureAllocation
    var accepts = 0
    var attempts = 0
    val p = 0.5
    val debug = false
    if ( debug ) println("i: "+i)
    if ( debug ) println("State:\n"+state)
    repeat(20) {
      if ( debug ) println("*")
      var proposal = FeatureAllocation(state.matrix.map(_.clone))
      for ( j <- 0 until state.nFeatures; if ! state.isSingleton(i,j) ) {
        proposal = if ( rdg.nextUniform(0,1) <= p ) proposal.add(i,j) else proposal.remove(i,j)
      }
      if ( debug ) println("State:\n"+state)
      if ( debug ) println("Proposal:\n"+proposal)
      val diff = logPosterior0(i, proposal, featureAllocationPrior, lglfm) - logPosterior0(i, state, featureAllocationPrior, lglfm)
      if ( debug ) println("Diff: "+diff)
      attempts += 1
      if ( rdg.nextUniform(0.0,1.0) < math.exp(diff) ) {
        if ( debug ) println("Accept")
        accepts += 1
        state = proposal
      }
    }
    (state, accepts, attempts)
  }

  def updateFeatureAllocationOfExistingWholeRowBroken2(i: Int, featureAllocation: FeatureAllocation, featureAllocationPrior: FeatureAllocationDistribution, lglfm: LinearGaussianLatentFeatureModel, rdg: RandomDataGenerator): (FeatureAllocation, Int, Int) = {
    var state = featureAllocation
    var accepts = 0
    var attempts = 0
    val p = 0.5
    val debug = false
    if ( debug ) println("i: "+i)
    if ( debug ) println("State:\n"+state)
    repeat(20) {
      if ( debug ) println("*")
      var proposal = state.removeRow(i)
      for ( j <- 0 until proposal.nFeatures ) {
        proposal = if ( rdg.nextUniform(0,1) <= p ) proposal.add(i,j) else proposal.remove(i,j)
      }
      repeat(rdg.nextInt(0,6)) {
        proposal = proposal.add(i)
      }
      if ( debug ) println("State:\n"+state)
      if ( debug ) println("Proposal:\n"+proposal)
      val diff = logPosterior0(i, proposal, featureAllocationPrior, lglfm) - logPosterior0(i, state, featureAllocationPrior, lglfm)
      if ( debug ) println("Diff: "+diff)
      attempts += 1
      if ( rdg.nextUniform(0.0,1.0) < math.exp(diff) ) {
        if ( debug ) println("Accept")
        accepts += 1
        state = proposal
      }
    }
    (state, accepts, attempts)
  }

  def updateFeatureAllocationOfExistingSimplyTies(i: Int, featureAllocation: FeatureAllocation, featureAllocationPrior: FeatureAllocationDistribution, lglfm: LinearGaussianLatentFeatureModel, rdg: RandomDataGenerator): (FeatureAllocation, Int, Int) = {
    val state = featureAllocation
    var accepts = 0
    var attempts = 0
    val existing = state.asList.filter { case (feature,size,count) => ! ( ( size == 1 ) || feature.contains(i) ) }
    val nExistingFeatures = existing.size
    if ( nExistingFeatures > 0 ) {
      val which = rdg.nextInt(0,nExistingFeatures-1)
      val current = existing(which)
    }

    /*
    for ( pair <- stateMap.keys; if ! ( ( pair._2 == 1 ) && ( pair._1.contains(i) ) ) ) {
      val proposal = if ( pair._1.contains(i) ) state.remove(i,j) else state.add(i,j)
      state.asMap
      val diff = logPosterior0(i, proposal, featureAllocationPrior, lglfm) - logPosterior0(i, state, featureAllocationPrior, lglfm)
      attempts += 1
      if ( ( diff >= 0 ) || ( log(rdg.nextUniform(0.0,1.0)) < diff ) ) {
        accepts += 1
        state = proposal
      }
    }
    */
    (state, accepts, attempts)
  }

}

