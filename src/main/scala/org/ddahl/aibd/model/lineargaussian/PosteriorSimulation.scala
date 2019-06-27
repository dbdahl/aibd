package org.ddahl.aibd.model.lineargaussian

import org.ddahl.aibd._
import distribution._
import org.apache.commons.math3.distribution.GammaDistribution
import util._
import util.Functions.{harmonicNumber, logOnInt}
import org.ddahl.commonsmath.RandomDataGeneratorImprovements
import org.apache.commons.math3.random.RandomDataGenerator
import org.apache.commons.math3.util.FastMath.{log, sqrt}
import scala.collection.parallel.immutable.ParVector

object PosteriorSimulation {

  def update4AIBD(featureAllocation: FeatureAllocation, featureAllocationPrior: FeatureAllocationDistribution, lglfm: LinearGaussianLatentFeatureModel, massPriorShape: Double, massPriorRate: Double, nPerShuffle: Int, temperaturePriorShape: Double, temperaturePriorRate: Double, maxStandardDeviationX: Double, maxStandardDeviationW: Double, sdProposedTemperature: Double, sdProposedStandardDeviationX: Double, sdProposedStandardDeviationW: Double, corProposedSdXSdW: Double, nOtherUpdatesPerAllocationUpdate: Int, nSamples: Int, thin: Int, progressWidth: Int, rdg: RandomDataGenerator, rankOneUpdates: Boolean, newFeaturesTruncationDivisor: Double = 1000): (Array[FeatureAllocation], Array[Array[Int]], Array[Array[Double]]) = {
    var stateFA = featureAllocation
    var stateFAPrior = featureAllocationPrior
    var stateLGLFM = lglfm
    val monitorFA = MCMCAcceptanceMonitor1()
    val monitorFAMass = MCMCAcceptanceMonitor1()
    val monitorFAPermutation = MCMCAcceptanceMonitor1()
    val monitorFATemperature = MCMCAcceptanceMonitor1()
    val monitorLGLFM = MCMCAcceptanceMonitor1()
    val tmAll = TimeMonitor()
    val tmAllocation = TimeMonitor()
    val tmMass = TimeMonitor()
    val tmPermutation = TimeMonitor()
    val tmTemperature = TimeMonitor()
    val tmLGLFM = TimeMonitor()
    val nIterations = thin*nSamples
    val (width,rate) = if ( progressWidth <= 0 ) (0,1)
    else {
      val r = nSamples / progressWidth
      if ( r == 0 ) (nSamples, 1) else (nSamples / r, r)
    }
    if ( width > 0 ) print("[" + (" " * width) + "]" + ("\b" * (width + 1)) + "   ")
    val resultFA = Array.ofDim[FeatureAllocation](nSamples)
    val resultPermutation = Array.ofDim[Int](nSamples, stateFA.nItems)
    val resultOthers = Array.ofDim[Double](nSamples,4)
    var b = 1
    while (b <= nIterations) {
      tmAll {
        stateFA = monitorFA(tmAllocation(updateFeatureAllocation(stateFA, stateFAPrior, stateLGLFM, rdg, rankOneUpdates, newFeaturesTruncationDivisor)))
        Misc.repeat(nOtherUpdatesPerAllocationUpdate) {
          stateFAPrior = stateFAPrior match {
            case faPrior: FeatureAllocationDistribution with HasMass[_] =>
              if ( ( massPriorShape <= 0 ) || ( massPriorRate <= 0 ) ) faPrior
              else monitorFAMass(tmMass(updateMass(stateFA, faPrior, rdg, massPriorShape, massPriorRate)))
            case _ =>
              stateFAPrior
          }
          stateFAPrior = stateFAPrior match {
            case faPrior: AttractionIndianBuffetDistribution =>
              if (nPerShuffle < 2) faPrior
              else monitorFAPermutation(tmPermutation(updatePermutation(stateFA, faPrior, rdg, nPerShuffle)))
            case _ =>
              stateFAPrior
          }
          stateFAPrior = stateFAPrior match {
            case faPrior: AttractionIndianBuffetDistribution =>
              if ((temperaturePriorShape <= 0) || (temperaturePriorRate <= 0)) faPrior
              else monitorFATemperature(tmTemperature(updateTemperature(stateFA, faPrior, rdg, temperaturePriorShape, temperaturePriorRate, sdProposedTemperature)))
            case _ =>
              stateFAPrior
          }
          stateLGLFM = if ((sdProposedStandardDeviationX <= 0.0) || (sdProposedStandardDeviationW <= 0.0)) stateLGLFM
          else monitorLGLFM(tmLGLFM(updateSamplingModel(stateFA, stateFAPrior, stateLGLFM, rdg, maxStandardDeviationX, maxStandardDeviationW, sdProposedStandardDeviationX, sdProposedStandardDeviationW, corProposedSdXSdW)))
        }
        if (b % thin == 0) {
          val index = (b-1)/thin
          resultFA(index) = stateFA
          resultPermutation(index) = stateFAPrior match {
            case faPrior: AttractionIndianBuffetDistribution => faPrior.permutation.toArray
            case _ => Array.ofDim[Int](stateFA.nItems)
          }
          resultOthers(index)(0) = stateFAPrior match {
            case faPrior: FeatureAllocationDistribution with HasMass[_] => faPrior.mass
            case _ => 0.0
          }
          resultOthers(index)(1) = stateFAPrior match {
            case faPrior: AttractionIndianBuffetDistribution =>
              faPrior.similarity match {
                case hasTemperature: Similarity with HasTemperature[_] => hasTemperature.temperature
                case _ => 0.0
              }
            case _ => 0.0
          }
          resultOthers(index)(2) = stateLGLFM.standardDeviationX
          resultOthers(index)(3) = stateLGLFM.standardDeviationW
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
      println("Rank-one updates: " + rankOneUpdates)
      println("Lapse times:")
      println("  Total: " + tmAll)
      println("  Allocation: " + tmAllocation)
      println("  Mass: " + tmMass)
      println("  Permutation: " + tmPermutation)
      println("  Temperature: " + tmTemperature)
      println("  Sampling model: " + tmLGLFM)
      println("Acceptance rates:")
      println("  Allocation: " + monitorFA.rate)
      println("  Permutation: " + monitorFAPermutation.rate)
      println("  Temperature: " + monitorFATemperature.rate)
      println("  Sampling model: " + monitorLGLFM.rate)
    }
    (resultFA, resultPermutation, resultOthers)
  }

  def updateMass[T](featureAllocation: FeatureAllocation, hasMassPrior: FeatureAllocationDistribution with HasMass[T], rdg: RandomDataGenerator, priorShape: Double, priorRate: Double): (FeatureAllocationDistribution with HasMass[T], Int, Int) = {
    val posteriorShape = priorShape + featureAllocation.nFeatures
    val posteriorRate = priorRate + harmonicNumber(featureAllocation.nItems)
    (hasMassPrior.updateMass(rdg.nextGamma(posteriorShape, 1/posteriorRate)), 1, 1)
  }

  def updatePermutation(featureAllocation: FeatureAllocation, aibdPrior: AttractionIndianBuffetDistribution, rdg: RandomDataGenerator, nPerShuffle: Int): (AttractionIndianBuffetDistribution, Int, Int) = {
    val proposedPermutation = aibdPrior.permutation.nPerShuffle(nPerShuffle).shuffle(rdg)
    val proposedAIBDPrior = aibdPrior.updatePermutation(proposedPermutation)
    val logMHRatio = proposedAIBDPrior.logProbability(featureAllocation) - aibdPrior.logProbability(featureAllocation)
    if ( ( logMHRatio > 0 ) || ( log(rdg.nextUniform(0,1)) < logMHRatio ) ) (proposedAIBDPrior,1,1) else (aibdPrior,0,1)
  }

  def updateTemperature(featureAllocation: FeatureAllocation, aibdPrior: AttractionIndianBuffetDistribution, rdg: RandomDataGenerator, temperatureShape: Double, temperatureRate: Double, sdProposedTemperature: Double): (FeatureAllocationDistribution, Int, Int) = {
    val (currentTemperature, proposedTemperature, proposedAIBDPrior) = aibdPrior.similarity match {
      case x: Similarity with HasTemperature[_] =>
        val proposedTemperature = rdg.nextGaussian(x.temperature, sdProposedTemperature)
        if ( proposedTemperature < 0 ) return (aibdPrior,0,1)
        (x.temperature, proposedTemperature, aibdPrior.updateSimilarity(x.updateTemperature(proposedTemperature)))
      case _ => return (aibdPrior, 0, 0)
    }
    val gammaDistribution = new GammaDistribution(null, temperatureShape, 1/temperatureRate)
    val logMHRatio = proposedAIBDPrior.logProbability(featureAllocation) + gammaDistribution.logDensity(proposedTemperature) - aibdPrior.logProbability(featureAllocation) - gammaDistribution.logDensity(currentTemperature)
    if ( ( logMHRatio > 0 ) || ( log(rdg.nextUniform(0,1)) < logMHRatio ) ) (proposedAIBDPrior,1,1) else (aibdPrior,0,1)
  }

  def updateSamplingModel(featureAllocation: FeatureAllocation, featureAllocationPrior: FeatureAllocationDistribution, lglfm: LinearGaussianLatentFeatureModel, rdg: RandomDataGenerator, maxStandardDeviationX: Double, maxStandardDeviationW: Double, sdProposedStandardDeviationX: Double, sdProposedStandardDeviationW: Double, corProposedSdXSdW: Double): (LinearGaussianLatentFeatureModel, Int, Int) = {
    val m1 = lglfm.standardDeviationX
    val m2 = lglfm.standardDeviationW
    val proposedStandardDeviationX = rdg.nextGaussian(m1, sdProposedStandardDeviationX)
    if ( ( proposedStandardDeviationX < 0 ) || ( proposedStandardDeviationX > maxStandardDeviationX ) ) return (lglfm,0,1)
    val proposedStandardDeviationW = rdg.nextGaussian(m2 + (sdProposedStandardDeviationW/sdProposedStandardDeviationX) * corProposedSdXSdW * (proposedStandardDeviationX - m1), sqrt((1 - corProposedSdXSdW*corProposedSdXSdW)*sdProposedStandardDeviationW*sdProposedStandardDeviationW))
    if ( ( proposedStandardDeviationW < 0 ) || ( proposedStandardDeviationW > maxStandardDeviationW ) ) return (lglfm,0,1)
    val proposedLGLFM = lglfm.updateStandardDeviations(proposedStandardDeviationX, proposedStandardDeviationW)
    val logMHRatio = proposedLGLFM.logLikelihood(featureAllocation) - lglfm.logLikelihood(featureAllocation)
    if ( ( logMHRatio > 0 ) || ( log(rdg.nextUniform(0,1)) < logMHRatio ) ) (proposedLGLFM,1,1) else (lglfm,0,1)
  }

  def updateFeatureAllocation(featureAllocation: FeatureAllocation, featureAllocationPrior: FeatureAllocationDistribution, lglfm: LinearGaussianLatentFeatureModel, rdg: RandomDataGenerator, rankOneUpdates: Boolean, newFeaturesTruncationDivisor: Double = 1000): (FeatureAllocation, Int, Int) = {
    val nItems = lglfm.N
    val logNewFeaturesTruncationDivisor = log(newFeaturesTruncationDivisor)
    var state = featureAllocation
    var accepts = 0
    var attempts = 0
    for (i <- 0 until nItems) {
      //val (stateNew, n, d) = updateFeatureAllocationOfExistingByEnumeration(i, state, featureAllocationPrior, lglfm, rdg, true, rankOneUpdates)
      val (stateNew, n, d) = updateFeatureAllocationOfExistingOneByOne(i, state, featureAllocationPrior, lglfm, rdg, rankOneUpdates)
      state = stateNew
      accepts += n
      attempts += d
      state = updateFeatureAllocationOfSingletons(i, state, featureAllocationPrior, lglfm, rdg, logNewFeaturesTruncationDivisor)
    }
    (state, accepts, attempts)
  }

  private def logPosterior0(i: Int, fa: FeatureAllocation, likelihoodComponentsOption: Option[LikelihoodComponents], featureAllocationPrior: FeatureAllocationDistribution, lglfm: LinearGaussianLatentFeatureModel): Double = {
    featureAllocationPrior.logProbability(i, fa) + { if ( likelihoodComponentsOption.isDefined ) lglfm.logLikelihood(likelihoodComponentsOption.get) else lglfm.logLikelihood(fa) }
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
      val newCumState = weights.head._1.add(i)
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
        if (proposals.isEmpty) Vector[(FeatureAllocation, Double)]()
        else {
          val lc1 = lglfm.computeLikelihoodComponents(proposals.head)
          val lc2 = lglfm.deallocateFeaturesFor(i, lc1)
          if (parallel) ParVector(proposals:_*).map(logPosterior2(i, _, featureAllocationPrior, lglfm, lc2)).toVector else proposals.map(logPosterior2(i, _, featureAllocationPrior, lglfm, lc2))
        }
      } else {
        if (parallel) ParVector(proposals:_*).map(logPosterior1(i, _, featureAllocationPrior, lglfm)).toVector else proposals.map(logPosterior1(i, _, featureAllocationPrior, lglfm))
      }
    }
    (rdg.nextItem(logWeights, onLogScale = true)._1, 1, 1)
  }

  def updateFeatureAllocationOfExistingOneByOne(i: Int, featureAllocation: FeatureAllocation, featureAllocationPrior: FeatureAllocationDistribution, lglfm: LinearGaussianLatentFeatureModel, rdg: RandomDataGenerator, rankOneUpdates: Boolean): (FeatureAllocation, Int, Int) = {
    if ( featureAllocation.nFeatures == 0 ) return (featureAllocation,0,0)
    var accepts = 0
    var attempts = 0
    var state = featureAllocation
    var stateMap = state.asCountMap
    var stateLikelihoodComponentsOption = if ( rankOneUpdates ) Some(lglfm.computeLikelihoodComponents(state)) else None
    var stateLogPosterior = logPosterior0(i, state, stateLikelihoodComponentsOption, featureAllocationPrior, lglfm)
    for ( j <- rdg.nextPermutation(state.nFeatures, state.nFeatures); if ! state.isSingleton(i,j) ) {
      val proposal = if ( state.features(j).contains(i) ) state.remove(i,j) else state.add(i,j)
      val proposalMap = proposal.asCountMap
      val proposalLikelihoodComponentsOption = if ( rankOneUpdates ) Some(lglfm.reallocateFeaturesFor(i, stateLikelihoodComponentsOption.get, proposal.matrixRow(i))) else None
      val proposalLogPosterior = logPosterior0(i, proposal, proposalLikelihoodComponentsOption, featureAllocationPrior, lglfm)
      val logMHRatio = proposalLogPosterior - stateLogPosterior - logOnInt(stateMap((state.features(j),state.sizes(j)))) + logOnInt(proposalMap((proposal.features(j),proposal.sizes(j))))
      attempts += 1
      if ( ( logMHRatio >= 0 ) || ( log(rdg.nextUniform(0.0,1.0)) < logMHRatio ) ) {
        accepts += 1
        state = proposal
        stateMap = proposalMap
        stateLikelihoodComponentsOption = proposalLikelihoodComponentsOption
        stateLogPosterior = proposalLogPosterior
      }
    }
    (state, accepts, attempts)
  }

}

