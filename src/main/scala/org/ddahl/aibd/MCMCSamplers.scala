package org.ddahl.aibd

import org.ddahl.aibd.Utils.logOnInt
import org.ddahl.commonsmath._
import org.apache.commons.math3.random.RandomDataGenerator
import org.apache.commons.math3.distribution.EnumeratedDistribution
import org.apache.commons.math3.util.CombinatoricsUtils.factorialLog
import org.apache.commons.math3.util.Pair
import org.apache.commons.math3.util.FastMath.{log, exp}
import scala.collection.JavaConverters._

object MCMCSamplers {

  // Original method of Griffiths & Ghahramani (2005), but better explained in Section 3 of "Accelerated Sampling for the Indian Buffet Process" by Doshi-Velez, Ghahramani (2009)
  def updateFeatureAllocationIBP[A](nScans: Int, fa: FeatureAllocation[A], ibp: IndianBuffetProcess[A], logLikelihood: (Int, FeatureAllocation[A]) => Double, maxNewFeatures: Int, rdg: RandomDataGenerator): FeatureAllocation[A] = {
    import ibp.nItems
    val logRate = log(ibp.mass / nItems)
    var faCurrent = fa
    repeat(nScans) {
      for (i <- 0 until nItems) {
        // Work on existing features
        for (f <- faCurrent) {
          if (f.contains(i) && f.size == 1) faCurrent = faCurrent.remove(f)
          else {
            val (m, fa0, fa1) = if (f.contains(i)) {
              val fM = f.remove(i)
              val faM = faCurrent.replace(f, fM)
              (fM.size, faM, faCurrent)
            } else {
              val fP = f.add(i)
              val faP = faCurrent.replace(f, fP)
              (f.size, faCurrent, faP)
            }
            val logWeight0 = logOnInt(m) + logLikelihood(i, fa0)
            val logWeight1 = logOnInt(m) + logLikelihood(i, fa1)
            val numerator = exp(logWeight1)
            val denominator = exp(logWeight0) + numerator
            if (rdg.nextUniform(0.0, denominator) <= numerator) faCurrent = fa1
            else faCurrent = fa0
          }
        }
        // Work on new features
        val candidates = Array.fill(maxNewFeatures)(Feature(ibp.parameterDistribution.sample(rdg), i))
          .scanLeft((faCurrent, exp(0.0 + logLikelihood(i, faCurrent)), 0.0, 1)) { (state, f) =>
            val newFA = state._1.add(f)
            val logIncrement = state._3 + logRate - logOnInt(state._4)
            val weight = exp(logIncrement + logLikelihood(i, newFA))
            (newFA, weight, logIncrement, state._4 + 1)
          }
        val weights = candidates.map(fa => new Pair[FeatureAllocation[A], java.lang.Double](fa._1, fa._2)).toList.asJava
        val dist = new EnumeratedDistribution(rdg.getRandomGenerator, weights)
        faCurrent = dist.sample
      }
    }
    faCurrent
  }

  def updateFeatureAllocationInOrOut[A](nScans: Int, fa: FeatureAllocation[A], faDistribution: FeatureAllocationDistribution[A], logLikelihood: (Int, FeatureAllocation[A]) => Double, rdg: RandomDataGenerator, parallel: Boolean): (FeatureAllocation[A], Int, Int) = {
    import faDistribution.{nItems, permutation}
    var faCurrent = fa
    var attempts = 0
    var acceptances = 0
    repeat(nScans) {
      for (index <- 0 until nItems) {
        val i = permutation(index)
        var faCurrentMH = faDistribution.logDensityStartingFromIndex(index, faCurrent, parallel) + logLikelihood(i, faCurrent)
        for (f <- faCurrent.filterNot(f => f.contains(i) && (f.size == 1))) {
          val faProposal = if (f.contains(i)) faCurrent.remove(i, f)
          else faCurrent.add(i, f)
          val faProposalMH = faDistribution.logDensityStartingFromIndex(index, faProposal, parallel) + logLikelihood(i, faProposal)
          attempts += 1
          if ((faProposalMH >= faCurrentMH) || (log(rdg.nextUniform(0.0, 1.0)) < faProposalMH - faCurrentMH)) {
            acceptances += 1
            faCurrent = faProposal
            faCurrentMH = faProposalMH
          }
        }
      }
    }
    (faCurrent, acceptances, attempts)
  }

  private def neighbors[A](i: Int, fa: FeatureAllocation[A]): Vector[FeatureAllocation[A]] = {
    def engine(fa: FeatureAllocation[A], features: List[Feature[A]]): List[FeatureAllocation[A]] = {
      if (features.isEmpty) List(fa)
      else {
        val fa2 = fa.add(i, features.head)
        engine(fa, features.tail) ++ engine(fa2, features.tail)
      }
    }

    engine(fa, fa.toList).toSet.toVector
  }

  def updateFeatureAllocationNeighborhoods[A](nScans: Int, nRepeatsCalculation: (Int) => Int, fa: FeatureAllocation[A], faDistribution: FeatureAllocationDistribution[A], logLikelihood: (Int, FeatureAllocation[A]) => Double, rdg: RandomDataGenerator, parallel: Boolean): (FeatureAllocation[A], Int, Int, Int, Int) = {
    import faDistribution.{nItems, permutation}
    var faCurrent = fa
    var attempts = 0
    var acceptances = 0
    var hitAttempts = 0
    var hitAcceptances = 0
    repeat(nScans) {
      for (index <- 0 until nItems) {
        val i = permutation(index)
        val hashMap = new scala.collection.mutable.HashMap[FeatureAllocation[A], Double]()
        var faCurrentMH = hashMap.getOrElseUpdate(faCurrent, faDistribution.logDensityStartingFromIndex(index, faCurrent, parallel) + logLikelihood(i, faCurrent))
        val newFeatures = faCurrent.filter(f => f.contains(i) && (f.size == 1))
        val faCurrentWithoutI = faCurrent.remove(i)
        val faCandidates = neighbors(i, faCurrentWithoutI).map(_.add(newFeatures))
        val nRepeats = nRepeatsCalculation(faCurrentWithoutI.nFeatures)
        hitAttempts += 1
        var hit = false
        repeat(nRepeats) {
          val faProposal = faCandidates(rdg.nextInt(0, faCandidates.size - 1))
          val faProposalMH = hashMap.getOrElseUpdate(faProposal, faDistribution.logDensityStartingFromIndex(index, faProposal, parallel) + logLikelihood(i, faProposal))
          attempts += 1
          if ((faProposalMH >= faCurrentMH) || (log(rdg.nextUniform(0.0, 1.0)) < faProposalMH - faCurrentMH)) {
            if (!hit) {
              hitAcceptances += 1
              hit = true
            }
            acceptances += 1
            faCurrent = faProposal
            faCurrentMH = faProposalMH
          }
        }
      }
    }
    (faCurrent, acceptances, attempts, hitAcceptances, hitAttempts)
  }

  def updateFeatureAllocationGibbs[A](fa: FeatureAllocation[A], faDistribution: FeatureAllocationDistribution[A], logLikelihood: FeatureAllocation[A] => Double, nSamples: Int, thin: Int, rdg: RandomDataGenerator, parallel: Boolean): Seq[FeatureAllocation[A]] = {
    val nIterations = thin*nSamples
    var state = fa
    var results = List[FeatureAllocation[A]]()
    var b = 1
    val logLike = (i: Int, fa: FeatureAllocation[A]) => logLikelihood(fa)
    while (b <= nIterations) {
      state = updateFeatureAllocationGibbs(thin, state, faDistribution, logLike, rdg, parallel)
      if (b % thin == 0) results = state :: results
      b += 1
    }
    results.reverse
  }

  def updateFeatureAllocationGibbs[A](nScans: Int, fa: FeatureAllocation[A], faDistribution: FeatureAllocationDistribution[A], logLikelihood: (Int, FeatureAllocation[A]) => Double, rdg: RandomDataGenerator, parallel: Boolean): FeatureAllocation[A] = {
    import faDistribution.{nItems, permutation}
    var faCurrent = fa
    repeat(nScans) {
      for (index <- 0 until nItems) {
        val i = permutation(index)
        val newFeatures = faCurrent.filter(f => f.contains(i) && (f.size == 1))
        val faCandidates = neighbors(i, faCurrent.remove(i)).map(_.add(newFeatures))
        val iterator = if (parallel) faCandidates.par else faCandidates
        val weights = iterator.map(fa => {
          new Pair[FeatureAllocation[A], java.lang.Double](fa, exp(faDistribution.logDensityStartingFromIndex(index, fa, parallel) + logLikelihood(i, fa)))
        }).toList.asJava
        val dist = new EnumeratedDistribution(rdg.getRandomGenerator, weights)
        faCurrent = dist.sample
      }
    }
    faCurrent
  }

  // Demonstrates that the MH ratio is 1.
  def updateFeatureAllocationGibbsSlow[A](nScans: Int, fa: FeatureAllocation[A], faDistribution: FeatureAllocationDistribution[A], logLikelihood: (Int, FeatureAllocation[A]) => Double, rdg: RandomDataGenerator, parallel: Boolean): FeatureAllocation[A] = {
    import faDistribution.{mass, nItems, permutation}
    var faCurrent = fa
    var attempts = 0
    var acceptances = 0
    repeat(nScans) {
      for (index <- 0 until nItems) {
        val i = permutation(index)
        val alpha = mass / (index + 1)
        val neighborhood = neighbors(i, faCurrent.remove(i))
        val newFeatures = faCurrent.filter(f => f.contains(i) && (f.size == 1))

        def engine(target: Option[FeatureAllocation[A]]): (FeatureAllocation[A], Double, Double) = {
          val faCandidates = neighborhood.map(_.add(newFeatures))
          val iterator = if (parallel) faCandidates.par else faCandidates
          val logWeights = iterator.map(fa => {
            fa -> (faDistribution.logDensityStartingFromIndex(index, fa, false) + logLikelihood(i, fa))
          })
          val weights = logWeights.map(x => x._1 -> exp(x._2))
          val sumOfWeights = weights.toList.map(_._2).sum
          val probabilities = weights.map(x => x._1 -> x._2 / sumOfWeights)
          val faProposal = target.getOrElse({
            val probabilitiesForJava = probabilities.toList.map(x => new Pair[FeatureAllocation[A], java.lang.Double](x._1, x._2)).asJava
            val dist = new EnumeratedDistribution(rdg.getRandomGenerator, probabilitiesForJava)
            dist.sample
          })
          val newFeaturesMap = newFeatures.groupBy(identity).mapValues(_.size)
          val logMHRatioContributionNewFeatures = newFeaturesMap.map(x => {
            -factorialLog(x._2) + x._2 * faDistribution.parameterDistribution.logDensity(x._1.parameter)
          }).sum + factorialLog(newFeatures.size)
          val logMHRatioContribution = newFeatures.size * log(alpha) - alpha - factorialLog(newFeatures.size) +
            logMHRatioContributionNewFeatures +
            log(probabilities.find(x => x._1 == faProposal).get._2)
          (faProposal, logWeights.find(x => x._1 == faProposal).get._2, logMHRatioContribution)
        }

        val proposalTuple = engine(None)
        val currentTuple = engine(Some(faCurrent))
        val logMHRatio = proposalTuple._2 - currentTuple._2 - proposalTuple._3 + currentTuple._3
        assert(-0.0000001 <= logMHRatio && logMHRatio <= 0.0000001)
        attempts += 1
        if ((logMHRatio >= 0.0) || (log(rdg.nextUniform(0.0, 1.0)) < logMHRatio)) {
          acceptances += 1
          faCurrent = proposalTuple._1
        }
      }
    }
    faCurrent
  }

  def updateFeatureAllocationSingletons[A](nScans: Int, fa: FeatureAllocation[A], faDistribution: FeatureAllocationDistribution[A], logLikelihood: (Int, FeatureAllocation[A]) => Double, rdg: RandomDataGenerator): (FeatureAllocation[A], Int, Int) = {
    import faDistribution.{mass, nItems, permutation}
    var faCurrent = fa
    var attempts = 0
    var acceptances = 0
    repeat(nScans) {
      for (index <- 0 until nItems) {
        val i = permutation(index)

        def engine(fa: FeatureAllocation[A], features: Vector[Feature[A]]): Double = {
          val logDensity = faDistribution.logDensityStartingFromIndex(index, fa, false) + logLikelihood(i, fa)
          val featuresMap = features.groupBy(identity).mapValues(_.size)
          val logMHRatioContribution = featuresMap.map(x => {
            -factorialLog(x._2) + x._2 * faDistribution.parameterDistribution.logDensity(x._1.parameter)
          }).sum + factorialLog(features.size)
          logDensity - logMHRatioContribution
        }

        val set = Set(i)
        val alpha = mass / (index + 1)
        val oldFeatures = faCurrent.filter(f => f.contains(i) && (f.size == 1))
        val newFeatures = Vector.fill(rdg.nextPoisson(alpha).toInt) {
          Feature(faDistribution.parameterDistribution.sample(rdg), set)
        }
        val faProposal = faCurrent.remove(oldFeatures).add(newFeatures)
        val logMHRatioContributionSizes = if (oldFeatures.size == newFeatures.size) 0.0 else (oldFeatures.size - newFeatures.size) * log(alpha) - factorialLog(oldFeatures.size) + factorialLog(newFeatures.size)
        val logMHRatio = engine(faProposal, newFeatures) - engine(faCurrent, oldFeatures.toVector) + logMHRatioContributionSizes
        attempts += 1
        if ((logMHRatio >= 0.0) || (log(rdg.nextUniform(0.0, 1.0)) < logMHRatio)) {
          acceptances += 1
          faCurrent = faProposal
        }
      }
    }
    (faCurrent, acceptances, attempts)
  }

  //-----

  def updateFeatureAllocationGibbsPriorOnly(fa: FeatureAllocation[Null], ibp: IndianBuffetProcess[Null], nSamples: Int, thin: Int, rdg: RandomDataGenerator, newFeaturesTruncationDivisor: Double = 1000): Seq[FeatureAllocation[Null]] = {
    val nItems = fa.nItems
    val rate = ibp.mass/nItems
    var results = List[FeatureAllocation[Null]]()
    var state = fa
    for (b <- 1 to nSamples) {
      for (i <- 0 until nItems) {
        for (feature <- state.features) {
          val m = feature.size - (if (feature.contains(i)) 1 else 0)
          state = if (m == 0) state.remove(feature)
          else if (rdg.nextUniform(0, 1) < m.toDouble / nItems) {
            if (feature.contains(i)) state else state.add(i, feature)
          } else {
            if (feature.contains(i)) state.remove(i, feature) else state
          }
        }
        @scala.annotation.tailrec
        def engine(weights: List[(Int,Double)], k: Int, max: Double): List[(Int,Double)] = {
          val newWeight = weights.head._2 * rate/k
          val expanded = (k,newWeight) :: weights
          if ( newWeight < max/newFeaturesTruncationDivisor ) expanded
          else engine(expanded, k+1, if ( newWeight > max ) newWeight else max)
        }
        val weights = engine(List((0,rate)),1,0.0)
        val nNewFeatures = rdg.nextItem(weights.toIndexedSeq)._1
        val featureWithOnlyI = Feature(null, i)
        for (r <- 0 until nNewFeatures) state = state.add(featureWithOnlyI)
      }
      if (b % thin == 0) results = state :: results
    }
    results.reverse
  }

  def updateFeatureAllocationGibbsWithLikelihood(fa: FeatureAllocation[Null], ibp: IndianBuffetProcess[Null], logLikelihood: FeatureAllocation[Null] => Double, nSamples: Int, thin: Int, rdg: RandomDataGenerator, newFeaturesTruncationDivisor: Double = 1000): Seq[FeatureAllocation[Null]] = {
    val nItems = fa.nItems
    val rate = ibp.mass/nItems
    val nIterations = thin*nSamples
    val logRate = log(rate)
    val logNewFeaturesTruncationDivisor = log(newFeaturesTruncationDivisor)
    var state = fa
    var results = List[FeatureAllocation[Null]]()
    var b = 1
    while (b <= nIterations) {
      for (i <- 0 until nItems) {
        for (feature <- state.features) {
          val contains = feature.contains(i)
          val m = feature.size - (if (contains) 1 else 0)
          state = if (m == 0) state.remove(feature) else {
            val (state0, state1) = if ( contains ) (state.remove(i,feature),state) else (state,state.add(i,feature))
            val weight0 = exp(logLikelihood(state0) + logOnInt(nItems - m))
            val weight1 = exp(logLikelihood(state1) + logOnInt(m))
            if (rdg.nextUniform(0, 1) < weight1 / ( weight0 + weight1 )) state1 else state0
          }
        }
        val featureWithOnlyI = Feature(null, i)
        @scala.annotation.tailrec
        def engine(weights: List[(FeatureAllocation[Null],Double)], cumProduct: Double, k: Int, max: Double): List[(FeatureAllocation[Null],Double)] = {
          val newCumProduct = cumProduct + logRate - logOnInt(k)
          val newCumState = weights.head._1.add(featureWithOnlyI)
          val newWeight = newCumProduct + logLikelihood(newCumState)
          val expanded = (newCumState,newWeight) :: weights
          if ( newWeight < max - logNewFeaturesTruncationDivisor ) expanded
          else engine(expanded, newCumProduct, k+1, if ( newWeight > max ) newWeight else max)
        }
        val weights = engine((state,logRate + logLikelihood(state)) :: Nil,logRate,1,Double.NegativeInfinity).toIndexedSeq
        state = rdg.nextItem(weights, onLogScale = true)._1
      }
      if (b % thin == 0) results = state :: results
      b += 1
    }
    results.reverse
  }

  def updateFeatureAllocationIndependence(fa: FeatureAllocation[Null], ibp: IndianBuffetProcess[Null], logLikelihood: FeatureAllocation[Null] => Double, nSamples: Int, thin: Int, rdg: RandomDataGenerator): Seq[FeatureAllocation[Null]] = {
    val nItems = fa.nItems
    val nIterations = thin*nSamples
    var state = fa
    var stateLogLikelihood = logLikelihood(state)
    var results = List[FeatureAllocation[Null]]()
    val monitor = MCMCAcceptanceMonitor1()
    var b = 1
    while (b <= nIterations) {
      for (i <- 0 until nItems) {
        val proposal = ibp.sample(i,state,rdg)
        val proposalLogLikelihood = logLikelihood(proposal)
        state = monitor {
          if ( log(rdg.nextUniform(0.0,1.0)) < proposalLogLikelihood - stateLogLikelihood ) {
            stateLogLikelihood = proposalLogLikelihood
            (proposal,1,1)
          } else {
            (state,0,1)
          }
        }
      }
      if (b % thin == 0) results = state :: results
      b += 1
    }
    //println("Monitor: "+monitor.rate)
    results.reverse
  }

  def allPossibleConfigurationsAmongExistingKeepingSingletons[A](i: Int, fa: FeatureAllocation[A]): List[FeatureAllocation[A]] = {
    val singletonFeatures = fa.featuresOf(i).filter(_.size == 1)
    val faWithoutI = fa.remove(i)
    var faList = List[FeatureAllocation[A]]()
    def engine(availableFeatures: List[Feature[A]], occupiedFeatures: List[Feature[A]]): Unit = {
      if ( availableFeatures.isEmpty ) {
        faList = occupiedFeatures.foldLeft(faWithoutI)((cumFA,f) => cumFA.add(i,f)) :: faList
      } else {
        engine(availableFeatures.tail,                           occupiedFeatures)
        engine(availableFeatures.tail, availableFeatures.head :: occupiedFeatures)
      }
    }
    //engine(faWithoutI, faWithoutI.toList).map(fa => fa.add(singletonFeatures).leftOrderedForm).toSet.toList  // Remove duplicates
    engine(faWithoutI.toList, Nil)
    faList.map(_.add(singletonFeatures).leftOrderedForm).toSet.toList
  }

  def updateFeatureAllocationViaNeighborhoods(fa: FeatureAllocation[Null], ibp: IndianBuffetProcess[Null], logLikelihood: FeatureAllocation[Null] => Double, nSamples: Int, thin: Int, rdg: RandomDataGenerator, newFeaturesTruncationDivisor: Double = 1000): Seq[FeatureAllocation[Null]] = {
    val nItems = fa.nItems
    val nIterations = thin*nSamples
    val logNewFeaturesTruncationDivisor = log(newFeaturesTruncationDivisor)
    var state = fa
    var results = List[FeatureAllocation[Null]]()
    var b = 1
    while (b <= nIterations) {
      for (i <- 0 until nItems) {
        val proposals = allPossibleConfigurationsAmongExistingKeepingSingletons(i, state)
        val logWeights = proposals.toIndexedSeq.map(fa => (fa, ibp.logDensity(fa) + logLikelihood(fa)))
        state = rdg.nextItem(logWeights, onLogScale = true)._1
        state = FeatureAllocation(state.nItems, state.filterNot(f => (f.size == 1 ) && (f.contains(i))).toVector)
        val featureWithOnlyI = Feature(null, i)
        @scala.annotation.tailrec
        def engine(weights: List[(FeatureAllocation[Null],Double)], max: Double): List[(FeatureAllocation[Null],Double)] = {
          val newCumState = weights.head._1.add(featureWithOnlyI)
          val newLogWeight = ibp.logDensity(newCumState) + logLikelihood(newCumState)
          val expanded = (newCumState,newLogWeight) :: weights
          if ( newLogWeight < max - logNewFeaturesTruncationDivisor ) expanded
          else engine(expanded, if ( newLogWeight > max ) newLogWeight else max)
        }
        val weights = engine((state,ibp.logDensity(state) + logLikelihood(state)) :: Nil,Double.NegativeInfinity).toIndexedSeq
        state = rdg.nextItem(weights, onLogScale = true)._1
      }
      if (b % thin == 0) results = state :: results
      b += 1
    }
    results.reverse
  }

}

