package org.ddahl.aibd

import org.ddahl.aibd.parameter.{NullParameterDistribution, ParameterDistribution}
import org.ddahl.commonsmath._
import org.apache.commons.math3.util.CombinatoricsUtils.{binomialCoefficientLog, factorialLog}
import org.apache.commons.math3.util.FastMath.log
import org.apache.commons.math3.random.RandomDataGenerator

class AttractionIndianBuffetDistribution[A] protected (val mass: Double, override val permutation: Permutation, val similarity: Similarity, val parameterDistribution: ParameterDistribution[A]) extends FeatureAllocationDistribution[A] {

  val nItems = permutation.nItems
  if (similarity.nItems != nItems) throw new RuntimeException("Inconsistent number of items.")
  private val faEmpty = FeatureAllocation.empty[A](nItems)

  def dropParameter = new AttractionIndianBuffetDistribution[Null](mass,permutation,similarity,NullParameterDistribution)

  def sample(rdg: RandomDataGenerator): FeatureAllocation[A] = {
    var fa = faEmpty
    var index = 0
    while (index < nItems) {
      fa = sampleByIndex(index, fa, rdg)
      index += 1
    }
    fa
  }

  def sampleByIndex(index: Int, fa: FeatureAllocation[A], rdg: RandomDataGenerator): FeatureAllocation[A] = {
    var _fa = fa
    val i = permutation(index)
    val const = if (similarity.isUniform) 1.0
    else {
      index / (0 until index).map(jj => similarity(i, permutation(jj))).sum
    }
    fa.foreach { f =>
      val numerator: Double = if (similarity.isUniform) f.size.toDouble
      else const * f.set.foldLeft(0.0)((sum, j) => sum + similarity(i, j))
      val p = numerator / (index + 1)
      if (rdg.nextUniform(0.0, 1.0) <= p) {
        _fa = _fa.add(i, f)
      }
    }
    repeat(rdg.nextPoisson(mass / (index + 1))) {
      _fa = _fa.add(i, parameterDistribution.sample(rdg))
    }
    _fa
  }

  def logDensity(fa: FeatureAllocation[A], parallel: Boolean): Double = logDensityStartingFromIndex(0, fa, parallel)

  override def logDensityStartingFromIndex(index: Int, fa: FeatureAllocation[A], parallel: Boolean): Double = {
    var unclaimedFeatures = fa.toVector
    var claimedFeatures = Vector[Feature[A]]()
    val cache = Array.tabulate(nItems)(_index => {
      def strip(features: Vector[Feature[A]]) = {
        features.map(f => {
          val items = f.set.filter(j => permutation.inverse(j) < _index)
          Feature(f.parameter, items)
        })
      }
      val i = permutation(_index)
      val (old, not) = claimedFeatures.partition(_.contains(i))
      val (fresh, stillUnclaimedFeatures) = unclaimedFeatures.partition(_.contains(i))
      claimedFeatures = fresh ++ claimedFeatures
      unclaimedFeatures = stillUnclaimedFeatures
      if (_index < index) null
      else (fresh.map(f => Feature(f.parameter, i)), strip(old), strip(not))
    })
    val range = index until nItems
    if (parallel) range.par.map(index => logDensityByIndex(index, fa, cache(index)._1, cache(index)._2, cache(index)._3)).sum
    else range.map(index => logDensityByIndex(index, fa, cache(index)._1, cache(index)._2, cache(index)._3)).sum
  }

  private def logDensityByIndex(index: Int, fa: FeatureAllocation[A], newFeatures: Vector[Feature[A]], yesFeatures: Vector[Feature[A]], noFeatures: Vector[Feature[A]]): Double = {
    val i = permutation(index)
    var logD = 0.0
    val const = if (similarity.isUniform) 1.0
    else index / (0 until index).foldLeft(0.0)((sum, jj) => sum + similarity(i, permutation(jj)))
    val set = yesFeatures.toSet | noFeatures.toSet
    set.foreach(f => {
      val numerator = if (similarity.isUniform) f.size.toDouble
      else const * f.set.foldLeft(0.0)((sum, j) => sum + similarity(i, j))
      val p = numerator / (index + 1)
      val kYes = yesFeatures.count(_ == f)
      val kNo = noFeatures.count(_ == f)
      logD += binomialCoefficientLog(kYes + kNo, kYes) + kYes * log(p) + kNo * log(1 - p)
    })
    val nNewFeatures = newFeatures.size
    val alpha = mass / (index + 1)
    logD += nNewFeatures * log(alpha) - alpha - factorialLog(nNewFeatures)
    val xx = newFeatures.groupBy(identity).mapValues(_.size).map(_._2)
    logD += factorialLog(nNewFeatures) - xx.map(factorialLog).sum
    newFeatures.foreach(f => {
      logD += parameterDistribution.logDensity(f.parameter)
    })
    logD
  }

  override def update(nScans: Int, fa: FeatureAllocation[A], rdg: RandomDataGenerator, parallel: Boolean): (FeatureAllocationDistribution[A], Int, Int) = {
    updatePermutation(nScans, fa, rdg, parallel)
  }

  def updatePermutation(nScans: Int, fa: FeatureAllocation[A], rdg: RandomDataGenerator, parallel: Boolean): (AttractionIndianBuffetDistribution[A], Int, Int) = {
    var current = this
    var attempts = 0
    var acceptances = 0
    if ( permutation.nPerShuffle > 0 ) {
      repeat(nScans) {
        val proposal = AttractionIndianBuffetDistribution(mass, permutation.shuffle(rdg), similarity, parameterDistribution)
        val logMHRatio = proposal.logDensity(fa, parallel) - logDensity(fa, parallel)
        attempts += 1
        if ((logMHRatio >= 0.0) || (log(rdg.nextUniform(0.0, 1.0)) < logMHRatio)) {
          acceptances += 1
          current = proposal
        }
      }
    }
    (current, acceptances, attempts)
  }

  def toEmpericalDistributionViaMCMC(rdg: RandomDataGenerator, nSamples: Int, nCores: Int = Runtime.getRuntime.availableProcessors) = {
    val dist = this.dropParameter
    val nSamplesPerCore = (nSamples - 1) / nCores + 1
    val efad = rdg.nextRandomDataGenerators(nCores).map(rdg => {
      var efad2 = EmpericalFeatureAllocationDistribution()
      var fa = dist.sample(rdg)
      repeat(nSamplesPerCore) {
        // fa = MCMCSamplers.updateFeatureAllocationGibbs(1, fa, dist, (i: Int, fa: FeatureAllocation[Null]) => 0.0, rdg, false)
        fa = MCMCSamplers.updateFeatureAllocationNeighborhoods(1, (x: Int) => x, fa, dist, (i: Int, fa: FeatureAllocation[Null]) => 0.0, rdg, false)._1
        // fa = MCMCSamplers.updateFeatureAllocationInOrOut(1, fa, dist, (i: Int, fa: FeatureAllocation[Null]) => 0.0, rdg, false)._1
        fa = MCMCSamplers.updateFeatureAllocationSingletons(1, fa, dist, (i: Int, fa: FeatureAllocation[Null]) => 0.0, rdg)._1
        efad2 = efad2.tally(fa)
      }
      efad2
    }).fold(EmpericalFeatureAllocationDistribution())((x, y) => {
      x ++ y
    })
    efad
  }

}

object AttractionIndianBuffetDistribution {

  def apply(mass: Double, permutation: Permutation, similarity: Similarity) = {
    new AttractionIndianBuffetDistribution(mass, permutation, similarity, NullParameterDistribution)
  }

  def apply[A](mass: Double, permutation: Permutation, similarity: Similarity, parameterDistribution: ParameterDistribution[A]) = {
    new AttractionIndianBuffetDistribution(mass, permutation, similarity, parameterDistribution)
  }

}
