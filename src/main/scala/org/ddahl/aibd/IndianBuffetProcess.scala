package org.ddahl.aibd

import org.ddahl.aibd.Utils.{harmonicNumber, logFactorial}
import org.ddahl.aibd.parameter.{NullParameterDistribution, ParameterDistribution}
import org.ddahl.commonsmath._
import org.apache.commons.math3.random.RandomDataGenerator
import org.apache.commons.math3.util.FastMath.log

class IndianBuffetProcess[A] protected(val mass: Double, val nItems: Int, val parameterDistribution: ParameterDistribution[A]) extends FeatureAllocationDistribution[A] {

  private val faEmpty = FeatureAllocation.empty[A](nItems)

  def dropParameter = new IndianBuffetProcess[Null](mass, nItems, NullParameterDistribution)

  def sample(rdg: RandomDataGenerator): FeatureAllocation[A] = {
    var fa = faEmpty
    var index = 0
    while (index < nItems) {
      val indexPlusOne = index + 1
      fa.foreach { f =>
        if (rdg.nextUniform(0.0, indexPlusOne) < f.size) fa = fa.add(index, f)
      }
      repeat(rdg.nextPoisson(mass / indexPlusOne)) {
        fa = fa.add(index, parameterDistribution.sample(rdg))
      }
      index += 1
    }
    fa
  }

  private val const1 = -mass * harmonicNumber(nItems)
  private val const2 = log(mass) - logFactorial(nItems)

  def logDensity(fa: FeatureAllocation[A], parallel: Boolean): Double = {
    var sum = const1 + fa.nFeatures * const2
    sum -= fa.foldLeft(Map.empty[Feature[Null], Int])((map, f) => {
      val ff = f.dropParameter
      map + ((ff, map.getOrElse(ff,0) + 1))
    }).values.foldLeft(0.0)((s, x) => s + logFactorial(x))
    val fs = if (parallel) fa.par else fa
    sum += fs.aggregate(0.0)((s, f) => {
      val mk = f.size
      s + logFactorial(nItems - mk) + logFactorial(mk - 1) + parameterDistribution.logDensity(f.parameter)
    }, _ + _)
    sum
  }

  def toEmpericalDistributionViaMCMC(rdg: RandomDataGenerator, nSamples: Int, nCores: Int = Runtime.getRuntime.availableProcessors) = {
    val dist = this.dropParameter
    val nSamplesPerCore = (nSamples - 1) / nCores + 1
    val efad = rdg.nextRandomDataGenerators(nCores).map(rdg => {
      var efad2 = EmpericalFeatureAllocationDistribution()
      var fa = dist.sample(rdg)
      repeat(nSamplesPerCore) {
        // fa = MCMCSamplers.updateFeatureAllocationGibbs(1, fa, dist, (i: Int, fa: FeatureAllocation[Null]) => 0.0, rdg, false)
        // fa = MCMCSamplers.updateFeatureAllocationNeighborhoods(1, (x: Int) => x, fa, dist, (i: Int, fa: FeatureAllocation[Null]) => 0.0, rdg, false)._1
        // fa = MCMCSamplers.updateFeatureAllocationInOrOut(1, fa, dist, (i: Int, fa: FeatureAllocation[Null]) => 0.0, rdg, false)._1
        fa = MCMCSamplers.updateFeatureAllocationIBP(1, fa, dist, (i: Int, fa: FeatureAllocation[Null]) => 0.0, 4, rdg)
        // fa = MCMCSamplers.updateFeatureAllocationSingletons(1, fa, dist, (i: Int, fa: FeatureAllocation[Null]) => 0.0, rdg)._1
        efad2 = efad2.tally(fa)
      }
      efad2
    }).fold(EmpericalFeatureAllocationDistribution())((x, y) => {
      x ++ y
    })
    efad
  }

}

object IndianBuffetProcess {

  def apply(mass: Double, nItems: Int) = new IndianBuffetProcess(mass, nItems, NullParameterDistribution)

  def apply[A](mass: Double, nItems: Int, parameterDistribution: ParameterDistribution[A]) = new IndianBuffetProcess(mass, nItems, parameterDistribution)

}

