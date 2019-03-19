package org.ddahl.aibd.model.lineargaussian

import org.ddahl.matrix._
import org.apache.commons.math3.random.RandomDataGenerator
import org.apache.commons.math3.util.FastMath.log
import org.ddahl.commonsmath.RandomDataGeneratorImprovements
import org.ddahl.aibd.{Feature, FeatureAllocation, IndianBuffetProcess}

object PosteriorSimulation {

  def fa2Z(fa: FeatureAllocation[Null]): Matrix = wrap(fa.leftOrderedForm.toMatrix.map(_.map(_.toDouble)))

  def Z2fa(Z: Matrix, nItems: Int): FeatureAllocation[Null] = {
    if ( Z == null ) FeatureAllocation.empty[Null](nItems)
    else {
      FeatureAllocation(Z.rows, (0 until Z.cols).map { i =>
        val col = Z(::, i)
        Feature(null, col.zipWithIndex.filter(x => x._1 != 0.0).map(_._2))
      }.toVector).leftOrderedForm
    }
  }

  def updateFeatureAllocationViaNeighborhoods(Z: Matrix, ibp: IndianBuffetProcess[Null], lglfm: LinearGaussianLatentFeatureModel, nSamples: Int, thin: Int, rdg: RandomDataGenerator, newFeaturesTruncationDivisor: Double = 1000): Seq[FeatureAllocation[Null]] = {
    val nItems = lglfm.N
    var state = Z
    val nIterations = thin*nSamples
    val logNewFeaturesTruncationDivisor = log(newFeaturesTruncationDivisor)
    var results = List[FeatureAllocation[Null]]()
    var b = 1
    while (b <= nIterations) {
      for (i <- 0 until nItems) {
        val (singletons, existing) = FeatureAllocationUtilities.partitionBySingletonsOf(i, state)
        val proposals = if (existing == null) Array(singletons)
        else FeatureAllocationUtilities.enumerateCombinationsFor(i, existing).map(_ | singletons)
        val proposals2 = proposals.map(Z => Z2fa(Z,nItems).leftOrderedForm).toSet.toArray.map(fa2Z)
        val logWeights = proposals2.map(Z => (Z, ibp.logDensity(Z2fa(Z,nItems)) + lglfm.logLikelihood(Z)))
        state = rdg.nextItem(logWeights, onLogScale = true)._1
        state = FeatureAllocationUtilities.partitionBySingletonsOf(i, state)._2
        val featureWithOnlyI = Array.ofDim[Double](nItems)
        featureWithOnlyI(i) = 1.0
        @scala.annotation.tailrec
        def engine(weights: List[(Matrix,Double)], max: Double): List[(Matrix,Double)] = {
          val newCumState = weights.head._1 | featureWithOnlyI
          val newLogWeight = ibp.logDensity(Z2fa(newCumState,nItems)) + lglfm.logLikelihood(newCumState)
          val expanded = (newCumState,newLogWeight) :: weights
          if ( newLogWeight < max - logNewFeaturesTruncationDivisor ) expanded
          else engine(expanded, if ( newLogWeight > max ) newLogWeight else max)
        }
        val setup = (state,ibp.logDensity(Z2fa(state,nItems)) + lglfm.logLikelihood(state)) :: Nil
        val weights = engine(setup,Double.NegativeInfinity).toIndexedSeq
        state = rdg.nextItem(weights, onLogScale = true)._1
      }
      if (b % thin == 0) results = Z2fa(state, nItems) :: results
      b += 1
    }
    results.reverse
  }

}
