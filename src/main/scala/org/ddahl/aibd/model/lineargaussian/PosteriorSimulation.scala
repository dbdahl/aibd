package org.ddahl.aibd.model.lineargaussian

import org.ddahl.matrix._
import org.apache.commons.math3.random.RandomDataGenerator
import org.apache.commons.math3.util.FastMath.log
import org.ddahl.commonsmath.RandomDataGeneratorImprovements
import org.ddahl.aibd.{Feature, FeatureAllocation, IndianBuffetProcess}
import org.ddahl.aibd.MCMCSamplers.allPossibleConfigurationsAmongExistingKeepingSingletons

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

  def updateFeatureAllocationViaNeighborhoods(fa: Matrix, ibp: IndianBuffetProcess[Null], lglfm: LinearGaussianLatentFeatureModel, nSamples: Int, thin: Int, rdg: RandomDataGenerator, newFeaturesTruncationDivisor: Double = 1000): Seq[FeatureAllocation[Null]] = {
    val nItems = lglfm.N
    var state = fa
    val nIterations = thin*nSamples
    val logNewFeaturesTruncationDivisor = log(newFeaturesTruncationDivisor)
    var results = List[FeatureAllocation[Null]]()
    var b = 1
    while (b <= nIterations) {
      for (i <- 0 until nItems) {
        val proposals = allPossibleConfigurationsAmongExistingKeepingSingletons(i, Z2fa(state, nItems))
        val logWeights = proposals.toIndexedSeq.map(fa => (fa, ibp.logDensity(fa) + lglfm.logLikelihood(fa2Z(fa))))
        state = fa2Z(rdg.nextItem(logWeights, onLogScale = true)._1)
        state = fa2Z(FeatureAllocation(nItems, Z2fa(state, nItems).filterNot(f => (f.size == 1 ) && (f.contains(i))).toVector))
        val featureWithOnlyI = Feature(null, i)
        @scala.annotation.tailrec
        def engine(weights: List[(FeatureAllocation[Null],Double)], max: Double): List[(FeatureAllocation[Null],Double)] = {
          val newCumState = weights.head._1.add(featureWithOnlyI)
          val newLogWeight = ibp.logDensity(newCumState) + lglfm.logLikelihood(fa2Z(newCumState))
          val expanded = (newCumState,newLogWeight) :: weights
          if ( newLogWeight < max - logNewFeaturesTruncationDivisor ) expanded
          else engine(expanded, if ( newLogWeight > max ) newLogWeight else max)
        }
        val faTmp = Z2fa(state,nItems)
        val setup1 = (faTmp,ibp.logDensity(faTmp) + lglfm.logLikelihood(state)) :: Nil
        val weights = engine(setup1,Double.NegativeInfinity).toIndexedSeq
        state = fa2Z(rdg.nextItem(weights, onLogScale = true)._1)
      }
      if (b % thin == 0) results = Z2fa(state, nItems) :: results
      b += 1
    }
    results.reverse
  }

}
