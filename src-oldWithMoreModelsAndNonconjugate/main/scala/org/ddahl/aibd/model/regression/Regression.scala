package org.ddahl.aibd.model.regression

import org.apache.commons.math3.random.RandomDataGenerator
import org.apache.commons.math3.util.FastMath.exp
import org.ddahl.aibd._
import org.ddahl.aibd.parameter.MultivariateNormalParameterDistribution
import org.ddahl.commonsmath._

object Regression extends App {

  val rdg = new RandomDataGenerator()
  //val seeds = (912445L, 13455422L)
  val seeds = ((math.random * 10000).toLong, (math.random * 10000).toLong)
  rdg.reSeed(seeds._1)
  val R = org.ddahl.rscala.RClient()
  R.eval(s"set.seed(${seeds._2})")
  R.eval(
    """
      data(airquality)
      #airquality <- airquality[sample(nrow(airquality),5),]
      missing <- apply(airquality,1,function(x) any(is.na(x)))
      # t.test(airquality$Ozone[missing],airquality$Ozone[!missing])
      # t.test(airquality$Solar.R[missing],airquality$Solar.R[!missing])
      airquality <- airquality[!missing,]
      doy <- as.Date(sprintf("%s-%s-%s",1973,airquality$Month,airquality$Day))
      dd <- as.numeric(doy)
      y <- airquality[,"Ozone"]
      X <- cbind(1,scale(as.matrix(airquality[,c("Solar.R","Wind","Temp")])))
      distance <- as.matrix(dist(scale(doy)))
      nItems <- nrow(X)
      nCovariates <- ncol(X)
    """.stripMargin)
  val nItems = R.evalI0("nItems")
  val stdDevOfBetas = 20.0
  val precisionShape = 4.0
  val precisionRate = 400.0
  val mass = 1.0
  val temperature = 1.0
  val nPerShuffle = nItems
  val nCovariates = R.evalI0("nCovariates")
  val distance = MatrixFactory(R.evalD2("distance"), false)
  val similarity = Similarity.fromDistance(distance.getDataRef, (x: Double) => exp(-temperature * x))
  println(similarity.toText("%4.2f"))
  val permutation = Permutation.random(nItems, rdg).nPerShuffle(nPerShuffle)
  val parameterDistribution = MultivariateNormalParameterDistribution.usingCovariance(MatrixFactory(Array.fill(nCovariates, 1)(0.0), false), MatrixFactory.identity(nCovariates) :* (stdDevOfBetas * stdDevOfBetas))
  var aibd = AttractionIndianBuffetDistribution(mass, permutation, similarity, parameterDistribution)
  val rsm = RegressionSamplingModel(R.evalD1("y"), R.evalD2("X"))
  var fa = aibd.sample(rdg)
  var nFeatures = 0
  val covariates = MatrixFactory(rsm.covariates.map(_.toArray).toArray, false)
  val result = Array.ofDim[Double](50, nItems)
  var precision = precisionShape / precisionRate
  val monitorFAExisting = MCMCAcceptanceMonitor2()
  val monitorFASingletons = MCMCAcceptanceMonitor1()
  val monitorPermutation = MCMCAcceptanceMonitor1()
  //val sweeper = RandomSweep(rdg)
  val sweeper = SystematicSweep()
  sweeper.add("precision", 1) {
    precision = rsm.updatePrecision(fa, rdg, precisionShape, precisionRate)
  }
  sweeper.add("coefficients", 1) {
    fa = rsm.updateCoefficients(fa, rdg, precision, parameterDistribution)
  }
  sweeper.add("allocation", 1) {
    fa = monitorFAExisting {
      MCMCSamplers.updateFeatureAllocationNeighborhoods(1, (x: Int) => x, fa, aibd, rsm.mkLogLikelihood(precision), rdg, false)
    }
  }
  sweeper.add("singletons", 1) {
    fa = monitorFASingletons {
      MCMCSamplers.updateFeatureAllocationSingletons(1, fa, aibd, rsm.mkLogLikelihood(precision), rdg)
    }
  }
  sweeper.add("permutation", 1) {
    aibd = monitorPermutation {
      aibd.updatePermutation(1, fa, rdg, false)
    }
  }
  var i = 0
  println(fa.nFeatures)
  time {
    while (i < result.length) {
      sweeper(1)
      nFeatures += fa.nFeatures
      //println(monitorFASingletons.rate)
      println(sweeper)
      //println(1 / sqrt(precision))
      //println(fa)
      println(fa.nFeatures)
      //println(fa.sharedFeaturesCountMatrix(nItems).toPrettyString("0"))
      //println(aibd.permutation)
      result(i) = rsm.mean(fa).toArray
      i += 1
    }
  }
  println(nFeatures / result.length.toDouble)
  println(monitorFAExisting.rate1 + " " + monitorFAExisting.rate2)
  println(monitorFASingletons.rate)
  println(monitorPermutation.rate)
  R.eval("result <- %-", result)
  println(R.evalD1("apply(result,2,mean)").mkString(" "))
  println(R.evalD1("y").mkString(" "))
}
