package org.ddahl.aibd.model.regression

import org.apache.commons.math3.random.RandomDataGenerator
import org.apache.commons.math3.util.FastMath.{exp, pow, sqrt}
import org.ddahl.aibd.parameter.MultivariateNormalParameterDistribution
import org.ddahl.aibd._
import org.ddahl.commonsmath._

object RegressionSimulationStudySetup extends App {

  //val Array(nItems, seed1, seed2) = args.map(_.toInt)
  val Array(nItems, seed1, seed2) = Array(64, 1027078439, 1409452087)
  println(seed1 + " " + seed2)
  val rdg = new RandomDataGenerator()
  rdg.reSeed(seed1)
  val R = org.ddahl.rscala.RClient()
  R.eval(s"set.seed(${seed2})")
  R.eval(
    s"""
       |nItems <- $nItems
       |loc <- t(sapply(1:nItems,function(i) matrix(c(0,0),nrow=1) + rnorm(2) %*% chol(matrix(c(1,0.9,0.9,1),nrow=2))))
       |distance <- as.matrix(dist(scale(loc)))
       |X <- cbind(1,sample(9:16,nItems,prob=c(1,2,3,6,5,3,3,6),replace=TRUE))
    """.stripMargin)
  val nCovariates = R.evalI0("ncol(X)")
  val temperature = 2.0
  val mass = 1.0
  val pdMean = Array(Array(0.0), Array(0.0))
  val pdCov = Array(Array(0.1, 0.0), Array(0.0, pow(0.025, 2)))
  val precisionShape = 60.0
  val precisionRate = 0.3
  val similarity = Similarity.fromDistance(R.evalD2("distance"), d => exp(-temperature * d))
  val parameterDistribution = MultivariateNormalParameterDistribution.usingCovariance(MatrixFactory(pdMean), MatrixFactory(pdCov))
  var aibd = AttractionIndianBuffetDistribution(mass, Permutation.random(nItems, rdg), similarity, parameterDistribution)
  val originalFA = aibd.sample(rdg)
  if (originalFA.nFeatures != 4) {
    println("I still haven't found what I'm looking for.  -- U2")
    sys.exit()
  }
  val originalPrecision = rdg.nextGamma(precisionShape, 1.0 / precisionRate)
  R.eval("temperature <- %-", temperature)
  R.eval("mass <- %-", mass)
  R.eval("pdMean <- %-", pdMean)
  R.eval("pdCov <- %-", pdCov)
  R.eval("precisionShape <- %-", precisionShape)
  R.eval("precisionRate <- %-", precisionRate)
  R.eval("precision <- %-", originalPrecision)
  R.eval("Z <- %-", originalFA.toMatrix)
  R.eval(
    """
      |labels <- apply(Z,1,function(x) sum(2^(0:(length(x)-1))*x))
      |table <- table(labels)
      |freq <- apply(Z,2,sum)
      |print(Z)
      |print(freq)
      |print(length(table))
      |print(table)
    """.stripMargin)
  println(originalFA.nFeatures)
  println(originalFA.rate)
  R.eval("features <- %-", originalFA.map(f => f.parameter.toArray).toArray)
  R.eval(
    s"""
       |print(features)
       |B <- if ( ncol(Z) == 0 ) matrix(0,nrow=nrow(X),ncol=ncol(X)) else Z %*% features
       |# overall <- matrix(rep(c(2.7,0.08),each=nItems),ncol=ncol(X))
       |# y <- apply(X*(overall+B),1,sum) + rnorm(nItems,sd=1/sqrt(precision))
       |y <- apply(X*B,1,sum) + rnorm(nItems,sd=1/sqrt(precision))
       |save(list=ls(),file="${nItems}_${seed1}_${seed2}.Rbin")
    """.stripMargin)
  R.eval(
    """
      |plot(loc[,1],loc[,2],type="n")
      |text(loc[,1],loc[,2],labels=labels)
      |locator(1)
    """.stripMargin)

}

object RegressionSimulationStudyRun extends App {

  var cla = Array("64_1027078439_1409452087.Rbin", "5", "10", "234523", "0", "aibd").toList
  //var cla = args.toList
  val rbin = cla.head
  cla = cla.tail
  val List(nBurnIn, nDraws, seed3, missingObservation) = cla.take(4).map(_.toInt)
  cla = cla.drop(4)
  val model = cla.head
  cla = cla.tail
  val rdg = new RandomDataGenerator()
  rdg.reSeed(seed3)
  val R = org.ddahl.rscala.RClient()
  R.eval("setwd('~/docs/devel/aibd/demonstration/regression')")
  R.eval(s"""load('$rbin')""")
  val nCovariates = R.evalI0("ncol(X)")
  val temperature = R.evalD0("temperature")
  val nItems = R.evalI0("nItems")
  val similarity = model match {
    case "aibd" => Similarity.fromDistance(R.evalD2("distance"), d => exp(-temperature * d))
    case "ibp" => Similarity.uniform(nItems)
    case _ => sys.error("Unrecognized distribution")
  }
  val mass = R.evalD0("mass")
  val parameterDistribution = MultivariateNormalParameterDistribution.usingCovariance(MatrixFactory(R.evalD2("pdMean")), MatrixFactory(R.evalD2("pdCov")))
  val precisionShape = R.evalD0("precisionShape")
  val precisionRate = R.evalD0("precisionRate")
  val monitorFANeighborhoods = MCMCAcceptanceMonitor2()
  val monitorFASingletons = MCMCAcceptanceMonitor1()
  val monitorPermutation = MCMCAcceptanceMonitor1()
  var aibd = AttractionIndianBuffetDistribution(mass, Permutation.random(nItems, rdg), similarity, parameterDistribution)
  var rsm = RegressionSamplingModel(R.evalD1("y"), R.evalD2("X"))
  var precision = precisionShape / precisionRate
  val drawsUnivariateQuantities = Array.ofDim[Double](nDraws, 5)
  val drawsCoefficients = Array.ofDim[Double](nDraws, nItems, nCovariates)
  rdg.reSeed(seed3)
  var fa = aibd.sample(rdg)
  //val sweeper = RandomSweep(rdg)
  val sweeper = SystematicSweep()
  sweeper.add("coefficients", 1) {
    fa = rsm.updateCoefficients(fa, rdg, precision, parameterDistribution)
  }
  sweeper.add("allocation", 1) {
    fa = monitorFANeighborhoods {
      MCMCSamplers.updateFeatureAllocationNeighborhoods(1, (x: Int) => x, fa, aibd, rsm.mkLogLikelihood(precision), rdg, false)
    }
  }
  sweeper.add("singletons", 1) {
    fa = monitorFASingletons {
      MCMCSamplers.updateFeatureAllocationSingletons(1, fa, aibd, rsm.mkLogLikelihood(precision), rdg)
    }
  }
  sweeper.add("permutation", 4) {
    aibd = monitorPermutation {
      aibd.updatePermutation(1, fa, rdg, false)
    }
  }
  sweeper.add("precision", 1) {
    precision = rsm.updatePrecision(fa, rdg, precisionShape, precisionRate)
  }
  sweeper(nBurnIn)
  var i = 0
  while (i < nDraws) {
    if (i % 1 == 0) {
      println(fa.nFeatures)
      println(sweeper)
      println(rsm.timer)
    }
    val (mean, currentY) = if (missingObservation >= 0) {
      val mean = rsm.mean(missingObservation, fa)
      val currentY = rdg.nextGaussian(mean, 1 / sqrt(precision))
      //rsm = rsm.updatedResponse(missingObservation, currentY)
      (mean, currentY)
    } else {
      (0.0, 0.0)
    }
    sweeper(1)
    drawsUnivariateQuantities(i) = Array(currentY, mean, precision, fa.nFeatures, fa.rate)
    drawsCoefficients(i) = rsm.coefficients(fa).map(_.toArray).toArray
    i += 1
  }
  R.eval("intercepts <- %-", drawsCoefficients.map(_.map(_ (0))))
  R.eval("slopes <- %-", drawsCoefficients.map(_.map(_ (1))))
  R.eval("univariates <- %-", drawsUnivariateQuantities)
  R.eval("mcmcRates <- %-", Array(monitorPermutation.rate, monitorFASingletons.rate, monitorFANeighborhoods.rate1, monitorFANeighborhoods.rate2))
  R.eval("sweeperStats <- %-", sweeper.toString)
  R.eval("args <- %-", args)
  R.eval(
    s"""
       |outdir <- tools::file_path_sans_ext('${rbin}')
       |dir.create(outdir, showWarnings=FALSE)
       |missingObservation <- $missingObservation+1
       |save(intercepts,slopes,univariates,mcmcRates,missingObservation,sweeperStats,args,file=paste(outdir,"/${model}_${missingObservation}_$seed3.Rbin",sep=""))
    """.stripMargin)
  println(monitorPermutation.rate + " " + monitorFASingletons.rate + " " + monitorFANeighborhoods.rate1 + " " + monitorFANeighborhoods.rate2)
  println(sweeper)

}
