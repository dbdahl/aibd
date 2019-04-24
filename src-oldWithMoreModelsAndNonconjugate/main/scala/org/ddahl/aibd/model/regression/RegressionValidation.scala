package org.ddahl.aibd.model.regression

import org.apache.commons.math3.random.RandomDataGenerator
import org.apache.commons.math3.util.FastMath.{exp, pow}
import org.ddahl.aibd.parameter.MultivariateNormalParameterDistribution
import org.ddahl.aibd._
import org.ddahl.commonsmath._

object RegressionValidation extends App {

  val Array(nReps, nBurnIn, nDraws, nItems) = args.init.map(_.toInt)
  val alpha = args.last.toDouble
  val rdg = new RandomDataGenerator()
  //val seeds = (912445L, 13455422L)
  val seeds = ((math.random * 10000).toLong, (math.random * 10000).toLong)
  rdg.reSeed(seeds._1)
  val R = org.ddahl.rscala.RClient()
  R.eval(s"set.seed(${seeds._2})")
  R.eval(
    s"""
       |nItems <- $nItems
       |loc <- t(sapply(1:nItems,function(i) matrix(c(0,0),nrow=1) + rnorm(2) %*% chol(matrix(c(1,0.8,0.8,1),nrow=2))))
       |distance <- as.matrix(dist(scale(loc)))
       |X <- cbind(1,sample(9:16,nItems,prob=c(1,2,3,6,5,3,3,6),replace=TRUE))
       |hit <- matrix(NA,nrow=$nReps,ncol=3)
       |probs <- c($alpha/2,1-$alpha/2)
    """.stripMargin)
  val nCovariates = R.evalI0("ncol(X)")
  val temperature = 1.0
  val similarity = Similarity.fromDistance(R.evalD2("distance"), d => exp(-temperature * d))
  val mass = 1.0
  val parameterDistribution = MultivariateNormalParameterDistribution.usingCovariance(MatrixFactory(Array(Array(0.0), Array(0.0))), MatrixFactory(Array(Array(0.1, 0.0), Array(0.0, pow(0.02, 2)))))
  val precisionShape = 30.0
  val precisionRate = 0.3
  val monitorFANeighborhoods = MCMCAcceptanceMonitor2()
  val monitorFASingletons = MCMCAcceptanceMonitor1()
  val monitorPermutation = MCMCAcceptanceMonitor1()
  for (rep <- 0 until nReps) {
    println(rep.toDouble / nReps)
    var aibd = AttractionIndianBuffetDistribution(mass, Permutation.random(nItems, rdg), similarity, parameterDistribution)
    val originalFA = aibd.sample(rdg)
    val originalPrecision = rdg.nextGamma(precisionShape, 1.0 / precisionRate)
    R.eval("precision <- %-", originalPrecision)
    R.eval("Z <- %-", originalFA.toMatrix)
    R.eval("features <- %-", originalFA.map(f => f.parameter.toArray).toArray)
    R.eval(
      """
        |B <- if ( ncol(Z) == 0 ) matrix(0,nrow=nrow(X),ncol=ncol(X)) else Z %*% features
        |# overall <- matrix(rep(c(2.8,0.08),each=nItems),ncol=ncol(X))
        |# y <- apply(X*(overall+B),1,sum) + rnorm(nItems,sd=1/sqrt(precision))
        |y <- apply(X*B,1,sum) + rnorm(nItems,sd=1/sqrt(precision))
      """.
        stripMargin)
    val rsm = RegressionSamplingModel(R.evalD1("y"), R.evalD2("X"))
    var precision = precisionShape / precisionRate
    var fa = aibd.sample(rdg)
    val drawsUnivariateQuantities = Array.ofDim[Double](nDraws, 1)
    val drawsCoefficients = Array.ofDim[Double](drawsUnivariateQuantities.length, nItems, nCovariates)
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
    while (i < drawsUnivariateQuantities.length) {
      sweeper(1)
      drawsUnivariateQuantities(i) = Array(precision)
      drawsCoefficients(i) = rsm.coefficients(fa).map(_.toArray).toArray
      i += 1
    }
    R.eval("rate <- %-", originalFA.rate)
    R.eval("intercepts <- %-", drawsCoefficients.map(_.map(_ (0))))
    R.eval("slopes <- %-", drawsCoefficients.map(_.map(_ (1))))
    R.eval("univariates <- %-", drawsUnivariateQuantities)
    R.eval(
      s"""
         |ci <- t(apply(intercepts,2,function(x) quantile(x,probs=probs)))
         |hit[$rep+1,1] <- mean(ci[,1] <= B[,1] & B[,1] <= ci[,2])
         |ci <- t(apply(slopes,2,function(x) quantile(x,probs=probs)))
         |hit[$rep+1,2] <- mean(ci[,1] <= B[,2] & B[,2] <= ci[,2])
         |hit[$rep+1,3] <- prod(quantile(univariates[,1],probs=probs)-precision)<0
    """.stripMargin)
  }
  R.eval(
    s"""
       |mean <- apply(hit,2,mean)
       |sd <- apply(hit,2,sd)
       |lower <- mean - 1.96*sd/sqrt($nReps)
       |upper <- mean + 1.96*sd/sqrt($nReps)
       |print(cbind(lower,mean,upper))
    """.stripMargin)
  println(monitorPermutation.rate + " " + monitorFASingletons.rate + " " + monitorFANeighborhoods.rate1 + " " + monitorFANeighborhoods.rate2)

}
