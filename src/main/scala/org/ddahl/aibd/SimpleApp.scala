package org.ddahl.aibd

import org.ddahl.aibd.parameter.{GaussianParameterDistribution, NullParameterDistribution}
import org.ddahl.commonsmath._
import org.apache.commons.math3.linear.RealMatrix
import org.apache.commons.math3.random.RandomDataGenerator
import org.apache.commons.math3.util.FastMath.{exp, pow}

object EnumerationIntegration extends App {
  val R = org.ddahl.rscala.RClient()
  R.eval(
    """
      which <- c("New Hampshire","Iowa","Wisconsin","California","Nevada")
      which <- c("New Hampshire","Iowa","California","Nevada")
      data <- scale(USArrests[which,])
      distance <- as.matrix(dist(data))
      #distance <- matrix(1,nrow=4,ncol=4)
      #distance <- matrix(1,nrow=5,ncol=5)
    """.stripMargin)
  val maxNFeatures = 8
  val mass = 1.0
  val temperature = 1.0
  val similarity = Similarity.fromDistance(R.evalD2("distance"), (x: Double) => exp(-temperature * x))
  println(similarity.toText("%4.2f"))
  val aibd = MarginalizedAttractionIndianBuffetDistribution(mass, similarity, true, NullParameterDistribution)
  val y = time {
    aibd.integrateEnumeration(
      (fa,p) => MatrixFactory(fa.pairwiseAllocationMatrix) :* p,
      (x1: RealMatrix, x2: RealMatrix) => x1 + x2,
      maxNFeatures, false)
  }
  println(y)
}

object MonteCarloIntegration extends App {
  val R = org.ddahl.rscala.RClient()
  R.eval(
    """
      which <- c("New Hampshire","Iowa","Wisconsin","California","Nevada")
      #which <- c("New Hampshire","Iowa","California","Nevada")
      data <- scale(USArrests[which,])
      distance <- as.matrix(dist(data))
      #distance <- matrix(1,nrow=4,ncol=4)
      #distance <- matrix(1,nrow=5,ncol=5)
    """.stripMargin)
  val mass = 1.0
  val temperature = 1.0
  val nCores = Runtime.getRuntime.availableProcessors
  val rdg = new RandomDataGenerator()
  val similarity = Similarity.fromDistance(R.evalD2("distance"), (x: Double) => exp(-temperature * x))
  println(similarity.toText("%4.2f"))
  val aibd = MarginalizedAttractionIndianBuffetDistribution(mass, similarity, false, NullParameterDistribution)
  val exponent = 6
  val nSamplesIdeally = pow(10.0, exponent).toInt
  val nSamples = (nSamplesIdeally / nCores) * nCores
  println(exponent)
  val y = time {
    aibd.integrateMonteCarloMatrix(fa => MatrixFactory(fa.pairwiseAllocationMatrix), rdg, nSamples) :/ nSamples
  }
  println(y.toPrettyString())
}

object Simple extends App {
  val R = org.ddahl.rscala.RClient()
  R.eval(
    """
      which <- c("New Hampshire","Iowa","Wisconsin","California","Nevada")
      data <- scale(USArrests[which,])
      distance <- as.matrix(dist(data))
    """.stripMargin)
  val mass = 1.0
  val temperature = 1.0
  val rdg = new RandomDataGenerator()
  val distance = MatrixFactory(R.evalD2("distance"), false)
  println("Distance Matrix:")
  println(distance.toPrettyString("0.00"))
  println
  println("Similarity Matrix:")
  val similarity = Similarity.fromDistance(distance.getDataRef, (x: Double) => exp(-temperature * x))
  println(similarity.toText("%4.2f"))
  val aibd = MarginalizedAttractionIndianBuffetDistribution(mass, similarity, false, new GaussianParameterDistribution(0, 1))
  for (i <- 1 to 10) {
    val fa = aibd.sample(rdg)
    println
    println(s"Sampled features (draw $i):")
    println("... in set notation:")
    println(fa)
    println("... in matrix notation:")
  }
}
