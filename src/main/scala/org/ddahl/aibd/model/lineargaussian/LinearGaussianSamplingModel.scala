package org.ddahl.aibd.model.lineargaussian

import org.ddahl.aibd._
import org.ddahl.aibd.parameter.{IndependentNormalsParameterDistribution, MultivariateNormalParameterDistribution}
import org.ddahl.commonsmath._
import org.apache.commons.math3.random.RandomDataGenerator
import org.apache.commons.math3.util.FastMath._

class LinearGaussianSamplingModel private (private val response: Array[Array[Double]], val maxNewFeatures: Int) {

  val nItems: Int = response.size
  require(nItems > 0,"Number of observations must be at least one.")

  val nResponses: Int = response(0).size

  val timer = TimeMonitor()

  private def sse(mean: Array[Array[Double]]) = {
    var sum = 0.0
    var i = 0
    while ( i < nItems ) {
      val responsei = response(i)
      val meani = mean(i)
      var j = 0
      while ( j < nResponses ) {
        val r = responsei(j) - meani(j)
        sum += r * r
        j += 1
      }
      i += 1
    }
    sum
  }

  private def logLikelihood(precision: Double, mean: Array[Array[Double]]): Double = {
    nItems*nResponses*(0.5*log(precision) + LinearGaussianSamplingModel.logNormalizingConstant) - 0.5*sse(mean)*precision
  }

  def logLikelihood(precision: Double, fa: FeatureAllocation[Vector[Double]]): Double = {
    logLikelihood(precision, LinearGaussianSamplingModel.mean(fa, nResponses))
  }

  def logLikelihood(precision: Array[Double], fas: Seq[FeatureAllocation[Vector[Double]]]): Array[Double] = {
    precision.zip(fas).par.map(x => logLikelihood(x._1, x._2)).toArray
  }

  private def sse(i: Int, mean: Array[Double]) = {
    var sum = 0.0
    var j = 0
    while ( j < nResponses ) {
      val r = response(i)(j) - mean(j)
      sum += r * r
      j += 1
    }
    sum
  }

  private def logLikelihood(i: Int, precision: Double, mean: Array[Double]): Double = {
    nResponses*(0.5*log(precision) + LinearGaussianSamplingModel.logNormalizingConstant) - 0.5*sse(i, mean)*precision
  }

  def logLikelihood(i: Int, precision: Double, fa: FeatureAllocation[Vector[Double]]): Double = {
    timer {
      logLikelihood(i, precision, LinearGaussianSamplingModel.mean(i, fa, nResponses))
    }
  }

  def maximumLikelihoodEstimate(precision: Array[Double], fas: Seq[FeatureAllocation[Vector[Double]]]): (Double, FeatureAllocation[Vector[Double]]) = {
    precision.zip(fas).par.maxBy(x => logLikelihood(x._1, x._2))
  }

  def maximumAPosterioriEstimate(precision: Array[Double], fas: Seq[FeatureAllocation[Vector[Double]]], featureAllocationDistribution: FeatureAllocationDistribution[Vector[Double]]) = {
    val parameterPrior = featureAllocationDistribution.parameterDistribution
    precision.zip(fas).par.maxBy( x => {
      logLikelihood(x._1, x._2) + featureAllocationDistribution.logDensityWithParameters(x._2,false)
    })
  }

  def mkLogLikelihood(precision: Double): (Int, FeatureAllocation[Vector[Double]]) => Double = logLikelihood(_, precision, _)

  def updatePrecision(fa: FeatureAllocation[Vector[Double]], rdg: RandomDataGenerator, shape: Double, rate: Double): Double = {
    val newShape = shape + nItems * nResponses / 2.0
    val newRate = rate + 0.5 * sse(LinearGaussianSamplingModel.mean(fa, nResponses))
    rdg.nextGamma(newShape, 1.0 / newRate)
  }

  def updateCoefficients(fa: FeatureAllocation[Vector[Double]], rdg: RandomDataGenerator, precision: Double, parameterPrior: IndependentNormalsParameterDistribution): FeatureAllocation[Vector[Double]] = {
    val priorPrecisionTimesMean = Array.tabulate(nResponses) { i => parameterPrior.precision(i) * parameterPrior.mean(i) }
    var faCurrent = fa
    for (fCurrent <- faCurrent) {
      faCurrent = faCurrent.remove(fCurrent)
      val n = fCurrent.size
      val demeaned = Array.ofDim[Double](nResponses)
      fCurrent.set.foreach { i =>
        val responsei = response(i)
        val mean = LinearGaussianSamplingModel.mean(i,faCurrent,nResponses)
        for ( j <- 0 until nResponses ) demeaned(j) += ( responsei(j) - mean(j) )
      }
      val posteriorPrecision = Array.tabulate(nResponses) { i =>
        parameterPrior.precision(i) + n * precision
      }
      val posteriorMean = Array.tabulate(nResponses) { i =>
        ( priorPrecisionTimesMean(i) + precision * demeaned(i) ) / posteriorPrecision(i)
      }
      val inpd = IndependentNormalsParameterDistribution.usingPrecision(posteriorMean,posteriorPrecision)
      val result = inpd.sample(rdg)
      faCurrent = faCurrent.add(fCurrent.replace(result))
    }
    faCurrent
  }

  def updateCoefficients(fa: FeatureAllocation[Vector[Double]], rdg: RandomDataGenerator, precision: Double, parameterPrior: MultivariateNormalParameterDistribution): FeatureAllocation[Vector[Double]] = {
    val priorPrecisionTimesMean = parameterPrior.precisionMatrix * parameterPrior.meanMatrix
    val precisionMatrix = precision *: MatrixFactory.identity(nResponses)
    val zero = MatrixFactory(nResponses,1)
    var faCurrent = fa
    for (fCurrent <- faCurrent) {
      faCurrent = faCurrent.remove(fCurrent)
      val items = fCurrent.set.toArray
      val demeaned = MatrixFactory(items.map { i =>
        val responsei = response(i)
        val meani = LinearGaussianSamplingModel.mean(i,faCurrent,nResponses)
        Array.tabulate(nResponses) { j => responsei(j) - meani(j) }
      },false)
      val sum = ( MatrixFactory.ones(1,items.length) * demeaned ).transpose
      val posteriorPrecisionMatrix = parameterPrior.precisionMatrix + demeaned.rows *: precisionMatrix
      val mvn = MultivariateNormalParameterDistribution.usingPrecision(zero, posteriorPrecisionMatrix)
      val posteriorMean = mvn.precisionCholesky.getSolver.solve( priorPrecisionTimesMean + precisionMatrix * sum )
      val result = mvn.sample(rdg,Some(posteriorMean))
      faCurrent = faCurrent.add(fCurrent.replace(result))
    }
    faCurrent
  }

}

object LinearGaussianSamplingModel {

  private final val logNormalizingConstant = -0.5 * log(2 * PI)

  def apply(response: Array[Array[Double]], maxNewFeatures: Int = 5): LinearGaussianSamplingModel = {
    new LinearGaussianSamplingModel(response.map(_.clone), maxNewFeatures)
  }

  def sse(response: Array[Array[Double]], fa: FeatureAllocation[Vector[Double]]): Double = {
    apply(response).sse(mean(fa,response(0).length))
  }

  def mean(i: Int, fa: FeatureAllocation[Vector[Double]], nResponses: Int): Array[Double] = {
    val mean = new Array[Double](nResponses)
    val nFeatures = fa.nFeatures
    val features = fa.toVector
    var k = 0
    while ( k < nFeatures ) {
      val f = features(k)
      if ( f.contains(i) ) {
        val param = f.parameter
        var m = 0
        while ( m < nResponses ) {
          mean(m) += param(m)
          m += 1
        }
      }
      k += 1
    }
    mean
  }

  def mean(fa: FeatureAllocation[Vector[Double]], nResponses: Int): Array[Array[Double]] = {
    val mean = Array.ofDim[Double](fa.nItems, nResponses)
    val nFeatures = fa.nFeatures
    val features = fa.toVector
    var k = 0
    while ( k < nFeatures ) {
      val f = features(k)
      val param = f.parameter
      var i = 0
      while ( i < fa.nItems ) {
        if ( f.contains(i) ) {
          var m = 0
          while ( m < nResponses ) {
            mean(i)(m) += param(m)
            m += 1
          }
        }
        i += 1
      }
      k += 1
    }
    mean
  }

  def sample(precision: Double, fa: FeatureAllocation[Vector[Double]], M: Int, rdg: RandomDataGenerator): Array[Array[Double]] = {
    val sd = 1/math.sqrt(precision)
    val result = Array.ofDim[Double](fa.nItems,M)
    for ( i <- 0 until fa.nItems ) {
      val mean = fa.featuresOf(i).foldLeft(new Array[Double](M)) { (sum,fa) =>
        var m = 0
        while ( m < M ) {
          sum(m) += fa.parameter(m)
          m += 1
        }
        sum
      }
      result(i) = Array.tabulate(M)(m => rdg.nextGaussian(mean(m),sd))
    }
    result
  }

  def sample(nSamples: Int, state: Option[(FeatureAllocation[Vector[Double]], Double, FeatureAllocationDistribution[Vector[Double]])], samplingModel: LinearGaussianSamplingModel, precisionShape: Double, precisionRate: Double, featureAllocationModel: FeatureAllocationDistribution[Vector[Double]], rdg: RandomDataGenerator, progressCallback: Int => Unit) = {
    var (fa, precision, faDistribution) = if ( state.isDefined ) state.get
    else ( featureAllocationModel.sample(rdg), precisionShape / precisionRate, featureAllocationModel)
    val monitorFAExisting = MCMCAcceptanceMonitor2()
    val monitorFASingletons = MCMCAcceptanceMonitor1()
    val monitorPermutation = MCMCAcceptanceMonitor1()
    val sweeper = SystematicSweep()
    sweeper.add("precision", 1) {
      precision = samplingModel.updatePrecision(fa, rdg, precisionShape, precisionRate)
    }
    sweeper.add("coefficients", 1) {
      faDistribution.parameterDistribution match {
        case pd: MultivariateNormalParameterDistribution =>
          fa = samplingModel.updateCoefficients(fa, rdg, precision, pd)
        case pd: IndependentNormalsParameterDistribution =>
          fa = samplingModel.updateCoefficients(fa, rdg, precision, pd)
      }
    }
    if ( faDistribution.isInstanceOf[IndianBuffetProcess[Vector[Double]]] ) {
      val ibp = faDistribution.asInstanceOf[IndianBuffetProcess[Vector[Double]]]
      sweeper.add("allocation", 1) {
        fa = MCMCSamplers.updateFeatureAllocationIBP(1, fa, ibp, samplingModel.mkLogLikelihood(precision), samplingModel.maxNewFeatures, rdg)
      }
    }
    else if ( faDistribution.isInstanceOf[AttractionIndianBuffetDistribution[Vector[Double]]]) {
      var aibd = faDistribution.asInstanceOf[AttractionIndianBuffetDistribution[Vector[Double]]]
      sweeper.add("allocation", 1) {
        fa = monitorFAExisting {
          MCMCSamplers.updateFeatureAllocationNeighborhoods(1, (x: Int) => x, fa, aibd, samplingModel.mkLogLikelihood(precision), rdg, false)
        }
      }
      sweeper.add("singletons", 1) {
        fa = monitorFASingletons {
          MCMCSamplers.updateFeatureAllocationSingletons(1, fa, aibd, samplingModel.mkLogLikelihood(precision), rdg)
        }
      }
      sweeper.add("permutation", 1) {
        aibd = monitorPermutation {
          aibd.updatePermutation(1, fa, rdg, false)
        }
      }
      faDistribution = aibd
    }
    val fas = Array.ofDim[FeatureAllocation[Vector[Double]]](nSamples)
    val precisions = Array.ofDim[Double](nSamples)
    val nFeatures = Array.ofDim[Int](nSamples)
    var i = 0
    while ( i < nSamples ) {
      sweeper(1)
      fas(i) = fa
      nFeatures(i) = fa.nFeatures
      precisions(i) = precision
      i += 1
      progressCallback(i)
    }
    val result = (
      (fa, precision, faDistribution),
      fas,
      nFeatures,
      precisions,
      Array(monitorFAExisting.rate1, monitorFAExisting.rate2, monitorFASingletons.rate, monitorPermutation.rate),
      sweeper.toString
    )
    result
  }

}

