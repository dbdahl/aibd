import org.ddahl.aibd._
import org.ddahl.aibd.parameter.{BernoulliParameterDistribution, NullParameterDistribution}
import org.apache.commons.math3.random.RandomDataGenerator
import org.apache.commons.math3.util.FastMath._
import org.scalatest.FlatSpec

object Utils {

  def sumOverAll(mass: Double, nItems: Int, maxNFeatures: Int): Double = {
    val aibd = IndianBuffetProcess(mass, nItems, NullParameterDistribution)
    val list = FeatureAllocation.enumerate(nItems, List[Null](null), maxNFeatures)
    list.map(fa => {
      exp(aibd.logDensity(fa, true))
    }).sum
  }

  def between(x: Double, lower: Double, upper: Double) = {
    assert(x >= lower && x <= upper)
  }

  def almostOne(x: Double, lower: Double) = between(x, lower, 1.0)

  def compareDistribution[A](fad: FeatureAllocationDistribution[A], efad: EmpericalFeatureAllocationDistribution, tol1: Double = 5.0E-3, tol2: Double = 5.0E-3, tol3: Double = 0.995) = {
    val logNSamples = log(efad.size)
    val dist = fad.dropParameter
    val (jensenShannonDivergence, totalVarianceDistance, totalProbability) = {
      var s1 = 0.0
      var s2 = 0.0
      var s3 = 0.0
      efad.map.foreach { case (fa, count) =>
        val p1 = dist.logDensity(fa, true)
        val p2 = log(count) - logNSamples
        s1 += (0.5 * ((exp(p1) - exp(p2)) * (p1 - p2)))
        s2 = max(s2, abs(exp(p1) - exp(p2)))
        s3 += exp(p1)
      }
      (s1, s2, s3)
    }
    println(s"$totalVarianceDistance $jensenShannonDivergence $totalProbability")
    between(totalVarianceDistance, 0.0, tol1)
    between(jensenShannonDivergence, 0.0, tol2)
    almostOne(totalProbability, tol3)
  }

}

class IndianBuffetTests extends FlatSpec {

  import Utils._

 "IndianBuffetProcess with NullParameterDistribution" should "have probability mass function that sums to one" in {
    almostOne(sumOverAll(1.0, 2, 5), 0.995)
    almostOne(sumOverAll(2.0, 2, 8), 0.995)
    almostOne(sumOverAll(1.0, 3, 6), 0.995)
    almostOne(sumOverAll(2.0, 3, 9), 0.995)
  }

  "IndianBuffetProcess with BernoulliParameterDistribution" should "have probability mass function that sums to one" in {
    almostOne(sumOverAll(0.5, 2, 4), 0.995)
    almostOne(sumOverAll(0.5, 3, 4), 0.995)
  }

  "AttractionIndianBuffetDistribution with Implicit and Explicit Uniform Similarity" should "be the same" in {
    val pd = NullParameterDistribution
    val mass = 1.0
    val nItems = 4
    val maxNFeatures = 4
    val uniformMatrix = Array.fill(nItems, nItems)(1.2314)
    val similarity = Similarity.fromDistance(uniformMatrix, identity)
    val permutation = Permutation.natural(nItems)
    val aibd1 = AttractionIndianBuffetDistribution(mass, permutation, similarity, pd)
    val aibd2 = AttractionIndianBuffetDistribution(mass, permutation, Similarity.uniform(nItems), pd)
    assert(!aibd1.similarity.isUniform)
    assert(aibd2.similarity.isUniform)
    FeatureAllocation.foreach(nItems, pd.discreteSupport.get, maxNFeatures) { fa =>
      assert(aibd1.logDensity(fa, true) == aibd2.logDensity(fa, true))
    }
  }

  "AttractionIndianBuffetDistribution with Uniform Distances and IndianBuffetProcess" should "be the same" in {
    val pd = NullParameterDistribution
    val mass = 1.0
    val nItems = 4
    val maxNFeatures = 4
    val uniformMatrix = Array.fill(nItems, nItems)(1.2314)
    val similarity = Similarity.fromDistance(uniformMatrix, identity)
    val permutation = Permutation.natural(nItems)
    val aibd1 = AttractionIndianBuffetDistribution(mass, permutation, similarity, pd)
    val aibd2 = IndianBuffetProcess(mass, nItems, pd)
    assert(!aibd1.similarity.isUniform)
    FeatureAllocation.foreach(nItems, pd.discreteSupport.get, maxNFeatures) { fa =>
      assert((aibd1.logDensity(fa, true) - aibd2.logDensity(fa, true)).abs <= 0.00000001)
    }
  }

  "AttractionIndianBuffetDistribution with Nonuniform Distances" should "have probability mass function that sums to one" in {
    val R = org.ddahl.rscala.RClient()
    R.eval(
      """
      which <- c("New Hampshire","Iowa","California","Nevada")
      data <- scale(USArrests[which,])
      distance <- as.matrix(dist(data))
      """.stripMargin)
    val temperature = 1.0
    val mass = 1.0
    val similarity = Similarity.fromDistance(R.evalD2("distance"), (x: Double) => exp(-temperature * x))
    val rdg = new RandomDataGenerator()
    val perm = Permutation.random(similarity.nItems, rdg)
    val parameterDistribution = NullParameterDistribution
    val aibd = AttractionIndianBuffetDistribution(mass, perm, similarity, parameterDistribution)
    val maxNFeatures = 7
    var one = 0.0
    FeatureAllocation.foreach(similarity.nItems, parameterDistribution.discreteSupport.get, maxNFeatures) { fa =>
      one += exp(aibd.logDensity(fa, true))
    }
    almostOne(one, 0.995)
  }

  "MarginalizedAttractionIndianBuffetDistribution with Nonuniform Distances" should "have probability mass function that sums to one" in {
    val R = org.ddahl.rscala.RClient()
    R.eval(
      """
      which <- c("New Hampshire","Iowa","California","Nevada")
      data <- scale(USArrests[which,])
      distance <- as.matrix(dist(data))
      """.stripMargin)
    val temperature = 1.0
    val mass = 1.0
    val similarity = Similarity.fromDistance(R.evalD2("distance"), (x: Double) => exp(-temperature * x))
    val parameterDistribution = NullParameterDistribution
    val maibd = MarginalizedAttractionIndianBuffetDistribution(mass, similarity, true, parameterDistribution)
    val maxNFeatures = 5
    var one = 0.0
    FeatureAllocation.foreach(similarity.nItems, parameterDistribution.discreteSupport.get, maxNFeatures) { fa =>
      one += exp(maibd.logDensity(fa, true))
    }
    almostOne(one, 0.98)
  }

 "IndianBuffetProcess with BernoulliParameterDistribution" should "sample from the entire sample space and its empirical distribution should match the probability mass function" in {
    val nItems: Int = 3
    val mass: Double = 1.0
    val prob: Double = 0.3
    val nSamples = 1000000
    val parameterDistribution = BernoulliParameterDistribution(prob)
    val ibp = IndianBuffetProcess(mass, nItems, parameterDistribution)
    val rdg = new RandomDataGenerator()
    time {
      compareDistribution(ibp, ibp.toEmpericalDistributionViaDirect(rdg, nSamples), 0.015, 0.015, 0.99)
    }
  }

  "AttractionIndianBuffetDistribution with BernoulliParameterDistribution" should "sample from the entire sample space and its empirical distribution should match the probability mass function" in {
    val mass: Double = 0.75
    val nSamples = 1000000
    // val parameterDistribution = NullParameterDistribution
    val parameterDistribution = BernoulliParameterDistribution(0.4)
    val R = org.ddahl.rscala.RClient()
    R.eval(
      """
      which <- c("New Hampshire","Iowa","Wisconsin","California","Nevada")
      which <- c("New Hampshire","Iowa","California","Nevada")
      which <- c("New Hampshire","Iowa","Nevada")
      data <- scale(USArrests[which,])
      distance <- as.matrix(dist(data))
      """.stripMargin)
    val temperature = 1.0
    val similarity = Similarity.fromDistance(R.evalD2("distance"), (x: Double) => exp(-temperature * x))
    val rdg = new RandomDataGenerator()
    val perm = Permutation.random(similarity.nItems, rdg)
    val aibd = AttractionIndianBuffetDistribution(mass, perm, similarity, parameterDistribution)
    time {
      compareDistribution(aibd, aibd.toEmpericalDistributionViaDirect(rdg, nSamples), 0.007, 0.007)
    }
  }

  "MarginalizedAttractionIndianBuffetDistribution with NullParameterDistribution" should "sample from the entire sample space and its empirical distribution should match the probability mass function" in {
    val mass: Double = 0.75
    val nSamples = 100000
    val parameterDistribution = NullParameterDistribution
    val R = org.ddahl.rscala.RClient()
    R.eval(
      """
      which <- c("New Hampshire","Iowa","Wisconsin","California","Nevada")
      which <- c("New Hampshire","Iowa","California","Nevada")
      which <- c("New Hampshire","Iowa","Nevada")
      data <- scale(USArrests[which,])
      distance <- as.matrix(dist(data))
      """.stripMargin)
    val temperature = 1.0
    val similarity = Similarity.fromDistance(R.evalD2("distance"), (x: Double) => exp(-temperature * x))
    val maibd = MarginalizedAttractionIndianBuffetDistribution(mass, similarity, true, parameterDistribution)
    val rdg = new RandomDataGenerator()
    time {
      compareDistribution(maibd, maibd.toEmpericalDistributionViaDirect(rdg, nSamples), 5.5E-3, 5.5E-3, 0.995)
    }
  }

  "Direct Sampler for IBP with BernoulliParameterDistribution" should "matches theory" in {
    val nItems: Int = 3
    val mass: Double = 0.75
    val nSamples = 100000
    val parameterDistribution = BernoulliParameterDistribution(0.4)
    val ibp = IndianBuffetProcess(mass, nItems, parameterDistribution)
    val rdg = new RandomDataGenerator()
    time {
      compareDistribution(ibp, ibp.toEmpericalDistributionViaDirect(rdg, nSamples), 0.006, 0.027, 0.982)
    }
  }

  "MCMC Sampler for Attraction IB Distribution with BernoulliParameterDistribution" should "matches theory" in {
    val mass: Double = 0.75
    val nSamples = 100000
    val parameterDistribution = BernoulliParameterDistribution(0.4)
    val R = org.ddahl.rscala.RClient()
    R.eval(
      """
      which <- c("New Hampshire","Iowa","California","Nevada")
      which <- c("New Hampshire","Iowa","Nevada")
      data <- scale(USArrests[which,])
      distance <- as.matrix(dist(data))
      """.stripMargin)
    val temperature = 1.0
    val similarity = Similarity.fromDistance(R.evalD2("distance"), (x: Double) => exp(-temperature * x))
    val rdg = new RandomDataGenerator()
    val perm = Permutation.random(similarity.nItems, rdg)
    val aibd = AttractionIndianBuffetDistribution(mass, perm, similarity, parameterDistribution)
    time {
      compareDistribution(aibd, aibd.toEmpericalDistributionViaMCMC(rdg, nSamples), 0.0045, 0.026, 0.982)
    }
  }

}

class IndianBuffetTests2 extends FlatSpec {

  import Utils._

  "MCMC Sampler for IBP with BernoulliParameterDistribution" should "matches theory" in {
    val mass: Double = 0.75
    val nSamples = 100000
    val parameterDistribution = BernoulliParameterDistribution(0.4)
    val rdg = new RandomDataGenerator()
    val ibp = IndianBuffetProcess(mass, 3, parameterDistribution)
    time {
      compareDistribution(ibp, ibp.toEmpericalDistributionViaMCMC(rdg, nSamples), 0.0045, 0.026, 0.982)
    }
  }

}
