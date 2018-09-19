package org.ddahl.aibd

class EmpericalFeatureAllocationDistribution private (val count: Int, val map: Map[FeatureAllocation[Null], Int]) {

  def size = count

  def tally(fa: FeatureAllocation[Null]): EmpericalFeatureAllocationDistribution = {
    val fac = fa.leftOrderedForm
    new EmpericalFeatureAllocationDistribution(count+1, map + ((fac,map.getOrElse(fac,0)+1)))
  }

  def ++(that: EmpericalFeatureAllocationDistribution): EmpericalFeatureAllocationDistribution = {
    val keys = this.map.keys.toSet ++ that.map.keys.toSet
    val map = keys.map(k => (k,this.map.getOrElse(k,0) + that.map.getOrElse(k,0))).toMap
    new EmpericalFeatureAllocationDistribution(this.count + that.count, map)
  }

}

object EmpericalFeatureAllocationDistribution {

  def apply() = {
    val m = Map.empty[FeatureAllocation[Null], Int]
    new EmpericalFeatureAllocationDistribution(0, m)
  }

}
