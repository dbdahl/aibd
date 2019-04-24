package org.ddahl.aibd.distribution

trait HasMass[T] {

  val mass: Double
  val logMass: Double

  def updateMass(mass: Double): FeatureAllocationDistribution with HasMass[T]

}

