import Dependencies._

lazy val root = (project in file(".")).
  settings(
    name := "aibd2",
    organization := "org.ddahl",
    scalaVersion := "2.12.8",
    version      := "0.1.0.2-SNAPSHOT",
    //version      := "0.1.0.2",
    crossScalaVersions := Seq("2.11.12", "2.12.8"),
    scalacOptions ++= List("-feature", "-deprecation", "-unchecked", "-Xlint"),
    libraryDependencies += scalaTest % Test,
    libraryDependencies += "org.ddahl" %% "commonsmath" % "1.2.2.6",
    libraryDependencies += "org.ddahl" %% "rscala" % "3.2.6",
    libraryDependencies += "org.ddahl" %% "sdols" % "1.7.3.2",
    libraryDependencies += "org.apache.commons" % "commons-math3" % "3.6.1",
    resolvers += Resolver.bintrayRepo("dahl", "maven"),
    licenses := List(("Apache-2.0",url("https://www.apache.org/licenses/LICENSE-2.0")))
  )

