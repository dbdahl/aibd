import Dependencies._

lazy val root = (project in file(".")).
  settings(
    name := "aibd2",
    organization := "org.ddahl.aibd2",
    scalaVersion := "2.12.6",
    version      := "0.1.0-SNAPSHOT",
    crossScalaVersions := Seq("2.11.12", "2.12.6"),
    scalacOptions ++= List("-feature", "-deprecation", "-unchecked", "-Xlint"),
    libraryDependencies += scalaTest % Test,
    libraryDependencies += "org.ddahl" %% "commonsmath" % "1.2-SNAPSHOT",
    libraryDependencies += "org.ddahl" %% "rscala" % "3.2.0-SNAPSHOT",
    libraryDependencies += "org.ddahl" %% "sdols" % "1.6-SNAPSHOT",
    libraryDependencies += "org.apache.commons" % "commons-math3" % "3.6.1"
  )

