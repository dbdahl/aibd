name := "aibd"

organization := "org.ddahl"

version := "0.1.10"
//version := "0.1.10-SNAPSHOT"

scalaVersion := "2.13.6"
crossScalaVersions := Seq("2.11.12", "2.12.14", "2.13.6")
scalacOptions ++= Seq( "-deprecation", "-unchecked", "-feature" )

libraryDependencies ++= Seq(
  "org.apache.commons" % "commons-math3" % "3.6.1" withSources(),
  "org.scalactic" %% "scalactic" % "3.2.9",
  "org.scalatest" %% "scalatest" % "3.2.9" % "test",
  "org.scalacheck" %% "scalacheck" % "1.14.1" % "test"
)

libraryDependencies ++= {
  CrossVersion.partialVersion(scalaVersion.value) match {
    case Some((2, major)) if major >= 13 =>
      Seq("org.scala-lang.modules" %% "scala-parallel-collections" % "1.0.3")
    case _ =>
      Seq()
  }
}

Compile / unmanagedJars := {
  val rPackages = Seq("commonsMath")
  rPackages.flatMap { p =>
    import scala.sys.process._
    import java.io.File
    val exe = if ( sys.env.getOrElse("R_HOME","") == "" ) "R" else {
      Seq(sys.env("R_HOME"),"bin","R").mkString(File.separator)
    }
    val output = Seq(exe,"--slave","-e",s"writeLines(rscala:::jarsOfPackage('${p}','${scalaBinaryVersion.value}'))").!!
    val cells = output.split(sys.props("line.separator")).toSeq
    println(s"JARs from '${p}' package:")
    println(cells.mkString("\n"))
    cells.map { path => Attributed.blank(new File(path)) }
  }
}

