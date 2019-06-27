name := "aibd"

organization := "org.ddahl"

//version := "0.1.1"
version := "0.1.0-SNAPSHOT"

scalaVersion := "2.12.8"
//crossScalaVersions := Seq("2.11.12", "2.12.8", "2.13.0")
crossScalaVersions := Seq("2.11.12", "2.12.8")
scalacOptions ++= Seq( "-deprecation", "-unchecked", "-feature" )

libraryDependencies ++= Seq(
  "org.apache.commons" % "commons-math3" % "3.6.1" withSources(),
  "org.scalactic" %% "scalactic" % "3.0.8",
  "org.scalatest" %% "scalatest" % "3.0.8" % "test",
  "org.scalacheck" %% "scalacheck" % "1.14.0" % "test"
)

libraryDependencies ++= {
  CrossVersion.partialVersion(scalaVersion.value) match {
    case Some((2, major)) if major >= 13 =>
      Seq("org.scala-lang.modules" %% "scala-parallel-collections" % "0.2.0")
    case _ =>
      Seq()
  }
}

Compile / unmanagedJars := {
  val rPackages = Seq("commonsMath","sdols")
  rPackages.flatMap { p =>
    import scala.sys.process._
    import scala.language.postfixOps
    import java.io.File
    val exe = if ( sys.env.getOrElse("R_HOME","") == "" ) "R" else {
      Seq(sys.env("R_HOME"),"bin","R").mkString(File.separator)
    }
    val output = Seq(exe,"--slave","-e",s"writeLines(rscala:::jarsOfPackage('${p}','${scalaBinaryVersion.value}'))") !!
    val cells = output.split("\n").toSeq
    println(cells)
    cells.map { path => Attributed.blank(new File(path)) }
  }
}

