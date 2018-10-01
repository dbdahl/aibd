installJARs <- function() {
  unlink("java",recursive=TRUE)
  unlink("inst",recursive=TRUE)
  filename <- tempfile()
  on.exit(function() delete(filename))
  download.file("https://dahl.byu.edu/public/aibd2-JARs.zip",destfile=filename)
  unzip(filename)
}
