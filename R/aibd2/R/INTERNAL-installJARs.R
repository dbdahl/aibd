#' @importFrom utils download.file unzip
#'
installJARs <- function() {
  unlink("java",recursive=TRUE)
  unlink("inst",recursive=TRUE)
  filename <- tempfile()
  on.exit(function() unlink(filename))
  download.file("https://dahl.byu.edu/public/aibd2-JARs.zip",destfile=filename)
  unzip(filename)
}
