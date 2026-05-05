# Package-wide global variables
.XOMultinomEnv <- new.env()

.onLoad <- function(lib, pkg) {
  if (getRversion() >= "2.15.1") {
    utils::globalVariables(c(".", ""))
  }
  .XOMultinomEnv$path.to.me <- tools::file_path_as_absolute(lib)
}

.onAttach <- function(libname, pkgname) {
  packageStartupMessage(sprintf(
    "Package %s (%s) loaded.\nTo cite, type citation(\"%s\")",
    pkgname, utils::packageDescription(pkgname)$Version, pkgname
  ))
  packageStartupMessage(
    "Please report improvements and bugs to: ",
    "https://github.com/sergioventurini/XOMultinom/issues"
  )
}

.onUnload <- function(libpath) {
  library.dynam.unload("XOMultinom", libpath)
}