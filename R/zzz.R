.onLoad <- function(libname, pkgname) {
  PROJ_DATA <- Sys.getenv("PROJ_DATA")
  proj_data <- "/usr/share/proj"
  Sys.setenv(PROJ_DATA = proj_data)
  if (!PROJ_DATA == proj_data) {
    message(sprintf(".onLoad %s PROJ_DATA has been changed from '%s' to '%s'", pkgname, PROJ_DATA, proj_data))
  }
}
