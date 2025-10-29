#' Launch seneR Shiny App
#' @export
runApp <- function() {
  appDir <- system.file("shiny", package = "seneR")
  if (appDir == "") stop("Shiny app not found in this package.")
  shiny::runApp(appDir, display.mode = "normal")
}