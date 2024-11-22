#' shiny_bluntPred
#' Launch an interactive web interface.
#'
#' @param ... other params sent to shiny::runApp().
#' @return Nothing
#' @import shiny
#' @import shinydashboard
#' @import DT
#' @export
shiny_bluntPred <- function(...) {
  shiny::runApp(appDir=system.file("shiny_bluntPred", package="bluntPred"), ...)
}
