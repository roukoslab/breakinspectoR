#' shiny_bluntPred
#' Launch an interactive web interface.
#'
#' @param ... other params sent to shiny::runApp().
#' @return Nothing
#' @export
shiny_bluntPred <- function(...) {
  message("Loading dependencies...")
  if (requireNamespace("shiny"         , quietly=TRUE) &
      requireNamespace("shinydashboard", quietly=TRUE) &
      requireNamespace("h2o"     , quietly=TRUE)) {
    message("Launching app")
    shiny::runApp(appDir=system.file("shiny_bluntPred", package="breakinspectoR"), ...)
  } else {
    stop("The following CRAN packages are required to run the shiny app:\n",
         "  - shiny\n",
         "  - shinydashboard\n",
         "  - h2o (version *exactly* 3.36.1.2)\n")
  }
}
