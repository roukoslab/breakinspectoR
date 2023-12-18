#' shiny_BTmotif
#' Launch an interactive web interface.
#'
#' @param ... other params sent to shiny::runApp().
#' @return Nothing
#' @export
shiny_BTmotif <- function(...) {
  message("Loading dependencies...")
  if (requireNamespace("shiny"         , quietly=TRUE) &
      requireNamespace("shinydashboard", quietly=TRUE) &
      requireNamespace("h2o"           , quietly=TRUE) &
      requireNamespace("ggplot2"       , quietly=TRUE) &
      requireNamespace("ggseqlogo"     , quietly=TRUE)) {
    message("Launching app")
    shiny::runApp(appDir=system.file("shiny_BTmotif", package="breakinspectoR"), ...)
  } else {
    stop("The following CRAN packages are required to run the shiny app:\n",
         "  - shiny\n",
         "  - shinydashboard\n",
         "  - h2o\n",
         "  - ggplot2\n",
         "  - ggseqlogo\n")
  }
}
