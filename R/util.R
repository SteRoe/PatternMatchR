#' starts the app
#' @description function to start the app
#' @export
PatternMatchRApp <- function() {
  shiny::shinyApp(ui, server)
}
