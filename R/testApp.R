ui <- shiny::shinyUI(
  shiny::fluidPage(
    InteractiveComplexHeatmap::InteractiveComplexHeatmapOutput(
      "heatmap_2",
      height1 = 100,
      width1 = 100,
      height2 = 100,
      width2 = 100,
      inline = FALSE
    )
  )
)

server <- function(input, output, session) {
  base::print(paste0(Sys.time(), " ver: 0.2.5."))
  #draw empty HM, without the github version won't work. Why?
  InteractiveComplexHeatmap::makeInteractiveComplexHeatmap(
    input = input,
    output = output,
    session = session,
    ht_list = emptyHM(),
    heatmap_id = "heatmap_2",
    show_layer_fun = FALSE,
    click_action = NULL,
    brush_action = NULL,
    hover_action = NULL
  )
}

#' emptyHM
#' creates an empty heatmap
#' @return empty heatmap
# examples emptyHM()
emptyHM <- function() {
  mat <- base::matrix(c(1, 2, 3, 4), nrow = 2, ncol = 2)
  ht <- ComplexHeatmap::Heatmap(mat)
  if (grDevices::dev.cur() > 1) {
    grDevices::dev.off()
  }
  grDevices::pdf(NULL)
  grDevices::dev.off()
  ht <- ComplexHeatmap::draw(ht)
  #return(ht)
}

#' starts the app
#' @description function to start the app
#' @export
PatternMatchRApp <- function() {
  shiny::shinyApp(ui, server)
}
