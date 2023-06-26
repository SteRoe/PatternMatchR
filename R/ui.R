ui <- shiny::shinyUI(
  shiny::fluidPage(
    shinyjs::useShinyjs(),

    shinyjs::inlineCSS(list(.red = "background: lightcoral",
                   .green = "background: lightgreen",
                   .blue = "background: lightblue")),

    #smaller font for preformatted text
    shiny::tags$head(shiny::tags$style(
      shiny::HTML("      pre, table.table {        font-size: smaller;      }    ")
    )),

    InteractiveComplexHeatmap::InteractiveComplexHeatmapOutput(
      "heatmap_2",
      height1 = 10,
      width1 = 10,
      height2 = 10,
      width2 = 10,
      inline = FALSE
    )
  )
)
