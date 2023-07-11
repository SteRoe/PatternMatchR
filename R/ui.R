generate_ui <- function() {
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
      shiny::tabsetPanel(
        shiny::tabPanel(
          "PatternMatchR",
          shinyBS::bsCollapse(
            id = "clpMain",
            open = c("session", "folders", "merge", "reduce traits by clustering", "HeatMap"),
            multiple = TRUE,
            shinyBS::bsCollapsePanel(
              "session",
              shiny::fluidRow(
                shiny::tags$table(
                  style = "width: 100%",
                  shiny::tags$tr(
                    shiny::tags$td(
                      align = "center",
                      shiny::actionButton(inputId = "btnExportCombinedData", label = "Export session data to <CombinedHM.RDS>"),
                      shinyFiles::shinySaveButton(id = "save", label = "Export session data",
                                                  title = "Save session data as file",
                                                  filetype = list(RDS = "RDS", text = "txt"),
                                                  viewtype = "detail")
                    ),
                    shiny::tags$td(
                      align = "center",
                      shiny::actionButton(inputId = "btnImportCombinedData", label = "Import session data from <CombinedHM.RDS>"),
                      shinyFiles::shinyFilesButton(id = "file", label = "Import session data",
                                                   title = "Load session data from file",
                                                   multiple = FALSE, viewtype = "detail")
                    )
                  )
                )
              )
            ),
            shinyBS::bsCollapsePanel(
              "folders",
              shiny::fluidRow(
                shiny::column(
                  id = "colRed",
                  width = 4,
                  "red traits",
                  DT::dataTableOutput("trait1DirList"),
                  shiny::actionButton("btnLoadDir1", label = "Load trait 1 (red) dir")
                ),
                shiny::column(
                  id = "colGreen",
                  width = 4,
                  "green traits",
                  DT::dataTableOutput("trait2DirList"),
                  shiny::actionButton("btnLoadDir2", label = "Load trait 2 (green) dir")
                ),
                shiny::column(
                  id = "colBlue",
                  width = 4,
                  "blue traits",
                  DT::dataTableOutput("trait3DirList"),
                  shiny::actionButton("btnLoadDir3", label = "Load trait 3 (blue) dir")
                )
              ),
              shiny::actionButton("btnLoadDirAll", label = "Step 1: Load all trait dirs"),
              shiny::verbatimTextOutput("txtLoadOut", placeholder = TRUE),
              shiny::fluidRow(
              ),
              shiny::actionButton("btnMerge", label = "Step 2: Merge data from all folders"),
              shiny::verbatimTextOutput("txtMergeOut", placeholder = TRUE)
            ),
            shinyBS::bsCollapsePanel(
              "count borders",
              shiny::tabsetPanel(
                shiny::tabPanel("P_VAL border vs. # of probes",
                  shiny::fluidRow(
                    shiny::actionButton("btnCountP_ValProbes", label = "Count Probes for p-values (may take a long time)"),
                    shiny::actionButton("btnCountProbesP_ValParallel", label = "Count Probes for p-values parallel (may take less time)")
                  ),
                  shiny::fluidRow(
                    #shiny::column(
                    #width = 5,
                    DT::dataTableOutput("DTP_VALborder"),
                    plotly::plotlyOutput("plotDendrogramP_VALborder", height = "80%")
                    #)
                  )
                ),
                shiny::tabPanel("Delta Methylation border vs. # of probes",
                  shiny::fluidRow(
                    shiny::actionButton("btnCountProbesDeltaMethParallel", label = "Count Probes for delta methylation parallel")
                  ),
                  shiny::fluidRow(
                    #shiny::column(
                    #width = 5,
                    DT::dataTableOutput("DTDMborder"),
                    plotly::plotlyOutput("plotDendrogramDMborder", height = "80%")
                    #)
                  )
                ),
                shiny::tabPanel("n border vs. # of probes",
                  shiny::fluidRow(
                    shiny::actionButton("btnCountProbesNParallel", label = "Count Probes for n parallel")
                  ),
                  shiny::fluidRow(
                    #shiny::column(
                    #width = 5,
                    DT::dataTableOutput("DTNborder"),
                    plotly::plotlyOutput("plotDendrogramNborder", height = "80%")
                    #)
                  )
                )
              )
            ),

            shinyBS::bsCollapsePanel(
              "merge",
              shiny::fluidRow(
                shiny::column(
                  width = 4,
                  shiny::sliderInput(
                    "sldP_Val",
                    "maximum (left slider) and minimum (right slider) p-val, 5e-x",
                    min = 3,
                    max = 200,
                    step = -1,
                    value = c(3, 199)
                  )
                ),
                shiny::column(
                  width = 4,
                  shiny::sliderInput(
                    "sldDM",
                    "maximum (left slider) and minimum (right slider) delta methylation",
                    min = 0,
                    max = 200,
                    step = .01,
                    value = c(3, 199)
                  )
                ),
                shiny::column(
                  width = 4,
                  shiny::sliderInput(
                    "sldN",
                    "maximum (left slider) and minimum (right slider) n",
                    min = 0,
                    max = 200,
                    step = 1,
                    value = c(3, 199)
                  )
                )
              ),
              shiny::actionButton("btnReduce", label = "Step 3: Reduce data (omit CpGs) by p-value"),
              shiny::verbatimTextOutput("txtPReduceOut", placeholder = TRUE)
            ),

            shinyBS::bsCollapsePanel(
              "reduce traits by clustering",
              shiny::fluidRow(
                shiny::sliderInput(
                  "sldNumClusters",
                  "number of clusters (omit traits)",
                  min = 1,
                  max = 1,
                  step = 1,
                  value = 1 # value = c(1, 10)
                ),
                shiny::tabsetPanel(
                  shiny::tabPanel(
                    "Dendrogram",
                    shiny::verbatimTextOutput("txtDendrogramTraitsLong", placeholder = TRUE),
                    shiny::plotOutput("plotDendrogramTraitsLong")
                  ),
                  shiny::tabPanel(
                    "Clustergram",
                    shiny::verbatimTextOutput("txtClustergramTraitsLong", placeholder = TRUE),
                    shiny::plotOutput("plotClustergramTraitsLong")
                  ),
                  shiny::tabPanel(
                    "DT Cluster Medoids",
                    shiny::verbatimTextOutput("txtDTTraitsMedoids", placeholder = TRUE),
                    DT::dataTableOutput("DTTraitsMedoids")
                  ),
                  shiny::tabPanel(
                    "DT Cluster Assignment",
                    shiny::verbatimTextOutput("txtDTTraitsClusters", placeholder = TRUE),
                    DT::dataTableOutput("DTTraitsClusters")
                  )
                ),
                shiny::actionButton("btnCluster", label = "Step 4: Cluster & reduce trait data"),
                shiny::verbatimTextOutput("txtClusterOut", placeholder = TRUE)
              )
            ),

            shinyBS::bsCollapsePanel(
              "HeatMap",
              shiny::fluidRow(
                "Original data",
                shiny::tabsetPanel(
                  shiny::tabPanel(
                    "HeatMap",
                    shiny::tabsetPanel(
                      shiny::tabPanel(
                        "Histogram",
                        "Histogram of all p-values in heatmap",
                        plotly::plotlyOutput("histP_Val", inline = TRUE)
                      ),
                      shiny::tabPanel(
                        "HeatMap",
                        shiny::actionButton("plotCombinedHM", label = "Step 5: Plot Heatmap"),
                        shiny::verbatimTextOutput("txtHMDescription", placeholder = TRUE),
                        InteractiveComplexHeatmap::InteractiveComplexHeatmapOutput(
                          "heatmap_1",
                          height1 = 2000,
                          width1 = 1500,
                          height2 = 2000,
                          width2 = 1500,
                          inline = FALSE,
                          output_ui = shiny::htmlOutput("info1")
                        )
                      ),
                      shiny::tabPanel(
                        "SPLOM",
                        shiny::selectizeInput("markingVar",
                          label = "select variable for color marking",
                          choices = NULL,
                          width = "100%"
                        ),
                        "SPLOM from selected area in heatmap",
                        plotly::plotlyOutput("SPLOM",
                                             height = 2000, #height = "100%",
                                             width = 1500, #width = "100%",
                                             inline = TRUE)
                      )
                    )
                  ),
                  shiny::tabPanel(
                    "Table P_VAL",
                    # DT tableout
                    DT::dataTableOutput("DTP_VAL")
                  ),
                  shiny::tabPanel(
                    "Table Delta Methylation",
                    # DT tableout
                    DT::dataTableOutput("DTDM")
                  ),
                  shiny::tabPanel(
                    "Table N",
                    # DT tableout
                    DT::dataTableOutput("DTN")
                  ),
                  shiny::tabPanel(
                    "Dendrogram Probes",
                    plotly::plotlyOutput("plotDendrogramProbes", height = "80%")
                  ),
                  shiny::tabPanel(
                    "Annotated Table Probes",
                    DT::dataTableOutput("DTProbes")
                  ),
                  shiny::tabPanel(
                    "Dendrogram Traits",
                    plotly::plotlyOutput("plotDendrogramTraits", height = "80%")
                  ),
                  shiny::tabPanel(
                    "Table Traits",
                    DT::dataTableOutput("DTTraits")
                  )
                )
              )
            )
          )
        ),
        shiny::tabPanel(
          "Settings",
          shiny::fluidRow(
            shiny::column(
              width = 12,
            shiny::checkboxInput("chkDebug", "DEBUG mode", TRUE),
            shiny::actionButton(inputId = "btnBrowser", label = "Break to Browser()"),
            shiny::actionButton(inputId = "btnDebug", label = "Debug"),
            shiny::verbatimTextOutput("txtDebugOut", placeholder = TRUE)
            )
          ),
          shinyBS::bsCollapse(
            id = "clpSettings",
            multiple = TRUE,
            shinyBS::bsCollapsePanel(
              "methylation file",
              shiny::fluidRow(
                shiny::textInput(inputId = "inpDNAmFile", label = "File with DNAm", value = "FileName", width = NULL, placeholder = TRUE),
              )
            ),
            shinyBS::bsCollapsePanel(
              "minimum cases for trait",
              shiny::fluidRow(
                shiny::textInput(inputId = "txtCases", label = "minimum # cases for each trait", value = 100, width = NULL, placeholder = TRUE),
              )
            ),
            shinyBS::bsCollapsePanel(
              "p-val settings",
              shiny::fluidRow(
                shiny::column(
                  width = 2,
                  shiny::textInput(inputId = "txtMaxProbes", label = "maximum Probes",
                                   value = 500000, width = NULL, placeholder = TRUE),
                  # sqrt(2^31)=46340
                ),
                shiny::column(
                  width = 2,
                  shiny::textInput(inputId = "txtMaxClassesProbes",
                    label = "maximum colored Classes Probes (in dendrogram)",
                    value = 10,
                    placeholder = TRUE
                  ),
                  shiny::textInput(inputId = "txtMaxClassesTraits",
                    label = "maximum colored Classes Traits (in dendrogram)",
                    value = 10,
                    placeholder = TRUE
                  )
                ),
                shiny::column(
                  width = 4,
                  shiny::tags$html("n"),
                  shiny::verbatimTextOutput("txtResultingN", placeholder = TRUE),
                  shiny::tags$html("minimum p-value"),
                  shiny::verbatimTextOutput("txtMinP_Val", placeholder = TRUE),
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
            )
          )
        )
      )
    )
  )
}
