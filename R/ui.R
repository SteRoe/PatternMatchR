generate_ui <- function() {
  shiny::shinyUI(
    shiny::fluidPage(
      # waiter::useWaiter(), # dependencies
      # waiter::transparent(alpha = 0.1),
      # waiter::waiterShowOnLoad(), # shows before anything else
      shinyjs::useShinyjs(),
      shiny::tags$script('
              $(document).on("keyup", function (e) {
              Shiny.onInputChange("keypressed", e.which);
              });
              '),
      shinyjs::inlineCSS(list(.red = "background: lightcoral",
                     .green = "background: lightgreen",
                     .blue = "background: lightblue")),

      #smaller font for preformatted text
      shiny::tags$head(shiny::tags$style(
        shiny::HTML("      pre, table.table {        font-size: smaller;      }    ")
      )),
      "Hostname / PID",
      shiny::verbatimTextOutput(outputId = "Sys.PID", placeholder = TRUE),
      shiny::tabsetPanel(
        shiny::tabPanel(
          "PatternMatchR",
          #shiny::tagList(
          shinyBS::bsCollapse(
          id = "clpPreprocess",
          open = c("Folders", "Merge", "Reduce Data"),
            multiple = TRUE,
            shinyBS::bsCollapsePanel(
              "Folders",
              shiny::fluidRow(
                shiny::column(
                  id = "colRed",
                  width = 4,
                  "red traits",
                  DT::DTOutput("trait1DirList"), #DT::dataTableOutput("trait1DirList"),
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
              shiny::verbatimTextOutput("txtLoadOut", placeholder = TRUE)
            ), #end bsCollapsePanel
            shinyBS::bsCollapsePanel(
              "Merge",
              shiny::fluidRow(
                shiny::actionButton("btnMerge", label = "Step 2: Merge data from all folders"),
                shiny::verbatimTextOutput("txtMergeOut", placeholder = TRUE)
              )
            ), #end bsCollapsePanel
            shinyBS::bsCollapsePanel(
              "Count Borders",
              shiny::tabsetPanel(
                shiny::tabPanel("P_VAL border vs. # of probes",
                  shiny::fluidRow(
                    # shinyjs::disabled(shiny::actionButton("btnCountP_ValProbes", label = "Count Probes for p-values (may take a long time)")),
                    shinyjs::disabled(shiny::actionButton("btnCountProbesP_ValParallel", label = "Count Probes for p-values parallel"))
                  ),
                  shiny::fluidRow(
                    shiny::tabsetPanel(
                      shiny::tabPanel("Dendrogram",
                        plotly::plotlyOutput("plotDendrogramP_VALborder", height = "80%")
                      ),
                      shiny::tabPanel("Table",
                        DT::dataTableOutput("DTP_VALborder")
                      )
                    )
                  )
                ),
                shiny::tabPanel("Delta Methylation border vs. # of probes",
                  shiny::fluidRow(
                    shinyjs::disabled(shiny::actionButton("btnCountProbesDeltaMethParallel", label = "Count Probes for delta methylation parallel"))
                  ),
                  shiny::fluidRow(
                    shiny::tabsetPanel(
                      shiny::tabPanel("Dendrogram",
                        plotly::plotlyOutput("plotDendrogramDMborder", height = "80%")
                      ),
                      shiny::tabPanel("Table",
                        DT::dataTableOutput("DTDMborder")
                      )
                    )
                  )
                ),
                shiny::tabPanel("n border vs. # of probes",
                  shiny::fluidRow(
                    shinyjs::disabled(shiny::actionButton("btnCountProbesNParallel", label = "Count Probes for n parallel"))
                  ),
                  shiny::fluidRow(
                    shiny::tabsetPanel(
                      shiny::tabPanel("Dendrogram",
                        plotly::plotlyOutput("plotDendrogramNborder", height = "80%")
                      ),
                      shiny::tabPanel("Table",
                        DT::dataTableOutput("DTNborder")
                      )
                    )
                  )
                )
              )
            ), #end bsCollapsePanel
            shinyBS::bsCollapsePanel(
              "Reduce Data",
#              shiny::verbatimTextOutput("txtTest", placeholder = TRUE),
              shiny::fluidRow(
                shiny::column(
                  width = 4,
                  shinyjs::disabled(shiny::sliderInput(
                    "sldP_Val",
                    "maximum (left slider) and minimum (right slider) p-val, 5e-x",
                    min = 0,
                    max = 0,
                    step = 0, #-1, #1
                    value = c(0, 0)
                  ))
                ),
                shiny::column(
                  width = 4,
                  shinyjs::disabled(shiny::sliderInput(
                    "sldDM",
                    "minimum (left slider) and maximum (right slider) delta methylation",
                    min = 0,
                    max = 0,
                    step = 0, #.01,
                    value = c(0, 0)
                  ))
                ),
                shiny::column(
                  width = 4,
                  shinyjs::disabled(shiny::sliderInput(
                    "sldN",
                    "minimum (left slider) and maximum (right slider) n",
                    min = 0,
                    max = 0,
                    step = 0,
                    value = c(0, 0)
                  ))
                )
              ),
              shinyjs::disabled(shiny::actionButton("btnReduce", label = "Step 3: Reduce data (omit CpGs) by applying thresholds for p-value, DM or n limit")),
              shiny::verbatimTextOutput("txtPReduceOut", placeholder = TRUE)
            ) #end bsCollapsePanel "Reduce Data"
          ), #end bsCollapse "clpPreprocess"
          shinyBS::bsCollapse(
            id = "clpClustering",
            open = c("Clustering"),
            multiple = TRUE,
            shinyBS::bsCollapsePanel(
              "Clustering",
              shiny::fluidRow(
                shiny::column(
                  width = 4,
                  shinyjs::disabled(shiny::sliderInput(
                    "sldNumClusters",
                    "number of clusters (omit traits based on p-value)",
                    min = 0,
                    max = 0,
                    step = 0,
                    value = 0 # value = c(1, 10)
                  ))
                )
              ), #end fluidRow
              shiny::fluidRow(
                shiny::column(
                  width = 4,
                  shinyjs::disabled(shiny::sliderInput(
                    "sld_NumNeighbours",
                    "distance to look for neighbours",
                    min = 10,
                    max = 10000,
                    step = 10,
                    value = 10
                  ))
                )
              ), #end fluidRow
              shiny::tabsetPanel(
                shiny::tabPanel(
                  "Clustering based on p-value",
                  Clustering_UI("P_Val")
                ),
                shiny::tabPanel(
                  "Clustering based on log(FC)",
                  Clustering_UI("LogFC")
                )
              ) #end tabSetPanel
            ) #end bsCollapsePanel "Clustering"
          ), #end bsCollapse "clpClustering"
          shinyBS::bsCollapse(
            id = "clpHeatmap",
            shinyBS::bsCollapsePanel(
              "Heatmaps",
              shiny::fluidRow(
                shiny::tabsetPanel(
                  shiny::tabPanel( ##full non-modified data p-val start
                    "Non-modified Data (p-val)",
                    shiny::tabsetPanel(
                      shiny::tabPanel(
                        "HeatMap P_Val",
                        shinyjs::disabled(
                          shiny::actionButton("btnPlotCombinedHM_P_Val", label = "Step 5a: Plot Heatmap P_Val")
                        ),
                        shinyjs::disabled(
                          shiny::verbatimTextOutput("txtHMDescription_P_Val", placeholder = TRUE)
                        ),
                        InteractiveComplexHeatmap::InteractiveComplexHeatmapOutput(
                          "Heatmap_P_Val",
                          height1 = '95vh', #1200,
                          width1 = 950,
                          height2 = '95vh', #1200,
                          width2 = 950,
                          inline = FALSE
                        )
                      ),
                      shiny::tabPanel(
                        "HeatMap P_Val Details",
                        HeatMap_UI("HeatMap_Full_DetailsPval")
                      ),
                      shiny::tabPanel(
                        "Dendrogram Probes",
                        plotly::plotlyOutput("traitReducedPlotDendrogramProbesPval", height = "80%")
                      ),
                      shiny::tabPanel(
                        "Annotated Table Probes",
                        DT::dataTableOutput("traitReducedDTProbesPval")
                      ),
                      shiny::tabPanel(
                        "Dendrogram Traits",
                        plotly::plotlyOutput("traitReducedPlotDendrogramTraitsPval", height = "80%")
                      ),
                      shiny::tabPanel(
                        "Table Traits",
                        DT::dataTableOutput("traitReducedDTTraitsPval")
                      ),
                      shiny::tabPanel(
                        "Histogram",
                        "Histogram of all p-values in full heatmap (number of p-values = number probes * number traits)",
                        plotly::plotlyOutput("traitReducedHistP_Val", inline = TRUE)
                      ) #end tabPanel
                    ) #end tabSetPanel
                  ), #end tabPanel ##full non-modified data p-val end
                  shiny::tabPanel( ##full non-modified data log(FC) start
                    "Non-modified Data (log(FC))",
                    shiny::tabsetPanel(
                      shiny::tabPanel(
                        "HeatMap log(FC)",
                        shinyjs::disabled(
                          shiny::actionButton("btnPlotCombinedHM_LogFC", label = "Step 5b: Plot Heatmap log(FC)")
                        ),
                        shinyjs::disabled(
                          shiny::verbatimTextOutput("txtHMDescription_LogFC", placeholder = TRUE)
                        ),
                        InteractiveComplexHeatmap::InteractiveComplexHeatmapOutput(
                          "Heatmap_LogFC",
                          height1 = '95vh', #1200,
                          width1 = 950,
                          height2 = '95vh', #1200,
                          width2 = 950,
                          inline = FALSE
                        )
                      ),
                      shiny::tabPanel(
                        "HeatMap log(FC) Details",
                        HeatMap_UI("HeatMap_Full_DetailsLogFC")
                      ),
                      shiny::tabPanel(
                        "Dendrogram Probes",
                        plotly::plotlyOutput("traitReducedPlotDendrogramProbesLogFC", height = "80%")
                      ),
                      shiny::tabPanel(
                        "Annotated Table Probes",
                        DT::dataTableOutput("traitReducedDTProbesLogFC")
                      ),
                      shiny::tabPanel(
                        "Dendrogram Traits",
                        plotly::plotlyOutput("traitReducedPlotDendrogramTraitsLogFC", height = "80%")
                      ),
                      shiny::tabPanel(
                        "Table Traits",
                        DT::dataTableOutput("traitReducedDTTraitsLogFC")
                      ),
                      shiny::tabPanel(
                        "Histogram",
                        "Histogram of all log(FC) in full heatmap (number of log(FC) = number probes * number traits)",
                        plotly::plotlyOutput("traitReducedHistLogFC", inline = TRUE)
                      )
                    )
                  ), ## full non-modified data log(FC) end
                  shiny::tabPanel(
                    ##full DW data start
                    "Distance weighted Data (negative p-values due to distance weighting) - experimental",
                    shiny::tabsetPanel(
                      shiny::tabPanel(
                        "HeatMap P_Val",
                        shiny::actionButton("btnPlotCombinedDWHM_P_Val", label = "Step 5c: Plot distance weighted Heatmap P_Val"),
                        shiny::verbatimTextOutput("txtDWHMDescription_P_Val", placeholder = TRUE),
                        # shiny::column(
                        #   width = 6,
                        #   shinyjs::disabled(shiny::numericInput(inputId = "numDWHMHSize", label = "Width", value = 4000, min = 1000, max = 10000))
                        # ),
                        # shiny::column(
                        #   width = 6,
                        #   shinyjs::disabled(shiny::numericInput(inputId = "numDWHMVSize", label = "Height", value = 4000, min = 1000, max = 10000))
                        # ),
                        InteractiveComplexHeatmap::InteractiveComplexHeatmapOutput(
                          "DWHeatmap_P_Val",
                          height1 = '95vh',
                          width1 = 950,
                          height2 = '95vh',
                          width2 = 950,
                          inline = FALSE
                        )
                      ),
                      shiny::tabPanel(
                        "Table P_Val",
                        "Table of p-value; clustering order comes from clustering of p-values.",
                        DT::dataTableOutput("fullDWDTP_VALPval")
                      ),
                      shiny::tabPanel(
                        "Table Delta Methylation",
                        "Table of delta methylation; clustering order comes from clustering of p-values.",
                        DT::dataTableOutput("fullDWDTDMPval")
                      ),
                      shiny::tabPanel(
                        "Table Delta Methylation log(FC)",
                        "Table of log fold change(delta methylation); clustering order comes from clustering of p-values.",
                        DT::dataTableOutput("fullDWDTlogFCPval")
                      ),
                      shiny::tabPanel(
                        "Table N",
                        "Table of n; clustering order comes from clustering of p-values.",
                        DT::dataTableOutput("fullDWDTNPval")
                      ),
                      shiny::tabPanel(
                        "Dendrogram Probes",
                        plotly::plotlyOutput("fullDWPlotDendrogramProbesPval", height = "80%")
                      ),
                      shiny::tabPanel(
                        "Annotated Table Probes",
                        DT::dataTableOutput("fullDWDTProbesPval")
                      ),
                      shiny::tabPanel(
                        "Dendrogram Traits",
                        plotly::plotlyOutput("fullDWPlotDendrogramTraitsPval", height = "80%")
                      ),
                      shiny::tabPanel(
                        "Table Traits",
                        DT::dataTableOutput("fullDWDTTraitsPval")
                      ),
                      shiny::tabPanel(
                        "Histogram P_Val",
                        "Histogram of all p-values in condensed heatmap (number of p-values = number probes * number traits)",
                        plotly::plotlyOutput("fullDWHistP_Val", inline = TRUE)
                      )
                    )
                  ) ##full DW data end
                )
              )
            ) #end bsCollapsePanel "Heatmap"
            # )
          ),
          shinyBS::bsCollapse(
          id = "clpGlobals",
            shinyBS::bsCollapsePanel(
              "Global Search",
              Search_Full_UI("Search")
            ), #end bsCollapsePanel
            shinyBS::bsCollapsePanel(
              "Global Selection",
              GlobalSelection_UI("GlobalSelection")
            ), #end bsCollapsePanel
            shinyBS::bsCollapsePanel(
              "Selection Visualization",
              shiny::tabsetPanel(
                shiny::tabPanel(
                  "VolcanoPlot Delta Methylation log(FC)",
                  "P-values and log fold change (delta methylation)",
                  VolcanoPlot_UI("VolcanoPlot")
                ),
                shiny::tabPanel(
                  "Near Range Methylation Profile Near Range DMR PC Plot",
                  "Near Range DMR PC Plot",
                  PCPlot_UI("PCPlot")
                )
              )
            ) #end bsCollapsePanel "Selection Visualization"
          ) #end bsCollapse "clpGlobals"
        ), #end tabPanel "PatternMatchR"
        shiny::tabPanel(
          "Settings",
          shiny::fluidRow(
            shiny::column(
              width = 12,
              shiny::checkboxInput("chkDebug", "DEBUG mode", FALSE),
              shiny::actionButton(inputId = "btnBrowser", label = "Break to Browser()"),
              shiny::actionButton(inputId = "btnDebug", label = "Debug"),
              shiny::verbatimTextOutput("txtDebugOut", placeholder = TRUE)
            ) #end column
          ), #end fluidRow
          shinyBS::bsCollapse(
            id = "clpSettings",
            multiple = TRUE,
            shinyBS::bsCollapsePanel(
              "methylation file",
              shiny::fluidRow(
                shiny::textInput(inputId = "inpDNAmFile", label = "File with DNAm", value = "FileName", width = NULL, placeholder = TRUE)
              ) #end fluidRow
            ), #end bsCollapsePanel "methylation file"
            shinyBS::bsCollapsePanel(
              "minimum cases for trait",
              shiny::fluidRow(
                shiny::textInput(inputId = "txtCases", label = "minimum # cases for each trait", value = 100, width = NULL, placeholder = TRUE)
              ) #end fluidRow
            ), #end bsCollapsePanel "minimum cases for trait"
            shinyBS::bsCollapsePanel(
              "p-val settings",
              shiny::fluidRow(
                shiny::column(
                  width = 2,
                  shiny::textInput(inputId = "txtMaxProbes", label = "maximum Probes",
                                   value = 500000, width = NULL, placeholder = TRUE)
                ), #end column
                shiny::column(
                  width = 2,
                  shiny::textInput(inputId = "txtMaxClassesProbes",
                    label = "maximum colored Classes Probes (in dendrogram)",
                    value = 7,
                    placeholder = TRUE
                  ),
                  shiny::textInput(inputId = "txtMaxClassesTraits",
                    label = "maximum colored Classes Traits (in dendrogram)",
                    value = 7,
                    placeholder = TRUE
                  )
                ) #end column
                # shiny::column(
                #   width = 4,
                #   shiny::tags$html("n"),
                #   shiny::verbatimTextOutput("txtResultingN", placeholder = TRUE),
                #   shiny::tags$html("minimum p-value"),
                #   shiny::verbatimTextOutput("txtMinP_Val", placeholder = TRUE),
                #   InteractiveComplexHeatmap::InteractiveComplexHeatmapOutput(
                #     "heatmap_2",
                #     height1 = 10,
                #     width1 = 10,
                #     height2 = 10,
                #     width2 = 10,
                #     inline = FALSE
                #   )
                # )
              ) #end fluidRow
            ) #end bsCollapsePanel "p-val settings"
          ) #end bsCollapse "clpSettings"
        ) #end tabPanel "Settings"
      ) # end tabSetPanel
    ) # end fluidPage
  ) # end shinyUI
} #end function
