generate_ui <- function() {
  shiny::shinyUI(
    shiny::fluidPage(
      shiny::useBusyIndicators(),
      shiny::busyIndicatorOptions(spinner_type = "ring"),
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
          open = c("Folders", "Merge Data", "Reduce Data"),
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
              "Merge Data",
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
              ReduceData_UI("Reduce")
            ) #end bsCollapsePanel "Reduce Data"
          ), #end bsCollapse "clpPreprocess"
          shinyBS::bsCollapse(
            id = "clpClustering",
            open = c("Clustering for Traits", "Clustering for Probes and Filtering for Neighbours"),
            multiple = TRUE,
            shinyBS::bsCollapsePanel(
              "Clustering for Traits",
              Clustering_UI("Traits"),
              shiny::tabsetPanel(
                shiny::tabPanel(
                  "Clustering based on p-value",
                  ClusteringTraits_UI("PVal")
                ),
                shiny::tabPanel(
                  "Clustering based on log(FC)",
                  ClusteringTraits_UI("LogFC")
                )
              ) #end tabSetPanel
            ), #end bsCollapsePanel "Clustering Traits"
            shinyBS::bsCollapsePanel(
              "Clustering for Probes and Filtering for Neighbours",
              ClusteringProbesGeneral_UI("Probes"),
              shiny::tabsetPanel(
                shiny::tabPanel(
                  "Clustering based on p-value",
                  ClusteringProbes_UI("PVal")
                ),
                shiny::tabPanel(
                  "Clustering based on log(FC)",
                  ClusteringProbes_UI("LogFC")
                )
              ) #end tabSetPanel
            )
          ), #end bsCollapse "clpClusteringTraits"

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
                          shiny::actionButton("btnPlotCombinedHM_P_Val", label = "Step 6a: Plot Heatmap (p-val)")
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
                      ), #end tabPanel
                      shiny::tabPanel(
                        "HeatMap P_Val Details",
                        HeatMap_UI("PVal")
                      ) #end tabPanel
                    ) #end tabSetPanel
                  ), #end tabPanel ##full non-modified data p-val end
                  shiny::tabPanel(
                    "Data (p-val) w/o Gaps",
                    shiny::tabsetPanel(
                      shiny::tabPanel(
                        "HeatMap p-val",
                        shinyjs::disabled(
                          shiny::actionButton("btnPlotCombinedHM_P_ValWOGap", label = "Step 6b: Plot Heatmap (p-val) w/o Gaps")
                        ),
                        shinyjs::disabled(
                          shiny::verbatimTextOutput("txtHMDescription_P_ValWOGap", placeholder = TRUE)
                        ),
                        InteractiveComplexHeatmap::InteractiveComplexHeatmapOutput(
                          "Heatmap_P_ValWOGap",
                          height1 = '95vh', #1200,
                          width1 = 950,
                          height2 = '95vh', #1200,
                          width2 = 950,
                          inline = FALSE
                        )
                      ), #end tabPanel
                      shiny::tabPanel(
                        "HeatMap P_Val Details",
                        HeatMap_UI("PValWOGap")
                      ) #end tabPanel
                    ) #end tabSetPanel
                  ), #end tabPanel
                  shiny::tabPanel( ##full non-modified data log(FC) start
                    "Non-modified Data (log(FC))",
                    shiny::tabsetPanel(
                      shiny::tabPanel(
                        "HeatMap log(FC)",
                        shinyjs::disabled(
                          shiny::actionButton("btnPlotCombinedHM_LogFC", label = "Step 6c: Plot Heatmap log(FC)")
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
                        HeatMap_UI("LogFC")
                      ) #end tabPanel
                    ) #end tabSetPanel
                  ), ## full non-modified data log(FC) end
                  shiny::tabPanel( ##full non-modified data log(FC) start
                    "Data (log(FC)) w/o Gaps",
                    shiny::tabsetPanel(
                      shiny::tabPanel(
                        "HeatMap log(FC)",
                        shinyjs::disabled(
                          shiny::actionButton("btnPlotCombinedHM_LogFCWOGap", label = "Step 6d: Plot Heatmap log(FC) w/o Gaps")
                        ),
                        shinyjs::disabled(
                          shiny::verbatimTextOutput("txtHMDescription_LogFCWOGap", placeholder = TRUE)
                        ),
                        InteractiveComplexHeatmap::InteractiveComplexHeatmapOutput(
                          "Heatmap_LogFCWOGap",
                          height1 = '95vh', #1200,
                          width1 = 950,
                          height2 = '95vh', #1200,
                          width2 = 950,
                          inline = FALSE
                        )
                      ),
                      shiny::tabPanel(
                        "HeatMap log(FC) Details",
                        HeatMap_UI("LogFCWOGap")
                      ) #end tabPanel
                    ) #end tabSetPanel
                  )
                  # shiny::tabPanel(
                  #   ##full DW data start
                  #   "Distance weighted Data (negative p-values due to distance weighting) - experimental",
                  #   shiny::tabsetPanel(
                  #     shiny::tabPanel(
                  #       "HeatMap P_Val",
                  #       shiny::actionButton("btnPlotCombinedDWHM_P_Val", label = "Step 5c: Plot distance weighted Heatmap P_Val"),
                  #       shiny::verbatimTextOutput("txtDWHMDescription_P_Val", placeholder = TRUE),
                  #       # shiny::column(
                  #       #   width = 6,
                  #       #   shinyjs::disabled(shiny::numericInput(inputId = "numDWHMHSize", label = "Width", value = 4000, min = 1000, max = 10000))
                  #       # ),
                  #       # shiny::column(
                  #       #   width = 6,
                  #       #   shinyjs::disabled(shiny::numericInput(inputId = "numDWHMVSize", label = "Height", value = 4000, min = 1000, max = 10000))
                  #       # ),
                  #       InteractiveComplexHeatmap::InteractiveComplexHeatmapOutput(
                  #         "DWHeatmap_P_Val",
                  #         height1 = '95vh',
                  #         width1 = 950,
                  #         height2 = '95vh',
                  #         width2 = 950,
                  #         inline = FALSE
                  #       )
                  #     ),
                  #     shiny::tabPanel(
                  #       "Table P_Val",
                  #       "Table of p-value; clustering order comes from clustering of p-values.",
                  #       DT::dataTableOutput("fullDWDTP_VALPval")
                  #     ),
                  #     shiny::tabPanel(
                  #       "Table Delta Methylation",
                  #       "Table of delta methylation; clustering order comes from clustering of p-values.",
                  #       DT::dataTableOutput("fullDWDTDMPval")
                  #     ),
                  #     shiny::tabPanel(
                  #       "Table Delta Methylation log(FC)",
                  #       "Table of log fold change(delta methylation); clustering order comes from clustering of p-values.",
                  #       DT::dataTableOutput("fullDWDTlogFCPval")
                  #     ),
                  #     shiny::tabPanel(
                  #       "Table N",
                  #       "Table of n; clustering order comes from clustering of p-values.",
                  #       DT::dataTableOutput("fullDWDTNPval")
                  #     ),
                  #     shiny::tabPanel(
                  #       "Dendrogram Probes",
                  #       plotly::plotlyOutput("fullDWPlotDendrogramProbesPval", height = "80%")
                  #     ),
                  #     shiny::tabPanel(
                  #       "Annotated Table Probes",
                  #       DT::dataTableOutput("fullDWDTProbesPval")
                  #     ),
                  #     shiny::tabPanel(
                  #       "Dendrogram Traits",
                  #       plotly::plotlyOutput("fullDWPlotDendrogramTraitsPval", height = "80%")
                  #     ),
                  #     shiny::tabPanel(
                  #       "Table Traits",
                  #       DT::dataTableOutput("fullDWDTTraitsPval")
                  #     ),
                  #     shiny::tabPanel(
                  #       "Histogram P_Val",
                  #       "Histogram of all p-values in condensed heatmap (number of p-values = number probes * number traits)",
                  #       plotly::plotlyOutput("fullDWHistP_Val", inline = TRUE)
                  #     )
                  #   )
                  # ) ##full DW data end
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
                  "VolcanoPlot Delta Methylation ~ log(FC)",
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
