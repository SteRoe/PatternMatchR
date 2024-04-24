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
      "Hostname / PID",
      shiny::verbatimTextOutput(outputId = "Sys.PID", placeholder = TRUE),
      shiny::tabsetPanel(
        shiny::tabPanel(
          "PatternMatchR",
          shinyBS::bsCollapse(
            id = "clpMain",
            open = c("Folders", "Merge", "Reduce Data", "Omit Traits", "Reduce Traits by Clustering", "Trait Reduced Data", "Condensed Distance Weighted Data (contains only CpG with nearby neighbours)", "Full Distance Weighted Data"),
            multiple = TRUE,
            # shinyBS::bsCollapsePanel(
            #   "session",
            #   shiny::fluidRow(
            #     shiny::tags$table(
            #       style = "width: 100%",
            #       shiny::tags$tr(
            #         shiny::tags$td(
            #           align = "center",
            #           shiny::actionButton(inputId = "btnExportCombinedData", label = "Export session data to <CombinedHM.RDS>"),
            #           shinyFiles::shinySaveButton(id = "save", label = "Export session data",
            #                                       title = "Save session data as file",
            #                                       filetype = list(RDS = "RDS", text = "txt"),
            #                                       viewtype = "detail")
            #         ),
            #         shiny::tags$td(
            #           align = "center",
            #           shiny::actionButton(inputId = "btnImportCombinedData", label = "Import session data from <CombinedHM.RDS>"),
            #           shinyFiles::shinyFilesButton(id = "file", label = "Import session data",
            #                                        title = "Load session data from file",
            #                                        multiple = FALSE, viewtype = "detail")
            #         )
            #       )
            #     )
            #   )
            # ),
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
            ),
            shinyjs::disabled(shinyBS::bsCollapsePanel(
              "Merge",
              shiny::fluidRow(
                shiny::actionButton("btnMerge", label = "Step 2: Merge data from all folders"),
                shiny::verbatimTextOutput("txtMergeOut", placeholder = TRUE)
              )
            )),
            shinyjs::disabled(shinyBS::bsCollapsePanel(
              "Count Borders",
              shiny::tabsetPanel(
                shiny::tabPanel("P_VAL border vs. # of probes",
                  shiny::fluidRow(
                    shinyjs::disabled(shiny::actionButton("btnCountP_ValProbes", label = "Count Probes for p-values (may take a long time)")),
                    shinyjs::disabled(shiny::actionButton("btnCountProbesP_ValParallel", label = "Count Probes for p-values parallel (may take less time)"))
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
            )),
            shinyjs::disabled(shinyBS::bsCollapsePanel(
              "Reduce Data",
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
              shinyjs::disabled(shiny::actionButton("btnReduce", label = "Step 3: Reduce data (omit CpGs) by p-value, DM or n limit")),
              shiny::verbatimTextOutput("txtPReduceOut", placeholder = TRUE)
            )),
            shinyjs::disabled(shinyBS::bsCollapsePanel(
              "Omit Traits",
              shiny::fluidRow(
                shiny::column(
                  width = 4,
                  shinyjs::disabled(shiny::sliderInput(
                    "sldNumClusters",
                    "number of clusters (omit traits)",
                    min = 0,
                    max = 0,
                    step = 0,
                    value = 0 # value = c(1, 10)
                  ))
                )
              ),
              shinyjs::disabled(shiny::actionButton("btnOmitTraits", label = "Step 4: Omit Traits")),
              shiny::verbatimTextOutput("txtOmitOut", placeholder = TRUE),
              shiny::fluidRow(
                shiny::tabsetPanel(
                  shiny::tabPanel(
                    "Dendrogram Traits",
                    shiny::verbatimTextOutput("txtDendrogramTraitsLong", placeholder = TRUE),
                    shiny::plotOutput("plotDendrogramTraitsLong")
                  ),
                  shiny::tabPanel(
                    "Clustergram Traits",
                    shiny::verbatimTextOutput("txtClustergramTraitsLong", placeholder = TRUE),
                    shiny::plotOutput("plotClustergramTraitsLong")
                  ),
                  shiny::tabPanel(
                    "DT Cluster Medoids Traits",
                    shiny::verbatimTextOutput("txtDTTraitsMedoids", placeholder = TRUE),
                    DT::dataTableOutput("DTTraitsMedoids")
                  ),
                  shiny::tabPanel(
                    "DT Cluster Assignment Traits",
                    shiny::verbatimTextOutput("txtDTTraitsClusters", placeholder = TRUE),
                    DT::dataTableOutput("DTTraitsClusters")
                  ),
                  shiny::tabPanel(
                    "Histograms/DT on CpG Distances of Clustering Results",
                    shiny::tabsetPanel(
                      shiny::tabPanel(
                        "Mean Distance Probes = 10 CpG up/down",
                        shiny::tabsetPanel(
                          shiny::tabPanel(
                            "Histogram",
                            shiny::tabsetPanel(
                              shiny::tabPanel(
                                "min",
                                plotly::plotlyOutput("histMinDistance10", inline = TRUE)
                              ),
                              shiny::tabPanel(
                                "mean",
                                plotly::plotlyOutput("histMeanDistance10", inline = TRUE)
                              ),
                              shiny::tabPanel(
                                "max",
                                plotly::plotlyOutput("histMaxDistance10", inline = TRUE)
                              )
                            )
                          ),
                          shiny::tabPanel(
                            "Table",
                            shiny::tabsetPanel(
                              shiny::tabPanel(
                                "reduced",
                                DT::dataTableOutput("DTDistance10reduced")
                              ),
                              shiny::tabPanel(
                                "full",
                                #table with histogram values
                                DT::dataTableOutput("DTDistance10")
                              )
                            )
                          )
                        )
                      ),
                      shiny::tabPanel(
                        "Mean Distance Probes = 100 CpG up/down",
                        shiny::tabsetPanel(
                          shiny::tabPanel(
                            "Histogram",
                            shiny::tabsetPanel(
                              shiny::tabPanel(
                                "min",
                                plotly::plotlyOutput("histMinDistance100", inline = TRUE)
                              ),
                              shiny::tabPanel(
                                "mean",
                                plotly::plotlyOutput("histMeanDistance100", inline = TRUE)
                              ),
                              shiny::tabPanel(
                                "max",
                                plotly::plotlyOutput("histMaxDistance100", inline = TRUE)
                              )
                            )
                          ),
                          shiny::tabPanel(
                            "Table",
                            shiny::tabsetPanel(
                              shiny::tabPanel(
                                "reduced",
                                DT::dataTableOutput("DTDistance100reduced")
                              ),
                              shiny::tabPanel(
                                "full",
                                #table with histogram values
                                DT::dataTableOutput("DTDistance100")
                              )
                            )
                          )
                        )
                      ),
                      shiny::tabPanel(
                        "Mean Distance Probes = 1000 CpG up/down",
                        shiny::tabsetPanel(
                          shiny::tabPanel(
                            "Histogram",
                            shiny::tabsetPanel(
                              shiny::tabPanel(
                                "min",
                                plotly::plotlyOutput("histMinDistance1000", inline = TRUE)
                              ),
                              shiny::tabPanel(
                                "mean",
                                plotly::plotlyOutput("histMeanDistance1000", inline = TRUE)
                              ),
                              shiny::tabPanel(
                                "max",
                                plotly::plotlyOutput("histMaxDistance1000", inline = TRUE)
                              )
                            )
                          ),
                          shiny::tabPanel(
                            "Table",
                            shiny::tabsetPanel(
                              shiny::tabPanel(
                                "reduced",
                                DT::dataTableOutput("DTDistance1000reduced")
                              ),
                              shiny::tabPanel(
                                "full",
                                #table with histogram values
                                DT::dataTableOutput("DTDistance1000")
                              )
                            )
                          )
                        )
                      ),
                      shiny::tabPanel(
                        "Mean Distance Probes = 10000 CpG up/down",
                        shiny::tabsetPanel(
                          shiny::tabPanel(
                            "Histogram",
                            shiny::tabsetPanel(
                              shiny::tabPanel(
                                "min",
                                plotly::plotlyOutput("histMinDistance10000", inline = TRUE)
                              ),
                              shiny::tabPanel(
                                "mean",
                                plotly::plotlyOutput("histMeanDistance10000", inline = TRUE)
                              ),
                              shiny::tabPanel(
                                "max",
                                plotly::plotlyOutput("histMaxDistance10000", inline = TRUE)
                              )
                            )
                          ),
                          shiny::tabPanel(
                            "Table",
                            shiny::tabsetPanel(
                              shiny::tabPanel(
                                "reduced",
                                DT::dataTableOutput("DTDistance10000reduced")
                              ),
                              shiny::tabPanel(
                                "full",
                                #table with histogram values
                                DT::dataTableOutput("DTDistance10000")
                              )
                            )
                          )
                        )
                      )
                    )
                  )
                )
              )
            )),
            shinyjs::disabled(shinyBS::bsCollapsePanel(
              "Trait Reduced Data",
              shiny::fluidRow(
                shiny::tabsetPanel(
                  shiny::tabPanel(
                    "HeatMap P_Val",
                    shiny::tabsetPanel(
                      shiny::tabPanel(
                        "HeatMap P_Val",
                        shinyjs::disabled(shiny::actionButton("btnPlotCombinedHM_P_Val", label = "Step 5: Plot Heatmap P_Val")),
                        shiny::verbatimTextOutput("txtHMDescription_P_Val", placeholder = TRUE),
                        shiny::column(
                          width = 6,
                          shinyjs::disabled(shiny::numericInput(inputId = "numHMHSize", label = "Width", value = 4000, min = 1000, max = 10000))
                        ),
                        shiny::column(
                          width = 6,
                          shinyjs::disabled(shiny::numericInput(inputId = "numHMVSize", label = "Height", value = 4000, min = 1000, max = 10000))
                        ),
                        InteractiveComplexHeatmap::InteractiveComplexHeatmapOutput(
                          "Heatmap_P_Val",
                          height1 = 4000,
                          width1 = 4000,
                          height2 = 4000,
                          width2 = 4000,
                          inline = FALSE
                        )
                      ),
                      shiny::tabPanel(
                        "Search for CpG in Heatmap",
                        shiny::textAreaInput(inputId = "txtSearchCpG", label = "Search CpG", value = ""),
                        shiny::verbatimTextOutput(outputId = "txtSearchResultCpG"),
                        shiny::actionButton("btnSearchCpGHM", label = "Search CpG"),
                        "Search for Trait in Heatmap",
                        shiny::textAreaInput(inputId = "txtSearchTrait", label = "Search Trait", value = ""),
                        shiny::verbatimTextOutput(outputId = "txtSearchResultTrait"),
                        shiny::actionButton("btnSearchTraitHM", label = "Search Trait")
                      ),
                      shiny::tabPanel(
                        "SPLOM of Original Data",
                        shiny::selectizeInput("markingVar",
                          label = "select variable for color marking",
                          choices = NULL,
                          width = "100%"
                        ),
                        shiny::column(
                          width = 6,
                          shiny::numericInput(inputId = "numSPLOMHSize", label = "Width (checked max 3000)", value = 2500, min = 1000, max = 3000)
                        ),
                        shiny::column(
                          width = 6,
                          shiny::numericInput(inputId = "numSPLOMVSize", label = "Height (checked max 12000)", value = 8000, min = 1000, max = 12000)
                        ),
                        shiny::tabsetPanel(
                          shiny::tabPanel(
                            "SPLOM from selected area in heatmap",
                            plotly::plotlyOutput("SPLOM",
                                             height = 2000, #height = "100%",
                                             width = 1500, #width = "100%",
                                             inline = TRUE)
                          ),
                          shiny::tabPanel(
                            "SPLOM trait/trait",
                            plotly::plotlyOutput("SPLOMTrait",
                                                 height = 2000, #height = "100%",
                                                 width = 1500, #width = "100%",
                                                 inline = TRUE)
                          ),
                          shiny::tabPanel(
                            "SPLOM probe/probe",
                            plotly::plotlyOutput("SPLOMProbe",
                                                 height = 2000, #height = "100%",
                                                 width = 1500, #width = "100%",
                                                 inline = TRUE)
                          )
                        )
                      ),
                      shiny::tabPanel(
                        "Selected CpG",
                        DT::dataTableOutput("DTSelectedCpG")
                      ),
                      shiny::tabPanel(
                        "Selected trait",
                        DT::dataTableOutput("DTSelectedTrait")
                      ),
                      shiny::tabPanel(
                        "Selected p-value",
                        DT::dataTableOutput("DTSelectedP_Val")
                      ),
                      shiny::tabPanel(
                        "Selected Histogram on Delta Methylation",
                        shiny::verbatimTextOutput("txtselectDMTraits", placeholder = TRUE),
                        plotly::plotlyOutput("selectHistDM", inline = TRUE)
                      )
                    )
                  ),
                  shiny::tabPanel(
                    "Table P_VAL",
                    # DT tableout
                    DT::dataTableOutput("traitReducedDTP_VAL")
                  ),
                  shiny::tabPanel(
                    "Table Delta Methylation",
                    # DT tableout
                    DT::dataTableOutput("traitReducedDTDM")
                  ),
                  shiny::tabPanel(
                    "Table N",
                    # DT tableout
                    DT::dataTableOutput("traitReducedDTN")
                  ),
                  shiny::tabPanel(
                    "Dendrogram Probes",
                    plotly::plotlyOutput("traitReducedPlotDendrogramProbes", height = "80%")
                  ),
                  shiny::tabPanel(
                    "Annotated Table Probes",
                    DT::dataTableOutput("traitReducedDTProbes")
                  ),
                  shiny::tabPanel(
                    "Dendrogram Traits",
                    plotly::plotlyOutput("traitReducedPlotDendrogramTraits", height = "80%")
                  ),
                  shiny::tabPanel(
                    "Table Traits",
                    DT::dataTableOutput("traitReducedDTTraits")
                  ),
                  shiny::tabPanel(
                    "Histogram P_Val",
                    "Histogram of all p-values in heatmap",
                    plotly::plotlyOutput("traitReducedHistP_Val", inline = TRUE)
                  )
                )
              )
            )),
            shinyjs::disabled(shinyBS::bsCollapsePanel(
              "Condensed Distance Weighted Data (contains only CpG with nearby neighbours)",
              shiny::fluidRow(
                shiny::tabsetPanel(
                  shiny::tabPanel(
                    "HeatMap P_Val",
                    shiny::actionButton("btnPlotCombinedCondDWHM_P_Val", label = "Step 6a: Plot condensed distance weighted Heatmap P_Val"),
                    shiny::verbatimTextOutput("txtCondDWHMDescription_P_Val", placeholder = TRUE),
                    shiny::column(
                      width = 6,
                      shinyjs::disabled(shiny::numericInput(inputId = "numCondDWHMHSize", label = "Width", value = 4000, min = 1000, max = 10000))
                    ),
                    shiny::column(
                      width = 6,
                      shinyjs::disabled(shiny::numericInput(inputId = "numCondDWHMVSize", label = "Height", value = 4000, min = 1000, max = 10000))
                    ),
                    InteractiveComplexHeatmap::InteractiveComplexHeatmapOutput(
                      "condDWHeatmap_P_Val",
                      height1 = 4000,
                      width1 = 4000,
                      height2 = 4000,
                      width2 = 4000,
                      inline = FALSE
                    )
                  ),
                  shiny::tabPanel(
                    "Table P_Val",
                    DT::dataTableOutput("condDWDTP_VAL")
                  ),
                  shiny::tabPanel(
                    "Table Delta Methylation",
                    DT::dataTableOutput("condDWDTDM")
                  ),
                  shiny::tabPanel(
                    "Table N",
                    DT::dataTableOutput("condDWDTN")
                  ),
                  shiny::tabPanel(
                    "Dendrogram Probes",
                    plotly::plotlyOutput("condDWPlotDendrogramProbes", height = "80%")
                  ),
                  shiny::tabPanel(
                    "Annotated Table Probes",
                    DT::dataTableOutput("condDWDTProbes")
                  ),
                  shiny::tabPanel(
                    "Dendrogram Traits",
                    plotly::plotlyOutput("condDWPlotDendrogramTraits", height = "80%")
                  ),
                  shiny::tabPanel(
                    "Table Traits",
                    DT::dataTableOutput("condDWDTTraits")
                  ),
                  shiny::tabPanel(
                    "Histogram P_Val",
                    "Histogram of all p-values in heatmap",
                    plotly::plotlyOutput("condDWHistP_Val", inline = TRUE)
                  )
                )
              )
            )),
            shinyjs::disabled(shinyBS::bsCollapsePanel(
              "Full Distance Weighted Data",
              shiny::fluidRow(
                shiny::tabsetPanel(
                  shiny::tabPanel(
                    "HeatMap P_Val",
                    shiny::actionButton("btnPlotCombinedDWHM_P_Val", label = "Step 6b: Plot distance weighted Heatmap P_Val"),
                    shiny::verbatimTextOutput("txtDWHMDescription_P_Val", placeholder = TRUE),
                    shiny::column(
                      width = 6,
                      shinyjs::disabled(shiny::numericInput(inputId = "numDWHMHSize", label = "Width", value = 4000, min = 1000, max = 10000))
                    ),
                    shiny::column(
                      width = 6,
                      shinyjs::disabled(shiny::numericInput(inputId = "numDWHMVSize", label = "Height", value = 4000, min = 1000, max = 10000))
                    ),
                    InteractiveComplexHeatmap::InteractiveComplexHeatmapOutput(
                      "DWHeatmap_P_Val",
                      height1 = 4000,
                      width1 = 4000,
                      height2 = 4000,
                      width2 = 4000,
                      inline = FALSE
                    )
                  ),
                  shiny::tabPanel(
                    "Table P_Val",
                    DT::dataTableOutput("fullDWDTP_VAL")
                  ),
                  shiny::tabPanel(
                    "Table Delta Methylation",
                    DT::dataTableOutput("fullDWDTDM")
                  ),
                  shiny::tabPanel(
                    "Table N",
                    DT::dataTableOutput("fullDWDTN")
                  ),
                  shiny::tabPanel(
                    "Dendrogram Probes",
                    plotly::plotlyOutput("fullDWPlotDendrogramProbes", height = "80%")
                  ),
                  shiny::tabPanel(
                    "Annotated Table Probes",
                    DT::dataTableOutput("fullDWDTProbes")
                  ),
                  shiny::tabPanel(
                    "Dendrogram Traits",
                    plotly::plotlyOutput("fullDWPlotDendrogramTraits", height = "80%")
                  ),
                  shiny::tabPanel(
                    "Table Traits",
                    DT::dataTableOutput("fullDWDTTraits")
                  ),
                  shiny::tabPanel(
                    "Histogram P_Val",
                    plotly::plotlyOutput("fullDWHistP_Val", inline = TRUE)
                  )
                )
              )
            )
          ))
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
                shiny::textInput(inputId = "inpDNAmFile", label = "File with DNAm", value = "FileName", width = NULL, placeholder = TRUE)
              )
            ),
            shinyBS::bsCollapsePanel(
              "minimum cases for trait",
              shiny::fluidRow(
                shiny::textInput(inputId = "txtCases", label = "minimum # cases for each trait", value = 100, width = NULL, placeholder = TRUE)
              )
            ),
            shinyBS::bsCollapsePanel(
              "p-val settings",
              shiny::fluidRow(
                shiny::column(
                  width = 2,
                  shiny::textInput(inputId = "txtMaxProbes", label = "maximum Probes",
                                   value = 500000, width = NULL, placeholder = TRUE)
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
                )
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
              )
            )
          )
        )
      )
    )
  )
}
