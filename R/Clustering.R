Clustering_UI <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    shiny::fluidRow(
      shiny::column(
        width = 4,
        shinyjs::disabled(shiny::sliderInput(
          ns("sldNumClustersTraits"),
          "number of trait clusters",
          min = 0,
          max = 0,
          step = 0,
          value = 0 # value = c(1, 10)
        ))
      )
    ) #end fluidRow
  )
}

ClusteringTraits_UI <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    shinyjs::disabled(
      if (id == "PVal") {
        shiny::actionButton(ns("btnOmitTraits"), label = "Step 4a: Omit traits (p-val)")
      }
      else {
        shiny::actionButton(ns("btnOmitTraits"), label = "Step 4b: Omit traits (log(FC))")
      }
    ),
    shiny::verbatimTextOutput(ns("txtOmitTraitsOut"), placeholder = TRUE),
    shiny::fluidRow(
      shiny::tabsetPanel(#tabset
        shiny::tabPanel( #tab non-distance weighted data
          "Non-Distance Weighted Trait Clustering Results",
          shiny::tabsetPanel(
            shiny::tabPanel(
              "Dendrogram Traits",
              shiny::plotOutput(outputId = ns("plotDendrogramTraitsLong"))
            #  plotly::plotlyOutput(ns("plotDendrogramTraitsLong")) #, inline = TRUE
            ),
            shiny::tabPanel(
              "Clustergram Traits",
              shiny::plotOutput(outputId = ns("plotClustergramTraitsLong"))
            #  plotly::plotlyOutput(outputId = ns("plotClustergramTraitsLong"), width = "100%")
            ),
            shiny::tabPanel(
              "DT Cluster Medoids Traits",
              #DT::dataTableOutput(ns("DTTraitsMedoids")),
              DT::DTOutput(outputId = ns("DTTraitsMedoids"))
            ),
            shiny::tabPanel(
              "DT Cluster Assignment Traits",
              DT::DTOutput(outputId = ns("DTTraitsClusters"))
            )
          ) #end tabsetPanel
        ), #end tab non-distance weighted data
        # shiny::tabPanel( #tab distance weighted data
        #   "Distance Weighted Clustering Results - experimental",
        #   shiny::tabsetPanel(
        #     shiny::tabPanel(
        #       "Dendrogram Traits",
        #       shiny::plotOutput(outputId = ns("plotDWDendrogramTraitsLong"))
        #       #  plotly::plotlyOutput(ns("plotDWDendrogramTraitsLong")) #, inline = TRUE
        #     ),
        #     shiny::tabPanel(
        #       "Clustergram Traits",
        #       shiny::plotOutput(outputId = ns("plotDWClustergramTraitsLong"))
        #       #  plotly::plotlyOutput(outputId = ns("plotDWClustergramTraitsLong"), width = "100%")
        #     ),
        #     shiny::tabPanel(
        #       "DT Cluster Medoids Traits",
        #       #DT::dataTableOutput(ns("DTDWTraitsMedoids")),
        #       DT::DTOutput(outputId = ns("DTDWTraitsMedoids"))
        #     ),
        #     shiny::tabPanel(
        #       "DT Cluster Assignment Traits",
        #       DT::DTOutput(outputId = ns("DTDWTraitsClusters"))
        #     )
        #   )
        # )#end tab distance weighted data
      )#end tabset
    ) #end fluidRow
  ) #end tagList
}

ClusteringProbesGeneral_UI <- function(id) {
  ns <- shiny::NS(id)
  #shinyjs::disabled(
  shiny::tagList(
    # shiny::sliderInput(
    #   "sldNumTest",
    #   "test",
    #   min = 0,
    #   max = 10,
    #   step = 1,
    #   value = 8),
    shiny::fluidRow(
      shiny::column(
        width = 4,
        shiny::sliderInput(
          ns("sldNumClustersProbes"),
          "number of probe clusters",
          min = 0,
          max = 0,
          step = 0,
          value = 0 # value = c(1, 10)
        )
      ),
    # ), #end fluidRow
    # shiny::fluidRow(
      shiny::column(
        width = 4,
        shiny::sliderInput(
          ns("sldNumNeighbours"),
          "distance to look for neighbours",
          min = 10,
          max = 10000,
          step = 10,
          value = 5000
        )
      )
    ) #end fluidRow
  ) #end tagList
  #) #end disabled
}

ClusteringProbes_UI <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    shinyjs::disabled(
      if (id == "PVal") {
        shiny::actionButton(ns("btnOmitProbes"), label = "Step 5a: Omit probes (p-val)")
      }
      else {
        shiny::actionButton(ns("btnOmitProbes"), label = "Step 5b: Omit probes (log(FC))")
      }
    ),
    shiny::verbatimTextOutput(ns("txtOmitProbesOut"), placeholder = TRUE),
    shiny::fluidRow(
      shiny::tabsetPanel(#tabset
        shiny::tabPanel( #tab non-distance weighted data
          "Non-Distance Weighted Probes Clustering Results",
          shiny::tabsetPanel(
            shiny::tabPanel(
              "Dendrogram Probes",
              shiny::plotOutput(outputId = ns("plotDendrogramProbesLong"))
              #  plotly::plotlyOutput(ns("plotDendrogramTraitsLong")) #, inline = TRUE
            ),
            shiny::tabPanel(
              "Clustergram Probes",
              shiny::plotOutput(outputId = ns("plotClustergramProbesLong"))
              #  plotly::plotlyOutput(outputId = ns("plotClustergramTraitsLong"), width = "100%")
            ),
            shiny::tabPanel(
              "DT Cluster Medoids Probes",
              #DT::dataTableOutput(ns("DTTraitsMedoids")),
              DT::DTOutput(outputId = ns("DTProbesMedoids"))
            ),
            shiny::tabPanel(
              "DT Cluster Assignment Probes",
              DT::DTOutput(outputId = ns("DTProbesClusters"))
            ),
            shiny::tabPanel(
              "DT Annotated Probes",
              DT::DTOutput(outputId = ns("DTAnnotatedProbes"))
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
                          plotly::plotlyOutput(outputId = ns("histMinDistance10"), inline = TRUE)
                        ),
                        shiny::tabPanel(
                          "mean",
                          plotly::plotlyOutput(outputId = ns("histMeanDistance10"), inline = TRUE)
                        ),
                        shiny::tabPanel(
                          "max",
                          plotly::plotlyOutput(outputId = ns("histMaxDistance10"), inline = TRUE)
                        )
                      )
                    ),
                    shiny::tabPanel(
                      "Table",
                      shiny::tabsetPanel(
                        shiny::tabPanel(
                          "reduced",
                          DT::dataTableOutput(outputId = ns("DTDistance10reduced"))
                        ),
                        shiny::tabPanel(
                          "full",
                          #table with histogram values
                          DT::dataTableOutput(outputId = ns("DTDistance10"))
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
                          plotly::plotlyOutput(outputId = ns("histMinDistance100"), inline = TRUE)
                        ),
                        shiny::tabPanel(
                          "mean",
                          plotly::plotlyOutput(outputId = ns("histMeanDistance100"), inline = TRUE)
                        ),
                        shiny::tabPanel(
                          "max",
                          plotly::plotlyOutput(outputId = ns("histMaxDistance100"), inline = TRUE)
                        )
                      )
                    ),
                    shiny::tabPanel(
                      "Table",
                      shiny::tabsetPanel(
                        shiny::tabPanel(
                          "reduced",
                          DT::dataTableOutput(outputId = ns("DTDistance100reduced"))
                        ),
                        shiny::tabPanel(
                          "full",
                          #table with histogram values
                          DT::dataTableOutput(outputId = ns("DTDistance100"))
                        ) #end tabPanel
                      ) #end tabsetPanel
                    ) #end tabPanel
                  ) #end tabsetPanel
                ), #end tabPanel
                shiny::tabPanel(
                  "Mean Distance Probes = 1000 CpG up/down",
                  shiny::tabsetPanel(
                    shiny::tabPanel(
                      "Histogram",
                      shiny::tabsetPanel(
                        shiny::tabPanel(
                          "min",
                          plotly::plotlyOutput(outputId = ns("histMinDistance1000"), inline = TRUE)
                        ),
                        shiny::tabPanel(
                          "mean",
                          plotly::plotlyOutput(outputId = ns("histMeanDistance1000"), inline = TRUE)
                        ),
                        shiny::tabPanel(
                          "max",
                          plotly::plotlyOutput(outputId = ns("histMaxDistance1000"), inline = TRUE)
                        )
                      )
                    ),
                    shiny::tabPanel(
                      "Table",
                      shiny::tabsetPanel(
                        shiny::tabPanel(
                          "reduced",
                          DT::dataTableOutput(outputId = ns("DTDistance1000reduced"))
                        ),
                        shiny::tabPanel(
                          "full",
                          #table with histogram values
                          DT::dataTableOutput(outputId = ns("DTDistance1000"))
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
                          plotly::plotlyOutput(outputId = ns("histMinDistance10000"), inline = TRUE)
                        ),
                        shiny::tabPanel(
                          "mean",
                          plotly::plotlyOutput(outputId = ns("histMeanDistance10000"), inline = TRUE)
                        ),
                        shiny::tabPanel(
                          "max",
                          plotly::plotlyOutput(outputId = ns("histMaxDistance10000"), inline = TRUE)
                        )
                      )
                    ),
                    shiny::tabPanel(
                      "Table",
                      shiny::tabsetPanel(
                        shiny::tabPanel(
                          "reduced",
                          DT::dataTableOutput(outputId = ns("DTDistance10000reduced"))
                        ),
                        shiny::tabPanel(
                          "full",
                          #table with histogram values
                          DT::dataTableOutput(outputId = ns("DTDistance10000"))
                        ) #end tabPanel
                      ) #end tabsetPanel
                    ) #end tabPanel
                  ) #end tabsetPanel
                ) #end tabPanel
              ) #end tabsetPanel
            ) #end tabPanel
          ) #end tabsetPanel
        ), #end tab non-distance weighted data
        # shiny::tabPanel( #tab distance weighted data
        #   "Distance Weighted Clustering Results - experimental",
        #   shiny::tabsetPanel(
        #     shiny::tabPanel(
        #       "Dendrogram Traits",
        #       shiny::plotOutput(outputId = ns("plotDWDendrogramTraitsLong"))
        #       #  plotly::plotlyOutput(ns("plotDWDendrogramTraitsLong")) #, inline = TRUE
        #     ),
        #     shiny::tabPanel(
        #       "Clustergram Traits",
        #       shiny::plotOutput(outputId = ns("plotDWClustergramTraitsLong"))
        #       #  plotly::plotlyOutput(outputId = ns("plotDWClustergramTraitsLong"), width = "100%")
        #     ),
        #     shiny::tabPanel(
        #       "DT Cluster Medoids Traits",
        #       #DT::dataTableOutput(ns("DTDWTraitsMedoids")),
        #       DT::DTOutput(outputId = ns("DTDWTraitsMedoids"))
        #     ),
        #     shiny::tabPanel(
        #       "DT Cluster Assignment Traits",
        #       DT::DTOutput(outputId = ns("DTDWTraitsClusters"))
        #     )
        #   )
        # )#end tab distance weighted data
      )#end tabset
    ) #end fluidRow
  ) #end tagList
}

Clustering_SERVER <- function(id, session) {
  shiny::moduleServer(id, function(input, output, session) {
    base::tryCatch({
      shiny::observe({
        #modify sldNumClustersTraits...
        if (!is.valid(session$userData$sessionVariables$pReducedDataStructure()$combinedDFP_Val_Labels$dfP_Val)) {
          shinyjs::disable("sldNumClustersTraits")
        }
        else {
          result <- session$userData$sessionVariables$pReducedDataStructure()$combinedDFP_Val_Labels
          maxTraits <- ncol(result$dfP_Val)
          #value <- input$sldNumClustersTraits #maxTraits
          #if (value == 0) {value <- maxTraits}
          value <- maxTraits
          shiny::updateSliderInput(session = session, inputId = "sldNumClustersTraits", max = maxTraits, min = 1, value = value, step = 1)
          shinyjs::enable("sldNumClustersTraits")
        }
      })

      session$userData$sessionVariables$numClustersTraits <- shiny::reactive({input$sldNumClustersTraits})
    },
    error = function(e) {
      base::message("An error occurred in moduleServer in ClusteringTraits_SERVER:\n", e)
      browser() #should not happen
    },
    warning = function(w) {
      base::message("An error occurred in moduleServer in ClusteringTraits_SERVER:\n", w)
      browser() #should not happen
    },
    finally = {
      base::print(base::paste0(sysTimePID(), " finished moduleServer in ClusteringTraits_SERVER."))
    })
  }) #end moduleServer
}

ClusteringTraits_SERVER <- function(id, pReducedDataStructure, traitReducedDataStructure, session) {
  shiny::moduleServer(id, function(input, output, session) {
    base::tryCatch({
      shiny::observe({
        if (!is.valid(session$userData$sessionVariables$pReducedDataStructure()$combinedDFP_Val_Labels$dfP_Val)) {
          shinyjs::disable("sldNumClustersTraits")
        }
        else {
          shinyjs::enable("sldNumClustersTraits")
        }
      })

      shiny::observe({
        if (id == "PVal") {
          if (!is.valid(pReducedDataStructure()$combinedDFP_Val_Labels$dfP_Val)) {
            shinyjs::disable("btnOmitTraits")
          }
          else {
            shinyjs::enable("btnOmitTraits")
          }
        }
      })

      shiny::observe({
        if (id == "LogFC") {
          if (!is.valid(pReducedDataStructure()$combinedDFP_Val_Labels$dfLogFC)) {
            shinyjs::disable("btnOmitTraits")
          }
          else {
            shinyjs::enable("btnOmitTraits")
          }
        }
      })

      shiny::observeEvent(input$btnOmitTraits,
        ignoreInit = TRUE,
        {
          base::tryCatch(
            {
              if (id == "PVal") {
                result <- getTraitReducedDataStructure(id, input, output, session)
                #session$userData$sessionVariables$traitReducedDataStructurePVal(result)
                session$userData$sessionVariables$traitReducedDataPVal(result)
              }
              else if (id == "LogFC") {
                result <- getTraitReducedDataStructure(id, input, output, session)
                #result <- getTraitReducedDataStructureLogFC(input, output, session)
                session$userData$sessionVariables$traitReducedDataLogFC(result)
              }
              else {
                browser() #should not happen
              }
            },
            error = function(e) {
              base::message("An error occurred in shiny::observeEvent(input$btnOmitTraits):\n", e)
            },
            warning = function(w) {
              base::message("A warning occurred in shiny::observeEvent(input$btnOmitTraits):\n", w)
            },
            finally = {
              base::print(base::paste0(sysTimePID(), " finished omitting traits."))
            }
          )
        },
        ignoreNULL = FALSE
      ) #end observeEvent

      output$txtOmitTraitsOut <- shiny::reactive({
        result <- NULL
        base::tryCatch(
          {
            base::print(base::paste0(sysTimePID(), " start generating output$txtOmitTraitsOut."))
            if (is.valid(traitReducedDataStructure()$combinedDFP_Val_Labels$dfP_Val)) {
              result <- updateTxtOmitTraitsOut(id, traitReducedDataStructure())
              if (id == "PVal") {
                result <- base::paste0("Clustering by p-val: ", result)
              }
              else {
                result <- base::paste0("Clustering by log(FC): ", result)
              }
            }
            else {
              result <- NULL
            }
          },
          error = function(e) {
            base::message("An error occurred in shiny::reactive(output$txtOmitTraitsOut):\n", e)
          },
          warning = function(w) {
            base::message("A warning occurred in shiny::reactive(output$txtOmitTraitsOut):\n", w)
          },
          finally = {
            base::print(base::paste0(sysTimePID(), " finished generating output$txtOmitTraitsOut"))
            return(result)
          }
        )
      }) #end output$txtOmitTraitsOut

      shiny::observe({
          output$plotDendrogramTraitsLong <- shiny::renderPlot(getPlot(traitReducedDataStructure()$traitDendrogram))
          #output$plotDendrogramTraitsLong <- plotly::renderPlotly(getPlot(traitReducedDataStructure()$traitDendrogram))
      })

      shiny::observe({
          output$plotClustergramTraitsLong <- shiny::renderPlot(getPlot(traitReducedDataStructure()$traitClustergram))
          #output$plotClustergramTraitsLong <- plotly::renderPlotly(getPlot(traitReducedDataStructure()$traitClustergram))
      })

      shiny::observe({
          output$DTTraitsMedoids <- DT::renderDataTable(as.data.frame(getMedoidsTable(id, traitReducedDataStructure()$traitClusterMedoids)))
      })

      shiny::observe({
        output$DTTraitsClusters <- DT::renderDataTable(as.data.frame(getClustersTable(id, traitReducedDataStructure()$traitClusters,
                                                                                      traitReducedDataStructure()$traitClusterMedoids)))
      })
    },
    error = function(e) {
      base::message("An error occurred in moduleServer in ClusteringTraits_SERVER:\n", e)
      browser() #should not happen
    },
    warning = function(w) {
      base::message("An error occurred in moduleServer in ClusteringTraits_SERVER:\n", w)
      browser() #should not happen
    },
    finally = {
      base::print(base::paste0(sysTimePID(), " finished moduleServer in ClusteringTraits_SERVER."))
    })
  }) #end moduleServer
} #end Clustering_TraitsSERVER

ClusteringProbesGeneral_SERVER <- function(id, session) {
  shiny::moduleServer(id, function(input, output, session) {
    base::tryCatch({
      shiny::observe({
        if (!is.valid(session$userData$sessionVariables$traitReducedDataStructurePVal()$combinedDFP_Val_Labels$dfP_Val)
            && !is.valid(session$userData$sessionVariables$traitReducedDataStructureLogFC()$combinedDFP_Val_Labels$dfLogFC)) {
          shinyjs::disable("sldNumClustersProbes")
          shinyjs::disable("sldNumNeighbours")
        }
        else {
          resultPVal <- session$userData$sessionVariables$traitReducedDataStructurePVal()$combinedDFP_Val_Labels
          resultLogFC <- session$userData$sessionVariables$traitReducedDataStructureLogFC()$combinedDFP_Val_Labels
#browser()
          maxProbesPVal <- nrow(resultPVal$dfP_Val)
          maxProbesLogFC <- nrow(resultLogFC$dfLogFC)
          maxProbes <- base::max(maxProbesPVal, maxProbesLogFC)
          value <- maxProbes
          shiny::updateSliderInput(session = session, inputId = "sldNumClustersProbes", max = maxProbes, min = 1, value = value, step = 1)
          shinyjs::enable("sldNumClustersProbes")
          shinyjs::enable("sldNumNeighbours")
        }
      })
      session$userData$sessionVariables$numClustersProbes <- shiny::reactive({input$sldNumClustersProbes})
      session$userData$sessionVariables$numNeighbours <- shiny::reactive({input$sldNumNeighbours})
    },
    error = function(e) {
      base::message("An error occurred in moduleServer in ClusteringProbesGeneral_SERVER:\n", e)
      browser() #should not happen
    },
    warning = function(w) {
      base::message("An error occurred in moduleServer in ClusteringProbesGeneral_SERVER:\n", w)
      browser() #should not happen
    },
    finally = {
      base::print(base::paste0(sysTimePID(), " finished moduleServer in ClusteringProbesGeneral_SERVER"))
    })
  }) #end moduleServer
} # end ClusteringProbesGeneral_SERVER

ClusteringProbes_SERVER <- function(id, traitReducedDataStructure, probeReducedDataStructure, session) {
  shiny::moduleServer(id, function(input, output, session) {
    base::tryCatch({
      shiny::observe({
        if (id == "PVal") {
          if (!is.valid(traitReducedDataStructure()$combinedDFP_Val_Labels$dfP_Val)) {
            shinyjs::disable("btnOmitProbes")
          }
          else {
            shinyjs::enable("btnOmitProbes")
          }
        }
      })

      shiny::observe({
        if (id == "LogFC") {
          if (!is.valid(traitReducedDataStructure()$combinedDFP_Val_Labels$dfLogFC)) { #session$userData$sessionVariables$traitReducedDataLogFC
            shinyjs::disable("btnOmitProbes")
          }
          else {
            shinyjs::enable("btnOmitProbes")
          }
        }
      })

      shiny::observeEvent(input$btnOmitProbes,
        ignoreInit = TRUE,
        {
          base::tryCatch(
            {
              result <- NULL
              if (id == "PVal") {
                result <- getProbeReducedDataStructure(id, input, output, session)
                session$userData$sessionVariables$probeReducedDataPVal(result)
              }
              else if (id == "LogFC") {
                result <- getProbeReducedDataStructure(id, input, output, session)
                session$userData$sessionVariables$probeReducedDataLogFC(result)
              }
              else {
                browser() #should not happen
              }
            },
            error = function(e) {
              base::message("An error occurred in shiny::observeEvent(input$btnOmitProbes):\n", e)
            },
            warning = function(w) {
              base::message("A warning occurred in shiny::observeEvent(input$btnOmitProbes):\n", w)
            },
            finally = {
              base::print(base::paste0(sysTimePID(), " finished omitting probes."))
            }
          )
        },
        ignoreNULL = FALSE
      ) #end observeEvent

      output$txtOmitProbesOut <- shiny::reactive({
        base::tryCatch(
          {
            base::print(base::paste0(sysTimePID(), " start generating output$txtOmitProbesOut."))
            #if (is.valid(probeReducedDataStructure())) {
            if (is.valid(probeReducedDataStructure()$combinedDFP_Val_Labels$dfP_Val)) {
              result <- updateTxtOmitProbesOut(id, probeReducedDataStructure())
              if (id == "PVal") {
                result <- base::paste0("Clustering by p-val: ", result)
              }
              else {
                result <- base::paste0("Clustering by log(FC): ", result)
              }
            }
            else {
              result <- NULL
            }
          },
          error = function(e) {
            base::message("An error occurred in shiny::reactive(output$txtOmitProbesOut):\n", e)
          },
          warning = function(w) {
            base::message("A warning occurred in shiny::reactive(output$txtOmitProbesOut):\n", w)
          },
          finally = {
            base::print(base::paste0(sysTimePID(), " finished generating output$txtOmitProbesOut"))
            return(result)
          }
        )
      }) #end output$txtOmitProbesOut

      distNeigboursProbes10 <- shiny::reactive(calculateDistanceNeigboursProbes(wd = session$userData$packageWd, clustResProbes = probeReducedDataStructure()$clustResProbes, annotation = session$userData$annotation, distanceToLook = 10, numCores = session$userData$numCores))
      output$DTDistance10 <- DT::renderDataTable(distNeigboursProbes10(), escape = FALSE, extensions = c("Buttons"), options = list(dom = "Bfrtip", buttons = c("csv"), pageLength = 10000))
      distNeigboursProbes10reduced <- shiny::reactive(na.omit(distNeigboursProbes10()))
      output$DTDistance10reduced <- DT::renderDataTable(distNeigboursProbes10reduced(), escape = FALSE, extensions = c("Buttons"), options = list(dom = "Bfrtip", buttons = c("csv"), pageLength = 10000))
      MinDistance <- shiny::reactive(distNeigboursProbes10()[, 2])
      histMinDistance10 <- shiny::reactive(plotly::plot_ly(x = MinDistance(), type = "histogram", name = "histMinDistance10"))
      output$histMinDistance10 <- plotly::renderPlotly(histMinDistance10())
      MeanDistance <- shiny::reactive(distNeigboursProbes10()[, 3])
      histMeanDistance10 <- shiny::reactive(plotly::plot_ly(x = MeanDistance(), type = "histogram", name = "histMeanDistance10"))
      output$histMeanDistance10 <- plotly::renderPlotly(histMeanDistance10())
      MaxDistance <- shiny::reactive(distNeigboursProbes10()[, 4])
      histMaxDistance10 <- shiny::reactive(plotly::plot_ly(x = MaxDistance(), type = "histogram", name = "histMaxDistance10"))
      output$histMaxDistance10 <- plotly::renderPlotly(histMaxDistance10())

      distNeigboursProbes100 <- shiny::reactive(calculateDistanceNeigboursProbes(wd = session$userData$packageWd, clustResProbes = probeReducedDataStructure()$clustResProbes, annotation = session$userData$annotation, distanceToLook = 100, numCores = session$userData$numCores))
      output$DTDistance100 <- DT::renderDataTable(distNeigboursProbes100(), escape = FALSE, extensions = c("Buttons"), options = list(dom = "Bfrtip", buttons = c("csv"), pageLength = 10000))
      distNeigboursProbes100reduced = shiny::reactive(na.omit(distNeigboursProbes100()))
      output$DTDistance100reduced <- DT::renderDataTable(distNeigboursProbes100reduced(), escape = FALSE, extensions = c("Buttons"), options = list(dom = "Bfrtip", buttons = c("csv"), pageLength = 10000))
      MinDistance <- shiny::reactive(distNeigboursProbes100()[, 2])
      histMinDistance100 <- shiny::reactive(plotly::plot_ly(x = MinDistance(), type = "histogram", name = "histMinDistance100"))
      output$histMinDistance100 <- plotly::renderPlotly(histMinDistance100())
      MeanDistance <- shiny::reactive(distNeigboursProbes100()[, 3])
      histMeanDistance100 <- shiny::reactive(plotly::plot_ly(x = MeanDistance(), type = "histogram", name = "histMeanDistance100"))
      output$histMeanDistance100 <- plotly::renderPlotly(histMeanDistance100())
      MaxDistance <- shiny::reactive(distNeigboursProbes100()[, 4])
      histMaxDistance100 <- shiny::reactive(plotly::plot_ly(x = MaxDistance(), type = "histogram", name = "histMaxDistance100"))
      output$histMaxDistance100 <- plotly::renderPlotly(histMaxDistance100())

      distNeigboursProbes1000 <- shiny::reactive(calculateDistanceNeigboursProbes(wd = session$userData$packageWd, clustResProbes = probeReducedDataStructure()$clustResProbes, annotation = session$userData$annotation, distanceToLook = 1000, numCores = session$userData$numCores))
      output$DTDistance1000 <- DT::renderDataTable(distNeigboursProbes1000(), escape = FALSE, extensions = c("Buttons"), options = list(dom = "Bfrtip", buttons = c("csv"), pageLength = 10000))
      distNeigboursProbes1000reduced = shiny::reactive(na.omit(distNeigboursProbes1000()))
      output$DTDistance1000reduced <- DT::renderDataTable(distNeigboursProbes1000reduced(), escape = FALSE, extensions = c("Buttons"), options = list(dom = "Bfrtip", buttons = c("csv"), pageLength = 10000))
      MinDistance <- shiny::reactive(distNeigboursProbes1000()[, 2])
      histMinDistance1000 <- shiny::reactive(plotly::plot_ly(x = MinDistance(), type = "histogram", name = "histMinDistance1000"))
      output$histMinDistance1000 <- plotly::renderPlotly(histMinDistance1000())
      MeanDistance <- shiny::reactive(distNeigboursProbes1000()[, 3])
      histMeanDistance1000 <- shiny::reactive(plotly::plot_ly(x = MeanDistance(), type = "histogram", name = "histMeanDistance1000"))
      output$histMeanDistance1000 <- plotly::renderPlotly(histMeanDistance1000())
      MaxDistance <- shiny::reactive(distNeigboursProbes1000()[, 4])
      histMaxDistance1000 <- shiny::reactive(plotly::plot_ly(x = MaxDistance(), type = "histogram", name = "histMaxDistance1000"))
      output$histMaxDistance1000 <- plotly::renderPlotly(histMaxDistance1000())

      distNeigboursProbes10000 <- shiny::reactive(calculateDistanceNeigboursProbes(wd = session$userData$packageWd, clustResProbes = probeReducedDataStructure()$clustResProbes, annotation = session$userData$annotation, distanceToLook = 10000, numCores = session$userData$numCores))
      output$DTDistance10000 <- DT::renderDataTable(distNeigboursProbes10000(), escape = FALSE, extensions = c("Buttons"), options = list(dom = "Bfrtip", buttons = c("csv"), pageLength = 10000))
      distNeigboursProbes10000reduced = shiny::reactive(na.omit(distNeigboursProbes10000()))
      output$DTDistance10000reduced <- DT::renderDataTable(distNeigboursProbes10000reduced(), escape = FALSE, extensions = c("Buttons"), options = list(dom = "Bfrtip", buttons = c("csv"), pageLength = 10000))
      MinDistance <- shiny::reactive(distNeigboursProbes10000()[, 2])
      histMinDistance10000 <- shiny::reactive(plotly::plot_ly(x = MinDistance(), type = "histogram", name = "histMinDistance10000"))
      output$histMinDistance10000 <- plotly::renderPlotly(histMinDistance10000())
      MeanDistance <- shiny::reactive(distNeigboursProbes10000()[, 3])
      histMeanDistance10000 <- shiny::reactive(plotly::plot_ly(x = MeanDistance(), type = "histogram", name = "histMeanDistance10000"))
      output$histMeanDistance10000 <- plotly::renderPlotly(histMeanDistance10000())
      MaxDistance <- shiny::reactive(distNeigboursProbes10000()[, 4])
      histMaxDistance10000 <- shiny::reactive(plotly::plot_ly(x = MaxDistance(), type = "histogram", name = "histMaxDistance10000"))
      output$histMaxDistance10000 <- plotly::renderPlotly(histMaxDistance10000())

      shiny::observe({
          output$plotDendrogramProbesLong <- shiny::renderPlot(getPlot(probeReducedDataStructure()$probeDendrogram))
      })
#      #works not:
#      output$plotDendrogramProbesLong <- shiny::reactive(shiny::renderPlot(getPlot(traitReducedDataStructure()$traitDendrogram)))

      shiny::observe({
          output$plotClustergramProbesLong <- shiny::renderPlot(getPlot(probeReducedDataStructure()$probeClustergram))
          #output$plotClustergramTraitsLong <- plotly::renderPlotly(getPlot(probeReducedDataStructure()$probeClustergram))
      })

      shiny::observe({
          output$DTProbesMedoids <- DT::renderDataTable(as.data.frame(getMedoidsTable(id, probeReducedDataStructure()$probeClusterMedoids)))
          #output$DTTraitsMedoids <- DT::renderDataTable(as.data.frame(mtcars))
      })

      shiny::observe({
        output$DTProbesClusters <- DT::renderDataTable(as.data.frame(getClustersTable(id, probeReducedDataStructure()$probeClusters,
                                                                                      probeReducedDataStructure()$probeClusterMedoids)))
      })
      shiny::observe({
        output$DTAnnotatedProbes <- DT::renderDataTable(as.data.frame(getAnnotatedProbeTable(id, probeReducedDataStructure()$probeClusterMedoids, session)))
      })

    },
    error = function(e) {
      base::message("An error occurred in moduleServer in ClusteringProbes_SERVER:\n", e)
      browser() #should not happen
    },
    warning = function(w) {
      base::message("An error occurred in moduleServer in ClusteringProbes_SERVER:\n", w)
      browser() #should not happen
    },
    finally = {
      base::print(base::paste0(sysTimePID(), " finished moduleServer in ClusteringProbes_SERVER"))
      #return(result)
    })
  } #end moduleServer
)} #end ClusteringProbes_SERVER

#' getTraitReducedDataStructure
#' delivers clustering results for traits
#' @param id either "PVal" or "LogFC"
#' @param input shiny input object
#' @param output shiny output object
#' @param session shiny session object
#' @return nested list with results of trait clustering
#' examples getTraitReducedDataStructure(id, input, output, session)
getTraitReducedDataStructure <- function(id, input, output, session) {
  shinyId <- shiny::showNotification("Creating trait reduced data structure...", duration = NULL, closeButton = FALSE)
  on.exit(shiny::removeNotification(shinyId), add = TRUE)
  base::tryCatch(
    {
      if (is.valid(session$userData$sessionVariables$pReducedData())) {
        result <- base::list(combinedDFP_Val_Labels = NULL,
                             #matP_Val = session$userData$sessionVariables$matP_Val(), #is part of combinedDFP_Val_Labels()
                             matP_Val.t = NULL,
                             distMatTraits = NULL,
                             clustResTraits = NULL,
                             traitClusters = NULL,
                             traitClusterMedoids = NULL,
                             traitDendrogram = NULL,
                             traitClustergram = NULL
        )
        result$combinedDFP_Val_Labels <- session$userData$sessionVariables$pReducedData()
        #check whats going on with dfVolcano... everything is fine, except too few cases from debug mode
        if (id == "PVal") {
          mat <- as.matrix(result$combinedDFP_Val_Labels$dfP_Val)
        }
        else if (id == "LogFC") {
          mat <- as.matrix(result$combinedDFP_Val_Labels$dfLogFC)
        }
        mat <- na.omit(mat)
        gc()
        result$matP_Val.t<- t(mat)

        numberCores <- session$userData$numCores
        base::print(paste0(sysTimePID(), " (traitReducedDataStructure) before distance matrix for n(reduced traits) = ", base::nrow(result$matP_Val.t), " (takes some time). Using n(cores) = ", numberCores, "."))
        if (is.valid(result$matP_Val.t)) {
          result$distMatTraits <- getDistMat(numberCores = numberCores, matrix = result$matP_Val.t)
        }
        else {
          result$distMatTraits <- NULL
        }
        base::print(paste0(sysTimePID(), " (traitReducedDataStructure) after distance matrix for reduced traits."))
        #for unknown reason getClustResFast() crashes, if executed without Sys.sleep in advance...
        Sys.sleep(1)
        base::print(paste0(sysTimePID(), " (traitReducedDataStructure) before clustering for traits.", base::nrow(result$matP_Val.t)))
        if (is.valid(result$distMatTraits)) {
          result$clustResTraits <- getClustResFast(result$distMatTraits)
        }
        else {
          result$clustResTraits <- NULL
        }

        base::print(paste0(sysTimePID(), " (traitReducedDataStructure) after clustering results for traits."))
        base::print(base::paste0(sysTimePID(), " (traitReducedDataStructure) start generating traitClusters."))
        maxTraits <- dendextend::nleaves(result$clustResTraits)
        numClusters <- session$userData$sessionVariables$numClustersTraits()
        numClusters <- min(numClusters, maxTraits)
        if (is.valid(result$clustResTraits) && numClusters > 1) {
          result$traitClusters <- dendextend::cutree(result$clustResTraits,
                                                     k = numClusters)
        }
        else {
          result$traitClusters <- NULL
        }
        base::print(base::paste0(sysTimePID(), " (traitReducedDataStructure) start generating traitClusterMedoids."))
        if (is.valid(result$clustResTraits) && is.valid(result$distMatTraits)) {
          result$traitClusterMedoids <- getClusterMedoids(distMat = result$distMatTraits,
                                                          clustRes = result$clustResTraits,
                                                          numClusters = numClusters)
        }
        else {
          result$traitClusterMedoids <- NULL
        }

        #clustergram needs to be generated before reducing to the medoids in order to show all traits
        base::print(base::paste0(sysTimePID(), " (traitReducedDataStructure) start generating traitClustergram."))
        if (is.valid(result$matP_Val.t) && is.valid(result$clustResTraits)) {
          result$traitClustergram <- getClustergramLong(mat = result$matP_Val.t,
                                                        clustRes = result$clustResTraits,
                                                        nClusters = numClusters)
        }
        else {
          result$traitClustergram <- NULL
        }
        if (is.valid(result$traitClusterMedoids)) {
          #select only medoid CpG for further analyses...
          medoids <- which(colnames(result$combinedDFP_Val_Labels$dfP_Val) %in% unlist(result$traitClusterMedoids))
          result$combinedDFP_Val_Labels$dfP_Val <- result$combinedDFP_Val_Labels$dfP_Val[, medoids]
          result$combinedDFP_Val_Labels$dfDM <- result$combinedDFP_Val_Labels$dfDM[, medoids]
          result$combinedDFP_Val_Labels$dfN <- result$combinedDFP_Val_Labels$dfN[, medoids]
          result$combinedDFP_Val_Labels$dfLogFC <- result$combinedDFP_Val_Labels$dfLogFC[, medoids]

          #also shorten:
          result$combinedDFP_Val_Labels$traitID <- result$combinedDFP_Val_Labels$traitID[medoids]
          result$combinedDFP_Val_Labels$mergedOriginDF <- result$combinedDFP_Val_Labels$mergedOriginDF[medoids]
          result$combinedDFP_Val_Labels$mergedColnames <- result$combinedDFP_Val_Labels$mergedColnames[medoids]
          result$combinedDFP_Val_Labels$mergedOriginalColnames <- result$combinedDFP_Val_Labels$mergedOriginalColnames[medoids]
          result$combinedDFP_Val_Labels$mergedOriginTrait <- result$combinedDFP_Val_Labels$mergedOriginTrait[medoids]

          if (length(result$traitClusterMedoids) < ncol(mat)) { #do this only if we have less medoids than before
            #clustering again
            if (id == "PVal") {
              mat <- as.matrix(result$combinedDFP_Val_Labels$dfP_Val)
            }
            else if (id == "LogFC") {
              mat <- as.matrix(result$combinedDFP_Val_Labels$dfP_Val)
            }
            mat <- na.omit(mat)
            gc()
            result$matP_Val.t<- t(mat)
            # clustering for traits
            base::print(base::paste0(sysTimePID(), " (traitReducedDataStructure) before distance matrix for medoids n(traits) = ", base::ncol(mat), " (takes some time). Using n(cores) = ", numberCores, "."))
            if (is.valid(result$matP_Val.t)) {
              result$distMatTraits <- getDistMat(numberCores = numberCores, matrix = result$matP_Val.t)
            }
            else {
              result$distMatTraits <- NULL
            }
            base::print(base::paste0(sysTimePID(), " (traitReducedDataStructure start clustResTraits for medoids."))
            if (is.valid(result$distMatTraits)) {
              result$clustResTraits <- getClustResFast(result$distMatTraits)
            }
            else {
              result$clustResTraits <- NULL
              browser() # should not happen
            }
          }
        }
        else {
          result$clustResTraits <- NULL
          browser() # should not happen
        }
        base::print(base::paste0(sysTimePID(), " (traitReducedDataStructure)  start generating traitDendrogram."))
        #do this only, if a dataset is already loaded
        if (is.valid(result$clustResTraits)) {
          result$traitDendrogram <- getDend(clustRes = result$clustResTraits, nClusters = numClusters)
          #check labels: labels(result$traitDendrogram)
          base::print(base::paste0(sysTimePID(), " (traitReducedDataStructure) after making traitDendrogram."))
        }
        else {
          result$traitDendrogram <- NULL
        }

        #keep original order for all matrices to plot HM (orders data itself again by dendrogram)
        result$combinedDFP_Val_Labels$dfP_Val_Original <- result$combinedDFP_Val_Labels$dfP_Val
        result$combinedDFP_Val_Labels$dfDM_Original <- result$combinedDFP_Val_Labels$dfDM
        result$combinedDFP_Val_Labels$dfN_Original <- result$combinedDFP_Val_Labels$dfN
        result$combinedDFP_Val_Labels$dfLogFC_Original <- result$combinedDFP_Val_Labels$dfLogFC
        if (is.valid(result$clustResTraits)) {
          # add "number" and reorder columns; order comes from result$clustResTraits
          result$combinedDFP_Val_Labels$dfP_Val <- result$combinedDFP_Val_Labels$dfP_Val[, result$clustResTraits$order]
          result$combinedDFP_Val_Labels$dfDM <- result$combinedDFP_Val_Labels$dfP_Val[, result$clustResTraits$order]
          result$combinedDFP_Val_Labels$dfN <- result$combinedDFP_Val_Labels$dfN[, result$clustResTraits$order]
          result$combinedDFP_Val_Labels$dfLogFC <- result$combinedDFP_Val_Labels$dfLogFC[, result$clustResTraits$order]
        }
        else {
          browser() #should not happen
        }
      }
    },
    error = function(e) {
      base::message("An error occurred in getTraitReducedDataStructure():\n", e)
      browser() #should not happen
    },
    warning = function(w) {
      base::message("A warning occurred in getTraitReducedDataStructure():\n", w)
      browser() #should not happen
    },
    finally = {
      base::print(base::paste0(sysTimePID(), " finished getTraitReducedDataStructure()."))
      return(result)
    }
  )
}

#' getProbeReducedDataStructure
#' delivers clustering results for probes as well as probes that are within a defined range (DNAdistances) around each selected probe
#' @param id either "PVal" or "LogFC"
#' @param input shiny input object
#' @param output shiny output object
#' @param session shiny session object
#' @return nested list with results of probe clustering
#' examples getProbeReducedDataStructure(id, input, output, session)
getProbeReducedDataStructure <- function(id, input, output, session) {
  shinyId <- shiny::showNotification("Creating probe reduced data structure...", duration = NULL, closeButton = FALSE)
  on.exit(shiny::removeNotification(shinyId), add = TRUE)
  base::tryCatch({
    if (id == "PVal") {
      p <- session$userData$sessionVariables$traitReducedDataStructurePVal()
    }
    else if (id == "LogFC") {
      p <- session$userData$sessionVariables$traitReducedDataStructureLogFC()
    }
    #if (is.valid(session$userData$sessionVariables$traitReducedDataStructurePVal()$combinedDFP_Val_Labels$dfP_Val)) {
    if (TRUE) {
      # p <- session$userData$sessionVariables$traitReducedDataStructurePVal()
      result <- base::list(combinedDFP_Val_Labels = p$combinedDFP_Val_Labels,
                           matP_Val.t = p$matP_Val.t,
                           distMatTraits = p$distMatTraits,
                           clustResTraits = p$clustResTraits,
                           traitClusters = p$traitClusters,
                           traitClusterMedoids = p$traitClusterMedoids,
                           traitDendrogram = p$traitDendrogram,
                           traitClustergram = p$traitClustergram,
                           distMatProbes = NULL,
                           clustResProbes = NULL,
                           probeClusters = NULL,
                           probeClusterMedoids = NULL,
                           probeDendrogram = NULL,
                           probeClustergram = NULL,
                           DNAdistances = NULL,
                           dfVolcano = NULL,
                           dfKeyShadow = NULL
      )
      p <- NULL
      #calculate clustering results
      base::print(base::paste0(sysTimePID(), " (probeReducedDataStructure) start generating distMatProbes."))
      if (id == "PVal") {
        mat <- as.matrix(result$combinedDFP_Val_Labels$dfP_Val)
      }
      else if (id == "LogFC") {
        mat <- as.matrix(result$combinedDFP_Val_Labels$dfLogFC)
      }
      mat <- na.omit(mat)
      # read slider for number probes
      nProbes <- session$userData$sessionVariables$numClustersProbes()
      #if (is.valid(dfP_Val) && base::nrow(dfP_Val) >= 5) {
      if (is.valid(mat) && base::nrow(mat) >= 5) {
        # dflogFC[dfP_Val > 0.05] <- NA # 1
        base::print(base::paste0(sysTimePID(), " (probeReducedDataStructure) calculating distance matrix with rows= ", nrow(mat), " cols= ", ncol(mat)))
        base::print(base::class(mat))
        # base::print(base::paste0(sysTimePID(), " set missing p-values to 1."))
        # dflogFC[base::is.na(dflogFC)] <- 1 # set missing P_VAL to 1
        base::print(Cstack_info())
        numberCores <- session$userData$numCores
        # clustering for probes
        base::print(base::paste0(sysTimePID(), " (probeReducedDataStructure) before distance matrix for n(probes) = ", base::nrow(mat), " (takes some time). Using n(cores) = ", numberCores, "."))
        gc()
        result$distMatProbes <- getDistMat(numberCores = numberCores, matrix = mat)
        base::print(base::paste0(sysTimePID(), " (probeReducedDataStructure) after distance matrix for probes.", base::nrow(mat)))
        base::print(base::paste0(sysTimePID(), " (probeReducedDataStructure) start clustResProbes."))
        distMat <- result$distMatProbes
        if (is.valid(distMat)) {
          result$clustResProbes <- getClustResFast(distMat)
        }
        else {
          result$clustResProbes <- NULL
          browser() # should not happen
        }
        maxProbes <- dendextend::nleaves(result$clustResProbes)
        numClusters <- session$userData$sessionVariables$numClustersProbes()
        numClusters <- min(numClusters, maxProbes)

        if (is.valid(result$clustResProbes) && numClusters > 1) {
          result$probeClusters <- dendextend::cutree(result$clustResProbes,
                                                     k = numClusters)
        }
        else {
          result$probeClusters <- NULL
        }

        base::print(base::paste0(sysTimePID(), " (probeReducedDataStructure) start generating probeClusterMedoids."))
        if (is.valid(result$clustResProbes) && is.valid(result$distMatProbes)) {
          result$probeClusterMedoids <- getClusterMedoids(distMat = result$distMatProbes,
                                                          clustRes = result$clustResProbes,
                                                          numClusters = numClusters)
        }
        else {
          result$probeClusterMedoids <- NULL
        }

        base::print(base::paste0(sysTimePID(), " (probeReducedDataStructure) start generating probeClustergram."))
        #if (is.valid(result$matP_Val.t) && is.valid(result$clustResProbes)) {
        if (is.valid(result$clustResProbes)) {
          result$probeClustergram <- getClustergramLong(mat = mat,
                                                        clustRes = result$clustResProbes,
                                                        nClusters = numClusters)
        }
        else {
          result$probeClustergram <- NULL
        }
        if (is.valid(result$probeClusterMedoids)) {
          #select only medoid CpG for further analyses...
          medoids <- which(rownames(result$combinedDFP_Val_Labels$dfP_Val) %in% unlist(result$probeClusterMedoids))
          result$combinedDFP_Val_Labels$dfP_Val <- result$combinedDFP_Val_Labels$dfP_Val[medoids, ]
          result$combinedDFP_Val_Labels$dfDM <- result$combinedDFP_Val_Labels$dfDM[medoids, ]
          result$combinedDFP_Val_Labels$dfN <- result$combinedDFP_Val_Labels$dfN[medoids, ]
          result$combinedDFP_Val_Labels$dfLogFC <- result$combinedDFP_Val_Labels$dfLogFC[medoids, ]
          #do this also for the unordered _Original data structures
          result$combinedDFP_Val_Labels$dfP_Val_Original <- result$combinedDFP_Val_Labels$dfP_Val_Original[medoids, ]
          result$combinedDFP_Val_Labels$dfDM_Original <- result$combinedDFP_Val_Labels$dfDM_Original[medoids, ]
          result$combinedDFP_Val_Labels$dfN_Original <- result$combinedDFP_Val_Labels$dfN_Original[medoids, ]
          result$combinedDFP_Val_Labels$dfLogFC_Original <- result$combinedDFP_Val_Labels$dfLogFC_Original[medoids, ]

          if (length(result$probeClusterMedoids) < nrow(mat)) { #do this only if we have less medoids than before
            #clustering again
            if (id == "PVal") {
              mat <- as.matrix(result$combinedDFP_Val_Labels$dfP_Val)
            }
            else if (id == "LogFC") {
              mat <- as.matrix(result$combinedDFP_Val_Labels$dfLogFC)
            }
            mat <- na.omit(mat)
            # clustering for probes
            base::print(base::paste0(sysTimePID(), " (probeReducedDataStructure) before distance matrix for medoids n(probes) = ", base::nrow(mat), " (takes some time). Using n(cores) = ", numberCores, "."))
            gc()
            result$distMatProbes <- getDistMat(numberCores = numberCores, matrix = mat)

            base::print(base::paste0(sysTimePID(), " (probeReducedDataStructure) start clustResProbes for medoids."))
            distMat <- result$distMatProbes
            if (is.valid(distMat)) {
              result$clustResProbes <- getClustResFast(distMat)
            }
            else {
              result$clustResProbes <- NULL
              browser() # should not happen
            }
          }
        }
        else {
          result$clustResProbes <- NULL
          browser() # should not happen
        }
        base::print(base::paste0(sysTimePID(), " (probeReducedDataStructure) start generating probeDendrogram."))
        if (is.valid(result$clustResProbes)) {
          #reorder columns; order comes from result$clustResProbes
          result$combinedDFP_Val_Labels$dfP_Val <- result$combinedDFP_Val_Labels$dfP_Val[result$clustResProbes$order, ]
          result$combinedDFP_Val_Labels$dfDM <- result$combinedDFP_Val_Labels$dfP_Val[result$clustResProbes$order, ]
          result$combinedDFP_Val_Labels$dfN <- result$combinedDFP_Val_Labels$dfN[result$clustResProbes$order, ]
          result$combinedDFP_Val_Labels$dfLogFC <- result$combinedDFP_Val_Labels$dfLogFC[result$clustResProbes$order, ]

          result$probeDendrogram <- getDend(clustRes = result$clustResProbes, nClusters = numClusters)
          #check labels: labels(result$probeDendrogram)
          base::print(base::paste0(sysTimePID(), " (probeReducedDataStructure) after making probeDendrogram."))
        }
        else {
          result$probeDendrogram <- NULL
        }

        result$combinedDFP_Val_Labels$dfP_Val_w_number <- result$combinedDFP_Val_Labels$dfP_Val
        nprobes <- nrow(result$combinedDFP_Val_Labels$dfP_Val_w_number)
        result$combinedDFP_Val_Labels$dfP_Val_w_number$number <- seq(1:nprobes)
        col_order <- c("number", colnames(result$combinedDFP_Val_Labels$dfP_Val_w_number))
        result$combinedDFP_Val_Labels$dfP_Val_w_number <- result$combinedDFP_Val_Labels$dfP_Val_w_number[, col_order]
        result$combinedDFP_Val_Labels$dfP_Val_w_number <- result$combinedDFP_Val_Labels$dfP_Val_w_number[ , -which(colnames(result$combinedDFP_Val_Labels$dfP_Val_w_number) %in% "number.1")]

        result$combinedDFP_Val_Labels$dfDM_w_number <- result$combinedDFP_Val_Labels$dfDM
        nprobes <- nrow(result$combinedDFP_Val_Labels$dfDM_w_number)
        result$combinedDFP_Val_Labels$dfDM_w_number$number <- seq(1:nprobes)
        col_order <- c("number", colnames(result$combinedDFP_Val_Labels$dfDM_w_number))
        result$combinedDFP_Val_Labels$dfDM_w_number <- result$combinedDFP_Val_Labels$dfDM_w_number[, col_order]
        result$combinedDFP_Val_Labels$dfDM_w_number <- result$combinedDFP_Val_Labels$dfDM_w_number[ , -which(colnames(result$combinedDFP_Val_Labels$dfDM_w_number) %in% "number.1")]

        result$combinedDFP_Val_Labels$dfN_w_number <- result$combinedDFP_Val_Labels$dfN
        nprobes <- nrow(result$combinedDFP_Val_Labels$dfN_w_number)
        result$combinedDFP_Val_Labels$dfN_w_number$number <- seq(1:nprobes)
        col_order <- c("number", colnames(result$combinedDFP_Val_Labels$dfN_w_number))
        result$combinedDFP_Val_Labels$dfN_w_number <- result$combinedDFP_Val_Labels$dfN_w_number[, col_order]
        result$combinedDFP_Val_Labels$dfN_w_number <- result$combinedDFP_Val_Labels$dfN_w_number[ , -which(colnames(result$combinedDFP_Val_Labels$dfN_w_number) %in% "number.1")]

        result$combinedDFP_Val_Labels$dfLogFC_w_number <- result$combinedDFP_Val_Labels$dfLogFC
        nprobes <- nrow(result$combinedDFP_Val_Labels$dfLogFC_w_number)
        result$combinedDFP_Val_Labels$dfLogFC_w_number$number <- seq(1:nprobes)
        col_order <- c("number", colnames(result$combinedDFP_Val_Labels$dfLogFC_w_number))
        result$combinedDFP_Val_Labels$dfLogFC_w_number <- result$combinedDFP_Val_Labels$dfLogFC_w_number[, col_order]
        result$combinedDFP_Val_Labels$dfLogFC_w_number <- result$combinedDFP_Val_Labels$dfLogFC_w_number[ , -which(colnames(result$combinedDFP_Val_Labels$dfLogFC_w_number) %in% "number.1")]

        Distance <- session$userData$sessionVariables$numNeighbours()
        DNAdistances <- calculateDistanceNeigboursProbes(wd = session$userData$packageWd, clustResProbes = result$clustResProbes, annotation = session$userData$annotation, distanceToLook = Distance, numCores = numberCores)
        result$DNAdistances <- DNAdistances
        dfLogFC <- result$combinedDFP_Val_Labels$dfLogFC
        if (is.valid(dfLogFC)) {
          dfP_Val <- result$combinedDFP_Val_Labels$dfP_Val
          if (is.valid(dfP_Val)) {
            #take everything into one table with columns p-val and logFC ...
            dfLogFC$probe <- row.names(dfLogFC)
            dfLogFC <- tidyr::pivot_longer(dfLogFC, cols  = -probe, names_to = c("trait"))
            colnames(dfLogFC)[3] <- "LogFC"
            dfP_Val$probe <- row.names(dfP_Val) #feed in traitLabels here
            dfP_Val <- tidyr::pivot_longer(dfP_Val, cols  = -probe, names_to = c("trait"))
            colnames(dfP_Val)[3] <- "P_Val"
            dfVolcano <- base::merge(dfP_Val, dfLogFC, by.x = c("probe","trait"), by.y = c("probe","trait"), all.x = TRUE, all.y = FALSE)
            #merge chr and position
            annotation <- base::subset(session$userData$annotation, select = c("name", "chromosome", "position", "gene.symbol"))
            dfVolcano <- base::merge(dfVolcano, annotation, by.x = "probe", by.y = "name", all.x = FALSE, all.y = FALSE)
            #add distances to dfVolcano
            DNAdistances <- result$DNAdistances
            row.names(DNAdistances) <- DNAdistances$ID
            DNAdistances$cg <- DNAdistances$ID #row.names(DNAdistances)
            DNAdistancesNumber <- DNAdistances[,c("cg", "number")]
            DNAdistancesNumber <- tidyr::pivot_longer(DNAdistancesNumber, cols  = -cg, names_to = c("number"))
            DNAdistancesNumber <- DNAdistancesNumber[,c("cg", "value")]
            colnames(DNAdistancesNumber)[2] <- "DistanceNumber"
            DNAdistancesMin <- DNAdistances[,c("cg", "minDistance")]
            DNAdistancesMin <- tidyr::pivot_longer(DNAdistancesMin, cols  = -cg, names_to = c("min"))
            DNAdistancesMin <- DNAdistancesMin[,c("cg", "value")]
            colnames(DNAdistancesMin)[2] <- "DistanceMin"

            DNAdistancesMean <- DNAdistances[,c("cg", "meanDistance")]
            DNAdistancesMean <- tidyr::pivot_longer(DNAdistancesMean, cols  = -cg, names_to = c("mean"))
            DNAdistancesMean <- DNAdistancesMean[,c("cg", "value")]
            colnames(DNAdistancesMean)[2] <- "DistanceMean"

            DNAdistancesMax <- DNAdistances[,c("cg", "maxDistance")]
            DNAdistancesMax <- tidyr::pivot_longer(DNAdistancesMax, cols  = -cg, names_to = c("max"))
            DNAdistancesMax <- DNAdistancesMax[,c("cg", "value")]
            colnames(DNAdistancesMax)[2] <- "DistanceMax"

            dfVolcano <- base::merge(dfVolcano, DNAdistancesNumber, by.x = c("probe"), by.y = c("cg"), all.x = TRUE, all.y = FALSE)
            dfVolcano <- base::merge(dfVolcano, DNAdistancesMin, by.x = c("probe"), by.y = c("cg"), all.x = TRUE, all.y = FALSE)
            dfVolcano <- base::merge(dfVolcano, DNAdistancesMean, by.x = c("probe"), by.y = c("cg"), all.x = TRUE, all.y = FALSE)
            dfVolcano <- base::merge(dfVolcano, DNAdistancesMax, by.x = c("probe"), by.y = c("cg"), all.x = TRUE, all.y = FALSE)
            #create shadow table for key assignment
            #add probe trait and traitsource to key
            originTrait <- result$combinedDFP_Val_Labels$mergedOriginTrait
            originTrait <- rep(originTrait, nprobes)
            keys <- seq(1:nrow(dfVolcano))
            dfKeyShadow <- base::data.frame(key = keys)
            dfKeyShadow$probe <- dfVolcano$probe
            traitLabels <- result$combinedDFP_Val_Labels$mergedOriginalColnames
            traitLabels <- rep(traitLabels, nprobes)
            dfKeyShadow$trait <- traitLabels #dfVolcano$traitLabel
            traitID <- result$combinedDFP_Val_Labels$traitID
            traitID <- rep(traitID, nprobes)
            dfKeyShadow$traitID <- traitID
            dfKeyShadow$traitSource <- originTrait #this one works for clustering by p-value (ln 1580), but not for clustering by log(FC), originTrait passt, dfVolcano ist zu klein 75268...
            dfVolcano$key <- dfKeyShadow$key
            rownames(dfVolcano) <- dfVolcano$key
            rownames(dfKeyShadow) <- dfKeyShadow$key
            #sort by p-val and logFC
            dfVolcano <- dfVolcano[base::order(dfVolcano$P_Val, dfVolcano$LogFC, decreasing = c(FALSE, TRUE), na.last = c(TRUE,TRUE)),]
            #take only first 2^16 entities to be able to plot volcano plot using plotly
            if(nrow(dfVolcano)>2^16) {
              dfVolcano <- dfVolcano[1:2^16-1,]
              dfKeyShadow <- dfKeyShadow[1:2^16-1,]
            }
            result$dfVolcano <- dfVolcano
            result$dfKeyShadow <- dfKeyShadow
          }
          else {
            result$dfVolcano <- NULL
            browser() # should not happen
          }
        }
        else {
          result$dfVolcano <- NULL
          browser() # should not happen
        }
# browser() #limit results to those whith distances within
#         dfP_Val <- result$combinedDFP_Val_Labels$dfP_Val
#         dfP_Val <- dfP_Val[which(rownames(dfP_Val) %in% DNAdistances$ID), ]
#         result$combinedDFP_Val_Labels$dfP_Val <- dfP_Val
#         rm(dfP_Val)
#         dfDM <- result$combinedDFP_Val_Labels$dfDM
#         dfDM <- dfDM[which(rownames(dfDM) %in% DNAdistances$ID), ]
#         result$combinedDFP_Val_Labels$dfDM <- dfDM
#         rm(dfDM)
#         dfLogFC <- result$combinedDFP_Val_Labels$dfLogFC
#         #check for missings (should be there) and negative values (too), yes we have negatives and missings here
#         dfLogFC <- dfLogFC[which(rownames(dfLogFC) %in% DNAdistances$ID), ]
#         result$combinedDFP_Val_Labels$dfLogFC <- dfLogFC
#         rm(dfLogFC)
#         dfN <- result$combinedDFP_Val_Labels$dfN
#         dfN <- dfN[which(rownames(dfN) %in% DNAdistances$ID), ]
#         result$combinedDFP_Val_Labels$dfN <- dfN
#         rm(dfN)
      }
      else {
        result <- NULL
        base::message(base::paste0(sysTimePID(), " (probeReducedDataStructure) number of rows in p-val table too low: ", base::nrow(mat)))
        browser() # should not happen
      }
      }
    },
    error = function(e) {
      base::message("An error occurred in getProbeReducedDataStructure():\n", e)
      browser() #should not happen
    },
    warning = function(w) {
      base::message("A warning occurred in getProbeReducedDataStructure():\n", w)
      browser() #should not happen
    },
    finally = {
      base::print(base::paste0(sysTimePID(), " finished getProbeReducedDataStructure()."))
      return(result)
  })
}

#' getDistMat
#' calculates distance matrix using parallelDist
#' @param numberCores number of CPU cores to use
#' @param matrix matrix to calculate distance matrix for
#' @return distance matrix
#' examples getDistMat(numberCores, matrix)
getDistMat <- function(numberCores, matrix) {
  shinyId <- shiny::showNotification("Calculating distance matrix parallel...", duration = NULL, closeButton = FALSE)
  base::on.exit(shiny::removeNotification(shinyId), add = TRUE)
  #distMat <- spam::as.spam.dist(
  if (length(matrix) > 100000) {
    distMat <- stats::as.dist(
      parallelDist::parallelDist(
        base::as.matrix(matrix),
        method = "euclidean",
        labels = TRUE,
        threads = numberCores
      )
    )
    if (all(is.na(distMat))) {
      browser() #this should not happen, if dist(matrix) delivers such a result, then something is wrong with parallelDist
    }
  }
  else {
    distMat <- stats::dist(base::as.matrix(matrix), method = "euclidean")
  }
  return(distMat)
}

#' getClustResFast
#' calculates hierarchical clustering using fastcluster
#' @param distanceMatrix distance matrix
#' @return hclust object
#' examples getClustResFast(distanceMatrix)
getClustResFast <- function(distanceMatrix) {
  shinyId <- shiny::showNotification("Calculating clustering results...", duration = NULL, closeButton = FALSE)
  base::on.exit(shiny::removeNotification(shinyId), add = TRUE)
  base::print(base::paste0(sysTimePID(), " start getClustResFast()"))
  #check size of distanceMatrix
  if (is.valid(distanceMatrix)) {
    base::print(base::paste0(sysTimePID(), " clustering trait data."))
    #startTime <- Sys.time()
    gc()
    ClustRes <- fastcluster::hclust(stats::as.dist(distanceMatrix), method = "ward.D2")
  }
  else {
    ClustRes <- NULL
    base::message(base::paste0(sysTimePID(), " is.valid(distanceMatrix) == FALSE"))
  }
  return(ClustRes)
}

#' #' updateTxtClusterOut
#' #' generates summary text after clustering
#' #' @param traitReducedcombinedDFP_Val_Labels data structure with trait reduced results
#' #' @param minP_Val minimum p_value for model to use for clustering
#' #' @param maxP_Val maximum p_value for model to use for clustering
#' #' @param minN minimum n for model to use for clustering
#' #' @param sldNumClasses number of classes to use for clustering
#' #' @return text
#' #' examples updateTxtClusterOut(traitReducedcombinedDFP_Val_Labels, minP_Val, maxP_Val, minN, sldNumClasses)
#' updateTxtClusterOut <- function(traitReducedcombinedDFP_Val_Labels, minP_Val, maxP_Val, minN, sldNumClasses) {
#'   base::tryCatch(
#'    {
#'       result <- NULL
#'       if (is.valid(traitReducedcombinedDFP_Val_Labels)) {
#'         maxClasses <- length(traitReducedcombinedDFP_Val_Labels$mergedColnames)
#'         numRow <- nrow(traitReducedcombinedDFP_Val_Labels$dfP_Val)
#'         numCol <- ncol(traitReducedcombinedDFP_Val_Labels$dfP_Val)
#'         minClasses <- 1 #dendextend::min_depth(session$userData$sessionVariables$dendTraits)
#'         result <- base::paste0("finding trait clusters successful. found minClusters = ",
#'                             minClasses, "; maxClusters: ", maxClasses, "; Clustering made for numClasses = ", sldNumClasses, ".\n",
#'                             "size of result df: nrow(CpG)=", numRow, ", ncol(trait)=", numCol, ".")
#'       }
#'     },
#'     error = function(e) {
#'       base::message("An error occurred in updateTxtClusterOut():\n", e)
#'     },
#'     warning = function(w) {
#'       base::message("A warning occurred in updateTxtClusterOut():\n", w)
#'     },
#'     finally = {
#'       return(shiny::HTML(result))
#'     }
#'   )
#' }

#' calculateDistanceNeigboursProbes
#' calculate distance from each probe to its neigbours and gives back data frame with distance metrics
#' @param wd working directory
#' @param clustResProbes data structure with clustering result
#' @param annotation annotation of CpG (names, location etc.)
#' @param distanceToLook maximum distance to look for
#' @param numCores number of cores to use for distance calculation
#' @return data.frame with min, mean, max distance and nuber of CpG in distanceToLook window
#' examples calculateDistanceNeigboursProbes(wd, clustResProbes, annotation, distanceToLook, numCores)
calculateDistanceNeigboursProbes <- function(wd, clustResProbes, annotation, distanceToLook, numCores) {
  shinyId <- shiny::showNotification("Calculating distance neighbours probes...", duration = NULL, closeButton = FALSE)
  base::on.exit(shiny::removeNotification(shinyId), add = TRUE)
  base::tryCatch(
    {
      if (!is.null(clustResProbes)){
        base::print(base::paste0(sysTimePID(), " start calculateDistanceNeigboursProbes() with max distance (distanceToLook) = ", distanceToLook, "."))
        #get chr and location from annotation
        maxDistanceToLook <- distanceToLook
        annotation <- subset(annotation, select = c("name", "chromosome", "position"))
        #merge annotation
        CpG <- data.table::as.data.table(clustResProbes$labels[clustResProbes$order]) #as.data.frame(clustResProbes$labels[clustResProbes$order])
        colnames(CpG)[1] <- "label"
        #CpG$order <- seq(1:nrow(CpG))
        CpG$order <- seq_len(base::nrow(CpG))
        distances <- base::merge(CpG, annotation, by.x = "label", by.y = "name")
        #DNAdistancesUp <- data.table::as.data.table(seq_along(distances[, 2]), 2) #base::data.frame(seq_along(distances[, 2]), 2)
        DNAdistances <- data.table::as.data.table(seq_along(distances[, 2]), 5) #base::data.frame(seq_along(distances[, 2]), 5)
        #sort order given by clustering
        distances <- distances[base::order(distances$order),]
        #library(future) #we have this already in DESCRIPTION file, but without "library(future)" here, it won't work. Strange.
        #library(doFuture)
        #future::plan(strategy = future::multisession, workers = numCores)
        #calculate mean distance to distanceToLook next probes in given order, omit probes with different chr
        #      foreach::foreach(i = seq_along(distances[, 2]), .combine = rbind, .verbose = TRUE) %dofuture% { #for all objects in distances
        foreach::foreach(i = seq_along(distances[, 2])) %do% { #for all objects in distances
        #for(i in seq_along(distances[,2])) { #for all objects in distances
          #base::source(paste0(wd, "/R/Clustering.R")) #this is necessary for foreach %dopar% to run properly
          currentCpG <- distances[i, ]
          #cut out defined number of probes (up- and downstream from CpG)
          lowerBorder <- i-distanceToLook
          if (lowerBorder < 1) {lowerBorder <- 1}#
          upperBorder <- i+distanceToLook
          if (upperBorder > nrow(distances)) {upperBorder <- nrow(distances)}
          distancesToLook <- distances [lowerBorder:upperBorder, ]
          #do not iterate over all probes, but exclude those from different chr first
          chr <- currentCpG$chr
          distancesToLook <- distancesToLook[distancesToLook$chr == chr,]
          distanceToLook <- nrow(distancesToLook)
          if (distanceToLook > 1) {
            DNAdistancesUp <- data.table::as.data.table(seq_along(1:distanceToLook), 2) #base::data.frame(seq_along(1:distanceToLook), 2)
            #upstream
            foreach::foreach(j = 1:distanceToLook) %do% { #max. distance given by distanceToLook
            #for(j in 1:distanceToLook) { #max. distance given by distanceToLook
              CpG <- distances[j, ]
              if (currentCpG$label != CpG$label) {
                DNAdistancesUp[j, 2] <- base::abs(currentCpG$position - CpG$position)
              }
              else {
                DNAdistancesUp[j, 2] <- NA
              }
              DNAdistancesUp[j, 1] <- currentCpG$label
            }
            #check number of nearby CpG
            if (is.numeric(DNAdistancesUp[j, 2])) {
              DNAdistances[i, 2] <- base::min((DNAdistancesUp[, 2]), na.rm = TRUE)
              DNAdistances[i, 3] <- base::mean((DNAdistancesUp[, 2]), na.rm = TRUE)
              DNAdistances[i, 4] <- base::max((DNAdistancesUp[, 2]), na.rm = TRUE)
              DNAdistances[i, 5] <- base::length(na.omit(DNAdistancesUp[, 2]))
            }
            else {
              DNAdistances[i, 2] <- NA
              DNAdistances[i, 3] <- NA
              DNAdistances[i, 4] <- NA
              DNAdistances[i, 5] <- NA
            }
            DNAdistances[i, 1] <- currentCpG$label
          }
          else {
            #no near CpG on the same chromosome found
            DNAdistances[i, 2] <- NA
            DNAdistances[i, 3] <- NA
            DNAdistances[i, 4] <- NA
            DNAdistances[i, 5] <- NA
          }
          DNAdistances[i, 1] <- currentCpG$label
        }
        colnames(DNAdistances) <- c("ID", "minDistance", "meanDistance", "maxDistance", "number")
        distances <- na.omit(DNAdistances)
        base::message(base::paste0(sysTimePID(), " found n = ", nrow(distances), " neigbouring CpG with distance <=", maxDistanceToLook, ""))
#tbc(): list works, but tibble not
#       DNAdistances <- tibble::rownames_to_column(as.data.frame(DNAdistances), var = "rowname")
      }
      else {
        DNAdistances <- NULL
      }
    },
    error = function(e) {
      base::message("An error occurred in calculateDistanceNeigboursProbes():\n", e)
      browser() #should not happen
    },
    warning = function(w) {
      base::message("A warning occurred in calculateDistanceNeigboursProbes():\n", w)
      browser() #should not happen
    },
    finally = {
      base::print(base::paste0(sysTimePID(), " end calculateDistanceNeigboursProbes()."))
      return(DNAdistances)
    }
  )
}

#' getClustergramLong
#' produces a clustergram
#' @param mat transposed matrix with result values
#' @param clustRes (reduced) clustering results
#' @param nClusters number of cluster medoids
#' @return plot with clustergram
#' examples getClustergramLong(mat, clustRes, nClusters)
getClustergramLong <- function(mat, clustRes, nClusters) {
  shinyId <- shiny::showNotification("Creating clustergram long...", duration = NULL, closeButton = FALSE)
  on.exit(shiny::removeNotification(shinyId), add = TRUE)
  base::print(base::paste0(sysTimePID(), " start plotting clustergram."))
  base::tryCatch(
    {
      if (nClusters > 1) {
        Clusters <- dendextend::cutree(clustRes, k = nClusters)
        mycolors <- colorRampPalette(RColorBrewer::brewer.pal(8, "Set2"))(nClusters)
        base::tryCatch(
          {
            p <- factoextra::fviz_cluster(list(data = mat, cluster = Clusters),
                                          palette = mycolors, #palette = "jco",
                                          ggtheme = ggplot2::theme_classic())
          },
          error = function(e) {
            base::message("An error occurred in p <- factoextra::fviz_cluster():\n", e)
            browser() #should not happen
            return(NULL)
          },
          warning = function(w) {
            base::message("An warning occurred in p <- factoextra::fviz_cluster():\n", w)
            browser() #should not happen
            return(NULL)
          },
          finally = {
            return(p)
          }
        )
      }
      else {
        p <- getEmptyPlot()
      }
    },
    error = function(e) {
      base::message("An error occurred in getClustergramLong():\n", e)
    },
    warning = function(w) {
      base::message("A warning occurred in getClustergramLong():\n", w)
    },
    finally = {
      base::print(base::paste0(sysTimePID(), " end plotting clustergram."))
      return(p)
    }
  )
}

#' getClusterMedoids
#' gets Medoids of earlier produced clusters; select trait with lowest (square) sum of distances as cluster medoid
#' @param distMat distance matrix
#' @param clustRes (reduced) clustering results
#' @param numClusters number of clusters to produce
#' @return list with cluster medoids
#' examples getClusterMedoids(distMat, clustRes, numClusters)
getClusterMedoids <- function(distMat, clustRes, numClusters) {
  shinyId <- shiny::showNotification("Getting trait cluster medoids...", duration = NULL, closeButton = FALSE)
  on.exit(shiny::removeNotification(shinyId), add = TRUE)
  base::tryCatch(
    {
      base::print(base::paste0(sysTimePID(), " start find cluster medoids."))
      if (is.valid(clustRes) && numClusters > 1) {
        i <- NULL
        numClusters <- numClusters
        ClustRes <- clustRes
        base::print(base::paste0(sysTimePID(), " cutting tree."))
        Clusters <- cutree(ClustRes, k = numClusters)
        numClusters <- max(Clusters)
        base::print(base::paste0(sysTimePID(), " defined ", numClusters, " clusters."))
        ClusterMedoids <- base::list()
        base::print(base::paste0(sysTimePID(), " finding medoids of n = ", numClusters, " clusters."))
        #iterate over each cluster:
        foreach::foreach(i = 1:numClusters) %do% {
          #select objects in cluster from DistMat
          idx <- which(Clusters == i)
          if (length(idx) > 1) {
            distMatCluster <- usedist::dist_subset(distMat, idx)
            dM <- as.matrix(distMatCluster)
            CS <- colSums(dM^2)
            min(CS)
            #select trait with lowest (square) sum of distances as cluster medoid
            ClusterMedoid <- which(CS == min(CS))
            ClusterMedoid <- labels(ClusterMedoid)[[1]] #take names from ClusterMedoid
          }
          else {
            ClusterMedoid <- ClustRes$labels[[i]] #take name from ClustRes
          }
          ClusterMedoids[[i]] <- ClusterMedoid
        }
      }
      else {
        ClusterMedoids <- NULL
      }
    },
    error = function(e) {
      base::message("An error occurred in getClusterMedoids():\n", e)
    },
    warning = function(w) {
      base::message("A warning occurred in getClusterMedoids():\n", w)
    },
    finally = {
      return(ClusterMedoids)
      base::print(base::paste0(sysTimePID(), " end find cluster medoids."))
    }
  )
}

getAnnotatedProbeTable <- function(id, listProbeMedoids, session) {
  shinyId <- shiny::showNotification("Getting cluster trait table...", duration = NULL, closeButton = FALSE)
  on.exit(shiny::removeNotification(shinyId), add = TRUE)
  base::tryCatch({
    if (is.valid(listProbeMedoids)) {
      result <- listProbeMedoids
      result <- as.data.frame(unlist(result))
      colnames(result) <- NULL
      colnames(result)[[1]] <- "probeID"
      #merge Annotation:
      result$order <- base::seq_len(base::nrow(result))
      rownames(result) <- result$probeID
      # add annotation
      result <-
        base::merge(
          x = result,
          y = session$userData$annotation,
          by.x = "probeID",
          by.y = "name",
          all.x = TRUE,
          all.y = FALSE
        )
      # sort
      result <- result[base::order(result$order), ]
      rownames(result) <- result$probeID
      result <- addLinkToEWASDataHubShort(result, session$userData$config$baseURL_EWASDataHub, session$userData$config$probeAttribut)
      result <- addLinkToMRCEWASCatalogShort(result, session$userData$config$baseURL_MRCEWASCatalog, session$userData$config$probeAttribut)
      result$probeID <- NULL
    }
  },
  error = function(e) {
    base::message("An error occurred in getAnnotatedProbeTable():\n", e)
    browser() #should not happen
  },
  warning = function(w) {
    base::message("A warning occurred in getAnnotatedProbeTable():\n", w)
    browser() #should not happen
  },
  finally = {
    return(result)
    base::print(base::paste0(sysTimePID(), " end getAnnotatedProbeTable()."))
  })
}

#' getClustersTable
#' gets a table earlier produced clusters
#' @param listClusters list with clusters
#' @param listMedoids list with medoids
#' @return data.frame with information on clusters and its medoids
#' examples getClustersTable(listClusters, listMedoids)
getClustersTable <- function(id, listClusters, listMedoids) {
  shinyId <- shiny::showNotification("Getting cluster table...", duration = NULL, closeButton = FALSE)
  on.exit(shiny::removeNotification(shinyId), add = TRUE)
  base::tryCatch({
    if (is.valid(listClusters)) {
      #result <- as.data.frame(listClusters)
      result <- as.data.frame(unlist(listClusters))
      colnames(result)[[1]] <- "Cluster #"
      result["Medoid name"] <- NA
      i <- NULL
      foreach::foreach(i = seq_along(listClusters)) %do% {
        result[i, 2] <- listMedoids[[result[i, 1]]]
      }
    }
  },
  error = function(e) {
    base::message("An error occurred in getClustersTable():\n", e)
    browser() #should not happen
  },
  warning = function(w) {
    base::message("A warning occurred in getClustersTable():\n", w)
    browser() #should not happen
  },
  finally = {
    return(result)
    base::print(base::paste0(sysTimePID(), " end getClustersTable()."))
  })
}

#' getMedoidsTable
#' gets a table with Medoids of earlier produced clusters
#' @param listMedoids list with medoids
#' @return data.frame with information on medoids
#' examples getMedoidsTable(listMedoids)
getMedoidsTable <- function(id, listMedoids) {
  shinyId <- shiny::showNotification("Getting medoids table...", duration = NULL, closeButton = FALSE)
  on.exit(shiny::removeNotification(shinyId), add = TRUE)
  base::tryCatch({
      if (is.valid(listMedoids)) {
        result <- listMedoids
        result <- as.data.frame(unlist(result))
        colnames(result) <- NULL
        colnames(result)[[1]] <- "Medoid"
      }
    },
    error = function(e) {
      base::message("An error occurred in getMedoidsTable():\n", e)
      browser() #should not happen
    },
    warning = function(w) {
      base::message("A warning occurred in getMedoidsTable():\n", w)
      browser() #should not happen
    },
    finally = {
      return(result)
      base::print(base::paste0(sysTimePID(), " end getMedoidsTable()."))
    })
}

getPlot <- function(plotObject) {
  if (is.valid(plotObject)) {
    return(plot(plotObject))
  }
}

#' getDend
#' gets a dendrogram with clustering results
#' @param clustRes clustering results for traits
#' @param nClusters list with traits for coloring
#' @return dendrogram
#' examples getDendTraits(clustRes, nClusters)
getDend <- function(clustRes, nClusters) {
  shinyId <- shiny::showNotification("Getting dendrogram...", duration = NULL, closeButton = FALSE)
  on.exit(shiny::removeNotification(shinyId), add = TRUE)
  base::tryCatch({
    if (is.valid(clustRes) && nClusters > 1) {
      base::print(base::paste0(sysTimePID(), " start getDend()."))
      dend <- stats::as.dendrogram(clustRes)
      #the following fails, if number of objects inside clustResTraits is not between 1 and 120:
      #Error in stats::cutree(tree, k = k, h = h, ...) : elements of 'k' must be between 1 and 120
      #so check
      if (nClusters < 120 && nClusters >= 1) {
        dendTraits <- dendextend::color_branches(dend, k = nClusters)
        dendTraits <- dendextend::color_labels(dend, k = nClusters)
      }
      else {
        message("Warning: number of clusters is not 1 < nClusters < 120\n")
      }
    }
    else {
      dend <- getEmptyPlot()
      #        dend <- plot.new() #plot.window(xlim = c(0 , 1), ylim = c( 5, 10)) #NULL
    }
  },
  error = function(e) {
    base::message("An error occurred in getDend():\n", e)
    browser() #should not happen
  },
  warning = function(w) {
    base::message("A warning occurred in getDend():\n", w)
    browser() #should not happen
  },
  finally = {
    return(dend)
    base::print(base::paste0(sysTimePID(), " end getDend()."))
  })
}

# updateTxtOmitTraitsOut <- function(combinedDFP_Val_Labels) {
#   base::tryCatch(
#     {
#       result <- NULL
#       if (is.valid(combinedDFP_Val_Labels)) {
#         result <- base::paste0("reduce data successful. result table is: nrow (CpG): ",
#                                nrow(combinedDFP_Val_Labels$dfP_Val),
#                                "; ncol (trait): ", ncol(combinedDFP_Val_Labels$dfP_Val), " (if # is less than defined number of clusters, then duplicate variable names have been omitted)")
#       }
#       else {
#         base::message(base::paste0(sysTimePID(), " is.valid(combinedDFP_Val_Labels) == FALSE."))
#       }
#     },
#     error = function(e) {
#       base::message("An error occurred in updateTxtOmitTraitsOut():\n", e)
#     },
#
#     warning = function(w) {
#       base::message("A warning occurred in updateTxtOmitTraitsOut():\n", w)
#     },
#     finally = {
#       return(shiny::HTML(result))
#     }
#   )
# }
#
updateTxtOmitTraitsOut <- function(id, traitReducedDataStructure) {
  base::tryCatch(
    {
      result <- NULL

      if (is.valid(traitReducedDataStructure$combinedDFP_Val_Labels$dfP_Val)) {
        result <- base::paste0("reduce data successful. result table is: nrow (CpG): ",
                               nrow(traitReducedDataStructure$combinedDFP_Val_Labels$dfP_Val),
                               "; ncol (trait): ", length(levels(as.factor(traitReducedDataStructure$traitClusters))), " (if # is less than defined number of clusters, then duplicate variable names have been omitted)")
      }
      else {
        base::message(base::paste0(sysTimePID(), " is.valid(traitReducedDataStructure()) == FALSE."))
      }
    },
    error = function(e) {
      base::message("An error occurred in updateTxtOmitTraitsOut():\n", e)
    },

    warning = function(w) {
      base::message("A warning occurred in updateTxtOmitTraitsOut():\n", w)
    },
    finally = {
      return(shiny::HTML(result))
    }
  )
}

updateTxtOmitProbesOut <- function(id, probeReducedDataStructure) {
  base::tryCatch(
    {
      result <- NULL
      if (is.valid(probeReducedDataStructure$combinedDFP_Val_Labels$dfP_Val)) {
        result <- base::paste0("reduce probes successful. n(probe): ", length(levels(as.factor(probeReducedDataStructure$probeClusters))), " (if # is less than defined number of clusters, then duplicate probe names have been omitted)")
      }
      else {
        base::message(base::paste0(sysTimePID(), " is.valid(probeReducedDataStructure()) == FALSE."))
      }
    },
    error = function(e) {
      base::message("An error occurred in updateTxtOmitProbesOut():\n", e)
    },

    warning = function(w) {
      base::message("A warning occurred in updateTxtOmitProbesOut():\n", w)
    },
    finally = {
      return(shiny::HTML(result))
    }
  )
}
