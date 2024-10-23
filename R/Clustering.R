Clustering_UI <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    shinyjs::disabled(
      if (id == "P_Val") {
        shiny::actionButton(ns("btnOmitTraits"), label = "Step 4a: Omit Traits (p-val)")
      }
      else {
        shiny::actionButton(ns("btnOmitTraits"), label = "Step 4b: Omit Traits (log(FC))")
      }
    ),
    shiny::verbatimTextOutput(ns("txtOmitOut"), placeholder = TRUE),
    shiny::fluidRow(
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
    ) #end fluidRow
  ) #end tagList
}

Clustering_SERVER <- function(id, pReducedDataStructure, traitReducedDataStructure, session) {
  shiny::moduleServer(id, function(input, output, session) {
    base::tryCatch({
      shiny::observe({
        if (!is.valid(pReducedDataStructure()$combinedDFP_Val_Labels$dfP_Val)) {
          shinyjs::disable("btnOmitTraits")
        }
        else {
          shinyjs::enable("btnOmitTraits")
        }
      })

      shiny::observe({
        if (!is.valid(pReducedDataStructure()$combinedDFP_Val_Labels$dfLogFC)) {
          shinyjs::disable("btnOmitTraits")
        }
        else {
          shinyjs::enable("btnOmitTraits")
        }
      })

      shiny::observeEvent(input$btnOmitTraits,
                          ignoreInit = TRUE,
                          {
                            base::tryCatch(
                              {
                                traitClusterMedoids <- NULL
                                result <- NULL
                                if (is.valid(pReducedDataStructure()$combinedDFP_Val_Labels$dfP_Val)) {
                                  if (is.valid(session$userData$sessionVariables$numClusters())) {
                                    result <- pReducedDataStructure()
                                    base::print(base::paste0(sysTimePID(), " (btnOmitTraitsPval) start generating traitClusterMedoids."))
                                    if (is.valid(result$clustResTraits) && is.valid(result$distMatTraits)) {
                                      traitClusterMedoids <- getTraitClusterMedoids(clustResTraits = result$clustResTraits,
                                                                                    distMatTraits = result$distMatTraits,
                                                                                    numClusters = session$userData$sessionVariables$numClusters())
                                    }
                                    else {
                                      traitClusterMedoids <- NULL
                                    }
                                  }
                                }
                                if (is.valid(pReducedDataStructure()$combinedDFP_Val_Labels$dfP_Val) &&
                                    is.valid(traitClusterMedoids)) {
                                  keys <- session$userData$config$keyAttributes
                                  result <- getTraitReducedData(pReducedDataStructure()$combinedDFP_Val_Labels,
                                                                traitClusterMedoids, keys)
                                }
                                else{
                                  result <- NULL
                                }
                              },
                              error = function(e) {
                                base::message("An error occurred in shiny::observeEvent(input$btnOmitTraitsPval):\n", e)
                              },
                              warning = function(w) {
                                base::message("A warning occurred in shiny::observeEvent(input$btnOmitTraitsPval):\n", w)
                              },
                              finally = {
                                if (id == "P_Val") {
                                  session$userData$sessionVariables$traitReducedDataPVal(result)
                                  #traitReducedDataStructure <- result
                                  #dataStructure(result) #this  wont work: change data itself to a reactive structure
                                }
                                else if (id == "LogFC") {
                                  session$userData$sessionVariables$traitReducedDataLogFC(result)
                                  #traitReducedDataStructure <- result
                                }
                                else {
                                  browser() #should not happen
                                }
                                base::print(base::paste0(sysTimePID(), " finished omitting traits."))
                              }
                            )
                          },
                          ignoreNULL = FALSE
      ) #end observeEvent

      output$txtOmitOut <- shiny::reactive({
        base::tryCatch(
          {
            base::print(base::paste0(sysTimePID(), " start generating output$txtOmitOut."))
            if (is.valid(traitReducedDataStructure()$combinedDFP_Val_Labels)) {
              result <- updateTxtOmitTraitsOut(traitReducedDataStructure()$combinedDFP_Val_Labels)
              if (id == "P_Val") {
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
            base::message("An error occurred in shiny::reactive(output$txtOmitOut):\n", e)
          },
          warning = function(w) {
            base::message("A warning occurred in shiny::reactive(output$txtOmitOut):\n", w)
          },
          finally = {
            base::print(base::paste0(sysTimePID(), " finished generating output$txtOmitOut"))
            return(result)
          }
        )
      }) #end output$txtOmitOut

      distNeigboursProbes10 <- shiny::reactive(calculateDistanceNeigboursProbes(wd = session$userData$packageWd, clustResProbes = traitReducedDataStructure()$clustResProbes, annotation = session$userData$annotation, distanceToLook = 10, numCores = session$userData$numCores))
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

      distNeigboursProbes100 <- shiny::reactive(calculateDistanceNeigboursProbes(wd = session$userData$packageWd, clustResProbes = traitReducedDataStructure()$clustResProbes, annotation = session$userData$annotation, distanceToLook = 100, numCores = session$userData$numCores))
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

      distNeigboursProbes1000 <- shiny::reactive(calculateDistanceNeigboursProbes(wd = session$userData$packageWd, clustResProbes = traitReducedDataStructure()$clustResProbes, annotation = session$userData$annotation, distanceToLook = 1000, numCores = session$userData$numCores))
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

      distNeigboursProbes10000 <- shiny::reactive(calculateDistanceNeigboursProbes(wd = session$userData$packageWd, clustResProbes = traitReducedDataStructure()$clustResProbes, annotation = session$userData$annotation, distanceToLook = 10000, numCores = session$userData$numCores))
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
        if(is.valid(traitReducedDataStructure()$traitDendrogram)) {
          output$plotDendrogramTraitsLong <- shiny::renderPlot(getPlot(traitReducedDataStructure()$traitDendrogram))
          #output$plotDendrogramTraitsLong <- plotly::renderPlotly(getPlot(traitReducedDataStructure()$traitDendrogram))
        }
      })

      shiny::observe({
        if(is.valid(traitReducedDataStructure()$traitClustergram)) {
          output$plotClustergramTraitsLong <- shiny::renderPlot(getPlot(traitReducedDataStructure()$traitClustergram))
          #output$plotClustergramTraitsLong <- plotly::renderPlotly(getPlot(traitReducedDataStructure()$traitClustergram))
        }
      })

      shiny::observe({
        if(is.valid(traitReducedDataStructure()$traitClusterMedoids)) {
          output$DTTraitsMedoids <- DT::renderDataTable(as.data.frame(getMedoidsTable(traitReducedDataStructure()$traitClusterMedoids)))
          #output$DTTraitsMedoids <- DT::renderDataTable(as.data.frame(mtcars))
        }
      })

      shiny::observe({
        output$DTTraitsClusters <- DT::renderDataTable(as.data.frame(getClustersTable(traitReducedDataStructure()$traitClusters,
                                                                                                     traitReducedDataStructure()$traitClusterMedoids)))
      })
    },
    error = function(e) {
      base::message("An error occurred in moduleServer in Clustering_SERVER:\n", e)
      browser() #should not happen
    },
    warning = function(w) {
      base::message("An error occurred in moduleServer in Clustering_SERVER:\n", w)
      browser() #should not happen
    },
    finally = {
      base::print(base::paste0(sysTimePID(), " finished moduleServer in Clustering_SERVER."))
      #return(result)
    })
  } #end moduleServer
)} #end Clustering_SERVER

#' getDistMat
#' calculates distance matrix using parallelDist
#' @param numberCores number of CPU cores to use
#' @param matrix matrix to calculate distance matrix for
#' @return distance matrix
#' examples getDistMat(numberCores, matrix)
getDistMat <- function(numberCores, matrix) {
  id <- shiny::showNotification("Calculating distance matrix parallel...", duration = NULL, closeButton = FALSE)
  base::on.exit(shiny::removeNotification(id), add = TRUE)
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
    browser() #this should not happen, if dist(matrix) delivers a result, then something is wrong with parallelDist
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
  id <- shiny::showNotification("Calculating clustering results...", duration = NULL, closeButton = FALSE)
  base::on.exit(shiny::removeNotification(id), add = TRUE)
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
  id <- shiny::showNotification("Calculating distance neighbours probes...", duration = NULL, closeButton = FALSE)
  base::on.exit(shiny::removeNotification(id), add = TRUE)
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

#' getplotClustergramTraitsLong
#' produces a clustergram for traits
#' @param mat.t transposed matrix with result values
#' @param clustResTraits (reduced) clustering results for traits
#' @param traitClusters trait cluster medoids
#' @return plot with clustergram
#' examples getplotClustergramTraitsLong(mat.t, clustResTraits, traitClusters)
getplotClustergramTraitsLong <- function(mat.t, clustResTraits, traitClusters) {
  id <- shiny::showNotification("Creating trait clustergram long...", duration = NULL, closeButton = FALSE)
  on.exit(shiny::removeNotification(id), add = TRUE)
  base::print(base::paste0(sysTimePID(), " start plotting clustergram."))
  base::tryCatch(
    {
      if (traitClusters > 1) {
        Clusters <- dendextend::cutree(clustResTraits, k = traitClusters) # Clusters <- cutree(clustResTraits, k = traitClusters)
        mycolors <- colorRampPalette(RColorBrewer::brewer.pal(8, "Set2"))(traitClusters)
        base::tryCatch(
          {
            p <- factoextra::fviz_cluster(list(data = mat.t, cluster = Clusters),
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
        # p <- factoextra::fviz_mclust(list(data = matP_Val.t, cluster = Clusters),
        #                               palette = "jco",
        #                               ggtheme = ggplot2::theme_classic())
      }
      else {
        p <- getEmptyPlot()
      }
    },
    error = function(e) {
      base::message("An error occurred in getplotClustergramTraitsLong():\n", e)
    },
    warning = function(w) {
      base::message("A warning occurred in getplotClustergramTraitsLong():\n", w)
    },
    finally = {
      base::print(base::paste0(sysTimePID(), " end plotting clustergram."))
      return(p)
    }
  )
}
