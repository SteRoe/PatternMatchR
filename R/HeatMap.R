HeatMap_UI <- function(id) {
  ns <- shiny::NS(id)
  htmltools::tagList(
    shiny::tabsetPanel(id = ns("tabsetHeatMap"),
      shiny::tabPanel(
       "Table P_VAL",
       if (id == "PVal") {
         "Table of p-value; clustering order comes from clustering of p-values."
       }
       else if (id == "LogFC") {
         "Table of p-value; clustering order comes from clustering of log(FC)."
       }
       else if (id == "PValWOGap") {
         "Table of p-value; clustering order comes from clustering of p-values."
       }
       else if (id == "LogFCWOGap") {
         "Table of p-value; clustering order comes from clustering of log(FC)."
       },
       DT::dataTableOutput(ns("DTP_VAL"))
      ),
      shiny::tabPanel(
       "Table Delta Methylation",
       if (id == "PVal") {
        "Table of delta methylation; clustering order comes from clustering of p-values."
       }
       else if (id == "LogFC") {
         "Table of delta methylation; clustering order comes from clustering of log(FC)."
       }
       else if (id == "PValWOGap") {
         "Table of delta methylation; clustering order comes from clustering of p-values."
       }
       else if (id == "LogFCWOGap") {
         "Table of delta methylation; clustering order comes from clustering of log(FC)."
       },
       DT::dataTableOutput(ns("DTDM"))
      ),
      shiny::tabPanel(
       "Table N",
       if (id == "PVal") {
        "Table of n; clustering order comes from clustering of p-values."
       }
       else if (id == "LogFC") {
         "Table of n; clustering order comes from clustering of log(FC)."
       }
       else if (id == "PValWOGap") {
         "Table of n; clustering order comes from clustering of p-values."
        }
       else if (id == "LogFCWOGap") {
         "Table of n; clustering order comes from clustering of log(FC)."
       },
       DT::dataTableOutput(ns("DTN"))
      ),
      shiny::tabPanel(
        "Table Delta Methylation log(FC)",
        if (id == "PVal") {
          "Table of log fold change(delta methylation); clustering order comes from clustering of p-values."
        }
        else if (id == "LogFC") {
          "Table of log fold change(delta methylation); clustering order comes from clustering of log(FC)."
        }
        else if (id == "PValWOGap") {
          "Table of log fold change(delta methylation); clustering order comes from clustering of p-values."
        }
        else if (id == "LogFCWOGap") {
          "Table of log fold change(delta methylation); clustering order comes from clustering of log(FC)."
        },
        DT::dataTableOutput(ns("DTLogFC"))
      )
    )
  )
}

HeatMap_SERVER <- function(id, session) {
  shiny::moduleServer(id, function(input, output, session) {
    base::tryCatch({
      shiny::observe({
        if (id == "PVal") {
          probeReducedDataStructure <- session$userData$sessionVariables$probeReducedDataStructurePVal()
        }
        else if (id == "LogFC") {
          probeReducedDataStructure <- session$userData$sessionVariables$probeReducedDataStructureLogFC()
        }
        else if (id == "PValWOGap") {
          probeReducedDataStructure <- session$userData$sessionVariables$probeReducedDataStructureWOGapPVal()
        }
        else if (id == "LogFCWOGap") {
          probeReducedDataStructure <- session$userData$sessionVariables$probeReducedDataStructureWOGapLogFC()
        }
        else{
          browser() #should not happen
        }
        output$DTP_VAL <- DT::renderDataTable(as.data.frame(probeReducedDataStructure$combinedDFP_Val_Labels$dfP_Val_w_number))
        output$DTDM <- DT::renderDataTable(as.data.frame(probeReducedDataStructure$combinedDFP_Val_Labels$dfDM_w_number))
        output$DTLogFC <- DT::renderDataTable(as.data.frame(probeReducedDataStructure$combinedDFP_Val_Labels$dfLogFC_w_number))
        output$DTN <- DT::renderDataTable(as.data.frame(probeReducedDataStructure$combinedDFP_Val_Labels$dfN_w_number))
      })
    },
    error = function(e) {
      if (attributes(e)$class[1] != "shiny.silent.error") {
        base::message("An error occurred in moduleServer in HeatMap_SERVER:\n", e)
        browser() #should not happen
      }
    },
    warning = function(w) {
      base::message("An error occurred in moduleServer in HeatMap_SERVER:\n", w)
      browser() #should not happen
    },
    finally = {
      base::print(base::paste0(sysTimePID(), " finished moduleServer in HeatMap_SERVER."))

    })
  }) #end shiny::moduleServer
}

plotCombinedHM <- function(id, input, output, session) {
  shinyId <- shiny::showNotification("Plotting heatmap", duration = NULL, closeButton = FALSE)
  base::on.exit(shiny::removeNotification(shinyId), add = TRUE)
  base::print(base::paste0(sysTimePID(), " start plotting heatmap."))
  output$txtHMDescription_P_Val <-
    shiny::renderText(base::paste0("calculating heatmap..., current plot is not valid"))
  while (!is.null(grDevices::dev.list())) {
    grDevices::dev.off()
  }
  if (id == "PVal") {
    #combinedDFP_Val_Labels <- session$userData$sessionVariables$probeReducedDataStructurePVal()$combinedDFP_Val_Labels
    probeReducedDataStructure <- session$userData$sessionVariables$probeReducedDataStructurePVal()
    traitReducedDataStructure <- session$userData$sessionVariables$traitReducedDataStructurePVal()
#    dfOriginal <- combinedDFP_Val_Labels$dfP_Val_Original #we need the unordered original matrix here, because ComplexHeatmap orders the data by itself according to the order from the dendrogram. That is important!!!!
#  browser() #if step 3 was omitted, we see an error here...
#    dfOriginal[dfOriginal > 0.05] <- NA # 1
  }
  else if (id == "LogFC") {
    #combinedDFP_Val_Labels <- session$userData$sessionVariables$probeReducedDataStructureLogFC()$combinedDFP_Val_Labels
    probeReducedDataStructure <- session$userData$sessionVariables$probeReducedDataStructureLogFC()
    traitReducedDataStructure <- session$userData$sessionVariables$traitReducedDataStructureLogFC()
#    dfOriginal <- combinedDFP_Val_Labels$dfLogFC_Original #we need the unordered original matrix here, because ComplexHeatmap orders the data by itself according to the order from the dendrogram. That is important!!!!
  }
  else if (id == "PValWOGap") {
    probeReducedDataStructure <- session$userData$sessionVariables$probeReducedDataStructureWOGapPVal()
    traitReducedDataStructure <- session$userData$sessionVariables$traitReducedDataStructurePVal()
    #combinedDFP_Val_Labels <- session$userData$sessionVariables$probeReducedDataStructureWOGapPVal()$combinedDFP_Val_Labels
  }
  else if (id == "LogFCWOGap") {
    probeReducedDataStructure <- session$userData$sessionVariables$probeReducedDataStructureWOGapLogFC()
    traitReducedDataStructure <- session$userData$sessionVariables$traitReducedDataStructureLogFC()
    #combinedDFP_Val_Labels <- session$userData$sessionVariables$probeReducedDataStructureWOGapLogFC()$combinedDFP_Val_Labels
  }
  else {
    browser # should not happen
  }
  base::print(base::paste0(sysTimePID(), " calculating combined heatmap."))
  if (TRUE) { # if (nrow(dfOriginal) > 5) {
    startTime <- Sys.time()
    base::tryCatch({
      base::print(base::paste0(sysTimePID(), " gc()"))
      gc()
      base::options(expressions = 500000)
      combinedDFP_Val_Labels <- probeReducedDataStructure$combinedDFP_Val_Labels
      dendProbes <- probeReducedDataStructure$probeDendrogram
      # check clustResProbes > 8
      maxClassesProbes <- 7
      if (maxClassesProbes <= 7) {
        dendProbes <- dendextend::color_branches(dendProbes, maxClassesProbes)
      }
      dendTraits <- traitReducedDataStructure$traitDendrogram
      Distances <- probeReducedDataStructure$DNAdistances
      base::print(base::paste0(sysTimePID(), " before calculating heatmap"))

      selectedRowIndicesYellow <- NULL
      if(is.valid(session$userData$sessionVariables$selectedProbe())) {
        selectedRowIndicesYellow <- session$userData$sessionVariables$selectedProbe()
      }
      # if (is.valid(input$txtSearchFullCpG)) {
      #   selectedRowIndicesYellow <- unlist(strsplit(input$txtSearchFullCpGPVal, split = " ")) #is a list of cg-numbers from search field "txtSearchCpG"
      # }
      # if (is.valid(input$txtSearchFullTrait)) {
      #   selectedColIndices <- unlist(strsplit(input$txtSearchFullTrait, split = " "))
      # }
      #selectedRowIndicesOrange <- session$userData$sessionVariables$distancesBelowThreshold()
      selectedRowIndicesOrange <- NULL
      if (id == "PVal") {
        id2 <- "PVal"
      }
      else if (id == "LogFC") {
        id2 <- "LogFC"
      }
      else if (id == "PValWOGap") {
        id2 <- "PVal"
      }
      else if (id == "LogFCWOGap") {
        id2 <- "LogFC"
      }
      base::print(base::paste0(sysTimePID(), " before l <- combinedDFInteractiveHeatMapP_Val(combinedDFP_Val_Labels, dendProbes, dendTraits, selectedRowIndices, selectedColIndices)"))
      l <-
        combinedDFInteractiveHeatMap(id2, combinedDFP_Val_Labels, dendProbes, dendTraits, Distances, selectedRowIndicesYellow, selectedColIndices, selectedRowIndicesOrange, session)
      base::print(base::paste0(sysTimePID(), " before combinedHM <- l$combinedHM"))
      combinedHM <- l$combinedHM
      endTime <- Sys.time()
      elapsedTime <- endTime - startTime
      base::print(base::paste0(sysTimePID(), " after calculating heatmap. Elapsed time: ", elapsedTime, " sec."))
      base::print(base::paste0(sysTimePID(), " before plotting heatmap; InteractiveComplexHeatmap::makeInteractiveComplexHeatmap()"))
      while (!base::is.null(grDevices::dev.list())) {
        grDevices::dev.off()
      }
      if (id == "PVal") {
        hm_id <- "Heatmap_P_Val"
      }
      else if (id == "LogFC") {
        hm_id <- "Heatmap_LogFC"
      }
      else if (id == "PValWOGap") {
        hm_id <- "Heatmap_P_ValWOGap"
      }
      else if (id == "LogFCWOGap") {
        hm_id <- "Heatmap_LogFCWOGap"
      }
      startTime <- Sys.time()
      InteractiveComplexHeatmap::makeInteractiveComplexHeatmap(
        input = input,
        output = output,
        session = session,
        ht_list = combinedHM,
        heatmap_id = hm_id,
        show_layer_fun = TRUE,
        click_action = click_action_fullHM,
        brush_action = brush_action_fullHM,
        hover_action = hover_action_fullHM
      )
    },
    error = function(e) {
      if (attributes(e)$class[1] != "shiny.silent.error") {
        base::message("An error occurred in plotCombinedHM():\n", e)
        Cstack_info()
        browser() #should not happen
      }
    },
    warning = function(w) {
      base::message("A warning occurred in plotCombinedHM():\n", w)
      browser() #should not happen
    },
    finally = {
      endTime <- Sys.time()
      elapsedTime <- endTime - startTime
      base::print(base::paste0(sysTimePID(), " after plotting heatmap; InteractiveComplexHeatmap::makeInteractiveComplexHeatmap(). Elapsed time: ", elapsedTime, " sec."))
      HMDescription <-
        shiny::renderText(base::paste0(sysTimePID(), " done plotting heatmap..., current plot is valid; elapsed time: ", elapsedTime, " sec."))
      if (id == "PVal") {
        output$txtHMDescription_P_Val <- HMDescription
      }
      else if (id == "LogFC") {
        output$txtHMDescription_LogFC <- HMDescription
      }
      else if (id == "PValWOGap") {
        output$txtHMDescription_P_ValWOGap <- HMDescription
      }
      else if (id == "LogFCWOGap") {
        output$txtHMDescription_LogFCWOGap <- HMDescription
      }
    })
  }
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

# emptyHM <- compiler::cmpfun(emptyHM)

#' creates a regular heatmap
#' @param combinedDF_Labels list of data.frame and labels generated from function mergeDFP_Val_Labels()
#' @param dendProbes dendrogram (without labels) for probes (rows), generated externally and providing information for clustering of heatmap
#' @param dendTraits dendrogram (without labels) result for traits (columns), generated externally and providing information for clustering of heatmap
#' @param selectedRowIndicesYellow indicies of HM rows to mark in yellow color
#' @param selectedColIndices indicies of HM cols to mark (not used so far)
#' @param selectedRowIndicesOrange indicies of HM rows to mark in orange color
#' @param session session object for reference
#' @return heatmap object for InteractiveComplexHeatmap::makeInteractiveComplexHeatmap
# examples combinedDFInteractiveHeatMapP_Val(combinedDF_Labels, dendProbes, dendTraits, Distances, selectedRowIndicesYellow, selectedColIndices, selectedRowIndicesOrange)
combinedDFInteractiveHeatMap <-
  function(id,
           combinedDF_Labels,
           dendProbes = NA,
           dendTraits = NA,
           Distances = NA,
           selectedRowIndicesYellow = NA,
           selectedColIndices = NA,
           selectedRowIndicesOrange = NA,
           session = session) {
    #function is not called twice (check caller function), if so, check here
    shinyId <- shiny::showNotification("Creating heatmap...", duration = NULL, closeButton = FALSE)
    base::on.exit(shiny::removeNotification(shinyId), add = TRUE)
    base::tryCatch(
      {
        startTime <- Sys.time()
        base::print(base::paste0(sysTimePID(), " start preparing HM; combinedDFInteractiveHeatMap()"))
        matP_Val <- base::as.matrix(combinedDF_Labels$dfP_Val_Original)
        matDM <- base::as.matrix(combinedDF_Labels$dfDM_Original)
        matN <- base::as.matrix(combinedDF_Labels$dfN_Original)
        matLogFC <- base::as.matrix(combinedDF_Labels$dfLogFC_Original)
        if (id == "PVal") {
          matHM <- matP_Val
        }
        else if (id == "LogFC") {
          matHM <- matLogFC
        }
        # use rasterization like described in
        # https://jokergoo.github.io/2020/06/30/rasterization-in-complexheatmap/
        base::print(base::paste0(sysTimePID(), " making labels"))
        labelsDF1 <- combinedDF_Labels$mergedOriginalColnames[which(combinedDF_Labels$mergedOriginTrait == 1)] # labelsDF1 <- combinedDF_Labels$labelsDF1
        labelsDF2 <- combinedDF_Labels$mergedOriginalColnames[which(combinedDF_Labels$mergedOriginTrait == 2)] # labelsDF2 <- combinedDF_Labels$labelsDF2
        labelsDF3 <- combinedDF_Labels$mergedOriginalColnames[which(combinedDF_Labels$mergedOriginTrait == 3)] # labelsDF3 <- combinedDF_Labels$labelsDF3
        l1 <- base::rep("trait 1", base::length(labelsDF1))
        l2 <- base::rep("trait 2", base::length(labelsDF2))
        l3 <- base::rep("trait 3", base::length(labelsDF3))
        labels <- base::c(l1, l2, l3)
        while (!base::is.null(grDevices::dev.list())) {
          grDevices::dev.off()
        }
        base::print(base::paste0(sysTimePID(), " making annotations"))
        ha <- ComplexHeatmap::columnAnnotation(
          classes = labels,
          col = base::list(
            classes = base::c(
              "trait 1" = "red",
              "trait 2" = "green",
              "trait 3" = "blue"
            )
          )
        )

        base::print(base::paste0(sysTimePID(), " making colors"))
        max.Col1 <- 0.05
        min.Col1 <- base::min(matHM, na.rm = TRUE)
        max.Col2 <- 0.05
        min.Col2 <- base::min(matHM, na.rm = TRUE)
        max.Col3 <- 0.05
        min.Col3 <- base::min(matHM, na.rm = TRUE)
        col1 <-
          circlize::colorRamp2(
            base::seq(max.Col1, min.Col1, length = 9),
            RColorBrewer::brewer.pal(9, "OrRd")
          )
        col2 <-
          circlize::colorRamp2(
            base::seq(max.Col2, min.Col2, length = 9),
            RColorBrewer::brewer.pal(9, "YlGn")
          )
        col3 <-
          circlize::colorRamp2(
            base::seq(max.Col3, min.Col3, length = 9),
            RColorBrewer::brewer.pal(9, "GnBu")
          )
        base::print(base::paste0(sysTimePID(), " preparing heatmap n(probes)=", base::dim(matHM)[1], " x n(traits)=", base::dim(matHM)[2]))
        if ("dendrogram" %in% class(dendProbes) && "dendrogram" %in% class(dendTraits)) {
          base::print(base::paste0(sysTimePID(), " length(unlist(dendProbes)): ", length(unlist(dendProbes))))
          base::print(base::paste0(sysTimePID(), " length(unlist(dendTraits)): ", length(unlist(dendTraits))))
          ht <-
            ComplexHeatmap::Heatmap(
              matHM,
              rect_gp = grid::gpar(type = "none"),
              cluster_rows = dendProbes,
              cluster_columns = dendTraits,
              top_annotation = ha,
              layer_fun = function(j, i, x, y, w, h, fill) {
                if (length(i) == base::nrow(matHM) * base::ncol(matHM)) {
                  # we are in main HM
                  subHM <- FALSE
                } else {
                  # we are in sub HM
                  subHM <- TRUE
                }
                l <- labels[j] == "trait 1"
                if (any(l)) {
                  grid::grid.rect(x[l], y[l], w[l], h[l], gp = grid::gpar(
                    fill = col1(ComplexHeatmap::pindex(matHM, i[l], j[l])),
                    col = col1(ComplexHeatmap::pindex(matHM, i[l], j[l]))
                  ))
                }
                l <- labels[j] == "trait 2"
                if (any(l)) {
                  grid::grid.rect(x[l], y[l], w[l], h[l], gp = grid::gpar(
                    fill = col2(ComplexHeatmap::pindex(matHM, i[l], j[l])),
                    col = col2(ComplexHeatmap::pindex(matHM, i[l], j[l]))
                  ))
                }
                l <- labels[j] == "trait 3"
                if (any(l)) {
                  grid::grid.rect(x[l], y[l], w[l], h[l], gp = grid::gpar(
                    fill = col3(ComplexHeatmap::pindex(matHM, i[l], j[l])),
                    col = col3(ComplexHeatmap::pindex(matHM, i[l], j[l]))
                  ))
                }
                #mark selected row indices
                labelsProbes  <- labels(dendProbes)
                if (is.valid(selectedRowIndicesYellow)) {
                  rowsToMarkYellow <- labelsProbes %in% selectedRowIndicesYellow
                  if (any(rowsToMarkYellow)) {
                    grid::grid.rect(x, y[rowsToMarkYellow], w, h[rowsToMarkYellow], gp = grid::gpar(
                      fill = "yellow",
                      col = "yellow"
                    ))
                  }
                }
                if (is.valid(selectedRowIndicesOrange)) {
                  rowsToMarkOrange <- labelsProbes %in% selectedRowIndicesOrange
                  if (any(rowsToMarkOrange)) {
                    grid::grid.rect(x, y[rowsToMarkOrange], w, h[rowsToMarkOrange], gp = grid::gpar(
                      fill = "orange",
                      col = "orange"
                    ))
                  }
                }
                if (subHM == TRUE) {
                  grid::grid.text(
                    paste0(
                      "p:",
                      sprintf("%.G", ComplexHeatmap::pindex(matP_Val, i, j)),
                      "\n",
                      "d:",
                      sprintf("%.G", ComplexHeatmap::pindex(matDM, i, j)),
                      "\n",
                      "n:",
                      ComplexHeatmap::pindex(matN, i, j),
                      "log(fc):",
                      sprintf("%.G", ComplexHeatmap::pindex(matLogFC, i, j))
                    ),
                    x,
                    y,
                    gp = grid::gpar(fontsize = 8)
                  )
                }
              },
              show_heatmap_legend = FALSE,
              use_raster = TRUE,
              raster_by_magick = TRUE
            )
        } else {
          base::print(base::paste0(sysTimePID(), " at least one distance matrix is not of class \"dendrogram\""))
        }
        base::print(base::paste0(sysTimePID(), " making legend"))
        lgd <- list(
          ComplexHeatmap::Legend(
            title = "trait 1",
            col_fun = col1,
            at = c(max.Col1, min.Col1),
            labels = c(max.Col1, paste0(
              extractMantissaExponent(min.Col1)$exponent
            ))
          ),
          ComplexHeatmap::Legend(
            title = "trait 2",
            col_fun = col2,
            at = c(max.Col2, min.Col2),
            labels = c(max.Col2, paste0(
              extractMantissaExponent(min.Col2)$exponent
            ))
          ),
          ComplexHeatmap::Legend(
            title = "trait 3",
            col_fun = col3,
            at = c(max.Col3, min.Col3),
            labels = c(max.Col3, paste0(
              extractMantissaExponent(min.Col3)$exponent
            ))
          )
        )
      },
      error = function(e) {
        if (attributes(e)$class[1] != "shiny.silent.error") {
          base::message("An error occurred in combinedDFInteractiveHeatMapP_Val():\n", e)
        }
      },
      warning = function(w) {
        base::message("A warning occurred in combinedDFInteractiveHeatMapP_Val():\n", w)
      },
      finally = {
        endTime <- Sys.time()
        elapsedTime <- endTime - startTime
        base::print(base::paste0(sysTimePID(), " end preparing heatmap. Elapsed time: ", elapsedTime, "."))
      }
    )
    base::tryCatch(
      {
        matDistances <- cbind(Distances$minDistance, Distances$meanDistance)
        htDistances <-
          ComplexHeatmap::Heatmap(
            matDistances,
            cluster_rows = dendProbes
          )
        minDist <- base::min(matDistances, na.rm = TRUE)
        maxDist <- base::max(matDistances, na.rm = TRUE)
        colDist <-
          circlize::colorRamp2(
            base::seq(minDist, maxDist, length = 9),
            RColorBrewer::brewer.pal(9, "OrRd")
          )
        lgdHMDist <- list(
          ComplexHeatmap::Legend(
            title = "distances",
            col_fun = colDist,
            at = c(maxDist, minDist),
            labels = c(maxDist, minDist)
          )
        )
      },
      error = function(e) {
        if (attributes(e)$class[1] != "shiny.silent.error") {
          base::message("An error occurred in making htDistances():\n", e)
        }
      },
      warning = function(w) {
        base::message("A warning occurred in making htDistances():\n", w)
      },
      finally = {
        base::print(base::paste0(sysTimePID(), " end making htDistances()."))
      }
    )
    base::tryCatch(
      {
        base::print(base::paste0(sysTimePID(), " clearing grDevices()"))
        if (grDevices::dev.cur() > 1) {
          grDevices::dev.off()
        }
        grDevices::pdf(NULL)
      },
      error = function(e) {
        if (attributes(e)$class[1] != "shiny.silent.error") {
          base::message("An error occurred in clearing grDevices():\n", e)
        }
      },
      warning = function(w) {
        base::message("A warning occurred in clearing grDevices():\n", w)
      },
      finally = {
        base::print(base::paste0(sysTimePID(), " end clearing grDevices()."))
      }
    )
    base::tryCatch(
      {
        base::print(base::paste0(sysTimePID(), " start drawing heatmap (takes some time). (step before ComplexHeatmap::draw())"))
        # with huge heatmaps, the following error occurs:
        # Error in Cairo: Failed to create Cairo backend!
        ht <- ComplexHeatmap::draw(htDistances + ht, annotation_legend_list = lgd)
      },
      error = function(e) {
        if (attributes(e)$class[1] != "shiny.silent.error") {
          base::message(base::paste0(sysTimePID(), " Error: unable to draw HM. ", e$message))
        }
      },
      warning = function(w) {
        base::message(base::paste0(sysTimePID(), " unable to draw HM. ", w$message))
      },
      finally = {
        base::print(base::paste0(sysTimePID(), " grDevices::dev.off()"))
        grDevices::dev.off()
        base::print(base::paste0(sysTimePID(), " l <- base::list()"))
        l <- base::list()
        base::print(base::paste0(sysTimePID(), " l$combinedHM <- ht"))
        l$combinedHM <- ht
        endTime <- Sys.time()
        elapsedTime <- endTime - startTime
        base::print(base::paste0(sysTimePID(), " end combinedDFInteractiveHeatMap()"))
        base::print(base::paste0(sysTimePID(), " finished drawing heatmap (takes some time). (step after ComplexHeatmap::draw(); Elapsed time: ", elapsedTime, "."))
        base::return(l)
      }
    )
  }

#' creates a regular heatmap based on p-values
#' @param combinedDF_Labels list of data.frame and labels generated from function mergeDFP_Val_Labels()
#' @param dendProbes dendrogram (without labels) for probes (rows), generated externally and providing information for clustering of heatmap
#' @param dendTraits dendrogram (without labels) result for traits (columns), generated externally and providing information for clustering of heatmap
#' @param selectedRowIndicesYellow indicies of HM rows to mark in yellow color
#' @param selectedColIndices indicies of HM cols to mark (not used so far)
#' @param selectedRowIndicesOrange indicies of HM rows to mark in orange color
#' @param session session object for reference
#' @return heatmap object for InteractiveComplexHeatmap::makeInteractiveComplexHeatmap
# examples combinedDFInteractiveHeatMapP_Val(combinedDF_Labels, dendProbes, dendTraits, Distances, selectedRowIndicesYellow, selectedColIndices, selectedRowIndicesOrange)
combinedDFInteractiveHeatMapLogFC <-
  function(combinedDF_Labels,
           dendProbes = NA,
           dendTraits = NA,
           Distances = NA,
           selectedRowIndicesYellow = NA,
           selectedColIndices = NA,
           selectedRowIndicesOrange = NA,
           session = session) {
    #function is not called twice (check caller function), if so, check here
    shinyId <- shiny::showNotification("Creating heatmap for (log(FC))...", duration = NULL, closeButton = FALSE)
    base::on.exit(shiny::removeNotification(shinyId), add = TRUE)
    base::tryCatch(
      {
        #check length(unlist(dendProbes))
        startTime <- Sys.time()
        base::print(base::paste0(sysTimePID(), " start preparing HM; combinedDFInteractiveHeatMapLogFC()"))
        matP_Val <- base::as.matrix(combinedDF_Labels$dfP_Val)
        matDM <- base::as.matrix(combinedDF_Labels$dfDM)
        matN <- base::as.matrix(combinedDF_Labels$dfN)
        matLogFC <- base::as.matrix(combinedDF_Labels$dfLogFC)
        # use rasterization like described in
        # https://jokergoo.github.io/2020/06/30/rasterization-in-complexheatmap/
        base::print(base::paste0(sysTimePID(), " making labels"))
        labelsDF1 <- combinedDF_Labels$mergedOriginalColnames[which(combinedDF_Labels$mergedOriginTrait == 1)] # labelsDF1 <- combinedDF_Labels$labelsDF1
        labelsDF2 <- combinedDF_Labels$mergedOriginalColnames[which(combinedDF_Labels$mergedOriginTrait == 2)] # labelsDF2 <- combinedDF_Labels$labelsDF2
        labelsDF3 <- combinedDF_Labels$mergedOriginalColnames[which(combinedDF_Labels$mergedOriginTrait == 3)] # labelsDF3 <- combinedDF_Labels$labelsDF3
        l1 <- base::rep("trait 1", base::length(labelsDF1))
        l2 <- base::rep("trait 2", base::length(labelsDF2))
        l3 <- base::rep("trait 3", base::length(labelsDF3))
        labels <- base::c(l1, l2, l3)
        while (!base::is.null(grDevices::dev.list())) {
          grDevices::dev.off()
        }
        base::print(base::paste0(sysTimePID(), " making annotations"))
        ha <- ComplexHeatmap::columnAnnotation(
          classes = labels,
          col = base::list(
            classes = base::c(
              "trait 1" = "red",
              "trait 2" = "green",
              "trait 3" = "blue"
            )
          )
        )

        base::print(base::paste0(sysTimePID(), " making colors"))
        max.Col1 <- 0.05
        min.Col1 <- base::min(matLogFC, na.rm = TRUE)
        max.Col2 <- 0.05
        min.Col2 <- base::min(matLogFC, na.rm = TRUE)
        max.Col3 <- 0.05
        min.Col3 <- base::min(matLogFC, na.rm = TRUE)
        col1 <-
          circlize::colorRamp2(
            base::seq(max.Col1, min.Col1, length = 9),
            RColorBrewer::brewer.pal(9, "OrRd")
          )
        col2 <-
          circlize::colorRamp2(
            base::seq(max.Col2, min.Col2, length = 9),
            RColorBrewer::brewer.pal(9, "YlGn")
          )
        col3 <-
          circlize::colorRamp2(
            base::seq(max.Col3, min.Col3, length = 9),
            RColorBrewer::brewer.pal(9, "GnBu")
          )
        base::print(base::paste0(sysTimePID(), " preparing heatmap n(probes)=", base::dim(matLogFC)[1], " x n(traits)=", base::dim(matLogFC)[2]))
        if ("dendrogram" %in% class(dendProbes) && "dendrogram" %in% class(dendTraits)) {
          base::print(base::paste0(sysTimePID(), " length(unlist(dendProbes)): ", length(unlist(dendProbes))))
          base::print(base::paste0(sysTimePID(), " length(unlist(dendTraits)): ", length(unlist(dendTraits))))
          length(unlist(dendTraits)) == base::dim(matLogFC)[2]
          length(unlist(dendProbes)) == base::dim(matLogFC)[1]
          ht <-
            ComplexHeatmap::Heatmap(
              matLogFC,
              rect_gp = grid::gpar(type = "none"),
              cluster_rows = dendProbes,
              cluster_columns = dendTraits,
              top_annotation = ha,
              layer_fun = function(j, i, x, y, w, h, fill) {
                if (length(i) == base::nrow(matLogFC) * base::ncol(matLogFC)) {
                  # we are in main HM
                  subHM <- FALSE
                } else {
                  # we are in sub HM
                  subHM <- TRUE
                }
                l <- labels[j] == "trait 1"
                if (any(l)) {
                  grid::grid.rect(x[l], y[l], w[l], h[l], gp = grid::gpar(
                    fill = col1(ComplexHeatmap::pindex(matLogFC, i[l], j[l])),
                    col = col1(ComplexHeatmap::pindex(matLogFC, i[l], j[l]))
                  ))
                }
                l <- labels[j] == "trait 2"
                if (any(l)) {
                  grid::grid.rect(x[l], y[l], w[l], h[l], gp = grid::gpar(
                    fill = col2(ComplexHeatmap::pindex(matLogFC, i[l], j[l])),
                    col = col2(ComplexHeatmap::pindex(matLogFC, i[l], j[l]))
                  ))
                }
                l <- labels[j] == "trait 3"
                if (any(l)) {
                  grid::grid.rect(x[l], y[l], w[l], h[l], gp = grid::gpar(
                    fill = col3(ComplexHeatmap::pindex(matLogFC, i[l], j[l])),
                    col = col3(ComplexHeatmap::pindex(matLogFC, i[l], j[l]))
                  ))
                }
                #mark selected row indices
                labelsProbes  <- labels(dendProbes)
                if (is.valid(selectedRowIndicesYellow)) {
                  rowsToMarkYellow <- labelsProbes %in% selectedRowIndicesYellow
                  if (any(rowsToMarkYellow)) {
                    grid::grid.rect(x, y[rowsToMarkYellow], w, h[rowsToMarkYellow], gp = grid::gpar(
                      fill = "yellow",
                      col = "yellow"
                    ))
                  }
                }
                if (is.valid(selectedRowIndicesOrange)) {
                  rowsToMarkOrange <- labelsProbes %in% selectedRowIndicesOrange
                  if (any(rowsToMarkOrange)) {
                    grid::grid.rect(x, y[rowsToMarkOrange], w, h[rowsToMarkOrange], gp = grid::gpar(
                      fill = "orange",
                      col = "orange"
                    ))
                  }
                }
                if (subHM == TRUE) {
                  grid::grid.text(
                    paste0(
                      "p:",
                      sprintf("%.G", ComplexHeatmap::pindex(matP_Val, i, j)),
                      "\n",
                      "d:",
                      sprintf("%.G", ComplexHeatmap::pindex(matDM, i, j)),
                      "\n",
                      "n:",
                      ComplexHeatmap::pindex(matN, i, j),
                      "\n",
                      "n:",
                      sprintf("%.G", ComplexHeatmap::pindex(matLogFC, i, j))
                    ),
                    x,
                    y,
                    gp = grid::gpar(fontsize = 8)
                  )
                }
              },
              show_heatmap_legend = FALSE,
              use_raster = TRUE,
              raster_by_magick = TRUE
            )
        } else {
          base::print(base::paste0(sysTimePID(), " at least one distance matrix is not of class \"dendrogram\""))
        }
        base::print(base::paste0(sysTimePID(), " making legend"))
        lgd <- list(
          ComplexHeatmap::Legend(
            title = "trait 1",
            col_fun = col1,
            at = c(max.Col1, min.Col1),
            labels = c(max.Col1, paste0(
              extractMantissaExponent(min.Col1)$exponent
            ))
          ),
          ComplexHeatmap::Legend(
            title = "trait 2",
            col_fun = col2,
            at = c(max.Col2, min.Col2),
            labels = c(max.Col2, paste0(
              extractMantissaExponent(min.Col2)$exponent
            ))
          ),
          ComplexHeatmap::Legend(
            title = "trait 3",
            col_fun = col3,
            at = c(max.Col3, min.Col3),
            labels = c(max.Col3, paste0(
              extractMantissaExponent(min.Col3)$exponent
            ))
          )
        )
      },
      error = function(e) {
        if (attributes(e)$class[1] != "shiny.silent.error") {
          base::message("An error occurred in combinedDFInteractiveHeatMapLogFC():\n", e)
        }
      },
      warning = function(w) {
        base::message("A warning occurred in combinedDFInteractiveHeatMapLogFC():\n", w)
      },
      finally = {
        endTime <- Sys.time()
        elapsedTime <- endTime - startTime
        base::print(base::paste0(sysTimePID(), " end preparing heatmap. Elapsed time: ", elapsedTime, "."))
      }
    )
    base::tryCatch(
      {
        matDistances <- cbind(Distances$minDistance, Distances$meanDistance)
        htDistances <-
          ComplexHeatmap::Heatmap(
            matDistances,
            cluster_rows = dendProbes
          )
        minDist <- base::min(matDistances, na.rm = TRUE)
        maxDist <- base::max(matDistances, na.rm = TRUE)
        colDist <-
          circlize::colorRamp2(
            base::seq(minDist, maxDist, length = 9),
            RColorBrewer::brewer.pal(9, "OrRd")
          )
        lgdHMDist <- list(
          ComplexHeatmap::Legend(
            title = "distances",
            col_fun = colDist,
            at = c(maxDist, minDist),
            labels = c(maxDist, minDist)
          )
        )
      },
      error = function(e) {
        if (attributes(e)$class[1] != "shiny.silent.error") {
          base::message("An error occurred in making htDistances():\n", e)
        }
      },
      warning = function(w) {
        base::message("A warning occurred in making htDistances():\n", w)
      },
      finally = {
        base::print(base::paste0(sysTimePID(), " end making htDistances()."))
      }
    )
    base::tryCatch(
      {
        base::print(base::paste0(sysTimePID(), " clearing grDevices()"))
        if (grDevices::dev.cur() > 1) {
          grDevices::dev.off()
        }
        grDevices::pdf(NULL)
      },
      error = function(e) {
        if (attributes(e)$class[1] != "shiny.silent.error") {
          base::message("An error occurred in clearing grDevices():\n", e)
        }
      },
      warning = function(w) {
        base::message("A warning occurred in clearing grDevices():\n", w)
      },
      finally = {
        base::print(base::paste0(sysTimePID(), " end clearing grDevices()."))
      }
    )
    base::tryCatch(
      {
        base::print(base::paste0(sysTimePID(), " start drawing heatmap (takes some time). (step before ComplexHeatmap::draw())"))
        # with huge heatmaps, the following error occurs:
        # Error in Cairo: Failed to create Cairo backend!
        # ht <- ComplexHeatmap::draw(ht + ht2 + ht3, annotation_legend_list = lgd)
        # ht <- ComplexHeatmap::draw(ht, annotation_legend_list = lgd)
        ht <- ComplexHeatmap::draw(htDistances + ht, annotation_legend_list = lgd)
      },
      error = function(e) {
        if (attributes(e)$class[1] != "shiny.silent.error") {
          base::message(base::paste0(sysTimePID(), " Error: unable to draw HM. ", e$message))
        }
      },
      warning = function(w) {
        base::message(base::paste0(sysTimePID(), " unable to draw HM. ", w$message))
      },
      finally = {
        base::print(base::paste0(sysTimePID(), " grDevices::dev.off()"))
        grDevices::dev.off()
        base::print(base::paste0(sysTimePID(), " l <- base::list()"))
        l <- base::list()
        base::print(base::paste0(sysTimePID(), " l$combinedHM <- ht"))
        l$combinedHM <- ht
        endTime <- Sys.time()
        elapsedTime <- endTime - startTime
        base::print(base::paste0(sysTimePID(), " end combinedDFInteractiveHeatMapLogFC()"))
        base::print(base::paste0(sysTimePID(), " finished drawing heatmap (takes some time). (step after ComplexHeatmap::draw(); Elapsed time: ", elapsedTime, "."))
        base::return(l)
      }
    )
  }

HeatMapDistances <-
  function(Distances, dendProbes = NA, session = session) {
    shinyId <- shiny::showNotification("Creating heatmap for probe distances...", duration = NULL, closeButton = FALSE)
    base::on.exit(shiny::removeNotification(shinyId), add = TRUE)
    base::tryCatch(
      {
        startTime <- Sys.time()
        base::print(base::paste0(sysTimePID(), " start preparing HM; HeatMapDistances()"))
        #matDistances <- cbind(Distances$minDistance, Distances$meanDistance, Distances$maxDistance, Distances$number)
        matDistances <- cbind(Distances$minDistance, Distances$meanDistance)
        # N <- length(Distances$minDistance)
        # M <- 5
        # matDistances <- matrix( rnorm(N*M,mean=0,sd=1), N, M)
        while (!base::is.null(grDevices::dev.list())) {
          grDevices::dev.off()
        }
        ht <-
          ComplexHeatmap::Heatmap(
            matDistances,
            #rect_gp = grid::gpar(type = "none"),
            cluster_rows = dendProbes
          )
      },
      error = function(e) {
        if (attributes(e)$class[1] != "shiny.silent.error") {
          base::message("An error occurred in HeatMapDistances():\n", e)
        }
      },
      warning = function(w) {
        base::message("A warning occurred in HeatMapDistances():\n", w)
      },
      finally = {
        endTime <- Sys.time()
        elapsedTime <- endTime - startTime
        base::print(base::paste0(sysTimePID(), " end preparing heatmap. Elapsed time: ", elapsedTime, "."))
      }
    )

    base::tryCatch(
      {
        base::print(base::paste0(sysTimePID(), " clearing grDevices()"))
        if (grDevices::dev.cur() > 1) {
          grDevices::dev.off()
        }
        grDevices::pdf(NULL)
      },
      error = function(e) {
        if (attributes(e)$class[1] != "shiny.silent.error") {
          base::message("An error occurred in clearing grDevices():\n", e)
        }
      },
      warning = function(w) {
        base::message("A warning occurred in clearing grDevices():\n", w)
      },
      finally = {
        base::print(base::paste0(sysTimePID(), " end clearing grDevices()."))
      }
    )
    base::tryCatch(
      {
        base::print(base::paste0(sysTimePID(), " start drawing heatmap (takes some time). (step before ComplexHeatmap::draw())"))
        # with huge heatmaps, the following error occurs:
        # Error in Cairo: Failed to create Cairo backend!
        ht <- ComplexHeatmap::draw(ht)
      },
      error = function(e) {
        if (attributes(e)$class[1] != "shiny.silent.error") {
          base::message(base::paste0(sysTimePID(), " Error: unable to draw HM. ", e$message))
        }
      },
      warning = function(w) {
        base::message(base::paste0(sysTimePID(), " unable to draw HM. ", w$message))
      },
      finally = {
        grDevices::dev.off()
        result <- base::list()
        base::print(base::paste0(sysTimePID(), " result$HeatMapDistances <- ht"))
        result$HeatMapDistances <- ht
        endTime <- Sys.time()
        elapsedTime <- endTime - startTime
        base::print(base::paste0(sysTimePID(), " end HeatMapDistances()"))
        base::print(base::paste0(sysTimePID(), " finished drawing heatmap (takes some time). (step after HeatMapDistances::draw(); Elapsed time: ", elapsedTime, "."))
        base::return(result)
      }
    )
  }

getSearchResultCpGHMPositions <- function(txtSearchCpG, dataStructure) {
  base::tryCatch(
    {
      base::print(base::paste0(sysTimePID(), " searching CpG"))
      #look into clustResProbes and find position of CpG
      CpG <- dataStructure$clustResProbes$labels[dataStructure$clustResProbes$order]
      positions <- base::which(CpG %in% unlist(base::strsplit(base::trimws(txtSearchCpG), " ")))
    },
    error = function(e) {
      if (attributes(e)$class[1] != "shiny.silent.error") {
        base::message("An error occurred in getSearchResultCpG():\n", e)
      }
    },
    warning = function(w) {
      base::message("A warning occurred in getSearchResultCpG():\n", w)
    },
    finally = {
      base::print(base::paste0(sysTimePID(), " end getSearchResultCpG()."))
      base::return(positions)
    }
  )
}

getSearchResultTraitHMPositions <- function(txtSearchTrait, dataStructure) {
  base::tryCatch(
    {
      base::print(base::paste0(sysTimePID(), " searching Trait"))
      #look into clustResTraits and find position of Trait
      Trait <- dataStructure$clustResTraits$labels[dataStructure$clustResTraits$order]
      positions <- base::which(Trait %in% unlist(base::strsplit(base::trimws(txtSearchTrait), " ")))
    },
    error = function(e) {
      if (attributes(e)$class[1] != "shiny.silent.error") {
        base::message("An error occurred in getSearchResultTrait():\n", e)
      }
    },
    warning = function(w) {
      base::message("A warning occurred in getSearchResultTrait():\n", w)
    },
    finally = {
      base::print(base::paste0(sysTimePID(), " end getSearchResultTrait()."))
      base::return(positions)
    }
  )
}

# plotHMDNADistances <- function(input, output, session) {
#   DNAdistances <- session$userData$sessionVariables$traitReducedDataStructurePVal()$DNAdistances
#   dendProbes <- session$userData$sessionVariables$traitReducedDataStructurePVal()$probeDendrogram
#   dendProbes <-
#     dendextend::color_branches(dendProbes, as.integer(input$txtMaxClassesProbes))
#   dendTraits <- session$userData$sessionVariables$traitReducedDataStructurePVal()$traitDendrogram
#
#   base::print(base::paste0(sysTimePID(), " before calculating heatmap"))
#
#   base::print(base::paste0(sysTimePID(), " length(unlist(dendProbes)): ", length(unlist(dendProbes))))
#   base::print(base::paste0(sysTimePID(), " length(unlist(dendTraits)): ", length(unlist(dendTraits))))
#
#   HMDistances <- HeatMapDistances(DNAdistances, dendProbes = dendProbes, session = session)
#   HMDistances <- HMDistances$HeatMapDistances
#   InteractiveComplexHeatmap::makeInteractiveComplexHeatmap(
#     input = input,
#     output = output,
#     session = session,
#     ht_list = HMDistances,
#     heatmap_id = "Heatmap_DNADistances",
#     show_layer_fun = FALSE
#   )
# }

#' click_action_fullHM
#' @param df data.frame containing data which is enclosed by brush action
#' @param input shiny input object
#' @param output shiny output object
#' @param session shiny session object
#' @return nothing; function info label for click action in HM
# examples click_action_fullHM_P_Val(df, input, output, session)
click_action_fullHM <- function(df, input, output, session) {
  base::tryCatch(
    {
      base::print(base::paste0(sysTimePID(), " start click_action_fullHM()."))
      # output[["info_HM_P_Val"]] <- shiny::renderUI({
      #   if (!is.null(df)) {
      #     htmltools::HTML(
      #       GetoptLong::qq(
      #         "<p style='background-color:#FF8080;color:white;padding:5px;'>
      #         row_label: @{df$row_label}, col_label: @{df$column_label},
      #         row: @{df$row_index}, column: @{df$column_index}</p>"
      #       )
      #     )
      #   }
      # })
    },
    error = function(e) {
      if (attributes(e)$class[1] != "shiny.silent.error") {
        base::message("An error occurred in click_action_fullHM():\n", e)
      }
    },
    warning = function(w) {
      base::message("A warning occurred in click_action_fullHM():\n", w)
    },
    finally = {
      base::print(base::paste0(sysTimePID(), " end click_action_fullHM()."))
    }
  )
}

#' brush_action_fullHM
#' @param df data.frame containing data which is enclosed by brush action
#' @param input shiny input object
#' @param output shiny output object
#' @param session shiny session object
#' @return nothing; function creates SPLOM
# examples brush_action_fullHM(df, input, output, session)
brush_action_fullHM <- function(df, input, output, session) {
  base::tryCatch(
    if (!is.null(df)) {
      row_index <- collapse::funique(unlist(df$row_index))
      column_index <- collapse::funique(unlist(df$column_index))
      #feed in selected CpG here

      selectedProbe <- rownames(session$userData$sessionVariables$probeReducedDataStructurePVal()$combinedDFP_Val_Labels$dfP_Val)[row_index]
      session$userData$sessionVariables$selectedProbe(selectedProbe)

      originTrait <- session$userData$sessionVariables$probeReducedDataStructurePVal()$combinedDFP_Val_Labels$mergedOriginTrait[column_index]
      traitLabels <- session$userData$sessionVariables$probeReducedDataStructurePVal()$combinedDFP_Val_Labels$mergedOriginalColnames[column_index]
      selectedTrait <- cbind(traitLabels, originTrait)
      colnames(selectedTrait) <- c("traitName", "traitSource")
      session$userData$sessionVariables$selectedTrait(selectedTrait)
      selectedTraitID <- session$userData$sessionVariables$probeReducedDataStructurePVal()$combinedDFP_Val_Labels$traitID[column_index]
      session$userData$sessionVariables$selectedTraitID(selectedTraitID)
      base::print(base::paste0(sysTimePID(), " last selection was from full heatmap brush event."))
      output$txtGlobalSelectOut <- shiny::renderText("last selection was from full heatmap brush event.")
    },
    error = function(e) {
      if (attributes(e)$class[1] != "shiny.silent.error") {
        base::message("An error occurred in brush_action_fullHM():\n", e)
      }
    },
    warning = function(w) {
      base::message("A warning occurred in brush_action_fullHM():\n", w)
    },
    finally = {
      base::print(base::paste0(sysTimePID(), " end brush_action_fullHM()."))
    }
  )
}

#' hover_action_fullHM
#' @param df data.frame containing data which is enclosed by brush action
#' @param input shiny input object
#' @param output shiny output object
#' @param session shiny session object
#' @return nothing; function info label for click action in HM
# examples click_action_fullHM(df, input, output, session)
hover_action_fullHM <- function(df, input, output, session) {
  base::tryCatch(
    {
      base::print(base::paste0(sysTimePID(), " hover_action_fullHM.", as.character(head(df))))
    },
    error = function(e) {
      if (attributes(e)$class[1] != "shiny.silent.error") {
        base::message("An error occurred in hover_action_fullHM():\n", e)
      }
    },
    warning = function(w) {
      base::message("A warning occurred in hover_action_fullHM():\n", w)
    },
    finally = {
      base::print(base::paste0(sysTimePID(), " end hover_action_fullHM()."))
    }
  )
}

#' click_action_fullHM_LogFC
#' @param df data.frame containing data which is enclosed by brush action
#' @param input shiny input object
#' @param output shiny output object
#' @param session shiny session object
#' @return nothing; function info label for click action in HM
# examples click_action_fullHM_LogFC(df, input, output, session)
click_action_fullHM_LogFC <- function(df, input, output, session) {
  base::tryCatch(
    {
      base::print(base::paste0(sysTimePID(), " start click_action_fullHM_LogFC()."))
    },
    error = function(e) {
      if (attributes(e)$class[1] != "shiny.silent.error") {
        base::message("An error occurred in click_action_fullHM_LogFC():\n", e)
      }
    },
    warning = function(w) {
      base::message("A warning occurred in click_action_fullHM_LogFC():\n", w)
    },
    finally = {
      base::print(base::paste0(sysTimePID(), " end click_action_fullHM_LogFC()."))
    }
  )
}

#' brush_action_fullHM_LogFC
#' @param df data.frame containing data which is enclosed by brush action
#' @param input shiny input object
#' @param output shiny output object
#' @param session shiny session object
#' @return nothing; function creates SPLOM
# examples brush_action_fullHM_LogFC(df, input, output, session)
brush_action_fullHM_LogFC <- function(df, input, output, session) {
  base::tryCatch(
    if (!is.null(df)) {
      row_index <- collapse::funique(unlist(df$row_index))
      column_index <- collapse::funique(unlist(df$column_index))
      #feed in selected CpG here
      session$userData$sessionVariables$selectedProbe(rownames(session$userData$sessionVariables$probeReducedDataStructureLogFC()$combinedDFP_Val_Labels$dfLogFC)[row_index])
      originTrait <- session$userData$sessionVariables$probeReducedDataStructureLogFC()$combinedDFP_Val_Labels$mergedOriginTrait[column_index]
      traitLabels <- session$userData$sessionVariables$probeReducedDataStructureLogFC()$combinedDFP_Val_Labels$mergedOriginalColnames[column_index]
      selectedTrait <- cbind(traitLabels, originTrait)
      colnames(selectedTrait) <- c("traitName", "traitSource")
      session$userData$sessionVariables$selectedTrait(selectedTrait)
      selectedTraitID <- session$userData$sessionVariables$probeReducedDataStructureLogFC()$combinedDFP_Val_Labels$traitID[column_index]
      session$userData$sessionVariables$selectedTraitID(selectedTraitID)
      base::print(base::paste0(sysTimePID(), " last selection was from full heatmap log(FC) brush event."))
      output$txtGlobalSelectOut <- shiny::renderText("last selection was from full heatmap log(FC) brush event.")
    },
    error = function(e) {
      if (attributes(e)$class[1] != "shiny.silent.error") {
        base::message("An error occurred in brush_action_fullHM_LogFC():\n", e)
      }
    },
    warning = function(w) {
      base::message("A warning occurred in brush_action_fullHM_LogFC():\n", w)
    },
    finally = {
      base::print(base::paste0(sysTimePID(), " end brush_action_fullHM_LogFC()."))
    }
  )
}

#' hover_action_fullHM_LogFC
#' @param df data.frame containing data which is enclosed by brush action
#' @param input shiny input object
#' @param output shiny output object
#' @param session shiny session object
#' @return nothing; function info label for click action in HM
# examples hover_action_fullHM_LogFC(df, input, output, session)
hover_action_fullHM_LogFC <- function(df, input, output, session) {
  base::tryCatch(
    {
      base::print(base::paste0(sysTimePID(), " hover_action_fullHM_LogFC", as.character(head(df))))
    },
    error = function(e) {
      if (attributes(e)$class[1] != "shiny.silent.error") {
        base::message("An error occurred in hover_action_fullHM_LogFC():\n", e)
      }
    },
    warning = function(w) {
      base::message("A warning occurred in hover_action_fullHM_LogFC():\n", w)
    },
    finally = {
      base::print(base::paste0(sysTimePID(), " end hover_action_fullHM_LogFC()."))
    }
  )
}

click_action_condHM_P_Val <- function(df, input, output, session) {
  base::tryCatch(
    {
      base::print(base::paste0(sysTimePID(), " start click_action_condHM_P_Val()."))
      # only placeholder at the moment
    },
    error = function(e) {
      if (attributes(e)$class[1] != "shiny.silent.error") {
        base::message("An error occurred in click_action_condHM_P_Val():\n", e)
      }
    },
    warning = function(w) {
      base::message("A warning occurred in click_action_condHM_P_Val():\n", w)
    },
    finally = {
      base::print(base::paste0(sysTimePID(), " end click_action_condHM_P_Val()."))
    }
  )
}

#' brush_action_condHM_P_Val
#' @param df data.frame containing data which is enclosed by brush action
#' @param input shiny input object
#' @param output shiny output object
#' @param session shiny session object
#' @return nothing; function creates SPLOM
# examples brush_action_CondHM_P_Val(df, input, output, session)
brush_action_condHM_P_Val <- function(df, input, output, session) {
  base::tryCatch(
    if (!is.null(df)) {
      row_index <- collapse::funique(unlist(df$row_index)) #row_index <- unique(unlist(df$row_index))
      column_index <- collapse::funique(unlist(df$column_index))
      #feed in selected CpG here
      session$userData$sessionVariables$selectedProbe(rownames(session$userData$sessionVariables$probeReducedDataStructure()$combinedDFP_Val_Labels$dfP_Val)[row_index])
      originTrait <- session$userData$sessionVariables$probeReducedDataStructure()$combinedDFP_Val_Labels$mergedOriginTrait[column_index]
      traitLabels <- colnames(session$userData$sessionVariables$probeReducedDataStructure()$combinedDFP_Val_Labels$dfP_Val)[column_index]
      selectedTrait <- cbind(traitLabels, originTrait)
      colnames(selectedTrait) <- c("traitName", "traitSource")
      session$userData$sessionVariables$selectedTrait(selectedTrait)
      base::print(base::paste0(sysTimePID(), " last selection was from condensed heatmap p-val brush event."))
      output$txtGlobalSelectOut <- shiny::renderText("last selection was from condensed heatmap p-val brush event.")
    },
    error = function(e) {
      if (attributes(e)$class[1] != "shiny.silent.error") {
        base::message("An error occurred in brush_action_condHM_P_Val():\n", e)
      }
    },
    warning = function(w) {
      base::message("A warning occurred in brush_action_condHM_P_Val():\n", w)
    },
    finally = {
      base::print(base::paste0(sysTimePID(), " end brush_action_condHM_P_Val()."))
    }
  )
}

hover_action_condHM_P_Val <- function(df, input, output, session) {
  base::tryCatch(
    {
      base::print(base::paste0(sysTimePID(), " hover_action_condHM_P_Val.", as.character(head(df))))
      # only placeholder at the moment
    },
    error = function(e) {
      if (attributes(e)$class[1] != "shiny.silent.error") {
        base::message("An error occurred in hover_action_condHM_P_Val():\n", e)
      }
    },
    warning = function(w) {
      base::message("A warning occurred in hover_action_condHM_P_Val():\n", w)
    },
    finally = {
      base::print(base::paste0(sysTimePID(), " end hover_action_condHM_P_Val()."))
    }
  )
}

click_action_condHeatmap_LogFC <- function(df, input, output, session) {
  base::tryCatch(
    {
      base::print(base::paste0(sysTimePID(), " start click_action_condHeatmap_LogFC()."))
      # only placeholder at the moment
    },
    error = function(e) {
      if (attributes(e)$class[1] != "shiny.silent.error") {
        base::message("An error occurred in click_action_condHeatmap_LogFC():\n", e)
      }
    },
    warning = function(w) {
      base::message("A warning occurred in click_action_condHeatmap_LogFC():\n", w)
    },
    finally = {
      base::print(base::paste0(sysTimePID(), " end click_action_condHeatmap_LogFC()."))
    }
  )
}

brush_action_condHeatmap_LogFC <- function(df, input, output, session) {
  base::tryCatch(
    if (!is.null(df)) {
      row_index <- collapse::funique(unlist(df$row_index)) #row_index <- unique(unlist(df$row_index))
      column_index <- collapse::funique(unlist(df$column_index))
      #load results into session global data structure #feed in selected CpG here
      session$userData$sessionVariables$selectedProbe(rownames(session$userData$sessionVariables$probeReducedDataStructure()$combinedDFP_Val_Labels$dfLogFC)[row_index])
      originTrait <- session$userData$sessionVariables$probeReducedDataStructure()$combinedDFP_Val_Labels$mergedOriginTrait[column_index]
      traitLabels <- colnames(session$userData$sessionVariables$probeReducedDataStructure()$combinedDFP_Val_Labels$dfP_Val)[column_index]
      selectedTrait <- cbind(traitLabels, originTrait)
      colnames(selectedTrait) <- c("traitName", "traitSource")
      session$userData$sessionVariables$selectedTrait(selectedTrait)
      base::print(base::paste0(sysTimePID(), " last selection was from condensed heatmap logFC brush event."))
      output$txtGlobalSelectOut <- shiny::renderText("last selection was from condensed heatmap logFC brush event.")
    },
    error = function(e) {
      if (attributes(e)$class[1] != "shiny.silent.error") {
        base::message("An error occurred in brush_action_condHeatmap_LogFC():\n", e)
      }
    },
    warning = function(w) {
      base::message("A warning occurred in brush_action_condHeatmap_LogFC():\n", w)
    },
    finally = {
      base::print(base::paste0(sysTimePID(), " end brush_action_condHeatmap_LogFC()."))
    }
  )
}

hover_action_condHeatmap_LogFC <- function(df, input, output, session) {
  base::tryCatch(
    {
      base::print(base::paste0(sysTimePID(), " hover_action_condHeatmap_LogFC", as.character(head(df))))
      # only placeholder at the moment
    },
    error = function(e) {
      if (attributes(e)$class[1] != "shiny.silent.error") {
        base::message("An error occurred in hover_action_condHeatmap_LogFC():\n", e)
      }
    },
    warning = function(w) {
      base::message("A warning occurred in hover_action_condHeatmap_LogFC():\n", w)
    },
    finally = {
      base::print(base::paste0(sysTimePID(), " end hover_action_condHeatmap_LogFC()."))
    }
  )
}
