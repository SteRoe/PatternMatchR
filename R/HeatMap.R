HeatMap_UI <- function(id) {
  ns <- shiny::NS(id)
  htmltools::tagList(
    shiny::tabsetPanel(id = ns("tabsetHeatMap"),
      shiny::tabPanel(
       "Table P_VAL",
       if (id == "HeatMap_Full_DetailsPval") {
         "Table of p-value; clustering order comes from clustering of p-values."
       }
       else if (id == "HeatMap_Full_DetailsLogFC") {
         "Table of p-value; clustering order comes from clustering of log(FC)."
       },
       DT::dataTableOutput(ns("traitReducedDTP_VAL"))
      ),
      shiny::tabPanel(
       "Table Delta Methylation",
       if (id == "HeatMap_Full_DetailsPval") {
        "Table of delta methylation; clustering order comes from clustering of p-values."
       }
       else if (id == "HeatMap_Full_DetailsLogFC") {
         "Table of delta methylation; clustering order comes from clustering of log(FC)."
       },
       DT::dataTableOutput(ns("traitReducedDTDM"))
      ),
      shiny::tabPanel(
       "Table N",
       if (id == "HeatMap_Full_DetailsPval") {
        "Table of n; clustering order comes from clustering of p-values."
       }
       else if (id == "HeatMap_Full_DetailsLogFC") {
         "Table of n; clustering order comes from clustering of log(FC)."
       },
       DT::dataTableOutput(ns("traitReducedDTN"))
      ),
      shiny::tabPanel(
       "Table Delta Methylation log(FC)",
       if (id == "HeatMap_Full_DetailsPval") {
        "Table of log fold change(delta methylation); clustering order comes from clustering of p-values."
       }
       else if (id == "HeatMap_Full_DetailsLogFC") {
         "Table of log fold change(delta methylation); clustering order comes from clustering of log(FC)."
       },
       DT::dataTableOutput(ns("traitReducedDTLogFC"))
      )
    )
  )
}

HeatMap_SERVER <- function(id, session) {
  shiny::moduleServer(id, function(input, output, session) {
    if (id == "HeatMap_Full_DetailsPval") {
      output$traitReducedDTP_VAL <- DT::renderDataTable(as.data.frame(session$userData$sessionVariables$traitReducedDataStructurePVal()$combinedDFP_Val_Labels$dfP_Val_w_number))
      output$traitReducedDTDM <- DT::renderDataTable(as.data.frame(session$userData$sessionVariables$traitReducedDataStructurePVal()$combinedDFP_Val_Labels$dfDM_w_number))
      output$traitReducedDTLogFC <- DT::renderDataTable(as.data.frame(session$userData$sessionVariables$traitReducedDataStructurePVal()$combinedDFP_Val_Labels$dfLogFC_w_number))
      output$traitReducedDTN <- DT::renderDataTable(as.data.frame(session$userData$sessionVariables$traitReducedDataStructurePVal()$combinedDFP_Val_Labels$dfN_w_number))
    }
    else if (id == "HeatMap_Full_DetailsLogFC") {
      output$traitReducedDTP_VAL <- DT::renderDataTable(as.data.frame(session$userData$sessionVariables$traitReducedDataStructureLogFC()$combinedDFP_Val_Labels$dfP_Val_w_number))
      output$traitReducedDTDM <- DT::renderDataTable(as.data.frame(session$userData$sessionVariables$traitReducedDataStructureLogFC()$combinedDFP_Val_Labels$dfDM_w_number))
      output$traitReducedDTLogFC <- DT::renderDataTable(as.data.frame(session$userData$sessionVariables$traitReducedDataStructureLogFC()$combinedDFP_Val_Labels$dfLogFC_w_number))
      output$traitReducedDTN <- DT::renderDataTable(as.data.frame(session$userData$sessionVariables$traitReducedDataStructureLogFC()$combinedDFP_Val_Labels$dfN_w_number))
    }
    else{
      browser() #should not happen
    }
  }) #end shiny::moduleServer
}

plotCombinedHM_P_Val <- function(input, output, session) {
  id <- shiny::showNotification("Plotting heatmap for p-val...", duration = NULL, closeButton = FALSE)
  base::on.exit(shiny::removeNotification(id), add = TRUE)
  base::print(base::paste0(sysTimePID(), " start plotting heatmap for P_Val."))
  output$txtHMDescription_P_Val <-
    shiny::renderText(base::paste0("calculating heatmap..., current plot is not valid"))
  while (!is.null(grDevices::dev.list())) {
    grDevices::dev.off()
  }
  combinedDFP_Val_Labels <- session$userData$sessionVariables$traitReducedDataStructurePVal()$combinedDFP_Val_Labels
  dfP_Val <- combinedDFP_Val_Labels$dfP_Val
  #browser() #if step 3 was omitted, we see an error here...
  dfP_Val[dfP_Val > 0.05] <- NA # 1
  base::print(
    base::paste0(
      sysTimePID(),
      " calculating combined heatmap with rows= ",
      nrow(dfP_Val),
      " cols= ",
      ncol(dfP_Val)
    )
  )
  if (nrow(dfP_Val) > 5) {
    startTime <- Sys.time()
    base::tryCatch({
      base::print(base::paste0(sysTimePID(), " gc()"))
      gc()
      # check clustResProbes > 8
      base::options(expressions = 500000)
      dendProbes <- session$userData$sessionVariables$traitReducedDataStructurePVal()$probeDendrogram
      maxClassesProbes <- 7 #as.integer(input$txtMaxClassesProbes)
      if (maxClassesProbes <= 7) {
        dendProbes <- dendextend::color_branches(dendProbes, maxClassesProbes)
      }
      dendTraits <- session$userData$sessionVariables$traitReducedDataStructurePVal()$traitDendrogram
      Distances <- session$userData$sessionVariables$traitReducedDataStructurePVal()$DNAdistances

      base::print(base::paste0(sysTimePID(), " before calculating heatmap"))

      base::print(base::paste0(sysTimePID(), " length(unlist(dendProbes)): ", length(unlist(dendProbes))))
      base::print(base::paste0(sysTimePID(), " length(unlist(dendTraits)): ", length(unlist(dendTraits))))
      length(unlist(dendTraits)) == base::dim(combinedDFP_Val_Labels$dfP_Val)[2]
      length(unlist(dendProbes)) == base::dim(combinedDFP_Val_Labels$dfP_Val)[1]
      selectedRowIndicesYellow <- NULL
      if(is.valid(session$userData$sessionVariables$selectedProbe())) {
        selectedRowIndicesYellow <- session$userData$sessionVariables$selectedProbe()
      }
      # if (is.valid(input$txtSearchFullCpG)) {
      #   selectedRowIndicesYellow <- unlist(strsplit(input$txtSearchFullCpG, split = " ")) #is a list of cg-numbers from search field "txtSearchCpG"
      # }
      # if (is.valid(input$txtSearchFullTrait)) {
      #   selectedColIndices <- unlist(strsplit(input$txtSearchFullTrait, split = " "))
      # }
      #selectedRowIndicesOrange <- session$userData$sessionVariables$distancesBelowThreshold()
      selectedRowIndicesOrange <- NULL
      base::print(base::paste0(sysTimePID(), " before l <- combinedDFInteractiveHeatMapP_Val(combinedDFP_Val_Labels, dendProbes, dendTraits, selectedRowIndices, selectedColIndices)"))
      l <-
        combinedDFInteractiveHeatMapP_Val(combinedDFP_Val_Labels, dendProbes, dendTraits, Distances, selectedRowIndicesYellow, selectedColIndices, selectedRowIndicesOrange, session)
      base::print(base::paste0(sysTimePID(), " before combinedHM <- l$combinedHM"))
      combinedHM <- l$combinedHM
      endTime <- Sys.time()
      elapsedTime <- endTime - startTime
      base::print(base::paste0(sysTimePID(), " after calculating heatmap. Elapsed time: ", elapsedTime, " sec."))
      base::print(base::paste0(sysTimePID(), " before plotting heatmap; InteractiveComplexHeatmap::makeInteractiveComplexHeatmap()"))
      while (!base::is.null(grDevices::dev.list())) {
        grDevices::dev.off()
      }
      startTime <- Sys.time()
      InteractiveComplexHeatmap::makeInteractiveComplexHeatmap(
        input = input,
        output = output,
        session = session,
        ht_list = combinedHM,
        heatmap_id = "Heatmap_P_Val",
        show_layer_fun = TRUE,
        click_action = click_action_fullHM_P_Val,
        brush_action = brush_action_fullHM_P_Val,
        hover_action = hover_action_fullHM_P_Val
      )
    },
    error = function(e) {
      base::message("An error occurred in plotCombinedHM_P_Val():\n", e)
      Cstack_info()
      browser() #should not happen
    },
    warning = function(w) {
      base::message("A warning occurred in plotCombinedHM_P_Val():\n", w)
      browser() #should not happen
    },
    finally = {
      endTime <- Sys.time()
      elapsedTime <- endTime - startTime
      base::print(base::paste0(sysTimePID(), " after plotting heatmap for P_Val; InteractiveComplexHeatmap::makeInteractiveComplexHeatmap(). Elapsed time: ", elapsedTime, " sec."))
      output$txtHMDescription_P_Val <-
        shiny::renderText(
          base::paste0(
            sysTimePID(),
            " done plotting heatmap for P_Val..., current plot is valid. n(probe) = ",
            base::nrow(base::as.matrix(combinedDFP_Val_Labels[[1]])),
            "; n(trait) = ",
            base::ncol(base::as.matrix(combinedDFP_Val_Labels[[1]])),
            "; elapsed time: ",
            elapsedTime, " sec."
          )
        )
    })
  }
}

plotCombinedHM_LogFC <- function(input, output, session) {
  id <- shiny::showNotification("Plotting heatmap for log(FC)...", duration = NULL, closeButton = FALSE)
  base::on.exit(shiny::removeNotification(id), add = TRUE)
  base::print(base::paste0(sysTimePID(), " start plotting heatmap for log(FC)."))
  output$txtHMDescription_LogFC <-
    shiny::renderText(base::paste0("calculating heatmap..., current plot is not valid"))
  while (!is.null(grDevices::dev.list())) {
    grDevices::dev.off()
  }
  combinedDFP_Val_Labels <- session$userData$sessionVariables$traitReducedDataStructureLogFC()$combinedDFP_Val_Labels
  dfLogFC <- combinedDFP_Val_Labels$dfLogFC
  #browser() #if step 3 was omitted, we see an error here...
  dfLogFC[dfLogFC > 0.05] <- NA # 1
  base::print(
    base::paste0(
      sysTimePID(),
      " calculating combined heatmap with rows= ",
      nrow(dfLogFC),
      " cols= ",
      ncol(dfLogFC)
    )
  )
  if (nrow(dfLogFC) > 5) {
    startTime <- Sys.time()
    base::tryCatch({
      base::print(base::paste0(sysTimePID(), " gc()"))
      gc()
      # check clustResProbes > 8
      base::options(expressions = 500000)
      dendProbes <- session$userData$sessionVariables$traitReducedDataStructureLogFC()$probeDendrogram
      maxClassesProbes <- 7 #as.integer(input$txtMaxClassesProbes)
      if (maxClassesProbes <= 7) {
        dendProbes <- dendextend::color_branches(dendProbes, maxClassesProbes)
      }
      dendTraits <- session$userData$sessionVariables$traitReducedDataStructureLogFC()$traitDendrogram
      Distances <- session$userData$sessionVariables$traitReducedDataStructureLogFC()$DNAdistances

      base::print(base::paste0(sysTimePID(), " before calculating heatmap"))

      base::print(base::paste0(sysTimePID(), " length(unlist(dendProbes)): ", length(unlist(dendProbes))))
      base::print(base::paste0(sysTimePID(), " length(unlist(dendTraits)): ", length(unlist(dendTraits))))
      length(unlist(dendTraits)) == base::dim(combinedDFP_Val_Labels$dfLogFC)[2]
      length(unlist(dendProbes)) == base::dim(combinedDFP_Val_Labels$dfLogFC)[1]
      selectedRowIndicesYellow <- NULL
      if(is.valid(session$userData$sessionVariables$selectedProbe())) {
        selectedRowIndicesYellow <- session$userData$sessionVariables$selectedProbe()
      }
      selectedRowIndicesOrange <- NULL
      base::print(base::paste0(sysTimePID(), " before l <- combinedDFInteractiveHeatMapLogFC(combinedDFP_Val_Labels, dendProbes, dendTraits, selectedRowIndices, selectedColIndices)"))
      l <-
        combinedDFInteractiveHeatMapLogFC(combinedDFP_Val_Labels, dendProbes, dendTraits, Distances, selectedRowIndicesYellow, selectedColIndices, selectedRowIndicesOrange, session)
      base::print(base::paste0(sysTimePID(), " before combinedHMLogFC <- l$combinedHM"))
      combinedHM <- l$combinedHM
      endTime <- Sys.time()
      elapsedTime <- endTime - startTime
      base::print(base::paste0(sysTimePID(), " after calculating heatmap. Elapsed time: ", elapsedTime, " sec."))
      base::print(base::paste0(sysTimePID(), " before plotting heatmap; InteractiveComplexHeatmap::makeInteractiveComplexHeatmap()"))
      while (!base::is.null(grDevices::dev.list())) {
        grDevices::dev.off()
      }
      startTime <- Sys.time()
      InteractiveComplexHeatmap::makeInteractiveComplexHeatmap(
        input = input,
        output = output,
        session = session,
        ht_list = combinedHM,
        heatmap_id = "Heatmap_LogFC",
        show_layer_fun = TRUE,
        click_action = click_action_fullHM_LogFC,
        brush_action = brush_action_fullHM_LogFC,
        hover_action = hover_action_fullHM_LogFC
      )
    },
    error = function(e) {
      base::message("An error occurred in plotCombinedHM_LogFC():\n", e)
      Cstack_info()
      browser() #should not happen
    },
    warning = function(w) {
      base::message("A warning occurred in plotCombinedHM_LogFC():\n", w)
      browser() #should not happen
    },
    finally = {
      endTime <- Sys.time()
      elapsedTime <- endTime - startTime
      base::print(base::paste0(sysTimePID(), " after plotting heatmap for log(FC); InteractiveComplexHeatmap::makeInteractiveComplexHeatmap(). Elapsed time: ", elapsedTime, " sec."))
      output$txtHMDescription_LogFC <-
        shiny::renderText(
          base::paste0(
            sysTimePID(),
            " done plotting heatmap for log(FC)..., current plot is valid. n(probe) = ",
            base::nrow(base::as.matrix(combinedDFP_Val_Labels[[1]])),
            "; n(trait) = ",
            base::ncol(base::as.matrix(combinedDFP_Val_Labels[[1]])),
            "; elapsed time: ",
            elapsedTime, " sec."
          )
        )
    })
  }
}

getSelectedTraitReducedcombinedDFP_Val_Labels <- function(combinedDFP_Val_Labels, row_index, column_index, session) {
  id <- shiny::showNotification("Getting trait reduced combined data structure...", duration = NULL, closeButton = FALSE)
  base::on.exit(shiny::removeNotification(id), add = TRUE)
  base::print(base::paste0(sysTimePID(), " start getSelectedTraitReducedcombinedDFP_Val_Labels()"))
  base::tryCatch(
    {
      result <- list()
      result$dfP_Val <- combinedDFP_Val_Labels$dfP_Val[row_index, ]
      result$dfDM <- combinedDFP_Val_Labels$dfDM[row_index, ]
      result$dfN <- combinedDFP_Val_Labels$dfN[row_index, ]
      result$labelsDF1 <- combinedDFP_Val_Labels$labelsDF1
      result$labelsDF2 <- combinedDFP_Val_Labels$labelsDF2
      result$labelsDF3 <- combinedDFP_Val_Labels$labelsDF3
      result$mergedOriginDF <- combinedDFP_Val_Labels$mergedOriginDF
      result$mergedColnames <- combinedDFP_Val_Labels$mergedColnames
      result$mergedOriginTrait <- combinedDFP_Val_Labels$mergedOriginTrait
      result$mergedDFList <- combinedDFP_Val_Labels$mergedDFList
    },
    error = function(e) {
      base::message("An error occurred in getSelectedTraitReducedcombinedDFP_Val_Labels():\n", e)
    },
    warning = function(w) {
      base::message("A warning occurred in getSelectedTraitReducedcombinedDFP_Val_Labels():\n", w)
    },
    finally = {
      base::print(base::paste0(sysTimePID(), " end getSelectedTraitReducedcombinedDFP_Val_Labels()"))
      return(result)
    })
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

#' creates a regular heatmap based on log fold change values
#' @param combinedDF_Labels list of data.frame and labels generated from function mergeDFP_Val_Labels()
#' @param dendProbes dendrogram (without labels) for probes (rows), generated externally and providing information for clustering of heatmap
#' @param dendTraits dendrogram (without labels) result for traits (columns), generated externally and providing information for clustering of heatmap
#' @param selectedRowIndicesYellow indicies of HM rows to mark in yellow color
#' @param selectedColIndices indicies of HM cols to mark (not used so far)
#' @param selectedRowIndicesOrange indicies of HM rows to mark in orange color
#' @param session session object for reference
#' @return heatmap object for InteractiveComplexHeatmap::makeInteractiveComplexHeatmap
# examples combinedDFInteractiveHeatMapDMlogFC(combinedDF_Labels, dendProbes, dendTraits, Distances, selectedRowIndicesYellow, selectedColIndices, selectedRowIndicesOrange)
combinedDFInteractiveHeatMapDMLogFC <-
  function(combinedDF_Labels,
           dendProbes = NA,
           dendTraits = NA,
           Distances = NA,
           selectedRowIndicesYellow = NA,
           selectedColIndices = NA,
           selectedRowIndicesOrange = NA,
           session = session) {
    id <- shiny::showNotification("Plotting heatmap for log(FC)...", duration = NULL, closeButton = FALSE)
    base::on.exit(shiny::removeNotification(id), add = TRUE)
    base::tryCatch(
      {
        startTime <- Sys.time()
        base::print(base::paste0(sysTimePID(), " start preparing HM; combinedDFInteractiveHeatMapP_Val()"))
        matP_Val <- base::as.matrix(combinedDF_Labels$dfP_Val)
        matDM <- base::as.matrix(combinedDF_Labels$dfDM)
        matLogFC <- base::as.matrix(combinedDF_Labels$dfLogFC)
        matN <- base::as.matrix(combinedDF_Labels$dfN)
        # use rasterization like described in
        # https://jokergoo.github.io/2020/06/30/rasterization-in-complexheatmap/
        base::print(base::paste0(sysTimePID(), " making labels"))
        labelsDF1 <- combinedDF_Labels$labelsDF1
        labelsDF2 <- combinedDF_Labels$labelsDF2
        labelsDF3 <- combinedDF_Labels$labelsDF3
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
        max.Col1 <- base::max(matLogFC, na.rm = TRUE)
        min.Col1 <- base::min(matLogFC, na.rm = TRUE)
        max.Col2 <- base::max(matLogFC, na.rm = TRUE)
        min.Col2 <- base::min(matLogFC, na.rm = TRUE)
        max.Col3 <- base::max(matLogFC, na.rm = TRUE)
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
                rowsToMarkYellow <- labelsProbes %in% selectedRowIndicesYellow
                if (any(rowsToMarkYellow)) {
                  grid::grid.rect(x, y[rowsToMarkYellow], w, h[rowsToMarkYellow], gp = grid::gpar(
                    fill = "yellow",
                    col = "yellow"
                  ))
                }
                rowsToMarkOrange <- labelsProbes %in% selectedRowIndicesOrange
                if (any(rowsToMarkOrange)) {
                  grid::grid.rect(x, y[rowsToMarkOrange], w, h[rowsToMarkOrange], gp = grid::gpar(
                    fill = "orange",
                    col = "orange"
                  ))
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
                      "dLogFC:",
                      sprintf("%.G", ComplexHeatmap::pindex(matLogFC, i, j)),
                      "\n",
                      "n:",
                      ComplexHeatmap::pindex(matN, i, j)
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
        base::message("An error occurred in combinedDFInteractiveHeatMapDMLogFC():\n", e)
      },
      warning = function(w) {
        base::message("A warning occurred in combinedDFInteractiveHeatMapDMLogFC():\n", w)
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
        #lgd <- list(lgdHMDist[1], lgd[1], lgd[2], lgd[3])
      },
      error = function(e) {
        base::message("An error occurred in making htDistances():\n", e)
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
        base::message("An error occurred in clearing grDevices():\n", e)
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
        base::print(base::paste0(sysTimePID(), " start drawing heatmap (takes some time). (step before ComplexHeatmap::draw()"))
        # with huge heatmaps, the following error occurs:
        # Error in Cairo: Failed to create Cairo backend!
        # ht <- ComplexHeatmap::draw(ht + ht2 + ht3, annotation_legend_list = lgd)
        # ht <- ComplexHeatmap::draw(ht, annotation_legend_list = lgd)
        ht <- ComplexHeatmap::draw(htDistances + ht, annotation_legend_list = lgd)
      },
      error = function(err) {
        base::message(base::paste0(sysTimePID(), " Error: unable to draw HM. ", err$message))
      },
      warning = function(w) {
        base::message(base::paste0(sysTimePID(), " unable to draw HM. ", w$message))
      },
      finally = {
        base::print(base::paste0(sysTimePID(), " grDevices::dev.off()"))
        grDevices::dev.off()
        base::print(base::paste0(sysTimePID(), " l <- base::list()"))
        l <- base::list()
        base::print(base::paste0(sysTimePID(), " l$combinedHMDMLogFC <- ht"))
        l$combinedHMDMLogFC <- ht
        endTime <- Sys.time()
        elapsedTime <- endTime - startTime
        base::print(base::paste0(sysTimePID(), " end combinedDFInteractiveHeatMapDMLogFC()"))
        base::print(base::paste0(sysTimePID(), " finished drawing heatmap (takes some time). (step after ComplexHeatmap::draw(); Elapsed time: ", elapsedTime, "."))
        base::return(l)
      })
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
combinedDFInteractiveHeatMapP_Val <-
  function(combinedDF_Labels,
           dendProbes = NA,
           dendTraits = NA,
           Distances = NA,
           selectedRowIndicesYellow = NA,
           selectedColIndices = NA,
           selectedRowIndicesOrange = NA,
           session = session) {
    #function is not called twice (check caller function), if so, check here
    id <- shiny::showNotification("Creating heatmap for (p-val)...", duration = NULL, closeButton = FALSE)
    base::on.exit(shiny::removeNotification(id), add = TRUE)
    base::tryCatch(
      {
        startTime <- Sys.time()
        base::print(base::paste0(sysTimePID(), " start preparing HM; combinedDFInteractiveHeatMapP_Val()"))
        matP_Val <- base::as.matrix(combinedDF_Labels$dfP_Val)
        matDM <- base::as.matrix(combinedDF_Labels$dfDM)
        matN <- base::as.matrix(combinedDF_Labels$dfN)
        # use rasterization like described in
        # https://jokergoo.github.io/2020/06/30/rasterization-in-complexheatmap/
        base::print(base::paste0(sysTimePID(), " making labels"))
        labelsDF1 <- combinedDF_Labels$labelsDF1
        labelsDF2 <- combinedDF_Labels$labelsDF2
        labelsDF3 <- combinedDF_Labels$labelsDF3
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
        min.Col1 <- base::min(matP_Val, na.rm = TRUE)
        max.Col2 <- 0.05
        min.Col2 <- base::min(matP_Val, na.rm = TRUE)
        max.Col3 <- 0.05
        min.Col3 <- base::min(matP_Val, na.rm = TRUE)
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
        base::print(base::paste0(sysTimePID(), " preparing heatmap n(probes)=", base::dim(matP_Val)[1], " x n(traits)=", base::dim(matP_Val)[2]))
        if ("dendrogram" %in% class(dendProbes) && "dendrogram" %in% class(dendTraits)) {
          base::print(base::paste0(sysTimePID(), " length(unlist(dendProbes)): ", length(unlist(dendProbes))))
          base::print(base::paste0(sysTimePID(), " length(unlist(dendTraits)): ", length(unlist(dendTraits))))
          length(unlist(dendTraits)) == base::dim(matP_Val)[2]
          length(unlist(dendProbes)) == base::dim(matP_Val)[1]
          ht <-
            ComplexHeatmap::Heatmap(
              matP_Val,
              rect_gp = grid::gpar(type = "none"),
              cluster_rows = dendProbes,
              cluster_columns = dendTraits,
              top_annotation = ha,
              layer_fun = function(j, i, x, y, w, h, fill) {
                if (length(i) == base::nrow(matP_Val) * base::ncol(matP_Val)) {
                  # we are in main HM
                  subHM <- FALSE
                } else {
                  # we are in sub HM
                  subHM <- TRUE
                }
                l <- labels[j] == "trait 1"
                if (any(l)) {
                  grid::grid.rect(x[l], y[l], w[l], h[l], gp = grid::gpar(
                    fill = col1(ComplexHeatmap::pindex(matP_Val, i[l], j[l])),
                    col = col1(ComplexHeatmap::pindex(matP_Val, i[l], j[l]))
                  ))
                }
                l <- labels[j] == "trait 2"
                if (any(l)) {
                  grid::grid.rect(x[l], y[l], w[l], h[l], gp = grid::gpar(
                    fill = col2(ComplexHeatmap::pindex(matP_Val, i[l], j[l])),
                    col = col2(ComplexHeatmap::pindex(matP_Val, i[l], j[l]))
                  ))
                }
                l <- labels[j] == "trait 3"
                if (any(l)) {
                  grid::grid.rect(x[l], y[l], w[l], h[l], gp = grid::gpar(
                    fill = col3(ComplexHeatmap::pindex(matP_Val, i[l], j[l])),
                    col = col3(ComplexHeatmap::pindex(matP_Val, i[l], j[l]))
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
                      ComplexHeatmap::pindex(matN, i, j)
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
        base::message("An error occurred in combinedDFInteractiveHeatMapP_Val():\n", e)
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
        #lgd <- list(lgdHMDist[1], lgd[1], lgd[2], lgd[3])
      },
      error = function(e) {
        base::message("An error occurred in making htDistances():\n", e)
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
        base::message("An error occurred in clearing grDevices():\n", e)
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
        base::print(base::paste0(sysTimePID(), " start drawing heatmap (takes some time). (step before ComplexHeatmap::draw()"))
        # with huge heatmaps, the following error occurs:
        # Error in Cairo: Failed to create Cairo backend!
        # ht <- ComplexHeatmap::draw(ht + ht2 + ht3, annotation_legend_list = lgd)
        # ht <- ComplexHeatmap::draw(ht, annotation_legend_list = lgd)
        ht <- ComplexHeatmap::draw(htDistances + ht, annotation_legend_list = lgd)
      },
      error = function(err) {
        base::message(base::paste0(sysTimePID(), " Error: unable to draw HM. ", err$message))
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
        base::print(base::paste0(sysTimePID(), " end combinedDFInteractiveHeatMapP_Val()"))
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
    id <- shiny::showNotification("Creating heatmap for (log(FC))...", duration = NULL, closeButton = FALSE)
    base::on.exit(shiny::removeNotification(id), add = TRUE)
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
        labelsDF1 <- combinedDF_Labels$labelsDF1
        labelsDF2 <- combinedDF_Labels$labelsDF2
        labelsDF3 <- combinedDF_Labels$labelsDF3
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
        base::message("An error occurred in combinedDFInteractiveHeatMapLogFC():\n", e)
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
        base::message("An error occurred in making htDistances():\n", e)
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
        base::message("An error occurred in clearing grDevices():\n", e)
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
        base::print(base::paste0(sysTimePID(), " start drawing heatmap (takes some time). (step before ComplexHeatmap::draw()"))
        # with huge heatmaps, the following error occurs:
        # Error in Cairo: Failed to create Cairo backend!
        # ht <- ComplexHeatmap::draw(ht + ht2 + ht3, annotation_legend_list = lgd)
        # ht <- ComplexHeatmap::draw(ht, annotation_legend_list = lgd)
        ht <- ComplexHeatmap::draw(htDistances + ht, annotation_legend_list = lgd)
      },
      error = function(err) {
        base::message(base::paste0(sysTimePID(), " Error: unable to draw HM. ", err$message))
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
    id <- shiny::showNotification("Creating heatmap for probe distances...", duration = NULL, closeButton = FALSE)
    base::on.exit(shiny::removeNotification(id), add = TRUE)
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
        base::message("An error occurred in HeatMapDistances():\n", e)
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
        base::message("An error occurred in clearing grDevices():\n", e)
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
        base::print(base::paste0(sysTimePID(), " start drawing heatmap (takes some time). (step before ComplexHeatmap::draw()"))
        # with huge heatmaps, the following error occurs:
        # Error in Cairo: Failed to create Cairo backend!
        ht <- ComplexHeatmap::draw(ht)
      },
      error = function(err) {
        base::message(base::paste0(sysTimePID(), " Error: unable to draw HM. ", err$message))
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
      base::message("An error occurred in getSearchResultCpG():\n", e)
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
      base::message("An error occurred in getSearchResultTrait():\n", e)
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

plotCombinedHM_DMLogFC <- function(input, output, session) {
  id <- shiny::showNotification("Plotting heatmap for log(FC)...", duration = NULL, closeButton = FALSE)
  base::on.exit(shiny::removeNotification(id), add = TRUE)
  base::print(base::paste0(sysTimePID(), " start plotting heatmap for log(FC)."))
  output$txtCondHMDescription_DM <-
    shiny::renderText(base::paste0("calculating heatmap..., current plot is not valid"))
  while (!is.null(grDevices::dev.list())) {
    grDevices::dev.off()
  }
  combinedDFP_Val_Labels <- session$userData$sessionVariables$probeReducedDataStructure()$combinedDFP_Val_Labels
  dfLogFC <- combinedDFP_Val_Labels$dfLogFC
  #leave out low logFC?

  if (nrow(dfLogFC) > 5) {
    startTime <- Sys.time()
  }
  base::tryCatch({
    base::print(base::paste0(sysTimePID(), " gc()"))
    gc()
    base::options(expressions = 500000)
    dendProbes <- session$userData$sessionVariables$probeReducedDataStructure()$probeDendrogram
    maxClassesProbes <- as.integer(input$txtMaxClassesProbes)
    if (maxClassesProbes <= 7) {
      dendProbes <- dendextend::color_branches(dendProbes, maxClassesProbes)
    }
    dendTraits <- session$userData$sessionVariables$probeReducedDataStructure()$traitDendrogram
    Distances <- session$userData$sessionVariables$probeReducedDataStructure()$DNAdistances

    base::print(base::paste0(sysTimePID(), " before calculating heatmap"))

    selectedRowIndicesYellow <- unlist(strsplit(input$txtSearchFullCpG, split = " ")) #is a list of cg-numbers from search field "txtSearchCpG"
    selectedColIndices <- unlist(strsplit(input$txtSearchFullTrait, split = " "))
    #selectedRowIndicesOrange <- session$userData$sessionVariables$distancesBelowThreshold()
    selectedRowIndicesOrange <- NULL
    base::print(base::paste0(sysTimePID(), " before l <- combinedDFInteractiveHeatMapP_Val(combinedDFP_Val_Labels, dendProbes, dendTraits, selectedRowIndices, selectedColIndices)"))
    l <-
      combinedDFInteractiveHeatMapDMLogFC(combinedDFP_Val_Labels, dendProbes, dendTraits, Distances, selectedRowIndicesYellow, selectedColIndices, selectedRowIndicesOrange, session)
    base::print(base::paste0(sysTimePID(), " before combinedHMDMLogFC <- l$combinedHMDMLogFC"))
    combinedHMDMLogFC <- l$combinedHMDMLogFC

    endTime <- Sys.time()
    elapsedTime <- endTime - startTime
    base::print(base::paste0(sysTimePID(), " after calculating heatmap. Elapsed time: ", elapsedTime, " sec."))
    base::print(base::paste0(sysTimePID(), " before plotting heatmap; InteractiveComplexHeatmap::makeInteractiveComplexHeatmap()"))
    while (!base::is.null(grDevices::dev.list())) {
      grDevices::dev.off()
    }
    startTime <- Sys.time()
    InteractiveComplexHeatmap::makeInteractiveComplexHeatmap(
      input = input,
      output = output,
      session = session,
      ht_list = combinedHMDMLogFC,
      heatmap_id = "condHeatmap_LogFC",
      show_layer_fun = TRUE,
      click_action = click_action_condHeatmap_LogFC,
      brush_action = brush_action_condHeatmap_LogFC,
      hover_action = hover_action_condHeatmap_LogFC
    )
  },
  error = function(e) {
    base::message("An error occurred in plotCombinedHM_LogFC():\n", e)
    Cstack_info()
    browser() #should not happen
  },
  warning = function(w) {
    base::message("A warning occurred in plotCombinedHM_LogFC():\n", w)
    browser() #should not happen
  },
  finally = {
    endTime <- Sys.time()
    elapsedTime <- endTime - startTime
    base::print(base::paste0(sysTimePID(), " after plotting heatmap for DM; InteractiveComplexHeatmap::makeInteractiveComplexHeatmap(). Elapsed time: ", elapsedTime, " sec."))
    output$txtCondHMDescription_DM <-
      shiny::renderText(
        base::paste0(
          sysTimePID(),
          " done plotting heatmap for DM..., current plot is valid. n(probe) = ",
          base::nrow(base::as.matrix(combinedDFP_Val_Labels[[1]])),
          "; n(trait) = ",
          base::ncol(base::as.matrix(combinedDFP_Val_Labels[[1]])),
          "; elapsed time: ",
          elapsedTime, " sec."
        )
      )
  })
}

#' #' plotCombinedDWHM_P_Val
#' #' plots distance weighted heatmap
#' #' @param input shiny input
#' #' @param output shiny output
#' #' @param session shiny session
#' #' @return nothing
#' #' examples plotCombinedDWHM_P_Val(input, output, session)
#' plotCombinedDWHM_P_Val <- function(input, output, session) {
#'   id <- shiny::showNotification("Plotting distance weighted heatmap for p-val...", duration = NULL, closeButton = FALSE)
#'   base::on.exit(shiny::removeNotification(id), add = TRUE)
#'   base::print(base::paste0(sysTimePID(), " start plotting distance weighted heatmap for P_Val."))
#'   output$txtDWHMDescription_P_Val <-
#'     shiny::renderText(base::paste0("calculating DW heatmap..., current plot is not valid"))
#'   while (!is.null(grDevices::dev.list())) {
#'     grDevices::dev.off()
#'   }
#'   combinedDFP_Val_Labels <- session$userData$sessionVariables$distanceMultipliedTraitReducedDataStructurePVal()$combinedDFP_Val_Labels
#'
#'   dfP_Val <- combinedDFP_Val_Labels$dfP_Val
#'   #browser() #if step 3 was omitted, we see an error here...
#'   #  dfP_Val[dfP_Val > 0.05] <- NA # 1
#'   base::print(
#'     base::paste0(
#'       sysTimePID(),
#'       " calculating combined heatmap with rows= ",
#'       nrow(dfP_Val),
#'       " cols= ",
#'       ncol(dfP_Val)
#'     )
#'   )
#'   base::print(base::class(dfP_Val))
#'   if (nrow(dfP_Val) > 5) {
#'     startTime <- Sys.time()
#'     base::print(base::paste0(sysTimePID(), " gc()"))
#'     gc()
#'     base::tryCatch({
#'       # check clustResProbes > 8
#'       #base::length(session$userData$sessionVariables$traitReducedDataStructurePVal()$clustResProbes)
#'       #base::options(expression = 500000)
#'       base::options(expressions = 500000)
#'       dendProbes <- session$userData$sessionVariables$distanceMultipliedTraitReducedDataStructurePVal()$probeDendrogram
#'       maxClassesProbes <- as.integer(input$txtMaxClassesProbes)
#'       if (maxClassesProbes <= 7) {
#'         dendProbes <- dendextend::color_branches(dendProbes, maxClassesProbes)
#'       }
#'       dendTraits <- session$userData$sessionVariables$distanceMultipliedTraitReducedDataStructurePVal()$traitDendrogram
#'
#'       Distances <- session$userData$sessionVariables$distanceMultipliedTraitReducedDataStructurePVal()$DNAdistances
#'
#'       base::print(base::paste0(sysTimePID(), " before calculating heatmap"))
#'
#'       base::print(base::paste0(sysTimePID(), " length(unlist(dendProbes)): ", length(unlist(dendProbes))))
#'       base::print(base::paste0(sysTimePID(), " length(unlist(dendTraits)): ", length(unlist(dendTraits))))
#'       length(unlist(dendTraits)) == base::dim(combinedDFP_Val_Labels$dfP_Val)[2]
#'       length(unlist(dendProbes)) == base::dim(combinedDFP_Val_Labels$dfP_Val)[1]
#'       selectedRowIndicesYellow <- unlist(strsplit(input$txtSearchFullCpG, split = " "))
#'       base::message(base::paste0(sysTimePID(), " length(selectedRowIndicesYellow): ", length(selectedRowIndicesYellow)))
#'       selectedColIndices <- unlist(strsplit(input$txtSearchFullTrait, split = " "))
#'       selectedRowIndicesOrange <- session$userData$sessionVariables$distancesBelowThreshold()
#'       base::message(base::paste0(sysTimePID(), " length(selectedRowIndicesOrange): ", length(selectedRowIndicesOrange)))
#'       base::print(base::paste0(sysTimePID(), " before l <- combinedDFInteractiveHeatMapP_Val(combinedDFP_Val_Labels, dendProbes, dendTraits, selectedRowIndices, selectedColIndices)"))
#'       #browser() #check, whether this is called twice
#'       l <-
#'         combinedDFInteractiveHeatMapP_Val(combinedDFP_Val_Labels, dendProbes, dendTraits, Distances, selectedRowIndicesYellow, selectedColIndices, selectedRowIndicesOrange, session)
#'       base::print(base::paste0(sysTimePID(), " before combinedHM <- l$combinedHM"))
#'       combinedHM <- l$combinedHM
#'
#'       endTime <- Sys.time()
#'       elapsedTime <- endTime - startTime
#'       base::print(base::paste0(sysTimePID(), " after calculating heatmap. Elapsed time: ", elapsedTime, " sec."))
#'       base::print(base::paste0(sysTimePID(), " before plotting heatmap; InteractiveComplexHeatmap::makeInteractiveComplexHeatmap()"))
#'       while (!base::is.null(grDevices::dev.list())) {
#'         grDevices::dev.off()
#'       }
#'       startTime <- Sys.time()
#'       InteractiveComplexHeatmap::makeInteractiveComplexHeatmap(
#'         input = input,
#'         output = output,
#'         session = session,
#'         ht_list = combinedHM,
#'         heatmap_id = "DWHeatmap_P_Val",
#'         show_layer_fun = TRUE
#'         #        click_action = click_action_HM_P_Val,
#'         #        brush_action = brush_action_HM_P_Val,
#'         #        hover_action = hover_action_HM_P_Val
#'       )
#'     },
#'     error = function(e) {
#'       base::message("An error occurred in plotCombinedDWHM_P_Val():\n", e)
#'       Cstack_info()
#'       browser() #should not happen
#'     },
#'     warning = function(w) {
#'       base::message("A warning occurred in plotCombinedDWHM_P_Val():\n", w)
#'       browser() #should not happen
#'     },
#'     finally = {
#'       endTime <- Sys.time()
#'       elapsedTime <- endTime - startTime
#'       base::print(base::paste0(sysTimePID(), " after plotting heatmap for DWP_Val; InteractiveComplexHeatmap::makeInteractiveComplexHeatmap(). Elapsed time: ", elapsedTime, " sec."))
#'       output$txtDWHMDescription_P_Val <-
#'         shiny::renderText(
#'           base::paste0(
#'             sysTimePID(),
#'             " done plotting heatmap for DWP_Val..., current plot is valid. n(probe) = ",
#'             base::nrow(base::as.matrix(combinedDFP_Val_Labels[[1]])),
#'             "; n(trait) = ",
#'             base::ncol(base::as.matrix(combinedDFP_Val_Labels[[1]])),
#'             "; elapsed time: ",
#'             elapsedTime, " sec."
#'           )
#'         )
#'     })
#'   }
#' }

#' plotCombinedCondHM_P_Val
#' plots condensed heatmap
#' @param input shiny input
#' @param output shiny output
#' @param session shiny session
#' @return nothing
#' examples plotCombinedCondHM_P_Val(input, output, session)
plotCombinedCondHM_P_Val <- function(input, output, session) {
  id <- shiny::showNotification("Plotting condensed heatmap for p-val...", duration = NULL, closeButton = FALSE)
  base::on.exit(shiny::removeNotification(id), add = TRUE)
  base::tryCatch({
    base::print(base::paste0(sysTimePID(), " start plotting condensed heatmap for P_Val."))
    output$txtCondHMDescription_P_Val <-
      shiny::renderText(base::paste0("calculating condensed heatmap..., current plot is not valid"))
    while (!is.null(grDevices::dev.list())) {
      grDevices::dev.off()
    }
    #combinedDFP_Val_Labels <- session$userData$sessionVariables$probeReducedDataStructure(numNeighbours = 10)$combinedDFP_Val_Labels
    combinedDFP_Val_Labels <- session$userData$sessionVariables$probeReducedDataStructure()$combinedDFP_Val_Labels
    dfP_Val <- combinedDFP_Val_Labels$dfP_Val
    #browser() #if step 3 was omitted, we see an error here...

    base::print(base::paste0(sysTimePID(), " calculating combined heatmap with rows= ", nrow(dfP_Val), " cols= ", ncol(dfP_Val)))
    if (nrow(dfP_Val) > 5) {
      startTime <- Sys.time()
      base::print(base::paste0(sysTimePID(), " gc()"))
      gc()
      # check clustResProbes > 8
      base::options(expressions = 500000)

      #      if (is.valid(combinedDF_Labels$dfP_Val)) {
      dendProbes <- session$userData$sessionVariables$probeReducedDataStructure()$probeDendrogram
      #browser() #we can either try to make a subset of dendrogram or to create a new dendrogram from the base data from the heatmap... we then need also a new clustering for the subset...
      maxClassesProbes <- as.integer(input$txtMaxClassesProbes)
      if (maxClassesProbes <= 7) {
        dendProbes <- dendextend::color_branches(dendProbes, maxClassesProbes)
      }
      dendTraits <- session$userData$sessionVariables$probeReducedDataStructure()$traitDendrogram

      Distances <- session$userData$sessionVariables$probeReducedDataStructure()$DNAdistances

      base::print(base::paste0(sysTimePID(), " before calculating condensed heatmap"))
      base::print(base::paste0(sysTimePID(), " length(unlist(dendProbes)): ", length(unlist(dendProbes))))
      base::print(base::paste0(sysTimePID(), " length(unlist(dendTraits)): ", length(unlist(dendTraits))))
      selectedRowIndicesYellow <- unlist(strsplit(input$txtSearchFullCpG, split = " "))
      base::message(base::paste0(sysTimePID(), " length(selectedRowIndicesYellow): ", length(selectedRowIndicesYellow)))
      selectedColIndices <- unlist(strsplit(input$txtSearchFullTrait, split = " "))
      selectedRowIndicesOrange <- session$userData$sessionVariables$distancesBelowThreshold()
      base::message(base::paste0(sysTimePID(), " length(selectedRowIndicesOrange): ", length(selectedRowIndicesOrange)))
      base::print(base::paste0(sysTimePID(), " before l <- combinedDFInteractiveHeatMapP_Val(combinedDFP_Val_Labels, dendProbes, dendTraits, selectedRowIndices, selectedColIndices)"))
      #browser() #check, whether this is called twice
      l <-
        combinedDFInteractiveHeatMapP_Val(combinedDFP_Val_Labels, dendProbes, dendTraits, Distances, selectedRowIndicesYellow, selectedColIndices, selectedRowIndicesOrange, session)
      base::print(base::paste0(sysTimePID(), " before combinedHM <- l$combinedHM"))
      combinedHM <- l$combinedHM

      endTime <- Sys.time()
      elapsedTime <- endTime - startTime
      base::print(base::paste0(sysTimePID(), " after calculating condensed heatmap. Elapsed time: ", elapsedTime, " sec."))
      base::print(base::paste0(sysTimePID(), " before plotting condensed heatmap; InteractiveComplexHeatmap::makeInteractiveComplexHeatmap()"))
      while (!base::is.null(grDevices::dev.list())) {
        grDevices::dev.off()
      }
      startTime <- Sys.time()
      InteractiveComplexHeatmap::makeInteractiveComplexHeatmap(
        input = input,
        output = output,
        session = session,
        ht_list = combinedHM,
        heatmap_id = "condHeatmap_P_Val",
        show_layer_fun = TRUE,
        click_action = click_action_condHM_P_Val,
        brush_action = brush_action_condHM_P_Val,
        hover_action = hover_action_condHM_P_Val
      )
    }
  },
  error = function(e) {
    base::message("An error occurred in plotCombinedHM_P_Val():\n", e)
    Cstack_info()
    browser() #should not happen
  },
  warning = function(w) {
    base::message("A warning occurred in plotCombinedHM_P_Val():\n", w)
    browser() #should not happen
  },
  finally = {
    endTime <- Sys.time()
    elapsedTime <- endTime - startTime
    base::print(base::paste0(sysTimePID(), " after plotting heatmap for P_Val; InteractiveComplexHeatmap::makeInteractiveComplexHeatmap(). Elapsed time: ", elapsedTime, " sec."))
    output$txtCondHMDescription_P_Val <-
      shiny::renderText(
        base::paste0(
          sysTimePID(),
          " done plotting heatmap for condensed P_Val..., current plot is valid. n(probe) = ",
          base::nrow(base::as.matrix(combinedDFP_Val_Labels[[1]])),
          "; n(trait) = ",
          base::ncol(base::as.matrix(combinedDFP_Val_Labels[[1]])),
          "; elapsed time: ",
          elapsedTime, " sec."
        )
      )
  })
}

#' #' plotCombinedCondDWHM_P_Val
#' #' plots condensed distance weighted heatmap
#' #' @param input shiny input
#' #' @param output shiny output
#' #' @param session shiny session
#' #' @return nothing
#' #' examples plotCombinedCondDWHM_P_Val(input, output, session)
#' plotCombinedCondDWHM_P_Val <- function(input, output, session) {
#'   id <- shiny::showNotification("Plotting condensed distance weighted heatmap for p-val...", duration = NULL, closeButton = FALSE)
#'   base::on.exit(shiny::removeNotification(id), add = TRUE)
#'   base::tryCatch({
#'     base::print(base::paste0(sysTimePID(), " start plotting condensed distance weighted heatmap for P_Val."))
#'     output$txtCondDWHMDescription_P_Val <-
#'       shiny::renderText(base::paste0("calculating condensedDW heatmap..., current plot is not valid"))
#'     while (!is.null(grDevices::dev.list())) {
#'       grDevices::dev.off()
#'     }
#'     #combinedDFP_Val_Labels <- session$userData$sessionVariables$probeReducedDataStructure(numNeighbours = 10)$combinedDFP_Val_Labels
#'     combinedDFP_Val_Labels <- session$userData$sessionVariables$distanceMultipliedProbeReducedDataStructure()$combinedDFP_Val_Labels
#'     dfP_Val <- combinedDFP_Val_Labels$dfP_Val
#'     #browser() #if step 3 was omitted, we see an error here...
#'
#'     base::print(base::paste0(sysTimePID(), " calculating combined heatmap with rows= ", nrow(dfP_Val), " cols= ", ncol(dfP_Val)))
#'     base::print(base::class(dfP_Val))
#'     if (nrow(dfP_Val) > 5) {
#'       startTime <- Sys.time()
#'       base::print(base::paste0(sysTimePID(), " gc()"))
#'       gc()
#'       # check clustResProbes > 8
#'       base::options(expressions = 500000)
#'
#'       #      if (is.valid(combinedDF_Labels$dfP_Val)) {
#'       dendProbes <- session$userData$sessionVariables$distanceMultipliedProbeReducedDataStructure()$probeDendrogram
#'       #browser() #we can either try to make a subset of dendrogram or to create a new dendrogram from the base data from the heatmap... we then need also a new clustering for the subset...
#'       maxClassesProbes <- as.integer(input$txtMaxClassesProbes)
#'       if (maxClassesProbes <= 7) {
#'         dendProbes <- dendextend::color_branches(dendProbes, as.integer(maxClassesProbes))
#'       }
#'       dendTraits <- session$userData$sessionVariables$distanceMultipliedProbeReducedDataStructure()$traitDendrogram
#'
#'       Distances <- session$userData$sessionVariables$distanceMultipliedProbeReducedDataStructure()$DNAdistances
#'
#'       base::print(base::paste0(sysTimePID(), " before calculating condensed DW heatmap"))
#'       base::print(base::paste0(sysTimePID(), " length(unlist(dendProbes)): ", length(unlist(dendProbes))))
#'       base::print(base::paste0(sysTimePID(), " length(unlist(dendTraits)): ", length(unlist(dendTraits))))
#'       selectedRowIndicesYellow <- unlist(strsplit(input$txtSearchFullCpG, split = " "))
#'       base::message(base::paste0(sysTimePID(), " length(selectedRowIndicesYellow): ", length(selectedRowIndicesYellow)))
#'       selectedColIndices <- unlist(strsplit(input$txtSearchFullTrait, split = " "))
#'       selectedRowIndicesOrange <- session$userData$sessionVariables$distancesBelowThreshold()
#'       base::message(base::paste0(sysTimePID(), " length(selectedRowIndicesOrange): ", length(selectedRowIndicesOrange)))
#'       base::print(base::paste0(sysTimePID(), " before l <- combinedDFInteractiveHeatMapP_Val(combinedDFP_Val_Labels, dendProbes, dendTraits, selectedRowIndices, selectedColIndices)"))
#'       #browser() #check, whether this is called twice
#'       l <-
#'         combinedDFInteractiveHeatMapP_Val(combinedDFP_Val_Labels, dendProbes, dendTraits, Distances, selectedRowIndicesYellow, selectedColIndices, selectedRowIndicesOrange, session)
#'       base::print(base::paste0(sysTimePID(), " before combinedHM <- l$combinedHM"))
#'       combinedHM <- l$combinedHM
#'
#'       endTime <- Sys.time()
#'       elapsedTime <- endTime - startTime
#'       base::print(base::paste0(sysTimePID(), " after calculating condensed DW heatmap. Elapsed time: ", elapsedTime, " sec."))
#'       base::print(base::paste0(sysTimePID(), " before plotting condensed DW heatmap; InteractiveComplexHeatmap::makeInteractiveComplexHeatmap()"))
#'       while (!base::is.null(grDevices::dev.list())) {
#'         grDevices::dev.off()
#'       }
#'       startTime <- Sys.time()
#'       InteractiveComplexHeatmap::makeInteractiveComplexHeatmap(
#'         input = input,
#'         output = output,
#'         session = session,
#'         ht_list = combinedHM,
#'         heatmap_id = "condDWHeatmap_P_Val",
#'         show_layer_fun = TRUE
#'       )
#'     }
#'   },
#'   error = function(e) {
#'     base::message("An error occurred in plotCombinedDWHM_P_Val():\n", e)
#'     Cstack_info()
#'     browser() #should not happen
#'   },
#'   warning = function(w) {
#'     base::message("A warning occurred in plotCombinedDWHM_P_Val():\n", w)
#'     browser() #should not happen
#'   },
#'   finally = {
#'     endTime <- Sys.time()
#'     elapsedTime <- endTime - startTime
#'     base::print(base::paste0(sysTimePID(), " after plotting heatmap for DWP_Val; InteractiveComplexHeatmap::makeInteractiveComplexHeatmap(). Elapsed time: ", elapsedTime, " sec."))
#'     output$txtCondDWHMDescription_P_Val <-
#'       shiny::renderText(
#'         base::paste0(
#'           sysTimePID(),
#'           " done plotting heatmap for condensed DWP_Val..., current plot is valid. n(probe) = ",
#'           base::nrow(base::as.matrix(combinedDFP_Val_Labels[[1]])),
#'           "; n(trait) = ",
#'           base::ncol(base::as.matrix(combinedDFP_Val_Labels[[1]])),
#'           "; elapsed time: ",
#'           elapsedTime, " sec."
#'         )
#'       )
#'   })
#' }

#' click_action_fullHM_P_Val
#' @param df data.frame containing data which is enclosed by brush action
#' @param input shiny input object
#' @param output shiny output object
#' @param session shiny session object
#' @return nothing; function info label for click action in HM
# examples click_action_fullHM_P_Val(df, input, output, session)
click_action_fullHM_P_Val <- function(df, input, output, session) {
  base::tryCatch(
    {
      base::print(base::paste0(sysTimePID(), " start click_action_fullHM_P_Val()."))
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
      base::message("An error occurred in click_action_fullHM_P_Val():\n", e)
    },
    warning = function(w) {
      base::message("A warning occurred in click_action_fullHM_P_Val():\n", w)
    },
    finally = {
      base::print(base::paste0(sysTimePID(), " end click_action_fullHM_P_Val()."))
    }
  )
}

#' brush_action_fullHM_P_Val
#' @param df data.frame containing data which is enclosed by brush action
#' @param input shiny input object
#' @param output shiny output object
#' @param session shiny session object
#' @return nothing; function creates SPLOM
# examples brush_action_fullHM_P_Val(df, input, output, session)
brush_action_fullHM_P_Val <- function(df, input, output, session) {
  base::tryCatch(
    if (!is.null(df)) {
      row_index <- collapse::funique(unlist(df$row_index))
      column_index <- collapse::funique(unlist(df$column_index))
      #feed in selected CpG here
      session$userData$sessionVariables$selectedProbe(rownames(session$userData$sessionVariables$traitReducedDataStructurePVal()$combinedDFP_Val_Labels$dfP_Val)[row_index])
      originTrait <- session$userData$sessionVariables$traitReducedDataStructurePVal()$combinedDFP_Val_Labels$mergedOriginTrait[column_index]
      traitLabels <- session$userData$sessionVariables$traitReducedDataStructurePVal()$combinedDFP_Val_Labels$mergedOriginalColnames[column_index]
      selectedTrait <- cbind(traitLabels, originTrait)
      colnames(selectedTrait) <- c("traitName", "traitSource")
      session$userData$sessionVariables$selectedTrait(selectedTrait)
      selectedTraitID <- session$userData$sessionVariables$traitReducedDataStructurePVal()$combinedDFP_Val_Labels$traitID[column_index]
      session$userData$sessionVariables$selectedTraitID(selectedTraitID)
      base::print(base::paste0(sysTimePID(), " last selection was from full heatmap p-val brush event."))
      output$txtGlobalSelectOut <- shiny::renderText("last selection was from full heatmap p-val brush event.")
    },
    error = function(e) {
      base::message("An error occurred in brush_action_fullHM_P_Val():\n", e)
    },
    warning = function(w) {
      base::message("A warning occurred in brush_action_fullHM_P_Val():\n", w)
    },
    finally = {
      base::print(base::paste0(sysTimePID(), " end brush_action_fullHM_P_Val()."))
    }
  )
}

#' hover_action_fullHM_P_Val
#' @param df data.frame containing data which is enclosed by brush action
#' @param input shiny input object
#' @param output shiny output object
#' @param session shiny session object
#' @return nothing; function info label for click action in HM
# examples click_action_fullHM_P_Val(df, input, output, session)
hover_action_fullHM_P_Val <- function(df, input, output, session) {
  base::tryCatch(
    {
      base::print(base::paste0(sysTimePID(), " hover_action_fullHM_P_Val.", as.character(head(df))))
    },
    error = function(e) {
      base::message("An error occurred in hover_action_fullHM_P_Val():\n", e)
    },
    warning = function(w) {
      base::message("A warning occurred in hover_action_fullHM_P_Val():\n", w)
    },
    finally = {
      base::print(base::paste0(sysTimePID(), " end hover_action_fullHM_P_Val()."))
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
      base::message("An error occurred in click_action_fullHM_LogFC():\n", e)
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
      session$userData$sessionVariables$selectedProbe(rownames(session$userData$sessionVariables$traitReducedDataStructureLogFC()$combinedDFP_Val_Labels$dfLogFC)[row_index])
      originTrait <- session$userData$sessionVariables$traitReducedDataStructureLogFC()$combinedDFP_Val_Labels$mergedOriginTrait[column_index]
      traitLabels <- session$userData$sessionVariables$traitReducedDataStructureLogFC()$combinedDFP_Val_Labels$mergedOriginalColnames[column_index]
      selectedTrait <- cbind(traitLabels, originTrait)
      colnames(selectedTrait) <- c("traitName", "traitSource")
      session$userData$sessionVariables$selectedTrait(selectedTrait)
      selectedTraitID <- session$userData$sessionVariables$traitReducedDataStructureLogFC()$combinedDFP_Val_Labels$traitID[column_index]
      session$userData$sessionVariables$selectedTraitID(selectedTraitID)
      base::print(base::paste0(sysTimePID(), " last selection was from full heatmap log(FC) brush event."))
      output$txtGlobalSelectOut <- shiny::renderText("last selection was from full heatmap log(FC) brush event.")
    },
    error = function(e) {
      base::message("An error occurred in brush_action_fullHM_LogFC():\n", e)
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
      base::message("An error occurred in hover_action_fullHM_LogFC():\n", e)
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
      base::message("An error occurred in click_action_condHM_P_Val():\n", e)
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
browser() #tbc()
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
      base::message("An error occurred in brush_action_condHM_P_Val():\n", e)
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
      base::message("An error occurred in hover_action_condHM_P_Val():\n", e)
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
      base::message("An error occurred in click_action_condHeatmap_LogFC():\n", e)
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
browser() #tbc()
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
      base::message("An error occurred in brush_action_condHeatmap_LogFC():\n", e)
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
      base::message("An error occurred in hover_action_condHeatmap_LogFC():\n", e)
    },
    warning = function(w) {
      base::message("A warning occurred in hover_action_condHeatmap_LogFC():\n", w)
    },
    finally = {
      base::print(base::paste0(sysTimePID(), " end hover_action_condHeatmap_LogFC()."))
    }
  )
}
