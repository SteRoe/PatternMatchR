#' structure of sessionVariables:

#' result data for/from clustering
#' session$userData$sessionVariables$distMatProbes - distance matrix for probes out of p-value results of regression
#' session$userData$sessionVariables$clustResProbes - clustering reusult for probes out of p-value results of regression
#' session$userData$sessionVariables$distMatTraits - distance matrix for traits out of p-value results of regression
#' session$userData$sessionVariables$clustResTraits - clustering reusult for traits out of p-value results of regression
#' session$userData$sessionVariables$combinedHMP_Val - data structure for combined heatmap

#' getSelectedOriginalData
#' @param combinedDFP_Val_Labels list with datastructure pointing to original data
#' @param row_index indicies of selection
#' @param column_index indicies of selection
#' @param session shiny session object
#' @return df with merged original data
# examples getSelectedOriginalData(combinedDFP_Val_Labels, row_index, column_index)
getSelectedOriginalData <- function(combinedDFP_Val_Labels, row_index, column_index, session) {
  base::print(base::paste0(Sys.time(), " start getSelectedOriginalData()"))
  tryCatch(
    {
      selectedColnames <- combinedDFP_Val_Labels$mergedColnames[column_index] # we here have the colnames of the selected traits
      selectedColnames <- removeAdjFromColname(selectedColnames)
      selectedTraitSources <- combinedDFP_Val_Labels$mergedOriginTrait[column_index] # we here have the trait indicies of the selected traits
      # selectedTraitOriginDF <- combinedDFP_Val_Labels$mergedOriginDF[column_index] #we here have the DF indicies of the selected traits
      selectedColnamesTrait1 <- selectedColnames[selectedTraitSources == 1]
      selectedColnamesTrait2 <- selectedColnames[selectedTraitSources == 2]
      selectedColnamesTrait3 <- selectedColnames[selectedTraitSources == 3]
      #to be sure we select only colnames, which are within PHENODF:
      selectedColnamesTrait1 <- intersect(colnames(session$userData$sessionVariables$resultDFListTrait1()$listPHENOdata[[1]]$PHENODF), selectedColnamesTrait1)
      selectedColnamesTrait2 <- intersect(colnames(session$userData$sessionVariables$resultDFListTrait2()$listPHENOdata[[1]]$PHENODF), selectedColnamesTrait2)
      selectedColnamesTrait3 <- intersect(colnames(session$userData$sessionVariables$resultDFListTrait3()$listPHENOdata[[1]]$PHENODF), selectedColnamesTrait3)
      # get selected original data from trait data
      selectedDFTrait1 <- session$userData$sessionVariables$resultDFListTrait1()$listPHENOdata[[1]]$PHENODF[selectedColnamesTrait1]
      selectedDFTrait2 <- session$userData$sessionVariables$resultDFListTrait2()$listPHENOdata[[1]]$PHENODF[selectedColnamesTrait2]
      selectedDFTrait3 <- session$userData$sessionVariables$resultDFListTrait3()$listPHENOdata[[1]]$PHENODF[selectedColnamesTrait3]
      # merge all trait data together by Kind_ID or rowname (better)
      rn <- rownames(selectedDFTrait1)
      selectedDFTrait1$Row.names <- rn
      rn <- rownames(selectedDFTrait2)
      selectedDFTrait2$Row.names <- rn
      rn <- rownames(selectedDFTrait3)
      selectedDFTrait3$Row.names <- rn
      selectedDF <- NULL
      if (!base::is.null(selectedDFTrait1) && !base::is.null(selectedDFTrait2)) {
        selectedDF <- merge(selectedDFTrait1, selectedDFTrait2, by = "Row.names", all.x = FALSE, all.y = FALSE)
      }
      else {
        if (!base::is.null(selectedDFTrait2)) {
          selectedDF <- selectedDFTrait2
        }
      }
      if (!base::is.null(selectedDF) && !base::is.null(selectedDFTrait3)) {
        selectedDF <- merge(selectedDF, selectedDFTrait3, by = "Row.names", all.x = FALSE, all.y = FALSE)
      }
      else {
        if (!base::is.null(selectedDFTrait3)) {
          selectedDF <- selectedDFTrait3
        }
      }
      rownames(selectedDF) <- selectedDF$Row.names

      # get selected methylation data...
      selectedRownames <- rownames(combinedDFP_Val_Labels$dfP_Val)[row_index]
      #subset selectedRownames to only keep those, that are loaded in $Beta_tDF
      selectedRownames <- intersect(colnames(session$userData$Beta_tDF), selectedRownames)
      selectedBeta <- session$userData$Beta_tDF[, selectedRownames] #if error
      #"nicht definierte Spalten gewählt" occurs, this is due to debug mode,
      #where most columns in Beta_tDF are not loaded.
      #... and merge with already merged trait data
      rn <- rownames(selectedBeta)
      if (is.valid(rn)) {
        selectedBeta$Row.names <- rn
        #selectedDF$Row.names <- rownames(selectedDF)
        selectedDF_Beta <- merge(selectedDF, selectedBeta, by = "Row.names", all.x = FALSE, all.y = FALSE)
        rownames(selectedDF_Beta) <- selectedDF_Beta$Row.names
        selectedDF_Beta$Row.names <- NULL
      }
      else {
        message("We miss rownames in selectedDF_Beta here... (in getSelectedOriginalData()).\n
                Reason might be, that the beta data set was not loaded in full length (debugMode == TRUE?).\n")
      }
      if (nrow(selectedDF_Beta) > 256 || ncol(selectedDF_Beta) > 256) {
        base::message(base::paste0(Sys.time(), "nrow(selectedDF) = ",
                                   nrow(selectedDF_Beta),
                                   "&& ncol(selectedDF) = ",
                                   ncol(selectedDF_Beta),
                                   " that might be too much for fast processing!"))
      }
      # result gives 3d-data structure: multiple methylation profiles!!!
      # create SPLOM... each variable a dimension
      # data structure needed by plotly is described here: https://plotly.com/r/splom/

      # get involved original data sources
      # dataSources <- combinedDFP_Val_Labels$

      # iterate over those data sources
      # append data to result data frame

      # return df with merged original data from selected area
    },
    error = function(e) {
      message("An error occurred in getSelectedOriginalData():\n", e)
    },
    warning = function(w) {
      message("A warning occurred in getSelectedOriginalData():\n", w)
    },
    finally = {
      base::print(base::paste0(Sys.time(), " end getSelectedOriginalData()"))
      return(selectedDF_Beta)
    }
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
}

# emptyHM <- compiler::cmpfun(emptyHM)

#' creates a regular heatmap
#' @param combinedDF_Labels list of data.frame and labels generated from function mergeDFP_Val_Labels()
#' @param dendProbes dendrogram (without labels) for probes (rows), generated externally and providing information for clustering of heatmap
#' @param dendTraits dendrogram (without labels) result for traits (columns), generated externally and providing information for clustering of heatmap
#' @return heatmap object for InteractiveComplexHeatmap::makeInteractiveComplexHeatmap
# examples combinedDFInteractiveHeatMapP_Val(combinedDF_Labels, dendProbes, dendTraits)
combinedDFInteractiveHeatMapP_Val <-
  function(combinedDF_Labels,
           #           dendProbes = NA,
           #          dendTraits = NA, # {
           dendProbes = NA,
           dendTraits = NA) {
    tryCatch(
      {
        #browser() #check for usage of already merged data
        base::print(base::paste0(Sys.time(), " start making HM"))
        mat <- base::as.matrix(combinedDF_Labels$dfP_Val)
        matDM <- base::as.matrix(combinedDF_Labels$dfDM)
        matN <- base::as.matrix(combinedDF_Labels$dfN)
        # ((((((((((((((((((((((((()))))))))))))))))))))))))
        # reduce matrix dimensions before plotting like in
        # https://gdevailly.netlify.app/post/plotting-big-matrices-in-r/
        # ((((((((((((((((((((((((()))))))))))))))))))))))))
        # and use rasterization like described in
        # https://jokergoo.github.io/2020/06/30/rasterization-in-complexheatmap/
        base::print(base::paste0(Sys.time(), " making labels"))
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
        base::print(base::paste0(Sys.time(), " making annotations"))
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

        base::print(base::paste0(Sys.time(), " making colors"))
        max.Col1 <- 0.05
        min.Col1 <- base::min(mat, na.rm = TRUE)
        max.Col2 <- 0.05
        min.Col2 <- base::min(mat, na.rm = TRUE)
        max.Col3 <- 0.05
        min.Col3 <- base::min(mat, na.rm = TRUE)
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
        base::print(base::paste0(Sys.time(), " making heatmap n(probes)=", base::dim(mat)[1], " x n(traits)=", base::dim(mat)[2]))
        if ("dendrogram" %in% class(dendProbes) && "dendrogram" %in% class(dendTraits)) {
          base::print(base::paste0(Sys.time(), " length(unlist(dendProbes)): ", length(unlist(dendProbes))))
          base::print(base::paste0(Sys.time(), " length(unlist(dendTraits)): ", length(unlist(dendTraits))))
          length(unlist(dendTraits)) == base::dim(mat)[2]
          length(unlist(dendProbes)) == base::dim(mat)[1]
#browser()
          ht <-
            ComplexHeatmap::Heatmap(
              mat,
              rect_gp = grid::gpar(type = "none"),
              cluster_rows = dendProbes,
              cluster_columns = dendTraits,
              top_annotation = ha,
              layer_fun = function(j, i, x, y, w, h, fill) {
                if (length(i) == base::nrow(mat) * base::ncol(mat)) {
                  # we are in main HM
                  subHM <- FALSE
                } else {
                  # we are in sub HM
                  subHM <- TRUE
                }
                l <- labels[j] == "trait 1"
                if (any(l)) {
                  grid::grid.rect(x[l], y[l], w[l], h[l], gp = grid::gpar(
                    fill = col1(ComplexHeatmap::pindex(mat, i[l], j[l])),
                    col = col1(ComplexHeatmap::pindex(mat, i[l], j[l]))
                  ))
                }
                l <- labels[j] == "trait 2"
                if (any(l)) {
                  grid::grid.rect(x[l], y[l], w[l], h[l], gp = grid::gpar(
                    fill = col2(ComplexHeatmap::pindex(mat, i[l], j[l])),
                    col = col2(ComplexHeatmap::pindex(mat, i[l], j[l]))
                  ))
                }
                l <- labels[j] == "trait 3"
                if (any(l)) {
                  grid::grid.rect(x[l], y[l], w[l], h[l], gp = grid::gpar(
                    fill = col3(ComplexHeatmap::pindex(mat, i[l], j[l])),
                    col = col3(ComplexHeatmap::pindex(mat, i[l], j[l]))
                  ))
                }
                if (subHM == TRUE) {
                  grid::grid.text(
                    paste0(
                      "p:",
                      sprintf("%.G", ComplexHeatmap::pindex(mat, i, j)),
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
          base::print(base::paste0(Sys.time(), " at least one distance matrix is not of class \"dendrogram\""))
        }
        base::print(base::paste0(Sys.time(), " making legend"))
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
        message("An error occurred in combinedDFInteractiveHeatMapP_Val():\n", e)
      },
      warning = function(w) {
        message("A warning occurred in combinedDFInteractiveHeatMapP_Val():\n", w)
      },
      finally = {
        base::print(base::paste0(Sys.time(), " end preparing heatmap."))
      }
    )

    tryCatch(
      {
        base::print(base::paste0(Sys.time(), " clearing grDevices()"))
        if (grDevices::dev.cur() > 1) {
          grDevices::dev.off()
        }
        grDevices::pdf(NULL)
      },
      error = function(e) {
        message("An error occurred in combinedDFInteractiveHeatMapP_Val():\n", e)
      },
      warning = function(w) {
        message("A warning occurred in combinedDFInteractiveHeatMapP_Val():\n", w)
      },
      finally = {
        base::print(base::paste0(Sys.time(), " end clearing grDevices()."))
      }
    )

    tryCatch(
      {
        base::print(base::paste0(Sys.time(), " plotting heatmap (takes some time)"))
        # with huge heatmaps, the following error occurs:
        # Error in Cairo: Failed to create Cairo backend!
        ht <- ComplexHeatmap::draw(ht, annotation_legend_list = lgd)
      },
      error = function(err) {
        base::print(base::paste0(Sys.time(), " Error: unable to plot HM. ", err$message))
      }
    )

    grDevices::dev.off()
    l <- base::list()
    l$combinedHMP_VAL <- ht
    base::print(base::paste0(Sys.time(), " end plotting heatmap"))
    return(l)
  }

#' click_action_HM
#' @param df data.frame containing data which is enclosed by brush action
#' @param input shiny input object
#' @param output shiny output object
#' @param session shiny session object
#' @return nothing; function info label for click action in HM
# examples click_action_HM(df, input, output, session)
click_action_HM <- function(df, input, output, session) {
  # browser()
  output[["info1"]] <- shiny::renderUI({
    if (!is.null(df)) {
      htmltools::HTML(
        GetoptLong::qq(
          "<p style='background-color:#FF8080;color:white;padding:5px;'>
          row_label: @{df$row_label}, col_label: @{df$column_label},
          row: @{df$row_index}, column: @{df$column_index}</p>"
        )
      )
    }
  })
}

#' brush_action_HM
#' @param df data.frame containing data which is enclosed by brush action
#' @param input shiny input object
#' @param output shiny output object
#' @param session shiny session object
#' @return nothing; function creates SPLOM
# examples brush_action_HM(df, input, output, session)
brush_action_HM <- function(df, input, output, session) {
  tryCatch(
    if (!is.null(df)) {
      session$userData$sessionVariables$SPLOM <- FALSE
      row_index <- collapse::funique(unlist(df$row_index)) #row_index <- unique(unlist(df$row_index))
      column_index <- collapse::funique(unlist(df$column_index))
      # session$userData$sessionVariables$selected_row_index <-
      #   row_index
      # session$userData$sessionVariables$selected_column_index <-
      #   column_index
      session$userData$sessionVariables$selectedOriginalData <-
        getSelectedOriginalData(session$userData$sessionVariables$pReducedcombinedDFP_Val_Labels(), row_index, column_index, session)
      if (!is.null(session$userData$sessionVariables$selectedOriginalData)) {
        FactorialVars <- getBinaryFactorialVars(session$userData$sessionVariables$selectedOriginalData)
        if (!is.null(FactorialVars)) {
          shiny::updateSelectizeInput(
            session = session,
            inputId = "markingVar",
            choices = FactorialVars,
            server = TRUE
          )
        } else {
          shiny::updateSelectizeInput(
            session = session,
            inputId = "markingVar",
            choices = NULL,
            server = TRUE
          )
        }
      } else {
        shiny::updateSelectizeInput(
          session = session,
          inputId = "markingVar",
          choices = NULL,
          server = TRUE
        )
      }
      if (session$userData$sessionVariables$SPLOM == FALSE) { # SPLOM empty
        fig <- getSPLOM(session$userData$sessionVariables$selectedOriginalData, session$userData$sessionVariables$markingVar)
        output$SPLOM <- plotly::renderPlotly(fig)
        session$userData$sessionVariables$SPLOM <- TRUE
      }
    },
    error = function(e) {
      message("An error occurred in brush_action_HM():\n", e)
    },
    warning = function(w) {
      message("A warning occurred in brush_action_HM():\n", w)
    },
    finally = {
    }
  )
}
