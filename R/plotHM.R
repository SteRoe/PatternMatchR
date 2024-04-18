getSelectedTraitReducedcombinedDFP_Val_Labels <- function(combinedDFP_Val_Labels, row_index, column_index, session) {
  base::print(base::paste0(sysTimePID(), " start getSelectedTraitReducedcombinedDFP_Val_Labels()"))
  base::tryCatch(
    {
      colInd <- which(combinedDFP_Val_Labels$mergedColnames %in% column_index)
      rowInd <- which(combinedDFP_Val_Labels$mergedRownames %in% row_index)
      #selectedColnames <- combinedDFP_Val_Labels$mergedColnames[colInd] # we here have the colnames of the selected traits
      selectedColnames <- combinedDFP_Val_Labels$mergedOriginalColnames[colInd] # we here have the colnames of the selected traits
      selectedColnames <- removeAdjFromColname(selectedColnames)
      result <- list()
      result$dfP_Val <- combinedDFP_Val_Labels$dfP_Val[row_index,]
      result$dfDM <- combinedDFP_Val_Labels$dfDM[row_index,]
      result$dfN <- combinedDFP_Val_Labels$dfN[row_index,]
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

#' getSelectedOriginalData
#' @param combinedDFP_Val_Labels list with datastructure pointing to original data
#' @param row_index indicies of selection
#' @param column_index indicies of selection
#' @param session shiny session object
#' @return df with merged original data
# examples getSelectedOriginalData(combinedDFP_Val_Labels, row_index, column_index)
getSelectedOriginalData <- function(combinedDFP_Val_Labels, row_index, column_index, session) {
  base::print(base::paste0(sysTimePID(), " start getSelectedOriginalData()"))
  base::tryCatch(
    {
      colInd <- which(combinedDFP_Val_Labels$mergedColnames %in% column_index)
      #selectedColnames <- combinedDFP_Val_Labels$mergedColnames[colInd] # we here have the colnames of the selected traits
      selectedColnames <- combinedDFP_Val_Labels$mergedOriginalColnames[colInd]
      selectedColnames <- removeAdjFromColname(selectedColnames)
      selectedTraitSources <- combinedDFP_Val_Labels$mergedOriginTrait[colInd] # we here have the trait indicies of the selected traits
      selectedColnamesTrait1 <- selectedColnames[selectedTraitSources == 1]
      selectedColnamesTrait2 <- selectedColnames[selectedTraitSources == 2]
      selectedColnamesTrait3 <- selectedColnames[selectedTraitSources == 3]
      #to be sure we select only colnames, which are within PHENODF:
      selectedColnamesTrait1 <- intersect(colnames(session$userData$sessionVariables$resultDFListTrait1()$listPHENOdata[[1]]$PHENODF), selectedColnamesTrait1)
      selectedColnamesTrait2 <- intersect(colnames(session$userData$sessionVariables$resultDFListTrait2()$listPHENOdata[[1]]$PHENODF), selectedColnamesTrait2)
      selectedColnamesTrait3 <- intersect(colnames(session$userData$sessionVariables$resultDFListTrait3()$listPHENOdata[[1]]$PHENODF), selectedColnamesTrait3)
      if (!is.valid(selectedColnamesTrait1)) {
        base::message(base::paste0(sysTimePID(), " file names in trait 1 folder do not match colnames in pheno file! SPLOM will not work."))
      }
      if (!is.valid(selectedColnamesTrait2)) {
        base::message(base::paste0(sysTimePID(), " file names in trait 2 folder do not match colnames in pheno file! SPLOM will not work."))
      }
      if (!is.valid(selectedColnamesTrait3)) {
        base::message(base::paste0(sysTimePID(), " file names in trait 3 folder do not match colnames in pheno file! SPLOM will not work."))
      }
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
      rowInd <- which(rownames(combinedDFP_Val_Labels$dfP_Val) %in% row_index)
      selectedRownames <- rownames(combinedDFP_Val_Labels$dfP_Val)[rowInd]
      #subset selectedRownames to only keep those, that are loaded in $Beta_tDF
      selectedRownames <- intersect(colnames(session$userData$Beta_tDF), selectedRownames)
      selectedBeta <- session$userData$Beta_tDF[, selectedRownames] #if error
      #"nicht definierte Spalten gewählt" occurs, this is due to debug mode,
      #where most columns in Beta_tDF are not loaded.
      #... and merge with already merged trait data
      rn <- rownames(selectedBeta)
      if (is.valid(rn)) {
        selectedBeta$Row.names <- rn
        selectedDF_Beta <- merge(selectedDF, selectedBeta, by = "Row.names", all.x = FALSE, all.y = FALSE)
        #selectedDF_Beta <- merge(selectedBeta, selectedDF, by = "Row.names", all.x = FALSE, all.y = FALSE) #beta first, then traits
        rownames(selectedDF_Beta) <- selectedDF_Beta$Row.names
        selectedDF_Beta$Row.names <- NULL
      }
      else {
        message("We miss rownames in selectedDF_Beta here... (in getSelectedOriginalData()).\n
                Reason might be, that the beta data set was not loaded in full length (debugMode == TRUE?).\n")
      }
      if (nrow(selectedDF_Beta) > 256 || ncol(selectedDF_Beta) > 256) {
        base::message(base::paste0(sysTimePID(), " Warning: nrow(selectedDF) = ",
                                   nrow(selectedDF_Beta),
                                   " || ncol(selectedDF) = ",
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
      base::message("An error occurred in getSelectedOriginalData():\n", e)
    },
    warning = function(w) {
      base::message("A warning occurred in getSelectedOriginalData():\n", w)
    },
    finally = {
      base::print(base::paste0(sysTimePID(), " end getSelectedOriginalData()"))
      return(selectedDF_Beta)
    }
  )
}

getSelectedOriginalDataTraits <- function(combinedDFP_Val_Labels, row_index, column_index, session) {
  base::print(base::paste0(sysTimePID(), " start getSelectedOriginalDataTraits()"))
  base::tryCatch(
    {
      colInd <- which(combinedDFP_Val_Labels$mergedColnames %in% column_index)
      #selectedColnames <- combinedDFP_Val_Labels$mergedColnames[colInd] # we here have the colnames of the selected traits
      selectedColnames <- combinedDFP_Val_Labels$mergedOriginalColnames[colInd] # we here have the colnames of the selected traits
      selectedColnames <- removeAdjFromColname(selectedColnames)
      selectedTraitSources <- combinedDFP_Val_Labels$mergedOriginTrait[colInd] # we here have the trait indicies of the selected traits
      selectedColnamesTrait1 <- selectedColnames[selectedTraitSources == 1]
      selectedColnamesTrait2 <- selectedColnames[selectedTraitSources == 2]
      selectedColnamesTrait3 <- selectedColnames[selectedTraitSources == 3]
      #to be sure we select only colnames, which are within PHENODF:
      selectedColnamesTrait1 <- intersect(colnames(session$userData$sessionVariables$resultDFListTrait1()$listPHENOdata[[1]]$PHENODF), selectedColnamesTrait1)
      selectedColnamesTrait2 <- intersect(colnames(session$userData$sessionVariables$resultDFListTrait2()$listPHENOdata[[1]]$PHENODF), selectedColnamesTrait2)
      selectedColnamesTrait3 <- intersect(colnames(session$userData$sessionVariables$resultDFListTrait3()$listPHENOdata[[1]]$PHENODF), selectedColnamesTrait3)
      if (!is.valid(selectedColnamesTrait1)) {
        base::message(base::paste0(sysTimePID(), " file names in trait 1 folder do not match colnames in pheno file! SPLOM will not work."))
      }
      if (!is.valid(selectedColnamesTrait2)) {
        base::message(base::paste0(sysTimePID(), " file names in trait 2 folder do not match colnames in pheno file! SPLOM will not work."))
      }
      if (!is.valid(selectedColnamesTrait3)) {
        base::message(base::paste0(sysTimePID(), " file names in trait 3 folder do not match colnames in pheno file! SPLOM will not work."))
      }
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

      if (length(selectedColnames) > 256) {
        base::message(base::paste0(sysTimePID(), "length(selectedColnames) = ",
                                   length(selectedColnames),
                                   " that might be too much for fast processing!"))
      }
      #remove row.names from selectedDF
      if("Row.names" %in% colnames(selectedDF)) {
        selectedDF$Row.names <- NULL
      }
    },
    error = function(e) {
      base::message("An error occurred in getSelectedOriginalDataTraits():\n", e)
    },
    warning = function(w) {
      base::message("A warning occurred in getSelectedOriginalDataTraits():\n", w)
    },
    finally = {
      base::print(base::paste0(sysTimePID(), " end getSelectedOriginalDataTraits()"))
      return(selectedDF)
    }
  )
}

getSelectedOriginalDataProbes <- function(combinedDFP_Val_Labels, row_index, column_index, session) {
  base::print(base::paste0(sysTimePID(), " start getSelectedOriginalDataProbes()"))
  base::tryCatch(
    {
      # get selected methylation data...
      rowInd <- which(rownames(combinedDFP_Val_Labels$dfP_Val) %in% row_index)
      selectedRownames <- rownames(combinedDFP_Val_Labels$dfP_Val)[rowInd]
#      selectedRownames <- rownames(combinedDFP_Val_Labels$dfP_Val)[row_index] #which(column_index %in% combinedDFP_Val_Labels$mergedColnames)
      #subset selectedRownames to only keep those, that are loaded in $Beta_tDF
      selectedRownames <- intersect(colnames(session$userData$Beta_tDF), selectedRownames)
      selectedBeta <- session$userData$Beta_tDF[, selectedRownames] #if error
      #"nicht definierte Spalten gewählt" occurs, this is due to debug mode,
      #where most columns in Beta_tDF are not loaded.
      #... and merge with trait data from markingVar
      #get Values of markingVar from traits table associated with ID#
#browser()
      markingVar <- session$userData$sessionVariables$markingVar()
      traits <- session$userData$sessionVariables$selectedOriginalDataTraits()
      #select only ID# and markingVar
      traits$id <- rownames(traits)
      Vars <- c("id", markingVar)
      if (all(Vars %in% colnames(traits))) {
        traits <- traits[,Vars]
        traits$markingVar <- traits[,markingVar]
        traits[,markingVar] <- NULL
        #merge
        selectedBeta$id <- rownames(selectedBeta)
        selectedBeta <- merge(selectedBeta, traits, by.x = "id", by.y = "id")
        rownames(selectedBeta) <- selectedBeta$id
        selectedBeta$id <- NULL
      }
      rn <- rownames(selectedBeta)
      if (!is.valid(rn)) {
        message("We miss rownames in selectedDF_Beta here... (in getSelectedOriginalData()).\n
                Reason might be, that the beta data set was not loaded in full length (debugMode == TRUE?).\n")
      }
      if (length(selectedRownames) > 256) {
        base::message(base::paste0(sysTimePID(), "length(selectedRownames) = ",
                                   length(selectedRownames),
                                   " that might be too much for fast processing!"))
      }
      #remove row.names from selectedDF
      if("Row.names" %in% colnames(selectedBeta)) {
        selectedBeta$Row.names <- NULL
      }
    },
    error = function(e) {
      base::message("An error occurred in getSelectedOriginalDataProbes():\n", e)
    },
    warning = function(w) {
      base::message("A warning occurred in getSelectedOriginalDataProbes():\n", w)
    },
    finally = {
      base::print(base::paste0(sysTimePID(), " end getSelectedOriginalDataProbes()"))
      return(selectedBeta)
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
# examples combinedDFInteractiveHeatMapP_Val(combinedDF_Labels, dendProbes, dendTraits, selectedRowIndicesYellow, selectedColIndices, selectedRowIndicesOrange)
combinedDFInteractiveHeatMapP_Val <-
  function(combinedDF_Labels,
           dendProbes = NA,
           dendTraits = NA,
           selectedRowIndicesYellow = NA,
           selectedColIndices = NA,
           selectedRowIndicesOrange = NA,
           session) {
    #function is not called twice (check caller function), if so, check here
    base::tryCatch(
      {
        startTime <- Sys.time()
        base::print(base::paste0(sysTimePID(), " start preparing HM; combinedDFInteractiveHeatMapP_Val()"))
        matP_Val <- base::as.matrix(combinedDF_Labels$dfP_Val)
        matDM <- base::as.matrix(combinedDF_Labels$dfDM)
        matN <- base::as.matrix(combinedDF_Labels$dfN)
        # ((((((((((((((((((((((((()))))))))))))))))))))))))
        # reduce matrix dimensions before plotting like in
        # https://gdevailly.netlify.app/post/plotting-big-matrices-in-r/
        # ((((((((((((((((((((((((()))))))))))))))))))))))))
        # and use rasterization like described in
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
    #add two right second heatmap for new measures here:
    base::tryCatch(
      {
        #numrows <- base::nrow(matP_Val)
        #vec <- stats::rnorm(numrows)
        #tbc()
        #tbc()
#browser()
        vec <- session$userData$sessionVariables$distNeigboursProbes100()
        vec <- vec[,2] #take means
        #normalize vec to -1...1
        hist(vec, plot = TRUE)
        vec <-  scales::rescale(vec, to = c(-1, 1))
        #invert vec
        vec <- vec * -1
        ht2 <- ComplexHeatmap::Heatmap(vec, name = "distances on DNA (100); low value: short distance.")
        #vec <- stats::rnorm(numrows)
        #ht3 <- ComplexHeatmap::Heatmap(vec, name = "new measure2")
        vec <- session$userData$sessionVariables$distNeigboursProbes1000()
        vec <- vec[,2] #take means
        hist(vec, plot = TRUE)
        #normalize vec to -1...1
        vec <-  scales::rescale(vec, to = c(-1, 1))
        #invert vec
        vec <- vec * -1
        ht3 <- ComplexHeatmap::Heatmap(vec, name = "distances on DNA (1000); low value: short distance.")
        #fill nA with 1, so distances of far away CpG dont become modified
        #other variant: fill with 0, so those become zeroes
        #other variant: leave NA, so those become omitted
        vec <- data.table::nafill(vec, type = "const", fill = 1)
        library(data.table)
browser() #check dimensions of dt and vec!!!
        dt <- data.table(matP_Val)
        DWP_Val <- data.table(t(t(dt)*vec))
        # for (i in seq_along(dt[,2])) {
        #   dt[, (i) := dt[[i]] * vec]
        # }
        # DWP_Val <- dt
#browser()
#vec has Distances of neighbouring CpG dist = 1000 up/down and becomes multiplied with p_val
        # we here try to find trait pairs, that differ in trait colour (red, green, blue), but share high similarity
        DWP_Val <- as.matrix(DWP_Val)
        DWht <-
          ComplexHeatmap::Heatmap(
            DWP_Val,
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
      },
      error = function(e) {
        base::message("An error occurred in add new measure:\n", e)
      },
      warning = function(w) {
        base::message("A warning occurred in add new measure:\n", w)
      },
      finally = {
        base::print(base::paste0(sysTimePID(), " end add new measure."))
      }
    )
    base::tryCatch(
      {
        base::print(base::paste0(sysTimePID(), " start drawing heatmap (takes some time). (step before ComplexHeatmap::draw()"))
        # with huge heatmaps, the following error occurs:
        # Error in Cairo: Failed to create Cairo backend!
        ht <- ComplexHeatmap::draw(ht + ht2 + ht3, annotation_legend_list = lgd)
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
        base::print(base::paste0(sysTimePID(), " l$combinedHMP_VAL <- ht"))
        l$combinedHMP_VAL <- ht
        l$DWcombinedHMP_VAL <- DWht
        endTime <- Sys.time()
        elapsedTime <- endTime - startTime
        base::print(base::paste0(sysTimePID(), " end combinedDFInteractiveHeatMapP_Val()"))
        base::print(base::paste0(sysTimePID(), " finished drawing heatmap (takes some time). (step after ComplexHeatmap::draw(); Elapsed time: ", elapsedTime, "."))
        base::return(l)
      }
    )
  }

getSearchResultCpG <- function(txtSearchCpG, session) {
  base::tryCatch(
    {
      base::print(base::paste0(sysTimePID(), " searching CpG"))
      #look into clustResProbes and find position of CpG
      CpG <- session$userData$sessionVariables$traitReducedDataStructure()$clustResProbes$labels[session$userData$sessionVariables$traitReducedDataStructure()$clustResProbes$order]
      positions <- base::which(CpG %in% unlist(strsplit(txtSearchCpG, " ")))
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

getSearchResultTrait <- function(txtSearchTrait, session) {
  base::tryCatch(
    {
      base::print(base::paste0(sysTimePID(), " searching Trait"))
      #look into clustResTraits and find position of Trait
      Trait <- session$userData$sessionVariables$clustResTraits()$labels[session$userData$sessionVariables$traitReducedDataStructure()$clustResTraits()$order]
      positions <- base::which(Trait %in% unlist(strsplit(txtSearchTrait, " ")))
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

plotCombinedHM_P_Val <- function(input, output, session) {
  base::print(base::paste0(sysTimePID(), " start plotting heatmap for P_Val."))
  output$txtHMDescription_P_Val <-
    shiny::renderText(base::paste0("calculating heatmap..., current plot is not valid"))
  while (!is.null(grDevices::dev.list())) {
    grDevices::dev.off()
  }
  base::print(base::paste0(sysTimePID(), " creating empty heatmap."))
  # combinedHMP_VAL <- emptyHM()
  # InteractiveComplexHeatmap::makeInteractiveComplexHeatmap(input = input, output = output,
  #                                                          session = session, ht_list = combinedHMP_VAL,
  #                                                          heatmap_id = "heatmap_1")
  combinedDFP_Val_Labels <- session$userData$sessionVariables$traitReducedDataStructure()$combinedDFP_Val_Labels
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
  base::print(base::class(dfP_Val))
  if (nrow(dfP_Val) > 5) {
    startTime <- Sys.time()
    output$txtResultingN <-
      shiny::renderText(paste0("number of resulting probes: ", nrow(dfP_Val)))
    base::print(base::paste0(sysTimePID(), " gc()"))
    base::tryCatch({
      # check clustResProbes > 8
      base::options(expressions = 500000)
      dendProbes <- session$userData$sessionVariables$traitReducedDataStructure()$probeDendrogram
      dendProbes <-
        dendextend::color_branches(dendProbes, as.integer(input$txtMaxClassesProbes))
      #dendTraits <- session$userData$sessionVariables$traitReducedDendTraits()
      dendTraits <- session$userData$sessionVariables$traitReducedDataStructure()$traitDendrogram

      base::print(base::paste0(sysTimePID(), " before calculating heatmap"))

      base::print(base::paste0(sysTimePID(), " length(unlist(dendProbes)): ", length(unlist(dendProbes))))
      base::print(base::paste0(sysTimePID(), " length(unlist(dendTraits)): ", length(unlist(dendTraits))))
      length(unlist(dendTraits)) == base::dim(combinedDFP_Val_Labels$dfP_Val)[2]
      length(unlist(dendProbes)) == base::dim(combinedDFP_Val_Labels$dfP_Val)[1]
      selectedRowIndicesYellow <- unlist(strsplit(input$txtSearchCpG, split = " ")) #is a list of cg-numbers from search field "txtSearchCpG"
      selectedColIndices <- unlist(strsplit(input$txtSearchTrait, split = " "))
      #selectedRowIndicesOrange <- session$userData$sessionVariables$distancesBelowThreshold()
      selectedRowIndicesOrange <- NULL
      base::print(base::paste0(sysTimePID(), " before l <- combinedDFInteractiveHeatMapP_Val(combinedDFP_Val_Labels, dendProbes, dendTraits, selectedRowIndices, selectedColIndices)"))
#browser() #check, whether this is called twice
      l <-
        combinedDFInteractiveHeatMapP_Val(combinedDFP_Val_Labels, dendProbes, dendTraits, selectedRowIndicesYellow, selectedColIndices, selectedRowIndicesOrange, session)
      base::print(base::paste0(sysTimePID(), " before combinedHMP_VAL <- l$combinedHMP_VAL"))
      combinedHMP_VAL <- l$combinedHMP_VAL

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
        ht_list = combinedHMP_VAL,
        heatmap_id = "Heatmap_P_Val",
        show_layer_fun = TRUE,
        click_action = click_action_HM_P_Val,
        brush_action = brush_action_HM_P_Val,
        hover_action = hover_action_HM_P_Val
      )
#      session$userData$sessionVariables$SPLOM <- FALSE
    },
    error = function(e) {
      base::message("An error occurred in plotCombinedHM_P_Val():\n", e)
      Cstack_info()
      browser()
    },
    warning = function(w) {
      base::message("A warning occurred in plotCombinedHM_P_Val():\n", w)
      browser()
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

#' plotCombinedDWHM_P_Val
#' plots distance weighted heatmap
#' @param input shiny input
#' @param output shiny output
#' @param session shiny session
#' @return nothing
#' examples plotCombinedDWHM_P_Val(input, output, session)
plotCombinedDWHM_P_Val <- function(input, output, session) {
  base::print(base::paste0(sysTimePID(), " start plotting distance weighted heatmap for P_Val."))
  output$txtDWHMDescription_P_Val <-
    shiny::renderText(base::paste0("calculating DW heatmap..., current plot is not valid"))
  while (!is.null(grDevices::dev.list())) {
    grDevices::dev.off()
  }
  base::print(base::paste0(sysTimePID(), " creating empty heatmap."))
  combinedDFP_Val_Labels <- session$userData$sessionVariables$traitReducedDataStructure()$combinedDFP_Val_Labels

  dfP_Val <- combinedDFP_Val_Labels$dfP_Val
  #browser() #if step 3 was omitted, we see an error here...
#  dfP_Val[dfP_Val > 0.05] <- NA # 1
  base::print(
    base::paste0(
      sysTimePID(),
      " calculating combined heatmap with rows= ",
      nrow(dfP_Val),
      " cols= ",
      ncol(dfP_Val)
    )
  )
  base::print(base::class(dfP_Val))
  if (nrow(dfP_Val) > 5) {
    startTime <- Sys.time()
    base::print(base::paste0(sysTimePID(), " gc()"))
    gc()
    base::tryCatch({
      # check clustResProbes > 8
      #base::length(session$userData$sessionVariables$traitReducedDataStructure()$clustResProbes)
      #base::options(expression = 500000)
      base::options(expressions = 500000)
      dendProbes <- session$userData$sessionVariables$traitReducedDataStructure()$probeDendrogram
      dendProbes <-
        dendextend::color_branches(dendProbes, as.integer(input$txtMaxClassesProbes))
      dendTraits <- session$userData$sessionVariables$traitReducedDataStructure()$traitDendrogram

      base::print(base::paste0(sysTimePID(), " before calculating heatmap"))

      base::print(base::paste0(sysTimePID(), " length(unlist(dendProbes)): ", length(unlist(dendProbes))))
      base::print(base::paste0(sysTimePID(), " length(unlist(dendTraits)): ", length(unlist(dendTraits))))
      length(unlist(dendTraits)) == base::dim(combinedDFP_Val_Labels$dfP_Val)[2]
      length(unlist(dendProbes)) == base::dim(combinedDFP_Val_Labels$dfP_Val)[1]
      selectedRowIndicesYellow <- unlist(strsplit(input$txtSearchCpG, split = " "))
      base::message(base::paste0(sysTimePID(), " length(selectedRowIndicesYellow): ", length(selectedRowIndicesYellow)))
      selectedColIndices <- unlist(strsplit(input$txtSearchTrait, split = " "))
      selectedRowIndicesOrange <- session$userData$sessionVariables$distancesBelowThreshold()
      base::message(base::paste0(sysTimePID(), " length(selectedRowIndicesOrange): ", length(selectedRowIndicesOrange)))
      base::print(base::paste0(sysTimePID(), " before l <- combinedDFInteractiveHeatMapP_Val(combinedDFP_Val_Labels, dendProbes, dendTraits, selectedRowIndices, selectedColIndices)"))
      #browser() #check, whether this is called twice
      l <-
        combinedDFInteractiveHeatMapP_Val(combinedDFP_Val_Labels, dendProbes, dendTraits, selectedRowIndicesYellow, selectedColIndices, selectedRowIndicesOrange, session)
      base::print(base::paste0(sysTimePID(), " before combinedHMP_VAL <- l$combinedHMP_VAL"))
      combinedHMP_VAL <- l$DWcombinedHMP_VAL #l$combinedHMP_VAL

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
        ht_list = combinedHMP_VAL,
        heatmap_id = "DWHeatmap_P_Val",
        show_layer_fun = TRUE
#        click_action = click_action_HM_P_Val,
#        brush_action = brush_action_HM_P_Val,
#        hover_action = hover_action_HM_P_Val
      )
      #      session$userData$sessionVariables$SPLOM <- FALSE
    },
    error = function(e) {
      base::message("An error occurred in plotCombinedDWHM_P_Val():\n", e)
      Cstack_info()
      browser()
    },
    warning = function(w) {
      base::message("A warning occurred in plotCombinedDWHM_P_Val():\n", w)
      browser()
    },
    finally = {
      endTime <- Sys.time()
      elapsedTime <- endTime - startTime
      base::print(base::paste0(sysTimePID(), " after plotting heatmap for DWP_Val; InteractiveComplexHeatmap::makeInteractiveComplexHeatmap(). Elapsed time: ", elapsedTime, " sec."))
      output$txtDWHMDescription_P_Val <-
        shiny::renderText(
          base::paste0(
            sysTimePID(),
            " done plotting heatmap for DWP_Val..., current plot is valid. n(probe) = ",
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

#' plotCombinedCondDWHM_P_Val
#' plots condensed distance weighted heatmap
#' @param input shiny input
#' @param output shiny output
#' @param session shiny session
#' @return nothing
#' examples plotCombinedCondDWHM_P_Val(input, output, session)
plotCombinedCondDWHM_P_Val <- function(input, output, session) {
  base::tryCatch({
    base::print(base::paste0(sysTimePID(), " start plotting condensed distance weighted heatmap for P_Val."))
    output$txtCondDWHMDescription_P_Val <-
      shiny::renderText(base::paste0("calculating condensedDW heatmap..., current plot is not valid"))
    while (!is.null(grDevices::dev.list())) {
      grDevices::dev.off()
    }
    base::print(base::paste0(sysTimePID(), " creating empty heatmap."))
    combinedDFP_Val_Labels <- session$userData$sessionVariables$traitReducedDataStructure()$combinedDFP_Val_Labels
    dfP_Val <- combinedDFP_Val_Labels$dfP_Val
    #browser() #if step 3 was omitted, we see an error here...
###tbc()
    ### delete CpG with no neighbours here...
    ### see, whether dfP_Val has CpG names
    ### take only those, which have counterparts in DNAdistances:
#    DNAdistances <- session$userData$sessionVariables$distNeigboursProbes10()
#    DNAdistances <- na.omit(DNAdistances)
browser() # for whatever reason there is no reduction here....
#we need either to write this back to session$userData$sessionVariables$traitReducedDataStructure()
#or write the following function based on dfP_Val...
#we need to create dendProbes and dendTraits also new for the reduced variant...


    base::print(
      base::paste0(
        sysTimePID(),
        " calculating combined heatmap with rows= ",
        nrow(dfP_Val),
        " cols= ",
        ncol(dfP_Val)
      )
    )
    base::print(base::class(dfP_Val))
    if (nrow(dfP_Val) > 5) {
      startTime <- Sys.time()
      base::print(base::paste0(sysTimePID(), " gc()"))
      gc()
      # check clustResProbes > 8
      base::options(expressions = 500000)
      dendProbes <- session$userData$sessionVariables$traitReducedDataStructure()$probeDendrogram
browser() #we can either try to make a subset of dendrogram or to create a new dendrogram from the base data from the heatmap... we then need also a new clustering for the subset...
      dendProbes <-
        dendextend::color_branches(dendProbes, as.integer(input$txtMaxClassesProbes))
      dendTraits <- session$userData$sessionVariables$traitReducedDataStructure()$traitDendrogram
browser()
#dendTraits <- dendTraits[which(rownames(dendTraits) %in% DNAdistances$ID),]

      base::print(base::paste0(sysTimePID(), " before calculating condensed DW heatmap"))

      base::print(base::paste0(sysTimePID(), " length(unlist(dendProbes)): ", length(unlist(dendProbes))))
      base::print(base::paste0(sysTimePID(), " length(unlist(dendTraits)): ", length(unlist(dendTraits))))
#      length(unlist(dendTraits)) == base::dim(combinedDFP_Val_Labels$dfP_Val)[2]
#      length(unlist(dendProbes)) == base::dim(combinedDFP_Val_Labels$dfP_Val)[1]
      selectedRowIndicesYellow <- unlist(strsplit(input$txtSearchCpG, split = " "))
      base::message(base::paste0(sysTimePID(), " length(selectedRowIndicesYellow): ", length(selectedRowIndicesYellow)))
      selectedColIndices <- unlist(strsplit(input$txtSearchTrait, split = " "))
      selectedRowIndicesOrange <- session$userData$sessionVariables$distancesBelowThreshold()
      base::message(base::paste0(sysTimePID(), " length(selectedRowIndicesOrange): ", length(selectedRowIndicesOrange)))
      base::print(base::paste0(sysTimePID(), " before l <- combinedDFInteractiveHeatMapP_Val(combinedDFP_Val_Labels, dendProbes, dendTraits, selectedRowIndices, selectedColIndices)"))
      #browser() #check, whether this is called twice
      l <-
        combinedDFInteractiveHeatMapP_Val(combinedDFP_Val_Labels, dendProbes, dendTraits, selectedRowIndicesYellow, selectedColIndices, selectedRowIndicesOrange, session)
      base::print(base::paste0(sysTimePID(), " before combinedHMP_VAL <- l$combinedHMP_VAL"))
      combinedHMP_VAL <- l$DWcombinedHMP_VAL #l$combinedHMP_VAL

      endTime <- Sys.time()
      elapsedTime <- endTime - startTime
      base::print(base::paste0(sysTimePID(), " after calculating condensed DW heatmap. Elapsed time: ", elapsedTime, " sec."))
      base::print(base::paste0(sysTimePID(), " before plotting condensed DW heatmap; InteractiveComplexHeatmap::makeInteractiveComplexHeatmap()"))
      while (!base::is.null(grDevices::dev.list())) {
        grDevices::dev.off()
      }
      startTime <- Sys.time()
      InteractiveComplexHeatmap::makeInteractiveComplexHeatmap(
        input = input,
        output = output,
        session = session,
        ht_list = combinedHMP_VAL,
        heatmap_id = "condDWHeatmap_P_Val",
        show_layer_fun = TRUE
      )
    }
  },
  error = function(e) {
    base::message("An error occurred in plotCombinedDWHM_P_Val():\n", e)
    Cstack_info()
    browser()
  },
  warning = function(w) {
    base::message("A warning occurred in plotCombinedDWHM_P_Val():\n", w)
    browser()
  },
  finally = {
    endTime <- Sys.time()
    elapsedTime <- endTime - startTime
    base::print(base::paste0(sysTimePID(), " after plotting heatmap for DWP_Val; InteractiveComplexHeatmap::makeInteractiveComplexHeatmap(). Elapsed time: ", elapsedTime, " sec."))
    output$txtCondDWHMDescription_P_Val <-
      shiny::renderText(
        base::paste0(
          sysTimePID(),
          " done plotting heatmap for condensed DWP_Val..., current plot is valid. n(probe) = ",
          base::nrow(base::as.matrix(combinedDFP_Val_Labels[[1]])),
          "; n(trait) = ",
          base::ncol(base::as.matrix(combinedDFP_Val_Labels[[1]])),
          "; elapsed time: ",
          elapsedTime, " sec."
        )
      )
  })
}

#' click_action_HM_P_Val
#' @param df data.frame containing data which is enclosed by brush action
#' @param input shiny input object
#' @param output shiny output object
#' @param session shiny session object
#' @return nothing; function info label for click action in HM
# examples click_action_HM_P_Val(df, input, output, session)
click_action_HM_P_Val <- function(df, input, output, session) {
  base::tryCatch(
    {
      base::print(base::paste0(sysTimePID(), " start click_action_HM_P_Val()."))
      output[["info_HM_P_Val"]] <- shiny::renderUI({
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
    },
  error = function(e) {
    base::message("An error occurred in click_action_HM_P_Val():\n", e)
  },
  warning = function(w) {
    base::message("A warning occurred in click_action_HM_P_Val():\n", w)
  },
  finally = {
    base::print(base::paste0(sysTimePID(), " end click_action_HM_P_Val()."))
  }
  )
}

#' brush_action_HM_P_Val
#' @param df data.frame containing data which is enclosed by brush action
#' @param input shiny input object
#' @param output shiny output object
#' @param session shiny session object
#' @return nothing; function creates SPLOM
# examples brush_action_HM_P_Val(df, input, output, session)
brush_action_HM_P_Val <- function(df, input, output, session) {
  base::tryCatch(
    if (!is.null(df)) {
      row_index <- collapse::funique(unlist(df$row_index)) #row_index <- unique(unlist(df$row_index))
      column_index <- collapse::funique(unlist(df$column_index))
      #feed in selected CpG here
      session$userData$sessionVariables$selectedCpG(rownames(session$userData$sessionVariables$traitReducedDataStructure()$combinedDFP_Val_Labels$dfP_Val)[row_index])
      #add annotation
      rownames(session$userData$annotation) <- session$userData$annotation$name
      selectedAnnotation <- session$userData$annotation[session$userData$sessionVariables$selectedCpG(),]
      nprobes <- nrow(selectedAnnotation)
      selectedAnnotation$number <- seq(1:nprobes)
      selectedAnnotation$probeID <- selectedAnnotation$name
      col_order <- c("number", "probeID", "type", "target",	"name", "chromosome",	"position", "meth.dye", "gene.symbol", "gene.accession", "gene.region", "cpg.island.name", "relation.to.island", "snp.exclude", "450k", "common", "epic", "epic2")
      selectedAnnotation <- selectedAnnotation[, col_order]
      #add links to EWAS data hub
      selectedAnnotation <- addLinkToEWASDataHubShort(selectedAnnotation, session$userData$config$baseURL_EWASDataHub, session$userData$config$probeAttribut)
      selectedAnnotation <- addLinkToMRCEWASCatalogShort(selectedAnnotation, session$userData$config$baseURL_MRCEWASCatalog, session$userData$config$probeAttribut)

      selectedAnnotation <- addLinkToEWASDataHub(selectedAnnotation, session$userData$config$baseURL_EWASDataHub, session$userData$config$probeAttribut)
      selectedAnnotation <- addLinkToMRCEWASCatalog(selectedAnnotation, session$userData$config$baseURL_MRCEWASCatalog, session$userData$config$probeAttribut)
      selectedAnnotation$probeID <- NULL
      mergedOriginDF <- session$userData$sessionVariables$traitReducedDataStructure()$combinedDFP_Val_Labels$mergedOriginDF[column_index]
      originTrait <- session$userData$sessionVariables$traitReducedDataStructure()$combinedDFP_Val_Labels$mergedOriginTrait[column_index]
      traitLabels <- colnames(session$userData$sessionVariables$traitReducedDataStructure()$combinedDFP_Val_Labels$dfP_Val)[column_index]
      selectedTrait <- cbind(traitLabels, originTrait)
      colnames(selectedTrait) <- c("traitName", "traitSource")
      session$userData$sessionVariables$selectedTrait(selectedTrait)
      #create DT from selectedAnnotation
      output$DTSelectedCpG <- DT::renderDataTable(as.data.frame(selectedAnnotation), escape = FALSE, extensions = c("Buttons"), options = list(dom = "Bfrtip", buttons = c("csv"), pageLength = 10000))
      #output$DTSelectedTrait <- DT::renderDataTable(as.data.frame(session$userData$sessionVariables$selectedTrait()), escape = FALSE, extensions = c("Buttons"), options = list(dom = "Bfrtip", buttons = c("csv"), pageLength = 10000))

      #session$userData$sessionVariables$selectedOriginalData(getSelectedOriginalData(session$userData$sessionVariables$pReducedcombinedDFP_Val_Labels(), row_index, column_index, session))
      session$userData$sessionVariables$selectedTraitReducedcombinedDFP_Val_Labels(getSelectedTraitReducedcombinedDFP_Val_Labels(session$userData$sessionVariables$traitReducedDataStructure()$dfP_Val_Labels, session$userData$sessionVariables$selectedCpG(), session$userData$sessionVariables$selectedTrait(), session))
      session$userData$sessionVariables$selectedOriginalData(getSelectedOriginalData(session$userData$sessionVariables$traitReducedDataStructure()$dfP_Val_Labels, session$userData$sessionVariables$selectedCpG(), session$userData$sessionVariables$selectedTrait(), session))

      #session$userData$sessionVariables$selectedOriginalDataTraits(getSelectedOriginalDataTraits(session$userData$sessionVariables$pReducedcombinedDFP_Val_Labels(), row_index, column_index, session))
      session$userData$sessionVariables$selectedOriginalDataTraits(getSelectedOriginalDataTraits(session$userData$sessionVariables$traitReducedDataStructure()$dfP_Val_Labels, session$userData$sessionVariables$selectedCpG(), session$userData$sessionVariables$selectedTrait(), session))
      #session$userData$sessionVariables$selectedOriginalDataProbes(getSelectedOriginalDataProbes(session$userData$sessionVariables$pReducedcombinedDFP_Val_Labels(), row_index, column_index, session))
      session$userData$sessionVariables$selectedOriginalDataProbes(getSelectedOriginalDataProbes(session$userData$sessionVariables$traitReducedDataStructure()$dfP_Val_Labels, session$userData$sessionVariables$selectedCpG(), session$userData$sessionVariables$selectedTrait(), session))
      if (!is.null(session$userData$sessionVariables$selectedOriginalData())) {
        FactorialVars <- getBinaryFactorialVars(session$userData$sessionVariables$selectedOriginalData())
        #if (!is.null(FactorialVars)) {
        if (is.valid(FactorialVars)) {
          shiny::updateSelectizeInput(
            session = session,
            inputId = "markingVar",
            choices = FactorialVars,
            server = TRUE
          )
          #session$userData$sessionVariables$markingVar <- unlist(FactorialVars[1])
          #session$userData$sessionVariables$markingVar(unlist(FactorialVars[1]))
          message(session$userData$sessionVariables$markingVar())
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
    },
    error = function(e) {
      base::message("An error occurred in brush_action_HM_P_Val():\n", e)
    },
    warning = function(w) {
      base::message("A warning occurred in brush_action_HM_P_Val():\n", w)
    },
    finally = {
    }
  )
}

#' hover_action_HM_P_Val
#' @param df data.frame containing data which is enclosed by brush action
#' @param input shiny input object
#' @param output shiny output object
#' @param session shiny session object
#' @return nothing; function info label for click action in HM
# examples click_action_HM(df, input, output, session)
hover_action_HM_P_Val <- function(df, input, output, session) {
  base::tryCatch(
    {
      base::print(base::paste0(sysTimePID(), " hover_action_HM_P_Val.", as.character(head(df))))
    },
    error = function(e) {
      base::message("An error occurred in hover_action_HM():\n", e)
    },
    warning = function(w) {
      base::message("A warning occurred in hover_action_HM():\n", w)
    },
    finally = {
      base::print(base::paste0(sysTimePID(), " end hover_action_HM_P_Val()."))
    }
  )
}

# click_action_HM_DM <- function(df, input, output, session) {}
# brush_action_HM_DM <- function(df, input, output, session) {}
# hover_action_HM_DM <- function(df, input, output, session) {}
