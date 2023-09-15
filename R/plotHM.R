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
  base::print(base::paste0(sysTimePID(), " start getSelectedOriginalData()"))
  base::tryCatch(
    {
      colInd <- which(combinedDFP_Val_Labels$mergedColnames %in% column_index)
      selectedColnames <- combinedDFP_Val_Labels$mergedColnames[colInd] # we here have the colnames of the selected traits
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
        base::message(base::paste0(sysTimePID(), "nrow(selectedDF) = ",
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
      message("An error occurred in getSelectedOriginalData():\n", e)
    },
    warning = function(w) {
      message("A warning occurred in getSelectedOriginalData():\n", w)
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
      selectedColnames <- combinedDFP_Val_Labels$mergedColnames[colInd] # we here have the colnames of the selected traits
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
      message("An error occurred in getSelectedOriginalDataTraits():\n", e)
    },
    warning = function(w) {
      message("A warning occurred in getSelectedOriginalDataTraits():\n", w)
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
      message("An error occurred in getSelectedOriginalDataProbes():\n", e)
    },
    warning = function(w) {
      message("A warning occurred in getSelectedOriginalDataProbes():\n", w)
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
#' @return heatmap object for InteractiveComplexHeatmap::makeInteractiveComplexHeatmap
# examples combinedDFInteractiveHeatMapP_Val(combinedDF_Labels, dendProbes, dendTraits)
combinedDFInteractiveHeatMapP_Val <-
  function(combinedDF_Labels,
           #           dendProbes = NA,
           #          dendTraits = NA, # {
           dendProbes = NA,
           dendTraits = NA,
           selectedRowIndices = NA,
           selectedColIndices = NA) {
    #function is not called twice (check caller function), if so, check here
    base::tryCatch(
      {
        startTime <- Sys.time()
        base::print(base::paste0(sysTimePID(), " start preparing HM; combinedDFInteractiveHeatMapP_Val()"))
        mat <- base::as.matrix(combinedDF_Labels$dfP_Val)
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
        base::print(base::paste0(sysTimePID(), " preparing heatmap n(probes)=", base::dim(mat)[1], " x n(traits)=", base::dim(mat)[2]))
        if ("dendrogram" %in% class(dendProbes) && "dendrogram" %in% class(dendTraits)) {
          base::print(base::paste0(sysTimePID(), " length(unlist(dendProbes)): ", length(unlist(dendProbes))))
          base::print(base::paste0(sysTimePID(), " length(unlist(dendTraits)): ", length(unlist(dendTraits))))
          length(unlist(dendTraits)) == base::dim(mat)[2]
          length(unlist(dendProbes)) == base::dim(mat)[1]
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
                #mark selected row indices
                labelsProbes  <- labels(dendProbes)
                rowsToMark <- labelsProbes %in% selectedRowIndices
                if (any(rowsToMark)) {
                  grid::grid.rect(x, y[rowsToMark], w, h[rowsToMark], gp = grid::gpar( #grid::grid.rect(x[rowsToMark], y[rowsToMark], w[rowsToMark], h[rowsToMark], gp = grid::gpar(
                    fill = "yellow",
                    col = "yellow"
                  ))
                }
#                 #mark selected column indices
#                 labelsTraits  <- labels(dendTraits)
# #                selectedColIndices <- c("MeP_in_ug_g","MeP_in_ug_l")
#                 colsToMark <- labelsTraits %in% selectedColIndices
#                 if (any(colsToMark)) {
#                   grid::grid.rect(x[colsToMark], y, w[colsToMark], h, gp = grid::gpar(
#                     fill = "yellow",
#                     col = "yellow"
#                   ))
#                 }
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
        message("An error occurred in combinedDFInteractiveHeatMapP_Val():\n", e)
      },
      warning = function(w) {
        message("A warning occurred in combinedDFInteractiveHeatMapP_Val():\n", w)
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
        message("An error occurred in clearing grDevices():\n", e)
      },
      warning = function(w) {
        message("A warning occurred in clearing grDevices():\n", w)
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
        ht <- ComplexHeatmap::draw(ht, annotation_legend_list = lgd)
      },
      error = function(err) {
        base::print(base::paste0(sysTimePID(), " Error: unable to draw HM. ", err$message))
      },
      warning = function(w) {
        base::print(base::paste0(sysTimePID(), " Error: unable to wraw HM. ", w$message))
      },
      finally = {
        base::print(base::paste0(sysTimePID(), " grDevices::dev.off()"))
        grDevices::dev.off()
        base::print(base::paste0(sysTimePID(), " l <- base::list()"))
        l <- base::list()
        base::print(base::paste0(sysTimePID(), " l$combinedHMP_VAL <- ht"))
        l$combinedHMP_VAL <- ht
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
      CpG <- session$userData$sessionVariables$clustResProbes()$labels[session$userData$sessionVariables$clustResProbes()$order]
      positions <- base::which(CpG %in% unlist(strsplit(txtSearchCpG, " ")))
    },
    error = function(e) {
      message("An error occurred in getSearchResultCpG():\n", e)
    },
    warning = function(w) {
      message("A warning occurred in getSearchResultCpG():\n", w)
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
      Trait <- session$userData$sessionVariables$clustResTraits()$labels[session$userData$sessionVariables$clustResTraits()$order]
      positions <- base::which(Trait %in% unlist(strsplit(txtSearchTrait, " ")))
    },
    error = function(e) {
      message("An error occurred in getSearchResultTrait():\n", e)
    },
    warning = function(w) {
      message("A warning occurred in getSearchResultTrait():\n", w)
    },
    finally = {
      base::print(base::paste0(sysTimePID(), " end getSearchResultTrait()."))
      base::return(positions)
    }
  )
}

plotCombinedHM <- function(input, output, session) {
  base::print(base::paste0(sysTimePID(), " start plotting heatmap."))
  output$txtHMDescription <-
    shiny::renderText(base::paste0("calculating heatmap..., current plot is not valid"))
  while (!is.null(grDevices::dev.list())) {
    grDevices::dev.off()
  }
  base::print(base::paste0(sysTimePID(), " creating empty heatmap."))
  # combinedHMP_VAL <- emptyHM()
  # InteractiveComplexHeatmap::makeInteractiveComplexHeatmap(input = input, output = output,
  #                                                          session = session, ht_list = combinedHMP_VAL,
  #                                                          heatmap_id = "heatmap_1")
  #          combinedDFP_Val_Labels <- session$userData$sessionVariables$pReducedcombinedDFP_Val_Labels()
  combinedDFP_Val_Labels <- session$userData$sessionVariables$traitReducedcombinedDFP_Val_Labels()

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
      base::length(session$userData$sessionVariables$clustResProbes())
      base::options(expression = 500000)
      dendProbes <- session$userData$sessionVariables$dendProbes()
      dendProbes <-
        dendextend::color_branches(dendProbes, as.integer(input$txtMaxClassesProbes))
      dendTraits <- session$userData$sessionVariables$traitReducedDendTraits()

      base::print(base::paste0(sysTimePID(), " before calculating heatmap"))

      base::print(base::paste0(sysTimePID(), " length(unlist(dendProbes)): ", length(unlist(dendProbes))))
      base::print(base::paste0(sysTimePID(), " length(unlist(dendTraits)): ", length(unlist(dendTraits))))
      length(unlist(dendTraits)) == base::dim(combinedDFP_Val_Labels$dfP_Val)[2]
      length(unlist(dendProbes)) == base::dim(combinedDFP_Val_Labels$dfP_Val)[1]
      selectedRowIndices <- unlist(strsplit(input$txtSearchCpG, split = " "))
      selectedColIndices <- unlist(strsplit(input$txtSearchTrait, split = " "))
      base::print(base::paste0(sysTimePID(), " before l <- combinedDFInteractiveHeatMapP_Val(combinedDFP_Val_Labels, dendProbes, dendTraits, selectedRowIndices, selectedColIndices)"))
#browser() #check, whether this is called twice
      l <-
        combinedDFInteractiveHeatMapP_Val(combinedDFP_Val_Labels, dendProbes, dendTraits, selectedRowIndices, selectedColIndices)
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
        heatmap_id = "heatmap_1",
        show_layer_fun = TRUE,
        click_action = click_action_HM,
        brush_action = brush_action_HM,
        hover_action = hover_action_HM
      )
#      session$userData$sessionVariables$SPLOM <- FALSE
    },
    error = function(e) {
      message("An error occurred in plotCombinedHM():\n", e)
      browser()
    },
    warning = function(w) {
      message("A warning occurred in plotCombinedHM():\n", w)
      browser()
    },
    finally = {
      endTime <- Sys.time()
      elapsedTime <- endTime - startTime
      base::print(base::paste0(sysTimePID(), " after plotting heatmap; InteractiveComplexHeatmap::makeInteractiveComplexHeatmap(). Elapsed time: ", elapsedTime, " sec."))
      output$txtHMDescription <-
        shiny::renderText(
          base::paste0(
            sysTimePID(),
            " done plotting heatmap..., current plot is valid. n(probe) = ",
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

#' click_action_HM
#' @param df data.frame containing data which is enclosed by brush action
#' @param input shiny input object
#' @param output shiny output object
#' @param session shiny session object
#' @return nothing; function info label for click action in HM
# examples click_action_HM(df, input, output, session)
click_action_HM <- function(df, input, output, session) {
  base::tryCatch(
    {
      base::print(base::paste0(sysTimePID(), " start click_action_HM()."))
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
    },
  error = function(e) {
    message("An error occurred in click_action_HM():\n", e)
  },
  warning = function(w) {
    message("A warning occurred in click_action_HM():\n", w)
  },
  finally = {
    base::print(base::paste0(sysTimePID(), " end click_action_HM()."))
  }
  )
}

#' brush_action_HM
#' @param df data.frame containing data which is enclosed by brush action
#' @param input shiny input object
#' @param output shiny output object
#' @param session shiny session object
#' @return nothing; function creates SPLOM
# examples brush_action_HM(df, input, output, session)
brush_action_HM <- function(df, input, output, session) {
  base::tryCatch(
    if (!is.null(df)) {
      row_index <- collapse::funique(unlist(df$row_index)) #row_index <- unique(unlist(df$row_index))
      column_index <- collapse::funique(unlist(df$column_index))
      #feed in selected CpG here
      session$userData$sessionVariables$selectedCpG(rownames(session$userData$sessionVariables$pReducedcombinedDFP_Val_Labels()$dfP_Val)[row_index])
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
      session$userData$sessionVariables$selectedTrait(colnames(session$userData$sessionVariables$pReducedcombinedDFP_Val_Labels()$dfP_Val)[column_index])
      #create DT from selectedAnnotation
      output$DTSelectedCpG <- DT::renderDataTable(as.data.frame(selectedAnnotation), escape = FALSE, extensions = c("Buttons"), options = list(dom = "Bfrtip", buttons = c("csv"), pageLength = 10000))
      output$DTSelectedTrait <- DT::renderDataTable(as.data.frame(session$userData$sessionVariables$selectedTrait()), escape = FALSE, extensions = c("Buttons"), options = list(dom = "Bfrtip", buttons = c("csv"), pageLength = 10000))
      #session$userData$sessionVariables$selectedOriginalData(getSelectedOriginalData(session$userData$sessionVariables$pReducedcombinedDFP_Val_Labels(), row_index, column_index, session))
      session$userData$sessionVariables$selectedOriginalData(getSelectedOriginalData(session$userData$sessionVariables$pReducedcombinedDFP_Val_Labels(), session$userData$sessionVariables$selectedCpG(), session$userData$sessionVariables$selectedTrait(), session))

      #session$userData$sessionVariables$selectedOriginalDataTraits(getSelectedOriginalDataTraits(session$userData$sessionVariables$pReducedcombinedDFP_Val_Labels(), row_index, column_index, session))
      session$userData$sessionVariables$selectedOriginalDataTraits(getSelectedOriginalDataTraits(session$userData$sessionVariables$pReducedcombinedDFP_Val_Labels(), session$userData$sessionVariables$selectedCpG(), session$userData$sessionVariables$selectedTrait(), session))
      #session$userData$sessionVariables$selectedOriginalDataProbes(getSelectedOriginalDataProbes(session$userData$sessionVariables$pReducedcombinedDFP_Val_Labels(), row_index, column_index, session))
      session$userData$sessionVariables$selectedOriginalDataProbes(getSelectedOriginalDataProbes(session$userData$sessionVariables$pReducedcombinedDFP_Val_Labels(), session$userData$sessionVariables$selectedCpG(), session$userData$sessionVariables$selectedTrait(), session))
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
      message("An error occurred in brush_action_HM():\n", e)
    },
    warning = function(w) {
      message("A warning occurred in brush_action_HM():\n", w)
    },
    finally = {
    }
  )
}

#' hover_action_HM
#' @param df data.frame containing data which is enclosed by brush action
#' @param input shiny input object
#' @param output shiny output object
#' @param session shiny session object
#' @return nothing; function info label for click action in HM
# examples click_action_HM(df, input, output, session)
hover_action_HM <- function(df, input, output, session) {
  base::tryCatch(
    {
      base::print(base::paste0(sysTimePID(), " hover_action_HM.", as.character(head(df))))
    },
    error = function(e) {
      message("An error occurred in hover_action_HM():\n", e)
    },
    warning = function(w) {
      message("A warning occurred in hover_action_HM():\n", w)
    },
    finally = {
      base::print(base::paste0(sysTimePID(), " end hover_action_HM()."))
    }
  )
}
