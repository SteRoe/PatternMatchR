#' combinedDFP_Val_Labels$dfP_Val - data frame with p-values from regression results
#' combinedDFP_Val_Labels$dfDM - data frame with delta methylations from regression results
#' combinedDFP_Val_Labels$dfN - data frame with n from regression results
#' combinedDFP_Val_Labels$labelsDF1 - variable labels from data frame for trait 1
#' combinedDFP_Val_Labels$labelsDF2 - variable labels from data frame for trait 2
#' combinedDFP_Val_Labels$labelsDF3 - variable labels from data frame for trait 3
#' combinedDFP_Val_Labels$labelsDF <- c(LabelsDF1, LabelsDF2, LabelsDF3) - combined variable labels for all three traits
#' combinedDFP_Val_Labels$mergedOriginDF - number (label) of original data frame
#' combinedDFP_Val_Labels$mergedColnames <- merged original Colnames
#' combinedDFP_Val_Labels$mergedOriginTrait <- number of trait (1,2 or 3), a particular variable belongs to

#' getResultDfP_D_N
#' @param listOfResultDF data.frame containing a list of data.frames to be merged
#' @param P_D_N scalar value "P", "D" or "N" describing whether to merge P_VAL, DeltaMeth or N
#' @return merged data.frame
# examples getResultDfP_D_N(listDF, "P")
getResultDfP_D_N <- function(listOfResultDF, P_D_N) {
  base::tryCatch(
    {
      base::print(base::paste0(sysTimePID(), " start getResultDfP_D_N(): ", P_D_N, "."))
      i <- NULL
      if (base::length(listOfResultDF) != 0) {
        #foreach::foreach(i = 1:length(listOfResultDF)) %do% {
        for (i in seq_along(listOfResultDF)) {
          base::print(base::paste0(
            sysTimePID(),
            " processing trait ",
            i,
            " of ",
            base::length(listOfResultDF)
          ))
          if (P_D_N == "P") {
            DF <- base::subset(listOfResultDF[[i]], select = c("probeID", "P_VAL"))
          } else if (P_D_N == "D") {
            DF <-
              base::subset(listOfResultDF[[i]], select = c("probeID", "DeltaMeth"))
          } else if (P_D_N == "N") {
            DF <- base::subset(listOfResultDF[[i]], select = c("probeID", "N"))
          }
          if (i == 1) {
            merged <- DF
            rownames(merged) <- DF$probeID
          } else {
            base::print(base::paste0(sysTimePID(), " merge"))
            # merge
            merged <-
              base::merge(
                x = merged,
                y = DF,
                by.x = "probeID",
                by.y = "probeID",
                all.x = TRUE,
                all.y = TRUE
              )
          }
          colnames(merged)[i + 1] <- base::names(listOfResultDF)[i]
        }
      }
      rownames(merged) <- merged$probeID
      merged$probeID <-
        NULL # here we have probeIDs as rownames, therefore this variable is no longer needed
      base::print(base::paste0(sysTimePID(), " finished getResultDfP_D_N()"))
    },
    error = function(e) {
      base::message("An error occurred in getResultDfP_D_N():\n", e)
    },
    warning = function(w) {
      base::message("A warning occurred in getResultDfP_D_N():\n", w)
    },
    finally = {
      base::print(base::paste0(sysTimePID(), " end getResultDfP_D_N(): ", P_D_N, "."))
      return(merged)
    }
  )
}

#' mergeDFP_Val_Labels
#' merges data.frames out of sessionVariable containing p-values, delta methylation values and N by probeID
#' @param resultDFListTrait1 trait structure1 to merge
#' @param resultDFListTrait2 trait structure2 to merge
#' @param resultDFListTrait3 trait structure3 to merge
#' @param minN minimum n
#' @return named list of data.frames, one df for merged P_Val, one for merged DeltaMethylation, one for merged N as well as labels:
#' @return result$dfP_Val for p-values
#' @return result$dfDM for delta methylation values
#' @return result$dfN for n
#' @return result$labelsDF1 for labels belonging to original df1
#' @return result$labelsDF2 for labels belonging to original df2
#' @return result$labelsDF3 for labels belonging to original df3
# examples mergeDFP_Val_Labels(resultDFListTrait1, resultDFListTrait2, resultDFListTrait3, minN)
mergeDFP_Val_Labels <- function(resultDFListTrait1, resultDFListTrait2, resultDFListTrait3, minN) {
  base::tryCatch(
    {
      base::print(base::paste0(sysTimePID(), " start mergeDFP_Val_Labels()."))
      # merge three df
      if (!base::is.null(resultDFListTrait1$resultDFP_Val)) {
        dfList <- base::list(
          dfP_Val = NULL,
          dfDM = NULL,
          dfN = NULL,
          resultOriginDF = NULL,
          resultColnames = NULL,
          resultOriginalColnames = NULL,
          listPHENOdata = NULL
        )
        dfList$dfP_Val <- resultDFListTrait1$resultDFP_Val
        dfList$dfDM <- resultDFListTrait1$resultDFDM
        dfList$dfN <- resultDFListTrait1$resultDFN
        dfList$resultOriginDF <- resultDFListTrait1$resultOriginDF
        dfList$resultColnames <- resultDFListTrait1$resultColnames
        dfList$resultOriginalColnames <- resultDFListTrait1$resultOriginalColnames
        dfList$listPHENOdata <- resultDFListTrait1$listPHENOdata
        dfP_Val1 <- dfList$dfP_Val
        if (base::exists("dfP_Val1")) {
          checkResultP_Val_cg(dfP_Val1)
        }
        dfDM1 <- dfList$dfDM
        dfN1 <- dfList$dfN
        OriginDF1 <- dfList$resultOriginDF
        Colnames1 <- dfList$resultColnames
        OriginalColnames1 <- dfList$resultOriginalColnames
        DFList1 <- dfList$listPHENOdata
        originTrait1 <- rep("1", length(Colnames1))
        if (!((base::nrow(dfP_Val1) > 0) && (base::ncol(dfP_Val1) > 0))) {
          base::message(base::paste0(sysTimePID(), "nrow(DFP_Val1) or ncol(DFP_Val1) == 0"))
        }
      } else {
        base::message(base::paste0(sysTimePID(), "DF1 is not valid"))
      }
      if (!base::is.null(resultDFListTrait2$resultDFP_Val)) {
        dfList <- base::list()
        dfList$dfP_Val <- resultDFListTrait2$resultDFP_Val
        dfList$dfDM <- resultDFListTrait2$resultDFDM
        dfList$dfN <- resultDFListTrait2$resultDFN
        dfList$resultOriginDF <- resultDFListTrait2$resultOriginDF
        dfList$resultColnames <- resultDFListTrait2$resultColnames
        dfList$resultOriginalColnames <- resultDFListTrait2$resultOriginalColnames
        dfList$listPHENOdata <- resultDFListTrait2$listPHENOdata
        dfP_Val2 <- dfList$dfP_Val
        if (base::exists("dfP_Val2")) {
          checkResultP_Val_cg(dfP_Val2)
        }
        dfDM2 <- dfList$dfDM
        dfN2 <- dfList$dfN
        OriginDF2 <- dfList$resultOriginDF
        Colnames2 <- dfList$resultColnames
        OriginalColnames2 <- dfList$resultOriginalColnames
        DFList2 <- dfList$listPHENOdata
        originTrait2 <- rep("2", length(Colnames2))
        if (!((base::nrow(dfP_Val2) > 0) && (base::ncol(dfP_Val2) > 0))) {
          base::message(base::paste0(sysTimePID(), "nrow(DF2) or ncol(DF2) == 0"))
        }
      } else {
        base::message(base::paste0(sysTimePID(), "DF2 is not valid"))
      }
      if (!base::is.null(resultDFListTrait3$resultDFP_Val)) {
        dfList <- base::list()
        dfList$dfP_Val <- resultDFListTrait3$resultDFP_Val
        dfList$dfDM <- resultDFListTrait3$resultDFDM
        dfList$dfN <- resultDFListTrait3$resultDFN
        dfList$resultOriginDF <- resultDFListTrait3$resultOriginDF
        dfList$resultColnames <- resultDFListTrait3$resultColnames
        dfList$resultOriginalColnames <- resultDFListTrait3$resultOriginalColnames
        dfList$listPHENOdata <- resultDFListTrait3$listPHENOdata
        dfP_Val3 <- dfList$dfP_Val
        if (base::exists("dfP_Val3")) {
          checkResultP_Val_cg(dfP_Val3)
        }
        dfDM3 <- dfList$dfDM
        dfN3 <- dfList$dfN
        OriginDF3 <- dfList$resultOriginDF
        Colnames3 <- dfList$resultColnames
        OriginalColnames3 <- dfList$resultOriginalColnames
        DFList3 <- dfList$listPHENOdata
        originTrait3 <- rep("3", length(Colnames3))
        if (!((base::nrow(dfP_Val3) > 0) && (base::ncol(dfP_Val3) > 0))) {
          base::message(base::paste0(sysTimePID(), "nrow(DF3) or ncol(DF3) == 0"))
        }
      } else {
        base::message(base::paste0(sysTimePID(), "DF3 is not valid"))
      }
      if (base::exists("dfP_Val1")) {
        rn <- base::rownames(dfP_Val1)
        dfP_Val1 <- base::as.data.frame(dfP_Val1)
        dfP_Val1$Row.names <- rn
        base::row.names(dfP_Val1) <- rn
        rn <- base::rownames(dfDM1)
        dfDM1 <- base::as.data.frame(dfDM1)
        dfDM1$Row.names <- rn
        base::rownames(dfDM1) <- rn
        rn <- base::rownames(dfN1)
        dfN1 <- base::as.data.frame(dfN1)
        dfN1$Row.names <- rn
        base::rownames(dfN1) <- rn

        rn <- base::rownames(dfP_Val1)
        mergedDFP_Val <- base::as.data.frame(dfP_Val1)
        mergedDFP_Val$Row.names <- rn
        # crazy error message here: left side converting to list? due to dfP_Val1 was of class matrix
        # mergedDFP_Val$Row.names <- rn
        rn <- base::rownames(dfDM1)
        mergedDFDM <- base::as.data.frame(dfDM1)
        mergedDFDM$Row.names <- rn
        rn <- base::rownames(dfN1)
        mergedDFN <- base::as.data.frame(dfN1)
        mergedDFN$Row.names <- rn
        mergedOriginDF <- OriginDF1
        mergedColnames <- Colnames1
        mergedOriginalColnames <- OriginalColnames1
        mergedOriginTrait <- originTrait1
        mergedDFList <- DFList1
      }
      if (base::exists("dfP_Val2")) {
        if (base::exists("mergedDFP_Val")) {
          rn <- base::rownames(dfP_Val2)
          dfP_Val2 <- base::as.data.frame(dfP_Val2)
          dfP_Val2$Row.names <- rn
          base::rownames(dfP_Val2) <- rn
          rn <- base::rownames(dfDM2)
          dfDM2 <- base::as.data.frame(dfDM2)
          dfDM2$Row.names <- rn
          base::rownames(dfDM2) <- rn
          rn <- base::rownames(dfN2)
          dfN2 <- base::as.data.frame(dfN2)
          dfN2$Row.names <- rn
          base::rownames(dfN2) <- rn

          mergedDFP_Val <-
            base::merge(
              mergedDFP_Val,
              dfP_Val2,
              by.x = "Row.names",
              by.y = "Row.names",
              all.x = FALSE,
              all.y = FALSE
            )
          mergedDFDM <-
            base::merge(
              mergedDFDM,
              dfDM2,
              by.x = "Row.names",
              by.y = "Row.names",
              all.x = FALSE,
              all.y = FALSE
            )
          mergedDFN <-
            base::merge(
              mergedDFN,
              dfN2,
              by.x = "Row.names",
              by.y = "Row.names",
              all.x = FALSE,
              all.y = FALSE
            )
          mergedOriginDF <- c(mergedOriginDF, OriginDF2)
          mergedColnames <- c(mergedColnames, Colnames2)
          mergedOriginalColnames <- c(mergedOriginalColnames, OriginalColnames2)
          mergedOriginTrait <- c(mergedOriginTrait, originTrait2)
          mergedDFList <- c(mergedDFList, DFList2)
          if ("Row.names" %in% base::colnames(dfP_Val2)) {
            rownames(dfP_Val2) <- dfP_Val2$Row.names
            dfP_Val2$Row.names <- NULL
          }
          if ("Row.names" %in% base::colnames(dfDM2)) {
            rownames(dfDM2) <- dfDM2$Row.names
            dfDM2$Row.names <- NULL
          }
          if ("Row.names" %in% base::colnames(dfN2)) {
            rownames(dfN2) <- dfN2$Row.names
            dfN2$Row.names <- NULL
          }
        } else {
          rn <- base::rownames(dfP_Val2)
          mergedDFP_Val <- base::as.data.frame(dfP_Val2)
          mergedDFP_Val$Row.names <- rn
          base::rownames(mergedDFP_Val) <- rn
          rn <- base::rownames(dfDM2)
          mergedDFDM <- base::as.data.frame(dfDM2)
          mergedDFDM$Row.names <- rn
          base::rownames(mergedDFDM) <- rn
          rn <- base::rownames(dfN2)
          mergedDFN <- base::as.data.frame(dfN2)
          mergedDFN$Row.names <- rn
          base::rownames(mergedDFN) <- rn
          mergedOriginDF <- OriginDF2
          mergedColnames <- Colnames2
          mergedOriginalColnames <- OriginalColnames2
          mergedOriginTrait <- originTrait2
          mergedDFList <- DFList2
        }
      }
      if (base::exists("dfP_Val3")) {
        if (base::exists("mergedDFP_Val")) {
          rn <- base::rownames(dfP_Val3)
          dfP_Val3 <- base::as.data.frame(dfP_Val3)
          dfP_Val3$Row.names <- rn
          base::rownames(dfP_Val3) <- rn
          rn <- base::rownames(dfDM3)
          dfDM3 <- base::as.data.frame(dfDM3)
          dfDM3$Row.names <- rn
          base::rownames(dfDM3) <- rn
          rn <- base::rownames(dfN3)
          dfN3 <- base::as.data.frame(dfN3)
          dfN3$Row.names <- rn
          base::rownames(dfN3) <- rn
          mergedDFP_Val <-
            base::merge(
              mergedDFP_Val,
              dfP_Val3,
              by.x = "Row.names",
              by.y = "Row.names",
              all.x = FALSE,
              all.y = FALSE
            )
          mergedDFDM <-
            base::merge(
              mergedDFDM,
              dfDM3,
              by.x = "Row.names",
              by.y = "Row.names",
              all.x = FALSE,
              all.y = FALSE
            )
          mergedDFN <-
            base::merge(
              mergedDFN,
              dfN3,
              by.x = "Row.names",
              by.y = "Row.names",
              all.x = FALSE,
              all.y = FALSE
            )
          mergedOriginDF <- c(mergedOriginDF, OriginDF3)
          mergedColnames <- c(mergedColnames, Colnames3)
          mergedOriginalColnames <- c(mergedOriginalColnames, Colnames3)
          mergedOriginTrait <- c(mergedOriginTrait, originTrait3)
          mergedDFList <- c(mergedDFList, DFList3)
          if ("Row.names" %in% base::colnames(dfP_Val3)) {
            rownames(dfP_Val3) <- dfP_Val3$Row.names
            dfP_Val3$Row.names <- NULL
          }
          if ("Row.names" %in% base::colnames(dfDM3)) {
            rownames(dfDM3) <- dfDM3$Row.names
            dfDM3$Row.names <- NULL
          }
          if ("Row.names" %in% base::colnames(dfN3)) {
            rownames(dfN3) <- dfN3$Row.names
            dfN3$Row.names <- NULL
          }
        } else {
          rn <- base::rownames(dfP_Val3)
          mergedDFP_Val <- base::as.data.frame(dfP_Val3)
          mergedDFP_Val$Row.names <- rn
          base::rownames(mergedDFP_Val) <- rn
          rn <- base::rownames(dfDM3)
          mergedDFDM <- base::as.data.frame(dfDM3)
          mergedDFDM$Row.names <- rn
          base::rownames(mergedDFDM) <- rn
          rn <- base::rownames(dfN3)
          mergedDFN <- base::as.data.frame(dfN3)
          mergedDFN$Row.names <- rn
          base::rownames(mergedDFN) <- rn
          mergedOriginDF <- OriginDF3
          mergedColnames <- Colnames3
          mergedOriginalColnames <- OriginalColnames3
          mergedOriginTrait <- originTrait3
          mergedDFList <- DFList3
        }
      }
#browser() #check mergedDFList
      if (base::exists("mergedDFP_Val")) {
        if ("Row.names" %in% base::colnames(mergedDFP_Val)) {
          rownames(mergedDFP_Val) <- mergedDFP_Val$Row.names
          mergedDFP_Val$Row.names <- NULL
        }
        if ("Row.names" %in% base::colnames(mergedDFDM)) {
          rownames(mergedDFDM) <- mergedDFDM$Row.names
          mergedDFDM$Row.names <- NULL
        }
        if ("Row.names" %in% base::colnames(mergedDFN)) {
          rownames(mergedDFN) <- mergedDFN$Row.names
          mergedDFN$Row.names <- NULL
        }
      }
      if (base::exists("dfP_Val1")) {
        if ("Row.names" %in% base::colnames(dfP_Val1)) {
          dfP_Val1$Row.names <- NULL
        }
        splitPointStart <- 1
        splitPointEnd <- base::ncol(dfP_Val1)
        LabelsDF1 <-
          base::colnames(mergedDFP_Val)[splitPointStart:splitPointEnd]
      }
      if (base::exists("dfP_Val2")) {
        if ("Row.names" %in% base::colnames(dfP_Val2)) {
          dfP_Val2$Row.names <- NULL
        }
        if (base::exists("splitPointEnd")) {
          splitPointStart <- splitPointEnd + 1
        } else {
          splitPointStart <- 1
        }
        splitPointEnd <- splitPointStart + base::ncol(dfP_Val2) - 1
        LabelsDF2 <-
          base::colnames(mergedDFP_Val)[splitPointStart:splitPointEnd]
      }
      if (base::exists("dfP_Val3")) {
        if ("Row.names" %in% base::colnames(dfP_Val3)) {
          dfP_Val3$Row.names <- NULL
        }
        if (base::exists("splitPointEnd")) {
          splitPointStart <- splitPointEnd + 1
        } else {
          splitPointStart <- 1
        }
        splitPointEnd <- splitPointStart + base::ncol(dfP_Val3) - 1
        LabelsDF3 <-
          base::colnames(mergedDFP_Val)[splitPointStart:splitPointEnd]
      }
      if (base::exists("mergedDFP_Val")) {
        result <- base::list(dfP_Val = NULL, dfDM = NULL, dfN = NULL,
                             labelsDF1 = NULL, labelsDF2 = NULL,
                             labelsDF3 = NULL, mergedOriginDF = NULL,
                             mergedColnames = NULL, mergedOriginalColnames = NULL, mergedOriginTrait = NULL,
                             mergedDFList = NULL)
        result$dfP_Val <- mergedDFP_Val # matP_Val
        mergedDFDM <- base::abs(mergedDFDM) # all Values to positive values
        result$dfDM <- mergedDFDM # matDM
        result$dfN <- mergedDFN # matN
        if (base::exists("LabelsDF1")) {
          result$labelsDF1 <- LabelsDF1
        }
        if (base::exists("LabelsDF2")) {
          result$labelsDF2 <- LabelsDF2
        }
        if (base::exists("LabelsDF3")) {
          result$labelsDF3 <- LabelsDF3
        }
        result$mergedOriginDF <- mergedOriginDF
        result$mergedColnames <- mergedColnames
        result$mergedOriginalColnames <- mergedOriginalColnames
        result$mergedOriginTrait <- mergedOriginTrait
        result$mergedDFList <- mergedDFList
      }
      else {
        result <- NULL
      }
    },
    error = function(e) {
      base::message("An error occurred in mergeDFP_Val_Labels():\n", e)
    },
    warning = function(w) {
      base::message("A warning occurred in mergeDFP_Val_Labels():\n", w)
    },
    finally = {
      base::print(base::paste0(sysTimePID(), " size of merged data.frame: ", dim(result$mergedOriginDF), " ."))
      base::print(base::paste0(sysTimePID(), " end mergeDFP_Val_Labels()."))
      return(result)
    }
  )
}

#' loadDir
#' @param session session info
#' @param traitDirList list containing directories to load into PatternMatchR
#' @return data.frame with contents of traitDirList
# examples loadDir(session, traitDirList)
loadDir <- function(session, traitDirList) {
  base::tryCatch(
    {
      #load all trait folders
      base::print(base::paste0(sysTimePID(), " before loading '", traitDirList, "'."))
      traitDFs <-
        base::lapply(traitDirList,
                     FUN = loadResultDF,
                     session = session, loadRDS = TRUE
        )
      #merge all loaded folders
      base::print(base::paste0(sysTimePID(), " before merge folders()"))
      resultDFList <- loadtraitDFs(traitDFs)
      if (base::exists("resultDFList$resultDFP_Val")) {
        #browser()
        checkResultP_Val_cg(resultDFList$resultDFP_Val)
      }
    },
    error = function(e) {
      base::message("An error occurred in loadDir():\n", e)
    },
    warning = function(w) {
      base::message("A warning occurred in loadDir():\n", w)
    },
    finally = {
      base::print(base::paste0(sysTimePID(), " '", traitDirList, "' successfully loaded."))
      return(resultDFList)
    }
  )
}

#' updateTxtLoadOut
#' @param resultDFListTrait1 data.frame containing results from trait1 (red)
#' @param resultDFListTrait2 data.frame containing results from trait2 (green)
#' @param resultDFListTrait3 data.frame containing results from trait3 (blue)
#' @return HTML to show in result info line after data load section of PatternMatchR
# examples updateTxtLoadOut(resultDFListTrait1, resultDFListTrait2, resultDFListTrait3)
updateTxtLoadOut <- function(session, resultDFListTrait1, resultDFListTrait2, resultDFListTrait3) {
  base::tryCatch(
    {
      i <- NULL
      result <- NULL
      if (is.valid(resultDFListTrait1)) {
        listPHENOData <- resultDFListTrait1$listPHENOdata
        DFP_Val <- resultDFListTrait1$resultDFP_Val
        res <- NULL
        foreach::foreach(i = seq_along(listPHENOData)) %do% {
          res <- base::paste0(res, listPHENOData[[i]]$PHENOFileName, "; ")
        }
        folder <- resultDFListTrait1$folder
        result <- base::paste0("loaded data from trait 1 list (red) from folder ", folder , " with pheno file: ", res, "nrow=", nrow(DFP_Val), "; ncol=", ncol(DFP_Val), ".\n")
      }
      if (is.valid(resultDFListTrait2)) {
        listPHENOData <- resultDFListTrait2$listPHENOdata
        DFP_Val <- resultDFListTrait2$resultDFP_Val
        res <- NULL
        foreach::foreach(i = seq_along(listPHENOData)) %do% {
          res <- base::paste0(res, listPHENOData[[i]]$PHENOFileName, "; ")
        }
        folder <- resultDFListTrait2$folder
        result <- base::paste0(result, "loaded data from trait 2 list (green) from folder ", folder, " with pheno file: ", res, "nrow=", nrow(DFP_Val), "; ncol=", ncol(DFP_Val), ".\n")
      }
      if (is.valid(resultDFListTrait3)) {
        listPHENOData <- resultDFListTrait3$listPHENOdata
        DFP_Val <- resultDFListTrait3$resultDFP_Val
        res <- NULL
        foreach::foreach(i = seq_along(listPHENOData)) %do% {
          res <- base::paste0(res, listPHENOData[[i]]$PHENOFileName, "; ")
        }
        folder <- resultDFListTrait3$folder
        result <- base::paste0(result, "loaded data from trait 3 list (blue) from folder ", folder, " with pheno file: ", res, "nrow=", nrow(DFP_Val), "; ncol=", ncol(DFP_Val), ".\n")
      }
    },
    error = function(e) {
      base::message("An error occurred in updateTxtLoadOut():\n", e)
    },
    warning = function(w) {
      base::message("A warning occurred in updateTxtLoadOut():\n", w)
    },
    finally = {
      return(shiny::HTML(result))
    }
  )
}
#' updateTxtMergeOut
#' @param combinedDFP_Val_Labels data.frame containing merged results
#' @return HTML to show in result info line after data merge section of PatternMatchR
# examples updateTxtMergeOut(combinedDFP_Val_Labels)
updateTxtMergeOut <- function(combinedDataStructure) {
  base::tryCatch(
    {
      result <- NULL
      #if (is.valid(combinedDataStructure) && nrow(combinedDataStructure$combinedDFP_Val_Labels$dfP_Val) > 0) {
      if (is.valid(combinedDataStructure$combinedDFP_Val_Labels$dfP_Val)) {
        result <- base::paste0("merge successful. result table is: nrow ",
                               nrow(combinedDataStructure$combinedDFP_Val_Labels$dfP_Val),
                               "; ncol: ", ncol(combinedDataStructure$combinedDFP_Val_Labels$dfP_Val))
      }
      else {
        base::print(base::paste0(sysTimePID(), " is.valid(combinedDFP_Val_Labels) == FALSE."))
      }
    },
    error = function(e) {
      base::message("An error occurred in updateTxtMergeOut():\n", e)
    },
    warning = function(w) {
      base::message("A warning occurred in updateTxtMergeOut():\n", w)
    },
    finally = {
      return(shiny::HTML(result))
    }
  )
}

updateSliders <- function(session, combinedDFP_Val_Labels) {
  DF <- combinedDFP_Val_Labels$dfP_Val
  DF <- as.matrix(DF)
  minP <- base::apply(DF, 2, FUN = function(x) {base::min(x[x > 0], na.rm = TRUE)})
  minP <- base::min(minP)
  minP <- extractMantissaExponent(minP)$exponent #base::round(extractMantissaExponent(minP)$exponent, 5)
  maxP <- base::apply(DF, 2, FUN = function(x) {base::max(x[x > 0], na.rm = TRUE)})
  maxP <- base::max(maxP)
  maxP <- extractMantissaExponent(maxP)$exponent #base::round(extractMantissaExponent(maxP)$exponent, 5)
  shiny::updateSliderInput(session = session, inputId = "sldP_Val", min = minP, max = maxP, value = c(minP, maxP))
  DF <- combinedDFP_Val_Labels$dfDM
  DF <- as.matrix(DF)
  minDM <- base::apply(DF, 2, FUN = function(x) {base::min(x[x > 0], na.rm = TRUE)})
  minDM <- base::min(minDM)
  if (minDM < 0) {
    base::message(base::paste0(sysTimePID(), " Warning: minDM < 0. Please check your data.")) #that should not be the case, please check data!
    minDM <- 0
  }
  maxDM <- base::apply(DF, 2, FUN = function(x) {base::max(x[x > 0], na.rm = TRUE)})
  maxDM <- base::max(maxDM)
  if (maxDM > 1) {
    base::message(base::paste0(sysTimePID(), "Warning: maxDM > 1. Please check your data.")) #that should not be the case, please check data!
    maxDM <- 1
  }
  #shiny::updateSliderInput(session = session, inputId = "sldDM", min = minDM, max = maxDM, value = c(minDM, maxDM), step = NULL)
  shiny::updateSliderInput(session = session, inputId = "sldDM", min = minDM, max = maxDM, value = c(minDM, maxDM), step = 0.001)
  DF <- combinedDFP_Val_Labels$dfN
  DF <- as.matrix(DF)
  minN <- base::apply(DF, 2, FUN = function(x) {base::min(as.integer(x[x > 0]), na.rm = TRUE)})
  minN <- base::min(minN)
  if (minN < 1) {
    base::message(base::paste0(sysTimePID(), "Warning: minN < 1. Please check your data.")) #that should not be the case, please check data!
    minN <- 1
  }
  if (minN != as.integer(minN)) {
    base::message(base::paste0(sysTimePID(), "Warning: minN != as.integer(minN). Please check your data.")) #that should not be the case, please check data!
    minN <- as.integer(minN)
  }
  maxN <- base::apply(DF, 2, FUN = function(x) {base::max(as.integer(x[x > 0]), na.rm = TRUE)})
  maxN <- base::max(maxN)
  if (maxN != as.integer(maxN)) {
    base::message(base::paste0(sysTimePID(), "Warning: maxN != as.integer(maxN). Please check your data.")) #that should not be the case, please check data!
    maxN <- as.integer(maxN)
  }
  if (maxN < 1) {
    base::message(base::paste0(sysTimePID(), "Warning: maxN < 1. Please check your data.")) #that should not be the case, please check data!
    browser()
  }
  shiny::updateSliderInput(session = session, inputId = "sldN", min = minN, max = maxN, value = c(minN, maxN))
}

# mergeDFP_Val_Labels <- compiler::cmpfun(mergeDFP_Val_Labels)

