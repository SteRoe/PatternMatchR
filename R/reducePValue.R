#' reduceTraitData
#' reduces data in structure combinedDFP_Val_Labels by p-value range and minimum n; entire CpGs will be removed
#' @param combinedDFP_Val_Labels data structure which should become p_val-reduced
#' @param minP_Val minimum p-value to retain
#' @param maxP_Val maximum p-value to retain
#' @param minN minimum n for model in regression result
#' @param debugMode return smaller data structure (only session$userData$sessionVariables$debugNumber records) for faster debugging
#' @return result list()
getPReducedTraitData <- function(session, combinedDFP_Val_Labels, minP_Val, maxP_Val, minDM, maxDM, minN, maxN, debugMode) {
  if (maxN < 1) {
    base::print(base::paste0(Sys.time(), "Warning: maxN < 1. Please check your data.")) #that should not be the case, please check data!
    browser()
  }
  base::tryCatch(
    {
      base::print(base::paste0(Sys.time(), " start pReduceTraitData()."))
      if (is.valid(combinedDFP_Val_Labels)) {
        result <- combinedDFP_Val_Labels
        dfP_Val <- result$dfP_Val
        dfDM <- result$dfDM
        dfN <- result$dfN
        base::print(base::paste0(Sys.time(),
                                 " result matrix before reduce has N(row) probes=",
                                 base::nrow(dfP_Val),
                                 " and n(col) traits=",
                                 base::ncol(dfP_Val), "."))
        LabelsDF1 <- result$labelsDF1
        LabelsDF2 <- result$labelsDF2
        LabelsDF3 <- result$labelsDF3
        mergedOriginDF <- result$mergedOriginDF
        mergedColnames <- result$mergedColnames
        mergedOriginTrait <- result$mergedOriginTrait
        mergedDFList <- result$mergedDFList
        # omit p_values from dfs - max PVal
        #if we obtain 0 hits later, check maxP_Val and minP_Val
        cgsToRetainMaxP <- dfP_Val < maxP_Val
        cgsToRetainMaxP <- unique(rownames(which(cgsToRetainMaxP == TRUE, arr.ind = TRUE)))

        cgsToRetainMinP <- dfP_Val > minP_Val
        cgsToRetainMinP <- unique(rownames(which(cgsToRetainMinP == TRUE, arr.ind = TRUE)))
        #here we obtain 0 hits
        if (base::exists("cgsToRetainMaxP") && base::exists("cgsToRetainMinP")) {
          cgsToRetainP <- base::intersect(cgsToRetainMaxP, cgsToRetainMinP)
        }
        else if (base::exists("cgsToRetainMaxP")) {
          cgsToRetainP <- cgsToRetainMaxP
        }
        else if (base::exists("cgsToRetainMinP")) {
          cgsToRetainP <- cgsToRetainMinP
        }
        if (length(cgsToRetainP) == 0) {
          base::print(base::paste0(Sys.time(), "Warning: length(cgsToRetainP) == 0. Please check your data.")) #that should not be the case, please check data!
          browser()
        }
        #take only DM outside slider defined range in cgsToRetainDM
#tbc() think on outside definition of DM range
        cgsToRetainMaxDM <- dfDM < maxDM #dfDM > maxDM
        cgsToRetainMaxDM <- unique(rownames(which(cgsToRetainMaxDM == TRUE, arr.ind = TRUE)))
        cgsToRetainMinDM <- dfDM > minDM #dfDM < minDM
        cgsToRetainMinDM <- unique(rownames(which(cgsToRetainMinDM == TRUE, arr.ind = TRUE)))
        cgsToRetainDM <- unique(c(cgsToRetainMaxDM, cgsToRetainMinDM))
        if (base::exists("cgsToRetainMaxDM") && base::exists("cgsToRetainMinDM")) {
          cgsToRetainDM <- base::intersect(cgsToRetainMaxDM, cgsToRetainMinDM)
        }
        else if (base::exists("cgsToRetainMaxDM")) {
          cgsToRetainDM <- cgsToRetainMaxDM
        }
        else if (base::exists("cgsToRetainMinDM")) {
          cgsToRetainDM <- cgsToRetainMinDM
        }
        if (length(cgsToRetainDM) == 0) {
          base::print(base::paste0(Sys.time(), "Warning: length(cgsToRetainDM) == 0. Please check your data.")) #that should not be the case, please check data!
          browser()
        }
        cgsToRetainMaxN <- dfN < maxN
        cgsToRetainMaxN <- unique(rownames(which(cgsToRetainMaxN == TRUE, arr.ind = TRUE)))
        cgsToRetainMinN <- dfN > minN
        cgsToRetainMinN <- unique(rownames(which(cgsToRetainMinN == TRUE, arr.ind = TRUE)))
        if (base::exists("cgsToRetainMaxN") && base::exists("cgsToRetainMinN")) {
          cgsToRetainN <- base::intersect(cgsToRetainMaxN, cgsToRetainMinN)
        }
        else if (base::exists("cgsToRetainMaxN")) {
          cgsToRetainN <- cgsToRetainMaxN
        }
        else if (base::exists("cgsToRetainMinN")) {
          cgsToRetainN <- cgsToRetainMinN
        }
        #intersect of all three cgsToRetain here:
        cgsToRetain <- intersect(cgsToRetainP, cgsToRetainDM)
        cgsToRetain <- intersect(cgsToRetain, cgsToRetainN)
        if (!base::exists("cgsToRetain") || length(cgsToRetain) == 0) {
          base::print(base::paste0(Sys.time(), "Warning: length(cgsToRetain) == 0. Please check your data.")) #that should not be the case, please check data!
          browser()
          base::print(base::paste0(Sys.time(), " max p-val border too low: ",
                                   maxP_Val, "; no remaining CpG."))
        }
        dfP_Val <- dfP_Val[cgsToRetain, ]
        dfDM <- dfDM[cgsToRetain, ]
        dfN <- dfN[cgsToRetain, ]
        base::print(base::paste0(Sys.time(),
                                 " result matrix after reduce has N(row) probes=",
                                 base::nrow(dfP_Val),
                                 " and n(col) traits=", base::ncol(dfP_Val),
                                 "."))
        #        dfP_Val[dfP_Val > 0.05] <- NA # 1
        base::print(base::class(dfP_Val))
        base::print(Cstack_info())
        if (base::nrow(dfP_Val) >= 5) {
          if (debugMode == TRUE && base::nrow(dfP_Val) > session$userData$sessionVariables$debugNumber) {
            base::print(base::paste0(Sys.time(), " debug mode n probes=session$userData$sessionVariables$debugNumber"))
            #dfP_Val <- dfP_Val[1:session$userData$sessionVariables$debugNumber, ]
            dfP_Val <- head(dfP_Val,session$userData$sessionVariables$debugNumber)
            dfDM <- head(dfDM,session$userData$sessionVariables$debugNumber)
            dfN <- head(dfN,session$userData$sessionVariables$debugNumber)
          }
          base::print(base::paste0(
            Sys.time(),
            " set infinite p-values to .Machine$integer.max."
          ))
          indexdfP_Val_infinite <- is.infinite(as.matrix(dfP_Val))
          dfP_Val[indexdfP_Val_infinite] <- .Machine$integer.max
          if (!base::is.data.frame(dfDM)) {
            dfDM <- base::as.data.frame(dfDM)
          }
          base::print(base::paste0(Sys.time(), " shortening dfDM"))
          dfDM <- dfDM[rownames(dfP_Val), colnames(dfP_Val)]
          if (!base::is.data.frame(dfN)) {
            dfN <- base::as.data.frame(dfN)
          }
          base::print(base::paste0(Sys.time(), " shortening dfN"))
          dfN <- dfN[rownames(dfP_Val), colnames(dfP_Val)]
          combinedDFP_Val_Labels <- base::list(dfP_Val = NULL, dfDM = NULL,
                                               dfN = NULL, labelsDF1 = NULL,
                                               labelsDF2 = NULL, labelsDF3 = NULL,
                                               mergedOriginDF = NULL, mergedColnames = NULL,
                                               mergedOriginTrait = NULL, mergedDFList = NULL)
          combinedDFP_Val_Labels$dfP_Val <- dfP_Val
          combinedDFP_Val_Labels$dfDM <- dfDM
          combinedDFP_Val_Labels$dfN <- dfN

          combinedDFP_Val_Labels$labelsDF1 <- LabelsDF1
          combinedDFP_Val_Labels$labelsDF2 <- LabelsDF2
          combinedDFP_Val_Labels$labelsDF3 <- LabelsDF3
          combinedDFP_Val_Labels$mergedOriginDF <- mergedOriginDF
          combinedDFP_Val_Labels$mergedColnames <- mergedColnames
          combinedDFP_Val_Labels$mergedOriginTrait <- mergedOriginTrait
          combinedDFP_Val_Labels$mergedDFList <- mergedDFList #(consists of mergedDFList$PHENODF and mergedDFList$PHENOFileName for each source file)
        }
        else {
          base::message(base::paste0(Sys.time(), " less than 5 probes remained. Please check your maximum and minimum p-val settings."))
        }
      }
    },
    error = function(e) {
      message("An error occurred in getPReducedTraitData():\n", e)
    },
    warning = function(w) {
      message("A warning occurred in getPReducedTraitData():\n", w)
    },
    finally = {
      #base::print(base::paste0(Sys.time(), " size of reduced data.frame: ", dim(combinedDFP_Val_Labels$dfP_Val), " ."))
      base::print(base::paste0(Sys.time(), " size of reduced data.frame: nrow (CpG)=",
                               nrow(combinedDFP_Val_Labels$dfP_Val), "; ncol (trait)=",
                               ncol(combinedDFP_Val_Labels$dfP_Val), " ."))
      base::print(base::paste0(Sys.time(), " pReduceTraitData() finished."))
      return(combinedDFP_Val_Labels)
    }
  )
}

#' updateTxtpReduceOut
#' generates summary text after p-value reduction
#' @param pReducedcombinedDFP_Val_Labels data structure with trait reduced results
#' @return text
#' examples updateTxtpReduceOut(pReducedcombinedDFP_Val_Labels)
updateTxtpReduceOut <- function(pReducedcombinedDFP_Val_Labels) {
  base::tryCatch(
    {
      result <- NULL
      if (is.valid(pReducedcombinedDFP_Val_Labels)) {
        result <- base::paste0("p reduce successful. result table is: nrow (CpG): ",
                               nrow(pReducedcombinedDFP_Val_Labels$dfP_Val),
                               "; ncol (trait): ", ncol(pReducedcombinedDFP_Val_Labels$dfP_Val))
      }
      else {
        base::message(base::paste0(Sys.time(), " is.valid(pReducedcombinedDFP_Val_Labels) == FALSE."))
      }
    },
    error = function(e) {
      message("An error occurred in updateTxtpReduceOut():\n", e)
    },
    warning = function(w) {
      message("A warning occurred in updateTxtpReduceOut():\n", w)
    },
    finally = {
      return(shiny::HTML(result))
    }
  )
}
